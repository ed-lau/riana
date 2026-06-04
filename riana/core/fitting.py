"""Kinetic curve fitting for D₂O time series (M3 Week 4).

Replaces the legacy ``riana_fit.py``'s analytic ``calculate_fs_m0`` path
(which uses a fixed per-peptide site count from sequence-count heuristics
and the analytic ``a_max = (1-ria_max)^n`` formula — biases k_deg ≈ −0.5
on the Week 0 baseline) with the IsoSpec forward-model FS solver lifted
in Phase F1.

For each peptide:

1. **Spep from coefficients** (`algorithms.isotope_dist.spep_from_coefficients`):
   ``Spep = Σ aa_coefficient[c] * count(c, sequence)`` using a pre-learned
   per-cell-line coefficient table (the M2 frozen
   ``d2o_aa_coefficients_<line>.csv``) or a literature default. This is
   the **production** path — Spep is deterministic from sequence given
   the coefficient table, not iteratively re-fit per-peptide.
2. **Per-timepoint FS** via :func:`algorithms.isotope_dist.solve_fs_d2o`
   using the computed Spep — full-envelope least-squares, not
   ``iso0 / Σ iso``.
3. **Kinetic curve fit** of ``(t, FS)`` to the selected model
   (simple / guan / fornasiero).
4. **Bootstrap CI** on ``k_deg`` by resampling ``(t, FS)`` pairs with
   replacement; replaces the legacy ``sqrt(diag(pcov))`` heuristic.

PROJECT_REVIEW.md §2b fixes landing here:

  - ``label == 4`` dispatch (the legacy string-check ``label == 'aa'``
    never fired because ``label`` is an int)
  - FS-denominator drift: gone — IsoSpec forward FS uses the full
    envelope shape, not ``iso0 / colsums``
  - Bootstrap CI replaces ``sqrt(diag(pcov))[0]``
  - Fixed-site-count model replaced: ``calculate_label_n`` is no longer
    called — Spep comes from the cell-line coefficient table

The legacy ``riana_fit.py`` stays as the ``--engine legacy`` regression
gate through Week 4; Phase F5 wires ``--engine new`` to dispatch here.
"""

from __future__ import annotations

import logging
import re
from concurrent import futures
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Callable, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from riana.algorithms.isotope_dist import solve_fs_d2o, spep_from_coefficients
from riana.algorithms.mass_calc import calculate_ion_mz
from riana.config import FitConfig
from riana.core import models

_LOGGER = logging.getLogger(__name__)

_MODELS: dict[str, Callable] = {
    "simple": models.one_exponent,
    "guan": models.two_compartment_guan,
    "fornasiero": models.two_compartment_fornasiero,
}

_K_DEG_INIT = 0.5
_K_DEG_BOUNDS = ([1e-4], [10.0])

_DEFAULT_N_BOOT = 200
_DEFAULT_BOOT_CI_PCT = (5.0, 95.0)


@dataclass(frozen=True, slots=True)
class FitResult:
    """Per-peptide fit output."""

    concat: str
    k_deg: float
    r_squared: float
    sd: float          # bootstrap std of k_deg
    spep: float        # per-residue-coefficient sum (deterministic per peptide)
    ci_lo: float       # bootstrap CI lower (default 5th pct)
    ci_hi: float       # bootstrap CI upper (default 95th pct)
    t: list[float]
    fs: list[float]
    protein_id: str


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def load_aa_coefficients(path: str | Path) -> dict[str, float]:
    """Load a ``d2o_aa_coefficients_<line>.csv`` into the dict ``fit_run`` expects.

    The CSV is the format committed under
    ``tests/data/calibration_d2o_mixing/<line>/`` — two required columns
    ``amino_acid, coefficient``. Other columns (bootstrap stats etc.) are
    ignored.
    """
    df = pd.read_csv(path)
    return {str(aa): float(c) for aa, c in zip(df["amino_acid"], df["coefficient"])}


def fit_run(
    config: FitConfig,
    integrate_dfs: list[pd.DataFrame],
    aa_coefficients: Mapping[str, float],
    *,
    n_boot: int = _DEFAULT_N_BOOT,
    boot_ci_pct: tuple[float, float] = _DEFAULT_BOOT_CI_PCT,
    random_state: int = 1337,
) -> pd.DataFrame:
    """Fit kinetic constants across a D₂O time series.

    Args:
        config: :class:`riana.config.FitConfig` (model, label, depth,
            q-value, k_p / k_r / r_p, ria_max, threads).
        integrate_dfs: one ``pandas.DataFrame`` per timepoint, each in the
            ``_riana.txt`` schema (with ``concat`` / ``sample`` /
            ``percolator q-value`` / ``isoN`` columns).
        aa_coefficients: ``{aa_letter: coefficient}`` for per-peptide Spep
            computation. Typically loaded via :func:`load_aa_coefficients`
            from the M2 per-cell-line frozen table; literature default for
            other cell types.
        n_boot: bootstrap resamples for the k_deg CI.
        boot_ci_pct: percentile bounds (default 5, 95).
        random_state: seed for the bootstrap RNG.

    Returns:
        DataFrame indexed by ``concat`` with columns
        ``k_deg, R_squared, sd, spep, ci_lo, ci_hi, t, fs, protein id``.
    """
    if config.model not in _MODELS:
        raise ValueError(
            f"unknown kinetic model {config.model!r}; "
            f"expected one of {sorted(_MODELS)}"
        )
    model_fn = _MODELS[config.model]

    rdf = pd.concat(integrate_dfs, ignore_index=True)
    rdf = rdf[rdf["percolator q-value"] < config.q_value].copy()
    rdf = (
        rdf.groupby(["concat", "file_idx"], group_keys=False)
        .filter(lambda x: x["sample"].nunique() >= config.depth)
        .copy()
    )
    if rdf.empty:
        raise ValueError(
            "No peptides survive --q-value / --depth filtering. "
            "Relax thresholds or check inputs."
        )

    concat_list = sorted(rdf["concat"].unique())
    fit_partial = partial(
        _fit_one_concat,
        rdf=rdf,
        model_fn=model_fn,
        config=config,
        aa_coefficients=dict(aa_coefficients),
        n_boot=n_boot,
        boot_ci_pct=boot_ci_pct,
        base_seed=random_state,
    )

    if config.threads <= 1:
        results = [fit_partial(c) for c in concat_list]
    else:
        with futures.ThreadPoolExecutor(max_workers=config.threads) as ex:
            results = list(ex.map(fit_partial, concat_list))

    _LOGGER.info(
        "fit_run: %d peptides processed, %d converged",
        len(results),
        sum(1 for r in results if r is not None and not np.isnan(r.k_deg)),
    )
    return _build_output_df(results)


# ---------------------------------------------------------------------------
# Per-peptide internals
# ---------------------------------------------------------------------------


def _fit_one_concat(
    concat: str,
    *,
    rdf: pd.DataFrame,
    model_fn: Callable,
    config: FitConfig,
    aa_coefficients: dict[str, float],
    n_boot: int,
    boot_ci_pct: tuple[float, float],
    base_seed: int,
) -> FitResult | None:
    """Per-peptide fit: Spep from coefficients → per-timepoint FS → k_deg."""
    peptide_rows = (
        rdf[rdf["concat"] == concat]
        .drop_duplicates(subset=["sample", "file_idx"])
        .copy()
    )
    if peptide_rows.empty:
        return _null_result(concat, "")

    iso_cols = sorted(
        [c for c in peptide_rows.columns if re.match(r"^iso\d+$", c)],
        key=lambda c: int(c[3:]),
    )
    if not iso_cols:
        return _null_result(concat, "")

    obs_matrix = peptide_rows[iso_cols].to_numpy(dtype=np.float64)
    t_arr = np.array(
        [float(re.sub(r"[^0-9.]", "", s)) for s in peptide_rows["sample"]],
        dtype=np.float64,
    )
    protein_id = (
        str(peptide_rows["protein id"].iloc[0])
        if "protein id" in peptide_rows.columns
        else ""
    )

    seq = concat.rsplit("_", 1)[0]
    try:
        pep_mass = calculate_ion_mz(seq)
    except (KeyError, ValueError):
        return _null_result(concat, protein_id)

    # Spep from coefficients — deterministic per peptide given the table.
    spep_float = spep_from_coefficients(seq, aa_coefficients)
    spep_int = max(1, int(round(spep_float)))

    # Per-timepoint FS via solve_fs_d2o, at the experiment's precursor
    # enrichment ``config.ria_max``. KeyError fires for peptides whose
    # sequence carries a non-canonical residue (B/X/U/etc.) absent from
    # the production aa_atoms table — drop those peptides cleanly rather
    # than crashing the whole fit.
    sums = obs_matrix.sum(axis=1)
    valid = sums > 0
    fs_arr = np.full_like(t_arr, np.nan)
    try:
        for i in range(len(t_arr)):
            if valid[i]:
                fs_arr[i] = solve_fs_d2o(
                    seq, pep_mass, obs_matrix[i], spep_int,
                    ria_max=float(config.ria_max), n_iso=len(iso_cols),
                )
    except (KeyError, ValueError):
        return _null_result(concat, protein_id)

    fit_mask = ~np.isnan(fs_arr) & valid
    if int(fit_mask.sum()) < config.depth:
        return _null_result(concat, protein_id)

    # Kinetic-model asymptotes: FS goes 0 → 1 (full pool turned over).
    # a_max here is the *kinetic* model's saturation level — always 1.0
    # for FS ∈ [0, 1], independent of the experiment's precursor
    # enrichment (which is config.ria_max, passed to solve_fs_d2o above).
    # The legacy curve_fit also used a_max=1.0 (riana_fit.py:368).
    kinetic_kwargs = dict(
        a_0=0.0, a_max=1.0,
        k_p=config.k_p, k_r=config.k_r, r_p=config.r_p,
    )

    try:
        popt, _ = curve_fit(
            partial(model_fn, **kinetic_kwargs),
            t_arr[fit_mask], fs_arr[fit_mask],
            bounds=_K_DEG_BOUNDS,
            p0=[_K_DEG_INIT],
            maxfev=2000,
        )
    except (RuntimeError, ValueError):
        return _null_result(concat, protein_id)
    k_deg = float(popt[0])

    pred = np.array([
        model_fn(ti, k_deg=k_deg, **kinetic_kwargs) for ti in t_arr[fit_mask]
    ])
    residuals = fs_arr[fit_mask] - pred
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((fs_arr[fit_mask] - np.mean(fs_arr[fit_mask])) ** 2))
    r_squared = float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot

    # Bootstrap CI on k_deg.
    rng = np.random.default_rng(base_seed ^ (hash(concat) & 0xFFFFFFFF))
    boot_ks: list[float] = []
    t_fit = t_arr[fit_mask]
    fs_fit = fs_arr[fit_mask]
    for _ in range(n_boot):
        idx = rng.integers(0, len(t_fit), size=len(t_fit))
        try:
            popt_b, _ = curve_fit(
                partial(model_fn, **kinetic_kwargs),
                t_fit[idx], fs_fit[idx],
                bounds=_K_DEG_BOUNDS,
                p0=[k_deg],
                maxfev=2000,
            )
            boot_ks.append(float(popt_b[0]))
        except (RuntimeError, ValueError):
            continue

    if len(boot_ks) >= 10:
        ci_lo, ci_hi = (float(p) for p in np.percentile(boot_ks, boot_ci_pct))
        sd = float(np.std(boot_ks, ddof=1))
    else:
        ci_lo = ci_hi = sd = float("nan")

    return FitResult(
        concat=concat,
        k_deg=k_deg,
        r_squared=r_squared,
        sd=sd,
        spep=float(spep_float),
        ci_lo=ci_lo,
        ci_hi=ci_hi,
        t=t_arr[fit_mask].tolist(),
        fs=fs_arr[fit_mask].tolist(),
        protein_id=protein_id,
    )


def _null_result(concat: str, protein_id: str) -> FitResult:
    """Sentinel for peptides we couldn't fit — kept in the output for census."""
    return FitResult(
        concat=concat, k_deg=float("nan"), r_squared=float("nan"),
        sd=float("nan"), spep=float("nan"),
        ci_lo=float("nan"), ci_hi=float("nan"),
        t=[], fs=[], protein_id=protein_id,
    )


# ---------------------------------------------------------------------------
# Output DataFrame
# ---------------------------------------------------------------------------


def _build_output_df(results: list[FitResult | None]) -> pd.DataFrame:
    """DataFrame in the legacy riana_fit_peptides.txt schema + Phase F2 adds."""
    rows = [
        {
            "concat": r.concat,
            "t": r.t,
            "fs": r.fs,
            "k_deg": r.k_deg,
            "R_squared": r.r_squared,
            "sd": r.sd,
            "spep": r.spep,
            "ci_lo": r.ci_lo,
            "ci_hi": r.ci_hi,
            "protein id": r.protein_id,
        }
        for r in results if r is not None
    ]
    return pd.DataFrame(rows).set_index("concat")
