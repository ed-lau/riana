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

  - amino-acid ``a_0`` dispatch (the legacy string-check ``label == 'aa'``
    never fired because ``label`` was an int at the time; the v1.1.0
    taxonomy cleanup makes the labels the strings ``"D2O"``/``"O18"``/``"AA"``)
  - FS-denominator drift: gone — IsoSpec forward FS uses the full
    envelope shape, not ``iso0 / colsums``
  - Bootstrap CI replaces ``sqrt(diag(pcov))[0]``
  - Fixed-site-count model replaced: ``calculate_label_n`` is no longer
    called — Spep comes from the cell-line coefficient table

The legacy ``riana_fit.py`` stays as the ``--engine legacy`` regression
gate through Week 4; Phase F5 wires ``--engine new`` to dispatch here.
"""

from __future__ import annotations

import hashlib
import importlib.resources
import logging
import re
from concurrent import futures
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Callable, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from riana import constants
from riana.algorithms.isotope_dist import (
    init_channel_masses, init_envelope_width,
    solve_fs_d2o, solve_fs_d2o_ds, solve_fs_o18, solve_fs_o18_ds,
    spep_from_coefficients, spep_from_length_coefficients,
)
from riana.algorithms.mass_calc import calculate_ion_mz, parse_unimod_ids
from riana.config import FitConfig
from riana.core import models
from riana.utils import strip_concat

_LOGGER = logging.getLogger(__name__)

_MODELS: dict[str, Callable] = {
    "simple": models.one_exponent,
    "guan": models.two_compartment_guan,
    "fornasiero": models.two_compartment_fornasiero,
    # Calibration recovery (mixing-proportion runs): a through-origin line FS≈f
    # fit by the same curve_fit path; ``k_deg`` holds the recovered slope, R² the
    # recovery quality. Dispatched automatically for experiment_type=="calibration"
    # (see core.pipeline.fit_project) and selectable via --model / the GUI combo.
    "calibration": models.calibration_line,
}

_K_DEG_INIT = 0.5
_K_DEG_BOUNDS = ([1e-4], [10.0])

_DEFAULT_N_BOOT = 200
_DEFAULT_BOOT_CI_PCT = (5.0, 95.0)


def _concat_seed(base_seed: int, concat: str) -> int:
    """A stable bootstrap seed for one peptide ``concat``.

    Uses a content hash (``hashlib``, **not** the salted built-in ``hash``) so
    the per-peptide bootstrap stream is identical across runs *and* across the
    process-pool fit (:attr:`riana.config.FitConfig.workers`) — a given peptide
    seeds the same regardless of which worker draws it or how many there are.
    Mirrors :func:`riana.core.protein._group_rng`.
    """
    digest = hashlib.blake2b(concat.encode(), digest_size=8).digest()
    return (base_seed ^ int.from_bytes(digest, "big")) & 0xFFFFFFFFFFFFFFFF


def _fit_key(concat: str) -> str:
    """The chemical-mod-stripped peptidoform identity that GROUPS peptidoforms
    into one turnover curve (M7 tier 1b).

    Removes the ``[UNIMOD:N]`` tokens of purely chemical mods
    (:data:`riana.constants.CHEMICAL_MODS` — Met-Ox) so a peptide's oxidized and
    unoxidized forms — which share the FS-vs-time signature (oxidation is
    post-synthesis, it does not reset the D₂O clock) — pool into a single fit.
    They are still **integrated separately** at their own clean m/z; only the fit
    *groups* them. Biological mods (phospho) stay in the key so they remain
    distinct units. A peptidoform with no chemical mod returns its own ``concat``
    unchanged → byte-identical to the pre-merge fit.
    """
    parts = concat.rsplit("_", 1)
    if len(parts) != 2:
        return concat
    seq_with_mods, charge = parts
    for unimod_id in constants.CHEMICAL_MODS:
        seq_with_mods = seq_with_mods.replace(f"[UNIMOD:{unimod_id}]", "")
    return f"{seq_with_mods}_{charge}"

#: Isotopomer channels the D2O envelope solver currently needs in the integrate
#: output. :func:`algorithms.isotope_dist.solve_fs_d2o` matches the observed
#: envelope against the IsoSpec forward model over the contiguous m0-m5 channels,
#: so all of these must be present (the legacy ``--iso 0 6`` pair is not enough,
#: and would silently misalign the observed envelope against the model). Matches
#: the ``riana integrate --iso`` default. A future ``--fs`` may let the fit use a
#: subset, but integrate should still emit the full set.
_REQUIRED_D2O_ISOTOPOMERS = (0, 1, 2, 3, 4, 5)

# --fs auto (init-width-keyed limited-isotopomer scoring; FitConfig.fs_auto).
# Per peptidoform the fit scores the leading FS_AUTO_BASE channels (iso0-3) — the
# flat default that tightens recovery for the bulk — but WIDENS to all captured
# channels for peptides whose natural-abundance envelope reaches index
# >= FS_AUTO_INIT_W_THRESHOLD, i.e. iso4-5 sit *inside* that envelope and carry
# clean, model-predicted signal at every timepoint (so iso0-3 would truncate real
# signal). Keyed on the natural (θ=0) init width, NOT the integrate-time N_ISO:
# init width is purely compositional (RIA-/θ-independent), so the threshold needs
# no per-experiment retuning.
#
# Derivation (D2O calibration mixing series ac16/cm/ipsc, RIA 6%, ground-truth θ;
# report 2026-06-23_adaptive_niso_limited_isotopomer.md): MAE(iso0-3) vs MAE(all)
# crosses over at init width 6 on ALL three lines (Δ ≈ 0 at 5, clearly +ve at 6);
# the resulting policy is a strict Pareto win over flat iso0-3 (better overall
# within±0.05 AND recovers the wide-tail MAE). See bench_niso_crossover.py.
#
# REVISIT IF: (a) unusually high RIA — init width keys on the *natural* envelope,
# so at high enrichment the *labelling* can push real signal into iso4-5 for
# peptides whose natural width is < 6; init-width keying deliberately ignores that
# (conservative), so a higher-enrichment study with ground truth might justify a
# lower threshold or an N_ISO/init-width blend. (b) Samples of only very long
# peptides — nearly all would widen, so the base/threshold split stops
# discriminating. (c) The endgame is per-channel abundance weighting (soft
# matcher), which dissolves the hard threshold entirely. Threshold kept as a
# named constant (not a user dial) so a revisit is a one-line, documented change.
FS_AUTO_BASE = 4
FS_AUTO_INIT_W_THRESHOLD = 6

#: Mass-defect fs_ds (the orthogonal DeuteRater second estimate; v1.1.0 item 1b).
#: Inverted on the low, well-measured spacing channels only — the high channels carry
#: a per-peptide model-vs-centroid reference offset (report 2026-06-25); iso4/5 stay
#: in the GUI for QC but out of the number. iso0's spacing is ≡0 (no signal).
THETA_DS_CHANNELS = (1, 2, 3)
#: ¹⁸O Δspacing channels. ¹⁸O is a +2 Da label so iso1 is the *only* channel free of
#: an ¹⁸O-bearing isotopolog (iso2 = one ¹⁸O; iso3 = one ¹⁸O + one ¹³C; iso4 = two
#: ¹⁸O / one ¹⁸O + two ¹³C all carry signal) — score iso2–4, skipping the flat iso1.
#: NB the ¹⁸O-vs-¹³C mass-defect difference is small (~2.5 mDa at iso2), so o18 fs_ds
#: is a low-signal completeness estimate, not a reliable second estimate (see
#: :func:`~riana.algorithms.isotope_dist.solve_fs_o18_ds`).
THETA_DS_CHANNELS_O18 = (2, 3, 4)


def _fs_ds_points(
    dspacing_rows: list[list[float]],
    t_fit: np.ndarray,
    *,
    seq: str,
    pep_mass: float,
    spep: int,
    charge: int,
    ria_max: float,
    n_iso: int,
    mods: tuple[int, ...],
    label: str,
) -> list[float]:
    """Per-timepoint mass-defect fs_ds by **nonlinear spacing inversion**
    (:func:`solve_fs_d2o_ds`) — the spacing analog of the intensity FS solve.

    For each fitted point, the per-channel Δspacing (:data:`THETA_DS_CHANNELS`,
    iso1-3) is **anchored to the unlabelled (t/f == 0) point** when one is present
    among the fitted points (subtract its Δspacing — removes the per-peptide
    reference offset; theory-referenced fallback otherwise), then inverted against
    the init↔final mixture-spacing curve. Replaces the earlier linear
    ``ΔSₓ/ΔSₓmax`` ratio, which over-read mid-range (the mixture mass-shift is
    concave in f — report 2026-06-25). NaN for a point with < 2 usable channels.
    """
    if not dspacing_rows:
        return []
    zero_idx = next((i for i, tt in enumerate(t_fit) if tt == 0.0), None)
    anchor = dspacing_rows[zero_idx] if zero_idx is not None else None
    # ¹⁸O (+2 Da) scores iso2–4 (iso1 is flat); D₂O scores iso1–3.
    is_o18 = label == "O18"
    channels = THETA_DS_CHANNELS_O18 if is_o18 else THETA_DS_CHANNELS
    out: list[float] = []
    for row in dspacing_rows:
        obs_ds: dict[int, float] = {}
        for k in channels:
            if k >= len(row):
                continue
            v = row[k] - (anchor[k] if anchor is not None else 0.0)
            if np.isfinite(v):
                obs_ds[k] = float(v)
        if len(obs_ds) < 2:
            out.append(float("nan"))
            continue
        if is_o18:
            out.append(solve_fs_o18_ds(
                seq, pep_mass, obs_ds, spep, charge,
                ria_max=ria_max, n_iso=n_iso, mods=mods))
        else:
            out.append(solve_fs_d2o_ds(
                seq, pep_mass, obs_ds, spep, charge,
                ria_max=ria_max, n_iso=n_iso, mods=mods, label=label))
    return out


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
    #: M7 Stage B proteoform suffix (e.g. ``pS34476``), constant per peptidoform;
    #: empty when no biological mod. The rollup appends it to the resolved accession.
    mod_sites: str
    #: M5 per-timepoint prediction-interval bounds, aligned 1:1 with ``t``/``fs``
    #: (the residual-bootstrap band; ``nan`` when the bootstrap did not converge).
    fs_lo: list[float]
    fs_hi: list[float]
    #: Per-point biological replicate, aligned 1:1 with ``t``/``fs`` — the
    #: protein rollup groups within (protein, labeling time, biological replicate).
    bio_rep: list[int]
    #: Per-point provenance, aligned 1:1 with ``t``/``fs``: ``"q_value"`` (direct
    #: ID) or ``"mbr"`` (transferred). Lets the rollup/GUI weight or colour points.
    evidence: list[str] = field(default_factory=list)
    #: Per-point Met-Ox flag, aligned 1:1 with ``t``/``fs`` (the point came from a
    #: chemical-mod peptidoform merged at fit). Carries the breakdown to the rollup.
    metox: list[bool] = field(default_factory=list)
    #: v1.1.0 item 1a — the orthogonal mass-defect QC substrate, aligned 1:1 with
    #: ``t``/``fs`` (one inner list per fitted point, per isotopomer channel, in
    #: **mDa**). ``dmass`` is the absolute shift ``obs_mz[k] − init_ref[k]`` (drift-
    #: sensitive); ``dspacing`` is the M0-internal ``(obs_mz[k]−obs_mz[0]) −
    #: (init_ref[k]−init_ref[0])`` (DeuteRater ΔSₓ, drift-robust; index 0 ≡ 0). Empty
    #: when the integrate output carries no ``iso{N}_obs_mz`` columns (QC is opt-in
    #: on the data). NaN per channel when that channel was not observed in-window.
    dmass: list[list[float]] = field(default_factory=list)
    dspacing: list[list[float]] = field(default_factory=list)
    #: v1.1.0 item 1b — the mass-defect θ (θ_ΔS), aligned 1:1 with ``t``/``fs``: a
    #: drift-robust SECOND fraction-new estimate from the Δspacing (iso0–3, anchored,
    #: weighted median + MAD), to cross-check the intensity FS — never the primary θ.
    #: NaN per point when < 2 channels are usable; empty when no obs_mz / no charge.
    fs_ds: list[float] = field(default_factory=list)
    #: Census of the **fitted** points (a curve's composition): total, the count
    #: from MBR transfers, the count from Met-Ox peptidoforms merged at fit (M7
    #: tier 1b), and the count that are neither ("clean"). A point that is both MBR
    #: and Met-Ox is in n_mbr and n_metox but not n_clean.
    n_points: int = 0
    n_mbr: int = 0
    n_metox: int = 0
    n_clean: int = 0


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


#: Package directory holding the bundled per-AA coefficient presets.
_COEFF_PKG = "riana.data.coefficients"


def available_coefficient_presets() -> list[str]:
    """Names of the bundled ``--coefficients`` presets (CSV stems).

    These ship inside the package (``riana/data/coefficients/*.csv``):
    ``deberneh_2025_rss`` (N_EH by the recommended RSS nested-fit — the in-vivo D₂O
    **default**) and ``ilchenko_2019`` (M0-M1 N_aa), both LC-MS in-vivo mouse tables
    that tighten within-protein k agreement vs the ``commerford_1983`` tritium
    literature (report 2026-06-25_d2o_coefficient_tables); plus the calibration-derived
    cell-line tables ``alamillo_2025_{ac16,ipsc,cm}``. A user may also pass a
    filesystem path to their own re-derived table.
    """
    try:
        root = importlib.resources.files(_COEFF_PKG)
    except (ModuleNotFoundError, FileNotFoundError):
        return []
    return sorted(
        p.name[:-4] for p in root.iterdir()
        if p.name.endswith(".csv")
    )


def load_aa_coefficients(path: str | Path) -> dict[str, float]:
    """Load a per-AA coefficient table into the dict ``fit_run`` expects.

    ``path`` is resolved as either a **bundled preset name** (one of
    :func:`available_coefficient_presets`, e.g. ``"deberneh_2025_rss"``)
    or a **filesystem path** to a CSV in the same format — two required columns
    ``amino_acid, coefficient`` (the format committed under
    ``tests/data/calibration_d2o_mixing/<line>/``). Extra columns (bootstrap
    stats etc.) are ignored.
    """
    name = str(path)
    if name in available_coefficient_presets():
        src = importlib.resources.files(_COEFF_PKG).joinpath(f"{name}.csv")
        with importlib.resources.as_file(src) as real_path:
            df = pd.read_csv(real_path)
    else:
        df = pd.read_csv(path)
    return {str(aa): float(c) for aa, c in zip(df["amino_acid"], df["coefficient"])}


def load_o18_coefficients(path: str | Path) -> dict[str, float]:
    """Load an ¹⁸O length-model coefficient table for ``--label o18``.

    Like :func:`load_aa_coefficients`, ``path`` is a **bundled preset name**
    (e.g. ``"juber_2026_o18_ac16"``) or a **filesystem path** to a CSV — here with columns
    ``feature, coefficient`` (the format written by
    ``tests/benchmark/bench_o18_coefficients.py``; features =
    :data:`riana.algorithms.isotope_dist.O18_LENGTH_FEATURES`). The returned
    ``{feature: coefficient}`` dict feeds
    :func:`riana.algorithms.isotope_dist.spep_from_length_coefficients`.
    """
    name = str(path)
    if name in available_coefficient_presets():
        src = importlib.resources.files(_COEFF_PKG).joinpath(f"{name}.csv")
        with importlib.resources.as_file(src) as real_path:
            df = pd.read_csv(real_path)
    else:
        df = pd.read_csv(path)
    return {str(f): float(c) for f, c in zip(df["feature"], df["coefficient"])}


def fit_run(
    config: FitConfig,
    integrate_dfs: list[pd.DataFrame],
    aa_coefficients: Mapping[str, float],
    *,
    time_column: str | None = None,
    n_boot: int = _DEFAULT_N_BOOT,
    boot_ci_pct: tuple[float, float] = _DEFAULT_BOOT_CI_PCT,
    random_state: int = 1337,
) -> pd.DataFrame:
    """Fit kinetic constants across a D₂O time series.

    Args:
        config: :class:`riana.config.FitConfig` (model, label, depth,
            q-value, k_p / k_r / r_p, ria_max, workers).
        integrate_dfs: one or more ``pandas.DataFrame`` in the ``_riana.txt``
            schema (``concat`` / ``sample`` / ``percolator q-value`` / ``isoN``).
        aa_coefficients: ``{aa_letter: coefficient}`` for per-peptide Spep
            computation. Typically loaded via :func:`load_aa_coefficients`
            from the M2 per-cell-line frozen table; literature default for
            other cell types.
        time_column: when given (the M6a manifest path), the kinetic-curve
            x-axis is read from this numeric column (e.g. ``"labeling_time"``,
            from the SDRF identity) and each row is treated as one
            already-recombined data point — :func:`riana.core.pipeline` does the
            fraction merge and replicate-point assembly upstream. When ``None``
            (the legacy/single-run path), the timepoint is parsed out of the
            ``sample`` string and data points are deduped on ``(sample,
            file_idx)``, preserving the pre-M6a behavior.
        n_boot: bootstrap resamples for the k_deg CI.
        boot_ci_pct: percentile bounds (default 5, 95).
        random_state: seed for the bootstrap RNG.

    Returns:
        DataFrame indexed by ``concat`` with columns
        ``k_deg, R_squared, sd, spep, ci_lo, ci_hi, t, fs, protein id``.
    """
    # config.label is validated to {"hw", "o18"} by FitConfig.__post_init__. hw
    # uses the per-AA D₂O Spep + solve_fs_d2o; o18 uses the length-model Spep +
    # solve_fs_o18 (the 3-isotope-O envelope). The dispatch is in _fit_one_concat.
    if config.model not in _MODELS:
        raise ValueError(
            f"unknown kinetic model {config.model!r}; "
            f"expected one of {sorted(_MODELS)}"
        )
    model_fn = _MODELS[config.model]

    rdf = pd.concat(integrate_dfs, ignore_index=True)

    # The D2O envelope solver needs the m0-m5 channels present and aligned;
    # guard here so an integrate run with too few / the wrong isotopomers (e.g.
    # the legacy `--iso 0 6` pair) fails with a clear message instead of
    # silently misaligning the observed envelope against the IsoSpec model.
    present_isos = {int(c[3:]) for c in rdf.columns if re.match(r"^iso\d+$", c)}
    # Required channels (Guard 1, run-level): the default full-envelope fit needs
    # iso0-5; with --fs limited-isotopomer scoring it needs exactly the scored
    # leading channels — so `integrate --iso 0 1 2 3` + `fit --fs 0 1 2 3` is valid
    # (you only need to extract what you score). Error if any are absent as a
    # COLUMN. (A short peptidoform whose envelope ends early — its high channel is
    # NaN, not column-absent — is the separate per-peptide clamp in solve_fs_d2o.)
    if config.fs_auto:
        # auto scores iso0-3 by default and widens *within* the captured channels,
        # so it needs only the base present; the widen target is whatever was
        # integrated (len(iso_cols)). Require the base leading channels.
        required = tuple(range(FS_AUTO_BASE))
        hint = f"--iso '{' '.join(map(str, required))}'"
    elif config.score_channels is not None:
        required = tuple(range(config.score_channels))
        hint = f"--iso '{' '.join(map(str, required))}'"
    else:
        required = _REQUIRED_D2O_ISOTOPOMERS
        hint = "--iso '0 1 2 3 4 5'"
    missing = [i for i in required if i not in present_isos]
    if missing:
        label_name = "o18" if config.label == "o18" else "D2O"
        raise ValueError(
            f"The {label_name} fit needs isotopomers {list(required)} present in "
            f"the integrate output, but {['iso%d' % i for i in missing]} are absent "
            f"(found {sorted(present_isos)}). Re-run `riana integrate {hint}` "
            f"or `--iso auto`."
        )

    rdf = rdf[rdf["percolator q-value"] < config.q_value].copy()
    # --exclude-mbr: drop match-between-runs transfers before fitting (the
    # with/without-MBR A/B lever). MBR points are kept by default. No-op when the
    # column is absent (Percolator / pre-MBR integrate runs).
    if config.exclude_mbr and "evidence" in rdf.columns:
        rdf = rdf[rdf["evidence"] != "mbr"].copy()
    # M7 tier 1b: group peptidoforms by their chemical-mod-stripped key so a
    # peptide's Met-Ox and unoxidized forms pool into one turnover curve (they are
    # integrated separately, then merged here). For data with no chemical mod this
    # is identity — ``fit_key == concat`` — so the fit is byte-identical. The depth
    # gate is applied to the *merged* group, so two forms that are each too shallow
    # alone can still qualify once pooled.
    rdf = rdf.assign(fit_key=rdf["concat"].map(_fit_key))
    if time_column is not None:
        if time_column not in rdf.columns:
            raise ValueError(
                f"time_column {time_column!r} not in the integrate frames "
                f"(have {list(rdf.columns)}). The pipeline must add it."
            )
        # depth = distinct labeling timepoints per merged peptidoform — the
        # kinetic-identifiability quantity (a one-exponent curve needs >= depth
        # distinct x to constrain k). NOT raw PSM rows: a peptidoform seen many
        # times at one timepoint is a single kinetic point, and counting rows
        # would admit curves a fit can't identify. Robust to PSM / replicate /
        # peptidoform / MBR multiplicity; the raw count stays visible via
        # n_points. (Harmonizes with the legacy path's sample.nunique below.)
        rdf = (
            rdf.groupby("fit_key", group_keys=False)
            .filter(lambda x: x[time_column].nunique() >= config.depth)
            .copy()
        )
    else:
        rdf = (
            rdf.groupby(["fit_key", "file_idx"], group_keys=False)
            .filter(lambda x: x["sample"].nunique() >= config.depth)
            .copy()
        )
    if rdf.empty:
        raise ValueError(
            "No peptides survive --q-value / --depth filtering. "
            "Relax thresholds or check inputs."
        )

    concat_list = sorted(rdf["fit_key"].unique())
    fit_partial = partial(
        _fit_one_concat,
        rdf=rdf,
        model_fn=model_fn,
        config=config,
        aa_coefficients=dict(aa_coefficients),
        time_column=time_column,
        n_boot=n_boot,
        boot_ci_pct=boot_ci_pct,
        base_seed=random_state,
    )

    if config.workers > 1:
        # Process-level parallelism — the real lever for the GIL-bound fit
        # (IsoSpec FS + per-peptide residual bootstrap). The shared inputs (the
        # frame, coefficients, ...) are pickled to each worker ONCE via the
        # initializer; only the lightweight concat strings cross per task. The
        # per-concat seed makes the result independent of worker count.
        init_args = (rdf, model_fn, config, dict(aa_coefficients),
                     time_column, n_boot, boot_ci_pct, random_state)
        chunk = max(1, len(concat_list) // (config.workers * 8))
        with futures.ProcessPoolExecutor(
            max_workers=config.workers,
            initializer=_init_fit_worker,
            initargs=init_args,
        ) as ex:
            results = list(ex.map(_fit_one_concat_worker, concat_list, chunksize=chunk))
    else:
        # Serial: the per-peptide IsoSpec FS + residual bootstrap is GIL-bound, so
        # threading it gave no speedup; `workers` (processes) is the only lever.
        results = [fit_partial(c) for c in concat_list]

    _LOGGER.info(
        "fit_run: %d peptides processed, %d converged",
        len(results),
        sum(1 for r in results if r is not None and not np.isnan(r.k_deg)),
    )
    out = _build_output_df(results)
    # M5: the tidy per-timepoint fraction-new table rides alongside the wide
    # per-peptide frame so callers (CLI, pipeline) can serialize it without
    # re-deriving θ. ``.attrs`` survives the ``.copy()`` / column-add that
    # ``fit_project`` does per curve.
    out.attrs["fractions_long"] = build_fractions_long(results)
    return out


# ---------------------------------------------------------------------------
# Per-peptide internals
# ---------------------------------------------------------------------------


# --- process-pool plumbing ---------------------------------------------------
# Each worker holds the shared fit inputs in a module global, set ONCE by the
# pool initializer, so only concat strings cross the boundary per task (not the
# 100k-row frame, per task). Module-level (not closures) so they pickle for
# ``ProcessPoolExecutor`` on spawn-start platforms (macOS).
_FIT_WORKER_STATE: dict[str, object] = {}


def _init_fit_worker(
    rdf, model_fn, config, aa_coefficients,
    time_column, n_boot, boot_ci_pct, base_seed,
) -> None:
    _FIT_WORKER_STATE.update(
        rdf=rdf, model_fn=model_fn, config=config,
        aa_coefficients=aa_coefficients, time_column=time_column,
        n_boot=n_boot, boot_ci_pct=boot_ci_pct, base_seed=base_seed,
    )


def _fit_one_concat_worker(concat: str) -> "FitResult | None":
    s = _FIT_WORKER_STATE
    return _fit_one_concat(
        concat, rdf=s["rdf"], model_fn=s["model_fn"], config=s["config"],
        aa_coefficients=s["aa_coefficients"], time_column=s["time_column"],
        n_boot=s["n_boot"], boot_ci_pct=s["boot_ci_pct"], base_seed=s["base_seed"],
    )


def _fit_one_concat(
    concat: str,
    *,
    rdf: pd.DataFrame,
    model_fn: Callable,
    config: FitConfig,
    aa_coefficients: dict[str, float],
    time_column: str | None = None,
    n_boot: int,
    boot_ci_pct: tuple[float, float],
    base_seed: int,
) -> FitResult | None:
    """Per-peptide fit: Spep from coefficients → per-timepoint FS → k_deg.

    ``concat`` here is the **fit key** (:func:`_fit_key`): for ordinary peptides
    it is the concat itself, but for a peptide seen both oxidized and unoxidized
    it is the chemical-mod-stripped identity shared by both forms (M7 tier 1b), so
    their per-row FS points pool into one curve. Each row still carries its own
    ``concat`` → its own mods / precursor mass / IsoSpec envelope.
    """
    if time_column is not None:
        # M6a manifest path: rows are already one point per (peptide, biorep,
        # timepoint); the x-axis is the numeric identity column, not the sample
        # string.
        peptide_rows = rdf[rdf["fit_key"] == concat].copy()
    else:
        # Dedup includes ``concat`` so a peptide's Ox and non-Ox forms at the same
        # sample stay as two points (they are distinct precursors); for a single
        # peptidoform this collapses to the old (sample, file_idx) dedup.
        peptide_rows = (
            rdf[rdf["fit_key"] == concat]
            .drop_duplicates(subset=["concat", "sample", "file_idx"])
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
    if time_column is not None:
        t_arr = peptide_rows[time_column].to_numpy(dtype=np.float64)
    else:
        t_arr = np.array(
            [float(re.sub(r"[^0-9.]", "", s)) for s in peptide_rows["sample"]],
            dtype=np.float64,
        )
    protein_id = (
        str(peptide_rows["protein id"].iloc[0])
        if "protein id" in peptide_rows.columns
        else ""
    )
    mod_sites = (
        str(peptide_rows["mod sites"].iloc[0])
        if "mod sites" in peptide_rows.columns
        and not pd.isna(peptide_rows["mod sites"].iloc[0])
        else ""
    )
    # Per-point biological replicate (manifest path); 1 on the legacy path,
    # where bioreps are not a concept. Aligned to peptide_rows / t_arr order so
    # ``[fit_mask]`` selects the same points as ``t``/``fs``.
    bio_rep_arr = (
        peptide_rows["biological_replicate"].to_numpy()
        if "biological_replicate" in peptide_rows.columns
        else np.ones(len(peptide_rows), dtype=int)
    )

    # The bare backbone and Spep are shared by every form under this fit key
    # (chemical-mod stripping leaves the same residues + the same biological
    # mods). ``seq`` (brackets and charge stripped) feeds spep_from_coefficients
    # and solve_fs_d2o's residue iteration; it requires pure AA letters.
    seq = strip_concat(concat)
    # Label dispatch: o18 uses the length-model Spep (b·(L-1) + carboxyl/amide/
    # hydroxyl side chains); hw (D₂O) uses the per-AA labile-H table.
    is_o18 = config.label == "o18"
    spep_float = (spep_from_length_coefficients(seq, aa_coefficients) if is_o18
                  else spep_from_coefficients(seq, aa_coefficients))
    spep_int = max(1, int(round(spep_float)))

    # Per-ROW peptidoform identity — its mods (parsed off the row's own concat)
    # and precursor mass, so each form's IsoSpec envelope reflects its actual
    # atom composition (an Ox row carries the extra O and a +16 m/z). For a
    # single-peptidoform fit this is one entry == the old single-concat path.
    # KeyError/ValueError here = a non-canonical residue (B/X/U) absent from the
    # production aa_atoms table → drop the peptide cleanly.
    row_concats = peptide_rows["concat"].to_numpy()
    # Per-point provenance, aligned to peptide_rows / t_arr: MBR transfer vs direct
    # ID, and whether the row is a Met-Ox (chemical-mod) form (its own concat
    # differs from the stripped fit_key `concat`). Feeds the curve census below.
    is_mbr = (
        peptide_rows["evidence"].to_numpy() == "mbr"
        if "evidence" in peptide_rows.columns
        else np.zeros(len(peptide_rows), dtype=bool)
    )
    is_metox = row_concats != concat
    try:
        forms = {
            c: (tuple(parse_unimod_ids(c.rsplit("_", 1)[0])),
                calculate_ion_mz(c.rsplit("_", 1)[0]))
            for c in set(row_concats)
        }
    except (KeyError, ValueError):
        return _null_result(concat, protein_id, mod_sites)

    # --fs auto: per-peptidoform score width keyed on the natural-abundance (θ=0)
    # init width — iso0-3 for typical peptides, all captured channels for ones
    # whose natural envelope is broad enough that iso4-5 carry clean signal. One
    # cached init-width lookup per form. (Flat --fs / full-envelope: one scalar.)
    if config.fs_auto:
        form_score = {
            rc: (FS_AUTO_BASE
                 if init_envelope_width(seq, pm, mods=md) < FS_AUTO_INIT_W_THRESHOLD
                 else len(iso_cols))
            for rc, (md, pm) in forms.items()
        }

    # Per-timepoint FS via solve_fs_d2o, at the experiment's precursor
    # enrichment ``config.ria_max``.
    sums = obs_matrix.sum(axis=1)
    valid = sums > 0
    fs_arr = np.full_like(t_arr, np.nan)
    try:
        for i in range(len(t_arr)):
            if valid[i]:
                mods_i, pep_mass_i = forms[row_concats[i]]
                score_channels_i = (form_score[row_concats[i]] if config.fs_auto
                                    else config.score_channels)
                fs_arr[i] = (solve_fs_o18 if is_o18 else solve_fs_d2o)(
                    seq, pep_mass_i, obs_matrix[i], spep_int,
                    ria_max=float(config.ria_max), n_iso=len(iso_cols),
                    mods=mods_i, score_channels=score_channels_i,
                )
    except (KeyError, ValueError):
        return _null_result(concat, protein_id, mod_sites)

    fit_mask = ~np.isnan(fs_arr) & valid
    if int(fit_mask.sum()) < config.depth:
        return _null_result(concat, protein_id, mod_sites)

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
        return _null_result(concat, protein_id, mod_sites)
    k_deg = float(popt[0])

    pred = np.array([
        model_fn(ti, k_deg=k_deg, **kinetic_kwargs) for ti in t_arr[fit_mask]
    ])
    residuals = fs_arr[fit_mask] - pred
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((fs_arr[fit_mask] - np.mean(fs_arr[fit_mask])) ** 2))
    r_squared = float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot

    # Unified residual bootstrap (fixed t-design). Resample the fit residuals,
    # refit k_deg, and from each refit derive BOTH the k_deg CI and a
    # per-timepoint *prediction interval*: model(t_i; k*) plus a freshly
    # resampled residual, so fs_lower/fs_upper reflect this peptide's
    # measurement scatter (the substrate the protein rollup weights on). Keeping
    # every timepoint in the design (vs. a pairs bootstrap that can drop one) is
    # more robust on a sparse 3-5 point turnover curve; on noise-free data the
    # residuals — hence the band — collapse to ~0.
    rng = np.random.default_rng(_concat_seed(base_seed, concat))
    t_fit = t_arr[fit_mask]
    fs_fit = fs_arr[fit_mask]
    n_pts = len(t_fit)
    boot_ks: list[float] = []
    boot_obs: list[np.ndarray] = []
    for _ in range(n_boot):
        fs_star = pred + residuals[rng.integers(0, n_pts, size=n_pts)]
        try:
            popt_b, _ = curve_fit(
                partial(model_fn, **kinetic_kwargs),
                t_fit, fs_star,
                bounds=_K_DEG_BOUNDS,
                p0=[k_deg],
                maxfev=2000,
            )
        except (RuntimeError, ValueError):
            continue
        k_b = float(popt_b[0])
        boot_ks.append(k_b)
        pred_b = np.asarray(
            model_fn(t_fit, k_deg=k_b, **kinetic_kwargs), dtype=np.float64
        )
        boot_obs.append(pred_b + residuals[rng.integers(0, n_pts, size=n_pts)])

    if len(boot_ks) >= 10:
        ci_lo, ci_hi = (float(p) for p in np.percentile(boot_ks, boot_ci_pct))
        sd = float(np.std(boot_ks, ddof=1))
        lo_arr, hi_arr = np.percentile(np.vstack(boot_obs), boot_ci_pct, axis=0)
        fs_lo = [float(x) for x in lo_arr]
        fs_hi = [float(x) for x in hi_arr]
    else:
        ci_lo = ci_hi = sd = float("nan")
        fs_lo = [float("nan")] * n_pts
        fs_hi = [float("nan")] * n_pts

    # Census of the fitted points (aligned to fit_mask): MBR / Met-Ox / clean.
    mbr_fit = is_mbr[fit_mask]
    metox_fit = is_metox[fit_mask]

    # v1.1.0 item 1a — orthogonal mass-defect QC, per fitted point × channel (mDa).
    # The observed per-channel m/z (iso{k}_obs_mz, written by integrate) compared to
    # the init (θ=0) averaged-isotopolog reference recomputed HERE from IsoSpec —
    # path-independent, NOT the recorded ppm_error (which is init-anchored on
    # --iso auto but an analytic comb on the fixed path). ``dmass`` = absolute shift
    # obs−init_ref; ``dspacing`` = M0-internal (obs[k]−obs[0])−(ref[k]−ref[0]), the
    # drift-robust DeuteRater ΔSₓ. Empty when the integrate output predates the
    # obs_mz columns (QC is opt-in on the data); NaN per channel not seen in-window.
    obs_mz_cols = [f"{c}_obs_mz" for c in iso_cols]
    dmass_fit: list[list[float]] = []
    dspacing_fit: list[list[float]] = []
    if all(c in peptide_rows.columns for c in obs_mz_cols):
        obs_mz_matrix = peptide_rows[obs_mz_cols].to_numpy(dtype=np.float64)
        charge_arr = (
            peptide_rows["charge"].to_numpy(dtype=np.float64)
            if "charge" in peptide_rows.columns
            else np.full(len(peptide_rows), np.nan)
        )
        n_ch = len(iso_cols)
        nan_row = [float("nan")] * n_ch
        for i in np.nonzero(fit_mask)[0]:
            z = charge_arr[i]
            mods_i, pep_mass_i = forms[row_concats[i]]
            if not np.isfinite(z) or z <= 0:
                dmass_fit.append(list(nan_row))
                dspacing_fit.append(list(nan_row))
                continue
            ref_neutral = np.asarray(
                init_channel_masses(seq, pep_mass_i, n_ch, mods=mods_i),
                dtype=np.float64,
            )
            init_ref_mz = (ref_neutral + z * constants.PROTON_MASS) / z
            obs = obs_mz_matrix[i]
            dmass = (obs - init_ref_mz) * 1e3
            dspacing = ((obs - obs[0]) - (init_ref_mz - init_ref_mz[0])) * 1e3
            dmass_fit.append([float(x) for x in dmass])
            dspacing_fit.append([float(x) for x in dspacing])

    # v1.1.0 item 1b — the mass-defect fs_ds: an orthogonal, drift-robust SECOND
    # estimate of fraction-new from the per-channel Δspacing, cross-checking the
    # intensity FS (never displacing it). Anchored to the unlabelled point when
    # present, then NONLINEARLY inverted against the init↔final mixture spacing
    # curve (solve_fs_d2o_ds / solve_fs_o18_ds — the spacing analog of the intensity
    # solve; the old linear ΔS/ΔSₓmax over-read mid-range, report 2026-06-25).
    # Label-dispatched: D₂O scores iso1–3, ¹⁸O scores iso2–4 (only iso1 is free of an
    # ¹⁸O isotopolog). The ¹⁸O variant is a low-signal completeness estimate (small
    # ¹⁸O-vs-¹³C mass-defect difference) — see solve_fs_o18_ds.
    fs_ds_fit: list[float] = []
    if dspacing_fit:
        first_i = int(np.nonzero(fit_mask)[0][0])
        mods0, pep_mass0 = forms[row_concats[first_i]]
        z0 = int(charge_arr[first_i]) if np.isfinite(charge_arr[first_i]) else 0
        if z0 > 0:
            try:
                fs_ds_fit = _fs_ds_points(
                    dspacing_fit, t_fit, seq=seq, pep_mass=pep_mass0,
                    spep=spep_int, charge=z0, ria_max=float(config.ria_max),
                    n_iso=len(iso_cols), mods=mods0,
                    label=("O18" if is_o18 else "D2O"),
                )
            except (KeyError, ValueError):
                fs_ds_fit = []

    return FitResult(
        concat=concat,
        k_deg=k_deg,
        r_squared=r_squared,
        sd=sd,
        spep=float(spep_float),
        ci_lo=ci_lo,
        ci_hi=ci_hi,
        t=t_fit.tolist(),
        fs=fs_fit.tolist(),
        protein_id=protein_id,
        mod_sites=mod_sites,
        fs_lo=fs_lo,
        fs_hi=fs_hi,
        bio_rep=[int(b) for b in bio_rep_arr[fit_mask]],
        evidence=["mbr" if m else "q_value" for m in mbr_fit],
        metox=[bool(m) for m in metox_fit],
        dmass=dmass_fit,
        dspacing=dspacing_fit,
        fs_ds=fs_ds_fit,
        n_points=int(fit_mask.sum()),
        n_mbr=int(mbr_fit.sum()),
        n_metox=int(metox_fit.sum()),
        n_clean=int((~mbr_fit & ~metox_fit).sum()),
    )


def _null_result(concat: str, protein_id: str, mod_sites: str = "") -> FitResult:
    """Sentinel for peptides we couldn't fit — kept in the output for census."""
    return FitResult(
        concat=concat, k_deg=float("nan"), r_squared=float("nan"),
        sd=float("nan"), spep=float("nan"),
        ci_lo=float("nan"), ci_hi=float("nan"),
        t=[], fs=[], protein_id=protein_id, mod_sites=mod_sites,
        fs_lo=[], fs_hi=[], bio_rep=[],
    )


# ---------------------------------------------------------------------------
# Output DataFrame
# ---------------------------------------------------------------------------


def _build_output_df(results: list[FitResult | None]) -> pd.DataFrame:
    """DataFrame in the legacy riana_fit_peptides.txt schema + Phase F2 adds.

    ``fs_lower`` / ``fs_upper`` are the M5 per-timepoint prediction-interval
    bounds, carried here as list-cells aligned to ``t`` / ``fs`` so the wide
    file stays self-contained for the GUI curve view; the tidy one-row-per-point
    form is :func:`build_fractions_long`.
    """
    rows = [
        {
            "concat": r.concat,
            "t": r.t,
            "fs": r.fs,
            "fs_lower": r.fs_lo,
            "fs_upper": r.fs_hi,
            "evidence": r.evidence,
            "metox": r.metox,
            "dmass": r.dmass,
            "dspacing": r.dspacing,
            "fs_ds": r.fs_ds,
            "k_deg": r.k_deg,
            "R_squared": r.r_squared,
            "sd": r.sd,
            "spep": r.spep,
            "ci_lo": r.ci_lo,
            "ci_hi": r.ci_hi,
            "protein id": r.protein_id,
            "mod sites": r.mod_sites,
            "n_points": r.n_points,
            "n_mbr": r.n_mbr,
            "n_metox": r.n_metox,
            "n_clean": r.n_clean,
        }
        for r in results if r is not None
    ]
    return pd.DataFrame(rows).set_index("concat")


#: Column order for the M5 long-format per-timepoint fraction-new table.
_FRACTIONS_LONG_COLUMNS = [
    "concat", "protein id", "mod sites", "biological_replicate", "labeling_time",
    "fs", "fs_lower", "fs_upper", "fs_ds", "evidence", "metox",
]

#: Per-timepoint list-cell columns on the wide per-peptide frame. They duplicate
#: the tidy ``riana_fit_fractions.txt`` (and lose its biorep labels), so they are
#: dropped when *writing* ``riana_fit_peptides.txt`` but kept in-memory (the GUI
#: curve reads ``t``/``fs``/``evidence`` from the result frame, not the file).
_PER_TIMEPOINT_COLS = ("t", "fs", "fs_lower", "fs_upper", "evidence", "metox",
                       "dmass", "dspacing", "fs_ds")


def peptide_summary(result_df: pd.DataFrame) -> pd.DataFrame:
    """The wide per-peptide frame as a scalar summary for ``riana_fit_peptides.txt``.

    Drops the per-timepoint list-cells (:data:`_PER_TIMEPOINT_COLS`) — the
    per-timepoint detail (with biorep labels) lives in ``riana_fit_fractions.txt``.
    """
    return result_df.drop(columns=list(_PER_TIMEPOINT_COLS), errors="ignore")


def build_fractions_long(results: list[FitResult | None]) -> pd.DataFrame:
    """Explode per-peptide fits into the M5 long-format fraction-new table.

    One row per ``(concat, biological_replicate, labeling_time)`` point with the
    per-timepoint fraction-new ``fs`` and its prediction-interval bounds
    ``fs_lower`` / ``fs_upper``. This is the substrate the protein rollup
    consumes (PROJECT_REVIEW Track C); :func:`riana.core.pipeline.fit_project`
    tags it with ``experiment`` / ``condition``. Peptides that did not fit
    contribute no rows (their ``t`` list is empty).
    """
    rows = []
    # v1.1.0 item 1a — per-channel mass-defect columns (dmass_iso{k} / dspacing_iso{k},
    # mDa). Exploded here from the per-point inner lists; collected so the column set
    # is stable across peptides and omitted entirely when no peptide carried obs_mz.
    extra_cols: set[str] = set()
    for r in results:
        if r is None:
            continue
        ev = r.evidence if len(r.evidence) == len(r.t) else ["q_value"] * len(r.t)
        mx = r.metox if len(r.metox) == len(r.t) else [False] * len(r.t)
        dm = r.dmass if len(r.dmass) == len(r.t) else [None] * len(r.t)
        ds = r.dspacing if len(r.dspacing) == len(r.t) else [None] * len(r.t)
        tds = r.fs_ds if len(r.fs_ds) == len(r.t) else [float("nan")] * len(r.t)
        for ti, fsi, lo, hi, br, evi, mxi, dmi, dsi, tdi in zip(
            r.t, r.fs, r.fs_lo, r.fs_hi, r.bio_rep, ev, mx, dm, ds, tds
        ):
            row = {
                "concat": r.concat,
                "protein id": r.protein_id,
                "mod sites": r.mod_sites,
                "biological_replicate": int(br),
                "labeling_time": float(ti),
                "fs": float(fsi),
                "fs_lower": float(lo),
                "fs_upper": float(hi),
                "fs_ds": float(tdi),
                "evidence": evi,
                "metox": bool(mxi),
            }
            if dmi is not None:
                for k, v in enumerate(dmi):
                    col = f"dmass_iso{k}"
                    row[col] = float(v)
                    extra_cols.add(col)
            if dsi is not None:
                for k, v in enumerate(dsi):
                    col = f"dspacing_iso{k}"
                    row[col] = float(v)
                    extra_cols.add(col)
            rows.append(row)
    ordered_extra = sorted(
        extra_cols, key=lambda c: (c.rsplit("_iso", 1)[0], int(c.rsplit("iso", 1)[1]))
    )
    return pd.DataFrame(rows, columns=_FRACTIONS_LONG_COLUMNS + ordered_extra)
