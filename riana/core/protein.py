# -*- coding: utf-8 -*-

"""Protein-level turnover rollup (Track C).

Rolls the `riana fit` per-peptide outputs up to one turnover estimate per
protein, grouped by ``(experiment, condition, protein)`` so side-by-side groups
stay distinct. One ``k_deg`` per protein from the selected ``method``:

- **``"weighted"`` (default)** — the pseudoreplication-safe refit and the reason
  M5 was built. Within each ``(protein, biological_replicate, labeling_time)``
  the peptides' fraction-new θ are collapsed by an **inverse-variance weighted
  average** (σ from the M5 prediction interval, ``σ ≈ (fs_upper − fs_lower) /
  3.29`` for the 5–95 band), then **one** protein ``k_deg`` is fit to the
  collapsed ``(t, θ)`` points *across timepoints AND bioreps* — different
  biological replicates stay independent points, giving the refit honest
  degrees of freedom.
- **``"pooled"``** — fit one ``k_deg`` to **all** peptide×timepoint θ points
  with no collapse. The pseudoreplication-naive stance (peptides treated as
  independent replicates ⇒ over-confident CI); kept for comparison, **not
  recommended**.

CI via the residual bootstrap :func:`riana.core.fitting.fit_run` uses. A
near-free **``peptide_median_k``** column (median of the peptides' fitted k) is
carried for comparison regardless of method — point estimators over peptide k
are a trivial ``groupby`` the user can also do on the peptide file. (Harmonic
mean was dropped: DIY *and* outlier-sensitive on the low-k tail. The linearized
``log(1−θ) = −kt`` fit + a cross-sample Δk test are a deferred milestone.)

**Parsimony is a summarize-time decision** (a shared peptide's envelope blends
both proteins' turnover, so its signature can't be attributed):

- ``"unique"`` (default) — keep only single-accession peptides; a peptide whose
  ``protein id`` is multi-accession (comma/semicolon-joined) is dropped. The
  attributed ``protein`` is the bare UniProt accession (isoform suffix kept when
  the sole accession *is* an isoform).
- ``"isoform"`` — fold isoform-only-shared peptides into the canonical entry,
  but only when no isoform in the group carries its own unique peptide
  (a dataset-wide pass over the peptide↔protein map). Adapted from
  ``02_R_parsimony_reference.Rmd`` (Juber/Lau) with fixes: Riana's ``,``
  separator (not ``;``), a base-accession fallback where the R left
  ``collapsed_uniprot`` NA (an all-isoform group with no canonical present), and
  the standard UniProt ``-N`` suffix. (The R weights by ``log2(Int)``; we weight
  the θ collapse by the M5 inverse variance instead — see ``_weighted_theta``.)

Consumes ``riana_fit_peptides.txt`` (per-peptide ``k_deg`` + ``protein id``) and
the M5 ``riana_fit_fractions.txt`` (``concat, biological_replicate,
labeling_time, fs, fs_lower, fs_upper``). The wide and long frames carry
``experiment`` / ``condition`` when produced by the manifest fit path; the
legacy single-curve path has neither, so both default to ``""``.
"""

from __future__ import annotations

import hashlib
import logging
import re
from concurrent import futures
from functools import partial
from typing import Mapping

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from riana.core import models
from riana.exceptions import DataError

_LOGGER = logging.getLogger(__name__)

_MODELS = {
    "simple": models.one_exponent,
    "guan": models.two_compartment_guan,
    "fornasiero": models.two_compartment_fornasiero,
}
_K_DEG_INIT = 0.5
_K_DEG_BOUNDS = ([1e-4], [10.0])

#: A normal's 5–95 percentile span is ``2 · 1.6449 · σ``; invert to recover σ
#: from the M5 prediction-interval width.
_PI_SPAN_SIGMA = 2.0 * 1.6448536269514722  # ≈ 3.2897

_PARSIMONY = ("unique", "isoform")
#: Rollup estimators. ``weighted`` (default) = the biorep-aware per-timepoint
#: inverse-variance collapse then refit (pseudoreplication-safe). ``pooled`` =
#: fit one k to *all* peptide×timepoint θ points with no collapse — the
#: pseudoreplication-naive stance (kept for comparison, not recommended).
_METHODS = ("weighted", "pooled")
_GROUP_KEYS = ["experiment", "condition", "protein"]

#: Output column order for ``riana_protein.txt``. One ``k_deg`` (+ CI / R²) from
#: the selected ``method``; ``peptide_median_k`` is the near-free median of the
#: peptides' fitted k, carried for comparison regardless of method.
PROTEIN_COLUMNS = [
    "experiment", "condition", "protein", "method",
    "n_peptides", "n_replicates", "n_timepoints", "n_points",
    "k_deg", "ci_lo", "ci_hi", "R_squared", "peptide_median_k",
]


# --------------------------------------------------------------------------- #
# public entry point
# --------------------------------------------------------------------------- #
def rollup_proteins(
    peptides: pd.DataFrame,
    fractions: pd.DataFrame,
    *,
    model: str = "simple",
    kinetic_kwargs: Mapping[str, float] | None = None,
    parsimony: str = "unique",
    method: str = "weighted",
    min_peptides: int = 2,
    min_points: int = 3,
    min_r2: float | None = None,
    alt_k: float = 0.025,
    alt_se: float = 0.05,
    alt_r2: float = 0.0,
    threads: int = 1,
    n_boot: int = 200,
    boot_ci_pct: tuple[float, float] = (5.0, 95.0),
    random_state: int = 1337,
) -> pd.DataFrame:
    """Roll per-peptide fits up to one ``k_deg`` per ``(experiment, condition,
    protein)`` via the median and the biorep-aware weighted refit.

    Args:
        peptides: ``riana_fit_peptides.txt`` content — needs ``protein id`` and
            ``k_deg``.
        fractions: ``riana_fit_fractions.txt`` content — needs ``concat``,
            ``protein id``, ``biological_replicate``, ``labeling_time``, ``fs``,
            ``fs_lower``, ``fs_upper``.
        model: kinetic model for the refit (must match how peptides were fit).
        kinetic_kwargs: ``k_p`` / ``k_r`` / ``r_p`` for guan/fornasiero
            (ignored by ``simple``); ``a_0`` / ``a_max`` are pinned 0 / 1.
        parsimony: ``"unique"`` (default) or ``"isoform"`` (deferred).
        min_peptides: a protein needs at least this many attributed peptides to
            be reported.
        min_points: the refit needs at least this many collapsed ``(t, θ)``
            points.
        min_r2: optional peptide R² admission gate applied *before* rollup. When
            ``None`` (default) no gate is applied — the inverse-variance weighting
            already down-weights noisy peptides; pass a value (e.g. 0.8) to also
            hard-exclude peptides whose kinetic fit doesn't follow the model, for
            an A/B against the unfiltered result. A peptide is kept if
            ``R² ≥ min_r2`` OR — the JCI slow-turnover admit, so legitimately
            slow peptides (low R² only because θ barely moves) survive — if
            ``k ≤ alt_k and SE ≤ alt_se and R² ≥ alt_r2``.
        alt_k / alt_se / alt_r2: the slow-turnover admit thresholds (only used
            when ``min_r2`` is set). ``SE`` is the fit's ``sd`` (bootstrap k_deg
            std).
        method: ``"weighted"`` (default, the inverse-variance per-timepoint
            collapse) or ``"pooled"`` (all peptide×timepoint points, no collapse;
            pseudoreplication-naive).
        threads: worker threads for the per-protein refit (the expensive
            ``curve_fit`` × bootstrap). Each protein gets an independent RNG
            stream seeded from ``random_state``, so the result is **identical**
            regardless of ``threads`` (and of completion order).
        n_boot / boot_ci_pct / random_state: bootstrap CI controls.

    Returns:
        One row per ``(experiment, condition, protein)`` — a ``method`` tag, the
        selected estimator's ``k_deg`` / ``ci_lo`` / ``ci_hi`` / ``R_squared``
        (``NaN`` where it could not fit), ``n_peptides`` / ``n_points``, and the
        comparison ``peptide_median_k``. ``result.attrs["protein_points"]`` maps
        each protein to the ``(t_list, fs_list)`` its refit used (GUI curve).
    """
    if model not in _MODELS:
        raise DataError(
            f"unknown kinetic model {model!r}; expected one of {sorted(_MODELS)}")
    if parsimony not in _PARSIMONY:
        raise DataError(
            f"parsimony must be one of {list(_PARSIMONY)}, got {parsimony!r}")
    if method not in _METHODS:
        raise DataError(
            f"method must be one of {list(_METHODS)}, got {method!r}")
    model_fn = _MODELS[model]
    kk = dict(a_0=0.0, a_max=1.0, **dict(kinetic_kwargs or {}))

    peptides = _ensure_group_cols(peptides)
    fractions = _ensure_group_cols(fractions)
    if "concat" not in peptides.columns:
        raise DataError("peptides input is missing the 'concat' column.")
    # Decide attribution once over the peptide map, then apply to both frames.
    mapping = _resolve_parsimony(peptides, parsimony)
    peptides = _apply_parsimony(peptides, mapping)
    fractions = _apply_parsimony(fractions, mapping)

    # Optional peptide R² admission gate (off by default).
    if min_r2 is not None:
        admitted = _r2_admitted(peptides, min_r2, alt_k, alt_se, alt_r2)
        peptides = peptides[peptides["concat"].isin(admitted)].copy()
        fractions = fractions[fractions["concat"].isin(admitted)].copy()

    stats = _peptide_stats(peptides, min_peptides=min_peptides)
    refit, points = _refit_table(
        fractions, model_fn=model_fn, kinetic_kwargs=kk, method=method,
        min_peptides=min_peptides, min_points=min_points,
        n_boot=n_boot, boot_ci_pct=boot_ci_pct, random_state=random_state,
        threads=threads,
    )

    out = pd.merge(stats, refit, on=_GROUP_KEYS, how="outer")
    out["method"] = method
    for col in PROTEIN_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    result = (
        out[PROTEIN_COLUMNS]
        .sort_values(_GROUP_KEYS)
        .reset_index(drop=True)
    )
    # The collapsed (t, θ) points behind each refit, for the GUI curve view.
    # Keyed by (experiment, condition, protein).
    result.attrs["protein_points"] = points
    return result


def build_rollup_fractions(result: pd.DataFrame) -> pd.DataFrame:
    """Long/tidy table of the collapsed ``(t, θ)`` points behind each protein
    refit — the inverse-variance-weighted fraction-new the GUI curve plots,
    one row per ``(experiment, condition, protein, labeling_time)``. Built from
    ``result.attrs["protein_points"]`` (empty when none were attached).
    """
    rows = []
    for (exp, cond, prot), pts in result.attrs.get("protein_points", {}).items():
        t_list, fs_list = pts
        for t, fs in zip(t_list, fs_list):
            rows.append({
                "experiment": exp, "condition": cond, "protein": prot,
                "labeling_time": float(t), "fs": float(fs),
            })
    return pd.DataFrame(
        rows,
        columns=["experiment", "condition", "protein", "labeling_time", "fs"],
    )


# --------------------------------------------------------------------------- #
# parsimony
# --------------------------------------------------------------------------- #
#: ``sp|P12345-2|NAME`` / ``tr|ACC|NAME`` → the accession (group 1).
_SP_TR_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|")
#: An isoform suffix: Swiss-Prot ``P12345-2`` and JCAST ``P12345-J1`` → strip
#: the trailing ``-N`` / ``-JN``. A bare trailing ``-`` (no digits) is NOT a
#: suffix, so ``\d+`` is required.
_ISOFORM_SUFFIX_RE = re.compile(r"-J?\d+$")
#: Riana joins multi-protein ids with ``,``; FragPipe/MSFragger uses ``;``.
_ACC_SEP_RE = re.compile(r"[;,]")


def _accession(token: str) -> str:
    """``sp|P12345-2|CDV3_HUMAN`` → ``P12345-2``; bare tokens pass through."""
    token = token.strip()
    m = _SP_TR_RE.match(token)
    return m.group(1) if m else token


def _accessions(protein_id: object) -> list[str]:
    """Split a (``,``/``;``-joined) ``protein id`` into UniProt accessions."""
    return [_accession(t) for t in _ACC_SEP_RE.split(str(protein_id)) if t.strip()]


def _base(acc: str) -> str:
    """Drop the isoform suffix: ``P12345-2`` → ``P12345``."""
    return _ISOFORM_SUFFIX_RE.sub("", acc)


def _is_isoform(acc: str) -> bool:
    return bool(_ISOFORM_SUFFIX_RE.search(acc))


def _collapsed(accs: list[str]) -> str:
    """The protein this accession set rolls up to: the first canonical (no ``-N``,
    sorted) accession, else — when the group is all isoforms with no canonical
    present — the shared base accession (the R reference left this NA)."""
    canonical = sorted(a for a in accs if not _is_isoform(a))
    if canonical:
        return canonical[0]
    if len(accs) == 1:
        return accs[0]                       # unique-to-an-isoform: keep it
    return sorted({_base(a) for a in accs})[0]


def _resolve_parsimony(peptides: pd.DataFrame, parsimony: str) -> pd.DataFrame:
    """Decide, once over the whole peptide↔protein map, which peptides are
    attributable and to which protein.

    Returns ``DataFrame[concat, protein]`` of the kept peptides. Built from the
    peptide table (the authoritative map) and then applied to *both* the peptide
    and fraction frames, so the ``isoform`` rule's dataset-wide "does any isoform
    carry a unique peptide" evidence is consistent across them.
    """
    if "protein id" not in peptides.columns:
        raise DataError("peptides input is missing the 'protein id' column.")
    pep = peptides[["concat", "protein id"]].drop_duplicates("concat").copy()
    pep["accs"] = pep["protein id"].map(_accessions)

    if parsimony == "unique":
        keep = pep["accs"].map(len) == 1
        out = pep.loc[keep, ["concat"]].copy()
        out["protein"] = pep.loc[keep, "accs"].map(lambda a: a[0])
        return out.reset_index(drop=True)

    # parsimony == "isoform": dataset-wide isoform evidence.
    # Isoform accessions that are the SOLE accession of some peptide (i.e. that
    # isoform has its own unique peptide somewhere in the data).
    isoforms_with_unique = {
        a[0] for a in pep.loc[pep["accs"].map(len) == 1, "accs"]
        if _is_isoform(a[0])
    }
    rows = []
    for concat, accs in zip(pep["concat"], pep["accs"]):
        if len(accs) == 1:                            # single accession → keep
            rows.append((concat, _collapsed(accs)))
            continue
        if len({_base(a) for a in accs}) > 1:         # multiple genes → reject
            continue
        if any(_is_isoform(a) and a in isoforms_with_unique for a in accs):
            continue                                  # an isoform is real → reject
        rows.append((concat, _collapsed(accs)))       # fold into canonical
    return pd.DataFrame(rows, columns=["concat", "protein"])


def _apply_parsimony(df: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """Filter *df* to the attributable peptides and attach the ``protein`` label."""
    if "concat" not in df.columns:
        raise DataError("input is missing the 'concat' column.")
    return df.drop(columns=["protein"], errors="ignore").merge(
        mapping, on="concat", how="inner"
    )


def _ensure_group_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Guarantee ``experiment`` / ``condition`` exist (legacy path has neither)."""
    df = df.copy()
    for c in ("experiment", "condition"):
        df[c] = df[c].fillna("") if c in df.columns else ""
    return df


def _r2_admitted(
    peptides: pd.DataFrame,
    min_r2: float,
    alt_k: float,
    alt_se: float,
    alt_r2: float,
) -> set:
    """Concats passing the R² admission gate (with a JCI slow-turnover admit).

    Keep a peptide if ``R² ≥ min_r2``, OR — to retain legitimately slow-turnover
    peptides whose R² is low only because θ barely moves — if
    ``k ≤ alt_k and SE ≤ alt_se and R² ≥ alt_r2`` (``SE`` = the fit ``sd``
    column). Non-converged peptides (``NaN`` k/R²/sd) fail every comparison and
    are excluded, which is the intent.
    """
    need = {"concat", "R_squared", "k_deg", "sd"}
    missing = need - set(peptides.columns)
    if missing:
        raise DataError(
            f"--min-r2 needs columns {sorted(missing)} in the peptides input.")
    r2 = peptides["R_squared"].to_numpy(dtype=float)
    k = peptides["k_deg"].to_numpy(dtype=float)
    se = peptides["sd"].to_numpy(dtype=float)
    keep = (r2 >= min_r2) | ((k <= alt_k) & (se <= alt_se) & (r2 >= alt_r2))
    return set(peptides.loc[keep, "concat"])


def _group_rng(random_state: int, key) -> np.random.Generator:
    """An independent RNG stream for one protein group, stable across runs/threads.

    Seeds ``default_rng`` from ``random_state`` + a content hash of the group key
    (``hashlib``, not the salted built-in ``hash``), so each ``(experiment,
    condition, protein)`` bootstraps from a fixed, order-independent stream — the
    rollup result is identical whether it ran on 1 thread or N.
    """
    digest = hashlib.blake2b(repr(key).encode(), digest_size=8).digest()
    return np.random.default_rng(
        [int(random_state), int.from_bytes(digest, "little")])


# --------------------------------------------------------------------------- #
# per-protein peptide stats (n_peptides + the comparison median-k)
# --------------------------------------------------------------------------- #
def _peptide_stats(peptides: pd.DataFrame, *, min_peptides: int) -> pd.DataFrame:
    """Per-protein peptide count + the near-free ``peptide_median_k``.

    ``n_peptides`` is the distinct attributed peptides; ``peptide_median_k`` is
    the median of their fitted ``k_deg`` — a trivial comparison column (not the
    headline estimate; that comes from the refit ``method``).
    """
    if "k_deg" not in peptides.columns:
        raise DataError("peptides input is missing the 'k_deg' column.")
    rows = []
    for keys, grp in peptides.groupby(_GROUP_KEYS, sort=False):
        n_pep = grp["concat"].nunique()
        if n_pep < min_peptides:
            continue
        ks = grp["k_deg"].to_numpy(dtype=float)
        ks = ks[np.isfinite(ks)]
        exp, cond, prot = keys
        rows.append({
            "experiment": exp, "condition": cond, "protein": prot,
            "n_peptides": int(n_pep),
            "peptide_median_k": float(np.median(ks)) if len(ks) else float("nan"),
        })
    return pd.DataFrame(
        rows, columns=_GROUP_KEYS + ["n_peptides", "peptide_median_k"])


# --------------------------------------------------------------------------- #
# the protein refit (method = weighted collapse | pooled points)
# --------------------------------------------------------------------------- #
def _refit_table(
    fractions: pd.DataFrame,
    *,
    model_fn,
    kinetic_kwargs: dict,
    method: str,
    min_peptides: int,
    min_points: int,
    n_boot: int,
    boot_ci_pct: tuple[float, float],
    random_state: int,
    threads: int = 1,
) -> tuple[pd.DataFrame, dict]:
    """Returns ``(refit table, points)`` where ``points`` maps
    ``(experiment, condition, protein)`` → the ``(t_list, fs_list)`` the refit
    was fit on — the substrate for the GUI's per-protein curve.

    ``method="weighted"`` collapses peptides within each (biorep, timepoint) by
    inverse-variance before fitting; ``method="pooled"`` fits all peptide×
    timepoint points directly (pseudoreplication). Each protein is an
    independent work unit dispatched over ``threads`` workers; per-group RNG
    streams keep the result identical regardless of thread count.
    """
    need = {"concat", "biological_replicate", "labeling_time",
            "fs", "fs_lower", "fs_upper"}
    missing = need - set(fractions.columns)
    if missing:
        raise DataError(f"fractions input is missing columns {sorted(missing)}.")

    groups = [
        (keys, grp) for keys, grp in fractions.groupby(_GROUP_KEYS, sort=False)
        if grp["concat"].nunique() >= min_peptides
    ]

    def _one(item):
        keys, grp = item
        n_rep = int(grp["biological_replicate"].nunique())
        n_tp = int(grp["labeling_time"].nunique())
        if method == "pooled":
            # All peptide×timepoint θ points, no collapse (pseudoreplication).
            fs = grp["fs"].to_numpy(dtype=float)
            t = grp["labeling_time"].to_numpy(dtype=float)
            keep = np.isfinite(fs) & np.isfinite(t)
            t_list, fs_list = t[keep].tolist(), fs[keep].tolist()
        else:  # weighted: collapse peptides within each (biorep, timepoint).
            t_list, fs_list = [], []
            for (_br, t), cell in grp.groupby(
                ["biological_replicate", "labeling_time"], sort=False
            ):
                theta = _weighted_theta(
                    cell["fs"].to_numpy(), cell["fs_lower"].to_numpy(),
                    cell["fs_upper"].to_numpy(),
                )
                if np.isfinite(theta):
                    t_list.append(float(t))
                    fs_list.append(theta)
        if len(t_list) < min_points:
            return keys, None, None
        fit = _fit_kdeg(
            np.asarray(t_list), np.asarray(fs_list),
            model_fn=model_fn, kinetic_kwargs=kinetic_kwargs,
            n_boot=n_boot, boot_ci_pct=boot_ci_pct,
            rng=_group_rng(random_state, keys),
        )
        if fit is None:
            return keys, None, None
        exp, cond, prot = keys
        k, r2, lo, hi = fit
        row = {
            "experiment": exp, "condition": cond, "protein": prot,
            "n_replicates": n_rep, "n_timepoints": n_tp,
            "n_points": int(len(t_list)), "k_deg": k,
            "ci_lo": lo, "ci_hi": hi, "R_squared": r2,
        }
        return keys, row, (t_list, fs_list)

    if threads <= 1 or len(groups) <= 1:
        computed = [_one(g) for g in groups]
    else:
        with futures.ThreadPoolExecutor(max_workers=threads) as ex:
            computed = list(ex.map(_one, groups))

    rows = []
    points: dict = {}
    for keys, row, pts in computed:
        if row is None:
            continue
        rows.append(row)
        points[keys] = pts
    table = pd.DataFrame(
        rows, columns=_GROUP_KEYS + [
            "n_replicates", "n_timepoints", "n_points",
            "k_deg", "ci_lo", "ci_hi", "R_squared"]
    )
    return table, points


def _weighted_theta(
    fs: np.ndarray, lo: np.ndarray, hi: np.ndarray
) -> float:
    """Inverse-variance weighted mean of θ within one (protein, biorep, t) cell.

    σ comes from the M5 prediction-interval width (``(hi − lo) / 3.29``).
    Peptides with a missing/zero σ are filled with the cell's median valid σ
    (counted as a typical peptide rather than dropped); if no peptide in the
    cell has a usable σ, the mean is unweighted.
    """
    fs = np.asarray(fs, dtype=float)
    keep = np.isfinite(fs)
    if not keep.any():
        return float("nan")
    sigma = (np.asarray(hi, float) - np.asarray(lo, float)) / _PI_SPAN_SIGMA
    valid = np.isfinite(sigma) & (sigma > 0)
    if valid.any():
        sigma = np.where(valid, sigma, float(np.median(sigma[valid])))
        w = 1.0 / np.square(sigma)
    else:
        w = np.ones_like(fs)
    w = np.where(keep, w, 0.0)
    return float(np.sum(w * np.where(keep, fs, 0.0)) / np.sum(w))


def _fit_kdeg(
    t: np.ndarray,
    fs: np.ndarray,
    *,
    model_fn,
    kinetic_kwargs: dict,
    n_boot: int,
    boot_ci_pct: tuple[float, float],
    rng: np.random.Generator,
) -> tuple[float, float, float, float] | None:
    """Fit one ``k_deg`` to collapsed ``(t, θ)`` + a residual-bootstrap CI."""
    try:
        popt, _ = curve_fit(
            partial(model_fn, **kinetic_kwargs), t, fs,
            bounds=_K_DEG_BOUNDS, p0=[_K_DEG_INIT], maxfev=2000)
    except (RuntimeError, ValueError):
        return None
    k = float(popt[0])
    pred = np.asarray(model_fn(t, k_deg=k, **kinetic_kwargs), dtype=float)
    resid = fs - pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((fs - fs.mean()) ** 2))
    r2 = float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot
    n = len(t)
    boot: list[float] = []
    for _ in range(n_boot):
        fs_star = pred + resid[rng.integers(0, n, n)]
        try:
            pb, _ = curve_fit(
                partial(model_fn, **kinetic_kwargs), t, fs_star,
                bounds=_K_DEG_BOUNDS, p0=[k], maxfev=2000)
            boot.append(float(pb[0]))
        except (RuntimeError, ValueError):
            continue
    if len(boot) >= 10:
        lo, hi = (float(p) for p in np.percentile(boot, boot_ci_pct))
    else:
        lo = hi = float("nan")
    return k, r2, lo, hi
