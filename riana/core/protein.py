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

**Proteoform keys (M7 Stage B).** When the input carries a ``mod sites`` column
(a per-peptidoform biological-mod suffix in protein coordinates, e.g. ``pS34476``,
set by the IO reader for phospho — see :func:`riana.io.mztab._proteoform_sites`),
the resolved ``protein`` becomes ``accession_<suffix>`` so a phosphopeptidoform
rolls up as its own turnover unit (``A2ASS6_pS34476``) instead of collapsing into
the bare protein. Two peptides covering the same site share a key (the site is
protein-coordinate, not peptide-relative). Unmodified peptidoforms and
constitutive/artifactual mods (N-term Acetyl, etc. → empty suffix) stay on the
bare accession. Inputs without the column behave exactly as before.

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
from typing import Callable, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from riana.core import models
from riana.exceptions import DataError
from riana.progress import iter_progress
from riana.records import GROUP_KEY_COLUMNS, PROTEIN_KEY_COLUMNS

_LOGGER = logging.getLogger(__name__)

_MODELS = {
    "simple": models.one_exponent,
    "guan": models.two_compartment_guan,
    "fornasiero": models.two_compartment_fornasiero,
}
#: The linearized cross-sample model — not a ``curve_fit`` ODE like the others, but
#: an OLS in φ = log(1−θ) space that fits a protein's conditions jointly to yield a
#: Δk test (:mod:`riana.core.linear_model`). A ``model`` choice, mutually exclusive
#: with ``_MODELS``.
LINEAR_MODEL = "linear simple"
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
#: The rollup grouping key — ``(experiment, condition, protein)``, defined once
#: in :mod:`riana.records` (:data:`PROTEIN_KEY_COLUMNS`) and shared with the fit
#: recombination so the stages can't drift. (Distinct from ``PROTEIN_COLUMNS``
#: below, which is the *output* schema.)
_GROUP_KEYS = list(PROTEIN_KEY_COLUMNS)

#: Output column order for ``riana_protein.txt``. One ``k_deg`` (+ CI / R²) from
#: the selected ``method``; ``peptide_median_k`` is the near-free median of the
#: peptides' fitted k, carried for comparison regardless of method.
#: ``n_mbr`` / ``n_metox`` / ``n_clean`` are the **raw** (peptide × biorep ×
#: timepoint) fraction-point composition — how many measurements came from MBR
#: transfers / Met-Ox-merged forms / neither — so they are *not* on the same scale
#: as ``n_points`` (the collapsed per-timepoint refit points); they answer "how
#: many of this protein's data points are clean".
PROTEIN_COLUMNS = [
    *PROTEIN_KEY_COLUMNS, "method",
    "n_peptides", "n_replicates", "n_timepoints", "n_points",
    "n_mbr", "n_metox", "n_clean",
    "k_deg", "ci_lo", "ci_hi", "R_squared", "peptide_median_k",
]

#: Output schema for ``model="linear simple"``: the per-condition protein rows
#: (k_deg = −slope of φ=log(1−θ), its CI, joint R²) plus the protein-level
#: cross-condition Δk test (``delta_k`` = k(other) − k(reference), its SE, the
#: contrast p, and the Benjamini-Hochberg p across proteins).
PROTEIN_LINEAR_COLUMNS = [
    *PROTEIN_KEY_COLUMNS, "method", "n_peptides", "n_points",
    "n_mbr", "n_metox", "n_clean",
    "k_deg", "ci_lo", "ci_hi", "R_squared", "peptide_median_k",
    "delta_k", "delta_k_se", "delta_k_p", "delta_k_p_adj",
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
    min_spep: int | None = None,
    min_r2: float | None = None,
    k_cv_max: float = 0.2,
    rescue_r2: float = 0.6,
    min_fit_points: int | None = None,
    workers: int = 1,
    n_boot: int = 200,
    boot_ci_pct: tuple[float, float] = (5.0, 95.0),
    random_state: int = 1337,
    phi_limit: float = -4.0,
    reference_condition: str | None = None,
    test_condition: str | None = None,
    exclude_mbr: bool = False,
    progress_callback: "Callable[[int, int], None] | None" = None,
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
        min_spep: optional Spep (labelling-site) admission gate applied *before*
            rollup — drop peptides below it (defense-in-depth; the primary gate is
            at fit, so a manifest rollup already inherits it). Off (``None``) by
            default. See report 2026-06-28_spep_curation_gate.
        min_r2: optional peptide R² admission gate applied *before* rollup. When
            ``None`` (default) no gate is applied — the inverse-variance weighting
            already down-weights noisy peptides; pass a value (e.g. 0.8) to also
            hard-exclude peptides whose kinetic fit doesn't follow the model, for
            an A/B against the unfiltered result. A peptide is kept if
            ``R² ≥ min_r2`` OR — the flat-curve rescue, so well-measured but
            slow / low-dynamic-range peptides (low R² only because θ barely moves)
            survive — if ``R² ≥ rescue_r2 and k_cv < k_cv_max``.
        k_cv_max / rescue_r2: the relative-uncertainty rescue thresholds (only
            used when ``min_r2`` is set). ``k_cv`` = ``(ci_hi − ci_lo)/(2·|k|)``,
            the rate constant's scale-free CV. ``k_cv_max ≤ 0`` disables the
            rescue; the ``rescue_r2`` floor guards against degenerate k≈0 rail-hits
            (see :func:`_r2_admitted`).
        min_fit_points: peptide-level biological-replicate gate — keep only
            peptidoforms fit on ≥ this many points. ``None`` (default) auto-resolves
            to **2 for a single-timepoint experiment** (detected from one distinct
            labeling time) and **off otherwise**. For single-timepoint data R² is
            degenerate, so curation is this replicate floor + the ``k_cv`` gate
            (``k_cv_max``), with the R²/rescue machinery bypassed.
        method: ``"weighted"`` (default, the inverse-variance per-timepoint
            collapse) or ``"pooled"`` (all peptide×timepoint points, no collapse;
            pseudoreplication-naive).
        workers: worker *processes* for the per-protein refit — the parallelism
            lever for the GIL-bound refit (``curve_fit`` × bootstrap). Each protein
            gets an independent RNG stream seeded from ``random_state``, so the
            result is **identical** regardless of ``workers`` (and of completion
            order).
        n_boot / boot_ci_pct / random_state: bootstrap CI controls (ignored by
            ``model="linear simple"``, whose CIs are analytic).
        phi_limit / reference_condition / test_condition: only for
            ``model="linear simple"`` — the plateau-truncation threshold in φ-space
            (default −4 ≈ θ 0.98) and the Δk contrast conditions. ``reference`` is
            the baseline (default the alphabetically-first); naming ``test`` as well
            contrasts that specific pair even when >2 conditions are present (the
            multi-group interim — see :func:`riana.core.linear_model.fit_linear_deltak`).

    Returns:
        One row per ``(experiment, condition, protein)`` — a ``method`` tag, the
        selected estimator's ``k_deg`` / ``ci_lo`` / ``ci_hi`` / ``R_squared``
        (``NaN`` where it could not fit), ``n_peptides`` / ``n_points``, and the
        comparison ``peptide_median_k``. ``result.attrs["protein_points"]`` maps
        each protein to the ``(t_list, fs_list)`` its refit used (GUI curve).
        ``model="linear simple"`` instead returns the :data:`PROTEIN_LINEAR_COLUMNS`
        schema (per-condition φ-slope k + the protein-level Δk test).
    """
    if model not in _MODELS and model != LINEAR_MODEL:
        raise DataError(
            f"unknown kinetic model {model!r}; expected one of "
            f"{sorted(_MODELS) + [LINEAR_MODEL]}")
    if parsimony not in _PARSIMONY:
        raise DataError(
            f"parsimony must be one of {list(_PARSIMONY)}, got {parsimony!r}")
    if method not in _METHODS:
        raise DataError(
            f"method must be one of {list(_METHODS)}, got {method!r}")
    kk = dict(a_0=0.0, a_max=1.0, **dict(kinetic_kwargs or {}))

    peptides = _ensure_group_cols(peptides)
    fractions = _ensure_group_cols(fractions)
    if "concat" not in peptides.columns:
        raise DataError("peptides input is missing the 'concat' column.")
    # Decide attribution once over the peptide map, then apply to both frames.
    mapping = _resolve_parsimony(peptides, parsimony)
    peptides = _apply_parsimony(peptides, mapping)
    fractions = _apply_parsimony(fractions, mapping)

    # --exclude-mbr: drop match-between-runs data points before the rollup refit
    # (MBR points are used by default). No-op when the column is absent.
    if exclude_mbr and "evidence" in fractions.columns:
        fractions = fractions[fractions["evidence"] != "mbr"].copy()

    # Optional Spep admission gate (off by default; the primary gate runs at fit
    # so a manifest-driven rollup inherits it — this is the explicit-input /
    # belt-and-braces lever). Spep is sequence-derived, so a concat's value is the
    # same across conditions; filter on it.
    if min_spep is not None and min_spep > 0:
        if "spep" not in peptides.columns:
            raise DataError("peptides input is missing 'spep' (needed for --min-spep).")
        keep = peptides.loc[peptides["spep"] >= min_spep, "concat"].unique()
        peptides = peptides[peptides["concat"].isin(keep)].copy()
        fractions = fractions[fractions["concat"].isin(keep)].copy()

    # Single-timepoint detection (data-driven, needs no flag): an experiment with one
    # labeling timepoint has no kinetic curve, so R² is degenerate (≈0/NaN) and the
    # R²/rescue-floor gate does not apply. Curate on the rate-constant relative
    # uncertainty (k_cv) plus a biological-replicate floor (min_fit_points) instead.
    single_tp = (
        "labeling_time" in fractions.columns
        and int(fractions["labeling_time"].nunique(dropna=True)) <= 1
    )
    if min_fit_points is None:
        min_fit_points = 2 if single_tp else 0
    if single_tp:
        _LOGGER.info(
            "rollup: single labeling timepoint detected — R² is not applicable, "
            "curating on >=%d replicate fit points and k_cv < %s "
            "(tune with --min-fit-points / --k-cv; --min-fit-points 1 --k-cv 0 = off).",
            max(min_fit_points, 1), k_cv_max,
        )

    # Peptide-level replicate gate: keep peptidoforms fit on >= min_fit_points points
    # (distinct (biorep, timepoint) fit points) — the dominant curation lever for
    # single-timepoint data, where n_points IS the biological-replicate count.
    if min_fit_points and min_fit_points > 1:
        if "n_points" not in peptides.columns:
            raise DataError(
                "peptides input is missing 'n_points' (needed for --min-fit-points).")
        keep = peptides.loc[peptides["n_points"] >= min_fit_points, "concat"].unique()
        peptides = peptides[peptides["concat"].isin(keep)].copy()
        fractions = fractions[fractions["concat"].isin(keep)].copy()

    # Curation gate. Single-timepoint: k_cv only (R² bypassed). Multi-timepoint: the
    # R² gate + its flat-curve k_cv rescue, off by default (min_r2 is None).
    if single_tp:
        if k_cv_max is not None and k_cv_max > 0.0:
            admitted = _k_cv_admitted(peptides, k_cv_max)
            peptides = peptides[peptides["concat"].isin(admitted)].copy()
            fractions = fractions[fractions["concat"].isin(admitted)].copy()
    elif min_r2 is not None:
        admitted = _r2_admitted(peptides, min_r2, k_cv_max, rescue_r2)
        peptides = peptides[peptides["concat"].isin(admitted)].copy()
        fractions = fractions[fractions["concat"].isin(admitted)].copy()

    stats = _peptide_stats(peptides, min_peptides=min_peptides)
    # Per-protein data-point census from the (filtered) fraction points: how many
    # are MBR transfers / Met-Ox-merged / clean — so the protein k carries the same
    # composition columns as the per-peptide fit output. Merged into stats so both
    # the kinetic and linear paths emit them.
    stats = stats.merge(_breakdown_stats(fractions), on=_GROUP_KEYS, how="left")
    for col in ("n_mbr", "n_metox", "n_clean"):
        stats[col] = stats[col].fillna(0).astype(int)

    if model == LINEAR_MODEL:
        return _rollup_linear(
            stats, fractions, method=method, min_peptides=min_peptides,
            min_points=min_points, phi_limit=phi_limit,
            reference_condition=reference_condition,
            test_condition=test_condition,
            progress_callback=progress_callback)

    model_fn = _MODELS[model]
    refit, points = _refit_table(
        fractions, model_fn=model_fn, kinetic_kwargs=kk, method=method,
        min_peptides=min_peptides, min_points=min_points,
        n_boot=n_boot, boot_ci_pct=boot_ci_pct, random_state=random_state,
        workers=workers, progress_callback=progress_callback,
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


def _rollup_linear(
    stats: pd.DataFrame,
    fractions: pd.DataFrame,
    *,
    method: str,
    min_peptides: int,
    min_points: int,
    phi_limit: float,
    reference_condition: str | None,
    test_condition: str | None,
    progress_callback: "Callable[[int, int], None] | None" = None,
) -> pd.DataFrame:
    """The ``model="linear simple"`` path — φ-space OLS + cross-condition Δk.

    Collapses peptides to the same ``(condition, t, θ)`` points the weighted/pooled
    refit uses, then fits each protein's conditions **jointly** in φ = log(1−θ)
    space (:func:`riana.core.linear_model.fit_linear_deltak`) for a per-condition
    k and a Δk test. Returns the :data:`PROTEIN_LINEAR_COLUMNS` schema and attaches
    the collapsed points for the GUI/fractions output (plotted in φ-space).
    """
    from riana.core.linear_model import fit_linear_deltak

    long = _collapse_long(fractions, method=method, min_peptides=min_peptides)
    lin = fit_linear_deltak(
        long, phi_limit=phi_limit, min_points=min_points,
        reference_condition=reference_condition,
        test_condition=test_condition,
        progress_callback=progress_callback)

    out = pd.merge(stats, lin, on=_GROUP_KEYS, how="right")
    out["method"] = LINEAR_MODEL
    for col in PROTEIN_LINEAR_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    result = (
        out[PROTEIN_LINEAR_COLUMNS]
        .sort_values(_GROUP_KEYS)
        .reset_index(drop=True)
    )
    points = {
        keys: (g["labeling_time"].tolist(), g["theta"].tolist())
        for keys, g in long.groupby(_GROUP_KEYS, sort=False)
    }
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
        columns=[*PROTEIN_KEY_COLUMNS, "labeling_time", "fs"],
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

    # M7 Stage B: a per-concat biological-mod proteoform suffix (e.g. ``pS34476``)
    # appended to the resolved accession, so a phosphopeptidoform rolls up as its
    # own unit. Absent column / empty value → bare accession (the pre-B behaviour).
    if "mod sites" in peptides.columns:
        sites = peptides[["concat", "mod sites"]].drop_duplicates("concat")
        site_map = dict(zip(sites["concat"], sites["mod sites"].fillna("")))
    else:
        site_map = {}

    def key(concat: str, accession: str) -> str:
        suffix = str(site_map.get(concat, "") or "")
        return f"{accession}_{suffix}" if suffix else accession

    if parsimony == "unique":
        keep = pep["accs"].map(len) == 1
        out = pep.loc[keep, ["concat"]].copy()
        accs = pep.loc[keep, "accs"].map(lambda a: a[0])
        out["protein"] = [key(c, a) for c, a in zip(out["concat"], accs)]
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
            rows.append((concat, key(concat, _collapsed(accs))))
            continue
        if len({_base(a) for a in accs}) > 1:         # multiple genes → reject
            continue
        if any(_is_isoform(a) and a in isoforms_with_unique for a in accs):
            continue                                  # an isoform is real → reject
        rows.append((concat, key(concat, _collapsed(accs))))  # fold into canonical
    return pd.DataFrame(rows, columns=["concat", "protein"])


def _apply_parsimony(df: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """Filter *df* to the attributable peptides and attach the ``protein`` label."""
    if "concat" not in df.columns:
        raise DataError("input is missing the 'concat' column.")
    return df.drop(columns=["protein"], errors="ignore").merge(
        mapping, on="concat", how="inner"
    )


def _ensure_group_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Guarantee the group key columns exist (the legacy path has neither)."""
    df = df.copy()
    for c in GROUP_KEY_COLUMNS:
        df[c] = df[c].fillna("") if c in df.columns else ""
    return df


def _k_cv_admitted(peptides: pd.DataFrame, k_cv_max: float) -> set:
    """Concats whose rate-constant relative uncertainty ``k_cv < k_cv_max`` — the
    **single-timepoint** curation gate. At one labeling timepoint R² is degenerate,
    so it is bypassed and admission rides on ``k_cv = (ci_hi − ci_lo) / (2·|k|)``
    alone (matches :func:`riana.core.fitting._k_cv`). A single-point fit has ``k_cv``
    NaN (undefined uncertainty — no replication) and a k≈0 rail-hit gives NaN/inf, so
    both fail the ``<`` and are excluded, as intended.
    """
    need = {"concat", "k_deg", "ci_lo", "ci_hi"}
    missing = need - set(peptides.columns)
    if missing:
        raise DataError(f"single-timepoint k_cv gate needs columns {sorted(missing)}.")
    k = peptides["k_deg"].to_numpy(dtype=float)
    lo = peptides["ci_lo"].to_numpy(dtype=float)
    hi = peptides["ci_hi"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        k_cv = (hi - lo) / (2.0 * np.abs(k))
    return set(peptides.loc[k_cv < k_cv_max, "concat"])


def _r2_admitted(
    peptides: pd.DataFrame,
    min_r2: float,
    k_cv_max: float,
    rescue_r2: float,
) -> set:
    """Concats passing the R² admission gate (with a relative-uncertainty rescue).

    Keep a peptide if ``R² ≥ min_r2`` (the primary goodness-of-fit gate), OR — to
    rescue well-measured peptides whose R² is low only because the curve is flat
    (slow turnover / low dynamic range: the pathology Lau *Nat Commun* 2018 flags)
    — if ``R² ≥ rescue_r2 AND k_cv < k_cv_max``, where ``k_cv`` is the rate
    constant's relative uncertainty ``(ci_hi − ci_lo) / (2·|k|)`` (a scale-free CV
    of k̂; see :func:`riana.core.fitting._k_cv`). The ``rescue_r2`` floor is
    load-bearing: without it the rescue admits degenerate k≈0 rail-hits whose
    bootstrap CI collapses to a spuriously tight ``k_cv≈0`` (R² deeply negative),
    which on noisier data wrecks protein-level ranking (reports 2026-06-26 /
    2026-07-01 — any floor > the negative-R² band suffices; 0.6 is the validated
    default). ``k_cv_max ≤ 0`` disables the rescue (R²-only gate). Non-converged
    peptides (``NaN`` k/R²/CI) fail every comparison and are excluded, as intended.
    """
    rescue = k_cv_max is not None and k_cv_max > 0.0
    need = {"concat", "R_squared"} | (
        {"k_deg", "ci_lo", "ci_hi"} if rescue else set())
    missing = need - set(peptides.columns)
    if missing:
        raise DataError(
            f"--min-r2 needs columns {sorted(missing)} in the peptides input.")
    r2 = peptides["R_squared"].to_numpy(dtype=float)
    keep = r2 >= min_r2
    if rescue:
        k = peptides["k_deg"].to_numpy(dtype=float)
        lo = peptides["ci_lo"].to_numpy(dtype=float)
        hi = peptides["ci_hi"].to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            k_cv = (hi - lo) / (2.0 * np.abs(k))   # matches fitting._k_cv; NaN/inf at k≈0 fail <
        keep = keep | ((r2 >= rescue_r2) & (k_cv < k_cv_max))
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


def _breakdown_stats(fractions: pd.DataFrame) -> pd.DataFrame:
    """Per-protein fraction-point census: n_mbr / n_metox / n_clean.

    Counts the raw ``(concat, biorep, timepoint)`` points by their per-point
    ``evidence`` (``q_value`` | ``mbr``) and ``metox`` flag (both on
    ``riana_fit_fractions.txt``). A point that is both MBR and Met-Ox lands in
    n_mbr and n_metox but not n_clean (matching the per-peptide fit breakdown).
    Absent columns (pre-MBR / pre-M7 runs) default to clean.
    """
    f = fractions
    is_mbr = (f["evidence"].astype(str) == "mbr"
              if "evidence" in f.columns else pd.Series(False, index=f.index))
    if "metox" in f.columns and f["metox"].dtype == bool:
        is_metox = f["metox"]
    elif "metox" in f.columns:
        is_metox = f["metox"].astype(str).str.strip().str.lower().isin(("true", "1"))
    else:
        is_metox = pd.Series(False, index=f.index)
    g = f.assign(
        n_mbr=is_mbr.astype(int),
        n_metox=is_metox.astype(int),
        n_clean=(~is_mbr & ~is_metox).astype(int),
    ).groupby(_GROUP_KEYS, sort=False)[["n_mbr", "n_metox", "n_clean"]].sum()
    return g.reset_index()


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
    workers: int = 1,
    progress_callback: "Callable[[int, int], None] | None" = None,
) -> tuple[pd.DataFrame, dict]:
    """Returns ``(refit table, points)`` where ``points`` maps
    ``(experiment, condition, protein)`` → the ``(t_list, fs_list)`` the refit
    was fit on — the substrate for the GUI's per-protein curve.

    ``method="weighted"`` collapses peptides within each (biorep, timepoint) by
    inverse-variance before fitting; ``method="pooled"`` fits all peptide×
    timepoint points directly (pseudoreplication). Each protein is an independent
    work unit; the per-protein ``curve_fit`` × bootstrap is GIL-bound, so
    ``workers`` (process-level) is the parallelism lever. Per-group RNG streams
    keep the result identical regardless of worker count and completion order.
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

    refit_one = partial(
        _refit_one_group, model_fn=model_fn, kinetic_kwargs=kinetic_kwargs,
        method=method, min_points=min_points, n_boot=n_boot,
        boot_ci_pct=boot_ci_pct, random_state=random_state,
    )

    if workers > 1 and len(groups) > 1:
        # Process-level parallelism — the real lever for the GIL-bound refit
        # (curve_fit + per-protein residual bootstrap). The group frames are
        # pickled to each worker ONCE via the initializer; only an integer index
        # crosses per task. The per-group RNG makes the result worker-count-
        # independent. Mirrors the fit -W pattern (core/fitting.py).
        init_args = (groups, model_fn, kinetic_kwargs, method, min_points,
                     n_boot, boot_ci_pct, random_state)
        chunk = max(1, len(groups) // (workers * 8))
        with futures.ProcessPoolExecutor(
            max_workers=workers,
            initializer=_init_refit_worker,
            initargs=init_args,
        ) as ex:
            computed = list(iter_progress(
                ex.map(_refit_one_group_worker, range(len(groups)), chunksize=chunk),
                len(groups), progress_callback))
    else:
        # Serial: the per-protein curve_fit + bootstrap is GIL-bound, so threading
        # gave no speedup; `workers` (processes) is the only parallelism lever.
        computed = list(iter_progress(
            (refit_one(g) for g in groups), len(groups), progress_callback))

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


def _collapse_group(grp: pd.DataFrame, method: str) -> tuple[list, list]:
    """Collapse one protein group's peptide θ into ``(t_list, fs_list)`` points.

    ``method="pooled"`` keeps every peptide×timepoint θ (pseudoreplication);
    ``method="weighted"`` collapses peptides within each (biorep, timepoint) by
    inverse variance. Shared by the nonlinear refit (:func:`_refit_one_group`) and
    the linear collapse (:func:`_collapse_long`) so both see identical points.
    """
    if method == "pooled":
        fs = grp["fs"].to_numpy(dtype=float)
        t = grp["labeling_time"].to_numpy(dtype=float)
        keep = np.isfinite(fs) & np.isfinite(t)
        return t[keep].tolist(), fs[keep].tolist()
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
    return t_list, fs_list


def _collapse_long(
    fractions: pd.DataFrame, *, method: str, min_peptides: int
) -> pd.DataFrame:
    """Long table of collapsed ``(experiment, condition, protein, labeling_time,
    theta)`` points — the substrate the linear (φ-space) model fits, built with
    the *same* collapse the weighted/pooled refit uses (:func:`_collapse_group`).
    """
    rows = []
    for (exp, cond, prot), grp in fractions.groupby(_GROUP_KEYS, sort=False):
        if grp["concat"].nunique() < min_peptides:
            continue
        t_list, fs_list = _collapse_group(grp, method)
        for t, fs in zip(t_list, fs_list):
            rows.append({"experiment": exp, "condition": cond, "protein": prot,
                         "labeling_time": float(t), "theta": float(fs)})
    return pd.DataFrame(
        rows, columns=["experiment", "condition", "protein",
                       "labeling_time", "theta"])


def _refit_one_group(
    item,
    *,
    model_fn,
    kinetic_kwargs: dict,
    method: str,
    min_points: int,
    n_boot: int,
    boot_ci_pct: tuple[float, float],
    random_state: int,
):
    """Refit one ``(experiment, condition, protein)`` group → ``(keys, row, pts)``.

    Module-level (not a closure) so it pickles for the ``ProcessPoolExecutor``
    on spawn-start platforms (macOS). The per-group RNG (``_group_rng``) makes the
    result independent of worker/thread count and completion order.
    """
    keys, grp = item
    n_rep = int(grp["biological_replicate"].nunique())
    n_tp = int(grp["labeling_time"].nunique())
    t_list, fs_list = _collapse_group(grp, method)
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


# --- process-pool plumbing (mirrors core/fitting._init_fit_worker) -----------
# Each worker holds the group list + shared refit config in a module global, set
# ONCE by the pool initializer, so only an integer index crosses the boundary per
# task (not every group frame, per task).
_REFIT_WORKER_STATE: dict[str, object] = {}


def _init_refit_worker(
    groups, model_fn, kinetic_kwargs, method, min_points,
    n_boot, boot_ci_pct, random_state,
) -> None:
    _REFIT_WORKER_STATE.update(
        groups=groups, model_fn=model_fn, kinetic_kwargs=kinetic_kwargs,
        method=method, min_points=min_points, n_boot=n_boot,
        boot_ci_pct=boot_ci_pct, random_state=random_state,
    )


def _refit_one_group_worker(index: int):
    s = _REFIT_WORKER_STATE
    return _refit_one_group(
        s["groups"][index], model_fn=s["model_fn"],
        kinetic_kwargs=s["kinetic_kwargs"], method=s["method"],
        min_points=s["min_points"], n_boot=s["n_boot"],
        boot_ci_pct=s["boot_ci_pct"], random_state=s["random_state"],
    )


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
