"""DIA-NN parquet PSM intake — quantms-diann output → :class:`PSMRecord` (M6b).

The DIA counterpart of :mod:`riana.io.mztab`. quantms-diann runs DIA-NN
(≥ 2.2.0; validated against 2.5.0) and emits one ``diann_report.parquet``
covering every run, plus its own SDRF (the DIA variant — a few columns differ,
but the identity model is shared, so :mod:`riana.io.sdrf` backs both). This
reader disaggregates the parquet by its ``Run`` column and attaches the per-run
:class:`~riana.records.RunIdentity` from the SDRF, so every emitted
:class:`PSMRecord` carries its full identity — exactly like the mzTab path. The
integrator stays parser-agnostic.

**The one real difference from DDA: no MS2 scan to anchor on.** A DDA PSM names
its MS2 ``spectra_ref`` scan, and the integrator extracts the *preceding* MS1
precursor cycle. DIA has no precursor-specific MS2 scan; DIA-NN instead reports
an inferred chromatographic apex ``RT`` per precursor per run. So this reader
emits ``scan = -1`` (a "resolve me" sentinel) and carries the apex into
``PSMRecord.retention_time`` (seconds, the mzTab convention). The integrator
(:func:`riana.core.integration.resolve_rt_anchored_scans`) maps that RT to the
nearest MS1 scan in *this* mzML at integration time, then runs the same
scan-based extraction — DIA-NN is "just another ID + RT source", and Riana still
pulls the MS1 isotopologues from the mzML itself.

Field mapping (DIA-NN parquet → :class:`PSMRecord`):

- ``Run`` → ``file_name`` and the join key to the SDRF ``comment[data file]``
  stem / mzML stem; ``file_idx`` is assigned by sorted ``Run`` order.
- ``Stripped.Sequence`` → ``sequence``; ``peptide_mass`` is *recomputed* via
  :func:`riana.algorithms.mass_calc.calculate_ion_mz` so it lines up
  bit-for-bit with the mzTab / Percolator paths (Carbamidomethyl(C) always
  counted). ``Precursor.Mz`` is used only as a self-check, not the target.
- ``RT`` (minutes) → ``retention_time`` (× 60 → seconds).
- ``Precursor.Charge`` → ``charge``; ``Precursor.Mz`` → ``precursor_mz``.
- ``Q.Value`` → ``percolator_q_value`` (the field carries the q-value from
  whatever engine produced the file); ``PEP`` → ``percolator_pep``.
- ``Protein.Ids`` (``;``-joined) → ``protein_id`` (normalized to ``,`` to match
  the rest of Riana / the rollup parsimony).

**Variable modifications (M7).** Like the DDA mzTab path, DIA-NN's inline
``(UniMod:N)`` tokens are folded into the sequence as ``[UNIMOD:N]`` tokens by
:func:`_encode_peptidoform`, so a starter-set peptidoform (N-term Acetyl
``UniMod:1``, Phospho ``UniMod:21``) integrates at the *modified* m/z and the
envelope sees the mod's atoms. The fixed Carbamidomethyl(C) (``UniMod:4``) is
stripped, not tokenized — it is folded per-cysteine in ``mass_calc``. A
peptidoform carrying any mod outside the fixed + starter sets (e.g. Deamidation
``UniMod:7``) is dropped when ``drop_variable_mods`` is true (default). (Met-Ox
``UniMod:35`` *is* in the starter set, so it is kept and fit-merged onto the
unmodified curve, not dropped.)

Decoy filtering: rows with ``Decoy == 1`` are dropped by default. DIA-NN's
report is already FDR-filtered (Q.Value ≤ 0.01), so in practice every row is a
target; the guard is defensive.
"""

from __future__ import annotations

import logging
import os
import re
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from riana import constants
from riana.algorithms import mass_calc as accmass
from riana.exceptions import DataError
from riana.records import PSMRecord, RunIdentity

_LOGGER = logging.getLogger(__name__)

# An inline UniMod token in a DIA-NN Modified.Sequence, e.g. the ``(UniMod:4)``
# in ``AAC(UniMod:4)DEFK`` — positioned right after its residue (N-term mods
# lead the string), so the same string both locates and identifies the mod.
_UNIMOD_RE = re.compile(r"\(UniMod:(\d+)\)", re.IGNORECASE)

# The parquet columns this reader needs (a subset — DIA-NN emits ~70).
_COLUMNS = [
    "Run",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Precursor.Charge",
    "Precursor.Mz",
    "Decoy",
    "Q.Value",
    "PEP",
    "RT",
    "Protein.Ids",
]

# Optional phospho-proteoform columns (M7 Stage B on the DIA path). DIA-NN emits
# these only for a PTM-aware search; read them when present so an older / no-mod
# report still loads. ``Protein.Sites`` is the localized site(s) in protein
# coordinates (``[ACC:S287]`` / ``[ACC:C274,S287]`` — note it lists ALL modified
# residues incl. fixed Carbamidomethyl); ``PTM.Site.Confidence`` is the
# localization probability used to gate ambiguous calls.
_SITE_COLUMNS = ["Protein.Sites", "PTM.Site.Confidence"]

# Default localization gate: a peptidoform whose PTM.Site.Confidence is below this
# folds into the bare protein (its site is too ambiguous to define a proteoform).
# 0.75 = the common "class-I localized" cutoff (keeps ~73% of phospho rows on the
# cardiac DIA series; the 27% below are genuinely ambiguous).
_DEFAULT_MIN_SITE_CONFIDENCE = 0.75

# A DIA-NN modified residue: a letter immediately followed by ``(UniMod:N)``.
# A leading (N-terminal) mod has no preceding letter, so it does not match — it
# folds into the bare protein, matching io.mztab's ``pos < 1`` skip.
_DIANN_MOD_RE = re.compile(r"([A-Z])\(UniMod:(\d+)\)")
# A residue+position token inside Protein.Sites, e.g. ``S287`` / ``C274``.
_DIANN_SITE_RE = re.compile(r"([A-Z])(\d+)")


def _diann_proteoform_sites(modified_sequence: str, protein_sites: object,
                            confidence: object, min_confidence: float) -> str:
    """``_``-joined biological-mod proteoform suffix (e.g. ``pS287`` /
    ``pS1332_pS1333``) from DIA-NN's ``Protein.Sites``, byte-identical to the
    mzTab/DDA key format (:func:`riana.io.mztab._proteoform_sites`).

    Returns ``""`` — so the peptidoform folds into the bare protein — when it
    carries no biological mod, ``Protein.Sites`` is missing, or the localization
    ``confidence`` is below ``min_confidence``. DIA-NN lists *every* modified site
    (including the constitutive fixed Carbamidomethyl on C); we keep only the
    residues that carry a ``constants.BIOLOGICAL_MODS`` mod *in this peptidoform*
    (so the C sites drop out), with the prefix from ``constants.MOD_SITE_PREFIX``.
    """
    ps = str(protein_sites)
    if not ps or ps.lower() in ("nan", "none", "null"):
        return ""
    try:
        if float(confidence) < min_confidence:
            return ""
    except (TypeError, ValueError):
        return ""
    # residue letter -> site prefix, for the biological mods actually on this
    # peptidoform (phospho S/T/Y here; never the fixed Carbamidomethyl C).
    bio_prefix = {
        residue: constants.MOD_SITE_PREFIX.get(int(unimod), "")
        for residue, unimod in _DIANN_MOD_RE.findall(str(modified_sequence))
        if int(unimod) in constants.BIOLOGICAL_MODS
    }
    if not bio_prefix:
        return ""
    # Protein.Sites: "[ACC:S287,C274]" (a shared ';'-group -> first accession,
    # matching mzTab's first-start convention). Tokens after the ':' are the sites.
    first = ps.split(";")[0]
    inner = first[first.find(":") + 1:] if ":" in first else ""
    tags = [
        (int(pos), f"{bio_prefix[residue]}{residue}{pos}")
        for residue, pos in _DIANN_SITE_RE.findall(inner)
        if residue in bio_prefix
    ]
    return "_".join(tag for _, tag in sorted(tags))


def read_diann(
    path: str | os.PathLike[str],
    sample_map: Mapping[str, RunIdentity],
    *,
    drop_decoys: bool = True,
    drop_variable_mods: bool = True,
    min_site_confidence: float = _DEFAULT_MIN_SITE_CONFIDENCE,
) -> tuple[list[PSMRecord], dict[int, str]]:
    """Parse a DIA-NN ``report.parquet`` into typed records, keyed by the SDRF.

    Args:
        path: the quantms-diann ``diann_report.parquet``.
        sample_map: ``{data_file_stem: RunIdentity}`` from
            :func:`riana.io.sdrf.read_sdrf` (``SdrfTable.sample_map``). Each PSM
            is tagged with the identity of its ``Run`` (joined on the stem), and
            ``PSMRecord.sample`` is set to that run's ``source name``.
        drop_decoys: when true (default), rows with ``Decoy == 1`` are skipped.
        drop_variable_mods: when true (default), a peptidoform carrying a mod the
            v1 forward model can't account for (anything outside the fixed +
            starter UniMod sets, e.g. Deamidation) is dropped; starter-set mods are
            always folded into the sequence as ``[UNIMOD:N]`` tokens (see module
            docstring). With it off, such a peptidoform is kept bare.

    Returns:
        ``(records, file_index_map)`` — the PSMs and a ``{file_idx: Run}``
        mapping (file_idx assigned by sorted ``Run`` order).

    Raises:
        DataError: the parquet is missing / unreadable / empty, or references a
            ``Run`` the SDRF does not cover.
    """
    path = Path(path)
    if not path.exists():
        raise DataError(f"DIA-NN parquet not found: {path}")

    try:
        import pandas as pd
        import pyarrow.parquet as pq

        # Read the optional phospho-site columns only when the report carries them
        # (a PTM-aware search), so a no-mod / older report still loads.
        available = set(pq.ParquetFile(path).schema.names)
        read_cols = _COLUMNS + [c for c in _SITE_COLUMNS if c in available]
        df = pd.read_parquet(path, columns=read_cols)
    except ImportError as e:  # pragma: no cover - environment guard
        raise DataError(
            "reading DIA-NN parquet needs pyarrow — `pip install pyarrow` "
            "(or `pip install riana[dia]`)."
        ) from e
    except Exception as e:  # pyarrow/pandas raise a variety; surface uniformly
        raise DataError(f"failed to read DIA-NN parquet {path}: {e}") from e

    if df.empty:
        return [], {}

    # Rename the dotted DIA-NN columns to clean identifiers so row access is
    # explicit (``itertuples`` mangles ``Precursor.Charge`` → a positional
    # ``_N`` field, which silently breaks if the column order ever shifts).
    df = df.rename(columns={
        "Stripped.Sequence": "sequence",
        "Modified.Sequence": "modified_sequence",
        "Precursor.Charge": "charge",
        "Precursor.Mz": "precursor_mz",
        "Q.Value": "q_value",
        "Protein.Ids": "protein_ids",
        "Protein.Sites": "protein_sites",
        "PTM.Site.Confidence": "ptm_confidence",
    })
    has_sites = "protein_sites" in df.columns
    has_conf = "ptm_confidence" in df.columns

    if drop_decoys and "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]

    n_before = len(df)
    n_dropped = 0

    # file_idx by sorted Run order — stable and order-independent of the SDRF.
    runs = sorted(df["Run"].unique())
    file_index_map = {i: str(run) for i, run in enumerate(runs)}
    run_to_idx = {run: i for i, run in file_index_map.items()}

    unmatched = sorted({r for r in runs if r not in sample_map})
    if unmatched:
        raise DataError(
            f"DIA-NN parquet {path} has Run(s) {unmatched} not in the SDRF "
            f"(known: {sorted(sample_map)}). Check that `comment[data file]` "
            "stems match the DIA-NN Run names."
        )

    records: list[PSMRecord] = []
    proton = constants.PROTON_MASS
    mz_ppm_diffs: list[float] = []
    for row in df.itertuples(index=False):
        run = str(row.Run)
        identity = sample_map[run]
        # Fold variable mods into the sequence as [UNIMOD:N] tokens (or drop a
        # peptidoform the v1 forward model can't account for). ``drop_variable_mods``
        # gates only the drop; with it off, an unmodelable peptidoform is kept
        # bare (integrates at a partial m/z — an escape hatch, not recommended).
        encoded = _encode_peptidoform(row.modified_sequence)
        if encoded is None:
            if drop_variable_mods:
                n_dropped += 1
                continue
            sequence = str(row.sequence)
        else:
            sequence = encoded
        charge = int(row.charge)
        peptide_mass = float(
            accmass.calculate_ion_mz(sequence)
        )
        precursor_mz = _safe_float(row.precursor_mz, 0.0)
        if precursor_mz > 0 and charge > 0:
            recomputed_mz = (peptide_mass + charge * proton) / charge
            mz_ppm_diffs.append((recomputed_mz - precursor_mz) / precursor_mz * 1e6)
        rt_min = _safe_float(row.RT, 0.0)
        # M7 Stage B proteoform key on the DIA path: map DIA-NN's localized
        # ``Protein.Sites`` to the biological-mod suffix (e.g. ``pS287``), gated on
        # ``PTM.Site.Confidence``. Only for a modeled peptidoform (``encoded`` set);
        # an unmodelable one folds bare. Empty -> bare protein, exactly as before.
        mod_sites = ""
        if has_sites and encoded is not None:
            mod_sites = _diann_proteoform_sites(
                row.modified_sequence, row.protein_sites,
                row.ptm_confidence if has_conf else 1.0, min_site_confidence,
            )
        records.append(
            PSMRecord(
                scan=-1,  # DIA: no MS2 scan — resolved from RT at integrate time.
                charge=charge,
                sequence=sequence,
                peptide_mass=peptide_mass,
                sample=identity.sample,
                file_idx=run_to_idx[run],
                file_name=run,
                retention_time=rt_min * 60.0,  # minutes → seconds (mzTab convention)
                identity=identity,
                protein_id=_normalize_accessions(row.protein_ids),
                flanking_aa="",  # DIA-NN does not report flanking residues.
                precursor_mz=precursor_mz,
                neutral_mass=0.0,
                percolator_score=0.0,  # no Percolator score on the DIA path
                percolator_q_value=_safe_float(row.q_value, 1.0),
                percolator_pep=_safe_float(row.PEP, 1.0),
                distinct_matches=0,
                mod_sites=mod_sites,
            )
        )

    if n_dropped:
        _LOGGER.info(
            "io.diann: dropped %d/%d peptidoforms carrying a mod outside the M7 "
            "v1 set (kept %d); starter-set mods (N-term Acetyl, Phospho) are "
            "encoded as [UNIMOD:N] and integrated at the modified m/z.",
            n_dropped, n_before, len(records),
        )
    if mz_ppm_diffs:
        median_ppm = float(np.median(np.abs(mz_ppm_diffs)))
        if median_ppm > 50.0:
            _LOGGER.warning(
                "io.diann: recomputed peptide m/z differs from DIA-NN "
                "Precursor.Mz by a median |%.1f| ppm — a sequence/charge "
                "parsing mismatch? Expected ≪ search tolerance.", median_ppm,
            )
        else:
            _LOGGER.info(
                "io.diann: %d PSMs, recomputed m/z vs Precursor.Mz median "
                "|%.2f| ppm.", len(records), median_ppm,
            )
    return records, file_index_map


def _encode_peptidoform(modified_sequence: object) -> str | None:
    """Turn a DIA-NN ``Modified.Sequence`` into a ``[UNIMOD:N]``-tagged sequence.

    DIA-NN inlines its mods as ``(UniMod:N)`` right after the modified residue.
    Returns:

    - the bare residues when there are no mods, or only the fixed Carbamidomethyl
      (``UniMod:4``) — CAM is folded per-cysteine by ``mass_calc`` and the inline
      token is stripped so it is not double-counted;
    - the residues with ``[UNIMOD:N]`` tokens at each starter-set variable mod
      (N-term Acetyl ``1``, Phospho ``21``), in place, so it integrates at the
      modified m/z and the envelope sees the mod's atoms;
    - ``None`` to **drop** the peptidoform when it carries any mod outside
      ``FIXED_UNIMODS`` ∪ ``STARTER_VARIABLE_UNIMODS``. Mirrors
      :func:`riana.io.mztab._encode_peptidoform`.
    """
    s = str(modified_sequence)
    ids = [int(u) for u in _UNIMOD_RE.findall(s)]
    if any(u not in constants.FIXED_UNIMODS | constants.STARTER_VARIABLE_UNIMODS
           for u in ids):
        return None

    def _replace(match: re.Match) -> str:
        u = int(match.group(1))
        return "" if u in constants.FIXED_UNIMODS else f"[UNIMOD:{u}]"

    return _UNIMOD_RE.sub(_replace, s)


def _normalize_accessions(protein_ids: object) -> str:
    """``;``-joined DIA-NN ``Protein.Ids`` → Riana's ``,``-joined convention."""
    s = "" if protein_ids is None else str(protein_ids)
    if s.lower() in ("nan", "none"):
        return ""
    return ",".join(p for p in re.split(r"[;,]", s) if p.strip())


def _safe_float(value: object, default: float) -> float:
    """Coerce *value* to float, returning *default* on ``None`` / ``NaN``."""
    if value is None:
        return default
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    if f != f:  # NaN
        return default
    return f
