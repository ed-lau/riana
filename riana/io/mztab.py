"""mzTab PSM intake — quantms / OpenMS output → :class:`PSMRecord`.

The **primary** ID path as of M6a. quantms produces *one* mzTab covering every
run in the experiment; :func:`read_mztab` disaggregates it by ``ms_run`` and
attaches the per-run :class:`~riana.records.RunIdentity` from the SDRF
(:mod:`riana.io.sdrf`), so every emitted :class:`PSMRecord` carries its full
identity. The (now demoted) single-mzML :mod:`riana.io.percolator` path emits
the same record type, so the integrator stays parser-agnostic.

Field mapping (quantms PSM section → :class:`PSMRecord`):

- ``spectra_ref`` ``ms_run[N]:... scan=NNNN`` →
  ``file_idx = N-1`` (0-indexed by ms_run order), ``file_name`` from
  ``MTD ms_run[N]-location``, ``scan = NNNN``.
- ``sequence`` → ``sequence`` (bare amino-acid string). ``peptide_mass`` is
  *recomputed* via :func:`riana.accmass.calculate_ion_mz` so it lines up
  bit-for-bit with the Percolator path's recompute — Carbamidomethyl(C) is
  always counted (it is the only fixed mod the calibration set uses).

**Variable modifications (M7 caveat).** The envelope is keyed on the *stripped*
``sequence``, so a variable-mod peptidoform (Oxidation(M) ``UNIMOD:35``,
Phospho(S) ``UNIMOD:21``, N-term Acetyl ``UNIMOD:1``, …) would integrate at the
**unmodified** m/z. Until M7 threads variable mods into the envelope,
``drop_variable_mods`` (default true) filters those PSM rows out via the mzTab
``modifications`` column (``pos-UNIMOD:id``, comma-joined) rather than
mis-target them. The fixed Carbamidomethyl(C) (``UNIMOD:4``) is kept — it is
folded into the recomputed peptide mass. This mirrors the DIA-NN path's
``drop_variable_mods`` (:mod:`riana.io.diann`), so a variable-mod search
(e.g. ``timeseries_lve_atr``) integrates only its unmodified peptidoforms on
both the DDA and DIA surfaces.
- ``opt_global_q-value`` → ``percolator_q_value`` (named ``percolator_*``
  on :class:`PSMRecord` for compatibility — the field carries the q-value
  from whatever ID engine produced the file).
- ``opt_global_Posterior_Error_Probability_score`` → ``percolator_pep``.
- ``search_engine_score[1]`` → ``percolator_score``.
- ``retention_time`` → ``retention_time`` (seconds, mzTab convention) — the DIA
  RT-apex prior (M6b); the DDA integration path uses the MS2 scan, not this.

Decoy filtering: rows with ``opt_global_cv_MS:1002217_decoy_peptide == 1``
are dropped by default. quantms emits target+decoy PSMs together, so dropping
decoys here matches the target-only Percolator surface.

**Run identity (M6a).** Pass ``sample_map`` (``{mzML stem: RunIdentity}`` from
:func:`riana.io.sdrf.read_sdrf`) to attach the full per-run identity, keyed by
the ``ms_run[N]-location`` basename. The legacy single-label call (a bare
``sample=`` string, no SDRF) is kept for the synthetic/round-trip tests.
"""

from __future__ import annotations

import logging
import os
import re
from collections.abc import Mapping
from pathlib import Path

from pyteomics import mztab

from riana.algorithms import mass_calc as accmass
from riana.exceptions import DataError
from riana.records import PSMRecord, RunIdentity

_LOGGER = logging.getLogger(__name__)


_MS_RUN_RE = re.compile(r"ms_run\[(\d+)\]")
_SPECTRA_REF_RE = re.compile(r"ms_run\[(\d+)\]:.*?scan=(\d+)")
_LOCATION_RE = re.compile(r"^(?:file://)?(.*)$")
# A UniMod token inside an mzTab ``modifications`` cell, e.g. ``7-UNIMOD:4``
# (``pos-UNIMOD:id``, comma-joined for multiple sites). Matched case-insensitively
# so it also accepts the DIA-NN-style ``UniMod:`` spelling.
_UNIMOD_RE = re.compile(r"UNIMOD:(\d+)", re.IGNORECASE)
# Carbamidomethyl — the fixed cysteine mod the calibration/turnover searches use;
# folded into the recomputed peptide mass, so it is NOT a "variable" peptidoform.
_FIXED_UNIMODS = frozenset({"4"})


def read_mztab(
    path: str | os.PathLike[str],
    sample_map: Mapping[str, RunIdentity] | None = None,
    *,
    sample: str | None = None,
    drop_decoys: bool = True,
    drop_variable_mods: bool = True,
) -> tuple[list[PSMRecord], dict[int, str]]:
    """Parse a quantms-style mzTab into typed records.

    Args:
        path: ``*.mzTab`` from quantms / OpenMS.
        sample_map: ``{mzML stem: RunIdentity}`` from
            :func:`riana.io.sdrf.read_sdrf`. When given, each PSM is tagged with
            the identity of its ``ms_run`` (joined on the location basename), and
            ``PSMRecord.sample`` is set to that run's ``source name``. Mutually
            exclusive with ``sample``.
        sample: legacy single-label fallback (no SDRF) — written into every
            record's ``sample`` field with ``identity=None``. Mutually exclusive
            with ``sample_map``.
        drop_decoys: when true (default), rows flagged as decoys are skipped.
        drop_variable_mods: when true (default), PSM rows carrying a non-fixed
            UniMod (Oxidation, Phospho, N-term Acetyl, …) in the ``modifications``
            column are dropped — they would otherwise integrate at the unmodified
            m/z until M7 (see module docstring). Carbamidomethyl (``UNIMOD:4``)
            is fixed and kept.

    Returns:
        ``(records, file_index_map)`` — the PSMs and a ``{file_idx: file_name}``
        mapping (file_idx is 0-based ms_run index, file_name is the basename
        of the mzML location with its extension stripped).

    Raises:
        DataError: neither/both of ``sample_map`` / ``sample`` given, or (with a
            ``sample_map``) an ``ms_run`` location that the SDRF does not cover.
    """
    if (sample_map is None) == (sample is None):
        raise DataError(
            "read_mztab needs exactly one of `sample_map` (SDRF identity) or "
            "`sample` (legacy single-label)."
        )
    path = Path(path)
    if not path.exists():
        raise DataError(f"mzTab not found: {path}")

    try:
        table = mztab.MzTab(str(path))
    except Exception as e:  # pyteomics raises a variety; surface uniformly
        raise DataError(f"failed to parse mzTab {path}: {e}") from e

    file_index_map = _build_file_index_map(table.metadata, path)
    if not file_index_map:
        raise DataError(f"no ms_run entries in mzTab metadata: {path}")

    if sample_map is not None:
        unmatched = sorted(
            {stem for stem in file_index_map.values() if stem not in sample_map}
        )
        if unmatched:
            raise DataError(
                f"mzTab {path} references ms_run mzML(s) {unmatched} that are not "
                f"in the SDRF (known: {sorted(sample_map)}). Check that "
                "`comment[data file]` matches the mzML names."
            )

    psm_df = table.spectrum_match_table
    if psm_df is None or len(psm_df) == 0:
        return [], file_index_map

    decoy_col = "opt_global_cv_MS:1002217_decoy_peptide"
    if drop_decoys and decoy_col in psm_df.columns:
        psm_df = psm_df[psm_df[decoy_col] == 0]

    if drop_variable_mods and "modifications" in psm_df.columns:
        n_before = len(psm_df)
        psm_df = psm_df[psm_df["modifications"].map(_only_fixed_mods)]
        n_dropped = n_before - len(psm_df)
        if n_dropped:
            _LOGGER.info(
                "io.mztab: dropped %d/%d variable-mod peptidoforms "
                "(integrate at the unmodified m/z until M7); keeping %d.",
                n_dropped, n_before, len(psm_df),
            )

    records: list[PSMRecord] = []
    for row in psm_df.to_dict(orient="records"):
        ms_run, scan = _parse_spectra_ref(row["spectra_ref"])
        file_idx = ms_run - 1
        file_name = file_index_map.get(file_idx, "")
        identity = sample_map.get(file_name) if sample_map is not None else None
        record_sample = identity.sample if identity is not None else (sample or "")
        sequence = str(row["sequence"])
        records.append(
            PSMRecord(
                scan=scan,
                charge=int(row["charge"]),
                sequence=sequence,
                peptide_mass=float(accmass.calculate_ion_mz(sequence)),
                sample=record_sample,
                file_idx=file_idx,
                file_name=file_name,
                retention_time=_safe_float(row.get("retention_time"), 0.0),
                identity=identity,
                protein_id=str(row.get("accession", "") or ""),
                # mzTab spells pre/post as comma-joined residues across all
                # protein matches; collapse to the first for parity with the
                # Percolator ``flanking aa`` (which is a 2-char string).
                flanking_aa=_first_flanker(row.get("pre"), row.get("post")),
                precursor_mz=_safe_float(row.get("exp_mass_to_charge"), 0.0),
                # neutral mass = (m/z - proton) * charge, but the mzTab carries
                # exp_mass_to_charge directly; leave neutral_mass=0.0 so a
                # missing value is honest rather than fabricated.
                neutral_mass=0.0,
                percolator_score=_safe_float(row.get("search_engine_score[1]"), 1.0),
                percolator_q_value=_safe_float(row.get("opt_global_q-value"), 1.0),
                percolator_pep=_safe_float(
                    row.get("opt_global_Posterior_Error_Probability_score"), 1.0
                ),
                distinct_matches=0,
            )
        )
    return records, file_index_map


def _build_file_index_map(metadata: dict, source_path: Path) -> dict[int, str]:
    """Build ``{file_idx (0-based): file_name}`` from ``MTD ms_run[N]-location``.

    File names are stripped of any ``file://`` scheme and ``.mzML`` extension
    so they line up with the basenames the Percolator path emits.
    """
    out: dict[int, str] = {}
    for key, value in metadata.items():
        if not key.endswith("-location"):
            continue
        m = _MS_RUN_RE.match(key)
        if m is None:
            continue
        ms_run = int(m.group(1))
        loc = str(value)
        loc_m = _LOCATION_RE.match(loc)
        location = loc_m.group(1) if loc_m else loc
        basename = Path(os.path.basename(location)).stem
        out[ms_run - 1] = basename
    return out


def _only_fixed_mods(modifications: object) -> bool:
    """True when every UniMod in an mzTab ``modifications`` cell is a fixed mod.

    The cell is ``pos-UNIMOD:id`` tokens (comma-joined), or ``None`` / ``"null"``
    for an unmodified PSM (which trivially qualifies). Carbamidomethyl
    (``UNIMOD:4``) is fixed and folded into the recomputed peptide mass, so it is
    kept; any other UniMod (Oxidation, Phospho, N-term Acetyl, …) marks a
    variable peptidoform dropped until M7. Mirrors
    :func:`riana.io.diann._only_fixed_mods`.
    """
    if modifications is None:
        return True
    s = str(modifications)
    if s.lower() in ("nan", "none", "null", ""):
        return True
    return all(t in _FIXED_UNIMODS for t in _UNIMOD_RE.findall(s))


def _parse_spectra_ref(value: str) -> tuple[int, int]:
    m = _SPECTRA_REF_RE.search(str(value))
    if m is None:
        raise DataError(
            f"could not parse spectra_ref '{value}' "
            "(expected 'ms_run[N]:... scan=NNNN')"
        )
    return int(m.group(1)), int(m.group(2))


def _safe_float(value: object, default: float) -> float:
    """Coerce *value* to float, returning *default* on ``None`` / ``NaN``.

    Critically, ``0.0`` is preserved — a previous version used
    ``float(x or default)`` which mis-handled q-value=0.0 (the most confident
    PSMs) as the default sentinel. Pyteomics returns ``None`` for genuinely
    missing entries and ``float('nan')`` for unparseable ones; we treat both
    as missing.
    """
    if value is None:
        return default
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    if f != f:  # NaN check without numpy
        return default
    return f


def _first_flanker(pre: object, post: object) -> str:
    """Collapse comma-joined pre/post residues to a single 2-char string."""
    pre_s = "" if pre is None else str(pre).split(",")[0]
    post_s = "" if post is None else str(post).split(",")[0]
    if pre_s == "null":
        pre_s = ""
    if post_s == "null":
        post_s = ""
    return (pre_s + post_s)[:2]
