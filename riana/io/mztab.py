"""mzTab PSM intake — quantms / OpenMS output → :class:`PSMRecord`.

The Week 2 alternative to :mod:`riana.io.percolator`. quantms produces *one*
mzTab covering every fraction in the experiment; the Percolator path through
snakemake produces *one file per fraction*. Both adapters emit the same
:class:`PSMRecord` type so the Week 3 integrator is parser-agnostic.

Field mapping (quantms PSM section → :class:`PSMRecord`):

- ``spectra_ref`` ``ms_run[N]:... scan=NNNN`` →
  ``file_idx = N-1`` (0-indexed by ms_run order), ``file_name`` from
  ``MTD ms_run[N]-location``, ``scan = NNNN``.
- ``sequence`` → ``sequence`` (bare amino-acid string). ``peptide_mass`` is
  *recomputed* via :func:`riana.accmass.calculate_ion_mz` so it lines up
  bit-for-bit with the Percolator path's recompute — Carbamidomethyl(C) is
  always counted (it is the only fixed mod the calibration set uses;
  variable mods are not yet handled — peptidoform parsing lands when
  variable-mod calibration data appears).
- ``opt_global_q-value`` → ``percolator_q_value`` (named ``percolator_*``
  on :class:`PSMRecord` for compatibility — the field carries the q-value
  from whatever ID engine produced the file).
- ``opt_global_Posterior_Error_Probability_score`` → ``percolator_pep``.
- ``search_engine_score[1]`` → ``percolator_score``.

Decoy filtering: rows with ``opt_global_cv_MS:1002217_decoy_peptide == 1``
are dropped by default. quantms emits target+decoy PSMs together; the
Percolator path's snakemake pipeline already wrote out target-only files,
so dropping decoys here aligns the two surfaces.

**Multi-fraction TODO.** ``sample`` is a single string applied to every
emitted record because today every calibration mzML is a single,
unfractionated run keyed 1:1 to a time point. When a real multi-fraction
quantms run lands, this should grow a ``sample_map: dict[ms_run_idx, str]``
parameter sourced from the SDRF samplesheet (``ms_run → sample``); the
ms_run-keyed PSM emission already supports that. Deferred until the first
multi-fraction dataset.
"""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from pathlib import Path

from pyteomics import mztab

from riana.algorithms import mass_calc as accmass
from riana.exceptions import DataError
from riana.records import PSMRecord


_MS_RUN_RE = re.compile(r"ms_run\[(\d+)\]")
_SPECTRA_REF_RE = re.compile(r"ms_run\[(\d+)\]:.*?scan=(\d+)")
_LOCATION_RE = re.compile(r"^(?:file://)?(.*)$")


def read_mztab(
    path: str | os.PathLike[str],
    sample: str,
    ignored_mods: Sequence[str] = (),
    *,
    drop_decoys: bool = True,
) -> tuple[list[PSMRecord], dict[int, str]]:
    """Parse a quantms-style mzTab into typed records.

    Args:
        path: ``*.mzTab`` from quantms / OpenMS.
        sample: sample label written into every emitted record.
        ignored_mods: forwarded to :func:`accmass.calculate_ion_mz` when the
            peptide_mass is recomputed.
        drop_decoys: when true (default), rows flagged as decoys are skipped.

    Returns:
        ``(records, file_index_map)`` — the PSMs and a ``{file_idx: file_name}``
        mapping (file_idx is 0-based ms_run index, file_name is the basename
        of the mzML location with the ``.mzML`` extension stripped). The
        caller uses this map to bridge ms_run identity to the snakemake
        Percolator file_idx scheme (which is per-fraction).
    """
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

    psm_df = table.spectrum_match_table
    if psm_df is None or len(psm_df) == 0:
        return [], file_index_map

    decoy_col = "opt_global_cv_MS:1002217_decoy_peptide"
    if drop_decoys and decoy_col in psm_df.columns:
        psm_df = psm_df[psm_df[decoy_col] == 0]

    records: list[PSMRecord] = []
    for row in psm_df.to_dict(orient="records"):
        ms_run, scan = _parse_spectra_ref(row["spectra_ref"])
        file_idx = ms_run - 1
        file_name = file_index_map.get(file_idx, "")
        sequence = str(row["sequence"])
        records.append(
            PSMRecord(
                scan=scan,
                charge=int(row["charge"]),
                sequence=sequence,
                peptide_mass=float(
                    accmass.calculate_ion_mz(sequence, ignored_mods=ignored_mods)
                ),
                sample=sample,
                file_idx=file_idx,
                file_name=file_name,
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
