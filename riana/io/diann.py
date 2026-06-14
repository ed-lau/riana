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

**Variable modifications (M7 caveat).** Like the DDA mzTab path, the envelope is
keyed on the *stripped* sequence, so a variable-mod peptidoform (e.g.
Oxidation(M), ``UniMod:35``) would integrate at the **unmodified** m/z. Until M7
threads variable mods into the envelope, ``drop_variable_mods`` (default true)
filters those peptidoforms out rather than mis-target them. The fixed
Carbamidomethyl(C) (``UniMod:4``) is kept — it is folded into the recomputed
peptide mass. On the validated cardiac DIA set this drops ~1.8% of rows and
costs only the few precursors observed *exclusively* as an oxidized form.

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

# A UniMod token inside a DIA-NN Modified.Sequence, e.g. ``(UniMod:35)``.
_UNIMOD_RE = re.compile(r"\(UniMod:(\d+)\)")
# Carbamidomethyl — the fixed cysteine mod the calibration/turnover searches use;
# folded into the recomputed peptide mass, so it is NOT a "variable" peptidoform.
_FIXED_UNIMODS = frozenset({"4"})

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


def read_diann(
    path: str | os.PathLike[str],
    sample_map: Mapping[str, RunIdentity],
    *,
    drop_decoys: bool = True,
    drop_variable_mods: bool = True,
) -> tuple[list[PSMRecord], dict[int, str]]:
    """Parse a DIA-NN ``report.parquet`` into typed records, keyed by the SDRF.

    Args:
        path: the quantms-diann ``diann_report.parquet``.
        sample_map: ``{data_file_stem: RunIdentity}`` from
            :func:`riana.io.sdrf.read_sdrf` (``SdrfTable.sample_map``). Each PSM
            is tagged with the identity of its ``Run`` (joined on the stem), and
            ``PSMRecord.sample`` is set to that run's ``source name``.
        drop_decoys: when true (default), rows with ``Decoy == 1`` are skipped.
        drop_variable_mods: when true (default), peptidoforms carrying a
            non-fixed UniMod (e.g. Oxidation) are dropped — they would otherwise
            integrate at the unmodified m/z until M7 (see module docstring).

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

        df = pd.read_parquet(path, columns=_COLUMNS)
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
    })

    if drop_decoys and "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]

    n_before = len(df)
    if drop_variable_mods:
        df = df[df["modified_sequence"].map(_only_fixed_mods)]
        n_dropped = n_before - len(df)
        if n_dropped:
            _LOGGER.info(
                "io.diann: dropped %d/%d variable-mod peptidoforms "
                "(integrate at the unmodified m/z until M7); keeping %d.",
                n_dropped, n_before, len(df),
            )

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
        sequence = str(row.sequence)
        charge = int(row.charge)
        peptide_mass = float(
            accmass.calculate_ion_mz(sequence)
        )
        precursor_mz = _safe_float(row.precursor_mz, 0.0)
        if precursor_mz > 0 and charge > 0:
            recomputed_mz = (peptide_mass + charge * proton) / charge
            mz_ppm_diffs.append((recomputed_mz - precursor_mz) / precursor_mz * 1e6)
        rt_min = _safe_float(row.RT, 0.0)
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
            )
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


def _only_fixed_mods(modified_sequence: object) -> bool:
    """True when every UniMod in *modified_sequence* is a fixed mod we keep.

    A bare (unmodified) sequence trivially qualifies. Carbamidomethyl (UniMod:4)
    is fixed and folded into the recomputed peptide mass, so it is kept; any
    other UniMod (Oxidation, N-term acetyl, ...) marks a variable peptidoform we
    drop until M7.
    """
    toks = _UNIMOD_RE.findall(str(modified_sequence))
    return all(t in _FIXED_UNIMODS for t in toks)


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
