"""Percolator PSM intake — emits :class:`riana.records.PSMRecord`.

Ports the header-sniff Crux-vs-standalone logic from
:class:`riana.peptides.ReadPercolator`. Two changes vs. 0.9.0:

- Function-style API. The result is a plain ``list[PSMRecord]`` (no carrier
  class, no module-level state) so the consumer composes its own per-fraction
  grouping/filtering with stdlib comprehensions.
- The Crux/standalone dispatch is by explicit header inspection (the 0.9.0
  exception-as-control-flow path predates the header check that landed during
  M1 stabilization; standalone is still supported because MSFragger output
  uses it).

The mzTab adapter (:mod:`riana.io.mztab`) emits the same :class:`PSMRecord`
type so downstream integration is parser-agnostic — that is the whole point
of the Week 2 I/O layer.

**Multi-fraction TODO.** ``sample`` is applied to every emitted record as a
single string. The snakemake calibration runs are 1 fraction per sample so
this is fine today, but Riana's multi-fraction support already lives at the
``file_idx`` level — when a real multi-fraction Percolator output appears,
the right move is a ``sample_map: dict[file_idx, str]`` parameter (the
mzTab adapter has the same TODO; the SDRF reader that backs both is
deferred until the first multi-fraction dataset).
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

from riana import accmass
from riana.exceptions import DataError
from riana.records import PSMRecord


# Crux Percolator emits a fixed-width TSV with these columns; the absence of
# ``spectrum precursor m/z`` is the signal we are looking at standalone output.
_CRUX_SIGNATURE = ("spectrum precursor m/z", "sequence")

# Standalone Percolator PSMId is MSFragger-shaped: ``file.scan.scan.charge_index``.
_PSMID_RE = re.compile(r"^(?P<file>.+)\.(?P<scan>\d+)\.\d+\.(?P<charge>\d+)_\d+$")


def read_percolator(
    path: str | os.PathLike[str],
    sample: str,
    ignored_mods: Sequence[str] = (),
) -> list[PSMRecord]:
    """Parse a Percolator target-PSMs file into typed records.

    Args:
        path: ``percolator.target.psms.txt`` (Crux or standalone format).
        sample: sample label written into every emitted record.
        ignored_mods: mods to ignore when recomputing peptide mass. Forwarded
            to :func:`riana.accmass.calculate_ion_mz`. The recompute exists
            because the Crux ``peptide mass`` column omits cysteine-IAA mass.

    Returns:
        One :class:`PSMRecord` per PSM row. ``pep_id`` is left unassigned
        (``-1``); use :func:`assign_pep_ids_per_fraction` to set it.
    """
    path = Path(path)
    try:
        with open(path, "r") as f:
            header_line = f.readline()
    except OSError as e:
        raise DataError(f"Failed to load percolator file {path}: {e}") from e

    header_cols = header_line.rstrip("\n").split("\t")
    is_crux = all(col in header_cols for col in _CRUX_SIGNATURE)

    if is_crux:
        return _read_crux(path, sample, ignored_mods)
    return _read_standalone(path, sample, ignored_mods)


def _read_crux(
    path: Path, sample: str, ignored_mods: Sequence[str]
) -> list[PSMRecord]:
    try:
        df = pd.read_csv(path, sep="\t")
    except OSError as e:
        raise DataError(f"Failed to load percolator file {path}: {e}") from e

    # Recompute peptide mass so cysteine-IAA mass is always counted (Crux's
    # column omits it). Matches the 0.9.0 ReadPercolator behavior bit-for-bit.
    peptide_masses = [
        accmass.calculate_ion_mz(seq, ignored_mods=ignored_mods)
        for seq in df["sequence"]
    ]

    rows = df.to_dict(orient="records")
    records: list[PSMRecord] = []
    for row, peptide_mass in zip(rows, peptide_masses):
        records.append(
            PSMRecord(
                scan=int(row["scan"]),
                charge=int(row["charge"]),
                sequence=str(row["sequence"]),
                peptide_mass=float(peptide_mass),
                sample=sample,
                file_idx=int(row["file_idx"]),
                protein_id=str(row.get("protein id", "") or ""),
                flanking_aa=str(row.get("flanking aa", "") or ""),
                precursor_mz=float(row.get("spectrum precursor m/z", 0.0)),
                neutral_mass=float(row.get("spectrum neutral mass", 0.0)),
                percolator_score=float(row.get("percolator score", 0.0)),
                percolator_q_value=float(row.get("percolator q-value", 1.0)),
                percolator_pep=float(row.get("percolator PEP", 1.0)),
                distinct_matches=int(row.get("distinct matches/spectrum", 0)),
            )
        )
    return records


def _read_standalone(
    path: Path, sample: str, ignored_mods: Sequence[str]
) -> list[PSMRecord]:
    """Parse standalone Percolator (MSFragger) output."""
    with open(path, "r") as f:
        f_ln = f.readlines()
    if len(f_ln) <= 1:
        return []

    # Standalone-Percolator rows have a variable column count because protein
    # IDs are tab-separated. The first 5 columns are fixed.
    head_df = pd.DataFrame([ln.split("\t")[0:5] for ln in f_ln[1:]])
    head_df.columns = ["PSMId", "score", "qvalue", "pep", "peptide"]
    head_df["qvalue"] = head_df["qvalue"].astype(float)
    head_df["pep"] = head_df["pep"].astype(float)

    # "X.SEQUENCE.Y" -> "SEQUENCE"; strip a trailing bare digit that
    # MSFragger sometimes appends before the closing flanker.
    sequences: list[str] = []
    for pep in head_df["peptide"]:
        seq = pep[2:-2]
        if len(seq) >= 2 and seq[-1].isdigit() and not seq[-2].isdigit():
            seq = seq[:-1]
        sequences.append(seq)
    head_df["sequence"] = sequences
    head_df["flanking_aa"] = [pep[0] + pep[-1] for pep in head_df["peptide"]]

    # PSMId carries (file, scan, charge); parse via the documented regex
    # rather than positional split gymnastics.
    parsed = head_df["PSMId"].astype(str).str.extract(_PSMID_RE)
    if parsed.isnull().any().any():
        bad = head_df.loc[parsed.isnull().any(axis=1), "PSMId"].head(3).tolist()
        raise DataError(
            "Standalone Percolator PSMId did not match expected MSFragger "
            f"format 'filename.scan.scan.charge_index'. Examples: {bad}"
        )

    head_df["scan"] = parsed["scan"].astype(int)
    head_df["charge"] = parsed["charge"].astype(int)
    head_df["file_name"] = parsed["file"].apply(
        lambda p: Path(os.path.basename(p)).stem
    )
    sorted_files = sorted(head_df["file_name"].unique())
    file_idx_map = {name: i for i, name in enumerate(sorted_files)}
    head_df["file_idx"] = head_df["file_name"].map(file_idx_map)

    protein_ids = [",".join(ln.rstrip().split("\t")[5:]) for ln in f_ln[1:]]

    records: list[PSMRecord] = []
    for row, prot in zip(head_df.itertuples(index=False), protein_ids):
        peptide_mass = accmass.calculate_ion_mz(row.sequence, ignored_mods=ignored_mods)
        records.append(
            PSMRecord(
                scan=int(row.scan),
                charge=int(row.charge),
                sequence=row.sequence,
                peptide_mass=float(peptide_mass),
                sample=sample,
                file_idx=int(row.file_idx),
                file_name=row.file_name,
                protein_id=prot,
                flanking_aa=row.flanking_aa,
                # Standalone output does not carry these spectrum fields.
                precursor_mz=0.0,
                neutral_mass=0.0,
                percolator_score=float(row.score),
                percolator_q_value=float(row.qvalue),
                percolator_pep=float(row.pep),
                distinct_matches=0,
            )
        )
    return records


def file_indices(records: Iterable[PSMRecord]) -> list[int]:
    """Sorted list of distinct ``file_idx`` values present in *records*."""
    return sorted({r.file_idx for r in records})


def filter_by_q_value(
    records: Iterable[PSMRecord], threshold: float
) -> list[PSMRecord]:
    """Keep records with ``percolator_q_value < threshold``."""
    return [r for r in records if r.percolator_q_value < threshold]


def fraction_psms(
    records: Iterable[PSMRecord],
    file_idx: int,
    *,
    assign_pep_ids: bool = True,
) -> list[PSMRecord]:
    """Return PSMs for one fraction (``file_idx``), sorted by scan.

    When ``assign_pep_ids`` is true (the 0.9.0 default), the per-fraction
    sequential ``pep_id`` 0..N-1 is written back into the returned records.
    """
    from dataclasses import replace

    rows = sorted(
        (r for r in records if r.file_idx == file_idx), key=lambda r: r.scan
    )
    if not assign_pep_ids:
        return rows
    return [replace(r, pep_id=i) for i, r in enumerate(rows)]
