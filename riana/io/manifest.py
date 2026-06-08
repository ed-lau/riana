"""The stage-aware project manifest — ``riana_manifest.tsv`` (M6a).

The manifest is the **project index** that glues Riana's linear ``integrate →
fit → protein`` chain together (PROJECT_REVIEW.md §3, locked decision #2). Each
stage appends rows tagged with its ``stage`` and the run/group
:class:`~riana.records.RunIdentity`, so a later stage groups its inputs from the
manifest rather than re-reading the SDRF or parsing identity out of filenames.

Identity is *header-authoritative*: integrate freezes the full identity into each
output's provenance header (reproducible, SDRF-independent at fit time) **and**
records it here for indexing. Fit reads ``stage == "integrate"`` rows, groups by
:attr:`RunIdentity.curve_key`, and writes back its own ``stage == "fit"`` rows.

The file is schema-versioned (``# manifest_schema 1`` on the first line) so the
breaking output-granularity change and any future format change are detectable
downstream (PROJECT_REVIEW.md §3, additional recommendation #3).
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

from riana.exceptions import DataError
from riana.records import RunIdentity

MANIFEST_FILENAME = "riana_manifest.tsv"
SCHEMA_VERSION = 1
_SCHEMA_PREFIX = "# manifest_schema"
_STAGES = ("integrate", "fit", "protein")

# Column order on disk. The identity block mirrors RunIdentity's fields; the
# float fields serialize None as "" and round-trip back to None.
_FIELDS = (
    "stage",
    "output_path",
    "experiment",
    "sample",
    "data_file",
    "biological_replicate",
    "technical_replicate",
    "fraction",
    "labeling_time",
    "labeling_time_unit",
    "mixing_proportion",
    "condition",
    "acquisition",
    "precursor_enrichment",
    "experiment_type",
    "config_hash",
    "git_sha",
)


@dataclass(frozen=True, slots=True)
class ManifestRow:
    """One stage output: an :class:`RunIdentity` + where it was written."""

    stage: str
    output_path: str
    identity: RunIdentity
    config_hash: str = ""
    git_sha: str = ""

    def __post_init__(self) -> None:
        if self.stage not in _STAGES:
            raise ValueError(f"stage must be one of {_STAGES}, got {self.stage!r}")

    @property
    def key(self) -> tuple[str, str]:
        """Upsert key — one row per ``(stage, output_path)``."""
        return (self.stage, self.output_path)

    def to_record(self) -> dict[str, str]:
        i = self.identity
        return {
            "stage": self.stage,
            "output_path": self.output_path,
            "experiment": i.experiment,
            "sample": i.sample,
            "data_file": i.data_file,
            "biological_replicate": str(i.biological_replicate),
            "technical_replicate": str(i.technical_replicate),
            "fraction": str(i.fraction),
            "labeling_time": _fmt_opt(i.labeling_time),
            "labeling_time_unit": i.labeling_time_unit,
            "mixing_proportion": _fmt_opt(i.mixing_proportion),
            "condition": i.condition,
            "acquisition": i.acquisition,
            "precursor_enrichment": _fmt_opt(i.precursor_enrichment),
            "experiment_type": i.experiment_type,
            "config_hash": self.config_hash,
            "git_sha": self.git_sha,
        }

    @classmethod
    def from_record(cls, rec: dict[str, str]) -> "ManifestRow":
        identity = RunIdentity(
            experiment=rec.get("experiment", ""),
            sample=rec.get("sample", ""),
            data_file=rec.get("data_file", ""),
            biological_replicate=_to_int(rec.get("biological_replicate"), 1),
            technical_replicate=_to_int(rec.get("technical_replicate"), 1),
            fraction=_to_int(rec.get("fraction"), 1),
            labeling_time=_to_opt_float(rec.get("labeling_time")),
            labeling_time_unit=rec.get("labeling_time_unit", ""),
            mixing_proportion=_to_opt_float(rec.get("mixing_proportion")),
            condition=rec.get("condition", ""),
            acquisition=rec.get("acquisition", "DDA") or "DDA",
            precursor_enrichment=_to_opt_float(rec.get("precursor_enrichment")),
        )
        return cls(
            stage=rec["stage"],
            output_path=rec["output_path"],
            identity=identity,
            config_hash=rec.get("config_hash", ""),
            git_sha=rec.get("git_sha", ""),
        )


def append_manifest(
    path: str | os.PathLike[str], rows: list[ManifestRow]
) -> list[ManifestRow]:
    """Upsert *rows* into the manifest at *path*, returning the full contents.

    Existing rows with the same ``(stage, output_path)`` key are replaced (so a
    re-run is idempotent); everything else is preserved. Creates the file with
    the schema header when absent.
    """
    path = Path(path)
    existing = read_manifest(path) if path.exists() else []
    merged: dict[tuple[str, str], ManifestRow] = {r.key: r for r in existing}
    for r in rows:
        merged[r.key] = r
    ordered = sorted(merged.values(), key=lambda r: (r.stage, r.output_path))
    _write_manifest(path, ordered)
    return ordered


def read_manifest(
    path: str | os.PathLike[str], stage: str | None = None
) -> list[ManifestRow]:
    """Read the manifest, optionally filtered to one ``stage``."""
    path = Path(path)
    if not path.exists():
        raise DataError(f"manifest not found: {path}")
    with open(path, "r", newline="") as f:
        first = f.readline()
        _check_schema(first, path)
        reader = csv.DictReader(f, delimiter="\t")
        rows = [ManifestRow.from_record(rec) for rec in reader]
    if stage is not None:
        if stage not in _STAGES:
            raise ValueError(f"stage must be one of {_STAGES}, got {stage!r}")
        rows = [r for r in rows if r.stage == stage]
    return rows


# --------------------------------------------------------------------------- #
# internals
# --------------------------------------------------------------------------- #
def _write_manifest(path: Path, rows: list[ManifestRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        f.write(f"{_SCHEMA_PREFIX} {SCHEMA_VERSION}\n")
        writer = csv.DictWriter(f, fieldnames=list(_FIELDS), delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r.to_record())


def _check_schema(first_line: str, path: Path) -> None:
    line = first_line.strip()
    if not line.startswith(_SCHEMA_PREFIX):
        raise DataError(
            f"{path} is not a Riana manifest (missing '{_SCHEMA_PREFIX} N' header)."
        )
    try:
        version = int(line[len(_SCHEMA_PREFIX):].strip())
    except ValueError:
        raise DataError(f"{path}: unparseable manifest schema header {line!r}.") from None
    if version > SCHEMA_VERSION:
        raise DataError(
            f"{path} is manifest schema v{version}, but this Riana understands "
            f"v{SCHEMA_VERSION}. Upgrade Riana."
        )


def _fmt_opt(value: float | None) -> str:
    return "" if value is None else repr(value)


def _to_opt_float(value: str | None) -> float | None:
    if value is None or value.strip() == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _to_int(value: str | None, default: int) -> int:
    if value is None or value.strip() == "":
        return default
    try:
        return int(float(value))
    except ValueError:
        return default
