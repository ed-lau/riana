"""TSV / JSON writers with provenance headers.

Every Riana output gets a ``# riana ... | git ... | config_hash ...``
comment header on its first line(s), per the M3 cross-cutting recommendation
(PROJECT_REVIEW.md §4.1). The header makes a result file self-identifying:
re-running with the same git SHA and config_hash should reproduce the
numbers bit-for-bit (the integration core is deterministic).

Used by the Week 4 rewrite of ``riana integrate`` / ``riana fit``; the
Week 2 :mod:`bench_id_path` benchmark also writes through here so its
output is provenance-stamped from day one.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from riana import __version__


@dataclass(frozen=True, slots=True)
class Provenance:
    """Identifying metadata stamped into the header of every output file.

    ``config_hash`` is the SHA-256 of a canonicalised JSON dump of the input
    config — same inputs ⇒ same hash. ``id_source`` is the path to the PSM
    file the result depends on (Percolator or mzTab). ``extra`` is open-ended
    for benchmark-specific bookkeeping.
    """

    riana_version: str
    git_sha: str
    config_hash: str
    id_source: str = ""
    extra: Mapping[str, str] = ()  # type: ignore[assignment]

    def comment_lines(self) -> list[str]:
        """The lines (each starting ``# ``) to prepend to a TSV output."""
        lines = [
            f"# riana {self.riana_version}",
            f"# git {self.git_sha}",
            f"# config_hash {self.config_hash}",
        ]
        if self.id_source:
            lines.append(f"# id_source {self.id_source}")
        for k, v in dict(self.extra).items():
            lines.append(f"# {k} {v}")
        return lines


def make_provenance(
    config: Mapping[str, object],
    *,
    id_source: str | os.PathLike[str] = "",
    extra: Mapping[str, str] | None = None,
) -> Provenance:
    """Build a :class:`Provenance` for the current run."""
    return Provenance(
        riana_version=__version__,
        git_sha=_git_sha(),
        config_hash=hash_config(config),
        id_source=str(id_source),
        extra=dict(extra or {}),
    )


def hash_config(config: Mapping[str, object]) -> str:
    """SHA-256 (16 hex chars) of a canonical JSON dump of *config*."""
    payload = json.dumps(_canonicalise(dict(config)), sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _canonicalise(obj: object) -> object:
    """Recursively turn os.PathLike into str so json.dumps is stable."""
    if isinstance(obj, dict):
        return {k: _canonicalise(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_canonicalise(x) for x in obj]
    if isinstance(obj, os.PathLike):
        return os.fspath(obj)
    return obj


def _git_sha() -> str:
    """Return the current commit SHA (12 chars), or ``"unknown"`` if unavailable."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short=12", "HEAD"],
            cwd=Path(__file__).resolve().parent,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip() or "unknown"
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"


def write_tsv(
    path: str | os.PathLike[str],
    columns: Sequence[str],
    rows: Sequence[Mapping[str, object]],
    provenance: Provenance,
) -> None:
    """Write *rows* to *path* as TSV with the provenance header on top."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        for line in provenance.comment_lines():
            f.write(line + "\n")
        writer = csv.DictWriter(
            f, fieldnames=list(columns), delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(
    path: str | os.PathLike[str],
    payload: Mapping[str, object],
    provenance: Provenance,
) -> None:
    """Write *payload* as JSON; provenance is nested under ``"_provenance"``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    body: dict[str, object] = {
        "_provenance": {
            "riana_version": provenance.riana_version,
            "git_sha": provenance.git_sha,
            "config_hash": provenance.config_hash,
            "id_source": provenance.id_source,
            **dict(provenance.extra),
        },
        **dict(payload),
    }
    with open(path, "w") as f:
        json.dump(body, f, indent=2, default=str)
        f.write("\n")
