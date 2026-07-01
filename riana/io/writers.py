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
import functools
import hashlib
import json
import os
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from riana import __version__

#: Float format for the *estimate* outputs (fit peptides / fractions / protein).
#: ~6 significant figures — readable and well beyond the measurement precision of
#: k / θ / R² (reproducibility is guaranteed by the provenance header, not
#: bit-exact floats). NOT used for integrate ``_riana.txt``, whose ``isoN_obs_mz``
#: needs ppm-level digits (and which is parity-gated).
ESTIMATE_FLOAT_FORMAT = "%.6g"


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


def read_provenance_header(path: str | os.PathLike[str]) -> dict[str, str]:
    """Parse the ``# key value`` provenance header of a Riana output file.

    The inverse of :meth:`Provenance.comment_lines`: returns the leading
    ``# key value`` comment lines (``riana``, ``git``, ``config_hash``,
    ``id_source``, and any ``extra`` such as ``model`` / ``label`` / ``method``)
    as a dict, stopping at the first non-comment line. Used when *displaying* a
    prior result (GUI "Load results"): the estimates are in the table, but the
    model that shaped the fitted curve lives only here. Unknown / malformed lines
    are skipped; a value may contain spaces (only the first token is the key).
    """
    header: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            parts = line[1:].strip().split(None, 1)
            if len(parts) == 2:
                header[parts[0]] = parts[1]
    return header


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


@functools.lru_cache(maxsize=1)
def _git_sha() -> str:
    """Return the current commit SHA (12 chars), or ``"unknown"`` if unavailable.

    Cached for the process lifetime: the SHA is constant for a run, and the
    ``subprocess`` call ``fork``s. ``make_provenance`` runs once per output file,
    so without the cache a many-file ``integrate`` would fork ``git`` hundreds of
    times from the (multi-threaded, post-``ProcessPoolExecutor``) main process —
    on macOS that intermittently **deadlocks** in the child's ``pthread_atfork``
    handlers and freezes the run. Warm it once up front (single-threaded, before
    the pool) via :func:`warm_git_sha` so even the single remaining fork is safe.
    """
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


def warm_git_sha() -> str:
    """Prime the :func:`_git_sha` cache while the process is still single-threaded.

    Call this before spawning the integrate/fit worker pool so the one ``git``
    ``fork`` happens in a clean process, not concurrently with pool threads.
    """
    return _git_sha()


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


def write_dataframe_tsv(
    path: str | os.PathLike[str],
    df,
    provenance: Provenance,
    *,
    include_index: bool = True,
    float_format: str | None = None,
) -> None:
    """Write a ``pandas.DataFrame`` to *path* as TSV with the provenance header.

    The header is the same ``# riana | git | config_hash | id_source``
    comment block ``write_tsv`` emits; the body is whatever
    ``DataFrame.to_csv(sep='\\t', index=include_index)`` would produce.

    Bench readers should ``pd.read_csv(path, sep='\\t', comment='#')`` to
    skip the provenance lines — that's the convention M3 Week 4 introduces
    and the helper docstring on ``Provenance`` is the canonical reference.

    Args:
        path: destination path.
        df: pandas DataFrame.
        provenance: from :func:`make_provenance`.
        include_index: write ``df.index`` as the first column. Defaults
            true because the legacy ``_riana.txt`` and
            ``riana_fit_peptides.txt`` schemas both rely on the index.
        float_format: optional ``DataFrame.to_csv`` float format (e.g.
            :data:`ESTIMATE_FLOAT_FORMAT`). ``None`` keeps full float64 precision.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        for line in provenance.comment_lines():
            f.write(line + "\n")
        df.to_csv(f, sep="\t", index=include_index, float_format=float_format)


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
