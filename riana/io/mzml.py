"""Indexed / streaming mzML reader, backed by :mod:`pyteomics.mzml`.

Replaces the 0.9.0 :mod:`riana.spectra` reader, which used ``pymzml`` and
loaded every MS1 spectrum's peak arrays eagerly into ``msdata`` (a Python
list of dicts of ``np.ndarray``) — the source of Riana's memory blow-up on
2 GB mzMLs. Here we keep only the per-MS1 ``(scan, rt)`` index in memory and
fetch peaks lazily per scan, or stream them in order via :meth:`ms1_iter`.

The ``scan_idx`` / ``rt_idx`` attributes are kept compatible with the
:mod:`riana.riana_integrate` consumer (``np.searchsorted`` / boolean-mask
indexing patterns at riana_integrate.py:336–352), so a Week-3 swap-over is
a localized change.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import tempfile
from pathlib import Path
from typing import Iterator

import numpy as np
from pyteomics import mzml

from riana.exceptions import DataError


# Thermo nativeID: ``controllerType=0 controllerNumber=1 scan=1432``.
# WIFF/Sciex format uses ``cycle=N experiment=M``; we don't handle that here
# (none of Riana's targets ship Sciex mzML today).
_SCAN_RE = re.compile(r"scan=(\d+)")


class IndexedMzML:
    """Indexed mzML reader exposing the MS1 (scan, rt) table + lazy peak access.

    On construction, the file is scanned once to build the MS1 index — peak
    arrays are decoded but discarded so memory peaks at one spectrum, not
    the whole run.

    Attributes:
        path: source mzML path (``.mzML`` or ``.mzML.gz``).
        scan_idx: ``np.ndarray[int64]`` of MS1 scan numbers, in file order.
        rt_idx: ``np.ndarray[float64]`` of MS1 retention times, **in minutes**.
    """

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = Path(path)
        if not self.path.exists():
            raise DataError(f"mzML not found: {self.path}")

        # pyteomics' offset index needs seekable bytes; a gzip stream is not
        # randomly seekable. Decompress to a temp file when the input is .gz
        # so peak fetch stays O(1). Memory still peaks at one MS1 spectrum.
        self._tmpdir: tempfile.TemporaryDirectory | None = None
        self._reader_path = self._materialize_seekable(self.path)

        scans: list[int] = []
        rts_min: list[float] = []
        spec_ids: list[str] = []
        with mzml.MzML(str(self._reader_path), use_index=True) as reader:
            for spec in reader:
                if spec.get("ms level") != 1:
                    continue
                scans.append(_scan_from_id(spec["id"]))
                rts_min.append(_rt_minutes(spec))
                spec_ids.append(spec["id"])

        self.scan_idx = np.asarray(scans, dtype=np.int64)
        self.rt_idx = np.asarray(rts_min, dtype=np.float64)
        self._scan_to_spec_id: dict[int, str] = dict(zip(scans, spec_ids))

        # Long-lived random-access reader. pyteomics keeps an offset index;
        # ``get_by_id`` decodes one spectrum on demand.
        self._reader = mzml.MzML(str(self._reader_path), use_index=True)

    def _materialize_seekable(self, path: Path) -> Path:
        """Return a path to a seekable mzML; decompress .gz to a temp file."""
        if path.suffix.lower() != ".gz":
            return path
        self._tmpdir = tempfile.TemporaryDirectory(prefix="riana-mzml-")
        out = Path(self._tmpdir.name) / path.with_suffix("").name
        with gzip.open(path, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return out

    def peaks(self, scan: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (m/z, intensity) arrays for the MS1 scan number *scan*."""
        try:
            spec_id = self._scan_to_spec_id[scan]
        except KeyError as e:
            raise DataError(
                f"scan {scan} is not an MS1 in {self.path.name}"
            ) from e
        spec = self._reader.get_by_id(spec_id)
        return spec["m/z array"], spec["intensity array"]

    def ms1_iter(self) -> Iterator[tuple[int, float, np.ndarray, np.ndarray]]:
        """Yield ``(scan, rt_minutes, mz_array, intensity_array)`` per MS1.

        A fresh streaming pass — one spectrum decoded at a time. Use this
        instead of ``peaks()`` when the caller will touch most of the run.
        """
        with mzml.MzML(str(self._reader_path), use_index=True) as reader:
            for spec in reader:
                if spec.get("ms level") != 1:
                    continue
                yield (
                    _scan_from_id(spec["id"]),
                    _rt_minutes(spec),
                    spec["m/z array"],
                    spec["intensity array"],
                )

    def close(self) -> None:
        self._reader.close()
        if self._tmpdir is not None:
            self._tmpdir.cleanup()
            self._tmpdir = None

    def __enter__(self) -> "IndexedMzML":
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()

    def __repr__(self) -> str:
        return f"IndexedMzML(path={self.path}, ms1={len(self.scan_idx)})"


def _scan_from_id(spec_id: str) -> int:
    """Pull the integer scan number out of a Thermo nativeID-style id."""
    m = _SCAN_RE.search(spec_id)
    if m is None:
        raise DataError(
            f"could not extract scan number from spectrum id '{spec_id}' "
            "(non-Thermo nativeID formats are not yet supported)"
        )
    return int(m.group(1))


def _rt_minutes(spec: dict) -> float:
    """RT in minutes, converting from seconds when the mzML declares that unit."""
    scan = spec["scanList"]["scan"][0]
    rt = scan["scan start time"]
    # pyteomics exposes the value as a `unitfloat`; the unit cvParam comes
    # through as `.unit_info`. Files in the wild use either second or minute.
    unit = getattr(rt, "unit_info", "") or ""
    value = float(rt)
    if "second" in unit.lower():
        value /= 60.0
    return value
