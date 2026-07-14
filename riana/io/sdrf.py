"""SDRF intake — the carrier for Riana's run-identity data model (M6a).

A `Sample and Data Relationship Format <https://github.com/bigbio/proteomics-sample-metadata>`_
TSV is the single source of run identity for both acquisition modes: DDA via
quantms (mzTab) and DIA via quantms-diann (DIA-NN parquet, M6b). Only the
PSM/quant file format differs; the identity model is shared, so this one reader
backs both (PROJECT_REVIEW.md §3, Track A).

Riana reads a **documented subset** of SDRF columns and ignores everything else
quantms fills in:

==============================================  ===============================
SDRF column                                     :class:`~riana.records.RunIdentity`
==============================================  ===============================
``comment[data file]`` (ext stripped)           ``data_file`` — the mzML join key
``source name``                                 ``sample``
``characteristics[biological replicate]``       ``biological_replicate``
``comment[technical replicate]``                ``technical_replicate``
``comment[fraction identifier]``                ``fraction``
``characteristics[labeling time]``  *(turnover)*  ``labeling_time`` (+unit)
``characteristics[mixing proportion]`` *(calib)*  ``mixing_proportion``
``factor value[...]``                            ``condition``
``comment[proteomics data acquisition method]``  ``acquisition`` (DDA/DIA)
``characteristics[precursor enrichment]``        ``precursor_enrichment`` (RIA)
``comment[precursor mass tolerance]``            :attr:`SdrfTable.precursor_mass_tol_ppm`
``comment[modification parameters]`` (1..n)      :attr:`SdrfTable.modifications`
==============================================  ===============================

Locked decisions (PROJECT_REVIEW.md §3):

- The **independent-variable column declares the experiment type**: exactly one
  of ``characteristics[labeling time]`` (turnover — the kinetic-curve x-axis) or
  ``characteristics[mixing proportion]`` (calibration fixtures) must be present.
  Fit dispatches on :attr:`SdrfTable.experiment_type`.
- Labeling time is a *characteristic* (sample-intrinsic), deliberately not
  ``factor value[time]`` — that avoids colliding with a drug-treatment time
  course where the *treatment* time is the genuine factor value.
- ``comment[precursor mass tolerance]`` **is** read into
  :attr:`SdrfTable.precursor_mass_tol_ppm` and used as the integration mass
  tolerance (CLI ``--mass_tol`` overrides; dataclass default is the fallback).
  This **reverses** the earlier M3 rec4 decision (which treated the search window
  as too tight and used a deliberately wider ~50 ppm window). That was a
  profile-vs-centroid category error: on **centroid** mzML (the quantms/
  ThermoRawFileParser norm) each isotopomer is one line spread only by mass
  accuracy (~3–10 ppm), so the search tolerance *is* the right integration
  window; a ~50 ppm "peak-width" window just imports co-eluting interference
  (verified on LVE: out-of-range θ 27%→8%, R²med 0.69→0.86 tightening 50→10).
  Only a ``ppm`` unit is honored; a ``Da`` value is ignored (can't map to a ppm
  window) with a warning.

Variable modifications are parsed and exposed on :attr:`SdrfTable.modifications`
but not yet threaded into the isotope-envelope model — that is the PTM-aware
M7 work; today the FS solve uses the unmodified backbone and the fixed
Carbamidomethyl(C) handled by :func:`riana.algorithms.mass_calc.calculate_ion_mz`.
"""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, replace
from pathlib import Path

from riana import multiplex
from riana.exceptions import DataError
from riana.records import RunIdentity

_LOGGER = logging.getLogger(__name__)

# mzML/raw/vendor extensions stripped from ``comment[data file]`` to recover the
# stem that joins to the mzTab ``ms_run[N]-location`` basename and the on-disk
# mzML. ``.mzML.gz`` is handled before the single-extension fallback.
_DOUBLE_EXTS = (".mzml.gz",)
_SINGLE_EXTS = (".raw", ".mzml", ".d", ".wiff", ".wiff2", ".dia", ".gz")

_LABELING_TIME_RE = re.compile(r"^\s*([+-]?[0-9]*\.?[0-9]+)\s*(.*?)\s*$")
_NOT_APPLICABLE = {"", "not applicable", "not available", "na", "n/a", "null"}

_LABELING_TIME_COL = "characteristics[labeling time]"
_MIXING_PROPORTION_COL = "characteristics[mixing proportion]"
_DATA_FILE_COL = "comment[data file]"
_SOURCE_NAME_COL = "source name"
_LABEL_COL = "comment[label]"

# An isobaric (TMT / iTRAQ) channel label in ``comment[label]``. Its presence means
# one mzML multiplexes several samples into ONE MS1 cluster, so Riana collapses the
# per-channel rows to a single per-file run (``_collapse_isobaric_runs``) — it
# measures the samples' intensity-weighted-average turnover, not per channel.
_ISOBARIC_LABEL_RE = re.compile(r"^\s*(tmt|itraq)", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class Modification:
    """A parsed ``comment[modification parameters]`` value (UNIMOD spec)."""

    name: str  #: ``NT=`` modification name, e.g. ``Oxidation``.
    residue: str  #: ``TA=`` target amino acid(s), e.g. ``M`` (empty = any).
    accession: str  #: ``AC=`` accession, e.g. ``UNIMOD:35``.
    mod_type: str  #: ``MT=`` ``Fixed`` or ``Variable``.


@dataclass(frozen=True, slots=True)
class SdrfTable:
    """Parsed SDRF — one :class:`RunIdentity` per row plus shared metadata."""

    experiment: str
    runs: tuple[RunIdentity, ...]
    modifications: tuple[Modification, ...]
    #: ``"turnover"`` or ``"calibration"`` — consistent across all runs.
    experiment_type: str
    #: ``"DDA"`` or ``"DIA"`` — consistent across all runs.
    acquisition: str
    #: ``comment[precursor mass tolerance]`` in ppm (search tolerance), used as
    #: the integration mass window unless CLI ``--mass_tol`` overrides. ``None``
    #: when the column is absent or given in Da (not a ppm window).
    precursor_mass_tol_ppm: float | None = None
    #: For a **sample-axis-multiplexed** SDRF (dimethyl/SILAC — one mzML carries N
    #: channels at distinct precursor masses, kept as distinct samples), the
    #: ``{(data_file_stem, (label, channel)): RunIdentity}`` map the mzTab reader
    #: routes each PSM through, by its label mod. ``None`` for a normal one-run-per-
    #: file SDRF. See :mod:`riana.multiplex`.
    multiplex_channel_map: dict[tuple[str, tuple[str, str]], RunIdentity] | None = None

    @property
    def is_multiplexed(self) -> bool:
        """True when this SDRF declares sample-axis multiplexing (dimethyl/SILAC) —
        i.e. one mzML holds several channels the intake keeps as distinct samples."""
        return self.multiplex_channel_map is not None

    @property
    def sample_map(self) -> dict[str, RunIdentity]:
        """``{data_file_stem: RunIdentity}`` — the join key into mzTab/DIA-NN.

        Keyed by the (extension-stripped) data-file stem rather than a positional
        ms_run index so the join is order-independent (the SDRF row order need
        not match the mzTab ``ms_run[N]`` order).
        """
        return {r.data_file: r for r in self.runs}

    @property
    def variable_modifications(self) -> tuple[Modification, ...]:
        return tuple(m for m in self.modifications if m.mod_type.lower() == "variable")

    @property
    def fixed_modifications(self) -> tuple[Modification, ...]:
        return tuple(m for m in self.modifications if m.mod_type.lower() == "fixed")


def read_sdrf(path: str | os.PathLike[str], experiment: str | None = None) -> SdrfTable:
    """Parse an SDRF TSV into an :class:`SdrfTable`.

    Args:
        path: the ``*.sdrf.tsv`` file.
        experiment: the top-level experiment label stamped on every
            :class:`RunIdentity`. Defaults to the SDRF file stem (with a
            trailing ``.sdrf`` stripped).

    Raises:
        DataError: missing required columns, neither/both experiment-type
            columns present, a non-numeric independent value, or a duplicate
            data file (an ambiguous mzML join).
    """
    path = Path(path)
    if not path.exists():
        raise DataError(f"SDRF not found: {path}")
    if experiment is None:
        experiment = re.sub(r"\.sdrf$", "", path.stem, flags=re.IGNORECASE)

    header, rows = _read_tsv(path)
    col_index = _index_columns(header)  # name -> [col positions]

    _require_column(col_index, _SOURCE_NAME_COL, path)
    _require_column(col_index, _DATA_FILE_COL, path)

    has_time = _LABELING_TIME_COL in col_index
    has_mix = _MIXING_PROPORTION_COL in col_index
    if has_time == has_mix:
        which = "both" if has_time else "neither"
        raise DataError(
            f"SDRF {path} must carry exactly one of "
            f"'{_LABELING_TIME_COL}' (turnover) or '{_MIXING_PROPORTION_COL}' "
            f"(calibration); found {which}."
        )
    experiment_type = "turnover" if has_time else "calibration"

    factor_cols = [name for name in col_index if name.startswith("factor value[")]
    mod_cols = col_index.get("comment[modification parameters]", [])

    # Isobaric (TMT/iTRAQ) SDRF? Then one data file legitimately recurs once per
    # channel; the rows are collapsed to a single per-file run after the loop.
    label_cols = col_index.get(_LABEL_COL, [])

    def _label_cell(row: list[str]) -> str:
        return row[label_cols[0]] if label_cols and label_cols[0] < len(row) else ""

    is_isobaric = bool(label_cols) and any(
        _ISOBARIC_LABEL_RE.match(_label_cell(row)) for row in rows
    )
    # Sample-axis multiplexing (dimethyl/SILAC): one mzML carries several channels at
    # DISTINCT precursor masses, so — unlike isobaric TMT — the channel rows are KEPT
    # as separate samples (each fit at its own enrichment). Detected via the registry
    # ``comment[label]`` CV terms. Mutually exclusive with isobaric.
    is_multiplex = not is_isobaric and bool(label_cols) and any(
        multiplex.is_multiplex_cv(_label_cell(row)) for row in rows
    )
    channel_map: dict[tuple[str, tuple[str, str]], RunIdentity] = {}

    runs: list[RunIdentity] = []
    acquisitions: set[str] = set()
    mass_tol_ppms: set[float] = set()
    seen_data_files: dict[str, int] = {}
    modifications: list[Modification] = []
    seen_mods: set[Modification] = set()

    for line_no, row in enumerate(rows, start=2):  # line 1 is the header

        def cell(name: str) -> str:
            cols = col_index.get(name)
            if not cols:
                return ""
            return row[cols[0]].strip() if cols[0] < len(row) else ""

        data_file = _strip_data_ext(cell(_DATA_FILE_COL))
        if not data_file:
            raise DataError(f"{path} line {line_no}: empty '{_DATA_FILE_COL}'.")
        if data_file in seen_data_files and not (is_isobaric or is_multiplex):
            raise DataError(
                f"{path}: data file {data_file!r} appears on lines "
                f"{seen_data_files[data_file]} and {line_no} — the mzML join "
                "key must be unique per run (unless isobaric/TMT-multiplexed or "
                "sample-axis-multiplexed, e.g. dimethyl/SILAC)."
            )
        seen_data_files[data_file] = line_no

        labeling_time, time_unit, mixing = None, "", None
        if has_time:
            labeling_time, time_unit = _parse_labeling_time(
                cell(_LABELING_TIME_COL), path, line_no
            )
        else:
            mixing = _parse_float(
                cell(_MIXING_PROPORTION_COL), _MIXING_PROPORTION_COL, path, line_no
            )

        acquisition = _parse_acquisition(
            cell("comment[proteomics data acquisition method]")
        )
        acquisitions.add(acquisition)

        mt_ppm = _parse_mass_tol_ppm(cell("comment[precursor mass tolerance]"))
        if mt_ppm is not None:
            mass_tol_ppms.add(mt_ppm)

        condition = "|".join(
            v for c in factor_cols if (v := cell(c)) and v.lower() not in _NOT_APPLICABLE
        )

        identity = RunIdentity(
            experiment=experiment,
            sample=cell(_SOURCE_NAME_COL),
            data_file=data_file,
            biological_replicate=_parse_int(
                cell("characteristics[biological replicate]"), default=1
            ),
            technical_replicate=_parse_int(
                cell("comment[technical replicate]"), default=1
            ),
            fraction=_parse_int(cell("comment[fraction identifier]"), default=1),
            labeling_time=labeling_time,
            labeling_time_unit=time_unit,
            mixing_proportion=mixing,
            condition=condition,
            acquisition=acquisition,
            precursor_enrichment=_parse_optional_float(
                cell("characteristics[precursor enrichment]")
            ),
        )
        runs.append(identity)

        # Sample-axis multiplexing: map (file stem, channel) → this row's identity so
        # the mzTab reader can route each PSM to its channel by the peptidoform's mod.
        if is_multiplex:
            channel = multiplex.label_channel_for_cv(_label_cell(row))
            if channel is None:
                raise DataError(
                    f"{path} line {line_no}: comment[label] "
                    f"{_label_cell(row)!r} is not a recognised multiplexing channel "
                    "(a sample-axis-multiplexed sheet must use registry CV terms, "
                    "e.g. NT=DIMETHYL0/DIMETHYL8, on every row)."
                )
            key = (data_file, channel)
            if key in channel_map:
                raise DataError(
                    f"{path} line {line_no}: channel {channel} of file "
                    f"{data_file!r} is declared more than once."
                )
            channel_map[key] = identity

        # Modifications are shared across rows in practice; dedup but keep order.
        for col in mod_cols:
            raw = row[col].strip() if col < len(row) else ""
            mod = _parse_modification(raw)
            if mod is not None and mod not in seen_mods:
                seen_mods.add(mod)
                modifications.append(mod)

    if not runs:
        raise DataError(f"SDRF {path} has a header but no data rows.")
    if len(acquisitions) > 1:
        raise DataError(
            f"SDRF {path} mixes acquisition methods {sorted(acquisitions)}; "
            "Riana expects one method per SDRF."
        )

    if is_isobaric:
        n_channels = len(runs)
        runs = _collapse_isobaric_runs(runs)
        _LOGGER.info(
            "SDRF %s: isobaric (TMT/iTRAQ) — collapsed %d channel rows to %d "
            "per-file runs (one MS1 measurement per file = the multiplexed "
            "average).", path, n_channels, len(runs),
        )
    elif is_multiplex:
        # Sample-axis multiplexing: KEEP every channel row (opposite of the isobaric
        # collapse) — each channel is a distinct sample fit at its own enrichment.
        _LOGGER.info(
            "SDRF %s: sample-axis multiplexed (dimethyl/SILAC) — %d channel rows "
            "across %d files kept as distinct samples (routed to channels %s).",
            path, len(runs), len(seen_data_files),
            sorted({ch for _, ch in channel_map}),
        )

    # Precursor mass tolerance (search window) → integration tolerance. Uniform
    # across runs in practice; if rows disagree, use the smallest (tightest) and
    # warn. A present-but-Da column resolves to no ppm value (warn once).
    prec_mass_tol_ppm: float | None = None
    if "comment[precursor mass tolerance]" in col_index:
        if mass_tol_ppms:
            prec_mass_tol_ppm = min(mass_tol_ppms)
            if len(mass_tol_ppms) > 1:
                _LOGGER.warning(
                    "SDRF %s: multiple precursor mass tolerances %s ppm; using "
                    "tightest (%g ppm).", path, sorted(mass_tol_ppms),
                    prec_mass_tol_ppm,
                )
        else:
            _LOGGER.warning(
                "SDRF %s: 'comment[precursor mass tolerance]' present but no "
                "usable ppm value (Da given?); falling back to the integration "
                "default.", path,
            )

    return SdrfTable(
        experiment=experiment,
        runs=tuple(runs),
        modifications=tuple(modifications),
        experiment_type=experiment_type,
        acquisition=acquisitions.pop(),
        precursor_mass_tol_ppm=prec_mass_tol_ppm,
        multiplex_channel_map=(channel_map if is_multiplex else None),
    )


def _collapse_isobaric_runs(runs: list[RunIdentity]) -> list[RunIdentity]:
    """Collapse an isobaric SDRF's per-channel rows into one run per data file.

    The multiplexed samples co-elute in a SINGLE MS1 cluster, so Riana measures
    their intensity-weighted-average D2O turnover — one measurement per file, not
    per channel. The merged run's ``condition`` is the ``|``-joined **set** of the
    channels' conditions:

    - a file whose channels are all ONE treatment keeps that clean label — so a
      design that splits treatments into separate TMT batches (all-control files vs
      all-treatment files) flows straight into the linear-simple 2-sample Δk;
    - a file that POOLS treatments into one plex (this dataset: control + nocodazole
      together) carries the combined label and is not separable at MS1.

    All other identity fields are file-level (they agree across channels — same
    replicate / fraction / labeling time / enrichment); the per-channel source name
    is replaced by the data-file stem. First-seen file order is preserved.
    """
    groups: dict[str, list[RunIdentity]] = {}
    for r in runs:
        groups.setdefault(r.data_file, []).append(r)
    collapsed: list[RunIdentity] = []
    for data_file, group in groups.items():
        condition = "|".join(sorted({r.condition for r in group if r.condition}))
        collapsed.append(replace(group[0], condition=condition, sample=data_file))
    return collapsed


# --------------------------------------------------------------------------- #
# parsing helpers
# --------------------------------------------------------------------------- #
def _read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    try:
        with open(path, "r", newline="") as f:
            lines = f.read().splitlines()
    except OSError as e:
        raise DataError(f"failed to read SDRF {path}: {e}") from e
    if not lines:
        raise DataError(f"SDRF {path} is empty.")
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:] if ln.strip()]
    return header, rows


def _index_columns(header: list[str]) -> dict[str, list[int]]:
    """Map each (lower-cased) column name to its position(s).

    SDRF legitimately repeats column names — notably
    ``comment[modification parameters]`` (one per fixed/variable mod) — so a
    name maps to a *list* of positions, not one.
    """
    out: dict[str, list[int]] = {}
    for i, name in enumerate(header):
        out.setdefault(name.strip().lower(), []).append(i)
    return out


def _require_column(col_index: dict[str, list[int]], name: str, path: Path) -> None:
    if name not in col_index:
        raise DataError(f"SDRF {path} is missing the required column '{name}'.")


def _strip_data_ext(name: str) -> str:
    base = os.path.basename(name.strip())
    low = base.lower()
    for ext in _DOUBLE_EXTS:
        if low.endswith(ext):
            return base[: -len(ext)]
    for ext in _SINGLE_EXTS:
        if low.endswith(ext):
            return base[: -len(ext)]
    return base


def _parse_labeling_time(value: str, path: Path, line_no: int) -> tuple[float, str]:
    m = _LABELING_TIME_RE.match(value)
    if m is None:
        raise DataError(
            f"{path} line {line_no}: could not parse '{_LABELING_TIME_COL}' "
            f"value {value!r} (expected e.g. '14 days')."
        )
    return float(m.group(1)), m.group(2).strip()


def _parse_float(value: str, col: str, path: Path, line_no: int) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        raise DataError(
            f"{path} line {line_no}: '{col}' value {value!r} is not numeric."
        ) from None


def _parse_optional_float(value: str) -> float | None:
    if value.strip().lower() in _NOT_APPLICABLE:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_mass_tol_ppm(value: str) -> float | None:
    """Parse a precursor-mass-tolerance cell to ppm.

    Accepts ``"10 ppm"``, ``"10ppm"``, or a bare ``"10"`` (assumed ppm). Returns
    ``None`` for empty/NA *or* a non-ppm unit (e.g. ``"0.02 Da"``) — a Da window
    can't be applied as a ppm tolerance, so the caller falls back to the default.
    """
    v = value.strip().lower()
    if v in _NOT_APPLICABLE:
        return None
    m = re.match(r"^([0-9]*\.?[0-9]+)\s*([a-z/]+)?$", v)
    if m is None:
        return None
    unit = (m.group(2) or "ppm").strip()
    if unit != "ppm":
        return None
    return float(m.group(1))


def _parse_int(value: str, default: int) -> int:
    v = value.strip().lower()
    if v in _NOT_APPLICABLE:
        return default
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _parse_acquisition(value: str) -> str:
    return "DIA" if "independent" in value.lower() else "DDA"


def _parse_modification(raw: str) -> Modification | None:
    """Parse ``NT=Oxidation;TA=M;AC=UNIMOD:35;MT=Variable`` into a tuple."""
    if not raw or raw.strip().lower() in _NOT_APPLICABLE:
        return None
    fields: dict[str, str] = {}
    for part in raw.split(";"):
        if "=" in part:
            k, _, v = part.partition("=")
            fields[k.strip().upper()] = v.strip()
    if "NT" not in fields:
        return None
    return Modification(
        name=fields.get("NT", ""),
        residue=fields.get("TA", ""),
        accession=fields.get("AC", ""),
        mod_type=fields.get("MT", ""),
    )
