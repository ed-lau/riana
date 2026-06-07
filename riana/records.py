# -*- coding: utf-8 -*-

""" Typed data records for the M3 (1.0.0) pipeline.

These replace the untyped ``pandas`` rows / positional lists the 0.9.0 pipeline
passes around. They are *skeleton* definitions: introduced in M3 Week 1 so later
weeks have stable types to build against, and refined when their consumers are
written — ``core/integration.py`` (Week 3) and ``core/fitting.py`` (Week 4).

Field names are deliberately anchored to the existing ``*_riana.txt`` column
schema (and the Percolator ``id_df`` columns) so the M2/Week-0 benchmarks, which
parse those files, keep working unchanged through the rewrite. New columns the
rewrite adds (mass-accuracy outputs) are strict additions. See PROJECT_REVIEW.md
§3 and the M2 findings.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class RunIdentity:
    """The full identity of a single MS run (one mzML / one PSM-file row).

    The M6a "run-identity data model" (PROJECT_REVIEW.md §3, Track A). Today
    identity is positional and conventional — ``file_idx`` is a sort index and
    the timepoint is parsed out of the ``sample`` string — which breaks the
    moment an experiment has biological replicates, multiple conditions, or one
    mzTab spanning many runs. This record is sourced from the SDRF
    (:mod:`riana.io.sdrf`) at intake, attached to every :class:`PSMRecord`, and
    carried *header-authoritatively* through the stage-aware
    ``riana_manifest.tsv`` to fit and protein rollup.

    The experiment type is *declared by which independent-variable column is
    present*: ``labeling_time`` (turnover — the kinetic-curve x-axis) XOR
    ``mixing_proportion`` (calibration fixtures). Fit dispatches on
    :attr:`experiment_type`.
    """

    #: Top-level experiment grouping. Distinguishes studies when several
    #: manifests/SDRFs are combined; constant within one SDRF. Defaults to the
    #: SDRF file stem in :func:`riana.io.sdrf.read_sdrf`.
    experiment: str
    #: SDRF ``source name`` — the biological sample this run came from.
    sample: str
    #: ``comment[data file]`` with its extension stripped — the join key to the
    #: mzML stem (and to the mzTab ``ms_run[N]-location`` basename).
    data_file: str
    #: ``characteristics[biological replicate]`` — genuine replicates (distinct
    #: animals/cultures). Load-bearing for the biorep-aware protein refit.
    biological_replicate: int = 1
    #: ``comment[technical replicate]``.
    technical_replicate: int = 1
    #: ``comment[fraction identifier]`` — fractions of one sample are merged at
    #: peptide level at fit time.
    fraction: int = 1
    #: ``characteristics[labeling time]``, numeric (turnover x-axis). ``None``
    #: for calibration fixtures. Mutually exclusive with ``mixing_proportion``.
    labeling_time: float | None = None
    #: Unit token parsed alongside ``labeling_time`` (e.g. ``"day"``), metadata
    #: only — fitting uses the numeric value.
    labeling_time_unit: str = ""
    #: ``characteristics[mixing proportion]`` (calibration x-axis). ``None`` for
    #: turnover. Mutually exclusive with ``labeling_time``.
    mixing_proportion: float | None = None
    #: ``factor value[...]`` (condition / group) for side-by-side comparison and
    #: future cross-group stats. Empty when no factor-value column is present.
    condition: str = ""
    #: ``comment[proteomics data acquisition method]`` collapsed to ``"DDA"`` or
    #: ``"DIA"``.
    acquisition: str = "DDA"
    #: ``characteristics[precursor enrichment]`` (RIA), optional. ``None`` falls
    #: back to the global ``--ria`` default at fit time.
    precursor_enrichment: float | None = None

    @property
    def experiment_type(self) -> str:
        """``"calibration"`` if a mixing proportion is set, else ``"turnover"``."""
        return "calibration" if self.mixing_proportion is not None else "turnover"

    @property
    def independent_value(self) -> float | None:
        """The fit x-axis value — ``mixing_proportion`` or ``labeling_time``."""
        return (
            self.mixing_proportion
            if self.experiment_type == "calibration"
            else self.labeling_time
        )

    @property
    def curve_key(self) -> tuple[str, str, int]:
        """The grouping key for one kinetic curve: ``(experiment, sample,
        biological_replicate)``. Fractions and timepoints of the same curve
        share this key; different bioreps are genuine replicates (kept apart)."""
        return (self.experiment, self.sample, self.biological_replicate)


@dataclass(frozen=True, slots=True)
class PSMRecord:
    """A single peptide-spectrum match from an identification search.

    Maps to the Percolator ``id_df`` columns produced by
    :class:`riana.peptides.ReadPercolator` and carried into ``*_riana.txt``.
    The ``io/percolator.py`` and ``io/mztab.py`` parsers (Week 2) will emit
    these; both ID paths must populate the same record so the downstream
    pipeline is parser-agnostic.
    """

    scan: int
    charge: int
    sequence: str
    #: riana-recalculated monoisotopic neutral mass (Da) — the ``peptide mass``
    #: column. Recomputed from sequence so cysteine-IAA mass is always included.
    peptide_mass: float
    sample: str
    file_idx: int = 0
    #: basename (no extension) of the source mzML. Populated by the standalone
    #: Percolator path (parsed from PSMId) and the mzTab path (from
    #: ``MTD ms_run[N]-location``); empty for Crux Percolator output, which
    #: only carries ``file_idx``. Required join key for ``bench_id_path.py``.
    file_name: str = ""
    protein_id: str = ""
    flanking_aa: str = ""
    #: ``spectrum precursor m/z`` — 0.0 when the standalone-Percolator path
    #: cannot supply it.
    precursor_mz: float = 0.0
    #: ``spectrum neutral mass`` — 0.0 when unavailable (standalone Percolator).
    neutral_mass: float = 0.0
    percolator_score: float = 0.0
    percolator_q_value: float = 1.0
    percolator_pep: float = 1.0
    #: ``distinct matches/spectrum``.
    distinct_matches: int = 0
    #: per-fraction sequential id assigned by the parser (the ``pep_id`` column).
    pep_id: int = -1
    #: how the PSM entered the result: ``"q_value"`` (direct ID) or ``"mbr"``
    #: (match-between-runs, when MBR returns post-0.9.0).
    evidence: str = "q_value"
    #: PSM retention time in seconds (mzTab convention). 0.0 when unavailable
    #: (standalone Percolator). Added in M6a as the DIA RT-apex prior (locked
    #: decision #3) so the M6b DIA path is designed in, not bolted on; the DDA
    #: integration path uses the MS2 scan, not this field.
    retention_time: float = 0.0
    #: The full SDRF-sourced run identity (:class:`RunIdentity`). ``None`` on the
    #: demoted single-mzML Percolator path, which has no SDRF.
    identity: RunIdentity | None = None

    @property
    def concat(self) -> str:
        """``sequence_charge`` identifier — the ``concat`` column the benchmarks
        group on. Derived here so it cannot drift from sequence/charge."""
        return f"{self.sequence}_{self.charge}"


@dataclass(frozen=True, slots=True)
class Chromatogram:
    """An extracted-ion chromatogram (XIC) for one isotopomer of one peptide.

    The per-scan intensity trace across the MS1 scans inside the RT window.
    This is the intermediate the 0.9.0 pipeline integrates with ``np.trapezoid``;
    Week 3 peak detection consumes it to find boundaries before integrating.
    """

    #: 0-based isotopomer index (0 = monoisotopic).
    isotopomer: int
    #: theoretical m/z this trace was extracted at.
    target_mz: float
    #: ±ppm half-width used for extraction (0.9.0 ``-m`` semantic).
    mass_tol_ppm: float
    scans: tuple[int, ...] = ()
    #: retention time per scan, in minutes — parallel to ``scans``.
    rt: tuple[float, ...] = ()
    #: summed centroid intensity per scan — parallel to ``scans``.
    intensity: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        if not (len(self.scans) == len(self.rt) == len(self.intensity)):
            raise ValueError(
                "Chromatogram scans, rt, and intensity must be equal length "
                f"(got {len(self.scans)}, {len(self.rt)}, {len(self.intensity)})"
            )


@dataclass(frozen=True, slots=True)
class IsotopomerPeak:
    """The integrated result for a single isotopomer of one peptide.

    ``area`` is the integrated intensity reported today as the ``isoN`` column.
    The mass-accuracy fields are M3 additions (PROJECT_REVIEW.md, "mass-accuracy
    output"): they are populated by Week 3 mass-domain refinement and stay
    ``None`` under the legacy fixed-window path, so a missing value is an honest
    "not computed" rather than a fabricated zero.
    """

    #: 0-based isotopomer index.
    isotopomer: int
    #: integrated intensity (the ``isoN`` column).
    area: float
    #: observed (intensity-weighted) m/z — ``isoN_obs_mz``.
    obs_mz: float | None = None
    #: (observed − theoretical) m/z error in ppm — ``isoN_ppm_error``.
    ppm_error: float | None = None
    #: signal-to-noise ratio from the baseline residual — ``isoN_snr``.
    snr: float | None = None
    #: composite peak-quality score (S/N + symmetry) — ``isoN_quality``.
    quality: float | None = None
