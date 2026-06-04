# -*- coding: utf-8 -*-

""" Frozen configuration dataclasses for the M3 (1.0.0) pipeline.

One typed config object per subcommand, intended as the single source of truth
shared by the CLI and (M4) GUI so the two surfaces cannot drift — the 0.9.0 code
duplicates validation across ``main.py`` and ``riana_ui/`` with non-identical
types (PROJECT_REVIEW.md §2d, §4.2).

These are *skeleton* definitions (M3 Week 1). Field defaults mirror the 0.9.0
``riana integrate`` / ``riana fit`` argparse defaults; each field is annotated
with its CLI flag. The constructors that build these from parsed CLI args land
with ``cli.py`` in Week 4 — deliberately not here, since Week 4 replaces argparse
with typer/click.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class IntegrationConfig:
    """Configuration for ``riana integrate``.

    The mass-window semantic is the 0.9.0 one: ``mass_tol_ppm`` is the ±N ppm
    half-width around the theoretical m/z (window is 2·N ppm wide). The M2 v0.9.0
    baseline was generated at ``mass_tol_ppm=15``; M3 integration benchmarks must
    use the same value or the regression comparison is invalid (M2 notes).
    """

    #: -s / --sample. Must end in a number (encodes the time point).
    sample: str = "time0"
    #: -i / --iso. Isotopomers to integrate; 0.9.0 default is the m0/m6 pair.
    isotopomers: tuple[int, ...] = (0, 6)
    #: -m / --mass_tol. ±ppm half-width (see class docstring).
    mass_tol_ppm: int = 50
    #: -r / --r_time. RT tolerance in minutes, applied in both directions.
    r_time: float = 1.0
    #: -q / --q_value. Integrate only PSMs with q-value below this.
    q_value: float = 1e-2
    #: -u / --unique. Restrict to peptides mapping to a single protein.
    unique_only: bool = False
    #: -w / --write_intensities. Also emit the pre-integration intensity trace.
    write_intensities: bool = False
    #: -S / --smoothing. Savitzky-Golay window size; None disables smoothing.
    smoothing: int | None = None
    #: -D / --mass_difference. Mass step between isotopomers (C13 default).
    mass_difference: float = 1.003354835
    #: -X / --ignored_mods. Modification masses excluded from the peptide mass.
    ignored_mods: tuple[float, ...] = ()
    #: -F / --forced_mods. Modification masses always added (e.g. SILAC).
    forced_mods: tuple[float, ...] = (0.0,)
    #: -t / --thread. Worker count for the per-peptide integration map.
    threads: int = 1
    #: -o / --out. Output directory.
    out_dir: str = "."
    #: No CLI flag in 0.9.0 — the pipeline hardcodes ``use_range=True`` (span the
    #: RT window across all PSM scans of a peptide-charge). Exposed as a field so
    #: the Week 3 rewrite can make it configurable without a signature change.
    use_range: bool = True

    # --- M3 Week 3 (Phase C/D) science fields. Defaults are the new
    # canonical pipeline; the legacy fixed-window path stays available via
    # ``peak_method="fixed_window"`` for Phase A parity / regression gating.

    #: Peak boundary method on the iso0 XIC.  ``"fixed_window"`` is the
    #: current default — it reproduces the 0.9.0 rectangle integration.
    #: ``"detected"`` runs :func:`algorithms.peaks.detect_peak` +
    #: co-elution check, falling back per-peptide to fixed-window when
    #: detection fails.
    #:
    #: **Why fixed_window is still the default in Week 3**: the Phase C
    #: regression bench across ac16/ipsc/cm showed that scipy.signal.find_peaks
    #: + peak_widths produces narrow per-peptide-per-fraction boundaries that
    #: are *not stable across the 9 D₂O proportions* — iso5 area drops ~70%
    #: under detection because the boundary lands at slightly different scan
    #: indices at different proportions, breaking R²>0.95 curation (curated
    #: peptide count fell 24–33% across lines; uncurated m0_rmse worsened
    #: ~65%). Baseline tuning (linear → noise-floor) didn't recover.
    #: Conclusion: opt-in only until a cross-proportion-stable peak picker
    #: lands (e.g. consensus boundaries fitted once per peptide, Skyline-
    #: grade noise-aware detection).
    peak_method: str = "fixed_window"
    #: Chromatographic baseline subtraction inside the integrated window.
    #: ``"none"`` matches 0.9.0; ``"noise_floor"`` is the recommended
    #: default — flat p10-of-trace baseline, robust to noise spikes at
    #: the boundary endpoints (the failure mode that broke the first cut
    #: of Phase C; see PROJECT_REVIEW.md Week 3 notes). ``"linear"`` is
    #: the Skyline classic (brittle on our boundaries); ``"snip"`` /
    #: ``"asls"`` are pybaselines fallbacks.
    baseline_method: str = "noise_floor"
    #: Savitzky–Golay polynomial order. ``2`` is the §2c fix (the 0.9.0
    #: polyorder=1 path is mathematically a moving average); only used when
    #: ``smoothing`` is set.
    smoothing_polyorder: int = 2
    #: Per-fraction calibration drift alert threshold, in ppm. Phase D logs
    #: a warning when the per-fraction median ppm error exceeds this.
    ppm_alert: float = 20.0

    def __post_init__(self) -> None:
        if not 1 <= self.mass_tol_ppm <= 500:
            raise ValueError(f"mass_tol_ppm must be in [1, 500], got {self.mass_tol_ppm}")
        if not 0.0 <= self.q_value <= 1.0:
            raise ValueError(f"q_value must be in [0, 1], got {self.q_value}")
        if self.smoothing is not None and (self.smoothing < 3 or self.smoothing % 2 == 0):
            raise ValueError(f"smoothing must be an odd integer >= 3, got {self.smoothing}")
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads}")
        if self.peak_method not in ("detected", "fixed_window"):
            raise ValueError(
                f"peak_method must be 'detected' or 'fixed_window', got {self.peak_method!r}"
            )
        if self.baseline_method not in ("none", "noise_floor", "linear", "snip", "asls"):
            raise ValueError(
                f"baseline_method must be one of 'none', 'noise_floor', 'linear', "
                f"'snip', 'asls'; got {self.baseline_method!r}"
            )
        if self.smoothing_polyorder < 2:
            raise ValueError(
                f"smoothing_polyorder must be >= 2 (the §2c fix); got {self.smoothing_polyorder}"
            )
        if self.ppm_alert <= 0:
            raise ValueError(f"ppm_alert must be > 0, got {self.ppm_alert}")


@dataclass(frozen=True, slots=True)
class FitConfig:
    """Configuration for ``riana fit``.

    Most exposed to change in M3: Week 4 rewrites ``core/fitting.py`` to apply the
    §2b science fixes — AA ``a_max`` dispatch, FS-denominator drift, bootstrap
    kinetic-fit CIs, and replacing the m0/mA-analytic FS calculation with the
    IsoSpec forward/solve model. Treat this as a living definition until then.
    """

    #: -m / --model. One of "simple", "guan", "fornasiero".
    model: str = "simple"
    #: -l / --label. 1=2H in vivo, 2=2H in vitro, 3=18O, 4=amino-acid labeling.
    label: int = 1
    #: -a / --aa. Label-carrying residue(s), for label=4 (e.g. "K", "KR").
    aa: str = "K"
    #: --kp. Precursor rate constant for the two-compartment models.
    k_p: float = 0.5
    #: --kr. Reutilization rate constant for the Fornasiero model.
    k_r: float = 0.05
    #: --rp. Bound/free precursor ratio for the Fornasiero model.
    r_p: float = 10.0
    #: -q / --q_value. Fit only data points with q-value below this.
    q_value: float = 1e-2
    #: -d / --depth. Fit only peptides seen in at least this many samples.
    depth: int = 3
    #: -r / --ria. Precursor enrichment level (RIA max) — the asymptotic
    #: D₂O fraction in body water / culture media (e.g. 0.06 ≈ 6% v/v
    #: D₂O). Used by the IsoSpec forward model
    #: (:func:`algorithms.isotope_dist.solve_fs_d2o`) and by the legacy
    #: analytic FS path (``core/fsynthesis.py``). User-controlled per
    #: experiment because metabolic-water dilution and protocol
    #: variability push this around.
    ria_max: float = 0.06
    #: -f / --fs. Fine-structure FS formula; None uses the m0 analytic path.
    fs_formula: str | None = None
    #: -p / --plotcurves. Emit per-peptide fitted-curve plots.
    plot_curves: bool = False
    #: -t / --thread. Worker count for the per-peptide fit map.
    threads: int = 1
    #: -o / --out. Output directory.
    out_dir: str = "."

    def __post_init__(self) -> None:
        if self.model not in ("simple", "guan", "fornasiero"):
            raise ValueError(f"model must be simple/guan/fornasiero, got {self.model!r}")
        if self.label not in (1, 2, 3, 4):
            raise ValueError(f"label must be 1, 2, 3, or 4, got {self.label}")
        if not 0.0 <= self.q_value <= 1.0:
            raise ValueError(f"q_value must be in [0, 1], got {self.q_value}")
        if self.depth < 1:
            raise ValueError(f"depth must be >= 1, got {self.depth}")
