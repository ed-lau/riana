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
    half-width around the theoretical m/z (window is 2·N ppm wide). Default
    **10 ppm** — the centroid mass-accuracy / search-tolerance norm (the value the
    calibration + animal SDRFs carry, and what both the LVE turnover and the
    ac16/iPSC/CM calibration were re-validated at). The integration window should
    track *mass accuracy* on **centroid** data — one centroid line per isotopomer,
    spread only by calibration drift (~3–10 ppm) — i.e. ≈ the search tolerance,
    NOT a profile *peak width* (~30 ppm at 700 m/z on a 60k Orbitrap). The earlier
    50-ppm default conflated the two: on centroid mzML it imported co-eluting
    near-isobar interference into the heavy isotopomer channels (out-of-range θ
    27%→8% and R²med 0.69→0.86 on LVE when tightened 50→10). On the SDRF path the
    value is read from ``comment[precursor mass tolerance]`` (the search
    tolerance); a CLI ``--mass_tol`` overrides; this dataclass default is the
    last-resort fallback. (The committed v0.9.0 baselines used 15 ppm; the
    integration_port parity tests pin that explicitly.)
    """

    #: -s / --sample. Must end in a number (encodes the time point).
    sample: str = "time0"
    #: -i / --iso. Isotopomers to integrate. Default is the contiguous m0-m5
    #: envelope the D2O fit consumes (``solve_fs_d2o`` matches the observed
    #: envelope against the IsoSpec forward model over these channels); the
    #: legacy ``(0, 6)`` m0/m6 pair is non-contiguous and is NOT fittable by the
    #: new engine (see ``core.fitting._REQUIRED_D2O_ISOTOPOMERS``). Override for
    #: other workflows (e.g. SILAC cluster extraction via ``forced_mods``).
    isotopomers: tuple[int, ...] = (0, 1, 2, 3, 4, 5)
    #: -m / --mass_tol. ±ppm half-width (see class docstring). Default 10 ppm
    #: (centroid mass-accuracy / search-tolerance norm); SDRF
    #: ``comment[precursor mass tolerance]`` supplies it on the SDRF path, CLI
    #: overrides.
    mass_tol_ppm: int = 10
    #: -r / --extraction_half_width (legacy alias --r_time). EXTRACTION
    #: half-width in RT minutes (how much XIC to pull, both directions). Must
    #: be ≥ ``integration_half_width`` (plus room for the apex offset when
    #: ``peak_rt="apex"``). For the ms2/fixed path the runner sets it equal to
    #: ``integration_half_width`` (0.9.0 parity). Renamed from ``r_time`` —
    #: that name conflated extraction with integration; they're distinct now.
    #: Default 0.5 suits the apex default (covers the apex offset + a 0.15
    #: window); the runner/CLI derive it from ``integration_half_width``.
    extraction_half_width: float = 0.5
    #: -q / --q_value. Integrate only PSMs with q-value below this.
    q_value: float = 1e-2
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

    # --- M3 peak-integration fields. DEFAULT since the 2026-06 peak-detection
    # spike: an apex-centred narrow window (peak_rt="apex",
    # integration_half_width=0.15, apex_selection="tallest", baseline none),
    # which beat the 0.9.0 fixed rectangle on every calibration line (ac16/
    # ipsc/cm, 0% & 100%) and an in-vivo mouse set. ms2 + whole-window
    # reproduces 0.9.0 (parity/regression escape hatch); consensus is the
    # high-D₂O alternative; auto detects boundaries. See PROJECT_REVIEW §2c.

    #: WHERE the integration window is anchored. **Default ``"apex"``** (spike
    #: winner): locate the chromatographic apex on the iso0 XIC
    #: (:func:`algorithms.peaks.find_apex`; prominence-gated, ``apex_selection``)
    #: and integrate apex ± ``integration_half_width``. ``"ms2"`` reproduces
    #: 0.9.0 — centre on the PSM/MS2 RT and integrate the whole
    #: ±``extraction_half_width`` (parity escape hatch). ``"consensus"`` takes
    #: the **median** apex over m0..m{``apex_n_consensus``-1}
    #: (:func:`algorithms.peaks.consensus_apex`) — labelling-independent and
    #: contamination-robust; a marginal edge at high D₂O, so prefer it there.
    #: ``integration_half_width="auto"`` detects the boundaries (implies apex).
    peak_rt: str = "apex"
    #: HOW WIDE to integrate — a float half-width in RT minutes (window =
    #: centre ± this), or ``"auto"`` to detect boundaries from the iso0 shape
    #: (``peak_widths`` at ``width_rel_height``; implies apex). **Default 0.15**
    #: (the spike optimum across lines); the best value tracks the
    #: chromatographic peak width — **dial it to your gradient** (≈0.1 for sharp
    #: UPLC, ≈0.2–0.33 for broad peaks). In the ms2 path the window IS the whole
    #: ±``extraction_half_width`` (runner sets it = this for 0.9.0 parity); for
    #: apex/consensus the window is a sub-interval of the wider extraction.
    integration_half_width: float | str = 0.15
    #: Chromatographic baseline subtraction inside the integrated window.
    #: ``"none"`` (default) — the narrow apex window removes background by
    #: *exclusion*, which beat subtraction across all three lines in the spike
    #: (see PROJECT_REVIEW). ``"noise_floor"`` — flat p10-of-(full-trace)
    #: baseline; competitive but ~no-op at low labelling, kept for future
    #: tuning (its quality depends on the extraction width — see the
    #: baseline.py off-peak TODO). ``"snip"`` / ``"asls"`` are pybaselines
    #: fallbacks. (Skyline-style ``"linear"`` was tested and discarded — it
    #: over-subtracts on our narrow on-peak boundaries; see PROJECT_REVIEW.)
    baseline_method: str = "none"
    #: Savitzky–Golay polynomial order. ``2`` is the §2c fix (the 0.9.0
    #: polyorder=1 path is mathematically a moving average); only used when
    #: ``smoothing`` is set.
    smoothing_polyorder: int = 2
    #: Per-fraction calibration drift alert threshold, in ppm. Phase D logs
    #: a warning when the per-fraction median ppm error exceeds this.
    ppm_alert: float = 20.0
    #: Local-MAD prominence multiplier for the apex finder
    #: (:func:`algorithms.peaks.find_apex` / ``detect_peak``). A candidate apex
    #: must clear ``prominence_k · 1.4826 · MAD(trace)``; higher ⇒ stricter.
    #: Used by ``peak_rt="apex"`` and ``integration_half_width="auto"``.
    prominence_k: float = 3.0
    #: Apex-height fraction at which ``integration_half_width="auto"`` measures
    #: the peak width (via ``scipy.signal.peak_widths``). 0.05 = 5% of apex
    #: (Skyline-classic; boundary down in the noisy flanks → wanders across
    #: proportions). 0.5 = FWHM (on the steep near-apex flank → far more
    #: stable, and amplitude-robust so iso0 suppression at high D₂O doesn't
    #: move it). Only used when ``integration_half_width="auto"``.
    width_rel_height: float = 0.05
    #: Apex selection rule for ``peak_rt="apex"`` / ``"consensus"``.
    #: **Default ``"tallest"``** (most intense prominent candidate — won the
    #: cross-dataset screen); ``"nearest"`` picks the candidate closest to the
    #: MS2 RT prior (better when a tall co-eluting neighbour is a risk).
    apex_selection: str = "tallest"
    #: Half-width (RT min) bounding the apex search around the PSM RT prior
    #: (the best-q PSM's scan for a peptide-charge with multiple PSMs — see
    #: :func:`riana.core.integration.integrate_run`). ``0`` = the whole
    #: extracted trace, which (with ``use_range=True`` spanning every PSM scan)
    #: lets the apex roam to the tallest peak anywhere in the window and grab a
    #: co-eluting isobar — the failure mode found on the animal D₂O series.
    #: **Default ``0.25``**: keep the apex within ±0.25 min of the confident ID,
    #: so the integration window (apex ± ``integration_half_width``) still nests
    #: inside the extraction (``extraction_half_width`` ≥ 0.25 + 0.15 = 0.40).
    apex_search_half_width: float = 0.25
    #: Number of leading isotopomer channels (m0..m{n-1}) the
    #: ``peak_rt="consensus"`` median-apex pools over (co-elution consensus).
    apex_n_consensus: int = 4

    # --- Intake scan↔RT guard (Track A). -------------------------------------
    #: When true (default), :func:`riana.core.integration.integrate_run` verifies
    #: per run that each PSM's ``spectra_ref`` scan → this mzML's MS1 RT
    #: reconciles with the mzTab-reported ``retention_time``, and **errors** if
    #: the per-run *median* offset exceeds :attr:`scan_rt_tol_min`. This catches
    #: the quantms filename-prefix scan-scramble (mzML basenames that are
    #: prefixes of one another) and wrong mzML↔mzTab pairings — a class of bug
    #: that was previously silent (PROJECT_REVIEW Track A). No-ops on the
    #: Percolator path (no ``retention_time``). Disable with ``--no-rt-check``
    #: only for a run you know is correctly paired.
    check_scan_rt: bool = True
    #: Median scan↔RT offset (RT minutes) above which :attr:`check_scan_rt`
    #: errors. **Default 2.0** — clears the run-dependent ProteomicsLFQ alignment
    #: offset (≤~0.9 min measured on real output) with margin, while a
    #: scan-scrambled run sits tens of minutes off (~25× the threshold).
    scan_rt_tol_min: float = 2.0

    def __post_init__(self) -> None:
        if not 1 <= self.mass_tol_ppm <= 500:
            raise ValueError(f"mass_tol_ppm must be in [1, 500], got {self.mass_tol_ppm}")
        if not 0.0 <= self.q_value <= 1.0:
            raise ValueError(f"q_value must be in [0, 1], got {self.q_value}")
        if self.smoothing is not None and (self.smoothing < 3 or self.smoothing % 2 == 0):
            raise ValueError(f"smoothing must be an odd integer >= 3, got {self.smoothing}")
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads}")
        if self.peak_rt not in ("ms2", "apex", "consensus"):
            raise ValueError(
                f"peak_rt must be 'ms2', 'apex' or 'consensus', got {self.peak_rt!r}")
        if self.apex_selection not in ("nearest", "tallest"):
            raise ValueError(
                f"apex_selection must be 'nearest' or 'tallest', got {self.apex_selection!r}")
        if self.apex_search_half_width < 0:
            raise ValueError(
                f"apex_search_half_width must be >= 0, got {self.apex_search_half_width}")
        if self.apex_n_consensus < 1:
            raise ValueError(
                f"apex_n_consensus must be >= 1, got {self.apex_n_consensus}")
        if self.integration_half_width != "auto":
            w = self.integration_half_width
            if isinstance(w, bool) or not isinstance(w, (int, float)) or w <= 0:
                raise ValueError(
                    "integration_half_width must be a positive number or 'auto', "
                    f"got {self.integration_half_width!r}"
                )
        if self.prominence_k <= 0:
            raise ValueError(f"prominence_k must be > 0, got {self.prominence_k}")
        if not 0.0 < self.width_rel_height < 1.0:
            raise ValueError(
                f"width_rel_height must be in (0, 1), got {self.width_rel_height}"
            )
        if self.baseline_method not in ("none", "noise_floor", "snip", "asls"):
            raise ValueError(
                f"baseline_method must be one of 'none', 'noise_floor', "
                f"'snip', 'asls'; got {self.baseline_method!r}"
            )
        if self.smoothing_polyorder < 2:
            raise ValueError(
                f"smoothing_polyorder must be >= 2 (the §2c fix); got {self.smoothing_polyorder}"
            )
        if self.ppm_alert <= 0:
            raise ValueError(f"ppm_alert must be > 0, got {self.ppm_alert}")
        if self.scan_rt_tol_min <= 0:
            raise ValueError(
                f"scan_rt_tol_min must be > 0, got {self.scan_rt_tol_min}")


@dataclass(frozen=True, slots=True)
class FitConfig:
    """Configuration for ``riana fit``.

    Heavy-water (D₂O) fitting via the IsoSpec forward/solve model: per-peptide
    Spep from a ``--coefficients`` table, per-timepoint FS by full-envelope
    least-squares, bootstrap kinetic-fit CIs (the M3 Week 4 science fixes). The
    per-AA coefficient table is the single source of cell/tissue specificity —
    ``label`` only selects the labeling chemistry (``hw``; ``o18`` is a
    post-M4 placeholder).
    """

    #: -m / --model. One of "simple", "guan", "fornasiero".
    model: str = "simple"
    #: -l / --label. The labeling chemistry. ``"hw"`` (heavy water / D₂O,
    #: default) is the only path the M4 fit engine implements — cell/tissue
    #: specificity comes from the per-AA ``--coefficients`` table, not the
    #: label. ``"o18"`` (¹⁸O) is recognized but its fit path is being
    #: reimplemented post-M4 (``fit_run`` raises a clear error). Amino-acid /
    #: SILAC labeling was dropped from fitting — integrate still extracts SILAC
    #: peaks via ``forced_mods``; do the L/(H+L) curve fit downstream.
    label: str = "hw"
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
    #: -f / --fs. Reserved (post-M4): restrict the envelope SSE in
    #: :func:`algorithms.isotope_dist.solve_fs_d2o` to a subset of isotopomer
    #: channels (e.g. m0-m2) to reduce sensitivity to co-eluting contaminants
    #: in the higher isotopomers. Currently **ignored** — the full integrated
    #: envelope is used. Kept on the config so the post-M4 wiring is a localized
    #: change. (Not the legacy fine-structure-ratio meaning.)
    fs_formula: str | None = None
    #: -t / --thread. Worker count for the per-peptide fit map.
    threads: int = 1
    #: -o / --out. Output directory.
    out_dir: str = "."

    def __post_init__(self) -> None:
        if self.model not in ("simple", "guan", "fornasiero"):
            raise ValueError(f"model must be simple/guan/fornasiero, got {self.model!r}")
        if self.label not in ("hw", "o18"):
            raise ValueError(f"label must be 'hw' or 'o18', got {self.label!r}")
        if not 0.0 <= self.q_value <= 1.0:
            raise ValueError(f"q_value must be in [0, 1], got {self.q_value}")
        if self.depth < 1:
            raise ValueError(f"depth must be >= 1, got {self.depth}")
