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

    def __post_init__(self) -> None:
        if not 1 <= self.mass_tol_ppm <= 500:
            raise ValueError(f"mass_tol_ppm must be in [1, 500], got {self.mass_tol_ppm}")
        if not 0.0 <= self.q_value <= 1.0:
            raise ValueError(f"q_value must be in [0, 1], got {self.q_value}")
        if self.smoothing is not None and (self.smoothing < 3 or self.smoothing % 2 == 0):
            raise ValueError(f"smoothing must be an odd integer >= 3, got {self.smoothing}")
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads}")


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
    #: -r / --ria. Final isotope enrichment level (RIA max).
    ria_max: float = 0.5
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
