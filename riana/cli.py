# -*- coding: utf-8 -*-

"""Riana command-line interface (Typer).

This is the M4 CLI rewrite. It replaces the 0.9.0 argparse ``main.py`` and the
``--engine legacy``/``--engine new`` split: the typed pipeline
(``riana.core`` + ``riana.io``, driven by the frozen
:class:`riana.config.IntegrationConfig` / :class:`~riana.config.FitConfig`) is
now the *only* engine. Those frozen dataclasses are the single source of truth
shared with the (M4 Phase 2) Qt GUI — their ``__post_init__`` carries the domain
validation, so both surfaces validate identically (PROJECT_REVIEW §2d, §4.2).

The CLI layer here only does the *CLI-shaped* checks Typer/click can't express
on the dataclass (paths exist, thread ≤ cpu count, sample ends in a digit) and
the arg→config marshalling; everything numeric is the dataclass's job.

Defaults reproduce the 2026-06 peak-detection spike winner (apex-centred narrow
window). Reproduce 0.9.0 integration with ``--peak-rt ms2
--integration-half-width 1.0``.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import List, Optional

import typer

from riana import __version__

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Riana — integrate isotopomer abundances from MS1 data and fit "
    "protein-turnover kinetics. https://github.com/ed-lau/riana",
)

# Rich-help panels that segregate the power-user dials from the everyday options
# in ``riana integrate --help`` (the GUI mirrors this with a collapsed "Advanced"
# group). The everyday flags stay in Typer's default "Options" panel.
_ADV = "Advanced integration (tuning dials)"
_MBR = "Match-between-runs (MBR)"


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"riana {__version__}")
        raise typer.Exit()


def _parse_number_list(value: str | None, cast, *, sort: bool, unique: bool):
    """Parse a comma/space-separated numeric list (``"0 6"`` or ``"0,6"``).

    Replaces the 0.9.0 argparse ``nargs='+'`` (``-i 0 6``); a single token keeps
    the flag Snakemake- and shell-quote-friendly.
    """
    if value is None or not value.strip():
        return ()
    parts = [p for p in re.split(r"[,\s]+", value.strip()) if p]
    out = [cast(p) for p in parts]
    if unique and len(out) != len(set(out)):
        raise typer.BadParameter(f"duplicate values not allowed: {value!r}")
    if sort:
        out.sort()
    return tuple(out)


@app.callback()
def _main(
    version: Optional[bool] = typer.Option(
        None, "--version", "-v",
        callback=_version_callback, is_eager=True,
        help="Show the riana version and exit.",
    ),
) -> None:
    """Riana CLI root."""


# --------------------------------------------------------------------------- #
# integrate
# --------------------------------------------------------------------------- #
@app.command()
def integrate(
    mzml_path: Path = typer.Argument(
        ..., exists=True, file_okay=False, dir_okay=True, readable=True,
        help="Folder containing the mzML file(s).",
    ),
    id_path: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True,
        help="The search-ID file: a Percolator target psms.txt (single-mzML "
        "path), or — with --sdrf — the quantms mzTab (DDA) / DIA-NN "
        "report.parquet (DIA) covering every run.",
    ),
    sdrf: Optional[Path] = typer.Option(
        None, "--sdrf", exists=True, dir_okay=False, readable=True,
        help="SDRF samplesheet (the primary path). Drives identity-keyed "
        "intake: id_path is read as an mzTab (DDA) or DIA-NN report.parquet "
        "(when the SDRF declares DIA), one <mzml_stem>_riana.txt is written per "
        "run with its full identity in the header, and a riana_manifest.tsv is "
        "written/updated. Without --sdrf, id_path is a Percolator file (the "
        "demoted single-mzML tier).",
    ),
    sample: str = typer.Option(
        "time0", "-s", "--sample",
        help="Sample name for the Percolator (no-SDRF) path; must end with a "
        "number encoding the time point, e.g. time1. Ignored with --sdrf "
        "(identity comes from the SDRF).",
    ),
    iso: str = typer.Option(
        "5", "-i", "--iso",
        help="Isotopomers to integrate. Give a single index N for the m0..mN "
        "envelope (default '5' = m0-m5, what the D2O fit consumes). An explicit "
        "comma/space list is still accepted for a non-contiguous set (e.g. the "
        "o18 '0 6' pair). Pass 'auto' for adaptive N_ISO: per-peptide channels "
        "from the IsoSpec init+final envelope (≥1% abundance), extracted at "
        "averaged-isotopolog accurate mass, padded to the run-wide max width.",
    ),
    ria: Optional[float] = typer.Option(
        None, "--ria", metavar="FRAC",
        help="Precursor enrichment (RIA max), e.g. 0.06 for 6% v/v D2O. Shapes the "
        "final envelope under --iso auto. If omitted, taken from the SDRF "
        "characteristics[precursor enrichment], else the 0.06 default; given here, "
        "overrides the SDRF.",
    ),
    out: Path = typer.Option(
        Path("."), "-o", "--out", help="Output directory [default: .]."),
    q_value: float = typer.Option(
        1e-2, "-q", "--q_value", metavar="FDR",
        help="Integrate only PSMs with q-value below this [default: 1e-2]."),
    extraction_half_width: Optional[float] = typer.Option(
        None, "-r", "--extraction_half_width", "--r_time", rich_help_panel=_ADV,
        help="Extraction half-width (RT min, both directions): how much XIC to "
        "pull around the PSM. If omitted, derived from --integration-half-width "
        "/ --peak-rt. (alias: --r_time)",
    ),
    peak_rt: str = typer.Option(
        "apex", "--peak-rt",
        help="Window anchor: 'apex' (default, spike winner), 'ms2' (0.9.0 "
        "parity), or 'consensus' (median apex over m0..m3; prefer at high D2O).",
    ),
    integration_half_width: str = typer.Option(
        "0.15", "--integration-half-width", metavar="MIN|auto",
        help="Integration half-width in RT min (default 0.15; dial to your "
        "chromatographic peak width), or 'auto' to detect boundaries.",
    ),
    baseline_method: str = typer.Option(
        "none", "--baseline", rich_help_panel=_ADV,
        help="In-window baseline subtraction: none (default), noise_floor, "
        "snip, asls.",
    ),
    apex_selection: str = typer.Option(
        "tallest", "--apex-selection", rich_help_panel=_ADV,
        help="Apex pick rule for apex/consensus: tallest (default) or nearest.",
    ),
    apex_search_half_width: float = typer.Option(
        0.25, "--apex-search-half-width", metavar="MIN", rich_help_panel=_ADV,
        help="Half-width (RT min) bounding the apex search around the PSM/MBR RT "
        "prior (default 0.25; 0 = whole extraction). Tighter keeps the apex on the "
        "confident RT anchor — relevant for cross-proportion / MBR stability.",
    ),
    write_intensities: bool = typer.Option(
        False, "-w", "--write_intensities", rich_help_panel=_ADV,
        help="Also write the pre-integration intensity trace."),
    mass_tol: Optional[int] = typer.Option(
        None, "-m", "--mass_tol", metavar="PPM",
        help="Mass tolerance half-width in ppm (±N ppm). If omitted, taken from "
        "the SDRF comment[precursor mass tolerance], else the 10 ppm default. "
        "Given here, overrides the SDRF."),
    smoothing: Optional[int] = typer.Option(
        None, "-S", "--smoothing", rich_help_panel=_ADV,
        help="Savitzky-Golay smoothing window (odd int ≥ 3)."),
    smoothing_polyorder: int = typer.Option(
        2, "--smoothing-polyorder", metavar="ORDER", rich_help_panel=_ADV,
        help="Savitzky-Golay polynomial order; only used with --smoothing "
        "[default: 2]."),
    mass_difference: float = typer.Option(
        1.003354835, "-D", "--mass_difference", rich_help_panel=_ADV,
        help="Mass difference between isotopomers [default: 1.003354835]."),
    ppm_alert: float = typer.Option(
        20.0, "--ppm-alert", metavar="PPM", rich_help_panel=_ADV,
        help="Per-fraction drift-alert threshold: warn when the median ppm error "
        "exceeds this [default: 20]."),
    prominence_k: float = typer.Option(
        3.0, "--prominence-k", metavar="K", rich_help_panel=_ADV,
        help="Apex-finder strictness — a candidate apex must clear "
        "K·1.4826·MAD(trace); higher is stricter. Used by --peak-rt apex/consensus "
        "and --integration-half-width auto [default: 3.0]."),
    width_rel_height: float = typer.Option(
        0.05, "--width-rel-height", metavar="FRAC", rich_help_panel=_ADV,
        help="Apex-height fraction at which --integration-half-width auto measures "
        "the peak width: 0.05 = 5% of apex (Skyline-classic), 0.5 = FWHM (more "
        "stable). Only used with --integration-half-width auto [default: 0.05]."),
    apex_n_consensus: int = typer.Option(
        4, "--apex-n-consensus", metavar="N", rich_help_panel=_ADV,
        help="Leading isotopomer channels (m0..m{N-1}) the consensus apex pools "
        "over. Only used with --peak-rt consensus [default: 4]."),
    workers: int = typer.Option(
        1, "-W", "--workers", metavar="N",
        help="Runs to integrate concurrently on the --sdrf path (one mzML in "
        "memory per worker; 2-4 suits a many-timepoint time series) [default: "
        "1]. The parallelism lever — per-run extraction is GIL-bound serial."),
    no_id_check: bool = typer.Option(
        False, "--no-id-check", rich_help_panel=_ADV,
        help="Disable the intake scan↔precursor guard — the per-run check that "
        "the mzTab spectra_ref scans point at the matching precursor m/z in this "
        "mzML (catches a wrong mzML↔mzTab pairing / quantms filename-prefix "
        "scan-scramble; mass-based, so immune to OpenMS RT alignment). Only "
        "disable for a run you know is correctly paired."),
    precursor_tol_ppm: float = typer.Option(
        10.0, "--precursor-tol-ppm", metavar="PPM", rich_help_panel=_ADV,
        help="Per-scan precursor-m/z match tolerance (ppm) for the intake "
        "scan↔precursor guard [default: 10.0]. Lenient to precursor-refinement "
        "drift; a wrong file lands hundreds of ppm off. Ignored with --no-id-check."),
    resume: bool = typer.Option(
        False, "--resume",
        help="On the --sdrf path, skip runs already integrated in the output's "
        "manifest at the current settings (each run's <stem>_riana.txt is written "
        "as it finishes, so a crashed run resumes where it stopped). Assumes the "
        "same SDRF/mzTab inputs."),
    mbr: bool = typer.Option(
        False, "--mbr", rich_help_panel=_MBR,
        help="Match-between-runs (mzTab/DDA, --sdrf path): transfer a confidently "
        "identified precursor into the runs of its (experiment, condition) curve "
        "that missed it, recovering points lost to stochastic MS2 sampling. "
        "Transferred rows are flagged evidence='mbr'; one with no detectable apex "
        "is dropped. No-op on DIA (DIA-NN already propagates)."),
    mbr_min_donor_runs: int = typer.Option(
        2, "--mbr-min-donor-runs", metavar="N", rich_help_panel=_MBR,
        help="MBR donor gate: transfer a precursor only if it is confidently "
        "identified in at least N runs of the curve group [default: 2]."),
    mbr_donor_q: float = typer.Option(
        1e-2, "--mbr-donor-q", metavar="Q", rich_help_panel=_MBR,
        help="MBR donor-confidence q-value threshold [default: 0.01]."),
    mbr_min_snr: float = typer.Option(
        4.0, "--mbr-min-snr", metavar="SNR", rich_help_panel=_MBR,
        help="Drop MBR transfers whose apex SNR (prominence/local-noise) is below "
        "this (inf SNR — a sparse MAD=0 trace — also fails). Default 4 (0 = ungated): "
        "the fit A/B showed ungated MBR is harmful but SNR≥4 makes it net-positive "
        "at in-vivo R² gates with no clean-curve pollution."),
    mbr_min_scans: int = typer.Option(
        3, "--mbr-min-scans", metavar="N", rich_help_panel=_MBR,
        help="Drop MBR transfers with fewer than N nonzero scans in the integration "
        "window (a sparse XIC can't define a reliable peak). Default 3 (0 = off); "
        "keeps ~98% of real-quality peaks. The companion to --mbr-min-snr."),
) -> None:
    """Integrate isotopomer abundance over retention time."""
    import dataclasses
    import json

    from riana.config import IntegrationConfig
    from riana.core.integration import integrate_run
    from riana.exceptions import DataError
    from riana.io.mzml import IndexedMzML, list_mzml_files, mzml_stem
    from riana.io.percolator import file_indices, fraction_psms, read_percolator
    from riana.io.writers import make_provenance, write_dataframe_tsv
    from riana.logger import get_logger

    # --- CLI-shaped validation (the dataclass does the numeric domain checks) -
    if workers > (os.cpu_count() or 1):
        raise typer.BadParameter(
            f"--workers {workers} exceeds CPU count ({os.cpu_count()}).")
    # The --sample digit convention is only how the no-SDRF Percolator path
    # encodes the timepoint; with --sdrf the timepoint is an SDRF column.
    if sdrf is None and (not sample or not sample[-1].isdigit()):
        raise typer.BadParameter(
            f"--sample must end with a number (got {sample!r}).")

    # --iso auto = adaptive N_ISO (per-peptide envelope-driven channels). The
    # explicit `isotopomers` tuple is then only a fallback for the few fixed-path
    # branches; keep the m0-m5 default so anything reading it still sees a sane set.
    adaptive_iso = iso.strip().lower() == "auto"
    if adaptive_iso:
        isotopomers = (0, 1, 2, 3, 4, 5)
    else:
        parsed = _parse_number_list(iso, int, sort=True, unique=True)
        if not parsed:
            raise typer.BadParameter("--iso must be 'auto', a single index N "
                                     "(= iso0..isoN), or an explicit list.")
        # A single int N is the easy form for the usual contiguous m0..mN capture
        # (e.g. '5' = the m0-m5 D2O envelope). A multi-value list stays explicit,
        # so a genuinely non-contiguous set (the o18 '0 6' pair) is preserved.
        isotopomers = (tuple(range(parsed[0] + 1)) if len(parsed) == 1
                       else tuple(parsed))

    ihw: float | str = ("auto" if integration_half_width == "auto"
                        else float(integration_half_width))
    # Derive the extraction half-width unless -r given: ms2 integrates the whole
    # extraction (= ihw); apex/consensus need room for the apex offset (+0.33);
    # 'auto' uses a 0.33 base. (Lifted from the M3 _integrate_new adapter.)
    if extraction_half_width is not None:
        ehw = float(extraction_half_width)
    elif peak_rt == "ms2" and ihw != "auto":
        ehw = float(ihw)
    else:
        ehw = (0.33 if ihw == "auto" else float(ihw)) + 0.33

    os.makedirs(out, exist_ok=True)

    # Read the SDRF up front (primary path) so its precursor mass tolerance can
    # feed the integration window. Resolution: explicit --mass_tol > SDRF
    # comment[precursor mass tolerance] > dataclass default. Mirrors the
    # per-sample RIA resolution (SDRF -> --ria default).
    sdrf_table = None
    if sdrf is not None:
        from riana.io.sdrf import read_sdrf
        try:
            sdrf_table = read_sdrf(sdrf)
        except DataError as exc:
            raise typer.BadParameter(str(exc)) from exc
    _default_mass_tol = IntegrationConfig.__dataclass_fields__["mass_tol_ppm"].default
    if mass_tol is not None:
        mass_tol_ppm, mass_tol_src = int(mass_tol), "--mass_tol"
    elif sdrf_table is not None and sdrf_table.precursor_mass_tol_ppm is not None:
        mass_tol_ppm = int(round(sdrf_table.precursor_mass_tol_ppm))
        mass_tol_src = "SDRF comment[precursor mass tolerance]"
    else:
        mass_tol_ppm, mass_tol_src = _default_mass_tol, "default"

    # Precursor enrichment (RIA max): --ria > SDRF characteristics[precursor
    # enrichment] (experiment-level — physically one value per experiment) > the
    # 0.06 default. Only load-bearing under --iso auto (it shapes the final
    # envelope), but resolved uniformly so the header records it either way.
    _default_ria = IntegrationConfig.__dataclass_fields__["ria_max"].default
    if ria is not None:
        ria_max, ria_src = float(ria), "--ria"
    elif sdrf_table is not None:
        _enrich = sorted({
            round(float(r.precursor_enrichment), 6) for r in sdrf_table.runs
            if r.precursor_enrichment is not None
        })
        if len(_enrich) == 1:
            ria_max, ria_src = _enrich[0], "SDRF characteristics[precursor enrichment]"
        elif len(_enrich) > 1:
            # Physically constant within an experiment; if rows disagree take the
            # max (the conservative — widest-envelope — choice for adaptive N_ISO).
            ria_max, ria_src = _enrich[-1], (
                f"SDRF characteristics[precursor enrichment] (max of {_enrich})")
        else:
            ria_max, ria_src = _default_ria, "default"
    else:
        ria_max, ria_src = _default_ria, "default"

    # The frozen config's __post_init__ owns the numeric domain validation
    # (mass_tol range, q_value range, peak_rt enum, ...); surface it as a clean
    # CLI error rather than a traceback.
    try:
        config = IntegrationConfig(
            sample=sample,
            isotopomers=isotopomers,
            adaptive_iso=adaptive_iso,
            ria_max=ria_max,
            mass_tol_ppm=mass_tol_ppm,
            extraction_half_width=ehw,
            peak_rt=peak_rt,
            integration_half_width=ihw,
            baseline_method=baseline_method,
            apex_selection=apex_selection,
            apex_search_half_width=float(apex_search_half_width),
            apex_n_consensus=int(apex_n_consensus),
            prominence_k=float(prominence_k),
            width_rel_height=float(width_rel_height),
            q_value=float(q_value),
            write_intensities=bool(write_intensities),
            smoothing=smoothing,
            smoothing_polyorder=int(smoothing_polyorder),
            mass_difference=float(mass_difference),
            ppm_alert=float(ppm_alert),
            out_dir=str(out),
            check_scan_id=not no_id_check,
            scan_precursor_tol_ppm=float(precursor_tol_ppm),
            mbr=bool(mbr),
            mbr_min_donor_runs=int(mbr_min_donor_runs),
            mbr_donor_q=float(mbr_donor_q),
            mbr_min_snr=float(mbr_min_snr),
            mbr_min_scans=int(mbr_min_scans),
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    logger = get_logger(__name__, str(out))
    logger.info(f"riana {__version__}")
    logger.info("integrate (typed pipeline)")
    logger.info(f"mass tolerance: ±{mass_tol_ppm} ppm (from {mass_tol_src})")
    if adaptive_iso:
        logger.info(
            f"isotopomers: adaptive (--iso auto, ≥{config.iso_abundance_floor:.0%} "
            f"abundance, RIA max {ria_max:g} from {ria_src})")
    else:
        logger.info(f"isotopomers: {list(isotopomers)} (fixed)")

    # Prime the provenance git-SHA cache now, while this process is still
    # single-threaded — before the integrate worker pool spawns. ``_git_sha``
    # ``fork``s ``git``, and forking from the multi-threaded pool main process
    # deadlocks on macOS (``pthread_atfork``); per-file ``make_provenance`` then
    # reuses the cached value instead of forking once per output file.
    from riana.io.writers import warm_git_sha
    warm_git_sha()

    # --- SDRF path (primary): identity-keyed mzTab intake via the shared
    # pipeline — one <stem>_riana.txt per run + a manifest. ---------------------
    if sdrf_table is not None:
        from riana.core.pipeline import integrate_project

        try:
            logger.info(
                f"SDRF {sdrf}: {len(sdrf_table.runs)} runs, "
                f"{sdrf_table.experiment_type}, {sdrf_table.acquisition}"
            )
            integrate_project(
                config, sdrf_table, mzml_path, id_path, out,
                max_workers=int(workers), resume=bool(resume), logger=logger,
            )
        except DataError as exc:
            raise typer.BadParameter(str(exc)) from exc
        logger.info("integrate: done")
        logger.handlers.clear()
        return

    # --- Percolator path (demoted single-mzML testing/legacy tier). -----------
    all_psms = read_percolator(str(id_path), sample=sample)

    # mzML directory layout: sort by name, accept .mzML / .mzML.gz (mirrors the
    # 0.9.0 fraction-index assignment when no percolator.log.txt is present).
    mzml_files = list_mzml_files(mzml_path)
    if not mzml_files:
        raise typer.BadParameter(f"No mzML files in {mzml_path}.")
    indices = file_indices(all_psms)
    if len(mzml_files) != len(indices):
        raise typer.BadParameter(
            f"mzML count ({len(mzml_files)}) != distinct file_idx count "
            f"({len(indices)}) in {id_path}.")

    for idx in indices:
        mzml_basename = mzml_stem(mzml_files[idx])
        mzml_file = os.path.join(mzml_path, mzml_files[idx])
        logger.info(f"integrating fraction {idx}: {mzml_basename}")

        fraction = fraction_psms(all_psms, idx)
        with IndexedMzML(mzml_file) as mzml:
            df = integrate_run(config, fraction, mzml, file_label=mzml_basename)

        out_file = Path(out) / f"{sample}_riana.txt"
        provenance = make_provenance(
            dataclasses.asdict(config),
            id_source=str(id_path),
            extra={"mzml": mzml_basename},
        )
        write_dataframe_tsv(out_file, df, provenance, include_index=True)

        drift = df.attrs.get("drift_summary")
        if drift is not None:
            drift_path = out_file.with_suffix(".drift.json")
            with drift_path.open("w") as fh:
                json.dump(dataclasses.asdict(drift), fh, indent=2)
        logger.info(f"wrote {out_file} (+ drift sidecar)")

    logger.info("integrate: done")
    logger.handlers.clear()


# --------------------------------------------------------------------------- #
# fit
# --------------------------------------------------------------------------- #
@app.command()
def fit(
    riana_path: Optional[List[Path]] = typer.Argument(
        None, exists=True, dir_okay=False,
        help="Integrate output _riana.txt files (legacy/single-curve path; the "
        "sample field encodes the time point, e.g. time0, time6). Omit when "
        "using --manifest.",
    ),
    manifest: Optional[Path] = typer.Option(
        None, "--manifest", exists=True, dir_okay=False, readable=True,
        help="riana_manifest.tsv from `integrate --sdrf` (the primary path). "
        "Runs are grouped into kinetic curves by (experiment, condition) with "
        "the timepoint taken from the SDRF identity, not the filename.",
    ),
    coefficients: Optional[str] = typer.Option(
        None, "--coefficients",
        help="Labeling-site coefficient table — REQUIRED; the format depends on "
        "--label. hw: a per-AA D2O table (preset deberneh_2025_rss [recommended] | "
        "ilchenko_2019 | commerford_1983 | alamillo_2025_{ac16,ipsc,cm}, or an "
        "(amino_acid, coefficient) CSV). o18: a length-model table (preset "
        "juber_2026_o18_ac16 [in-vitro] | rachdaoui_2009_o18 [in-vivo mouse], "
        "or a (feature, coefficient) CSV).",
    ),
    model: str = typer.Option(
        "simple", "-m", "--model",
        help="Kinetic model: simple (default), guan, fornasiero. Or 'calibration' "
             "— a through-origin FS-vs-mixing-proportion recovery line (R² = "
             "recovery quality); auto-selected on --manifest calibration runs."),
    label: str = typer.Option(
        "hw", "-l", "--label",
        help="Labeling chemistry: 'hw' (heavy water / D2O, default) or 'o18' "
        "(metabolic H2-18O; length-model Spep + 3-isotope-O envelope). Each "
        "needs its matching --coefficients table. (Amino-acid / SILAC fitting "
        "was dropped — integrate SILAC peaks, fit L/(H+L) downstream.)",
    ),
    kp: float = typer.Option(
        0.5, "--kp", help="Precursor rate constant (two-compartment models)."),
    kr: float = typer.Option(
        0.05, "--kr", help="Reutilization rate constant (Fornasiero)."),
    rp: float = typer.Option(
        10.0, "--rp", help="Bound/free precursor ratio (Fornasiero)."),
    q_value: float = typer.Option(
        1e-2, "-q", "--q_value", metavar="FDR",
        help="Fit only data points with q-value below this [default: 1e-2]."),
    depth: int = typer.Option(
        3, "-d", "--depth",
        help="Fit only peptidoforms seen at this many distinct labeling "
        "timepoints (per condition) [default: 3]. Counts distinct timepoints, "
        "not raw PSM rows — repeats at one timepoint don't count toward depth."),
    ria: Optional[float] = typer.Option(
        None, "-r", "--ria",
        help="Precursor enrichment (asymptotic labeled fraction, e.g. 0.06 for "
        "6%% v/v). Default: per-experiment from the manifest's "
        "precursor_enrichment (the SDRF value); given here, it overrides that "
        "for every curve. The legacy non-manifest path falls back to 0.06."),
    out: Path = typer.Option(
        Path("."), "-o", "--out", help="Output directory [default: .]."),
    fs: Optional[str] = typer.Option(
        None, "-f", "--fs", metavar="N|auto",
        help="Limited-isotopomer scoring: fit the FS on the leading channels "
        "iso0..isoN (give a single index N, e.g. '3' = iso0-iso3) to dodge "
        "co-eluting contaminants in the higher channels — integrate wide, fit "
        "narrow. Or 'auto': per-peptide width keyed on the natural-abundance "
        "envelope (iso0-3 for typical peptides, wider for broad ones; "
        "RIA-invariant). Default: score the full envelope."),
    workers: int = typer.Option(
        1, "-W", "--workers", metavar="N",
        help="Worker *processes* for the per-peptide fit [default: 1]. The "
        "parallelism lever: dispatches over a process pool to sidestep the GIL "
        "(the per-peptide fit is GIL-bound). Results are identical regardless of N "
        "(per-peptide deterministic seed)."),
    exclude_mbr: bool = typer.Option(
        False, "--exclude-mbr",
        help="Drop match-between-runs data points (evidence='mbr') before "
        "fitting. MBR points are used by default; this is the with/without-MBR "
        "A/B lever. The n_mbr / n_clean output columns report the split either way."),
    min_spep: Optional[int] = typer.Option(
        None, "--min-spep", metavar="N",
        help="Curation floor on a peptidoform's labelling-site count (Spep): drop "
        "those below N before fitting, so under-powered curves (too few sites for "
        "the envelope to shift measurably) never reach the results or rollup. "
        "Default is label-aware (8 for hw/D2O, 6 for o18 — o18's +2 Da/site shifts "
        "more per site); 0 disables. Tune up for short time series / slow turnover."),
    fraction_collapse: str = typer.Option(
        "sum", "--fraction-collapse", metavar="sum|anchor",
        help="How to combine LC fractions / technical replicates of one "
        "(peptidoform, charge, biorep, timepoint) before fitting (manifest path; "
        "SDRF comment[fraction identifier] defines the fractions). 'sum' (default): "
        "sum each isoN channel across fractions and intensity-weight the mass/QC "
        "columns. 'anchor': keep only the highest-total-intensity fraction (legacy "
        "parity). Intensities are always combined before a single FS is solved — "
        "fraction FS values are never averaged."),
    fs_rail_drop: bool = typer.Option(
        True, "--fs-rail-drop/--no-fs-rail-drop",
        help="Drop per-timepoint FS points beyond a 0.05 margin of the physical [0,1] "
        "range (over-labelled ≥1.05 or below-natural ≤−0.05) before counting fit "
        "points / depth — such a point is a failed solve, not a real measurement. On by "
        "default for all fits (benchmark-picked rails; see the 2026-07-05 report); "
        "--no-fs-rail-drop reproduces the pre-1.2.0 numbers (rail-hits then fall to the "
        "R² gate)."),
) -> None:
    """Fit kinetic models to a D2O-labeling integrate time series."""
    import dataclasses

    import pandas as pd

    from riana.config import FitConfig
    from riana.core.fitting import (
        available_coefficient_presets, fit_run, load_aa_coefficients,
        load_o18_coefficients, peptide_summary,
    )
    from riana.io.writers import (
        ESTIMATE_FLOAT_FORMAT, make_provenance, write_dataframe_tsv,
    )
    from riana.logger import get_logger

    if workers > (os.cpu_count() or 1):
        raise typer.BadParameter(
            f"--workers {workers} exceeds CPU count ({os.cpu_count()}).")

    if (manifest is None) == (not riana_path):
        raise typer.BadParameter(
            "provide either positional _riana.txt file(s) (legacy path) or "
            "--manifest (the SDRF path), not both / neither.")

    # --fs: one of three forms — 'auto' (per-peptide init-width-keyed widening); a
    # single int N meaning "score iso0..isoN" (the easy form, N is the highest
    # channel index); or, for back-compat, the explicit leading-contiguous list
    # iso0..isoN. The solver scores LEADING channels, so a gap or a non-zero start
    # is a user error, not a silent reinterpretation. All resolve to a channel
    # COUNT = N+1 = len(list).
    score_channels: Optional[int] = None
    fs_auto = False
    if fs is not None and fs.strip().lower() == "auto":
        fs_auto = True
    elif fs is not None:
        # list(...) — _parse_number_list returns a tuple; compare as a list below.
        fs_list = list(_parse_number_list(fs, int, sort=True, unique=True))
        if len(fs_list) == 1:  # single int N -> iso0..isoN
            score_channels = fs_list[0] + 1
        elif fs_list == list(range(len(fs_list))):  # explicit leading list
            score_channels = len(fs_list)
        else:
            raise typer.BadParameter(
                f"--fs must be 'auto', a single channel index N (e.g. '3' = "
                f"iso0-iso3), or a leading-contiguous list from iso0 (e.g. "
                f"'0 1 2 3'); got {fs!r}.")
        if score_channels < 2:
            raise typer.BadParameter(
                f"--fs needs >=2 channels (iso0 + a labelled one), e.g. '1' = "
                f"iso0-iso1; got {fs!r}.")

    # --coefficients is required for both labels — hw uses a per-AA D₂O table,
    # o18 uses a length-model table; the two CSV formats differ, so pick the
    # preset matching the label.
    if not coefficients:
        presets = " | ".join(available_coefficient_presets())
        fmt = ("(feature, coefficient) — e.g. juber_2026_o18_ac16" if label == "o18"
               else "(amino_acid, coefficient) — e.g. deberneh_2025_rss | ilchenko_2019 "
               "| commerford_1983 | alamillo_2025_{ac16,ipsc,cm}")
        raise typer.BadParameter(
            f"--coefficients is required for --label {label}. Pass a bundled "
            f"preset ({presets}) or a path to a {fmt} CSV.")

    os.makedirs(out, exist_ok=True)
    try:
        config = FitConfig(
            model=model,
            label=label,
            k_p=float(kp),
            k_r=float(kr),
            r_p=float(rp),
            q_value=float(q_value),
            depth=int(depth),
            ria_max=(float(ria) if ria is not None
                     else FitConfig.__dataclass_fields__["ria_max"].default),
            score_channels=score_channels,
            fs_auto=fs_auto,
            workers=int(workers),
            out_dir=str(out),
            exclude_mbr=bool(exclude_mbr),
            fraction_collapse=fraction_collapse,
            min_spep=min_spep,
            fs_rail_drop=bool(fs_rail_drop),
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    logger = get_logger(__name__, str(out))
    logger.info(f"riana {__version__}")
    logger.info("fit (typed pipeline)")
    logger.info(
        "Spep gate: --min-spep %d%s",
        config.min_spep,
        " (label default)" if min_spep is None else "")
    if fs_auto:
        logger.info("limited-isotopomer scoring: --fs auto (per-peptide "
                    "init-width-keyed widening)")
    elif score_channels is not None:
        logger.info(
            f"limited-isotopomer scoring: --fs iso0-iso{score_channels - 1} "
            f"({score_channels} channels)"
        )

    _load_coeffs = load_o18_coefficients if label == "o18" else load_aa_coefficients
    coeffs = _load_coeffs(coefficients) if coefficients else {}
    if coeffs:
        logger.info(f"loaded {len(coeffs)} {label} coefficients from {coefficients}")

    from riana.progress import ProgressReporter
    progress = ProgressReporter(0, "fit", logger)

    if manifest is not None:
        from riana.core.pipeline import fit_project, resolve_manifest_write

        # Default / same-folder -o updates the project in place; a *different* -o
        # forks a self-contained derived project there, leaving this manifest
        # untouched — so pointing -o elsewhere never clobbers the input run.
        out_dir, target_manifest = resolve_manifest_write(manifest, out, "fit")
        if out_dir != Path(manifest).resolve().parent:
            logger.info(
                "--manifest + -o: forking a derived project into %s "
                "(input manifest left untouched)", out_dir)
        out = out_dir
        logger.info(f"fitting from manifest {manifest}")
        result_df = fit_project(config, manifest, coeffs, logger=logger,
                                ria_override=ria, progress_callback=progress)
        progress.close()
        id_source = str(manifest)
    else:
        dfs = [pd.read_table(p, comment="#") for p in riana_path]
        logger.info(f"read {len(dfs)} timepoint files; fitting ...")
        # The explicit-files path does NOT collapse fractions — it has no SDRF
        # identity, so rows sharing a (peptide, sample) are fit as independent
        # points (pseudo-replication). Warn if that is the case; the manifest path
        # is the one that collapses LC fractions / technical replicates.
        if fraction_collapse != "sum":
            logger.warning(
                "--fraction-collapse is ignored on the explicit-files path; it "
                "only applies to --manifest fits.")
        _combined = pd.concat(dfs, ignore_index=True)
        if {"concat", "sample"}.issubset(_combined.columns) and \
                _combined.duplicated(subset=["concat", "sample"]).any():
            logger.warning(
                "multiple rows per (peptide, sample) found across the input files "
                "— these are fit as independent points (pseudo-replication). For "
                "LC-fraction / technical-replicate data, fit via the SDRF "
                "--manifest path so fractions collapse per (peptide, biorep, "
                "timepoint) first.")
        result_df = fit_run(config, dfs, coeffs, progress_callback=progress)
        progress.close()
        id_source = ",".join(str(p) for p in riana_path)

    os.makedirs(out, exist_ok=True)
    out_path = Path(out) / "riana_fit_peptides.txt"
    provenance = make_provenance(
        dataclasses.asdict(config),
        id_source=id_source,
        extra={"model": model, "label": label, "coefficients": str(coefficients)},
    )
    # Write the scalar per-peptide summary (the per-timepoint detail lives in the
    # fractions file); estimates rounded to ~6 sig figs.
    write_dataframe_tsv(out_path, peptide_summary(result_df), provenance,
                        include_index=True, float_format=ESTIMATE_FLOAT_FORMAT)
    logger.info(f"wrote {out_path}")
    written = [out_path]

    # M5: per-timepoint fraction-new (long format, one row per peptide-timepoint
    # with prediction-interval bounds) — the substrate for the protein rollup.
    fractions = result_df.attrs.get("fractions_long")
    if fractions is not None and not fractions.empty:
        frac_path = Path(out) / "riana_fit_fractions.txt"
        write_dataframe_tsv(frac_path, fractions, provenance,
                            include_index=False,
                            float_format=ESTIMATE_FLOAT_FORMAT)
        logger.info(
            f"wrote {frac_path} ({len(fractions)} peptide-timepoints)")
        written.append(frac_path)

    # Record stage="fit" rows so `rollup --manifest` can find these outputs and
    # the manifest is the single project index (the integrate→fit→rollup chain).
    # On a fork this is the derived project's manifest, not the input one.
    if manifest is not None:
        from riana.core.pipeline import record_stage_rows

        record_stage_rows(target_manifest, "fit", written, result_df, provenance)
        logger.info(f"recorded {len(written)} fit rows in {target_manifest}")

    n_fitted = int(result_df["k_deg"].notna().sum())
    n_well = int((result_df["R_squared"] >= 0.9).sum())
    logger.info(f"{n_fitted} peptides converged; {n_well} R²≥0.9")
    logger.handlers.clear()


# --------------------------------------------------------------------------- #
# rollup
# --------------------------------------------------------------------------- #
@app.command()
def rollup(
    fit_dir: Optional[Path] = typer.Argument(
        None, exists=True, file_okay=False, dir_okay=True, readable=True,
        help="Directory holding riana_fit_peptides.txt + riana_fit_fractions.txt "
        "from `riana fit`. Omit when using --manifest.",
    ),
    manifest: Optional[Path] = typer.Option(
        None, "--manifest", exists=True, dir_okay=False, readable=True,
        help="riana_manifest.tsv (the SDRF/project path). Finds the fit outputs "
        "from its stage='fit' rows, writes riana_rollup_proteins.txt next to the "
        "manifest, and records a stage='protein' row. -o is ignored here."),
    model: str = typer.Option(
        "simple", "-m", "--model",
        help="Kinetic model for the protein refit — match how the peptides were "
        "fit (simple, guan, fornasiero). Or 'linear simple': the linearized "
        "φ=log(1−θ) cross-sample model — fits a protein's conditions jointly for a "
        "per-condition k and a Δk test (writes delta_k / delta_k_p / "
        "delta_k_p_adj). Mutually exclusive with the ODE models."),
    phi_limit: float = typer.Option(
        -4.0, "--phi-limit", metavar="PHI",
        help="['linear simple' only] Plateau-truncation threshold in φ-space: "
        "points with φ=log(1−θ) at/below this are dropped per curve (saturated "
        "tail = measurement noise, not slope). −4 ≈ θ 0.98, −3 ≈ θ 0.95."),
    linear_weights: str = typer.Option(
        "wls", "--linear-weights", metavar="SCHEME",
        help="['linear simple' only] Weighting of the φ-space fit: 'wls' (default) "
        "weights by the delta-method inverse variance (1−θ)², taken from the fitted "
        "value — φ=log(1−θ) makes FS-scale noise heteroscedastic, so an unweighted "
        "fit is anti-conservative (~28% false positives at α=0.05, and k biased ~14% "
        "low in the fast tail). 'ols' restores the old unweighted fit (audit only)."),
    reference_condition: str = typer.Option(
        None, "--reference-condition", metavar="COND",
        help="['linear simple' only] Baseline condition of the Δk contrast — "
        "delta_k = k(test) − k(reference). Default: alphabetically first."),
    test_condition: str = typer.Option(
        None, "--test-condition", metavar="COND",
        help="['linear simple' only] Comparison condition of the Δk contrast "
        "(requires --reference-condition). Naming a pair contrasts exactly those "
        "two even when >2 conditions are present — an interim for multi-group "
        "projects before full all-pairwise. NOTE: the joint fit still pools the "
        "residual variance over ALL conditions in the data, so scope the SDRF "
        "conditions to the ones you mean to compare."),
    method: str = typer.Option(
        "weighted", "--method",
        help="Protein estimator. 'weighted' (default): biorep-aware per-timepoint "
        "inverse-variance collapse, then one refit (pseudoreplication-safe). "
        "'pooled': fit all peptide×timepoint points directly (pseudoreplication; "
        "for comparison). The median of peptide k is always carried as "
        "peptide_median_k for reference."),
    parsimony: str = typer.Option(
        "unique", "--parsimony",
        help="Protein attribution at summarize time. 'unique' (default): only "
        "single-accession peptides contribute (a shared peptide's envelope "
        "blends both proteins' turnover, so it can't be attributed). 'isoform': "
        "also fold isoform-only-shared peptides into the canonical entry, unless "
        "an isoform in the group carries its own unique peptide."),
    kp: float = typer.Option(
        0.5, "--kp", help="Precursor rate constant (two-compartment models)."),
    kr: float = typer.Option(
        0.05, "--kr", help="Reutilization rate constant (Fornasiero)."),
    rp: float = typer.Option(
        10.0, "--rp", help="Bound/free precursor ratio (Fornasiero)."),
    min_peptides: int = typer.Option(
        2, "--min-peptides",
        help="Min attributed peptides for a protein to be reported [default: 2]."),
    min_points: int = typer.Option(
        3, "--min-points",
        help="Min collapsed (t, theta) points for the refit [default: 3]."),
    min_fit_points: Optional[int] = typer.Option(
        None, "--min-fit-points", metavar="N",
        help="Peptide-level biological-replicate gate: keep only peptidoforms fit "
        "on >= N points. Default auto = 2 for a SINGLE-TIMEPOINT experiment "
        "(detected from one distinct labeling time — R2 is degenerate there, so "
        "curation rides on this replicate floor + --k-cv, with R2/rescue bypassed), "
        "off otherwise. Set 3+ to require more replicates, or 1 to disable."),
    min_spep: Optional[int] = typer.Option(
        None, "--min-spep", metavar="N",
        help="Optional Spep (labelling-site) admission gate before rollup — drop "
        "peptides below N. Off by default (the primary gate is at `fit`, so a "
        "manifest rollup already inherits it); set this for explicit-file inputs "
        "or belt-and-braces."),
    min_r2: Optional[float] = typer.Option(
        0.8, "--min-r2", metavar="R2",
        help="Peptide R² admission gate before rollup (default 0.8). Keep a "
        "peptide if R² ≥ this, OR (flat-curve rescue) R² ≥ --rescue-r2 and "
        "k_cv < --k-cv. The gate is needed for good within-protein geom-CV (the "
        "inverse-variance weighting alone under-curates; see reports/). Pass 0 "
        "(or any value ≤ 0) to disable it entirely."),
    k_cv: float = typer.Option(
        0.2, "--k-cv",
        help="Flat-curve rescue: max relative uncertainty of k̂ — "
        "k_cv = (ci_hi−ci_lo)/(2·|k|), a scale-free CV of the rate constant — for "
        "a low-R² peptide to still be admitted (only with --min-r2). Being "
        "scale-free it needs no retuning across time-series ranges or k units "
        "(/day vs /h), unlike an absolute SE. Set ≤ 0 to disable (R²-only gate)."),
    rescue_r2: float = typer.Option(
        0.6, "--rescue-r2",
        help="Flat-curve rescue: min R² floor for the --k-cv rescue (only with "
        "--min-r2). Guards against degenerate k≈0 fits whose CI collapses to a "
        "spuriously tight k_cv; the exact value barely matters (any floor above "
        "the negative-R² band works), 0.6 is validated."),
    workers: int = typer.Option(
        1, "-W", "--workers", metavar="N",
        help="Worker *processes* for the per-protein refit [default: 1]. The "
        "parallelism lever for the GIL-bound rollup: dispatches proteins over a "
        "process pool. Results are identical regardless of N (per-protein "
        "deterministic seed)."),
    exclude_mbr: bool = typer.Option(
        False, "--exclude-mbr",
        help="Drop match-between-runs fraction points (evidence='mbr') before the "
        "protein refit. MBR points are used by default; the n_mbr / n_clean output "
        "columns report the composition either way."),
    out: Path = typer.Option(
        Path("."), "-o", "--out", help="Output directory [default: .]."),
) -> None:
    """Roll per-peptide fits up to protein turnover.

    Two estimates per (experiment, condition, protein): the median of the
    peptides' k_deg, and a biorep-aware per-timepoint inverse-variance weighted
    refit over the M5 fraction-new substrate. Reads the `riana fit` outputs —
    from *fit_dir*, or via --manifest (the project path) — and writes
    ``riana_rollup_proteins.txt`` (+ stage='rollup' manifest rows on the --manifest path).
    """
    import pandas as pd

    from riana.core.pipeline import (
        fit_outputs_from_manifest,
        record_stage_rows,
        resolve_manifest_write,
    )
    from riana.core.protein import build_rollup_fractions, rollup_proteins
    from riana.exceptions import DataError
    from riana.io.writers import (
        ESTIMATE_FLOAT_FORMAT, make_provenance, write_dataframe_tsv,
    )
    from riana.logger import get_logger

    if (manifest is None) == (fit_dir is None):
        raise typer.BadParameter(
            "provide either FIT_DIR (the fit output folder) or --manifest, "
            "not both / neither.")
    if workers > (os.cpu_count() or 1):
        raise typer.BadParameter(
            f"--workers {workers} exceeds CPU count ({os.cpu_count()}).")
    # --min-r2 0 (or any value ≤ 0) disables the gate; the default 0.8 is on.
    min_r2 = min_r2 if (min_r2 is not None and min_r2 > 0) else None

    forked_into = None
    if manifest is not None:
        # Locate the fit outputs from the manifest's stage='fit' rows. Default /
        # same-folder -o writes the rollup next to the manifest (in place); a
        # *different* -o forks a self-contained derived project there (upstream
        # integrate+fit rows carried over), leaving this manifest untouched.
        try:
            pep_path, frac_path = (Path(p) for p in
                                   fit_outputs_from_manifest(manifest))
        except DataError as exc:
            raise typer.BadParameter(str(exc)) from exc
        out_dir, target_manifest = resolve_manifest_write(manifest, out, "rollup")
        if out_dir != Path(manifest).resolve().parent:
            forked_into = out_dir
        out = out_dir
        id_source = str(manifest)
    else:
        pep_path = fit_dir / "riana_fit_peptides.txt"
        frac_path = fit_dir / "riana_fit_fractions.txt"
        for p in (pep_path, frac_path):
            if not p.exists():
                raise typer.BadParameter(
                    f"{p.name} not found in {fit_dir}. Run `riana fit` first.")
        id_source = str(fit_dir)

    if test_condition and not reference_condition:
        raise typer.BadParameter(
            "--test-condition requires --reference-condition (the baseline of the "
            "Δk contrast).", param_hint="--test-condition")
    if test_condition and reference_condition and test_condition == reference_condition:
        raise typer.BadParameter(
            "--test-condition must differ from --reference-condition — a condition "
            "compared to itself has no Δk.", param_hint="--test-condition")

    os.makedirs(out, exist_ok=True)
    logger = get_logger(__name__, str(out))
    logger.info(f"riana {__version__}")
    logger.info(f"rollup (method={method}, parsimony={parsimony}, model={model})")
    if forked_into is not None:
        logger.info(
            "--manifest + -o: forking a derived project into %s "
            "(input manifest left untouched)", forked_into)

    peptides = pd.read_table(pep_path, comment="#")
    fractions = pd.read_table(frac_path, comment="#")
    from riana.progress import ProgressReporter
    progress = ProgressReporter(0, "rollup", logger)
    try:
        result = rollup_proteins(
            peptides, fractions, model=model, method=method,
            kinetic_kwargs=dict(k_p=kp, k_r=kr, r_p=rp),
            parsimony=parsimony, min_peptides=int(min_peptides),
            min_points=int(min_points), min_spep=min_spep, min_r2=min_r2,
            k_cv_max=float(k_cv), rescue_r2=float(rescue_r2),
            min_fit_points=min_fit_points, workers=int(workers),
            phi_limit=float(phi_limit), linear_weights=str(linear_weights),
            reference_condition=reference_condition,
            test_condition=test_condition,
            exclude_mbr=bool(exclude_mbr), progress_callback=progress,
        )
        progress.close()
    except (DataError, NotImplementedError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    out_path = Path(out) / "riana_rollup_proteins.txt"
    provenance = make_provenance(
        {"model": model, "parsimony": parsimony, "kp": kp, "kr": kr, "rp": rp,
         "min_peptides": min_peptides, "min_points": min_points,
         "min_fit_points": min_fit_points,
         "min_r2": min_r2, "k_cv": k_cv, "rescue_r2": rescue_r2, "method": method,
         "phi_limit": phi_limit, "linear_weights": linear_weights,
         "reference_condition": reference_condition,
         "test_condition": test_condition},
        id_source=id_source,
        extra={"method": method, "parsimony": parsimony, "model": model},
    )
    write_dataframe_tsv(out_path, result, provenance, include_index=False,
                        float_format=ESTIMATE_FLOAT_FORMAT)
    logger.info(f"wrote {out_path}")
    written = [out_path]

    # The collapsed inverse-variance-weighted (t, θ) the refit used (the GUI
    # curve substrate) — long/tidy, for a readable record of the weighting.
    rollup_fractions = build_rollup_fractions(result)
    if not rollup_fractions.empty:
        frac_path = Path(out) / "riana_rollup_fractions.txt"
        write_dataframe_tsv(frac_path, rollup_fractions, provenance,
                            include_index=False,
                            float_format=ESTIMATE_FLOAT_FORMAT)
        logger.info(f"wrote {frac_path} ({len(rollup_fractions)} points)")
        written.append(frac_path)

    # Record the stage='rollup' rows so the manifest indexes the whole chain (the
    # derived project's manifest on a fork, not the input one).
    if manifest is not None:
        record_stage_rows(target_manifest, "rollup", written, result, provenance)
        logger.info(f"recorded {len(written)} rollup rows in {target_manifest}")

    n_fit = int(result["k_deg"].notna().sum())
    logger.info(
        f"{len(result)} proteins ({method}); {n_fit} with a fitted k_deg")
    logger.handlers.clear()


# --------------------------------------------------------------------------- #
# gui
# --------------------------------------------------------------------------- #
@app.command()
def gui() -> None:
    """Launch the PySide6 GUI (requires the optional ``[gui]`` extra).

    The GUI forms build the *same* frozen IntegrationConfig/FitConfig as the CLI,
    so both surfaces validate identically. PySide6/qasync/pyqtgraph are imported
    lazily here so they never load on the core CLI path.
    """
    # Import from riana.gui.app (not the lazy riana.gui wrapper): this is the
    # import that actually pulls in PySide6/qasync/pyqtgraph, so a missing
    # [gui] extra is caught here with a clean message instead of a traceback.
    try:
        from riana.gui.app import run_gui
    except ImportError as exc:
        raise typer.BadParameter(
            "The GUI needs the optional dependencies. Install them with:\n"
            "    pip install 'riana[gui]'\n"
            f"(import failed: {exc})"
        ) from exc
    raise typer.Exit(code=run_gui())


def main() -> None:
    """Console-script entry point (``riana = riana.cli:main``)."""
    app()


if __name__ == "__main__":
    main()
