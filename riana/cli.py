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
        "0 1 2 3 4 5", "-i", "--iso",
        help="Isotopomers to integrate, comma/space separated. Default "
        "'0 1 2 3 4 5' is the m0-m5 envelope the D2O fit consumes; pick a "
        "custom set for other workflows (e.g. SILAC cluster extraction via -F).",
    ),
    out: Path = typer.Option(
        Path("."), "-o", "--out", help="Output directory [default: .]."),
    q_value: float = typer.Option(
        1e-2, "-q", "--q_value", metavar="FDR",
        help="Integrate only PSMs with q-value below this [default: 1e-2]."),
    extraction_half_width: Optional[float] = typer.Option(
        None, "-r", "--extraction_half_width", "--r_time",
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
        "none", "--baseline",
        help="In-window baseline subtraction: none (default), noise_floor, "
        "snip, asls.",
    ),
    apex_selection: str = typer.Option(
        "tallest", "--apex-selection",
        help="Apex pick rule for apex/consensus: tallest (default) or nearest.",
    ),
    write_intensities: bool = typer.Option(
        False, "-w", "--write_intensities",
        help="Also write the pre-integration intensity trace."),
    mass_tol: Optional[int] = typer.Option(
        None, "-m", "--mass_tol", metavar="PPM",
        help="Mass tolerance half-width in ppm (±N ppm). If omitted, taken from "
        "the SDRF comment[precursor mass tolerance], else the 10 ppm default. "
        "Given here, overrides the SDRF."),
    smoothing: Optional[int] = typer.Option(
        None, "-S", "--smoothing",
        help="Savitzky-Golay smoothing window (odd int ≥ 3)."),
    mass_difference: float = typer.Option(
        1.003354835, "-D", "--mass_difference",
        help="Mass difference between isotopomers [default: 1.003354835]."),
    workers: int = typer.Option(
        1, "-W", "--workers", metavar="N",
        help="Runs to integrate concurrently on the --sdrf path (one mzML in "
        "memory per worker; 2-4 suits a many-timepoint time series) [default: "
        "1]. The parallelism lever — per-run extraction is GIL-bound serial."),
    no_rt_check: bool = typer.Option(
        False, "--no-rt-check",
        help="Disable the intake scan↔RT guard — the per-run check that the "
        "mzTab spectra_ref scans reconcile with this mzML's retention times "
        "(catches the quantms filename-prefix scan-scramble / wrong mzML↔mzTab "
        "pairing). Only disable for a run you know is correctly paired."),
    resume: bool = typer.Option(
        False, "--resume",
        help="On the --sdrf path, skip runs already integrated in the output's "
        "manifest at the current settings (each run's <stem>_riana.txt is written "
        "as it finishes, so a crashed run resumes where it stopped). Assumes the "
        "same SDRF/mzTab inputs."),
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

    isotopomers = _parse_number_list(iso, int, sort=True, unique=True)
    if not isotopomers:
        raise typer.BadParameter("--iso must list at least one isotopomer.")

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

    # The frozen config's __post_init__ owns the numeric domain validation
    # (mass_tol range, q_value range, peak_rt enum, ...); surface it as a clean
    # CLI error rather than a traceback.
    try:
        config = IntegrationConfig(
            sample=sample,
            isotopomers=isotopomers,
            mass_tol_ppm=mass_tol_ppm,
            extraction_half_width=ehw,
            peak_rt=peak_rt,
            integration_half_width=ihw,
            baseline_method=baseline_method,
            apex_selection=apex_selection,
            q_value=float(q_value),
            write_intensities=bool(write_intensities),
            smoothing=smoothing,
            mass_difference=float(mass_difference),
            out_dir=str(out),
            check_scan_rt=not no_rt_check,
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    logger = get_logger(__name__, str(out))
    logger.info(f"riana {__version__}")
    logger.info("integrate (typed pipeline)")
    logger.info(f"mass tolerance: ±{mass_tol_ppm} ppm (from {mass_tol_src})")

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
        help="Per-AA D2O labeling-site table — REQUIRED for --label hw. Either "
        "a bundled preset name (commerford | ac16 | ipsc | cm) or a path to a "
        "CSV with columns (amino_acid, coefficient).",
    ),
    model: str = typer.Option(
        "simple", "-m", "--model",
        help="Kinetic model: simple (default), guan, fornasiero."),
    label: str = typer.Option(
        "hw", "-l", "--label",
        help="Labeling chemistry: 'hw' (heavy water / D2O, default). 'o18' is "
        "recognized but its fit is being reimplemented post-M4. (Amino-acid / "
        "SILAC fitting was dropped — integrate SILAC peaks, fit L/(H+L) "
        "downstream.)",
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
        help="Fit only peptides seen in at least this many samples [default: 3]."),
    ria: float = typer.Option(
        0.06, "-r", "--ria",
        help="Precursor enrichment level (asymptotic D2O fraction, e.g. 0.06 "
        "for 6%% v/v) [default: 0.06]."),
    out: Path = typer.Option(
        Path("."), "-o", "--out", help="Output directory [default: .]."),
    fs: Optional[str] = typer.Option(
        None, "-f", "--fs",
        help="Reserved (post-M4): restrict the envelope SSE to a subset of "
        "isotopomer channels to reduce co-eluting-contaminant sensitivity. "
        "Currently ignored — the full integrated envelope is used."),
    workers: int = typer.Option(
        1, "-W", "--workers", metavar="N",
        help="Worker *processes* for the per-peptide fit [default: 1]. The "
        "parallelism lever: dispatches over a process pool to sidestep the GIL "
        "(the per-peptide fit is GIL-bound). Results are identical regardless of N "
        "(per-peptide deterministic seed)."),
) -> None:
    """Fit kinetic models to a D2O-labeling integrate time series."""
    import dataclasses

    import pandas as pd

    from riana.config import FitConfig
    from riana.core.fitting import (
        available_coefficient_presets, fit_run, load_aa_coefficients,
        peptide_summary,
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

    # --coefficients is required for the hw path; o18 errors in fit_run anyway.
    if label == "hw" and not coefficients:
        presets = " | ".join(available_coefficient_presets())
        raise typer.BadParameter(
            "--coefficients is required for --label hw. Pass a bundled preset "
            f"({presets}) or a path to a (amino_acid, coefficient) CSV.")

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
            ria_max=float(ria),
            fs_formula=fs,
            workers=int(workers),
            out_dir=str(out),
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    logger = get_logger(__name__, str(out))
    logger.info(f"riana {__version__}")
    logger.info("fit (typed pipeline)")
    if fs:
        logger.warning(
            "--fs is reserved for a post-M4 feature (channel-subset envelope "
            "SSE) and is currently ignored; the full envelope is used."
        )

    coeffs = load_aa_coefficients(coefficients) if coefficients else {}
    if coeffs:
        logger.info(f"loaded {len(coeffs)} AA coefficients from {coefficients}")

    if manifest is not None:
        from riana.core.pipeline import fit_project

        # The manifest's folder IS the project: write outputs next to it and
        # update that manifest, so a single --manifest drives integrate→fit→
        # rollup. -o is ignored on this path (warn if it was set elsewhere).
        proj_dir = Path(manifest).resolve().parent
        if Path(out).resolve() != proj_dir and str(out) != ".":
            logger.warning(
                "--manifest: outputs go next to the manifest (%s); ignoring -o %s",
                proj_dir, out)
        out = proj_dir
        logger.info(f"fitting from manifest {manifest}")
        result_df = fit_project(config, manifest, coeffs, logger=logger)
        id_source = str(manifest)
    else:
        dfs = [pd.read_table(p, comment="#") for p in riana_path]
        logger.info(f"read {len(dfs)} timepoint files; fitting ...")
        result_df = fit_run(config, dfs, coeffs)
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
    if manifest is not None:
        from riana.core.pipeline import record_stage_rows

        record_stage_rows(manifest, "fit", written, result_df, provenance)
        logger.info(f"recorded {len(written)} fit rows in {manifest}")

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
    reference_condition: str = typer.Option(
        None, "--reference-condition", metavar="COND",
        help="['linear simple' only] Baseline condition for the Δk contrast — "
        "delta_k = k(other) − k(reference). Default: alphabetically first."),
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
    min_r2: Optional[float] = typer.Option(
        None, "--min-r2", metavar="R2",
        help="Optional peptide R² admission gate before rollup (off by default — "
        "the inverse-variance weighting already down-weights noisy peptides). "
        "When set, keep a peptide if R² ≥ this, OR (slow-turnover admit) "
        "k ≤ --alt-k and SE ≤ --alt-se. Pass e.g. 0.8 to A/B against the "
        "unfiltered result."),
    alt_k: float = typer.Option(
        0.025, "--alt-k",
        help="Slow-turnover admit: max k_deg for a low-R² peptide to still be "
        "kept (only with --min-r2)."),
    alt_se: float = typer.Option(
        0.05, "--alt-se",
        help="Slow-turnover admit: max k_deg bootstrap SE (the 'sd' column) for "
        "a low-R² peptide to still be kept (only with --min-r2)."),
    workers: int = typer.Option(
        1, "-W", "--workers", metavar="N",
        help="Worker *processes* for the per-protein refit [default: 1]. The "
        "parallelism lever for the GIL-bound rollup: dispatches proteins over a "
        "process pool. Results are identical regardless of N (per-protein "
        "deterministic seed)."),
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

    from riana.core.pipeline import fit_outputs_from_manifest, record_stage_rows
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

    ignored_out = None
    if manifest is not None:
        # The manifest's folder is the project: locate the fit outputs from its
        # stage='fit' rows, write next to it, ignore -o.
        try:
            pep_path, frac_path = (Path(p) for p in
                                   fit_outputs_from_manifest(manifest))
        except DataError as exc:
            raise typer.BadParameter(str(exc)) from exc
        proj_dir = Path(manifest).resolve().parent
        if Path(out).resolve() != proj_dir and str(out) != ".":
            ignored_out = out
        out = proj_dir
        id_source = str(manifest)
    else:
        pep_path = fit_dir / "riana_fit_peptides.txt"
        frac_path = fit_dir / "riana_fit_fractions.txt"
        for p in (pep_path, frac_path):
            if not p.exists():
                raise typer.BadParameter(
                    f"{p.name} not found in {fit_dir}. Run `riana fit` first.")
        id_source = str(fit_dir)

    os.makedirs(out, exist_ok=True)
    logger = get_logger(__name__, str(out))
    logger.info(f"riana {__version__}")
    logger.info(f"rollup (method={method}, parsimony={parsimony}, model={model})")
    if ignored_out is not None:
        logger.warning(
            "--manifest: riana_rollup_proteins.txt goes next to the manifest (%s); "
            "ignoring -o %s", out, ignored_out)

    peptides = pd.read_table(pep_path, comment="#")
    fractions = pd.read_table(frac_path, comment="#")
    try:
        result = rollup_proteins(
            peptides, fractions, model=model, method=method,
            kinetic_kwargs=dict(k_p=kp, k_r=kr, r_p=rp),
            parsimony=parsimony, min_peptides=int(min_peptides),
            min_points=int(min_points), min_r2=min_r2,
            alt_k=float(alt_k), alt_se=float(alt_se), workers=int(workers),
            phi_limit=float(phi_limit), reference_condition=reference_condition,
        )
    except (DataError, NotImplementedError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    out_path = Path(out) / "riana_rollup_proteins.txt"
    provenance = make_provenance(
        {"model": model, "parsimony": parsimony, "kp": kp, "kr": kr, "rp": rp,
         "min_peptides": min_peptides, "min_points": min_points,
         "min_r2": min_r2, "alt_k": alt_k, "alt_se": alt_se, "method": method,
         "phi_limit": phi_limit, "reference_condition": reference_condition},
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

    # Record the stage='rollup' rows so the manifest indexes the whole chain.
    if manifest is not None:
        record_stage_rows(manifest, "rollup", written, result, provenance)
        logger.info(f"recorded {len(written)} rollup rows in {manifest}")

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
