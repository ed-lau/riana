"""Run the M3 Week 3 (v1.0) integrator across the calibration set.

Mirrors :mod:`run_integrate_v0_9_0` but calls :func:`core.integration.integrate_run`
directly through the typed-records pipeline. Outputs land at
``tests/data/calibration_d2o_mixing/<line>/integrate_outputs/v1.0.0/<sample>_riana.txt``
and feed ``bench_peak_boundary.py --method``, ``bench_m0_ma_recovery.py
--inputs``, etc.

Config matches the Week 3 *shipped* defaults: ``peak_rt='ms2'``
and ``baseline_method='none'`` (see PROJECT_REVIEW.md §3 and commits
``3d8c715`` for the rationale — opt-in detection until a cross-proportion-
stable picker lands). Versus ``v0.9.0/`` this dir adds the Phase D
mass-accuracy columns (``iso{N}_obs_mz``, ``iso{N}_ppm_error``) and a
per-fraction ``<sample>_riana.drift.json`` sidecar; the legacy area
columns remain numerically faithful to ``v0.9.0/`` within 1e-3 rel
(Phase A parity).

Phase E folds the same dispatch into ``riana integrate --engine new``;
this standalone runner stays for benchmark convenience.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from riana.config import IntegrationConfig  # noqa: E402
from riana.core.integration import integrate_run  # noqa: E402
from riana.io.mzml import IndexedMzML  # noqa: E402
from riana.io.percolator import read_percolator  # noqa: E402

VERSION_LABEL = "v1.0.0"
PINNED_CONFIG = dict(
    isotopomers=(0, 1, 2, 3, 4, 5),
    q_value=0.01,
    extraction_half_width=0.33,
    mass_tol_ppm=15,
       
    peak_rt="ms2",
    integration_half_width=0.33,
    baseline_method="none",
)


def run_one_fraction(line: str, row: dict, output_dir: Path,
                     config_overrides: dict) -> None:
    proportion = float(row["nominal_proportion"])
    mzml_basename = row["mzml_filename"]
    sample = f"time{proportion:g}"

    mzml_src = REPO_ROOT / "data" / f"calibration_{line}" / "mzml" / mzml_basename
    psms_path = (
        REPO_ROOT / "data" / f"calibration_{line}" / "snakemake_results"
        / sample / "percolator" / "percolator.target.psms.txt"
    )
    if not mzml_src.exists():
        raise FileNotFoundError(f"mzML not found: {mzml_src}")
    if not psms_path.exists():
        raise FileNotFoundError(f"percolator PSMs not found: {psms_path}")

    out_file = output_dir / f"{sample}_riana.txt"
    if out_file.exists():
        print(f"[{line} {sample}] already done, skipping ({out_file})")
        return

    print(f"[{line} {sample}] loading PSMs + mzML ...", flush=True)
    psms = read_percolator(psms_path, sample=sample)
    config = IntegrationConfig(sample=sample, **{**PINNED_CONFIG, **config_overrides})

    t0 = time.time()
    with IndexedMzML(mzml_src) as mzml:
        df = integrate_run(config, psms, mzml)
    print(
        f"[{line} {sample}] integrated {len(df)} PSMs in {time.time()-t0:.1f}s",
        flush=True,
    )

    # Phase F3: provenance header (riana version + git SHA + config hash +
    # id source) on _riana.txt. Bench readers use comment='#' to skip.
    from riana.io.writers import make_provenance, write_dataframe_tsv

    provenance = make_provenance(
        dataclasses.asdict(config),
        id_source=str(psms_path),
        extra={"line": line, "mzml": mzml_basename},
    )
    write_dataframe_tsv(out_file, df, provenance, include_index=True)

    # Phase D sidecar: per-fraction calibration drift summary as a small JSON.
    drift = df.attrs.get("drift_summary")
    if drift is not None:
        drift_path = out_file.with_suffix(".drift.json")
        with drift_path.open("w") as f:
            json.dump(dataclasses.asdict(drift), f, indent=2)


def run_line(line: str, out_label: str, config_overrides: dict) -> Path:
    gt_path = (
        REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing" / line
        / "ground_truth.csv"
    )
    output_dir = (
        REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing" / line
        / "integrate_outputs" / out_label
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    with gt_path.open() as f:
        rows = list(csv.DictReader(f))
    print(f"\n=== {line.upper()}: {len(rows)} fractions -> {output_dir} ===")
    for row in rows:
        run_one_fraction(line, row, output_dir, config_overrides)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--line", choices=["ac16", "ipsc", "cm", "all"], default="all")
    parser.add_argument("--peak-rt", choices=["ms2", "apex"],
                        default=PINNED_CONFIG["peak_rt"])
    parser.add_argument("--baseline-method",
                        choices=["none", "noise_floor", "snip", "asls"],
                        default=PINNED_CONFIG["baseline_method"])
    parser.add_argument("--integration-half-width", default="0.33",
                        help="RT half-width in min, or 'auto' to detect boundaries")
    parser.add_argument("--width-rel-height", type=float, default=0.05,
                        help="apex-height fraction when width=auto (0.05=5%%, 0.5=FWHM)")
    parser.add_argument("--prominence-k", type=float, default=3.0,
                        help="apex prominence floor multiplier (apex/auto)")
    parser.add_argument("--out-label", default=VERSION_LABEL,
                        help="integrate_outputs/<label>/ subdir to write into")
    args = parser.parse_args()
    targets = ["ac16", "ipsc", "cm"] if args.line == "all" else [args.line]
    ihw = ("auto" if args.integration_half_width == "auto"
           else float(args.integration_half_width))
    # Extraction half-width: for ms2/fixed it equals the window; for apex/auto
    # it must be wider so the detected apex can sit off the MS2 RT and still get
    # its full window — be generous (+0.33 min apex-offset allowance).
    if args.peak_rt == "ms2" and ihw != "auto":
        extraction_half_width = float(ihw)
    else:
        w = 0.33 if ihw == "auto" else float(ihw)
        extraction_half_width = w + 0.33
    overrides = dict(peak_rt=args.peak_rt,
                     integration_half_width=ihw,
                     width_rel_height=args.width_rel_height,
                     prominence_k=args.prominence_k,
                     baseline_method=args.baseline_method,
                     extraction_half_width=extraction_half_width)
    for line in targets:
        run_line(line, args.out_label, overrides)


if __name__ == "__main__":
    main()
