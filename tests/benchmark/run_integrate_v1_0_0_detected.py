"""Run the M3 Week 3 detected-peak integrator across the calibration set.

Mirrors :mod:`run_integrate_v0_9_0` but calls :func:`core.integration.integrate_run`
directly with the new science fields on (peak_method='detected',
baseline_method='linear', smoothing=None, polyorder=2). Outputs land at
``tests/data/calibration_d2o_mixing/<line>/integrate_outputs/v1.0.0_detected/<sample>_riana.txt``
ready to be picked up by ``bench_peak_boundary.py --method``,
``bench_m0_ma_recovery.py --inputs``, etc.

Phase E will fold the same dispatch into ``riana integrate --engine new``.
This standalone exists so Phase C's regression benchmarks have a runnable
input without waiting on the CLI rewrite.
"""

from __future__ import annotations

import argparse
import csv
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

VERSION_LABEL = "v1.0.0_detected"
PINNED_CONFIG = dict(
    isotopomers=(0, 1, 2, 3, 4, 5),
    q_value=0.01,
    r_time=0.33,
    mass_tol_ppm=15,
    threads=4,
    forced_mods=(0.0,),
    peak_method="detected",
    baseline_method="linear",
)


def run_one_fraction(line: str, row: dict, output_dir: Path) -> None:
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
    config = IntegrationConfig(sample=sample, **PINNED_CONFIG)

    t0 = time.time()
    with IndexedMzML(mzml_src) as mzml:
        df = integrate_run(config, psms, mzml)
    print(
        f"[{line} {sample}] integrated {len(df)} PSMs in {time.time()-t0:.1f}s",
        flush=True,
    )

    df.to_csv(out_file, sep="\t")


def run_line(line: str) -> Path:
    gt_path = (
        REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing" / line
        / "ground_truth.csv"
    )
    output_dir = (
        REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing" / line
        / "integrate_outputs" / VERSION_LABEL
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    with gt_path.open() as f:
        rows = list(csv.DictReader(f))
    print(f"\n=== {line.upper()}: {len(rows)} fractions -> {output_dir} ===")
    for row in rows:
        run_one_fraction(line, row, output_dir)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--line", choices=["ac16", "ipsc", "cm", "all"], default="all")
    args = parser.parse_args()
    targets = ["ac16", "ipsc", "cm"] if args.line == "all" else [args.line]
    for line in targets:
        run_line(line)


if __name__ == "__main__":
    main()
