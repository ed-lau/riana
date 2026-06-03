"""Run the M3 Week 3 (v1.0) integrator across the calibration set.

Mirrors :mod:`run_integrate_v0_9_0` but calls :func:`core.integration.integrate_run`
directly through the typed-records pipeline. Outputs land at
``tests/data/calibration_d2o_mixing/<line>/integrate_outputs/v1.0.0/<sample>_riana.txt``
and feed ``bench_peak_boundary.py --method``, ``bench_m0_ma_recovery.py
--inputs``, etc.

Config matches the Week 3 *shipped* defaults: ``peak_method='fixed_window'``
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
    r_time=0.33,
    mass_tol_ppm=15,
    threads=4,
    forced_mods=(0.0,),
    peak_method="fixed_window",
    baseline_method="none",
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
    # Phase D sidecar: per-fraction calibration drift summary as a small JSON
    # alongside the TSV. Kept separate from _riana.txt so existing bench
    # loaders (pd.read_csv without comment handling) keep working.
    drift = df.attrs.get("drift_summary")
    if drift is not None:
        drift_path = out_file.with_suffix(".drift.json")
        with drift_path.open("w") as f:
            json.dump(dataclasses.asdict(drift), f, indent=2)


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
