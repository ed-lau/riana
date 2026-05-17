"""
Run `riana integrate` against the calibration dataset, producing per-fraction
*_riana.txt outputs that feed bench_aa_coefficients.py and bench_fs_recovery.py.

Inputs read (per cell line):
  - data/calibration_<line>/mzml/<basename>.mzML.gz
  - data/calibration_<line>/snakemake_results/time<N>/percolator/percolator.target.psms.txt
  - tests/data/calibration_d2o_mixing/<line>/ground_truth.csv (proportion -> mzML basename)

Outputs (per cell line):
  - tests/data/calibration_d2o_mixing/<line>/integrate_outputs/v0.9.0/time<N>_riana.txt

Pinned parameters (match the snakemake-era config_template.yaml except for the
0.9.0 mass-tolerance semantic correction: `-m N` now means ±N ppm, was ±N/2):
  -i 0..5         isotopomers 0-5
  -q 0.01         FDR <= 1%
  -r 0.33         retention window ±0.33 min
  -m 15           mass tolerance ±15 ppm
  -t 4            4 threads

`-D` (mass-step between isotopomers) is intentionally omitted: the current CLI
takes a float (default 1.003354835 = 13C-12C step), and at low-to-mid D2O
enrichment the +1 envelope is dominated by 13C contributions on tryptic
peptides anyway. The snakemake-era config_template.yaml passed `-D D` (a
label) but that no longer parses in 0.9.0+.

Each fraction is run against a tempdir holding a single symlinked mzML, because
`riana integrate` expects a folder and would otherwise pick an arbitrary file
based on file_idx ordering.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
VERSION_LABEL = 'v0.9.0'

PINNED_ARGS = [
    '-i', '0', '1', '2', '3', '4', '5',
    '-q', '0.01',
    '-r', '0.33',
    '-m', '15',
    '-t', '4',
]


def run_one_fraction(line: str, row: dict, output_dir: Path,
                     extra_args: list[str], dry_run: bool) -> None:
    proportion = float(row['nominal_proportion'])
    mzml_basename = row['mzml_filename']
    sample = f'time{proportion:g}'

    mzml_src = REPO_ROOT / 'data' / f'calibration_{line}' / 'mzml' / mzml_basename
    psms = (REPO_ROOT / 'data' / f'calibration_{line}' / 'snakemake_results'
            / sample / 'percolator' / 'percolator.target.psms.txt')

    if not mzml_src.exists():
        raise FileNotFoundError(f'mzML not found: {mzml_src}')
    if not psms.exists():
        raise FileNotFoundError(f'percolator PSMs not found: {psms}')

    with tempfile.TemporaryDirectory(prefix=f'riana_integrate_{line}_{sample}_') as td:
        mzml_dir = Path(td) / 'mzml'
        mzml_dir.mkdir()
        os.symlink(mzml_src, mzml_dir / mzml_basename)

        cmd = [
            sys.executable, '-m', 'riana', 'integrate',
            str(mzml_dir),
            str(psms),
            '-s', sample,
            '-o', str(output_dir),
            *extra_args,
        ]
        print(f'[{line} {sample}] {" ".join(cmd)}', flush=True)
        if dry_run:
            return

        t0 = time.time()
        subprocess.run(cmd, check=True)
        print(f'[{line} {sample}] done in {time.time() - t0:.1f}s', flush=True)

        # riana integrate writes to <out>/<sample>_riana.txt; verify it landed.
        expected = output_dir / f'{sample}_riana.txt'
        if not expected.exists():
            # Some versions of riana also write into a sample subdir.
            alt = output_dir / sample / f'{sample}_riana.txt'
            if alt.exists():
                shutil.move(str(alt), str(expected))
            else:
                raise FileNotFoundError(f'expected output not produced: {expected}')


def run_line(line: str, dry_run: bool, output_dir: Path | None = None,
             extra_args: list[str] | None = None) -> Path:
    """Run integrate for all 9 fractions of a cell line. Returns the output dir."""
    gt_path = (REPO_ROOT / 'tests' / 'data' / 'calibration_d2o_mixing' / line
               / 'ground_truth.csv')
    if output_dir is None:
        output_dir = (REPO_ROOT / 'tests' / 'data' / 'calibration_d2o_mixing' / line
                      / 'integrate_outputs' / VERSION_LABEL)
    if extra_args is None:
        extra_args = PINNED_ARGS
    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)
    with gt_path.open() as f:
        rows = list(csv.DictReader(f))
    print(f'\n=== {line.upper()}: {len(rows)} fractions -> {output_dir} ===')
    for row in rows:
        run_one_fraction(line, row, output_dir, extra_args=extra_args,
                         dry_run=dry_run)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--line', choices=['ac16', 'ipsc', 'all'], default='all')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print commands without executing')
    args = parser.parse_args()

    targets = ['ac16', 'ipsc'] if args.line == 'all' else [args.line]
    for line in targets:
        run_line(line, dry_run=args.dry_run)


if __name__ == '__main__':
    main()
