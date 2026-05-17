"""
Smoothing sweep: how does `riana integrate --smoothing` change the calibration
benchmark?

For each smoothing window S in {none,3,5,7,9}, this:
  1. Runs `riana integrate` with `-S S` (S=none omits the flag) into
     integrate_outputs/v0.9.0_smoothing_<S>/  (cached — skipped if already
     present; S=none reuses the existing v0.9.0 baseline outputs).
  2. Re-runs the bench_aa_coefficients pipeline (Spep fit + per-AA regression).
  3. Re-runs the bench_fs_recovery pipeline.
  4. Records n_peptides, train/test R^2, median Spep, median FS bias, and the
     20 per-AA coefficients.

Outputs smoothing_sweep_summary.csv (one row per S) and prints per-AA
coefficient deltas vs the no-smoothing baseline. This quantifies §2c point 3
of PROJECT_REVIEW.md — whether the polyorder-1 Savitzky-Golay smoothing helps
or distorts the integrated areas that feed the calibration.

WARNING: each non-cached smoothing value triggers a full 9-fraction integrate
run (~30-90 min per cell line). Use --skip-integrate to run only the bench
phase on already-computed integrate outputs, or --dry-run to preview.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _helpers import forward_model as fm  # noqa: E402
from bench_aa_coefficients import (  # noqa: E402
    fit_aa_coefficients,
    fit_spep_per_peptide,
    load_and_curate,
)
from bench_fs_recovery import run_recovery, summarize  # noqa: E402
from run_integrate_v0_9_0 import PINNED_ARGS, REPO_ROOT, run_line  # noqa: E402

N_ISO = 4
DIAG_MIN_SPEP = 6
RANDOM_STATE = 1337


def integrate_dir_for(line: str, smoothing: str) -> Path:
    base = REPO_ROOT / 'tests' / 'data' / 'calibration_d2o_mixing' / line / 'integrate_outputs'
    if smoothing == 'none':
        return base / 'v0.9.0'
    return base / f'v0.9.0_smoothing_{smoothing}'


def ensure_integrate(line: str, smoothing: str, skip_integrate: bool,
                     dry_run: bool) -> Path:
    out_dir = integrate_dir_for(line, smoothing)
    n_existing = len(list(out_dir.glob('*_riana.txt'))) if out_dir.exists() else 0
    if n_existing >= 9:
        print(f'[{line} S={smoothing}] integrate outputs present ({n_existing}) '
              f'-> {out_dir}')
        return out_dir
    if skip_integrate:
        raise FileNotFoundError(
            f'--skip-integrate set but {out_dir} has only {n_existing}/9 outputs'
        )
    extra_args = list(PINNED_ARGS)
    if smoothing != 'none':
        extra_args += ['-S', smoothing]
    run_line(line, dry_run=dry_run, output_dir=out_dir, extra_args=extra_args)
    return out_dir


def bench_one(line: str, smoothing: str, integrate_dir: Path,
              ground_truth: pd.DataFrame) -> dict:
    riana_df = load_and_curate(integrate_dir, ground_truth, r2_min=0.95)
    spep_df = fit_spep_per_peptide(riana_df, n_iso=N_ISO)
    result = fit_aa_coefficients(spep_df, random_state=RANDOM_STATE)
    coeff_dict = dict(zip(result['coeff_df']['amino_acid'],
                          result['coeff_df']['coefficient']))
    rec_df = run_recovery(riana_df, coeff_dict, n_iso=N_ISO)
    summary = summarize(rec_df, min_spep=DIAG_MIN_SPEP)

    row = {
        'smoothing': smoothing,
        'n_peptides': int(len(spep_df)),
        'train_r2': result['train_r2'],
        'test_r2': result['test_r2'],
        'med_spep': float(spep_df['spep_continuous'].median()),
        'med_fs_bias': summary['diag_filtered']['median_bias'],
    }
    for aa in fm.AA_LIST:
        row[f'coef_{aa}'] = float(coeff_dict[aa])
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--line', choices=['ac16', 'ipsc'], required=True)
    parser.add_argument('--smoothing-values', default='none,3,5,7,9',
                        help='comma-separated; "none" omits -S, others map to -S N')
    parser.add_argument('--skip-integrate', action='store_true',
                        help='Only run the bench phase on existing integrate outputs')
    parser.add_argument('--dry-run', action='store_true',
                        help='Preview integrate commands without running')
    args = parser.parse_args()

    smoothing_values = [s.strip() for s in args.smoothing_values.split(',')]
    gt_path = (REPO_ROOT / 'tests' / 'data' / 'calibration_d2o_mixing'
               / args.line / 'ground_truth.csv')
    ground_truth = pd.read_csv(gt_path)
    output_dir = (REPO_ROOT / 'tests' / 'data' / 'calibration_d2o_mixing'
                  / args.line / 'benchmark_results' / 'v0.9.0')
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for smoothing in smoothing_values:
        print(f'\n========== {args.line} smoothing={smoothing} ==========')
        integrate_dir = ensure_integrate(args.line, smoothing,
                                         skip_integrate=args.skip_integrate,
                                         dry_run=args.dry_run)
        if args.dry_run:
            continue
        row = bench_one(args.line, smoothing, integrate_dir, ground_truth)
        rows.append(row)
        print(f'[{args.line} S={smoothing}] n={row["n_peptides"]}  '
              f'train_r2={row["train_r2"]:.4f}  test_r2={row["test_r2"]:.4f}  '
              f'med_fs_bias={row["med_fs_bias"]:+.4f}')

    if args.dry_run:
        print('\n[dry-run] no bench phase executed')
        return

    sweep_df = pd.DataFrame(rows)
    out_csv = output_dir / 'smoothing_sweep_summary.csv'
    sweep_df.to_csv(out_csv, index=False)

    print('\n=== Smoothing sweep summary ===')
    show = ['smoothing', 'n_peptides', 'train_r2', 'test_r2', 'med_spep', 'med_fs_bias']
    print(sweep_df[show].to_string(index=False, float_format='{:.4f}'.format))

    # Per-AA coefficient delta vs the no-smoothing baseline.
    if 'none' in sweep_df['smoothing'].values:
        base = sweep_df[sweep_df['smoothing'] == 'none'].iloc[0]
        coef_cols = [f'coef_{aa}' for aa in fm.AA_LIST]
        print('\n=== max |coef delta| vs smoothing=none ===')
        for _, r in sweep_df.iterrows():
            if r['smoothing'] == 'none':
                continue
            max_d = max(abs(r[c] - base[c]) for c in coef_cols)
            print(f'  S={r["smoothing"]}: max |Δcoef| = {max_d:.4f}  '
                  f'Δtest_r2 = {r["test_r2"] - base["test_r2"]:+.4f}  '
                  f'Δmed_fs_bias = {r["med_fs_bias"] - base["med_fs_bias"]:+.4f}')

    with (output_dir / 'smoothing_sweep_summary.json').open('w') as f:
        json.dump({'line': args.line, 'rows': rows}, f, indent=2)
    print(f'\n[done] wrote {out_csv}')


if __name__ == '__main__':
    main()
