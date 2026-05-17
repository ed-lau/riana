"""
Port of NB87a's N_ISO sweep.

D2O's +1 Da label modifies the m1 peak directly, so iso0:1 already carries
Spep information at N_ISO=2. This sweep re-runs the whole bench_aa_coefficients
+ bench_fs_recovery pipeline at N_ISO in {2,3,4,5,6} and reports how the per-AA
coefficients, train/test R^2, median Spep, and FS-recovery bias respond.

Low N_ISO is more robust to co-eluting interference (low-mass peaks have higher
SNR) but discards label signal in iso2+. The sweep quantifies that trade-off
and tells us whether RIANA needs to emit iso6+ at all.

Output (under --output-dir):
  - n_iso_sweep_summary.csv  (one row per N_ISO; train_r2, test_r2, med_spep,
                              med_fs_bias, coef_<AA> x 20 — matches 87a_out schema)

With --reference <csv>, diffs the produced summary against a reference
n_iso_sweep_summary.csv and prints the max coefficient delta (port-validation).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _helpers import forward_model as fm  # noqa: E402
from bench_aa_coefficients import (  # noqa: E402
    fit_aa_coefficients,
    fit_spep_per_peptide,
    load_and_curate,
)
from bench_fs_recovery import run_recovery, summarize  # noqa: E402


def run_sweep(riana_df: pd.DataFrame, n_iso_values: list[int],
              random_state: int, diag_min_spep: int) -> pd.DataFrame:
    rows = []
    for n_iso in n_iso_values:
        print(f'\n[N_ISO={n_iso}] fitting Spep ...', flush=True)
        spep_df = fit_spep_per_peptide(riana_df, n_iso=n_iso)
        result = fit_aa_coefficients(spep_df, random_state=random_state)
        coeff_df = result['coeff_df']
        coeff_dict = dict(zip(coeff_df['amino_acid'], coeff_df['coefficient']))

        rec_df = run_recovery(riana_df, coeff_dict, n_iso=n_iso)
        summary = summarize(rec_df, min_spep=diag_min_spep)
        med_fs_bias = summary['diag_filtered']['median_bias']

        row = {
            'n_iso': n_iso,
            'train_r2': result['train_r2'],
            'test_r2': result['test_r2'],
            'med_spep': float(spep_df['spep_continuous'].median()),
            'med_fs_bias': med_fs_bias,
        }
        for aa in fm.AA_LIST:
            row[f'coef_{aa}'] = float(coeff_dict[aa])
        rows.append(row)
        print(f'[N_ISO={n_iso}] train_r2={row["train_r2"]:.4f}  '
              f'test_r2={row["test_r2"]:.4f}  med_spep={row["med_spep"]:.2f}  '
              f'med_fs_bias={row["med_fs_bias"]:+.4f}', flush=True)
    return pd.DataFrame(rows)


def compare_to_reference(sweep_df: pd.DataFrame, ref_path: Path) -> None:
    ref = pd.read_csv(ref_path)
    coef_cols = [c for c in sweep_df.columns if c.startswith('coef_')]
    merged = sweep_df.merge(ref, on='n_iso', suffixes=('_new', '_ref'))
    max_delta = 0.0
    for c in coef_cols:
        d = (merged[f'{c}_new'] - merged[f'{c}_ref']).abs().max()
        max_delta = max(max_delta, float(d))
    r2_delta = max(
        (merged['train_r2_new'] - merged['train_r2_ref']).abs().max(),
        (merged['test_r2_new'] - merged['test_r2_ref']).abs().max(),
    )
    print(f'\n[reference] {ref_path}')
    print(f'  max |coef delta| across all N_ISO = {max_delta:.2e}')
    print(f'  max |R² delta|   across all N_ISO = {r2_delta:.2e}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--ground-truth', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--n-iso-values', default='2,3,4,5,6',
                        help='comma-separated N_ISO values')
    parser.add_argument('--random-state', type=int, default=1337)
    parser.add_argument('--r2-min', type=float, default=0.95)
    parser.add_argument('--diag-min-spep', type=int, default=6)
    parser.add_argument('--reference', type=Path, default=None,
                        help='reference n_iso_sweep_summary.csv to diff against')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    n_iso_values = [int(x) for x in args.n_iso_values.split(',')]

    ground_truth = pd.read_csv(args.ground_truth)
    riana_df = load_and_curate(args.inputs, ground_truth, r2_min=args.r2_min)
    print(f'[curate] {riana_df["concat"].nunique()} peptides, {len(riana_df)} rows')

    sweep_df = run_sweep(riana_df, n_iso_values,
                         random_state=args.random_state,
                         diag_min_spep=args.diag_min_spep)
    out_path = args.output_dir / 'n_iso_sweep_summary.csv'
    sweep_df.to_csv(out_path, index=False)

    print('\n=== N_ISO sweep summary ===')
    show = ['n_iso', 'train_r2', 'test_r2', 'med_spep', 'med_fs_bias']
    print(sweep_df[show].to_string(index=False, float_format='{:.4f}'.format))
    print(f'\n[done] wrote {out_path}')

    if args.reference is not None:
        compare_to_reference(sweep_df, args.reference)


if __name__ == '__main__':
    main()
