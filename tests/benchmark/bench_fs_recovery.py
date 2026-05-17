"""
Port of NB87a step 4: FS recovery.

Given a per-AA coefficient table produced by bench_aa_coefficients.py and the
same RIANA integrate outputs, solve fractional synthesis `fs` from each
(peptide, proportion) envelope and compare to the nominal mixing proportion.

Pipeline:
  1. Load + curate riana_df (identical to bench_aa_coefficients.py).
  2. Load per-AA coefficients CSV.
  3. For each row: estimate Spep = round(Sum(coeff_aa * n_aa)) from the
     sequence, then solve fs via forward_model.solve_fs_d2o using the first
     N_ISO isotopomers and bounds (-0.1, 1.2).
  4. Emit per-row CSV + summary statistics.

Outputs:
  - fs_recovery.csv      (one row per (peptide, proportion))
  - fs_recovery_summary.json  (median bias, per-proportion median, % outside [0,1])
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _helpers import forward_model as fm  # noqa: E402
from bench_aa_coefficients import ISO_COLS, load_and_curate  # noqa: E402

DIAG_PROP_EXCLUDE = [0.0, 100.0]
DIAG_MIN_SPEP_DEFAULT = 6


def estimate_spep_from_coeffs(sequence: str, coeff_dict: dict[str, float]) -> int:
    clean = re.sub(r'\[.*?\]', '', sequence)
    return max(1, round(sum(coeff_dict.get(aa, 0.0) for aa in clean)))


def run_recovery(riana_df: pd.DataFrame, coeff_dict: dict[str, float],
                 n_iso: int) -> pd.DataFrame:
    fm.clear_envelope_cache()
    rows = []
    t0 = time.time()
    n = len(riana_df)
    for i, row in riana_df.iterrows():
        seq = row['sequence']
        charge = int(row['charge'])
        pep_mass = row['pep_mass']
        spep = estimate_spep_from_coeffs(seq, coeff_dict)
        obs = [row[c] for c in ISO_COLS]
        fs = fm.solve_fs_d2o(seq, charge, pep_mass, obs, spep, n_iso=n_iso)
        rows.append({
            'concat': row['concat'],
            'sequence': seq,
            'charge': charge,
            'proportion': row['proportion'],
            'spep': spep,
            'fs_d2o': fs,
            'total_intensity': float(sum(obs)),
        })
        if (i + 1) % 2000 == 0:
            print(f'  recovered {i + 1}/{n} rows in {time.time() - t0:.1f}s',
                  flush=True)
    print(f'  recovered {n}/{n} rows in {time.time() - t0:.1f}s', flush=True)
    df = pd.DataFrame(rows).dropna(subset=['fs_d2o'])
    df['log2_intensity'] = np.log2(df['total_intensity'].clip(lower=1.0))
    df['fs_bias'] = df['fs_d2o'] - df['proportion'] / 100.0
    return df


def summarize(rec_df: pd.DataFrame, min_spep: int) -> dict:
    full_med_bias = float(rec_df['fs_bias'].median())
    full_per_prop = (rec_df.groupby('proportion')['fs_d2o']
                     .median()
                     .round(6)
                     .to_dict())
    n_total = len(rec_df)
    n_below = int((rec_df['fs_d2o'] < 0.0).sum())
    n_above = int((rec_df['fs_d2o'] > 1.0).sum())

    diag = rec_df[
        (~rec_df['proportion'].isin(DIAG_PROP_EXCLUDE))
        & (rec_df['spep'] >= min_spep)
    ]
    diag_med_bias = float(diag['fs_bias'].median()) if len(diag) else float('nan')
    diag_per_prop = (diag.groupby('proportion')['fs_d2o']
                     .median()
                     .round(6)
                     .to_dict()) if len(diag) else {}

    return {
        'all_rows': {
            'n': n_total,
            'median_bias': full_med_bias,
            'median_fs_by_proportion': full_per_prop,
            'frac_fs_below_0': n_below / n_total if n_total else None,
            'frac_fs_above_1': n_above / n_total if n_total else None,
        },
        'diag_filtered': {
            'n': int(len(diag)),
            'exclude_proportions': DIAG_PROP_EXCLUDE,
            'min_spep': min_spep,
            'median_bias': diag_med_bias,
            'median_fs_by_proportion': diag_per_prop,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputs', type=Path, required=True,
                        help='Directory of *_riana.txt files')
    parser.add_argument('--ground-truth', type=Path, required=True)
    parser.add_argument('--coefficients', type=Path, required=True,
                        help='d2o_amino_acid_coefficients.csv from bench_aa_coefficients.py')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--n-iso', type=int, default=4)
    parser.add_argument('--r2-min', type=float, default=0.95)
    parser.add_argument('--diag-min-spep', type=int, default=DIAG_MIN_SPEP_DEFAULT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f'[load] inputs={args.inputs}  coefficients={args.coefficients}')
    ground_truth = pd.read_csv(args.ground_truth)
    coeff_df = pd.read_csv(args.coefficients)
    coeff_dict = dict(zip(coeff_df['amino_acid'], coeff_df['coefficient']))

    riana_df = load_and_curate(args.inputs, ground_truth, r2_min=args.r2_min)
    print(f'[curate] {riana_df["concat"].nunique()} peptides, {len(riana_df)} rows')

    print(f'[recover] solving fs, N_ISO={args.n_iso}, fs bounds {fm.FS_BOUNDS}')
    rec_df = run_recovery(riana_df, coeff_dict, n_iso=args.n_iso)
    rec_df.to_csv(args.output_dir / 'fs_recovery.csv', index=False)

    summary = summarize(rec_df, min_spep=args.diag_min_spep)
    summary['args'] = {
        'inputs': str(args.inputs),
        'ground_truth': str(args.ground_truth),
        'coefficients': str(args.coefficients),
        'n_iso': args.n_iso,
        'r2_min': args.r2_min,
        'diag_min_spep': args.diag_min_spep,
    }
    with (args.output_dir / 'fs_recovery_summary.json').open('w') as f:
        json.dump(summary, f, indent=2)

    print(f'\n[done] n={summary["all_rows"]["n"]}  '
          f'median_bias={summary["all_rows"]["median_bias"]:+.4f}  '
          f'(diag-filtered n={summary["diag_filtered"]["n"]}, '
          f'median_bias={summary["diag_filtered"]["median_bias"]:+.4f})')
    print(f'       outputs -> {args.output_dir}')


if __name__ == '__main__':
    main()
