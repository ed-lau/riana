"""
M3 Week 0: build a frozen per-cell-line D2O amino-acid coefficient table.

M2 finding 3 (see PROJECT_REVIEW.md §3): per-cell-line frozen coefficient
tables are viable as a *constant* reference. Scoring two integrations against
the same frozen predicted m0/mA preserves their relative ranking, which is
what regression-gating M3 needs — so the table's absolute bias cancels.

`bench_aa_coefficients.py` emits a coefficient table from a single 80/20
split (random_state-dependent). That is fine for an R²/test-R² readout but
too noisy to freeze. This script instead bootstraps the per-AA non-negative
regression: it resamples peptides with replacement N_BOOT times, refits
`LinearRegression(fit_intercept=False, positive=True)` each time, and reports
the per-AA mean coefficient plus a bootstrap std and 95% percentile CI. The
out-of-bag R² across resamples is an honest generalization estimate (no
held-out split needed).

Input is the `spep_per_peptide.csv` produced by `bench_aa_coefficients.py`
(the expensive per-peptide Spep fit is already done and committed for
v0.9.0). In M3, re-run `bench_aa_coefficients.py` on the improved integration,
then re-run this script on the new spep table to re-derive the frozen table.

Outputs (written under --output-dir):
  - d2o_aa_coefficients_<line>.csv          frozen table; `coefficient`
      column is a drop-in for bench_fs_recovery.py / bench_m0_ma_recovery.py
      `--coefficients`.
  - d2o_aa_coefficients_<line>.provenance.json   source spep, integration,
      git SHA, n_peptides, n_boot, OOB R², date — so the freeze is traceable.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

sys.path.insert(0, str(Path(__file__).parent))
from _helpers import forward_model as fm  # noqa: E402


def count_aa_matrix(sequences: pd.Series) -> np.ndarray:
    """(n_pep, 20) residue-count matrix; modification brackets stripped."""
    rows = []
    for seq in sequences:
        clean = re.sub(r'\[.*?\]', '', seq)
        rows.append([clean.count(aa) for aa in fm.AA_LIST])
    return np.asarray(rows, dtype=float)


def bootstrap_coefficients(
    X: np.ndarray,
    y: np.ndarray,
    n_boot: int,
    random_state: int,
) -> dict:
    """Resample peptides with replacement, refit non-neg regression each time."""
    rng = np.random.default_rng(random_state)
    n = len(y)
    n_aa = X.shape[1]
    coefs = np.zeros((n_boot, n_aa))
    oob_r2 = np.full(n_boot, np.nan)
    all_idx = np.arange(n)

    t0 = time.time()
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        reg = LinearRegression(fit_intercept=False, positive=True)
        reg.fit(X[idx], y[idx])
        coefs[b] = reg.coef_
        oob = np.setdiff1d(all_idx, idx, assume_unique=False)
        if len(oob) > 1:
            oob_r2[b] = r2_score(y[oob], reg.predict(X[oob]))
        if (b + 1) % 500 == 0:
            print(f'  bootstrap {b + 1}/{n_boot} in {time.time() - t0:.1f}s',
                  flush=True)
    print(f'  bootstrap {n_boot}/{n_boot} in {time.time() - t0:.1f}s', flush=True)

    return {
        'mean': coefs.mean(axis=0),
        'std': coefs.std(axis=0, ddof=1),
        'ci_lo': np.percentile(coefs, 2.5, axis=0),
        'ci_hi': np.percentile(coefs, 97.5, axis=0),
        'frac_nonzero': (coefs > 0).mean(axis=0),
        'oob_r2_median': float(np.nanmedian(oob_r2)),
        'oob_r2_lo': float(np.nanpercentile(oob_r2, 2.5)),
        'oob_r2_hi': float(np.nanpercentile(oob_r2, 97.5)),
    }


def git_sha() -> str:
    try:
        out = subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'],
            capture_output=True, text=True, cwd=Path(__file__).parent,
        )
        return out.stdout.strip() or 'unknown'
    except Exception:  # noqa: BLE001 - provenance is best-effort
        return 'unknown'


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--spep', type=Path, required=True,
                        help='spep_per_peptide.csv from bench_aa_coefficients.py')
    parser.add_argument('--line', required=True,
                        help='cell-line tag, e.g. ac16 / ipsc')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--source', default='',
                        help='free-text description of the source integration '
                             '(e.g. "v0.9.0, -m 15, N_ISO=4")')
    parser.add_argument('--n-boot', type=int, default=2000)
    parser.add_argument('--random-state', type=int, default=1337)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f'[load] spep={args.spep}')
    spep_df = pd.read_csv(args.spep)
    X = count_aa_matrix(spep_df['sequence'])
    y = spep_df['spep_continuous'].to_numpy(dtype=float)
    print(f'[regress] bootstrapping {args.n_boot}x non-negative regression '
          f'over {len(y)} peptides')
    boot = bootstrap_coefficients(X, y, args.n_boot, args.random_state)

    coeff_df = pd.DataFrame({
        'amino_acid': fm.AA_LIST,
        'coefficient': boot['mean'],
        'boot_std': boot['std'],
        'ci_lo': boot['ci_lo'],
        'ci_hi': boot['ci_hi'],
        'boot_frac_nonzero': boot['frac_nonzero'],
    })
    out_csv = args.output_dir / f'd2o_aa_coefficients_{args.line}.csv'
    coeff_df.to_csv(out_csv, index=False)

    provenance = {
        'line': args.line,
        'source_spep': str(args.spep),
        'source_integration': args.source,
        'git_sha': git_sha(),
        'date': date.today().isoformat(),
        'n_peptides': int(len(y)),
        'n_boot': args.n_boot,
        'random_state': args.random_state,
        'oob_r2_median': boot['oob_r2_median'],
        'oob_r2_ci': [boot['oob_r2_lo'], boot['oob_r2_hi']],
    }
    out_json = args.output_dir / f'd2o_aa_coefficients_{args.line}.provenance.json'
    with out_json.open('w') as f:
        json.dump(provenance, f, indent=2)

    print(f'\n[done] frozen table -> {out_csv}')
    print(f'       n_peptides={provenance["n_peptides"]}  '
          f'OOB R²={boot["oob_r2_median"]:.4f} '
          f'[{boot["oob_r2_lo"]:.4f}, {boot["oob_r2_hi"]:.4f}]')


if __name__ == '__main__':
    main()
