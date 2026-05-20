"""
M3 Week 0: observed-vs-predicted m0 (and full-envelope) mass-abundance.

The headline regression-gating metric for M3. Given a frozen per-cell-line
coefficient table from `build_frozen_tables.py`, this script:

  1. Estimates per-peptide Spep = round(Σ coeff_aa · n_aa), as in
     bench_fs_recovery.estimate_spep_from_coeffs (Σ is over the residue
     composition of the modification-stripped sequence; coeff_aa comes from
     the frozen `coefficient` column).
  2. Builds the predicted envelope at each nominal mixing fraction
     f = proportion/100:
         pred = (1 - f) · init_norm + f · final_norm
     where init/final come from forward_model's IsoSpec cache (natural
     abundance vs. Spep-many H-sites at RIA_D2O enrichment).
  3. Normalizes the observed iso0..iso(N_ISO-1) row-wise and scores:
         m0_err   = m0_obs  - m0_pred                     (per row)
         env_rmse = sqrt(mean((obs_norm - pred)^2))       (per row, N_ISO bins)

Why this and not bench_fs_recovery's median FS bias: M2 finding 1 showed
median FS bias is insensitive to integration quality (smoothing sweep moved
it ±0.002). m0_rmse / env_rmse are per-row and so resolve the per-peptide
spread that peak detection should move.

Why both populations: M2 finding 2 — the R²>0.95 curation gate discards
exactly the co-eluting / low-SNR peptides where peak detection helps most.
Reporting on the uncurated population is the only honest way to score
integration improvements that target those peptides. To get both from one
pass, we call `load_and_curate(r2_min=-1)` (no-op filter) and stamp each
row with `is_curated = r2 > 0.95`.

The frozen table is a *constant*, so its absolute bias cancels when scoring
two integrations against it — what regression-gating M3 needs.

Outputs (written under --output-dir):
  - m0_ma_recovery.csv          (one row per (peptide, proportion))
  - m0_ma_recovery_summary.json (all/curated/uncurated × overall and
                                 per-proportion metrics)
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


def estimate_spep(sequence: str, coeff_dict: dict[str, float]) -> int:
    clean = re.sub(r'\[.*?\]', '', sequence)
    return max(1, round(sum(coeff_dict.get(aa, 0.0) for aa in clean)))


def predict_envelope(
    sequence: str,
    charge: int,
    pep_mass: float,
    spep: int,
    proportion: float,
    n_iso: int,
) -> np.ndarray:
    """Mixing model: (1 - f)·init_norm + f·final_norm, both row-normalized."""
    init_env = fm._get_init_env(sequence, charge, pep_mass, n=n_iso)
    final_env = fm._get_final_env(sequence, charge, pep_mass, spep, n=n_iso)
    i_sum = init_env.sum()
    f_sum = final_env.sum()
    if i_sum == 0 or f_sum == 0:
        return np.full(n_iso, np.nan)
    init_norm = init_env / i_sum
    final_norm = final_env / f_sum
    f = proportion / 100.0
    return (1.0 - f) * init_norm + f * final_norm


def run_recovery(
    riana_df: pd.DataFrame,
    coeff_dict: dict[str, float],
    n_iso: int,
) -> pd.DataFrame:
    fm.clear_envelope_cache()
    rows = []
    t0 = time.time()
    n = len(riana_df)
    for i, row in riana_df.iterrows():
        seq = row['sequence']
        charge = int(row['charge'])
        pep_mass = row['pep_mass']
        spep = estimate_spep(seq, coeff_dict)

        obs = np.array([row[c] for c in ISO_COLS][:n_iso], dtype=float)
        obs_sum = obs.sum()
        if obs_sum == 0:
            continue
        obs_norm = obs / obs_sum

        pred = predict_envelope(seq, charge, pep_mass, spep,
                                row['proportion'], n_iso)
        env_rmse = float(np.sqrt(np.mean((obs_norm - pred) ** 2)))
        m0_obs = float(obs_norm[0])
        m0_pred = float(pred[0])

        rows.append({
            'concat': row['concat'],
            'sequence': seq,
            'charge': charge,
            'proportion': row['proportion'],
            'r2': float(row['r2']),
            'is_curated': bool(row['r2'] > 0.95),
            'spep_est': spep,
            'm0_obs': m0_obs,
            'm0_pred': m0_pred,
            'm0_err': m0_obs - m0_pred,
            'env_rmse': env_rmse,
            'total_intensity': float(obs_sum),
        })
        if (i + 1) % 2000 == 0:
            print(f'  scored {i + 1}/{n} rows in {time.time() - t0:.1f}s',
                  flush=True)
    print(f'  scored {n}/{n} rows in {time.time() - t0:.1f}s', flush=True)
    return pd.DataFrame(rows)


def _agg(df: pd.DataFrame) -> dict:
    if not len(df):
        return {
            'n': 0, 'n_peptides': 0,
            'm0_rmse': float('nan'), 'm0_mae': float('nan'),
            'm0_bias_median': float('nan'),
            'env_rmse_median': float('nan'), 'env_rmse_p95': float('nan'),
            'per_proportion': [],
        }
    per_prop = df.groupby('proportion').agg(
        m0_rmse=('m0_err', lambda x: float(np.sqrt(np.mean(x ** 2)))),
        m0_iqr=('m0_err', lambda x: float(np.percentile(x, 75) - np.percentile(x, 25))),
        m0_bias=('m0_err', lambda x: float(np.median(x))),
        env_rmse_median=('env_rmse', lambda x: float(np.median(x))),
        n=('m0_err', 'size'),
    ).round(6).reset_index()
    return {
        'n': int(len(df)),
        'n_peptides': int(df['concat'].nunique()),
        'm0_rmse': float(np.sqrt(np.mean(df['m0_err'].values ** 2))),
        'm0_mae': float(np.mean(np.abs(df['m0_err'].values))),
        'm0_bias_median': float(np.median(df['m0_err'].values)),
        'env_rmse_median': float(np.median(df['env_rmse'].values)),
        'env_rmse_p95': float(np.percentile(df['env_rmse'].values, 95)),
        'per_proportion': per_prop.to_dict(orient='records'),
    }


def summarize(rec_df: pd.DataFrame) -> dict:
    return {
        'all': _agg(rec_df),
        'curated': _agg(rec_df[rec_df['is_curated']]),
        'uncurated': _agg(rec_df[~rec_df['is_curated']]),
    }


def compare_integrations(
    methods: dict[str, Path],
    ground_truth: pd.DataFrame,
    coeff_dict: dict[str, float],
    n_iso: int,
) -> tuple[pd.DataFrame, dict]:
    """Score several named integration output dirs and tabulate by method.

    `methods` maps a method label -> a directory of *_riana.txt files. Each is
    scored exactly as a single bench_m0_ma_recovery run (curated + uncurated).
    Returns a long-format comparison DataFrame (method x population) and the
    raw per-method summaries.

    This is the shared engine behind the bench_peak_boundary / bench_baseline
    stubs: in M3 Week 0 only one method exists (fixed-window v0.9.0), so the
    comparison has a single method; Week 3 registers detected-boundary /
    baseline-subtracted runs as additional methods and the same table then
    shows whether the algorithm change moved m0_rmse.
    """
    rows = []
    summaries = {}
    for name, inputs_dir in methods.items():
        print(f'[compare] method {name!r} <- {inputs_dir}')
        riana_df = load_and_curate(inputs_dir, ground_truth, r2_min=-1.0)
        rec_df = run_recovery(riana_df, coeff_dict, n_iso=n_iso)
        summ = summarize(rec_df)
        summaries[name] = summ
        for pop in ('all', 'curated', 'uncurated'):
            a = summ[pop]
            rows.append({
                'method': name,
                'population': pop,
                'n': a.get('n', 0),
                'n_peptides': a.get('n_peptides', 0),
                'm0_rmse': a.get('m0_rmse', float('nan')),
                'm0_mae': a.get('m0_mae', float('nan')),
                'm0_bias_median': a.get('m0_bias_median', float('nan')),
                'env_rmse_median': a.get('env_rmse_median', float('nan')),
                'env_rmse_p95': a.get('env_rmse_p95', float('nan')),
            })
    return pd.DataFrame(rows), summaries


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--inputs', type=Path, required=True,
                        help='Directory of *_riana.txt files')
    parser.add_argument('--ground-truth', type=Path, required=True)
    parser.add_argument('--coefficients', type=Path, required=True,
                        help='Frozen d2o_aa_coefficients_<line>.csv from '
                             'build_frozen_tables.py (the constant reference)')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--n-iso', type=int, default=4)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f'[load] inputs={args.inputs}  coefficients={args.coefficients}')
    ground_truth = pd.read_csv(args.ground_truth)
    coeff_df = pd.read_csv(args.coefficients)
    coeff_dict = dict(zip(coeff_df['amino_acid'], coeff_df['coefficient']))

    # r2_min=-1 keeps the full uncurated population (M2 finding 2). The
    # downstream split into curated/uncurated uses the per-peptide r2 column.
    riana_df = load_and_curate(args.inputs, ground_truth, r2_min=-1.0)
    n_curated = int((riana_df['r2'] > 0.95).sum())
    print(f'[curate] {riana_df["concat"].nunique()} peptides, {len(riana_df)} rows '
          f'(of which {n_curated} curated rows with r2>0.95)')

    print(f'[recover] predicting envelopes, N_ISO={args.n_iso}')
    rec_df = run_recovery(riana_df, coeff_dict, n_iso=args.n_iso)
    rec_df.to_csv(args.output_dir / 'm0_ma_recovery.csv', index=False)

    summary = summarize(rec_df)
    summary['args'] = {
        'inputs': str(args.inputs),
        'ground_truth': str(args.ground_truth),
        'coefficients': str(args.coefficients),
        'n_iso': args.n_iso,
    }
    with (args.output_dir / 'm0_ma_recovery_summary.json').open('w') as f:
        json.dump(summary, f, indent=2)

    a = summary['all']; c = summary['curated']; u = summary['uncurated']
    print(f'\n[done] all       n={a["n"]:5d}  '
          f'm0_rmse={a["m0_rmse"]:.4f}  env_rmse(med)={a["env_rmse_median"]:.4f}')
    print(f'       curated   n={c["n"]:5d}  '
          f'm0_rmse={c["m0_rmse"]:.4f}  env_rmse(med)={c["env_rmse_median"]:.4f}')
    print(f'       uncurated n={u["n"]:5d}  '
          f'm0_rmse={u["m0_rmse"]:.4f}  env_rmse(med)={u["env_rmse_median"]:.4f}')
    print(f'       outputs -> {args.output_dir}')


if __name__ == '__main__':
    main()
