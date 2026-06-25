"""
Port of NB90c: learn the ¹⁸O Spep length-model coefficients via the IsoSpec
forward model, mirroring bench_aa_coefficients.py (the D2O / NB87a trainer).

Pipeline (matches data/notebook/90c_O18_LengthModel_IsoSpec_AC16.ipynb):

  1. Load all <inputs>/*_riana.txt, attach `nominal_proportion` from
     ground_truth.csv by riana filename (same convention as the ac16/ipsc/cm
     D2O lines). iso0..iso5 are read directly (current pipeline output); a
     legacy m0..m5 schema is renamed for robustness.
  2. Curate (NB90c): dedupe (concat, proportion); keep peptides observed at all
     proportions; m0_ma = iso0 / sum(iso0..iso5); filter R²(m0_ma~prop) >
     --r2-min (NB90c default 0.95); parse SEQ_charge; compute pep_mass.
  3. Per-peptide Spep: minimize_scalar on the o18 forward-model SSE loss over
     [0.5, n_O_atoms], using the first --n-iso isotopomers (NB90c default 6 —
     the +2 Da ¹⁸O label spreads signal into m2/m4/m6, so iso0:5 all inform).
  4. **Length-model regression** (NB90c, NOT D2O's 20-per-AA — that has a K/R
     non-identifiability): Spep = b·(L-1) + c_D·D + c_E·E + c_N·N + c_Q·Q,
     intercept fixed 0, bounded lsq_linear (bvls), 80/20 split.

Outputs (under --output-dir):
  - o18_length_coefficients.csv  (feature, coefficient, std_error + metadata)
  - spep_per_peptide.csv         (concat, sequence, charge, pep_mass, spep_*)
  - summary.json                 (n_peptides, train_r2, test_r2, coef, args)

NB: validation is data-blocked — point --inputs at the o18 calibration series
once it is re-searched through the mzTab path (the harness is ready now).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, lsq_linear
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split

sys.path.insert(0, str(Path(__file__).parent))
from _helpers import o18_forward_model as o18  # noqa: E402

ISO_COLS = ['iso0', 'iso1', 'iso2', 'iso3', 'iso4', 'iso5']
_LEGACY_RENAME = {'m0': 'iso0', 'm1': 'iso1', 'm2': 'iso2',
                  'm3': 'iso3', 'm4': 'iso4', 'm5': 'iso5'}


def load_and_curate(
    inputs_dir: Path,
    ground_truth: pd.DataFrame,
    r2_min: float,
    drop_proportions: tuple[float, ...] = (),
) -> pd.DataFrame:
    """Load *_riana.txt, join proportions, apply NB90c curation."""
    gt_by_filename = ground_truth.set_index('riana_filename')['nominal_proportion']
    riana_files = sorted(inputs_dir.glob('*_riana.txt'))
    if not riana_files:
        raise FileNotFoundError(f'No *_riana.txt found in {inputs_dir}')

    parts = []
    for f in riana_files:
        if f.name not in gt_by_filename.index:
            raise KeyError(
                f'{f.name} not in ground_truth.csv. Known: {list(gt_by_filename.index)}'
            )
        proportion = float(gt_by_filename[f.name])
        if proportion in drop_proportions:
            continue
        df = pd.read_csv(f, sep='\t', comment='#')  # '#' skips the provenance header
        df = df.rename(columns=_LEGACY_RENAME)       # m0..m5 -> iso0..iso5 if legacy
        df['proportion'] = proportion
        parts.append(df)
    if not parts:
        raise ValueError(f'No input files left after drop_proportions={drop_proportions}')
    riana_df = pd.concat(parts, ignore_index=True)

    missing = [c for c in ISO_COLS if c not in riana_df.columns]
    if missing:
        raise KeyError(f'_riana.txt missing isotopomer columns {missing} '
                       f'(need iso0..iso5 for the o18 fit; have {list(riana_df.columns)})')

    riana_df = riana_df.drop_duplicates(subset=['concat', 'proportion'])
    riana_df = riana_df.sort_values('proportion')

    expected = len(parts)
    vc = riana_df.groupby('concat')['proportion'].nunique()
    riana_df = riana_df[riana_df['concat'].isin(vc[vc == expected].index)].copy()
    riana_df = riana_df[['concat', 'proportion'] + ISO_COLS].reset_index(drop=True)
    if riana_df.empty:
        raise ValueError(
            f'No peptide observed at all {expected} proportions — coverage curation '
            'emptied the set. Check the ground_truth.csv proportions and that the '
            'concat keys are stable across runs.')

    iso_sum = riana_df[ISO_COLS].sum(axis=1).replace(0, np.nan)
    riana_df['m0_ma'] = riana_df['iso0'] / iso_sum

    # R² per peptide for m0_ma vs proportion (explicit loop — pandas-version-stable)
    r2_map: dict[str, float] = {}
    for concat, group in riana_df.groupby('concat'):
        x = group['proportion'].values
        y = group['m0_ma'].values
        if len(x) < 3 or np.var(x) == 0 or np.var(y) == 0:
            r2_map[concat] = 0.0
        else:
            r2_map[concat] = float(np.corrcoef(x, y)[0, 1] ** 2)
    riana_df['r2'] = riana_df['concat'].map(r2_map)
    riana_df = riana_df[riana_df['r2'] > r2_min].copy()
    if riana_df.empty:
        raise ValueError(f'No peptide passed the R² > {r2_min} curation gate.')

    riana_df['sequence'] = riana_df['concat'].str.rsplit('_', n=1).str[0]
    riana_df['charge'] = riana_df['concat'].str.rsplit('_', n=1).str[1].astype(int)
    riana_df['pep_mass'] = riana_df.apply(
        lambda r: o18.calculate_ion_mz(r['sequence'], charge=r['charge']), axis=1)
    riana_df = riana_df[riana_df[ISO_COLS].sum(axis=1) > 0].copy()
    return riana_df.reset_index(drop=True)


def fit_spep_per_peptide(riana_df: pd.DataFrame, n_iso: int) -> pd.DataFrame:
    """Per-peptide continuous Spep via the o18 forward-model SSE loss."""
    o18.clear_envelope_cache()
    uniques = (riana_df[['concat', 'sequence', 'charge', 'pep_mass']]
               .drop_duplicates('concat').reset_index(drop=True))
    rows = []
    n = len(uniques)
    t0 = time.time()
    for i, row in uniques.iterrows():
        concat, seq = row['concat'], row['sequence']
        charge, pep_mass = int(row['charge']), row['pep_mass']
        pep_rows = riana_df[riana_df['concat'] == concat].sort_values('proportion')
        obs_raw = pep_rows[ISO_COLS].values[:, :n_iso].astype(float)
        row_sums = obs_raw.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        obs_norm = obs_raw / row_sums
        proportions_frac = pep_rows['proportion'].values / 100.0
        max_spep = float(o18.count_atoms(seq, charge=0)[2])   # total O atoms

        def loss(s: float) -> float:
            return o18.peptide_spep_loss(s, seq, charge, pep_mass, obs_norm, proportions_frac)

        res = minimize_scalar(loss, bounds=(0.5, max(1.0, max_spep)),
                              method='bounded', options={'xatol': 1e-3})
        rows.append({
            'concat': concat, 'sequence': seq, 'charge': charge, 'pep_mass': pep_mass,
            'spep_continuous': float(res.x),
            'spep_isospec': max(1, int(round(res.x))),
            'loss': float(res.fun),
        })
        if (i + 1) % 200 == 0:
            print(f'  fit {i + 1}/{n} ({(i + 1) / n * 100:.1f}%) in {time.time() - t0:.1f}s',
                  flush=True)
    print(f'  fit {n}/{n} peptides in {time.time() - t0:.1f}s', flush=True)
    return pd.DataFrame(rows)


def _build_design_matrix(spep_df: pd.DataFrame) -> pd.DataFrame:
    df = spep_df.copy()
    clean = df['sequence'].apply(o18.clean_seq)
    df['length_minus1'] = clean.str.len() - 1
    for aa in ['D', 'E', 'N', 'Q']:
        df[aa] = clean.str.count(aa)
    return df


def fit_length_model(spep_df: pd.DataFrame, random_state: int) -> dict:
    """5-param bounded length model (NB90c), intercept fixed 0."""
    df = _build_design_matrix(spep_df)
    X = df[o18.FEATURE_COLS].values.astype(float)
    y = df['spep_continuous'].values

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=random_state)

    res = lsq_linear(X_train, y_train,
                     bounds=(o18.FEATURE_BOUNDS_LOW, o18.FEATURE_BOUNDS_HIGH),
                     method='bvls')
    coef = res.x
    train_r2 = float(r2_score(y_train, X_train @ coef))
    test_r2 = float(r2_score(y_test, X_test @ coef))

    mse = float(np.mean((y_train - X_train @ coef) ** 2))
    try:
        xtx_inv = np.linalg.inv(X_train.T @ X_train)
        std_errors = np.sqrt(np.abs(mse * np.diag(xtx_inv)))
    except np.linalg.LinAlgError:
        std_errors = np.full(len(coef), np.nan)

    coeff_df = pd.DataFrame({
        'feature': o18.FEATURE_COLS,
        'coefficient': coef,
        'std_error': std_errors,
        'isotope': '18O — IsoSpec forward model (NB90c length model)',
        'intercept': 0.0,
        'train_r2': train_r2,
        'test_r2': test_r2,
    })
    return {
        'coeff_df': coeff_df, 'coef': coef.tolist(),
        'train_r2': train_r2, 'test_r2': test_r2,
        'n_train': int(len(X_train)), 'n_test': int(len(X_test)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputs', type=Path, required=True,
                        help='Directory of o18 calibration *_riana.txt files')
    parser.add_argument('--ground-truth', type=Path, required=True,
                        help='ground_truth.csv with riana_filename, nominal_proportion')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--n-iso', type=int, default=o18.N_ISO,
                        help=f'isotopomers used in the fit (NB90c default {o18.N_ISO})')
    parser.add_argument('--random-state', type=int, default=12345)
    parser.add_argument('--r2-min', type=float, default=0.95)
    parser.add_argument('--drop-proportion', type=float, nargs='+', default=[], metavar='PCT')
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f'[load] inputs={args.inputs}  ground_truth={args.ground_truth}')
    ground_truth = pd.read_csv(args.ground_truth)
    riana_df = load_and_curate(args.inputs, ground_truth, r2_min=args.r2_min,
                               drop_proportions=tuple(args.drop_proportion))
    print(f'[curate] {riana_df["concat"].nunique()} peptides, {len(riana_df)} rows '
          f'after R² > {args.r2_min}')

    print(f'[spep] per-peptide Spep via o18 forward model, N_ISO={args.n_iso}')
    spep_df = fit_spep_per_peptide(riana_df, n_iso=args.n_iso)
    spep_df.to_csv(args.output_dir / 'spep_per_peptide.csv', index=False)

    print('[regress] 5-param bounded length model (NB90c)')
    result = fit_length_model(spep_df, random_state=args.random_state)
    result['coeff_df'].to_csv(args.output_dir / 'o18_length_coefficients.csv', index=False)

    b, c_D, c_E, c_N, c_Q = result['coef']
    summary = {
        'n_peptides': int(len(spep_df)),
        'n_train': result['n_train'], 'n_test': result['n_test'],
        'train_r2': result['train_r2'], 'test_r2': result['test_r2'],
        'coef': {'b': b, 'c_D': c_D, 'c_E': c_E, 'c_N': c_N, 'c_Q': c_Q},
        'RIA_O18': o18.RIA_O18,
        'args': {'inputs': str(args.inputs), 'ground_truth': str(args.ground_truth),
                 'n_iso': args.n_iso, 'random_state': args.random_state, 'r2_min': args.r2_min},
    }
    with (args.output_dir / 'summary.json').open('w') as f:
        json.dump(summary, f, indent=2)

    print(f'\n[done] n_peptides={summary["n_peptides"]}  '
          f'train_r2={summary["train_r2"]:.4f}  test_r2={summary["test_r2"]:.4f}')
    print(f'       Spep = {b:.4f}·(L-1) + {c_D:.4f}·D + {c_E:.4f}·E '
          f'+ {c_N:.4f}·N + {c_Q:.4f}·Q')
    print(f'       outputs -> {args.output_dir}')


if __name__ == '__main__':
    main()
