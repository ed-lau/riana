"""
Port of NB87a steps 1-3: per-peptide Spep via IsoSpec forward model, then
per-AA non-negative linear regression.

Pipeline (matches data/notebook/87a_D2O_LearnAALabelingSites_IsoSpec_AC16.ipynb):

  1. Load all <riana_dir>/*_riana.txt, attach `nominal_proportion` from
     ground_truth.csv by riana filename.
  2. Curate:
     - dedupe (concat, proportion)
     - keep peptides observed at all 9 proportions
     - rename m0..m5 -> iso0..iso5
     - compute m0_ma = iso0 / sum(iso0..iso5), filter R^2(m0_ma~prop) > --r2-min
     - parse sequence + charge from concat (SEQ_charge)
     - compute pep_mass via forward_model.calculate_ion_mz
  3. Per-peptide Spep fit: minimize_scalar on peptide_spep_loss over
     [0.5, max_H_sites] for each peptide, using first N_ISO isotopomers.
  4. Per-AA non-neg regression: Spep_continuous ~ Sum(coeff_aa * n_aa),
     fit_intercept=False, positive=True, 80/20 split with --random-state.

Outputs (written under --output-dir):
  - d2o_amino_acid_coefficients.csv  (matches 87a_out schema)
  - spep_per_peptide.csv             (concat, sequence, charge, pep_mass, spep_*)
  - summary.json                     (n_peptides, train_r2, test_r2, args)
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
from scipy.optimize import minimize_scalar
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split

sys.path.insert(0, str(Path(__file__).parent))
from _helpers import forward_model as fm  # noqa: E402

ISO_COLS = ['iso0', 'iso1', 'iso2', 'iso3', 'iso4', 'iso5']


def load_and_curate(
    inputs_dir: Path,
    ground_truth: pd.DataFrame,
    r2_min: float,
    drop_proportions: tuple[float, ...] = (),
) -> pd.DataFrame:
    """Load *_riana.txt files, join proportions, apply NB87a curation.

    drop_proportions: nominal proportions to exclude *entirely* from curation
    (e.g. a known-weak acquisition such as cm's time50). The "observed at all
    proportions" coverage requirement then applies only to the surviving
    proportions — so dropping one weak fraction recovers peptides that would
    otherwise be discarded for missing in that single run. The dropped
    fraction's *_riana.txt is still scored downstream where the caller keeps
    the full population (bench_m0_ma_recovery runs with r2_min=-1).
    """
    gt_by_filename = ground_truth.set_index('riana_filename')['nominal_proportion']
    riana_files = sorted(inputs_dir.glob('*_riana.txt'))
    if not riana_files:
        raise FileNotFoundError(f'No *_riana.txt found in {inputs_dir}')

    parts = []
    for f in riana_files:
        if f.name not in gt_by_filename.index:
            raise KeyError(
                f'{f.name} not in ground_truth.csv. '
                f'Known: {list(gt_by_filename.index)}'
            )
        proportion = float(gt_by_filename[f.name])
        if proportion in drop_proportions:
            continue
        df = pd.read_csv(f, sep='\t')
        df['proportion'] = proportion
        df['file_name'] = f.name
        parts.append(df)
    if not parts:
        raise ValueError(
            f'No input files left after drop_proportions={drop_proportions}'
        )
    riana_df = pd.concat(parts, ignore_index=True)

    riana_df = riana_df.drop_duplicates(subset=['concat', 'proportion'])
    riana_df = riana_df.sort_values('proportion')

    expected = len(parts)
    vc = riana_df.groupby('concat')['proportion'].nunique()
    complete = vc[vc == expected].index
    riana_df = riana_df[riana_df['concat'].isin(complete)].copy()

    riana_df = riana_df.rename(columns={
        'm0': 'iso0', 'm1': 'iso1', 'm2': 'iso2',
        'm3': 'iso3', 'm4': 'iso4', 'm5': 'iso5',
    })
    riana_df = riana_df[['concat', 'proportion'] + ISO_COLS].reset_index(drop=True)

    iso_sum = riana_df[ISO_COLS].sum(axis=1).replace(0, np.nan)
    riana_df['m0_ma'] = riana_df['iso0'] / iso_sum

    def _r2(group: pd.DataFrame) -> pd.DataFrame:
        x = group['proportion'].values
        y = group['m0_ma'].values
        if len(x) < 3 or np.var(x) == 0:
            group['r2'] = 0.0
            return group
        r = np.corrcoef(x, y)[0, 1]
        group['r2'] = 0.0 if np.isnan(r) else r ** 2
        return group

    riana_df = riana_df.groupby('concat', group_keys=False).apply(_r2)
    riana_df = riana_df[riana_df['r2'] > r2_min].copy()

    riana_df['sequence'] = riana_df['concat'].str.rsplit('_', n=1).str[0]
    riana_df['charge'] = riana_df['concat'].str.rsplit('_', n=1).str[1].astype(int)
    riana_df['pep_mass'] = riana_df.apply(
        lambda r: fm.calculate_ion_mz(r['sequence'], charge=r['charge']),
        axis=1,
    )
    riana_df = riana_df[riana_df[ISO_COLS].sum(axis=1) > 0].copy()
    return riana_df.reset_index(drop=True)


def fit_spep_per_peptide(riana_df: pd.DataFrame, n_iso: int) -> pd.DataFrame:
    """For each unique peptide, fit continuous Spep across all proportions."""
    fm.clear_envelope_cache()
    uniques = (riana_df[['concat', 'sequence', 'charge', 'pep_mass']]
               .drop_duplicates('concat')
               .reset_index(drop=True))
    rows = []
    n = len(uniques)
    t0 = time.time()
    for i, row in uniques.iterrows():
        concat = row['concat']
        seq = row['sequence']
        charge = int(row['charge'])
        pep_mass = row['pep_mass']
        pep_rows = riana_df[riana_df['concat'] == concat].sort_values('proportion')
        obs_raw = pep_rows[ISO_COLS].values[:, :n_iso].astype(float)
        row_sums = obs_raw.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        obs_norm = obs_raw / row_sums
        proportions_frac = pep_rows['proportion'].values / 100.0
        max_spep = float(fm.count_atoms(seq, charge=0)[1])  # total H atoms

        def loss(s: float) -> float:
            return fm.peptide_spep_loss(
                s, seq, charge, pep_mass, obs_norm, proportions_frac
            )

        res = minimize_scalar(loss, bounds=(0.5, max_spep), method='bounded',
                              options={'xatol': 1e-3})
        rows.append({
            'concat': concat,
            'sequence': seq,
            'charge': charge,
            'pep_mass': pep_mass,
            'spep_continuous': float(res.x),
            'spep_isospec': max(1, int(round(res.x))),
            'loss': float(res.fun),
        })
        if (i + 1) % 200 == 0:
            print(f'  fit {i + 1}/{n} peptides ({(i + 1) / n * 100:.1f}%) '
                  f'in {time.time() - t0:.1f}s', flush=True)
    print(f'  fit {n}/{n} peptides in {time.time() - t0:.1f}s', flush=True)
    return pd.DataFrame(rows)


def fit_aa_coefficients(spep_df: pd.DataFrame, random_state: int) -> dict:
    """Per-AA non-negative LinearRegression, no intercept. 80/20 split."""
    def count_aas(seq: str) -> dict[str, int]:
        clean = re.sub(r'\[.*?\]', '', seq)
        return {aa: clean.count(aa) for aa in fm.AA_LIST}

    X = pd.DataFrame([count_aas(s) for s in spep_df['sequence']], columns=fm.AA_LIST)
    y = spep_df['spep_continuous'].values

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=random_state
    )
    reg = LinearRegression(fit_intercept=False, positive=True)
    reg.fit(X_train, y_train)
    train_r2 = float(reg.score(X_train, y_train))
    test_r2 = float(r2_score(y_test, reg.predict(X_test)))

    y_pred = reg.predict(X_train)
    mse = float(np.mean((y_train - y_pred) ** 2))
    try:
        xtx_inv = np.linalg.inv(X_train.values.T @ X_train.values)
        var_coef = mse * np.diag(xtx_inv)
        std_errors = np.sqrt(np.abs(var_coef))
    except np.linalg.LinAlgError:
        std_errors = np.full(len(fm.AA_LIST), np.nan)

    coeff_df = pd.DataFrame({
        'amino_acid': fm.AA_LIST,
        'coefficient': reg.coef_,
        'std_error': std_errors,
        'isotope': 'D2O (H2) — IsoSpec forward model',
        'train_r2': train_r2,
        'test_r2': test_r2,
    })
    return {
        'coeff_df': coeff_df,
        'train_r2': train_r2,
        'test_r2': test_r2,
        'n_train': int(len(X_train)),
        'n_test': int(len(X_test)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputs', type=Path, required=True,
                        help='Directory of *_riana.txt files')
    parser.add_argument('--ground-truth', type=Path, required=True,
                        help='ground_truth.csv with riana_filename, nominal_proportion')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--n-iso', type=int, default=4)
    parser.add_argument('--random-state', type=int, default=1337)
    parser.add_argument('--r2-min', type=float, default=0.95)
    parser.add_argument('--drop-proportion', type=float, nargs='+', default=[],
                        metavar='PCT',
                        help='nominal proportion(s) to exclude from curation, '
                             'e.g. --drop-proportion 50 to drop a weak '
                             'acquisition; coverage then requires the peptide '
                             'at every surviving proportion')
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f'[load] inputs={args.inputs}  ground_truth={args.ground_truth}')
    if args.drop_proportion:
        print(f'[load] dropping proportion(s) from curation: {args.drop_proportion}')
    ground_truth = pd.read_csv(args.ground_truth)
    riana_df = load_and_curate(args.inputs, ground_truth, r2_min=args.r2_min,
                               drop_proportions=tuple(args.drop_proportion))
    print(f'[curate] {riana_df["concat"].nunique()} peptides, {len(riana_df)} rows '
          f'after R² > {args.r2_min} filter')

    print(f'[spep] fitting per-peptide Spep, N_ISO={args.n_iso}')
    spep_df = fit_spep_per_peptide(riana_df, n_iso=args.n_iso)
    spep_df.to_csv(args.output_dir / 'spep_per_peptide.csv', index=False)

    print('[regress] per-AA non-negative LinearRegression')
    result = fit_aa_coefficients(spep_df, random_state=args.random_state)
    result['coeff_df'].to_csv(
        args.output_dir / 'd2o_amino_acid_coefficients.csv', index=False
    )

    summary = {
        'n_peptides': int(len(spep_df)),
        'n_train': result['n_train'],
        'n_test': result['n_test'],
        'train_r2': result['train_r2'],
        'test_r2': result['test_r2'],
        'args': {
            'inputs': str(args.inputs),
            'ground_truth': str(args.ground_truth),
            'n_iso': args.n_iso,
            'random_state': args.random_state,
            'r2_min': args.r2_min,
        },
    }
    with (args.output_dir / 'summary.json').open('w') as f:
        json.dump(summary, f, indent=2)

    print(f'\n[done] n_peptides={summary["n_peptides"]}  '
          f'train_r2={summary["train_r2"]:.4f}  test_r2={summary["test_r2"]:.4f}')
    print(f'       outputs -> {args.output_dir}')


if __name__ == '__main__':
    main()
