"""
M3 Week 0: fit-module recovery benchmark via the pseudo-time trick.

The regression gate for the M3 Week 4 `core/fitting.py` rewrite (and the §2b
scientific fixes that come with it). Pulled forward from post-M3 per the
roadmap reconciliation, so the fit rewrite ships verified.

THE PSEUDO-TIME TRICK
The calibration series is a D2O *mixing* series, not a time series — sample i
is a fixed fraction f_i of fully-labeled lysate. But the simple kinetic model
is fs(t) = 1 - exp(-k_deg·t) (one_exponent with a_0=0, a_max=1). Invert it:
choose a k_deg0, and relabel sample i with a pseudo-time

    t_i = -ln(1 - f_i) / k_deg0

Then a peptide whose recovered fractional synthesis equals its mixing
fraction (fs_i ≈ f_i) must, when fit against t_i, return exactly k_deg0 —
for *every* peptide, regardless of biology. So the fit module is tested with
one source of variability (integration + FS math) instead of the usual
biological/kinetic cocktail. Deviation of the fitted k_deg from k_deg0 is the
metric: `k_rel_err = (k_deg - k_deg0) / k_deg0`.

The 100% fraction maps to t = -ln(0)/k_deg0 = ∞ and is excluded; the other 8
fractions (0%..87.5%) form the pseudo-time series. depth=3 (riana fit default)
keeps peptides seen in ≥3 of them.

WHAT A NON-ZERO v0.9.0 BASELINE MEANS
If fs_i ≠ f_i — e.g. the §2b "fit uses an incorrect a_0/a_max baseline"
defect, or integration noise — the fit will not return k_deg0. A systematic
v0.9.0 k bias is therefore expected and is exactly the regression the Week 4
§2b fixes should remove. This script records that baseline; it does not
assume v0.9.0 is correct.

This runs the real `riana fit` CLI (subprocess) on pseudo-time-relabeled
copies of the integrate outputs — it tests the shipped fit path end to end.

Outputs (written under --output-dir):
  - fit_recovery.csv          (one row per fitted peptide)
  - fit_recovery_summary.json (all peptides + well-fit R²≥0.9 subset)
  - pseudotime_map.csv is (re)generated at --pseudotime-map if absent.
"""
from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

# RIA_D2O used to come from _helpers.forward_model (M2 oracle, frozen). The
# bench's --ria default now lives inline at 0.06 (≈ 6% v/v D₂O) so it
# matches FitConfig.ria_max and the new CLI default in main.py.
_DEFAULT_RIA_D2O = 0.06

WELL_FIT_R2 = 0.9


def build_pseudotime_map(k_deg0: float, proportions: list[float]) -> pd.DataFrame:
    """t = -ln(1 - f) / k_deg0; f=100% maps to ∞ and is excluded."""
    rows = []
    for prop in sorted(proportions):
        f = prop / 100.0
        if f >= 1.0:
            rows.append({'nominal_proportion': prop, 'k_deg0': k_deg0,
                         'pseudo_time': np.nan, 'included': False})
        else:
            t = -np.log(1.0 - f) / k_deg0
            rows.append({'nominal_proportion': prop, 'k_deg0': k_deg0,
                         'pseudo_time': round(float(t), 6) or 0.0,
                         'included': True})
    return pd.DataFrame(rows)


def load_or_build_map(map_path: Path, k_deg0: float,
                      proportions: list[float]) -> pd.DataFrame:
    if map_path.exists():
        print(f'[map] reading {map_path}')
        return pd.read_csv(map_path)
    print(f'[map] {map_path} absent — building with k_deg0={k_deg0}')
    df = build_pseudotime_map(k_deg0, proportions)
    map_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(map_path, index=False)
    return df


def relabel_inputs(inputs_dir: Path, gt: pd.DataFrame, pt_map: pd.DataFrame,
                   dest: Path) -> list[Path]:
    """Copy each included fraction's *_riana.txt, rewriting `sample` -> pseudo-time.

    riana fit reads the time of a sample from its `sample` column via
    re.sub('[^0-9.]', '', sample) — so writing `time<pseudo_time>` there is all
    that is needed to turn the mixing series into a pseudo-time series.
    """
    fname_by_prop = dict(zip(gt['nominal_proportion'].round(4),
                             gt['riana_filename']))
    out = []
    for _, row in pt_map[pt_map['included']].iterrows():
        prop = round(float(row['nominal_proportion']), 4)
        t = float(row['pseudo_time'])
        fname = fname_by_prop.get(prop)
        if fname is None:
            raise KeyError(f'proportion {prop} not in ground_truth.csv')
        src = inputs_dir / fname
        if not src.exists():
            raise FileNotFoundError(f'integrate output missing: {src}')
        df = pd.read_csv(src, sep='\t', index_col=0)
        df['sample'] = f'time{t:.6f}'
        dst = dest / fname
        df.to_csv(dst, sep='\t')
        out.append(dst)
    return out


def run_riana_fit(files: list[Path], out_dir: Path, ria: float, label: int,
                  depth: int, q_value: float, threads: int) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, '-m', 'riana', 'fit',
        *[str(f) for f in files],
        '-m', 'simple',
        '-l', str(label),
        '-r', str(ria),
        '-q', str(q_value),
        '-d', str(depth),
        '-t', str(threads),
        '-o', str(out_dir),
    ]
    print(f'[fit] {" ".join(cmd)}')
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + '\n' + proc.stderr + '\n')
        raise RuntimeError(f'riana fit failed (exit {proc.returncode})')
    result = out_dir / 'riana_fit_peptides.txt'
    if not result.exists():
        raise FileNotFoundError(f'riana fit produced no output: {result}')
    return result


def _n_points(t_field) -> int:
    try:
        val = ast.literal_eval(t_field) if isinstance(t_field, str) else t_field
        return len(val) if hasattr(val, '__len__') else 1
    except (ValueError, SyntaxError):
        return 0


def score(fit_df: pd.DataFrame, k_deg0: float) -> pd.DataFrame:
    df = fit_df.rename(columns={fit_df.columns[0]: 'concat'}).copy()
    df['k_deg0'] = k_deg0
    df['n_points'] = df['t'].map(_n_points)
    df['k_rel_err'] = (df['k_deg'] - k_deg0) / k_deg0
    df['k_abs_err'] = df['k_deg'] - k_deg0
    keep = ['concat', 'k_deg', 'k_deg0', 'k_rel_err', 'k_abs_err',
            'R_squared', 'sd', 'n_points', 'protein id']
    return df[[c for c in keep if c in df.columns]]


def _agg(df: pd.DataFrame) -> dict:
    fitted = df.dropna(subset=['k_deg'])
    if not len(fitted):
        return {'n': 0}
    rel = fitted['k_rel_err'].to_numpy()
    return {
        'n': int(len(fitted)),
        'median_k_deg': float(fitted['k_deg'].median()),
        'median_k_rel_err': float(np.median(rel)),
        'k_rel_err_rmse': float(np.sqrt(np.mean(rel ** 2))),
        'k_rel_err_iqr': float(np.percentile(rel, 75) - np.percentile(rel, 25)),
        'frac_within_10pct': float(np.mean(np.abs(rel) <= 0.10)),
        'frac_within_25pct': float(np.mean(np.abs(rel) <= 0.25)),
    }


def summarize(scored: pd.DataFrame, k_deg0: float) -> dict:
    well_fit = scored[scored['R_squared'] >= WELL_FIT_R2]
    return {
        'k_deg0': k_deg0,
        'n_total': int(len(scored)),
        'n_failed': int(scored['k_deg'].isna().sum()),
        'all': _agg(scored),
        'well_fit': {'r2_min': WELL_FIT_R2, **_agg(well_fit)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--inputs', type=Path, required=True,
                        help='Directory of *_riana.txt integrate outputs')
    parser.add_argument('--ground-truth', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--pseudotime-map', type=Path,
                        default=Path('tests/data/calibration_d2o_mixing/'
                                     'pseudotime_map.csv'),
                        help='read if present, else generated here')
    parser.add_argument('--k-deg0', type=float, default=0.5,
                        help='target rate constant; only used when the '
                             'pseudotime map must be built [default: 0.5]')
    parser.add_argument('--ria', type=float, default=_DEFAULT_RIA_D2O,
                        help='ria_max for riana fit [default: 6%% D2O RIA]')
    parser.add_argument('--label', type=int, default=1)
    parser.add_argument('--depth', type=int, default=3)
    parser.add_argument('--q-value', type=float, default=0.01)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    gt = pd.read_csv(args.ground_truth)
    pt_map = load_or_build_map(args.pseudotime_map, args.k_deg0,
                               gt['nominal_proportion'].tolist())
    k_deg0 = float(pt_map['k_deg0'].iloc[0])
    included = pt_map[pt_map['included']]
    print(f'[map] k_deg0={k_deg0}  {len(included)} pseudo-time points '
          f'(100% excluded as t=∞)')

    with tempfile.TemporaryDirectory(prefix='bench_fit_recovery_') as td:
        relabel_dir = Path(td) / 'relabeled'
        relabel_dir.mkdir()
        files = relabel_inputs(args.inputs, gt, pt_map, relabel_dir)
        print(f'[relabel] {len(files)} fractions -> pseudo-time')
        fit_out = Path(td) / 'fit_out'
        result_path = run_riana_fit(files, fit_out, ria=args.ria,
                                    label=args.label, depth=args.depth,
                                    q_value=args.q_value, threads=args.threads)
        fit_df = pd.read_csv(result_path, sep='\t')

    scored = score(fit_df, k_deg0)
    scored.to_csv(args.output_dir / 'fit_recovery.csv', index=False)

    summary = summarize(scored, k_deg0)
    summary['args'] = {
        'inputs': str(args.inputs),
        'ground_truth': str(args.ground_truth),
        'pseudotime_map': str(args.pseudotime_map),
        'ria': args.ria,
        'label': args.label,
        'depth': args.depth,
        'q_value': args.q_value,
    }
    with (args.output_dir / 'fit_recovery_summary.json').open('w') as f:
        json.dump(summary, f, indent=2)

    a = summary['all']; w = summary['well_fit']
    print(f'\n[done] k_deg0={k_deg0}  fitted={a.get("n", 0)}  '
          f'failed={summary["n_failed"]}')
    print(f'       all      median_k={a.get("median_k_deg", float("nan")):.4f}  '
          f'rel_err(med)={a.get("median_k_rel_err", float("nan")):+.4f}  '
          f'rel_err_rmse={a.get("k_rel_err_rmse", float("nan")):.4f}  '
          f'within±10%={a.get("frac_within_10pct", float("nan")):.2f}')
    print(f'       well-fit median_k={w.get("median_k_deg", float("nan")):.4f}  '
          f'rel_err(med)={w.get("median_k_rel_err", float("nan")):+.4f}  '
          f'rel_err_rmse={w.get("k_rel_err_rmse", float("nan")):.4f}  '
          f'within±10%={w.get("frac_within_10pct", float("nan")):.2f}  '
          f'(n={w.get("n", 0)})')
    print(f'       outputs -> {args.output_dir}')


if __name__ == '__main__':
    main()
