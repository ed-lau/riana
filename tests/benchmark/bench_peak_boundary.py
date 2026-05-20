"""
M3 Week 0 stub: peak-boundary integration comparison.

Compares observed-vs-predicted m0 / full-envelope RMSE (the
bench_m0_ma_recovery scoring) across integration methods that differ in how
the chromatographic peak boundary is determined.

WEEK 0 STATUS — STUB. The only integration method that exists today is
fixed-window integration (v0.9.0: every MS1 scan within ±r_time of the PSM
span is summed, see PROJECT_REVIEW.md §2c point 1). So this script currently
records a one-method baseline. It is committed now, before the rewrite, so
that the moment Week 3 lands `algorithms/peaks.py` there is an unchanged
gate to compare against.

WEEK 3 — register the new integration runs as extra `--method` entries:
  detected_boundary   scipy.signal.find_peaks + peak_widths
  skyline             Skyline-style boundary + co-elution grouping
The comparison table then shows whether peak detection moves m0_rmse. Per
M2 finding 2, watch the *uncurated* population row most closely — the
R²>0.95 gate discards exactly the co-eluting / low-SNR peptides peak
detection is meant to rescue.

The frozen per-cell-line coefficient table (build_frozen_tables.py) is the
constant reference: scoring every method against the same predicted m0/mA
preserves the relative ranking, so the freeze bias cancels.

Usage:
  python bench_peak_boundary.py \
    --method fixed_window=<integrate_dir> [--method detected=<dir> ...] \
    --ground-truth <ground_truth.csv> \
    --coefficients <d2o_aa_coefficients_<line>.csv> \
    --output-dir <dir>

Outputs (written under --output-dir):
  - peak_boundary_comparison.csv   (method x population: m0/env RMSE)
  - peak_boundary_comparison.json  (per-method per-proportion detail)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from bench_m0_ma_recovery import compare_integrations  # noqa: E402


def parse_method(spec: str) -> tuple[str, Path]:
    """Parse a `name=path` --method argument."""
    if '=' not in spec:
        raise argparse.ArgumentTypeError(
            f'--method must be name=path, got {spec!r}')
    name, _, path = spec.partition('=')
    if not name or not path:
        raise argparse.ArgumentTypeError(
            f'--method must be name=path, got {spec!r}')
    return name, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--method', type=parse_method, action='append',
                        required=True, metavar='NAME=DIR',
                        help='named integration output dir; repeatable')
    parser.add_argument('--ground-truth', type=Path, required=True)
    parser.add_argument('--coefficients', type=Path, required=True,
                        help='frozen d2o_aa_coefficients_<line>.csv')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--n-iso', type=int, default=4)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    methods = dict(args.method)
    if len(methods) != len(args.method):
        parser.error('duplicate --method name')
    if len(methods) == 1:
        print('[stub] one method only — recording the fixed-window baseline. '
              'Week 3 adds detected-boundary / skyline methods.')

    print(f'[load] ground_truth={args.ground_truth}  '
          f'coefficients={args.coefficients}')
    ground_truth = pd.read_csv(args.ground_truth)
    coeff_df = pd.read_csv(args.coefficients)
    coeff_dict = dict(zip(coeff_df['amino_acid'], coeff_df['coefficient']))

    comparison, summaries = compare_integrations(
        methods, ground_truth, coeff_dict, n_iso=args.n_iso)
    comparison.to_csv(args.output_dir / 'peak_boundary_comparison.csv',
                      index=False)

    out = {
        'benchmark': 'peak_boundary',
        'status': 'week0_stub_fixed_window_only' if len(methods) == 1
        else 'multi_method',
        'methods': {name: str(path) for name, path in methods.items()},
        'args': {
            'ground_truth': str(args.ground_truth),
            'coefficients': str(args.coefficients),
            'n_iso': args.n_iso,
        },
        'summaries': summaries,
    }
    with (args.output_dir / 'peak_boundary_comparison.json').open('w') as f:
        json.dump(out, f, indent=2)

    print('\n[done] method x population m0_rmse:')
    pivot = comparison.pivot(index='method', columns='population',
                             values='m0_rmse')
    print(pivot.round(4).to_string())
    print(f'       outputs -> {args.output_dir}')


if __name__ == '__main__':
    main()
