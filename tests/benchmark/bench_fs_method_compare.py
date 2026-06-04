"""bench_fs_method_compare.py — legacy iso0-only vs new full-envelope FS solver.

Both methods receive the same per-peptide Spep (computed via
``algorithms.isotope_dist.spep_from_coefficients`` from a per-cell-line
``d2o_aa_coefficients_<line>.csv``). The only difference between
``legacy`` and ``new`` is the FS-extraction step:

- ``legacy``:  ``fs = (mi - a_0) / (a_max - a_0)``
  where ``mi = iso0/Σ iso``, ``a_0 = calculate_a0(seq, label=1)`` (natural
  abundance iso0 fraction), and ``a_max = a_0 * (1 - ria_max)**Spep``.
  Uses iso0 alone — the §2b "FS-denominator drift" path.

- ``new``: :func:`algorithms.isotope_dist.solve_fs_d2o` — full-envelope
  least-squares between ``obs_norm`` and
  ``(1-fs)·init_norm + fs·final_norm``.

Reports per-peptide / per-proportion FS for both methods, plus aggregate
recovery stats vs. the ground-truth ``nominal_proportion``. The recovery
delta isolates the **FS extraction method** itself, with Spep held
constant — answers "is the IsoSpec forward solver actually better than
iso0-only, given identical Spep input?"

Usage:
  python -m tests.benchmark.bench_fs_method_compare \\
    --inputs tests/data/calibration_d2o_mixing/ac16/integrate_outputs/v1.0.0 \\
    --ground-truth tests/data/calibration_d2o_mixing/ac16/ground_truth.csv \\
    --coefficients tests/data/calibration_d2o_mixing/ac16/d2o_aa_coefficients_ac16.csv \\
    --output-dir tests/data/calibration_d2o_mixing/ac16/benchmark_results/fs_method_compare
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from riana.algorithms.isotope_dist import solve_fs_d2o, spep_from_coefficients
from riana.algorithms.mass_calc import calculate_ion_mz
from riana.core.fitting import load_aa_coefficients
from riana.core.fsynthesis import calculate_a0, calculate_fs_m0


def _per_peptide_fs(
    df: pd.DataFrame,
    seq: str,
    spep: int,
    pep_mass: float,
    proportions: np.ndarray,
    ria_max: float,
    label: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute legacy and new FS for one peptide across all proportions.

    Returns (fs_legacy, fs_new), each shape (n_prop,) with NaN for
    timepoints where the iso0/Σ row is zero or the solver could not
    return a finite FS.
    """
    iso_cols = sorted(
        [c for c in df.columns if c.startswith("iso") and c[3:].isdigit()],
        key=lambda c: int(c[3:]),
    )
    obs = df[iso_cols].to_numpy(dtype=np.float64)
    sums = obs.sum(axis=1)

    # Legacy: iso0 / Σ iso → analytic (a - a_0) / (a_max - a_0).
    mi = np.where(sums > 0, obs[:, 0] / np.where(sums > 0, sums, 1), 0.0)
    fs_legacy = calculate_fs_m0(
        a=mi, seq=seq, label=label, ria_max=ria_max,
        num_labeling_sites=spep,
    )
    fs_legacy = np.asarray(fs_legacy, dtype=np.float64)
    # calculate_fs_m0 returns zeros when a_max == a_0; treat as NaN so the
    # recovery stats don't get pinned to 0 spuriously.
    fs_legacy = np.where(sums > 0, fs_legacy, np.nan)

    # New: full-envelope solve_fs_d2o per row.
    fs_new = np.array([
        solve_fs_d2o(seq, pep_mass, obs[i], spep, ria_max=ria_max, n_iso=len(iso_cols))
        if sums[i] > 0 else np.nan
        for i in range(len(obs))
    ])

    return fs_legacy, fs_new


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--inputs", type=Path, required=True,
                        help="Directory of *_riana.txt files")
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--coefficients", type=Path, required=True,
                        help="d2o_aa_coefficients_<line>.csv")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--q-value", type=float, default=0.01)
    parser.add_argument("--label", type=int, default=1)
    parser.add_argument("--ria", type=float, default=0.06,
                        help="precursor enrichment (default: 0.06 for 6%% v/v D2O)")
    parser.add_argument("--drop-proportion", type=float, nargs="+", default=[])
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    gt = pd.read_csv(args.ground_truth)
    if args.drop_proportion:
        gt = gt[~gt["nominal_proportion"].isin(args.drop_proportion)].copy()
        print(f"[load] dropped proportions: {args.drop_proportion}")
    coeffs = load_aa_coefficients(args.coefficients)
    print(f"[load] {len(coeffs)} AA coefficients from {args.coefficients}")

    # Concat all timepoints, attach the nominal proportion via riana_filename.
    frames = []
    for _, row in gt.iterrows():
        rfile = args.inputs / row["riana_filename"]
        if not rfile.exists():
            print(f"[skip] missing {rfile}")
            continue
        # M3 Week 4: _riana.txt now carries a provenance header; comment='#' skips it.
        f = pd.read_csv(rfile, sep="\t", index_col=0, comment="#")
        f = f[f["percolator q-value"] < args.q_value].copy()
        f["nominal_proportion"] = float(row["nominal_proportion"]) / 100.0
        frames.append(f)
    df = pd.concat(frames, ignore_index=True)
    print(f"[load] {len(df)} PSM-rows across {len(frames)} fractions "
          f"(q < {args.q_value})")

    # Per-peptide: compute both FS series across the fractions it was seen in.
    out_rows = []
    by_concat = df.groupby("concat", sort=False)
    n_done = 0
    for concat, peptide_rows in by_concat:
        seq = concat.rsplit("_", 1)[0]
        try:
            pep_mass = calculate_ion_mz(seq)
        except (KeyError, ValueError):
            continue

        spep_float = spep_from_coefficients(seq, coeffs)
        spep_int = max(1, int(round(spep_float)))
        proportions = peptide_rows["nominal_proportion"].to_numpy(dtype=np.float64)

        try:
            fs_legacy, fs_new = _per_peptide_fs(
                peptide_rows, seq, spep_int, pep_mass, proportions,
                ria_max=args.ria, label=args.label,
            )
        except (KeyError, ValueError):
            continue

        for i in range(len(peptide_rows)):
            out_rows.append({
                "concat": concat,
                "sample": peptide_rows["sample"].iloc[i],
                "nominal_proportion": proportions[i],
                "spep": spep_float,
                "fs_legacy": fs_legacy[i],
                "fs_new": fs_new[i],
            })
        n_done += 1
        if n_done % 500 == 0:
            print(f"  processed {n_done} peptides ...")

    out = pd.DataFrame(out_rows)
    out["err_legacy"] = out["fs_legacy"] - out["nominal_proportion"]
    out["err_new"] = out["fs_new"] - out["nominal_proportion"]
    out.to_csv(args.output_dir / "fs_method_compare.csv", index=False)
    print(f"[done] wrote {len(out)} rows to fs_method_compare.csv")

    def _agg(err: pd.Series) -> dict:
        finite = err.dropna()
        if finite.empty:
            return {"n": 0}
        return {
            "n": int(len(finite)),
            "median": float(finite.median()),
            "rmse": float(np.sqrt(np.mean(finite ** 2))),
            "iqr": float(finite.quantile(0.75) - finite.quantile(0.25)),
            "frac_within_0.05": float((finite.abs() <= 0.05).mean()),
            "frac_within_0.10": float((finite.abs() <= 0.10).mean()),
        }

    summary = {
        "n_peptides": int(out["concat"].nunique()),
        "n_rows": int(len(out)),
        "legacy_overall": _agg(out["err_legacy"]),
        "new_overall": _agg(out["err_new"]),
        "by_proportion": {
            f"{p:g}": {
                "legacy": _agg(out[out["nominal_proportion"] == p]["err_legacy"]),
                "new":    _agg(out[out["nominal_proportion"] == p]["err_new"]),
            }
            for p in sorted(out["nominal_proportion"].unique())
        },
        "args": {
            "inputs": str(args.inputs),
            "ground_truth": str(args.ground_truth),
            "coefficients": str(args.coefficients),
            "q_value": args.q_value,
            "label": args.label,
            "ria_max": args.ria,
            "drop_proportion": list(args.drop_proportion),
        },
    }
    with (args.output_dir / "fs_method_compare_summary.json").open("w") as f:
        json.dump(summary, f, indent=2)

    legacy = summary["legacy_overall"]
    new = summary["new_overall"]
    print()
    print(f"FS recovery vs. nominal proportion (Spep from {args.coefficients.name}):")
    print(f"  metric           legacy            new")
    print(f"  n                {legacy['n']:>10}     {new['n']:>10}")
    print(f"  median err       {legacy['median']:+10.4f}     {new['median']:+10.4f}")
    print(f"  RMSE             {legacy['rmse']:>10.4f}     {new['rmse']:>10.4f}")
    print(f"  IQR              {legacy['iqr']:>10.4f}     {new['iqr']:>10.4f}")
    print(f"  within ±0.05     {legacy['frac_within_0.05']*100:>9.1f}%     {new['frac_within_0.05']*100:>9.1f}%")
    print(f"  within ±0.10     {legacy['frac_within_0.10']*100:>9.1f}%     {new['frac_within_0.10']*100:>9.1f}%")


if __name__ == "__main__":
    main()
