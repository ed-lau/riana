"""bench_niso_crossover.py — derive the N_ISO crossover for `--fs` widening.

Track B follow-on (the per-peptide `--fs` decision). The flat global `--fs 0 1 2 3`
(iso0-3) tightens recovery for the bulk but **over-truncates genuinely wide-envelope
peptides** (report 2026-06-23 §5: N_ISO≥12 pay +0.08 |fs−f|). This script re-derives,
at *per-integer* N_ISO granularity, where iso0-3 stops helping and starts hurting —
i.e. the crossover that keys an "extend, don't shift" widening
``score_channels = 4 if N_ISO ≤ T else 6`` (6 = all channels the fixed v1.0.0 integrate
captured, so this is a pure fit-side change, no re-integrate).

Per (peptide, proportion) row it computes, against the **same** per-line Spep used by
``bench_fs_method_compare`` (the new-method ``solve_fs_d2o``):
  - N_ISO  = ``len(adaptive_channel_masses(seq, pep_mass, ria_max))`` — the IsoSpec
             init∪final >1% envelope width (B1), the same quantity ``--iso auto`` uses.
  - err4   = solve(score_channels=4)  − f   (iso0-3, the shipped default)
  - err6   = solve(score_channels=None) − f  (all populated channels = iso0-5)

Outputs (under --output-dir):
  - niso_crossover.csv      per-row (concat, sample, nominal_proportion, n_iso, fs4, fs6)
  - niso_crossover.json     per-N_ISO |err| table + crossover + heuristic A/B vs flat

Usage:
  python tests/benchmark/bench_niso_crossover.py --line ac16
  python tests/benchmark/bench_niso_crossover.py --line all --threshold 11
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from riana.algorithms.isotope_dist import (  # noqa: E402
    adaptive_channel_masses,
    solve_fs_d2o,
    spep_from_coefficients,
)
from riana.algorithms.mass_calc import calculate_ion_mz  # noqa: E402
from riana.core.fitting import load_aa_coefficients  # noqa: E402

CALIB_ROOT = REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing"
# Same per-line recovery config as run_calibration_benchmark (cm drops its 50% outlier).
LINE_RECOVERY = {
    "ac16": {"coeff": "d2o_aa_coefficients_ac16.csv", "drop": []},
    "ipsc": {"coeff": "d2o_aa_coefficients_ipsc.csv", "drop": []},
    "cm": {"coeff": "d2o_aa_coefficients_cm_drop50.csv", "drop": [50.0]},
}


def _load_concat(line: str, q_value: float, drop: list[float]) -> pd.DataFrame:
    integrate_dir = CALIB_ROOT / line / "integrate_outputs" / "v1.0.0"
    gt = pd.read_csv(CALIB_ROOT / line / "ground_truth.csv")
    if drop:
        gt = gt[~gt["nominal_proportion"].isin(drop)].copy()
    frames = []
    for _, row in gt.iterrows():
        rfile = integrate_dir / row["riana_filename"]
        if not rfile.exists():
            continue
        f = pd.read_csv(rfile, sep="\t", index_col=0, comment="#")
        f = f[f["percolator q-value"] < q_value].copy()
        f["nominal_proportion"] = float(row["nominal_proportion"]) / 100.0
        frames.append(f)
    return pd.concat(frames, ignore_index=True)


def _per_peptide(df: pd.DataFrame, coeffs: dict, ria: float) -> pd.DataFrame:
    iso_cols = sorted(
        [c for c in df.columns if c.startswith("iso") and c[3:].isdigit()],
        key=lambda c: int(c[3:]),
    )
    n_cols = len(iso_cols)
    rows = []
    n_done = 0
    for concat, pep in df.groupby("concat", sort=False):
        seq = concat.rsplit("_", 1)[0]
        try:
            pep_mass = calculate_ion_mz(seq)
            n_iso = len(adaptive_channel_masses(seq, pep_mass, ria_max=ria))
        except (KeyError, ValueError):
            continue
        spep = max(1, int(round(spep_from_coefficients(seq, coeffs))))
        obs = pep[iso_cols].to_numpy(dtype=np.float64)
        sums = np.nansum(obs, axis=1)
        props = pep["nominal_proportion"].to_numpy(dtype=np.float64)
        samples = pep["sample"].to_numpy()
        for i in range(len(obs)):
            if sums[i] <= 0:
                continue
            fs4 = solve_fs_d2o(seq, pep_mass, obs[i], spep, ria_max=ria,
                               n_iso=n_cols, score_channels=4)
            fs6 = solve_fs_d2o(seq, pep_mass, obs[i], spep, ria_max=ria,
                               n_iso=n_cols, score_channels=None)
            rows.append({"concat": concat, "sample": samples[i],
                         "nominal_proportion": props[i], "n_iso": n_iso,
                         "fs4": fs4, "fs6": fs6})
        n_done += 1
        if n_done % 2000 == 0:
            print(f"  {n_done} peptides ...", flush=True)
    return pd.DataFrame(rows)


def _mae(s: pd.Series) -> float:
    s = s.dropna()
    return float(s.abs().mean()) if len(s) else float("nan")


def _within(s: pd.Series) -> float:
    s = s.dropna()
    return float((s.abs() <= 0.05).mean()) if len(s) else float("nan")


def analyse(out: pd.DataFrame, threshold: int) -> dict:
    out = out.copy()
    out["err4"] = out["fs4"] - out["nominal_proportion"]
    out["err6"] = out["fs6"] - out["nominal_proportion"]
    # Per-integer N_ISO stratification (tail collapsed at 12+ for sample size).
    out["niso_bin"] = np.where(out["n_iso"] >= 12, 12, out["n_iso"])
    strata = {}
    for b, g in out.groupby("niso_bin"):
        strata[int(b)] = {
            "n_rows": int(len(g)),
            "n_pep": int(g["concat"].nunique()),
            "mae_iso03": _mae(g["err4"]),
            "mae_all": _mae(g["err6"]),
            "delta_iso03_minus_all": _mae(g["err4"]) - _mae(g["err6"]),
        }
    # Crossover = smallest N_ISO at/above which iso0-3 MAE exceeds all-channel MAE
    # (and stays worse), i.e. where widening starts to pay.
    crossover = None
    bins = sorted(strata)
    for b in bins:
        if strata[b]["delta_iso03_minus_all"] > 0 and b >= 6:
            crossover = b
            break
    # Heuristic A/B: per-row pick iso0-3 below T, all-channels at/above T.
    out["err_heur"] = np.where(out["n_iso"] >= threshold, out["err6"], out["err4"])
    wide = out[out["n_iso"] >= 12]
    ab = {
        "flat_iso03": {"within_0.05": _within(out["err4"]), "mae": _mae(out["err4"]),
                       "wide_mae": _mae(wide["err4"])},
        "flat_all": {"within_0.05": _within(out["err6"]), "mae": _mae(out["err6"]),
                     "wide_mae": _mae(wide["err6"])},
        "heuristic": {"within_0.05": _within(out["err_heur"]),
                      "mae": _mae(out["err_heur"]), "wide_mae": _mae(wide["err_heur"]),
                      "threshold": threshold,
                      "n_pep_widened": int(out[out["n_iso"] >= threshold]
                                           ["concat"].nunique())},
    }
    return {"n_pep": int(out["concat"].nunique()), "n_rows": int(len(out)),
            "crossover_niso": crossover, "by_niso": strata, "ab": ab}


def run_line(line: str, threshold: int, out_root: Path) -> dict:
    rec = LINE_RECOVERY[line]
    ria = 0.0598
    print(f"\n=== {line.upper()} ===", flush=True)
    df = _load_concat(line, q_value=0.01, drop=rec["drop"])
    coeffs = load_aa_coefficients(CALIB_ROOT / line / rec["coeff"])
    out = _per_peptide(df, coeffs, ria)
    res = analyse(out, threshold)
    out_dir = out_root / line
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / "niso_crossover.csv", index=False)
    with (out_dir / "niso_crossover.json").open("w") as f:
        json.dump(res, f, indent=2)
    # Console table.
    print(f"  N_ISO  n_pep   MAE(iso0-3)  MAE(all)   Δ(iso03-all)")
    for b in sorted(res["by_niso"]):
        s = res["by_niso"][b]
        tag = "12+" if b == 12 else str(b)
        flag = "  <- iso0-3 worse" if s["delta_iso03_minus_all"] > 0 else ""
        print(f"  {tag:>5}  {s['n_pep']:>5}   {s['mae_iso03']:>10.4f}  "
              f"{s['mae_all']:>9.4f}  {s['delta_iso03_minus_all']:>+9.4f}{flag}")
    print(f"  crossover N_ISO = {res['crossover_niso']}")
    ab = res["ab"]
    print(f"  A/B within±0.05 | overall  wide-MAE")
    for k in ("flat_iso03", "flat_all", "heuristic"):
        print(f"    {k:11s} {ab[k]['within_0.05']*100:5.1f}%   {ab[k]['wide_mae']:.4f}")
    return res


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--line", choices=["ac16", "cm", "ipsc", "all"], default="ac16")
    parser.add_argument("--threshold", type=int, default=11,
                        help="N_ISO at/above which to widen iso0-3 -> all [default 11]")
    parser.add_argument("--output-dir", type=Path,
                        default=REPO_ROOT / "runs" / "niso_crossover")
    args = parser.parse_args()
    lines = ["ac16", "cm", "ipsc"] if args.line == "all" else [args.line]
    for line in lines:
        run_line(line, args.threshold, args.output_dir)


if __name__ == "__main__":
    main()
