"""bench_niso_crossover.py — derive the crossover for `--fs` widening, N_ISO vs init-width.

Track B follow-on (the per-peptide `--fs` decision). The flat global `--fs 0 1 2 3`
(iso0-3) tightens recovery for the bulk but **over-truncates genuinely wide-envelope
peptides** (report 2026-06-23 §5: N_ISO≥12 pay +0.08 |fs−f|). This script re-derives,
at *per-integer* granularity, where iso0-3 stops helping and starts hurting — keyed on
two candidate criteria:

  - **N_ISO**  = ``len(adaptive_channel_masses(seq, pep_mass, ria_max))`` — the IsoSpec
                 init∪final >1% envelope width (the ``--iso auto`` quantity). Folds in
                 RIA, but **inflates with labelling** so it decouples from the natural
                 width at high enrichment.
  - **init_w** = last isotopomer index >1% in the **natural-abundance (θ=0) envelope**
                 alone — a purely compositional, **RIA- and θ-invariant** width. The
                 measured driver of the clean widen (report 2026-06-23: N_ISO≥11
                 peptides have iso4/iso5 natural-populated 100%/95% — that, not
                 labelling, is why widening is θ-robust at 6%).

Per (peptide, proportion) row it solves ``solve_fs_d2o`` at score_channels ∈ {4,5,6}
(iso0-3 / iso0-4 / iso0-5; 6 = all the fixed v1.0.0 integrate captured), against the
per-line Spep. Heuristic arms A/B'd:

  - flat iso0-3 (shipped) / flat all
  - **niso**          : 4 if N_ISO ≤ niso-threshold else all          (binary)
  - **initw_binary**  : 4 if init_w ≤ initw-threshold-1 else all       (binary, RIA-stable)
  - **initw_graded**  : clamp(init_w+1, 4, 6) — *score exactly the natural width* (RIA-stable)

Outputs (under --output-dir/<line>): niso_crossover.csv (per-row), niso_crossover.json
(per-criterion strata + crossover + A/B).

Usage:
  python tests/benchmark/bench_niso_crossover.py --line all
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
    _binned_envelope,
    adaptive_channel_masses,
    get_peptide_distribution,
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


def _init_last_idx(seq: str, pep_mass: float, n: int = 16) -> int:
    """Last isotopomer index >1% in the natural-abundance (θ=0) envelope — the
    RIA-/θ-invariant compositional width. Channels 0..init_w are natural-populated."""
    init = get_peptide_distribution(seq, label="D2O")
    _im, ip = _binned_envelope(init, pep_mass, n)
    tot = sum(ip) or 1.0
    last = 0
    for i in range(n):
        if ip[i] / tot >= 0.01:
            last = i
    return last


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
            init_w = _init_last_idx(seq, pep_mass)
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
            fs = {
                k: solve_fs_d2o(seq, pep_mass, obs[i], spep, ria_max=ria,
                                n_iso=n_cols, score_channels=(None if k == 6 else k))
                for k in (4, 5, 6)
            }
            rows.append({"concat": concat, "sample": samples[i],
                         "nominal_proportion": props[i], "n_iso": n_iso,
                         "init_w": init_w, "fs4": fs[4], "fs5": fs[5], "fs6": fs[6]})
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


def _strata(out: pd.DataFrame, key: str, cap: int) -> tuple[dict, int | None]:
    """Per-integer |err| for iso0-3 vs all, by `key` (binned at `cap`+). Crossover =
    smallest key value (≥ floor) where iso0-3 MAE first exceeds all-channel MAE."""
    binned = np.where(out[key] >= cap, cap, out[key])
    strata, floor = {}, 4 if key == "init_w" else 6
    for b, g in out.groupby(binned):
        strata[int(b)] = {
            "n_rows": int(len(g)), "n_pep": int(g["concat"].nunique()),
            "mae_iso03": _mae(g["err4"]), "mae_all": _mae(g["err6"]),
            "delta_iso03_minus_all": _mae(g["err4"]) - _mae(g["err6"]),
        }
    crossover = next((b for b in sorted(strata)
                      if strata[b]["delta_iso03_minus_all"] > 0 and b >= floor), None)
    return strata, crossover


def analyse(out: pd.DataFrame, niso_t: int, initw_t: int) -> dict:
    out = out.copy()
    out["err4"] = out["fs4"] - out["nominal_proportion"]
    out["err5"] = out["fs5"] - out["nominal_proportion"]
    out["err6"] = out["fs6"] - out["nominal_proportion"]

    niso_strata, niso_x = _strata(out, "n_iso", 12)
    initw_strata, initw_x = _strata(out, "init_w", 7)

    # Heuristic arms.
    out["err_niso"] = np.where(out["n_iso"] >= niso_t, out["err6"], out["err4"])
    out["err_initw_bin"] = np.where(out["init_w"] >= initw_t, out["err6"], out["err4"])
    graded = np.clip(out["init_w"] + 1, 4, 6)  # score exactly the natural width
    out["err_initw_grad"] = np.select(
        [graded == 4, graded == 5, graded == 6],
        [out["err4"], out["err5"], out["err6"]])

    wide = out[out["n_iso"] >= 12]  # the genuinely-wide tail, fixed reference
    arms = {
        "flat_iso03": "err4", "flat_all": "err6", "niso": "err_niso",
        "initw_binary": "err_initw_bin", "initw_graded": "err_initw_grad",
    }
    ab = {name: {"within_0.05": _within(out[col]), "mae": _mae(out[col]),
                 "wide_mae": _mae(wide[col])} for name, col in arms.items()}
    ab["niso"]["threshold"] = niso_t
    ab["initw_binary"]["threshold"] = initw_t
    ab["initw_binary"]["n_pep_widened"] = int(out[out["init_w"] >= initw_t]
                                              ["concat"].nunique())

    return {"n_pep": int(out["concat"].nunique()), "n_rows": int(len(out)),
            "crossover_niso": niso_x, "crossover_initw": initw_x,
            "by_niso": niso_strata, "by_initw": initw_strata, "ab": ab}


def run_line(line: str, niso_t: int, initw_t: int, out_root: Path) -> dict:
    rec = LINE_RECOVERY[line]
    ria = 0.0598
    print(f"\n=== {line.upper()} ===", flush=True)
    df = _load_concat(line, q_value=0.01, drop=rec["drop"])
    coeffs = load_aa_coefficients(CALIB_ROOT / line / rec["coeff"])
    out = _per_peptide(df, coeffs, ria)
    res = analyse(out, niso_t, initw_t)
    out_dir = out_root / line
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / "niso_crossover.csv", index=False)
    with (out_dir / "niso_crossover.json").open("w") as f:
        json.dump(res, f, indent=2)

    for key, strata, x in (("N_ISO", res["by_niso"], res["crossover_niso"]),
                           ("init_w", res["by_initw"], res["crossover_initw"])):
        cap = 12 if key == "N_ISO" else 7
        print(f"  {key:>6}  n_pep   MAE(iso0-3)  MAE(all)   Δ(iso03-all)")
        for b in sorted(strata):
            s = strata[b]
            tag = f"{cap}+" if b == cap else str(b)
            flag = "  <- iso0-3 worse" if s["delta_iso03_minus_all"] > 0 else ""
            print(f"  {tag:>6}  {s['n_pep']:>5}   {s['mae_iso03']:>10.4f}  "
                  f"{s['mae_all']:>9.4f}  {s['delta_iso03_minus_all']:>+9.4f}{flag}")
        print(f"  crossover {key} = {x}")
    ab = res["ab"]
    print(f"  A/B                within±0.05   wide-MAE")
    for k in ("flat_iso03", "flat_all", "niso", "initw_binary", "initw_graded"):
        print(f"    {k:14s} {ab[k]['within_0.05']*100:8.1f}%   {ab[k]['wide_mae']:.4f}")
    return res


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--line", choices=["ac16", "cm", "ipsc", "all"], default="ac16")
    parser.add_argument("--niso-threshold", type=int, default=11,
                        help="N_ISO at/above which to widen iso0-3 -> all [default 11]")
    parser.add_argument("--initw-threshold", type=int, default=6,
                        help="init (θ=0) width at/above which to widen [default 6, "
                        "the derived cross-line crossover]")
    parser.add_argument("--output-dir", type=Path,
                        default=REPO_ROOT / "runs" / "niso_crossover")
    args = parser.parse_args()
    lines = ["ac16", "cm", "ipsc"] if args.line == "all" else [args.line]
    for line in lines:
        run_line(line, args.niso_threshold, args.initw_threshold, args.output_dir)


if __name__ == "__main__":
    main()
