"""bench_mbr_quality.py — do MBR-transferred peaks carry real, sensible signal?

Run after ``riana integrate --mbr`` (e.g. ``runs/lve_atr_mbr``). Validates the
match-between-runs transfers *before* they reach the fit, on three questions:

1. **Survival.** How many transfers found a real MS1 apex (kept) vs were dropped
   by the graceful no-apex guard? (Dropped counts are parsed from the run log.)
2. **MS1 signal.** Do surviving MBR rows carry real intensity, or noise-floor
   junk? Compares the iso0 intensity + monoisotopic fraction m0 = iso0/Σiso of
   MBR rows against the directly-identified (q_value) rows.
3. **Does θ make sense?** In D₂O labeling m0 declines monotonically with labeling
   time, so a *good* transferred point should sit on its precursor's real
   trajectory. Two checks, each **benchmarked against held-out real points** (drop
   a real point, predict it from its real neighbours) so we learn whether MBR
   points are as trustworthy as real ones, not just "small":
   - **interp residual** — |observed m0 − linear-interp from the bracketing real
     points| for MBR vs held-out-real interior points.
   - **monotone corridor** — fraction of points that fall within the [next, prev]
     real m0 bracket (± tol), i.e. respect the decline, for MBR vs held-out-real.

m0 = iso0/Σiso0..5 is the integrate-stage proxy for θ (θ rises as m0 falls); the
labeling time is the ``timeNN`` token (= the day: 0,1,2,3,4,6,8,10,15,20,25,30).

Usage:
    python tests/benchmark/bench_mbr_quality.py --run runs/lve_atr_mbr
"""

from __future__ import annotations

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

_ISO = [f"iso{i}" for i in range(6)]


def _load(run_dir: str) -> pd.DataFrame:
    rows = []
    for path in sorted(glob.glob(os.path.join(run_dir, "*_riana.txt"))):
        base = os.path.basename(path)
        m = re.search(r"_(LVE|ATR)_time(\d+)_riana\.txt$", base)
        if not m:
            continue
        chamber, t = m.group(1), int(m.group(2))
        df = pd.read_csv(path, sep="\t", comment="#")
        cols = ["concat", "evidence", *_ISO]
        for c in ("apex_snr", "n_scans"):
            if c in df.columns:
                cols.append(c)
        df = df[cols].copy()
        for c in ("apex_snr", "n_scans"):
            if c not in df.columns:
                df[c] = np.nan  # pre-gate run
        df["chamber"] = chamber
        df["time"] = t
        rows.append(df)
    if not rows:
        raise SystemExit(f"no *_riana.txt under {run_dir}")
    out = pd.concat(rows, ignore_index=True)
    iso_sum = out[_ISO].sum(axis=1)
    out["m0"] = np.where(iso_sum > 0, out["iso0"] / iso_sum, np.nan)
    return out


def _parse_drops(run_dir: str) -> int:
    total = 0
    for log in glob.glob(os.path.join(run_dir, "*.log")):
        with open(log, errors="ignore") as fh:
            for line in fh:
                m = re.search(r"dropped (\d+) MBR transfer", line)
                if m:
                    total += int(m.group(1))
    return total


def _pct(x):
    return f"{100 * x:.1f}%"


def survival(df: pd.DataFrame, run_dir: str, planned: int = 0) -> None:
    print("\n" + "=" * 70)
    print("1. SURVIVAL — transfers that found a real MS1 apex vs were dropped")
    print("=" * 70)
    for ch, g in df.groupby("chamber"):
        n_real = int((g["evidence"] == "q_value").sum())
        n_mbr = int((g["evidence"] == "mbr").sum())
        print(f"  {ch}: {n_real:,} real rows | {n_mbr:,} MBR rows survived")
    surv = int((df["evidence"] == "mbr").sum())
    # The per-run no-apex drop count is logged inside worker processes, so it
    # doesn't reach the main logfile (_parse_drops will usually find 0 until that
    # is surfaced on the result object). Prefer planned-minus-surviving: pass
    # --planned (the total transfer count from the plan step).
    planned = planned or _parse_drops(run_dir)
    if planned:
        dropped = planned - surv
        print(f"  planned transfers: {planned:,}  |  survived: {surv:,}  |  "
              f"dropped (no apex): {dropped:,} = {_pct(dropped / planned)}")
    else:
        print(f"  survived: {surv:,}  (pass --planned <N> from the plan step for "
              "the drop rate; worker-logged drops don't reach the main logfile)")


def ms1_signal(df: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("2. MS1 SIGNAL — surviving MBR rows vs real rows")
    print("=" * 70)
    for ev in ("q_value", "mbr"):
        g = df[df["evidence"] == ev]
        i0 = g["iso0"].to_numpy(dtype=float)
        frac_zero = float(np.mean(i0 <= 0))
        q = np.percentile(i0[i0 > 0], [10, 50, 90]) if (i0 > 0).any() else [0, 0, 0]
        m0 = g["m0"].dropna()
        print(f"  {ev:<8} n={len(g):>7,}  iso0 p10/50/90 = "
              f"{q[0]:,.0f} / {q[1]:,.0f} / {q[2]:,.0f}  "
              f"frac iso0<=0 = {_pct(frac_zero)}  "
              f"m0 median = {m0.median():.3f}")
    print("  (MBR rows are expected lower-intensity — they were missed because "
          "scarce — but should NOT be ~0; frac iso0<=0 ≈ 0 confirms the drop "
          "kept only real peaks.)")


def _bracketed_residuals(traj: pd.DataFrame):
    """For one precursor trajectory, return (mbr_resid, real_loo_resid, corridor)
    lists: interp residual + monotone-corridor membership for MBR points and for
    held-out (leave-one-out) interior real points."""
    traj = traj.sort_values("time")
    real = traj[traj["evidence"] == "q_value"]
    if len(real) < 3:
        return [], [], [], []
    rt = real["time"].to_numpy(float)
    rm = real["m0"].to_numpy(float)

    def at(point_t, point_m0, donor_t, donor_m0):
        # bracketing real donors around point_t
        lo = donor_t < point_t
        hi = donor_t > point_t
        if not lo.any() or not hi.any():
            return None
        tp, mp = donor_t[lo][-1], donor_m0[lo][-1]   # nearest before
        tn, mn = donor_t[hi][0], donor_m0[hi][0]     # nearest after
        pred = mp + (mn - mp) * (point_t - tp) / (tn - tp)
        resid = abs(point_m0 - pred)
        tol = 0.05
        corridor = (min(mp, mn) - tol) <= point_m0 <= (max(mp, mn) + tol)
        return resid, corridor

    mbr_res, mbr_cor = [], []
    for _, r in traj[traj["evidence"] == "mbr"].iterrows():
        got = at(r["time"], r["m0"], rt, rm)
        if got and np.isfinite(r["m0"]):
            mbr_res.append(got[0]); mbr_cor.append(got[1])

    loo_res, loo_cor = [], []
    for i in range(len(real)):
        mask = np.ones(len(real), bool); mask[i] = False
        got = at(rt[i], rm[i], rt[mask], rm[mask])
        if got and np.isfinite(rm[i]):
            loo_res.append(got[0]); loo_cor.append(got[1])
    return mbr_res, loo_res, mbr_cor, loo_cor


def trajectory_sense(df: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("3. DOES θ MAKE SENSE — MBR points vs held-out real points")
    print("   (m0 declines with labeling time; a good MBR point sits on the curve)")
    print("=" * 70)
    for ch, g in df.groupby("chamber"):
        mbr_res, loo_res, mbr_cor, loo_cor = [], [], [], []
        for _, traj in g.groupby("concat"):
            a, b, c, d = _bracketed_residuals(traj)
            mbr_res += a; loo_res += b; mbr_cor += c; loo_cor += d
        if not mbr_res:
            print(f"  {ch}: no MBR points with bracketing real neighbours.")
            continue
        print(f"\n  {ch}:")
        print(f"    interp |Δm0|  MBR   median={np.median(mbr_res):.3f} "
              f"p90={np.percentile(mbr_res, 90):.3f}  (n={len(mbr_res):,})")
        print(f"    interp |Δm0|  REAL  median={np.median(loo_res):.3f} "
              f"p90={np.percentile(loo_res, 90):.3f}  (n={len(loo_res):,}, "
              f"held-out baseline)")
        print(f"    monotone corridor  MBR={_pct(np.mean(mbr_cor))}  "
              f"REAL={_pct(np.mean(loo_cor))}")
    print("\n  Read: MBR |Δm0| ≈ REAL |Δm0| and corridor% ≈ REAL ⇒ transferred "
          "points are as trustworthy as real ones. MBR ≫ REAL ⇒ noise/mis-transfer.")


def _mbr_points_with_neighbours(df: pd.DataFrame) -> pd.DataFrame:
    """Per surviving MBR point bracketed by real neighbours: run, snr, iso0,
    monotone-corridor membership, interp residual."""
    recs = []
    for (ch, concat), g in df.groupby(["chamber", "concat"]):
        g = g.sort_values("time")
        real = g[g["evidence"] == "q_value"]
        if len(real) < 3:
            continue
        rt = real["time"].to_numpy(float)
        rm = real["m0"].to_numpy(float)
        for r in g[g["evidence"] == "mbr"].itertuples():
            lo, hi = rt < r.time, rt > r.time
            if not lo.any() or not hi.any() or not np.isfinite(r.m0):
                continue
            tp, mp = rt[lo][-1], rm[lo][-1]
            tn, mn = rt[hi][0], rm[hi][0]
            pred = mp + (mn - mp) * (r.time - tp) / (tn - tp)
            cor = (min(mp, mn) - 0.05) <= r.m0 <= (max(mp, mn) + 0.05)
            ns = int(r.n_scans) if np.isfinite(r.n_scans) else 0
            recs.append((f"{ch}_t{int(r.time):02d}", float(r.apex_snr),
                         float(r.iso0), ns, bool(cor), abs(float(r.m0) - pred)))
    return pd.DataFrame(recs, columns=["run", "snr", "iso0", "n_scans", "cor", "resid"])


_SCAN_GRID = (0, 3, 5, 7)
_SNR_GRID = (0, 4, 6, 8)


def gate_sweep(df: pd.DataFrame) -> None:
    """Two-part MBR gate sweep: min nonzero scans (N) × min apex-SNR (T).

    An **inf** apex_snr (sparse MAD=0 trace, no noise floor) is treated as FAIL —
    it is not a defined SNR. Each cell is keep% / corridor% on MBR points with
    bracketing real neighbours (corridor = θ-recovery proxy; real held-out
    baseline ~94%). N>=0, T=0 is the no-gate baseline.
    """
    print("\n" + "=" * 70)
    print("5. TWO-PART GATE SWEEP — min nonzero scans (N) × min apex-SNR (T)")
    print("   inf-SNR (sparse, no noise floor) = FAIL; cell = keep% / corridor%")
    print("=" * 70)
    R = _mbr_points_with_neighbours(df)
    if R["snr"].isna().all() or R["n_scans"].fillna(0).eq(0).all():
        print("  no apex_snr / n_scans columns — re-run integrate so the gate "
              "diagnostics are written.")
        return
    R["eff_snr"] = np.where(np.isfinite(R["snr"]), R["snr"], 0.0)  # inf -> fail
    n_tot = len(R)
    print(f"  {n_tot:,} MBR points w/ neighbours | inf-SNR "
          f"{np.isinf(R['snr']).mean():.0%} | baseline corridor {R['cor'].mean():.0%}")
    print("  " + "N\\T".rjust(6) + "".join(("T>=" + str(t)).rjust(13) for t in _SNR_GRID))
    cand = []
    for n in _SCAN_GRID:
        cells = []
        for t in _SNR_GRID:
            k = R[(R["n_scans"] >= n) & (R["eff_snr"] >= t)]
            cells.append(f"{len(k) / n_tot:3.0%}/{k['cor'].mean():3.0%}"
                         if len(k) else "-/-")
            if len(k):
                cand.append((len(k) / n_tot, float(k["cor"].mean()), n, t,
                             float(k["resid"].median())))
        print("  " + ("N>=" + str(n)).rjust(6) + "".join(c.rjust(13) for c in cells))

    # conservative pick: closest to ~20% keep, then best corridor
    near = [c for c in cand if 0.12 <= c[0] <= 0.28]
    if near:
        best = max(near, key=lambda c: c[1])
        print(f"\n  conservative pick (~20% keep): --mbr-min-scans {best[2]} "
              f"--mbr-min-snr {best[3]} → keep {best[0]:.0%}, corridor {best[1]:.0%}, "
              f"|Δm0|med {best[4]:.3f}  (real ~94% / 0.013)")


def real_vs_mbr_gate(df: pd.DataFrame) -> None:
    """Calibration: would the gate discard genuine (q-value) IDs too? A gate that
    fails many *real* points is an abundance filter, not an MBR-junk filter; the
    gap (MBR fail% ≫ real fail%) is the discrimination. Applied over ALL rows (no
    neighbour requirement) — real rows are never gated in production, this is the
    hypothetical comparison."""
    print("\n" + "=" * 70)
    print("6. GATE vs REAL POINTS — does the gate also discard confident IDs?")
    print("=" * 70)
    if df["apex_snr"].isna().all() or df["n_scans"].fillna(0).eq(0).all():
        print("  no apex_snr / n_scans — re-run integrate so the columns exist.")
        return
    print(f"  {'pop':<6}{'n':>10}{'inf-SNR':>10}"
          + "".join(f"{f'scans<{n}':>11}" for n in _SCAN_GRID[1:]))
    for ev, label in (("q_value", "real"), ("mbr", "MBR")):
        g = df[df["evidence"] == ev]
        inf = float(np.isinf(g["apex_snr"]).mean())
        cells = "".join(f"{np.mean(g['n_scans'] < n):>11.0%}" for n in _SCAN_GRID[1:])
        print(f"  {label:<6}{len(g):>10,}{_pct(inf):>10}{cells}")
    print("\n  combined-gate fail% (drop if n_scans<N OR inf-SNR OR apex_snr<T):")
    print(f"  {'cut':<16}{'real fail':>12}{'MBR fail':>12}{'ratio':>9}")
    for n, t in ((3, 4), (5, 6), (7, 8)):
        fr = {}
        for ev in ("q_value", "mbr"):
            g = df[df["evidence"] == ev]
            eff = np.where(np.isfinite(g["apex_snr"]), g["apex_snr"], 0.0)
            fr[ev] = float(((g["n_scans"] < n) | (eff < t)).mean())
        ratio = fr["mbr"] / fr["q_value"] if fr["q_value"] else float("inf")
        print(f"  N>={n},T>={t:<10}{_pct(fr['q_value']):>12}{_pct(fr['mbr']):>12}"
              f"{ratio:>8.1f}x")
    print("  (MBR fail% ≫ real fail% ⇒ the gate separates transfers from genuine "
          "signal; if close, it is mostly an abundance filter.)")


def examples(df: pd.DataFrame, n: int = 5) -> None:
    print("\n" + "=" * 70)
    print(f"4. EXAMPLE TRAJECTORIES (m0 by time; * = MBR)")
    print("=" * 70)
    g = df[df["chamber"] == "LVE"]
    # precursors with >=2 MBR points and >=4 real, for a visible curve
    cand = []
    for concat, traj in g.groupby("concat"):
        n_mbr = int((traj["evidence"] == "mbr").sum())
        n_real = int((traj["evidence"] == "q_value").sum())
        if n_mbr >= 2 and n_real >= 4:
            cand.append(concat)
    for concat in cand[:n]:
        traj = g[g["concat"] == concat].sort_values("time")
        cells = [f"{int(r.time):>2}:{r.m0:.2f}{'*' if r.evidence == 'mbr' else ' '}"
                 for r in traj.itertuples()]
        print(f"  {concat:<22} " + " ".join(cells))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default="runs/lve_atr_mbr")
    ap.add_argument("--planned", type=int, default=0,
                    help="total transfers from the plan step, for the drop rate")
    args = ap.parse_args()
    df = _load(args.run)
    survival(df, args.run, args.planned)
    ms1_signal(df)
    trajectory_sense(df)
    gate_sweep(df)
    real_vs_mbr_gate(df)
    examples(df)


if __name__ == "__main__":
    main()
