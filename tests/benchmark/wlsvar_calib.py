"""WP3 calibration targets: the same statistics on a real fit and on a simulated one.

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md`` WP3. The simulation
(``bench_rising_null.py --design atf6``) is tuned until these match ATF6. Every statistic
is computed from the frame the shipped linear rollup collapses (its gate, Spep floor,
admission and parsimony), so real and simulated data go through the same estimator.

- *Design and gating*: tests (protein × experiment with a Δk), admitted precursors per
  protein × condition (median, 90th percentile), protein k (median, 10th / 90th).
- *Completeness*: per protein × condition with ≥ 2 precursors, the share of t > 0 cells
  holding every admitted precursor (median), and the share of those groups complete in
  every cell.
- *Precision*: median band σ (= band / 3.29) at θ̄ ≤ 0.05 and at θ̄ 0.4–0.7.
- *Peptide spread* (WP1b): DerSimonian–Laird τ per θ̄ bin, as τ_rel of an apparent-k
  offset, and Q/df at baseline.
- *Animal* (WP1d, tissues of the same animals): cross-tissue correlation of the per-cell
  deviation (φ scale, centred within the timepoint), its variance and covariance on the
  relative-k scale, the correlation of the proteome-wide animal shift across tissues, and
  the SD of that shift over the last three timepoints.
- *Timepoint deviation shared across tissues* (cross-sectional designs): per protein ×
  tissue × day, the mean relative-k residual over the day's animals (both conditions);
  the share of its variance that the tissues have in common, after removing the animals'
  own cross-tissue covariance ÷ n (``day_shared_frac``; ATF6 0.85).

python tests/benchmark/wlsvar_calib.py runs/wlsvar_eval/atf6/fit --ref control --test ATF6-KO
python tests/benchmark/wlsvar_calib.py runs/wlsvar_calib/try1/ratio1_0/theta --ref A --test B
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from riana.core import protein  # noqa: E402
from tests.benchmark.eval_rising_variance import collapse, dk_shipped  # noqa: E402
from tests.benchmark.wlsvar_wp1 import (animal_frame, capture, residuals,  # noqa: E402
                                        tau_frame)


def completeness(frame: pd.DataFrame) -> tuple[float, float]:
    x = frame[(frame.labeling_time > 0) & np.isfinite(frame.fs)]
    rows = []
    for _, g in x.groupby(["experiment", "condition", "protein"]):
        if g.concat.nunique() < 2:
            continue
        m = g.pivot_table(index=["labeling_time", "biological_replicate"], columns="concat",
                          values="fs", aggfunc="size").notna()
        rows.append(m.all(axis=1).mean())
    c = pd.Series(rows)
    return float(c.median()), float((c == 1).mean())


def animal_shift_sd(seen: dict) -> float:
    r = residuals(collapse(seen, False))
    r = r[r.individual.astype(str) != ""].copy()
    r["rk"] = (r.phi - r.fit) / r.fit
    lo, hi = r.rk.quantile([0.01, 0.99])
    r = r[(r.rk >= lo) & (r.rk <= hi)]
    cell = ["protein", "condition", "labeling_time", "experiment"]
    r["rk_c"] = r.rk - r.groupby(cell).rk.transform("mean")
    sh = r.groupby(["condition", "labeling_time", "individual"]).rk_c.median()
    sd = sh.groupby("labeling_time").std()
    return float(sd.iloc[-3:].mean())


def day_shared_frac(seen: dict) -> float:
    """Share of the per-protein × tissue × day deviation that the tissues have in common.

    Per cell (protein, condition, day, tissue, animal): relative-k residual rk = (φ − fit) / fit
    from the shipped per-condition wls line. The day mean d(protein, tissue, day) over the
    animals of both conditions carries the tissue-shared day term plus the tissue-shared animal
    term ÷ n; the animal part is removed with the within-cell cross-tissue covariance of the
    animal deviations. Mean over tissue pairs."""
    r = residuals(collapse(seen, False))
    r = r[r.individual.astype(str) != ""].copy()
    r["rk"] = (r.phi - r.fit) / r.fit
    lo, hi = r.rk.quantile([0.01, 0.99])
    r = r[(r.rk >= lo) & (r.rk <= hi)]
    day = ["protein", "experiment", "labeling_time"]
    r["a"] = r.rk - r.groupby(day + ["condition"]).rk.transform("mean")
    r = r[r.groupby(day).rk.transform("size") >= 4]
    d = r.groupby(day).agg(d=("rk", "mean"), n=("rk", "size")).reset_index()
    wd = d.pivot_table(index=["protein", "labeling_time"], columns="experiment", values="d")
    wa = r.pivot_table(index=["protein", "condition", "labeling_time", "individual"],
                       columns="experiment", values="a")
    nbar = float(d.n.mean())
    m = nbar / 2                      # within-cell centring removes 1/m of the animal variance
    tissues = sorted(r.experiment.unique())
    fr = []
    for i, a in enumerate(tissues):
        for b in tissues[i + 1:]:
            x, y = wd[[a, b]].dropna(), wa[[a, b]].dropna()
            cov_d, va, vb = (x[a] * x[b]).mean(), (x[a] ** 2).mean(), (x[b] ** 2).mean()
            cov_a, vaa, vab = (float((y[u] * y[v]).mean()) * m / (m - 1)
                               for u, v in ((a, b), (a, a), (b, b)))
            fr.append((cov_d - cov_a / nbar) / np.sqrt((va - vaa / nbar) * (vb - vab / nbar)))
    return float(np.mean(fr))


def stats(fit_dir: Path, ref: str, test: str) -> dict:
    seen = capture(fit_dir, ref, test)
    f = seen["frame"]
    out: dict = {}
    long = collapse(seen, False)
    dk = dk_shipped(long, "wls", ref, test)
    out["tests"] = len(dk)
    out["k_med"], out["k_q10"], out["k_q90"] = (float(x) for x in
                                                dk.k_A.quantile([.5, .1, .9]))
    npre = f.groupby(["experiment", "condition", "protein"]).concat.nunique()
    npre = npre[npre >= 2]
    out["precursors_med"], out["precursors_q90"] = float(npre.median()), float(npre.quantile(.9))
    out["complete_cells_med"], out["fully_complete"] = completeness(f)
    sig = (f.fs_upper - f.fs_lower) / protein._PI_SPAN_SIGMA
    th = f.fs
    out["sigma_base"] = float(sig[(th <= 0.05) & (f.labeling_time > 0)].median())
    out["sigma_late"] = float(sig[(th >= 0.4) & (th <= 0.7)].median())
    tb, ts = tau_frame(seen)
    out["Q_df_base"] = float(tb.Q_over_df.iloc[0])
    for lo, name in ((0.2, "tau_rel_25"), (0.4, "tau_rel_45"), (0.6, "tau_rel_65")):
        row = tb[tb.theta_bin == f"({lo}, {round(lo + 0.1, 1)}]"]
        out[name] = float(row.tau_rel_k.iloc[0]) if len(row) else np.nan
    out["tau2_over_v_late"] = float(ts.tau2_over_v_late.iloc[0])
    if f.experiment.nunique() > 1 and "individual" in f and (f.individual.astype(str) != "").any():
        pairs, prot = animal_frame(seen)
        out["corr_e"] = float(pairs.corr_e.mean())
        out["var_r"] = float(pairs[["var_r_first", "var_r_second"]].to_numpy().mean())
        out["cov_r"] = float(pairs.cov_r.mean())
        out["corr_shift"] = float(prot.corr_proteome_shift.mean())
        out["shift_sd_late"] = animal_shift_sd(seen)
        out["day_shared_frac"] = day_shared_frac(seen)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("fit_dirs", type=Path, nargs="+")
    p.add_argument("--ref", required=True)
    p.add_argument("--test", required=True)
    p.add_argument("-o", "--out", type=Path)
    a = p.parse_args()
    rows = [dict(fit=str(d), **stats(d, a.ref, a.test)) for d in a.fit_dirs]
    df = pd.DataFrame(rows).set_index("fit").T
    with pd.option_context("display.width", 220, "display.float_format", "{:.4g}".format):
        print(df.to_string())
    if a.out:
        df.to_csv(a.out, sep="\t")


if __name__ == "__main__":
    main()
