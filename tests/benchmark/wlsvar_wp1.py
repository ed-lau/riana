"""WP1 diagnostics for the wls-var / D / E evaluation (no new estimators).

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md`` WP1.

``q1 FIT [--ref R --test T]`` — does the cell variance Ṽ that ``wls-var`` weights on track
the real scatter? On the shipped collapse (and the shared-proportion collapse C when there
are two conditions), each protein is fitted as shipped (``wls``: through-origin slope per
condition, OLS then one step with weights (1 − θ̂)² from the fitted value). Each cell's
residual is put on the θ scale, u = (φ − φ̂)(1 − θ̂), and divided by (1 − h) (its leverage),
so under ``wls`` E[u²/(1 − h)] is the same for every cell of a protein.

- Slope of log(u²/(1 − h)) on log V within protein (protein fixed effects): 1 = the variance
  tracks the scatter, 0 = uninformative. Given for the raw cell variance V = 1/Σ(1/σ²)
  (attenuated by its own noise) and for the moderated Ṽ that ``wls-var`` uses (eBayes
  posterior; a calibrated Ṽ has slope 1).
- The same within-protein regressions on log(peptides in the cell): Ṽ falls as 1/n
  (slope −1) by construction; if the scatter falls more slowly, the cells have a floor
  (peptide spread that Ṽ ignores).
- Mean scatter by Ṽ decile (within-protein, relative to the protein mean), for the shape.

``animal FIT --ref R --test T`` (WP1d) — the animal term, from tissues of the same animal
(experiments = tissues, ``individual`` = the animal). Residuals of the shipped ``wls`` fit are
centred within (protein, condition, timepoint, tissue), so curve shape and the timepoint mean
drop out and what is left is each animal's deviation. Per tissue pair, pooled over proteins
and timepoints: the correlation and covariance of the deviations on two scales,
e = φ − φ̂ (additive in φ) and r = e / φ̂ (relative, i.e. a shift in log k: proportional).
Under a proportional animal effect the covariance is the same for every pair on the r scale
and grows with the two tissues' k on the e scale; under an additive one, the reverse. Also the
proteome-wide part: each animal's median r over proteins, correlated across its tissues.

``tau FIT [--ref R --test T]`` (WP1b–c) — D's between-peptide variance τ². Per collapsed cell
with ≥ 2 peptides, Cochran's Q = Σ w (θ_i − θ̄)² with w = 1/σ² from the peptides' bands (the
collapse's own weights). Pooled by θ̄ bin, DerSimonian–Laird: τ² = Σ(Q − (n − 1)) / Σ(Σw − Σw²/Σw).
Its shape is compared with an apparent-k offset per peptide, τ ∝ g(θ) = (1 − θ)(−ln(1 − θ)), and
with a θ-proportional offset, τ ∝ θ: the implied τ_rel = τ / g (or τ / θ) is flat across bins
for the right one. (c) τ² against the cell variance Ṽ (median per bin): where τ² ≪ Ṽ at the
late, slope-carrying points, D reduces to C plus a cell variance.

python tests/benchmark/wlsvar_wp1.py q1 runs/wlsvar_eval/boomi/fit
python tests/benchmark/wlsvar_wp1.py q1 runs/wlsvar_eval/lve_atr/fit --ref control --test atrium
python tests/benchmark/wlsvar_wp1.py tau runs/wlsvar_eval/boomi/fit
python tests/benchmark/wlsvar_wp1.py animal runs/wlsvar_eval/atf6/fit --ref control --test ATF6-KO
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from riana.core import protein  # noqa: E402
from riana.core.linear_model import _fit_fdist  # noqa: E402
from tests.benchmark.eval_rising_variance import collapse  # noqa: E402

PHI_LIMIT = -4.0
CELL = ["experiment", "condition", "protein", "labeling_time", "biological_replicate"]


def capture(fit_dir: Path, ref: str | None, test: str | None) -> dict:
    """Run the shipped linear rollup and keep the gated frame it collapses."""
    from typer.testing import CliRunner
    from riana.cli import app

    seen: dict = {}
    orig = protein._collapse_long

    def spy(fractions, *, method, min_peptides):
        seen.update(frame=fractions.copy(), method=method, min_peptides=min_peptides)
        return orig(fractions, method=method, min_peptides=min_peptides)

    args = ["rollup", str(fit_dir), "--model", "linear", "--linear-weights", "wls"]
    if ref and test:
        args += ["--reference-condition", ref, "--test-condition", test]
    protein._collapse_long = spy
    try:
        with tempfile.TemporaryDirectory() as tmp:
            res = CliRunner().invoke(app, [*args, "-o", tmp])
            if res.exit_code != 0:
                raise RuntimeError(f"rollup failed on {fit_dir}:\n{res.output[-2000:]}")
    finally:
        protein._collapse_long = orig
    return seen


def residuals(long: pd.DataFrame) -> pd.DataFrame:
    """Shipped ``wls`` fit per (experiment, protein); θ-scale residual and leverage per cell."""
    c = long.copy()
    c["phi"] = np.log1p(-c.theta.clip(upper=1 - 1e-12))
    c = c[(c.labeling_time > 0) & np.isfinite(c.phi) & (c.phi > PHI_LIMIT)].copy()
    out = []
    for _, g in c.groupby(["experiment", "protein"], sort=False):
        t, y = g.labeling_time.to_numpy(float), g.phi.to_numpy(float)
        cond = g.condition.to_numpy()
        fit = np.empty_like(y)
        h = np.empty_like(y)
        ok = True
        for cc in np.unique(cond):  # per-condition through-origin slopes are separable
            m = cond == cc
            tt, yy = t[m], y[m]
            if len(tt) < 2:
                ok = False
                break
            b = (tt @ yy) / (tt @ tt)                       # OLS
            w = np.exp(2 * np.maximum(b * tt, PHI_LIMIT))   # (1 − θ̂)² from the fitted value
            b = (w * tt) @ yy / ((w * tt) @ tt)             # one WLS step
            fit[m] = b * tt
            h[m] = w * tt ** 2 / ((w * tt) @ tt)
        if not ok or len(y) - len(np.unique(cond)) < 2:
            continue
        g = g.assign(fit=fit, h=h)
        g["u"] = (g.phi - g.fit) * np.exp(np.maximum(g.fit, PHI_LIMIT))
        out.append(g)
    return pd.concat(out, ignore_index=True)


def _fe_slope(y: pd.Series, x: pd.Series, grp: pd.Series) -> tuple[float, float, int]:
    """Within-group OLS slope of y on x, its SE, and the number of points used."""
    d = pd.DataFrame({"y": y, "x": x, "g": grp}).replace([np.inf, -np.inf], np.nan).dropna()
    d = d[d.groupby("g").x.transform("size") > 1]
    yd = d.y - d.groupby("g").y.transform("mean")
    xd = d.x - d.groupby("g").x.transform("mean")
    sxx = float(xd @ xd)
    if sxx <= 0:
        return float("nan"), float("nan"), len(d)
    b = float(xd @ yd) / sxx
    e = yd - b * xd
    dof = len(d) - d.g.nunique() - 1
    return b, float(np.sqrt((e @ e) / dof / sxx)), len(d)


def q1(fit_dir: Path, ref: str | None, test: str | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    seen = capture(fit_dir, ref, test)
    f = seen["frame"]
    npep = f.groupby(CELL).concat.nunique().rename("n_pep").reset_index()
    two = f.condition.nunique() >= 2
    rows, deciles = [], []
    for name, shared in [("A", False)] + ([("C", True)] if two else []):
        long = collapse(seen, shared)
        v, dfv = long.theta_var.to_numpy(float), long.theta_df.to_numpy(float)
        ok = np.isfinite(v) & (v > 0) & np.isfinite(dfv) & (dfv > 0)
        d0, s0 = _fit_fdist(v[ok], dfv[ok])
        long["v_mod"] = (d0 * s0 + dfv * v) / (d0 + dfv)
        r = residuals(long).merge(npep, on=CELL, how="left")
        r = r[np.isfinite(r.theta_var) & (r.theta_var > 0) & (r.h < 1)]
        r["s2"] = r.u ** 2 / (1 - r.h)
        grp = r.experiment.astype(str) + "|" + r.protein.astype(str)
        ly = np.log(r.s2.where(r.s2 > 0))
        for lab, x in [("log V (raw)", np.log(r.theta_var)), ("log Ṽ (moderated)", np.log(r.v_mod)),
                       ("log n_pep", np.log(r.n_pep))]:
            b, se, n = _fe_slope(ly, x, grp)
            rows.append(dict(collapse=name, outcome="log scatter", regressor=lab, slope=b, se=se,
                             cells=n, d0=d0))
        if two:
            # within a condition (protein × condition fixed effects) vs between conditions:
            # per protein, the test − reference difference of condition means
            pc = grp + "|" + r.condition.astype(str)
            b, se, n = _fe_slope(ly, np.log(r.v_mod), pc)
            rows.append(dict(collapse=name, outcome="log scatter", regressor="log Ṽ within condition",
                             slope=b, se=se, cells=n, d0=d0))
            m = pd.DataFrame({"g": grp, "c": r.condition.astype(str), "y": ly,
                              "x": np.log(r.v_mod)}).dropna().groupby(["g", "c"])[["y", "x"]].mean()
            m = m.unstack("c").dropna()
            dy = m["y"][test] - m["y"][ref]
            dx = m["x"][test] - m["x"][ref]
            X = np.column_stack([np.ones(len(dx)), dx])
            beta, *_ = np.linalg.lstsq(X, dy, rcond=None)
            e = dy - X @ beta
            cov = np.linalg.inv(X.T @ X) * float(e @ e) / (len(dy) - 2)
            rows.append(dict(collapse=name, outcome="Δ log scatter (test − ref)",
                             regressor="Δ log Ṽ between conditions", slope=beta[1],
                             se=float(np.sqrt(cov[1, 1])), cells=len(dy), d0=d0))
        for lab, x in [("log n_pep", np.log(r.n_pep))]:
            b, se, n = _fe_slope(np.log(r.v_mod), x, grp)
            rows.append(dict(collapse=name, outcome="log Ṽ", regressor=lab, slope=b, se=se,
                             cells=n, d0=d0))
            b, se, n = _fe_slope(np.log(r.theta_var), x, grp)
            rows.append(dict(collapse=name, outcome="log V (raw)", regressor=lab, slope=b, se=se,
                             cells=n, d0=d0))
        # shape: within-protein relative Ṽ decile vs relative scatter (means, not logs)
        rel_v = r.v_mod / r.groupby(grp).v_mod.transform("mean")
        rel_s = r.s2 / r.groupby(grp).s2.transform("mean")
        q = pd.qcut(np.log(rel_v), 10, labels=False, duplicates="drop")
        dd = pd.DataFrame({"decile": q, "rel_v": rel_v, "rel_scatter": rel_s}).groupby(
            "decile").agg(cells=("rel_v", "size"), rel_v=("rel_v", "mean"),
                          rel_scatter=("rel_scatter", "mean"))
        deciles.append(dd.reset_index().assign(collapse=name))
    return pd.DataFrame(rows), pd.concat(deciles, ignore_index=True)


def tau(fit_dir: Path, ref: str | None, test: str | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    seen = capture(fit_dir, ref, test)
    f = seen["frame"].copy()
    f = f[np.isfinite(f.fs)]
    sig = (f.fs_upper - f.fs_lower) / protein._PI_SPAN_SIGMA
    sig = sig.where(np.isfinite(sig) & (sig > 0))
    f["sig"] = sig.fillna(sig.groupby([f[c] for c in CELL]).transform("median"))
    f = f[np.isfinite(f.sig)]
    f["w"] = 1.0 / f.sig ** 2
    f["wx"] = f.w * f.fs
    f["w2"] = f.w ** 2
    g = f.groupby(CELL)
    c = pd.DataFrame({"n": g.size(), "sw": g.w.sum(), "swx": g.wx.sum(), "sw2": g.w2.sum()})
    c["theta"] = c.swx / c.sw
    f = f.join(c.theta.rename("cell_theta"), on=CELL)
    f["q"] = f.w * (f.fs - f.cell_theta) ** 2
    c["Q"] = f.groupby(CELL).q.sum()
    c = c[c.n >= 2].reset_index()
    c["den"] = c.sw - c.sw2 / c.sw
    c["v"] = 1.0 / c.sw
    bins = [-1, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2]
    c["bin"] = pd.cut(c.theta, bins)
    rows = []
    for b, d in c.groupby("bin", observed=True):
        t2 = max(0.0, float((d.Q - (d.n - 1)).sum() / d.den.sum()))
        th = float(d.theta.median())
        gk = (1 - th) * -np.log1p(-min(max(th, 1e-6), 1 - 1e-6))
        rows.append(dict(theta_bin=str(b), cells=len(d), theta_median=th,
                         Q_over_df=float(d.Q.sum() / (d.n - 1).sum()), tau=np.sqrt(t2),
                         tau_rel_k=np.sqrt(t2) / gk if gk > 0 else np.nan,
                         tau_rel_theta=np.sqrt(t2) / th if th > 0 else np.nan,
                         v_cell_median=float(d.v.median()),
                         tau2_over_v=t2 / float(d.v.median()),
                         peptides_median=float(d.n.median())))
    late = c[c.labeling_time >= c.groupby(["experiment", "condition", "protein"]).labeling_time
             .transform(lambda x: x.quantile(0.5))]
    t2l = max(0.0, float((late.Q - (late.n - 1)).sum() / late.den.sum()))
    summ = pd.DataFrame([dict(cells=len(c), late_cells=len(late), tau2_late=t2l,
                              v_late_median=float(late.v.median()),
                              tau2_over_v_late=t2l / float(late.v.median()),
                              share_cells_tau2_gt_v=float((t2l > late.v).mean()))])
    return pd.DataFrame(rows), summ


def animal(fit_dir: Path, ref: str, test: str, shared: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    seen = capture(fit_dir, ref, test)
    r = residuals(collapse(seen, shared))
    r = r[r.individual.astype(str) != ""].copy()
    r["e"] = r.phi - r.fit
    r["r"] = r.e / r.fit
    cell = ["protein", "condition", "labeling_time", "experiment"]
    for v in ("e", "r"):
        r[v] = r[v] - r.groupby(cell)[v].transform("mean")
    n_t = r.groupby(cell).individual.transform("size")
    r = r[n_t >= 2]
    k = (-r.fit / r.labeling_time).groupby([r.experiment, r.protein]).median()
    tissues = sorted(r.experiment.unique())
    wide = {v: r.pivot_table(index=["protein", "condition", "labeling_time", "individual"],
                             columns="experiment", values=v) for v in ("e", "r")}
    rows = []
    for i, a in enumerate(tissues):
        for b in tissues[i + 1:]:
            row = dict(pair=f"{a}–{b}")
            kk = pd.concat([k[a], k[b]], axis=1, keys=["a", "b"]).dropna()
            row["k_ratio_median"] = float(np.median(kk.b / kk.a))
            for v in ("e", "r"):
                w = wide[v][[a, b]].dropna()
                if v == "r":  # relative deviations explode where φ̂ ≈ 0: trim to the central 98 %
                    lo, hi = w.stack().quantile([0.01, 0.99])
                    w = w[(w >= lo).all(axis=1) & (w <= hi).all(axis=1)]
                row[f"n_{v}"] = len(w)
                row[f"corr_{v}"] = float(w[a].corr(w[b]))
                row[f"cov_{v}"] = float(np.mean(w[a] * w[b]))
                row[f"var_{v}_{'first'}"] = float(np.mean(w[a] ** 2))
                row[f"var_{v}_{'second'}"] = float(np.mean(w[b] ** 2))
            rows.append(row)
    # proteome-wide animal shift: median relative deviation over proteins, per animal × tissue
    shift = r.groupby(["individual", "experiment"]).r.median().unstack()
    c = shift.corr()
    prows = [dict(pair=f"{a}–{b}", animals=int(shift[[a, b]].dropna().shape[0]),
                  corr_proteome_shift=float(c.loc[a, b]))
             for i, a in enumerate(tissues) for b in tissues[i + 1:]]
    return pd.DataFrame(rows), pd.DataFrame(prows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("test", choices=["q1", "animal", "tau"])
    p.add_argument("--shared", action="store_true", help="animal: on collapse C")
    p.add_argument("fit_dir", type=Path)
    p.add_argument("--ref")
    p.add_argument("--test-condition", "--test", dest="cond")
    p.add_argument("-o", "--out", type=Path, help="write the tables here (TSV prefix)")
    a = p.parse_args()
    if a.test == "animal":
        reg, dec = animal(a.fit_dir, a.ref, a.cond, a.shared)
    elif a.test == "tau":
        reg, dec = tau(a.fit_dir, a.ref, a.cond)
    else:
        reg, dec = q1(a.fit_dir, a.ref, a.cond)
    with pd.option_context("display.width", 200, "display.max_rows", 200,
                           "display.float_format", "{:.3f}".format):
        print(reg.to_string(index=False))
        print()
        print(dec.to_string(index=False))
    if a.out:
        reg.to_csv(f"{a.out}_slopes.tsv", sep="\t", index=False)
        dec.to_csv(f"{a.out}_deciles.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
