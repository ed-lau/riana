"""Offline evaluator for the Δk Monte Carlo (``run_deltak_mc.sh`` → ``bench_rising_null``).

For every simulated dataset it re-tests the SAME gated peptides under the alternative peptide
collapses and weights, with the shipped functions wherever they exist:

1. *Gated peptides.* The shipped ``riana rollup`` runs in-process with the bench's arguments and
   ``riana.core.protein._collapse_long`` wrapped to record its input — so the R² gate, Spep
   floor, peptide admission and parsimony are exactly the shipped ones.
2. *Cells.* ``per-condition`` = the shipped collapse (``_collapse_long``; checked against the
   CLI's ``riana_rollup_fractions.txt``). ``shared`` = the same collapse after each peptide's
   band is rescaled so its *level* is pooled over both conditions (``shared_bands``): the
   peptides then enter a cell in the same proportions in both conditions, so a peptide's
   systematic offset cannot masquerade as a Δk.
3. *Tests on each cell table.* The shipped ``fit_linear_deltak`` (``wls``, ``wls-var``) and the
   shipped ``fit_interaction`` (``wls``); plus ``pool2``, a research variance function
   reimplemented here (``fit_pool2``) — Var(φ) = c·Ṽ/(1 − θ̂)² + s·φ̂², c and s pooled per
   experiment — whose ``wls-var`` twin is checked against the shipped p-values.
4. *Scoring* against ``truth.tsv``: rejection rate at α 0.05 / 0.01 (overall and by true k),
   95 % CI coverage of Δk, the SD of the z-scores, mean log₂ k_B/k_A; for the interaction the
   mean-Δk, additive and ratio-scale tests' rejection rates.

Writes ``<root>/evaluation.tsv`` (one row per dataset × route × collapse × test) and
``<root>/evaluation_audit.tsv`` (the reproduction checks), and prints the arm-level table.

``--real REF TEST`` instead takes a real fit directory (no truth): calls per collapse × weights,
their overlap with the shipped default, each protein's Δk against its peptides' own Δk (from the
fit's per-peptide curves, no collapse), the size of the mixing shift against the SE, and the
spread of the per-peptide band-variance ratio between the conditions.

    python -m tests.benchmark.eval_rising_variance runs/deltak_mc
    python -m tests.benchmark.eval_rising_variance runs/atf6_individual/fit --real control ATF6-KO

Reference: ``reports/2026-11-12_deltak_inference_manual.md``.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls

import riana.core.protein as protein
from riana.core.interaction_model import fit_interaction
from riana.core.linear_model import _fit_fdist, fit_linear_deltak

PHI_LIMIT = -4.0
REF, TEST = "A", "B"                                         # the bench's conditions


def rd(path: Path) -> pd.DataFrame:
    return pd.read_table(path, comment="#", low_memory=False)


# --- 1. gated peptides, from the shipped rollup ----------------------------------------------

def capture_gated(fit_dir: Path, ref: str = REF, test: str = TEST) -> dict:
    """Run the shipped ``riana rollup`` on ``fit_dir`` and return the frame it collapses."""
    from typer.testing import CliRunner
    from riana.cli import app

    seen: dict = {}
    orig = protein._collapse_long

    def spy(fractions, *, method, min_peptides):
        seen.update(frame=fractions.copy(), method=method, min_peptides=min_peptides)
        return orig(fractions, method=method, min_peptides=min_peptides)

    protein._collapse_long = spy
    try:
        with tempfile.TemporaryDirectory() as tmp:
            res = CliRunner().invoke(app, [
                "rollup", str(fit_dir), "--model", "linear simple", "--reference-condition",
                ref, "--test-condition", test, "--linear-weights", "wls", "-o", tmp])
            if res.exit_code != 0:
                raise RuntimeError(f"rollup failed on {fit_dir}:\n{res.output[-2000:]}")
    finally:
        protein._collapse_long = orig
    return seen


# --- 2. cells ----------------------------------------------------------------------------------

def shared_bands(frame: pd.DataFrame) -> pd.DataFrame:
    """Rescale each peptide's band so its variance level is shared by the conditions.

    A peptide's weight in a cell is 1/σ², σ = band/3.29. Per (experiment, peptide, condition)
    the σ² over its points has mean m_c; over both conditions, m. Each point's band is scaled
    by √(m/m_c): the time profile within a condition is kept, the level becomes the
    peptide's pooled one. Points without a usable band are left alone (the collapse fills
    them with the cell median, as shipped).
    """
    f = frame.copy()
    s2 = ((f.fs_upper - f.fs_lower) / protein._PI_SPAN_SIGMA) ** 2
    s2 = s2.where(np.isfinite(s2) & (s2 > 0))
    m_c = s2.groupby([f.experiment, f.concat, f.condition]).transform("mean")
    m = s2.groupby([f.experiment, f.concat]).transform("mean")
    r = np.sqrt(m / m_c).fillna(1.0)
    f["fs_lower"] = f.fs - (f.fs - f.fs_lower) * r
    f["fs_upper"] = f.fs + (f.fs_upper - f.fs) * r
    return f


def collapse(seen: dict, shared: bool) -> pd.DataFrame:
    frame = shared_bands(seen["frame"]) if shared else seen["frame"]
    long = protein._collapse_long(frame, method=seen["method"],
                                  min_peptides=seen["min_peptides"])
    ind = protein._point_individuals(seen["frame"])
    long["individual"] = [ind.get((e, c, float(t), int(b)), "") for e, c, t, b in zip(
        long.experiment, long.condition, long.labeling_time, long.biological_replicate)]
    return long


# --- 3. tests ----------------------------------------------------------------------------------

def dk_shipped(long: pd.DataFrame, weights: str, ref: str = REF,
               test: str = TEST) -> pd.DataFrame:
    """The shipped Δk test, reshaped to one row per (experiment, protein)."""
    out = fit_linear_deltak(long, phi_limit=PHI_LIMIT, weights=weights,
                            reference_condition=ref, test_condition=test)
    k = out.pivot_table(index=["experiment", "protein"], columns="condition",
                        values="k_deg").rename(columns={ref: "k_A", test: "k_B"})
    n = out.groupby(["experiment", "protein"]).n_points.sum().rename("n")
    dk = out.dropna(subset=["delta_k_p"]).drop_duplicates(["experiment", "protein"]).set_index(
        ["experiment", "protein"])[["delta_k", "delta_k_se", "delta_k_p", "delta_k_p_adj"]]
    j = dk.join(k).join(n).reset_index()
    return j.rename(columns={"delta_k": "dk", "delta_k_se": "se", "delta_k_p": "p",
                             "delta_k_p_adj": "p_adj"}).assign(df=lambda d: d.n - 2)


def interaction_shipped(long: pd.DataFrame, weights: str, workers: int, ref: str = REF,
                        test: str = TEST) -> pd.DataFrame:
    out = fit_interaction(long, reference_condition=ref, test_condition=test,
                          phi_limit=PHI_LIMIT, weights=weights, workers=workers)
    return out.drop_duplicates("protein")


def _wls(t, cond, y, w):
    X = np.column_stack([t * (cond == 0), t * (cond == 1)])
    sw = np.sqrt(w)
    Xw, yw = X * sw[:, None], y * sw
    try:
        inv = np.linalg.inv(Xw.T @ Xw)
    except np.linalg.LinAlgError:
        return None
    b = inv @ (Xw.T @ yw)
    e = y - X @ b
    n = len(y)
    if n <= 2:
        return None
    h = np.einsum("ij,jk,ik->i", Xw, inv, Xw)
    return b, inv * float(np.sum(w * e ** 2) / (n - 2)), e, h, X


def _prepare(long: pd.DataFrame) -> pd.DataFrame:
    """φ, the eBayes-moderated cell variance Ṽ (as the shipped wls-var) and the t > 0 /
    plateau filters, for the reimplemented fits."""
    c = long.copy()
    ok = (np.isfinite(c.theta_var) & (c.theta_var > 0) & np.isfinite(c.theta_df)
          & (c.theta_df > 0))
    d0, s0 = _fit_fdist(c.theta_var[ok].to_numpy(), c.theta_df[ok].to_numpy())
    vm = (d0 * s0 + c.theta_df * c.theta_var) / (d0 + c.theta_df)
    c["vmod"] = np.where(np.isfinite(vm) & (vm > 0), vm, s0)
    c["phi"] = np.log1p(-np.clip(c.theta.to_numpy(float), 0.001, 0.999))
    return c[(c.labeling_time > 0) & (c.phi > PHI_LIMIT)].copy()


def _fit_all(c: pd.DataFrame, variant: str, params: dict | None = None, test: str = TEST):
    rows, pieces = [], []
    for (exp, prot), g in c.groupby(["experiment", "protein"], sort=False):
        cond = (g.condition.astype(str) == test).to_numpy().astype(int)
        if min((cond == 0).sum(), (cond == 1).sum()) < 2 or len(g) < 3:
            continue
        t, y = g.labeling_time.to_numpy(float), g.phi.to_numpy(float)
        v = g.vmod.to_numpy(float)
        ols = _wls(t, cond, y, np.ones_like(y))
        if ols is None:
            continue
        phat = np.maximum(ols[4] @ ols[0], PHI_LIMIT)
        A, B = v / np.exp(2 * phat), phat ** 2          # measurement / animal-turnover bases
        if variant == "wlsvar":
            w = np.exp(2 * phat) / v
        elif variant == "pool2":
            w = 1.0 / (params["c"][exp] * A + params["s"][exp] * B)
        else:
            raise ValueError(variant)
        res = _wls(t, cond, y, w)
        if res is None:
            continue
        b, V, e, h, X = res
        se = float(np.sqrt(V[0, 0] + V[1, 1] - 2 * V[0, 1]))
        dk = -(b[1] - b[0])
        rows.append(dict(experiment=exp, protein=prot, k_A=-b[0], k_B=-b[1], dk=dk, se=se,
                         df=len(y) - 2,
                         p=float(2 * stats.t.sf(abs(dk / se), len(y) - 2)) if se > 0 else np.nan))
        pieces.append(pd.DataFrame(dict(experiment=exp, e2=e ** 2 / np.maximum(1 - h, 0.05),
                                        A=A, B=B)))
    return pd.DataFrame(rows), (pd.concat(pieces) if pieces else pd.DataFrame())


def fit_pool2(long: pd.DataFrame, iters: int = 4, test: str = TEST) -> tuple[pd.DataFrame, dict]:
    """Research variance function Var(φ) = c·Ṽ/(1 − θ̂)² + s·φ̂², (c, s) per experiment from the
    leverage-corrected squared residuals by weighted NNLS, iterated with the refit (starts
    from wls-var). The s·φ̂² term is variance growing as (k·t)², the signature of
    between-animal turnover variation."""
    c = _prepare(long)
    params = {"c": {e: 1.0 for e in c.experiment.unique()},
              "s": {e: 0.0 for e in c.experiment.unique()}}
    variant = "wlsvar"
    for _ in range(iters):
        _, P = _fit_all(c, variant, params, test)
        for exp, g in P.groupby("experiment"):
            M = np.column_stack([g.A.to_numpy(), g.B.to_numpy()])
            wt = 1.0 / np.maximum(M @ np.array([params["c"][exp], params["s"][exp]]), 1e-12)
            coef, _ = nnls(M * wt[:, None], g.e2.to_numpy() * wt)
            params["c"][exp], params["s"][exp] = max(float(coef[0]), 1e-12), float(coef[1])
        variant = "pool2"
    fits, P = _fit_all(c, "pool2", params, test)
    # How much of the modeled variance is the time-growing (animal) term, per experiment.
    params["share"] = {exp: float(np.mean(params["s"][exp] * g.B / (params["c"][exp] * g.A
                                                                     + params["s"][exp] * g.B)))
                       for exp, g in P.groupby("experiment")}
    return fits, params


# --- 4. scoring --------------------------------------------------------------------------------

K_BINS = [0, 0.01, 0.03, 0.08, np.inf]
K_LABELS = ["k<.01", "k.01-.03", "k.03-.08", "k>.08"]


def score_dk(fits: pd.DataFrame, truth: pd.DataFrame, ratio: float) -> dict:
    tw = truth.pivot_table(index=["experiment", "protein"], columns="condition",
                           values="k_true").rename(columns={REF: "tA", TEST: "tB"})
    j = fits.merge(tw.reset_index(), on=["experiment", "protein"])
    dk_true = j.tB - j.tA
    q = stats.t.ppf(0.975, j.df)
    pos = (j.k_A > 0) & (j.k_B > 0)
    bins = pd.cut(j.tA, K_BINS, labels=K_LABELS)
    by_k = {f"rej05_{b}": float((g.p < .05).mean()) for b, g in j.groupby(bins, observed=True)}
    return dict(tests=len(j), reject05=float((j.p < .05).mean()),
                reject01=float((j.p < .01).mean()), **by_k,
                dk_cover95=float(((j.dk - q * j.se <= dk_true) & (dk_true <= j.dk + q * j.se)).mean()),
                sd_z=float(((j.dk - dk_true) / j.se).std()),
                dk_rel_bias=float(np.median((j.dk - dk_true) / dk_true)) if ratio != 1 else np.nan,
                log2_ratio=float(np.log2(j.k_B[pos] / j.k_A[pos]).mean()),
                log2_ratio_true=float(np.log2(ratio)))


def score_interaction(inter: pd.DataFrame) -> dict:
    rate = lambda p: float((p.dropna() < .05).mean()) if p.notna().any() else np.nan
    return dict(tests=len(inter), mean_rej05=rate(inter.mean_delta_k_p),
                additive_rej05=rate(inter.interaction_additive_p),
                ratio_rej05=rate(inter.interaction_p),
                sigma_animal_med=float(inter.sigma_animal.median()))


# --- real data (no truth: how far the calls move) ----------------------------------------------

def _bh_per_experiment(fits: pd.DataFrame) -> pd.Series:
    from statsmodels.stats.multitest import multipletests
    out = pd.Series(np.nan, index=fits.index)
    for _, g in fits.dropna(subset=["p"]).groupby("experiment"):
        out[g.index] = multipletests(g.p.to_numpy(), method="fdr_bh")[1]
    return out


def peptide_dk(fit_dir: Path, frame: pd.DataFrame, ref: str, test: str) -> pd.DataFrame:
    """Each admitted peptide's own Δk (its k in ``test`` minus its k in ``ref``, from the fit's
    per-peptide curves), summarised per (experiment, protein): the median, and the protein's
    within-protein peptide k CV in ``ref`` — independent of any collapse."""
    pep = rd(fit_dir / "riana_fit_peptides.txt")
    adm = frame[["experiment", "protein", "concat"]].drop_duplicates()
    kk = adm.merge(pep.pivot_table(index=["experiment", "concat"], columns="condition",
                                   values="k_deg").reset_index(), on=["experiment", "concat"])
    kk["dk_pep"] = kk[test] - kk[ref]
    return kk.groupby(["experiment", "protein"]).agg(
        dk_pep_med=("dk_pep", "median"),
        pep_k_cv=(ref, lambda s: s.std() / s.mean())).reset_index()


def evaluate_real(fit_dir: Path, ref: str, test: str, workers: int) -> pd.DataFrame:
    """Calls (BH < 0.05, per experiment as shipped) under each collapse × weights on a real
    fit; their overlap with the shipped default (per-condition collapse, ``wls``); and how
    close each protein Δk sits to its peptides' own Δk."""
    seen = capture_gated(fit_dir, ref, test)
    pdk = peptide_dk(fit_dir, seen["frame"], ref, test)
    print(f"within-protein peptide k CV ({ref}): median {pdk.pep_k_cv.median():.3f}, "
          f"90th percentile {pdk.pep_k_cv.quantile(.9):.3f}")
    rows, calls = [], {}
    for shared in (False, True):
        coll = "shared" if shared else "per-condition"
        long = collapse(seen, shared)
        for w in ("wls", "wls-var", "pool2"):
            if w == "pool2":
                fits, params = fit_pool2(long, test=test)
                fits["p_adj"] = _bh_per_experiment(fits)
                extra = {f"s_share_{e}": v for e, v in params["share"].items()}
            else:
                fits, extra = dk_shipped(long, w, ref, test), {}
            hit = set(map(tuple, fits.loc[fits.p_adj < .05, ["experiment", "protein"]].to_numpy()))
            calls[(coll, w)] = hit
            m = fits.merge(pdk, on=["experiment", "protein"]).dropna(subset=["dk", "dk_pep_med"])
            h = m[m.p_adj < .05]
            extra.update(r_vs_peptide=float(np.corrcoef(m.dk, m.dk_pep_med)[0, 1]),
                         gap_vs_peptide=float((m.dk - m.dk_pep_med).abs().median()),
                         calls_abs_dk=float(h.dk.abs().median()) if len(h) else np.nan,
                         calls_abs_peptide_dk=float(h.dk_pep_med.abs().median()) if len(h) else np.nan,
                         p01=int((fits.p < .01).sum()))
            rows.append(dict(collapse=coll, test="dk", weights=w, tests=len(fits),
                             calls=len(hit), **extra))
        if long.experiment.nunique() < 2:              # one experiment: no interaction
            continue
        inter = interaction_shipped(long, "wls", workers, ref, test)
        for col, name in (("mean_delta_k_p_adj", "interaction: mean Δk"),
                          ("interaction_p_adj", "interaction: ratio scale"),
                          ("interaction_additive_p_adj", "interaction: additive")):
            hit = set(inter.loc[inter[col] < .05, "protein"])
            calls[(coll, name)] = hit
            rows.append(dict(collapse=coll, test=name, weights="wls", tests=len(inter),
                             calls=len(hit), sigma_animal_med=float(inter.sigma_animal.median())))
            if shared and name == "interaction: ratio scale":
                print(f"ratio-scale interaction hits, shared collapse: {sorted(hit)}")
    # Size of the mixing effect against the SE and the remaining Δk, and the chance spread of
    # the per-peptide band-variance ratio between conditions that drives it.
    pc = dk_shipped(collapse(seen, False), "wls", ref, test)
    sh = dk_shipped(collapse(seen, True), "wls", ref, test)
    m = pc.merge(sh, on=["experiment", "protein"], suffixes=("", "_s"))
    shift = (m.dk - m.dk_s).abs()
    f = seen["frame"]
    s2 = ((f.fs_upper - f.fs_lower) / protein._PI_SPAN_SIGMA) ** 2
    lvl = s2.groupby([f.experiment, f.concat, f.condition]).mean().unstack()
    ratio = (lvl[test] / lvl[ref]).dropna()
    npts = f.groupby(["experiment", "concat", "condition"]).size().unstack().median()
    print(f"peptides per protein (median) {f.groupby(['experiment', 'protein']).concat.nunique().median():.0f}; "
          f"points per peptide per condition {npts.round(0).to_dict()}")
    print(f"Δk SE median {m.se.median():.4f} | |mixing shift| median {shift.median():.4f}, "
          f"90th pct {shift.quantile(.9):.4f} | |shift|/SE median {(shift / m.se).median():.2f}, "
          f"share > 2 SE {(shift > 2 * m.se).mean():.3f} | |Δk shared|/SE median "
          f"{(m.dk_s.abs() / m.se_s).median():.2f}")
    print(f"peptide band-variance ratio {test}/{ref}: median {ratio.median():.2f}, IQR "
          f"{ratio.quantile(.25):.2f}–{ratio.quantile(.75):.2f}, 5–95 % "
          f"{ratio.quantile(.05):.2f}–{ratio.quantile(.95):.2f}")
    R = pd.DataFrame(rows)
    base = {"dk": calls[("per-condition", "wls")]}
    R["shared_with_default"] = [
        len(calls[(r.collapse, r.weights if r.test == "dk" else r.test)]
            & (base["dk"] if r.test == "dk" else calls[("per-condition", r.test)]))
        for r in R.itertuples()]
    return R


# --- driver ------------------------------------------------------------------------------------

def evaluate_dataset(d: Path, arm: str, kp: float, workers: int) -> tuple[list, list]:
    ratio = float(d.name[len("ratio"):].rsplit("_", 1)[0])
    truth = rd(d / "truth.tsv")
    rows, audit = [], []
    base = dict(arm=arm, dataset=d.name, ratio=ratio, kp=kp)
    routes = ("theta", "f") if kp < 50 else ("theta",)      # instant arms: f ≡ θ
    for route in routes:
        seen = capture_gated(d / route)
        for shared in (False, True):
            coll = "shared" if shared else "per-condition"
            long = collapse(seen, shared)
            if not shared:                                    # audit: shipped cells / p-values
                cli = rd(d / route / "rollup_wls" / "riana_rollup_fractions.txt")
                m = long.merge(cli, on=["experiment", "condition", "protein", "labeling_time",
                                        "biological_replicate"])
                audit.append(dict(**base, route=route, check="cells vs CLI fractions",
                                  n=len(m), cli_rows=len(cli),
                                  max_abs=float(max((m.theta - m.fs).abs().max(),
                                                    (m.theta_var - m.fs_var).abs().max()))))
            for w in ("wls", "wls-var"):
                fits = dk_shipped(long, w)
                rows.append(dict(**base, route=route, collapse=coll, test="dk", weights=w,
                                 **score_dk(fits, truth, ratio)))
                if not shared:
                    cli = rd(d / route / f"rollup_{w}" / "riana_rollup_proteins.txt").dropna(
                        subset=["delta_k_p"]).drop_duplicates(["experiment", "protein"])
                    m = fits.merge(cli, on=["experiment", "protein"])
                    audit.append(dict(**base, route=route, check=f"dk p vs CLI ({w})",
                                      n=len(m), cli_rows=len(cli),
                                      max_abs=float((m.p - m.delta_k_p).abs().max())))
                if w == "wls-var":
                    mine, _ = _fit_all(_prepare(long), "wlsvar")
                    m = mine.merge(fits[["experiment", "protein", "p"]],
                                   on=["experiment", "protein"], suffixes=("", "_ship"))
                    audit.append(dict(**base, route=route,
                                      check=f"reimplemented wls-var p vs shipped ({coll})",
                                      n=len(m), cli_rows=len(fits),
                                      max_abs=float((m.p - m.p_ship).abs().max())))
            fits, params = fit_pool2(long)
            rows.append(dict(**base, route=route, collapse=coll, test="dk", weights="pool2",
                             pool2_s=float(np.mean(list(params["s"].values()))),
                             **score_dk(fits, truth, ratio)))
            inter = interaction_shipped(long, "wls", workers)
            rows.append(dict(**base, route=route, collapse=coll, test="interaction",
                             weights="wls", **score_interaction(inter)))
        print(f"  {arm}/{d.name}/{route} done", flush=True)
    return rows, audit


def summarise(R: pd.DataFrame) -> pd.DataFrame:
    keys = ["arm", "ratio", "route", "test", "collapse", "weights"]
    num = [c for c in R.columns if c not in keys + ["dataset", "kp", "tests"]
           and pd.api.types.is_numeric_dtype(R[c])]

    def agg(g: pd.DataFrame) -> pd.Series:
        out = {c: (np.average(g[c].dropna(), weights=g.tests[g[c].notna()])
                   if g[c].notna().any() else np.nan) for c in num}
        out["tests"] = g.tests.sum()
        out["datasets"] = len(g)
        out["reject05_range"] = (f"{g.reject05.min():.3f}–{g.reject05.max():.3f}"
                                 if "reject05" in g and g.reject05.notna().any() else "")
        return pd.Series(out)

    return R.groupby(keys, sort=False).apply(agg, include_groups=False)


def main(root: str, workers: int = 16) -> None:
    rootp = Path(root)
    rows, audit = [], []
    for arm_dir in sorted(p for p in rootp.iterdir() if p.is_dir()):
        datasets = sorted(arm_dir.glob("ratio*_*"))
        datasets = [d for d in datasets if d.is_dir()]
        for d in datasets:
            pfile = arm_dir / f"{d.name.rsplit('_', 1)[0]}_params.json"
            kp = float(json.loads(pfile.read_text())["kp"]) if pfile.exists() else np.nan
            r, a = evaluate_dataset(d, arm_dir.name, kp, workers)
            rows += r
            audit += a
    R, A = pd.DataFrame(rows), pd.DataFrame(audit)
    R.to_csv(rootp / "evaluation.tsv", sep="\t", index=False)
    A.to_csv(rootp / "evaluation_audit.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print("\naudit (max |Δ| between this evaluator and the shipped outputs):")
    print(A.groupby("check").agg(n=("n", "sum"), max_abs=("max_abs", "max")).to_string())
    print()
    print(summarise(R).round(3).to_string())


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root", nargs="?", default="runs/deltak_mc",
                    help="Monte Carlo root (run_deltak_mc.sh's OUT), or a fit directory with --real")
    ap.add_argument("-W", "--workers", type=int, default=16)
    ap.add_argument("--real", nargs=2, metavar=("REFERENCE", "TEST"),
                    help="real-data mode: ROOT is a fit directory; compare these conditions")
    a = ap.parse_args()
    if a.real:
        pd.set_option("display.width", 250)
        R = evaluate_real(Path(a.root), *a.real, a.workers)
        print(R.round(3).to_string())
    else:
        main(a.root, a.workers)
