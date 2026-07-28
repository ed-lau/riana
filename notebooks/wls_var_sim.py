"""Peptide-level MC for the wls-var scenarios: single-peptide proteins + collapse-vs-pooled.

Key modeling: each peptide has an intrinsic k-offset SHARED across conditions (a peptide is
measured in both control and treatment), so the offset cancels in Δk — the realistic paired design.
"""
import numpy as np, pandas as pd, collections
from scipy.special import digamma, polygamma
from scipy.optimize import brentq
from tests.benchmark.bench_linear_weights import DAYS, PHI_LIMIT, to_phi, fit, ebayes_var

_T = DAYS[DAYS > 0]


def simulate_protein(kA, kB, n_pep, rng, *, n_biorep=3, tau=0.3, sig=0.05, sig_spread=0.5, var_df=6):
    """One protein, 2 conditions, PAIRED peptides. offset_j ~ N(0,tau) shared across conditions
    (peptide intrinsic k-bias → within-condition between-peptide scatter = the dispersion, but it
    cancels in Δk). θ = 1-exp(-k_cond·e^offset·t)+N(0,σ_pep); σ_pep heteroscedastic; σ̂ noisy at
    df=var_df. Rows: cond, biorep, t, pep, θ, σ̂, df."""
    offset = rng.normal(0, tau, n_pep)                               # SHARED across conditions
    sigp = sig * np.exp(rng.normal(0, sig_spread, n_pep) - sig_spread ** 2 / 2)
    rows = []
    for cond, kc in ((0, kA), (1, kB)):
        kpep = kc * np.exp(offset)
        for j in range(n_pep):
            for br in range(n_biorep):
                for t in _T:
                    theta = 1 - np.exp(-kpep[j] * t) + rng.normal(0, sigp[j])
                    sighat = sigp[j] * np.sqrt(rng.chisquare(var_df) / var_df)
                    rows.append((cond, br, t, j, theta, sighat, var_df))
    return np.array(rows, float)


def collapse(arr):
    """Inverse-variance collapse per (cond, biorep, t) cell → 1 point + Var + Satterthwaite df."""
    cells = collections.defaultdict(list)
    for cond, br, t, pep, th, sg, df in arr:
        cells[(cond, br, t)].append((th, sg, df))
    out = []
    for (cond, br, t), items in cells.items():
        th = np.array([x[0] for x in items]); sg = np.array([x[1] for x in items]); d = np.array([x[2] for x in items])
        w = 1 / sg ** 2
        out.append((cond, t, float(np.sum(w * th) / np.sum(w)), float(1 / np.sum(w)),
                    float(np.sum(w) ** 2 / np.sum(w ** 2 / np.maximum(d, 1)))))
    return np.array(out, float)                                      # cond, t, theta, var, df


def pooled(arr):
    """Every (cond, biorep, t, peptide) is its own point; var = σ̂² (per-peptide, unbiased)."""
    return np.column_stack([arr[:, 0], arr[:, 2], arr[:, 4], arr[:, 5] ** 2, arr[:, 6]])


def fit_pts(pts, how, d0=2.0):
    cond, t, theta, var, df = pts[:, 0], pts[:, 1], pts[:, 2], pts[:, 3], pts[:, 4]
    phi = to_phi(theta)
    v = var
    if how == "ebayes":
        v = ebayes_var(var, df, float(np.nanmean(var)), d0); how = "wls_fit_var"
    elif how == "wls_fit":
        v = None
    return fit(t, cond, phi, theta, how, PHI_LIMIT, var=v)


def run(kA, kB, n_pep, nsim, seed, d0=2.0, methods=("wls_fit", "wls_fit_var", "ebayes")):
    """MC over `nsim` proteins. `reject` = the Δk-test rejection rate — Type-I when kA==kB,
    power otherwise. (k_A absolute coverage is intentionally NOT reported: it is contaminated
    by the peptide-sampling bias that cancels in the paired Δk, so it is not the inferential
    target.) Returns a per-(method, rollup) summary."""
    rng = np.random.default_rng(seed)
    acc = {(m, r): {"dk": [], "rej": []} for m in methods for r in ("collapse", "pooled")}
    for _ in range(nsim):
        arr = simulate_protein(kA, kB, n_pep, rng)
        for rname, pts in (("collapse", collapse(arr)), ("pooled", pooled(arr))):
            for m in methods:
                r = fit_pts(pts, m, d0)
                if r is None:
                    continue
                acc[(m, rname)]["dk"].append(r["dk"])
                acc[(m, rname)]["rej"].append(r["dk_p"] < 0.05)
    rows = []
    for (m, rname), a in acc.items():
        dk = np.array(a["dk"])
        rows.append(dict(method=m, rollup=rname, n=len(dk),
                         dk_bias=np.mean(dk) - (kB - kA),
                         dk_rmse=np.sqrt(np.mean((dk - (kB - kA)) ** 2)),
                         reject=np.mean(a["rej"])))
    return pd.DataFrame(rows)


def d0_grid(kA, kB, n_pep, nsim, seed, d0s=(0.5, 1, 2, 4, 8), rollup="collapse"):
    """Sweep the eBayes prior df d0 for one rollup method — where does Type-I land vs gain?
    Returns a table of (d0, Type-I-or-power, dk_rmse) for the `ebayes` arm."""
    rows = []
    for d0 in d0s:
        r = run(kA, kB, n_pep, nsim, seed, d0=d0, methods=("ebayes",))
        rr = r[r.rollup == rollup].iloc[0]
        rows.append(dict(d0=d0, reject=rr["reject"], dk_rmse=rr["dk_rmse"]))
    return pd.DataFrame(rows)


def fit_fdist(s2, df, trim=0.0):
    """Smyth (2004) fitFDist → (d0, s0²) from sample variances s2 and their df. `trim` > 0
    Winsorizes the log-variance deviations before the moment-match — a limma `robust=TRUE`
    (Phipson 2016) flavour, needed because proteomics variances are heavy-tailed and the
    outliers bias d0 LOW (under-shrinkage)."""
    s2, df = np.asarray(s2, float), np.asarray(df, float)
    ok = np.isfinite(s2) & (s2 > 0) & (df > 0); s2, df = s2[ok], df[ok]
    e = np.log(s2) - digamma(df / 2) + np.log(df / 2)
    if trim > 0:
        lo, hi = np.quantile(e, [trim, 1 - trim]); e = np.clip(e, lo, hi)
    evar = np.var(e, ddof=1) - np.mean(polygamma(1, df / 2))
    if evar <= 0:
        return np.inf, float(np.exp(np.mean(e)))
    try:
        half = brentq(lambda x: polygamma(1, x) - evar, 1e-4, 1e5); d0 = 2 * half
    except Exception:
        d0 = np.inf
    s0 = float(np.exp(np.mean(e) + digamma(d0 / 2) - np.log(d0 / 2))) if np.isfinite(d0) else float(np.exp(np.mean(e)))
    return d0, s0


# --------------------------------------------------------------------------- #
# Real-data validation — biorep-split stability of per-protein k (no ground truth)
# --------------------------------------------------------------------------- #
def collapse_cells(ff, pi_span=3.29):
    """Peptide `riana_fit_fractions` → per (protein, condition, biorep, t) cells with θ,
    Var(θ)=1/Σ(1/σ²), and a Satterthwaite df (per-peptide σ-df ≈ its fit-point count)."""
    ff = ff[ff.labeling_time > 0].copy()
    ff["sig"] = (ff.fs_upper - ff.fs_lower) / pi_span
    ff = ff[np.isfinite(ff.sig) & (ff.sig > 0) & np.isfinite(ff.fs)].copy()
    ff["npts"] = ff.groupby(["protein id", "concat"])["labeling_time"].transform("nunique")
    ff["w"] = 1 / ff.sig ** 2; ff["wth"] = ff.w * ff.fs; ff["w2d"] = ff.w ** 2 / np.maximum(ff.npts, 1)
    g = ff.groupby(["protein id", "condition", "biological_replicate", "labeling_time"]).agg(
        sw=("w", "sum"), swth=("wth", "sum"), sw2d=("w2d", "sum"), npep=("concat", "nunique")).reset_index()
    g["theta"] = g.swth / g.sw; g["var"] = 1 / g.sw; g["df"] = g.sw ** 2 / g.sw2d
    g["phi"] = to_phi(g.theta.to_numpy())
    return g


def _fit_k(t, phi, var=None, phi_lim=-4.0):
    """Single-condition k via one IRLS step: φ = −k·t through origin, weight (1−θ̂)²[/var]."""
    keep = np.isfinite(phi) & (phi > phi_lim) & (t > 0); t, phi = t[keep], phi[keep]
    if len(t) < 3:
        return np.nan
    b0 = np.sum(t * phi) / np.sum(t * t)                     # OLS pilot slope
    w = np.exp(2 * np.maximum(b0 * t, phi_lim))
    if var is not None:
        w = w / np.clip(var[keep], 1e-12, None)
    return -np.sum(w * t * phi) / np.sum(w * t * t)


def biorep_split(run_dir, d0_trim=0.1):
    """Biorep-split stability of per-protein k: fit k per (protein, condition, biorep) under
    `wls` vs `wls-var` (eBayes with a robust-fitFDist d0), and measure |log(k_b1/k_b2)| within
    each (protein, condition). Lower = more consistent across biological replicates = better.
    Returns (summary_df, d0) or (None, nan) if the run has < 2 bioreps."""
    ff = pd.read_table(f"{run_dir}/riana_fit_fractions.txt", comment="#")
    c = collapse_cells(ff)
    brs = sorted(c.biological_replicate.unique())
    if len(brs) < 2:
        return None, float("nan")
    d0, _ = fit_fdist(c["var"].to_numpy(), c["df"].to_numpy(), trim=d0_trim)
    pool = float(np.nanmean(c["var"]))
    rows = []
    for (prot, cond, br), g in c.groupby(["protein id", "condition", "biological_replicate"]):
        t, phi, var, df = g.labeling_time.to_numpy(), g.phi.to_numpy(), g["var"].to_numpy(), g["df"].to_numpy()
        rows.append(dict(prot=prot, cond=cond, br=br, maxpep=int(g.npep.max()),
                         k_wls=_fit_k(t, phi), k_var=_fit_k(t, phi, ebayes_var(var, df, pool, d0))))
    K = pd.DataFrame(rows)
    mp = K.groupby(["prot", "cond"])["maxpep"].max()
    out = {}
    for col, lab in [("k_wls", "wls"), ("k_var", "wls-var")]:
        p = K.pivot_table(index=["prot", "cond"], columns="br", values=col)[brs[:2]].dropna()
        p = p[(p[brs[0]] > 0) & (p[brs[1]] > 0)]
        d = np.abs(np.log(p[brs[0]] / p[brs[1]])); m = mp.reindex(d.index)
        out[lab] = dict(all=d.median(), multi_pep=d[m >= 2].median(), single_pep=d[m < 2].median(), n=len(d))
    return pd.DataFrame(out).T, d0


if __name__ == "__main__":
    pd.set_option("display.width", 120)
    fmt = lambda x: f"{x:8.3f}"
    print("SCENARIO 2 — collapse vs pooled, multi-peptide (n_pep=4).  TYPE-I: kA=kB=0.1\n")
    print(run(0.1, 0.1, 4, 500, 1).to_string(index=False, float_format=fmt))
    print("\n  POWER: kA=0.1, kB=0.14\n")
    print(run(0.1, 0.14, 4, 500, 3).to_string(index=False, float_format=fmt))
    print("\nSCENARIO 1 — single-peptide proteins (n_pep=1).  TYPE-I: kA=kB=0.1\n")
    print(run(0.1, 0.1, 1, 500, 2).to_string(index=False, float_format=fmt))
    rng = np.random.default_rng(9); s2, dfs = [], []
    for _ in range(400):
        c = collapse(simulate_protein(0.1, 0.1, 4, rng)); s2 += c[:, 3].tolist(); dfs += c[:, 4].tolist()
    d0, s0 = fit_fdist(np.array(s2), np.array(dfs))
    print(f"\nfitFDist on collapsed-cell variances → d0 = {d0:.2f}  (bench frontier liked d0≈0.5–2)")
