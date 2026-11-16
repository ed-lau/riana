"""WP4 prototypes: D (random-effects collapse) and E (peptide-level mixed model).

``wls-var`` was removed from Riana on 2026-11-15: D+var need Riana at commit
``e12980c`` or earlier (the wls-var / D / E evaluation).

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md`` WP4. Evaluator code only, not a
CLI option. Each variant maps the gated frame the shipped rollup collapses (``seen``, from
``eval_rising_variance.capture_gated``) to one Δk test per (experiment, protein), with the
columns of ``eval_rising_variance.dk_shipped``: dk, se, p, p_adj (BH per experiment), k_A,
k_B, df.

**D, random-effects collapse.** C's shared band levels (``shared_bands``), then each
peptide's weight in a cell becomes 1/(σᵢ² + τ²(θ̄)), with τ(θ̄) = τ_rel·g(θ̄) and
g(θ) = (1 − θ)(−ln(1 − θ)), the apparent-k offset form WP1b found. τ_rel is estimated per
experiment, pooled over both conditions, by moments (DerSimonian–Laird with the g scaling):
τ_rel² = Σ(Q − (n − 1)) / Σ g²(Σw − Σw²/Σw) over cells with ≥ 2 peptides, from the shared
bands. The shipped collapse then runs on the inflated bands, so its cell variance is
Ṽ_RE = 1/Σ(1/(σᵢ² + τ²)).

- ``D``: the shipped ``wls`` test on D's cells;
- ``D+var``: the shipped ``wls-var`` test on D's cells (Ṽ_RE in place of Ṽ, moderated as
  shipped).

**E, peptide-level mixed model**, per (experiment, protein), on the admitted peptides'
points (t > 0, φ above the plateau limit)::

    φ = −k_c·t − b_p·t + u_animal·t + ε

- k_c: one slope per condition (Δk = k_test − k_ref); b_p: fixed peptide slope offsets with
  sum-to-zero coding, shared by the conditions (an apparent-k offset per peptide, WP1b), so
  k_c is the unweighted peptide-average k and a peptide's offset cannot enter Δk;
- u_animal: a random slope per animal (``individual``; else the condition × time × replicate
  composite, which is the cell). When every animal sits in one condition (genotype), it is
  one slope. When animals carry both conditions (tissues of one animal, e.g. T6), it is a
  slope per condition with a free 2 × 2 covariance, so the shared part cancels from Δk (the
  paired model);
- Var(ε) ∝ 1/w, from a first OLS pass as the shipped ``wls``: ``E``: w = (1 − θ̂)²;
  ``E+band``: w = (1 − θ̂)²/σᵢ², σᵢ from the shared band level (C).

REML, statsmodels ``MixedLM`` on the √w-scaled data (as ``riana.core.interaction_model``),
BFGS then Nelder–Mead; when the fit fails or the animal variance is at 0, the WLS fit with
no animal term. df: animals − 2 (between animals) or animals − 1 (paired). φ = ln(1 − θ)
with θ capped at 0.999 but not floored: a peptide θ below 0 gives φ > 0 (the collapse's
floor at 0.001 would bias single noisy points upward).
"""
from __future__ import annotations

import os
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

import riana.core.protein as protein
from tests.benchmark.eval_rising_variance import (REF, TEST, _bh_per_experiment, dk_shipped,
                                                  shared_bands)

PHI_LIMIT = -4.0
CELL = ["experiment", "condition", "protein", "labeling_time", "biological_replicate"]


def _sigma(f: pd.DataFrame) -> pd.Series:
    """Band σ per point, missing ones filled with the cell median (as ``_weighted_theta``)."""
    s = (f.fs_upper - f.fs_lower) / protein._PI_SPAN_SIGMA
    s = s.where(np.isfinite(s) & (s > 0))
    return s.fillna(s.groupby([f[c] for c in CELL]).transform("median"))


def _g(theta: np.ndarray) -> np.ndarray:
    th = np.clip(theta, 0.0, 0.999)
    return (1 - th) * -np.log1p(-th)


# --- D -----------------------------------------------------------------------------------------

def tau_rel(f: pd.DataFrame) -> dict[str, float]:
    """τ_rel per experiment, from the bands already in ``f`` (pass C's shared bands)."""
    f = f[np.isfinite(f.fs)].copy()
    f["sig"] = _sigma(f)
    f = f[np.isfinite(f.sig)]
    f["w"] = 1 / f.sig ** 2
    f["wx"] = f.w * f.fs
    f["w2"] = f.w ** 2
    c = f.groupby(CELL).agg(n=("w", "size"), sw=("w", "sum"), swx=("wx", "sum"),
                            sw2=("w2", "sum"))
    c["theta"] = c.swx / c.sw
    f = f.join(c.theta.rename("cth"), on=CELL)
    c["Q"] = (f.w * (f.fs - f.cth) ** 2).groupby([f[k] for k in CELL]).sum()
    c = c[(c.n >= 2) & (c.theta > 0)].reset_index()
    c["den"] = _g(c.theta.to_numpy()) ** 2 * (c.sw - c.sw2 / c.sw)
    out = {}
    for exp, d in c.groupby("experiment"):
        out[exp] = float(np.sqrt(max(0.0, (d.Q - (d.n - 1)).sum() / d.den.sum())))
    return out


def collapse_D(seen: dict) -> tuple[pd.DataFrame, dict]:
    f = shared_bands(seen["frame"])
    tr = tau_rel(f)
    sig = _sigma(f)
    w = 1 / sig ** 2
    key = [f[c] for c in CELL]
    cth = (w * f.fs).groupby(key).transform("sum") / w.groupby(key).transform("sum")
    tau = f.experiment.map(tr).fillna(0.0).to_numpy() * _g(cth.to_numpy())
    s2 = np.sqrt(sig ** 2 + tau ** 2)
    half = s2 * protein._PI_SPAN_SIGMA / 2
    f = f.assign(fs_lower=f.fs - half, fs_upper=f.fs + half)
    long = protein._collapse_long(f, method=seen["method"], min_peptides=seen["min_peptides"])
    ind = protein._point_individuals(seen["frame"])
    long["individual"] = [ind.get((e, c, float(t), int(b)), "") for e, c, t, b in zip(
        long.experiment, long.condition, long.labeling_time, long.biological_replicate)]
    return long, tr


def variant_D(seen: dict, ref: str = REF, test: str = TEST, weights: str = "wls"):
    long, _ = collapse_D(seen)
    return dk_shipped(long, weights, ref, test)


# --- E -----------------------------------------------------------------------------------------

def _animal_keys(g: pd.DataFrame) -> np.ndarray:
    comp = ("c:" + g.condition.astype(str) + "|" + g.labeling_time.astype(str) + "|"
            + g.biological_replicate.astype(str))
    if "individual" in g:
        named = g.individual.map(protein.individual_label).astype(str)   # empty g: keep str
        return np.where(named != "", "i:" + named, comp)
    return comp.to_numpy()


def fit_E_protein(g: pd.DataFrame, ref: str, test: str, band: bool,
                  mult: bool = False) -> dict | None:
    """Peptide-level mixed model for one protein (see the module docstring). ``mult``: the
    peptide offsets scale with the condition's k (slope k_c·(1 + m_p) instead of k_c + b_p);
    the offset columns are rescaled by k̂_c / k̄ over a few WLS passes, then the same
    MixedLM is fit with that scaling held fixed."""
    th = np.minimum(g.fs.to_numpy(float), 0.999)
    t = g.labeling_time.to_numpy(float)
    phi = np.log1p(-np.maximum(th, -0.5))
    keep = np.isfinite(phi) & (t > 0) & (phi > PHI_LIMIT)
    g, phi, t = g[keep], phi[keep], t[keep]
    cond = g.condition.astype(str).to_numpy()
    isA, isB = (cond == ref), (cond == test)
    animal = _animal_keys(g)
    nA, nB = len(np.unique(animal[isA])), len(np.unique(animal[isB]))
    if nA < 2 or nB < 2:
        return None
    pep = pd.factorize(g.concat)[0]
    P = int(pep.max()) + 1
    C = np.vstack([np.eye(P - 1), -np.ones(P - 1)]) if P > 1 else np.zeros((1, 0))
    X = np.column_stack([t * isA, t * isB] + [t * C[pep, j] for j in range(P - 1)])
    y = phi
    n, p = X.shape
    if n <= p + 2:
        return None
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    w = np.exp(2 * np.maximum(X @ b0, PHI_LIMIT))
    if band:
        w = w / g.sig.to_numpy(float) ** 2
    w = w / w.mean()
    sw = np.sqrt(w)
    if mult and P > 1:
        Xa = X[:, 2:]
        for _ in range(4):
            bb = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
            kc = -bb[:2]
            kbar = kc.mean()
            if not (kc.min() > 0 and kbar > 0):
                break
            scale = np.where(isA, kc[0], kc[1]) / kbar
            X = np.column_stack([X[:, :2], Xa * scale[:, None]])
    Xt, yt = X * sw[:, None], y * sw
    codes, uniq = pd.factorize(pd.Series(animal))
    paired = bool((pd.Series(cond).groupby(codes).nunique() > 1).any())
    # random-effect designs to try in turn: paired → a slope per condition (free 2 × 2
    # covariance), then one shared slope; between animals → one slope; finally no animal term
    designs = ([np.column_stack([t * isA, t * isB]), t[:, None]] if paired else [t[:, None]])
    n_anim = len(uniq)
    df = n_anim - 1 if paired else n_anim - 2
    b = V = None
    s2a = 0.0
    from statsmodels.regression.mixed_linear_model import MixedLM
    for Z in designs:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = MixedLM(yt, Xt, groups=codes, exog_re=Z * sw[:, None]).fit(
                    reml=True, method=["bfgs", "nm"])
            cre = np.asarray(res.cov_re)
            s2a = float(np.trace(cre))
            if not (res.converged and np.all(np.isfinite(res.fe_params))
                    and s2a > 1e-8 * res.scale):
                continue
            bb = np.asarray(res.fe_params, float)
            VV = np.asarray(res.cov_params(), float)[:p, :p]
            if VV[0, 0] + VV[1, 1] - 2 * VV[0, 1] > 0:
                b, V = bb, VV
                break
        except Exception:  # noqa: BLE001
            continue
    if b is None:                        # no animal term: WLS
        try:
            inv = np.linalg.inv(Xt.T @ Xt)
        except np.linalg.LinAlgError:
            return None
        b = inv @ (Xt.T @ yt)
        r = yt - Xt @ b
        V = inv * float(r @ r) / (n - p)
        s2a = 0.0
    dk = -(b[1] - b[0])
    se = float(np.sqrt(V[0, 0] + V[1, 1] - 2 * V[0, 1]))
    if not (se > 0) or df < 1:
        return None
    return dict(k_A=-b[0], k_B=-b[1], dk=dk, se=se, df=df, n=n, n_pep=P, animals=n_anim,
                paired=paired, sigma_animal=float(np.sqrt(max(s2a, 0.0))),
                p=float(2 * stats.t.sf(abs(dk / se), df)))


def _fit_E_chunk(args):
    items, ref, test, band, time, *rest = args
    mult = bool(rest and rest[0])
    powell = rest[1] if len(rest) > 1 else None
    out = []
    for (exp, prot), g in items:
        r = (fit_E_time_protein(g, ref, test, band, powell) if time
             else fit_E_protein(g, ref, test, band, mult))
        if r is not None:
            out.append(dict(experiment=exp, protein=prot, **r))
    return out


def variant_E(seen: dict, ref: str = REF, test: str = TEST, band: bool = False,
              workers: int | None = None, time: bool = False,
              mult: bool = False, powell: bool | None = None) -> pd.DataFrame:
    f = shared_bands(seen["frame"]) if band else seen["frame"].copy()
    f = f[np.isfinite(f.fs) & f.condition.astype(str).isin([ref, test])].copy()
    f["sig"] = _sigma(f)
    if band:
        f = f[np.isfinite(f.sig)]
    # the same proteins as the collapse tests: ≥ min_peptides precursors per condition
    npc = f.groupby(["experiment", "condition", "protein"]).concat.transform("nunique")
    f = f[npc >= seen["min_peptides"]]
    items = list(f.groupby(["experiment", "protein"], sort=False))
    workers = (workers or int(os.environ.get("WLSVAR_WORKERS", "0"))
               or max(1, (os.cpu_count() or 2) - 4))
    chunks = [items[i::workers * 4] for i in range(workers * 4)]
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for part in ex.map(_fit_E_chunk, [(c, ref, test, band, time, mult, powell)
                                                     for c in chunks]) for r in part]
    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj"] = _bh_per_experiment(out)
    return out


VARIANTS = {
    "D": lambda s, ref=REF, test=TEST: variant_D(s, ref, test, "wls"),
    "D+var": lambda s, ref=REF, test=TEST: variant_D(s, ref, test, "wls-var"),
    "E": lambda s, ref=REF, test=TEST: variant_E(s, ref, test, band=False),
    "E+band": lambda s, ref=REF, test=TEST: variant_E(s, ref, test, band=True),
    "E+time": lambda s, ref=REF, test=TEST: variant_E(s, ref, test, band=False, time=True),
    "E+band+time": lambda s, ref=REF, test=TEST: variant_E(s, ref, test, band=True, time=True),
    # E+time as an implementation would fit it: Powell + the GLS covariance (§8 of the report)
    "E+time-powell": lambda s, ref=REF, test=TEST: variant_E(s, ref, test, band=False,
                                                             time=True, powell=True),
    "Em": lambda s, ref=REF, test=TEST: variant_E(s, ref, test, band=False, mult=True),
}


# --- E, joint over tissues: the interaction --------------------------------------------------

def _gls_cov(X: np.ndarray, z: np.ndarray, groups: np.ndarray, vc: dict, scale: float,
             s2a: float, s2c: float) -> np.ndarray | None:
    """(Σ_a X_a' V_a⁻¹ X_a)⁻¹ at the fitted variance components, V_a = scale·I + s2a·z z' +
    s2c·Z_cell Z_cell'. ``res.cov_params()`` inverts the joint Hessian over the fixed effects
    and the variance parameters, which is ill-conditioned when a component sits near zero
    (negative variances, a near-singular Wald block)."""
    p = X.shape[1]
    M = np.zeros((p, p))
    for a in np.unique(groups):
        r = np.flatnonzero(groups == a)
        Va = (scale * np.eye(len(r)) + s2a * np.outer(z[r], z[r])
              + s2c * vc[a] @ vc[a].T)
        M += X[r].T @ np.linalg.solve(Va, X[r])
    try:
        V = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        return None
    return V if np.all(np.isfinite(V)) and np.all(np.diag(V) > 0) else None


def _gls_cov_multi(X: np.ndarray, groups: np.ndarray, comps: list, scale: float
                   ) -> np.ndarray | None:
    """``_gls_cov`` for any set of components: V_g = scale·I + Σ_k s2_k Z_k[g] Z_k[g]',
    ``comps`` = [(s2_k, {g: Z_k[g]})], Z_k[g] on group g's rows in data order."""
    p = X.shape[1]
    M = np.zeros((p, p))
    for gi in np.unique(groups):
        r = np.flatnonzero(groups == gi)
        Vg = scale * np.eye(len(r))
        for s2, Z in comps:
            Vg += s2 * Z[gi] @ Z[gi].T
        M += X[r].T @ np.linalg.solve(Vg, X[r])
    try:
        V = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        return None
    return V if np.all(np.isfinite(V)) and np.all(np.diag(V) > 0) else None


def _joint_fit_animal(Xt, yt, zt, animal, cell, n_anim):
    """E, joint: groups = animal, random animal slope + the cell variance component."""
    from statsmodels.regression.mixed_linear_model import MixedLM
    # variance component "cell": per animal, one column per tissue (t·√w on that cell's rows)
    vc = {}
    for a in range(n_anim):
        rows = np.flatnonzero(animal == a)
        cells = np.unique(cell[rows])
        vc[a] = np.column_stack([np.where(cell[rows] == c, zt[rows], 0.0) for c in cells])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            order = np.argsort(animal, kind="stable")    # statsmodels wants groups contiguous
            res = MixedLM(yt[order], Xt[order], groups=animal[order],
                          exog_re=zt[order][:, None],
                          exog_vc={"cell": vc}).fit(reml=True, method=["bfgs", "nm"])
        s2a = float(np.asarray(res.cov_re)[0, 0])
        s2c = float(np.asarray(res.vcomp)[0])
        if not (res.converged and np.all(np.isfinite(res.fe_params))):
            return None
        b = np.asarray(res.fe_params, float)
        V = _gls_cov(Xt, zt, animal, vc, float(res.scale), max(s2a, 0.0), max(s2c, 0.0))
    except Exception:  # noqa: BLE001
        return None
    return None if V is None else dict(b=b, V=V, s2a=s2a, s2c=s2c, s2t=0.0, s2d=0.0)


def _joint_fit_time(Xt, yt, zt, animal, cell, sidx, tp, n_tp):
    """EJ+time: groups = labeling time, random time slope; components animal, cell and
    tissue × time, all nested in the time."""
    from statsmodels.regression.mixed_linear_model import MixedLM
    tissue_day = sidx * n_tp + tp

    def by_day(unit):
        out = {}
        for d in range(n_tp):
            rows = np.flatnonzero(tp == d)
            out[d] = np.column_stack([np.where(unit[rows] == u, zt[rows], 0.0)
                                      for u in np.unique(unit[rows])])
        return out

    vcs = {"animal": by_day(animal), "cell": by_day(cell), "tday": by_day(tissue_day)}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            order = np.argsort(tp, kind="stable")
            # BFGS stalls at its start with three nested components and ≈ 5 groups; Powell
            # converges and agrees with a long Nelder–Mead run (ATF6 spot check)
            res = MixedLM(yt[order], Xt[order], groups=tp[order], exog_re=zt[order][:, None],
                          exog_vc=vcs).fit(reml=True, method=["powell", "nm"], maxiter=2000)
        s2t = float(np.asarray(res.cov_re)[0, 0])
        s2a, s2c, s2d = (float(v) for v in np.asarray(res.vcomp))     # sorted names
        if not (res.converged and np.all(np.isfinite(res.fe_params))):
            return None
        b = np.asarray(res.fe_params, float)
        zday = {d: zt[np.flatnonzero(tp == d)][:, None] for d in range(n_tp)}
        V = _gls_cov_multi(Xt, tp, [(max(s2t, 0.0), zday), (max(s2a, 0.0), vcs["animal"]),
                                    (max(s2c, 0.0), vcs["cell"]), (max(s2d, 0.0), vcs["tday"])],
                           float(res.scale))
    except Exception:  # noqa: BLE001
        return None
    return None if V is None else dict(b=b, V=V, s2a=s2a, s2c=s2c, s2t=s2t, s2d=s2d)


def fit_E_joint_protein(g: pd.DataFrame, ref: str, test: str,
                        time: bool = False) -> dict | None:
    """The interaction model of ``riana.core.interaction_model`` at the peptide level::

        φ = Σ_s β_s·t·[s] + β_shift·t·[test] + Σ_j γ_j·t·[test]·C[s, j] − b_{p,s}·t
            + u_animal·t + v_{animal,s}·t + ε

    b_{p,s}: fixed peptide slope offsets, sum-to-zero within each tissue; u: the animal slope
    shared by its tissues; v: the cell (animal × tissue) slope, one variance; weights
    (1 − θ̂)² from an OLS pass. Tests as shipped (mean Δk, additive γ = 0, ratio scale by the
    delta method). df: mean Δk animals − 1 (between animals); additive and ratio
    cells − animals − 2(S − 1) (cells within animals), the containment rule at the cell level.

    ``time`` (EJ+time): plus a random slope per labeling time shared by every sample taken then
    (all tissues, both conditions) and a tissue × time slope (the part the tissues do not
    share); time ⊃ animal ⊃ cell. Cross-sectional designs with ≥ 3 times; otherwise, or if
    that fit fails, as without (``time_term`` False). df: mean Δk animals − times − 1;
    additive and ratio cells − animals − (S − 1)·times − (S − 1) (the tissue × condition ×
    time stratum pooled with the cells; balanced-design counts).
    """
    th = np.minimum(g.fs.to_numpy(float), 0.999)
    t = g.labeling_time.to_numpy(float)
    phi = np.log1p(-np.maximum(th, -0.5))
    keep = np.isfinite(phi) & (t > 0) & (phi > PHI_LIMIT)
    g, phi, t = g[keep].copy(), phi[keep], t[keep]
    g["animal"] = _animal_keys(g)
    ok = (g.groupby(["experiment", "condition"]).animal.nunique().unstack()
          .reindex(columns=[ref, test]).fillna(0))     # a missing condition counts as 0
    exps = sorted(e for e in ok.index if ok.loc[e].min() >= 2)
    if len(exps) < 2:
        return None
    m = g.experiment.isin(exps).to_numpy()
    g, phi, t = g[m], phi[m], t[m]
    S = len(exps)
    sidx = g.experiment.map({e: i for i, e in enumerate(exps)}).to_numpy()
    is_test = (g.condition.astype(str) == test).to_numpy(float)
    C = np.vstack([np.eye(S - 1), -np.ones(S - 1)])
    cols = ([t * (sidx == s) for s in range(S)] + [t * is_test]
            + [t * is_test * C[sidx, j] for j in range(S - 1)])
    for s in range(S):
        ms = sidx == s
        pep = pd.factorize(g.concat[ms])[0]
        P = int(pep.max()) + 1
        if P > 1:
            Cp = np.vstack([np.eye(P - 1), -np.ones(P - 1)])
            for j in range(P - 1):
                col = np.zeros(len(t))
                col[ms] = t[ms] * Cp[pep, j]
                cols.append(col)
    X = np.column_stack(cols)
    y = phi
    n, p = X.shape
    animal, uniq = pd.factorize(g.animal)
    cell = pd.factorize(g.animal + "|" + g.experiment.astype(str))[0]
    n_anim, n_cell = len(uniq), int(cell.max()) + 1
    if n <= p + 3 or n_cell <= n_anim:
        return None
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    w = np.exp(2 * np.maximum(X @ b0, PHI_LIMIT))
    sw = np.sqrt(w / w.mean())
    Xt, yt, zt = X * sw[:, None], y * sw, t * sw
    tp, tp_uniq = pd.factorize(pd.Series(t))
    n_tp = len(tp_uniq)
    if time and (n_tp < 3 or (pd.Series(tp).groupby(animal).nunique() > 1).any()):
        time = False
    fit = _joint_fit_time(Xt, yt, zt, animal, cell, sidx, tp, n_tp) if time else None
    if fit is None:
        time = False
        fit = _joint_fit_animal(Xt, yt, zt, animal, cell, n_anim)
    if fit is None:
        return None
    b, V = fit["b"], fit["V"]
    ref_slope = b[:S]
    test_slope = ref_slope + b[S] + C @ b[S + 1:2 * S]
    if time:
        df_mean = n_anim - n_tp - 1
        df_within = n_cell - n_anim - (S - 1) * n_tp - (S - 1)
    else:
        df_mean = n_anim - 1
        df_within = n_cell - n_anim - 2 * (S - 1)
    if df_within < 1 or df_mean < 1:
        return None
    se_shift = float(np.sqrt(V[S, S]))
    mean_p = float(2 * stats.t.sf(abs(b[S] / se_shift), df_mean)) if se_shift > 0 else np.nan
    gi = slice(S + 1, 2 * S)

    def wald(est, cov):
        try:
            return float(est @ np.linalg.solve(cov, est) / len(est))
        except np.linalg.LinAlgError:
            return float("nan")

    add_F = wald(b[gi], V[gi, gi])
    add_p = float(stats.f.sf(add_F, S - 1, df_within)) if np.isfinite(add_F) else np.nan
    rat_p = np.nan
    if np.all(ref_slope < 0) and np.all(test_slope < 0):
        J = np.zeros((S, p))
        for s in range(S):
            J[s, s] = 1.0 / test_slope[s] - 1.0 / ref_slope[s]
            J[s, S] = 1.0 / test_slope[s]
            J[s, S + 1:2 * S] = C[s] / test_slope[s]
        D = np.hstack([np.eye(S - 1), -np.ones((S - 1, 1))])
        r = np.log(test_slope / ref_slope)
        rat_F = wald(D @ r, D @ J @ V @ J.T @ D.T)
        rat_p = float(stats.f.sf(rat_F, S - 1, df_within)) if np.isfinite(rat_F) else np.nan
    out = dict(n_experiments=S, n_animals=n_anim,
               sigma_animal=float(np.sqrt(max(fit["s2a"], 0))),
               sigma_cell=float(np.sqrt(max(fit["s2c"], 0))), mean_delta_k=-float(b[S]),
               mean_delta_k_p=mean_p, interaction_additive_p=add_p, interaction_p=rat_p)
    if time or fit["s2t"]:
        out.update(sigma_time=float(np.sqrt(max(fit["s2t"], 0))),
                   sigma_tissue_time=float(np.sqrt(max(fit["s2d"], 0))))
    return out


def _fit_EJ_chunk(args):
    items, ref, test, *rest = args          # (items, ref, test[, time]): old callers pass 3
    time = bool(rest and rest[0])
    out = []
    for prot, g in items:
        r = fit_E_joint_protein(g, ref, test, time=time)
        if r is not None:
            out.append(dict(protein=prot, **r))
    return out


def interaction_E(seen: dict, ref: str = REF, test: str = TEST,
                  workers: int | None = None, time: bool = False) -> pd.DataFrame:
    """Peptide-level interaction per protein (``time``: EJ+time); BH per test across
    proteins."""
    from statsmodels.stats.multitest import multipletests
    f = seen["frame"]
    f = f[np.isfinite(f.fs) & f.condition.astype(str).isin([ref, test])].copy()
    npc = f.groupby(["experiment", "condition", "protein"]).concat.transform("nunique")
    f = f[npc >= seen["min_peptides"]]
    items = list(f.groupby("protein", sort=False))
    workers = (workers or int(os.environ.get("WLSVAR_WORKERS", "0"))
               or max(1, (os.cpu_count() or 2) - 4))
    chunks = [items[i::workers * 4] for i in range(workers * 4)]
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for part in ex.map(_fit_EJ_chunk, [(c, ref, test, time) for c in chunks])
                for r in part]
    out = pd.DataFrame(rows)
    for col in ("mean_delta_k_p", "interaction_p", "interaction_additive_p"):
        out[col + "_adj"] = np.nan
        ok = out[col].notna()
        if ok.any():
            out.loc[ok, col + "_adj"] = multipletests(out.loc[ok, col], method="fdr_bh")[1]
    return out


# --- E with a timepoint term ---------------------------------------------------------------

def fit_E_time_protein(g: pd.DataFrame, ref: str, test: str, band: bool,
                       powell: bool | None = None) -> dict | None:
    """E plus a random slope per labeling time, shared by every animal (both conditions)
    sampled then: the timepoint's common deviation from the line (curve shape, a shared
    precursor history, a batch), which cancels from Δk when both conditions are sampled at
    each time. Nested: time ⊃ animal (⊃ cell, when animals carry both conditions; cross-
    sectional designs only). Falls back to E when the time term cannot be fitted.
    df: animals − timepoints − 1 (between animals), animals − 1 (paired)."""
    th = np.minimum(g.fs.to_numpy(float), 0.999)
    t = g.labeling_time.to_numpy(float)
    phi = np.log1p(-np.maximum(th, -0.5))
    keep = np.isfinite(phi) & (t > 0) & (phi > PHI_LIMIT)
    g, phi, t = g[keep], phi[keep], t[keep]
    cond = g.condition.astype(str).to_numpy()
    isA, isB = (cond == ref), (cond == test)
    animal = _animal_keys(g)
    if len(np.unique(animal[isA])) < 2 or len(np.unique(animal[isB])) < 2:
        return None
    a_codes, a_uniq = pd.factorize(pd.Series(animal))
    tp_codes, tp_uniq = pd.factorize(pd.Series(t))
    # cross-sectional only: every animal at one time
    if (pd.Series(tp_codes).groupby(a_codes).nunique() > 1).any() or len(tp_uniq) < 3:
        return fit_E_protein(g, ref, test, band)
    pep = pd.factorize(g.concat)[0]
    P = int(pep.max()) + 1
    C = np.vstack([np.eye(P - 1), -np.ones(P - 1)]) if P > 1 else np.zeros((1, 0))
    X = np.column_stack([t * isA, t * isB] + [t * C[pep, j] for j in range(P - 1)])
    y = phi
    n, p = X.shape
    if n <= p + 3:
        return None
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    w = np.exp(2 * np.maximum(X @ b0, PHI_LIMIT))
    if band:
        w = w / g.sig.to_numpy(float) ** 2
    sw = np.sqrt(w / w.mean())
    Xt, yt, zt = X * sw[:, None], y * sw, t * sw
    paired = bool((pd.Series(cond).groupby(a_codes).nunique() > 1).any())
    cell = pd.factorize(pd.Series(animal) + "|" + pd.Series(cond))[0]
    n_anim, n_tp = len(a_uniq), len(tp_uniq)

    def vc_of(unit):
        out = {}
        for gi in range(n_tp):
            rows = np.flatnonzero(tp_codes == gi)
            units = np.unique(unit[rows])
            out[gi] = np.column_stack([np.where(unit[rows] == u, zt[rows], 0.0) for u in units])
        return out

    vc = {"animal": vc_of(a_codes)}
    if paired:
        vc["cell"] = vc_of(cell)
    from statsmodels.regression.mixed_linear_model import MixedLM
    # WLSVAR_ET_OPT=powell: Powell + the GLS covariance, as EJ+time (the implementation
    # proposal); the default (BFGS + the Hessian covariance) is what the cached results used
    if powell is None:
        powell = os.environ.get("WLSVAR_ET_OPT", "") == "powell"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if powell:
                res = MixedLM(yt, Xt, groups=tp_codes, exog_re=zt[:, None],
                              exog_vc=vc).fit(reml=True, method=["powell", "nm"], maxiter=2000)
            else:
                res = MixedLM(yt, Xt, groups=tp_codes, exog_re=zt[:, None],
                              exog_vc=vc).fit(reml=True, method=["bfgs", "nm"])
        b = np.asarray(res.fe_params, float)
        if powell:
            s2v = np.asarray(res.vcomp, float)                      # sorted: animal, cell
            zday = {d: zt[np.flatnonzero(tp_codes == d)][:, None] for d in range(n_tp)}
            comps = [(max(float(np.asarray(res.cov_re)[0, 0]), 0.0), zday)]
            comps += [(max(float(v), 0.0), vc[k]) for v, k in zip(s2v, sorted(vc))]
            V = _gls_cov_multi(Xt, tp_codes, comps, float(res.scale))
            if V is None:
                return fit_E_protein(g, ref, test, band)
        else:
            V = np.asarray(res.cov_params(), float)[:p, :p]
        var_se = V[0, 0] + V[1, 1] - 2 * V[0, 1]
        if not (res.converged and np.all(np.isfinite(b)) and var_se > 0):
            return fit_E_protein(g, ref, test, band)
        s2t = float(np.asarray(res.cov_re)[0, 0])
        vcomp = np.asarray(res.vcomp, float)
    except Exception:  # noqa: BLE001
        return fit_E_protein(g, ref, test, band)
    dk = -(b[1] - b[0])
    se = float(np.sqrt(var_se))
    df = n_anim - 1 if paired else n_anim - n_tp - 1
    if df < 1:
        return None
    return dict(k_A=-b[0], k_B=-b[1], dk=dk, se=se, df=df, n=n, n_pep=P, animals=n_anim,
                paired=paired, sigma_time=float(np.sqrt(max(s2t, 0))),
                sigma_animal=float(np.sqrt(max(vcomp[0], 0))),
                p=float(2 * stats.t.sf(abs(dk / se), df)))
