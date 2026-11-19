"""Peptide-level linear model — ``rollup --model linear --level peptide``.

The protein-level linear model (:mod:`riana.core.linear_model`) first collapses a protein's
peptides into one θ per (biological replicate, timepoint) and then fits φ = log(1 − θ)
against time. This module fits the peptides' own points instead, in one mixed model per
protein, so that the spread between a protein's peptides and the animals' shared deviations
enter the standard errors rather than the collapsed values. It is the candidate **E+time**
(and **EJ+time** for ``--interaction``) of ``reports/2026-11-14_wls_var_D_E_evaluation.md``,
chosen 2026-11-16 as the opt-in peptide-level method (report §5, §8).

**Δk model (E+time)**, per (experiment, protein), on the admitted peptides' points of the
compared conditions (t > 0, φ above ``phi_limit``)::

    φ = −k_c·t − b_p·t + u_time·t + v_animal·t [+ w_cell·t] + ε

- ``k_c``: one slope per condition; Δk = k_test − k_ref.
- ``b_p``: fixed peptide slope offsets, sum-to-zero, shared by the conditions — an apparent-k
  offset per peptide. k_c is then the peptides' unweighted average k, and a peptide's own
  offset cannot enter Δk.
- ``u_time``: a random slope per labeling time, shared by every sample taken then (the
  timepoint's common deviation from the line: curve shape, a shared precursor history, a
  batch). It cancels from Δk when both conditions are sampled at each time.
- ``v_animal``: a random slope per animal, nested in its time. The animal is the
  ``individual`` column (SDRF ``characteristics[individual]``); without one, the condition ×
  time × replicate composite, which is the cell.
- ``w_cell``: when animals carry both conditions (tissues of one animal), an animal ×
  condition slope; the animal's shared part then cancels from Δk (the paired model). Pairing
  is detected from the data, not set by a flag.
- ``Var(ε) ∝ 1/w``, ``w = (1 − θ̂)²`` from a first OLS pass, as the protein-level ``wls``.

REML with statsmodels ``MixedLM`` on the √w-scaled data, groups = labeling time, Powell then
Nelder–Mead (BFGS stalls with nested components and few groups). The fixed-effect covariance
is the GLS one at the fitted variance components (:func:`_gls_cov_multi`); statsmodels'
``cov_params`` inverts the joint Hessian, which is ill-conditioned when a component sits near
zero. Denominator df: animals − 1 when paired, else animals − times − 1.

**Fallback (E)** when the design is not cross-sectional (an animal at several times), has
fewer than three times, or the time model fails: the random slope is per animal only
(groups = animal; paired: a slope per condition with a free covariance, then one shared
slope), BFGS then Nelder–Mead, df animals − 1 (paired) or animals − conditions. When that
fails too, or the animal variance is at 0: weighted least squares with no random term.
``random_effects`` records the path (``time``, ``animal`` or ``none``).

**Conditions.** The Δk contrast pair: ``--reference-condition`` / ``--test-condition``, or
the experiment's two conditions. A condition of the pair that a protein lacks (or that has
fewer than two animals left) leaves that protein with the other condition's k only, fitted
by the same model on one condition. An experiment with more than two conditions needs the
pair named. Δk rests on the peptides seen in both conditions (each one's offset cancels);
a peptide in one condition only sharpens that condition's k and the variance components. A
protein whose two conditions share no peptide has no Δk (k per condition only), and a
design whose fixed effects are not all identified is not fitted.

**Per-condition k CI.** From the same covariance, at df = times − 1 when the time slope is in
the model (the time stratum carries k's variance) and at the Δk df otherwise. The Δk test is
what the evaluation calibrated; the per-condition interval is the matching containment
interval, not separately validated.

φ = ln(1 − θ) with θ capped at 0.999 and floored at −0.5 (not at the protein model's 0.001:
a peptide θ below 0 gives φ > 0, and flooring single noisy points would bias them upward).
"""

from __future__ import annotations

import logging
import warnings
from concurrent import futures
from typing import Callable

import numpy as np
import pandas as pd
from scipy import stats

from riana.multiplex import peptide_identity
from riana.progress import iter_progress

_LOGGER = logging.getLogger(__name__)

#: Columns ``--level peptide`` adds to the linear output (``riana_rollup_proteins.txt``):
#: the Δk df, the animals in the fit, whether animals carried both conditions, the
#: random-effects path (``time`` / ``animal`` / ``none``) and the fitted SDs.
PEPTIDE_LEVEL_COLUMNS = [
    "delta_k_df", "n_animals", "paired", "random_effects", "sigma_time", "sigma_animal",
]
#: Columns the peptide-level interaction adds to :data:`riana.core.interaction_model.
#: INTERACTION_COLUMNS`.
PEPTIDE_INTERACTION_COLUMNS = ["random_effects", "sigma_cell", "sigma_time",
                               "sigma_tissue_time"]
#: Per-condition k interval: 95 %, as the protein-level model's ``conf_int``.
_CI_ALPHA = 0.05
_THETA_CEILING = 0.999
_THETA_FLOOR = -0.5


def peptide_phi(fs: np.ndarray) -> np.ndarray:
    """φ = ln(1 − θ) for peptide points: θ capped at 0.999, floored at −0.5."""
    th = np.minimum(np.asarray(fs, dtype=float), _THETA_CEILING)
    return np.log1p(-np.maximum(th, _THETA_FLOOR))


def _animal_keys(g: pd.DataFrame) -> np.ndarray:
    """The animal of each point: ``individual`` where named, else the condition × time ×
    biological-replicate composite (the cell)."""
    from riana.core.protein import individual_label
    comp = ("c:" + g["condition"].astype(str) + "|" + g["labeling_time"].astype(str) + "|"
            + g["biological_replicate"].astype(str))
    if "individual" in g:
        named = g["individual"].map(individual_label).astype(str)
        return np.where(named != "", "i:" + named, comp)
    return comp.to_numpy()


def _shares_peptide(g: pd.DataFrame, a: str, b: str) -> bool:
    """Whether conditions ``a`` and ``b`` hold at least one peptide in common — what a
    peptide-level Δk rests on (a peptide seen in one condition only enters that
    condition's k and the variance components, not Δk)."""
    cond = g["condition"].astype(str)
    return bool(set(g["concat"][cond == a]) & set(g["concat"][cond == b]))


def _pos_def(V: np.ndarray) -> bool:
    return bool(np.all(np.isfinite(V)) and np.all(np.linalg.eigvalsh((V + V.T) / 2) > 0))


def _by_identity(fractions: pd.DataFrame) -> pd.DataFrame:
    """``concat`` → the peptide's identity across samples (channel-label tokens removed,
    :func:`riana.multiplex.peptide_identity`): the peptide offsets, the shared-peptide
    check and the peptide counts then treat a dimethyl light and heavy form as one."""
    return fractions.assign(concat=peptide_identity(fractions["concat"]))


def _full_rank(X: np.ndarray) -> bool:
    """Every fixed effect identified (a guard against designs the checks above miss)."""
    return bool(np.linalg.matrix_rank(X) == X.shape[1])


def _sum_to_zero(pep: np.ndarray) -> tuple[np.ndarray, int]:
    """Sum-to-zero contrast matrix for ``pep`` codes and the number of peptides."""
    P = int(pep.max()) + 1
    C = np.vstack([np.eye(P - 1), -np.ones(P - 1)]) if P > 1 else np.zeros((1, 0))
    return C, P


def _gls_cov(X: np.ndarray, z: np.ndarray, groups: np.ndarray, vc: dict, scale: float,
             s2a: float, s2c: float) -> np.ndarray | None:
    """(Σ_a X_a' V_a⁻¹ X_a)⁻¹ at the fitted variance components, V_a = scale·I + s2a·z z' +
    s2c·Z_cell Z_cell'."""
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
    """GLS covariance for any set of components: V_g = scale·I + Σ_k s2_k Z_k[g] Z_k[g]',
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


# --------------------------------------------------------------------------- #
# Δk: one protein
# --------------------------------------------------------------------------- #
def _slope_var(V: np.ndarray, K: int) -> float:
    """Variance of the Δk contrast (K = 2), or of the one slope (K = 1)."""
    return float(V[0, 0] + V[1, 1] - 2 * V[0, 1]) if K == 2 else float(V[0, 0])


def _fit_animal(g: pd.DataFrame, conds: list[str], phi_limit: float) -> dict | None:
    """E: random slope per animal (the fallback; see the module docstring)."""
    phi = peptide_phi(g["fs"].to_numpy(float))
    t = g["labeling_time"].to_numpy(float)
    cond = g["condition"].astype(str).to_numpy()
    K = len(conds)
    is_c = [cond == c for c in conds]
    animal = _animal_keys(g)
    if any(len(np.unique(animal[m])) < 2 for m in is_c):
        return None
    pep = pd.factorize(g["concat"])[0]
    C, P = _sum_to_zero(pep)
    X = np.column_stack([t * m for m in is_c] + [t * C[pep, j] for j in range(P - 1)])
    y = phi
    n, p = X.shape
    if n <= p + 2 or not _full_rank(X):
        return None
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    w = np.exp(2 * np.maximum(X @ b0, phi_limit))
    w = w / w.mean()
    sw = np.sqrt(w)
    Xt, yt = X * sw[:, None], y * sw
    codes, uniq = pd.factorize(pd.Series(animal))
    paired = bool((pd.Series(cond).groupby(codes).nunique() > 1).any())
    # random-effect designs in turn: paired → a slope per condition (free covariance), then
    # one shared slope; between animals → one slope; finally no animal term
    designs = ([np.column_stack([t * m for m in is_c]), t[:, None]] if paired
               else [t[:, None]])
    n_anim = len(uniq)
    df = n_anim - 1 if paired else n_anim - K
    b = V = None
    s2a = 0.0
    random = "none"
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
            # statsmodels' inverse Hessian can be indefinite near the variance boundary
            # (a negative per-condition variance, a Δk SE from a broken matrix): such a
            # fit falls through to the next design, finally WLS.
            if _slope_var(VV, K) > 0 and _pos_def(VV[:K, :K]):
                b, V, random = bb, VV, "animal"
                break
        except Exception:  # noqa: BLE001 - statsmodels raises various
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
    if not (_slope_var(V, K) > 0) or df < 1:
        return None
    return dict(b=b[:K], V=V[:K, :K], df=df, df_k=df, n=[int(m.sum()) for m in is_c],
                n_pep=P, n_animals=n_anim, paired=paired, random_effects=random,
                sigma_time=0.0, sigma_animal=float(np.sqrt(max(s2a, 0.0))))


def _fit_time(g: pd.DataFrame, conds: list[str], phi_limit: float) -> dict | None:
    """E+time (see the module docstring); falls back to :func:`_fit_animal`."""
    phi = peptide_phi(g["fs"].to_numpy(float))
    t = g["labeling_time"].to_numpy(float)
    keep = np.isfinite(phi) & (t > 0) & (phi > phi_limit)
    g, phi, t = g[keep], phi[keep], t[keep]
    cond = g["condition"].astype(str).to_numpy()
    K = len(conds)
    is_c = [cond == c for c in conds]
    animal = _animal_keys(g)
    if any(len(np.unique(animal[m])) < 2 for m in is_c):
        return None
    a_codes, a_uniq = pd.factorize(pd.Series(animal))
    tp_codes, tp_uniq = pd.factorize(pd.Series(t))
    # cross-sectional only: every animal at one time
    if (pd.Series(tp_codes).groupby(a_codes).nunique() > 1).any() or len(tp_uniq) < 3:
        return _fit_animal(g, conds, phi_limit)
    pep = pd.factorize(g["concat"])[0]
    C, P = _sum_to_zero(pep)
    X = np.column_stack([t * m for m in is_c] + [t * C[pep, j] for j in range(P - 1)])
    y = phi
    n, p = X.shape
    if n <= p + 3 or not _full_rank(X):
        return None
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    w = np.exp(2 * np.maximum(X @ b0, phi_limit))
    sw = np.sqrt(w / w.mean())
    Xt, yt, zt = X * sw[:, None], y * sw, t * sw
    paired = bool((pd.Series(cond).groupby(a_codes).nunique() > 1).any())
    cell = pd.factorize(pd.Series(animal) + "|" + pd.Series(cond))[0]
    n_anim, n_tp = len(a_uniq), len(tp_uniq)

    def vc_of(unit: np.ndarray) -> dict:
        out = {}
        for gi in range(n_tp):
            rows = np.flatnonzero(tp_codes == gi)
            units = np.unique(unit[rows])
            out[gi] = np.column_stack([np.where(unit[rows] == u, zt[rows], 0.0)
                                       for u in units])
        return out

    vc = {"animal": vc_of(a_codes)}
    if paired:
        vc["cell"] = vc_of(cell)
    from statsmodels.regression.mixed_linear_model import MixedLM
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = MixedLM(yt, Xt, groups=tp_codes, exog_re=zt[:, None],
                          exog_vc=vc).fit(reml=True, method=["powell", "nm"], maxiter=2000)
        b = np.asarray(res.fe_params, float)
        s2v = np.asarray(res.vcomp, float)                      # sorted: animal, cell
        zday = {d: zt[np.flatnonzero(tp_codes == d)][:, None] for d in range(n_tp)}
        comps = [(max(float(np.asarray(res.cov_re)[0, 0]), 0.0), zday)]
        comps += [(max(float(v), 0.0), vc[k]) for v, k in zip(s2v, sorted(vc))]
        V = _gls_cov_multi(Xt, tp_codes, comps, float(res.scale))
        if V is None:
            return _fit_animal(g, conds, phi_limit)
        if not (res.converged and np.all(np.isfinite(b)) and _slope_var(V, K) > 0):
            return _fit_animal(g, conds, phi_limit)
        s2t = float(np.asarray(res.cov_re)[0, 0])
        vcomp = np.asarray(res.vcomp, float)
    except Exception:  # noqa: BLE001 - statsmodels raises various
        return _fit_animal(g, conds, phi_limit)
    df = n_anim - 1 if paired else n_anim - n_tp - 1
    if K == 2 and df < 1:
        return None
    return dict(b=b[:K], V=V[:K, :K], df=df, df_k=n_tp - 1,
                n=[int(m.sum()) for m in is_c], n_pep=P, n_animals=n_anim, paired=paired,
                random_effects="time", sigma_time=float(np.sqrt(max(s2t, 0))),
                sigma_animal=float(np.sqrt(max(vcomp[0], 0))))


def _rows_for(exp: str, prot: str, g: pd.DataFrame, pair: list[str],
              phi_limit: float) -> list[dict]:
    """Fit one protein's points of ``pair`` (the contrast, or one condition) → output rows."""
    phi = peptide_phi(g["fs"].to_numpy(float))
    t = g["labeling_time"].to_numpy(float)
    keep = np.isfinite(phi) & (t > 0) & (phi > phi_limit)
    animal = _animal_keys(g)
    cond = g["condition"].astype(str).to_numpy()
    # A condition needs two animals left after the point filters; else it is dropped.
    conds = [c for c in pair if len(np.unique(animal[keep & (cond == c)])) >= 2]
    if not conds:
        return []
    if len(conds) == 2 and not _shares_peptide(g[keep], *conds):
        # Δk is a within-peptide contrast (each peptide's offset is shared by the two
        # conditions): with no peptide in both, it is not identified — k per condition only.
        return [r for c in conds
                for r in _rows_for(exp, prot, g[cond == c], [c], phi_limit)]
    sub = g[np.isin(cond, conds)]
    fit = _fit_time(sub, conds, phi_limit)
    if fit is None and len(conds) == 2:          # no Δk: each condition's k alone
        return [r for c in conds
                for r in _rows_for(exp, prot, g[cond == c], [c], phi_limit)]
    if fit is None:
        return []
    b, V, K = fit["b"], fit["V"], len(conds)
    dk = dse = dp = np.nan
    if K == 2:
        dk = float(-(b[1] - b[0]))
        dse = float(np.sqrt(_slope_var(V, 2)))
        dp = float(2 * stats.t.sf(abs(dk / dse), fit["df"]))
    tq = (float(stats.t.ppf(1 - _CI_ALPHA / 2, fit["df_k"])) if fit["df_k"] >= 1
          else np.nan)
    rows = []
    for i, c in enumerate(conds):
        k, se = float(-b[i]), float(np.sqrt(V[i, i]))
        rows.append({
            "experiment": exp, "condition": c, "protein": prot,
            "n_points": fit["n"][i], "k_deg": k,
            "ci_lo": k - tq * se, "ci_hi": k + tq * se, "R_squared": np.nan,
            "delta_k": dk, "delta_k_se": dse, "delta_k_p": dp,
            "delta_k_df": float(fit["df"]) if K == 2 else np.nan,
            "n_animals": fit["n_animals"], "paired": fit["paired"],
            "random_effects": fit["random_effects"],
            "sigma_time": fit["sigma_time"], "sigma_animal": fit["sigma_animal"],
        })
    return rows


# --------------------------------------------------------------------------- #
# Δk: driver
# --------------------------------------------------------------------------- #
_WORKER_STATE: dict = {}


def _init_worker(items, phi_limit) -> None:
    _WORKER_STATE["items"] = items
    _WORKER_STATE["phi_limit"] = phi_limit


def _deltak_worker(index: int) -> list[dict]:
    exp, prot, g, pair = _WORKER_STATE["items"][index]
    return _rows_for(exp, prot, g, pair, _WORKER_STATE["phi_limit"])


def _check_pair_conditions(fractions: pd.DataFrame, reference: str | None,
                           test: str | None) -> None:
    from riana.exceptions import DataError
    present = set(fractions["condition"].astype(str).unique())
    for role, cond in (("reference", reference), ("test", test)):
        if cond is not None and str(cond) not in present:
            raise DataError(
                f"{role} condition {cond!r} is not among the conditions in the data "
                f"({sorted(present)}). Check the spelling / the SDRF condition values.")
    if reference is not None and test is not None and str(reference) == str(test):
        raise DataError(
            f"reference and test condition are the same ({reference!r}) — pick two "
            "different conditions.")


def _experiment_pairs(fractions: pd.DataFrame, reference: str | None,
                      test: str | None) -> dict[str, list[str]]:
    """The conditions fitted per experiment: the contrast pair (or the one condition)."""
    from riana.exceptions import DataError
    out: dict[str, list[str]] = {}
    for exp, conds in fractions.groupby("experiment")["condition"].agg(
            lambda s: sorted(set(s.astype(str)))).items():
        if reference is not None and test is not None:
            out[str(exp)] = [str(reference), str(test)]
        elif len(conds) == 1:
            out[str(exp)] = conds
        elif len(conds) == 2:
            ref = str(reference) if reference is not None and reference in conds else conds[0]
            out[str(exp)] = [ref, next(c for c in conds if c != ref)]
        else:
            raise DataError(
                f"--level peptide fits the two conditions of the Δk contrast; experiment "
                f"{exp!r} has {len(conds)} ({conds}). Name the pair with "
                "--reference-condition and --test-condition.")
    return out


def _admitted(fractions: pd.DataFrame, conds_by_exp: dict[str, list[str]],
              min_peptides: int, keys: list[str]) -> pd.DataFrame:
    """Finite points of the fitted conditions, from (experiment, condition, protein)
    groups with at least ``min_peptides`` peptides (as the collapse requires)."""
    f = fractions[np.isfinite(fractions["fs"].to_numpy(float))].copy()
    f["condition"] = f["condition"].astype(str)
    f["experiment"] = f["experiment"].astype(str)
    wanted = f["experiment"].map(conds_by_exp)
    f = f[[c in w if isinstance(w, list) else False
           for c, w in zip(f["condition"], wanted)]]
    npc = f.groupby(keys)["concat"].transform("nunique")
    return f[npc >= min_peptides]


def fit_peptide_deltak(
    fractions: pd.DataFrame,
    *,
    reference_condition: str | None = None,
    test_condition: str | None = None,
    phi_limit: float = -4.0,
    min_peptides: int = 1,
    workers: int = 1,
    progress_callback: "Callable[[int, int], None] | None" = None,
) -> pd.DataFrame:
    """Per-protein peptide-level k per condition and the Δk test (E+time).

    Args:
        fractions: the admitted peptide points (``riana_fit_fractions.txt`` after the
            rollup's attribution and curation gates): ``experiment``, ``condition``,
            ``protein``, ``concat``, ``labeling_time``, ``biological_replicate``, ``fs``
            and optionally ``individual``.
        reference_condition / test_condition: the contrast, as for
            :func:`riana.core.linear_model.fit_linear_deltak`; Δk = k(test) − k(reference).
            Without ``test_condition`` an experiment needs at most two conditions.
        phi_limit: plateau limit; points with φ at or below it are dropped.
        min_peptides: a protein's condition needs this many peptides.
        workers: processes (statsmodels is GIL-bound). Results do not depend on it.
        progress_callback: ``callback(done, total)`` per protein.

    Returns:
        One row per (experiment, protein, condition): the linear output columns
        (:data:`riana.core.linear_model.LINEAR_COLUMNS`, ``R_squared`` NaN) plus
        :data:`PEPTIDE_LEVEL_COLUMNS`. ``delta_k_p_adj`` is BH per experiment.
    """
    from riana.core.linear_model import LINEAR_COLUMNS
    from riana.exceptions import DataError

    need = {"experiment", "condition", "protein", "concat", "labeling_time",
            "biological_replicate", "fs"}
    missing = need - set(fractions.columns)
    if missing:
        raise DataError(f"peptide-level points missing columns {sorted(missing)}.")
    columns = [*LINEAR_COLUMNS, *PEPTIDE_LEVEL_COLUMNS]
    if fractions.empty:
        return pd.DataFrame(columns=columns)
    fractions = _by_identity(fractions)
    _check_pair_conditions(fractions, reference_condition, test_condition)
    pairs = _experiment_pairs(fractions, reference_condition, test_condition)
    f = _admitted(fractions, pairs, min_peptides, ["experiment", "condition", "protein"])
    items = [(str(exp), str(prot), g, pairs[str(exp)])
             for (exp, prot), g in f.groupby(["experiment", "protein"], sort=False)]
    if workers > 1 and len(items) > 1:
        chunk = max(1, len(items) // (workers * 8))
        with futures.ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker, initargs=(items, phi_limit),
        ) as ex:
            computed = list(iter_progress(
                ex.map(_deltak_worker, range(len(items)), chunksize=chunk),
                len(items), progress_callback))
    else:
        computed = list(iter_progress(
            (_rows_for(e, p, g, pr, phi_limit) for e, p, g, pr in items),
            len(items), progress_callback))
    out = pd.DataFrame([r for rows in computed for r in rows],
                       columns=[c for c in columns if c != "delta_k_p_adj"])
    out["delta_k_p_adj"] = np.nan
    if out.empty:
        return out[columns]
    tested = out.dropna(subset=["delta_k_p"]).drop_duplicates(["experiment", "protein"])
    if not tested.empty:
        paths = tested["random_effects"].value_counts().to_dict()
        _LOGGER.info(
            "--level peptide: %d Δk tests; random effects by path %s; %d paired "
            "(animals under both conditions).", len(tested), paths,
            int(tested["paired"].sum()))
        from statsmodels.stats.multitest import multipletests
        adj: dict[tuple, float] = {}
        for exp, grp in tested.groupby("experiment"):
            padj = multipletests(grp["delta_k_p"].to_numpy(), method="fdr_bh")[1]
            adj.update(zip(zip(grp["experiment"], grp["protein"]), padj))
        out["delta_k_p_adj"] = [adj.get((e, p), np.nan)
                                for e, p in zip(out["experiment"], out["protein"])]
    return out[columns]


# --------------------------------------------------------------------------- #
# Interaction (EJ+time): one protein
# --------------------------------------------------------------------------- #
def _joint_fit_animal(Xt, yt, zt, animal, cell, n_anim) -> dict | None:
    """EJ: groups = animal, random animal slope + the cell variance component."""
    from statsmodels.regression.mixed_linear_model import MixedLM
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
    except Exception:  # noqa: BLE001 - statsmodels raises various
        return None
    return None if V is None else dict(b=b, V=V, s2a=s2a, s2c=s2c, s2t=0.0, s2d=0.0)


def _joint_fit_time(Xt, yt, zt, animal, cell, sidx, tp, n_tp) -> dict | None:
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
    except Exception:  # noqa: BLE001 - statsmodels raises various
        return None
    return None if V is None else dict(b=b, V=V, s2a=s2a, s2c=s2c, s2t=s2t, s2d=s2d)


def _interaction_rows(prot: str, g: pd.DataFrame, ref: str, test: str,
                      phi_limit: float) -> list[dict]:
    """The experiment × condition interaction for one protein, at the peptide level::

        φ = Σ_s β_s·t·[s] + β_shift·t·[test] + Σ_j γ_j·t·[test]·C[s, j] − b_{p,s}·t
            + u_time·t + v_animal·t + w_cell·t + x_{tissue×time}·t + ε

    b_{p,s}: fixed peptide slope offsets, sum-to-zero within each experiment; u: a slope per
    labeling time shared by every sample taken then; v: the animal slope shared by its
    experiments; w: the cell (animal × experiment) slope; x: the experiment × time slope.
    Tests as :func:`riana.core.interaction_model.fit_interaction` (mean Δk, additive γ = 0,
    ratio scale by the delta method) at containment df: mean Δk animals − times − 1;
    additive and ratio cells − animals − (S − 1)·times − (S − 1). Without the time term (an
    animal at several times, fewer than three times, or no convergence): groups = animal,
    df animals − 1 and cells − animals − 2(S − 1).
    """
    phi = peptide_phi(g["fs"].to_numpy(float))
    t = g["labeling_time"].to_numpy(float)
    keep = np.isfinite(phi) & (t > 0) & (phi > phi_limit)
    g, phi, t = g[keep].copy(), phi[keep], t[keep]
    g["animal"] = _animal_keys(g)
    ok = (g.groupby(["experiment", "condition"])["animal"].nunique().unstack()
          .reindex(columns=[ref, test]).fillna(0))     # a missing condition counts as 0
    # an experiment enters when both conditions have two animals and share a peptide
    # (its Δk is a within-peptide contrast, as in the Δk model)
    exps = sorted(e for e in ok.index if ok.loc[e].min() >= 2
                  and _shares_peptide(g[g["experiment"] == e], ref, test))
    if len(exps) < 2:
        return []
    m = g["experiment"].isin(exps).to_numpy()
    g, phi, t = g[m], phi[m], t[m]
    S = len(exps)
    sidx = g["experiment"].map({e: i for i, e in enumerate(exps)}).to_numpy()
    is_test = (g["condition"].astype(str) == test).to_numpy(float)
    C = np.vstack([np.eye(S - 1), -np.ones(S - 1)])
    cols = ([t * (sidx == s) for s in range(S)] + [t * is_test]
            + [t * is_test * C[sidx, j] for j in range(S - 1)])
    for s in range(S):
        ms = sidx == s
        pep = pd.factorize(g["concat"][ms])[0]
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
    animal, uniq = pd.factorize(g["animal"])
    cell = pd.factorize(g["animal"] + "|" + g["experiment"].astype(str))[0]
    n_anim, n_cell = len(uniq), int(cell.max()) + 1
    if n <= p + 3 or n_cell <= n_anim or not _full_rank(X):
        return []
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    w = np.exp(2 * np.maximum(X @ b0, phi_limit))
    sw = np.sqrt(w / w.mean())
    Xt, yt, zt = X * sw[:, None], y * sw, t * sw
    tp, tp_uniq = pd.factorize(pd.Series(t))
    n_tp = len(tp_uniq)
    time = not (n_tp < 3 or (pd.Series(tp).groupby(animal).nunique() > 1).any())
    fit = _joint_fit_time(Xt, yt, zt, animal, cell, sidx, tp, n_tp) if time else None
    if fit is None:
        time = False
        fit = _joint_fit_animal(Xt, yt, zt, animal, cell, n_anim)
    if fit is None:
        return []
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
        return []
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
    rat_F = rat_p = np.nan
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
    protein_level = {
        "n_experiments": S, "n_animals": n_anim,
        "sigma_animal": float(np.sqrt(max(fit["s2a"], 0))),
        "mean_delta_k": -float(b[S]), "mean_delta_k_se": se_shift,
        "mean_delta_k_df": float(df_mean), "mean_delta_k_p": mean_p,
        "interaction_F": rat_F, "interaction_df1": float(S - 1),
        "interaction_df2": float(df_within), "interaction_p": rat_p,
        "interaction_additive_F": add_F, "interaction_additive_df2": float(df_within),
        "interaction_additive_p": add_p,
        "random_effects": "time" if time else "animal",
        "sigma_cell": float(np.sqrt(max(fit["s2c"], 0))),
        "sigma_time": float(np.sqrt(max(fit["s2t"], 0))) if time else np.nan,
        "sigma_tissue_time": float(np.sqrt(max(fit["s2d"], 0))) if time else np.nan,
    }
    rows = []
    for s, name in enumerate(exps):
        l_d = np.zeros(p)
        l_d[S] = 1.0
        l_d[S + 1:2 * S] = C[s]
        rows.append({
            "protein": prot, "experiment": name, "n_points": int((sidx == s).sum()),
            "k_reference": -float(ref_slope[s]), "k_test": -float(test_slope[s]),
            "delta_k": -float(l_d @ b), "delta_k_se": float(np.sqrt(l_d @ V @ l_d)),
            "k_ratio": float(test_slope[s] / ref_slope[s]),
            **protein_level,
        })
    return rows


def _init_inter_worker(items, args) -> None:
    _WORKER_STATE["inter_items"] = items
    _WORKER_STATE["inter_args"] = args


def _inter_worker(index: int) -> list[dict]:
    prot, g = _WORKER_STATE["inter_items"][index]
    return _interaction_rows(prot, g, *_WORKER_STATE["inter_args"])


def _animal_key(f: pd.DataFrame) -> str:
    """How animals were identified (for provenance), as the protein-level interaction."""
    if "individual" not in f.columns:
        return "composite"
    from riana.core.protein import individual_label
    named = f["individual"].map(individual_label) != ""
    if not named.any():
        return "composite"
    return "individual" if named.all() else "individual+composite"


def fit_peptide_interaction(
    fractions: pd.DataFrame,
    *,
    reference_condition: str,
    test_condition: str,
    phi_limit: float = -4.0,
    min_peptides: int = 1,
    workers: int = 1,
    progress_callback: "Callable[[int, int], None] | None" = None,
) -> pd.DataFrame:
    """Per-protein experiment × condition interaction at the peptide level (EJ+time).

    Takes the admitted peptide points (as :func:`fit_peptide_deltak`) and returns one row
    per (protein, experiment) in :data:`riana.core.interaction_model.INTERACTION_COLUMNS`
    plus :data:`PEPTIDE_INTERACTION_COLUMNS`; BH across proteins per test.
    ``attrs["animal_key"]`` records how animals were identified.
    """
    from riana.core.interaction_model import INTERACTION_COLUMNS
    from riana.exceptions import DataError

    need = {"experiment", "condition", "protein", "concat", "labeling_time",
            "biological_replicate", "fs"}
    missing = need - set(fractions.columns)
    if missing:
        raise DataError(f"peptide-level points missing columns {sorted(missing)}.")
    ref, test = str(reference_condition), str(test_condition)
    fractions = _by_identity(fractions)
    _check_pair_conditions(fractions, ref, test)
    columns = [*INTERACTION_COLUMNS, *PEPTIDE_INTERACTION_COLUMNS]
    pairs = {str(e): [ref, test] for e in fractions["experiment"].astype(str).unique()}
    f = _admitted(fractions, pairs, min_peptides, ["experiment", "condition", "protein"])
    items = [(str(prot), g) for prot, g in f.groupby("protein", sort=False)]
    args = (ref, test, phi_limit)
    if workers > 1 and len(items) > 1:
        chunk = max(1, len(items) // (workers * 8))
        with futures.ProcessPoolExecutor(
            max_workers=workers, initializer=_init_inter_worker, initargs=(items, args),
        ) as ex:
            computed = list(iter_progress(
                ex.map(_inter_worker, range(len(items)), chunksize=chunk),
                len(items), progress_callback))
    else:
        computed = list(iter_progress(
            (_interaction_rows(prot, g, *args) for prot, g in items),
            len(items), progress_callback))
    adjusted = {"mean_delta_k_p_adj": "mean_delta_k_p", "interaction_p_adj": "interaction_p",
                "interaction_additive_p_adj": "interaction_additive_p"}
    out = pd.DataFrame([r for rows in computed for r in rows],
                       columns=[c for c in columns if c not in adjusted])
    for adj in adjusted:
        out[adj] = np.nan
    out.attrs["animal_key"] = _animal_key(f)
    if out.empty:
        return out[columns]
    from statsmodels.stats.multitest import multipletests
    per_protein = out.drop_duplicates("protein").set_index("protein")
    for adj, raw in adjusted.items():
        pv = per_protein[raw].dropna()
        if not pv.empty:
            padj = dict(zip(pv.index, multipletests(pv.to_numpy(), method="fdr_bh")[1]))
            out[adj] = out["protein"].map(padj)
    return out[columns]
