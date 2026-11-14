"""Stage-2-only refit for relabeled designs (WP0 of the wls-var / D / E evaluation).

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md`` §4 WP0. A permutation (T2) or a
replicate split (T1) changes only which curve each point belongs to; the per-point fs (the
IsoSpec envelope solve, stage 1) does not depend on the labels. This tool re-runs only stage 2
of ``riana.core.fitting._fit_one_concat`` (depth gate, ``curve_fit``, residual bootstrap → k,
R², CI and the per-point prediction band) on the points stored by a fit, regrouped by the new
labels, and writes a fit directory the shipped ``riana rollup`` reads.

- **Point store.** Use a fit made at ``-d 1``: a ``-d 3`` fit stores only the points of
  peptides that passed depth in their *original* group, and a relabeled group can need others.
  Rail-dropped points are never stored (as in a full fit).
- **Exactness.** The code mirrors the shipped stage 2 line by line (same model, bounds, initial
  value, ``n_boot``, percentiles and per-peptide seed ``_concat_seed(1337, fit key)``, same
  point order). Checked by ``--validate``: on the unpermuted labels it must reproduce the
  shipped depth-3 fit's k, R², CI and bands.
- **Shortcut.** The bootstrap is skipped when R² < ``--r2-floor`` (default 0.6 = the rollup's
  ``--rescue-r2`` default): such a peptide cannot be admitted (R² ≥ 0.8, or R² ≥ 0.6 with
  k_cv), so its CI and bands are never read. Its row is written with NaN CI.

- **``--fast``** (model ``simple`` only). The 200 bootstrap refits of a curve are solved
  together: a vectorized, bounded Newton (Gauss–Newton where the Hessian is not positive,
  step-halving on any rise in the SSE) on f = 1 − e^(−kt), started at k̂ as the shipped
  ``curve_fit`` is. The random draws stay in the shipped order (two index draws per resample,
  in a loop), so only the optimizer differs: it converges to 1e-13 instead of TRF's 1e-8
  tolerance. ``--validate`` again decides whether it reproduces the shipped fit.

- **Unfitted curves.** The shipped fit also writes a NaN row (``n_points`` 0) for every curve
  that has ≥ depth distinct timepoints among its raw rows (q-value filtered, before rail-drop)
  but < depth after it. Under ``--peptide-admission all`` such a row counts as the peptide
  appearing in that condition, so it changes admission. Rail-dropped rows are not in the store,
  so the raw presence comes from the integrate files of the source manifest (cached in
  ``<store>/raw_presence.pkl``); the Spep-passed peptides and their protein ids come from
  the depth-1 store's peptide table.

Labels come from two manifests: the source (original labels) and the relabeled one
(``wlsvar_manifests.py t1|t2``), joined on the data file.

python tests/benchmark/wlsvar_stage2.py STORE_DIR SRC_MANIFEST NEW_MANIFEST OUT_DIR [-d 3] [-W 16]
python tests/benchmark/wlsvar_stage2.py STORE_DIR SRC NEW OUT --validate REFERENCE_FIT_DIR
"""
from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from riana.core.fitting import (_DEFAULT_BOOT_CI_PCT, _DEFAULT_N_BOOT, _K_DEG_BOUNDS,
                                _K_DEG_INIT, _MODELS, _concat_seed, _fit_key, _k_cv)
from riana.io.manifest import read_manifest

BASE_SEED = 1337
KIN = dict(a_0=0.0, a_max=1.0, k_p=0.5, k_r=0.05, r_p=10.0)
CURVE = ["experiment", "condition", "fit_key"]


def label_map(src: Path, new: Path) -> dict:
    """(experiment, condition, t, biological_replicate) → (new condition, new replicate)."""
    old = {r.identity.data_file: r.identity for r in read_manifest(src, stage="integrate")}
    rel = {r.identity.data_file: r.identity for r in read_manifest(new, stage="integrate")}
    out = {}
    for f, i in old.items():
        j = rel[f]
        key = (i.experiment, i.condition, float(i.labeling_time), int(i.biological_replicate))
        val = (j.condition, int(j.biological_replicate))
        if out.setdefault(key, val) != val:
            raise SystemExit(f"point key {key} maps to two labels")
    return out


def stage2(job: tuple, model: str, depth: int, r2_floor: float) -> tuple:
    """One curve: mirrors riana.core.fitting._fit_one_concat after the fs solve."""
    key, t_fit, fs_fit = job
    model_fn = _MODELS[model]
    if np.unique(t_fit).size < depth:
        return key, None
    try:
        popt, _ = curve_fit(partial(model_fn, **KIN), t_fit, fs_fit, bounds=_K_DEG_BOUNDS,
                            p0=[_K_DEG_INIT], maxfev=2000)
    except (RuntimeError, ValueError):
        return key, None
    k_deg = float(popt[0])
    pred = np.array([model_fn(ti, k_deg=k_deg, **KIN) for ti in t_fit])
    residuals = fs_fit - pred
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((fs_fit - np.mean(fs_fit)) ** 2))
    r2 = float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot
    n_pts = len(t_fit)
    nan_band = [float("nan")] * n_pts
    if not (r2 >= r2_floor):
        return key, (k_deg, r2, float("nan"), float("nan"), float("nan"), nan_band, nan_band)
    rng = np.random.default_rng(_concat_seed(BASE_SEED, key[2]))
    boot_ks, boot_obs = [], []
    for _ in range(_DEFAULT_N_BOOT):
        fs_star = pred + residuals[rng.integers(0, n_pts, size=n_pts)]
        try:
            popt_b, _ = curve_fit(partial(model_fn, **KIN), t_fit, fs_star,
                                  bounds=_K_DEG_BOUNDS, p0=[k_deg], maxfev=2000)
        except (RuntimeError, ValueError):
            continue
        k_b = float(popt_b[0])
        boot_ks.append(k_b)
        pred_b = np.asarray(model_fn(t_fit, k_deg=k_b, **KIN), dtype=np.float64)
        boot_obs.append(pred_b + residuals[rng.integers(0, n_pts, size=n_pts)])
    if n_pts >= 2 and len(boot_ks) >= 10:
        ci_lo, ci_hi = (float(p) for p in np.percentile(boot_ks, _DEFAULT_BOOT_CI_PCT))
        sd = float(np.std(boot_ks, ddof=1))
        lo, hi = np.percentile(np.vstack(boot_obs), _DEFAULT_BOOT_CI_PCT, axis=0)
        return key, (k_deg, r2, ci_lo, ci_hi, sd, list(lo), list(hi))
    return key, (k_deg, r2, float("nan"), float("nan"), float("nan"), nan_band, nan_band)


_LO, _HI = float(_K_DEG_BOUNDS[0][0]), float(_K_DEG_BOUNDS[1][0])


def fit_k_vec(t: np.ndarray, Y: np.ndarray, k0: float, tol: float = 1e-13,
              max_iter: int = 200) -> np.ndarray:
    """Least-squares k of f = 1 − e^(−kt) for each row of Y (B × n), bounded, from k0."""
    k = np.full(Y.shape[0], float(k0))

    def sse(kk, yy):
        return np.sum((1.0 - np.exp(-np.outer(kk, t)) - yy) ** 2, axis=1)

    s = sse(k, Y)
    todo = np.ones(k.size, bool)
    for _ in range(max_iter):
        if not todo.any():
            break
        kk, yy = k[todo], Y[todo]
        e = np.exp(-np.outer(kk, t))
        r = (1.0 - e) - yy
        J = t * e                                   # ∂f/∂k
        g = np.sum(J * r, axis=1)                   # ½ ∂SSE/∂k
        gn = np.sum(J * J, axis=1)
        h = gn - np.sum(r * t * t * e, axis=1)      # + r·∂²f/∂k²  (∂²f/∂k² = −t²e)
        h = np.where(h > 0, h, gn)
        step = np.where(h > 1e-300, -g / np.where(h > 1e-300, h, 1.0), 0.0)
        lam = np.ones_like(kk)
        k_new = np.clip(kk + step, _LO, _HI)
        s_old = s[todo]
        s_new = sse(k_new, yy)
        rise = s_new > s_old
        for _h in range(60):
            if not rise.any():
                break
            lam[rise] *= 0.5
            k_new[rise] = np.clip(kk[rise] + lam[rise] * step[rise], _LO, _HI)
            s_new[rise] = sse(k_new[rise], yy[rise])
            rise = s_new > s_old
        k_new = np.where(rise, kk, k_new)
        s_new = np.where(rise, s_old, s_new)
        moved = np.abs(k_new - kk)
        idx = np.flatnonzero(todo)
        k[idx], s[idx] = k_new, s_new
        todo[idx[moved <= tol * np.maximum(kk, _LO)]] = False
    return k


def stage2_fast(job: tuple, model: str, depth: int, r2_floor: float) -> tuple:
    """As :func:`stage2`, with the bootstrap refits solved together by :func:`fit_k_vec`."""
    key, t_fit, fs_fit = job
    model_fn = _MODELS[model]
    if np.unique(t_fit).size < depth:
        return key, None
    try:
        popt, _ = curve_fit(partial(model_fn, **KIN), t_fit, fs_fit, bounds=_K_DEG_BOUNDS,
                            p0=[_K_DEG_INIT], maxfev=2000)
    except (RuntimeError, ValueError):
        return key, None
    k_deg = float(popt[0])
    pred = np.array([model_fn(ti, k_deg=k_deg, **KIN) for ti in t_fit])
    residuals = fs_fit - pred
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((fs_fit - np.mean(fs_fit)) ** 2))
    r2 = float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot
    n_pts = len(t_fit)
    nan_band = [float("nan")] * n_pts
    if not (r2 >= r2_floor):
        return key, (k_deg, r2, float("nan"), float("nan"), float("nan"), nan_band, nan_band)
    rng = np.random.default_rng(_concat_seed(BASE_SEED, key[2]))
    i1 = np.empty((_DEFAULT_N_BOOT, n_pts), np.int64)
    i2 = np.empty((_DEFAULT_N_BOOT, n_pts), np.int64)
    for b in range(_DEFAULT_N_BOOT):                # the shipped draw order
        i1[b] = rng.integers(0, n_pts, size=n_pts)
        i2[b] = rng.integers(0, n_pts, size=n_pts)
    boot_ks = fit_k_vec(t_fit, pred + residuals[i1], k_deg)
    boot_obs = 1.0 - np.exp(-np.outer(boot_ks, t_fit)) + residuals[i2]
    if n_pts >= 2:
        ci_lo, ci_hi = (float(p) for p in np.percentile(boot_ks, _DEFAULT_BOOT_CI_PCT))
        sd = float(np.std(boot_ks, ddof=1))
        lo, hi = np.percentile(boot_obs, _DEFAULT_BOOT_CI_PCT, axis=0)
        return key, (k_deg, r2, ci_lo, ci_hi, sd, list(lo), list(hi))
    return key, (k_deg, r2, float("nan"), float("nan"), float("nan"), nan_band, nan_band)


def raw_presence(store: Path, src: Path, q_value: float) -> pd.DataFrame:
    """Unique (experiment, condition, labeling_time, biological_replicate, fit_key) over the
    raw integrate rows (q < q_value, MBR kept as in the fit), under the source labels."""
    cache = store / "raw_presence.pkl"
    if cache.exists():
        return pd.read_pickle(cache)
    parts = []
    for r in read_manifest(src, stage="integrate"):
        path = Path(r.output_path)
        if not path.exists():
            path = src.parent / r.output_path
        d = pd.read_table(path, comment="#", usecols=["concat", "percolator q-value"],
                          low_memory=False)
        keys = d.loc[d["percolator q-value"] < q_value, "concat"].map(_fit_key).unique()
        i = r.identity
        parts.append(pd.DataFrame({"experiment": i.experiment, "condition": i.condition,
                                   "labeling_time": float(i.labeling_time),
                                   "biological_replicate": int(i.biological_replicate),
                                   "fit_key": keys}))
    out = pd.concat(parts, ignore_index=True).drop_duplicates()
    out.to_pickle(cache)
    return out


def _chunk(jobs: list, model: str, depth: int, r2_floor: float, fast: bool = False) -> list:
    fn = stage2_fast if fast else stage2
    return [fn(j, model, depth, r2_floor) for j in jobs]


def header(path: Path, depth: int) -> str:
    lines = [ln for ln in open(path) if ln.startswith("#")]
    lines = [f"# depth {depth}\n" if ln.startswith("# depth ") else ln for ln in lines]
    return "".join(lines) + "# stage2_refit tests/benchmark/wlsvar_stage2.py\n"


def run(store: Path, labels: dict, out: Path, depth: int, workers: int, r2_floor: float,
        fast: bool = False, src: Path | None = None) -> None:
    t0 = time.time()
    frac = pd.read_table(store / "riana_fit_fractions.txt", comment="#", low_memory=False)
    pep = pd.read_table(store / "riana_fit_peptides.txt", comment="#", low_memory=False)
    model = next(ln.split()[2] for ln in open(store / "riana_fit_peptides.txt")
                 if ln.startswith("# model "))
    if fast and model != "simple":
        raise SystemExit(f"--fast supports model simple only (store has {model})")
    keys = list(zip(frac.experiment, frac.condition, frac.labeling_time.astype(float),
                    frac.biological_replicate.astype(int)))
    new = [labels[k] for k in keys]
    frac["condition"] = [c for c, _ in new]
    frac["biological_replicate"] = [b for _, b in new]
    # curves the shipped fit would list unfitted: ≥ depth raw timepoints under the new labels
    unfit = []
    if src is not None:
        q = float(next(ln.split()[2] for ln in open(store / "riana_fit_peptides.txt")
                       if ln.startswith("# q_value ")))
        raw = raw_presence(store, src, q)
        raw = raw[raw.fit_key.isin(set(pep.concat))]          # passed the Spep gate
        nl = [labels[k] for k in zip(raw.experiment, raw.condition, raw.labeling_time,
                                     raw.biological_replicate)]
        raw = raw.assign(condition=[c for c, _ in nl])
        nt = raw.groupby(CURVE).labeling_time.nunique()
        unfit = list(nt.index[nt >= depth])
    frac["fit_key"] = frac.concat.map(_fit_key)
    groups = frac.groupby(CURVE, sort=False).indices
    jobs = [(k, frac.labeling_time.to_numpy(float)[ix], frac.fs.to_numpy(float)[ix])
            for k, ix in groups.items()]
    size = max(1, len(jobs) // (workers * 8))
    chunks = [jobs[i:i + size] for i in range(0, len(jobs), size)]
    fn = partial(_chunk, model=model, depth=depth, r2_floor=r2_floor, fast=fast)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = dict(r for part in ex.map(fn, chunks) for r in part)
    lo = np.full(len(frac), np.nan)
    hi = np.full(len(frac), np.nan)
    keep = np.zeros(len(frac), bool)
    rows = []
    first = frac.groupby(CURVE, sort=False).first()
    spep = pep.drop_duplicates("concat").set_index("concat").spep
    for k, ix in groups.items():
        r = res.get(k)
        if r is None:
            continue
        k_deg, r2, ci_lo, ci_hi, sd, blo, bhi = r
        keep[ix], lo[ix], hi[ix] = True, blo, bhi
        g = frac.iloc[ix]
        mbr = (g.evidence.astype(str) == "mbr").to_numpy()
        metox = g.metox.astype(bool).to_numpy()
        f0 = first.loc[k]
        rows.append({"concat": k[2], "k_deg": k_deg, "R_squared": r2, "sd": sd,
                     "spep": spep.get(k[2], np.nan), "ci_lo": ci_lo, "ci_hi": ci_hi,
                     "k_cv": _k_cv(ci_lo, ci_hi, k_deg), "protein id": f0["protein id"],
                     "mod sites": f0["mod sites"], "n_points": len(ix),
                     "n_mbr": int(mbr.sum()), "n_metox": int(metox.sum()),
                     "n_clean": int((~mbr & ~metox).sum()), "experiment": k[0],
                     "condition": k[1]})
    fitted = {(r["experiment"], r["condition"], r["concat"]) for r in rows}
    prot = pep.drop_duplicates("concat").set_index("concat")
    for e, c, fk in unfit:
        if (e, c, fk) in fitted:
            continue
        rows.append({"concat": fk, "k_deg": np.nan, "R_squared": np.nan, "sd": np.nan,
                     "spep": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "k_cv": np.nan,
                     "protein id": prot.at[fk, "protein id"], "mod sites": np.nan,
                     "n_points": 0, "n_mbr": 0, "n_metox": 0, "n_clean": 0,
                     "experiment": e, "condition": c})
    frac["fs_lower"], frac["fs_upper"] = lo, hi
    out.mkdir(parents=True, exist_ok=True)
    for name, df in [("riana_fit_peptides.txt", pd.DataFrame(rows)),
                     ("riana_fit_fractions.txt", frac[keep].drop(columns="fit_key"))]:
        with open(out / name, "w") as fh:
            fh.write(header(store / name, depth))
            df.to_csv(fh, sep="\t", index=False)
    print(f"{out}: {len(fitted):,} curves fitted from {len(jobs):,} groups, "
          f"{len(rows) - len(fitted):,} unfitted rows, in {(time.time() - t0) / 60:.1f} min")


def validate(out: Path, ref: Path, r2_floor: float) -> None:
    on = ["experiment", "condition", "concat"]
    a = pd.read_table(out / "riana_fit_peptides.txt", comment="#")
    b = pd.read_table(ref / "riana_fit_peptides.txt", comment="#")
    ua = a[a.k_deg.isna()].set_index(on)
    ub = b[b.k_deg.isna()].set_index(on)
    both = ua.index.intersection(ub.index)
    print(f"unfitted rows: stage2 {len(ua):,}, reference {len(ub):,}, both {len(both):,}; "
          f"protein id equal {bool((ua.loc[both, 'protein id'] == ub.loc[both, 'protein id']).all())}")
    a = a.dropna(subset=["k_deg"]).set_index(on)
    b = b.dropna(subset=["k_deg"]).set_index(on)
    print(f"curves: stage2 {len(a):,}, reference {len(b):,}, both {len(a.index.intersection(b.index)):,}")
    c = a.index.intersection(b.index)
    for col in ["k_deg", "R_squared"]:
        print(f"  {col}: max |diff| {np.nanmax(np.abs(a.loc[c, col] - b.loc[c, col])):.3g}")
    g = c[(b.loc[c, "R_squared"] >= r2_floor).to_numpy()]
    for col in ["ci_lo", "ci_hi"]:
        print(f"  {col} (R² ≥ {r2_floor}, {len(g):,}): max |diff| "
              f"{np.nanmax(np.abs(a.loc[g, col] - b.loc[g, col])):.3g}")
    fa = pd.read_table(out / "riana_fit_fractions.txt", comment="#", low_memory=False)
    fb = pd.read_table(ref / "riana_fit_fractions.txt", comment="#", low_memory=False)
    k = ["experiment", "condition", "concat", "labeling_time", "biological_replicate"]
    for f in (fa, fb):  # a key can hold several points (e.g. charge states): number them in order
        f["nth"] = f.groupby(k).cumcount()
    m = fa.merge(fb, on=k + ["nth"], suffixes=("_a", "_b"))
    print(f"  fs of matched points: max |diff| {np.nanmax(np.abs(m.fs_a - m.fs_b)):.3g}")
    m = m[m.fs_lower_a.notna()]
    print(f"  points: stage2 {len(fa):,}, reference {len(fb):,}, compared {len(m):,}; band max |diff| "
          f"{np.nanmax(np.abs(m.fs_lower_a - m.fs_lower_b)):.3g} / "
          f"{np.nanmax(np.abs(m.fs_upper_a - m.fs_upper_b)):.3g}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("store", type=Path)
    p.add_argument("src", type=Path)
    p.add_argument("new", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("-d", "--depth", type=int, default=3)
    p.add_argument("-W", "--workers", type=int, default=8)
    p.add_argument("--r2-floor", type=float, default=0.6)
    p.add_argument("--validate", type=Path, help="reference fit directory to compare against")
    p.add_argument("--fast", action="store_true", help="vectorized bootstrap (model simple)")
    a = p.parse_args()
    run(a.store, label_map(a.src, a.new), a.out, a.depth, a.workers, a.r2_floor, a.fast,
        a.src)
    if a.validate:
        validate(a.out, a.validate, a.r2_floor)


if __name__ == "__main__":
    sys.exit(main())
