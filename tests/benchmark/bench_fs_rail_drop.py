"""FS rail-drop A/B — the multi-point rail-drop re-validation harness.

Fits each turnover project both with and without the FS rail-drop
(``FitConfig.fs_rail_drop``) and reports the deltas that decided the 1.2.0 default:
peptide yield after the production rollup gate, R² of the fitted population,
within-protein geom-CV, and — the unconfounded metric — the **matched** geom-CV over
peptides admitted in BOTH arms (does rail-drop degrade the peptides already there, or
is any CV change just the newly-rescued admits?).

The load-bearing detail is ``--depth``: at the loose ``--depth 3`` default rail-drop
appears to worsen CV (it thins a curve to the 2-point floor), but at the production
``--depth 6`` curation the curves stay well-sampled and the effect reverses. Always
run at the curation depth actually used. See
``reports/2026-07-05_multipoint_rail_drop.md``.

Usage:
    python -m tests.benchmark.bench_fs_rail_drop            # all 5 sets, depth 6
    python -m tests.benchmark.bench_fs_rail_drop --runs boomi_ipsc_d2o --depth 3
    python -m tests.benchmark.bench_fs_rail_drop -o /tmp/rail_ab.tsv
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from riana.config import FitConfig
from riana.core.fitting import load_aa_coefficients, load_o18_coefficients
from riana.core.pipeline import fit_project

#: name -> (coefficient table, label). The five standard turnover benchmark sets.
RUNS = {
    "boomi_ipsc_d2o": ("alamillo_2025_ipsc", "hw"),
    "juber_ac16_d2o": ("alamillo_2025_ac16", "hw"),
    "lve_atr_clean":  ("deberneh_2025_rss",  "hw"),
    "boomi_ipsc_o18": ("juber_2026_o18_ac16", "o18"),
    "juber_ac16_o18": ("juber_2026_o18_ac16", "o18"),
}
# Production rollup admission gate (rollup --min-r2 0.8 + the k_cv rescue).
MIN_R2, RESCUE_R2, K_CV = 0.8, 0.6, 0.2


def admit_mask(pep: pd.DataFrame) -> np.ndarray:
    r2 = pep["R_squared"].to_numpy(float)
    kcv = pep["k_cv"].to_numpy(float) if "k_cv" in pep else np.full(len(pep), np.nan)
    return (r2 >= MIN_R2) | ((r2 >= RESCUE_R2) & (kcv < K_CV))


def geomcv(pep: pd.DataFrame, min_pep: int = 3) -> tuple[float, int]:
    """Median over proteins (>= min_pep unique-accession peptides) of the geometric
    CV of peptide k_deg. Lower = peptides of one protein agree more."""
    d = pep[np.isfinite(pep["k_deg"]) & (pep["k_deg"] > 0)]
    d = d[~d["protein id"].astype(str).str.contains(r"[;,]")]
    cvs = [np.sqrt(np.expm1(np.var(np.log(g["k_deg"].to_numpy()), ddof=1)))
           for _, g in d.groupby("protein id") if len(g) >= min_pep]
    return (float(np.median(cvs)) if cvs else float("nan")), len(cvs)


def fit_arm(manifest, coeffs, label, rail_drop, depth, workers, n_boot):
    cfg = FitConfig(model="simple", label=label, q_value=1e-2, depth=depth,
                    workers=workers, fs_rail_drop=rail_drop)
    out = fit_project(cfg, manifest, coeffs, n_boot=n_boot).reset_index()
    out = out[np.isfinite(out["k_deg"])].copy()
    out["admit"] = admit_mask(out)
    return out


def run(names, depth, workers, n_boot, out_path):
    rows = []
    for name in names:
        coef, label = RUNS[name]
        load = load_o18_coefficients if label == "o18" else load_aa_coefficients
        coeffs = load(coef)
        manifest = f"runs/{name}/riana_manifest.tsv"
        arms = {}
        for rail in (False, True):
            t0 = time.time()
            arms[rail] = fit_arm(manifest, coeffs, label, rail, depth, workers, n_boot)
            print(f"[{name} rail={'ON ' if rail else 'OFF'}] "
                  f"fitted={len(arms[rail])} admit={int(arms[rail].admit.sum())} "
                  f"({round(time.time()-t0)}s)", flush=True)
        off, on = arms[False], arms[True]
        oa, na = off[off.admit], on[on.admit]
        cv_off, _ = geomcv(oa)
        cv_on, _ = geomcv(na)
        both = set(oa.concat) & set(na.concat)
        cv_m_off, npm = geomcv(oa[oa.concat.isin(both)])
        cv_m_on, _ = geomcv(na[na.concat.isin(both)])
        rows.append(dict(
            run=name, admit_off=int(off.admit.sum()), admit_on=int(on.admit.sum()),
            r2all_off=float(np.nanmedian(off["R_squared"])),
            r2all_on=float(np.nanmedian(on["R_squared"])),
            kmed_off=float(np.nanmedian(oa["k_deg"])),
            kmed_on=float(np.nanmedian(na["k_deg"])),
            geomcv_off=cv_off, geomcv_on=cv_on,
            matchedcv_off=cv_m_off, matchedcv_on=cv_m_on, n_matched=len(both)))
        r = rows[-1]
        print(f"  admit {r['admit_off']}->{r['admit_on']} "
              f"({100*(r['admit_on']-r['admit_off'])/max(r['admit_off'],1):+.1f}%)  "
              f"r2all {r['r2all_off']:.3f}->{r['r2all_on']:.3f}  "
              f"matchedcv {cv_m_off:.4f}->{cv_m_on:.4f} (Δ{cv_m_on-cv_m_off:+.4f})  "
              f"Δkmed={r['kmed_on']-r['kmed_off']:+.4f}", flush=True)
    df = pd.DataFrame(rows)
    if out_path:
        df.to_csv(out_path, sep="\t", index=False)
        print(f"\nwrote {out_path}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs", nargs="*", default=list(RUNS), choices=list(RUNS))
    ap.add_argument("--depth", type=int, default=6)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("-o", "--output-dir", dest="out", default=None,
                    help="Optional TSV path for the per-run result table.")
    a = ap.parse_args()
    run(a.runs, a.depth, a.workers, a.n_boot, a.out)


if __name__ == "__main__":
    main()
