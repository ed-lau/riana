"""FS rail-drop threshold sweep — is the drop bound ideal, or should it be tighter?

The clamp-only rails (the solver's ±``FS_BOUNDS``) drop only points the solver *pinned*
to its bound. A solved-but-implausible FS (1.15 over-labelled, −0.07 sub-natural) is
also not a measurement. This sweeps ``FitConfig.fs_rail_hi`` / ``fs_rail_lo`` (rail-drop
ON) and, for each threshold, reports admitted yield, R², and the **matched** geom-CV
vs the clamp baseline (same peptides).

Findings (reports/2026-07-05_multipoint_rail_drop.md):
  - ``-0.1 / 1.1`` is an artifact: its lower −0.1 sits *below* the solver's −0.09999
    low-clamp, so it stops dropping the lower rail. The lower threshold must stay above
    the clamp convergence.
  - ``0.0 / 1.0`` (the exact physical bound) eats AC16's near-zero t0 anchors → −16%
    yield. Dataset-dangerous.
  - ``-0.05 / 1.05`` (a 0.05 margin around [0,1]) is the robust optimum and the 1.2.0
    default — cleaner R²/CV on every set, +5–9% yield on the noisy ones, a no-op on ¹⁸O.

Run the ¹⁸O confirmation with ``--runs boomi_ipsc_o18 juber_ac16_o18`` (it is a no-op:
no ¹⁸O FS lands in the (1.05, 1.199] band).

Usage:
    python -m tests.benchmark.bench_fs_rail_threshold                 # 3 D₂O sets
    python -m tests.benchmark.bench_fs_rail_threshold --runs boomi_ipsc_o18 juber_ac16_o18
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from riana.config import FitConfig
from riana.core.fitting import load_aa_coefficients, load_o18_coefficients
from riana.core.pipeline import fit_project

RUNS = {
    "boomi_ipsc_d2o": ("alamillo_2025_ipsc", "hw"),
    "juber_ac16_d2o": ("alamillo_2025_ac16", "hw"),
    "lve_atr_clean":  ("deberneh_2025_rss",  "hw"),
    "boomi_ipsc_o18": ("juber_2026_o18_ac16", "o18"),
    "juber_ac16_o18": ("juber_2026_o18_ac16", "o18"),
}
D2O = ["boomi_ipsc_d2o", "juber_ac16_d2o", "lve_atr_clean"]
# (label, hi, lo). None = the shipped default rails (1.05 / -0.05).
CONFIGS = [
    ("clamp",         1.199, -0.099),
    ("hi1.1_lo-.1",   1.1,   -0.1),
    ("hi1.05_lo-.05", 1.05,  -0.05),
    ("hi1.0_lo0.0",   1.0,   0.0),
]
MIN_R2, RESCUE_R2, K_CV = 0.8, 0.6, 0.2


def admit_mask(p):
    r2 = p["R_squared"].to_numpy(float)
    kcv = p["k_cv"].to_numpy(float) if "k_cv" in p else np.full(len(p), np.nan)
    return (r2 >= MIN_R2) | ((r2 >= RESCUE_R2) & (kcv < K_CV))


def geomcv(p, mp=3):
    d = p[np.isfinite(p["k_deg"]) & (p["k_deg"] > 0)]
    d = d[~d["protein id"].astype(str).str.contains(r"[;,]")]
    cv = [np.sqrt(np.expm1(np.var(np.log(g["k_deg"].to_numpy()), ddof=1)))
          for _, g in d.groupby("protein id") if len(g) >= mp]
    return float(np.median(cv)) if cv else float("nan")


def fit_cfg(manifest, coeffs, label, hi, lo, depth, workers, n_boot):
    cfg = FitConfig(model="simple", label=label, q_value=1e-2, depth=depth,
                    workers=workers, fs_rail_drop=True, fs_rail_hi=hi, fs_rail_lo=lo)
    out = fit_project(cfg, manifest, coeffs, n_boot=n_boot).reset_index()
    out = out[np.isfinite(out["k_deg"])].copy()
    out["admit"] = admit_mask(out)
    return out


def run(names, depth, workers, n_boot):
    for name in names:
        coef, label = RUNS[name]
        load = load_o18_coefficients if label == "o18" else load_aa_coefficients
        coeffs = load(coef)
        manifest = f"runs/{name}/riana_manifest.tsv"
        clamp_adm = None
        base = 0
        for lbl, hi, lo in CONFIGS:
            t0 = time.time()
            out = fit_cfg(manifest, coeffs, label, hi, lo, depth, workers, n_boot)
            adm = out[out.admit]
            cv, r2all = geomcv(adm), float(np.nanmedian(out["R_squared"]))
            extra = ""
            if lbl == "clamp":
                base, clamp_adm = int(out.admit.sum()), adm
            else:
                both = set(clamp_adm.concat) & set(adm.concat)
                dcv = geomcv(adm[adm.concat.isin(both)]) - geomcv(clamp_adm[clamp_adm.concat.isin(both)])
                d = int(out.admit.sum()) - base
                extra = f" Δadmit={d:+d} ({100*d/max(base,1):+.1f}%) Δmatchedcv={dcv:+.4f}"
            print(f"[{name:15s} {lbl:14s}] admit={int(out.admit.sum()):5d} "
                  f"r2all={r2all:+.3f} geomcv={cv:.4f}{extra} ({round(time.time()-t0)}s)",
                  flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs", nargs="*", default=D2O, choices=list(RUNS))
    ap.add_argument("--depth", type=int, default=6)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()
    run(a.runs, a.depth, a.workers, a.n_boot)


if __name__ == "__main__":
    main()
