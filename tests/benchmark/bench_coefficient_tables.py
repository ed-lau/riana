# -*- coding: utf-8 -*-
"""Compare D₂O labelling-site coefficient tables on real in-vivo turnover data.

Two table-quality metrics, per coefficient preset, from a manifest fit (intensity
FS, the production θ_ΔI path):

1. **Within-protein robust geometric CV of k** — a protein's peptides should agree
   on one turnover rate, so the spread of per-peptide k *within* a (protein,
   condition) is a coefficient-table quality signal (a better Spep table → tighter
   agreement). CV = ``1.4826·MAD(ln k)`` over curated (R²≥``--min-r2``) peptides;
   reported as the median over (protein, condition) groups with ≥``--min-pep``
   peptides. Lower = better.
2. **fs vs fs_ds agreement** — median(fs_ds − fs) and within-0.1 fraction (the
   mass-defect cross-check; a better table should also tighten this).

Usage:
  python -m tests.benchmark.bench_coefficient_tables \
      runs/lve_fixed_ab/riana_manifest.tsv --ria 0.045 \
      --tables commerford_1983 ilchenko_2019 deberneh_2025_rss
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from riana.config import FitConfig
from riana.core.fitting import load_aa_coefficients
from riana.core.pipeline import fit_project


def geo_cv_ln(k: np.ndarray) -> float:
    k = k[np.isfinite(k) & (k > 0)]
    if len(k) < 2:
        return np.nan
    lnk = np.log(k)
    return 1.4826 * np.median(np.abs(lnk - np.median(lnk)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--ria", type=float, default=0.045)
    ap.add_argument("--tables", nargs="+",
                    default=["commerford_1983", "ilchenko_2019", "deberneh_2025_rss"])
    ap.add_argument("--min-r2", type=float, default=0.9)
    ap.add_argument("--min-pep", type=int, default=3)
    a = ap.parse_args()

    print(f"manifest={a.manifest}  ria={a.ria}  min-r2={a.min_r2}  min-pep={a.min_pep}")
    print(f"\n{'table':18} {'med within-prot k geoCV':>24} {'n groups':>9} "
          f"{'fs_ds−fs bias':>14} {'|Δ|<0.1':>8} {'n pep':>7}")
    for tbl in a.tables:
        cfg = FitConfig(model="simple", label="hw", q_value=0.01, depth=4,
                        ria_max=a.ria)
        out = fit_project(cfg, a.manifest, load_aa_coefficients(tbl),
                          n_boot=0, random_state=1)
        df = out.reset_index()
        cur = df[df["R_squared"] >= a.min_r2].copy()
        grp_cols = [c for c in ("protein id", "condition") if c in cur.columns]
        cvs = []
        for _, g in cur.groupby(grp_cols):
            if len(g) >= a.min_pep:
                cv = geo_cv_ln(g["k_deg"].to_numpy(dtype=float))
                if np.isfinite(cv):
                    cvs.append(cv)
        frac = out.attrs["fractions_long"]
        m = frac.dropna(subset=["fs_ds", "fs"])
        diff = (m["fs_ds"] - m["fs"]).to_numpy()
        print(f"{tbl:18} {np.median(cvs):>24.4f} {len(cvs):>9} "
              f"{np.median(diff):>+14.3f} {100*np.mean(np.abs(diff)<0.1):>7.1f}% "
              f"{int(df['k_deg'].notna().sum()):>7}")


if __name__ == "__main__":
    main()
