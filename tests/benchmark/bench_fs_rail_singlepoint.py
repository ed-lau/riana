"""FS rail-drop on a SINGLE-timepoint set — is the 1.05/−0.05 rail harmful there?

Single-timepoint curation (auto-detected in ``rollup_proteins``) rides on the replicate
floor (``min_fit_points`` auto-2) + ``k_cv``, with R² bypassed — so a tighter rail could
thin bioreps below the 2-replicate floor. This fits the 16-plex TMTpro AC16 24 h set
three ways (rail-drop OFF; ON with the old clamp; ON with the 1.05/−0.05 default) and
runs the real single-tp rollup.

Finding (reports/2026-07-05_multipoint_rail_drop.md): the tighter rail is **neutral** vs
clamp (admit −0.5%, proteins −2, median k identical). Rail-drop ON vs OFF drops the
count more (removing fake-replicate peptides whose railed FS manufactured a spurious
k_cv = 0) — the intended behaviour, not a regression.

Usage:
    python -m tests.benchmark.bench_fs_rail_singlepoint
    python -m tests.benchmark.bench_fs_rail_singlepoint --run splatd2o_tmt16 \
        --coefficients alamillo_2025_ac16
"""
from __future__ import annotations

import argparse
import time

import numpy as np

from riana.config import FitConfig
from riana.core.fitting import load_aa_coefficients
from riana.core.pipeline import fit_project
from riana.core.protein import rollup_proteins

# (label, fs_rail_drop, hi, lo). None hi/lo = the shipped default (1.05 / -0.05).
ARMS = [
    ("OFF",          False, None,  None),
    ("ON clamp",     True,  1.199, -0.099),
    ("ON 1.05/-.05", True,  None,  None),
]


def run(run_name, coefficients, workers, n_boot):
    coeffs = load_aa_coefficients(coefficients)
    manifest = f"runs/{run_name}/riana_manifest.tsv"
    for lbl, drop, hi, lo in ARMS:
        t0 = time.time()
        cfg = FitConfig(model="simple", label="hw", q_value=1e-2, depth=1, ria_max=0.06,
                        min_spep=8, workers=workers, fs_rail_drop=drop,
                        fs_rail_hi=hi, fs_rail_lo=lo)
        out = fit_project(cfg, manifest, coeffs, n_boot=n_boot)
        frac, pep = out.attrs["fractions_long"], out.reset_index()
        fit = pep[np.isfinite(pep["k_deg"])]
        reps2 = fit[fit["n_points"] >= 2]
        admit = reps2[reps2["k_cv"] < 0.2]  # single-tp gate: >=2 reps AND k_cv<0.2
        prot = rollup_proteins(pep, frac, model="simple", method="weighted",
                               parsimony="unique", min_peptides=2, k_cv_max=0.2,
                               workers=workers)
        print(f"[{lbl:13s}] fit={len(fit):5d}  n>=2reps={len(reps2):5d}  "
              f"admit(k_cv<0.2)={len(admit):5d}  proteins={len(prot):4d}  "
              f"k_med={np.nanmedian(admit['k_deg']):.4f}  ({round(time.time()-t0)}s)",
              flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", default="splatd2o_tmt16")
    ap.add_argument("--coefficients", default="alamillo_2025_ac16")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()
    run(a.run, a.coefficients, a.workers, a.n_boot)


if __name__ == "__main__":
    main()
