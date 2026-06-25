# -*- coding: utf-8 -*-
"""PROTOTYPE (v1.1.0 item 1b, NOT wired into production): the mass-defect θ_ΔS
second estimate vs the intensity FS θ_ΔI.

Reads a ``riana_fit_fractions.txt`` (which already carries the per-(peptide,
sample), M0-internal, theory-referenced ``dspacing_iso{k}``) and builds θ_ΔS:

  θ_ΔS(f) = weighted_median_{k∈1..3}[ ΔSₓ(f,k) / ΔSₓmax(k) ]   (MAD outlier guard)

where ΔSₓ(f,k) is **f0/t0-anchored** (subtract the peptide's own unlabeled-row
dspacing) when an unlabeled row exists, else the theory-referenced dspacing
(fallback); and ΔSₓmax(k) is the IsoSpec theoretical init→final averaged-spacing
change (the normalizer). iso0–3 only (iso4/5 are QC-only). It then correlates θ_ΔS
against θ_ΔI (the ``fs`` column) and — for calibration — the known mixing
proportion (``labeling_time``).

Usage:
  python -m tests.benchmark.bench_mass_defect_theta \
      runs/calib_ac16_v1/riana_fit_fractions.txt --coefficients ac16 --ria 0.06
  ... add --no-anchor to use the theory fallback everywhere (A/B the anchor).
"""
from __future__ import annotations

import argparse
import re

import numpy as np
import pandas as pd

from riana.algorithms.isotope_dist import (
    _binned_envelope, get_peptide_distribution, spep_from_coefficients,
)
from riana.algorithms.mass_calc import calculate_ion_mz
from riana.core.fitting import load_aa_coefficients

CHANNELS = (1, 2, 3)          # iso0–3 spacing channels (iso0 spacing ≡ 0)
DSMAX_FLOOR = 0.2             # mDa; skip channels with negligible theoretical range


def delta_s_max(seq: str, charge: int, coeffs: dict, ria: float, n: int = 4) -> dict:
    """Theoretical ΔSₓmax(k) in m/z mDa — init→final averaged-spacing change."""
    pm = calculate_ion_mz(seq)
    spep = max(1, int(round(spep_from_coefficients(seq, coeffs))))
    init = get_peptide_distribution(seq, label=1)
    final = get_peptide_distribution(
        seq, deuterium_enrichment_level=ria, label=1, num_labeling_sites=spep)
    im, _ = _binned_envelope(init, pm, n + 1)
    fm, _ = _binned_envelope(final, pm, n + 1)
    return {k: ((fm[k] - fm[0]) - (im[k] - im[0])) / charge * 1e3
            for k in range(1, n + 1)}


def weighted_median(vals: np.ndarray, wts: np.ndarray) -> float:
    order = np.argsort(vals)
    v, w = np.asarray(vals)[order], np.asarray(wts)[order]
    c = np.cumsum(w)
    return float(v[np.searchsorted(c, w.sum() / 2.0)])


def theta_ds(ds_by_k: dict, dsmax_by_k: dict) -> float:
    ks = [k for k in CHANNELS
          if abs(dsmax_by_k.get(k, 0.0)) > DSMAX_FLOOR
          and not np.isnan(ds_by_k.get(k, np.nan))]
    if len(ks) < 2:
        return float("nan")
    th = np.array([ds_by_k[k] / dsmax_by_k[k] for k in ks])
    w = np.array([abs(dsmax_by_k[k]) for k in ks])           # signal-range weights
    med = np.median(th)
    mad = np.median(np.abs(th - med)) or 1e-9
    keep = np.abs(th - med) <= 3 * mad
    if keep.sum() < 2:
        keep = np.ones_like(th, dtype=bool)
    return weighted_median(th[keep], w[keep])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("fractions")
    ap.add_argument("--coefficients", default="ac16")
    ap.add_argument("--ria", type=float, default=0.06)
    ap.add_argument("--no-anchor", action="store_true",
                    help="use theory-referenced dspacing everywhere (A/B the f0 anchor)")
    a = ap.parse_args()

    coeffs = load_aa_coefficients(a.coefficients)
    d = pd.read_table(a.fractions, comment="#")
    f0 = d["labeling_time"].min()
    # Average the unlabeled anchor per concat (replicate t0/f0 rows — bioreps —
    # collapse to one baseline spacing per peptidoform).
    base = (d[d["labeling_time"] == f0]
            .groupby("concat")[[f"dspacing_iso{k}" for k in CHANNELS]]
            .mean())

    dsmax = {}
    for c in d["concat"].unique():
        m = re.match(r"(.+)_(\d+)$", c)
        try:
            dsmax[c] = delta_s_max(re.sub(r"\[[^\]]*\]", "", m.group(1)),
                                   int(m.group(2)), coeffs, a.ria)
        except Exception:
            dsmax[c] = None

    rows = []
    for _, r in d.iterrows():
        c = r["concat"]
        dm = dsmax.get(c)
        if dm is None:
            continue
        ds = {}
        anchored = (not a.no_anchor) and (c in base.index)
        for k in CHANNELS:
            theory = r.get(f"dspacing_iso{k}", np.nan)
            ds[k] = theory - base.loc[c, f"dspacing_iso{k}"] if anchored else theory
        rows.append({"concat": c, "f": r["labeling_time"], "fs": r["fs"],
                     "theta_ds": theta_ds(ds, dm)})

    o = pd.DataFrame(rows).dropna(subset=["theta_ds", "fs"])
    o = o[np.isfinite(o["theta_ds"]) & np.isfinite(o["fs"])]
    mode = "theory-only" if a.no_anchor else "f0-anchored (theory fallback)"
    print(f"=== θ_ΔS prototype [{mode}] — {len(o)} points ===")
    print(f"corr(θ_ΔS, θ_ΔI=fs):  r={np.corrcoef(o['theta_ds'], o['fs'])[0,1]:.3f}")
    print(f"corr(θ_ΔS, known f):  r={np.corrcoef(o['theta_ds'], o['f'])[0,1]:.3f}")
    print(f"corr(θ_ΔI, known f):  r={np.corrcoef(o['fs'], o['f'])[0,1]:.3f}")
    diff = o["theta_ds"] - o["fs"]
    print(f"θ_ΔS − θ_ΔI: median={diff.median():+.3f}  "
          f"MAD={np.median(np.abs(diff - diff.median())):.3f}")
    print("\n  f     median θ_ΔS   median θ_ΔI   (n)")
    for f, g in o.groupby("f"):
        print(f"{f:5.3f}     {g['theta_ds'].median():6.3f}       "
              f"{g['fs'].median():6.3f}    ({len(g)})")


if __name__ == "__main__":
    main()
