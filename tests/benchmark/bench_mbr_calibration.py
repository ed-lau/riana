"""bench_mbr_calibration.py — ground-truth θ recovery for MBR on the D₂O mixing series.

The calibration series mixes fully-labelled (≥10 doublings in 6% D₂O) and unlabelled
lysate at *known* proportions f ∈ {0, 0.125, …, 1.0}. ``solve_fs_d2o`` should recover
θ = f, so for an MBR-transferred peak **|θ − f| is a TRUE ground-truth accuracy
readout** — unlike the LVE within-protein corridor (a proxy). MBR transfers a peptide
across proportions (same peptide, same RT; the mix changes the envelope, not the
elution), so a transferred point at proportion f should still recover θ ≈ f.

Run ``integrate --mbr`` on the calibration mzTab first (the gate defaults apply)::

    riana integrate data/calibration_ac16/mzml <…>.mzTab \\
        --sdrf data/calibration_ac16/samplesheet_ac16.sdrf.tsv --mbr -o runs/cal_ac16_mbr
    python tests/benchmark/bench_mbr_calibration.py --run runs/cal_ac16_mbr --coefficients ac16

Reports, per proportion and overall, the θ-recovery error |θ − f| and bias (θ − f)
for MBR transfers vs direct (q-value) IDs — the direct test of whether transferred
peaks land on the right θ.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

from riana.algorithms.isotope_dist import clear_envelope_cache, solve_fs_d2o
from riana.core.fitting import load_aa_coefficients, spep_from_coefficients
from riana.io.manifest import read_manifest

_ISO = [f"iso{i}" for i in range(6)]
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _theta_for_run(path: str, coeffs: dict, ria: float, q_value: float) -> pd.DataFrame:
    """θ per surviving row of one calibration run's `_riana.txt`."""
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df[df["percolator q-value"] < q_value].copy()
    if df.empty:
        return pd.DataFrame(columns=["evidence", "theta"])
    seqs = df["sequence"].astype(str).to_numpy()
    masses = df["peptide mass"].to_numpy(dtype=float)
    iso = df[_ISO].to_numpy(dtype=float)
    spep_cache: dict[str, int] = {}
    theta = np.full(len(df), np.nan)
    for i in range(len(df)):
        s = seqs[i]
        sp = spep_cache.get(s)
        if sp is None:
            sp = max(1, int(round(spep_from_coefficients(s, coeffs))))
            spep_cache[s] = sp
        try:
            theta[i] = solve_fs_d2o(s, masses[i], iso[i], sp, ria_max=ria, n_iso=6)
        except (KeyError, ValueError):
            pass
    return pd.DataFrame({"evidence": df["evidence"].to_numpy(), "theta": theta})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True, help="integrate --mbr output dir (manifest)")
    ap.add_argument("--coefficients", default="alamillo_2025_ac16", help="per-AA Spep table [alamillo_2025_ac16]")
    ap.add_argument("--ria", type=float, default=0.06, help="D₂O enrichment [0.06]")
    ap.add_argument("--q-value", type=float, default=1e-2)
    args = ap.parse_args()

    coeffs = load_aa_coefficients(args.coefficients)
    rows = read_manifest(os.path.join(args.run, "riana_manifest.tsv"), stage="integrate")
    clear_envelope_cache()

    parts = []
    for r in rows:
        f = r.identity.mixing_proportion
        if f is None:
            continue
        path = r.output_path if os.path.isabs(r.output_path) else os.path.join(_REPO, r.output_path)
        t = _theta_for_run(path, coeffs, args.ria, args.q_value)
        t["f"] = float(f)
        parts.append(t)
    df = pd.concat(parts, ignore_index=True).dropna(subset=["theta"])

    print("=" * 74)
    print("MBR θ-RECOVERY vs ground truth (θ should = mixing proportion f)")
    print("=" * 74)
    print(f"  {'f':>6}{'n_real':>9}{'n_mbr':>8}"
          f"{'|θ−f| real':>12}{'|θ−f| mbr':>11}{'bias real':>11}{'bias mbr':>10}")
    for f, g in df.groupby("f"):
        real = g[g.evidence == "q_value"]; mbr = g[g.evidence == "mbr"]
        def med_abs(x):
            return f"{np.median(np.abs(x.theta - f)):.3f}" if len(x) else "  -"
        def bias(x):
            return f"{np.median(x.theta - f):+.3f}" if len(x) else "  -"
        print(f"  {f:>6.3f}{len(real):>9,}{len(mbr):>8,}"
              f"{med_abs(real):>12}{med_abs(mbr):>11}{bias(real):>11}{bias(mbr):>10}")
    # overall (exclude f=0 and f=1 endpoints from the |θ−f| summary? keep all)
    real = df[df.evidence == "q_value"]; mbr = df[df.evidence == "mbr"]
    print("-" * 74)
    print(f"  overall |θ−f| median: real {np.median(np.abs(real.theta - real.f)):.3f}"
          f"  mbr {np.median(np.abs(mbr.theta - mbr.f)):.3f}   "
          f"(n real {len(real):,} / mbr {len(mbr):,})")
    print("  Read: mbr |θ−f| ≈ real ⇒ transferred peaks recover the right θ; "
          "mbr ≫ real ⇒ MBR mis-quantifies.")


if __name__ == "__main__":
    main()
