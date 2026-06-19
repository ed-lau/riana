"""bench_mbr_ab.py — does MBR help the fit? include-MBR vs --exclude-mbr A/B.

Compares two ``riana_fit_peptides.txt`` from the *same* integrate run — one fit
with MBR (default), one with ``--exclude-mbr``. The exclude arm is real-only, so
it is identical across MBR gates; only the include arm changes when you tighten
``--mbr-min-snr`` / ``--mbr-min-scans`` at integrate time. Produce the two files:

    riana fit --manifest <run>/riana_manifest.tsv --coefficients <c>              # include
    riana fit --manifest <run>/riana_manifest.tsv --coefficients <c> --exclude-mbr

then::

    python tests/benchmark/bench_mbr_ab.py --include <incl>.txt --exclude <excl>.txt

Readouts: peptide yield, R²>0.95 count (the quality gate), k_deg stability, the
quality of curves that fit *only* because of MBR, and — for curves fit both ways —
whether appending MBR points to an existing curve **degrades** its R² (off-curve
transfers) or not. Keyed on ``(concat, condition)`` (a concat recurs across
conditions).
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

KEY = ["concat", "condition"]


def _converged(path: str) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", comment="#")
    return d[d["R_squared"].notna()].set_index(KEY)


def _iqr(s: pd.Series) -> float:
    return float(s.quantile(0.75) - s.quantile(0.25))


def ab(include_path: str, exclude_path: str) -> None:
    ci = _converged(include_path)   # real + MBR
    ce = _converged(exclude_path)   # real only

    print("=" * 72)
    print("MBR A/B — include (real+MBR) vs exclude (real only)")
    print("=" * 72)
    n_i, n_e = len(ci), len(ce)
    r_i = int((ci["R_squared"] > 0.95).sum())
    r_e = int((ce["R_squared"] > 0.95).sum())
    print(f"  {'metric':<22}{'exclude':>12}{'include':>12}{'Δ':>12}")
    print(f"  {'converged peptides':<22}{n_e:>12,}{n_i:>12,}{n_i - n_e:>+12,}")
    print(f"  {'R²>0.95 (quality)':<22}{r_e:>12,}{r_i:>12,}{r_i - r_e:>+12,}")
    print(f"  {'k_deg median':<22}{ce.k_deg.median():>12.4f}"
          f"{ci.k_deg.median():>12.4f}{ci.k_deg.median() - ce.k_deg.median():>+12.4f}")
    print(f"  {'k_deg IQR':<22}{_iqr(ce.k_deg):>12.4f}{_iqr(ci.k_deg):>12.4f}"
          f"{_iqr(ci.k_deg) - _iqr(ce.k_deg):>+12.4f}")

    # R²-cutoff sweep — does MBR gain net peptides at a looser (in-vivo) gate?
    sweep = []
    for c in (0.95, 0.9, 0.8, 0.7, 0.0):
        e = int((ce["R_squared"] > c).sum())
        i = int((ci["R_squared"] > c).sum())
        sweep.append(f">{c}: {i - e:+,}")
    print("  R²-cutoff net (include−exclude): " + "  ".join(sweep))

    # Curves that converge ONLY with MBR (MBR filled them to depth).
    enabled = ci.loc[ci.index.difference(ce.index)]
    if len(enabled):
        print(f"\n  MBR-ENABLED curves (fit only with MBR): {len(enabled):,}")
        print(f"    R²>0.95: {int((enabled.R_squared > 0.95).sum()):,} "
              f"({100 * (enabled.R_squared > 0.95).mean():.0f}%)  "
              f"median R² {enabled.R_squared.median():.3f}  "
              f"(all-include median {ci.R_squared.median():.3f})")

    # Curves fit BOTH ways that gained MBR points — does MBR degrade them?
    shared = ci.index.intersection(ce.index)
    si, se = ci.loc[shared], ce.loc[shared]
    gained = si["n_mbr"] > 0
    dR2 = (si["R_squared"] - se["R_squared"])[gained]
    dk = (si["k_deg"] - se["k_deg"])[gained]
    if len(dR2):
        worse = int((dR2 < -0.01).sum())
        better = int((dR2 > 0.01).sum())
        same = int((dR2.abs() <= 0.01).sum())
        print(f"\n  SHARED curves that gained ≥1 MBR point: {int(gained.sum()):,}")
        # median, not mean: R² is unbounded below for a broken fit, so the mean is
        # outlier-dominated and meaningless here.
        print(f"    ΔR² (incl−excl) median {dR2.median():+.4f}  "
              f"|Δk| median {dk.abs().median():.4f}")
        print(f"    R² worse {worse:,} / better {better:,} / ~same {same:,}  "
              f"→ {'DEGRADES' if worse > 2 * better else 'neutral/helps'}")
        # Dose: ΔR² scales with how much MBR DOMINATES the curve, not raw count —
        # the case for a fill-gaps-only / cap-fraction policy (low dose is ~free).
        dose = pd.DataFrame({
            "dR2": dR2.to_numpy(),
            "frac": (si["n_mbr"] / si["n_points"])[gained].to_numpy(),
        })
        dose["fb"] = pd.cut(dose["frac"], [0, 0.1, 0.25, 0.5, 1.01],
                            labels=["<10%", "10-25%", "25-50%", ">50%"])
        cells = [f"{b}: {g.dR2.median():+.3f}(n={len(g):,})"
                 for b, g in dose.groupby("fb", observed=True)]
        print("    ΔR² by MBR-fraction of curve: " + "  ".join(cells))

    # Verdict weighs the real harm (pollution of clean curves) and the realistic
    # in-vivo gate (R²>0.8); a noise-level change at the strict 0.95 gate is not a
    # "hurt". (The strict gate rarely moves because MBR rescues mid-quality curves.)
    e80 = int((ce["R_squared"] > 0.8).sum())
    i80 = int((ci["R_squared"] > 0.8).sum())
    net95, net80 = r_i - r_e, i80 - e80
    pollutes = len(dR2) > 0 and float(dR2.median()) < -0.02
    noise95 = abs(net95) <= max(30, int(0.02 * r_e))
    print("\n  VERDICT:", end=" ")
    if pollutes:
        print(f"MBR HURTS — pollutes clean curves (shared ΔR² median "
              f"{float(dR2.median()):+.3f}); tighten the gate.")
    elif net80 > 0 and (net95 >= 0 or noise95):
        print(f"MBR net-POSITIVE at R²>0.8 ({net80:+,}), ~neutral at strict 0.95 "
              f"({net95:+,}), no clean-curve pollution — keep.")
    elif noise95 and abs(net80) <= max(30, int(0.02 * e80)):
        print(f"MBR ~NEUTRAL ({net95:+,} at 0.95, {net80:+,} at 0.8) — adds peptides "
              "at no quality cost.")
    else:
        print(f"MBR net {net95:+,} at R²>0.95 / {net80:+,} at R²>0.8.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--include", required=True, help="include-MBR riana_fit_peptides.txt")
    ap.add_argument("--exclude", required=True, help="--exclude-mbr riana_fit_peptides.txt")
    args = ap.parse_args()
    ab(args.include, args.exclude)


if __name__ == "__main__":
    main()
