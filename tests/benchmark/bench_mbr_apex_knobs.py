"""bench_mbr_apex_knobs.py — sweep apex/extraction knobs on the ac16 calibration MBR.

FINDING (2026-06-20, see reports/2026-06-17_mbr_v1_design.md § Update 2026-06-20b):
this sweep was INCONCLUSIVE / chased the wrong mechanism. MBR θ ≈ 0 for every
proportion and apex_selection nearest = tallest, because the real cause is an
mzTab↔mzML **RT-axis mismatch** (the calibration mzTab is .raw-searched → 2–4 min off
the local .mzML; MBR is RT-anchored, direct is scan-based hence immune), NOT an apex
mis-pick. Re-run only AFTER the maintainer re-searches quantms on the exact .mzML (then
drop --no-rt-check). Kept as the auditable negative result + it exposed the consensus+MBR
gate bug (consensus arms returned n_mbr=0: consensus_apex emits a spread, not an SNR).

MBR mis-quantifies θ on the calibration (a high-label D₂O *mixing* series) but not on
the in-vivo Track D series. Hypothesis: the apex search keys on **iso0**, which D₂O
suppresses, so for a transferred (MBR) row it can latch onto the unlabelled-dominated
peak and bias θ low. This sweep re-integrates the ac16 calibration with ``--mbr`` under
a knob grid and measures the ground-truth error |θ − f| (θ should equal the mixing
proportion f) for **MBR transfers vs direct IDs** — ranking the candidate fixes
empirically rather than by argument.

θ uses the **ac16** per-AA coefficients (the calibration cell line), NOT the in-vivo
(Commerford) set used for Track D, and the ac16 enrichment ria = 0.0598.

Usage: python tests/benchmark/bench_mbr_apex_knobs.py [--force] [--only NAME]
Re-uses each arm's integrate output if present (pass --force to re-integrate).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_mbr_calibration import _theta_for_run  # noqa: E402

from riana.algorithms.isotope_dist import clear_envelope_cache  # noqa: E402
from riana.core.fitting import load_aa_coefficients  # noqa: E402
from riana.io.manifest import read_manifest  # noqa: E402

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_MZML = "data/calibration_ac16/mzml"
_MZTAB = ("data/calibration_ac16/quantms_results/quant_tables/"
          "samplesheet_ac16_alpine.sdrf_openms_design_openms.mzTab")
_SDRF = "data/calibration_ac16/samplesheet_ac16.sdrf.tsv"
_RIA = 0.0598
_COEFFS = "ac16"
_OUT = "runs/cal_ac16_knobs"

# (name, integrate flags). First arm is the shipped default = baseline.
_ARMS = [
    ("apex/tallest hw.25 (default)",
     ["--peak-rt", "apex", "--apex-selection", "tallest", "--apex-search-half-width", "0.25"]),
    ("consensus/tallest",
     ["--peak-rt", "consensus", "--apex-selection", "tallest", "--apex-search-half-width", "0.25"]),
    ("apex/nearest",
     ["--peak-rt", "apex", "--apex-selection", "nearest", "--apex-search-half-width", "0.25"]),
    ("consensus/nearest",
     ["--peak-rt", "consensus", "--apex-selection", "nearest", "--apex-search-half-width", "0.25"]),
    ("apex/tallest hw.10 (tight)",
     ["--peak-rt", "apex", "--apex-selection", "tallest", "--apex-search-half-width", "0.10"]),
    ("apex/tallest 5ppm",
     ["--peak-rt", "apex", "--apex-selection", "tallest", "--apex-search-half-width", "0.25",
      "--mass_tol", "5"]),
]


def _slug(name: str) -> str:
    out = name.replace("/", "-")
    for ch in " ().":
        out = out.replace(ch, "_")
    return out.strip("_").replace("__", "_")


def _integrate(armdir: str, flags: list[str], force: bool) -> None:
    manifest = os.path.join(_REPO, armdir, "riana_manifest.tsv")
    if os.path.exists(manifest) and not force:
        print(f"  [cached] {armdir}")
        return
    # --no-rt-check: the calibration mzTab was searched on the .raw (OpenMS-aligned
    # RT), so the scan↔RT guard trips (up to ~3.9 min) even though the scans truly
    # correspond — confirmed benign earlier by 1.36 ppm median mass accuracy.
    # Gate OFF (--mbr-min-snr/scans 0): isolate the apex *method* effect on θ from
    # the quality gate — and the gate's SNR floor would otherwise drop EVERY
    # consensus row (consensus_apex yields a spread, not an SNR → apex_snr NaN).
    cmd = [sys.executable, "-m", "riana.cli", "integrate", _MZML, _MZTAB,
           "--sdrf", _SDRF, "--mbr", "--no-rt-check",
           "--mbr-min-snr", "0", "--mbr-min-scans", "0", "-o", armdir] + flags
    print(f"  integrating → {armdir} ...", flush=True)
    r = subprocess.run(cmd, cwd=_REPO, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])
        raise SystemExit(f"integrate FAILED for {armdir}")


def _theta_df(armdir: str, coeffs: dict) -> pd.DataFrame:
    rows = read_manifest(os.path.join(_REPO, armdir, "riana_manifest.tsv"), stage="integrate")
    clear_envelope_cache()
    parts = []
    for r in rows:
        f = r.identity.mixing_proportion
        if f is None:
            continue
        path = r.output_path if os.path.isabs(r.output_path) else os.path.join(_REPO, r.output_path)
        t = _theta_for_run(path, coeffs, _RIA, 1e-2)
        t["f"] = float(f)
        parts.append(t)
    return pd.concat(parts, ignore_index=True).dropna(subset=["theta"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-integrate even if cached")
    ap.add_argument("--only", default=None, help="run a single arm by name substring")
    args = ap.parse_args()

    coeffs = load_aa_coefficients(_COEFFS)
    arms = [(n, f) for n, f in _ARMS if args.only is None or args.only in n]

    print("Integrating arms (ac16 calibration, --mbr):")
    results = {}
    for name, flags in arms:
        armdir = os.path.join(_OUT, _slug(name))
        _integrate(armdir, flags, args.force)
        results[name] = _theta_df(armdir, coeffs)

    base = arms[0][0]
    base_mbr = results[base][results[base].evidence == "mbr"]
    base_abs = float(np.median(np.abs(base_mbr.theta - base_mbr.f))) if len(base_mbr) else float("nan")

    print("\n" + "=" * 86)
    print("θ-RECOVERY |θ−f| BY ARM  (ac16 coeffs, ria 0.0598; θ should = f)")
    print("=" * 86)
    print(f"  {'arm':<30}{'n_mbr':>7}{'mbr|θ−f|':>10}{'real|θ−f|':>11}"
          f"{'mbr bias':>10}{'Δ vs base':>11}")
    for name, _ in arms:
        df = results[name]
        mbr = df[df.evidence == "mbr"]
        real = df[df.evidence == "q_value"]
        m_abs = float(np.median(np.abs(mbr.theta - mbr.f))) if len(mbr) else float("nan")
        r_abs = float(np.median(np.abs(real.theta - real.f))) if len(real) else float("nan")
        m_bias = float(np.median(mbr.theta - mbr.f)) if len(mbr) else float("nan")
        delta = m_abs - base_abs
        print(f"  {name:<30}{len(mbr):>7,}{m_abs:>10.3f}{r_abs:>11.3f}"
              f"{m_bias:>+10.3f}{delta:>+11.3f}")

    # Per-proportion MBR |θ−f|, arms as columns — the high-f cells are where it breaks.
    fs = sorted(results[base]["f"].unique())
    print("\nMBR |θ−f| per proportion f (rows) × arm (cols):")
    header = "  " + f"{'f':>6}" + "".join(f"{_slug(n):>20}" for n, _ in arms)
    print(header)
    for f in fs:
        cells = []
        for name, _ in arms:
            g = results[name]
            g = g[(g.evidence == "mbr") & (g.f == f)]
            cells.append(f"{np.median(np.abs(g.theta - f)):>20.3f}" if len(g) else f"{'-':>20}")
        print(f"  {f:>6.3f}" + "".join(cells))
    print("\nRead: lower mbr|θ−f| = better θ recovery; compare each arm's mbr to its own real.")


if __name__ == "__main__":
    main()
