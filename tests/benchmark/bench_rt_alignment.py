"""bench_rt_alignment.py — is the quantms mzTab RT aligned across runs?

MBR (PROJECT_REVIEW Track A) transfers a donor precursor's identity + retention
time into an acceptor run where it was not picked for MS2. That only works if we
know the acceptor-run RT. Two regimes:

  - **Aligned** — ProteomicsLFQ ran RT alignment and the PSM ``retention_time``
    is on a common frame. Co-identified precursors then sit on the RT identity
    line across runs, and MBR can use the donor RT directly.
  - **Raw** — the PSM RT is the original per-run spectrum RT. Co-identified
    precursors show a smooth off-diagonal warp (slope != 1 / nonzero
    intercept / separation-growing offset), and MBR must *learn* the alignment
    from shared IDs before transferring.

This probe reads the raw quantms mzTab PSM section, keeps target PSMs
(q <= 0.01, non-decoy), takes the median RT per (run, precursor=SEQUENCE_charge),
and for each run pair reports the RT identity-line residual (median |dRT|) plus a
robust line fit (slope/intercept/r) on representative pairs. A compact pairwise
median-|dRT| matrix over one chamber's timepoints shows whether the offset grows
with timepoint separation (the raw signature).

Usage:
    python tests/benchmark/bench_rt_alignment.py            # LVE chamber
    python tests/benchmark/bench_rt_alignment.py --chamber ATR
"""

from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MZTAB = os.path.join(
    REPO,
    "data/timeseries_lve_atr/quantms_results/quant_tables",
    "samplesheet_lve_atr.sdrf_openms_design_openms.mzTab",
)

# 0-based PSM column indices (verified against the PSH header)
C_SEQ, C_RT, C_Z, C_REF, C_Q, C_DECOY = 1, 10, 11, 14, 20, 21
_RUN_RE = re.compile(r"ms_run\[(\d+)\]")


def parse_mztab(path: str, q_max: float = 0.01):
    """Return (run_stem_by_idx, rt[run_idx][precursor] -> median RT seconds)."""
    run_stem: dict[int, str] = {}
    # per (run_idx, precursor) collect RTs, then take the median
    rts: dict[int, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            if line.startswith("MTD"):
                if "-location" in line:
                    parts = line.rstrip("\n").split("\t")
                    m = _RUN_RE.search(parts[1])
                    stem = os.path.basename(parts[2]).replace(".mzML", "")
                    run_stem[int(m.group(1))] = stem
                continue
            if not line.startswith("PSM\t"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[C_DECOY] == "1":
                continue
            try:
                if float(p[C_Q]) > q_max:
                    continue
                rt = float(p[C_RT])
            except (ValueError, IndexError):
                continue
            m = _RUN_RE.search(p[C_REF])
            if not m:
                continue
            run = int(m.group(1))
            prec = f"{p[C_SEQ]}_{p[C_Z]}"
            rts[run][prec].append(rt)

    med: dict[int, dict[str, float]] = {}
    for run, d in rts.items():
        med[run] = {prec: float(np.median(v)) for prec, v in d.items()}
    return run_stem, med


def chamber_runs(run_stem: dict[int, str], chamber: str):
    """Ordered [(timepoint_label, run_idx)] for a chamber, by the timeNN token."""
    out = []
    for idx, stem in run_stem.items():
        if f"_{chamber}_" not in stem:
            continue
        m = re.search(r"_time(\d+)", stem)
        out.append((int(m.group(1)), f"t{m.group(1)}", idx))
    out.sort()
    return [(lbl, idx) for _, lbl, idx in out]


def pair_stats(a: dict[str, float], b: dict[str, float]):
    """RT comparison for precursors co-identified in runs a and b."""
    shared = a.keys() & b.keys()
    if len(shared) < 20:
        return None
    xa = np.array([a[p] for p in shared])
    xb = np.array([b[p] for p in shared])
    d = xb - xa
    # robust line fit b ~ slope*a + intercept (least squares on shared medians)
    slope, intercept = np.polyfit(xa, xb, 1)
    resid_fit = xb - (slope * xa + intercept)
    return dict(
        n=len(shared),
        med_abs_drt=float(np.median(np.abs(d))),
        p90_abs_drt=float(np.percentile(np.abs(d), 90)),
        median_drt=float(np.median(d)),
        slope=float(slope),
        intercept=float(intercept),
        r=float(np.corrcoef(xa, xb)[0, 1]),
        rmse_identity=float(np.sqrt(np.mean(d**2))),
        rmse_fit=float(np.sqrt(np.mean(resid_fit**2))),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chamber", default="LVE", choices=["LVE", "ATR"])
    ap.add_argument("--mztab", default=MZTAB)
    args = ap.parse_args()

    print(f"parsing {os.path.basename(args.mztab)} ...")
    run_stem, med = parse_mztab(args.mztab)
    runs = chamber_runs(run_stem, args.chamber)
    labels = [lbl for lbl, _ in runs]
    idx_by_label = dict(runs)
    print(f"{args.chamber}: {len(runs)} runs  ({', '.join(labels)})")

    # pairwise median |dRT| matrix (seconds)
    print("\n" + "=" * 72)
    print(f"PAIRWISE median |dRT| (s) across {args.chamber} timepoints")
    print("identity-line residual for co-identified precursors")
    print("=" * 72)
    print("      " + "".join(f"{l:>6}" for l in labels))
    for la in labels:
        row = f"{la:>5} "
        for lb in labels:
            if la == lb:
                row += f"{'·':>6}"
                continue
            s = pair_stats(med[idx_by_label[la]], med[idx_by_label[lb]])
            row += f"{s['med_abs_drt']:>6.1f}" if s else f"{'-':>6}"
        print(row)

    # detailed fit on representative pairs (adjacent / mid / far from t0)
    print("\n" + "=" * 72)
    print("REPRESENTATIVE PAIRS vs t0  (slope~1 & intercept~0 & small RMSE => aligned)")
    print("=" * 72)
    t0 = labels[0]
    picks = [labels[1], labels[len(labels) // 2], labels[-1]]
    hdr = f"{'pair':<14}{'n':>7}{'med|dRT|':>10}{'p90':>8}{'slope':>8}{'intcpt':>9}{'r':>8}{'rmse_id':>9}{'rmse_fit':>9}"
    print(hdr)
    for lb in picks:
        s = pair_stats(med[idx_by_label[t0]], med[idx_by_label[lb]])
        if not s:
            continue
        print(
            f"{t0}->{lb:<9}{s['n']:>7}{s['med_abs_drt']:>10.1f}{s['p90_abs_drt']:>8.1f}"
            f"{s['slope']:>8.3f}{s['intercept']:>9.2f}{s['r']:>8.4f}"
            f"{s['rmse_identity']:>9.1f}{s['rmse_fit']:>9.1f}"
        )
    print(
        "\nmed|dRT| = median |RT_b - RT_a| on the identity line (no fit).\n"
        "rmse_id vs rmse_fit: if rmse_fit << rmse_id, a linear realignment helps\n"
        "  (=> RTs are NOT already on a common frame)."
    )


if __name__ == "__main__":
    main()
