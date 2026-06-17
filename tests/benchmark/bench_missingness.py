"""bench_missingness.py — per-precursor missingness across a labeling time series.

Motivation (PROJECT_REVIEW Track A — "Re-explore match-between-runs for the
mzTab/DDA path"). DDA picks precursors for MS2 stochastically, so a peptide can
drop out of individual runs and leave holes in its turnover curve; DIA-NN, by
contrast, propagates IDs across runs internally (built-in MBR). Before building
any MBR for the mzTab/DDA path, the roadmap says to *measure the gap first*:

  - How incomplete are DDA turnover curves (per precursor, across timepoints)?
  - How much of that missingness is **recoverable** by MBR — i.e. the precursor
    is observed in *some* runs (a donor exists) but missing in others — vs
    genuinely absent everywhere (MBR cannot help)?
  - Does DIA's internal propagation already make it complete enough to *not*
    need MBR, confirming DDA is where MBR would pay?

It reads `riana integrate` ``*_riana.txt`` outputs (one per run), builds a
precursor x run presence matrix per curve, and reports completeness +
recoverable-gap metrics. No fit step is required.

Precursor identity is the ``concat`` column (``SEQUENCE_charge``) — the
precursor-level key MBR would transfer (integration needs a specific m/z, so
charge states are distinct). A sequence-level (charge-collapsed) view is also
reported because curve-fitting rolls charge states up per peptide.

Usage:
    python tests/benchmark/bench_missingness.py            # all known sets
    python tests/benchmark/bench_missingness.py --csv out.csv
"""

from __future__ import annotations

import argparse
import glob
import itertools
import os
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@dataclass
class RunSet:
    """One turnover curve: an ordered list of (slot_label, [run files])."""

    name: str
    acquisition: str  # "DDA" | "DIA"
    # ordered list of (timepoint_label, list_of_run_paths); >1 path == replicates
    slots: list[tuple[str, list[str]]]


def _read_precursors(path: str) -> set[str]:
    """Distinct precursor ids (``concat``) present in one ``_riana.txt`` run."""
    df = pd.read_csv(path, sep="\t", comment="#", usecols=lambda c: c in ("concat",))
    return set(df["concat"].dropna().astype(str))


def _read_sequences(path: str) -> set[str]:
    df = pd.read_csv(path, sep="\t", comment="#", usecols=lambda c: c in ("sequence",))
    return set(df["sequence"].dropna().astype(str))


def _build_dda_set(name: str, run_dir: str, chamber: str) -> RunSet:
    """DDA: one acquisition per timepoint, ordered by the ``timeNN`` token."""
    files = sorted(glob.glob(os.path.join(run_dir, f"*_{chamber}_time*_riana.txt")))
    slots: list[tuple[str, list[str]]] = []
    for f in files:
        m = re.search(r"_time(\d+)_riana\.txt$", os.path.basename(f))
        slots.append((f"t{m.group(1)}", [f]))
    return RunSet(name=name, acquisition="DDA", slots=slots)


def _build_dia_set(name: str, run_dir: str) -> RunSet:
    """DIA: group the per-run files by labeling_time from the manifest (reps)."""
    man = pd.read_csv(
        os.path.join(run_dir, "riana_manifest.tsv"), sep="\t", comment="#"
    )
    man = man[man["stage"] == "integrate"].copy()
    man["labeling_time"] = man["labeling_time"].astype(float)
    slots: list[tuple[str, list[str]]] = []
    for day, grp in sorted(man.groupby("labeling_time")):
        paths = [os.path.join(REPO, p) for p in grp["output_path"]]
        slots.append((f"d{day:g}", sorted(paths)))
    return RunSet(name=name, acquisition="DIA", slots=slots)


def _presence_matrix(runset: RunSet, level: str) -> tuple[pd.DataFrame, list[str]]:
    """Boolean precursor x slot matrix (slot present if seen in >=1 of its reps)."""
    reader = _read_precursors if level == "precursor" else _read_sequences
    slot_sets: dict[str, set[str]] = {}
    for label, paths in runset.slots:
        ids: set[str] = set()
        for p in paths:
            ids |= reader(p)
        slot_sets[label] = ids
    labels = [lbl for lbl, _ in runset.slots]
    universe = sorted(set().union(*slot_sets.values()))
    mat = pd.DataFrame(
        {lbl: [pid in slot_sets[lbl] for pid in universe] for lbl in labels},
        index=universe,
    )
    return mat, labels


def analyze(runset: RunSet, level: str = "precursor") -> dict:
    mat, labels = _presence_matrix(runset, level)
    n_slots = len(labels)
    seen = mat.sum(axis=1).to_numpy()  # how many slots each id is present in
    n_ids = len(seen)

    # completeness
    complete = int((seen == n_slots).sum())
    singletons = int((seen == 1).sum())
    hist = {k: int((seen == k).sum()) for k in range(1, n_slots + 1)}

    # recoverable gap: ids present in >=2 slots have a donor for MBR.
    has_donor = seen >= 2
    total_slots_donor = int(has_donor.sum()) * n_slots
    filled_slots_donor = int(seen[has_donor].sum())
    recoverable_gap = (
        (total_slots_donor - filled_slots_donor) / total_slots_donor
        if total_slots_donor
        else 0.0
    )

    # overall empty-slot rate across ALL observed ids (upper bound on holes)
    total_slots = n_ids * n_slots
    overall_gap = (total_slots - int(seen.sum())) / total_slots

    # anchor (t0) loss — DDA curves need the first timepoint as the m0 baseline.
    t0_present = bool(mat[labels[0]].any())
    # of ids seen at >=1 *later* slot, how many miss t0?
    later_any = mat[labels[1:]].any(axis=1).to_numpy()
    miss_t0 = int(((~mat[labels[0]].to_numpy()) & later_any).sum())
    miss_t0_frac = miss_t0 / max(int(later_any.sum()), 1)

    # "fittable" yield at coverage thresholds (fraction of timepoints present)
    thresholds = sorted({2, 3, n_slots // 2, max(n_slots - 2, 2), n_slots})
    fittable = {f">={k}/{n_slots}": int((seen >= k).sum()) for k in thresholds}

    return dict(
        name=runset.name,
        acquisition=runset.acquisition,
        level=level,
        n_slots=n_slots,
        n_ids=n_ids,
        complete=complete,
        complete_frac=complete / n_ids,
        singletons=singletons,
        singleton_frac=singletons / n_ids,
        median_seen=float(np.median(seen)),
        mean_seen=float(seen.mean()),
        overall_gap=overall_gap,
        recoverable_gap=recoverable_gap,
        n_with_donor=int(has_donor.sum()),
        miss_t0_frac=miss_t0_frac,
        hist=hist,
        fittable=fittable,
        t0_present=t0_present,
    )


def _slot_completeness(slot_sets: list[set[str]]) -> float:
    """Fraction of the union present in *all* slots (the DIA completeness metric)."""
    union = set().union(*slot_sets)
    if not union:
        return 0.0
    full = set.intersection(*slot_sets)
    return len(full) / len(union)


def control_three_timepoint(dda_sets: list[RunSet], dia: RunSet) -> None:
    """Geometry control: DIA has only 3 timepoints; DDA has 12. Completeness over
    3 slots is far easier than over 12, so re-measure DDA over 3 timepoints to
    isolate per-acquisition stochastic loss from the curve-length effect.

    For each DDA chamber, completeness is averaged over *every* 3-timepoint
    combination drawn from the non-t0 timepoints (excludes the special t0 anchor,
    per the request). DIA is shown as-is (3 tp x 3 reps) and at 1 rep/timepoint
    (matches DDA's one acquisition per slot) so geometry *and* shots-per-slot are
    controlled.
    """
    print("\n" + "=" * 78)
    print("GEOMETRY CONTROL — completeness over 3 timepoints (DIA-matched)")
    print("=" * 78)

    for rs in dda_sets:
        # per-timepoint precursor sets, excluding t0 (slot 0)
        sets = [set().union(*[_read_precursors(p) for p in paths]) for _, paths in rs.slots]
        labels = [lbl for lbl, _ in rs.slots]
        non_t0 = list(range(1, len(sets)))
        comps = [
            _slot_completeness([sets[i] for i in combo])
            for combo in itertools.combinations(non_t0, 3)
        ]
        comps = np.array(comps)
        print(
            f"  {rs.name:<12} 3-of-{len(non_t0)} non-t0 tp  "
            f"({len(comps)} combos): complete = "
            f"mean {_fmt_pct(comps.mean())}  median {_fmt_pct(np.median(comps))}  "
            f"range [{_fmt_pct(comps.min())}, {_fmt_pct(comps.max())}]"
        )

    # DIA as-is (3 reps) and at 1 rep/timepoint (first sorted run per slot = biorep 1)
    dia_3rep = [set().union(*[_read_precursors(p) for p in paths]) for _, paths in dia.slots]
    dia_1rep = [_read_precursors(paths[0]) for _, paths in dia.slots]
    print(
        f"  {'LV (DIA)':<12} 3 tp x 3 reps        : complete = "
        f"{_fmt_pct(_slot_completeness(dia_3rep))}"
    )
    print(
        f"  {'LV (DIA)':<12} 3 tp x 1 rep (biorep1): complete = "
        f"{_fmt_pct(_slot_completeness(dia_1rep))}   <- matches DDA's 1 acq/slot"
    )
    print(
        "\ninterpretation: if DDA-3tp completeness is still well below DIA-3tp,\n"
        "the gap is acquisition/ID-pipeline (stochastic MS2 + no propagation),\n"
        "not curve length -> MBR is warranted independent of the timepoint-count confound."
    )


def _fmt_pct(x: float) -> str:
    return f"{100 * x:5.1f}%"


def print_report(results: list[dict]) -> None:
    print("\n" + "=" * 78)
    print("PER-PRECURSOR MISSINGNESS  (identity = concat = SEQUENCE_charge)")
    print("=" * 78)
    hdr = (
        f"{'set':<14}{'acq':<5}{'slots':>6}{'precursors':>12}"
        f"{'complete':>10}{'median':>8}{'overall':>9}{'recover':>9}"
    )
    print(hdr)
    print(
        f"{'':14}{'':5}{'(tp)':>6}{'(union)':>12}{'(all tp)':>10}"
        f"{'#seen':>8}{'gap':>9}{'gap':>9}"
    )
    print("-" * 78)
    for r in results:
        if r["level"] != "precursor":
            continue
        print(
            f"{r['name']:<14}{r['acquisition']:<5}{r['n_slots']:>6}"
            f"{r['n_ids']:>12}{_fmt_pct(r['complete_frac']):>10}"
            f"{r['median_seen']:>8.0f}{_fmt_pct(r['overall_gap']):>9}"
            f"{_fmt_pct(r['recoverable_gap']):>9}"
        )
    print("-" * 78)
    print(
        "complete gap = empty (precursor,run) slots / all slots, over every "
        "observed precursor\nrecover gap  = empty slots among precursors seen "
        "in >=2 runs (a donor exists) -> MBR-addressable"
    )

    for r in results:
        if r["level"] != "precursor":
            continue
        print(f"\n--- {r['name']} ({r['acquisition']}) ---")
        print(
            f"  precursors (union)        : {r['n_ids']:,}\n"
            f"  complete curves (all {r['n_slots']} tp): "
            f"{r['complete']:,} ({_fmt_pct(r['complete_frac'])})\n"
            f"  seen in only 1 tp         : {r['singletons']:,} "
            f"({_fmt_pct(r['singleton_frac'])})  <- no donor, MBR can't help\n"
            f"  precursors w/ a donor (>=2): {r['n_with_donor']:,}\n"
            f"  median timepoints / precursor: {r['median_seen']:.0f} / {r['n_slots']}\n"
            f"  missing the t0 anchor      : {_fmt_pct(r['miss_t0_frac'])} of "
            f"precursors seen at a later tp"
        )
        print("  coverage yield (precursors seen in >= k of N timepoints):")
        for k, v in r["fittable"].items():
            print(f"      {k:<9}: {v:,}")
        print("  completeness histogram (precursors seen in exactly k tp):")
        line = "      "
        for k, v in r["hist"].items():
            line += f"{k}:{v:,}  "
        print(line)

    # sequence-level addendum
    print("\n" + "=" * 78)
    print("SEQUENCE-LEVEL (charge-collapsed) completeness, for reference")
    print("=" * 78)
    for r in results:
        if r["level"] != "sequence":
            continue
        print(
            f"  {r['name']:<12} {r['acquisition']:<4} "
            f"complete {_fmt_pct(r['complete_frac'])}  "
            f"overall gap {_fmt_pct(r['overall_gap'])}  "
            f"recover gap {_fmt_pct(r['recoverable_gap'])}  "
            f"(n={r['n_ids']:,})"
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", help="write the per-set summary rows to this CSV")
    args = ap.parse_args()

    runs = os.path.join(REPO, "runs")
    sets = [
        _build_dda_set("LVE (DDA)", os.path.join(runs, "lve_atr"), "LVE"),
        _build_dda_set("ATR (DDA)", os.path.join(runs, "lve_atr"), "ATR"),
        _build_dia_set("LV (DIA)", os.path.join(runs, "lve_dia")),
    ]

    results: list[dict] = []
    for rs in sets:
        for level in ("precursor", "sequence"):
            results.append(analyze(rs, level))

    print_report(results)
    control_three_timepoint(sets[:2], sets[2])

    if args.csv:
        rows = [
            {k: v for k, v in r.items() if k not in ("hist", "fittable")}
            for r in results
        ]
        pd.DataFrame(rows).to_csv(args.csv, index=False)
        print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
