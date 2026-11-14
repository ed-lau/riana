"""Relabeled manifests for the wls-var / D / E evaluation (WP0).

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md``. Each builder reads the
``integrate`` rows of an existing run's manifest and writes a new manifest whose
rows point at the same ``*_riana.txt`` files with relabeled identities, so
``riana fit --manifest <new> -o <dir>`` refits without re-integrating:

- ``individuals``: attach ``characteristics[individual]`` from an SDRF (for runs
  integrated before the column existed, e.g. ATF6).
- ``t1``: replicate split. At every (experiment, condition, timepoint) the
  biological replicates are dealt at random to pseudo-conditions A and B (true
  Δk = 0). Needs exactly 2 replicates per timepoint.
- ``t2``: permutation within timepoint. Condition labels are shuffled between
  individuals at each timepoint, keeping the per-condition counts; an
  individual's runs (all its tissues) move together. Each timepoint's new groups
  mix the original ones (no timepoint keeps the true contrast whole).
- ``t3``: time split. Timepoints after 0 are ranked; odd ranks go to one
  manifest, even ranks to the other; t = 0 goes to both.

Usage::

    python tests/benchmark/wlsvar_manifests.py individuals SRC OUT --sdrf SDRF
    python tests/benchmark/wlsvar_manifests.py t1 SRC OUT --seed 1
    python tests/benchmark/wlsvar_manifests.py t2 SRC OUT --seed 1
    python tests/benchmark/wlsvar_manifests.py t3 SRC OUT_ODD OUT_EVEN
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import numpy as np

from riana.io.manifest import ManifestRow, _write_manifest, read_manifest
from riana.io.sdrf import read_sdrf


def _relabel(row: ManifestRow, **identity) -> ManifestRow:
    return replace(row, identity=replace(row.identity, **identity))


def individuals(rows: list[ManifestRow], sdrf: Path) -> list[ManifestRow]:
    ind = {r.data_file: r.individual for r in read_sdrf(sdrf).runs}
    missing = {r.identity.data_file for r in rows} - ind.keys()
    if missing:
        raise SystemExit(f"{len(missing)} runs not in the SDRF, e.g. {sorted(missing)[:3]}")
    return [_relabel(r, individual=ind[r.identity.data_file]) for r in rows]


def t1(rows: list[ManifestRow], seed: int) -> list[ManifestRow]:
    rng = np.random.default_rng(seed)
    reps = defaultdict(set)
    for r in rows:
        i = r.identity
        reps[(i.experiment, i.condition, i.labeling_time)].add(i.biological_replicate)
    deal = {}
    for key, rs in reps.items():
        if len(rs) != 2:
            raise SystemExit(f"T1 needs 2 replicates per timepoint; {key} has {sorted(rs)}")
        a, b = rng.permutation(sorted(rs))
        deal[key + (a,)], deal[key + (b,)] = "A", "B"
    out = []
    for r in rows:
        i = r.identity
        out.append(_relabel(r, condition=deal[(i.experiment, i.condition, i.labeling_time,
                                               i.biological_replicate)],
                            biological_replicate=1))
    return out


def t2(rows: list[ManifestRow], seed: int, mixed: bool = True) -> list[ManifestRow]:
    """Shuffle condition labels between individuals at each timepoint.

    ``mixed``: redraw a timepoint's shuffle until every new condition holds at least one
    individual of every original condition, so no timepoint carries the true contrast
    intact (with 3 + 3 animals, 18 of the 20 splits qualify).
    """
    rng = np.random.default_rng(seed)
    if any(not r.identity.individual for r in rows):
        raise SystemExit("T2 needs an individual on every run (build `individuals` first)")
    animal = {}  # (labeling_time, individual) -> condition
    for r in rows:
        key = (r.identity.labeling_time, r.identity.individual)
        if animal.setdefault(key, r.identity.condition) != r.identity.condition:
            raise SystemExit(f"individual {key} carries two conditions")
    new = {}
    for t in {k[0] for k in animal}:
        keys = sorted(k for k in animal if k[0] == t)
        orig = [animal[k] for k in keys]
        for _ in range(1000):
            perm = rng.permutation(orig)
            if not mixed or all({o for o, c in zip(orig, perm) if c == cc} == set(orig)
                                for cc in set(orig)):
                break
        else:
            raise SystemExit(f"t = {t}: no mixed shuffle of {orig}")
        for k, c in zip(keys, perm):
            new[k] = str(c)
    # biological_replicate must stay unique per (condition, timepoint): renumber
    reps = defaultdict(dict)
    for (t, ind), c in sorted(new.items()):
        reps[(c, t)].setdefault(ind, len(reps[(c, t)]) + 1)
    out = []
    for r in rows:
        k = (r.identity.labeling_time, r.identity.individual)
        out.append(_relabel(r, condition=new[k], biological_replicate=reps[(new[k], k[0])][k[1]]))
    return out


def t3(rows: list[ManifestRow]) -> tuple[list[ManifestRow], list[ManifestRow]]:
    times = sorted({r.identity.labeling_time for r in rows
                    if r.identity.labeling_time and r.identity.labeling_time > 0})
    odd = set(times[0::2])
    zero = [r for r in rows if not r.identity.labeling_time]
    pos = [r for r in rows if r.identity.labeling_time and r.identity.labeling_time > 0]
    return (zero + [r for r in pos if r.identity.labeling_time in odd],
            zero + [r for r in pos if r.identity.labeling_time not in odd])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("test", choices=["individuals", "t1", "t2", "t3"])
    p.add_argument("src", type=Path, help="source riana_manifest.tsv")
    p.add_argument("out", type=Path, nargs="+", help="output manifest(s); t3 takes two")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--sdrf", type=Path)
    a = p.parse_args()
    rows = read_manifest(a.src, stage="integrate")
    if a.test == "individuals":
        outs = [individuals(rows, a.sdrf)]
    elif a.test == "t1":
        outs = [t1(rows, a.seed)]
    elif a.test == "t2":
        outs = [t2(rows, a.seed)]
    else:
        outs = list(t3(rows))
    if len(outs) != len(a.out):
        raise SystemExit(f"{a.test} writes {len(outs)} manifest(s); got {len(a.out)} paths")
    for path, rs in zip(a.out, outs):
        _write_manifest(path, rs)
        conds = sorted({r.identity.condition for r in rs})
        print(f"{path}: {len(rs)} runs, conditions {conds}")


if __name__ == "__main__":
    main()
