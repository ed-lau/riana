"""--depth semantics spike: rows vs distinct-timepoints qualification.

The fit's depth gate (``core/fitting.py``) counts **raw PSM rows** per
``(experiment, condition, fit_key)`` — `len(x) >= depth` on the modern/manifest
path — with no recombine before it. Because a peptidoform can have several PSM
rows in a *single* run (up to ~9 in LVE/ATR), ``depth=3`` on rows can admit a
curve seen at a single labeling timepoint, which a one-exponent kinetic fit
cannot identify. This bench quantifies how many qualified curves would change
under the proposed harmonization to **distinct labeling timepoints**
(``x[time].nunique() >= depth``), which is <= rows always (it can only reject).

Usage: python bench_depth_semantics.py [manifest.tsv]
Reads the integrate frames the manifest points at and replicates the fit's
assembly (concat per group, q-value filter, fit_key grouping) — no fitting.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from riana.core.fitting import _fit_key

Q_VALUE = 0.01
DEPTHS = (3, 4, 5)


def _load_group(rows: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for r in rows.itertuples():
        df = pd.read_csv(r.output_path, sep="\t", comment="#")
        keep = ["concat", "percolator q-value"]
        if "evidence" in df.columns:
            keep.append("evidence")
        df = df[keep].copy()
        df["labeling_time"] = float(r.labeling_time)
        frames.append(df)
    rdf = pd.concat(frames, ignore_index=True)
    rdf = rdf[rdf["percolator q-value"] < Q_VALUE].copy()
    rdf["fit_key"] = rdf["concat"].map(_fit_key)
    return rdf


def _qualify(rdf: pd.DataFrame) -> pd.DataFrame:
    g = rdf.groupby("fit_key")
    return pd.DataFrame({"rows": g.size(), "tp": g["labeling_time"].nunique()})


def main(manifest: Path) -> None:
    m = pd.read_csv(manifest, sep="\t", comment="#")
    m = m[m["stage"] == "integrate"]
    print(f"manifest: {manifest}  ({len(m)} integrate runs, "
          f"{m.groupby(['experiment', 'condition']).ngroups} fit groups)\n")

    allq, allq_nombr = [], []
    for _, g in m.groupby(["experiment", "condition"]):
        rdf = _load_group(g)
        allq.append(_qualify(rdf))
        if "evidence" in rdf.columns:
            allq_nombr.append(_qualify(rdf[rdf["evidence"] != "mbr"]))

    q = pd.concat(allq, ignore_index=True)
    print(f"curves (fit_key × group): {len(q):,}")
    print(f"  PSM-rows per curve: median {q.rows.median():.0f}, "
          f"p90 {q.rows.quantile(0.9):.0f}, max {q.rows.max():.0f}")
    print(f"  distinct timepoints per curve: median {q.tp.median():.0f}, "
          f"max {q.tp.max():.0f}")
    print(f"  rows == distinct-tp for all? {(q.rows == q.tp).all()}\n")

    print(f"{'depth':>6} {'qual(rows)':>11} {'qual(tp)':>10} {'flip→fail':>10} "
          f"{'flip %':>7}")
    for d in DEPTHS:
        qual_rows = int((q.rows >= d).sum())
        qual_tp = int((q.tp >= d).sum())
        flip = int(((q.rows >= d) & (q.tp < d)).sum())
        print(f"{d:>6} {qual_rows:>11,} {qual_tp:>10,} {flip:>10,} "
              f"{100 * flip / max(qual_rows, 1):>6.1f}%")

    if allq_nombr:
        qn = pd.concat(allq_nombr, ignore_index=True)
        print("\nwithout MBR (clean rows only):")
        for d in DEPTHS:
            qr = int((qn.rows >= d).sum())
            flip = int(((qn.rows >= d) & (qn.tp < d)).sum())
            print(f"  depth={d}: qual(rows) {qr:,}, flip→fail {flip:,} "
                  f"({100 * flip / max(qr, 1):.1f}%)")


if __name__ == "__main__":
    mf = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "runs/lve_atr_mbr_gated/riana_manifest.tsv")
    main(mf)
