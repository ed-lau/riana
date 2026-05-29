"""bench_id_path.py — ID-path concordance for the M3 I/O layer.

Compares two PSM sources for the same calibration mzMLs:

  * Percolator: per-fraction ``percolator.target.psms.txt`` from the snakemake
    Crux/Percolator pipeline at
    ``data/calibration_<line>/snakemake_results/time<N>/percolator/``.
  * mzTab: the combined ``samplesheet_<line>_alpine.sdrf_openms_design_openms.mzTab``
    from the quantms run at ``data/calibration_<line>/quantms_results/quant_tables/``.

Reports run in two phases — the user asked for both PSM-level overlap and
downstream m0/mA agreement (M3 Week 2 plan, PROJECT_REVIEW.md §3 step 2 +
verification list under M3):

Phase 1 — **PSM overlap.** Per nominal proportion, intersect on
``(scan, sequence, charge)`` after filtering each source to q<0.01. Output:
``psm_overlap.csv`` with per-proportion Jaccard, q-value Spearman rho on the
intersection, count breakdowns. Fast, runs in <1 min on a 9-fraction line.

Phase 2 — **Downstream m0/mA agreement.** Project the mzTab PSMs into Crux
Percolator-shaped TSVs (one per fraction) and run ``riana integrate`` against
each, producing a parallel ``integrate_outputs/v0.9.0_mztab/`` tree. Then call
:mod:`bench_m0_ma_recovery` against both trees (with the same frozen
coefficients) and diff the summary JSONs. ~minutes per line because integrate
walks every MS1 in 9 mzMLs.

The two sources are not expected to be identical — quantms applies a combined
target-decoy FDR across all 9 fractions while snakemake searches each
fraction independently. The point of the bench is to quantify the divergence.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from riana.io.mztab import read_mztab
from riana.io.percolator import filter_by_q_value, read_percolator
from riana.io.writers import make_provenance, write_tsv
from riana.records import PSMRecord


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT / "data"
BENCH_DATA = REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing"

# Pinned integrate args — must match run_integrate_v0_9_0.py exactly, or the
# m0/mA comparison is between mismatched configs (PROJECT_REVIEW.md M3
# verification: "must be run with the same -m 15 mass window as the committed
# v0.9.0 baseline").
INTEGRATE_ARGS = [
    "-i", "0", "1", "2", "3", "4", "5",
    "-q", "0.01",
    "-r", "0.33",
    "-m", "15",
    "-t", "4",
]

# The Crux Percolator column order ``ReadPercolator`` expects when it reads
# back a TSV. ``write_mztab_as_percolator_tsv`` writes these exact columns.
PERCOLATOR_COLUMNS = [
    "file_idx", "scan", "charge",
    "spectrum precursor m/z", "spectrum neutral mass", "peptide mass",
    "percolator score", "percolator q-value", "percolator PEP",
    "distinct matches/spectrum",
    "sequence", "protein id", "flanking aa",
]


# ---------------------------------------------------------------------------
# Loading + indexing
# ---------------------------------------------------------------------------


def _proportion_to_sample(prop: float) -> str:
    """0.0 -> 'time0', 12.5 -> 'time12.5'."""
    return f"time{prop:g}"


def _mzml_basename(name: str) -> str:
    """'20230731_AC16_WL_0.mzML.gz' -> '20230731_AC16_WL_0'."""
    stem = Path(name).name
    if stem.endswith(".mzML.gz"):
        return stem[: -len(".mzML.gz")]
    return Path(stem).stem


def load_percolator_by_proportion(
    line: str, q_threshold: float | None = 0.01
) -> dict[float, list[PSMRecord]]:
    """Return ``{proportion: [PSMRecord]}`` from the snakemake Percolator tree.

    ``file_name`` on each record is back-filled from ``ground_truth.csv``
    (Crux Percolator output does not carry the mzML basename, but the bench
    needs it as the join key against mzTab).
    """
    truth = pd.read_csv(BENCH_DATA / line / "ground_truth.csv")
    out: dict[float, list[PSMRecord]] = {}
    for row in truth.itertuples(index=False):
        prop = float(row.nominal_proportion)
        sample = _proportion_to_sample(prop)
        psms_path = (
            DATA_ROOT / f"calibration_{line}" / "snakemake_results" / sample
            / "percolator" / "percolator.target.psms.txt"
        )
        if not psms_path.exists():
            raise FileNotFoundError(f"Percolator file not found: {psms_path}")
        records = read_percolator(psms_path, sample=sample)
        mzml_base = _mzml_basename(row.mzml_filename)
        records = [dataclasses.replace(r, file_name=mzml_base) for r in records]
        if q_threshold is not None:
            records = filter_by_q_value(records, q_threshold)
        out[prop] = records
    return out


def load_mztab_by_proportion(
    line: str, q_threshold: float | None = 0.01
) -> tuple[dict[float, list[PSMRecord]], Path]:
    """Return ``({proportion: [PSMRecord]}, mztab_path)``.

    The single mzTab is split per-proportion by matching each PSM's
    ``file_name`` to ``ground_truth.csv``. Records whose ``file_name`` does
    not appear in the ground-truth map are dropped (they would be acquisitions
    quantms saw but the snakemake side did not).
    """
    truth = pd.read_csv(BENCH_DATA / line / "ground_truth.csv")
    basename_to_prop = {
        _mzml_basename(row.mzml_filename): float(row.nominal_proportion)
        for row in truth.itertuples(index=False)
    }
    quant_dir = DATA_ROOT / f"calibration_{line}" / "quantms_results" / "quant_tables"
    mztab_candidates = sorted(quant_dir.glob("*.mzTab"))
    if not mztab_candidates:
        raise FileNotFoundError(f"no mzTab in {quant_dir}")
    mztab_path = mztab_candidates[0]
    records, _file_map = read_mztab(mztab_path, sample=line)

    by_prop: dict[float, list[PSMRecord]] = defaultdict(list)
    for r in records:
        prop = basename_to_prop.get(r.file_name)
        if prop is None:
            continue
        if q_threshold is not None and r.percolator_q_value >= q_threshold:
            continue
        by_prop[prop].append(r)
    return dict(by_prop), mztab_path


# ---------------------------------------------------------------------------
# Phase 1 — PSM overlap
# ---------------------------------------------------------------------------


def psm_overlap_table(
    percolator: dict[float, list[PSMRecord]],
    mztab: dict[float, list[PSMRecord]],
) -> pd.DataFrame:
    """Per-proportion Jaccard + q-value rank-correlation on the intersection."""
    rows: list[dict] = []
    all_props = sorted(set(percolator) | set(mztab))
    for prop in all_props:
        p_recs = percolator.get(prop, [])
        m_recs = mztab.get(prop, [])
        # Identity key: (scan, peptide-sequence-bare, charge). Both sources
        # carry an unmodified-residue ``sequence`` field, so this is a clean
        # join even though the underlying ID engines differ.
        p_keys = {(r.scan, r.sequence, r.charge): r for r in p_recs}
        m_keys = {(r.scan, r.sequence, r.charge): r for r in m_recs}
        shared = p_keys.keys() & m_keys.keys()
        union = p_keys.keys() | m_keys.keys()

        if shared:
            pq = np.asarray([p_keys[k].percolator_q_value for k in shared])
            mq = np.asarray([m_keys[k].percolator_q_value for k in shared])
            # Pearson on -log10(q+eps) is more interpretable than raw q (which
            # piles at zero); both q sources hit q-min ties so ranks are noisy.
            eps = 1e-12
            with np.errstate(divide="ignore"):
                pq_log = -np.log10(pq + eps)
                mq_log = -np.log10(mq + eps)
            if pq_log.std() > 0 and mq_log.std() > 0:
                q_corr = float(np.corrcoef(pq_log, mq_log)[0, 1])
            else:
                q_corr = float("nan")
        else:
            q_corr = float("nan")

        rows.append({
            "proportion": prop,
            "n_percolator": len(p_keys),
            "n_mztab": len(m_keys),
            "n_intersection": len(shared),
            "n_only_percolator": len(p_keys) - len(shared),
            "n_only_mztab": len(m_keys) - len(shared),
            "n_union": len(union),
            "jaccard": len(shared) / len(union) if union else float("nan"),
            "neg_log_q_corr": q_corr,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Phase 2 — write Crux-shaped percolator TSVs + run integrate
# ---------------------------------------------------------------------------


def write_mztab_as_percolator_tsv(
    records: Iterable[PSMRecord], path: Path
) -> None:
    """Project mzTab-derived records into the Crux Percolator TSV schema.

    ``ReadPercolator`` will read the result and the integrator behaves the
    same as on a real Percolator file. The integrator does not consume
    ``protein id`` or ``flanking aa`` directly, but we keep them populated
    so the file is greppable.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    # file_idx is reset to 0 because integrate runs against a single-mzML
    # tempdir per fraction (matching run_integrate_v0_9_0.py's pattern).
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PERCOLATOR_COLUMNS, delimiter="\t")
        writer.writeheader()
        for r in records:
            writer.writerow({
                "file_idx": 0,
                "scan": r.scan,
                "charge": r.charge,
                "spectrum precursor m/z": r.precursor_mz,
                "spectrum neutral mass": r.neutral_mass,
                "peptide mass": r.peptide_mass,
                "percolator score": r.percolator_score,
                "percolator q-value": r.percolator_q_value,
                "percolator PEP": r.percolator_pep,
                "distinct matches/spectrum": r.distinct_matches,
                "sequence": r.sequence,
                "protein id": r.protein_id,
                "flanking aa": r.flanking_aa or "--",
            })


def run_integrate_from_mztab(
    line: str, *, output_dir: Path, dry_run: bool = False
) -> Path:
    """Run ``riana integrate`` against mzTab-derived PSMs, one tempfile per fraction."""
    mztab_path = next(
        (DATA_ROOT / f"calibration_{line}" / "quantms_results" / "quant_tables")
        .glob("*.mzTab")
    )
    print(f"[mztab] loading {mztab_path}")
    records, _ = read_mztab(mztab_path, sample=line)

    truth = pd.read_csv(BENCH_DATA / line / "ground_truth.csv")
    by_basename: dict[str, list[PSMRecord]] = defaultdict(list)
    for r in records:
        by_basename[r.file_name].append(r)

    output_dir.mkdir(parents=True, exist_ok=True)
    for row in truth.itertuples(index=False):
        prop = float(row.nominal_proportion)
        sample = _proportion_to_sample(prop)
        mzml_basename = _mzml_basename(row.mzml_filename)
        mzml_src = DATA_ROOT / f"calibration_{line}" / "mzml" / row.mzml_filename
        if not mzml_src.exists():
            raise FileNotFoundError(f"mzML not found: {mzml_src}")

        psms_for_fraction = by_basename.get(mzml_basename, [])
        if not psms_for_fraction:
            print(f"[{line} {sample}] WARNING: no mzTab PSMs for {mzml_basename}; skipping")
            continue

        with tempfile.TemporaryDirectory(prefix=f"riana_mztab_{line}_{sample}_") as td:
            mzml_dir = Path(td) / "mzml"
            mzml_dir.mkdir()
            os.symlink(mzml_src, mzml_dir / row.mzml_filename)
            tsv_path = Path(td) / "percolator.target.psms.txt"
            write_mztab_as_percolator_tsv(psms_for_fraction, tsv_path)

            cmd = [
                sys.executable, "-m", "riana", "integrate",
                str(mzml_dir),
                str(tsv_path),
                "-s", sample,
                "-o", str(output_dir),
                *INTEGRATE_ARGS,
            ]
            print(f"[{line} {sample}] {' '.join(cmd)}", flush=True)
            if dry_run:
                continue
            t0 = time.time()
            subprocess.run(cmd, check=True)
            print(f"[{line} {sample}] done in {time.time() - t0:.1f}s", flush=True)

            expected = output_dir / f"{sample}_riana.txt"
            if not expected.exists():
                alt = output_dir / sample / f"{sample}_riana.txt"
                if alt.exists():
                    shutil.move(str(alt), str(expected))
                else:
                    raise FileNotFoundError(f"missing expected output: {expected}")
    return output_dir


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--line", required=True, choices=["ac16", "ipsc", "cm"])
    parser.add_argument(
        "--phase", choices=["psm", "integrate", "all"], default="psm",
        help="psm: PSM-level overlap only (default, fast). "
             "integrate: also run riana integrate on the mzTab-derived TSVs. "
             "all: run integrate, then m0/mA agreement (calls bench_m0_ma_recovery).",
    )
    parser.add_argument("--q-threshold", type=float, default=0.01)
    parser.add_argument(
        "--label", default="v0.9.0_mztab",
        help="Subdirectory name under benchmark_results/ and integrate_outputs/. "
             "Override to A/B two quantms searches without clobbering a "
             "committed baseline, e.g. --label v0.9.0_mztab_isoerr.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Default: tests/data/calibration_d2o_mixing/<line>/benchmark_results/<label>",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    line = args.line
    if args.output_dir is None:
        args.output_dir = BENCH_DATA / line / "benchmark_results" / args.label
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1 — always.
    print(f"\n=== phase 1: PSM overlap ({line}, q<{args.q_threshold}) ===")
    percolator = load_percolator_by_proportion(line, q_threshold=args.q_threshold)
    mztab, mztab_path = load_mztab_by_proportion(line, q_threshold=args.q_threshold)
    overlap = psm_overlap_table(percolator, mztab)
    print(overlap.to_string(index=False))

    prov = make_provenance(
        {"line": line, "q_threshold": args.q_threshold, "phase": args.phase},
        id_source=str(mztab_path),
    )
    write_tsv(
        args.output_dir / "psm_overlap.tsv",
        columns=list(overlap.columns),
        rows=overlap.to_dict(orient="records"),
        provenance=prov,
    )
    summary_path = args.output_dir / "psm_overlap_summary.json"
    summary = {
        "totals": {
            "percolator_psms": int(overlap["n_percolator"].sum()),
            "mztab_psms": int(overlap["n_mztab"].sum()),
            "intersection": int(overlap["n_intersection"].sum()),
            "jaccard_pooled": float(
                overlap["n_intersection"].sum() / overlap["n_union"].sum()
            ),
        },
        "per_proportion": overlap.to_dict(orient="records"),
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"[phase 1] wrote {args.output_dir}/psm_overlap.tsv + summary.json")

    if args.phase == "psm":
        return

    # Phase 2 — integrate against mzTab IDs.
    integrate_out = BENCH_DATA / line / "integrate_outputs" / args.label
    print(f"\n=== phase 2: integrate mzTab IDs -> {integrate_out} ===")
    run_integrate_from_mztab(line, output_dir=integrate_out, dry_run=args.dry_run)

    if args.phase == "integrate":
        print("[phase 2] integrate-only mode; m0/mA comparison skipped")
        return

    # Phase 3 — m0/mA recovery against both integrate outputs, same coefficients.
    print(f"\n=== phase 3: m0/mA recovery diff ===")
    coeff_csv = BENCH_DATA / line / f"d2o_aa_coefficients_{line}.csv"
    if not coeff_csv.exists():
        raise FileNotFoundError(
            f"frozen coefficients not found: {coeff_csv}. Re-run "
            "build_frozen_tables.py first."
        )
    for label, inputs in [
        ("percolator", BENCH_DATA / line / "integrate_outputs" / "v0.9.0"),
        ("mztab", integrate_out),
    ]:
        out_subdir = args.output_dir / f"m0_ma_recovery_{label}"
        out_subdir.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable, "-m", "tests.benchmark.bench_m0_ma_recovery",
            "--inputs", str(inputs),
            "--ground-truth", str(BENCH_DATA / line / "ground_truth.csv"),
            "--coefficients", str(coeff_csv),
            "--output-dir", str(out_subdir),
        ]
        print(f"[{label}] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"[{label}] wrote {out_subdir}")

    # Diff the two summary JSONs into a single comparison table.
    p_sum = json.loads(
        (args.output_dir / "m0_ma_recovery_percolator" / "m0_ma_recovery_summary.json").read_text()
    )
    m_sum = json.loads(
        (args.output_dir / "m0_ma_recovery_mztab" / "m0_ma_recovery_summary.json").read_text()
    )
    diff = {
        section: {
            metric: {
                "percolator": p_sum[section].get(metric),
                "mztab": m_sum[section].get(metric),
            }
            for metric in ("n", "n_peptides", "m0_rmse", "m0_mae", "m0_bias_median",
                           "env_rmse_median", "env_rmse_p95")
        }
        for section in ("all", "curated", "uncurated")
    }
    diff_path = args.output_dir / "m0_ma_diff.json"
    diff_path.write_text(json.dumps(diff, indent=2))
    print(f"[phase 3] wrote {diff_path}")
    for section in ("all", "curated", "uncurated"):
        p = p_sum[section]; m = m_sum[section]
        print(f"  {section:<10}  percolator: n={p['n']:5d} m0_rmse={p['m0_rmse']:.4f}   "
              f"mztab: n={m['n']:5d} m0_rmse={m['m0_rmse']:.4f}")


if __name__ == "__main__":
    main()
