"""Track D: within-protein fractional-labeling spread on the animal D₂O series.

The mixing-series benches score observed-vs-predicted envelopes against a *known*
heavy proportion (bench_m0_ma_recovery). The animal in-vivo series has **no
fractional-pool ground truth** — labeling is a turnover process, not a bench
mixture — so it is scored by an *internal-consistency* metric (Hammond et al.
2022, [[m3-animal-within-protein-metric]]):

    A good integrator MINIMIZES the spread of fractional labeling (θ) among the
    peptides of one protein, within a single time point — because every peptide
    of a protein turns over together, so any within-protein θ spread is noise the
    integration injected. Tighter spread ⇒ better integration.

**Terminology:** a **"cell"** in this bench = a *(protein, timepoint)*
combination — one protein's peptides at one labeling timepoint — **not** a
biological cell or cell type. ``n_cells`` in the outputs is the number of scored
(protein, timepoint) combinations (≥ ``min_pep`` peptides each), i.e. the sample
size behind each within-protein-θ median.

θ per (peptide, timepoint) is solved with the **production** IsoSpec forward
model — :func:`solve_fs_d2o` against an Spep from the literature mammalian
coefficient table (``commerford``; Spep is trusted input here, not fit) at the
sample's precursor enrichment. This is exactly what ``riana fit`` computes per
timepoint (``FitResult.fs``), so the long-format θ this bench materializes is
also a preview of the M5 per-timepoint-FS substrate.

**Two populations, like bench_m0_ma_recovery.** Curation (here, kinetic-fit
R²(θ~t) > 0.95 under the simple model) discards exactly the low-SNR peptides
where integration improvements show up ([[m3-curation-gate-revisit]]), so the
spread is reported for ``all`` / ``curated`` / ``uncurated`` separately. It is
also reported **per abundance stratum** (hi/mid/lo from the frozen set) — we
stratify rather than exclude the low-abundance tail.

**Bias control.** Membership comes only from the frozen, integrator-independent
set built by ``build_lve_bench_set.py`` (proteotypic peptides of ≥3-peptide
proteins at q<0.01). Every method is scored on that byte-identical set, so the
``compare_methods`` table reflects the integrator, not the panel.

Usage (one or more methods; each is an integrate output dir's manifest):

    python bench_within_protein_theta.py \\
        --method baseline=runs/lve_fixed/riana_manifest.tsv \\
        --bench-set tests/benchmark/bench_sets/bench_lve_peptides.tsv \\
        --output-dir runs/lve_fixed/bench_theta

Outputs (under --output-dir):
  - within_protein_theta_cells.csv     (per method × (condition, protein, t) cell)
  - within_protein_theta_summary.json  (method × {all,curated,uncurated} × strata)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from riana.algorithms.isotope_dist import (  # noqa: E402
    clear_envelope_cache,
    solve_fs_d2o,
    spep_from_coefficients,
)
from riana.config import FitConfig  # noqa: E402
from riana.core.fitting import fit_run, load_aa_coefficients  # noqa: E402
from riana.core.pipeline import recombine_for_fit  # noqa: E402
from riana.io.manifest import read_manifest  # noqa: E402

STRATA = ("lo", "mid", "hi")
CURATION_R2 = 0.95

#: A **"cell"** here is a *(protein, timepoint)* combination — one protein's set
#: of peptides at one labeling timepoint — NOT a biological cell / cell type.
#: (The full key also carries condition + abundance stratum, which for the LVE
#: single-condition series are determined by the protein.) ``n_cells`` in the
#: outputs is therefore the count of scored (protein, timepoint) combinations,
#: i.e. the sample size behind each within-protein-θ median.
_CELL_KEYS = ["condition", "protein_id", "stratum", "labeling_time"]


# --------------------------------------------------------------------------- #
# per-timepoint θ (production forward model) + kinetic-fit curation
# --------------------------------------------------------------------------- #
def _iso_cols(frame: pd.DataFrame) -> list[str]:
    return sorted(
        (c for c in frame.columns if re.fullmatch(r"iso\d+", c)),
        key=lambda c: int(c[3:]),
    )


def compute_theta_long(
    frame: pd.DataFrame,
    coeffs: dict[str, float],
    ria_max: float,
    q_value: float,
) -> pd.DataFrame:
    """θ per (peptide, timepoint) via the production ``solve_fs_d2o``.

    One row per surviving ``(concat, biological_replicate, labeling_time)`` — the
    *full* q-passing population (no kinetic-fit / depth filter), so the uncurated
    tail is preserved for honest two-population reporting.
    """
    iso_cols = _iso_cols(frame)
    n_iso = len(iso_cols)
    df = frame[frame["percolator q-value"] < q_value].copy()
    if df.empty:
        return pd.DataFrame(
            columns=["concat", "sequence", "biological_replicate",
                     "labeling_time", "theta"]
        )

    seqs = df["sequence"].astype(str).to_numpy()
    masses = df["peptide mass"].to_numpy(dtype=float)
    iso = df[iso_cols].to_numpy(dtype=float)
    times = df["labeling_time"].to_numpy(dtype=float)
    concats = df["concat"].astype(str).to_numpy()
    bioreps = (
        df["biological_replicate"].to_numpy()
        if "biological_replicate" in df.columns
        else np.ones(len(df), dtype=int)
    )

    clear_envelope_cache()
    spep_cache: dict[str, int] = {}
    thetas = np.full(len(df), np.nan)
    for i in range(len(df)):
        seq = seqs[i]
        spep = spep_cache.get(seq)
        if spep is None:
            spep = max(1, int(round(spep_from_coefficients(seq, coeffs))))
            spep_cache[seq] = spep
        try:
            thetas[i] = solve_fs_d2o(
                seq, masses[i], iso[i], spep, ria_max=ria_max, n_iso=n_iso
            )
        except (KeyError, ValueError):
            thetas[i] = np.nan  # non-canonical residue etc.; drop cleanly

    return pd.DataFrame({
        "concat": concats,
        "sequence": seqs,
        "biological_replicate": bioreps,
        "labeling_time": times,
        "theta": thetas,
    })


def kinetic_r2(
    frame: pd.DataFrame, config: FitConfig, coeffs: dict[str, float]
) -> dict[str, float]:
    """``{concat: R²(θ~t)}`` from the production kinetic fit (curation signal).

    ``n_boot=0`` skips the (expensive, here-unused) k_deg bootstrap. Peptides
    dropped by the fit's depth filter are simply absent ⇒ treated as uncurated.
    """
    fitted = fit_run(
        config, [frame], coeffs, time_column="labeling_time", n_boot=0
    )
    return {str(c): float(r) for c, r in zip(fitted.index, fitted["R_squared"])}


def score_method(
    manifest_path: Path,
    frozen_pep: pd.DataFrame,
    config: FitConfig,
    coeffs: dict[str, float],
    ria_default: float,
) -> pd.DataFrame:
    """Manifest → per-(condition, protein, stratum, sequence, timepoint) θ table.

    Charges and biological replicates of the same peptide are collapsed to one θ
    per (sequence, timepoint) by median; ``is_curated`` is OR-ed across them.
    """
    rows = read_manifest(manifest_path, stage="integrate")
    if not rows:
        raise SystemExit(f"no integrate rows in {manifest_path}")
    # Per-sample RIA: SDRF identity → global default (locked decision #5).
    ria_max = next(
        (r.identity.precursor_enrichment for r in rows
         if r.identity.precursor_enrichment is not None),
        ria_default,
    )
    curves = recombine_for_fit(rows)

    parts = []
    for (_experiment, condition), frame in curves.items():
        theta_long = compute_theta_long(frame, coeffs, ria_max, config.q_value)
        if theta_long.empty:
            continue
        r2 = kinetic_r2(frame, config, coeffs)
        theta_long["condition"] = condition or ""
        theta_long["is_curated"] = (
            theta_long["concat"].map(r2).fillna(-1.0) > CURATION_R2
        )
        parts.append(theta_long)
    if not parts:
        raise SystemExit(f"no θ computed for {manifest_path}")
    theta_long = pd.concat(parts, ignore_index=True)

    # Join the frozen, integrator-independent membership (proteotypic peptides of
    # ≥3-peptide proteins). This both restricts the population and assigns the
    # canonical single protein_id + abundance stratum.
    merged = theta_long.merge(
        frozen_pep[["sequence", "protein_id", "stratum"]],
        on="sequence", how="inner",
    )
    if merged.empty:
        raise SystemExit(
            f"no frozen-set peptides found in {manifest_path} outputs "
            "(sequence join empty — check the bench set matches this mzTab)."
        )

    seq_cell = (
        merged.dropna(subset=["theta"])
        .groupby(
            ["condition", "protein_id", "stratum", "sequence", "labeling_time"],
            as_index=False, observed=True,
        )
        .agg(theta=("theta", "median"), is_curated=("is_curated", "max"))
    )
    return seq_cell


# --------------------------------------------------------------------------- #
# spread aggregation
# --------------------------------------------------------------------------- #
def cell_spreads(seq_cell: pd.DataFrame, min_pep: int) -> pd.DataFrame:
    """Per (condition, protein, stratum, timepoint) within-protein θ spread.

    A cell qualifies only with ≥ *min_pep* distinct peptides bearing a θ at that
    timepoint. ``theta_robust_sd`` = 1.4826·MAD (the headline; outlier-robust);
    ``theta_iqr`` and plain SD are reported alongside.
    """
    empty = pd.DataFrame(columns=_CELL_KEYS + [
        "n_pep", "theta_median", "theta_iqr", "theta_robust_sd", "theta_sd",
        "frac_curated",
    ])
    if seq_cell.empty:
        return empty
    grp = seq_cell.groupby(_CELL_KEYS, observed=True)
    cell = grp["theta"].agg(
        n_pep="size",
        theta_median="median",
        q1=lambda v: float(np.percentile(v, 25)),
        q3=lambda v: float(np.percentile(v, 75)),
        theta_sd=lambda v: float(np.std(v, ddof=1)) if len(v) > 1 else np.nan,
        mad=lambda v: float(np.median(np.abs(v - np.median(v)))),
    )
    cell["frac_curated"] = grp["is_curated"].mean()
    cell = cell.reset_index()
    cell["theta_iqr"] = cell["q3"] - cell["q1"]
    cell["theta_robust_sd"] = 1.4826 * cell["mad"]
    cell = cell[cell["n_pep"] >= min_pep].drop(columns=["q1", "q3", "mad"])
    return cell.reset_index(drop=True)


def _agg_cells(cell: pd.DataFrame) -> dict:
    if cell.empty:
        return {"n_cells": 0, "n_proteins": 0, "n_timepoints": 0,
                "theta_robust_sd_median": float("nan"),
                "theta_iqr_median": float("nan"),
                "theta_robust_sd_p90": float("nan"), "per_stratum": []}
    per_stratum = []
    for s in STRATA:
        sub = cell[cell["stratum"] == s]
        per_stratum.append({
            "stratum": s,
            "n_cells": int(len(sub)),
            "n_proteins": int(sub["protein_id"].nunique()),
            "theta_robust_sd_median": (
                float(np.median(sub["theta_robust_sd"])) if len(sub)
                else float("nan")
            ),
            "theta_iqr_median": (
                float(np.median(sub["theta_iqr"])) if len(sub) else float("nan")
            ),
        })
    return {
        "n_cells": int(len(cell)),
        "n_proteins": int(cell["protein_id"].nunique()),
        "n_timepoints": int(cell["labeling_time"].nunique()),
        "theta_robust_sd_median": float(np.median(cell["theta_robust_sd"])),
        "theta_iqr_median": float(np.median(cell["theta_iqr"])),
        "theta_robust_sd_p90": float(np.percentile(cell["theta_robust_sd"], 90)),
        "per_stratum": per_stratum,
    }


def summarize(seq_cell: pd.DataFrame, min_pep: int) -> dict:
    """``{all,curated,uncurated}`` × overall + per-stratum spread aggregates."""
    return {
        "all": _agg_cells(cell_spreads(seq_cell, min_pep)),
        "curated": _agg_cells(
            cell_spreads(seq_cell[seq_cell["is_curated"]], min_pep)
        ),
        "uncurated": _agg_cells(
            cell_spreads(seq_cell[~seq_cell["is_curated"]], min_pep)
        ),
    }


# --------------------------------------------------------------------------- #
# multi-method comparison (mirrors bench_m0_ma_recovery.compare_integrations)
# --------------------------------------------------------------------------- #
def compare_methods(
    methods: dict[str, Path],
    frozen_pep: pd.DataFrame,
    config: FitConfig,
    coeffs: dict[str, float],
    ria_default: float,
    min_pep: int,
) -> tuple[pd.DataFrame, dict]:
    """Score each named method's manifest on the frozen set; tabulate by method.

    Returns ``(all-population cell table with a ``method`` column, per-method
    summaries)``. In Track D only the baseline exists; a Track B integrator
    variant registers as another method and the same table shows whether the
    within-protein spread tightened.
    """
    cells, summaries = [], {}
    for name, manifest in methods.items():
        t0 = time.time()
        print(f"[score] method {name!r} <- {manifest}", flush=True)
        seq_cell = score_method(manifest, frozen_pep, config, coeffs, ria_default)
        summaries[name] = summarize(seq_cell, min_pep)
        all_cell = cell_spreads(seq_cell, min_pep)
        all_cell.insert(0, "method", name)
        cells.append(all_cell)
        a = summaries[name]["all"]
        print(
            f"       {a['n_cells']} cells / {a['n_proteins']} proteins / "
            f"{a['n_timepoints']} timepoints  "
            f"θ_robust_sd(med)={a['theta_robust_sd_median']:.4f}  "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )
    return pd.concat(cells, ignore_index=True), summaries


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--method", action="append", required=True, metavar="LABEL=MANIFEST",
        help="A method to score, e.g. baseline=runs/lve/riana_manifest.tsv. "
        "Repeat to compare integrators on the same frozen set.",
    )
    parser.add_argument(
        "--bench-set", type=Path, required=True,
        help="bench_lve_peptides.tsv from build_lve_bench_set.py.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--coefficients", default="commerford",
        help="Per-AA Spep coefficient table — preset name or CSV path "
        "[default: commerford (literature mammalian)].",
    )
    parser.add_argument(
        "--ria", type=float, default=0.06,
        help="Fallback precursor enrichment (RIA max) if absent from the "
        "manifest identity [default: 0.06].",
    )
    parser.add_argument("--q-value", type=float, default=0.01)
    parser.add_argument(
        "--depth", type=int, default=3,
        help="Min timepoints for a peptide's kinetic-fit curation R² "
        "[default: 3].",
    )
    parser.add_argument(
        "--min-peptides", type=int, default=3,
        help="Min distinct peptides for a (protein, timepoint) cell to score "
        "[default: 3].",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    methods: dict[str, Path] = {}
    for spec in args.method:
        if "=" not in spec:
            raise SystemExit(f"--method must be LABEL=PATH, got {spec!r}")
        label, path = spec.split("=", 1)
        methods[label] = Path(path)

    frozen_pep = pd.read_csv(args.bench_set, sep="\t")
    frozen_pep["stratum"] = pd.Categorical(
        frozen_pep["stratum"], categories=list(STRATA), ordered=True
    )
    coeffs = load_aa_coefficients(args.coefficients)
    config = FitConfig(q_value=args.q_value, depth=args.depth)

    print(
        f"[load] bench set: {frozen_pep['protein_id'].nunique()} proteins, "
        f"{len(frozen_pep)} peptides; coefficients={args.coefficients}"
    )
    cells, summaries = compare_methods(
        methods, frozen_pep, config, coeffs, args.ria, args.min_peptides
    )

    cells_path = args.output_dir / "within_protein_theta_cells.csv"
    cells.to_csv(cells_path, index=False)
    summary_path = args.output_dir / "within_protein_theta_summary.json"
    with summary_path.open("w") as f:
        json.dump(
            {"methods": summaries,
             "args": {"bench_set": str(args.bench_set),
                      "coefficients": args.coefficients, "ria": args.ria,
                      "q_value": args.q_value, "depth": args.depth,
                      "min_peptides": args.min_peptides}},
            f, indent=2,
        )

    print("\n[done] within-protein θ spread (lower = tighter = better):")
    for name, summ in summaries.items():
        for pop in ("all", "curated", "uncurated"):
            a = summ[pop]
            strata = "  ".join(
                f"{s['stratum']}={s['theta_robust_sd_median']:.4f}"
                for s in a["per_stratum"]
            )
            print(
                f"  {name:>10} {pop:>9}  n={a['n_cells']:5d}  "
                f"θ_robust_sd(med)={a['theta_robust_sd_median']:.4f}   [{strata}]"
            )
    print(f"       cells   -> {cells_path}")
    print(f"       summary -> {summary_path}")


if __name__ == "__main__":
    main()
