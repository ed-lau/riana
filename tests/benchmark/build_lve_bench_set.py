"""Freeze the canonical within-protein-θ bench set for the animal D₂O series.

Track D needs a protein/peptide subset to score the within-protein fractional-
labeling spread metric ([[m3-animal-within-protein-metric]]). The cardinal rule
(PROJECT_REVIEW.md §3, Track D discussion 2026-06-07) is that **membership must
be integrator-independent** — computable from the search/ID output (the mzTab)
*before* any Riana integration — so every candidate integrator is later scored on
the byte-identical set and relative comparisons stay apples-to-apples. Anything
the integrator influences (θ, FS, R², integrated intensity, SNR) must NOT define
membership, or each integrator would redraw its own favourable panel.

This builder reads the quantms mzTab and emits that frozen set. Selection uses
only ID-side features:

  1. **q < --q-value** — search-engine ID confidence (default 0.01).
  2. **Proteotypic only** — a PSM whose ``accession`` names a single protein
     (no comma-joined group). A peptide shared across proteins cannot be
     attributed to one protein's spread, so razor/shared peptides are dropped.
  3. **≥ --min-peptides distinct sequences per protein** (default 3) — you need
     ≥3 peptides to form a non-degenerate within-protein spread.
  4. **Abundance stratum** — proteins are split into hi/mid/lo tertiles by an
     integrator-independent abundance proxy (spectral count = number of q-passing
     PSMs per protein). We *stratify and report per tertile*, we do NOT
     subsample to balance — low-abundance peptides are precisely the population
     the Track B integration work targets ([[m3-curation-gate-revisit]]), so the
     bench must keep them and surface their behaviour, not design it away.

Because Riana never changes which peptides are *identified* (quantms owns the
search), freezing on this set cannot disadvantage a future Riana integrator — it
only changes how well a candidate integrates a fixed input. The one axis that
*would* change the ID set is a different quantms search config, pinned separately
([[m3-quantms-search-config]]).

Outputs (written next to --out, default ``tests/benchmark/bench_sets/``):
  - ``bench_lve_peptides.tsv``  (protein_id, sequence, n_psms, stratum) — the
    frozen membership the bench joins on ``sequence``.
  - ``bench_lve_proteins.tsv``  (protein_id, n_peptides, n_psms, stratum) — the
    protein-level census used for stratified reporting.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Reuse the production mzTab parser so q-value / decoy / accession handling is
# exactly what `riana integrate` sees (no second, drifting reader).
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from riana.io.mztab import read_mztab  # noqa: E402
from riana.utils import strip_concat  # noqa: E402

STRATA = ("lo", "mid", "hi")


def build_bench_set(
    mztab_path: Path,
    q_value: float,
    min_peptides: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(peptides_df, proteins_df)`` — the frozen, ID-side bench set.

    Pure function over the mzTab (no Riana integration), so the result depends
    only on the search output and the two thresholds.
    """
    # Legacy single-label call: we only need accession / sequence / q-value, not
    # per-run identity, so a bare `sample=` (no SDRF) is the right surface here.
    records, _ = read_mztab(mztab_path, sample="lve")
    if not records:
        raise SystemExit(f"no PSMs parsed from {mztab_path}")

    rows = []
    for r in records:
        if r.percolator_q_value >= q_value:
            continue
        acc = (r.protein_id or "").strip()
        # Proteotypic only: a single accession, no comma-joined protein group.
        if not acc or "," in acc:
            continue
        rows.append({"protein_id": acc, "sequence": strip_concat(r.sequence)})
    if not rows:
        raise SystemExit(
            f"no proteotypic PSMs at q < {q_value} in {mztab_path}"
        )
    psm_df = pd.DataFrame(rows)

    # Peptide = distinct (protein, stripped sequence). n_psms across all charges
    # / scans is the spectral-count abundance signal at the peptide level.
    pep = (
        psm_df.groupby(["protein_id", "sequence"], as_index=False)
        .size()
        .rename(columns={"size": "n_psms"})
    )

    # Belt-and-braces proteotypicity: a sequence reported under a single
    # accession in some PSMs but a *different* single accession in others is not
    # truly protein-unique — drop it so no peptide contributes to two proteins'
    # spreads.
    seq_protein_count = pep.groupby("sequence")["protein_id"].transform("nunique")
    pep = pep[seq_protein_count == 1].copy()

    # Protein census: distinct peptides + total spectral count.
    prot = pep.groupby("protein_id", as_index=False).agg(
        n_peptides=("sequence", "nunique"),
        n_psms=("n_psms", "sum"),
    )
    prot = prot[prot["n_peptides"] >= min_peptides].copy()
    if prot.empty:
        raise SystemExit(
            f"no proteins with ≥ {min_peptides} proteotypic peptides at "
            f"q < {q_value}"
        )

    # Abundance tertiles by spectral count. `rank` breaks the heavy ties spectral
    # counts produce so qcut sees a (near-)continuous distribution; duplicates=
    # 'drop' degrades gracefully if a tertile edge still collapses.
    ranked = prot["n_psms"].rank(method="first")
    try:
        prot["stratum"] = pd.qcut(ranked, 3, labels=list(STRATA))
    except ValueError:
        prot["stratum"] = pd.cut(
            ranked, bins=3, labels=list(STRATA), include_lowest=True
        )
    prot = prot.sort_values(
        ["stratum", "n_psms"], ascending=[True, False]
    ).reset_index(drop=True)

    # Push the protein stratum down to its peptides for the join key the bench
    # uses.
    pep = pep.merge(
        prot[["protein_id", "stratum"]], on="protein_id", how="inner"
    )
    pep = pep.sort_values(["protein_id", "sequence"]).reset_index(drop=True)
    return pep, prot


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--mztab", type=Path, required=True,
        help="quantms/OpenMS *.mzTab for the animal time series.",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).resolve().parent / "bench_sets",
        help="Output directory for the frozen TSVs "
        "[default: tests/benchmark/bench_sets/].",
    )
    parser.add_argument("--q-value", type=float, default=0.01)
    parser.add_argument("--min-peptides", type=int, default=3)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"[load] mztab={args.mztab}")
    pep, prot = build_bench_set(args.mztab, args.q_value, args.min_peptides)

    pep_path = args.out / "bench_lve_peptides.tsv"
    prot_path = args.out / "bench_lve_proteins.tsv"
    pep.to_csv(pep_path, sep="\t", index=False)
    prot.to_csv(prot_path, sep="\t", index=False)

    counts = prot["stratum"].value_counts().reindex(STRATA, fill_value=0)
    print(
        f"[done] {len(prot)} proteins (≥{args.min_peptides} proteotypic peptides "
        f"at q<{args.q_value}), {len(pep)} peptides"
    )
    print(
        f"       strata  lo={counts['lo']}  mid={counts['mid']}  hi={counts['hi']}"
        f"   (spectral-count tertiles)"
    )
    print(f"       -> {pep_path}")
    print(f"       -> {prot_path}")


if __name__ == "__main__":
    main()
