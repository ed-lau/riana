# Dimethyl 0/8 duplex — light→heavy spillover go/no-go

**Date:** 2026-07-06 · **Branch:** `1.2.0` · **Bench:**
`tests/benchmark/bench_dimethyl_spillover.py` · **Data:** `data/singlepoint_dimethyld2o`

## Question

The next 1.2.0 item is **true (sample-axis) multiplexing** for the dimethyl-D₂O
duplex: fit each channel separately, combine at rollup as experiments. Before
building the channel→sample intake we need to know whether the two channels are
cleanly separable at MS1, or whether the near-term build must also ship a
light→heavy **demux** (spillover subtraction).

The duplex tags **control liver = light dimethyl** (UNIMOD:36, +28.0313/site) and
**rapamycin liver = +8 heavy dimethyl** (UNIMOD:330, +36.0757/site). Heavy is
**+8.0444 Da per dimethyl site** above light, and a peptide carries
**S = 1 (N-term) + #K** sites:

| S | dimethyl sites | heavy cluster offset | where heavy iso0 lands |
|---|----------------|----------------------|------------------------|
| **1** | N-term only (C-term R, **no internal K**) | **+8.04 Da** | **light iso8** ← closest |
| 2 | N-term + one K | +16.09 Da | light iso16 |
| 3 | N-term + two K | +24.13 Da | light iso24 |

So the **only** population at risk is **S=1** (R-terminated, no K). Its heavy cluster
iso0 sits on light iso8, and the concern is that the light cluster — broadened by
D₂O labelling — spills into the heavy FS-scoring window. Every K-bearing peptide
(S≥2, **72.5 %** of the tryptic population) has its heavy cluster ≥16 Da away and is
**completely non-overlapping** (p99 spill into the heavy window ≈ 0.00 %).

## Method

RIANA's **production** forward model (`isotope_dist.get_peptide_distribution` /
`get_envelope`) over an in-silico tryptic digest of the reviewed mouse proteome
(0–2 missed cleavages, len 7–35). No search results or fit output exist for this
dataset yet — and the intake needed to produce them is the very thing this test
informs — so we characterise across the realistic peptide space and sweep the
fraction-new (FS) instead of reading one dataset's measured FS. The observed light
envelope at fraction FS **is** the production mixture `(1-FS)·init + FS·final`, built
at the dataset's actual **light** enrichment (**4.614 % D₂O**, from the SDRF; note
this is *below* the 6 % cell-culture default → less broadening), with per-peptide
labile-site counts from the mouse in-vivo `deberneh_2025_rss` table.

**Decision metric** — for an S=1 peptide at balanced abundance (L=H), the fraction of
the **heavy iso0..iso5 window** (light indices 8..13, what the heavy FS solve scores)
that is actually **light** contamination:
`light_in_window / (light_in_window + heavy_in_window)`. It scales **linearly** with
the true L/H ratio (a heavy peptide more abundant than its light partner is
proportionally *less* corrupted, and vice-versa).

FS=1.0 is the conservative worst case (maximum broadening); day-8 mouse liver — fast
turnover — is a **high-FS** regime, so the high end of the sweep is the operative one.

## Result — the spill is entirely a function of peptide length

S=1 heavy-window contamination (% of the heavy iso0..5 signal that is light, L=H):

| length | FS 0.5 med / p90 | FS 0.9 med / p90 | FS 1.0 med / p90 | FS 1.0: >5% / >10% |
|--------|-----------------|-----------------|-----------------|--------------------|
| **7–12**  | 0.00 / 0.04 | 0.00 / 0.08 | 0.00 / 0.09 | 0 % / 0 % |
| **13–18** | 0.19 / 0.58 | 0.34 / 1.08 | 0.38 / 1.22 | 0 % / 0 % |
| **19–25** | 1.37 / 3.35 | 2.61 / 6.65 | 2.94 / 7.61 | 25 % / 4 % |
| **26–35** | 6.41 / 12.4 | 13.1 / 26.1 | 15.2 / 30.7 | **96 % / 75 %** |

The driver is natural-abundance envelope width: a long peptide's ¹³C envelope already
reaches iso6–8 *before* any D₂O broadening, so at high FS its light tail lands
squarely on the heavy iso0. Short peptides have no tail to spill.

Length-weighted over a realistic DDA-detectability profile, the S=1 population as a
whole: at FS=1.0, **~11 %** of heavy windows are >5 % light-contaminated, **~6 %**
>10 %, **~2 %** >20 %. S=1 is 27.5 % of the tryptic population, so the badly-corrupted
set (>10 %) is on the order of **~1–2 % of all quantifiable peptides**, concentrated
in the long-and-R-terminated corner that DDA quantifies least reliably anyway.

## Decision — **GO** for the near-term duplex build, with a cheap S=1 guard

The 0/8 duplex is genuinely separable at MS1. **Full demux is NOT a near-term
blocker.** Concretely, no subtraction is needed for:

- **all S≥2 peptides** (K-bearing, 72.5 % of tryptic) — heavy at ≥+16 Da, zero overlap;
- **all S=1 peptides ≤ ~18 residues** — spill is negligible at any FS (p99 < 2.5 %);
- **most S=1 peptides 19–25 residues** — median ~3 %, with a modest 5–8 % tail.

The **only** compromised population is **long S=1 peptides (≥ ~26 residues, R-terminal,
no K)** at high FS — badly corrupted (median 15 %, tail >30 %). Rather than block the
build on demux, ship the channel→sample duplex now and add a **deterministic S=1
spillover guard**: for a peptide with S=1, the forward-model light-tail into iso8+ is
computable at fit time from its sequence + spep + enrichment alone, so the at-risk
peptides can be gated *exactly* (not by a blunt length cutoff) and their heavy-channel
FS dropped or flagged. (A plain length gate at ~20–22 residues for S=1 is the trivial
first approximation.)

**Full light→heavy subtraction demux** (Track C "proper demultiplexing") only earns
its keep if we later want to *rescue* those long S=1 peptides rather than curate them
out — a clean future item, not a prerequisite. It becomes mandatory for the tighter
spacings SILAC/dimethyl will eventually bring (DIMETHYL2/4/6, or SILAC medium
channels), where even short peptides overlap.

## Caveats

- **L=H assumption.** Contamination scales linearly with the per-peptide light:heavy
  total-abundance ratio. A control-vs-rapamycin liver duplex is ~balanced by design,
  but individual proteins vary; the guard should key on the forward-model spill, which
  is abundance-independent, and treat L=H as the nominal severity.
- **Forward-model, not measured.** Once the channel→sample intake exists and the data
  is fit, re-run against the **actual** measured FS and identified S=1 population — but
  the physics bound here is decisive for the separable majority regardless.
- **4.614 % enrichment is on the low side.** A higher-enrichment (e.g. 6 %) duplex would
  broaden more and push the length threshold down; the guard should scale with the
  SDRF `precursor enrichment`, not a hard-coded length.
