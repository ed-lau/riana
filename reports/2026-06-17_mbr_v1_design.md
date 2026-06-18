# MBR v1 — design (mzTab/DDA path)

- **Date:** 2026-06-17
- **Status:** v1 integrate-core built; first real-data validation done (2026-06-18) → an MBR SNR/intensity floor is needed before fit (see [§ v1 validation](#v1-validation--first-real-data-run-2026-06-18))
- **Precursor:** [missingness + RT-alignment measurement](2026-06-17_mbr_dda_feasibility.md) (GO for DDA)
- **Roadmap:** `PROJECT_REVIEW.md` → Track A

Match-between-runs for the quantms **mzTab/DDA** path: transfer a confidently
identified precursor's identity + retention time into the runs of its turnover
curve that missed it, so curve points lost to stochastic MS2 sampling are
recovered. The [measurement](2026-06-17_mbr_dda_feasibility.md) established this is
worth doing (DDA curves 8–15% complete, ~40–45% recoverable gaps, 70% LVE
t0-anchor loss; DIA needs none) and feasible (quantms RT is aligned to ~4–9 s with
15–25 s run-specific residuals → a light per-run refinement suffices).

## Decisions

- **Ship v1 = pure RT-transfer.** Confident donors only; no mzTab changes.
- **Sub-threshold "rescue" tier is a hypothetical future improvement, not v1.**
  Promoting borderline acceptor IDs (q∈(0.01, 0.05]) would be the highest-confidence
  transfer (real MS2 at that RT), but the current mzTab is pre-filtered at ~1%
  PSM-FDR (247,723 PSMs ≤0.01; **179** in (0.01,0.02]; **0** above) — there is
  nothing to promote without a looser-FDR quantms re-export. Parked; revisit only if
  a looser export ever exists.
- **MBR FDR is its own future study/report,** not a v1 gate. v1 leans on the
  existing downstream R²>0.95 curation gate + an `evidence="mbr"` flag, and is
  validated empirically (below).

## Architecture

The mzTab is a whole-experiment file, so the cross-run donor assembly is free.
Hook in [`core/pipeline.plan_integration`](../riana/core/pipeline.py) after
`all_psms` is read and before the per-run split:

1. Group `all_psms` by `identity.group_key` = `(experiment, condition)` — one
   turnover curve. Foundation: a peptide elutes at ~the same RT across labeling
   times (D₂O changes the isotope envelope, not elution), confirmed by the
   alignment measurement.
2. `core/mbr.augment(all_psms, config)` appends synthetic MBR `PSMRecord`s.
3. The existing flow handles them: `resolve_rt_anchored_scans` already RT-anchors
   any `scan<0` row (acquisition-agnostic; selective — leaves real DDA scans
   untouched), the apex finder re-centers, extraction is unchanged.

**Integration-layer change (one, surgical):** the scan↔RT prefix-scramble guard
currently skips the *whole run* when any `scan<0` exists; once MBR rows are mixed
into a DDA run that would drop the guard for the real PSMs. Fix: guard only the
*originally* directly-scanned subset (captured before `resolve_rt_anchored_scans`),
so DIA stays a natural no-op and DDA+MBR keeps full protection. Behaviour is
identical until MBR rows exist, so it lands first as a safe prep refactor.

## Transfer algorithm (`core/mbr.py`)

Per curve group, per fraction (transfers never cross fractions):

1. **Observed sets** — per run, the precursor set (`concat = SEQUENCE_charge`).
2. **Donor gate** — a precursor is a donor if identified at **q≤`mbr_donor_q`** in
   **≥`mbr_min_donor_runs`** runs of the group (defaults 0.01 / 2; both tunable).
   Keep each donor's sequence/charge/`peptide_mass`/`mod_sites`/`protein_id`, its
   per-donor-run RT, and best q.
3. **Per-run RT alignment** — for each acceptor run, fit a **robust** RT map from
   the donor frame to that run from their co-identified precursors (median /
   Theil-Sen offset; the measurement showed near-constant offsets, so a global
   linear warp is unnecessary). Aligns the anchor to *this* run's frame.
4. **Emit** — for each (donor precursor, acceptor run that lacks it): a synthetic
   `PSMRecord(scan=−1, retention_time=aligned RT, evidence="mbr", q=0.0,
   …donor seq/charge/mass/mod_sites…, identity=acceptor identity)`.
5. Append to `all_psms`. The `evidence="mbr"` flag rides through to fit/rollup.

**Graceful failure (load-bearing, not an afterthought).** A transferred precursor
may have *no real peak* in the acceptor run (it was genuinely below detection, not
just unsequenced). The transfer must then **drop the row**, never integrate
baseline as signal. Mechanism: after the apex search for an MBR row, require a
detected apex within the search window; if none, discard the MBR row (counted +
logged, not written). This keeps MBR additive to yield without injecting noise.
(Real q-value PSMs are unaffected — they keep today's behaviour.)

*Open question — is the prominence gate strict enough? (spike, see below).* The
existing apex finder (`algorithms/peaks.detect_peak`) gates on
`prominence ≥ max(prominence_k · 1.4826·MAD, 1.0)` with `prominence_k=3.0` — a ~3σ
floor relative to the **trace's own** MAD noise, returning `None` if nothing clears
it. That floor is self-referential: on a near-pure-noise trace (precursor truly
absent) MAD is small, so a modest noise bump can still clear 3×MAD and produce a
*spurious* apex. So for MBR's absent-precursor case the prominence gate is likely a
necessary-but-insufficient first filter. v1 implements the drop on the existing
gate; **how often it false-fires must be measured** (spike), and an absolute
SNR/intensity floor for MBR rows may be needed. Backstop: a noise extraction won't
match the IsoSpec forward envelope, so the downstream R²>0.95 curation gate rejects
it at fit time regardless.

## Config surface

- **Integrate** (`IntegrationConfig` + `--mbr` …): `mbr: bool = False` (opt-in),
  `mbr_min_donor_runs: int = 2`, `mbr_donor_q: float = 1e-2`. (RT-alignment method
  is internal for v1: robust per-run offset.)
- **Fit / rollup opt-out** (`FitConfig` + `--exclude-mbr`): MBR points are *used by
  default* (the whole point), but `--exclude-mbr` filters `evidence=="mbr"` rows so
  the with/without A/B is a one-flag rerun and a cautious user can drop them.

## Output, marking, and data-point accounting

- **Export tables** already carry the `evidence` column (`q_value` | `mbr`); keep it
  prominent in `_riana.txt` and the fit/rollup outputs.
- **GUI** marks MBR points distinctly — a different point color in the Model curve
  and Protein curve views, and the evidence value shown on row select. The
  chromatogram view should indicate an MBR (RT-anchored) extraction.
- **Data-point breakdown columns (scope bundled per request).** Expand the fit/rollup
  "number of data points" reporting from a single `n` to a breakdown so the count of
  *clean* points is explicit:
  - `n_points` (total), `n_mbr` (MBR-transferred), `n_metox` (came from a Met-Ox
    peptidoform merged at fit per M7 `CHEMICAL_MODS={35}`), `n_both`, `n_clean`
    (neither). Applies per peptide (fit) and per protein (rollup).
  - This needs the fit aggregation to read `evidence` (MBR) and the merged
    peptidoform's mod state (Met-Ox) per data point — a fit-layer change downstream
    of the integrate core.

## Validation — preregistered readouts

Run each metric **with vs without** MBR (`--exclude-mbr` toggles it):

- **Primary accuracy — calibration ground truth.** The held-out approach is
  *intensity-biased* (real IDs were picked because they were abundant, so they
  transfer too easily to represent the genuinely-missing low-abundance peaks). Use
  the D₂O mixing series instead (`data/calibration_{ac16,cm,ipsc}`): it has a known
  θ at every mixing proportion *independent* of MS2 picking, so score MBR'd peaks
  against ground-truth θ / predicted m0 across the intensity range — including the
  hard low ones. Report bias + RMSE for MBR'd vs directly-IDed peaks.
- **Secondary — within-protein-θ variance** (Track D `bench_within_protein_theta`)
  on the in-vivo LVE/ATR sets: the guardrail. Should stay flat or *tighten*; if MBR
  injects noise it widens.
- **Yield / completeness:** # precursors past the R²>0.95 gate, # protein curves, and
  the `bench_missingness` curve-completeness + t0-anchor recovery (expect ↑,
  especially t0).
- **k_deg stability:** median/IQR overall, and specifically on curves that gained a
  t0 anchor — should sharpen, not shift systematically.
- **Held-out (demoted, caveated):** drop/restore real IDs as a *sanity* check only,
  acknowledging the intensity bias.
- **Breakdown** by donor count and post-alignment RT residual — to locate where
  transfer degrades.

## Open spikes (surfaced 2026-06-17; resolve within the relevant phase)

### `--depth` semantics with MBR + chemical mods
`depth` is the curve-qualification gate (minimum data before a peptide is fit), but
what it *counts* is already inconsistent and gets muddier with MBR + Met-Ox:

- **Today:** the manifest/SDRF path counts `len(group)` = **rows** (so biological
  replicates *and* Met-Ox-merged peptidoforms inflate it — `fitting.py:282`), while
  the legacy path counts `sample.nunique()` = **distinct samples** (`fitting.py:288`).
  Two different meanings.
- **With MBR:** MBR rows add to the count, and (in the 1-rep LVE case) add genuine new
  timepoints. Whether they *should* count toward `depth` is governed cleanly by
  `--exclude-mbr`: excluded → depth on clean rows; included → recovered timepoints
  (e.g. a restored t0) count, which is the intended benefit.
- **With Met-Ox:** the gate is on the *merged* fit_key group, so oxidized + unoxidized
  forms already pool into depth; the `n_metox`/`n_clean` breakdown exposes how much of
  a curve is mod- or MBR-derived.

**Recommendation:** harmonize `depth` to **distinct labeling timepoints**
(`labeling_time` nunique) on both paths — it is the kinetic-identifiability quantity
(a curve needs enough *distinct x* to estimate k), robust to replicate / peptidoform /
MBR multiplicity, and it fixes the rows-vs-nunique inconsistency. Keep raw point count
visible via `n_points`. A future **min-clean-depth** gate (≥N non-MBR / non-Met-Ox
timepoints) is the refinement if the breakdown shows curves over-reliant on transferred
or merged points.

**Why a spike, not a snap change:** flipping the manifest path from rows→timepoints is
a behaviour change (stricter for replicate-heavy data). Quantify first — how many
curves change qualification under rows vs distinct-timepoints, with/without MBR, on
LVE/ATR + the DIA (3 tp × 3 rep) set — then flip the default. (Per maintainer: this is
its own spike; Met-Ox-counting reactivity rides along.)

### Apex false-peak rate on absent precursors
See *graceful failure* above. Measure how often `detect_peak` returns a (spurious) apex
when a precursor is genuinely absent, to decide whether the 3×MAD prominence gate needs
an absolute SNR/intensity floor for MBR rows. Measurement: extract a set of precursors
at runs where they are confidently absent (present in ≪ donor count, far below donor
intensity), tabulate how many yield an apex and the prominence/intensity distribution
vs real peaks. Backstop already exists (downstream R²>0.95 envelope gate).

## Phasing

1. **Prep refactor** — scan↔RT guard on the directly-scanned subset (safe, behaviour-
   preserving). *(this turn)*
2. **Integrate core** — `IntegrationConfig` knobs + `core/mbr.py` (donor assembly,
   robust per-run RT offset, synthetic records) + `plan_integration` hook + graceful
   failure + `--mbr` CLI + tests. Deliverable: `integrate --mbr` emits flagged,
   RT-anchored MBR rows.
3. **Fit / rollup** — `--exclude-mbr` filter + the `n_points/n_mbr/n_metox/n_clean`
   breakdown columns.
4. **GUI** — MBR point coloring + evidence on row select.
5. **Validation** — calibration ground-truth bench (primary) + the Track D /
   completeness / k_deg readouts; A/B report.

## v1 validation — first real-data run (2026-06-18)

Auditable record of the first `integrate --mbr` on real data.

- **Run:** `riana integrate data/timeseries_lve_atr/mzml <mztab> --sdrf <…> --mbr -W 10 -o runs/lve_atr_mbr` — LVE+ATR, 24 runs, ±10 ppm (from SDRF), the apex defaults.
- **Analysis:** `python tests/benchmark/bench_mbr_quality.py --run runs/lve_atr_mbr`.
- **Metric note.** `m0 = iso0 / Σ(iso0..iso5)` is the **monoisotopic fraction**, the
  integrate-stage observable (θ is a fit-stage quantity). In D₂O labeling m0
  **declines** monotonically with labeling time (inverse of θ, which rises), so a
  *good* transferred point sits on its precursor's m0 decline. The "monotone
  corridor" flags points outside their `[next, prev]` real-neighbour bracket (±0.05);
  the "interp |Δm0|" is the residual from a linear interp of the bracketing real
  points. Both are **benchmarked against held-out real points** (drop a real point,
  predict from its real neighbours) so the bar is "as good as a real point," not
  "perfect" (real points themselves are ~94% in-corridor, |Δm0|~0.013, due to noise).

**Survival (graceful no-apex drop).** 73,337 / ~131,023 planned transfers survived
(**56%**); **44% dropped** for no detectable apex. (Drops = planned − surviving: the
per-run drop log is emitted in worker processes so it does not reach the main
logfile — to be surfaced via the result object.)

**MS1 signal (surviving MBR vs directly-identified rows).**

| evidence | n | iso0 p10/50/90 | frac iso0≤0 | m0 median |
|----------|---|----------------|------------|-----------|
| q_value | 247,723 | 87,245 / 710,014 / 12,639,672 | 0.4% | 0.317 |
| mbr | 73,337 | 1,297 / 42,722 / 420,978 | **0.0%** | 0.336 |

Survivors carry real intensity (≈17× lower median than real — expected, they were
missed for being scarce) and none are empty (0.0% iso0≤0) — the drop removed the
truly-absent traces.

**Trajectory sense (m0 decline vs held-out real).**

| chamber | interp \|Δm0\| MBR (med/p90) | interp \|Δm0\| real | corridor MBR | corridor real |
|---------|------------------------------|---------------------|--------------|---------------|
| LVE | 0.055 / 0.345 | 0.013 / 0.064 | 63.5% | 94.0% |
| ATR | 0.095 / 0.498 | 0.013 / 0.063 | 50.0% | 94.6% |

**Quality scales steeply with intensity** (MBR points with bracketing real
neighbours, n=28,681, by iso0 quintile):

| iso0 quintile | iso0 median | corridor | interp \|Δm0\| median |
|---------------|-------------|----------|----------------------|
| Q1 (low) | 2,041 | 22.8% | 0.228 |
| Q2 | 20,347 | 45.7% | 0.112 |
| Q3 | 62,581 | 60.7% | 0.067 |
| Q4 | 156,587 | 72.1% | 0.046 |
| Q5 (high) | 590,009 | 82.7% | 0.025 |
| *real (held-out)* | *710,014* | *~94%* | *0.013* |

**Verdict.** The no-apex drop is necessary but **insufficient**: ~40% of survivors
are off-trajectory (corridor 50–64% vs 94% real, |Δm0| ~7× worse), and quality
climbs monotonically with intensity (Q1 23% → Q5 83% ≈ real). The relative 3×MAD
prominence gate passes too many wrong-peak (co-eluting / noise) picks. **v1 needs an
MBR SNR/intensity floor** (the apex-spike's answer — yes), a tunable yield-vs-quality
dial, with the downstream envelope-fit R²>0.95 gate as the backstop. **Do not fit
floor-less MBR.** Next: add the floor, re-run, and re-measure corridor% toward real.
