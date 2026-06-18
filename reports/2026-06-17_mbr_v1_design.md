# MBR v1 — design (mzTab/DDA path)

- **Date:** 2026-06-17
- **Status:** design agreed, implementing
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
detected apex passing the normal prominence/SNR test within the search window; if
none, discard the MBR row (counted + logged, not written). This keeps MBR additive
to yield without injecting noise. (Real q-value PSMs are unaffected — they keep
today's behaviour.)

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
