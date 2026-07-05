# Riana — Roadmap & TODO

> **Status (2026-07-02):** version `1.1.0` **fully released** — PyPI (`riana==1.1.0`),
> GitHub release + Zenodo DOI, and `master` merged back. The experimental-science + QoL
> line is shipped; active development continues on the `1.1.0` branch, and the next line
> is **1.2.0** (TMT / multiplexing, deamidation, the adaptive-N_ISO + robust-envelope
> rework). The shipped record lives in `CHANGELOG.md`; **this document is the
> forward-looking roadmap** — open work is §3 (the Tracks) and Known Limitations (§6).
> §2 and M1–M4 are kept as collapsed pointers because code docstrings reference their
> anchors (§2b/c/d, §3, Tracks A–E).
>
> Maintainer: Edward Lau.

## 1. Project status

Riana is a single-author scientific Python tool for extracting and modelling
isotopomer time-series from MS1 data, used in protein-turnover research. The
scientific core (accurate-mass, kinetic models, fractional-synthesis math, IsoSpec)
is lifted-unchanged and proven; the 1.0.0 rewrite replaced the surrounding software
(typed `core/` + `algorithms/` + `io/`, streaming mzML, SDRF/manifest intake, peak
detection, a PySide6 GUI). The package layout is documented in the **README**; the
as-built feature record is **`CHANGELOG.md`**; the open work is §3 + §6.

## 2. Critical findings (post-evaluation)

### 2b. Scientific defects deferred to 1.0.0

All **fixed in 1.0.0** (see `CHANGELOG.md [1.0.0]`): the unreachable amino-acid
`a_0`/`label` dispatch (`core/fsynthesis`), the `iso0/colsums` FS-denominator
drift, the heuristic kinetic-CI (now a residual bootstrap), and the fixed
site-count FS bias — FS now comes from the IsoSpec per-peptide-Spep forward/solve
model, closing the ≈ −0.5 pseudo-time `k_deg` bias.

### 2c. Algorithmic feature gaps — addressed in 1.0.0 (peak detection, baseline, mass-domain)

The 1.0.0 integration rewrite addressed these (see `CHANGELOG.md` + the
`m3-peak-detection-spike` memory): chromatographic apex/consensus peak detection
(the apex-narrow window is now the default), in-window baseline options, and
per-isotopomer observed-mass / drift tracking. SG smoothing is opt-in (it distorts
areas → off by default). **Still open / not fully resolved:** robust in-window
**baseline subtraction** (the `linear` baseline was discarded; `none` is the
default) and a **cross-proportion-stable peak picker** — both live
integration-fidelity items under Track B / Known Limitations.

### 2d. Architecture findings (addressed in 1.0.0 rewrite)

All addressed by the 1.0.0 rewrite (see `CHANGELOG.md`): CLI/GUI no longer
duplicate validation (one frozen `IntegrationConfig`/`FitConfig`), the monolithic
`integrate_all` is gone (typed `core/` modules + async GUI off the main thread),
the broken Tkinter GUI is replaced by PySide6, and the global-state modules that
blocked isolated test runs are restructured.

## 3. Roadmap

Sequencing: M1 (0.9.0), M2 (calibration dataset + benchmarks), M3 (restructure
+ peak-detection spike), and M4 **Phase 1** (Typer CLI + `--engine legacy`
removal) have all shipped — see CHANGELOG `[1.0.0]`. M4 **Phase 2** (PySide6 +
async Qt GUI) is shipped too: the `riana gui` Integrate and Model tabs both run
end-to-end. The post-M4 planning round is done — the next batch of work is the
five-track "Post-M4 roadmap" at the end of this section (it supersedes the old
record-only M5–M8).

### M1 — Stabilization → 0.9.0 — DONE

Smallest viable bugfix release on the existing layout. Fixes B1–B13,
modernizes packaging, adds CI. Behavior-preserving except for the
documented mass-tolerance semantic correction. See `CHANGELOG.md`.

### M2 — Calibration test dataset (acquisition + dual-ID processing) — DONE

The per-cell-line D₂O + ¹⁸O calibration **datasets, the cm-drop50 decision, and the
current-defaults numbers** (OOB R², curation yields, m0 recovery RMSE for ac16 / ipsc /
cm / ac16-¹⁸O on `runs/calib_*_v1`) now live in
`reports/2026-06-23_calibration_benchmark_harness.md` (§ "Datasets & refreshed current
state"). The standing harness (`run_calibration_benchmark.py`) shipped; the NB87a
coefficient-stability approach + the M2→M3 metric findings are recorded there and in
`CHANGELOG.md`. **Open:** none (the curation-gate revisit it fed is tracked under Track C).

### M3 — Aggressive restructure → 1.0.0 — DONE

Delivered Weeks 0–4 + a pre-M4 peak-detection spike; itemized record in
`CHANGELOG.md [1.0.0]`. Outcome: the new `core/` + `algorithms/` + `io/` package
layout with typed records and frozen configs (**layout is documented in the
README**), streaming mzML + dual Percolator/mzTab intake, the integration rewrite
with apex/consensus peak detection (the **apex-narrow window is now the default**,
benchmark-gated on the calibration series; 0.9.0 reproducible via `--peak-rt ms2
--integration-half-width 1.0`), and the fitting rewrite (IsoSpec forward/solve FS +
bootstrap CIs — the §2b fixes). Peak-detection rationale incl. the discarded
`linear` baseline is in §2c; full record in the `m3-peak-detection-spike` memory.
**Open refinements** carried into the tracks below: a cross-proportion-stable peak
picker (Track B) and robust baseline subtraction (§2c / Known Limitations).

### M4 — Qt + CLI rewrite, legacy removal — DONE

Both phases shipped (see `CHANGELOG.md [1.0.0]`). **Phase 1 (2026-06-06):** the
Typer `riana/cli.py` replaced argparse; the entire legacy pipeline
(`riana_integrate`/`riana_fit` + shims, `spectra`/`peptides`/`project`) and the
broken Tkinter `riana_ui/` were deleted — the typed pipeline is the only engine
(`--engine` gone; 0.9.0 reproduced via `--peak-rt ms2 --integration-half-width 1.0`).
`riana fit` requires `--coefficients`; `--label` collapsed to `{hw, o18}`; AA/SILAC
fitting dropped. **Phase 2 (2026-06-07):** the PySide6 + `qasync` GUI (`riana/gui/`,
`[gui]` extra) on a process pool, with a Qt-free worker layer that provably matches
the CLI numerics. The GUI now has Integrate / Model / Protein tabs.

### Status & deferred work

**1.0.0** (M1–M4) and the **1.1.0** experimental-science + QoL line are shipped — the
itemized record is in `CHANGELOG.md`. 1.1.0 added the ¹⁸O rewrite + kinetic fit,
mass-defect θ / `fs_ds`, the D₂O coefficient tables (new default `deberneh_2025_rss`),
pyteomics 5.x, the label-taxonomy strings, LC-fraction collapse + winner-fraction MBR,
the Spep curation gate, CLI progress bars, and the GUI changes (SDRF-only, isotopomer
bar chart, fixed hint area). The remaining open work lives in the Tracks below.

**Deferred (named milestones):**
- **1.1.1 — folded into 1.1.0 (shipped):** its two items — GUI display of prior
  fit/rollup/integrate results from a manifest without refitting, and determinate
  fit/rollup progress bars — both shipped inside 1.1.0 (see `CHANGELOG.md`), so the
  milestone is empty.
- **1.2.0 (in progress on the `1.2.0` branch — see `CHANGELOG.md [Unreleased]`):**
  **shipped** — masking-vectorization perf, parallel test suite + `slow` tier, the
  selectable-condition-pair linear-simple Δk, and the **TMT / TMTpro + single-timepoint
  labeling** line (pinned-isotope envelopes, isobaric SDRF collapse, single-timepoint
  auto-detect + `--min-fit-points` + FS rail-drop). **Remaining** — SILAC/dimethyl
  multiplexing + the multi-point FS rail-drop (Track C), deamidation (own side
  project), the adaptive-N_ISO + robust-envelope rework.

The **Locked decisions** below are design invariants (kept for reference); the Tracks
A–E are the open research/engineering clusters.

#### Locked decisions (2026-06-07)

1. **Labeling time is a required SDRF column** — `characteristics[labeling
   time]`, not parsed from `source name` or the filename. Deliberately a
   *characteristic* (sample-intrinsic), not `factor value[time]`: it is the
   kinetic-curve x-axis, and keying on `characteristics[labeling time]` avoids a
   real collision with a drug-treatment time course (where the *treatment* time
   is the genuine `factor value[time]` while the metabolic-labeling duration is a
   separate property). The current `data/timeseries_lve/samplesheet_lve_sdrf.tsv`
   does *not* carry it (time is only in `LVE_d0` / `…time0`), so adopting Riana's
   SDRF convention means adding this column. Open sub-decision: the value/unit
   format (bare number with a documented unit vs. a unit token — LVE is days).
   Riana reads its own documented subset, so it does not depend on generic SDRF
   tooling treating the studied variable as a factor value. The calibration
   mixing series declares its type with the sibling `characteristics[mixing
   proportion]` instead; fit dispatches on which column is present (kinetic
   models vs. the `calibration` 1:1 recovery model — see Track C).
2. **Header-authoritative identity + a manifest.** Integrate freezes the full
   identity into each output's provenance header (a frozen SDRF snapshot →
   reproducible, SDRF-independent at fit time); a stage-aware
   `riana_manifest.tsv` is the index. Fit groups runs from the manifest rather
   than re-reading the SDRF.
3. **DIA-NN parquet intake is a fast-follow** (its own milestone), but
   `PSMRecord` grows `retention_time` now so the DIA RT-prior path is designed
   in, not bolted on.
4. **Protein comparison ships as a side-by-side view first**; in-app
   cross-sample statistics are deferred. The bridge is the linearized
   simple-model fit (`log(1−θ) = −kt`): once that supports a shared-variance
   two-sample linear model (marginal-mean k per sample, Δk test), the stats
   layer becomes reachable.
5. **Per-sample RIA is an SDRF column** (optional `characteristics[precursor
   enrichment]`; resolution order SDRF → global `--ria` default). It is
   load-bearing at the per-timepoint theta-solve — it builds `_get_final_env`
   in the IsoSpec forward model, which is what recovers θ in the first place —
   so it cannot stay a single global scalar once one mzTab spans multiple
   animals.

#### Track A — I/O & run-identity data model (the spine)

**Shipped** (see CHANGELOG): the `RunIdentity` model + SDRF/manifest intake (M6a),
DIA-NN parquet intake (M6b), per-run concurrency + `--resume`, the manifest project
chain (`integrate → fit → rollup`), the scan↔precursor intake guard, and MBR (mzTab/DDA)
with the winner-fraction policy. **Open:** DIA-NN multi-fraction intake (no data yet
— guard + document). The `RunIdentity` model itself is documented authoritatively in
`riana/records.py`'s docstring.

#### Track B — integration and fit fidelity (science research cluster)

**Shipped** (see CHANGELOG): adaptive N_ISO at integrate (opt-in `--iso auto`), the
H4′ mix-then-truncate FS solve, and limited-isotopomer scoring — flat `--fs` plus
**`--fs auto`**, which *already* widens per peptidoform: it keys on the θ=0
natural-abundance envelope width (`init_envelope_width`) and scores iso0-3 below the
threshold, the full captured envelope (e.g. iso0-5) at width ≥ 6 (`FS_AUTO_BASE` /
`FS_AUTO_INIT_W_THRESHOLD` in `core/fitting`). **Open:**
- **Cross-proportion-stable peak picker** (Phase C v2) — the `apex_search_half_width`
  / `consensus` levers for label-invariant boundary stability are untested in
  production (the §2c / M3 carry-over).
- **Robust observed-vs-IsoSpec matcher** — soft-trim / per-channel-SNR / Huber
  weighting for high-channel contamination; the research-grade endgame, deferred.
- **A genuine per-peptide channel optimizer** — `--fs auto` is a binary init-width
  threshold; the open work is an outlier-aware per-peptide channel choice that
  *balances the bias/variance tradeoff* (more channels = more signal but more
  isobaric-contamination risk), rather than the all-or-iso0-3 switch. Gated on the
  Track D animal benchmark.
- **`noise_floor` baseline subtraction** — implemented (`algorithms/baseline.noise_floor`,
  a flat low-quantile floor) but **off by default** (`--baseline none`): it was
  detrimental in every test so far. Needs more exhaustive testing to decide keep vs
  remove (alongside the still-open robust in-window baseline — §2c).
- **Integrate performance — masking vectorization** (follow-up to the shipped MS1
  cache). The per-run **MS1 peak precache shipped** (2026-06-30; CHANGELOG): each MS1
  is decoded once instead of re-decoded for every overlapping per-PSM RT window — ~8×
  faster `integrate` on dense fractionated runs, byte-identical. The residual per-PSM
  cost is the **masking** (`np.abs(mz−target)≤δ` per scan × iso), which is large
  because `use_range=True` windows span the peptide's whole concat scan range. **Open:**
  vectorize it via `np.searchsorted` on the m/z-sorted centroids (O(log n) per channel,
  byte-identical), and/or revisit whether the apex path needs the full-concat window vs
  the narrower `anchor ± extraction_half_width` (a science decision — changes results).
- **MS2 level integration for DIA-NN path** Explore the use of MS2 fragment isotopomer
  information for estimating theta/FS. First check how much the isotopomer envelope is
  truncated in MS2
- **Other rollup options** Explore other linear simple model rollup option than the 
  current inverse-variance weighted average and pooling, e.g., using mixture models to
  account for peptide-level and biological-level variances.

#### Track C — fitting / modeling science

**Shipped** (see CHANGELOG): M5 per-timepoint fraction-new, the fit-model set
(simple/guan/fornasiero/calibration), protein rollup (`riana rollup`), the
`linear simple` cross-condition Δk model, the ¹⁸O rewrite, M7 PTM-aware envelope
(phospho / N-term-Ac / K-ac / Met-Ox), the 2D-LC fraction collapse + mass-merge, and
the within-protein-θ animal benchmark. **Open:**
- **>2-condition Δk — full multi-group comparison.** The **interim selectable pair
  shipped 2026-07-03** (`rollup --test-condition` + `--reference-condition`, GUI
  dropdowns auto-populated from the manifest; see CHANGELOG): the user picks any two
  conditions and the named pair is contrasted from the **joint all-condition fit**,
  so a multi-group project gets a targeted Δk now (both names validated vs the data;
  the joint fit still pools variance over all conditions — a documented caveat).
  **Open:** the full **all-pairwise / Tukey** extension (all pairs + multiplicity
  correction), which builds *additively* on the same joint fit — blocked on a good
  ≥3-condition dataset.
- **Deamidation** (chemical fit-merge) — its own side project: the +0.984 / C13-M+1
  isobaric overlap needs joint envelope + deamidation-proportion modelling.
- **TMT / TMTpro — SHIPPED 2026-07-05** (see `CHANGELOG.md`). NOT the sample-axis
  "multiplexing" it was first framed as: the multiplexed samples' D₂O signatures are
  inseparable at MS1, so TMT is a **chemical fit-merge mod** whose built-in ¹³C/¹⁵N are
  pinned single-isotope pseudo-elements, and RIANA reports the intensity-weighted-
  average turnover of the plex. Includes isobaric SDRF collapse + single-timepoint
  fitting/curation. **Open follow-ups:**
  - **Multi-point FS rail-drop** — single-timepoint fitting drops per-point FS
    rail-hits (at ±`FS_BOUNDS`) before counting fit points/depth; extend the *identical*
    criterion to multi-point fits (widen the `single_timepoint` scope in
    `core.fitting`). Deferred pending **re-validation** — it shifts the established
    multi-timepoint numbers (rail-contaminated peptides the R² gate currently drops
    would be re-fit on their physical points).
  - **True multiplexing (SILAC / dimethyl)** — the genuine sample-axis mods: a mod
    marks a *different sample* (fit separately, combined at rollup as experiments),
    with heavy ¹³C/¹⁵N/²H as pseudo-elements — the **same machinery TMT now uses**.
    SDRF channel→sample mapping TBD; needs Sadygov dimethyl-D₂O reprocessing.
- **Proper demultiplexing** (Beyond initial 1.2.0) - correct the spillover from
light cluster into heavy cluster (e.g., iso_6 from light cluster overlaps with
iso_9 of the SILAC heavy D2O cluster)
- **GG-remnant (UNIMOD:121)** — Tier-2 PTM, the most turnover-relevant, but blocked
  on anti-K-ε-GG enriched D₂O data (none exists).
- **Cross-fraction RT-correlation MBR** — winner-fraction MBR ships the conservative
  policy; following a peptide that drifts fractions across timepoints is deferred.

#### Track E — GUI/UX & responsiveness

**Shipped** (see CHANGELOG): sortable tables + PNG export, GUI `-W/--workers`,
advanced-knob exposure, the Δmass/Δspacing (`fs_ds`) Model-tab overlay, the
isotopomer bar chart, the fixed hint area + tooltip audit, CLI progress bars, and
the SDRF/manifest-only narrowing. **Open:**
- **Surface per-run SDRF sample / fraction info in the Integrate view** — the
  "display prior results from a manifest" and "determinate fit/rollup progress" items
  this list used to carry **both shipped in 1.1.0** (Load-results on all three tabs +
  the Manager-queue determinate bars — see `CHANGELOG.md`). The residual open piece is
  showing each Integrate row's SDRF sample / fraction. (The integrate-side
  cache-vs-re-read-mzML / smoothing-reproducibility question is the separate
  `--save-traces` item below.)
- **Faithful-to-smoothing chromatogram trace** — the GUI re-extracts the raw XIC on
  row-click, so the trace doesn't reflect the S-G smoothing integration applies. Fix
  by re-applying at click, or an optional **`--save-traces`**
  `<stem>_riana_traces.parquet` sidecar (`(concat, isotopomer) → (rt[], intensity[])`)
  that also removes the ~10 s re-read.
- **Package `riana gui` as a standalone app** — a PyInstaller / py2app bundle (a
  `.app` on macOS, an `.exe` on Windows) so non-Python users can launch the GUI
  without a pip install, and macOS gets a real Dock icon / `.app` (the 2026-06-21
  note that an unbundled `python` process can't own the Dock icon). A nearer-term,
  lower-effort alternative to the Electron rewrite below.
- **Explore other GUI frameworks** - Replace PyQT with a modern Electron app.

#### Cross-cutting chores

**Done:** v1.0.0 release + Zenodo DOI, the 2026-06-21 repo-hygiene audit, the
Snakefile retirement, manifest schema-versioning, and the parallel-fit ≡ serial-fit
reproducibility test. **Open:**
- **Test-suite runtime** — mark the heavy benches / integration tests
  (`@pytest.mark.slow`) and adopt `pytest-xdist` (`-n`) so CI parallelizes.
- **`mypy --strict` rollout** — start on `riana/algorithms/` (smallest blast radius).
- **User-facing docs refresh** — the prose docs are post-M3 stale; non-obvious
  contracts now live in docstrings, but a full pass is its own chunk.
- **Benchmark CI smoke tier** — a subsampled fast bench tier to catch science
  regressions (e.g. within-protein-θ drift) that unit tests miss.
- **SDRF as a partial config source** — auto-load `comment[modification parameters]`
  into the fit config (mass tolerance stays a separate Riana parameter).

## 4. Cross-cutting recommendations

All shipped in 1.0.0 (see `CHANGELOG.md`): provenance-header reproducibility, the
frozen `IntegrationConfig`/`FitConfig` single source of truth, `from __future__
import annotations` throughout, the 2026-06-21 repo-hygiene audit, and the
`workflow/Snakefile` retirement (Riana stays a linear `integrate → fit → rollup`
chain glued by the manifest — orchestration-agnostic; see §5).

## 5. What this project / roadmap deliberately does not include

- **AA / SILAC fitting** — dropped from `riana fit` (the buggy `--label 4` `a_0`
  path is gone); Riana fits metabolic-water labelling (D₂O / ¹⁸O) only.
- **A bundled workflow engine** — the `workflow/Snakefile` is retired; SDRF +
  quantms / DIA-NN own search + ID, and Riana is a linear `integrate → fit → rollup`
  chain glued by the manifest (orchestration-agnostic). The SDRF/manifest path is
  the canonical intake; the per-timepoint Percolator path is kept for benchmarks /
  CLI dev and will be gradually deprecated.
- **Complex / Bayesian kinetic models** — nothing beyond the existing three
  (simple / guan / fornasiero); the focus is honest uncertainty on those, not more
  models.
- **Multi-omics integration, a web interface, or a plugin system** — not right for a
  single-developer scientific tool.

## 6. Known limitations

- **Peak fidelity / baseline:** integration uses an apex-centred narrow window (the
  calibration-gated 1.0 default); robust in-window baseline subtraction and a
  cross-proportion-stable peak picker are still open (§2c / Track B).
- **Linear Δk model:** `linear simple` is an ordinary (unweighted) OLS in
  φ = log(1−θ) space — homoscedastic over per-timepoint points whose precision
  actually varies (peptide depth, φ non-linearity); see
  `reports/2026-06-28_lve_atr_d2o.md`.
- **Memory** is bounded — streaming/indexed mzML, one fraction at a time
  (`io/mzml.py`, `IndexedMzML`).
