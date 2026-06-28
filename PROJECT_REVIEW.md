# Riana — Project Review & Roadmap

> **Status (version `1.0.0` in code; `v1.0.0` git tag + Zenodo DOI still pending —
> a release action; branch `m3-rewrite`).** The M1–M4 rewrite is complete
> (typed pipeline, streaming I/O, dual ID intake, apex-default peak integration,
> IsoSpec forward-model FS, Typer CLI, PySide6 GUI). The post-M4 planning round
> (2026-06-07) retired the record-only M5–M8 numbering into **five tracks** (A
> I/O & run-identity, B integration fidelity, C fitting science, D validation
> infra, E GUI/UX) — see §3.
>
> **Shipped since (full detail in `CHANGELOG.md` `[1.0.0]`; one-liners here):**
> - **Track A** — M6a run-identity model (`RunIdentity` + `io/sdrf.py` +
>   `io/manifest.py` + `core/pipeline.py`; `integrate --sdrf` / `fit --manifest`);
>   M6b DIA-NN parquet intake; intake scan↔RT guard; manifest project chain
>   (`integrate → fit → rollup`); **MBR for the mzTab/DDA path** (gated RT-transfer,
>   2026-06-18, feature-complete 06-20).
> - **Track C** — M5 per-timepoint fraction-new (`riana_fit_fractions.txt` with
>   bootstrap PIs); protein rollup (`riana rollup`, `--parsimony {unique,isoform}`,
>   `--method {weighted,pooled}`); **`linear simple` cross-sample Δk model** (φ-space
>   OLS + Δk test + BH); **M7 PTM-aware envelope** (atom-vector `+P`, mods threaded
>   IO→fit, proteoform rollup keys, Met-Ox fit-merge).
> - **Track D** — within-protein-θ benchmark + frozen LVE bench set (surfaced +
>   fixed the quantms filename-prefix scan-scramble and the mass_tol 50→10 ppm
>   centroid-window fix).
> - **Track E** — GUI rewired onto `core/pipeline`; Protein tab + φ-space plots;
>   `-W/--workers` everywhere (`-t/--thread` removed — GIL-bound, no speedup).
>
> **Next (unblocked):** expose hidden integrate knobs as marked *advanced* options;
> Tier-1 PTM mods (K-acetyl, K/R-methylation — cheap, machinery exists); cut a real
> `1.0.0` tag + repo hygiene; user-facing docs refresh. **Blocked/demand-driven:**
> DIA-NN phospho sites (user's variable-phospho rerun), >2-condition Δk (needs ≥3
> conditions), Track B fidelity (gateable by Track D), o18 (needs NB90b table),
> deamidation (side project). See §3 "Handoff — ordered priorities".
>
> Maintainer: Edward Lau. Last reviewed: 2026-06-20.

This document consolidates and supersedes the prior `documentation/` folder
(`PROJECT_EVALUATION.md`, `ROADMAP.md`, `MASS_ACCURACY_SPECIFICATION.md`). The
old contents are folded in below where still load-bearing; the rest is cut.

## 1. Project status

> **Update (2026-06, version `1.0.0`):** the three-step plan below is complete
> through M4 — 0.9.0 stabilization, the calibration dataset, and the 1.0.0
> rewrite (typed pipeline, streaming I/O, peak detection, mzTab intake, and a
> PySide6 GUI) have all shipped. The narrative below is the original pre-rewrite
> evaluation that motivated the plan, kept for context; the as-built record is
> CHANGELOG `[1.0.0]` + §3. The post-M4 planning round is done — the next
> batch of work is organized into five tracks in §3 "Post-M4 roadmap".

Riana is a single-author scientific Python tool for extracting and modeling
isotopomer time-series from MS1 data, used in protein-turnover research. The
core scientific functions — accurate-mass calculation, kinetic models,
fractional-synthesis math, IsoSpec wrapping — are correct and proven. The
surrounding software (data parsing, integration pipeline, Tkinter GUI,
packaging, tests) shows typical scientist-coder accumulation: tight coupling,
weak tests, dead code, brittle string parsing.

Rather than incrementally patch the structure, the plan is:

1. **0.9.0 — Stabilization.** Fix correctness defects in place. Modernize
   packaging and add CI. No restructure.
2. **Calibration dataset** (parallel work). A 9-level D₂O mixing series
   becomes the regression check for every future change.
3. **1.0.0 — Aggressive rewrite.** New package structure, typed records,
   streaming mzML, async pipeline, Skyline-style peak detection, mzTab
   intake, PySide6 GUI. The numerically-sensitive scientific functions are
   **lifted unchanged** rather than rewritten.

The 0.9.0 → 1.0.0 boundary is a deliberate breaking change. The user base is
small and scientific software is allowed to evolve.

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

### Post-M4 roadmap (planning round, 2026-06-07)

M1–M4 shipped the rewrite (typed pipeline, streaming I/O, peak detection,
mzTab intake, Typer CLI, PySide6 GUI). This round retires the record-only
M5–M8 numbering and the scattered deferrals (Phase C v2, the animal
within-protein benchmark, the `--fs` subset, the o18 rewrite) and reorganizes
everything — plus substantial new feature work surfaced in the session — into
**five tracks** with named milestones and a sequence.

**The through-line.** Nearly every I/O pain point is downstream of one missing
abstraction: a **run-identity data model**

```
RunIdentity = (experiment, source/sample, condition/group,
               biological_replicate, technical_replicate, fraction,
               labeling_time, acquisition, precursor_enrichment)
```

attached to every PSM at intake and carried — *header-authoritatively* — to
fit and protein rollup. Today identity is positional and conventional
(`file_idx` as a sort index, timepoint encoded in the filename,
one-Percolator-file-per-sample by habit); that is the fragility. The SDRF is the
carrier for both acquisition modes — DDA via quantms and DIA via quantms-diann,
which emits its own SDRF (a near-identical variant, a few columns different) —
so the identity model is shared; only the PSM/quant file format differs (mzTab
vs DIA-NN parquet).

#### Handoff — 1.1.0 planning round (as of 2026-06-24)

**1.0.0 RELEASED 2026-06-24** — GitHub release + Zenodo DOI; `v1.0.0` annotated tag at
`master` `c2c3d7b`. Shipped across the whole 1.0 line (now in CHANGELOG; cleared from
this list): the M3 rewrite, M7 v1 + Met-Ox + K-acetyl proteoform keys, the `linear
simple` cross-sample Δk milestone, MBR (mzTab/DDA), advanced-knob CLI+GUI exposure,
the docs refresh, and the **Track B N_ISO line** — adaptive `--iso auto`, the H4′ FS
solve, `--fs` limited-isotopomer scoring, the `run_calibration_benchmark.py` driver,
**`--fs auto`** per-peptide widening (keyed on the RIA-invariant natural-abundance
*init width*; crossover init-width 6 cross-line; `core.fitting.FS_AUTO_*`), single-int
`--fs`/`--iso`, **GUI parity** for those knobs, **DIA-NN phospho proteoforms**
(`io.diann` `Protein.Sites`→`pS###`, conf≥0.75), and the release-hardening dependency
fixes (**`lxml`** declared, **`pyteomics<5`** pinned — a fresh install reads mzML).

**`1.1.0` branch cut 2026-06-24** (off `master` `c2c3d7b`; version bumped to `1.1.0`).
This is the experimental-science line. A fresh planning round is pending; the agreed
ordered items + the user's framing (2026-06-24):

1. **Δmass-over-time QC (GUI) — FIRST; lays the groundwork for mass-defect→θ.** Plot
   per-peptide observed Δmass (from the `iso{N}_ppm_error` / orthogonal-θ substrate
   adaptive N_ISO already emits) vs labeling time. Build this diagnostic view first,
   then the **mass-defect → θ estimator** on the *same* substrate (θ from the labelled
   envelope's mass shift, orthogonal to the abundance-ratio FS). Track E → Track B.
2. **Pyteomics 5.x upgrade — SECOND.** Unpin `pyteomics<5` (the 1.0.0 release fix) and
   add `psims`; 5.0 brings the PSI-MS CV resolver and **potential per-mzML
   multithreading** worth evaluating for integrate. Verify mzML parity vs 4.7.5 first.
3. **o18 rewrite — blocked on data (re-run in progress, ~this week).** `fit --label
   o18` errors until: (a) the o18 calibration data is **re-run through the mzTab path**,
   and (b) a **new coefficient-type structure** is set up and trained with the
   **reverse model** — the same calibration→coefficients pipeline that produced the
   ac16/cm/ipsc D₂O tables, but for ¹⁸O (the NB90b frozen-coefficient artifact).
   **Cleanup folded into this rewrite (user 2026-06-24):** collapse the *internal*
   `get_peptide_distribution`/`fsynthesis` integer `label` {1 = D₂O in-vivo, 2 = D₂O
   in-vitro, 3 = ¹⁸O} to clear string labels `"D2O"` / `"O18"` (+ future double-water,
   post-1.2.0). The 1-vs-2 in-vivo/in-vitro split is obsolete — D₂O cell-specificity is
   the **coefficient table** (Commerford mammalian → ac16/ipsc/cm), not the label, just
   as the CLI already collapsed `--label` to `hw` (consider aligning `hw`→`D2O`). NB the
   `label == 3` branch in `isotope_dist.get_peptide_distribution` is currently a **stub**
   — it extends `atom_count` + isotope mass but NOT `isotope_probability_list` (and reuses
   the ¹⁶O mass), so the envelope can't run for ¹⁸O until the rewrite completes it.
4. **TMT — most invasive; deferred. Sets up a THIRD mod-handling type: multiplexing.**
   The mod taxonomy: **chemical** (Met-Ox, CAM) folds at the *peptide* level
   (`_fit_key` merge); **biological** (phospho, K-acetyl) gets a distinct *proteoform*
   rollup key; **multiplexing** (TMT, dimethyl, SILAC) tags **sample identity** — fit
   each channel's peptides *separately*, then combine at rollup **as if from different
   experiments — which they are** (each channel = a different biological sample). That
   needs an **SDRF channel→sample mapping** (TMT has a documented way to record per-
   channel biological samples; the dimethyl/SILAC SDRF syntax is TBD — to look into).
   **NB: TMT is a *chemical* mod as far as D₂O is concerned — there is one D₂O peptide
   cluster** (isobaric across channels); it multiplexes *other* information, not the D₂O
   envelope. Envelope mechanism: inject TMT's heavy ¹³C/¹⁵N as **pseudo-elements** the
   IsoSpec envelope tracks (like the D₂O/¹⁸O labels), *not* a mass override. Blocked on
   data: user has **TMT-D₂O** to process; **dimethyl-D₂O** needs reprocessing public
   **Sadygov** data. See `track_c_tmt_chemical_mod`.

Smaller / deferred items are in the lists below (deamidation side-project, non-universal
cysteine CAM, K/R-methylation Tier 2, MBR-on-calibration RT correction, chromatogram-
smoothing trace, matplotlib export, modernize look).

**Deferred GUI/UX (Track E) — sortable tables + graph export DONE 2026-06-21.**
Remaining: faithful-to-smoothing chromatogram trace, the **Δmass-over-time QC** (couple
it to adaptive N_ISO's `iso{N}_ppm_error` / orthogonal-θ substrate), matplotlib
static/SVG export, fit progress + real parallelism, modernize look. See Track E.

**Demand-driven / blocked on data or a decision (do when unblocked):**
- **K/R-methylation — relegated to Tier 2 (2026-06-21).** `UNIMOD:34/36/37`,
  biological keys (`meK###`/`me2R###`); composition is already in `mod_atoms` so it
  stays cheap, but demand is low and it is no longer a near-term priority. (K-acetyl
  already shipped.)
- **DIA-NN phospho proteoform sites — DONE 2026-06-23.** `io.diann` maps the
  localized `Protein.Sites` → biological-mod proteoform suffix (`P35486_pS293`),
  gated on `PTM.Site.Confidence` (default ≥ 0.75); Carbamidomethyl-C and Met-Ox-M
  sites are excluded (only `BIOLOGICAL_MODS` define a proteoform), key format
  byte-identical to the DDA/mzTab path. Validated on the cardiac variable-phospho
  series (577 sited proteoforms). The two real searches are complementary
  (`diann_results` = Met-Ox only, `diann_results_mods` = phospho only — quantms
  library generation ran out of memory with both), so each path is exercised by
  one dataset.
- **Deamidation** — its own **side project** (the +0.984 / C13-M+1 isobaric overlap
  needs joint envelope + deamidation-proportion modeling). Not started.
- **Non-universal cysteine alkylation (far-future, deamidation-tier priority).** CAM
  (UNIMOD:4) is today a fixed mod on *every* cysteine via the `mass_calc` `iaa=True`
  flag (composition already the unified `mod_atoms[4]`), which is correct for standard
  IAA data and is the only common denominator (Crux/Percolator never *declares* CAM —
  its sequences are bare; mzTab/DIA-NN declare but strip it). Revisit only with data
  that needs *partial* alkylation / alternative reagents (NEM) / free cysteines: then
  CAM becomes a per-cysteine **variable** mod threaded via `count_atoms(mods=…)`, not
  the blanket flag. Tokenizing it before that would just bloat every cysteine peptide's
  `concat` identity for zero information. See the `_count_residue_atoms` comment.
- **MBR-on-calibration RT correction (NOT a re-search — corrected 2026-06-23).** The
  mzml re-search **landed and does NOT fix the offset**: measured median |mzML-scan-RT −
  mzTab-RT| is *identical* for the `.raw`- and `.mzML`-searched mzTabs (ac16 0% file 2.16
  vs 2.15 min; 100% file 0.56 vs 0.58), i.e. the offset is **OpenMS ProteomicsLFQ RT
  alignment**, applied by both pipelines, not a search-input artifact (and it is *larger*
  at low θ, so not a high-θ effect). To MBR the calibration, subtract the measured
  per-run RT offset before `resolve_rt_anchored_scans` (or re-anchor MBR in scan space) —
  then re-run `bench_mbr_apex_knobs.py` without `--no-rt-check`. For the Percolator-path
  0.9.0 bench, irrelevant. See report `2026-06-23_adaptive_niso_limited_isotopomer.md` §4
  + memory `mbr_rt_axis_dependency`.
- **GG-remnant (tier 2 PTM)** — most turnover-relevant, but blocked on anti-K-ε-GG
  enriched D₂O data (none exists).
- **Heads-up:** the user will rename the run `atf6small` → `atf6_lve`.

**Parked / blocked (do when unblocked, not ahead of the above):**
- **>2-condition pairwise Δk contrasts** — `fit_linear_deltak` handles exactly two
  qualifying conditions today; the pairwise/Tukey extension is **blocked on a good
  ≥3-condition test dataset** (user to provide). Forward-compatible (per-condition
  k for any N; Δk left NaN when ≠2).
- **Track B integration fidelity** + the `calibration` 4th model + guan/fornasiero
  validation — demand-driven, gateable by the Track D within-protein-θ bench.
- **o18 rewrite** — needs the NB90b frozen coefficient table; `fit --label o18`
  errors until then.
- **GUI framework feasibility (Electron/Tauri/web)** — long-term decision doc, only
  if Qt becomes a concrete graphing ceiling.

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

- **M6a — identity model + intake refactor. DONE (2026-06-07, branch
  `m3-rewrite`).** Shipped as `RunIdentity` (`records.py`) + `io/sdrf.py` +
  `io/manifest.py` (schema-versioned, stage-aware) + `core/pipeline.py`
  (`integrate_project` one-file-per-run with identity-in-header + manifest;
  `recombine_for_fit`/`fit_project` grouping curves by `(experiment,
  condition)`), wired as `integrate --sdrf` / `fit --manifest`. Curve grouping
  was user-confirmed: one curve per condition, `characteristics[biological
  replicate]` are independent replicate points (not separate curves). Percolator
  is the demoted single-mzML tier. Variable mods are parsed/exposed but not yet
  threaded into the envelope (that stays M7). A parity test pins the manifest
  path to the legacy path. **Deferred to follow-ups:** DIA-NN intake (M6b),
  bounded file-parallelism (below), GUI rewiring onto `core/pipeline.py`
  (Track E), threading SDRF mods into integration. Original spec below.
- **M6a (original spec).** `io/sdrf.py` + a documented
  Riana-read column subset (`comment[data file]` → mzML join key; `source
  name`; `characteristics[biological replicate]`; `comment[technical
  replicate]`; `comment[fraction identifier]`; `characteristics[labeling time]`
  (required for turnover, the kinetic-curve x-axis) OR `characteristics[mixing
  proportion]` (calibration fixtures) — the independent-variable column declares
  experiment type and fit dispatches on it; `factor value[...]` → condition/group
  for side-by-side + future cross-group stats (key on whichever `factor
  value[...]` column(s) exist — `factor value[condition]` the canonical case,
  so `factor value[disease]`/`[genotype]` work too);
  `comment[proteomics data acquisition method]` → DDA/DIA;
  `characteristics[precursor enrichment]` → RIA (optional)); validation with
  clear errors; everything else quantms fills is ignored. Identity +
  `retention_time` on `PSMRecord`. mzTab becomes the primary path
  (`read_mztab(path, sample_map)` where `sample_map` comes from the SDRF);
  Percolator is demoted to a single-mzML testing/legacy tier (don't rescue
  multi-fraction Percolator). **Output: one `<mzml_stem>_riana.txt` per run**
  (no timepoint-name collisions), full identity + RIA frozen into the header.
  **Stage-aware `riana_manifest.tsv`** (`stage = integrate | fit | protein`) as
  the project index. **Fit-time recombination** — group runs by `(experiment,
  sample, bio_rep)`, merge fractions of the same `(sample, rep, timepoint)` at
  peptide level, order by time — lands in a shared **`core/pipeline.py`** (the
  extraction the GUI "mirror, don't refactor" note flagged), consumed by both
  CLI and GUI. Read `comment[modification parameters]` from the SDRF so variable
  mods need not be set on the CLI; **and read `comment[precursor mass tolerance]`
  as the integration window** (CLI `--mass_tol` overrides; default 10 ppm).
  *(Reverses the original rec4: it had said NOT to inherit the search tolerance,
  using a deliberately wider ~50 ppm "signal-capture" window. That was a
  profile-vs-centroid category error — on **centroid** mzML, ~50 ppm imports
  co-eluting interference into the heavy isotopomer channels. Verified 2026-06-09
  on LVE: tightening 50→10 ppm lifted R²med 0.69→0.87, ≥0.8 41→61%, within-protein
  θ 0.16→0.10, out-of-range θ 27%→8%; all current mzMLs are MS1 centroid. The
  search tolerance IS the right integration window for centroid data. A
  profile-mode mzML triggers an intake warning.)*
- **M6b — DIA-NN parquet intake. SHIPPED (2026-06-11, branch `m3-rewrite`).**
  `io/diann.py` reads the DIA-NN `report.parquet` (validated on **2.5.0** with the
  `diann` parquet output, ≥ 2.2.0 as emitted by quantms-diann); `plan_integration`
  dispatches the reader on `SdrfTable.acquisition` (DIA → parquet, else mzTab), so
  the same `io/sdrf.py` spine and both CLI + GUI surfaces get DIA for free — only
  the parquet reader is new. DIA has **no MS2 anchor**, so the reader emits
  `scan = -1` + carries DIA-NN's apex `RT`, and `core/integration.resolve_rt_anchored_scans`
  maps it to the nearest MS1 scan in the mzML at integrate time (+ an RT-in-bounds
  wrong-mzML check); the scan↔RT scramble guard is skipped for DIA (circular
  there). Riana still extracts the MS1 isotopologues itself — DIA-NN is "just
  another ID + RT source." **Variable-mod peptidoforms are dropped** (the M7
  caveat — they'd integrate at the unmodified m/z; ~1.8% on the cardiac set).
  Needs `pyarrow` (the `[dia]` extra, lazily imported). Validated end-to-end
  (`integrate → fit → rollup`) on the cardiac in-vivo DIA set under
  `data/timeseries_dia` → `runs/lve_dia`: **98.5% m0 coverage, ~1–2 ppm mass
  accuracy** (extraction is on-target), fit **R²med ~0.42 / ≥0.8 ~22%** (lower
  than the DDA LVE ~0.87 / ~61% — DIA's wide-window MS1 imports co-eluting
  interference into m1–m5, and the curve is 3 timepoints at RIA 4.6%), but
  **k_deg median 0.08/day** is biologically sound. *Open follow-ups:* a stronger
  DIA-specific QC than RT-in-bounds; if DIA turnover matters, an MS1-interference-
  robust envelope matcher (overlaps the Track B adaptive-N_ISO matcher).
- **Integrate concurrency + GUI rewiring (Track E, SHIPPED 2026-06-10).**
  `integrate_project` was split into `plan_integration` (SDRF+mzTab → per-run
  `RunTask`s) / `_integrate_results` (yields each frame as it finishes) /
  `finalize_run`. **Crash-resilient (2026-06-11):** the main process writes each
  run's `_riana.txt` + appends its manifest row *as that run completes* (workers
  only integrate), so an interruption keeps finished runs; `integrate --resume`
  skips runs already in the manifest at the current `config_hash`. The GUI
  Integrate tab now drives the *same* plan + per-run unit over its **own** shared
  pool via `asyncio.as_completed` bounded by a *Workers* spinbox (one mzML per
  concurrent run) — cross-file parallelism with **no nested process pools** —
  and gained an **SDRF** field that routes the search-ID file through
  `core/pipeline` (identity-stamped `<stem>_riana.txt` per run + manifest), on
  top of the per-run *Threads*. The GUI Model tab gained a **Manifest** field
  (`fit_project` via `tasks.run_fit_manifest`) and now also writes the M5
  `riana_fit_fractions.txt`. So integrate/fit go through `core/pipeline` on both
  surfaces. **Rollup threading shipped** (2026-06-10): the per-protein refit is
  dispatched over `--thread` / the Protein-tab spinbox with per-protein RNG
  streams (`_group_rng`), so the result is identical regardless of thread count.
  *Still TODO:* the SDRF mass-tolerance resolution the CLI does (the GUI uses the
  spinbox value explicitly).
- **Manifest project chain (SHIPPED 2026-06-10).** `fit --manifest` /
  `rollup --manifest` now root outputs at the **manifest's folder** (the project
  dir; `-o` ignored there with a warning) and write back their `stage="fit"` /
  `stage="rollup"` rows (`record_stage_rows` + `fit_outputs_from_manifest` in
  `core/pipeline.py`), so one `--manifest` drives `integrate → fit → rollup` and
  the manifest indexes every stage. `rollup`'s `FIT_DIR` argument is now optional
  (XOR `--manifest`). The aggregate fit/protein rows carry a coarse
  experiment-level identity (`aggregate_identity`) — indexing/provenance only.
  (Closes the gap where `fit` never updated the manifest despite the M6a plan.)

> **⚠️ quantms input gotcha — mzML filenames MUST NOT be prefixes of one another
> (verified 2026-06-08, quantms/OpenMS ~1.7.0).** ConsensusID/ProteomicsLFQ
> matches mzML/spectra by a **filename-prefix** rule, so when one basename is a
> prefix of another (`…_time1` vs `…_time10`/`…_time15`; `…_time3` vs `…_time30`)
> the mzTab `spectra_ref` **scan is pulled from the sibling file** — the scan↔file
> association silently scrambles for every entangled run. On the LVE series this
> wrecked integration on 8 of 12 runs (their `spectra_ref` scans indexed a
> sibling's mzML; e.g. `time30`'s scans → `time3`'s mzML to 0.13 min). Proof was
> single-variable: **renaming to zero-padded non-prefix names (`time00…time30`),
> same 2023 conversion, restored scan-match to ~90% on all 12.** Also previously
> seen as the reason Sage couldn't be added to ConsensusID. **Mitigations:**
> (1) zero-pad / otherwise de-prefix mzML basenames before a quantms run;
> (2) **anchor mzTab integration on `spectra_ref` scan, NOT the reported
> `retention_time`** — RT carries an alignment-frame offset (~0.1 min on clean
> runs, but the mzTab RT is OpenMS-aligned, not raw) so scan is the correct,
> confound-free key *once filenames are clean*; (3) **Track A intake guard
> (SHIPPED 2026-06-10):** `integrate_run` reconciles, per run, each PSM's
> `spectra_ref scan → mzML RT` against the mzTab `retention_time` and raises
> `DataError` when the **per-run median** offset exceeds `scan_rt_tol_min`
> (default **2.0 min**). The median (not per-PSM) gate is robust to the ~10% of
> PSMs that legitimately mismatch and to the run-dependent ProteomicsLFQ
> alignment offset — measured at ≤~0.9 min on the de-prefixed LVE runs (time00
> median 0.92, time03 0.22, time30 0.36; all 100% within 2 min), vs tens of
> minutes on a scrambled run (~25× separation). Errors loudly with the
> de-prefix/zero-pad remedy in the message; `--no-rt-check` (config
> `check_scan_rt=False`) overrides. No-ops on the Percolator path (no
> `retention_time`). **DIA caveat:** this is scan-anchored, so it does not
> transfer to M6b — DIA-NN reports an inferred per-run apex RT with no MS2-scan
> anchor to reconcile against; the DIA analog (reported RT within mzML bounds +
> run-column matches the paired mzML) is a separate, weaker check for M6b.

- **Match-between-runs (MBR) for the mzTab/DDA path — SHIPPED 2026-06-18,
  feature-complete 2026-06-20.** Full record in CHANGELOG `[1.0.0]` + the design
  report `reports/2026-06-17_mbr_v1_design.md` (commits `6b68991`→`47bb806`);
  benches `bench_missingness` / `bench_rt_alignment` / `bench_mbr_quality` /
  `bench_mbr_ab`. The load-bearing findings worth keeping in the roadmap:
  - **DDA is where MBR pays, not DIA.** DDA turnover curves are badly gappy (only
    8–15% of precursors span all 12 timepoints, ~50% of slots empty, t0 the sharpest
    loss at 70% LVE); DIA-NN's internal propagation already fills DIA. The gap
    survives a matched-geometry control (DIA 60% vs DDA 32–44% at 3 timepoints), so
    it is acquisition + ID-pipeline, not curve length.
  - **quantms aligns RT but imperfectly** — co-IDs sit within ~4–9 s median |ΔRT|
    but run-specific residuals reach 15–25 s (> the ±9 s window), near-constant
    per-run → the shipped fix is a robust per-run-pair RT offset (not a global
    realign), on top of M6b's `resolve_rt_anchored_scans` substrate.
  - **Verdict:** gated MBR is neutral at strict R²>0.95 and net-positive at the
    in-vivo gates (+180 at R²>0.8) with no pollution → shipped uncapped. Eval-only
    yield-for-consistency tradeoff is reasonable (+37 proteins at R²>0.8 for ~+3.5%
    within-protein scatter).
  - **Open:** the high-label calibration concern was **retracted** as an
    mzTab↔mzML RT-axis artifact (`.raw`-searched mzTab + `--no-rt-check` bypass) —
    maintainer TODO is to re-search on the exact `.mzML` and re-run without
    `--no-rt-check`. **Parked:** sub-threshold *rescue* tier / MBR-FDR is its own
    future study (mzTab is 1%-FDR pre-filtered: 179 PSMs in (0.01,0.02], 0 above).

#### Track B — integration fidelity (research cluster)

Plan these together — they are complementary integration-fidelity levers, and
the high-isotopomer channels are exactly where peak detection, baseline, and
N_ISO all interact:

- **Phase C v2 — cross-proportion-stable peak picker.** The one untested lever
  from the M3 spike: tight `apex_search_half_width` (~0.15) and/or
  `consensus_apex` for boundaries stable across labelling proportions (the
  boundary-stability finding).
- **Adaptive N_ISO via IsoSpec-at-integrate — the NEXT integration-fidelity item
  (prioritized 2026-06-21, ahead of TMT/SILAC/dimethyl).** Today the forward model
  runs **only at fit**; `integrate` extracts a *fixed* `--iso` channel set at
  analytic m/z (`calculate_ion_mz`: residue+mod masses, `base + iso·mass_diff/z`).
  **Promote the envelope to integrate time:** per peptidoform run the IsoSpec
  **init (0%)** and **final (RIA%)** envelopes and set **N_ISO per peptide** from
  the union of significant channels — a short peptide gets m0–m3, a long/
  heavily-labeled one m0–m8 — instead of the one-size `--iso 0…5/6`. Three things
  fall out:
  1. **`iso0` IS the precursor m0 mass.** One computation, inherently consistent
     with the envelope, so the separate mod-mass path collapses into it. TMT /
     heavy-dimethyl then need **no** analytic `unimod_mass + pinned-isotope` mass
     extension — the pinned-isotope *envelope* support (`mod_atoms` light +
     `mod_fixed_isotopes`, appended in `get_peptide_distribution`) is still the
     shared prerequisite, but the **mass comes from iso0**. (This RESOLVES the TMT
     A/B in favour of B — do this first, get the mass for free; see the TMT box in
     Track C.)
  2. **Adopt the averaged-isotopolog per-isotopomer accurate mass** as the
     extraction target (the `use_nominal_masses=True` weighted-average mass per
     nominal bin), replacing today's `base + iso·1.00335` spacing. **Accepted
     behavior change** (user 2026-06-21): it is *more* accurate — the true centroid
     of each isotopomer cluster, accounting for mixed C/H/N/O/S/²H contributions —
     and is distinct from the MS literature's "nominal mass" (averaged *across*
     isotopomers).
  3. **Orthogonal mass-defect θ** (subsumes the old "dual-mode FS" item). Because
     the ²H−¹H mass defect differs from ¹³C−¹²C, as the profile mixes with more D
     each isotopomer's **accurate mass shifts measurably from its θ=0 position**.
     That shift is an FS estimator *orthogonal* to the abundance-ratio one — solve
     θ both ways and cross-validate; for low-abundance peptides mass accuracy can
     beat spectral accuracy. Integrate already emits `iso{N}_obs_mz` /
     `iso{N}_ppm_error`, so the substrate exists (the near-term QC half is the
     Track E Δmass-over-time display).
  - **Choose limited isotopomers at FIT, not integrate** (subsumes the reserved
    `--fs` channel-subset SSE). Integrate **wide** (adaptive N_ISO captures the
    full envelope) and let `fit` solve θ over a **chosen subset** (e.g. iso0+iso1
    only) to dodge co-eluting contaminants in the high channels. Wide capture,
    narrow scoring — you cannot subset at fit what you did not integrate.
  - **Robust observed-vs-IsoSpec matcher** (the open research question): more
    channels means more contamination risk, so pair adaptive N_ISO with a
    contaminant-robust matcher (Huber / soft-trim / per-channel-SNR weighting) so
    the extra channels don't import co-eluting isobars.
  - **⚠️ HARD CONSTRAINT — normalize the FULL IsoSpec envelope BEFORE truncating
    to N_ISO** (the "H4′" finding, verified 2026-06-08). `solve_fs_d2o` today
    normalizes init/final *separately over the truncated channels* then mixes:
    `(1−f)·norm(init_trunc) + f·norm(final_trunc)`. That equals the physical
    `norm_trunc[(1−f)·init_full + f·final_full]` **only when init and final carry
    the same fraction of mass inside the window.** The labeled envelope spills
    past the window (more, the higher the RIA / Spep), so normalize-then-mix drops
    the θ-dependent denominator and FS↔θ goes **non-linear** — the math won't work.
    Negligible at LVE's RIA 4.6% (in-window mass fraction ΔS ≤ 0.04 ⇒
    current == correct to 3 dp), which is why it is *deferred*; but it becomes
    load-bearing the moment adaptive N_ISO widens/varies the window, or at high
    RIA / high Spep / ¹⁸O. **So build the mix-then-normalize FS solve as part of
    this item** (mix full-envelope *abundances*, THEN truncate + normalize), not
    after. The same trap kills a naive iso0/iso1-ratio interpolation — you must
    mix the abundances, then take the ratio.
  - **Cost note:** integrate gains per-peptidoform IsoSpec (today fit-only), but it
    is cached by `(sequence, mods)` and the distinct-peptidoform count ≪ PSM count,
    so it is bounded and parallel-safe under `-W`. The `use_nominal_masses` envelope
    cache (`algorithms.isotope_dist`) already exists.
  - **RESULT (SHIPPED 2026-06-23; report `2026-06-23_adaptive_niso_limited_isotopomer.md`).**
    Built B0–B4 and ran the 2×2 capture×scoring control on the calibration mixing series
    (ac16/cm/ipsc, ground-truth |θ−f|). **The science win is the fit-side scoring choice,
    not adaptive capture.** `--fs 0 1 2 3` (iso0-3) tightens recovery on all three lines
    (within-±0.05 +1.8–2.8 pp, lower IQR/bias) and is a *pure improvement over 0.9.0*
    (`v1.0.0/all == v0.9.0/all` exactly — byte-faithful port — no regression). Adaptive
    capture is neutral-to-slightly-negative *and* ~3–5× slower, with `fix·iso0-3 ≈
    adapt·iso0-3` everywhere. **Shipped:** `riana fit --fs` (`FitConfig.score_channels`,
    H4′ mix-then-normalize in `solve_fs_d2o`, run-level + per-peptide guards);
    `riana integrate --iso auto` stays **opt-in** (its live justifications are UX —
    parameter-free integrate — and the TMT precursor-mass-from-`iso0` path, not N_ISO
    accuracy).
  - **NEXT — per-peptide `--fs` (parked, designed; brainstorm in the report's Future
    Work).** The flat global `--fs iso0-3` over-truncates the genuinely-wide-envelope
    peptides (N_ISO ≥ 12: +0.08 recovery error — crossover at ≈ N_ISO 11). A **per-peptide
    `score_channels = f(N_ISO)`** (iso0-3 for short/medium, wider for very long peptides)
    captures Currie et al.'s site-keyed heuristic (<15 / 15-35 / >35 D-sites → higher
    channels) in the whole-cluster-RMSE paradigm. Two notes: (i) **extend, don't shift** —
    unlike Currie's 2-channel ratio (which drops a suppressed iso0 at high Spep), the RMSE
    abundance-down-weights iso0 automatically, so keep it; (ii) key on **N_ISO** (already
    computed by `adaptive_channel_masses` — this is where `--iso auto` finally earns its
    keep, as the per-peptide scoring-width *supplier*) or the cheaper Spep proxy. The
    research-grade endgame is the original **soft robust matcher** (per-channel SNR/Huber
    weight) rather than a hard per-peptide cutoff.

Gated by both the mixing-series benchmarks and the new animal benchmark
(Track D).

#### Track C — fitting / modeling science

- **M5 — persist per-timepoint fraction-new (SHIPPED 2026-06-10; CHANGELOG).**
  `riana fit` writes `riana_fit_fractions.txt` (long-format `(concat,
  biological_replicate, labeling_time)` + `fs` and prediction-interval bounds
  `fs_lower`/`fs_upper`). Load-bearing design: the bounds are a **prediction
  interval** from a unified residual bootstrap (`model(t_i; k*) + resampled
  residual`), so they capture each peptide's measurement scatter — the quantity
  the protein rollup inverse-variance-weights. **This is the substrate the rollup
  consumes.**
- **2D-LC / technical-replicate fraction collapse (NEW 2026-06-20; collapse policy +
  mass-merge BUILT 2026-06-27).** The depth gate counts **distinct labeling
  timepoints**, so multi-file multiplicity at one timepoint never inflated
  qualification. The collapse itself runs on the **manifest fit path**:
  `recombine_for_fit` → `_merge_fractions` collapses fractions of one `(concat =
  peptidoform+charge, biological_replicate, labeling_time)` into a single point —
  combining intensities *before* one FS is solved (never averaging per-fraction FS).
  Now configurable via **`fit --fraction-collapse sum|anchor`** (sum each `isoN`, the
  default, vs keep the single highest-total-intensity fraction), and the merge now
  **intensity-weights the per-channel mass/QC columns** (`iso{N}_obs_mz` /
  `_ppm_error` / `apex_snr`) instead of taking the first fraction's — a real `fs_ds`
  bug now fixed. The explicit-files path (no `--manifest`) does **not** collapse and
  now warns. Validated against the **real fractionated iPSC D₂O SDRF**
  (`data/timeseries_lauren_9_ipsc_d2o`, 192 files = 12 tp × 2 biorep × 8 fractions)
  through `read_sdrf`; **end-to-end fit validation awaits the mzML/quantms search**.
  MBR is now **winner-fraction-restricted**: a peptide's holes are filled only in the
  single LC fraction where it has the most IDs (ties → best q, then fraction number;
  `mbr._winner_fractions`), the conservative minimal policy. **Still open:**
  cross-fraction RT-correlation MBR (following a peptide that drifts fractions —
  deferred by the maintainer as non-minimal); DIA-NN multi-fraction intake (no data).
  Memory `track_c_fraction_collapse_gap`.
- **Fit-model set + the calibration model.** The kinetic models `{simple, guan,
  fornasiero}` are *all wired end-to-end* already (models math lifted unchanged;
  `_MODELS` dispatch → `curve_fit`; CLI `--model`; `FitConfig.model` validation;
  GUI combo + k_p/k_r/r_p spinboxes) — they fit `k_deg` with `k_p/k_r/r_p` held
  *fixed* at config values. **Caveat (verified 2026-06-07):** only `simple` is
  validated — guan/fornasiero have **no test or benchmark** (every gate, incl.
  the −0.5 bias close, ran on `simple`) and their `k_p/k_r/r_p` defaults
  (0.5/0.05/10.0) are placeholders. Add **`calibration` as a 4th model** — a 1:1
  line fit of observed FS vs. known mixing proportion (slope ≈ 1, bias, R² = the
  recovery metric; the productized "option b" recovery mode), dispatched
  identically and reusing the GUI curve view. Its x-axis is `characteristics[mixing
  proportion]`, so **fit dispatches on the experiment-type SDRF column**:
  `characteristics[labeling time]` → kinetic models; `characteristics[mixing
  proportion]` → `calibration`. **Future science (as-needed):** validate
  guan/fornasiero against the in-vivo animal labeling data (hard to model) and
  decide whether the precursor parameters (k_p / k_r / r_p) are *fitted* or
  *supplied* — both are fixed today, so meaningful two-compartment use needs real
  precursor priors. Unlikely near-term.
- **Protein rollup (SHIPPED 2026-06-10; CHANGELOG).** `riana rollup`
  (`core/protein.py`) writes `riana_rollup_proteins.txt` + `_fractions.txt`
  grouped by `(experiment, condition, protein)`, with `--parsimony
  {unique,isoform}`, `--method {weighted,pooled}`, an optional `--min-r2` gate,
  `-W` refit, and a GUI Protein tab with a per-protein refit curve. Load-bearing
  locked decisions (the design rationale that survives the shipping):
  - **Parsimony is a summarize-time decision, not integrate-time** (locked
    2026-06-10). A shared peptide's envelope blends both proteins' turnover and
    can't be attributed, so `integrate --unique` was **removed** — integrate
    extracts all peptides, `rollup --parsimony` attributes. `unique` drops
    multi-accession peptides; `isoform` (from `02_R_parsimony_reference.Rmd`)
    folds isoform-shared peptides onto the canonical entry unless an isoform
    carries its own unique peptide, rejecting cross-gene groups. A dataset-wide
    `_resolve_parsimony` pass; θ collapse weights by the **M5 inverse variance**.
  - **Biorep-aware collapse** (the `weighted` default): the per-timepoint
    weighted-average is *within* `(protein, labeling time, biological replicate)`
    — peptides of the same protein in the *same animal* are pseudoreplicates,
    different `characteristics[biological replicate]` are genuine replicate points
    → honest degrees of freedom. `pooled` is the pseudoreplication-naive
    comparator; point estimators were dropped (carried as a `peptide_median_k`
    column instead).
- **`linear simple` cross-sample Δk model (SHIPPED 2026-06-12; CHANGELOG).** A
  φ-space OLS model choice (mutually exclusive with the nonlinear ODE models):
  φ = `log(1−θ)` through-origin OLS `φ ~ 0 + day + day:condition` (statsmodels) →
  per-condition k (= −slope) + a pairwise **Δk** test + Benjamini-Hochberg; CIs
  are **analytic**, not bootstrap. Validated on `runs/lve_atr` (780 proteins, 410
  sig at BH p_adj<0.05) + GUI φ-space plotting. Two decisions that still bind:
  - **⚠️ Plateau truncation is required and linear-only** (`--phi-limit`, default
    −4 ≈ θ 0.98). Once φ saturates, later timepoints are floor-noise that drag the
    through-origin slope flat; the nonlinear models instead *need* the plateau
    (their asymptote fits it). Inseparable from `linear simple`.
  - **statsmodels is a dependency** (chosen over hand-rolled contrasts — wanted
    anyway for the deferred mixed/`limma`-style models). Reference notebook
    `data/notebook/03_R_linearmodel_reference.Rmd` (study-specific aspects do
    **not** carry over).
  - **Parked — >2-condition pairwise contrasts:** Δk needs exactly two qualifying
    conditions; `fit_linear_deltak` already emits per-condition k for any N and
    leaves Δk NaN when ≠2 (forward-compatible). Blocked on a ≥3-condition dataset.
- **o18 rewrite.** Same SSE-vs-IsoSpec → fit machinery as D₂O, but FS is
  computed with o18 isotope mass/proportion **and a different Spep model**: the
  labelling sites are only the oxygen-bearing residues, so the coefficient
  schema is a regression on `length, #D, #E, #N, #Q, #S` (possibly #T/#Y) — *not*
  a per-AA-over-20 dict. This needs a distinct `--label o18` dispatch (different
  coefficient format) and a frozen o18 coefficient table built via
  `data/notebook/90b_O18_LearnAALabelingSites_IsoSpec_AC16.ipynb` (the analog of
  87a for D₂O). Until then `riana fit --label o18` errors; o18 *integration* is
  unaffected.
- **M7 — PTM-aware envelope. v1 SHIPPED 2026-06-13 → 2026-06-17 (CHANGELOG).**
  Atom-vector `[C,H,O,N,S]` → `[C,H,O,N,S,P]`; a curated UNIMOD-id-keyed
  `mod_atoms` table (CAM as `UNIMOD:4`, retiring the `iaa` flag); variable mods
  threaded IO→integrate→fit (each `[UNIMOD:N]` form a distinct `concat` at its own
  m/z + envelope); proteoform rollup keys (`A2ASS6_pS34476`) driven by a
  `BIOLOGICAL_MODS={21}` (phospho gets its own key; N-term Ac / CAM / Met-Ox fold
  into the bare accession); `-X/-F` retired. Baseline `runs/lve_atr_m7/`. Two
  design notes that still bind the follow-on work:
  - **Two separate fixes, both required:** the integrate-side extraction target
    m/z *and* the fit-side IsoSpec `formula` must both carry the mod atoms (the
    pre-M7 drop sidestepped both). Mod hydrogens stay out of `num_labeling_sites`
    (mod D₂O enrichment unknown). The frozen 5-element M2 oracle is **not** touched
    — production asserts byte-identity against it on unmodified peptides.
  - **Site localization → proteoform key, no FASTA needed.** *mzTab* gives
    `start`/`end` + peptide-relative `pos-UNIMOD:id` → protein site `start+pos−1`
    (DONE). *DIA-NN parquet* gives `Protein.Sites` pre-formatted `[acc:res+pos]` +
    `PTM.Site.Confidence` (the **deferred** DIA phospho path — gate the key on
    confidence, filter to variable mods of interest). *Percolator* has no
    protein-coordinate site → no PTM support there.

  **Far-future — per-experiment chemical-vs-biological mod tagging.** Today the
  chemical/biological/encode policy is *global* id-sets in `constants`
  (`BIOLOGICAL_MODS`, `CHEMICAL_MODS`, …). Two cases break that: (a) the **same
  UniMod is chemical in one experiment, biological in another** — dimethyl
  (`UNIMOD:36`) is a chemical duplex label (Sadygov/Deberneh) but biological for
  histone methylation; (b) the **same UniMod means different things by site** —
  Acetyl (`UNIMOD:1`) is constitutive at the protein N-terminus but regulated K-ac
  on a side chain (handled today only by the crude pos≥1 site guard). The real fix
  is an **experiment-scoped policy**: read the modified-residue type from the SDRF
  UniMod and/or a CLI/GUI tag declaring, per run, which UniMod is chemical vs
  biological (vs a separate channel). Deferred — YAGNI until a user needs the
  conflicting interpretation; recorded so the global-set design isn't mistaken for
  the final word.

  **Roadmap — which mods come next, and the binding constraint.** The hard gate is
  **identifiability in a search over *un*enriched data**: no PTM-enrichment
  D₂O-labeling dataset exists yet, so we can only measure turnover of PTM forms
  detectable in the ordinary global-proteome runs. That, not envelope difficulty,
  is what sequences the list.
  - *Tier 0 (v1):* **phospho-STY**, **protein N-term Acetyl** — high value, reliably
    found unenriched (N-term Ac is near-universal/high-stoichiometry; abundant
    phosphosites do show up without enrichment, just fewer).
  - *Tier 1 — biological, own key, acceptable unenriched yield:* **Lysine acetylation
    (K-ac, `UNIMOD:1`) — SHIPPED 2026-06-21** (own `_acK###` key; N-term Ac folds via
    the `pos<1` guard). **K/R methylation** (mono/di/tri, `UNIMOD:34/36/37` =
    `[1,2,0,0,0,0]` / `[2,4,…]` / `[3,6,…]`) is **relegated to Tier 2** (2026-06-21,
    low demand) — still cheap (composition already in `mod_atoms`, no atom-vector
    work) but no longer near-term.
  - *Tier 1b — artifactual / chemical, fold-into-bare, envelope fidelity:* **Met-Ox
    (`UNIMOD:35` = `[0,0,1,0,0,0]`)** and **deamidation N/Q (`UNIMOD:7` =
    `[0,-1,1,-1,0,0]`)** — ubiquitous, trivially identified, not turnover units, but
    accounting their atoms **recovers the abundant modified peptidoforms we
    currently drop**, raising bare-protein peptide counts. **Chemical mods need a
    peptide-level merge, not just envelope accounting** — see the boxed design
    below; deamidation is the harder of the two (isobaric-overlap regime).

  **Chemical-mod handling — integrate-separate, fit-merge (Met-Ox SHIPPED
  2026-06-17, CHANGELOG; deamidation deferred).** The pattern future chemical mods
  follow: a purely chemical mod happens *post-synthesis*, so it does **not** reset
  the D₂O clock — the oxidized and unoxidized forms share the same FS-vs-time
  signature and must fold onto **one** curve, not two underpowered ones. Mechanism:
  **integrate the forms separately** (each at its own clean m/z — A2 already makes
  each `[UNIMOD:N]` a distinct `concat`), then **merge at the fit level** via a
  chemical-mod-stripped `_fit_key` (`CHEMICAL_MODS={35}`) that keeps biological
  mods distinct. The fit key is the layer *between* the integrate `concat` and the
  Stage-B proteoform key — and it sidesteps the "which m/z to integrate?" problem
  (never integrate a blended m/z). LVE_ATR motivation: of 545 oxidized Met
  peptidoforms, 451 also appear non-Ox → previously double-fit; merging gives the
  consolidation/power win (raw "fittable-series count" *drops* under dedup, so it
  is the wrong metric).
  - *Deamidation is the hard case — conditionally tractable, mass-gated.* +0.98401 Da
    (N→D / Q→E) sits only **0.0193 Da below the C13 M+1** (`1.00335`), i.e. a
    separation of **≈ 19340 / M_neutral ppm** (~19 ppm at 1000 Da → ~6 ppm at
    3000 Da; charge-independent in ppm, but higher z lands the peak at lower m/z
    where the Orbitrap resolves better). Because deamidation is **non-stoichiometric
    and time-correlated** (old/unlabeled protein is most deamidated), an unresolved
    deam-M0 bleeds into the non-deam M+1 *in lockstep with the labeling state* — a
    smooth, systematic bias that can warp `k` and pass the R²>0.95 curation silently.
    Two gates decide separability: (1) MS1 resolution actually centroided them apart
    (R ≈ M/0.0193 — ~52k at 1000 Da, ~104k at 2000 Da, lost by ~3000 Da or at fast-DIA
    15–30k MS1), and (2) a tight extraction window (±3–5 ppm, not the default ±10).
    **Where both hold** (small/mid peptides, high-res MS1, tight window) → it reduces
    to the Met-Ox case (integrate-separate + fit-merge). **Where they don't** (large
    peptides, low-res MS1, wide window) → merged centroid, irrecoverable → **flag /
    exclude the peptide (both forms)**. The *rigorous* alternative for the
    unresolvable regime — **jointly modelling both forms' FS/D₂O IsoSpec envelopes
    AND the (unknown, time-varying) deamidation proportion** — is a **standalone
    side project**, deferred (user 2026-06-13). So the near-term policy is a
    per-peptide **mass/resolution gate** (`19340/M_neutral` vs the window ppm), not
    blanket merge or blanket exclude. The bio-vs-chemical distinction is **moot** — both are
    the identical +0.984 mass shift, unactionable; FS-merge is valid for either, so
    the only question is separability. When deamidation lands, first measure the MS1
    resolution + deam-peptide mass distribution (e.g. on LVE) to size the separable
    fraction. This is a **fit-stage** concern; none of it is in v1 (Ox/deam are
    dropped), so it does not affect Stage B.
  - *Tier 2 — highest biological value, BLOCKED on data:* **Ubiquitin/ISG15
    GG-remnant (`UNIMOD:121`, GlyGly = `[4,6,2,2,0,0]`)** — it literally *is* the
    degradation tag, so it is the most turnover-relevant PTM imaginable, **but**
    diGly-remnant peptides are essentially undetectable without anti-K-ε-GG
    enrichment, and no enriched D₂O dataset exists. Deferred until such data lands
    (same blocker for the acyl class — succinyl / malonyl / crotonyl — that needs
    enrichment). Worth naming now precisely because the science case is strongest:
    the moment an enriched D₂O-diGly dataset appears, GG jumps the queue.

- **TMT — chemical isobaric label as a fit-merge mod (NEAR-TERM, designed
  2026-06-21; not built; SEQUENCED AFTER adaptive N_ISO — see Track B). Memory
  `track_c_tmt_chemical_mod`.** Support **TMT10plex
  (`UNIMOD:737`, identical chemistry to TMT6plex)** and **TMTpro 16/18-plex
  (`UNIMOD:2016`)** as **chemical fit-merge** mods (the Met-Ox bucket): searched
  *variable* to catch incomplete labeling, a peptide shows up as 1-tag vs 2-tag
  peptidoforms (different MS1 masses, integrated separately), and since TMT is
  post-harvest it shares one D₂O curve → strip in `_fit_key` (add to
  `CHEMICAL_MODS`). NOT a biological key; NOT a separate channel (dimethyl/SILAC are
  the separate-sample TODO below).
  - **⚠️ The trap — TMT carries fixed heavy isotopes.** TMT6plex is
    `C8 ¹³C4 H20 N1 ¹⁵N1 O2` (Δ 229.1629), TMTpro `C8 ¹³C7 H25 N1 ¹⁵N2 O3` (Δ
    304.2071), so the light-composition mass ≠ the true mass. "Total atoms in
    `mod_atoms`" gives a ~5 / ~9 Da-light precursor (extraction misses the peak);
    "light atoms + a separate mass override" is **also broken** — the FS solver
    matches **by channel index** (so a global m/z offset is harmless), *but*
    `get_envelope` bins the IsoSpec distribution at `round(pep_mass) ± 0.5`, so a
    light distribution sitting ~5 Da below the correct `pep_mass` anchor falls
    outside every bin → empty envelope → NaN on every TMT peptide.
  - **The correct fix — pinned single-isotope pseudo-elements.**
    `get_peptide_distribution` already appends a pseudo-element for deuterium the
    same way; append the fixed heavies as `count, masses=[¹³C/¹⁵N], probs=[1.0]`. A
    100%-abundance isotope adds **mass but zero broadening** — physically exact.
    Data model: `mod_atoms[id]` = light broadening atoms; new `mod_fixed_isotopes`
    table. This **pinned-isotope envelope support is the shared prerequisite** and
    is needed regardless.
  - **Precursor mass — resolved 2026-06-21: comes from `iso0`, not a separate mass
    path.** With **adaptive N_ISO sequenced first** (Track B — IsoSpec runs at
    integrate), the envelope's `iso0` *is* the precursor m0, so TMT needs **no**
    analytic `unimod_mass + pinned` extension at all — the bridge (Option A) is
    dropped. (The A/B was: A = ship TMT now with the analytic mass; B = do
    IsoSpec-at-integrate first and ride `iso0`. **B chosen** because adaptive N_ISO
    + fit-time isotopomer subsetting is higher-value than TMT — user 2026-06-21.)
    Until adaptive N_ISO lands, TMT is simply not built. ~tests + the `constants`
    tables (`mod_atoms` light entries + `mod_fixed_isotopes` + add 737/2016 to
    `STARTER_VARIABLE_UNIMODS` ∪ `CHEMICAL_MODS`); frozen M2 oracle stays
    5-element/light (assert unmodified byte-identity). **GUI fold-points are
    pre-wired for `tmt`** (built 2026-06-21), so the merged points colour the moment
    this lands.
- **Dimethyl duplex & SILAC — "multiplex/channel labels, fit separately" (TODO).**
  *Distinct from the fit-merge bucket above* (user correction 2026-06-21): a chemical
  duplex/multiplex label (reductive **dimethylation** — Sadygov/Deberneh; **SILAC**)
  marks a *different sample/channel*, so its light/heavy forms must be fit as
  **separate curves**, not merged. That makes the mod a **sample/condition axis**
  (conceptually an SDRF channel dimension) — a bigger architectural feature than the
  `_fit_key` merge, hence deferred. SILAC is biological (metabolic); dimethyl is
  chemical-but-still-separate. Note the per-experiment ambiguity (dimethyl `UNIMOD:36`
  is *chemical* in a duplex labeling study but *biological* for histone methylation)
  — resolved by the far-future per-experiment chemical-vs-biological UniMod tagging
  (Track C M7 box).

- **Animal within-protein-θ benchmark — BUILT (2026-06-10).** The 4th dataset
  (`data/timeseries_lve`, mouse in-vivo D₂O, mzTab + SDRF) has no fractional-pool
  ground truth, so it is scored by *within-protein θ/FS variance per timepoint*
  (Hammond 2022): a better integrator minimizes the spread among peptides of the
  same protein. Shipped as `build_lve_bench_set.py` (frozen, integrator-
  independent set: proteotypic peptides of ≥3-peptide proteins at q<0.01,
  abundance tertiles) + `bench_within_protein_theta.py` (per-timepoint θ via the
  production `solve_fs_d2o`; overall / per-stratum / curated-vs-uncurated;
  `compare_methods` so a Track B variant slots in on the identical set). Baseline
  at the adopted defaults (apex@0.15 + best-q anchor + 10 ppm): within-protein-θ
  robust-SD median **0.10** (curated 0.044), correct abundance ordering. Building
  it surfaced and fixed **two integration bugs**: (1) the quantms **filename-prefix
  scan-scramble** (Track A gotcha above) and (2) **mass_tol 50→10** (centroid
  search-tolerance, not profile width — LVE R²med 0.69→0.87). Now the orthogonal
  regression gate for Track B.
- **Curation-gate revisit.** Treat both knobs as tunable — the R²>0.95 threshold
  *and* the "observed at every proportion" coverage rule (it should tolerate a
  known-bad fraction rather than discard peptides wholesale). Report the R²
  *distribution* shift before/after a fidelity change, not just curated counts.
- **DIA test dataset (TODO, enables M6b).** A small quantms-diann run lands under
  `data/timeseries_dia` (DIA-NN `report.parquet` ≥ 2.2.0 + the DIA-variant SDRF),
  the fixture the `io/diann.py` adapter and DDA/DIA auto-detect are built and
  tested against.

#### Track E — GUI/UX & responsiveness (surfaces the engines above)

- **Sortable tables — DONE 2026-06-21.** All three tabs (Integrate/Model/Protein)
  set `setSortingEnabled(True)`; sorting is implemented as `DataFrameTableModel.sort()`
  reordering the backing frame in **pandas** (not a `QSortFilterProxyModel`). Two
  reasons the model-internal sort beat the proxy here: numeric columns sort by
  value, not by the `:.4g` display string (lexicographic `"100" < "9"`), and the
  row-click handlers keep using `dataframe.iloc[row]` directly with no
  `mapToSource` translation — the proxy would have silently mis-mapped every
  selection. The vectorised pandas sort is also the **large-table responsiveness**
  win: a proxy comparing cells in Python doesn't scale to large result frames.
- **Save/export graphs — DONE 2026-06-21 (PNG).** A "Save graph…" button under the
  chromatogram + fitted-curve views exports the live `PlotItem` via pyqtgraph's
  `ImageExporter` (`riana/gui/export.py`); replaces the removed `--plotcurves`.
  Disabled in the placeholder state. **Raster only:** pyqtgraph's `SVGExporter`
  throws on plots with scatter symbols (the observed-point / fold series these
  views always draw), so SVG was dropped rather than ship a button that crashes on
  the common case. A faithful **matplotlib static** export remains a future option
  if a vector figure is needed.
- **Δmass-over-time QC display** — the near-term half of dual-mode FS: surface
  the per-isotopomer accurate-mass shift over the init envelope across the time
  series, a high-level QC for advanced users.
- **Make the displayed chromatogram trace faithful to smoothing.** *(The form
  control half is DONE 2026-06-20 — smoothing window + poly-order are in the
  Integrate tab's Advanced group.)* What remains: the chromatogram re-extracts the
  *raw* XIC from the mzML at row-click, so the plotted trace does not reflect the
  S-G smoothing integration actually applied. Fix by either re-applying the run's
  `IntegrationConfig.smoothing` at click or — better — reading the persisted
  smoothed trace from the `--save-traces` sidecar (which also removes the ~10 s
  re-read).
- **Fit progress + real parallelism** (engine is Track C / `core/pipeline.py`,
  surfaced here). Today the GUI submits one opaque `run_fit` to a single worker
  with an indeterminate "busy" bar, and the inner peptide parallelism is a
  GIL-bound `ThreadPoolExecutor`. Make `fit_run` chunk-aware (peptide chunks),
  dispatch chunks across the **ProcessPool**, emit a progress callback the bar
  maps to. The CLI gets the same speedup. **Reproducibility caveat:** chunking
  across processes requires a deterministic per-peptide bootstrap seed (e.g. from
  `concat`), or output becomes scheduling-dependent — add a parallel == serial
  test (rec below).
- **Modernize look** — low priority, ready-made only (a QSS theme or
  `qt-material` / `qtmodern`).
- **`-W/--workers` in the GUI fit + rollup forms — DONE 2026-06-13.** Integrate
  already had it (the "Workers (files)" spin + semaphore dispatch over the shared
  pool). Added a **Workers** spin to the Model (fit) and Protein (rollup) tabs.
  **Nesting guard (verified the hard way — a `ProcessPoolExecutor` inside a pool
  worker raises `BrokenProcessPool` here):** when `workers > 1` the tab dispatches
  the run on a **main-process thread** (`run_in_executor(None, …)`) so `fit_run` /
  `rollup_proteins` create their pool from the *main* process, not nested inside a
  shared-pool worker; `workers == 1` keeps the shared pool (one CPU task off the
  main process). Per-peptide/-protein deterministic seeds keep the result
  worker-count-independent.
- **Expose the hidden integrate knobs in CLI + GUI as clearly-marked *advanced*
  options — DONE 2026-06-20 (CHANGELOG).** Every `IntegrationConfig` dial is now
  reachable on both surfaces: 6 were frozen at default with no flag/widget
  anywhere (`prominence_k`, `width_rel_height`, `apex_n_consensus`,
  `scan_rt_tol_min`, `smoothing_polyorder`, plus the GUI-only `ppm_alert` → added
  `--ppm-alert`). CLI groups them in "Advanced integration" + "Match-between-runs
  (MBR)" `--help` panels; the GUI grew a collapsed "Advanced…" group (the form is
  now scrollable) with the full dial set + a checkable MBR sub-group, and
  `build_config()` passes every field. The smoothing-in-GUI item above is
  subsumed (smoothing window + poly-order are in the Advanced group); the
  *faithful-displayed-trace* half (re-apply S-G at row-click / read a
  `--save-traces` sidecar) is still open.
- **GUI-framework feasibility report — is Qt a ceiling? (spiked 2026-06-12,
  long-term.)** If pyqtgraph/Qt limits graphing or interactivity as the protein /
  Δk / φ-space views grow, write a feasibility note comparing alternatives
  (Electron + JS charting, Tauri, a web/server split, Dash/Streamlit for the
  analysis surfaces) against the cost of leaving the native PySide6 app. Decision
  doc only — no migration implied; revisit only if a concrete Qt limit blocks a
  needed view.

#### Persist the integration signal (optional sidecar)

The GUI re-extracts a peptide's XIC from the mzML on every row click (~10 s). Add
an optional `--save-traces` → `<stem>_riana_traces.parquet` (keyed `(concat,
isotopomer) → (rt[], intensity[])`), never in the `.txt`; default off for batch
(size), on for GUI/diagnostic workflows. This is purely a GUI-speed/diagnostics
feature — the protein refit consumes per-(peptide, t) θ (the M5 intermediate),
not the raw XIC, so the two are decoupled.

*2026-06-28 spiked:* GUI should display results (from fit/rollup and optionally
integrate when manifest is present)
- Remove legacy Percolator path from GUI (keep in CLI for testing/dev and bench 
scripts, GUI users won't need that as most features require SDRF path) (2026-06-28 shipped)
- GUI integrate view should display sample/fraction information when SDRF is present


#### Suggested sequence

1. **Pre-1.0.0 chores** — repo hygiene + cut a real `1.0.0` tag (below). A clean
   line in the sand before piling on features.
2. **M6a** — the spine; unblocks Track D and the protein layer; pulls in
   `core/pipeline.py` and integrate concurrency.
3. **Track D** — stand up the orthogonal animal gate while M6a is fresh.
4. **M5** — small, but gates the protein layer.
5. **Track C protein rollup + Track E polish** — the protein tab needs M5 + the
   M6a manifest; GUI polish anytime.
6. **Track B** — the fidelity research, now with two gates ready.
7. **o18 / M7 / M6b / dual-mode FS** — slot in by demand and bandwidth.

#### Cross-cutting / pre-1.0.0 chores

- **Cut a real `1.0.0` tag + Zenodo *code* DOI (RELEASE ACTION — pending).**
  `riana/__init__.py` is **already `1.0.0`** (the version bump is done; this header
  prose elsewhere still says `1.0.0.dev1` — stale). Last tag is `v0.9.0`. What
  remains is the **release decision**: whether to merge `m3-rewrite` → `master`
  first, then `git tag v1.0.0` + push + mint the Zenodo code DOI. Left to the
  maintainer (branch/timing + an external service).
- **Repo hygiene — largely already clean (audited 2026-06-21).** `data/` is
  **fully gitignored** (27 GB local-only — a personal-disk concern, not a git one);
  there are **no stray *tracked* root outputs**; `docs/` has **no committed Quarto
  HTML**; `riana_website/` is legitimate Quarto **source** at normal `0755` perms
  (the old "mode 0700" note is stale — keep it). Only gitignored junk remains
  locally (`.DS_Store`, `.coverage`, `logfile.log`). No action needed beyond the
  optional local `data/` disk cleanup, which is the maintainer's.
- **User-facing docs refresh — post-M3, still stale (large TODO).** The narrative
  docs were not redone after the M3/M4 rewrite, so they lag the current CLI/GUI,
  the SDRF/manifest model, the rollup, and now `linear simple`. Until that
  happens, the **docstrings are the source of truth** and should carry the
  non-obvious contracts explicitly — e.g. `linear simple` CIs are analytic, not
  bootstrap (now noted in `core/linear_model.py`). Track a full docs pass as its
  own chunk of work, not a side-effect of feature PRs.
- **Snakemake — retire the bundled `workflow/Snakefile`** (supersedes §4.4). With
  Percolator demoted and quantms / DIA-NN handling search + ID end-to-end
  upstream, Riana's own pipeline collapses to a linear `integrate → fit →
  protein` chain whose glue is the manifest (each step reads the prior's manifest
  rows). Snakemake's DAG/partial-rerun value largely evaporates; ship the three
  CLI subcommands (optionally a thin `riana run` convenience wrapper) and stay
  orchestration-agnostic so Riana composes into whatever workflow already runs
  quantms.
- **mypy --strict** rollout (§4.3); `py.typed` already ships.
- **Test-suite runtime (~931s; profiled 2026-06-13).** ~845s lives in ~7 real-mzML
  tests, all single-fraction `integrate_run` calls (`test_ac16_time0_matches_
  committed_baseline` alone is **415s**; six `sample1`/cli integrate tests are
  51–114s). **`-W` does not help these** — `max_workers` is an `integrate_project`
  cross-*run* lever, but these integrate one fraction each and `integrate_run`'s
  per-PSM loop has no internal pool, so `n_parallel = min(workers, 1) = 1`. The
  effective levers (deferred — recorded, not done): (1) a **`slow` pytest marker**
  on the heavy real-data regression gates so local dev runs `pytest -m "not slow"`
  (the 415s ac16 test is already `skipif`-guarded on gitignored heavy inputs, so it
  only runs where those exist); (2) **`pytest-xdist` (`-n`)** for GIL-free
  cross-test process parallelism (~2× — the 415s ac16 test becomes the wall-clock
  pole). Single-fraction `integrate_run` is intentionally left serial (cross-run
  `-W` covers the multi-file production case).

#### Additional recommendations (added this round)

1. **Benchmark CI smoke tier.** The `bench_*` suite runs manually; add a
   subsampled fast tier to CI so a science regression (e.g. a fidelity change that
   worsens within-protein θ) fails the build.
2. **Reproducibility as an explicit invariant.** Add a parallel-fit == serial-fit
   test (forces the deterministic per-peptide seeding); protects the
   header-authoritative frozen-output story.
3. **Schema-version the manifest/header** so the breaking output-granularity
   change (timepoint-named → mzml-named + manifest) and future format changes are
   detectable downstream.
4. **SDRF as a partial config source** (scoped): read mods from the SDRF; keep
   integration mass tolerance a separate, wider Riana parameter (not the search
   tolerance).
5. **Sequence the breaking-change migration with the Snakefile retirement** so
   there is never a broken intermediate workflow state.


## 4. Cross-cutting recommendations

Mostly shipped in 1.0.0 (see `CHANGELOG.md`): provenance-header reproducibility,
the frozen `IntegrationConfig`/`FitConfig` single source of truth, `from __future__
import annotations` throughout, and the 2026-06-21 repo-hygiene audit. **Open:**
retire the bundled `workflow/Snakefile` — with quantms/DIA-NN owning search + ID,
Riana's pipeline is a linear `integrate → fit → rollup` chain glued by the manifest,
so ship the CLI subcommands and stay orchestration-agnostic.

## 5. What this plan deliberately does not include

- No "AI/ML peak detection," no cloud, no plugin architecture, no
  database backend, no web UI. Earlier evaluations suggested all of
  these; none are right for a single-author academic tool with this user
  base.
- No new kinetic models beyond the existing three. Focus on getting the
  existing three to produce reliable numbers with honest uncertainty.
- No multi-omics integration, no Kubernetes, no "plugin ecosystem with
  >10 community plugins." Same reason.

## 6. Known limitations (honest list)

- Memory: fixed in 1.0 via streaming/indexed reads (`io/mzml.py`,
  `IndexedMzML` — one fraction in memory). The 0.9.0 list-of-arrays
  loader (`spectra.py`) is removed.
- Peak fidelity: see §2c. The calibration dataset quantified the impact;
  the M3 spike adopted an apex-centred narrow window as the default.
- AA / SILAC fitting: dropped from `riana fit` in M4 (the buggy
  `--label 4` `a_0` path is gone). The SILAC dual-channel extraction knobs
  (`-F/--forced_mods`, `-X/--ignored_mods`) were retired in M7 Stage A3 — mods
  are now handled uniformly via UniMod; integrate extracts one channel per
  peptidoform at its own m/z.
- o18 fitting: recognized (`--label o18`) but errors pending a post-M4
  rewrite (length + selected-residue coefficients, not a per-AA dict).
  o18 *integration* is unaffected.
- GUI: the broken Tkinter `riana_ui/` is deleted (M4 Phase 1); the PySide6
  replacement shipped in M4 Phase 2 (`riana gui` — Integrate + Model tabs).
- The bundled `workflow/Snakefile` is an example, not a tested pipeline.
  Treat as a starting point.
