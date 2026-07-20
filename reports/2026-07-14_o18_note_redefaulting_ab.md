# ¹⁸O technical note — re-defaulting A/B: how the 1.2.0 fit/rollup changes move the note's numbers

**Date:** 2026-07-14 · **Status: analysis complete, note update pending** · **Note repo:** `../riana-o18-note`

## Why

The ¹⁸O technical note (`../riana-o18-note`) was built at **RIANA 1.1.0**. Three (really
**four**) results-affecting changes have since shipped as `1.2.0` **defaults**, so the note
currently documents superseded behaviour — and in two places cites flags that no longer exist:

1. **FS rail-drop** (`fit`, default ON, rails `1.05 / −0.05`) — drops per-timepoint FS points
   the solver pinned to its bound, before depth/points counting. `dea65f1`.
2. **k_cv rescue gate** (`rollup`, `--k-cv 0.2 --rescue-r2 0.6`) — a scale-free relative-uncertainty
   admit that **replaced** the old `--alt-k`/`--alt-se` slow-turnover admit. `8312374`.
3. **linear-simple WLS** (`rollup`, `--linear-weights wls`, default) — fitted-value-weighted LS
   on φ, fixing OLS's ~28 % Type-I inflation. `4ae5e50`, report `2026-07-13_linear_model_wls.md`.
4. **t = 0 exclusion** in the linear-φ fit — bundled into (3) but *unconditional* (see Bugs below).

**Data / method.** Each of the note's six runs (`boomi_ipsc_{o18,d2o}`, `juber_ac16_{o18,d2o}`,
`timeseries_lauren5_7_ipsc_mesoderm_o18`, `timeseries_lauren9`) was **re-fit twice from its
existing integrate outputs** — `--no-fs-rail-drop` (legacy) and default (`--fs-rail-drop`) — at the
note's uniform **depth 6** + label-aware Spep gate. No mzML re-integration. The cheap rollup
factorial (gate × linear-model) then runs on top of each fit arm, giving the cumulative ladder and
every single-factor isolate. The **L0 baseline uses the note's *actual* rollup code** (worktree at
`fa72d50`) because the current code cannot reproduce the retired `--alt-k` gate at any flag setting.

- Driver: `runs/note_ab_ladder.sh` (12 fits + 48 rollups) · `runs/note_ab_legacy.sh` (legacy rung)
- Analysis: `runs/note_ab_analyze.py` → `reports/_note_ab_tables.md` (full tables)
- Metric definitions mirror the **note's notebooks** (spep≥6 o18 / ≥8 D₂O, robust geom-CV =
  1.4826·MAD(ln k), Spearman on shared peptidoforms), so numbers are directly comparable to the note.

**Reproduction is exact.** The re-fit `railoff` arm regenerates the note's published fit counts to
the peptide (ac16 ¹⁸O: 9,401 converged), and the legacy rollup on it regenerates §6's headline
**to the protein: 869 tested / 519 sig / 253 faster / 266 slower** — identical to the note. Both
stages of the baseline are verified reproductions, so the deltas below are trustworthy.

## Headline: §6 (biological application) — the one place a *conclusion* moves

§6's Δk analysis (mesoderm − iPSC, `timeseries_lauren5_7`) is where the change bites. Cumulative
ladder, each rung adding one shipped default on top of the last:

| rung | tested | **sig (q<0.05)** | faster | slower |
|---|---|---|---|---|
| **L0** note today (legacy: alt-k gate, OLS, t0 kept) | 869 | **519** (59.7 %) | 253 | 266 |
| L1  + FS rail-drop | 750 | 419 (55.9 %) | 208 | 211 |
| L2  + k_cv gate (replaces alt-k) | 652 | 369 (56.6 %) | 196 | 173 |
| **L3  + WLS = current RIANA default** | 652 | **225** (34.5 %) | **84** | **141** |

**Two conclusion-level shifts:**

1. **The significant count more than halves, 519 → 225**, and **WLS is the dominant mover**
   (369 → 225 at fixed L2 curation — a 39 % cut). This is exactly the OLS Type-I inflation the WLS
   report predicted (OLS ~28 % false-positive on RIANA's own data); the note's 519 was inflated by
   OLS's anti-conservative p-values, not by more real biology.
2. **The "near-symmetric split" claim breaks.** §6.2 states *"253 faster and 266 slower … the split
   is near-symmetric and there is no global scale shift."* At the new default it is **84 faster /
   141 slower** — slower-in-mesoderm now dominates ~63 %. WLS preferentially culls the *faster*
   program, which §6.4 calls the sharp, high-confidence translation-apparatus signal. That framing
   needs rewriting.

**The protein-by-protein thesis survives; the pathway-level story does not.** Of L0's 30
most-significant proteins, **22 still significant at L3**; among *all* surviving significant proteins
only **5 sign-flip**. Named-mover survival (measured directly, OLD backup vs NEW canonical):

| direction | named movers (§6.3) | survive at WLS default |
|---|---|---|
| **slower** (pluripotent metabolic) | PGK1, RRM2, PHGDH, GLDC, EPRS1, UBA1, UGP2 | **6/7** — PGK1/PHGDH/EPRS1/UBA1/UGP2 strong (q≤1e-10); RRM2 weaker (q 6e-3); **GLDC drops out** of the admitted set |
| **faster** (translation/other) | EEF2, HNRNPH1, IDI1, GART | HNRNPH1/IDI1 hold; **EEF2 collapses** (Δk 0.030→0.011); **GART flips sign** (+0.029→−0.005) and dies |

**Fig 8 pathway enrichment DISSOLVES at the honest gate (verified live, 2026-07-14).** Re-running the
note's exact g:Profiler ORA (custom background = the tested proteins, BH-FDR<0.05, term size 5–500,
generic roots dropped) on the new sets:

- **faster (84):** **0** enriched terms vs the measured-proteome background. (Whole-genome background
  returns 217 terms topped by "cytoplasmic translation" — but that is exactly the abundance bias the
  note's custom background exists to reject; the translation genes are *present* but not
  *over-represented* relative to the measured proteome.)
- **slower (141):** only 4 terms, all generic roots (biological_process, cellular process, KEGG root,
  Metabolic pathways) — all dropped by the note's own GENERIC filter → **0** meaningful terms.

So **§6.4's entire two-program pathway result was riding on OLS's Type-I inflation**: the 253-protein
faster set cleared FDR, the honest 84-protein set does not. This inverts the note's asymmetry claim
(it called faster-translation the sharp signal) — the durable signal is the *slower* metabolic
individual movers, and the *faster* translation program is the artifact. **§6.4 / Fig 8 cannot stand
as written.** Options for handling this are a pending user decision (drop / report-null / rank-based
GSEA replacement).

## §5.2 (orthogonal-label validation) — yields rise, conclusions stable

The cross-label story is **robust**. Adopting the shipped default (railon + k_cv) as primary, per
the user's decision:

| pair | curation | yield ¹⁸O / D₂O | Spearman(k) pep | within-prot geom-CV ¹⁸O / D₂O |
|---|---|---|---|---|
| **AICS52 (boomi)** | L0 note (railoff, R²≥0.8) | 46.7 / 36.9 % | 0.679 | 0.142 / 0.140 |
| | **new default (railon, k_cv)** | **61.0 / 53.9 %** | **0.662** | 0.168 / 0.168 |
| **SCVI480 (lauren)** | L0 note | 10.7 / 15.3 % | 0.668 | 0.173 / 0.139 |
| | **new default** | **25.3 / 30.6 %** | **0.646** | 0.201 / 0.187 |
| **AC16 (juber)** | L0 note | 9.0 / 8.5 % | 0.673 | 0.159 / 0.136 |
| | **new default** | **21.7 / 20.3 %** | **0.649** | 0.172 / 0.179 |

- **Yield rises everywhere** (rail-drop cleans the fitted population — median R² of *all* fitted
  peptides jumps, e.g. boomi ¹⁸O 0.778→0.837, ac16 ¹⁸O 0.071→0.314 — and k_cv rescues flat-but-precise
  curves). The note's §5.2 "candidate default" CI-gate is now the shipped default.
- **Cross-label Spearman is essentially unchanged** (drops ≤ 0.02 at every pair) — the note's
  central validation claim (¹⁸O tracks D₂O on ranking) is untouched.
- **within-protein geom-CV rises modestly** (~+0.02–0.03) — the honest cost of admitting more,
  lower-R² peptides; still well within the note's "≈ 0.14" regime for the clean lines and far below
  any threshold that would change the "identical within-protein consistency" reading.
- **Absolute scale offset is stable**: median log₂(k_D₂O/k_¹⁸O) moves < 0.01 (boomi −0.464 → −0.458).
  The ~1.4× ¹⁸O/D₂O offset — a §5.2 headline — is unchanged, as expected (a multiplicative scale is
  invariant to curation).

**The old §5.2 over-rescue worry does not reproduce.** §5.2 reported a CI-rescue (`relunc<0.25`,
*no* R²-floor) that "over-rescues, geomCV → 0.64". The shipped gate's `--rescue-r2 0.6` floor is
exactly what tames it: at `k_cv<0.2 ∧ R²≥0.6` the geom-CV lands at 0.17–0.20, not 0.64.

## §4 (methods) — text is now factually wrong in two places

- **§4.5 cites `--alt-k` / `--alt-se`** ([`04_fs_kinetic_fitting.md:93`]) — **removed** in `8312374`.
  Replace with the k_cv gate (`--k-cv 0.2 --rescue-r2 0.6`).
- **§4.5 says the rollup fits "by OLS"** ([`04_fs_kinetic_fitting.md:87`]) — now **WLS** by default.
- **§4.3 depth/Spep curation** is unchanged and still correct; add the rail-drop as a fit default.

## Per-regime notes

- **¹⁸O rail-drop is *not* the strict no-op the multi-point report found.** That report's "¹⁸O no-op"
  is `1.05/−0.05` **vs the clamp** (both rail-ON). Here the comparison is rail-**ON vs OFF**, which is
  large for ¹⁸O: boomi ¹⁸O admitted 5704→6131 (r2 gate), ac16 ¹⁸O 468→594. The note's runs never had
  rail-drop at all, so they get the full ON-vs-OFF gain.
- **Rail-drop thins the depth-6 *fitted* population hard** (ac16 ¹⁸O 9,401 fitted → 3,009 railon)
  because railed points push curves below the depth-6 floor — yet **curated** yield *rises*, because
  what's dropped were failed solves. This is the intended mechanism (multi-point report), reproduced.
- **boomi/lauren D₂O** move the same direction as ¹⁸O (yield up, ρ stable) — the note's D₂O reference
  arms are consistent under the new defaults.

## Bugs found — 1 & 2 FIXED (2026-07-17), 3 cosmetic

1. **`--linear-weights ols` did not restore the pre-2026-07 fit — FIXED.** The t=0 drop at
   [`core/linear_model.py`] (`keep = t > 0`) was **unconditional**, not gated on the weights scheme —
   contradicting the module docstring, the CLI `--linear-weights` help, and the WLS report, which all
   say `ols` "restores the old unweighted fit for audit." Gated it on `weights != "ols"`, so `ols` now
   keeps t=0 (the genuine pre-2026-07 fit); plateau truncation stays common to both. Regression test
   `test_ols_keeps_t0_while_wls_drops_it`. *Note:* this changes what current-code `--linear-weights ols`
   produces (t=0 now kept), so the ladder's **OLS** rungs above (pre-fix, t=0 dropped) would shift by
   ≤ 6 proteins (349→355 / 230→235) if reproduced today; the **WLS default** rungs and every note number
   are unaffected (WLS still drops t=0).
2. **`rollup -o <fork>` / `fit -o <fork>` wrote a doubled `integrate` path — FIXED.** Fork manifests got
   `integrate` rows like `<absdir>/runs/<run>/runs/<run>/..._riana.txt` (non-existent): the seed resolved
   cwd-relative upstream paths against `proj_dir` (the manifest's folder) instead of the cwd the reader
   uses. Anchored the seed (and the fork-conflict check) at `Path.cwd()`. Harmless to rollup (reads `fit`
   rows) but would break any fork consumer that reads `integrate` rows (e.g. GUI chromatogram).
   Regression test `test_fork_seeds_cwd_relative_upstream_paths` (verified it fails on the old code).
3. **Report inconsistency (cosmetic).** `2026-07-05_multipoint_rail_drop.md` gives `juber_ac16_o18`
   rail-ON admit as both "+79 (+14.0 %)" (→ 645) in the A/B table and 719 in the ¹⁸O threshold table;
   my ladder gives 719 (matching the threshold table).

## Incident note

While validating the legacy baseline I ran the `fa72d50` rollup pointed at the **canonical**
`runs/timeseries_lauren5_7_ipsc_mesoderm_o18/riana_manifest.tsv`. The legacy code **ignores `-o`**
(the fork behaviour postdates it) and **overwrote** that run's `riana_rollup_{proteins,fractions}.txt`
in place (files are gitignored, so not git-recoverable). **No data was lost:** the legacy rollup is
deterministic and reproduced the note's published 869/519/253/266 exactly, so the overwritten file
is purpose-identical to what was there. All subsequent legacy runs write to isolated
`rollup_legacy/` dirs (`runs/note_ab_legacy.sh`) and never touch a canonical manifest.

## Decision taken (per user, 2026-07-14)

- **Primary curation for the note → RIANA's shipped default gate** (R²≥0.8 OR (R²≥0.6 ∧ k_cv<0.2)),
  with the strict R²≥0.8 numbers retained as a conservative sensitivity table.
- **Full clean re-fit, 6 runs × 2 arms** (done) — both arms code-matched, so the rail-drop delta is
  unconfounded.

## Note-repo update — DONE (2026-07-16)

**Fits promoted.** The railon (default) fits + `kcv_wls` rollups were copied into the note's
canonical `runs/<run>/` dirs (1.1.0 outputs backed up to `runs/_note_1.1.0_backup/`, 24 files).
Notebooks read those paths unchanged.

**Six figures regenerated** (fig1/2/3/6/M untouched — no fit/rollup dependency):

| fig | change |
|---|---|
| fig4 | curation → default gate; curated 6,131 → **6,671**, median k 0.0472/h, t½ 15 h |
| fig5 | default gate both levels; AICS52 pep ρ **0.662** (n=3950) / prot **0.689**; SCVI480 pep **0.646** / prot **0.647** (OLS-era prot 0.55 → WLS 0.65) |
| fig7 | WLS Δk: **652/225/84/141** (was 869/519/253/266); cross-state ρ **0.59** (was 0.41); volcano now slower-dominant |
| **fig8** | **rewritten to a negative result** — ORA 0 terms + GSEA 0/1129 (BH<0.25) vs the tested-proteome universe; figure shows the whole-genome bias as the teaching point |
| figS1 | promoted fits; t0 noise floor 0.11/0.13 (SCVI480) vs 0.07 (AICS52) |
| figS2 | D₂O degron GSEA **strengthens 5 → 8 motifs** (SPOP×2, FBXW7×2, β-TrCP, MIDN); ¹⁸O 1 → 0 (still underpowered) |

**Text rewritten** (all reflect 1.2.0 defaults): §4.3 (default gate + rail-drop), §4.5
(alt-k → k_cv, OLS → WLS), §5.2 table + robustness paragraph + cross-line paragraph + [^ac16] /
[^yieldgap] / [^degron] footnotes, §5.3, §6.2 (counts + drop "near-symmetric" + ρ 0.41→0.59),
§6.3 (slower-robust reorder; EEF2 collapse / GART flip noted), §6.4 (pathway null), §6.5, and
00_overview.md. Rank-based GSEA was tried first per the user's decision and also came back null,
so the fallback (report the null honestly) was taken.

**Two RIANA bugs (filed above) remain open** in the main repo — not blocking the note.
```
