# Multi-point FS rail-drop — re-validation A/B (boomi / juber / lve, D₂O + ¹⁸O)

- **Status:** **SHIPPED — default ON.** The single-timepoint FS rail-drop
  (drop per-timepoint FS points the solver clamped to its ±`FS_BOUNDS` diagnostic rail
  before counting fit points / depth) is extended to **multi-timepoint** fits with the
  *identical* criterion, gated by a new `FitConfig.fs_rail_drop` (CLI
  `--fs-rail-drop / --no-fs-rail-drop`, default on). Re-validated against the five
  standard turnover sets at the production **depth-6** curation.
- **Data:** `runs/{boomi_ipsc_d2o, juber_ac16_d2o, lve_atr_clean, boomi_ipsc_o18,
  juber_ac16_o18}` (coefficients: `alamillo_2025_{ipsc,ac16}`, `deberneh_2025_rss`,
  `juber_2026_o18_ac16`). Harness: `tests/benchmark/bench_fs_rail_drop.py` (in-process `fit_project`,
  `n_boot=200`, production rollup gate R²≥0.8 OR (R²≥0.6 ∧ k_cv<0.2), geomCV = median
  over proteins with ≥3 admitted unique-accession peptides).

## TL;DR

- **Default ON is correct at the real curation depth.** At `--depth 6`, across all five
  sets: admitted-peptide **yield rises** where rail-hits are common (+7–14% on AC16/iPSC),
  is neutral on the cleanest in-vivo set (lve −2.5%); **R² of the fitted population is
  universally cleaner** (+0.02 to +0.36); the **matched within-protein geom-CV**
  (same peptides admitted in both arms) is **better or exactly neutral on all five**; and
  the **median k is unbiased** everywhere (|Δk| ≤ 0.001 /unit). A rail-hit genuinely is a
  failed solve, not a measurement — dropping it recovers real peptides and tightens the
  ones already there without moving the answer.
- **Depth matters — the depth-3 "regression" was an artifact.** At the loose `--depth 3`
  default, rail-drop *appeared* to worsen within-protein CV (up to +0.048 on AC16). That
  was the loose floor letting rail-drop thin a curve to 2–3 points → noisy k. At the
  production depth-6 curation the curves stay well-sampled (matched median n_points
  17→16, 8→8), so the CV effect disappears and reverses. **Always re-validate curation
  changes at the curation depth actually used.**
- **The matched-CV is the honest metric.** The *raw* admitted geomCV is confounded — the
  admitted set grows, so which proteins clear the ≥3-peptide bar shifts. On boomi ¹⁸O the
  raw geomCV rose +0.008 while the matched-CV *improved* −0.0045. Compare the same
  peptides.

## Depth-6 results (ON − OFF, production gate)

| run | Δ admitted yield | Δ R²(all fitted) | Δ matched geom-CV | Δ median k |
|---|---|---|---|---|
| boomi_ipsc_d2o | +484 (+7.2%) | +0.073 | **−0.0064** | −0.0004 |
| juber_ac16_d2o | +146 (+13.9%) | +0.340 | **0.0000** | −0.0004 |
| lve_atr_clean  | −246 (−2.5%) | +0.023 | **−0.0057** | −0.0010 |
| boomi_ipsc_o18 | +339 (+5.0%) | +0.050 | **−0.0045** | −0.0004 |
| juber_ac16_o18 | +79 (+14.0%) | +0.355 | **−0.0000** | −0.0001 |

Reading: **yield up on 4/5, neutral on lve; R² of the whole fitted population cleaner on
5/5; matched within-protein agreement better on 3, exactly unchanged on 2; k unbiased on
5/5.** The newly-rescued admits are not noise — their standalone geom-CV is *tighter* than
the matched set on boomi (0.195 vs 0.210) and iPSC ¹⁸O (0.199 vs 0.209).

### The depth-3 vs depth-6 contrast (why depth is load-bearing)

At `--depth 3` the same A/B showed a within-protein CV rise (boomi +0.010, AC16 +0.048).
Fitting at the production depth-6 floor:

| run | Δ matched-CV @ depth 3 (est. from raw) | Δ matched-CV @ depth 6 |
|---|---|---|
| juber_ac16_d2o | +0.048 (raw, worst) | **0.0000** |
| boomi_ipsc_d2o | +0.010 (raw) | **−0.0064** |
| lve_atr_clean  | −0.010 (raw) | **−0.0057** |

The mechanism: rail-drop removes ~1 median point from a well-sampled curve (17→16, 8→8),
which is harmless at depth 6 but can push a depth-3 curve to the 2-point floor, where k is
unstable. The pathology was the loose floor, not the rail-drop.

## Decision

- **`fs_rail_drop` defaults ON** for all fits (single- and multi-timepoint), unified under
  one flag. `--no-fs-rail-drop` reproduces the pre-1.2.0 multi-timepoint numbers.
- The k estimate is unbiased, so downstream Δk / linear-φ contrasts are unaffected in
  central value; they gain the extra admitted peptides and the cleaner R² population.

## Drop-threshold sweep (is clamp-only ideal, or a tighter bound?)

The clamp-only drop (the solver's ±`FS_BOUNDS`) only removes points the solver *pinned*
to its bound. A solved-but-implausible FS (1.15 over-labelled, −0.07 sub-natural) is also
not a real measurement. Sweep of the drop threshold (rail-drop ON, `--depth 6`,
production gate; `tests/benchmark/bench_fs_rail_threshold.py`):

| threshold (lo / hi) | boomi Δyield / r2all | juber_ac16 Δyield / r2all | lve Δyield / r2all |
|---|---|---|---|
| clamp (−0.099 / 1.199) | — / 0.760 | — / 0.343 | — / 0.925 |
| −0.1 / 1.1 | −6.4% / 0.690 | −12.7% / 0.023 | +1.2% / 0.910 |
| **−0.05 / 1.05** | **+4.7% / 0.789** | **+9.4% / 0.361** | **−4.5% / 0.933** |
| 0.0 / 1.0 | +5.7% / 0.803 | −16.0% / −0.101 | −11.0% / 0.935 |

- **`−0.1 / 1.1` is an artifact, not a real config:** its `lo=−0.1` sits *below* the
  solver's low-clamp convergence (≈ −0.09999), so it stops dropping the lower rail *and*
  cuts real upper plateau — worse on 2 of 3. It only exists here to show the lower
  threshold must stay **above** the clamp convergence.
- **`0.0 / 1.0` (the exact physical bound) is dataset-dangerous:** best on clean boomi,
  but its `lo=0.0` eats AC16's real near-zero t0/early anchors (incomplete turnover, low
  FS scattered slightly negative) → **−16% yield, R² −0.101**. Rejected.
- **`−0.05 / 1.05` (a 0.05 margin around [0,1]) is the robust optimum:** cleaner R² and
  matched within-protein CV on all three, +5–9% yield on the noisy sets (AC16, iPSC),
  and on clean lve it trims the 4.5% *noisiest* peptides (matched-CV −0.0012, i.e. the
  retained peptides agree *better*). It removes solved-but-implausible points without
  eating genuine plateau / anchors.

**Decision: default rails = `1.05 / −0.05`** (`_FS_RAIL_HI` / `_FS_RAIL_LO`), overridable
per fit via `FitConfig.fs_rail_hi` / `fs_rail_lo`. `fs_rail_drop=False` disables the drop
entirely.

### o18 and single-timepoint (TMT) confirmation

The sweep above is D₂O + multi-timepoint. Before adopting `1.05 / −0.05` as the default
it was confirmed non-harmful on the other two regimes (`tests/benchmark/bench_fs_rail_threshold.py` for ¹⁸O,
`tests/benchmark/bench_fs_rail_singlepoint.py` for TMT):

- **¹⁸O (multi-timepoint) — strict no-op.** `1.05 / −0.05` reproduces the clamp result
  *exactly* on both ¹⁸O sets (boomi admit 7241 = 7241; juber admit 719 = 719; Δmatched-CV
  0.0000). The ¹⁸O solve lands either cleanly below 1.05 or fully railed at the 1.2 clamp
  (caught by both thresholds), with nothing in the intermediate (1.05, 1.199] band — so
  the tighter rail changes nothing for ¹⁸O.

  | ¹⁸O run | clamp admit | 1.05/−0.05 admit | Δ |
  |---|---|---|---|
  | boomi_ipsc_o18 | 7241 | 7241 | 0 |
  | juber_ac16_o18 | 719 | 719 | 0 |

- **Single-timepoint TMT (splatd2o_tmt16, 16-plex AC16 @ 24 h) — neutral.** Single-tp
  curation rides on the replicate floor (`min_fit_points` auto-2) + `k_cv` (R² bypassed),
  so a tighter rail *could* thin bioreps below the 2-replicate floor. It barely does:
  vs clamp, admitted peptides −0.5 % (5299 → 5271), proteins −2 (914 → 912), and median k
  is **identical** (0.0194 → 0.0194) — a tiny trim of a few over-labelled bioreps, no
  yield hit and no k bias. (Rail-drop *on* vs *off* on this set drops the count more —
  5504 → 5299 admitted — which is the intended removal of fake-replicate peptides whose
  railed FS manufactured a spurious `k_cv = 0`; that is the point of the drop, not a
  regression.)

  | TMT arm | admit (k_cv<0.2) | proteins | median k |
  |---|---|---|---|
  | rail-drop OFF | 5504 | 955 | 0.0196 |
  | ON, clamp | 5299 | 914 | 0.0194 |
  | ON, 1.05/−0.05 | 5271 | 912 | 0.0194 |


## Reproduce

```
python -m tests.benchmark.bench_fs_rail_drop                 # depth-6 A/B + matched-CV, all 5 sets
python -m tests.benchmark.bench_fs_rail_threshold            # drop-threshold sweep, 3 D₂O sets
python -m tests.benchmark.bench_fs_rail_threshold --runs boomi_ipsc_o18 juber_ac16_o18  # ¹⁸O no-op
python -m tests.benchmark.bench_fs_rail_singlepoint          # single-timepoint TMT
```
