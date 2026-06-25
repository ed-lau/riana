# Mass-defect θ (Δspacing) — a drift-robust second turnover estimate

- **Date:** 2026-06-25
- **Branch:** `1.1.0`
- **Status:** investigation complete — **1a Δmass QC shipped & verified; 1b design decided + prototyped.** θ_ΔS is a real, drift-robust, well-calibrated-in-aggregate but **low-precision** estimator → keep it as a **cross-check / second estimate**, never displacing the intensity FS (θ_ΔI). Production wiring (θ_ΔS column + GUI anchor toggle) next.
- **Prototype:** `tests/benchmark/bench_mass_defect_theta.py` (standalone, not yet production-wired)
- **Inputs:** `runs/calib_ac16_v1` (AC16 D₂O mixing, ground-truth proportion), `runs/lve_fixed_ab` (in-vivo LVE turnover); o18 cross-check `runs/calib_o18`
- **Roadmap / memory:** PROJECT_REVIEW Track C; memory `mass_defect_theta_design` (item #1)

## Question

D₂O (and ¹⁸O) labelling shifts each neutromer's **accurate mass**, not just its
abundance — a second, orthogonal turnover signal (the DeuteRater method, Naylor/Price
*Bioinformatics* 2017) that Riana recorded (`iso{N}_obs_mz`) but never used. Build the
QC display (1a), then decide whether the mass-defect θ (θ_ΔS) earns a place beside the
intensity FS (θ_ΔI) — and how to compute it (1b).

## TL;DR

- **θ_ΔS works, but as a *cross-check*, not a replacement.** On ground-truth calibration its
  population *median* tracks the known mixing proportion near-perfectly, yet per-point it is
  **~2.4× noisier than θ_ΔI** (r 0.35 vs 0.82). So: θ_ΔI stays primary; θ_ΔS is a drift-robust
  second estimate whose *disagreement* with θ_ΔI is the useful product.
- **Empirical f0/t0 anchoring is essential** — it removes a per-peptide model-vs-measurement
  reference offset and **2–4× the usable correlation** (calib r 0.075 → 0.348; LVE r 0.193 →
  0.343). Theory-referenced ΔS alone is near-useless per point.
- **Use iso0–3 only** (weighted median + MAD); keep iso4/5 in the GUI for QC. The high channels'
  apparent noise is **not** low SNR (the labelled envelope broadens 4–6× — iso4 to ~11% at f=1)
  but a per-channel reference offset that grows with index.
- **The Δspacing metric is verified correct**: computed `dmass_iso0` matches the recorded
  `iso0_ppm_error` to <0.001 mDa, and the o18 cross-check reproduces the expected physics
  (iso1 flat ≈0, iso2 negative — ¹⁸O is ~2.5 mDa *lighter* than the ¹³C₂ it displaces).

## What 1a ships

Per fitted (peptide, timepoint) × channel, in mDa (`core/fitting._fit_one_concat`,
emitted to `riana_fit_fractions.txt`; GUI Model-tab **Fit / Δspacing / Δmass** toggle):

- **Δmass** = `obs_mz(k) − init_ref(k)` — absolute shift (drift-sensitive).
- **Δspacing (ΔSₓ)** = `(obs(k)−obs(0)) − (ref(k)−ref(0))` — M0-internal (drift-robust); `ref`
  = IsoSpec natural-abundance averaged-isotopolog mass, recomputed at fit (path-independent).

## Results

### 1. The metric is correct (and label-general)

`dmass_iso0` vs the recorded `iso0_ppm_error` (independent integrate-side computation):
agreement **< 0.001 mDa** over 50 peptides. The **o18** calibration (`runs/calib_o18`,
`o18_ac16`) reproduces the expected ¹⁸O physics: ΔS_iso1 flat ≈0 (the +1 channel is pure
¹³C, no ¹⁸O; Spearman 0.0), ΔS_iso2 develops a clear **negative** shift with labelling
(¹⁸O at +2.0042 Da is ~2.5 mDa lighter than ¹³C₂ at +2.0067) — opposite-signed to D₂O and
~20× smaller (~0.1 mDa). Reproducing the right per-channel physics on a second chemistry
is strong evidence the computation is sound.

### 2. M0-internal cancels per-spectrum drift; a per-peptide reference offset remains

Comparing the Δmass vs Δspacing views of the same peptide: in Δmass all channels move
**together** (per-run m/z calibration drift); Δspacing removes that common mode (M0-internal),
leaving the genuine per-channel signal. But Δspacing does **not** sit at 0 at f=0 for every
peptide — a per-peptide, channel-growing **reference offset** persists (calib_ac16, ΔS at f=0,
4776 peptides):

| channel | median (mDa) | MAD |
|---|---|---|
| iso1 | −0.03 | 0.28 |
| iso2 | +0.05 | 0.49 |
| iso3 | +0.16 | 0.92 |
| iso4 | +0.49 | 1.58 |
| iso5 | +1.60 | 2.57 |

iso1/iso2 do centre on 0 (intuition holds); iso3–5 carry a systematic offset. It is a stable
**model-vs-measurement** bias — the IsoSpec averaged-isotopolog centroid over a ±0.5 bin vs the
instrument's intensity-weighted centroid — **not** a drift artifact (a fit pooling 9 samples,
each with its own file drift, recovers an intercept matching the single measured f0 offset to
~0.2–0.37 mDa) and **not** low SNR (the labelled envelope broadens 4–6×: Spep-18 peptide iso4
2.8%→**11.3%**, iso5 0.8%→**5.1%** at f=1; f0-anchored scatter does *not* shrink as those
channels gain intensity).

### 3. Empirical f0/t0 anchoring is essential

Anchoring each peptide's ΔS to its own unlabelled (f0/t0) spacing subtracts the offset exactly
(the theoretical reference cancels: `ΔSₓ_anchored = [obs(f,k)−obs(f,0)] − [obs(f0,k)−obs(f0,0)]`).
Prototype θ_ΔS = weighted-median over iso1–3 of `ΔSₓ/ΔSₓmax` (MAD guard):

| dataset | metric | theory-only | **f0-anchored** |
|---|---|---|---|
| calib_ac16 (known f) | r(θ_ΔS, f) | 0.075 | **0.348** |
| calib_ac16 | r(θ_ΔI, f) | — | 0.822 |
| lve_fixed_ab (turnover) | r(θ_ΔS, θ_ΔI) | 0.193 | **0.343** |
| lve_fixed_ab | agreement MAD | 0.30 | 0.26 |

Anchoring needs a measured unlabelled point — free for calibration, and an argument to **acquire
t0** for turnover (where it is canonically skipped). When t0 is absent, fall back to theory-
referenced ΔS on iso0–3 (the offset there is small: median +0.05/+0.16 mDa at iso2/iso3).

### 4. θ_ΔS is calibrated-in-aggregate but low-precision (→ cross-check)

calib_ac16 median θ by known proportion (f0-anchored):

| f | median θ_ΔS | median θ_ΔI (fs) |
|---|---|---|
| 0.00 | 0.00 | −0.03 |
| 0.25 | 0.30 | 0.22 |
| 0.50 | 0.58 | 0.46 |
| 0.75 | 0.82 | 0.73 |
| 1.00 | 0.98 | 0.99 |

θ_ΔS's median tracks the diagonal at least as well as θ_ΔI — in fact **less biased mid-range**
(θ_ΔI undershoots, 0.46 at f=0.5). But per-point precision is the opposite story (r 0.35 vs 0.82):
even a single clean channel's linear-R² vs known f is only 0.6–0.7 (whole-range signal ~1.3–1.6
mDa ≈ per-point noise). So θ_ΔS is **only usable as the iso0–3 weighted-median**, and its value is
as an orthogonal cross-check — its mid-range disagreement with θ_ΔI is itself a finding.

## Decision (1b)

- **θ_main = intensity FS (θ_ΔI)** — primary, written, unchanged.
- **θ_ΔS = experimental second estimate / drift-robust cross-check** — per-sample, M0-internal,
  **iso0–3**, weighted median + MAD, normalized by IsoSpec ΔSₓmax.
- **Reference = empirical f0/t0 anchor when available (toggle), theory-referenced fallback** —
  chosen so the *written per-timepoint* θ_ΔS is unbiased and point-comparable to θ_ΔI (fit-time
  intercept absorption would fix only the slope and force two fit models).
- **GUI**: keep iso0–5 in the Δ views (iso4/5 = QC); add an anchor-to-f0/t0 toggle (off = the
  absolute mass-accuracy QC that surfaced the offset).

## Reproduce

```bash
# fit a calibration manifest (auto-dispatches the calibration recovery model) → fractions
riana fit --manifest runs/calib_ac16_v1/riana_manifest.tsv --coefficients ac16 --ria 0.06
#   (turnover: riana fit runs/lve_fixed_ab/*_riana.txt --coefficients commerford --ria 0.045)

# θ_ΔS prototype vs θ_ΔI and known f; add --no-anchor to A/B the f0 anchor
python -m tests.benchmark.bench_mass_defect_theta \
  runs/calib_ac16_v1/riana_fit_fractions.txt --coefficients ac16 --ria 0.06
# calibration mzML/PSMs are gitignored; regenerate integrate via the SDRF path.
```
