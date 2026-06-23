# reports/

Human-readable write-ups of investigations, benchmarks, and design decisions —
the narrative trail behind the terse bullets in `PROJECT_REVIEW.md`. Each report
is a dated, self-contained snapshot of what was measured/decided and why, so the
reasoning survives even when the code or roadmap moves on.

Convention: `YYYY-MM-DD_short_topic.md`. Inputs often live under the gitignored
`runs/` (regenerable); reports quote the numbers so they stand alone.

| date | report | outcome |
|------|--------|---------|
| 2026-06-17 | [MBR DDA feasibility](2026-06-17_mbr_dda_feasibility.md) | DDA curves badly gappy (8–15% complete, 40–45% recoverable, 70% LVE t0-anchor loss); DIA needs none; quantms RT aligned-but-imperfect (15–25 s run residuals) → **GO** for mzTab/DDA MBR with a light per-run RT refinement. |
| 2026-06-17 | [MBR v1 design + validation](2026-06-17_mbr_v1_design.md) | **SHIPPED — gated MBR** (`--mbr-min-snr 4 --mbr-min-scans 3` default, uncapped). Pure RT-transfer + a two-part quality gate (apex-SNR floor with inf=fail + min nonzero scans). Fit A/B (runs/lve_atr_mbr): **ungated MBR HURTS** (R²>0.95 −30%, pollutes clean curves), but **gated MBR is net-positive at the in-vivo R² gates (+180 at R²>0.8), neutral at strict 0.95, no pollution**; cap-fraction guards would hurt (the value is in the high-fraction rescues). `--exclude-mbr` opts out; `n_mbr`/`n_clean` columns audit it. |
| 2026-06-23 | [Adaptive N_ISO & limited-isotopomer scoring](2026-06-23_adaptive_niso_limited_isotopomer.md) | **Fit-side `--fs` limited-isotopomer scoring is the keeper** — a 2×2 capture×scoring control (ac16 calibration) shows iso0-3 scoring tightens recovery (within-±0.05 +2-3pp, lower IQR/bias) while integrate-side adaptive `--iso auto` adds **nothing** and is ~3-5× slower → kept opt-in. H4′ mix-then-normalize (B3) shipped (no-op at low θ, load-bearing for narrow scoring); `--fs` productionized with run-level + per-peptide guards. Aside: the RT↔scan offset is **OpenMS alignment**, not a raw-vs-mzml artifact (mzml re-search does not fix it). |
| 2026-06-23 | [Robust CV for turnover rate constants](2026-06-23_robust_turnover_cv.md) | **Reference / derivation.** `1.4826·MAD(ln k)` is a Fisher-consistent, outlier-resistant estimate of σ_log, which equals the geometric CV to first order (bias <2.5% for CV ≲ 30%, always under-reports; exact `√(e^σ²−1)` past that). Composition of two standard results (Hampel/Rousseeuw-Croux + Limpert) with primary citations + validity bounds. Verified to match `bench_mbr_ab.py`'s within-protein k-CV. |
