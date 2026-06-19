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
