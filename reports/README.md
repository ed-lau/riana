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
