#!/usr/bin/env bash
# Regenerates every Monte Carlo number in reports/2026-11-12_deltak_inference_manual.md.
#
#   bash tests/benchmark/run_deltak_mc.sh            # ~1.5 h on 24 cores
#   W=8 OUT=runs/my_mc bash tests/benchmark/run_deltak_mc.sh
#
# Each arm writes <arm>/ratio<r>_params.json (settings, commit, tree state) and
# <arm>/ratio<r>_results.tsv (the shipped CLI's rejection rates); the evaluator then
# re-tests the same collapsed cells under the collapse / weight variants and writes
# $OUT/evaluation.tsv. Seeds are fixed (bench default 20261112 + dataset + 1000·ratio).
set -euo pipefail
W=${W:-22}
OUT=${OUT:-runs/deltak_mc}
PY=${PY:-python}
B() { "$PY" -m tests.benchmark.bench_rising_null -W "$W" --proteins 250 "$@"; }

# A. Slow precursor (k_p 0.08/day), full nuisance: null + the ratio grid.
B --work "$OUT/slow" --ratio 1 --datasets 3
for r in 1.25 1.5 2; do B --work "$OUT/slow" --ratio "$r" --datasets 1; done

# B. Instant precursor (k_p 100/day ≈ priming dose), full nuisance: null + one effect size.
B --work "$OUT/instant" --kp 100 --ratio 1 --datasets 3
B --work "$OUT/instant" --kp 100 --ratio 1.5 --datasets 1

# C. Instant precursor, null, 2×2 knock-out of the two nuisance families
#    (per-peptide labeling-site error; per-animal turnover + enrichment). Interference and
#    per-peptide precision stay on in every arm.
B --work "$OUT/instant_nosite" --kp 100 --site-err 0 --datasets 1
B --work "$OUT/instant_noanimal" --kp 100 --animal-k-sd 0 --ria-cv 0 --datasets 1
B --work "$OUT/instant_neither" --kp 100 --site-err 0 --animal-k-sd 0 --ria-cv 0 --datasets 1

"$PY" -m tests.benchmark.eval_rising_variance "$OUT"
