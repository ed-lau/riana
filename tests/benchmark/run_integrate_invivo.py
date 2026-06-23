"""
Integrate one in-vivo acquisition from a quantms mzTab → a ``_riana.txt`` that
serves as the subsample source for ``bench_zero_sweep.py --mzml`` (Part B of the
generalizability check). Used for ``timeseries_lve`` time0 (= day-0, unlabelled
→ natural-abundance IsoSpec ground truth, no Spep, like proportion 0).

PSMs come from the OpenMS/quantms mzTab, filtered to the ms_run whose location
matches the given mzML (via ``bench_zero_sweep._load_psms``). Config is the
spike default (ms2 / 0.15 / baseline none) — the integrated areas only seed the
stratified subsample; ``bench_zero_sweep`` re-extracts the XICs itself.

Usage:
  python run_integrate_invivo.py \
    --mzml  data/timeseries_lve/mzml/20230907_JC_Boulder_LVE_time0.mzML \
    --mztab data/timeseries_lve/quantms_results/quant_tables/samplesheet_lve_sdrf_openms_design_openms.mzTab \
    --out   tests/data/timeseries_lve/time0_riana.txt
"""
from __future__ import annotations

import argparse
import dataclasses
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).parent))

from riana.config import IntegrationConfig  # noqa: E402
from riana.core.integration import integrate_run  # noqa: E402
from riana.io.mzml import IndexedMzML  # noqa: E402
from riana.io.writers import make_provenance, write_dataframe_tsv  # noqa: E402
from bench_zero_sweep import _load_psms  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mzml", type=Path, required=True)
    ap.add_argument("--mztab", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sample", default="time0", help="must end in a number")
    args = ap.parse_args()

    psms = _load_psms(args.mztab, "mztab", args.sample, args.mzml.name)
    print(f"[psms] {len(psms)} for run {args.mzml.name}", flush=True)
    if not psms:
        raise SystemExit("no PSMs matched the mzML run in the mzTab")

    cfg = IntegrationConfig(
        sample=args.sample, isotopomers=(0, 1, 2, 3, 4, 5), q_value=0.01,
        mass_tol_ppm=15,
        peak_rt="ms2", integration_half_width=0.15,
        extraction_half_width=0.15, baseline_method="none",
    )
    t0 = time.time()
    with IndexedMzML(args.mzml) as mzml:
        df = integrate_run(cfg, psms, mzml)
    print(f"[integrate] {len(df)} PSMs in {time.time()-t0:.0f}s", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    prov = make_provenance(dataclasses.asdict(cfg), id_source=str(args.mztab),
                           extra={"dataset": "timeseries_lve", "run": args.mzml.name})
    write_dataframe_tsv(args.out, df, prov, include_index=True)
    print(f"[done] -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
