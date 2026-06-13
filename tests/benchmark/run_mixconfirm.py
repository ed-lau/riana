"""
Mixing-confirm runner: integrate several peak-detection configs on the full
9-proportion calibration series, loading each proportion's mzML **once** and
running all configs against it before moving on (≈ N_proportions mzML loads
instead of N_configs × N_proportions — the load is the bottleneck).

Configs compared (all baseline=none, the spike default), for the cross-proportion
mixing confirm of the 0%/100% screen winners:
  - ms2_w015        : peak_rt=ms2,        integration_half_width=0.15
  - apex_tall_w015  : peak_rt=apex,       apex_selection=tallest, ihw=0.15
  - consensus_w015  : peak_rt=consensus,  apex_n_consensus=4,     ihw=0.15

Outputs land at integrate_outputs/<config>/<sample>_riana.txt; score with
  bench_mixing_linearity.py --method ms2=... --method apex_tall=... --method consensus=...
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from riana.config import IntegrationConfig  # noqa: E402
from riana.core.integration import integrate_run  # noqa: E402
from riana.io.mzml import IndexedMzML  # noqa: E402
from riana.io.percolator import read_percolator  # noqa: E402
from riana.io.writers import make_provenance, write_dataframe_tsv  # noqa: E402

BASE = dict(
    isotopomers=(0, 1, 2, 3, 4, 5),
    q_value=0.01,
    mass_tol_ppm=15,
        forced_mods=(0.0,),
    baseline_method="none",
)

# name -> per-config overrides (extraction_half_width derived per the spike rule:
# = ihw for ms2; = ihw + 0.33 for apex/consensus so the apex can sit off-RT).
CONFIGS = {
    "ms2_w015": dict(peak_rt="ms2", integration_half_width=0.15,
                     extraction_half_width=0.15),
    "apex_tall_w015": dict(peak_rt="apex", apex_selection="tallest",
                           integration_half_width=0.15, extraction_half_width=0.48),
    "consensus_w015": dict(peak_rt="consensus", apex_n_consensus=4,
                           apex_selection="nearest",
                           integration_half_width=0.15, extraction_half_width=0.48),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--line", default="ac16")
    args = parser.parse_args()
    line = args.line

    gt_path = (REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing" / line
               / "ground_truth.csv")
    out_root = (REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing" / line
                / "integrate_outputs")
    out_dirs = {name: out_root / name for name in CONFIGS}
    for d in out_dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    with gt_path.open() as f:
        rows = list(csv.DictReader(f))
    print(f"=== {line.upper()}: {len(rows)} proportions x {len(CONFIGS)} configs "
          f"(load-once-per-proportion) ===", flush=True)

    for row in rows:
        proportion = float(row["nominal_proportion"])
        sample = f"time{proportion:g}"
        mzml_src = REPO_ROOT / "data" / f"calibration_{line}" / "mzml" / row["mzml_filename"]
        psms_path = (REPO_ROOT / "data" / f"calibration_{line}" / "snakemake_results"
                     / sample / "percolator" / "percolator.target.psms.txt")
        pending = {n: out_dirs[n] / f"{sample}_riana.txt" for n in CONFIGS
                   if not (out_dirs[n] / f"{sample}_riana.txt").exists()}
        if not pending:
            print(f"[{sample}] all configs done, skipping", flush=True)
            continue

        print(f"[{sample}] loading mzML + PSMs once ...", flush=True)
        psms = read_percolator(psms_path, sample=sample)
        t0 = time.time()
        with IndexedMzML(mzml_src) as mzml:
            print(f"[{sample}] loaded in {time.time()-t0:.0f}s; integrating "
                  f"{len(pending)} configs", flush=True)
            for name, out_file in pending.items():
                cfg = IntegrationConfig(sample=sample, **{**BASE, **CONFIGS[name]})
                t1 = time.time()
                df = integrate_run(cfg, psms, mzml)
                prov = make_provenance(dataclasses.asdict(cfg), id_source=str(psms_path),
                                       extra={"line": line, "config": name})
                write_dataframe_tsv(out_file, df, prov, include_index=True)
                print(f"    {name:16s} {len(df)} PSMs in {time.time()-t1:.0f}s "
                      f"-> {out_file.name}", flush=True)


if __name__ == "__main__":
    main()
