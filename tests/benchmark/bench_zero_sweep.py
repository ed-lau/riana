"""
Fast, model-free integration-setting sweep on the 0% (unlabelled) proportion.

The exploration tier for the peak-detection spike (planning round 2026-06-05).
The full 9-proportion mixing metric (bench_mixing_linearity) is the *confirm*
step at ~100 min/setting; this is the *screen*: at 0% the peptide is natural
abundance, so the ground-truth isotopomer envelope is IsoSpec's theoretical
distribution from the formula — **no Spep, no D2O model, no solve_fs**
(forward_model._get_init_env, the f=0 term of bench_m0_ma_recovery's mixing
prediction). Observed-vs-IsoSpec envelope RMSE is then a clean absolute
integration-accuracy metric.

Speed comes from **load-once, integrate-many**: the mzML load + per-PSM XIC
extraction (the real cost, ~minutes) happens ONCE into a cache of generous
traces; each grid setting then only re-runs the cheap window+baseline+integrate
step on the cache. A few-hundred-cell grid costs about one integration.

BLIND SPOT (by design): 0% is the *easiest* regime — strong clean iso0, no
labelled spread, high channels are pure background. So it scores baseline /
contaminant / low-channel accuracy well, but is blind to high-D2O envelope
*capture* (a too-narrow window that would clip the spread labelled envelope at
high D2O looks perfect here). Use this to rank+prune the grid; CONFIRM the
shortlist on bench_mixing_linearity (full 9 proportions) which sees clipping.

Subsample: ~1000 *unmodified* peptides (IsoSpec is exact only on the backbone —
the M7 caveat), stratified length-tertile x abundance-half, fixed + paired
across every setting.

Usage:
  python bench_zero_sweep.py \
    --line ac16 \
    [--n-subsample 1000] [--seed 0] [--n-iso 6] \
    --output-dir <dir>
"""
from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).parent))

from riana.config import IntegrationConfig  # noqa: E402
from riana.core.integration import (  # noqa: E402
    _extract_per_psm, _baseline_corrected_slice, _g,
)
from riana.algorithms import peaks as pk  # noqa: E402
from riana.io.mzml import IndexedMzML  # noqa: E402
from riana.io.percolator import read_percolator  # noqa: E402
from riana.io.mztab import read_mztab  # noqa: E402
from _helpers import forward_model as fm  # noqa: E402
from bench_m0_ma_recovery import estimate_spep  # noqa: E402

# Generous extraction so every grid window fits inside the cached trace.
EXTRACT_HALF_WIDTH = 1.0
# Fallback width (min) when apex/auto detection fails on a PSM.
FALLBACK_HALF_WIDTH = 0.33


def _basename_stem(name: str) -> str:
    """Filename without the MS extension, for matching mzTab runs to an mzML."""
    s = Path(name).name
    for ext in (".mzML.gz", ".mzml.gz", ".mzML", ".mzml", ".raw", ".gz"):
        if s.endswith(ext):
            return s[: -len(ext)]
    return s


def _load_psms(psms_path: Path, fmt: str, sample: str, mzml_name: str):
    """PSMRecords for one run. mzTab is filtered to the ms_run whose location
    basename matches ``mzml_name`` (the in-vivo / quantms path)."""
    if fmt == "mztab":
        records, file_map = read_mztab(psms_path, sample=sample)
        stem = _basename_stem(mzml_name)
        if not any(_basename_stem(n) == stem for n in file_map.values()):
            raise SystemExit(
                f"no mzTab ms_run matches {mzml_name!r}; runs present: "
                f"{[_basename_stem(n) for n in file_map.values()]}")
        return [r for r in records if _basename_stem(r.file_name) == stem]
    return list(read_percolator(psms_path, sample=sample))


@dataclass(frozen=True)
class Setting:
    center: str             # "ms2" | "apex" (iso0) | "consensus" (m0..m{n-1})
    ihw: float | str        # half-width minutes, or "auto"
    baseline: str           # none | noise_floor
    selection: str = "nearest"          # nearest | tallest (apex/consensus)
    apex_search_half_width: float = 0.0  # 0 = whole trace (wide default)
    n_consensus: int = 4                 # channels for consensus (m0..m{n-1})
    width_rel_height: float = 0.05
    prominence_k: float = 3.0

    @property
    def label(self) -> str:
        w = "auto" if self.ihw == "auto" else f"{self.ihw:g}"
        tag = self.center + (str(self.n_consensus) if self.center == "consensus" else "")
        sel = f"_{self.selection}" if self.center in ("apex", "consensus") else ""
        extra = (f"_h{self.width_rel_height:g}_p{self.prominence_k:g}"
                 if self.ihw == "auto" else "")
        return f"{tag}_w{w}{sel}{extra}_{self.baseline}"


def build_grid() -> list[Setting]:
    widths = [0.1, 0.15, 0.2, 0.33]
    grid: list[Setting] = []
    # ms2: fixed window on the PSM RT (no apex selection).
    for w, b in itertools.product(widths, ["none", "noise_floor"]):
        grid.append(Setting("ms2", w, b))
    # apex (iso0) and consensus (m0..m3), nearest vs tallest, baseline none.
    # apex_search wide (0.0) for this first pass — tightening it is a later axis.
    for c, w, sel in itertools.product(["apex", "consensus"], widths,
                                       ["nearest", "tallest"]):
        grid.append(Setting(c, w, "none", selection=sel))
    # auto (iso0 detect_peak) kept for reference: 5% vs FWHM.
    for h in (0.05, 0.5):
        grid.append(Setting("apex", "auto", "none", width_rel_height=h))
    return grid


def build_subsample(zero_riana: Path, n: int, seed: int) -> pd.DataFrame:
    """Stratified ~n unmodified peptides: length-tertile x abundance-half."""
    df = pd.read_csv(zero_riana, sep="\t", comment="#")
    iso_cols = [c for c in ("iso0", "iso1", "iso2", "iso3", "iso4", "iso5")
                if c in df.columns]
    df = df.rename(columns={f"m{i}": f"iso{i}" for i in range(6)})
    iso_cols = ["iso0", "iso1", "iso2", "iso3", "iso4", "iso5"]
    df["total_int"] = df[iso_cols].sum(axis=1)
    df = df[df["total_int"] > 0].copy()
    df["sequence"] = df["concat"].str.rsplit("_", n=1).str[0]
    # unmodified only — exact IsoSpec ground truth (M7 caveat)
    df = df[~df["sequence"].str.contains(r"[\[\(]", regex=True)]
    df = df.drop_duplicates("concat")
    df["pep_len"] = df["sequence"].str.len()
    df["len_bin"] = pd.qcut(df["pep_len"], 3, labels=["short", "mid", "long"])
    df["abund_bin"] = pd.qcut(df["total_int"], 2, labels=["low", "high"])
    rng = np.random.default_rng(seed)
    per_cell = max(1, n // 6)
    picks = []
    for _, g in df.groupby(["len_bin", "abund_bin"], observed=True):
        take = min(per_cell, len(g))
        picks.append(g.sample(take, random_state=rng.integers(1 << 31)))
    sub = pd.concat(picks, ignore_index=True)
    return sub[["concat", "sequence", "pep_len", "len_bin", "abund_bin",
                "total_int"]]


def _window(rt_arr, traces, psm_rt, s: Setting):
    """Return (lo, hi) integration bounds for one setting on one cached PSM."""
    n = rt_arr.size
    iso0 = traces[0]
    if s.ihw == "auto":
        b = pk.detect_peak(rt_arr, iso0, scan_prior_rt=psm_rt,
                           rel_height=s.width_rel_height, prominence_k=s.prominence_k)
        if b is not None:
            return b.lo, b.hi
        centre, half = psm_rt, FALLBACK_HALF_WIDTH      # detection failed
    elif s.center == "apex":
        a = pk.find_apex(rt_arr, iso0, scan_prior_rt=psm_rt,
                         prominence_k=s.prominence_k,
                         apex_search_half_width=s.apex_search_half_width,
                         selection=s.selection)
        centre = rt_arr[a] if a is not None else psm_rt
        half = float(s.ihw)
    elif s.center == "consensus":
        res = pk.consensus_apex(rt_arr, traces[:s.n_consensus], scan_prior_rt=psm_rt,
                                prominence_k=s.prominence_k,
                                apex_search_half_width=s.apex_search_half_width,
                                selection=s.selection)
        centre = rt_arr[res[0]] if res is not None else psm_rt
        half = float(s.ihw)
    else:  # ms2 fixed
        centre, half = psm_rt, float(s.ihw)
    win = np.where(np.abs(rt_arr - centre) <= half)[0]
    if win.size == 0:
        return 0, n - 1
    return int(win[0]), int(win[-1])


def score_setting(cache: list[tuple], iso_cols, s: Setting, n_iso: int) -> pd.DataFrame:
    rows = []
    for concat, meta, rt_arr, traces, psm_rt, init_norm in cache:
        lo, hi = _window(rt_arr, traces, psm_rt, s)
        areas = np.empty(n_iso)
        for k in range(n_iso):
            corrected = _baseline_corrected_slice(rt_arr, traces[k], lo, hi, s.baseline)
            areas[k] = max(0.0, float(np.trapezoid(corrected, x=rt_arr[lo:hi + 1])))
        tot = areas.sum()
        if tot <= 0:
            continue
        obs = areas / tot
        err = obs - init_norm
        rows.append({
            "concat": concat, "len_bin": meta["len_bin"],
            "abund_bin": meta["abund_bin"], "pep_len": meta["pep_len"],
            "env_rmse": float(np.sqrt(np.mean(err ** 2))),
            **{f"bias_m{k}": float(err[k]) for k in range(n_iso)},
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--line", default="ac16")
    parser.add_argument("--n-subsample", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--n-iso", type=int, default=6)
    parser.add_argument("--proportion", type=float, default=0.0,
                        help="D2O mixing proportion to sweep. 0 = natural (no "
                             "Spep); >0 builds the labelled ground truth via the "
                             "mixing model and needs --coefficients.")
    parser.add_argument("--coefficients", type=Path, default=None,
                        help="frozen d2o_aa_coefficients_<line>.csv (Spep source "
                             "when --proportion > 0)")
    parser.add_argument("--mzml", type=Path, default=None,
                        help="explicit mzML path (e.g. in-vivo); bypasses the "
                             "calibration-structure resolution")
    parser.add_argument("--psms", type=Path, default=None,
                        help="explicit PSM file (percolator .txt or mzTab)")
    parser.add_argument("--psm-format", choices=["percolator", "mztab"],
                        default="percolator")
    parser.add_argument("--subsample-riana", type=Path, default=None,
                        help="riana.txt to draw the subsample from (required with --mzml)")
    parser.add_argument("--sample", default=None,
                        help="sample label for explicit mode (must end in a number)")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    n_iso = args.n_iso

    prop = args.proportion
    f_frac = prop / 100.0
    if args.mzml is not None:                  # explicit-path mode (in-vivo / mzTab)
        if args.psms is None or args.subsample_riana is None:
            parser.error("--mzml mode requires --psms and --subsample-riana")
        sample = args.sample or "time0"
        mzml_src = args.mzml
        psms_path = args.psms
        prop_riana = args.subsample_riana
    else:                                      # calibration-structure mode
        sample = f"time{prop:g}"
        line_dir = REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing" / args.line
        prop_riana = line_dir / "integrate_outputs" / "v1.0.0" / f"{sample}_riana.txt"
        gt = pd.read_csv(line_dir / "ground_truth.csv")
        row = gt[gt["nominal_proportion"] == prop].iloc[0]
        mzml_src = (REPO_ROOT / "data" / f"calibration_{args.line}" / "mzml"
                    / row["mzml_filename"])
        psms_path = (REPO_ROOT / "data" / f"calibration_{args.line}"
                     / "snakemake_results" / sample / "percolator"
                     / "percolator.target.psms.txt")
    coeff_dict = None
    if f_frac > 0:
        if args.coefficients is None:
            parser.error("--coefficients is required when --proportion > 0")
        cdf = pd.read_csv(args.coefficients)
        coeff_dict = dict(zip(cdf["amino_acid"], cdf["coefficient"]))

    print(f"[subsample] proportion={prop:g} fmt={args.psm_format}  from {prop_riana}",
          flush=True)
    sub = build_subsample(prop_riana, args.n_subsample, args.seed)
    want = set(sub["concat"])
    meta_by_concat = sub.set_index("concat").to_dict("index")
    print(f"[subsample] {len(sub)} unmodified peptides "
          f"({sub.groupby(['len_bin','abund_bin'], observed=True).size().to_dict()})",
          flush=True)

    # --- one-time load + extraction into the cache ---
    extract_cfg = IntegrationConfig(
        sample=sample, isotopomers=tuple(range(n_iso)), q_value=0.01,
        mass_tol_ppm=15, extraction_half_width=EXTRACT_HALF_WIDTH, use_range=False,
        smoothing=None, forced_mods=(0.0,),
        peak_rt="ms2", integration_half_width=EXTRACT_HALF_WIDTH,
        baseline_method="none",
    )
    iso_cols = [f"mod{_g(0.0)}_iso{n}" for n in range(n_iso)]
    psms = [p for p in _load_psms(psms_path, args.psm_format, sample, mzml_src.name)
            if p.concat in want]
    seen: set[str] = set()
    cache = []
    t0 = time.time()
    with IndexedMzML(mzml_src) as mzml:
        for psm in psms:
            if psm.concat in seen:
                continue
            seen.add(psm.concat)
            idf, _ma = _extract_per_psm(psm, {}, (0.0,), tuple(range(n_iso)),
                                        extract_cfg, mzml)
            rt_arr = idf["rt"].to_numpy(dtype=np.float64)
            traces = [idf[c].to_numpy(dtype=np.float64) for c in iso_cols]
            psm_rt = float(
                mzml.rt_idx[np.searchsorted(mzml.scan_idx, psm.scan, "left") - 1])
            # Anchor neutral (charge=0) to match psm.peptide_mass (neutral
            # monoisotopic) — same convention as the new engine's isotope_dist.
            # Mixing a charged dist with a neutral pep_mass shifts the bins by
            # +z Da (the harness bug that inflated RMSE to ~0.3). The envelope
            # *shape* is charge-independent, so charge=0 is exact here.
            init = fm._get_init_env(psm.sequence, 0, psm.peptide_mass, n=n_iso)
            isum = init.sum()
            if isum <= 0:
                continue
            pred = init / isum
            if f_frac > 0:
                # Labelled component via the mixing model (= _get_final_env at
                # proportion 100). Spep from the frozen coefficient table; neutral
                # (charge=0) anchoring as for init.
                spep = estimate_spep(psm.sequence, coeff_dict)
                final = fm._get_final_env(psm.sequence, 0, psm.peptide_mass,
                                          spep, n=n_iso)
                fsum = final.sum()
                if fsum <= 0:
                    continue
                pred = (1.0 - f_frac) * pred + f_frac * (final / fsum)
            cache.append((psm.concat, meta_by_concat[psm.concat], rt_arr,
                          traces, psm_rt, pred))
            if len(cache) % 200 == 0:
                print(f"  extracted {len(cache)} in {time.time()-t0:.0f}s", flush=True)
    print(f"[cache] {len(cache)} peptides extracted in {time.time()-t0:.0f}s", flush=True)

    # --- sweep on the cache ---
    grid = build_grid()
    print(f"[sweep] {len(grid)} settings", flush=True)
    summary_rows = []
    per_channel_rows = []
    for i, s in enumerate(grid):
        sc = score_setting(cache, iso_cols, s, n_iso)
        if not len(sc):
            continue
        rmse = sc["env_rmse"].to_numpy()
        row = {"setting": s.label, "center": s.center, "selection": s.selection,
               "ihw": s.ihw, "baseline": s.baseline,
               "width_rel_height": s.width_rel_height, "prominence_k": s.prominence_k,
               "n": len(sc), "rmse_median": float(np.median(rmse)),
               "rmse_mean": float(np.mean(rmse)),
               "rmse_p90": float(np.percentile(rmse, 90))}
        # per-stratum median RMSE
        for (lb, ab), g in sc.groupby(["len_bin", "abund_bin"], observed=True):
            row[f"rmse_{lb}_{ab}"] = float(np.median(g["env_rmse"]))
        summary_rows.append(row)
        per_channel_rows.append({"setting": s.label, **{
            f"bias_m{k}": float(np.median(sc[f"bias_m{k}"])) for k in range(n_iso)}})
        if (i + 1) % 10 == 0:
            print(f"  scored {i+1}/{len(grid)} settings", flush=True)

    summ = pd.DataFrame(summary_rows).sort_values("rmse_median")
    summ.to_csv(args.output_dir / "zero_sweep_summary.csv", index=False)
    pd.DataFrame(per_channel_rows).to_csv(
        args.output_dir / "zero_sweep_per_channel.csv", index=False)
    with (args.output_dir / "zero_sweep.json").open("w") as f:
        json.dump({"line": args.line, "proportion": prop,
                   "n_peptides": len(cache),
                   "n_settings": len(grid), "seed": args.seed,
                   "top": summ.head(15).to_dict("records")}, f, indent=2)

    print("\n[done] top 12 settings by median envelope RMSE vs IsoSpec natural:")
    cols = ["setting", "n", "rmse_median", "rmse_p90",
            "rmse_short_low", "rmse_long_low"]
    cols = [c for c in cols if c in summ.columns]
    print(summ[cols].head(12).to_string(index=False))
    print(f"\n       baseline of record (ms2_w0.33_none): "
          f"{summ.loc[summ['setting']=='ms2_w0.33_none','rmse_median'].values}")
    print(f"       outputs -> {args.output_dir}")


if __name__ == "__main__":
    main()
