"""
M3 pre-M4 peak-detection spike: model-free mixing-linearity sweep metric.

The Tier-1 gate for the peak-detection / baseline revisit (planning round
2026-06-04). Unlike bench_m0_ma_recovery, this metric uses **no IsoSpec
forward model, no Spep, no frozen coefficient table** — only the integrated
isotopomer areas and the known nominal mixing proportion f.

Why it is clean (derivation in the planning notes): in a two-pool mix of a
fully-labelled lysate (heavy fraction f) with an unlabelled lysate, the
observed unnormalised channel intensities are exactly linear superpositions
of the two fixed endpoint envelopes, so each normalised fraction

    m_i:mA(f) = iso_i / sum(iso0..iso5)

is a smooth *monotonic* function of f (a Mobius function, very nearly linear
at 6% D2O where little envelope spills past the iso5 truncation). A perfect
integrator therefore puts every peptide's m_i:mA on a clean line vs f; the
per-peptide line-fit R^2 and residual measure integration noise. The exact
Mobius curvature is a per-peptide *constant*, so it cancels when comparing
two integration methods on the same peptides — which is all this sweep does.

GATE vs PROBE (the role split):
  - GATE channels m0:mA, m1:mA, m2:mA — high-SNR, low-contaminant. This is
    where the real science (m0/m1, m0/m2 -> IsoSpec) lives, so it decides
    GO/NO-GO. Reported as mean/median R^2, fraction crossing R^2>0.95, and
    residual RMSE, on curated and uncurated populations.
  - PROBE channels m4:mA, m5:mA — low-abundance, where a ~constant background
    and a wandering boundary have the *largest* fractional effect, so they
    are the most sensitive readout of whether peak detection / baseline does
    anything. They are *contaminated* (co-eluting isobars) so they do NOT
    gate; read them with ROBUST stats (median R^2, IQR) so per-peptide
    contaminants show up as outliers, not as aggregate poison.
  - f=0 baseline readout: at the 0% proportion the labelled contribution to
    every channel is ~0, so observed m4:mA / m5:mA there is background +
    tiny natural abundance. Correct baseline subtraction pulls it down and
    tightens it across peptides. Reported as median + IQR of the value at
    proportion 0 — a near-direct probe of baseline-subtraction quality.

Coverage: this metric deliberately does NOT require a peptide observed at all
9 proportions (the NB87a coverage gate). A line needs only >= --min-points
proportions, so the partial-coverage / low-abundance peptides — exactly the
population peak detection is meant to rescue — are kept. Curated/uncurated is
the m0:mA R^2>0.95 split, reported separately, not a filter.

Usage:
  python bench_mixing_linearity.py \
    --method fixed=<dir> [--method detected=<dir> ...] \
    --ground-truth <ground_truth.csv> \
    --output-dir <dir>

Outputs (under --output-dir):
  - mixing_linearity_peptides.csv   one row per (method, concat, quantity)
  - mixing_linearity_comparison.csv method x quantity x population summary
  - mixing_linearity_summary.json   full nested summary + f=0 readout
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from bench_aa_coefficients import ISO_COLS  # noqa: E402

GATE = ["m0_ma", "m1_ma", "m2_ma"]
PROBE = ["m4_ma", "m5_ma"]
QUANTITIES = ["m0_ma", "m1_ma", "m2_ma", "m3_ma", "m4_ma", "m5_ma"]
R2_GATE = 0.95


def parse_method(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise argparse.ArgumentTypeError(f"--method must be name=path, got {spec!r}")
    name, _, path = spec.partition("=")
    if not name or not path:
        raise argparse.ArgumentTypeError(f"--method must be name=path, got {spec!r}")
    return name, Path(path)


def load_fractions(
    inputs_dir: Path,
    ground_truth: pd.DataFrame,
    min_points: int,
) -> pd.DataFrame:
    """Load *_riana.txt, attach proportion, compute m_i:mA per row.

    Unlike bench_aa_coefficients.load_and_curate this keeps partial-coverage
    peptides (>= min_points proportions) — see module docstring.
    """
    gt = ground_truth.set_index("riana_filename")["nominal_proportion"]
    files = sorted(inputs_dir.glob("*_riana.txt"))
    if not files:
        raise FileNotFoundError(f"No *_riana.txt found in {inputs_dir}")

    frames = []
    for f in files:
        if f.name not in gt.index:
            raise KeyError(f"{f.name} not in ground_truth.csv")
        df = pd.read_csv(f, sep="\t", comment="#")
        # v0.9.0 files name the channels m0..m5; v1.0.0 already uses iso0..iso5.
        df = df.rename(columns={f"m{i}": f"iso{i}" for i in range(6)})
        if not set(ISO_COLS).issubset(df.columns):
            raise KeyError(f"{f.name} missing iso columns: have {list(df.columns)}")
        sub = df[["concat", *ISO_COLS]].copy()
        sub["proportion"] = float(gt[f.name])
        frames.append(sub)

    big = pd.concat(frames, ignore_index=True)
    big = big.drop_duplicates(subset=["concat", "proportion"])

    iso_sum = big[ISO_COLS].sum(axis=1).replace(0, np.nan)
    for i in range(6):
        big[f"m{i}_ma"] = big[f"iso{i}"] / iso_sum
    big = big.dropna(subset=["m0_ma"])

    counts = big.groupby("concat")["proportion"].nunique()
    keep = counts[counts >= min_points].index
    big = big[big["concat"].isin(keep)].copy()
    big["pep_len"] = big["concat"].str.rsplit("_", n=1).str[0].str.replace(
        r"\[.*?\]", "", regex=True).str.len()
    return big


def _linefit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Return (r2, resid_rmse, slope) of an ordinary least-squares line."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3 or np.var(x) == 0:
        return (np.nan, np.nan, np.nan)
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    resid = y - yhat
    resid_rmse = float(np.sqrt(np.mean(resid ** 2)))
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return (r2, resid_rmse, float(slope))


def per_peptide(big: pd.DataFrame) -> pd.DataFrame:
    """One row per (concat, quantity): line-fit r2/resid + value at f=0."""
    rows = []
    for concat, g in big.groupby("concat"):
        x = g["proportion"].to_numpy(dtype=float)
        at0 = g.loc[g["proportion"] == 0]
        pep_len = int(g["pep_len"].iloc[0])
        # m0:mA r2 defines the curated/uncurated split (the existing gate).
        m0_r2, _, _ = _linefit(x, g["m0_ma"].to_numpy(dtype=float))
        is_curated = bool(np.isfinite(m0_r2) and m0_r2 > R2_GATE)
        for q in QUANTITIES:
            r2, resid, slope = _linefit(x, g[q].to_numpy(dtype=float))
            rows.append({
                "concat": concat,
                "quantity": q,
                "pep_len": pep_len,
                "n_points": int(g["proportion"].nunique()),
                "r2": r2,
                "resid_rmse": resid,
                "slope": slope,
                "val_at_f0": float(at0[q].iloc[0]) if len(at0) else np.nan,
                "is_curated": is_curated,
            })
    return pd.DataFrame(rows)


def _agg(df: pd.DataFrame) -> dict:
    r2 = df["r2"].to_numpy(dtype=float)
    r2 = r2[np.isfinite(r2)]
    resid = df["resid_rmse"].to_numpy(dtype=float)
    resid = resid[np.isfinite(resid)]
    if r2.size == 0:
        return {"n_peptides": 0}
    return {
        "n_peptides": int(df["concat"].nunique()),
        "mean_r2": float(np.mean(r2)),
        "median_r2": float(np.median(r2)),
        "frac_r2_gt_gate": float(np.mean(r2 > R2_GATE)),
        "median_resid": float(np.median(resid)) if resid.size else float("nan"),
        "p95_resid": float(np.percentile(resid, 95)) if resid.size else float("nan"),
        "resid_iqr": float(np.percentile(resid, 75) - np.percentile(resid, 25))
        if resid.size else float("nan"),
    }


def summarize_method(pp: pd.DataFrame) -> dict:
    out: dict = {"quantities": {}}
    for q in QUANTITIES:
        qdf = pp[pp["quantity"] == q]
        out["quantities"][q] = {
            "all": _agg(qdf),
            "curated": _agg(qdf[qdf["is_curated"]]),
            "uncurated": _agg(qdf[~qdf["is_curated"]]),
        }
    # f=0 baseline readout on the probe channels (robust).
    out["f0_readout"] = {}
    for q in PROBE:
        v = pp.loc[pp["quantity"] == q, "val_at_f0"].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        out["f0_readout"][q] = {
            "n": int(v.size),
            "median": float(np.median(v)) if v.size else float("nan"),
            "iqr": float(np.percentile(v, 75) - np.percentile(v, 25))
            if v.size else float("nan"),
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--method", type=parse_method, action="append",
                        required=True, metavar="NAME=DIR",
                        help="named integration output dir; repeatable")
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-points", type=int, default=4,
                        help="min proportions a peptide must appear at (default 4)")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    methods = dict(args.method)
    if len(methods) != len(args.method):
        parser.error("duplicate --method name")
    ground_truth = pd.read_csv(args.ground_truth)

    all_pp = []
    summaries = {}
    comp_rows = []
    for name, inputs_dir in methods.items():
        print(f"[load] {name!r} <- {inputs_dir}", flush=True)
        big = load_fractions(inputs_dir, ground_truth, args.min_points)
        pp = per_peptide(big)
        pp.insert(0, "method", name)
        all_pp.append(pp)
        summ = summarize_method(pp)
        summaries[name] = summ
        n_pep = pp["concat"].nunique()
        n_cur = pp.loc[pp["quantity"] == "m0_ma", "is_curated"].sum()
        print(f"       {n_pep} peptides ({n_cur} curated by m0:mA R^2>{R2_GATE})",
              flush=True)
        for q in QUANTITIES:
            for pop in ("all", "curated", "uncurated"):
                a = summ["quantities"][q][pop]
                comp_rows.append({
                    "method": name, "quantity": q, "population": pop,
                    "n_peptides": a.get("n_peptides", 0),
                    "mean_r2": a.get("mean_r2", float("nan")),
                    "median_r2": a.get("median_r2", float("nan")),
                    "frac_r2_gt_gate": a.get("frac_r2_gt_gate", float("nan")),
                    "median_resid": a.get("median_resid", float("nan")),
                })

    pd.concat(all_pp, ignore_index=True).to_csv(
        args.output_dir / "mixing_linearity_peptides.csv", index=False)
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(args.output_dir / "mixing_linearity_comparison.csv", index=False)
    with (args.output_dir / "mixing_linearity_summary.json").open("w") as f:
        json.dump({"methods": {k: str(v) for k, v in methods.items()},
                   "min_points": args.min_points, "summaries": summaries}, f, indent=2)

    # ---- console report ----
    print("\n[GATE] median R^2 (uncurated population — where detection should win):")
    piv = comp[(comp["population"] == "uncurated") & (comp["quantity"].isin(GATE))]
    print(piv.pivot(index="quantity", columns="method",
                    values="median_r2").round(4).to_string())
    print("\n[GATE] fraction crossing R^2>%.2f (uncurated):" % R2_GATE)
    print(piv.pivot(index="quantity", columns="method",
                    values="frac_r2_gt_gate").round(4).to_string())
    print("\n[PROBE] median R^2 (robust; all population):")
    pivp = comp[(comp["population"] == "all") & (comp["quantity"].isin(PROBE))]
    print(pivp.pivot(index="quantity", columns="method",
                     values="median_r2").round(4).to_string())
    print("\n[PROBE] f=0 baseline readout (median | IQR; lower+tighter = cleaner):")
    for name in methods:
        rd = summaries[name]["f0_readout"]
        cells = "  ".join(
            f"{q}={rd[q]['median']:.4f}|{rd[q]['iqr']:.4f}" for q in PROBE)
        print(f"   {name:18s} {cells}")
    print(f"\n[done] outputs -> {args.output_dir}")


if __name__ == "__main__":
    main()
