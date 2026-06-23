"""run_calibration_benchmark.py — the standing calibration A/B harness.

A thin orchestrator (no new science) over the tools that already exist:
:func:`run_integrate_v1_0_0.run_line` (Percolator-path integrate) and
``bench_fs_method_compare`` (|θ−f| recovery vs the ground-truth mixing
proportion). Per line it integrates (or **reuses** an existing
``integrate_outputs/<label>/`` when complete), scores recovery, writes the
canonical standing layout, and appends a one-line summary to ``BASELINE.md`` so
any future integration-knob / MBR tuning is a one-command A/B against a recorded
baseline rather than an ad-hoc re-run. See
``reports/2026-06-23_calibration_benchmark_harness.md`` for the design.

Output layout (per line)::

    runs/calib_<line>/                      # gitignored, regenerable
      <label>/
        recovery/fs_method_compare_summary.json   # the recovery anchor
        recovery/fs_method_compare.csv
        config.json                               # exact integrate + scoring config
      BASELINE.md                                 # one frozen line per <label>

Integrate output itself stays in its established home
``tests/data/calibration_d2o_mixing/<line>/integrate_outputs/<label>/`` (where
every other bench reads it); ``config.json`` records the path.

The default integrate config + ``--fs 0 1 2 3`` reproduces the committed
``benchmark_results/<line>/v1.0.0_fs0123/`` anchor (within±0.05: ac16 24.8 %,
cm 21.5 %, ipsc 28.9 %).

Examples::

    # Reproduce / refresh the fs0123 anchor for all three lines (reuses integrate):
    python tests/benchmark/run_calibration_benchmark.py --label v1.0.0 --fs 0 1 2 3

    # The proposed production-defaults baseline (re-integrates):
    python tests/benchmark/run_calibration_benchmark.py --label baseline \\
        --peak-rt apex --mass-tol 10 --integration-half-width 0.15 --fs 0 1 2 3

    # An adaptive-N_ISO A/B on one line:
    python tests/benchmark/run_calibration_benchmark.py --line ac16 \\
        --label adaptive --adaptive --fs 0 1 2 3

Within-protein-θ scoring (Track D) is intentionally *not* wired here: the mixing
series has ground-truth θ, so recovery is the right metric; the internal-
consistency θ-spread metric is for the no-ground-truth in-vivo series (use
``bench_within_protein_theta.py`` directly there).
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCH_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH_DIR))
sys.path.insert(0, str(REPO_ROOT))

from run_integrate_v1_0_0 import run_line  # noqa: E402

CALIB_ROOT = REPO_ROOT / "tests" / "data" / "calibration_d2o_mixing"
RUNS_ROOT = REPO_ROOT / "runs"
LINES = ("ac16", "cm", "ipsc")

# Per-line recovery config: the coefficient table and the proportions dropped
# from scoring. cm's 50 % fraction is a known outlier in its coefficient fit, so
# the recorded baseline trains on the drop50 table and drops proportion 50 (see
# the harness report + 2026-06-23_adaptive_niso_limited_isotopomer.md).
LINE_RECOVERY = {
    "ac16": {"coeff": "d2o_aa_coefficients_ac16.csv", "drop": []},
    "ipsc": {"coeff": "d2o_aa_coefficients_ipsc.csv", "drop": []},
    "cm": {"coeff": "d2o_aa_coefficients_cm_drop50.csv", "drop": [50.0]},
}


def _n_fractions(line: str) -> int:
    gt = CALIB_ROOT / line / "ground_truth.csv"
    with gt.open() as f:
        return sum(1 for _ in csv.DictReader(f))


def _integrate_complete(integrate_dir: Path, n_fractions: int) -> bool:
    if not integrate_dir.is_dir():
        return False
    return len(list(integrate_dir.glob("*_riana.txt"))) >= n_fractions


def _integrate_overrides(args: argparse.Namespace) -> dict:
    ihw = ("auto" if args.integration_half_width == "auto"
           else float(args.integration_half_width))
    # Mirror run_integrate_v1_0_0: apex/auto picking needs a wider extraction
    # window than the integration window so an off-MS2 apex still gets its full
    # width; ms2/fixed extracts exactly the window.
    if args.peak_rt == "ms2" and ihw != "auto":
        extraction_half_width = float(ihw)
    else:
        w = 0.33 if ihw == "auto" else float(ihw)
        extraction_half_width = w + 0.33
    return {
        "peak_rt": args.peak_rt,
        "mass_tol_ppm": args.mass_tol,
        "integration_half_width": ihw,
        "baseline_method": args.baseline_method,
        "extraction_half_width": extraction_half_width,
        "adaptive_iso": args.adaptive,
        "ria_max": args.ria,
    }


def _score_channels(fs: list[int]) -> int | None:
    """Translate an --fs channel list into bench_fs_method_compare's leading-N
    count. The bench scores the *leading* N channels only, so --fs must be a
    contiguous 0..N-1 prefix; an empty list means score all populated channels."""
    if not fs:
        return None
    if fs != list(range(len(fs))):
        raise SystemExit(
            f"--fs must be a contiguous leading prefix (0..N-1); got {fs}. "
            "The recovery bench scores the leading N channels only."
        )
    return len(fs)


def _run_recovery(line: str, integrate_dir: Path, out_dir: Path,
                  score_channels: int | None, ria: float) -> dict:
    rec = LINE_RECOVERY[line]
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "tests.benchmark.bench_fs_method_compare",
        "--inputs", str(integrate_dir),
        "--ground-truth", str(CALIB_ROOT / line / "ground_truth.csv"),
        "--coefficients", str(CALIB_ROOT / line / rec["coeff"]),
        "--output-dir", str(out_dir),
        "--ria", str(ria),
    ]
    if score_channels is not None:
        cmd += ["--score-channels", str(score_channels)]
    if rec["drop"]:
        cmd += ["--drop-proportion", *[str(p) for p in rec["drop"]]]
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)
    with (out_dir / "fs_method_compare_summary.json").open() as f:
        return json.load(f)


def _append_baseline(line: str, label: str, args: argparse.Namespace,
                     summary: dict) -> None:
    new = summary["new_overall"]
    md_path = RUNS_ROOT / f"calib_{line}" / "BASELINE.md"
    md_path.parent.mkdir(parents=True, exist_ok=True)
    if not md_path.exists():
        md_path.write_text(
            f"# Calibration recovery baselines — {line}\n\n"
            "Standing `bench_fs_method_compare` anchors (new-method FS, |θ−f| vs "
            "the ground-truth mixing proportion). One line per config `<label>`; "
            "re-running a label appends a fresh line (newest last).\n\n"
        )
    fs = "all" if not args.fs else "".join(str(c) for c in args.fs)
    with md_path.open("a") as f:
        f.write(
            f"- `{label}` (fs={fs}, ria={args.ria}, peak={args.peak_rt}, "
            f"mt={args.mass_tol}, ihw={args.integration_half_width}"
            f"{', adaptive' if args.adaptive else ''}): "
            f"within±0.05 **{new['frac_within_0.05'] * 100:.1f}%** · "
            f"within±0.10 {new['frac_within_0.10'] * 100:.1f}% · "
            f"IQR {new['iqr']:.3f} · median {new['median']:+.3f} · "
            f"n={new['n']} ({summary['n_peptides']} pep)\n"
        )


def run_one(line: str, args: argparse.Namespace) -> dict:
    label = args.label
    integrate_dir = CALIB_ROOT / line / "integrate_outputs" / label
    overrides = _integrate_overrides(args)
    n_fractions = _n_fractions(line)

    if args.force or not _integrate_complete(integrate_dir, n_fractions):
        print(f"\n=== {line.upper()} [{label}]: integrating ===")
        run_line(line, label, overrides)
    else:
        print(f"\n=== {line.upper()} [{label}]: reusing integrate "
              f"({integrate_dir}) ===")

    analysis_home = RUNS_ROOT / f"calib_{line}" / label
    score_channels = _score_channels(args.fs)
    summary = _run_recovery(line, integrate_dir, analysis_home / "recovery",
                            score_channels, args.ria)

    config = {
        "line": line,
        "label": label,
        "integrate_dir": str(integrate_dir),
        "reused_integrate": not args.force
        and _integrate_complete(integrate_dir, n_fractions),
        "integrate_overrides": overrides,
        "fs": args.fs,
        "score_channels": score_channels,
        "ria": args.ria,
        "recovery": LINE_RECOVERY[line],
    }
    (analysis_home / "config.json").write_text(json.dumps(config, indent=2))
    _append_baseline(line, label, args, summary)

    new = summary["new_overall"]
    print(f"[{line} {label}] within±0.05 {new['frac_within_0.05'] * 100:.1f}% · "
          f"IQR {new['iqr']:.3f} · n={new['n']} -> {analysis_home}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--line", choices=["ac16", "cm", "ipsc", "all"],
                        default="all")
    parser.add_argument("--label", required=True,
                        help="config subdir name (e.g. v1.0.0, baseline, mt10, "
                        "adaptive) — both the integrate_outputs/<label> dir and "
                        "the runs/calib_<line>/<label> analysis home")
    # Integrate knobs (passthrough -> run_line overrides).
    parser.add_argument("--peak-rt", choices=["ms2", "apex"], default="ms2")
    parser.add_argument("--mass-tol", type=float, default=15.0,
                        help="mass tolerance in ppm [default 15, the v1.0.0 anchor]")
    parser.add_argument("--integration-half-width", default="0.33",
                        help="RT half-width in min, or 'auto' to detect boundaries")
    parser.add_argument("--baseline-method",
                        choices=["none", "noise_floor", "snip", "asls"],
                        default="none")
    parser.add_argument("--adaptive", action="store_true",
                        help="adaptive N_ISO (--iso auto): per-peptide envelope channels")
    parser.add_argument("--ria", type=float, default=0.0598,
                        help="precursor enrichment (RIA max) [default 0.0598]")
    # Fit / scoring knob.
    parser.add_argument("--fs", type=int, nargs="*", default=[0, 1, 2, 3],
                        metavar="N",
                        help="limited-isotopomer scoring channels (leading prefix, "
                        "e.g. 0 1 2 3). Bare '--fs' scores all populated channels. "
                        "[default 0 1 2 3]")
    parser.add_argument("--force", action="store_true",
                        help="re-integrate even when integrate_outputs/<label> is complete")
    args = parser.parse_args()

    targets = list(LINES) if args.line == "all" else [args.line]
    results = {}
    for line in targets:
        results[line] = run_one(line, args)

    print("\n=== summary (within±0.05, new-method FS) ===")
    for line in targets:
        new = results[line]["new_overall"]
        print(f"  {line:5s} {args.label:16s} {new['frac_within_0.05'] * 100:5.1f}%  "
              f"IQR {new['iqr']:.3f}  median {new['median']:+.3f}")


if __name__ == "__main__":
    main()
