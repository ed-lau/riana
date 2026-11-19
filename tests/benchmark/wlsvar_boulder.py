"""Four-tissue helpers for the wls-var / D / E evaluation (EDL / SOL beside ventricle and atrium).

Scope: ``reports/2026-11-14_wls_var_D_E_evaluation.md`` §6. The combined ATR / LVE / EDL / SOL
search (``data/timeseries_boulder_4tissues``): one mouse per labeling day, all four tissues from
that mouse (``characteristics[individual]`` = ``mouse_Ddd``), the tissues as the conditions of one
experiment.

``pairs FIT OUT [--no-ids]`` — one fit directory per tissue pair (``OUT/<A>_<B>``) holding the
fit's rows for those two tissues only, so peptide admission (``all``) spans the pair as it
would in a two-tissue project. Conditions are renamed LVE / ATR / EDL / SOL. ``--no-ids``
blanks ``individual``, so E's animal becomes the tissue sample (as on lve_atr, which has no IDs).

``tissues FIT`` — the shared day deviation across tissues (WP1d for one animal per day). Per
tissue (admission ``own``), each protein's shipped ``wls`` line; each (protein, tissue, day)
cell's deviation on two scales: d = −(φ − φ̂)/t, an absolute shift in k, and r = (φ − φ̂)/φ̂ =
d/k̂, a relative one. On the cells held by all four tissues, a one-factor model per scale gives
each tissue's loading on the shared deviation, λᵢ² = mean over j, l ≠ i of cov(i, j)·cov(i, l) /
cov(j, l). A proportional deviation (δk ∝ k) has equal loadings on the r scale and loadings ∝ k
on the d scale; an additive one (δk the same in every tissue) the reverse. The same for the
proteome-wide day shift (median r over proteins per tissue × day).

python tests/benchmark/wlsvar_boulder.py pairs runs/deltak_eval/wlsvar_eval/boulder4/fit runs/deltak_eval/wlsvar_eval/boulder4/pairs
python tests/benchmark/wlsvar_boulder.py tissues runs/deltak_eval/wlsvar_eval/boulder4/fit
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from riana.core import protein  # noqa: E402
from tests.benchmark.eval_rising_variance import collapse  # noqa: E402
from tests.benchmark.wlsvar_eval import FIT_FILES  # noqa: E402
from tests.benchmark.wlsvar_wp1 import residuals  # noqa: E402

SHORT = {"left ventricle": "LVE", "left atrium": "ATR",
         "extensor digitorum longus": "EDL", "soleus": "SOL"}
ORDER = ["LVE", "ATR", "EDL", "SOL"]


def _read(path: Path) -> tuple[list[str], pd.DataFrame]:
    with open(path) as fh:
        head = []
        for line in fh:
            if not line.startswith("#"):
                break
            head.append(line)
    return head, pd.read_csv(path, sep="\t", skiprows=len(head), low_memory=False)


def pairs(fit_dir: Path, out: Path, no_ids: bool) -> list[Path]:
    tables = {name: _read(fit_dir / name) for name in FIT_FILES}
    dirs = []
    for a, b in combinations(ORDER, 2):
        d = out / f"{a}_{b}"
        d.mkdir(parents=True, exist_ok=True)
        for name, (head, df) in tables.items():
            s = df.assign(condition=df.condition.map(SHORT).fillna(df.condition))
            s = s[s.condition.isin([a, b])].copy()
            if no_ids and "individual" in s:
                s["individual"] = np.nan
            with open(d / name, "w") as fh:
                fh.writelines(head)
                s.to_csv(fh, sep="\t", index=False)
        dirs.append(d)
    return dirs


def capture_own(fit_dir: Path) -> dict:
    """The shipped linear rollup's gated frame, each tissue with its own peptides."""
    from typer.testing import CliRunner
    from riana.cli import app

    seen: dict = {}
    orig = protein._collapse_long

    def spy(fractions, *, method, min_peptides):
        seen.update(frame=fractions.copy(), method=method, min_peptides=min_peptides)
        return orig(fractions, method=method, min_peptides=min_peptides)

    protein._collapse_long = spy
    try:
        with tempfile.TemporaryDirectory() as tmp:
            res = CliRunner().invoke(app, ["rollup", str(fit_dir), "--model", "linear",
                                           "--peptide-admission", "own", "-o", tmp])
            if res.exit_code != 0:
                raise RuntimeError(f"rollup failed on {fit_dir}:\n{res.output[-2000:]}")
    finally:
        protein._collapse_long = orig
    return seen


def _loadings(w: pd.DataFrame) -> pd.Series:
    c = w.cov()
    out = {}
    for i in w.columns:
        rest = [j for j in w.columns if j != i]
        v = [c.loc[i, j] * c.loc[i, l] / c.loc[j, l] for j, l in combinations(rest, 2)]
        m = float(np.mean(v))
        out[i] = np.sqrt(m) if m > 0 else np.nan
    return pd.Series(out)


def tissues(fit_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    seen = capture_own(fit_dir)
    r = residuals(collapse(seen, False))
    r["tissue"] = r.condition.map(SHORT).fillna(r.condition)
    r["k"] = -r.fit / r.labeling_time
    r["d"] = -(r.phi - r.fit) / r.labeling_time
    r["r"] = (r.phi - r.fit) / r.fit
    lo, hi = r.r.quantile([0.01, 0.99])
    r = r[(r.r >= lo) & (r.r <= hi)]
    tis = [t for t in ORDER if t in set(r.tissue)]
    idx = ["protein", "labeling_time"]
    wide = {v: r.pivot_table(index=idx, columns="tissue", values=v)[tis].dropna()
            for v in ("d", "r", "k")}
    common = wide["d"].index.intersection(wide["r"].index)
    k_med = wide["k"].loc[wide["k"].index.intersection(common)].median()
    # per pair: correlation of the cell deviations
    prow = []
    for a, b in combinations(tis, 2):
        w = r.pivot_table(index=idx, columns="tissue", values="r")[[a, b]].dropna()
        prow.append(dict(pair=f"{a}–{b}", cells=len(w), corr_r=float(w[a].corr(w[b])),
                         k_ratio=float(k_med[b] / k_med[a])))
    # one-factor loadings on the cells every tissue holds, per scale
    lrow = pd.DataFrame({"k_median": k_med,
                         "load_d": _loadings(wide["d"].loc[common]),
                         "load_r": _loadings(wide["r"].loc[common]),
                         "sd_d": wide["d"].loc[common].std(),
                         "sd_r": wide["r"].loc[common].std()})
    lrow["load_d_over_k"] = lrow.load_d / lrow.k_median
    lrow.index.name = "tissue"
    lrow.attrs["cells"] = len(common)
    # proteome-wide day shift
    shift = r.groupby(["labeling_time", "tissue"]).r.median().unstack()[tis]
    srow = pd.DataFrame({"shift_sd_r": shift.std(), "load_shift_r": _loadings(shift)})
    srow["shift_corr_mean"] = [shift.corr().loc[t].drop(t).mean() for t in tis]
    srow.index.name = "tissue"
    return pd.DataFrame(prow), lrow.reset_index(), srow.reset_index()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["pairs", "tissues"])
    p.add_argument("fit_dir", type=Path)
    p.add_argument("out", type=Path, nargs="?")
    p.add_argument("--no-ids", action="store_true")
    a = p.parse_args()
    if a.cmd == "pairs":
        for d in pairs(a.fit_dir, a.out, a.no_ids):
            print(d)
        return
    pr, ld, sh = tissues(a.fit_dir)
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(pr.to_string(index=False))
        print(f"\none-factor loadings on {ld.attrs.get('cells', '?')} cells held by all tissues")
        print(ld.to_string(index=False))
        print("\nproteome-wide day shift (median r over proteins)")
        print(sh.to_string(index=False))
    if a.out:
        a.out.mkdir(parents=True, exist_ok=True)
        pr.to_csv(a.out / "tissues_pairs.tsv", sep="\t", index=False)
        ld.to_csv(a.out / "tissues_loadings.tsv", sep="\t", index=False)
        sh.to_csv(a.out / "tissues_shift.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
