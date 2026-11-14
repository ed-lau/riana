"""Design and gated-completeness audit of fit runs (wls-var / D / E scope, §2).

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md``. For each fit
directory (``riana_fit_fractions.txt`` + ``riana_fit_peptides.txt``):

- design: per (experiment, condition) the timepoints and replicates per timepoint;
- gated completeness: peptides passing an approximation of the rollup gate
  (R² ≥ 0.8, or R² ≥ 0.6 with k_cv < 0.2; Spep ≥ 8; admission ``all``). Per
  protein × condition with ≥ 2 admitted peptides, ``presence`` is the share of
  (peptide, t > 0 cell) pairs observed and ``complete`` the share of cells
  holding every admitted peptide. Identical mixing in both conditions (C, D)
  is only possible in complete cells; E does not need it.

python tests/benchmark/wlsvar_audit.py <fit_dir> [<fit_dir> ...]
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

CELL = ["labeling_time", "biological_replicate"]


def design(x: pd.DataFrame) -> None:
    for (e, c), g in x.groupby(["experiment", "condition"]):
        reps = g.groupby("labeling_time").biological_replicate.nunique()
        tps = sorted(g.labeling_time.unique().tolist())
        more = "…" if len(tps) > 12 else ""
        line = f"  {e} | {c}: {len(tps)} tp {tps[:12]}{more}; reps/tp {reps.min()}–{reps.max()}"
        if "individual" in g:
            line += f"; individuals {g.individual.nunique()}"
        print(line)


def admitted(pep: pd.DataFrame) -> set[tuple[str, str]]:
    ok = pep.R_squared >= 0.8
    if "k_cv" in pep:
        ok |= (pep.R_squared >= 0.6) & (pep.k_cv < 0.2)
    if "spep" in pep:
        ok &= pep.spep >= 8
    pep = pep[ok]
    ncond = pep.groupby("experiment").condition.nunique()
    if ncond.max() > 1:  # admission `all`: passes in every condition of its experiment
        cnt = pep.groupby(["experiment", "concat"]).condition.nunique()
        keep = cnt[cnt == cnt.index.get_level_values(0).map(ncond)].index
        pep = pep.set_index(["experiment", "concat"]).loc[keep].reset_index()
    return set(zip(pep.experiment, pep.concat))


def completeness(x: pd.DataFrame, adm: set[tuple[str, str]]) -> pd.DataFrame:
    x = x[x.labeling_time > 0]
    x = x[[k in adm for k in zip(x.experiment, x.concat)]]
    x = x[~x["protein id"].astype(str).str.contains(";")]
    rows = []
    for _, g in x.groupby(["experiment", "condition", "protein id"]):
        peps = g.concat.unique()
        if len(peps) < 2:
            rows.append(dict(n=len(peps)))
            continue
        cells = pd.MultiIndex.from_frame(g[CELL].drop_duplicates())
        m = g.pivot_table(index=CELL, columns="concat", values="fs",
                          aggfunc="size").reindex(cells).notna()
        rows.append(dict(n=len(peps), presence=m.values.mean(), complete=m.all(axis=1).mean()))
    return pd.DataFrame(rows)


def main(dirs: list[Path]) -> None:
    for d in dirs:
        f = d / "riana_fit_fractions.txt"
        hdr = pd.read_table(f, comment="#", nrows=0).columns
        cols = ["concat", "protein id", "fs", "experiment", "condition", *CELL, "individual"]
        x = pd.read_table(f, comment="#", usecols=[c for c in cols if c in hdr],
                          low_memory=False)
        x = x[x.fs.notna()]
        print(f"\n===== {d}  rows={len(x):,}")
        design(x)
        pep = pd.read_table(d / "riana_fit_peptides.txt", comment="#", low_memory=False)
        r = completeness(x, admitted(pep))
        multi = r[r.n >= 2]
        print(f"  protein x condition {len(r):,}; single-peptide {(r.n < 2).mean():.0%}; "
              f"multi-peptide: median n {multi.n.median():.0f}, presence "
              f"{multi.presence.median():.2f}, complete cells {multi.complete.median():.2f}, "
              f"fully complete {(multi.complete == 1).mean():.0%}")


if __name__ == "__main__":
    main([Path(a) for a in sys.argv[1:]])
