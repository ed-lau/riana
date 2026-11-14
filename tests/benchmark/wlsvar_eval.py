"""Truth-free yardsticks for the wls-var / D / E evaluation (WP2).

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md`` §3. Every fit directory is
re-tested on the SAME gated peptides under each collapse × weights with the shipped
functions, via ``eval_rising_variance`` (``capture_gated`` runs the shipped rollup and
records the frame it collapses; ``collapse(shared=True)`` is candidate C):

- A = per-condition collapse + ``wls`` (shipped default)
- B = per-condition collapse + ``wls-var``
- C = shared-proportion collapse + ``wls``  (also C + ``wls-var`` for completeness)

Tests:

- ``t3 ODD EVEN --ref R --test T``: time split. z = (Δk_odd − Δk_even) / √(SE²_odd + SE²_even)
  per protein; reports its median, SD, robust SD (1.4826 × MAD) and the share with |z| > 1.96,
  overall and by admitted peptides.
- ``null DIR [DIR ...] --ref R --test T``: replicate split (T1) or permutation (T2) refits;
  rejection rate of the Δk test at α 0.05 / 0.01 pooled over the refits, the SD of Δk/SE,
  and the BH < 0.05 calls per refit.
- ``inter DIR [DIR ...] --ref R --test T [-W N]``: the shipped ``--interaction`` model on the
  same refits (T2 permutations), under A and C (``wls``): rejection rate at α 0.05 of its
  mean-Δk, ratio-scale (primary) and additive tests, and their BH < 0.05 calls per refit.
- ``t6 FIT --ref LV --test RV [--split control ATF6-KO]``: animal-set split. The fit's
  experiment (chamber) and condition (genotype) swap roles: one directory per genotype,
  holding only the two chambers compared (admission ``all`` spans every condition of an
  experiment), then the T3 z between the two genotypes' chamber Δk. No refit: the curves
  are per (chamber, genotype) either way.

python tests/benchmark/wlsvar_eval.py t3 runs/wlsvar_eval/lve_atr/t3_odd runs/wlsvar_eval/lve_atr/t3_even --ref control --test atrium
python tests/benchmark/wlsvar_eval.py null runs/wlsvar_eval/boomi/t1/t1_* --ref A --test B
python tests/benchmark/wlsvar_eval.py t6 runs/wlsvar_eval/atf6/fit --ref LV --test RV
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.benchmark.eval_rising_variance import (  # noqa: E402
    capture_gated, collapse, dk_shipped, interaction_shipped)

VARIANTS = [("A", False, "wls"), ("B", False, "wls-var"), ("C", True, "wls"),
            ("C+wls-var", True, "wls-var")]
KEY = ["experiment", "protein"]


def dk_tables(fit_dir: Path, ref: str, test: str) -> dict[str, pd.DataFrame]:
    """Δk per (experiment, protein) under each variant, with the admitted peptide count."""
    seen = capture_gated(fit_dir, ref, test)
    f = seen["frame"]
    npep = f.groupby(KEY).concat.nunique().rename("n_pep")
    out = {}
    for name, shared, w in VARIANTS:
        t = dk_shipped(collapse(seen, shared), w, ref, test)
        out[name] = t.merge(npep.reset_index(), on=KEY, how="left")
    return out


def _rsd(v: pd.Series) -> float:
    return float(1.4826 * (v - v.median()).abs().median())


def t3(odd: Path, even: Path, ref: str, test: str) -> pd.DataFrame:
    a, b = dk_tables(odd, ref, test), dk_tables(even, ref, test)
    rows = []
    for name, *_ in VARIANTS:
        m = a[name].merge(b[name], on=KEY, suffixes=("_o", "_e")).dropna(
            subset=["dk_o", "dk_e", "se_o", "se_e"])
        z = (m.dk_o - m.dk_e) / np.sqrt(m.se_o ** 2 + m.se_e ** 2)
        npep = np.minimum(m.n_pep_o, m.n_pep_e)
        for lab, sel in [("all", npep > 0), ("1", npep == 1), ("2-3", npep.between(2, 3)),
                         ("4+", npep >= 4)]:
            zz = z[sel]
            if len(zz) < 5:
                continue
            rows.append(dict(variant=name, peptides=lab, n=len(zz), z_median=zz.median(),
                             z_sd=zz.std(), z_robust_sd=_rsd(zz),
                             share_abs_z_gt_196=(zz.abs() > 1.96).mean(),
                             dk_diff_robust_sd=_rsd((m.dk_o - m.dk_e)[sel]),
                             se_comb_median=np.sqrt(m.se_o ** 2 + m.se_e ** 2)[sel].median()))
    return pd.DataFrame(rows)


def null(dirs: list[Path], ref: str, test: str) -> pd.DataFrame:
    rows = []
    for d in dirs:
        for name, t in dk_tables(d, ref, test).items():
            t = t.dropna(subset=["p"])
            z = t.dk / t.se
            kq = pd.qcut((t.k_A + t.k_B) / 2, 4, labels=False)  # turnover quartile, per refit
            for lab, sel in [("all", t.n_pep > 0), ("1", t.n_pep == 1),
                             ("2+", t.n_pep >= 2)] + [(f"k Q{q + 1}", kq == q) for q in range(4)]:
                if not sel.any():
                    continue
                rows.append(dict(refit=d.name, variant=name, peptides=lab, tests=int(sel.sum()),
                                 rej05=(t.p[sel] < .05).mean(), rej01=(t.p[sel] < .01).mean(),
                                 z_sd=z[sel].std(), calls_bh05=int((t.p_adj[sel] < .05).sum())))
    r = pd.DataFrame(rows)
    pooled = r.groupby(["variant", "peptides"]).apply(lambda g: pd.Series(dict(
        refits=len(g), tests=g.tests.mean(),
        rej05=np.average(g.rej05, weights=g.tests), rej01=np.average(g.rej01, weights=g.tests),
        rej05_min=g.rej05.min(), rej05_max=g.rej05.max(),
        z_sd=g.z_sd.mean(), calls_bh05_mean=g.calls_bh05.mean())), include_groups=False)
    return pooled.reset_index()


INTER_TESTS = [("mean Δk", "mean_delta_k_p"), ("ratio (primary)", "interaction_p"),
               ("additive", "interaction_additive_p")]


def inter(dirs: list[Path], ref: str, test: str, workers: int) -> pd.DataFrame:
    rows = []
    for d in dirs:
        seen = capture_gated(d, ref, test)
        for name, shared in [("A", False), ("C", True)]:
            t = interaction_shipped(collapse(seen, shared), "wls", workers, ref, test)
            for lab, col in INTER_TESTS:
                p = t[col].dropna()
                rows.append(dict(refit=d.name, variant=name, test=lab, tests=len(p),
                                 rej05=float((p < .05).mean()),
                                 calls_bh05=int((t[col + "_adj"] < .05).sum())))
    r = pd.DataFrame(rows)
    pooled = r.groupby(["variant", "test"]).apply(lambda g: pd.Series(dict(
        refits=len(g), tests=g.tests.mean(), rej05=np.average(g.rej05, weights=g.tests),
        rej05_min=g.rej05.min(), rej05_max=g.rej05.max(),
        calls_bh05_mean=g.calls_bh05.mean())), include_groups=False)
    return pooled.reset_index()


FIT_FILES = ("riana_fit_peptides.txt", "riana_fit_fractions.txt")


def t6_dirs(fit_dir: Path, out: Path, chambers: tuple[str, str],
            split: tuple[str, str]) -> list[Path]:
    """One fit directory per animal set (genotype), experiment ← "heart", condition ← chamber."""
    dirs = [out / g for g in split]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
    for name in FIT_FILES:
        with open(fit_dir / name) as fh:
            head = []
            for line in fh:
                if not line.startswith("#"):
                    break
                head.append(line)
        df = pd.read_csv(fit_dir / name, sep="\t", skiprows=len(head), low_memory=False)
        for g, d in zip(split, dirs):
            s = df[(df.condition == g) & df.experiment.isin(chambers)].copy()
            s["condition"], s["experiment"] = s.experiment, "heart"
            with open(d / name, "w") as fh:
                fh.writelines(head)
                s.to_csv(fh, sep="\t", index=False)
    return dirs


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("test", choices=["t3", "null", "t6", "inter"])
    p.add_argument("-W", "--workers", type=int, default=8, help="inter: worker processes")
    p.add_argument("dirs", type=Path, nargs="+")
    p.add_argument("--ref", required=True)
    p.add_argument("--test-condition", "--test", dest="cond", required=True)
    p.add_argument("--split", nargs=2, default=["control", "ATF6-KO"],
                   help="t6: the two animal sets (conditions of the fit)")
    p.add_argument("--work", type=Path, help="t6: where to write the swapped fit dirs")
    p.add_argument("-o", "--out", type=Path, help="write the table here (TSV)")
    a = p.parse_args()
    if a.test == "t3":
        if len(a.dirs) != 2:
            raise SystemExit("t3 takes ODD and EVEN fit directories")
        r = t3(*a.dirs, a.ref, a.cond)
    elif a.test == "inter":
        r = inter(a.dirs, a.ref, a.cond, a.workers)
    elif a.test == "t6":
        work = a.work or a.dirs[0].parent / f"t6_{a.ref}_{a.cond}"
        r = t3(*t6_dirs(a.dirs[0], work, (a.ref, a.cond), tuple(a.split)), a.ref, a.cond)
    else:
        r = null(a.dirs, a.ref, a.cond)
    with pd.option_context("display.width", 200, "display.max_rows", 200,
                           "display.float_format", "{:.3f}".format):
        print(r.to_string(index=False))
    if a.out:
        r.to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
