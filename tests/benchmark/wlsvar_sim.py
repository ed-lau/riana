"""Calibrated simulation for the wls-var / D / E evaluation (WP3, then WP4's scoring).

Scope: ``reports/2026-11-13_wls_var_and_collapse_plan.md`` WP3. The manual's Monte Carlo
(``bench_rising_null.py``) stays as it is, so its numbers remain reproducible; this bench
adds the design and the nuisance that WP1–WP2 measured on real data, and runs the same
envelope → ``riana fit`` → shipped-rollup chain.

Design ``atf6`` (the shape of the ATF6 study): two conditions × four tissues sampled from
every animal, three animals per condition × day, days 1–14 (no day 0), RIA 0.046,
instant precursor (θ route, fitted as ATF6 was: depth 3, ``deberneh_2025_rss``, Spep ≥ 8).

Nuisance, each a ``--set`` knob (defaults are the WP3 calibration; see ``PARAMS``):

- protein k: log-normal (``k_med``, ``k_sd``); tissue factors from ATF6's chamber k ratios;
- peptides per protein: 2 + a geometric count (mean ``pep_extra``), at most ``pep_max``;
- envelope noise, relative + absolute, per-peptide multiplier (``noise_rel``, ``noise_abs``,
  ``noise_spread``); labeling-site error (``site_err``); interference (``interfere``);
- **per-peptide apparent-k offset** (``pep_k_sd``, log scale), shared by conditions, animals
  and tissues: the systematic peptide spread of WP1b (τ_rel ≈ 0.17 on ATF6);
- **per-animal proteome-wide turnover shift** (``shift_sd``, log scale), shared by the
  animal's tissues: WP1d's ≈ 8 % factor;
- per-animal × protein turnover (``animal_k_sd``), shared by the animal's tissues, and per
  animal × protein × tissue (``animal_tissue_k_sd``); a per-animal × tissue proteome-wide
  shift (``tissue_shift_sd``);
- per-animal enrichment (``ria_cv``);
- **missingness**: each peptide is missing from a run with probability
  ``miss`` × scale^(−``miss_slope``) (low-intensity peptides drop more);
- **timepoint-shared deviation** (off by default): per protein × tissue × day, shared by every
  animal sampled that day in both conditions, a turnover factor (``time_k_sd``, log scale), an
  absolute turnover offset (``time_kabs_sd``, per day; ATF6's σ_time is ≈ √(0.0072² +
  (0.088·k)²)) and an additive θ offset (``time_theta_sd``). It cancels from Δk and is what E+time models;
  ATF6 and boomi both show it. ``time_shared`` is the share of its variance common to the
  tissues of a day (ATF6 ≈ 0.85; 0 = independent per tissue). Drawn from a separate stream,
  so datasets made without it are unchanged;
- **per-peptide constant θ offset** (``pep_theta_sd``, off by default): shared by every run (a
  reference-envelope error). A through-origin line cannot absorb it, so it shows up as lack of
  fit shared by the conditions at each time, largest for slow proteins.

``truth.tsv`` holds k0 × tissue × ratio and ``M``, the protein's estimand factor: the mean over
its peptides of e^koff × Spep_true / Spep_tab. A protein-level estimate targets the
peptide average, so its true Δk is (ratio − 1) × k_true,A × M (``pooled`` reports coverage
against both).

The true ratio k_B / k_A (``--ratio``) applies to every protein, tissue and peptide, so
under ``--ratio 1`` every Δk test and every interaction test is a false positive, and under
``--ratio 1.5`` the ratio-scale interaction stays null.

Each dataset is scored under the variants in ``VARIANTS`` (A, B, C, C + ``wls-var``; D and E
join in WP4) and the shipped interaction under A and C. ``--calib`` also prints the WP3
calibration statistics (``wlsvar_calib.py``) for the dataset.

    python tests/benchmark/wlsvar_sim.py --work runs/wlsvar_calib/try1 --proteins 300 --calib
    python tests/benchmark/wlsvar_sim.py --work runs/wlsvar_sim/null --datasets 10
    python tests/benchmark/wlsvar_sim.py --work runs/wlsvar_sim/r15 --ratio 1.5 --datasets 4
    python tests/benchmark/wlsvar_sim.py --work runs/wlsvar_sim/null --score-only
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.benchmark.bench_rising_null import _random_peptide  # noqa: E402

COEFFS = "deberneh_2025_rss"
PARAMS = dict(                     # WP3 calibration to ATF6 (runs/wlsvar_calib/try6)
    days=[1.0, 3.0, 7.0, 10.0, 14.0],
    tissues={"LV": 1.0, "RV": 1.0, "LA": 1.17, "RA": 1.15},
    animals=3, ria=0.046, kp=100.0,
    k_med=0.08, k_sd=0.55,
    pep_extra=4.0, pep_max=30,
    noise_rel=0.085, noise_abs=0.006, noise_spread=0.5,
    site_err=0.10, interfere=0.08,
    pep_k_sd=0.17, shift_sd=0.085, tissue_shift_sd=0.02, animal_k_sd=0.06,
    animal_tissue_k_sd=0.05, ria_cv=0.02,
    miss=0.08, miss_slope=1.0,
    time_k_sd=0.0, time_kabs_sd=0.0, time_theta_sd=0.0, pep_theta_sd=0.0, time_shared=0.0,
    depth=3,
)
COND = ("A", "B")


# --- simulation --------------------------------------------------------------------------------

def _proteins(rng: np.random.Generator, n_proteins: int, P: dict) -> list[dict]:
    """The proteins and their peptides: the first draws of a dataset's stream."""
    from riana.algorithms.isotope_dist import spep_from_coefficients
    from riana.core.fitting import load_aa_coefficients

    coeffs = load_aa_coefficients(COEFFS)
    proteins = []
    for i in range(n_proteins):
        k0 = float(np.clip(P["k_med"] * math.exp(rng.normal(0, P["k_sd"])), 0.003, 1.5))
        n_pep = int(min(P["pep_max"], 2 + rng.geometric(1 / (1 + P["pep_extra"])) - 1))
        peps = []
        while len(peps) < n_pep:
            seq = _random_peptide(rng)
            spep_tab = spep_from_coefficients(seq, coeffs)
            if round(spep_tab) < 9:
                continue
            scale = float(math.exp(rng.normal(0, 1.0)))
            peps.append(dict(
                seq=seq, spep_tab=float(spep_tab),
                spep=max(1, int(round(spep_tab * math.exp(rng.normal(0, P["site_err"]))))),
                noise=float(math.exp(rng.normal(0, P["noise_spread"])
                                     - P["noise_spread"] ** 2 / 2)),
                scale=scale, koff=float(rng.normal(0, P["pep_k_sd"])),
                p_miss=float(min(0.95, P["miss"] * scale ** (-P["miss_slope"])))))
        M = float(np.mean([math.exp(q["koff"]) * q["spep"] / q["spep_tab"] for q in peps]))
        proteins.append(dict(pid=f"sp|SIM{i:05d}|SIM_{i}", k0=k0, peps=peps, M=M))
    return proteins


def estimand_factors(seed: int, n_proteins: int, P: dict) -> pd.DataFrame:
    """``M`` per protein, regenerated from a dataset's seed (datasets made before
    ``truth.tsv`` carried it)."""
    pr = _proteins(np.random.default_rng(seed), n_proteins, P)
    return pd.DataFrame({"protein": [p["pid"].split("|")[1] for p in pr],
                         "M": [p["M"] for p in pr]})


def simulate(work: Path, seed: int, n_proteins: int, ratio: float, P: dict,
             workers: int) -> Path:
    from riana.io.manifest import ManifestRow, append_manifest
    from riana.records import RunIdentity

    rng = np.random.default_rng(seed)
    proteins = _proteins(rng, n_proteins, P)
    animals = []
    for cond in COND:
        for day in P["days"]:
            for j in range(P["animals"]):
                shift = rng.normal(0, P["shift_sd"])
                u = rng.normal(0, P["animal_k_sd"], n_proteins)
                v = {t: rng.normal(0, P["animal_tissue_k_sd"], n_proteins)
                     + rng.normal(0, P["tissue_shift_sd"]) for t in P["tissues"]}
                k = {t: [p["k0"] * tf * (ratio if cond == "B" else 1.0)
                         * math.exp(shift + u[i] + v[t][i]) for i, p in enumerate(proteins)]
                     for t, tf in P["tissues"].items()}
                animals.append(dict(cond=cond, day=day, rep=j + 1,
                                    individual=f"{cond}_D{int(day)}_a{j + 1}",
                                    ria=P["ria"] * (1 + rng.normal(0, P["ria_cv"])),
                                    k=k, seed=int(rng.integers(2**31))))
    rng_t = np.random.default_rng(seed + 7919)       # its own stream: see the docstring
    tk = {t: {d: rng_t.normal(0, P["time_k_sd"], n_proteins) for d in P["days"]}
          for t in P["tissues"]}
    tth = {t: {d: rng_t.normal(0, P["time_theta_sd"], n_proteins) for d in P["days"]}
           for t in P["tissues"]}
    for pr in proteins:
        for pep in pr["peps"]:
            pep["toff"] = float(rng_t.normal(0, P["pep_theta_sd"]))
    tka = {t: {d: rng_t.normal(0, P["time_kabs_sd"], n_proteins) for d in P["days"]}
           for t in P["tissues"]}
    f = float(P.get("time_shared", 0.0))
    if f > 0:                      # drawn last, so datasets made without it are unchanged
        for comp, sd in ((tk, "time_k_sd"), (tth, "time_theta_sd"), (tka, "time_kabs_sd")):
            for d in P["days"]:
                s = rng_t.normal(0, P[sd], n_proteins)
                for t in P["tissues"]:
                    comp[t][d] = math.sqrt(1 - f) * comp[t][d] + math.sqrt(f) * s
    for a in animals:
        a["k"] = {t: [max(1e-4, kv * math.exp(tk[t][a["day"]][i]) + tka[t][a["day"]][i])
                      for i, kv in enumerate(ks)]
                  for t, ks in a["k"].items()}
        a["tadd"] = {t: tth[t][a["day"]] for t in P["tissues"]}
    work.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(workers) as ex:
        written = list(ex.map(_write_animal, [(work, a, proteins, P) for a in animals]))
    pd.DataFrame([dict(experiment=t, protein=p["pid"].split("|")[1], condition=c,
                       k_true=p["k0"] * tf * (ratio if c == "B" else 1.0), M=p["M"])
                  for p in proteins for t, tf in P["tissues"].items() for c in COND]
                 ).to_csv(work / "truth.tsv", sep="\t", index=False)
    mf = work / "riana_manifest.tsv"
    if mf.exists():
        mf.unlink()
    append_manifest(mf, [ManifestRow("integrate", str(path), RunIdentity(
        experiment=tissue, sample=path.stem, data_file=path.stem,
        biological_replicate=a["rep"], individual=a["individual"], labeling_time=a["day"],
        labeling_time_unit="day", condition=a["cond"], precursor_enrichment=P["ria"]))
        for a, paths in zip(animals, written) for tissue, path in paths])
    return mf


def _write_animal(arg) -> list[tuple[str, Path]]:
    """The animal's runs, one per tissue. Instant precursor: the envelope at fraction new f
    is (1 − f)·natural + f·plateau(RIA of this animal)."""
    from riana.algorithms.mass_calc import calculate_ion_mz
    from riana.algorithms.precursor import RisingEnvelope

    work, a, proteins, P = arg
    rng = np.random.default_rng(a["seed"])
    env = {}
    for p in proteins:
        for pep in p["peps"]:
            if pep["seq"] not in env:
                e0 = RisingEnvelope(pep["seq"], calculate_ion_mz(pep["seq"]), pep["spep"],
                                    a["ria"], 0.0, P["kp"])
                env[pep["seq"]] = (e0.natural[:6], e0.plateau[:6])
    out = []
    for tissue in P["tissues"]:
        rows = []
        for i, p in enumerate(proteins):
            for pep in p["peps"]:
                miss = rng.random() < pep["p_miss"]
                f = (-math.expm1(-a["k"][tissue][i] * math.exp(pep["koff"]) * a["day"])
                     + float(a["tadd"][tissue][i]) + pep["toff"])
                nat, pla = env[pep["seq"]]
                e = (1 - f) * nat + f * pla
                s = pep["noise"]
                obs = (e * (1 + s * P["noise_rel"] * rng.normal(size=6))
                       + s * P["noise_abs"] * rng.normal(size=6))
                if rng.random() < P["interfere"]:
                    obs[int(rng.integers(3, 6))] *= 1 + rng.uniform(0.2, 1.0)
                if miss:          # draws above are kept so the stream does not depend on it
                    continue
                obs = np.clip(obs, 1e-6, None) * 1e6 * pep["scale"]
                rows.append({"file_idx": 0, "scan": 1000 + len(rows), "charge": 2,
                             "concat": f"{pep['seq']}_2", "sequence": pep["seq"],
                             "sample": "sim", "percolator q-value": 1e-4,
                             "protein id": p["pid"], **{f"iso{c}": obs[c] for c in range(6)}})
        path = work / f"{tissue}_{a['individual']}_riana.txt"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        out.append((tissue, path))
    return out


def fit(mf: Path, out: Path, P: dict, workers: int) -> None:
    exe = Path(sys.executable).with_name("riana")
    res = subprocess.run([str(exe), "fit", "--manifest", str(mf), "--coefficients", COEFFS,
                          "-d", str(P["depth"]), "-W", str(workers), "-o", str(out)],
                         capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"riana fit failed:\n{res.stdout[-2000:]}\n{res.stderr[-2000:]}")


# --- scoring -----------------------------------------------------------------------------------

def _variants(only: list[str] | None = None):
    """name → function(seen) → Δk table (experiment, protein, dk, se, p, p_adj, k_A, k_B)."""
    from tests.benchmark.eval_rising_variance import collapse, dk_shipped
    v = {
        "A": lambda s: dk_shipped(collapse(s, False), "wls"),
        "B": lambda s: dk_shipped(collapse(s, False), "wls-var"),
        "C": lambda s: dk_shipped(collapse(s, True), "wls"),
        "C+wls-var": lambda s: dk_shipped(collapse(s, True), "wls-var"),
    }
    if not only or set(only) - set(v):           # WP4 prototypes (wlsvar_de.py)
        from tests.benchmark.wlsvar_de import VARIANTS as de
        v.update({k: (lambda s, fn=fn: fn(s)) for k, fn in de.items()})
    return v


def score(fit_dir: Path, ratio: float, workers: int, only: list[str] | None = None,
          interaction: bool = True) -> list[dict]:
    from scipy import stats

    from tests.benchmark.eval_rising_variance import (capture_gated, collapse,
                                                      interaction_shipped)
    seen = capture_gated(fit_dir)
    f = seen["frame"]
    npep = f.groupby(["experiment", "protein"]).concat.nunique().rename("n_pep")
    rows = []
    tests_dir = fit_dir.parent / "tests"
    tests_dir.mkdir(exist_ok=True)
    for name, fn in _variants([v for v in only if not v.startswith("EJ")] if only
                              else None).items():
        if only and name not in only:
            continue
        t0 = time.time()
        t = fn(seen).dropna(subset=["p"]).merge(npep.reset_index(), on=["experiment", "protein"],
                                                how="left")
        dt = time.time() - t0
        cols = [c for c in ("experiment", "protein", "k_A", "k_B", "dk", "se", "df", "p", "p_adj",
                            "n_pep", "sigma_animal", "sigma_time") if c in t]
        t[cols].to_csv(tests_dir / f"{name}.tsv", sep="\t", index=False)
        df = t["df"] if "df" in t else pd.Series(np.inf, index=t.index)
        q = stats.t.ppf(0.975, df)
        z = t.dk / t.se
        lr = np.log2(t.k_B / t.k_A).replace([np.inf, -np.inf], np.nan)
        kq = pd.qcut((t.k_A + t.k_B) / 2, 4, labels=False)
        row = dict(test="dk", variant=name, tests=len(t), rej05=float((t.p < .05).mean()),
                   rej01=float((t.p < .01).mean()), calls_bh05=int((t.p_adj < .05).sum()),
                   z_sd=float(z.std()), z_median=float(z.median()),
                   log2_ratio_median=float(lr.median()), seconds=round(dt, 1),
                   **{f"rej05_kQ{i + 1}": float((t.p[kq == i] < .05).mean()) for i in range(4)})
        if ratio == 1:
            row["cover95"] = float((t.dk.abs() <= q * t.se).mean())
        rows.append(row)
    if interaction:
        for name, shared in (("A", False), ("C", True), ("EJ", None), ("EJ+time", None)):
            if only and name not in only:
                continue
            if name.startswith("EJ"):            # WP4: the peptide-level joint model
                from tests.benchmark.wlsvar_de import interaction_E
                it = interaction_E(seen, workers=workers, time=name == "EJ+time")
            else:
                it = interaction_shipped(collapse(seen, shared), "wls", workers)
            it.to_csv(tests_dir / f"inter_{name}.tsv", sep="\t", index=False)
            for lab, col in (("inter mean Δk", "mean_delta_k_p"), ("inter ratio", "interaction_p"),
                             ("inter additive", "interaction_additive_p")):
                p = it[col].dropna()
                rows.append(dict(test=lab, variant=name, tests=len(p),
                                 rej05=float((p < .05).mean()), rej01=float((p < .01).mean()),
                                 calls_bh05=int((it[col + "_adj"] < .05).sum())))
    return rows


# --- driver ------------------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--datasets", type=int, default=1)
    ap.add_argument("--first", type=int, default=0, help="index of the first dataset")
    ap.add_argument("--proteins", type=int, default=300)
    ap.add_argument("--ratio", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=20261114)
    ap.add_argument("-W", "--workers", type=int, default=16)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override a PARAMS entry (repeatable)")
    ap.add_argument("--calib", action="store_true", help="print the WP3 calibration statistics")
    ap.add_argument("--score-only", action="store_true",
                    help="re-score existing datasets under --work (no simulation or fit)")
    ap.add_argument("--variants", nargs="*", help="score only these variants")
    ap.add_argument("--no-interaction", action="store_true")
    ap.add_argument("--pooled", action="store_true",
                    help="only pool the per-test tables already under --work")
    ap.add_argument("--size-adjust", type=Path, metavar="NULL_WORK",
                    help="power of the datasets under --work at each variant's own p threshold "
                         "for a 5 %% pooled null on NULL_WORK (ratio 1)")
    a = ap.parse_args()
    if a.size_adjust:
        size_adjusted(a.size_adjust, a.work, a.ratio)
        return
    if a.pooled:
        pooled(a.work, a.ratio)
        return
    P = dict(PARAMS)
    for kv in a.set:
        k, v = kv.split("=", 1)
        if k not in P:
            raise SystemExit(f"unknown parameter {k}")
        P[k] = type(P[k])(v) if not isinstance(P[k], (list, dict)) else json.loads(v)
    a.work.mkdir(parents=True, exist_ok=True)
    git = lambda *c: subprocess.run(["git", *c], capture_output=True, text=True).stdout.strip()
    if not a.score_only:
        (a.work / f"ratio{a.ratio:g}_params.json").write_text(json.dumps(dict(
            P, proteins=a.proteins, seed=a.seed, argv=sys.argv[1:],
            commit=git("rev-parse", "--short", "HEAD"),
            code_dirty=bool(git("status", "--porcelain", "--", "riana", "tests/benchmark"))),
            indent=2))
    all_rows = []
    for d in range(a.first, a.first + a.datasets):
        ds = a.work / f"ratio{a.ratio:g}_{d}"
        fit_dir = ds / "theta"
        t0 = time.time()
        if not a.score_only:
            if not (fit_dir / "riana_fit_peptides.txt").exists():
                mf = simulate(ds, a.seed + d + int(1000 * a.ratio), a.proteins, a.ratio, P,
                              a.workers)
                t1 = time.time()
                fit(mf, fit_dir, P, a.workers)
                print(f"dataset {d}: simulate {(t1 - t0) / 60:.1f} min, fit "
                      f"{(time.time() - t1) / 60:.1f} min", flush=True)
            if a.calib:
                from tests.benchmark.wlsvar_calib import stats
                s = stats(fit_dir, "A", "B")
                pd.Series(s).to_csv(ds / "calib.tsv", sep="\t", header=False)
                print(pd.Series(s).to_string(float_format="{:.4g}".format), flush=True)
        rows = score(fit_dir, a.ratio, a.workers, a.variants, not a.no_interaction)
        for r in rows:
            r.update(dataset=d, ratio=a.ratio)
        all_rows += rows
        tag = "_".join(a.variants) if a.variants else "all"
        pd.DataFrame(rows).to_csv(ds / f"scores_{tag}.tsv", sep="\t", index=False)
        print(f"dataset {d} scored ({(time.time() - t0) / 60:.1f} min)", flush=True)
    R = pd.DataFrame(all_rows)
    tag = "_".join(a.variants) if a.variants else "all"
    R.to_csv(a.work / f"ratio{a.ratio:g}_scores_{tag}_{a.first}-{a.first + a.datasets - 1}.tsv",
             sep="\t", index=False)
    summarise(R)


def pooled(work: Path, ratio: float) -> pd.DataFrame:
    """Pool the per-test tables of every dataset under ``work`` (``ratio<r>_<d>/tests/``).

    Under the null, the per-dataset z SD leaves out what the datasets do not share (each has
    its own animals, so its own proteome-wide shift); the pooled z SD and rejection rate are
    the marginal calibration a test must hold."""
    from scipy import stats
    rows = []
    prm = work / f"ratio{ratio:g}_params.json"
    meta = json.loads(prm.read_text()) if prm.exists() else None
    for ds in sorted(work.glob(f"ratio{ratio:g}_*")):
        if not (ds / "tests").is_dir():
            continue
        truth = None
        if (ds / "truth.tsv").exists():
            tr = pd.read_table(ds / "truth.tsv")
            if "M" not in tr and meta:             # made before truth.tsv carried M
                d = int(ds.name.rsplit("_", 1)[1])
                P = {k: meta.get(k, v) for k, v in PARAMS.items()}
                tr = tr.merge(estimand_factors(meta["seed"] + d + int(1000 * ratio),
                                               meta["proteins"], P), on="protein")
            if "M" in tr:
                w = tr.pivot_table(index=["experiment", "protein"], columns="condition",
                                   values="k_true")
                truth = pd.DataFrame({"dk_true": w["B"] - w["A"]}).join(
                    tr.groupby(["experiment", "protein"]).M.first()).reset_index()
        for f in sorted((ds / "tests").glob("*.tsv")):
            if f.stem.startswith("inter_"):
                continue
            t = pd.read_table(f).dropna(subset=["p"])
            if truth is not None:
                t = t.merge(truth, on=["experiment", "protein"], how="left")
            t["variant"], t["dataset"] = f.stem, ds.name
            rows.append(t)
    T = pd.concat(rows, ignore_index=True)
    T["z"] = T.dk / T.se
    T["df"] = T["df"].fillna(np.inf) if "df" in T else np.inf
    T["kq"] = T.groupby(["variant", "dataset"], group_keys=False).apply(
        lambda d: pd.qcut((d.k_A + d.k_B) / 2, 4, labels=False), include_groups=False)
    out = []
    for v, d in T.groupby("variant", sort=False):
        per = d.groupby("dataset")
        r = dict(variant=v, datasets=per.ngroups, tests=len(d) / per.ngroups,
                 rej05=float((d.p < .05).mean()), rej01=float((d.p < .01).mean()),
                 rej05_min=float(per.p.apply(lambda x: (x < .05).mean()).min()),
                 rej05_max=float(per.p.apply(lambda x: (x < .05).mean()).max()),
                 calls_bh05=float(per.p_adj.apply(lambda x: (x < .05).sum()).mean()),
                 z_sd_pooled=float(d.z.std()), z_sd_within=float(per.z.std().mean()),
                 log2_ratio=float(np.log2(d.k_B / d.k_A).replace([np.inf, -np.inf],
                                                                  np.nan).median()),
                 **{f"rej05_kQ{i + 1}": float((d.p[d.kq == i] < .05).mean()) for i in range(4)})
        q = stats.t.ppf(0.975, d["df"])
        if ratio == 1:
            r["cover95"] = float((d.dk.abs() <= q * d.se).mean())
        elif "dk_true" in d:
            r["cover95_k0"] = float(((d.dk - d.dk_true).abs() <= q * d.se).mean())
            r["cover95_estimand"] = float(((d.dk - d.dk_true * d.M).abs() <= q * d.se).mean())
        out.append(r)
    out = pd.DataFrame(out).set_index("variant")
    with pd.option_context("display.width", 250, "display.max_columns", 30,
                           "display.float_format", "{:.3f}".format):
        print(out.to_string())
    out.to_csv(work / f"ratio{ratio:g}_pooled.tsv", sep="\t")
    return out


INTER_TESTS = {"mean Δk": "mean_delta_k_p", "ratio": "interaction_p",
               "additive": "interaction_additive_p"}


def _tables(work: Path, ratio: float) -> dict[str, pd.DataFrame]:
    out: dict[str, list[pd.DataFrame]] = {}
    for ds in sorted(work.glob(f"ratio{ratio:g}_*")):
        if not (ds / "tests").is_dir():
            continue
        for f in (ds / "tests").glob("*.tsv"):
            t = pd.read_csv(f, sep="\t")
            if f.stem.startswith("inter_"):
                t = t.drop_duplicates("protein")
            out.setdefault(f.stem, []).append(t.assign(dataset=ds.name))
    return {k: pd.concat(v, ignore_index=True) for k, v in out.items()}


def size_adjusted(null_work: Path, alt_work: Path, ratio: float) -> pd.DataFrame:
    """Power at each variant's own p threshold p* that rejects 5 % of the pooled null.

    The null is pooled over its datasets, so each dataset's proteome-wide animal shift counts
    as the population null does. Also printed: the null's median rate over datasets and the
    SE of its pooled rate (between-dataset spread / √datasets), which the shift dominates."""
    null, alt = _tables(null_work, 1.0), _tables(alt_work, ratio)
    rows = []
    for v in alt:
        if v not in null:
            continue
        tests = INTER_TESTS.items() if v.startswith("inter_") else [("dk", "p")]
        for tn, col in tests:
            n, a = null[v][["dataset", col]].dropna(), alt[v][col].dropna()
            per = n.groupby("dataset")[col].apply(lambda x: (x < .05).mean())
            pstar = float(np.quantile(n[col], 0.05))
            rows.append(dict(variant=v.removeprefix("inter_"), test=tn,
                             null_rej05=float((n[col] < .05).mean()),
                             null_se=float(per.std(ddof=1) / np.sqrt(len(per))),
                             null_median_ds=float(per.median()), null_worst_ds=float(per.max()),
                             power=float((a < .05).mean()), p_star=pstar,
                             power_size_adj=float((a < pstar).mean())))
    r = pd.DataFrame(rows)
    for tn, g in r.groupby("test"):
        if "A" in set(g.variant):
            ref = g[g.variant == "A"].iloc[0]
            r.loc[g.index, "power_vs_A"] = g.power / ref.power
            r.loc[g.index, "power_size_adj_vs_A"] = g.power_size_adj / ref.power_size_adj
    r = r.sort_values(["test", "power_size_adj"])
    with pd.option_context("display.width", 250, "display.float_format", "{:.3f}".format):
        print(r.to_string(index=False))
    r.to_csv(alt_work / f"ratio{ratio:g}_size_adjusted.tsv", sep="\t", index=False)
    return r


def summarise(R: pd.DataFrame) -> pd.DataFrame:
    num = [c for c in R.columns if c not in ("test", "variant", "dataset", "ratio", "tests")
           and pd.api.types.is_numeric_dtype(R[c])]
    g = R.groupby(["test", "variant"], sort=False)
    out = g.apply(lambda d: pd.Series({
        "datasets": len(d), "tests": d.tests.mean(),
        **{c: (np.average(d[c], weights=d.tests) if c.startswith(("rej", "cover"))
               else d[c].mean()) for c in num if d[c].notna().all()},
        "rej05_range": f"{d.rej05.min():.3f}–{d.rej05.max():.3f}"}), include_groups=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30,
                           "display.float_format", "{:.3f}".format):
        print(out.to_string())
    return out


if __name__ == "__main__":
    main()
