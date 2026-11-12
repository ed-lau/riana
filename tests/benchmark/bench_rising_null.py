"""False positives of the Δk and interaction tests under a rising precursor — Monte Carlo
through the shipped pipeline (``riana fit`` → ``riana rollup``), from simulated envelopes.

The earlier calibration work (``bench_linear_weights.py``, ``notebooks/wls_var_sim.py``)
simulated θ directly with an instant precursor. Under a rising precursor the variance the
rollup weights on comes from a chain those sims do not have: the per-point envelope solve
(``riana.algorithms.precursor``), the peptide's label-content bootstrap band, the gain ∂A/∂f
that maps it to f, the inverse-variance collapse and the ``wls-var`` moderation. So this
bench starts from **envelopes** and runs the real CLI on them.

Design (generic, the shape of an unprimed in-vivo D₂O study): two conditions, three tissues
sampled from the same animals, one animal per condition × day, days 0–84, k_p = 0.08/day,
RIA 0.04. Proteins have 2–6 random tryptic peptides; true k log-uniform 0.004–0.15/day with
fixed tissue factors. Nuisance that real data have:

- envelope noise, relative + absolute per channel, with a per-peptide noise multiplier
  (log-normal) — peptides differ in precision;
- per-peptide labeling-site error (the true site count differs from the table's, ±10 %);
- per-animal enrichment (RIA × (1 + N(0, 5 %)), shared by the animal's tissues) and
  per-animal × protein turnover (k × e^{N(0, 0.15)}, shared by its tissues);
- interference: 5 % of points get +20–100 % on one of m3–m5.

``--ratio``: true k_B / k_A in every tissue. 1 is the null (every Δk and interaction test is a
false positive at rate α if calibrated); above 1 gives power and ratio recovery, and the
ratio-scale interaction stays null (the ratio is the same in every tissue).

Each dataset is fitted twice — θ (``--kp inf``) and f (k_p from the manifest) — and rolled up
with ``wls`` and ``wls-var`` (interaction on). Reported: rejection rate at p < 0.05 / 0.01 of
the per-tissue Δk test and of the interaction's mean-Δk / additive / ratio tests, plus the
k recovery against truth.

    python -m tests.benchmark.bench_rising_null --datasets 4 --proteins 250 -W 14
    python -m tests.benchmark.bench_rising_null --ratio 1.5 --datasets 1

Outputs (large, regenerable) go to ``--work`` (default ``runs/bench_rising_null``, gitignored).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

DAYS = (0.0, 1.0, 3.0, 7.0, 14.0, 21.0, 35.0, 56.0, 84.0)
KP, RIA = 0.08, 0.04
TISSUES = {"tissueA": 1.0, "tissueB": 0.5, "tissueC": 2.0}
COND = ("A", "B")
NOISE_REL, NOISE_ABS, NOISE_SPREAD = 0.035, 0.0025, 0.5
SITE_ERR, RIA_CV, ANIMAL_K_SD = 0.10, 0.05, 0.15
P_INTERFERE = 0.05
P_INTERFERE_RUN = [P_INTERFERE]   # passed to workers in the run dict (they re-import the module)
COEFFS = "deberneh_2025_rss"
_AA = "ACDEFGHIKLMNPQRSTVWY"
_AA_FREQ = np.array([8.3, 1.4, 5.5, 6.8, 3.9, 7.1, 2.3, 5.9, 5.8, 9.7, 2.4, 4.1, 4.7, 3.9,
                     5.4, 6.6, 5.3, 6.9, 1.1, 2.9])


def _random_peptide(rng: np.random.Generator) -> str:
    n = int(rng.integers(9, 21))
    body = rng.choice(list(_AA.replace("K", "").replace("R", "")), size=n - 1,
                      p=_freq_no_kr())
    return "".join(body) + ("K" if rng.random() < 0.5 else "R")


def _freq_no_kr() -> np.ndarray:
    f = np.array([p for a, p in zip(_AA, _AA_FREQ) if a not in "KR"])
    return f / f.sum()


def simulate(work: Path, seed: int, n_proteins: int, ratio: float, workers: int) -> Path:
    """Write per-run integrate files + a manifest; return the manifest path."""
    from riana.algorithms.isotope_dist import spep_from_coefficients
    from riana.core.fitting import load_aa_coefficients
    from riana.io.manifest import ManifestRow, append_manifest
    from riana.records import RunIdentity

    rng = np.random.default_rng(seed)
    coeffs = load_aa_coefficients(COEFFS)
    proteins = []
    for i in range(n_proteins):
        k0 = float(np.exp(rng.uniform(np.log(0.004), np.log(0.15))))
        peps = []
        while len(peps) < int(min(6, 2 + rng.poisson(2))):
            seq = _random_peptide(rng)
            spep_tab = spep_from_coefficients(seq, coeffs)
            if round(spep_tab) < 9:
                continue
            peps.append(dict(seq=seq, spep=max(1, int(round(spep_tab * math.exp(
                rng.normal(0, SITE_ERR))))), noise=float(math.exp(
                    rng.normal(0, NOISE_SPREAD) - NOISE_SPREAD ** 2 / 2)),
                scale=float(math.exp(rng.normal(0, 1.0)))))
        proteins.append(dict(pid=f"sp|SIM{i:05d}|SIM_{i}", k0=k0, peps=peps))
    runs = []
    for cond in COND:
        for day in DAYS:
            animal = f"{cond}_D{int(day)}"
            ria_a = RIA * (1 + rng.normal(0, RIA_CV))
            u = rng.normal(0, ANIMAL_K_SD, n_proteins)       # shared by the animal's tissues
            for tissue, tf in TISSUES.items():
                runs.append(dict(cond=cond, day=day, animal=animal, tissue=tissue,
                                 ria=ria_a, kp=KP, p_int=P_INTERFERE_RUN[0], k=[p["k0"] * tf * (ratio if cond == "B" else 1.0)
                                               * math.exp(u[i]) for i, p in enumerate(proteins)],
                                 seed=int(rng.integers(2**31))))
    work.mkdir(parents=True, exist_ok=True)
    args = [(work, r, proteins) for r in runs]
    with ProcessPoolExecutor(workers) as ex:
        paths = list(ex.map(_write_run, args))
    truth = pd.DataFrame([dict(experiment=t, protein=p["pid"].split("|")[1], condition=c,
                               k_true=p["k0"] * tf * (ratio if c == "B" else 1.0))
                          for p in proteins for t, tf in TISSUES.items() for c in COND])
    truth.to_csv(work / "truth.tsv", sep="\t", index=False)
    mf = work / "riana_manifest.tsv"
    if mf.exists():
        mf.unlink()
    append_manifest(mf, [ManifestRow("integrate", str(path), RunIdentity(
        experiment=r["tissue"], sample=path.stem, data_file=path.stem,
        biological_replicate=1, individual=r["animal"], labeling_time=r["day"],
        labeling_time_unit="day", condition=r["cond"], precursor_enrichment=RIA,
        precursor_rate_constant=KP)) for r, path in zip(runs, paths)])
    return mf


def _write_run(arg) -> Path:
    from riana.algorithms.mass_calc import calculate_ion_mz
    from riana.algorithms.precursor import RisingEnvelope

    work, r, proteins = arg
    rng = np.random.default_rng(r["seed"])
    rows = []
    for i, p in enumerate(proteins):
        f = -math.expm1(-r["k"][i] * r["day"])
        for j, pep in enumerate(p["peps"]):
            env = RisingEnvelope(pep["seq"], calculate_ion_mz(pep["seq"]), pep["spep"],
                                 r["ria"], r["day"], r["kp"])   # explicit: workers re-import
            e = env.predict(f)[:6]
            s = pep["noise"]
            obs = e * (1 + s * NOISE_REL * rng.normal(size=6)) + s * NOISE_ABS * rng.normal(size=6)
            if rng.random() < r["p_int"]:
                obs[int(rng.integers(3, 6))] *= 1 + rng.uniform(0.2, 1.0)
            obs = np.clip(obs, 1e-6, None) * 1e6 * pep["scale"]
            rows.append({"file_idx": 0, "scan": 1000 + len(rows), "charge": 2,
                         "concat": f"{pep['seq']}_2", "sequence": pep["seq"],
                         "sample": "sim", "percolator q-value": 1e-4, "protein id": p["pid"],
                         **{f"iso{c}": obs[c] for c in range(6)}})
    path = work / f"{r['tissue']}_{r['animal']}_riana.txt"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _riana(*args: str) -> None:
    exe = Path(sys.executable).with_name("riana")
    res = subprocess.run([str(exe), *args], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"riana {' '.join(args)} failed:\n{res.stdout[-2000:]}\n{res.stderr[-2000:]}")


def run_dataset(work: Path, seed: int, n_proteins: int, ratio: float, workers: int) -> pd.DataFrame:
    mf = simulate(work, seed, n_proteins, ratio, workers)
    rows = []
    for route, kp_args in (("theta", ["--kp", "inf"]), ("f", [])):
        fit_dir = work / route
        _riana("fit", "--manifest", str(mf), "--coefficients", COEFFS, "--fs", "auto",
               "-W", str(workers), "-o", str(fit_dir), *kp_args)
        for weights in ("wls", "wls-var"):
            out = fit_dir / f"rollup_{weights}"
            _riana("rollup", str(fit_dir), "--model", "linear simple", "--linear-weights",
                   weights, "--reference-condition", "A", "--test-condition", "B",
                   "--interaction", "-W", str(workers), "-o", str(out))
            rows.append(_summarise(work, out, route, weights, seed, ratio))
    return pd.DataFrame(rows)


def _summarise(work: Path, out: Path, route: str, weights: str, seed: int, ratio: float) -> dict:
    rd = lambda p: pd.read_table(p, comment="#", low_memory=False)
    prot = rd(out / "riana_rollup_proteins.txt")
    dk = prot.dropna(subset=["delta_k_p"]).drop_duplicates(["experiment", "protein"])
    inter = rd(out / "riana_rollup_interaction.txt").drop_duplicates("protein")
    truth = rd(work / "truth.tsv")
    ref = prot[prot.condition == "A"].merge(truth, on=["experiment", "protein", "condition"])
    rec = np.log(ref.k_deg / ref.k_true).replace([np.inf, -np.inf], np.nan).dropna()
    rate = lambda p, a: float((p.dropna() < a).mean()) if p.notna().any() else float("nan")
    return dict(seed=seed, ratio=ratio, route=route, weights=weights,
                dk_tests=len(dk), dk_rej05=rate(dk.delta_k_p, .05), dk_rej01=rate(dk.delta_k_p, .01),
                int_proteins=len(inter), mean_rej05=rate(inter.mean_delta_k_p, .05),
                add_rej05=rate(inter.interaction_additive_p, .05),
                ratio_rej05=rate(inter.interaction_p, .05),
                k_over_true=float(np.exp(rec.median())) if len(rec) else float("nan"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--datasets", type=int, default=4)
    ap.add_argument("--proteins", type=int, default=250)
    ap.add_argument("--ratio", type=float, default=1.0,
                    help="true k_B / k_A in every tissue (1 = the null)")
    ap.add_argument("--seed", type=int, default=20261112)
    ap.add_argument("-W", "--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--work", default="runs/bench_rising_null")
    ap.add_argument("--animal-k-sd", type=float, default=ANIMAL_K_SD,
                    help="per-animal × protein turnover SD (log scale); 0 turns it off")
    ap.add_argument("--ria-cv", type=float, default=RIA_CV,
                    help="per-animal enrichment CV; 0 turns it off")
    ap.add_argument("--kp", type=float, default=KP,
                    help="true precursor rate constant (also written to the manifest); a large "
                    "value (e.g. 100) is an effectively instant precursor")
    ap.add_argument("--site-err", type=float, default=SITE_ERR,
                    help="per-peptide labeling-site error SD (log scale); 0 turns it off")
    ap.add_argument("--interfere", type=float, default=P_INTERFERE,
                    help="fraction of points with high-channel interference; 0 turns it off")
    ap.add_argument("--noise-spread", type=float, default=NOISE_SPREAD,
                    help="per-peptide noise-multiplier SD (log scale); 0 = equal precision")
    a = ap.parse_args()
    globals().update(ANIMAL_K_SD=a.animal_k_sd, RIA_CV=a.ria_cv, KP=a.kp, SITE_ERR=a.site_err,
                     NOISE_SPREAD=a.noise_spread)
    P_INTERFERE_RUN[0] = a.interfere
    # Provenance for audit: the exact settings and code version behind this arm's numbers.
    git = lambda *c: subprocess.run(["git", *c], capture_output=True, text=True).stdout.strip()
    Path(a.work).mkdir(parents=True, exist_ok=True)
    (Path(a.work) / f"ratio{a.ratio:g}_params.json").write_text(json.dumps(dict(
        vars(a), argv=sys.argv[1:], commit=git("rev-parse", "--short", "HEAD"),
        code_dirty=bool(git("status", "--porcelain", "--", "riana", "tests/benchmark")),
        days=DAYS, ria=RIA, tissues=TISSUES, noise_rel=NOISE_REL, noise_abs=NOISE_ABS,
        coefficients=COEFFS), indent=2))
    res = []
    for d in range(a.datasets):
        work = Path(a.work) / f"ratio{a.ratio:g}_{d}"
        res.append(run_dataset(work, a.seed + d + int(1000 * a.ratio), a.proteins, a.ratio,
                               a.workers))
        print(f"dataset {d} done", flush=True)
    df = pd.concat(res, ignore_index=True)
    df.to_csv(Path(a.work) / f"ratio{a.ratio:g}_results.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    w = lambda g: pd.Series({
        "dk_tests": g.dk_tests.sum(),
        "dk_rej05": np.average(g.dk_rej05, weights=g.dk_tests),
        "dk_rej01": np.average(g.dk_rej01, weights=g.dk_tests),
        "int_proteins": g.int_proteins.sum(),
        "mean_rej05": np.average(g.mean_rej05, weights=g.int_proteins),
        "add_rej05": np.average(g.add_rej05, weights=g.int_proteins),
        "ratio_rej05": np.average(g.ratio_rej05, weights=g.int_proteins),
        "k_over_true": g.k_over_true.median()})
    print(f"\nratio {a.ratio:g}: rejection rates (Type-I where the null holds)")
    print(df.groupby(["route", "weights"]).apply(w, include_groups=False).round(3).to_string())


if __name__ == "__main__":
    main()
