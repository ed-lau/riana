"""Dimethyl 0/8 duplex — light→heavy spillover go/no-go (1.2.0 multiplexing).

The dimethyl-duplex (`data/singlepoint_dimethyld2o`) tags control liver with light
dimethyl (UNIMOD:36, +28.0313/site) and rapamycin liver with the +8 heavy dimethyl
(UNIMOD:330, +36.0757/site). Heavy is **+8.0444 Da per dimethyl site** above light,
and a peptide carries **S = 1 (N-term) + #K** sites. So for an **S=1 peptide**
(C-terminal R, no internal K) the heavy cluster's iso0 lands on **light iso8**, and
the question is whether the light cluster — broadened by D₂O labelling — spills
non-negligible intensity into the heavy FS-scoring window (light indices 8..13).

Low spill → the 0/8 duplex is separable and the near-term duplex build
(channel→sample, fit each channel separately) is sound as-is. High spill → it must
ship the light→heavy subtraction (demux; Track C "proper demultiplexing") before the
heavy FS is trustworthy.

This uses RIANA's **production** forward model
(`isotope_dist.get_peptide_distribution` / `get_envelope`) over an in-silico tryptic
digest of the mouse proteome, at the dataset's actual **light** enrichment
(4.614 % D₂O), swept across day-8 fraction-new (FS). FS=1.0 is the conservative worst
case (maximum broadening); day-8 mouse liver is a high-FS regime. The observed light
envelope at fraction FS is the production mixture `(1-FS)·init + FS·final`.

Decision metric: for an S=1 peptide at balanced abundance (L=H), the fraction of the
HEAVY iso0..iso5 window that is actually **light** contamination,
`light_in_window / (light_in_window + heavy_in_window)`. It scales linearly with the
true L/H ratio. See `reports/2026-07-06_dimethyl_duplex_spillover.md`.

Usage:
    python -m tests.benchmark.bench_dimethyl_spillover
    python -m tests.benchmark.bench_dimethyl_spillover --coefficients ilchenko_2019
    python -m tests.benchmark.bench_dimethyl_spillover --fasta /path/to/proteome.fasta
"""
from __future__ import annotations

import argparse

import numpy as np

from riana.algorithms.isotope_dist import (
    get_envelope,
    get_peptide_distribution,
    spep_from_coefficients,
)
from riana.algorithms.mass_calc import _calc_atom_mass, count_atoms
from riana.core.fitting import load_aa_coefficients

RIA_LIGHT, RIA_HEAVY = 0.04614, 0.05611          # per-channel D₂O enrichment (SDRF)
DIMETHYL_LIGHT = 36                              # UNIMOD:36, mod_atoms[36]=C2H4
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
N_WIDE = 28                                      # full-cluster normalization width
HEAVY_WIN = 6                                    # heavy FS solver scores iso0..iso5
LEN_BINS = [(7, 12), (13, 18), (19, 25), (26, 35)]
FS_SWEEP = [0.5, 0.7, 0.9, 1.0]
# approximate DDA-detectability weights per length bin (short peptides dominate)
LEN_WEIGHTS = {(7, 12): 0.40, (13, 18): 0.35, (19, 25): 0.18, (26, 35): 0.07}
DEFAULT_FASTA = "data/timeseries_lve/uniprot_mouse_reviewed.fasta"


def read_fasta(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def digest(protein, missed=2):
    """Trypsin: cut after K/R not before P; 0..`missed` missed cleavages."""
    sites = ([0] + [i + 1 for i in range(len(protein) - 1)
                    if protein[i] in "KR" and protein[i + 1] != "P"] + [len(protein)])
    peps = set()
    for a in range(len(sites) - 1):
        for b in range(a + 1, min(a + 2 + missed, len(sites))):
            peps.add(protein[sites[a]:sites[b]])
    return peps


def light_mixture(seq, ria, fs, coeffs):
    """Normalized (sum=1) mixture envelope for the S=1 dimethylated peptidoform on
    the light nominal comb — the production `(1-fs)·init + fs·final`."""
    mods = [DIMETHYL_LIGHT]
    pm = _calc_atom_mass(count_atoms(seq, mods=mods))
    spep = max(1, int(round(spep_from_coefficients(seq, coeffs))))
    init = np.array(get_envelope(
        get_peptide_distribution(seq, label="D2O", mods=mods), pm, n=N_WIDE), float)
    fin = np.array(get_envelope(get_peptide_distribution(
        seq, deuterium_enrichment_level=ria, label="D2O",
        num_labeling_sites=spep, mods=mods), pm, n=N_WIDE), float)
    init /= init.sum() or 1.0
    fin /= fin.sum() or 1.0
    return (1.0 - fs) * init + fs * fin


def q(a, p):
    return float(np.quantile(a, p)) if len(a) else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", default=DEFAULT_FASTA)
    ap.add_argument("--coefficients", default="deberneh_2025_rss",
                    help="bundled preset or CSV path (mouse in-vivo D₂O default)")
    ap.add_argument("--sample-per-bin", type=int, default=1200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    coeffs = load_aa_coefficients(args.coefficients)
    rng = np.random.default_rng(args.seed)

    print("Digesting proteome ...", flush=True)
    allpeps = set()
    for prot in read_fasta(args.fasta):
        allpeps |= digest("".join(c for c in prot if c in STD_AA))
    allp = [p for p in allpeps if 7 <= len(p) <= 35 and set(p) <= STD_AA]
    s1 = [p for p in allp if p.count("K") == 0]      # S=1: heavy at +8 Da (worst case)
    print(f"  tryptic len7-35: {len(allp):,}   S=1 (no K): {len(s1):,} "
          f"({100*len(s1)/len(allp):.1f}%)   S≥2: {len(allp)-len(s1):,} "
          f"({100*(len(allp)-len(s1))/len(allp):.1f}%, heavy at +16 Da → no overlap)")

    data = {}
    for lo, hi in LEN_BINS:
        bp = [p for p in s1 if lo <= len(p) <= hi]
        idx = rng.choice(len(bp), size=min(args.sample_per_bin, len(bp)), replace=False)
        bp = [bp[i] for i in idx]
        for fs in FS_SWEEP:
            contam = []
            for seq in bp:
                try:
                    lmix = light_mixture(seq, RIA_LIGHT, fs, coeffs)   # light cluster
                    hmix = light_mixture(seq, RIA_HEAVY, fs, coeffs)   # heavy shape
                except Exception:
                    continue
                lin = lmix[8:8 + HEAVY_WIN].sum()      # light contam in heavy window
                hin = hmix[0:HEAVY_WIN].sum()          # heavy signal in its window
                denom = lin + hin
                contam.append(100 * lin / denom if denom > 0 else 0.0)
            data[(lo, hi, fs)] = (np.array(contam), len(bp))

    print("\n" + "=" * 78)
    print("S=1 — % of the HEAVY iso0..5 window that is LIGHT contamination (L=H)")
    print("=" * 78)
    for lo, hi in LEN_BINS:
        print(f"\n--- length {lo}-{hi}  (n={data[(lo, hi, FS_SWEEP[0])][1]}) ---")
        print("  FS    med    p75    p90    p99    | >5%    >10%   >20%")
        for fs in FS_SWEEP:
            c = data[(lo, hi, fs)][0]
            print(f"  {fs:.1f}  {q(c,.5):5.2f}  {q(c,.75):5.2f}  {q(c,.9):5.2f}  "
                  f"{q(c,.99):5.2f}  | " + "  ".join(
                      f"{100*np.mean(c > t):4.0f}%" for t in (5, 10, 20)))

    print("\n" + "=" * 78)
    print("Population S=1 (length-weighted ~DDA): % heavy windows >X% light")
    print("=" * 78)
    print("   FS    >5%     >10%    >20%")
    for fs in FS_SWEEP:
        row = [sum(LEN_WEIGHTS[(lo, hi)] * np.mean(data[(lo, hi, fs)][0] > t)
                   for lo, hi in LEN_BINS) for t in (5, 10, 20)]
        print(f"  {fs:.1f}   " + "  ".join(f"{100*r:5.1f}%" for r in row))


if __name__ == "__main__":
    main()
