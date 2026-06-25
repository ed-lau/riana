"""
Forward-model helpers ported from
data/notebook/90c_O18_LengthModel_IsoSpec_AC16.ipynb (the ¹⁸O reverse model).

NB90c supersedes the deprecated NB90b: it fixes the per-peptide Spep optimizer
(continuous loss via envelope interpolation between integer-Spep envelopes) and
replaces the 20-per-AA coefficient model — which had a K/R non-identifiability
(every tryptic peptide has K+R=1) — with a 5-parameter biochemically-motivated
**length model**::

    Spep = b·(L-1) + c_D·D + c_E·E + c_N·N + c_Q·Q

This module mirrors :mod:`_helpers.forward_model` (the D2O frozen oracle) but for
¹⁸O. Two ¹⁸O-specific differences from D2O:

1. **3-isotope enriched-O pseudo-element** (¹⁶O/¹⁷O/¹⁸O) rather than the 2-isotope
   H/D swap — H₂¹⁸O moves the labeled-O probability mass to ¹⁸O (a +2 Da shift).
2. **N_ISO default 6** (vs 4 for D2O): the +2 Da label spreads probability into
   m2/m4/m6 (singly/doubly/triply-labeled), so iso0:5 all carry Spep information;
   capping lower plateaus the high-Spep fit.

Kept deliberately self-contained (a frozen benchmark oracle), like
:mod:`_helpers.forward_model` — it must not drift when production code changes.
"""
from __future__ import annotations

import re
from typing import Optional

import numpy as np
from IsoSpecPy import IsoParamsFromDict, IsoTotalProb

# ── Physical constants ────────────────────────────────────────────────────────
C_MASS = 12.0000000
H_MASS = 1.00782503223
O_MASS = 15.99491461957
N_MASS = 14.00307400443
S_MASS = 31.9720711744

# ── RIA_O18 — AC16, 6% v/v H2^18O (mole-fraction enrichment of the O pool) ────
_molarity_h2o = 0.997 * 1000 / 18.015
_molarity_o18 = 1.11 * 1000 / 20.015
_o18_vol = 0.06 * 0.97
RIA_O18 = (
    _o18_vol * _molarity_o18
    / ((1 - _o18_vol) * _molarity_h2o + _o18_vol * _molarity_o18)
)

# ── Natural and enriched oxygen isotope abundances ───────────────────────────
# Natural ¹⁶O/¹⁷O/¹⁸O. The enriched pool keeps the ¹⁶O:¹⁷O ratio but lifts ¹⁸O by
# the labeling fraction — H₂¹⁸O dilutes the natural pool and adds ¹⁸O.
NAT_16O, NAT_17O, NAT_18O = 0.99757, 0.00038, 0.00205
P16_LAB = NAT_16O * (1 - RIA_O18)
P17_LAB = NAT_17O * (1 - RIA_O18)
P18_LAB = NAT_18O + (1 - NAT_18O) * RIA_O18

# ── AA elemental composition [C, H, O, N, S] — C = carbamidomethyl Cys ───────
AA_ATOMS = {
    'A': [3, 5, 1, 1, 0],  'C': [5, 8, 2, 2, 1],  'D': [4, 5, 3, 1, 0],
    'E': [5, 7, 3, 1, 0],  'F': [9, 9, 1, 1, 0],  'G': [2, 3, 1, 1, 0],
    'H': [6, 7, 1, 3, 0],  'I': [6, 11, 1, 1, 0], 'K': [6, 12, 1, 2, 0],
    'L': [6, 11, 1, 1, 0], 'M': [5, 9, 1, 1, 1],  'N': [4, 6, 2, 2, 0],
    'P': [5, 7, 1, 1, 0],  'Q': [5, 8, 2, 2, 0],  'R': [6, 12, 1, 4, 0],
    'S': [3, 5, 2, 1, 0],  'T': [4, 7, 2, 1, 0],  'V': [5, 9, 1, 1, 0],
    'W': [11, 10, 1, 2, 0], 'Y': [9, 9, 2, 1, 0],
    'U': [0, 0, 0, 0, 0],  'X': [0, 0, 0, 0, 0],  'B': [0, 0, 0, 0, 0],
}
AA_LIST = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
]

# Default number of isotopomers used in fits and recovery. NB90c sets this to 6
# for ¹⁸O (the N_ISO sweep keeps improving through 6). Call sites accept an
# explicit `n_iso` and only fall back to this default when unset.
N_ISO = 6

# ── Length-model feature spec (NB90c) ────────────────────────────────────────
# Spep = b·(L-1) + c_D·D + c_E·E + c_N·N + c_Q·Q ; intercept fixed at 0.
# Bounds from oxygen counts: b ∈ [0,1] (one backbone carbonyl O per residue),
# c_D,c_E ∈ [0,2] (Asp/Glu carry 2 extra side-chain O), c_N,c_Q ∈ [0,1].
FEATURE_COLS = ['length_minus1', 'D', 'E', 'N', 'Q']
FEATURE_BOUNDS_LOW = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
FEATURE_BOUNDS_HIGH = np.array([1.0, 2.0, 2.0, 1.0, 1.0])


# ── Sequence cleaning ─────────────────────────────────────────────────────────

def clean_seq(seq: str) -> str:
    """Strip ``[mod]`` brackets and a leading ``n`` so length/atom counts match
    what IsoSpec sees — a carbamidomethyl ``[57.0215]`` would otherwise inflate
    ``len()`` by ~9 chars and bias the length-model feature."""
    return re.sub(r'\[.*?\]', '', re.sub('^n', '', seq))


# ── Atom count & ion m/z helpers ──────────────────────────────────────────────

def _count_residue_atoms(seq: str) -> list[int]:
    tot = [0, 0, 0, 0, 0]
    for ch in seq:
        if ch in AA_ATOMS:
            tot = [tot[i] + AA_ATOMS[ch][i] for i in range(5)]
    return tot


def count_atoms(sequence: str, charge: int = 0) -> list[int]:
    """Residue sum + H2O (peptide-bond termini) + charge-state protons."""
    res = _count_residue_atoms(clean_seq(sequence))
    atoms = [res[i] + [0, 2, 1, 0, 0][i] for i in range(5)]  # +H2O
    atoms[1] += charge   # one extra H per charge-state proton
    return atoms


def calculate_ion_mz(seq: str, charge: int = 0) -> float:
    """Monoisotopic neutral mass including H2O and charge-state protons."""
    atoms = count_atoms(seq, charge=charge)
    return (
        atoms[0] * C_MASS + atoms[1] * H_MASS + atoms[2] * O_MASS
        + atoms[3] * N_MASS + atoms[4] * S_MASS
    )


# ── IsoSpecPy envelope functions for O18 ─────────────────────────────────────

def get_peptide_distribution_o18(
    peptide: str,
    charge: int = 0,
    o18_enrichment_level: Optional[float] = None,
    num_o18_labeling_sites: int = 0,
):
    """
    IsoSpecPy distribution for an ¹⁸O-labeled (or natural-abundance) peptide.
    When ``o18_enrichment_level`` is None the natural-abundance distribution is
    returned. ``num_o18_labeling_sites`` is the integer Spep — count of oxygen
    sites exchanged into the enriched (¹⁶O/¹⁷O/¹⁸O) pool. Unlike the D2O 2-isotope
    H/D swap, the labeled O is a 3-isotope pseudo-element.
    """
    peptide_atoms = count_atoms(peptide, charge=charge)
    ac, im, ip, _ = IsoParamsFromDict(formula={
        "C": peptide_atoms[0], "H": peptide_atoms[1],
        "O": peptide_atoms[2], "N": peptide_atoms[3], "S": peptide_atoms[4],
    })

    if num_o18_labeling_sites > 0 and o18_enrichment_level is not None:
        ac = list(ac); im = list(im); ip = list(ip)
        ac[2] -= num_o18_labeling_sites          # reduce natural O count
        ac.append(num_o18_labeling_sites)        # add the labeled-O element
        im.append(im[2])                         # labeled O shares O's mass table
        ip.append((P16_LAB, P17_LAB, P18_LAB))   # 3-isotope enriched probabilities

    return IsoTotalProb(prob_to_cover=0.999, atomCounts=ac,
                        isotopeMasses=im, isotopeProbabilities=ip,
                        use_nominal_masses=True)


def get_envelope(dist, pep_mass: float, n: int = 8) -> list[float]:
    """
    Bin IsoSpecPy distribution (use_nominal_masses=True) into n isotope peaks
    around round(pep_mass)+iso with ±0.5 tolerance, so all peaks are captured
    regardless of mass defect.
    """
    nom_mass_0 = round(pep_mass)
    return [
        sum(p for m, p in zip(dist.masses, dist.probs)
            if abs(m - (nom_mass_0 + iso)) <= 0.5)
        for iso in range(n)
    ]


# ── Cached envelopes + per-peptide Spep loss ─────────────────────────────────

_envelope_cache: dict = {}


def clear_envelope_cache() -> None:
    _envelope_cache.clear()


def _get_init_env(sequence: str, charge: int, pep_mass: float,
                  n: int = N_ISO) -> np.ndarray:
    key = ('init', sequence, charge, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution_o18(sequence, charge=charge)
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return _envelope_cache[key]


def _get_final_env(sequence: str, charge: int, pep_mass: float,
                   spep: int, n: int = N_ISO) -> np.ndarray:
    key = ('final', sequence, charge, spep, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution_o18(
            sequence, charge=charge,
            o18_enrichment_level=RIA_O18,
            num_o18_labeling_sites=spep,
        )
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return _envelope_cache[key]


def peptide_spep_loss(
    spep_float: float,
    sequence: str,
    charge: int,
    pep_mass: float,
    obs_matrix: np.ndarray,        # shape (n_prop, n_iso), normalized per row
    proportions_frac: np.ndarray,  # shape (n_prop,)
) -> float:
    """
    SSE across all proportions between normalized observed and predicted
    envelopes. The predicted final envelope is linearly interpolated between
    integer-Spep envelopes at floor(spep_float) and floor(spep_float)+1 — the
    loss is continuous in spep_float (required for minimize_scalar / Brent) and
    physically equals the population-average envelope when a fraction `frac` of
    molecules carry ceil(s) sites and (1-frac) carry floor(s). Only the first
    n_iso isotopomers are compared (limited-isotopomer method).
    """
    spep_float = max(0.0, spep_float)
    spep_lo = int(np.floor(spep_float))
    spep_hi = spep_lo + 1
    frac = spep_float - spep_lo
    n_iso = obs_matrix.shape[1]

    init_env = _get_init_env(sequence, charge, pep_mass, n=n_iso)
    fenv_lo = _get_final_env(sequence, charge, pep_mass, spep_lo, n=n_iso)
    fenv_hi = _get_final_env(sequence, charge, pep_mass, spep_hi, n=n_iso)
    final_env = (1.0 - frac) * fenv_lo + frac * fenv_hi

    init_sum = init_env.sum()
    final_sum = final_env.sum()
    if init_sum == 0 or final_sum == 0:
        return 1e6

    init_norm = init_env / init_sum
    final_norm = final_env / final_sum

    total_sse = 0.0
    for fs, obs_norm in zip(proportions_frac, obs_matrix):
        pred = (1.0 - fs) * init_norm + fs * final_norm
        total_sse += float(np.sum((obs_norm - pred) ** 2))
    return total_sse


# ── Length-model Spep evaluator + FS recovery solver ─────────────────────────

def estimate_spep_from_length(sequence: str, coef) -> int:
    """Spep from the 5-param length model coefficients ``(b, c_D, c_E, c_N, c_Q)``.

    Floored at 1. ``sequence`` is cleaned of mods first so the length feature is
    consistent with the design matrix used at training.
    """
    b, c_D, c_E, c_N, c_Q = coef
    seq = clean_seq(sequence)
    L = len(seq)
    raw = (b * (L - 1) + c_D * seq.count('D') + c_E * seq.count('E')
           + c_N * seq.count('N') + c_Q * seq.count('Q'))
    return max(1, round(raw))


FS_BOUNDS = (-0.1, 1.2)


def solve_fs_o18(
    sequence: str,
    charge: int,
    pep_mass: float,
    observed_iso,
    spep: int,
    n_iso: int = N_ISO,
) -> float:
    """
    Solve fractional synthesis fs from one observed envelope, given the
    per-peptide Spep. Bounds widened to (-0.1, 1.2) — unphysical fs is preserved
    (not clipped) to expose model misspecification at p=0 / p=100.
    """
    from scipy.optimize import minimize_scalar  # local import — fast cold path

    init_env = np.array(_get_init_env(sequence, charge, pep_mass, n=n_iso))
    final_env = np.array(_get_final_env(sequence, charge, pep_mass, spep, n=n_iso))
    obs = np.array(observed_iso[:n_iso], dtype=float)
    obs_total = obs.sum()
    if obs_total == 0:
        return float('nan')
    obs_norm = obs / obs_total
    i_sum = init_env.sum()
    f_sum = final_env.sum()
    if i_sum == 0 or f_sum == 0:
        return float('nan')
    init_norm = init_env / i_sum
    final_norm = final_env / f_sum

    def sse(fs):
        pred = (1 - fs) * init_norm + fs * final_norm
        return float(np.sum((obs_norm - pred) ** 2))

    return float(minimize_scalar(sse, bounds=FS_BOUNDS, method='bounded').x)
