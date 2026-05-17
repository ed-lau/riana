"""
Forward-model helpers ported from
data/notebook/87a_D2O_LearnAALabelingSites_IsoSpec_AC16.ipynb.

The NB defines an IsoSpecPy-backed forward model for D2O-labeled peptide
isotope envelopes and a continuous Spep loss for per-peptide fitting.
This module exposes the same functions as importable callables.
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

# ── RIA_D2O — 6% v/v D2O ─────────────────────────────────────────────────────
# Mole-fraction enrichment of D in the H pool for 6% v/v D2O in H2O.
_molarity_h2o = 0.997 * 1000 / 18.015
_molarity_d2o = 1.1056 * 1000 / 20.027
_d2o_vol = 0.06 * 0.9996
RIA_D2O = (
    _d2o_vol * _molarity_d2o
    / ((1 - _d2o_vol) * _molarity_h2o + _d2o_vol * _molarity_d2o)
)

# ── RIA_O18 — 6% v/v H2^18O (kept for parity with the NB) ────────────────────
_molarity_o18 = 1.11 * 1000 / 20.015
_o18_vol = 0.06 * 0.97
RIA_O18 = (
    _o18_vol * _molarity_o18
    / ((1 - _o18_vol) * _molarity_h2o + _o18_vol * _molarity_o18)
)

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

# Default number of isotopomers used in fits and recovery. The NB sets this to 4
# and parameterizes the N_ISO sweep around it. Most call sites accept an
# explicit `n_iso` argument and only fall back to this default when unset.
N_ISO = 4


# ── Atom count & ion m/z helpers ──────────────────────────────────────────────

def _count_residue_atoms(seq: str) -> list[int]:
    tot = [0, 0, 0, 0, 0]
    for ch in seq:
        if ch in AA_ATOMS:
            tot = [tot[i] + AA_ATOMS[ch][i] for i in range(5)]
    return tot


def count_atoms(sequence: str, charge: int = 0) -> list[int]:
    """Residue sum + H2O (peptide-bond termini) + charge-state protons."""
    clean = re.sub(r'\[.*?\]', '', re.sub('^n', '', sequence))
    res = _count_residue_atoms(clean)
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


# ── IsoSpecPy envelope functions for D2O ─────────────────────────────────────

def get_peptide_distribution_d2o(
    peptide: str,
    charge: int = 0,
    deuterium_enrichment_level: Optional[float] = None,
    num_deuterium_labeling_sites: int = 0,
):
    """
    IsoSpecPy distribution for a D2O-labeled (or natural-abundance) peptide.
    When deuterium_enrichment_level is None the natural-abundance distribution
    is returned. num_deuterium_labeling_sites is the integer Spep — count of
    labile H sites exchanged with the D pool.
    """
    peptide_atoms = count_atoms(peptide, charge=charge)
    ac, im, ip, _ = IsoParamsFromDict(formula={
        "C": peptide_atoms[0], "H": peptide_atoms[1],
        "O": peptide_atoms[2], "N": peptide_atoms[3], "S": peptide_atoms[4],
    })

    if num_deuterium_labeling_sites > 0 and deuterium_enrichment_level is not None:
        ac = list(ac); im = list(im); ip = list(ip)
        ac[1] -= num_deuterium_labeling_sites
        ac.append(num_deuterium_labeling_sites)
        im.append(im[1])          # D shares H's exact-mass table
        p_H = 1.0 - deuterium_enrichment_level
        p_D = deuterium_enrichment_level
        ip.append((p_H, p_D))

    return IsoTotalProb(prob_to_cover=0.999, atomCounts=ac,
                        isotopeMasses=im, isotopeProbabilities=ip,
                        use_nominal_masses=True)


def get_envelope(dist, pep_mass: float, n: int = 8) -> list[float]:
    """
    Bin IsoSpecPy distribution (use_nominal_masses=True) into n isotope peaks.
    use_nominal_masses=True returns probability-weighted average exact masses
    per nominal integer bin; we bin around round(pep_mass)+iso with ±0.5
    tolerance so all peaks are captured regardless of mass defect.
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
        dist = get_peptide_distribution_d2o(sequence, charge=charge)
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return _envelope_cache[key]


def _get_final_env(sequence: str, charge: int, pep_mass: float,
                   spep: int, n: int = N_ISO) -> np.ndarray:
    key = ('final', sequence, charge, spep, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution_d2o(
            sequence, charge=charge,
            deuterium_enrichment_level=RIA_D2O,
            num_deuterium_labeling_sites=spep,
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
    loss is continuous in spep_float (required for minimize_scalar / Brent).
    Only the first n_iso isotopomers are compared (limited-isotopomer method).
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
    if init_sum == 0:
        return 1e6
    if final_sum == 0:
        return 1e6

    init_norm = init_env / init_sum
    final_norm = final_env / final_sum

    total_sse = 0.0
    for fs, obs_norm in zip(proportions_frac, obs_matrix):
        pred = (1.0 - fs) * init_norm + fs * final_norm
        total_sse += float(np.sum((obs_norm - pred) ** 2))
    return total_sse


# ── FS recovery solver ────────────────────────────────────────────────────────
FS_BOUNDS = (-0.1, 1.2)


def solve_fs_d2o(
    sequence: str,
    charge: int,
    pep_mass: float,
    observed_iso,
    spep: int,
    n_iso: int = N_ISO,
) -> float:
    """
    Solve fractional synthesis fs from a single observed envelope, given the
    per-peptide Spep. Bounds are widened to (-0.1, 1.2) — physically unphysical
    fs outside [0, 1] is preserved (not clipped) to expose model misspec.
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
