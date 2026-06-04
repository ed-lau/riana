# -*- coding: utf-8 -*-

"""Theoretical isotope distributions and the per-peptide Spep / FS solver.

Two layers:

1. ``get_peptide_distribution`` — production forward model, lifted unchanged
   from ``riana.utils.get_peptide_distribution`` in M3. Supports label types
   1 (D₂O in vivo), 2 (D₂O in vitro), 3 (¹⁸O). Used by ``riana fit``.
2. M3 Week 4 solver layer: ``get_envelope``, ``peptide_spep_loss``,
   ``solve_fs_d2o`` — lifted from the M2 benchmark's
   ``tests/benchmark/_helpers/forward_model.py``. They compose on top of
   ``get_peptide_distribution(label=1)`` and use the integer-nominal ±0.5
   envelope-binning convention (the M2 oracle's choice; validated by
   ``bench_fs_recovery``).

The M2 benchmark copy in ``_helpers/forward_model.py`` **stays frozen** as
the regression oracle — it must not import this module. The two are
independent implementations of the same forward model; that independence
is the whole point of the M2 calibration gate.

The legacy fixed-site-count fsynthesis (``core/fsynthesis.py``,
``H_MASS``-step ±0.1 binning) is what Week 4 replaces; the new
``core/fitting.py`` calls into this solver layer instead.
"""

import numpy as np
import IsoSpecPy
from IsoSpecPy import IsoTotalProb

from riana.algorithms.mass_calc import count_atoms
from riana.constants import RIA_D2O


def get_peptide_distribution(peptide: str,
                             deuterium_enrichment_level: float = None,
                             label: int = 1,
                             num_labeling_sites: int = 0,
                             ) -> IsoSpecPy.Iso:

    """
    Calculates the total isotope distribution of a peptide given the peptide sequence and deuterium enrichment level

    :param peptide:                     the peptide sequence
    :param deuterium_enrichment_level:  the deuterium enrichment level of the sample
    :param label:       int: 1=2H_in_vivo, 2=2H_in_vitro, 3=18O, 4=AA, if AA, return 1 assuming no heavy prior to labeling
    :param num_labeling_sites:          the number of labeling sites
    :return:                            IsoSpecPy Distribution of atom counts, isotope masses, and isotope probabilities
    """

    # Check that label must be one of hw, hw_cell, or o18
    assert label in [1, 2, 3], 'Label must be one of 1 (2H_in_vivo), 2 (2H_in_vitro), or 3 (18O)'

    if deuterium_enrichment_level is not None:
        assert 0 < deuterium_enrichment_level <= 1, 'Deuterium enrichment level must be greater than 0 and no greater than 1'

    # Get C, H, O, N, S count using the Riana count_atoms function
    peptide_atoms = count_atoms(peptide)
    # print(peptide_atoms)

    # Supply atom counts to IsoSpecPy.IsoParamsFromDict and unpack to get atom counts, isotope masses. and probabilities
    atom_count_list, isotope_mass_list, isotope_probability_list, _ = IsoSpecPy.IsoParamsFromDict(formula={"C": peptide_atoms[0],
                                                                                                           "H": peptide_atoms[1],
                                                                                                           "O": peptide_atoms[2],
                                                                                                           "N": peptide_atoms[3],
                                                                                                           "S": peptide_atoms[4]})

    if label == 1 or label == 2:
        # Subtract the number of labeling sites from hydrogen, extend the atom count list with accessible deuterium count
        atom_count_list[1] = atom_count_list[1] - num_labeling_sites
        atom_count_list.extend([num_labeling_sites])
        # print(f'Atom count list: {atom_count_list}')
        # Extend the isotope mass list for deuterium, which is the same as hydrogen
        isotope_mass_list.extend([isotope_mass_list[1]])
        # print(f'Isotope mass list: {isotope_mass_list}')

        # Extend the isotope probabilities for labelable hydrogen sites, which is the isotope enrichment level
        # For the unlabeled samples, we should use the background level of 0.0001157
        if deuterium_enrichment_level is None:
            isotope_probability_list.extend([isotope_probability_list[1]])
        else:
            isotope_probability_list.extend([(1-deuterium_enrichment_level, deuterium_enrichment_level)])
            # TODO: include the background deuterium level here too?

    elif label == 3:
        atom_count_list[2] = atom_count_list[2] - num_labeling_sites
        atom_count_list.extend([num_labeling_sites])
        # Extend the isotope mass list for O18, which is the same as O16
        isotope_mass_list.extend([isotope_mass_list[2]])

    # print(f'Isotope probability list: {isotope_probability_list}')

    isotope_dist = IsoTotalProb(prob_to_cover=.999,
                       atomCounts=atom_count_list,
                       isotopeMasses= isotope_mass_list,
                       isotopeProbabilities=isotope_probability_list,
                       use_nominal_masses = True)

    return isotope_dist


# ---------------------------------------------------------------------------
# M3 Week 4 — D₂O Spep / FS solver (lifted from M2 forward_model.py)
# ---------------------------------------------------------------------------

# Default envelope length for solver calls. M2's N_ISO sweep landed on 4 as
# the optimum; callers pass an explicit ``n_iso`` for the M3 N_ISO sweep.
_DEFAULT_N_ISO = 4

#: Spep optimization bounds (Brent on a scalar). Realistic peptides have
#: 5-30 labile H sites; wider bounds let the loss surface determine the
#: physical answer without clipping.
_SPEP_BOUNDS = (0.0, 100.0)

#: FS solver bounds — kept slightly wider than [0, 1] so unphysical
#: integrations (e.g. iso0 area drift) are exposed at fit time rather
#: than silently clipped.
FS_BOUNDS = (-0.1, 1.2)

# Module-level memoization. Same peptide / charge / spep produces the same
# envelope; caching avoids re-running IsoSpec inside the Spep optimization
# loop (which evaluates the loss at ~10s of spep candidates per peptide).
# Thread-safe enough for ThreadPoolExecutor — concurrent writes of the same
# key are idempotent at CPython dict level; worst case is duplicated work,
# not corruption. Call ``clear_envelope_cache()`` between independent runs.
_envelope_cache: dict = {}


def clear_envelope_cache() -> None:
    """Drop the memoized envelopes — for tests / repeated independent runs."""
    _envelope_cache.clear()


def get_envelope(dist, pep_mass: float, n: int = 8) -> list[float]:
    """Bin an IsoSpecPy distribution into ``n`` integer-nominal isotope peaks.

    ``use_nominal_masses=True`` in :func:`get_peptide_distribution` returns
    probability-weighted average exact masses per nominal integer bin, so
    we bin around ``round(pep_mass) + iso`` with ±0.5 tolerance — captures
    every peak regardless of mass defect.

    This is the M2 convention (validated by ``bench_fs_recovery``). The
    older ``core/fsynthesis.py`` used ``H_MASS``-step ±0.1, which broke for
    peptides where the mass defect pushes the apex outside the narrow
    tolerance. The new solver standardizes on integer-nominal ±0.5.
    """
    nom_mass_0 = round(pep_mass)
    return [
        sum(p for m, p in zip(dist.masses, dist.probs)
            if abs(m - (nom_mass_0 + iso)) <= 0.5)
        for iso in range(n)
    ]


def _get_init_env(sequence: str, pep_mass: float,
                  n: int = _DEFAULT_N_ISO) -> np.ndarray:
    """Natural-abundance envelope, cached by (sequence, n)."""
    key = ('init', sequence, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution(sequence, label=1)
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return _envelope_cache[key]


def _get_final_env(sequence: str, pep_mass: float, spep: int,
                   n: int = _DEFAULT_N_ISO) -> np.ndarray:
    """Fully-labeled envelope at ``RIA_D2O`` enrichment with ``spep`` sites.

    Cached by (sequence, spep, n). Spep is rounded to integer here because
    the Spep optimization loop linearly interpolates between adjacent
    integer-Spep envelopes (see :func:`peptide_spep_loss`).
    """
    key = ('final', sequence, spep, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution(
            sequence,
            deuterium_enrichment_level=RIA_D2O,
            label=1,
            num_labeling_sites=spep,
        )
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return _envelope_cache[key]


def peptide_spep_loss(
    spep_float: float,
    sequence: str,
    pep_mass: float,
    obs_matrix: np.ndarray,        # shape (n_prop, n_iso), normalized per row
    proportions_frac: np.ndarray,  # shape (n_prop,)
) -> float:
    """SSE across all proportions between normalized observed and predicted
    envelopes — the per-peptide Spep fit objective.

    Predicted envelope at a given proportion ``f`` is the linear mixture
    ``(1-f)·init_norm + f·final_norm``. ``final_env`` is linearly
    interpolated between integer-Spep envelopes at ``floor(spep_float)``
    and ``floor(spep_float)+1`` so the loss stays continuous in
    ``spep_float`` — required for :func:`scipy.optimize.minimize_scalar`
    (Brent) to converge cleanly.
    """
    spep_float = max(0.0, spep_float)
    spep_lo = int(np.floor(spep_float))
    spep_hi = spep_lo + 1
    frac = spep_float - spep_lo
    n_iso = obs_matrix.shape[1]

    init_env = _get_init_env(sequence, pep_mass, n=n_iso)
    fenv_lo = _get_final_env(sequence, pep_mass, spep_lo, n=n_iso)
    fenv_hi = _get_final_env(sequence, pep_mass, spep_hi, n=n_iso)
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


def solve_fs_d2o(
    sequence: str,
    pep_mass: float,
    observed_iso,
    spep: int,
    n_iso: int = _DEFAULT_N_ISO,
) -> float:
    """Per-timepoint fractional synthesis from one observed envelope.

    Given a peptide's fitted Spep (from :func:`peptide_spep_loss`-minimized
    optimization), this solves ``obs_norm ≈ (1-fs)·init_norm + fs·final_norm``
    for ``fs`` via 1-D minimization on the SSE.

    Bounds are widened to ``[-0.1, 1.2]`` so unphysical FS (from
    integration noise / co-eluting interference) surfaces in the output
    instead of being clipped — a diagnostic signal the downstream
    R² curation gate can act on.
    """
    from scipy.optimize import minimize_scalar  # local import keeps cold path fast

    init_env = np.asarray(_get_init_env(sequence, pep_mass, n=n_iso), dtype=float)
    final_env = np.asarray(_get_final_env(sequence, pep_mass, spep, n=n_iso), dtype=float)
    obs = np.asarray(observed_iso[:n_iso], dtype=float)
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

    def sse(fs: float) -> float:
        pred = (1.0 - fs) * init_norm + fs * final_norm
        return float(np.sum((obs_norm - pred) ** 2))

    return float(minimize_scalar(sse, bounds=FS_BOUNDS, method='bounded').x)


def fit_peptide_spep(
    sequence: str,
    pep_mass: float,
    obs_matrix: np.ndarray,
    proportions_frac: np.ndarray,
) -> float:
    """Fit the effective per-peptide labeling-site count Spep.

    **Calibration-only.** This function requires **known proportions** —
    the per-sample nominal heavy fractions from a D₂O mixing-calibration
    experiment, where each ``proportions_frac[i]`` is independently
    measured (nominal mixing ratio), not derived from the kinetic model.

    Used by :mod:`tests.benchmark.bench_aa_coefficients` to learn the
    per-amino-acid coefficient table (``d2o_aa_coefficients_<line>.csv``)
    from the mixing series. The production ``core/fitting.py`` does
    **not** call this — for real user time-series data, ``proportions``
    are unknown a priori (they're what we're fitting). Production
    instead computes ``Spep = Σ aa_coefficient[c] * count(c, sequence)``
    using a pre-learned coefficient table.

    Args:
        sequence: bare amino-acid sequence.
        pep_mass: peptide neutral monoisotopic mass (for nominal-bin centering).
        obs_matrix: per-(proportion, isotope) integrated areas, normalized
            row-wise to sum to 1. Shape ``(n_prop, n_iso)``.
        proportions_frac: per-proportion nominal heavy fraction in
            ``[0, 1]`` from the mixing experiment ground truth.
    """
    from scipy.optimize import minimize_scalar  # local import

    result = minimize_scalar(
        peptide_spep_loss,
        args=(sequence, pep_mass, obs_matrix, proportions_frac),
        bounds=_SPEP_BOUNDS,
        method='bounded',
    )
    return float(result.x)


def spep_from_coefficients(
    sequence: str,
    aa_coefficients: dict[str, float],
    default_per_residue: float = 0.0,
) -> float:
    """Per-peptide Spep from a pre-learned per-amino-acid coefficient table.

    This is the **production** path used by ``core/fitting.py`` for real
    user samples (real D₂O time series with no ground-truth proportions).

    The coefficient table is typically the per-cell-line frozen table
    learned by :mod:`tests.benchmark.bench_aa_coefficients` and committed
    as ``d2o_aa_coefficients_<line>.csv`` (M2 work). For data from a cell
    type without its own calibration, fall back to literature mammalian
    values (Commerford 1983 et al.) — supply that table here.

    Args:
        sequence: bare amino-acid sequence (no modifications encoded).
        aa_coefficients: ``{aa_letter: coefficient}`` — per-residue
            expected labile-H labeling-site count.
        default_per_residue: fallback coefficient for residues absent
            from the table (e.g. when the table omits selenocysteine).
            Defaults to ``0.0`` (effectively ignore unknown residues).

    Returns:
        Sum of per-residue coefficients across the peptide sequence.
    """
    return float(sum(
        aa_coefficients.get(aa, default_per_residue) for aa in sequence
    ))
