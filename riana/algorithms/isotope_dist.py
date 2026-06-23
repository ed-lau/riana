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

import math

import numpy as np
import IsoSpecPy
from IsoSpecPy import IsoTotalProb

from riana import constants
from riana.algorithms.mass_calc import count_atoms


def get_peptide_distribution(peptide: str,
                             deuterium_enrichment_level: float = None,
                             label: int = 1,
                             num_labeling_sites: int = 0,
                             mods: list = (),
                             ) -> IsoSpecPy.Iso:

    """
    Calculates the total isotope distribution of a peptide given the peptide sequence and deuterium enrichment level

    :param peptide:                     the peptide sequence
    :param deuterium_enrichment_level:  the deuterium enrichment level of the sample
    :param label:       int: 1=2H_in_vivo, 2=2H_in_vitro, 3=18O, 4=AA, if AA, return 1 assuming no heavy prior to labeling
    :param num_labeling_sites:          the number of labeling sites
    :param mods:        iterable of UniMod accession ids for variable mods on this
                        peptidoform (M7); their atom compositions shape the
                        envelope. Empty by default → bare-backbone envelope.
    :return:                            IsoSpecPy Distribution of atom counts, isotope masses, and isotope probabilities
    """

    # Check that label must be one of hw, hw_cell, or o18
    assert label in [1, 2, 3], 'Label must be one of 1 (2H_in_vivo), 2 (2H_in_vitro), or 3 (18O)'

    if deuterium_enrichment_level is not None:
        assert 0 < deuterium_enrichment_level <= 1, 'Deuterium enrichment level must be greater than 0 and no greater than 1'

    # Get C, H, O, N, S, P count using the Riana count_atoms function
    peptide_atoms = count_atoms(peptide, mods=mods)
    # print(peptide_atoms)

    # Supply atom counts to IsoSpecPy.IsoParamsFromDict and unpack to get atom counts, isotope masses. and probabilities.
    # P (phosphorus) is monoisotopic so a count of 0 leaves the distribution unchanged.
    atom_count_list, isotope_mass_list, isotope_probability_list, _ = IsoSpecPy.IsoParamsFromDict(formula={"C": peptide_atoms[0],
                                                                                                           "H": peptide_atoms[1],
                                                                                                           "O": peptide_atoms[2],
                                                                                                           "N": peptide_atoms[3],
                                                                                                           "S": peptide_atoms[4],
                                                                                                           "P": peptide_atoms[5]})

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
#: Memoized natural-abundance (θ=0) envelope WIDTH (init_w), keyed by
#: (sequence, mods). One small int per peptidoform — see ``init_envelope_width``.
_init_width_cache: dict = {}


def clear_envelope_cache() -> None:
    """Drop the memoized envelopes — for tests / repeated independent runs."""
    _envelope_cache.clear()
    _adaptive_cache.clear()
    _init_width_cache.clear()


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
                  n: int = _DEFAULT_N_ISO,
                  mods: tuple[int, ...] = ()) -> np.ndarray:
    """Natural-abundance envelope, cached by (sequence, mods, n)."""
    key = ('init', sequence, mods, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution(sequence, label=1, mods=mods)
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return _envelope_cache[key]


def init_envelope_width(sequence: str, pep_mass: float, *,
                        mods: tuple[int, ...] = (),
                        abundance_floor: float = 0.01, n: int = 12) -> int:
    """Last isotopomer index whose relative abundance clears ``abundance_floor``
    in the **natural-abundance (unlabeled, θ=0) envelope** — i.e. channels
    ``0..init_w`` are populated by the peptide's own isotope pattern alone.

    This is a *purely compositional* width: it depends only on the peptide's
    atoms (via the same `_get_init_env` the FS solver already builds), with **no
    RIA, no labelling sites, and no Commerford union** — unlike the integrate-time
    "N_ISO" from ``adaptive_channel_masses``, which unions this init envelope with
    the fully-labelled final and so *grows with enrichment*. That RIA-independence
    is exactly why it is the criterion used to widen `--fs` scoring (see
    ``core.fitting.FS_AUTO_*``): a channel inside the natural envelope carries
    clean, model-predicted signal at *every* timepoint, so scoring it is safe
    regardless of θ or RIA.

    One init-envelope computation per (sequence, mods) — cached as a single int,
    so it is amortized over the per-peptide bootstrap's many solve calls. ``n``
    only needs to exceed the widen threshold; the exact value past it is unused.
    """
    key = (sequence, tuple(mods), round(float(abundance_floor), 6), int(n))
    cached = _init_width_cache.get(key)
    if cached is not None:
        return cached
    env = _get_init_env(sequence, pep_mass, n=int(n), mods=tuple(mods))
    tot = float(env.sum()) or 1.0
    last = max((i for i in range(len(env)) if env[i] / tot >= abundance_floor),
               default=0)
    _init_width_cache[key] = last
    return last


def _get_final_env(sequence: str, pep_mass: float, spep: int,
                   ria_max: float,
                   n: int = _DEFAULT_N_ISO,
                   mods: tuple[int, ...] = ()) -> np.ndarray:
    """Fully-labeled envelope at precursor enrichment ``ria_max``
    with ``spep`` labile sites.

    Cached by (sequence, mods, spep, ria_max-rounded, n). ``ria_max`` is the
    experiment's precursor enrichment (~0.06 for 6% v/v D₂O culture
    media, possibly different for in-vivo / metabolic-water cases).
    """
    # Round ria_max for the cache key so very close values share an
    # envelope; the IsoSpec calc is insensitive to 1e-7 changes anyway.
    ria_key = round(float(ria_max), 6)
    key = ('final', sequence, mods, spep, ria_key, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution(
            sequence,
            deuterium_enrichment_level=ria_max,
            label=1,
            num_labeling_sites=spep,
            mods=mods,
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
    ria_max: float = 0.06,
    mods: tuple[int, ...] = (),
) -> float:
    """SSE across all proportions between normalized observed and predicted
    envelopes — the per-peptide Spep fit objective.

    Predicted envelope at a given proportion ``f`` is the linear mixture
    ``(1-f)·init_norm + f·final_norm``. ``final_env`` is linearly
    interpolated between integer-Spep envelopes at ``floor(spep_float)``
    and ``floor(spep_float)+1`` so the loss stays continuous in
    ``spep_float`` — required for :func:`scipy.optimize.minimize_scalar`
    (Brent) to converge cleanly.

    ``ria_max`` is the experiment's precursor enrichment (default 0.06,
    i.e. 6% v/v D₂O); callers in production fitting and benchmark code
    pass the user-controlled value.
    """
    spep_float = max(0.0, spep_float)
    spep_lo = int(np.floor(spep_float))
    spep_hi = spep_lo + 1
    frac = spep_float - spep_lo
    n_iso = obs_matrix.shape[1]

    init_env = _get_init_env(sequence, pep_mass, n=n_iso, mods=mods)
    fenv_lo = _get_final_env(sequence, pep_mass, spep_lo, ria_max, n=n_iso, mods=mods)
    fenv_hi = _get_final_env(sequence, pep_mass, spep_hi, ria_max, n=n_iso, mods=mods)
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


#: Channel count treated as the "full cluster" when normalizing the theoretical
#: init / final envelopes before mixing (the H4′ basis). Wide enough that both
#: the natural-abundance and the fully-labelled (≤ ria_max precursor enrichment)
#: envelopes carry ~all their mass for any realistic peptide — so a convex
#: mixture is itself full-cluster-normalized and per-channel values are true
#: full-cluster fractions. iso_max (adaptive capture cap) is 15 < this.
_FULL_CLUSTER_N = 24


def solve_fs_d2o(
    sequence: str,
    pep_mass: float,
    observed_iso,
    spep: int,
    ria_max: float = 0.06,
    n_iso: int = _DEFAULT_N_ISO,
    mods: tuple[int, ...] = (),
    score_channels: int | None = None,
) -> float:
    """Per-timepoint fractional synthesis from one observed envelope.

    Given a peptide's fitted Spep, solve ``fs`` in ``obs ≈ mixture(fs)`` by
    1-D SSE minimization, using the **H4′ normalization order** (the only order
    that stays linear under truncation):

        mix in the FULL-cluster basis → truncate to the SCORING channels →
        renormalize → compare to the observed renormalized over the same channels.

    IsoSpec returns each envelope already normalized over the whole cluster, so
    a convex combination ``(1-fs)·init_full + fs·final_full`` is itself full-
    cluster-normalized and each channel is a true full-cluster fraction. We then
    truncate that mixture to the scoring channels and renormalize, matching the
    observed cluster (raw intensities) renormalized over the same channels. This
    is exact under any truncation, unlike normalize-each-then-mix (which is only
    correct when init and final share the in-window mass fraction — the H4′
    finding). It is the prerequisite for ``score_channels`` (limited-isotopomer
    scoring): the narrower the subset, the more the two orders diverge.

    Args:
        score_channels: number of leading channels (m0..m{k-1}) to SCORE the fit
            on. ``None`` (default) scores on the peptide's full populated channel
            set. A small value (e.g. 2 = iso0+iso1) dodges co-eluting contaminants
            in the high isotopomers — "integrate wide, fit narrow" (Track B / B4;
            Sadygov & Currie JPR 2025). Capture (integrate) is unaffected; this is
            a fit-time choice.

    Bounds are widened to ``[-0.1, 1.2]`` so unphysical FS (from integration
    noise / co-eluting interference) surfaces in the output instead of being
    clipped — a diagnostic signal the downstream R² curation gate can act on.
    """
    from scipy.optimize import minimize_scalar  # local import keeps cold path fast

    obs = np.asarray(observed_iso, dtype=float)[:n_iso]
    # Adaptive N_ISO: a peptidoform integrated to fewer channels than the run-wide
    # output width carries **trailing NaN** padding. Use only its populated leading
    # channels. On the fixed path ``obs`` has no NaN ⇒ ``n_real == n_iso``.
    valid = ~np.isnan(obs)
    if not valid.any():
        return float('nan')
    n_real = int(np.max(np.nonzero(valid)[0])) + 1
    obs = obs[:n_real]
    if np.isnan(obs).any():
        # Interior NaN (a gap, not trailing padding) — unexpected; bail honestly.
        return float('nan')

    # Scoring width (B4): the chosen subset, clamped per-peptide to what this
    # peptidoform actually populated. A peptide so short its envelope ends before
    # the requested subset (e.g. --fs iso0-3 but only iso0-2 clear 1%) is scored on
    # the channels it has. Need ≥ 2 channels (m0 + a labelled one) to separate init
    # from final; a 1-channel peptidoform cannot be fit → NaN (honest, not a crash).
    k = n_real if score_channels is None else min(int(score_channels), n_real)
    if k < 2:
        return float('nan')
    obs_score = obs[:k]
    obs_total = obs_score.sum()
    if obs_total == 0:
        return float('nan')
    obs_norm = obs_score / obs_total

    # Full-cluster-normalized init / final (mix in this basis, THEN truncate).
    n_full = max(_FULL_CLUSTER_N, n_real)
    init_full = np.asarray(_get_init_env(sequence, pep_mass, n=n_full, mods=mods), dtype=float)
    final_full = np.asarray(
        _get_final_env(sequence, pep_mass, spep, ria_max, n=n_full, mods=mods), dtype=float,
    )
    i_sum = init_full.sum()
    f_sum = final_full.sum()
    if i_sum == 0 or f_sum == 0:
        return float('nan')
    init_full = init_full / i_sum
    final_full = final_full / f_sum

    def sse(fs: float) -> float:
        pred_full = (1.0 - fs) * init_full + fs * final_full
        pred_score = pred_full[:k]
        ps = pred_score.sum()
        if ps <= 0:
            return 1e6
        return float(np.sum((obs_norm - pred_score / ps) ** 2))

    return float(minimize_scalar(sse, bounds=FS_BOUNDS, method='bounded').x)


def fit_peptide_spep(
    sequence: str,
    pep_mass: float,
    obs_matrix: np.ndarray,
    proportions_frac: np.ndarray,
    ria_max: float = 0.06,
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
        args=(sequence, pep_mass, obs_matrix, proportions_frac, ria_max),
        bounds=_SPEP_BOUNDS,
        method='bounded',
    )
    return float(result.x)


# ---------------------------------------------------------------------------
# Adaptive N_ISO — IsoSpec-at-integrate channel selection (Track B)
# ---------------------------------------------------------------------------

#: Adaptive channel-mass cache, keyed by (sequence, mods, ria_key, floor_key,
#: iso_max). Distinct peptidoforms ≪ PSM count, so this is bounded; one entry per
#: peptidoform per run. Per-process under the ``-W`` ProcessPool (each worker
#: rebuilds its own). Cleared by :func:`clear_envelope_cache`.
_adaptive_cache: dict = {}


def _spep_upper_bound(sequence: str) -> int:
    """Conservative labile-H site count for the adaptive *final* envelope.

    Uses the Commerford, Carsten & Cronkite 1983 literature table
    (``constants.label_deuterium_commerford``); its per-residue coefficients
    exceed the calibration/DE-learned ones, so the resulting final envelope is
    the widest plausible and never truncates a real channel. **Width
    determination only** — the exact per-line Spep is fit downstream
    (:func:`spep_from_coefficients`). Rounded up and floored at 1 so even an
    all-low-coefficient peptide carries a labelled channel.
    """
    s = sum(constants.label_deuterium_commerford.get(aa, 0.0) for aa in sequence)
    return max(1, int(math.ceil(s)))


def _binned_envelope(dist, pep_mass: float, n: int) -> tuple[list[float], list[float]]:
    """``(avg_mass, abundance)`` per integer-nominal bin ``round(pep_mass)+iso``.

    ``avg_mass`` is the abundance-weighted average exact mass within the ±0.5 bin
    (the averaged-isotopolog accurate mass — the true centroid of the isotopomer
    cluster); ``abundance`` is the summed probability. An empty bin gives
    ``(nan, 0.0)``. Companion to :func:`get_envelope`, which returns only the
    abundances.
    """
    masses = list(dist.masses)
    probs = list(dist.probs)
    nom0 = round(pep_mass)
    out_m: list[float] = []
    out_p: list[float] = []
    for iso in range(n):
        center = nom0 + iso
        wsum = 0.0
        msum = 0.0
        for m, p in zip(masses, probs):
            if abs(m - center) <= 0.5:
                wsum += p
                msum += m * p
        if wsum > 0:
            out_m.append(msum / wsum)
            out_p.append(wsum)
        else:
            out_m.append(float("nan"))
            out_p.append(0.0)
    return out_m, out_p


def adaptive_channel_masses(
    sequence: str,
    pep_mass: float,
    *,
    ria_max: float = 0.06,
    mods: tuple[int, ...] = (),
    abundance_floor: float = 0.01,
    iso_max: int = 15,
) -> tuple[float, ...]:
    """Per-channel averaged-isotopolog NEUTRAL masses for one peptidoform's
    adaptive isotopomer set (Track B "adaptive N_ISO").

    The channel **count** is the highest isotopomer index whose relative
    abundance clears ``abundance_floor`` in *either* the natural-abundance
    (init) or the fully-labelled (final, at ``ria_max`` with the conservative
    Commerford upper-bound site count) envelope, capped at ``iso_max``. Index 0
    (the precursor m0) is always present. Short peptides get m0-m3, long/
    heavily-labelled ones m0-m8 — instead of a one-size fixed set.

    Each channel's **mass is the init (unlabeled) averaged-isotopolog mass** of
    the nominal bin ``round(pep_mass)+iso`` — the θ=0 reference. (The init/final
    midpoint was considered as a drift-robust target but rejected for v1: it
    overcorrects when labelling is light — low coefficients / low θ keep the real
    peak near init — and anchoring at init makes the recorded ``iso{N}_ppm_error``
    a clean drift-from-unlabelled signal, the orthogonal mass-defect θ estimator.)
    The final envelope is used **only** for the width union here, not the mass.

    Returns a tuple of neutral masses, length = N_channels (= last index + 1).
    Cached by (sequence, mods, ria_max, abundance_floor, iso_max).
    """
    ria_key = round(float(ria_max), 6)
    floor_key = round(float(abundance_floor), 6)
    key = (sequence, tuple(mods), ria_key, floor_key, int(iso_max))
    cached = _adaptive_cache.get(key)
    if cached is not None:
        return cached

    n = int(iso_max) + 1
    spep = _spep_upper_bound(sequence)
    init = get_peptide_distribution(sequence, label=1, mods=tuple(mods))
    # Consistency guard: the envelope is built from sequence + parsed mod atoms,
    # but ``pep_mass`` (the binning anchor) comes from the record. If a mod is in
    # the mass but not the atoms — legacy ``[mass]`` Percolator brackets, which
    # ``parse_unimod_ids`` cannot turn into a composition — the envelope sits a
    # full mod-mass below the anchor and would collapse to iso0. Refuse rather
    # than emit a garbage set: the caller falls back to the fixed analytic
    # channels at the correct ``pep_mass``. (Real [UNIMOD:N] data passes — atoms
    # and mass agree.) Same root cause as the get_envelope binning note in
    # track_c_tmt_chemical_mod.
    env_m0 = min(init.masses)
    if abs(env_m0 - pep_mass) > 0.5:
        raise ValueError(
            f"adaptive envelope m0 {env_m0:.4f} disagrees with pep_mass "
            f"{pep_mass:.4f} by >0.5 Da — unaccounted modification "
            f"(mods={tuple(mods)}); use the fixed channel set."
        )
    final = get_peptide_distribution(
        sequence, deuterium_enrichment_level=ria_max, label=1,
        num_labeling_sites=spep, mods=tuple(mods),
    )
    init_m, init_p = _binned_envelope(init, pep_mass, n)
    _final_m, final_p = _binned_envelope(final, pep_mass, n)

    init_tot = sum(init_p) or 1.0
    final_tot = sum(final_p) or 1.0
    last = 0
    for iso in range(n):
        rel = max(init_p[iso] / init_tot, final_p[iso] / final_tot)
        if rel >= abundance_floor:
            last = iso

    out: list[float] = []
    for iso in range(last + 1):
        im = init_m[iso]
        # A channel within 0..last is populated in the init envelope for any
        # real D₂O distribution (no interior gaps); the analytic-spacing fallback
        # only guards a pathological empty init bin.
        out.append(im if not math.isnan(im) else pep_mass + iso * 1.003354835)
    result = tuple(out)
    _adaptive_cache[key] = result
    return result


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
