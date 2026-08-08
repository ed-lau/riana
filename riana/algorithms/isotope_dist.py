# -*- coding: utf-8 -*-

"""Theoretical isotope distributions and the per-peptide Spep / FS solver.

Two layers:

1. ``get_peptide_distribution`` — production forward model, lifted unchanged
   from ``riana.utils.get_peptide_distribution`` in M3. Labeling chemistry is
   selected by the string ``label``: ``"D2O"`` (heavy water; subsumes the old
   in-vivo/in-vitro split) or ``"O18"`` (¹⁸O). Used by ``riana fit``.
2. M3 Week 4 solver layer: ``get_envelope``, ``peptide_spep_loss``,
   ``solve_fs_d2o`` — lifted from the M2 benchmark's
   ``tests/benchmark/_helpers/forward_model.py``. They compose on top of
   ``get_peptide_distribution(label="D2O")`` and use the integer-nominal ±0.5
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

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any, cast

import numpy as np
import numpy.typing as npt
import IsoSpecPy
from IsoSpecPy import IsoTotalProb

from riana import constants
from riana.algorithms.mass_calc import count_atoms


def get_peptide_distribution(peptide: str,
                             deuterium_enrichment_level: float | None = None,
                             label: str = "D2O",
                             num_labeling_sites: int = 0,
                             mods: Sequence[int] = (),
                             ) -> IsoSpecPy.Iso:

    """
    Calculates the total isotope distribution of a peptide given the peptide sequence and deuterium enrichment level

    :param peptide:                     the peptide sequence
    :param deuterium_enrichment_level:  the labeling enrichment (RIA) of the sample
                                        (the ²H fraction for ``"D2O"``, the ¹⁸O
                                        fraction for ``"O18"``)
    :param label:       str: the labeling chemistry — ``"D2O"`` (heavy water, the
                        ²H H/D swap; subsumes the old in-vivo/in-vitro split) or
                        ``"O18"`` (¹⁸O metabolic labeling, the 3-isotope-O shift)
    :param num_labeling_sites:          the number of labeling sites
    :param mods:        iterable of UniMod accession ids for variable mods on this
                        peptidoform (M7); their atom compositions shape the
                        envelope. Empty by default → bare-backbone envelope.
    :return:                            IsoSpecPy Distribution of atom counts, isotope masses, and isotope probabilities
    """

    assert label in ("D2O", "O18"), 'Label must be "D2O" (heavy water / ²H) or "O18" (¹⁸O)'

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

    # Fixed-isotope modifications (TMT/TMTpro): built-in ¹³C/¹⁵N that are ~100%
    # heavy by synthesis (NOT natural abundance), so they shift both mass and
    # envelope shape. Append each as a single-isotope pseudo-element (prob 1.0) —
    # mass added, zero combinatorial broadening — the same construction the D2O/O18
    # label uses below. ``pep_mass`` (the caller's nominal-bin anchor) includes
    # these via ``mass_calc.unimod_mass``, so the envelope and anchor stay in
    # lockstep. Identical isotope masses are aggregated across all mod instances on
    # the peptidoform (e.g. TMTpro on the N-term AND each lysine).
    pinned_isotopes: dict[float, int] = {}
    for mod_id in mods:
        for count, iso_mass in constants.mod_fixed_isotopes.get(mod_id, ()):
            pinned_isotopes[iso_mass] = pinned_isotopes.get(iso_mass, 0) + count
    for iso_mass, count in pinned_isotopes.items():
        atom_count_list.append(count)
        isotope_mass_list.append((iso_mass,))
        isotope_probability_list.append((1.0,))

    if label == "D2O":
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

    elif label == "O18":
        # ¹⁸O metabolic labeling (H₂¹⁸O): move ``num_labeling_sites`` oxygens into
        # an enriched-O pseudo-element. Unlike the D₂O 2-isotope H/D swap, labeled
        # O is a 3-isotope element (¹⁶O/¹⁷O/¹⁸O) — the enrichment dilutes the
        # natural pool and lifts ¹⁸O by the labeling fraction (the +2 Da shift).
        # ``deuterium_enrichment_level`` carries the o18 RIA here (it is the
        # generic labeling-enrichment level, despite the D₂O-flavoured name).
        # Matches the NB90c reverse model + tests/benchmark/_helpers/o18_forward_model.
        if num_labeling_sites > 0 and deuterium_enrichment_level is not None:
            atom_count_list[2] = atom_count_list[2] - num_labeling_sites
            atom_count_list.extend([num_labeling_sites])
            # Labeled O shares oxygen's exact-mass table (¹⁶O/¹⁷O/¹⁸O).
            isotope_mass_list.extend([isotope_mass_list[2]])
            nat16, nat17, nat18 = isotope_probability_list[2]
            ria = deuterium_enrichment_level
            isotope_probability_list.extend([(
                nat16 * (1 - ria),
                nat17 * (1 - ria),
                nat18 + (1 - nat18) * ria,
            )])

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
_envelope_cache: dict[Any, Any] = {}
#: Memoized natural-abundance (θ=0) envelope WIDTH (init_w), keyed by
#: (sequence, mods). One small int per peptidoform — see ``init_envelope_width``.
_init_width_cache: dict[Any, int] = {}


def clear_envelope_cache() -> None:
    """Drop the memoized envelopes — for tests / repeated independent runs."""
    _envelope_cache.clear()
    _adaptive_cache.clear()
    _init_width_cache.clear()


def get_envelope(dist: IsoSpecPy.Iso, pep_mass: float, n: int = 8) -> list[float]:
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
                  mods: tuple[int, ...] = ()) -> npt.NDArray[np.float64]:
    """Natural-abundance envelope, cached by (sequence, mods, n)."""
    key = ('init', sequence, mods, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution(sequence, label="D2O", mods=mods)
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return cast("npt.NDArray[np.float64]", _envelope_cache[key])


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
                   mods: tuple[int, ...] = ()) -> npt.NDArray[np.float64]:
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
            label="D2O",
            num_labeling_sites=spep,
            mods=mods,
        )
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return cast("npt.NDArray[np.float64]", _envelope_cache[key])


def peptide_spep_loss(
    spep_float: float,
    sequence: str,
    pep_mass: float,
    obs_matrix: npt.NDArray[np.float64],        # shape (n_prop, n_iso), normalized per row
    proportions_frac: npt.NDArray[np.float64],  # shape (n_prop,)
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
    observed_iso: npt.ArrayLike,
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
    obs_matrix: npt.NDArray[np.float64],
    proportions_frac: npt.NDArray[np.float64],
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
_adaptive_cache: dict[Any, tuple[float, ...]] = {}


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


def _binned_envelope(dist: IsoSpecPy.Iso, pep_mass: float, n: int) -> tuple[list[float], list[float]]:
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
    init = get_peptide_distribution(sequence, label="D2O", mods=tuple(mods))
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
        sequence, deuterium_enrichment_level=ria_max, label="D2O",
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


def init_channel_masses(
    sequence: str,
    pep_mass: float,
    n: int,
    mods: tuple[int, ...] = (),
) -> tuple[float, ...]:
    """Per-channel **init (unlabeled, θ=0) averaged-isotopolog NEUTRAL masses**.

    The θ=0 reference for the orthogonal mass-defect θ (DeuteRater ΔS / Δmass,
    v1.1.0 item 1). Each entry is the abundance-weighted average exact mass of the
    nominal bin ``round(pep_mass)+iso`` in the natural-abundance envelope — i.e.
    the *centroid* the labelled peak walks away from as deuterium incorporates.

    This is the "init half" of :func:`adaptive_channel_masses` (which documents
    these as the drift-robust θ=0 anchor) but with the channel count ``n`` given
    by the caller — at fit the integrated width is already known, so no Commerford
    final-envelope width union is needed. An empty init bin (pathological) falls
    back to the analytic neutron-spacing comb. Cached in the shared
    :data:`_envelope_cache` (cleared by :func:`clear_envelope_cache`).

    Companion to :func:`_get_init_env` (abundances): same distribution, but the
    averaged masses rather than the per-channel probabilities.
    """
    key = ('init_mass', sequence, mods, int(n))
    cached = _envelope_cache.get(key)
    if cached is not None:
        return cast("tuple[float, ...]", cached)
    dist = get_peptide_distribution(sequence, label="D2O", mods=tuple(mods))
    init_m, _init_p = _binned_envelope(dist, pep_mass, int(n))
    out = tuple(
        m if not math.isnan(m) else pep_mass + iso * 1.003354835
        for iso, m in enumerate(init_m)
    )
    _envelope_cache[key] = out
    return out


def _spacing_components(
    sequence: str, pep_mass: float, spep: int, ria_max: float, n: int,
    mods: tuple[int, ...], label: str,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64],
           npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Cached per-channel (mass, prob) of the init and final envelopes — the pieces
    the mixture-spacing curve is built from. NaN-bin masses fall back to the analytic
    neutron comb. Returns ``(init_m, init_p, final_m, final_p)`` as np arrays."""
    key = ('spacing_comp', sequence, tuple(mods), int(spep),
           round(float(ria_max), 6), int(n), label)
    cached = _envelope_cache.get(key)
    if cached is None:
        init = get_peptide_distribution(sequence, label="D2O", mods=tuple(mods))
        final = get_peptide_distribution(
            sequence, deuterium_enrichment_level=ria_max, label=label,
            num_labeling_sites=spep, mods=tuple(mods),
        )
        im_l, ip_l = _binned_envelope(init, pep_mass, int(n))
        fm_l, fp_l = _binned_envelope(final, pep_mass, int(n))
        im: npt.NDArray[np.float64] = np.array(
            [m if not math.isnan(m) else pep_mass + j * 1.003354835
             for j, m in enumerate(im_l)])
        fm: npt.NDArray[np.float64] = np.array(
            [m if not math.isnan(m) else pep_mass + j * 1.003354835
             for j, m in enumerate(fm_l)])
        ip: npt.NDArray[np.float64] = np.asarray(ip_l, float)
        fp: npt.NDArray[np.float64] = np.asarray(fp_l, float)
        ip = ip / ip.sum() if ip.sum() > 0 else ip
        fp = fp / fp.sum() if fp.sum() > 0 else fp
        cached = (im, ip, fm, fp)
        _envelope_cache[key] = cached
    return cast(
        "tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], "
        "npt.NDArray[np.float64], npt.NDArray[np.float64]]",
        cached,
    )


def _mixture_dspacing(f: float,
                      im: npt.NDArray[np.float64], ip: npt.NDArray[np.float64],
                      fm: npt.NDArray[np.float64], fp: npt.NDArray[np.float64],
                      charge: int, k: int) -> float:
    """Predicted **ΔSₓ(f, k)** in m/z mDa for the init↔final mixture at fraction
    ``f``: the M0-internal spacing of the per-channel intensity-weighted mixture
    centroid, minus its f=0 value. The nonlinear (concave) curve fs_ds inverts."""
    def cz(j: int) -> float:
        w = (1.0 - f) * ip[j] + f * fp[j]
        return float(((1.0 - f) * ip[j] * im[j] + f * fp[j] * fm[j]) / w) if w > 0 else np.nan
    s_f = cz(k) - cz(0)
    s_0 = im[k] - im[0]          # mixture at f=0 is the init envelope
    return float((s_f - s_0) / max(1, int(charge)) * 1e3)


def solve_fs_d2o_ds(
    sequence: str,
    pep_mass: float,
    obs_dspacing: dict[int, float],
    spep: int,
    charge: int,
    ria_max: float = 0.06,
    n_iso: int = _DEFAULT_N_ISO,
    mods: tuple[int, ...] = (),
    label: str = "D2O",
) -> float:
    """Mass-defect fraction-new from the **per-channel Δspacing** — the spacing
    analog of :func:`solve_fs_d2o` (v1.1.0 item 1b).

    Solves ``f`` in ``obs_ΔSₓ(k) ≈ ΔSₓ(f, k)`` by 1-D SSE minimization over the
    given channels, where ``ΔSₓ(f, k)`` is the **nonlinear** init↔final mixture
    spacing curve (:func:`_mixture_dspacing`) — NOT the linear ``ΔSₓ/ΔSₓmax``, which
    over-reads mid-range because the mixture mass-shift is concave in f (Price 2017;
    report 2026-06-25). ``obs_dspacing`` is the **t0/f0-anchored** observed
    M0-internal Δspacing per channel (m/z mDa), keyed by isotopomer index (iso0 is
    ≡0 and excluded by the caller). Needs ≥ 2 channels → else NaN. Bounds match
    :data:`FS_BOUNDS`. ``label`` selects the labeled-envelope chemistry
    (``"D2O"``); the ¹⁸O entry point is :func:`solve_fs_o18_ds` (``label="O18"``).
    """
    from scipy.optimize import minimize_scalar  # local import keeps cold path fast

    ks = [k for k, v in obs_dspacing.items() if v is not None and np.isfinite(v)]
    if len(ks) < 2:
        return float("nan")
    im, ip, fm, fp = _spacing_components(
        sequence, pep_mass, spep, ria_max, int(n_iso), tuple(mods), label)

    def sse(f: float) -> float:
        return float(sum(
            (obs_dspacing[k] - _mixture_dspacing(f, im, ip, fm, fp, charge, k)) ** 2
            for k in ks))

    return float(minimize_scalar(sse, bounds=FS_BOUNDS, method="bounded").x)


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


# ---------------------------------------------------------------------------
# o18 (¹⁸O) Spep / FS solver — the NB90c reverse model (v1.1.0)
# ---------------------------------------------------------------------------

#: o18 length-model feature spec, aligned to the trained coefficient table
#: (``riana/data/coefficients/o18_*.csv``). Spep is
#: ``b·(L-1) + c_D·D + c_E·E + c_N·N + c_Q·Q + c_S·S`` (intercept 0). Backbone is
#: ``L-1`` — the two C-terminal carboxyl oxygens back-exchange with the H₂¹⁶O
#: digest (Previs), leaving the (L-1) internal peptide-bond carbonyls as the
#: stable backbone sites; serine's hydroxyl O is the one labile side-chain site
#: beyond NB90c's D/E/N/Q (threonine/tyrosine were tested and are null).
O18_LENGTH_FEATURES = ("length_minus1", "D", "E", "N", "Q", "S")


def spep_from_length_coefficients(
    sequence: str,
    coefficients: dict[str, float],
) -> float:
    """Per-peptide ¹⁸O Spep from the length-model coefficient table.

    The o18 production analogue of :func:`spep_from_coefficients` (the per-AA
    D₂O model). ``coefficients`` maps each feature in
    :data:`O18_LENGTH_FEATURES` to its learned value. ``sequence`` must be a bare
    AA string (mods stripped) so length and residue counts match the design
    matrix used at training.
    """
    n = len(sequence)
    total = 0.0
    for feat in O18_LENGTH_FEATURES:
        x = (n - 1) if feat == "length_minus1" else sequence.count(feat)
        total += coefficients.get(feat, 0.0) * x
    return float(total)


def _get_o18_final_env(sequence: str, pep_mass: float, spep: int,
                       ria_max: float,
                       n: int = _DEFAULT_N_ISO,
                       mods: tuple[int, ...] = ()) -> npt.NDArray[np.float64]:
    """Fully-labeled ¹⁸O envelope at precursor enrichment ``ria_max`` with
    ``spep`` labile oxygen sites — the o18 analogue of :func:`_get_final_env`.

    The init (natural-abundance) envelope is label-independent, so
    :func:`_get_init_env` is shared with the D₂O path; only the labeled envelope
    differs (¹⁸O 3-isotope O vs D₂O H/D). Cached under a ``'final_o18'`` key so it
    never collides with the D₂O final-envelope cache.
    """
    ria_key = round(float(ria_max), 6)
    key = ('final_o18', sequence, mods, spep, ria_key, n)
    if key not in _envelope_cache:
        dist = get_peptide_distribution(
            sequence,
            deuterium_enrichment_level=ria_max,
            label="O18",
            num_labeling_sites=spep,
            mods=mods,
        )
        env = np.array(get_envelope(dist, pep_mass, n=n + 2))[:n]
        _envelope_cache[key] = env
    return cast("npt.NDArray[np.float64]", _envelope_cache[key])


def solve_fs_o18(
    sequence: str,
    pep_mass: float,
    observed_iso: npt.ArrayLike,
    spep: int,
    ria_max: float = 0.06,
    n_iso: int = _DEFAULT_N_ISO,
    mods: tuple[int, ...] = (),
    score_channels: int | None = None,
) -> float:
    """Per-timepoint fractional synthesis from one observed ¹⁸O envelope.

    The o18 analogue of :func:`solve_fs_d2o`, with the identical H4′
    normalization order (mix full-cluster → truncate to scoring channels →
    renormalize) and widened FS bounds; only the labeled (final) envelope
    differs. See :func:`solve_fs_d2o` for the full rationale. The body is kept
    parallel rather than refactored so the heavily-tested D₂O path stays
    untouched.
    """
    from scipy.optimize import minimize_scalar  # local import keeps cold path fast

    obs = np.asarray(observed_iso, dtype=float)[:n_iso]
    valid = ~np.isnan(obs)
    if not valid.any():
        return float('nan')
    n_real = int(np.max(np.nonzero(valid)[0])) + 1
    obs = obs[:n_real]
    if np.isnan(obs).any():
        return float('nan')
    k = n_real if score_channels is None else min(int(score_channels), n_real)
    if k < 2:
        return float('nan')
    obs_score = obs[:k]
    obs_total = obs_score.sum()
    if obs_total == 0:
        return float('nan')
    obs_norm = obs_score / obs_total

    n_full = max(_FULL_CLUSTER_N, n_real)
    init_full = np.asarray(_get_init_env(sequence, pep_mass, n=n_full, mods=mods), dtype=float)
    final_full = np.asarray(
        _get_o18_final_env(sequence, pep_mass, spep, ria_max, n=n_full, mods=mods), dtype=float,
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


def solve_fs_o18_ds(
    sequence: str,
    pep_mass: float,
    obs_dspacing: dict[int, float],
    spep: int,
    charge: int,
    ria_max: float = 0.06,
    n_iso: int = _DEFAULT_N_ISO,
    mods: tuple[int, ...] = (),
) -> float:
    """¹⁸O mass-defect fraction-new from the per-channel Δspacing — the ¹⁸O analog
    of :func:`solve_fs_d2o_ds`. Mechanically identical (the same nonlinear
    init↔final mixture-spacing inversion over the t0/f0-anchored ``obs_dspacing``),
    but with the **3-isotope enriched-O reverse model** (``label="O18"``: the
    weighted-average mass shift under ¹⁸O masses at ``ria_max`` with ``spep`` labile
    O sites). Thin delegation to :func:`solve_fs_d2o_ds` with ``label="O18"`` so the
    heavily-tested spacing core is shared.

    **Channels.** ¹⁸O is a +2 Da label, so only **iso1** is free of any ¹⁸O-bearing
    isotopolog; **iso2** (one ¹⁸O), **iso3** (one ¹⁸O + one ¹³C, +2+1) and **iso4**
    (two ¹⁸O, or one ¹⁸O + two ¹³C) all carry signal — the caller
    (:func:`riana.core.fitting._fs_ds_points`) scores iso2–4.

    **Low-signal caveat (implemented for completeness/symmetry).** The Δspacing reads
    the mass-defect *difference* between the heavy isotopolog and the ¹³C peak it
    displaces. For ¹⁸O that difference is small — the ¹⁸O isotopolog at iso2 is only
    ~2.5 mDa lighter than 2×¹³C (≈half the per-mass-unit defect of D₂O's D-vs-¹³C) —
    and ¹⁸O carries few labile sites, so the ¹⁸O Δspacing holds far less information
    than D₂O's. This is a symmetry estimator, **not** a reliable second estimate; the
    intensity :func:`solve_fs_o18` stays primary.
    """
    return solve_fs_d2o_ds(
        sequence, pep_mass, obs_dspacing, spep, charge,
        ria_max=ria_max, n_iso=n_iso, mods=mods, label="O18",
    )
