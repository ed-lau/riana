"""Rising precursor: the per-point envelope solve for the fraction new under a slowly
rising D₂O precursor (PROJECT_REVIEW "Slow precursor" box; report
``reports/2026-11-11_non_first_order_linearization_feasibility.md``).

Without a priming dose the body-water enrichment rises as ``P(s) = 1 − e^{−k_p·s}``
(fraction of plateau). A molecule made at time ``s`` carries the enrichment of that time,
so at labeling time ``t`` the peptide pool is the old (natural) molecules plus new ones
spread over synthesis times:

    E_pred(t; k) = e^{−kt}·E(natural) + ∫₀ᵗ k·e^{−k(t−s)}·E(p(s)) ds,
    p(s) = p_nat + (RIA_max − p_nat)·P(s)

The fraction made since labeling began, ``f = 1 − e^{−kt}``, is still first-order under
any shared precursor curve, so each point is solved for the ``k_i`` (equivalently
``f_i``) whose predicted envelope matches the observed one, and everything downstream
(linear Δk, interaction) runs on ``f``. With ``k_p → ∞`` the integrand is the plateau
envelope throughout and the solve reduces to today's θ.

**Label content** ``A`` — the pool's mean enrichment ÷ the plateau enrichment, to first
order the envelope's mean mass shift whatever the precursor model — is the measured
scale: ``A(t; k) = ∫₀ᵗ k·e^{−k(t−s)}·P(s) ds`` (the Guan curve), ``≤ P(t)`` for any k.
The fit's rails and its peptide curve fit (R², ``k_cv``) work on ``A``, not on ``f``:
early in a series ``A`` spans only ``[0, P(t)]``, so the envelope → ``f`` map is nearly
flat there and noise in the envelope becomes large swings in ``f``.

**Envelopes by convolution.** The labeled envelope at enrichment ``p`` is the
peptide's envelope with its labeling-site hydrogens fixed at ¹H (``E0′``, IsoSpec, once
per peptidoform) convolved with ``Binomial(spep, p)`` in isotopomer steps. That matches
IsoSpec within ~10⁻⁴ and needs no IsoSpec call per enrichment, which the synthesis-time
integral would otherwise make per grid point. D₂O only (``--label hw``).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln

from riana.algorithms.isotope_dist import (
    FS_BOUNDS, _FULL_CLUSTER_N, _get_init_env, get_envelope, get_peptide_distribution,
)

FloatArray = npt.NDArray[np.float64]

#: Synthesis-time intervals over ``[0, t]`` for the integral (exact interval weights,
#: midpoint enrichment). 100 keeps the k bias of the simulation benches ≤ 1 %.
N_SYNTH = 100
#: Upper bound of the per-point solve in ``f``. A read at or above the precursor
#: ceiling ``P(t)`` (no ``f`` < 1 reproduces it) solves to this bound.
F_UPPER = 0.999
#: Solves at or above this ``f`` are failed solves (the read is at/above the ceiling);
#: the fit's rails drop them. Beyond the φ plateau truncation (−4 ≈ f 0.98) anyway.
F_FAILED = 0.99
#: Lower bound of the per-point solve, in **label content**: the ``f`` bound is set per
#: timepoint so that ``A`` can read down to here (at early t that is ``f`` ≪ 0 — the
#: gain is tiny). The fit's below-natural rail (−0.05 on ``A``) sits inside it.
LABEL_FLOOR = -0.1

_SITE_FREE_CACHE: dict[Any, FloatArray] = {}
_F_LOWER_CACHE: dict[tuple[float, float], float] = {}
_NATURAL_D: list[float] = []


def clear_precursor_cache() -> None:
    """Drop the memoized site-free envelopes and solve bounds (tests)."""
    _SITE_FREE_CACHE.clear()
    _F_LOWER_CACHE.clear()


def natural_deuterium() -> float:
    """IsoSpec's natural ²H abundance (what the natural-abundance envelope uses)."""
    if not _NATURAL_D:
        import IsoSpecPy

        _, _, probs, _ = IsoSpecPy.IsoParamsFromDict(formula={"H": 1})
        _NATURAL_D.append(float(probs[0][1]))
    return _NATURAL_D[0]


def precursor_level(s: npt.ArrayLike, k_p: float) -> FloatArray:
    """``P(s) = 1 − e^{−k_p·s}``, the precursor enrichment as a fraction of plateau."""
    arr: FloatArray = -np.expm1(-float(k_p) * np.asarray(s, dtype=np.float64))
    return arr


def label_content(t: npt.ArrayLike, k: npt.ArrayLike, k_p: float) -> FloatArray:
    """Label content ``A(t; k) = ∫₀ᵗ k·e^{−k(t−s)}·(1 − e^{−k_p·s}) ds`` — the Guan curve
    ``1 − (k_p·e^{−kt} − k·e^{−k_p·t})/(k_p − k)``, written pole-free at ``k = k_p``.

    Vectorized over ``t`` and ``k``; ``A = 0`` at ``t ≤ 0``; ``k < 0`` continues the
    formula smoothly (a below-natural read); ``k_p = ∞`` gives ``1 − e^{−kt}``.
    """
    t_a, k_a = np.broadcast_arrays(np.asarray(t, dtype=np.float64),
                                   np.asarray(k, dtype=np.float64))
    f = -np.expm1(-k_a * t_a)
    if math.isinf(k_p):
        out_inf: FloatArray = np.where(t_a > 0, f, 0.0)
        return out_inf
    x = (k_a - k_p) * t_a
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        # Near the pole (|x| small) the series-safe form; elsewhere the direct form, which
        # stays finite when k_p·t is large (e^{−k_p·t} underflows to 0 there, and the
        # series form would multiply that 0 by an overflowed e^{|x|}).
        near = k_a * t_a * np.exp(-k_p * t_a) * np.where(
            np.abs(x) < 1e-8, 1.0 - x / 2.0, -np.expm1(-x) / np.where(x == 0, 1.0, x))
        far = k_a * (np.exp(-k_p * t_a) - np.exp(-k_a * t_a)) / np.where(
            k_a == k_p, 1.0, k_a - k_p)
        lag = np.where(np.abs(x) < 1e-3, near, far)
    out: FloatArray = np.where(t_a > 0, f - lag, 0.0)
    return out


def label_gain(t: npt.ArrayLike, k: npt.ArrayLike, k_p: float) -> FloatArray:
    """``∂A/∂f`` at ``(t, k)``: how far label content moves per unit of ``f`` at that
    time — ``(∂A/∂k)/(t·e^{−kt})``. 1 under an instant precursor and at ``t = 0``; small
    early in a slow-precursor series. ``Var(f) = Var(A)/gain²``."""
    t_a, k_a = np.broadcast_arrays(np.asarray(t, dtype=np.float64),
                                   np.asarray(k, dtype=np.float64))
    h = 1e-6 * np.maximum(np.abs(k_a), 1e-3)
    dA = (label_content(t_a, k_a + h, k_p) - label_content(t_a, k_a - h, k_p)) / (2.0 * h)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        g = dA / (t_a * np.exp(-k_a * t_a))
    out: FloatArray = np.where(t_a > 0, g, 1.0)
    return out


def f_lower_bound(t: float, k_p: float) -> float:
    """The ``f`` at which label content reads :data:`LABEL_FLOOR` at time ``t`` — the
    per-timepoint lower bound of the solve. ``FS_BOUNDS[0]`` at ``t ≤ 0``."""
    if t <= 0:
        return float(FS_BOUNDS[0])
    key = (float(t), float(k_p))
    hit = _F_LOWER_CACHE.get(key)
    if hit is not None:
        return hit
    lo = -0.1
    while float(label_content(t, lo, k_p)) > LABEL_FLOOR and lo > -1e6:
        lo *= 2.0   # bracket by doubling (label content falls with k, overflow-safe)
    k_lo = float(brentq(lambda k: float(label_content(t, k, k_p)) - LABEL_FLOOR, lo, 0.0))
    val = float(-math.expm1(-k_lo * t))
    _F_LOWER_CACHE[key] = val
    return val


def _binomial_rows(n: int, p: FloatArray) -> FloatArray:
    """Binomial(n, p) pmf over 0..n, one row per enrichment in ``p``."""
    j = np.arange(n + 1, dtype=np.float64)
    pc = np.clip(p, 1e-300, 1.0 - 1e-16)[:, None]
    log_c = gammaln(n + 1.0) - gammaln(j + 1.0) - gammaln(n - j + 1.0)
    rows: FloatArray = np.exp(log_c + j * np.log(pc) + (n - j) * np.log1p(-pc))
    return rows


def _site_free_envelope(sequence: str, pep_mass: float, spep: int, n_full: int,
                        mods: tuple[int, ...]) -> FloatArray:
    """``E0′``: the envelope with the ``spep`` labeling-site hydrogens fixed at ¹H,
    normalized over ``n_full`` channels (cached per peptidoform)."""
    key = (sequence, tuple(mods), round(float(pep_mass), 6), int(spep), int(n_full))
    hit = _SITE_FREE_CACHE.get(key)
    if hit is not None:
        return hit
    dist = get_peptide_distribution(sequence, deuterium_enrichment_level=1e-15,
                                    label="D2O", num_labeling_sites=int(spep),
                                    mods=tuple(mods))
    env = np.asarray(get_envelope(dist, pep_mass, n=n_full + 2), dtype=np.float64)[:n_full]
    total = float(env.sum())
    env = env / total if total > 0 else env
    _SITE_FREE_CACHE[key] = env
    return env


def labeled_envelopes(site_free: FloatArray, spep: int, p: FloatArray) -> FloatArray:
    """Envelopes at each site enrichment in ``p`` (rows), by convolution of ``E0′`` with
    ``Binomial(spep, p)``, truncated and renormalized to ``len(site_free)`` channels."""
    n_full = len(site_free)
    binom = _binomial_rows(int(spep), np.atleast_1d(np.asarray(p, dtype=np.float64)))
    out = np.zeros((binom.shape[0], n_full))
    for j in range(min(binom.shape[1], n_full)):
        out[:, j:] += binom[:, [j]] * site_free[: n_full - j]
    sums = out.sum(axis=1, keepdims=True)
    rows: FloatArray = out / np.where(sums > 0, sums, 1.0)
    return rows


class RisingEnvelope:
    """The forward model at one labeling time for one peptidoform: the predicted
    full-cluster envelope as a function of ``f``, and the per-point solve.

    Built per (peptidoform, t) inside one peptide's fit and discarded after — the
    synthesis-time basis is ``N_SYNTH × n_full`` floats, too many to cache run-wide.
    """

    def __init__(self, sequence: str, pep_mass: float, spep: int, ria_max: float,
                 t: float, k_p: float, *, mods: Sequence[int] = (),
                 n_full: int = _FULL_CLUSTER_N) -> None:
        self.t = float(t)
        self.k_p = float(k_p)
        mods_t = tuple(int(m) for m in mods)
        init = np.asarray(_get_init_env(sequence, pep_mass, n=n_full, mods=mods_t),
                          dtype=np.float64)
        self.natural: FloatArray = init / init.sum()
        site_free = _site_free_envelope(sequence, pep_mass, spep, n_full, mods_t)
        p_nat = natural_deuterium()
        if self.t > 0:
            self.edges: FloatArray = np.linspace(0.0, self.t, N_SYNTH + 1)
            mid = 0.5 * (self.edges[:-1] + self.edges[1:])
            p = p_nat + (float(ria_max) - p_nat) * precursor_level(mid, self.k_p)
            self.basis: FloatArray = labeled_envelopes(site_free, spep, p)
        else:
            # No history at t = 0: the plateau mixture (today's θ), the only defined solve.
            self.plateau: FloatArray = labeled_envelopes(site_free, spep,
                                                         np.array([float(ria_max)]))[0]

    def predict(self, f: float) -> FloatArray:
        """Predicted full-cluster envelope (sums to 1) at fraction new ``f``."""
        if self.t > 0:
            f = min(f, 1.0 - 1e-12)   # f = 1 (k → ∞): every molecule made at the end
            k = -math.log1p(-f) / self.t
            w = (np.exp(-k * (self.t - self.edges[1:]))
                 - np.exp(-k * (self.t - self.edges[:-1])))
            pred: FloatArray = (1.0 - f) * self.natural + w @ self.basis
            return pred
        mix: FloatArray = (1.0 - f) * self.natural + f * self.plateau
        return mix

    def label_content_of(self, f: float) -> float:
        """Label content of the pool at fraction new ``f`` (``= f`` at ``t = 0``)."""
        if self.t <= 0:
            return float(f)
        k = -math.log1p(-f) / self.t
        return float(label_content(self.t, k, self.k_p))

    def solve(self, observed_iso: npt.ArrayLike, n_iso: int,
              score_channels: int | None = None) -> tuple[float, float]:
        """``(f, label content)`` of one observed envelope; ``(nan, nan)`` when it cannot
        be solved. Observed-channel handling and the H4′ order (mix in the full-cluster
        basis → truncate to the scored channels → renormalize) as
        :func:`riana.algorithms.isotope_dist.solve_fs_d2o`."""
        obs = np.asarray(observed_iso, dtype=np.float64)[:n_iso]
        valid = ~np.isnan(obs)
        if not valid.any():
            return float("nan"), float("nan")
        n_real = int(np.max(np.nonzero(valid)[0])) + 1
        obs = obs[:n_real]
        if np.isnan(obs).any():
            return float("nan"), float("nan")
        kk = n_real if score_channels is None else min(int(score_channels), n_real)
        if kk < 2:
            return float("nan"), float("nan")
        total = float(obs[:kk].sum())
        if total == 0:
            return float("nan"), float("nan")
        obs_norm = obs[:kk] / total

        def sse(f: float) -> float:
            pred = self.predict(f)[:kk]
            ps = float(pred.sum())
            if not np.isfinite(ps) or ps <= 0:
                return 1e6
            return float(np.sum((obs_norm - pred / ps) ** 2))

        if self.t > 0:
            bounds = (f_lower_bound(self.t, self.k_p), F_UPPER)
        else:
            bounds = (float(FS_BOUNDS[0]), float(FS_BOUNDS[1]))
        f = float(minimize_scalar(sse, bounds=bounds, method="bounded").x)
        return f, self.label_content_of(f)
