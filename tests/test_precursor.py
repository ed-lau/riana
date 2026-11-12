"""Rising-precursor envelope solve (``riana.algorithms.precursor``).

The forward model is checked against independent constructions: the label-content
closed form against numerical integration, the convolution envelopes against IsoSpec
at each enrichment, the instant-precursor limit against today's plateau mixture, and
the per-point solve against envelopes integrated from IsoSpec over synthesis times.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad

from riana.algorithms.isotope_dist import (
    _get_final_env, _get_init_env, get_envelope, get_peptide_distribution,
)
from riana.algorithms.mass_calc import calculate_ion_mz
from riana.algorithms.precursor import (
    F_UPPER, LABEL_FLOOR, RisingEnvelope, f_lower_bound, label_content, label_gain,
    labeled_envelopes, natural_deuterium, precursor_level, _site_free_envelope,
)

SEQ, SPEP, RIA, KP = "VLSPADKTNVKAAWGK", 28, 0.04, 0.08
N = 24


def _pm() -> float:
    return calculate_ion_mz(SEQ)


@pytest.mark.parametrize("t,k", [(1.0, 0.02), (3.0, 0.3), (10.0, KP), (10.0, KP + 1e-9),
                                 (31.0, 1.0), (79.0, 0.01), (5.0, -0.01)])
def test_label_content_matches_the_integral(t, k):
    num = quad(lambda s: k * math.exp(-k * (t - s)) * (1 - math.exp(-KP * s)), 0, t)[0]
    assert float(label_content(t, k, KP)) == pytest.approx(num, rel=1e-7, abs=1e-12)


def test_label_content_limits():
    t = np.array([0.0, 1.0, 10.0])
    np.testing.assert_allclose(label_content(t, 0.2, math.inf), -np.expm1(-0.2 * t))
    assert float(label_content(0.0, 0.5, KP)) == 0.0
    # Any protein's label content is capped by the precursor reached by then.
    assert float(label_content(5.0, 1e3, KP)) == pytest.approx(
        float(precursor_level(5.0, KP)), rel=1e-3)


def test_label_gain_is_dA_df():
    t, k = 3.0, 0.05
    f = -math.expm1(-k * t)
    h = 1e-6
    k_plus, k_minus = (-math.log1p(-(f + h)) / t, -math.log1p(-(f - h)) / t)
    fd = (float(label_content(t, k_plus, KP)) - float(label_content(t, k_minus, KP))) / (2 * h)
    assert float(label_gain(t, k, KP)) == pytest.approx(fd, rel=1e-5)
    assert float(label_gain(t, k, math.inf)) == pytest.approx(1.0)
    assert float(label_gain(0.0, k, KP)) == 1.0
    # Early in a slow-precursor series the envelope barely moves with f.
    assert float(label_gain(1.0, 0.01, KP)) < 0.05


def test_convolution_envelopes_match_isospec():
    pm = _pm()
    site_free = _site_free_envelope(SEQ, pm, SPEP, N, ())
    for p in (natural_deuterium(), 0.01, RIA):
        conv = labeled_envelopes(site_free, SPEP, np.array([p]))[0]
        iso = np.asarray(get_envelope(get_peptide_distribution(
            SEQ, deuterium_enrichment_level=p, label="D2O", num_labeling_sites=SPEP),
            pm, n=N + 2))[:N]
        np.testing.assert_allclose(conv, iso / iso.sum(), atol=2e-4)


def test_instant_precursor_limit_is_the_plateau_mixture():
    pm = _pm()
    init = _get_init_env(SEQ, pm, n=N)
    fin = _get_final_env(SEQ, pm, SPEP, RIA, n=N)
    env = RisingEnvelope(SEQ, pm, SPEP, RIA, 10.0, 1e6)
    for f in (0.1, 0.6):
        mix = (1 - f) * init / init.sum() + f * fin / fin.sum()
        np.testing.assert_allclose(env.predict(f), mix, atol=2e-4)


def test_mean_mass_shift_is_label_content():
    pm = _pm()
    init = _get_init_env(SEQ, pm, n=N)
    idx = np.arange(N)
    for t, k in ((1.0, 0.02), (10.0, 0.05), (40.0, 0.3)):
        env = RisingEnvelope(SEQ, pm, SPEP, RIA, t, KP)
        shift = env.predict(-math.expm1(-k * t)) @ idx - (init / init.sum()) @ idx
        read = shift / (SPEP * (RIA - natural_deuterium()))
        assert read == pytest.approx(float(label_content(t, k, KP)), rel=3e-3, abs=1e-5)


def _isospec_integrated(t: float, k: float, ns: int = 200) -> np.ndarray:
    """Independent forward model: IsoSpec envelope per synthesis interval."""
    pm = _pm()
    init = _get_init_env(SEQ, pm, n=N)
    out = math.exp(-k * t) * init / init.sum()
    edges = np.linspace(0, t, ns + 1)
    w = np.exp(-k * (t - edges[1:])) - np.exp(-k * (t - edges[:-1]))
    nat = natural_deuterium()
    for wj, s in zip(w, 0.5 * (edges[:-1] + edges[1:])):
        p = nat + (RIA - nat) * float(precursor_level(s, KP))
        e = np.asarray(get_envelope(get_peptide_distribution(
            SEQ, deuterium_enrichment_level=p, label="D2O", num_labeling_sites=SPEP),
            pm, n=N + 2))[:N]
        out = out + wj * e / e.sum()
    return out


@pytest.mark.parametrize("t,k", [(1.0, 0.02), (5.0, 0.05), (20.0, 0.05), (40.0, 0.01)])
def test_solve_recovers_f_from_an_independent_forward_model(t, k):
    obs = _isospec_integrated(t, k)[:6] * 1e6
    f, a = RisingEnvelope(SEQ, _pm(), SPEP, RIA, t, KP).solve(obs, n_iso=6)
    assert f == pytest.approx(-math.expm1(-k * t), abs=1e-3)
    assert a == pytest.approx(float(label_content(t, k, KP)), abs=2e-4)


def test_below_natural_reads_solve_on_the_label_scale():
    """Early, a read slightly lighter than natural is a small negative label content
    (f far below −0.05, since the gain is tiny) — not clipped at an f bound."""
    pm = _pm()
    env = RisingEnvelope(SEQ, pm, SPEP, RIA, 1.0, KP)
    f_lo = f_lower_bound(1.0, KP)
    k_lo = -math.log1p(-f_lo) / 1.0
    assert float(label_content(1.0, k_lo, KP)) == pytest.approx(LABEL_FLOOR, abs=1e-9)
    from scipy.optimize import brentq
    k_neg = brentq(lambda k: float(label_content(1.0, k, KP)) + 0.03, -50.0, 0.0)
    f_neg = -math.expm1(-k_neg)
    f, a = env.solve(env.predict(f_neg)[:6] * 1e6, n_iso=6)
    assert f < -0.5
    assert a == pytest.approx(-0.03, abs=1e-4)


def test_read_above_the_ceiling_solves_to_the_bound():
    pm = _pm()
    env = RisingEnvelope(SEQ, pm, SPEP, RIA, 1.0, KP)
    fin = _get_final_env(SEQ, pm, SPEP, RIA, n=N)
    init = _get_init_env(SEQ, pm, n=N)
    heavy = 0.8 * init / init.sum() + 0.2 * fin / fin.sum()      # A = 0.2 > P(1) ≈ 0.07
    f, _ = env.solve(heavy[:6] * 1e6, n_iso=6)
    assert f == pytest.approx(F_UPPER, abs=1e-4)


def test_solve_observed_channel_handling():
    env = RisingEnvelope(SEQ, _pm(), SPEP, RIA, 3.0, KP)
    assert all(math.isnan(x) for x in env.solve(np.full(6, np.nan), n_iso=6))
    assert all(math.isnan(x) for x in env.solve(np.array([1.0, np.nan, 1.0]), n_iso=3))
    assert all(math.isnan(x) for x in env.solve(np.array([1.0, 0.5]), n_iso=2,
                                                score_channels=1))
    # Trailing NaN padding (adaptive N_ISO) is ignored, as in solve_fs_d2o.
    obs = env.predict(0.2)[:6] * 1e6
    padded = np.concatenate([obs[:4], [np.nan, np.nan]])
    f, _ = env.solve(padded, n_iso=6)
    assert f == pytest.approx(0.2, abs=1e-3)


def test_label_content_is_finite_for_a_fast_precursor():
    """k_p·t large (a primed or near-instant precursor): no 0·∞ — A → 1 − e^{−kt}."""
    t = np.array([1.0, 84.0])
    for k in (-0.1, 0.01, 0.3):
        a = label_content(t, k, 100.0)
        assert np.isfinite(a).all()
        np.testing.assert_allclose(a, -np.expm1(-k * t), rtol=0.02, atol=1e-3)
    assert math.isfinite(f_lower_bound(84.0, 100.0))
    assert float(label_content(10.0, KP + 2e-4, KP)) == pytest.approx(
        float(label_content(10.0, KP + 2e-3, KP)), rel=0.05)     # smooth across the branch switch
