"""Unit tests for Phase B algorithm modules.

Each module under :mod:`riana.algorithms` is a set of pure numpy functions
that Phase C composes into the integrator. These tests verify each piece
in isolation on synthetic inputs (gaussian peaks + known baselines + known
mass shifts) so the composition stage in Phase C has a stable floor.
"""

from __future__ import annotations

import numpy as np
import pytest

from riana.algorithms import baseline as ba
from riana.algorithms import calibration as cal
from riana.algorithms import isotope_dist as iso
from riana.algorithms import peaks as pk
from riana.algorithms import smoothing as sm
from riana.algorithms.mass_calc import calculate_ion_mz


# ---------------------------------------------------------------------------
# Synthetic helpers — used by multiple test groups
# ---------------------------------------------------------------------------


def _gaussian(rt: np.ndarray, center: float, sigma: float, height: float) -> np.ndarray:
    return height * np.exp(-0.5 * ((rt - center) / sigma) ** 2)


# ---------------------------------------------------------------------------
# smoothing
# ---------------------------------------------------------------------------


def test_savgol_polyorder_1_rejected():
    """The §2c fix: polyorder=1 (moving average) is no longer permitted."""
    with pytest.raises(ValueError, match="polyorder must be"):
        sm.savgol(np.arange(20.0), window=5, polyorder=1)


def test_savgol_preserves_gaussian_height():
    rt = np.linspace(0, 1, 101)
    sig = _gaussian(rt, 0.5, 0.05, 1000.0)
    out = sm.savgol(sig, window=7, polyorder=2)
    # polyorder=2 SG with a narrow window keeps peak height within a few %.
    assert abs(out.max() - sig.max()) / sig.max() < 0.05


def test_savgol_short_or_zero_trace_passthrough():
    assert np.array_equal(sm.savgol(np.zeros(5), window=3), np.zeros(5))
    assert np.array_equal(sm.savgol(np.array([1.0, 2.0]), window=5), np.array([1.0, 2.0]))


def test_savgol_validates_window():
    with pytest.raises(ValueError, match="odd"):
        sm.savgol(np.zeros(10), window=4, polyorder=2)
    with pytest.raises(ValueError, match="window"):
        sm.savgol(np.zeros(10), window=3, polyorder=3)


# ---------------------------------------------------------------------------
# baseline
# ---------------------------------------------------------------------------


def test_noise_floor_recovers_constant_noise_under_peak():
    rng = np.random.default_rng(0)
    rt = np.linspace(0, 1, 201)
    noise = 50.0 + 5.0 * rng.standard_normal(rt.size)  # mean 50, σ 5
    sig = noise + _gaussian(rt, 0.5, 0.03, 800.0)
    bl = ba.noise_floor(sig)
    # p10 of (noise + occasional peak) is approximately noise's p10.
    # noise mean=50, σ=5 ⇒ p10 ≈ 50 - 1.28*5 ≈ 43.6.
    assert bl[0] == pytest.approx(bl[-1])  # flat
    assert 35 < bl[0] < 55


def test_noise_floor_handles_zero_or_empty_trace():
    assert ba.noise_floor(np.zeros(10))[0] == 0.0
    assert ba.noise_floor(np.array([], dtype=np.float64)).size == 0


def test_noise_floor_robust_to_endpoint_spike():
    """The whole point of noise_floor vs local_linear: an endpoint spike
    that would inflate local_linear's interpolation does not move the
    p10 noise estimate."""
    rt = np.linspace(0, 1, 51)
    sig = 10.0 + _gaussian(rt, 0.5, 0.04, 100.0)
    spike = sig.copy()
    spike[0] = 500.0  # noise spike at endpoint
    bl_clean = ba.noise_floor(sig)
    bl_spike = ba.noise_floor(spike)
    # p10 is unaffected by a single high spike.
    assert abs(bl_spike[0] - bl_clean[0]) < 1.0


def test_local_linear_subtracts_known_linear_drift():
    rt = np.linspace(0, 1, 51)
    drift = 100.0 + 50.0 * rt
    peak = _gaussian(rt, 0.5, 0.05, 500.0)
    sig = drift + peak
    bl = ba.local_linear(rt, sig, lo=10, hi=40)
    # Inside the peak window, baseline matches the linear drift's endpoints
    # (peak-free at the boundaries).
    np.testing.assert_allclose(bl[[10, 40]], sig[[10, 40]], atol=1e-9)
    # And subtracting it leaves a near-pure gaussian (centre intensity
    # within a small fraction of the original peak height).
    centred = sig - bl
    assert centred.max() / 500.0 > 0.95


def test_local_linear_flat_extends_outside():
    rt = np.linspace(0, 1, 21)
    sig = np.full(21, 7.0)
    bl = ba.local_linear(rt, sig, lo=5, hi=15)
    assert np.all(bl[:5] == sig[5])
    assert np.all(bl[16:] == sig[15])


def test_local_linear_rejects_bad_indices():
    with pytest.raises(ValueError):
        ba.local_linear(np.arange(5.0), np.arange(5.0), lo=3, hi=1)


def test_snip_baseline_below_signal():
    rt = np.linspace(0, 1, 201)
    sig = 50.0 * np.sin(np.pi * rt) + _gaussian(rt, 0.5, 0.04, 300.0)
    bl = ba.snip(sig, n_iter=40)
    # SNIP should track the slow hump (~50) without consuming the peak.
    assert bl.max() < sig.max()
    assert (bl <= sig + 1e-6).all()


def test_asls_returns_smooth_baseline():
    rt = np.linspace(0, 1, 201)
    sig = 30.0 + 200.0 * _gaussian(rt, 0.5, 0.03, 1.0)  # broad hump + sharp peak
    bl = ba.asls(sig, lam=1e5, p=0.01)
    # AsLS baseline stays below the peak apex and is smooth (small variance
    # in first differences relative to the signal).
    assert bl.max() < sig.max()
    sig_diff_var = np.var(np.diff(sig))
    bl_diff_var = np.var(np.diff(bl))
    assert bl_diff_var < sig_diff_var * 0.1


# ---------------------------------------------------------------------------
# peaks
# ---------------------------------------------------------------------------


def test_detect_peak_finds_apex_near_psm_prior():
    rt = np.linspace(0, 1, 201)
    sig = _gaussian(rt, 0.6, 0.03, 1000.0) + 10.0 * np.random.default_rng(0).standard_normal(rt.size)
    result = pk.detect_peak(rt, sig, scan_prior_rt=0.6)
    assert result is not None
    assert abs(rt[result.apex_idx] - 0.6) < 0.02
    assert result.lo < result.apex_idx < result.hi


def test_detect_peak_picks_nearest_when_multiple():
    """With two peaks, the one closer to the PSM-scan prior wins."""
    rt = np.linspace(0, 1, 201)
    sig = _gaussian(rt, 0.3, 0.02, 500.0) + _gaussian(rt, 0.7, 0.02, 500.0)
    # Same-height peaks: nearer prior wins.
    res_left = pk.detect_peak(rt, sig, scan_prior_rt=0.25)
    res_right = pk.detect_peak(rt, sig, scan_prior_rt=0.75)
    assert rt[res_left.apex_idx] < 0.5 < rt[res_right.apex_idx]


def test_detect_peak_returns_none_on_flat_or_zero():
    rt = np.linspace(0, 1, 50)
    assert pk.detect_peak(rt, np.zeros(50), scan_prior_rt=0.5) is None
    # Pure noise should not pass the prominence floor.
    rng = np.random.default_rng(0)
    noise = 10.0 + rng.standard_normal(50)
    assert pk.detect_peak(rt, noise, scan_prior_rt=0.5) is None


def test_coelution_ok_passes_aligned_apices():
    rt = np.linspace(0, 1, 201)
    iso0 = pk.detect_peak(rt, _gaussian(rt, 0.5, 0.02, 1000.0), scan_prior_rt=0.5)
    # iso1 with apex one MS1 cycle later — should pass.
    iso1 = pk.detect_peak(rt, _gaussian(rt, 0.505, 0.02, 600.0), scan_prior_rt=0.5)
    assert pk.coelution_ok(iso0, iso1, rt) is True


def test_coelution_ok_rejects_distant_apices():
    rt = np.linspace(0, 1, 201)
    iso0 = pk.detect_peak(rt, _gaussian(rt, 0.3, 0.02, 1000.0), scan_prior_rt=0.3)
    iso1 = pk.detect_peak(rt, _gaussian(rt, 0.7, 0.02, 600.0), scan_prior_rt=0.7)
    assert pk.coelution_ok(iso0, iso1, rt) is False


def test_coelution_ok_tolerates_missing_iso1():
    rt = np.linspace(0, 1, 201)
    iso0 = pk.detect_peak(rt, _gaussian(rt, 0.5, 0.02, 1000.0), scan_prior_rt=0.5)
    # iso1=None must NOT reject — genuine low-iso1 envelopes at high D2O
    # are legitimate and should keep the iso0 integration.
    assert pk.coelution_ok(iso0, None, rt) is True


def test_snr_reflects_baseline_quality():
    rt = np.linspace(0, 1, 201)
    sig = _gaussian(rt, 0.5, 0.03, 1000.0) + 10.0
    bl = np.full_like(sig, 10.0)
    # Clean peak above flat baseline → S/N is large.
    assert pk.snr(sig, bl) > 100


def test_symmetry_returns_one_for_gaussian():
    rt = np.linspace(0, 1, 201)
    sig = _gaussian(rt, 0.5, 0.05, 1000.0)
    apex = int(sig.argmax())
    assert pk.symmetry(rt, sig, apex) > 0.95


def test_quality_blends_snr_and_symmetry():
    # High S/N + perfect symmetry ⇒ near 1.0.
    assert pk.quality(snr_value=20.0, symmetry_value=1.0) == pytest.approx(1.0, abs=0.05)
    # Zero on either side ⇒ 0.
    assert pk.quality(snr_value=0.0, symmetry_value=1.0) == 0.0
    assert pk.quality(snr_value=20.0, symmetry_value=0.0) == 0.0


# ---------------------------------------------------------------------------
# calibration
# ---------------------------------------------------------------------------


def test_weighted_obs_mz_returns_centroid_under_tol():
    rng = np.random.default_rng(0)
    target = 500.0
    # Three centroids inside ±15 ppm of 500.0 (delta ≈ 0.0075 Da).
    mzs = np.array([499.9999, 500.0010, 500.0030])
    ints = np.array([1.0, 4.0, 1.0])
    obs, ppm = cal.weighted_obs_mz(mzs, ints, target_mz=target, half_width_ppm=15.0)
    # Intensity-weighted centroid lands near the dominant centroid (500.0010).
    assert abs(obs - 500.0010) < 1e-3
    assert abs(ppm - (obs - target) / target * 1e6) < 1e-3


def test_weighted_obs_mz_returns_none_when_no_centroid_in_window():
    obs, ppm = cal.weighted_obs_mz(
        np.array([400.0, 600.0]), np.array([1.0, 1.0]),
        target_mz=500.0, half_width_ppm=10.0,
    )
    assert obs is None and ppm is None


def test_drift_summary_aggregates_ppm_errors():
    ppm = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
    s = cal.drift_summary(ppm)
    assert s.n == 5
    assert s.median_ppm == 4.0
    assert s.suggested_shift_ppm == -4.0
    # MAD of [2,3,4,5,6] is 1; 1.4826 * 1 ≈ 1.4826.
    assert abs(s.mad_ppm - 1.4826) < 1e-6


def test_drift_summary_ignores_nan_entries():
    ppm = np.array([2.0, np.nan, 4.0, np.nan, 6.0])
    s = cal.drift_summary(ppm)
    assert s.n == 3
    assert s.median_ppm == 4.0


def test_drift_summary_empty_safe():
    s = cal.drift_summary(np.array([], dtype=np.float64))
    assert s.n == 0 and s.suggested_shift_ppm == 0.0


# ---------------------------------------------------------------------------
# isotope_dist — Spep / FS solver (M3 Week 4 lift)
# ---------------------------------------------------------------------------


# Fixed test peptide. ~9 labile H sites (Ala/Val/Leu/Asp/Ile/Glu/Lys/Arg
# all carry α-H; this is the M2 default test peptide pattern).
_TEST_SEQ = "VAPEPTIDEK"


def _peptide_mass(seq: str) -> float:
    """Neutral monoisotopic mass via the production mass_calc."""
    return calculate_ion_mz(seq)


def test_get_envelope_returns_n_bins_summing_to_distribution_prob():
    pep_mass = _peptide_mass(_TEST_SEQ)
    dist = iso.get_peptide_distribution(_TEST_SEQ, label=1)
    env = iso.get_envelope(dist, pep_mass, n=4)
    assert len(env) == 4
    # Natural-abundance envelope is monotone-decreasing past iso0.
    assert env[0] > env[1] > env[2]
    # The 4 captured bins should cover most of the (prob_to_cover=0.999) mass.
    assert sum(env) > 0.95


def test_envelope_cache_idempotent_under_clear():
    iso.clear_envelope_cache()
    pep_mass = _peptide_mass(_TEST_SEQ)
    env1 = iso._get_init_env(_TEST_SEQ, pep_mass, n=4)
    env2 = iso._get_init_env(_TEST_SEQ, pep_mass, n=4)
    # Same key returns the SAME cached array (identity check) until cleared.
    assert env1 is env2
    iso.clear_envelope_cache()
    env3 = iso._get_init_env(_TEST_SEQ, pep_mass, n=4)
    # After clear, a fresh array is returned (still equal content).
    assert env3 is not env1
    np.testing.assert_allclose(env3, env1)


def test_peptide_spep_loss_minimized_at_true_spep():
    """Synthetic: build observed envelopes from a known Spep, then verify
    the SSE loss has its minimum near the true Spep."""
    iso.clear_envelope_cache()
    pep_mass = _peptide_mass(_TEST_SEQ)
    proportions = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    true_spep = 8

    init = iso._get_init_env(_TEST_SEQ, pep_mass, n=4)
    final = iso._get_final_env(_TEST_SEQ, pep_mass, true_spep, ria_max=0.06, n=4)
    init_norm = init / init.sum()
    final_norm = final / final.sum()
    obs_matrix = np.stack([
        (1 - f) * init_norm + f * final_norm for f in proportions
    ])

    # Loss at true Spep is essentially zero; at any other integer it grows.
    loss_true = iso.peptide_spep_loss(
        true_spep, _TEST_SEQ, pep_mass, obs_matrix, proportions,
    )
    loss_off_lo = iso.peptide_spep_loss(
        true_spep - 3, _TEST_SEQ, pep_mass, obs_matrix, proportions,
    )
    loss_off_hi = iso.peptide_spep_loss(
        true_spep + 3, _TEST_SEQ, pep_mass, obs_matrix, proportions,
    )
    assert loss_true < 1e-6
    assert loss_off_lo > loss_true
    assert loss_off_hi > loss_true


def test_fit_peptide_spep_recovers_synthetic_spep():
    """Brent optimization recovers the true Spep within 0.5 sites."""
    iso.clear_envelope_cache()
    pep_mass = _peptide_mass(_TEST_SEQ)
    proportions = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    true_spep = 8
    init = iso._get_init_env(_TEST_SEQ, pep_mass, n=4)
    final = iso._get_final_env(_TEST_SEQ, pep_mass, true_spep, ria_max=0.06, n=4)
    obs = np.stack([
        (1 - f) * init / init.sum() + f * final / final.sum() for f in proportions
    ])

    fitted = iso.fit_peptide_spep(_TEST_SEQ, pep_mass, obs, proportions)
    assert abs(fitted - true_spep) < 0.5, f"fitted Spep {fitted} far from true {true_spep}"


def test_solve_fs_d2o_recovers_known_fractions():
    """Per-timepoint FS solver recovers each known proportion within ±0.02."""
    iso.clear_envelope_cache()
    pep_mass = _peptide_mass(_TEST_SEQ)
    spep = 8
    init = iso._get_init_env(_TEST_SEQ, pep_mass, n=4)
    final = iso._get_final_env(_TEST_SEQ, pep_mass, spep, ria_max=0.06, n=4)
    init_norm = init / init.sum()
    final_norm = final / final.sum()

    for true_fs in (0.0, 0.125, 0.5, 0.75, 1.0):
        obs = (1 - true_fs) * init_norm + true_fs * final_norm
        # The solver expects un-normalized intensities (it normalizes
        # internally); scale by a synthetic peak area to test that path.
        obs_unnorm = obs * 1e6
        recovered = iso.solve_fs_d2o(_TEST_SEQ, pep_mass, obs_unnorm, spep, n_iso=4)
        assert abs(recovered - true_fs) < 0.02, (
            f"FS recovery {recovered:.3f} far from true {true_fs:.3f}"
        )


def test_solve_fs_d2o_returns_nan_on_zero_intensity():
    iso.clear_envelope_cache()
    pep_mass = _peptide_mass(_TEST_SEQ)
    result = iso.solve_fs_d2o(
        _TEST_SEQ, pep_mass, observed_iso=[0.0, 0.0, 0.0, 0.0],
        spep=8, n_iso=4,
    )
    assert np.isnan(result)
