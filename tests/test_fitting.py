"""Synthetic-data tests for ``core/fitting.py`` (M3 Week 4 Phase F2).

Builds per-timepoint observed envelopes from a known (k_deg, Spep) pair via
the IsoSpec forward model, runs the new ``fit_run`` with a synthetic
coefficient table calibrated to give the right Spep for our test peptides,
and verifies that:

1. The fitted k_deg recovers the true value within ±10% on at least
   90% of test peptides.
2. The output ``spep`` column matches the coefficient-table sum.
3. The output DataFrame has the legacy columns plus the Phase F2 additions.

These are synthetic-noise-free tests; they verify the *fit pipeline math*,
not robustness to integration noise (that's ``bench_fit_recovery``'s job
in Phase F4). Bootstrap CI width is meaningless on noise-free data
(every bootstrap sample yields the same fit) — that's tested at the
``bench_fit_recovery`` level, not here.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from riana.algorithms.isotope_dist import (
    _get_final_env,
    _get_init_env,
    clear_envelope_cache,
    init_channel_masses,
)
from riana.algorithms.mass_calc import calculate_ion_mz
from riana.constants import PROTON_MASS
from riana.config import FitConfig
from riana.core.fitting import fit_run
from riana.core import models


_TEST_PEPTIDES = [
    ("VAPEPTIDEK", 2),
    ("LSHGLNVR", 2),
    ("AEFVEVTK", 2),
    ("STDNAATR", 2),
    ("HVELFK", 2),
]


def _coefficients_for_target_spep(peptides, spep_target: int) -> dict[str, float]:
    """Synthesize an AA-coefficient table so each test peptide's ``Σ
    coefficient[c]`` rounds to ``spep_target`` exactly.

    This is the production "Spep from sequence + table" path; the test
    just makes the table cooperate. For each unique AA across the test
    peptides, the coefficient is ``spep_target / mean_peptide_length`` so
    the sums end up near ``spep_target``.
    """
    lens = [len(seq) for seq, _ in peptides]
    avg_len = sum(lens) / len(lens)
    per_aa = spep_target / avg_len
    # Cover all 20 standard AAs so any peptide gets a defined coefficient.
    return {aa: per_aa for aa in "ACDEFGHIKLMNPQRSTVWY"}

# Synthetic experiment: 9 D2O proportions (mimicking the calibration series)
# at known pseudo-time under a simple exponential model with k_deg₀=0.5.
_K_DEG_TRUE = 0.5
_PROPORTIONS = np.array([0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875])
_TIMES = -np.log(1.0 - _PROPORTIONS) / _K_DEG_TRUE   # pseudo-times


def _make_synthetic_dfs(
    peptides: list[tuple[str, int]], spep_by_seq: dict[str, int],
) -> list[pd.DataFrame]:
    """Build one integrate-output DataFrame per timepoint with the known
    mixing-fraction envelopes for each peptide. Each peptide's Spep is
    looked up in ``spep_by_seq`` so the synthetic data matches whatever
    the production coefficient table will compute."""
    clear_envelope_cache()
    dfs = []
    for ti, prop in zip(_TIMES, _PROPORTIONS):
        rows = []
        for pep_idx, (seq, charge) in enumerate(peptides):
            pep_mass = calculate_ion_mz(seq)
            spep = spep_by_seq[seq]
            # The m0-m5 envelope the production fit consumes (matches the
            # `riana integrate --iso` default and core.fitting's required set).
            init = _get_init_env(seq, pep_mass, n=6)
            final = _get_final_env(seq, pep_mass, spep, ria_max=0.06, n=6)
            init_norm = init / init.sum()
            final_norm = final / final.sum()
            mix = (1.0 - prop) * init_norm + prop * final_norm
            # Scale to a synthetic intensity so the fit code path matches the
            # real "iso{N} = area" semantic (areas, not normalized fractions).
            scaled = mix * 1e6
            rows.append({
                "file_idx": 0,
                "scan": 1000 + pep_idx,
                "charge": charge,
                "concat": f"{seq}_{charge}",
                "sequence": seq,
                "sample": f"time{ti:.6f}",
                "percolator q-value": 1e-4,
                "protein id": f"sp|P0000{pep_idx}|TEST_HUMAN",
                "iso0": scaled[0], "iso1": scaled[1], "iso2": scaled[2],
                "iso3": scaled[3], "iso4": scaled[4], "iso5": scaled[5],
            })
        dfs.append(pd.DataFrame(rows))
    return dfs


def _spep_by_seq_from_coefficients(
    peptides, coefficients: dict[str, float]
) -> dict[str, int]:
    """Per-peptide integer Spep matching what fit_run will compute."""
    return {
        seq: max(1, int(round(sum(coefficients.get(c, 0.0) for c in seq))))
        for seq, _ in peptides
    }


def test_fit_run_recovers_k_deg_on_synthetic_data():
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    config = FitConfig(
        model="simple", label="hw", q_value=0.05, depth=3,
        ria_max=0.06,
    )

    result = fit_run(config, dfs, coeffs, n_boot=50, random_state=42)

    # Output schema: legacy + Phase F2 additions.
    for col in ("k_deg", "R_squared", "sd", "spep", "ci_lo", "ci_hi",
                "t", "fs", "protein id"):
        assert col in result.columns, f"missing {col}"

    assert len(result) == len(_TEST_PEPTIDES)

    # k_deg recovery — pure synthetic data should be tight.
    k_deg_arr = result["k_deg"].dropna().to_numpy()
    assert len(k_deg_arr) >= int(0.9 * len(_TEST_PEPTIDES)), (
        f"only {len(k_deg_arr)}/{len(_TEST_PEPTIDES)} peptides converged"
    )
    rel_err = (k_deg_arr - _K_DEG_TRUE) / _K_DEG_TRUE
    assert (np.abs(rel_err) < 0.10).mean() >= 0.9, (
        f"k_deg recovery too loose: rel_err = {rel_err}"
    )

    # Spep is deterministic from sequence × coefficient table.
    spep_arr = result["spep"].dropna().to_numpy()
    expected = [sum(coeffs.get(c, 0.0) for c in seq) for seq, _ in _TEST_PEPTIDES]
    np.testing.assert_allclose(np.sort(spep_arr), np.sort(expected), atol=1e-9)


def test_fit_run_fs_score_channels_runs_and_guards_missing_channels():
    """--fs limited-isotopomer scoring: fit_run accepts score_channels and still
    recovers k; and Guard 1 errors if the integrate output lacks a requested
    channel (user set integrate --iso too narrow)."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)

    cfg = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                    ria_max=0.06, score_channels=4)
    # (a) full iso0-5 integrate, score iso0-3 → runs and converges.
    result = fit_run(cfg, dfs, coeffs, n_boot=0, random_state=42)
    assert result["k_deg"].notna().sum() >= int(0.9 * len(_TEST_PEPTIDES))

    # (b) the relaxation: a MINIMAL integrate of exactly the scored channels
    # (iso0-3) is valid — you only need to extract what you score.
    minimal = [d.drop(columns=["iso4", "iso5"]) for d in dfs]
    res_min = fit_run(cfg, minimal, coeffs, n_boot=0, random_state=42)
    assert res_min["k_deg"].notna().sum() >= int(0.9 * len(_TEST_PEPTIDES))

    # (c) Guard 1: integrate too narrow (iso0-2) but --fs asks iso0-3 → clear error.
    narrow = [d.drop(columns=["iso3", "iso4", "iso5"]) for d in dfs]
    with pytest.raises(ValueError, match=r"iso3|--iso"):
        fit_run(cfg, narrow, coeffs, n_boot=0, random_state=42)


def test_fit_run_fs_auto_runs_and_guards_base_channels():
    """--fs auto: per-peptide init-width-keyed scoring runs and converges, and
    Guard 1 needs only the base (iso0-3) present (auto widens *within* whatever
    was captured)."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)

    cfg = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                    ria_max=0.06, fs_auto=True)
    result = fit_run(cfg, dfs, coeffs, n_boot=0, random_state=42)
    assert result["k_deg"].notna().sum() >= int(0.9 * len(_TEST_PEPTIDES))

    # A minimal iso0-3 integrate is valid under auto (only the base is required).
    minimal = [d.drop(columns=["iso4", "iso5"]) for d in dfs]
    res_min = fit_run(cfg, minimal, coeffs, n_boot=0, random_state=42)
    assert res_min["k_deg"].notna().sum() >= int(0.9 * len(_TEST_PEPTIDES))

    # Guard 1: integrate narrower than the base (iso0-2) → clear error.
    narrow = [d.drop(columns=["iso3", "iso4", "iso5"]) for d in dfs]
    with pytest.raises(ValueError, match=r"iso3|--iso"):
        fit_run(cfg, narrow, coeffs, n_boot=0, random_state=42)


def test_fit_config_fs_auto_and_score_channels_mutually_exclusive():
    """--fs auto and an explicit --fs channel count cannot both be set."""
    with pytest.raises(ValueError, match="mutually exclusive"):
        FitConfig(fs_auto=True, score_channels=4)


def test_fit_run_workers_deterministic():
    """The process-pool fit (``workers>1``) must match the serial path exactly.

    The per-peptide bootstrap is seeded from a content hash of the ``concat``
    (not the worker), so ``-W`` changes throughput, never the numbers. Noise is
    injected so residuals — hence the bootstrap CIs — are non-zero and the seed
    is actually observable.
    """
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    iso_cols = ["iso0", "iso1", "iso2", "iso3", "iso4", "iso5"]
    rng = np.random.default_rng(7)
    for df in dfs:
        noisy = df[iso_cols].to_numpy() * (1 + rng.normal(0, 0.02, df[iso_cols].shape))
        df[iso_cols] = noisy

    serial = fit_run(
        FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                  ria_max=0.06, workers=1),
        dfs, coeffs, n_boot=50, random_state=42,
    ).sort_index()
    pooled = fit_run(
        FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                  ria_max=0.06, workers=2),
        dfs, coeffs, n_boot=50, random_state=42,
    ).sort_index()

    assert list(serial.index) == list(pooled.index)
    for col in ("k_deg", "R_squared", "ci_lo", "ci_hi", "spep"):
        np.testing.assert_allclose(
            pd.to_numeric(serial[col], errors="coerce").to_numpy(),
            pd.to_numeric(pooled[col], errors="coerce").to_numpy(),
            rtol=1e-9, atol=1e-9, equal_nan=True,
            err_msg=f"{col} differs between workers=1 and workers=2",
        )


def test_fit_run_filters_by_depth():
    """Peptides with fewer than `depth` timepoints are filtered out."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES[:1], coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES[:1], spep_by_seq=spep_by_seq)
    # Only feed 2 timepoints — fewer than depth=3.
    short_dfs = dfs[:2]
    config = FitConfig(
        model="simple", label="hw", q_value=0.05, depth=3,
        ria_max=0.06,
    )
    with pytest.raises(ValueError, match="No peptides survive"):
        fit_run(config, short_dfs, coeffs, n_boot=10)


def test_depth_counts_distinct_timepoints_not_psm_rows():
    """Modern path: depth is distinct labeling timepoints, not raw PSM rows.

    A peptidoform with >= depth rows but < depth distinct timepoints (repeated in
    2D-LC fractions / technical replicates at one timepoint) must be filtered —
    counting rows would admit a curve a one-exponent fit can't identify.
    """
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES[:1], coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES[:1], spep_by_seq=spep_by_seq)
    for ti, df in zip(_TIMES, dfs):
        df["labeling_time"] = float(ti)
        df["biological_replicate"] = 1
    cfg = FitConfig(model="simple", label="hw", q_value=0.05, depth=3, ria_max=0.06)

    # 3 PSM rows but only 2 distinct timepoints (t0 repeated) -> filtered out.
    repeated = [dfs[0], dfs[1], dfs[0].copy()]
    with pytest.raises(ValueError, match="No peptides survive"):
        fit_run(cfg, repeated, coeffs, n_boot=10, time_column="labeling_time")

    # 3 distinct timepoints -> qualifies.
    out = fit_run(cfg, dfs[:3], coeffs, n_boot=10, time_column="labeling_time")
    assert len(out) == 1


def test_fit_run_rejects_unknown_model():
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES[:1], coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES[:1], spep_by_seq=spep_by_seq)
    bad_config = FitConfig(
        model="simple", label="hw", q_value=0.05, depth=3,
        ria_max=0.06,
    )
    object.__setattr__(bad_config, "model", "nonexistent")
    with pytest.raises(ValueError, match="unknown kinetic model"):
        fit_run(bad_config, dfs, coeffs, n_boot=10)


def test_fit_run_requires_canonical_isotopomers():
    """Integrate output missing m0-m5 (e.g. the legacy `--iso 0 6` pair) is
    rejected with a clear message rather than silently misaligning the envelope."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    # Drop the higher isotopomers so only m0-m3 remain (mimics too-narrow --iso).
    dfs = [df.drop(columns=["iso4", "iso5"]) for df in dfs]
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)
    with pytest.raises(ValueError, match="D2O fit needs isotopomers"):
        fit_run(config, dfs, coeffs, n_boot=10)
    # The guard names the active label (o18 vs D2O), not always "D2O".
    o18_cfg = FitConfig(model="simple", label="o18", q_value=0.05, depth=3,
                        ria_max=0.06)
    with pytest.raises(ValueError, match="o18 fit needs isotopomers"):
        fit_run(o18_cfg, dfs, coeffs, n_boot=10)


def test_fit_run_handles_bracketed_modification_strings():
    """Concat IDs from search engines often carry [mass] mod annotations
    (phospho, oxidation, etc.). fit_run must not crash on them — the
    bracketed-mass mod is included in pep_mass (via calculate_ion_mz),
    the envelope is computed on the stripped backbone for now (proper
    PTM forward-modelling is a post-M4 planning item)."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    # Inject a bracketed mod into one peptide's concat across all timepoints,
    # mimicking what a real Crux Percolator output would contain after a
    # variable-mod search.
    for df in dfs:
        m = df["concat"] == "VAPEPTIDEK_2"
        df.loc[m, "concat"] = "VAPEPTIDES[79.9663]K_2"
        df.loc[m, "sequence"] = "VAPEPTIDES[79.9663]K"

    config = FitConfig(
        model="simple", label="hw", q_value=0.05, depth=3,
        ria_max=0.06,
    )
    result = fit_run(config, dfs, coeffs, n_boot=10, random_state=42)

    # All 5 peptides should appear in the output (no crash; bracketed
    # peptide may or may not converge depending on synthetic-envelope
    # internal consistency, but should not raise).
    assert len(result) == len(_TEST_PEPTIDES)
    # Bracketed concat survives in the output index.
    assert any("[79.9663]" in c for c in result.index)


def test_fit_run_emits_fractions_long_with_prediction_intervals():
    """M5: fit_run attaches a tidy per-timepoint fraction-new table with PI bounds."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)
    result = fit_run(config, dfs, coeffs, n_boot=200, random_state=42)

    long = result.attrs["fractions_long"]
    assert list(long.columns) == [
        "concat", "protein id", "mod sites", "biological_replicate", "labeling_time",
        "fs", "fs_lower", "fs_upper", "fs_ds", "evidence", "metox",
    ]
    # Synthetic data has no MBR -> every point is a direct ID.
    assert (long["evidence"] == "q_value").all()
    # One row per (peptide, timepoint); on clean data all 5 peptides converge
    # over all 8 timepoints, and the long count matches the wide t list-cells.
    assert len(long) == sum(len(t) for t in result["t"])
    assert len(long) == 8 * len(_TEST_PEPTIDES)
    # The band is an interval (lower <= upper) and, on noise-free data, tight.
    assert (long["fs_lower"] <= long["fs_upper"] + 1e-9).all()
    width = (long["fs_upper"] - long["fs_lower"]).to_numpy()
    assert np.nanmedian(width) < 0.05
    # No SDRF identity on this path -> biorep defaults to 1.
    assert (long["biological_replicate"] == 1).all()
    # The wide frame also carries the PI list-cells for the GUI curve view.
    assert {"fs_lower", "fs_upper"} <= set(result.columns)


def _add_obs_mz(dfs, peptides, shift_per_channel):
    """Inject ``iso{k}_obs_mz`` columns at the init reference + a per-channel m/z
    shift, so the v1.1.0 item-1a Δmass path has known-answer inputs.

    ``shift_per_channel`` is a sequence of per-channel m/z offsets (length = number
    of iso channels) added to the init (θ=0) averaged-isotopolog reference m/z.
    The same shift is used at every timepoint — enough to check the arithmetic.
    """
    n = len(shift_per_channel)
    for df in dfs:
        cols = {f"iso{k}_obs_mz": [] for k in range(n)}
        for _, row in df.iterrows():
            seq = row["sequence"]
            z = int(row["charge"])
            pep_mass = calculate_ion_mz(seq)
            ref_neutral = init_channel_masses(seq, pep_mass, n)
            for k in range(n):
                ref_mz = (ref_neutral[k] + z * PROTON_MASS) / z
                cols[f"iso{k}_obs_mz"].append(ref_mz + shift_per_channel[k])
        for c, vals in cols.items():
            df[c] = vals
    return dfs


def test_init_channel_masses_reference():
    """The θ=0 reference: iso0 is the monoisotopic mass, channels step by ~1 Da."""
    clear_envelope_cache()
    seq, n = "VAPEPTIDEK", 6
    pep_mass = calculate_ion_mz(seq)
    masses = init_channel_masses(seq, pep_mass, n)
    assert len(masses) == n
    # iso0's averaged-isotopolog mass is dominated by the monoisotopic peak.
    assert masses[0] == pytest.approx(pep_mass, abs=2e-3)
    # Successive channels are ~one neutron apart (slightly > 1 from the averaging).
    steps = np.diff(masses)
    assert np.all(steps > 0.99) and np.all(steps < 1.02)


def test_dmass_dspacing_known_answer_and_global_drift_cancels():
    """Item 1a: dmass = obs − init_ref; dspacing is M0-internal so a uniform m/z
    drift cancels out of it while dmass records it."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    cfg = FitConfig(model="simple", label="hw", q_value=0.05, depth=3, ria_max=0.06)

    # (a) Uniform +5 mDa drift on every channel -> dmass ≈ 5 mDa everywhere;
    #     dspacing ≈ 0 (the M0-internal subtraction removes the global offset).
    drift = 5e-3  # m/z
    dfs = _add_obs_mz(
        _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq),
        _TEST_PEPTIDES, [drift] * 6)
    long = fit_run(cfg, dfs, coeffs, n_boot=0, random_state=42).attrs["fractions_long"]
    for k in range(6):
        assert f"dmass_iso{k}" in long.columns
        assert f"dspacing_iso{k}" in long.columns
    # charges in _TEST_PEPTIDES are all 2 -> the m/z drift is the same number.
    np.testing.assert_allclose(long["dmass_iso3"].to_numpy(), 5.0, atol=1e-6)
    np.testing.assert_allclose(long["dspacing_iso3"].to_numpy(), 0.0, atol=1e-6)
    assert (long["dspacing_iso0"] == 0.0).all()  # iso0 spacing is 0 by definition

    # (b) A per-channel ramp k*2 mDa -> dspacing_iso{k} ≈ k*2 mDa (relative to iso0).
    ramp = [k * 2e-3 for k in range(6)]
    dfs2 = _add_obs_mz(
        _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq),
        _TEST_PEPTIDES, ramp)
    long2 = fit_run(cfg, dfs2, coeffs, n_boot=0, random_state=42).attrs["fractions_long"]
    for k in range(6):
        np.testing.assert_allclose(long2[f"dspacing_iso{k}"].to_numpy(), k * 2.0,
                                   atol=1e-6)
        np.testing.assert_allclose(long2[f"dmass_iso{k}"].to_numpy(), k * 2.0,
                                   atol=1e-6)


def test_dmass_is_additive_no_effect_on_k_deg():
    """Regression: the Δmass path is purely additive — adding obs_mz columns leaves
    the FS / k_deg output identical (and omitting them omits the Δ columns)."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    cfg = FitConfig(model="simple", label="hw", q_value=0.05, depth=3, ria_max=0.06)

    base = fit_run(cfg, _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq),
                   coeffs, n_boot=0, random_state=42)
    with_mz = fit_run(
        cfg,
        _add_obs_mz(_make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq),
                    _TEST_PEPTIDES, [1e-3] * 6),
        coeffs, n_boot=0, random_state=42)

    # The FS solve is byte-identical; only floating-point summation order differs.
    np.testing.assert_allclose(base["k_deg"].to_numpy(), with_mz["k_deg"].to_numpy(),
                               rtol=1e-9, atol=1e-12, equal_nan=True)
    # No obs_mz -> no Δ columns; obs_mz present -> Δ columns appear.
    base_long = base.attrs["fractions_long"]
    assert not any(c.startswith(("dmass_", "dspacing_")) for c in base_long.columns)
    assert any(c.startswith("dmass_") for c in with_mz.attrs["fractions_long"].columns)


def _make_calibration_dfs(peptides, spep_by_seq):
    """Mixing-calibration series: ``sample`` encodes the heavy fraction f (0..1)
    directly (``time{f}``), so the legacy x-parse yields the proportion and a
    perfect integrator recovers FS == f — a 1:1 line."""
    clear_envelope_cache()
    dfs = []
    for f in _PROPORTIONS:
        rows = []
        for k, (seq, charge) in enumerate(peptides):
            pep_mass = calculate_ion_mz(seq)
            init = _get_init_env(seq, pep_mass, n=6)
            final = _get_final_env(seq, pep_mass, spep_by_seq[seq], ria_max=0.06, n=6)
            mix = (1 - f) * (init / init.sum()) + f * (final / final.sum())
            scaled = mix * 1e6
            rows.append({
                "file_idx": 0, "charge": charge, "concat": f"{seq}_{charge}",
                "sequence": seq, "sample": f"time{f:.6f}",
                "percolator q-value": 1e-4, "protein id": f"sp|P{k}|TEST",
                **{f"iso{j}": scaled[j] for j in range(6)},
            })
        dfs.append(pd.DataFrame(rows))
    return dfs


def test_theta_delta_s_known_answer_and_anchor():
    """_theta_delta_s: weighted-median of ΔSₓ/ΔSₓmax over iso1-3, t0-anchored."""
    from riana.core.fitting import _theta_delta_s
    dsmax = (0.0, 1.0, 2.0, 3.0)          # mDa per channel (iso0 unused)
    t = np.array([0.0, 0.5, 1.0])
    # Perfect signal: row k = f · dsmax[k]  →  θ = f after /ΔSₓmax.
    rows = [[0.0, 0.0, 0.0, 0.0],
            [0.0, 0.5, 1.0, 1.5],
            [0.0, 1.0, 2.0, 3.0]]
    np.testing.assert_allclose(_theta_delta_s(rows, t, dsmax), [0.0, 0.5, 1.0], atol=1e-9)
    # Anchoring removes a per-peptide constant offset on every channel (here +0.4):
    rows_off = [[0.0, 0.4, 0.4, 0.4],
                [0.0, 0.9, 1.4, 1.9],
                [0.0, 1.4, 2.4, 3.4]]
    np.testing.assert_allclose(_theta_delta_s(rows_off, t, dsmax), [0.0, 0.5, 1.0],
                               atol=1e-9)
    # A channel with negligible ΔSₓmax is dropped (needs ≥2 usable channels).
    assert np.isnan(_theta_delta_s([[0.0, 1.0, 0.0, 0.0]], np.array([1.0]),
                                   (0.0, 1.0, 0.0, 0.0)))


def test_delta_spacing_max_is_positive_monotone_for_d2o():
    """ΔSₓmax(k): iso0≡0, grows with channel (D heavier than ¹³C → positive)."""
    from riana.algorithms.isotope_dist import delta_spacing_max, clear_envelope_cache
    clear_envelope_cache()
    seq = "VAPEPTIDEK"
    dsmax = delta_spacing_max(seq, calculate_ion_mz(seq), spep=8, charge=2,
                              n=6, ria_max=0.06, label_int=1)
    assert dsmax[0] == 0.0
    assert dsmax[2] > dsmax[1] > 0          # positive, growing for D₂O
    assert all(np.isfinite(dsmax))


def test_fit_config_accepts_calibration_model():
    assert FitConfig(model="calibration").model == "calibration"
    with pytest.raises(ValueError, match="calibration"):
        FitConfig(model="nonsense")


def test_calibration_model_recovers_unit_slope():
    """--model calibration: a through-origin line of recovered FS vs known mixing
    proportion. On clean synthetic mixing data the slope ≈ 1 and R² ≈ 1."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_calibration_dfs(_TEST_PEPTIDES, spep_by_seq)
    cfg = FitConfig(model="calibration", label="hw", q_value=0.05, depth=3,
                    ria_max=0.06)
    out = fit_run(cfg, dfs, coeffs, n_boot=20, random_state=1)
    slopes = out["k_deg"].dropna().to_numpy()
    r2 = out["R_squared"].dropna().to_numpy()
    assert len(slopes) == len(_TEST_PEPTIDES)
    np.testing.assert_allclose(slopes, 1.0, atol=0.05)
    assert (r2 > 0.999).all()


def test_breakdown_columns_and_exclude_mbr():
    """n_points/n_mbr/n_metox/n_clean census + the --exclude-mbr filter."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    # Mark the first two timepoints' rows as MBR transfers (the rest direct IDs).
    for i, df in enumerate(dfs):
        df["evidence"] = "mbr" if i < 2 else "q_value"

    cfg = FitConfig(model="simple", label="hw", q_value=0.05, depth=3, ria_max=0.06)
    inc = fit_run(cfg, dfs, coeffs, n_boot=50, random_state=42)
    # 8 timepoints, 2 marked MBR -> every peptide curve is 2 MBR + 6 clean.
    assert (inc["n_points"] == 8).all()
    assert (inc["n_mbr"] == 2).all()
    assert (inc["n_metox"] == 0).all()
    assert (inc["n_clean"] == 6).all()
    assert inc.attrs["fractions_long"]["evidence"].eq("mbr").sum() == 2 * len(_TEST_PEPTIDES)

    cfg_excl = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                         ria_max=0.06, exclude_mbr=True)
    exc = fit_run(cfg_excl, dfs, coeffs, n_boot=50, random_state=42)
    assert (exc["n_points"] == 6).all()
    assert (exc["n_mbr"] == 0).all()
    assert (exc["n_clean"] == 6).all()
    assert exc.attrs["fractions_long"]["evidence"].eq("q_value").all()


def test_prediction_interval_widens_with_scatter():
    """The residual-bootstrap band tracks measurement scatter: injecting noise
    into the envelopes widens the per-timepoint prediction interval."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    clean = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)

    rng = np.random.default_rng(0)
    iso_cols = ["iso0", "iso1", "iso2", "iso3", "iso4", "iso5"]
    noisy = []
    for df in clean:
        d = df.copy()
        d[iso_cols] = d[iso_cols].to_numpy() * rng.lognormal(0.0, 0.05, d[iso_cols].shape)
        noisy.append(d)

    w_clean = (lambda L: np.nanmedian((L["fs_upper"] - L["fs_lower"]).to_numpy()))(
        fit_run(config, clean, coeffs, n_boot=200, random_state=7).attrs["fractions_long"]
    )
    w_noisy = (lambda L: np.nanmedian((L["fs_upper"] - L["fs_lower"]).to_numpy()))(
        fit_run(config, noisy, coeffs, n_boot=200, random_state=7).attrs["fractions_long"]
    )
    assert w_noisy > w_clean


def test_peptide_summary_drops_per_timepoint_list_cells():
    """The written peptides summary is scalar (no t/fs lists); the in-memory
    result keeps them for the GUI curve."""
    from riana.core.fitting import peptide_summary

    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES, coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES, spep_by_seq=spep_by_seq)
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)
    result = fit_run(config, dfs, coeffs, n_boot=10, random_state=42)

    summary = peptide_summary(result)
    assert not ({"t", "fs", "fs_lower", "fs_upper"} & set(summary.columns))
    assert {"k_deg", "R_squared", "spep", "protein id"} <= set(summary.columns)
    # The in-memory result is untouched (GUI curve reads t/fs from it).
    assert {"t", "fs"} <= set(result.columns)


def test_load_aa_coefficients_reads_csv(tmp_path):
    """load_aa_coefficients reads d2o_aa_coefficients CSVs."""
    csv = tmp_path / "coeffs.csv"
    csv.write_text("amino_acid,coefficient,boot_std\nA,0.5,0.05\nK,1.2,0.1\n")
    from riana.core.fitting import load_aa_coefficients
    out = load_aa_coefficients(csv)
    assert out == {"A": 0.5, "K": 1.2}


def test_fit_merges_met_ox_with_unoxidized_into_one_curve():
    """M7 tier 1b: a peptide seen unoxidized at some timepoints and Met-oxidized
    at others pools into ONE turnover curve (the chemical-mod-stripped fit key),
    so a pair that is each too shallow alone clears the depth gate together and
    recovers the planted k."""
    bare, charge, spep, k = "SAMMLPEPTIDEK", 2, 12, 0.5
    coeffs = _coefficients_for_target_spep([(bare, charge)], spep)
    clear_envelope_cache()
    mass_bare = calculate_ion_mz(bare)
    mass_ox = calculate_ion_mz("SAM[UNIMOD:35]MLPEPTIDEK")
    # 3 unoxidized timepoints + 3 oxidized = 6 pooled; each form alone is < depth 4.
    plan = [
        (0.5, "SAMMLPEPTIDEK_2", (), mass_bare),
        (1.0, "SAMMLPEPTIDEK_2", (), mass_bare),
        (1.5, "SAMMLPEPTIDEK_2", (), mass_bare),
        (2.0, "SAM[UNIMOD:35]MLPEPTIDEK_2", (35,), mass_ox),
        (3.0, "SAM[UNIMOD:35]MLPEPTIDEK_2", (35,), mass_ox),
        (4.0, "SAM[UNIMOD:35]MLPEPTIDEK_2", (35,), mass_ox),
    ]
    dfs = []
    for ti, concat, mods, pep_mass in plan:
        theta = 1.0 - np.exp(-k * ti)
        init = _get_init_env(bare, pep_mass, n=6, mods=mods)
        final = _get_final_env(bare, pep_mass, spep, ria_max=0.06, n=6, mods=mods)
        mix = (1 - theta) * (init / init.sum()) + theta * (final / final.sum())
        scaled = mix * 1e6
        dfs.append(pd.DataFrame([{
            "file_idx": 0, "scan": 1000, "charge": charge,
            "concat": concat, "sequence": bare, "sample": f"time{ti:.6f}",
            "percolator q-value": 1e-4, "protein id": "sp|P1|TEST",
            **{f"iso{n}": scaled[n] for n in range(6)},
        }]))
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=4,
                       ria_max=0.06)
    out = fit_run(config, dfs, coeffs, n_boot=50, random_state=42)
    # ONE merged result, keyed by the stripped fit key, pooling all 6 points.
    assert list(out.index) == ["SAMMLPEPTIDEK_2"]
    row = out.loc["SAMMLPEPTIDEK_2"]
    assert len(row["t"]) == 6
    assert row["k_deg"] == pytest.approx(k, abs=0.05)
