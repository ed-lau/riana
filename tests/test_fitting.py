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
)
from riana.algorithms.mass_calc import calculate_ion_mz
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
            init = _get_init_env(seq, pep_mass, n=4)
            final = _get_final_env(seq, pep_mass, spep, ria_max=0.06, n=4)
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
                "iso0": scaled[0], "iso1": scaled[1],
                "iso2": scaled[2], "iso3": scaled[3],
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
        model="simple", label=1, q_value=0.05, depth=3,
        ria_max=0.06, threads=1,
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


def test_fit_run_filters_by_depth():
    """Peptides with fewer than `depth` timepoints are filtered out."""
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES[:1], coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES[:1], spep_by_seq=spep_by_seq)
    # Only feed 2 timepoints — fewer than depth=3.
    short_dfs = dfs[:2]
    config = FitConfig(
        model="simple", label=1, q_value=0.05, depth=3,
        ria_max=0.06, threads=1,
    )
    with pytest.raises(ValueError, match="No peptides survive"):
        fit_run(config, short_dfs, coeffs, n_boot=10)


def test_fit_run_rejects_unknown_model():
    coeffs = _coefficients_for_target_spep(_TEST_PEPTIDES, 8)
    spep_by_seq = _spep_by_seq_from_coefficients(_TEST_PEPTIDES[:1], coeffs)
    dfs = _make_synthetic_dfs(_TEST_PEPTIDES[:1], spep_by_seq=spep_by_seq)
    bad_config = FitConfig(
        model="simple", label=1, q_value=0.05, depth=3,
        ria_max=0.06, threads=1,
    )
    object.__setattr__(bad_config, "model", "nonexistent")
    with pytest.raises(ValueError, match="unknown kinetic model"):
        fit_run(bad_config, dfs, coeffs, n_boot=10)


def test_load_aa_coefficients_reads_csv(tmp_path):
    """load_aa_coefficients reads d2o_aa_coefficients CSVs."""
    csv = tmp_path / "coeffs.csv"
    csv.write_text("amino_acid,coefficient,boot_std\nA,0.5,0.05\nK,1.2,0.1\n")
    from riana.core.fitting import load_aa_coefficients
    out = load_aa_coefficients(csv)
    assert out == {"A": 0.5, "K": 1.2}
