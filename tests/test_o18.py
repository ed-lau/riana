"""¹⁸O (label=o18) fit-side regression guards.

Covers the production o18 path added in v1.1.0: the 3-isotope enriched-O
envelope (``get_peptide_distribution`` ``label="O18"``), ``solve_fs_o18``, the
length-model Spep (``spep_from_length_coefficients``), and the bundled
``juber_2026_o18_ac16`` coefficient preset. The science is validated against the NB90c
reverse model (``tests/benchmark/_helpers/o18_forward_model.py``); these tests
just pin the production wiring so it cannot silently regress.
"""
from __future__ import annotations

import numpy as np
import pytest

from riana.algorithms.isotope_dist import (
    O18_LENGTH_FEATURES,
    get_envelope,
    get_peptide_distribution,
    solve_fs_o18,
    spep_from_length_coefficients,
)
from riana.algorithms.mass_calc import calculate_ion_mz
from riana.core.fitting import load_o18_coefficients

_SEQ = "TDLEKDIISDTSGDFR"
_RIA_O18 = 0.058315  # AC16 6% v/v H2(18)O
_SPEP = 7


def _envelope(label_kwargs: dict) -> np.ndarray:
    mass = calculate_ion_mz(_SEQ)
    dist = get_peptide_distribution(_SEQ, **label_kwargs)
    env = np.array(get_envelope(dist, mass, n=8))[:6]
    return env / env.sum()


def test_o18_envelope_is_three_isotope_and_shifts_mass():
    """The labeled o18 envelope (label="O18") moves probability off iso0 (+2 Da/site)."""
    init = _envelope(dict(label="D2O"))  # natural abundance (label-independent)
    final = _envelope(dict(label="O18", deuterium_enrichment_level=_RIA_O18,
                           num_labeling_sites=_SPEP))
    assert np.isclose(init.sum(), 1.0) and np.isclose(final.sum(), 1.0)
    # ¹⁸O incorporation depletes the monoisotopic peak and fills higher channels.
    assert final[0] < init[0]
    assert final[2] > init[2]


@pytest.mark.parametrize("fs_true", [0.0, 0.25, 0.5, 0.75, 1.0])
def test_o18_solve_fs_recovers_known_mix(fs_true):
    """solve_fs_o18 recovers fs from an init/final mixture at known fs."""
    init = _envelope(dict(label="D2O"))
    final = _envelope(dict(label="O18", deuterium_enrichment_level=_RIA_O18,
                           num_labeling_sites=_SPEP))
    obs = (1 - fs_true) * init + fs_true * final
    fs = solve_fs_o18(_SEQ, calculate_ion_mz(_SEQ), obs, _SPEP,
                      ria_max=_RIA_O18, n_iso=6)
    assert abs(fs - fs_true) < 0.05


def test_o18_spep_length_model():
    """Spep = b·(L-1) + c_D·D + c_E·E + c_N·N + c_Q·Q + c_S·S, intercept 0."""
    coef = {"length_minus1": 0.16, "D": 1.5, "E": 1.2, "N": 0.38, "Q": 0.03, "S": 0.52}
    seq = "DDES"  # L=4 → 0.16·3 + 1.5·2 + 1.2·1 + 0.52·1
    expected = 0.16 * 3 + 1.5 * 2 + 1.2 * 1 + 0.52 * 1
    assert spep_from_length_coefficients(seq, coef) == pytest.approx(expected)
    # an unknown feature key contributes nothing
    assert spep_from_length_coefficients("AAAA", {"D": 9.9}) == 0.0


def test_o18_preset_loads_all_features():
    """The bundled juber_2026_o18_ac16 preset carries every length-model feature."""
    coef = load_o18_coefficients("juber_2026_o18_ac16")
    assert set(coef) == set(O18_LENGTH_FEATURES)
    # carboxyl (D/E) dominate the in-vitro signal; backbone is small; Q ~ 0.
    assert coef["D"] > coef["S"] > coef["length_minus1"]
    assert coef["Q"] < 0.1


def test_o18_previs_preset_matches_paper_worked_example():
    """The in-vivo rachdaoui_2009_o18 (mouse) preset reproduces Previs' own worked
    example: LGEYGFQNAILVR → 12 peptide bonds + N + Q + 2·E = 16
    (Rachdaoui/Previs, Mol Cell Proteomics 2009; 8(12):2653). Backbone is 1 ¹⁸O/peptide
    bond in vivo, vs ~0.16 in the in-vitro AC16 table — the regime difference."""
    prev = load_o18_coefficients("rachdaoui_2009_o18")
    assert prev == {"length_minus1": 1.0, "D": 0.0, "E": 2.0,
                    "N": 1.0, "Q": 1.0, "S": 0.0}
    assert spep_from_length_coefficients("LGEYGFQNAILVR", prev) == 16.0


def test_o18_dspacing_iso1_flat_iso2_carries_signal():
    """¹⁸O (+2 Da) leaves iso1 flat — no ¹⁸O-bearing isotopolog reaches +1 — while
    iso2 carries the mass-defect shift (one ¹⁸O is ~2.5 mDa lighter than 2×¹³C)."""
    from riana.algorithms.isotope_dist import _spacing_components, _mixture_dspacing
    pm = calculate_ion_mz(_SEQ)
    im, ip, fm, fp = _spacing_components(_SEQ, pm, _SPEP, _RIA_O18, 6, (), "O18")
    s_iso1 = _mixture_dspacing(0.8, im, ip, fm, fp, 2, 1)
    s_iso2 = _mixture_dspacing(0.8, im, ip, fm, fp, 2, 2)
    assert abs(s_iso1) < abs(s_iso2)


def test_o18_solve_fs_ds_recovers_known_mix():
    """solve_fs_o18_ds inverts the ¹⁸O Δspacing (the label="O18" envelope) back to the
    known fraction-new on the iso2-4 channels — the spacing analog of solve_fs_o18.
    (Self-consistency against its own forward model: signal is small but nonzero, so
    the noiseless minimum is unique.)"""
    from riana.algorithms.isotope_dist import (
        _spacing_components, _mixture_dspacing, solve_fs_o18_ds)
    pm = calculate_ion_mz(_SEQ)
    z = 2
    im, ip, fm, fp = _spacing_components(_SEQ, pm, _SPEP, _RIA_O18, 6, (), "O18")
    for f_true in (0.25, 0.5, 0.75):
        obs = {k: _mixture_dspacing(f_true, im, ip, fm, fp, z, k) for k in (2, 3, 4)}
        f = solve_fs_o18_ds(_SEQ, pm, obs, _SPEP, z, ria_max=_RIA_O18, n_iso=6)
        assert abs(f - f_true) < 0.05
