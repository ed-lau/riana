"""M7 Stage A2 — variable-mod (PTM) accounting through the forward model.

Covers the pieces that turn a variable-mod peptidoform into a correctly
integrated + fitted measurement, without needing a heavy mzML/parquet fixture:

- ``[UNIMOD:N]`` mass resolution in ``calculate_ion_mz`` (precursor m/z) and the
  ``unimod_mass`` / ``parse_unimod_ids`` helpers;
- the mzTab and DIA-NN ``_encode_peptidoform`` folders (encode starter-set mods,
  fold/strip the fixed CAM, drop unmodelable peptidoforms);
- ``solve_fs_d2o`` envelope threading (the mod reshapes the envelope; the
  bare-sequence path stays byte-identical).
"""

import numpy as np
import pytest

from riana.algorithms.isotope_dist import (
    _get_init_env,
    clear_envelope_cache,
    solve_fs_d2o,
)
from riana.algorithms.mass_calc import (
    calculate_ion_mz,
    parse_unimod_ids,
    unimod_mass,
)
from riana.io.diann import _encode_peptidoform as diann_encode
from riana.io.mztab import _encode_peptidoform as mztab_encode
from riana.io.mztab import _proteoform_sites


# --- mass resolution --------------------------------------------------------

def test_unimod_mass_matches_known_monoisotopic_shifts():
    # Phospho HPO3 and Acetyl C2H2O monoisotopic deltas.
    assert unimod_mass(21) == pytest.approx(79.96633, abs=1e-4)
    assert unimod_mass(1) == pytest.approx(42.01057, abs=1e-4)


def test_calculate_ion_mz_resolves_unimod_brackets():
    bare = calculate_ion_mz("SAMPLERK")
    assert calculate_ion_mz("SAMPLER[UNIMOD:21]K") - bare == pytest.approx(79.96633, abs=1e-4)
    assert calculate_ion_mz("[UNIMOD:1]SAMPLERK") - bare == pytest.approx(42.01057, abs=1e-4)


def test_calculate_ion_mz_still_handles_legacy_bracket_mass():
    # The demoted Percolator path encodes mods as raw [mass] floats.
    bare = calculate_ion_mz("SAMPLERK")
    assert calculate_ion_mz("SAM[15.99]PLERK") - bare == pytest.approx(15.99, abs=1e-6)
    assert calculate_ion_mz("SAM[0]PLERK") == pytest.approx(bare, abs=1e-9)


def test_calculate_ion_mz_rejects_non_numeric_non_unimod_bracket():
    with pytest.raises(ValueError):
        calculate_ion_mz("SAM[oops]PLERK")


def test_parse_unimod_ids_preserves_order_and_ignores_mass_brackets():
    assert parse_unimod_ids("[UNIMOD:1]SAMPLER[UNIMOD:21]K") == [1, 21]
    assert parse_unimod_ids("SAM[15.99]PLERK") == []
    assert parse_unimod_ids("SAMPLERK") == []


# --- mzTab encoder ----------------------------------------------------------

@pytest.mark.parametrize("modifications,expected", [
    (None, "SAMPLERK"),                              # unmodified
    ("null", "SAMPLERK"),                            # explicit null
    ("3-UNIMOD:4", "SAMPLERK"),                      # CAM only → folded, bare
    ("0-UNIMOD:1", "[UNIMOD:1]SAMPLERK"),            # N-term Acetyl
    ("0-UNIMOD:1,1-UNIMOD:4", "[UNIMOD:1]SAMPLERK"),  # Acetyl + CAM → encode Ac only
    ("3-UNIMOD:35", "SAM[UNIMOD:35]PLERK"),          # Met-Ox (tier 1b) → encoded
    ("0-UNIMOD:1,3-UNIMOD:35", "[UNIMOD:1]SAM[UNIMOD:35]PLERK"),  # Ac + Ox both kept
    ("2-UNIMOD:7", None),                            # deamidation → drop (unsupported)
])
def test_mztab_encode_peptidoform(modifications, expected):
    assert mztab_encode("SAMPLERK", modifications) == expected


def test_mztab_phospho_token_lands_after_the_modified_residue():
    # RVKSPEPVTSHPK, phospho at residue 4 (the S) → token immediately after it.
    assert mztab_encode("RVKSPEPVTSHPK", "4-UNIMOD:21") == "RVKS[UNIMOD:21]PEPVTSHPK"


# --- mzTab proteoform site (Stage B) ----------------------------------------

@pytest.mark.parametrize("sequence,modifications,start,expected", [
    ("RVKSPEPVTSHPK", "4-UNIMOD:21", 34473, "pS34476"),   # TITIN S34476
    ("IGHHSTSDDSSAYR", "5-UNIMOD:21", 330, "pS334"),
    ("GEIEHHCSGLHR", "8-UNIMOD:21,7-UNIMOD:4", 50, "pS57"),  # CAM ignored
    ("SPSPK", "1-UNIMOD:21,3-UNIMOD:21", 200, "pS200_pS202"),  # two sites
    ("SAMPLERK", "0-UNIMOD:1,3-UNIMOD:4", 100, ""),       # N-term Ac + CAM → bare
    ("PEPTIDEK", "null", 10, ""),                          # unmodified → bare
    ("PEPTIDEK", "2-UNIMOD:21", None, ""),                 # no start → bare (no fabrication)
])
def test_mztab_proteoform_sites(sequence, modifications, start, expected):
    assert _proteoform_sites(sequence, modifications, start) == expected


# --- chemical-mod fit key (tier 1b) -----------------------------------------

@pytest.mark.parametrize("concat,expected", [
    ("PEPM[UNIMOD:35]TIDEK_2", "PEPMTIDEK_2"),                  # Met-Ox stripped
    ("PEPS[UNIMOD:21]M[UNIMOD:35]TIDEK_2", "PEPS[UNIMOD:21]MTIDEK_2"),  # phospho kept
    ("PEPMTIDEK_2", "PEPMTIDEK_2"),                            # no chemical mod → unchanged
    ("[UNIMOD:1]PEPM[UNIMOD:35]TIDEK_3", "[UNIMOD:1]PEPMTIDEK_3"),  # N-term Ac kept, Ox stripped
])
def test_fit_key_strips_only_chemical_mods(concat, expected):
    from riana.core.fitting import _fit_key
    assert _fit_key(concat) == expected


# --- DIA-NN encoder ---------------------------------------------------------

@pytest.mark.parametrize("modified_sequence,expected", [
    ("AAAFEQLQK", "AAAFEQLQK"),                              # unmodified
    ("AAADEWTTC(UniMod:4)TPPSGLQGK", "AAADEWTTCTPPSGLQGK"),  # CAM stripped (folded)
    ("(UniMod:1)AACDEFK", "[UNIMOD:1]AACDEFK"),              # N-term Acetyl
    ("AAS(UniMod:21)PEPK", "AAS[UNIMOD:21]PEPK"),            # phospho in place
    ("AAC(UniMod:4)S(UniMod:21)PEPK", "AACS[UNIMOD:21]PEPK"),  # CAM stripped, phospho kept
    ("AAM(UniMod:35)PEPK", "AAM[UNIMOD:35]PEPK"),           # Met-Ox (tier 1b) → encoded
    ("AAN(UniMod:7)PEPK", None),                            # deamidation → drop (unsupported)
])
def test_diann_encode_peptidoform(modified_sequence, expected):
    assert diann_encode(modified_sequence) == expected


# --- envelope threading -----------------------------------------------------

def test_mods_reshape_the_envelope_but_bare_path_is_unchanged():
    clear_envelope_cache()
    seq, pep_mass = "SAMPLERSTK", calculate_ion_mz("SAMPLERSTK")
    env_bare = np.asarray(_get_init_env(seq, pep_mass, n=6, mods=()))
    env_phos = np.asarray(_get_init_env(seq, pep_mass, n=6, mods=(21,)))
    # Phospho's 3 oxygens reshape the M+1/M+2 channels.
    assert float(np.max(np.abs(env_phos - env_bare))) > 1e-3
    # The empty-mods call must equal the legacy no-mods default exactly.
    env_default = np.asarray(_get_init_env(seq, pep_mass, n=6))
    assert np.array_equal(env_bare, env_default)


def test_solve_fs_d2o_accepts_mods_and_recovers_endpoints():
    # On a noise-free fully-labeled-vs-initial mixture for a phospho peptidoform,
    # the FS solve should recover ~0 and ~1 when fed the matching envelopes.
    clear_envelope_cache()
    seq, mods = "SAMPLERSTK", (21,)
    pep_mass = calculate_ion_mz("SAMPLER[UNIMOD:21]STK")
    init = np.asarray(_get_init_env(seq, pep_mass, n=4, mods=mods))
    final = np.asarray(_get_final_env_via_solver(seq, pep_mass, mods))
    fs0 = solve_fs_d2o(seq, pep_mass, init, spep=10, n_iso=4, mods=mods)
    fs1 = solve_fs_d2o(seq, pep_mass, final, spep=10, n_iso=4, mods=mods)
    assert fs0 == pytest.approx(0.0, abs=0.05)
    assert fs1 == pytest.approx(1.0, abs=0.05)


def _get_final_env_via_solver(seq, pep_mass, mods):
    from riana.algorithms.isotope_dist import _get_final_env
    return _get_final_env(seq, pep_mass, spep=10, ria_max=0.06, n=4, mods=mods)
