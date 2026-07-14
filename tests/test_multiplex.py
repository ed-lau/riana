"""Sample-axis multiplexing (dimethyl duplex) — registry + channel→sample intake.

Phase 1-2 of the 1.2.0 multiplexing build: the ``riana.multiplex`` label registry,
the heavy-dimethyl (UNIMOD:330) / DIMETHYL4 (199) pinned-isotope constants, and the
SDRF/mzTab channel intake that keeps light and heavy as distinct samples fit at their
own precursor enrichment. See the plan + ``reports/2026-07-06_dimethyl_duplex_spillover.md``.
"""
from __future__ import annotations

import textwrap

import pytest

from riana import multiplex
from riana.algorithms.isotope_dist import adaptive_channel_masses, get_envelope, get_peptide_distribution
from riana.algorithms.mass_calc import calculate_ion_mz, parse_unimod_ids, unimod_mass
from riana.core.pipeline import plan_integration
from riana.config import IntegrationConfig
from riana.exceptions import DataError
from riana.io.sdrf import read_sdrf
from riana.io.mztab import read_mztab
from riana.utils import strip_concat


# --- the label registry ------------------------------------------------------

def test_channel_of_routes_by_dimethyl_mod():
    assert multiplex.channel_of([36]) == ("dimethyl", "light")
    assert multiplex.channel_of([36, 36]) == ("dimethyl", "light")   # N-term + K, one channel
    assert multiplex.channel_of([199]) == ("dimethyl", "medium")
    assert multiplex.channel_of([330]) == ("dimethyl", "heavy")
    # non-multiplex mods (CAM/Met-Ox/phospho) are ignored; a bare peptide → None
    assert multiplex.channel_of([4, 35, 21]) is None
    assert multiplex.channel_of([]) is None
    # a channel mod alongside a biological mod still routes by the channel
    assert multiplex.channel_of([330, 21]) == ("dimethyl", "heavy")


def test_channel_of_rejects_conflicting_channels():
    with pytest.raises(ValueError, match="conflicting"):
        multiplex.channel_of([36, 330])  # a PSM cannot be both light and heavy


def test_cv_term_detection_and_mapping():
    assert multiplex.is_multiplex_cv("NT=DIMETHYL0;AC=PRIDE:0000848")
    assert multiplex.is_multiplex_cv("NT=DIMETHYL8;AC=PRIDE:0000852")
    assert not multiplex.is_multiplex_cv("label free sample")
    assert not multiplex.is_multiplex_cv("NT=TMT126;AC=MS:...")   # isobaric ≠ multiplex
    assert not multiplex.is_multiplex_cv("")
    assert multiplex.label_channel_for_cv("NT=DIMETHYL0;AC=x") == ("dimethyl", "light")
    assert multiplex.label_channel_for_cv("NT=DIMETHYL8;AC=x") == ("dimethyl", "heavy")
    assert multiplex.label_channel_for_cv("label free sample") is None


def test_silac_registered_as_geometry_only():
    silac = next(l for l in multiplex.REGISTRY if l.name == "silac")
    assert silac.intake_wired is False          # no CV/mod intake yet
    assert silac.labels_n_term is False and silac.site_residues == ("K", "R")
    heavy = next(c for c in silac.channels if c.name == "heavy")
    # heterogeneous per-residue shift (the point of registering it): K ≠ R
    assert heavy.shift_per_site["K"] == pytest.approx(6.0201, abs=1e-3)
    assert heavy.shift_per_site["R"] == pytest.approx(10.0083, abs=1e-3)


# --- pinned-isotope constants (mass + envelope lockstep) ---------------------

@pytest.mark.parametrize("uid, expect", [(36, 28.0313), (199, 32.0564), (330, 36.0757)])
def test_dimethyl_unimod_masses(uid, expect):
    assert unimod_mass(uid) == pytest.approx(expect, abs=1e-3)


def test_heavy_minus_light_is_plus_eight_per_site():
    assert unimod_mass(330) - unimod_mass(36) == pytest.approx(8.0444, abs=1e-3)
    assert unimod_mass(199) - unimod_mass(36) == pytest.approx(4.0251, abs=1e-3)


@pytest.mark.parametrize("concat", ["[UNIMOD:36]WCALSHLER",
                                    "[UNIMOD:199]WCALSHLER",
                                    "[UNIMOD:330]WCALSHLER"])
def test_envelope_anchor_matches_precursor_mass(concat):
    """The pinned-isotope pseudo-elements must keep the precursor mass
    (``calculate_ion_mz``) and the IsoSpec envelope m0 in lockstep — else the
    adaptive-N_ISO guard rejects the peptidoform and the fit falls back."""
    mods = parse_unimod_ids(concat)
    pep_mass = calculate_ion_mz(concat, ion="M", charge=0)
    bare = strip_concat(concat)
    dist = get_peptide_distribution(bare, label="D2O", mods=mods)
    assert min(dist.masses) == pytest.approx(pep_mass, abs=1e-3)
    # the envelope is populated (not collapsed to a wrong bin) and the guard passes
    assert sum(get_envelope(dist, pep_mass, n=6)) > 0.99
    assert len(adaptive_channel_masses(bare, pep_mass, mods=tuple(mods))) >= 3


# --- SDRF channel→sample intake ----------------------------------------------

_DM_SDRF = textwrap.dedent("""\
    source name\tcharacteristics[biological replicate]\tcharacteristics[precursor enrichment]\tcharacteristics[labeling time]\tcomment[label]\tcomment[data file]\tcomment[proteomics data acquisition method]\tfactor value[condition]
    ctrl_r1\t1\t0.04614\t8 day\tNT=DIMETHYL0;AC=PRIDE:0000848\tfileA.mzML\tDDA\tcontrol
    rapa_r1\t1\t0.05611\t8 day\tNT=DIMETHYL8;AC=PRIDE:0000852\tfileA.mzML\tDDA\trapamycin
    ctrl_r2\t2\t0.04614\t8 day\tNT=DIMETHYL0;AC=PRIDE:0000848\tfileB.mzML\tDDA\tcontrol
    rapa_r2\t2\t0.05611\t8 day\tNT=DIMETHYL8;AC=PRIDE:0000852\tfileB.mzML\tDDA\trapamycin
    """)


def _write_dm_sdrf(tmp_path):
    p = tmp_path / "dm.sdrf.tsv"
    p.write_text(_DM_SDRF)
    return read_sdrf(p)


def test_multiplex_sdrf_keeps_channels_as_distinct_samples(tmp_path):
    """A dimethyl SDRF (2 rows/file) is NOT collapsed (unlike TMT): every channel
    stays a distinct RunIdentity with its own condition + precursor enrichment."""
    t = _write_dm_sdrf(tmp_path)
    assert t.is_multiplexed
    assert len(t.runs) == 4  # 2 files × 2 channels, kept
    cmap = t.multiplex_channel_map
    light = cmap[("fileA", ("dimethyl", "light"))]
    heavy = cmap[("fileA", ("dimethyl", "heavy"))]
    assert light.condition == "control" and light.precursor_enrichment == 0.04614
    assert heavy.condition == "rapamycin" and heavy.precursor_enrichment == 0.05611
    assert light.sample == "ctrl_r1" and heavy.sample == "rapa_r1"


def test_non_multiplex_sdrf_has_no_channel_map(tmp_path):
    src = _DM_SDRF.replace("NT=DIMETHYL0;AC=PRIDE:0000848", "label free sample")
    src = src.replace("NT=DIMETHYL8;AC=PRIDE:0000852", "label free sample")
    # now fileA/fileB each appear twice with a non-channel label → duplicate guard fires
    p = tmp_path / "lf.sdrf.tsv"
    p.write_text(src)
    with pytest.raises(DataError, match="unique per run"):
        read_sdrf(p)


# --- mzTab channel routing + per-channel planning ----------------------------

_DM_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t1.0.0
    MTD\tms_run[1]-location\tfile://fileA.mzML
    MTD\tms_run[2]-location\tfile://fileB.mzML

    PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score[1]\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tspectra_ref\tpre\tpost\tstart\tend\topt_global_Posterior_Error_Probability_score\topt_global_q-value\topt_global_cv_MS:1002217_decoy_peptide\topt_global_cv_MS:1000889_peptidoform_sequence
    PSM\tSAMPLER\t0\tsp|P1|X\t1\tdb\tnull\tse\t0.001\t0-UNIMOD:36\t500.0\t2\t400.0\t400.0\tms_run[1]:controllerType=0 controllerNumber=1 scan=10\tK\tR\t1\t7\t0.01\t0.001\t0\t.(Dimethyl)SAMPLER
    PSM\tSAMPLER\t1\tsp|P1|X\t1\tdb\tnull\tse\t0.001\t0-UNIMOD:330\t500.0\t2\t404.0\t404.0\tms_run[1]:controllerType=0 controllerNumber=1 scan=11\tK\tR\t1\t7\t0.01\t0.001\t0\t.(Dimethyl:2H(6)13C(2))SAMPLER
    PSM\tSAMPLER\t2\tsp|P1|X\t1\tdb\tnull\tse\t0.001\t0-UNIMOD:36\t500.0\t2\t400.0\t400.0\tms_run[2]:controllerType=0 controllerNumber=1 scan=12\tK\tR\t1\t7\t0.01\t0.001\t0\t.(Dimethyl)SAMPLER
    PSM\tSAMPLER\t3\tsp|P1|X\t1\tdb\tnull\tse\t0.001\t0-UNIMOD:330\t500.0\t2\t404.0\t404.0\tms_run[2]:controllerType=0 controllerNumber=1 scan=13\tK\tR\t1\t7\t0.01\t0.001\t0\t.(Dimethyl:2H(6)13C(2))SAMPLER
    """)


def test_mztab_routes_each_psm_to_its_channel_identity(tmp_path):
    sdrf = _write_dm_sdrf(tmp_path)
    mzt = tmp_path / "dm.mzTab"
    mzt.write_text(_DM_MZTAB)
    psms, fidx = read_mztab(mzt, channel_map=sdrf.multiplex_channel_map)
    assert len(psms) == 4
    by_cond = {}
    for p in psms:
        by_cond.setdefault(p.identity.condition, []).append(p)
    assert set(by_cond) == {"control", "rapamycin"}
    # light [UNIMOD:36] → control@0.04614; heavy [UNIMOD:330] → rapamycin@0.05611
    for p in by_cond["control"]:
        assert 36 in parse_unimod_ids(p.sequence) and p.identity.precursor_enrichment == 0.04614
    for p in by_cond["rapamycin"]:
        assert 330 in parse_unimod_ids(p.sequence) and p.identity.precursor_enrichment == 0.05611


def test_plan_integration_emits_one_task_per_channel(tmp_path):
    sdrf = _write_dm_sdrf(tmp_path)
    mzt = tmp_path / "dm.mzTab"
    mzt.write_text(_DM_MZTAB)
    mzml_dir = tmp_path / "mzml"
    mzml_dir.mkdir()
    for f in ("fileA.mzML", "fileB.mzML"):
        (mzml_dir / f).write_text("")  # plan only checks existence
    tasks = plan_integration(IntegrationConfig(), sdrf, mzml_dir, mzt)
    assert len(tasks) == 4  # 2 files × 2 channels
    # two channels share a physical mzML but have distinct output stems + identities
    by_stem = {t.stem: t for t in tasks}
    assert set(by_stem) == {"ctrl_r1", "rapa_r1", "ctrl_r2", "rapa_r2"}
    assert by_stem["ctrl_r1"].mzml_path == by_stem["rapa_r1"].mzml_path  # same file
    assert by_stem["ctrl_r1"].identity.condition == "control"
    assert by_stem["rapa_r1"].identity.condition == "rapamycin"
    assert len(by_stem["ctrl_r1"].psms) == 1 and len(by_stem["rapa_r1"].psms) == 1
