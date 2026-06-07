"""Smoke tests for ``riana.io.*`` (Week 2 I/O layer).

Sample1 BSA fixture is the CI-cheap regression set; the heavyweight parity
checks against the calibration dataset run only in benchmarks. The mzTab
test uses a synthetic 2-PSM file written into a tempdir so CI does not
depend on the gigabyte-scale quantms output that lives outside the repo.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import numpy as np
import pytest

from riana.io import mzml as iomzml
from riana.io import mztab as iomztab
from riana.io import percolator as ioperc
from riana.io import writers as iowriters

FIXTURE_DIR = Path("tests/data/sample1")
PSMS = FIXTURE_DIR / "percolator.target.psms.txt"
MZML_GZ = FIXTURE_DIR / "20180216_BSA.mzML.gz"


# --- percolator --------------------------------------------------------------


def test_percolator_reads_crux_format():
    records = ioperc.read_percolator(PSMS, sample="sample1")
    assert len(records) > 0
    first = records[0]
    # Fields the downstream pipeline reads on every PSM.
    assert first.scan > 0
    assert first.charge > 0
    assert first.sequence
    assert first.peptide_mass > 0
    assert first.sample == "sample1"
    assert first.concat == f"{first.sequence}_{first.charge}"


def test_percolator_matches_readpercolator_parity():
    """The typed parser must be a drop-in for the 0.9.0 ReadPercolator.

    Compared against ``percolator_parity_v0_9_0.csv`` — the legacy
    ``ReadPercolator.id_df`` (peptide mass / scan / sequence) frozen before the
    legacy module was removed in M4. Anchored to the same parity we verified on
    the ac16 calibration data (peptide_mass diff 0.0); catches field-mapping
    regressions in the new parser.
    """
    import pandas as pd

    old = pd.read_csv(FIXTURE_DIR / "percolator_parity_v0_9_0.csv")
    new = ioperc.read_percolator(PSMS, sample="sample1")
    assert len(new) == len(old)
    # peptide mass is a recomputed float; the golden round-trips through CSV, so
    # compare at full float precision rather than bit-exact ==.
    np.testing.assert_allclose(
        np.array([r.peptide_mass for r in new]),
        old["peptide mass"].to_numpy(),
        rtol=1e-9, atol=1e-6,
    )
    assert (old["scan"].to_numpy() == np.array([r.scan for r in new])).all()
    assert (
        old["sequence"].to_numpy() == np.array([r.sequence for r in new])
    ).all()


def test_percolator_fraction_psms_assigns_pep_ids():
    records = ioperc.read_percolator(PSMS, sample="sample1")
    idx = ioperc.file_indices(records)[0]
    rows = ioperc.fraction_psms(records, idx)
    assert [r.pep_id for r in rows] == list(range(len(rows)))
    # Sorted by scan.
    assert [r.scan for r in rows] == sorted(r.scan for r in rows)


# --- mzml --------------------------------------------------------------------


def test_indexed_mzml_matches_pymzml_ms1_count():
    """pyteomics-backed reader sees the same MS1 scans as the 0.9.0 pymzml one.

    Compared against ``mzml_ms1_v0_9_0.npz`` — the legacy ``riana.spectra.Mzml``
    ``scan_idx`` / ``rt_idx`` and a mid-scan peak count, frozen before the
    legacy reader was removed in M4.
    """
    golden = np.load(FIXTURE_DIR / "mzml_ms1_v0_9_0.npz")
    old_scan_idx = golden["scan_idx"]
    old_rt_idx = golden["rt_idx"]
    mid_scan = int(golden["mid_scan"])
    mid_npeaks = int(golden["mid_npeaks"])
    with iomzml.IndexedMzML(MZML_GZ) as new:
        assert (new.scan_idx == old_scan_idx).all()
        # rt_idx is in minutes for both readers; tolerate 1e-9 for float quirks.
        assert np.max(np.abs(new.rt_idx - old_rt_idx)) < 1e-9
        # Random-access peaks fetch round-trips against the eager pymzml load.
        new_mz, _ = new.peaks(mid_scan)
        assert len(new_mz) == mid_npeaks


def test_indexed_mzml_ms1_iter_yields_full_set():
    with iomzml.IndexedMzML(MZML_GZ) as r:
        scans = [scan for scan, _, _, _ in r.ms1_iter()]
    assert scans
    assert sorted(scans) == scans  # MS1 emitted in file order


# --- mztab -------------------------------------------------------------------


_MINIMAL_MZTAB = textwrap.dedent(
    """\
    MTD\tmzTab-version\t1.0.0
    MTD\tmzTab-mode\tSummary
    MTD\tmzTab-type\tIdentification
    MTD\tdescription\tsynthetic test fixture
    MTD\tpsm_search_engine_score[1]\t[MS, MS:1003115, OpenMS:Target-decoy PSM q-value, ]
    MTD\tms_run[1]-format\t[MS, MS:1000584, mzML file, ]
    MTD\tms_run[1]-location\tfile://run_A.mzML
    MTD\tms_run[1]-id_format\t[MS, MS:1000768, Thermo nativeID format, ]
    MTD\tms_run[2]-format\t[MS, MS:1000584, mzML file, ]
    MTD\tms_run[2]-location\tfile://run_B.mzML
    MTD\tms_run[2]-id_format\t[MS, MS:1000768, Thermo nativeID format, ]

    PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score[1]\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tspectra_ref\tpre\tpost\tstart\tend\topt_global_Posterior_Error_Probability_score\topt_global_q-value\topt_global_cv_MS:1002217_decoy_peptide\topt_global_cv_MS:1000889_peptidoform_sequence
    PSM\tPEPTIDEK\t0\tsp|P00001|TEST_HUMAN\t1\tdb\tnull\t[, , dummy, 1]\t0.001\tnull\t10.5\t2\t472.74\t472.73\tms_run[1]:controllerType=0 controllerNumber=1 scan=101\tK\tR\t10\t17\t0.01\t0.001\t0\tPEPTIDEK
    PSM\tELVISLIVES\t1\tsp|P00002|TEST2_HUMAN\t1\tdb\tnull\t[, , dummy, 1]\t0.5\tnull\t20.1\t3\t372.21\t372.20\tms_run[2]:controllerType=0 controllerNumber=1 scan=202\tK\tA\t1\t10\t0.5\t0.4\t0\tELVISLIVES
    PSM\tDECYDECDR\t2\tDECOY_sp|P00003|TEST3_HUMAN\t1\tdb\tnull\t[, , dummy, 1]\t0.9\tnull\t30.0\t2\t555.00\t555.00\tms_run[1]:controllerType=0 controllerNumber=1 scan=303\tK\tR\t1\t10\t0.9\t0.8\t1\tDECYDECDR
    """
)


def test_mztab_parses_synthetic_fixture(tmp_path):
    fixture = tmp_path / "minimal.mzTab"
    fixture.write_text(_MINIMAL_MZTAB)

    records, file_map = iomztab.read_mztab(fixture, sample="syn")

    assert file_map == {0: "run_A", 1: "run_B"}
    assert len(records) == 2  # decoy dropped by default
    by_scan = {r.scan: r for r in records}
    assert by_scan[101].file_idx == 0
    assert by_scan[101].file_name == "run_A"
    assert by_scan[101].sequence == "PEPTIDEK"
    assert by_scan[101].charge == 2
    assert by_scan[202].file_idx == 1
    assert by_scan[202].file_name == "run_B"
    # peptide_mass is recomputed via accmass; just confirm it is finite and >0.
    assert all(r.peptide_mass > 0 for r in records)


def test_mztab_keeps_decoys_when_flag_off(tmp_path):
    fixture = tmp_path / "minimal.mzTab"
    fixture.write_text(_MINIMAL_MZTAB)
    records, _ = iomztab.read_mztab(fixture, sample="syn", drop_decoys=False)
    assert len(records) == 3


def test_mztab_preserves_zero_qvalue(tmp_path):
    """Regression: a q-value of exactly 0.0 must stay 0.0, not flip to 1.0.

    A previous version used ``float(value or default)``, which treats 0.0 as
    falsy and silently substitutes the missing-value default. On real quantms
    output that mis-classified ~80% of PSMs (the most confident ones) as
    completely non-confident.
    """
    body = _MINIMAL_MZTAB + (
        "PSM\tTPCNFER\t3\tsp|P00004|TEST4_HUMAN\t1\tdb\tnull\t"
        "[, , dummy, 1]\t0.0\tnull\t40.0\t2\t400.0\t400.0\t"
        "ms_run[1]:controllerType=0 controllerNumber=1 scan=404\t"
        "K\tR\t1\t7\t0.0\t0.0\t0\tTPCNFER\n"
    )
    fixture = tmp_path / "zeroq.mzTab"
    fixture.write_text(body)
    records, _ = iomztab.read_mztab(fixture, sample="syn")
    top = next(r for r in records if r.sequence == "TPCNFER")
    assert top.percolator_q_value == 0.0
    assert top.percolator_pep == 0.0
    assert top.percolator_score == 0.0


# --- writers -----------------------------------------------------------------


def test_provenance_hash_is_order_invariant():
    h1 = iowriters.hash_config({"a": 1, "b": [2, 3]})
    h2 = iowriters.hash_config({"b": [2, 3], "a": 1})
    assert h1 == h2


def test_write_tsv_emits_header(tmp_path):
    prov = iowriters.make_provenance({"mass_tol_ppm": 15}, id_source="x.txt")
    out = tmp_path / "result.tsv"
    iowriters.write_tsv(
        out,
        columns=["concat", "iso0"],
        rows=[{"concat": "AAAK_2", "iso0": 1.5e6}],
        provenance=prov,
    )
    lines = out.read_text().splitlines()
    assert lines[0].startswith("# riana ")
    assert any(ln.startswith("# git ") for ln in lines)
    assert any(ln.startswith("# config_hash ") for ln in lines)
    # Header row immediately after comments.
    header_idx = next(i for i, ln in enumerate(lines) if not ln.startswith("#"))
    assert lines[header_idx].split("\t") == ["concat", "iso0"]


def test_write_json_nests_provenance(tmp_path):
    prov = iowriters.make_provenance({"k": 1})
    out = tmp_path / "result.json"
    iowriters.write_json(out, {"rmse": 0.02}, provenance=prov)
    body = json.loads(out.read_text())
    assert "_provenance" in body
    assert body["_provenance"]["riana_version"]
    assert body["rmse"] == pytest.approx(0.02)
