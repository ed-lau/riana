"""Tests for the DIA-NN parquet intake (M6b, :mod:`riana.io.diann`) and the
RT→MS1-scan resolver it depends on (:func:`riana.core.integration.resolve_rt_anchored_scans`).

Synthetic fixtures only — the real cardiac DIA mzMLs/parquet are gigabyte-scale
and live outside the repo, so CI builds a tiny parquet in a tempdir and a stub
mzML index.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from riana.core.integration import resolve_rt_anchored_scans
from riana.exceptions import DataError
from riana.io.diann import read_diann
from riana.records import PSMRecord, RunIdentity


def _identity(stem: str, sample: str, biorep: int, t: float) -> RunIdentity:
    return RunIdentity(
        experiment="dia_exp", sample=sample, data_file=stem,
        biological_replicate=biorep, labeling_time=t, labeling_time_unit="days",
        condition="control", acquisition="DIA", precursor_enrichment=0.046,
    )


def _write_parquet(path: Path, rows: list[dict]) -> None:
    cols = ["Run", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
            "Precursor.Mz", "Decoy", "Q.Value", "PEP", "RT", "Protein.Ids"]
    pd.DataFrame(rows, columns=cols).to_parquet(path, index=False)


@pytest.fixture()
def sample_map():
    return {
        "843_LV": _identity("843_LV", "LV_843", 1, 14.0),
        "930_LV": _identity("930_LV", "LV_930", 1, 3.0),
    }


def _row(run, seq, mod_seq, z, decoy=0, q=0.001, rt=40.0, prot="Q3TW96"):
    return dict(zip(
        ["Run", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
         "Precursor.Mz", "Decoy", "Q.Value", "PEP", "RT", "Protein.Ids"],
        [run, seq, mod_seq, z, 800.0, decoy, q, q, rt, prot],
    ))


_SITE_COLS = ["Run", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
              "Precursor.Mz", "Decoy", "Q.Value", "PEP", "RT", "Protein.Ids",
              "Protein.Sites", "PTM.Site.Confidence"]


def _site_row(run, seq, mod_seq, z, sites, conf, prot="Q3TW96", rt=40.0):
    """A row with DIA-NN's PTM-aware columns (Protein.Sites + confidence)."""
    return dict(zip(_SITE_COLS,
                    [run, seq, mod_seq, z, 800.0, 0, 0.001, 0.001, rt, prot,
                     sites, conf]))


def _write_parquet_sites(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows, columns=_SITE_COLS).to_parquet(path, index=False)


def test_read_diann_basic(tmp_path, sample_map):
    p = tmp_path / "diann_report.parquet"
    _write_parquet(p, [
        _row("843_LV", "AAAAGALAPGPLPDLAAR", "AAAAGALAPGPLPDLAAR", 2, rt=56.9),
        _row("930_LV", "AAAAGALAPGPLPDLAAR", "AAAAGALAPGPLPDLAAR", 2, rt=50.9),
        _row("843_LV", "SAMPLEDPEPTIDEK", "SAMPLEDPEPTIDEK", 2, rt=30.0,
             prot="A0A0G2JDV3;Q61107;Q91Z40"),
    ])
    records, file_index_map = read_diann(p, sample_map)

    assert len(records) == 3
    assert file_index_map == {0: "843_LV", 1: "930_LV"}  # sorted Run order

    r = next(x for x in records if x.file_name == "930_LV")
    assert r.scan == -1  # DIA: unresolved sentinel
    assert r.charge == 2
    assert r.sequence == "AAAAGALAPGPLPDLAAR"
    assert r.peptide_mass > 0
    assert r.sample == "LV_930"
    assert r.identity is not None and r.identity.acquisition == "DIA"
    # RT carried in seconds (mzTab convention): 50.9 min → 3054 s.
    assert r.retention_time == pytest.approx(50.9 * 60.0)
    assert r.percolator_q_value == pytest.approx(0.001)

    # ';'-joined Protein.Ids normalized to ','.
    multi = next(x for x in records if x.sequence == "SAMPLEDPEPTIDEK")
    assert multi.protein_id == "A0A0G2JDV3,Q61107,Q91Z40"


def test_read_diann_encodes_supported_mods_drops_unsupported_and_decoys(tmp_path, sample_map):
    p = tmp_path / "report.parquet"
    _write_parquet(p, [
        # kept bare: fixed Carbamidomethyl folded per-cysteine
        _row("843_LV", "AAADEWTTCTPPSGLQGK", "AAADEWTTC(UniMod:4)TPPSGLQGK", 2),
        # kept + encoded: Met-Ox (tier 1b) integrates separately, merges at fit
        _row("843_LV", "MLSEDQVK", "M(UniMod:35)LSEDQVK", 2),
        # dropped: deamidation is not in the supported set
        _row("843_LV", "NLSEDQVK", "N(UniMod:7)LSEDQVK", 2),
        # dropped: decoy
        _row("843_LV", "DECOYPEPTIDEK", "DECOYPEPTIDEK", 2, decoy=1),
    ])
    records, _ = read_diann(p, sample_map)
    seqs = {r.sequence for r in records}
    assert seqs == {"AAADEWTTCTPPSGLQGK", "M[UNIMOD:35]LSEDQVK"}

    # opting out keeps the unsupported-mod form too (still drops the decoy)
    records2, _ = read_diann(p, sample_map, drop_variable_mods=False)
    assert {r.sequence for r in records2} == {
        "AAADEWTTCTPPSGLQGK", "M[UNIMOD:35]LSEDQVK", "NLSEDQVK"}


def test_read_diann_phospho_proteoform_sites(tmp_path, sample_map):
    """M7 Stage B on the DIA path: DIA-NN's localized Protein.Sites -> the
    biological-mod proteoform suffix (pS###), gated on PTM.Site.Confidence and
    filtered to phospho (the constitutive Carbamidomethyl C site is dropped)."""
    p = tmp_path / "report.parquet"
    _write_parquet_sites(p, [
        # confident single phospho -> pS9
        _site_row("843_LV", "AAADGDDSLYPIAVLIDELR",
                  "AAADGDDS(UniMod:21)LYPIAVLIDELR", 2, "[Q76MZ3:S9]", 1.0,
                  prot="Q76MZ3"),
        # phospho + carbamidomethyl: the C100 site is constitutive -> dropped, pS98
        _site_row("843_LV", "ADLSNCLYKDMPATIDSVFAR",
                  "ADLS(UniMod:21)NC(UniMod:4)LYKDMPATIDSVFAR", 2,
                  "[Q0II04:S98,C100]", 0.93, prot="Q0II04"),
        # low confidence (< 0.75) -> folds into the bare protein
        _site_row("930_LV", "AAAESSAIQSISHV",
                  "AAAES(UniMod:21)SAIQSISHV", 2, "[P14404:S1223]", 0.60,
                  prot="P14404"),
        # two confident phospho sites -> pS1332_pS1333 (sorted, _-joined)
        _site_row("930_LV", "AIRLGKELSSAESQLHDTQELLQEETR",
                  "AIRLGKELS(UniMod:21)S(UniMod:21)AESQLHDTQELLQEETR", 3,
                  "[Q6URW6:S1332,S1333]", 0.99, prot="Q6URW6"),
    ])
    records, _ = read_diann(p, sample_map)
    by_seq = {r.sequence: r.mod_sites for r in records}
    assert by_seq["AAADGDDS[UNIMOD:21]LYPIAVLIDELR"] == "pS9"
    assert by_seq["ADLS[UNIMOD:21]NCLYKDMPATIDSVFAR"] == "pS98"   # C100 dropped
    assert by_seq["AAAES[UNIMOD:21]SAIQSISHV"] == ""              # gated -> bare
    assert by_seq["AIRLGKELS[UNIMOD:21]S[UNIMOD:21]AESQLHDTQELLQEETR"] == \
        "pS1332_pS1333"

    # A custom gate keeps the 0.60 row (proves the threshold is honoured).
    recs2, _ = read_diann(p, sample_map, min_site_confidence=0.5)
    assert {r.mod_sites for r in recs2 if r.sequence ==
            "AAAES[UNIMOD:21]SAIQSISHV"} == {"pS1223"}


def test_read_diann_metox_excluded_from_phospho_key(tmp_path, sample_map):
    """A Met-Ox + phospho peptidoform keys on the phospho site ONLY — Met-Ox (a
    CHEMICAL_MODS member, merged at fit) never defines a proteoform, just like
    Carbamidomethyl. (The two real searches don't co-occur, so this is synthetic.)"""
    p = tmp_path / "report.parquet"
    _write_parquet_sites(p, [
        _site_row("843_LV", "MLSPEPTIDESK", "M(UniMod:35)LS(UniMod:21)PEPTIDESK",
                  2, "[Q9Z1P6:M87,S89]", 1.0, prot="Q9Z1P6"),
    ])
    records, _ = read_diann(p, sample_map)
    assert len(records) == 1
    assert records[0].sequence == "M[UNIMOD:35]LS[UNIMOD:21]PEPTIDESK"  # Ox kept
    assert records[0].mod_sites == "pS89"  # M87 (Met-Ox) excluded from the key


def test_read_diann_no_site_columns_still_loads(tmp_path, sample_map):
    """A report without the PTM-aware columns (older / no-mod search) loads with
    empty mod_sites — graceful degradation, no hard failure."""
    p = tmp_path / "report.parquet"
    _write_parquet(p, [_row("843_LV", "SAMPLEDPEPTIDEK", "SAMPLEDPEPTIDEK", 2)])
    records, _ = read_diann(p, sample_map)
    assert len(records) == 1 and records[0].mod_sites == ""


def test_read_diann_unmatched_run_errors(tmp_path, sample_map):
    p = tmp_path / "report.parquet"
    _write_parquet(p, [_row("UNKNOWN_RUN", "PEPTIDEK", "PEPTIDEK", 2)])
    with pytest.raises(DataError, match="not in the SDRF"):
        read_diann(p, sample_map)


class _StubMzML:
    """Minimal IndexedMzML stand-in for the resolver (scan/rt index + path)."""

    def __init__(self, scans, rts_min):
        self.scan_idx = np.asarray(scans, dtype=np.int64)
        self.rt_idx = np.asarray(rts_min, dtype=np.float64)
        self.path = Path("stub_run.mzML")


def _psm(scan, rt_s):
    return PSMRecord(
        scan=scan, charge=2, sequence="PEPTIDEK", peptide_mass=900.0,
        sample="s", retention_time=rt_s,
    )


def test_resolve_rt_anchored_scans_maps_to_nearest_ms1():
    # MS1 scans 1, 21, 41 at RT 10.0, 20.0, 30.0 min.
    mz = _StubMzML([1, 21, 41], [10.0, 20.0, 30.0])
    psms = [
        _psm(-1, 20.4 * 60),  # nearest 20.0 → scan 21
        _psm(-1, 29.9 * 60),  # nearest 30.0 → scan 41
        _psm(-1, 10.1 * 60),  # nearest 10.0 → scan 1
    ]
    out = resolve_rt_anchored_scans(psms, mz)
    assert [p.scan for p in out] == [21, 41, 1]


def test_resolve_rt_anchored_scans_noop_for_dda():
    mz = _StubMzML([1, 21, 41], [10.0, 20.0, 30.0])
    psms = [_psm(17, 20.0 * 60)]  # already has a real scan
    out = resolve_rt_anchored_scans(psms, mz)
    assert out[0].scan == 17  # untouched


def test_resolve_rt_anchored_scans_wrong_mzml_raises():
    # mzML spans 10–30 min; the DIA RTs are all way outside → wrong pairing.
    mz = _StubMzML([1, 21, 41], [10.0, 20.0, 30.0])
    psms = [_psm(-1, 80.0 * 60), _psm(-1, 90.0 * 60), _psm(-1, 100.0 * 60)]
    with pytest.raises(DataError, match="wrong mzML"):
        resolve_rt_anchored_scans(psms, mz)
