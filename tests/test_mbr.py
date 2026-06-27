"""Unit tests for match-between-runs donor assembly + RT alignment (core/mbr)."""

from __future__ import annotations

import pytest

from riana.config import IntegrationConfig
from riana.core.mbr import augment_with_mbr
from riana.records import PSMRecord, RunIdentity


def _ident(condition: str, data_file: str, experiment: str = "exp",
           fraction: int = 1) -> RunIdentity:
    return RunIdentity(
        experiment=experiment, sample=data_file, data_file=data_file,
        condition=condition, fraction=fraction,
    )


def _psm(seq, charge, q, rt, file_idx, ident) -> PSMRecord:
    return PSMRecord(
        scan=file_idx * 1000 + 1, charge=charge, sequence=seq,
        peptide_mass=1000.0 + len(seq), sample=ident.sample, file_idx=file_idx,
        file_name=ident.data_file, percolator_q_value=q, retention_time=rt,
        identity=ident, mod_sites="", protein_id="P1",
    )


# Three runs of one curve group (same experiment+condition), distinct file_idx.
R0 = _ident("ctrl", "run0")
R1 = _ident("ctrl", "run1")
R2 = _ident("ctrl", "run2")

# A: located in all 3 runs (no acceptor). B: in runs 0,1 (→ transfer to run 2).
# C: in run 0 only (< 2 donor runs → no transfer). Run RTs carry a per-run shift:
# run1 ≈ +10, run2 ≈ +22 vs run0, so the consensus/offset alignment is testable.
_GROUP = [
    _psm("PEPA", 2, 0.0, 100.0, 0, R0),
    _psm("PEPA", 2, 0.0, 110.0, 1, R1),
    _psm("PEPA", 2, 0.0, 122.0, 2, R2),
    _psm("PEPB", 2, 0.0, 200.0, 0, R0),
    _psm("PEPB", 2, 0.0, 210.0, 1, R1),
    _psm("PEPC", 2, 0.0, 300.0, 0, R0),
]

_CFG = IntegrationConfig(mbr=True, mbr_min_donor_runs=2, mbr_donor_q=1e-2)


def _transfers(records):
    return [p for p in records if p.evidence == "mbr"]


def test_transfers_only_to_missing_runs():
    out = augment_with_mbr(_GROUP, _CFG)
    t = _transfers(out)
    # only PEPB → run 2 qualifies (A is everywhere, C has one donor run)
    assert len(t) == 1
    rec = t[0]
    assert rec.concat == "PEPB_2"
    assert rec.file_idx == 2


def test_mbr_restricted_to_winner_fraction():
    """A precursor split across fractions is matched only in the fraction with the
    most IDs (the winner); holes in its minor fraction are left untouched."""
    f5 = {r: _ident("ctrl", f"f5_r{r}", fraction=5) for r in range(4)}
    f6 = {r: _ident("ctrl", f"f6_r{r}", fraction=6) for r in range(4)}
    psms = []
    # PEPA anchors every run of both fractions, so all runs exist (no PEPA holes).
    for r in range(4):
        psms.append(_psm("PEPA", 2, 0.0, 100.0, r, f5[r]))
        psms.append(_psm("PEPA", 2, 0.0, 100.0, 10 + r, f6[r]))
    # PEPB: located in 3 runs of fraction 5 (hole at r3) and 2 of fraction 6
    # (holes at r2, r3). Unrestricted MBR would fill all three holes; the winner
    # policy (fraction 5 has more IDs) restricts it to fraction 5's single hole.
    for r in range(3):
        psms.append(_psm("PEPB", 2, 0.0, 200.0, r, f5[r]))
    for r in range(2):
        psms.append(_psm("PEPB", 2, 0.0, 200.0, 10 + r, f6[r]))
    b = [p for p in _transfers(augment_with_mbr(psms, _CFG))
         if p.concat == "PEPB_2"]
    assert len(b) == 1
    assert b[0].identity.fraction == 5
    assert b[0].file_idx == 3


def test_no_transfer_for_singleton_donor():
    out = augment_with_mbr(_GROUP, _CFG)
    assert all(p.concat != "PEPC_2" for p in _transfers(out))


def test_no_transfer_when_present_in_all_runs():
    out = augment_with_mbr(_GROUP, _CFG)
    assert all(p.concat != "PEPA_2" for p in _transfers(out))


def test_transfer_record_fields():
    rec = _transfers(augment_with_mbr(_GROUP, _CFG))[0]
    # RT-anchored sentinel + flag
    assert rec.scan == -1
    assert rec.evidence == "mbr"
    assert rec.percolator_q_value == 0.0
    # precursor identity from the donor template (PEPB)
    assert rec.sequence == "PEPB"
    assert rec.charge == 2
    assert rec.peptide_mass == pytest.approx(1000.0 + len("PEPB"))
    # run identity from the acceptor (run 2)
    assert rec.identity is R2
    assert rec.file_name == "run2"
    assert rec.sample == "run2"


def test_transfer_rt_uses_consensus_plus_run_offset():
    # consensus_rt[B] = median(200, 210) = 205; run-2 offset from its only shared
    # precursor A: 122 - consensus_rt[A]=median(100,110,122)=110 → +12.
    rec = _transfers(augment_with_mbr(_GROUP, _CFG))[0]
    assert rec.retention_time == pytest.approx(205.0 + 12.0)


def test_min_donor_runs_knob():
    cfg3 = IntegrationConfig(mbr=True, mbr_min_donor_runs=3)
    # PEPB now has only 2 donor runs (< 3) → no transfer at all.
    assert _transfers(augment_with_mbr(_GROUP, cfg3)) == []


def test_transfers_do_not_cross_groups():
    other = [
        _psm("PEPB", 2, 0.0, 200.0, 10, _ident("treat", "run10")),
        _psm("PEPB", 2, 0.0, 210.0, 11, _ident("treat", "run11")),
    ]
    out = augment_with_mbr(_GROUP + other, _CFG)
    # ctrl PEPB → run2 (1), treat PEPB present in both its runs → 0. Still just 1.
    assert len(_transfers(out)) == 1
    assert _transfers(out)[0].identity is R2


def test_records_without_identity_pass_through_untouched():
    no_id = PSMRecord(scan=5, charge=2, sequence="ZZ", peptide_mass=900.0,
                      sample="s", file_idx=0, percolator_q_value=0.0)
    out = augment_with_mbr([no_id, *_GROUP], _CFG)
    assert no_id in out
    assert len(_transfers(out)) == 1


def test_low_quality_donor_excluded():
    # PEPB identified at q=0.5 in run 1 → only 1 confident donor run → no transfer.
    group = [
        _psm("PEPA", 2, 0.0, 100.0, 0, R0),
        _psm("PEPA", 2, 0.0, 110.0, 1, R1),
        _psm("PEPB", 2, 0.0, 200.0, 0, R0),
        _psm("PEPB", 2, 0.5, 210.0, 1, R1),
    ]
    assert _transfers(augment_with_mbr(group, _CFG)) == []
