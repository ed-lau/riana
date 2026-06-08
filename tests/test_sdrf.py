"""Tests for ``riana.io.sdrf`` (M6a run-identity intake).

Uses tiny committed SDRF fixtures under ``tests/data/sdrf/`` (the working SDRFs
live under the gitignored ``data/`` tree, so the test suite carries its own).
Covers the documented column subset, the turnover-vs-calibration discriminant,
duplicate ``comment[modification parameters]`` handling, and the validation
errors.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from riana.exceptions import DataError
from riana.io.sdrf import read_sdrf

SDRF_DIR = Path("tests/data/sdrf")
TURNOVER = SDRF_DIR / "turnover.sdrf.tsv"
CALIBRATION = SDRF_DIR / "calibration.sdrf.tsv"
TWO_CONDITION = SDRF_DIR / "two_condition.sdrf.tsv"


def test_turnover_sdrf_parses_identity():
    t = read_sdrf(TURNOVER)
    assert t.experiment == "turnover"
    assert t.experiment_type == "turnover"
    assert t.acquisition == "DDA"
    assert len(t.runs) == 3

    r0 = t.runs[0]
    assert r0.sample == "mouse_t0"
    assert r0.data_file == "run_t0"  # ".raw" stripped
    assert r0.labeling_time == 0.0
    assert r0.labeling_time_unit == "day"
    assert r0.mixing_proportion is None
    assert r0.condition == "control"
    assert r0.precursor_enrichment == pytest.approx(0.046)
    assert r0.biological_replicate == 1
    assert r0.independent_value == 0.0

    # x-axis values across the series
    assert sorted(r.labeling_time for r in t.runs) == [0.0, 1.0, 2.0]


def test_sample_map_is_keyed_by_mzml_stem():
    t = read_sdrf(TURNOVER)
    assert set(t.sample_map) == {"run_t0", "run_t1", "run_t2"}
    assert t.sample_map["run_t1"].sample == "mouse_t1"


def test_duplicate_modification_columns_are_both_read():
    t = read_sdrf(TURNOVER)
    fixed = [(m.name, m.residue) for m in t.fixed_modifications]
    variable = [(m.name, m.residue) for m in t.variable_modifications]
    assert ("Carbamidomethyl", "C") in fixed
    assert ("Oxidation", "M") in variable


def test_calibration_sdrf_declares_calibration_type():
    t = read_sdrf(CALIBRATION)
    assert t.experiment_type == "calibration"
    r = t.runs[1]
    assert r.mixing_proportion == 0.5
    assert r.labeling_time is None
    assert r.independent_value == 0.5
    assert sorted(x.mixing_proportion for x in t.runs) == [0.0, 0.5, 1.0]


def test_two_condition_sdrf_groups_and_replicates():
    t = read_sdrf(TWO_CONDITION)
    assert t.acquisition == "DIA"
    conditions = sorted({r.condition for r in t.runs})
    assert conditions == ["control", "knockout"]
    # control has a t14 biorep-1 and biorep-2 (replicate animals at one time).
    ctrl_t14 = [
        r for r in t.runs
        if r.condition == "control" and r.labeling_time == 14.0
    ]
    assert sorted(r.biological_replicate for r in ctrl_t14) == [1, 2]


def test_experiment_label_defaults_to_stem(tmp_path):
    src = (TWO_CONDITION).read_text()
    p = tmp_path / "myproject.sdrf.tsv"
    p.write_text(src)
    t = read_sdrf(p)
    assert t.experiment == "myproject"  # trailing ".sdrf" stripped
    assert read_sdrf(p, experiment="override").experiment == "override"


# --- validation --------------------------------------------------------------


def _write(tmp_path, text: str) -> Path:
    p = tmp_path / "bad.sdrf.tsv"
    p.write_text(textwrap.dedent(text))
    return p


def test_error_when_neither_experiment_column(tmp_path):
    p = _write(tmp_path, """\
        source name\tcomment[data file]
        s0\trun0.raw
    """)
    with pytest.raises(DataError, match="exactly one"):
        read_sdrf(p)


def test_error_when_both_experiment_columns(tmp_path):
    p = _write(tmp_path, """\
        source name\tcomment[data file]\tcharacteristics[labeling time]\tcharacteristics[mixing proportion]
        s0\trun0.raw\t0 day\t0
    """)
    with pytest.raises(DataError, match="exactly one"):
        read_sdrf(p)


def test_error_on_missing_required_column(tmp_path):
    p = _write(tmp_path, """\
        source name\tcharacteristics[labeling time]
        s0\t0 day
    """)
    with pytest.raises(DataError, match="comment\\[data file\\]"):
        read_sdrf(p)


def test_error_on_duplicate_data_file(tmp_path):
    p = _write(tmp_path, """\
        source name\tcomment[data file]\tcharacteristics[labeling time]
        s0\trun0.raw\t0 day
        s1\trun0.raw\t1 day
    """)
    with pytest.raises(DataError, match="unique per run"):
        read_sdrf(p)


def test_error_on_non_numeric_labeling_time(tmp_path):
    p = _write(tmp_path, """\
        source name\tcomment[data file]\tcharacteristics[labeling time]
        s0\trun0.raw\tearly
    """)
    with pytest.raises(DataError, match="labeling time"):
        read_sdrf(p)


def test_error_on_mixed_acquisition(tmp_path):
    p = _write(tmp_path, """\
        source name\tcomment[data file]\tcharacteristics[labeling time]\tcomment[proteomics data acquisition method]
        s0\trun0.raw\t0 day\tNT=data-dependent acquisition;AC=MS:1003221
        s1\trun1.raw\t1 day\tNT=data-independent acquisition;AC=MS:1003215
    """)
    with pytest.raises(DataError, match="acquisition"):
        read_sdrf(p)
