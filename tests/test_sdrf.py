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
from riana.io.sdrf import _parse_mass_tol_ppm, read_sdrf

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


def test_precursor_mass_tolerance_read_as_ppm():
    """comment[precursor mass tolerance] feeds the integration window (rec4
    reversed): the search tolerance IS the right window for centroid mzML."""
    t = read_sdrf(TURNOVER)
    assert t.precursor_mass_tol_ppm == 10.0


@pytest.mark.parametrize(
    "value, expected",
    [
        ("10 ppm", 10.0),
        ("10ppm", 10.0),
        ("10", 10.0),        # bare number assumed ppm
        ("4.5 ppm", 4.5),
        ("0.02 Da", None),   # Da window can't be a ppm tolerance
        ("not applicable", None),
        ("", None),
    ],
)
def test_parse_mass_tol_ppm_units(value, expected):
    assert _parse_mass_tol_ppm(value) == expected


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


def test_isobaric_tmt_collapses_channels_per_file(tmp_path):
    """A TMT SDRF (comment[label]=TMT*) collapses its per-channel rows to one run
    per data file — the multiplexed samples share one MS1 cluster. A file whose
    channels are all one treatment keeps that clean condition (so split-batch
    designs flow into the linear-simple 2-sample Δk); a file pooling treatments
    carries the combined ``|``-joined label."""
    cols = ["source name", "characteristics[organism]",
            "characteristics[biological replicate]",
            "characteristics[labeling time]", "comment[fraction identifier]",
            "comment[label]", "comment[data file]",
            "comment[proteomics data acquisition method]",
            "characteristics[precursor enrichment]", "factor value[treatment]"]
    rows = [
        ["s1", "Homo sapiens", "1", "24 hours", "1", "TMT126", "fileA.mzML", "DDA", "0.06", "control"],
        ["s2", "Homo sapiens", "1", "24 hours", "1", "TMT127N", "fileA.mzML", "DDA", "0.06", "control"],
        ["s3", "Homo sapiens", "1", "24 hours", "1", "TMT128N", "fileA.mzML", "DDA", "0.06", "control"],
        ["s4", "Homo sapiens", "2", "24 hours", "1", "TMT126", "fileB.mzML", "DDA", "0.06", "control"],
        ["s5", "Homo sapiens", "2", "24 hours", "1", "TMT127N", "fileB.mzML", "DDA", "0.06", "nocodazole"],
        ["s6", "Homo sapiens", "2", "24 hours", "1", "TMT128N", "fileB.mzML", "DDA", "0.06", "control"],
    ]
    p = tmp_path / "tmt.sdrf.tsv"
    p.write_text("\n".join("\t".join(r) for r in [cols, *rows]) + "\n")

    t = read_sdrf(p)
    assert len(t.runs) == 2, "6 channel rows across 2 files -> 2 collapsed runs"
    by_file = {r.data_file: r for r in t.runs}
    assert by_file["fileA"].condition == "control"            # single-treatment plex
    assert by_file["fileB"].condition == "control|nocodazole"  # pooled plex
    assert by_file["fileA"].biological_replicate == 1
    assert by_file["fileB"].biological_replicate == 2
    assert by_file["fileA"].sample == "fileA"  # per-channel source name -> file stem


def test_isobaric_channels_disagreeing_on_labeling_time_raise(tmp_path):
    """Isobaric collapse keeps one channel's file-level fields — labeling_time is the
    fit x-axis — so it must reject an SDRF whose channels within a file disagree on it
    (a typo, or an unsupported pulsed-labeling-in-channels layout) rather than silently
    stamp an arbitrary channel's time."""
    cols = ["source name", "characteristics[organism]",
            "characteristics[biological replicate]",
            "characteristics[labeling time]", "comment[fraction identifier]",
            "comment[label]", "comment[data file]",
            "comment[proteomics data acquisition method]",
            "characteristics[precursor enrichment]", "factor value[treatment]"]
    rows = [
        ["s1", "Homo sapiens", "1", "24 hours", "1", "TMT126", "fileA.mzML", "DDA", "0.06", "control"],
        ["s2", "Homo sapiens", "1", "48 hours", "1", "TMT127N", "fileA.mzML", "DDA", "0.06", "control"],
    ]
    p = tmp_path / "bad.sdrf.tsv"
    p.write_text("\n".join("\t".join(r) for r in [cols, *rows]) + "\n")
    with pytest.raises(DataError, match="labeling_time"):
        read_sdrf(p)


def test_experiment_label_defaults_to_stem(tmp_path):
    src = (TWO_CONDITION).read_text()
    p = tmp_path / "myproject.sdrf.tsv"
    p.write_text(src)
    t = read_sdrf(p)
    assert t.experiment == "myproject"  # trailing ".sdrf" stripped
    assert read_sdrf(p, experiment="override").experiment == "override"


def _stratified_sdrf(tmp_path) -> Path:
    """A 2-tissue × 2-condition turnover SDRF for --experiment-column tests."""
    cols = ["source name", "characteristics[organism]",
            "characteristics[organism part]",
            "characteristics[biological replicate]",
            "characteristics[labeling time]", "comment[data file]",
            "comment[proteomics data acquisition method]",
            "characteristics[precursor enrichment]", "factor value[condition]"]
    rows = [
        ["s1", "Mus musculus", "Left atrium", "1", "1 day", "a1.mzML", "DIA", "0.046", "control"],
        ["s2", "Mus musculus", "Left atrium", "1", "1 day", "a2.mzML", "DIA", "0.046", "knockout"],
        ["s3", "Mus musculus", "Left ventricle", "1", "1 day", "b1.mzML", "DIA", "0.046", "control"],
        ["s4", "Mus musculus", "Left ventricle", "1", "1 day", "b2.mzML", "DIA", "0.046", "knockout"],
    ]
    p = tmp_path / "strat.sdrf.tsv"
    p.write_text("\n".join("\t".join(r) for r in [cols, *rows]) + "\n")
    return p


def test_experiment_column_stratifies_runs_per_row(tmp_path):
    p = _stratified_sdrf(tmp_path)
    # default: one experiment (the stem) on every run
    assert {r.experiment for r in read_sdrf(p).runs} == {"strat"}
    # --experiment-column: per-row experiment from the named column
    t = read_sdrf(p, experiment_column="characteristics[organism part]")
    by_file = {r.data_file: r.experiment for r in t.runs}
    assert by_file == {"a1": "Left atrium", "a2": "Left atrium",
                       "b1": "Left ventricle", "b2": "Left ventricle"}
    # condition (the contrast axis) is unaffected by the stratifier
    assert {r.condition for r in t.runs} == {"control", "knockout"}


def test_experiment_column_missing_raises(tmp_path):
    p = _stratified_sdrf(tmp_path)
    with pytest.raises(DataError, match="not found"):
        read_sdrf(p, experiment_column="characteristics[nope]")


def test_experiment_column_empty_cell_raises(tmp_path):
    cols = ["source name", "characteristics[organism]",
            "characteristics[organism part]",
            "characteristics[biological replicate]",
            "characteristics[labeling time]", "comment[data file]",
            "comment[proteomics data acquisition method]",
            "characteristics[precursor enrichment]", "factor value[condition]"]
    rows = [
        ["s1", "Mus musculus", "Left atrium", "1", "1 day", "a1.mzML", "DIA", "0.046", "control"],
        ["s2", "Mus musculus", "not available", "1", "1 day", "a2.mzML", "DIA", "0.046", "control"],
    ]
    p = tmp_path / "gap.sdrf.tsv"
    p.write_text("\n".join("\t".join(r) for r in [cols, *rows]) + "\n")
    with pytest.raises(DataError, match="empty"):
        read_sdrf(p, experiment_column="characteristics[organism part]")


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
