"""Tests for ``riana.io.manifest`` (M6a stage-aware project index)."""

from __future__ import annotations

from pathlib import Path

import pytest

from riana.exceptions import DataError
from riana.io.manifest import (
    SCHEMA_VERSION,
    ManifestRow,
    append_manifest,
    read_manifest,
)
from riana.records import RunIdentity


def _turnover(sample, data_file, t, rep=1):
    return RunIdentity(
        experiment="exp", sample=sample, data_file=data_file,
        biological_replicate=rep, labeling_time=t, labeling_time_unit="day",
        condition="control", precursor_enrichment=0.046,
    )


def test_write_read_roundtrip(tmp_path):
    mf = tmp_path / "riana_manifest.tsv"
    rows = [
        ManifestRow("integrate", "out/t0_riana.txt", _turnover("s0", "t0", 0.0),
                    config_hash="abc", git_sha="g1"),
        ManifestRow("integrate", "out/t1_riana.txt", _turnover("s1", "t1", 1.0),
                    config_hash="abc", git_sha="g1"),
    ]
    append_manifest(mf, rows)

    back = read_manifest(mf)
    assert len(back) == 2
    r0 = next(r for r in back if r.output_path == "out/t0_riana.txt")
    assert r0.identity.sample == "s0"
    assert r0.identity.labeling_time == 0.0
    assert r0.identity.labeling_time_unit == "day"
    assert r0.identity.precursor_enrichment == pytest.approx(0.046)
    assert r0.identity.experiment_type == "turnover"
    assert r0.config_hash == "abc"


def test_none_optionals_roundtrip(tmp_path):
    """Calibration runs have labeling_time=None; turnover have mixing=None."""
    mf = tmp_path / "riana_manifest.tsv"
    cal = RunIdentity(experiment="c", sample="cal", data_file="c0",
                      mixing_proportion=0.5)
    append_manifest(mf, [ManifestRow("integrate", "out/c0.txt", cal)])
    back = read_manifest(mf)[0]
    assert back.identity.mixing_proportion == 0.5
    assert back.identity.labeling_time is None
    assert back.identity.precursor_enrichment is None
    assert back.identity.experiment_type == "calibration"


def test_upsert_is_idempotent(tmp_path):
    mf = tmp_path / "riana_manifest.tsv"
    row = ManifestRow("integrate", "out/t0.txt", _turnover("s0", "t0", 0.0),
                      config_hash="v1")
    append_manifest(mf, [row])
    # Re-run with the same (stage, output_path) but new content -> replace.
    append_manifest(mf, [ManifestRow("integrate", "out/t0.txt",
                                     _turnover("s0", "t0", 0.0), config_hash="v2")])
    back = read_manifest(mf)
    assert len(back) == 1
    assert back[0].config_hash == "v2"


def test_stage_filter(tmp_path):
    mf = tmp_path / "riana_manifest.tsv"
    ident = _turnover("s0", "t0", 0.0)
    append_manifest(mf, [
        ManifestRow("integrate", "out/t0.txt", ident),
        ManifestRow("fit", "out/fit.txt", ident),
    ])
    assert len(read_manifest(mf, stage="integrate")) == 1
    assert len(read_manifest(mf, stage="fit")) == 1
    assert len(read_manifest(mf)) == 2


def test_schema_header_written(tmp_path):
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, [ManifestRow("integrate", "o", _turnover("s", "d", 0.0))])
    first = mf.read_text().splitlines()[0]
    assert first == f"# manifest_schema {SCHEMA_VERSION}"


def test_rejects_non_manifest_file(tmp_path):
    bad = tmp_path / "not_a_manifest.tsv"
    bad.write_text("stage\toutput_path\nintegrate\tx\n")
    with pytest.raises(DataError, match="manifest"):
        read_manifest(bad)


def test_rejects_future_schema(tmp_path):
    fut = tmp_path / "future.tsv"
    fut.write_text("# manifest_schema 999\nstage\toutput_path\n")
    with pytest.raises(DataError, match="schema"):
        read_manifest(fut)


def test_invalid_stage_rejected():
    with pytest.raises(ValueError, match="stage"):
        ManifestRow("bogus", "o", _turnover("s", "d", 0.0))


def test_created_at_is_stamped_and_roundtrips(tmp_path):
    mf = tmp_path / "riana_manifest.tsv"
    row = ManifestRow("rollup", "out/riana_rollup_proteins.txt",
                      _turnover("s", "d", 0.0))
    assert row.created_at and "T" in row.created_at        # ISO timestamp set
    append_manifest(mf, [row])
    back = read_manifest(mf)[0]
    assert back.created_at == row.created_at               # preserved on read


def test_legacy_protein_stage_reads_as_rollup(tmp_path):
    """Pre-rename manifests (stage='protein') still load, mapped to 'rollup'."""
    mf = tmp_path / "riana_manifest.tsv"
    # Hand-write a legacy row (no created_at column, old 'protein' stage).
    mf.write_text(
        "# manifest_schema 1\n"
        "stage\toutput_path\texperiment\tsample\tdata_file\n"
        "protein\triana_protein.txt\texp\t\t\n"
    )
    rows = read_manifest(mf)
    assert len(rows) == 1
    assert rows[0].stage == "rollup"          # aliased
    assert rows[0].created_at == ""           # absent in the old file
    assert len(read_manifest(mf, stage="rollup")) == 1
