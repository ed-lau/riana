"""Tests for ``riana.core.pipeline`` (M6a shared orchestration).

Three things to nail down:

1. **Recombination** groups runs into curves by ``(experiment, condition)`` and
   merges fractions of the same ``(biorep, timepoint)`` at peptide level, keeping
   different biological replicates as independent points.
2. **Parity**: the manifest/identity fit path recovers the *same* ``k_deg`` as
   the legacy positional path on trivial identity — i.e. the new spine doesn't
   change the science, only where the timepoint comes from.
3. **integrate_project** runs end-to-end on the real BSA mzML via a synthetic
   SDRF + mzTab, writing one ``<stem>_riana.txt`` per run with identity in the
   header plus a manifest row.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from riana.algorithms.isotope_dist import (
    _get_final_env,
    _get_init_env,
    clear_envelope_cache,
)
from riana.algorithms.mass_calc import calculate_ion_mz
from riana.config import FitConfig, IntegrationConfig
from riana.core.fitting import fit_run
from riana.core.pipeline import (
    fit_project,
    identity_to_extra,
    integrate_project,
    recombine_for_fit,
)
from riana.io.manifest import ManifestRow, append_manifest, read_manifest
from riana.io.sdrf import read_sdrf
from riana.records import RunIdentity

SAMPLE1 = Path("tests/data/sample1")
BSA_MZML = SAMPLE1 / "20180216_BSA.mzML.gz"

# --- synthetic D2O envelope series (mirrors tests/test_fitting.py) ------------
_TEST_PEPTIDES = [("VAPEPTIDEK", 2), ("LSHGLNVR", 2), ("AEFVEVTK", 2),
                  ("STDNAATR", 2), ("HVELFK", 2)]
_K_DEG_TRUE = 0.5
_PROPORTIONS = np.array([0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875])
_TIMES = -np.log(1.0 - _PROPORTIONS) / _K_DEG_TRUE


def _coeffs():
    avg_len = sum(len(s) for s, _ in _TEST_PEPTIDES) / len(_TEST_PEPTIDES)
    per_aa = 8 / avg_len
    return {aa: per_aa for aa in "ACDEFGHIKLMNPQRSTVWY"}


def _spep_by_seq(coeffs):
    return {s: max(1, int(round(sum(coeffs.get(c, 0.0) for c in s))))
            for s, _ in _TEST_PEPTIDES}


def _make_timepoint_dfs(coeffs):
    """One integrate-output DataFrame per timepoint (sample encodes the time)."""
    spep_by = _spep_by_seq(coeffs)
    clear_envelope_cache()
    dfs = []
    for ti, prop in zip(_TIMES, _PROPORTIONS):
        rows = []
        for k, (seq, charge) in enumerate(_TEST_PEPTIDES):
            pep_mass = calculate_ion_mz(seq)
            init = _get_init_env(seq, pep_mass, n=6)
            final = _get_final_env(seq, pep_mass, spep_by[seq], ria_max=0.06, n=6)
            mix = (1 - prop) * (init / init.sum()) + prop * (final / final.sum())
            scaled = mix * 1e6
            rows.append({
                "file_idx": 0, "scan": 1000 + k, "charge": charge,
                "concat": f"{seq}_{charge}", "sequence": seq,
                "sample": f"time{ti:.6f}", "percolator q-value": 1e-4,
                "protein id": f"sp|P{k}|TEST",
                **{f"iso{j}": scaled[j] for j in range(6)},
            })
        dfs.append(pd.DataFrame(rows))
    return dfs


# --- recombination -----------------------------------------------------------


def _integrate_rows_from_dfs(tmp_path, dfs, *, condition="control"):
    """Write each timepoint df as a per-run _riana.txt and return manifest rows."""
    rows = []
    for ti, df in zip(_TIMES, dfs):
        stem = f"run_t{ti:.4f}"
        path = tmp_path / f"{stem}_riana.txt"
        df.to_csv(path, sep="\t", index=False)
        ident = RunIdentity(
            experiment="syn", sample=stem, data_file=stem,
            labeling_time=float(ti), labeling_time_unit="au", condition=condition,
        )
        rows.append(ManifestRow("integrate", str(path), ident))
    return rows


def _calibration_rows_from_dfs(tmp_path, dfs, *, condition="control"):
    """Per-run _riana.txt + manifest rows for a mixing-calibration series — each
    run carries ``mixing_proportion`` (the heavy fraction) instead of a labeling
    time, so RunIdentity.experiment_type == 'calibration'."""
    rows = []
    for prop, df in zip(_PROPORTIONS, dfs):
        stem = f"calrun_f{prop:.4f}"
        path = tmp_path / f"{stem}_riana.txt"
        df.to_csv(path, sep="\t", index=False)
        ident = RunIdentity(
            experiment="syn", sample=stem, data_file=stem,
            mixing_proportion=float(prop), condition=condition,
        )
        rows.append(ManifestRow("integrate", str(path), ident))
    return rows


def _make_calibration_dfs(coeffs):
    """Mixing-calibration envelopes (FS == proportion); proportion comes from the
    manifest identity, so the sample string is just a label here."""
    spep_by = _spep_by_seq(coeffs)
    clear_envelope_cache()
    dfs = []
    for prop in _PROPORTIONS:
        rows = []
        for k, (seq, charge) in enumerate(_TEST_PEPTIDES):
            pep_mass = calculate_ion_mz(seq)
            init = _get_init_env(seq, pep_mass, n=6)
            final = _get_final_env(seq, pep_mass, spep_by[seq], ria_max=0.06, n=6)
            mix = (1 - prop) * (init / init.sum()) + prop * (final / final.sum())
            scaled = mix * 1e6
            rows.append({
                "file_idx": 0, "charge": charge, "concat": f"{seq}_{charge}",
                "sequence": seq, "sample": f"cal{prop:.4f}",
                "percolator q-value": 1e-4, "protein id": f"sp|P{k}|TEST",
                **{f"iso{j}": scaled[j] for j in range(6)},
            })
        dfs.append(pd.DataFrame(rows))
    return dfs


def test_fit_project_auto_dispatches_calibration_model(tmp_path):
    """A manifest whose runs carry mixing_proportion (experiment_type==calibration)
    is fit with the calibration recovery line — even when config.model is the
    kinetic default — so the recovered slope ≈ 1 and R² ≈ 1, and a default-model
    run matches an explicit --model calibration run."""
    coeffs = _coeffs()
    rows = _calibration_rows_from_dfs(tmp_path, _make_calibration_dfs(coeffs))
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)

    # config.model is the kinetic default; the manifest's calibration identity
    # must override it (the decided experiment-type dispatch).
    auto = fit_project(FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                                 ria_max=0.06), mf, coeffs, n_boot=0, random_state=1)
    slopes = auto["k_deg"].dropna().to_numpy()
    np.testing.assert_allclose(slopes, 1.0, atol=0.05)
    assert (auto["R_squared"].dropna().to_numpy() > 0.999).all()

    # An explicit --model calibration produces the identical result — proof the
    # default-model run dispatched to calibration, not the simple exponential.
    explicit = fit_project(FitConfig(model="calibration", label="hw", q_value=0.05,
                                     depth=3, ria_max=0.06), mf, coeffs,
                           n_boot=0, random_state=1)
    np.testing.assert_allclose(auto.sort_index()["k_deg"].to_numpy(dtype=float),
                               explicit.sort_index()["k_deg"].to_numpy(dtype=float),
                               rtol=0, atol=1e-9, equal_nan=True)


def test_fit_project_resolves_ria_per_experiment(tmp_path, monkeypatch):
    """Each curve is shaped with the manifest's precursor_enrichment (the SDRF
    RIA), and an explicit ria_override wins over it — the ¹⁸O head-to-head needs
    AC16=0.0583 vs iPSC=0.0897 honored, not the 0.06 config default."""
    import riana.core.fitting as fitting_mod

    coeffs = _coeffs()
    dfs = _make_timepoint_dfs(coeffs)
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)

    # Manifest rows carrying an experiment-level enrichment ≠ the 0.06 default.
    rows = []
    for ti, df in zip(_TIMES, dfs):
        stem = f"riaE_t{ti:.4f}"
        path = tmp_path / f"{stem}_riana.txt"
        df.to_csv(path, sep="\t", index=False)
        rows.append(ManifestRow("integrate", str(path), RunIdentity(
            experiment="syn", sample=stem, data_file=stem,
            labeling_time=float(ti), labeling_time_unit="au", condition="control",
            precursor_enrichment=0.0897)))
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)

    # Spy on the ria_max each curve's fit_run receives. fit_project re-imports
    # fit_run from riana.core.fitting per call, so patching the source attr works.
    seen: list[float] = []
    real = fitting_mod.fit_run

    def _spy(cfg, *a, **k):
        seen.append(cfg.ria_max)
        return real(cfg, *a, **k)

    monkeypatch.setattr(fitting_mod, "fit_run", _spy)

    fit_project(config, mf, coeffs, n_boot=0, random_state=0)
    assert seen == [pytest.approx(0.0897)]            # manifest enrichment, not 0.06

    seen.clear()
    fit_project(config, mf, coeffs, n_boot=0, random_state=0, ria_override=0.05)
    assert seen == [pytest.approx(0.05)]              # explicit override wins


def test_recombine_groups_by_condition(tmp_path):
    coeffs = _coeffs()
    rows = _integrate_rows_from_dfs(tmp_path, _make_timepoint_dfs(coeffs),
                                    condition="control")
    rows += _integrate_rows_from_dfs(tmp_path, _make_timepoint_dfs(coeffs),
                                     condition="knockout")
    curves = recombine_for_fit(rows)
    assert set(curves) == {("syn", "control"), ("syn", "knockout")}
    frame = curves[("syn", "control")]
    assert "labeling_time" in frame.columns
    # 8 timepoints x 5 peptides
    assert len(frame) == 8 * len(_TEST_PEPTIDES)


def test_recombine_merges_fractions(tmp_path):
    """Two fractions of the same (condition, biorep, time) sum at peptide level."""
    coeffs = _coeffs()
    df = _make_timepoint_dfs(coeffs)[3]  # one timepoint
    ti = float(_TIMES[3])
    rows = []
    for frac in (1, 2):
        stem = f"run_f{frac}"
        (tmp_path / f"{stem}_riana.txt").write_text(
            df.assign(file_idx=frac - 1).to_csv(sep="\t", index=False)
        )
        rows.append(ManifestRow("integrate", str(tmp_path / f"{stem}_riana.txt"),
            RunIdentity(experiment="syn", sample="s", data_file=stem,
                        labeling_time=ti, fraction=frac, condition="control")))
    frame = recombine_for_fit(rows)[("syn", "control")]
    # Fractions merged -> one row per peptide, iso channels summed (2x).
    assert len(frame) == len(_TEST_PEPTIDES)
    one = frame[frame["concat"] == "VAPEPTIDEK_2"].iloc[0]
    expected = 2 * df[df["concat"] == "VAPEPTIDEK_2"]["iso0"].iloc[0]
    assert one["iso0"] == pytest.approx(expected)


def test_bioreps_stay_independent_points(tmp_path):
    """Same condition+time, different biorep -> two points, not merged."""
    coeffs = _coeffs()
    df = _make_timepoint_dfs(coeffs)[3]
    ti = float(_TIMES[3])
    rows = []
    for rep in (1, 2):
        stem = f"run_r{rep}"
        (tmp_path / f"{stem}_riana.txt").write_text(df.to_csv(sep="\t", index=False))
        rows.append(ManifestRow("integrate", str(tmp_path / f"{stem}_riana.txt"),
            RunIdentity(experiment="syn", sample=f"s{rep}", data_file=stem,
                        labeling_time=ti, biological_replicate=rep,
                        condition="control")))
    frame = recombine_for_fit(rows)[("syn", "control")]
    # Two independent points per peptide (one per biorep).
    assert len(frame) == 2 * len(_TEST_PEPTIDES)


# --- parity: manifest path == legacy path ------------------------------------


def test_manifest_fit_matches_legacy_kdeg(tmp_path):
    coeffs = _coeffs()
    dfs = _make_timepoint_dfs(coeffs)
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)

    # Legacy positional path: time parsed from the sample string.
    legacy = fit_run(config, dfs, coeffs, n_boot=10, random_state=42)

    # Manifest path: time from the SDRF identity.
    rows = _integrate_rows_from_dfs(tmp_path, dfs)
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)
    manifested = fit_project(config, mf, coeffs, n_boot=10, random_state=42)

    # Point estimate (k_deg) and R^2 are order-independent -> must match exactly.
    legacy = legacy.sort_index()
    manifested = manifested.sort_index()
    assert list(legacy.index) == list(manifested.index)
    # rtol 1e-6 proves the science is identical; the residual is float
    # round-trip through the per-run _riana.txt CSV, not a method difference.
    np.testing.assert_allclose(
        legacy["k_deg"].to_numpy(dtype=float),
        manifested["k_deg"].to_numpy(dtype=float),
        rtol=1e-6, atol=1e-9, equal_nan=True,
    )
    np.testing.assert_allclose(
        legacy["R_squared"].to_numpy(dtype=float),
        manifested["R_squared"].to_numpy(dtype=float),
        rtol=1e-6, atol=1e-9, equal_nan=True,
    )
    # fit_project tags each curve with its group identity.
    assert set(manifested["condition"]) == {"control"}
    assert set(manifested["experiment"]) == {"syn"}

    # M5: the per-timepoint long table rides on .attrs, tagged by curve identity.
    long = manifested.attrs["fractions_long"]
    assert {"concat", "labeling_time", "fs", "fs_lower", "fs_upper",
            "biological_replicate", "experiment", "condition"} <= set(long.columns)
    assert set(long["condition"]) == {"control"}
    assert set(long["experiment"]) == {"syn"}
    # One row per (peptide, timepoint) of the fitted curve.
    assert len(long) == 8 * len(_TEST_PEPTIDES)


def test_fit_project_long_keeps_bioreps_separate(tmp_path):
    """Two bioreps at each timepoint -> distinct long rows per (peptide, biorep)."""
    coeffs = _coeffs()
    dfs = _make_timepoint_dfs(coeffs)
    rows = _integrate_rows_from_dfs(tmp_path, dfs, condition="control")
    # A second biorep of the same curve (distinct data_file stems).
    for ti, df in zip(_TIMES, dfs):
        stem = f"run2_t{ti:.4f}"
        path = tmp_path / f"{stem}_riana.txt"
        df.to_csv(path, sep="\t", index=False)
        rows.append(ManifestRow("integrate", str(path), RunIdentity(
            experiment="syn", sample=stem, data_file=stem,
            labeling_time=float(ti), labeling_time_unit="au",
            biological_replicate=2, condition="control")))
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)
    long = fit_project(config, mf, coeffs, n_boot=20).attrs["fractions_long"]
    assert set(long["biological_replicate"]) == {1, 2}
    # Each (peptide, biorep) contributes its own 8 timepoints.
    assert len(long) == 2 * 8 * len(_TEST_PEPTIDES)


def test_fit_project_two_conditions_concat(tmp_path):
    """Regression: a two-condition fit must not crash at the final concat.

    Each per-curve result frame carries its fractions_long table on ``.attrs``;
    ``pd.concat`` → ``__finalize__`` reconciles ``.attrs`` by equality-comparing
    values across frames, which raised "Can only compare identically-labeled"
    on the differently-shaped DataFrames once there was >1 curve (the Δk path).
    Single-condition fits never reconciled attrs across frames, so this slipped
    through until the two-condition lve_atr fixture.
    """
    coeffs = _coeffs()
    rows = []
    # Distinct file stems per condition (real LVE/ATR runs are distinct mzMLs);
    # the manifest keys on output_path, so a shared stem would collapse to one.
    for cond in ("control", "atrium"):
        for ti, df in zip(_TIMES, _make_timepoint_dfs(coeffs)):
            stem = f"{cond}_t{ti:.4f}"
            path = tmp_path / f"{stem}_riana.txt"
            df.to_csv(path, sep="\t", index=False)
            rows.append(ManifestRow("integrate", str(path), RunIdentity(
                experiment="syn", sample=stem, data_file=stem,
                labeling_time=float(ti), labeling_time_unit="au",
                condition=cond)))
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)
    config = FitConfig(model="simple", label="hw", q_value=0.05, depth=3,
                       ria_max=0.06)
    out = fit_project(config, mf, coeffs, n_boot=20)
    # Both curves survive the concat, tagged distinctly.
    assert set(out["condition"]) == {"control", "atrium"}
    long = out.attrs["fractions_long"]
    assert set(long["condition"]) == {"control", "atrium"}
    # Each condition contributes its own 8 timepoints × peptides.
    assert len(long) == 2 * 8 * len(_TEST_PEPTIDES)


# --- integrate_project end-to-end on the real BSA mzML -----------------------

# retention_time values (541.1 s @ scan 4408, 732.6 s @ scan 6838) are the BSA
# mzML's real MS1 RTs at those scans, so the intake scan↔RT guard reconciles.
_BSA_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t1.0.0
    MTD\tms_run[1]-location\tfile://20180216_BSA.mzML

    PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score[1]\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tspectra_ref\tpre\tpost\tstart\tend\topt_global_Posterior_Error_Probability_score\topt_global_q-value\topt_global_cv_MS:1002217_decoy_peptide\topt_global_cv_MS:1000889_peptidoform_sequence
    PSM\tRHPEYAVSVLLR\t0\tsp|P02769|ALBU_BOVIN\t1\tdb\tnull\t[, , dummy, 1]\t0.001\tnull\t541.1\t3\t470.6\t470.6\tms_run[1]:controllerType=0 controllerNumber=1 scan=4408\tK\tR\t1\t12\t0.01\t0.001\t0\tRHPEYAVSVLLR
    PSM\tHPYFYAPELLYYANK\t1\tsp|P02769|ALBU_BOVIN\t1\tdb\tnull\t[, , dummy, 1]\t0.001\tnull\t732.6\t3\t620.3\t620.3\tms_run[1]:controllerType=0 controllerNumber=1 scan=6838\tK\tR\t1\t15\t0.01\t0.001\t0\tHPYFYAPELLYYANK
    """)

_BSA_SDRF = textwrap.dedent("""\
    source name\tcharacteristics[biological replicate]\tcharacteristics[precursor enrichment]\tcharacteristics[labeling time]\tcomment[data file]\tcomment[fraction identifier]\tcomment[technical replicate]\tcomment[proteomics data acquisition method]\tfactor value[condition]
    bsa_t0\t1\t0.06\t0 day\t20180216_BSA.mzML\t1\t1\tNT=data-dependent acquisition;AC=MS:1003221\tcontrol
    """)


@pytest.mark.skipif(not BSA_MZML.exists(), reason="sample1 BSA mzML missing")
def test_integrate_project_end_to_end(tmp_path):
    mztab = tmp_path / "bsa.mzTab"
    mztab.write_text(_BSA_MZTAB)
    sdrf_path = tmp_path / "bsa.sdrf.tsv"
    sdrf_path.write_text(_BSA_SDRF)
    sdrf = read_sdrf(sdrf_path)

    config = IntegrationConfig(
        isotopomers=(0, 1, 2, 3, 4, 5), mass_tol_ppm=25,
        peak_rt="ms2", integration_half_width=1.0, extraction_half_width=1.0,
        out_dir=str(tmp_path),
    )
    rows = integrate_project(config, sdrf, SAMPLE1, mztab, tmp_path)

    # One <stem>_riana.txt per run.
    out_file = tmp_path / "20180216_BSA_riana.txt"
    assert out_file.exists()
    assert len(rows) == 1 and rows[0].identity.sample == "bsa_t0"

    # Identity is frozen into the provenance header.
    header = "\n".join(
        ln for ln in out_file.read_text().splitlines() if ln.startswith("#")
    )
    assert "# experiment bsa" in header           # SDRF stem
    assert "# experiment_type turnover" in header
    assert "# sample bsa_t0" in header
    assert "# labeling_time 0.0" in header
    assert "# condition control" in header

    # The manifest indexes the run.
    mf_rows = read_manifest(tmp_path / "riana_manifest.tsv", stage="integrate")
    assert len(mf_rows) == 1
    assert mf_rows[0].identity.data_file == "20180216_BSA"
    assert mf_rows[0].output_path == str(out_file)

    # The data body has the integrated isotopomers for the 2 PSMs.
    body = pd.read_table(out_file, comment="#")
    assert {"iso0", "iso5", "concat", "sample"} <= set(body.columns)
    assert (body["sample"] == "bsa_t0").all()


@pytest.mark.skipif(not BSA_MZML.exists(), reason="sample1 BSA mzML missing")
def test_integrate_project_resume_skips_done_runs(tmp_path):
    """`resume=True` keeps an already-integrated run (file not rewritten)."""
    import time

    mztab = tmp_path / "bsa.mzTab"
    mztab.write_text(_BSA_MZTAB)
    (tmp_path / "bsa.sdrf.tsv").write_text(_BSA_SDRF)
    sdrf = read_sdrf(tmp_path / "bsa.sdrf.tsv")
    config = IntegrationConfig(
        isotopomers=(0, 1, 2, 3, 4, 5), mass_tol_ppm=25, peak_rt="ms2",
        integration_half_width=1.0, extraction_half_width=1.0, out_dir=str(tmp_path),
    )
    integrate_project(config, sdrf, SAMPLE1, mztab, tmp_path)
    out_file = tmp_path / "20180216_BSA_riana.txt"
    mtime = out_file.stat().st_mtime_ns

    time.sleep(0.01)
    rows = integrate_project(config, sdrf, SAMPLE1, mztab, tmp_path, resume=True)
    assert len(rows) == 1                              # the kept run is returned
    assert out_file.stat().st_mtime_ns == mtime        # not rewritten


def test_resume_partition_logic(tmp_path):
    """The resume filter: keep done runs (file + matching config_hash), run rest."""
    import logging

    from riana.core.pipeline import RunTask, _resume_partition
    from riana.io.manifest import ManifestRow, append_manifest

    ident0 = RunIdentity(experiment="e", sample="s0", data_file="r0")
    t0 = RunTask(0, "r0", ident0, "a.mzML", [])
    t1 = RunTask(1, "r1", RunIdentity(experiment="e", sample="s1", data_file="r1"),
                 "b.mzML", [])
    tasks = [t0, t1]
    (tmp_path / "r0_riana.txt").write_text("done")     # r0 output exists
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, [ManifestRow(
        "integrate", str(tmp_path / "r0_riana.txt"), ident0, config_hash="H")])
    log = logging.getLogger("t")

    # resume + matching hash -> r0 kept, r1 runs.
    to_run, kept = _resume_partition(tasks, mf, tmp_path, "H", True, log)
    assert [t.stem for t in to_run] == ["r1"] and set(kept) == {0}
    # different hash (settings changed) -> both run.
    to_run2, kept2 = _resume_partition(tasks, mf, tmp_path, "OTHER", True, log)
    assert {t.stem for t in to_run2} == {"r0", "r1"} and kept2 == {}
    # resume=False -> everything runs.
    to_run3, _ = _resume_partition(tasks, mf, tmp_path, "H", False, log)
    assert len(to_run3) == 2
    # output file missing -> not kept even if the manifest row matches.
    (tmp_path / "r0_riana.txt").unlink()
    to_run4, kept4 = _resume_partition(tasks, mf, tmp_path, "H", True, log)
    assert {t.stem for t in to_run4} == {"r0", "r1"} and kept4 == {}


def test_record_stage_rows_and_fit_outputs_roundtrip(tmp_path):
    """fit/rollup register their outputs as stage rows; rollup finds them back."""
    from riana.core.pipeline import (
        aggregate_identity,
        fit_outputs_from_manifest,
        record_stage_rows,
    )
    from riana.io.writers import make_provenance

    coeffs = _coeffs()
    rows = _integrate_rows_from_dfs(tmp_path, _make_timepoint_dfs(coeffs),
                                    condition="control")
    mf = tmp_path / "riana_manifest.tsv"
    append_manifest(mf, rows)

    result = pd.DataFrame(
        {"k_deg": [0.5], "experiment": ["syn"], "condition": ["control"]})
    ident = aggregate_identity(result)
    assert ident.experiment == "syn" and ident.condition == "control"

    prov = make_provenance({"x": 1})
    pep = tmp_path / "riana_fit_peptides.txt"
    pep.write_text("x")
    frac = tmp_path / "riana_fit_fractions.txt"
    frac.write_text("x")
    record_stage_rows(mf, "fit", [pep, frac], result, prov)

    got_pep, got_frac = fit_outputs_from_manifest(mf)
    assert Path(got_pep).name == "riana_fit_peptides.txt"
    assert Path(got_frac).name == "riana_fit_fractions.txt"
    # The integrate rows are preserved alongside the new fit rows.
    assert len(read_manifest(mf, stage="integrate")) == len(_TIMES)
    assert len(read_manifest(mf, stage="fit")) == 2


def test_aggregate_identity_blanks_mixed_groups():
    from riana.core.pipeline import aggregate_identity
    mixed = pd.DataFrame({"experiment": ["a", "b"], "condition": ["x", "x"]})
    ident = aggregate_identity(mixed)
    assert ident.experiment == ""          # >1 experiment -> blank
    assert ident.condition == "x"          # single condition -> kept


def test_identity_to_extra_omits_none():
    cal = RunIdentity(experiment="c", sample="s", data_file="d",
                      mixing_proportion=0.5)
    extra = identity_to_extra(cal)
    assert extra["experiment_type"] == "calibration"
    assert extra["mixing_proportion"] == "0.5"
    assert "labeling_time" not in extra  # None omitted
    assert "precursor_enrichment" not in extra
