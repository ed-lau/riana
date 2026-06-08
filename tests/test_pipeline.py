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
                       ria_max=0.06, threads=1)

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


# --- integrate_project end-to-end on the real BSA mzML -----------------------

_BSA_MZTAB = textwrap.dedent("""\
    MTD\tmzTab-version\t1.0.0
    MTD\tms_run[1]-location\tfile://20180216_BSA.mzML

    PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score[1]\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tspectra_ref\tpre\tpost\tstart\tend\topt_global_Posterior_Error_Probability_score\topt_global_q-value\topt_global_cv_MS:1002217_decoy_peptide\topt_global_cv_MS:1000889_peptidoform_sequence
    PSM\tRHPEYAVSVLLR\t0\tsp|P02769|ALBU_BOVIN\t1\tdb\tnull\t[, , dummy, 1]\t0.001\tnull\t300.0\t3\t470.6\t470.6\tms_run[1]:controllerType=0 controllerNumber=1 scan=4408\tK\tR\t1\t12\t0.01\t0.001\t0\tRHPEYAVSVLLR
    PSM\tHPYFYAPELLYYANK\t1\tsp|P02769|ALBU_BOVIN\t1\tdb\tnull\t[, , dummy, 1]\t0.001\tnull\t400.0\t3\t620.3\t620.3\tms_run[1]:controllerType=0 controllerNumber=1 scan=6838\tK\tR\t1\t15\t0.01\t0.001\t0\tHPYFYAPELLYYANK
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


def test_identity_to_extra_omits_none():
    cal = RunIdentity(experiment="c", sample="s", data_file="d",
                      mixing_proportion=0.5)
    extra = identity_to_extra(cal)
    assert extra["experiment_type"] == "calibration"
    assert extra["mixing_proportion"] == "0.5"
    assert "labeling_time" not in extra  # None omitted
    assert "precursor_enrichment" not in extra
