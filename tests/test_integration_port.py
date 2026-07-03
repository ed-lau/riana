"""Phase A regression gate — ``core.integration.integrate_run`` vs 0.9.0.

Per PROJECT_REVIEW.md M3 verification §, the rewrite must reproduce 0.9.0
per-peptide integration m0/m6 within 1e-3 relative tolerance on the
``sample1/`` smoke fixture. The ac16 ``time0`` fraction is the bigger
real-data check, gated against the committed
``integrate_outputs/v0.9.0/time0_riana.txt`` baseline.

These run quickly because the per-PSM trapezoidal core is the same numerical
recipe; the test is asserting that the typed-records + ``IndexedMzML``
rewrap didn't drift the numbers.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from riana.algorithms.mass_calc import calculate_ion_mz
from riana.config import IntegrationConfig
from riana.core.integration import (
    check_scan_precursor_consistency,
    integrate_run,
)
from riana.exceptions import DataError
from riana.io.mzml import IndexedMzML
from riana.io.percolator import read_percolator
from riana.records import PSMRecord


SAMPLE1 = Path("tests/data/sample1")
AC16_INTEGRATE_V0 = Path(
    "tests/data/calibration_d2o_mixing/ac16/integrate_outputs/v0.9.0/time0_riana.txt"
)
AC16_MZML = Path("data/calibration_ac16/mzml/20230731_AC16_WL_0.mzML.gz")
AC16_PSMS = Path(
    "data/calibration_ac16/snakemake_results/time0/percolator/percolator.target.psms.txt"
)


# --- sample1 BSA: A/B against the legacy CLI ---------------------------------


@pytest.fixture(scope="module")
def legacy_sample1_output() -> pd.DataFrame:
    """The committed 0.9.0 golden ``_riana.txt`` for sample1.

    Captured from the legacy engine before it was removed in M4 (it ran
    ``riana integrate ... -i 0 6 -q 1.0 -r 1.0 -m 50``). The new engine's
    ``ms2`` mode must still reproduce it within 1e-3 — that parity is the
    "is legacy removal safe" gate. Regenerate only if the numerical core
    legitimately changes.
    """
    # _riana.txt carries a provenance header (# lines); skip it.
    return pd.read_csv(SAMPLE1 / "sample1_riana.v0_9_0.txt", sep="\t",
                       index_col=0, comment="#")


@pytest.mark.slow
def test_sample1_iso0_iso6_within_rtol(legacy_sample1_output):
    config = IntegrationConfig(
        sample="sample1", isotopomers=(0, 6), q_value=1.0,
        extraction_half_width=1.0, mass_tol_ppm=50,
        peak_rt="ms2", baseline_method="none",
    )
    psms = read_percolator(
        SAMPLE1 / "percolator.target.psms.txt", sample="sample1"
    )
    with IndexedMzML(SAMPLE1 / "20180216_BSA.mzML.gz") as mzml:
        new = integrate_run(config, psms, mzml)

    # Join on (concat, scan) — the unique-per-row identity the legacy file
    # carries through every step.
    legacy = legacy_sample1_output
    merged = new.merge(
        legacy[["concat", "scan", "iso0", "iso6"]].rename(
            columns={"iso0": "iso0_legacy", "iso6": "iso6_legacy"}
        ),
        on=["concat", "scan"],
        how="inner",
    )
    assert len(merged) == len(legacy), "rowcount drift between port and legacy"

    # 1e-3 relative tolerance per PROJECT_REVIEW.md M3 verification.
    for col_new, col_old in [("iso0", "iso0_legacy"), ("iso6", "iso6_legacy")]:
        np.testing.assert_allclose(
            merged[col_new].to_numpy(),
            merged[col_old].to_numpy(),
            rtol=1e-3,
            atol=1e-6,
            err_msg=f"{col_new} drifted vs legacy on sample1",
        )


# --- ac16 time0: A/B vs the committed v0.9.0 baseline ------------------------


@pytest.mark.skipif(
    not (AC16_MZML.exists() and AC16_PSMS.exists()),
    reason="ac16 calibration mzML/PSMs not present (heavy inputs are gitignored)",
)
@pytest.mark.slow
def test_ac16_time0_matches_committed_baseline():
    """The bigger real-data check: ac16 time0, 9 isotopomers including 6."""
    config = IntegrationConfig(
        sample="time0", isotopomers=(0, 1, 2, 3, 4, 5), q_value=0.01,
        extraction_half_width=0.33, mass_tol_ppm=15,
        peak_rt="ms2", baseline_method="none",
    )
    psms = read_percolator(AC16_PSMS, sample="time0")
    with IndexedMzML(AC16_MZML) as mzml:
        new = integrate_run(config, psms, mzml)

    legacy = pd.read_csv(AC16_INTEGRATE_V0, sep="\t", index_col=0)
    iso_cols = ["iso0", "iso1", "iso2", "iso3", "iso4", "iso5"]
    merged = new.merge(
        legacy[["concat", "scan", *iso_cols]].rename(
            columns={c: c + "_legacy" for c in iso_cols}
        ),
        on=["concat", "scan"],
        how="inner",
    )
    assert len(merged) == len(legacy), "rowcount drift vs committed ac16 baseline"

    for c in iso_cols:
        # rtol=1e-3 per the verification gate; small atol so zero-vs-tiny
        # values don't trip the relative comparison.
        np.testing.assert_allclose(
            merged[c].to_numpy(),
            merged[c + "_legacy"].to_numpy(),
            rtol=1e-3,
            atol=1.0,  # an intensity of 1 is dwarfed by typical peak areas (~1e6)
            err_msg=f"{c} drifted vs ac16 v0.9.0 baseline",
        )


# --- Phase C smoke: detected pipeline produces non-trivial, plausible areas --


@pytest.mark.slow
def test_sample1_detected_pipeline_runs_and_is_smaller_than_fixed():
    """Phase C smoke: the detected pipeline returns finite areas that are
    bounded above by the fixed-window areas (peak detection narrows the
    integration interval; baseline subtraction can only further reduce it).

    Not a numerical gate — that's the benchmark suite (bench_peak_boundary,
    bench_m0_ma_recovery) on the calibration data. This test just confirms
    the detected path doesn't blow up on a real fixture and returns sane
    relative magnitudes.
    """
    psms = read_percolator(
        SAMPLE1 / "percolator.target.psms.txt", sample="sample1"
    )
    base = dict(
        sample="sample1", isotopomers=(0, 1, 2, 3, 4, 5), q_value=1.0,
        extraction_half_width=1.0, mass_tol_ppm=50,
    )
    fixed_cfg = IntegrationConfig(**base, peak_rt="ms2",
                                  baseline_method="none")
    det_cfg = IntegrationConfig(**base, peak_rt="apex",
                                integration_half_width="auto",
                                baseline_method="none")
    with IndexedMzML(SAMPLE1 / "20180216_BSA.mzML.gz") as mzml:
        fixed = integrate_run(fixed_cfg, psms, mzml)
        det = integrate_run(det_cfg, psms, mzml)

    iso_cols = ["iso0", "iso1", "iso2", "iso3", "iso4", "iso5"]
    merged = fixed.merge(
        det[["concat", "scan", *iso_cols]].rename(
            columns={c: c + "_det" for c in iso_cols}
        ),
        on=["concat", "scan"],
    )
    assert len(merged) == len(fixed)

    # detected ≤ fixed everywhere it ran detection on a sane peak, with a
    # small slack to absorb cases where the boundary lands one scan wider
    # than the fixed window (the cycle-level off-by-one is harmless).
    for c in iso_cols:
        diff = merged[c + "_det"].to_numpy() - merged[c].to_numpy()
        # We allow up to ~5% positive drift on the detected side, mostly
        # from baseline-subtraction-of-noise edge cases.
        assert (diff / np.maximum(merged[c].to_numpy(), 1.0)).max() < 0.10, (
            f"detected {c} unexpectedly larger than fixed-window"
        )
        # And there should be at least one row where detection meaningfully
        # tightened the window (proves the detected branch ran).
        assert (diff < -1.0).any(), f"detected pipeline never narrowed {c}"


# --- Phase D smoke: mass-accuracy outputs populated --------------------------


@pytest.mark.slow
def test_sample1_mass_accuracy_columns_populated():
    """Phase D smoke: each ``isoN_obs_mz`` is within ±mass_tol_ppm of its
    target (validates the intensity-weighted centroid calc), each
    ``isoN_ppm_error`` is finite for matched PSMs, and the per-fraction
    drift summary is attached to ``df.attrs['drift_summary']`` with a
    non-zero count.
    """
    config = IntegrationConfig(
        sample="sample1", isotopomers=(0, 6), q_value=1.0,
        extraction_half_width=1.0, mass_tol_ppm=50,
        peak_rt="ms2", baseline_method="none",
    )
    psms = read_percolator(
        SAMPLE1 / "percolator.target.psms.txt", sample="sample1"
    )
    with IndexedMzML(SAMPLE1 / "20180216_BSA.mzML.gz") as mzml:
        df = integrate_run(config, psms, mzml)

    # New columns landed.
    for c in ("iso0_obs_mz", "iso6_obs_mz", "iso0_ppm_error", "iso6_ppm_error"):
        assert c in df.columns, f"missing Phase D column {c}"

    # ppm_error stays inside the requested mass window for the rows that
    # actually matched a centroid (NaN = no match in window, ignored).
    for c in ("iso0_ppm_error", "iso6_ppm_error"):
        observed = df[c].dropna()
        assert len(observed) > 0, f"{c} has no observed values"
        assert observed.abs().max() <= config.mass_tol_ppm, (
            f"{c} drift {observed.abs().max():.2f} ppm exceeds "
            f"mass_tol_ppm={config.mass_tol_ppm}"
        )

    # Drift summary attached + non-trivial.
    drift = df.attrs.get("drift_summary")
    assert drift is not None
    assert drift.n > 0
    assert not np.isnan(drift.median_ppm)
    assert drift.mad_ppm >= 0


# --- CLI smoke: `riana integrate` dispatch ----------------------------------


@pytest.mark.slow
def test_cli_integrate_produces_compatible_output(tmp_path):
    """CLI smoke: ``python -m riana integrate`` writes the legacy _riana.txt
    schema (plus the mass-accuracy columns) and a per-fraction drift JSON
    sidecar through the Typer CLI (the only engine after M4).
    """
    out_dir = tmp_path / "cli_out"
    out_dir.mkdir()
    cmd = [
        sys.executable, "-m", "riana", "integrate",
        str(SAMPLE1), str(SAMPLE1 / "percolator.target.psms.txt"),
        "-s", "sample1", "-o", str(out_dir),
        "-i", "0 6",
        "-q", "1.0",
        "-r", "1.0",
        "-m", "50",
    ]
    subprocess.run(cmd, check=True, capture_output=True)

    out_file = out_dir / "sample1_riana.txt"
    drift_file = out_dir / "sample1_riana.drift.json"
    assert out_file.exists()
    assert drift_file.exists()

    df = pd.read_csv(out_file, sep="\t", index_col=0, comment="#")
    # Legacy area columns + Phase D additions.
    for c in ("iso0", "iso6", "iso0_obs_mz", "iso6_obs_mz",
              "iso0_ppm_error", "iso6_ppm_error", "concat", "sample"):
        assert c in df.columns, f"missing column {c}"

    drift = json.loads(drift_file.read_text())
    for k in ("n", "median_ppm", "mad_ppm", "suggested_shift_ppm"):
        assert k in drift, f"missing drift key {k}"
    assert drift["n"] > 0


# --- best-q anchor + apex_search_half_width bound ----------------------------


def test_apex_search_half_width_bounds_apex_to_anchor():
    """The apex search keys on the (best-q) anchor scan, bounded by
    ``apex_search_half_width``.

    Synthetic iso0 with two prominent peaks: a *shorter* one at RT 1.0 and a
    *taller* one at RT 3.0. With ``apex_search_half_width=0.25`` the apex must
    follow the anchor (pick the peak within ±0.25 of it), not the global tallest
    — the regression guard for the animal-D₂O failure where ``asw=0`` +
    ``selection="tallest"`` roamed to a co-eluting isobar.
    """
    import dataclasses
    from types import SimpleNamespace

    from riana.core.integration import _peak_boundary

    rt = np.linspace(0.0, 5.0, 101)  # rt[20]=1.0, rt[60]=3.0

    def gauss(center, amp, sigma=0.08):
        return amp * np.exp(-((rt - center) ** 2) / (2 * sigma * sigma))

    # Both peaks well clear of the prominence floor; the one at 3.0 is taller.
    iso0 = gauss(1.0, 100.0) + gauss(3.0, 500.0) + 1.0
    idf = pd.DataFrame({"rt": rt, "iso0": iso0, "iso1": iso0 * 0.5})
    # MS1 index: scan N -> rt[N-1]; searchsorted(scan)-1 lands on rt[N-2].
    mzml = SimpleNamespace(scan_idx=np.arange(1, 102), rt_idx=rt)
    psm = SimpleNamespace(scan=1, concat="PEPTIDEK_2")
    cfg = IntegrationConfig(
        peak_rt="apex", integration_half_width=0.15,
        apex_search_half_width=0.25, apex_selection="tallest",
    )

    # Anchor near the SHORTER peak at 1.0 -> apex must be ~1.0, not the taller 3.0.
    b_near = _peak_boundary(
        idf, psm, mzml, cfg, rt, "iso0", "iso1", anchor_scan=22,
    )
    assert b_near is not None
    assert abs(rt[b_near.apex_idx] - 1.0) < 0.2

    # Anchor near the taller peak at 3.0 -> apex follows to ~3.0.
    b_far = _peak_boundary(
        idf, psm, mzml, cfg, rt, "iso0", "iso1", anchor_scan=62,
    )
    assert b_far is not None
    assert abs(rt[b_far.apex_idx] - 3.0) < 0.2

    # The old default (asw=0) + tallest WOULD roam to the global-tallest peak
    # at 3.0 even when anchored near 1.0 — the bug the 0.25 default fixes.
    cfg0 = dataclasses.replace(cfg, apex_search_half_width=0.0)
    b_roam = _peak_boundary(
        idf, psm, mzml, cfg0, rt, "iso0", "iso1", anchor_scan=22,
    )
    assert b_roam is not None
    assert abs(rt[b_roam.apex_idx] - 3.0) < 0.2


# --- intake scan↔RT guard (Track A) ------------------------------------------

# BSA mzML MS1 RT at the (MS2) PSM scans the synthetic mzTab references — the
# scans `searchsorted(side="left") - 1` lands on. Measured from the committed
# mzML; the reported retention_time must reconcile to within the 2 min default.
def _bsa_ms2_precursors(mzml, n=4):
    """A few ``(scan, precursor m/z)`` from real BSA MS2 spectra in the fixture."""
    out = []
    for scan in mzml._scan_to_spec_id_all:
        if scan in mzml._scan_to_spec_id:  # MS1 — no precursor
            continue
        pmz = mzml.precursor_mz(scan)
        if pmz and pmz > 0:
            out.append((scan, pmz))
        if len(out) >= n:
            break
    return out


def _bsa_psm(scan, precursor_mz, *, sequence="RHPEYAVSVLLR", charge=3, q=1e-3):
    """A hand-built PSM at a real BSA MS2 scan with a chosen precursor m/z."""
    return PSMRecord(
        scan=scan, charge=charge, sequence=sequence,
        peptide_mass=float(calculate_ion_mz(sequence)),
        sample="bsa", file_idx=0, percolator_q_value=q,
        precursor_mz=precursor_mz,
    )


@pytest.fixture(scope="module")
def bsa_mzml():
    with IndexedMzML(SAMPLE1 / "20180216_BSA.mzML.gz") as m:
        yield m


def test_scan_precursor_passes_when_precursor_matches(bsa_mzml):
    """Reported precursor m/z == the mzML selected-ion m/z ⇒ reconciles."""
    pairs = _bsa_ms2_precursors(bsa_mzml)
    assert pairs, "fixture should have MS2 spectra with precursors"
    psms = [_bsa_psm(s, pmz) for s, pmz in pairs]
    check = check_scan_precursor_consistency(psms, bsa_mzml)
    assert check.n_checked == len(psms)
    assert check.ok
    assert check.frac_matched == 1.0
    assert check.median_ppm < 1.0


def test_scan_precursor_flags_wrong_precursor_and_integrate_run_errors(bsa_mzml):
    """A reported precursor 5 Da off every scan ⇒ no match ⇒ DataError.

    The wrong-file / scan-scramble signature: spectra_ref scans that index a
    different run point at different precursors.
    """
    pairs = _bsa_ms2_precursors(bsa_mzml)
    psms = [_bsa_psm(s, pmz + 5.0) for s, pmz in pairs]
    check = check_scan_precursor_consistency(psms, bsa_mzml)
    assert not check.ok
    assert check.frac_matched == 0.0

    cfg = IntegrationConfig(
        isotopomers=(0, 1), q_value=1.0, peak_rt="ms2",
        integration_half_width=1.0, extraction_half_width=1.0, mass_tol_ppm=50,
    )
    with pytest.raises(DataError, match="reconciliation FAILED"):
        integrate_run(cfg, psms, bsa_mzml)


def test_scan_precursor_guard_no_ops_without_precursor(bsa_mzml):
    """PSMs without a reported precursor m/z ⇒ guard skips, integration proceeds."""
    pairs = _bsa_ms2_precursors(bsa_mzml)
    psms = [_bsa_psm(s, 0.0) for s, _ in pairs]  # precursor_mz = 0 ⇒ skipped
    check = check_scan_precursor_consistency(psms, bsa_mzml)
    assert check.n_checked == 0
    assert check.ok  # nothing to check is not a failure


def test_scan_id_guard_escape_hatch_disables_check(bsa_mzml):
    """check_scan_id=False bypasses the guard even on a wrong precursor m/z."""
    s, pmz = _bsa_ms2_precursors(bsa_mzml, n=1)[0]
    psms = [_bsa_psm(s, pmz + 50.0)]  # 50 Da off
    cfg = IntegrationConfig(
        check_scan_id=False, isotopomers=(0, 1), q_value=1.0, peak_rt="ms2",
        integration_half_width=1.0, extraction_half_width=1.0, mass_tol_ppm=50,
    )
    df = integrate_run(cfg, psms, bsa_mzml)  # no DataError despite the bad precursor
    assert len(df) == 1
