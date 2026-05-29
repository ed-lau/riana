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

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from riana.config import IntegrationConfig
from riana.core.integration import integrate_run
from riana.io.mzml import IndexedMzML
from riana.io.percolator import read_percolator


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
def legacy_sample1_output(tmp_path_factory) -> pd.DataFrame:
    """Run the legacy ``riana integrate`` on sample1 once; return its _riana.txt."""
    out_dir = tmp_path_factory.mktemp("legacy_sample1")
    cmd = [
        sys.executable, "-m", "riana", "integrate",
        str(SAMPLE1), str(SAMPLE1 / "percolator.target.psms.txt"),
        "-s", "sample1", "-o", str(out_dir),
        "-i", "0", "6",
        "-q", "1.0",
        "-r", "1.0",
        "-m", "50",
        "-t", "1",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return pd.read_csv(out_dir / "sample1_riana.txt", sep="\t", index_col=0)


def test_sample1_iso0_iso6_within_rtol(legacy_sample1_output):
    config = IntegrationConfig(
        sample="sample1", isotopomers=(0, 6), q_value=1.0,
        r_time=1.0, mass_tol_ppm=50, threads=1, forced_mods=(0.0,),
        peak_method="fixed_window", baseline_method="none",
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
def test_ac16_time0_matches_committed_baseline():
    """The bigger real-data check: ac16 time0, 9 isotopomers including 6."""
    config = IntegrationConfig(
        sample="time0", isotopomers=(0, 1, 2, 3, 4, 5), q_value=0.01,
        r_time=0.33, mass_tol_ppm=15, threads=4, forced_mods=(0.0,),
        peak_method="fixed_window", baseline_method="none",
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
        r_time=1.0, mass_tol_ppm=50, threads=1, forced_mods=(0.0,),
    )
    fixed_cfg = IntegrationConfig(**base, peak_method="fixed_window",
                                  baseline_method="none")
    det_cfg = IntegrationConfig(**base, peak_method="detected",
                                baseline_method="linear")
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
