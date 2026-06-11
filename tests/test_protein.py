"""Tests for ``riana.core.protein`` — the Track C protein rollup.

Two layers: (1) unit tests of the summarize-time parsimony resolver
(``unique`` / ``isoform`` attribution over the peptide↔protein map); (2) the
rollup estimators on synthetic, deterministic frames whose per-timepoint θ
follow ``1 − exp(−k·t)`` exactly, so the median-of-k and the biorep-aware
inverse-variance weighted refit both recover the planted k.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from riana.core.protein import (
    _r2_admitted,
    _resolve_parsimony,
    _weighted_theta,
    rollup_proteins,
)
from riana.exceptions import DataError


def _make_frames(
    k_by_protein: dict[str, float],
    *,
    n_pep: int = 3,
    times=(0.5, 1.0, 2.0, 3.0, 4.0),
    bioreps=(1,),
    add_shared: bool = False,
):
    """Build (peptides, fractions) frames with planted per-protein kinetics.

    Protein ids are ``sp|ACC|NAME`` so the parsimony resolver maps them to the
    bare accession ``ACC`` (the value that lands in the output ``protein`` col).
    """
    pep_rows, frac_rows = [], []
    for pi, (prot, k) in enumerate(k_by_protein.items()):
        for j in range(n_pep):
            concat = f"PEP{pi}_{j}_2"
            pep_rows.append({"concat": concat, "protein id": prot, "k_deg": k})
            for br in bioreps:
                for t in times:
                    theta = 1.0 - np.exp(-k * t)
                    frac_rows.append({
                        "concat": concat, "protein id": prot,
                        "biological_replicate": br, "labeling_time": t,
                        "fs": theta, "fs_lower": theta - 0.01,
                        "fs_upper": theta + 0.01,
                    })
    if add_shared:
        concat = "SHARED_2"
        pep_rows.append({"concat": concat, "protein id": "sp|A|X,sp|B|Y",
                         "k_deg": 0.9})
        for br in bioreps:
            for t in times:
                theta = 1.0 - np.exp(-0.9 * t)
                frac_rows.append({
                    "concat": concat, "protein id": "sp|A|X,sp|B|Y",
                    "biological_replicate": br, "labeling_time": t,
                    "fs": theta, "fs_lower": theta - 0.01, "fs_upper": theta + 0.01,
                })
    return pd.DataFrame(pep_rows), pd.DataFrame(frac_rows)


def _peps(rows):
    return pd.DataFrame(rows, columns=["concat", "protein id"])


# --- parsimony resolver (unit) ----------------------------------------------- #
def test_unique_resolver_strips_to_accession_and_drops_shared():
    pep = _peps([
        ("u_2", "sp|P12345|G_HUMAN"),
        ("iso_2", "sp|P12345-2|G_HUMAN"),       # unique to an isoform -> kept as-is
        ("shared_2", "sp|A|X,sp|B|Y"),          # multi-accession -> dropped
    ])
    m = dict(zip(*_resolve_parsimony(pep, "unique")[["concat", "protein"]].values.T))
    assert m == {"u_2": "P12345", "iso_2": "P12345-2"}


def test_isoform_single_accession_kept_canonical_or_isoform():
    pep = _peps([
        ("canon_2", "sp|P12345|G_HUMAN"),
        ("iso_2", "sp|P12345-2|G_HUMAN"),
    ])
    m = dict(zip(*_resolve_parsimony(pep, "isoform")[["concat", "protein"]].values.T))
    assert m == {"canon_2": "P12345", "iso_2": "P12345-2"}


def test_isoform_shared_folds_into_canonical_when_no_isoform_unique():
    pep = _peps([
        ("shared_2", "sp|P12345|G,sp|P12345-2|G"),
        ("canon_only_2", "sp|P12345|G"),
    ])
    m = dict(zip(*_resolve_parsimony(pep, "isoform")[["concat", "protein"]].values.T))
    assert m["shared_2"] == "P12345"      # no isoform-unique evidence -> canonical
    assert m["canon_only_2"] == "P12345"


def test_isoform_shared_excluded_when_an_isoform_has_a_unique_peptide():
    pep = _peps([
        ("shared_2", "sp|P12345|G,sp|P12345-2|G"),
        ("iso_unique_2", "sp|P12345-2|G"),     # the isoform is physically real
    ])
    m = dict(zip(*_resolve_parsimony(pep, "isoform")[["concat", "protein"]].values.T))
    assert "shared_2" not in m             # ambiguous -> excluded
    assert m["iso_unique_2"] == "P12345-2"


def test_isoform_multi_gene_shared_is_rejected():
    pep = _peps([("para_2", "sp|P13541|MYH3,sp|Q02566|MYH6")])
    m = _resolve_parsimony(pep, "isoform")
    assert "para_2" not in set(m["concat"])


def test_isoform_all_isoforms_no_canonical_falls_back_to_base():
    # The R reference left collapsed_uniprot NA here; we fall back to the base.
    pep = _peps([("shared_2", "sp|P12345-2|G,sp|P12345-3|G")])
    m = dict(zip(*_resolve_parsimony(pep, "isoform")[["concat", "protein"]].values.T))
    assert m["shared_2"] == "P12345"


def test_isoform_handles_jcast_dash_J_suffix():
    # JCAST emits '-J1' / '-J2' isoforms; they must strip to the same base.
    pep = _peps([
        ("shared_2", "sp|P12345|G,sp|P12345-J1|G"),  # canonical + JCAST isoform
        ("canon_only_2", "sp|P12345|G"),
    ])
    m = dict(zip(*_resolve_parsimony(pep, "isoform")[["concat", "protein"]].values.T))
    assert m["shared_2"] == "P12345"      # one gene -> folds into canonical
    # And when the -J1 isoform has its own unique peptide, the shared one drops.
    pep2 = _peps([
        ("shared_2", "sp|P12345|G,sp|P12345-J1|G"),
        ("jiso_unique_2", "sp|P12345-J1|G"),
    ])
    m2 = dict(zip(*_resolve_parsimony(pep2, "isoform")[["concat", "protein"]].values.T))
    assert "shared_2" not in m2
    assert m2["jiso_unique_2"] == "P12345-J1"


# --- rollup estimators ------------------------------------------------------- #
def test_unique_parsimony_drops_shared_peptides():
    pep, frac = _make_frames({"sp|P0|X": 0.5}, add_shared=True)
    out = rollup_proteins(pep, frac, n_boot=20)
    assert set(out["protein"]) == {"P0"}


def test_peptide_median_k_column_recovers_planted_k():
    pep, frac = _make_frames({"sp|P0|X": 0.3, "sp|P1|Y": 0.8}, n_pep=4)
    out = rollup_proteins(pep, frac, n_boot=20).set_index("protein")
    assert out.loc["P0", "peptide_median_k"] == pytest.approx(0.3, abs=1e-9)
    assert out.loc["P1", "peptide_median_k"] == pytest.approx(0.8, abs=1e-9)
    assert out.loc["P0", "n_peptides"] == 4


def test_weighted_refit_recovers_planted_k():
    pep, frac = _make_frames({"sp|P0|X": 0.5})
    out = rollup_proteins(pep, frac, n_boot=50).set_index("protein")
    assert (out["method"] == "weighted").all()
    assert out.loc["P0", "k_deg"] == pytest.approx(0.5, abs=0.02)
    assert out.loc["P0", "R_squared"] > 0.999
    assert out.loc["P0", "n_points"] == 5         # collapsed (biorep, t) cells


def test_pooled_method_uses_all_points():
    pep, frac = _make_frames({"sp|P0|X": 0.5}, n_pep=3)  # 3 peptides x 5 t
    out = rollup_proteins(pep, frac, method="pooled", n_boot=50).set_index("protein")
    assert (out["method"] == "pooled").all()
    assert out.loc["P0", "k_deg"] == pytest.approx(0.5, abs=0.02)
    assert out.loc["P0", "n_points"] == 15        # all peptide×timepoint points


def test_bioreps_are_independent_refit_points():
    pep, frac = _make_frames({"sp|P0|X": 0.5}, bioreps=(1, 2))
    out = rollup_proteins(pep, frac, n_boot=20).set_index("protein")
    assert out.loc["P0", "n_points"] == 10
    assert out.loc["P0", "k_deg"] == pytest.approx(0.5, abs=0.02)


def test_weighted_theta_favours_tighter_ci():
    fs = np.array([0.5, 0.6])
    lo = np.array([0.49, 0.40])
    hi = np.array([0.51, 0.80])
    assert _weighted_theta(fs, lo, hi) == pytest.approx(0.5, abs=0.01)


def test_weighted_theta_falls_back_to_unweighted_when_no_ci():
    fs = np.array([0.4, 0.6])
    nan = np.array([np.nan, np.nan])
    assert _weighted_theta(fs, nan, nan) == pytest.approx(0.5)


def test_min_peptides_filters_small_proteins():
    pep, frac = _make_frames({"sp|P0|X": 0.5}, n_pep=1)
    out = rollup_proteins(pep, frac, min_peptides=2, n_boot=20)
    assert out.empty or "P0" not in set(out["protein"])


def test_condition_groups_stay_separate():
    pep, frac = _make_frames({"sp|P0|X": 0.5})
    pep2, frac2 = _make_frames({"sp|P0|X": 0.5})
    for f in (pep, frac):
        f["condition"] = "control"
    for f in (pep2, frac2):
        f["condition"] = "knockout"
    out = rollup_proteins(
        pd.concat([pep, pep2], ignore_index=True),
        pd.concat([frac, frac2], ignore_index=True),
        n_boot=20,
    )
    assert set(out["condition"]) == {"control", "knockout"}
    assert len(out) == 2  # same protein, two conditions


def test_isoform_rollup_runs_end_to_end():
    # Canonical + an isoform with no unique peptide -> both fold to canonical,
    # so the protein has 2x the peptides rolling up.
    pep, frac = _make_frames({"sp|P12345|G_HUMAN": 0.5})
    # Relabel half the peptides as shared canonical+isoform (no isoform-unique).
    shared = pep["concat"].isin(["PEP0_0_2"])
    for df in (pep, frac):
        df.loc[df["concat"] == "PEP0_0_2", "protein id"] = \
            "sp|P12345|G_HUMAN,sp|P12345-2|G_HUMAN"
    out = rollup_proteins(pep, frac, parsimony="isoform", n_boot=20)
    assert set(out["protein"]) == {"P12345"}
    assert out.iloc[0]["n_peptides"] == 3  # the shared peptide still counts


def test_r2_admit_gate_excludes_noisy_but_keeps_slow_turnover():
    pep = pd.DataFrame({
        "concat": ["good_2", "noisy_2", "slow_2"],
        "R_squared": [0.95, 0.40, 0.30],   # noisy + slow both below 0.8
        "k_deg": [0.50, 0.50, 0.010],      # slow_2 turns over very slowly
        "sd": [0.02, 0.30, 0.020],         # slow_2 is well-determined (low SE)
    })
    admitted = _r2_admitted(pep, min_r2=0.8, alt_k=0.025, alt_se=0.05, alt_r2=0.0)
    assert admitted == {"good_2", "slow_2"}   # noisy excluded; slow admitted via alt


def test_r2_gate_off_by_default_keeps_everything():
    # A deliberately low-R² peptide stays in when min_r2 is None (default).
    pep, frac = _make_frames({"sp|P0|X": 0.5}, n_pep=3)
    pep["R_squared"] = [0.2, 0.95, 0.95]
    pep["sd"] = 0.3
    out_off = rollup_proteins(pep, frac, n_boot=20).set_index("protein")
    assert out_off.loc["P0", "n_peptides"] == 3
    out_on = rollup_proteins(pep, frac, min_r2=0.8, n_boot=20).set_index("protein")
    assert out_on.loc["P0", "n_peptides"] == 2  # the R²=0.2 peptide gated out


def test_rollup_threads_give_identical_result():
    """Per-protein RNG streams make the rollup independent of thread count."""
    pep, frac = _make_frames(
        {"sp|P0|X": 0.3, "sp|P1|Y": 0.6, "sp|P2|Z": 0.9}, n_pep=4)
    a = rollup_proteins(pep, frac, n_boot=100, threads=1)
    b = rollup_proteins(pep, frac, n_boot=100, threads=4)
    pd.testing.assert_frame_equal(a, b)


def test_unknown_model_parsimony_method_raise():
    pep, frac = _make_frames({"sp|P0|X": 0.5})
    with pytest.raises(DataError, match="model"):
        rollup_proteins(pep, frac, model="nope")
    with pytest.raises(DataError, match="parsimony"):
        rollup_proteins(pep, frac, parsimony="razor")
    with pytest.raises(DataError, match="method"):
        rollup_proteins(pep, frac, method="bogus")
