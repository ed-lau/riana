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
    _k_cv_admitted,
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


# --- M7 Stage B proteoform keys ---------------------------------------------- #
def test_proteoform_key_appends_biological_mod_site():
    """A biological-mod ('mod sites') suffix is appended to the resolved accession;
    unmodified and constitutive-mod (N-term Ac → empty) peptidoforms stay bare."""
    pep = pd.DataFrame([
        {"concat": "PEPS[UNIMOD:21]K_2", "protein id": "sp|A2ASS6|TITIN_MOUSE",
         "mod sites": "pS34476"},
        {"concat": "PEPSK_2", "protein id": "sp|A2ASS6|TITIN_MOUSE", "mod sites": ""},
        {"concat": "[UNIMOD:1]PEPK_2", "protein id": "sp|A2ASS6|TITIN_MOUSE",
         "mod sites": ""},
    ])
    m = _resolve_parsimony(pep, "unique").set_index("concat")["protein"].to_dict()
    assert m["PEPS[UNIMOD:21]K_2"] == "A2ASS6_pS34476"
    assert m["PEPSK_2"] == "A2ASS6"
    assert m["[UNIMOD:1]PEPK_2"] == "A2ASS6"


def test_proteoform_key_is_bare_when_mod_sites_column_absent():
    """Back-compat: an input without a 'mod sites' column rolls up to the bare
    accession exactly as before Stage B."""
    pep = pd.DataFrame([{"concat": "PEPK_2", "protein id": "sp|P1|X"}])
    m = _resolve_parsimony(pep, "unique").set_index("concat")["protein"].to_dict()
    assert m["PEPK_2"] == "P1"


def test_rollup_separates_phospho_proteoform_from_bare_protein():
    """Two peptidoforms of one accession — a phospho form and the bare form —
    roll up as distinct units (P1_pS100 vs P1), not collapsed together."""
    k, times = 0.3, (0.5, 1.0, 2.0, 3.0, 4.0)
    pep_rows, frac_rows = [], []
    for concat, sites in [("BAREPEP_2", ""), ("PHOSPEP_2", "pS100")]:
        pep_rows.append({"concat": concat, "protein id": "sp|P1|X",
                         "mod sites": sites, "k_deg": k})
        for t in times:
            theta = 1.0 - np.exp(-k * t)
            frac_rows.append({"concat": concat, "protein id": "sp|P1|X",
                              "mod sites": sites, "biological_replicate": 1,
                              "labeling_time": t, "fs": theta,
                              "fs_lower": theta - 0.01, "fs_upper": theta + 0.01})
    out = rollup_proteins(pd.DataFrame(pep_rows), pd.DataFrame(frac_rows),
                          n_boot=20, min_peptides=1)
    proteins = set(out["protein"])
    assert "P1" in proteins         # the unmodified form → bare protein
    assert "P1_pS100" in proteins   # the phosphopeptidoform → its own unit


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


def test_k_cv_admitted_single_timepoint_gate():
    """The single-timepoint k_cv gate admits k_cv < threshold and rejects a wide CI
    or a NaN CI (single-point / non-converged)."""
    pep = pd.DataFrame([
        {"concat": "A", "k_deg": 0.02, "ci_lo": 0.019, "ci_hi": 0.021},  # k_cv=0.05 keep
        {"concat": "B", "k_deg": 0.02, "ci_lo": 0.005, "ci_hi": 0.05},   # k_cv=1.1 reject
        {"concat": "C", "k_deg": 0.02, "ci_lo": np.nan, "ci_hi": np.nan},  # NaN reject
    ])
    assert _k_cv_admitted(pep, 0.2) == {"A"}


def test_rollup_single_timepoint_auto_detects_and_curates_on_replicates_and_kcv():
    """One distinct labeling time auto-detects single-timepoint: R² is bypassed and
    curation = min_fit_points (auto 2) + the k_cv gate. A single-replicate peptide
    (n_points=1, k_cv NaN) drops with no flags; --min-fit-points 1 --k-cv 0 keeps it."""
    T = 24.0
    pep = pd.DataFrame([
        {"concat": "A_2", "protein id": "sp|P1|X", "k_deg": 0.02, "R_squared": 0.0,
         "ci_lo": 0.019, "ci_hi": 0.021, "n_points": 3},                   # keep
        {"concat": "B_2", "protein id": "sp|P1|X", "k_deg": 0.02, "R_squared": np.nan,
         "ci_lo": np.nan, "ci_hi": np.nan, "n_points": 1},                 # single -> drop
    ])
    rows = []
    for c, brs in (("A_2", [1, 2, 3]), ("B_2", [1])):
        for br in brs:
            rows.append({"concat": c, "protein id": "sp|P1|X",
                         "biological_replicate": br, "labeling_time": T,
                         "fs": 0.38, "fs_lower": 0.36, "fs_upper": 0.40})
    frac = pd.DataFrame(rows)
    auto = rollup_proteins(pep, frac, model="simple", min_points=1,
                           min_peptides=1, n_boot=20).set_index("protein")
    assert int(auto.loc["P1", "n_peptides"]) == 1     # only the replicated A survives
    off = rollup_proteins(pep, frac, model="simple", min_points=1, min_peptides=1,
                          min_fit_points=1, k_cv_max=0.0, n_boot=20).set_index("protein")
    assert int(off.loc["P1", "n_peptides"]) == 2      # gates disabled -> A and B


def test_rollup_min_spep_gate_drops_low_site_peptides():
    """--min-spep at rollup drops peptides below the floor before the refit
    (defense-in-depth for explicit-file inputs that weren't gated at fit)."""
    pep, frac = _make_frames({"sp|P0|X": 0.5}, n_pep=3)
    pep = pep.assign(spep=[10.0, 10.0, 4.0])      # PEP0_2_2 is under the floor
    gated = rollup_proteins(pep, frac, min_spep=6, n_boot=20).set_index("protein")
    open_ = rollup_proteins(pep, frac, n_boot=20).set_index("protein")
    assert gated.loc["P0", "n_peptides"] == 2     # the Spep-4 peptide dropped
    assert open_.loc["P0", "n_peptides"] == 3     # off by default


def test_rollup_progress_callback_ticks_to_total():
    """rollup_proteins drives the progress callback once per protein group."""
    pep, frac = _make_frames({"sp|P0|X": 0.5, "sp|P1|Y": 0.8}, n_pep=3)
    seen: list[tuple[int, int]] = []
    rollup_proteins(pep, frac, n_boot=0,
                    progress_callback=lambda d, t: seen.append((d, t)))
    assert seen, "callback never fired"
    done, total = seen[-1]
    assert done == total == len(seen)          # one tick per group, ends at total


def test_pooled_method_uses_all_points():
    pep, frac = _make_frames({"sp|P0|X": 0.5}, n_pep=3)  # 3 peptides x 5 t
    out = rollup_proteins(pep, frac, method="pooled", n_boot=50).set_index("protein")
    assert (out["method"] == "pooled").all()
    assert out.loc["P0", "k_deg"] == pytest.approx(0.5, abs=0.02)
    assert out.loc["P0", "n_points"] == 15        # all peptide×timepoint points


def test_rollup_breakdown_and_exclude_mbr():
    """n_mbr/n_metox/n_clean census + the --exclude-mbr filter at the rollup."""
    pep, frac = _make_frames({"sp|P1|X": 0.5}, n_pep=3)  # 3 pep × 5 tp = 15 points
    frac["evidence"] = "q_value"
    frac["metox"] = False
    frac.loc[frac["labeling_time"] == 0.5, "evidence"] = "mbr"   # 3 MBR points
    frac.loc[frac["labeling_time"] == 1.0, "metox"] = True       # 3 Met-Ox points

    inc = rollup_proteins(pep, frac, n_boot=20).set_index("protein")
    assert inc.loc["P1", ["n_mbr", "n_metox", "n_clean"]].tolist() == [3, 3, 9]

    exc = rollup_proteins(pep, frac, n_boot=20, exclude_mbr=True).set_index("protein")
    # the 3 MBR points (t=0.5) are dropped; the 3 Met-Ox (t=1.0) remain.
    assert exc.loc["P1", "n_mbr"] == 0
    assert exc.loc["P1", ["n_metox", "n_clean"]].tolist() == [3, 9]


def test_bioreps_are_independent_refit_points():
    pep, frac = _make_frames({"sp|P0|X": 0.5}, bioreps=(1, 2))
    out = rollup_proteins(pep, frac, n_boot=20).set_index("protein")
    assert out.loc["P0", "n_points"] == 10        # 2 bioreps x 5 timepoints
    assert out.loc["P0", "n_replicates"] == 2
    assert out.loc["P0", "n_timepoints"] == 5
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


def test_r2_admit_gate_rescues_tight_ci_and_floors_railed_fits():
    # The k_cv rescue: keep a low-R² peptide only if its rate constant is tightly
    # determined (k_cv < max) AND R² clears the rescue floor. k_cv = (hi-lo)/(2|k|):
    #   good  R²≥0.8 (primary)                      -> keep
    #   noisy k_cv=(.80-.20)/(2·.5)=0.60 too wide   -> drop
    #   slow  R²=0.65≥floor, k_cv=.002/.02=0.10<0.2 -> keep (flat but well-measured)
    #   railed R²=-0.5 below floor (degenerate k≈0, spuriously tight k_cv=0.10)
    #          -> drop; the floor is what excludes it.
    pep = pd.DataFrame({
        "concat": ["good_2", "noisy_2", "slow_2", "railed_2"],
        "R_squared": [0.95, 0.40, 0.65, -0.50],
        "k_deg": [0.50, 0.50, 0.010, 0.001],
        "ci_lo": [0.48, 0.20, 0.009, 0.0009],
        "ci_hi": [0.52, 0.80, 0.011, 0.0011],
    })
    admitted = _r2_admitted(pep, min_r2=0.8, k_cv_max=0.2, rescue_r2=0.6)
    assert admitted == {"good_2", "slow_2"}


def test_r2_admit_gate_disables_rescue_when_k_cv_max_not_positive():
    # k_cv_max <= 0 turns the rescue off -> only the primary R² ≥ min_r2 survives
    # (and the CI columns aren't even required).
    pep = pd.DataFrame({
        "concat": ["good_2", "slow_2"],
        "R_squared": [0.95, 0.65],
        "k_deg": [0.50, 0.010],
    })
    assert _r2_admitted(pep, min_r2=0.8, k_cv_max=0.0, rescue_r2=0.6) == {"good_2"}


def test_r2_gate_off_by_default_keeps_everything():
    # A deliberately low-R² peptide stays in when min_r2 is None (default).
    pep, frac = _make_frames({"sp|P0|X": 0.5}, n_pep=3)
    pep["R_squared"] = [0.2, 0.95, 0.95]
    # Tight CI on every peptide (k_cv = 0.05): the default rescue is ON, yet the
    # R²=0.2 peptide is still gated out because it sits below the rescue R² floor.
    pep["ci_lo"] = pep["k_deg"] * 0.95
    pep["ci_hi"] = pep["k_deg"] * 1.05
    out_off = rollup_proteins(pep, frac, n_boot=20).set_index("protein")
    assert out_off.loc["P0", "n_peptides"] == 3
    out_on = rollup_proteins(pep, frac, min_r2=0.8, n_boot=20).set_index("protein")
    assert out_on.loc["P0", "n_peptides"] == 2  # the R²=0.2 peptide gated out


def test_rollup_workers_give_identical_result():
    """Process workers (-W) must give the same result as serial — the per-protein
    RNG stream is seeded from the group key, not worker count or order."""
    pep, frac = _make_frames(
        {"sp|P0|X": 0.3, "sp|P1|Y": 0.6, "sp|P2|Z": 0.9}, n_pep=4)
    a = rollup_proteins(pep, frac, n_boot=100, workers=1)
    b = rollup_proteins(pep, frac, n_boot=100, workers=2)
    pd.testing.assert_frame_equal(a, b)


def _make_two_condition_frames(k_by_prot_cond, *, n_pep=3,
                               times=(0, 1, 2, 3, 4, 6, 8, 10)):
    """(peptides, fractions) with a `condition` column and per-(protein,condition)
    kinetics — the substrate for the linear-simple Δk path."""
    pep_rows, frac_rows = [], []
    for pi, (prot, by_cond) in enumerate(k_by_prot_cond.items()):
        for cond, k in by_cond.items():
            for j in range(n_pep):
                concat = f"PEP{pi}_{cond}_{j}_2"
                pep_rows.append({"concat": concat, "protein id": prot,
                                 "k_deg": k, "condition": cond})
                for t in times:
                    theta = 1.0 - np.exp(-k * t)
                    frac_rows.append({
                        "concat": concat, "protein id": prot, "condition": cond,
                        "biological_replicate": 1, "labeling_time": float(t),
                        "fs": theta, "fs_lower": theta - 0.01,
                        "fs_upper": theta + 0.01})
    return pd.DataFrame(pep_rows), pd.DataFrame(frac_rows)


def test_rollup_linear_simple_delta_k():
    """The `linear simple` model: per-condition φ-slope k + a cross-condition Δk
    test, in the PROTEIN_LINEAR_COLUMNS schema."""
    from riana.core.protein import PROTEIN_LINEAR_COLUMNS
    pep, frac = _make_two_condition_frames({
        "sp|P0|X": {"control": 0.05, "atrium": 0.10},   # atrium faster
        "sp|P1|Y": {"control": 0.08, "atrium": 0.08},   # no difference
    })
    res = rollup_proteins(pep, frac, model="linear simple",
                          reference_condition="control", min_peptides=2)
    assert list(res.columns) == PROTEIN_LINEAR_COLUMNS
    assert (res["method"] == "linear simple").all()
    k = res.set_index(["protein", "condition"])["k_deg"]
    assert k[("P0", "control")] == pytest.approx(0.05, abs=0.01)
    assert k[("P0", "atrium")] == pytest.approx(0.10, abs=0.01)
    dk = res.dropna(subset=["delta_k"]).drop_duplicates("protein").set_index("protein")
    # P0: atrium − control ≈ +0.05 and significant; P1: ≈ 0.
    assert dk.loc["P0", "delta_k"] == pytest.approx(0.05, abs=0.01)
    assert dk.loc["P0", "delta_k_p_adj"] < 0.05
    assert abs(dk.loc["P1", "delta_k"]) < 0.01
    # BH p_adj present for every protein that has a Δk.
    assert dk["delta_k_p_adj"].notna().all()


def test_rollup_linear_simple_is_a_valid_model():
    """`linear simple` must not be rejected as an unknown model."""
    pep, frac = _make_two_condition_frames(
        {"sp|P0|X": {"control": 0.05, "atrium": 0.09}})
    rollup_proteins(pep, frac, model="linear simple", min_peptides=2)  # no raise


def test_unknown_model_parsimony_method_raise():
    pep, frac = _make_frames({"sp|P0|X": 0.5})
    with pytest.raises(DataError, match="model"):
        rollup_proteins(pep, frac, model="nope")
    with pytest.raises(DataError, match="parsimony"):
        rollup_proteins(pep, frac, parsimony="razor")
    with pytest.raises(DataError, match="method"):
        rollup_proteins(pep, frac, method="bogus")
