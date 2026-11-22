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
    _PI_SPAN_SIGMA,
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
    assert (out["level"] == "protein").all()
    assert out.loc["P0", "k_deg"] == pytest.approx(0.5, abs=0.02)
    assert out.loc["P0", "R_squared"] > 0.999
    assert out.loc["P0", "n_points"] == 5         # collapsed (biorep, t) cells


def test_rollup_fractions_carry_biological_replicate():
    """`build_rollup_fractions` writes `biological_replicate` per collapsed cell, so
    (experiment, condition, protein, labeling_time, biological_replicate) uniquely
    keys a row — enabling replicate subsetting / biorep-split modelling from the file."""
    from riana.core.protein import build_rollup_fractions
    pep, frac = _make_frames({"sp|P0|X": 0.5}, bioreps=(1, 2, 3))
    out = rollup_proteins(pep, frac, n_boot=20)
    rf = build_rollup_fractions(out)
    assert "biological_replicate" in rf.columns
    assert set(rf["biological_replicate"]) == {1, 2, 3}
    # one row per (biorep, timepoint) cell → the 5-tuple key is unique
    key = ["experiment", "condition", "protein", "labeling_time", "biological_replicate"]
    assert not rf.duplicated(key).any()
    assert len(rf) == 3 * 5   # 3 bioreps × 5 timepoints


def test_rollup_fractions_biorep_on_exponential_and_linear_paths():
    """Both the exponential collapse and the linear path carry biorep."""
    from riana.core.protein import build_rollup_fractions
    pep, frac = _make_frames({"sp|P0|X": 0.5}, bioreps=(1, 2, 3))
    expo = build_rollup_fractions(rollup_proteins(pep, frac, n_boot=20))
    assert set(expo["biological_replicate"]) == {1, 2, 3}
    # linear: two conditions so a Δk substrate exists
    pep2, frac2 = _make_frames({"sp|P0|X": 0.5}, bioreps=(1, 2, 3))
    for df in (pep2, frac2):
        df["condition"] = "control"
    pepK, fracK = _make_frames({"sp|P0|X": 0.9}, bioreps=(1, 2, 3))
    for df in (pepK, fracK):
        df["condition"] = "ko"
    pep_all = pd.concat([pep2, pepK], ignore_index=True)
    frac_all = pd.concat([frac2, fracK], ignore_index=True)
    lin = build_rollup_fractions(
        rollup_proteins(pep_all, frac_all, model="linear",
                        reference_condition="control", n_boot=20))
    assert "biological_replicate" in lin.columns
    assert set(lin["biological_replicate"]) == {1, 2, 3}


def test_k_cv_admitted_single_timepoint_gate():
    """The single-timepoint k_cv gate admits k_cv < threshold and rejects a wide CI
    or a NaN CI (single-point / non-converged)."""
    pep = pd.DataFrame([
        {"experiment": "e", "condition": "c",
         "concat": "A", "k_deg": 0.02, "ci_lo": 0.019, "ci_hi": 0.021},  # k_cv=0.05 keep
        {"experiment": "e", "condition": "c",
         "concat": "B", "k_deg": 0.02, "ci_lo": 0.005, "ci_hi": 0.05},   # k_cv=1.1 reject
        {"experiment": "e", "condition": "c",
         "concat": "C", "k_deg": 0.02, "ci_lo": np.nan, "ci_hi": np.nan},  # NaN reject
    ])
    assert _k_cv_admitted(pep, 0.2) == {("e", "c", "A")}


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
    auto = rollup_proteins(pep, frac, model="exponential", min_points=1,
                           min_peptides=1, n_boot=20).set_index("protein")
    assert int(auto.loc["P1", "n_peptides"]) == 1     # only the replicated A survives
    off = rollup_proteins(pep, frac, model="exponential", min_points=1, min_peptides=1,
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
    out = rollup_proteins(pep, frac, model="exponential", n_boot=20).set_index("protein")
    assert out.loc["P0", "n_points"] == 10        # 2 bioreps x 5 timepoints
    assert out.loc["P0", "n_replicates"] == 2
    assert out.loc["P0", "n_timepoints"] == 5
    assert out.loc["P0", "k_deg"] == pytest.approx(0.5, abs=0.02)


def test_weighted_theta_favours_tighter_ci():
    fs = np.array([0.5, 0.6])
    lo = np.array([0.49, 0.40])
    hi = np.array([0.51, 0.80])
    theta, var, df = _weighted_theta(fs, lo, hi)
    assert theta == pytest.approx(0.5, abs=0.01)      # the tight-CI peptide dominates
    # var = 1/Σ(1/σ²) is smaller than the tighter peptide's own σ² (info adds)
    sig_tight = (0.51 - 0.49) / _PI_SPAN_SIGMA
    assert 0.0 < var < sig_tight ** 2
    assert df == pytest.approx(2.0)                    # n_pep fallback (dfs not given)


def test_weighted_theta_falls_back_to_unweighted_when_no_ci():
    fs = np.array([0.4, 0.6])
    nan = np.array([np.nan, np.nan])
    theta, var, df = _weighted_theta(fs, nan, nan)
    assert theta == pytest.approx(0.5)                # unweighted mean, θ unchanged
    assert np.isnan(var) and np.isnan(df)             # no σ info → undefined variance


def test_weighted_theta_satterthwaite_df():
    """eff_df = (Σw)²/Σ(w²/d): a lone peptide → its own df; equal peptides → n·d."""
    fs = np.array([0.5, 0.5])
    lo = np.array([0.49, 0.49])
    hi = np.array([0.51, 0.51])            # identical σ → equal weights
    _, _, df = _weighted_theta(fs, lo, hi, dfs=np.array([6.0, 6.0]))
    assert df == pytest.approx(12.0)       # two equal peptides at df 6 → 12
    _, _, df1 = _weighted_theta(fs[:1], lo[:1], hi[:1], dfs=np.array([6.0]))
    assert df1 == pytest.approx(6.0)       # single peptide → its own df


def test_min_peptides_filters_small_proteins():
    pep, frac = _make_frames({"sp|P0|X": 0.5}, n_pep=1)
    out = rollup_proteins(pep, frac, min_peptides=2, n_boot=20)
    assert out.empty or "P0" not in set(out["protein"])
    # default (since 2026-11-17): one good peptide reports the protein
    out = rollup_proteins(pep, frac, n_boot=20)
    assert "P0" in set(out["protein"])


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
        "experiment": ["e"] * 4,
        "condition": ["c"] * 4,
        "concat": ["good_2", "noisy_2", "slow_2", "railed_2"],
        "R_squared": [0.95, 0.40, 0.65, -0.50],
        "k_deg": [0.50, 0.50, 0.010, 0.001],
        "ci_lo": [0.48, 0.20, 0.009, 0.0009],
        "ci_hi": [0.52, 0.80, 0.011, 0.0011],
    })
    admitted = _r2_admitted(pep, min_r2=0.8, k_cv_max=0.2, rescue_r2=0.6)
    assert admitted == {("e", "c", "good_2"), ("e", "c", "slow_2")}


def test_r2_admit_gate_disables_rescue_when_k_cv_max_not_positive():
    # k_cv_max <= 0 turns the rescue off -> only the primary R² ≥ min_r2 survives
    # (and the CI columns aren't even required).
    pep = pd.DataFrame({
        "experiment": ["e", "e"],
        "condition": ["c", "c"],
        "concat": ["good_2", "slow_2"],
        "R_squared": [0.95, 0.65],
        "k_deg": [0.50, 0.010],
    })
    assert _r2_admitted(pep, min_r2=0.8, k_cv_max=0.0, rescue_r2=0.6) == {("e", "c", "good_2")}


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


def test_rollup_gates_are_per_condition_not_leaked_across_conditions():
    """A curation gate is keyed on (experiment, condition, concat): a peptide that
    clears the R² gate in one condition but FAILS it in another is dropped only in
    the failing condition — its good-condition pass must not drag its bad-condition
    rows into the other condition's rollup (which would corrupt the per-condition k
    and the two-condition Δk contrast). Regression for the condition-blind ``.isin``
    on bare ``concat`` (which admitted a peptide in every condition once it passed
    in one)."""
    k = 0.1
    times = (0, 1, 2, 3, 4, 6, 8, 10)
    # PEP passes in ctrl (R²=0.95) but fails in drug (R²=0.30, below the 0.6 rescue
    # floor even at tight k_cv); QEP passes in both.
    r2_by = {("PEP_2", "ctrl"): 0.95, ("PEP_2", "drug"): 0.30,
             ("QEP_2", "ctrl"): 0.95, ("QEP_2", "drug"): 0.95}
    pep_rows, frac_rows = [], []
    for (concat, cond), r2 in r2_by.items():
        pep_rows.append({"concat": concat, "protein id": "sp|P0|X", "condition": cond,
                         "k_deg": k, "R_squared": r2,
                         "ci_lo": k * 0.98, "ci_hi": k * 1.02})   # k_cv=0.02 (tight)
        for t in times:
            theta = 1.0 - np.exp(-k * t)
            frac_rows.append({"concat": concat, "protein id": "sp|P0|X",
                              "condition": cond, "biological_replicate": 1,
                              "labeling_time": float(t), "fs": theta,
                              "fs_lower": theta - 0.01, "fs_upper": theta + 0.01})
    pep, frac = pd.DataFrame(pep_rows), pd.DataFrame(frac_rows)
    out = rollup_proteins(pep, frac, min_r2=0.8, min_peptides=1, min_points=3,
                          n_boot=20, peptide_admission="own").set_index("condition")
    assert int(out.loc["ctrl", "n_peptides"]) == 2   # both peptides pass in ctrl
    assert int(out.loc["drug", "n_peptides"]) == 1   # only QEP; PEP must NOT leak in


def _two_condition_gate_frames():
    """PEP passes the R² gate in ctrl only (0.95 vs 0.30); QEP passes in both."""
    k, times = 0.1, (0, 1, 2, 3, 4, 6, 8, 10)
    r2_by = {("PEP_2", "ctrl"): 0.95, ("PEP_2", "drug"): 0.30,
             ("QEP_2", "ctrl"): 0.95, ("QEP_2", "drug"): 0.95}
    pep_rows, frac_rows = [], []
    for (concat, cond), r2 in r2_by.items():
        pep_rows.append({"concat": concat, "protein id": "sp|P0|X", "condition": cond,
                         "k_deg": k, "R_squared": r2, "ci_lo": k * 0.98, "ci_hi": k * 1.02})
        for t in times:
            theta = 1.0 - np.exp(-k * t)
            frac_rows.append({"concat": concat, "protein id": "sp|P0|X", "condition": cond,
                              "biological_replicate": 1, "labeling_time": float(t),
                              "fs": theta, "fs_lower": theta - 0.01, "fs_upper": theta + 0.01})
    return pd.DataFrame(pep_rows), pd.DataFrame(frac_rows)


def test_peptide_admission_policies_any_own_all():
    """--peptide-admission maps a peptide's per-condition gate result across conditions:
    'own' = strictly per-condition, 'any' = passed-in-one kept-in-all (the old leak, now
    opt-in), 'all' = kept only if it passes in every condition it appears in."""
    pep, frac = _two_condition_gate_frames()

    def n_by_cond(policy):
        out = rollup_proteins(pep, frac, min_r2=0.8, min_peptides=1, min_points=3,
                              n_boot=20, peptide_admission=policy).set_index("condition")
        return {c: int(out.loc[c, "n_peptides"]) for c in ("ctrl", "drug")}

    assert n_by_cond("own") == {"ctrl": 2, "drug": 1}   # PEP kept only where it passes
    assert n_by_cond("any") == {"ctrl": 2, "drug": 2}   # PEP passed ctrl -> kept in both
    assert n_by_cond("all") == {"ctrl": 1, "drug": 1}   # PEP fails drug -> dropped in both


def test_fit_kdeg_single_point_returns_nan_ci():
    """_fit_kdeg on a single collapsed (t, θ) point returns a NaN CI, not a zero-width
    one: the 1-point fit passes through the point exactly, so every bootstrap resample
    gives the identical k (undefined uncertainty, not zero — which would give k_cv=0 and
    sail through the --k-cv gate). k itself is still solved directly."""
    from riana.core import models
    from riana.core.protein import _fit_kdeg
    k, r2, lo, hi = _fit_kdeg(
        np.array([24.0]), np.array([0.38]),
        model_fn=models.one_exponent, kinetic_kwargs=dict(a_0=0.0, a_max=1.0),
        n_boot=50, boot_ci_pct=(5.0, 95.0), rng=np.random.default_rng(0))
    assert np.isfinite(k)          # k solved directly at the single timepoint
    assert np.isnan(r2)            # R² N/A at one timepoint
    assert np.isnan(lo) and np.isnan(hi)   # undefined CI, not zero-width


def test_unknown_peptide_admission_raises():
    pep, frac = _two_condition_gate_frames()
    with pytest.raises(DataError, match="peptide_admission"):
        rollup_proteins(pep, frac, min_r2=0.8, peptide_admission="bogus")


@pytest.mark.parametrize("model, level, expected", [
    ("exponential", "protein", "own"), ("calibration", "protein", "own"),
    ("linear", "protein", "all"), ("linear", "peptide", "own")])
def test_auto_admission_resolves_by_model_and_level(model, level, expected):
    from riana.core.protein import resolve_admission
    assert resolve_admission("auto", model, level) == expected
    assert resolve_admission("any", model, level) == "any"     # explicit passes through


@pytest.mark.parametrize("model, ctrl_peptides", [("exponential", 2), ("linear", 1)])
def test_default_admission_is_own_for_exponential_and_all_for_linear(model, ctrl_peptides):
    """The default 'auto': the linear Δk (protein level) keeps only peptides passing in
    both conditions; the exponential k (no Δk) keeps each condition's own passes, so a
    condition's k does not depend on the other conditions in the table."""
    pep, frac = _two_condition_gate_frames()   # PEP passes ctrl only; QEP passes both
    kw = dict(model=model, min_r2=0.8, min_peptides=1, n_boot=20)
    if model == "linear":
        kw["reference_condition"] = "ctrl"
    default = rollup_proteins(pep, frac, **kw).set_index("condition")
    assert int(default.loc["ctrl", "n_peptides"]) == ctrl_peptides
    assert int(default.loc["drug", "n_peptides"]) == 1


def test_linear_rollup_uses_only_the_named_pair():
    """With a third condition in the table, a linear Δk equals the one made from the
    table cut to the pair: the third condition's peptides do not enter admission (PEP
    fails only in 'other'), the shared bands or the joint fit, and its rows are dropped."""
    pep, frac = _two_condition_gate_frames()
    k, times = 0.1, (0, 1, 2, 3, 4, 6, 8, 10)
    rng = np.random.default_rng(3)
    for concat, r2 in (("PEP_2", 0.30), ("QEP_2", 0.95)):
        pep = pd.concat([pep, pd.DataFrame([{
            "concat": concat, "protein id": "sp|P0|X", "condition": "other",
            "k_deg": k, "R_squared": r2, "ci_lo": k * 0.98, "ci_hi": k * 1.02}])])
        theta = 1.0 - np.exp(-0.2 * np.array(times)) + rng.normal(0, 0.03, len(times))
        frac = pd.concat([frac, pd.DataFrame({
            "concat": concat, "protein id": "sp|P0|X", "condition": "other",
            "biological_replicate": 1, "labeling_time": np.array(times, float),
            "fs": theta, "fs_lower": theta - 0.05, "fs_upper": theta + 0.05})])
    pep.loc[(pep.concat == "PEP_2") & (pep.condition == "drug"), "R_squared"] = 0.95
    kw = dict(model="linear", min_r2=0.8, reference_condition="ctrl",
              test_condition="drug")
    full = rollup_proteins(pep, frac, **kw)
    pair = rollup_proteins(pep[pep.condition != "other"], frac[frac.condition != "other"],
                           **kw)
    assert set(full["condition"]) == {"ctrl", "drug"}
    pd.testing.assert_frame_equal(full, pair)
    assert int(full["n_peptides"].iloc[0]) == 2          # PEP kept: 'other' is not asked


def test_linear_rollup_with_three_conditions_needs_the_pair():
    pep, frac = _two_condition_gate_frames()
    frac = pd.concat([frac, frac.assign(condition="other")])
    pep = pd.concat([pep, pep.assign(condition="other")])
    with pytest.raises(DataError, match="Name the pair"):
        rollup_proteins(pep, frac, model="linear", min_r2=0.8)
    with pytest.raises(DataError, match="Name the pair"):
        rollup_proteins(pep, frac, model="linear", min_r2=0.8, reference_condition="ctrl")
    with pytest.raises(DataError, match="not among the conditions"):
        rollup_proteins(pep, frac, model="linear", min_r2=0.8, reference_condition="ctrl",
                        test_condition="drgu")
    with pytest.raises(DataError, match="--model exponential"):
        rollup_proteins(pep, frac, min_r2=0.8)          # the default (linear) refuses too
    assert set(rollup_proteins(pep, frac, model="exponential", min_r2=0.8,
                               n_boot=10)["condition"]) == {
        "ctrl", "drug", "other"}                         # exponential: no contrast


def test_dimethyl_channel_forms_are_one_peptide_across_conditions():
    """A dimethyl light form (control channel) and heavy form (treatment channel) of a
    peptide are the same peptide to the rollup: admission 'all' pairs them (the heavy
    form failing drops the light one), as for an unlabelled peptide."""
    pep, frac = _two_condition_gate_frames()
    tag = {"ctrl": "[UNIMOD:36]", "drug": "[UNIMOD:330]"}
    for df in (pep, frac):
        df["concat"] = [f"{tag[c]}{s[:-2]}K{tag[c]}_2"
                        for s, c in zip(df["concat"], df["condition"])]
    out = rollup_proteins(pep, frac, model="linear", min_r2=0.8,
                          reference_condition="ctrl").set_index("condition")
    assert int(out.loc["ctrl", "n_peptides"]) == 1       # PEP's heavy form fails drug
    own = rollup_proteins(pep, frac, model="linear", min_r2=0.8, reference_condition="ctrl",
                          peptide_admission="own").set_index("condition")
    assert int(own.loc["ctrl", "n_peptides"]) == 2


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
                concat = f"PEP{pi}_{j}_2"   # the same peptide in every condition
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
    """The `linear` model: per-condition φ-slope k + a cross-condition Δk
    test, in the PROTEIN_LINEAR_COLUMNS schema."""
    from riana.core.protein import PROTEIN_LINEAR_COLUMNS
    pep, frac = _make_two_condition_frames({
        "sp|P0|X": {"control": 0.05, "atrium": 0.10},   # atrium faster
        "sp|P1|Y": {"control": 0.08, "atrium": 0.08},   # no difference
    })
    res = rollup_proteins(pep, frac, model="linear",
                          reference_condition="control", min_peptides=2)
    assert list(res.columns) == PROTEIN_LINEAR_COLUMNS
    assert (res["level"] == "protein").all()
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


def _composition_frames():
    """One protein, two peptides that differ in k (0.05, 0.20) with the SAME curve in
    both conditions (true Δk = 0). Their bands swap between the conditions: in X the
    slow peptide is tight and the fast one wide, in Y the reverse — the per-condition
    precision difference that makes a collapse mix the peptides in different proportions."""
    pep_rows, frac_rows = [], []
    for cond in ("X", "Y"):
        for concat, k in (("SLOW_2", 0.05), ("FAST_2", 0.20)):
            tight = (concat == "SLOW_2") == (cond == "X")
            half = 0.005 if tight else 0.05
            pep_rows.append({"concat": concat, "protein id": "sp|P0|X", "k_deg": k,
                             "condition": cond})
            for t in (1, 2, 3, 4, 6, 8, 10):
                theta = 1.0 - np.exp(-k * t)
                frac_rows.append({
                    "concat": concat, "protein id": "sp|P0|X", "condition": cond,
                    "biological_replicate": 1, "labeling_time": float(t),
                    "fs": theta, "fs_lower": theta - half, "fs_upper": theta + half})
    return pd.DataFrame(pep_rows), pd.DataFrame(frac_rows)


def test_linear_collapse_shares_band_levels_so_composition_gives_no_delta_k(monkeypatch):
    """Candidate C (reports/2026-11-14 §5): the linear collapse pools each peptide's band
    level over the conditions, so swapping which peptide is precise in which condition
    no longer moves the collapsed θ — Δk stays at the true 0. Without the shared level
    the same data give a large spurious Δk."""
    import riana.core.protein as protein_mod
    pep, frac = _composition_frames()

    def delta_k() -> float:
        res = rollup_proteins(pep, frac, model="linear", reference_condition="X",
                              min_peptides=2)
        return float(res.dropna(subset=["delta_k"])["delta_k"].iloc[0])

    assert abs(delta_k()) < 1e-9
    monkeypatch.setattr(protein_mod, "_shared_bands", lambda f: f)
    assert abs(delta_k()) > 0.05


def test_shared_bands_pools_the_level_and_keeps_the_time_profile():
    from riana.core.protein import _shared_bands
    _, frac = _composition_frames()
    frac["experiment"] = "e"
    # A peptide seen in X only, with a time-varying band: left exactly as is.
    solo = frac[(frac.concat == "SLOW_2") & (frac.condition == "X")].assign(concat="SOLO_2")
    solo = solo.assign(fs_upper=solo.fs + 0.002 * solo.labeling_time)
    out = _shared_bands(pd.concat([frac, solo], ignore_index=True))
    s2 = ((out.fs_upper - out.fs_lower) / _PI_SPAN_SIGMA) ** 2
    lvl = s2.groupby([out.concat, out.condition]).mean()
    for concat in ("SLOW_2", "FAST_2"):          # one level in both conditions
        assert lvl[(concat, "X")] == pytest.approx(lvl[(concat, "Y")], rel=1e-12)
    pd.testing.assert_frame_equal(
        out[out.concat == "SOLO_2"].reset_index(drop=True), solo.reset_index(drop=True))
    assert (out.fs == pd.concat([frac, solo], ignore_index=True).fs).all()


def test_shared_bands_leaves_a_single_condition_rollup_unchanged():
    from riana.core.protein import _shared_bands
    _, frac = _composition_frames()
    one = frac[frac.condition == "X"].assign(experiment="e")
    assert _shared_bands(one) is one


def test_rollup_linear_simple_is_a_valid_model():
    """`linear` must not be rejected as an unknown model."""
    pep, frac = _make_two_condition_frames(
        {"sp|P0|X": {"control": 0.05, "atrium": 0.09}})
    rollup_proteins(pep, frac, model="linear", min_peptides=2)  # no raise


def test_unknown_model_parsimony_level_raise():
    pep, frac = _make_frames({"sp|P0|X": 0.5})
    with pytest.raises(DataError, match="model"):
        rollup_proteins(pep, frac, model="nope")
    with pytest.raises(DataError, match="parsimony"):
        rollup_proteins(pep, frac, parsimony="razor")
    with pytest.raises(DataError, match="level"):
        rollup_proteins(pep, frac, level="bogus")
    with pytest.raises(DataError, match="needs model 'linear'"):
        rollup_proteins(pep, frac, model="exponential", level="peptide")  # no peptide level


def test_rollup_fractions_carry_individual_when_the_fit_did():
    """``individual`` (the animal) is mapped onto each collapsed point by
    (experiment, condition, labeling_time, biological_replicate); an ID column a
    TSV round-trip turned into floats (gaps → NaN, 843 → 843.0) reads back as text."""
    from riana.core.protein import build_rollup_fractions
    pep, frac = _make_frames({"sp|P0|X": 0.5}, bioreps=(1, 2, 3))
    assert "individual" not in build_rollup_fractions(
        rollup_proteins(pep, frac, n_boot=20)).columns
    ids = {1: 843.0, 2: 844.0, 3: np.nan}
    frac = frac.assign(individual=frac["biological_replicate"].map(ids))
    for kw in ({}, {"model": "linear"}):
        rf = build_rollup_fractions(rollup_proteins(pep, frac, n_boot=20, **kw))
        cols = list(rf.columns)
        assert cols[cols.index("biological_replicate") + 1] == "individual"
        assert (rf["individual"] == rf["biological_replicate"].map(
            {1: "843", 2: "844", 3: ""})).all()


# --------------------------------------------------------------------------- #
# Model name + what the fit's fs holds (rising precursor, 2026-11-13)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("old, new", [("simple", "exponential"), ("linear simple", "linear")])
def test_pre_1_2_model_names_are_refused_with_the_new_name(old, new):
    """No aliases since 1.2.0 (breaking change): the old names are refused, the error
    names the new one."""
    pep, frac = _make_frames({"sp|P0|X": 0.5})
    with pytest.raises(DataError, match=f"renamed '{new}' in Riana 1.2.0"):
        rollup_proteins(pep, frac, model=old)


def test_model_from_header_reads_files_written_before_1_2():
    """Reading old outputs (the app's Load results) still maps the old names."""
    from riana.constants import model_from_header

    assert model_from_header("simple") == "exponential"
    assert model_from_header("linear simple") == "linear"
    assert model_from_header("linear") == "linear"


@pytest.mark.parametrize("model", ["guan", "fornasiero"])
def test_removed_models_are_refused_with_a_pointer_to_kp(model):
    pep, frac = _make_frames({"sp|P0|X": 0.5})
    with pytest.raises(DataError, match=r"removed in Riana 1\.2\.0.*--kp"):
        rollup_proteins(pep, frac, model=model)


def test_check_fs_basis_refuses_an_unknown_basis():
    from riana.core.protein import check_fs_basis

    with pytest.raises(DataError, match="unknown fs_basis"):
        check_fs_basis("phi")
    for basis in ("theta", "f", "mixed"):   # every rollup model applies to each
        check_fs_basis(basis)


@pytest.mark.parametrize("fit_model, refused", [
    ("guan", True), ("fornasiero", True), ("simple", False),
    ("exponential", False), (None, False)])
def test_check_fit_model_refuses_fits_from_removed_models(tmp_path, fit_model, refused):
    from riana.core.protein import check_fit_model

    f = tmp_path / "riana_fit_fractions.txt"
    head = "# riana 1.1.0\n" + (f"# model {fit_model}\n" if fit_model else "")
    f.write_text(head + "concat\tfs\n")
    if refused:
        with pytest.raises(DataError, match=f"made with --model {fit_model}"):
            check_fit_model(f)
    else:
        check_fit_model(f)


def test_fit_basis_header_reads_the_fit_header(tmp_path):
    from riana.core.protein import fit_basis_header

    old = tmp_path / "old.txt"
    old.write_text("# riana 1.1.0\nconcat\tfs\n")
    assert fit_basis_header(old) == {"fs_basis": "theta"}
    new = tmp_path / "new.txt"
    new.write_text("# riana 1.2.0\n# fs_basis f\n# k_p_applied A/ctrl=0.1\nconcat\tfs\n")
    assert fit_basis_header(new) == {"fs_basis": "f", "k_p_applied": "A/ctrl=0.1"}


def test_protein_level_delta_k_needs_a_shared_peptide():
    """A protein whose two conditions share no peptide gets k per condition but no Δk
    (it would compare different peptides), as at the peptide level; BH runs over the
    proteins that have one. The interaction input drops that experiment for it."""
    from riana.core.protein import build_rollup_fractions, interaction_points
    pep, frac = _make_two_condition_frames({
        "sp|P0|X": {"control": 0.05, "atrium": 0.10},
        "sp|P1|Y": {"control": 0.08, "atrium": 0.12},
    })
    split = lambda df: df.assign(concat=np.where(   # P1: different peptides per condition
        df["protein id"].str.contains("P1"), df["concat"] + "_" + df["condition"], df["concat"]))
    pep, frac = split(pep), split(frac)
    res = rollup_proteins(pep, frac, model="linear", reference_condition="control")
    p1 = res[res["protein"] == "P1"]
    assert len(p1) == 2 and p1["k_deg"].notna().all()
    assert p1["delta_k"].isna().all() and p1["delta_k_p_adj"].isna().all()
    p0 = res[res["protein"] == "P0"]
    assert p0["delta_k"].iloc[0] == pytest.approx(0.05, abs=1e-3)
    pts = interaction_points(res, build_rollup_fractions(res))
    assert set(pts["protein"]) == {"P0"}


@pytest.mark.parametrize("model", ["linear", "exponential"])
def test_rollup_attrs_are_shared_not_copied_into_slices(model):
    """pandas deep-copies ``attrs`` into every slice, row and group; the per-protein
    points must be shared instead, or per-protein loops over a large rollup (the
    interaction filter, --plot) spend ~1 s per protein copying them."""
    pep, frac = _make_two_condition_frames({
        "sp|P0|X": {"control": 0.05, "atrium": 0.10},
        "sp|P1|Y": {"control": 0.08, "atrium": 0.12},
    })
    res = rollup_proteins(pep, frac, model=model, reference_condition="control")
    points = res.attrs["protein_points"]
    assert points
    _, row = next(res.iterrows())
    for derived in (res[res["k_deg"].notna()], row,
                    next(iter(res.groupby("protein")))[1]):
        assert derived.attrs["protein_points"] is points
        assert derived.attrs["point_individuals"] is res.attrs["point_individuals"]
