"""Peptide-level linear model (``rollup --model linear --level peptide``): E+time / EJ+time.

The numerics are checked against the evaluation's cached results on real data (report
``2026-11-14_wls_var_D_E_evaluation.md`` §5 decision note); these tests pin the behaviour on
small synthetic designs: recovery, the random-effects path and df, pairing, the fallbacks,
the condition handling and the output schema.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from riana.core.linear_model import LINEAR_COLUMNS
from riana.core.peptide_model import (
    PEPTIDE_INTERACTION_COLUMNS, PEPTIDE_LEVEL_COLUMNS, fit_peptide_deltak,
    fit_peptide_interaction,
)
from riana.exceptions import DataError

_TIMES = (1.0, 2.0, 4.0, 7.0, 10.0, 14.0)


def _points(k_by_cond, *, n_pep=4, n_rep=3, times=_TIMES, paired=False,
            longitudinal=False, experiments=("e",), exp_ratio=None, seed=0,
            protein="P0", noise=0.01):
    """Peptide points for one protein: per condition a k, per peptide a slope offset,
    per animal and per labeling time a random slope. ``paired``: each animal carries
    every condition; ``longitudinal``: each animal is sampled at every time;
    ``experiments``: strata of the same animals (tissues), ``exp_ratio`` scaling the
    test condition's k per experiment."""
    rng = np.random.default_rng(seed)
    offsets = rng.normal(0, 0.01, n_pep)
    offsets -= offsets.mean()
    u_time = {t: rng.normal(0, 0.004) for t in times}
    conds = list(k_by_cond)
    rows = []
    for t in times:
        for r in range(n_rep):
            shared = rng.normal(0, 0.003)
            for c in conds:
                animal = (f"a{r}" if longitudinal
                          else f"t{t:g}-r{r}" if paired else f"{c}-t{t:g}-r{r}")
                u_anim = shared if (paired or longitudinal) else rng.normal(0, 0.003)
                for e in experiments:
                    k_c = k_by_cond[c] * ((exp_ratio or {}).get(e, 1.0)
                                          if c == conds[-1] else 1.0)
                    for j in range(n_pep):
                        k = k_c + offsets[j] + u_anim + u_time[t]
                        th = 1 - np.exp(-k * t) + rng.normal(0, noise)
                        rows.append({
                            "experiment": e, "condition": c, "protein": protein,
                            "concat": f"PEP{j}_2", "labeling_time": t,
                            "biological_replicate": r + 1, "individual": animal,
                            "fs": th, "fs_lower": th - 0.02, "fs_upper": th + 0.02})
    return pd.DataFrame(rows)


def test_recovers_delta_k_with_the_time_model_and_the_output_schema():
    f = _points({"ctrl": 0.10, "ko": 0.15})
    out = fit_peptide_deltak(f, reference_condition="ctrl", test_condition="ko")
    assert list(out.columns) == [*LINEAR_COLUMNS, *PEPTIDE_LEVEL_COLUMNS]
    assert set(out["condition"]) == {"ctrl", "ko"}
    row = out.iloc[0]
    assert row["delta_k"] == pytest.approx(0.05, abs=0.01)
    assert row["delta_k_p_adj"] < 0.05
    assert row["random_effects"] == "time" and not row["paired"]
    n_anim = 2 * len(_TIMES) * 3
    assert row["n_animals"] == n_anim
    assert row["delta_k_df"] == n_anim - len(_TIMES) - 1      # animals − times − 1
    k = out.set_index("condition")
    assert k.loc["ctrl", "ci_lo"] < 0.10 < k.loc["ctrl", "ci_hi"]
    assert out["R_squared"].isna().all()


def test_pairing_is_detected_from_shared_animals():
    f = _points({"LV": 0.08, "RV": 0.08}, paired=True, seed=1)
    row = fit_peptide_deltak(f, reference_condition="LV", test_condition="RV").iloc[0]
    assert row["paired"]
    assert row["delta_k_df"] == row["n_animals"] - 1            # paired: animals − 1
    assert abs(row["delta_k"]) < 0.01


def test_longitudinal_design_falls_back_to_the_animal_model():
    """An animal sampled at several times is not the cross-sectional time model."""
    f = _points({"ctrl": 0.10, "ko": 0.12}, longitudinal=True, seed=2)
    row = fit_peptide_deltak(f, reference_condition="ctrl", test_condition="ko").iloc[0]
    assert row["random_effects"] in ("animal", "none")


def test_a_protein_with_one_condition_gets_k_without_delta_k():
    f = _points({"ctrl": 0.10, "ko": 0.15})
    f = pd.concat([f, _points({"ctrl": 0.07}, protein="P1", seed=3)], ignore_index=True)
    out = fit_peptide_deltak(f, reference_condition="ctrl", test_condition="ko")
    p1 = out[out["protein"] == "P1"]
    assert list(p1["condition"]) == ["ctrl"]
    assert p1["delta_k"].isna().all()
    assert p1["k_deg"].iloc[0] == pytest.approx(0.07, abs=0.01)


def test_more_than_two_conditions_need_the_pair_named():
    f = _points({"a": 0.10, "b": 0.12, "c": 0.14})
    with pytest.raises(DataError, match="Name the pair"):
        fit_peptide_deltak(f)
    out = fit_peptide_deltak(f, reference_condition="a", test_condition="c")
    assert set(out["condition"]) == {"a", "c"}                   # only the pair is fitted
    assert out["delta_k"].iloc[0] == pytest.approx(0.04, abs=0.01)
    with pytest.raises(DataError, match="not among the conditions"):
        fit_peptide_deltak(f, reference_condition="a", test_condition="zz")


def test_min_peptides_and_workers():
    f = _points({"ctrl": 0.10, "ko": 0.15})
    one = f[f["concat"] == "PEP0_2"]
    assert fit_peptide_deltak(one, min_peptides=2).empty
    f2 = pd.concat([f, _points({"ctrl": 0.06, "ko": 0.06}, protein="P1", seed=4)],
                   ignore_index=True)
    serial = fit_peptide_deltak(f2, reference_condition="ctrl")
    parallel = fit_peptide_deltak(f2, reference_condition="ctrl", workers=2)
    pd.testing.assert_frame_equal(serial, parallel)


def test_rollup_proteins_peptide_level():
    from riana.core.protein import PROTEIN_PEPTIDE_COLUMNS, rollup_proteins
    f = _points({"ctrl": 0.10, "ko": 0.15}).rename(columns={"protein": "protein id"})
    f["protein id"] = "sp|P0|X"
    pep = (f[["experiment", "condition", "concat", "protein id"]].drop_duplicates()
           .assign(k_deg=0.1))
    res = rollup_proteins(pep, f, model="linear", level="peptide",
                          reference_condition="ctrl")
    assert list(res.columns) == PROTEIN_PEPTIDE_COLUMNS
    assert (res["level"] == "peptide").all()
    assert res["delta_k"].iloc[0] == pytest.approx(0.05, abs=0.01)
    assert res.attrs["protein_points"]          # the collapsed display points
    with pytest.raises(DataError, match="needs model 'linear'"):
        rollup_proteins(pep, f, level="peptide")


def test_peptide_interaction_on_tissues_of_the_same_animals():
    f = _points({"ctrl": 0.10, "ko": 0.13}, experiments=("LA", "LV"),
                exp_ratio={"LA": 1.0, "LV": 1.0}, seed=5)
    out = fit_peptide_interaction(f, reference_condition="ctrl", test_condition="ko")
    from riana.core.interaction_model import INTERACTION_COLUMNS
    assert list(out.columns) == [*INTERACTION_COLUMNS, *PEPTIDE_INTERACTION_COLUMNS]
    assert set(out["experiment"]) == {"LA", "LV"}
    row = out.iloc[0]
    assert row["random_effects"] == "time"
    assert row["mean_delta_k"] == pytest.approx(0.03, abs=0.01)
    assert row["n_experiments"] == 2
    assert out.attrs["animal_key"] == "individual"
