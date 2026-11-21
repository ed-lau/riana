"""Tests for ``riana.core.rollup_plot`` (the CLI ``rollup --plot`` PDF export)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from riana.core.rollup_plot import render_linear_pdf


def _frames(experiments):
    """Synthetic linear-simple rollup tables over the given experiments."""
    prot_rows, frac_rows = [], []
    for exp in experiments:
        for cond, k in (("control", 0.03), ("knockout", 0.06)):
            prot_rows.append(dict(
                experiment=exp, condition=cond, protein="P1", method="weighted",
                k_deg=k, ci_lo=k * 0.9, ci_hi=k * 1.1, R_squared=0.99,
                delta_k=0.03, delta_k_se=0.005, delta_k_p=1e-4, delta_k_p_adj=1e-3))
            for t in (0, 1, 3, 7, 14):
                theta = 1 - np.exp(-k * t)
                frac_rows.append(dict(experiment=exp, condition=cond, protein="P1",
                                      labeling_time=t, fs=theta,
                                      fs_var=1e-4, fs_df=10.0))
    return pd.DataFrame(prot_rows), pd.DataFrame(frac_rows)


def _is_pdf(path):
    return path.exists() and path.read_bytes()[:4] == b"%PDF"


def test_render_single_experiment_writes_one_unsuffixed_pdf(tmp_path):
    proteins, fractions = _frames(["exp1"])
    written = render_linear_pdf(proteins, fractions, tmp_path)
    assert written == [tmp_path / "riana_rollup_curves.pdf"]
    assert _is_pdf(written[0])


def test_render_multi_experiment_writes_one_pdf_each(tmp_path):
    proteins, fractions = _frames(["LA", "LV"])
    written = render_linear_pdf(proteins, fractions, tmp_path)
    names = sorted(p.name for p in written)
    assert names == ["riana_rollup_curves_LA.pdf", "riana_rollup_curves_LV.pdf"]
    assert all(_is_pdf(p) for p in written)


def test_render_sanitizes_experiment_names_with_spaces(tmp_path):
    proteins, fractions = _frames(["Left atrium", "Left ventricle"])
    written = render_linear_pdf(proteins, fractions, tmp_path)
    names = sorted(p.name for p in written)
    assert names == ["riana_rollup_curves_Left_atrium.pdf",
                     "riana_rollup_curves_Left_ventricle.pdf"]


def test_render_requires_delta_k_column(tmp_path):
    proteins, fractions = _frames(["exp1"])
    with pytest.raises(ValueError, match="not 'linear'"):
        render_linear_pdf(proteins.drop(columns=["delta_k"]), fractions, tmp_path)


@pytest.mark.parametrize("basis, symbol", [("theta", "θ"), ("f", "f"), ("mixed", "θ or f")])
def test_phi_axis_is_labelled_by_the_fit_basis(basis, symbol):
    """On a rising-precursor fit (fs = fraction new f) the φ axis says f, not θ."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from riana.core.linear_model import phi_axis_label
    from riana.core.rollup_plot import _plot_one

    assert phi_axis_label(basis) == f"φ = log(1 − {symbol})"
    fig, ax = plt.subplots()
    _plot_one(ax, "P1", {"control": ([1, 3, 7], [0.1, 0.3, 0.5], 0.1, 0.09, 0.11)},
              phi_limit=-4.0, fs_basis=basis)
    assert ax.get_ylabel() == f"Clearance  φ = log(1 − {symbol})"
    plt.close(fig)


def test_render_passes_the_fit_basis(tmp_path, monkeypatch):
    from riana.core import rollup_plot

    seen = []
    real = rollup_plot._plot_one
    monkeypatch.setattr(rollup_plot, "_plot_one",
                        lambda *a, **kw: (seen.append(kw["fs_basis"]), real(*a, **kw)))
    proteins, fractions = _frames(["exp1"])
    render_linear_pdf(proteins, fractions, tmp_path, fs_basis="f")
    assert seen == ["f"]
