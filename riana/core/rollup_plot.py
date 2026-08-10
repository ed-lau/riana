"""Headless PDF export of ``linear simple`` rollup curves — the batch/CLI analog of
the GUI's one-at-a-time φ-space view (``riana.gui.curve_view.CurveView.plot_linear``).

``riana rollup --model "linear simple" --plot`` renders every protein's control-vs-
condition Δk comparison to a multi-page PDF, one PDF per experiment (Riana groups on
``(experiment, condition)``, so each experiment is a self-contained set of contrasts).
matplotlib is a core dependency, so this adds nothing new; it runs on the Agg backend.

Reads the two rollup tables the GUI's ``load_rollup_results`` also reconstructs from:
``riana_rollup_proteins.txt`` (per-condition ``k_deg`` / ``ci_lo`` / ``ci_hi`` / ``delta_k``
/ ``delta_k_p_adj``) and ``riana_rollup_fractions.txt`` (per-point collapsed ``fs`` = θ).

Performance: per-experiment PDFs are embarrassingly parallel, so ``workers > 1``
dispatches one experiment per process (the ``-W`` lever); within a PDF the pages
render serially, reusing a single figure and pre-grouping the frames by protein
(both matter — the naive per-protein DataFrame filter + a fresh figure + tight_layout
per page is ~2× slower).
"""
from __future__ import annotations

import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from riana.core.linear_model import to_phi

# Match the GUI palette (gui/curve_view.py:_CONDITION_COLOURS) so CLI and GUI plots
# read the same colour-per-condition.
_CONDITION_COLOURS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]


def _plot_one(ax, title: str, per_condition: dict, *, phi_limit: float,
              delta_k=None, delta_k_p_adj=None) -> None:
    """Render one protein's φ-space condition comparison onto *ax* (mirrors the GUI)."""
    t_max = 0.0
    for i, cond in enumerate(sorted(per_condition)):
        t, theta, k, ci_lo, ci_hi = per_condition[cond]
        t = np.asarray(t, dtype=float)
        phi = to_phi(np.asarray(theta, dtype=float))
        colour = _CONDITION_COLOURS[i % len(_CONDITION_COLOURS)]
        kept = phi > phi_limit
        t_max = max(t_max, float(t.max()) if t.size else 0.0)
        ax.scatter(t[kept], phi[kept], s=32, facecolors=colour, edgecolors=colour,
                   label=str(cond), zorder=3)                      # kept (fit) points
        ax.scatter(t[~kept], phi[~kept], s=32, facecolors="none", edgecolors=colour,
                   zorder=3)                                        # truncated (hollow)
        if k is not None and np.isfinite(k):
            grid = np.linspace(0.0, t_max if t_max > 0 else 1.0, 100)
            ax.plot(grid, -k * grid, color=colour, lw=2, label=f"{cond}  k={k:.3g}")
            if ci_lo is not None and ci_hi is not None and np.isfinite(ci_lo) and np.isfinite(ci_hi):
                # higher k ⇒ more-negative φ, so ci_hi is the lower edge
                ax.fill_between(grid, -ci_hi * grid, -ci_lo * grid, color=colour,
                                alpha=0.15, lw=0)
    ax.axhline(phi_limit, ls="--", color="grey", lw=0.8, alpha=0.7)
    ax.set_xlabel("Time")
    ax.set_ylabel("Clearance  φ = log(1 − θ)")
    if delta_k is not None and np.isfinite(delta_k):
        title += f"   Δk={delta_k:+.3g}"
    if delta_k_p_adj is not None and np.isfinite(delta_k_p_adj):
        title += f"  (p_adj={delta_k_p_adj:.2g})"
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6, loc="upper right", framealpha=0.6)


def _points_from_groups(frac_g: dict, prot_g: dict, prot: str) -> dict:
    """``{condition: (t, θ, k, ci_lo, ci_hi)}`` for one protein, from pre-grouped
    per-protein slices (O(1) lookup instead of a full-frame filter per protein)."""
    fr = frac_g.get(prot)
    if fr is None:
        return {}
    pr = prot_g.get(prot)
    kci = {}
    if pr is not None:
        for _, r in pr.iterrows():
            kci[str(r["condition"])] = (float(r["k_deg"]), float(r["ci_lo"]), float(r["ci_hi"]))
    per_condition = {}
    for cond, sub in fr.groupby("condition"):
        sub = sub.sort_values("labeling_time")
        k, ci_lo, ci_hi = kci.get(str(cond), (float("nan"),) * 3)
        per_condition[str(cond)] = (sub["labeling_time"].to_numpy(),
                                    sub["fs"].to_numpy(), k, ci_lo, ci_hi)
    return per_condition


def _safe(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_") or "experiment"


def _render_experiment(args) -> tuple[str, int]:
    """Render one experiment's proteins to one PDF (a picklable process-pool task).

    *args* = ``(exp, proteins_exp, fractions_exp, out_path, phi_limit)`` — pre-sliced
    to this experiment so nothing else crosses the process boundary.
    """
    exp, proteins_exp, fractions_exp, out_path, phi_limit = args
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    pairs = (proteins_exp.drop_duplicates(["experiment", "protein"])
             .sort_values("delta_k_p_adj", na_position="last"))
    prots = list(pairs["protein"])
    frac_g = {p: g for p, g in fractions_exp.groupby("protein")}
    prot_g = {p: g for p, g in proteins_exp.groupby("protein")}

    n = 0
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    fig.subplots_adjust(left=0.12, right=0.97, top=0.92, bottom=0.12)  # fixed → skip tight_layout
    dk = dict(zip(pairs["protein"], pairs.get("delta_k", pd.Series(dtype=float))))
    padj = dict(zip(pairs["protein"], pairs.get("delta_k_p_adj", pd.Series(dtype=float))))
    with PdfPages(out_path) as pdf:
        for prot in prots:
            per_condition = _points_from_groups(frac_g, prot_g, prot)
            if not per_condition:
                continue
            ax.cla()
            _plot_one(ax, f"{exp} · {prot}", per_condition, phi_limit=phi_limit,
                      delta_k=dk.get(prot), delta_k_p_adj=padj.get(prot))
            pdf.savefig(fig)
            n += 1
    plt.close(fig)
    return str(out_path), n


def render_linear_pdf(proteins: pd.DataFrame, fractions: pd.DataFrame, out_dir: Path,
                      *, phi_limit: float = -4.0, workers: int = 1,
                      progress: "Callable[[str, int], None] | None" = None) -> list[Path]:
    """Write ``linear simple`` Δk curves to one PDF per experiment.

    One page per protein (sorted by ``delta_k_p_adj``, most significant first),
    filled kept points vs hollow plateau-truncated points at ``φ > phi_limit``, a
    through-origin −k·t line + CI ribbon per condition, and Δk / p_adj in the title.
    A single unnamed experiment writes ``riana_rollup_curves.pdf``; multiple write
    ``riana_rollup_curves_<experiment>.pdf``. ``workers > 1`` renders experiments in
    parallel (one process per experiment). Returns the paths written.
    """
    if "delta_k" not in proteins.columns:
        raise ValueError("rollup is not 'linear simple' (no delta_k) — nothing to plot")

    out_dir = Path(out_dir)
    experiments = [e for e in proteins["experiment"].dropna().unique()]
    multi = len(experiments) > 1
    tasks = []
    for exp in experiments:
        pexp = proteins[proteins["experiment"] == exp]
        fexp = fractions[fractions["experiment"] == exp]
        if pexp.drop_duplicates(["experiment", "protein"]).empty:
            continue
        fname = f"riana_rollup_curves_{_safe(exp)}.pdf" if multi else "riana_rollup_curves.pdf"
        tasks.append((exp, pexp, fexp, out_dir / fname, phi_limit))

    written: list[Path] = []
    n_par = min(int(workers), len(tasks))
    if n_par > 1:
        with ProcessPoolExecutor(max_workers=n_par) as ex:
            results = list(ex.map(_render_experiment, tasks))
    else:
        results = [_render_experiment(t) for t in tasks]
    for path, n in results:
        written.append(Path(path))
        if progress is not None:
            progress(path, n)
    return written
