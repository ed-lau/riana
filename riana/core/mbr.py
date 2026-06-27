"""Match-between-runs (MBR) for the mzTab/DDA path.

DDA picks precursors for MS2 stochastically, so a peptide drops out of individual
runs and leaves holes in its turnover curve (measured: DDA curves 8–15% complete,
~40–45% recoverable, 70% LVE t0-anchor loss — ``reports/2026-06-17_mbr_dda_feasibility.md``).
This module fills those holes: for each ``(experiment, condition)`` curve group it
finds **donor** precursors (confidently identified in ≥ ``mbr_min_donor_runs`` runs)
and, for each run of the group that *lacks* the precursor, emits a synthetic
:class:`~riana.records.PSMRecord` carrying the donor's identity and an
RT-frame-aligned retention time, flagged ``evidence="mbr"`` with ``scan=-1``.

The synthetic records are appended to the PSM list at intake
(:func:`riana.core.pipeline.plan_integration`) and flow through the *existing*
RT-anchored extraction path (``resolve_rt_anchored_scans`` → apex re-detect), the
same one DIA uses. A transferred precursor that has no detectable apex in the
acceptor run is dropped at integration time (``core.integration.integrate_run``),
never integrated as baseline.

**RT alignment.** quantms RT is aligned across runs but with run-specific residuals
up to ~15–25 s (``reports/2026-06-17_mbr_dda_feasibility.md``) — larger than the
narrow integration window. We remove that with a robust per-run offset to a
consensus frame: ``consensus_rt[p]`` = median RT of precursor *p* over the runs
where it is located; ``δ[run]`` = median over *that run's* located precursors of
``rt - consensus_rt``; the transfer RT for (p, run) is ``consensus_rt[p] + δ[run]``.
Offsets (not a global linear warp) match the measured near-constant shifts, and the
medians are robust to the heavy precursor-RT tail. The anchor only needs to land
inside the apex-search window — the apex finder does the final placement.

Design + validation plan: ``reports/2026-06-17_mbr_v1_design.md``.
"""

from __future__ import annotations

import dataclasses
import logging
from collections import defaultdict
from collections.abc import Sequence

import numpy as np

from riana.config import IntegrationConfig
from riana.records import PSMRecord

_LOGGER = logging.getLogger(__name__)


def augment_with_mbr(
    psms: Sequence[PSMRecord], config: IntegrationConfig
) -> list[PSMRecord]:
    """Return ``psms`` plus match-between-runs transfer records.

    Records are grouped into curves by ``identity.group_key`` (``(experiment,
    condition)``). Within a curve, each precursor is matched in a **single winner
    fraction** — the LC fraction in which it has the most identifications (ties
    broken by best q-value, then fraction number — :func:`_winner_fractions`).
    Transfers are emitted only among *that* fraction's runs; the precursor's
    incidental appearances in other fractions are left alone. This is the
    conservative minimal policy: it fills holes where a peptide reliably elutes and
    deliberately does **not** model cross-fraction drift (no RT-correlation matching
    across fractions). For single-fraction data every winner is the lone fraction,
    so this is a structural no-op.

    Records without an SDRF identity (the demoted Percolator path) pass through
    untouched: MBR needs the run-identity grouping. The input order is preserved and
    the synthetic records are appended.
    """
    by_group: dict[tuple, list[PSMRecord]] = defaultdict(list)
    n_no_identity = 0
    for p in psms:
        if p.identity is None:
            n_no_identity += 1
            continue
        by_group[p.identity.group_key].append(p)

    transfers: list[PSMRecord] = []
    for group_psms in by_group.values():
        winner = _winner_fractions(group_psms, config.q_value)
        by_fraction: dict[int, list[PSMRecord]] = defaultdict(list)
        for p in group_psms:
            by_fraction[p.identity.fraction].append(p)
        for frac, frac_psms in by_fraction.items():
            eligible = {c for c, wf in winner.items() if wf == frac}
            if not eligible:
                continue
            transfers.extend(
                _transfer_within_group(frac_psms, config, eligible))

    if transfers:
        _LOGGER.info(
            "MBR: %d transfer records added across %d curve group(s) "
            "(donor q<=%.3g in >=%d runs; per-precursor winner fraction).",
            len(transfers), len(by_group), config.mbr_donor_q,
            config.mbr_min_donor_runs,
        )
    elif n_no_identity == 0:
        _LOGGER.info("MBR: no eligible transfers found.")
    return list(psms) + transfers


def _winner_fractions(
    group_psms: Sequence[PSMRecord], q_value: float
) -> dict[str, int]:
    """Map each precursor to its **winner fraction** within one curve group.

    The winner is the fraction with the most distinct runs carrying a located
    (``q <= q_value``) ID — where the peptide most reliably elutes — with ties
    broken by best (lowest) q-value, then lowest fraction number (deterministic).
    Precursors with no located ID are absent (they cannot seed MBR).
    """
    runs_by: dict[str, dict[int, set[int]]] = defaultdict(lambda: defaultdict(set))
    best_q: dict[str, dict[int, float]] = defaultdict(dict)
    for p in group_psms:
        if p.percolator_q_value > q_value:
            continue
        frac = p.identity.fraction
        runs_by[p.concat][frac].add(p.file_idx)
        prev = best_q[p.concat].get(frac)
        if prev is None or p.percolator_q_value < prev:
            best_q[p.concat][frac] = p.percolator_q_value
    winner: dict[str, int] = {}
    for concat, by_frac in runs_by.items():
        winner[concat] = min(
            by_frac,
            key=lambda f: (-len(by_frac[f]), best_q[concat][f], f),
        )
    return winner


def _transfer_within_group(
    group_psms: Sequence[PSMRecord],
    config: IntegrationConfig,
    eligible_concats: set[str],
) -> list[PSMRecord]:
    """Emit MBR records for one ``(group_key, fraction)`` curve group, restricted to
    the precursors whose winner fraction is this one (``eligible_concats``)."""
    # "Located" = confidently extracted (q <= q_value): defines where a precursor
    # already has a real peak (so it is NOT an acceptor there) and supplies the
    # RT anchors. The donor gate is the stricter-or-equal mbr_donor_q.
    located = [p for p in group_psms if p.percolator_q_value <= config.q_value]
    if not located:
        return []

    run_rep: dict[int, PSMRecord] = {}          # run (file_idx) → a record, for acceptor identity
    present: dict[str, set[int]] = defaultdict(set)   # concat → runs with a located ID
    donor_runs: dict[str, set[int]] = defaultdict(set)  # concat → runs at q<=mbr_donor_q
    donor_tpl: dict[str, PSMRecord] = {}         # concat → best-q located record (template)
    rt_lists: dict[int, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))

    for p in located:
        run_rep.setdefault(p.file_idx, p)
        present[p.concat].add(p.file_idx)
        rt_lists[p.file_idx][p.concat].append(p.retention_time)
        if (p.concat not in donor_tpl
                or p.percolator_q_value < donor_tpl[p.concat].percolator_q_value):
            donor_tpl[p.concat] = p
        if p.percolator_q_value <= config.mbr_donor_q:
            donor_runs[p.concat].add(p.file_idx)

    all_runs = set(run_rep)
    if len(all_runs) < config.mbr_min_donor_runs:
        return []  # not enough runs in the group to corroborate a donor

    # Median RT per (run, precursor) → consensus RT per precursor → per-run offset.
    rt_rc = {
        fidx: {c: float(np.median(v)) for c, v in precs.items()}
        for fidx, precs in rt_lists.items()
    }
    rt_by_concat: dict[str, list[float]] = defaultdict(list)
    for precs in rt_rc.values():
        for concat, rt in precs.items():
            rt_by_concat[concat].append(rt)
    consensus_rt = {c: float(np.median(v)) for c, v in rt_by_concat.items()}
    run_offset = {
        fidx: float(np.median([rt - consensus_rt[c] for c, rt in precs.items()]))
        for fidx, precs in rt_rc.items()
        if precs
    }

    out: list[PSMRecord] = []
    for concat, runs_seen in donor_runs.items():
        if concat not in eligible_concats:
            continue  # this precursor's winner fraction is a different one
        if len(runs_seen) < config.mbr_min_donor_runs:
            continue
        acceptors = all_runs - present[concat]  # runs with no located ID for this precursor
        if not acceptors:
            continue
        tpl = donor_tpl[concat]
        for fidx in acceptors:
            rep = run_rep[fidx]
            out.append(
                dataclasses.replace(
                    tpl,
                    scan=-1,  # RT-anchored: resolved to the nearest MS1 scan at integrate time
                    retention_time=consensus_rt[concat] + run_offset.get(fidx, 0.0),
                    evidence="mbr",
                    pep_id=-1,
                    file_idx=fidx,
                    file_name=rep.file_name,
                    sample=rep.sample,
                    identity=rep.identity,
                    percolator_q_value=0.0,  # sentinel: transferred, not q-scored
                    percolator_score=0.0,
                    percolator_pep=1.0,
                    distinct_matches=0,
                )
            )
    return out
