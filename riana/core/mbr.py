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

    PSMs are grouped by ``(identity.group_key, identity.fraction)`` — one turnover
    curve, one fraction (transfers never cross fractions, which carry different
    peptides). Records without an SDRF identity (the demoted Percolator path) are
    passed through untouched: MBR needs the run-identity grouping. The input order
    is preserved and the synthetic records are appended.
    """
    groups: dict[tuple, list[PSMRecord]] = defaultdict(list)
    n_no_identity = 0
    for p in psms:
        if p.identity is None:
            n_no_identity += 1
            continue
        groups[(p.identity.group_key, p.identity.fraction)].append(p)

    transfers: list[PSMRecord] = []
    for group_psms in groups.values():
        transfers.extend(_transfer_within_group(group_psms, config))

    if transfers:
        _LOGGER.info(
            "MBR: %d transfer records added across %d curve group(s) "
            "(donor q<=%.3g in >=%d runs).",
            len(transfers), len(groups), config.mbr_donor_q,
            config.mbr_min_donor_runs,
        )
    elif n_no_identity == 0:
        _LOGGER.info("MBR: no eligible transfers found.")
    return list(psms) + transfers


def _transfer_within_group(
    group_psms: Sequence[PSMRecord], config: IntegrationConfig
) -> list[PSMRecord]:
    """Emit MBR records for one ``(group_key, fraction)`` curve group."""
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
