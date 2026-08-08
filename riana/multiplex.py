# -*- coding: utf-8 -*-

"""Multiplexing (sample-axis) metabolic labels — the label registry (1.2.0).

A *multiplexing* label marks a **different sample** that is co-injected but stays
**separable at MS1** because the label shifts the whole precursor mass: dimethyl
(light UNIMOD:36 / medium :2H(4) 199 / heavy :2H(6)13C(2) 330) and SILAC (heavy
Lys/Arg). This is the opposite of isobaric TMT/iTRAQ, whose channels share one MS1
cluster and are *merged* into an average (``io.sdrf._collapse_isobaric_runs``,
``constants.CHEMICAL_MODS``). Here the channels are **expanded** into distinct
samples: each is fit separately at its own precursor enrichment and combined at
rollup as separate conditions (a two-condition Δk).

This module is the single source of truth for *which* labels multiplex and how their
channels map to (a) the SDRF ``comment[label]`` CV term and (b) the UNIMOD mod a
peptidoform carries. The ``io.sdrf`` / ``io.mztab`` channel intake routes each PSM to
its channel by :func:`channel_of`; the SDRF reader recognises a multiplexed sheet by
:func:`is_multiplex_cv` and keys each channel's identity by :func:`label_channel_for_cv`.

The per-channel ``shift_per_site`` geometry is declared here too (used by the
forthcoming forward-model spillover gate — the sibling-cluster offset is
``Σ shift_per_site[site]``): dimethyl is **uniform** (+8.0444 Da on the N-term and each
K for the heavy channel); SILAC is **heterogeneous** (K vs R differ), which the model
already accommodates. The heavy UNIMOD masses (``constants.mod_fixed_isotopes``) are
the pinned-isotope pseudo-elements TMT introduced.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Channel:
    """One channel of a multiplexing label.

    ``mod`` is the UNIMOD id a peptidoform in this channel carries (``None`` for a
    no-mod light channel, e.g. SILAC light). ``cv_term`` is the ``NT=`` value of the
    SDRF ``comment[label]`` cell (``None`` when intake for this label is not yet
    wired). ``shift_per_site`` is the heavy mass shift (Da) **relative to the light
    channel**, keyed by site type (``"n_term"`` or a residue letter) — empty for light.
    """

    name: str
    mod: int | None
    cv_term: str | None = None
    shift_per_site: Mapping[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class MultiplexLabel:
    """A multiplexing label: its channels + the site rule the geometry needs.

    ``labels_n_term`` / ``site_residues`` are which positions carry the label (with
    each channel's ``shift_per_site`` they give the sibling-cluster offset).
    ``intake_wired`` is ``False`` for a label registered as **geometry only** (no CV /
    mod intake yet — e.g. SILAC, pending data + heavy-UNIMOD constants).
    """

    name: str
    channels: tuple[Channel, ...]
    labels_n_term: bool
    site_residues: tuple[str, ...]
    intake_wired: bool = True


#: Dimethyl duplex/triplex. Fully wired for intake. Heavy channel sits +8.0444 Da/site
#: above light (medium +4.0251); values equal the ``unimod_mass`` deltas (330−36, 199−36).
DIMETHYL = MultiplexLabel(
    name="dimethyl",
    labels_n_term=True,
    site_residues=("K",),
    channels=(
        Channel("light", 36, "DIMETHYL0"),
        Channel("medium", 199, "DIMETHYL4", {"n_term": 4.0251, "K": 4.0251}),
        Channel("heavy", 330, "DIMETHYL8", {"n_term": 8.0444, "K": 8.0444}),
    ),
)

#: SILAC (K+6 / R+10). **Geometry only** — no CV terms / heavy-UNIMOD constants yet
#: (they land with SILAC-D₂O data). Registered now to prove the heterogeneous-shift
#: model: Lys6 (+6.0201 on K) and Arg10 (+10.0083 on R) differ, so the sibling offset
#: is per-residue, not uniform.
SILAC = MultiplexLabel(
    name="silac",
    labels_n_term=False,
    site_residues=("K", "R"),
    intake_wired=False,
    channels=(
        Channel("light", None),
        Channel("heavy", None, None, {"K": 6.0201, "R": 10.0083}),
    ),
)

#: All registered multiplexing labels.
REGISTRY: tuple[MultiplexLabel, ...] = (DIMETHYL, SILAC)

# Derived lookups (built once). Only intake-wired channels contribute a CV/mod key.
_CV_TO_CHANNEL: dict[str, tuple[str, str]] = {}
_MOD_TO_CHANNEL: dict[int, tuple[str, str]] = {}
for _lbl in REGISTRY:
    if not _lbl.intake_wired:
        continue  # geometry-only label (e.g. SILAC): contributes no CV/mod intake key
    for _ch in _lbl.channels:
        if _ch.cv_term is not None:
            _CV_TO_CHANNEL[_ch.cv_term.upper()] = (_lbl.name, _ch.name)
        if _ch.mod is not None:
            _MOD_TO_CHANNEL[_ch.mod] = (_lbl.name, _ch.name)


def _cv_name(cell: str) -> str | None:
    """The ``NT=`` value of an SDRF ``comment[label]`` cell, e.g.
    ``"NT=DIMETHYL0;AC=PRIDE:0000848"`` → ``"DIMETHYL0"``. ``None`` if absent."""
    if not cell:
        return None
    for part in cell.split(";"):
        key, _, val = part.partition("=")
        if key.strip().upper() == "NT":
            return val.strip()
    return None


def is_multiplex_cv(cell: str) -> bool:
    """Does this SDRF ``comment[label]`` cell name a registered multiplexing channel?

    Used by :func:`io.sdrf.read_sdrf` to detect a sample-axis-multiplexed sheet (the
    sibling of the ``tmt|itraq`` isobaric check) and switch to per-channel intake.
    """
    name = _cv_name(cell)
    return name is not None and name.upper() in _CV_TO_CHANNEL


def label_channel_for_cv(cell: str) -> tuple[str, str] | None:
    """``(label_name, channel_name)`` for an SDRF ``comment[label]`` cell, or ``None``.

    e.g. ``"NT=DIMETHYL8;AC=PRIDE:0000852"`` → ``("dimethyl", "heavy")``.
    """
    name = _cv_name(cell)
    return _CV_TO_CHANNEL.get(name.upper()) if name else None


def channel_of(mods: Iterable[int]) -> tuple[str, str] | None:
    """``(label_name, channel_name)`` for a peptidoform, from its UNIMOD mod ids.

    Non-multiplexing mods (CAM, Met-Ox, phospho, …) are ignored. Returns ``None`` when
    the peptidoform carries no multiplexing mod (route to the default/unlabelled path).
    Raises :class:`ValueError` if it carries mods from **two different** channels (a
    mixed-channel ID — the caller should drop that PSM), which cannot happen for a
    correctly labelled peptide (all its dimethyls are one channel).
    """
    found = {_MOD_TO_CHANNEL[m] for m in mods if m in _MOD_TO_CHANNEL}
    if not found:
        return None
    if len(found) > 1:
        raise ValueError(
            f"peptidoform carries conflicting multiplexing channels {sorted(found)} "
            "— a PSM cannot belong to two channels at once"
        )
    return next(iter(found))
