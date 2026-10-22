"""Ladder complementation: place every ladder's fragments and base calls on a
shared position axis so several ladders / isoforms can be read side by side."""
from __future__ import annotations

import logging
from typing import NamedTuple

import pandas as pd

from . import constants as K
from .basecalling import base_calling

log = logging.getLogger(__name__)


class Ladder(NamedTuple):
    dataframe: pd.DataFrame
    fullmass: float
    ori: int  # 5 or 3


def _position(mass: float, seats: int, orientation: int, for_base: bool) -> int:
    pos = int(mass // K.POSITION_MASS)
    if orientation == 3:
        # original code: fragments map to seats+1-pos, base calls to seats-pos
        pos = seats - pos if for_base else seats + 1 - pos
    return pos


def place_ladder(calls_df: pd.DataFrame, pairs, full_mass: float, orientation: int,
                 positions: pd.Index) -> tuple[dict, dict]:
    """Return {position: [masses]} and {position: [bases]} for one ladder."""
    seats = int(full_mass // K.POSITION_MASS)
    masses: dict[int, list[float]] = {}
    bases: dict[int, list[str]] = {}

    for m in calls_df.Mass.fillna(0):
        if m <= 0:
            continue
        pos = _position(m, seats, orientation, for_base=False)
        if pos not in positions:
            log.warning("fragment %.4f maps to position %d outside 1..%d, skipped", m, pos, positions[-1])
            continue
        masses.setdefault(pos, []).append(float(m))

    for heavy, light, name in pairs:
        pair = calls_df[calls_df.Mass.isin([heavy, light])]
        if pair.empty:
            continue
        anchor = pair.Mass.max() if orientation == 5 else pair.Mass.min()
        pos = _position(anchor, seats, orientation, for_base=True)
        if pos not in positions:
            log.warning("base call %s at %.4f maps outside the tRNA, skipped", name, anchor)
            continue
        bases.setdefault(pos, []).append(name)
    return masses, bases


def ladder_complementation(ladders: list[Ladder], ppm: float = 10) -> pd.DataFrame:
    """One row per position; columns Base{i} / Mass{i} for each ladder i.

    Cells hold comma-separated values when several fragments/bases land on one
    position. Results still need expert review (as in the original).
    """
    if not ladders:
        raise ValueError("no ladders given")
    seats = int(ladders[0].fullmass // K.POSITION_MASS)
    positions = pd.RangeIndex(1, seats + 2, name="position")
    out = pd.DataFrame(index=positions)
    for i, ladder in enumerate(ladders):
        calls = base_calling(ladder.dataframe, ppm=ppm)
        masses, bases = place_ladder(calls.compounds, calls.pairs, ladder.fullmass,
                                     ladder.ori, positions)
        out[f"Base{i}"] = [", ".join(bases.get(p, [])) for p in positions]
        out[f"Mass{i}"] = [", ".join(map(str, masses.get(p, []))) for p in positions]
    return out
