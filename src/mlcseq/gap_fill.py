"""GapFill: recover ladder fragments that MassSum missed.

Walks along a ladder skeleton (sorted by mass). For each gap between two kept
fragments whose mass difference can be explained by some nucleotide composition,
it looks inside a box in (mass, retention time) space between the two ends and
keeps compounds that are themselves explainable from both ends.

Logic, thresholds and tie-breaks are ported 1:1 from the original notebook.
"""
from __future__ import annotations

import logging
from typing import NamedTuple

import numpy as np
import pandas as pd

from . import constants as K
from .composition import has_composition, soft_connected

log = logging.getLogger(__name__)

MIN_STEP = 600            # skeleton dots closer than this to the last kept dot are skipped
GAP_LIMIT_MAJOR = 3200
GAP_LIMIT_DEFAULT = 10000


class End(NamedTuple):
    Mass: float
    RT: float
    Vol: float


def _end(row) -> End:
    return End(float(row["Mass"]), float(row["RT"]), float(row["Vol"]))


def gap_box(left: End, right: End) -> tuple[float, float, float, float]:
    """Search window (mass_lo, mass_hi, rt_lo, rt_hi) between two ladder ends."""
    slope, icpt = np.polyfit([left.Mass, right.Mass], [left.RT, right.RT], 1)
    lo, hi = left.Mass + 200, right.Mass - 200
    rt_lo, rt_hi = slope * (lo - 150) + icpt, slope * (hi + 150) + icpt
    if lo > 10000:
        rt_lo, rt_hi = rt_lo - 0.2, rt_hi + 0.2
    return lo, hi, rt_lo, rt_hi


def standalone_dots(df: pd.DataFrame, left: End, right: End) -> pd.DataFrame:
    """Compounds inside the gap box that are explainable from BOTH ends."""
    lo, hi, rt_lo, rt_hi = gap_box(left, right)
    box = df[(df.Mass > lo) & (df.Mass < hi) & (df.RT > rt_lo) & (df.RT < rt_hi)]
    keep = [has_composition(abs(m - left.Mass)) and has_composition(abs(m - right.Mass))
            for m in box.Mass]
    found = box[np.array(keep, dtype=bool)] if len(box) else box
    log.debug("gap %.3f-%.3f: %d in box, %d explainable", left.Mass, right.Mass, len(box), len(found))
    return found.drop_duplicates()


def filter_dots_in_mid_gap(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only the best-connected compounds inside one gap.

    Builds a graph (node = mass, edge = explainable mass difference), then
    repeatedly removes the least-connected node (ties: lowest intensity, then
    first seen) until every remaining node has the same number of neighbours.
    """
    if df.empty or len(df) == 1:
        return df
    masses = df.Mass.to_numpy(dtype=float)
    vols = df.Vol.to_numpy(dtype=float)

    nbrs: dict[float, set] = {}
    vol: dict[float, float] = {}
    for i in range(len(masses)):
        for j in range(len(masses)):
            if soft_connected(masses[i] - masses[j]):
                a, b = masses[i], masses[j]
                nbrs.setdefault(a, set()).add(b)
                nbrs.setdefault(b, set()).add(a)
                vol[a], vol[b] = vols[i], vols[j]
    if not nbrs:
        return df.iloc[0:0]

    while nbrs:
        ranked = sorted(nbrs, key=lambda n: (len(nbrs[n]), vol[n]))  # stable sort
        if len(nbrs[ranked[0]]) >= len(nbrs[ranked[-1]]):
            break
        drop = ranked[0]
        for other in nbrs.pop(drop):
            nbrs[other].discard(drop)
    return df[df.Mass.isin(list(nbrs))]


def gap_fill(df: pd.DataFrame, skeleton: pd.DataFrame, fullmass_dot: pd.DataFrame,
             orientation: int = 5, major: bool = False,
             start_at: float = 0, stop_at: float = 0) -> pd.DataFrame:
    """Return compounds from `df` that fill gaps in the ladder `skeleton`.

    df            : all compounds of this ladder zone (index labels are kept)
    skeleton      : MassSum result for this ladder
    fullmass_dot  : row(s) of `df` at the intact mass; first row is used as the end
    orientation   : 5 or 3
    major         : if True only gaps < 3200 Da are explored (faster); else < 10000 Da
    """
    if orientation not in (3, 5):
        raise ValueError("orientation must be 3 or 5")
    limit = GAP_LIMIT_MAJOR if major else GAP_LIMIT_DEFAULT
    skeleton = skeleton.sort_values("Mass")
    last_label = skeleton.index[-1] if len(skeleton) else None
    start_tag = K.TAG_5P + K.TAG_3OH if orientation == 5 else K.TAG_3OH - K.TAG_5P

    kept: list = []
    gap_labels: list = []
    for label, row in skeleton.iterrows():
        mass = float(row.Mass)
        if start_at > 0 and mass < start_at:
            continue
        if stop_at > 0 and mass > stop_at:
            break
        left = End(start_tag, 0.001, 0.0) if not kept else _end(skeleton.loc[kept[-1]])
        pre = left.Mass
        is_last = label == last_label

        if not is_last and mass - pre < MIN_STEP:
            continue
        if mass - pre > limit:
            kept.append(label)
            continue

        if has_composition(mass - pre):
            kept.append(label)
            found = filter_dots_in_mid_gap(standalone_dots(df, left, _end(row)))
            gap_labels.extend(found.index)

        if is_last and mass - pre < limit and not fullmass_dot.empty and kept:
            left = _end(skeleton.loc[kept[-1]])
            full = _end(fullmass_dot.iloc[0])
            if stop_at > 0 and full.Mass > stop_at:
                break
            shift = K.TAG_5P if orientation == 5 else -K.TAG_5P
            right = full._replace(Mass=full.Mass + shift)
            gap_labels.extend(standalone_dots(df, left, right).index)
            gap_labels.append(fullmass_dot.index[0])

    gaps = df.loc[gap_labels]
    return gaps[~gaps.index.isin(skeleton.index)]
