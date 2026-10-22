"""Base calling from mass differences, and homology search on intact masses."""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pandas as pd

from . import constants as K


class BaseCalls(NamedTuple):
    compounds: pd.DataFrame          # compounds that take part in at least one call
    pairs: list[tuple[float, float, str]]  # (heavier mass, lighter mass, base name)


def _select(bases: dict, names) -> dict:
    if not names:
        return dict(bases)
    missing = [n for n in names if n not in bases]
    if missing:
        raise KeyError(f"Unknown base name(s): {missing}")
    return {n: bases[n] for n in names}


def _diff_pairs(heavy: np.ndarray, light: np.ndarray, delta: float, ppm: float):
    """Index pairs (i, j) with heavy[i] - light[j] == delta within `ppm` (row-major order)."""
    target = light[np.newaxis, :] + delta
    err = np.abs(heavy[:, np.newaxis] - target) * 1e6 / target
    return np.argwhere(err < ppm)


def base_calling(df: pd.DataFrame, bases: dict | None = None, names=(), ppm: float = 10,
                 acid_labile: bool = False) -> BaseCalls:
    """Find every pair of compounds whose mass difference equals a known building block.

    Replaces `base_calling_random`. `bases` defaults to the standard nucleotide table
    (constants.BASES_DEFAULT). With acid_labile=True the data must have a boolean
    `isAD` column and pairs are searched between non-degraded (heavy) and degraded rows.
    """
    bases = _select(K.BASES_ACID_LABILE if acid_labile and bases is None
                    else (bases or K.BASES_DEFAULT), names)
    df = df.sort_index()
    if acid_labile:
        heavy_df, light_df = df[~df.isAD.astype(bool)], df[df.isAD.astype(bool)]
    else:
        heavy_df = light_df = df
    heavy = heavy_df.Mass.to_numpy(dtype=float)
    light = light_df.Mass.to_numpy(dtype=float)

    pairs, positions = [], set()
    for name, delta in bases.items():
        for i, j in _diff_pairs(heavy, light, delta, ppm):
            pairs.append((heavy[i], light[j], name))
            positions.update((i, j))
    pairs = [p for p in pairs if p[0] != p[1]]

    if acid_labile:
        used = {m for p in pairs for m in p[:2]}
        compounds = df[df.Mass.isin(used)]
    else:
        compounds = df.iloc[sorted(positions)]
    return BaseCalls(compounds, pairs)


def homology_search(df: pd.DataFrame, names=(), ppm: float = 10,
                    min_mass: float = K.INTACT_MIN_MASS) -> BaseCalls:
    """Group intact tRNAs (mass > min_mass) that differ by C, A, +O, +Me or +2Me."""
    return base_calling(df[df.Mass > min_mass], bases=K.BASES_HOMOLOGY, names=names, ppm=ppm)
