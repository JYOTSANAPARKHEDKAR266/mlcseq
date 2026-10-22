"""Which nucleotide compositions can explain a mass gap?

Replaces `components()` from the original gap_fill notebook. The original rebuilt a
brute-force `itertools.product` table on EVERY call (millions of tuples), which is
why GapFill took minutes. Here the table for each gap size is built once with numpy,
sorted, cached, and queried with a binary search. The search space and tolerances
are identical to the original, so the answers are identical.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from . import constants as K

BLOCK_NAMES = list(K.GAPFILL_BLOCKS)            # A C G U D y Methyl
BLOCK_MASSES = np.array(list(K.GAPFILL_BLOCKS.values()))
METHYL = K.GAPFILL_BLOCKS["Methyl"]
MAX_GAP = 10000.0
TOLERANCE = 0.1


def _compositions(total: int, parts: int, cap: int) -> np.ndarray:
    """All non-negative integer vectors of length `parts`, each <= cap, summing to total."""
    rows = np.zeros((1, 0), dtype=np.int16)
    remaining = np.array([total])
    for _ in range(parts - 1):
        new_rows, new_rem = [], []
        for v in range(min(cap, total) + 1):
            ok = remaining >= v
            new_rows.append(np.hstack([rows[ok], np.full((ok.sum(), 1), v, np.int16)]))
            new_rem.append(remaining[ok] - v)
        rows, remaining = np.vstack(new_rows), np.concatenate(new_rem)
    ok = remaining <= cap
    return np.hstack([rows[ok], remaining[ok, None].astype(np.int16)])


@lru_cache(maxsize=None)
def _table(size: int) -> tuple[np.ndarray, np.ndarray]:
    """(counts[n, 7], masses[n]) for a gap of roughly `size` nucleotides, sorted by mass.

    Search space (same as original): `size` units split over A, C, G, U, D, y
    (each at most size//2 once size >= 20), plus 0..max_methyl-1 methyl groups,
    where a composition may carry i methyls only if it has >= i-1 A/C/G/U.
    """
    cap = size if size < 20 else size // 2
    base = _compositions(size, 6, cap)
    max_methyl = 3 if size < 3 else size
    acgu = base[:, :4].sum(axis=1)
    blocks = [np.hstack([base[acgu >= i - 1], np.full(((acgu >= i - 1).sum(), 1), i, np.int16)])
              for i in range(max_methyl)]
    counts = np.vstack(blocks)
    masses = counts.astype(float) @ BLOCK_MASSES
    order = np.argsort(masses, kind="stable")
    return counts[order], masses[order]


def _candidates(mass: float) -> tuple[np.ndarray, np.ndarray]:
    mass = abs(mass)
    if mass > MAX_GAP:
        return np.empty((0, 7), np.int16), np.empty(0)
    size = round(mass / K.POSITION_MASS)
    if 200 < mass < 320:
        size = 1
    elif size < 1:
        # very small gap: only whole methyl groups are considered
        n = int(round(mass / METHYL))
        counts = np.zeros((1, 7), np.int16)
        counts[0, -1] = n
        return counts, np.array([n * METHYL])
    counts, masses = _table(size)
    hi = np.searchsorted(masses, mass + 1, side="right")   # original kept Mass <= mass + 1
    return counts[:hi], masses[:hi]


def _hits(mass: float) -> tuple[np.ndarray, np.ndarray]:
    counts, masses = _candidates(mass)
    diff = abs(mass) - masses
    ok = (diff > -TOLERANCE) & (diff < TOLERANCE)
    return counts[ok], masses[ok]


def has_composition(mass: float) -> bool:
    """True if some nucleotide/methyl composition explains `mass` within 0.1 Da."""
    return _hits(mass)[1].size > 0


def gap_compositions(mass: float) -> pd.DataFrame:
    """All compositions explaining `mass` (one row each), for inspection."""
    counts, masses = _hits(mass)
    out = pd.DataFrame(counts, columns=BLOCK_NAMES)
    out["Mass"] = masses
    out["MassDiff"] = abs(mass) - masses
    return out


def soft_connected(mass_diff: float) -> bool:
    """Original helper: a positive mass gap that some composition explains."""
    return mass_diff > 0 and has_composition(mass_diff)
