"""Mass matching helpers (vectorised replacements for match_dfs / match_dfs_v2)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def drop_dups(df: pd.DataFrame, digit: int = 0) -> pd.DataFrame:
    """Keep the most intense compound per rounded mass, then sort by mass."""
    out = df.copy()
    out["RoundedMass"] = np.round(out.Mass, digit)
    out = out.sort_values("Vol", ascending=False).drop_duplicates("RoundedMass")
    return out.sort_values("Mass")


def find_candidates(query: np.ndarray, reference: np.ndarray, ppm: float) -> list[np.ndarray]:
    """For every query mass m, return positions in `reference` (ascending, i.e. in
    reference order) whose mass d satisfies |d - m| < 1 Da and |d - m| / m * 1e6 < ppm.
    """
    query = np.asarray(query, dtype=float)
    reference = np.asarray(reference, dtype=float)
    order = np.argsort(reference, kind="stable")
    ref_sorted = reference[order]
    out = []
    for m in query:
        lo = np.searchsorted(ref_sorted, m - 1, side="right")
        hi = np.searchsorted(ref_sorted, m + 1, side="left")
        pos = order[lo:hi]
        if pos.size:
            ppm_err = np.abs(1e6 * (reference[pos] - m)) / m
            pos = np.sort(pos[ppm_err < ppm])
        out.append(pos)
    return out


def match_masses(source: pd.DataFrame, target: pd.DataFrame, ppm: float = 10,
                 shift: float = 0.0, copy_cols=()) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Find compounds in `target` whose mass matches any (shifted) mass in `source`.

    Replaces the original `match_dfs`. Returns:
      matched_target - rows of `target` (index reset to 0..n-1) that matched
      source_annot   - copy of `source` with Match / MatchedMass / PPM columns
    """
    target = target.reset_index(drop=True)
    src = source.copy()
    masses = src.Mass.to_numpy(dtype=float) + shift
    cands = find_candidates(masses, target.Mass.to_numpy(), ppm)

    hit = np.array([c.size > 0 for c in cands], dtype=bool)
    first = np.array([c[0] if c.size else -1 for c in cands])
    src["Mass"] = masses
    src["Match"] = hit
    if hit.any():
        matched = target.Mass.to_numpy()[first[hit]]
        _set_where(src, "MatchedMass", hit, matched)
        for col in copy_cols:
            _set_where(src, col, hit, target[col].to_numpy()[first[hit]])
        _set_where(src, "PPM", hit, np.abs(masses[hit] - matched) * 1e6 / matched)

    all_pos = np.unique(np.concatenate(cands)) if cands else np.array([], dtype=int)
    matched_target = target.iloc[all_pos.astype(int)].copy()
    return matched_target, src


def annotate_matches(target: pd.DataFrame, reference: pd.DataFrame, ppm: float = 30,
                     shift: float = 0.0, copy_cols=()) -> pd.DataFrame:
    """Annotate each row of `target` with its first matching (shifted) reference mass.

    Replaces the side-effect part of the original `match_dfs_v2(..., inplace=True)`.
    Behaviour kept from the original: `Match` is reset on every call, but
    MatchedMass / Sft / PPM / copy_cols from EARLIER calls are left in place for rows
    that do not match this time. Returns a new DataFrame (no hidden mutation).
    """
    out = target.copy()
    ref = reference.reset_index(drop=True)
    ref_mass = ref.Mass.to_numpy(dtype=float) + shift
    cands = find_candidates(out.Mass.to_numpy(), ref_mass, ppm)

    hit = np.array([c.size > 0 for c in cands], dtype=bool)
    first = np.array([c[0] if c.size else -1 for c in cands])
    out["Match"] = hit
    if not hit.any():
        return out
    matched = ref_mass[first[hit]]
    _set_where(out, "MatchedMass", hit, matched)
    _set_where(out, "Sft", hit, np.full(hit.sum(), shift))
    for col in copy_cols:
        _set_where(out, col, hit, ref[col].to_numpy()[first[hit]])
    _set_where(out, "PPM", hit,
               np.abs(out.Mass.to_numpy()[hit] - matched) * 1e6 / matched)
    return out


def _set_where(df: pd.DataFrame, col: str, mask: np.ndarray, values) -> None:
    """Positional (not label-based) column update, so duplicate index labels are safe."""
    values = np.asarray(values)
    if col in df.columns:
        arr = df[col].to_numpy(copy=True)
    else:
        arr = np.full(len(df), np.nan, dtype=object if values.dtype == object else float)
    if arr.dtype != object and values.dtype == object:
        arr = arr.astype(object)
    arr[mask] = values
    df[col] = arr
