"""Loading deconvoluted mass lists (Thermo BioPharma Finder exports, Excel or CSV)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

THERMO_COLUMNS = {
    "Monoisotopic Mass": "Mass",
    "Apex RT": "RT",
    "Sum Intensity": "Vol",
    "Relative Abundance": "RA",
    "Fractional Abundance": "FA",
}
FULL_COLUMNS = ["Mass", "RT", "Vol", "RA", "FA"]
CORE_COLUMNS = ["Mass", "RT", "Vol"]


def normalize_columns(df: pd.DataFrame, key_columns_only: bool = True) -> pd.DataFrame:
    """Rename Thermo export columns to Mass/RT/Vol/RA/FA.

    If all five columns exist, rows with missing values are dropped. If RA/FA are
    absent, only Mass/RT/Vol are kept (no rows dropped) - same as the original.
    The row index is preserved: manual exclusion lists refer to it.
    """
    df = df.rename(columns=THERMO_COLUMNS)
    if not key_columns_only:
        return df
    if set(FULL_COLUMNS).issubset(df.columns):
        df = df[FULL_COLUMNS].dropna()
    else:
        missing = set(CORE_COLUMNS) - set(df.columns)
        if missing:
            raise ValueError(f"Input is missing required columns: {sorted(missing)}")
        df = df[CORE_COLUMNS]
    return df.astype("float64")


def load_sample(path: str | Path, key_columns_only: bool = True) -> pd.DataFrame:
    """Load an .xlsx/.xls or .csv mass list into a Mass/RT/Vol(/RA/FA) DataFrame."""
    path = Path(path)
    if path.suffix.lower() == ".csv":
        raw = pd.read_csv(path)
    elif path.suffix.lower() in {".xlsx", ".xls"}:
        raw = pd.read_excel(path)
    else:
        raise ValueError(f"Unsupported file type: {path.suffix} (use .xlsx or .csv)")
    return normalize_columns(raw, key_columns_only)


def load_many(paths) -> list[pd.DataFrame]:
    """Load several files, sorted by file name (as the original notebooks do)."""
    return [load_sample(p) for p in sorted(map(str, paths))]
