"""MassSum: a 5' fragment + its complementary 3' fragment = intact mass + H2O."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import constants as K


def mass_sum(df3p: pd.DataFrame, df5p: pd.DataFrame, full_mass: float, error: float = 0.1,
             tag: float | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the (3' rows, 5' rows) whose masses add up to the intact mass.

    Pair sums and the target are rounded to 0.1 Da before comparing with +/- `error`
    (original behaviour). If the sample carries a 5' tag, pass its mass as `tag`.
    """
    m3 = df3p.Mass.to_numpy(dtype=float)
    m5 = df5p.Mass.to_numpy(dtype=float)
    sums = np.round(m3[:, np.newaxis] + m5, 1)
    if tag is not None:
        target = round(full_mass + tag + K.MASS_SUM_P - K.H2O, 1)
    else:
        target = round(full_mass + K.H2O, 1)
    i3, i5 = np.nonzero((sums >= target - error) & (sums <= target + error))
    return df3p.iloc[np.unique(i3)], df5p.iloc[np.unique(i5)]
