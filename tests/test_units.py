"""Fast unit tests on small synthetic inputs (no data files needed)."""
from itertools import filterfalse, product

import numpy as np
import pandas as pd
import pytest

from mlcseq import constants as K
from mlcseq.basecalling import base_calling, homology_search
from mlcseq.complementation import Ladder, ladder_complementation
from mlcseq.composition import has_composition
from mlcseq.io import normalize_columns
from mlcseq.mass_sum import mass_sum
from mlcseq.matching import annotate_matches, drop_dups, match_masses


# ---------------------------------------------------------------- helpers
def ladder_from_seq(seq: str, start: float = 0.0) -> np.ndarray:
    m = {"A": K.A, "C": K.C, "G": K.G, "U": K.U}
    return start + np.cumsum([m[b] for b in seq])


def frame(masses, rt=None, vol=None):
    n = len(masses)
    return pd.DataFrame({"Mass": masses,
                         "RT": rt if rt is not None else np.linspace(1, 10, n),
                         "Vol": vol if vol is not None else np.arange(n, 0, -1) * 1e5})


# ---------------------------------------------------------------- io
def test_normalize_thermo_columns():
    raw = pd.DataFrame({"Monoisotopic Mass": [1.0, 2.0], "Apex RT": [3, 4],
                        "Sum Intensity": [5, 6], "Other": ["x", "y"]})
    out = normalize_columns(raw)
    assert list(out.columns) == ["Mass", "RT", "Vol"]
    assert out.dtypes.eq("float64").all()


def test_normalize_missing_required_column():
    with pytest.raises(ValueError):
        normalize_columns(pd.DataFrame({"Monoisotopic Mass": [1.0]}))


# ---------------------------------------------------------------- base calling
def test_base_calling_reads_a_simple_ladder():
    seq = "GCAU"
    df = frame(ladder_from_seq(seq, start=1000.0))
    calls = base_calling(df)
    called = [p[2] for p in sorted(calls.pairs, key=lambda p: p[0])]
    assert called == list(seq[1:])          # differences between consecutive rungs


def test_homology_finds_methylation():
    df = frame([24000.0, 24000.0 + K.METHYL_BASECALL, 5000.0])
    calls = homology_search(df, names=["Me"])
    assert len(calls.pairs) == 1 and calls.pairs[0][2] == "Me"
    assert 5000.0 not in calls.compounds.Mass.values


# ---------------------------------------------------------------- mass sum
def test_mass_sum_pairs_complementary_fragments():
    full = 10000.0
    five = frame([3000.0, 4000.0, 1234.0])
    three = frame([full + K.H2O - 3000.0, full + K.H2O - 4000.0, 999.0])
    ms3, ms5 = mass_sum(three, five, full_mass=full)
    assert sorted(ms5.Mass) == [3000.0, 4000.0]
    assert 999.0 not in ms3.Mass.values


# ---------------------------------------------------------------- matching
def test_match_masses_ppm_window():
    src = pd.DataFrame({"Mass": [1000.0, 2000.0]})
    tgt = pd.DataFrame({"Mass": [1000.005, 2000.5, 5.0]})
    hits, annot = match_masses(src, tgt, ppm=10)
    assert hits.Mass.tolist() == [1000.005]
    assert annot.Match.tolist() == [True, False]


def test_annotate_keeps_previous_matches_like_original():
    tgt = pd.DataFrame({"Mass": [1000.0, 2000.0]})
    ref = pd.DataFrame({"Mass": [1000.0], "tRNA": ["X"]})
    a = annotate_matches(tgt, ref, ppm=10, copy_cols=["tRNA"])
    b = annotate_matches(a, ref, ppm=10, shift=50, copy_cols=["tRNA"])
    assert b.Match.tolist() == [False, False]
    assert b.tRNA.tolist()[0] == "X"        # stale annotation kept (original behaviour)


def test_drop_dups_keeps_most_intense():
    df = frame([100.2, 100.4, 200.0], vol=[1.0, 5.0, 2.0])
    out = drop_dups(df)
    assert out.Mass.tolist() == [100.4, 200.0]


# ---------------------------------------------------------------- composition
def _original_components_nonempty(mass):
    """Verbatim port of the original brute-force search (slow, used as an oracle)."""
    names = list(K.GAPFILL_BLOCKS)
    vals = np.array(list(K.GAPFILL_BLOCKS.values()))
    mass = abs(mass)
    if mass > 10000:
        return False
    size = round(mass / 320)
    if 200 < mass < 320:
        size = 1
    elif size < 1:
        n = int(round(mass / vals[-1]))
        if abs(n * vals[-1] - mass) < 0.2:
            return -0.1 < mass - n * vals[-1] < 0.1
        return False
    m, n = size, len(names)
    r = m if m < 20 else int(m // 2)
    base = list(filterfalse(lambda i: sum(i) != m, product(range(r + 1), repeat=n - 1)))
    max_methyl = 3 if m < 3 else m
    perms = [b + (i,) for i in range(max_methyl) for b in base if sum(b[:4]) >= i - 1]
    masses = np.matmul(perms, vals)
    masses = masses[masses <= mass + 1]
    d = mass - masses
    return bool(((d > -0.1) & (d < 0.1)).any())


@pytest.mark.parametrize("seed", range(3))
def test_composition_matches_original_bruteforce(seed):
    rng = np.random.default_rng(seed)
    exact = [K.A + K.C, 2 * K.G + K.U + K.GAPFILL_BLOCKS["Methyl"], K.GAPFILL_BLOCKS["Methyl"] * 3]
    masses = np.concatenate([rng.uniform(0, 2600, 40), exact])
    for m in masses:
        assert has_composition(m) == _original_components_nonempty(m), m


# ---------------------------------------------------------------- complementation
def test_complementation_places_bases():
    seq = "GCAUUAGC"
    five = frame(ladder_from_seq(seq, start=K.TAG_5P + K.H2O))
    full = five.Mass.iloc[-1] + 20 * 320
    out = ladder_complementation([Ladder(five, full, 5)])
    called = "".join(b for b in out["Base0"] if b)
    assert called == seq[1:]
