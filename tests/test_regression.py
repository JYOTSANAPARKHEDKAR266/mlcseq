"""End-to-end checks: the refactor must reproduce the ORIGINAL code's outputs.

Reference files in tests/reference/ were produced by tools/run_original.py, which
runs the original notebooks unchanged under their pinned 2021 environment
(Python 3.8, pandas 1.2.3, numpy 1.20.1).

Needs the original repository for input data:
    git clone https://github.com/rnamodifications/MLC-Seq
    MLCSEQ_DATA_ROOT=/path/to/MLC-Seq pytest -m slow
"""
import ast
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from mlcseq.workflows import load_config, run_phe, run_total

HERE = Path(__file__).parent
REF = HERE / "reference"
CONFIGS = HERE.parent / "configs"
DATA = os.environ.get("MLCSEQ_DATA_ROOT")

pytestmark = [pytest.mark.slow,
              pytest.mark.skipif(not DATA, reason="set MLCSEQ_DATA_ROOT to the original repo")]


def ref(name):
    return pd.read_csv(REF / f"{name}.csv", index_col=0)


def same_rows(expected: pd.DataFrame, actual: pd.DataFrame):
    """Same compounds (index labels), order-insensitive."""
    assert sorted(expected.index) == sorted(actual.index)


def cell_set(v):
    if pd.isna(v) or v == "":
        return ()
    return tuple(sorted(str(v).split(", ")))


def same_table(expected: pd.DataFrame, actual: pd.DataFrame):
    """Cell-by-cell, treating 'a, b' cells as unordered lists (the original's order
    came from Python set iteration and carries no meaning)."""
    assert expected.shape == actual.shape
    expected = expected.set_axis(expected.index.astype(int))
    for pos in expected.index:
        for col in expected.columns:
            e = expected.loc[pos, col]
            a = actual.loc[pos, col]
            assert cell_set(e) == cell_set(a), (pos, col, e, a)


@pytest.fixture(scope="module")
def lite():
    return run_phe(DATA, load_config(CONFIGS / "phe_lite.json"))


@pytest.fixture(scope="module")
def full():
    return run_phe(DATA, load_config(CONFIGS / "phe_full.json"))


@pytest.fixture(scope="module")
def total():
    return run_total(DATA)


def test_lite_homology(lite):
    same_rows(ref("lite_homology_ctrl"), lite.homology_control.compounds)
    same_rows(ref("lite_homology_deg"), lite.homology_degraded.compounds)
    expected = json.loads((REF / "lite_homology_ctrl_pairs.json").read_text())
    assert sorted(map(tuple, expected)) == sorted(lite.homology_control.pairs)


def test_lite_masssum_and_gapfill(lite):
    same_rows(ref("lite_masssum_5p"), lite.masssum_5p[0])
    same_rows(ref("lite_masssum_3p"), lite.masssum_3p[0])
    same_rows(ref("lite_gap_5p"), lite.gap_5p[1])
    same_rows(ref("lite_gap_3p"), lite.gap_3p[1])


def test_lite_complementation(lite):
    same_table(pd.read_csv(REF / "lite_ladder_comp.csv", index_col=0, dtype=str, keep_default_na=False),
               lite.complementation)


def test_full_every_stage(full):
    for i in range(5):
        same_rows(ref(f"full_masssum_5p_iso{i+1}"), full.masssum_5p[i])
        same_rows(ref(f"full_masssum_3p_iso{i+1}"), full.masssum_3p[i])
    for i in (1, 2):
        same_rows(ref(f"full_gap_5p_iso{i}"), full.gap_5p[i])
        same_rows(ref(f"full_gap_3p_iso{i}"), full.gap_3p[i])
    for i in range(1, 5):
        same_rows(ref(f"full_final_5p_iso{i}"), full.ladders_5p[i])
        same_rows(ref(f"full_final_3p_iso{i}"), full.ladders_3p[i])
    same_table(pd.read_csv(REF / "full_ladder_comp.csv", index_col=0, dtype=str, keep_default_na=False),
               full.complementation)


def test_total_matching(total):
    exp = ref("total_ctrl_annotated")
    cols = ["Mass", "RT", "Vol", "MatchedMass", "Sft", "tRNA", "PPM", "Isoform", "Match"]
    pd.testing.assert_frame_equal(exp[cols], total.control[cols], check_dtype=False)
    assert sorted(ref("total_dfms_unique").Mass) == sorted(total.matched_unique.Mass)
    same_rows(ref("total_homology_me"), total.methyl_homology.compounds)


def test_total_coverage(total):
    cov = total.coverage
    for line in (REF / "total_coverage.txt").read_text().splitlines():
        n, trna, rate, sites = line.split(" ", 3)
        row = cov.loc[int(n)]
        assert row.tRNA == trna
        assert f"{row.coverage:.2f}" == rate
        assert sorted(int(s) for s in ast.literal_eval(sites)) == \
            [int(s) for s in str(row.positions).split()]
