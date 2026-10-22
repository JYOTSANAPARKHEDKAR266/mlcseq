# mlcseq
[![tests](https://github.com/JYOTSANAPARKHEDKAR266/mlcseq/actions/workflows/tests.yml/badge.svg)](https://github.com/JYOTSANAPARKHEDKAR266/mlcseq/actions/workflows/tests.yml)

A cleaned-up, tested Python package of the **MLC-Seq** algorithms from
*Yuan et al., "Mass Spectrometry-Based Direct Sequencing of tRNAs De Novo and
Quantitative Mapping of Multiple RNA Modifications", JACS 2024*.
The original code is at https://github.com/rnamodifications/MLC-Seq (GPL-3.0).

**The results are the same as the original code, it runs 3–60× faster, and it
installs on current Python.**

## Results compared with the original

The original notebooks were run unchanged in their own pinned 2021 environment
(Python 3.8, pandas 1.2.3, numpy 1.20.1) and every intermediate result was saved
(`tools/run_original.py` → `tests/reference/`). The test suite checks that this
package gives the same output at every stage:

| Workflow | Stages compared | Same as original? | Original | This package |
|---|---|---|---|---|
| tRNA-Phe lite | homology, MassSum, GapFill, complementation | yes | 30 s | 1.2 s |
| tRNA-Phe full (5 isoforms, all manual filters) | all 5 isoforms, every stage | yes | 69 s | 1.2 s |
| Yeast total tRNA | database matching, methyl homology, coverage of 34 tRNAs | yes | 34 s | 10 s |

One small, intentional difference: when several fragments or bases land on the
same position, the original listed them in Python-`set` order (effectively
random). This package lists them in a fixed order. The tests compare those cells
as unordered lists.

On a current install (pandas 3), the original code ran for several minutes and
then crashed in the last step (ladder complementation).

## Install

```bash
git clone https://github.com/rnamodifications/MLC-Seq          # input data
cd mlcseq
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"            # add ,plot for plotly figures
```

## Run

```bash
# tRNA-Phe, short version (same as trna_phe_analysis_lite.ipynb)
mlcseq phe --data-root ../MLC-Seq --config phe_lite -o out/lite

# tRNA-Phe, full version with all manual decisions from the paper
mlcseq phe --data-root ../MLC-Seq --config phe_full -o out/full

# yeast total tRNA
mlcseq total --data-root ../MLC-Seq -o out/total

# single steps on your own files (.xlsx or .csv)
mlcseq homology sample.xlsx --top 10
mlcseq masssum zone5p.csv zone3p.csv --full-mass 24252.311 -o out/ms
```

From Python:

```python
from mlcseq import load_sample, mass_sum, gap_fill, ladder_complementation, Ladder
```

## Tests

```bash
pytest                                              # fast unit tests
MLCSEQ_DATA_ROOT=../MLC-Seq pytest                  # + full regression vs original
```

## What changed

| Original | Now |
|---|---|
| Jupyter notebooks loading each other with `%run` | Python package with modules, a command-line tool and `pip install` |
| 80 pinned 2021 packages; doesn't install on new Python | 3 dependencies (numpy, pandas, openpyxl); plotly optional |
| Manual filtering hard-coded in notebook cells | Every manual decision in a JSON config (`configs/phe_full.json`) |
| GapFill rebuilt a brute-force table of millions of tuples on every call | Table built once per gap size with numpy, cached, binary search (same search space) |
| Row-by-row `iterrows` mass matching | Sorted arrays + binary search |
| `networkx` only to count neighbours | Plain dictionaries |
| Functions silently modify their inputs (`inplace=True`) | Functions return new tables |
| Relied on pandas chained assignment (`df.iloc[1]['Mass'] += ...`), which silently does nothing in pandas 3 | Explicit values |
| Excel input only | Excel or CSV |
| No tests | 12 unit tests (including a check of the composition search against the original brute force) + 6 end-to-end regression tests |

## Quirks kept on purpose (to match the paper's results)

These are kept so the output matches the paper. Each one is a candidate for a
later improvement.

1. **Position = mass // 320.** Every fragment's position uses an average
   nucleotide mass of 320 Da. On a 76-nt tRNA with modifications, this can be
   off by one or more positions.
2. **Three different methyl masses** (14.01, 14.0106, 14.0156) are used in
   different steps. The correct CH₂ value is 14.01565. See `constants.py`.
3. **First match, not best match.** Database matching takes the first reference
   within tolerance, not the closest one.
4. **Old matches are not cleared.** In the total-tRNA workflow, `Match` is reset
   for each shift tried, but `MatchedMass`/`tRNA` from earlier shifts are kept.
5. **The GapFill `major` flag is described backwards** in the original
   docstring. `major=True` limits gaps to 3200 Da; the default is 10000 Da.
6. **Ladder separation (5′ vs 3′) is still manual.** The lasso widget is in
   `mlcseq.plotting.lasso_selector`. The examples use the authors' saved
   selections (`examples/data/phe5p.xlsx`, `phe3p.xlsx`).

## Layout

```
src/mlcseq/
  constants.py        masses and tolerances
  io.py               load Thermo exports (.xlsx/.csv)
  matching.py         ppm matching, de-duplication
  basecalling.py      base calling and homology search
  mass_sum.py         MassSum
  composition.py      which nucleotide sets explain a mass gap
  gap_fill.py         GapFill
  complementation.py  ladder complementation
  workflows.py        the three example workflows
  plotting.py         optional plotly figures and lasso selector
  cli.py              command-line tool
configs/              workflow settings and manual decisions
tests/                unit and regression tests, reference outputs
tools/run_original.py runs the original notebooks headlessly
```

## License and credit

GPL-3.0, the same as the original. All algorithms are by the original authors
(Shenglong Zhang lab). This repository only re-implements and tests their code.
