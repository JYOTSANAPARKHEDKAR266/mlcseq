"""Run the ORIGINAL MLC-Seq notebook code (unchanged algorithms) headlessly and
save every intermediate result, so the refactor can be checked against it.

Usage: python run_original.py <repo_dir> <out_dir>
Plotting libraries are replaced by no-op mocks; algorithm code is executed as-is.
"""
import io
import json
import os
import re
import sys
import contextlib
from unittest import mock

for name in ["plotly", "plotly.express", "plotly.graph_objects", "matplotlib",
             "matplotlib.pyplot", "seaborn"]:
    sys.modules[name] = mock.MagicMock()

import numpy as np
import pandas as pd

REPO, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)


def nb_code(path):
    nb = json.load(open(path))
    cells = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
    return cells


def clean(src):
    src = re.sub(r"^%%capture.*$", "", src, flags=re.M)
    src = re.sub(r"^%run.*$", "", src, flags=re.M)
    return src


ns = {"pd": pd, "np": np, "__name__": "orig"}
os.chdir(os.path.join(REPO, "examples"))
for mod in ["utils", "ladder_separation", "homology_search", "mass_sum", "gap_fill",
            "ladder_complementation"]:
    for cell in nb_code(os.path.join(REPO, "modules", mod + ".ipynb")):
        exec(clean(cell), ns)
# no-op plotting
for fn in ["plotly_zone", "plotly_zones", "plotly_multi_zones", "plotly_basecalling",
           "plot_zones", "plot_basecalling"]:
    ns[fn] = lambda *a, **k: (mock.MagicMock(), mock.MagicMock())
ns["plt"] = mock.MagicMock()


def save(df, name):
    df.to_csv(os.path.join(OUT, name + ".csv"))


def run_notebook(nb_name, skip_cells=()):
    cells = nb_code(os.path.join(REPO, "examples", nb_name))
    log = io.StringIO()
    for i, cell in enumerate(cells):
        if i in skip_cells:
            continue
        src = clean(cell)
        # never write into the repo's outputs folder
        src = src.replace("'../outputs/", f"'{OUT}/orig_")
        with contextlib.redirect_stdout(log):
            exec(src, ns)
    return log.getvalue()


which = sys.argv[3] if len(sys.argv) > 3 else "all"

if which in ("all", "lite"):
    run_notebook("trna_phe_analysis_lite.ipynb")
    save(ns["bcr"][0], "lite_homology_ctrl")
    json.dump([list(p) for p in ns["bcr"][1]], open(f"{OUT}/lite_homology_ctrl_pairs.json", "w"))
    save(ns["df_homo_deg"], "lite_homology_deg")
    save(ns["df_masssum_5p"], "lite_masssum_5p")
    save(ns["df_masssum_3p"], "lite_masssum_3p")
    save(ns["df_gap_5p"], "lite_gap_5p")
    save(ns["df_gap_3p"], "lite_gap_3p")
    save(ns["df_ladder_comp"], "lite_ladder_comp")

if which in ("all", "full"):
    # skip the interactive lasso-selection cells (widget)
    cells = nb_code(os.path.join(REPO, "examples", "trna_phe_analysis.ipynb"))
    skip = {i for i, c in enumerate(cells) if "zone_selection(" in c or "df_sample.iloc[idxs]" in c}
    run_notebook("trna_phe_analysis.ipynb", skip_cells=skip)
    for i in range(5):
        save(ns["df_masssum_5ps"][i], f"full_masssum_5p_iso{i+1}")
        save(ns["df_masssum_3ps"][i], f"full_masssum_3p_iso{i+1}")
    for i in range(4):
        save(ns["df_ms_gf_5ps_refined"][i], f"full_final_5p_iso{i+1}")
        save(ns["df_ms_gf_3ps_refined"][i], f"full_final_3p_iso{i+1}")
    save(ns["df_gap_1st_5p"], "full_gap_5p_iso1")
    save(ns["df_gap_1st_3p"], "full_gap_3p_iso1")
    save(ns["df_gap_2nd_5p"], "full_gap_5p_iso2")
    save(ns["df_gap_2nd_3p"], "full_gap_3p_iso2")
    save(ns["df_ladder_comp"], "full_ladder_comp")

if which in ("all", "total"):
    log = run_notebook("trna_total_analysis.ipynb")
    save(ns["df_ctrl"], "total_ctrl_annotated")
    save(ns["dfms_unique"], "total_dfms_unique")
    save(ns["bcr"][0], "total_homology_me")
    json.dump([list(p) for p in ns["bcr"][1]], open(f"{OUT}/total_homology_me_pairs.json", "w"))
    coverage = [l for l in log.splitlines() if re.match(r"^\d+ ", l)]
    open(f"{OUT}/total_coverage.txt", "w").write("\n".join(coverage) + "\n")

print("done", which)
