"""End-to-end workflows that reproduce the three example notebooks of the
original repository, with every manual decision moved into a JSON config."""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path

import pandas as pd

from . import constants as K
from .basecalling import BaseCalls, homology_search
from .complementation import Ladder, ladder_complementation
from .gap_fill import gap_fill
from .io import load_many, load_sample
from .mass_sum import mass_sum
from .matching import annotate_matches, drop_dups, match_masses

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------- helpers
def top_by_intensity(df: pd.DataFrame, n: int) -> pd.DataFrame:
    return df.sort_values("Vol", ascending=False).iloc[:n]


def intact_dot(df: pd.DataFrame, full_mass: float, tol: float = 0.1) -> pd.DataFrame:
    return df[(df.Mass > full_mass - tol) & (df.Mass < full_mass + tol)]


def apply_exclusions(df: pd.DataFrame, rule: dict | None) -> pd.DataFrame:
    """Drop rows by index label and/or mass cut-offs.

    rule = {"ids": [...], "mass_below": x, "mass_above": y}  (all optional)
    """
    if not rule:
        return df
    drop = df.index.isin(rule.get("ids", []))
    if "mass_below" in rule:
        drop |= (df.Mass < rule["mass_below"]).to_numpy()
    if "mass_above" in rule:
        drop |= (df.Mass > rule["mass_above"]).to_numpy()
    missing = set(rule.get("ids", [])) - set(df.index)
    if missing:
        log.warning("exclusion ids not present in ladder: %s", sorted(missing))
    return df[~drop]


def _write(df: pd.DataFrame, out_dir: Path | None, name: str) -> None:
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_dir / f"{name}.csv")


# --------------------------------------------------------------------------- tRNA-Phe
@dataclass
class PheResult:
    homology_control: BaseCalls
    homology_degraded: BaseCalls
    masssum_5p: list[pd.DataFrame] = field(default_factory=list)
    masssum_3p: list[pd.DataFrame] = field(default_factory=list)
    gap_5p: dict[int, pd.DataFrame] = field(default_factory=dict)
    gap_3p: dict[int, pd.DataFrame] = field(default_factory=dict)
    ladders_5p: dict[int, pd.DataFrame] = field(default_factory=dict)
    ladders_3p: dict[int, pd.DataFrame] = field(default_factory=dict)
    complementation: pd.DataFrame | None = None


def run_phe(data_root: str | Path, config: dict, out_dir: str | Path | None = None) -> PheResult:
    """Yeast tRNA-Phe workflow (both the 'lite' and the full notebook).

    Steps: homology search -> MassSum per isoform -> manual exclusions (from config)
    -> GapFill -> manual exclusions -> ladder complementation.
    """
    root = Path(data_root)
    out = Path(out_dir) if out_dir else None

    ctrl = load_sample(root / config["control"])
    homo_ctrl = homology_search(ctrl)

    deg = load_sample(root / config["degraded"])
    deg_top = top_by_intensity(deg[deg.Mass > K.INTACT_MIN_MASS], config.get("top_intact", 10))
    homo_deg = homology_search(deg_top)

    zone5 = load_sample(root / config["zone_5p"])
    zone3 = load_sample(root / config["zone_3p"])
    top5 = top_by_intensity(zone5, config.get("top_n", 1000))
    top3 = top_by_intensity(zone3, config.get("top_n", 1000))

    isoforms = config.get("isoforms")
    if not isoforms:  # lite mode: take the most intense intact compound
        isoforms = [float(homo_deg.compounds.sort_values("Vol", ascending=False).iloc[0].Mass)]
    res = PheResult(homo_ctrl, homo_deg)

    ms_excl = config.get("masssum_exclude", {})
    gf_excl = config.get("gapfill_exclude", {})
    gapfill_isoforms = config.get("gapfill_isoforms", list(range(1, len(isoforms) + 1)))
    ladder_isoforms = config.get("ladder_isoforms", gapfill_isoforms)

    for n, full in enumerate(isoforms, start=1):
        ms3, ms5 = mass_sum(top3, top5, full_mass=full)
        res.masssum_5p.append(ms5)
        res.masssum_3p.append(ms3)
        _write(ms5, out, f"masssum_5p_iso{n}")
        _write(ms3, out, f"masssum_3p_iso{n}")

        rules = ms_excl.get(str(n), {})
        lad5 = apply_exclusions(ms5, rules.get("5p"))
        lad3 = apply_exclusions(ms3, rules.get("3p"))

        if n in gapfill_isoforms:
            log.info("GapFill isoform %d (%.4f)", n, full)
            g5 = gap_fill(zone5, lad5, intact_dot(zone5, full), orientation=5, major=True)
            g3 = gap_fill(zone3, lad3, intact_dot(zone3, full), orientation=3, major=True)
            res.gap_5p[n], res.gap_3p[n] = g5, g3
            _write(g5, out, f"gapfill_5p_iso{n}")
            _write(g3, out, f"gapfill_3p_iso{n}")
            rules = gf_excl.get(str(n), {})
            lad5 = apply_exclusions(pd.concat([lad5, g5]), rules.get("5p"))
            lad3 = apply_exclusions(pd.concat([lad3, g3]), rules.get("3p"))

        res.ladders_5p[n], res.ladders_3p[n] = lad5, lad3
        _write(lad5, out, f"ladder_5p_iso{n}")
        _write(lad3, out, f"ladder_3p_iso{n}")

    ladders = []
    for n in ladder_isoforms:
        full = isoforms[n - 1]
        ladders += [Ladder(res.ladders_5p[n], full, 5), Ladder(res.ladders_3p[n], full, 3)]
    res.complementation = ladder_complementation(ladders)
    _write(res.complementation, out, "ladder_complementation")
    if out is not None:
        res.complementation.to_excel(out / "ladder_complementation.xlsx")
    return res


# --------------------------------------------------------------------------- total tRNA
SEQBASE_COLUMNS = ["Monoisotopic Mass (full length)", "Monoisotopic Mass (CC-tailed)",
                   "Monoisotopic Mass (C-tailed)", "AminoAcid & CCA"]


@dataclass
class TotalResult:
    control: pd.DataFrame         # merged control compounds with match annotations
    matched_unique: pd.DataFrame  # compounds matching any reference intact mass
    methyl_homology: BaseCalls
    coverage: pd.DataFrame        # ladder coverage per reference tRNA


def merge_controls(files) -> pd.DataFrame:
    """Load control runs, keep intact masses, align RT/intensity scales, de-duplicate."""
    dfs = [d[d.Mass > K.INTACT_MIN_MASS].copy() for d in load_many(files)]
    max_rt = max(d.RT.max() for d in dfs)
    max_vol = max(d.Vol.max() for d in dfs)
    for d in dfs:
        d["RT"] = d.RT * (max_rt / d.RT.max())
        d["Vol"] = d.Vol * (max_vol / d.Vol.max())
    return drop_dups(pd.concat(dfs))


def ladder_coverage(df_ad: pd.DataFrame, library: pd.DataFrame, ppm: float = 30) -> pd.DataFrame:
    """For each reference tRNA, which positions are supported by a matching fragment."""
    ladder_cols = ["5prim_df", "3prim_df", "3prim_CCtailed_df", "3prim_Ctailed_df"]
    full_cols = [None, "MassCCA", "MassCC", "MassC"]
    rows = []
    for i, (_, rec) in enumerate(library.iterrows(), start=1):
        seq = [c for c in rec.iloc[1:93] if c != "-"]
        trna, seq = seq[0], seq[1:]
        sites: set[float] = set()
        for col, full_col in zip(ladder_cols, full_cols):
            ladder = pd.read_csv(StringIO(rec[col]), sep="\t")
            hits, _ = match_masses(ladder, df_ad, ppm=ppm)
            hits = drop_dups(hits)
            mass = hits.Mass if full_col is None else rec[full_col] - hits.Mass
            sites.update(float(s) for s in (mass // K.POSITION_MASS))
        sites = sorted(s for s in sites if s > 0)
        rows.append({"row": i, "tRNA": trna, "length": len(seq),
                     "coverage": len(sites) / len(seq),
                     "positions": " ".join(str(int(s)) for s in sites)})
    return pd.DataFrame(rows).set_index("row")


def run_total(data_root: str | Path, out_dir: str | Path | None = None) -> TotalResult:
    """Yeast total-tRNA workflow (database matching instead of de novo reading)."""
    root = Path(data_root)
    out = Path(out_dir) if out_dir else None

    ctrl = merge_controls((root / "samples/total_tRNA/Control").glob("*.xls*"))

    seqbase_file = root / "statics/total_tRNA_seqs_base.xlsx"
    seqbase = pd.read_excel(seqbase_file)
    matched = []
    for col in SEQBASE_COLUMNS:
        ref = seqbase.assign(Mass=seqbase[col])
        hits, _ = match_masses(ref, ctrl, ppm=50)
        matched.append(hits.sort_values("Mass").drop(columns="RoundedMass"))
    matched_unique = drop_dups(pd.concat(matched))

    refs = {name: pd.read_excel(seqbase_file, sheet_name=name) for name in ("CCA", "CC", "C")}
    for shift in [2 * K.METHYL_UTILS, K.METHYL_UTILS, K.NA, K.K,
                  -2 * K.METHYL_UTILS, -1 * K.METHYL_UTILS, 0]:
        for name in ("C", "CCA", "CC"):   # same order as the original notebook
            ctrl = annotate_matches(ctrl, refs[name], ppm=30, shift=shift, copy_cols=["tRNA"])
            ctrl.loc[ctrl.Match.to_numpy(), "Isoform"] = name

    report = ctrl.sort_values("tRNA").reset_index(drop=True)
    report.index += 1
    report = report[["Mass", "RT", "Vol", "MatchedMass", "tRNA", "PPM", "Isoform", "Sft"]]

    methyl = homology_search(ctrl[ctrl.MatchedMass > 0], names=["Me"], ppm=10)

    df_ad = pd.concat(load_many((root / "samples/total_tRNA/AcidDeg").glob("*.xls*")))
    library = pd.read_excel(root / "statics/total_tRNA_mass_ladders.xlsx", skiprows=1)
    coverage = ladder_coverage(df_ad, library)

    _write(ctrl, out, "control_annotated")
    _write(report, out, "matched_control_report")
    _write(matched_unique, out, "matched_unique")
    _write(methyl.compounds, out, "methyl_homology")
    _write(coverage, out, "ladder_coverage")
    return TotalResult(ctrl, matched_unique, methyl, coverage)


def load_config(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())
