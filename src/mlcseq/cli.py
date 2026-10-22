"""Command line interface:  mlcseq <command> --help"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

from .basecalling import homology_search
from .io import load_sample
from .mass_sum import mass_sum
from .workflows import load_config, run_phe, run_total

CONFIG_DIR = Path(__file__).resolve().parents[2] / "configs"


def _resolve_config(value: str) -> dict:
    path = Path(value)
    if not path.exists() and (CONFIG_DIR / f"{value}.json").exists():
        path = CONFIG_DIR / f"{value}.json"
    return load_config(path)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="mlcseq", description="MLC-Seq RNA sequencing tools")
    p.add_argument("-v", "--verbose", action="store_true", help="show progress messages")
    sub = p.add_subparsers(dest="cmd", required=True)

    h = sub.add_parser("homology", help="homology search on intact masses of one sample")
    h.add_argument("sample", help=".xlsx or .csv mass list")
    h.add_argument("--top", type=int, default=0, help="only use the N most intense intact compounds")
    h.add_argument("--ppm", type=float, default=10)
    h.add_argument("-o", "--out", help="CSV for the base-call pairs")

    m = sub.add_parser("masssum", help="MassSum on separated 5'/3' zones")
    m.add_argument("zone_5p")
    m.add_argument("zone_3p")
    m.add_argument("--full-mass", type=float, required=True)
    m.add_argument("--top", type=int, default=1000)
    m.add_argument("-o", "--out", required=True, help="output folder")

    ph = sub.add_parser("phe", help="tRNA-Phe workflow (config: phe_lite, phe_full or a .json path)")
    ph.add_argument("--data-root", required=True, help="folder of the original MLC-Seq repo")
    ph.add_argument("--config", default="phe_lite")
    ph.add_argument("-o", "--out", required=True)

    t = sub.add_parser("total", help="yeast total-tRNA workflow")
    t.add_argument("--data-root", required=True)
    t.add_argument("-o", "--out", required=True)

    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO if a.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    t0 = time.perf_counter()

    if a.cmd == "homology":
        df = load_sample(a.sample)
        if a.top:
            df = df[df.Mass > 23000].sort_values("Vol", ascending=False).iloc[:a.top]
        calls = homology_search(df, ppm=a.ppm)
        import pandas as pd
        pairs = pd.DataFrame(calls.pairs, columns=["HeavyMass", "LightMass", "Difference"])
        print(pairs.to_string(index=False) if len(pairs) else "no homologous pairs found")
        if a.out:
            pairs.to_csv(a.out, index=False)
    elif a.cmd == "masssum":
        z5 = load_sample(a.zone_5p).sort_values("Vol", ascending=False).iloc[:a.top]
        z3 = load_sample(a.zone_3p).sort_values("Vol", ascending=False).iloc[:a.top]
        ms3, ms5 = mass_sum(z3, z5, full_mass=a.full_mass)
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        ms5.to_csv(out / "masssum_5p.csv")
        ms3.to_csv(out / "masssum_3p.csv")
        print(f"5' ladder: {len(ms5)} fragments, 3' ladder: {len(ms3)} fragments -> {out}")
    elif a.cmd == "phe":
        res = run_phe(a.data_root, _resolve_config(a.config), a.out)
        filled = (res.complementation.filter(like="Base") != "").any(axis=1).sum()
        print(f"positions with a base call: {filled}/{len(res.complementation)} -> {a.out}")
    elif a.cmd == "total":
        res = run_total(a.data_root, a.out)
        n = res.control.MatchedMass.notna().sum()
        print(f"{n}/{len(res.control)} intact masses matched the reference -> {a.out}")

    print(f"done in {time.perf_counter() - t0:.1f}s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
