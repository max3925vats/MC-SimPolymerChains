#!/usr/bin/env python3
"""
Plot runt results: overlay the bead number-density profile near the central
sphere for a chosen swept parameter, holding the others fixed. General over any
generated cases -- not tied to a particular study.

Examples
--------
    # density vs sphere radius, at fixed chain length and packing (athermal):
    python3 tools/plot_results.py --vary sphere --fix n=8 --fix eta=0.2 \
        --fix eps_ff=0 --fix eps_sf=0

    # effect of sphere-chain attraction, fixed geometry:
    python3 tools/plot_results.py --vary eps_sf --fix n=8 --fix sphere=1.5 \
        --fix eta=0.2 --fix eps_ff=0
"""

from __future__ import annotations
import argparse
import os

import plotlib

_HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(_HERE, "..", "runs")
FIGS = os.path.join(_HERE, "..", "figures")

# observable -> (output file, column index, is_density, y-axis label)
OBSERVABLES = {
    "density": ("denrunt.out", 2, True, r"$\rho(r)/\rho_0$"),
}


def main():
    ap = argparse.ArgumentParser(description="Overlay runt density profiles")
    plotlib.add_common_args(ap, OBSERVABLES)
    ap.add_argument("--outdir", default=FIGS)
    args = ap.parse_args()

    out = plotlib.overlay(
        RUNS, "runt.ic", OBSERVABLES,
        observable=args.observable, vary=args.vary,
        fix=plotlib.parse_fix(args.fix),
        xmin=args.xmin, xmax=args.xmax, outdir=args.outdir)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
