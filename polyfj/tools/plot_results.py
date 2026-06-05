#!/usr/bin/env python3
"""
Plot polyfj results: overlay a chosen observable against r/sigma for a swept
parameter, holding the others fixed. General over any generated cases -- not
tied to a particular study.

Observables: the total bead density and its end-/mid-bead resolved components;
the centre-of-mass-resolved radius-of-gyration and end-to-end mean squares; the
segmental order parameter; and the chain-shape semi-axes.

Examples
--------
    # density vs sphere radius, fixed chain length and packing:
    python3 tools/plot_results.py --vary sphere --fix n=8 --fix eta=0.2

    # density vs chain length, fixed geometry:
    python3 tools/plot_results.py --vary n --fix sphere=1.5 --fix eta=0.2

    # radius-of-gyration profile vs chain length:
    python3 tools/plot_results.py --observable rg2 --vary n \
        --fix sphere=1.5 --fix eta=0.2
"""

from __future__ import annotations
import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "runt", "tools"))
import plotlib  # noqa: E402

RUNS = os.path.join(_HERE, "..", "runs")
FIGS = os.path.join(_HERE, "..", "figures")

# observable -> (output file, column index, is_density, y-axis label)
OBSERVABLES = {
    "density":      ("polyden.out", 2, True,  r"$\rho(r)/\rho_0$"),
    "density_ends": ("polyden.out", 3, True,  r"$\rho_{\mathrm{ends}}(r)/\rho_0$"),
    "density_mid":  ("polyden.out", 4, True,  r"$\rho_{\mathrm{mid}}(r)/\rho_0$"),
    # rgfj.out layout (modern driver): idx dist rho_com Rg2 Re2  (5 cols, 0-based)
    "rho_com":      ("rgfj.out",    2, False, r"$\rho_{\mathrm{CoM}}(r)$"),
    "rg2":          ("rgfj.out",    3, False, r"$\langle R_g^2\rangle/\sigma^2$"),
    "re2":          ("rgfj.out",    4, False, r"$\langle R_e^2\rangle/\sigma^2$"),
    # seg.out: idx dist S  (3 cols, 0-based col 2)
    "seg":          ("seg.out",     2, False, r"$S(r)$"),
    # g2fj.out layout (modern driver): idx dist a b c  (5 cols, 0-based cols 2/3/4)
    "shape_a":      ("g2fj.out",    2, False, r"$a/\sigma$"),
    "shape_b":      ("g2fj.out",    3, False, r"$b/\sigma$"),
    "shape_c":      ("g2fj.out",    4, False, r"$c/\sigma$"),
}


def main():
    ap = argparse.ArgumentParser(description="Overlay polyfj profiles")
    plotlib.add_common_args(ap, OBSERVABLES)
    ap.add_argument("--outdir", default=FIGS)
    args = ap.parse_args()

    out = plotlib.overlay(
        RUNS, "polyfj.ic", OBSERVABLES,
        observable=args.observable, vary=args.vary,
        fix=plotlib.parse_fix(args.fix),
        xmin=args.xmin, xmax=args.xmax, outdir=args.outdir)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
