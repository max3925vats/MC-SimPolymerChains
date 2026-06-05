#!/usr/bin/env python3
"""
Generate input decks (polyfj.inp + polyfj.ic) for the hard-chain /
hard-sphere Monte Carlo program `polyfj`.

The tool sweeps a Cartesian product over the physical parameters and writes one
self-contained case directory per combination. It is not tied to any particular
study: pass whatever chain lengths, sphere sizes and packing fractions you want.

polyfj is athermal (hard core only, no attractions), so the sweep covers chain
length, central sphere radius and packing fraction. In addition to the density
profile, polyfj outputs structural data: the radius-of-gyration and end-to-end
tensors, the segmental order parameter, and the chain-shape semi-axes.

Conventions match the runt tooling (reduced units, sigma = 1; sphere_radius is
the central sphere radius in units of sigma; nchains = round(target_beads / n);
box length follows from eta). The geometry / initial-configuration generator is
shared -- see ../../runt/tools/setup_runs.py. The only format differences here:
polyfj.ic carries the sphere radius RSP, and polyfj.inp has 7 values (no
temperature or attraction parameters).
"""

from __future__ import annotations
import argparse
import itertools
import os
import sys

# reuse the validated geometry / generator from the runt tooling
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "runt", "tools"))
from setup_runs import build_config, chains_for  # noqa: E402


def case_name(n: int, sphere_radius: float, eta: float) -> str:
    return f"n{n}_s{sphere_radius:g}_e{eta:g}"


def write_ic(path: str, eta: float, nchains: int, n: int, rsp: float,
             al, coords) -> None:
    """polyfj.ic layout expected by polyfj.f:
       PFC / blank / NMOL1 / N / blank / RSP / AL / blank / (I J X Y Z) per bead."""
    with open(path, "w") as f:
        f.write(f"{eta:.4f}\n\n{nchains}\n{n}\n\n{rsp:.4f}\n{al:.10E}\n\n")
        k = 0
        for i in range(1, nchains + 1):
            for j in range(1, n + 1):
                x, y, z = coords[k]
                f.write(f"{i} {j} {x:.10E} {y:.10E} {z:.10E}\n")
                k += 1


def write_inp(path: str, ncon: int, nskip: int, bsz: float, dlr: float,
              dint: float, fdick: float, frept: float) -> None:
    """polyfj.inp: 7 values, one per line, in polyfj.f's READ order."""
    vals = [str(ncon), str(nskip), f"{bsz:.4f}", f"{dlr:.4f}", f"{dint:.4f}",
            f"{fdick:.6f}", f"{frept:.6f}"]
    with open(path, "w") as f:
        f.write("\n".join(vals) + "\n")


def main() -> None:
    default_runs = os.path.join(os.path.dirname(_HERE), "runs")

    ap = argparse.ArgumentParser(
        description="Generate polyfj input decks over a parameter sweep")
    ap.add_argument("--n", type=int, nargs="+", default=[4, 8, 16],
                    help="chain lengths (beads per chain)")
    ap.add_argument("--sphere", type=float, nargs="+", default=[1.0, 1.5, 2.0],
                    help="central sphere radii (units of sigma)")
    ap.add_argument("--eta", type=float, nargs="+", default=[0.1, 0.2, 0.4],
                    help="packing fractions")
    ap.add_argument("--target-beads", type=int, default=1600,
                    help="approximate total bead count (sets nchains per length)")
    ap.add_argument("--ncon", type=int, default=2_000_000,
                    help="total MC moves per run (smoke default; raise for converged runs)")
    ap.add_argument("--nskip", type=int, default=1000, help="moves per accumulation")
    ap.add_argument("--bsz", type=float, default=0.05, help="density bin size")
    ap.add_argument("--dlr", type=float, default=0.1, help="max chain translation")
    ap.add_argument("--dint", type=float, default=0.1, help="max internal displacement")
    ap.add_argument("--outdir", default=default_runs, help="output runs directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    third = 1.0 / 3.0

    print(f"{'case':22s} {'n':>3} {'Nchain':>6} {'AL':>8} {'RSP':>4}")
    count = 0
    for n, s, eta in itertools.product(args.n, args.sphere, args.eta):
        nchains = chains_for(n, args.target_beads)
        name = case_name(n, s, eta)
        case_dir = os.path.join(args.outdir, name)
        os.makedirs(case_dir, exist_ok=True)

        al, coords = build_config(n=n, nchains=nchains, eta=eta, sphere_radius=s)
        write_ic(os.path.join(case_dir, "polyfj.ic"), eta, nchains, n, s, al, coords)
        write_inp(os.path.join(case_dir, "polyfj.inp"),
                  ncon=args.ncon, nskip=args.nskip, bsz=args.bsz,
                  dlr=args.dlr, dint=args.dint, fdick=third, frept=third)

        print(f"{name:22s} {n:>3} {nchains:>6} {al:8.3f} {s:4.1f}")
        count += 1

    print(f"\n{count} cases written under {args.outdir}  (NCON={args.ncon})")


if __name__ == "__main__":
    main()
