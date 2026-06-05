#!/usr/bin/env python3
"""
Generate input decks (runt.inp + runt.ic) for the hard-core + Yukawa chain /
central-sphere Monte Carlo program `runt`.

The tool sweeps a Cartesian product over the physical parameters and writes one
self-contained case directory per combination. It is not tied to any particular
study: pass whatever chain lengths, sphere sizes, packing fractions and
attraction strengths you want.

Physical model and conventions
------------------------------
* Reduced units: bead diameter sigma = 1.0, bead radius = 0.5.
* sphere_radius (RSP) is the central sphere radius in units of sigma. The
  bead/sphere contact (centre-to-centre) is RSP + 0.5.
* Attraction strength: the site-site / sphere-site potentials are
  beta*u = -eps * exp[-kappa (r - r_contact)] / r, while the code computes
  BATT = BEPS * exp[...] / r. Hence BEPS = -eps: an attractive eps (e.g. 0.5)
  is written as a negative BEPS/BBEPS. eps = 0 is athermal (hard core only).
* Total bead count is held roughly constant across chain lengths by choosing
  nchains = round(target_beads / n); the box length then follows from the
  packing fraction, eta = (pi/6) * n_beads / AL^3.

Initial configuration
---------------------
Beads are laid on a cubic lattice with spacing >= sigma (so no pair is closer
than contact), walked in serpentine order so consecutive beads are bonded
neighbours, with sites inside the central sphere removed. The first n_beads
surviving sites are chunked into chains. The few stretched bonds at lattice
seams / the sphere edge relax during MC equilibration -- the same role the
original move-and-insert generator played. The central sphere sits at the origin
(and its periodic images at the cell corners), matching runt's density
measurement which folds positions toward 0 via DNINT.
"""

from __future__ import annotations
import argparse
import itertools
import math
import os

import numpy as np

SIGMA = 1.0            # bead diameter (reduced units)
BEAD_R = 0.5           # bead radius


# ---------------------------------------------------------------------------
# Geometry / initial-configuration generator (shared with the polyfj tooling)
# ---------------------------------------------------------------------------

def box_length(n_beads: int, eta: float) -> float:
    """Cubic box side from packing fraction: eta = (pi/6) * n_beads / AL^3."""
    return ((math.pi / 6.0) * n_beads / eta) ** (1.0 / 3.0)


def chains_for(n: int, target_beads: int) -> int:
    """Number of chains giving ~target_beads total beads for chain length n."""
    return max(1, round(target_beads / n))


def _min_image(d: np.ndarray, al: float) -> np.ndarray:
    """Nearest-image displacement under cubic PBC."""
    return d - al * np.round(d / al)


def _serpentine_sites(m: int) -> list[tuple[int, int, int]]:
    """Boustrophedon (snake) ordering of an m x m x m integer lattice, so that
    consecutive entries are nearest-neighbour steps almost everywhere (jumps
    occur only at occasional layer/row seams)."""
    sites = []
    for iz in range(m):
        ys = range(m) if iz % 2 == 0 else range(m - 1, -1, -1)
        for k, iy in enumerate(ys):
            xs = range(m) if k % 2 == 0 else range(m - 1, -1, -1)
            for ix in xs:
                sites.append((ix, iy, iz))
    return sites


def build_config(n: int, nchains: int, eta: float, sphere_radius: float,
                 seed: int = 0) -> tuple[float, np.ndarray]:
    """Return (AL, coords[n_beads, 3]) for `nchains` chains of `n` beads on a
    sigma-spaced serpentine lattice with the central sphere excised. `seed` is
    accepted for API stability but unused (the lattice start is deterministic).
    """
    del seed
    n_beads = n * nchains
    al = box_length(n_beads, eta)

    # lattice spacing >= sigma, with a small margin so no pair lands exactly at
    # contact (which floating-point rounding could read as an overlap).
    m = int(al / (SIGMA * 1.0001))
    spacing = al / m
    excl = sphere_radius + BEAD_R          # min centre-to-centre, bead <-> sphere

    order = _serpentine_sites(m)
    pts = (np.asarray(order, dtype=float) + 0.5) * spacing   # centre in cells
    dd = _min_image(pts, al)
    pts = pts[np.einsum("ij,ij->i", dd, dd) >= excl * excl]

    if len(pts) < n_beads:
        raise RuntimeError(
            f"lattice has {len(pts)} usable sites < {n_beads} beads needed "
            f"(n={n}, nchains={nchains}, eta={eta}, sphere_radius={sphere_radius}). "
            f"This eta is too dense to place on a sigma-spaced lattice.")

    return al, np.mod(pts[:n_beads], al)


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def case_name(n: int, sphere_radius: float, eta: float,
              eps_ff: float, eps_sf: float) -> str:
    return (f"n{n}_s{sphere_radius:g}_e{eta:g}"
            f"_ff{eps_ff:g}_sf{eps_sf:g}")


def write_ic(path: str, eta: float, nchains: int, n: int,
             al: float, coords: np.ndarray) -> None:
    """runt.ic layout expected by runt.f:
       PFC / blank / NMOL / N / blank / blank / AL / then (I J X Y Z) per bead."""
    with open(path, "w") as f:
        f.write(f"{eta:.4f}\n\n{nchains}\n{n}\n\n\n{al:.10E}\n")
        k = 0
        for i in range(1, nchains + 1):
            for j in range(1, n + 1):
                x, y, z = coords[k]
                f.write(f"{i} {j} {x:.10E} {y:.10E} {z:.10E}\n")
                k += 1


def write_inp(path: str, ncon: int, nskip: int, bsz: float, dlr: float,
              dint: float, fdick: float, frept: float, beps: float,
              bbeps: float, rsp: float,
              rcut: float | None = None) -> None:
    """runt.inp: one value per line, in the order runt.f READs them.

    Lines 1-10 are always written (exact-energy mode).  When `rcut` is
    provided and > 0, an 11th line containing the cutoff radius is appended,
    enabling the cutoff-energy path in the modern runt binary.
    """
    vals = [str(ncon), str(nskip), f"{bsz:.4f}", f"{dlr:.4f}", f"{dint:.4f}",
            f"{fdick:.6f}", f"{frept:.6f}", f"{beps:.4f}", f"{bbeps:.4f}",
            f"{rsp:.4f}"]
    if rcut is not None and rcut > 0.0:
        vals.append(f"{rcut:.4f}")
    with open(path, "w") as f:
        f.write("\n".join(vals) + "\n")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    default_runs = os.path.join(os.path.dirname(here), "runs")

    ap = argparse.ArgumentParser(
        description="Generate runt input decks over a parameter sweep")
    ap.add_argument("--n", type=int, nargs="+", default=[4, 8, 16],
                    help="chain lengths (beads per chain)")
    ap.add_argument("--sphere", type=float, nargs="+", default=[1.0, 1.5, 2.0],
                    help="central sphere radii (units of sigma)")
    ap.add_argument("--eta", type=float, nargs="+", default=[0.1, 0.2, 0.4],
                    help="packing fractions")
    ap.add_argument("--eps-ff", type=float, nargs="+", default=[0.0],
                    help="chain-chain attraction strengths (>=0 attractive; 0 = athermal)")
    ap.add_argument("--eps-sf", type=float, nargs="+", default=[0.0],
                    help="sphere-chain attraction strengths (>=0 attractive; 0 = athermal)")
    ap.add_argument("--target-beads", type=int, default=1600,
                    help="approximate total bead count (sets nchains per length)")
    ap.add_argument("--ncon", type=int, default=2_000_000,
                    help="total MC moves per run (smoke default; raise for converged runs)")
    ap.add_argument("--nskip", type=int, default=1000, help="moves per accumulation")
    ap.add_argument("--bsz", type=float, default=0.05, help="density bin size")
    ap.add_argument("--dlr", type=float, default=0.1, help="max chain translation")
    ap.add_argument("--dint", type=float, default=0.1, help="max internal displacement")
    ap.add_argument("--outdir", default=default_runs, help="output runs directory")
    ap.add_argument("--rcut", type=float, default=None,
                    help="cutoff radius for the Yukawa potential (modern binary only). "
                         "When provided (> 0), writes an 11th line to runt.inp enabling "
                         "cutoff-energy mode.  Omit (default) for exact-energy mode, "
                         "which is faithful to the legacy binary.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    third = 1.0 / 3.0

    print(f"{'case':40s} {'n':>3} {'Nchain':>6} {'AL':>8} {'BEPS':>6} {'BBEPS':>6} {'RSP':>4}")
    count = 0
    for n, s, eta, ff, sf in itertools.product(
            args.n, args.sphere, args.eta, args.eps_ff, args.eps_sf):
        nchains = chains_for(n, args.target_beads)
        name = case_name(n, s, eta, ff, sf)
        case_dir = os.path.join(args.outdir, name)
        os.makedirs(case_dir, exist_ok=True)

        al, coords = build_config(n=n, nchains=nchains, eta=eta, sphere_radius=s)
        write_ic(os.path.join(case_dir, "runt.ic"), eta, nchains, n, al, coords)
        write_inp(os.path.join(case_dir, "runt.inp"),
                  ncon=args.ncon, nskip=args.nskip, bsz=args.bsz, dlr=args.dlr,
                  dint=args.dint, fdick=third, frept=third,
                  beps=-ff, bbeps=-sf, rsp=s,    # sign flip: eps -> code BEPS
                  rcut=args.rcut)

        print(f"{name:40s} {n:>3} {nchains:>6} {al:8.3f} "
              f"{-ff:6.2f} {-sf:6.2f} {s:4.1f}")
        count += 1

    print(f"\n{count} cases written under {args.outdir}  (NCON={args.ncon})")


if __name__ == "__main__":
    main()
