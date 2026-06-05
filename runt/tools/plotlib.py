#!/usr/bin/env python3
"""
Shared plotting core for the polymer-chain / central-sphere MC tooling.

Reads the per-case output directories produced by run_all.sh and draws an
overlay of a chosen observable against the radial coordinate r/sigma: one curve
per value of a chosen "vary" parameter, with the remaining parameters pinned by
"fix" constraints. It is agnostic to any particular study -- it plots whatever
cases exist on disk.

Case directories are named by their physical parameters, e.g.
    n8_s1.5_e0.2            (polyfj: chain length, sphere radius, packing)
    n8_s1.5_e0.2_ff0_sf0.5  (runt: + chain-chain and sphere-chain attraction)
which this module parses back into a parameter dict with canonical keys
n, sphere, eta, eps_ff, eps_sf.
"""

from __future__ import annotations
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# token prefix in a case-dir name -> canonical parameter key
_PREFIX = {"n": "n", "s": "sphere", "e": "eta", "ff": "eps_ff", "sf": "eps_sf"}
_TOKEN = re.compile(r"^(ff|sf|n|s|e)([0-9.]+)$")


def parse_params(name: str) -> dict | None:
    """Parse a case-dir name into {n, sphere, eta, eps_ff, eps_sf} (those present)."""
    params: dict = {}
    for tok in name.split("_"):
        m = _TOKEN.match(tok)
        if not m:
            return None
        key = _PREFIX[m.group(1)]
        val = m.group(2)
        params[key] = int(val) if key == "n" else float(val)
    return params or None


def read_ic_header(case_dir: str, ic_name: str):
    """Return (nmol, n, al) from a .ic header. NMOL/N/AL sit on lines 3/4/7
    (0-indexed 2/3/6) in both runt.ic and polyfj.ic."""
    path = os.path.join(case_dir, ic_name)
    if not os.path.isfile(path):
        return None
    lines = open(path).read().splitlines()
    try:
        return int(lines[2]), int(lines[3]), float(lines[6])
    except (IndexError, ValueError):
        return None


def load_observable(case_dir: str, ic_name: str, fname: str, col: int,
                    is_density: bool):
    """Return (r, y) for an observable column, or None. Density columns are
    normalised by the bulk number density rho0 = n_beads / AL^3 and clipped at
    r <= AL/2 (beyond the inscribed sphere the shell normalisation overcounts)."""
    path = os.path.join(case_dir, fname)
    hdr = read_ic_header(case_dir, ic_name)
    if not os.path.isfile(path) or hdr is None:
        return None
    nmol, n, al = hdr
    d = np.genfromtxt(path)
    if d.ndim != 2 or len(d) == 0 or d.shape[1] <= col:
        return None
    r, y = d[:, 1], d[:, col]
    mask = np.isfinite(y) & (r <= al / 2.0)
    if is_density:
        rho0 = (nmol * n) / al ** 3
        y = y / rho0
        mask &= y > 0
    return r[mask], y[mask]


def _coerce(val: str):
    try:
        return int(val)
    except ValueError:
        return float(val)


def overlay(runs_dir: str, ic_name: str, observables: dict, *,
            observable: str, vary: str, fix: dict, xmin: float, xmax: float,
            outdir: str) -> str:
    """Draw one overlay figure and return the output path.

    observables: name -> (filename, column, is_density, ylabel).
    vary: canonical parameter key to sweep (one curve per value).
    fix:  canonical key -> value constraints selecting the cases.
    """
    if observable not in observables:
        raise SystemExit(f"unknown observable '{observable}'; "
                         f"choose from {', '.join(observables)}")
    fname, col, is_density, ylabel = observables[observable]

    # collect matching cases
    series = []
    for name in sorted(os.listdir(runs_dir)):
        p = parse_params(name)
        if not p or vary not in p:
            continue
        if any(key not in p or abs(p[key] - val) > 1e-9 for key, val in fix.items()):
            continue
        res = load_observable(os.path.join(runs_dir, name), ic_name,
                              fname, col, is_density)
        if res is not None and len(res[0]):
            series.append((p[vary], res))
    if not series:
        raise SystemExit(
            f"no cases match observable={observable}, vary={vary}, fix={fix} "
            f"(have you run the simulations in {runs_dir}?)")

    series.sort(key=lambda s: s[0])
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for vval, (r, y) in series:
        ax.plot(r, y, lw=1.2, label=f"{vary}={vval:g}")
    fixed_txt = ", ".join(f"{k}={v:g}" for k, v in sorted(fix.items()))
    ax.set_xlabel(r"$r/\sigma$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(xmin, xmax)
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False, fontsize=9, title=vary)
    ax.set_title(f"{observable} vs r" + (f"  ({fixed_txt})" if fixed_txt else ""))
    fig.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    tag = "_".join(f"{k}{v:g}" for k, v in sorted(fix.items()))
    stem = f"{observable}_vary-{vary}" + (f"_{tag}" if tag else "")
    out = os.path.join(outdir, stem + ".png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def add_common_args(ap, observables: dict, default_observable: str = "density"):
    """Attach the standard --observable/--vary/--fix/--xmin/--xmax/--outdir args."""
    ap.add_argument("--observable", default=default_observable,
                    choices=list(observables),
                    help="quantity to plot on the y-axis")
    ap.add_argument("--vary", required=True,
                    choices=["n", "sphere", "eta", "eps_ff", "eps_sf"],
                    help="parameter swept within the figure (one curve per value)")
    ap.add_argument("--fix", action="append", default=[], metavar="KEY=VAL",
                    help="pin a parameter (repeatable), e.g. --fix n=8 --fix eta=0.2")
    ap.add_argument("--xmin", type=float, default=1.0)
    ap.add_argument("--xmax", type=float, default=7.0)


def parse_fix(items: list[str]) -> dict:
    fix = {}
    for it in items:
        if "=" not in it:
            raise SystemExit(f"--fix expects KEY=VAL, got '{it}'")
        k, v = it.split("=", 1)
        fix[k] = _coerce(v)
    return fix
