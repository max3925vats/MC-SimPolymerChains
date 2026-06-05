#!/usr/bin/env python3
"""
compare_profiles.py — Density-profile comparator for runt/polyfj MC runs.

Compares two denrunt.out (or polyden.out) files produced by the modern fpm
runt and the legacy gfortran runt.  The files have 3 numeric columns:

    idx   dist   density

where  dist  is in real sigma units and  density  is number density.

Usage
-----
python3 tools/compare_profiles.py FILE_A FILE_B [OPTIONS]

Options
-------
--naver N        Number of accumulation samples (ncon/nskip).  When given,
                 per-bin Poisson error is estimated and the comparison uses
                 sigma-based flagging.  When omitted, relative-difference
                 flagging is used instead.
--tol-sigmas F   Flag bins where |dA-dB| > F * combined_sigma  [default 3.0]
--rtol F         Relative-tolerance for no-naver mode  [default 0.05]
--rmin F         Skip bins with dist < rmin  (sigma units)
--rmax F         Skip bins with dist > rmax  (sigma units)
--selftest       Run built-in PASS+FAIL self-tests and exit (no file args needed)
"""

import argparse
import math
import os
import sys
import tempfile

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_profile(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load a 3-column density profile; skip blank / short lines robustly.

    Returns (idx, dist, density) as 1-D float arrays (idx cast to int).
    """
    rows = []
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, 1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 3:
                # Too few fields — skip (could be a header in some variants)
                continue
            try:
                row = [float(p) for p in parts[:3]]
            except ValueError:
                # Non-numeric first token → skip
                continue
            rows.append(row)

    if not rows:
        raise ValueError(f"No valid data rows found in '{path}'")

    data = np.array(rows, dtype=float)
    idx  = data[:, 0].astype(int)
    dist = data[:, 1]
    dens = data[:, 2]
    return idx, dist, dens


def align_profiles(
    distA: np.ndarray,
    densA: np.ndarray,
    distB: np.ndarray,
    densB: np.ndarray,
    atol: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (dist, dA, dB) restricted to bins present in both files.

    If the grids are identical (same length + all |dA_dist - dB_dist| < atol)
    no reindexing is needed.  Otherwise the intersection on dist is used and
    a warning is printed.
    """
    if len(distA) == len(distB) and np.all(np.abs(distA - distB) < atol):
        return distA.copy(), densA.copy(), densB.copy()

    # Intersect on dist values rounded to grid precision
    keysA = np.round(distA / atol).astype(int)
    keysB = np.round(distB / atol).astype(int)
    common = np.array(sorted(set(keysA) & set(keysB)))

    if len(common) == 0:
        raise ValueError(
            "Files have no overlapping dist values (within 1e-6). "
            "Check that both runs used the same bin size (bsz)."
        )

    # Build index maps
    mapA = {k: i for i, k in enumerate(keysA)}
    mapB = {k: i for i, k in enumerate(keysB)}

    idxA = np.array([mapA[k] for k in common])
    idxB = np.array([mapB[k] for k in common])

    print(
        f"  [align] Grids differ — found {len(common)} common bins "
        f"(file_A had {len(distA)}, file_B had {len(distB)}).",
        file=sys.stderr,
    )

    return distA[idxA], densA[idxA], densB[idxB]


# ---------------------------------------------------------------------------
# Core comparison logic
# ---------------------------------------------------------------------------

def compare(
    file_a: str,
    file_b: str,
    naver: int | None = None,
    tol_sigmas: float = 3.0,
    rtol: float = 0.05,
    rmin: float | None = None,
    rmax: float | None = None,
) -> bool:
    """Compare two density profiles.

    Returns True if PASS (no flagged bins), False if FAIL.
    """
    # --- Load ---
    _, distA, densA = load_profile(file_a)
    _, distB, densB = load_profile(file_b)

    # --- Align on common dist grid ---
    dist, densA, densB = align_profiles(distA, densA, distB, densB)

    # --- Spatial filter ---
    mask = np.ones(len(dist), dtype=bool)
    if rmin is not None:
        mask &= dist >= rmin
    if rmax is not None:
        mask &= dist <= rmax

    # --- Floor filter: skip bins where both profiles are essentially zero ---
    # Uses a small absolute floor (1e-9 in sigma^{-3}) to avoid 0/0 noise
    floor = 1e-9
    mask &= (densA > floor) | (densB > floor)

    dist_c  = dist[mask]
    densA_c = densA[mask]
    densB_c = densB[mask]
    n_bins  = len(dist_c)

    if n_bins == 0:
        print("  [warn] No bins remain after filtering — nothing to compare.")
        print("PASS")
        return True

    # --- Per-bin difference ---
    diff = np.abs(densA_c - densB_c)

    if naver is not None:
        # Poisson error estimate -------------------------------------------
        # Shell volumes: bin spacing from dist grid (assume uniform)
        bsz = dist_c[1] - dist_c[0] if len(dist_c) > 1 else (dist_c[0] * 0.1)
        rl = dist_c - bsz / 2.0
        ru = dist_c + bsz / 2.0
        # Clip to non-negative radii
        rl = np.clip(rl, 0.0, None)
        shell_vol = (4.0 / 3.0) * math.pi * (ru**3 - rl**3)

        # Approximate counts: density * shell_volume * naver
        countA = densA_c * shell_vol * naver
        countB = densB_c * shell_vol * naver

        # sigma_density = density / sqrt(count) = 1/(sqrt(count)*naver*vol_bin)
        # but equivalently: sigma_dens = sqrt(density / (shell_vol * naver))
        sigA = np.sqrt(densA_c / (shell_vol * naver + 1e-300))
        sigB = np.sqrt(densB_c / (shell_vol * naver + 1e-300))
        combined_sig = np.sqrt(sigA**2 + sigB**2)

        metric = diff / np.where(combined_sig > 0, combined_sig, np.inf)
        threshold = tol_sigmas

    else:
        # Relative difference ---------------------------------------------
        mean_dens = 0.5 * (densA_c + densB_c)
        rel_diff  = diff / np.where(mean_dens > 0, mean_dens, np.inf)
        metric    = rel_diff
        combined_sig = None   # not used in rtol mode
        countA = countB = None
        threshold = rtol

    # --- Flagging ---
    flagged   = metric > threshold
    n_flagged = int(np.sum(flagged))

    # Worst bin
    worst = int(np.argmax(metric))
    worst_dist  = dist_c[worst]
    worst_dA    = densA_c[worst]
    worst_dB    = densB_c[worst]
    worst_diff  = diff[worst]
    worst_metric = metric[worst]

    # --- Report ---
    print(f"  File A : {file_a}")
    print(f"  File B : {file_b}")
    print(f"  Bins compared  : {n_bins}")
    if naver is not None:
        print(f"  Mode   : Poisson sigma  (naver={naver}, tol={tol_sigmas} sigma)")
    else:
        print(f"  Mode   : relative diff  (rtol={rtol})")

    if rmin is not None or rmax is not None:
        lo = f"{rmin:.3f}" if rmin is not None else "-inf"
        hi = f"{rmax:.3f}" if rmax is not None else "+inf"
        print(f"  Range  : dist in [{lo}, {hi}] sigma")

    print()
    print(f"  Worst bin:")
    print(f"    dist   = {worst_dist:.4f} sigma")
    print(f"    dens_A = {worst_dA:.6e}")
    print(f"    dens_B = {worst_dB:.6e}")
    print(f"    |diff| = {worst_diff:.6e}")
    if naver is not None:
        worst_sig = float(np.sqrt(sigA[worst]**2 + sigB[worst]**2))
        worst_count = 0.5 * (float(countA[worst]) + float(countB[worst]))
        print(f"    sigma  = {worst_sig:.6e}  =>  {worst_metric:.2f} sigma")
        print(f"    ~count = {worst_count:.1f}  (approx mean bin count)")
    else:
        print(f"    rel    = {worst_metric:.4f}")

    print()
    print(f"  Flagged bins   : {n_flagged} / {n_bins}")
    print()

    passed = n_flagged == 0
    status = "PASS" if passed else "FAIL"
    print(status)
    return passed


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _write_profile(path: str, dist: np.ndarray, dens: np.ndarray) -> None:
    """Write a synthetic 3-column profile (idx dist density)."""
    with open(path, "w") as fh:
        for i, (d, rho) in enumerate(zip(dist, dens), 1):
            fh.write(f"  {i:4d}  {d:8.4f}   {rho:.6e}\n")


def run_selftest() -> None:
    """Synthetic self-test: verify PASS and FAIL paths.

    Writes temporary profile files, calls compare(), asserts exit codes.
    Prints the full comparator output for both cases.
    """
    import subprocess

    # Synthetic profile parameters
    bsz  = 0.05
    dist = np.arange(bsz / 2, 8.0, bsz)  # ~160 bins
    # Smooth Gaussian-like density peak at r=2.5
    dens_base = 0.5 * np.exp(-((dist - 2.5) ** 2) / (2 * 0.4**2)) + 0.05

    script = os.path.abspath(__file__)

    with tempfile.TemporaryDirectory() as tmpdir:
        fa = os.path.join(tmpdir, "profile_a.out")
        fb_pass = os.path.join(tmpdir, "profile_b_pass.out")
        fb_fail = os.path.join(tmpdir, "profile_b_fail.out")

        # PASS case: identical profiles
        _write_profile(fa, dist, dens_base)
        _write_profile(fb_pass, dist, dens_base)

        # FAIL case: profile with one bin perturbed by a huge amount
        dens_fail = dens_base.copy()
        spike_idx = int(len(dist) * 0.4)  # bin near r=3.2
        dens_fail[spike_idx] *= 5.0       # 5× spike — far outside any tolerance

        _write_profile(fb_fail, dist, dens_fail)

        # --- PASS path ---
        print("=" * 60)
        print("SELF-TEST 1: identical profiles  →  expect PASS")
        print("=" * 60)
        naver = 50000
        result_pass = subprocess.run(
            [
                sys.executable, script,
                fa, fb_pass,
                "--naver", str(naver),
                "--tol-sigmas", "3.0",
                "--rmin", "1.5",
            ],
            capture_output=False,
        )
        print()
        assert result_pass.returncode == 0, (
            f"SELF-TEST 1 FAILED: expected exit 0 (PASS), got {result_pass.returncode}"
        )
        print(f"[selftest] exit code = {result_pass.returncode}  ✓ PASS as expected\n")

        # --- FAIL path ---
        print("=" * 60)
        print("SELF-TEST 2: profiles with large spike in one bin  →  expect FAIL")
        print("=" * 60)
        result_fail = subprocess.run(
            [
                sys.executable, script,
                fa, fb_fail,
                "--naver", str(naver),
                "--tol-sigmas", "3.0",
                "--rmin", "1.5",
            ],
            capture_output=False,
        )
        print()
        assert result_fail.returncode == 1, (
            f"SELF-TEST 2 FAILED: expected exit 1 (FAIL), got {result_fail.returncode}"
        )
        print(f"[selftest] exit code = {result_fail.returncode}  ✓ FAIL as expected\n")

        # --- rtol path (no --naver) ---
        print("=" * 60)
        print("SELF-TEST 3: identical profiles, rtol mode (no --naver)  →  expect PASS")
        print("=" * 60)
        result_rtol = subprocess.run(
            [
                sys.executable, script,
                fa, fb_pass,
                "--rtol", "0.05",
                "--rmin", "1.5",
            ],
            capture_output=False,
        )
        print()
        assert result_rtol.returncode == 0, (
            f"SELF-TEST 3 FAILED: expected exit 0 (PASS), got {result_rtol.returncode}"
        )
        print(f"[selftest] exit code = {result_rtol.returncode}  ✓ PASS as expected\n")

    print("=" * 60)
    print("All self-tests passed.")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compare two density-profile files from runt/polyfj MC runs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("file_a", nargs="?", help="Path to profile A (denrunt.out)")
    p.add_argument("file_b", nargs="?", help="Path to profile B (denrunt.out)")
    p.add_argument(
        "--naver", type=int, default=None,
        help="Accumulation samples (ncon/nskip). Enables Poisson error estimate.",
    )
    p.add_argument(
        "--tol-sigmas", type=float, default=3.0, dest="tol_sigmas",
        help="Flag bins where |diff| > tol_sigmas * combined_sigma  [default 3.0]",
    )
    p.add_argument(
        "--rtol", type=float, default=0.05,
        help="Relative tolerance for no-naver mode  [default 0.05]",
    )
    p.add_argument("--rmin", type=float, default=None, help="Min dist (sigma) to include")
    p.add_argument("--rmax", type=float, default=None, help="Max dist (sigma) to include")
    p.add_argument(
        "--selftest", action="store_true",
        help="Run built-in PASS+FAIL self-tests and exit (no file args needed).",
    )
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.selftest:
        run_selftest()
        sys.exit(0)

    if not args.file_a or not args.file_b:
        parser.error("Two positional file arguments are required unless --selftest is used.")

    passed = compare(
        file_a=args.file_a,
        file_b=args.file_b,
        naver=args.naver,
        tol_sigmas=args.tol_sigmas,
        rtol=args.rtol,
        rmin=args.rmin,
        rmax=args.rmax,
    )
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
