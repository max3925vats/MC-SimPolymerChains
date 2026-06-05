#!/usr/bin/env bash
#
# EXAMPLE (not part of the general tooling): reproduce the athermal analyses
# from the accompanying SRFP-2016 report using the general polyfj tools. This is
# just one particular choice of parameter sweeps; the underlying scripts accept
# any parameters (see tools/setup_runs.py --help and tools/plot_results.py --help).
#
# Run from the polyfj/ directory:  ./examples/reproduce_report.sh
#
# For converged (report-quality) curves, append e.g. --ncon 1000000000 to the
# setup commands below; the defaults are a fast smoke setting.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

# --- generate the decks for each sweep (overlapping cases are written once) ---
# effect of central sphere size (8-mers, two packings, three sphere radii)
python3 tools/setup_runs.py --n 8 --sphere 1.0 1.5 2.0 --eta 0.2 0.4
# effect of chain length (fixed packing, two sphere radii, three lengths)
python3 tools/setup_runs.py --n 4 8 16 --sphere 1.5 2.0 --eta 0.2
# effect of packing fraction (8-mers, one sphere radius, three packings)
python3 tools/setup_runs.py --n 8 --sphere 1.5 --eta 0.1 0.2 0.4

# --- run every generated case ---
./run_all.sh

# --- draw the overlays ---
# sphere-size effect, one figure per packing
python3 tools/plot_results.py --vary sphere --fix n=8 --fix eta=0.2
python3 tools/plot_results.py --vary sphere --fix n=8 --fix eta=0.4
# chain-length effect, one figure per sphere radius
python3 tools/plot_results.py --vary n --fix sphere=1.5 --fix eta=0.2
python3 tools/plot_results.py --vary n --fix sphere=2.0 --fix eta=0.2
# packing-fraction effect
python3 tools/plot_results.py --vary eta --fix n=8 --fix sphere=1.5
# structural profiles (chain-length effect on Rg^2 and segmental order)
python3 tools/plot_results.py --observable rg2 --vary n --fix sphere=1.5 --fix eta=0.2
python3 tools/plot_results.py --observable seg --vary n --fix sphere=1.5 --fix eta=0.2

echo "Figures written to polyfj/figures/"
