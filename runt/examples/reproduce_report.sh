#!/usr/bin/env bash
#
# EXAMPLE (not part of the general tooling): reproduce the attraction analysis
# from the accompanying SRFP-2016 report using the general runt tools. This is
# just one particular choice of parameter sweeps; the underlying scripts accept
# any parameters (see tools/setup_runs.py --help and tools/plot_results.py --help).
#
# Run from the runt/ directory:  ./examples/reproduce_report.sh
#
# For converged (report-quality) curves, append e.g. --ncon 1000000000 to the
# setup command below; the default is a fast smoke setting.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

# --- generate the decks: 8-mers, two sphere radii, fixed packing, the 2x2
#     grid of chain-chain and sphere-chain attraction strengths {0, 0.5} ---
python3 tools/setup_runs.py --n 8 --sphere 1.5 2.0 --eta 0.2 \
    --eps-ff 0 0.5 --eps-sf 0 0.5

# --- run every generated case ---
./run_all.sh

# --- draw the overlays: attraction effect, one figure per sphere radius ---
# vary the sphere-chain attraction (chain-chain off)
python3 tools/plot_results.py --vary eps_sf --fix n=8 --fix sphere=1.5 --fix eta=0.2 --fix eps_ff=0
python3 tools/plot_results.py --vary eps_sf --fix n=8 --fix sphere=2.0 --fix eta=0.2 --fix eps_ff=0
# vary the chain-chain attraction (sphere-chain off)
python3 tools/plot_results.py --vary eps_ff --fix n=8 --fix sphere=1.5 --fix eta=0.2 --fix eps_sf=0
python3 tools/plot_results.py --vary eps_ff --fix n=8 --fix sphere=2.0 --fix eta=0.2 --fix eps_sf=0

echo "Figures written to runt/figures/"
