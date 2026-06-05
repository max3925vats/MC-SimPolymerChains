#!/usr/bin/env bash
#
# Compile the (BBATT-fixed, CVERT-added) runt once, then run every generated
# case under runs/. Each case has its own runt.inp / runt.ic and writes its
# outputs (runt.out, denrunt.out, runt.fc) into its own directory.
#
# Usage:
#   ./run_all.sh                 # run every case in runs/
#   ./run_all.sh Fig4 t_n8       # run only cases whose dir name matches a filter
#
# Notes:
#   * Generate the cases first:  python3 tools/setup_runs.py
#   * runt seeds its RNG from the wall clock, so runs are not bit-reproducible;
#     the report block-averaged 4 runs per condition. Re-run as needed.
#   * Runtime scales with NCON (set in runt.inp). The smoke default (2e6) is
#     minutes per case; the report's 1e9 is hours -- regenerate with
#     `python3 tools/setup_runs.py --ncon 1000000000` when ready.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_DIR="$ROOT/runs"
BIN="$ROOT/runt_bin"
SRC=("$ROOT/runt.f" "$ROOT/runt-moves.f")

if [[ ! -d "$RUNS_DIR" ]]; then
  echo "No runs/ directory. Generate cases first:  python3 tools/setup_runs.py" >&2
  exit 1
fi

echo ">> Compiling runt ..."
gfortran -std=legacy -O2 -o "$BIN" "${SRC[@]}"
echo ">> Built $BIN"

# Optional name filters from the command line.
filters=("$@")
match() {
  [[ ${#filters[@]} -eq 0 ]] && return 0
  local name="$1"
  for f in "${filters[@]}"; do [[ "$name" == *"$f"* ]] && return 0; done
  return 1
}

n_ok=0; n_fail=0; failed=()
for d in "$RUNS_DIR"/*/; do
  name="$(basename "$d")"
  match "$name" || continue
  if [[ ! -f "$d/runt.inp" || ! -f "$d/runt.ic" ]]; then
    echo "!! $name: missing runt.inp/runt.ic, skipping" >&2
    continue
  fi
  echo ">> Running $name ..."
  start=$(date +%s)
  if ( cd "$d" && "$BIN" ) >"$d/run.log" 2>&1; then
    elapsed=$(( $(date +%s) - start ))
    echo "   done ($name) in ${elapsed}s -> $d/runt.out"
    n_ok=$((n_ok + 1))
  else
    echo "!! FAILED: $name (see $d/run.log)" >&2
    n_fail=$((n_fail + 1)); failed+=("$name")
  fi
done

echo
echo ">> Summary: $n_ok succeeded, $n_fail failed"
if [[ $n_fail -gt 0 ]]; then
  printf '   failed: %s\n' "${failed[@]}"
  exit 1
fi
