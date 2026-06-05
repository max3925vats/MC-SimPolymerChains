#!/usr/bin/env bash
#
# Compile the (RNG-fixed) polyfj, then run every generated athermal case under
# runs/. Each case has its own polyfj.inp / polyfj.ic and writes its outputs
# (polyfj.out, polyden.out, seg.out, rgfj.out, g2fj.out, polyfj.fc) into its
# own directory.
#
# Usage:
#   ./run_all.sh                          # run every case (legacy engine, default)
#   ./run_all.sh Fig2 n16                # run only cases matching a filter
#   POLYMC_ENGINE=modern ./run_all.sh    # use fpm-built modern binary instead
#
# Notes:
#   * Generate the cases first:  python3 tools/setup_runs.py
#   * polyfj seeds its RNG from the wall clock; runs are not bit-reproducible.
#   * Runtime scales with NCON (set in polyfj.inp). The smoke default (2e6) is
#     minutes per case; the report's 1e9 is hours -- regenerate with
#     `python3 tools/setup_runs.py --ncon 1000000000` when ready.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_DIR="$ROOT/runs"

# Engine selection: default is legacy (compile + run local gfortran binary).
# Set POLYMC_ENGINE=modern to use the fpm-built binary from the repo build/ dir.
ENGINE="${POLYMC_ENGINE:-legacy}"

if [[ ! -d "$RUNS_DIR" ]]; then
  echo "No runs/ directory. Generate cases first:  python3 tools/setup_runs.py" >&2
  exit 1
fi

if [[ "$ENGINE" == "modern" ]]; then
  # Repo root is one level above this script's directory (polyfj/).
  REPO_ROOT="$(cd "$ROOT/.." && pwd)"
  echo ">> Engine: modern  (fpm build)"
  echo ">> Building via fpm ..."
  ( cd "$REPO_ROOT" && fpm build )
  BIN="$(find "$REPO_ROOT/build" -type f -path '*/app/polyfj' | head -1)"
  if [[ -z "$BIN" ]]; then
    echo "!! Could not locate modern polyfj binary under $REPO_ROOT/build" >&2
    exit 1
  fi
  echo ">> Modern binary: $BIN"
else
  # Legacy path: compile from source, unchanged behaviour.
  BIN="$ROOT/polyfj_bin"
  SRC=("$ROOT/polyfj.f" "$ROOT/polyfj-moves.f")
  echo ">> Engine: legacy"
  echo ">> Compiling polyfj ..."
  gfortran -std=legacy -O2 -o "$BIN" "${SRC[@]}"
  echo ">> Built $BIN"
fi

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
  if [[ ! -f "$d/polyfj.inp" || ! -f "$d/polyfj.ic" ]]; then
    echo "!! $name: missing polyfj.inp/polyfj.ic, skipping" >&2
    continue
  fi
  echo ">> Running $name ..."
  start=$(date +%s)
  if ( cd "$d" && "$BIN" ) >"$d/run.log" 2>&1; then
    elapsed=$(( $(date +%s) - start ))
    echo "   done ($name) in ${elapsed}s -> $d/polyfj.out"
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
