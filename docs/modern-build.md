# Modern Fortran Build — polymc

This document covers building and running the modern (fpm) rewrite alongside the
legacy Fortran source.  The legacy path is the default everywhere; the modern
path is opt-in via `POLYMC_ENGINE=modern`.

## Prerequisites

| Tool | Install (macOS) |
|------|----------------|
| gfortran | `brew install gcc` |
| fpm (Fortran Package Manager) | `brew install fpm` |

Verify: `gfortran --version` and `fpm --version`.

## Build

```bash
# from repo root
fpm build          # compiles src/ and app/ targets
fpm test           # runs the test suite under test/
```

Makefile shortcut (wraps fpm):

```bash
make build         # equivalent to fpm build
```

Built binaries are placed under `build/<compiler-hash>/app/runt` and
`build/<compiler-hash>/app/polyfj`.

## Running via the harness

Both `runt/run_all.sh` and `polyfj/run_all.sh` accept an environment variable
`POLYMC_ENGINE` (values: `legacy` [default] or `modern`).

### runt

```bash
# Default — compiles legacy runt.f with gfortran and runs it (unchanged behaviour)
cd runt
./run_all.sh

# Opt-in modern engine — fpm build, then runs the modern binary
POLYMC_ENGINE=modern ./run_all.sh

# With case-name filters (works in both modes)
POLYMC_ENGINE=modern ./run_all.sh Fig4 n8
```

### polyfj

```bash
cd polyfj

# Legacy (default)
./run_all.sh

# Modern
POLYMC_ENGINE=modern ./run_all.sh

# With filters
POLYMC_ENGINE=modern ./run_all.sh Fig2 n16
```

When `POLYMC_ENGINE=modern` the script runs `fpm build` from the repo root,
locates the binary under `build/`, and drops it into the same per-case loop.
The output files (`runt.out`, `denrunt.out`, `runt.fc` / `polyfj.out`, etc.)
are identical in layout to the legacy versions.

## Running a single case manually

```bash
# Build once
fpm build

# Locate the binary
BIN=$(find build -type f -path '*/app/runt' | head -1)

# Run in-place in a case directory
cd runt/runs/n4_s1.0_e0.1_ff0_sf0
"$BIN"
```

The modern binary reads `runt.inp` and `runt.ic` from the current working
directory, just like the legacy binary.

## runt energy modes

Line count in `runt.inp` controls the energy calculation:

| Lines | Mode | When to use |
|-------|------|-------------|
| 10 (default) | Exact | Faithful to legacy; use for production / validation |
| 11 (rcut appended) | Cutoff | Approximation — faster; use for performance benchmarking |

Generate exact-mode decks (default):

```bash
python3 tools/setup_runs.py
```

Generate cutoff-mode decks with rcut = 5.0 σ:

```bash
python3 tools/setup_runs.py --rcut 5.0 --outdir /tmp/runt_cutoff_runs
```

Or append the cutoff radius manually as an 11th line to any existing `runt.inp`.

Note: the cutoff path is only implemented in the modern binary.  Running a
cutoff deck with the legacy binary will silently ignore the 11th line (legacy
Fortran stops reading after line 10).

## polyfj observable styles

The modern polyfj binary accepts an optional flag controlling how density
profiles are normalised:

| Flag | Mode |
|------|------|
| (none — default) | Corrected observables |
| `--legacy-obs` | Legacy-style observables (for direct comparison with old output) |

This flag is passed directly to the binary, not through the harness.

## Validation

See `docs/validation.md` for the legacy-vs-modern regression test protocol and
acceptance criteria.
