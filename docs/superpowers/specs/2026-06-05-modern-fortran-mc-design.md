# Modern Fortran Rewrite of the Polymer-Chain / Central-Sphere MC Codes

**Date:** 2026-06-05
**Status:** Approved design, pending implementation plan
**Scope:** Rewrite the `runt` (hard-core + Yukawa, thermal) and `polyfj` (hard-chain,
athermal) Monte Carlo programs as a modern Fortran codebase: two executables over a
shared module library, with a real cell list, a modern RNG, and a test + regression
harness. The legacy `.f` code and the existing Python tooling are preserved.

---

## 1. Goals

Driven by three priorities the user selected:

- **Performance & scale** — implement the neighbour cell list that the legacy code
  declared but never used, so per-move cost drops from O(N_beads) to O(local).
- **Extend the science** — a clean, modular base for new potentials, observables,
  geometries, and move types.
- **Maintainability** — `implicit none`, derived types, `allocatable` arrays,
  `intent`, and unit tests, replacing `COMMON` blocks and fixed-size arrays.

Non-goals for v1: GPU, OpenMP/threaded kernels (conditions already run as separate
processes), changing the physics model, or replacing the Python analysis tooling.

## 2. Architecture

A new **fpm project at the repository root**, alongside the untouched legacy code:

```
fpm.toml              # primary build (Fortran Package Manager)
Makefile              # thin gfortran fallback build
src/                  # shared module library (the kernel)
app/runt.f90          # runt driver (thin: orchestration only)
app/polyfj.f90        # polyfj driver (thin)
test/                 # test-drive unit tests
```

Legacy `runt/` and `polyfj/` (the `.f` sources, decks, and Python tooling) remain in
the tree as the validation reference until the modern version passes regression.

### 2.1 Shared core modules (`src/`)

| Module | Responsibility |
|---|---|
| `polymc_kinds` | `dp` precision parameter |
| `polymc_rng` | xoshiro256** generator (splitmix64 seeding), `uniform()`, `random_unit_vector()` |
| `polymc_box` | `box_t` type; minimum-image, reduced<->real scaling |
| `polymc_cell_list` | cell-list build/update + neighbour iteration (the corrected `OVER1`) |
| `polymc_overlap` | hard-core overlap test (sphere + inter/intra) via the cell list |
| `polymc_chain` | chain helpers (in-place reversal a la `CVERT`, bond ops) |
| `polymc_config` | `params_t`/`system_t` types; read `.inp`/`.ic`, write `.fc`/`.out` |
| `polymc_binning` | radial histogram / profile accumulation |
| `polymc_linalg` | symmetric 3x3 eigensolver (Jacobi) for the shape tensor |
| `polymc_io` | formatted output helpers |

### 2.2 Per-program modules

- **runt**: `runt_potential` (Yukawa pair + sphere, unit-correct), `runt_energy`
  (total + incremental), `runt_moves` (Dickman / reptation / CCB with Metropolis
  acceptance), density observable. The energy supports **two selectable modes**:
  - **exact** (default) — full inter-molecular sum, no cutoff, faithful to legacy
    `NEWEN1`/`OLDEN1`; O(N_beads) per move.
  - **cutoff** (opt-in via input `rcut > 0`) — Yukawa truncated at `rcut`, summed via
    the cell list; O(local) per move. This is a deliberate approximation, provided so
    the user can benchmark it against the exact result.
- **polyfj**: `polyfj_moves` (hard-core accept-on-no-overlap), `polyfj_observables`
  (total/end/mid density, Rg & Re tensors, segmental order parameter, shape semi-axes
  via `polymc_linalg`).

## 3. Key design decisions

- **RNG upgrade (behavioural change).** Replace Park-Miller (period ~2.1e9, too short
  for 1e9 moves x several draws each -> sequence reuse and correlations) with
  xoshiro256** (period 2^256). Deterministically seedable, so runs become reproducible
  (the legacy time-seeded code was not). Seed comes from the input (with a documented
  default); a fixed seed gives bit-identical reruns of the *modern* code.
- **Cell list accelerates hard-core overlap (polyfj) and the optional cutoff energy
  path (runt).** Hard-core overlap queries go through `polymc_cell_list`; correctness is
  guaranteed by a unit test asserting cell-list overlap equals brute-force O(N^2)
  overlap on random configurations.
- **Exact by default; approximation is opt-in.** `runt`'s energy defaults to the exact
  full sum (no approximation beyond the original), so `runt` retains the legacy
  O(N_beads) per-move cost. A finite `rcut` enables the truncated cell-list path
  (O(local)) purely so the user can quantify the approximation against the exact run.
  `polyfj` (athermal) always gets the cell-list overlap speedup.
- **I/O format compatibility.** Read the existing `runt.inp`/`runt.ic`/`polyfj.inp`/
  `polyfj.ic` and write the existing `*.out`/`*.fc` outputs byte-compatibly, so the
  Python `setup_runs`/`run_all`/plotters and all generated decks work unchanged, and
  legacy-vs-modern regression is direct.
- **State via derived types**, not `COMMON`. `allocatable` arrays sized from input;
  the fixed `NBMAX`/`NMAX` ceilings are removed.

## 4. Validation strategy (two tiers)

1. **Unit tests (deterministic, exact)** under `test/`:
   - PBC minimum-image correctness.
   - `random_unit_vector` distribution (unit length; mean ~ 0 over many draws).
   - hard-core overlap detection on hand-built configs.
   - **cell-list overlap equals brute-force** on random configs (the trust anchor for
     the optimization). `runt`'s energy is exact-by-construction (full sum), checked by
     an incremental-equals-full-recompute test rather than against a cutoff.
   - Yukawa potential values at known separations.
   - Jacobi eigensolver vs analytically known symmetric matrices.
2. **Statistical regression against the legacy binaries** (the user's explicit
   requirement): run the existing committed decks through both the legacy `.f` build
   and the modern build (in **exact** mode) at matched conditions; the density profiles
   (and, for polyfj, the structural profiles) must agree within combined Monte-Carlo
   error. Because the RNGs differ, agreement is statistical, not bit-exact. A small
   comparison script loads both `denrunt.out`/`polyden.out` sets and checks per-bin
   agreement within a tolerance scaled by the sampling error.
3. **Cutoff-vs-exact comparison** (runt): run the same deck in exact and cutoff modes
   and report the per-bin density difference, so the size of the cutoff approximation is
   measured, not assumed.

## 5. Development methodology

**Subagent-driven development with two-stage review** (user's explicit choice):

- The implementation plan (next step, via the writing-plans skill) decomposes the work
  into independent, ordered tasks (e.g. one module or module-group per task), each with
  its own tests and acceptance criteria.
- Each task is implemented by a dedicated sub-agent against the task brief.
- **Two-stage review per task:** (1) an automated/code-review pass for correctness,
  style, and contract adherence; (2) an independent verification pass that runs the
  task's tests (and, where relevant, the cell-list-vs-brute-force and legacy-regression
  checks) before the task is accepted.
- Work proceeds on the `modern-fortran-rewrite` branch; the legacy code stays intact so
  every stage can be validated against it.

## 6. Phasing

1. **Foundation** — `polymc_kinds/rng/box/cell_list/overlap` + unit tests, including
   cell-list-vs-brute-force. Establish fpm + Makefile builds and the test target.
2. **runt** on the new core — potentials, energy, moves, density; validate vs legacy
   `runt` (regression on existing decks).
3. **polyfj** on the new core — moves + structural observables + Jacobi; validate vs
   legacy `polyfj`.
4. **Cutover & polish** — point `run_all.sh` at the new binaries (keeping the legacy
   path available), finalize docs. Legacy `.f` retained in-tree.

## 7. Risks / open points

- gfortran is the assumed compiler (the legacy build already requires it). fpm install
  is via Homebrew.
- Statistical regression needs runs long enough to shrink MC error below the comparison
  tolerance; short smoke runs will be too noisy to validate against. The plan will
  specify the regression run length.
- Cutover of `run_all.sh` is deferred until phase 4 so the legacy pipeline keeps working
  throughout.
