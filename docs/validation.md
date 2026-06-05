# Validating the modern runt/polyfj against legacy

This runbook documents the exact commands to verify that the modern fpm-built
`runt` binary produces density profiles statistically consistent with the
legacy gfortran `runt`.  Run these after any significant change to the modern
code, or when benchmarking a new deck.

> **Do not run this document's commands inside the repo working tree.**
> All simulation output goes to `/tmp/val_*` to keep the tree clean.

---

## 1. Build both binaries

**Legacy** (single-file gfortran compile, fixed standard):

```bash
gfortran -std=legacy -O2 \
  -o /tmp/runt_legacy \
  runt/runt.f runt/runt-moves.f
```

**Modern** (fpm build from repo root):

```bash
fpm build
MODERN_RUNT=$(find build -type f -path '*/app/runt' | head -1)
echo "Modern binary: $MODERN_RUNT"
```

---

## 2. Set up run directories

Choose the deck with attraction on (ff0 = athermal sphere, sf0.5 = attractive
chain segments, n=8, sigma=1.5):

```bash
DECK=runt/runs/n8_s1.5_e0.2_ff0_sf0.5
```

Create two run directories and copy the deck files into each:

```bash
mkdir -p /tmp/val_legacy /tmp/val_modern_exact
cp "$DECK/runt.ic"  /tmp/val_legacy/
cp "$DECK/runt.ic"  /tmp/val_modern_exact/
cp "$DECK/runt.inp" /tmp/val_legacy/
cp "$DECK/runt.inp" /tmp/val_modern_exact/
```

Edit `runt.inp` in both directories for a production-length run.
The deck's `runt.inp` already has good move-type fractions; only `ncon` and
`nskip` need to change.  The file has one line per parameter:

- Line 1 = `ncon` (total MC moves)
- Line 2 = `nskip` (moves per accumulation)
- Line 11 = `rcut` (optional cutoff; blank/absent = exact mode)

Set `ncon=200000000` and `nskip=1000`.  Exact mode requires the deck to have
**no 11th line** (rcut absent).  The committed decks have exactly 10 lines, so
they already satisfy this — no edit is needed for exact mode beyond lines 1-2:

```bash
# In-place edit with sed (lines 1 and 2 only, leave all other lines untouched).
# Exact mode = no rcut line; the 10-line deck already has none, so we do not
# add or blank an 11th line here.
for DIR in /tmp/val_legacy /tmp/val_modern_exact; do
    sed -i '' '1s/.*/200000000/' "$DIR/runt.inp"
    sed -i '' '2s/.*/1000/'      "$DIR/runt.inp"
done
```

Verify the edits are correct:

```bash
echo "=== legacy ===" && cat -n /tmp/val_legacy/runt.inp
echo "=== modern ===" && cat -n /tmp/val_modern_exact/runt.inp
```

---

## 3. Run both simulations

These take on the order of minutes to hours depending on hardware.
Run them sequentially (or in separate terminals in parallel):

**Legacy** (wall-clock seed; non-reproducible between runs):

```bash
( cd /tmp/val_legacy && /tmp/runt_legacy )
```

**Modern exact** (fixed seed 12345; fully reproducible):

```bash
( cd /tmp/val_modern_exact && "$MODERN_RUNT" )
```

Both write `denrunt.out` with 3 columns: `idx  dist  density`.

---

## 4. Compare modern-exact vs legacy

`naver = ncon / nskip = 200000000 / 1000 = 200000`.

The comparison uses Poisson-estimated per-bin error (from the `--naver` flag),
so the criterion is statistical agreement at ±3 sigma.  The `--rmin 1.5` flag
skips the sphere interior (which is all zeros) and `--rmax 7` avoids the
low-count outer tail.

```bash
python3 tools/compare_profiles.py \
    /tmp/val_legacy/denrunt.out \
    /tmp/val_modern_exact/denrunt.out \
    --naver 200000 \
    --tol-sigmas 3 \
    --rmin 1.5 \
    --rmax 7
```

**Expected result:** `PASS` — no bins flagged at 3 sigma.  The two codes
implement the same physics with the same move-acceptance criterion; the only
difference is the RNG sequence (different seeds).  At 200 000 accumulations
the statistical error per bin is small enough that any systematic bug would be
clearly visible.

If this run prints `FAIL`, check the worst-bin output: a few marginal bins at
3.0–3.5 sigma are within expected fluctuations; re-run the legacy once more
(it re-seeds from the wall clock) to confirm the discrepancy is real.

---

## 5. Cutoff-vs-exact study

Quantify the error introduced by the pairwise cutoff approximation.
Create a third directory with the same deck but with an explicit `rcut=5.0`
**appended** as line 11 of `runt.inp`.

The committed deck has exactly 10 lines (no 11th line), so a `sed` substitute on
line 11 would be a silent no-op and the run would fall back to exact mode — the
cutoff study must **append** the line, not substitute it.

```bash
mkdir -p /tmp/val_modern_cut
cp "$DECK/runt.ic" /tmp/val_modern_cut/
cp "$DECK/runt.inp" /tmp/val_modern_cut/   # 10-line deck copied as-is

# Set ncon and nskip (lines 1-2 exist in the deck), then APPEND rcut as line 11
sed -i '' '1s/.*/200000000/' /tmp/val_modern_cut/runt.inp
sed -i '' '2s/.*/1000/'      /tmp/val_modern_cut/runt.inp
echo "5.0" >> /tmp/val_modern_cut/runt.inp   # append → file is now exactly 11 lines

# Verify: line 11 must read 5.0 and the file must have exactly 11 lines
echo "=== cutoff inp ===" && cat -n /tmp/val_modern_cut/runt.inp
```

Run it:

```bash
( cd /tmp/val_modern_cut && "$MODERN_RUNT" )
```

Compare exact vs cutoff:

```bash
python3 tools/compare_profiles.py \
    /tmp/val_modern_exact/denrunt.out \
    /tmp/val_modern_cut/denrunt.out \
    --naver 200000 \
    --tol-sigmas 3 \
    --rmin 1.5 \
    --rmax 7
```

**Expected result:** `PASS` for a good cutoff distance (5 sigma is well beyond
the Yukawa decay length).  The worst-bin line reports the maximum discrepancy
in units of the statistical error; record this value.  If it is << 1 sigma,
the cutoff approximation is negligible for this deck.  If it exceeds ~2 sigma,
consider increasing `rcut` or using exact mode.

---

## 6. Improving statistics with multiple seeds

The legacy binary re-seeds from the wall clock on every invocation, so running
it multiple times produces independent trajectories.  For tighter uncertainty
estimates, run N replicas and average the density profiles before comparing:

```bash
for i in 1 2 3 4; do
    RUNDIR="/tmp/val_legacy_$i"
    mkdir -p "$RUNDIR"
    cp "$DECK/runt.ic"  "$RUNDIR/"
    cp /tmp/val_legacy/runt.inp "$RUNDIR/"   # reuse edited inp
    ( cd "$RUNDIR" && /tmp/runt_legacy ) &
done
wait
```

Average the profiles in Python (load all four `denrunt.out` files, stack
along axis 0, take the mean along axis 0) and write the averaged file, then
pass it to `compare_profiles.py` with `--naver $((4 * 200000))`.

---

## 7. polyfj (athermal) validation

`polyfj` is the athermal hard-chain / hard-sphere program. It has a central
hard sphere (radius `rsp`, from `polyfj.ic`), so density starts beyond the
contact distance `rsp + 0.5`. The comparator reads the first three columns of
`polyden.out` (`idx  dist  total_density`), so it works on `polyden.out`
directly (the extra end/mid columns are ignored).

### 7.1 Build both
```bash
# legacy polyfj
gfortran -std=legacy -O2 -o /tmp/polyfj_legacy polyfj/polyfj.f polyfj/polyfj-moves.f
# modern polyfj (fpm)
fpm build
MODERN_POLYFJ=$(find build -type f -path '*/app/polyfj' | head -1)
```

### 7.2 Set up two run dirs from a deck (e.g. n8, sphere 1.5, eta 0.2)
```bash
DECK=polyfj/runs/n8_s1.5_e0.2
mkdir -p /tmp/val_polyfj_legacy /tmp/val_polyfj_modern
for d in /tmp/val_polyfj_legacy /tmp/val_polyfj_modern; do
  cp "$DECK/polyfj.ic" "$DECK/polyfj.inp" "$d/"
  sed -i '' '1s/.*/200000000/' "$d/polyfj.inp"   # ncon (converged)
  sed -i '' '2s/.*/1000/'      "$d/polyfj.inp"   # nskip
done
```

### 7.3 Run both (legacy is wall-clock seeded; modern is fixed-seed)
```bash
( cd /tmp/val_polyfj_legacy && /tmp/polyfj_legacy )
( cd /tmp/val_polyfj_modern && "$MODERN_POLYFJ" )
```

### 7.4 Compare the density profile (naver = ncon/nskip = 200000)
```bash
python3 tools/compare_profiles.py \
    /tmp/val_polyfj_legacy/polyden.out \
    /tmp/val_polyfj_modern/polyden.out \
    --naver 200000 --tol-sigmas 3 --rmin 2.0 --rmax 7
```
Expected: **PASS** — the modern athermal density reproduces legacy within
Monte-Carlo error.

### 7.5 Structural diagnostics (seg.out / g2fj.out) — note
The modern driver writes **corrected** segmental-order and shape semi-axis
diagnostics by default, which intentionally do NOT match the legacy
`seg.out`/`g2fj.out` (the legacy formulas carry reduced-unit artifacts — an
`AL^2` factor in seg, an inertia-tensor formula in shape). To emit
legacy-faithful values for a direct comparison, run the modern driver with the
`--legacy-obs` flag:
```bash
( cd /tmp/val_polyfj_modern && "$MODERN_POLYFJ" --legacy-obs )
```
The total density in `polyden.out` is identical either way (the flag only
affects `seg.out` and `g2fj.out`).
