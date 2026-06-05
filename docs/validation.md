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

Set `ncon=200000000` and `nskip=1000`, keep line 11 **blank** (exact mode):

```bash
# In-place edit with sed (line 1 and 2 only, leave all other lines)
for DIR in /tmp/val_legacy /tmp/val_modern_exact; do
    sed -i '' '1s/.*/200000000/' "$DIR/runt.inp"
    sed -i '' '2s/.*/1000/'      "$DIR/runt.inp"
    # Ensure line 11 is blank (exact mode, no cutoff)
    # macOS sed: replace the 11th line with empty
    sed -i '' '11s/.*//' "$DIR/runt.inp"
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
Create a third directory with the same deck but with an explicit `rcut=5.0` on
line 11 of `runt.inp`:

```bash
mkdir -p /tmp/val_modern_cut
cp "$DECK/runt.ic" /tmp/val_modern_cut/
cp "$DECK/runt.inp" /tmp/val_modern_cut/

# Set ncon, nskip, and rcut
sed -i '' '1s/.*/200000000/' /tmp/val_modern_cut/runt.inp
sed -i '' '2s/.*/1000/'      /tmp/val_modern_cut/runt.inp
sed -i '' '11s/.*/5.0/'      /tmp/val_modern_cut/runt.inp

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

## 7. polyfj equivalent (pending Task 19)

Once the modern `polyfj` driver exists, the same workflow applies using
`polyden.out` instead of `denrunt.out`.  The column layout is:

```
idx   dist   total_density   ...
```

Use `--rmin` and `--rmax` appropriate for the polyfj geometry (no central
sphere — start comparison from `dist = 0.5` sigma).

Command template (update paths after Task 19):

```bash
python3 tools/compare_profiles.py \
    /tmp/val_polyfj_legacy/polyden.out \
    /tmp/val_polyfj_modern/polyden.out \
    --naver 200000 \
    --tol-sigmas 3
```
