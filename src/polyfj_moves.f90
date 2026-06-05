! polyfj_moves.f90 — athermal (hard-core) MC moves for polyfj
!
! Faithful port of the legacy polyfj/polyfj-moves.f (DICK, REPT, CCB) to
! modern Fortran with real (sigma) units.
!
! Key athermal differences from runt_moves:
!   - No energy / Boltzmann weights.  A trial bead is valid (weight 1) iff it
!     has zero hard overlaps; invalid (weight 0) otherwise.
!   - CCB acceptance is purely WO/WN (Rosenbluth ratio of overlap counts),
!     with no exp(EOLD-ENEW) factor.
!   - Coordinates are in real sigma units (bond length = 1 sigma exactly).
!
! Reversal semantics:
!   The legacy CVERT subroutine reversed the chain in the global position array
!   AND immediately updated the linked-cell data structure.  In the modern port
!   we must NOT modify sys%x/y/z before acceptance, because the cell list is
!   only rebuilt by the driver after each accepted move.  Premature modification
!   of sys%x/y/z would leave the cell list stale, causing missed overlap
!   detections in subsequent moves.  The fix: copy the chain into a local
!   working array (xn/yn/zn), apply the reversal there, and only write back to
!   sys on acceptance.
!
! Public interface:
!   subroutine move_dick(sys, p, cl, rng, imol, accepted)
!   subroutine move_rept(sys, p, cl, rng, imol, accepted)
!   subroutine move_ccb (sys, p, cl, rng, imol, accepted)
!
! Caller contract:
!   - cl is up-to-date for the current config; the caller rebuilds cl after
!     each accepted move.  Do NOT rebuild cl inside a move.
!   - Overlap checks use exclude_mol=imol so the moving chain's own old beads
!     are not counted as inter-molecular neighbours.
!   - sys%rsp is the central hard-sphere radius (read from polyfj.ic).
!     If rsp==0, sphere_overlaps always returns .false. (no sphere present).

module polyfj_moves
  use polymc_kinds,     only: dp
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t
  use polymc_rng,       only: rng_t, rng_uniform, rng_unit_vector
  use polymc_overlap,   only: bead_overlaps, sphere_overlaps
  use polymc_config,    only: system_t, params_t
  implicit none
  private

  public :: move_dick, move_rept, move_ccb

  ! Number of Rosenbluth trial directions per bead (matches legacy NSAMP=15)
  integer, parameter :: NSAMP = 15

  ! Branch probability constants matching legacy literal values
  real(dp), parameter :: PROB_BOND_PERTURB  = 0.25_dp   ! DICK: bond perturbation probability
  real(dp), parameter :: PROB_CHAIN_REVERSE = 0.5_dp    ! REPT/CCB: chain reversal probability

contains

  ! ------------------------------------------------------------------
  !> Minimum-image squared distance between two points (local helper).
  ! ------------------------------------------------------------------
  pure real(dp) function mi_dist2(ax, ay, az, bx, by, bz, L)
    real(dp), intent(in) :: ax, ay, az, bx, by, bz, L
    real(dp) :: dx, dy, dz
    dx = min_image(ax - bx, L)
    dy = min_image(ay - by, L)
    dz = min_image(az - bz, L)
    mi_dist2 = dx*dx + dy*dy + dz*dz
  end function mi_dist2

  ! ------------------------------------------------------------------
  !> Reverse a local chain array in place (n beads).
  ! ------------------------------------------------------------------
  pure subroutine reverse_local(x, y, z, n)
    real(dp), intent(inout) :: x(:), y(:), z(:)
    integer,  intent(in)    :: n
    integer  :: j
    real(dp) :: tx, ty, tz
    do j = 1, n / 2
      tx = x(j);     x(j) = x(n+1-j); x(n+1-j) = tx
      ty = y(j);     y(j) = y(n+1-j); y(n+1-j) = ty
      tz = z(j);     z(j) = z(n+1-j); z(n+1-j) = tz
    end do
  end subroutine reverse_local

  ! ====================================================================
  !> move_dick — Dickman regrowth move (athermal).
  !
  ! Port of legacy DICK in polyfj/polyfj-moves.f.
  !
  ! Algorithm:
  !   1. Displace bead 1 by (rng_uniform-0.5)*dlr on each axis.
  !   2. Hard-reject on sphere overlap or inter-molecular bead overlap.
  !   3. For beads ii=2..N: take old bond direction; with prob 0.25 perturb
  !      by dint*unit_vector then renormalise to unit length; step from
  !      previous trial bead.
  !   4. Hard-reject each trial bead on: sphere overlap, inter overlap, or
  !      intra overlap vs trial beads 1..ii-2 (non-bonded, dist<1 sigma).
  !   5. If all N beads placed without overlap → accept, write trial into sys.
  !
  ! No chain reversal in DICK (only REPT and CCB reverse).
  ! ====================================================================
  subroutine move_dick(sys, p, cl, rng, imol, accepted)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol     ! 1-based molecule index
    logical,           intent(out)   :: accepted

    integer  :: ibase, ii, jj
    real(dp) :: tx(sys%n), ty(sys%n), tz(sys%n)   ! trial chain
    real(dp) :: d1, d2, d3, cnm, cn
    real(dp) :: vx, vy, vz, rdis

    accepted = .false.
    ibase    = (imol - 1) * sys%n

    ! ------------------------------------------------------------------
    ! Step 1: build trial bead 1 — displace by random vector in [-dlr/2, +dlr/2]^3
    ! Legacy: ZS1=(RANF-0.5); XITR(1) = X(IBD) + ZS1*DLR
    ! ------------------------------------------------------------------
    tx(1) = sys%x(ibase + 1) + (rng_uniform(rng) - 0.5_dp) * p%dlr
    ty(1) = sys%y(ibase + 1) + (rng_uniform(rng) - 0.5_dp) * p%dlr
    tz(1) = sys%z(ibase + 1) + (rng_uniform(rng) - 0.5_dp) * p%dlr

    ! Hard reject: overlap with central sphere
    if (sphere_overlaps(sys%box, sys%rsp, tx(1), ty(1), tz(1))) return

    ! Hard reject: inter-molecular bead overlap
    if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                      sys%mol_of, tx(1), ty(1), tz(1), imol)) return

    ! ------------------------------------------------------------------
    ! Steps 2-4: build trial beads ii=2..N
    ! Legacy DO 200 II=2,N
    ! ------------------------------------------------------------------
    do ii = 2, sys%n
      ! Old bond vector (unit length in stored config)
      d1 = sys%x(ibase + ii) - sys%x(ibase + ii - 1)
      d2 = sys%y(ibase + ii) - sys%y(ibase + ii - 1)
      d3 = sys%z(ibase + ii) - sys%z(ibase + ii - 1)

      ! With probability 0.25, perturb bond direction and renormalise.
      ! Legacy: IF(RANF.GT.0.25)GOTO 211  → perturb when RANF <= 0.25
      if (rng_uniform(rng) <= PROB_BOND_PERTURB) then
        call rng_unit_vector(rng, vx, vy, vz)
        d1 = d1 + p%dint * vx
        d2 = d2 + p%dint * vy
        d3 = d3 + p%dint * vz
        ! Renormalise to unit length (sigma units)
        cnm = d1*d1 + d2*d2 + d3*d3
        cn  = 1.0_dp / sqrt(cnm)
        d1  = d1 * cn
        d2  = d2 * cn
        d3  = d3 * cn
      end if

      ! Step from previous trial bead
      tx(ii) = tx(ii - 1) + d1
      ty(ii) = ty(ii - 1) + d2
      tz(ii) = tz(ii - 1) + d3

      ! Hard reject: overlap with central sphere
      if (sphere_overlaps(sys%box, sys%rsp, tx(ii), ty(ii), tz(ii))) return

      ! Hard reject: inter-molecular bead overlap
      if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                        sys%mol_of, tx(ii), ty(ii), tz(ii), imol)) return

      ! Hard reject: intra-molecular overlap vs non-bonded earlier trial beads
      ! Legacy DO 215 JJ=1,II-2 (only runs for II>=3)
      if (ii > 2) then
        do jj = 1, ii - 2
          rdis = mi_dist2(tx(jj), ty(jj), tz(jj), tx(ii), ty(ii), tz(ii), &
                          sys%box%L)
          if (rdis < 1.0_dp) return
        end do
      end if
    end do

    ! ------------------------------------------------------------------
    ! All N trial beads accepted with no hard overlap → accept move.
    ! ------------------------------------------------------------------
    accepted = .true.
    sys%x(ibase + 1 : ibase + sys%n) = tx(1:sys%n)
    sys%y(ibase + 1 : ibase + sys%n) = ty(1:sys%n)
    sys%z(ibase + 1 : ibase + sys%n) = tz(1:sys%n)

  end subroutine move_dick

  ! ====================================================================
  !> move_rept — Reptation move (athermal).
  !
  ! Port of legacy REPT in polyfj/polyfj-moves.f.
  !
  ! Algorithm:
  !   1. Copy chain into local working array; with prob 0.5, reverse it.
  !   2. Grow new head at xw(1) + unit_vector (1-sigma bond).
  !   3. Hard-reject: sphere overlap, inter overlap, intra overlap vs
  !      surviving beads 2..N-1.
  !   4. Accept → new chain = [new_head, xw(1), xw(2), ..., xw(N-1)].
  !
  ! IMPORTANT: The reversal is done in the LOCAL working array, NOT in
  ! sys%x/y/z.  This keeps sys and the cell list consistent at all times.
  ! The legacy CVERT modified sys in place AND updated the cell data
  ! structure; here we avoid that by working in a copy.
  ! ====================================================================
  subroutine move_rept(sys, p, cl, rng, imol, accepted)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol
    logical,           intent(out)   :: accepted

    integer  :: ibase, ii
    real(dp) :: xw(sys%n), yw(sys%n), zw(sys%n)   ! local working copy
    real(dp) :: vx, vy, vz
    real(dp) :: tnx, tny, tnz   ! new head bead position
    real(dp) :: rdis

    accepted = .false.
    ibase    = (imol - 1) * sys%n

    ! ------------------------------------------------------------------
    ! Step 1: copy chain into working array; optionally reverse.
    ! Legacy: IF(RANF.GT.0.5) CALL CVERT(I)
    ! We reverse in the local copy to avoid modifying sys before acceptance.
    ! ------------------------------------------------------------------
    do ii = 1, sys%n
      xw(ii) = sys%x(ibase + ii)
      yw(ii) = sys%y(ibase + ii)
      zw(ii) = sys%z(ibase + ii)
    end do
    if (rng_uniform(rng) > PROB_CHAIN_REVERSE) then
      call reverse_local(xw, yw, zw, sys%n)
    end if

    ! ------------------------------------------------------------------
    ! Step 2: grow new head bead at xw(1) + unit_vector (1-sigma bond).
    ! ------------------------------------------------------------------
    call rng_unit_vector(rng, vx, vy, vz)
    tnx = xw(1) + vx
    tny = yw(1) + vy
    tnz = zw(1) + vz

    ! Hard reject: overlap with central sphere
    if (sphere_overlaps(sys%box, sys%rsp, tnx, tny, tnz)) return

    ! Hard reject: inter-molecular bead overlap
    ! (sys%x/y/z still holds original un-reversed positions; cl is valid)
    if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                      sys%mol_of, tnx, tny, tnz, imol)) return

    ! Hard reject: intra overlap vs surviving beads 2..N-1 of working chain.
    ! Legacy DO 320 II=2,N-1
    do ii = 2, sys%n - 1
      rdis = mi_dist2(tnx, tny, tnz, xw(ii), yw(ii), zw(ii), sys%box%L)
      if (rdis < 1.0_dp) return
    end do

    ! ------------------------------------------------------------------
    ! Move accepted.  New chain = [tnx, xw(1), xw(2), ..., xw(N-1)].
    ! Write directly into sys (shift xw(1..N-1) → positions 2..N, head at 1).
    ! ------------------------------------------------------------------
    accepted = .true.
    do ii = sys%n, 2, -1
      sys%x(ibase + ii) = xw(ii - 1)
      sys%y(ibase + ii) = yw(ii - 1)
      sys%z(ibase + ii) = zw(ii - 1)
    end do
    sys%x(ibase + 1) = tnx
    sys%y(ibase + 1) = tny
    sys%z(ibase + 1) = tnz

  end subroutine move_rept

  ! ====================================================================
  !> move_ccb — Configurational-bias MC (athermal Rosenbluth) move.
  !
  ! Port of legacy CCB in polyfj/polyfj-moves.f.
  !
  ! In this athermal version the Rosenbluth weights are purely overlap-based:
  !   ET(K) = 1.0  if trial direction K has no hard overlap
  !   ET(K) = 0.0  if trial direction K has a hard overlap
  ! There is NO Boltzmann factor (no exp(-energy)).
  !
  ! Algorithm:
  !   1. Copy chain into xn; with prob 0.5, reverse xn in place.
  !      (Never touch sys%x/y/z so the cell list stays valid.)
  !   2. ICUT = INT((N-1)*uniform) + 2  [2..N]
  !   3. NEW (forward) Rosenbluth: regrow beads ICUT..N.
  !        For each bead J: generate NSAMP random unit-bond trials from xn(J-1).
  !        ET(K) = 1 if no sphere/inter/intra overlap; 0 otherwise.
  !        SUM = sum(ET); if SUM < tiny → hard reject entire move.
  !        Normalise ET/=SUM; select trial proportionally; WN *= ET_selected.
  !        (ET_selected = 1/num_valid_trials.)
  !   4. Store new chain as stx/sty/stz; reload xn from sys (un-reversed original).
  !   5. OLD (backward) Rosenbluth: for beads ICUT..N.
  !        K=1: existing bead (ET(1)=1 always, since config is valid).
  !        K=2..NSAMP: fresh random trials, ET=1/0.
  !        ET(1) /= SUM (normalise ET(1) only); WO *= ET(1).
  !        (ET(1) = 1/total_valid_trials_including_old.)
  !   6. Accept iff WO/WN > rng_uniform (legacy line 313).
  !   7. On accept: write stx/sty/stz into sys.
  !
  ! IMPORTANT: The reversal is applied to the local copy xn only.  sys%x/y/z
  ! is NOT modified before acceptance, keeping the cell list valid throughout.
  !
  ! CCB weight bookkeeping (athermal):
  !   WN = product over regrown beads of (1/num_valid_new_directions).
  !   WO = product over regrown beads of (1/num_valid_including_old_direction).
  !   Acceptance probability = WO/WN.
  ! ====================================================================
  subroutine move_ccb(sys, p, cl, rng, imol, accepted)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol
    logical,           intent(out)   :: accepted

    integer  :: ibase, j, k, ii, icut, i1
    real(dp) :: xn(sys%n), yn(sys%n), zn(sys%n)    ! working chain (may be reversed)
    real(dp) :: xbase(sys%n), ybase(sys%n), zbase(sys%n) ! base chain after reversal
    real(dp) :: stx(sys%n), sty(sys%n), stz(sys%n) ! stored new chain
    real(dp) :: xt(NSAMP), yt(NSAMP), zt(NSAMP)    ! trial positions
    real(dp) :: et(NSAMP)                           ! trial weights (0 or 1, then normalised)
    real(dp) :: wn, wo, wsum, xran, s
    real(dp) :: vx, vy, vz, rdis

    accepted = .false.
    wn       = 1.0_dp
    wo       = 1.0_dp
    ibase    = (imol - 1) * sys%n

    ! ------------------------------------------------------------------
    ! Step 1: copy chain into working array; optionally reverse in place.
    !
    ! Legacy: CVERT modifies X1 in place (reversal) AND updates the cell
    ! list. Then XN is initialised from X1. After NEW phase, XN is reloaded
    ! from X1 (still reversed) for the OLD Rosenbluth step.
    !
    ! Modern fix: never modify sys%x/y/z before acceptance (keeps the cell
    ! list valid). Instead we copy sys → xbase, optionally reverse xbase,
    ! and use xbase as the starting point for BOTH the NEW and OLD
    ! Rosenbluth phases. sys%x/y/z is only written on acceptance.
    ! ------------------------------------------------------------------
    do j = 1, sys%n
      i1       = ibase + j
      xbase(j) = sys%x(i1)
      ybase(j) = sys%y(i1)
      zbase(j) = sys%z(i1)
    end do
    if (rng_uniform(rng) > PROB_CHAIN_REVERSE) then
      call reverse_local(xbase, ybase, zbase, sys%n)
    end if

    ! Initialise working array from (possibly reversed) base.
    do j = 1, sys%n
      xn(j) = xbase(j)
      yn(j) = ybase(j)
      zn(j) = zbase(j)
    end do

    ! ------------------------------------------------------------------
    ! Step 2: choose cut point ICUT in [2, N].
    ! Legacy: ICUT = INT((N-1)*RANF(NSEED)) + 2
    ! ------------------------------------------------------------------
    icut = int(real(sys%n - 1, dp) * rng_uniform(rng)) + 2

    ! ------------------------------------------------------------------
    ! Step 3: NEW (forward) Rosenbluth — regrow beads ICUT..N.
    ! ------------------------------------------------------------------
    do j = icut, sys%n

      wsum = 0.0_dp

      do k = 1, NSAMP
        et(k) = 1.0_dp   ! default: valid

        ! Generate random unit bond from working bead J-1
        call rng_unit_vector(rng, vx, vy, vz)
        xt(k) = xn(j - 1) + vx
        yt(k) = yn(j - 1) + vy
        zt(k) = zn(j - 1) + vz

        ! Check sphere overlap — set ET(K)=0, accumulate, go to next k
        if (sphere_overlaps(sys%box, sys%rsp, xt(k), yt(k), zt(k))) then
          et(k) = 0.0_dp
          wsum  = wsum + et(k)
          cycle
        end if

        ! Check inter-molecular bead overlap (using unmodified sys%x/y/z)
        if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                          sys%mol_of, xt(k), yt(k), zt(k), imol)) then
          et(k) = 0.0_dp
          wsum  = wsum + et(k)
          cycle
        end if

        ! Check intra-molecular overlap vs non-bonded beads 1..J-2.
        ! Legacy DO 16 II=1,J-2 (only runs for J>=3).
        if (j >= 3) then
          do ii = 1, j - 2
            rdis = mi_dist2(xt(k), yt(k), zt(k), xn(ii), yn(ii), zn(ii), &
                            sys%box%L)
            if (rdis < 1.0_dp) then
              et(k) = 0.0_dp
              exit   ! found intra overlap; et(k) is now 0
            end if
          end do
        end if

        wsum = wsum + et(k)
      end do  ! k = 1..NSAMP

      ! Hard reject: no valid trial position found.
      if (wsum < 1.0d-10) return

      ! Normalise ET(K) by SUM.
      wsum = 1.0_dp / wsum
      do k = 1, NSAMP
        et(k) = et(k) * wsum
      end do

      ! Select trial proportionally to ET(K).
      ! Fallback k=NSAMP is safe because wsum > 0 ensures at least one ET > 0.
      xran = rng_uniform(rng)
      s    = 0.0_dp
      k    = NSAMP
      block
        integer :: kk
        do kk = 1, NSAMP
          s = s + et(kk)
          if (xran < s) then
            k = kk
            exit
          end if
        end do
      end block

      ! Accumulate Rosenbluth weight for new chain.
      wn    = wn * et(k)
      xn(j) = xt(k)
      yn(j) = yt(k)
      zn(j) = zt(k)

    end do  ! j = icut..n  (NEW chain grown)

    ! ------------------------------------------------------------------
    ! Step 4: store new chain; set up OLD chain for weight calculation.
    ! The OLD chain uses the ORIGINAL (un-reversed) chain positions.
    ! Legacy: after DO 50 loop, copy XN→STX; reload XN from X1 (the original).
    ! Here, xorig holds the original positions regardless of reversal.
    ! ------------------------------------------------------------------
    do j = 1, sys%n
      stx(j) = xn(j)
      sty(j) = yn(j)
      stz(j) = zn(j)
    end do

    ! Reload xn from the base chain (possibly reversed) for OLD Rosenbluth.
    ! Legacy: "DO 65 J=1,N: XN(J) = X1(IMOL+J)" where X1 was already
    ! modified by CVERT. Here xbase holds the same post-reversal state.
    do j = 1, sys%n
      xn(j) = xbase(j)
      yn(j) = ybase(j)
      zn(j) = zbase(j)
    end do

    ! ------------------------------------------------------------------
    ! Step 5: OLD (backward) Rosenbluth — compute WO.
    ! For each bead J=ICUT..N:
    !   K=1: use existing original bead (ET(1)=1 always; valid by assumption).
    !   K=2..NSAMP: generate fresh random trials, ET=1/0.
    !   ET(1) /= SUM; WO *= ET(1).
    !
    ! NOTE: xn here holds the BASE (possibly reversed) chain. This matches
    ! legacy which reloads XN from X1 (which was already reversed by CVERT)
    ! after storing the new chain.
    ! ------------------------------------------------------------------
    do j = icut, sys%n

      wsum = 0.0_dp

      do k = 1, NSAMP

        if (k == 1) then
          ! Existing bead: always valid in a non-overlapping configuration.
          et(k) = 1.0_dp
          wsum  = wsum + et(k)
          cycle
        end if

        ! K >= 2: generate fresh random trial from old chain position J-1.
        call rng_unit_vector(rng, vx, vy, vz)
        xt(k) = xn(j - 1) + vx
        yt(k) = yn(j - 1) + vy
        zt(k) = zn(j - 1) + vz

        et(k) = 1.0_dp   ! default: valid

        ! Sphere overlap
        if (sphere_overlaps(sys%box, sys%rsp, xt(k), yt(k), zt(k))) then
          et(k) = 0.0_dp
          wsum  = wsum + et(k)
          cycle
        end if

        ! Inter-molecular bead overlap
        if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                          sys%mol_of, xt(k), yt(k), zt(k), imol)) then
          et(k) = 0.0_dp
          wsum  = wsum + et(k)
          cycle
        end if

        ! Intra-molecular overlap vs non-bonded beads 1..J-2 (original chain)
        if (j >= 3) then
          do ii = 1, j - 2
            rdis = mi_dist2(xt(k), yt(k), zt(k), xn(ii), yn(ii), zn(ii), &
                            sys%box%L)
            if (rdis < 1.0_dp) then
              et(k) = 0.0_dp
              exit
            end if
          end do
        end if

        wsum = wsum + et(k)
      end do  ! k = 1..NSAMP

      ! Normalise ET(1) by SUM and accumulate WO.
      ! Legacy: ET(1) = ET(1)/SUM; WO = WO*ET(1)
      et(1) = et(1) / wsum
      wo    = wo * et(1)

    end do  ! j = icut..n  (OLD chain weights)

    ! ------------------------------------------------------------------
    ! Step 6: Acceptance criterion (athermal Rosenbluth).
    ! Legacy: BOLTZ = WO/WN; IF(BOLTZ.GT.RANF(...)) accept.
    ! ------------------------------------------------------------------
    if (wo / wn <= rng_uniform(rng)) return

    ! ------------------------------------------------------------------
    ! Step 7: accept — write new chain into sys.
    ! ------------------------------------------------------------------
    accepted = .true.
    do j = 1, sys%n
      sys%x(ibase + j) = stx(j)
      sys%y(ibase + j) = sty(j)
      sys%z(ibase + j) = stz(j)
    end do

  end subroutine move_ccb

end module polyfj_moves
