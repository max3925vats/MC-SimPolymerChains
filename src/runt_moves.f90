! runt_moves.f90 — Monte Carlo moves for the runt polymer model
!
! Implements three moves faithful to the legacy runt-moves.f (DICK / REPT / CCB),
! ported to real (sigma) units.  The legacy code worked in reduced units divided
! by AL; here positions are in sigma so:
!   - dlr and dint are used directly (no /AL division)
!   - bond renormalisation uses CN = 1/sqrt(CNM) giving unit-length bonds
!   - rng_unit_vector returns a unit vector; adding it gives a 1-sigma bond
!
! Reference: runt/runt-moves.f (committed version with RANF + corrected BBATT
! + CVERT).
!
! Public interface:
!   move_dick(sys, p, cl, rng, imol, accepted, de)
!   move_rept(sys, p, cl, rng, imol, accepted, de)
!   move_ccb (sys, p, cl, rng, imol, accepted, de)
!
! Caller contract:
!   - cl reflects the CURRENT config (rebuilt by driver after each accept)
!   - within a move, cl is used with exclude_mol=imol so old beads of the
!     moving chain are skipped; do NOT rebuild cl inside a move
!   - on accept: sys%x/y/z updated; de = ENEW - EOLD; accepted = .true.
!   - on reject: sys unchanged; de = 0; accepted = .false.

module runt_moves
  use polymc_kinds,     only: dp
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t
  use polymc_rng,       only: rng_t, rng_uniform, rng_unit_vector
  use polymc_overlap,   only: bead_overlaps, sphere_overlaps
  use polymc_config,    only: system_t, params_t
  use runt_potential,   only: u_pair, u_sphere
  use runt_energy,      only: bead_inter_energy
  implicit none
  private

  public :: move_dick, move_rept, move_ccb

  ! Number of Rosenbluth trials per bead in CCB (matches legacy NSAMP=15)
  integer, parameter :: NSAMP = 15

  ! Branch probability constants (match legacy literal values exactly)
  real(dp), parameter :: PROB_BOND_PERTURB  = 0.25_dp   ! DICK: chance to perturb a bond
  real(dp), parameter :: PROB_CHAIN_REVERSE = 0.5_dp    ! REPT/CCB: chance to reverse chain first

contains

  ! ------------------------------------------------------------------
  !> Helper: minimum-image squared distance between two points.
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
  !> Helper: minimum-image distance from point to origin.
  ! ------------------------------------------------------------------
  pure real(dp) function dist_to_origin(px, py, pz, L)
    real(dp), intent(in) :: px, py, pz, L
    real(dp) :: dx, dy, dz
    dx = min_image(px, L)
    dy = min_image(py, L)
    dz = min_image(pz, L)
    dist_to_origin = sqrt(dx*dx + dy*dy + dz*dz)
  end function dist_to_origin

  ! ------------------------------------------------------------------
  !> Reverse a local chain array in place (n beads, 1-based).
  !
  ! Used by REPT/CCB to apply the 50% chain reversal to a LOCAL working
  ! copy instead of sys.  See module note on reversal semantics.
  ! ------------------------------------------------------------------
  pure subroutine reverse_local(x, y, z, n)
    real(dp), intent(inout) :: x(:), y(:), z(:)
    integer,  intent(in)    :: n
    integer  :: j
    real(dp) :: tx, ty, tz
    do j = 1, n / 2
      tx = x(j);  x(j) = x(n+1-j);  x(n+1-j) = tx
      ty = y(j);  y(j) = y(n+1-j);  y(n+1-j) = ty
      tz = z(j);  z(j) = z(n+1-j);  z(n+1-j) = tz
    end do
  end subroutine reverse_local

  ! ====================================================================
  !> move_dick — Dickman regrowth move.
  !
  ! Faithful port of legacy DICK in runt-moves.f.
  !
  ! Algorithm:
  !   1. Translate bead 1 by a random displacement in [-dlr/2, +dlr/2]^3
  !   2. For beads 2..N: keep old bond direction; with prob 0.25 perturb by
  !      dint*unit_vector then renormalise to unit length; step from previous bead
  !   3. At each bead: sphere overlap => reject; inter overlap => reject;
  !      intra overlap (vs non-bonded partners) => reject; else accumulate ENEW
  !   4. After all N trial beads built: compute EOLD for the current chain
  !   5. Metropolis acceptance
  !
  ! Unit note: all displacements are in sigma; dlr, dint already in sigma.
  ! ====================================================================
  subroutine move_dick(sys, p, cl, rng, imol, accepted, de)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol    ! 1-based molecule index
    logical,           intent(out)   :: accepted
    real(dp),          intent(out)   :: de

    integer  :: ibase, ii, jj
    real(dp) :: tx(sys%n), ty(sys%n), tz(sys%n)   ! trial chain
    real(dp) :: d1, d2, d3, cnm, cn
    real(dp) :: vx, vy, vz, rdis, r
    real(dp) :: enew, eold

    accepted = .false.
    de       = 0.0_dp
    ibase    = (imol - 1) * sys%n

    ! ------------------------------------------------------------------
    ! Build trial bead 1: displace by random vector in [-dlr/2, +dlr/2]^3
    ! ------------------------------------------------------------------
    tx(1) = sys%x(ibase + 1) + (rng_uniform(rng) - 0.5_dp) * p%dlr
    ty(1) = sys%y(ibase + 1) + (rng_uniform(rng) - 0.5_dp) * p%dlr
    tz(1) = sys%z(ibase + 1) + (rng_uniform(rng) - 0.5_dp) * p%dlr

    ! Sphere overlap check (hard reject)
    if (sphere_overlaps(sys%box, sys%rsp, tx(1), ty(1), tz(1))) return

    ! Inter-molecular bead overlap check (hard reject)
    if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                      sys%mol_of, tx(1), ty(1), tz(1), imol)) return

    ! Accumulate inter + sphere energy for bead 1
    r    = dist_to_origin(tx(1), ty(1), tz(1), sys%box%L)
    enew = bead_inter_energy(sys, p, cl, imol, tx(1), ty(1), tz(1)) &
         + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

    ! ------------------------------------------------------------------
    ! Build trial beads 2..N
    ! ------------------------------------------------------------------
    do ii = 2, sys%n
      ! Old bond vector (real units; bonds are unit length in stored config)
      d1 = sys%x(ibase + ii) - sys%x(ibase + ii - 1)
      d2 = sys%y(ibase + ii) - sys%y(ibase + ii - 1)
      d3 = sys%z(ibase + ii) - sys%z(ibase + ii - 1)

      ! With probability 0.25, perturb bond direction then renormalise
      ! (legacy: IF(RANF(NSEED).GT.0.25)GOTO 211 -> i.e., perturb with prob 0.25)
      if (rng_uniform(rng) <= PROB_BOND_PERTURB) then
        call rng_unit_vector(rng, vx, vy, vz)
        d1 = d1 + p%dint * vx
        d2 = d2 + p%dint * vy
        d3 = d3 + p%dint * vz
        ! Renormalise to unit length (sigma units: target length = 1.0)
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

      ! Sphere overlap check (hard reject)
      if (sphere_overlaps(sys%box, sys%rsp, tx(ii), ty(ii), tz(ii))) return

      ! Inter-molecular bead overlap check (hard reject)
      if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                        sys%mol_of, tx(ii), ty(ii), tz(ii), imol)) return

      ! Accumulate inter + sphere energy for this bead
      r    = dist_to_origin(tx(ii), ty(ii), tz(ii), sys%box%L)
      enew = enew + bead_inter_energy(sys, p, cl, imol, tx(ii), ty(ii), tz(ii)) &
                  + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

      ! Intra-molecular overlap check vs non-bonded partners (j..ii-2)
      ! and accumulate intra energy (faithful to legacy: checks ii-2 partners)
      if (ii > 2) then
        do jj = 1, ii - 2
          rdis = mi_dist2(tx(jj), ty(jj), tz(jj), tx(ii), ty(ii), tz(ii), &
                          sys%box%L)
          ! Hard overlap: dist^2 < 1 (sigma^2) — reject
          if (rdis < 1.0_dp) return
          ! Soft energy contribution
          enew = enew + u_pair(sqrt(rdis), p%beps, p%rcut)
        end do
      end if
    end do

    ! ------------------------------------------------------------------
    ! All N trial beads accepted; compute EOLD for the current chain.
    ! Legacy computes:
    !   - intra pairs: j=3..N, k=1..j-2 (non-bonded)
    !   - per-bead: inter (OLDEN1) + sphere
    ! ------------------------------------------------------------------
    eold = 0.0_dp

    ! Intra-molecular contribution
    do jj = 3, sys%n
      do ii = 1, jj - 2
        rdis = mi_dist2(sys%x(ibase + jj), sys%y(ibase + jj), sys%z(ibase + jj), &
                        sys%x(ibase + ii), sys%y(ibase + ii), sys%z(ibase + ii), &
                        sys%box%L)
        eold = eold + u_pair(sqrt(rdis), p%beps, p%rcut)
      end do
    end do

    ! Inter-molecular + sphere per bead
    do jj = 1, sys%n
      r    = dist_to_origin(sys%x(ibase + jj), sys%y(ibase + jj), &
                            sys%z(ibase + jj), sys%box%L)
      eold = eold + bead_inter_energy(sys, p, cl, imol, &
                      sys%x(ibase + jj), sys%y(ibase + jj), sys%z(ibase + jj)) &
                  + u_sphere(r, p%bbeps, sys%rsp, p%rcut)
    end do

    ! ------------------------------------------------------------------
    ! Metropolis acceptance criterion (legacy: accept if EOLD > ENEW always,
    ! if EOLD <= ENEW accept with exp(EOLD-ENEW))
    ! ------------------------------------------------------------------
    if (eold <= enew) then
      if (exp(eold - enew) <= rng_uniform(rng)) return
    end if

    ! Move accepted
    accepted = .true.
    de       = enew - eold
    sys%x(ibase + 1 : ibase + sys%n) = tx(1:sys%n)
    sys%y(ibase + 1 : ibase + sys%n) = ty(1:sys%n)
    sys%z(ibase + 1 : ibase + sys%n) = tz(1:sys%n)

  end subroutine move_dick

  ! ====================================================================
  !> move_rept — Reptation move.
  !
  ! Faithful port of legacy REPT in runt-moves.f.
  !
  ! Algorithm:
  !   1. With prob 0.5, reverse the chain in sys (CVERT equivalent)
  !   2. Grow new head bead at pos(1) + unit_vector (length 1 sigma bond)
  !   3. Check sphere overlap, inter overlap, intra overlap vs beads 2..N-1
  !   4. Compute ENEW for the new head; compute EOLD for the removed tail (bead N)
  !   5. Metropolis acceptance
  !   6. On accept: new chain = [new_head, old_1, old_2, ..., old_{N-1}]
  !
  ! Reversal semantics (fix): the 50% chain reversal is applied to a LOCAL
  ! working copy (xw/yw/zw), NEVER to sys%x/y/z.  The driver only rebuilds cl
  ! after an accepted move, so mutating sys before acceptance would leave cl
  ! stale and corrupt later overlap detection on a rejected reversal.  All of
  ! the moving chain's own positions (head-growth base, surviving beads, the
  ! removed tail, and the EOLD intra computation) come from the local copy.
  ! Inter-molecular / sphere terms still read sys with exclude_mol=imol, which
  ! is correct: the moving chain's own beads are excluded, so its stale-in-cl
  ! old positions are never consulted.  sys is written only on acceptance, so a
  ! rejected move leaves sys bitwise unchanged.
  ! ====================================================================
  subroutine move_rept(sys, p, cl, rng, imol, accepted, de)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol
    logical,           intent(out)   :: accepted
    real(dp),          intent(out)   :: de

    integer  :: ibase, ii
    real(dp) :: xw(sys%n), yw(sys%n), zw(sys%n)   ! local working copy
    real(dp) :: vx, vy, vz
    real(dp) :: tnx, tny, tnz   ! new head bead position
    real(dp) :: rdis, r
    real(dp) :: enew, eold

    accepted = .false.
    de       = 0.0_dp
    ibase    = (imol - 1) * sys%n

    ! ------------------------------------------------------------------
    ! Step 1: copy chain into local working array; with prob 0.5 reverse
    ! the LOCAL copy (legacy CVERT, but never touching sys before accept).
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
    ! Step 2: grow new head bead at xw(1) + unit_vector
    ! Legacy: XITR(1) = X(IMOL+1) + DELX/AL  where |DELX/AL| = 1/AL = 1 sigma
    ! Modern: directly add unit vector (length 1 sigma)
    ! ------------------------------------------------------------------
    call rng_unit_vector(rng, vx, vy, vz)
    tnx = xw(1) + vx
    tny = yw(1) + vy
    tnz = zw(1) + vz

    ! Sphere overlap check (hard reject)
    if (sphere_overlaps(sys%box, sys%rsp, tnx, tny, tnz)) return

    ! Inter-molecular bead overlap check (hard reject)
    ! (sys still holds original positions; cl is valid; own chain excluded)
    if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                      sys%mol_of, tnx, tny, tnz, imol)) return

    ! Accumulate ENEW for new head
    r    = dist_to_origin(tnx, tny, tnz, sys%box%L)
    enew = bead_inter_energy(sys, p, cl, imol, tnx, tny, tnz) &
         + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

    ! Intra overlap check: new head vs working beads 2..N-1 (surviving beads)
    ! and accumulate intra energy with beads 2..N-1
    ! Legacy DO 320 II=2,N-1
    do ii = 2, sys%n - 1
      rdis = mi_dist2(tnx, tny, tnz, xw(ii), yw(ii), zw(ii), sys%box%L)
      if (rdis < 1.0_dp) return
      enew = enew + u_pair(sqrt(rdis), p%beps, p%rcut)
    end do

    ! ------------------------------------------------------------------
    ! Step 4: compute EOLD for the removed tail bead (working bead N)
    ! Legacy REPT: intra pairs of bead N vs beads 1..N-2, inter, sphere.
    ! Yukawa inter/sphere energy is invariant under the reversal relabeling,
    ! so reading sys for the inter term is equivalent; we use the working
    ! copy for the tail position and intra pairs for consistency.
    ! ------------------------------------------------------------------
    eold = 0.0_dp

    ! Intra: working bead N vs working beads 1..N-2
    do ii = 1, sys%n - 2
      rdis = mi_dist2(xw(ii), yw(ii), zw(ii), &
                      xw(sys%n), yw(sys%n), zw(sys%n), sys%box%L)
      eold = eold + u_pair(sqrt(rdis), p%beps, p%rcut)
    end do

    ! Inter + sphere for removed working bead N
    r    = dist_to_origin(xw(sys%n), yw(sys%n), zw(sys%n), sys%box%L)
    eold = eold + bead_inter_energy(sys, p, cl, imol, &
                    xw(sys%n), yw(sys%n), zw(sys%n)) &
                + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

    ! ------------------------------------------------------------------
    ! Metropolis acceptance
    ! ------------------------------------------------------------------
    if (eold <= enew) then
      if (exp(eold - enew) <= rng_uniform(rng)) return
    end if

    ! ------------------------------------------------------------------
    ! Move accepted: new chain = [new_head, xw(1), xw(2), ..., xw(N-1)]
    ! Written into sys only now (reject leaves sys bitwise unchanged).
    ! ------------------------------------------------------------------
    accepted = .true.
    de       = enew - eold

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
  !> move_ccb — Configurational-bias MC (Rosenbluth) move.
  !
  ! Faithful port of legacy CCB in runt-moves.f.
  !
  ! Algorithm:
  !   1. With prob 0.5, reverse the chain in sys
  !   2. Copy current chain into XN/YN/ZN (the "working" chain for growth)
  !   3. Pick ICUT = INT((N-1)*uniform)+2  (2..N)
  !   4. Forward (NEW) Rosenbluth growth from ICUT..N:
  !        For each bead J: generate NSAMP random trial positions (unit bond from J-1)
  !        For each trial K: compute energy eni; ET(K)=exp(-eni), or 0 on hard overlap
  !        [sphere overlap on ANY trial => hard reject entire move]
  !        Select one trial proportionally; WN *= ET_selected / SUM_ET
  !   5. Store new chain as STX/STY/STZ; reload XN from original sys
  !   6. Backward (OLD) Rosenbluth: for each bead J=ICUT..N, K=1 uses old pos,
  !        K=2..NSAMP uses random trials; [sphere on any random trial => hard reject]
  !        WO *= ET(1) / SUM_ET
  !   7. Energy verification: compute ENEW (XITR=new chain), EOLD (old X) fully
  !      including intra and inter+sphere for all beads
  !   8. Acceptance: BOLTZ = WO/WN * exp(EOLD-ENEW); accept if BOLTZ >= uniform
  !   9. On accept: write STX/STY/STZ into sys
  !
  ! "Hard reject entire move" on sphere overlap: legacy CCB returns with ISUC=1
  ! (reject) if ANY of the NSAMP trials hits the sphere during the growth phase.
  ! This is reproduced here.
  !
  ! Reversal semantics (fix): the 50% chain reversal is applied to a LOCAL base
  ! copy (xbase/ybase/zbase), NEVER to sys%x/y/z.  The cut point ICUT indexes
  ! into this (possibly reversed) base chain — growing from either end is the
  ! whole point of the reversal.  BOTH the NEW and OLD Rosenbluth phases and the
  ! full EOLD/ENEW energy verification use the base copy for the moving chain's
  ! own positions; inter/sphere terms read sys with exclude_mol=imol (correct:
  ! own beads excluded, so stale-in-cl old positions are never consulted).
  ! sys is written only on acceptance, so a rejected move leaves sys bitwise
  ! unchanged.  (The Yukawa inter/sphere energy is invariant under the reversal
  ! relabeling, but we read it consistently with the base copy positions.)
  ! ====================================================================
  subroutine move_ccb(sys, p, cl, rng, imol, accepted, de)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol
    logical,           intent(out)   :: accepted
    real(dp),          intent(out)   :: de

    integer  :: ibase, j, k, ii, icut
    real(dp) :: xn(sys%n), yn(sys%n), zn(sys%n)   ! working/trial chain
    real(dp) :: xbase(sys%n), ybase(sys%n), zbase(sys%n) ! base chain (post-reversal)
    real(dp) :: stx(sys%n), sty(sys%n), stz(sys%n) ! stored new chain
    real(dp) :: xt(NSAMP), yt(NSAMP), zt(NSAMP)    ! trial positions
    real(dp) :: et(NSAMP)                           ! trial weights
    real(dp) :: wn, wo, wsum, xran, s
    real(dp) :: eni, eno, rdis
    real(dp) :: vx, vy, vz, r
    real(dp) :: enew, eold, boltz

    accepted = .false.
    de       = 0.0_dp
    ibase    = (imol - 1) * sys%n

    wn = 1.0_dp
    wo = 1.0_dp

    ! ------------------------------------------------------------------
    ! Step 1: copy current chain into a LOCAL base array; with prob 0.5
    ! reverse the LOCAL copy.  sys%x/y/z is never touched before accept,
    ! so the cell list stays consistent on a rejected move.
    ! ------------------------------------------------------------------
    do j = 1, sys%n
      xbase(j) = sys%x(ibase + j)
      ybase(j) = sys%y(ibase + j)
      zbase(j) = sys%z(ibase + j)
    end do
    if (rng_uniform(rng) > PROB_CHAIN_REVERSE) then
      call reverse_local(xbase, ybase, zbase, sys%n)
    end if

    ! ------------------------------------------------------------------
    ! Step 2: initialise working array XN from the (possibly reversed) base
    ! ------------------------------------------------------------------
    do j = 1, sys%n
      xn(j) = xbase(j)
      yn(j) = ybase(j)
      zn(j) = zbase(j)
    end do

    ! ------------------------------------------------------------------
    ! Step 3: choose cut point ICUT in [2, N]
    ! Legacy: ICUT = INT((N-1)*RANF(NSEED)) + 2
    ! ------------------------------------------------------------------
    icut = int(real(sys%n - 1, dp) * rng_uniform(rng)) + 2

    ! ------------------------------------------------------------------
    ! Step 4: NEW (forward) Rosenbluth growth — regrow beads ICUT..N
    ! ------------------------------------------------------------------
    do j = icut, sys%n

      wsum = 0.0_dp

      do k = 1, NSAMP
        eni  = 0.0_dp
        et(k) = 1.0_dp   ! default (will be overwritten)

        ! Generate random unit bond from bead J-1
        call rng_unit_vector(rng, vx, vy, vz)
        xt(k) = xn(j - 1) + vx
        yt(k) = yn(j - 1) + vy
        zt(k) = zn(j - 1) + vz

        ! Sphere overlap: hard reject entire move (legacy: RETURN with ISUC=1)
        if (sphere_overlaps(sys%box, sys%rsp, xt(k), yt(k), zt(k))) return

        ! Inter-molecular bead overlap: ET(K)=0, continue to next trial
        if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                          sys%mol_of, xt(k), yt(k), zt(k), imol)) then
          et(k) = 0.0_dp
          wsum  = wsum + et(k)
          cycle
        end if

        ! Inter + sphere energy for this trial bead
        r   = dist_to_origin(xt(k), yt(k), zt(k), sys%box%L)
        eni = bead_inter_energy(sys, p, cl, imol, xt(k), yt(k), zt(k)) &
            + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

        ! Intra overlap check vs non-bonded partners j-2 and below
        ! and accumulate intra energy
        block
          logical :: intra_overlap
          intra_overlap = .false.
          if (j >= 3) then
            do ii = 1, j - 2
              rdis = mi_dist2(xt(k), yt(k), zt(k), xn(ii), yn(ii), zn(ii), &
                              sys%box%L)
              if (rdis < 1.0_dp) then
                intra_overlap = .true.
                exit
              end if
              eni = eni + u_pair(sqrt(rdis), p%beps, p%rcut)
            end do
          end if
          ! On hard overlap: et(k)=0; otherwise accumulate Boltzmann weight.
          ! wsum always incremented (preserves goto fall-through semantics).
          if (intra_overlap) then
            et(k) = 0.0_dp
          else
            et(k) = exp(-eni)
          end if
        end block
        wsum = wsum + et(k)
      end do  ! k = 1..NSAMP

      ! If all weights zero: hard reject (no valid position)
      if (wsum < 1.0d-10) return

      ! Normalise
      wsum = 1.0_dp / wsum
      do k = 1, NSAMP
        et(k) = et(k) * wsum
      end do

      ! Select trial proportionally to ET(K)
      xran = rng_uniform(rng)
      s    = 0.0_dp
      k    = NSAMP    ! fallback to last trial
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

      ! Accumulate normalised Rosenbluth weight: WN *= ET_selected/sum_ET (matches legacy)
      wn    = wn * et(k)
      ! Update working chain
      xn(j) = xt(k)
      yn(j) = yt(k)
      zn(j) = zt(k)

    end do  ! j = icut..n  (NEW chain)

    ! ------------------------------------------------------------------
    ! Step 5: store new chain; reload XN from the (possibly reversed) base
    ! chain for the OLD Rosenbluth phase.  Legacy reloads XN from X1, which
    ! CVERT had already reversed in place; xbase holds that same state.
    ! ------------------------------------------------------------------
    do j = 1, sys%n
      stx(j) = xn(j)
      sty(j) = yn(j)
      stz(j) = zn(j)
    end do

    do j = 1, sys%n
      xn(j) = xbase(j)
      yn(j) = ybase(j)
      zn(j) = zbase(j)
    end do

    ! ------------------------------------------------------------------
    ! Step 6: OLD (backward) Rosenbluth — compute WO
    ! For each bead J=ICUT..N: K=1 uses the existing bead position;
    ! K=2..NSAMP uses fresh random trials.
    ! WO *= ET(1) / SUM_ET  (normalised weight of the old position)
    ! ------------------------------------------------------------------
    do j = icut, sys%n

      wsum = 0.0_dp

      do k = 1, NSAMP

        eno = 0.0_dp

        if (k == 1) then
          ! Use the existing (base-chain) bead position J, already verified
          ! non-overlapping.  xn currently holds the (possibly reversed) base.
          r   = dist_to_origin(xn(j), yn(j), zn(j), sys%box%L)
          eno = bead_inter_energy(sys, p, cl, imol, xn(j), yn(j), zn(j)) &
              + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

          ! Intra contribution for existing bead vs base beads 1..J-2
          if (j >= 3) then
            do ii = 1, j - 2
              rdis = mi_dist2(xn(j), yn(j), zn(j), &
                              xn(ii), yn(ii), zn(ii), sys%box%L)
              ! Legacy: IF(RDIS.LT.DEFF) WRITE error (not a reject) — config is valid
              eno = eno + u_pair(sqrt(rdis), p%beps, p%rcut)
            end do
          end if

          et(k) = exp(-eno)
          wsum  = wsum + et(k)
          cycle  ! go to next k
        end if

        ! k >= 2: generate random trial
        call rng_unit_vector(rng, vx, vy, vz)
        xt(k) = xn(j - 1) + vx
        yt(k) = yn(j - 1) + vy
        zt(k) = zn(j - 1) + vz

        ! Sphere overlap: hard reject entire move
        if (sphere_overlaps(sys%box, sys%rsp, xt(k), yt(k), zt(k))) return

        ! Inter-molecular bead overlap: ET(K)=0
        if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                          sys%mol_of, xt(k), yt(k), zt(k), imol)) then
          et(k) = 0.0_dp
          wsum  = wsum + et(k)
          cycle
        end if

        ! Inter + sphere
        r   = dist_to_origin(xt(k), yt(k), zt(k), sys%box%L)
        eno = bead_inter_energy(sys, p, cl, imol, xt(k), yt(k), zt(k)) &
            + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

        ! Intra
        block
          logical :: intra_overlap
          intra_overlap = .false.
          if (j >= 3) then
            do ii = 1, j - 2
              rdis = mi_dist2(xt(k), yt(k), zt(k), xn(ii), yn(ii), zn(ii), &
                              sys%box%L)
              if (rdis < 1.0_dp) then
                intra_overlap = .true.
                exit
              end if
              eno = eno + u_pair(sqrt(rdis), p%beps, p%rcut)
            end do
          end if
          ! On hard overlap: et(k)=0; otherwise accumulate Boltzmann weight.
          ! wsum always incremented (preserves goto fall-through semantics).
          if (intra_overlap) then
            et(k) = 0.0_dp
          else
            et(k) = exp(-eno)
          end if
        end block
        wsum = wsum + et(k)

      end do  ! k = 1..NSAMP

      ! WO gets the normalised weight of the existing (K=1) position
      et(1) = et(1) / wsum
      wo    = wo * et(1)

    end do  ! j = icut..n  (OLD chain)

    ! ------------------------------------------------------------------
    ! Step 7: full energy verification (legacy CCB "CALCULATE AND VERIFY ENERGIES")
    ! ENEW = full energy of new chain (STX/STY/STZ stored in XN temporarily via XITR)
    ! EOLD = full energy of old chain (sys%x/y/z)
    ! Both include inter+sphere per bead AND intra non-bonded pairs.
    ! ------------------------------------------------------------------
    eold = 0.0_dp
    enew = 0.0_dp

    do j = 1, sys%n

      ! NEW chain: inter + sphere for bead J
      r    = dist_to_origin(stx(j), sty(j), stz(j), sys%box%L)
      enew = enew + bead_inter_energy(sys, p, cl, imol, stx(j), sty(j), stz(j)) &
                  + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

      ! OLD chain: inter + sphere for (possibly reversed) base bead J
      r    = dist_to_origin(xbase(j), ybase(j), zbase(j), sys%box%L)
      eold = eold + bead_inter_energy(sys, p, cl, imol, &
                      xbase(j), ybase(j), zbase(j)) &
                  + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

      ! Intra pairs: J vs beads 1..J-2
      if (j >= 3) then
        do k = 1, j - 2
          ! New chain intra
          rdis = mi_dist2(stx(j), sty(j), stz(j), stx(k), sty(k), stz(k), &
                          sys%box%L)
          enew = enew + u_pair(sqrt(rdis), p%beps, p%rcut)

          ! Old chain intra (base copy)
          rdis = mi_dist2(xbase(j), ybase(j), zbase(j), &
                          xbase(k), ybase(k), zbase(k), sys%box%L)
          eold = eold + u_pair(sqrt(rdis), p%beps, p%rcut)
        end do
      end if

    end do

    ! ------------------------------------------------------------------
    ! Step 8: acceptance criterion (Rosenbluth-weighted Metropolis)
    ! Legacy: BOLTZ = WO/WN * exp(EOLD-ENEW); accept if BOLTZ >= RANF
    ! ------------------------------------------------------------------
    boltz = wo / wn * exp(eold - enew)
    if (boltz < rng_uniform(rng)) return

    ! ------------------------------------------------------------------
    ! Step 9: accept — write new chain into sys
    ! ------------------------------------------------------------------
    accepted = .true.
    de       = enew - eold

    do j = 1, sys%n
      sys%x(ibase + j) = stx(j)
      sys%y(ibase + j) = sty(j)
      sys%z(ibase + j) = stz(j)
    end do

  end subroutine move_ccb

end module runt_moves
