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
  use polymc_chain,     only: chain_reverse
  use polymc_config,    only: system_t, params_t
  use runt_potential,   only: u_pair, u_sphere
  use runt_energy,      only: bead_inter_energy
  implicit none
  private

  public :: move_dick, move_rept, move_ccb

  ! Number of Rosenbluth trials per bead in CCB (matches legacy NSAMP=15)
  integer, parameter :: NSAMP = 15

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
      if (rng_uniform(rng) <= 0.25_dp) then
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
  ! After accept, the chain_reverse may have changed the bead order in sys;
  ! the caller rebuilds cl from the updated sys.
  ! ====================================================================
  subroutine move_rept(sys, p, cl, rng, imol, accepted, de)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol
    logical,           intent(out)   :: accepted
    real(dp),          intent(out)   :: de

    integer  :: ibase, ii, i1, i2
    real(dp) :: vx, vy, vz
    real(dp) :: tnx, tny, tnz   ! new head bead position
    real(dp) :: rdis, r
    real(dp) :: enew, eold

    accepted = .false.
    de       = 0.0_dp
    ibase    = (imol - 1) * sys%n

    ! ------------------------------------------------------------------
    ! Step 1: randomly reverse chain (with prob 0.5)
    ! NOTE: chain_reverse modifies sys%x/y/z permanently for this molecule.
    ! This is faithful to legacy REPT which calls CVERT unconditionally on
    ! the sys arrays.  The driver rebuilds cl after accept (which includes
    ! the reversal).  On reject the reversal stays in sys; this is also
    ! faithful to legacy (CVERT is not undone on reject).
    ! ------------------------------------------------------------------
    if (rng_uniform(rng) > 0.5_dp) then
      call chain_reverse(sys%x, sys%y, sys%z, ibase, sys%n)
    end if

    ! ------------------------------------------------------------------
    ! Step 2: grow new head bead at pos(1) + unit_vector
    ! Legacy: XITR(1) = X(IMOL+1) + DELX/AL  where |DELX/AL| = 1/AL = 1 sigma
    ! Modern: directly add unit vector (length 1 sigma)
    ! ------------------------------------------------------------------
    call rng_unit_vector(rng, vx, vy, vz)
    tnx = sys%x(ibase + 1) + vx
    tny = sys%y(ibase + 1) + vy
    tnz = sys%z(ibase + 1) + vz

    ! Sphere overlap check (hard reject)
    if (sphere_overlaps(sys%box, sys%rsp, tnx, tny, tnz)) return

    ! Inter-molecular bead overlap check (hard reject)
    if (bead_overlaps(sys%box, cl, sys%x, sys%y, sys%z, sys%nbeads, &
                      sys%mol_of, tnx, tny, tnz, imol)) return

    ! Accumulate ENEW for new head
    r    = dist_to_origin(tnx, tny, tnz, sys%box%L)
    enew = bead_inter_energy(sys, p, cl, imol, tnx, tny, tnz) &
         + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

    ! Intra overlap check: new head vs beads 2..N-1 (old beads that survive)
    ! and accumulate intra energy with beads 2..N-1
    ! Legacy DO 320 II=2,N-1 uses IMOL+II (0-based offset gives beads 2..N-1)
    do ii = 2, sys%n - 1
      i1   = ibase + ii
      rdis = mi_dist2(tnx, tny, tnz, sys%x(i1), sys%y(i1), sys%z(i1), &
                      sys%box%L)
      if (rdis < 1.0_dp) return
      enew = enew + u_pair(sqrt(rdis), p%beps, p%rcut)
    end do

    ! ------------------------------------------------------------------
    ! Step 4: compute EOLD for the removed tail bead (bead N)
    ! Legacy REPT: intra pairs of bead N vs beads 1..N-2, inter, sphere
    ! ------------------------------------------------------------------
    eold = 0.0_dp
    i2   = ibase + sys%n   ! global index of bead N (tail to be removed)

    ! Intra: bead N vs beads 1..N-2
    do ii = 1, sys%n - 2
      i1   = ibase + ii
      rdis = mi_dist2(sys%x(i1), sys%y(i1), sys%z(i1), &
                      sys%x(i2), sys%y(i2), sys%z(i2), sys%box%L)
      eold = eold + u_pair(sqrt(rdis), p%beps, p%rcut)
    end do

    ! Inter + sphere for removed bead N
    r    = dist_to_origin(sys%x(i2), sys%y(i2), sys%z(i2), sys%box%L)
    eold = eold + bead_inter_energy(sys, p, cl, imol, &
                    sys%x(i2), sys%y(i2), sys%z(i2)) &
                + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

    ! ------------------------------------------------------------------
    ! Metropolis acceptance
    ! ------------------------------------------------------------------
    if (eold <= enew) then
      if (exp(eold - enew) <= rng_uniform(rng)) return
    end if

    ! ------------------------------------------------------------------
    ! Move accepted: new chain = [new_head, old_1, old_2, ..., old_{N-1}]
    ! Legacy reset loop: XITR(J)=X(IMOL+J-1) for J=2..N, then writes XITR to X
    ! ------------------------------------------------------------------
    accepted = .true.
    de       = enew - eold

    ! Shift old beads 1..N-1 into positions 2..N (copy backwards to be safe)
    ! Then write new head into position 1.
    do ii = sys%n, 2, -1
      sys%x(ibase + ii) = sys%x(ibase + ii - 1)
      sys%y(ibase + ii) = sys%y(ibase + ii - 1)
      sys%z(ibase + ii) = sys%z(ibase + ii - 1)
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
  ! ====================================================================
  subroutine move_ccb(sys, p, cl, rng, imol, accepted, de)
    type(system_t),    intent(inout) :: sys
    type(params_t),    intent(in)    :: p
    type(cell_list_t), intent(in)    :: cl
    type(rng_t),       intent(inout) :: rng
    integer,           intent(in)    :: imol
    logical,           intent(out)   :: accepted
    real(dp),          intent(out)   :: de

    integer  :: ibase, j, k, ii, icut, i1, i2
    real(dp) :: xn(sys%n), yn(sys%n), zn(sys%n)   ! working/trial chain
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
    ! Step 1: randomly reverse chain (with prob 0.5)
    ! ------------------------------------------------------------------
    if (rng_uniform(rng) > 0.5_dp) then
      call chain_reverse(sys%x, sys%y, sys%z, ibase, sys%n)
    end if

    ! ------------------------------------------------------------------
    ! Step 2: copy current chain into working array XN
    ! ------------------------------------------------------------------
    do j = 1, sys%n
      i1    = ibase + j
      xn(j) = sys%x(i1)
      yn(j) = sys%y(i1)
      zn(j) = sys%z(i1)
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
        if (j >= 3) then
          do ii = 1, j - 2
            rdis = mi_dist2(xt(k), yt(k), zt(k), xn(ii), yn(ii), zn(ii), &
                            sys%box%L)
            if (rdis < 1.0_dp) then
              et(k) = 0.0_dp
              goto 19  ! jump to wsum accumulation
            end if
            eni = eni + u_pair(sqrt(rdis), p%beps, p%rcut)
          end do
        end if

        et(k) = exp(-eni)
19      continue
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

      ! Accumulate Rosenbluth weight WN (selected normalised weight)
      wn    = wn * et(k)
      ! Update working chain
      xn(j) = xt(k)
      yn(j) = yt(k)
      zn(j) = zt(k)

    end do  ! j = icut..n  (NEW chain)

    ! ------------------------------------------------------------------
    ! Step 5: store new chain; reload XN from original sys positions
    ! ------------------------------------------------------------------
    do j = 1, sys%n
      stx(j) = xn(j)
      sty(j) = yn(j)
      stz(j) = zn(j)
    end do

    do j = 1, sys%n
      i1    = ibase + j
      xn(j) = sys%x(i1)
      yn(j) = sys%y(i1)
      zn(j) = sys%z(i1)
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
          ! Use existing bead position (already verified non-overlapping)
          i1  = ibase + j
          r   = dist_to_origin(sys%x(i1), sys%y(i1), sys%z(i1), sys%box%L)
          eno = bead_inter_energy(sys, p, cl, imol, &
                                  sys%x(i1), sys%y(i1), sys%z(i1)) &
              + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

          ! Intra contribution for existing bead
          if (j >= 3) then
            do ii = 1, j - 2
              i2   = ibase + ii
              rdis = mi_dist2(sys%x(i1), sys%y(i1), sys%z(i1), &
                              sys%x(i2), sys%y(i2), sys%z(i2), sys%box%L)
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
        if (j >= 3) then
          do ii = 1, j - 2
            rdis = mi_dist2(xt(k), yt(k), zt(k), xn(ii), yn(ii), zn(ii), &
                            sys%box%L)
            if (rdis < 1.0_dp) then
              et(k) = 0.0_dp
              goto 79
            end if
            eno = eno + u_pair(sqrt(rdis), p%beps, p%rcut)
          end do
        end if

        et(k) = exp(-eno)
79      continue
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
      i1 = ibase + j

      ! NEW chain: inter + sphere for bead J
      r    = dist_to_origin(stx(j), sty(j), stz(j), sys%box%L)
      enew = enew + bead_inter_energy(sys, p, cl, imol, stx(j), sty(j), stz(j)) &
                  + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

      ! OLD chain: inter + sphere for bead J
      r    = dist_to_origin(sys%x(i1), sys%y(i1), sys%z(i1), sys%box%L)
      eold = eold + bead_inter_energy(sys, p, cl, imol, &
                      sys%x(i1), sys%y(i1), sys%z(i1)) &
                  + u_sphere(r, p%bbeps, sys%rsp, p%rcut)

      ! Intra pairs: J vs beads 1..J-2
      if (j >= 3) then
        do k = 1, j - 2
          ! New chain intra
          rdis = mi_dist2(stx(j), sty(j), stz(j), stx(k), sty(k), stz(k), &
                          sys%box%L)
          enew = enew + u_pair(sqrt(rdis), p%beps, p%rcut)

          ! Old chain intra
          i2   = ibase + k
          rdis = mi_dist2(sys%x(i1), sys%y(i1), sys%z(i1), &
                          sys%x(i2), sys%y(i2), sys%z(i2), sys%box%L)
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
