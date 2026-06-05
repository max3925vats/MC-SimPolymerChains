! test_runt_moves.f90 — unit tests for runt_moves (move_dick, move_rept, move_ccb)
!
! Uses test-drive framework (v0.6.0).
! NOTE: gfortran-15 allocatable-finalizer workaround:
!   allocate(testsuite(N)) + element-wise assignment; no array constructor.
!
! System: nmol=8 chains of n=4 in L=8.0, rsp=1.0
! Non-overlapping lattice initial config ensures no hard overlaps at start.
!
! Tests:
!   (a) energy_bookkeeping — track erun via de; periodically compare to
!       total_energy recompute; tolerance 1e-6.
!   (b) no_overlap_after_accept — every accepted move leaves no bead-bead
!       hard overlaps in the whole system.
!   (c) athermal_limit — with beps=bbeps=0, every accepted move has de==0
!       and no overlaps; some moves accepted, some rejected.

module test_runt_moves_suite
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,     only: dp, i8
  use polymc_rng,       only: rng_t, rng_seed, rng_uniform
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t, cl_init, cl_build
  use polymc_config,    only: system_t, params_t
  use runt_energy,      only: total_energy
  use runt_moves,       only: move_dick, move_rept, move_ccb
  implicit none
  private
  public :: collect_runt_moves

  ! Box / chain parameters used throughout
  integer,  parameter :: NMOL_TEST  = 8
  integer,  parameter :: N_TEST     = 4
  real(dp), parameter :: L_TEST     = 10.0_dp
  real(dp), parameter :: RSP_TEST   = 1.0_dp

contains

  !> Register all tests
  subroutine collect_runt_moves(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(3))
    testsuite(1) = new_unittest("energy_bookkeeping",    test_energy_bookkeeping)
    testsuite(2) = new_unittest("no_overlap_after_accept", test_no_overlap_after_accept)
    testsuite(3) = new_unittest("athermal_limit",         test_athermal_limit)
  end subroutine collect_runt_moves

  ! ------------------------------------------------------------------
  !> Build a non-overlapping lattice initial configuration.
  !
  ! nmol=8 chains of n=4 in a cubic box of L=10.
  ! Molecules are placed on a 2x2x2 lattice (spacing 5.0).
  ! Within each molecule, beads are separated by 1.0 sigma along x
  ! (unit bond length, matching what MC moves expect).
  ! Chain spans 3 sigma (beads at -1.5, -0.5, +0.5, +1.5 from center).
  ! Min inter-molecular bead distance = 5.0 - 1.5 - 1.5 = 2.0 sigma >> 1.
  ! Min distance from origin (sphere center): closest mol center is 2.5,
  ! closest bead = 2.5 - 1.5 = 1.0 sigma from center.
  ! Sphere contact = rsp + 0.5 = 1.5; bead at 1.0 < 1.5 would overlap sphere.
  ! So offset mol centers to (3.5, ...) so min bead dist to origin = 2.0 > 1.5.
  ! ------------------------------------------------------------------
  subroutine build_lattice_system(sys, p, cl)
    type(system_t),    intent(out) :: sys
    type(params_t),    intent(out) :: p
    type(cell_list_t), intent(out) :: cl

    integer  :: imol, ibead, ix, iy, iz
    real(dp) :: molx, moly, molz
    real(dp) :: r_cl

    sys%nmol   = NMOL_TEST
    sys%n      = N_TEST
    sys%nbeads = NMOL_TEST * N_TEST
    sys%box%L  = L_TEST
    sys%rsp    = RSP_TEST

    allocate(sys%x(sys%nbeads), sys%y(sys%nbeads), sys%z(sys%nbeads))
    allocate(sys%mol_of(sys%nbeads))

    ! Place 8 molecules on a 2x2x2 grid with cell spacing 5.0, L=10.
    ! Centres at (2.5+ix*5, 2.5+iy*5, 2.5+iz*5) for ix,iy,iz in {0,1}.
    ! Bead j (j=1..4) placed along x at centre + (j-2.5)*1.0; y=z=centre.
    ! Bond length = 1.0 sigma.
    ! Min bead distance from sphere centre (origin): min-image of bead position;
    !   closest bead is at x=2.5-1.5=1.0 with y,z=2.5 -> dist = sqrt(1+6.25+6.25)=3.67 > 1.5 ✓
    ! Min inter-mol bead dist: 5.0 - 1.5 - 1.5 = 2.0 sigma >> 1 ✓
    imol = 0
    do iz = 0, 1
      do iy = 0, 1
        do ix = 0, 1
          imol  = imol + 1
          molx  = 2.5_dp + real(ix, dp) * 5.0_dp
          moly  = 2.5_dp + real(iy, dp) * 5.0_dp
          molz  = 2.5_dp + real(iz, dp) * 5.0_dp

          ! Beads along x with unit bond length (1.0 sigma)
          do ibead = 1, N_TEST
            sys%x((imol-1)*N_TEST + ibead) = molx + (real(ibead,dp) - 2.5_dp) * 1.0_dp
            sys%y((imol-1)*N_TEST + ibead) = moly
            sys%z((imol-1)*N_TEST + ibead) = molz
            sys%mol_of((imol-1)*N_TEST + ibead) = imol
          end do
        end do
      end do
    end do

    ! Parameters: attractive interactions, exact mode
    p%dlr   = 0.5_dp
    p%dint  = 0.5_dp
    p%beps  = -0.5_dp
    p%bbeps = -0.5_dp
    p%rsp   = RSP_TEST
    p%rcut  = -1.0_dp   ! exact (no cutoff)

    ! Build cell list (r_cut = 1.0 is the hard-core distance; use that)
    r_cl = 1.0_dp
    call cl_init(cl, sys%box, r_cl)
    call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

  end subroutine build_lattice_system

  ! ------------------------------------------------------------------
  !> Rebuild cell list after an accepted move (matches driver behavior).
  ! ------------------------------------------------------------------
  subroutine rebuild_cl(sys, cl)
    type(system_t),    intent(in)    :: sys
    type(cell_list_t), intent(inout) :: cl
    call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)
  end subroutine rebuild_cl

  ! ------------------------------------------------------------------
  !> Check for bead-bead hard overlaps in the whole system.
  !  Returns .true. if any pair of beads from DIFFERENT molecules are
  !  within sigma=1 of each other (min-image), OR if any intra non-bonded
  !  pair (|i-j|>=2 within same molecule) has dist < 1.
  !
  !  Bonded beads (adjacent in chain) can be at distance 1 by construction,
  !  so we only check non-bonded pairs (separated by >= 2 in chain, or
  !  different molecules).
  ! ------------------------------------------------------------------
  logical function system_has_overlap(sys)
    type(system_t), intent(in) :: sys

    integer  :: i, j
    real(dp) :: dx, dy, dz, d2
    integer  :: mol_i, mol_j, pos_in_mol_i, pos_in_mol_j

    system_has_overlap = .false.

    do i = 1, sys%nbeads
      do j = i + 1, sys%nbeads
        mol_i = sys%mol_of(i)
        mol_j = sys%mol_of(j)

        ! For same molecule: only check non-bonded pairs (skip bonded: |diff|<2)
        if (mol_i == mol_j) then
          pos_in_mol_i = i - (mol_i - 1) * sys%n
          pos_in_mol_j = j - (mol_j - 1) * sys%n
          if (abs(pos_in_mol_i - pos_in_mol_j) < 2) cycle   ! bonded: skip
        end if

        dx = min_image(sys%x(j) - sys%x(i), sys%box%L)
        dy = min_image(sys%y(j) - sys%y(i), sys%box%L)
        dz = min_image(sys%z(j) - sys%z(i), sys%box%L)
        d2 = dx*dx + dy*dy + dz*dz

        if (d2 < 1.0_dp) then
          system_has_overlap = .true.
          return
        end if
      end do
    end do
  end function system_has_overlap

  ! ==================================================================
  !> Test (a): energy bookkeeping.
  !
  ! Run 2000 mixed moves, tracking erun by accumulating de.
  ! Every 200 moves assert |erun - total_energy(...)| < 1e-6.
  ! Move cycle: imol = mod(step-1, nmol) + 1; move_type = mod(step-1, 3)
  !   0 -> move_dick, 1 -> move_rept, 2 -> move_ccb
  ! ==================================================================
  subroutine test_energy_bookkeeping(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl
    type(rng_t)       :: rng

    real(dp)  :: erun, etot, de
    logical   :: acc
    integer   :: step, imol, move_type
    real(dp), parameter :: tol = 1.0e-6_dp

    call build_lattice_system(sys, p, cl)
    call rng_seed(rng, 42_i8)

    ! Initial energy
    erun = total_energy(sys, p, cl)

    do step = 1, 2000
      imol      = mod(step - 1, sys%nmol) + 1
      move_type = mod(step - 1, 3)

      select case (move_type)
      case (0)
        call move_dick(sys, p, cl, rng, imol, acc, de)
      case (1)
        call move_rept(sys, p, cl, rng, imol, acc, de)
      case (2)
        call move_ccb(sys, p, cl, rng, imol, acc, de)
      end select

      if (acc) then
        erun = erun + de
        call rebuild_cl(sys, cl)
      end if

      ! Every 200 steps compare running energy to recompute
      if (mod(step, 200) == 0) then
        etot = total_energy(sys, p, cl)
        call check(error, abs(erun - etot) < tol, &
                   "energy bookkeeping: erun drifts from total_energy recompute")
        if (allocated(error)) then
          deallocate(sys%x, sys%y, sys%z, sys%mol_of)
          return
        end if
      end if
    end do

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_energy_bookkeeping

  ! ==================================================================
  !> Test (b): no bead-bead hard overlaps after every accepted move,
  !> AND sys is bitwise unchanged after every rejected move.
  !
  ! Run 5000 mixed moves (cycling imol; many rept+ccb so the 50% chain
  ! reversal fires often).  After each ACCEPTED move, rebuild cl (as the
  ! driver does) and assert the whole system has no hard overlap.
  !
  ! After each REJECTED move, assert sys%x/y/z are BITWISE identical to the
  ! pre-move snapshot.  This directly guards the reversal-on-reject bug:
  ! a rejected chain reversal must not mutate sys (which would also leave
  ! the cell list stale and corrupt subsequent overlap detection).
  ! ==================================================================
  subroutine test_no_overlap_after_accept(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl
    type(rng_t)       :: rng

    logical  :: acc
    real(dp) :: de
    integer  :: step, imol, move_type
    real(dp), allocatable :: x0(:), y0(:), z0(:)   ! pre-move snapshot

    call build_lattice_system(sys, p, cl)
    call rng_seed(rng, 137_i8)

    allocate(x0(sys%nbeads), y0(sys%nbeads), z0(sys%nbeads))

    do step = 1, 5000
      imol      = mod(step - 1, sys%nmol) + 1
      ! Weight the move mix toward rept+ccb so reversals fire often:
      ! cycle dick, rept, ccb, rept, ccb (3 of 5 moves are reversal-capable).
      select case (mod(step - 1, 5))
      case (0)
        move_type = 0
      case (1, 3)
        move_type = 1
      case default
        move_type = 2
      end select

      ! Snapshot sys before the move (to verify reject leaves it unchanged).
      x0 = sys%x
      y0 = sys%y
      z0 = sys%z

      select case (move_type)
      case (0)
        call move_dick(sys, p, cl, rng, imol, acc, de)
      case (1)
        call move_rept(sys, p, cl, rng, imol, acc, de)
      case (2)
        call move_ccb(sys, p, cl, rng, imol, acc, de)
      end select

      if (acc) then
        call rebuild_cl(sys, cl)
        call check(error, .not. system_has_overlap(sys), &
                   "no_overlap: hard overlap found after accepted move")
        if (allocated(error)) then
          deallocate(x0, y0, z0)
          deallocate(sys%x, sys%y, sys%z, sys%mol_of)
          return
        end if
      else
        ! Rejected move: sys must be BITWISE unchanged (no premature reversal).
        call check(error, all(sys%x == x0) .and. all(sys%y == y0) &
                          .and. all(sys%z == z0), &
                   "reject: sys positions mutated by a rejected move")
        if (allocated(error)) then
          deallocate(x0, y0, z0)
          deallocate(sys%x, sys%y, sys%z, sys%mol_of)
          return
        end if
      end if
    end do

    deallocate(x0, y0, z0)
    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_no_overlap_after_accept

  ! ==================================================================
  !> Test (c): athermal limit (beps=bbeps=0).
  !
  ! Every accepted move must have de==0 (within floating-point tolerance)
  ! and leave no hard overlaps.  Run 500 moves; assert at least 1 accepted
  ! and at least 1 rejected (sanity that the move machinery works).
  ! ==================================================================
  subroutine test_athermal_limit(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl
    type(rng_t)       :: rng

    logical  :: acc
    real(dp) :: de
    integer  :: step, imol, move_type, n_acc, n_rej
    real(dp), parameter :: tol_de = 1.0e-12_dp

    call build_lattice_system(sys, p, cl)

    ! Override to zero interactions
    p%beps  = 0.0_dp
    p%bbeps = 0.0_dp

    call rng_seed(rng, 271_i8)

    n_acc = 0
    n_rej = 0

    do step = 1, 500
      imol      = mod(step - 1, sys%nmol) + 1
      move_type = mod(step - 1, 3)

      select case (move_type)
      case (0)
        call move_dick(sys, p, cl, rng, imol, acc, de)
      case (1)
        call move_rept(sys, p, cl, rng, imol, acc, de)
      case (2)
        call move_ccb(sys, p, cl, rng, imol, acc, de)
      end select

      if (acc) then
        n_acc = n_acc + 1
        call rebuild_cl(sys, cl)

        ! de must be 0 in athermal limit
        call check(error, abs(de) < tol_de, &
                   "athermal: accepted move has non-zero de")
        if (allocated(error)) then
          deallocate(sys%x, sys%y, sys%z, sys%mol_of)
          return
        end if

        ! No overlaps after accept
        call check(error, .not. system_has_overlap(sys), &
                   "athermal: hard overlap found after accepted move")
        if (allocated(error)) then
          deallocate(sys%x, sys%y, sys%z, sys%mol_of)
          return
        end if
      else
        n_rej = n_rej + 1
      end if
    end do

    ! Sanity: some moves were accepted and some were rejected
    call check(error, n_acc > 0, "athermal: no moves were accepted (500 steps)")
    if (allocated(error)) then
      deallocate(sys%x, sys%y, sys%z, sys%mol_of)
      return
    end if

    call check(error, n_rej > 0, "athermal: no moves were rejected (500 steps)")

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_athermal_limit

end module test_runt_moves_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_runt_moves
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_runt_moves_suite, only: collect_runt_moves
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  testsuites = [new_testsuite("runt_moves", collect_runt_moves)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_runt_moves
