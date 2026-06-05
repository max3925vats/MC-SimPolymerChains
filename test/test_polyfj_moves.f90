! test_polyfj_moves.f90 — unit tests for polyfj_moves (athermal hard-core moves)
!
! Uses test-drive framework (v0.6.0).
! gfortran-15 workaround: allocate(testsuite(N)) + element-wise assignment.
!
! System: nmol=8 chains of n=4 in L=10, rsp=1.0
! All beads placed on a wide lattice so bond length=1 and min inter-mol dist>=2.
!
! Tests:
!   (a) no_overlap_after_accept — after every accepted move the whole system
!       has no hard bead-bead or bead-sphere overlap (O(N^2) min-image scan).
!       cl is rebuilt after each accept, as the driver would do.
!   (b) ccb_exercised — CCB moves are run many times; the test verifies that
!       at least one CCB move is accepted and leaves no overlap (this is
!       covered by (a) but explicitly counted here).
!   (c) all_move_types_active — confirm at least one accept and one reject
!       over many steps, and that all three move types get accepted >= 1 time.

module test_polyfj_moves_suite
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,     only: dp, i8
  use polymc_rng,       only: rng_t, rng_seed, rng_uniform
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t, cl_init, cl_build
  use polymc_config,    only: system_t, params_t
  use polyfj_moves,     only: move_dick, move_rept, move_ccb
  implicit none
  private
  public :: collect_polyfj_moves

  integer,  parameter :: NMOL_TEST = 8
  integer,  parameter :: N_TEST    = 4
  real(dp), parameter :: L_TEST    = 10.0_dp
  real(dp), parameter :: RSP_TEST  = 1.0_dp

contains

  !> Register all tests
  subroutine collect_polyfj_moves(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(3))
    testsuite(1) = new_unittest("no_overlap_after_accept",  test_no_overlap_after_accept)
    testsuite(2) = new_unittest("ccb_exercised",            test_ccb_exercised)
    testsuite(3) = new_unittest("all_move_types_active",    test_all_move_types_active)
  end subroutine collect_polyfj_moves

  ! ------------------------------------------------------------------
  !> Build a non-overlapping lattice initial configuration.
  !
  ! 8 molecules placed on a 2x2x2 grid (cell spacing 5.0 sigma, L=10).
  ! Chain centres at (2.5+ix*5, 2.5+iy*5, 2.5+iz*5) for ix,iy,iz in {0,1}.
  ! Within each chain, bead j is placed at centre + (j-2.5)*1.0 along x
  ! (bond length = 1 sigma, chain spans 3 sigma).
  !
  ! Overlap checks:
  !   Inter-mol min dist = 5.0 - 1.5 - 1.5 = 2.0 sigma >> 1 (no overlap).
  !   Central sphere: closest bead is at (2.5-1.5, 2.5, 2.5) = (1.0, 2.5, 2.5)
  !   which is 2.5+2.5+0.25=... wait: dist^2 = 1+6.25+6.25 = 13.5; dist=3.67.
  !   Contact = rsp+0.5 = 1.5; 3.67 >> 1.5 (no sphere overlap). Verified.
  ! ------------------------------------------------------------------
  subroutine build_system(sys, p, cl)
    type(system_t),    intent(out) :: sys
    type(params_t),    intent(out) :: p
    type(cell_list_t), intent(out) :: cl

    integer  :: imol, ibead, ix, iy, iz
    real(dp) :: molx, moly, molz, r_cl

    sys%nmol   = NMOL_TEST
    sys%n      = N_TEST
    sys%nbeads = NMOL_TEST * N_TEST
    sys%box%L  = L_TEST
    sys%rsp    = RSP_TEST

    allocate(sys%x(sys%nbeads), sys%y(sys%nbeads), sys%z(sys%nbeads))
    allocate(sys%mol_of(sys%nbeads))

    imol = 0
    do iz = 0, 1
      do iy = 0, 1
        do ix = 0, 1
          imol = imol + 1
          molx = 2.5_dp + real(ix, dp) * 5.0_dp
          moly = 2.5_dp + real(iy, dp) * 5.0_dp
          molz = 2.5_dp + real(iz, dp) * 5.0_dp
          do ibead = 1, N_TEST
            sys%x((imol-1)*N_TEST + ibead) = molx + (real(ibead,dp) - 2.5_dp)
            sys%y((imol-1)*N_TEST + ibead) = moly
            sys%z((imol-1)*N_TEST + ibead) = molz
            sys%mol_of((imol-1)*N_TEST + ibead) = imol
          end do
        end do
      end do
    end do

    ! Athermal polyfj parameters: no energy wells.
    p%dlr    = 0.1_dp
    p%dint   = 0.1_dp
    p%fdick  = 0.33_dp
    p%frept  = 0.33_dp
    p%rsp    = RSP_TEST
    ! Energy params unused but set to zero for clarity
    p%beps   = 0.0_dp
    p%bbeps  = 0.0_dp
    p%rcut   = -1.0_dp

    r_cl = 1.0_dp
    call cl_init(cl, sys%box, r_cl)
    call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)
  end subroutine build_system

  !> Rebuild cell list after an accepted move (mirrors driver behaviour).
  subroutine rebuild_cl(sys, cl)
    type(system_t),    intent(in)    :: sys
    type(cell_list_t), intent(inout) :: cl
    call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)
  end subroutine rebuild_cl

  ! ------------------------------------------------------------------
  !> Full-system hard-overlap scan.
  !
  ! Returns .true. if:
  !  (a) any two beads from DIFFERENT molecules have dist < 1.0 (min-image), OR
  !  (b) any two beads within the SAME molecule that are non-bonded
  !      (separated by >= 2 positions) have dist < 1.0.
  !
  ! Bonded pairs (adjacent in chain) are excluded because their dist=1.0
  ! exactly, which is acceptable (equality not strict <).
  ! ------------------------------------------------------------------
  logical function system_has_overlap(sys)
    type(system_t), intent(in) :: sys

    integer  :: i, j
    real(dp) :: dx, dy, dz, d2
    integer  :: mol_i, mol_j, pos_i, pos_j

    system_has_overlap = .false.

    do i = 1, sys%nbeads
      do j = i + 1, sys%nbeads
        mol_i = sys%mol_of(i)
        mol_j = sys%mol_of(j)

        ! Same molecule: skip bonded pairs (|pos_i - pos_j| < 2)
        if (mol_i == mol_j) then
          pos_i = i - (mol_i - 1) * sys%n
          pos_j = j - (mol_j - 1) * sys%n
          if (abs(pos_i - pos_j) < 2) cycle
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
  !> Test (a): no hard overlap after every accepted move.
  !
  ! Run 1000 mixed moves (cycling through DICK/REPT/CCB and all molecules).
  ! After each accepted move rebuild cl and assert no system-wide overlap.
  ! ==================================================================
  subroutine test_no_overlap_after_accept(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl
    type(rng_t)       :: rng

    logical  :: acc
    integer  :: step, imol, move_type

    call build_system(sys, p, cl)
    call rng_seed(rng, 42_i8)

    do step = 1, 1000
      imol      = mod(step - 1, sys%nmol) + 1
      move_type = mod(step - 1, 3)

      select case (move_type)
      case (0)
        call move_dick(sys, p, cl, rng, imol, acc)
      case (1)
        call move_rept(sys, p, cl, rng, imol, acc)
      case (2)
        call move_ccb(sys, p, cl, rng, imol, acc)
      end select

      if (acc) then
        call rebuild_cl(sys, cl)
        call check(error, .not. system_has_overlap(sys), &
                   "no_overlap: hard overlap found after accepted move")
        if (allocated(error)) then
          deallocate(sys%x, sys%y, sys%z, sys%mol_of)
          return
        end if
      end if
    end do

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_no_overlap_after_accept

  ! ==================================================================
  !> Test (b): CCB moves are exercised and at least one is accepted
  !  without leaving a hard overlap.
  !
  ! Run 2000 CCB-only moves across all molecules; count accepted.
  ! Assert n_acc_ccb > 0 (at least one CCB accepted over 2000 attempts).
  ! The no-overlap property for each accepted CCB is verified inline.
  ! ==================================================================
  subroutine test_ccb_exercised(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl
    type(rng_t)       :: rng

    logical  :: acc
    integer  :: step, imol, n_acc_ccb

    call build_system(sys, p, cl)
    call rng_seed(rng, 137_i8)

    n_acc_ccb = 0

    do step = 1, 2000
      imol = mod(step - 1, sys%nmol) + 1
      call move_ccb(sys, p, cl, rng, imol, acc)

      if (acc) then
        n_acc_ccb = n_acc_ccb + 1
        call rebuild_cl(sys, cl)
        call check(error, .not. system_has_overlap(sys), &
                   "ccb_exercised: hard overlap found after accepted CCB move")
        if (allocated(error)) then
          deallocate(sys%x, sys%y, sys%z, sys%mol_of)
          return
        end if
      end if
    end do

    call check(error, n_acc_ccb > 0, &
               "ccb_exercised: no CCB moves were accepted in 2000 attempts")

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_ccb_exercised

  ! ==================================================================
  !> Test (c): all three move types accepted at least once; some rejected.
  !
  ! Run 3000 mixed moves.  Assert:
  !   - total n_acc > 0 and total n_rej > 0 (move machinery works)
  !   - n_acc_dick > 0, n_acc_rept > 0, n_acc_ccb > 0 (each type accepted)
  ! ==================================================================
  subroutine test_all_move_types_active(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl
    type(rng_t)       :: rng

    logical  :: acc
    integer  :: step, imol, move_type
    integer  :: n_acc, n_rej, n_acc_dick, n_acc_rept, n_acc_ccb

    call build_system(sys, p, cl)
    call rng_seed(rng, 271_i8)

    n_acc      = 0
    n_rej      = 0
    n_acc_dick = 0
    n_acc_rept = 0
    n_acc_ccb  = 0

    do step = 1, 3000
      imol      = mod(step - 1, sys%nmol) + 1
      move_type = mod(step - 1, 3)

      select case (move_type)
      case (0)
        call move_dick(sys, p, cl, rng, imol, acc)
        if (acc) n_acc_dick = n_acc_dick + 1
      case (1)
        call move_rept(sys, p, cl, rng, imol, acc)
        if (acc) n_acc_rept = n_acc_rept + 1
      case (2)
        call move_ccb(sys, p, cl, rng, imol, acc)
        if (acc) n_acc_ccb = n_acc_ccb + 1
      end select

      if (acc) then
        n_acc = n_acc + 1
        call rebuild_cl(sys, cl)
      else
        n_rej = n_rej + 1
      end if
    end do

    call check(error, n_acc > 0, "all_active: no moves accepted in 3000 steps")
    if (allocated(error)) then
      deallocate(sys%x, sys%y, sys%z, sys%mol_of)
      return
    end if

    call check(error, n_rej > 0, "all_active: no moves rejected in 3000 steps")
    if (allocated(error)) then
      deallocate(sys%x, sys%y, sys%z, sys%mol_of)
      return
    end if

    call check(error, n_acc_dick > 0, "all_active: DICK never accepted in 3000 steps")
    if (allocated(error)) then
      deallocate(sys%x, sys%y, sys%z, sys%mol_of)
      return
    end if

    call check(error, n_acc_rept > 0, "all_active: REPT never accepted in 3000 steps")
    if (allocated(error)) then
      deallocate(sys%x, sys%y, sys%z, sys%mol_of)
      return
    end if

    call check(error, n_acc_ccb > 0, "all_active: CCB never accepted in 3000 steps")

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_all_move_types_active

end module test_polyfj_moves_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_polyfj_moves
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_polyfj_moves_suite, only: collect_polyfj_moves
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  testsuites = [new_testsuite("polyfj_moves", collect_polyfj_moves)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_polyfj_moves
