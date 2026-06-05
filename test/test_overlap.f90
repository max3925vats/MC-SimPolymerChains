! test_overlap.f90 — unit tests for polymc_overlap
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran-15 finalizer bug: do NOT use array constructors for
! testsuite(:). Use allocate + element-by-element assignment instead.
!
! Tests:
!   1. pair_overlap      — query at dist 0.9 overlaps; at 1.1 does not.
!   2. exclude_self      — bead of same mol is excluded; of different mol counted.
!   3. matches_bruteforce — cell-list result matches full O(N) scan for ~400 beads.
!   4. sphere            — sphere_overlaps at rsp+0.4 → .true.; rsp+0.6 → .false.

module test_overlap_suite
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,     only: dp, i8
  use polymc_rng,       only: rng_t, rng_seed, rng_uniform
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t, cl_init, cl_build, cl_neighbours
  use polymc_overlap,   only: bead_overlaps, sphere_overlaps
  implicit none
  private
  public :: collect_overlap

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_overlap(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(4))
    testsuite(1) = new_unittest("pair_overlap",       test_pair_overlap)
    testsuite(2) = new_unittest("exclude_self",       test_exclude_self)
    testsuite(3) = new_unittest("matches_bruteforce", test_matches_bruteforce)
    testsuite(4) = new_unittest("sphere",             test_sphere)
  end subroutine collect_overlap

  ! ------------------------------------------------------------------
  ! One bead at known location (molecule 2).  Query from molecule 1.
  !   dist 0.9 centre-to-centre -> overlap (.true.)
  !   dist 1.1 centre-to-centre -> no overlap (.false.)
  ! ------------------------------------------------------------------
  subroutine test_pair_overlap(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t)       :: box
    type(cell_list_t) :: cl
    integer, parameter :: nb = 1
    real(dp) :: bx(nb), by(nb), bz(nb)
    integer  :: mol_of(nb)
    logical  :: res

    box%L = 10.0_dp
    call cl_init(cl, box, 1.5_dp)

    ! Bead of molecule 2 placed at (5.0, 5.0, 5.0)
    bx(1) = 5.0_dp;  by(1) = 5.0_dp;  bz(1) = 5.0_dp
    mol_of(1) = 2
    call cl_build(cl, box, bx, by, bz, nb)

    ! Query point at distance 0.9 along x from the bead (molecule 1)
    res = bead_overlaps(box, cl, bx, by, bz, nb, mol_of, &
                        5.9_dp, 5.0_dp, 5.0_dp, exclude_mol=1)
    call check(error, res, "dist 0.9: expected overlap (.true.)")
    if (allocated(error)) return

    ! Query point at distance 1.1 along x (no overlap expected)
    res = bead_overlaps(box, cl, bx, by, bz, nb, mol_of, &
                        6.1_dp, 5.0_dp, 5.0_dp, exclude_mol=1)
    call check(error, .not. res, "dist 1.1: expected no overlap (.false.)")
  end subroutine test_pair_overlap

  ! ------------------------------------------------------------------
  ! A bead of molecule 1 is at distance 0.5 from the query point.
  !   exclude_mol=1 -> same mol, skip -> no overlap (.false.)
  !   exclude_mol=2 -> different mol, count -> overlap (.true.)
  ! ------------------------------------------------------------------
  subroutine test_exclude_self(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t)       :: box
    type(cell_list_t) :: cl
    integer, parameter :: nb = 1
    real(dp) :: bx(nb), by(nb), bz(nb)
    integer  :: mol_of(nb)
    logical  :: res

    box%L = 10.0_dp
    call cl_init(cl, box, 1.5_dp)

    ! Bead of molecule 1 at (3.0, 3.0, 3.0)
    bx(1) = 3.0_dp;  by(1) = 3.0_dp;  bz(1) = 3.0_dp
    mol_of(1) = 1
    call cl_build(cl, box, bx, by, bz, nb)

    ! Query point at dist 0.5 along x — same molecule excluded
    res = bead_overlaps(box, cl, bx, by, bz, nb, mol_of, &
                        3.5_dp, 3.0_dp, 3.0_dp, exclude_mol=1)
    call check(error, .not. res, &
               "exclude_self: same mol at dist 0.5 should NOT count as overlap")
    if (allocated(error)) return

    ! Same geometry, but now exclude_mol=2 -> bead 1 (mol 1) IS counted
    res = bead_overlaps(box, cl, bx, by, bz, nb, mol_of, &
                        3.5_dp, 3.0_dp, 3.0_dp, exclude_mol=2)
    call check(error, res, &
               "exclude_self: different mol at dist 0.5 SHOULD give overlap")
  end subroutine test_exclude_self

  ! ------------------------------------------------------------------
  ! Brute-force oracle vs. bead_overlaps for ~400 beads, 200 queries
  ! ------------------------------------------------------------------
  subroutine test_matches_bruteforce(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t)       :: box
    type(cell_list_t) :: cl
    type(rng_t)       :: rng
    integer, parameter :: nb  = 400
    integer, parameter :: nq  = 200
    real(dp), parameter :: L  = 8.0_dp
    real(dp), parameter :: r_cut = 1.0_dp

    real(dp) :: bx(nb), by(nb), bz(nb)
    integer  :: mol_of(nb)
    real(dp) :: qx, qy, qz, dx, dy, dz, d2
    logical  :: cl_res, bf_res
    integer  :: k, q, excl

    box%L = L
    call rng_seed(rng, 77_i8)

    ! Random bead positions
    do k = 1, nb
      bx(k) = rng_uniform(rng) * L
      by(k) = rng_uniform(rng) * L
      bz(k) = rng_uniform(rng) * L
      ! Chains of 4: mol_of(k) = ceil(k/4)
      mol_of(k) = (k - 1) / 4 + 1
    end do

    call cl_init(cl, box, r_cut)
    call cl_build(cl, box, bx, by, bz, nb)

    excl = 3  ! exclude molecule 3 throughout

    do q = 1, nq
      qx = rng_uniform(rng) * L
      qy = rng_uniform(rng) * L
      qz = rng_uniform(rng) * L

      ! Cell-list result
      cl_res = bead_overlaps(box, cl, bx, by, bz, nb, mol_of, &
                              qx, qy, qz, exclude_mol=excl)

      ! Independent brute-force oracle: min-image dist < sigma(=1)
      bf_res = .false.
      do k = 1, nb
        if (mol_of(k) == excl) cycle
        dx = min_image(bx(k) - qx, L)
        dy = min_image(by(k) - qy, L)
        dz = min_image(bz(k) - qz, L)
        d2 = dx*dx + dy*dy + dz*dz
        if (d2 < 1.0_dp) then
          bf_res = .true.
          exit
        end if
      end do

      call check(error, cl_res .eqv. bf_res, &
                 "bead_overlaps disagrees with brute-force oracle")
      if (allocated(error)) return
    end do
  end subroutine test_matches_bruteforce

  ! ------------------------------------------------------------------
  ! sphere_overlaps: rsp=1.5
  !   Point along x at radius 1.9 (= rsp+0.4) -> .true.
  !   Point along x at radius 2.1 (= rsp+0.6) -> .false.
  ! Also test a near-boundary along y inside primary cell.
  ! ------------------------------------------------------------------
  subroutine test_sphere(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t) :: box
    real(dp), parameter :: rsp = 1.5_dp
    logical :: res

    box%L = 10.0_dp

    ! Point at (1.9, 0, 0) — within primary cell, dist from origin = 1.9 < 2.0
    res = sphere_overlaps(box, rsp, 1.9_dp, 0.0_dp, 0.0_dp)
    call check(error, res, &
               "sphere: r=1.9 < rsp+0.5=2.0 should overlap (.true.)")
    if (allocated(error)) return

    ! Point at (2.1, 0, 0) — dist from origin = 2.1 > 2.0
    res = sphere_overlaps(box, rsp, 2.1_dp, 0.0_dp, 0.0_dp)
    call check(error, .not. res, &
               "sphere: r=2.1 > rsp+0.5=2.0 should not overlap (.false.)")
    if (allocated(error)) return

    ! Test min-image path: point placed near the far edge of the box
    ! (9.0, 0, 0) in L=10 box; min-image distance from origin = min(9, 1) = -1
    ! |min_image| = 1.0 < 2.0 -> overlap
    res = sphere_overlaps(box, rsp, 9.0_dp, 0.0_dp, 0.0_dp)
    call check(error, res, &
               "sphere min-image: x=9 in L=10 maps to dist 1.0 < 2.0 (.true.)")
  end subroutine test_sphere

end module test_overlap_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_overlap
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_overlap_suite, only: collect_overlap
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element constructor is safe from the gfortran-15 finalizer bug.
  testsuites = [new_testsuite("overlap", collect_overlap)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_overlap
