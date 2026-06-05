! test_cell_list.f90 — unit tests for polymc_cell_list
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran 15 finalizer bug: do NOT use array constructors for
! testsuite(:). Use allocate + element-by-element assignment instead.
!
! Strategy:
!   1. Generate a random config (L=8, r_cut=1, nbeads=500).
!   2. cl_init, cl_build.
!   3. For 200 random query points:
!        a. cl_neighbours -> candidates; keep those within r_cut (min-image).
!        b. Brute-force: scan all beads for those within r_cut.
!        c. Sort both id lists; assert identical.
!   4. Assert every returned id is in 1..nbeads.
!   5. Assert mc >= 3 after cl_init.

module test_cell_list_suite
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,     only: dp, i8
  use polymc_rng,       only: rng_t, rng_seed, rng_uniform
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t, cl_init, cl_build, cl_neighbours
  implicit none
  private
  public :: collect_cell_list

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_cell_list(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(3))
    testsuite(1) = new_unittest("mc_at_least_3",      test_mc_at_least_3)
    testsuite(2) = new_unittest("ids_in_range",       test_ids_in_range)
    testsuite(3) = new_unittest("matches_bruteforce", test_matches_bruteforce)
  end subroutine collect_cell_list

  ! ------------------------------------------------------------------
  ! After cl_init with reasonable L and r_cut, mc >= 3
  ! ------------------------------------------------------------------
  subroutine test_mc_at_least_3(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t)       :: box
    type(cell_list_t) :: cl
    box%L = 8.0_dp
    call cl_init(cl, box, 1.0_dp)
    call check(error, cl%mc >= 3, "cl_init: mc should be >= 3")
  end subroutine test_mc_at_least_3

  ! ------------------------------------------------------------------
  ! Every bead id returned by cl_neighbours is in 1..nbeads
  ! ------------------------------------------------------------------
  subroutine test_ids_in_range(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t)       :: box
    type(cell_list_t) :: cl
    type(rng_t)       :: rng
    integer, parameter :: nbeads = 500
    real(dp) :: bx(nbeads), by(nbeads), bz(nbeads)
    real(dp) :: qx, qy, qz
    integer  :: idx(nbeads), nidx, k, q
    integer, parameter :: nq = 200

    box%L = 8.0_dp
    call rng_seed(rng, 11_i8)
    ! Random bead positions in [0, L)
    do k = 1, nbeads
      bx(k) = rng_uniform(rng) * box%L
      by(k) = rng_uniform(rng) * box%L
      bz(k) = rng_uniform(rng) * box%L
    end do

    call cl_init(cl, box, 1.0_dp)
    call cl_build(cl, box, bx, by, bz, nbeads)

    do q = 1, nq
      qx = rng_uniform(rng) * box%L
      qy = rng_uniform(rng) * box%L
      qz = rng_uniform(rng) * box%L
      call cl_neighbours(cl, box, qx, qy, qz, idx, nidx)
      do k = 1, nidx
        call check(error, idx(k) >= 1 .and. idx(k) <= nbeads, &
                   "cl_neighbours returned id outside 1..nbeads")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_ids_in_range

  ! ------------------------------------------------------------------
  ! cl_neighbours result matches brute-force for 200 random query points
  ! ------------------------------------------------------------------
  subroutine test_matches_bruteforce(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t)       :: box
    type(cell_list_t) :: cl
    type(rng_t)       :: rng
    integer, parameter :: nbeads = 500
    integer, parameter :: nq     = 200
    real(dp), parameter :: r_cut = 1.0_dp
    real(dp) :: bx(nbeads), by(nbeads), bz(nbeads)
    real(dp) :: qx, qy, qz, dx, dy, dz, dist2
    integer  :: cl_ids(nbeads), bf_ids(nbeads)
    integer  :: ncl, nbf, k, q

    box%L = 8.0_dp
    call rng_seed(rng, 42_i8)

    ! Random bead positions in [0, L)
    do k = 1, nbeads
      bx(k) = rng_uniform(rng) * box%L
      by(k) = rng_uniform(rng) * box%L
      bz(k) = rng_uniform(rng) * box%L
    end do

    call cl_init(cl, box, r_cut)
    call cl_build(cl, box, bx, by, bz, nbeads)

    do q = 1, nq
      qx = rng_uniform(rng) * box%L
      qy = rng_uniform(rng) * box%L
      qz = rng_uniform(rng) * box%L

      ! --- cell-list path: gather candidates, distance-filter ---
      call cl_neighbours(cl, box, qx, qy, qz, cl_ids, ncl)
      ! Filter to those within r_cut (in-place compact into cl_ids)
      nbf = 0  ! temporarily reuse for counting
      do k = 1, ncl
        dx = min_image(bx(cl_ids(k)) - qx, box%L)
        dy = min_image(by(cl_ids(k)) - qy, box%L)
        dz = min_image(bz(cl_ids(k)) - qz, box%L)
        dist2 = dx*dx + dy*dy + dz*dz
        if (dist2 <= r_cut*r_cut) then
          nbf = nbf + 1
          cl_ids(nbf) = cl_ids(k)  ! compact in-place
        end if
      end do
      ncl = nbf

      ! --- brute-force path ---
      nbf = 0
      do k = 1, nbeads
        dx = min_image(bx(k) - qx, box%L)
        dy = min_image(by(k) - qy, box%L)
        dz = min_image(bz(k) - qz, box%L)
        dist2 = dx*dx + dy*dy + dz*dz
        if (dist2 <= r_cut*r_cut) then
          nbf = nbf + 1
          bf_ids(nbf) = k
        end if
      end do

      ! Sort both and compare
      call isort(cl_ids, ncl)
      call isort(bf_ids, nbf)

      call check(error, ncl == nbf, &
                 "cl_neighbours count does not match brute-force")
      if (allocated(error)) return

      do k = 1, ncl
        call check(error, cl_ids(k) == bf_ids(k), &
                   "cl_neighbours id list does not match brute-force")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_matches_bruteforce

  ! ------------------------------------------------------------------
  ! Simple insertion sort for small integer arrays
  ! ------------------------------------------------------------------
  subroutine isort(a, n)
    integer, intent(inout) :: a(:)
    integer, intent(in)    :: n
    integer :: i, j, tmp
    do i = 2, n
      tmp = a(i)
      j = i - 1
      do while (j >= 1 .and. a(j) > tmp)
        a(j+1) = a(j)
        j = j - 1
      end do
      a(j+1) = tmp
    end do
  end subroutine isort

end module test_cell_list_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_cell_list
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_cell_list_suite, only: collect_cell_list
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element constructor is safe from the gfortran-15 finalizer bug.
  testsuites = [new_testsuite("cell_list", collect_cell_list)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_cell_list
