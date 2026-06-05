! test_box.f90 — unit tests for polymc_box (box type, min_image, cell_index)
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran 15 finalizer bug: do NOT use array constructors for
! testsuite(:). Use allocate + element-by-element assignment instead.

module test_box_suite
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use polymc_kinds, only: dp
  use polymc_box,   only: box_t, min_image, cell_index
  implicit none
  private
  public :: collect_box

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_box(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(5))
    testsuite(1) = new_unittest("min_image_pos",    test_min_image_pos)
    testsuite(2) = new_unittest("min_image_neg",    test_min_image_neg)
    testsuite(3) = new_unittest("min_image_zero",   test_min_image_zero)
    testsuite(4) = new_unittest("cell_index_basic", test_cell_index_basic)
    testsuite(5) = new_unittest("cell_index_clamp", test_cell_index_clamp)
  end subroutine collect_box

  ! ------------------------------------------------------------------
  ! min_image: d=0.6, L=1.0  ->  -0.4
  ! ------------------------------------------------------------------
  subroutine test_min_image_pos(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: m
    m = min_image(0.6_dp, 1.0_dp)
    call check(error, abs(m - (-0.4_dp)) < 1.0e-12_dp, &
               "min_image(0.6, 1.0) should be -0.4")
  end subroutine test_min_image_pos

  ! ------------------------------------------------------------------
  ! min_image: d=-0.6, L=1.0  ->  +0.4
  ! ------------------------------------------------------------------
  subroutine test_min_image_neg(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: m
    m = min_image(-0.6_dp, 1.0_dp)
    call check(error, abs(m - 0.4_dp) < 1.0e-12_dp, &
               "min_image(-0.6, 1.0) should be +0.4")
  end subroutine test_min_image_neg

  ! ------------------------------------------------------------------
  ! min_image: d=0.0, L=10.0  ->  0.0
  ! ------------------------------------------------------------------
  subroutine test_min_image_zero(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: m
    m = min_image(0.0_dp, 10.0_dp)
    call check(error, m == 0.0_dp, &
               "min_image(0.0, 10.0) should be exactly 0.0")
  end subroutine test_min_image_zero

  ! ------------------------------------------------------------------
  ! cell_index basic: L=10, mc=5
  !   (0.0,0.0,0.0)   -> ix=iy=iz=1
  !   (9.99,5.0,2.5)  -> ix=5, iy in 1..5, iz in 1..5
  ! ------------------------------------------------------------------
  subroutine test_cell_index_basic(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t) :: box
    integer :: ix, iy, iz

    box%L = 10.0_dp

    ! Origin maps to cell (1,1,1)
    call cell_index(box, 5, 0.0_dp, 0.0_dp, 0.0_dp, ix, iy, iz)
    call check(error, ix == 1, "origin: ix should be 1")
    if (allocated(error)) return
    call check(error, iy == 1, "origin: iy should be 1")
    if (allocated(error)) return
    call check(error, iz == 1, "origin: iz should be 1")
    if (allocated(error)) return

    ! x=9.99 is just below L, should fall in cell 5 (last)
    call cell_index(box, 5, 9.99_dp, 5.0_dp, 2.5_dp, ix, iy, iz)
    call check(error, ix == 5, "x=9.99: ix should be 5")
    if (allocated(error)) return
    call check(error, iy >= 1 .and. iy <= 5, "y=5.0: iy out of 1..5")
    if (allocated(error)) return
    call check(error, iz >= 1 .and. iz <= 5, "z=2.5: iz out of 1..5")
    if (allocated(error)) return
  end subroutine test_cell_index_basic

  ! ------------------------------------------------------------------
  ! cell_index clamp: every index stays in 1..mc across a range of coords
  ! including values at 0.0, near L, and mid-box.
  ! ------------------------------------------------------------------
  subroutine test_cell_index_clamp(error)
    type(error_type), allocatable, intent(out) :: error
    type(box_t) :: box
    integer :: ix, iy, iz, i
    real(dp) :: coords(6)
    ! Test positions: 0, near-zero, quarter, half, three-quarter, just-below-L
    coords = [0.0_dp, 0.01_dp, 2.5_dp, 5.0_dp, 7.5_dp, 9.9999_dp]
    box%L  = 10.0_dp

    do i = 1, size(coords)
      call cell_index(box, 5, coords(i), coords(i), coords(i), ix, iy, iz)
      call check(error, ix >= 1 .and. ix <= 5, "ix out of range")
      if (allocated(error)) return
      call check(error, iy >= 1 .and. iy <= 5, "iy out of range")
      if (allocated(error)) return
      call check(error, iz >= 1 .and. iz <= 5, "iz out of range")
      if (allocated(error)) return
    end do
  end subroutine test_cell_index_clamp

end module test_box_suite

! ==================================================================
!  Driver program
! ==================================================================
program test_box
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_box_suite, only: collect_box
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element constructor is safe from the gfortran-15 finalizer bug.
  testsuites = [new_testsuite("box", collect_box)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_box
