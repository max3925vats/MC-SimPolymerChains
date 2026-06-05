! test_chain.f90 — unit tests for polymc_chain
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran-15 finalizer bug: do NOT use array constructors for
! testsuite(:). Use allocate + element-by-element assignment instead.
!
! Tests:
!   1. reverse_even   — chain of n=4 beads; reversed coords match expected.
!   2. untouched      — chain 1 (indices 1..4) unchanged after reversing chain 2.
!   3. reverse_odd    — chain of n=3; middle bead stays in place.
!   4. reverse_single — n=1; no-op (no crash, coord unchanged).
!   5. permutation    — reversal is a permutation (set of coords preserved).

module test_chain_suite
  use testdrive,     only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,  only: dp
  use polymc_chain,  only: chain_reverse
  implicit none
  private
  public :: collect_chain

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_chain(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(5))
    testsuite(1) = new_unittest("reverse_even",   test_reverse_even)
    testsuite(2) = new_unittest("untouched",       test_untouched)
    testsuite(3) = new_unittest("reverse_odd",     test_reverse_odd)
    testsuite(4) = new_unittest("reverse_single",  test_reverse_single)
    testsuite(5) = new_unittest("permutation",     test_permutation)
  end subroutine collect_chain

  ! ------------------------------------------------------------------
  ! Chain 2 (ibase=4, n=4): beads at (1,0,0),(2,0,0),(3,0,0),(4,0,0).
  ! After reverse: x(5..8) == [4,3,2,1]; y,z all still 0.
  ! ------------------------------------------------------------------
  subroutine test_reverse_even(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x(8), y(8), z(8)
    integer :: j

    x = 0.0_dp;  y = 0.0_dp;  z = 0.0_dp
    ! Chain 2: indices 5..8
    do j = 1, 4
      x(4+j) = real(j, dp)   ! x = 1,2,3,4
    end do

    call chain_reverse(x, y, z, 4, 4)

    ! Expected: x(5..8) = 4,3,2,1
    call check(error, x(5) == 4.0_dp, "reverse_even: x(5) should be 4")
    if (allocated(error)) return
    call check(error, x(6) == 3.0_dp, "reverse_even: x(6) should be 3")
    if (allocated(error)) return
    call check(error, x(7) == 2.0_dp, "reverse_even: x(7) should be 2")
    if (allocated(error)) return
    call check(error, x(8) == 1.0_dp, "reverse_even: x(8) should be 1")
    if (allocated(error)) return

    ! y and z all zero throughout
    do j = 5, 8
      call check(error, y(j) == 0.0_dp, "reverse_even: y should stay 0")
      if (allocated(error)) return
      call check(error, z(j) == 0.0_dp, "reverse_even: z should stay 0")
      if (allocated(error)) return
    end do
  end subroutine test_reverse_even

  ! ------------------------------------------------------------------
  ! Chain 1 (ibase=0, n=4) untouched when we reverse chain 2.
  ! ------------------------------------------------------------------
  subroutine test_untouched(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x(8), y(8), z(8)
    integer :: j

    x = 0.0_dp;  y = 0.0_dp;  z = 0.0_dp
    ! Chain 1: indices 1..4 — distinct values so any clobber is detectable
    do j = 1, 4
      x(j) = real(j, dp) * 10.0_dp
      y(j) = real(j, dp) * 20.0_dp
      z(j) = real(j, dp) * 30.0_dp
    end do
    ! Chain 2: indices 5..8
    do j = 1, 4
      x(4+j) = real(j, dp)
    end do

    call chain_reverse(x, y, z, 4, 4)  ! reverse chain 2 only

    ! Chain 1 must be untouched
    do j = 1, 4
      call check(error, x(j) == real(j, dp)*10.0_dp, "untouched: x chain 1 changed")
      if (allocated(error)) return
      call check(error, y(j) == real(j, dp)*20.0_dp, "untouched: y chain 1 changed")
      if (allocated(error)) return
      call check(error, z(j) == real(j, dp)*30.0_dp, "untouched: z chain 1 changed")
      if (allocated(error)) return
    end do
  end subroutine test_untouched

  ! ------------------------------------------------------------------
  ! Odd chain n=3: positions a,b,c -> c,b,a; middle bead unchanged.
  ! ------------------------------------------------------------------
  subroutine test_reverse_odd(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x(3), y(3), z(3)

    x(1) = 1.0_dp;  x(2) = 5.0_dp;  x(3) = 9.0_dp
    y(1) = 2.0_dp;  y(2) = 6.0_dp;  y(3) = 10.0_dp
    z(1) = 3.0_dp;  z(2) = 7.0_dp;  z(3) = 11.0_dp

    call chain_reverse(x, y, z, 0, 3)

    ! First and last swapped
    call check(error, x(1) == 9.0_dp,  "reverse_odd: x(1) should be 9")
    if (allocated(error)) return
    call check(error, x(3) == 1.0_dp,  "reverse_odd: x(3) should be 1")
    if (allocated(error)) return
    ! Middle unchanged
    call check(error, x(2) == 5.0_dp,  "reverse_odd: x(2) (middle) should stay 5")
    if (allocated(error)) return
    call check(error, y(2) == 6.0_dp,  "reverse_odd: y(2) (middle) should stay 6")
    if (allocated(error)) return
    call check(error, z(2) == 7.0_dp,  "reverse_odd: z(2) (middle) should stay 7")
    if (allocated(error)) return
    ! y, z endpoints also swapped
    call check(error, y(1) == 10.0_dp, "reverse_odd: y(1) should be 10")
    if (allocated(error)) return
    call check(error, z(1) == 11.0_dp, "reverse_odd: z(1) should be 11")
  end subroutine test_reverse_odd

  ! ------------------------------------------------------------------
  ! n=1: single bead — must be a no-op (no crash, value unchanged).
  ! ------------------------------------------------------------------
  subroutine test_reverse_single(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x(1), y(1), z(1)

    x(1) = 42.0_dp;  y(1) = -1.0_dp;  z(1) = 3.14_dp

    call chain_reverse(x, y, z, 0, 1)

    call check(error, x(1) == 42.0_dp,  "reverse_single: x should be unchanged")
    if (allocated(error)) return
    call check(error, y(1) == -1.0_dp,  "reverse_single: y should be unchanged")
    if (allocated(error)) return
    call check(error, z(1) == 3.14_dp,  "reverse_single: z should be unchanged")
  end subroutine test_reverse_single

  ! ------------------------------------------------------------------
  ! The set (multiset) of coordinates is preserved after reversal.
  ! Use n=5 with distinct x-values; sort-and-compare via sum/product.
  ! ------------------------------------------------------------------
  subroutine test_permutation(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x(5), y(5), z(5)
    real(dp) :: sum_before, prod_before, sum_after, prod_after
    integer  :: j

    do j = 1, 5
      x(j) = real(j, dp) * 1.7_dp + 0.3_dp
      y(j) = real(j, dp) * 0.5_dp
      z(j) = real(j, dp) * 2.3_dp - 1.0_dp
    end do

    sum_before  = sum(x) + sum(y) + sum(z)
    prod_before = product(x)

    call chain_reverse(x, y, z, 0, 5)

    sum_after  = sum(x) + sum(y) + sum(z)
    prod_after = product(x)

    call check(error, abs(sum_after  - sum_before)  < 1.0e-12_dp, &
               "permutation: sum changed after reversal")
    if (allocated(error)) return
    call check(error, abs(prod_after - prod_before) < 1.0e-10_dp, &
               "permutation: product changed after reversal")
  end subroutine test_permutation

end module test_chain_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_chain
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_chain_suite, only: collect_chain
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  testsuites = [new_testsuite("chain", collect_chain)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_chain
