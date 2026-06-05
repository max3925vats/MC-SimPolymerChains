! test_rng.f90 — unit tests for polymc_rng (xoshiro256** generator)
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran 15 has a bug where finalizers on derived types with
! allocatable character components crash when those types appear as
! temporaries in array constructors ([a, b, c]).  Both unittest_type and
! testsuite_type carry such finalizers.  Workaround: allocate the output
! array explicitly and assign elements one by one.

module test_rng_suite
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use polymc_kinds, only: dp, i8
  use polymc_rng, only: rng_t, rng_seed, rng_uniform, rng_unit_vector
  implicit none
  private
  public :: collect_rng

contains

  !> Register all tests in this suite (gfortran-15-safe: no array ctor)
  subroutine collect_rng(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(4))
    testsuite(1) = new_unittest("range",        test_range)
    testsuite(2) = new_unittest("mean",         test_mean)
    testsuite(3) = new_unittest("unitvec",      test_unitvec)
    testsuite(4) = new_unittest("reproducible", test_reproducible)
  end subroutine collect_rng

  ! ------------------------------------------------------------------
  ! range: 100 000 draws all satisfy 0 <= u < 1
  ! ------------------------------------------------------------------
  subroutine test_range(error)
    type(error_type), allocatable, intent(out) :: error
    type(rng_t) :: rng
    real(dp)    :: u
    integer     :: i
    call rng_seed(rng, 42_i8)
    do i = 1, 100000
      u = rng_uniform(rng)
      call check(error, u >= 0.0_dp, "rng_uniform returned negative value")
      if (allocated(error)) return
      call check(error, u < 1.0_dp, "rng_uniform returned value >= 1")
      if (allocated(error)) return
    end do
  end subroutine test_range

  ! ------------------------------------------------------------------
  ! mean: mean of 1 000 000 draws within 1e-3 of 0.5
  ! ------------------------------------------------------------------
  subroutine test_mean(error)
    type(error_type), allocatable, intent(out) :: error
    type(rng_t) :: rng
    real(dp)    :: s
    integer     :: i
    integer, parameter :: n = 1000000
    call rng_seed(rng, 123_i8)
    s = 0.0_dp
    do i = 1, n
      s = s + rng_uniform(rng)
    end do
    s = s / real(n, dp)
    call check(error, abs(s - 0.5_dp) < 1.0e-3_dp, &
               "mean of rng_uniform deviates from 0.5 by more than 1e-3")
  end subroutine test_mean

  ! ------------------------------------------------------------------
  ! unitvec: 100 000 vectors all have |v| within 1e-12 of 1;
  !          mean x-component within 1e-2 of 0 (isotropy)
  ! ------------------------------------------------------------------
  subroutine test_unitvec(error)
    type(error_type), allocatable, intent(out) :: error
    type(rng_t) :: rng
    real(dp)    :: vx, vy, vz, len, sx
    integer     :: i
    integer, parameter :: n = 100000
    call rng_seed(rng, 7_i8)
    sx = 0.0_dp
    do i = 1, n
      call rng_unit_vector(rng, vx, vy, vz)
      len = sqrt(vx*vx + vy*vy + vz*vz)
      call check(error, abs(len - 1.0_dp) < 1.0e-12_dp, &
                 "rng_unit_vector length deviates from 1 by more than 1e-12")
      if (allocated(error)) return
      sx = sx + vx
    end do
    sx = sx / real(n, dp)
    call check(error, abs(sx) < 1.0e-2_dp, &
               "mean x-component of unit vectors deviates from 0 by more than 1e-2")
  end subroutine test_unitvec

  ! ------------------------------------------------------------------
  ! reproducible: two RNGs seeded identically produce identical first draw
  ! ------------------------------------------------------------------
  subroutine test_reproducible(error)
    type(error_type), allocatable, intent(out) :: error
    type(rng_t) :: rng1, rng2
    real(dp)    :: u1, u2
    call rng_seed(rng1, 999_i8)
    call rng_seed(rng2, 999_i8)
    u1 = rng_uniform(rng1)
    u2 = rng_uniform(rng2)
    call check(error, u1 == u2, &
               "two identically-seeded RNGs produced different first draws")
  end subroutine test_reproducible

end module test_rng_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_rng
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_rng_suite, only: collect_rng
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element constructor is safe from the gfortran-15 finalizer bug.
  testsuites = [new_testsuite("rng", collect_rng)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_rng
