! test_runt_potential.f90 — unit tests for runt_potential (Yukawa potentials)
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran-15 allocatable-finalizer bug workaround:
!   allocate(testsuite(N)) + element-wise assignment; no array ctor.
!
! Tests:
!   1. u_pair_contact        — at r=1, exp(0)=1, /1 → beps (exact)
!   2. u_pair_monotonic      — magnitude decays: |u(2)| < |u(1.5)| < |u(1)|
!   3. u_pair_cutoff         — beyond rcut → 0; inside rcut → equals exact
!   4. u_pair_exact_value    — u(2, -0.5, -1) ≈ -0.5*exp(-2.5)/2
!   5. u_sphere_contact      — at contact rsp+0.5=2, bbeps/r = -0.5/2
!   6. u_sphere_inside_core  — r < rsp+0.5 → 0 (hard core)
!   7. u_sphere_cutoff       — beyond rcut → 0

module test_runt_potential_suite
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,     only: dp
  use runt_potential,   only: u_pair, u_sphere
  implicit none
  private
  public :: collect_runt_potential

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_runt_potential(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(7))
    testsuite(1) = new_unittest("u_pair_contact",       test_u_pair_contact)
    testsuite(2) = new_unittest("u_pair_monotonic",     test_u_pair_monotonic)
    testsuite(3) = new_unittest("u_pair_cutoff",        test_u_pair_cutoff)
    testsuite(4) = new_unittest("u_pair_exact_value",   test_u_pair_exact_value)
    testsuite(5) = new_unittest("u_sphere_contact",     test_u_sphere_contact)
    testsuite(6) = new_unittest("u_sphere_inside_core", test_u_sphere_inside_core)
    testsuite(7) = new_unittest("u_sphere_cutoff",      test_u_sphere_cutoff)
  end subroutine collect_runt_potential

  ! ------------------------------------------------------------------
  ! u_pair at contact r=1: exp(-kappa*(1-1))/1 = 1/1 = 1 -> beps*1 = beps
  ! beps = -0.5, rcut = -1.0 (no truncation)
  ! ------------------------------------------------------------------
  subroutine test_u_pair_contact(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: u
    u = u_pair(1.0_dp, -0.5_dp, -1.0_dp)
    call check(error, abs(u - (-0.5_dp)) < 1.0e-12_dp, &
               "u_pair: at contact r=1 should equal beps=-0.5")
  end subroutine test_u_pair_contact

  ! ------------------------------------------------------------------
  ! Monotonic magnitude decay: |u(2)| < |u(1.5)| < |u(1)|
  ! (beps negative → u negative, magnitude is |u|)
  ! ------------------------------------------------------------------
  subroutine test_u_pair_monotonic(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: u1, u15, u2
    u1  = u_pair(1.0_dp, -0.5_dp, -1.0_dp)
    u15 = u_pair(1.5_dp, -0.5_dp, -1.0_dp)
    u2  = u_pair(2.0_dp, -0.5_dp, -1.0_dp)
    call check(error, abs(u2) < abs(u15), &
               "u_pair monotonic: |u(2)| should be < |u(1.5)|")
    if (allocated(error)) return
    call check(error, abs(u15) < abs(u1), &
               "u_pair monotonic: |u(1.5)| should be < |u(1)|")
  end subroutine test_u_pair_monotonic

  ! ------------------------------------------------------------------
  ! Cutoff tests: rcut=3.0
  !   r=3.1 > rcut → 0
  !   r=2.9 < rcut → equals the exact (no-cutoff) value
  ! ------------------------------------------------------------------
  subroutine test_u_pair_cutoff(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: u_beyond, u_inside, u_exact
    u_beyond = u_pair(3.1_dp, -0.5_dp, 3.0_dp)
    u_inside = u_pair(2.9_dp, -0.5_dp, 3.0_dp)
    u_exact  = u_pair(2.9_dp, -0.5_dp, -1.0_dp)
    call check(error, u_beyond == 0.0_dp, &
               "u_pair cutoff: r=3.1 > rcut=3.0 should return 0")
    if (allocated(error)) return
    call check(error, u_inside == u_exact, &
               "u_pair cutoff: r=2.9 < rcut=3.0 should equal exact value")
  end subroutine test_u_pair_cutoff

  ! ------------------------------------------------------------------
  ! Exact value: u_pair(2.0, -0.5, -1.0) = -0.5 * exp(-2.5*1.0) / 2.0
  ! r-1 = 1.0, so exp(-2.5*(r-1)) = exp(-2.5)
  ! ------------------------------------------------------------------
  subroutine test_u_pair_exact_value(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: u, u_expected
    u          = u_pair(2.0_dp, -0.5_dp, -1.0_dp)
    u_expected = -0.5_dp * exp(-2.5_dp * 1.0_dp) / 2.0_dp
    call check(error, abs(u - u_expected) < 1.0e-12_dp, &
               "u_pair exact: u(2.0,-0.5,-1) should be -0.5*exp(-2.5)/2")
  end subroutine test_u_pair_exact_value

  ! ------------------------------------------------------------------
  ! u_sphere at contact: rsp=1.5, r=rsp+0.5=2.0
  ! exp(-kappa*(2.0-2.0))/2.0 = 1/2.0 → bbeps/2.0 = -0.5/2.0 = -0.25
  ! ------------------------------------------------------------------
  subroutine test_u_sphere_contact(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: u, u_expected
    u          = u_sphere(2.0_dp, -0.5_dp, 1.5_dp, -1.0_dp)
    u_expected = -0.5_dp / 2.0_dp
    call check(error, abs(u - u_expected) < 1.0e-12_dp, &
               "u_sphere: at contact r=rsp+0.5=2.0, should be bbeps/r")
  end subroutine test_u_sphere_contact

  ! ------------------------------------------------------------------
  ! u_sphere inside hard core: r=1.9 < rsp+0.5=2.0 → should return 0
  ! ------------------------------------------------------------------
  subroutine test_u_sphere_inside_core(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: u
    u = u_sphere(1.9_dp, -0.5_dp, 1.5_dp, -1.0_dp)
    call check(error, u == 0.0_dp, &
               "u_sphere: r=1.9 inside hard core (rsp+0.5=2.0) should return 0")
  end subroutine test_u_sphere_inside_core

  ! ------------------------------------------------------------------
  ! u_sphere cutoff: rcut=5.0, r=5.1 > rcut → 0
  ! ------------------------------------------------------------------
  subroutine test_u_sphere_cutoff(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: u
    u = u_sphere(5.1_dp, -0.5_dp, 1.5_dp, 5.0_dp)
    call check(error, u == 0.0_dp, &
               "u_sphere cutoff: r=5.1 > rcut=5.0 should return 0")
  end subroutine test_u_sphere_cutoff

end module test_runt_potential_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_runt_potential
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_runt_potential_suite, only: collect_runt_potential
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element constructor is safe from the gfortran-15 finalizer bug.
  testsuites = [new_testsuite("runt_potential", collect_runt_potential)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_runt_potential
