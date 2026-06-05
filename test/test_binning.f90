! test_binning.f90 — unit tests for polymc_binning (radial histogram + shell-volume density)
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran-15 finalizer bug: do NOT use array constructors for
! testsuite(:) or unittest(:). Use allocate + element-by-element assignment.

module test_binning_suite
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,   only: dp, i8
  use polymc_rng,     only: rng_t, rng_seed, rng_uniform
  use polymc_binning, only: hist_t, hist_init, hist_add, hist_density
  implicit none
  private
  public :: collect_binning

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_binning(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(3))
    testsuite(1) = new_unittest("counts",          test_counts)
    testsuite(2) = new_unittest("edge",            test_edge)
    testsuite(3) = new_unittest("uniform_density", test_uniform_density)
  end subroutine collect_binning

  ! ------------------------------------------------------------------
  ! counts: bsz=0.5, rmax=5.0 => nbin=10
  !   r=0.1  => k=int(0.1/0.5)+1 = 0+1 = 1  (bin 1)
  !   r=0.6  => k=int(0.6/0.5)+1 = 1+1 = 2  (bin 2)
  !   r=2.4  => k=int(2.4/0.5)+1 = 4+1 = 5  (bin 5)
  !   r=4.9  => k=int(4.9/0.5)+1 = 9+1 = 10 (bin 10)
  !   r=5.5  => k=int(5.5/0.5)+1 = 11+1= 12 (ignored, > nbin)
  ! ------------------------------------------------------------------
  subroutine test_counts(error)
    type(error_type), allocatable, intent(out) :: error
    type(hist_t) :: h

    call hist_init(h, 0.5_dp, 5.0_dp)

    ! Verify nbin
    call check(error, h%nbin == 10, "nbin should be 10 for rmax=5.0, bsz=0.5")
    if (allocated(error)) return

    ! Verify all counts start at zero
    call check(error, all(h%count == 0_i8), "all counts should start at 0")
    if (allocated(error)) return

    ! Add points
    call hist_add(h, 0.1_dp)   ! -> bin 1
    call hist_add(h, 0.6_dp)   ! -> bin 2
    call hist_add(h, 2.4_dp)   ! -> bin 5
    call hist_add(h, 4.9_dp)   ! -> bin 10
    call hist_add(h, 5.5_dp)   ! -> ignored (out of range)

    call check(error, h%count(1)  == 1_i8, "bin 1 should have count 1")
    if (allocated(error)) return
    call check(error, h%count(2)  == 1_i8, "bin 2 should have count 1")
    if (allocated(error)) return
    call check(error, h%count(5)  == 1_i8, "bin 5 should have count 1")
    if (allocated(error)) return
    call check(error, h%count(10) == 1_i8, "bin 10 should have count 1")
    if (allocated(error)) return

    ! Bins 3,4,6,7,8,9 should remain zero
    call check(error, h%count(3) == 0_i8, "bin 3 should be 0")
    if (allocated(error)) return
    call check(error, h%count(9) == 0_i8, "bin 9 should be 0")
    if (allocated(error)) return

    ! Total: exactly 4 counts accumulated
    call check(error, sum(h%count) == 4_i8, "total count should be 4")
  end subroutine test_counts

  ! ------------------------------------------------------------------
  ! edge: r=0.0 lands in bin 1; r<0 is ignored
  ! ------------------------------------------------------------------
  subroutine test_edge(error)
    type(error_type), allocatable, intent(out) :: error
    type(hist_t) :: h

    call hist_init(h, 0.5_dp, 5.0_dp)

    ! r=0.0 exactly: int(0.0/0.5)+1 = 0+1 = 1 => bin 1
    call hist_add(h, 0.0_dp)
    call check(error, h%count(1) == 1_i8, "r=0.0 should land in bin 1")
    if (allocated(error)) return

    ! r=-0.1: k = int(-0.1/0.5)+1; negative, should be ignored
    call hist_add(h, -0.1_dp)
    call check(error, sum(h%count) == 1_i8, "r<0 should be ignored, total stays 1")
  end subroutine test_edge

  ! ------------------------------------------------------------------
  ! uniform_density: generate 2e6 points uniformly in a sphere of
  ! radius rmax, call hist_density, verify the estimated number density
  ! is approximately constant across mid-bins (3..8) — i.e. the
  ! shell-volume normalization has removed the r^2 growth.
  !
  ! Strategy: uniform sampling in a sphere via rejection (u,v,w uniform
  ! in [-1,1]^3, accept if u^2+v^2+w^2 < rmax^2). Each accepted point
  ! contributes its radius r = sqrt(u^2+v^2+w^2) to the histogram.
  ! True number density rho = N_accepted / V_sphere is constant.
  ! The variation in bins 3..8 should be within 10% of the mean.
  ! ------------------------------------------------------------------
  subroutine test_uniform_density(error)
    type(error_type), allocatable, intent(out) :: error

    type(hist_t) :: h
    type(rng_t)  :: rng

    real(dp), parameter :: bsz  = 0.25_dp
    real(dp), parameter :: rmax = 5.0_dp
    integer,  parameter :: ntry = 4000000   ! rejection loop draws; ~52% acceptance
    integer             :: i, nsamp
    real(dp) :: u, v, w, r2, rmax2
    real(dp), allocatable :: dist(:), dens(:)
    real(dp) :: dens_mid(6), mean_mid, dev

    call rng_seed(rng, 42_i8)
    call hist_init(h, bsz, rmax)
    rmax2 = rmax * rmax
    nsamp = 0

    ! Generate uniform points inside the sphere
    do i = 1, ntry
      u = (rng_uniform(rng) - 0.5_dp) * 2.0_dp * rmax
      v = (rng_uniform(rng) - 0.5_dp) * 2.0_dp * rmax
      w = (rng_uniform(rng) - 0.5_dp) * 2.0_dp * rmax
      r2 = u*u + v*v + w*w
      if (r2 < rmax2) then
        call hist_add(h, sqrt(r2))
        nsamp = nsamp + 1
      end if
    end do

    allocate(dist(h%nbin), dens(h%nbin))
    call hist_density(h, nsamp, dist, dens)

    ! Check that mid-bins 3..8 have approximately constant density
    ! (within 10% of their own mean)
    dens_mid = dens(3:8)
    mean_mid = sum(dens_mid) / 6.0_dp

    call check(error, mean_mid > 0.0_dp, "mean mid-bin density should be positive")
    if (allocated(error)) return

    do i = 1, 6
      dev = abs(dens_mid(i) - mean_mid) / mean_mid
      call check(error, dev < 0.10_dp, &
        "mid-bin density deviates more than 10% from mean — r^2 normalization may be wrong")
      if (allocated(error)) return
    end do

  end subroutine test_uniform_density

end module test_binning_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_binning
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_binning_suite, only: collect_binning
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element array constructor is safe from the gfortran-15 finalizer bug.
  testsuites = [new_testsuite("binning", collect_binning)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_binning
