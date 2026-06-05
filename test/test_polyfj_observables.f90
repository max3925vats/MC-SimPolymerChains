! test_polyfj_observables.f90 — unit tests for polyfj_observables
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran-15 allocatable-finalizer bug workaround:
!   allocate(testsuite(N)) + element-wise assignment; no array ctor.
!
! Tests:
!   1. density_counts    — 1 chain of 4 beads; verify h_all=4, h_end=2, h_mid=2
!   2. rg2_straight      — straight chain n=3 along x; Rg^2 = 2/3 ≈ 0.6667
!   3. re2_known         — same chain; Re^2 = 4.0
!   4. shape_known       — same chain; a≈sqrt(2/3), b≈0, c≈0
!   5. segorder_rod      — chain along +x far from origin; <cos>≈1

module test_polyfj_observables_suite
  use testdrive,              only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,           only: dp, i8
  use polymc_config,          only: system_t
  use polymc_box,             only: box_t
  use polyfj_observables,     only: polyfj_obs_t, obs_init, obs_sample
  implicit none
  private
  public :: collect_polyfj_obs

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_polyfj_obs(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(5))
    testsuite(1) = new_unittest("density_counts",  test_density_counts)
    testsuite(2) = new_unittest("rg2_straight",    test_rg2_straight)
    testsuite(3) = new_unittest("re2_known",       test_re2_known)
    testsuite(4) = new_unittest("shape_known",     test_shape_known)
    testsuite(5) = new_unittest("segorder_rod",    test_segorder_rod)
  end subroutine collect_polyfj_obs

  ! ------------------------------------------------------------------
  ! Helper: build a minimal system_t with one chain of n beads.
  !
  ! Positions (real sigma units) are supplied in xs/ys/zs(1..n).
  ! Box side length L is chosen large enough that min-image of all
  ! coordinates equals the coordinate itself (|coord| << L/2).
  ! ------------------------------------------------------------------
  subroutine make_system(sys, n, xs, ys, zs, L)
    type(system_t), intent(out)  :: sys
    integer,        intent(in)   :: n
    real(dp),       intent(in)   :: xs(n), ys(n), zs(n)
    real(dp),       intent(in)   :: L

    integer :: j

    sys%nmol   = 1
    sys%n      = n
    sys%nbeads = n
    sys%box%L  = L
    sys%rsp    = 0.0_dp

    allocate(sys%x(n), sys%y(n), sys%z(n), sys%mol_of(n))
    do j = 1, n
      sys%x(j) = xs(j)
      sys%y(j) = ys(j)
      sys%z(j) = zs(j)
      sys%mol_of(j) = 1
    end do
  end subroutine make_system

  ! ------------------------------------------------------------------
  ! 1. density_counts
  !
  ! Chain of n=4 beads.  Positions chosen so that:
  !   - all 4 beads fall in distinct bins far from overflow
  !   - j=1 and j=4 are end beads (h_end)
  !   - j=2 and j=3 are middle beads (for n=4: nm1=2, nm2=3) (h_mid)
  !
  ! Positions (placed along z-axis at distances 0.25, 0.75, 1.25, 1.75
  ! with bsz=0.5): bins 1, 2, 3, 4 respectively.
  !   bead 1 (j=1): dist=0.25 → bin 1  (end)
  !   bead 2 (j=2): dist=0.75 → bin 2  (mid, nm1=2 for n=4)
  !   bead 3 (j=3): dist=1.25 → bin 3  (mid, nm2=3 for n=4)
  !   bead 4 (j=4): dist=1.75 → bin 4  (end)
  !
  ! L=20 so min_image(x,L)=x for all coords here.
  ! ------------------------------------------------------------------
  subroutine test_density_counts(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)      :: sys
    type(polyfj_obs_t)  :: obs
    integer, parameter  :: n = 4
    real(dp), parameter :: bsz  = 0.5_dp
    real(dp), parameter :: bszc = 0.1_dp
    real(dp), parameter :: L    = 20.0_dp
    real(dp), parameter :: rmax = L / 2.0_dp

    real(dp) :: xs(n), ys(n), zs(n)
    integer(i8) :: total_all

    ! Place beads at distances 0.25, 0.75, 1.25, 1.75 along z
    xs = 0.0_dp;  ys = 0.0_dp
    zs(1) = 0.25_dp   ! bin 1, end
    zs(2) = 0.75_dp   ! bin 2, mid (nm1=2)
    zs(3) = 1.25_dp   ! bin 3, mid (nm2=3)
    zs(4) = 1.75_dp   ! bin 4, end

    call make_system(sys, n, xs, ys, zs, L)
    call obs_init(obs, bsz, bszc, rmax)
    call obs_sample(obs, sys)

    ! Total count in h_all must equal n=4
    total_all = sum(obs%h_all%count)
    call check(error, total_all == 4_i8, "density_counts: h_all total /= 4")
    if (allocated(error)) return

    ! Each bead should be in a distinct bin
    call check(error, obs%h_all%count(1) == 1_i8, "density_counts: h_all(1) /= 1")
    if (allocated(error)) return
    call check(error, obs%h_all%count(2) == 1_i8, "density_counts: h_all(2) /= 1")
    if (allocated(error)) return
    call check(error, obs%h_all%count(3) == 1_i8, "density_counts: h_all(3) /= 1")
    if (allocated(error)) return
    call check(error, obs%h_all%count(4) == 1_i8, "density_counts: h_all(4) /= 1")
    if (allocated(error)) return

    ! h_end: j=1 (bin 1) and j=4 (bin 4)
    call check(error, sum(obs%h_end%count) == 2_i8, "density_counts: h_end total /= 2")
    if (allocated(error)) return
    call check(error, obs%h_end%count(1) == 1_i8, "density_counts: h_end(1) /= 1")
    if (allocated(error)) return
    call check(error, obs%h_end%count(4) == 1_i8, "density_counts: h_end(4) /= 1")
    if (allocated(error)) return

    ! h_mid: j=2 (nm1=2, bin 2) and j=3 (nm2=3, bin 3)
    call check(error, sum(obs%h_mid%count) == 2_i8, "density_counts: h_mid total /= 2")
    if (allocated(error)) return
    call check(error, obs%h_mid%count(2) == 1_i8, "density_counts: h_mid(2) /= 1")
    if (allocated(error)) return
    call check(error, obs%h_mid%count(3) == 1_i8, "density_counts: h_mid(3) /= 1")
    if (allocated(error)) return

  end subroutine test_density_counts

  ! ------------------------------------------------------------------
  ! 2. rg2_straight
  !
  ! Chain n=3 at (0,0,0),(1,0,0),(2,0,0).
  ! CoM = (1,0,0).  Rg^2 = (1/3)*[(-1)^2 + 0^2 + 1^2] = 2/3 ≈ 0.6667.
  !
  ! Translate so CoM is away from origin a bit (add offset to all):
  ! beads at (3,0,0),(4,0,0),(5,0,0) → CoM=(4,0,0), Rg^2 unchanged = 2/3.
  !
  ! bszc=0.5 → CoM distance xcmb = min_image(4,L)
  ! Use L=20 so min_image(4,20) = 4; bin = 1+int(4/0.5) = 9.
  ! rg2_sum(9)/nc(9) should ≈ 2/3.
  ! ------------------------------------------------------------------
  subroutine test_rg2_straight(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)      :: sys
    type(polyfj_obs_t)  :: obs
    integer, parameter  :: n = 3
    real(dp), parameter :: bsz  = 0.5_dp
    real(dp), parameter :: bszc = 0.5_dp
    real(dp), parameter :: L    = 20.0_dp
    real(dp), parameter :: rmax = L / 2.0_dp
    real(dp), parameter :: tol  = 1.0e-12_dp
    real(dp), parameter :: rg2_exact = 2.0_dp / 3.0_dp

    real(dp) :: xs(n), ys(n), zs(n)
    real(dp) :: com_dist, rg2_mean
    integer  :: kbin

    xs(1) = 3.0_dp;  ys(1) = 0.0_dp;  zs(1) = 0.0_dp
    xs(2) = 4.0_dp;  ys(2) = 0.0_dp;  zs(2) = 0.0_dp
    xs(3) = 5.0_dp;  ys(3) = 0.0_dp;  zs(3) = 0.0_dp

    call make_system(sys, n, xs, ys, zs, L)
    call obs_init(obs, bsz, bszc, rmax)
    call obs_sample(obs, sys)

    ! CoM = (4,0,0) → distance 4.0; min_image(4, 20)=4
    ! bin = 1 + int(4.0 / 0.5) = 1 + 8 = 9
    com_dist = 4.0_dp
    kbin = 1 + int(com_dist / bszc)

    call check(error, kbin <= obs%nbinc, "rg2_straight: kbin out of range")
    if (allocated(error)) return

    call check(error, obs%nc(kbin) == 1_i8, "rg2_straight: nc(kbin) /= 1")
    if (allocated(error)) return

    rg2_mean = obs%rg2_sum(kbin) / real(obs%nc(kbin), dp)
    call check(error, abs(rg2_mean - rg2_exact) < tol, &
               "rg2_straight: Rg^2 /= 2/3")
    if (allocated(error)) return

  end subroutine test_rg2_straight

  ! ------------------------------------------------------------------
  ! 3. re2_known
  !
  ! Same chain as rg2_straight: (3,0,0),(4,0,0),(5,0,0).
  ! Re^2 = (5-3)^2 = 4.0.
  ! ------------------------------------------------------------------
  subroutine test_re2_known(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)      :: sys
    type(polyfj_obs_t)  :: obs
    integer, parameter  :: n = 3
    real(dp), parameter :: bsz  = 0.5_dp
    real(dp), parameter :: bszc = 0.5_dp
    real(dp), parameter :: L    = 20.0_dp
    real(dp), parameter :: rmax = L / 2.0_dp
    real(dp), parameter :: tol  = 1.0e-12_dp

    real(dp) :: xs(n), ys(n), zs(n)
    real(dp) :: com_dist, re2_mean
    integer  :: kbin

    xs(1) = 3.0_dp;  ys(1) = 0.0_dp;  zs(1) = 0.0_dp
    xs(2) = 4.0_dp;  ys(2) = 0.0_dp;  zs(2) = 0.0_dp
    xs(3) = 5.0_dp;  ys(3) = 0.0_dp;  zs(3) = 0.0_dp

    call make_system(sys, n, xs, ys, zs, L)
    call obs_init(obs, bsz, bszc, rmax)
    call obs_sample(obs, sys)

    com_dist = 4.0_dp
    kbin = 1 + int(com_dist / bszc)

    call check(error, kbin <= obs%nbinc, "re2_known: kbin out of range")
    if (allocated(error)) return

    call check(error, obs%nc(kbin) == 1_i8, "re2_known: nc(kbin) /= 1")
    if (allocated(error)) return

    re2_mean = obs%re2_sum(kbin) / real(obs%nc(kbin), dp)
    call check(error, abs(re2_mean - 4.0_dp) < tol, "re2_known: Re^2 /= 4")
    if (allocated(error)) return

  end subroutine test_re2_known

  ! ------------------------------------------------------------------
  ! 4. shape_known
  !
  ! Same chain (3,0,0),(4,0,0),(5,0,0); CoM=(4,0,0); deviations:(-1,0,0),(0,0,0),(1,0,0).
  ! Gyration tensor G = (1/3)*diag(2,0,0) = diag(2/3, 0, 0).
  ! Eigenvalues (ascending): 0, 0, 2/3.
  ! sqrt(evals) sorted descending: a=sqrt(2/3), b=0, c=0.
  !
  ! Test: a_sum(kbin)/nc(kbin) ≈ sqrt(2/3) ≈ 0.8165, b≈0, c≈0 within 1e-10.
  ! ------------------------------------------------------------------
  subroutine test_shape_known(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)      :: sys
    type(polyfj_obs_t)  :: obs
    integer, parameter  :: n = 3
    real(dp), parameter :: bsz  = 0.5_dp
    real(dp), parameter :: bszc = 0.5_dp
    real(dp), parameter :: L    = 20.0_dp
    real(dp), parameter :: rmax = L / 2.0_dp
    real(dp), parameter :: tol  = 1.0e-10_dp

    real(dp) :: xs(n), ys(n), zs(n)
    real(dp) :: com_dist, a_mean, b_mean, c_mean
    integer  :: kbin

    xs(1) = 3.0_dp;  ys(1) = 0.0_dp;  zs(1) = 0.0_dp
    xs(2) = 4.0_dp;  ys(2) = 0.0_dp;  zs(2) = 0.0_dp
    xs(3) = 5.0_dp;  ys(3) = 0.0_dp;  zs(3) = 0.0_dp

    call make_system(sys, n, xs, ys, zs, L)
    call obs_init(obs, bsz, bszc, rmax)
    call obs_sample(obs, sys)

    com_dist = 4.0_dp
    kbin = 1 + int(com_dist / bszc)

    call check(error, kbin <= obs%nbinc, "shape_known: kbin out of range")
    if (allocated(error)) return

    call check(error, obs%nc(kbin) == 1_i8, "shape_known: nc(kbin) /= 1")
    if (allocated(error)) return

    a_mean = obs%a_sum(kbin) / real(obs%nc(kbin), dp)
    b_mean = obs%b_sum(kbin) / real(obs%nc(kbin), dp)
    c_mean = obs%c_sum(kbin) / real(obs%nc(kbin), dp)

    call check(error, abs(a_mean - sqrt(2.0_dp/3.0_dp)) < tol, &
               "shape_known: a /= sqrt(2/3)")
    if (allocated(error)) return

    call check(error, abs(b_mean) < tol, "shape_known: b /= 0")
    if (allocated(error)) return

    call check(error, abs(c_mean) < tol, "shape_known: c /= 0")
    if (allocated(error)) return

  end subroutine test_shape_known

  ! ------------------------------------------------------------------
  ! 5. segorder_rod
  !
  ! A chain perfectly aligned along the +x axis, far from origin.
  ! Consecutive position vectors are nearly parallel → cos ≈ 1.
  !
  ! Chain of n=5 at x = 100, 101, 102, 103, 104 (y=0, z=0).
  ! Position vectors r_j and r_{j-1} point almost exactly in +x;
  ! their dot product / (|r1||r2|) = r_{j-1}*r_j / (r_{j-1}*r_j) ≈ 1
  ! for large offsets.
  !
  ! Exact cos for two consecutive beads at (d, 0, 0) and (d+1, 0, 0):
  !   cos = d*(d+1) / (d*(d+1)) = 1.0 exactly.
  ! Wait: cos = dot(r1,r2)/(|r1||r2|) = (d*(d+1)+0+0)/(d*(d+1)) = 1 exactly.
  !
  ! Actually for a rod exactly along x-axis:
  !   r_j = (x_j, 0, 0);  dot(r_{j-1}, r_j) = x_{j-1}*x_j
  !   |r_{j-1}| = x_{j-1},  |r_j| = x_j
  !   cos = (x_{j-1}*x_j)/(x_{j-1}*x_j) = 1.0 exactly.
  !
  ! XM = |r_j - r_{j-1}|/2 = 0.5; bsz=0.5 → bin = 1+int(0.5/0.5)=2.
  ! So all pairs land in bin 2.
  ! seg_cos_sum(2)/seg_cnt(2) should = 1.0 within 1e-14.
  ! ------------------------------------------------------------------
  subroutine test_segorder_rod(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)      :: sys
    type(polyfj_obs_t)  :: obs
    integer, parameter  :: n = 5
    real(dp), parameter :: bsz  = 0.5_dp
    real(dp), parameter :: bszc = 0.5_dp
    real(dp), parameter :: L    = 300.0_dp   ! large box so min_image of CoM≈coords
    real(dp), parameter :: rmax = L / 2.0_dp
    real(dp), parameter :: tol  = 1.0e-14_dp

    real(dp) :: xs(n), ys(n), zs(n)
    real(dp) :: cos_mean
    integer  :: j, kbin_seg

    do j = 1, n
      xs(j) = 100.0_dp + real(j-1, dp)  ! 100, 101, 102, 103, 104
      ys(j) = 0.0_dp
      zs(j) = 0.0_dp
    end do

    call make_system(sys, n, xs, ys, zs, L)
    call obs_init(obs, bsz, bszc, rmax)
    call obs_sample(obs, sys)

    ! XM = 0.5 for all consecutive pairs → bin = 1+int(0.5/0.5) = 2
    kbin_seg = 2

    call check(error, obs%seg_cnt(kbin_seg) == int(n-1, i8), &
               "segorder_rod: seg_cnt(2) /= n-1")
    if (allocated(error)) return

    cos_mean = obs%seg_cos_sum(kbin_seg) / real(obs%seg_cnt(kbin_seg), dp)
    call check(error, abs(cos_mean - 1.0_dp) < tol, &
               "segorder_rod: <cos> /= 1")
    if (allocated(error)) return

  end subroutine test_segorder_rod

end module test_polyfj_observables_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_polyfj_observables
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive,                          only: run_testsuite, new_testsuite, testsuite_type
  use test_polyfj_observables_suite,      only: collect_polyfj_obs
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element array constructor is safe from gfortran-15 finalizer bug.
  testsuites = [new_testsuite("polyfj_observables", collect_polyfj_obs)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_polyfj_observables
