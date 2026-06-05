! test_runt_energy.f90 — unit tests for runt_energy (bead_inter_energy, total_energy)
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran-15 allocatable-finalizer bug workaround:
!   allocate(testsuite(N)) + element-wise assignment; no array ctor.
!
! System: nmol=10 chains, n=4 beads each, nbeads=40, L=8.0, rsp=1.5
! Parameters: beps=-0.5, bbeps=-0.5, rcut=-1 (exact) or 100 (cutoff)
!
! Tests:
!   a) exact_inter_matches_oracle  — bead_inter_energy (exact) == direct double loop
!   b) cutoff_huge_matches_exact   — bead_inter_energy with rcut=100 == exact value
!   c) total_energy_consistency    — 0.5*sum_k inter_k + sphere + intra == total_energy
!   d) zero_coefficients_zero_energy — beps=bbeps=0 => total_energy == 0

module test_runt_energy_suite
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,     only: dp, i8
  use polymc_rng,       only: rng_t, rng_seed, rng_uniform
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t, cl_init, cl_build, cl_neighbours
  use polymc_config,    only: system_t, params_t
  use runt_potential,   only: u_pair, u_sphere
  use runt_energy,      only: bead_inter_energy, total_energy
  implicit none
  private
  public :: collect_runt_energy

  ! Yukawa inverse-range parameter (must match runt_potential.f90)
  real(dp), parameter :: kappa = 2.5_dp

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_runt_energy(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(4))
    testsuite(1) = new_unittest("exact_inter_matches_oracle",  test_exact_inter_matches_oracle)
    testsuite(2) = new_unittest("cutoff_huge_matches_exact",   test_cutoff_huge_matches_exact)
    testsuite(3) = new_unittest("total_energy_consistency",    test_total_energy_consistency)
    testsuite(4) = new_unittest("zero_coefficients_zero_energy", test_zero_coefficients_zero_energy)
  end subroutine collect_runt_energy

  ! ------------------------------------------------------------------
  ! Build a small deterministic test system in-place.
  ! nmol=10, n=4, nbeads=40, L=8.0, rsp=1.5
  ! Positions: xoshiro256** with given seed, in [0, L)
  ! mol_of(k) = (k-1)/n + 1
  ! ------------------------------------------------------------------
  subroutine build_test_system(sys, rng_seed_val)
    type(system_t), intent(out) :: sys
    integer(i8),    intent(in)  :: rng_seed_val

    type(rng_t) :: rng
    integer :: k

    sys%nmol   = 10
    sys%n      = 4
    sys%nbeads = sys%nmol * sys%n    ! = 40
    sys%box%L  = 8.0_dp
    sys%rsp    = 1.5_dp

    allocate(sys%x(sys%nbeads), sys%y(sys%nbeads), sys%z(sys%nbeads))
    allocate(sys%mol_of(sys%nbeads))

    call rng_seed(rng, rng_seed_val)
    do k = 1, sys%nbeads
      sys%x(k)      = rng_uniform(rng) * sys%box%L
      sys%y(k)      = rng_uniform(rng) * sys%box%L
      sys%z(k)      = rng_uniform(rng) * sys%box%L
      sys%mol_of(k) = (k - 1) / sys%n + 1
    end do
  end subroutine build_test_system

  ! ------------------------------------------------------------------
  ! Oracle for inter-molecular Yukawa energy of bead at (px,py,pz)
  ! belonging to molecule imol.
  ! Computes SUM over all j with mol_of(j)/=imol of u_pair(r, beps, rcut_exact)
  ! using an INDEPENDENT direct loop — not using bead_inter_energy at all.
  ! ------------------------------------------------------------------
  real(dp) function oracle_inter(sys, beps, imol, px, py, pz)
    type(system_t), intent(in) :: sys
    real(dp),       intent(in) :: beps, px, py, pz
    integer,        intent(in) :: imol

    integer  :: j
    real(dp) :: dx, dy, dz, r

    oracle_inter = 0.0_dp
    do j = 1, sys%nbeads
      if (sys%mol_of(j) /= imol) then
        dx = min_image(sys%x(j) - px, sys%box%L)
        dy = min_image(sys%y(j) - py, sys%box%L)
        dz = min_image(sys%z(j) - pz, sys%box%L)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        ! Direct Yukawa formula (no cutoff) — written by hand, independent of u_pair
        oracle_inter = oracle_inter + beps * exp(-kappa * (r - 1.0_dp)) / r
      end if
    end do
  end function oracle_inter

  ! ------------------------------------------------------------------
  ! Test (a): exact inter == direct double loop (oracle)
  ! Test a few trial beads at existing positions with their own imol.
  ! ------------------------------------------------------------------
  subroutine test_exact_inter_matches_oracle(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t) :: sys
    type(params_t) :: p
    type(cell_list_t) :: cl   ! unused in exact mode, but required by interface
    real(dp) :: e_fn, e_oracle
    integer  :: k, test_beads(5)
    real(dp), parameter :: tol = 1.0e-10_dp

    call build_test_system(sys, 42_i8)

    p%beps  = -0.5_dp
    p%bbeps = -0.5_dp
    p%rcut  = -1.0_dp   ! exact mode

    ! Build a dummy cell list (unused in exact mode but interface requires it)
    call cl_init(cl, sys%box, max(1.0_dp, -p%rcut))  ! r_cut dummy, just init
    call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

    ! Test 5 beads at different positions/molecules
    test_beads = [1, 7, 15, 22, 38]

    do k = 1, 5
      e_fn = bead_inter_energy(sys, p, cl, &
               sys%mol_of(test_beads(k)), &
               sys%x(test_beads(k)), sys%y(test_beads(k)), sys%z(test_beads(k)))
      e_oracle = oracle_inter(sys, p%beps, &
                   sys%mol_of(test_beads(k)), &
                   sys%x(test_beads(k)), sys%y(test_beads(k)), sys%z(test_beads(k)))
      call check(error, abs(e_fn - e_oracle) < tol, &
                 "exact inter: bead_inter_energy does not match oracle")
      if (allocated(error)) return
    end do

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_exact_inter_matches_oracle

  ! ------------------------------------------------------------------
  ! Test (b): cutoff with huge rcut == exact value
  !
  ! rcut=100 forces mc=3 (L/rcut < 1, clamped to 3).
  ! With mc=3 the 3x3x3 stencil = ALL 27 cells = entire box, so every
  ! bead is a candidate.  After molecule-filtering, result must equal exact.
  ! ------------------------------------------------------------------
  subroutine test_cutoff_huge_matches_exact(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p_exact, p_cutoff
    type(cell_list_t) :: cl_exact, cl_big
    real(dp) :: e_exact, e_cutoff
    integer  :: k, test_beads(3)
    real(dp), parameter :: tol = 1.0e-10_dp

    call build_test_system(sys, 99_i8)

    p_exact%beps  = -0.5_dp
    p_exact%bbeps = -0.5_dp
    p_exact%rcut  = -1.0_dp   ! exact mode

    p_cutoff%beps  = -0.5_dp
    p_cutoff%bbeps = -0.5_dp
    p_cutoff%rcut  = 100.0_dp  ! huge cutoff

    ! Build cell lists
    ! Exact mode: use a modest r_cut for the cell list (unused in exact path)
    call cl_init(cl_exact, sys%box, 1.0_dp)
    call cl_build(cl_exact, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

    ! Cutoff mode: r_cut=100 => mc=max(3, int(8/100))=3; 3x3x3=all cells
    call cl_init(cl_big, sys%box, p_cutoff%rcut)
    call cl_build(cl_big, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

    ! Self-check: with mc=3 the 3x3x3 stencil covers the whole box, so
    ! cl_neighbours must return ALL nbeads as candidates for any query point.
    ! This guards against future cell-list changes silently narrowing the stencil.
    block
      integer :: idx_chk(sys%nbeads), nidx_chk
      call cl_neighbours(cl_big, sys%box, &
                         sys%x(1), sys%y(1), sys%z(1), idx_chk, nidx_chk)
      call check(error, nidx_chk == sys%nbeads, &
                 "cutoff huge rcut: cl_neighbours did not return all beads (mc=3 stencil broken)")
      if (allocated(error)) return
    end block

    test_beads = [3, 18, 35]

    do k = 1, 3
      e_exact = bead_inter_energy(sys, p_exact, cl_exact, &
                  sys%mol_of(test_beads(k)), &
                  sys%x(test_beads(k)), sys%y(test_beads(k)), sys%z(test_beads(k)))
      e_cutoff = bead_inter_energy(sys, p_cutoff, cl_big, &
                   sys%mol_of(test_beads(k)), &
                   sys%x(test_beads(k)), sys%y(test_beads(k)), sys%z(test_beads(k)))
      call check(error, abs(e_exact - e_cutoff) < tol, &
                 "cutoff huge rcut: does not match exact value")
      if (allocated(error)) return
    end do

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_cutoff_huge_matches_exact

  ! ------------------------------------------------------------------
  ! Test (c): total_energy consistency
  !
  ! Verify: total_energy = inter + sphere + intra
  ! where inter = 0.5 * sum_k bead_inter_energy(mol_of(k), pos_k)
  ! and sphere, intra are computed here directly.
  ! Equivalently: 2*(total_energy - sphere - intra) == sum_k inter_k
  ! ------------------------------------------------------------------
  subroutine test_total_energy_consistency(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl

    real(dp) :: e_total
    real(dp) :: inter_sum   ! sum_k bead_inter_energy(k)
    real(dp) :: sphere_sum  ! sum_k u_sphere(r_k, ...)
    real(dp) :: intra_sum   ! sum over chains, non-bonded intra pairs

    integer  :: k, i, j, ij, ik
    real(dp) :: dx, dy, dz, r, r_k
    real(dp), parameter :: tol = 1.0e-9_dp

    call build_test_system(sys, 77_i8)

    p%beps  = -0.5_dp
    p%bbeps = -0.5_dp
    p%rcut  = -1.0_dp   ! exact mode

    ! Build cell list (unused in exact mode)
    call cl_init(cl, sys%box, 1.0_dp)
    call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

    ! Compute total_energy
    e_total = total_energy(sys, p, cl)

    ! Compute sphere sum: r_k = min-image distance of bead k to origin
    sphere_sum = 0.0_dp
    do k = 1, sys%nbeads
      dx = min_image(sys%x(k), sys%box%L)
      dy = min_image(sys%y(k), sys%box%L)
      dz = min_image(sys%z(k), sys%box%L)
      r_k = sqrt(dx*dx + dy*dy + dz*dz)
      sphere_sum = sphere_sum + u_sphere(r_k, p%bbeps, sys%rsp, p%rcut)
    end do

    ! Compute intra sum: for each chain i, pairs (j, k) with j=3..n, k=1..j-2
    intra_sum = 0.0_dp
    do i = 1, sys%nmol
      do j = 3, sys%n
        do k = 1, j - 2
          ij = (i - 1)*sys%n + j
          ik = (i - 1)*sys%n + k
          dx = min_image(sys%x(ij) - sys%x(ik), sys%box%L)
          dy = min_image(sys%y(ij) - sys%y(ik), sys%box%L)
          dz = min_image(sys%z(ij) - sys%z(ik), sys%box%L)
          r = sqrt(dx*dx + dy*dy + dz*dz)
          intra_sum = intra_sum + u_pair(r, p%beps, p%rcut)
        end do
      end do
    end do

    ! Compute inter sum = sum_k bead_inter_energy(k)
    inter_sum = 0.0_dp
    do k = 1, sys%nbeads
      inter_sum = inter_sum + bead_inter_energy(sys, p, cl, &
                    sys%mol_of(k), sys%x(k), sys%y(k), sys%z(k))
    end do

    ! Check: total_energy == 0.5*inter_sum + sphere_sum + intra_sum
    call check(error, abs(e_total - (0.5_dp*inter_sum + sphere_sum + intra_sum)) < tol, &
               "total_energy inconsistent: does not match 0.5*inter + sphere + intra")

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_total_energy_consistency

  ! ------------------------------------------------------------------
  ! Test (d): zero beps and bbeps => total_energy == 0
  ! ------------------------------------------------------------------
  subroutine test_zero_coefficients_zero_energy(error)
    type(error_type), allocatable, intent(out) :: error

    type(system_t)    :: sys
    type(params_t)    :: p
    type(cell_list_t) :: cl
    real(dp) :: e
    real(dp), parameter :: tol = 1.0e-12_dp

    call build_test_system(sys, 13_i8)

    p%beps  = 0.0_dp
    p%bbeps = 0.0_dp
    p%rcut  = -1.0_dp

    call cl_init(cl, sys%box, 1.0_dp)
    call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

    e = total_energy(sys, p, cl)
    call check(error, abs(e) < tol, &
               "zero coefficients: total_energy should be 0")

    deallocate(sys%x, sys%y, sys%z, sys%mol_of)
  end subroutine test_zero_coefficients_zero_energy

end module test_runt_energy_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_runt_energy
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_runt_energy_suite, only: collect_runt_energy
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element constructor is safe from the gfortran-15 finalizer bug.
  testsuites = [new_testsuite("runt_energy", collect_runt_energy)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_runt_energy
