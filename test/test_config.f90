! test_config.f90 — unit tests for polymc_config
!
! Tests: read_runt_inp, read_polyfj_inp, read_ic, write_fc (round-trip)
!
! gfortran-15 workaround: allocate(testsuite(N)) + element-wise assignment.
! No array constructors for testsuite(:).

module test_config_suite
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use polymc_kinds, only: dp, i8
  use polymc_config, only: system_t, params_t, &
                            read_runt_inp, read_polyfj_inp, read_ic, write_fc
  implicit none
  private
  public :: collect_config

contains

  !> Register all config tests
  subroutine collect_config(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(7))
    testsuite(1) = new_unittest("runt_ic_read",      test_runt_ic_read)
    testsuite(2) = new_unittest("runt_inp_read",     test_runt_inp_read)
    testsuite(3) = new_unittest("polyfj_inp_read",   test_polyfj_inp_read)
    testsuite(4) = new_unittest("polyfj_ic_read",    test_polyfj_ic_read)
    testsuite(5) = new_unittest("rcut_present",      test_rcut_present)
    testsuite(6) = new_unittest("rcut_absent",       test_rcut_absent)
    testsuite(7) = new_unittest("round_trip",        test_round_trip)
  end subroutine collect_config

  ! ------------------------------------------------------------------
  ! Read runt.ic — check nmol, n, nbeads, box L, coords in [0,L),
  ! mol_of(1)==1, mol_of(9)==2.
  ! ------------------------------------------------------------------
  subroutine test_runt_ic_read(error)
    type(error_type), allocatable, intent(out) :: error
    type(system_t) :: sys
    integer :: k

    call read_ic('runt/runs/n8_s1.5_e0.2_ff0_sf0/runt.ic', sys, .false.)

    call check(error, sys%nmol == 200,     "runt.ic: nmol should be 200")
    if (allocated(error)) return
    call check(error, sys%n == 8,          "runt.ic: n should be 8")
    if (allocated(error)) return
    call check(error, sys%nbeads == 1600,  "runt.ic: nbeads should be 1600")
    if (allocated(error)) return
    call check(error, abs(sys%box%L - 16.12_dp) < 0.05_dp, &
               "runt.ic: box L should be ~16.12")
    if (allocated(error)) return

    ! All coords in [0, L)
    do k = 1, sys%nbeads
      call check(error, sys%x(k) >= 0.0_dp .and. sys%x(k) < sys%box%L, &
                 "runt.ic: x out of [0,L)")
      if (allocated(error)) return
      call check(error, sys%y(k) >= 0.0_dp .and. sys%y(k) < sys%box%L, &
                 "runt.ic: y out of [0,L)")
      if (allocated(error)) return
      call check(error, sys%z(k) >= 0.0_dp .and. sys%z(k) < sys%box%L, &
                 "runt.ic: z out of [0,L)")
      if (allocated(error)) return
    end do

    call check(error, sys%mol_of(1) == 1, "runt.ic: mol_of(1) should be 1")
    if (allocated(error)) return
    call check(error, sys%mol_of(9) == 2, "runt.ic: mol_of(9) should be 2")
    if (allocated(error)) return
  end subroutine test_runt_ic_read

  ! ------------------------------------------------------------------
  ! Read runt.inp — check fdick, rsp, and rcut==−1 (no 11th line).
  ! ------------------------------------------------------------------
  subroutine test_runt_inp_read(error)
    type(error_type), allocatable, intent(out) :: error
    type(params_t) :: p

    call read_runt_inp('runt/runs/n8_s1.5_e0.2_ff0_sf0/runt.inp', p)

    call check(error, abs(p%fdick - 1.0_dp/3.0_dp) < 1.0e-4_dp, &
               "runt.inp: fdick should be ~1/3")
    if (allocated(error)) return
    call check(error, abs(p%rsp - 1.5_dp) < 1.0e-6_dp, &
               "runt.inp: rsp should be 1.5")
    if (allocated(error)) return
    ! No 11th line in this fixture → rcut must remain the default -1.0
    call check(error, p%rcut == -1.0_dp, &
               "runt.inp: rcut should be -1.0 (no 11th line)")
    if (allocated(error)) return
  end subroutine test_runt_inp_read

  ! ------------------------------------------------------------------
  ! Read polyfj.inp — check ncon>0, frept~1/3.
  ! ------------------------------------------------------------------
  subroutine test_polyfj_inp_read(error)
    type(error_type), allocatable, intent(out) :: error
    type(params_t) :: p

    call read_polyfj_inp('polyfj/runs/n8_s1.5_e0.2/polyfj.inp', p)

    call check(error, p%ncon > 0_i8, &
               "polyfj.inp: ncon should be > 0")
    if (allocated(error)) return
    call check(error, abs(p%frept - 1.0_dp/3.0_dp) < 1.0e-4_dp, &
               "polyfj.inp: frept should be ~1/3")
    if (allocated(error)) return
  end subroutine test_polyfj_inp_read

  ! ------------------------------------------------------------------
  ! Read polyfj.ic — check rsp, nmol, n.
  ! ------------------------------------------------------------------
  subroutine test_polyfj_ic_read(error)
    type(error_type), allocatable, intent(out) :: error
    type(system_t) :: sys

    call read_ic('polyfj/runs/n8_s1.5_e0.2/polyfj.ic', sys, .true.)

    call check(error, abs(sys%rsp - 1.5_dp) < 1.0e-6_dp, &
               "polyfj.ic: rsp should be 1.5")
    if (allocated(error)) return
    call check(error, sys%nmol == 200, "polyfj.ic: nmol should be 200")
    if (allocated(error)) return
    call check(error, sys%n == 8,      "polyfj.ic: n should be 8")
    if (allocated(error)) return
  end subroutine test_polyfj_ic_read

  ! ------------------------------------------------------------------
  ! rcut fixture: write a synthetic runt.inp with an 11th line = 5.0,
  ! read it back, assert rcut == 5.0.
  ! ------------------------------------------------------------------
  subroutine test_rcut_present(error)
    type(error_type), allocatable, intent(out) :: error
    type(params_t) :: p
    integer :: iu, ios
    character(len=256) :: tmppath

    tmppath = 'build/test_rcut.inp'

    ! Write the temp file (10 standard lines + rcut line)
    open(newunit=iu, file=trim(tmppath), status='replace', action='write', &
         iostat=ios)
    call check(error, ios == 0, "rcut_present: could not open temp file for writing")
    if (allocated(error)) return

    write(iu, '(g0)') 2000000    ! ncon
    write(iu, '(g0)') 1000       ! nskip
    write(iu, '(g0)') 0.05       ! bsz
    write(iu, '(g0)') 0.10       ! dlr
    write(iu, '(g0)') 0.10       ! dint
    write(iu, '(g0)') 0.333333   ! fdick
    write(iu, '(g0)') 0.333333   ! frept
    write(iu, '(g0)') 0.0        ! beps
    write(iu, '(g0)') 0.0        ! bbeps
    write(iu, '(g0)') 1.5        ! rsp
    write(iu, '(g0)') 5.0        ! rcut — 11th line present
    close(iu)

    call read_runt_inp(trim(tmppath), p)

    call check(error, abs(p%rcut - 5.0_dp) < 1.0e-9_dp, &
               "rcut_present: rcut should be 5.0")
    if (allocated(error)) return
  end subroutine test_rcut_present

  ! ------------------------------------------------------------------
  ! rcut absent: read a file with only 10 lines; rcut stays -1.0.
  ! Confirmed by the committed runt.inp fixture (which has a trailing
  ! blank line 11 — list-directed read of a blank line gives EOF/EOR,
  ! not a value, so rcut should remain -1.0).
  ! ------------------------------------------------------------------
  subroutine test_rcut_absent(error)
    type(error_type), allocatable, intent(out) :: error
    type(params_t) :: p

    call read_runt_inp('runt/runs/n8_s1.5_e0.2_ff0_sf0/runt.inp', p)

    call check(error, p%rcut == -1.0_dp, &
               "rcut_absent: rcut should remain -1.0 when 11th line is absent/blank")
    if (allocated(error)) return
  end subroutine test_rcut_absent

  ! ------------------------------------------------------------------
  ! Round-trip: read runt.ic, write_fc to a temp path, read back,
  ! assert nmol/n/L/first-bead positions match.
  ! ------------------------------------------------------------------
  subroutine test_round_trip(error)
    type(error_type), allocatable, intent(out) :: error
    type(system_t) :: sys_orig, sys_back
    real(dp) :: pfc
    integer :: k
    character(len=256) :: tmppath

    tmppath = 'build/test_roundtrip.ic'
    pfc = 0.2_dp

    call read_ic('runt/runs/n8_s1.5_e0.2_ff0_sf0/runt.ic', sys_orig, .false.)
    call write_fc(trim(tmppath), sys_orig, pfc)
    call read_ic(trim(tmppath), sys_back, .false.)

    call check(error, sys_back%nmol == sys_orig%nmol, "round-trip: nmol mismatch")
    if (allocated(error)) return
    call check(error, sys_back%n == sys_orig%n,       "round-trip: n mismatch")
    if (allocated(error)) return
    call check(error, abs(sys_back%box%L - sys_orig%box%L) < 1.0e-9_dp, &
               "round-trip: box L mismatch")
    if (allocated(error)) return

    ! Check all bead positions
    do k = 1, sys_orig%nbeads
      call check(error, abs(sys_back%x(k) - sys_orig%x(k)) < 1.0e-6_dp, &
                 "round-trip: x mismatch")
      if (allocated(error)) return
      call check(error, abs(sys_back%y(k) - sys_orig%y(k)) < 1.0e-6_dp, &
                 "round-trip: y mismatch")
      if (allocated(error)) return
      call check(error, abs(sys_back%z(k) - sys_orig%z(k)) < 1.0e-6_dp, &
                 "round-trip: z mismatch")
      if (allocated(error)) return
    end do
  end subroutine test_round_trip

end module test_config_suite

! ==================================================================
!  Driver program
! ==================================================================
program test_config
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_config_suite, only: collect_config
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  testsuites = [new_testsuite("config", collect_config)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_config
