! app/polyfj.f90 — driver for the polyfj (hard-chain near a sphere) MC simulation.
!
! Mirrors the legacy polyfj/polyfj.f main-loop structure in modern Fortran,
! using real (sigma) units throughout.  Output files are numeric and parseable
! by numpy.loadtxt; the observable style (corrected vs legacy) is switchable
! at runtime via the --legacy-obs command-line flag.
!
! Build via fpm:  fpm build
! Run:            cd <deck_dir> && /path/to/build/.../polyfj
!                 cd <deck_dir> && /path/to/build/.../polyfj --legacy-obs
!
! Expected inputs (current working directory):
!   polyfj.inp   — run parameters (read_polyfj_inp)
!   polyfj.ic    — initial configuration (read_ic, has_rsp=.true.)
!
! Outputs:
!   polyden.out  — 5 cols: idx dist total_density end_density mid_density
!   seg.out      — 3 cols: idx dist S  (segmental order parameter)
!   rgfj.out     — 5 cols: idx dist rho_com Rg2 Re2
!   g2fj.out     — 5 cols: idx dist a b c  (chain shape semi-axes)
!   polyfj.out   — human-readable summary with move statistics
!   polyfj.fc    — restart configuration (write_fc)
!
! Observable style (--legacy-obs flag):
!   corrected (default):
!     seg:    S = (3*mean_cos - 1)/2
!     g2fj:   a/b/c = mean of per-chain sqrt(gyration eigenvalue)   (a_sum/nc)
!   legacy:
!     seg:    S = (3*mean_cos*L^2 - 1)/2  (faithful to AL^2 artifact in polyfj.f)
!     g2fj:   a/b/c = sqrt(acig_sum/nc)   (sqrt-of-mean, legacy polyfj.f:496-498)

program polyfj_driver
  use polymc_kinds,        only: dp, i8
  use polymc_rng,          only: rng_t, rng_seed, rng_uniform
  use polymc_box,          only: box_t
  use polymc_cell_list,    only: cell_list_t, cl_init, cl_build
  use polymc_config,       only: system_t, params_t, read_polyfj_inp, read_ic, write_fc
  use polymc_binning,      only: hist_t, hist_density
  use polyfj_moves,        only: move_dick, move_rept, move_ccb
  use polyfj_observables,  only: polyfj_obs_t, obs_init, obs_sample
  implicit none

  real(dp), parameter :: pi = 3.141592653589793_dp
  real(dp), parameter :: bszcm_default = 0.10_dp  ! legacy BSZCM

  ! --- Simulation objects ---
  type(system_t)       :: sys
  type(params_t)       :: p
  type(cell_list_t)    :: cl
  type(rng_t)          :: rng
  type(polyfj_obs_t)   :: obs

  ! --- Observable style switch ---
  logical :: legacy_obs

  ! --- Counters ---
  integer(i8) :: icon
  integer :: naver
  integer :: ntd, ntr, ntj   ! move attempts: Dickman, Reptation, CCB
  integer :: nsd, nsr, nsj   ! move successes
  real(dp) :: fsd, fsr, fsj  ! success fractions

  ! --- Loop / move temporaries ---
  integer  :: imol, k
  real(dp) :: xran
  logical  :: accepted

  ! --- Geometry / binning ---
  real(dp) :: rmax, bszc
  real(dp) :: pfc, vol_sphere

  ! --- Density output arrays ---
  integer :: nbin_den
  real(dp), allocatable :: dist_arr(:), dens_all(:), dens_end(:), dens_mid(:)

  ! --- Output unit variables ---
  integer :: iu_out, iu_den, iu_seg, iu_rg, iu_g2

  ! --- Temporaries for output ---
  real(dp) :: rl, ru, vol_shell, dist_k
  real(dp) :: rho_com, rg2_mean, re2_mean
  real(dp) :: a_out, b_out, c_out
  real(dp) :: mean_cos, s_param

  ! ==================================================================
  ! 0. Parse command-line arguments for observable style flag
  ! ==================================================================
  legacy_obs = .false.
  block
    integer :: iarg, nargs
    character(len=256) :: arg
    nargs = command_argument_count()
    do iarg = 1, nargs
      call get_command_argument(iarg, arg)
      if (trim(arg) == '--legacy-obs') then
        legacy_obs = .true.
        exit
      end if
    end do
  end block

  if (legacy_obs) then
    write(*, '(A)') 'observable style: LEGACY  (seg: S=(3*cos*L^2-1)/2; g2fj: sqrt-of-mean)'
  else
    write(*, '(A)') 'observable style: CORRECTED  (seg: S=(3*cos-1)/2; g2fj: mean-then-output)'
  end if

  ! ==================================================================
  ! 1. Read inputs
  ! ==================================================================
  call read_polyfj_inp('polyfj.inp', p)
  ! polyfj.ic carries RSP — set by read_ic (has_rsp=.true.)
  call read_ic('polyfj.ic', sys, .true.)

  ! Guard: nskip must be positive
  if (p%nskip <= 0) error stop "polyfj: nskip must be > 0"

  ! ==================================================================
  ! 2. Seed RNG (deterministic; default seed 12345 from params_t)
  !    Runs are reproducible by design: same seed → same trajectory.
  !    (Legacy seeded from the wall clock, making reproduction impossible;
  !    this rewrite deliberately avoids that.)
  ! ==================================================================
  call rng_seed(rng, p%seed)

  ! ==================================================================
  ! 3. Build cell list (athermal: overlap detection only; cutoff = 1.0 sigma)
  ! ==================================================================
  call cl_init(cl, sys%box, 1.0_dp)
  call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

  ! ==================================================================
  ! 4. Initialise observables
  !    rmax = half box body-diagonal, consistent with legacy INT(AL/2/BSZ)
  !    bszc = legacy BSZCM = 0.10 sigma
  ! ==================================================================
  rmax = sys%box%L * sqrt(3.0_dp) / 2.0_dp
  bszc = bszcm_default
  call obs_init(obs, p%bsz, bszc, rmax)

  ! ==================================================================
  ! 5. Initialise counters
  ! ==================================================================
  ntd   = 0;  ntr  = 0;  ntj  = 0
  nsd   = 0;  nsr  = 0;  nsj  = 0
  naver = 0

  ! ==================================================================
  ! 6. Main MC loop
  !    Move-selection order mirrors legacy polyfj.f lines 193-211:
  !      REPT first (xran < frept), then DICK (xran < fdick+frept), else CCB.
  ! ==================================================================
  write(*, '(A)') 'Simulation begun'

  do icon = 1_i8, p%ncon

    ! Pick a molecule at random (1-based)
    imol = int(sys%nmol * rng_uniform(rng)) + 1

    ! Choose move type (legacy order: REPT, then DICK, then CCB)
    xran = rng_uniform(rng)

    if (xran < p%frept) then
      ! REPTATION
      ntr = ntr + 1
      call move_rept(sys, p, cl, rng, imol, accepted)
      if (accepted) nsr = nsr + 1

    else if (xran < p%fdick + p%frept) then
      ! DICKMAN regrowth
      ntd = ntd + 1
      call move_dick(sys, p, cl, rng, imol, accepted)
      if (accepted) nsd = nsd + 1

    else
      ! CCB (configurational-bias)
      ntj = ntj + 1
      call move_ccb(sys, p, cl, rng, imol, accepted)
      if (accepted) nsj = nsj + 1
    end if

    ! Rebuild cell list after each accepted move
    if (accepted) then
      call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)
    end if

    ! Accumulation step
    if (mod(icon, int(p%nskip, i8)) == 0_i8) then
      naver = naver + 1
      call obs_sample(obs, sys)
    end if

  end do

  write(*, '(A)') 'Ending simulation'

  ! ==================================================================
  ! 7. Compute move-fraction statistics
  ! ==================================================================
  fsd = 0.0_dp;  fsr = 0.0_dp;  fsj = 0.0_dp
  if (ntd > 0) fsd = real(nsd, dp) / real(ntd, dp)
  if (ntr > 0) fsr = real(nsr, dp) / real(ntr, dp)
  if (ntj > 0) fsj = real(nsj, dp) / real(ntj, dp)

  ! ==================================================================
  ! 8. Compute packing fraction (using sphere-excluded volume — legacy)
  !    pfc = (pi/6) * nmol * n / L^3
  !    (Note: legacy PFCMC uses VOL = L^3 - (4/3)*pi*rsp^3 denominator
  !     for the run summary, but polyfj.fc header uses the simpler PFC
  !     read from .ic; we use the same formula for write_fc here.)
  ! ==================================================================
  pfc = (pi / 6.0_dp) * real(sys%nmol, dp) * real(sys%n, dp) / sys%box%L**3

  ! ==================================================================
  ! 9. Write polyden.out — 5 cols: idx dist total end mid
  !    Uses hist_density for all three density histograms.
  ! ==================================================================
  nbin_den = obs%h_all%nbin
  allocate(dist_arr(nbin_den), dens_all(nbin_den), dens_end(nbin_den), dens_mid(nbin_den))

  call hist_density(obs%h_all, naver, dist_arr, dens_all)
  call hist_density(obs%h_end, naver, dist_arr, dens_end)
  call hist_density(obs%h_mid, naver, dist_arr, dens_mid)

  open(newunit=iu_den, file='polyden.out', status='replace', action='write')
  do k = 1, nbin_den
    write(iu_den, '(I10,1X,F12.5,3(3X,E14.6))') &
      k, dist_arr(k), dens_all(k), dens_end(k), dens_mid(k)
  end do
  close(iu_den)

  ! ==================================================================
  ! 10. Write seg.out — 3 cols: idx dist S
  !     Observable style controls the S formula:
  !       corrected: S = (3*mean_cos - 1)/2
  !       legacy:    S = (3*mean_cos*L^2 - 1)/2  (faithful AL^2 artifact)
  ! ==================================================================
  open(newunit=iu_seg, file='seg.out', status='replace', action='write')
  do k = 1, obs%nbin
    if (obs%seg_cnt(k) > 0_i8) then
      mean_cos = obs%seg_cos_sum(k) / real(obs%seg_cnt(k), dp)
      if (legacy_obs) then
        ! Faithful to legacy polyfj.f line 285: SR(K) += COSO*AL*AL,
        ! then SRI = SR(I)/NR(I) and SRI = (3*SRI - 1)/2.
        ! mean_cos here has no AL^2 factor (obs_sample accumulates plain coso).
        ! Multiply by L^2 here to reproduce legacy behavior.
        s_param = (3.0_dp * mean_cos * sys%box%L**2 - 1.0_dp) / 2.0_dp
      else
        s_param = (3.0_dp * mean_cos - 1.0_dp) / 2.0_dp
      end if
      dist_k = (real(k, dp) - 0.5_dp) * obs%bsz
      write(iu_seg, '(I10,1X,F12.5,3X,E14.6)') k, dist_k, s_param
    end if
  end do
  close(iu_seg)

  ! ==================================================================
  ! 11. Write rgfj.out — 5 cols: idx dist rho_com Rg2 Re2
  !     rho_com = nc / (shell_vol * naver)  where shell uses bszc bins.
  !     (This is a number density of chain CoMs per unit volume.)
  ! ==================================================================
  open(newunit=iu_rg, file='rgfj.out', status='replace', action='write')
  do k = 1, obs%nbinc
    if (obs%nc(k) > 0_i8) then
      rl = real(k-1, dp) * obs%bszc
      ru = real(k,   dp) * obs%bszc
      vol_shell = (4.0_dp / 3.0_dp) * pi * (ru**3 - rl**3) * real(naver, dp)
      rho_com   = real(obs%nc(k), dp) / vol_shell
      rg2_mean  = obs%rg2_sum(k) / real(obs%nc(k), dp)
      re2_mean  = obs%re2_sum(k) / real(obs%nc(k), dp)
      dist_k    = (real(k, dp) - 0.5_dp) * obs%bszc
      write(iu_rg, '(I10,1X,F12.5,3(3X,E14.6))') k, dist_k, rho_com, rg2_mean, re2_mean
    end if
  end do
  close(iu_rg)

  ! ==================================================================
  ! 12. Write g2fj.out — 5 cols: idx dist a b c
  !     Observable style controls finalisation:
  !       corrected: a = a_sum/nc  (mean of per-chain sqrt(gyration eigenvalue))
  !       legacy:    a = sqrt(acig_sum/nc)  (sqrt-of-mean, polyfj.f lines 496-498)
  ! ==================================================================
  open(newunit=iu_g2, file='g2fj.out', status='replace', action='write')
  do k = 1, obs%nbinc
    if (obs%nc(k) > 0_i8) then
      dist_k = (real(k, dp) - 0.5_dp) * obs%bszc
      if (legacy_obs) then
        ! Legacy: ACIGI = sqrt(ACIG(K) / NC(K)), matching polyfj.f lines 496-498.
        a_out = sqrt(max(obs%acig_sum(k) / real(obs%nc(k), dp), 0.0_dp))
        b_out = sqrt(max(obs%bcig_sum(k) / real(obs%nc(k), dp), 0.0_dp))
        c_out = sqrt(max(obs%ccig_sum(k) / real(obs%nc(k), dp), 0.0_dp))
      else
        a_out = obs%a_sum(k) / real(obs%nc(k), dp)
        b_out = obs%b_sum(k) / real(obs%nc(k), dp)
        c_out = obs%c_sum(k) / real(obs%nc(k), dp)
      end if
      write(iu_g2, '(I10,1X,F12.5,3(3X,E14.6))') k, dist_k, a_out, b_out, c_out
    end if
  end do
  close(iu_g2)

  ! ==================================================================
  ! 13. Write polyfj.out — human-readable summary
  ! ==================================================================
  open(newunit=iu_out, file='polyfj.out', status='replace', action='write')

  write(iu_out, '(A)') 'RESULTS OF CHAIN NEAR A SPHERE RUN (modern-fortran-rewrite)'
  write(iu_out, *)
  write(iu_out, '(A,F10.4)') 'PACKING FRACTION        ', pfc
  write(iu_out, '(A,I6)')    'NUMBER OF BEADS         ', sys%n
  write(iu_out, '(A,I8)')    'NUMBER OF CHAINS        ', sys%nmol
  write(iu_out, '(A,F10.4)') 'RADIUS OF LARGER SPHERE ', sys%rsp
  write(iu_out, '(A,F10.4)') 'BOX LENGTH (sigma)      ', sys%box%L
  write(iu_out, *)
  write(iu_out, '(A,I0)')    'NO. TOTAL MOVES   ', p%ncon
  write(iu_out, '(A,I0)')    'NUMBER OF SKIPS   ', p%nskip
  write(iu_out, '(A,F10.4)') 'BIN SIZE (sigma)  ', p%bsz
  write(iu_out, '(A,F10.4)') 'BSZC (sigma)      ', bszc
  write(iu_out, '(A,I0)')    'NBIN (density)    ', obs%nbin
  write(iu_out, '(A,I0)')    'NBINC (CoM)       ', obs%nbinc
  write(iu_out, *)
  if (legacy_obs) then
    write(iu_out, '(A)') 'OBSERVABLE STYLE: LEGACY (seg: S=(3*cos*L^2-1)/2; g2fj: sqrt-of-mean)'
  else
    write(iu_out, '(A)') 'OBSERVABLE STYLE: CORRECTED (seg: S=(3*cos-1)/2; g2fj: mean-then-output)'
  end if
  write(iu_out, *)
  write(iu_out, '(A)') 'Move statistics:'
  write(iu_out, '(A)')          '              REPTATION      DICKMAN         CCB'
  write(iu_out, '(A,I10,I13,I13)') 'ATTEMPTED ', ntr, ntd, ntj
  write(iu_out, '(A,I10,I13,I13)') 'SUCCESSFUL', nsr, nsd, nsj
  write(iu_out, '(A,F10.5,F13.5,F13.5)') 'FRACTION  ', fsr, fsd, fsj
  write(iu_out, *)
  write(iu_out, '(A)') 'DENSITY PROFILES NEAR A SPHERE (idx  dist  total_dens  end_dens  mid_dens)'
  do k = 1, nbin_den
    if (naver > 0) then
      write(iu_out, '(I10,1X,F12.5,3(3X,E14.6))') &
        k, dist_arr(k), dens_all(k), dens_end(k), dens_mid(k)
    end if
  end do

  close(iu_out)

  ! ==================================================================
  ! 14. Write polyfj.fc — restart configuration
  ! ==================================================================
  call write_fc('polyfj.fc', sys, pfc)

  ! ==================================================================
  ! 15. Print move statistics to stdout
  ! ==================================================================
  write(*, *)
  write(*, '(A)') '--- Move statistics ---'
  write(*, '(A)')          '              REPTATION      DICKMAN         CCB'
  write(*, '(A,I10,I13,I13)') 'ATTEMPTED ', ntr, ntd, ntj
  write(*, '(A,I10,I13,I13)') 'SUCCESSFUL', nsr, nsd, nsj
  write(*, '(A,F10.5,F13.5,F13.5)') 'FRACTION  ', fsr, fsd, fsj
  write(*, *)
  write(*, '(A,I0,A)') 'naver = ', naver, ' snapshots accumulated'
  write(*, '(A)') 'outputs: polyden.out  seg.out  rgfj.out  g2fj.out  polyfj.out  polyfj.fc'

end program polyfj_driver
