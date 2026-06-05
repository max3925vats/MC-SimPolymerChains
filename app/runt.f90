! app/runt.f90 — driver for the runt (big-bead chain) MC simulation.
!
! Mirrors the legacy runt/runt.f main-loop logic in modern Fortran, using
! real (sigma) units throughout.  Output files follow legacy FORMAT statements
! closely enough that existing Python tooling / numpy.loadtxt can parse them.
!
! Build via fpm:  fpm build
! Run:            cd <deck_dir> && /path/to/build/.../runt
!
! Expected inputs (current working directory):
!   runt.inp   — run parameters (read_runt_inp)
!   runt.ic    — initial configuration (read_ic, has_rsp=.false.)
!
! Outputs:
!   runt.fc    — restart configuration (write_fc)
!   runt.out   — results summary with move statistics and density profile
!   denrunt.out — 3-column density profile: bin_idx  dist  density

program runt_driver
  use polymc_kinds,     only: dp, i8
  use polymc_rng,       only: rng_t, rng_seed, rng_uniform
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t, cl_init, cl_build
  use polymc_config,    only: system_t, params_t, read_runt_inp, read_ic, write_fc
  use polymc_binning,   only: hist_t, hist_init, hist_add, hist_density
  use runt_energy,      only: total_energy
  use runt_moves,       only: move_dick, move_rept, move_ccb
  implicit none

  real(dp), parameter :: pi = 3.141592653589793_dp

  ! --- Simulation objects ---
  type(system_t)    :: sys
  type(params_t)    :: p
  type(cell_list_t) :: cl
  type(rng_t)       :: rng
  type(hist_t)      :: hden

  ! --- Counters ---
  integer(i8) :: icon
  integer :: naver
  integer :: ntd, ntr, ntj   ! move attempts: Dickman, Reptation, CCB
  integer :: nsd, nsr, nsj   ! move successes
  real(dp) :: fsd, fsr, fsj  ! success fractions

  ! --- Energy / accumulation ---
  real(dp) :: energy, de, avener

  ! --- Loop / move temporaries ---
  integer  :: imol
  real(dp) :: xran
  logical  :: accepted

  ! --- Histogram / density ---
  real(dp) :: rmax
  integer  :: nbin, ibin
  real(dp), allocatable :: dist_arr(:), dens_arr(:)

  ! --- Output / geometry ---
  real(dp) :: pfc, vol, dist_k
  integer  :: k, iu_out, iu_den

  ! ==================================================================
  ! 1. Read inputs
  ! ==================================================================
  call read_runt_inp('runt.inp', p)
  call read_ic('runt.ic', sys, .false.)

  ! runt's sphere radius comes from the .inp, not the .ic header
  sys%rsp = p%rsp

  ! Guard: nskip must be positive — mod(icon, nskip) in the accumulation
  ! step would divide by zero for nskip <= 0, producing wrong results or
  ! a runtime fault on some compilers.
  if (p%nskip <= 0) error stop "runt: nskip must be > 0"

  ! ==================================================================
  ! 2. Seed RNG (deterministic; seed default 12345)
  !    Runs are reproducible by design: the same seed always produces
  !    the same trajectory.  To vary the run, change the 'seed' value
  !    in runt.inp (or expose it via command-line / environment).
  !    (The legacy code seeded from the wall clock, making exact
  !    reproduction impossible — this rewrite deliberately avoids that.)
  ! ==================================================================
  call rng_seed(rng, p%seed)

  ! ==================================================================
  ! 3. Report energy mode
  ! ==================================================================
  if (p%rcut <= 0.0_dp) then
    write(*, '(A)') 'energy mode: EXACT (no cutoff)'
  else
    write(*, '(A,F8.4)') 'energy mode: CUTOFF rcut=', p%rcut
  end if

  ! ==================================================================
  ! 4. Build cell list
  ! ==================================================================
  call cl_init(cl, sys%box, max(1.0_dp, p%rcut))
  call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)

  ! ==================================================================
  ! 5. Initial energy
  ! ==================================================================
  energy = total_energy(sys, p, cl)
  write(*, '(A,ES20.10)') 'INITIAL ENERGY', energy

  ! ==================================================================
  ! 6. Set up radial histogram
  !    rmax = half the box body-diagonal = L*sqrt(3)/2
  ! ==================================================================
  rmax = sys%box%L * sqrt(3.0_dp) / 2.0_dp
  call hist_init(hden, p%bsz, rmax)
  nbin = hden%nbin

  ! ==================================================================
  ! 7. Initialise counters
  ! ==================================================================
  ntd   = 0;  ntr  = 0;  ntj  = 0
  nsd   = 0;  nsr  = 0;  nsj  = 0
  naver = 0
  avener = 0.0_dp

  ! ==================================================================
  ! 8. Main MC loop
  ! ==================================================================
  write(*, '(A)') 'Simulation begun'

  do icon = 1_i8, p%ncon

    ! Pick a molecule at random (1-based)
    imol = int(sys%nmol * rng_uniform(rng)) + 1

    ! Choose move type by fraction: Dickman / Reptation / CCB
    xran = rng_uniform(rng)

    if (xran < p%fdick) then
      ntd = ntd + 1
      call move_dick(sys, p, cl, rng, imol, accepted, de)
      if (accepted) nsd = nsd + 1

    else if (xran < p%fdick + p%frept) then
      ntr = ntr + 1
      call move_rept(sys, p, cl, rng, imol, accepted, de)
      if (accepted) nsr = nsr + 1

    else
      ntj = ntj + 1
      call move_ccb(sys, p, cl, rng, imol, accepted, de)
      if (accepted) nsj = nsj + 1
    end if

    ! Update energy and rebuild cell list after accepted move
    if (accepted) then
      energy = energy + de
      call cl_build(cl, sys%box, sys%x, sys%y, sys%z, sys%nbeads)
    end if

    ! Accumulation step
    if (mod(icon, int(p%nskip, i8)) == 0_i8) then
      naver  = naver + 1
      avener = avener + energy

      ! Accumulate radial histogram: min-image distance of each bead to origin
      do k = 1, sys%nbeads
        dist_k = sqrt(  min_image(sys%x(k), sys%box%L)**2 &
                      + min_image(sys%y(k), sys%box%L)**2 &
                      + min_image(sys%z(k), sys%box%L)**2 )
        call hist_add(hden, dist_k)
      end do
    end if

  end do

  write(*, '(A)') 'Ending simulation'

  ! ==================================================================
  ! 9. Compute summary statistics
  ! ==================================================================
  if (ntd > 0) then; fsd = real(nsd, dp) / real(ntd, dp)
  else;              fsd = 0.0_dp
  end if
  if (ntr > 0) then; fsr = real(nsr, dp) / real(ntr, dp)
  else;              fsr = 0.0_dp
  end if
  if (ntj > 0) then; fsj = real(nsj, dp) / real(ntj, dp)
  else;              fsj = 0.0_dp
  end if

  ! Average energy per bead (guard against naver=0 or nbeads=0)
  if (naver > 0 .and. sys%nbeads > 0) then
    avener = avener / real(naver, dp) / real(sys%nbeads, dp)
  else
    avener = 0.0_dp
  end if
  write(*, '(A,ES20.10)') 'AVG ENERGY PER BEAD', avener

  ! ==================================================================
  ! 10a. Write restart file runt.fc
  !      pfc = (pi/6) * nmol * n / L^3  (packing fraction)
  ! ==================================================================
  vol = sys%box%L**3
  pfc = (pi / 6.0_dp) * real(sys%nmol, dp) * real(sys%n, dp) / vol
  call write_fc('runt.fc', sys, pfc)

  ! ==================================================================
  ! 10b. Compute density profile arrays
  ! ==================================================================
  allocate(dist_arr(nbin), dens_arr(nbin))
  call hist_density(hden, naver, dist_arr, dens_arr)

  ! ==================================================================
  ! 10c. Write runt.out — results summary
  ! ==================================================================
  open(newunit=iu_out, file='runt.out', status='replace', action='write')

  write(iu_out, '(A)') 'RESULTS OF HARD-YUKAWA CHAIN BIG BEAD RUN'
  write(iu_out, *)

  ! Packing fraction, chain geometry, MC parameters
  write(iu_out, 113) pfc, sys%n, sys%nmol, p%dlr, p%dint
  write(iu_out, 114) p%ncon, p%nskip
  write(iu_out, 115) p%bsz, nbin
  write(iu_out, 116) p%beps, p%bbeps

  ! Energy mode line (extra info vs legacy)
  if (p%rcut <= 0.0_dp) then
    write(iu_out, '(1X,A)') 'energy mode: EXACT (no cutoff)'
  else
    write(iu_out, '(1X,A,F8.4)') 'energy mode: CUTOFF rcut=', p%rcut
  end if

  ! Move statistics
  write(iu_out, '(1X,A)') 'Move statistics'
  write(iu_out, 117) ntd, ntr, ntj, nsd, nsr, nsj, fsd, fsr, fsj

  ! Box geometry
  write(iu_out, 118) sys%box%L, sys%box%L

  ! Density profile header + per-bin lines
  write(iu_out, '(1X,A)') 'DENSITY PROFILES NEAR THE BIG BEAD'
  do ibin = 1, nbin
    if (naver > 0) then
      write(iu_out, 1140) ibin, dist_arr(ibin), dens_arr(ibin)
    end if
  end do

  close(iu_out)

  ! ==================================================================
  ! 10d. Write denrunt.out — 3-column file for Python loadtxt
  !      Columns: bin_idx  dist  density   (FORMAT 1140)
  ! ==================================================================
  open(newunit=iu_den, file='denrunt.out', status='replace', action='write')
  do ibin = 1, nbin
    if (naver > 0) then
      write(iu_den, 1140) ibin, dist_arr(ibin), dens_arr(ibin)
    end if
  end do
  close(iu_den)

  write(*, '(A)') 'files created'

  ! --- FORMAT statements (match legacy runt.f) ---
113 format(1X,'PACKING FRACTION ',4X,F6.4 / &
           1X,' NUMBER OF BEADS ',6X,I4   / &
           1X,'NUMBER OF CHAINS ',6X,I4   / &
           1X,'     DELTA CHAIN ',5X,F5.3 / &
           1X,'  DELTA INTERNAL ',5X,F5.3 )
114 format(1X / 1X,'NO. TOTAL MOVES',2X,I10 / &
                1X,'NUMBER OF SKIPS',3X,I7 /)
115 format(1X,'      BIN SIZE',3X,F5.3 / &
           1X,'NUMBER OF BINS',4X,I4 /)
116 format(1X,'      BEPS',3X,F10.4 / &
           1X,'     BBEPS',3X,F10.4 /)
117 format(1X / 13X,'DICKMAN  REPTATION       CCB' / &
           1X,'ATTEMPTED',3X,I7,4X,I7,3X,I7 / &
           1X,'SUCCESSFUL',2X,I7,4X,I7,3X,I7 / &
           1X,'FRACTION',5X,F6.4,5X,F6.4,4X,F6.4 /)
118 format(1X,'BOX LENGTH',3X,E16.8 / &
           1X,'PERIODIC LENGTH',4X,E16.8 /)
! 1140: matches legacy FORMAT (1X,I4,1X,F6.3,3X,...) — 3 parseable numeric columns
1140 format(1X,I4,1X,F8.4,3X,ES14.6)

end program runt_driver
