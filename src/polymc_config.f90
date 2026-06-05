! polymc_config.f90 — simulation configuration types and deck I/O
!
! Public interface:
!   system_t        — particle/box state (positions in real sigma units)
!   params_t        — MC run parameters
!   read_runt_inp   — parse runt.inp (10 required values + optional rcut)
!   read_polyfj_inp — parse polyfj.inp (7 values)
!   read_ic         — read initial configuration (.ic file)
!   write_fc        — write restart configuration (.fc / .ic-compatible file)
!
! Coordinate convention:
!   Positions are stored in real (sigma) units in [0, L).
!   The legacy code divided by AL internally; this module does NOT — coords
!   are read and stored as-is.
!
! Blank-line handling:
!   Use read(unit,'(A)') dummy to consume exactly one record (blank or not).
!   Use read(unit,*) for value lines (list-directed).

module polymc_config
  use polymc_kinds, only: dp, i8
  use polymc_box,   only: box_t
  implicit none
  private
  public :: system_t, params_t, read_runt_inp, read_polyfj_inp, read_ic, write_fc

  !> Full system state: molecules, beads, box geometry, and positions.
  type :: system_t
    integer :: nmol  = 0   ! number of molecules
    integer :: n     = 0   ! beads per molecule
    integer :: nbeads = 0  ! total beads = nmol * n
    type(box_t) :: box     ! simulation box (L in sigma)
    real(dp) :: rsp = 0.0_dp  ! sphere radius (polyfj only; 0 for runt)
    real(dp), allocatable :: x(:), y(:), z(:)  ! bead positions [0, L)
    integer,  allocatable :: mol_of(:)          ! molecule index of each bead
  end type system_t

  !> MC run parameters parsed from .inp files.
  type :: params_t
    integer(i8) :: ncon  = 0          ! number of MC steps
    integer(i8) :: seed  = 12345_i8   ! RNG seed (not in .inp; kept for driver)
    integer  :: nskip    = 0          ! output frequency
    real(dp) :: bsz      = 0.0_dp    ! max bead displacement
    real(dp) :: dlr      = 0.0_dp    ! reptation step
    real(dp) :: dint     = 0.0_dp    ! internal-pivot step
    real(dp) :: fdick    = 0.0_dp    ! fraction of Dickman moves
    real(dp) :: frept    = 0.0_dp    ! fraction of reptation moves
    real(dp) :: beps     = 0.0_dp    ! bead–bead well depth (runt only)
    real(dp) :: bbeps    = 0.0_dp    ! bead–box well depth  (runt only)
    real(dp) :: rsp      = 0.0_dp    ! sphere radius        (runt only)
    real(dp) :: rcut     = -1.0_dp   ! cutoff radius: <=0 → exact; >0 → cutoff
  end type params_t

contains

  ! ------------------------------------------------------------------
  ! read_runt_inp — parse runt.inp
  !
  ! File layout (one value per line, list-directed):
  !   1  ncon
  !   2  nskip
  !   3  bsz
  !   4  dlr
  !   5  dint
  !   6  fdick
  !   7  frept
  !   8  beps
  !   9  bbeps
  !  10  rsp
  !  11  rcut  (OPTIONAL — if absent/EOF leave rcut = -1.0_dp)
  ! ------------------------------------------------------------------
  subroutine read_runt_inp(path, p)
    character(len=*), intent(in)  :: path
    type(params_t),   intent(out) :: p

    integer  :: iu, ios
    real(dp) :: rtmp
    character(len=256) :: msg

    p%rcut = -1.0_dp   ! default: exact (no cutoff)

    open(newunit=iu, file=trim(path), status='old', action='read', &
         iostat=ios, iomsg=msg)
    if (ios /= 0) error stop "polymc_config: cannot open " // trim(path)

    read(iu, *) p%ncon
    read(iu, *) p%nskip
    read(iu, *) p%bsz
    read(iu, *) p%dlr
    read(iu, *) p%dint
    read(iu, *) p%fdick
    read(iu, *) p%frept
    read(iu, *) p%beps
    read(iu, *) p%bbeps
    read(iu, *) p%rsp

    ! Optional 11th line: rcut.
    ! A blank line produces iostat=EOR (nonzero); EOF also nonzero.
    ! Either way we keep the default.
    read(iu, *, iostat=ios) rtmp
    if (ios == 0) p%rcut = rtmp

    close(iu)
  end subroutine read_runt_inp

  ! ------------------------------------------------------------------
  ! read_polyfj_inp — parse polyfj.inp
  !
  ! File layout (one value per line):
  !   1  ncon
  !   2  nskip
  !   3  bsz
  !   4  dlr
  !   5  dint
  !   6  fdick
  !   7  frept
  ! (No beps / bbeps / rsp / rcut in polyfj runs.)
  ! ------------------------------------------------------------------
  subroutine read_polyfj_inp(path, p)
    character(len=*), intent(in)  :: path
    type(params_t),   intent(out) :: p

    integer :: iu, ios
    character(len=256) :: msg

    p%rcut = -1.0_dp

    open(newunit=iu, file=trim(path), status='old', action='read', &
         iostat=ios, iomsg=msg)
    if (ios /= 0) error stop "polymc_config: cannot open " // trim(path)

    read(iu, *) p%ncon
    read(iu, *) p%nskip
    read(iu, *) p%bsz
    read(iu, *) p%dlr
    read(iu, *) p%dint
    read(iu, *) p%fdick
    read(iu, *) p%frept

    close(iu)
  end subroutine read_polyfj_inp

  ! ------------------------------------------------------------------
  ! read_ic — read an initial/restart configuration file
  !
  ! Common header (both runt and polyfj):
  !   Record 1 (A):  blank — consumed with (A) format  [spec says PFC value]
  !     -> Actually PFC is a value on line 1; we read it with (*) then...
  !
  ! Exact layout observed in fixture files:
  !
  !   runt.ic:
  !     line 1: PFC (e.g. 0.2000)        — read(*) into dummy_real
  !     line 2: blank                    — consume with ('(A)')
  !     line 3: NMOL                     — read(*) sys%nmol
  !     line 4: N                        — read(*) sys%n
  !     line 5: blank                    — consume with ('(A)')
  !     line 6: blank                    — consume with ('(A)')
  !     line 7: AL                       — read(*) sys%box%L
  !     lines 8..: bead records I J X Y Z
  !
  !   polyfj.ic:
  !     line 1: PFC                      — read(*) into dummy_real
  !     line 2: blank                    — consume with ('(A)')
  !     line 3: NMOL                     — read(*) sys%nmol
  !     line 4: N                        — read(*) sys%n
  !     line 5: blank                    — consume with ('(A)')
  !     line 6: RSP                      — read(*) sys%rsp
  !     line 7: AL                       — read(*) sys%box%L
  !     line 8: blank                    — consume with ('(A)')
  !     lines 9..: bead records I J X Y Z
  !
  ! has_rsp=.true.  → polyfj layout (RSP present, blank after AL)
  ! has_rsp=.false. → runt layout    (two blanks before AL, no blank after AL)
  ! ------------------------------------------------------------------
  subroutine read_ic(path, sys, has_rsp)
    character(len=*), intent(in)  :: path
    type(system_t),   intent(out) :: sys
    logical,          intent(in)  :: has_rsp

    integer  :: iu, ios, k, ibead_mol, ibead_chain
    real(dp) :: dummy_real
    character(len=256) :: dummy_str, msg

    open(newunit=iu, file=trim(path), status='old', action='read', &
         iostat=ios, iomsg=msg)
    if (ios /= 0) error stop "polymc_config: cannot open " // trim(path)

    ! Line 1: PFC (read value, discard)
    read(iu, *) dummy_real

    ! Line 2: blank separator
    read(iu, '(A)') dummy_str

    ! Line 3: NMOL
    read(iu, *) sys%nmol

    ! Line 4: N (beads per molecule)
    read(iu, *) sys%n

    ! Line 5: blank separator
    read(iu, '(A)') dummy_str

    if (has_rsp) then
      ! polyfj layout
      ! Line 6: RSP
      read(iu, *) sys%rsp
      ! Line 7: AL (box side length)
      read(iu, *) sys%box%L
      ! Line 8: blank separator before beads
      read(iu, '(A)') dummy_str
    else
      ! runt layout: two blanks before AL (lines 5+6 consumed above/below)
      ! Line 6: another blank
      read(iu, '(A)') dummy_str
      ! Line 7: AL
      read(iu, *) sys%box%L
      ! No blank between AL and bead records in runt.ic
    end if

    ! Allocate storage for all beads
    sys%nbeads = sys%nmol * sys%n
    allocate(sys%x(sys%nbeads), sys%y(sys%nbeads), sys%z(sys%nbeads))
    allocate(sys%mol_of(sys%nbeads))

    ! Read bead records: I  J  X  Y  Z
    ! I = molecule index, J = bead-within-molecule index (1..n)
    do k = 1, sys%nbeads
      read(iu, *) ibead_mol, ibead_chain, &
                  sys%x(k), sys%y(k), sys%z(k)
      sys%mol_of(k) = ibead_mol
    end do

    close(iu)
  end subroutine read_ic

  ! ------------------------------------------------------------------
  ! write_fc — write a restart configuration file
  !
  ! Output layout (re-readable by read_ic with has_rsp=.false.):
  !   line 1: pfc
  !   line 2: blank
  !   line 3: nmol
  !   line 4: n
  !   line 5: blank
  !   line 6: blank
  !   line 7: box%L
  !   lines 8..: I  J  X  Y  Z  (one per bead, 1-based J within molecule)
  !
  ! Spacing is simple and consistent; values use full double-precision
  ! exponential notation so they can be read back without loss.
  ! ------------------------------------------------------------------
  subroutine write_fc(path, sys, pfc)
    character(len=*), intent(in) :: path
    type(system_t),   intent(in) :: sys
    real(dp),         intent(in) :: pfc

    integer :: iu, k, mol_idx, bead_in_mol

    open(newunit=iu, file=trim(path), status='replace', action='write')

    ! Header
    write(iu, '(ES16.10)') pfc
    write(iu, *)                ! blank line 2
    write(iu, '(I0)')  sys%nmol
    write(iu, '(I0)')  sys%n
    write(iu, *)                ! blank line 5
    write(iu, *)                ! blank line 6
    write(iu, '(ES22.10)') sys%box%L

    ! Bead records: I J X Y Z
    ! mol_of(k) holds molecule index; bead_in_mol is 1..n within that molecule
    bead_in_mol = 0
    do k = 1, sys%nbeads
      mol_idx = sys%mol_of(k)

      ! Assumption: beads belonging to the same molecule are stored
      ! contiguously in the arrays (i.e., all beads of molecule i appear
      ! before any bead of molecule i+1).  The J counter below resets
      ! whenever mol_of(k) /= mol_of(k-1); it will produce wrong J values
      ! if this contiguity invariant is violated.
      if (k == 1) then
        bead_in_mol = 1
      else if (sys%mol_of(k) /= sys%mol_of(k-1)) then
        bead_in_mol = 1
      else
        bead_in_mol = bead_in_mol + 1
      end if

      write(iu, '(I6, 1X, I6, 3(1X, ES22.10))') &
            mol_idx, bead_in_mol, sys%x(k), sys%y(k), sys%z(k)
    end do

    close(iu)
  end subroutine write_fc

end module polymc_config
