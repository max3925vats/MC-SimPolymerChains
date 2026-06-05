! polymc_overlap.f90 — hard-core overlap checks for beads and central sphere
!
! Provides:
!   bead_overlaps  — is a trial point within sigma(=1) of any bead in a
!                    different molecule? Uses the cell list for O(1) candidate
!                    gathering, then applies min-image distance filter.
!   sphere_overlaps — is a trial point within (rsp + 0.5) of the origin?
!                     Uses min-image convention so the box wrapping is respected.

module polymc_overlap
  use polymc_kinds,     only: dp
  use polymc_box,       only: box_t, min_image
  use polymc_cell_list, only: cell_list_t, cl_neighbours
  implicit none
  private
  public :: bead_overlaps, sphere_overlaps

contains

  ! ------------------------------------------------------------------
  !> Check whether the trial point (px, py, pz) overlaps any bead that
  !  does NOT belong to exclude_mol.
  !
  !  Two beads "overlap" when their centre-to-centre distance (min-image)
  !  is < sigma = 1.0.  We compare squared distances against 1.0_dp.
  !
  !  Candidate beads are gathered from the 27-cell stencil via cl_neighbours
  !  (which returns all beads in those cells without distance filtering).
  !  The caller-side buffer `idx` is declared as an automatic array of size
  !  nbeads — safe for nbeads up to ~1e4 on the stack.
  !
  !  Early-returns .true. as soon as the first true overlap is found.
  ! ------------------------------------------------------------------
  logical function bead_overlaps(box, cl, x, y, z, nbeads, mol_of, &
                                  px, py, pz, exclude_mol)
    type(box_t),        intent(in) :: box
    type(cell_list_t),  intent(in) :: cl
    real(dp),           intent(in) :: x(:), y(:), z(:)  ! stored bead positions
    integer,            intent(in) :: nbeads
    integer,            intent(in) :: mol_of(nbeads)     ! molecule id per bead
    real(dp),           intent(in) :: px, py, pz         ! trial point
    integer,            intent(in) :: exclude_mol        ! molecule to skip

    integer :: idx(nbeads)  ! automatic stack buffer for candidate ids
    integer :: nidx         ! number of candidates returned
    integer :: k, j
    real(dp) :: dx, dy, dz, d2

    bead_overlaps = .false.

    ! Gather candidates from the 27-cell neighbourhood
    call cl_neighbours(cl, box, px, py, pz, idx, nidx)

    do k = 1, nidx
      j = idx(k)
      ! Skip beads belonging to the excluded molecule
      if (mol_of(j) == exclude_mol) cycle

      ! Minimum-image displacement components
      dx = min_image(x(j) - px, box%L)
      dy = min_image(y(j) - py, box%L)
      dz = min_image(z(j) - pz, box%L)

      d2 = dx*dx + dy*dy + dz*dz

      ! Overlap: centre-to-centre distance < sigma = 1  =>  d^2 < 1
      if (d2 < 1.0_dp) then
        bead_overlaps = .true.
        return  ! early exit on first overlap
      end if
    end do

  end function bead_overlaps

  ! ------------------------------------------------------------------
  !> Check whether the trial point (px, py, pz) overlaps the central
  !  hard sphere of radius rsp.
  !
  !  Bead radius is 0.5, so contact occurs at distance rsp + 0.5 from
  !  the origin.  Overlap if the min-image distance is strictly less
  !  than (rsp + 0.5).
  ! ------------------------------------------------------------------
  logical function sphere_overlaps(box, rsp, px, py, pz)
    type(box_t), intent(in) :: box
    real(dp),    intent(in) :: rsp       ! radius of the central hard sphere
    real(dp),    intent(in) :: px, py, pz

    real(dp) :: dx, dy, dz, d2, contact2

    ! Min-image displacement from the origin (origin is inside the box)
    dx = min_image(px, box%L)
    dy = min_image(py, box%L)
    dz = min_image(pz, box%L)

    d2       = dx*dx + dy*dy + dz*dz
    contact2 = (rsp + 0.5_dp)**2

    sphere_overlaps = (d2 < contact2)

  end function sphere_overlaps

end module polymc_overlap
