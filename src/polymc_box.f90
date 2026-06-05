! polymc_box.f90 — simulation box geometry
!
! Provides:
!   box_t       — derived type holding the cubic box side length L
!   min_image   — elemental minimum-image convention: d - L*nint(d/L)
!   cell_index  — maps a position into a 1-based cell coordinate in 1..mc

module polymc_box
  use polymc_kinds, only: dp
  implicit none
  private
  public :: box_t, min_image, cell_index

  !> Cubic simulation box; L is in the same units as particle positions (sigma)
  type :: box_t
    real(dp) :: L = 0.0_dp
  end type box_t

contains

  !> Minimum-image displacement: maps d into (-L/2, L/2]
  elemental function min_image(d, L) result(m)
    real(dp), intent(in) :: d   ! raw displacement component
    real(dp), intent(in) :: L   ! box side length
    real(dp) :: m
    m = d - L * nint(d / L)
  end function min_image

  !> Map a position (x,y,z) into 1-based cell indices ix,iy,iz in 1..mc.
  !
  !  Uses modulo to fold positions back into [0, L) before binning, then
  !  clamps to mc to guard against the exact-L floating-point edge case.
  pure subroutine cell_index(box, mc, x, y, z, ix, iy, iz)
    type(box_t), intent(in)  :: box
    integer,     intent(in)  :: mc           ! number of cells per dimension
    real(dp),    intent(in)  :: x, y, z      ! position (0 <= coord < L)
    integer,     intent(out) :: ix, iy, iz   ! 1-based cell indices

    ! int(modulo(coord,L)/L * mc) gives 0..mc-1; add 1 for 1-based indexing
    ix = 1 + int(modulo(x, box%L) / box%L * real(mc, dp))
    iy = 1 + int(modulo(y, box%L) / box%L * real(mc, dp))
    iz = 1 + int(modulo(z, box%L) / box%L * real(mc, dp))

    ! Clamp: floating-point rounding can push the index one above mc
    ix = min(ix, mc)
    iy = min(iy, mc)
    iz = min(iz, mc)
  end subroutine cell_index

end module polymc_box
