! polymc_chain.f90 — chain-level coordinate operations
!
! Provides in-place reversal of a polymer chain's bead ordering.
! Coordinates are stored in flat arrays x(:),y(:),z(:); chain i of
! length n occupies indices ibase+1 .. ibase+n where ibase=(i-1)*n.
!
! For a symmetric freely-jointed chain, reversal is a pure relabelling:
! all pairwise distances (and thus energy / overlap / density) are
! unchanged.  This lets reptation / CCB algorithms grow from either end.

module polymc_chain
  use polymc_kinds, only: dp
  implicit none
  private

  public :: chain_reverse

contains

  !> Reverse the n beads belonging to one chain, in place.
  !>
  !> Args:
  !>   x,y,z  - flat bead coordinate arrays (all chains); intent(inout)
  !>   ibase  - index offset: chain occupies ibase+1 .. ibase+n
  !>   n      - number of beads in this chain; n<=1 is a no-op
  !>
  !> Algorithm: swap bead ibase+j with bead ibase+(n+1-j) for j = 1 .. n/2.
  subroutine chain_reverse(x, y, z, ibase, n)
    real(dp), intent(inout) :: x(:), y(:), z(:)
    integer,  intent(in)    :: ibase, n

    integer  :: j, lo, hi
    real(dp) :: tmp_x, tmp_y, tmp_z

    do j = 1, n / 2
      lo = ibase + j
      hi = ibase + (n + 1 - j)
      ! Swap x
      tmp_x  = x(lo);  x(lo) = x(hi);  x(hi) = tmp_x
      ! Swap y
      tmp_y  = y(lo);  y(lo) = y(hi);  y(hi) = tmp_y
      ! Swap z
      tmp_z  = z(lo);  z(lo) = z(hi);  z(hi) = tmp_z
    end do
  end subroutine chain_reverse

end module polymc_chain
