! polymc_cell_list.f90 — standard periodic linked-cell neighbour search
!
! Provides:
!   cell_list_t   — holds the linked-cell data structure
!   cl_init       — initialise cell grid from box and cutoff
!   cl_build      — (re)build head/nextp from current bead positions
!   cl_neighbours — gather all bead ids whose cell lies within the 27-cell
!                   stencil around a query point; caller applies distance filter

module polymc_cell_list
  use polymc_kinds, only: dp
  use polymc_box,   only: box_t, cell_index
  implicit none
  private
  public :: cell_list_t, cl_init, cl_build, cl_neighbours

  !> Linked-cell list data structure.
  !
  !  head(c)    — id of the first bead in cell c (0 = cell is empty).
  !  nextp(i)   — id of the next bead after bead i in the same cell (0 = none).
  !  Linear cell index: id = ix + (iy-1)*mc + (iz-1)*mc*mc  (1-based)
  type :: cell_list_t
    integer  :: mc        = 0       ! cells per dimension
    real(dp) :: cell_size = 0.0_dp  ! L / mc
    integer, allocatable :: head(:)   ! size mc**3; 0 = empty
    integer, allocatable :: nextp(:)  ! size nbeads;  0 = end of chain
  end type cell_list_t

contains

  ! ------------------------------------------------------------------
  !> Initialise the cell grid.
  !  mc = max(3, floor(L/r_cut));  cell_size = L/mc;  allocate head(mc**3).
  ! ------------------------------------------------------------------
  subroutine cl_init(cl, box, r_cut)
    type(cell_list_t), intent(inout) :: cl
    type(box_t),       intent(in)    :: box
    real(dp),          intent(in)    :: r_cut

    cl%mc        = max(3, int(box%L / r_cut))   ! floor via int()
    cl%cell_size = box%L / real(cl%mc, dp)

    if (allocated(cl%head)) deallocate(cl%head)
    allocate(cl%head(cl%mc**3))
    cl%head = 0
  end subroutine cl_init

  ! ------------------------------------------------------------------
  !> Build the linked-cell list from bead positions.
  !  Allocates/re-sizes nextp(nbeads); resets head to 0; then inserts
  !  each bead at the head of its cell's chain (O(nbeads)).
  ! ------------------------------------------------------------------
  subroutine cl_build(cl, box, x, y, z, nbeads)
    type(cell_list_t), intent(inout) :: cl
    type(box_t),       intent(in)    :: box
    real(dp),          intent(in)    :: x(:), y(:), z(:)
    integer,           intent(in)    :: nbeads
    integer :: i, ix, iy, iz, cid

    ! (Re-)allocate nextp for this bead count
    if (allocated(cl%nextp)) deallocate(cl%nextp)
    allocate(cl%nextp(nbeads))

    ! Reset all cell heads to empty
    cl%head  = 0
    cl%nextp = 0

    ! Insert beads: prepend each bead to its cell's linked list
    do i = 1, nbeads
      call cell_index(box, cl%mc, x(i), y(i), z(i), ix, iy, iz)
      cid = ix + (iy - 1)*cl%mc + (iz - 1)*cl%mc*cl%mc
      cl%nextp(i) = cl%head(cid)   ! point to current head
      cl%head(cid) = i             ! become new head
    end do
  end subroutine cl_build

  ! ------------------------------------------------------------------
  !> Gather bead ids from the 27-cell stencil around query point (px,py,pz).
  !
  !  idx   — caller-preallocated integer array (size >= nbeads)
  !  nidx  — on exit: number of valid entries in idx(1:nidx)
  !
  !  Periodic wrap: cell coordinate 0  -> mc;  mc+1 -> 1.
  !  No distance filtering; caller is responsible for that.
  ! ------------------------------------------------------------------
  subroutine cl_neighbours(cl, box, px, py, pz, idx, nidx)
    type(cell_list_t), intent(in)    :: cl
    type(box_t),       intent(in)    :: box
    real(dp),          intent(in)    :: px, py, pz
    integer,           intent(inout) :: idx(:)
    integer,           intent(out)   :: nidx
    integer :: qix, qiy, qiz      ! cell of query point
    integer :: dix, diy, diz      ! stencil offsets (-1, 0, +1)
    integer :: nx, ny, nz         ! neighbour cell indices (before wrap)
    integer :: wx, wy, wz         ! wrapped cell indices
    integer :: cid, bead

    nidx = 0

    ! Cell containing the query point
    call cell_index(box, cl%mc, px, py, pz, qix, qiy, qiz)

    ! Loop over the 3x3x3 stencil
    do diz = -1, 1
      nz = qiz + diz
      ! Periodic wrap: map 0 -> mc,  mc+1 -> 1
      if (nz < 1)     then; wz = nz + cl%mc
      else if (nz > cl%mc) then; wz = nz - cl%mc
      else;                       wz = nz
      end if

      do diy = -1, 1
        ny = qiy + diy
        if (ny < 1)     then; wy = ny + cl%mc
        else if (ny > cl%mc) then; wy = ny - cl%mc
        else;                       wy = ny
        end if

        do dix = -1, 1
          nx = qix + dix
          if (nx < 1)     then; wx = nx + cl%mc
          else if (nx > cl%mc) then; wx = nx - cl%mc
          else;                       wx = nx
          end if

          ! Linear cell id (1-based)
          cid = wx + (wy - 1)*cl%mc + (wz - 1)*cl%mc*cl%mc

          ! Walk the linked list for this cell
          bead = cl%head(cid)
          do while (bead /= 0)
            nidx = nidx + 1
            idx(nidx) = bead
            bead = cl%nextp(bead)
          end do

        end do
      end do
    end do
  end subroutine cl_neighbours

end module polymc_cell_list
