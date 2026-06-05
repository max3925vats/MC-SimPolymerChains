! runt_energy.f90 — inter-molecular, sphere, and intra-molecular Yukawa energies
!
! Two public functions:
!
!   bead_inter_energy(sys, p, cl, imol, px, py, pz)
!     — Yukawa energy of a trial bead at (px,py,pz) in molecule imol
!       due to all beads in OTHER molecules.
!     — EXACT mode (p%rcut <= 0): loops ALL nbeads, skips same-molecule,
!       calls u_pair with rcut=-1 (no truncation).
!     — CUTOFF mode (p%rcut > 0): uses cl_neighbours to get the 27-cell
!       candidate set, skips same-molecule, distance-filters to rcut, sums u_pair.
!
!   total_energy(sys, p, cl)
!     — Full system energy:
!         inter   = 0.5 * sum_k bead_inter_energy(mol_of(k), pos_k)
!         sphere  = sum_k u_sphere(r_k, bbeps, rsp, rcut)  [r_k = min-image dist to origin]
!         intra   = sum over chains, non-bonded pairs (j,k) with j>=3, k<=j-2
!       return inter + sphere + intra
!
! Faithful to legacy ENERCAL in runt-moves.f (exact mode).
! In CUTOFF mode every term is consistently truncated at p%rcut.

module runt_energy
  use polymc_kinds,     only: dp
  use polymc_box,       only: min_image
  use polymc_cell_list, only: cell_list_t, cl_neighbours
  use polymc_config,    only: system_t, params_t
  use runt_potential,   only: u_pair, u_sphere
  implicit none
  private
  public :: bead_inter_energy, total_energy

contains

  ! ------------------------------------------------------------------
  !> Inter-molecular Yukawa energy of a trial bead at (px,py,pz).
  !
  !  @param sys   System state (positions, mol_of, box, etc.)
  !  @param p     Run parameters (beps, rcut)
  !  @param cl    Cell list (used only when p%rcut > 0)
  !  @param imol  Molecule index the trial bead belongs to
  !  @param px,py,pz  Trial bead position (sigma units)
  !  @return  Sum of u_pair(r, beps, rcut) over beads j with mol_of(j) /= imol
  ! ------------------------------------------------------------------
  real(dp) function bead_inter_energy(sys, p, cl, imol, px, py, pz)
    type(system_t),    intent(in) :: sys
    type(params_t),    intent(in) :: p
    type(cell_list_t), intent(in) :: cl
    integer,           intent(in) :: imol
    real(dp),          intent(in) :: px, py, pz

    integer  :: j, nidx
    real(dp) :: dx, dy, dz, r
    ! Stack array for cell-list candidate indices
    integer :: idx(sys%nbeads)

    bead_inter_energy = 0.0_dp

    if (p%rcut <= 0.0_dp) then
      ! EXACT mode: loop over all beads, skip same molecule
      do j = 1, sys%nbeads
        if (sys%mol_of(j) /= imol) then
          dx = min_image(sys%x(j) - px, sys%box%L)
          dy = min_image(sys%y(j) - py, sys%box%L)
          dz = min_image(sys%z(j) - pz, sys%box%L)
          r  = sqrt(dx*dx + dy*dy + dz*dz)
          ! Pass rcut=-1 so u_pair does no truncation
          bead_inter_energy = bead_inter_energy + u_pair(r, p%beps, -1.0_dp)
        end if
      end do

    else
      ! CUTOFF mode: gather 27-cell stencil candidates, then filter
      call cl_neighbours(cl, sys%box, px, py, pz, idx, nidx)
      do j = 1, nidx
        if (sys%mol_of(idx(j)) /= imol) then
          dx = min_image(sys%x(idx(j)) - px, sys%box%L)
          dy = min_image(sys%y(idx(j)) - py, sys%box%L)
          dz = min_image(sys%z(idx(j)) - pz, sys%box%L)
          r  = sqrt(dx*dx + dy*dy + dz*dz)
          ! u_pair handles the distance cutoff: returns 0 if r > rcut
          bead_inter_energy = bead_inter_energy + u_pair(r, p%beps, p%rcut)
        end if
      end do

    end if

  end function bead_inter_energy

  ! ------------------------------------------------------------------
  !> Full system potential energy.
  !
  !  Components (faithful to legacy ENERCAL):
  !    inter  = 0.5 * sum_k bead_inter_energy(mol_of(k), pos_k)
  !    sphere = sum_k u_sphere(r_k, bbeps, rsp, rcut)
  !             where r_k = sqrt(sum min_image(coord_k, L)^2)  [dist to origin]
  !    intra  = sum over mol i, j=3..n, k=1..j-2  of u_pair(r_jk, beps, rcut)
  !             global indices: bead j in mol i -> (i-1)*n + j
  !
  !  @param sys  System state
  !  @param p    Run parameters (beps, bbeps, rcut)
  !  @param cl   Cell list (used only in cutoff mode)
  !  @return  inter + sphere + intra
  ! ------------------------------------------------------------------
  real(dp) function total_energy(sys, p, cl)
    type(system_t),    intent(in) :: sys
    type(params_t),    intent(in) :: p
    type(cell_list_t), intent(in) :: cl

    integer  :: k, i, j, ij, ik
    real(dp) :: e_inter, e_sphere, e_intra
    real(dp) :: dx, dy, dz, r, r_k

    ! --- Inter-molecular contribution (factor 0.5 avoids double-counting) ---
    e_inter = 0.0_dp
    do k = 1, sys%nbeads
      e_inter = e_inter + bead_inter_energy(sys, p, cl, &
                            sys%mol_of(k), sys%x(k), sys%y(k), sys%z(k))
    end do
    e_inter = 0.5_dp * e_inter

    ! --- Sphere (central bead) contribution ---
    ! r_k = minimum-image distance from bead k to the origin (sphere centre)
    e_sphere = 0.0_dp
    do k = 1, sys%nbeads
      dx  = min_image(sys%x(k), sys%box%L)
      dy  = min_image(sys%y(k), sys%box%L)
      dz  = min_image(sys%z(k), sys%box%L)
      r_k = sqrt(dx*dx + dy*dy + dz*dz)
      e_sphere = e_sphere + u_sphere(r_k, p%bbeps, sys%rsp, p%rcut)
    end do

    ! --- Intra-molecular non-bonded contribution ---
    ! Pairs (j, k) within each chain that are separated by >= 2 bonds:
    !   j = 3..n, k = 1..j-2
    ! Global bead index for position j in chain i: (i-1)*n + j
    e_intra = 0.0_dp
    do i = 1, sys%nmol
      do j = 3, sys%n
        do k = 1, j - 2
          ij = (i - 1)*sys%n + j   ! global index of bead j in chain i
          ik = (i - 1)*sys%n + k   ! global index of bead k in chain i
          dx = min_image(sys%x(ij) - sys%x(ik), sys%box%L)
          dy = min_image(sys%y(ij) - sys%y(ik), sys%box%L)
          dz = min_image(sys%z(ij) - sys%z(ik), sys%box%L)
          r  = sqrt(dx*dx + dy*dy + dz*dz)
          e_intra = e_intra + u_pair(r, p%beps, p%rcut)
        end do
      end do
    end do

    total_energy = e_inter + e_sphere + e_intra

  end function total_energy

end module runt_energy
