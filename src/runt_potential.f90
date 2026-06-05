! runt_potential.f90 — Yukawa (screened-Coulomb) pair and sphere-bead potentials
!
! Physical setup (reduced units: sigma = 1, bead radius = 0.5):
!   u_pair:   bead–bead Yukawa, measured from contact at r = 1.
!   u_sphere: bead–sphere Yukawa, measured from contact at r = rsp + 0.5.
!
! Both functions accept an optional hard truncation rcut:
!   rcut <= 0  →  exact (no truncation; faithful to the legacy runt behaviour)
!   rcut >  0  →  return 0.0 if r > rcut
!
! The sign convention follows the legacy BATT/BBATT coefficients:
!   beps  = -eps_ff   (negative for an attractive fluid-fluid interaction)
!   bbeps = -eps_sf   (negative for an attractive sphere-fluid interaction)
! so the Yukawa energy u = beps * exp(-kappa*(r - r_contact)) / r is negative
! (attractive) when beps < 0.
!
! kappa = 2.5 (inverse range, hard-coded to match legacy constant 2.5D0).

module runt_potential
  use polymc_kinds, only: dp
  implicit none
  private
  public :: u_pair, u_sphere

  ! Yukawa inverse-range parameter — matches legacy 2.5D0
  real(dp), parameter :: kappa = 2.5_dp

contains

  ! ------------------------------------------------------------------
  !> Yukawa pair potential between two beads.
  !
  !  Contact distance between bead centres: sigma = 1 (bead radius 0.5 each).
  !  Functional form:  u = beps * exp(-kappa*(r - 1)) / r
  !
  !  @param r     Centre-to-centre distance (must be >= 0)
  !  @param beps  Code coefficient; typically beps = -eps_ff < 0 for attraction
  !  @param rcut  Truncation radius: rcut <= 0 means no truncation;
  !               rcut > 0 returns 0 if r > rcut
  ! ------------------------------------------------------------------
  elemental pure real(dp) function u_pair(r, beps, rcut)
    real(dp), intent(in) :: r
    real(dp), intent(in) :: beps
    real(dp), intent(in) :: rcut

    ! Apply optional hard truncation
    if (rcut > 0.0_dp .and. r > rcut) then
      u_pair = 0.0_dp
      return
    end if

    ! Yukawa measured from contact distance 1.0 (= sigma)
    u_pair = beps * exp(-kappa * (r - 1.0_dp)) / r

  end function u_pair

  ! ------------------------------------------------------------------
  !> Yukawa potential between a bead and the central hard sphere.
  !
  !  Contact occurs when the bead centre is at distance rsp + 0.5 from
  !  the sphere centre (sphere radius rsp, bead radius 0.5).
  !  Inside the hard core (r < rsp + 0.5) the potential is 0 (infinite wall,
  !  but the caller is expected not to call here with overlapping configurations;
  !  we return 0 rather than infinity for defensive correctness).
  !
  !  Functional form (outside core): u = bbeps * exp(-kappa*(r-(rsp+0.5))) / r
  !
  !  @param r      Bead-centre to sphere-centre distance
  !  @param bbeps  Code coefficient; typically bbeps = -eps_sf < 0 for attraction
  !  @param rsp    Radius of the central hard sphere
  !  @param rcut   Truncation radius: rcut <= 0 means no truncation;
  !                rcut > 0 returns 0 if r > rcut
  ! ------------------------------------------------------------------
  elemental pure real(dp) function u_sphere(r, bbeps, rsp, rcut)
    real(dp), intent(in) :: r
    real(dp), intent(in) :: bbeps
    real(dp), intent(in) :: rsp
    real(dp), intent(in) :: rcut

    real(dp) :: r_contact

    r_contact = rsp + 0.5_dp

    ! Return 0 inside the hard core
    if (r < r_contact) then
      u_sphere = 0.0_dp
      return
    end if

    ! Apply optional hard truncation
    if (rcut > 0.0_dp .and. r > rcut) then
      u_sphere = 0.0_dp
      return
    end if

    ! Yukawa measured from contact distance rsp + 0.5
    u_sphere = bbeps * exp(-kappa * (r - r_contact)) / r

  end function u_sphere

end module runt_potential
