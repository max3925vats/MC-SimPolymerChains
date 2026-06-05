! polyfj_observables.f90 — structural observable accumulators for polyfj
!
! Faithfully ports the accumulation math from legacy polyfj/polyfj.f,
! lines ~259-394.  Positions are expected in real sigma units (unlike
! the legacy code which stored reduced coords divided by AL).
!
! Public interface:
!   polyfj_obs_t  — derived type holding all accumulator arrays
!   obs_init      — allocate and zero all accumulators
!   obs_sample    — accumulate one snapshot over all chains
!
! Unit-conversion notes vs. legacy:
!   Legacy stores coords in reduced units (x_red = x_sigma / AL).
!   Here coords are already in sigma, so:
!     - density binning uses bsz directly (not bsz/AL)
!     - segmental-order accumulates plain cos(theta) with NO AL^2 factor
!     - Rg^2, Re^2 are in sigma^2 (no AL^2 multiplication needed)
!     - shape semi-axes: see note below
!
! Shape semi-axis accumulators — two parallel sets:
!
!   CORRECTED (default, a_sum / b_sum / c_sum):
!     Diagonalises the GYRATION tensor G_{ab} = (1/N)*sum_j dr_a*dr_b
!     and accumulates sqrt(eigenvalue) sorted descending (a>=b>=c).
!     The driver finalises as a_out = a_sum/nc (proper mean of per-chain values).
!
!   LEGACY-FAITHFUL (acig_sum / bcig_sum / ccig_sum):
!     Builds the INERTIA tensor I_{ab} (off-diagonal terms negated),
!     diagonalises it, and applies the legacy formula (polyfj.f lines 372-374):
!       cigai = sqrt(max(2.5*(eig2+eig3-eig1)/n, 0))   ! eig1<=eig2<=eig3
!       cigbi = sqrt(max(2.5*(eig1+eig3-eig2)/n, 0))
!       cigci = sqrt(max(2.5*(eig2+eig1-eig3)/n, 0))
!     The driver finalises as acig_out = sqrt(acig_sum/nc) to match legacy's
!     sqrt-of-mean convention (ACIG(K) = sum cigai*AL; output = sqrt(ACIG/NC)*AL).
!     The max(...,0) guard prevents NaN on degenerate (rod-like) chains where
!     legacy would have taken the sqrt of a negative value.
!     Note: the legacy *AL scaling factor is dropped — coords are in sigma here.
!
! Segmental order:
!   Accumulates plain cos(theta) (no AL^2 factor).
!   The driver computes S = (3*<cos>-1)/2 at write time.

module polyfj_observables
  use polymc_kinds,   only: dp, i8
  use polymc_box,     only: min_image
  use polymc_config,  only: system_t
  use polymc_binning, only: hist_t, hist_init, hist_add
  use polymc_linalg,  only: jacobi3
  implicit none
  private
  public :: polyfj_obs_t, obs_init, obs_sample

  type :: polyfj_obs_t
    ! ----------------------------------------------------------------
    ! Density profiles (binned by distance-to-origin of each bead)
    ! ----------------------------------------------------------------
    type(hist_t) :: h_all   ! all beads
    type(hist_t) :: h_end   ! chain-end beads  (j==1 or j==n)
    type(hist_t) :: h_mid   ! middle beads     (j==nm1 or j==nm2)

    ! ----------------------------------------------------------------
    ! Segmental order parameter (binned by half-bond-midpoint distance)
    ! ----------------------------------------------------------------
    integer :: nbin = 0
    real(dp) :: bsz = 0.0_dp
    real(dp), allocatable :: seg_cos_sum(:)  ! sum of cos(theta) per bin
    integer(i8), allocatable :: seg_cnt(:)   ! count per bin

    ! ----------------------------------------------------------------
    ! CoM-resolved chain shape / size (binned by CoM distance-to-origin)
    ! ----------------------------------------------------------------
    integer :: nbinc = 0
    real(dp) :: bszc = 0.0_dp
    integer(i8), allocatable :: nc(:)        ! chain count per bin
    real(dp), allocatable :: rg2_sum(:)      ! sum of Rg^2 [sigma^2]
    real(dp), allocatable :: re2_sum(:)      ! sum of Re^2 [sigma^2]
    real(dp), allocatable :: a_sum(:)        ! sum of largest  sqrt(gyration eval)  [corrected]
    real(dp), allocatable :: b_sum(:)        ! sum of middle   sqrt(gyration eval)  [corrected]
    real(dp), allocatable :: c_sum(:)        ! sum of smallest sqrt(gyration eval)  [corrected]

    ! Legacy inertia-tensor semi-axes (for legacy-faithful g2fj.out emission).
    ! These mirror ACIG/BCIG/CCIG from polyfj.f lines 390-392.
    ! Driver finalises as sqrt(acig_sum/nc), matching legacy's sqrt-of-mean convention.
    real(dp), allocatable :: acig_sum(:)     ! sum of cigai (legacy inertia semi-axis, sigma)
    real(dp), allocatable :: bcig_sum(:)     ! sum of cigbi
    real(dp), allocatable :: ccig_sum(:)     ! sum of cigci
  end type polyfj_obs_t

contains

  ! ------------------------------------------------------------------
  !> obs_init — allocate and zero all accumulator arrays.
  !>
  !> @param obs   [out] observable accumulator
  !> @param bsz   [in]  bin width for density/segorder histograms [sigma]
  !> @param bszc  [in]  bin width for CoM-resolved histograms [sigma]
  !> @param rmax  [in]  maximum relevant distance = box%L/2 [sigma]
  ! ------------------------------------------------------------------
  subroutine obs_init(obs, bsz, bszc, rmax)
    type(polyfj_obs_t), intent(inout) :: obs
    real(dp), intent(in) :: bsz
    real(dp), intent(in) :: bszc
    real(dp), intent(in) :: rmax

    integer :: nb, nbc

    ! Number of bins (matches legacy: NBIN = INT(AL/2/BSZ), NBINC = INT(AL/2/BSZCM))
    nb  = int(rmax / bsz)
    nbc = int(rmax / bszc)

    ! Store for use in obs_sample
    obs%bsz  = bsz
    obs%bszc = bszc
    obs%nbin  = nb
    obs%nbinc = nbc

    ! Density histograms (hist_init handles allocation and zeroing)
    call hist_init(obs%h_all, bsz, rmax)
    call hist_init(obs%h_end, bsz, rmax)
    call hist_init(obs%h_mid, bsz, rmax)

    ! Segmental order arrays
    if (allocated(obs%seg_cos_sum)) deallocate(obs%seg_cos_sum)
    if (allocated(obs%seg_cnt))     deallocate(obs%seg_cnt)
    allocate(obs%seg_cos_sum(nb))
    allocate(obs%seg_cnt(nb))
    obs%seg_cos_sum = 0.0_dp
    obs%seg_cnt     = 0_i8

    ! CoM-resolved arrays
    if (allocated(obs%nc))        deallocate(obs%nc)
    if (allocated(obs%rg2_sum))   deallocate(obs%rg2_sum)
    if (allocated(obs%re2_sum))   deallocate(obs%re2_sum)
    if (allocated(obs%a_sum))     deallocate(obs%a_sum)
    if (allocated(obs%b_sum))     deallocate(obs%b_sum)
    if (allocated(obs%c_sum))     deallocate(obs%c_sum)
    if (allocated(obs%acig_sum))  deallocate(obs%acig_sum)
    if (allocated(obs%bcig_sum))  deallocate(obs%bcig_sum)
    if (allocated(obs%ccig_sum))  deallocate(obs%ccig_sum)
    allocate(obs%nc(nbc))
    allocate(obs%rg2_sum(nbc))
    allocate(obs%re2_sum(nbc))
    allocate(obs%a_sum(nbc))
    allocate(obs%b_sum(nbc))
    allocate(obs%c_sum(nbc))
    allocate(obs%acig_sum(nbc))
    allocate(obs%bcig_sum(nbc))
    allocate(obs%ccig_sum(nbc))
    obs%nc        = 0_i8
    obs%rg2_sum   = 0.0_dp
    obs%re2_sum   = 0.0_dp
    obs%a_sum     = 0.0_dp
    obs%b_sum     = 0.0_dp
    obs%c_sum     = 0.0_dp
    obs%acig_sum  = 0.0_dp
    obs%bcig_sum  = 0.0_dp
    obs%ccig_sum  = 0.0_dp

  end subroutine obs_init

  ! ------------------------------------------------------------------
  !> obs_sample — accumulate one snapshot over all chains in sys.
  !>
  !> Ported faithfully from legacy polyfj.f lines 259-394.
  !> Coords in sys are real sigma units; legacy reduced-unit factors
  !> (BSZR, AL^2, etc.) are absorbed into the real-unit equivalents.
  ! ------------------------------------------------------------------
  subroutine obs_sample(obs, sys)
    type(polyfj_obs_t), intent(inout) :: obs
    type(system_t),     intent(in)    :: sys

    integer   :: i, j, ij, ij1, ij2, k
    integer   :: nm1, nm2              ! middle-bead in-chain indices
    real(dp)  :: xi, yi, zi            ! min-imaged bead position
    real(dp)  :: dist                  ! distance from origin
    real(dp)  :: x11, x22, xm, coso   ! segmental order temporaries

    ! CoM temporaries
    real(dp)  :: xcm, ycm, zcm        ! raw CoM (not min-imaged)
    real(dp)  :: x1cm, y1cm, z1cm     ! min-imaged CoM
    real(dp)  :: xcmb                 ! |CoM| (min-imaged)
    real(dp)  :: rjx, rjy, rjz        ! bead deviation from raw CoM
    real(dp)  :: rg2                  ! Rg^2 for this chain
    real(dp)  :: xe, ye, ze           ! end-to-end vector (raw)
    real(dp)  :: re2                  ! Re^2 for this chain

    ! Gyration tensor and eigendecomposition
    real(dp)  :: gyr(3,3), eval(3), evec(3,3)
    real(dp)  :: ea, eb, ec            ! semi-axes a>=b>=c (corrected)

    ! Legacy inertia tensor and eigendecomposition
    real(dp)  :: ai(3,3)               ! inertia tensor (legacy AI)
    real(dp)  :: ai_eval(3), ai_evec(3,3)
    real(dp)  :: eig1, eig2, eig3      ! eigenvalues sorted ascending
    real(dp)  :: cigai, cigbi, cigci   ! legacy semi-axes
    real(dp)  :: etmp

    integer   :: n, nmol
    real(dp)  :: L

    n    = sys%n
    nmol = sys%nmol
    L    = sys%box%L

    ! Middle-bead indices (1-based within chain), matching legacy lines 74-80:
    !   if (MOD(N,2)==0): NM1=N/2,   NM2=N/2+1
    !   if (MOD(N,2)/=0): NM1=NM2=(N-1)/2+1
    if (mod(n, 2) == 0) then
      nm1 = n / 2
      nm2 = nm1 + 1
    else
      nm1 = (n - 1) / 2 + 1
      nm2 = nm1
    end if

    ! ==============================================================
    ! Pass 1: density + segmental order (legacy lines 260-288)
    ! Each bead is visited once in the inner loop over j=1..n.
    ! ==============================================================
    do i = 1, nmol
      do j = 1, n
        ij = (i-1)*n + j

        ! Min-image bead position (legacy: X1I = X1(IJ) - DNINT(X1(IJ))
        ! then DIST = sqrt(X1I^2+...) with reduced coords; in sigma units
        ! the box centre is at L/2 but the min-image convention maps to (-L/2, L/2])
        ! Legacy: coords stored in [0,L) reduced, DNINT(x_red) gives nearest integer,
        ! so x_red - DNINT(x_red) is equivalent to min_image(x - L/2, L) ... actually
        ! legacy reduced coords are in [0,1) and centre is 0 (box spans [-0.5, 0.5]),
        ! so x_red - DNINT(x_red) is the min-image from the box centre.
        ! In our real units, min_image(x, L) already gives the signed displacement
        ! from the box centre at 0 (i.e. for coords stored in [0,L) we need
        ! min_image(x - L/2 offset ... ).
        ! Actually min_image(d, L) = d - L*nint(d/L), which for d in [0,L) gives
        ! d for d in [0,L/2) and d-L for d in [L/2,L). So min_image(x, L) correctly
        ! gives the displacement from nearest periodic image of the origin. That is
        ! the legacy intent.
        xi = min_image(sys%x(ij), L)
        yi = min_image(sys%y(ij), L)
        zi = min_image(sys%z(ij), L)
        dist = sqrt(xi*xi + yi*yi + zi*zi)

        ! Density accumulation (legacy lines 267-271)
        call hist_add(obs%h_all, dist)
        if (j == 1 .or. j == n)     call hist_add(obs%h_end, dist)
        if (j == nm1 .or. j == nm2) call hist_add(obs%h_mid, dist)

        ! Segmental order (legacy lines 274-286): for j >= 2
        ! Uses RAW (non-min-imaged) stored coords, matching legacy which uses
        ! X1(IJ1) and X1(IJ2) directly (the commented-out min-image block
        ! at line 300-305 confirms raw coords are intentional here too).
        if (j >= 2) then
          ij1 = (i-1)*n + j       ! = ij
          ij2 = (i-1)*n + j - 1

          ! Magnitudes of raw position vectors
          x11 = sqrt(sys%x(ij1)**2 + sys%y(ij1)**2 + sys%z(ij1)**2)
          x22 = sqrt(sys%x(ij2)**2 + sys%y(ij2)**2 + sys%z(ij2)**2)

          ! Half-separation distance (legacy: XM = |r2-r1|/2 / AL, then K=1+INT(XM/BSZR))
          ! BSZR = BSZ/AL, so K = 1 + INT((|r2-r1|/2 / AL) / (BSZ/AL))
          !                     = 1 + INT(|r2-r1|/2 / BSZ)
          ! In real sigma units: xm = |r2-r1|/2
          xm = sqrt((sys%x(ij2)-sys%x(ij1))**2 + &
                    (sys%y(ij2)-sys%y(ij1))**2 + &
                    (sys%z(ij2)-sys%z(ij1))**2) * 0.5_dp

          ! Cosine of angle between position vectors r1 and r2
          ! (legacy: COSO = dot(r1,r2)/(|r1|*|r2|) in reduced units)
          ! In real units: COSO is dimensionless — identical formula
          if (x11 > 0.0_dp .and. x22 > 0.0_dp) then
            coso = (sys%x(ij1)*sys%x(ij2) + &
                    sys%y(ij1)*sys%y(ij2) + &
                    sys%z(ij1)*sys%z(ij2)) / (x11 * x22)
          else
            coso = 0.0_dp
          end if

          ! Bin and accumulate (legacy: K=1+INT(XM/BSZR), NR(K)+=1, SR(K)+=COSO*AL*AL)
          ! In real units, xm is already in sigma, bin width is bsz in sigma.
          ! We accumulate plain coso (no AL^2 factor — see unit-conversion notes above).
          k = 1 + int(xm / obs%bsz)
          if (k >= 1 .and. k <= obs%nbin) then
            obs%seg_cnt(k)     = obs%seg_cnt(k) + 1_i8
            obs%seg_cos_sum(k) = obs%seg_cos_sum(k) + coso
          end if
        end if

      end do ! j
    end do ! i (density + segorder pass)

    ! ==============================================================
    ! Pass 2: CoM-resolved Rg / Re / shape (legacy lines 290-394)
    ! Uses RAW (non-min-imaged) stored coords for CoM and gyration tensor.
    ! Only the final CoM position itself is min-imaged for the bin index.
    ! ==============================================================
    do i = 1, nmol

      ! --- Compute raw CoM (legacy lines 292-315) ---
      xcm = 0.0_dp
      ycm = 0.0_dp
      zcm = 0.0_dp
      do j = 1, n
        ij = (i-1)*n + j
        ! Use raw coords, NOT min-imaged (legacy comment says so explicitly
        ! at lines 300-305 where it comments out the min-image version)
        xcm = xcm + sys%x(ij)
        ycm = ycm + sys%y(ij)
        zcm = zcm + sys%z(ij)
      end do
      xcm = xcm / real(n, dp)
      ycm = ycm / real(n, dp)
      zcm = zcm / real(n, dp)

      ! Min-image of CoM for distance-to-origin and bin assignment
      ! (legacy: X1CM = XCM - DNINT(XCM) in reduced units → min_image in sigma)
      x1cm = min_image(xcm, L)
      y1cm = min_image(ycm, L)
      z1cm = min_image(zcm, L)
      xcmb = sqrt(x1cm*x1cm + y1cm*y1cm + z1cm*z1cm)

      ! Bin index (legacy: K = 1 + INT(XCMB/BSZC), BSZC = BSZCM/AL)
      ! XCMB legacy is in reduced units; here xcmb is in sigma, bszc in sigma.
      k = 1 + int(xcmb / obs%bszc)
      if (k < 1 .or. k > obs%nbinc) cycle

      ! --- Rg^2 via sum of squared deviations from raw CoM ---
      ! Legacy uses <x^2> - xcm^2 component-wise (König's theorem).
      ! Equivalently: Rg^2 = (1/N) * sum_j |r_j - CoM|^2
      ! Both use RAW coords; the final value is in sigma^2.
      rg2 = 0.0_dp
      do j = 1, n
        ij  = (i-1)*n + j
        rjx = sys%x(ij) - xcm
        rjy = sys%y(ij) - ycm
        rjz = sys%z(ij) - zcm
        rg2 = rg2 + rjx*rjx + rjy*rjy + rjz*rjz
      end do
      rg2 = rg2 / real(n, dp)  ! sigma^2

      ! Legacy accumulation:
      !   RG(K) += (RG1I+RG2I+RG3I)*AL^2  where RGiI are in reduced^2
      !   In sigma units: rg2 is already in sigma^2, no AL^2 factor.
      obs%rg2_sum(k) = obs%rg2_sum(k) + rg2

      ! --- Re^2: end-to-end squared (raw coords, sigma^2) ---
      ! Legacy (lines 324-331): XE = X1(IN) - X1(I1) in reduced units,
      ! RE(K) += (RE1I+RE2I+RE3I)*AL^2
      ij  = (i-1)*n + 1   ! first bead
      ij2 = (i-1)*n + n   ! last bead
      xe  = sys%x(ij2) - sys%x(ij)
      ye  = sys%y(ij2) - sys%y(ij)
      ze  = sys%z(ij2) - sys%z(ij)
      re2 = xe*xe + ye*ye + ze*ze  ! sigma^2

      obs%re2_sum(k) = obs%re2_sum(k) + re2

      ! --- Gyration tensor + eigendecomposition (shape semi-axes) ---
      ! Build G_{ab} = (1/N) * sum_j (r_j - CoM)_a * (r_j - CoM)_b
      ! using raw coords; symmetric upper triangle then fill lower.
      gyr = 0.0_dp
      do j = 1, n
        ij  = (i-1)*n + j
        rjx = sys%x(ij) - xcm
        rjy = sys%y(ij) - ycm
        rjz = sys%z(ij) - zcm
        gyr(1,1) = gyr(1,1) + rjx*rjx
        gyr(2,2) = gyr(2,2) + rjy*rjy
        gyr(3,3) = gyr(3,3) + rjz*rjz
        gyr(1,2) = gyr(1,2) + rjx*rjy
        gyr(1,3) = gyr(1,3) + rjx*rjz
        gyr(2,3) = gyr(2,3) + rjy*rjz
      end do
      gyr = gyr / real(n, dp)
      gyr(2,1) = gyr(1,2)
      gyr(3,1) = gyr(1,3)
      gyr(3,2) = gyr(2,3)

      ! Diagonalise; eigenvalues returned sorted ascending by jacobi3
      call jacobi3(gyr, eval, evec)

      ! sqrt of eigenvalues, sorted descending (a>=b>=c)
      ! eval(1)<=eval(2)<=eval(3) from jacobi3, so:
      ea = sqrt(max(eval(3), 0.0_dp))   ! largest
      eb = sqrt(max(eval(2), 0.0_dp))   ! middle
      ec = sqrt(max(eval(1), 0.0_dp))   ! smallest

      ! Accumulate corrected semi-axes (driver: mean = a_sum/nc)
      obs%nc(k)    = obs%nc(k) + 1_i8
      obs%a_sum(k) = obs%a_sum(k) + ea
      obs%b_sum(k) = obs%b_sum(k) + eb
      obs%c_sum(k) = obs%c_sum(k) + ec

      ! --- Legacy inertia tensor (polyfj.f lines 333-374) ---
      ! Build I_{ab} = sum_j (delta_ab * |rj|^2 - rja*rjb) with the
      ! sign convention used by legacy (off-diagonal terms negated):
      !   I(1,1) = sum(rjy^2 + rjz^2);  I(1,2) = -sum(rjx*rjy); etc.
      ai = 0.0_dp
      do j = 1, n
        ij  = (i-1)*n + j
        rjx = sys%x(ij) - xcm
        rjy = sys%y(ij) - ycm
        rjz = sys%z(ij) - zcm
        ai(1,1) = ai(1,1) + rjy*rjy + rjz*rjz
        ai(2,2) = ai(2,2) + rjx*rjx + rjz*rjz
        ai(3,3) = ai(3,3) + rjx*rjx + rjy*rjy
        ai(1,2) = ai(1,2) - rjx*rjy
        ai(1,3) = ai(1,3) - rjx*rjz
        ai(2,3) = ai(2,3) - rjy*rjz
      end do
      ai(2,1) = ai(1,2)
      ai(3,1) = ai(1,3)
      ai(3,2) = ai(2,3)

      ! Diagonalise; jacobi3 returns eigenvalues sorted ascending
      call jacobi3(ai, ai_eval, ai_evec)

      ! ai_eval(1) <= ai_eval(2) <= ai_eval(3)
      eig1 = ai_eval(1)
      eig2 = ai_eval(2)
      eig3 = ai_eval(3)

      ! Legacy semi-axis formulae (polyfj.f lines 372-374).
      ! max(..., 0) guards sqrt against negative values that arise for
      ! rod-like chains (legacy did not guard and could produce NaN).
      cigai = sqrt(max(2.5_dp*(eig2 + eig3 - eig1) / real(n, dp), 0.0_dp))
      cigbi = sqrt(max(2.5_dp*(eig1 + eig3 - eig2) / real(n, dp), 0.0_dp))
      cigci = sqrt(max(2.5_dp*(eig2 + eig1 - eig3) / real(n, dp), 0.0_dp))

      ! Accumulate into the same CoM-distance bin k used above.
      ! Driver finalises as sqrt(acig_sum/nc) to match legacy's sqrt-of-mean.
      obs%acig_sum(k) = obs%acig_sum(k) + cigai
      obs%bcig_sum(k) = obs%bcig_sum(k) + cigbi
      obs%ccig_sum(k) = obs%ccig_sum(k) + cigci

    end do ! i (Rg/Re/shape pass)

  end subroutine obs_sample

end module polyfj_observables
