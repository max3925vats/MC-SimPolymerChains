! polymc_binning.f90 — radial histogram accumulator with shell-volume density normalization.
!
! hist_t accumulates counts in equal-width radial bins.  hist_density converts
! counts to a number density by dividing by the spherical shell volume of each
! bin and the number of samples, mirroring the legacy polyfj normalization.

module polymc_binning
  use polymc_kinds, only: dp, i8
  implicit none
  private
  public :: hist_t, hist_init, hist_add, hist_density

  ! pi is used in shell-volume calculation
  real(dp), parameter :: pi = 3.141592653589793_dp

  type :: hist_t
    real(dp)              :: bsz  = 0.0_dp   ! bin width
    integer               :: nbin = 0         ! number of bins
    integer(i8), allocatable :: count(:)      ! per-bin counts (1-indexed)
  end type

contains

  ! ------------------------------------------------------------------
  ! hist_init — set bin width, compute nbin = int(rmax/bsz),
  !             allocate count(nbin) and zero it.
  ! ------------------------------------------------------------------
  subroutine hist_init(h, bsz, rmax)
    type(hist_t), intent(inout) :: h
    real(dp),     intent(in)    :: bsz
    real(dp),     intent(in)    :: rmax

    h%bsz  = bsz
    h%nbin = int(rmax / bsz)

    if (allocated(h%count)) deallocate(h%count)
    allocate(h%count(h%nbin))
    h%count = 0_i8
  end subroutine hist_init

  ! ------------------------------------------------------------------
  ! hist_add — bin the radius r.
  !   k = int(r / bsz) + 1
  !   Increment count(k) if 1 <= k <= nbin; otherwise ignore.
  ! ------------------------------------------------------------------
  subroutine hist_add(h, r)
    type(hist_t), intent(inout) :: h
    real(dp),     intent(in)    :: r

    integer :: k

    ! Guard against negative r: int(-x/bsz) truncates toward zero giving k>=1,
    ! so we must explicitly reject r < 0 before binning.
    if (r < 0.0_dp) return

    k = int(r / h%bsz) + 1
    if (k >= 1 .and. k <= h%nbin) then
      h%count(k) = h%count(k) + 1_i8
    end if
  end subroutine hist_add

  ! ------------------------------------------------------------------
  ! hist_density — compute bin-centre positions (dist) and the
  !   number density (dens) normalized by shell volume and sample count.
  !
  !   For bin k (1-indexed):
  !     rl = (k-1)*bsz          (inner shell radius)
  !     ru =  k   *bsz          (outer shell radius)
  !     dist(k) = rl + 0.5*bsz  (bin centre)
  !     v = (4/3)*pi*(ru^3 - rl^3)   (shell volume)
  !     dens(k) = count(k) / (v * nsamp)
  !
  !   Caller must pass dist and dens sized to at least nbin.
  !   nsamp <= 0 gives dens = 0 everywhere.
  ! ------------------------------------------------------------------
  subroutine hist_density(h, nsamp, dist, dens)
    type(hist_t), intent(in)  :: h
    integer,      intent(in)  :: nsamp
    real(dp),     intent(out) :: dist(:)
    real(dp),     intent(out) :: dens(:)

    integer  :: k
    real(dp) :: rl, ru, v

    do k = 1, h%nbin
      rl = real(k-1, dp) * h%bsz
      ru = real(k,   dp) * h%bsz
      dist(k) = rl + 0.5_dp * h%bsz

      if (nsamp > 0) then
        v = (4.0_dp / 3.0_dp) * pi * (ru**3 - rl**3)
        dens(k) = real(h%count(k), dp) / (v * real(nsamp, dp))
      else
        dens(k) = 0.0_dp
      end if
    end do
  end subroutine hist_density

end module polymc_binning
