module polymc_rng
  use polymc_kinds, only: dp, i8
  implicit none
  private
  public :: rng_t, rng_seed, rng_uniform, rng_unit_vector

  type :: rng_t
    integer(i8) :: s(4) = 0_i8
  end type

contains

  ! xoshiro256** rotation helper
  pure function rotl(x, k) result(r)
    integer(i8), intent(in) :: x
    integer,     intent(in) :: k
    integer(i8)             :: r
    r = ior(ishft(x, k), ishft(x, k-64))
  end function

  ! Seed the RNG state via splitmix64 from a single 64-bit integer
  subroutine rng_seed(rng, seed)
    type(rng_t),    intent(out) :: rng
    integer(i8),    intent(in)  :: seed
    integer(i8) :: z
    integer     :: i
    z = seed
    do i = 1, 4                          ! splitmix64 fill
      z = z + (-7046029254386353131_i8)
      z = ieor(z, ishft(z,-30)) * (-4658895280553007687_i8)
      z = ieor(z, ishft(z,-27)) * (-7723592293110705685_i8)
      rng%s(i) = ieor(z, ishft(z,-31))
    end do
  end subroutine

  ! Raw xoshiro256** output — advances state and returns 64-bit integer
  function next(rng) result(r)
    type(rng_t), intent(inout) :: rng
    integer(i8) :: r, t
    r = rotl(rng%s(2) * 5_i8, 7) * 9_i8
    t = ishft(rng%s(2), 17)
    rng%s(3) = ieor(rng%s(3), rng%s(1))
    rng%s(4) = ieor(rng%s(4), rng%s(2))
    rng%s(2) = ieor(rng%s(2), rng%s(3))
    rng%s(1) = ieor(rng%s(1), rng%s(4))
    rng%s(3) = ieor(rng%s(3), t)
    rng%s(4) = rotl(rng%s(4), 45)
  end function

  ! Uniform [0,1) draw using top 53 bits of the raw output
  function rng_uniform(rng) result(u)
    type(rng_t), intent(inout) :: rng
    real(dp) :: u
    integer(i8) :: x
    x = ishft(next(rng), -11)            ! top 53 bits
    u = real(x, dp) * (1.0_dp / 9007199254740992.0_dp)   ! / 2^53
  end function

  ! Uniform unit vector on S^2 via Marsaglia (1972) rejection method
  subroutine rng_unit_vector(rng, vx, vy, vz)
    type(rng_t), intent(inout) :: rng
    real(dp),    intent(out)   :: vx, vy, vz
    real(dp) :: b1, b2, bsq, bh
    do
      b1 = 1.0_dp - 2.0_dp*rng_uniform(rng)
      b2 = 1.0_dp - 2.0_dp*rng_uniform(rng)
      bsq = b1*b1 + b2*b2
      if (bsq < 1.0_dp) exit
    end do
    bh = sqrt(1.0_dp - bsq)
    vx = 2.0_dp*b1*bh; vy = 2.0_dp*b2*bh; vz = 1.0_dp - 2.0_dp*bsq
  end subroutine

end module polymc_rng
