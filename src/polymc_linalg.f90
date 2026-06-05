! polymc_linalg.f90 — linear algebra utilities for polymc
!
! Provides a non-destructive symmetric 3×3 eigensolver using the
! classical cyclic Jacobi rotation method (ported from the legacy
! JACOBI subroutine in polyfj/polyfj.f, specialised to fixed 3×3).
!
! Convention: eigenvalues are returned sorted ascending; eigenvector
! columns are reordered to match (evec(:,k) belongs to eval(k)).
!
! Reference: Numerical Recipes §11.1 (Jacobi transformations of a
! symmetric matrix), also used verbatim in the legacy polyfj code.

module polymc_linalg
  use polymc_kinds, only: dp
  implicit none
  private
  public :: jacobi3

contains

  ! ----------------------------------------------------------------
  !> jacobi3 — symmetric 3×3 eigensolver via cyclic Jacobi rotations
  !>
  !> @param a     [in]  3×3 symmetric matrix (NOT modified)
  !> @param eval  [out] eigenvalues, sorted ascending
  !> @param evec  [out] eigenvector columns: evec(:,k) <-> eval(k)
  ! ----------------------------------------------------------------
  subroutine jacobi3(a, eval, evec)
    real(dp), intent(in)  :: a(3,3)
    real(dp), intent(out) :: eval(3)
    real(dp), intent(out) :: evec(3,3)

    ! Internal working copy; we never touch the caller's a.
    real(dp) :: aw(3,3)

    ! Workspace arrays matching the legacy JACOBI signature
    real(dp) :: b(3), z(3)
    real(dp) :: sm, tresh, g, h, t, theta, c, s, tau
    integer  :: i, j, ip, iq
    logical  :: converged

    ! --- initialise from the input (non-destructive) ---
    aw = a

    ! eigenvector matrix starts as identity
    evec = 0.0_dp
    evec(1,1) = 1.0_dp;  evec(2,2) = 1.0_dp;  evec(3,3) = 1.0_dp

    ! b and eval hold the diagonal; z accumulates off-diagonal drift
    do ip = 1, 3
      b(ip)    = aw(ip,ip)
      eval(ip) = b(ip)
      z(ip)    = 0.0_dp
    end do

    converged = .false.

    ! --- up to 50 sweeps (convergence is guaranteed long before that) ---
    do i = 1, 50

      ! sum of |off-diagonal| elements (upper triangle)
      sm = 0.0_dp
      do ip = 1, 2
        do iq = ip+1, 3
          sm = sm + abs(aw(ip,iq))
        end do
      end do
      if (sm == 0.0_dp) then
        converged = .true.
        exit
      end if

      ! threshold: coarse for first 3 sweeps, zero thereafter
      if (i < 4) then
        tresh = 0.2_dp * sm / 9.0_dp   ! 9 = 3**2
      else
        tresh = 0.0_dp
      end if

      ! loop over all off-diagonal pairs (upper triangle only)
      do ip = 1, 2
        do iq = ip+1, 3

          g = 100.0_dp * abs(aw(ip,iq))

          ! skip this element if it is already negligible
          if ( (i > 4) .and. &
               (abs(eval(ip)) + g == abs(eval(ip))) .and. &
               (abs(eval(iq)) + g == abs(eval(iq))) ) then
            aw(ip,iq) = 0.0_dp

          else if (abs(aw(ip,iq)) > tresh) then

            h = eval(iq) - eval(ip)

            ! compute rotation angle t = tan(angle); choose sign to
            ! minimise |rotation angle|
            if (abs(h) + g == abs(h)) then
              t = aw(ip,iq) / h
            else
              theta = 0.5_dp * h / aw(ip,iq)
              t = 1.0_dp / (abs(theta) + sqrt(1.0_dp + theta**2))
              if (theta < 0.0_dp) t = -t
            end if

            c   = 1.0_dp / sqrt(1.0_dp + t**2)
            s   = t * c
            tau = s / (1.0_dp + c)
            h   = t * aw(ip,iq)

            z(ip) = z(ip) - h
            z(iq) = z(iq) + h
            eval(ip) = eval(ip) - h
            eval(iq) = eval(iq) + h
            aw(ip,iq) = 0.0_dp

            ! update rows 1 .. ip-1
            do j = 1, ip-1
              g = aw(j,ip)
              h = aw(j,iq)
              aw(j,ip) = g - s*(h + g*tau)
              aw(j,iq) = h + s*(g - h*tau)
            end do

            ! update rows ip+1 .. iq-1
            do j = ip+1, iq-1
              g = aw(ip,j)
              h = aw(j,iq)
              aw(ip,j) = g - s*(h + g*tau)
              aw(j,iq) = h + s*(g - h*tau)
            end do

            ! update rows iq+1 .. 3
            do j = iq+1, 3
              g = aw(ip,j)
              h = aw(iq,j)
              aw(ip,j) = g - s*(h + g*tau)
              aw(iq,j) = h + s*(g - h*tau)
            end do

            ! accumulate the rotation into the eigenvector matrix
            do j = 1, 3
              g = evec(j,ip)
              h = evec(j,iq)
              evec(j,ip) = g - s*(h + g*tau)
              evec(j,iq) = h + s*(g - h*tau)
            end do

          end if

        end do   ! iq
      end do     ! ip

      ! update b and reset z
      do ip = 1, 3
        b(ip)    = b(ip) + z(ip)
        eval(ip) = b(ip)
        z(ip)    = 0.0_dp
      end do

    end do   ! sweep i

    ! Guard: symmetric 3×3 always converges in practice, but a silent
    ! non-converged return would be a hard-to-debug failure mode.
    if (.not. converged) error stop "jacobi3: failed to converge"

    ! --- sort eigenvalues ascending; reorder eigenvector columns to match ---
    call insertion_sort3(eval, evec)

  end subroutine jacobi3

  ! ----------------------------------------------------------------
  !> insertion_sort3 — sort 3-element array d ascending, applying
  !> the same permutation to the columns of v (3×3).
  !>
  !> Insertion sort is O(n^2) but trivially fast for n = 3.
  ! ----------------------------------------------------------------
  subroutine insertion_sort3(d, v)
    real(dp), intent(inout) :: d(3)
    real(dp), intent(inout) :: v(3,3)

    real(dp) :: key_d, key_v(3)
    integer  :: i, j

    do i = 2, 3
      key_d = d(i)
      key_v = v(:,i)
      j = i - 1
      do while (j >= 1 .and. d(j) > key_d)
        d(j+1)   = d(j)
        v(:,j+1) = v(:,j)
        j = j - 1
      end do
      d(j+1)   = key_d
      v(:,j+1) = key_v
    end do

  end subroutine insertion_sort3

end module polymc_linalg
