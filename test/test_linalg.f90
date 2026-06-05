! test_linalg.f90 — unit tests for polymc_linalg (symmetric 3×3 Jacobi eigensolver)
!
! Uses test-drive framework (v0.6.0).
!
! NOTE: gfortran-15 allocatable-finalizer bug workaround:
!   allocate(testsuite(N)) + element-wise assignment; no array ctor.
!
! Convention: jacobi3 returns eigenvalues sorted ascending.
! Eigenvector columns are ordered to match the sorted eigenvalues.

module test_linalg_suite
  use testdrive,      only: new_unittest, unittest_type, error_type, check
  use polymc_kinds,   only: dp
  use polymc_linalg,  only: jacobi3
  implicit none
  private
  public :: collect_linalg

contains

  !> Register all tests (gfortran-15-safe: allocate + assign, no array ctor)
  subroutine collect_linalg(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    allocate(testsuite(5))
    testsuite(1) = new_unittest("diagonal",       test_diagonal)
    testsuite(2) = new_unittest("known_symmetric", test_known_symmetric)
    testsuite(3) = new_unittest("eigen_relation",  test_eigen_relation)
    testsuite(4) = new_unittest("reconstruction",  test_reconstruction)
    testsuite(5) = new_unittest("input_unchanged", test_input_unchanged)
  end subroutine collect_linalg

  ! ------------------------------------------------------------------
  ! Helper: build a 3×3 diagonal matrix from three values
  ! ------------------------------------------------------------------
  pure function make_diag(d1, d2, d3) result(m)
    real(dp), intent(in) :: d1, d2, d3
    real(dp) :: m(3,3)
    m = 0.0_dp
    m(1,1) = d1;  m(2,2) = d2;  m(3,3) = d3
  end function make_diag

  ! ------------------------------------------------------------------
  ! Helper: 3×3 identity matrix
  ! ------------------------------------------------------------------
  pure function eye3() result(m)
    real(dp) :: m(3,3)
    m = 0.0_dp
    m(1,1) = 1.0_dp;  m(2,2) = 1.0_dp;  m(3,3) = 1.0_dp
  end function eye3

  ! ------------------------------------------------------------------
  ! diagonal: a = diag(3,5,7)
  !   sorted eval = [3, 5, 7]  within 1e-10
  !   evec columns are unit axis vectors (up to sign)
  ! ------------------------------------------------------------------
  subroutine test_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), eval(3), evec(3,3)
    real(dp), parameter :: tol = 1.0e-10_dp

    a = make_diag(3.0_dp, 5.0_dp, 7.0_dp)
    call jacobi3(a, eval, evec)

    ! eigenvalues sorted ascending
    call check(error, abs(eval(1) - 3.0_dp) < tol, &
               "diagonal: eval(1) /= 3")
    if (allocated(error)) return
    call check(error, abs(eval(2) - 5.0_dp) < tol, &
               "diagonal: eval(2) /= 5")
    if (allocated(error)) return
    call check(error, abs(eval(3) - 7.0_dp) < tol, &
               "diagonal: eval(3) /= 7")
    if (allocated(error)) return

    ! each eigenvector column must be a unit vector
    call check(error, abs(norm2(evec(:,1)) - 1.0_dp) < tol, &
               "diagonal: evec(:,1) not a unit vector")
    if (allocated(error)) return
    call check(error, abs(norm2(evec(:,2)) - 1.0_dp) < tol, &
               "diagonal: evec(:,2) not a unit vector")
    if (allocated(error)) return
    call check(error, abs(norm2(evec(:,3)) - 1.0_dp) < tol, &
               "diagonal: evec(:,3) not a unit vector")
    if (allocated(error)) return

    ! each column must be aligned with one of the axis vectors (|dot| ~ 1)
    call check(error, abs(abs(dot_product(evec(:,1), [1.0_dp,0.0_dp,0.0_dp])) &
               + abs(dot_product(evec(:,1), [0.0_dp,1.0_dp,0.0_dp])) &
               + abs(dot_product(evec(:,1), [0.0_dp,0.0_dp,1.0_dp])) - 1.0_dp) < tol, &
               "diagonal: evec(:,1) not aligned with any axis")
  end subroutine test_diagonal

  ! ------------------------------------------------------------------
  ! known_symmetric: a = [[2,1,0],[1,2,0],[0,0,3]]
  !   analytic eigenvalues: 1, 3, 3  (sorted: 1, 3, 3)
  ! ------------------------------------------------------------------
  subroutine test_known_symmetric(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), eval(3), evec(3,3)
    real(dp), parameter :: tol = 1.0e-9_dp

    a = 0.0_dp
    a(1,1) = 2.0_dp;  a(1,2) = 1.0_dp
    a(2,1) = 1.0_dp;  a(2,2) = 2.0_dp
    a(3,3) = 3.0_dp

    call jacobi3(a, eval, evec)

    call check(error, abs(eval(1) - 1.0_dp) < tol, &
               "known_symmetric: eval(1) /= 1")
    if (allocated(error)) return
    call check(error, abs(eval(2) - 3.0_dp) < tol, &
               "known_symmetric: eval(2) /= 3")
    if (allocated(error)) return
    call check(error, abs(eval(3) - 3.0_dp) < tol, &
               "known_symmetric: eval(3) /= 3")
  end subroutine test_known_symmetric

  ! ------------------------------------------------------------------
  ! eigen_relation: a = [[4,1,2],[1,3,0],[2,0,5]]
  !   verify: A * v_k ≈ λ_k * v_k  within 1e-8  for each k
  !   verify: V^T V ≈ I             within 1e-8
  ! ------------------------------------------------------------------
  subroutine test_eigen_relation(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), eval(3), evec(3,3)
    real(dp) :: av(3), lv(3), diff(3), vtv(3,3), idiff(3,3)
    real(dp), parameter :: tol = 1.0e-8_dp
    integer :: k, i, j

    a = 0.0_dp
    a(1,1) = 4.0_dp;  a(1,2) = 1.0_dp;  a(1,3) = 2.0_dp
    a(2,1) = 1.0_dp;  a(2,2) = 3.0_dp;  a(2,3) = 0.0_dp
    a(3,1) = 2.0_dp;  a(3,2) = 0.0_dp;  a(3,3) = 5.0_dp

    call jacobi3(a, eval, evec)

    ! A * evec(:,k) == eval(k) * evec(:,k)
    do k = 1, 3
      av = matmul(a, evec(:,k))
      lv = eval(k) * evec(:,k)
      diff = av - lv
      do i = 1, 3
        call check(error, abs(diff(i)) < tol, &
                   "eigen_relation: A*v_k /= lambda_k*v_k")
        if (allocated(error)) return
      end do
    end do

    ! V^T V should be identity (orthonormality)
    vtv = matmul(transpose(evec), evec)
    idiff = vtv - eye3()
    do i = 1, 3
      do j = 1, 3
        call check(error, abs(idiff(i,j)) < tol, &
                   "eigen_relation: eigenvectors not orthonormal")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_eigen_relation

  ! ------------------------------------------------------------------
  ! reconstruction: V * diag(eval) * V^T ≈ a  within 1e-8
  ! ------------------------------------------------------------------
  subroutine test_reconstruction(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), eval(3), evec(3,3), d(3,3), recon(3,3), diff(3,3)
    real(dp), parameter :: tol = 1.0e-8_dp
    integer :: k, i, j

    a = 0.0_dp
    a(1,1) = 4.0_dp;  a(1,2) = 1.0_dp;  a(1,3) = 2.0_dp
    a(2,1) = 1.0_dp;  a(2,2) = 3.0_dp;  a(2,3) = 0.0_dp
    a(3,1) = 2.0_dp;  a(3,2) = 0.0_dp;  a(3,3) = 5.0_dp

    call jacobi3(a, eval, evec)

    ! build diagonal matrix from eval
    d = 0.0_dp
    do k = 1, 3
      d(k,k) = eval(k)
    end do

    recon = matmul(matmul(evec, d), transpose(evec))
    diff  = recon - a
    do i = 1, 3
      do j = 1, 3
        call check(error, abs(diff(i,j)) < tol, &
                   "reconstruction: V*D*V^T /= A")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_reconstruction

  ! ------------------------------------------------------------------
  ! input_unchanged: the caller's matrix is identical before and after
  ! ------------------------------------------------------------------
  subroutine test_input_unchanged(error)
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), a_before(3,3), eval(3), evec(3,3)
    real(dp), parameter :: tol = 0.0_dp   ! exact bitwise equality
    integer :: i, j

    a = 0.0_dp
    a(1,1) = 4.0_dp;  a(1,2) = 1.0_dp;  a(1,3) = 2.0_dp
    a(2,1) = 1.0_dp;  a(2,2) = 3.0_dp;  a(2,3) = 0.0_dp
    a(3,1) = 2.0_dp;  a(3,2) = 0.0_dp;  a(3,3) = 5.0_dp

    a_before = a
    call jacobi3(a, eval, evec)

    do i = 1, 3
      do j = 1, 3
        call check(error, a(i,j) == a_before(i,j), &
                   "input_unchanged: a was modified by jacobi3")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_input_unchanged

end module test_linalg_suite

! ==================================================================
!  Driver program — fpm auto-tests compiles this as an executable
! ==================================================================
program test_linalg
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive,          only: run_testsuite, new_testsuite, testsuite_type
  use test_linalg_suite,  only: collect_linalg
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  ! Single-element array constructor is safe from gfortran-15 finalizer bug.
  testsuites = [new_testsuite("linalg", collect_linalg)]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if
end program test_linalg
