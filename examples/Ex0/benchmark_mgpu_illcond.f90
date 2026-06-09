!------------------------------------------------------------------------------
! benchmark_mgpu_illcond.f90
! Ill-conditioning test of the mixed-precision multi-GPU solver.
!
! Builds a dense matrix with an EXACTLY KNOWN condition number kappa via a
! Householder similarity A = H D H, where H = I - 2 v v^H / (v^H v) is unitary
! (Hermitian) and D = diag(d_i) has geometrically spaced positive entries with
! max/min = kappa. Since H is unitary, the singular values of A are exactly the
! d_i, so cond_2(A) = kappa is exact (no LAPACK / ZGECON needed). The exact
! solution is x = 1, giving the right-hand side b = A 1 in closed form. The test
! reports the forward error max|x_solved - 1| of the FP32-factor + FP64-refinement
! multi-GPU solver as a function of kappa and the number of refinement steps,
! to characterize where "full double precision" holds and where it degrades.
!
! Usage:  benchmark_mgpu_illcond  N  nch  kappa  max_refine
!------------------------------------------------------------------------------
program benchmark_mgpu_illcond
  use precision
#ifdef GPU_ENABLED
  use gpu_solver_interface
#endif
  implicit none

  integer :: n, nch, max_refine, i, j, ich, info
  real(dp) :: kappa, tol, c, alpha, errfwd, re, im
  complex(dp), allocatable :: A(:,:), X(:,:), v(:), bvec(:)
  real(dp), allocatable :: d(:)
  complex(dp) :: S1, S2, twoc
  character(len=32) :: arg

  n = 8000; nch = 64; kappa = 1.0d6; max_refine = 3; tol = 1.0d-13
  if (command_argument_count() >= 4) then
    call get_command_argument(1, arg); read(arg,*) n
    call get_command_argument(2, arg); read(arg,*) nch
    call get_command_argument(3, arg); read(arg,*) kappa
    call get_command_argument(4, arg); read(arg,*) max_refine
  end if

  allocate(A(n,n), X(n,nch), v(n), bvec(n), d(n))

  ! Singular values geometrically spaced, balanced around 1: d_i in [kappa^-1/2, kappa^+1/2],
  ! so max/min = kappa exactly.
  do i = 1, n
    d(i) = kappa**( real(i-1,dp)/real(n-1,dp) - 0.5_dp )
  end do

  ! Random complex Householder vector v.
  call random_seed()
  do i = 1, n
    call random_number(re); call random_number(im)
    v(i) = cmplx(re - 0.5_dp, im - 0.5_dp, dp)
  end do
  c = sum(abs(v)**2)                       ! v^H v  (real)
  alpha = 0.0_dp
  do i = 1, n
    alpha = alpha + abs(v(i))**2 * d(i)    ! v^H D v  (real, D real)
  end do
  twoc = cmplx(2.0_dp/c, 0.0_dp, dp)

  ! A_ij = d_i delta_ij - (2/c) v_i conj(v_j) (d_i + d_j - 2 alpha/c)
  do j = 1, n
    do i = 1, n
      A(i,j) = - twoc * v(i)*conjg(v(j)) * cmplx(d(i)+d(j) - 2.0_dp*alpha/c, 0.0_dp, dp)
    end do
    A(j,j) = A(j,j) + cmplx(d(j), 0.0_dp, dp)
  end do

  ! Exact solution x = 1  ->  b_i = sum_j A_ij
  !   = d_i - (2/c) v_i [ (d_i - 2 alpha/c) S1 + S2 ],  S1 = sum_j conj(v_j), S2 = sum_j conj(v_j) d_j
  S1 = sum(conjg(v))
  S2 = (0.0_dp, 0.0_dp)
  do j = 1, n
    S2 = S2 + conjg(v(j)) * cmplx(d(j), 0.0_dp, dp)
  end do
  do i = 1, n
    bvec(i) = cmplx(d(i),0.0_dp,dp) - twoc * v(i) * &
              ( cmplx(d(i) - 2.0_dp*alpha/c, 0.0_dp, dp)*S1 + S2 )
  end do
  do ich = 1, nch
    X(:,ich) = bvec
  end do

#ifdef GPU_ENABLED
  call gpu_solve_multi_mixed(A, X, n, nch, max_refine, tol, info)
#else
  block
    integer, allocatable :: ipiv(:)
    complex(dp), allocatable :: Acpy(:,:)
    allocate(ipiv(n), Acpy(n,n)); Acpy = A
    call ZGESV(n, nch, Acpy, n, ipiv, X, n, info)
    deallocate(ipiv, Acpy)
  end block
#endif

  errfwd = maxval(abs(X - (1.0_dp, 0.0_dp)))
  write(*,'(A,I7,A,ES9.1,A,I2,A,I3,A,ES12.4)') &
       'N=', n, '  kappa=', kappa, '  refine=', max_refine, &
       '  info=', info, '  max|x-1|=', errfwd

  deallocate(A, X, v, bvec, d)
end program benchmark_mgpu_illcond
