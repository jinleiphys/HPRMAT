!------------------------------------------------------------------------------
! benchmark_mgpu.f90
! Large-N multi-GPU mixed-precision (cusolverMg FP32 factorization + FP64 host
! refinement) capacity/scaling benchmark.
!
! Solves a dense complex system A X = B of dimension N with nch right-hand sides
! using the cusolverMg backend (solver_type 6), distributing A across the GPUs
! made visible by CUDA_VISIBLE_DEVICES. Accuracy is checked against a
! manufactured exact solution (X_exact = 1, so b is the row sum of A): the
! forward error max|X - 1| needs no reference solve or explicit inverse, so
! arbitrarily large N can be validated.
!
! Usage:  benchmark_mgpu  N  nch  [nrep]
!   N    : matrix dimension
!   nch  : number of right-hand sides (channels)
!   nrep : timed repetitions (default 1; one untimed warmup always precedes)
!------------------------------------------------------------------------------
program benchmark_mgpu
  use precision
#ifdef GPU_ENABLED
  use gpu_solver_interface
#endif
  implicit none

  integer :: n, nch, nrep, i, j, ich, irep, info, max_refine
  integer(8) :: t1, t2, rate
  complex(dp), allocatable :: A(:,:), X(:,:), B0(:,:), rowsum(:)
  real(dp), allocatable :: rebuf(:), imbuf(:)
  character(len=32) :: arg
  real(dp) :: elapsed, resid, tol

  n = 12800
  nch = 64
  nrep = 1
  max_refine = 2
  tol = 1.0d-12
  if (command_argument_count() >= 2) then
    call get_command_argument(1, arg); read(arg,*) n
    call get_command_argument(2, arg); read(arg,*) nch
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg); read(arg,*) nrep
  end if
  if (command_argument_count() >= 4) then
    call get_command_argument(4, arg); read(arg,*) max_refine
  end if

  write(*,'(A)') '=========================================================='
  write(*,'(A)') '  Multi-GPU (cusolverMg, FP64) large-N benchmark'
  write(*,'(A)') '=========================================================='
  write(*,'(A,I0,A,I0,A,I0,A,I0)') '  N = ', n, ',  nch = ', nch, ',  nrep = ', nrep, &
       ',  max_refine = ', max_refine

  allocate(A(n,n), X(n,nch), B0(n,nch), rowsum(n), rebuf(n), imbuf(n))

  ! Diagonally-dominant random complex matrix (well-conditioned). Column-buffered
  ! RANDOM_NUMBER fill keeps initialization cheap even at N ~ 5e4.
  call random_seed()
  do j = 1, n
    call random_number(rebuf)
    call random_number(imbuf)
    do i = 1, n
      A(i,j) = cmplx(rebuf(i), imbuf(i), dp)
    end do
  end do
  do i = 1, n
    A(i,i) = A(i,i) + cmplx(2.0_dp*real(n,dp), 0.0_dp, dp)
  end do
  deallocate(rebuf, imbuf)

  ! Manufactured solution for a BLAS-free accuracy check: take the exact solution
  ! X_exact = 1 for every entry, so the right-hand side is the row sum of A,
  ! b_i = sum_j A_ij. After the solve, the error is simply max|X - 1|. Row sums are
  ! formed in one pass (intrinsic SUM), so this scales to N ~ 5e4 with no GEMM.
  do i = 1, n
    rowsum(i) = sum(A(i,:))
  end do
  do ich = 1, nch
    X(:,ich) = rowsum
  end do
  B0 = X

#ifdef GPU_ENABLED
  ! Untimed warmup: pays cusolverMg handle/grid creation and first device alloc.
  call gpu_solve_multi_mixed(A, X, n, nch, max_refine, tol, info)
  if (info /= 0) then
    write(*,'(A,I0)') '  WARMUP FAILED, info = ', info
  end if

  call system_clock(count_rate=rate)
  call system_clock(t1)
  do irep = 1, nrep
    X = B0
    call gpu_solve_multi_mixed(A, X, n, nch, max_refine, tol, info)
  end do
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(rate, dp) / real(nrep, dp)
#else
  call system_clock(count_rate=rate)
  call system_clock(t1)
  block
    integer, allocatable :: ipiv(:)
    complex(dp), allocatable :: Acpy(:,:)
    allocate(ipiv(n), Acpy(n,n))
    Acpy = A
    call ZGESV(n, nch, Acpy, n, ipiv, X, n, info)
    deallocate(ipiv, Acpy)
  end block
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(rate, dp)
#endif

  ! Accuracy vs the manufactured solution X_exact = 1.
  resid = maxval(abs(X - (1.0_dp, 0.0_dp)))

  write(*,'(A)') '----------------------------------------------------------'
  write(*,'(A,I0)')      '  info            = ', info
  write(*,'(A,F12.4,A)') '  solve time      = ', elapsed, ' s'
  write(*,'(A,E12.4)')   '  max|X - X_exact|= ', resid
  write(*,'(A)') '=========================================================='

  deallocate(A, X, B0, rowsum)
end program benchmark_mgpu
