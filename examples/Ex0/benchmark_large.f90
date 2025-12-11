!------------------------------------------------------------------------------
! benchmark_large.f90
! Benchmark large matrix: Pierre vs HPRMAT Type 1, 2, 3
!------------------------------------------------------------------------------
program benchmark_large
  use precision
  use rmat_solvers
  implicit none

  integer :: nch, nlag, ntotal, i, j, k, irep, ir, nrep
  real(dp) :: normfac, elapsed, t_pierre
  integer(8) :: t1, t2, clock_rate
  complex(dp), allocatable :: cmat(:,:), cmat_copy(:,:), B_vector(:), Rmat(:,:)
  complex(dp), allocatable :: Rmat_ref(:,:), X_vector(:,:), X_result(:,:), WORK(:)
  integer, allocatable :: IPIV(:)
  real(dp) :: max_diff
  integer :: LWORK
  integer :: INFO
  character(len=32) :: arg

  ! Default: 12800 = 128 channels x 100 basis
  nch = 128
  nlag = 100
  nrep = 1

  ! Read command line arguments
  if (command_argument_count() >= 2) then
    call get_command_argument(1, arg)
    read(arg, *) nch
    call get_command_argument(2, arg)
    read(arg, *) nlag
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) nrep
  end if

  ntotal = nch * nlag
  normfac = 1.0_dp

  write(*,*) '========================================================'
  write(*,*) '  Large Matrix Benchmark: Pierre vs HPRMAT'
  write(*,*) '========================================================'
  write(*,'(A,I0,A,I0)') '  Matrix size: ', ntotal, ' x ', ntotal
  write(*,'(A,I0,A,I0)') '  nch = ', nch, ', nlag = ', nlag
  write(*,'(A,I0)') '  Repetitions: ', nrep
  write(*,*) ''

  LWORK = ntotal * 64  ! Workspace for ZGETRI
  allocate(cmat(ntotal, ntotal))
  allocate(cmat_copy(ntotal, ntotal))
  allocate(B_vector(nlag))
  allocate(X_vector(ntotal, nch))
  allocate(X_result(ntotal, nch))
  allocate(Rmat(nch, nch))
  allocate(Rmat_ref(nch, nch))
  allocate(IPIV(ntotal))
  allocate(WORK(LWORK))

  ! Initialize test matrix (diagonally dominant for stability)
  write(*,*) 'Initializing matrix...'
  call random_seed()
  do j = 1, ntotal
    do i = 1, ntotal
      cmat(i, j) = cmplx(rand(), rand(), dp)
    end do
  end do
  do i = 1, ntotal
    cmat(i, i) = cmat(i, i) + cmplx(ntotal * 2.0_dp, 0.0_dp, dp)
  end do

  ! Initialize B vector
  do i = 1, nlag
    B_vector(i) = cmplx(1.0_dp / sqrt(real(i, dp)), 0.0_dp, dp)
  end do

  write(*,*) ''
  write(*,*) '| Method | Time (s) | Speedup | Max Diff |'
  write(*,*) '|--------|----------|---------|----------|'

  ! Get clock rate for wall time measurement
  call system_clock(count_rate=clock_rate)

  !---------------------------------------------------------------------------
  ! Pierre's method: Matrix inversion (ZGETRF + ZGETRI)
  ! Then multiply: X = C^{-1} * B_columns, then extract R
  !---------------------------------------------------------------------------
  cmat_copy = cmat
  call system_clock(t1)
  do irep = 1, nrep
    cmat_copy = cmat
    call ZGETRF(ntotal, ntotal, cmat_copy, ntotal, IPIV, INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: ZGETRF failed'
      stop
    end if
    call ZGETRI(ntotal, cmat_copy, ntotal, IPIV, WORK, LWORK, INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: ZGETRI failed'
      stop
    end if
    ! Setup RHS matrix (same as HPRMAT)
    X_vector = (0.0_dp, 0.0_dp)
    do i = 1, nch
      do ir = 1, nlag
        X_vector((i-1)*nlag + ir, i) = B_vector(ir)
      end do
    end do
    ! X_result = C^{-1} * B
    call ZGEMM('N', 'N', ntotal, nch, ntotal, (1.0_dp, 0.0_dp), cmat_copy, ntotal, &
               X_vector, ntotal, (0.0_dp, 0.0_dp), X_result, ntotal)
    ! Extract R-matrix: R_ij = sum_r B(r) * X_{(i-1)*nlag+r, j}
    Rmat_ref = (0.0_dp, 0.0_dp)
    do j = 1, nch
      do i = 1, nch
        do ir = 1, nlag
          Rmat_ref(i, j) = Rmat_ref(i, j) + B_vector(ir) * X_result(ir + (i-1)*nlag, j)
        end do
      end do
    end do
  end do
  call system_clock(t2)
  t_pierre = real(t2 - t1, dp) / real(clock_rate, dp) / nrep
  write(*,'(A,F8.3,A)') '| Pierre | ', t_pierre, ' | 1.0x    | (ref)    |'

  !---------------------------------------------------------------------------
  ! Type 1: Dense LAPACK ZGESV
  !---------------------------------------------------------------------------
  call system_clock(t1)
  do irep = 1, nrep
    call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=1)
  end do
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(clock_rate, dp) / nrep
  max_diff = maxval(abs(Rmat - Rmat_ref))
  write(*,'(A,F8.3,A,F5.1,A,E8.1,A)') '| Type 1 | ', elapsed, ' | ', t_pierre/elapsed, 'x   | ', max_diff, ' |'

  !---------------------------------------------------------------------------
  ! Type 2: Mixed Precision
  !---------------------------------------------------------------------------
  call system_clock(t1)
  do irep = 1, nrep
    call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=2)
  end do
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(clock_rate, dp) / nrep
  max_diff = maxval(abs(Rmat - Rmat_ref))
  write(*,'(A,F8.3,A,F5.1,A,E8.1,A)') '| Type 2 | ', elapsed, ' | ', t_pierre/elapsed, 'x   | ', max_diff, ' |'

  !---------------------------------------------------------------------------
  ! Type 3: Woodbury
  !---------------------------------------------------------------------------
  call system_clock(t1)
  do irep = 1, nrep
    call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=3)
  end do
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(clock_rate, dp) / nrep
  max_diff = maxval(abs(Rmat - Rmat_ref))
  write(*,'(A,F8.3,A,F5.1,A,E8.1,A)') '| Type 3 | ', elapsed, ' | ', t_pierre/elapsed, 'x   | ', max_diff, ' |'

#ifdef GPU_ENABLED
  !---------------------------------------------------------------------------
  ! Type 4: GPU cuSOLVER
  !---------------------------------------------------------------------------
  call system_clock(t1)
  do irep = 1, nrep
    call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=4)
  end do
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(clock_rate, dp) / nrep
  max_diff = maxval(abs(Rmat - Rmat_ref))
  write(*,'(A,F8.3,A,F5.1,A,E8.1,A)') '| Type 4 | ', elapsed, ' | ', t_pierre/elapsed, 'x   | ', max_diff, ' |'
#endif

  write(*,*) ''
  write(*,*) 'Benchmark completed!'

  deallocate(cmat, cmat_copy, B_vector, X_vector, X_result, Rmat, Rmat_ref, IPIV, WORK)

end program benchmark_large
