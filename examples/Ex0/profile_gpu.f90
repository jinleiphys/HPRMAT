!------------------------------------------------------------------------------
! profile_gpu.f90
! GPU-only driver for stage profiling of the Type-4 (cuSOLVER mixed-precision)
! solver. Reuses the same dense diagonally-dominant synthetic matrix as
! benchmark_large.f90, does one untimed warmup solve (pays GPU context creation
! and host-matrix page-locking), then nrep timed steady-state Type-4 solves.
! Run under nsys to decompose the steady-state solve into H2D / getrf / getrs /
! convert / D2H and reconcile the per-stage sum against this wall time.
!   usage: profile_gpu <nch> <nlag> <nrep>
!------------------------------------------------------------------------------
program profile_gpu
  use precision
  use rmat_solvers
  implicit none
  integer :: nch, nlag, ntotal, i, j, irep, nrep
  real(dp) :: normfac, elapsed
  integer(8) :: t1, t2, clock_rate
  complex(dp), allocatable :: cmat(:,:), B_vector(:), Rmat(:,:)
  character(len=32) :: arg

  nch = 256; nlag = 100; nrep = 3
  if (command_argument_count() >= 2) then
    call get_command_argument(1, arg); read(arg,*) nch
    call get_command_argument(2, arg); read(arg,*) nlag
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg); read(arg,*) nrep
  end if
  ntotal = nch * nlag; normfac = 1.0_dp

  allocate(cmat(ntotal,ntotal), B_vector(nlag), Rmat(nch,nch))
  call random_seed()
  do j = 1, ntotal
    do i = 1, ntotal
      cmat(i,j) = cmplx(rand(), rand(), dp)
    end do
  end do
  do i = 1, ntotal
    cmat(i,i) = cmat(i,i) + cmplx(ntotal * 2.0_dp, 0.0_dp, dp)
  end do
  do i = 1, nlag
    B_vector(i) = cmplx(1.0_dp / sqrt(real(i, dp)), 0.0_dp, dp)
  end do

  call system_clock(count_rate=clock_rate)
  write(*,'(A,I0,A,I0,A,I0)') 'N=', ntotal, '  nch=', nch, '  nrep=', nrep

#ifdef GPU_ENABLED
  ! Untimed warmup: GPU context, device buffers, host-matrix pinning (cached).
  call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=4)
  call system_clock(t1)
  do irep = 1, nrep
    call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=4)
  end do
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(clock_rate, dp) / nrep
  write(*,'(A,I0,A,F9.4)') 'RESULT N=', ntotal, '  Type4_wall_s=', elapsed

  ! Hybrid mode (max_refine=2): GPU FP32 factorization + host FP64 refinement.
  call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=4, max_refine=2)
  call system_clock(t1)
  do irep = 1, nrep
    call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=4, max_refine=2)
  end do
  call system_clock(t2)
  elapsed = real(t2 - t1, dp) / real(clock_rate, dp) / nrep
  write(*,'(A,I0,A,F9.4)') 'RESULT N=', ntotal, '  Hybrid_wall_s=', elapsed
#else
  write(*,*) 'built without GPU_ENABLED'
#endif

  deallocate(cmat, B_vector, Rmat)
end program profile_gpu
