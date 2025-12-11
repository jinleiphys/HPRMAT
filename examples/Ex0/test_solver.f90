!------------------------------------------------------------------------------
! test_solver.f90
! Simple test program to benchmark different HPRMAT solver types
!------------------------------------------------------------------------------
program test_solver
  use precision
  use rmat_solvers
  implicit none

  integer, parameter :: nch = 10    ! Number of channels
  integer, parameter :: nlag = 50   ! Lagrange mesh points
  integer :: ntotal, i, j, stype
  real(dp) :: normfac, t1, t2, elapsed

  complex(dp), allocatable :: cmat(:,:), B_vector(:), Rmat(:,:)
  complex(dp), allocatable :: Rmat_ref(:,:)
  real(dp) :: max_diff

  ntotal = nch * nlag
  normfac = 1.0_dp

  write(*,*) '================================================'
  write(*,*) 'HPRMAT Solver Benchmark'
  write(*,*) '================================================'
  write(*,'(A,I0,A,I0,A,I0)') ' Matrix size: ', ntotal, ' x ', ntotal, ' (nch=', nch
  write(*,'(A,I0,A)') '              nlag=', nlag, ')'
  write(*,*) ''

  allocate(cmat(ntotal, ntotal))
  allocate(B_vector(nlag))
  allocate(Rmat(nch, nch))
  allocate(Rmat_ref(nch, nch))

  ! Initialize test matrix (symmetric positive definite)
  call random_seed()
  do j = 1, ntotal
    do i = 1, ntotal
      cmat(i, j) = cmplx(rand(), rand(), dp)
    end do
  end do
  ! Make it diagonally dominant for stability
  do i = 1, ntotal
    cmat(i, i) = cmat(i, i) + cmplx(ntotal * 2.0_dp, 0.0_dp, dp)
  end do

  ! Initialize B vector
  do i = 1, nlag
    B_vector(i) = cmplx(1.0_dp / sqrt(real(i, dp)), 0.0_dp, dp)
  end do

  ! Test each solver
  write(*,*) 'Solver Type   Time (s)    Max Diff from Ref'
  write(*,*) '------------------------------------------------'

  ! Solver 1: Dense LAPACK (reference)
  call cpu_time(t1)
  call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat_ref, solver_type=1)
  call cpu_time(t2)
  elapsed = t2 - t1
  write(*,'(A,I1,A,F10.4,A)') ' solver_type=', 1, '   ', elapsed, '     (reference)'

  ! Solver 2: Mixed Precision
  call cpu_time(t1)
  call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=2)
  call cpu_time(t2)
  elapsed = t2 - t1
  max_diff = maxval(abs(Rmat - Rmat_ref))
  write(*,'(A,I1,A,F10.4,A,E10.3)') ' solver_type=', 2, '   ', elapsed, '     ', max_diff

  ! Solver 3: Woodbury-Kinetic
  call cpu_time(t1)
  call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=3)
  call cpu_time(t2)
  elapsed = t2 - t1
  max_diff = maxval(abs(Rmat - Rmat_ref))
  write(*,'(A,I1,A,F10.4,A,E10.3)') ' solver_type=', 3, '   ', elapsed, '     ', max_diff

  ! Solver 4: GPU (will fallback to CPU if not available)
  call cpu_time(t1)
  call solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type=4)
  call cpu_time(t2)
  elapsed = t2 - t1
  max_diff = maxval(abs(Rmat - Rmat_ref))
  write(*,'(A,I1,A,F10.4,A,E10.3)') ' solver_type=', 4, '   ', elapsed, '     ', max_diff

  write(*,*) '------------------------------------------------'
  write(*,*) ''
  write(*,*) 'Test completed successfully!'

  deallocate(cmat, B_vector, Rmat, Rmat_ref)

end program test_solver
