!------------------------------------------------------------------------------
! HPRMAT C Interface
!
! This module provides C-compatible bindings for the HPRMAT R-matrix solver
! using ISO_C_BINDING for interoperability with C, C++, Rust, and other
! languages that need a C ABI.
!
! Author: Jin Lei
! Date: December 2025
!------------------------------------------------------------------------------
module hprmat_c_interface
  use, intrinsic :: iso_c_binding
  use precision
  use rmat_hp_mod, only: solver_type
  implicit none

contains

  !----------------------------------------------------------------------------
  ! hprmat_init: Initialize R-matrix calculation
  !
  ! C signature:
  !   void hprmat_init(int nr, int ns, double rmax, double* zrma);
  !----------------------------------------------------------------------------
  subroutine hprmat_init(nr, ns, rmax, zrma) bind(C, name='hprmat_init')
    integer(c_int), intent(in), value :: nr, ns
    real(c_double), intent(in), value :: rmax
    real(c_double), intent(out) :: zrma(nr*ns)

    call rmat_ini(nr, ns, rmax, zrma)
  end subroutine hprmat_init

  !----------------------------------------------------------------------------
  ! hprmat_set_solver: Set the solver type
  !
  ! C signature:
  !   void hprmat_set_solver(int solver);
  !----------------------------------------------------------------------------
  subroutine hprmat_set_solver(solver) bind(C, name='hprmat_set_solver')
    integer(c_int), intent(in), value :: solver
    solver_type = solver
  end subroutine hprmat_set_solver

  !----------------------------------------------------------------------------
  ! hprmat_get_solver: Get the current solver type
  !
  ! C signature:
  !   int hprmat_get_solver(void);
  !----------------------------------------------------------------------------
  function hprmat_get_solver() result(solver) bind(C, name='hprmat_get_solver')
    integer(c_int) :: solver
    solver = solver_type
  end function hprmat_get_solver

  !----------------------------------------------------------------------------
  ! hprmat_solve: Main R-matrix solver (simplified interface)
  !
  ! C signature:
  !   void hprmat_solve(int nch, int* lval, double* qk, double* eta,
  !                     double rmax, int nr, int ns,
  !                     double complex* cpot, int cpot_dim1,
  !                     double complex* cu, int* nopen,
  !                     int solver);
  !----------------------------------------------------------------------------
  subroutine hprmat_solve(nch, lval, qk, eta, rmax, nr, ns, &
                          cpot, cpot_dim1, cu, nopen, solver) &
                          bind(C, name='hprmat_solve')
    integer(c_int), intent(in), value :: nch, nr, ns, cpot_dim1, solver
    integer(c_int), intent(in) :: lval(nch)
    real(c_double), intent(in) :: qk(nch), eta(nch)
    real(c_double), intent(in), value :: rmax
    complex(c_double_complex), intent(in) :: cpot(cpot_dim1, nch, nch)
    complex(c_double_complex), intent(out) :: cu(nch, nch)
    integer(c_int), intent(out) :: nopen

    ! Local variables
    integer :: lval_f(nch), nvc_f(1), nopen_f, local_solver
    real*8 :: qk_f(nch), eta_f(nch)
    complex*16 :: cpot_f(cpot_dim1, nch, nch), cu_f(nch, nch)
    complex*16 :: cf_dummy(1, 1, 1), cpnl_dummy(1, 1, 1)
    logical :: twf

    ! Copy to Fortran types
    lval_f = lval
    qk_f = qk
    eta_f = eta
    cpot_f = cpot
    twf = .false.
    nvc_f(1) = 1

    ! Determine solver
    if (solver > 0) then
      local_solver = solver
    else
      local_solver = solver_type
    end if

    ! Call main rmatrix routine
    call rmatrix(nch, lval_f, qk_f, eta_f, rmax, nr, ns, cpot_f, cu_f, &
                 cpot_dim1, nch, nopen_f, twf, cf_dummy, 1, 1, 1, nvc_f, &
                 0, cpnl_dummy, local_solver)

    ! Copy output
    cu = cu_f
    nopen = nopen_f

  end subroutine hprmat_solve

  !----------------------------------------------------------------------------
  ! hprmat_solve_full: Full R-matrix solver with wave function output
  !
  ! C signature:
  !   void hprmat_solve_full(int nch, int* lval, double* qk, double* eta,
  !                          double rmax, int nr, int ns,
  !                          double complex* cpot, int cpot_dim1,
  !                          double complex* cu,
  !                          double complex* cf, int cf_dim1, int cf_dim2,
  !                          int nc_entrance, int* nvc,
  !                          double complex* cpnl, int cpnl_dim1,
  !                          int* nopen, int solver);
  !----------------------------------------------------------------------------
  subroutine hprmat_solve_full(nch, lval, qk, eta, rmax, nr, ns, &
                               cpot, cpot_dim1, cu, &
                               cf, cf_dim1, cf_dim2, nc_entrance, nvc, &
                               cpnl, cpnl_dim1, nopen, solver) &
                               bind(C, name='hprmat_solve_full')
    integer(c_int), intent(in), value :: nch, nr, ns, cpot_dim1, cf_dim1, cf_dim2
    integer(c_int), intent(in), value :: nc_entrance, cpnl_dim1, solver
    integer(c_int), intent(in) :: lval(nch), nvc(nc_entrance)
    real(c_double), intent(in) :: qk(nch), eta(nch)
    real(c_double), intent(in), value :: rmax
    complex(c_double_complex), intent(in) :: cpot(cpot_dim1, nch, nch)
    complex(c_double_complex), intent(in) :: cpnl(max(1,cpnl_dim1), nch, nch)
    complex(c_double_complex), intent(out) :: cu(nch, nch)
    complex(c_double_complex), intent(out) :: cf(cf_dim1, cf_dim2, nc_entrance)
    integer(c_int), intent(out) :: nopen

    ! Local variables
    integer :: lval_f(nch), nvc_f(nc_entrance), nopen_f, local_solver
    real*8 :: qk_f(nch), eta_f(nch)
    complex*16, allocatable :: cpot_f(:,:,:), cu_f(:,:), cf_f(:,:,:), cpnl_f(:,:,:)
    logical :: twf

    allocate(cpot_f(cpot_dim1, nch, nch))
    allocate(cu_f(nch, nch))
    allocate(cf_f(cf_dim1, cf_dim2, nc_entrance))

    ! Copy to Fortran types
    lval_f = lval
    nvc_f = nvc
    qk_f = qk
    eta_f = eta
    cpot_f = cpot
    twf = .true.

    ! Determine solver
    if (solver > 0) then
      local_solver = solver
    else
      local_solver = solver_type
    end if

    if (cpnl_dim1 > 0) then
      allocate(cpnl_f(cpnl_dim1, nch, nch))
      cpnl_f = cpnl(1:cpnl_dim1, :, :)
      call rmatrix(nch, lval_f, qk_f, eta_f, rmax, nr, ns, cpot_f, cu_f, &
                   cpot_dim1, nch, nopen_f, twf, cf_f, cf_dim1, cf_dim2, &
                   nc_entrance, nvc_f, cpnl_dim1, cpnl_f, local_solver)
      deallocate(cpnl_f)
    else
      allocate(cpnl_f(1, nch, nch))
      call rmatrix(nch, lval_f, qk_f, eta_f, rmax, nr, ns, cpot_f, cu_f, &
                   cpot_dim1, nch, nopen_f, twf, cf_f, cf_dim1, cf_dim2, &
                   nc_entrance, nvc_f, 0, cpnl_f, local_solver)
      deallocate(cpnl_f)
    end if

    ! Copy output
    cu = cu_f
    cf = cf_f
    nopen = nopen_f

    deallocate(cpot_f, cu_f, cf_f)

  end subroutine hprmat_solve_full

  !----------------------------------------------------------------------------
  ! hprmat_wavefunction: Calculate wave function on uniform mesh
  !
  ! C signature:
  !   void hprmat_wavefunction(int nch, int* lval, double* qk, double* eta,
  !                            double rmax, int nr, int ns,
  !                            double complex* cu, int nopen,
  !                            double complex* cf, int cf_dim1, int cf_dim2,
  !                            double* zrma, int iv, int nom,
  !                            int npoin, double h,
  !                            double complex* cwftab);
  !----------------------------------------------------------------------------
  subroutine hprmat_wavefunction(nch, lval, qk, eta, rmax, nr, ns, &
                                 cu, nopen_in, cf, cf_dim1, cf_dim2, &
                                 zrma, iv, nom, npoin, h, cwftab) &
                                 bind(C, name='hprmat_wavefunction')
    integer(c_int), intent(in), value :: nch, nr, ns, nopen_in, cf_dim1, cf_dim2
    integer(c_int), intent(in), value :: iv, nom, npoin
    integer(c_int), intent(in) :: lval(nch)
    real(c_double), intent(in) :: qk(nch), eta(nch), zrma(nr*ns)
    real(c_double), intent(in), value :: rmax, h
    complex(c_double_complex), intent(in) :: cu(nch, nch), cf(cf_dim1, cf_dim2, nom)
    complex(c_double_complex), intent(out) :: cwftab(npoin)

    ! Local variables
    integer :: lval_f(nch)
    real*8 :: qk_f(nch), eta_f(nch), zrma_f(nr*ns)
    complex*16 :: cu_f(nch, nch), cf_f(cf_dim1, cf_dim2, nom), cwftab_f(npoin)

    lval_f = lval
    qk_f = qk
    eta_f = eta
    zrma_f = zrma
    cu_f = cu
    cf_f = cf

    call wf_print(nch, lval_f, qk_f, eta_f, rmax, nr, ns, cu_f, &
                  nch, nopen_in, cf_f, cf_dim1, cf_dim2, zrma_f, &
                  iv, nom, npoin, h, cwftab_f)

    cwftab = cwftab_f

  end subroutine hprmat_wavefunction

  !----------------------------------------------------------------------------
  ! Low-level solvers (direct access for advanced users)
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! hprmat_solve_linear: Solve linear system for R-matrix
  !
  ! C signature:
  !   void hprmat_solve_linear(double complex* cmat, double complex* B_vector,
  !                            int nch, int nlag, double normfac,
  !                            double complex* Rmat, int solver);
  !----------------------------------------------------------------------------
  subroutine hprmat_solve_linear(cmat, B_vector, nch, nlag, normfac, Rmat, solver) &
                                 bind(C, name='hprmat_solve_linear')
    use rmat_solvers, only: solve_rmatrix
    integer(c_int), intent(in), value :: nch, nlag, solver
    real(c_double), intent(in), value :: normfac
    complex(c_double_complex), intent(in) :: cmat(nch*nlag, nch*nlag)
    complex(c_double_complex), intent(in) :: B_vector(nlag)
    complex(c_double_complex), intent(out) :: Rmat(nch, nch)

    complex(dp) :: cmat_f(nch*nlag, nch*nlag), B_f(nlag), Rmat_f(nch, nch)

    cmat_f = cmat
    B_f = B_vector

    call solve_rmatrix(cmat_f, B_f, nch, nlag, real(normfac, dp), Rmat_f, solver)

    Rmat = Rmat_f

  end subroutine hprmat_solve_linear

  !----------------------------------------------------------------------------
  ! Utility functions
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! hprmat_version: Get version string
  !
  ! C signature:
  !   void hprmat_version(char* version, int maxlen);
  !----------------------------------------------------------------------------
  subroutine hprmat_version(version, maxlen) bind(C, name='hprmat_version')
    character(kind=c_char), intent(out) :: version(*)
    integer(c_int), intent(in), value :: maxlen

    character(len=32) :: ver_str
    integer :: i

    ver_str = 'HPRMAT v1.0.0'

    do i = 1, min(len_trim(ver_str), maxlen-1)
      version(i) = ver_str(i:i)
    end do
    version(min(len_trim(ver_str)+1, maxlen)) = c_null_char

  end subroutine hprmat_version

  !----------------------------------------------------------------------------
  ! hprmat_solver_info: Get solver description
  !
  ! C signature:
  !   void hprmat_solver_info(int solver, char* info, int maxlen);
  !----------------------------------------------------------------------------
  subroutine hprmat_solver_info(solver, info, maxlen) bind(C, name='hprmat_solver_info')
    integer(c_int), intent(in), value :: solver
    character(kind=c_char), intent(out) :: info(*)
    integer(c_int), intent(in), value :: maxlen

    character(len=64) :: info_str
    integer :: i

    select case (solver)
    case (1)
      info_str = 'Dense LAPACK ZGESV (reference)'
    case (2)
      info_str = 'Mixed Precision (single + double refinement)'
    case (3)
      info_str = 'Woodbury-Kinetic (CPU optimized)'
    case (4)
      info_str = 'GPU cuSOLVER (NVIDIA GPU)'
    case default
      info_str = 'Unknown solver type'
    end select

    do i = 1, min(len_trim(info_str), maxlen-1)
      info(i) = info_str(i:i)
    end do
    info(min(len_trim(info_str)+1, maxlen)) = c_null_char

  end subroutine hprmat_solver_info

end module hprmat_c_interface
