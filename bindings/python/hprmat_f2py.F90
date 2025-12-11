!------------------------------------------------------------------------------
! HPRMAT f2py Interface (Simplified)
!
! This module provides f2py-compatible wrappers for HPRMAT.
! Only includes core routines that f2py can parse cleanly.
!
! Author: Jin Lei
! Date: December 2025
!------------------------------------------------------------------------------

subroutine py_rmat_init(nr, ns, rmax, zrma)
  !f2py intent(in) nr, ns, rmax
  !f2py intent(out) zrma
  !f2py depend(nr, ns) zrma
  implicit none
  integer, intent(in) :: nr, ns
  real*8, intent(in) :: rmax
  real*8, intent(out) :: zrma(nr*ns)

  call rmat_ini(nr, ns, rmax, zrma)
end subroutine py_rmat_init

subroutine py_set_solver(solver)
  !f2py intent(in) solver
  use rmat_hp_mod, only: solver_type
  implicit none
  integer, intent(in) :: solver
  solver_type = solver
end subroutine py_set_solver

subroutine py_get_solver(solver)
  !f2py intent(out) solver
  use rmat_hp_mod, only: solver_type
  implicit none
  integer, intent(out) :: solver
  solver = solver_type
end subroutine py_get_solver

subroutine py_rmatrix(nch, ntot, lval, qk, eta, rmax, nr, ns, cpot, cu, nopen, isolver)
  !f2py intent(in) nch, ntot, lval, qk, eta, rmax, nr, ns, cpot, isolver
  !f2py intent(out) cu, nopen
  !f2py depend(nch) lval, qk, eta, cu
  !f2py depend(ntot, nch) cpot
  use rmat_hp_mod
  implicit none

  integer, intent(in) :: nch, nr, ns, ntot, isolver
  integer, intent(in) :: lval(nch)
  real*8, intent(in) :: qk(nch), eta(nch), rmax
  complex*16, intent(in) :: cpot(ntot, nch, nch)
  complex*16, intent(out) :: cu(nch, nch)
  integer, intent(out) :: nopen

  ! Local variables
  complex*16, allocatable :: cf_dummy(:,:,:), cpnl_dummy(:,:,:)
  integer :: nvc_dummy(1)
  logical :: twf

  allocate(cf_dummy(1, 1, 1), cpnl_dummy(1, 1, 1))

  twf = .false.
  nvc_dummy(1) = 1

  call rmatrix(nch, lval, qk, eta, rmax, nr, ns, cpot, cu, &
               ntot, nch, nopen, twf, cf_dummy, 1, 1, 1, nvc_dummy, &
               0, cpnl_dummy, isolver)

  deallocate(cf_dummy, cpnl_dummy)

end subroutine py_rmatrix

subroutine py_rmatrix_wf(nch, ntot, lval, qk, eta, rmax, nr, ns, cpot, &
                         cu, cf, nc_entrance, nvc, nopen, isolver)
  !f2py intent(in) nch, ntot, lval, qk, eta, rmax, nr, ns, cpot, nc_entrance, nvc, isolver
  !f2py intent(out) cu, cf, nopen
  !f2py depend(nch) lval, qk, eta, cu
  !f2py depend(ntot, nch) cpot
  !f2py depend(ntot, nch, nc_entrance) cf
  !f2py depend(nc_entrance) nvc
  use rmat_hp_mod
  implicit none

  integer, intent(in) :: nch, nr, ns, nc_entrance, isolver, ntot
  integer, intent(in) :: lval(nch), nvc(nc_entrance)
  real*8, intent(in) :: qk(nch), eta(nch), rmax
  complex*16, intent(in) :: cpot(ntot, nch, nch)
  complex*16, intent(out) :: cu(nch, nch)
  complex*16, intent(out) :: cf(ntot, nch, nc_entrance)
  integer, intent(out) :: nopen

  ! Local variables
  complex*16, allocatable :: cpnl_dummy(:,:,:)
  logical :: twf

  allocate(cpnl_dummy(1, 1, 1))
  twf = .true.

  call rmatrix(nch, lval, qk, eta, rmax, nr, ns, cpot, cu, &
               ntot, nch, nopen, twf, cf, ntot, nch, nc_entrance, nvc, &
               0, cpnl_dummy, isolver)

  deallocate(cpnl_dummy)

end subroutine py_rmatrix_wf
