!------------------------------------------------------------------------------
! Module: constants
! Purpose: Physical constants for nuclear physics calculations
!------------------------------------------------------------------------------
module constants
  use precision
  implicit none

  ! Physical constants (CODATA values)
  real(dp), parameter :: hbarc = 197.3269718_dp   ! MeV.fm
  real(dp), parameter :: finec = 137.03599_dp     ! Fine structure constant inverse
  real(dp), parameter :: amu = 931.49432_dp       ! Atomic mass unit in MeV
  real(dp), parameter :: e2 = 1.43997_dp          ! e^2 in MeV.fm

end module constants
