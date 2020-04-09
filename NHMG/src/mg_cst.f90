!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !===================================================================!
  !- All constants of the program have to be declared in this module -!
  !===================================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mg_cst

  implicit none

  !======================================!
  !- DECLARATIONS of GLOBAL VARIABLE(S) -!
  !======================================!

  !NG For the moment define as global variable in mg_namelist.f90 file

  integer(kind = 4), parameter :: ip = 4
  integer(kind = 4), parameter :: rp = 8

  real(kind=rp), parameter :: two  = 2._rp
  real(kind=rp), parameter :: one  = 1._rp
  real(kind=rp), parameter :: hlf  = 0.5_rp
  real(kind=rp), parameter :: qrt  = 0.25_rp
  real(kind=rp), parameter :: eighth  = 0.125_rp
  real(kind=rp), parameter :: zero = 0._rp

end module mg_cst
