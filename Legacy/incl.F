! This is to test all .h include files to make sure that everything is
! declared properly: every object declared in every .h-file is either a
! parameter or is a variable inside a common block. Basically the idea
! is to include everything and then compile it relying on compiler
! capability to warn about unused variables, e.g.,
!
!         f95 -Wall ....
!
! (where f95 is gfortran compiler from GCC). The compiler issues a
! warning about unused variable if it is on the "loose" (not in common
! or parameter statement).  The sequence of the include statements is
! arbitrary, except the first three lines with "cppdefs.h", "param.h",
! and ALL_DATA which should be on top.! Defining ALL_DATA forces
! declaration of restricted-access netcdf-data related variables in
! several h-files related to forcing.

#include "cppdefs.h"
#define ALL_DATA

      end

