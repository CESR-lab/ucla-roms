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
#include "param.h"

#include "averages.h"
#include "boundary.h"
#include "buffer.h"
#include "climat.h"
#include "coupling.h"
#include "diagnostics.h"
#include "floats.h"
#include "forces.h"
#include "grid.h"
#include "init_floats.h"
#include "mixing.h"
#include "ncvars_floats.h"
#include "ncvars.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "private_scratch.h"
#include "scalars.h"
#include "scoord.h"
#include "sources.h"
#include "strings.h"
#include "tides.h"
#include "work.h"
      end

