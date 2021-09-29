      module param

      ! Master module containing key roms variables used throughout the code
      ! and by other modules. Direct copy of param.h into module format.

      implicit none

! DevinD: normally this was above param.opt in subroutines but since compiled alone need it here
#include "cppdefs.opt"

! Need this here above param.opt for BGC to work:
      integer, parameter ::
#ifdef SOLVE3D
     &         itemp=1
# ifdef SALINITY
     &       , isalt=2
# endif
#endif

#include "param.opt"

! Array dimensions and bounds of the used portions of sub-arrays:

#ifdef MPI
      integer, parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA

      integer ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
# ifndef EW_PERIODIC
      logical west_exchng,  east_exchng
# endif
# ifndef NS_PERIODIC
      logical south_exchng, north_exchng
# endif
#else
      integer, parameter :: Lm=LLm, Mm=MMm
#endif

! Derived dimension parameters, number of tracers and tracer
! identification indices:

      integer, parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2

      end module
