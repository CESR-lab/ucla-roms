      module param

      ! Master module containing key roms variables used throughout the code
      ! and by other modules. Direct copy of param.h into module format.

      implicit none

! DevinD: normally this was above param.opt in subroutines but since compiled alone need it here
#include "cppdefs.opt"

! Need this here above param.opt for BGC to work:
      integer, parameter :: itemp=1,isalt=2

#include "param.opt"

! Array dimensions and bounds of the used portions of sub-arrays:

#ifdef MPI
      integer, parameter :: nnodes=NP_XI*NP_ETA

      integer ::
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA

      integer ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
#else
      integer, parameter :: Lm=LLm, Mm=MMm
#endif

! Derived dimension parameters, number of tracers and tracer
! identification indices:

!     integer, parameter :: padd_X=(Lm_old+2)/2-(Lm_old+1)/2,
!    &                      padd_E=(Mm_old+2)/2-(Mm_old+1)/2
      integer, parameter :: padd_X= 0,padd_E= 0

# ifdef EW_PERIODIC
      logical,parameter :: ew_periodic =.true.
# else
      logical,parameter :: ew_periodic =.false.
# endif
# ifdef NS_PERIODIC
      logical,parameter :: ns_periodic =.true.
# else
      logical,parameter :: ns_periodic =.false.
# endif

# ifdef OBC_WEST
      logical,parameter :: obc_west =.true.
# else
      logical,parameter :: obc_west =.false.
# endif
# ifdef OBC_EAST
      logical,parameter :: obc_east =.true.
# else
      logical,parameter :: obc_east =.false.
# endif
# ifdef OBC_NORTH
      logical,parameter :: obc_north=.true.
# else
      logical,parameter :: obc_north=.false.
# endif
# ifdef OBC_SOUTH
      logical,parameter :: obc_south=.true.
# else
      logical,parameter :: obc_south=.false.
# endif


      end module
