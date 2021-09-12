      module param

      ! Master module containing key roms variables used throughout the code
      ! and by other modules. Direct copy of param.h into module format.

! Dimensions of Physical Grid and array dimensions:
!----------- -- -------- ---- --- ----- -----------
! LLm   Number of the internal points of the PHYSICAL grid in XI-
! MMm   and ETA-directions, excluding physical side boundary points,
!       peroodic ghost points, and MPI-margins (if any).
!
! Lm    Number of the internal points [see above] of array covering
! Mm    a single MPI-subdomain.  These are identical to LLm, MMm if
!       there is no MPI-partitioning.

      implicit none

! DevinD: normally this was above param.h in subroutines but since compiled alone need here
#include "cppdefs.h"

      integer, parameter ::
#if defined ANA_WEC_FRC
     &               LLm=1400, MMm=2, N=2    ! DevinD added
#elif defined USWC_WEC || defined USWC_sample
c     &               LLm=1600, MMm=800, N=50       ! US West Coast L2 DD
c     &               LLm=1700, MMm=850, N=50       ! US West Coast L3 DD
     &               LLm=199, MMm=99, N=50       ! test Devin L3 WEC DD

#elif defined PACIFIC_PD
     &               LLm=930, MMm=480, N=60    ! Incorrect grid
C     &               LLm=920, MMm=480, N=60    ! Pacific model PierreD 25km
c    &               LLm=1840, MMm=960, N=100  ! Pacific model PierreD 12.5km
#elif defined ANA_RIVER_USWC
     &               LLm=100, MMm=100, N=10 ! DevinD created USWC_sample
#else
     &                LLm=??, MMm=??, N=??
#endif


! Domain subdivision parameters:
!------- ----------- -----------
! NNODES             total number of MPI processes (nodes);
! NP_XI,  NP_ETA     number of MPI subdomains in XI-, ETA-directions;
! NSUB_X, NSUB_E     number of shared memory subdomains (tiles) in
!                                             XI- and ETA-directions;
      integer, parameter ::
#ifdef MPI
# if defined ANA_WEC_FRC
     &      NP_XI=8, NP_ETA=1, NSUB_X=1, NSUB_E=1  ! DevinD - analytical WEC
# elif defined USWC_WEC || defined USWC_sample || defined ANA_RIVER_USWC
     &      NP_XI=3, NP_ETA=2,  NSUB_X=1, NSUB_E=1 ! DevinD - WEC
# elif defined ANA_RIVER_USWC
	 &      NP_XI=3, NP_ETA=2,  NSUB_X=1, NSUB_E=1
# elif defined PACIFIC_PD
c     &      NP_XI=10, NP_ETA=5, NSUB_X=1, NSUB_E=1
     &      NP_XI=8, NP_ETA=6, NSUB_X=1, NSUB_E=1
# endif
#else
     &      NSUB_X=??, NSUB_E=??
#endif

! Array dimensions and bounds of the used portions of sub-arrays

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
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
# ifndef EW_PERIODIC
     &                , west_exchng,  east_exchng
# endif
# ifndef NS_PERIODIC
     &                , south_exchng, north_exchng
# endif
#else
      integer, parameter :: Lm=LLm, Mm=MMm
#endif

! Derived dimension parameters, number of tracers and tracer
! identification indices:

      integer, parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
#ifdef SOLVE3D
     &       , itemp=1
# ifdef SALINITY
     &       , isalt=2
#  if defined BIOLOGY_BEC2
     &       , ntrc_pas=0
     &       , ntrc_salt=1
     &       , itrc_bio=isalt+1  ! itemp+ntrc_salt+ntrc_pas+1
     &       , ntrc_bio_base=26

     &       , ntrc_bio_ncycle=
#   ifdef Ncycle_SY
     &  3
#   else /* not Ncycle_SY */
     &  0
#   endif /* Ncycle_SY */

     &       , ntrc_bio=ntrc_bio_base+ntrc_bio_ncycle	 ! +ntrc_bio_cocco
     &       , NT=isalt+ntrc_pas+ntrc_bio
#  else /* not BIOLOGY_BEC2 */
     &       , NT=2
#  endif /* BIOLOGY_BEC2 */

# else /* SALINITY */
     &       , NT=1
# endif
#endif /* SOLVE3D */

      end module
