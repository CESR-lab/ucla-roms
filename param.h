! Dimensions of Physical Grid and array dimensions:
!----------- -- -------- ---- --- ----- -----------
! LLm   Number of the internal points of the PHYSICAL grid in XI-
! MMm   and ETA-directions, excluding physical side boundary points,
!       peroodic ghost points, and MPI-margins (if any).
!
! Lm    Number of the internal points [see above] of array covering
! Mm    a single MPI-subdomain.  These are identical to LLm, MMm if
!       there is no MPI-partitioning.

      integer, parameter ::
#if defined ANA_WEC_FRC
     &               LLm=1400, MMm=2, N=2    ! DevinD added
#elif defined USWC_WEC || defined USWC_sample
c     &               LLm=1600, MMm=800, N=50       ! US West Coast L2 DD
c     &               LLm=1700, MMm=850, N=50       ! US West Coast L3 DD
     &               LLm=199, MMm=99, N=50       ! test Devin L3 WEC DD

#elif defined PACIFIC_PD
c     &               LLm=930, MMm=480, N=60    ! Incorrect grid
     &               LLm=920, MMm=480, N=60    ! Pacific model PierreD 25km
c    &               LLm=1840, MMm=960, N=100  ! Pacific model PierreD 12.5km
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
# elif defined USWC_WEC || defined USWC_sample
     &      NP_XI=3, NP_ETA=2,  NSUB_X=1, NSUB_E=1 ! DevinD - WEC
#elif defined PACIFIC_PD
     &      NP_XI=10, NP_ETA=5, NSUB_X=1, NSUB_E=1
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
#  ifdef BIOLOGY
     &       , NT=7, iNO3_=3, iNH4_=4, iDet_=5, iPhyt=6, iZoo_=7
#  elif defined USWC_sample
     &       , NT=2
#  else
     &       , NT=2
#  endif
# else
     &       , NT=1
# endif
#endif

#ifdef PSOURCE
     &       , Msrc=10   ! Number of point sources
#endif

! Tides
! -----
!#if defined TIDES
!      integer Ntides   ! Number of tides
!      parameter (Ntides=15)
!#endif
