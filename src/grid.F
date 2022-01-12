      module grid

#include "cppdefs.opt"

      use param
      implicit none

! This is include file "grid": Environmental two-dimensional arrays
! associated with curvilinear horizontal coordinate system.
!
! h         Model bottom topography (depth [m] at RHO-points.)
! f, fomn   Coriolis parameter [1/s] and compound term f/[pm*pn]
!                                                   at RHO points.
! angler      Angle [radians] between XI-axis and the direction
!                                       to the EAST at RHO-points.
! latr, lonr  Latitude (degrees north) and Longitude (degrees east)
!                                                  at RHO-points.
! xr, xp      XI-coordinates [m] at RHO- and PSI-points.
! yr, yp      ETA-coordinates [m] at RHO- and PSI-points.
!
! pm, pm  Coordinate transformation metric "m" and "n" associated
!         with the differential distances in XI- and ETA-directions.
!
! dm_u, dm_r  Grid spacing [meters] in the XI-direction
! dm_v, dm_p       at U-, RHO-,  V- and vorticity points.
! dn_u, dn_r  Grid spacing [meters] in the ETA-direction
! dn_v, dn_p      at U-, RHO-,  V- and vorticity points.
!
! dmde     ETA-derivative of inverse metric factor "m" d(1/M)/d(ETA)
! dndx     XI-derivative  of inverse metric factor "n" d(1/N)/d(XI)
!
! pmon_u   Compound term, pm/pn at U-points.
! pnom_v   Compound term, pn/pm at V-points.
!
! umask, rmask  Land-sea masking arrays at RHO-,U-,V- and PSI-points
! pmask, vmask      (rmask,umask,vmask) = (0=Land, 1=Sea);
!                    pmask = (0=Land, 1=Sea, 1-gamma2 =boundary).
!
      real,allocatable,dimension(:,:) :: h
      real,allocatable,dimension(:,:) :: hinv
      real,allocatable,dimension(:,:) :: f
      real,allocatable,dimension(:,:) :: fomn
# ifdef NON_TRADITIONAL
      real,allocatable,dimension(:,:) :: f_XI
      real,allocatable,dimension(:,:) :: f_ETA
# endif

# ifdef CURVGRID
      real,allocatable,dimension(:,:) :: angler
# endif

#ifdef SPHERICAL
      real,allocatable,dimension(:,:) :: latr
      real,allocatable,dimension(:,:) :: lonr
#else
      real,allocatable,dimension(:,:) :: xp
      real,allocatable,dimension(:,:) :: xr
      real,allocatable,dimension(:,:) :: yp
      real,allocatable,dimension(:,:) :: yr
#endif

      real,allocatable,dimension(:,:) :: pm
      real,allocatable,dimension(:,:) :: pn
      real,allocatable,dimension(:,:) :: dm_r
      real,allocatable,dimension(:,:) :: dn_r
      real,allocatable,dimension(:,:) :: pn_u
      real,allocatable,dimension(:,:) :: dm_u
      real,allocatable,dimension(:,:) :: dn_u
      real,allocatable,dimension(:,:) :: dm_v
      real,allocatable,dimension(:,:) :: pm_v
      real,allocatable,dimension(:,:) :: dn_v
      real,allocatable,dimension(:,:) :: dm_p
      real,allocatable,dimension(:,:) :: dn_p

      real,allocatable,dimension(:,:) :: iA_u
      real,allocatable,dimension(:,:) :: iA_v

#if (defined CURVGRID && defined UV_ADV)
      real,allocatable,dimension(:,:) :: dmde
      real,allocatable,dimension(:,:) :: dndx
#endif
      real,allocatable,dimension(:,:) :: pmon_u
      real,allocatable,dimension(:,:) :: pnom_v
      real,allocatable,dimension(:,:) :: grdscl

#ifdef MASKING
      real,allocatable,dimension(:,:) :: rmask
      real,allocatable,dimension(:,:) :: pmask
      real,allocatable,dimension(:,:) :: umask
      real,allocatable,dimension(:,:) :: vmask
#endif

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_grid  ![
      implicit none

      allocate( h(GLOBAL_2D_ARRAY) )     ! Potential first touch issue. Before only rmask was set in init_arrays...
      allocate( hinv(GLOBAL_2D_ARRAY) )
      allocate( f(GLOBAL_2D_ARRAY) )
      allocate( fomn(GLOBAL_2D_ARRAY) )
# ifdef NON_TRADITIONAL
      allocate( f_XI(GLOBAL_2D_ARRAY) )
      allocate( f_ETA(GLOBAL_2D_ARRAY) )
# endif
 
# ifdef CURVGRID
      allocate( angler(GLOBAL_2D_ARRAY) )
# endif
 
#ifdef SPHERICAL
      allocate( latr(GLOBAL_2D_ARRAY) )
      allocate( lonr(GLOBAL_2D_ARRAY) )
#else
      allocate( xp(GLOBAL_2D_ARRAY) )
      allocate( xr(GLOBAL_2D_ARRAY) )
      allocate( yp(GLOBAL_2D_ARRAY) )
      allocate( yr(GLOBAL_2D_ARRAY) )
#endif
 
      allocate( pm(GLOBAL_2D_ARRAY) )
      allocate( pn(GLOBAL_2D_ARRAY) )
      allocate( dm_r(GLOBAL_2D_ARRAY) )
      allocate( dn_r(GLOBAL_2D_ARRAY) )
      allocate( pn_u(GLOBAL_2D_ARRAY) )
      allocate( dm_u(GLOBAL_2D_ARRAY) )
      allocate( dn_u(GLOBAL_2D_ARRAY) )
      allocate( dm_v(GLOBAL_2D_ARRAY) )
      allocate( pm_v(GLOBAL_2D_ARRAY) )
      allocate( dn_v(GLOBAL_2D_ARRAY) )
      allocate( dm_p(GLOBAL_2D_ARRAY) )
      allocate( dn_p(GLOBAL_2D_ARRAY) )

      allocate( iA_u(GLOBAL_2D_ARRAY) )
      allocate( iA_v(GLOBAL_2D_ARRAY) )
 
#if (defined CURVGRID && defined UV_ADV)
      allocate( dmde(GLOBAL_2D_ARRAY) )
      allocate( dndx(GLOBAL_2D_ARRAY) )
#endif
      allocate( pmon_u(GLOBAL_2D_ARRAY) )
      allocate( pnom_v(GLOBAL_2D_ARRAY) )
      allocate( grdscl(GLOBAL_2D_ARRAY) )
 
#ifdef MASKING
      allocate( rmask(GLOBAL_2D_ARRAY) )
      allocate( pmask(GLOBAL_2D_ARRAY) )
      allocate( umask(GLOBAL_2D_ARRAY) )
      allocate( vmask(GLOBAL_2D_ARRAY) )
#endif

      ! WARNING: "rmask" MUST BE initialized to all-one state in order to
      ! read grid variables (coordinates, metric, topography), which should
      ! not be masked.

#ifdef MASKING
      rmask=1.
#endif

      end subroutine init_arrays_grid  !]

      !----------------------------------------------------------------

      end module grid
