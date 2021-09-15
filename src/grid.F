      module grid

#include "cppdefs.h"

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
      real h(GLOBAL_2D_ARRAY)
      real hinv(GLOBAL_2D_ARRAY)
      real f(GLOBAL_2D_ARRAY)
      real fomn(GLOBAL_2D_ARRAY)
# ifdef NON_TRADITIONAL
      real f_XI(GLOBAL_2D_ARRAY)
      real f_ETA(GLOBAL_2D_ARRAY)
# endif
 
# ifdef CURVGRID
      real angler(GLOBAL_2D_ARRAY)
# endif
 
#ifdef SPHERICAL
      real latr(GLOBAL_2D_ARRAY)
      real lonr(GLOBAL_2D_ARRAY)
#else
      real xp(GLOBAL_2D_ARRAY)
      real xr(GLOBAL_2D_ARRAY)
      real yp(GLOBAL_2D_ARRAY)
      real yr(GLOBAL_2D_ARRAY)
#endif
 
      real pm(GLOBAL_2D_ARRAY)
      real pn(GLOBAL_2D_ARRAY)
      real dm_r(GLOBAL_2D_ARRAY)
      real dn_r(GLOBAL_2D_ARRAY)
      real pn_u(GLOBAL_2D_ARRAY)
      real dm_u(GLOBAL_2D_ARRAY)
      real dn_u(GLOBAL_2D_ARRAY)
      real dm_v(GLOBAL_2D_ARRAY)
      real pm_v(GLOBAL_2D_ARRAY)
      real dn_v(GLOBAL_2D_ARRAY)
      real dm_p(GLOBAL_2D_ARRAY)
      real dn_p(GLOBAL_2D_ARRAY)

      real iA_u(GLOBAL_2D_ARRAY)
      real iA_v(GLOBAL_2D_ARRAY)
 
#if (defined CURVGRID && defined UV_ADV)
      real dmde(GLOBAL_2D_ARRAY)
      real dndx(GLOBAL_2D_ARRAY)
#endif
      real pmon_u(GLOBAL_2D_ARRAY)
      real pnom_v(GLOBAL_2D_ARRAY)
      real grdscl(GLOBAL_2D_ARRAY)
 
#ifdef MASKING
      real rmask(GLOBAL_2D_ARRAY)
      real pmask(GLOBAL_2D_ARRAY)
      real umask(GLOBAL_2D_ARRAY)
      real vmask(GLOBAL_2D_ARRAY)
#endif

      end module grid
