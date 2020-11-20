! This is include file "grid.h": Environmental two-dimensional arrays
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
CSDISTRIBUTE_RESHAPE h(BLOCK_PATTERN) BLOCK_CLAUSE
      real hinv(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE h(BLOCK_PATTERN) BLOCK_CLAUSE
      real f(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE f(BLOCK_PATTERN) BLOCK_CLAUSE
      real fomn(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE fomn(BLOCK_PATTERN) BLOCK_CLAUSE
      common /grd_h/h /grd_hinv/hinv /grd_f/f /grd_fomn/fomn
# ifdef NON_TRADITIONAL
      real f_XI(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE f_XI(BLOCK_PATTERN) BLOCK_CLAUSE
      real f_ETA(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE f_ETA(BLOCK_PATTERN) BLOCK_CLAUSE
      common /grd_fXI/f_XI /grd_fETA/f_ETA
# endif
 
# ifdef CURVGRID
      real angler(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE angler(BLOCK_PATTERN) BLOCK_CLAUSE
      common /grd_angler/angler
# endif
 
#ifdef SPHERICAL
      real latr(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE   latr(BLOCK_PATTERN) BLOCK_CLAUSE
      real lonr(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE   lonr(BLOCK_PATTERN) BLOCK_CLAUSE
      common /grd_latr/latr /grd_lonr/lonr
#else
      real xp(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE xp(BLOCK_PATTERN) BLOCK_CLAUSE
      real xr(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE xr(BLOCK_PATTERN) BLOCK_CLAUSE
      real yp(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE yp(BLOCK_PATTERN) BLOCK_CLAUSE
      real yr(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE yr(BLOCK_PATTERN) BLOCK_CLAUSE
      common /grd_xr/xr /grd_xp/xp /grd_yp/yp /grd_yr/yr
#endif
 
      real pm(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE pm(BLOCK_PATTERN) BLOCK_CLAUSE
      real pn(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE pn(BLOCK_PATTERN) BLOCK_CLAUSE
      real dm_r(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dm_r(BLOCK_PATTERN) BLOCK_CLAUSE
      real dn_r(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dn_r(BLOCK_PATTERN) BLOCK_CLAUSE
      real pn_u(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE pn_u(BLOCK_PATTERN) BLOCK_CLAUSE
      real dm_u(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dm_u(BLOCK_PATTERN) BLOCK_CLAUSE
      real dn_u(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dn_u(BLOCK_PATTERN) BLOCK_CLAUSE
      real dm_v(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dm_v(BLOCK_PATTERN) BLOCK_CLAUSE
      real pm_v(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE pm_v(BLOCK_PATTERN) BLOCK_CLAUSE
      real dn_v(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dn_v(BLOCK_PATTERN) BLOCK_CLAUSE
      real dm_p(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dm_p(BLOCK_PATTERN) BLOCK_CLAUSE
      real dn_p(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dn_p(BLOCK_PATTERN) BLOCK_CLAUSE
      common /mtrix_pm/pm     /mtrix_pn/pn
     &       /mtrix_dm_r/dm_r /mtrix_dn_r/dn_r
     &       /mtrix_dm_u/dm_u /mtrix_dn_u/dn_u
     &       /mtrix_dm_v/dm_v /mtrix_dn_v/dn_v
     &       /mtrix_pn_u/pn_u /mtrix_pm_v/pm_v
     &       /mtrix_dm_p/dm_p /mtrix_dn_p/dn_p

      real iA_u(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE iA_u(BLOCK_PATTERN) BLOCK_CLAUSE
      real iA_v(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE iA_v(BLOCK_PATTERN) BLOCK_CLAUSE
      common /mtrix_iAu/iA_u  /mtrix_iAv/iA_v
 
#if (defined CURVGRID && defined UV_ADV)
      real dmde(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dmde(BLOCK_PATTERN) BLOCK_CLAUSE
      real dndx(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE dndx(BLOCK_PATTERN) BLOCK_CLAUSE
      common /mtrix_dmde/dmde   /mtrix_dndx/dndx
#endif
      real pmon_u(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE pmon_u(BLOCK_PATTERN) BLOCK_CLAUSE
      real pnom_v(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE pnom_v(BLOCK_PATTERN) BLOCK_CLAUSE
      real grdscl(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE grdscl(BLOCK_PATTERN) BLOCK_CLAUSE
      common /mtrix_pmon_u/pmon_u /mtrix_pnom_v/pnom_v
     &                            /mtrix_grdscl/grdscl
 
#ifdef MASKING
      real rmask(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE rmask(BLOCK_PATTERN) BLOCK_CLAUSE
      real pmask(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE pmask(BLOCK_PATTERN) BLOCK_CLAUSE
      real umask(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE umask(BLOCK_PATTERN) BLOCK_CLAUSE
      real vmask(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE vmask(BLOCK_PATTERN) BLOCK_CLAUSE
      common /mask_r/rmask /mask_p/pmask
     &       /mask_u/umask /mask_v/vmask
#endif
