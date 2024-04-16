#ifdef FLOATS 
! This is include file "ncvars_floats.h".
! ==== == ======= ==== ============
! indices in character array "vname", which holds variable names
!                                                and attributes.
! indxfltGrd      Float grid level
! indxfltTemp     Temperature at float location
! indxfltSalt     Salinity at float location
! indxfltRho      Density at float location
! indxfltVel      Averaged velocity module
!
      integer, parameter ::  fltfield=5, indxfltGrd=1, indxfltTemp=2,
     &                    indxfltSalt=3, indxfltRho=4,  indxfltVel=5

      integer ncidflt,  nrecflt,    nrpfflt,  fltGlevel,
     &        fltTstep, fltTime,    fltXgrd,  fltYgrd,   fltZgrd,
     &        fltVel,   rstnfloats, rstTinfo, rstfltgrd, rsttrack
#ifdef SPHERICAL
     &      , fltLon,   fltLat
#else
     &      , fltX,     fltY
#endif
#ifdef SOLVE3D
     &      , fltDepth, fltDen,     fltTemp
# ifdef SALINITY
     &                        ,     fltSal
# endif
#endif
      logical wrtflt(fltfield)

      common/incvars_floats/ ncidflt, nrecflt, nrpfflt, fltGlevel,
     &        fltTstep, fltTime,    fltXgrd,  fltYgrd,   fltZgrd, 
     &        fltVel,   rstnfloats, rstTinfo, rstfltgrd, rsttrack
#ifdef SPHERICAL
     &      , fltLon,   fltLat
#else
     &      , fltX,     fltY
#endif
#ifdef SOLVE3D
     &      , fltDepth, fltDen,     fltTemp
# ifdef SALINITY
     &                        ,     fltSal
# endif
#endif
     &      , wrtflt

      character*80  fltname,   fposnam
      common /cncvars_floats/ fltname,   fposnam
#endif /*FLOATS*/

