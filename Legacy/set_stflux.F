#include "cppdefs.h"
! DevinD added BULK_FLUX as set surface tracers accounted for in flux modules
! when BULK_FLUX is turned on
#if defined SOLVE3D && (!defined ANA_STFLUX || !defined ANA_SSFLUX) && !defined BULK_FLUX && !defined FLUX_FRC

      subroutine set_stflux_tile(istr,iend,jstr,jend, itrc, ierr)

! Interpolate and assemble surface tracer flux for current tile.

      use tracers, only: t

      implicit none
      integer istr,iend,jstr,jend, itrc, ierr, i,j, it1,it2
      real tmid, cff1,cff2, cff
      logical out_of_bounds

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"

# define STFLUX_DATA
# if defined QCORRECTION && !defined ANA_SST
      integer it3,it4, it5,it6
      real  cff3,cff4, cff5,cff6
#  define SST_DATA
#  define DQDT_DATA
# endif
# if defined SFLX_CORR && defined SALINITY
      integer  it7,it8
      real cff7,cff8
#  define SSS_DATA
# endif
# include "forces.h"

# ifdef SEA_ICE_NOFLUX
#  include "mixing.h"
# endif

# include "compute_extended_bounds.h"

! Set coefficients for interpolation.  Also check that for the next
! time step [when time=time+dt] both weights will still be positive --
! if not set synchro_flag to signal that a fresh timeslice of data
! should be read from an appropriate netCDF input file (master thread
! only). After that either load time-invariant data, or interpolate in
! time or complain about error and signal to quit, if interpolation
! is needed, but not possible.  Note unit conversion factors:

! for heat flux: 1./(rho0*Cp) to convert [W/(m^2*deg C)] into [m/sec]
! thus into units of kinematic flux; same applies to (T-SST)*dQdSST;

! for salinity(freshwater) flux: 0.01/day2sec to convert from [cm/day]
! into [PSU m/s]

      if (ntstf(itrc) > 1) then  !--> interpolate in time
        tmid=time+0.5*dt  ;  out_of_bounds=.false.
        it2=itstf(itrc)   ;  cff1=stf_time(it2,itrc)-tmid
        it1=3-it2         ;  cff2=tmid-stf_time(it1,itrc)
        if (cff1 < 0. .or. cff2 < 0.) out_of_bounds=.true.
        if (cff1 < dt .and. ntstf(itrc) > 1) synchro_flag=.true.
        if (itrc == itemp) then
          cff=1./(rho0*Cp*(cff1+cff2))
# ifdef SALINITY
        elseif (itrc == isalt) then
          cff=0.01/(day2sec*(cff1+cff2))
# endif
        else
          cff=1./(cff1+cff2)
        endif
        cff1=cff1*cff ; cff2=cff2*cff

# if defined QCORRECTION && !defined ANA_SST
        it4=itsst  ;   cff3=sst_time(it4)-tmid
        it3=3-it4  ;   cff4=tmid-sst_time(it3)
        if (cff3 < 0. .or. cff4 < 0.) out_of_bounds=.true.
        if (cff3 < dt .and. ntstf(itrc) > 1) synchro_flag=.true.
        cff=1./(cff3+cff4) ; cff3=cff3*cff ; cff4=cff4*cff

        it6=itdqdt  ;  cff5=dqdt_time(it6)-tmid
        it5=3-it6   ;  cff6=tmid-dqdt_time(it5)
        if (cff5 < 0. .or. cff6 < 0.) out_of_bounds=.true.
        if (cff5 < dt .and. ntstf(itrc) > 1) synchro_flag=.true.
        cff=1./(rho0*Cp*(cff5+cff6))
        cff5=cff5*cff ; cff6=cff6*cff
# endif

# if defined SFLX_CORR && defined SALINITY
        it8=itsss  ;   cff7=sss_time(it8)-tmid
        it7=3-it8  ;   cff8=tmid-sss_time(it7)
        if (cff7 < 0. .or. cff8 < 0.) out_of_bounds=.true.
        if (cff7 < dt .and. ntstf(itrc) > 1) synchro_flag=.true.
        cff=1./(cff7+cff8) ;  cff7=cff7*cff ; cff8=cff8*cff
# endif
        if (itrc == itemp) then
          do j=jstrR,jendR
            do i=istrR,iendR
              stflx(i,j,itemp)=cff1*stflxg(i,j,it1,itemp)
     &                        +cff2*stflxg(i,j,it2,itemp)

# if defined QCORRECTION && !defined ANA_SST
     &         +(cff5*dqdtg(i,j,it5)+cff6*dqdtg(i,j,it6))
     &         *( t(i,j,N,nrhs,itemp) -cff3*sstg(i,j,it3)
     &                                -cff4*sstg(i,j,it4)
     &                                                  )


! Restrict stflx to prevent surface temperature to go below -2
! degrees C.

#  ifdef SEA_ICE_NOFLUX
              stflx(i,j,itemp)=max( stflx(i,j,itemp),
#   ifdef LMD_KPP
     &                          max(hbls(i,j,nstp), Hz(i,j,N))
#   else
     &                                               Hz(i,j,N)
#   endif
     &                         *(-2.D0-t(i,j,N,nrhs,itemp))/dt
     &                                                       )
#  endif
# endif

            enddo
          enddo
# ifdef SALINITY
        elseif (itrc == isalt) then      ! convert fresh-water flux
          do j=jstrR,jendR
            do i=istrR,iendR
              stflx(i,j,isalt)=t(i,j,N,nrhs,isalt)*(
     &                     cff1*stflxg(i,j,it1,isalt)
     &                    +cff2*stflxg(i,j,it2,isalt)
     &                                              )

! Add relaxation of surface salinity back to climatological value to
! avoid long-term drift.  Note that dSSSdt below is "piston velocity"
! expressed in [m/s].

#  if defined SFLX_CORR && defined SALINITY

     &                 -dSSSdt*( t(i,j,N,nrhs,isalt)
c**  &                             -tclm(i,j,N,isalt)

     &         -cff7*sssg(i,j,it7)-cff8*sssg(i,j,it8)
     &                                              )
#  endif
            enddo
          enddo
# endif
        endif  !<-- itemp or isalt

        if (out_of_bounds) then
          if (ZEROTH_TILE) then
            write(*,'(/1x,3A/3(1x,A,F16.10)/)')       '### WARNING: ',
     &  'set_stflux_tile :: Current model time is outside bounds of ',
     &  '''stf_time''.',  'stf_tstart =', stf_time(it1,itrc)*sec2day,
     &                    'tdays =',      tdays,
     &                    'stf_tend =',   stf_time(it2,itrc)*sec2day
          endif
          ierr=ierr+1
        endif
      elseif (iic == ntstart) then        !--> load time-invariant
        if (itrc == itemp) then           !    surface tracer flux
          cff=1./(rho0*Cp)
          do j=jstrR,jendR
            do i=istrR,iendR
              stflx(i,j,itemp)=cff*stflxg(i,j,itstf(itemp),itemp)
            enddo
          enddo
# ifdef SALINITY
        elseif (itrc == isalt) then
          cff=0.01/86400.
          do j=jstrR,jendR
            do i=istrR,iendR
              stflx(i,j,isalt)=cff*stflxg(i,j,itstf(isalt),isalt)
     &                                       *t(i,j,N,nrhs,isalt)
            enddo
          enddo
# endif
        endif
      endif   !<-- ntstf(itrc) > 1

      ! DevinD debug
      if(mynode==0) then
        print *, 'stflx(20,20,itemp)= ', stflx(20,20,itemp)
        print *, 'stflx(20,20,isalt)= ', stflx(20,20,isalt)
      endif

      end
#else
      subroutine set_stflux_empty
      end
#endif /* SOLVE3D  && !ANA_STFLUX || (SALINITY && !ANA_SSFLUX) */
