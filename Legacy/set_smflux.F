#include "cppdefs.h"
! DevinD added BULK_FLUX as get wind accounted for in flux modules
! when BULK_FLUX is turned on
! DevinD flux_frc module. Remove this file completely once working
#if !defined ANA_SMFLUX && !defined BULK_FLUX && !defined FLUX_FRC

      subroutine set_smflux_tile(istr,iend,jstr,jend, ierr)
      implicit none
      integer istr,iend,jstr,jend, ierr, it1,it2, i,j
      real tmid, cff1,cff2, cff
# ifdef WIND_MAGN
      real cff3,cff4
      integer it3,it4
# endif
# define WIND_DATA
# include "param.h"
# ifdef CURVGRID
#  include "grid.h"
# endif
# include "forces.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible. Note that multiplication by 1/rho0 is to
! convert from [Newton/m^2] to [m^2/s^2], i.e to make it kinematic
! momentum flux.

      tmid=time+0.5*dt                 !--> interpolate wind in time
      if (ntwnd > 1) then
        it2=itwnd ; cff1=wnd_time(it2)-tmid
        it1=3-it2 ; cff2=tmid-wnd_time(it1)
        if (cff1 < dt) synchro_flag=.true.
        cff=1./(rho0*(cff1+cff2))
        cff1=cff1*cff ; cff2=cff2*cff

        do j=jstrR,jendR
          do i=istrR,iendR
            sustr(i,j)=cff1*uwind(i,j,it1)+cff2*uwind(i,j,it2)
            svstr(i,j)=cff1*vwind(i,j,it1)+cff2*vwind(i,j,it2)
          enddo
        enddo

        if (cff1 < 0. .or. cff2 < 0.) then
          if (ZEROTH_TILE) then
            write(*,'(/1x,3A/3(1x,A,F16.10)/)')      '### WARNING: ',
     &    'set_smflux_tile :: Current model time is outside bounds ',
     &    'of ''wnd_time''.', 'sms_tstart =', wnd_time(it1)*sec2day,
     &                        'tdays =',      tdays,
     &                        'sms_tend =',   wnd_time(it2)*sec2day
          endif
          ierr=ierr+1
        endif
      elseif (iic == ntstart) then           ! load time-invariant
        cff=1./rho0                          ! wind at initialization
        do j=jstrR,jendR                     ! only when there is only 
          do i=IstR,iendR                    ! one record, ntwnd == 1 
            sustr(i,j)=cff*uwind(i,j,itwnd)  
            svstr(i,j)=cff*vwind(i,j,itwnd)
          enddo
        enddo
      endif

# ifdef WIND_MAGN
      if (ntwmag > 1) then
        it4=itwmag ; cff3=wmag_time(it4)-tmid
        it3=3-it4  ; cff4=tmid-wmag_time(it3)
        if (cff3 < dt) synchro_flag=.true.
        cff=1./(rho0*(cff3+cff4))
        cff3=cff3*cff ; cff4=cff4*cff

        do j=jstrR,jendR
          do i=istrR,iendR
            wndmag(i,j)=cff3*windmag(i,j,it3)+cff4*windmag(i,j,it4)
          enddo
        enddo

        if (cff3 < 0. .or. cff4 < 0.) then
          if (ZEROTH_TILE) then
            write(*,'(/1x,3A/3(1x,A,F16.10)/)')      '### WARNING: ',
     &    'set_smflux_tile :: Current model time is outside bounds ',
     &    'of ''wmag_time''.','atau_tstart =',wmag_time(it3)*sec2day,
     &                        'tdays =',      tdays,
     &                        'atau_tend =',  wmag_time(it4)*sec2day
          endif
          ierr=ierr+1
        endif
      elseif (iic == ntstart) then           ! load time-invariant
        cff=1./rho0                          ! wind stress magnitude
        do j=jstrR,jendR                     ! at initialization only
          do i=IstR,iendR                    ! when  ntwmag == 1
            wndmag(i,j)=cff*windmag(i,j,itwmag)
          enddo
        enddo
      endif
# endif

      ! DevinD debug
      if(mynode==0) then
        print *, 'sustr(20,20)= ', sustr(20,20)
        print *, 'svstr(20,20)= ', svstr(20,20)
      endif

      end
#else
      subroutine get_smflux_empty
      end
#endif /* !ANA_SMFLUX */

