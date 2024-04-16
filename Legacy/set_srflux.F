#include "cppdefs.h"
! DevinD moved set_srflux to bulk_frc.F module when BULK_FLUX on. Avoid variable clash
#if defined SOLVE3D && !defined ANA_SRFLUX && !defined BULK_FLUX && !defined FLUX_FRC

c--#define DIURNAL_SRFLUX
c--#define ROBUST_DIURNAL_SRFLUX


      subroutine set_srflux_tile (istr,iend,jstr,jend, ierr)

! Set-up shortwave radiation flux for current tile.

      implicit none
      integer istr,iend,jstr,jend, ierr, i,j, it1,it2
      real tmid, cff, cff1, cff2
# ifdef DIURNAL_SRFLUX
      real Ampl, cos_h, dec,cos_d,sin_d, tan_d, phi, csph,snph, h0
      real, parameter :: year2day=365.25,  day2year=1.D0/year2day
# elif defined ROBUST_DIURNAL_SRFLUX
      real Ampl, dtau, tau1,tau2
# endif
# define SWRAD_DATA
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "forces.h"

# include "compute_extended_bounds.h"


# ifdef ROBUST_DIURNAL_SRFLUX
      dtau=2.D0*pi * dt*sec2day
      cff=4.D0*pi*(tdays-dble(int(tdays))-0.5D0)
      tau1=max(-pi, min(pi, cff))
      tau2=max(-pi, min(pi, cff+2.D0*dtau))
      cff=tau2-tau1
      Ampl=( cff + 2.D0*cos(0.5D0*(tau1+tau2))*sin(0.5D0*cff) )/dtau
# endif


! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data or interpolate in time.
! In either case apply factor 1./(rho0*Cp) to convert units from
! [Watts/m^2, data units] to [degC*m/sec, units of kinematic flux].
! Complain about error and signal to quit, if interpolation is
! needed, but not possible.

      if (ntsrf > 1) then                      ! interpolate shortwave
        tmid=time+0.5*dt                       ! radiation flux in time
        it2=itsrf ; cff1=srf_time(it2)-tmid
        it1=3-it2 ; cff2=tmid-srf_time(it1)
        if (cff1 < dt) synchro_flag=.true.
        cff=1./(rho0*Cp*(cff1+cff2))
        cff1=cff1*cff ; cff2=cff2*cff

# ifdef DIURNAL_SRFLUX

! DIURNAL CYCLE - USED IN BOTH PHYSICAL AND ECOSYSTEM MODELS
! Patrick Marchesiello - 1999: Modulate average dayly insolation
! to get diurnal cycle by:
!
!              cos(h)*cos(d)*cos(phi)  +  sin(d)*sin(phi)
!       pi * ---------------------------------------------
!             sin(h0)*cos(d)*cos(phi) + h0*sin(d)*sin(phi)
!
! where: h, d, phi -- are hour, declination, latitude angles;
!        h0 is hour angle at sunset and sunrise
!
! Yusuke Uchiyama, 2009: UTC correction based on lonr is added.
!                               ocean_time should be in UTC.
#  define UTC_CORRECTION

#  ifndef UTC_CORRECTION
        cos_h=cos(2.*pi*(tdays+0.5-int(tdays+0.5)))
#  endif
        dec=-0.406*cos(deg2rad*(tdays-int(tdays*day2year)*year2day))
        cos_d=cos(dec) ; sin_d=sin(dec) ; tan_d=tan(dec)
# endif

        do j=jstrR,jendR
          do i=istrR,iendR
# ifdef DIURNAL_SRFLUX
#  ifdef UTC_CORRECTION
            cos_h=cos( 2.*pi*(tdays+0.5 -int(tdays+0.5))
     &                               +deg2rad*lonr(i,j) )
#  endif
            phi=deg2rad*latr(i,j)
            h0=acos(-tan(phi)*tan_d)
            csph=cos_d*cos(phi) ; snph=sin_d*sin(phi)

            Ampl=max( 0.,  pi*(cos_h*csph +  snph)
     &                     /(sin(h0)*csph + h0*snph)
     &                                             )
# endif

            srflx(i,j)=cff1*swradg(i,j,it1)+cff2*swradg(i,j,it2)

# if defined DIURNAL_SRFLUX || defined ROBUST_DIURNAL_SRFLUX
#  ifdef BULK_FLUX
            srflx(i,j)=srflx(i,j)*Ampl         ! just modulate
#  else
            cff=stflx(i,j,itemp)-srflx(i,j)    ! subtract short-wave
            srflx(i,j)=srflx(i,j)*Ampl         ! radiatin from the net,
            stflx(i,j,itemp)= cff+srflx(i,j)   ! modulate and add back
#  endif
# endif
          enddo
        enddo

        if (cff1 < 0. .or. cff2 < 0.) then
          if (ZEROTH_TILE) then
            write(*,'(/1x,3A/3(1x,A,F16.10)/)')       '### WARNING: ',
     &     'set_srflux_tile :: Current model time is outside bounds ',
     &     'of ''srf_time''.', 'srf_tstart =', srf_time(it1)*sec2day,
     &                         'tdays =',      tdays,
     &                         'srf_tend =',   srf_time(it2)*sec2day
          endif
          ierr=ierr+1
        endif
      else     !--> ntsrf == 1                ! load time-invariant
        if (iic == ntstart) then              ! shortwave radiation
          cff=1/(rho0*Cp)                     ! flux
          do j=jstrR,jendR
            do i=istrR,iendR
              srflx(i,j)=cff*swradg(i,j,itsrf)
            enddo
          enddo
        endif
      endif

      ! DevinD debug
      if(mynode==0) then
        print *, 'srflx(20,20)= ', srflx(20,20)
      endif

      end
#else
      subroutine set_srflux_empty
      end
#endif /* SOLVE3D && !ANA_SRFLUX */

