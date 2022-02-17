      module flux_frc

      ! Direct flux forcing module

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)
      ! (only refactoring old code's get/set _st/sm/srflux into module)

#include "cppdefs.opt"

      ! Modules needed:
      use param
      use read_write, only: set_frc_var_tile, ncvs
      use tracers,    only: t
      use surf_flux ! sustr, svstr, stflx, srflx, sss, sst
      use scalars


      implicit none

      private ! All variables private to module unless public specified

      ! Includes:
#include "flux_frc.opt"


      public set_flux_frc
      public init_arrays_flux_frc

      contains

! ----------------------------------------------------------------------
      subroutine init_arrays_flux_frc  ![
      implicit none

      allocate( nc_sustr%data( GLOBAL_2D_ARRAY,2) )
      allocate( nc_svstr%data( GLOBAL_2D_ARRAY,2) )
      allocate( nc_shflx%data( GLOBAL_2D_ARRAY,2) )
      allocate( nc_swflux%data(GLOBAL_2D_ARRAY,2) )
      allocate( nc_swrad%data( GLOBAL_2D_ARRAY,2) )

      end subroutine init_arrays_flux_frc  !]

! ----------------------------------------------------------------------
      subroutine set_flux_frc(istr,iend,jstr,jend)  ![

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! local
      integer :: i, j

! scalars for rho0
! index bounds for do loop division by rho0
#include "compute_extended_bounds.h"

      ! 1) set surface momentum flux
      call set_smflux(istr,iend,jstr,jend)

      ! 2) set surface heat flux: stflx(itemp)
      call set_shflux(istr,iend,jstr,jend)

      ! 3) set short-wave radiation flux
      call set_srflux(istr,iend,jstr,jend)

      ! 4) set water flux: stflx(isalt)
#ifdef SALINITY
      call set_swflux(istr,iend,jstr,jend)
#endif


      end subroutine set_flux_frc  !]

! ----------------------------------------------------------------------
      subroutine set_smflux(istr,iend,jstr,jend)  ![
      ! set surface momentum flux, sustr, svstr

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j
#include "compute_extended_bounds.h"

      call set_frc_var_tile( istr, iend, jstr, jend, nc_sustr, sustr, interp_frc )
      call set_frc_var_tile( istr, iend, jstr, jend, nc_svstr, svstr, interp_frc )

      ! Combine this loop into one loop at the end
      do j=jstrR,jendR
        do i=istrR,iendR
          sustr(i,j)=sustr(i,j) / rho0
          svstr(i,j)=svstr(i,j) / rho0
        enddo
      enddo
      ! Move rho0 into arguement of set_var_frc_tile eventually?

      ! DevinD debug
!      if(mynode==0) then
!        print *, 'sustr(20,20)= ', sustr(20,20)
!        print *, 'svstr(20,20)= ', svstr(20,20)
!      endif

      end subroutine set_smflux  !]

! ----------------------------------------------------------------------
      subroutine set_shflux(istr,iend,jstr,jend)  ![
      ! Set surface heat flux

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j
      real rho0cp
#include "compute_extended_bounds.h"


      call set_frc_var_tile( istr, iend, jstr, jend, nc_shflx, stflx(:,:,itemp), interp_frc )
!     call set_frc_var_tile( istr, iend, jstr, jend, nc_shflx, btflx(:,:,itemp), interp_frc )


      rho0cp=rho0*Cp ! Set coefficient outside of loop for efficiency
      do j=jstrR,jendR
        do i=istrR,iendR
          stflx(i,j,itemp)=stflx(i,j,itemp)/rho0cp
!         btflx(i,j,itemp)=stflx(i,j,itemp)/rho0cp

# if defined QCORRECTION && !defined ANA_SST
     &         -dSSTdt*(t(i,j,N,nrhs,itemp) - sst(i,j))

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

      end subroutine set_shflux  !]

! ----------------------------------------------------------------------
#ifdef SALINITY
      subroutine set_swflux(istr,iend,jstr,jend)  ![
      ! set surface freshwater flux: stflx(isalt)

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j
      real cmday2ms
# include "compute_extended_bounds.h"

      call set_frc_var_tile( istr, iend, jstr, jend, nc_swflux, stflx(:,:,isalt), interp_frc )

      cmday2ms=0.01/day2sec ! Set coefficient outside of loop for efficiency
      do j=jstrR,jendR
        do i=istrR,iendR
          ! cm/day -> m/s
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N,nrhs,isalt)*cmday2ms

! Add relaxation of surface salinity back to climatological value to
! avoid long-term drift.  Note that dSSSdt below is "piston velocity"
! expressed in [m/s].
# if defined SFLX_CORR
     &                 -dSSSdt*( t(i,j,N,nrhs,isalt)-sss(i,j) )
# endif

        enddo
      enddo

      ! DevinD debug
!      if(mynode==0) then
!        print *, 'stflx(20,20,isalt)= ', stflx(20,20,isalt)
!      endif

      end subroutine set_swflux  !]
#endif /* SALINITY */

! ----------------------------------------------------------------------
      subroutine set_srflux(istr,iend,jstr,jend)  ![
      ! set short-wave radiation flux

      use grid

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j, it1, it2
      real tmid, cff, cff1, cff2
# ifdef DIURNAL_SRFLUX
      real Ampl, cos_h, dec,cos_d,sin_d, tan_d, phi, csph,snph, h0
      real, parameter :: year2day=365.25,  day2year=1.D0/year2day
# elif defined ROBUST_DIURNAL_SRFLUX
      real Ampl, dtau, tau1,tau2
# endif
# define SWRAD_DATA

#include "compute_extended_bounds.h"

      ! Do initial time interpolation only of swrad:
      call set_frc_var_tile( istr, iend, jstr, jend, nc_swrad, srflx, interp_frc )

      ! Note: the additional CPPFLAG option have been left in from old
      ! code but are untested here, might not compile.

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

                                               ! interpolate shortwave
!        tmid=time+0.5*dt                       ! radiation flux in time
!        it2=itsrf ; cff1=srf_time(it2)-tmid
!        it1=3-it2 ; cff2=tmid-srf_time(it1)
!        if (cff1 < dt) synchro_flag=.true.
!        cff=1./(rho0*Cp*(cff1+cff2))
!        cff1=cff1*cff ; cff2=cff2*cff

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

        ! This loop is not subject to any CPP flags
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

!            srflx(i,j)=cff1*swradg(i,j,it1)+cff2*swradg(i,j,it2)
            srflx(i,j)=srflx(i,j)/(rho0*Cp) ! DevinD equivalent

# if defined DIURNAL_SRFLUX || defined ROBUST_DIURNAL_SRFLUX
            cff=stflx(i,j,itemp)-srflx(i,j)    ! subtract short-wave
            srflx(i,j)=srflx(i,j)*Ampl         ! radiatin from the net,
            stflx(i,j,itemp)= cff+srflx(i,j)   ! modulate and add back
# endif
          enddo
        enddo

      end subroutine set_srflux  !]

! ----------------------------------------------------------------------

      end module flux_frc

