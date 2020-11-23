      module flux_frc

      ! Direct flux forcing module

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)
      ! (only refactoring old code's get/set _st/sm/srflux into module)

#include "cppdefs.h"

      ! Modules needed:
      use read_write, only: set_frc_var_tile
      use tracers,    only: t
      use surf_flux ! sustr, svstr, stflx, srflx, sss, sst


      implicit none

      private ! All variables private to module unless public specified

      ! Includes:
      ! param.h needed for GLOBAL_2D_array to work. NT = number tracer from param.h
#include "param.h"
#include "scalars.h"

      ! ****************************************************************
      ! User inputs

      integer, parameter :: nflx = 5 ! Number of flux variables

      ! Interpolate forcing from coarser input grid (=1) or not (=0)
      integer :: interp_frc = 0 ! Factor 2 only for now

      ! sustr - surface u-momentum stress flux (input data in N/m^2)
      character(len=5) :: sustr_name = 'sustr'
      character(len=8) :: sustr_time_name  = 'sms_time'
      integer, parameter :: sustr_grd_type = 1 ! u-point
      integer, parameter :: isustr = 1 ! Needed for its netcdf variables
                                       ! Sequential for 'nflx'
      ! svstr - surface v-momentum stress flux (input data in N/m^2)
      character(len=5) :: svstr_name = 'svstr'
      character(len=8) :: svstr_time_name  = 'sms_time'
      integer, parameter :: svstr_grd_type = 2 ! v-point
      integer, parameter :: isvstr = 2
      ! stflx(itemp) - surface heat flux (input data in W/m^2)
      character(len=6) :: shflx_name = 'shflux'
      character(len=8) :: shflx_time_name  = 'shf_time'
      integer, parameter :: ishflx = 3 ! Could swap to beginning in order to use itemp=1 & isalt=2?
      ! stflx(isalt) - surface freshwater flux (input data in cm/day)
      ! Might want to use #if def SALINITY?
      character(len=6) :: swflx_name = 'swflux'
      character(len=8) :: swflx_time_name  = 'swf_time'
      integer, parameter :: iswflx = 4
      ! swrad - surface short-wave radiation flux (input data in W/m^2)
      character(len=5) :: swrad_name = 'swrad'
      character(len=8) :: swrad_time_name  = 'srf_time'
      integer, parameter :: iswrad = 5

      ! End of user inputs
      ! ****************************************************************

      !-- Flux netcdf variables as arrays/matrices of 'nflx' length
      ! Abstracted from user.
      real                     :: flx_data(GLOBAL_2D_ARRAY,2,nflx) ! Could set to zero here?
      integer, dimension(nflx) :: flx_grd_type  = 0
      integer, dimension(nflx) :: flx_file_indx = 0
      integer, dimension(nflx) :: flx_irec      = 0
      integer, dimension(nflx) :: flx_it1 = 1, flx_it2 = 2
      real, dimension(2,nflx)  :: flx_times = -99 ! [-99,-99]


      public set_flux_frc

      contains
! ----------------------------------------------------------------------

      subroutine set_flux_frc(istr,iend,jstr,jend)

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j

! scalars.h for rho0
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


      end subroutine set_flux_frc
! ----------------------------------------------------------------------

      subroutine set_smflux(istr,iend,jstr,jend)
      ! set surface momentum flux, sustr, svstr

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j
#include "compute_extended_bounds.h"

      ! sustr
      call set_frc_var_tile(
     &  istr,iend,              jstr,jend,            ! Tile bounds indices
     &  sustr_name,             sustr_time_name,      ! Text name of array and time
     &  flx_data(:,:,:,isustr), sustr,                ! Variable arrays
     &  flx_times(:,isustr),    sustr_grd_type,       ! Input times & Grid type
     &  flx_file_indx(isustr),  flx_irec(isustr),     ! File indx & Current input record
     &  flx_it1(isustr),        flx_it2(isustr),      ! Time index placeholders
     &  interp_frc )                                  ! Online spatial interpolation flag

      ! svstr
      call set_frc_var_tile(
     &  istr,iend,              jstr,jend,            ! Tile bounds indices
     &  svstr_name,             svstr_time_name,      ! Text name of array and time
     &  flx_data(:,:,:,isvstr), svstr,                ! Variable arrays
     &  flx_times(:,isvstr),    svstr_grd_type,       ! Input times & Grid type
     &  flx_file_indx(isvstr),  flx_irec(isvstr),     ! File indx & Current input record
     &  flx_it1(isvstr),        flx_it2(isvstr),      ! Time index placeholders
     &  interp_frc )                                  ! Online spatial interpolation flag

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


      end subroutine set_smflux
! ----------------------------------------------------------------------

      subroutine set_shflux(istr,iend,jstr,jend)
      ! Set surface heat flux

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j
      real rho0cp
#include "compute_extended_bounds.h"

      ! stflx(itemp)
      call set_frc_var_tile(
     &  istr,iend,              jstr,jend,            ! Tile bounds indices
     &  shflx_name,             shflx_time_name,      ! Text name of array and time
     &  flx_data(:,:,:,ishflx), stflx(:,:,itemp),     ! Variable arrays
     &  flx_times(:,ishflx),    flx_grd_type(ishflx), ! Input times & Grid type
     &  flx_file_indx(ishflx),  flx_irec(ishflx),     ! File indx & Current input record
     &  flx_it1(ishflx),        flx_it2(ishflx),      ! Time index placeholders
     &  interp_frc )                                  ! Online spatial interpolation flag

      rho0cp=rho0*Cp ! Set coefficient outside of loop for efficiency
      do j=jstrR,jendR
        do i=istrR,iendR
          stflx(i,j,itemp)=stflx(i,j,itemp)/rho0cp

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

      ! DevinD debug
!      if(mynode==0) then
!        print *, 'stflx(20,20,itemp)= ', stflx(20,20,itemp)
!      endif

      end subroutine set_shflux
! ----------------------------------------------------------------------
#ifdef SALINITY
      subroutine set_swflux(istr,iend,jstr,jend)
      ! set surface freshwater flux: stflx(isalt)

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer i, j
      real cmday2ms
# include "compute_extended_bounds.h"

      call set_frc_var_tile(
     &  istr,iend,              jstr,jend,            ! Tile bounds indices
     &  swflx_name,             swflx_time_name,      ! Text name of array and time
     &  flx_data(:,:,:,iswflx), stflx(:,:,isalt),     ! Variable arrays
     &  flx_times(:,iswflx),    flx_grd_type(iswflx), ! Input times & Grid type
     &  flx_file_indx(iswflx),  flx_irec(iswflx),     ! File indx & Current input record
     &  flx_it1(iswflx),        flx_it2(iswflx),      ! Time index placeholders
     &  interp_frc )                                  ! Online spatial interpolation flag

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

      end subroutine set_swflux
#endif /* SALINITY */
! ----------------------------------------------------------------------

      subroutine set_srflux(istr,iend,jstr,jend)
      ! set short-wave radiation flux

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
!# include "param.h" ! Have in module preamble
!# include "scalars.h"! Have in module preamble
# include "grid.h"

#include "compute_extended_bounds.h"

      ! Do initial time interpolation only of swrad:
      call set_frc_var_tile(
     &  istr,iend,              jstr,jend,            ! Tile bounds indices
     &  swrad_name,             swrad_time_name,      ! Text name of array and time
     &  flx_data(:,:,:,iswrad), srflx,                ! Variable arrays
     &  flx_times(:,iswrad),    flx_grd_type(iswrad), ! Input times & Grid type
     &  flx_file_indx(iswrad),  flx_irec(iswrad),     ! File indx & Current input record
     &  flx_it1(iswrad),        flx_it2(iswrad),      ! Time index placeholders
     &  interp_frc )                                  ! Online spatial interpolation flag

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

      end subroutine set_srflux
! ----------------------------------------------------------------------

      end module flux_frc

