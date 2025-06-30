      module flux_frc

      ! Direct flux forcing module

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)
      ! (only refactoring old code's get/set _st/sm/srflux into module)

#include "cppdefs.opt"

      ! Modules needed:
      use param
      use dimensions
      use tracers,    only: t
      use surf_flux ! sustr, svstr, stflx, srflx, sss, sst
      use scalars
      use roms_read_write
      use grid


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

      allocate( nc_sustr%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_sustr%vdata = 0
      allocate( nc_svstr%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_svstr%vdata = 0
      allocate( nc_shflx%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_shflx%vdata = 0
      allocate( nc_swflux%vdata(GLOBAL_2D_ARRAY,2) );nc_swflux%vdata = 0
      allocate( nc_swrad%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_swrad%vdata = 0

      if (interp_frc==1) then
        nc_sustr%coarse=1
        nc_svstr%coarse=1
        nc_shflx%coarse=1
        nc_swflux%coarse=1
        nc_swrad%coarse=1
      endif

      ! Print user options (flux_frc.opt) to netcdf attributes
      ! Note: to turn flux forces off, edit cppdefs!
      flux_frc_opt = ''
      if (interp_frc==1) then
        call store_string_att(flux_frc_opt, 'Interpolation ON')
      else
        call store_string_att(flux_frc_opt, 'Interpolation OFF')
      endif
      flux_frc_opt = trim(adjustl(flux_frc_opt))
      end subroutine init_arrays_flux_frc  !]
! ----------------------------------------------------------------------
      subroutine set_flux_frc(istr,iend,jstr,jend)  ![

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! local
      integer :: i, j,ierr

      ! 1) set surface momentum flux
      call set_frc_data(nc_sustr,sustr,'u')
      call set_frc_data(nc_svstr,svstr,'v')
      sustr = sustr/rho0
      svstr = svstr/rho0

      ! 2) set surface heat flux: stflx(itemp)
      call set_frc_data(nc_shflx,stflx(:,:,itemp),'r')
      stflx(:,:,itemp) = stflx(:,:,itemp) /(rho0*Cp)
#ifdef SEA_ICE_NOFLUX
      do j=j0,j1
        do i=i0,i1
          if( t(i,j,nz,nrhs,itemp) .le. -1.8 ) then
            stflx(i,j,itemp)=0.
#   if defined LMD_KPP
            srflx(i,j)=0.
#    endif
          endif
        enddo
      enddo
#endif
# if defined QCORRECTION
      do j=j0,j1
        do i=i0,i1
          stflx(i,j,itemp)=
     &       -dSSTdt*(t(i,j,nz,nrhs,itemp) - sst(i,j))
        enddo
      enddo
#endif

      ! 2a) set bottom heat flux: stflx(itemp)
!     call set_frc_data(nc_bhflx,btflx(i0:i1,j0:j1,itemp),'r')

      ! 3) set short-wave radiation flux
      call set_srflux

#ifdef SALINITY
      ! 4) set water flux: stflx(isalt)
      call set_swflux
#endif

#ifdef EXCHANGE
!     call exchange_xxx(sustr,svstr)
!     call exchange_xxx(srflx,stflx(:,:,itemp))
#ifdef SALINITY
!     call exchange_xxx(stflx(:,:,isalt))
#endif
#endif
!     srflx = 0
!     print *,'maxval srflx: ', maxval(abs(srflx))

      end subroutine set_flux_frc  !]
! ----------------------------------------------------------------------
      subroutine set_swflux  ![
      ! set surface freshwater flux: stflx(isalt)

      implicit none

      ! local
      integer i,j
      real cmday2ms

      call set_frc_data(nc_swflux,stflx(:,:,isalt),'r')

      cmday2ms=0.01/day2sec ! Set coefficient outside of loop for efficiency
      do j=j0,j1
        do i=i0,i1
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

      end subroutine set_swflux  !]
! ----------------------------------------------------------------------
      subroutine set_srflux  ![
      ! set short-wave radiation flux

      implicit none

      ! local
      integer i, j, it1, it2
      real tmid, cff, cff1, cff2
# ifdef DIURNAL_SRFLUX
      real Ampl, cos_h, dec,cos_d,sin_d, tan_d, phi, csph,snph, h0
      real, parameter :: year2day=365.25,  day2year=1.D0/year2day
# endif

      call set_frc_data(nc_swrad,srflx,'r')
      srflx = srflx/(rho0*Cp)

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

        dec=-0.406*cos(deg2rad*(tdays-int(tdays*day2year)*year2day))
        cos_d=cos(dec) ; sin_d=sin(dec) ; tan_d=tan(dec)

        do j=j0:j1
          do i=i0:i1
            cos_h=cos( 2.*pi*(tdays+0.5 -int(tdays+0.5))
     &                               +deg2rad*lonr(i,j) )
            phi=deg2rad*latr(i,j)
            h0=acos(-tan(phi)*tan_d)
            csph=cos_d*cos(phi) ; snph=sin_d*sin(phi)

            Ampl=max( 0.,  pi*(cos_h*csph +  snph)
     &                     /(sin(h0)*csph + h0*snph)
     &                                             )

            cff=stflx(i,j,itemp)-srflx(i,j)    ! subtract short-wave
            srflx(i,j)=srflx(i,j)*Ampl         ! radiating from the net,
            stflx(i,j,itemp)= cff+srflx(i,j)   ! modulate and add back
          enddo
        enddo
# endif

      end subroutine set_srflux  !]
! ----------------------------------------------------------------------

      end module flux_frc

