      module flux_frc

      ! Direct flux forcing module

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)
      ! (only refactoring old code's get/set _st/sm/srflux into module)

#include "cppdefs.opt"

      ! Modules needed:
      use dimensions
      use tracers
      use surf_flux
      use roms_read_write
      use grid
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
      subroutine set_flux_frc  ![

      implicit none

      ! local
      integer :: i,j

      ! Surface momentum flux
      call set_frc_data(nc_sustr,sustr,'u')
      call set_frc_data(nc_svstr,svstr,'v')
      sustr = sustr/rho0
      svstr = svstr/rho0

      ! Surface heat flux
      call set_frc_data(nc_shflx,stflx(:,:,itemp),'r')
      stflx(:,:,itemp) = stflx(:,:,itemp) /(rho0*Cp)

      ! Surface salinity flux
      call set_frc_data(nc_swflux,stflx(:,:,isalt),'r')
!     stflx(:,:,isalt) = stflx(:,:,isalt)*t(:,:,nz,nrhs,isalt)*cmday2ms
      swflx = -stflx(:,:,isalt)*cmday2ms
      stflx(:,:,isalt) = 0

      ! Incoming net shortwave radiation
      call set_frc_data(nc_swrad,srflx,'r')
      srflx = srflx/(rho0*Cp)

      ! Bottom heat flux: stflx(itemp)
!     call set_frc_data(nc_bhflx,btflx(:,:,itemp),'r')

      ! Modulation of surface fluxes

      ! Sea ice blocks heat fluxs and incoming light
      if (sea_ice_noflux) then
        do j=0,ny+1
          do i=0,nx+1
            if( t(i,j,nz,nrhs,itemp) .le. -1.8 ) then
              stflx(i,j,itemp) = 0.
              srflx(i,j) = 0.
            endif
          enddo
        enddo
      endif

      ! SST restoring
      if (qcorrection) then
        stflx(:,:,itemp) = stflx(:,:,itemp)
     &        -dSSSdt*(t(:,:,nz,nrhs,isalt)-sss(:,:))
      endif

      ! SST restoring
      if (sflx_correction) then
        stflx(:,:,isalt) = stflx(:,:,isalt)
     &        -dSSSdt*(t(:,:,nz,nrhs,isalt)-sss(:,:))
      endif

      ! Diurnal modulation of averaged short wave radiation flux
      if (diurnal_srflux) then
        call diurnal_mod_srflx
      endif

      end subroutine set_flux_frc  !]
! ----------------------------------------------------------------------
      subroutine diurnal_mod_srflx  ![
      ! Daily modulation of short-wave radiation flux

      implicit none

      ! local
      integer ::  i,j
      real    :: cff

      real Ampl, cos_h, dec,cos_d,sin_d, tan_d, phi, csph,snph, h0
      real,parameter :: year2day=365.25
      real,parameter :: day2year=1.0/year2day

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

        do j=j0,j1
          do i=i0,i1
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

      end subroutine diurnal_mod_srflx  !]
! ----------------------------------------------------------------------

      end module flux_frc

