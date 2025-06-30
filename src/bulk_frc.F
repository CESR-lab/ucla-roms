      module bulk_frc

#include "cppdefs.opt"

      use param
      use dimensions
      use mpi_exchanges
      use roms_read_write

      implicit none

      private  ! make all variable private to this module unless public specified

#include "bulk_frc.opt"

      real,allocatable,dimension(:,:) ::  tair             ! air temp
      real,allocatable,dimension(:,:) ::     Q             ! SPECIFIC humidity (kg/kg) (not relative), qair in old code
      real,allocatable,dimension(:,:) :: prate             ! prate (precipitation rate) - 'rain' variable in old code
      real,allocatable,dimension(:,:) :: lwrad             ! Long wave radiation - get_bulk_rad.F & set_bulk_rad() in old code
      !real:: srflx(GLOBAL_2D_ARRAY)                       ! short wave radiation: srflx stored in surf_flux.F module as used in flux_frc too

#if defined TAU_CORRECTION
      real,allocatable,dimension(:,:) :: taux              ! The TAU_CORRECTION flag is used to active a correction term to
      real,allocatable,dimension(:,:) :: tauy              ! bring bulk force field closer to that of the measured data.
#endif

      real,public,allocatable,dimension(:,:) :: wspd_used  ! wind speed at rho point (also used in BGC routines)

      public init_arrays_bulk_frc
      public set_bulk_frc

      contains

! ----------------------------------------------------------------------
      subroutine init_arrays_bulk_frc  ![
      implicit none
      allocate(  tair(GLOBAL_2D_ARRAY) )
      allocate(     Q(GLOBAL_2D_ARRAY) )
      allocate( prate(GLOBAL_2D_ARRAY) )
      allocate( lwrad(GLOBAL_2D_ARRAY) )

      allocate( wspd_used(GLOBAL_2D_ARRAY) )

!     vdata needs to be of global_2d_array size
      allocate( nc_uwnd%vdata(GLOBAL_2D_ARRAY,2) )
      allocate( nc_vwnd%vdata(GLOBAL_2D_ARRAY,2) )
!
      allocate( nc_tair%vdata(GLOBAL_2D_ARRAY,2) )
      allocate( nc_Q%vdata(GLOBAL_2D_ARRAY,2) )
      allocate( nc_prate%vdata(GLOBAL_2D_ARRAY,2) )
      allocate( nc_lwrad%vdata(GLOBAL_2D_ARRAY,2) )
      allocate( nc_swrad%vdata(GLOBAL_2D_ARRAY,2) )

#if defined TAU_CORRECTION
      allocate( taux(GLOBAL_2D_ARRAY) )
      allocate( tauy(GLOBAL_2D_ARRAY) )

      allocate( nc_taux%vdata(GLOBAL_2D_ARRAY,2) )
      allocate( nc_tauy%vdata(GLOBAL_2D_ARRAY,2) )
#endif

      if (interp_frc==1) then
        nc_uwnd%coarse = 1
        nc_vwnd%coarse = 1
        nc_tair%coarse = 1
        nc_Q%coarse = 1
        nc_prate%coarse = 1
        nc_lwrad%coarse = 1
        nc_swrad%coarse = 1
      endif

      ! Print user options (bulk_frc.opt) to netcdf attributes
      ! Note: to turn bulk forces off, edit cppdefs!
      bulk_frc_opt = ''
      if (interp_frc==1) then
        call store_string_att(bulk_frc_opt, 'Interpolation ON')
      else
        call store_string_att(bulk_frc_opt, 'Interpolation OFF')
      endif

      end subroutine init_arrays_bulk_frc  !]
! ----------------------------------------------------------------------
      subroutine set_bulk_frc(istr,iend,jstr,jend) ![

      use surf_flux, only: srflx,uwnd,vwnd

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! Short/Long radiation:
      ! Thermodynamic forcing: Note that BULK_FLUX requires computing the
      ! short-wave radiation flux first because bulk flux routine performs
      ! the final assembly of everything.

      ! All variables (including wind) should be on rho-points in the surface forcing
      ! file. The ranges here need to match the sizes in the partitioned files
      call set_frc_data(nc_swrad,srflx,'r')
      call set_frc_data(nc_lwrad,lwrad,'r')

      call set_frc_data(nc_uwnd,uwnd,'r')
      call set_frc_data(nc_vwnd,vwnd,'r')

      call set_frc_data(nc_tair,tair,'r')
      call set_frc_data(nc_Q,Q,'r')

      ! Possibly check units of Q and others

      call set_frc_data(nc_prate,prate,'r')

#if defined TAU_CORRECTION
      call set_frc_data(nc_taux,taux,'r')
      call set_frc_data(nc_tauy,tauy,'r')
!     call exchange_xxx(taux,tauy)
#endif

!     call exchange_xxx(srflx,lwrad)
!     call exchange_xxx(uwnd,vwnd)
!     call exchange_xxx(tair,Q)
!     call exchange_xxx(prate)


      ! NEED TO ADD QCORRECTOR CALLS FOR SEA SURFACE TEMP AND/OR SALINITY:

      ! CALCULATE FLUX INTERACTIONS:
      call calc_all_bulk_forces(istr,iend,jstr,jend)

      if (do_check_units) then
        if(mynode==0) print *,'Checking Bulk force units'
        call nc_check_units(frcfile(nc_Q%ifile),nc_Q%vname,'kg/kg')
        call nc_check_units(
     &         frcfile(nc_prate%ifile),nc_prate%vname,'cm/day')
        do_check_units = .false.
      endif


      end subroutine set_bulk_frc  !]
!-----------------------------------------------------------------------
      subroutine calc_all_bulk_forces(istr,iend,jstr,jend)  ![

! Compute bulk parameterization for surface wind stress and net heat
! fluxes. Adapted from COARE code originally written by David Rutgers
! and Frank Bradley. EMINUSP option for equivalent salt fluxes added
!                                          by Paul Goodman (10/2004).
! Modified by Kate Hedstrom for COARE version 3.0 (03/2005).
! Modified by Jim Edson to correct specific hunidities.
! Modified for ROMS_AGRIF by Fred Marin (2008)
! Modified for ROMS_AGRIF by P. Marchesiello (2009)
! Modified and adapted for roms-ucla by P.Damien (12/2019)

! References:

!  Fairall, C. W., E. F. Bradley, D. P. Rogers, J. B. Edson, and
!      G. S. Young, 1996: Bulk parameterization of air-sea fluxes for
!      tropical ocean-global atmosphere Coupled-Ocean Atmosphere
!      Response Experiment, JGR, 101, 3747-3764.

!  Fairall, C. W., E. F. Bradley, J. E. Hare, A. A. Grachev, and
!      J. B. Edson, 2003: Bulk parameterization of air-sea fluxes:
!      Updates and verification for the COARE algorithm,
!      J. Climate, 16, 571-591.

!  Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B. Edson,
!      and G.S. Young, 1996:  Cool-skin and warm-layer effects on sea
!      surface temperature, JGR, 101, 1295-1308.

!  Liu, W. T., K. B. Katsaros, and J. A. Businger, 1979:  Bulk
!      parameterization of the air-sea exchange of heat and water
!      vapor including the molecular constraints at the interface,
!      J. Atmos. Sci, 36, 1722-1735.

!  Taylor, P. K., and M. A. Yelland, 2001: The dependence of sea
!      surface roughness on the height and steepness of the waves.
!      J. Phys. Oceanogr., 31, 572-590.

!  Oost, W. A., G. J. Komen, C. M. J. Jacobs, and C. van Oort, 2002:
!      New evidence for a relation between wind stress and wave age
!      from measurements during ASGAMAGE. Boundary-Layer Meteor.,
!      103, 409-438.

      ! Modules required
      ! ----------------
      use surf_flux
      use tracers, only: t
      use grid
      use mixing
      use ocean_vars
      use scalars

      implicit none

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend

      integer i,j, IterMax,iter
      real a
      real rho0i
      real TseaC,TseaK,Qsea
      real TairC,TairK,rhoAir ! ,Qair
      real VisAir,Hlv ! Q, RH,
      real delW,delT,delQ
      real u10,Zo10,Cd10,Ch10,Ct10,Cd
      real Ct,CC,Ri,Ribcu,Zetu,L10
      real Wstar,Tstar,Qstar
      real ZoW,ZoT,ZoT10,ZoQ,ZoL,L,Rr,Bff
      real Wpsi,Tpsi,Qpsi
      real wspd0,Wgus,charn

      ! DevinD commented out bulk_psiu,bulk_psit as wasn't linking?
      ! same as ncdf_read_mod in read_write.F module. Don't think you
      ! need to declare integer of function name when function lives in
      ! the same module.
!      real bulk_psiu,bulk_psit

      real hfsen,hflat,hflw,upvel,evap

      real blk_Rgas,blk_ZW,blk_ZT,blk_ZQ,blk_beta,blk_Cpa,blk_Zabl
      parameter (blk_Rgas=287.1)
      parameter (blk_ZW=10.0)  ! The height at which the atm. Wind is taken
      parameter (blk_ZT=10.0)  ! The height at which the atm. Temp is taken
      parameter (blk_ZQ=10.0)  ! The height at which the atm. Humi.is taken
      parameter (blk_Zabl=600.0)
      parameter (blk_beta=1.2)
      parameter (blk_Cpa=1004.67)
      real emiss_lw,SigmaSB
      parameter (emiss_lw=0.985)        ! longwave emissivity of ocean,
                                        ! ref. Dickey et al. (1994)
      parameter (SigmaSB=5.6697E-8)     ! Stefan-Boltzmann constant [m^2/K^4]
      real rhow,patm
      parameter (rhow=1000.0)
      parameter (patm=1010.0)
      real eps,r3
      parameter (eps=1.d-20)
      parameter (r3=1.0/3.0)

!!! CFB feedback
      ! wind-stress correction using wind speed:  rho0*sustr + s_tau*Uo
      !   s_tau = cfb_slope * wspd + cfb_offset [N.m^-3.s]
      !  (recommendended and default if BULK_FLUX - needs wspd data)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029)
      parameter (cfb_offset=0.008)
      real    Wspd_min, stau_ref

!!!      real aer, cer
      real aer(PRIVATE_2D_SCRATCH_ARRAY),
     &     cer(PRIVATE_2D_SCRATCH_ARRAY)

      ! DevinD commented out all scalar values of flux variables tair,
      ! prate, rhum, etc because in the reading process of set_bulk_frc_tile
      ! we have already done a temporal interpolation so don't need to do
      ! it here. E.g. we have tair(i,j) already.
      ! However, having extra array tair, etc for each variable will create
      ! a memory overhead, so might be better to revert to old system at
      ! some point of doing time interpolation here, in order to avoid
      ! storing the interpolated time in its own array.

      real wspd ! uwnd, vwnd, - DevinD commented
!      real tair ! DevinD commented
!      real prate ! DevinD commented
!      real rhum ! DevinD commented
      real radlw
!!!      real uwnd,vwnd, wspd,wspd0, szero, ss
!!!      real Tsea,Ssea,evap, radlw
!!!      real rho0i, prate
!!!      real ew,rw,qw,es, rhum, rs,ra,qa
!!!      real tair, vtair,rhoair,a,ce,cs,htevap
!!!      real hflat,hfsen


!  Input variables:
!
!      - radsw, srflx
!      - radlw        downward longwave radiation
!      - prate
!      - sst           t(N,itemp) or sstsk
!      - sss           t(N,isalt)
!      - tair
!      - rhum
!      - wspd
!      - uwnd
!      - vwnd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

      real, parameter ::
!!!     &      cpair=0.240*4186.8, ! specific heat for air in
!!!                                ! Joules/(kg*degC), where 4186.8 is
!!!                                ! to convert from cal/g to Joules/kg
     &      cpi=1./Cp    ! inverse of specific heat for seawater
                         ! [kg*degC/Joule], Cp is defined in "scalars"


      ! DevinD don't need these anymore!
      integer it1,it2, it3,it4, it5,it6, it7,it8, it9,it10
      real tmid, cff1,cff2, cff3,cff4, cff5,cff6, cff7,cff8,
     &                                 cff9,cff10, cff
# ifdef QCORRECTION
#  define SST_DATA
      integer it13,it14, it11,it12
      real cff13,cff14, cff11,cff12
# endif
# ifdef SFLX_CORR
#  define SSS_DATA
      integer it15,it16
      real cff15,cff16
# endif
# define WIND_DATA
# define TAIR_DATA
# define RHUM_DATA
# define PRATE_DATA
# define LWRAD_DATA

! Extended bounds needed since u2rho and rho2u averaging occurs, and
! lmd_kpp has SMOOTH_HBL which requires extended loop ranges for smoothing.
# include "compute_extended_bounds.h"


      ! SET_SRFLUX.F PLACED HERE FOR MODULE VERSION
      ! ===========================================

! Set-up shortwave radiation flux for current tile.

!      implicit none
!      integer istr,iend,jstr,jend, ierr, i,j, it1,it2
!      real tmid, cff, cff1, cff2
!# ifdef DIURNAL_SRFLUX
!      real Ampl, cos_h, dec,cos_d,sin_d, tan_d, phi, csph,snph, h0
!      real, parameter :: year2day=365.25,  day2year=1.D0/year2day
!# elif defined ROBUST_DIURNAL_SRFLUX
!      real Ampl, dtau, tau1,tau2
!# endif
!# include "param.h"
!# include "scalars"
!# include "grid"
!
!# include "compute_extended_bounds.h"
!
!
!
!! Set coefficients for interpolation. Check that for the next time
!! step [when time=time+dt] both weights will still be positive, and
!! if not, set synchro_flag to signal that new data should be read
!! from an appropriate netCDF input file (master thread only).
!! After that either load time-invariant data or interpolate in time.
!! In either case apply factor 1./(rho0*Cp) to convert units from
!! [Watts/m^2, data units] to [degC*m/sec, units of kinematic flux].
!! Complain about error and signal to quit, if interpolation is
!! needed, but not possible.
!
!      if (ntsrf > 1) then                      ! interpolate shortwave
!        tmid=time+0.5*dt                       ! radiation flux in time
!        it2=itsrf ; cff1=srf_time(it2)-tmid
!        it1=3-it2 ; cff2=tmid-srf_time(it1)
!        if (cff1 < dt) synchro_flag=.true.
!        cff=1./(rho0*Cp*(cff1+cff2))
!        cff1=cff1*cff ; cff2=cff2*cff
!
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
            ! DevinD - time interp not needed but convert
            ! convert [Watts/m^2] into kinematic flux [degC*m/s]
!            srflx(i,j)=cff1*swrad_data(i,j,it1)+cff2*swradg(i,j,it2)
            srflx(i,j)= srflx(i,j)/(rho0*Cp)

# if defined DIURNAL_SRFLUX || defined ROBUST_DIURNAL_SRFLUX
#  ifdef BULK_FRC
            srflx(i,j)=srflx(i,j)*Ampl         ! just modulate
#  else
            cff=stflx(i,j,itemp)-srflx(i,j)    ! subtract short-wave
            srflx(i,j)=srflx(i,j)*Ampl         ! radiatin from the net,
            stflx(i,j,itemp)= cff+srflx(i,j)   ! modulate and add back
#  endif
# endif
          enddo
        enddo
!
!        if (cff1 < 0. .or. cff2 < 0.) then
!          if (ZEROTH_TILE) then
!            write(*,'(/1x,3A/3(1x,A,F16.10)/)')       '### WARNING: ',
!     &     'set_srflux_tile :: Current model time is outside bounds ',
!     &     'of ''srf_time''.', 'srf_tstart =', srf_time(it1)*sec2day,
!     &                         'tdays =',      tdays,
!     &                         'srf_tend =',   srf_time(it2)*sec2day
!          endif
!          ierr=ierr+1
!        endif
!      else     !--> ntsrf == 1                ! load time-invariant
!        if (iic == ntstart) then              ! shortwave radiation
!          cff=1/(rho0*Cp)                     ! flux
!          do j=jstrR,jendR
!            do i=istrR,iendR
!              srflx(i,j)=cff*swradg(i,j,itsrf)
!            enddo
!          enddo
!        endif
!      endif
!      end


      ! SUBROUTINE FROM BULK_FLUX_COARE.F STARTS HERE
      ! =============================================


      rho0i=1.0/rho0


      ! All calculations in this loop at rho-points. Loop over extended
      ! bounds.
      do j=jstrR,jendR
        do i=istrR,iendR

!          uwnd=cff1*uwind(i,j,it1)+cff2*uwind(i,j,it2)
!          vwnd=cff1*vwind(i,j,it1)+cff2*vwind(i,j,it2)
          wspd_used(i,j)=sqrt(uwnd(i,j)*uwnd(i,j)+vwnd(i,j)*vwnd(i,j))

          ! convert [Watts/m^2] into kinematic flux [degC*m/s]
          radlw=lwrad(i,j) / (rho0*Cp)

          wspd0=wspd_used(i,j)
          TairC=tair(i,j) ! DevinD
          TairK=TairC+273.16

          TseaC=t(i,j,N,nrhs,itemp)
          TseaK=TseaC+273.16

          ! label(HUMIDITY) - chose between specific or relative humidity
!          RH=rhum(i,j)       ! Specific humidity (if input data is in specific hum.)
!          RH=rhum(i,j)/100.0 ! To get specific humidity (if input data is pdamien's pacific model L217 COAREnew)

!!!          Ssea=t(i,j,N,nrhs,isalt)

!
!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
          hflw=radlw  ! positive downward
     &        - emiss_lw*rho0i*cpi*SigmaSB*TseaK*TseaK*TseaK*TseaK
!
!-----------------------------------------------------------------------
!  Compute specific humidities (kg/kg).
!
!    note that Qair is the saturation specific humidity at Tair
!                 Q is the actual specific humidity
!              Qsea is the saturation specific humidity at Tsea

!          Saturation vapor pressure in mb is first computed and then
!          converted to specific humidity in kg/kg
!
!          The saturation vapor pressure is computed from Teten formula
!          using the approach of Buck (1981):
!
!          Esat(mb) = (1.0007+3.46E-6*patm)*6.1121*
!                  EXP(17.502*TairC(C)/(240.97+TairC(C)))
!
!          The ambient vapor is found from the definition of the
!          Relative humidity:
!
!          RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
!                                       E = RH*Esat     if RH fractional
!
!          The specific humidity is then found using the relationship:
!
!          Q = 0.622 E/(P + (0.622-1)e)
!
!          Q(kg/kg) = 0.62197*(E(mb)/(patm(mb)-0.378*E(mb)))
!
!-----------------------------------------------------------------------
!
! air saturation vapor pressure (mb), using Teten formula.

!          cff=(1.0007+3.46e-6*patm)*6.1121
!     &        *exp(17.502*TairC/(240.97+TairC))

! specific humidity at Saturation, Qair (kg/kg).

!          Qair=0.62197*(cff/(patm-0.378*cff))

! specific humidity, Q (kg/kg).
!
! The following is CROCO's version to automatically establish the units of the
! the input file. For ROMS we don't want the code to do this automatically, we
! want the user to be sure.
!
!          if (RH < 2.) then                          !RH fraction
!            cff=cff*RH                                 !Vapor pres (mb)
!            Q=0.62197*(cff/(patm-0.378*cff))           !Spec hum (kg/kg)
!          else          !RH input was actually specific humidity in g/kg
!            Q=RH/1000.0                                !Spec Hum (kg/kg)
!          endif

!   RH is already fraction of SPECIFIC humidity in kg/kg (DPD original)
!             Q=RH (DPD original)

!   RH input is currently actually specific humidity in g/kg.
!      DevinD: this should be corrected so that input is in kg/kg.
!             Q=RH/1000.0                                !Spec Hum (kg/kg)

! water saturation vapor pressure (mb), using Teten formula.

          cff=(1.0007+3.46e-6*patm)*6.1121*
     &            exp(17.502*TseaC/(240.97+TseaC))

! vapor pressure reduced for salinity (Kraus & Businger, 1994, pp 42).

          cff=cff*0.98

! Qsea (kg/kg) from vapor pressure.

          Qsea=0.62197*(cff/(patm-0.378*cff))
!
!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
          rhoAir=patm*100.0/(blk_Rgas*TairK*
     &                              (1.0+0.61*Q(i,j)))
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
          VisAir=1.326E-5*(1.0+TairC*(6.542E-3+TairC*
     &               (8.301e-6-4.84e-9*TairC)))
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
          Hlv=(2.501-0.00237*TseaC)*1.0e+6
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus=0.5
          delW=SQRT(wspd0*wspd0+Wgus*Wgus)
          delQ=Qsea-Q(i,j)
          delT=TseaC-TairC
!
!  Neutral coefficients.
!
          ZoW=0.0001
          u10=delW*LOG(10.0/ZoW)/LOG(blk_ZW/ZoW)
          Wstar=0.035*u10
          Zo10=0.011*Wstar*Wstar/g+0.11*VisAir/Wstar
          Cd10=(vonKar/LOG(10.0/Zo10))**2
          Ch10=0.00115
          Ct10=Ch10/sqrt(Cd10)
          ZoT10=10.0/exp(vonKar/Ct10)
          Cd=(vonKar/LOG(blk_ZW/Zo10))**2
!
!  Compute Richardson number.
!
          Ct=vonKar/LOG(blk_ZT/ZoT10)  ! T transfer coefficient
          CC=vonKar*Ct/Cd
          Ribcu=-blk_ZW/(blk_Zabl*0.004*blk_beta**3)
          Ri=-g*blk_ZW*(delT+0.61*TairK*delQ)/
     &          (TairK*delW*delW)
          if (Ri.lt.0.0) then
            Zetu=CC*Ri/(1.0+Ri/Ribcu)       ! Unstable
          else
            Zetu=CC*Ri/(1.0+3.0*Ri/CC)      ! Stable
          endif
          L10=blk_ZW/Zetu
          if (Zetu.gt.50.0) then
            IterMax=1
          else
            IterMax=3
          endif
!
!  First guesses for Monin-Obukhov similarity scales.
!
          Wstar=delW*vonKar/(LOG(blk_ZW/Zo10)-
     &                             bulk_psiu(blk_ZW/L10,pi))
          Tstar=-delT*vonKar/(LOG(blk_ZT/ZoT10)-
     &              bulk_psit(blk_ZT/L10,pi))
          Qstar=-delQ*vonKar/(LOG(blk_ZQ/ZoT10)-
     &              bulk_psit(blk_ZQ/L10,pi))
!
!  Modify Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!
          if (delW.gt.18.0) then
            charn=0.018
          elseif ((10.0.lt.delW).and.(delW.le.18.0)) then
            charn=0.011+0.125*(0.018-0.011)*(delW-10.)
          else
            charn=0.011
          endif
!
!  Iterate until convergence. It usually converges within four
!  iterations.
!
          do Iter=1,IterMax
            ZoW=charn*Wstar*Wstar/g+0.11*VisAir/(Wstar+eps)
            Rr=ZoW*Wstar/VisAir
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ=MIN(1.15e-4,5.5e-5/Rr**0.6)
            ZoT=ZoQ
            ZoL=vonKar*g*blk_ZW*
     &             (Tstar*(1.0+0.61*Q(i,j))+0.61*TairK*Qstar)/
     &             (TairK*Wstar*Wstar*(1.0+0.61*Q(i,j))+eps)
            L=blk_ZW/(ZoL+eps)
!
!  Evaluate stability functions at Z/L.
!
            Wpsi=bulk_psiu(ZoL,pi)
            Tpsi=bulk_psit(blk_ZT/L,pi)
            Qpsi=bulk_psit(blk_ZQ/L,pi)
!
!  Compute wind scaling parameters, Wstar.
!
            Wstar=MAX(eps,delW*vonKar/(LOG(blk_ZW/ZoW)-Wpsi))
            Tstar=-delT*vonKar/(LOG(blk_ZT/ZoT)-Tpsi)
            Qstar=-delQ*vonKar/(LOG(blk_ZQ/ZoQ)-Qpsi)
!
!  Compute gustiness in wind speed.
!
            Bff=-g/TairK*Wstar*(Tstar+0.61*TairK*Qstar)
            if (Bff.gt.0.0) then
              Wgus=blk_beta*(Bff*blk_Zabl)**r3
            else
              Wgus=0.2
            endif
            delW=SQRT(wspd0*wspd0+Wgus*Wgus)
          enddo

!
!-----------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
!
!  Compute transfer coefficients for momentum (Cd).
!
          wspd0=SQRT(wspd0*wspd0+Wgus*Wgus)
          Cd=Wstar*Wstar/(wspd0*wspd0+eps)
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
          hfsen=-blk_Cpa*rhoAir*Wstar*Tstar!

!  Compute turbulent latent heat flux (W/m2), Hl.
!
          hflat=-Hlv*rhoAir*Wstar*Qstar
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
          upvel=-1.61*Wstar*Qstar-(1.0+1.61*Q(i,j))*Wstar*Tstar/TairK
          hflat=hflat+rhoAir*Hlv*upvel*Q(i,j)
!
!=======================================================================
!  Convert to ROMS tracer and momentum surface fluxes
!=======================================================================
!
!  Compute kinematic, surface, net heat flux (degC m/s).  Notice that
!  the signs of latent and sensible fluxes are reversed because fluxes
!  calculated from the bulk formulations above are positive out of the
!  ocean. in ROMS, fluxes are positive downward when they warm the ocean.
!
!  For EMINUSP option,  EVAP = LHeat (W/m2) / Hlv (J/kg) = kg/m2/s
!                       PREC = rain = kg/m2/s
!
!  To convert these rates to m/s divide by freshwater density, rhow.
!
!  Note that when the air is undersaturated in water vapor (Q < Qsea)
!  the model will evaporate and LHeat > 0:
!
!                   LHeat positive out of the ocean
!                    evap positive out of the ocean
!
!  Note that if evaporating, the salt flux is positive
!        and if     raining, the salt flux is negative
!
!  Note that fresh water flux is positive out of the ocean and the
!  salt flux (stflx(isalt)) is positive into the ocean. It is converted
!  to (psu m/s) for stflx(isalt) in "set_vbc.F".
!----------------------------------------------------------------------
!
!  Latent and sensible heat fluxes
!
          hflat=-hflat*rho0i*cpi
          hfsen=-hfsen*rho0i*cpi
!
!  Total surface heat flux
!
          stflx(i,j,itemp)=srflx(i,j)+hflw+hflat+hfsen
!
# ifdef SALINITY
!  Total surface salt flux
!
          evap=-cp*hflat/Hlv
          ! Divide prate by 100.*3600.*24. to convert from [cm/day] to [m/s]
          stflx(i,j,isalt)=(evap- ( prate(i,j)/(100.*3600.*24.) ) )*t(i,j,N,nrhs,isalt)
#endif
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
# ifdef SALINITY
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
# endif
# endif
!
!  Save for the calculation of wind stress
!
          aer(i,j)=rhoAir*wspd0*rho0i
          cer(i,j)=Cd

! kinematic surface wind stress [m^2/s^2]
!---------- ------- ---- ------ ---------
! sustr = rhoA/rhoW * Cd * wspd*uwnd
! svstr = rhoA/rhoW * Cd * wspd*vwnd

! BULK_SM_UPDATE removed as it is always needed
          sustr_r(i,j)=aer(i,j)*cer(i,j)*uwnd(i,j)
#  ifdef MASKING
     &                                           *rmask(i,j)
#  endif
          svstr_r(i,j)=aer(i,j)*cer(i,j)*vwnd(i,j)
#  ifdef MASKING
     &                                           *rmask(i,j)
#  endif

!
!  Store fluxes for diagnostics
!  When storing in netCDF file, multiplied
!  by rho0.Cp to get W/m2
!
!          shflx_rsw(i,j)=srflx(i,j)
!          shflx_lat(i,j)=hflat
!          shflx_sen(i,j)=hfsen
!          shflx_rlw(i,j)=hflw
!
!---------------------------------------------------------------
! Flux correction to surface net heat flux.
!---------------------------------------------------------------
!
# ifdef QCORRECTION
          ! Not sure that this still works....
          cff=cff11*dqdtg(i,j,it11)+cff12*dqdtg(i,j,it12)
          stflx(i,j,itemp)=stflx(i,j,itemp) +cff*(
     &           t(i,j,N,nrhs,itemp) -cff13*sstg(i,j,it13)
     &                               -cff14*sstg(i,j,it14)
     &                                                  )
#  ifdef MASKING
     &                                           *rmask(i,j)
#  endif
# endif /* QCORRECTION */

!
!--------------------------------------------------------------
! Flux correction to surface salt flux.
!--------------------------------------------------------------
# ifdef SFLX_CORR
# ifdef SALINITY
          stflx(i,j,isalt)=stflx(i,j,isalt) -dSSSdt*(
     &           t(i,j,N,nrhs,isalt) - sss(i,j)  )
#  ifdef MASKING
     &                                           *rmask(i,j)
#  endif
# endif
# endif

!---------------------------------------------------------------
! Restrict stflx to prevent surface temperature to go below -1.8
! degrees C.
!---------------------------------------------------------------
#ifdef SEA_ICE_NOFLUX
              if( t(i,j,N,nrhs,itemp) .le. -1.8 ) then
                   stflx(i,j,itemp)=0.
#   if defined LMD_KPP
                   srflx(i,j)=0.
#    endif
              endif
#endif

        enddo
      enddo


!---------------------------------------------------------------------
!  Correct kinematic, surface wind stress (m2/s2) due to surface
!  current feedback     -----> CFB
!---------------------------------------------------------------------
!
!         sustr = sustr + s_tau * Uo,  (or equivalently for svstr)
!
! with  s_tau = cfb_slope  * Wspd + cfb_offset  [N.m^-3.s]  CFB_STRESS
!
!          wind wpeed  Wspd < 3 m/s,          CFB_STRESS
!            --> take s_tau mean value of -0.0027 N.m^-3.s
!
!      Caution: sustr is in [m^2/s^2], sms_scale=1/rho0 [m^3/kg]
!      must be used for conversion
!
!---------------------------------------------------------------------
!
!
!!!!!!!!!!!!! This exchange may probably not be necessary !!!!!
!#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!        call exchange_2_tile(istr,iend,jstr,jend,
!     &                       u(START_2D_ARRAY,1,nnew), N,
!     &                       v(START_2D_ARRAY,1,nnew), N)
!# endif
      ! cfb_offset=  0.008
      ! cfb_slope = -0.0029
      ! This seems discontinuous.... should be Wspd_min = 3.7 m/s

      Wspd_min=3.        ! [m/s]
      stau_ref=-0.0027   ! [N.m^-3.s]

      ! maybe move computation of ustar to here
      ! Note loop ranges: iend+1 is needed because lmd_kpp uses SMOOTH_HBL
      ! and smoothing needs iend+1 buffer information for ustar calc.
      ! iend+1 instead of iendR because rho2u sustr_r = u(i) + u(i+1) & since
      ! iendR=iend+2 for interior tiles (as per compute_extended_bounds.h)
      ! if we used i=iendR then u(i+1) = u(iend+3) is not a value.
      do j=jstrR,jendR
        do i=istrR,iend+1

          wspd=sqrt( uwnd(i,j)*uwnd(i,j) + vwnd(i,j)*vwnd(i,j) )

          if (wspd .gt. Wspd_min) then
            cff = cfb_slope*wspd + cfb_offset ! both cfb's are parameters
          else
            cff = stau_ref
          endif

          ! u values averaged to rho points as sustr still at rho-point
          ! Can't do this at rho2u step below because then sustr(i-1) doesn't
          ! include 'u' addition.
          sustr_r(i,j)=sustr_r(i,j) + cff*0.5*(u(i,j,N,nrhs)+u(i+1,j,N,nrhs))/rho0
          !! This rho0 should be worked into the coefficients
          !! Do a mask multiply at reading time for input vars ?

# ifdef TAU_CORRECTION
          ! Correction to bring bulk forcing towards observed fluxes. taux at rho-point.
          sustr_r(i,j)=sustr_r(i,j) + taux(i,j) / rho0
# endif
          ! rmask applies to tau_correction if on, else svstr_r above that.
          ! either way end result is that svstr_r is rmask'ed.

          ! Convert sustr at rho-point to u-point for Eq.System
          ! Need if (i >= istrR+1) because loop over istrR, but we use
          ! index i-1 below which is out of array bounds at tile's western boundary.
          if (i >= istrR+1) then
            sustr(i,j) = ( sustr_r(i-1,j) + sustr_r(i,j) ) / 2
# ifdef MASKING
          !! JM: I don't think that this is actually doing anything useful
     &                    *umask(i,j)
# endif
          end if

        enddo
      enddo

      ! jend+1 instead of jendR because rho2u svstr_r = v(j) + v(j+1) & since
      ! jendR=jend+2 for interior tiles (as per compute_extended_bounds.h)
      ! if we used j=jendR then v(j+1) = u(jend+3) is not a value.
      do j=jstrR,jend+1
        do i=istrR,iendR

          wspd=sqrt( uwnd(i,j)*uwnd(i,j) + vwnd(i,j)*vwnd(i,j) )

          if (wspd .gt. Wspd_min) then
            cff = cfb_slope*wspd + cfb_offset
          else
            cff = stau_ref
          endif

          ! v values averaged to rho points as svstr still at rho-point
          svstr_r(i,j)=svstr_r(i,j) + cff*0.5*(v(i,j,N,nrhs)+v(i,j+1,N,nrhs))/rho0

# ifdef TAU_CORRECTION
          ! Correction to bring bulk forcing towards observed fluxes.
          svstr_r(i,j)=svstr_r(i,j) + tauy(i,j) /rho0
# endif
          ! rmask applies to tau_correction if on, else svstr_r above that.
          ! either way end result is that svstr_r is rmask'ed.

          ! Convert svstr at rho-point to v-point for Eq.System
          ! Need if (j >= jstrR+1) because loop over jstrR, but we use
          ! index j-1 below which is out of array bounds at tile's southern boundary.
          if (j >= jstrR+1) then
            svstr(i,j) = ( svstr_r(i,j-1) + svstr_r(i,j) ) / 2
# ifdef MASKING
     &                    *vmask(i,j)
# endif
          end if
        enddo
      enddo

      end subroutine calc_all_bulk_forces  !]
!-----------------------------------------------------------------------
      function bulk_psiu (ZoL, pi)  ![
!
!=======================================================================
!                                                                      !
!  This function evaluates the stability function for  wind speed      !
!  by matching Kansas  and free convection forms.  The convective      !
!  form follows Fairall et al. (1996) with profile constants from      !
!  Grachev et al. (2000) BLM.  The  stable  form is from Beljaars      !
!  and Holtslag (1991).                                                !
!                                                                      !
!=======================================================================
!
!
!  Function result
!
      real bulk_psiu
!
!  Imported variable declarations.
!
      real ZoL, pi
!
!  Local variable declarations.
!
      real r3
      parameter (r3=1.0/3.0)
      real Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      if (ZoL<0.0) then
        x=(1.0-15.0*ZoL)**0.25
        psik=2.0*LOG(0.5*(1.0+x))+LOG(0.5*(1.0+x*x))-
     &       2.0*ATAN(x)+0.5*pi
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0)
        y=(1.0-10.15*ZoL)**r3
        psic=1.5*LOG(r3*(1.0+y+y*y))-
     &       cff*ATAN((1.0+2.0*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0+cff)
        bulk_psiu=(1.0-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      else
        cff=MIN(50.0,0.35*ZoL)
        bulk_psiu=-((1.0+ZoL)+0.6667*(ZoL-14.28)/
     &            EXP(cff)+8.525)
      endif

      return
      end function bulk_psiu  !]
!-----------------------------------------------------------------------
      function bulk_psit (ZoL, pi)  ![

!=======================================================================
!                                                                      !
!  This function evaluates the  stability function  for moisture and   !
!  heat by matching Kansas and free convection forms. The convective   !
!  form follows Fairall et al. (1996) with  profile  constants  from   !
!  Grachev et al. (2000) BLM.  The stable form is from  Beljaars and   !
!  and Holtslag (1991).                                                !
!
!=======================================================================
!
!  Function result

      real bulk_psit
!
!  Imported variable declarations.
!
      real ZoL, pi
!
!  Local variable declarations.
!
      real r3
      parameter (r3=1.0/3.0)
      real Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      if (ZoL<0.0) then
        x=(1.0-15.0*ZoL)**0.5
        psik=2.0*LOG(0.5*(1.0+x))
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0)
        y=(1.0-34.15*ZoL)**r3
        psic=1.5*LOG(r3*(1.0+y+y*y))-
     &       cff*ATAN((1.0+2.0*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0+cff)
        bulk_psit=(1.0-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      else
        cff=MIN(50.0,0.35*ZoL)
        bulk_psit=-((1.0+2.0*ZoL)**1.5+
     &            0.6667*(ZoL-14.28)/EXP(cff)+8.525)
      endif

      return
      end function bulk_psit  !]
!-----------------------------------------------------------------------

      end module bulk_frc
