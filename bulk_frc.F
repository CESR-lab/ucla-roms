      module bulk_frc

      ! ===================
      ! BULK FORCING MODULE
      ! ===================

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.h"
! CPP flag BULK_FRC_MOD covering entire module so that the module is empty if
! BULK_FRC_MOD is unused, and thus module will not take up unnecessary size in
! the executable.
#if defined BULK_FLUX /* for whole module! */

      ! Modules required
      ! ================

      use read_write, only: set_frc_var_tile

      implicit none

      ! Make all variable private to this module unless public specified
      private

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_array to work
#include "param.h"

      ! Variables
      ! =========

      ! wind
      ! ----

      ! Note in old code u-wind was uwnd but in new code it is
      ! uwind. Should be consistent as some point. Same for v-wind.

       ! u-wind (input data in m/s)
       ! ------
        ! stores raw data u-wind for 2 times: time_A < model_time < time_B
      real uwnd_data(GLOBAL_2D_ARRAY,2)
      real uwnd(GLOBAL_2D_ARRAY) ! Stores time interpolated u-wind
      ! grid type are needed for ncdf_read in read_write ncdf module
      integer :: uwnd_grd_type = 0 ! should be = r_var but need to include ncvars.h
      character(len=4) :: uwnd_name = 'uwnd'
      ! timestep variables for wind bulk force
      ! Notice time here is called something different wnd not uwnd.
      ! get_smth.F has a method to deal with different time names, might
      ! be worth using. Hard coded now for ease.
      character(len=8) :: uwnd_time_name = 'wnd_time'
      integer :: uwnd_file_indx = 0 ! Set file index in force file list as zero to begin with
      integer :: uwnd_irec      = 0 ! Set netcdf record number to zero to start
      ! These were only one variable (itwnd) in old code, but maybe it is clearer
      ! with 2. Although if it's one it's forcing you to realise that
      ! it jumps between values.
      integer :: uwnd_it1 = 1, uwnd_it2 = 2  ! Placeholder for read in time
      real, dimension(2) :: uwnd_times = [-99,-99] ! Store input file times

       ! v-wind (input data in m/s)
       ! ------
      real    :: vwnd_data(GLOBAL_2D_ARRAY,2)
      real    :: vwnd(GLOBAL_2D_ARRAY)
      integer :: vwnd_grd_type  = 0
      integer :: vwnd_file_indx = 0
      integer :: vwnd_irec      = 0
      integer :: vwnd_it1 = 1, vwnd_it2 = 2
      real, dimension(2) :: vwnd_times = [-99,-99]
      character(len=4)   :: vwnd_name = 'vwnd'
      character(len=8)   :: vwnd_time_name = 'wnd_time'

      ! Atmospheric tracers
      ! -------------------

       ! tair (air temp)
       ! ---------------
       ! Input data units: Degrees C
      real    :: tair_data(GLOBAL_2D_ARRAY,2)
      real    :: tair(GLOBAL_2D_ARRAY)
      integer :: tair_grd_type  = 0
      integer :: tair_file_indx = 0
      integer :: tair_irec      = 0
      integer :: tair_it1 = 1, tair_it2 = 2
      real, dimension(2) :: tair_times = [-99,-99]
      character(len=4)   :: tair_name = 'tair'
      character(len=8)   :: tair_time_name = 'tra_time'

       ! rhum (currently specific humidity! Not relative hum.) qair in old code
       ! ----
       ! Search for label(HUMIDITY) to convert to relative humidity instead
       ! if input data is in relative humidity. Currently set for input
       ! data in specific humidity
       ! Input data units: fraction
      real    :: rhum_data(GLOBAL_2D_ARRAY,2)
      real    :: rhum(GLOBAL_2D_ARRAY)
      integer :: rhum_grd_type  = 0
      integer :: rhum_file_indx = 0
      integer :: rhum_irec      = 0
      integer :: rhum_it1 = 1, rhum_it2 = 2
      real, dimension(2) :: rhum_times = [-99,-99]
      character(len=4)   :: rhum_name = 'qair'
      character(len=8)   :: rhum_time_name = 'tra_time'

       ! Precifitation rate (a.k.a. rain fall)
       ! ------------------
       ! Input data units: [kg/m^2/s^1], which is converted to m/s later.
      real    :: prate_data(GLOBAL_2D_ARRAY,2)
      real    :: prate(GLOBAL_2D_ARRAY)
      integer :: prate_grd_type  = 0
      integer :: prate_file_indx = 0
      integer :: prate_irec      = 0
      integer :: prate_it1 = 1, prate_it2 = 2
      real, dimension(2) :: prate_times = [-99,-99]
      character(len=4)   :: prate_name = 'rain'
      character(len=9)   :: prate_time_name = 'prec_time'

       ! net terrestrial longwave radiation [Watts/m^2]
       ! ----------------------------------
      real    :: lwrad_data(GLOBAL_2D_ARRAY,2)
      real    :: lwrad(GLOBAL_2D_ARRAY)
      integer :: lwrad_grd_type  = 0
      integer :: lwrad_file_indx = 0
      integer :: lwrad_irec      = 0
      integer :: lwrad_it1 = 1, lwrad_it2 = 2
      real, dimension(2) :: lwrad_times = [-99,-99]
      character(len=5)   :: lwrad_name = 'radlw'
      character(len=8)   :: lwrad_time_name = 'rad_time'

       ! net solar shortwave radiation [Watts meter-2]
       ! -----------------------------
      ! swrad_data (as swradg) and srflx needed in init_arrays.F to comply with old code
      real, public    :: swrad_data(GLOBAL_2D_ARRAY,2)
      real, public    :: srflx(GLOBAL_2D_ARRAY) ! Kept same name convention as old code
      integer :: swrad_grd_type  = 0
      integer :: swrad_file_indx = 0
      integer :: swrad_irec      = 0
      integer :: swrad_it1 = 1, swrad_it2 = 2
      real, dimension(2) :: swrad_times = [-99,-99]
      character(len=5)   :: swrad_name = 'radsw'
      character(len=8)   :: swrad_time_name = 'rad_time'


      ! OTHER VARIABLES
      ! ===============

      ! /* Pierre sustr-OUTPUT */
      real sustr_blk(GLOBAL_2D_ARRAY)
      real svstr_blk(GLOBAL_2D_ARRAY)
      real shflx_rsw(GLOBAL_2D_ARRAY)
      real shflx_lat(GLOBAL_2D_ARRAY)
      real shflx_sen(GLOBAL_2D_ARRAY)
      real shflx_rlw(GLOBAL_2D_ARRAY)
      real shflx_tot(GLOBAL_2D_ARRAY)

      public set_bulk_frc_tile

      contains


! ----------------------------------------------------------------------

      subroutine set_bulk_frc_tile(istr,iend,jstr,jend)

      ! ===================
      ! SET ALL BULK FORCES
      ! ===================

      ! Landing subroutine from set_forces
      ! where all bulk forces are set.

      ! All bulk variables need time interpolation only
      ! here so can use same generic routine

      ! ALL FORCING INPUT DATA NEEDS TO BE IN DAYS!

      implicit none

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend

      ! 1) Find & save force variable file_indx & correct record (irec)
      ! 2) Read in forcing data
      ! 3) Interpolate input forcing to correct model time

      ! Short/Long radiation
      ! --------------------

        ! Short wave radiation -
        ! this accounts for get_srflx.F & time interp of set_srflx.F in newer code
        ! and get_srflux.F & set_srflux() in old code
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          swrad_name,     swrad_time_name,! Text names
     &          swrad_data,     srflx,          ! Variable arrays
     &          swrad_times,    swrad_grd_type, ! Input times & Grid type
     &          swrad_file_indx,swrad_irec,     ! File indx & Current input record
     &          swrad_it1,      swrad_it2 )     ! Time index placeholders
!      mpi_master_only print *, 'srflx(20,20)',srflx(20,20)

        ! Long wave radiation - get_bulk_rad.F & set_bulk_rad() in old code
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          lwrad_name,     lwrad_time_name,! Text names
     &          lwrad_data,     lwrad,          ! Variable arrays
     &          lwrad_times,    lwrad_grd_type, ! Input times & Grid type
     &          lwrad_file_indx,lwrad_irec,     ! File indx & Current input record
     &          lwrad_it1,      lwrad_it2 )     ! Time index placeholders
!      mpi_master_only print *, 'lwrad(20,20)',lwrad(20,20)

      ! Wind (get_bulk_wnd.F & set_bulk_wnd() in old code)
      ! ----

        ! u-wind -
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          uwnd_name,      uwnd_time_name, ! Text names
     &          uwnd_data,      uwnd,           ! Variable arrays
     &          uwnd_times,     uwnd_grd_type,  ! Input times & Grid type
     &          uwnd_file_indx, uwnd_irec,      ! File indx & Current input record
     &          uwnd_it1,       uwnd_it2 )      ! Time index placeholders
!      mpi_master_only print *, 'uwnd(20,20)',uwnd(20,20)

        ! v-wind
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          vwnd_name,      vwnd_time_name, ! Text names
     &          vwnd_data,      vwnd,           ! Variable arrays
     &          vwnd_times,     vwnd_grd_type,  ! Input times & Grid type
     &          vwnd_file_indx, vwnd_irec,      ! File indx & Current input record
     &          vwnd_it1,       vwnd_it2 )      ! Time index placeholders
!      mpi_master_only print *, 'vwnd(20,20)',vwnd(20,20)

      ! Atmospheric tracers (get_bulk_tra & set_bulk_tra() in old code)
      ! -------------------

        ! tair (air temp)
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          tair_name,      tair_time_name, ! Text names
     &          tair_data,      tair,           ! Variable arrays
     &          tair_times,     tair_grd_type,  ! Input times & Grid type
     &          tair_file_indx, tair_irec,      ! File indx & Current input record
     &          tair_it1,       tair_it2 )      ! Time index placeholders
!      mpi_master_only print *, 'tair(20,20)',tair(20,20)

        ! rhum (currently specific (not relative) humidity) qair in old code
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          rhum_name,      rhum_time_name, ! Text names
     &          rhum_data,      rhum,           ! Variable arrays
     &          rhum_times,     rhum_grd_type,  ! Input times & Grid type
     &          rhum_file_indx, rhum_irec,      ! File indx & Current input record
     &          rhum_it1,       rhum_it2 )      ! Time index placeholders
!      mpi_master_only print *, 'rhum(20,20)',rhum(20,20)

        ! prate (precipitation rate) - 'rain' variable in old code
        ! get_bulk_prec.F & set_bulk_prec() in old code
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          prate_name,     prate_time_name,! Text names
     &          prate_data,     prate,          ! Variable arrays
     &          prate_times,    prate_grd_type, ! Input times & Grid type
     &          prate_file_indx,prate_irec,     ! File indx & Current input record
     &          prate_it1,      prate_it2 )     ! Time index placeholders
!      mpi_master_only print *, 'prate(20,20)',prate(20,20)


      ! NEED TO ADD QCORRECTOR CALLS FOR SEA SURFACE TEMP AND/OR SALINITY
      ! =================================================================

      ! CALCULATE FLUX INTERACTIONS
      ! ===========================
      call calc_all_bulk_forces(istr,iend,jstr,jend)


      end subroutine set_bulk_frc_tile


!-----------------------------------------------------------------------

      subroutine calc_all_bulk_forces(istr,iend,jstr,jend)

      ! =========================================
      ! CALCULATE INTERACTIVE BULK FORCE TRANSFER
      ! =========================================

      ! This started as a copy of bulk_flux_COARE.F
      ! from the new code


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
      use flux_all, only: sustr, svstr, stflx

      implicit none

! THIS WAS UNCOMMENTED!!!! DevinD - JM says wind must be at rho-points
!# ifndef WIND_AT_RHO_POINTS
!???? ERROR: bulk_flux requires WIND_AT_RHO_POINTS
!# endif

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend

!
!# include "param.h" ! included at the top of the module
# include "scalars.h"
# include "grid.h"
# include "ocean3d.h"

# ifdef SEA_ICE_NOFLUX
#  ifdef LMD_KPP
#   include "mixing.h"
#  endif
# endif

#define BULK_LW
#undef SST_SKIN
#define BULK_EP

      integer i,j, IterMax,iter
      real a
      real rho0i
      real TseaC,TseaK,Qsea
      real TairC,TairK,rhoAir,Qair
      real Q,RH,VisAir,Hlv
      real delW,delT,delQ
      real u10,Zo10,Cd10,Ch10,Ct10,Cd
      real Ct,CC,Ri,Ribcu,Zetu,L10
      real Wstar,Tstar,Qstar
      real ZoW,ZoT,ZoT10,ZoQ,ZoL,L,Rr,Bf
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
      parameter (blk_ZW=10.0)
      parameter (blk_ZT=10.0)
      parameter (blk_ZQ=10.0)
      parameter (blk_Zabl=600.0)
      parameter (blk_beta=1.2)
      parameter (blk_Cpa=1004.67)
# ifdef BULK_LW
      real emiss_lw,SigmaSB
      parameter (emiss_lw=0.985)        ! longwave emissivity of ocean,
                                        ! ref. Dickey et al. (1994)
      parameter (SigmaSB=5.6697E-8)     ! Stefan-Boltzmann constant [m^2/K^4]
# endif
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

      real wspd ! uwnd, vwnd, - DevinD
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
!      - radlw         #ifdef BULK_LW: downward longwave radiation)
!                      #else: net outgoing longwave radiation)
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
                         ! [kg*degC/Joule], Cp is defined in "scalars.h"

!!!     &      cpaircpi=cpair*cpi ! ratio of specific heat for air
!!!                               ! and seawater
!!!# ifdef BULK_LW
!!!      real, parameter ::
!!!     &      CtoK=273.16,       ! conversion from [C] to [K]
!!!# endif

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

! DevinD commented forces.h as everything should be module contained
!# include "forces.h"



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
!# define SWRAD_DATA
!# include "param.h"
!# include "scalars.h"
!# include "grid.h"
!# include "forces.h"
!
!# include "compute_extended_bounds.h"
!
!
# ifdef ROBUST_DIURNAL_SRFLUX
      dtau=2.D0*pi * dt*sec2day
      cff=4.D0*pi*(tdays-dble(int(tdays))-0.5D0)
      tau1=max(-pi, min(pi, cff))
      tau2=max(-pi, min(pi, cff+2.D0*dtau))
      cff=tau2-tau1
      Ampl=( cff + 2.D0*cos(0.5D0*(tau1+tau2))*sin(0.5D0*cff) )/dtau
# endif
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

! inverse seawater density, use mean value for seawater density.

      rho0i=1.0/rho0

!      tmid=time+0.5*dt ! DevinD commented
!> wind (m/s) - DevinD purely time interp, done already
!      it2=itwnd ; cff1=wnd_time(it2)-tmid
!      it1=3-it2 ; cff2=tmid-wnd_time(it1)
!      if (cff1 < dt) synchro_flag=.true.
!      cff=1./(cff1+cff2)
!      cff1=cff1*cff ; cff2=cff2*cff

!> tair (degC) - DevinD purely time interp, done already
!      it4=ittair ; cff3=tair_time(it4)-tmid
!      it3=3-it4  ; cff4=tmid-tair_time(it3)
!      if (cff3 < dt) synchro_flag=.true.
!      cff=1./(cff3+cff4)
!      cff3=cff3*cff ; cff4=cff4*cff

!> rhum (frac RELATIVE HUMIDITY)
! convert to fraction from percent
!      it6=itrhum ; cff5=rhum_time(it6)-tmid
!      it5=3-it6  ; cff6=tmid-rhum_time(it5)
!      if (cff5 < dt) synchro_flag=.true.
!      cff=1./(100.*(cff5+cff6))
!      cff5=cff5*cff ; cff6=cff6*cff

!> rhum (frac of SPECIFIC HUMIDITY)  - DevinD purely time interp, done already
!      it6=itrhum ; cff5=rhum_time(it6)-tmid
!      it5=3-it6  ; cff6=tmid-rhum_time(it5)
!      if (cff5 < dt) synchro_flag=.true.
!      cff=1./(cff5+cff6)
!      cff5=cff5*cff ; cff6=cff6*cff

!> prate (m/s)
!      it8=itprate ; cff7=prate_time(it8)-tmid
!      it7=3-it8   ; cff8=tmid-prate_time(it7)
!      if (cff7 < dt) synchro_flag=.true.

! Drakkar rain is in units of [kg/m^2/s^2]. convert it to [m/s]
!! UNITS?? isn't it [kg/m^2/s]?
!cff7 and cff8 units of 1/precip?
! using fresh water density rho=1000.

! SY New units of cm/day--convert to m/s.
!      cff=1./(100.*3600.*24.*(cff7+cff8))
!      cff7=cff7*cff ; cff8=cff8*cff
! units of [kg/m^2/s^1] --convert to m/s.
!      cff=1./(1000.*(cff7+cff8))
!      cff7=cff7*cff ; cff8=cff8*cff


!> lwrad (degC.m/s)
!      it10=itlrf ; cff9=lrf_time(it10)-tmid
!      it9=3-it10 ; cff10=tmid-lrf_time(it9)
!      if (cff9 < dt) synchro_flag=.true.
!      cff=1./(rho0*Cp*(cff9+cff10))     ! <-- convert [Watts/m^2] into
!      cff9=cff9*cff ; cff10=cff10*cff   !     kinematic flux [degC*m/s]

# ifdef QCORRECTION
!      it12=itdqdt ; cff11=dqdt_time(it12)-tmid  ! Flux correction for
!      it11=3-it12 ; cff12=tmid-dqdt_time(it11)  ! net surface heat flux
!      if (cff11 < dt) synchro_flag=.true.       ! (relaxation toward
!      cff=1./(rho0*Cp*(cff11+cff12))            ! observed SST data)
!      cff11=cff11*cff ; cff12=cff12*cff         ! and, optionally, ...
!
!      it14=itsst  ; cff13=sst_time(it14)-tmid
!      it13=3-it14 ; cff14=tmid-sst_time(it13)
!      if (cff13 < dt) synchro_flag=.true.
!      cff=1/(cff13+cff14)
!      cff13=cff13*cff ; cff14=cff14*cff
# endif
# ifdef SFLX_CORR
!      it16=itsss  ; cff15=sss_time(it16)-tmid   ! the same for surface
!      it15=3-it16 ; cff16=tmid-sss_time(it15)   ! freshwater flux
!      if (cff15 < dt) synchro_flag=.true.
!      cff=1./(cff15+cff16)
!      cff15=cff15*cff ; cff16=cff16*cff
# endif


      do j=jstrR,jendR
        do i=istrR,iendR

          ! DevinD uwnd & vwnd already interpolated???
!          uwnd=cff1*uwind(i,j,it1)+cff2*uwind(i,j,it2)
!          vwnd=cff1*vwind(i,j,it1)+cff2*vwind(i,j,it2)
          ! DevinD is wspd for rho-point or u/v-point calc?
          ! Currently setup for rho-point
          wspd=sqrt(uwnd(i,j)*uwnd(i,j)+vwnd(i,j)*vwnd(i,j))

!          tair=cff3*tairg(i,j,it3)  +cff4*tairg(i,j,it4) ! DevinD time only
!          rhum=cff5*rhumg(i,j,it5)  +cff6*rhumg(i,j,it6) ! DevinD time only
!          prate=cff7*prateg(i,j,it7)+cff8*prateg(i,j,it8) ! DevinD divide 1000 below

          ! convert [Watts/m^2] into kinematic flux [degC*m/s]
          radlw=lwrad(i,j) / (rho0*Cp)

          wspd0=wspd
          TairC=tair(i,j) ! DevinD
          TairK=TairC+273.16
# ifdef SST_SKIN
          TseaC=sst_skin(i,j)   ! sea surface temperature and salinity
# else
          TseaC=t(i,j,N,nrhs,itemp)
# endif
          TseaK=TseaC+273.16

          ! label(HUMIDITY) - chose between specific or relative humidity
          RH=rhum(i,j)        ! Specific humidity (if input data is in specific hum.)
!          RH=rhum(i,j)/100.0 ! Specific humidity (if input data is in relative hum.)

!!!          Ssea=t(i,j,N,nrhs,isalt)

# ifndef BULK_LW
          hflw=-radlw ! positive downward, i.e., warming the ocean
# else
!
!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
          hflw=radlw  ! positive downward
     &        - emiss_lw*rho0i*cpi*SigmaSB*TseaK*TseaK*TseaK*TseaK

# endif
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

          cff=(1.0007+3.46e-6*patm)*6.1121
     &        *exp(17.502*TairC/(240.97+TairC))

! specific humidity at Saturation, Qair (kg/kg).

          Qair=0.62197*(cff/(patm-0.378*cff))

! specific humidity, Q (kg/kg).
!
!          if (RH < 2.) then                          !RH fraction
!            cff=cff*RH                                 !Vapor pres (mb)
!            Q=0.62197*(cff/(patm-0.378*cff))           !Spec hum (kg/kg)
!          else          !RH input was actually specific humidity in g/kg
!            Q=RH/1000.0                                !Spec Hum (kg/kg)
!          endif

!   RH is already fraction of SPECIFIC humidity in kg/kg
             Q=RH

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
     &                              (1.0+0.61*Q))
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
          delQ=Qsea-Q
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
     &             (Tstar*(1.0+0.61*Q)+0.61*TairK*Qstar)/
     &             (TairK*Wstar*Wstar*(1.0+0.61*Q)+eps)
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
          Bf=-g/TairK*Wstar*(Tstar+0.61*TairK*Qstar)
          if (Bf.gt.0.0) then
            Wgus=blk_beta*(Bf*blk_Zabl)**r3
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
          hfsen=-blk_Cpa*rhoAir*Wstar*Tstar
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
          hflat=-Hlv*rhoAir*Wstar*Qstar
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
          upvel=-1.61*Wstar*Qstar-(1.0+1.61*Q)*Wstar*Tstar/TairK
          hflat=hflat+rhoAir*Hlv*upvel*Q
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
# ifdef BULK_EP
!
!  Total surface salt flux
!
          evap=-cp*hflat/Hlv
          ! DevinD divide prate by 1000: units of [kg/m^2/s^1] convert to m/s.
          stflx(i,j,isalt)=(evap- (prate(i,j)/1000) )*t(i,j,N,nrhs,isalt)
# endif

# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
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
          sustr(i,j)=aer(i,j)*cer(i,j)*uwnd(i,j)
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif
          svstr(i,j)=aer(i,j)*cer(i,j)*vwnd(i,j)
#  ifdef MASKING
     &                  *rmask(i,j)
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
! Pierre : This is shit, yu want to change this !!
# ifdef SFLX_CORR
          cff=cff11*dqdtg(i,j,it11)+cff12*dqdtg(i,j,it12)
          stflx(i,j,isalt)=stflx(i,j,isalt) +cff*(
     &           t(i,j,N,nrhs,isalt) -cff15*sssg(i,j,it15)
     &                               -cff16*sssg(i,j,it16)
     &                                                  )
#  ifdef MASKING
     &                                           *rmask(i,j)
#  endif
# endif

!!# ifdef SFLX_CORR
!!          stflx(i,j,isalt)=stflx(i,j,isalt)
!!# ifdef SFLX_CORR_COEF
!!     &  -Hz(i,j,N)/(dSdt*day2sec)*(t(i,j,N,nrhs,isalt)-sss(i,j))
!!# else
!!     &                 +dqdt(i,j)*(t(i,j,N,nrhs,isalt)-sss(i,j))
!!# endif
!!#  ifdef WET_DRY
!!     &                                           *rmask_wet(i,j)
!!#  endif
!!# endif /* SFLX_CORR */

!---------------------------------------------------------------
!
!---------------------------------------------------------------
!

! Restrict stflx to prevent surface temperature to go below -2
! degrees C.
!
!!# ifdef SEA_ICE_NOFLUX
!!              stflx(i,j,itemp)=max( stflx(i,j,itemp),
!!#  ifdef LMD_KPP
!!     &                          max(hbls(i,j,nstp), Hz(i,j,N))
!!#  else
!!     &                                               Hz(i,j,N)
!!#  endif
!!     &                         *(-2.D0-t(i,j,N,nrhs,itemp))/dt
!!     &                                                       )
!!# endif

#ifdef SEA_ICE_NOFLUX
              if( t(i,j,N,nrhs,itemp) .le. -1.8 ) then
                   stflx(i,j,itemp)=0.
#   if defined LMD_KPP
                   srflx(i,j)=0.
#    endif
              endif
#endif
!-------------------------------------------------
!
!  Store fluxes for diagnostics
!  When storing in netCDF file, multiplied
!  by rho0.Cp to get W/m2
!
          shflx_rsw(i,j)=srflx(i,j)*rho0*Cp
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif
          shflx_lat(i,j)=hflat*rho0*Cp
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif
          shflx_sen(i,j)=hfsen*rho0*Cp
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif
          shflx_rlw(i,j)=hflw*rho0*Cp
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif
          shflx_tot(i,j)=stflx(i,j,itemp)*rho0*Cp
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif


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

      do j=jstrR,jendR
        do i=istrR,iend
!          uwnd=cff1*uwind(i,j,it1)+cff2*uwind(i,j,it2)
!          vwnd=cff1*vwind(i,j,it1)+cff2*vwind(i,j,it2)
!          wspd=sqrt( uwnd(i,j)*uwnd(i,j) + vwnd(i,j)*vwnd(i,j) )
          ! Rho-2-u point averaging
          wspd=0.5*sqrt(  (uwnd(i,j)+uwnd(i-1,j))**2
     &                  + (vwnd(i,j)+vwnd(i-1,j))**2 )
          if (wspd .gt. Wspd_min) then
            cff = cfb_slope*wspd + cfb_offset
          else
            cff = stau_ref
          endif
!          sustr(i,j)=sustr(i,j) + cff*0.5*(u(i,j,N,nrhs)+u(i+1,j,N,nrhs))/rho0
          !! This rho0 should be worked into the coefficients
          ! NOTE: sustr is now at u-point not rho-point
      ! DevinD debug
!      if (mynode==0 .and. i==20 .and. j==20) then
!        print *, 'u(20,20,N,nrhs)',u(20,20,N,nrhs)
!        print *, 'wspd',wspd,'cff',cff
!      endif
          sustr(i,j)=sustr(i,j) + cff*u(i,j,N,nrhs)/rho0
          !! Do a mask multiply at reading time for input vars ?
# ifdef MASKING
          !! I don't think that this is actually doing anything usefull
          ! This would be improved with & *umask instead of it's own line.
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
        enddo
      enddo

      do j=jstrR,jend
        do i=istrR,iendR
!          uwnd=cff1*uwind(i,j,it1)+cff2*uwind(i,j,it2)
!          vwnd=cff1*vwind(i,j,it1)+cff2*vwind(i,j,it2)
!          wspd=sqrt(uwnd(i,j)*uwnd(i,j) + vwnd(i,j)*vwnd(i,j) )
          ! Rho-2-v point averaging
          wspd=0.5*sqrt(  (uwnd(i,j)+uwnd(i,j-1))**2
     &                  + (vwnd(i,j)+vwnd(i,j-1))**2 )
          if (wspd .gt. Wspd_min) then
            cff = cfb_slope*wspd + cfb_offset
          else
            cff = stau_ref
          endif
!          svstr(i,j)=svstr(i,j) + cff*0.5*(v(i,j,N,nrhs)+v(i,j+1,N,nrhs))/rho0
          ! NOTE: sustr is now at v-point not rho-point
          svstr(i,j)=svstr(i,j) + cff*v(i,j,N,nrhs)/rho0
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
        enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! DIAGNOSTICS !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j=jstrR,jendR
        do i=istrR,iendR
          sustr_blk(i,j)= sustr(i,j)*rho0  !<== diagnostics
          svstr_blk(i,j)= svstr(i,j)*rho0  !<== diagnostics
        enddo
      enddo


      ! DevinD debug
!      mpi_master_only print *, 'sustr(20,20)',sustr(20,20)
!      mpi_master_only print *, 'svstr(20,20)',svstr(20,20)
!      mpi_master_only print *, 'stflx(20,20,itemp)',stflx(20,20,itemp)
!      mpi_master_only print *, 'stflx(20,20,isalt)',stflx(20,20,isalt)


      end subroutine calc_all_bulk_forces


!-----------------------------------------------------------------------

      function bulk_psiu (ZoL, pi)
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
      end function bulk_psiu


!-----------------------------------------------------------------------

      function bulk_psit (ZoL, pi)

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
      end function bulk_psit


!-----------------------------------------------------------------------

#endif /* BULK_FLUX for whole module */

      end module bulk_frc
