#include "cppdefs.h"
#ifdef BULK_FLUX
# ifndef WIND_AT_RHO_POINTS
???? ERROR: bulk_flux requires WIND_AT_RHO_POINTS
# endif


      subroutine bulk_flux_tile(istr,iend,jstr,jend)

! Computes the bulk parameterization of surface wind stress and
! surface net heat fluxes using COAMPS formulation.

      implicit none
      integer istr,iend,jstr,jend, i,j

# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "ocean3d.h"
# ifdef SEA_ICE_NOFLUX
#  ifdef LMD_KPP
#   include "mixing.h"
#  endif
# endif



      real uwnd,vwnd, wspd,wspd0, szero, ss
      real Tsea,Ssea,evap, radlw
      real rho0i, prate
      real ew,rw,qw,es, rhum, rs,ra,qa
      real tair, vtair,rhoair,a,ce,cs,htevap
      real hflat,hfsen


#define BULK_LW
#undef SST_SKIN
#define BULK_EP



      real, parameter ::
     &      patm=1013.0,     ! use mean value for atmospheric pressure
                             !                          at sea surface

     &      cpair=0.240*4186.8, ! specific heat for air in
                                ! Joules/(kg*degC), where 4186.8 is
                                ! to convert from cal/g to Joules/kg
     &      cpi=1./Cp,   ! inverse of specific heat for seawater
                         ! [kg*degC/Joule], Cp is defined in "scalars.h"

     &      cpaircpi=cpair*cpi ! ratio of specific heat for air
                               ! and seawater
# ifdef BULK_LW
      real, parameter ::
     &      CtoK=273.16,       ! conversion from [C] to [K]
     &      emiss_lw=0.985,    ! longwave emissivity of ocean,
                               ! ref. Dickey et al. (1994)
     &      SigmaSB=5.6697e-8  ! Stefan-Boltzmann constant [m^2/K^4]
# endif

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
# include "forces.h"



# include "compute_extended_bounds.h"

! inverse seawater density, use mean value for seawater density.

      rho0i=1.0/rho0

      tmid=time+0.5*dt
!> wind
      it2=itwnd ; cff1=wnd_time(it2)-tmid 
      it1=3-it2 ; cff2=tmid-wnd_time(it1)
      if (cff1 < dt) synchro_flag=.true.
      cff=1./(cff1+cff2)
      cff1=cff1*cff ; cff2=cff2*cff

!> tair
      it4=ittair ; cff3=tair_time(it4)-tmid
      it3=3-it4  ; cff4=tmid-tair_time(it3)
      if (cff3 < dt) synchro_flag=.true.
      cff=1./(cff3+cff4)
      cff3=cff3*cff ; cff4=cff4*cff

!> rhum
      it6=itrhum ; cff5=rhum_time(it6)-tmid
      it5=3-it6  ; cff6=tmid-rhum_time(it5) 
      if (cff5 < dt) synchro_flag=.true.
      cff=1./(cff5+cff6)
      cff5=cff5*cff ; cff6=cff6*cff

!> prate
      it8=itprate ; cff7=prate_time(it8)-tmid
      it7=3-it8   ; cff8=tmid-prate_time(it7) 
      if (cff7 < dt) synchro_flag=.true.

! Drakkar rain is in units of [kg/m^2/s^2]. convert it to [m/s]
! using fresh water density rho=1000.

      cff=1./(1000.*(cff7+cff8))
      cff7=cff7*cff ; cff8=cff8*cff

!> lwrad
      it10=itlrf ; cff9=lrf_time(it10)-tmid 
      it9=3-it10 ; cff10=tmid-lrf_time(it9)
      if (cff9 < dt) synchro_flag=.true.
      cff=1./(rho0*Cp*(cff9+cff10))     ! <-- convert [Watts/m^2] into
      cff9=cff9*cff ; cff10=cff10*cff   !     kinematic flux [degC*m/s]

# ifdef QCORRECTION
      it12=itdqdt ; cff11=dqdt_time(it12)-tmid  ! Flux correction for
      it11=3-it12 ; cff12=tmid-dqdt_time(it11)  ! net surface heat flux
      if (cff11 < dt) synchro_flag=.true.       ! (relaxation toward
      cff=1./(rho0*Cp*(cff11+cff12))            ! observed SST data)
      cff11=cff11*cff ; cff12=cff12*cff         ! and, optionally, ...

      it14=itsst  ; cff13=sst_time(it14)-tmid
      it13=3-it14 ; cff14=tmid-sst_time(it13)
      if (cff13 < dt) synchro_flag=.true.
      cff=1/(cff13+cff14)
      cff13=cff13*cff ; cff14=cff14*cff
# endif
# ifdef SFLX_CORR
      it16=itsss  ; cff15=sss_time(it16)-tmid   ! the same for surface
      it15=3-it16 ; cff16=tmid-sss_time(it15)   ! freshwater flux
      if (cff15 < dt) synchro_flag=.true.
      cff=1./(cff15+cff16)
      cff15=cff15*cff ; cff16=cff16*cff
# endif


      do j=jstrR,jendR
        do i=istrR,iendR
          uwnd=cff1*uwind(i,j,it1)+cff2*uwind(i,j,it2)
          vwnd=cff1*vwind(i,j,it1)+cff2*vwind(i,j,it2)
          wspd=sqrt(uwnd*uwnd+vwnd*vwnd)

          tair=cff3*tairg(i,j,it3)  +cff4*tairg(i,j,it4)
          rhum=cff5*rhumg(i,j,it5)  +cff6*rhumg(i,j,it6)
          prate=cff7*prateg(i,j,it7)+cff8*prateg(i,j,it8)

          radlw=cff9*lwradg(i,j,it9)+cff10*lwradg(i,j,it10)


# ifdef SST_SKIN
          Tsea=sst_skin(i,j)   ! sea surface temperature and salinity
# else
          Tsea=t(i,j,N,nrhs,itemp)
# endif
          Ssea=t(i,j,N,nrhs,isalt)

! saturation vapor pressure (ew) at surface of water.
! use formula of Buck.
! note:  change due to salinity correction is about 1.7%.
! simplify for constant atm pressure (1013 mb).

          ew=6.1121*exp(17.502*Tsea/(240.97 + Tsea))
     &                  *1.0042*(1.0 - 5.37E-4*Ssea)

! saturation mixing ratio at surface of water (rw).

          rw=0.62197*ew/(patm-ew)

! saturation specific humidity at surface of water (qw).

          qw=rw/(1.0+rw)

! saturation vapor pressure (es) of air at temp tair.
! simplify for constant atm pressure (1013 mb) and zero salinity.

          es=6.1121*exp(17.502*tair/(240.97+tair))*1.0042

! saturation mixing ratio for air (rs).

          rs=0.62197*es/(patm-es)

! mixing ratio for air "ra" from relative humidity.

          ra=rhum*rs

! specific humidity for air.

          qa=ra/(1.0+ra)

! air density in kg/m3 (vtair = virtual temperature of air).
! formula is from Smithsonian Meteor Tables (List, 1951, p 290,295).
! note:  the correction for water vapor is ~ 1%, so neglect.

          vtair=(tair+273.16)*(1.0+1.60779*ra)/(1.0+ra)
          rhoair=0.34838*patm/vtair

! drag coefficients for latent and sensible heat exchange.

! drag coefficient calculation here is from Kondo, 1975.
! use min windspeed of 0.1 m/s.
! note: using neutral values of 0.0014 for "ce" and 0.0011 for "cs"
! results in mean latent and sensible heat flux magnitudes that
! match the values computed by COAMS fairly well.

          wspd0=max(0.1, wspd)
          szero=(Tsea-tair)/(wspd0*wspd0)
          ss=szero*abs(szero)/(abs(szero)+0.01)
          ss=max(ss, -3.3)
          if (ss <= 0.) then
            a=0.1+0.03*ss+0.9*exp(4.8*ss)
          else
            a=1.0+0.63*sqrt(ss)
          endif
          ce=0.0014*a
          cs=0.0011*a    !--> discard a,ss

! latent heat of evaporation in Joules/kg; the formula is
!          htevap = 597.8-0.573*T(degC) cal/g.
! to convert from cal/g to Joules/kg, muliply by 4186.8.

          htevap=2.5029e6-2.40e3*Tsea

! evaporation rate [m/s], latent and sensible heat fluxes [deg_C-m/s]
! sign convention: positive latent and sensible heat fluxes downward;

          a=rhoair*wspd0*rho0i
          evap=a*ce*(qw-qa)
          hflat=-evap*htevap*cpi
          hfsen=-a*cs*cpaircpi*(Tsea-tair)

! Net surface heat flux.
!---- ------- ---- -----

          stflx(i,j,itemp)=srflx(i,j)+hflat+hfsen
# ifdef BULK_LW
     &      +emiss_lw*( radlw-rho0i*cpi*sigmaSB*(Tsea+CtoK)**4 )
# else
     &                                -radlw
# endif
# ifdef BULK_EP
          stflx(i,j,isalt)=(evap-prate)*Ssea
# endif
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
# endif

! kinematic surface wind stress [m^2/s^2]
!---------- ------- ---- ------ ---------
! sustr = rhoA/rhoW * Cd * wspd*uwnd
! svstr = rhoA/rhoW * Cd * wspd*vwnd

# ifdef BULK_SM_UPDATE
          sustr(i,j)=a*ce*uwnd
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif
          svstr(i,j)=a*ce*vwnd
#  ifdef MASKING
     &                  *rmask(i,j)
#  endif
# endif

! NOTE THAT "dqdt" NEGATIVE.

# ifdef QCORRECTION
          cff=cff11*dqdtg(i,j,it11)+cff12*dqdtg(i,j,it12)

          stflx(i,j,itemp)=stflx(i,j,itemp) +cff*(
     &           t(i,j,N,nrhs,itemp) -cff13*sstg(i,j,it13)
     &                               -cff14*sstg(i,j,it14)
     &                                                  )
# endif
# ifdef SFLX_CORR
          stflx(i,j,isalt)=stflx(i,j,isalt) -dSSSdt*(
     &           t(i,j,N,nrhs,isalt) -cff15*sssg(i,j,it15)
     &                               -cff16*sssg(i,j,it16)
     &                                                  )
# endif

! Restrict stflx to prevent surface temperature to go below -2
! degrees C.

# ifdef SEA_ICE_NOFLUX
              stflx(i,j,itemp)=max( stflx(i,j,itemp),
#  ifdef LMD_KPP
     &                          max(hbls(i,j,nstp), Hz(i,j,N))
#  else
     &                                               Hz(i,j,N)
#  endif
     &                         *(-2.D0-t(i,j,N,nrhs,itemp))/dt
     &                                                       )
# endif
        enddo
      enddo
      end

#else
      subroutine bulk_flux_empty
      end
#endif
