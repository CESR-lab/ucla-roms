#include "cppdefs.h"
#ifdef BULK_FLUX
# ifndef WIND_AT_RHO_POINTS
???? ERROR: bulk_flux requires WIND_AT_RHO_POINTS
# endif

      subroutine bulk_flux_tile(istr,iend,jstr,jend, aer,cer)

! Compute bulk parameterization for surface wind stress and net heat
! fluxes. Adapted from COARE code originally written by David Rutgers
! and Frank Bradley. EMINUSP option for equivalent salt fluxes added
!                                          by Paul Goodman (10/2004).
! Modified by Kate Hedstrom for COARE version 3.0 (03/2005).
! Modified by Jim Edson to correct specific hunidities.
! Modified for ROMS_AGRIF by Fred Marin (2008)
! Modified for ROMS_AGRIF by P. Marchesiello (2009)

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



      implicit none
      integer istr,iend,jstr,jend
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: aer, cer
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
      integer i,j, IterMax,iter
      real a,cff
      real rho0i,cpi
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

      real bulk_psiu,bulk_psit

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
      parameter (emiss_lw=0.985)
      parameter (SigmaSB=5.6697E-8)
# endif
      real rhow,patm
      parameter (rhow=1000.0)
      parameter (patm=1010.0)
      real eps,r3
      parameter (eps=1.d-20)
      parameter (r3=1.0/3.0)

# ifdef QCORRECTION
#  define SST_DATA
      integer it7,it8, it9,it10
      real tmid, cff7,cff8, cff9,cff10, cff
#  ifdef SFLX_CORR
      integer it11,it12
      real cff11,cff12
#  endif
# endif
# include "forces.h"

#include "compute_extended_bounds.h"

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
!
!
!  define inverse seawater density, use mean value for seawater density.

      rho0i=1.0/rho0
!
!  set inverse of specific heat for seawater (kg-degC/Joule).
!  cp is defined in scalars.h

      cpi=1.0/cp


      do j=jstrR,jendR
        do i=istrR,iendR


!  Input bulk parameterization fields

          wspd0=wspd(i,j)
          TairC=tair(i,j)
          TairK=TairC+273.16
# ifdef SST_SKIN
          TseaC=sst_skin(i,j)
# else
          TseaC=t(i,j,N,nrhs,itemp)
# endif
          TseaK=TseaC+273.16
          RH=rhum(i,j)
# ifndef BULK_LW
          hflw=-radlw(i,j) ! positive downward, i.e., warming the ocean
# else
!
!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
          hflw=radlw(i,j)  ! positive downward
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

          if (RH < 2.) then                          !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q=0.62197*(cff/(patm-0.378*cff))           !Spec hum (kg/kg)
          else          !RH input was actually specific humidity in g/kg
            Q=RH/1000.0                                !Spec Hum (kg/kg)
          endif

! water saturation vapor pressure (mb), using Teten formula.

          cff=(1.0007+3.46e-6*patm)*6.1121*
     &            exp(17.502*TseaC/(240.97+TseaC))

! vapor pressure reduced for salinity (Kraus & Businger, 1994, pp 42).

          cff=cff*0.98

! Qsea (kg/kg) from vapor pressure.

          Qsea=0.62197*(cff/(patm-0.378*cff))

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
          if (Ri<0.0) then
            Zetu=CC*Ri/(1.0+Ri/Ribcu)       ! Unstable
          else
            Zetu=CC*Ri/(1.0+3.0*Ri/CC)      ! Stable
          endif
          L10=blk_ZW/Zetu
          if (Zetu>50.0) then
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
          if (delW>18.0) then
            charn=0.018
          elseif ((10.0<delW).and.(delW<=18.0)) then
            charn=0.011+0.125*(0.018-0.011)*(delW-10.)
          else
            charn=0.011
          endif
!
!  Iterate until convergence. It usually converges within four
!  iterations.
!
        do iter=1,IterMax
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
          if (Bf>0.0) then
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
          stflx(i,j,isalt)=(evap-prate(i,j))*t(i,j,N,nrhs,isalt)
# endif

# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask_wet(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
# endif
!
!  Save for the calculation of wind stress
!
          aer(i,j)=rhoAir*wspd0*rho0i
          cer(i,j)=Cd
        enddo
      enddo

# ifdef QCORRECTION
      tmid=time+0.5*dt                    ! Flux correction for
      it7=3-itsst   ; it8=itsst           ! net surface heat flux
      cff7=sst_time(it8)-tmid             ! (relaxation toward
      cff8=tmid-sst_time(it7)             ! observed SST) and,...
      if (cff7 < dt) synchro_flag=.true.
      cff=1/(cff7+cff8)
      cff7=cff7*cff ; cff8=cff8*cff

      it9=3-itdqdt ; it10=itdqdt           ! optionally, same for
      cff9=dqdt_time(it10)-tmid            ! surface freshwater flux
      cff10=tmid-dqdt_time(it9)
      if (cff9 < dt) synchro_flag=.true.
      cff=1./(rho0*Cp*(cff9+cff10))
      cff9=cff9*cff ; cff10=cff10*cff

#  ifdef SFLX_CORR
      it11=3-itsss ; it12=itsss
      cff11=sss_time(it12)-tmid
      cff12=tmid-sss_time(it11)
      if (cff11 < dt) synchro_flag=.true.
      cff=1./(cff11+cff12)
      cff11=cff11*cff ; cff12=cff12*cff
#  endif
      do j=jstrR,jendR
        do i=istrR,iendR
          cff=cff9*dqdtg(i,j,it9)+cff10*dqdtg(i,j,it10)
          stflx(i,j,itemp)=stflx(i,j,itemp) +cff*(
     &         t(i,j,N,nrhs,itemp) -cff7*sstg(i,j,it7)
     &                             -cff8*sstg(i,j,it8)
     &                                              )
#  ifdef SFLX_CORR
          stflx(i,j,isalt)=stflx(i,j,isalt) +cff*(
     &         t(i,j,N,nrhs,isalt) -cff11*sssg(i,j,it11)
     &                             -cff12*sssg(i,j,it12)
     &                                              )
#  endif
        enddo
      enddo
# endif /* QCORRECTION */

!
!--------------------------------------------------------------
!  Compute kinematic, surface wind stress (m2/s2).
!--------------------------------------------------------------
!
# ifdef BULK_SM_UPDATE
!
!     sustr = rhoA/rhoW * Cd * wspd*uwnd  --> u point
!
      do j=jstrR,jendR
        do i=istrR+1,iendR
          a=0.5*(aer(i-1,j)+aer(i,j))
          cff=0.5*(cer(i-1,j)+cer(i,j))
          sustr(i,j)=a*cff*uwnd(i,j)
#  ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
#  endif
#  ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
#  endif
        enddo
      enddo
!
!     svstr = rhoA/rhoW * Cd * wspd*vwnd  --> v point
!
      do j=jstrR+1,jendR
        do i=istrR,iendR
          a=0.5*(aer(i,j-1)+aer(i,j))
          cff=0.5*(cer(i,j-1)+cer(i,j))
          svstr(i,j)=a*cff*vwnd(i,j)
#  ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
#  endif
#  ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
#  endif
        enddo
      enddo
# endif
      end




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
      end function bulk_psiu

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
      end function bulk_psit

#else
      subroutine bulk_flux_empty
      end
#endif
