#include "cppdefs.h"
#ifdef BIOLOGY

      subroutine biology_tile (istr,iend,jstr,jend)
!
! Compute biological forcing functions as defined by the Fasham et
! al. [JMR, 48, 591-639, 1990].  This routine was originated by John
! Moisan and adpapted for 3D code by MANU Sept. 8 98.   It computes
! r.h.s. terms associated with biological conversions. In this
! particular implementation we have: NO3, NH4, Detritus,
! PHYTOplankton and ZOOplanknton.
!
      implicit none
      integer istr,iend,jstr,jend
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "scalars.h"
      real solar, albedo, trans, PAR0_max, kwater, kphyto, alpha,
     &     K_NO3, K_NH4,  phi,   mu_40,    mu_43,  gmax,   K_Phyt,
     &     beta,  mu_30,  mu_32, mu_50,    mu_52,  mu_53,  pi,
     &            deg2rad,       ccf1,     ccf2,   ccf3,   ccf4
      parameter (
     &   solar   = 1353.,! the solar max is from Brock, 1981
     &   albedo  = 0.04, ! albedo of the ocean surface.
     &   trans   = 0.8,  ! fraction of total radiation
                         ! transmitted through atmosphere
!
! Potosynthetic Available Radiation (PAR) at the sea surface (without
! correction due to ellipticity of the Earth orbit, absorbtion by the
! effect of solar altitude and atmospheric clouds, see below; 0.43
! here is the estimated PAR fraction of the total solar radiation
!
     &   PAR0_max = solar*trans*0.43*(1.-albedo),
!
! Parameters as in Table 1; Fasham et al. [JMR, 48, 591-639, 1990]
!
     &   kwater = 0.04, ! light attenuation due to sea water [m-1]
     &   kphyto = 0.03, ! light attenuation by Phytoplankton
                        !                           [(m^2 mMol N)-1]
     &   alpha  = 0.025,! initial slope of the P-I curve
                        !                            [(W m-2)-1 d-1]
     &   K_NO3  = 0.5,  ! half-saturation for Phytoplankton NO3
                        !                        uptake [mMol N m-3]
     &   K_NH4  = 0.5,  ! half-saturation for Phytoplankton NH4
                        !                        uptake [mMol N m-3]
     &   phi    = 1.5,  ! Phytoplankton ammonium inhibition
                        !                     parameter [(mMol N)-1]

     &   mu_40  = 0.018,! Phyto loss to sink rate[d-1]
     &   mu_43  = 0.072,! Phyto mortality to Detritus rate d-1]
     &   gmax   = 0.75, ! maximum Zooplankton growth rate [d-1]
     &   beta   = 0.75, ! Zooplankton assimilation efficiency of
                        !                         Zooplankton [n.d.]
     &   K_Phyt = 1.0,  ! Zooplankton half-saturation conts. for
     &                  !                            ingestion [d-1]
     &   mu_50  = 0.025,! Zooplankton loss to sink [d-1]
     &   mu_52  = 0.1,  ! Zooplankton specific excretion rate [d-1]
     &   mu_53  = 0.025,! Zooplankton mortality to Detritus [d-1]

     &   mu_30  = 0.02, ! Detrital loss to sink rate [d-1]
     &   mu_32  = 0.03, ! Detrital breakdown to NH4 rate [d-1]

     &   ccf1   = 0.6071538329,      ! Set OSW Papa CLOUD
     &   ccf2   = 1.187075734,       ! correction coefficients
     &   ccf3   = 0.7726212144,
     &   ccf4   =-0.2782480419,

     &   pi     = 3.14159265358979323846,
     &   deg2rad= 2.*pi/360.)

      integer i,j,k, ITER
      real dt_bio,  day,  cos_Day, PAR0,     NO3(N),  NO3_bak(N),
     &     cff,     dec,  cos_Thr, PAR0_ell, NH4(N),  NH4_bak(N),
     &     cloud,         cos_Znt, PARz,     Det(N),  Det_bak(N),
     &     daym(14),      cos_dec, Vp,       Phyt(N), Phyt_bak(N),
     &     coktas(14),    sin_dec, aJ(N),    Zoo(N),  Zoo_bak(N)

# define VERTICAL_SINKING
# ifdef VERTICAL_SINKING
      real wDet,wPhyt,wZoo
      parameter (wDet=10., wPhyt=1.,wZoo=1.)

      real FC(0:N),dDet(0:N),dPhyt(0:N),dZoo(0:N),
     &         aL,aR, cff1,cff2,cff6,cu
# endif
!
# include "compute_auxiliary_bounds.h"
!
      daym(1)=0.                 ! Days of the year
      daym(2)=16.
      daym(3)=46.
      daym(4)=75.
      daym(5)=105.
      daym(6)=136.
      daym(7)=166.
      daym(8)=197.
      daym(9)=228.
      daym(10)=258.
      daym(11)=289.
      daym(12)=319.
      daym(13)=350.
      daym(14)=360.

      coktas(1)=6.29             ! OWS Papa CLOUD
      coktas(2)=6.26             ! CLIMATOLOGY
      coktas(3)=6.31
      coktas(4)=6.31
      coktas(5)=6.32
      coktas(6)=6.70
      coktas(7)=7.12
      coktas(8)=7.26
      coktas(9)=6.93
      coktas(10)=6.25
      coktas(11)=6.19
      coktas(12)=6.23
      coktas(13)=6.31
      coktas(14)=6.29

      day=tdays-int(tdays/360.)*360.  ! convert in days 1-360

!
! Coefficient for the cloud correction algorithm from fitting
! to Bishop and Rossow data
!
      do j=1,13
       if (day>=daym(j) .and. day<=daym(j+1)) i=j
      enddo
      cloud=0.125*(coktas(i)*(daym(i+1)-day)
     &             +coktas(i+1)*(day-daym(i))
     &                  )/(daym(i+1)-daym(i))
!
! Calculate PAR0 by vertical integration assume the surface chlor
! is a good approx of chlor with depth [this is a rough initial est.]
!
      cos_Day=cos(deg2rad*day)                   ! seasonal cycle
      cos_Thr=cos(2.*pi*(tdays-int(tdays)-0.5))  ! dayly cycle

c??
      cos_Day=0.         ! replace astronomical parameters
      cos_Thr=1./pi      ! with their effective average values.
      cloud=0.125*6.25
c??

      dec=-0.406*cos_Day !in radians and from Evans and Parlsow, 1985
      cos_dec=cos(dec)
      sin_dec=sin(dec)
!
! Correction of solar constant due to the ellipticity of the Earth
! orbit, after [Duffie and Beckman, 1980]
!
      PAR0_ell=PAR0_max*(1.+0.033*cos_Day)

      dt_bio=dt/(24.*3600.)  ! time step as fraction of day.

!
! Since the following solver is iterative to achieve implicit
! discretization of the biological interaction, two time slices are
! required, BIO and BIO_bak, where BIO is understood as vector of
! biological state variables: BIO=[NO3,NH4,Det,Phyt,Zoo]. Assume
! that the iterations converge, the newly obtained state variables
! satisfy equations
!
!           BIO = BIO_bak + dt_bio * rhs(BIO)
!
! where rhs(BIO) is the vector of biological r.h.s. computed at
! the new time step. During the iterative procedure a series of
! fractional time steps is performed in a chained mode (splitting
! by different biological conversion processes) in sequence NO3 --
! NH4 -- Phyt -- Zoo -- Det, that is the main food chain. In all
! stages the concentration of the component being consumed is
! treated in fully implicit manner, so that the algorithm guarantees
! non-negative values, no matter how strong is the concentration of
! active consuming component (Phyto or Zoo).
!
! The overall algorithm, as well as any stage of it is formulated
! in conservative form (except explicit sinking) in sense that the
! sum of concentration of all five components is conserved.
!

#  ifdef EW_PERIODIC
#   define I_RANGE istr,iend
#  else
#   define I_RANGE istrR,iendR
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE jstr,jend
#  else
#   define J_RANGE jstrR,jendR
#  endif

      do j=J_RANGE
        do i=I_RANGE
!
! Extract biological variables from tracer arrays; place them into
! scratch variables; restrict their values to be positive definite.
!
          do k=1,N
            NO3_bak(k) =max(t(i,j,k,nnew,iNO3_),0.) ! Nitrate
            NH4_bak(k) =max(t(i,j,k,nnew,iNH4_),0.) ! Ammonium
            Det_bak(k) =max(t(i,j,k,nnew,iDet_),0.) ! Detritus
            Phyt_bak(k)=max(t(i,j,k,nnew,iPhyt),0.) ! Phytoplankton
            Zoo_bak(k) =max(t(i,j,k,nnew,iZoo_),0.) ! Zooplankton

            NO3(k)  = NO3_bak(k)
            NH4(k)  = NH4_bak(k)
            Det(k)  = Det_bak(k)
            Phyt(k) = Phyt_bak(k)
            Zoo(k)  = Zoo_bak(k)
          enddo
!
! Calulate aJ (here: cos_Znt -- cos of solar zenith angle)
!
          cff=deg2rad*latr(i,j)
          cos_Znt=cos_Thr*cos_dec*cos(cff)+sin_dec*sin(cff)
          if (cos_Znt>0.) then
            PAR0=PAR0_ell*cos_Znt*(1.-ccf1+ccf2*cos_Znt)
     &        *(1.-cloud*( ccf3+ccf4*sqrt(1.-cos_Znt*cos_Znt)))
            do k=1,N                              ! From Eppley, d-1:
              Vp=0.851*1.066**t(i,j,k,nnew,itemp) !   Vp=2.9124317 at
                                                  !   t=19.25 degrees
              PARz=PAR0*exp(-abs(z_r(i,j,k))*(kwater+kphyto*Phyt(k)))
              aJ(k)=Vp*alpha*PARz/sqrt(Vp*Vp+alpha*alpha*PARz*PARz)
            enddo
          else
            do k=1,N
              aJ(k)=0.           ! <-- during the night
            enddo
          endif

          DO ITER=1,3   !--> Start internal iterations to achieve
!                            nonlinear backward-implicit solution.
! NO3 uptake by Phyto                                           [1-4]
!
            do k=1,N
              cff=dt_bio*Phyt(k)*aJ(k)*exp(-phi*NH4(k))/(K_NO3+NO3(k))
              NO3(k)=NO3_bak(k)/(1.+cff)
              Phyt(k)=Phyt_bak(k)+cff*NO3(k)
            enddo
!
! NH4 uptake by Phyto                                           [2-4]
!
            do k=1,N
              cff=dt_bio*Phyt(k)*aJ(k)/(K_NH4+NH4(k))
              NH4(k)=NH4_bak(k)/(1.+cff)
              Phyt(k)=Phyt(k)+cff*NH4(k)
            enddo
!
! Phytoplankton grazing by Zooplankton;                         [4-5]
!               mortality to Detritus (at rate mu_43);          [4-3]
!                        loss to sink (at rate mu_40);          [4-0]
!
            do k=1,N
              cff=dt_bio*gmax*Zoo(k)/(K_Phyt+Phyt(k))

c              cff=dt_bio*gmax*Zoo(k)/(K_Phyt        )

              Phyt(k)=Phyt(k)/(1.+cff+dt_bio*(mu_40+mu_43))
c_no_sink            Phyt(k)=Phyt(k)/(1.+cff+dt_bio*mu_43)

              Zoo(k)=Zoo_bak(k)+Phyt(k)*cff*beta
              Det(k)=Det_bak(k)+Phyt(k)*(cff*(1.-beta)+dt_bio*mu_43)
            enddo
!
! Zooplankton mortality to sinking (rate mu_50);                [5-0]
!                           to Det (rate mu_53);                [5-3]
!                 excretion to NH4 (rate mu_52)                 [5-2]
!
            do k=1,N
              Zoo(k)=Zoo(k)/(1.+dt_bio*(mu_50+mu_52+mu_53))
c_no_sink           Zoo(k)=Zoo(k)/(1.+dt_bio*(mu_52+mu_53))

              NH4(k)=NH4(k)+dt_bio*mu_52*Zoo(k)
              Det(k)=Det(k)+dt_bio*mu_53*Zoo(k)
            enddo
!
! Detritus breakdown to NH4  (at rate mu_32)                    [3-2]
!               loss to sink (at rate mu_30)                    [3-0]
!
            do k=1,N
              Det(k)=Det(k)/(1.+dt_bio*(mu_30+mu_32))
c_no_sink             Det(k)=Det(k)/(1.+dt_bio*mu_32)

              NH4(k)=NH4(k)+dt_bio*mu_32*Det(k)
            enddo

# ifdef VERTICAL_SINKING
!
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
!
! Part (i): Construct parabilic splines: compute vertical derivatives
! of the felds Det,Phyt,Zoo. The derivatives are placed at W-points;
! Neumann boundary conditions are assumed on top and bottom.
!
      FC(0)=0.
      dDet(0)=0.
      dPhyt(0)=0.
      dZoo(0)=0.

      cff6=6./ds
      do k=1,N-1,+1
        cff=1./(2.*Hz(i,j,k+1)+Hz(i,j,k)*(2.-FC(k-1)))
        FC(k)=cff*Hz(i,j,k+1)
        dDet(k) =cff*(cff6*(Det(k+1) -Det(k)) -Hz(i,j,k)*dDet(k-1))
        dPhyt(k)=cff*(cff6*(Phyt(k+1)-Phyt(k))-Hz(i,j,k)*dPhyt(k-1))
        dZoo(k) =cff*(cff6*(Zoo(k+1) -Zoo(k)) -Hz(i,j,k)*dZoo(k-1))
      enddo

      dDet(N)=0.
      dPhyt(N)=0.
      dZoo(N)=0.

      do k=N-1,1,-1
        dDet(k) =dDet(k) -FC(k)*dDet(k+1)
        dPhyt(k)=dPhyt(k)-FC(k)*dPhyt(k+1)
        dZoo(k) =dZoo(k) -FC(k)*dZoo(k+1)
      enddo
!
! Part (ii): Convert dDet,dPhyt,dZoo [which are now derivatives of
! fields Det,Phyt,Zoo at the grid box interfaces] into field values
! at these interfaces, assuming parabolic profiles within each grid
! box. Restrict these values to lie between bounds determined from
! box-averaged values of grid boxes adjscent from above and below.
! (This restriction is part of PPM-like monotonization procedure.)
!
      cff=ds/3
      dDet(0) =Det(1)  ! -cff*Hz(i,j,1)*(dDet(0) +0.5*dDet(1))
      dPhyt(0)=Phyt(1) ! -cff*Hz(i,j,1)*(dPhyt(0)+0.5*dPhyt(1))
      dZoo(0) =Zoo(1)  ! -cff*Hz(i,j,1)*(dZoo(0) +0.5*dZoo(1))

      dDet(N) =Det(N)  ! +cff*Hz(i,j,N)*(dDet(N) +0.5*dDet(N-1))
      dPhyt(N)=Phyt(N) ! +cff*Hz(i,j,N)*(dPhyt(N)+0.5*dPhyt(N-1))
      dZoo(N) =Zoo(N)  ! +cff*Hz(i,j,N)*(dZoo(N) +0.5*dZoo(N-1))

      do k=2,N,+1        !<-- irreversible
        dDet(k-1) =Det(k) -cff*Hz(i,j,k)*(0.5*dDet(k)+dDet(k-1))
        dDet(k-1) =max(dDet(k-1), min(Det(k-1), Det(k)))
        dDet(k-1) =min(dDet(k-1), max(Det(k-1), Det(k)))

        dPhyt(k-1)=Phyt(k)-cff*Hz(i,j,k)*(0.5*dPhyt(k)+dPhyt(k-1))
        dPhyt(k-1)=max(dPhyt(k-1),min(Phyt(k-1),Phyt(k)))
        dPhyt(k-1)=min(dPhyt(k-1),max(Phyt(k-1),Phyt(k)))

        dZoo(k-1) =Zoo(k) -cff*Hz(i,j,k)*(0.5*dZoo(k)+dZoo(k-1))
        dZoo(k-1) =max(dZoo(k-1), min(Zoo(k-1), Zoo(k)))
        dZoo(k-1) =min(dZoo(k-1), max(Zoo(k-1), Zoo(k)))
      enddo
!
! Part (iii): Conver dDet,dPhyt,dZoo into flux-integrated values,
! complete PPM flux limiting. This procedure starts from assigning
! Left and Right (aR,aL) values of the interpolating parabolae, then
! monotonicity conditions are checked and aL,aR are modified to fit.
! Overall, it follows Colella-Woodward, 1984 almost exactly.
!
      do k=1,N,+1        !<-- irreversible
        FC(k)=dt_bio/(ds*Hz(i,j,k))

        aR=dDet(k)
        aL=dDet(k-1)
        cff1=(aR-aL)*6.*(Det(k)-0.5*(aR+aL))
        cff2=(aR-aL)**2
        if ((aR-Det(k))*(Det(k)-aL)<0.) then
          aL=Det(k)
          aR=Det(k)
        elseif (cff1>cff2) then
          aL=3.*Det(k)-2.*aR
        elseif (cff1<-cff2) then
          aR=3.*Det(k)-2.*aL
        endif
        cu=wDet*FC(k)
        dDet(k-1)=Det(k) - (1.-cu)*(0.5*(aR-aL)-( 0.5*(aR+aL)-
     &                                       Det(k))*(1.-2.*cu))

        aR=dPhyt(k)
        aL=dPhyt(k-1)
        cff1=(aR-aL)*6.*(Phyt(k)-0.5*(aR+aL))
        cff2=(aR-aL)**2
        if ((aR-Phyt(k))*(Phyt(k)-aL)<0.) then
          aL=Phyt(k)
          aR=Phyt(k)
        elseif (cff1>cff2) then
          aL=3.*Phyt(k)-2.*aR
        elseif (cff1<-cff2) then
          aR=3.*Phyt(k)-2.*aL
        endif
        cu=wPhyt*FC(k)
        dPhyt(k-1)=Phyt(k) - (1.-cu)*( 0.5*(aR-aL)-(0.5*(aR+aL)-
     &                                        Phyt(k))*(1.-2.*cu))

        aR=dZoo(k)
        aL=dZoo(k-1)
        cff1=(aR-aL)*6.*(Zoo(k)-0.5*(aR+aL))
        cff2=(aR-aL)**2
        if ((aR-Zoo(k))*(Zoo(k)-aL)<0.) then
          aL=Zoo(k)
          aR=Zoo(k)
        elseif (cff1>cff2) then
          aL=3.*Zoo(k)-2.*aR
        elseif (cff1<-cff2) then
          aR=3.*Zoo(k)-2.*aL
        endif
        cu=wZoo*FC(k)
        dZoo(k-1)=Zoo(k) - (1.-cu)*(0.5*(aR-aL)-( 0.5*(aR+aL)-
     &                                       Zoo(k))*(1.-2.*cu))
      enddo

      dDet(N)=0.   ! Set no-flux boundary conditions at top.
      dPhyt(N)=0.
      dZoo(N)=0.
!
! Apply fluxes:
!
      do k=1,N
        Det(k) =Det(k)  + wDet *FC(k)*(dDet(k) -dDet(k-1))
        Phyt(k)=Phyt(k) + wPhyt*FC(k)*(dPhyt(k)-dPhyt(k-1))
        Zoo(k) =Zoo(k)  + wZoo *FC(k)*(dZoo(k) -dZoo(k-1))
      enddo
# endif    /* VERTICAL_SINKING */

      ENDDO  ! <-- ITER

c          cff=NO3(1)+NH4(1)+Det(1)+Phyt(1)+Zoo(1)
c          write(*,'(F7.2,6F12.9)') tdays, NO3,NH4,Det,Phyt,Zoo,cff
!
! Write back
!
          do k=1,N
           t(i,j,k,nnew,iNO3_)=t(i,j,k,nnew,iNO3_)-NO3_bak(k) +NO3(k)
           t(i,j,k,nnew,iNH4_)=t(i,j,k,nnew,iNH4_)-NH4_bak(k) +NH4(k)
           t(i,j,k,nnew,iDet_)=t(i,j,k,nnew,iDet_)-Det_bak(k) +Det(k)
           t(i,j,k,nnew,iPhyt)=t(i,j,k,nnew,iPhyt)-Phyt_bak(k)+Phyt(k)
           t(i,j,k,nnew,iZoo_)=t(i,j,k,nnew,iZoo_)-Zoo_bak(k) +Zoo(k)
          enddo
        enddo      ! <-- i
      enddo      ! <-- j
#else
      subroutine biology_empty ()
#endif
      end


