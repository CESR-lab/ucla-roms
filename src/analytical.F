#include "cppdefs.opt"

! ANALYTICAL PACKAGE:   Set up various analytical forcing fields
!=========== ========   for the model.
!
! ana_bsedim_tile   Analytical bottom sediment grain size
!                          and density.
! ana_meanRHO_tile  Analytical mean density anomaly.
! ana_smflux_tile   Analytical kinematic surface momentum flux
!                          (wind stress).
! ana_srflux_tile   Analytical kinematic surface shortwave
!                          radiation.
! ana_ssh_tile      Analytical sea surface height climatology.
! ana_sst_tile      Analytical sea surface temperature and dQdSST
!                         which are used during heat flux correction.
! ana_stflux_tile   Analytical kinematic surface flux of tracer type
!                          variables.
! ana_tclima_tile   Analytical tracer climatology fields.
! ana_uclima_tile   Analytical tracer climatology fields.
! ana_wwave_tile    Analytical wind induced wave amplitude,
!                         direction and period.

#if defined ANA_BSEDIM && defined SG_BBL96
      subroutine ana_bsedim_tile(istr,iend,jstr,jend)

! Set initial bottom sediment grain diameter size [m]  and density
! used in the bottom boundary formulation [kg/m^3].

      use param
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j

# include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          Ssize(i,j)=???
          Sdens(i,j)=???
        enddo
      enddo
      end
#endif /* ANA_BSEDIM && SG_BBL96 */







#if defined SOLVE3D && defined ANA_MEANRHO
      subroutine ana_meanRHO_tile(istr,iend,jstr,jend)

!  Set mean density anomaly [kg/m^3] using an analytical expression.
      use param
      use grid
      use ocean3d
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real cff1, cff2 , cff3, cff4

# include "compute_extended_bounds.h"

!  Set mean density anomaly (kg/m^3) at horizontal and vertical
!  RHO-points.

# if defined CANYON_A || defined GRAV_ADJ || defined OVERFLOW
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=R0
          enddo
        enddo
      enddo
# elif defined CANYON_B
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=R0-3.488*exp(z_w(i,j,k)/800.)*
     &                  (1.-(2./3.)*tanh(z_r(i,j,k)/800.))
          enddo
        enddo
      enddo
# elif defined SEAMOUNT
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=28.-2.*exp(z_r(i,j,k)/1000.)
          enddo
        enddo
      enddo
# elif defined SHELFRONT
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=R0
          enddo
        enddo
      enddo
# elif defined UPWELLING
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=26.459472-2.24*exp(z_r(i,j,k)/50.)
          enddo
        enddo
      enddo
# elif defined JET
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
           cff1=TL1 + TL2*z_r(i,j,k) + TL3*exp(z_r(i,j,k)/ZL1)
           cff2=TC1 + TC2*z_r(i,j,k) + TC3*exp((z_r(i,j,k)-ZC1)/ZC2)
           cff3=YC1 +YC2* tanh(YC3*xr(i,j)+YC4)

!          cff4=0.5*( cff2*(1+tanh((cff3-yr(i,j))/LCUR))
!                    +cff1*(1+tanh((yr(i,j)-cff3)/LCUR)))
           cff4=cff1
            rhobar(i,j,k)=R0+Tcoef*cff4
          enddo
        enddo
      enddo
# else
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=???
          enddo
        enddo
      enddo
# endif
# ifdef MASKING
      do k=1,N                  ! Apply Land/Sea mask.
        do j=jstrR,jendR
          do i=istrR,iendR
            rhobar(i,j,k)=rhobar(i,j,k)*rmask(i,j)
          enddo
        enddo
      enddo
# endif
      end
#endif /* ANA_MEANRHO && SOLVE3D */


#ifdef ANA_SMFLUX
      subroutine ana_smflux_tile(istr,iend,jstr,jend)

      use param
      use surf_flux, only: sustr,svstr
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
      real TauX, TauY, cff

# include "compute_extended_bounds.h"

      if (FIRST_TIME_STEP) then  ! Set wind stress (kinematic surface
        do j=jstrR,jendR         ! momentum flux [m^2/s^2]) components
          do i=istrR,iendR       ! in XI- and ETA-directions [m^2/s^2])
            sustr(i,j)=0.        ! at horizontal U- and V-points.
            svstr(i,j)=0.        ! For code compactness, both are set
          enddo                  ! to zero at first (default values),
        enddo                    ! then either one of them, or both may
      endif                      ! be overwritten by nontrivial values.

# ifdef DOUBLE_GYRE
      if (FIRST_TIME_STEP) then    ! Wind is constant in time, so
        TauX=(0.1/1000.)           ! initialize it during the first
        cff=2.*pi/el               ! time step, and use thereafter. 
        do j=jstrR,jendR
          sustr(istr,j)=TauX*(1.-cos(cff*yr(istr,j)))
          do i=istr+1,iendR
            sustr(i,j)=sustr(istr,j)
          enddo
        enddo
      endif
# elif defined CANYON_A || defined CANYON_B
      do j=jstrR,jendR
        do i=istr,iendR
          sustr(i,j)=0.0001*0.5*sin(2.*pi*tdays/10.)*(
     &               1.-tanh((yr(i,j)-0.5*el)/10000.) )
        enddo
      enddo
# elif defined EKMAN_SPIRAL
      TauX=0.1/rho0
      do j=jstrR,jendR
        do i=istrR,iendR
          sustr(i,j)=TauX
        enddo
      enddo

# elif defined COLD_FILAMENT
      TauY=0.1/rho0   ;  TauX=-1.0*TauY ! This setup allows rotating
      cff=TauY/sqrt(TauY**2+TauX**2)    ! wind direction by setting the
      TauX=cff*TauX   ;  TauY=cff*TauY  ! initial value of TauX, then
      do j=jstrR,jendR                  ! adjacting both components to
        do i=istrR,iendR                ! kee the overall amplitute as
          sustr(i,j)=TauX               ! initially specified by TauY.
          svstr(i,j)=TauY
        enddo
      enddo

# elif defined UPWELLING
      if (tdays < 2.) then
        TauX=+0.1*sin(pi*tdays/4.)/rho0
      else
        TauX=+0.1/rho0
      endif

      do j=jstrR,jendR
        do i=istrR,iendR
          sustr(i,j)=TauX
c?        svstr(i,j)=0.2*TauX !<-- to make it ETALON_CHECK sensitive
        enddo
      enddo
# elif defined BALTIC

      do j=jstrR,jendR
        do i=istrR,iendR
          sustr(i,j)=0.1/rho0  ! Wind blowing West --> East
        enddo
      enddo


# elif defined USWEST
      do j=jstrR,jendR
        do i=istr,iendR
c         cff1=(latr(i,j)-latr(Lm/2,Mm/2))/20.
c         sustr(i,j)=1.d-4.*cff1
c         sustr(i,j)=-1.d-4
          sustr(i,j)=0.
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
c          svstr(i,j)=-1.0d-4
          svstr(i,j)=0.
        enddo
      enddo
# endif
      end
#endif /* ANA_SMFLUX */

#if defined ANA_SRFLUX && defined SOLVE3D

      subroutine ana_srflux_tile(istr,iend,jstr,jend)

! Set kinematic surface solar shortwave radiation flux "srflx"
! [degC m/s] using an analytical expression.
      use param
      use surf_flux, only: srflx
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j

# include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          srflx(i,j)=0.
        enddo
      enddo
      end
#endif /* ANA_SRFLUX */


#if defined ANA_SSH
      subroutine ana_ssh(tile)            ! Set analytical sea-surface
      use param
      use hidden_mpi_vars
      implicit none                       ! height climatology [m]
      integer tile
#include "compute_tile_bounds.h"
      call ana_ssh_tile(istr,iend,jstr,jend)
      end

      subroutine ana_ssh_tile(istr,iend,jstr,jend)

      use param
      use climat
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
#ifdef JET
      real TL0, TC0, cff1, cff2, cff3, cff4
#endif

# include "compute_extended_bounds.h"

# if defined JET
      TL0=-TL1 - 0.5*TL2*ZREF**2 - ZL1*TL3*exp(ZREF/ZL1)
      TC0=-TC1 - 0.5*TC2*ZREF**2 - ZC2*TC3*exp((ZREF-ZC1)/ZC2)
      do j=jstrR,jendR
        do i=istrR,iendR
         cff1=YC1+YC2*tanh(YC3*xr(i,j)+YC4)
         cff2=TL0+ZL1*TL3
         cff3=TC0+ZC2*TC3*exp(-ZC1/ZC2)
         cff4=-0.5*g*Tcoef*(cff2+cff3)
         ssh(i,j)=cff4/(rho0*g)
                +Tcoef*cff2*(1+tanh((cff1-yr(i,j))/LCUR))/(2*rho0)
                +Tcoef*cff3*(1+tanh((yr(i,j)-cff1)/LCUR))/(2*rho0)
        enddo
      enddo
# else
      do j=jstrR,jendR
        do i=istrR,iendR
          ssh(i,j)=???
        enddo
      enddo
#endif /* JET */
      end
#endif



#ifdef SOLVE3D
# if defined ANA_STFLUX || defined ANA_SSFLUX
      subroutine ana_stflux_tile(istr,iend,jstr,jend, itrc)

! Set kinematic surface flux of tracer type variables "stflx" (tracer
! units m/s) using analytical expressions.
      use param
      use surf_flux, only: stflx
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, itrc, i,j

# include "compute_extended_bounds.h"

      if (itrc==itemp) then

! Set kinematic surface heat flux [degC m/s] at horizontal
! RHO-points.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,itemp)=0.
          enddo
        enddo

#  ifdef SALINITY
      elseif (itrc == isalt) then

! Set kinematic surface freshwater flux (m/s) at horizontal
! RHO-points, scaling by surface salinity is done in STEP3D.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,isalt)=0.
          enddo
        enddo
#  endif
      else

!  Set kinematic surface flux of additional tracers, if any.

      endif
      end
# endif /* ANA_STFLUX || ANA_SSFLUX */



# if defined ANA_TCLIMA && defined TCLIMATOLOGY
      subroutine ana_tclima(tile)
      use param
      use hidden_mpi_vars
      implicit none                   ! Set up analytical tracer
      integer tile                    ! climatology fields.
# include "compute_tile_bounds.h"
      call ana_tclima_tile(istr,iend,jstr,jend)
      end

      subroutine ana_tclima_tile(istr,iend,jstr,jend)

      use param
      use climat
      use grid
      use ocean3d
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j,k
#  ifdef JET
      real cff1, cff2, cff3
#  endif

#  include "compute_extended_bounds.h"

#  if defined JET
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
           cff1=TL1+ TL2*z_r(i,j,k) +TL3*exp(z_r(i,j,k)/ZL1)
           cff2=TC1+ TC2*z_r(i,j,k) +TC3*exp((z_r(i,j,k)-ZC1)/ZC2)
           cff3=YC1+YC2*tanh(YC3*xr(i,j)+YC4)
           tclm(i,j,k,itemp)=0.5*(cff2*(1+tanh((cff3-yr(i,j))/LCUR))
     &                           +cff1*(1+tanh((yr(i,j)-cff3)/LCUR)))
          enddo
        enddo
      enddo

#  elif defined ISWAKE

! There is no code here because for this case climatological data is
! saved from the initial conditions and thereafter remaines unchanged.

#  else
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            tclm(i,j,k,itemp)=???
#   ifdef SALINITY
            tclm(i,j,k,isalt)=???
#   endif
          enddo
        enddo
      enddo
#  endif
      end
# endif /* ANA_TCLIMA && TCLIMATOLOGY */
#endif /* SOLVE3D */


#if defined ANA_UCLIMA && defined UCLIMATOLOGY
      subroutine ana_uclima(tile)
      use param
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer tile

#include "compute_tile_bounds.h"
      call ana_uclima_tile   (istr,iend,jstr,jend)!,A2d(1,1),A2d(1,2))
      end

      subroutine ana_uclima_tile(istr,iend,jstr,jend)!, wrk1,wrk2)

!---------------------------------------------------------------------
!  This routine sets analytical momentum climatology fields.
!---------------------------------------------------------------------

      use param
      use climat
      use grid
      use ocean3d

      implicit none
      integer istr,iend,jstr,jend, i,j,k
c      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk1, wrk2
# ifdef JET
      real TL0, TC0, cff1, cff2, cff3, cff4, cff5
# endif

# include "compute_extended_bounds.h"

! Set climatology fields for tracer type variables.
!---- ----------- ------ --- ------ ---- ----------

# if defined USWEST
      do j=jstrR,jendR
        do i=istrR,iendR
          ubclm(i,j)=0.
          vbclm(i,j)=0.
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            uclm(i,j,k)=0.
            vclm(i,j,k)=0.
          enddo
        enddo
      enddo
#  endif
# elif defined JET
      TL0=-TL1 - 0.5*TL2*ZREF**2 - ZL1*TL3*exp(ZREF/ZL1)
      TC0=-TC1 - 0.5*TC2*ZREF**2 - ZC2*TC3*exp((ZREF-ZC1)/ZC2)
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          cff1=TL0*h(i,j) - 0.5*TL1*h(i,j)**2 + TL2*(h(i,j)**3)/6.
     &                         - (ZL1**2)*TL3*(1-exp(-h(i,j)/ZL1))
          cff2=TC0*h(i,j) - 0.5*TC1*h(i,j)**2 + TC2*(h(i,j)**3)/6.
     &           - (ZC2**2)*TC3*exp(-ZC1/ZC2)*(1-exp(-h(i,j)/ZC2))
          cff3=YC1+YC2*tanh(YC3*xr(i,j)+YC4)
          cff4=-g*Tcoef*(1-(tanh((cff3-yr(i,j))/LCUR))**2)
     &                       /(2*rho0*f(i,j)*LCUR*h(i,j))
! ubar
          wrk1(i,j)=cff4*(cff1-cff2)
! vbar
          cff5=YC2*YC3*(1-(tanh(YC3*xr(i,j)+YC4))**2)
          wrk2(i,j)=cff4*cff5*(cff1-cff2)
        enddo
      enddo
      do j=jstrR,jendR
        do i=istr,iendR
          ubclm(i,j)=0.5*(wrk1(i,j)+wrk1(i-1,j))
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          vbclm(i,j)=0.5*(wrk2(i,j)+wrk2(i,j-1))
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            cff1=TL0 + TL1*z_r(i,j,k) + 0.5*TL2*(z_r(i,j,k)**2)
     &                            + ZL1*TL3*exp(z_r(i,j,k)/ZL1)

            cff2=TC0 + TC1*z_r(i,j,k) + 0.5*TC2*(z_r(i,j,k)**2)
     &                      + ZC2*TC3*exp((z_r(i,j,k)-ZC1)/ZC2)
            cff3=YC1 + YC2*tanh(YC3*xr(i,j)+YC4)
            cff4=g*Tcoef*(1.-(tanh((cff3-yr(i,j))/LCUR))**2)
     &                                 /(2*rho0*f(i,j)*LCUR)
! u
            wrk1(i,j)=cff4*(cff1-cff2)
! v
            cff5=YC2*YC3*(1.-(tanh(YC3*xr(i,j)+YC4))**2)
            wrk2(i,j)=cff4*cff5*(cff1-cff2)
          enddo
        enddo
        do j=jstrR,jendR
          do i=istr,iendR
            uclm(i,j,k)=0.5*(wrk1(i,j)+wrk1(i-1,j))
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
            vclm(i,j,k)=0.5*(wrk2(i,j)+wrk2(i,j-1))
          enddo
        enddo
      enddo
#  endif
# elif defined ISWAKE

! There is no code here because for this case climatological data is
! saved from the initial conditions and thereafter remaines unchanged.

# else
      do j=jstrR,jendR
        do i=istrR,iendR
          ubclm(i,j)=???
          vbclm(i,j)=???
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            uclm(i,j,k)=???
            vclm(i,j,k)=???
          enddo
        enddo
      enddo
#  endif
# endif
      end
#endif /* ANA_UCLIMA && UCLIMATOLOGY */
!
#if defined ANA_WWAVE && defined SG_BBL96
      subroutine ana_wwave_tile(istr,iend,jstr,jend)

! Set wind induced wave amplitude, direction and period for
! bottom boundary layer formulation.
      use param
      use grid

      implicit none
      integer istr,iend,jstr,jend, i,j

# include "compute_extended_bounds.h"

!  Set wind induced wave amplitude (m), direction (radians) and
!  period (s) at RHO-points.

      do j=jstrR,jendR
        do i=istrR,iendR
          Awave(i,j)=???
          Dwave(i,j)=???
          Pwave(i,j)=???
        enddo
      enddo
      end
#endif /* ANA_WWAVE && SG_BBL96 */

#if defined PSOURCE && defined ANA_PSOURCE
                              ! Set analytical tracer and mass point
      subroutine ana_psource  ! sources and sinks. River runoff can
      use param
      use scalars
      implicit none           ! be consider as a point source.
!
      integer is, k
      real cff, fac

! Set-up point Sources/Sink number (Nsrc), direction (Dsrc), I- and
! J-grid locations (Isrc,Jsrc).  Currently, the direction can be along
! XI- (Dsrc = 0) or along ETA-direction (Dsrc > 0).  The mass sources
! are located at U- or V-points so the grid locations should range
! from 1 =< Isrc =< L  and  1 =< Jsrc =< M.

      if (FIRST_TIME_STEP) then
# if defined IBERIA
        Nsrc=4
        Dsrc(1)=1            ! Minho
        Isrc(1)=109
        Jsrc(1)=137

        Dsrc(2)=1            ! Douro
        Isrc(2)=116
        Jsrc(2)=114

        Dsrc(3)=1            ! Mondego
        Isrc(3)=107
        Jsrc(3)=78

        Dsrc(4)=1            ! TEJO
        Isrc(4)=102
        Jsrc(4)=36
#  ifdef SOLVE3D
        do is=1,Nsrc                   ! Set-up nondimensional shape
          do k=1,9                     ! function to distribute mass
            Qshape(is,k)=0.1           ! point sources/sinks
          enddo                        ! vertically.
          Qshape(is,10)=0.06           ! These most add up to unity!
          Qshape(is,11)=0.03
          Qshape(is,12)=0.01
          do k=13,19
            Qshape(is,k)=0.0
          enddo
        enddo
#  endif

# elif defined RIVER
        Nsrc=3
        do is=1,Nsrc
          Dsrc(is)=0
          Isrc(is)=1
          Jsrc(is)=3*MMm/4 +is-1
        enddo
#  ifdef SOLVE3D
        cff=1./float(N)
        do k=1,N                         ! Uniform vertical
          do is=1,Nsrc                   ! distribution
            Qshape(is,k)=cff
          enddo
        enddo
#  endif
# else
#  error  "### ERROR:  CPP-key 'ANA_PSOURCE' is defined, but no code"
#  error  "### ERROR:  is provided to set up Dsrc, Isrc, Jsrc"
           xxx <-- to stop compilation here
# endif
      endif
!
! Set-up vertically integrated mass transport [m3/s] of point
! sources (these may time-dependent; positive in the positive U-
! or V-direction and vice-versa) and vertically distribute them
! according to mass transport profile chosen above.
!
# ifdef SOLVE3D
#  ifdef IBERIA
!     if (tdays<0.5) then
!       fac=1.0+TANH((time-43200.)/43200.)
!     else
!       fac=0.5 + 0.5*sin(tdays*pi/180. +45.)
!     endif
!     do is=1,Nsrc
!       Qbar(is)=fac*400.
!       Qbar(is)=400.
!     enddo
!
      Qbar(1)=300.
      Qbar(2)=600.
      Qbar(3)=200.
      Qbar(4)=800.
#  elif defined RIVER
!     if (tdays-dstart < 0.5) then
      if (tdays<0.5) then
        fac=1.0+TANH((time-43200.)/43200.)
      else
        fac=1.0
      endif
      do is=1,Nsrc
        if (is==2) then
          cff=13.
        else
          cff=7.
        endif
        Qbar(is)=fac*1500. * cff/27.
      enddo
#  else
#   error "### ERROR:  CPP-key 'ANA_PSOURCE' is defined, but no code"
#   error "### ERROR:  is provided to set up Qbar(is) analytically."
           xxx <-- to stop compilation here
#  endif
      do is=1,Nsrc
        do k=1,N
          Qsrc(is,k)=Qbar(is)*Qshape(is,k)
        enddo
      enddo
# endif /* SOLVE3D */
# if defined PSOURCE && defined SOLVE3D
!
!  Set-up tracer (tracer units) point Sources/Sinks.
!
#  ifdef IBERIA
      do k=1,N
!       do is=1,Nsrc
          Tsrc(1,k,itemp)=13.0
          Tsrc(1,k,isalt)=10.0
          Tsrc(2,k,itemp)=13.0
          Tsrc(2,k,isalt)=10.0
          Tsrc(3,k,itemp)=14.0
          Tsrc(3,k,isalt)=10.0
          Tsrc(4,k,itemp)=14.0
          Tsrc(4,k,isalt)=10.0
!       enddo
      enddo
#  elif defined RIVER
      do k=1,N
        do is=1,Nsrc
          Tsrc(is,k,itemp)=13.0
          Tsrc(is,k,isalt)=0.0
        enddo
      enddo
#  else
#   error "### ERROR:  CPP-key 'ANA_PSOURCE' is defined, but no code"
#   error "### ERROR:  is provided to set up Tsrc(is) analytically."
           xxx <-- to stop compilation here
#  endif
# endif
      end
#endif

      subroutine ana_empty
      end
