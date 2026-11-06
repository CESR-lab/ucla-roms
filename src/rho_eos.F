#include "cppdefs.opt"

      subroutine rho_eos(tidx)

!  ![
! Compute density anomaly from T,S via Equation Of State (EOS) for
!-------- ------- ------- ----  for seawater. Following Jackett and
! McDougall, 1995, physical EOS is assumed to have form
!
!                           rho0 + rho1(T,S)
!           rho(T,S,z) = ------------------------                 (1)
!                         1 - 0.1*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation[kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
!     K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2.  (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
!                           rho0
!                     ----------------                            (3)
!                      1 - 0.1|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.]  This results in
!
!                   rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1|z|]
!    rho1 + 0.1|z| -----------------------------------------------
!                        K00 + K01 + (K1-0.1)*|z| + K2*|z|^2
!                                                                 (4)
! which is suitable for pressure-gradient calculation.

!
! Optionally, if CPP-switch SPLIT_EOS is defined, term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
!              rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z|              (5)
!
! where
!                            rho1(T,S) - rho0*K01(T,S)/K00
!             qp1(T,S)= 0.1 -------------------------------       (6)
!                                   K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
!         d rho     d rho1           d qp1             d |z|
!        ------- = -------- + |z| * -------  +  qp1 * -------     (7)
!         d x,s      d x,s           d x,s             d x,s
!
!                  |<----- adiabatic ----->|   |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Vaisala frequency [1/s2] at
! horizontal RHO-points and vertical W-points,
!
!                          g    d rho  |
!             bvf^2 = - ------ ------- |                          (8)
!                        rho0    d z   | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! WARNING: Shared target arrays in the code below: "rho1",
!          "bvf" (if needed), and
!
!          SPLIT_EOS is defined: "qp1"  ["rho" does not exist]
!                   not defined  "rho"  ["qp1" does not exist]
!
!
! Reference:  Jackett, D. R. and T. J. McDougall, 1995, Minimal
!             Adjustment of Hydrostatic Profiles to Achieve Static
!             Stability. J. Atmos. Ocean. Tec., vol. 12, pp. 381-389.
!
! << This equation of state formulation has been derived by Jackett
! and McDougall (1992), unpublished manuscript, CSIRO, Australia. It
! computes in-situ density anomaly as a function of potential
! temperature (Celsius) relative to the surface, salinity (PSU),
! and depth (meters).  It assumes  no  pressure  variation along
! geopotential  surfaces,  that  is,  depth  and  pressure  are
! interchangeable. >>
!                                          John Wilkin, 29 July 92
!  !]

      use dimensions
      use roms_mpi
      use param
      use tracers, only: t
      use coupling
      use eos_vars
      use grid
      use mixing
      use ocean_vars
      use scalars

      implicit none

      integer, intent(in) :: tidx

      integer i,j,k
# ifdef NONLIN_EOS
      real Ts,Tt,sqrtTs, dr00, K0

! EOS coefficient notation: first index is for salinity, second for
! temperature powers, "S" for 3/2 power of salinity.

      real, parameter ::   r00=999.842594,     r01=6.793952E-2,
     &      r02=-9.095290E-3,   r03=1.001685E-4,    r04=-1.120083E-6,
     &                                              r05=6.536332E-9,
     &      r10=0.824493,       r11=-4.08990E-3,    r12=7.64380E-5,
     &                          r13=-8.24670E-7,    r14=5.38750E-9,
     &      rS0=-5.72466E-3,    rS1=1.02270E-4,     rS2=-1.65460E-6,
     &      r20=4.8314E-4,

     &      K00=19092.56,       K01=209.8925,       K02=-3.041638,
     &                          K03=-1.852732e-3,   K04=-1.361629e-5,
     &      K10=104.4077,       K11=-6.500517,      K12=0.1553190,
     &                                              K13=2.326469e-4,
     &      KS0=-5.587545,      KS1=+0.7390729,     KS2=-1.909078e-2

      real K0_Duk 

# endif /* NONLIN_EOS */
      real cff,cff1
      real dpth

# ifdef NONLIN_EOS
      Tt=3.8D0 ; Ts=34.5D0 ; sqrtTs=sqrt(Ts)
      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

      dr00=r00-rho0
# endif

      do j=1-bf,ny+bf
# ifdef NONLIN_EOS
        do k=1,nz                                      ! NONLINEAR
          do i=1-bf,nx+bf                              !  EQUATION
            Tt=t(i,j,k,tidx,itemp)                     !  OF STATE
            Ts=t(i,j,k,tidx,isalt) ; sqrtTs=sqrt(max(0.,Ts))

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                            r04+Tt*r05 ))))

     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                             r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(rS1+Tt*rS2))
     &                                                +Ts*r20 ) )
#  ifdef MASKING
     &                                              *rmask(i,j)
#  endif
            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

            qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)
     &                               /((K00+K0)*(K00+K0_Duk))
#   ifdef MASKING
     &                                            *rmask(i,j)
#   endif
#  if defined LMD_KPP || defined LMD_BKPP || defined LMD_RIMIX
          enddo
        enddo

        cff=g/rho0
        do k=1,nz-1
          do i=1-bf,nx+bf

c**         dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dpth=          -0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            bvf(i,j,k)=-cff*( rho1(i,j,k+1)-rho1(i,j,k)  ! Elementary
     &                        +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic
     &                            *dpth*(1.- qp2*dpth)   ! difference
     &                       )/(z_r(i,j,k+1)-z_r(i,j,k))
#   ifdef MASKING
     &                   *rmask(i,j)
#   endif
#  endif
          enddo
        enddo
#  if defined LMD_KPP || defined LMD_BKPP || defined LMD_RIMIX
        do i=1-bf,nx+bf
          bvf(i,j,nz)=bvf(i,j,nz-1)
          bvf(i,j,0 )=bvf(i,j,   1)
        enddo
#  endif

# else /* NONLIN_EOS */

       ! Linear equation of state

        cff=Tcoef*T0-Scoef*S0

        do k=1,nz
          do i=1-bf,nx+bf
#  ifdef SPLIT_EOS
            rho1(i,j,k)=cff -Tcoef*t(i,j,k,tidx,itemp) 
     &                      +Scoef*t(i,j,k,tidx,isalt)
            qp1(i,j,k)=0.
#  else
            rho(i,j,k)=cff  -Tcoef*t(i,j,k,tidx,itemp)
     &                      +Scoef*t(i,j,k,tidx,isalt)
#  endif

#  ifdef MASKING
#   ifdef SPLIT_EOS
            rho1(i,j,k)=rho1(i,j,k)*rmask(i,j)
#   else
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#   endif
#  endif
          enddo
        enddo
#  if defined LMD_KPP || defined LMD_BKPP || defined LMD_RIMIX
        cff=g/rho0
        do k=1,nz-1
          do i=1-bf,nx+bf
            bvf(i,j,k)=cff*(rho(i,j,k)-rho(i,j,k+1))
     &                     /(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo
        do i=1-bf,nx+bf
          bvf(i,j,nz)=bvf(i,j,nz-1)
          bvf(i,j,0 )=bvf(i,j,   1)
        enddo
#  endif
# endif /* NONLIN_EOS */

! Computation of 2-way averaged densities for barotropic mode.
! In the code segment below "rhoA" is vertically averaged density
! perturbation normalized by rho0, while "rhoS" is vertically
! integrated pressure normalized by (1/2)*total_depth^2/rho0,
! i.e., by vertically integrated pressure generated by water
! column with uniform density rho0.  Hence both "rhoA" and "rhoS"
! are nondimensional quantities of comparable value.

# ifdef VAR_RHO_2D
        do i=1-bf,nx+bf
#  ifdef SPLIT_EOS
c**       dpth=z_w(i,j,N)-z_r(i,j,N)
          dpth=-z_r(i,j,nz)

          cff=Hz(i,j,nz)*(rho1(i,j,nz)+qp1(i,j,nz)*dpth*(1.-qp2*dpth))
#  else
          cff=Hz(i,j,nz)*rho(i,j,nz)
#  endif
          rhoS(i,j)=0.5*cff*Hz(i,j,nz)
          rhoA(i,j)=cff
        enddo
        do k=nz-1,1,-1
          do i=1-bf,nx+bf
#  ifdef SPLIT_EOS
            dpth=-z_r(i,j,k)
            cff=Hz(i,j,k)*(rho1(i,j,k)+qp1(i,j,k)*dpth*(1.-qp2*dpth))
#  else
            cff=Hz(i,j,k)*rho(i,j,k)
#  endif
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5*cff)
            rhoA(i,j)=rhoA(i,j)+cff
          enddo
        enddo
        cff1=1./rho0
        do i=1-bf,nx+bf
          cff=1./(z_w(i,j,nz)-z_w(i,j,0))
          rhoA(i,j)=cff*cff1*rhoA(i,j)
          rhoS(i,j)=2.*cff*cff*cff1*rhoS(i,j)
        enddo
# endif         /* VAR_RHO_2D */
# ifdef NONTRAD_COR
        cff=0.5*rho0/g
        do k=1,nz
          do i=1-bf,nx+bf
#  if !defined SPLIT_EOS
            stop 'NONTRAD_COR needs SPLIT_EOS'
#  endif
            rho1(i,j,k)=rho1(i,j,k) +cff*(
     &            f_XI(i,j)*(v(i,j,k,tidx)+v(i,j+1,k,tidx))
     &          -f_ETA(i,j)*(u(i,j,k,tidx)+u(i+1,j,k,tidx))
     &                                                    )
          enddo
        enddo
# endif
      enddo    ! <-- j

      end
