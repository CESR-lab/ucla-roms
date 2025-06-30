#include "cppdefs.opt"
#ifdef SOLVE3D
# define DUKO_2001


      subroutine rho_eos(tidx)
      use param
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer, intent(in) :: tidx
      integer,save :: tile=0

# include "compute_tile_bounds.h"
# if defined NONLIN_EOS && !defined SPLIT_EOS
      call rho_eos_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2),tidx)
# else
      call rho_eos_tile (istr,iend,jstr,jend,tidx)
# endif
      end


# if defined NONLIN_EOS && !defined SPLIT_EOS
      subroutine rho_eos_tile (istr,iend,jstr,jend, K_up,K_dw,tidx)
# else
      subroutine rho_eos_tile (istr,iend,jstr,jend,tidx)
# endif

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

      integer istr,iend,jstr,jend, i,j,k
# ifdef NONLIN_EOS
      real Ts,Tt,sqrtTs, dr00, K0
#  ifndef SPLIT_EOS
      real K1, K2
      real dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: K_up,K_dw
#  endif

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

#  ifndef SPLIT_EOS
      real, parameter ::   B00=0.4721788,
     &      B01=0.01028859,     B02=-2.512549e-4,   B03=-5.939910e-7,
     &      B10=-0.01571896,    B11=-2.598241e-4,   B12=7.267926e-6,
     &                          BS1=2.042967e-3,

     &      E00=+1.045941e-5,   E01=-5.782165e-10,  E02=+1.296821e-7,
     &      E10=-2.595994e-7,   E11=-1.248266e-9,   E12=-3.508914e-9
#  endif
#  ifdef DUKO_2001
      real K0_Duk  !, rho1_0
#  endif
# endif /* NONLIN_EOS */
      real cff,cff1
# ifdef SPLIT_EOS
      real dpth
# endif

# include "compute_extended_bounds.h"

# ifdef NONLIN_EOS
      Tt=3.8D0 ; Ts=34.5D0 ; sqrtTs=sqrt(Ts)
#  ifdef DUKO_2001
      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

#  endif
      dr00=r00-rho0
# endif

      do j=jstrR,jendR
# ifdef NONLIN_EOS
        do k=1,N                                       ! NONLINEAR
          do i=istrR,iendR                             !  EQUATION
            Tt=t(i,j,k,tidx,itemp)                     !  OF STATE
#  ifdef SALINITY
            Ts=t(i,j,k,tidx,isalt) ; sqrtTs=sqrt(max(0.,Ts))
#  endif
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

#  ifdef SPLIT_EOS
#   ifdef DUKO_2001
            qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)
     &                               /((K00+K0)*(K00+K0_Duk))
#   else
            qp1(i,j,k)=( 0.000649D0 +0.1D0*( K00*rho1(i,j,k)
     &                             -rho0*K0)/(K00*(K00+K0)) )
#   endif
#   ifdef MASKING
     &                                            *rmask(i,j)
#   endif
#  else
            K1=B00+Tt*(B01+Tt*(B02+Tt*B03)) +Ts*( B10+Tt*( B11
     &                                    +Tt*B12 )+sqrtTs*BS1 )

            K2=E00+Tt*(E01+Tt*E02) +Ts*(E10+Tt*(E11+Tt*E12))

            dpth=-z_r(i,j,k)
            cff=K00-0.1*dpth
            cff1=K0+dpth*(K1+K2*dpth)
            rho(i,j,k)=( rho1(i,j,k)*cff*(K00+cff1)
     &                           -0.1*dpth*rho0*cff1
     &                            )/(cff*(cff+cff1))
#   ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#   endif
#  endif
#  if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP
#   ifndef SPLIT_EOS
            dpth=-z_w(i,j,k)
            K_up(i,k)=K0+dpth*(K1+K2*dpth)

            dpth=-z_w(i,j,k-1)
            K_dw(i,k)=K0+dpth*(K1+K2*dpth)
#   endif
          enddo
        enddo

        cff=g/rho0
        do k=1,N-1
          do i=istrR,iendR

#   ifdef SPLIT_EOS
c**         dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dpth=          -0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            bvf(i,j,k)=-cff*( rho1(i,j,k+1)-rho1(i,j,k)  ! Elementary
     &                        +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic
     &                            *dpth*(1.- qp2*dpth)   ! difference
     &                       )/(z_r(i,j,k+1)-z_r(i,j,k))
#   else
            cff1=0.1*(z_w(i,j,N)-z_w(i,j,k))
            bvf(i,j,k)=-cff*(   (rho1(i,j,k+1)-rho1(i,j,k))
     &                        *(K00+K_dw(i,k+1))*(K00+K_up(i,k))

     &                     -cff1*( rho0*(K_dw(i,k+1)-K_up(i,k))
     &                            +K00*(rho1(i,j,k+1)-rho1(i,j,k))
     &                                  +rho1(i,j,k+1)*K_dw(i,k+1)
     &                                      -rho1(i,j,k)*K_up(i,k)

     &          ) )/(  (K00+K_dw(i,k+1)-cff1)*(K00+K_up(i,k)-cff1)
     &                                  *(z_r(i,j,k+1)-z_r(i,j,k))
     &                                                          )

#   endif
#   ifdef MASKING
     &                   *rmask(i,j)
#   endif
#  endif
          enddo
        enddo
#  if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
   || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING
        do i=istrR,iendR
          bvf(i,j,N)=bvf(i,j,N-1)
          bvf(i,j,0)=bvf(i,j,  1)
        enddo
#  endif

# else /* NONLIN_EOS */
!-----------------> LINEAR EQUATION OF STATE
#  ifdef SALINITY
        cff=Tcoef*T0-Scoef*S0
#  else
        cff=Tcoef*T0
#  endif
        do k=1,N
          do i=istrR,iendR
#  ifdef SPLIT_EOS
            rho1(i,j,k)=cff -Tcoef*t(i,j,k,tidx,itemp)
#   ifdef SALINITY
     &                      +Scoef*t(i,j,k,tidx,isalt)
#   endif
            qp1(i,j,k)=0.
#  else
            rho(i,j,k)=cff  -Tcoef*t(i,j,k,tidx,itemp)
#   ifdef SALINITY
     &                      +Scoef*t(i,j,k,tidx,isalt)
#   endif
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
#  if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
   || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING
        cff=g/rho0
        do k=1,N-1
          do i=istrR,iendR
            bvf(i,j,k)=cff*(rho(i,j,k)-rho(i,j,k+1))
     &                     /(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo
        do i=istrR,iendR
          bvf(i,j,N)=bvf(i,j,N-1)
          bvf(i,j,0)=bvf(i,j,  1)
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
        do i=istrR,iendR
#  ifdef SPLIT_EOS
c**       dpth=z_w(i,j,N)-z_r(i,j,N)
          dpth=-z_r(i,j,N)

          cff=Hz(i,j,N)*(rho1(i,j,N)+qp1(i,j,N)*dpth*(1.-qp2*dpth))
#  else
          cff=Hz(i,j,N)*rho(i,j,N)
#  endif
          rhoS(i,j)=0.5*cff*Hz(i,j,N)
          rhoA(i,j)=cff
        enddo
        do k=N-1,1,-1
          do i=istrR,iendR
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
        do i=istrR,iendR
          cff=1./(z_w(i,j,N)-z_w(i,j,0))
          rhoA(i,j)=cff*cff1*rhoA(i,j)
          rhoS(i,j)=2.*cff*cff*cff1*rhoS(i,j)
        enddo
# endif         /* VAR_RHO_2D */
# ifdef NON_TRADITIONAL
        cff=0.5*rho0/g
        do k=1,N
          do i=istrR,iendR
            rho1(i,j,k)=rho1(i,j,k) +cff*(
     &            f_XI(i,j)*(v(i,j,k,tidx)+v(i,j+1,k,tidx))
     &          -f_ETA(i,j)*(u(i,j,k,tidx)+u(i+1,j,k,tidx))
     &                                                    )
          enddo
        enddo
# endif
      enddo    ! <-- j

      end
#else
      subroutine rho_eos_empty
      end
#endif /* SOLVE3D */
