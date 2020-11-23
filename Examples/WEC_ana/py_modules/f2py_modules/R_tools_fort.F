# 1 "R_tools_fort.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "R_tools_fort.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROMS ROUTINES
!!
!! copied from actual ROMS scripts
!!
!! compile with:
!! "cpp R_tools_fort.F R_tools_fort.f"
!! "f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m R_tools_fort R_tools_fort.f" for python use
!!
!! print R_tools_fort.rho_eos.__doc__
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



# 1 "set_global_definitions.h" 1
# 15 "set_global_definitions.h"
c--#define ALLOW_SINGLE_BLOCK_MODE
# 263 "set_global_definitions.h"
c-#ifdef
c-# define float dfloat
c-# define FLoaT dfloat
c-# define FLOAT dfloat
c-# define sqrt dsqrt
c-# define SQRT dsqrt
c-# define exp dexp
c-# define EXP dexp
c-# define dtanh dtanh
c-# define TANH dtanh
c-#endif
# 19 "R_tools_fort.F" 2
# 1 "cppdefs.h" 1
# 30 "cppdefs.h"
c--#ifndef MAX_GRID_LEVEL
c--# define MAX_GRID_LEVEL 2
c--# include "redefs.X"
c--#endif
# 43 "cppdefs.h"
c# define PSOURCE
c# define Q_PSOURCE
c# define PASSIVE_TRACER
c# define ANA_PSOURCE
# 56 "cppdefs.h"
c# define EXACT_RESTART







c------# define DIURNAL_SRFLUX
# 77 "cppdefs.h"
c# define LMD_RIMIX
# 773 "cppdefs.h"
# 1 "set_global_definitions.h" 1
# 15 "set_global_definitions.h"
c--#define ALLOW_SINGLE_BLOCK_MODE
# 263 "set_global_definitions.h"
c-#ifdef
c-# define float dfloat
c-# define FLoaT dfloat
c-# define FLOAT dfloat
c-# define sqrt dsqrt
c-# define SQRT dsqrt
c-# define exp dexp
c-# define EXP dexp
c-# define dtanh dtanh
c-# define TANH dtanh
c-#endif
# 774 "cppdefs.h" 2
# 20 "R_tools_fort.F" 2
# 1 "rho_eos.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute density anomaly (adapted from rho_eos.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine rho_eos(Lm,Mm,N, T,S, z_r,z_w,rho0, rho)



!
! Compute density anomaly from T,S via Equation Of State (EOS) for
!-------- ------- ------- ---- for seawater. Following Jackett and
! McDougall, 1995, physical EOS is assumed to have form
!
! rho0 + rho1(T,S)
! rho(T,S,z) = ------------------------ (1)
! 1 - 0.1*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation[kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
! K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2. (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
! rho0
! ---------------- (3)
! 1 - 0.1|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.] This results in
!
! rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1|z|]
! rho1 + 0.1|z| -----------------------------------------------
! K00 + K01 + (K1-0.1)*|z| + K2*|z|^2
! (4)
! which is suitable for pressure-gradient calculation.
!
! Optionally, if CPP-switch is defined, term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
! rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z| (5)
!
! where
! rho1(T,S) - rho0*K01(T,S)/K00
! qp1(T,S)= 0.1 ------------------------------- (6)
! K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
! d rho d rho1 d qp1 d |z|
! ------- = -------- + |z| * ------- + qp1 * ------- (7)
! d x,s d x,s d x,s d x,s
!
! |<----- adiabatic ----->| |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Vaisala frequency [1/s^2] at
! horizontal RHO-points and vertical W-points,
!
! g d rho |
! bvf^2 = - ------ ------- | (8)
! rho0 d z | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! WARNING: Shared target arrays in the code below: "rho1",
! "bvf" (if needed), and
!
! is defined: "qp1" ["rho" does not exist]
! not defined "rho" ["qp1" does not exist]
!
!
! Reference: Jackett, D. R. and T. J. McDougall, 1995, Minimal
! Adjustment of Hydrostatic Profiles to Achieve Static
! Stability. J. Atmos. Ocean. Tec., vol. 12, pp. 381-389.
!
! << This equation of state formulation has been derived by Jackett
! and McDougall (1992), unpublished manuscript, CSIRO, Australia. It
! computes in-situ density anomaly as a function of potential
! temperature (Celsius) relative to the surface, salinity (PSU),
! and depth (meters). It assumes no pressure variation along
! geopotential surfaces, that is, depth and pressure are
! interchangeable. >>
! John Wilkin, 29 July 92
!


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 T(0:Lm+1,0:Mm+1,N), S(0:Lm+1,0:Mm+1,N),
     & z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N),
     & rho(0:Lm+1,0:Mm+1,N),
     & rho1(0:Lm+1,0:Mm+1,N), qp1(0:Lm+1,0:Mm+1,N),
     & dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00

      real*8, parameter :: r00=999.842594, r01=6.793952E-2,
     & r02=-9.095290E-3, r03=1.001685E-4, r04=-1.120083E-6,
     & r05=6.536332E-9,
     & r10=0.824493, r11=-4.08990E-3, r12=7.64380E-5,
     & r13=-8.24670E-7, r14=5.38750E-9,
     & rS0=-5.72466E-3, rS1=1.02270E-4, rS2=-1.65460E-6,
     & r20=4.8314E-4,
     & K00=19092.56, K01=209.8925, K02=-3.041638,
     & K03=-1.852732e-3, K04=-1.361629e-5,
     & K10=104.4077, K11=-6.500517, K12=0.1553190,
     & K13=2.326469e-4,
     & KS0=-5.587545, KS1=+0.7390729, KS2=-1.909078e-2,
     & qp2=0.0000172, g=9.81





      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N, T,S, z_r,z_w, rho0
Cf2py intent(out) rho



!!!!!!!!!!!!!!!!!!
# 160 "rho_eos.F"
      dr00=r00-rho0

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunk_size -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )




      do j=jstr,jend
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)

            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     & r04+Tt*r05 ))))
     & +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     & r13+Tt*r14 )))
     & +sqrtTs*(rS0+Tt*(
     & rS1+Tt*rS2 ))+Ts*r20 ))

            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     & +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     & +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
# 206 "rho_eos.F"
            qp1(i,j,k)=( 0.000649D0 +0.1D0*( K00*rho1(i,j,k)
     & -rho0*K0)/(K00*(K00+K0)) )


         rho(i,j,k) = rho1(i,j,k) + qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k))


          enddo
        enddo



      enddo ! <-- j

      return
      end
# 21 "R_tools_fort.F" 2
# 1 "bvf_eos.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute density anomaly (adapted from rho_eos.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bvf_eos(Lm,Mm,N, T,S, z_r,z_w,rho0,bvf)



!
! Compute density anomaly from T,S via Equation Of State (EOS) for
!-------- ------- ------- ---- for seawater. Following Jackett and
! McDougall, 1995, physical EOS is assumed to have form
!
! rho0 + rho1(T,S)
! rho(T,S,z) = ------------------------ (1)
! 1 - 0.1*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation[kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
! K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2. (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
! rho0
! ---------------- (3)
! 1 - 0.1|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.] This results in
!
! rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1|z|]
! rho1 + 0.1|z| -----------------------------------------------
! K00 + K01 + (K1-0.1)*|z| + K2*|z|^2
! (4)
! which is suitable for pressure-gradient calculation.
!
! Optionally, if CPP-switch is defined, term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
! rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z| (5)
!
! where
! rho1(T,S) - rho0*K01(T,S)/K00
! qp1(T,S)= 0.1 ------------------------------- (6)
! K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
! d rho d rho1 d qp1 d |z|
! ------- = -------- + |z| * ------- + qp1 * ------- (7)
! d x,s d x,s d x,s d x,s
!
! |<----- adiabatic ----->| |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Vaisala frequency [1/s^2] at
! horizontal RHO-points and vertical W-points,
!
! g d rho |
! bvf^2 = - ------ ------- | (8)
! rho0 d z | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! WARNING: Shared target arrays in the code below: "rho1",
! "bvf" (if needed), and
!
! is defined: "qp1" ["rho" does not exist]
! not defined "rho" ["qp1" does not exist]
!
!
! Reference: Jackett, D. R. and T. J. McDougall, 1995, Minimal
! Adjustment of Hydrostatic Profiles to Achieve Static
! Stability. J. Atmos. Ocean. Tec., vol. 12, pp. 381-389.
!
! << This equation of state formulation has been derived by Jackett
! and McDougall (1992), unpublished manuscript, CSIRO, Australia. It
! computes in-situ density anomaly as a function of potential
! temperature (Celsius) relative to the surface, salinity (PSU),
! and depth (meters). It assumes no pressure variation along
! geopotential surfaces, that is, depth and pressure are
! interchangeable. >>
! John Wilkin, 29 July 92
!


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 T(0:Lm+1,0:Mm+1,N), S(0:Lm+1,0:Mm+1,N),
     & z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N),
     & rho(0:Lm+1,0:Mm+1,N), bvf(0:Lm+1,0:Mm+1,0:N),
     & rho1(0:Lm+1,0:Mm+1,N), qp1(0:Lm+1,0:Mm+1,N),
     & dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00

      real*8, parameter :: r00=999.842594, r01=6.793952E-2,
     & r02=-9.095290E-3, r03=1.001685E-4, r04=-1.120083E-6,
     & r05=6.536332E-9,
     & r10=0.824493, r11=-4.08990E-3, r12=7.64380E-5,
     & r13=-8.24670E-7, r14=5.38750E-9,
     & rS0=-5.72466E-3, rS1=1.02270E-4, rS2=-1.65460E-6,
     & r20=4.8314E-4,
     & K00=19092.56, K01=209.8925, K02=-3.041638,
     & K03=-1.852732e-3, K04=-1.361629e-5,
     & K10=104.4077, K11=-6.500517, K12=0.1553190,
     & K13=2.326469e-4,
     & KS0=-5.587545, KS1=+0.7390729, KS2=-1.909078e-2,
     & qp2=0.0000172, g=9.81





      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N, T,S, z_r,z_w, rho0
Cf2py intent(out) bvf



!!!!!!!!!!!!!!!!!!
# 160 "bvf_eos.F"
      dr00=r00-rho0

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunk_size -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )




      do j=jstr,jend
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)

            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     & r04+Tt*r05 ))))
     & +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     & r13+Tt*r14 )))
     & +sqrtTs*(rS0+Tt*(
     & rS1+Tt*rS2 ))+Ts*r20 ))

            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     & +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     & +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))






               qp1(i,j,k) = (0.000649D0 + 0.1D0 * (K00*rho1(i,j,k)
     & - rho0*K0)/(K00*(K00+K0)) )



          enddo
        enddo

        cff=g/rho0
        do k=1,N-1
          do i=imin,imax

            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            cff2=( rho1(i,j,k+1)-rho1(i,j,k) ! Elementary
     & +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic
     & *dpth*(1.-2.*qp2*dpth) ! difference
     & )


            bvf(i,j,k)=-cff*cff2 / (z_r(i,j,k+1)-z_r(i,j,k))

          enddo
        enddo




        do i=imin,imax
          bvf (i,j,N) = bvf (i,j,N-1)
          bvf (i,j,0) = bvf (i,j, 1)
        enddo




      enddo ! <-- j

      return
      end
# 22 "R_tools_fort.F" 2
# 1 "prsgrd.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute pressure gradient (adapted from prsgrd32AC1.F in ROMS code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine prsgrd(Lm,Mm,N, z_r, z_w, rho, rho0, rho1, qp1, Hz,
     & dn_u, dm_v, ru ,rv)
! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensianal coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid. These
! derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.
!
! In the code below, if CPP-switch is defined, the Equation
! of State (EOS) is assumed to have form
!
! rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1.-qp2*dpth]
!
! where rho1 is potential density at 1 atm and qp1 is compressibility
! coefficient, which does not depend on z, and dpth=zeta-z, and qp2
! is just a constant. In this case
!
! d rho d rho1 d qp1 d z
! ------- = ------ + ----- *dpth*[..] - qp1*[1.-2.*qp2*dpth]*------
! d s,x d s,x d s,x d s,x
!
! |<--- adiabatic part --->| |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach quarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! statification in physical sense as long as discrete values of
! density are positively stratified.
!
! This scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho)
! [with the exception of harmonic averaging algorithm in the case
! when CPP-switch is defined, see above]. If parameter
! OneFifth (see above) is set to zero, the scheme becomes identical
! to standard Jacobian.
!
! NOTE: This routine is an alternative form of prsgrd32 and it
! produces results identical to that if its prototype.


      implicit none
      integer Lm, Mm, N, imin, imax, jmin, jmax, i,j,k,
     & numthreads, trd, chunksize, margin, jstr,
     & jend, istr, iend, istrU, jstrV


      real*8 rho(0:Lm+1,0:Mm+1,N), z_r(0:Lm+1,0:Mm+1,N),
     & z_w(0:Lm+1,0:Mm+1,0:N), rho1(0:Lm+1,0:Mm+1,N),
     & qp1(0:Lm+1,0:Mm+1,N), Hz(0:Lm+1,0:Mm+1,N),
     & dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1), dpth
      real*8, dimension(0:Lm+1,0:Mm+1,N) :: ru,rv,P
      real*8, dimension(0:Lm+1,0:Mm+1) :: dR,dZ
      real*8, dimension(0:Lm+1,0:Mm+1) :: FC,dZx,rx,dRx
      real*8 GRho, HalfGRho, cff, cfr, rho0
      real*8, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.
     & ,g=9.81, qp2=0.0000172


!!!! INPUTS: Lm,Mm,N --> not inputted in python
!!!! !!!!!!!!!!!!! rho , rho1,qp1,qp2---> all returned from rho_eos_V2.F


Cf2py intent(in) Lm,Mm,N,z_r,z_w,rho,rho0, rho1,qp1,Hz, dn_u, dm_v
Cf2py intent(out) ru, rv


      ! Preliminary step (same for XI and ETA-components
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      numthreads=1
      trd=0
      chunksize=(jmax-jmin + numthreads)/numthreads
      margin=(chunksize*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunksize -margin, jmin )
      jend=min( (trd+1)*chunksize-1-margin, jmax )

      trd=0
      chunksize=(imax-imin + numthreads)/numthreads
      margin=(chunksize*numthreads -imax+imin-1)/2
      istr=max( trd *chunksize -margin, imin )
      iend=min( (trd+1)*chunksize-1-margin, imax )

      !!! not using this code for entire domain


      istrU = istr+1
      jstrV = jstr+1

      GRho=g/rho0
      HalfGRho=0.5 * GRho

      do j=jstr-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1) - z_r(i,j,k)

c** dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dpth= -0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k) ! elementary
     & +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ! difference



          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1 !---> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0
            endif

c** dpth=z_w(i,j,N)-z_r(i,j,k)
            dpth= -z_r(i,j,k)
            dR(i,k)=dR(i,k) -qp1(i,j,k)*dZ(i,k)*(1.-2.*qp2*dpth)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)

          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     & +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     & /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     & *(z_r(i,j,k+1)-z_r(i,j,k))

     & -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     & -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     & -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     & -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     & ))
          enddo
        enddo
      enddo !<-- j


      ! Compute XI-component of pressure gradient term:
      !------------------------------------------------

      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))

c** dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
c** & -z_r(i,j,k)-z_r(i-1,j,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i-1,j,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k) ! elementary
     & +(qp1(i,j,k)-qp1(i-1,j,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ) ! difference



          enddo
        enddo
        do j=jstr,jend
           do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif

            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     & *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))

            enddo

          do i=istrU,iend
            ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     & P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     & (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))

     & -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     & -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     & )))
          enddo
        enddo
!
!! ETA-component of pressure gradient term
!-----------------------------------------
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))

c** dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
c** & -z_r(i,j,k)-z_r(i,j-1,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i,j-1,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k) ! elementary
     & +(qp1(i,j,k)-qp1(i,j-1,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ) ! difference



          enddo
        enddo
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif

            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     & *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))

          enddo
          if (j.ge.jstrV) then
            do i=istr,iend
              rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     & P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     & (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))

     & -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     & -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     & )))
            enddo
          endif
        enddo
      enddo !<-- k
      return
      end
# 23 "R_tools_fort.F" 2
# 1 "rho_eos_v2.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute density anomaly (adapted from rho_eos.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c#define DUKO_2001

      subroutine rho_eos_v2(Lm,Mm,N, T,S, z_r,z_w,rho0, rho,rho1,qp1)




!
! Compute density anomaly from T,S via Equation Of State (EOS) for
!-------- ------- ------- ---- for seawater. Following Jackett and
! McDougall, 1995, physical EOS is assumed to have form
!
! rho0 + rho1(T,S)
! rho(T,S,z) = ------------------------ (1)
! 1 - 0.1*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation[kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
! K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2. (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
! rho0
! ---------------- (3)
! 1 - 0.1|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.] This results in
!
! rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1|z|]
! rho1 + 0.1|z| -----------------------------------------------
! K00 + K01 + (K1-0.1)*|z| + K2*|z|^2
! (4)
! which is suitable for pressure-gradient calculation.
!
! Optionally, if CPP-switch is defined, term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
! rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z| (5)
!
! where
! rho1(T,S) - rho0*K01(T,S)/K00
! qp1(T,S)= 0.1 ------------------------------- (6)
! K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
! d rho d rho1 d qp1 d |z|
! ------- = -------- + |z| * ------- + qp1 * ------- (7)
! d x,s d x,s d x,s d x,s
!
! |<----- adiabatic ----->| |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Vaisala frequency [1/s^2] at
! horizontal RHO-points and vertical W-points,
!
! g d rho |
! bvf^2 = - ------ ------- | (8)
! rho0 d z | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! WARNING: Shared target arrays in the code below: "rho1",
! "bvf" (if needed), and
!
! is defined: "qp1" ["rho" does not exist]
! not defined "rho" ["qp1" does not exist]
!
!
! Reference: Jackett, D. R. and T. J. McDougall, 1995, Minimal
! Adjustment of Hydrostatic Profiles to Achieve Static
! Stability. J. Atmos. Ocean. Tec., vol. 12, pp. 381-389.
!
! << This equation of state formulation has been derived by Jackett
! and McDougall (1992), unpublished manuscript, CSIRO, Australia. It
! computes in-situ density anomaly as a function of potential
! temperature (Celsius) relative to the surface, salinity (PSU),
! and depth (meters). It assumes no pressure variation along
! geopotential surfaces, that is, depth and pressure are
! interchangeable. >>
! John Wilkin, 29 July 92
!


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 T(0:Lm+1,0:Mm+1,N), S(0:Lm+1,0:Mm+1,N),
     & z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N),
     & rho(0:Lm+1,0:Mm+1,N),
     & rho1(0:Lm+1,0:Mm+1,N), qp1(0:Lm+1,0:Mm+1,N),
     & dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00

      real*8, parameter :: r00=999.842594, r01=6.793952E-2,
     & r02=-9.095290E-3, r03=1.001685E-4, r04=-1.120083E-6,
     & r05=6.536332E-9,
     & r10=0.824493, r11=-4.08990E-3, r12=7.64380E-5,
     & r13=-8.24670E-7, r14=5.38750E-9,
     & rS0=-5.72466E-3, rS1=1.02270E-4, rS2=-1.65460E-6,
     & r20=4.8314E-4,
     & K00=19092.56, K01=209.8925, K02=-3.041638,
     & K03=-1.852732e-3, K04=-1.361629e-5,
     & K10=104.4077, K11=-6.500517, K12=0.1553190,
     & K13=2.326469e-4,
     & KS0=-5.587545, KS1=+0.7390729, KS2=-1.909078e-2,
     & qp2=0.0000172, g=9.81





      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N, T,S, z_r,z_w, rho0
Cf2py intent(out) rho, rho1, qp1




!!!rho_eos_V2 is different than original
!!! in that it returns rho1,qp1,qp2 for prsgrd calculation


!!!!!!!!!!!!!!!!!!
# 167 "rho_eos_v2.F"
      dr00=r00-rho0

! qp2_O=qp2

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunk_size -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )




      do j=jstr,jend
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)

            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     & r04+Tt*r05 ))))
     & +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     & r13+Tt*r14 )))
     & +sqrtTs*(rS0+Tt*(
     & rS1+Tt*rS2 ))+Ts*r20 ))

            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     & +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     & +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
# 215 "rho_eos_v2.F"
            qp1(i,j,k)=( 0.000649D0 +0.1D0*( K00*rho1(i,j,k)
     & -rho0*K0)/(K00*(K00+K0)) )


         rho(i,j,k) = rho1(i,j,k) + qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k))


          enddo
        enddo



      enddo ! <-- j

      return
      end
# 24 "R_tools_fort.F" 2
# 1 "prsgrd_nosplit.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute pressure gradient (adapted from prsgrd32AC1.F in ROMS code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine prsgrd_nosplit(Lm,Mm,N, z_r, z_w, rho, rho0, Hz,
     & dn_u, dm_v, ru ,rv)
! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensianal coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid. These
! derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.
!
! In the code below, if CPP-switch is defined, the Equation
! of State (EOS) is assumed to have form
!
! rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1.-qp2*dpth]
!
! where rho1 is potential density at 1 atm and qp1 is compressibility
! coefficient, which does not depend on z, and dpth=zeta-z, and qp2
! is just a constant. In this case
!
! d rho d rho1 d qp1 d z
! ------- = ------ + ----- *dpth*[..] - qp1*[1.-2.*qp2*dpth]*------
! d s,x d s,x d s,x d s,x
!
! |<--- adiabatic part --->| |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach quarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! statification in physical sense as long as discrete values of
! density are positively stratified.
!
! This scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho)
! [with the exception of harmonic averaging algorithm in the case
! when CPP-switch is defined, see above]. If parameter
! OneFifth (see above) is set to zero, the scheme becomes identical
! to standard Jacobian.
!
! NOTE: This routine is an alternative form of prsgrd32 and it
! produces results identical to that if its prototype.


      implicit none
      integer Lm, Mm, N, imin, imax, jmin, jmax, i,j,k,
     & numthreads, trd, chunksize, margin, jstr,
     & jend, istr, iend, istrU, jstrV


      real*8 rho(0:Lm+1,0:Mm+1,N), z_r(0:Lm+1,0:Mm+1,N),
     & z_w(0:Lm+1,0:Mm+1,0:N),
     & Hz(0:Lm+1,0:Mm+1,N),
     & dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1), dpth
      real*8, dimension(0:Lm+1,0:Mm+1,N) :: ru,rv,P
      real*8, dimension(0:Lm+1,0:Mm+1) :: dR,dZ
      real*8, dimension(0:Lm+1,0:Mm+1) :: FC,dZx,rx,dRx
      real*8 GRho, HalfGRho, cff, cfr, rho0
      real*8, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.
     & ,g=9.81, qp2=0.0000172


!!!! INPUTS: Lm,Mm,N --> not inputted in python
!!!! !!!!!!!!!!!!! rho , rho1,qp1,qp2---> all returned from rho_eos_V2.F


Cf2py intent(in) Lm,Mm,N,z_r,z_w,rho,rho0,Hz, dn_u, dm_v
Cf2py intent(out) ru, rv


      ! Preliminary step (same for XI and ETA-components
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      numthreads=1
      trd=0
      chunksize=(jmax-jmin + numthreads)/numthreads
      margin=(chunksize*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunksize -margin, jmin )
      jend=min( (trd+1)*chunksize-1-margin, jmax )

      trd=0
      chunksize=(imax-imin + numthreads)/numthreads
      margin=(chunksize*numthreads -imax+imin-1)/2
      istr=max( trd *chunksize -margin, imin )
      iend=min( (trd+1)*chunksize-1-margin, imax )

      !!! not using this code for entire domain
      istrU = istr+1
      jstrV = jstr+1
      GRho=g/rho0
      HalfGRho=0.5 * GRho

      do j=jstr-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1) - z_r(i,j,k)
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)

          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1 !---> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0
            endif
         enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     & +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     & /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     & *(z_r(i,j,k+1)-z_r(i,j,k))

     & -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     & -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     & -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     & -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     & ))
          enddo
        enddo
      enddo !<-- j


      ! Compute XI-component of pressure gradient term:
      !------------------------------------------------

      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
            rx(i,j)=(rho(i,j,k)-rho(i-1,j,k))
         enddo
        enddo
        do j=jstr,jend
           do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif
           enddo

          do i=istrU,iend
            ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     & P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     & (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))

     & -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     & -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     & )))
          enddo
        enddo
!
!! ETA-component of pressure gradient term
!-----------------------------------------
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
            rx(i,j)=(rho(i,j,k)-rho(i,j-1,k))

          enddo
        enddo
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif
         enddo
          if (j.ge.jstrV) then
            do i=istr,iend
              rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     & P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     & (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))

     & -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     & -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     & )))
            enddo
          endif
        enddo
      enddo !<-- k
      return
      end
# 25 "R_tools_fort.F" 2
# 1 "prsgrd32ac1_momentum.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute pressure gradient (adapted from prsgrd32AC1.F in ROMS code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine prsgrd32ac1_momentum(Lm,Mm,N, z_r, z_w, rho, rho0, rho1
     & , qp1, pm,pn,ru ,rv)
! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensianal coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid. These
! derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.
!
! In the code below, if CPP-switch is defined, the Equation
! of State (EOS) is assumed to have form
!
! rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1.-qp2*dpth]
!
! where rho1 is potential density at 1 atm and qp1 is compressibility
! coefficient, which does not depend on z, and dpth=zeta-z, and qp2
! is just a constant. In this case
!
! d rho d rho1 d qp1 d z
! ------- = ------ + ----- *dpth*[..] - qp1*[1.-2.*qp2*dpth]*------
! d s,x d s,x d s,x d s,x
!
! |<--- adiabatic part --->| |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach quarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! statification in physical sense as long as discrete values of
! density are positively stratified.
!
! This scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho)
! [with the exception of harmonic averaging algorithm in the case
! when CPP-switch is defined, see above]. If parameter
! OneFifth (see above) is set to zero, the scheme becomes identical
! to standard Jacobian.
!
! NOTE: This routine is an alternative form of prsgrd32 and it
! produces results identical to that if its prototype.


      implicit none
      integer Lm, Mm, N, imin, imax, jmin, jmax, i,j,k,
     & numthreads, trd, chunksize, margin, jstr,
     & jend, istr, iend, istrU, jstrV


      real*8 rho(0:Lm+1,0:Mm+1,N), z_r(0:Lm+1,0:Mm+1,N),
     & z_w(0:Lm+1,0:Mm+1,0:N), rho1(0:Lm+1,0:Mm+1,N),
     & qp1(0:Lm+1,0:Mm+1,N),
     & pn(0:Lm+1,0:Mm+1), pm(0:Lm+1,0:Mm+1), dpth
      real*8, dimension(0:Lm+1,0:Mm+1,N) :: ru,rv,P
      real*8, dimension(0:Lm+1,0:Mm+1) :: dR,dZ
      real*8, dimension(0:Lm+1,0:Mm+1) :: FC,dZx,rx,dRx
      real*8 GRho, HalfGRho, cff, cfr, rho0
      real*8, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.
     & ,g=9.81, qp2=0.0000172


!!!! INPUTS: Lm,Mm,N --> not inputted in python
!!!! !!!!!!!!!!!!! rho , rho1,qp1,qp2---> all returned from rho_eos_V2.F


Cf2py intent(in) Lm,Mm,N,z_r,z_w,rho,rho0, rho1,qp1,pm,pn
Cf2py intent(out) ru, rv


      ! Preliminary step (same for XI and ETA-components
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      numthreads=1
      trd=0
      chunksize=(jmax-jmin + numthreads)/numthreads
      margin=(chunksize*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunksize -margin, jmin )
      jend=min( (trd+1)*chunksize-1-margin, jmax )

      trd=0
      chunksize=(imax-imin + numthreads)/numthreads
      margin=(chunksize*numthreads -imax+imin-1)/2
      istr=max( trd *chunksize -margin, imin )
      iend=min( (trd+1)*chunksize-1-margin, imax )

      !!! not using this code for entire domain


      istrU = istr+1
      jstrV = jstr+1

      GRho=g/rho0
      HalfGRho=0.5 * GRho

      do j=jstr-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1) - z_r(i,j,k)

c** dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dpth= -0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k) ! elementary
     & +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ! difference



          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1 !---> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0
            endif

c** dpth=z_w(i,j,N)-z_r(i,j,k)
            dpth= -z_r(i,j,k)
            dR(i,k)=dR(i,k) -qp1(i,j,k)*dZ(i,k)*(1.-2.*qp2*dpth)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)

          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     & +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     & /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     & *(z_r(i,j,k+1)-z_r(i,j,k))

     & -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     & -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     & -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     & -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     & ))
          enddo
        enddo
      enddo !<-- j


      ! Compute XI-component of pressure gradient term:
      !------------------------------------------------

      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))

c** dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
c** & -z_r(i,j,k)-z_r(i-1,j,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i-1,j,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k) ! elementary
     & +(qp1(i,j,k)-qp1(i-1,j,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ) ! difference



          enddo
        enddo
        do j=jstr,jend
           do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif

            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     & *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))

            enddo

          do i=istrU,iend
            ru(i,j,k)=0.5*(pm(i,j)+pm(i-1,j))*(
     & P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     & (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))

     & -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     & -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     & )))
          enddo
        enddo
!
!! ETA-component of pressure gradient term
!-----------------------------------------
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))

c** dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
c** & -z_r(i,j,k)-z_r(i,j-1,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i,j-1,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k) ! elementary
     & +(qp1(i,j,k)-qp1(i,j-1,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ) ! difference



          enddo
        enddo
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif

            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     & *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))

          enddo
          if (j.ge.jstrV) then
            do i=istr,iend
              rv(i,j,k)=0.5*(pn(i,j)+pn(i,j-1))*(
     & P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     & (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))

     & -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     & -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     & )))
            enddo
          endif
        enddo
      enddo !<-- k
      return
      end
# 26 "R_tools_fort.F" 2
# 1 "prsgrd32ac1_momentum_v2.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute pressure gradient (adapted from prsgrd32AC1.F in ROMS code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine prsgrd32ac1_momentum_v2(Lm,Mm,N, z_r, z_w, rho, rho0,
     & rho1 , qp1, pm,pn,ru ,rv)
! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensianal coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid. These
! derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.
!
! In the code below, if CPP-switch is defined, the Equation
! of State (EOS) is assumed to have form
!
! rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1.-qp2*dpth]
!
! where rho1 is potential density at 1 atm and qp1 is compressibility
! coefficient, which does not depend on z, and dpth=zeta-z, and qp2
! is just a constant. In this case
!
! d rho d rho1 d qp1 d z
! ------- = ------ + ----- *dpth*[..] - qp1*[1.-2.*qp2*dpth]*------
! d s,x d s,x d s,x d s,x
!
! |<--- adiabatic part --->| |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach quarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! statification in physical sense as long as discrete values of
! density are positively stratified.
!
! This scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho)
! [with the exception of harmonic averaging algorithm in the case
! when CPP-switch is defined, see above]. If parameter
! OneFifth (see above) is set to zero, the scheme becomes identical
! to standard Jacobian.
!
! NOTE: This routine is an alternative form of prsgrd32 and it
! produces results identical to that if its prototype.


      implicit none
      integer Lm, Mm, N, imin, imax, jmin, jmax, i,j,k,
     & numthreads, trd, chunksize, margin, jstr,
     & jend, istr, iend, istrU, jstrV


      real*8 rho(0:Lm+1,0:Mm+1,N), z_r(0:Lm+1,0:Mm+1,N),
     & z_w(0:Lm+1,0:Mm+1,0:N), rho1(0:Lm+1,0:Mm+1,N),
     & qp1(0:Lm+1,0:Mm+1,N),
     & pn(0:Lm+1,0:Mm+1), pm(0:Lm+1,0:Mm+1), dpth
      real*8, dimension(0:Lm+1,0:Mm+1,N) :: ru,rv,P
      real*8, dimension(0:Lm+1,0:Mm+1) :: dR,dZ
      real*8, dimension(0:Lm+1,0:Mm+1) :: FC,dZx,rx,dRx
      real*8 GRho, HalfGRho, cff, cfr, rho0
      real*8, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.
     & ,g=9.81, qp2=0.0000172


!!!! INPUTS: Lm,Mm,N --> not inputted in python
!!!! !!!!!!!!!!!!! rho , rho1,qp1,qp2---> all returned from rho_eos_V2.F


Cf2py intent(in) Lm,Mm,N,z_r,z_w,rho,rho0, rho1,qp1,pm,pn
Cf2py intent(out) ru, rv


      ! Preliminary step (same for XI and ETA-components
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      numthreads=1
      trd=0
      chunksize=(jmax-jmin + numthreads)/numthreads
      margin=(chunksize*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunksize -margin, jmin )
      jend=min( (trd+1)*chunksize-1-margin, jmax )

      trd=0
      chunksize=(imax-imin + numthreads)/numthreads
      margin=(chunksize*numthreads -imax+imin-1)/2
      istr=max( trd *chunksize -margin, imin )
      iend=min( (trd+1)*chunksize-1-margin, imax )

      !!! not using this code for entire domain


      istrU = istr+1
      jstrV = jstr+1

      GRho=g/rho0
      HalfGRho=0.5 * GRho

      do j=jstrV-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1) - z_r(i,j,k)

c** dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dpth= -0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k) ! elementary
     & +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ! difference



          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1 !---> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0
            endif

c** dpth=z_w(i,j,N)-z_r(i,j,k)
            dpth= -z_r(i,j,k)
            dR(i,k)=dR(i,k) -qp1(i,j,k)*dZ(i,k)*(1.-2.*qp2*dpth)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)

          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     & +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     & /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     & *(z_r(i,j,k+1)-z_r(i,j,k))

     & -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     & -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     & -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     & -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     & ))
          enddo
        enddo
      enddo !<-- j


      ! Compute XI-component of pressure gradient term:
      !------------------------------------------------

      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))

c** dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
c** & -z_r(i,j,k)-z_r(i-1,j,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i-1,j,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k) ! elementary
     & +(qp1(i,j,k)-qp1(i-1,j,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ) ! difference



          enddo
        enddo
        do j=jstr,jend
           do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif

            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     & *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))

            enddo

          do i=istrU,iend
            ru(i,j,k)=0.5*(pn(i,j)+pn(i-1,j))*(
     & P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     & (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))

     & -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     & -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     & )))
          enddo
        enddo
!
!! ETA-component of pressure gradient term
!-----------------------------------------
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))

c** dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
c** & -z_r(i,j,k)-z_r(i,j-1,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i,j-1,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k) ! elementary
     & +(qp1(i,j,k)-qp1(i,j-1,k)) ! adiabatic
     & *dpth*(1.-qp2*dpth) ) ! difference



          enddo
        enddo
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif

            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     & *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))

          enddo
          if (j.ge.jstrV) then
            do i=istr,iend
              rv(i,j,k)=0.5*(pm(i,j)+pm(i,j-1))*(
     & P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     & (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))

     & -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     & -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     & -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     & -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     & )))
            enddo
          endif
        enddo
      enddo !<-- k
      return
      end
# 27 "R_tools_fort.F" 2
# 1 "sigma_to_z_intr.F" 1
      subroutine sigma_to_z_intr (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     & z_lev, var_zlv, imin,jmin,kmin, FillValue)
c subroutine sigma_to_z_intr (Lm,Mm,N, nz, z_r, z_w, rmask, var,
c & z_lev, var_zlv,FillValue)


! Interpolate field "var" defined in sigma-space to a set of several
! fixed z-levels z_lev(1:nz).





      implicit none
      integer Lm,Mm,N, nz, imin,imax,jmin,jmax, kmin, i,j,k,m

      real*8 z_r(0:Lm+1,0:Mm+1,N), rmask(0:Lm+1,0:Mm+1),
     & z_w(0:Lm+1,0:Mm+1,0:N), zlev
      real*4 z_lev(nz), var(imin:Lm+1,jmin:Mm+1,kmin:N), FillValue,
     & var_zlv(imin:Lm+1,jmin:Mm+1,nz)

      integer km(0:Lm+1)
      real*8 zz(0:Lm+1,0:N+1), dpth

     & , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


!############ F2PY INPUT/OUT FOR PYTHON !##############
Cf2py intent(in) Lm,Mm,N,nz,z_r,z_w,rmask,var,z_lev,imin,jmin,kmin,FillValue
Cf2py intent(out) var_zlv
!#######################################################

      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max( trd *chunk_size -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )






      do j=jstr,jend
        if (kmin==1) then
          if (imin==0 .and. jmin==0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,j,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,j,0)
              zz(i,N+1)=z_w(i,j,N)
            enddo
          elseif (imin==1 .and. jmin==0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            enddo
          elseif (imin==0 .and. jmin==1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
              zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            enddo
          elseif (imin==1 .and. jmin==1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     & +z_r(i,j-1,k)+z_r(i-1,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     & +z_w(i,j-1,0)+z_w(i-1,j-1,0))

              zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     & +z_w(i,j-1,N)+z_w(i-1,j-1,N))
             enddo
          endif
        else
          if (imin==0 .and. jmin==0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,j,k)
              enddo
            enddo
          elseif (imin==1 .and. jmin==0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
              enddo
            enddo
          elseif (imin==0 .and. jmin==1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
              enddo
            enddo
          elseif (imin==1 .and. jmin==1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     & +z_w(i,j-1,k)+z_w(i-1,j-1,k))
              enddo
            enddo
          endif
        endif

        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1 !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff>0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     & *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo

        do m=1,nz
          zlev=z_lev(m)
c** write(*,*) 'm =',m, ' zlev =', zlev
          if (kmin==0) then !
            do i=imin,imax !
              dpth=zz(i,N)-zz(i,0)
              if (rmask(i,j)<0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(zlev-zz(i,N))>0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(i,0)-zlev)>0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)
              if (rmask(i,j)<0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(zlev-zz(i,N+1))>0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(zlev-zz(i,N))>0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-zlev)>0.) then
                km(i)=-2 !<-- below bottom
              elseif (dpth*(zz(i,1)-zlev)>0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif

          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i)==-1) then
                if((zz(i,k+1)-zlev)*(zlev-zz(i,k)) >= 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i)==-3) then
              var_zlv(i,j,m)=0. !<-- masked out
            elseif (km(i)==-2) then
              var_zlv(i,j,m)=FillValue !<-- below bottom
            elseif (km(i)==N+2) then
              var_zlv(i,j,m)=-FillValue !<-- above surface
            elseif (km(i)==N) then
              var_zlv(i,j,m)=var(i,j,N) !-> R-point, above z_r(N)

     & +FC(i,N)*(zlev-zz(i,N))




            elseif (km(i)==kmin-1) then !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin) ! but above bottom

     & -FC(i,kmin)*(zz(i,kmin)-zlev)




            else
              k=km(i)

              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=zlev-zz(i,k)
              q=zz(i,k+1)-zlev

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     & -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     & +p*FC(i,k+1) -q*FC(i,k) )
     & )
# 244 "sigma_to_z_intr.F"
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j
      return
      end
# 28 "R_tools_fort.F" 2
