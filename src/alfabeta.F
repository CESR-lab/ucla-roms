#include "cppdefs.opt"
#if defined SOLVE3D && defined LMD_KPP

      subroutine alfabeta_tile(istr,iend,jstr,jend,
     &                         imin,imax,jmin,jmax,
     &                                  alpha,beta)

! Compute thermal expansion "alpha" and saline contraction "beta"
! coefficients as functions of potential temperature and salinity
! using polynomial expressions of Jackett & McDougall (1992).
!
!       alpha(Ts,Tt,0)=-d(rho1(Ts,Tt,0))/d(Tt) / rho0
!       beta(Ts,Tt,0) = d(rho1(Ts,Tt,0))/d(Ts) / rho0
!
! Both are evaluated at the surface.  Adapted from original "rati"
! and "beta" routines.

      use param
      use tracers, only: t
      use eos_vars
      use grid
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: alpha,beta

# ifdef NONLIN_EOS
      real Tt, cff
      real, parameter ::      r01=6.793952E-2,  r02=-9.095290E-3,
     &     r03=+1.001685E-4,  r04=-1.120083E-6, r05=+6.536332E-9,
#  ifdef SALINITY
     &       r10=+0.824493,   r11=-4.08990E-3,  r12=+7.64380E-5,
     &                        r13=-8.24670E-7,  r14=+5.38750E-9,
     &       rS0=-5.72466E-3, rS1=+1.02270E-4,  rS2=-1.65460E-6,
     &       r20=+4.8314E-4
      real Ts, sqrtTs
#  endif
# else
!#  include "eos_vars"
# endif
# ifdef NONLIN_EOS
      cff=1./rho0
# endif
      do j=jmin,jmax
        do i=imin,imax
# ifdef NONLIN_EOS
          Tt=t(i,j,N,nstp,itemp)
#  ifdef SALINITY
          Ts=t(i,j,N,nstp,isalt) ; sqrtTs=sqrt(max(0.,Ts))
#  endif
          alpha(i,j)=-cff*( r01+Tt*( 2.*r02+Tt*( 3.*r03+Tt*(
     &                                  4.*r04 +Tt*5.*r05 )))
#  ifdef SALINITY
     &                         +Ts*( r11+Tt*( 2.*r12+Tt*(
     &                                   3.*r13 +Tt*4.*r14 ))
     &                              +sqrtTs*(rS1+Tt*2.*rS2) )
#  endif
     &                                                     )
#  ifdef SALINITY
          beta(i,j)= cff*( r10+Tt*(r11+Tt*(r12+Tt*(r13+Tt*r14)))
     &                         +1.5*(rS0+Tt*(rS1+Tt*rS2))*sqrtTs
     &                                              +2.*r20*Ts )
#  endif
# else
          alpha(i,j)=abs(Tcoef)  !--> for linear Equation of State
#  ifdef SALINITY
          beta(i,j)=abs(Scoef)
#  else
          beta(i,j)=0.
#  endif
# endif /* NONLIN_EOS */
        enddo
      enddo
      end
#else
      subroutine alfabeta_empty
      end
#endif /* LMD_KPP */
