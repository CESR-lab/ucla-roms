#include "cppdefs.h"
#ifdef LMD_KPP

      subroutine lmd_swfrac_tile (istr,iend,jstr,jend,
     &                     imin,imax,jmin,jmax, Z,swdk)
!
! Compute fraction of solar shortwave flux penetrating to the
! specified depth "Z" [meters, POSITIVE DOWNWARD] due to exponential
! decay in Jerlov water type using Paulson and Simpson (1977)
! two-wavelength band solar  absorption model.
!
! output: swdk    shortwave (radiation) fractional decay.
!
! Reference:  Paulson, C.A., and J.J. Simpson, 1977: Irradiance
! meassurements in the upper ocean, J. Phys. Oceanogr., 7, 952-956.
!
!  This routine was adapted from Bill Large 1995 code.
!
      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j, Jwt
      real Z(PRIVATE_2D_SCRATCH_ARRAY),
     &     swdk(PRIVATE_2D_SCRATCH_ARRAY),
     &     mu1(5),mu2(5), r1(5), attn1,attn2,
     &     frac1,frac2, xi1,xi2, swdk1,swdk2

      mu1(1)=0.35    !  Define reciprocal of the absorption
      mu1(2)=0.6     !  coefficient for each of two solar
      mu1(3)=1.0     !  wavelength bands as a function
      mu1(4)=1.5     !  of water type (Ref: Paulson and
      mu1(5)=1.4     !  Simpson, 1977).

      mu2(1)=23.0
      mu2(2)=20.0
      mu2(3)=17.0
      mu2(4)=14.0
      mu2(5)=7.9

      r1(1)=0.58     !  Define fraction of the total radiance
      r1(2)=0.62     !  for wavelength band 1 as a function of
      r1(3)=0.67     !  Jerlov water type. The fraction for
      r1(4)=0.77     !  wavelength band 2 is r2=1-r1.
      r1(5)=0.78
                     !  Set Jerlov water type to assign
      Jwt=1          !  everywhere; an integer from 1 to 5.

      attn1=-1./mu1(Jwt)
      attn2=-1./mu2(Jwt)
      frac1=r1(Jwt)
      frac2=(1.-r1(Jwt))

      do j=jmin,jmax
        do i=imin,imax
          xi1=attn1*Z(i,j)
          if (xi1>-20.) then      ! Quick fix to avoid
            swdk1=frac1*exp(xi1)     ! computing exp(-15000)
          else
            swdk1=0.
          endif

          xi2=attn2*Z(i,j)
          if (xi2>-20.) then
            swdk2=frac2*exp(xi2)
          else
            swdk2=0.
          endif

          swdk(i,j)=swdk1+swdk2
        enddo
      enddo
      return
      end
#else
      subroutine lmd_swfrac_empty
      end
#endif /* LMD_KPP */


