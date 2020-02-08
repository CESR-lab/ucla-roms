#include "cppdefs.h"
#define PP_MIXING
#ifdef PP_MIXING

****
****  WARNING: NOT PARALLELIZED, DO NOT USE
****

      subroutine pp_vmix

! Compute vertical mixing coefficients for momentum and tracers using
! Pacanowski and Philander  (1981)  mixing scheme wich is based on
! Richardson number.
!
!  Output:
!
!            vertical viscosity coefficient (m^2/s).
!            vertical diffusion coefficient for potential
!              temperature (m^2/s).
!            vertical diffusion coefficient for salinity (m^2/s).
!     bvf      Brunt-Vaisala frequency (1/sec).
!
!  Reference:
!
!   Pacanowski, R.C. and G.H. Philander, 1981: Parameterization of
!       vertical mixing in numerical models of tropical oceans, J.
!       Phys. Oceanog., 11, 1443-1451.

      implicit none
# include "param.h"
# include "mixing.h"

# include "scalars.h"

      integer pp_n
      real pp_alpha, pp_nu0, pp_numax
      parameter (pp_n=2, pp_alpha=5., pp_nu0=0.01, pp_numax=0.01)


      integer i, j, k
      real fac

   THESE TWO ARE TO BE RETAINED!!!

      real Ri_big,Ri_low ! Upper and lower limits for
      save Ri_big,Ri_low ! gradient Richardson number.

!  On first pass, set upper and lower gradient Richardson number
!  limiting values.

      if (iic==ntstart) then
        if (pp_n==2) then
          Ri_big=1./pp_alpha*( sqrt((pp_nu0+Akv_bak)/pp_numax)-1.)
          Ri_low=1./pp_alpha*(-sqrt((pp_nu0+Akv_bak)/pp_numax)-1.)
        else
          Ri_big=1./pp_alpha*( ((pp_nu0+Akv_bak)/pp_numax)
     &                                            **(-1./pp_n)-1.)
          Ri_low=-1.0e-6
        endif
      endif

!  Compute gradient Richardson number, Brunt-Vaisala frequency, and
!  horizontal velocity shear squared.
!-------------------------------------------------------------------

      call ri_number (Rig,bvf,shear2)

!  Compute vertical mixing coefficients for momentum and tracers.

      do k=1,N-1
        do j=jstr,jend
          do i=istr,iend
            if (Rig(i,j,k)<Ri_big .and. Rig(i,j,k)>Ri_low) then
              fac=1.+pp_alpha*Ri_big
              Akv(i,j,k)=pp_numax
            else
              fac=1.+pp_alpha*Rig(i,j,k)
              Akv(i,j,k)=Akv_bak+pp_nu0/(fac**pp_n)
            endif
            Akt(i,j,k,itemp)=Akt_bak(itemp)+Akv(i,j,k)/fac
# ifdef SALINITY
            Akt(i,j,k,isalt)=Akt_bak(isalt)+Akv(i,j,k)/fac
# endif
          enddo
        enddo
      enddo

# ifdef EXCHANGE
#  ifdef SALINITY
      call exchange_3_tile ( istr,iend,jstr,jend,
     &                       Akt(START_2D_ARRAY,0,isalt), N+1,
#  else
      call exchange_2_tile ( istr,iend,jstr,jend,
#  endif
     &            Akv, N+1,  Akt(START_2D_ARRAY,0,itemp), N+1)
# endif
      end
#else
      subroutine pp_vmix_empty
      end
#endif /* PP_MIXING */
