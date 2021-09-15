#include "cppdefs.h"
#if defined ANA_VMIX && defined SOLVE3D

      subroutine ana_vmix (tile)

      use param

      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call ana_vmix_tile (istr,iend,jstr,jend)
      end

      subroutine ana_vmix_tile (istr,iend,jstr,jend)

! Set vertical kinematic vertical viscosity "Akv" [m^2/s] for momentum
! and vertical mixing c efficients "Akt" [m^2/s] for tracer equations
! using analytical expressions.

      use param
      use grid

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real cff

# include "ocean3d.h"
# include "mixing.h"
# include "scalars.h"

# include "compute_auxiliary_bounds.h"

# if defined DOUBLE_GYRE
      if (FIRST_TIME_STEP) then
        cff=1./100.
        do k=0,N
          do j=jstr,jend
            do i=istr,iend
             Akv(i,j,k)=0.001+0.25*exp(cff*(z_w(i,j,k)-z_w(i,j,N)))

             Akt(i,j,k,itemp)=0.
#  ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
            enddo
          enddo
        enddo
      endif

# elif defined CANYON_B
      cff=1./50.                    ! Setup both top and bottom
      do k=0,N                      ! boundary layers.
        do j=jstr,jend
          do i=istr,iend
            Akv(i,j,k)=0.001+0.0095*(
     &                      exp(cff*(z_w(i,j,k)-z_w(i,j,N)))
     &                     +exp(cff*(z_w(i,j,0)-z_w(i,j,k)))
     &                                                     )
            Akt(i,j,k,itemp)=Akt_bak(itemp)
          enddo
        enddo
      enddo
# elif defined UPWELLING || defined IBERIA
      cff=1./100.
      do k=0,N
        do j=jstr,jend
          do i=istr,iend
           Akv(i,j,k)=0.005+0.015*exp(cff*(z_w(i,j,k)-z_w(i,j,N)))



c           Akv(i,j,k)=0.004+0.008*exp(cff*(z_w(i,j,k)-z_w(i,j,N)))
            Akt(i,j,k,itemp)=0.2*Akt_bak(itemp)
c??
c??
c             Akv(i,j,k)=0.
c             Akt(i,j,k,itemp)=0.
c??
c??
#  ifdef SALINITY
            Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
          enddo
        enddo
      enddo
# else
      do k=1,N-1
        do j=jstr,jend
          do i=istr,iend
            Akv(i,j,k)=???
            Akt(i,j,k,itemp)=???
            Akt(i,j,k,isalt)=???
          enddo
        enddo
      enddo
# endif /* CANYON_B */

# ifdef EXCHANGE
      call exchange_2_tile (istr,iend,jstr,jend,  Akv, N+1,
     &                    Akt(START_2D_ARRAY,0,itemp), N+1)
#  ifdef SALINITY
      call exchange_tile (istr,iend,jstr,jend,
     &                    Akt(START_2D_ARRAY,0,isalt), N+1)
#  endif
# endif
      end
#else
      subroutine ana_vmix_empty
      end
#endif /* ANA_VMIX */
