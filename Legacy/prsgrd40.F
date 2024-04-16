#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd40_tile (istr,iend,jstr,jend,   A3d(1,1), A3d(1,2),
     &                                 A3d(1,3), A3d(1,4), A2d(1,1))
      end

      subroutine prsgrd40_tile (istr,iend,jstr,jend, ru,rv, P,FX,FC)

! Finite-volime pressure gradient scheme by S. J. Lin, 1997
!------- ------ -------- -------- ------ -- -- -- ---- -----
! Reference: Lin, Shian-Jiann, 1997: A finite volume integration
! method for omputing pressure gradient force in general vertical
! coordinates. Quarterly Journal of Royal Meteorological Society.,
! vol. 123, pp. 1749-1762.

      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv, FX
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: P
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: FC
      real grho0, Half_g
# include "grid.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "scalars.h"

# include "compute_auxiliary_bounds.h"

      Half_g=0.5*g
      grho0=g/rho0
      do j=jstrV-1,jend
        do i=istrU-1,iend  !--> set pressure at free surface,
          P(i,j,N)=0.      !    initialize vertical integration
          FC(i,N)=0.
        enddo
        do k=N,1,-1
          do i=istrU-1,iend
            P(i,j,k-1)=P(i,j,k)+Hz(i,j,k)*rho(i,j,k)
            FX(i,j,k)=0.5*Hz(i,j,k)*( P(i,j,k)+P(i,j,k-1))
          enddo
        enddo

        if (j>=jstr) then
          do k=N,1,-1
            do i=istrU,iend
              FC(i,k-1)=0.5*(P(i,j,k-1)+P(i-1,j,k-1))
     &                        *(z_w(i,j,k-1)-z_w(i-1,j,k-1))

              ru(i,j,k)=( Half_g*(Hz(i-1,j,k)+Hz(i,j,k))
     &                             *(z_w(i-1,j,N)-z_w(i,j,N))

     &       +grho0*(FX(i-1,j,k)-FX(i,j,k)+FC(i,k)-FC(i,k-1))
     &                                            )*dn_u(i,j)
            enddo
          enddo
        endif

        if (j>=jstrV) then
          do k=N,1,-1
            do i=istr,iend
              FC(i,k-1)=0.5*(P(i,j,k-1)+P(i,j-1,k-1))
     &                        *(z_w(i,j,k-1)-z_w(i,j-1,k-1))

              rv(i,j,k)=( Half_g*(Hz(i,j-1,k)+Hz(i,j,k))
     &                             *(z_w(i,j-1,N)-z_w(i,j,N))

     &       +grho0*(FX(i,j-1,k)-FX(i,j,k)+FC(i,k)-FC(i,k-1))
     &                                            )*dm_v(i,j)
            enddo
          enddo
        endif
      enddo    !<-- j
      end

#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */
