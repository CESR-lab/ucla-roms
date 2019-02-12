#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd41_tile (istr,iend,jstr,jend,   A3d(1,1), A3d(1,2),
     &                                 A3d(1,3), A2d(1,1), A2d(1,2))
      end

      subroutine prsgrd41_tile (istr,iend,jstr,jend, ru,rv, FX,FC,P)

      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv, FX
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: FC,P

      real cff,cff1
      real, parameter ::  Half=0.5, Two=2., OneSixth=1./6.

# include "grid.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "scalars.h"

! Modified version of S. J. Lin, 1997 finite-volime pressure gradient
! ======== ======= == == == ==== ==== ============= ======== ========
! scheme. Original method of Lin uses trapezoidal integration rules
! for pressure. Unfortunately, piecewise-linear assumption about the
! structure of pressure field corresponds to only piecewise-constant
! density, and, as the result, original Lin's method fails to produce
! exact error cancellation if density is linear function of depth.
! The following code is a modified version, which is also cast into
! a finite volume framework (see ru,rv which appear in form of
! differences of two pairs of fluxes), but these fluxes are computed
! by assuming  bi-linear density distribution within the control
! volume with subsequent analytical integration.
!
! Reference: Lin, Shian-Jiann, 1997: A finite volume integration
!            method for omputing pressure gradient force in general
!            vertical coordinates. Quarterly Journal of Royal
!            Meteorological Society., vol. 123, pp. 1749-1762.
!
! NOTE: This method assumes that discretized density rho(i,j,k) is a
! set of grid box averages, rather than a set of values at locations
! z_r(i,j,k). Consequently, z_r(i,j,k) does not appear anywhere in
! this code.

# include "compute_auxiliary_bounds.h"

      do j=jstrV-1,jend         ! Vertical derivatives of density.
        do k=1,N-1
          do i=istrU-1,iend
             FC(i,k)=Two*(rho(i,j,k+1)-rho(i,j,k))
     &                     /(Hz(i,j,k+1)+Hz(i,j,k))
          enddo
        enddo                   ! FX are pressure fluxes across
                                ! side facets of controm volumes.
        do i=istrU-1,iend
          FX(i,j,N)=Half*Hz(i,j,N)*Hz(i,j,N)*( rho(i,j,N)
     &                      +OneSixth*Hz(i,j,N)*FC(i,N-1))
          P(i,N-1)=rho(i,j,N)*Hz(i,j,N)
        enddo
        do k=N-1,2,-1
          do i=istrU-1,iend
            FX(i,j,k)=Hz(i,j,k)*( P(i,k)+Half*Hz(i,j,k)*( rho(i,j,k)
     &        +OneSixth*Hz(i,j,k)*( FC(i,k  )*(Hz(i,j,k+1)+Hz(i,j,k))
     &                             +FC(i,k-1)*(Hz(i,j,k)+Hz(i,j,k-1))
     &                   )/(Hz(i,j,k+1)+Two*Hz(i,j,k)+Hz(i,j,k-1)) ))

            P(i,k-1)=P(i,k)+rho(i,j,k)*Hz(i,j,k)
          enddo
        enddo
        do i=istrU-1,iend
          FX(i,j,1)=Hz(i,j,1)*( P(i,1)+Half*Hz(i,j,1)*( rho(i,j,1)
     &                               +OneSixth*Hz(i,j,1)*FC(i,1) ))
        enddo  !--> discard P,FC

        if (j>=jstr) then     ! In the code segment below P(i,k-1)
          do i=istrU,iend       ! is pressure averaged over the lower
            P(i,N)=0.           ! horizontal (along sigma) facet of
            FC(i,N)=0.          ! the control volume around velocity
          enddo                 ! u(i,j,k), while FC is net pressure
          cff=0.5*g             ! force accross that facet. Lengthy
          cff1=g/rho0           ! expression after OneSixth is net
          do k=N,1,-1           ! water mass in the control volume.
            do i=istrU,iend
              P(i,k-1)=P(i,k) + OneSixth*( rho(i-1,j,k)*Hz(i-1,j,k)
     &            +(rho(i,j,k)+rho(i-1,j,k))*(Hz(i,j,k)+Hz(i-1,j,k))
     &                                        +rho(i,j,k)*Hz(i,j,k))

              FC(i,k-1)=P(i,k-1)*(z_w(i,j,k-1)-z_w(i-1,j,k-1))

              ru(i,j,k)=( cff*(Hz(i-1,j,k)+Hz(i,j,k))*(z_w(i-1,j,N)
     &                                                  -z_w(i,j,N))

     &               +cff1*(FX(i-1,j,k)-FX(i,j,k)+FC(i,k)-FC(i,k-1))
     &                                                   )*dn_u(i,j)
            enddo
          enddo   !--> discard P,FC
        endif

        if (j>=jstrV) then
          do i=istr,iend
            P(i,N)=0.
            FC(i,N)=0.
          enddo
          cff=0.5*g
          cff1=g/rho0
          do k=N,1,-1
            do i=istr,iend
              P(i,k-1)=P(i,k) + OneSixth*( rho(i,j-1,k)*Hz(i,j-1,k)
     &            +(rho(i,j,k)+rho(i,j-1,k))*(Hz(i,j,k)+Hz(i,j-1,k))
     &                                        +rho(i,j,k)*Hz(i,j,k))

              FC(i,k-1)=P(i,k-1)*(z_w(i,j,k-1)-z_w(i,j-1,k-1))

              rv(i,j,k)=( cff*(Hz(i,j-1,k)+Hz(i,j,k))*(z_w(i,j-1,N)
     &                                                  -z_w(i,j,N))

     &               +cff1*(FX(i,j-1,k)-FX(i,j,k)+FC(i,k)-FC(i,k-1))
     &                                                   )*dm_v(i,j)
            enddo
          enddo    !--> discard P,FC
        endif
      enddo    !<-- j
      end

#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */
