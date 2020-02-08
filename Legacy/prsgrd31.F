#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd31_tile (istr,iend,jstr,jend, A3d(1,1), A3d(1,2),
     &                                         A2d(1,1), A2d(1,2))
      end

      subroutine prsgrd31_tile (istr,iend,jstr,jend, ru,rv, rsurf,
     &                                                       pgrd)
      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: rsurf
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY) :: pgrd

! Compute pressure gradient term: STANDARD(POM) or WEIGHTED JACOBIAN
!
! input:  rho               ! density anomaly [kg/m^3].
! output: ru = - pgrd_XI    ! XI- and ETA-components of pressure
!         rv = - pgrd_ETA   !     -gradient term as it appears in
!                           ! the r.h.s. of 3D momentum equations.
!
! Reference: Song, Y.T., 1998: A general pressure gradient
! formutlation for numerical ocean models. Part I: Scheme design
! and diagnostic analysis. Mon. Wea. Rev., v. 126, pp. 3213-3230.
!
! Switches: WJ_GRADP:  WEIGHTED/STANDARD jacobian formulation.
c**
c**   Original coefficient by Tony was 0.25
c**
c--#define WJ_GRADP 0.125

! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----
! Computation starts with extrapolation of density field toward sea
! surface , after which compute pressure gradient at the topmost grid
! box around u(:,:,N) point, including the contribution due to free-
! surface elevation (barotropic part) and due to the density difference
! in the top-most grid box (baroclinic part). This operation
! initializes vertical integration. Once done, proceed to the grid
! points below throughout the vertical column using either Weighted or
! Standard Jacobian.  The standard jacobian is rewritten in diagonal
! form (lile in Lin, 1997), which minimizes the number of operations
! relatively to any other form.

      real cff
# ifdef WJ_GRADP
      real gamma
# endif
# include "grid.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "scalars.h"

# include "compute_auxiliary_bounds.h"

      do j=jstrV-1,jend
        do i=istrU-1,iend
          rsurf(i,j)=rho(i,j,N) + (rho(i,j,N)-rho(i,j,N-1))
     &                              *(z_w(i,j,N)-z_r(i,j,N))
     &                            /(z_r(i,j,N)-z_r(i,j,N-1))
        enddo

        if (j >= jstr) then
          cff=0.5*g/rho0
          do i=istrU,iend
            pgrd(i)=(g+cff*(rsurf(i-1,j)+rsurf(i,j)))*( z_w(i-1,j,N)
     &                                                   -z_w(i,j,N))

     &     +cff*( (rho(i-1,j,N)-rsurf(i,j))*(z_w(i-1,j,N)-z_r(i,j,N))
     &           +(rsurf(i-1,j)-rho(i,j,N))*(z_w(i,j,N)-z_r(i-1,j,N))
     &                                                              )

            ru(i,j,N)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))*dn_u(i,j)*pgrd(i)
          enddo
# ifdef WJ_GRADP
          cff=0.25*g/rho0
# else
          cff=0.5*g/rho0
# endif
          do k=N-1,1,-1
            do i=istrU,iend
# ifdef WJ_GRADP
              gamma=WJ_GRADP
     &         *(z_r(i,j,k+1)-z_r(i-1,j,k+1)+z_r(i,j,k)-z_r(i-1,j,k))
     &         *(z_r(i,j,k+1)-z_r(i-1,j,k+1)-z_r(i,j,k)+z_r(i-1,j,k))
     &     /((z_r(i,j,k+1)-z_r(i,j,k))*(z_r(i-1,j,k+1)-z_r(i-1,j,k)))
# endif
            pgrd(i)=pgrd(i)-cff*(
# ifdef WJ_GRADP
     &                    ( (1.+gamma)*(rho(i,j,k+1)-rho(i-1,j,k+1))
     &                     +(1.-gamma)*(rho(i,j,k  )-rho(i-1,j,k  )))
     &                               *( z_r(i,j,k+1)+z_r(i-1,j,k+1)
     &                                 -z_r(i,j,k  )-z_r(i-1,j,k  ))

     &                               -( rho(i,j,k+1)+rho(i-1,j,k+1)
     &                                 -rho(i,j,k  )-rho(i-1,j,k  ))
     &                   *( (1.+gamma)*(z_r(i,j,k+1)-z_r(i-1,j,k+1))
     &                     +(1.-gamma)*(z_r(i,j,k  )-z_r(i-1,j,k  )))
# else
     &                                   (rho(i,j,k+1)-rho(i-1,j,k))
     &                                  *(z_r(i-1,j,k+1)-z_r(i,j,k))
     &                                  +(rho(i,j,k)-rho(i-1,j,k+1))
     &                                  *(z_r(i,j,k+1)-z_r(i-1,j,k))
# endif
     &                                                             )

              ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*pgrd(i)
            enddo
          enddo     !--> discard pgrd
        endif

! ETA-component: same sequence as above.
!---- ---------- ---- -------- -- ------

        if (j >= jstrV) then
          cff=0.5*g/rho0
          do i=istr,iend
            pgrd(i)=(g+cff*(rsurf(i,j-1)+rsurf(i,j)))*( z_w(i,j-1,N)
     &                                                   -z_w(i,j,N))

     &     +cff*( (rho(i,j-1,N)-rsurf(i,j))*(z_w(i,j-1,N)-z_r(i,j,N))
     &           +(rsurf(i,j-1)-rho(i,j,N))*(z_w(i,j,N)-z_r(i,j-1,N))
     &                                                              )

            rv(i,j,N)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))*dm_v(i,j)*pgrd(i)
          enddo
# ifdef WJ_GRADP
          cff=0.25*g/rho0
# else
          cff=0.5*g/rho0
# endif
          do k=N-1,1,-1
            do i=istr,iend
# ifdef WJ_GRADP
              gamma=WJ_GRADP
     &         *(z_r(i,j,k+1)-z_r(i,j-1,k+1)+z_r(i,j,k)-z_r(i,j-1,k))
     &         *(z_r(i,j,k+1)-z_r(i,j-1,k+1)-z_r(i,j,k)+z_r(i,j-1,k))
     &     /((z_r(i,j,k+1)-z_r(i,j,k))*(z_r(i,j-1,k+1)-z_r(i,j-1,k)))
# endif
              pgrd(i)=pgrd(i)-cff*(
# ifdef WJ_GRADP
     &                    ( (1.+gamma)*(rho(i,j,k+1)-rho(i,j-1,k+1))
     &                     +(1.-gamma)*(rho(i,j,k  )-rho(i,j-1,k  )))
     &                               *( z_r(i,j,k+1)+z_r(i,j-1,k+1)
     &                                 -z_r(i,j,k  )-z_r(i,j-1,k  ))

     &                               -( rho(i,j,k+1)+rho(i,j-1,k+1)
     &                                 -rho(i,j,k  )-rho(i,j-1,k  ))
     &                   *( (1.+gamma)*(z_r(i,j,k+1)-z_r(i,j-1,k+1))
     &                     +(1.-gamma)*(z_r(i,j,k  )-z_r(i,j-1,k  )))
# else
     &                                   (rho(i,j,k+1)-rho(i,j-1,k))
     &                                  *(z_r(i,j-1,k+1)-z_r(i,j,k))
     &                                  +(rho(i,j,k)-rho(i,j-1,k+1))
     &                                  *(z_r(i,j,k+1)-z_r(i,j-1,k))
# endif
     &                                                             )

              rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*pgrd(i)
            enddo
          enddo    !--> discard pgrd
        endif
      enddo
      end
#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */
