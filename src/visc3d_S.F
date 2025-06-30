#include "cppdefs.opt"
# if defined SOLVE3D && defined UV_VIS2

      subroutine visc3d

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none
      integer,parameter :: tile=0

# include "compute_tile_bounds.h"
      call visc3d_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                                       A2d(1,3), A2d(1,4))
      end

      subroutine visc3d_tile (istr,iend,jstr,jend, UFx,UFe,VFx,VFe)

! Compute harmonic mixing of momentum, along constant S-surfaces as
! horizontal divergence of the stress tensor. A transverse isotropy is
! assumed so the stress tensor is splitted into vertical and horizontal
! subtensors.
!
! Reference: Wajsowicz, R. C., 1993: A consistent formulation of
!            the anisotropic stress tensor for use in models of
!            the large-scale ocean circulation, JCP, 105, 333-338.

      use dimensions
      use param
      use diagnostics
      use coupling
      use grid
      use mixing
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j,k, indx
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe
      real cff

# include "compute_auxiliary_bounds.h"

      indx=3-nstp    !--> time index for target arrays;

!  Compute horizontal harmonic viscosity along constant S-surfaces.
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.

      do k=1,N
        do j=jstrV-1,jend
          do i=istrU-1,iend
            cff=0.5*Hz(i,j,k)*visc2_r(i,j)*(

     &      dn_r(i,j)*pm(i,j)*( (pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nstp)
     &                         -(pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nstp)
     &                                                              )

     &     -dm_r(i,j)*pn(i,j)*( (pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nstp)
     &                         -(pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nstp)
     &                                                             ))
            UFx(i,j)=  cff*dn_r(i,j)*dn_r(i,j)
            VFe(i,j)= -cff*dm_r(i,j)*dm_r(i,j)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend+1
            cff=0.125*(Hz(i-1,j,k)+Hz(i,j,k)+Hz(i-1,j-1,k)
     &                         +Hz(i,j-1,k))*visc2_p(i,j)*(

     &      0.25*(pm(i-1,j)+pm(i,j)+pm(i-1,j-1)+pm(i,j-1))*dn_p(i,j)
     &                     *( (pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nstp)
     &                       -(pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nstp)
     &                                                              )

     &     +0.25*(pn(i-1,j)+pn(i,j)+pn(i-1,j-1)+pn(i,j-1))*dm_p(i,j)
     &                     *( (pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nstp)
     &                       -(pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nstp)
     &                                                             ))
# ifdef MASKING
     &                                                   *pmask(i,j)
# endif
            UFe(i,j)= cff*dm_p(i,j)*dm_p(i,j)
            VFx(i,j)= cff*dn_p(i,j)*dn_p(i,j)
          enddo
        enddo

! Apply viscous terms. Note that at this stage arrays u,v(...,3-nstp)
! contain Hz*U and Hz*V with units of [m2/s].   Also compute vertical
! integral of viscous terms and add it into coupling terms for the
! barotropic mode

        do j=jstr,jend
          do i=istrU,iend
            cff=0.125*(pm(i-1,j)+pm(i,j))*(pn(i-1,j) +pn(i,j))
     &             *( (pn(i-1,j)+pn(i,j))*(UFx(i,j)-UFx(i-1,j))
     &               +(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i,j))
     &                                                        )
            rufrc(i,j)=rufrc(i,j) + cff
            u(i,j,k,indx)=u(i,j,k,indx) + dt*cff
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            cff=0.125*(pm(i,j)+pm(i,j-1))*(pn(i,j) +pn(i,j-1))
     &             *( (pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j))
     &               +(pm(i,j-1)+pm(i,j))*(VFe(i,j)-VFe(i,j-1))
     &                                                        )
            rvfrc(i,j)=rvfrc(i,j) + cff
            v(i,j,k,indx)=v(i,j,k,indx) + dt*cff
          enddo
        enddo
# if defined NHMG
!   Do things
# endif
# ifdef DIAGNOSTICS
      if (diag_uv .and. calc_diag) then
        do j=1,ny
          do i=1,nx
            Udiag(i,j,k,ihmix) = 0.5*dxdyi_u(i,j)*
     &              ( (pn(i-1,j)+pn(i,j))*(UFx(i,j)-UFx(i-1,j))
     &               +(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i,j)) )
            Vdiag(i,j,k,ihmix) = 0.5*dxdyi_v(i,j)*
     &              ( (pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j))
     &               +(pm(i,j-1)+pm(i,j))*(VFe(i,j)-VFe(i,j-1)) )
          enddo
        enddo
      endif
# endif
      enddo ! k loop

      end

#else
      subroutine visc3d_empty
      end
#endif
