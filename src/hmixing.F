      module hmixing
!     Contains subroutines that add horizontal mixing to ru/rv

#include "cppdefs.opt"

      use dimensions
      use diagnostics
      use ocean_vars
      use grid
      use mixing
      use scalars

      implicit none

      private

      real,dimension(:,:),allocatable :: UFx,UFy,VFx,VFy

      public uv_visc3d

      contains
!------------------------------------------------------------------
      subroutine init_hmixing
      ! allocate work space for the stress components

      allocate( UFx(GLOBAL_2D_ARRAY) ); UFx=0.
      allocate( UFy(GLOBAL_2D_ARRAY) ); UFy=0.
      allocate( VFx(GLOBAL_2D_ARRAY) ); VFx=0.
      allocate( VFy(GLOBAL_2D_ARRAY) ); VFy=0.

      end subroutine init_hmixing
!------------------------------------------------------------------
      subroutine uv_visc3d(ru,rv,k)

! Compute harmonic mixing of momentum, along constant S-surfaces as
! horizontal divergence of the stress tensor. A transverse isotropy is
! assumed so the stress tensor is splitted into vertical and horizontal
! subtensors.
!
! Reference: Wajsowicz, R. C., 1993: A consistent formulation of
!            the anisotropic stress tensor for use in models of
!            the large-scale ocean circulation, JCP, 105, 333-338.
!
! Stresses evaluated at n+1/2 (nrhs). Called from step3d_uv1 only
! The flux divergence is accumulated into ru,rv in area-integrated
! units [m4/s2]; the barotropic mode receives the vertical integral
! through the rufrc = sum_k ru assembly in step3d_uv1.

      ! Import/export
      integer, intent(in) :: k
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv
      ! Local
      integer i,j
      real cff

      if (.not.allocated(UFx)) call init_hmixing

!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.

      do j=0,ny
        do i=0,nx
          cff=0.5*Hz(i,j,k)*visc2_r(i,j)*(

     &    dn_r(i,j)*pm(i,j)*( (pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)
     &                       -(pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs)
     &                                                            )

     &   -dm_r(i,j)*pn(i,j)*( (pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)
     &                       -(pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)
     &                                                           ))
          UFx(i,j)=  cff*dn_r(i,j)*dn_r(i,j)
          VFy(i,j)= -cff*dm_r(i,j)*dm_r(i,j)
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx+1
          cff=0.125*(Hz(i-1,j,k)+Hz(i,j,k)+Hz(i-1,j-1,k)
     &                       +Hz(i,j-1,k))*visc2_p(i,j)*(

     &    0.25*(pm(i-1,j)+pm(i,j)+pm(i-1,j-1)+pm(i,j-1))*dn_p(i,j)
     &                   *( (pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)
     &                     -(pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs)
     &                                                            )

     &   +0.25*(pn(i-1,j)+pn(i,j)+pn(i-1,j-1)+pn(i,j-1))*dm_p(i,j)
     &                   *( (pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)
     &                     -(pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)
     &                                                           ))
     &                                                 *pmask(i,j)
          UFy(i,j)= cff*dm_p(i,j)*dm_p(i,j)
          VFx(i,j)= cff*dn_p(i,j)*dn_p(i,j)
        enddo
      enddo

! Add the flux divergence to ru,rv

      do j=1,ny
        do i=1,nx
          ru(i,j,k) = ru(i,j,k)
     &      + 0.5*( (pn(i-1,j)+pn(i,j))*(UFx(i,j)  -UFx(i-1,j))
     &             +(pm(i-1,j)+pm(i,j))*(UFy(i,j+1)-UFy(i,j)) )

          rv(i,j,k) = rv(i,j,k)
     &      + 0.5*( (pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j))
     &             +(pm(i,j-1)+pm(i,j))*(VFy(i,j)  -VFy(i,j-1)) )
        enddo
      enddo

      if (diag_uv) then
        do j=1,ny
          do i=1,nx
            Udiag(i,j,k,ihmix) = 0.5*dxdyi_u(i,j)*
     &              ( (pn(i-1,j)+pn(i,j))*(UFx(i,j)  -UFx(i-1,j))
     &               +(pm(i-1,j)+pm(i,j))*(UFy(i,j+1)-UFy(i,j)) )
            Vdiag(i,j,k,ihmix) = 0.5*dxdyi_v(i,j)*
     &              ( (pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j))
     &               +(pm(i,j-1)+pm(i,j))*(VFy(i,j)  -VFy(i,j-1)) )
          enddo
        enddo
      endif

      end subroutine uv_visc3d

      end module hmixing
