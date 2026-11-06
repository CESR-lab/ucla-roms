#include "cppdefs.opt"

      subroutine t3dmix
!     Computes horizontal mixing of tracers along geopotential surfaces.
      use param
      use tracers,  only: t, nt
      use scalars,  only: nnew
      use roms_mpi, only: exchange_xxx
      use private_scratch

      implicit none

      integer :: itrc
      logical,parameter :: t3d_dif2_gp=.true.
      logical,parameter :: t3d_dif4_gp=.false.

      if (.not.t3d_dif2_gp.and..not.t3d_dif4_gp) then
        ! nothing to do
        return
      endif

      if (t3d_dif2_gp.and.t3d_dif4_gp) then
        if (mynode.eq.0)
     &    print *,'only one of t3d_dif2_gp or t3d_dif4_gp can be true'
        stop
      endif

      if (t3d_dif2_gp) then
        do itrc=1,nt
          call t3dmix_dif2_gp (itrc,A2d(1,1),A2d(1,2),
     &                              A2d(1,3),A2d(1,5),
     &                              A2d(1,7),A2d(1,9),
     &                              A2d(1,11),A2d(1,13))
        enddo
      elseif (t3d_dif4_gp) then
        do itrc=1,nt
          call t3dmix_dif4_gp (itrc,A3d(1,1),
     &                              A2d(1,1),A2d(1,2),
     &                              A2d(1,3),A2d(1,5),
     &                              A2d(1,7),A2d(1,9),
     &                              A2d(1,11),A2d(1,13))
        enddo
      endif

      do itrc=1,nt-1,2
        call exchange_xxx( t(:,:,:,nnew,itrc),
     &                     t(:,:,:,nnew,itrc+1) )
      enddo
      if (mod(nt,2) == 1) then
        call exchange_xxx( t(:,:,:,nnew,nt) )
      endif


      end subroutine t3dmix
!------------------------------------------------------------------------------
      subroutine t3dmix_dif2_gp (itrc,FX,FE,FSC,dTdz,
     &                                dTdx,dTde,dZdx,dZde) ![
      use param
      use grid
      use mixing
      use ocean_vars
      use scalars
      use tracers
      use diagnostics
      use roms_mpi

      implicit none

      ! Import/export
      integer,intent(in) :: itrc
      real, dimension(GLOBAL_2D_ARRAY)   :: FX,FE
      real, dimension(GLOBAL_2D_ARRAY,2) :: FSC,dTdz,
     &                                      dTdx,dTde,
     &                                      dZdx,dZde
      ! Local
      integer i,j,k,k1,k2
      real cff, cff1, cff2, cff3, cff4

                                             !!! WARNING: RECURSIVE
! Compute horizontal (along geopotential     !!!   BLOCKING SEQUENCE
! surfaces) Laplacian diffusion terms.       !!!

      k2=1            ! vertical   dTdx,dTde(:,:,k2) k+1   rho-points
      do k=0,nz,+1    ! placement    dTdz,FSC(:,:,k2) k+1/2   W-points
        k1=k2         ! of scratch dTdx,dTde(:,:,k1) k     rho-points
        k2=3-k1       ! variables    dTdz,FSC(:,:,k1) k-1/2   W-points
        if (k<nz) then
          do j=1,ny
            do i=1,nx+1
              cff=0.5*(pm(i,j)+pm(i-1,j))*umask(i,j)
              dZdx(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i-1,j,k+1,nrhs,itrc)
     &                                                      )
            enddo
          enddo
          do j=1,ny+1
            do i=1,nx
              cff=0.5*(pn(i,j)+pn(i,j-1))*vmask(i,j)
              dZde(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i,j-1,k+1))
              dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i,j-1,k+1,nrhs,itrc)
     &                                                      )
            enddo
          enddo
        endif

        if (k==0 .or. k==nz) then
          do j=0,ny+1
            do i=0,nx+1
              dTdz(i,j,k2)=0.0     ! Set no-flux bottom or
              FSC(i,j,k2)=0.0      ! top boundary conditions.
            enddo
          enddo
        else
          do j=0,ny+1
            do i=0,nx+1
              dTdz(i,j,k2)=(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc)
     &                                   )/(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif

        if (k>0) then
          do j=1,ny
            do i=1,nx+1
              FX(i,j)=0.25*(diff2(i,j,itrc)+diff2(i-1,j,itrc))
     &                     *(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             min(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &            +max(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=1,ny+1
            do i=1,nx
              FE(i,j)=0.25*(diff2(i,j,itrc)+diff2(i,j-1,itrc))
     &                     *(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &         dTde(i,j,k1)-0.5*(
     &             min(dZde(i,j,k1),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dZde(i,j,k1),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k<nz) then
            do j=1,ny
              do i=1,nx

                cff1=min(dZdx(i  ,j,k1),0.)
                cff2=min(dZdx(i+1,j,k2),0.)

                cff3=max(dZdx(i  ,j,k2),0.)
                cff4=max(dZdx(i+1,j,k1),0.)

                FSC(i,j,k2)=0.5*diff2(i,j,itrc)*(
     &                       cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))
     &                        +cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))

     &                        +cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))
     &                        +cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))
     &                                                               )
                cff1=min(dZde(i,j  ,k1),0.)
                cff2=min(dZde(i,j+1,k2),0.)

                cff3=max(dZde(i,j  ,k2),0.)
                cff4=max(dZde(i,j+1,k1),0.)

                FSC(i,j,k2)=FSC(i,j,k2)+0.5*diff2(i,j,itrc)*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))
     &                                                              )
              enddo
            enddo
          endif
          do j=1,ny
            do i=1,nx
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*( pm(i,j)*
     &               pn(i,j)*( FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j) )
     &                                       +FSC(i,j,k2)-FSC(i,j,k1)
     &                                                   )/Hz(i,j,k)
            enddo
          enddo
        endif  ! <-- k>0
      enddo   ! --> k

      if (diag_trc.and.wrt_t_dia(itrc) ) then
        print *,'diagnostics might not be working here'
        stop
!       call set_diags_t_h_mix(itrc)
      endif

      end subroutine t3dmix_dif2_gp !]
!------------------------------------------------------------------------------
      subroutine t3dmix_dif4_gp (itrc,LapT,FX,FE,FSC,dTdz,
     &                                dTdx,dTde,dZdx,dZde) ![

      use param
      use grid
      use mixing
      use ocean_vars
      use scalars
      use tracers
      use diagnostics
      use roms_mpi

      implicit none

      ! Import/export
      integer,intent(in) :: itrc
      real, dimension(GLOBAL_2D_ARRAY,0:nz) :: LapT
      real, dimension(GLOBAL_2D_ARRAY)      :: FX,FE
      real, dimension(GLOBAL_2D_ARRAY,2)    :: FSC,dTdz,
     &                                         dTdx,dTde,
     &                                         dZdx,dZde
      ! Local
      integer :: i,j,k,k1,k2
      real cff, cff1, cff2, cff3, cff4

! Compute horizontal biharmonic diffusion      !!! WARNING: RECURSIVE
! along geopotential surfaces. The biharmonic  !!!  BLOCKING SEQUENCE
! operator is computed by applying the rotated !!!
! Laplacian operator twice.

! THE FIRST LAPLACIAN: ! vertical   dTdx,dTde(:,:,k2)  k+1   rho-points
!                      ! placement    dTdz,FSC(:,:,k2) k+1/2   W-points
      k2=1             ! of scratch dTdx,dTde(:,:,k1)  k     rho-points
      do k=0,nz,+1     ! variables    dTdz,FSC(:,:,k1) k-1/2   W-points
        k1=k2
        k2=3-k1
        if (k<nz) then
          do j=0,ny+1
            do i=0,nx+2
              cff=0.5*(pm(i,j)+pm(i-1,j))*umask(i,j)
              dZdx(i,j,k2)=cff*(z_r(i,j,k+1)-z_r(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-
     &                               t(i-1,j,k+1,nrhs,itrc))
            enddo
          enddo
          do j=0,ny+2
            do i=0,nx+1
              cff=0.5*(pn(i,j)+pn(i,j-1))*vmask(i,j)
              dZde(i,j,k2)=cff*(z_r(i,j,k+1)-z_r(i,j-1,k+1))
              dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-
     &                               t(i,j-1,k+1,nrhs,itrc))
            enddo
          enddo
        endif

        if (k==0 .or. k==nz) then
          do j=-1,ny+2
            do i=-1,nx+2
              dTdz(i,j,k2)=0.0
              FSC(i,j,k2)=0.0
            enddo
          enddo
        else
          do j=-1,ny+2
            do i=-1,nx+2
              dTdz(i,j,k2)=(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc))
     &                                    /(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif

        if (k>0) then
          do j=0,ny+1
            do i=0,nx+2
              FX(i,j)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &                                                dTdx(i  ,j,k1)
     &      -0.5*( min(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &            +max(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=0,ny+2
            do i=0,nx+1
              FE(i,j)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &                                                  dTde(i,j,k1)
     &      -0.5*( min(dZde(i,j,k1),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dZde(i,j,k1),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k<nz) then
            do j=0,ny+1
              do i=0,nx+1
                cff1=min(dZdx(i  ,j,k1),0.)
                cff2=min(dZdx(i+1,j,k2),0.)

                cff3=max(dZdx(i  ,j,k2),0.)
                cff4=max(dZdx(i+1,j,k1),0.)

                FSC(i,j,k2)=0.5*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))
     &                                                              )
                cff1=min(dZde(i,j  ,k1),0.)
                cff2=min(dZde(i,j+1,k2),0.)

                cff3=max(dZde(i,j  ,k2),0.)
                cff4=max(dZde(i,j+1,k1),0.)

                FSC(i,j,k2)=FSC(i,j,k2)+0.5*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))
     &                                                              )
              enddo
            enddo
          endif
          do j=0,ny+1
            do i=0,nx+1
              LapT(i,j,k)=( pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                       +FE(i,j+1)-FE(i,j) )
     &                                   +FSC(i,j,k2)-FSC(i,j,k1)
     &                                               )/Hz(i,j,k)
            enddo
          enddo
        endif   ! <-- k>0
      enddo   ! --> k

! Apply lateral boundary conditions to the first Laplacian
! (closed or gradient depending on obc_*).

        if (west_bnd) then
          do k=1,nz
            do j=0,ny+1
              if (obc_west) then
                LapT(0,j,k)=LapT(1,j,k)
              else
                LapT(0,j,k)=0.
              endif
            enddo
          enddo
        endif
        if (east_bnd) then
          do k=1,nz
            do j=0,ny+1
              if (obc_east) then
                LapT(nx+1,j,k)=LapT(nx,j,k)
              else
                LapT(nx+1,j,k)=0.
              endif
            enddo
          enddo
        endif

        if (south_bnd) then
          do k=1,nz
            do i=0,nx+1
              if (obc_south) then
                LapT(i,0,k)=LapT(i,1,k)
              else
                LapT(i,0,k)=0.
              endif
            enddo
          enddo
        endif
        if (north_bnd) then
          do k=1,nz
            do i=0,nx+1
              if (obc_north) then
                LapT(i,ny+1,k)=LapT(i,ny,k)
              else
                LapT(i,ny+1,k)=0.
              endif
            enddo
          enddo
        endif

!
! THE SECOND LAPLACIAN
!
      k2=1
      do k=0,nz,+1
       k1=k2
       k2=3-k1

        if (k<nz) then
          do j=1,ny
            do i=1,nx+1
              cff=0.5*(pm(i,j)+pm(i-1,j))*umask(i,j)
              dZdx(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i-1,j,k+1))
            enddo
          enddo
          do j=1,ny+1
            do i=1,nx
              cff=0.5*(pn(i,j)+pn(i,j-1))*vmask(i,j)
              dZde(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i,j-1,k+1))
              dTde(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i,j-1,k+1))
            enddo
          enddo
        endif

        if (k==0 .or. k==nz) then
          do j=0,ny+1
            do i=0,nx+1
              dTdz(i,j,k2)=0.0
              FSC(i,j,k2)=0.0
            enddo
          enddo
        else
          do j=0,ny+1
            do i=0,nx+1
              dTdz(i,j,k2)=(LapT(i,j,k+1)-LapT(i,j,k))
     &                       /(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif

        if (k>0) then
          do j=1,ny
            do i=1,nx+1
              FX(i,j)=-0.25*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
     &                     *(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             min(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &            +max(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=1,ny+1
            do i=1,nx
              FE(i,j)=-0.25*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
     &                     *(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &         dTde(i,j,k1)-0.5*(
     &             min(dZde(i,j,k1),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dZde(i,j,k1),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k<nz) then
            do j=1,ny
              do i=1,nx
                cff1=min(dZdx(i  ,j,k1),0.)
                cff2=min(dZdx(i+1,j,k2),0.)

                cff3=max(dZdx(i  ,j,k2),0.)
                cff4=max(dZdx(i+1,j,k1),0.)

                FSC(i,j,k2)=-0.5*diff4(i,j,itrc)*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))
     &                                                              )
                cff1=min(dZde(i,j  ,k1),0.)
                cff2=min(dZde(i,j+1,k2),0.)

                cff3=max(dZde(i,j  ,k2),0.)
                cff4=max(dZde(i,j+1,k1),0.)

                FSC(i,j,k2)=FSC(i,j,k2)-0.5*diff4(i,j,itrc)*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))
     &                                                              )
              enddo
            enddo
          endif
          do j=1,ny
            do i=1,nx
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) + dt*( pm(i,j)
     &              *pn(i,j)*( FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j) )
     &                        +FSC(i,j,k2)-FSC(i,j,k1)    )/Hz(i,j,k)
            enddo
          enddo
        endif  ! <-- k>0
      enddo   ! --> k

      if (diag_trc.and.wrt_t_dia(itrc) ) then
        print *,'diagnostics might not be working here'
        stop
!       call set_diags_t_h_mix(itrc)
      endif

      end subroutine t3dmix_dif4_gp  !]
