#include "cppdefs.opt"

      subroutine t3dmix
!     Computes horizontal mixing of tracers along S-surfaces.
      use param
      use private_scratch
      use tracers
      use roms_mpi
      use scalars

      implicit none

      integer :: itrc
      logical,parameter :: t3d_dif2_s=.true.
      logical,parameter :: t3d_dif4_s=.false.

      if (.not.t3d_dif2_s.and..not.t3d_dif4_s) then
        ! nothing to do
        return
      endif

      if (t3d_dif2_s.and.t3d_dif4_s) then
        if (mynode.eq.0) 
     &    print *,'only one of t3d_dif2_s or t3d_dif4_s can be true'
        stop
      endif

      if (t3d_dif2_s) then
        do itrc=1,nt
          call t3dmix_dif2_s (itrc,A2d(1,1),A2d(1,2))
        enddo
      elseif (t3d_dif4_s) then
        do itrc=1,nt
          call t3dmix_dif4_s (itrc,A2d(1,1),A2d(1,2),A2d(1,3))
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
      subroutine t3dmix_dif2_s (itrc,FX,FE) ![
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
      real, dimension(GLOBAL_2D_ARRAY) :: FX,FE
      ! Local
      integer i,j,k


      do k=1,nz
! Add in horizontal Laplacian diffusion along constant S-surfaces.
! Compute x- and y-components of diffusive tracer flux.

        do j=1,ny
          do i=1,nx+1
            FX(i,j)=0.25*(diff2(i,j,itrc)+diff2(i-1,j,itrc))
     &                   *pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))*(
     &                     t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc)
     &                                                              )
     &                                                   *umask(i,j)
          enddo
        enddo
        do j=1,ny+1
          do i=1,nx
            FE(i,j)=0.25*(diff2(i,j,itrc)+diff2(i,j-1,itrc))
     &                     *pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))*(
     &                      t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc)
     &                                                              )
     &                                                   *vmask(i,j)
          enddo
        enddo

        if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
          do j=1,ny
            do i=1,nx+1
              Tdiag(i,j,k,tmixx,td(itrc)) =
     &           Tdiag(i,j,k,tmixx,td(itrc)) - FX(i,j)
            enddo
          enddo
          do j=1,ny+1
            do i=1,nx
              Tdiag(i,j,k,tmixy,td(itrc)) =
     &           Tdiag(i,j,k,tmixy,td(itrc)) - FE(i,j)
            enddo
          enddo
        endif

        ! Add in horizontal diffusion of tracer [T/s].
        do j=1,ny
          do i=1,nx
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*pm(i,j)*pn(i,j)
     &                         *(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
     &                                                     /Hz(i,j,k)
          enddo
        enddo

      enddo

      end subroutine t3dmix_dif2_s !]
!------------------------------------------------------------------------------
      subroutine t3dmix_dif4_s(itrc,FX,FE,LapT) ![

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
      real, dimension(GLOBAL_2D_ARRAY) :: FX,FE,LapT
      ! Local
      integer :: i,j,k

!  Compute/Add in horizontal biharmonic diffusion along constant
!  S-surfaces.
!-----------------------------------------------------------------
!  The biharmonic operator is computed by applying the Laplacian
!  operator twice.
!
!  Compute horizontal tracer flux in the XI-direction at U-points.
      do k= 1,nz

        do j=0,ny+1
          do i=0,nx+2
            FX(i,j)=pmon_u(i,j)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))*(
     &                      t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc)
     &                                                              )
     &                                                   *umask(i,j)
          enddo
        enddo

!  Compute horizontal tracer flux in the ETA-direction at V-points.

        do j=0,ny+2
          do i=0,nx+1
            FE(i,j)=pnom_v(i,j)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))*(
     &                      t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc)
     &                                                            )
     &                                                 *vmask(i,j)
          enddo
        enddo

!  Compute the first Laplacian, without mixing coefficient.
!  Multiply by the metrics of the second Laplacian.
!  Save into work array "LapT".

        do j=0,ny+1
          do i=0,nx+1
            LapT(i,j)=( FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j) )
     &                               *pm(i,j)*pn(i,j)/Hz(i,j,k)
          enddo
        enddo

!  Apply boundary conditions (except periodic; closed or gradient)
!  to the first Laplacian.

        if (west_bnd) then
          do j=0,ny+1
            if (obc_west) then
              LapT(0,j)=LapT(1,j)
            else
              LapT(0,j)=0.
            endif
          enddo
        endif
        if (east_bnd) then
          do j=0,ny+1
            if (obc_east) then
              LapT(nx+1,j)=LapT(nx,j)
            else
              LapT(nx+1,j)=0.
            endif
          enddo
        endif

        if (south_bnd) then
          do i=0,nx+1
            if (obc_south) then
              LapT(i,0)=LapT(i,1)
            else
              LapT(i,0)=0.
            endif
          enddo
        endif
        if (north_bnd) then
          do i=0,nx+1
            if (obc_north) then
              LapT(i,ny+1)=LapT(i,ny)
            else
              LapT(i,ny+1)=0.
            endif
          enddo
        endif

! Compute FX=d(LapT)/d(xi) and FE=d(LapT)/d(eta) terms.  After that
! compute and add in biharmonic mixing [T m^3/s]. Multiply by mixing
! coefficient.

        do j=1,ny
          do i=1,nx+1
            FX(i,j)=0.25*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
     &                   *pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))
     &                               *(LapT(i,j)-LapT(i-1,j))
     &                                           *umask(i,j)
          enddo
        enddo
        do j=1,ny+1
          do i=1,nx
            FE(i,j)=0.25*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
     &                   *pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))
     &                               *(LapT(i,j)-LapT(i,j-1))
     &                                           *vmask(i,j)
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*pm(i,j)*pn(i,j)
     &                         *(FX(i,j)-FX(i+1,j)+FE(i,j)-FE(i,j+1))
     &                                                     /Hz(i,j,k)
          enddo
        enddo
      enddo

      end subroutine t3dmix_dif4_s  !]
