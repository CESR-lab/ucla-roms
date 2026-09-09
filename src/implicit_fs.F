      module implicit_fs

#include "cppdefs.opt"

! Implicit free surface for non-hydrostatic runs (NH_IMPLICIT_FS):
! the surface condition of the NH pressure becomes the Robin condition
!
!     q + beta dq/dz = theta g dt^2 [ (1-theta) w_s(n) + theta w*_s ],
!     beta = theta^2 g dt^2,
!
! obtained from the dynamic condition q_s = theta g (zeta(n+1)-zeta(n))
! (the hydrostatic gradient is evaluated at zeta(n) in prsgrd, no zeta*
! shift) and the theta-weighted kinematic condition
! zeta(n+1) = zeta(n) + dt [(1-theta) w_s(n) + theta w_s(n+1)]:
! theta = 1 is backward Euler, theta = 1/2 Crank-Nicolson (2nd order).
! The projection then advances the surface waves implicitly at every
! scale, so the hydrostatic restoring force and its non-hydrostatic
! counterpart are never evaluated at different time levels; this
! removes the barotropic amplifier of the split-explicit Dirichlet
! coupling (Examples/NH/TANK_CHECKER).  beta -> 0 is Dirichlet,
! beta -> infinity is Neumann.
!
! Implementation: the library (NHMG, mg_namelist%robin_beta, sigtop)
! replaces the Dirichlet top-face factor 1 by dzw/(dzw+beta); the known
! part of the surface source is applied here by rescaling the provisional
! surface flux before the solve.  After the projection the final
! transports U(n+1) define zeta(n+1) = zeta(n) - dt div U_theta exactly
! with U_theta = theta U(n+1) + (1-theta) U(n) = DU_avg2 (tracer fluxes)
! and DU_avg1 = U(n+1), so step3d_uv2, set_depth and
! the tracer step work unchanged.  The split-explicit barotropic step is
! bypassed (main.F); set_HUV1 replaces the predictor's depth-mean
! transports by an AB2 extrapolation of the implicit ones (DU_half), as
! the standard code does with the fast-averaged ones; the predictor uses
! the extrapolated NH trend as before.
!
! Open boundaries follow the barotropic option of the split-explicit
! code.  OBC_M2FLATHER (radiating): the face transport at n+1 obeys
! Flather with the elevation of the adjacent cell,
!     F(n+1) = F_ext + n c dn (zeta_1(n+1) - zeta_ext),  c = sqrt(g h),
! n the outward normal sign.  Made local by the barotropic approximation
! of the boundary column, zeta_1(n+1)-zeta_1(n) = pbar/(theta g dt) with
! pbar the column mean of the library pressure: per level the face flux
! correction is du(k) = n Arx(k)/(theta c dt) p(k) (library obcrob_u/v,
! a symmetric diagonal leak in the masked face row) and the provisional
! face transport is set to the explicit Flather value with zeta_1(n)
! before the solve (set_face_transport, u3dbc profile kept).  The
! boundary cell keeps its Robin surface row; the ghost elevation is
! zero-gradient (diagnostic).  OBC_M2SPECIFIED ('clamped'): the adjacent
! cells get a PRESCRIBED surface flux (zeta_ext(n+1)-zeta(n))/dt instead
! of the Robin row (library sfcfac = 0 there) and the face transport is
! set to ubar_ext*D: incoming signals enter, outgoing waves reflect.

      ! module dependencies outside the option guard so that the
      ! dependency generator orders the compilation (the option only
      ! guards the contents)
      use param
      use dimensions
      use scalars
      use grid
      use ocean_vars
      use coupling
      use roms_mpi
      use boundary
#ifdef NH_IMPLICIT_FS
      use mg_namelist, only: robin_beta, surface_neumann, sfcfac,
     &                       obcrob_u, obcrob_v
      use mpi
#endif

      implicit none
      private
#ifdef NH_IMPLICIT_FS

      real, public :: theta_fs = 0.5    ! implicitness of the surface; env THETA_FS
      real, allocatable, dimension(:,:) :: Un, Vn   ! transports at time n
      real, allocatable, dimension(:,:), public :: DU_half, DV_half  ! predictor target transports

      public :: init_implicit_fs, implicit_fs_surface_flux,
     &          implicit_fs_finish

      contains

!----------------------------------------------------------------------
      subroutine init_implicit_fs  ![
      implicit none
      character(len=16) :: cenv
      integer :: ilen, istat, i, j
      call get_environment_variable('THETA_FS',cenv,ilen,istat)
      if (istat==0 .and. ilen>0) read(cenv,*) theta_fs
      if (surface_neumann) then
        if (mynode==0) write(*,*) 'implicit_fs: NH_IMPLICIT_FS needs',
     &    ' surface_neumann = .false. (the Robin surface replaces Dirichlet)'
        call MPI_Abort(MPI_COMM_WORLD,1,i)
      endif
      robin_beta = theta_fs*theta_fs*g*dt*dt
      allocate( Un(GLOBAL_2D_ARRAY) ) ; Un=0.
      allocate( Vn(GLOBAL_2D_ARRAY) ) ; Vn=0.
      allocate( DU_half(GLOBAL_2D_ARRAY) ) ; DU_half=DU_avg1
      allocate( DV_half(GLOBAL_2D_ARRAY) ) ; DV_half=DV_avg1
#ifdef OBC_M2FLATHER
      ! radiating open faces: leak coefficient n/(theta c dt) per column
      if (.not.allocated(obcrob_u)) allocate( obcrob_u(0:ny+1,0:nx+1) )
      if (.not.allocated(obcrob_v)) allocate( obcrob_v(0:ny+1,0:nx+1) )
      obcrob_u=0. ; obcrob_v=0.
      if (west_bnd.and.obc_west) then
        do j=1,ny
          obcrob_u(j,1)=-1./(theta_fs*dt*sqrt(g*0.5*(h(0,j)+h(1,j))))
        enddo
      endif
      if (east_bnd.and.obc_east) then
        do j=1,ny
          obcrob_u(j,nx+1)=1./(theta_fs*dt*sqrt(g*0.5*(h(nx,j)+h(nx+1,j))))
        enddo
      endif
      if (south_bnd.and.obc_south) then
        do i=1,nx
          obcrob_v(1,i)=-1./(theta_fs*dt*sqrt(g*0.5*(h(i,0)+h(i,1))))
        enddo
      endif
      if (north_bnd.and.obc_north) then
        do i=1,nx
          obcrob_v(ny+1,i)=1./(theta_fs*dt*sqrt(g*0.5*(h(i,ny)+h(i,ny+1))))
        enddo
      endif
      if (mynode==0 .and. (obc_west.or.obc_east.or.obc_south.or.obc_north))
     &  write(*,*) ' implicit_fs: open boundaries radiating (implicit Flather)'
#else
      ! clamped open boundaries: prescribed surface flux in the adjacent cells
      if (.not.allocated(sfcfac)) allocate( sfcfac(0:ny+1,0:nx+1) )
      sfcfac=1.
      if (west_bnd .and.obc_west ) sfcfac(:,1 )=0.
      if (east_bnd .and.obc_east ) sfcfac(:,nx)=0.
      if (south_bnd.and.obc_south) sfcfac(1 ,:)=0.
      if (north_bnd.and.obc_north) sfcfac(ny,:)=0.
      if (mynode==0 .and. (obc_west.or.obc_east.or.obc_south.or.obc_north))
     &  write(*,*) ' implicit_fs: open boundaries clamped (zeta_ext, ubar_ext)'
#endif
      if (mynode==0) write(*,'(A,F6.3,A,E12.4)')
     &   ' implicit_fs: theta_fs =',theta_fs,'  robin_beta [m] =',robin_beta
      end subroutine init_implicit_fs  !]

!----------------------------------------------------------------------
      subroutine implicit_fs_surface_flux  ![
      ! corrector, just before the projection: w(:,:,nz,nnew) holds the
      ! provisional surface volume flux W*_s.  The Robin source is
      ! f = theta g dt^2 [(1-theta) w_s(n) + theta w*_s]; its known part
      ! moves to the right-hand side:  W_rhs = sigma W*_s
      !        - theta (1-theta) g dt^2/(dzw+beta) * W_s(n),
      ! with sigma = dzw/(dzw+beta) and W_s(n) = -div U(n) (kinematic).
      implicit none
      integer :: i,j
      real :: dzwtop, th, wsn
      real, dimension(0:max(nx,ny)+1) :: ubt
      th=theta_fs
      Un=DU_avg1 ; Vn=DV_avg1              ! U(n), kept for finish
      do j=1,ny
        do i=1,nx
          dzwtop=0.5*Hz(i,j,nz)
          wsn=-( Un(i+1,j)-Un(i,j)+Vn(i,j+1)-Vn(i,j) )
          w(i,j,nz,nnew)=( w(i,j,nz,nnew)*dzwtop
     &                    -th*(1.-th)*g*dt*dt*wsn )/(dzwtop+robin_beta)
        enddo
      enddo
#ifdef OBC_M2FLATHER
      ! radiating: depth-mean of the provisional face flux = explicit
      ! Flather with zeta_1(n); the projection adds n kappa p (library),
      ! so the final transport is Flather at zeta_1(n+1)
      if (west_bnd.and.obc_west) then
        do j=1,ny
          ubt(j)=ubar_west(j) -sqrt(g/(0.5*(h(0,j)+h(1,j))))
     &                          *(zeta(1,j,kstp)-zeta_west(j))
        enddo
        call set_face_transport(1,1,ny,ubt,.true.)
      endif
      if (east_bnd.and.obc_east) then
        do j=1,ny
          ubt(j)=ubar_east(j) +sqrt(g/(0.5*(h(nx,j)+h(nx+1,j))))
     &                          *(zeta(nx,j,kstp)-zeta_east(j))
        enddo
        call set_face_transport(nx+1,1,ny,ubt,.true.)
      endif
      if (south_bnd.and.obc_south) then
        do i=1,nx
          ubt(i)=vbar_south(i) -sqrt(g/(0.5*(h(i,0)+h(i,1))))
     &                           *(zeta(i,1,kstp)-zeta_south(i))
        enddo
        call set_face_transport(1,1,nx,ubt,.false.)
      endif
      if (north_bnd.and.obc_north) then
        do i=1,nx
          ubt(i)=vbar_north(i) +sqrt(g/(0.5*(h(i,ny)+h(i,ny+1))))
     &                           *(zeta(i,ny,kstp)-zeta_north(i))
        enddo
        call set_face_transport(ny+1,1,nx,ubt,.false.)
      endif
#else
      ! clamped: prescribed surface flux in the adjacent cells and
      ! ubar_ext*D at the boundary faces (u,v hold fluxes here)
      if (west_bnd.and.obc_west) then
        do j=1,ny
          w(1,j,nz,nnew)=(zeta_west(j)-zeta(1,j,kstp))/(dt*pm(1,j)*pn(1,j))
        enddo
        call set_face_transport(1,1,ny,ubar_west,.true.)
      endif
      if (east_bnd.and.obc_east) then
        do j=1,ny
          w(nx,j,nz,nnew)=(zeta_east(j)-zeta(nx,j,kstp))
     &                                        /(dt*pm(nx,j)*pn(nx,j))
        enddo
        call set_face_transport(nx+1,1,ny,ubar_east,.true.)
      endif
      if (south_bnd.and.obc_south) then
        do i=1,nx
          w(i,1,nz,nnew)=(zeta_south(i)-zeta(i,1,kstp))/(dt*pm(i,1)*pn(i,1))
        enddo
        call set_face_transport(1,1,nx,vbar_south,.false.)
      endif
      if (north_bnd.and.obc_north) then
        do i=1,nx
          w(i,ny,nz,nnew)=(zeta_north(i)-zeta(i,ny,kstp))
     &                                        /(dt*pm(i,ny)*pn(i,ny))
        enddo
        call set_face_transport(ny+1,1,nx,vbar_north,.false.)
      endif
#endif
      end subroutine implicit_fs_surface_flux  !]

!----------------------------------------------------------------------
      subroutine implicit_fs_finish  ![
      ! corrector, after the projection: u,v(nnew) hold dz*u, dz*v.
      ! Final transports -> DU_avg1; U_theta -> DU_avg2, and
      ! zeta(knew) = zeta(kstp) - dt div(DU_avg2) exactly.
      implicit none
      integer :: i,j,k
      real :: th
      th=theta_fs

      DU_avg1=0. ; DV_avg1=0.        ! faces 1..nx+1: the physical east/north
      do k=1,nz                      ! faces hold the boundary values, the
        do j=1,ny                    ! interior-rank ghosts are exchanged below
          do i=1,nx+1
            DU_avg1(i,j)=DU_avg1(i,j)+u(i,j,k,nnew)*dn_u(i,j)
          enddo
        enddo
        do j=1,ny+1
          do i=1,nx
            DV_avg1(i,j)=DV_avg1(i,j)+v(i,j,k,nnew)*dm_v(i,j)
          enddo
        enddo
      enddo
      if (west_bnd .and..not.obc_west ) DU_avg1(1,:)=0.
      if (east_bnd .and..not.obc_east ) DU_avg1(nx+1,:)=0.
      if (south_bnd.and..not.obc_south) DV_avg1(:,1)=0.
      if (north_bnd.and..not.obc_north) DV_avg1(:,ny+1)=0.
      call exchange_xxx(DU_avg1,DV_avg1)
      DU_avg2=th*DU_avg1+(1.-th)*Un ; DV_avg2=th*DV_avg1+(1.-th)*Vn

      do j=1,ny
        do i=1,nx
          zeta(i,j,knew)=zeta(i,j,kstp) - dt*pm(i,j)*pn(i,j)*(
     &                       DU_avg2(i+1,j)-DU_avg2(i,j)
     &                      +DV_avg2(i,j+1)-DV_avg2(i,j) )
        enddo
      enddo
      call zetabc(zeta(:,:,knew))
#ifdef OBC_M2FLATHER
      ! ghost elevation is diagnostic under the radiating faces
      if (west_bnd .and.obc_west ) zeta(0   ,1:ny,knew)=zeta(1 ,1:ny,knew)
      if (east_bnd .and.obc_east ) zeta(nx+1,1:ny,knew)=zeta(nx,1:ny,knew)
      if (south_bnd.and.obc_south) zeta(1:nx,0   ,knew)=zeta(1:nx,1 ,knew)
      if (north_bnd.and.obc_north) zeta(1:nx,ny+1,knew)=zeta(1:nx,ny,knew)
#else
      if (west_bnd .and.obc_west ) zeta(0   ,1:ny,knew)=zeta_west(1:ny)
      if (east_bnd .and.obc_east ) zeta(nx+1,1:ny,knew)=zeta_east(1:ny)
      if (south_bnd.and.obc_south) zeta(1:nx,0   ,knew)=zeta_south(1:nx)
      if (north_bnd.and.obc_north) zeta(1:nx,ny+1,knew)=zeta_north(1:nx)
#endif
      call exchange_xxx(zeta(:,:,knew))
      Zt_avg1=zeta(:,:,knew)
      call set_depth            ! new z_w, Hz (done by step2d otherwise)
      ! target transports for the next predictor (set_HUV1): AB2 of the
      ! implicit transports, U(n+3/2) ~ 1.5 U(n+1) - 0.5 U(n)
      DU_half=1.5*DU_avg1-0.5*Un ; DV_half=1.5*DV_avg1-0.5*Vn

      end subroutine implicit_fs_finish  !]

!----------------------------------------------------------------------
      subroutine set_face_transport(ib,l0,l1,ubext,isu)  ![
      ! Replace the depth-mean of the provisional boundary-face flux by
      ! ubar_ext*D*dn (u faces, isu) or vbar_ext*D*dm (v faces), keeping
      ! the vertical profile; u,v(nnew) are in flux units (m^3/s) here.
      implicit none
      integer, intent(in) :: ib,l0,l1
      real, dimension(0:), intent(in) :: ubext
      logical, intent(in) :: isu
      integer :: l,k
      real :: Uf,D,tgt,dz
      do l=l0,l1
        if (isu) then
          D=0.5*(h(ib,l)+h(ib-1,l)+zeta(ib,l,kstp)+zeta(ib-1,l,kstp))
          tgt=ubext(l)*D*dn_u(ib,l)
          Uf=0.
          do k=1,nz
            Uf=Uf+u(ib,l,k,nnew)
          enddo
          do k=1,nz
            dz=0.5*(Hz(ib,l,k)+Hz(ib-1,l,k))
            u(ib,l,k,nnew)=u(ib,l,k,nnew)+(tgt-Uf)*dz/D
          enddo
        else
          D=0.5*(h(l,ib)+h(l,ib-1)+zeta(l,ib,kstp)+zeta(l,ib-1,kstp))
          tgt=ubext(l)*D*dm_v(l,ib)
          Uf=0.
          do k=1,nz
            Uf=Uf+v(l,ib,k,nnew)
          enddo
          do k=1,nz
            dz=0.5*(Hz(l,ib,k)+Hz(l,ib-1,k))
            v(l,ib,k,nnew)=v(l,ib,k,nnew)+(tgt-Uf)*dz/D
          enddo
        endif
      enddo
      end subroutine set_face_transport  !]

#endif /* NH_IMPLICIT_FS */
      end module implicit_fs
