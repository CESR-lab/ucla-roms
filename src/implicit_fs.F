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
! the extrapolated NH trend as before.  Closed domains
! only for now (open boundaries need implicit boundary rows).

#ifdef NH_IMPLICIT_FS
      use param
      use dimensions
      use scalars
      use grid
      use ocean_vars
      use coupling
      use roms_mpi
      use mg_namelist, only: robin_beta, surface_neumann
      use mpi

      implicit none
      private

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
      integer :: ilen, istat, i
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

      DU_avg1=0. ; DV_avg1=0.
      do k=1,nz
        do j=1,ny
          do i=1,nx
            DU_avg1(i,j)=DU_avg1(i,j)+u(i,j,k,nnew)*dn_u(i,j)
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            DV_avg1(i,j)=DV_avg1(i,j)+v(i,j,k,nnew)*dm_v(i,j)
          enddo
        enddo
      enddo
      if (west_bnd)  DU_avg1(1,:)=0.
      if (east_bnd)  DU_avg1(nx+1,:)=0.
      if (south_bnd) DV_avg1(:,1)=0.
      if (north_bnd) DV_avg1(:,ny+1)=0.
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
      call exchange_xxx(zeta(:,:,knew))
      Zt_avg1=zeta(:,:,knew)
      call set_depth            ! new z_w, Hz (done by step2d otherwise)
      ! target transports for the next predictor (set_HUV1): AB2 of the
      ! implicit transports, U(n+3/2) ~ 1.5 U(n+1) - 0.5 U(n)
      DU_half=1.5*DU_avg1-0.5*Un ; DV_half=1.5*DV_avg1-0.5*Vn

      end subroutine implicit_fs_finish  !]

#endif /* NH_IMPLICIT_FS */
      end module implicit_fs
