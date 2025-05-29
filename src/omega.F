#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine omega

! Compute S-coordinate vertical velocity, w=[Hz/(m*n)]*omega [m^3/s],
! which has the meaning of FINITE_VOLUME FLUX across MOVING grid-box
! interfaces of RHO-boxes. To compute it we first integrate divergence
! of horizontal mass fluxes from bottom up, starting with the no-flow
! boundary condition at the bottom (k=0); After this operation W(:,:,N)
! contains vertical velocity flux at the free surface, which is the
! time-tendency of the free surface, d_zeta/d_t multiplied by grid-box
! area as seen from above;  To convert W(:,:,:) into S-coordinate
! vertical velocity flux, one needs to subtract the vertical velocities
! of moving S-coordinate surfaces, which are proportional the product
! of d_zeta/d_t and the fraction of the distance from the point to the
! bottom divided by the total depth of water column, i.e. the whole
! S-coordinate system is "breathes" by linear in Z-space expansion
! and contraction set by variation in free surface.

! Parameter setting: "cu_min" is threshold value for Courant Number
! below which vertical advection is fully explicit; "cu_max" is the
! maximum CN which the explicit component "We" is allowed to reach.

      use param
      use pipe_frc
      use cdr_frc
      use grid
      use ocean_vars
      use scalars
      use private_scratch
      use hidden_mpi_vars

      implicit none

      ! local
      integer :: i,j,k
      real,dimension(:,:),pointer :: CX
      real,dimension(:),pointer   :: wrk
      real                        :: dtau,c2d,dh,cw,cff,
     &                               cw_min,cw_max,cw_max2
      real,parameter              :: cu_min=0.6D0, cu_max=1.0D0,
     &  cmnx_ratio=cu_min/cu_max,  cutoff=2.D0-cmnx_ratio,
     &                      r4cmx=0.25D0/(1.D0-cmnx_ratio)

      CX(1-bf:nx+bf,1-bf:nz+bf) => scratch_2d(:,1)
      wrk => scratch_1d(:,1)

      if (CORR_STAGE) then
        dtau=dt
      elseif (FIRST_TIME_STEP) then
        dtau=0.5*dt
      else
        dtau=0.6*dt
      endif

      Wi(1:nx,1:ny,0)=0.

      if (cdr_source.and.cdr_volume) then
        do cidx=1,cdr_nprf
          icdr = cdr_icdr(cidx)
          i = cdr_iloc(cidx)
          j = cdr_iloc(cidx)
          do k=1,nz
            Wi(i,j,k)=Wi(i,j,k-1)
     &         + cdr_vol(icdr)*cdr_prf(cidx,k)
          enddo
        enddo
      endif

      do j=jstr,jend

      !!! For NHMG pipes, the flux needs to be out of the bottom.
      !!! Additionally, modify, where needed, the We to Wi split 
      !!! for pipe sources
        do k=1,N,+1        !--> recursive
          do i=1,nx
            Wi(i,j,k)=Wi(i,j,k-1) -FlxU(i+1,j,k) +FlxU(i,j,k)
     &                            -FlxV(i,j+1,k) +FlxV(i,j,k)
            if (pipe_source) then
              if (pipe_idx(i,j)>0) then         ! pipe_source removes branching because it's a parameter
                pidx = pipe_idx(i,j)
                Wi(i,j,k)=Wi(i,j,k)
     &                   + pipe_flx(i,j)*pipe_prf(pidx,k)
              endif
            endif

            CX(i,k)=max(FlxU(i+1,j,k),0.)-min(FlxU(i,j,k),0.)
     &             +max(FlxV(i,j+1,k),0.)-min(FlxV(i,j,k),0.)
          enddo
        enddo
        do i=1,nx
          wrk(i)=Wi(i,j,N)/(z_w(i,j,N)-z_w(i,j,0))
          Wi(i,j,N)=0.
          We(i,j,N)=0.                      ! note that
          We(i,j,0)=0.                      ! CX(i,k)*CX(i,0)/Hz(i,j,k) 
          CX(i,0)=dtau*pm(i,j)*pn(i,j)      ! is horizontal 2D Courant
        enddo                               ! number within Hz(i,j,k)
        do k=N-1,1,-1
          do i=1,nx
            Wi(i,j,k)=Wi(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))  !! this removes the grid motion

c*          if (Wi(i,j,k) > 0.) then         ! Three different variants 
c*            c2d=CX(i,k)   ; dh=Hz(i,j,k)   ! for computing 2D Courant 
c*          else                             ! number at the interface: 
c*            c2d=CX(i,k+1) ; dh=Hz(i,j,k+1) ! (1) use value from the 
c*          endif                            !     grid box upstream in 
                                             !     vertical direction;
c>          c2d=0.5*(CX(i,k) +CX(i,k+1))
c>          dh=0.5*(Hz(i,j,k)+Hz(i,j,k+1))   ! (2) average the two; or

            c2d=max(CX(i,k),  CX(i,k+1))     ! (3) pick the maximum
            dh=min(Hz(i,j,k),Hz(i,j,k+1))    !     of the two.


            cw_max=cu_max*dh-c2d*CX(i,0)
            if (cw_max > 0.D0) then
              cw_max2=cw_max*cw_max
              cw_min=cw_max*cmnx_ratio

              cw=abs(Wi(i,j,k))*CX(i,0)  !<-- cw/dh = vertical Courant

              if (cw < cw_min) then
                cff=cw_max2
              elseif (cw < cutoff*cw_max) then
                cff=cw_max2 +r4cmx*(cw-cw_min)**2 
              else
                cff=cw_max*cw
              endif

              We(i,j,k)=cw_max2*Wi(i,j,k)/cff  !<-- NORMAL SPLITING
              Wi(i,j,k)=Wi(i,j,k)-We(i,j,k)    !<-- CODE TO BE USED
            else
              We(i,j,k)=0.          !--> Wi(i,j,k) remains unchanged
            endif

c**         We(i,j,k)=0.            !<-- fully implicit(BE), test only

c**         We(i,j,k)=Wi(i,j,k)     !<-- for testing only: revert
c**         Wi(i,j,k)=0.            !<-- back to fully explicit code
          enddo
        enddo
      enddo

      ! done with scratch space
      nullify(CX,wrk)

      ! lateral boundary conditions
      if (inode==0.and.(.not.ew_periodic)) then
        do k=0,N
          do j=1,ny
            We(0,j,k)=We(1,j,k)
            Wi(0,j,k)=Wi(1,j,k)
          enddo
        enddo
      endif
      if (inode==np_xi-1.and.(.not.ew_periodic)) then
        do k=0,N
          do j=1,ny
            We(nx+1,j,k)=We(nx,j,k)
            Wi(nx+1,j,k)=Wi(nx,j,k)
          enddo
        enddo
      endif
      if (jnode==0.and.(.not.ns_periodic)) then
        do k=0,N
          do i=1,nx
            We(i,0,k)=We(i,1,k)
            Wi(i,0,k)=Wi(i,1,k)
          enddo
        enddo
      endif
      if (jnode==np_eta-1.and.(.not.ns_periodic)) then
        do k=0,N
          do i=1,nx
            We(i,ny+1,k)=We(i,ny,k)
            Wi(i,ny+1,k)=Wi(i,ny,k)
          enddo
        enddo
      endif

      ! corners
      if ((.not.ew_periodic).and.(.not.ns_periodic)) then
        if (jnode==0.and.inode==0) then
          do k=0,N
            We(0,0,k)=We(1,1,k)
            Wi(0,0,k)=Wi(1,1,k)
          enddo
        endif
        if (jnode==0.and.inode==np_xi-1) then
          do k=0,N
            We(0,ny+1,k)=We(1,ny,k)
            Wi(0,ny+1,k)=Wi(1,ny,k)
          enddo
        endif
        if (jnode==np_eta-1.and.inode==0) then
          do k=0,N
            We(nx+1, 0,k)=We(nx,1,k)
            Wi(nx+1, 0,k)=Wi(nx,1,k)
          enddo
        endif
        if (jnode==np_eta-1.and.inode==np_xi-1) then
          do k=0,N
            We(nx+1,ny+1,k)=We(nx,ny,k)
            Wi(nx+1,ny+1,k)=Wi(nx,ny,k)
          enddo
        endif
      endif
# ifdef EXCHANGE
      call exchange_xxx(We,Wi)
# endif
      end
#else
      subroutine omega_empty
      end
#endif /* SOLVE3D */

