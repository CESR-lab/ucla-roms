#include "cppdefs.opt"

      subroutine step3d_uv1
      use param
      use roms_mpi
      use private_scratch
      implicit none

      call step3d_uv1_tile( A3d(1,1),A3d(1,2),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                              A2d(1,5),A2d(1,6)
# ifdef NHMG
     &                                             ,A3d(1,5)
# endif
     &                                                              )
      end

      subroutine step3d_uv1_tile(ru,rv,
     &                                            WC,FC,CF,DC,
     &                             UFx,UFe,VFx,VFe, wrk1,wrk2
# ifdef NHMG
     &                                            ,rw
# endif
     &                                                              )

      use param
# ifdef NHMG
      use nhmg
      use mg_grids, mggrid => grid
      use mg_tictoc, only : tic, toc
# endif
# ifdef WEC
      use wec_frc
# endif
      use diagnostics
      use dimensions
      use surf_flux, only: sustr, svstr
      use coupling
      use eos_vars
      use grid
      use mixing
      use ocean_vars
      use scalars
      use roms_mpi
      use advection
      use tri_solvers
      use hmixing

      implicit none
      integer  i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,nz) :: ru,rv
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:nz) :: WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
      real cff
# ifdef NHMG
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:nz) :: rw
      real  Flxw,Uflxw,Vflxw
      real  ub,ut
      real  div,dmax
#endif
      real,dimension(nx,0:nz) :: Av,om
      real,dimension(nx,nz) :: dz
      real,dimension(nx) :: drag


      if (diag_w) then
!       Wdiag(1:nx,1:ny,:,iwprsgr) =rw(1:nx,1:ny,:)
      endif

! set the implicit advection and mixing  to zero on the boundaries
      om(:,0)  = 0
      om(:,nz) = 0
      Av(:,0)  = 0
      Av(:,nz) = 0

# ifdef NHMG

      do k=1,nz
        call w_hadv_up3(rw,k)    ! upwinded in the corrector
      enddo

      do j=1,ny

        call w_vadv_up3(rw,j)    ! upwinded in the corrector
        do i=1,nx
          DC(i,0)=dt*pm(i,j)*pn(i,j)
        enddo
        do k=1,nz           !! w is dz*w after this
          do i=1,nx
            w(i,j,k,nnew)=w(i,j,k,nnew) +DC(i,0)*rw(i,j,k)
          enddo
        enddo

        do k=0,nz
          do i=1,nx
            om(i,k)= DC(i,0)*Wi(i,j,k)
            Av(i,k)= Akv(i,j,k)
          enddo
        enddo
        do k=1,nz
          do i=1,nx
            dz(i,k)= Hz(i,j,k)
          enddo
        enddo

        call vmix_solve_w(w(:,:,:,nnew),dt,dz,Av,om,j)

        ! w is expected to be as dz*w
        do k=1,nz-1
          do i=1,nx
            w(i,j,k,nnew) = 0.5*(Hz(i,j,k)+Hz(i,j,k+1))*w(i,j,k,nnew)
          enddo
        enddo
        do i=1,nx
          w(i,j,nz,nnew) = 0.5*Hz(i,j,nz)*w(i,j,nz,nnew)
        enddo

        !------- end computing of w(:,:,:,nnew) -----

      enddo  ! j-loop
# endif /* NHMG */


      do k=1,nz
        call uv_cor_curv(ru,rv,k)
        call uv_hadv_cor(ru,rv,k)
# ifdef UV_VIS2
        call uv_visc3d(ru,rv,k)
# endif
      enddo

# ifdef WEC
      call compute_wec_rhs_uv_terms(1,nx,1,ny,ru,rv,wrk1,wrk2)
# endif

      if (diag_uv) then
        call set_diags_u_4th_adv
        call set_diags_v_4th_adv
      endif

      do j=1,ny

        call uv_vadv_cor(ru,rv,j)

        if (diag_uv) then
          u_prev = u(1:nx,j,:,nnew)
          v_prev = v(1:nx,j,:,nnew)
        endif

        do i=1,nx
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        enddo
        do k=1,nz
          do i=1,nx
            u(i,j,k,nnew) = u(i,j,k,nnew) +DC(i,0)*ru(i,j,k)
          enddo
        enddo
        do i=1,nx
          u(i,j,nz,nnew)= u(i,j,nz,nnew) +dt*sustr(i,j)
        enddo

        do k=1,nz-1
          do i=1,nx
            om(i,k)= DC(i,0)*0.5*( Wi(i,j,k)+Wi(i-1,j,k) )
            Av(i,k)= 0.5*( Akv(i,j,k) + Akv(i-1,j,k) )
          enddo
        enddo
        do k=1,nz
          do i=1,nx
            dz(i,k)= 0.5*( Hz(i,j,k)+Hz(i-1,j,k) )
          enddo
        enddo
        do i=1,nx
          drag(i) = 0.5*(r_D(i,j) + r_D(i-1,j) )
        enddo 

        call vmix_solve(u(:,:,:,nnew),dt,dz,Av,om,j,drag)

        do i=1,nx
          rufrc(i,j)=dm_u(i,j)*dn_u(i,j)*(
     &               sustr(i,j) -drag(i)*u(i,j,1,nnew) )
        enddo
        do k=1,nz
          do i=1,nx
            u(i,j,k,nnew) = u(i,j,k,nnew)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo


! start of v solve
        do i=1,nx
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
        enddo
        do k=1,nz
          do i=1,nx
            v(i,j,k,nnew) = v(i,j,k,nnew) +DC(i,0)*rv(i,j,k)
          enddo
        enddo
        do i=1,nx
          v(i,j,nz,nnew)= v(i,j,nz,nnew) +dt*svstr(i,j)
        enddo

        do k=1,nz-1
          do i=1,nx
            om(i,k)= DC(i,0)*0.5*( Wi(i,j,k)+Wi(i,j-1,k) )
            Av(i,k)= 0.5*( Akv(i,j,k) + Akv(i,j-1,k) )
          enddo
        enddo
        do k=1,nz
          do i=1,nx
            dz(i,k)= 0.5*( Hz(i,j,k) + Hz(i,j-1,k) )
          enddo
        enddo
        do i=1,nx
          drag(i) = 0.5*(r_D(i,j) + r_D(i,j-1) )
        enddo

        call vmix_solve(v(:,:,:,nnew),dt,dz,Av,om,j,drag)

        do i=1,nx
          rvfrc(i,j)=dm_v(i,j)*dn_v(i,j)*(
     &               svstr(i,j) -drag(i)*v(i,j,1,nnew) )
        enddo
        do k=1,nz
          do i=1,nx
            v(i,j,k,nnew) = v(i,j,k,nnew)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
            rvfrc(i,j)=rvfrc(i,j) +rv(i,j,k)
          enddo
        enddo

        if (diag_uv) then
          do k=1,nz
           do i=1,nx
            Udiag(i,j,k,idiss) = ru(i,j,k)*dxdyi_u(i,j)*umask(i,j)
     &          - Udiag(i,j,k,ipgr)
     &          - Udiag(i,j,k,icori)
     &          - Udiag(i,j,k,iadv)
     &          - Udiag(i,j,k,ihmix)
            Udiag(i,j,k,ivmix) = (u(i,j,k,nnew)-u_prev(i,k))/dt -
     &              ru(i,j,k)*dxdyi_u(i,j)

            Vdiag(i,j,k,idiss) = rv(i,j,k)*dxdyi_v(i,j)*vmask(i,j)
     &          - Vdiag(i,j,k,ipgr)
     &          - Vdiag(i,j,k,icori)
     &          - Vdiag(i,j,k,iadv)
     &          - Vdiag(i,j,k,ihmix)
            Vdiag(i,j,k,ivmix) = (v(i,j,k,nnew)-v_prev(i,k))/dt -
     &              rv(i,j,k)*dxdyi_v(i,j)
           enddo
          enddo
        endif

      enddo  !<-- j
      !! velocities are in the form:  u*dz, v*dz, w*dz

      !! Translate to u,v so we can do boundary conditions
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k,nnew)=u(i,j,k,nnew)/(0.5*(Hz(i,j,k)+Hz(i-1,j,k)))
            v(i,j,k,nnew)=v(i,j,k,nnew)/(0.5*(Hz(i,j,k)+Hz(i,j-1,k)))
          enddo
        enddo
      enddo

      call u3dbc(wrk1)
      call v3dbc(wrk1)

      do k=1,nz
        do j=1,ny+1
          do i=1,nx+1
            u(i,j,k,nnew)=u(i,j,k,nnew)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            v(i,j,k,nnew)=v(i,j,k,nnew)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
          enddo
        enddo
      enddo

# if defined NHMG

      ! Translate to fluxes
      do k=1,nz
        do j=1,ny
          do i=1,nx+1
            u(i,j,k,nnew) = u(i,j,k,nnew)*dn_u(i,j)
          enddo
        enddo
        do j=1,ny+1
          do i=1,nx
            v(i,j,k,nnew) = v(i,j,k,nnew)*dm_v(i,j)
          enddo
        enddo
        if (k==nz) then
          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &            /(0.5*(            Hz(i,j,k)))
            enddo
          enddo
        else
          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &          /(0.5*(Hz(i,j,k+1)+Hz(i,j,k)))
            enddo
          enddo
        endif
      enddo

#  ifdef EXCHANGE
      !! Todo: we should be able to avoid this. We only need u and v at the boundaries
      !! Todo: computational ranges for u and v need to be extended to nx+1, ny+1
      call exchange_xxx(u(:,:,:,nnew),v(:,:,:,nnew) )
#  endif

      !! only do this for neumann nhmg conditions
      ! Compute a 'barotropic' correction to w such that it matches
      ! the current ubar,vbar divergence (see Molemaker et al., 2018)
      if (surface_neumann) then

        nh_ubar = 0.
        nh_vbar = 0.
        do k=1,nz
          do j=1,ny
            do i=1,nx+1
              nh_ubar(i,j) = nh_ubar(i,j) + u(i,j,k,nnew)
            enddo
          enddo
          do j=1,ny+1
            do i=1,nx
              nh_vbar(i,j) = nh_vbar(i,j) + v(i,j,k,nnew)
            enddo
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            nh_wcor(i,j) = w(i,j,nz,nnew) +
     &         (nh_ubar(i+1,j)-nh_ubar(i,j)+nh_vbar(i,j+1)-nh_vbar(i,j))
          enddo
        enddo
        do k=1,nz
          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = w(i,j,k,nnew) - nh_wcor(i,j)
     &                       * (z_w(i,j,k)-z_w(i,j,0))
     &                       / (z_w(i,j,nz)-z_w(i,j,0))
            enddo
          enddo
        enddo
      endif

      call nhmg_solve(nx,ny,nz,halo,padd_X,padd_E,
     &           u(:,:,:,nnew),v(:,:,:,nnew),w(:,:,:,nnew) )


      ! At this point u,v,w are fluxes with the dimension of m^3/s
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k,nnew) = (u(i,j,k,nnew) + mggrid(1)%du(k,j,i) )*pn_u(i,j)

            rufrc(i,j) = rufrc(i,j) + dn_u(i,j)*mggrid(1)%du(k,j,i)/dt
            nhdu(i,j,k,iprec1) = mggrid(1)%du(k,j,i)/dt !stored as a trend
          enddo
        enddo

        do j=1,ny
          do i=1,nx
            v(i,j,k,nnew) = (v(i,j,k,nnew) + mggrid(1)%dv(k,j,i) )*pm_v(i,j)

            rvfrc(i,j) = rvfrc(i,j) + dm_v(i,j)*mggrid(1)%dv(k,j,i)/dt
            nhdv(i,j,k,iprec1) = mggrid(1)%dv(k,j,i)/dt !stored as a trend
          enddo
        enddo

        if (k==nz) then
          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = (w(i,j,k,nnew) + mggrid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))*(0.5*(            Hz(i,j,k)))

              nhdw(i,j,k,iprec1) = mggrid(1)%dw(k+1,j,i)/dt !stored as a trend
            enddo
          enddo
        else
          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = (w(i,j,k,nnew) + mggrid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))*(0.5*(Hz(i,j,k+1)+Hz(i,j,k)))

              nhdw(i,j,k,iprec1) = mggrid(1)%dw(k+1,j,i)/dt !stored as a trend
            enddo
          enddo
        endif
      enddo

      ! u,v, and w contain again dz*u,dz*v, and dz*w

      iprec1 = 1+mod(iprec1,2)
      iprec2 = 1+mod(iprec2,2)

# endif /* NHMG */


       !!!! JM not clear what this does either
# undef r_D
# if defined ADV_ISONEUTRAL && defined EXCHANGE
       !!!! JM not clear what this does !!!
      call exchange_xxx(dRdx,dRde)
# endif

      end
