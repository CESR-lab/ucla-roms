#include "cppdefs.opt"
# define UPSTREAM_UV

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

      implicit none
      integer  i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,nz) :: ru,rv
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:nz) :: WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
      real cff
      real, parameter ::  delta=0.1666666666666667 ! delta=0.125
# ifdef UPSTREAM_UV
     &                  , gamma=0.3333333333333333 ! gamma=0.25
# endif
# ifdef NHMG
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:nz) :: rw
      real  Flxw,Uflxw,Vflxw
      real  ub,ut
      real  div,dmax
!     integer :: ierr !,imx,jmx,kmx
#endif

! Taken from rhs3d35S.F for WEC, used in compute_horiz_rhs_uv_terms.h
! but needs to be declared here in declarations section not in
! the compute_horiz_rhs_uv_terms.h file else wont compile.
# ifdef WEC
      real cff1
# endif

      if (diag_w) then
!       Wdiag(1:nx,1:ny,:,iwprsgr) =rw(1:nx,1:ny,:)
      endif

# ifdef NHMG

#  define FX UFx
#  define FE VFe

      do k=1,nz
#  include "compute_horiz_rhs_w_terms.h"
      enddo

#  undef FE
#  undef FX

      do j=1,ny

#  include "compute_vert_rhs_w_terms.h"
        do i=1,nx
          DC(i,0)=dt*pm(i,j)*pn(i,j)
        enddo
        do k=1,nz           !! w is dz*w after this
          do i=1,nx
            w(i,j,k,nnew)=w(i,j,k,nnew) +DC(i,0)*rw(i,j,k)
          enddo
        enddo

      enddo

      if (diag_w) call set_diags_w_at_uv1( 1, nx, 1, ny, DC, rw )
# endif  /* NHMG */

      do k=1,nz
# include "compute_horiz_rhs_uv_terms.h"
      enddo

# ifdef WEC
      call compute_wec_rhs_uv_terms(1,nx,1,ny,ru,rv,wrk1,wrk2)
# endif

      if (diag_uv) then
        call set_diags_u_4th_adv
        call set_diags_v_4th_adv
      endif

      do j=1,ny

# include "compute_vert_rhs_uv_terms.h"

! Apply wave bottom streaming stress as a body force in the bottom-most
! grid cells only (if !defined BODY_FRICTION)
!! JM move this bit to inside compute_wec_rhs_uv_terms...
# ifdef WECXXX
#  if defined WAVE_FRICTION && !defined BODY_FRICTION
        !! JM add this somehow like a bustr in a manner like below
        !! or move into compute_wec_rhs_uv.
        do i=1,nx
          ru(i,j,1) =ru(i,j,1) + dm_u(i,j)*dn_u(i,j)*frc2dx(i,j)
        enddo
        do i=1,nx
          rv(i,j,1) =rv(i,j,1) + dm_v(i,j)*dn_v(i,j)*frc2de(i,j)
        enddo
#  endif
# endif

        if (diag_uv) then
          u_prev = u(1:nx,j,:,nnew)
          v_prev = v(1:nx,j,:,nnew)
        endif

        do i=1,nx
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        enddo
        do k=1,nz
          do i=1,nx
            DC(i,k)= u(i,j,k,nnew) +DC(i,0)*ru(i,j,k)
            !! add ru into rufrc here as well
          enddo
        enddo

        do i=1,nx
          FC(i,nz-1)= 2.*dt*(Akv(i,j,nz-1)+Akv(i-1,j,nz-1))
     &                      /(  Hz(i,j,nz  )+Hz(i-1,j,nz  )
     &                         +Hz(i,j,nz-1)+Hz(i-1,j,nz-1))

          WC(i,nz-1)= DC(i,0)*0.5*(Wi(i,j,nz-1)+Wi(i-1,j,nz-1))

          cff=1./( 0.5*(Hz(i,j,nz)+Hz(i-1,j,nz))
     &                   +FC(i,nz-1)-min(WC(i,nz-1),0.) )

          CF(i,nz-1)=cff*( FC(i,nz-1)+max(WC(i,nz-1),0.) )

!         DC(i,nz)=cff*( u(i,j,nz,nnew) +DC(i,0)*ru(i,j,nz)
          DC(i,nz)=cff*( DC(i,nz) +dt*sustr(i,j) )
        enddo
        do k=nz-1,2,-1      !--> forward elimination
          do i=1,nx
            FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i-1,j,k-1))
     &                        /(  Hz(i,j,k  )+Hz(i-1,j,k  )
     &                           +Hz(i,j,k-1)+Hz(i-1,j,k-1))

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i-1,j,k-1))

            cff=1./( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                           +FC(i,k-1)-min(WC(i,k-1),0.)
     &                             +FC(i,k)+max(WC(i,k),0.)
     &                    -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                      )
            CF(i,k-1)=cff*(   FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*( u(i,j,k,nnew) +DC(i,0)*ru(i,j,k)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
          enddo
        enddo
        do i=1,nx
          DC(i,1)=( u(i,j,1,nnew) +DC(i,0)*ru(i,j,1)
     &                      +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                          )/( 0.5*(Hz(i,j,1)+Hz(i-1,j,1))
     &                            +0.5*dt*(r_D(i,j)+r_D(i-1,j))
     &                                +FC(i,1)+max(WC(i,1),0.)
     &                       -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                       )

          u(i,j,1,nnew)=DC(i,1) * 0.5*(Hz(i,j,1)+Hz(i-1,j,1))

          rufrc(i,j)=ru(i,j,1) +dm_u(i,j)*dn_u(i,j)*(
     &                                      sustr(i,j)
     &            -0.5*(r_D(i-1,j)+r_D(i,j))*DC(i,1) )
        enddo
        do k=2,nz,+1          !--> backsubstitution
          do i=1,nx
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

            u(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i-1,j,k))

            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo

        if (diag_uv) then
          do k=1,nz
           do i=1,nx
            Udiag(i,j,k,idiss) = ru(i,j,k)*dxdyi_u(i,j)*umask(i,j)
     &          - Udiag(i,j,k,ipgr)
     &          - Udiag(i,j,k,icori)
     &          - Udiag(i,j,k,iadv)
            Udiag(i,j,k,ivmix) = (u(i,j,k,nnew)-u_prev(i,k))/dt -
     &              ru(i,j,k)*dxdyi_u(i,j)
           enddo
          enddo
        endif

        do i=1,nx
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

          FC(i,nz-1)= 2.*dt*(Akv(i,j,nz-1)+Akv(i,j-1,nz-1))
     &                      /(  Hz(i,j,nz  )+Hz(i,j-1,nz  )
     &                         +Hz(i,j,nz-1)+Hz(i,j-1,nz-1))

          WC(i,nz-1)= DC(i,0)*0.5*(Wi(i,j,nz-1)+Wi(i,j-1,nz-1))

          cff=1./( 0.5*(Hz(i,j,nz)+Hz(i,j-1,nz))
     &                   +FC(i,nz-1)-min(WC(i,nz-1),0.) )

          CF(i,nz-1)=cff*( FC(i,nz-1)+max(WC(i,nz-1),0.) )

          DC(i,nz)=cff*( v(i,j,nz,nnew) +DC(i,0)*rv(i,j,nz)
     &                                     +dt*svstr(i,j)
     &                                                  )
        enddo
        do k=nz-1,2,-1      !--> forward elimination
          do i=1,nx
            FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                        /(  Hz(i,j,k  )+Hz(i,j-1,k  )
     &                           +Hz(i,j,k-1)+Hz(i,j-1,k-1))

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i,j-1,k-1))

            cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                           +FC(i,k-1)-min(WC(i,k-1),0.)
     &                             +FC(i,k)+max(WC(i,k),0.)
     &                    -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                      )
            CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*( v(i,j,k,nnew) +DC(i,0)*rv(i,j,k)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
          enddo
        enddo
        do i=1,nx
          DC(i,1)=( v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
     &                      +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                          )/( 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
     &                            +0.5*dt*(r_D(i,j)+r_D(i,j-1))
     &                                +FC(i,1)+max(WC(i,1),0.)
     &                       -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                       )
          v(i,j,1,nnew)=DC(i,1) * 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
          rvfrc(i,j)=rv(i,j,1) +dm_v(i,j)*dn_v(i,j)*(
     &                                      svstr(i,j)
     &            -0.5*(r_D(i,j-1)+r_D(i,j))*DC(i,1) )
        enddo
        do k=2,nz,+1          !--> backsubstitution
          do i=1,nx
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

            v(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
            rvfrc(i,j)=rvfrc(i,j) + rv(i,j,k)
          enddo
        enddo

        if (diag_uv) then
          do k=1,nz
           do i=1,nx
            Vdiag(i,j,k,idiss) = rv(i,j,k)*dxdyi_v(i,j)*vmask(i,j)
     &          - Vdiag(i,j,k,ipgr)
     &          - Vdiag(i,j,k,icori)
     &          - Vdiag(i,j,k,iadv)
            Vdiag(i,j,k,ivmix) = (v(i,j,k,nnew)-v_prev(i,k))/dt -
     &              rv(i,j,k)*dxdyi_v(i,j)
           enddo
          enddo
        endif

      enddo  !<-- j
      !! velocities are in the form:  u*dz, v*dz, w*dz

# if defined NHMG

      if (west_bnd) then
        do k=1,nz
          do j=1,ny
            u(1,j,k,nnew) = u(1,j,k,nnew)*
     &         0.5*(Hz(0,j,k)+Hz(1,j,k))
          enddo
        enddo
      endif
      if (east_bnd) then
        do k=1,nz
          do j=1,ny
            u(nx+1,j,k,nnew) = u(nx+1,j,k,nnew)*
     &         0.5*(Hz(nx+1,j,k)+Hz(nx,j,k))
          enddo
        enddo
      endif
      if (south_bnd) then
        do k=1,nz
          do i=1,nx
            v(i,1,k,nnew) = v(i,1,k,nnew)*
     &         0.5*(Hz(i,0,k)+Hz(i,1,k))
          enddo
        enddo
      endif
      if (north_bnd) then
        do k=1,nz
          do i=1,nx
            v(i,ny+1,k,nnew) = v(i,ny+1,k,nnew)*
     &         0.5*(Hz(i,ny+1,k)+Hz(i,ny,k))
          enddo
        enddo
      endif


      ! Translate to fluxes, different than in pre_step, we have dz*u
      ! instead of u

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
