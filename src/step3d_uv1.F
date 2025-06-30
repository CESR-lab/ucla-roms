#include "cppdefs.opt"
#ifdef SOLVE3D
# define UPSTREAM_UV

      subroutine step3d_uv1(tile)
      use param
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer tile

# include "compute_tile_bounds.h"
      call step3d_uv1_tile( istr,iend,jstr,jend,    A3d(1,1),A3d(1,2),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                              A2d(1,5),A2d(1,6)
# ifdef NHMG
     &                                             ,A3d(1,5)
# endif
     &                                                              )
      end

      subroutine step3d_uv1_tile( istr,iend,jstr,jend, ru,rv,
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

      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
      real cff
      real, parameter ::  delta=0.1666666666666667 ! delta=0.125
# ifdef UPSTREAM_UV
     &                  , gamma=0.3333333333333333 ! gamma=0.25
# endif
# ifdef NHMG
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: rw
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

# include "compute_auxiliary_bounds.h"

# ifdef DIAGNOSTICS
#  ifdef NHMGDIAG
      if (diag_uv.and.calc_diag) then
        Wdiag(istr:iend,jstr:jend,:,iwprsgr) =rw(istr:iend,jstr:jend,:)
      endif
#  endif
# endif

# ifdef NHMG

#  define FX UFx
#  define FE VFe

      do k=1,N
#  include "compute_horiz_rhs_w_terms.h"
      enddo

#  undef FE
#  undef FX

      do j=jstr,jend

#  include "compute_vert_rhs_w_terms.h"
        do i=istr,iend
          DC(i,0)=dt*pm(i,j)*pn(i,j)
        enddo
        do k=1,N           !! w is dz*w after this
          do i=istr,iend
            w(i,j,k,nnew)=w(i,j,k,nnew) +DC(i,0)*rw(i,j,k)
          enddo
        enddo

      enddo
#  ifdef DIAGNOSTICS_NHMG
      if (calc_diag) call set_diags_w_at_uv1( istr, iend, jstr, jend, DC, rw )
#  endif
# endif  /* NHMG */

      do k=1,N
# include "compute_horiz_rhs_uv_terms.h"
      enddo

# ifdef WEC
      call compute_wec_rhs_uv_terms(istr,iend,jstr,jend,ru,rv,wrk1,wrk2)
# endif

# ifdef DIAGNOSTICS
        if (diag_uv .and. calc_diag) then
          call set_diags_u_4th_adv
          call set_diags_v_4th_adv
        endif
# endif

      do j=jstr,jend

# include "compute_vert_rhs_uv_terms.h"

! Apply wave bottom streaming stress as a body force in the bottom-most
! grid cells only (if !defined BODY_FRICTION)
!! JM move this bit to inside compute_wec_rhs_uv_terms...
# ifdef WECXXX
#  if defined WAVE_FRICTION && !defined BODY_FRICTION
        !! JM add this somehow like a bustr in a manner like below
        !! or move into compute_wec_rhs_uv.
        do i=istrU,iend
          ru(i,j,1) =ru(i,j,1) + dm_u(i,j)*dn_u(i,j)*frc2dx(i,j)
        enddo
        do i=istr,iend
          rv(i,j,1) =rv(i,j,1) + dm_v(i,j)*dn_v(i,j)*frc2de(i,j)
        enddo
#  endif
# endif

#define IMPLICIT_BOTTOM_DRAG
# ifdef IMPLICIT_BOTTOM_DRAG
        if (diag_uv) then
          u_prev = u(1:nx,j,:,nnew)
          v_prev = v(1:nx,j,:,nnew)
        endif

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                      /(  Hz(i,j,N  )+Hz(i-1,j,N  )
     &                         +Hz(i,j,N-1)+Hz(i-1,j,N-1))

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i-1,j,N-1))

          cff=1./( 0.5*(Hz(i,j,N)+Hz(i-1,j,N))
     &                   +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*( u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                     +dt*sustr(i,j)
     &                                                  )
        enddo
        do k=N-1,2,-1      !--> forward elimination
          do i=istrU,iend
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
        do i=istrU,iend
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
        do k=2,N,+1          !--> backsubstitution
          do i=istrU,iend
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

            u(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i-1,j,k))

            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo
# ifdef DIAGNOSTICS
        if (diag_uv.and.calc_diag) then
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
# endif

        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                        /(  Hz(i,j,N  )+Hz(i,j-1,N  )
     &                           +Hz(i,j,N-1)+Hz(i,j-1,N-1))

            WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i,j-1,N-1))

            cff=1./( 0.5*(Hz(i,j,N)+Hz(i,j-1,N))
     &                     +FC(i,N-1)-min(WC(i,N-1),0.) )

            CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

            DC(i,N)=cff*( v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                       +dt*svstr(i,j)
     &                                                    )
          enddo
          do k=N-1,2,-1      !--> forward elimination
            do i=istr,iend
              FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                          /(  Hz(i,j,k  )+Hz(i,j-1,k  )
     &                             +Hz(i,j,k-1)+Hz(i,j-1,k-1))

              WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i,j-1,k-1))

              cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                              +FC(i,k-1)-min(WC(i,k-1),0.)
     &                                +FC(i,k)+max(WC(i,k),0.)
     &                       -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                        )
              CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

              DC(i,k)=cff*( v(i,j,k,nnew) +DC(i,0)*rv(i,j,k)
     &                   +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
            enddo
          enddo
          do i=istr,iend
            DC(i,1)=( v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
     &                       +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                            )/( 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
     &                              +0.5*dt*(r_D(i,j)+r_D(i,j-1))
     &                                  +FC(i,1)+max(WC(i,1),0.)
     &                         -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                          )
            v(i,j,1,nnew)=DC(i,1) * 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
            rvfrc(i,j)=rv(i,j,1) +dm_v(i,j)*dn_v(i,j)*(
     &                                        svstr(i,j)
     &              -0.5*(r_D(i,j-1)+r_D(i,j))*DC(i,1) )
          enddo
          do k=2,N,+1          !--> backsubstitution
            do i=istr,iend
              DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

              v(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              rvfrc(i,j)=rvfrc(i,j) + rv(i,j,k)
            enddo
          enddo

# ifdef DIAGNOSTICS
         if (diag_uv.and.calc_diag) then
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
# endif

        endif   !<-- j >= jstrV

# else  /* --> not IMPLICIT_BOTTOM_DRAG */

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          u(i,j,N,nnew)=u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                     +dt*sustr(i,j) !! add brk2dx
          u(i,j,1,nnew)=u(i,j,1,nnew)+DC(i,0)*ru(i,j,1)

          rufrc(i,j)=ru(i,j,1) +ru(i,j,N) +dn_u(i,j)*dm_u(i,j)*(
     &                                                 sustr(i,j)
     &                 -0.5*(r_D(i-1,j)+r_D(i,j))*u(i,j,1,nstp) )!! add brk2dx
        enddo
        do k=2,N-1
          do i=istrU,iend
            u(i,j,k,nnew)=u(i,j,k,nnew)+DC(i,0)*ru(i,j,k)  ! u*dz
            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo
# ifdef DIAGNOSTICS
        if (diag_uv.and.calc_diag) then
          do k=1,nz
            Udiag(:,j,k,idiss) = ru(1:nx,j,k)*dxdyi_u(:,j)*umask(i,j)
     &          - Udiag(:,j,k,ipgr)
     &          - Udiag(:,j,k,icori)
     &          - Udiag(:,j,k,iadv)
          enddo
        endif
# endif

        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            v(i,j,N,nnew)=v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                       +dt*svstr(i,j)
            v(i,j,1,nnew)=v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)

            rvfrc(i,j)=rv(i,j,1) +rv(i,j,N) +dm_v(i,j)*dn_v(i,j)*(
     &                                                   svstr(i,j) ! add brk2de
     &                   -0.5*(r_D(i,j-1)+r_D(i,j))*v(i,j,1,nstp) )
          enddo
          do k=2,N-1
            do i=istr,iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k)
              rvfrc(i,j)=rvfrc(i,j) +rv(i,j,k)
            enddo
          enddo
# ifdef DIAGNOSTICS
        if (diag_uv.and.calc_diag) then
          do k=1,nz
            Vdiag(:,j,k,idiss) = rv(1:nx,j,k)*dxdyi_v(:,j)*vmask(i,j)
     &          - Vdiag(:,j,k,ipgr)
     &          - Vdiag(:,j,k,icori)
     &          - Vdiag(:,j,k,iadv)
          enddo
        endif
# endif
        endif  !<-- if j>=jstrV

# endif /* IMPLICIT_BOTTOM_DRAG */
      enddo  !<-- j
      !! velocities are in the form:  u*dz, v*dz, w*dz

# if defined NHMG

#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            u(istr,j,k,nnew) = u(istr,j,k,nnew)*
     &         0.5*(Hz(istr-1,j,k)+Hz(istr,j,k))
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            u(iend+1,j,k,nnew) = u(iend+1,j,k,nnew)*
     &         0.5*(Hz(iend+1,j,k)+Hz(iend,j,k))
          enddo
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            v(i,jstr,k,nnew) = v(i,jstr,k,nnew)*
     &         0.5*(Hz(i,jstr-1,k)+Hz(i,jstr,k))
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            v(i,jend+1,k,nnew) = v(i,jend+1,k,nnew)*
     &         0.5*(Hz(i,jend+1,k)+Hz(i,jend,k))
          enddo
        enddo
      endif
#  endif


      ! Translate to fluxes, different than in pre_step, we have dz*u
      ! instead of u

      do k=1,N
        do j=jstr,jend
          do i=istr,iend+1
            u(i,j,k,nnew) = u(i,j,k,nnew)*dn_u(i,j)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            v(i,j,k,nnew) = v(i,j,k,nnew)*dm_v(i,j)
          enddo
        enddo
        if (k==N) then
          do j=Jstr,Jend
            do i=istr,iend
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &            /(0.5*(            Hz(i,j,k)))
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=istr,iend
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &          /(0.5*(Hz(i,j,k+1)+Hz(i,j,k)))
            enddo
          enddo
        endif
      enddo

#  ifdef EXCHANGE
      !! Todo: we should be able to avoid this. We only need u and v at the boundaries
      !! Todo: computational ranges for u and v need to be extended to Iend+1, Jend+1
      call exchange_xxx(u(:,:,:,nnew),v(:,:,:,nnew) )
#  endif

      !! only do this for neumann nhmg conditions
      ! Compute a 'barotropic' correction to w such that it matches
      ! the current ubar,vbar divergence (see Molemaker et al., 2018)
      if (surface_neumann) then

        nh_ubar = 0.
        nh_vbar = 0.
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend+1
              nh_ubar(i,j) = nh_ubar(i,j) + u(i,j,k,nnew)
            enddo
          enddo
          do j=Jstr,Jend+1
            do i=Istr,Iend
              nh_vbar(i,j) = nh_vbar(i,j) + v(i,j,k,nnew)
            enddo
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            nh_wcor(i,j) = w(i,j,N,nnew) +
     &         (nh_ubar(i+1,j)-nh_ubar(i,j)+nh_vbar(i,j+1)-nh_vbar(i,j))
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              w(i,j,k,nnew) = w(i,j,k,nnew) - nh_wcor(i,j)
     &                       * (z_w(i,j,k)-z_w(i,j,0))
     &                       / (z_w(i,j,N)-z_w(i,j,0))
            enddo
          enddo
        enddo
      endif

      call nhmg_solve(Lm,Mm,N,halo,padd_X,padd_E,
     &           u(:,:,:,nnew),v(:,:,:,nnew),w(:,:,:,nnew) )


      ! At this point u,v,w are fluxes with the dimension of m^3/s
      do k=1,N
        do j=jstr,jend
          do i=Istr,Iend
            u(i,j,k,nnew) = (u(i,j,k,nnew) + mggrid(1)%du(k,j,i) )*pn_u(i,j)

            rufrc(i,j) = rufrc(i,j) + dn_u(i,j)*mggrid(1)%du(k,j,i)/dt
            nhdu(i,j,k,iprec1) = mggrid(1)%du(k,j,i)/dt !stored as a trend
          enddo
        enddo

        do j=jstr,jend
          do i=Istr,Iend
            v(i,j,k,nnew) = (v(i,j,k,nnew) + mggrid(1)%dv(k,j,i) )*pm_v(i,j)

            rvfrc(i,j) = rvfrc(i,j) + dm_v(i,j)*mggrid(1)%dv(k,j,i)/dt
            nhdv(i,j,k,iprec1) = mggrid(1)%dv(k,j,i)/dt !stored as a trend
          enddo
        enddo

        if (k==N) then
          do j=jstr,jend
            do i=Istr,Iend
              w(i,j,k,nnew) = (w(i,j,k,nnew) + mggrid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))*(0.5*(            Hz(i,j,k)))

              nhdw(i,j,k,iprec1) = mggrid(1)%dw(k+1,j,i)/dt !stored as a trend
            enddo
          enddo
        else
          do j=jstr,jend
            do i=Istr,Iend
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

      subroutine check_step_uv1_switches(ierr)

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      use param
      use strings

      implicit none
      integer ierr, is,ie, lenstr

      ie=lenstr(cpps)
      is=ie+2 ; ie=is+13
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<step3d_uv1.F>'

# include "track_advec_switches.h"

      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end

#else
      subroutine step3d_uv1_empty
      end
#endif
