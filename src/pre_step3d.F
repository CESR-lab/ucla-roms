#include "cppdefs.opt"

                                    ! Note that arrays A2d(1,1:4) in
      subroutine pre_step3d        ! the call just below are repeated
      use param
      use private_scratch
      implicit none                 ! within the list of agruments to
                                    ! use the same memory for vertical
                                    ! and horizontal scratch arrays.

      call pre_step3d_tile (A3d(1,1),A3d(1,2),
     &         A3d(1,3),A3d(1,4), A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                             A2d(1,5),A2d(1,6)
# ifdef NHMG
     &                                            ,A3d(1,5)
# endif
     &          )
      end

      subroutine pre_step3d_tile (ru,rv, Hz_bak,
     &                                         Hz_fwd,   WC,FC,CF,DC,
     &                                    UFx,UFe,VFx,VFe, wrk1,wrk2
# ifdef NHMG
     &                                            ,rw
# endif
     &                                                              )

      use dimensions
      use param
      use river_frc
      use wec_frc
#ifdef NHMG
      use nhmg
      use mg_grids, mggrid => grid
      use mg_tictoc
#endif
      use surf_flux, only: sustr, svstr
      use tracers,   only: iTandS, t
      use diagnostics 
      use coupling
      use grid
      use mixing
      use ocean_vars
      use scalars
      use roms_mpi
      use advection
      use tri_solvers

      implicit none
      integer ::  i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,nz) :: ru,rv,Hz_bak,
     &                                                      Hz_fwd
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:nz) ::  WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
#ifdef NHMG
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:nz) :: rw
      real  Flxw,Uflxw,Vflxw
      real  div,dmax
      real  ub,ut
#endif

      real dtau, cf_stp,cf_bak,  cff, FlxDiv
      real, parameter :: AM3_crv=1./6., epsil=1.E-33
      integer indx, itrc

      real,dimension(nx,0:nz) :: Av,om
      real,dimension(nx,nz) :: dz
      real,dimension(nx) :: drag


! Preliminary step: initialize computations of the new time step
! 3D primitive variables.
!
! Start computation of the auxiliary tracer field.
!------ ----------- -- --- --------- ------ ------
! After this stage the resultant t(:,:,:,nnew,:) is time-centered
! halfway between steps n and n+1. A high spatial order, centered,
! non-conservative, but constancy preserving scheme is applied to
! accomplish it.  The constancy preservation property is achieved
! by introducing an artificial continuity equation [a''la Easter,
! 1993], so that the actual advection scheme is still in the flux
! form, however the new-time-step grid box height "Hz_fwd" (see
! code segment just below) has no relation with the true grid-box
! heights determined from the updated free surface (not available
! at this stage yet), but is rather chosen to absorb the
! 3D-divergence of mass fluxes FlxU, FlxV, and W consistently with
! time-stepping algorithm of this preliminary step (recall that
! computation of "Hz_fwd" mimics time step for tracers themselves).
! Later in this code the "Hz_fwd"-field plays the role of new-step
! "Hz" in the updates for tracer and 3D momenta to n+1/2, and it
! does not participate in any further computation.  Hence, division
! by "Hz_fwd" during computation of t(:,:,:,nnew,:) below is merely
! a mechanism to ensure constancy preservation, at the expense of
! loosing conservation property.
!
! This is acceptable because t(:,:,:,n+1/2,:) fields will be used
! exclussively to compute the tracer fluxes during subsequent
! step3d_t operation, and the final values of t(i,j,k,n+1,itrc)
! alfer step3d_t will be computed in a flux-conservative manner.
! The overall time step will be both conservative and constancy
! preserving.

      indx=3-nstp

      if (FIRST_TIME_STEP) then       ! Coefficients for alternative-
        dtau=0.5*dt                   ! form LF-AM3 stepping algorithm,
        cf_stp=1.                     ! see Fig. 8, Eqs. (2.38)-(2.39);
        cf_bak=0.                     ! also Eq. (4.8) from SM2005;
      else                            ! Here "dtau" is the actual time
        dtau=dt*(1.-AM3_crv)          ! increment of predictor substep;
        cf_stp=0.5 +AM3_crv
        cf_bak=0.5 -AM3_crv
      endif                           ! Construct artificial slow-time
                                      ! continuity equation for pseudo-
      do k=1,nz                        ! compressible predictor substep,
        cff=0.5*dtau                  ! Eq. (4.7) from SM2005.
        do j=0,ny
          do i=0,nx
            FlxDiv=cff*pm(i,j)*pn(i,j)*( FlxU(i+1,j,k)-FlxU(i,j,k)
     &                                  +FlxV(i,j+1,k)-FlxV(i,j,k)
     &               +We(i,j,k)+Wi(i,j,k) -We(i,j,k-1)-Wi(i,j,k-1)
     &                                                           )
            Hz_bak(i,j,k)=Hz(i,j,k) +FlxDiv
            Hz_fwd(i,j,k)=Hz(i,j,k) -FlxDiv
          enddo
        enddo

# define FX UFx
# define FE VFe
                                             ! Advance tracer fields
        do itrc=1,nt                         ! starting with applying
                                             ! horizontal fluxes. Also
          call t_hadv_pre(FX,FE,k,itrc)
          call river_tfluxes(FX,FE,k,itrc)
                                             ! INITIALIZE CORRECTOR
          do j=1,ny                          ! STEP by pre-multiplying
            do i=1,nx
              t(i,j,k,nnew,itrc)=Hz_bak(i,j,k)*(
     &                                    cf_stp*t(i,j,k,nstp,itrc)
     &                                   +cf_bak*t(i,j,k,indx,itrc) )
     &                    -dtau*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                           +FE(i,j+1)-FE(i,j) )

              t(i,j,k,indx,itrc)=Hz(i,j,k)*t(i,j,k,nstp,itrc)
            enddo
          enddo                              ! tracer at "nstp" by Hz
        enddo !<-- itrc                      ! also at "nstp" before
                                             ! it is owerwritten by

# ifdef NHMG
#  include "compute_horiz_rhs_w_terms.h"
# endif

# undef FE                                /* ! barotropic mode. */
# undef FX

        call uv_cor_curv(ru,rv,k)
        call uv_hadv_pre(ru,rv,k)

      enddo !<-- k

#ifdef WEC
      call compute_wec_rhs_uv_terms(1,nx,1,ny,ru,rv,wrk1,wrk2)
#endif


! Compute dynamic bottom drag coefficient.  Note that "rd" must be
! aliased to a horizontal array beyond #4 to avoid overlap with one
! of the vertical arrays used below.

# define rd wrk1
# include "compute_rd_bott_drag.h"

! Continue computation of the tracer fields: 
!  - Compute its change due to vertical advection.  
!  - Do an implicit solve for vertical mixing and implicit vertical advection

! set the implicit advection and mixing  to zero on the boundaries
      om(:,0)  = 0
      om(:,nz) = 0
      Av(:,0)  = 0
      Av(:,nz) = 0

      do j=1,ny  !! Start of the big j-loop

        do i=1,nx
          DC(i,0)=dtau*pm(i,j)*pn(i,j)
        enddo
        do k=1,nz-1
          do i=1,nx
           om(i,k)=DC(i,0)*Wi(i,j,k)
           Av(i,k)= Akt(i,j,k)
          enddo
        enddo
        do k=1,nz
          do i=1,nx
           dz(i,k)= Hz_fwd(i,j,k)
          enddo
        enddo

        do itrc=1,nt

          call t_vadv_pre(FC,j,itrc)

          do k=1,nz
            do i=1,nx   !! t(..) is Hz*T
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &            -dtau*pm(i,j)*pn(i,j)*(FC(i,k)-FC(i,k-1))
            enddo
          enddo 

          call vmix_solve(t(:,:,:,nnew,itrc),dtau,dz,Av,om,j)

        enddo   !<-- itrc



!                      ! At this moment arrays "ru", "rv" contain all
! Momentum equations:  ! terms of finite-volume r.h.s. for 3D momentum
!--------- ----------  ! equation EXCEPT the implicit part of vertical
!                      ! advection, vertical viscosity, bottom drag,
!                      ! and wind forcing at surface.
!                      ! ru,rv, and rw are volume integrated and [m4/s2]

        ! add the explicit part of avection to the rhs of u and v
        call uv_vadv_pre(ru,rv,j)

!-----  Solve implicitly for vertical mixing and vertical advection of u -----
        do i=1,nx
          DC(i,0)=dtau*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        enddo
        do k=1,nz
          do i=1,nx
            u(i,j,k,nnew)=0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*(
     &                   cf_stp*u(i,j,k,nstp)+cf_bak*u(i,j,k,indx) )
     &                                        + DC(i,0)*ru(i,j,k)

            !! dz*u^n. After index rotation this will be u(nnew)
            u(i,j,k,indx)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nstp)
          enddo
        enddo
        do i=1,nx
           u(i,j,nz,nnew) = u(i,j,nz,nnew) + dtau*sustr(i,j)
        enddo

        do k=1,nz-1
          do i=1,nx
            om(i,k)= DC(i,0)*0.5*( Wi(i,j,k)+Wi(i-1,j,k) )
            Av(i,k)= 0.5*( Akv(i,j,k) + Akv(i-1,j,k) )
          enddo
        enddo
        do k=1,nz
          do i=1,nx
            dz(i,k)= 0.5*( Hz_fwd(i,j,k) + Hz_fwd(i-1,j,k) )
          enddo
        enddo
        do i=1,nx
          drag(i) = 0.5*(rd(i,j) + rd(i-1,j) )
        enddo

        call vmix_solve(u(:,:,:,nnew),dtau,dz,Av,om,j,drag)

!-----  Solve implicitly for vertical mixing and vertical advection of v -----
        do i=1,nx
          DC(i,0)=dtau*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
        enddo
        do k=1,nz
          do i=1,nx
            v(i,j,k,nnew)=0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*(
     &                  cf_stp*v(i,j,k,nstp)+cf_bak*v(i,j,k,indx) )
     &                                          + DC(i,0)*rv(i,j,k)

            v(i,j,k,indx)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nstp)
          enddo
        enddo
        do i=1,nx
           v(i,j,nz,nnew) = v(i,j,nz,nnew) + dtau*svstr(i,j)
        enddo

        do k=1,nz-1
          do i=1,nx
            om(i,k)= DC(i,0)*0.5*( Wi(i,j,k)+Wi(i,j-1,k) )
            Av(i,k)= 0.5*( Akv(i,j,k) + Akv(i,j-1,k) )
          enddo
        enddo
        do k=1,nz
          do i=1,nx
            dz(i,k)= 0.5*( Hz_fwd(i,j,k) + Hz_fwd(i,j-1,k) )
          enddo
        enddo
        do i=1,nx
          drag(i) = 0.5*(rd(i,j) + rd(i,j-1) )
        enddo

        call vmix_solve(v(:,:,:,nnew),dtau,dz,Av,om,j,drag)

!-----  Solve implicitly for vertical mixing and vertical advection of w -----
# ifdef NHMG
![
#  include "compute_vert_rhs_w_terms.h"

        do i=1,nx
          DC(i,0) =dtau*pm(i,j)*pn(i,j)
        enddo
        !! w(indx),w(nstp) in m/s,
        !! w(0) is always zero
        do k=1,nz-1
          do i=1,nx
            w(i,j,k,nnew)=0.5*(Hz_bak(i,j,k+1)+Hz_bak(i,j,k))*(
     &                  cf_stp*w(i,j,k,nstp)+cf_bak*w(i,j,k,indx) )
     &                                           + DC(i,0)*rw(i,j,k)

            w(i,j,k,indx)=0.5*(Hz(i,j,k+1)+Hz(i,j,k))*w(i,j,k,nstp)
          enddo
        enddo
        !! w(nz) is weighted with 0.5*Hz(nz)
        do i=1,nx
          w(i,j,nz,nnew)=0.5*(  Hz_bak(i,j,nz))*(
     &                cf_stp*w(i,j,nz,nstp)+cf_bak*w(i,j,nz,indx) )
     &                                         + DC(i,0)*rw(i,j,nz)
          w(i,j,nz,indx)=0.5*(  Hz(i,j,nz))*w(i,j,nz,nstp)
        enddo

        do k=0,nz
          do i=1,nx
            om(i,k)= DC(i,0)*Wi(i,j,k)
            Av(i,k)= Akv(i,j,k)
          enddo
        enddo
        do k=1,nz
          do i=1,nx
            dz(i,k)= Hz_fwd(i,j,k)
          enddo
        enddo

        call vmix_solve_w(w(:,:,:,nnew),dtau,dz,Av,om,j)

        !------- end computing of w(:,:,:,nnew) -----
!]
# endif /* NHMG */

      enddo     !<-- j
# undef rd
      !! velocities are now in m/s

      if (river_source) then
        do j=1,ny
          do i=1,nx
            if (abs(riv_uflx(i,j)).gt.1e-3) then
              riv_depth = 0.5*( z_w(i-1,j,nz)-z_w(i-1,j,0)
     &                     +z_w(i  ,j,nz)-z_w(i  ,j,0) )
              iriver = nint(riv_uflx(i,j)/10)
              riv_uvel = riv_vol(iriver)*(riv_uflx(i,j)-10*iriver)/
     &                   ( dn_u(i,j)*riv_depth)
              do k= 1,nz
                u(i,j,k,nnew) = riv_uvel
              enddo
            endif
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            if (abs(riv_vflx(i,j)).gt.1e-3) then
              riv_depth = 0.5*( z_w(i,j-1,nz)-z_w(i,j-1,0)
     &                     +z_w(i,j  ,nz)-z_w(i,j  ,0) )
              iriver = nint(riv_vflx(i,j)/10)
              riv_vvel = riv_vol(iriver)*(riv_vflx(i,j)-10*iriver)/
     &                   ( dm_v(i,j)*riv_depth)
              do k= 1,nz
                v(i,j,k,nnew) = riv_vvel
              enddo
            endif
          enddo
        enddo
      endif  ! <-- river_source


      call u3dbc(wrk1)
      call v3dbc(wrk1)

! WARNING: Preliminary time step for 3D momentum equitions is not
! complete after this moment: the computed fields u,v(i,j,k,nnew)
! have wrong vertical integrals, which will be corrected later
! after computation of barotropic mode.


! Set PHYSICAL lateral boundary conditions for tracer fields.

      do itrc=1,nt
        call t3dbc(itrc, wrk1)
# ifdef EXCHANGE
        ! Another opportunity to pack more mpi_exchanges
        call exchange_xxx(t(:,:,:,nnew,itrc))
# endif
      enddo

# ifdef NHMG
!======================================================================
! Non-hydro projection of momentum
!======================================================================

!     if (iic.lt.-10) then ! should never happen, so always project twice
      if (iic.ge.(ntstart+2)) then ! Use AB2 extrapolation for the NH prsgrd
!       if  (mynode==0) print *, 'AB2 extrapolate'
        do k=1,nz
          do j=1,ny
            do i=1,nx
              u(i,j,k,nnew) = u(i,j,k,nnew)
     &            + (1.5*nhdu(i,j,k,iprec2)-0.5*nhdu(i,j,k,iprec1))*dtau
     &            /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i-1,j,k))) * pn_u(i,j)
            enddo
          enddo

          do j=1,ny
            do i=1,nx
              v(i,j,k,nnew) = v(i,j,k,nnew)
     &            + (1.5*nhdv(i,j,k,iprec2)-0.5*nhdv(i,j,k,iprec1))*dtau
     &            /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i,j-1,k))) * pm_v(i,j)
            enddo
          enddo

          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = w(i,j,k,nnew)
     &            + (1.5*nhdw(i,j,k,iprec2)-0.5*nhdw(i,j,k,iprec1))*dtau
     &            * (pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo

      else  !! instead of using AB2 to extrapolate the nh pressure gradient from corrector steps, we project here
      if (west_bnd) then
        do k=1,nz
          do j=1,ny
            Hz_fwd(0,j,k)=Hz_fwd(1,j,k)
          enddo
        enddo
      endif
      if (east_bnd) then
        do k=1,nz
          do j=1,ny
            Hz_fwd(nx+1,j,k)=Hz_fwd(nx,j,k)
          enddo
        enddo
      endif
      if (south_bnd) then
        do k=1,nz
          do i=1,nx
            Hz_fwd(i,0,k)=Hz_fwd(i,1,k)
          enddo
        enddo
      endif
      if (north_bnd) then
        do k=1,nz
          do i=1,nx
            Hz_fwd(i,ny+1,k)=Hz_fwd(i,ny,k)
          enddo
        enddo
      endif

!    At the end of the tridiag solves, the velocities are in m/s
!    Multiply u,v,w to get back fluxes for use inside mg_solve
        do k=1,nz
          do j=1,ny
            do i=1,nx
              u(i,j,k,nnew) = u(i,j,k,nnew)
     &            *(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i-1,j,k)))/pn_u(i,j)
            enddo
          enddo
          do j=1,ny
            do i=1,nx
              v(i,j,k,nnew) = v(i,j,k,nnew)
     &            *(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i,j-1,k)))/pm_v(i,j)
            enddo
          enddo
          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = w(i,j,k,nnew)/(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo
#  ifdef EXCHANGE
        !! try to see if we could avoid this. NHMG only!
        call exchange_xxx(u(:,:,:,nnew),v(:,:,:,nnew))
#  endif

        ! Compute a 'barotropic' correction to w such that it matches
        ! the current ubar,vbar divergence (see Molemaker et al., 2018)
        if (surface_neumann) then
!         print *,'neumann'
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
     &        (nh_ubar(i+1,j)-nh_ubar(i,j)+nh_vbar(i,j+1)-nh_vbar(i,j))
            enddo
          enddo
          do k=1,nz
            do j=1,ny
              do i=1,nx
                w(i,j,k,nnew) = w(i,j,k,nnew) - nh_wcor(i,j)
     &                         * (z_w(i,j,k)-z_w(i,j,0))
     &                         / (z_w(i,j,nz)-z_w(i,j,0))
              enddo
            enddo
          enddo
        endif  ! surface_neumann !

        call nhmg_solve(nx,ny,nz,halo,padd_X,padd_E,
     &          u(:,:,:,nnew),v(:,:,:,nnew),w(:,:,:,nnew) )

      dmax = 0
      do k=1,nz
        do j=1,ny
          do i=1,nx
            div= (u(i+1,j,k,nnew) + mggrid(1)%du(k,j,i+1))
     &         - (u(i  ,j,k,nnew) + mggrid(1)%du(k,j,i  ))
     &         + (v(i,j+1,k,nnew) + mggrid(1)%dv(k,j+1,i))
     &         - (v(i,j  ,k,nnew) + mggrid(1)%dv(k,j  ,i))
     &         + (w(i,j,k  ,nnew) + mggrid(1)%dw(k+1,j,i))
     &         - (w(i,j,k-1,nnew) + mggrid(1)%dw(k  ,j,i))
            dmax = max(dmax,abs(div))
          enddo
        enddo
      enddo

        ! add the non-hydro correction (du, dv, dw are in flux form)
        ! add the non-hydrostatic pressure gradient to rufrc, rvfrc
        ! Can't do to nx+1 because we dont have Hz_fwd(nx+1)
        do k=1,nz
          do j=1,ny
            do i=1,nx
              u(i,j,k,nnew) = (u(i,j,k,nnew) + mggrid(1)%du(k,j,i))
     &           /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i-1,j,k)))* pn_u(i,j)
            enddo
          enddo

          do j=1,ny
            do i=1,nx
              v(i,j,k,nnew) = (v(i,j,k,nnew) + mggrid(1)%dv(k,j,i))
     &           /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i,j-1,k)))* pm_v(i,j)
            enddo
          enddo

          do j=1,ny
            do i=1,nx
              w(i,j,k,nnew) = (w(i,j,k,nnew) + mggrid(1)%dw(k+1,j,i))
     &             * (pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo

      endif  !! AB2 or projection

!     ! velocities are now in m/s

      call w3dbc(wrk1)
      call exchange_xxx(w(:,:,:,nnew))

      !! u and v are exchanged in set_depth (set_HUV1)
      !! for neatness, we could move this exchange there as well


# endif /* NHMG */


      end

