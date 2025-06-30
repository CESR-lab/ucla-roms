#include "cppdefs.opt"
#ifdef SOLVE3D
                                    ! Note that arrays A2d(1,1:4) in
      subroutine pre_step3d (tile)  ! the call just below are repeated
      use param
      use hidden_mpi_vars
      use private_scratch
      implicit none                 ! within the list of agruments to
      integer tile                  ! use the same memory for vertical
                                    ! and horizontal scratch arrays.

# include "compute_tile_bounds.h"
      call pre_step3d_tile (istr,iend,jstr,jend,   A3d(1,1),A3d(1,2),
     &         A3d(1,3),A3d(1,4), A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                             A2d(1,5),A2d(1,6)
# ifdef NHMG
     &                                            ,A3d(1,5)
# endif
     &          )
      end

      subroutine pre_step3d_tile (istr,iend,jstr,jend, ru,rv, Hz_bak,
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
      use diagnostics                     ! needed for compute_horiz_uv_terms.h below
      use coupling
      use grid
      use mixing
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv,Hz_bak,
     &                                                      Hz_fwd
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) ::  WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
#ifdef NHMG
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: rw
      real  Flxw,Uflxw,Vflxw
      real  div,dmax
      real  ub,ut
!     integer :: ierr !,imx,jmx,kmx
#endif

      real dtau, cf_stp,cf_bak,  cff, FlxDiv
      real, parameter :: AM3_crv=1./6., epsil=1.E-33
     &                 , delta=0.1666666666666667 !! delta=0.125
# ifdef UPSTREAM_UV
     &                 , gamma=0.3333333333333333 !! gamma=0.25
# endif
      integer indx, itrc, iAkt

! Taken from rhs3d35S.F for WEC, used in compute_horiz_rhs_uv_terms.h
! but needs to be declared here in declarations section not in
! the compute_horiz_rhs_uv_terms.h file else wont compile.
# ifdef WEC
      real cff1
# endif

# include "compute_auxiliary_bounds.h"

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
      do k=1,N                        ! compressible predictor substep,
        cff=0.5*dtau                  ! Eq. (4.7) from SM2005.
        do j=jstrV-1,jend
          do i=istrU-1,iend
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
# include "compute_horiz_tracer_fluxes.h"
!! river included in above include file
                                             ! INITIALIZE CORRECTOR
          do j=jstr,jend                     ! STEP by pre-multiplying
            do i=istr,iend
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

# include "compute_horiz_rhs_uv_terms.h"

      enddo !<-- k

#ifdef WEC
      call compute_wec_rhs_uv_terms(istr,iend,jstr,jend,ru,rv,wrk1,wrk2)
#endif


! Compute dynamic bottom drag coefficient.  Note that "rd" must be
! aliased to a horizontal array beyond #4 to avoid overlap with one
! of the vertical arrays used below.

# define rd wrk1
# include "compute_rd_bott_drag.h"

! ...continue computation of the auxiliary tracer field: compute its
! change due to vertical advection.  Vertical advective fluxes require
! interpolation of tracer values to the verical grid-box interfaces
! (W-points). This can be is done by either using parabolic spline
! interpolation or, more simple local cubic polynomial [with or without
! monotonicity constraint; linear interpolation is considered obsolete,
! but the code is retained for reference].

      do j=jstr,jend  !! Start of the giant j-loop

        do itrc=1,nt
# include "compute_vert_tracer_fluxes.h"

          do k=1,n
            do i=istr,iend   !! t(..) is Hz*T
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &            -dtau*pm(i,j)*pn(i,j)*(FC(i,k)-FC(i,k-1))
# ifdef PIPE_SOURCEXXX /*only add this during the corrector step */
              if (pipe_idx(i,j)>0) then
                pidx = pipe_idx(i,j)
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &            +dtau*pm(i,j)*pn(i,j)*
     &              pipe_flx(i,j)*pipe_prf(pidx,k)*pipe_trc(pidx,itrc)
              endif
# endif

# ifdef CONST_TRACERS
              t(i,j,k,nnew,itrc)=t(i,j,k,nstp,itrc)
# endif
            enddo
          enddo            !--> discard FC

          iAkt=min(itrc,iTandS)

          !! Hz*T(k) - dt*2*Ak(k  )*(T(k+1)-T(k  ))/(Hz(k+1)+Hz(k  ))
          !!         + dt*2*Ak(k-1)*(T(k  )-T(k-1))/(Hz(k  )+Hz(k-1)) =
          !!                                                  Hz*T(k)
          !! See numerical recipes p40
          !! b(1) = Hz(1) + dt*2*Ak(1)/(Hz(2)+Hz(1))
          !! c(1) =       - dt*2*Ak(1)/(Hz(2)+Hz(1))
          !! U(k) --> DC(k)
          do i=istr,iend
            FC(i,1)=2.*dtau*Akt(i,j,1,iAkt)/( Hz_fwd(i,j,2) !!  dt*2*Ak(1)/(Hz(2)+Hz(1))
     &                                       +Hz_fwd(i,j,1))
            DC(i,0)=dtau*pm(i,j)*pn(i,j)
            WC(i,1)=DC(i,0)*Wi(i,j,1)

            cff=1./(Hz_fwd(i,j,1) +FC(i,1)+max(WC(i,1),0.))  !! 1./bet
            CF(i,1)=cff*(          FC(i,1)-min(WC(i,1),0.))  !! c(1)/bet
            DC(i,1)=cff*t(i,j,1,nnew,itrc)                 !! u(1) = r(1)/bet
          enddo

          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=2.*dtau*Akt(i,j,k,iAkt)/( Hz_fwd(i,j,k+1)
     &                                         +Hz_fwd(i,j,k) )
              WC(i,k)=DC(i,0)*Wi(i,j,k)

              cff=1./( Hz_fwd(i,j,k) +FC(i,k)+max(WC(i,k),0.)  !!  1./((b(k) -a(k)*c(k-1))/bet)
     &                             +FC(i,k-1)-min(WC(i,k-1),0.)
     &                  -CF(i,k-1)*(FC(i,k-1)+max(WC(i,k-1),0.))
     &                                                          )
              CF(i,k)=cff*(FC(i,k)-min(WC(i,k),0.))             !! c(k)/bet

              DC(i,k)=cff*( t(i,j,k,nnew,itrc) +DC(i,k-1)*(     !!  u(k)=(r(k)-u(j-1)*a(j))/bet
     &                          FC(i,k-1)+max(WC(i,k-1),0.) ))
            enddo
          enddo  !--> discard DC(:,0)

          do i=istr,iend
            t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc) +DC(i,N-1)*(
     &                                FC(i,N-1)+max(WC(i,N-1),0.) )
     &               )/( Hz_fwd(i,j,N) +FC(i,N-1)-min(WC(i,N-1),0.)
     &                      -CF(i,N-1)*(FC(i,N-1)+max(WC(i,N-1),0.))
     &                                                            )
          enddo

          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
            enddo
          enddo
        enddo   !<-- itrc  !--> discard DC,CF,FC




!                      ! At this moment arrays "ru", "rv" contain all
! Momentum equations:  ! terms of finite-volume r.h.s. for 3D momentum
!--------- ----------  ! equation EXCEPT the implicit part of vertical
!                      ! advection, vertical viscosity, bottom drag,
!                      ! and wind forcing at surface.
!                      ! ru,rv, and rw are volume integrated and [m4/s2]

# include "compute_vert_rhs_uv_terms.h"

# ifdef WEC
#  if defined WAVE_FRICTION && !defined BODY_FRICTION
        do i=istrU,iend
          ru(i,j,1) =ru(i,j,1) + dm_u(i,j)*dn_u(i,j)*frc2dx(i,j)
        enddo
        do i=istr,iend
          rv(i,j,1) =rv(i,j,1) + dm_v(i,j)*dn_v(i,j)*frc2de(i,j)
        enddo
#  endif
# endif   /* WEC */

# ifdef NHMG

#  include "compute_vert_rhs_w_terms.h"

        do i=istr,iend
          DC(i,0) =dtau*pm(i,j)*pn(i,j)
        enddo
        !! w(indx),w(nstp) in m/s,
        !! w(0) is always zero
        do k=1,N-1
          do i=istr,iend
            DC(i,k)=0.5*(Hz_bak(i,j,k+1)+Hz_bak(i,j,k))*(
     &                  cf_stp*w(i,j,k,nstp)+cf_bak*w(i,j,k,indx) )
     &                                           + DC(i,0)*rw(i,j,k)

            w(i,j,k,indx)=0.5*(Hz(i,j,k+1)+Hz(i,j,k))*w(i,j,k,nstp)
          enddo
        enddo
        k = N
        !! here is the special volume weighting for w(N)
        !! Still, write out the volume integrated (dz*w) vertical mixing eq.
        !! Just to check :)
        do i=istr,iend
          DC(i,k)=0.5*(   Hz_bak(i,j,k))*(
     &                cf_stp*w(i,j,k,nstp)+cf_bak*w(i,j,k,indx) )
     &                                         + DC(i,0)*rw(i,j,k)

          w(i,j,k,indx)=0.5*(  Hz(i,j,k))*w(i,j,k,nstp)
        enddo

        !! start of tri-diag solve
        !   FC(i,1)=0.5*dtau*(Akt(i,j,k+1)+Akt(i,j,k)/( Hz_fwd(i,j,k)
        do i=istr,iend
          FC(i,N-1)= 0.5*dtau*(-Akv(i,j,N)+Akv(i,j,N-1)) ! see notes, exception for k=N
     &                      /Hz_fwd(i,j,N)

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N)+Wi(i,j,N-1))

          cff = 0.5*Hz_fwd(i,j,N) +FC(i,N-1)-min(WC(i,N-1),0.)

          CF(i,N-1)= ( FC(i,N-1)+max(WC(i,N-1),0.) )/cff ! gam(n-1)=c(n)/bet

          DC(i,N)=  DC(i,N)/cff                          ! u(1) = r(1)/bet
        enddo

        do k=N-1,2,-1      !--> forward elimination for w
          do i=istr,iend
            FC(i,k-1)= 0.5*dtau*(Akv(i,j,k)+Akv(i,j,k-1))!! A(k-1)
     &                   /Hz_fwd(i,j,k)

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k)+Wi(i,j,k-1))

            cff=1./( 0.5*(Hz_fwd(i,j,k)+Hz_fwd(i,j,k-1))
     &                            +FC(i,k-1)-min(WC(i,k-1),0.)
     &                              +FC(i,k)+max(WC(i,k),0.)
     &                     -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                       )
            CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*(DC(i,k)+DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)))
          enddo
        enddo
        !! Use the fact that w(0) is always 0
        do i=istr,iend
          FC(i,0)= 0.5*dtau*(Akv(i,j,1)+Akv(i,j,0))!! A(0)
     &                   /Hz_fwd(i,j,1)

          WC(i,0)= DC(i,0)*0.5*(Wi(i,j,1)+Wi(i,j,0))

          w(i,j,1,nnew)=( DC(i,1)  +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                         )/( 0.5*(Hz_fwd(i,j,2)+Hz_fwd(i,j,1))
     &                                       +FC(i,0)-min(WC(i,0),0.)
     &                                       +FC(i,1)+max(WC(i,1),0.)
     &                              -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                               )
        enddo
        do k=2,N,+1          !--> backsubstitution for w
          do i=istr,iend
            w(i,j,k,nnew)=DC(i,k) +CF(i,k-1)*w(i,j,k-1,nnew)
          enddo
        enddo
        !------- end computing of w(:,:,:,nnew) -----
# endif /* NHMG */

        do i=istrU,iend
          DC(i,0)=dtau*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        enddo
        do k=1,N
          do i=istrU,iend
            !! DC is dz*u
            DC(i,k)=0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*(
     &                  cf_stp*u(i,j,k,nstp)+cf_bak*u(i,j,k,indx) )
     &                                           +DC(i,0)*ru(i,j,k)

            !! also dz*u
            u(i,j,k,indx)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nstp)
          enddo
        enddo

        do i=istrU,iend
          FC(i,N-1)= 2.*dtau*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                 /( Hz_fwd(i,j,N  )+Hz_fwd(i-1,j,N  )
     &                   +Hz_fwd(i,j,N-1)+Hz_fwd(i-1,j,N-1))

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i-1,j,N-1))

          cff=1./( 0.5*(Hz_fwd(i,j,N)+Hz_fwd(i-1,j,N))
     &                    +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*(  FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*( DC(i,N)
     &                                  +dtau*sustr(i,j)
#    if defined WEC && defined SURFACE_BREAK
     &                            + dtau*brk2dx(i,j)
#    endif
     &                                                 )
        enddo
        do k=N-1,2,-1      !--> forward elimination
          do i=istrU,iend
            FC(i,k-1)= 2.*dtau*(Akv(i,j,k-1)+Akv(i-1,j,k-1))
     &                   /( Hz_fwd(i,j,k  )+Hz_fwd(i-1,j,k  )
     &                     +Hz_fwd(i,j,k-1)+Hz_fwd(i-1,j,k-1))

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i-1,j,k-1))

            cff=1./( 0.5*(Hz_fwd(i,j,k)+Hz_fwd(i-1,j,k))
     &                            +FC(i,k-1)-min(WC(i,k-1),0.)
     &                              +FC(i,k)+max(WC(i,k),0.)
     &                     -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                       )
            CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*(DC(i,k)+DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)))
          enddo
        enddo
        do i=istrU,iend
          u(i,j,1,nnew)=( DC(i,1)  +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                         )/( 0.5*(Hz_fwd(i,j,1)+Hz_fwd(i-1,j,1))
# ifdef IMPLCT_NO_SLIP_BTTM_BC
     &                                   +0.5*dtau*(rd(i,j)+rd(i-1,j))
# endif
     &                                       +FC(i,1)+max(WC(i,1),0.)
     &                              -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                               )
        enddo
        do k=2,N,+1          !--> backsubstitution
          do i=istrU,iend
            u(i,j,k,nnew)=DC(i,k) +CF(i,k-1)*u(i,j,k-1,nnew)
          enddo
        enddo
        !------- end computing of u(:,:,:,nnew) -----

        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dtau*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          enddo
          do k=1,N
            do i=istr,iend
              DC(i,k)=0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*(
     &                    cf_stp*v(i,j,k,nstp)+cf_bak*v(i,j,k,indx) )
     &                                            + DC(i,0)*rv(i,j,k)

              v(i,j,k,indx)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nstp)
            enddo
          enddo
          do i=istr,iend
            FC(i,N-1)= 2.*dtau*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                   /( Hz_fwd(i,j,N  )+Hz_fwd(i,j-1,N  )
     &                     +Hz_fwd(i,j,N-1)+Hz_fwd(i,j-1,N-1))

            WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i,j-1,N-1))

            cff=1./( 0.5*(Hz_fwd(i,j,N)+Hz_fwd(i,j-1,N))
     &                      +FC(i,N-1)-min(WC(i,N-1),0.) )

            CF(i,N-1)=cff*(  FC(i,N-1)+max(WC(i,N-1),0.) )

            DC(i,N)=cff*( DC(i,N)
     &                                    +dtau*svstr(i,j)
#    if defined WEC && defined SURFACE_BREAK
     &                            + dtau*brk2de(i,j)
#    endif
     &                                                   )
          enddo
          do k=N-1,2,-1      !--> forward elimination
            do i=istr,iend
              FC(i,k-1)= 2.*dtau*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                     /( Hz_fwd(i,j,k  )+Hz_fwd(i,j-1,k  )
     &                       +Hz_fwd(i,j,k-1)+Hz_fwd(i,j-1,k-1))

              WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i,j-1,k-1))

              cff=1./( 0.5*(Hz_fwd(i,j,k)+Hz_fwd(i,j-1,k))
     &                             +FC(i,k-1)-min(WC(i,k-1),0.)
     &                               +FC(i,k)+max(WC(i,k),0.)
     &                      -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                        )
              CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

              DC(i,k)=cff*(DC(i,k)+DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)))
            enddo
          enddo
          do i=istr,iend
            v(i,j,1,nnew)=( DC(i,1)+DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                         )/( 0.5*(Hz_fwd(i,j,1)+Hz_fwd(i,j-1,1))
# ifdef IMPLCT_NO_SLIP_BTTM_BC
     &                                   +0.5*dtau*(rd(i,j)+rd(i,j-1))
# endif
     &                                       +FC(i,1)+max(WC(i,1),0.)
     &                              -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                               )
          enddo
          do k=2,N,+1          !--> backsubstitution
            do i=istr,iend
              v(i,j,k,nnew)=DC(i,k) +CF(i,k-1)*v(i,j,k-1,nnew)
            enddo
          enddo
        endif  !<-- j>=jstrV
      enddo     !<-- j
# undef rd
      !! velocities are now in m/s

      if (river_source) then
        do j=jstr,jend
          do i=istrU,iend
            if (abs(riv_uflx(i,j)).gt.1e-3) then
              riv_depth = 0.5*( z_w(i-1,j,N)-z_w(i-1,j,0)
     &                     +z_w(i  ,j,N)-z_w(i  ,j,0) )
              iriver = nint(riv_uflx(i,j)/10)
              riv_uvel = riv_vol(iriver)*(riv_uflx(i,j)-10*iriver)/
     &                   ( dn_u(i,j)*riv_depth)
              do k= 1,N
                u(i,j,k,nnew) = riv_uvel
              enddo
            endif
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            if (abs(riv_vflx(i,j)).gt.1e-3) then
              riv_depth = 0.5*( z_w(i,j-1,N)-z_w(i,j-1,0)
     &                     +z_w(i,j  ,N)-z_w(i,j  ,0) )
              iriver = nint(riv_vflx(i,j)/10)
              riv_vvel = riv_vol(iriver)*(riv_vflx(i,j)-10*iriver)/
     &                   ( dm_v(i,j)*riv_depth)
              do k= 1,N
                v(i,j,k,nnew) = riv_vvel
              enddo
            endif
          enddo
        enddo
      endif  ! <-- river_source

!     if (split_bc) then
!       call u3dbc(u_slow)
!       call u3dbc(u_fast)
!       call v3dbc(v_slow)
!       call v3dbc(v_fast)

!     else
!     endif

      call u3dbc_tile (istr,iend,jstr,jend, wrk1)
      call v3dbc_tile (istr,iend,jstr,jend, wrk1)

! WARNING: Preliminary time step for 3D momentum equitions is not
! complete after this moment: the computed fields u,v(i,j,k,nnew)
! have wrong vertical integrals, which will be corrected later
! after computation of barotropic mode.


! Set PHYSICAL lateral boundary conditions for tracer fields.

      do itrc=1,NT
        call t3dbc_tile (istr,iend,jstr,jend, itrc, wrk1)
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
        do k=1,N
          do j=jstr,jend
            do i=istrU,iend
              u(i,j,k,nnew) = u(i,j,k,nnew)
     &            + (1.5*nhdu(i,j,k,iprec2)-0.5*nhdu(i,j,k,iprec1))*dtau
     &            /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i-1,j,k))) * pn_u(i,j)
            enddo
          enddo

          do j=jstrV,jend
            do i=istr,iend
              v(i,j,k,nnew) = v(i,j,k,nnew)
     &            + (1.5*nhdv(i,j,k,iprec2)-0.5*nhdv(i,j,k,iprec1))*dtau
     &            /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i,j-1,k))) * pm_v(i,j)
            enddo
          enddo

          do j=jstr,jend
            do i=istr,iend
              w(i,j,k,nnew) = w(i,j,k,nnew)
     &            + (1.5*nhdw(i,j,k,iprec2)-0.5*nhdw(i,j,k,iprec1))*dtau
     &            * (pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo

      else  !! instead of using AB2 to extrapolate the nh pressure gradient from corrector steps, we project here
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            Hz_fwd(istr-1,j,k)=Hz_fwd(istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            Hz_fwd(iend+1,j,k)=Hz_fwd(iend,j,k)
          enddo
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            Hz_fwd(i,jstr-1,k)=Hz_fwd(i,jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            Hz_fwd(i,jend+1,k)=Hz_fwd(i,jend,k)
          enddo
        enddo
      endif
#  endif

!    At the end of the tridiag solves, the velocities are in m/s
!    Multiply u,v,w to get back fluxes for use inside mg_solve
        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              u(i,j,k,nnew) = u(i,j,k,nnew)
     &            *(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i-1,j,k)))/pn_u(i,j)
            enddo
          enddo
          do j=jstr,jend
            do i=istr,iend
              v(i,j,k,nnew) = v(i,j,k,nnew)
     &            *(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i,j-1,k)))/pm_v(i,j)
            enddo
          enddo
          do j=jstr,jend
            do i=istr,iend
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
          do k=1,N
            do j=jstr,jend
              do i=istr,iend+1
                nh_ubar(i,j) = nh_ubar(i,j) + u(i,j,k,nnew)
              enddo
            enddo
            do j=jstr,jend+1
              do i=istr,iend
                nh_vbar(i,j) = nh_vbar(i,j) + v(i,j,k,nnew)
              enddo
            enddo
          enddo

          do j=jstr,jend
            do i=istr,iend
              nh_wcor(i,j) = w(i,j,N,nnew) +
     &        (nh_ubar(i+1,j)-nh_ubar(i,j)+nh_vbar(i,j+1)-nh_vbar(i,j))
            enddo
          enddo
          do k=1,N
            do j=jstr,jend
              do i=istr,iend
                w(i,j,k,nnew) = w(i,j,k,nnew) - nh_wcor(i,j)
     &                         * (z_w(i,j,k)-z_w(i,j,0))
     &                         / (z_w(i,j,N)-z_w(i,j,0))
              enddo
            enddo
          enddo
        endif  ! surface_neumann !

        call nhmg_solve(Lm,Mm,N,halo,padd_X,padd_E,
     &          u(:,:,:,nnew),v(:,:,:,nnew),w(:,:,:,nnew) )

      dmax = 0
      do k=1,N
        do j=jstr,jend
          do i=istr,iend
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
        ! Can't do to iend+1 because we dont have Hz_fwd(iend+1)
        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              u(i,j,k,nnew) = (u(i,j,k,nnew) + mggrid(1)%du(k,j,i))
     &           /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i-1,j,k)))* pn_u(i,j)
            enddo
          enddo

          do j=jstr,jend
            do i=istr,iend
              v(i,j,k,nnew) = (v(i,j,k,nnew) + mggrid(1)%dv(k,j,i))
     &           /(0.5*(Hz_fwd(i,j,k) + Hz_fwd(i,j-1,k)))* pm_v(i,j)
            enddo
          enddo

          do j=jstr,jend
            do i=istr,iend
              w(i,j,k,nnew) = (w(i,j,k,nnew) + mggrid(1)%dw(k+1,j,i))
     &             * (pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo

      endif  !! AB2 or projection

!     ! velocities are now in m/s

!     call u3dbc_tile(istr,iend,jstr,jend, wrk1)
!     call v3dbc_tile(istr,iend,jstr,jend, wrk1)
!     call w3dbc_tile(istr,iend,jstr,jend, wrk1)

      !! u and v are exchanged in set_depth (set_HUV1)
      !! for neatness, we could move this exchange there as well
#  ifdef EXCHANGE
      call exchange_xxx(w(:,:,:,nnew))
#  endif

# endif /* NHMG */

# ifdef OBC_CHECK
      print *, 'nnew=', nnew
      print *, '   u(istr  ,10,10,nnew)=', u(istr,10,10,nnew)
      print *, '   u(iend+1,10,10,nnew)=', u(iend+1,10,10,nnew)
# endif

      end

!------------------------------------------------------------------------------------------
      subroutine check_pre_step_switches (ierr)

! This routine keeps track of the status of local CPP-settings in
! "pre_step3d4S.F".  It must be placed here rather than a separate
! file in order to be exposed to the relevant CPP-settings. It does
! not affect any model results, other than "CPPS" signature saved
! as a global attribute in output netCDF files.

      use param
      use strings
      implicit none
      integer ierr, is,ie, lenstr

      ie=lenstr(cpps)
      is=ie+2 ; ie=is+15
      if (ie>max_opt_size) goto 99
      cpps(is:ie)='<pre_step3d4S.F>'

# include "track_advec_switches.h"

      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      return
      end

#else
      subroutine pre_step3d_empty
      end
#endif  /* SOLVE3D */
