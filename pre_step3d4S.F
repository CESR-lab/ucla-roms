#include "cppdefs.h"
#ifdef SOLVE3D
                                    ! Note that arrays A2d(1,1:4) in
      subroutine pre_step3d (tile)  ! the call just below are repeated
      implicit none                 ! within the list of agruments to
      integer tile                  ! use the same memory for vertical
# include "param.h"              /* ! and horizontal scratch arrays. */
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call pre_step3d_tile (istr,iend,jstr,jend,   A3d(1,1),A3d(1,2),
     &         A3d(1,3),A3d(1,4), A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                             A2d(1,5),A2d(1,6))
      end

      subroutine pre_step3d_tile (istr,iend,jstr,jend, ru,rv, Hz_bak,
     &                                         Hz_fwd,   WC,FC,CF,DC,
     &                                    UFx,UFe,VFx,VFe, wrk1,wrk2)
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv,Hz_bak,
     &                                                      Hz_fwd
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) ::  WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
      real dtau, cf_stp,cf_bak,  cff, FlxDiv
      real, parameter :: AM3_crv=1./6., epsil=1.E-33
     &                 , delta=0.1666666666666667 !! delta=0.125
# ifdef UPSTREAM_UV
     &                 , gamma=0.3333333333333333 !! gamma=0.25
# endif
      integer indx, itrc, iAkt
# ifdef PSOURCE
     &                  , is
# endif
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# include "sources.h"

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
! computation of "Hz_fwd" mimics the time step for tracers themselves).
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
        do itrc=1,NT                         ! starting with applying
                                             ! horizontal fluxes. Also
# include "compute_horiz_tracer_fluxes.h"
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
# undef FE                                /* ! barotropic mode. */
# undef FX

# include "compute_horiz_rhs_uv_terms.h"

      enddo !<-- k

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

      do j=jstr,jend
        do itrc=1,NT

# include "compute_vert_tracer_fluxes.h"

          do k=1,N
            do i=istr,iend
             t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &            -dtau*pm(i,j)*pn(i,j)*(FC(i,k)-FC(i,k-1))

# ifdef CONST_TRACERS
              t(i,j,k,nnew,itrc)=t(i,j,k,nstp,itrc)
# endif
            enddo
          enddo            !--> discard FC

# ifdef BIOLOGY
#  ifdef SALINITY
          iAkt=min(itrc,isalt)
#  else
          iAkt=min(itrc,itemp)
#  endif
# else
          iAkt=itrc
# endif
          do i=istr,iend
            FC(i,1)=2.*dtau*Akt(i,j,1,iAkt)/( Hz_fwd(i,j,2)
     &                                       +Hz_fwd(i,j,1))
            DC(i,0)=dtau*pm(i,j)*pn(i,j)
            WC(i,1)=DC(i,0)*Wi(i,j,1)

            cff=1./(Hz_fwd(i,j,1) +FC(i,1)+max(WC(i,1),0.))
            CF(i,1)=cff*(          FC(i,1)-min(WC(i,1),0.))
            DC(i,1)=cff*t(i,j,1,nnew,itrc)
          enddo

          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=2.*dtau*Akt(i,j,k,iAkt)/( Hz_fwd(i,j,k+1)
     &                                         +Hz_fwd(i,j,k) )
              WC(i,k)=DC(i,0)*Wi(i,j,k)

              cff=1./( Hz_fwd(i,j,k) +FC(i,k)+max(WC(i,k),0.)
     &                             +FC(i,k-1)-min(WC(i,k-1),0.)
     &                  -CF(i,k-1)*(FC(i,k-1)+max(WC(i,k-1),0.))
     &                                                          )
              CF(i,k)=cff*(FC(i,k)-min(WC(i,k),0.))

              DC(i,k)=cff*( t(i,j,k,nnew,itrc) +DC(i,k-1)*(
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

# include "compute_vert_rhs_uv_terms.h"

        do i=istrU,iend
          DC(i,0)=dtau*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        enddo
        do k=1,N
          do i=istrU,iend
            DC(i,k)=0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*(
     &                  cf_stp*u(i,j,k,nstp)+cf_bak*u(i,j,k,indx) )
     &                                           +DC(i,0)*ru(i,j,k)

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
# ifdef WIND_AT_RHO_POINTS
     &               +0.5*dtau*(sustr(i-1,j)+sustr(i,j))
# else
     &                                  +dtau*sustr(i,j)
# endif
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
# ifdef PRED_COUPLED_MODE
          rufrc(i,j)=ru(i,j,1) +dm_u(i,j)*dn_u(i,j)*(
#  ifdef WIND_AT_RHO_POINTS
     &                   0.5*(sustr(i-1,j)+sustr(i,j))
#  else
     &                                      sustr(i,j)
#  endif
     &        -0.5*(rd(i-1,j)+rd(i,j))*u(i,j,1,nnew) )
# endif
        enddo
        do k=2,N,+1          !--> backsubstitution
          do i=istrU,iend
            u(i,j,k,nnew)=DC(i,k) +CF(i,k-1)*u(i,j,k-1,nnew)
# ifdef PRED_COUPLED_MODE
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k)
# endif
          enddo
        enddo

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
# ifdef WIND_AT_RHO_POINTS
     &                 +0.5*dtau*(svstr(i,j-1)+svstr(i,j))
# else
     &                                    +dtau*svstr(i,j)
# endif
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

# ifdef PRED_COUPLED_MODE
            rvfrc(i,j)=rv(i,j,1) +dm_v(i,j)*dn_v(i,j)*(
#  ifdef WIND_AT_RHO_POINTS
     &                     0.5*(svstr(i,j-1)+svstr(i,j))
#  else
     &                                        svstr(i,j)
#  endif
     &          -0.5*(rd(i,j-1)+rd(i,j))*v(i,j,1,nnew) )
# endif
          enddo
          do k=2,N,+1          !--> backsubstitution
            do i=istr,iend
              v(i,j,k,nnew)=DC(i,k) +CF(i,k-1)*v(i,j,k-1,nnew)
# ifdef PRED_COUPLED_MODE
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k)
# endif
            enddo
          enddo
        endif  !<-- j>=jstrV
      enddo     !<-- j
# undef rd

      call u3dbc_tile (istr,iend,jstr,jend, wrk1)
      call v3dbc_tile (istr,iend,jstr,jend, wrk1)

# ifdef PSOURCE
      do is=1,Nsrc           ! Apply point mass sources associated
        i=Isrc(is)           ! with river runoff simulations.
        j=Jsrc(is)
        if ( istrR <= i .and. i <= iendR .and.
     &       jstrR <= j .and. j <= jendR ) then
          if (Dsrc(is) == 0) then
            do k=1,N
              u(i,j,k,nnew)=2.*Qsrc(is,k)/( dn_u(i,j)*(
     &                       z_w(i-1,j,k)-z_w(i-1,j,k-1)
     &                      +z_w(i  ,j,k)-z_w(i  ,j,k-1)
     &                                                ))
            enddo
          else
            do k=1,N
              v(i,j,k,nnew)=2.*Qsrc(is,k)/( dm_v(i,j)*(
     &                       z_w(i,j-1,k)-z_w(i,j-1,k-1)
     &                      +z_w(i,j  ,k)-z_w(i,j  ,k-1)
     &                                                ))
            enddo
          endif
        endif
      enddo
# endif

! WARNING: Preliminary time step for 3D momentum equitions is not
! complete after this moment: the computed fields u,v(i,j,k,nnew)
! have wrong vertical integrals, which will be corrected later
! after computation of barotropic mode.


! Set PHYSICAL lateral boundary conditions for tracer fields.

      do itrc=1,NT
        call t3dbc_tile (istr,iend,jstr,jend, itrc, wrk1)
# ifdef EXCHANGE
        call exchange_tile (istr,iend,jstr,jend,
     &                      t(START_2D_ARRAY,1,nnew,itrc), N)
# endif
      enddo
      end


      subroutine check_pre_step_switches (ierr)

! This routine keeps track of the status of local CPP-settings in
! "pre_step3d4S.F".  It must be placed here rather than a separate
! file in order to be exposed to the relevant CPP-settings. It does
! not affect any model results, other than "CPPS" signature saved
! as a global attribute in output netCDF files.

      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2 ; ie=is+15
      if (ie>max_opt_size) goto 99
      cpps(is:ie)='<pre_step3d4S.F>'

# include "track_advec_switches.h"

      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      return
      end

#else
      subroutine pre_step3d_empty
      end
#endif  /* SOLVE3D */
