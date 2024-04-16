#include "cppdefs.h"
#ifdef SOLVE3D

# define UPSTREAM_TS
c--# define AKIMA
                                  ! Note that arrays A2d(1,1:4) within
      subroutine step3d_t(tile)   ! the list of agruments in call below
      implicit none               ! are repeated to use the same memory
      integer tile                ! for vertical and horizontal scratch
                                  ! arrays.
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"

      call step3d_t_tile( istr,iend,jstr,jend,
     &                    A2d(1,1), A2d(1,2), A2d(1,3), A2d(1,4),
     &                    A2d(1,1), A2d(1,2), A2d(1,3)         )
      end

      subroutine step3d_t_tile( istr,iend,jstr,jend, WC,FC,CF,DC,
     &                                                FX,FE,wrk1)
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY)     :: FX,FE,wrk1

      real cff
      real, parameter :: epsil=1.E-33
      integer itrc, iAkt

# include "grid.h"
# include "ocean3d.h"
# include "forces.h"
# include "mixing.h"
# include "climat.h"
# include "scalars.h"
# include "sources.h"

# include "compute_auxiliary_bounds.h"

! Compute and add in contribution due to horizontal advection fluxes.
! Note that at this moment the content of array t(:,:,:,nnew,itrc) is
! concentration of tracer "itrc" taken at time corresponding to "nstp"
! and already MULTIPLIED by Hz(:,:,:) at "nstp" in "pre_step3d_tile".

      do k=1,N
        do itrc=1,NT

# include "compute_horiz_tracer_fluxes.h"

          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) -dt*pm(i,j)
     &                              *pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                        +FE(i,j+1)-FE(i,j)
     &                                                         )
            enddo
          enddo                 !-->  discard FX,FE
        enddo
      enddo

      do j=jstr,jend                  ! Compute and add explicit part
        do itrc=1,NT                  ! of vertical advective fluxes.

# include "compute_vert_tracer_fluxes.h"

          do k=1,N
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) -dt*pm(i,j)*pn(i,j)
     &                                              *(FC(i,k)-FC(i,k-1))
            enddo
          enddo               !--> discard FC

! Add surface fluxes

          do i=istr,iend
            t(i,j,N,nnew,itrc)=t(i,j,N,nnew,itrc)+dt*stflx(i,j,itrc)
          enddo

# ifdef LMD_KPP
! Add the solar radiation flux in temperature equation. Also compute
! the nonlocal transport flux for unstable (convective) forcing
! conditions into matrix DC when using the Large et al. 1994 KPP
! scheme.

          if (itrc == itemp) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=srflx(i,j)*swr_frac(i,j,k)
#  ifdef LMD_NONLOCAL
     &                 -ghat(i,j,k)*(stflx(i,j,itemp)-srflx(i,j))
#  endif
                t(i,j,k+1,nnew,itemp)=t(i,j,k+1,nnew,itemp) -dt*cff
                t(i,j,k  ,nnew,itemp)=t(i,j,k  ,nnew,itemp) +dt*cff
              enddo
            enddo

c?? WARNING: the above implies that light (if any) reaching all the
c?? way to the bottom is entirely absorbed within the botom-most grid
c?? box, rather that reflected/scattered back to the water column. It
c?? is not clear, however, how to deal with this situation correctly
c??
c??         do i=istr,iend
c??           cff=srflx(i,j)*swr_frac(i,j,0)
c??#  ifdef LMD_NONLOCAL
c??  &               -ghat(i,j,0)*(stflx(i,j,itemp)-srflx(i,j))
c??#  endif
c??           t(i,j,1,nnew,itemp)=t(i,j,1,nnew,itemp) -dt*cff
c??         enddo

#  if defined LMD_NONLOCAL && defined SALINITY
          elseif (itrc == isalt) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=-dt*ghat(i,j,k)*stflx(i,j,isalt)
                t(i,j,k+1,nnew,isalt)=t(i,j,k+1,nnew,isalt) -cff
                t(i,j,k  ,nnew,isalt)=t(i,j,k  ,nnew,isalt) +cff
              enddo
            enddo
#  endif
          endif
# endif

! Perform an implicit time step for vertical diffusion
!
!   dq(k)     1     [         q(k+1)-q(k)             q(k)-q(k-1) ]
!  ------ = ----- * [ Akt(k)* ----------- - Akt(k-1)* ----------- ]
!    dt     Hz(k)   [            dz(k)                   dz(k-1)  ]
!
! where q(k) represents tracer field t(:,:,k,:,itrc). Doing so
! implies solution of a tri-diagonal system
!
!     -FC(k-1)*q_new(k-1) +[Hz(k)+FC(k-1)+FC(k)]*q_new(k)
!                       -FC(k)*q_new(k+1) = Hz(k)*q_old(k)
!
!                dt*Akt(k)
! where FC(k) = ----------- is normalized diffusivity coefficient
!                  dz(k)
!
! defined at W-points; q_new(k) is the new-time-step (unknown) tracer
! field; q_old(k) is old-time-step tracer (known).   As long as Akt(k)
! is nonnegative, the tri-diagonal matrix is diagonally dominant which
! guarantees stability of the Gaussian elimination algorithm.  Top and
! bottom boundary conditions are assumed to be no-flux, effectively
! Akt(N)=Akt(0)=0, hence FC(N)=FC(1)=0. This leads to equations for
! top and bottom grid boxes;
!
!   -FC(N-1)*q_new(N-1) +[Hz(N)+FC(N-1)]*q_new(N) = Hz(N)*q_old(N)
!
!          [Hz(1)+FC(1)]*q_new(1) -FC(1)*q_new(2) = Hz(1)*q_old(1)
!
! The FC(N)=FC(0)=0 boundary conditions does not mean that physical
! boundary conditions are no flux: the forcing fluxes have been
! applied explicitly above.  Instead, the no-flux condition should
! be interpreted as that the implicit step merely redistributes the
! tracer concentration throughout the water column. At this moment
! the content of array t(:,:,:,nnew,itrc) has meaning of Hz*tracer.
! After the implicit step it becomes just tracer.

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
            DC(i,0)=dt*pm(i,j)*pn(i,j) !<--horizontal metric

            FC(i,1)=2.*dt*Akt(i,j,1,iAkt)/(Hz(i,j,1)+Hz(i,j,2))
            WC(i,1)=DC(i,0)*Wi(i,j,1)

            cff=1./(Hz(i,j,1) +FC(i,1)+max(WC(i,1),0.))
            CF(i,1)=cff*(      FC(i,1)-min(WC(i,1),0.))
            DC(i,1)=cff*t(i,j,1,nnew,itrc)
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=2.*dt*Akt(i,j,k,iAkt)/(Hz(i,j,k)+Hz(i,j,k+1))
              WC(i,k)=Wi(i,j,k)*DC(i,0)

              cff=1./( Hz(i,j,k) +FC(i,k)+max(WC(i,k),0.)
     &                         +FC(i,k-1)-min(WC(i,k-1),0.)
     &              -CF(i,k-1)*(FC(i,k-1)+max(WC(i,k-1),0.))
     &                                                     )

              CF(i,k)=cff*(       FC(i,k)-min(WC(i,k),0.) )

              DC(i,k)=cff*( t(i,j,k,nnew,itrc) +DC(i,k-1)*(
     &                          FC(i,k-1)+max(WC(i,k-1),0.) ))
            enddo
          enddo          !--> discard DC(:,0)
          do i=istr,iend
            t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc) +DC(i,N-1)*(
     &                               FC(i,N-1)+max(WC(i,N-1),0.) )
     &                   )/( Hz(i,j,N) +FC(i,N-1)-min(WC(i,N-1),0.)
     &                      -CF(i,N-1)*(FC(i,N-1)+max(WC(i,N-1),0.))
     &                                                             )
# ifdef MASKING
     &                                                  *rmask(i,j)
# endif
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=(DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc))
# ifdef MASKING
     &                                                      *rmask(i,j)
# endif

            enddo
          enddo           !--> discard FC,WC,CF,DC

#ifdef CONST_TRACERS
          do k=1,N
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nstp,itrc)
            enddo
          enddo
#endif
        enddo      ! <-- itrc
      enddo      ! <-- j

! Set lateral boundary conditions; nudge toward tracer climatology;
! apply land-sea mask and exchange periodic boundary conditions.

      do itrc=1,NT
        call t3dbc_tile (istr,iend,jstr,jend, itrc, wrk1)

# if defined MOORE_BEC
      enddo                                    ! Multi-Group species
                                               ! biological Model (MGM)
      call ecosys_tile (istr,iend,jstr,jend)
      call mgm_avg_tile (istr,iend,jstr,jend)
      call ecosys_vst_tile (istr,iend,jstr,jend, nnew,dt)

      do itrc=1,NT
# elif defined BIOLOGY
      enddo                                    ! Insert interactive
      call biology_tile (istr,iend,jstr,jend)  ! part of the biological
      do itrc=1,NT                             ! model
# endif

# if defined MASKING || defined TCLIMATOLOGY
        do k=1,N
          do j=jstrR,jendR
            do i=istrR,iendR
#  ifdef TCLIMATOLOGY
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &           +dt*Tnudgcof(i,j,itrc)*( tclm(i,j,k,itrc)
     &                                 -t(i,j,k,nnew,itrc))
#  endif
#  ifdef MASKING
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
#  endif
            enddo
          enddo
        enddo
# endif
      enddo  ! <-- itrc

! Exchange computational margins and/or periodic boundaries:
!--------- ------------- -------- ------ -------- -----------
! This is done only here if the code is configured without explicit
! lateral diffusion. If there is one, then the exchange routine is
! called from the diffusive routine which immediately follows this
! step, therefore there is no need to call it here.

# if !defined TS_DIF2 && !defined TS_DIF4
#  ifdef EXCHANGE
      do itrc=1,NT
        call exchange_tile(istr,iend,jstr,jend,
     &        t(START_2D_ARRAY,1,nnew,itrc), N)
      enddo
#  endif
# endif
      end

      subroutine check_step_t_switches(ierr)

! The following routine keeps track of the status of CPP settings in
! "compute_horiz_tracer_fluxes.h" and "compute_horiz_vert_fluxes.h".
! This code must be placed here rather than in a separate file so it
! is exposed to the relevant CPP-settings.   It does not affect any
! model results, other than signature in global attribute "CPPS" in
! the output netCDF files.

      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2 ; ie=is+13
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<step3d_t3S.F>'

# include "track_advec_switches.h"

      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
#else
      subroutine step3d_t_empty
      end
#endif  /* SOLVE3D */
