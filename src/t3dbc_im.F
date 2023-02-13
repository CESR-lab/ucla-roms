#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine t3dbc_tile (istr,iend,jstr,jend, itrc, grad)

! Set lateral boundary conditions for tracer field t(:,:,:,itrc)

      use param
      use tracers, only: t
      use boundary
      use grid
      use ocean_vars
      use scalars
      use sponge_tune

      implicit none
      integer istr,iend,jstr,jend, itrc, i,j,k
      real grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cff, cext, dtfwd
c-->     &                                           , dft,dfx,dfy
      real, parameter :: eps=1.E-33

# include "compute_auxiliary_bounds.h"

      if (PRED_STAGE) then    ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

# define dt illegal

# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
#  if defined OBC_WEST && defined OBC_TORLANSKI
        do k=1,nz
          do j=jstr,jend+1
            grad(istr-1,j)=( t(istr-1,j  ,k,nstp,itrc)
     &                      -t(istr-1,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                 *vmask(istr-1,j)
#   endif
            grad(istr  ,j)=( t(istr  ,j  ,k,nstp,itrc)
     &                      -t(istr  ,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                   *vmask(istr,j)
#   endif
          enddo
          do j=jstr,jend
            cx=-dtfwd*u(istr,j,k,nrhs)*pm(istr-1,j)
            cy=0.5*dtfwd*(v(istr-1,j,k,nrhs)+v(istr-1,j+1,k,nrhs))
     &                                               *pn(istr-1,j)
            if (cx > 0.) then ! outward flow
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(istr-1,j,k,nnew,itrc)=(1.-cx)*( t(istr-1,j,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(istr-1,j  )
     &                                    -min(cy,0.)*grad(istr-1,j+1)
     &                                                              )
     &                                  +cx*(  t(istr,j,k,nstp,itrc)
     &                                      -max(cy,0.)*grad(istr,j  )
     &                                      -min(cy,0.)*grad(istr,j+1)
     &                                                              )

#   if defined T_FRC_BRY
#     if defined SPONGE_TUNE
            if (ub_tune) cext = max(cext,min(ub_west(j),1.0))
#     endif
            t(istr-1,j,k,nnew,itrc)=(1.-cext)*t(istr-1,j,k,nnew,itrc)
     &                                    +cext*t_west(j,k,itrc)
#   endif
#   ifdef MASKING
            t(istr-1,j,k,nnew,itrc)=t(istr-1,j,k,nnew,itrc)
     &                                      *rmask(istr-1,j)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do j=jstr,jend
#   if defined OBC_WEST && defined OBC_TSPECIFIED
            t(istr-1,j,k,nnew,itrc)=t_west(j,k,itrc)
#   else
            t(istr-1,j,k,nnew,itrc)=t(istr,j,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(istr-1,j)
#   endif
          enddo
        enddo
#  endif
      endif     ! <-- WESTERN_EDGE

      if (EASTERN_EDGE) then
#  if defined OBC_EAST && defined OBC_TORLANSKI

                                         !  Eastern edge radiation BC
        do k=1,N                         !  ======= ==== ========= ==
          do j=jstr,jend+1
           grad(iend  ,j)=( t(iend  ,j  ,k,nstp,itrc)
     &                     -t(iend  ,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                  *vmask(iend,j)
#   endif
           grad(iend+1,j)=( t(iend+1,j  ,k,nstp,itrc)
     &                     -t(iend+1,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                *vmask(iend+1,j)
#   endif
          enddo
          do j=jstr,jend

            cx=dtfwd*u(iend+1,j,k,nrhs)*pm(iend+1,j)
            cy=0.5*dtfwd*(v(iend+1,j,k,nrhs)+v(iend+1,j+1,k,nrhs))
     &                                               *pn(iend+1,j)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(iend+1,j,k,nnew,itrc)=(1.-cx)*( t(iend+1,j,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(iend+1,j  )
     &                                    -min(cy,0.)*grad(iend+1,j+1)
     &                                                              )
     &                                  +cx*(  t(iend,j,k,nnew,itrc)
     &                                      -max(cy,0.)*grad(iend,j  )
     &                                      -min(cy,0.)*grad(iend,j+1)
     &                                                              )

#   if defined T_FRC_BRY
#     if defined SPONGE_TUNE
            if (ub_tune) cext = max(cext,min(ub_east(j),1.0))
#     endif
            t(iend+1,j,k,nnew,itrc)=(1.-cext)*t(iend+1,j,k,nnew,itrc)
     &                                         +cext*t_east(j,k,itrc)
#   endif
#   ifdef MASKING
            t(iend+1,j,k,nnew,itrc)=t(iend+1,j,k,nnew,itrc)
     &                                     *rmask(iend+1,j)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do j=jstr,jend
#   if defined OBC_EAST && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            t(iend+1,j,k,nnew,itrc)=t_east(j,k,itrc)
#    else
            t(iend+1,j,k,nnew,itrc)=tclm(iend+1,j,k,itrc)
#    endif
#   else
            t(iend+1,j,k,nnew,itrc)=t(iend,j,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(iend+1,j)
#   endif
          enddo
        enddo
#  endif
      endif    ! <-- EASTERN_EDGE
# endif        /* !EW_PERIODIC */







# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  if defined OBC_SOUTH && defined OBC_TORLANSKI
        do k=1,N
          do i=istr,iend+1
            grad(i,jstr  )=( t(i  ,jstr  ,k,nstp,itrc)
     &                      -t(i-1,jstr  ,k,nstp,itrc))
#   ifdef MASKING
     &                                   *umask(i,jstr)
#   endif
            grad(i,jstr-1)=( t(i  ,jstr-1,k,nstp,itrc)
     &                      -t(i-1,jstr-1,k,nstp,itrc))
#   ifdef MASKING
     &                                *umask(i,jstr-1)
#   endif
          enddo
          do i=istr,iend
            cx=-dtfwd*v(i,jstr,k,nrhs)*pn(i,jstr-1)
            cy=0.5*dtfwd*(u(i,jstr-1,k,nrhs)+u(i+1,jstr-1,k,nrhs))
     &                                               *pm(i,jstr-1)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(i,jstr-1,k,nnew,itrc)=(1.-cx)*( t(i,jstr-1,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(i  ,jstr-1)
     &                                    -min(cy,0.)*grad(i+1,jstr-1)
     &                                                              )
     &                                   +cx*(  t(i,jstr,k,nstp,itrc)
     &                                      -max(cy,0.)*grad(i  ,jstr)
     &                                      -min(cy,0.)*grad(i+1,jstr)
     &                                                              )

#   if defined T_FRC_BRY
#     if defined SPONGE_TUNE
            if (ub_tune) cext = max(cext,min(ub_south(i),1.0))
#     endif
            t(i,jstr-1,k,nnew,itrc)=(1.-cext)*t(i,jstr-1,k,nnew,itrc)
     &                                        +cext*t_south(i,k,itrc)
#   endif
#   ifdef MASKING
            t(i,jstr-1,k,nnew,itrc)=t(i,jstr-1,k,nnew,itrc)
     &                                     *rmask(i,jstr-1)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=istr,iend
#   if defined OBC_SOUTH && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            t(i,jstr-1,k,nnew,itrc)=t_south(i,k,itrc)
#    else
            t(i,jstr-1,k,nnew,itrc)=tclm(i,jstr-1,k,itrc)
#    endif
#   else
            t(i,jstr-1,k,nnew,itrc)=t(i,jstr,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(i,jstr-1)
#   endif
          enddo
        enddo
#  endif
      endif    ! <-- SOUTHERN_EDGE


      if (NORTHERN_EDGE) then
#  if defined OBC_NORTH && defined OBC_TORLANSKI
        do k=1,N
          do i=istr,iend+1
            grad(i,jend  )=( t(i  ,jend  ,k,nstp,itrc)
     &                      -t(i-1,jend  ,k,nstp,itrc))
#   ifdef MASKING
     &                                   *umask(i,jend)
#   endif
            grad(i,jend+1)=( t(i  ,jend+1,k,nstp,itrc)
     &                      -t(i-1,jend+1,k,nstp,itrc))
#   ifdef MASKING
     &                                 *umask(i,jend+1)
#   endif
          enddo
          do i=istr,iend
            cx=dtfwd*v(i,jend+1,k,nrhs)*pn(i,jend+1)
            cy=0.5*dtfwd*(u(i,jend+1,k,nrhs)+u(i+1,jend+1,k,nrhs))
     &                                               *pm(i,jend+1)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(i,jend+1,k,nnew,itrc)=(1.-cx)*( t(i,jend+1,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(i  ,jend+1)
     &                                    -min(cy,0.)*grad(i+1,jend+1)
     &                                                              )
     &                                  +cx*(  t(i,jend,k,nnew,itrc)
     &                                      -max(cy,0.)*grad(i  ,jend)
     &                                      -min(cy,0.)*grad(i+1,jend)
     &                                                              )

#   if defined T_FRC_BRY
#     if defined SPONGE_TUNE
            if (ub_tune) cext = max(cext,min(ub_north(i),1.0))
#     endif
            t(i,jend+1,k,nnew,itrc)=(1.-cext)*t(i,jend+1,k,nnew,itrc)
     &                                        +cext*t_north(i,k,itrc)
#   endif
#   ifdef MASKING
            t(i,jend+1,k,nnew,itrc)=t(i,jend+1,k,nnew,itrc)
     &                                     *rmask(i,jend+1)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=istr,iend
#   if defined OBC_NORTH && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            t(i,jend+1,k,nnew,itrc)=t_north(i,k,itrc)
#    else
            t(i,jend+1,k,nnew,itrc)=tclm(i,jend+1,k,itrc)
#    endif
#   else
            t(i,jend+1,k,nnew,itrc)=t(i,jend,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(i,jend+1)
#   endif
          enddo
        enddo
#  endif
      endif    ! <-- NORTHERN_EDGE
# endif /* ! NS_PERIODIC */

! Corner points between adjacent boundaries. Note that because boundary
! conditions for tracers are of gradient nature -- either Neumann, or
! radiation, the corner points have to be set in all cases, even if the
! adjacent boundaries are closed.  This differs from setting boundaries
! for velocities, where there is anisotropy associated with normal and
! tangential components.   In the case when one/or both points near the
! corner is/are masked, the use of masked points is avoided.

# ifndef EW_PERIODIC
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(istr,jstr-1)+rmask(istr-1,jstr)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(istr-1,jstr-1,k,nnew,itrc)=cff*(
     &              rmask(istr,jstr-1)*t(istr,jstr-1,k,nnew,itrc)
     &             +rmask(istr-1,jstr)*t(istr-1,jstr,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(istr-1,jstr-1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(istr-1,jstr-1,k,nnew,itrc)=0.5*( t(istr,jstr-1,k,nnew,
     &                           itrc)+t(istr-1,jstr,k,nnew,itrc))
        enddo
#   endif
      endif

      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(iend,jstr-1)+rmask(iend+1,jstr)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(iend+1,jstr-1,k,nnew,itrc)=cff*(
     &              rmask(iend,jstr-1)*t(iend,jstr-1,k,nnew,itrc)
     &             +rmask(iend+1,jstr)*t(iend+1,jstr,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(iend+1,jstr-1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(iend+1,jstr-1,k,nnew,itrc)=0.5*(t(iend,jstr-1,k,nnew,
     &                           itrc)+t(iend+1,jstr,k,nnew,itrc))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(istr,jend+1)+rmask(istr-1,jend)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(istr-1,jend+1,k,nnew,itrc)=cff*(
     &              rmask(istr,jend+1)*t(istr,jend+1,k,nnew,itrc)
     &             +rmask(istr-1,jend)*t(istr-1,jend,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(istr-1,jend+1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(istr-1,jend+1,k,nnew,itrc)=0.5*( t(istr,jend+1,k,nnew,
     &                           itrc)+t(istr-1,jend,k,nnew,itrc))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(iend,jend+1)+rmask(iend+1,jend)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(iend+1,jend+1,k,nnew,itrc)=cff*(
     &              rmask(iend,jend+1)*t(iend,jend+1,k,nnew,itrc)
     &             +rmask(iend+1,jend)*t(iend+1,jend,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(iend+1,jend+1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(iend+1,jend+1,k,nnew,itrc)=0.5*( t(iend,jend+1,k,nnew,
     &                           itrc)+t(iend+1,jend,k,nnew,itrc))
        enddo
#   endif
      endif
#  endif
# endif
      end
#else
      subroutine t3dbc_empty
      end
#endif /* SOLVE3D */
