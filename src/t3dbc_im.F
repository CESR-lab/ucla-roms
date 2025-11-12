#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine t3dbc_tile (itrc, grad)

! Set lateral boundary conditions for tracer field t(:,:,:,itrc)

      use param
      use tracers, only: t
      use boundary
      use grid
      use ocean_vars
      use scalars
      use obc_tune

      implicit none
      integer :: itrc, i,j,k
      real grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cff, cext, dtfwd
c-->     &                                           , dft,dfx,dfy
      real, parameter :: eps=1.E-33

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
          do j=1,ny+1
            grad(0,j)=( t(0,j  ,k,nstp,itrc)
     &                      -t(0,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                 *vmask(0,j)
#   endif
            grad(1  ,j)=( t(1  ,j  ,k,nstp,itrc)
     &                      -t(1  ,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                   *vmask(1,j)
#   endif
          enddo
          do j=1,ny
            cx=-dtfwd*u(1,j,k,nrhs)*pm(0,j)
            cy=0.5*dtfwd*(v(0,j,k,nrhs)+v(0,j+1,k,nrhs))
     &                                               *pn(0,j)
            if (cx > 0.) then ! outward flow
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(0,j,k,nnew,itrc)=(1.-cx)*( t(0,j,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(0,j  )
     &                                    -min(cy,0.)*grad(0,j+1)
     &                                                              )
     &                                  +cx*(  t(1,j,k,nstp,itrc)
     &                                      -max(cy,0.)*grad(1,j  )
     &                                      -min(cy,0.)*grad(1,j+1)
     &                                                              )

#   if defined T_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_west(j),1.0))
            t(0,j,k,nnew,itrc)=(1.-cext)*t(0,j,k,nnew,itrc)
     &                                    +cext*t_west(j,k,itrc)
#   endif
#   ifdef MASKING
            t(0,j,k,nnew,itrc)=t(0,j,k,nnew,itrc)
     &                                      *rmask(0,j)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do j=1,ny
#   if defined OBC_WEST && defined OBC_TSPECIFIED
            t(0,j,k,nnew,itrc)=t_west(j,k,itrc)
#   else
            t(0,j,k,nnew,itrc)=t(1,j,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(0,j)
#   endif
          enddo
        enddo
#  endif
      endif     ! <-- WESTERN_EDGE

      if (EASTERN_EDGE) then
#  if defined OBC_EAST && defined OBC_TORLANSKI

                                         !  Eastern edge radiation BC
        do k=1,N                         !  ======= ==== ========= ==
          do j=1,ny+1
           grad(nx  ,j)=( t(nx  ,j  ,k,nstp,itrc)
     &                     -t(nx  ,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                  *vmask(nx,j)
#   endif
           grad(nx+1,j)=( t(nx+1,j  ,k,nstp,itrc)
     &                     -t(nx+1,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                *vmask(nx+1,j)
#   endif
          enddo
          do j=1,ny

            cx=dtfwd*u(nx+1,j,k,nrhs)*pm(nx+1,j)
            cy=0.5*dtfwd*(v(nx+1,j,k,nrhs)+v(nx+1,j+1,k,nrhs))
     &                                               *pn(nx+1,j)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(nx+1,j,k,nnew,itrc)=(1.-cx)*( t(nx+1,j,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(nx+1,j  )
     &                                    -min(cy,0.)*grad(nx+1,j+1)
     &                                                              )
     &                                  +cx*(  t(nx,j,k,nnew,itrc)
     &                                      -max(cy,0.)*grad(nx,j  )
     &                                      -min(cy,0.)*grad(nx,j+1)
     &                                                              )

#   if defined T_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_east(j),1.0))
            t(nx+1,j,k,nnew,itrc)=(1.-cext)*t(nx+1,j,k,nnew,itrc)
     &                                         +cext*t_east(j,k,itrc)
#   endif
#   ifdef MASKING
            t(nx+1,j,k,nnew,itrc)=t(nx+1,j,k,nnew,itrc)
     &                                     *rmask(nx+1,j)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do j=1,ny
#   if defined OBC_EAST && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            t(nx+1,j,k,nnew,itrc)=t_east(j,k,itrc)
#    else
            t(nx+1,j,k,nnew,itrc)=tclm(nx+1,j,k,itrc)
#    endif
#   else
            t(nx+1,j,k,nnew,itrc)=t(nx,j,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(nx+1,j)
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
          do i=1,nx+1
            grad(i,1  )=( t(i  ,1  ,k,nstp,itrc)
     &                      -t(i-1,1  ,k,nstp,itrc))
#   ifdef MASKING
     &                                   *umask(i,1)
#   endif
            grad(i,0)=( t(i  ,0,k,nstp,itrc)
     &                      -t(i-1,0,k,nstp,itrc))
#   ifdef MASKING
     &                                *umask(i,0)
#   endif
          enddo
          do i=1,nx
            cx=-dtfwd*v(i,1,k,nrhs)*pn(i,0)
            cy=0.5*dtfwd*(u(i,0,k,nrhs)+u(i+1,0,k,nrhs))
     &                                               *pm(i,0)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(i,0,k,nnew,itrc)=(1.-cx)*( t(i,0,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(i  ,0)
     &                                    -min(cy,0.)*grad(i+1,0)
     &                                                              )
     &                                   +cx*(  t(i,1,k,nstp,itrc)
     &                                      -max(cy,0.)*grad(i  ,1)
     &                                      -min(cy,0.)*grad(i+1,1)
     &                                                              )

#   if defined T_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_south(i),1.0))
            t(i,0,k,nnew,itrc)=(1.-cext)*t(i,0,k,nnew,itrc)
     &                                        +cext*t_south(i,k,itrc)
#   endif
#   ifdef MASKING
            t(i,0,k,nnew,itrc)=t(i,0,k,nnew,itrc)
     &                                     *rmask(i,0)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=1,nx
#   if defined OBC_SOUTH && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            t(i,0,k,nnew,itrc)=t_south(i,k,itrc)
#    else
            t(i,0,k,nnew,itrc)=tclm(i,0,k,itrc)
#    endif
#   else
            t(i,0,k,nnew,itrc)=t(i,1,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(i,0)
#   endif
          enddo
        enddo
#  endif
      endif    ! <-- SOUTHERN_EDGE


      if (NORTHERN_EDGE) then
#  if defined OBC_NORTH && defined OBC_TORLANSKI
        do k=1,N
          do i=1,nx+1
            grad(i,ny  )=( t(i  ,ny  ,k,nstp,itrc)
     &                      -t(i-1,ny  ,k,nstp,itrc))
#   ifdef MASKING
     &                                   *umask(i,ny)
#   endif
            grad(i,ny+1)=( t(i  ,ny+1,k,nstp,itrc)
     &                      -t(i-1,ny+1,k,nstp,itrc))
#   ifdef MASKING
     &                                 *umask(i,ny+1)
#   endif
          enddo
          do i=1,nx
            cx=dtfwd*v(i,ny+1,k,nrhs)*pn(i,ny+1)
            cy=0.5*dtfwd*(u(i,ny+1,k,nrhs)+u(i+1,ny+1,k,nrhs))
     &                                               *pm(i,ny+1)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            t(i,ny+1,k,nnew,itrc)=(1.-cx)*( t(i,ny+1,k,nstp,itrc)
     &                                    -max(cy,0.)*grad(i  ,ny+1)
     &                                    -min(cy,0.)*grad(i+1,ny+1)
     &                                                              )
     &                                  +cx*(  t(i,ny,k,nnew,itrc)
     &                                      -max(cy,0.)*grad(i  ,ny)
     &                                      -min(cy,0.)*grad(i+1,ny)
     &                                                              )

#   if defined T_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_north(i),1.0))
            t(i,ny+1,k,nnew,itrc)=(1.-cext)*t(i,ny+1,k,nnew,itrc)
     &                                        +cext*t_north(i,k,itrc)
#   endif
#   ifdef MASKING
            t(i,ny+1,k,nnew,itrc)=t(i,ny+1,k,nnew,itrc)
     &                                     *rmask(i,ny+1)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=1,nx
#   if defined OBC_NORTH && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            t(i,ny+1,k,nnew,itrc)=t_north(i,k,itrc)
#    else
            t(i,ny+1,k,nnew,itrc)=tclm(i,ny+1,k,itrc)
#    endif
#   else
            t(i,ny+1,k,nnew,itrc)=t(i,ny,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(i,ny+1)
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
        cff=rmask(1,0)+rmask(0,1)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(0,0,k,nnew,itrc)=cff*(
     &              rmask(1,0)*t(1,0,k,nnew,itrc)
     &             +rmask(0,1)*t(0,1,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(0,0,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(0,0,k,nnew,itrc)=0.5*( t(1,0,k,nnew,
     &                           itrc)+t(0,1,k,nnew,itrc))
        enddo
#   endif
      endif

      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(nx,0)+rmask(nx+1,1)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(nx+1,0,k,nnew,itrc)=cff*(
     &              rmask(nx,0)*t(nx,0,k,nnew,itrc)
     &             +rmask(nx+1,1)*t(nx+1,1,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(nx+1,0,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(nx+1,0,k,nnew,itrc)=0.5*(t(nx,0,k,nnew,
     &                           itrc)+t(nx+1,1,k,nnew,itrc))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(1,ny+1)+rmask(0,ny)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(0,ny+1,k,nnew,itrc)=cff*(
     &              rmask(1,ny+1)*t(1,ny+1,k,nnew,itrc)
     &             +rmask(0,ny)*t(0,ny,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(0,ny+1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(0,ny+1,k,nnew,itrc)=0.5*( t(1,ny+1,k,nnew,
     &                           itrc)+t(0,ny,k,nnew,itrc))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(nx,ny+1)+rmask(nx+1,ny)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            t(nx+1,ny+1,k,nnew,itrc)=cff*(
     &              rmask(nx,ny+1)*t(nx,ny+1,k,nnew,itrc)
     &             +rmask(nx+1,ny)*t(nx+1,ny,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(nx+1,ny+1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(nx+1,ny+1,k,nnew,itrc)=0.5*( t(nx,ny+1,k,nnew,
     &                           itrc)+t(nx+1,ny,k,nnew,itrc))
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
