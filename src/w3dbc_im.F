#include "cppdefs.opt"
#if defined SOLVE3D && defined NHMG

      subroutine w3dbc_tile (grad)

      use param
      use boundary
      use grid
      use ocean_vars
      use scalars

      implicit none
      real grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cff, cext, dtfwd
      real, parameter :: eps=1.E-33

      if (PRED_STAGE) then    ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

# define dt illegal

# define OBC_M3SPECIFIEDW

# ifndef EW_PERIODIC

      if (WESTERN_EDGE) then
#  if defined OBC_WEST 
#   if defined OBC_M3ORLANSKI
        do k=1,N
          do j=1,ny+1
            grad(0,j)= w(0,j,k,nstp)- w(0,j-1,k,nstp)
            grad(1  ,j)= w(1  ,j,k,nstp)- w(1  ,j-1,k,nstp)
          enddo

          do j=1,ny

            cx=-dtfwd*u(1,j,k,nrhs)*pm(0,j)
            cy=0.5*dtfwd*(v(0,j,k,nrhs)+v(0,j+1,k,nrhs))
     &                                               *pn(0,j)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            w(0,j,k,nnew)=(1.-cx)*( w(0,j,k,nstp)
     &                                -max(cy,0.)*grad(0,j  )
     &                                -min(cy,0.)*grad(0,j+1) )
     &                             +cx*( w(1  ,j,k,nstp)
     &                                -max(cy,0.)*grad(1  ,j  )
     &                                -min(cy,0.)*grad(1  ,j+1) )

#    if defined M3_FRC_BRY
            w(0,j,k,nnew)= (1.-cext)*w(0,j,k,nnew)
     &                                +cext*w_west(j,k) 
#    endif
          enddo
        enddo
#   else  /* not OBC_M3ORLANSKI */
        ! not orlanski, can be either no gradient, or specified
        do k=1,N
          do j=1,ny
#    if defined OBC_M3SPECIFIED
            w(0,j,k,nnew)=w_west(j,k)
#    else
            w(0,j,k,nnew)=w(1,j,k,nnew)  ! no gradient
#    endif
          enddo
        enddo
#   endif /* OBC_M3ORLANSKI or SPECIFIED/NO-GRAD*/

#  else /* !OBC_WEST (closed)*/
        do k=1,N                          ! Wall: slip: gamma= 1
          do j=1,ny                  !    no-slip: gamma=-1
            w(0,j,k,nnew)=gamma2*w(1,j,k,nnew)
          enddo
        enddo
#  endif
      endif     ! <-- WESTERN_EDGE


      if (EASTERN_EDGE) then
#  if defined OBC_EAST 
#   if defined OBC_M3ORLANSKI
        do k=1,N
          do j=1,ny+1
           grad(nx  ,j)= w(nx  ,j,k,nstp)-w(nx  ,j-1,k,nstp)
           grad(nx+1,j)= w(nx+1,j,k,nstp)-w(nx+1,j-1,k,nstp)
          enddo
          do j=1,ny
            cx=    dtfwd* u(nx+1,j,k,nrhs)*pm(nx+1,j)
            cy=0.5*dtfwd*(v(nx+1,j,k,nrhs)+v(nx+1,j+1,k,nrhs))
     &                                               *pn(nx+1,j)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            w(nx+1,j,k,nnew)=(1.-cx)*( w(nx+1,j,k,nstp)
     &                          -max(cy,0.)*grad(nx+1,j  )
     &                          -min(cy,0.)*grad(nx+1,j+1)  )
     &                            +cx *( w(nx  ,j,k,nnew)
     &                          -max(cy,0.)*grad(nx  ,j  )
     &                          -min(cy,0.)*grad(nx  ,j+1)  )

#    if defined M3_FRC_BRY
            w(nx+1,j,k,nnew)=((1.-cext)*w(nx+1,j,k,nnew)
     &                             +cext*w_east(j,k)  )
#    endif
          enddo
        enddo
#   else  /* not ORLANSKI  (Spec or no grad)*/
        do k=1,N
          do j=1,ny
#    if defined OBC_M3SPECIFIED
            w(nx+1,j,k,nnew)=w_east(j,k)
#    else
            w(nx+1,j,k,nnew)=w(nx,j,k,nnew)
#    endif
          enddo
        enddo
#   endif /* orlanski,spec, or no-grad */

#  else /* not OBC_EAST (closed) */
        do k=1,N                          ! Wall: slip: gamma= 1
          do j=1,ny                  !    no-slip: gamma=-1
            w(nx+1,j,k,nnew)=gamma2*w(nx,j,k,nnew)
          enddo
        enddo
#  endif        /* OBC_EAST or not */
      endif    ! <-- EASTERN_EDGE

# endif        /* not EW_PERIODIC */



# ifndef NS_PERIODIC

      if (SOUTHERN_EDGE) then
#  if defined OBC_SOUTH
#   if defined OBC_M3ORLANSKI
        do k=1,N
          do i=1,nx+1
            grad(i,1  )= w(i  ,1  ,k,nstp)-w(i-1,1  ,k,nstp)
            grad(i,0)= w(i  ,0,k,nstp)-w(i-1,0,k,nstp)
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
            w(i,0,k,nnew)=(1.-cx)*( w(i,0,k,nstp)
     &                               -max(cy,0.)*grad(i  ,0)
     &                               -min(cy,0.)*grad(i+1,0) )
     &                             +cx*( w(i,1  ,k,nstp)
     &                               -max(cy,0.)*grad(i  ,1  )
     &                               -min(cy,0.)*grad(i+1,1  ) )
#    if defined W_FRC_BRY
            w(i,0,k,nnew)=( (1.-cext)*w(i,0,k,nnew)
     &                              +cext*w_south(i,k) )
#    endif
          enddo
        enddo
#   else /* not orlanski, but spec or no-grad */
        do k=1,N
          do i=1,nx
#    if defined OBC_WSPECIFIED
            w(i,0,k,nnew)=w_south(i,k)
#    else
            w(i,0,k,nnew)=w(i,1,k,nnew)
#    endif
          enddo
        enddo
#   endif /* orlanski, spec, or no-grad */
#  else /* not open, but closed */
        do k=1,N                          ! Wall: slip: gamma= 1
          do j=1,ny                  !    no-slip: gamma=-1
            w(i,0,k,nnew)=gamma2*w(i,1,k,nnew)
          enddo
        enddo
#  endif /* OBC_SOUTH or not */
      endif    ! <-- SOUTHERN_EDGE


      if (NORTHERN_EDGE) then
#  if defined OBC_NORTH
#   if defined OBC_M3ORLANSKI
        do k=1,N
          do i=1,nx+1
            grad(i,ny  )= w(i  ,ny  ,k,nstp)-w(i-1,ny  ,k,nstp)
            grad(i,ny+1)= w(i  ,ny+1,k,nstp)-w(i-1,ny+1,k,nstp)
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

            w(i,ny+1,k,nnew)=(1.-cx)*( w(i,ny+1,k,nstp)
     &                             -max(cy,0.)*grad(i  ,ny+1)
     &                             -min(cy,0.)*grad(i+1,ny+1) )
     &                             +cx*( w(i,ny  ,k,nnew)
     &                             -max(cy,0.)*grad(i  ,ny  )
     &                             -min(cy,0.)*grad(i+1,ny  ) )
#    if defined W_FRC_BRY 
            w(i,ny+1,k,nnew)= ((1.-cext)*w(i,ny+1,k,nnew)
     &                             +cext*w_north(i,k) )
#    endif
          enddo
        enddo
#   else  /* not orlanski, spec or non-grad instead */
        do k=1,N
          do i=1,nx
#    if defined OBC_M3SPECIFIED
            w(i,ny+1,k,nnew)=w_north(i,k)
#    else
            w(i,ny+1,k,nnew)=w(i,ny,k,nnew)  
#    endif
          enddo
        enddo
#   endif  /* orlanski, spec, or non-grad */
#  else /* not open, but closed */
        do k=1,N                          ! Wall: slip: gamma= 1
          do j=1,ny                  !    no-slip: gamma=-1
            w(i,ny+1,k,nnew)=gamma2*w(i,ny,k,nnew)
          enddo
        enddo
#  endif /* OBC_NORTH or not */
      endif    ! <-- NORTHERN_EDGE

# endif /* not NS_PERIODIC */

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
            w(0,0,k,nnew)=cff*(
     &              rmask(1,0)*w(1,0,k,nnew)
     &             +rmask(0,1)*w(0,1,k,nnew))
          enddo
        else
          do k=1,N
            w(0,0,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(0,0,k,nnew)=0.5*( w(1,0,k,nnew)
     &                                + w(0,1,k,nnew))
        enddo
#   endif
      endif

      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(nx,0)+rmask(nx+1,1)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            w(nx+1,0,k,nnew)=cff*(
     &              rmask(nx,0)*w(nx,0,k,nnew)
     &             +rmask(nx+1,1)*w(nx+1,1,k,nnew))
          enddo
        else
          do k=1,N
            w(nx+1,0,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(nx+1,0,k,nnew)=0.5*( w(nx,0,k,nnew)
     &                                + w(nx+1,1,k,nnew) )
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(1,ny+1)+rmask(0,ny)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            w(0,ny+1,k,nnew)=cff*(
     &              rmask(1,ny+1)*w(1,ny+1,k,nnew)
     &             +rmask(0,ny)*w(0,ny,k,nnew))
          enddo
        else
          do k=1,N
            w(0,ny+1,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(0,ny+1,k,nnew)=0.5*( w(1,ny+1,k,nnew)
     &                                + w(0,ny,k,nnew))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(nx,ny+1)+rmask(nx+1,ny)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            w(nx+1,ny+1,k,nnew)=cff*(
     &              rmask(nx,ny+1)*w(nx,ny+1,k,nnew)
     &             +rmask(nx+1,ny)*w(nx+1,ny,k,nnew))
          enddo
        else
          do k=1,N
            w(nx+1,ny+1,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(nx+1,ny+1,k,nnew)= 0.5*( w(nx,ny+1,k,nnew)
     &                                 + w(nx+1,ny,k,nnew))
        enddo
#   endif
      endif
#  endif
# endif
      end

#else 

      subroutine w3dbc_empty
      end

#endif /* SOLVE3D && NHMG */
