#include "cppdefs.h"
#if defined SOLVE3D && defined NHMG

      subroutine w3dbc_tile (istr,iend,jstr,jend,grad)
      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cff, cext, dtfwd
c-->     &                                           , dft,dfx,dfy
      real, parameter :: eps=1.E-33
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"

# include "compute_auxiliary_bounds.h"

      if (PRED_STAGE) then    ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

# define dt illegal

# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
#  if defined OBC_WEST 
#   if defined OBC_M3ORLANSKI
        do k=1,N
          do j=jstr,jend+1
            grad(istr-1,j)=( w(istr-1,j  ,k,nstp)
     &                      -w(istr-1,j-1,k,nstp))
#    ifdef MASKING
     &                            *vmask(istr-1,j)
#    endif
            grad(istr  ,j)=( w(istr  ,j  ,k,nstp)
     &                      -w(istr  ,j-1,k,nstp))
#    ifdef MASKING
     &                              *vmask(istr,j)
#    endif
          enddo
          do j=jstr,jend

            cx=-dtfwd*u(istr,j,k,nrhs)*pm(istr-1,j)
            cy=0.5*dtfwd*(v(istr-1,j,k,nrhs)+v(istr-1,j+1,k,nrhs))
     &                                               *pn(istr-1,j)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            w(istr-1,j,k,nnew)=(1.-cx)*( w(istr-1,j,k,nstp)
     &                                -max(cy,0.)*grad(istr-1,j  )
     &                                -min(cy,0.)*grad(istr-1,j+1)
     &                                                            )
     &                             +cx*( w(istr  ,j,k,nstp)
     &                                -max(cy,0.)*grad(istr,j  )
     &                                -min(cy,0.)*grad(istr,j+1)
     &                                                            )


#    if defined M3_FRC_BRY
            w(istr-1,j,k,nnew)=(1.-cext)*w(istr-1,j,k,nnew)
     &                                    +cext*w_west(j,k)
#    endif
#    ifdef MASKING
            w(istr-1,j,k,nnew)=w(istr-1,j,k,nnew)*rmask(istr-1,j)
#    endif
          enddo
        enddo
#   else  /* !OBC_M3ORLANSKI */
        do k=1,N
          do j=jstr,jend
#    if defined OBC_M3SPECIFIED
            w(istr-1,j,k,nnew)=w_west(j,k)
#    else
            w(istr-1,j,k,nnew)=w(istr,j,k,nnew)  ! no gradient
#    endif
#    ifdef MASKING
     &                                   *rmask(istr-1,j)
#    endif
          enddo
        enddo
#   endif /* OBC_M3ORLANSKI */

#  else /* !OBC_WEST (closed)*/
        do k=1,N                          ! Wall: slip: gamma= 1
          do j=jstr,jend                  !    no-slip: gamma=-1
            w(istr-1,j,k,nnew)=gamma2*w(istr-1,j,k,nnew)
          enddo
        enddo
#  endif
      endif     ! <-- WESTERN_EDGE



      if (EASTERN_EDGE) then
#  if defined OBC_EAST 
#   if defined OBC_M3ORLANSKI
        do k=1,N
          do j=jstr,jend+1
           grad(iend  ,j)=( w(iend  ,j  ,k,nstp)
     &                     -w(iend  ,j-1,k,nstp))
#   ifdef MASKING
     &                                *vmask(iend  ,j)
#   endif
           grad(iend+1,j)=( w(iend+1,j  ,k,nstp)
     &                     -w(iend+1,j-1,k,nstp))
#   ifdef MASKING
     &                                *vmask(iend+1,j)
#   endif
          enddo
          do j=jstr,jend

            cx=    dtfwd* u(iend+1,j,k,nrhs)*pm(iend+1,j)
            cy=0.5*dtfwd*(v(iend+1,j,k,nrhs)+v(iend+1,j+1,k,nrhs))
     &                                               *pn(iend+1,j)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            w(iend+1,j,k,nnew)=(1.-cx)*( w(iend+1,j,k,nstp)
     &                             -max(cy,0.)*grad(iend+1,j  )
     &                             -min(cy,0.)*grad(iend+1,j+1)
     &                                                         )
     &                             +cx*( w(iend  ,j,k,nnew)
     &                               -max(cy,0.)*grad(iend,j  )
     &                               -min(cy,0.)*grad(iend,j+1)
     &                                                         )

#   if defined M3_FRC_BRY
            w(iend+1,j,k,nnew)=(1.-cext)*w(iend+1,j,k,nnew)
     &                                  +cext*w_east(j,k)
#   endif
#   ifdef MASKING
            w(iend+1,j,k,nnew)=w(iend+1,j,k,nnew)*rmask(iend+1,j)
#   endif
          enddo
        enddo
#  else  /* not ORLANSKI */
        do k=1,N
          do j=jstr,jend
#   if defined OBC_M3SPECIFIED
            w(iend+1,j,k,nnew)=w_east(j,k)
#   else
            w(iend+1,j,k,nnew)=w(iend,j,k,nnew)
#   endif
#   ifdef MASKING
     &                               *rmask(iend+1,j)
#   endif
          enddo
        enddo
#  endif

#  else /* !OBC_EAST (closed) */
        do k=1,N                          ! Wall: slip: gamma= 1
          do j=jstr,jend                  !    no-slip: gamma=-1
            w(iend+1,j,k,nnew)=gamma2*w(iend+1,j,k,nnew)
          enddo
        enddo
#  endif        /* OBC_EAST */
      endif    ! <-- EASTERN_EDGE
# endif        /* !EW_PERIODIC */





# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  if defined OBC_SOUTH && defined OBC_TORLANSKI
        do k=1,N
          do i=istr,iend+1
            grad(i,jstr  )=( w(i  ,jstr  ,k,nstp)
     &                      -w(i-1,jstr  ,k,nstp))
#   ifdef MASKING
     &                                   *umask(i,jstr)
#   endif
            grad(i,jstr-1)=( w(i  ,jstr-1,k,nstp)
     &                      -w(i-1,jstr-1,k,nstp))
#   ifdef MASKING
     &                                *umask(i,jstr-1)
#   endif
          enddo
          do i=istr,iend
#   ifdef LEGACY_ORLANSKI
            dft=w(i,jstr,k,nstp)-w(i,jstr  ,k,nnew)
            dfx=w(i,jstr,k,nnew)-w(i,jstr+1,k,nnew)

#    ifdef OBC_RAD_NORMAL
            cy=0. ; cff=max(dfx*dfx, eps)
#    else
            if (dft*(grad(i,jstr)+grad(i+1,jstr)) > 0.) then
              dfy=grad(i,jstr)
            else
              dfy=grad(i+1,jstr)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
#     ifdef OBC_RAD_NPO
            cy=0.
#     else
            cy=min(cff,max(dft*dfy,-cff))
#     endif
#    endif

            cx=dft*dfx
            if (cx < 0.) then     !--> cancel, if inflow
              cx=0. ; cy=0.
#    if defined T_FRC_BRY || defined TNUDGING
              cext=dtfwd*ubind * pn(i,jstr-1)
            else
              cext=0.
#    endif
            endif

            w(i,jstr-1,k,nnew)=( cff*w(i,jstr-1,k,nstp)
     &                                   +cx*w(i,jstr,k,nnew)
     &                                -max(cy,0.)*grad(i  ,jstr-1)
     &                                -min(cy,0.)*grad(i+1,jstr-1)
     &                                                 )/(cff+cx)

#   else
            cx=-dtfwd*v(i,jstr,k,nrhs)*pn(i,jstr-1)
            cy=0.5*dtfwd*(u(i,jstr-1,k,nrhs)+u(i+1,jstr-1,k,nrhs))
     &                                               *pm(i,jstr-1)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            w(i,jstr-1,k,nnew)=(1.-cx)*( w(i,jstr-1,k,nstp)
     &                                    -max(cy,0.)*grad(i  ,jstr-1)
     &                                    -min(cy,0.)*grad(i+1,jstr-1)
     &                                                              )
     &                                   +cx*(  w(i,jstr,k,nstp)
     &                                      -max(cy,0.)*grad(i  ,jstr)
     &                                      -min(cy,0.)*grad(i+1,jstr)
     &                                                              )

#   endif /* LEGACY_ORLANSKI */

#   if defined T_FRC_BRY || defined TNUDGING
            w(i,jstr-1,k,nnew)=(1.-cext)*w(i,jstr-1,k,nnew)
#    ifdef T_FRC_BRY
     &                                    +cext*w_south(i,k)
#    else
     &                                    +cext*wclm(i,jstr-1,k)
#    endif
#   endif
#   ifdef MASKING
            w(i,jstr-1,k,nnew)=w(i,jstr-1,k,nnew)
     &                                    *rmask(i,jstr-1)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=istr,iend
#   if defined OBC_SOUTH && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            w(i,jstr-1,k,nnew)=w_south(i,k)
#    else
            w(i,jstr-1,k,nnew)=wclm(i,jstr-1,k)
#    endif
#   else
            w(i,jstr-1,k,nnew)=w(i,jstr,k,nnew)
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
            grad(i,jend  )=( w(i  ,jend  ,k,nstp)
     &                      -w(i-1,jend  ,k,nstp))
#   ifdef MASKING
     &                                   *umask(i,jend)
#   endif
            grad(i,jend+1)=( w(i  ,jend+1,k,nstp)
     &                      -w(i-1,jend+1,k,nstp))
#   ifdef MASKING
     &                                 *umask(i,jend+1)
#   endif
          enddo
          do i=istr,iend
#   ifdef LEGACY_ORLANSKI
            dft=w(i,jend,k,nstp)-w(i,jend  ,k,nnew)
            dfx=w(i,jend,k,nnew)-w(i,jend-1,k,nnew)

#    ifdef OBC_RAD_NORMAL
            cy=0. ; cff=max(dfx*dfx, eps)
#    else
            if (dft*(grad(i,jend)+grad(i+1,jend)) > 0.) then
              dfy=grad(i,jend)
            else
              dfy=grad(i+1,jend)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
#     ifdef OBC_RAD_NPO
            cy=0.
#     else
            cy=min(cff,max(dft*dfy,-cff))
#     endif
#    endif

            cx=dft*dfx
            if (cx < 0.) then     !--> cancel, if inflow
              cx=0. ; cy=0.
#    if defined T_FRC_BRY || defined TNUDGING
              cext=dtfwd*ubind * pn(i,jend+1)
            else
              cext=0.
#    endif
            endif

            w(i,jend+1,k,nnew)=( cff*w(i,jend+1,k,nstp)
     &                                 +cx*w(i,jend  ,k,nnew)
     &                                -max(cy,0.)*grad(i  ,jend+1)
     &                                -min(cy,0.)*grad(i+1,jend+1)
     &                                                 )/(cff+cx)

#   else
            cx=dtfwd*v(i,jend+1,k,nrhs)*pn(i,jend+1)
            cy=0.5*dtfwd*(u(i,jend+1,k,nrhs)+u(i+1,jend+1,k,nrhs))
     &                                               *pm(i,jend+1)
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            w(i,jend+1,k,nnew)=(1.-cx)*( w(i,jend+1,k,nstp)
     &                                    -max(cy,0.)*grad(i  ,jend+1)
     &                                    -min(cy,0.)*grad(i+1,jend+1)
     &                                                              )
     &                                  +cx*(  w(i,jend,k,nnew)
     &                                      -max(cy,0.)*grad(i  ,jend)
     &                                      -min(cy,0.)*grad(i+1,jend)
     &                                                              )
#   endif /* LEGACY_ORLANSKI */

#   if defined T_FRC_BRY || defined TNUDGING
            w(i,jend+1,k,nnew)=(1.-cext)*w(i,jend+1,k,nnew)
#    ifdef T_FRC_BRY
     &                                        +cext*w_north(i,k)
#    else
     &                                    +cext*wclm(i,jend+1,k)
#    endif
#   endif
#   ifdef MASKING
            w(i,jend+1,k,nnew)=w(i,jend+1,k,nnew)
     &                                     *rmask(i,jend+1)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=istr,iend
#   if defined OBC_NORTH && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            w(i,jend+1,k,nnew)=w_north(i,k)
#    else
            w(i,jend+1,k,nnew)=wclm(i,jend+1,k)
#    endif
#   else
            w(i,jend+1,k,nnew)=w(i,jend,k,nnew)
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
            w(istr-1,jstr-1,k,nnew)=cff*(
     &              rmask(istr,jstr-1)*w(istr,jstr-1,k,nnew)
     &             +rmask(istr-1,jstr)*w(istr-1,jstr,k,nnew))
          enddo
        else
          do k=1,N
            w(istr-1,jstr-1,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(istr-1,jstr-1,k,nnew)=0.5*( w(istr,jstr-1,k,nnew)
     &                                + w(istr-1,jstr,k,nnew))
        enddo
#   endif
      endif

      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(iend,jstr-1)+rmask(iend+1,jstr)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            w(iend+1,jstr-1,k,nnew)=cff*(
     &              rmask(iend,jstr-1)*w(iend,jstr-1,k,nnew)
     &             +rmask(iend+1,jstr)*w(iend+1,jstr,k,nnew))
          enddo
        else
          do k=1,N
            w(iend+1,jstr-1,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(iend+1,jstr-1,k,nnew)=0.5*( w(iend,jstr-1,k,nnew)
     &                                + w(iend+1,jstr,k,nnew) )
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(istr,jend+1)+rmask(istr-1,jend)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            w(istr-1,jend+1,k,nnew)=cff*(
     &              rmask(istr,jend+1)*w(istr,jend+1,k,nnew)
     &             +rmask(istr-1,jend)*w(istr-1,jend,k,nnew))
          enddo
        else
          do k=1,N
            w(istr-1,jend+1,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(istr-1,jend+1,k,nnew)=0.5*( w(istr,jend+1,k,nnew)
     &                                + w(istr-1,jend,k,nnew))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(iend,jend+1)+rmask(iend+1,jend)
        if (cff > 0.) then
          cff=1./cff
          do k=1,N
            w(iend+1,jend+1,k,nnew)=cff*(
     &              rmask(iend,jend+1)*w(iend,jend+1,k,nnew)
     &             +rmask(iend+1,jend)*w(iend+1,jend,k,nnew))
          enddo
        else
          do k=1,N
            w(iend+1,jend+1,k,nnew)=0.
          enddo
        endif
#   else
        do k=1,N
          w(iend+1,jend+1,k,nnew)= 0.5*( w(iend,jend+1,k,nnew)
     &                                 + w(iend+1,jend,k,nnew))
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
