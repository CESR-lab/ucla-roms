#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine v3dbc_tile (grad)

! Set lateral boundary conditions for ETA-component velocity
! v(:,:,:,nnew)

      use param
      use wec_frc
      use boundary
      use grid
      use ocean_vars
      use scalars
      use obc_tune

      implicit none
      integer :: i,j,k
      real grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cext, cff,
     &                                     dtfwd, dft,dfx,dfy
      real, parameter :: eps=1.E-33

      if (PRED_STAGE) then    ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

# define dt illegal

# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  ifdef OBC_SOUTH
#   ifdef OBC_M3ORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=1,nx+1                  ! ======== ==== =========
            grad(i,1  )=(v(i,1  ,k,nstp)-v(i-1,1  ,k,nstp))
#    ifdef MASKING
     &                                                *pmask(i,1)
#    endif
            grad(i,2)=(v(i,2,k,nstp)-v(i-1,2,k,nstp))
#    ifdef MASKING
     &                                              *pmask(i,2)
#    endif
          enddo
          do i=1,nx
            dft=v(i,2,k,nstp)-v(i,2,k,nnew)
            dfx=v(i,2,k,nnew)-v(i,3,k,nnew)

#    ifdef OBC_RAD_NORMAL
            cy=0. ; cff=max(dfx*dfx, eps)
#    else
            if (dft*(grad(i,2)+grad(i+1,2)) > 0.) then
              dfy=grad(i,2)
            else
              dfy=grad(i+1,2)
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
#    if defined M3_FRC_BRY
              if (v_south(i,k) > 0.) then
                cext=v_south(i,k)
              else
                cext=ubind
              endif
              cext=cext * dtfwd * 0.5*(pn(i,0)+pn(i,1))
            else
              cext=0.
#    endif
            endif

            v(i,1,k,nnew)=( cff*v(i,1,k,nstp)
     &                        +cx*v(i,2,k,nnew)
     &                     -max(cy,0.)*grad(i  ,1)
     &                     -min(cy,0.)*grad(i+1,1)
     &                                    )/(cff+cx)
#    if defined M3_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_south(i),1.0))
            v(i,1,k,nnew)=(1.-cext)*v(i,1,k,nnew)
     &                               +cext*v_south(i,k)
#    endif
#    ifdef MASKING
            v(i,1,k,nnew)=v(i,1,k,nnew)*vmask(i,1)
#    endif
          enddo
        enddo

#   else           /* alternative open */
        do k=1,N
          do i=1,nx
#    ifdef OBC_M3SPECIFIED
            v(i,1,k,nnew)=v_south(i,k)        ! specified
#    else
            v(i,1,k,nnew)=v(i,2,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                           *vmask(i,1)
#    endif
          enddo
        enddo
#   endif
#  else
        do k=1,N                               ! Southern edge closed
          do i=1,nx                       ! ======== ==== ======
            v(i,1,k,nnew)=0.                ! (no-flux, default)
          enddo
        enddo
#  endif              /* OBC_SOUTH */
      endif         !<-- SOUTHERN_EDGE



      if (NORTHERN_EDGE) then
#  ifdef OBC_NORTH
#   ifdef OBC_M3ORLANSKI
        do k=1,N                            ! Northern edge radiation
          do i=1,nx+1                  ! ======== ==== =========
            grad(i,ny  )=(v(i,ny  ,k,nstp)-v(i-1,ny  ,k,nstp))
#    ifdef MASKING
     &                                                *pmask(i,ny)
#    endif
            grad(i,ny+1)=(v(i,ny+1,k,nstp)-v(i-1,ny+1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(i,ny+1)
#    endif
          enddo
          do i=1,nx
            dft=v(i,ny,k,nstp)-v(i,ny  ,k,nnew)
            dfx=v(i,ny,k,nnew)-v(i,ny-1,k,nnew)

#    ifdef OBC_RAD_NORMAL
            cy=0. ; cff=max(dfx*dfx, eps)
#    else
            if (dft*(grad(i,ny)+grad(i+1,ny)) > 0.) then
              dfy=grad(i,ny)
            else
              dfy=grad(i+1,ny)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
#     ifdef OBC_RAD_NPO
            cy=0.
#     else
            cy=min(cff, max(dft*dfy, -cff))
#     endif
#    endif

            cx=dft*dfx
            if (cx < 0.) then     !--> cancel, if inflow
              cx=0. ; cy=0.
#    if defined M3_FRC_BRY || defined M3NUDGING
              if (v_north(i,k) < 0.) then
                cext=-v_north(i,k)
              else
                cext=ubind
              endif
              cext=cext * dtfwd * 0.5*(pn(i,ny)+pn(i,ny+1))
            else
              cext=0.
#    endif
            endif

            v(i,ny+1,k,nnew)=( cff*v(i,ny+1,k,nstp)
     &                              +cx*v(i,ny,k,nnew)
     &                      -max(cy,0.)*grad(i  ,ny+1)
     &                      -min(cy,0.)*grad(i+1,ny+1)
     &                                      )/(cff+cx)
#    if defined M3_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_north(i),1.0))
            v(i,ny+1,k,nnew)=(1.-cext)*v(i,ny+1,k,nnew)
     &                                   +cext*v_north(i,k)
#    endif
#    ifdef MASKING
            v(i,ny+1,k,nnew)=v(i,ny+1,k,nnew)*vmask(i,ny+1)
#    endif
          enddo
        enddo
#   else                /* alternative open */
        do k=1,N
          do i=1,nx
#    ifdef OBC_M3SPECIFIED
#     ifdef M3_FRC_BRY
            v(i,ny+1,k,nnew)=v_north(i,k)      ! specified
#     else
            v(i,ny+1,k,nnew)=vclm(i,ny+1,k)
#     endif
#    else
            v(i,ny+1,k,nnew)=v(i,ny,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *vmask(i,ny+1)
#    endif
          enddo
        enddo
#   endif
#  else
        do k=1,N                               ! Northern edge closed
          do i=1,nx                       ! ======== ==== ======
#   ifdef WEC 
cc            v(i,ny+1,k,nnew)=-vst(i,ny+1,k) ! no Lagrangian flux
            v(i,ny+1,k,nnew)=v(i,ny,k,nnew)  ! gradient (default)
#    ifdef MASKING
     &                        *vmask(i,ny+1)
#    endif
#   else
            v(i,ny+1,k,nnew)=0.               ! (no-flux, default)
#   endif          
          enddo
        enddo
#  endif
      endif     !<--  NORTHERN_EDGE
# endif          /* !NS_PERIODIC */



# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
#  ifdef OBC_WEST
#   ifdef OBC_M3ORLANSKI
        do k=1,N                             ! Western edge radiation
          do j=0,ny                  ! ======= ==== =========
            grad(0,j)=v(0,j+1,k,nstp)-v(0,j,k,nstp)
            grad(1  ,j)=v(1  ,j+1,k,nstp)-v(1  ,j,k,nstp)
          enddo
          do j=1,ny

          cx=-0.125*dtfwd*(u(1,j,k,nrhs)+u(1,j-1,k,nrhs))
     &                            *( pm(0,j)+pm(0,j-1)
     &                                +pm(1,j)+pm(1,j-1) )

          cy= 0.125*dtfwd*(v(0,j,k,nrhs)+v(1,j,k,nrhs))
     &                            *( pn(0,j)+pn(0,j-1)
     &                                +pn(1,j)+pn(1,j-1) )
          if (cx > 0.) then
            cext=0.
          else
            cext=-cx ; cx=0.
          endif

          v(0,j,k,nnew)=(1.-cx)*( v(0,j,k,nstp)
     &                          -max(cy,0.)*grad(0,j-1)
     &                          -min(cy,0.)*grad(0,j  )
     &                                                     )
     &                       +cx*(        v(1,j,k,nstp)
     &                            -max(cy,0.)*grad(1,j-1)
     &                            -min(cy,0.)*grad(1,j  )
     &                                                     )

#    if defined M3_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_west(j),1.0))
            v(0,j,k,nnew)=(1.-cext)*v(0,j,k,nnew)
     &                                    +cext*v_west(j,k)
#    endif
#    ifdef MASKING
            v(0,j,k,nnew)=v(0,j,k,nnew)*vmask(0,j)
#    endif
          enddo
        enddo
#   else
        do k=1,N
          do j=1,ny
#    ifdef OBC_M3SPECIFIED
#     ifdef M3_FRC_BRY
            v(0,j,k,nnew)=v_west(j,k)       ! specified
#     else
            v(0,j,k,nnew)=vclm(0,j,k)
#     endif
#    else
            v(0,j,k,nnew)=v(1,j,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *vmask(0,j)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef NS_PERIODIC
#    define J_RANGE 1,ny
#   else
#    define J_RANGE 1,j1
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do j=J_RANGE                  ! =====   no-slip (gamma2=-1)
            v(0,j,k,nnew)=gamma2*v(1,j,k,nnew)
#   ifdef MASKING
     &                                *vmask(0,j)
#   endif
          enddo
        enddo
#   undef J_RANGE
#  endif
      endif          !<-- WESTERN_EDGE

      if (EASTERN_EDGE) then
#  ifdef OBC_EAST
#   ifdef OBC_M3ORLANSKI
        do k=1,N                             ! Eastern edge radiation
          do j=0,ny                  ! ======= ==== =========
            grad(nx  ,j)=v(nx  ,j+1,k,nstp)-v(nx  ,j,k,nstp)
            grad(nx+1,j)=v(nx+1,j+1,k,nstp)-v(nx+1,j,k,nstp)
          enddo
          do j=1,ny

            cx=0.125*dtfwd*(u(nx+1,j,k,nrhs)+u(nx+1,j-1,k,nrhs))
     &                               *( pm(nx+1,j)+pm(nx+1,j-1)
     &                                   +pm(nx,j)+pm(nx,j-1) )

            cy=0.125*dtfwd*(v(nx,j,k,nrhs)+v(nx+1,j,k,nrhs))
     &                               *( pn(nx+1,j)+pn(nx+1,j-1)
     &                                   +pn(nx,j)+pn(nx,j-1) )
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif

            v(nx+1,j,k,nnew)=(1.-cx)*(   v(nx+1,j,k,nstp)
     &                          -max(cy,0.)*grad(nx+1,j-1)
     &                          -min(cy,0.)*grad(nx+1,j  )
     &                                                     )
     &                       +cx*(         v(nx,j,k,nstp)
     &                            -max(cy,0.)*grad(nx,j-1)
     &                            -min(cy,0.)*grad(nx,j  )
     &                                                     )

#    if defined M3_FRC_BRY 
            if (obc_tuning) cext = max(cext,min(ub_east(j),1.0))
            v(nx+1,j,k,nnew)=(1.-cext)*v(nx+1,j,k,nnew)
     &                                    +cext*v_east(j,k)
#    endif
#    ifdef MASKING
            v(nx+1,j,k,nnew)=v(nx+1,j,k,nnew)*vmask(nx+1,j)
#    endif
          enddo
        enddo
#   else               /* alternative open */
        do k=1,N
          do j=1,ny
#    ifdef OBC_M3SPECIFIED
#     ifdef M3_FRC_BRY
            v(nx+1,j,k,nnew)=v_east(j,k)       ! specified
#     else
            v(nx+1,j,k,nnew)=vclm(nx+1,j,k)
#     endif
#    else
            v(nx+1,j,k,nnew)=v(nx,j,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *vmask(nx+1,j)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef NS_PERIODIC
#    define J_RANGE 1,ny
#   else
#    define J_RANGE 1,j1
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do j=J_RANGE                  ! ====    no-slip (gamma2=-1)
            v(nx+1,j,k,nnew)=gamma2*v(nx,j,k,nnew)
#   ifdef MASKING
     &                                *vmask(nx+1,j)
#   endif
          enddo
        enddo
#   undef J_RANGE
#  endif
      endif     !<-- EASTERN_EDGE
# endif          /* !EW_PERIODIC */

                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========

# if defined OBC_SOUTH && defined OBC_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          v(0,1,k,nnew)=0.5*( v(0,2,k,nnew)
     &                               +v(1  ,1  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_SOUTH && defined OBC_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          v(nx+1,1,k,nnew)=0.5*( v(nx+1,2,k,nnew)
     &                               +v(nx  ,1  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          v(0,ny+1,k,nnew)=0.5*( v(0,ny,k,nnew)
     &                                 +v(1,ny+1,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          v(nx+1,ny+1,k,nnew)=0.5*( v(nx+1,ny,k,nnew)
     &                                 +v(nx,ny+1,k,nnew))
        enddo
      endif
# endif

! Daniel add
!====================================================================
!       NO LAGRANGIAN FLUX CONDITION AT LAND-SEA INTERFACES
!====================================================================
# if defined WEC && defined MASKING
      do i=i0,i1
        do j=1,j1
          do k=1,N
            v(i,j,k,nnew)=v(i,j,k,nnew)*vmask(i,j)
     &                   +vst(i,j,k)*(vmask(i,j)-1.0)
          enddo
        enddo
      enddo
# endif
!     return
      end
#else
      subroutine v3dbc_empty
      end
#endif /* SOLVE3D */

