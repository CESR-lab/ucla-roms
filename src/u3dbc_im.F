#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine u3dbc_tile (grad)

! Set lateral boundary conditions for XI-component velocity
! u(:,:,:,nnew).

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

# ifndef EW_PERIODIC
      if (inode==0) then
#  ifdef OBC_WEST
#   ifdef OBC_M3ORLANSKI
        do k=1,nz                         ! Western edge radiation
          do j=1,ny+1                     ! ======= ==== =========
            grad(1,j)=(u(1,j,k,nstp)-u(1,j-1,k,nstp))
#    ifdef MASKING
     &                                                *pmask(1,j)
#    endif
            grad(2,j)=(u(2,j,k,nstp)-u(2,j-1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(2,j)
#    endif
          enddo
          do j=1,ny
            dft=u(2,j,k,nstp)-u(2,j,k,nnew)
            dfx=u(2,j,k,nnew)-u(3,j,k,nnew)

#    ifdef OBC_RAD_NORMAL
            cy=0. ; cff=max(dfx*dfx, eps)
#    else
            if (dft*(grad(2,j)+grad(2,j+1)) > 0.) then
              dfy=grad(2,j)
            else
              dfy=grad(2,j+1)
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
              if (u_west(j,k) > 0.) then 
                cext=u_west(j,k)
              else
                cext=ubind
              endif
              cext=cext * dtfwd * 0.5*(pm(0,j)+pm(1,j))
            else
              cext=0.
#    endif
            endif

            u(1,j,k,nnew)=( cff*u(1,j,k,nstp)
     &                         +cx*u(2,j,k,nnew)
     &                      -max(cy,0.)*grad(1,j  )
     &                      -min(cy,0.)*grad(1,j+1)
     &                                     )/(cff+cx)
#    if defined M3_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_west(j),1.0))
            u(1,j,k,nnew)=(1.-cext)*u(1,j,k,nnew)
     &                                +cext*u_west(j,k)
#    endif
#    ifdef MASKING
            u(1,j,k,nnew)=u(1,j,k,nnew)*umask(1,j)
#    endif
          enddo
        enddo
#   else            /* alternative open */
        do k=1,nz
          do j=1,ny
#    ifdef OBC_M3SPECIFIED
            u(1,j,k,nnew)=u_west(j,k)         ! specified
#    else
            u(1,j,k,nnew)=u(2,j,k,nnew)  ! Gradient: default
#    endif
#    ifdef MASKING
     &                           *umask(1,j)
#    endif
          enddo
        enddo
#   endif
#  else           /* alternative to open */
        do k=1,nz                              ! Western edge closed
          do j=1,ny                       ! ======= ==== ======
            u(1,j,k,nnew)=0.                !  (no-flux: default)
          enddo
        enddo
#  endif              /* OBC_WEST */
      endif         !<-- WESTERN_EDGE


      if (EASTERN_EDGE) then
#  ifdef OBC_EAST
#   ifdef OBC_M3ORLANSKI
        do k=1,nz                            ! Eastern edge radiation
          do j=1,ny+1                   ! ======= ==== =========
            grad(nx  ,j)=(u(nx  ,j,k,nstp)-u(nx  ,j-1,k,nstp))
#    ifdef MASKING
     &                                                *pmask(nx,j)
#    endif
            grad(nx+1,j)=(u(nx+1,j,k,nstp)-u(nx+1,j-1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(nx+1,j)
#    endif
          enddo
          do j=1,ny
            dft=u(nx,j,k,nstp)-u(nx  ,j,k,nnew)
            dfx=u(nx,j,k,nnew)-u(nx-1,j,k,nnew)

#    ifdef OBC_RAD_NORMAL
            cy=0. ; cff=max(dfx*dfx, eps)
#    else
            if (dft*(grad(nx,j)+grad(nx,j+1)) > 0.) then
              dfy=grad(nx,j)
            else
              dfy=grad(nx,j+1)
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
              if (u_east(j,k) < 0.) then
                cext=-u_east(j,k)
              else
                cext=ubind
              endif
              cext=cext * dtfwd *0.5*(pm(nx,j)+pm(nx+1,j))
            else
              cext=0.
#    endif
            endif

            u(nx+1,j,k,nnew)=( cff*u(nx+1,j,k,nstp)
     &                              +cx*u(nx,j,k,nnew)
     &                       -max(cy,0.)*grad(nx+1,j  )
     &                       -min(cy,0.)*grad(nx+1,j+1)
     &                                        )/(cff+cx)
#    if defined M3_FRC_BRY
            if (obc_tuning) cext = max(cext,min(ub_east(j),1.0))
            u(nx+1,j,k,nnew)=(1.-cext)*u(nx+1,j,k,nnew)
     &                                    +cext*u_east(j,k)
#    endif
#    ifdef MASKING
            u(nx+1,j,k,nnew)=u(nx+1,j,k,nnew)*umask(nx+1,j)
#    endif
          enddo
        enddo
#   else                /* alternative open */
        do k=1,nz
          do j=1,ny
#    ifdef OBC_M3SPECIFIED
            u(nx+1,j,k,nnew)=u_east(j,k)       ! specified
#    else
            u(nx+1,j,k,nnew)=u(nx,j,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *umask(nx+1,j)
#    endif
          enddo
        enddo
#   endif
#  else
        do k=1,N                                ! Eastern edge closed
          do j=1,ny                        ! ======= ==== ======
#   ifdef WEC
c            u(nx+1,j,k,nnew)=-ust(nx+1,j,k) ! no Lagrangian flux
            u(nx+1,j,k,nnew)=u(nx,j,k,nnew)  ! gradient (default)
#    ifdef MASKING
     &                         *umask(nx+1,j)
#    endif
#   else
            u(nx+1,j,k,nnew)=0.               !  (no-flux: default)
#   endif   
          enddo
        enddo
#  endif
      endif         !<-- EASTERN_EDGE
# endif             /* !EW_PERIODIC */



# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  ifdef OBC_SOUTH
#   ifdef OBC_M3ORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=0,nx                 ! ======== ==== =========
            grad(i,0)=u(i+1,0,k,nstp)-u(i,0,k,nstp)
            grad(i,1  )=u(i+1,1  ,k,nstp)-u(i,1  ,k,nstp)
          enddo
          do i=1,nx

          cx=-0.125*dtfwd*(v(i,1,k,nrhs)+v(i-1,1,k,nrhs))
     &                            *( pn(i,0)+pn(i-1,0)
     &                                +pn(i,1)+pn(i-1,1) )

          cy= 0.125*dtfwd*(u(i,0,k,nrhs)+u(i,1,k,nrhs))
     &                            *( pm(i,0)+pm(i-1,0)
     &                                +pm(i,1)+pm(i-1,1) )

          if (cx > 0.) then
            cext=0.
          else
            cext=-cx ; cx=0.
          endif

          u(i,0,k,nnew)=(1.-cx)*(   u(i,0,k,nstp)
     &                          -max(cy,0.)*grad(i-1,0)
     &                          -min(cy,0.)*grad(i  ,0)
     &                                                     )
     &                       +cx*(         u(i,1,k,nstp)
     &                            -max(cy,0.)*grad(i-1,1)
     &                            -min(cy,0.)*grad(i  ,1)
     &                                                     )

#    if defined M3_FRC_BRY
           if (obc_tuning) cext = max(cext,min(ub_south(i),1.0))
           u(i,0,k,nnew)=(1.-cext)*u(i,0,k,nnew)
     &                                  +cext*u_south(i,k)
#    endif
#    ifdef MASKING
            u(i,0,k,nnew)=u(i,0,k,nnew)*umask(i,0)
#    endif
          enddo
        enddo
#   else
        do k=1,N
          do i=1,nx
#    ifdef OBC_M3SPECIFIED
            u(i,0,k,nnew)=u_south(i,k)      ! specified
#    else
            u(i,0,k,nnew)=u(i,1,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *umask(i,0)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef EW_PERIODIC
#    define I_RANGE 1,nx
#   else
#    define I_RANGE 1,i1
#   endif
        do k=1,nz                       ! Wall: free-slip (gamma2=+1)
          do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
            u(i,0,k,nnew)=gamma2*u(i,1,k,nnew)
#   ifdef MASKING
     &                                *umask(i,0)
#   endif
          enddo
        enddo
#   undef I_RANGE
#  endif
      endif              !<-- SOUTHERN_EDGE



      if (NORTHERN_EDGE) then
#  ifdef OBC_NORTH
#   ifdef OBC_M3ORLANSKI
        do k=1,N                            ! Northern edge radiation
          do i=0,nx                 ! ======== ==== =========
            grad(i,ny  )=u(i+1,ny  ,k,nstp)-u(i,ny  ,k,nstp)
            grad(i,ny+1)=u(i+1,ny+1,k,nstp)-u(i,ny+1,k,nstp)
          enddo
          do i=1,nx

          cx=0.125*dtfwd*(v(i,ny+1,k,nrhs)+v(i-1,ny+1,k,nrhs))
     &                               *( pn(i,ny+1)+pn(i-1,ny+1)
     &                                   +pn(i,ny)+pn(i-1,ny) )

          cy=0.125*dtfwd*(u(i,ny,k,nrhs)+u(i,ny+1,k,nrhs))
     &                               *( pm(i,ny+1)+pm(i-1,ny+1)
     &                                   +pm(i,ny)+pm(i-1,ny) )

          if (cx > 0.) then
            cext=0.
          else
            cext=-cx ; cx=0.
          endif

          u(i,ny+1,k,nnew)=(1.-cx)*(  u(i,ny+1,k,nstp)
     &                          -max(cy,0.)*grad(i-1,ny+1)
     &                          -min(cy,0.)*grad(i  ,ny+1)
     &                                                     )
     &                       +cx*(         u(i,ny,k,nstp)
     &                            -max(cy,0.)*grad(i-1,ny)
     &                            -min(cy,0.)*grad(i  ,ny)
     &                                                     )

#    if defined M3_FRC_BRY 
            if (obc_tuning) cext = max(cext,min(ub_north(i),1.0))
            u(i,ny+1,k,nnew)=(1.-cext)*u(i,ny+1,k,nnew)
     &                                   +cext*u_north(i,k)
#    endif
#    ifdef MASKING
            u(i,ny+1,k,nnew)=u(i,ny+1,k,nnew)*umask(i,ny+1)
#    endif
          enddo
        enddo
#   else               /* alternative open */
        do k=1,nz
          do i=1,nx
#    ifdef OBC_M3SPECIFIED
            u(i,ny+1,k,nnew)=u_north(i,k)      ! specified
#    else
            u(i,ny+1,k,nnew)=u(i,ny,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *umask(i,ny+1)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef EW_PERIODIC
#    define I_RANGE 1,nx
#   else
#    define I_RANGE 1,i1
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
            u(i,ny+1,k,nnew)=gamma2*u(i,ny,k,nnew)
#   ifdef MASKING
     &                                *umask(i,ny+1)
#   endif
          enddo
        enddo
#   undef I_RANGE
#  endif
      endif   !<-- NORTHERN_EDGE
# endif        /* !NS_PERIODIC */
                           
! Corners between adjacent open boundaries

# if defined OBC_SOUTH && defined OBC_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          u(1,0,k,nnew)=0.5*( u(2,0,k,nnew)
     &                               +u(1  ,1  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_SOUTH && defined OBC_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          u(nx+1,0,k,nnew)=0.5*( u(nx,0,k,nnew)
     &                                 +u(nx+1,1,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          u(1,ny+1,k,nnew)=0.5*( u(2,ny+1,k,nnew)
     &                               +u(1  ,ny  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          u(nx+1,ny+1,k,nnew)=0.5*( u(nx,ny+1,k,nnew)
     &                                 +u(nx+1,ny,k,nnew))
        enddo
      endif
# endif

! Daniel add
!====================================================================
!       NO LAGRANGIAN FLUX CONDITION AT LAND-SEA INTERFACES
!====================================================================
# if defined WEC && defined MASKING
      do i=1,i1
        do j=j0,j1
          do k=1,N
            u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
     &                   +ust(i,j,k)*(umask(i,j)-1.0)
          enddo
        enddo
      enddo
# endif

      end
#else
      subroutine u3dbc_empty
      end
#endif /* SOLVE3D */
