#include "cppdefs.opt"

      subroutine u3dbc_new

!  Set lateral boundary conditions for x-component velocity
!  u(:,:,:,nnew).

      use param
      use dimensions
      use boundary
      use grid
      use ocean3d
      use scalars
      use hidden_mpi_vars
!     use private_scratch, only: grad => A2d10

      implicit none
      integer :: i,j,k
      real    :: cx,cy, cext, cff, dtfwd, dft,dfx,dfy
      real, parameter :: eps=1.E-33
      real, dimension(:,:),allocatable,save :: gradEW,gradNS

!     real,dimension(-1:nx+2,:),pointer :: grad
!     grad => work20(-1:nx+!

      if (.not.allocated(gradEW)) then
        allocate(gradEW(2,0:ny+1))
        allocate(gradNS(0:nx+1,2))
      endif
      ! gradEW 1 is on the boundary, 2 is first interior
      ! gradNS 1 is in the buffer,   2 is first interior


      if (PRED_STAGE) then    ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.


# ifndef EW_PERIODIC
      !! only for non-periodic in EW
      if (.not.west_msg_exch) then
#  ifdef OBC_WEST
#   ifdef OBC_M3ORLANSKI
        do k=1,nz                            ! Western edge radiation
          do j= 1,ny+1                       ! ======= ==== =========
            gradEW(1,j)=(u(1,j,k,nstp)-u(1,j-1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(1  ,j)
#    endif
            gradEW(2,j)=(u(2,j,k,nstp)-u(2,j-1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(1+1,j)
#    endif
          enddo
          do j=1,ny
            dft=u(2,j,k,nstp)-u(2,j,k,nnew)
            dfx=u(2,j,k,nnew)-u(3,j,k,nnew)

            if (dft*(gradEW(2,j)+gradEW(2,j+1)) > 0.) then
              dfy=gradEW(2,j)
            else
              dfy=gradEW(2,j+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cy=min(cff,max(dft*dfy,-cff))

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
     &                         +cx*u(1+1,j,k,nnew)
     &                      -max(cy,0.)*gradEW(1,j  )
     &                      -min(cy,0.)*gradEW(1,j+1)
     &                                     )/(cff+cx)
#    if defined M3_FRC_BRY
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
          do j= 1,ny                           ! ======= ==== ======
            u(1,j,k,nnew)=0.                !  (no-flux: default)
          enddo
        enddo
#  endif              /* OBC_WEST */
      endif         !<-- WESTERN_EDGE


      if (.not.east_msg_exch) then
#  ifdef OBC_EAST
#   ifdef OBC_M3ORLANSKI
        do k=1,nz                            ! Eastern edge radiation
          do j=1,ny+1                   ! ======= ==== =========
            gradEW(2,j)=(u(nx  ,j,k,nstp)-u(nx  ,j-1,k,nstp))
#    ifdef MASKING
     &                                                *pmask(nx,j)
#    endif
            gradEW(1,j)=(u(nx+1,j,k,nstp)-u(nx+1,j-1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(nx+1,j)
#    endif
          enddo
          do j=1,ny
            dft=u(nx,j,k,nstp)-u(nx  ,j,k,nnew)
            dfx=u(nx,j,k,nnew)-u(nx-1,j,k,nnew)

            if (dft*(gradEW(2,j)+gradEW(2,j+1)) > 0.) then
              dfy=gradEW(2,j)
            else
              dfy=gradEW(2,j+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cy=min(cff,max(dft*dfy,-cff))
            cx=dft*dfx
            if (cx < 0.) then     !--> cancel, if inflow
              cx=0. ; cy=0.

              if (u_east(j,k) < 0.) then
                cext=-u_east(j,k)
              else
                cext=ubind
              endif
              cext=cext * dtfwd *0.5*(pm(nx,j)+pm(nx+1,j))
            else
              cext=0.
            endif

            u(nx+1,j,k,nnew)=( cff*u(nx+1,j,k,nstp)
     &                              +cx*u(nx,j,k,nnew)
     &                       -max(cy,0.)*gradEW(1,j  )
     &                       -min(cy,0.)*gradEW(1,j+1)
     &                                        )/(cff+cx)

            u(nx+1,j,k,nnew)=(1.-cext)*u(nx+1,j,k,nnew)
     &                                    +cext*u_east(j,k)
#    ifdef MASKING
            u(nx+1,j,k,nnew)=u(nx+1,j,k,nnew)*umask(nx+1,j)
#    endif
           if (mynode==5.and.k==35.and.j==1) then
!           print *,'gradEW: ', gradEW(1,j),gradEW(1,j+1)
!           print *,'cff,cx: ', cff,cx,cy
!           print *,'cext:  ', cext
           endif
          enddo
        enddo
        if (mynode==5) then
!       print *,'u_east: ', u(nx+1,1,35,nstp),u(nx+1,1,35,nnew)
        endif
#   else                /* alternative open */
        do k=1,nz
          do j=1,ny
#    ifdef OBC_M3SPECIFIED
            u(nx+1,j,k,nnew)=u_east(j,k)     ! specified
#    else
            u(nx+1,j,k,nnew)=u(nx,j,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *umask(nx+1,j)
#    endif
          enddo
        enddo
#   endif
#  else    /* Not OBC_EAST */
        do k=1,nz                               ! Eastern edge closed
          do j=1,ny                             ! ======= ==== ======
            u(nx+1,j,k,nnew)=0.               !  (no-flux: default)
          enddo
        enddo
#  endif
      endif         !<-- EASTERN_EDGE
# endif             /* !EW_PERIODIC */



# ifndef NS_PERIODIC
      if (.not.south_msg_exch) then
#  ifdef OBC_SOUTH
#   ifdef OBC_M3ORLANSKI
        do k=1,nz                        ! Southern edge radiation
          do i=0,nx                      ! ======== ==== =========
            gradNS(i,1)=u(i+1,0,k,nstp)-u(i,0,k,nstp)
            gradNS(i,2)=u(i+1,1,k,nstp)-u(i,1,k,nstp)
          enddo
          do i=1,nx

          cx=-0.125*dtfwd*(v(i,1,k,nrhs)+v(i-1,1,k,nrhs))
     &                            *( pn(i,1-1)+pn(i-1,1-1)
     &                                +pn(i,1)+pn(i-1,1) )

          cy= 0.125*dtfwd*(u(i,1-1,k,nrhs)+u(i,1,k,nrhs))
     &                            *( pm(i,1-1)+pm(i-1,1-1)
     &                                +pm(i,1)+pm(i-1,1) )

          if (cx > 0.) then
            cext=0.
          else
            cext=-cx ; cx=0.
          endif

          u(i,0,k,nnew)=(1.-cx)*(   u(i,0,k,nstp)
     &                       -max(cy,0.)*gradNS(i-1,1)
     &                       -min(cy,0.)*gradNS(i  ,1)
     &                                                     )
     &                       +cx*(         u(i,1,k,nstp)
     &                        -max(cy,0.)*gradNS(i-1,2)
     &                        -min(cy,0.)*gradNS(i  ,2)
     &                                                     )
           u(i,0,k,nnew)=(1.-cext)*u(i,0,k,nnew)
     &                                  +cext*u_south(i,k)
#    ifdef MASKING
     &                        *umask(i,0)
#    endif
          enddo
        enddo
#   else
        do k=1,nz
          do i=1,nx
#    ifdef OBC_M3SPECIFIED
            u(i,0,k,nnew)=u_south(i,k)   ! specified
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
        do k=1,nz                       ! Wall: free-slip (gamma2=+1)
          do i=1,nx                     ! =====   no-slip (gamma2=-1)
            u(i,0,k,nnew)=gamma2*u(i,1,k,nnew)
#   ifdef MASKING
     &                          *umask(i,0)
#   endif
          enddo
        enddo
#  endif
      endif              !<-- SOUTHERN_EDGE



      if (.not. north_msg_exch) then
#  ifdef OBC_NORTH
#   ifdef OBC_M3ORLANSKI
        do k=1,nz                        ! Northern edge radiation
          do i=0,nx                   ! ======== ==== =========
            gradNS(i,1)=u(i+1,ny+1,k,nstp)-u(i,ny+1,k,nstp)
            gradNS(i,2)=u(i+1,ny  ,k,nstp)-u(i,ny  ,k,nstp)
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
     &                          -max(cy,0.)*gradNS(i-1,1)
     &                          -min(cy,0.)*gradNS(i  ,1)
     &                                                     )
     &                       +cx*(         u(i,ny,k,nstp)
     &                          -max(cy,0.)*gradNS(i-1,2)
     &                          -min(cy,0.)*gradNS(i  ,2)
     &                                                     )


            u(i,ny+1,k,nnew)=(1.-cext)*u(i,ny+1,k,nnew)
     &                                   +cext*u_north(i,k)
#    ifdef MASKING
            u(i,ny+1,k,nnew)=u(i,ny+1,k,nnew)*umask(i,ny+1)
#    endif
          enddo
        enddo
#   else               /* alternative open */
        do k=1,nz
          do i=1,nx
            u_fast(i,ny+1,k,nnew)=u_north_fast(i,k)  ! specified
#    ifdef OBC_M3SPECIFIED
            u(i,ny+1,k,nnew)=u_north(i,k)    ! specified
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
        do k=1,nz                       ! Wall: free-slip (gamma2=+1)
          do i=1,nx+1                   ! =====   no-slip (gamma2=-1)
            u(i,ny+1,k,nnew)=gamma2*u(i,ny,k,nnew)
#   ifdef MASKING
     &                                *umask(i,ny+1)
#   endif
          enddo
        enddo
#  endif
      endif   !<-- NORTHERN_EDGE
# endif        /* !NS_PERIODIC */

! Corners between adjacent open boundaries

# if defined OBC_SOUTH && defined OBC_WEST
      if (.not.west_msg_exch .and. .not.south_msg_exch) then
        do k=1,nz
          u(1,1-1,k,nnew)=0.5*( u(1+1,1-1,k,nnew)
     &                               +u(1  ,1  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_SOUTH && defined OBC_EAST
      if (.not.east_msg_exch .and. .not.south_msg_exch) then
        do k=1,N
          u(nx+1,1-1,k,nnew)=0.5*( u(nx,1-1,k,nnew)
     &                                 +u(nx+1,1,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_WEST
      if (.not.west_msg_exch .and. .not.north_msg_exch) then
        do k=1,N
          u(1,ny+1,k,nnew)=0.5*( u(1+1,ny+1,k,nnew)
     &                               +u(1  ,ny  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_EAST
      if (.not.east_msg_exch .and. .not.north_msg_exch) then
        do k=1,N
          u(nx+1,ny+1,k,nnew)=0.5*( u(nx,ny+1,k,nnew)
     &                                 +u(nx+1,ny,k,nnew))
        enddo
      endif
# endif

      end
