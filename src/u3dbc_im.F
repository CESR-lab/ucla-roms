#include "cppdefs.opt"

      subroutine u3dbc(grad)

! Set lateral boundary conditions for x-component velocity
! u(:,:,:,nnew).

      use param
      use boundary
      use grid
      use ocean_vars
      use scalars
      use obc_tune

      implicit none

      ! import/export
      real,dimension(GLOBAL_2D_ARRAY) :: grad

      ! Local
      integer  :: i,j,k
      real :: cx,cy, cext, cff, dtfwd, dft,dfx,dfy
      real, parameter :: eps=1.E-33

! For open baroclinic bc's: m3orlanski, m3specified, no-gradient

# if defined OBC_M3ORLANSKI
      logical,parameter ::obc_m3orlanski  = .true.
      logical,parameter ::obc_m3specified = .false.
# elif defined OBC_M3SPECIFIED
      logical,parameter ::obc_m3orlanski  = .false.
      logical,parameter ::obc_m3specified = .true.
# else   /* no-gradient at obc */
      logical,parameter ::obc_m3orlanski  = .false.
      logical,parameter ::obc_m3specified = .false.
# endif

      if (PRED_STAGE) then    ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

! Normal faces: east & west -- u is the normal component

      if (west_bnd) then
        if (obc_west) then
          if (obc_m3orlanski) then
            do k=1,nz
              do j=1,ny+1
                grad(1,j)=(u(1,j,k,nstp)-u(1,j-1,k,nstp))*pmask(1,j)
                grad(2,j)=(u(2,j,k,nstp)-u(2,j-1,k,nstp))*pmask(2,j)
              enddo
              do j=1,ny
                dft=u(2,j,k,nstp)-u(2,j,k,nnew)
                dfx=u(2,j,k,nnew)-u(3,j,k,nnew)

                if (dft*(grad(2,j)+grad(2,j+1)) > 0.) then
                  dfy=grad(2,j)
                else
                  dfy=grad(2,j+1)
                endif
                cff=max(dfx*dfx+dfy*dfy, eps)
                cy=min(cff,max(dft*dfy,-cff))

                cx=dft*dfx
                if (cx < 0.) then     !--> cancel, if inflow
                  cx=0. ; cy=0.
                  if (u_west(j,k) > 0.) then 
                    cext=u_west(j,k)
                  else
                    cext=ubind
                  endif
                  cext=cext * dtfwd * 0.5*(pm(0,j)+pm(1,j))
                else
                  cext=0.
                endif

                u(1,j,k,nnew)=( cff*u(1,j,k,nstp)
     &                             +cx*u(2,j,k,nnew)
     &                          -max(cy,0.)*grad(1,j  )
     &                          -min(cy,0.)*grad(1,j+1)
     &                                         )/(cff+cx)

                if (obc_tuning) cext = max(cext,min(ub_west(j),1.0))

                u(1,j,k,nnew)=(1.-cext)*u(1,j,k,nnew)
     &                                +cext*u_west(j,k)

                u(1,j,k,nnew)=u(1,j,k,nnew)*umask(1,j)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do j=1,ny
                u(1,j,k,nnew)=u_west(j,k)*umask(1,j)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do j=1,ny
                u(1,j,k,nnew)=u(2,j,k,nnew)*umask(1,j)
              enddo
            enddo
          endif
        else  ! closed boundary
          do k=1,nz    
            do j=1,ny
              u(1,j,k,nnew)=0.   
            enddo
          enddo
        endif  ! obc_west
      endif   ! west_bnd

      if (east_bnd) then
        if (obc_east) then
          if (obc_m3orlanski) then
            do k=1,nz    
              do j=1,ny+1 
                grad(nx  ,j)=(u(nx  ,j,k,nstp)-u(nx  ,j-1,k,nstp))*pmask(nx  ,j)
                grad(nx+1,j)=(u(nx+1,j,k,nstp)-u(nx+1,j-1,k,nstp))*pmask(nx+1,j)
              enddo
              do j=1,ny
                dft=u(nx,j,k,nstp)-u(nx  ,j,k,nnew)
                dfx=u(nx,j,k,nnew)-u(nx-1,j,k,nnew)

                if (dft*(grad(nx,j)+grad(nx,j+1)) > 0.) then
                  dfy=grad(nx,j)
                else
                  dfy=grad(nx,j+1)
                endif
                cff=max(dfx*dfx+dfy*dfy,eps)
                cy=min(cff,max(dft*dfy,-cff))
                cx=dft*dfx
                if (cx < 0.) then
                  cx=0. ; cy=0.
                  if (u_east(j,k) < 0.) then
                    cext=-u_east(j,k)
                  else
                    cext=ubind
                  endif
                  cext=cext*dtfwd*0.5*(pm(nx,j)+pm(nx+1,j))
                else
                  cext=0.
                endif

                u(nx+1,j,k,nnew)=( cff*u(nx+1,j,k,nstp)
     &                              +cx*u(nx,j,k,nnew)
     &                       -max(cy,0.)*grad(nx+1,j  )
     &                       -min(cy,0.)*grad(nx+1,j+1)
     &                                        )/(cff+cx)

                if (obc_tuning) cext = max(cext,min(ub_east(j),1.0))

                u(nx+1,j,k,nnew)=(1.-cext)*u(nx+1,j,k,nnew)
     &                                    +cext*u_east(j,k)
                u(nx+1,j,k,nnew)=u(nx+1,j,k,nnew)*umask(nx+1,j)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do j=1,ny
                u(nx+1,j,k,nnew)=u_east(j,k)*umask(nx+1,j)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do j=1,ny
                u(nx+1,j,k,nnew)=u(nx,j,k,nnew)*umask(nx+1,j)
              enddo
            enddo
          endif
        else  ! closed boundary
          do k=1,nz
            do j=1,ny
              u(nx+1,j,k,nnew)=0.
            enddo
          enddo
        endif  ! obc_east
      endif   ! east_bnd

! Tangential faces: south & north -- u is the tangential component

      if (south_bnd) then
        if (obc_south) then
          if (obc_m3orlanski) then
            do k=1,nz
              do i=0,nx
                grad(i,0)=u(i+1,0,k,nstp)-u(i,0,k,nstp)
                grad(i,1)=u(i+1,1,k,nstp)-u(i,1,k,nstp)
              enddo
              do i=1,nx

                cx=-0.125*dtfwd*(v(i,1,k,nrhs)+v(i-1,1,k,nrhs))
     &                                *( pn(i,0)+pn(i-1,0)
     &                                  +pn(i,1)+pn(i-1,1) )

                cy= 0.125*dtfwd*(u(i,0,k,nrhs)+u(i,1,k,nrhs))
     &                            *( pm(i,0)+pm(i-1,0)
     &                                +pm(i,1)+pm(i-1,1) )

                if (cx > 0.) then
                  cext=0.
                else
                  cext=-cx ; cx=0.
                endif

                u(i,0,k,nnew)=(1.-cx)*( u(i,0,k,nstp)
     &                                -max(cy,0.)*grad(i-1,0)
     &                                -min(cy,0.)*grad(i  ,0) )
     &                            +cx*( u(i,1,k,nstp)
     &                                -max(cy,0.)*grad(i-1,1)
     &                                -min(cy,0.)*grad(i  ,1) )

                if (obc_tuning) cext = max(cext,min(ub_south(i),1.0))

                u(i,0,k,nnew)=(1.-cext)*u(i,0,k,nnew)
     &                                  +cext*u_south(i,k)
                u(i,0,k,nnew)=u(i,0,k,nnew)*umask(i,0)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do i=1,nx
                u(i,0,k,nnew)=u_south(i,k)*umask(i,0)
              enddo
            enddo
          else   ! no gradient
            do k=1,nz
              do i=1,nx
                u(i,0,k,nnew)=u(i,1,k,nnew)*umask(i,0)
              enddo
            enddo
          endif
        else   ! slip or no-slip at closed boundary
          do k=1,nz
            do i=1,i1
              u(i,0,k,nnew)=gamma2*u(i,1,k,nnew)*umask(i,0)
            enddo
          enddo
        endif  ! obc_south
      endif   ! south_bnd

      if (north_bnd) then
        if (obc_north) then
          if (obc_m3orlanski) then
            do k=1,nz
              do i=0,nx
                grad(i,ny  )=u(i+1,ny  ,k,nstp)-u(i,ny  ,k,nstp)
                grad(i,ny+1)=u(i+1,ny+1,k,nstp)-u(i,ny+1,k,nstp)
              enddo
              do i=1,nx

                cx=0.125*dtfwd*(v(i,ny+1,k,nrhs)+v(i-1,ny+1,k,nrhs))
     &                                  *( pn(i,ny+1)+pn(i-1,ny+1)
     &                                    +pn(i,ny  )+pn(i-1,ny  ) )

                cy=0.125*dtfwd*(u(i,ny  ,k,nrhs)+u(i  ,ny+1,k,nrhs))
     &                                  *( pm(i,ny+1)+pm(i-1,ny+1)
     &                                    +pm(i,ny  )+pm(i-1,ny  ) )

                if (cx > 0.) then
                  cext=0.
                else
                  cext=-cx ; cx=0.
                endif

                u(i,ny+1,k,nnew)=(1.-cx)*(  u(i,ny+1,k,nstp)
     &                          -max(cy,0.)*grad(i-1,ny+1)
     &                          -min(cy,0.)*grad(i  ,ny+1) )
     &                       +cx*(         u(i,ny,k,nstp)
     &                            -max(cy,0.)*grad(i-1,ny)
     &                            -min(cy,0.)*grad(i  ,ny) )

                if (obc_tuning) cext = max(cext,min(ub_north(i),1.0))

                u(i,ny+1,k,nnew)=(1.-cext)*u(i,ny+1,k,nnew)
     &                                   +cext*u_north(i,k)

                u(i,ny+1,k,nnew)=u(i,ny+1,k,nnew)*umask(i,ny+1)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do i=1,nx
                u(i,ny+1,k,nnew)=u_north(i,k)*umask(i,ny+1)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do i=1,nx
                u(i,ny+1,k,nnew)=u(i,ny,k,nnew)*umask(i,ny+1)
              enddo
            enddo
          endif
        else ! slip or no-slip at closed boundary
          do k=1,nz
            do i=1,i1
              u(i,ny+1,k,nnew)=gamma2*u(i,ny,k,nnew)
     &                                *umask(i,ny+1)
            enddo
          enddo
        endif  ! obc_north
      endif   ! north_bnd
                           
! Corners between adjacent open boundaries

      if (west_bnd.and.south_bnd.and.obc_west.and.obc_south) then
        do k=1,nz
          u(1,0,k,nnew)=0.5*( u(2,0,k,nnew)+u(1,1,k,nnew))
        enddo
      endif
      if (east_bnd.and.south_bnd.and.obc_east.and.obc_south) then
        do k=1,nz
          u(nx+1,0,k,nnew)=0.5*( u(nx,0,k,nnew)+u(nx+1,1,k,nnew))
        enddo
      endif
      if (west_bnd.and.north_bnd.and.obc_west.and.obc_north) then
        do k=1,nz
          u(1   ,ny+1,k,nnew)=0.5*(u(2 ,ny+1,k,nnew)+u(1   ,ny,k,nnew))
        enddo
      endif
      if (east_bnd.and.north_bnd.and.obc_east.and.obc_north) then
        do k=1,nz
          u(nx+1,ny+1,k,nnew)=0.5*(u(nx,ny+1,k,nnew)+u(nx+1,ny,k,nnew))
        enddo
      endif

      end subroutine u3dbc
