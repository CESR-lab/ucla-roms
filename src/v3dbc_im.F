#include "cppdefs.opt"

      subroutine v3dbc(grad)

! Set lateral boundary conditions for y-component velocity
! v(:,:,:,nnew).

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

! Normal faces: south & north -- v is the normal component

      if (south_bnd) then
        if (obc_south) then
          if (obc_m3orlanski) then
            do k=1,nz
              do i=1,nx+1
                grad(i,1)=(v(i,1,k,nstp)-v(i-1,1,k,nstp))*pmask(i,1)
                grad(i,2)=(v(i,2,k,nstp)-v(i-1,2,k,nstp))*pmask(i,2)
              enddo
              do i=1,nx
                dft=v(i,2,k,nstp)-v(i,2,k,nnew)
                dfx=v(i,2,k,nnew)-v(i,3,k,nnew)

                if (dft*(grad(i,2)+grad(i+1,2)) > 0.) then
                  dfy=grad(i,2)
                else
                  dfy=grad(i+1,2)
                endif
                cff=max(dfx*dfx+dfy*dfy, eps)
                cy=min(cff,max(dft*dfy,-cff))

                cx=dft*dfx
                if (cx < 0.) then     !--> cancel, if inflow
                  cx=0. ; cy=0.
                  if (v_south(i,k) > 0.) then
                    cext=v_south(i,k)
                  else
                    cext=ubind
                  endif
                  cext=cext * dtfwd * 0.5*(pn(i,0)+pn(i,1))
                else
                  cext=0.
                endif

                v(i,1,k,nnew)=( cff*v(i,1,k,nstp)
     &                             +cx*v(i,2,k,nnew)
     &                          -max(cy,0.)*grad(i  ,1)
     &                          -min(cy,0.)*grad(i+1,1)
     &                                         )/(cff+cx)

                if (obc_tuning) cext = max(cext,min(ub_south(i),1.0))

                v(i,1,k,nnew)=(1.-cext)*v(i,1,k,nnew)
     &                                +cext*v_south(i,k)

                v(i,1,k,nnew)=v(i,1,k,nnew)*vmask(i,1)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do i=1,nx
                v(i,1,k,nnew)=v_south(i,k)*vmask(i,1)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do i=1,nx
                v(i,1,k,nnew)=v(i,2,k,nnew)*vmask(i,1)
              enddo
            enddo
          endif
        else  ! closed boundary
          do k=1,nz
            do i=1,nx
              v(i,1,k,nnew)=0.
            enddo
          enddo
        endif  ! obc_south
      endif   ! south_bnd

      if (north_bnd) then
        if (obc_north) then
          if (obc_m3orlanski) then
            do k=1,nz
              do i=1,nx+1
                grad(i,ny  )=(v(i,ny  ,k,nstp)-v(i-1,ny  ,k,nstp))*pmask(i,ny  )
                grad(i,ny+1)=(v(i,ny+1,k,nstp)-v(i-1,ny+1,k,nstp))*pmask(i,ny+1)
              enddo
              do i=1,nx
                dft=v(i,ny,k,nstp)-v(i,ny  ,k,nnew)
                dfx=v(i,ny,k,nnew)-v(i,ny-1,k,nnew)

                if (dft*(grad(i,ny)+grad(i+1,ny)) > 0.) then
                  dfy=grad(i,ny)
                else
                  dfy=grad(i+1,ny)
                endif
                cff=max(dfx*dfx+dfy*dfy,eps)
                cy=min(cff,max(dft*dfy,-cff))
                cx=dft*dfx
                if (cx < 0.) then
                  cx=0. ; cy=0.
                  if (v_north(i,k) < 0.) then
                    cext=-v_north(i,k)
                  else
                    cext=ubind
                  endif
                  cext=cext*dtfwd*0.5*(pn(i,ny)+pn(i,ny+1))
                else
                  cext=0.
                endif

                v(i,ny+1,k,nnew)=( cff*v(i,ny+1,k,nstp)
     &                              +cx*v(i,ny,k,nnew)
     &                       -max(cy,0.)*grad(i  ,ny+1)
     &                       -min(cy,0.)*grad(i+1,ny+1)
     &                                        )/(cff+cx)

                if (obc_tuning) cext = max(cext,min(ub_north(i),1.0))

                v(i,ny+1,k,nnew)=(1.-cext)*v(i,ny+1,k,nnew)
     &                                    +cext*v_north(i,k)
                v(i,ny+1,k,nnew)=v(i,ny+1,k,nnew)*vmask(i,ny+1)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do i=1,nx
                v(i,ny+1,k,nnew)=v_north(i,k)*vmask(i,ny+1)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do i=1,nx
                v(i,ny+1,k,nnew)=v(i,ny,k,nnew)*vmask(i,ny+1)
              enddo
            enddo
          endif
        else  ! closed boundary
          do k=1,nz
            do i=1,nx
              v(i,ny+1,k,nnew)=0.
            enddo
          enddo
        endif  ! obc_north
      endif   ! north_bnd

! Tangential faces: west & east -- v is the tangential component

      if (west_bnd) then
        if (obc_west) then
          if (obc_m3orlanski) then
            do k=1,nz
              do j=0,ny
                grad(0,j)=v(0,j+1,k,nstp)-v(0,j,k,nstp)
                grad(1,j)=v(1,j+1,k,nstp)-v(1,j,k,nstp)
              enddo
              do j=1,ny

                cx=-0.125*dtfwd*(u(1,j,k,nrhs)+u(1,j-1,k,nrhs))
     &                                *( pm(0,j)+pm(0,j-1)
     &                                  +pm(1,j)+pm(1,j-1) )

                cy= 0.125*dtfwd*(v(0,j,k,nrhs)+v(1,j,k,nrhs))
     &                            *( pn(0,j)+pn(0,j-1)
     &                                +pn(1,j)+pn(1,j-1) )

                if (cx > 0.) then
                  cext=0.
                else
                  cext=-cx ; cx=0.
                endif

                v(0,j,k,nnew)=(1.-cx)*( v(0,j,k,nstp)
     &                                -max(cy,0.)*grad(0,j-1)
     &                                -min(cy,0.)*grad(0,j  )
     &                                                     )
     &                            +cx*( v(1,j,k,nstp)
     &                                -max(cy,0.)*grad(1,j-1)
     &                                -min(cy,0.)*grad(1,j  )
     &                                                     )

                if (obc_tuning) cext = max(cext,min(ub_west(j),1.0))

                v(0,j,k,nnew)=(1.-cext)*v(0,j,k,nnew)
     &                                  +cext*v_west(j,k)
                v(0,j,k,nnew)=v(0,j,k,nnew)*vmask(0,j)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do j=1,ny
                v(0,j,k,nnew)=v_west(j,k)*vmask(0,j)
              enddo
            enddo
          else   ! no gradient
            do k=1,nz
              do j=1,ny
                v(0,j,k,nnew)=v(1,j,k,nnew)*vmask(0,j)
              enddo
            enddo
          endif
        else   ! slip or no-slip at closed boundary
          do k=1,nz
            do j=1,j1
              v(0,j,k,nnew)=gamma2*v(1,j,k,nnew)*vmask(0,j)
            enddo
          enddo
        endif  ! obc_west
      endif   ! west_bnd

      if (east_bnd) then
        if (obc_east) then
          if (obc_m3orlanski) then
            do k=1,nz
              do j=0,ny
                grad(nx  ,j)=v(nx  ,j+1,k,nstp)-v(nx  ,j,k,nstp)
                grad(nx+1,j)=v(nx+1,j+1,k,nstp)-v(nx+1,j,k,nstp)
              enddo
              do j=1,ny

                cx=0.125*dtfwd*(u(nx+1,j,k,nrhs)+u(nx+1,j-1,k,nrhs))
     &                                  *( pm(nx+1,j)+pm(nx+1,j-1)
     &                                    +pm(nx  ,j)+pm(nx  ,j-1) )

                cy=0.125*dtfwd*(v(nx  ,j,k,nrhs)+v(nx+1,j  ,k,nrhs))
     &                                  *( pn(nx+1,j)+pn(nx+1,j-1)
     &                                    +pn(nx  ,j)+pn(nx  ,j-1) )

                if (cx > 0.) then
                  cext=0.
                else
                  cext=-cx ; cx=0.
                endif

                v(nx+1,j,k,nnew)=(1.-cx)*(  v(nx+1,j,k,nstp)
     &                          -max(cy,0.)*grad(nx+1,j-1)
     &                          -min(cy,0.)*grad(nx+1,j  )
     &                                                     )
     &                       +cx*(         v(nx,j,k,nstp)
     &                            -max(cy,0.)*grad(nx,j-1)
     &                            -min(cy,0.)*grad(nx,j  )
     &                                                     )

                if (obc_tuning) cext = max(cext,min(ub_east(j),1.0))

                v(nx+1,j,k,nnew)=(1.-cext)*v(nx+1,j,k,nnew)
     &                                   +cext*v_east(j,k)

                v(nx+1,j,k,nnew)=v(nx+1,j,k,nnew)*vmask(nx+1,j)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do j=1,ny
                v(nx+1,j,k,nnew)=v_east(j,k)*vmask(nx+1,j)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do j=1,ny
                v(nx+1,j,k,nnew)=v(nx,j,k,nnew)*vmask(nx+1,j)
              enddo
            enddo
          endif
        else ! slip or no-slip at closed boundary
          do k=1,nz
            do j=1,j1
              v(nx+1,j,k,nnew)=gamma2*v(nx,j,k,nnew)
     &                                *vmask(nx+1,j)
            enddo
          enddo
        endif  ! obc_east
      endif   ! east_bnd

! Corners between adjacent open boundaries

      if (west_bnd.and.south_bnd.and.obc_west.and.obc_south) then
        do k=1,nz
          v(0   ,1,k,nnew)=0.5*(v(0   ,2,k,nnew)+v(1 ,1,k,nnew))
        enddo
      endif
      if (east_bnd.and.south_bnd.and.obc_east.and.obc_south) then
        do k=1,nz
          v(nx+1,1,k,nnew)=0.5*(v(nx+1,2,k,nnew)+v(nx,1,k,nnew))
        enddo
      endif
      if (west_bnd.and.north_bnd.and.obc_west.and.obc_north) then
        do k=1,nz
          v(0   ,ny+1,k,nnew)=0.5*(v(0   ,ny,k,nnew)+v(1 ,ny+1,k,nnew))
        enddo
      endif
      if (east_bnd.and.north_bnd.and.obc_east.and.obc_north) then
        do k=1,nz
          v(nx+1,ny+1,k,nnew)=0.5*(v(nx+1,ny,k,nnew)+v(nx,ny+1,k,nnew))
        enddo
      endif

      end subroutine v3dbc
