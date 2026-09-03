#include "cppdefs.opt"

      subroutine w3dbc(grad)

! Set lateral boundary conditions for vertical flux w

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
      integer :: i,j,k
      real :: cx,cy, cff,cext,dtfwd

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

      if (west_bnd) then
        if (obc_west) then
          if (obc_m3orlanski) then
            do k=1,nz
              do j=1,ny+1
                grad(0,j)= w(0,j,k,nstp)- w(0,j-1,k,nstp)
                grad(1,j)= w(1,j,k,nstp)- w(1,j-1,k,nstp)
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

                if (obc_tuning) cext = max(cext,min(ub_west(j),1.0))

                w(0,j,k,nnew)=((1.-cext)*w(0,j,k,nnew)
     &                                +cext*w_west(j,k))*rmask(0,j)

              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do j=1,ny
                w(0,j,k,nnew)=w_west(j,k)*rmask(0,j)
              enddo
            enddo
          else   ! no gradient
            do k=1,nz
              do j=1,ny
                w(0,j,k,nnew)=w(1,j,k,nnew)*rmask(0,j)
              enddo
            enddo
          endif
        else   ! slip or no-slip at closed boundary
          do k=1,nz
            do j=1,ny
              w(0,j,k,nnew)=gamma2*w(1,j,k,nnew)*rmask(0,j)
            enddo
          enddo
        endif  ! obc_west
      endif   ! west_bnd


      if (east_bnd) then
        if (obc_east) then
          if (obc_m3orlanski) then
            do k=1,nz
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

                if (obc_tuning) cext = max(cext,min(ub_east(j),1.0))

                w(nx+1,j,k,nnew)=((1.-cext)*w(nx+1,j,k,nnew)
     &                             +cext*w_east(j,k)  )*rmask(nx+1,j)

              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do j=1,ny
                w(nx+1,j,k,nnew)=w_east(j,k)*rmask(nx+1,j)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do j=1,ny
                w(nx+1,j,k,nnew)=w(nx,j,k,nnew)*rmask(nx+1,j)
              enddo
            enddo
          endif
        else ! slip or no-slip at closed boundary
          do k=1,nz
            do j=1,ny
              w(nx+1,j,k,nnew)=gamma2*w(nx,j,k,nnew)*rmask(nx+1,j)
            enddo
          enddo
        endif  ! obc_east
      endif   ! east_bnd

      if (south_bnd) then
        if (obc_south) then
          if (obc_m3orlanski) then
            do k=1,nz
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

                if (obc_tuning) cext = max(cext,min(ub_south(i),1.0))

                w(i,0,k,nnew)=((1.-cext)*w(i,0,k,nnew)
     &                             +cext*w_south(i,k) )*rmask(i,0)
              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do i=1,nx
                w(i,0,k,nnew)=w_south(i,k)*rmask(i,0)
              enddo
            enddo
          else   ! no gradient
            do k=1,nz
              do i=1,nx
                w(i,0,k,nnew)=w(i,1,k,nnew)*rmask(i,0)
              enddo
            enddo
          endif
        else   ! slip or no-slip at closed boundary
          do k=1,nz
            do i=1,nx
              w(i,0,k,nnew)=gamma2*w(i,1,k,nnew)*rmask(i,0)
            enddo
          enddo
        endif  ! obc_south
      endif   ! south_bnd


      if (north_bnd) then
        if (obc_north) then
          if (obc_m3orlanski) then
            do k=1,nz
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

                if (obc_tuning) cext = max(cext,min(ub_north(i),1.0))

                w(i,ny+1,k,nnew)= ((1.-cext)*w(i,ny+1,k,nnew)
     &                             +cext*w_north(i,k) )*rmask(i,ny+1)

              enddo
            enddo
          elseif (obc_m3specified) then
            do k=1,nz
              do i=1,nx
                w(i,ny+1,k,nnew)=w_north(i,k)*rmask(i,ny+1)
              enddo
            enddo
          else ! no gradient
            do k=1,nz
              do i=1,nx
                w(i,ny+1,k,nnew)=w(i,ny,k,nnew)*rmask(i,ny+1)
              enddo
            enddo
          endif
        else ! slip or no-slip at closed boundary
          do k=1,nz
            do i=1,nx
              w(i,ny+1,k,nnew)=gamma2*w(i,ny,k,nnew)*rmask(i,ny+1)
            enddo
          enddo
        endif  ! obc_north
      endif   ! north_bnd

! Corner points between adjacent boundaries.
! Unlike tracers, we're fine taking the average of points that
! are potentially masked, since the values are benign (zeros)

      if (west_bnd.and.south_bnd.and.obc_west.and.obc_south) then
        do k=1,nz
          w(0,0,k,nnew)=0.5*( w(1,0,k,nnew)+ w(0,1,k,nnew))
        enddo
      endif

      if (east_bnd.and.south_bnd.and.obc_east.and.obc_south) then
        do k=1,nz
          w(nx+1,0,k,nnew)=0.5*( w(nx,0,k,nnew)+w(nx+1,1,k,nnew) )
        enddo
      endif

      if (west_bnd.and.north_bnd.and.obc_west.and.obc_north) then
        do k=1,nz
          w(0,ny+1,k,nnew)=0.5*( w(1,ny+1,k,nnew)+ w(0,ny,k,nnew))
        enddo
      endif

      if (east_bnd.and.north_bnd.and.obc_east.and.obc_north) then
        do k=1,nz
          w(nx+1,ny+1,k,nnew)=0.5*(w(nx,ny+1,k,nnew)+w(nx+1,ny,k,nnew))
        enddo
      endif

      end subroutine w3dbc

