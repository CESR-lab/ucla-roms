#include "cppdefs.opt"

      subroutine t3dbc(itrc, grad)

! Set lateral boundary conditions for tracer t(:,:,:,itrc)

      use param
      use tracers, only: t
      use boundary
      use grid
      use ocean_vars
      use scalars
      use obc_tune

      implicit none

      ! import/export
      integer :: itrc
      real,dimension(GLOBAL_2D_ARRAY) :: grad

      ! Local
      integer :: i,j,k
      real :: cx,cy,cff,cext,dtfwd

# if defined OBC_TORLANSKI
      logical,parameter ::obc_torlanski  = .true.
      logical,parameter ::obc_tspecified = .false.
# elif defined OBC_TSPECIFIED
      logical,parameter ::obc_torlanski  = .false.
      logical,parameter ::obc_tspecified = .true.
# else   /* no-gradient at obc */
      logical,parameter ::obc_torlanski  = .false.
      logical,parameter ::obc_tspecified = .false.
# endif

      if (PRED_STAGE) then    ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

      if (west_bnd) then
        if (obc_west.and.obc_torlanski) then
          do k=1,nz
            do j=1,ny+1
              grad(0,j)=(t(0,j,k,nstp,itrc)-t(0,j-1,k,nstp,itrc))
     &                               *vmask(0,j)
              grad(1,j)=(t(1,j,k,nstp,itrc)-t(1,j-1,k,nstp,itrc))
     &                               *vmask(1,j)
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
     &                                  -max(cy,0.)*grad(0,j  )
     &                                  -min(cy,0.)*grad(0,j+1)
     &                                                              )
     &                                +cx*(  t(1,j,k,nstp,itrc)
     &                                  -max(cy,0.)*grad(1,j  )
     &                                  -min(cy,0.)*grad(1,j+1)
     &                                                              )

              if (obc_tuning) cext = max(cext,min(ub_west(j),1.0))

              t(0,j,k,nnew,itrc)=(1.-cext)*t(0,j,k,nnew,itrc)
     &                                    +cext*t_west(j,k,itrc)
              t(0,j,k,nnew,itrc)=t(0,j,k,nnew,itrc)*rmask(0,j)
            enddo
          enddo
        elseif (obc_west.and.obc_tspecified) then
          do k=1,nz
            do j=1,ny
              t(0,j,k,nnew,itrc)=t_west(j,k,itrc)*rmask(0,j)
            enddo
          enddo
        else  ! no gradient for both default obc and closed bnd
          do k=1,nz
            do j=1,ny
              t(0,j,k,nnew,itrc)=t(1,j,k,nnew,itrc)*rmask(0,j)
            enddo
          enddo
        endif  ! Orlanski, specified or no-gradient
      endif   ! west_bnd

      if (east_bnd) then
        if (obc_east.and.obc_torlanski) then
          do k=1,nz
            do j=1,ny+1
              grad(nx  ,j)=( t(nx  ,j  ,k,nstp,itrc)
     &                      -t(nx  ,j-1,k,nstp,itrc))
     &                                  *vmask(nx,j)
              grad(nx+1,j)=( t(nx+1,j  ,k,nstp,itrc)
     &                      -t(nx+1,j-1,k,nstp,itrc))
     &                                *vmask(nx+1,j)
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

              if (obc_tuning) cext = max(cext,min(ub_east(j),1.0))

              t(nx+1,j,k,nnew,itrc)=(1.-cext)*t(nx+1,j,k,nnew,itrc)
     &                                         +cext*t_east(j,k,itrc)

              t(nx+1,j,k,nnew,itrc)=t(nx+1,j,k,nnew,itrc)*rmask(nx+1,j)
            enddo
          enddo
        elseif (obc_east.and.obc_tspecified) then
          do k=1,nz
            do j=1,ny
              t(nx+1,j,k,nnew,itrc)=t_east(j,k,itrc)*rmask(nx+1,j)
            enddo
          enddo
        else  ! no gradient for both default obc and closed bnd
          do k=1,nz
            do j=1,ny
              t(nx+1,j,k,nnew,itrc)=t(nx,j,k,nnew,itrc)*rmask(nx+1,j)
            enddo
          enddo
        endif  ! Orlanski, specified or no-gradient
      endif ! east_bnd

      if (south_bnd) then
        if (obc_south.and.obc_torlanski) then
          do k=1,nz
            do i=1,nx+1
              grad(i,1)=(t(i,1,k,nstp,itrc)-t(i-1,1,k,nstp,itrc))
     &                                     *umask(i,1)
              grad(i,0)=(t(i,0,k,nstp,itrc)-t(i-1,0,k,nstp,itrc))
     &                                     *umask(i,0)
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
     &                                    -max(cy,0.)*grad(i  ,1)
     &                                    -min(cy,0.)*grad(i+1,1)
     &                                                              )

              if (obc_tuning) cext = max(cext,min(ub_south(i),1.0))

              t(i,0,k,nnew,itrc)=(1.-cext)*t(i,0,k,nnew,itrc)
     &                                        +cext*t_south(i,k,itrc)
              t(i,0,k,nnew,itrc)=t(i,0,k,nnew,itrc)*rmask(i,0)
            enddo
          enddo
        elseif (obc_south.and.obc_tspecified) then
          do k=1,nz
            do i=1,nx
              t(i,0,k,nnew,itrc)=t_south(i,k,itrc)*rmask(i,0)
            enddo
          enddo
        else  ! no gradient for both default obc and closed bnd
          do k=1,nz
            do i=1,nx
              t(i,0,k,nnew,itrc)=t(i,1,k,nnew,itrc)*rmask(i,0)
            enddo
          enddo
        endif  ! Orlanski, specified or no-gradient
      endif  ! south_bnd

      if (north_bnd) then
        if (obc_north.and.obc_torlanski) then
          do k=1,nz
            do i=1,nx+1
              grad(i,ny  )=(t(i,ny  ,k,nstp,itrc)-t(i-1,ny  ,k,nstp,itrc))
     &                                 *umask(i,ny  )
              grad(i,ny+1)=(t(i,ny+1,k,nstp,itrc)-t(i-1,ny+1,k,nstp,itrc))
     &                                 *umask(i,ny+1)
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

              if (obc_tuning) cext = max(cext,min(ub_north(i),1.0))

              t(i,ny+1,k,nnew,itrc)=(1.-cext)*t(i,ny+1,k,nnew,itrc)
     &                                        +cext*t_north(i,k,itrc)
 
              t(i,ny+1,k,nnew,itrc)=t(i,ny+1,k,nnew,itrc)*rmask(i,ny+1)
            enddo
          enddo
        elseif (obc_north.and.obc_tspecified) then
          do k=1,nz
            do i=1,nx
              t(i,ny+1,k,nnew,itrc)=t_north(i,k,itrc)*rmask(i,ny+1)
            enddo
          enddo
        else  ! no gradient for both default obc and closed bnd
          do k=1,nz
            do i=1,nx
              t(i,ny+1,k,nnew,itrc)=t(i,ny,k,nnew,itrc)*rmask(i,ny+1)
            enddo
          enddo
        endif  ! Orlanski, specified or no-gradient
      endif  ! north_bnd

! Corner points between adjacent boundaries. Note that because boundary
! conditions for tracers are of gradient nature -- either Neumann, or
! radiation, the corner points have to be set in all cases, even if the
! adjacent boundaries are closed.  This differs from setting boundaries
! for velocities, where there is anisotropy associated with normal and
! tangential components.   In the case when one/or both points near the
! corner is/are masked, the use of masked points is avoided.

      if (west_bnd.and.south_bnd.and.obc_west.and.obc_south) then
        cff=rmask(1,0)+rmask(0,1)
        if (cff > 0.) then
          cff=1./cff
          do k=1,nz
            t(0,0,k,nnew,itrc)=cff*(
     &              rmask(1,0)*t(1,0,k,nnew,itrc)
     &             +rmask(0,1)*t(0,1,k,nnew,itrc))
          enddo
        else
          do k=1,nz
            t(0,0,k,nnew,itrc)=0.
          enddo
        endif
      endif

      if (east_bnd.and.south_bnd.and.obc_east.and.obc_south) then
        cff=rmask(nx,0)+rmask(nx+1,1)
        if (cff > 0.) then
          cff=1./cff
          do k=1,nz
            t(nx+1,0,k,nnew,itrc)=cff*(
     &              rmask(nx,0)*t(nx,0,k,nnew,itrc)
     &             +rmask(nx+1,1)*t(nx+1,1,k,nnew,itrc))
          enddo
        else
          do k=1,nz
            t(nx+1,0,k,nnew,itrc)=0.
          enddo
        endif
      endif

      if (west_bnd.and.north_bnd.and.obc_west.and.obc_north) then
        cff=rmask(1,ny+1)+rmask(0,ny)
        if (cff > 0.) then
          cff=1./cff
          do k=1,nz
            t(0,ny+1,k,nnew,itrc)=cff*(
     &              rmask(1,ny+1)*t(1,ny+1,k,nnew,itrc)
     &             +rmask(0,ny)*t(0,ny,k,nnew,itrc))
          enddo
        else
          do k=1,nz
            t(0,ny+1,k,nnew,itrc)=0.
          enddo
        endif
      endif

      if (east_bnd.and.north_bnd.and.obc_east.and.obc_north) then
        cff=rmask(nx,ny+1)+rmask(nx+1,ny)
        if (cff > 0.) then
          cff=1./cff
          do k=1,nz
            t(nx+1,ny+1,k,nnew,itrc)=cff*(
     &              rmask(nx,ny+1)*t(nx,ny+1,k,nnew,itrc)
     &             +rmask(nx+1,ny)*t(nx+1,ny,k,nnew,itrc))
          enddo
        else
          do k=1,nz
            t(nx+1,ny+1,k,nnew,itrc)=0.
          enddo
        endif
      endif

      end subroutine t3dbc
