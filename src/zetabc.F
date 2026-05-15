#include "cppdefs.opt"

      subroutine zetabc(zeta_new)

! Set boundary conditions for free surface zeta_new

      use param
      use boundary
      use grid
      use ocean_vars
      use scalars

      implicit none

      ! import/export
      real,dimension(GLOBAL_2D_ARRAY) :: zeta_new

      ! local
      integer :: i,j
      real :: cx

# if defined OBC_M2FLATHER
      logical,parameter ::obc_m2flather   = .true.
# else
      logical,parameter ::obc_m2flather   = .false.
# endif

      if (west_bnd) then
        if (obc_west.and.obc_m2flather) then
          do j=0,ny
            cx=dtfast*pm(1,j)*sqrt(g*h(1,j))
            zeta_new(0,j)=( (1.-cx)*zeta(0,j,kstp)
     &                               +cx*zeta(1,j,kstp) )
     &                                    *rmask(0,j)
          enddo
        else  ! specified or closed
          do j=0,ny
            zeta_new(0,j)=zeta_new(1,j)*rmask(0,j)
          enddo
        endif
      endif

      if (east_bnd) then
        if (obc_east.and.obc_m2flather) then
          do j=0,ny
            cx=dtfast*pm(nx,j)*sqrt(g*h(nx,j))
            zeta_new(nx+1,j)=( (1.-cx)*zeta(nx+1,j,kstp)
     &                               +cx*zeta(nx,j,kstp) )
     &                                *rmask(nx+1,j)
          enddo
        else  ! specified or closed
          do j=0,ny
            zeta_new(nx+1,j)=zeta_new(nx,j)
     &                        *rmask(nx+1,j)
          enddo
        endif
      endif

      if (south_bnd) then
        if (obc_south.and.obc_m2flather) then
          do i=0,nx
            cx=dtfast*pn(i,1)*sqrt(g*h(i,1))
            zeta_new(i,0)=( (1.-cx)*zeta(i,0,kstp)
     &                               +cx*zeta(i,1,kstp) )
     &                                    *rmask(i,0)
          enddo
        else  ! specified or closed
          do i=0,nx
            zeta_new(i,0)=zeta_new(i,1)
     &                       *rmask(i,0)
          enddo
        endif
      endif


      if (north_bnd) then
        if (obc_north.and.obc_m2flather) then
          do i=0,nx
            cx=dtfast*pn(i,ny)*sqrt(g*h(i,ny))
            zeta_new(i,ny+1)=( (1.-cx)*zeta(i,ny+1,kstp)
     &                               +cx*zeta(i,ny,kstp) )
     &                                      *rmask(i,ny+1)
          enddo
        else  ! specified or closed
          do i=0,nx
            zeta_new(i,ny+1)=zeta_new(i,ny)
     &                       *rmask(i,ny+1)
          enddo
        endif
      endif

      ! Corners between adjacent open boundaries
      if (west_bnd.and.south_bnd.and.obc_west.and.obc_south) then
        zeta_new(   0,0)=0.5*(zeta_new( 1,0)+zeta_new(  0,1))
      endif
      if (east_bnd.and.south_bnd.and.obc_east.and.obc_south) then
        zeta_new(nx+1,0)=0.5*(zeta_new(nx,0)+zeta_new(nx+1,1))
      endif
      if (west_bnd.and.north_bnd.and.obc_west.and.obc_north) then
        zeta_new(   0,ny+1)=0.5*(zeta_new( 1,ny+1)+zeta_new(   0,ny))
      endif
      if (east_bnd.and.north_bnd.and.obc_east.and.obc_north) then
        zeta_new(nx+1,ny+1)=0.5*(zeta_new(nx,ny+1)+zeta_new(nx+1,ny))
      endif

      end subroutine zetabc

