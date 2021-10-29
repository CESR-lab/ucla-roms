! This is an example of analytical grid definition
! It defines topography h, coriolis f,  mask, etc
!
#include "cppdefs.opt"
#ifdef ANA_GRID

      subroutine ana_grid (tile)

      use param
      use hidden_mpi_vars

      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call ana_grid_tile (istr,iend,jstr,jend)
      end

      subroutine ana_grid_tile (istr,iend,jstr,jend)

! Set up model grid using analytical expressions:
!---- -- ----- ---- ----- ---------- ------------
! output: stored in common blocks, see files "scalars" "grid"

! xl,el      Physical dimensions of the computational domain [usually
!                                                measured in  meters];
! h          Model bathymetry [meters, positive] at RHO-points.
! hmin,hmax  Minimum and maximum values of depth of bathymetry [m].
! f          Coriolis parameter (1/seconds) at RHO-points.
! pm,pn      Coordinate transformation metric "m" [1/meters]
!                   associated with the differential distances in
!                   XI- and ETA-directions, both are at RHO-points.
! xp,xr      XI-coordinates [m] at PSI- and RHO-points.
! yp,yr      ETA-coordinates [m] at PSI- and RHO-points.
      use param
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j

! Grid configuration parameters:
!----- ------------- -----------
! Size_XI, Size_ETA  are the physical dimensions of the computational
!                                domain [usually measured in  meters];
! depth      the maximum depth [meters, positive];
! f0,beta    Coriolis parameters to setup a beta-plane [1/s, 1/(m*s)].

      real cff,y, x0,y0,dx,dy
      real dh, shelf, slope, land, coast

      real, parameter ::
     &       depth = HD       ! H set in param.opt


# include "compute_extended_bounds.h"

      dx=100.0                         ! Set grid spacings for
      dy=100.0                         ! Cartesian rectangular grid

C$OMP MASTER                           ! Copy physical dimensions of
      xl=dx*dble(LLm) 
      el=dy*dble(MMm)                  ! the grid into glabally vizible
C$OMP END MASTER                       ! variables (ONE THREAD ONLY).

# ifdef MPI
      x0=dx*dble(iSW_corn) - xl/2  ! Coordinates of south-west
      y0=dy*dble(jSW_corn)             ! corner of MPI subdomain
                                       ! -xl/2 to centre over axis
# else                                
      x0=0. ; y0=0.
# endif

      do j=jstrR,jendR                  ! Setup Cartezian grid
        do i=istrR,iendR                ! (XI,ETA) at PSI- and RHO-
          xp(i,j)=x0+dx* dble(i-1)      ! points and compute metric
          xr(i,j)=x0+dx*(dble(i)-0.5D0) ! transformation coefficients
          yp(i,j)=y0+dy* dble(j-1)      ! pm and pn, which are
          yr(i,j)=y0+dy*(dble(j)-0.5D0) ! uniform in this case.

          pm(i,j)=1./dx
          pn(i,j)=1./dy
        enddo
      enddo

! Set Coriolis parameter [1/s] at RHO-points.

      x0=xl/2.
      y0=el/2.
      do j=jstrR,jendR
        do i=istrR,iendR
!          f(i,j)=f0+beta*( yr(i,j)-y0 )
          f(i,j)=f0
        enddo
      enddo

      do j=jstrR,jendR
        do i=istrR,iendR
          h(i,j)=depth ! constant depth
        enddo
      enddo

      ! Set up land masking

      do j=jstrR,jendR
        do i=istrR,iendR
          rmask(i,j) = 1
        enddo
      enddo
  
      end
#else
      subroutine ana_grid_empty
      end
#endif /* ANA_GRID */
