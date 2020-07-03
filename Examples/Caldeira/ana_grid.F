#include "cppdefs.h"
#ifdef ANA_GRID

!! This is ana_grid for CALDEIRA NH test case

      subroutine ana_grid (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call ana_grid_tile (istr,iend,jstr,jend)
      end

      subroutine ana_grid_tile (istr,iend,jstr,jend)

! Set up model grid using analytical expressions:
!---- -- ----- ---- ----- ---------- ------------
! output: stored in common blocks, see files "scalars.h" "grid.h"

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

      implicit none
      integer istr,iend,jstr,jend, i,j
# include "param.h"

! Grid configuration parameters:
!----- ------------- -----------
! Size_XI, Size_ETA  are the physical dimensions of the computational
!                                domain [usually measured in  meters];
! depth      the maximum depth [meters, positive];
! f0,beta    Coriolis parameters to setup a beta-plane [1/s, 1/(m*s)].

      real, parameter :: depth = 1300.
      real, parameter :: f0 = 0.5e-4, beta = 0.
      real, parameter :: Size_XI =2.e3, Size_ETA= 2.0e3

      real x0,y0,dx,dy,cff

# include "grid.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

C$OMP MASTER                           ! Copy physical dimensions of
      xl=Size_XI ; el=Size_ETA         ! the grid into glabally visible
C$OMP END MASTER                       ! variables (ONE THREAD ONLY).

      dx=Size_XI/dble(LLm)             ! Set grid spacings for
      dy=Size_ETA/dble(MMm)            ! Cartesian rectangular grid
# ifdef MPI
      x0=dx*dble(iSW_corn)             ! Coordinates of south-west
      y0=dy*dble(jSW_corn)             ! corner of MPI subdomain
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

      x0=Size_XI/2.
      y0=Size_ETA/2.
      do j=jstrR,jendR
        do i=istrR,iendR
          f(i,j)=f0+beta*( yr(i,j)-y0 )
# if defined NONTRAD_COR
!         feta(i,j) = f0*cos(pi/4)
!         fxi(i,j)  = f0*sin(pi/4)
# endif
        enddo
      enddo

! Set bathymetry [meters; positive] at RHO-points.

      cff=(1./5.e2)**2
      do j=jstrR,jendR
        do i=istrR,iendR
          h(i,j)=depth*(1-0.1*exp(-cff*((xr(i,j)-0.5*xl)**2+
     &                                  (yr(i,j)-0.5*el)**2)))
        enddo
      enddo

      end
#else
      subroutine ana_grid_empty
      end
#endif /* ANA_GRID */
