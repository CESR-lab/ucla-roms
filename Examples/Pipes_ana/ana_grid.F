! This is an example of analytical grid definition
! It is part of the rivers_ana example and defines
! topography h, coriolis f,  mask, and river flux location rflx
!
#include "cppdefs.h"
#ifdef ANA_GRID

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

      use pipe_frc
      implicit none
      integer istr,iend,jstr,jend, i,j
# include "param.h"

! Grid configuration parameters:
!----- ------------- -----------
! Size_XI, Size_ETA  are the physical dimensions of the computational
!                                domain [usually measured in  meters];
! depth      the maximum depth [meters, positive];
! f0,beta    Coriolis parameters to setup a beta-plane [1/s, 1/(m*s)].

      real, parameter ::
     &       Size_XI = 1.0e4,   Size_ETA= 1.0e4,
     &       depth=5.,          max_depth=100.0,
     &       f0=0.0e-4,            beta=0.

      real cff,y, x0,y0,dx,dy
      real dh, shelf, slope, land, coast
      real riv_west, riv_east,riv_cells
      real psz,px,py,pipe_cells

# include "grid.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

C$OMP MASTER                           ! Copy physical dimensions of
      xl=Size_XI ; el=Size_ETA         ! the grid into glabally vizible
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
        enddo
      enddo

      shelf=size_eta/5 ! shelf location in meters from south
      slope=(max_depth-depth)/(size_eta*4/5) ! Similar triangles o/a=dh/pm=(max_depth-depth)/(MMm*4/5)
      do j=jstrR,jendR
        do i=istrR,iendR

          if(yr(i,j)<shelf .and. SOUTHERN_EDGE) then
            ! Constant shallow region 20% of domain in south.
            h(i,j)=depth
          else
            ! Uniform gradient from south (shallow) to north (deep).
            dh=(yr(i,j)-shelf)*slope
            h(i,j)=depth+dh
          endif

        enddo
      enddo

      ! Set up land masking for river channel
      land  = el*0.1  ! Land extends 10% of domain from south
      coast = el*0.02 ! Coast is not as far
      riv_west=xl*0.4 ! River west bank at 40% from west
      riv_east=xl*0.6 ! River west bank at 60% from west

      do j=jstrR,jendR
        do i=istrR,iendR
          ! default is water
          rmask(i,j) = 1

          if(yr(i,j)<land) then
            if (xr(i,j)<riv_west .or. xr(i,j)>riv_east) then
              rmask(i,j)=0.0
            endif
          endif
          if(yr(i,j)<coast) then !! All land in the far south
            rmask(i,j) = 0.0
          endif
        enddo
      enddo

#if defined PIPE_SOURCE
      psz = xl*0.02 ! Width of the pipe
      px  = xl*.5  ! x location pipe 
      py  = el*.5  ! y location pipe 
      pipe_cells = nint(psz/dx)**2 !number of cells in this pipe
      do j=jstrR,jendR
        do i=istrR,iendR
          pipe_fraction(i,j) = 0.0
          pipe_idx(i,j) = 0
          if (xr(i,j)> px-0.5*psz .and. xr(i,j)<px+0.5*psz) then
            if (yr(i,j)> py-0.5*psz .and. yr(i,j)<py+0.5*psz) then
               pipe_fraction(i,j) = 1.0/pipe_cells
               pipe_idx(i,j) = 1
            endif
          endif
        enddo
      enddo
#endif
   
      end
#else
      subroutine ana_grid_empty
      end
#endif /* ANA_GRID */
