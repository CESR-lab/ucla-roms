#include "cppdefs.h"
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
! output: stored in common blocks, see files "scalars.h" "grid"

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

      implicit none
      integer istr,iend,jstr,jend, i,j

! Grid configuration parameters:
!----- ------------- -----------
! Size_XI, Size_ETA  are the physical dimensions of the computational
!                                domain [usually measured in  meters];
! depth      the maximum depth [meters, positive];
! f0,beta    Coriolis parameters to setup a beta-plane [1/s, 1/(m*s)].

      real, parameter ::
# if defined ANA_RIVER_USWC
     &       Size_XI = 1.0e4,   Size_ETA= 1.0e4,
     &       depth=5.,          max_depth=100.0,
     &       f0=0.0,            beta=0.
# else
#  error ANA_GRID CPP-switch is set, but grid parameters are missing.
# endif

      real cff,y, x0,y0,dx,dy
# ifdef ANA_RIVER_USWC
     &  , dh, shelf, slope, land, riv_west, riv_east
# endif

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

! Set bathymetry [meters; positive] at RHO-points.

      ! Default version: Setup flat-bottom
!      do j=jstrR,jendR
!        do i=istrR,iendR         ! topography of prescribed depth.
!          h(i,j)=depth           ! This is used by some examples:
!        enddo                    ! Soliton and Wave Radiation test.
!      enddo

      ! Non-default version: ANA_RIVER_USWC (DevinD)
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
      land=el/10       ! Land extends 10% of domain from south
      riv_west=xl*4/10 ! River west bank at 40% from west
      riv_east=xl*6/10 ! River west bank at 60% from west
      do j=jstrR,jendR
        do i=istrR,iendR
          ! Doesn't extend north beyond land
          if(yr(i,j)<land .and. SOUTHERN_EDGE) then
            ! Isn't in river mouth, also need row (j<=jstr) of mask in
            ! river as river forcing can only be from land cells
            if(j<=jstr .or. xr(i,j)<riv_west .or. xr(i,j)>riv_east) then
              rmask(i,j)=0.0
            endif
          endif
        enddo
      enddo

      end
#else
      subroutine ana_grid_empty
      end
#endif /* ANA_GRID */
