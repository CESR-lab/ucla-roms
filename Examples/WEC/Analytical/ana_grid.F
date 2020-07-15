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
# if defined ANA_WEC_FRC
     &                  Size_XI = 1400.0,   Size_ETA= 2.0,
     &                  depth=10.,          f0=0.0,        beta=0.
# else
#  error ANA_GRID CPP-switch is set, but grid parameters are missing.
# endif

      real cff,y, x0,y0,dx,dy
# if defined ISWAKE
      real rd, rd_inner, rd_outer, ypL,yC,ypR
# elif defined SEAMOUNT
      real xpL,xC,xpR,  ypL,yC,ypR
# elif defined RIVER
      integer imax,jchn
# endif
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
CR      write(*,'(4(A3,pe15.9,1x),I3)') 'dx=',dx, 'dy=',dy,
CR     &                         'x0=',x0, 'y0=',y0, mynode


# ifdef SEAMOUNT
#  define STRTCH 0.333333333333
      do i=istrR-1,iendR,+1   !--> irreversible
        cff=x0+dx*dble(i)
        xpR=cff +STRTCH*(Size_XI/(2.*pi))*sin(2.*pi*cff/Size_XI)
        if (i > istrR-1) then
          cff=x0+dx*(dble(i)-0.5)
          xC=cff +STRTCH*(Size_XI/(2.*pi))*sin(2.*pi*cff/Size_XI)
          cff=1./(xpR-xpL)
          do j=jstrR,jendR            ! Generate rectangular grid
            xp(i,j)=xpL               ! with grid resolition refined
            xr(i,j)=xC                ! in the middle by approximately
            pm(i,j)=cff               ! factor of 2 relatively to
          enddo                       ! periphery (if STRTCH = 1/3)
        endif
        xpL=xpR
      enddo
      do j=jstrR-1,jendR,+1   !--> irreversible
        cff=y0+dy*dble(j)
        ypR=cff +STRTCH*(Size_ETA/(2.*pi))*sin(2.*pi*cff/Size_ETA)
        if (j > jstrR-1) then
          cff=y0+dy*(dble(j)-0.5)
          yC=cff +STRTCH*(Size_ETA/(2.*pi))*sin(2.*pi*cff/Size_ETA)
          cff=1./(ypR-ypL)
          do i=istrR,iendR
            yp(i,j)=ypL
            yr(i,j)=yC
            pn(i,j)=cff
          enddo
        endif
        ypL=ypR
      enddo
#  undef STRTCH
# elif defined ISWAKE
      cff=0.25*Size_ETA                    ! Note that y-coordinate
      do j=jstrR-1,jendR,+1                ! defined here is centered
        ypR=y0+dy*dble(j) -0.5*Size_ETA    ! at the channel middle axis
        if (ypR < -cff) then
          ypR=ypR -(ypR+cff)**2/(0.5*Size_ETA-cff)
        elseif (ypR > cff) then
          ypR=ypR +(ypR-cff)**2/(0.5*Size_ETA-cff)
        endif

        yC=y0+dy*(dble(j)-0.5D0) -0.5*Size_ETA
        if (yC < -cff) then
          yC=yC -(yC+cff)**2/(0.5*Size_ETA-cff)
        elseif (yC > cff) then
          yC=yC +(yC-cff)**2/(0.5*Size_ETA-cff)
        endif

        if (j > jstrR-1) then
          do i=istrR,iendR
            xp(i,j)=x0+dx* dble(i-1)
            xr(i,j)=x0+dx*(dble(i)-0.5D0)
            pm(i,j)=1./dx

            yp(i,j)=ypL
            yr(i,j)=0.5*yC+0.25*(ypL+ypR)
            pn(i,j)=1./(ypR-ypL)
          enddo
        endif
        ypL=ypR !<-- recursive
      enddo
C$OMP MASTER
      el=2.*(Size_ETA-cff)
C$OMP END MASTER
# else
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
# endif

! Set Coriolis parameter [1/s] at RHO-points.

      x0=Size_XI/2.
      y0=Size_ETA/2.
      do j=jstrR,jendR
        do i=istrR,iendR
          f(i,j)=f0+beta*( yr(i,j)-y0 )
c     &                    +(0.1)*(xr(i,j)-x0))
        enddo
      enddo

! Set bathymetry [meters; positive] at RHO-points.

# if defined CANYON_A || defined CANYON_B
      do j=jstrR,jendR
        do i=istrR,iendR
          cff=32000.-16000.*(sin(pi*xr(i,j)/size_XI))**24
          h(i,j)=20.+0.5*(depth-20.)*(1.+tanh((yr(i,j)-cff)/10000.))
        enddo
      enddo
# elif defined OVERFLOW
      do j=jstrR,jendR
        do i=istrR,iendR
          h(i,j)=20.+0.5*(depth-20.)*(1.+tanh((yr(i,j)-40000.)/5000.))
        enddo
      enddo
# elif defined SEAMOUNT
      y0=Size_ETA/2.          ! Seamount topography: the seamount is
      x0=y0                   ! in the middle of the domain in both
      cff=(1./40000.0)**2     ! directions, if it is a square (unforced
      do j=jstrR,jendR        ! PGF error test) or closer to the west
        do i=istrR,iendR
          h(i,j)=depth-4500.0*exp( -cff*( (xr(i,j)-x0)**2
     &                                   +(yr(i,j)-y0)**2 ))
        enddo
      enddo                   ! if the domain is rectangular (forced)

# elif defined ISWAKE

! Island topography and coastline: two versions are available:
! cylinder-shaped island with vertical walls and island with
! cylindrical above-surface part and gaussian-shaped shelf.
! Note that the above-surface area of the island is smaller in
! the second case because an extra space is needed to accommodate
! shelf areas.  For a 80km-wide domain maximum allowed "rd_inner"
! should not exceed 10 and 7.5 km for cylinder and shelf-shaped
! islands respectively.

      x0=Size_ETA/2. ; y0=0.     ! location of the center of island
c--#  define GAUSSIAN_SLOPE_ISLAND
#  ifdef GAUSSIAN_SLOPE_ISLAND
      rd_inner=7.5e+3           ! radius of the above-surface part
#  else
      rd_inner=10.e+3           ! radius of the cylinder-shaped island
#  endif
      rd_outer=rd_inner*1.
      do j=jstrR,jendR
        do i=istrR,iendR
          h(i,j)=depth         ! <-- maximum depth and mask to water

! Setup slopping bottom                                ! slope dh/dy

     &     *exp(-1.25*(yr(i,j)/Size_ETA +0.5))

c     &     *exp(1.5*(yr(i,j)/Size_ETA  -0.5))

#  ifdef MASKING
          rmask(i,j)=1.     ! everywhere, then modify it
#  endif
          rd=sqrt((xr(i,j)-x0)**2+(yr(i,j)-y0)**2)

          if (rd <= rd_inner) then
#  ifdef MASKING
            rmask(i,j)=0.
#  endif
#  ifdef GAUSSIAN_SLOPE_ISLAND
            h(i,j)=20.
          else
            h(i,j)=depth - (depth-20.)*exp( -((rd-rd_inner)
     &                                   /(2.*rd_inner))**2 )
#  endif
          endif
        enddo
      enddo
      if (SOUTHERN_EDGE) then
        do i=istrR,iendR
          rmask(i,jstr-1)=1.
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istrR,iendR
          rmask(i,jend+1)=1.
        enddo
      endif

# elif defined SHELFRONT
      do j=jstrR,jendR
        do i=istrR,iendR
          cff=yr(i,j)/1000.
          if (cff<50.) then
            h(i,j)=50.+2.*cff
          elseif (cff<60.) then
            h(i,j)=160.+1.5*(cff-50.)**2-0.1*(cff-60.0)**2
          elseif (cff<100.) then
            h(i,j)=310.+30.*(cff-60.)
          elseif (cff<110.) then
            h(i,j)=1660.-1.5*(cff-110.)**2
          else
            h(i,j)=1660.
          endif
        enddo
      enddo


# elif defined UPWELLING
      do j=jstrR,jendR
#  ifdef MPI
        y=dy*dble(j+jSW_corn)
#  else
        y=dy*dble(j)
#  endif
        if (y>Size_ETA/2.) y=Size_ETA-y+dy
        cff=min(depth,84.5+66.526*tanh(0.00015*(y-0.125*Size_ETA)))
        do i=istrR,iendR
          h(i,j)=cff
        enddo
      enddo

# elif defined RIVER
      do i=istrR,iendR
#  ifdef MPI
        cff=(dble(i +iSW_corn)-0.5D0)/dble(LLm)
#  else
        cff=(dble(i          )-0.5D0)/dble(LLm)
#  endif
        if (cff<0.05) then
          h(i,jstrR)=15.
        elseif (cff<0.15) then
          h(i,jstrR)=15.+843.75*(cff-0.05)**2
        elseif (cff<0.85) then
          h(i,jstrR)=15.+168.75*(cff-0.1)
        elseif (cff<0.95) then
          h(i,jstrR)=150.-843.75*(cff-0.95)**2
        else
          h(i,jstrR)=150.
        endif
      enddo
      do j=jstrR+1,jendR
        do i=istrR,iendR
          h(i,j)=h(i,jstrR)
        enddo
      enddo
#  ifdef MASKING
      do j=jstrR,jendR                ! Set mask to all-water status,
        do i=istrR,iendR              ! then mask out 5-point wide
          rmask(i,j)=1.               ! strip of land on the west, and,
        enddo                         ! finally carve a 1-point wide
      enddo                           ! channel through that strip.
#   ifdef MPI
      imax=min(LLm/15-iSW_corn,iendR) ! translate "absolute" index
      jchn=3*MMm/4 -jSW_corn          ! into index within the
#   else
      imax=min(LLm/15,iendR)          ! presend MPI subdomain.
      jchn=3*MMm/4
#   endif
      do j=jstrR,jendR
        do i=istrR,imax
          rmask(i,j)=0.               ! <-- strip of land
        enddo
      enddo
      do j=jchn-1,jchn+1
        if (j>=jstrR .and. j<=jendR) then
          do i=istrR,imax
            rmask(i,j)=1.             ! <-- channel
          enddo
        endif
      enddo
#  endif
# else
      do j=jstrR,jendR           ! Default version: Setup flat-bottom
        do i=istrR,iendR         ! topography of prescribed depth.
          h(i,j)=depth           ! This is used by some examples:
        enddo                    ! Soliton and Wave Radiation test.
      enddo
# endif
      end
#else
      subroutine ana_grid_empty
      end
#endif /* ANA_GRID */
