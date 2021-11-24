      module analytical

      ! all analytical routines contained in this module

#include "cppdefs.opt"

      use param
      use ocean2d
      use ocean3d
      use grid
      use scalars
      use tracers

      implicit none
      private

! ANALYTICAL PACKAGE:   Set up various analytical forcing fields
!=========== ========   for the model.
!
! ana_bsedim_tile   Analytical bottom sediment grain size
!                          and density.
! ana_meanRHO_tile  Analytical mean density anomaly.
! ana_smflux_tile   Analytical kinematic surface momentum flux
!                          (wind stress).
! ana_srflux_tile   Analytical kinematic surface shortwave
!                          radiation.
! ana_ssh_tile      Analytical sea surface height climatology.
! ana_sst_tile      Analytical sea surface temperature and dQdSST
!                         which are used during heat flux correction.
! ana_stflux_tile   Analytical kinematic surface flux of tracer type
!                          variables.
! ana_tclima_tile   Analytical tracer climatology fields.
! ana_uclima_tile   Analytical tracer climatology fields.
! ana_wwave_tile    Analytical wind induced wave amplitude,
!                         direction and period.

      ! user input:
      logical, parameter :: restore = .true.

! REMOVE ALL FLAGS LATER
#ifdef ANA_GRID
      public ana_grid
#endif
#ifdef ANA_INITIAL
      public ana_init
#endif
#ifdef ANA_VMIX
      public ana_vmix
#endif

#ifdef ANA_SMFLUX
      public ana_smflux_tile
#endif
#if defined ANA_SRFLUX && defined SOLVE3D
      public ana_srflux_tile
#endif
#if defined ANA_SSH
      public ana_ssh_tile
#endif
#ifdef SOLVE3D
# if defined ANA_STFLUX || defined ANA_SSFLUX
      public ana_stflux_tile
# endif
#endif

      ! Non-typical analytical vars:


      contains

! ----------------------------------------------------------------------
#ifdef ANA_GRID
      subroutine ana_grid (tile)

      use param
      use hidden_mpi_vars

      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call ana_grid_tile (istr,iend,jstr,jend)

      end subroutine ana_grid

! ----------------------------------------------------------------------
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

      use river_frc
      implicit none
      integer istr,iend,jstr,jend, i,j

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

#if defined RIVER_SOURCE
      riv_cells = nint( (riv_east - riv_west)/dx) !number of cells in this river
      do j=jstrR,jendR
        do i=istrR,iendR
          if (xr(i,j)>riv_west .and. xr(i,j)<riv_east) then
            ! find 'coastline' masked cells
            if (rmask(i,j)==0 .and. rmask(i,j+1)==1) then
              rflx(i,j) = 1.0+1.0/riv_cells
            endif
          endif
        enddo
      enddo
#endif

      end subroutine ana_grid_tile
#endif /* ANA_GRID */

! ----------------------------------------------------------------------
#ifdef ANA_INITIAL
      subroutine ana_init (tile)  ! Set initial conditions for momentum
                                  ! free surface, and tracer variables
      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none               ! using analytical expressions.
      integer tile

# include "compute_tile_bounds.h"

      call ana_init_tile( istr,iend,jstr,jend )

      end subroutine ana_init

! ----------------------------------------------------------------------
      subroutine ana_init_tile (istr,iend,jstr,jend)

      use param
      use tracers, only: t
      use boundary
      use climat
      use eos_vars
      use grid
      use mixing
      use ocean2d
      use ocean3d
      use scalars

      implicit none

      integer istr,iend,jstr,jend, i,j,k, itrc, i0,j0

      real Ampl, fac, x,y, x0,y0, cff, cff0,cff1,cff2,cff3
# if defined COLD_FILAMENT

      real, parameter :: f0=0.8E-4,              ! Coriolis parameter
     &     h_0=60., dh=100.,  wdth=5.e+3,
     &     b_0=7.3e-3, N_0=6.3e-3, B=0.02, lam=1./50.

      real xx,zz, hbl, sgm
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: hq

#  define VISC_TERM_WIND_BALANCE
#  ifdef VISC_TERM_WIND_BALANCE
      real  b11,b12,b21,b22, det, TauX,TauY
      real, dimension(0:N) :: Av,FC, c11,c12,c21,c22, d1,d2
#  endif
# elif defined ISWAKE
      integer jwake,mid_wake
      real AmpR,AmpU, cffS, Zthrm, dZthrm, xitm, alpha, csa,sna
# endif

# include "compute_auxiliary_bounds.h"

      do j=jstrR,jendR         ! Set everything (except temperature
        do i=istrR,iendR       ! and salinity) to all-zero state, then
          zeta(i,j,1)=0.       ! modify some of the variables, if a
          ubar(i,j,1)=0.       ! non-trivial initialization required.
          vbar(i,j,1)=0.       ! Note: A code to initialize T [and S]
        enddo                  ! must always be supplied for 3D
      enddo                    ! applications.
# ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            u(i,j,k,1)=0.
            u(i,j,k,2)=0.
            v(i,j,k,1)=0.
            v(i,j,k,2)=0.
          enddo
        enddo
      enddo
# endif

# if defined RIVER_SOURCE
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
!           t(i,j,k,1,itemp)=24.0
            t(i,j,k,1,itemp)=4.+10.*exp(z_r(i,j,k)/50.)
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif

# else
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=???
            t(i,j,k,2,itemp)=???
            t(i,j,k,1,isalt)=???
            t(i,j,k,2,isalt)=???
          enddo
        enddo
      enddo
#  endif
# endif  /* end of configuration selection switch */

! The rest of the code is independent of the specific configuration.

! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.


# if defined OBC_WEST  || defined OBC_EAST ||\
     defined OBC_SOUTH || defined OBC_NORTH
#  ifdef ANA_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR                     ! Save perimeter values
          zeta_west(j)=zeta(istr-1,j,1)      ! of the initial fields
          ubar_west(j)=ubar(istr  ,j,1)      ! to be used as external
          vbar_west(j)=vbar(istr-1,j,1)      ! forcing data for side
        enddo                                ! boundaries.
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_EAST
      if (EASTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR
          zeta_east(j)=zeta(iend+1,j,1)
          ubar_east(j)=ubar(iend+1,j,1)
          vbar_east(j)=vbar(iend+1,j,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_east(j,k)=u(iend+1,j,k,1)
            v_east(j,k)=v(iend+1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_east(j,k,itemp)=t(iend+1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif


#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_south(i)=zeta(i,jstr-1,1)
          ubar_south(i)=ubar(i,jstr-1,1)
          vbar_south(i)=vbar(i,jstr  ,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,jstr-1,k,1)
            v_south(i,k)=v(i,jstr  ,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,itemp)=t(i,jstr-1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_north(i)=zeta(i,jend+1,1)
          ubar_north(i)=ubar(i,jend+1,1)
          vbar_north(i)=vbar(i,jend+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,jend+1,k,1)
            v_north(i,k)=v(i,jend+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,itemp)=t(i,jend+1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  else /* alternative to ANA_BRY */
#   ifdef UCLIMATOLOGY
      do j=jstrR,jendR                   ! Save initial data into
        do i=istrR,iendR                 ! climatology for the
          ubclm(i,j)=ubar(i,j,1)         ! subsequent use as inflow
          vbclm(i,j)=vbar(i,j,1)         ! boundary conditions
        enddo
      enddo
#   endif
#   ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#    ifdef UCLIMATOLOGY
            uclm(i,j,k)=u(i,j,k,1)
            vclm(i,j,k)=v(i,j,k,1)
#    endif
#    ifdef TCLIMATOLOGY
            tclm(i,j,k,itemp)=t(i,j,k,1,itemp)
#    endif
          enddo
        enddo
      enddo
#   endif /* SOLVE3D */
#  endif /* ANA_BRY */
# endif /* at least one of OBC_XXXXX */


# ifdef EXCHANGE
      call exchange2d_3_tile(istr,iend,jstr,jend,
     &                       zeta(START_2D_ARRAY,1),
     &                       ubar(START_2D_ARRAY,1),
     &                       vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
      call exchange_2_tile(istr,iend,jstr,jend,
     &                     u(START_2D_ARRAY,1,1), N,
     &                     v(START_2D_ARRAY,1,1), N)
      do itrc=1,NT
        call exchange_tile(istr,iend,jstr,jend,
     &                t(START_2D_ARRAY,1,1,itrc), N)
      enddo
#  endif
#  ifdef ISWAKE
#   ifdef UCLIMATOLOGY
      call exchange2d_2_tile(istr,iend,jstr,jend, ubclm,vbclm)
#   endif
#   ifdef SOLVE3D
#    ifdef UCLIMATOLOGY
      call exchange_2_tile(istr,iend,jstr,jend, uclm,N, vclm,N)
#    endif
#    ifdef TCLIMATOLOGY
      call exchange_tile(istr,iend,jstr,jend,
     &                    tclm(START_2D_ARRAY,1,itemp), N)
#    endif
#   endif
#  elif defined COLD_FILAMENT
      call exchange_tile(istr,iend,jstr,jend, Akv,N+1)

#  endif
# endif  /* EXCHANGE */


# if defined OBC_WEST && (defined NS_PERIODIC || defined MPI)
      if (WESTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_west,1)
        call exch_bry_EW_tile (jstr,jend, ubar_west,1)
        call exch_bry_EW_tile (jstr,jend, vbar_west,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_west,N)
        call exch_bry_EW_tile (jstr,jend, v_west,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_west(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_EAST && (defined NS_PERIODIC || defined MPI)
      if (EASTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_east,1)
        call exch_bry_EW_tile (jstr,jend, ubar_east,1)
        call exch_bry_EW_tile (jstr,jend, vbar_east,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_east,N)
        call exch_bry_EW_tile (jstr,jend, v_east,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_east(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_SOUTH && (defined EW_PERIODIC || defined MPI)
      if (SOUTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_south,1)
        call exch_bry_NS_tile (istr,iend, ubar_south,1)
        call exch_bry_NS_tile (istr,iend, vbar_south,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_south,N)
        call exch_bry_NS_tile (istr,iend, v_south,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_south(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_NORTH && (defined EW_PERIODIC || defined MPI)
      if (NORTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_north,1)
        call exch_bry_NS_tile (istr,iend, ubar_north,1)
        call exch_bry_NS_tile (istr,iend, vbar_north,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_north,N)
        call exch_bry_NS_tile (istr,iend, v_north,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_north(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

      end subroutine ana_init_tile
#endif /* ANA_INITIAL */

! ----------------------------------------------------------------------
#if defined ANA_VMIX && defined SOLVE3D
      subroutine ana_vmix (tile)

      use param
      use hidden_mpi_vars

      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call ana_vmix_tile (istr,iend,jstr,jend)
      end subroutine ana_vmix

! ----------------------------------------------------------------------
      subroutine ana_vmix_tile (istr,iend,jstr,jend)

! Set vertical kinematic vertical viscosity "Akv" [m^2/s] for momentum
! and vertical mixing c efficients "Akt" [m^2/s] for tracer equations
! using analytical expressions.

      use param
      use grid
      use mixing
      use ocean3d
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real cff

# include "compute_auxiliary_bounds.h"

      if (FIRST_TIME_STEP) then
        do k=0,N
          do j=jstr,jend
            do i=istr,iend
             Akv(i,j,k)=0.0

             Akt(i,j,k,itemp)=0.
#  ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
            enddo
          enddo
        enddo
      endif


# ifdef EXCHANGE
      call exchange_2_tile (istr,iend,jstr,jend,  Akv, N+1,
     &                    Akt(START_2D_ARRAY,0,itemp), N+1)
#  ifdef SALINITY
      call exchange_tile (istr,iend,jstr,jend,
     &                    Akt(START_2D_ARRAY,0,isalt), N+1)
#  endif
# endif

      end subroutine ana_vmix_tile
#endif /* ANA_VMIX */

! ----------------------------------------------------------------------

! ONLY THE _TILE ROUTINES ARE USED in set_forces! GET RID OF THE NON-TILE VERSIONS

! ----------------------------------------------------------------------
#ifdef ANA_SMFLUX
      subroutine ana_smflux_tile(istr,iend,jstr,jend)

      use param
      use surf_flux !, only: sustr,svstr
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
      real TauX, TauY, cff

# include "compute_extended_bounds.h"

      if (FIRST_TIME_STEP) then  ! Set wind stress (kinematic surface
        do j=jstrR,jendR         ! momentum flux [m^2/s^2]) components
          do i=istrR,iendR       ! in XI- and ETA-directions [m^2/s^2])
            sustr(i,j)=0.        ! at horizontal U- and V-points.
            svstr(i,j)=0.        ! For code compactness, both are set
          enddo                  ! to zero at first (default values),
        enddo                    ! then either one of them, or both may
      endif                      ! be overwritten by nontrivial values.

      end
#endif /* ANA_SMFLUX */

! ----------------------------------------------------------------------
#if defined ANA_SRFLUX && defined SOLVE3D
      subroutine ana_srflux_tile(istr,iend,jstr,jend)

! Set kinematic surface solar shortwave radiation flux "srflx"
! [degC m/s] using an analytical expression.
      use param
      use surf_flux, only: srflx
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j

# include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          srflx(i,j)=0.
        enddo
      enddo
      end
#endif /* ANA_SRFLUX */

! ----------------------------------------------------------------------
#if defined ANA_SSH
      subroutine ana_ssh(tile)            ! Set analytical sea-surface
      use param
      use hidden_mpi_vars
      implicit none                       ! height climatology [m]
      integer tile
#include "compute_tile_bounds.h"
      call ana_ssh_tile(istr,iend,jstr,jend)
      end

! ----------------------------------------------------------------------
      subroutine ana_ssh_tile(istr,iend,jstr,jend)

      use param
      use climat
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j

# include "compute_extended_bounds.h"

      end
#endif

! ----------------------------------------------------------------------
#ifdef SOLVE3D
# if defined ANA_STFLUX || defined ANA_SSFLUX
      subroutine ana_stflux_tile(istr,iend,jstr,jend, itrc)

! Set kinematic surface flux of tracer type variables "stflx" (tracer
! units m/s) using analytical expressions.
      use param
      use surf_flux, only: stflx
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, itrc, i,j

# include "compute_extended_bounds.h"

      if (itrc==itemp) then

! Set kinematic surface heat flux [degC m/s] at horizontal
! RHO-points.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,itemp)=0.
          enddo
        enddo

#  ifdef SALINITY
      elseif (itrc == isalt) then

! Set kinematic surface freshwater flux (m/s) at horizontal
! RHO-points, scaling by surface salinity is done in STEP3D.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,isalt)=0.
          enddo
        enddo
#  endif
      else

!  Set kinematic surface flux of additional tracers, if any.

      endif
      end
# endif /* ANA_STFLUX || ANA_SSFLUX */

#endif /* SOLVE3D */
! ----------------------------------------------------------------------

      end module analytical
