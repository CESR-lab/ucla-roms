      module grid

#include "cppdefs.opt"

      use param
      use dimensions
      use roms_read_write
      implicit none

!
! h         Model bottom topography (depth [m] at RHO-points.)
! f, fomn   Coriolis parameter [1/s] and compound term f/[pm*pn]
!                                                   at RHO points.
! angler      Angle [radians] between XI-axis and the direction
!                                       to the EAST at RHO-points.
! latr, lonr  Latitude (degrees north) and Longitude (degrees east)
!                                                  at RHO-points.
! xr, xp      XI-coordinates [m] at RHO- and PSI-points.
! yr, yp      ETA-coordinates [m] at RHO- and PSI-points.
!
! pm, pm  Coordinate transformation metric "m" and "n" associated
!         with the differential distances in XI- and ETA-directions.
!
! dm_u, dm_r  Grid spacing [meters] in the XI-direction
! dm_v, dm_p       at U-, RHO-,  V- and vorticity points.
! dn_u, dn_r  Grid spacing [meters] in the ETA-direction
! dn_v, dn_p      at U-, RHO-,  V- and vorticity points.
!
! dmde     ETA-derivative of inverse metric factor "m" d(1/M)/d(ETA)
! dndx     XI-derivative  of inverse metric factor "n" d(1/N)/d(XI)
!
! pmon_u   Compound term, pm/pn at U-points.
! pnom_v   Compound term, pn/pm at V-points.
!
! umask, rmask  Land-sea masking arrays at RHO-,U-,V- and PSI-points
! pmask, vmask      (rmask,umask,vmask) = (0=Land, 1=Sea);
!                    pmask = (0=Land, 1=Sea, 1-gamma2 =boundary).
!
      real,allocatable,dimension(:,:) :: h
      real,allocatable,dimension(:,:) :: hinv
      real,allocatable,dimension(:,:) :: f
      real,allocatable,dimension(:,:) :: fomn
# ifdef NON_TRADITIONAL
      real,allocatable,dimension(:,:) :: f_XI
      real,allocatable,dimension(:,:) :: f_ETA
# endif

# ifdef CURVGRID
      real,allocatable,dimension(:,:) :: angler
# endif

      real,allocatable,dimension(:,:) :: latr
      real,allocatable,dimension(:,:) :: lonr
      real,allocatable,dimension(:,:) :: xp
      real,allocatable,dimension(:,:) :: xr
      real,allocatable,dimension(:,:) :: yp
      real,allocatable,dimension(:,:) :: yr

      real,allocatable,dimension(:,:) :: pm
      real,allocatable,dimension(:,:) :: pn
      real,allocatable,dimension(:,:) :: dm_r
      real,allocatable,dimension(:,:) :: dn_r
      real,allocatable,dimension(:,:) :: pn_u
      real,allocatable,dimension(:,:) :: dm_u
      real,allocatable,dimension(:,:) :: dn_u
      real,allocatable,dimension(:,:) :: dm_v
      real,allocatable,dimension(:,:) :: pm_v
      real,allocatable,dimension(:,:) :: dn_v
      real,allocatable,dimension(:,:) :: dm_p
      real,allocatable,dimension(:,:) :: dn_p

      real,allocatable,dimension(:,:) :: iA_u
      real,allocatable,dimension(:,:) :: iA_v

#if (defined CURVGRID && defined UV_ADV)
      real,allocatable,dimension(:,:) :: dmde
      real,allocatable,dimension(:,:) :: dndx
#endif
      real,allocatable,dimension(:,:) :: pmon_u
      real,allocatable,dimension(:,:) :: pnom_v
      real,allocatable,dimension(:,:) :: grdscl

#ifdef MASKING
      real,allocatable,dimension(:,:) :: rmask
      real,allocatable,dimension(:,:) :: pmask
      real,allocatable,dimension(:,:) :: umask
      real,allocatable,dimension(:,:) :: vmask
#endif

      public :: get_grid

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_grid  ![
      use scalars
      implicit none

      ! WARNING: "rmask" MUST BE initialized to all-one (=1) state in order to
      ! read grid variables (coordinates, metric, topography), which should
      ! not be masked.

#ifdef MASKING
      allocate( rmask(GLOBAL_2D_ARRAY) ) ; rmask=1.
      allocate( pmask(GLOBAL_2D_ARRAY) ) ; pmask=init
      allocate( umask(GLOBAL_2D_ARRAY) ) ; umask=init
      allocate( vmask(GLOBAL_2D_ARRAY) ) ; vmask=init
#endif

      allocate( h(GLOBAL_2D_ARRAY) )    ; h=init           ! potential first touch issue.
      allocate( hinv(GLOBAL_2D_ARRAY) ) ; hinv=init        ! before only rmask was set in init_arrays...
      allocate( f(GLOBAL_2D_ARRAY) )    ; f=init
      allocate( fomn(GLOBAL_2D_ARRAY) ) ; fomn=init
# ifdef NON_TRADITIONAL
      allocate( f_XI(GLOBAL_2D_ARRAY) )  ; f_XI=init
      allocate( f_ETA(GLOBAL_2D_ARRAY) ) ; f_ETA=init
# endif
 
# ifdef CURVGRID
      allocate( angler(GLOBAL_2D_ARRAY) ) ; angler=init
# endif
 
#ifdef SPHERICAL
      allocate( latr(GLOBAL_2D_ARRAY) ) ; latr=init
      allocate( lonr(GLOBAL_2D_ARRAY) ) ; lonr=init
#else
      allocate( xp(GLOBAL_2D_ARRAY) ) ; xp=init
      allocate( xr(GLOBAL_2D_ARRAY) ) ; xr=init
      allocate( yp(GLOBAL_2D_ARRAY) ) ; yp=init
      allocate( yr(GLOBAL_2D_ARRAY) ) ; yr=init
#endif
 
      allocate( pm(GLOBAL_2D_ARRAY) )   ; pm=init
      allocate( pn(GLOBAL_2D_ARRAY) )   ; pn=init
      allocate( dm_r(GLOBAL_2D_ARRAY) ) ; dm_r=init
      allocate( dn_r(GLOBAL_2D_ARRAY) ) ; dn_r=init
      allocate( pn_u(GLOBAL_2D_ARRAY) ) ; pn_u=init
      allocate( dm_u(GLOBAL_2D_ARRAY) ) ; dm_u=init
      allocate( dn_u(GLOBAL_2D_ARRAY) ) ; dn_u=init
      allocate( dm_v(GLOBAL_2D_ARRAY) ) ; dm_v=init
      allocate( pm_v(GLOBAL_2D_ARRAY) ) ; pm_v=init
      allocate( dn_v(GLOBAL_2D_ARRAY) ) ; dn_v=init
      allocate( dm_p(GLOBAL_2D_ARRAY) ) ; dm_p=init
      allocate( dn_p(GLOBAL_2D_ARRAY) ) ; dn_p=init

      allocate( iA_u(GLOBAL_2D_ARRAY) ) ; iA_u=init
      allocate( iA_v(GLOBAL_2D_ARRAY) ) ; iA_v=init
 
#if (defined CURVGRID && defined UV_ADV)
      allocate( dmde(GLOBAL_2D_ARRAY) ) ; dmde=init
      allocate( dndx(GLOBAL_2D_ARRAY) ) ; dndx=init
#endif
      allocate( pmon_u(GLOBAL_2D_ARRAY) ) ; pmon_u=init
      allocate( pnom_v(GLOBAL_2D_ARRAY) ) ; pnom_v=init
      allocate( grdscl(GLOBAL_2D_ARRAY) ) ; grdscl=init

      end subroutine init_arrays_grid  !]
!----------------------------------------------------------------------
#ifdef ANA_GRID
! ----------------------------------------------------------------------
      subroutine ana_grid (tile)  ![

      use param
      use hidden_mpi_vars

      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call ana_grid_tile (istr,iend,jstr,jend)

      call wrt_ana_grid

      end subroutine ana_grid  !]
! ----------------------------------------------------------------------
      subroutine ana_grid_tile (istr,iend,jstr,jend)  ![

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
      use scalars

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
      real psz,px,py,pipe_cells
      real :: riv_east,riv_west

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
      do j=0,ny+1
        do i=0,nx+1

          if(yr(i,j)<shelf) then
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

      do j=jstr-1,jend+1
        do i=istr-1,iend+1
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

!     if (river_source) then
!        rflx is defined in river_frc.F (init_river_frc)
!     endif

      end subroutine ana_grid_tile  !]
! ----------------------------------------------------------------------
      subroutine wrt_ana_grid     ! write analytical grid to file  ![
      use param
      use mpi
      use netcdf
      implicit none               ! into a separate netCDF file.

      character(len=99),save  :: fname
      integer ncid,ierr

      call create_file('_grd',fname,.true.) ! last argument is to skip date label

      ierr=nf90_open(fname,nf90_write,ncid)
      call def_grid(ncid)
      call wrt_grid(ncid)
      ierr=nf90_close(ncid)

      end subroutine wrt_ana_grid  !]
#endif /* ANA_GRID */

!----------------------------------------------------------------------
      subroutine def_grid(ncid)  ![

! Define grid variables in output NetCDF file, which may be
! restart, history, averages, etc...
!
! Arguments: ncid    NetCDF unit-ID of NetCDF file, which must be
!                            already opened and in definition mode;

      use param
      use netcdf
      use nc_read_write
      implicit none

      ! input
      integer :: ncid
      ! local
      integer :: nf_ftype, varid, ierr

! Grid type switch: Spherical or Cartesian.
      ierr=nf90_def_var(ncid, 'spherical', nf90_char, varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'grid type logical switch')     ! need to be new line <72 else big space in .nc
      ierr=nf90_put_att(ncid, varid, 'option_T', 'spherical')
      ierr=nf90_put_att(ncid, varid, 'option_F', 'cartesian')

! Physical dimensions of model domain, xl,el (Cartesian grid only).
#ifndef SPHERICAL
      ierr=nf90_def_var(ncid, 'xl', nf90_double, varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                           'domain length in the XI-direction')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_def_var(ncid, 'el', nf90_double, varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                           'domain length in the ETA-direction')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
#endif

! Bathymetry.
      varid = nccreate(ncid,'h',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'bathymetry at RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'bath, scalar')

! Coriolis Parameter.
      varid = nccreate(ncid,'f',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'Coriolis parameter at RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'second-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'coriolis, scalar')

! Curvilinear coordinate metric coefficients pm,pn.
      varid = nccreate(ncid,'pm',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                       'curvilinear coordinate metric in XI')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'pm, scalar')
      varid = nccreate(ncid,'pn',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                       'curvilinear coordinate metric in ETA')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'pn, scalar')

! Longitude-latitude or cartezian coordinates of RHO-points.
#ifdef SPHERICAL
      varid = nccreate(ncid,'lon_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'longitude of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'degree_east')
      ierr=nf90_put_att(ncid, varid, 'field', 'lon_rho, scalar')
      varid = nccreate(ncid,'lat_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid,varid,'long_name',
     &                               'latitude of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'degree_north')
      ierr=nf90_put_att(ncid, varid, 'field', 'lat_rho, scalar')
#else
      varid = nccreate(ncid,'x_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'x-locations of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'x_rho, scalar')
      varid = nccreate(ncid,'y_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'y-locations of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'y_rho, scalar')
#endif
#ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points
      varid = nccreate(ncid,'angle',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'angle between XI-axis and EAST')
      ierr=nf90_put_att(ncid, varid, 'units', 'radians')
      ierr=nf90_put_att(ncid, varid, 'field', 'angle, scalar')
#endif
#ifdef MASKING
! Land-Sea mask at RHO-points.
      varid = nccreate(ncid,'mask_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'mask on RHO-points')
      ierr=nf90_put_att(ncid, varid, 'option_0', 'land' )
      ierr=nf90_put_att(ncid, varid, 'option_1', 'water')
#endif

      end subroutine def_grid  !]
!----------------------------------------------------------------------
      subroutine wrt_grid(ncid)  ![
      ! Write grid variables in output NetCDF file, which may be restart,
      ! history, averages, etc. All variables are assumed to be previously
      ! defined by def_grid.
      !
      ! Arguments: ncid    netCDF unit-ID of NetCDF file, which must be
      !                            already opened and in definition mode

      use param
      use scalars
      use netcdf
      use nc_read_write

      implicit none
      ! input
      integer,intent(in) :: ncid
      ! local
      integer :: varid, ierr, ncdf_write

#if defined MPI && !defined PARALLEL_FILES
      if (mynode == 0) then
#endif

! Grid type switch: Spherical or Cartesian.

        ierr=nf90_inq_varid(ncid, 'spherical', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_put_var(ncid, varid,
#ifdef SPHERICAL
     &                                         'T')
#else
     &                                         'F')
#endif
          if (ierr /= nf90_noerr) then
            call handle_ierr(ierr, 'wrt_grid spherical')
          endif
        else
          call handle_ierr(ierr, 'wrt_grid spherical')
        endif
        if (ierr /= nf90_noerr)  goto 99                   !--> ERROR

#ifndef SPHERICAL
! Physical Dimensions of Model Domain, xl,el (Cartesian grid only).

        ierr=nf90_inq_varid(ncid, 'xl', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_put_var(ncid, varid, xl)
          if (ierr /= nf90_noerr) then
            call handle_ierr(ierr, 'wrt_grid ')
          endif
        else
          call handle_ierr(ierr, 'wrt_grid ')
        endif
        if (ierr /= nf90_noerr)  goto 99                   !--> ERROR

        ierr=nf90_inq_varid(ncid, 'el', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_put_var(ncid, varid, el)
          if (ierr /= nf90_noerr) then
            call handle_ierr(ierr, 'wrt_grid ')
          endif
        else
          call handle_ierr(ierr, 'wrt_grid ')
        endif
        if (ierr /= nf90_noerr)  goto 99                   !--> ERROR
#endif


#if defined MPI && !defined PARALLEL_FILES
      endif
#endif

      call ncwrite(ncid,'h', h(i0:i1,j0:j1))                    ! Bathymetry.
      call ncwrite(ncid,'f', f(i0:i1,j0:j1))                    ! Coriolis parameter.

      call ncwrite(ncid,'pm', pm(i0:i1,j0:j1))                  ! pm=1./dx, pn=1./dy
      call ncwrite(ncid,'pn', pn(i0:i1,j0:j1))       

#ifdef SPHERICAL
      call ncwrite(ncid,'lon_rho', lonr(i0:i1,j0:j1))           ! longitude-latitude
      call ncwrite(ncid,'lat_rho', latr(i0:i1,j0:j1))           !
#else
      call ncwrite(ncid,'x_rho', xr(i0:i1,j0:j1))               ! Cartesian coordinates
      call ncwrite(ncid,'y_rho', yr(i0:i1,j0:j1))
#endif
#ifdef CURVGRID
      call ncwrite(ncid,'angle', angler(i0:i1,j0:j1))           ! angle between XI-axis and EAST at RHO-points
#endif
#ifdef MASKING
      call ncwrite(ncid,'mask_rho', rmask(i0:i1,j0:j1))         ! masking fields at RHO-points.
#endif
#ifdef MPI_SILENT_MODE
      if (mynode == 0) then
#endif
        write(*,'(6x,A)')  'wrt_grid :: wrote grid data '
#ifdef MPI_SILENT_MODE
      endif
#endif

      return
 99   may_day_flag=3
      return

  1   format(/1x,'### ERROR: wrt_grid :: Cannot find variable ''',
     &             A, ''' within netCDF' /25x, 'file ''', A,'''.'/)
  2   format(/1x,'### ERROR: wrt_grid :: Cannot write variable ''',
     &               A, ''' into netCDF' /25x, 'file ''', A,'''.'/)

      end subroutine wrt_grid  !]
!----------------------------------------------------------------------
      subroutine get_grid  ![
      ! read in grid data from NetCDF file or define analytical grid

      use roms_read_write, only: nc_check_units
      use scalars
      use netcdf
      use nc_read_write, only: ncread
      use dimensions

      implicit none

#ifdef ANA_GRID
      call ana_grid(0)
#else 

      character*1 char1
      integer :: ierr, ncid, varid, lstr, lenstr, checkdims
      integer,dimension(4) :: start

! Open grid netCDF file for reading. Check that dimensions in the file
! are consistent with the model, then read all necessary variables.

!      lstr=lenstr(grdname)
      ierr=nf90_open(grdname, nf90_nowrite, ncid)
      if (ierr == nf90_noerr) then
        ierr=checkdims (ncid, grdname, varid)
        if (ierr. ne. nf90_noerr) goto 99
      else
        write(*,'(/1x,4A/12x,A)')  '### ERROR: get_grid :: Cannot ',
     &          'open input NetCDF file ''', trim(grdname), '''.',
     &                                            nf90_strerror(ierr)
        goto 99                                          !--> ERROR
      endif

! Logical switch for spherical grid configuration:

      ierr=nf90_inq_varid (ncid, 'spherical', varid)
      if (ierr == nf90_noerr) then
        ierr=nf90_get_var(ncid, varid, char1)
        if (ierr /= nf90_noerr) then
          write(*,2) 'spherical', trim(grdname), nf90_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'spherical', trim(grdname)
        goto 99                                          !--> ERROR
      endif

      if (char1=='t' .or. char1=='T') then
# ifdef SPHERICAL
        mpi_master_only write(*,'(/1x,A/)') 'Spherical grid detected.'
# else
        write(*,'(/1x,2A/24x,A/)')   '### ERROR: get_grid :: ',
     &          'Spherical grid file detected, but CPP-switch',
     &                                 'SPHERICAL is not set.'
        goto 99                                          !--> ERROR
# endif
      else

! Physical dimensions of the basin in XI- and ETA-directions:

        ierr=nf90_inq_varid (ncid, 'xl', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_get_var(ncid, varid, xl)
          if (ierr /= nf90_noerr) then
            write(*,2) 'xl', trim(grdname), nf90_strerror(ierr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'xl', trim(grdname)
          goto 99                                        !--> ERROR
        endif

        ierr=nf90_inq_varid (ncid, 'el', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_get_var(ncid, varid, el)
          if (ierr /= nf90_noerr) then
            write(*,2) 'el', trim(grdname), nf90_strerror(ierr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'el', trim(grdname)
          goto 99                                        !--> ERROR
        endif
      endif

! Read in grid arrays.
!===== == ==== =======

      start=1; start(3)=0                                            ! 2D vars

      call ncread(ncid,'h',h(i0:i1,j0:j1),start)                     ! Bathymetry:
      call exchange_tile(1,nx,1,ny, h(START_2D_ARRAY),1)             ! h needed an exchange

      call ncread(ncid,'f',f(i0:i1,j0:j1),start)                     ! Coriolis parameter.
      call exchange_tile(1,nx,1,ny, f(START_2D_ARRAY),1)

! Coordinate transfomation metrics (m,n) associated with
! differential distances in XI- and ETA-directions.

      call ncread(ncid,'pm',pm(i0:i1,j0:j1),start)
      call exchange_tile(1,nx,1,ny, pm(START_2D_ARRAY),1)

      call ncread(ncid,'pn',pn(i0:i1,j0:j1),start)
      call exchange_tile(1,nx,1,ny, pn(START_2D_ARRAY),1)

# ifdef SPHERICAL
      call ncread(ncid,'lon_rho',lonr(i0:i1,j0:j1),start)  ! Coordinates (lon,lat [degrees]) or (x,y [meters]) at RHO-points.
      call exchange_tile(1,nx,1,ny, lonr(START_2D_ARRAY),1)
      call ncread(ncid,'lat_rho',latr(i0:i1,j0:j1),start)
      call exchange_tile(1,nx,1,ny, latr(START_2D_ARRAY),1)
# else
      call ncread(ncid,'x_rho',xr(i0:i1,j0:j1),start)
      call exchange_tile(1,nx,1,ny, xr(START_2D_ARRAY),1)
      call ncread(ncid,'y_rho',yr(i0:i1,j0:j1),start)
      call exchange_tile(1,nx,1,ny, yr(START_2D_ARRAY),1)
# endif
# ifdef CURVGRID
      call ncread(ncid,'angle',angler(i0:i1,j0:j1),start)  ! Angle (radians) between XI-axis and EAST at RHO-points.
      call exchange_tile(1,nx,1,ny, angler(START_2D_ARRAY),1)
# endif
# ifdef MASKING
      call ncread(ncid,'mask_rho',rmask(i0:i1,j0:j1),start)! Mask at RHO-points.
      call exchange_tile(1,nx,1,ny, rmask(START_2D_ARRAY),1)
# endif
      ierr=nf90_close(ncid)
      if (ierr==nf90_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
          write(*,'(6x,4A,1x,A,I4)') 'get_grid :: read grid data ',
     &            'from file ''',  trim(grdname), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(1x,4A/12x,A,1x,A,I4)')   '### ERROR: get_grid :: ',
     &                'Cannot close file ''', trim(grdname), '''.',
     &                 nf90_strerror(ierr) MYID
        goto 99
      endif
      return                                     !--> NORMAL RETURN

   1  format(/1x, '### ERROR: get_grid :: Cannot find variable ''', A,
     &                          ''' within netCDF file ''', A, '''.'/)
   2  format(/1x, '### ERROR: get_grid :: Cannot read variable ''', A,
     &                     ''' from netCDF file ''', A, '''.' /12x,A/)
  99  may_day_flag=2
      return                                             !--> ERROR

#endif /* ifdef ANA_GRID */

      end subroutine get_grid  !]

!----------------------------------------------------------------------

      end module grid
