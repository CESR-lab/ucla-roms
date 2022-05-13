      module grid

#include "cppdefs.opt"

      use param
      implicit none

! This is include file "grid": Environmental two-dimensional arrays
! associated with curvilinear horizontal coordinate system.
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

#ifdef SPHERICAL
      real,allocatable,dimension(:,:) :: latr
      real,allocatable,dimension(:,:) :: lonr
#else
      real,allocatable,dimension(:,:) :: xp
      real,allocatable,dimension(:,:) :: xr
      real,allocatable,dimension(:,:) :: yp
      real,allocatable,dimension(:,:) :: yr
#endif

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
                                  ! Sometimes it is convenient for post
      subroutine wrt_ana_grid     ! processing to save analytical grid  ![
      use param
      use ncvars
      use mpi
      use netcdf
      implicit none               ! into a separate netCDF file.

# include "netcdf.inc"
      character(len=32) fname
      integer ncid,ierr, lncn, lenstr, my_nf_def_dim

# if defined MPI && !defined PARALLEL_FILES
      integer status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv(blank, 1, MPI_INTEGER,
     &        mynode-1, 1, ocean_grid_comm, status, ierr)
# endif

      fname='grid.nc' ; lncn=lenstr(fname)
# ifdef MPI
#  ifdef PARALLEL_FILES
      ierr=0
      call insert_node(fname, lncn, mynode, NNODES, ierr)
#  else
      if (mynode == 0) then
#  endif
# endif
        ierr=nf90_create(fname, NF90_NETCDF4, ncid)
        if (ierr == nf90_noerr) then
          call put_global_atts (ncid, ierr)

          call def_grid(ncid)
          ierr=nf90_enddef(ncid)
        else
          write(*,'(/1x,4A/12x,A/)')     '### ERROR: wrt_ana_grid :: ',
     &            'Cannot create netCDF file ''', fname(1:lncn), ''',',
     &                                              nf90_strerror(ierr)
        endif
# if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf90_open(fname, nf_write, ncid)
        if (ierr/=nf90_noerr) then
          write(*,'(/1x,4A/12x,A/)')    '### ERROR: wrt_ana_grid :: ',
     &             'Cannot open netCDF file ''', fname(1:lncn), ''',',
     &                                             nf90_strerror(ierr)
        endif
      endif
# endif

      if (ierr == nf90_noerr) then
        call wrt_grid(ncid, fname, lncn)
        ierr=nf90_close(ncid)
        if (ierr /= nf90_noerr) then
           write(*,'(/1x,4A/12x,A/)')   '### ERROR: wrt_ana_grid :: ',
     &            'Cannot close netCDF file ''', fname(1:lncn), ''',',
     &                                             nf90_strerror(ierr)
        endif
      endif

# if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, ocean_grid_comm, ierr)
      endif
# endif
      end subroutine wrt_ana_grid  !]
#endif /* ANA_GRID */

!----------------------------------------------------------------------
      subroutine def_grid(ncid)  ![

! Define grid variables in output NetCDF file, which may be
! restart, history, averages, etc...
!
! Arguments: ncid    NetCDF unit-ID of NetCDF file, which must be
!                            already opened and in definition mode;
!            r2dgrd  integer array of size 2, which contains NetCDF
!                        IDs for dimensions 'xi_rho' and 'eta_rho'.

      use param
      use ncvars
      use netcdf
      use nc_read_write
      implicit none

      ! input
      integer :: ncid
      ! local
      integer :: nf_ftype, varid, ierr

      if (ncid == ncrst) then                              ! decide about precision:
        nf_ftype=NF_FTYPE
      else
        nf_ftype=NF_FOUT
      endif

! Grid type switch: Spherical or Cartesian.
      ierr=nf90_def_var(ncid, 'spherical', nf90_char, varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'grid type logical switch')     ! need to be new line <72 else big space in .nc
      ierr=nf90_put_att(ncid, varid, 'option_T', 'spherical')
      ierr=nf90_put_att(ncid, varid, 'option_F', 'cartesian')

! Physical dimensions of model domain, xl,el (Cartesian grid only).
#ifndef SPHERICAL
      ierr=nf90_def_var(ncid, 'xl', nf_ftype, varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                           'domain length in the XI-direction')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_def_var(ncid, 'el', nf_ftype, varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                           'domain length in the ETA-direction')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
#endif

! Bathymetry.
      varid = nccreate(ncid,'h',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'bathymetry at RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'bath, scalar')

! Coriolis Parameter.
      varid = nccreate(ncid,'f',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'Coriolis parameter at RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'second-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'coriolis, scalar')

! Curvilinear coordinate metric coefficients pm,pn.
      varid = nccreate(ncid,'pm',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                       'curvilinear coordinate metric in XI')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'pm, scalar')
      varid = nccreate(ncid,'pn',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                       'curvilinear coordinate metric in ETA')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'pn, scalar')

! Longitude-latitude or cartezian coordinates of RHO-points.
#ifdef SPHERICAL
      varid = nccreate(ncid,'lon_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'longitude of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'degree_east')
      ierr=nf90_put_att(ncid, varid, 'field', 'lon_rho, scalar')
      varid = nccreate(ncid,'lat_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid,varid,'long_name',
     &                               'latitude of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'degree_north')
      ierr=nf90_put_att(ncid, varid, 'field', 'lat_rho, scalar')
#else
      varid = nccreate(ncid,'x_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'x-locations of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'x_rho, scalar')
      varid = nccreate(ncid,'y_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'y-locations of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'y_rho, scalar')
#endif
#ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points
      varid = nccreate(ncid,'angle',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'angle between XI-axis and EAST')
      ierr=nf90_put_att(ncid, varid, 'units', 'radians')
      ierr=nf90_put_att(ncid, varid, 'field', 'angle, scalar')
#endif
#ifdef MASKING
! Land-Sea mask at RHO-points.
      varid = nccreate(ncid,'mask_rho',(/'xi_rho','eta_rho'/),(/xi_rho,eta_rho/), nf_ftype)
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
      use ncvars
      use scalars
      use netcdf
      use nc_read_write

      implicit none
      integer ncid, varid, ierr, ncdf_write

!#include "netcdf.inc"

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
          ierr=nf90_put_var(ncid, varid, 1, xl)
          if (ierr /= nf90_noerr) then
            call handle_ierr(ierr, 'wrt_grid ')
          endif
        else
          call handle_ierr(ierr, 'wrt_grid ')
        endif
        if (ierr /= nf90_noerr)  goto 99                   !--> ERROR

        ierr=nf90_inq_varid(ncid, 'el', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_put_var(ncid, varid, 1, el)
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

! NOTE: In the case when masking of output data is activated by
! CPP-switch MASK_LAND_DATA, ncdf_write uses sign of its argument
! "varid" as switch to perform on not perform masking.  Because
! none of the grid variables should actually be masked, in the
! code below varid is passed into ncdf_write with negative sign.

! Bathymetry.

      ierr=nf90_inq_varid(ncid, 'h', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0, r_var, h(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR

! Coriolis parameter.

      ierr=nf90_inq_varid(ncid, 'f', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, f(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR

! Curvilinear transformation metrics.

      ierr=nf90_inq_varid(ncid, 'pm', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, pm(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR

      ierr=nf90_inq_varid(ncid, 'pn', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0, r_var, pn(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR

! Longitude-latitude or cartezian coordinates of RHO-points

#ifdef SPHERICAL
      ierr=nf90_inq_varid(ncid, 'lon_rho', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write( ncid, -varid, 0,r_var, lonr(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR

      ierr=nf90_inq_varid(ncid, 'lat_rho', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, latr(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR
#else
      ierr=nf90_inq_varid(ncid, 'x_rho', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, xr(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR

      ierr=nf90_inq_varid(ncid, 'y_rho', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, yr(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR
#endif
#ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points

      ierr=nf90_inq_varid(ncid, 'angle', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, angler(START_2D_ARRAY),1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR

#endif
#ifdef MASKING
! Masking fields at RHO-points.

      ierr=nf90_inq_varid(ncid, 'mask_rho', varid)
      if (ierr == nf90_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, rmask(START_2D_ARRAY), 1)
        if (ierr /= nf90_noerr) then
          call handle_ierr(ierr, 'wrt_grid ')
        endif
      else
        call handle_ierr(ierr, 'wrt_grid ')
      endif
      if (ierr /= nf90_noerr)  goto 99                     !--> ERROR
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

      end module grid
