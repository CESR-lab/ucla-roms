#include "cppdefs.h"
#ifdef ANA_GRID
                                  ! Sometimes it is convenient for post
      subroutine wrt_ana_grid     ! processing to save analytical grid
      use param
      implicit none               ! into a separate netCDF file.

# include "ncvars.h"
# include "netcdf.inc"
      character(len=32) fname
      integer ncid,ierr, r2dgrd(2), lncn, lenstr, my_nf_def_dim

# if defined MPI && !defined PARALLEL_FILES
#  include "mpif.h"
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
        ierr=nf_create(fname, nf_clobber+nf_64bit_offset, ncid)
        if (ierr == nf_noerr) then
          call put_global_atts (ncid, ierr)

          ierr=my_nf_def_dim(ncid, 'xi_rho',  xi_rho,  r2dgrd(1))
          ierr=my_nf_def_dim(ncid, 'eta_rho', eta_rho, r2dgrd(2))

          call def_grid(ncid, r2dgrd)
          ierr=nf_enddef(ncid)
        else
          write(*,'(/1x,4A/12x,A/)')     '### ERROR: wrt_ana_grid :: ',
     &            'Cannot create netCDF file ''', fname(1:lncn), ''',',
     &                                              nf_strerror(ierr)
        endif
# if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf_open(fname, nf_write, ncid)
        if (ierr/=nf_noerr) then
          write(*,'(/1x,4A/12x,A/)')    '### ERROR: wrt_ana_grid :: ',
     &             'Cannot open netCDF file ''', fname(1:lncn), ''',',
     &                                             nf_strerror(ierr)
        endif
      endif
# endif

      if (ierr == nf_noerr) then
        call wrt_grid(ncid, fname, lncn)
        ierr=nf_close(ncid)
        if (ierr /= nf_noerr) then
           write(*,'(/1x,4A/12x,A/)')   '### ERROR: wrt_ana_grid :: ',
     &            'Cannot close netCDF file ''', fname(1:lncn), ''',',
     &                                             nf_strerror(ierr)
        endif
      endif

# if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, ocean_grid_comm, ierr)
      endif
# endif
      end
#endif /* ANA_GRID */



      subroutine def_grid(ncid, r2dgrd)

! Define grid variables in output NetCDF file, which may be
! restart, history, averages, etc...
!
! Arguments: ncid    NetCDF unit-ID of NetCDF file, which must be
!                            already opened and in definition mode;
!            r2dgrd  integer array of size 2, which contains NetCDF
!                        IDs for dimensions 'xi_rho' and 'eta_rho'.

      use param
      implicit none
      integer ncid, r2dgrd(2), nf_ftype, varid, ierr
#include "ncvars.h"
#include "netcdf.inc"

! Decide about precision:

      if (ncid == ncrst) then
        nf_ftype=NF_FTYPE
      else
        nf_ftype=NF_FOUT
      endif

! Grid type switch: Spherical or Cartesian.

      ierr=nf_def_var(ncid, 'spherical', nf_char, 0, 0, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name',24,
     &                                  'grid type logical switch')
      ierr=nf_put_att_text(ncid, varid, 'option_T', 9, 'spherical')
      ierr=nf_put_att_text(ncid, varid, 'option_F', 9, 'cartesian')

! Physical dimensions of model domain, xl,el (Cartesian grid only).

#ifndef SPHERICAL
      ierr=nf_def_var(ncid, 'xl', nf_ftype, 0, 0, varid)
      ierr=nf_put_att_text (ncid, varid, 'long_name',   33,
     &                  'domain length in the XI-direction')
      ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')

      ierr=nf_def_var(ncid, 'el', nf_ftype, 0, 0, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name',   34,
     &                 'domain length in the ETA-direction')
      ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
#endif

! Bathymetry.

      ierr=nf_def_var(ncid, 'h', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name',  24,
     &                           'bathymetry at RHO-points')
      ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
      ierr=nf_put_att_text(ncid, varid, 'field', 12,
     &                                       'bath, scalar')

! Coriolis Parameter.

      ierr=nf_def_var(ncid,'f', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name',   32,
     &                       'Coriolis parameter at RHO-points')
      ierr=nf_put_att_text(ncid, varid, 'units',  8, 'second-1')
      ierr=nf_put_att_text(ncid, varid, 'field', 16,
     &                                       'coriolis, scalar')

! Curvilinear coordinate metric coefficients pm,pn.

      ierr=nf_def_var(ncid, 'pm', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 35,
     &                      'curvilinear coordinate metric in XI')
      ierr=nf_put_att_text(ncid, varid, 'units',  7,    'meter-1')
      ierr=nf_put_att_text(ncid, varid, 'field', 10, 'pm, scalar')

      ierr=nf_def_var(ncid, 'pn', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 36,
     &                     'curvilinear coordinate metric in ETA')
      ierr=nf_put_att_text(ncid, varid, 'units',  7,    'meter-1')
      ierr=nf_put_att_text(ncid, varid, 'field', 10, 'pn, scalar')

! Longitude-latitude or cartezian coordinates of RHO-points.

#ifdef SPHERICAL
      ierr=nf_def_var(ncid, 'lon_rho', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 23,
     &                                    'longitude of RHO-points')
      ierr=nf_put_att_text(ncid, varid, 'units', 11,  'degree_east')
      ierr=nf_put_att_text(ncid, varid, 'field', 15,
     &                                            'lon_rho, scalar')

      ierr=nf_def_var(ncid, 'lat_rho', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid,varid,'long_name',22,
     &                                     'latitude of RHO-points')
      ierr=nf_put_att_text(ncid, varid, 'units', 12,
     &                                               'degree_north')
      ierr=nf_put_att_text(ncid, varid, 'field', 15,
     &                                            'lat_rho, scalar')
#else
      ierr=nf_def_var(ncid, 'x_rho', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 25,
     &                                  'x-locations of RHO-points')
      ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
      ierr=nf_put_att_text(ncid, varid, 'field',13, 'x_rho, scalar')

      ierr=nf_def_var (ncid, 'y_rho', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text (ncid, varid, 'long_name',25,
     &                                   'y-locations of RHO-points')
      ierr=nf_put_att_text (ncid, varid, 'units', 5, 'meter')
      ierr=nf_put_att_text (ncid, varid, 'field',13, 'y_rho, scalar')
#endif
#ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points

      ierr=nf_def_var(ncid, 'angle', nf_ftype, 2, r2dgrd,varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name',30,
     &                             'angle between XI-axis and EAST')
      ierr=nf_put_att_text(ncid, varid, 'units', 7, 'radians')
      ierr=nf_put_att_text(ncid, varid, 'field',13, 'angle, scalar')
#endif
#ifdef MASKING

! Land-Sea mask at RHO-points.

      ierr=nf_def_var(ncid, 'mask_rho', nf_ftype, 2, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name',18,
     &                                         'mask on RHO-points')
      ierr=nf_put_att_text(ncid, varid, 'option_0',    4,   'land' )
      ierr=nf_put_att_text(ncid, varid, 'option_1',    5,   'water')
#endif
      end




      subroutine wrt_grid(ncid, ncname, lncn)

! Write grid variables in output NetCDF file, which may be restart,
! history, averages, etc. All variables are assumed to be previously
! defined by def_grid.
!
! Arguments: ncid    netCDF unit-ID of NetCDF file, which must be
!                            already opened and in definition mode;
!            ncname  netCDF file name (used only in error messages)
!            lncn    length of ncname

      use param
      use grid

      implicit none
      character(len=*) ncname
      integer ncid, lncn, varid, ierr, ncdf_write

#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"

#if defined MPI && !defined PARALLEL_FILES
      if (mynode == 0) then
#endif

! Grid type switch: Spherical or Cartesian.

        ierr=nf_inq_varid(ncid, 'spherical', varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_var1_text(ncid, varid, 1,
#ifdef SPHERICAL
     &                                         'T')
#else
     &                                         'F')
#endif
          if (ierr /= nf_noerr) then
            write(*,2) 'spherical', ncname(1:lncn)
          endif
        else
          write(*,1) 'spherical', ncname(1:lncn)
        endif
        if (ierr /= nf_noerr)  goto 99                   !--> ERROR

#ifndef SPHERICAL
! Physical Dimensions of Model Domain, xl,el (Cartesian grid only).

        ierr=nf_inq_varid(ncid, 'xl', varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_var1_FTYPE(ncid, varid, 1, xl)
          if (ierr /= nf_noerr) then
            write(*,2) 'xl', ncname(1:lncn)
          endif
        else
          write(*,1) 'xl', ncname(1:lncn)
        endif
        if (ierr /= nf_noerr)  goto 99                   !--> ERROR

        ierr=nf_inq_varid(ncid, 'el', varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_var1_FTYPE(ncid, varid, 1, el)
          if (ierr /= nf_noerr) then
            write(*,2) 'el', ncname(1:lncn)
          endif
        else
          write(*,1) 'el', ncname(1:lncn)
        endif
        if (ierr /= nf_noerr)  goto 99                   !--> ERROR
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

      ierr=nf_inq_varid(ncid, 'h', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0, r_var, h(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'h', ncname(1:lncn)
        endif
      else
        write(*,1) 'h', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR

! Coriolis parameter.

      ierr=nf_inq_varid(ncid, 'f', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, f(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'f', ncname(1:lncn)
        endif
      else
        write(*,1) 'f', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR

! Curvilinear transformation metrics.

      ierr=nf_inq_varid(ncid, 'pm', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, pm(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'pm', ncname(1:lncn)
        endif
      else
        write(*,1) 'pm', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR

      ierr=nf_inq_varid(ncid, 'pn', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0, r_var, pn(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'pn', ncname(1:lncn)
        endif
      else
        write(*,1) 'pn', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR

! Longitude-latitude or cartezian coordinates of RHO-points

#ifdef SPHERICAL
      ierr=nf_inq_varid(ncid, 'lon_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write( ncid, -varid, 0,r_var, lonr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'lon_rho', ncname(1:lncn)
        endif
      else
        write(*,1) 'lon_rho', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR

      ierr=nf_inq_varid(ncid, 'lat_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, latr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'lat_rho', ncname(1:lncn)
        endif
      else
        write(*,1) 'lat_rho', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR
#else
      ierr=nf_inq_varid(ncid, 'x_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, xr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'x_rho', ncname(1:lncn)
        endif
      else
        write(*,1) 'x_rho', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR

      ierr=nf_inq_varid(ncid, 'y_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, yr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'y_rho', ncname(1:lncn)
        endif
      else
        write(*,1) 'y_rho', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR
#endif
#ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points

      ierr=nf_inq_varid(ncid, 'angle', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, angler(START_2D_ARRAY),1)
        if (ierr /= nf_noerr) then
          write(*,2) 'angle', ncname(1:lncn)
        endif
      else
        write(*,1) 'angle', ncname(1:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR

#endif
#ifdef MASKING
! Masking fields at RHO-points.

      ierr=nf_inq_varid(ncid, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncid, -varid, 0,r_var, rmask(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'mask_rho', ncname(1:lncn)
        endif
      else
        write(*,1) 'mask_rho', ncname(2:lncn)
      endif
      if (ierr /= nf_noerr)  goto 99                     !--> ERROR
#endif
#ifdef MPI_SILENT_MODE
      if (mynode == 0) then
#endif
        write(*,'(6x,4A,2x,A,I4)')  'wrt_grid :: wrote grid data ',
     &          'into file ''', ncname(1:lncn), '''.' MYID
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
      end
