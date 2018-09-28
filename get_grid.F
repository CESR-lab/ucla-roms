#include "cppdefs.h"
#ifndef ANA_GRID

      subroutine get_grid         ! Read in grid data from NetCDF file.
      implicit none
      character*1 char1
      integer ierr, ncid, varid, lstr, lenstr, checkdims, ncdf_read
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "ncvars.h"
# include "netcdf.inc"

! Open grid netCDF file for reading. Check that dimensions in the file
! are consistent with the model, then read all necessary variables.

      lstr=lenstr(grdname)
      ierr=nf_open(grdname(1:lstr), nf_nowrite, ncid)
      if (ierr == nf_noerr) then
        ierr=checkdims (ncid, grdname, varid)
        if (ierr. ne. nf_noerr) goto 99
      else
        write(*,'(/1x,4A/12x,A)')  '### ERROR: get_grid :: Cannot ',
     &          'open input NetCDF file ''', grdname(1:lstr), '''.',
     &                                            nf_strerror(ierr)
        goto 99                                          !--> ERROR
      endif

! Logical switch for spherical grid configuration:

      ierr=nf_inq_varid (ncid, 'spherical', varid)
      if (ierr == nf_noerr) then
        ierr=nf_get_var1_text (ncid, varid, 1, char1)
        if (ierr /= nf_noerr) then
          write(*,2) 'spherical', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'spherical', grdname(1:lstr)
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

        ierr=nf_inq_varid (ncid, 'xl', varid)
        if (ierr == nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, varid, 1, xl)
          if (ierr /= nf_noerr) then
            write(*,2) 'xl', grdname(1:lstr), nf_strerror(ierr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'xl', grdname(1:lstr)
          goto 99                                        !--> ERROR
        endif

        ierr=nf_inq_varid (ncid, 'el', varid)
        if (ierr == nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, varid, 1, el)
          if (ierr /= nf_noerr) then
            write(*,2) 'el', grdname(1:lstr), nf_strerror(ierr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'el', grdname(1:lstr)
          goto 99                                        !--> ERROR
        endif
      endif

! Read in grid arrays.
!===== == ==== =======
! Bathymetry:

      ierr=nf_inq_varid (ncid, 'h', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, h(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'h', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'h', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif

! Coriolis parameter.

      ierr=nf_inq_varid (ncid, 'f', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, f(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'f', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'f', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif

! Coordinate transfomation metrics (m,n) associated with
! differential distances in XI- and ETA-directions.

      ierr=nf_inq_varid (ncid, 'pm', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, pm(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'pm', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'pm', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'pn', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, pn(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'pn', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'pn', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif

! Coordinates (lon,lat [degrees]) or (x,y [meters]) at RHO-points.

# ifdef SPHERICAL
      ierr=nf_inq_varid (ncid, 'lon_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, lonr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'lon_rho', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'lon_rho', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'lat_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, latr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'lat_rho', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'lat_rho', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif
# else
      ierr=nf_inq_varid (ncid, 'x_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncf_fread(ncid, varid, 0, r_var, xr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'x_rho', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'x_rho', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif

      ierr=nf_inq_varid (ncid, 'y_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, yr(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'y_rho', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'y_rho', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif
# endif
# ifdef CURVGRID

! Angle (radians) between XI-axis and EAST at RHO-points.

      ierr=nf_inq_varid (ncid, 'angle', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, angler(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'angle', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'angle', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif
# endif
# ifdef MASKING

! Mask at RHO-points.

      ierr=nf_inq_varid (ncid, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        ierr=ncdf_read(ncid, varid, 0, r_var, rmask(START_2D_ARRAY), 1)
        if (ierr /= nf_noerr) then
          write(*,2) 'mask_rho', grdname(1:lstr), nf_strerror(ierr)
          goto 99                                        !--> ERROR
        endif
      else
        write(*,1) 'mask_rho', grdname(1:lstr)
        goto 99                                          !--> ERROR
      endif
# endif
      ierr=nf_close(ncid)
      if (ierr==nf_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
          write(*,'(6x,4A,1x,A,I4)') 'get_grid :: read grid data ',
     &            'from file ''',  grdname(1:lstr), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(1x,4A/12x,A,1x,A,I4)')   '### ERROR: get_grid :: ',
     &                'Cannot close file ''', grdname(1:lstr), '''.',
     &                 nf_strerror(ierr) MYID
        goto 99
      endif
      return                                     !--> NORMAL RETURN

   1  format(/1x, '### ERROR: get_grid :: Cannot find variable ''', A,
     &                          ''' within netCDF file ''', A, '''.'/)
   2  format(/1x, '### ERROR: get_grid :: Cannot read variable ''', A,
     &                     ''' from netCDF file ''', A, '''.' /12x,A/)
  99  may_day_flag=2
      return                                             !--> ERROR
      end
#else
      subroutine get_grid_empty
      end
#endif /* !ANA_GRID */

