! $Id: def_floats.F,v 2.1 2004/04/09 21:23:33 capet Exp $

#include "cppdefs.h"
#ifdef FLOATS
                                 ! Create float data NetCDF file:
                                 ! Define its variables, dimensions,
                                 ! and attributes
      subroutine def_floats(ncid, total_rec, ierr,res)

          implicit none
# include "param.h"
# include "mixing.h"
# include "ncvars.h"
# include "ncvars_floats.h"
# include "floats.h"
# include "scalars.h"
# include "strings.h"
# include "netcdf.inc"

      logical create_new_file, res
      integer ncid, total_rec, ierr
      integer lstr, lvar, fltdim, ftimedim, twodim
      integer pgrd(2), temp(2)
      integer lenstr, my_nf_def_dim
# ifdef SOLVE3D
      integer trcdim
# endif /* SOLVE3D */
      character*65 vinfo(4)
      if (may_day_flag/=0) return      !-->  EXIT

      ierr=0
      lstr=lenstr(fltname)
      if (nrpfflt>0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfflt))
        call insert_time_index (fltname, lstr, lvar, ierr)
        if (ierr /= 0) goto 99
      endif

!  Create a new float data file.
!---------------------------------

      create_new_file=ldefflt
      if (ncid/=-1) create_new_file=.false.

#if defined MPI & !defined PARALLEL_FILES
      if (mynode>0) create_new_file=.false.
#endif
 10   if (create_new_file) then
        lstr=lenstr(fltname)
        ierr=nf_create(fltname(1:lstr),nf_clobber,ncid)
        if (ierr==nf_noerr) then
          call put_global_atts (ncid, ierr)
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_floats :: ',
     &                      'Cannot create floats netCDF file ''',
     &                  fltname(1:lstr), '''.', nf_strerror(ierr)
          may_day_flag=3
          return                         !-->  EXIT
        endif

!  Define the dimensions of staggered fields.
!--------------------------------------------

        ierr=my_nf_def_dim(ncid,'drifter',nfloats,fltdim)
# ifdef SOLVE3D
        ierr=my_nf_def_dim(ncid,'tracer',NT,trcdim)
# endif
        ierr=my_nf_def_dim(ncid,'ftime',nf_unlimited,ftimedim)
        ierr=my_nf_def_dim(ncid,'two',2,twodim)

!  Define dimension vectors for point variables.

        pgrd(1)=fltdim
        pgrd(2)=ftimedim

        temp(1)=twodim
        temp(2)=ftimedim

!  Define running parameters.
!----------------------------

!  Time stepping parameters.



!  Define variables and their attributes.
!-----------------------------------------------------------

!  Define time step and model time.

        vinfo(1)='time_step'
        vinfo(2)='time step and record numbers from initialization'
        vinfo(3)='nondimensionnal'
        vinfo(4)='time step/record number, vector, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),nf_int,
     &                    2,temp,fltTstep)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltTstep,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltTstep,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltTstep,'field',lvar,
     &                         vinfo(4)(1:lvar))



! Time

        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var(ncid,vname(1,indxTime)(1:lvar),NF_FTYPE,
     &                    1,ftimedim,fltTime)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text(ncid,fltTime,'long_name',lvar,
     &                         vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text(ncid,fltTime,'units',lvar,
     &                         vname(3,indxTime)(1:lvar))
!
!  Define floats (lon,lat) or (x,y) locations.
!

# ifdef SPHERICAL

        vinfo(1)='lon'
        vinfo(2)='longitude of floats trajectories'
        vinfo(3)='degree_east'
        vinfo(4)='lon, scalar, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltLon)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltLon,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltLon,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltLon,'field',lvar,
     &                         vinfo(4)(1:lvar))

        vinfo(1)='lat'
        vinfo(2)='latitude of floats trajectories'
        vinfo(3)='degree_north'
        vinfo(4)='lat, scalar, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltLat)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltLat,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltLat,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltLat,'field',lvar,
     &                         vinfo(4)(1:lvar))


# else

        vinfo(1)='x'
        vinfo(2)='x-location of floats trajectories'
        vinfo(3)='meter'
        vinfo(4)='x, scalar, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltX)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltX,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltX,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltX,'field',lvar,
     &                         vinfo(4)(1:lvar))


        vinfo(1)='y'
        vinfo(2)='y-location of floats trajectories'
        vinfo(3)='meter'
        vinfo(4)='y, scalar, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltY)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltY,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltY,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltY,'field',lvar,
     &                         vinfo(4)(1:lvar))

# endif  /* SPHERICAL */
        if (wrtflt(indxfltGrd)) then
!  Define grid level

        vinfo(1)='grid_level'
        vinfo(2)='grid level in nested grid hierarchy'
        vinfo(3)='nondimensionnal'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),nf_int,
     &                    2,pgrd,fltGlevel)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltGlevel,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltGlevel,'units',lvar,
     &                         vinfo(3)(1:lvar))


! Define float X-position in the grid

        vinfo(1)='Xgrid'
        vinfo(2)='x-grid floats locations'
        vinfo(3)='nondimensional'
        vinfo(4)='Xgrid, scalar, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltXgrd)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltXgrd,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltXgrd,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltXgrd,'field',lvar,
     &                         vinfo(4)(1:lvar))



! Define float Y-position in the grid

        vinfo(1)='Ygrid'
        vinfo(2)='y-grid floats locations'
        vinfo(3)='nondimensional'
        vinfo(4)='Ygrid, scalar, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltYgrd)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltYgrd,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltYgrd,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltYgrd,'field',lvar,
     &                         vinfo(4)(1:lvar))

#  ifdef SOLVE3D
!
! Define float Z-position in the grid
!

        vinfo(1)='Zgrid'
        vinfo(2)='z-grid floats locations'
        vinfo(3)='nondimensional'
        vinfo(4)='Zgrid, scalar, series'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltZgrd)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltZgrd,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltZgrd,'units',lvar,
     &                         vinfo(3)(1:lvar))
        lvar=lenstr(vinfo(4))
        ierr=nf_put_att_text(ncid,fltZgrd,'field',lvar,
     &                         vinfo(4)(1:lvar))
#  endif /* SOLVE3D */

      endif  ! wrtflt(indxfltGrd)
!!
# ifdef SOLVE3D
!  Define depth.


        vinfo(1)='depth'
        vinfo(2)='depth of floats trajectories'
        vinfo(3)='meter'

        lvar=lenstr(vinfo(1))
        ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltDepth)
        lvar=lenstr(vinfo(2))
        ierr=nf_put_att_text(ncid,fltDepth,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
        lvar=lenstr(vinfo(3))
        ierr=nf_put_att_text(ncid,fltDepth,'units',lvar,
     &                         vinfo(3)(1:lvar))

        if (wrtflt(indxfltTemp)) then
!  Define temperature

          vinfo(1)='temp'
          vinfo(2)='temperature'
          vinfo(3)='degrees Celsius'

          lvar=lenstr(vinfo(1))
          ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltTemp)
          lvar=lenstr(vinfo(2))
          ierr=nf_put_att_text(ncid,fltTemp,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
          lvar=lenstr(vinfo(3))
          ierr=nf_put_att_text(ncid,fltTemp,'units',lvar,
     &                         vinfo(3)(1:lvar))
        endif

#  ifdef SALINITY
!  Salinity

        if (wrtflt(indxfltSalt)) then
          vinfo(1)='salt'
          vinfo(2)='salinity'
          vinfo(3)='PSU'

          lvar=lenstr(vinfo(1))
          ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltSal)
          lvar=lenstr(vinfo(2))
          ierr=nf_put_att_text(ncid,fltSal,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
          lvar=lenstr(vinfo(3))
          ierr=nf_put_att_text(ncid,fltSal,'units',lvar,
     &                         vinfo(3)(1:lvar))
        endif
#  endif
!  Density anomaly

        if (wrtflt(indxfltRho)) then
          vinfo(1)='rho'
          vinfo(2)='density anomaly'
          vinfo(3)='kilogram meter-3'

          lvar=lenstr(vinfo(1))
          ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltDen)
          lvar=lenstr(vinfo(2))
          ierr=nf_put_att_text(ncid,fltDen,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
          lvar=lenstr(vinfo(3))
          ierr=nf_put_att_text(ncid,fltDen,'units',lvar,
     &                         vinfo(3)(1:lvar))

        endif

# endif  /* SOLVE3D */
! Mean velocity (module)

        if (wrtflt(indxfltVel)) then
          vinfo(1)='vel'
          vinfo(2)='mean module velocity'
          vinfo(3)='meter/s'

          lvar=lenstr(vinfo(1))
          ierr=nf_def_var(ncid,vinfo(1)(1:lvar),NF_FTYPE,
     &                    2,pgrd,fltVel)
          lvar=lenstr(vinfo(2))
          ierr=nf_put_att_text(ncid,fltVel,'long_name',lvar,
     &                         vinfo(2)(1:lvar))
          lvar=lenstr(vinfo(3))
          ierr=nf_put_att_text(ncid,fltVel,'units',lvar,
     &                         vinfo(3)(1:lvar))
        endif


!  Leave definition mode.
!---------------------------

        ierr=nf_enddef(ncid)
        res=.true. ! marker that a file has been created
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
          write(*,'(4x,4A,1x,A,i4)')  'def_floats :: created ',
     &            'new file ''', fltname(1:lstr), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif




!  Open an existing float file, check its contents,
!  and prepare for appending data.
!===================================================
!  Inquire about the contents of stations NetCDF file:
!  Inquire about the dimensions and variables. Check for consistency.
!--------------------------------------------------------------------

      elseif (ncid==-1) then

c       write(*,*) 'I TRY TO OPEN THE FILE (XA)'

        ierr=nf_open (fltname(1:lstr), nf_write, ncid)

! for history file, there is a call to checkdims at this point
! think to add one for floats. xa

        if (ierr. ne. nf_noerr) then
#if defined MPI & !defined PARALLEL_FILES
          if (mynode==0) then
            create_new_file=.true.
            goto 10
          else
            write(*,'(/1x,4A,2x,A,I4/)') '### ERROR: def_floats :: ',
     &                  'Cannot open file ''', fltname(1:lstr), '''.'
     &                   MYID
            goto 99                                       !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif



!  Scan variable list from input NetCDF using switches for
!  floats variables. Get variable IDs.

! Time step indices:
!
        ierr=nf_inq_varid (ncid, 'time_step', fltTstep)
        if (ierr /= nf_noerr) then
          write(*,1) 'time_step', fltname(1:lstr)
          goto 99                                         !--> ERROR
        endif


! Time.

        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),fltTime)
        if (ierr /= nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), fltname(1:lstr)
          goto 99                                         !--> ERROR
        endif


! Floats (lon,lat) or (x,y) locations.

# ifdef SPHERICAL

        vinfo(1)='lon'
        lvar=lenstr(vinfo(1))
        ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltLon)
        if (ierr /= nf_noerr) then
          write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        vinfo(1)='lat'
        lvar=lenstr(vinfo(1))
        ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltLat)
        if (ierr /= nf_noerr) then
          write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# else
        vinfo(1)='x'
        lvar=lenstr(vinfo(1))
        ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltX)
        if (ierr /= nf_noerr) then
          write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        vinfo(1)='y'
        lvar=lenstr(vinfo(1))
        ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltY)
        if (ierr /= nf_noerr) then
          write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
          goto 99                                         !--> ERROR
        endif

# endif  /* SPHERICAL */

        if (wrtflt(indxfltGrd)) then

!  Grid level

          vinfo(1)='grid_level'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltGlevel)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                       !--> ERROR
          endif


! float X-position in the grid

          vinfo(1)='Xgrid'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltXgrd)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                         !--> ERROR
          endif


! float Y-position in the grid

          vinfo(1)='Ygrid'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltYgrd)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                         !--> ERROR
          endif


#  ifdef SOLVE3D

! float Z-position in the grid

          vinfo(1)='Zgrid'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltZgrd)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                         !--> ERROR
          endif
#  endif
        endif  !wrtflt(indxfltGrd)

# ifdef SOLVE3D


! Depth

        vinfo(1)='depth'
        lvar=lenstr(vinfo(1))
        ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltDepth)
        if (ierr /= nf_noerr) then
          write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
          goto 99                                         !--> ERROR
        endif

! Temperature

        if (wrtflt(indxfltTemp)) then
          vinfo(1)='temp'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltTemp)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                         !--> ERROR
          endif
        endif

#  ifdef SALINITY
!  Salinity
        if (wrtflt(indxfltSalt)) then
          vinfo(1)='salt'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltSal)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                         !--> ERROR
          endif
         endif
#  endif
!  Density anomaly.

        if (wrtflt(indxfltRho)) then
          vinfo(1)='rho'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltDen)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                         !--> ERROR
          endif
         endif

# endif /* SOLVE3D */

! Mean velocity (module)

        if (wrtflt(indxfltVel)) then
          vinfo(1)='vel'
          lvar=lenstr(vinfo(1))
          ierr=nf_inq_varid (ncid,vinfo(1)(1:lvar),fltVel)
          if (ierr /= nf_noerr) then
            write(*,1) vinfo(1)(1:lvar), fltname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif

!  Set unlimited time record dimension to current value.

        ierr=nf_inq_dimid(ncid,'ftime',ftimedim)
        ierr=nf_inq_dimlen(ncid,ftimedim,nrecflt)
        nrecflt=nrecflt+1
c       write(*,*)'nrecflt = ',nrecflt
        if (ierr /= nf_noerr) then
         write(*,*) '### ERROR: def_floats :: Cannot determine nrecflt.'
        endif

         write(*,'(6x,2A,i4,1x,A,i4)') 'def_floats :: Opened ',
     &                  'existing file from record = ', nrecflt

      endif  ! create or open

   1  format(/1x,'### ERROR: def_floats :: Cannot find variable ''',
     &                        A, ''' in netCDF file ''', A, '''.'/)

 99   return
      end
#else
      subroutine def_floats_empty
      end
#endif /* FLOATS */
