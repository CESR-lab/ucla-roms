#include "cppdefs.h"
#ifdef STATIONS
c--# define VERBOSE

      subroutine def_statn (ncid, total_rec, ierr) ! Create or open
                                                   ! output netCDF file
      implicit none                                ! file for stations
      integer ncid, total_rec, ierr
# include "param.h"
# include "scoord.h"
# include "grid.h"
# include "mixing.h"
# include "ncvars.h"
# include "netcdf.inc"
# include "scalars.h"
# include "strings.h"

      integer  k, lncn, lvar, stadim, timedim, varid
      integer pgrd(2)
      integer lenstr, my_nf_def_dim
# ifdef SOLVE3D
      integer itrc, srdim, swdim, trcdim
      integer rgrd(3), wgrd(3)
# endif
      real wrk(NS)
      logical create_new_file

# ifdef VERBOSE
      write(*,*) 'entered def_statn, total_rec =', total_rec
# endif
      ierr=0         ! <-- initialise: in may never be touched below
      lncn=lenstr(staname)
      if (nrpfstn>0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfstn))
        call insert_time_index (staname, lncn, lvar, ierr)
        if (ierr /= 0) goto 99
      endif

! Decide whether to create a new file, or open an existing one.
! Overall the code below is globally organized as a 3-way switch,
!
! 10  if (create_new_file) then
!        ... create new file, save netCDF ids for all variables;
!     elseif (ncid==-1) then
!        ... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be
!        .... already known (MPI single file output only).
!     endif
!
! which is designed to implement a flexible opening policy: if
! ldefhis=.true., it forces creation of a new file [hence, if the
! file already exists, it will be overwritten]; conversely, if
! ldefhis=.false. it is assumed that the file already exists and an
! attempt to open it is made; if the attempt is successful, the file
! is prepared for appending hew data; if it fails, a new file will
! be created.

      create_new_file=ldefhis            !<-- initialize
      if (ncid /= -1)  create_new_file=.false.
#if defined MPI && !defined PARALLEL_FILES
      if (mynode > 0) create_new_file=.false.
#endif


!  Create a new station data file.
!---------------------------------

  10  if (create_new_file) then
        ierr=nf_create(staname(1:lncn), nf_clobber+nf_64bit_offset,
     &                                                         ncid)
        if (ierr /= nf_noerr) then
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_stn :: Cannot ',
     &                'create netCDF file ''', staname(1:lncn), '''.',
     &                                              nf_strerror(ierr)
          goto 99                             !--> ERROR
# ifdef VERBOSE
        else
          write(*,*) 'def_statn ::  created ''', staname(1:lncn), ''''
# endif
        endif

! Put global attributes.
! --- ------ -----------

        call put_global_atts (ncid, ierr)

! Define the dimensions of staggered fields.
!--------------------------------------------

        ierr=my_nf_def_dim (ncid, 'station', nstation, stadim)
# ifdef SOLVE3D
        ierr=my_nf_def_dim (ncid, 's_rho', N,  srdim)
        ierr=my_nf_def_dim (ncid, 's_w',   N,  swdim)
        ierr=my_nf_def_dim (ncid, 'tracer',NT, trcdim)
# endif
        ierr=my_nf_def_dim (ncid, 'time', nf_unlimited, timedim)

!  Define dimension vectors for point variables.

        pgrd(1)=stadim
        pgrd(2)=timedim
# ifdef SOLVE3D

!  Define dimension vectors for cast variables at vertical RHO-points.

        rgrd(1)=srdim
        rgrd(2)=stadim
        rgrd(3)=timedim

!  Define dimension vectors for cast variables at vertical W-points.

        wgrd(1)=swdim
        wgrd(2)=stadim
        wgrd(3)=timedim
# endif


! Define time-invariant environmental variables:
!------- -------------- ------------- ----------

! Grid type switch: Spherical or Cartesian.

        ierr=nf_def_var (ncid, 'spherical', nf_char, 0, 0, varid)
        ierr=nf_put_att_text (ncid,varid, 'long_name', 24,
     &                                     'grid type logical switch')
        ierr=nf_put_att_text (ncid,varid, 'option_T',9,'spherical')
        ierr=nf_put_att_text (ncid,varid, 'option_F',9,'cartesian')

#ifndef SPHERICAL
! Physical dimensions of model domstaain, xl,el (Cartesian grid only).
        ierr=nf_def_var (ncid, 'xl', NF_FTYPE, 0, 0, varid)
        ierr=nf_put_att_text (ncid, varid, 'long_name',   33,
     &                   'domain length in the XI-direction')
        ierr=nf_put_att_text (ncid, varid, 'units', 5, 'meter')

        ierr=nf_def_var (ncid, 'el', NF_FTYPE, 0, 0, varid)
        ierr=nf_put_att_text (ncid, varid, 'long_name',   34,
     &                  'domain length in the ETA-direction')
        ierr=nf_put_att_text (ncid, varid, 'units', 5, 'meter')
#endif

! ipos,jpos indices for station positions

        ierr=nf_def_var (ncid, 'ipos', nf_int, 1, stadim, varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',25,
     &                                  'I-index station positions')
        ierr=nf_put_att_text(ncid,varid,'units',14,'nondimensional')
        ierr=nf_put_att_int(ncid,varid,'valid_min',NF_INT,1,2)
        ierr=nf_put_att_int(ncid,varid,'valid_max',NF_INT,1,LLm)

        ierr=nf_def_var(ncid,'jpos', nf_int, 1,stadim,varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',25,
     &                                  'J-index station positions')
        ierr=nf_put_att_text(ncid,varid,'units',14,'nondimensional')
        ierr=nf_put_att_int(ncid,varid,'valid_min',NF_INT,1,2)
        ierr=nf_put_att_int(ncid,varid,'valid_max',NF_INT,1,MMm)

! Longitude-latitude or cartezian coordinates of RHO-points

# ifdef SPHERICAL
        ierr=nf_def_var(ncid,'lon_rho',NF_FTYPE,1,stadim,varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',23,
     &                                     'longitude of RHO-points')
        ierr=nf_put_att_text(ncid,varid,'units',11, 'degree_east')

        ierr=nf_def_var(ncid,'lat_rho',NF_FTYPE,1,stadim,varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',22,
     &                                      'latitude of RHO-points')
        ierr=nf_put_att_text(ncid,varid,'units',12,'degree_north')

# else
        ierr=nf_def_var(ncid,'x_rho',NF_FTYPE,1,stadim,varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',25,
     &                           'x-locations of RHO-points')
        ierr=nf_put_att_text(ncid,varid,'units',5,'meter')

        ierr=nf_def_var(ncid,'y_rho',NF_FTYPE,1,stadim,varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',25,
     &                           'y-locations of RHO-points')
        ierr=nf_put_att_text(ncid,varid,'units',5,'meter')
# endif
# ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points, if appropriate.

        ierr=nf_def_var(ncid,'angle',NF_FTYPE,1,stadim,varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',30,
     &                           'angle between XI-axis and EAST')
        ierr=nf_put_att_text(ncid,varid,'units',7, 'radians')
# endif
! Bathymetry.

        ierr=nf_def_var(ncid,'h',NF_FTYPE, 1,stadim,varid)
        ierr=nf_put_att_text(ncid,varid,'long_name',24,
     &                                   'bathymetry at RHO-points')
        ierr=nf_put_att_text(ncid,varid,'units',5,'meter')


! Define evolving model variables:
!------- -------- ----- ----------

! Time.

        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                         NF_FTYPE, 1, timedim, stnTime)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, stnTime, 'long_name',
     &                       lvar, vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, stnTime, 'units',
     &                       lvar, vname(3,indxTime)(1:lvar))

! Free-surface and 2D barotropic velocity components

        if (wrtsta(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &                               NF_FOUT, 2, pgrd, stnZ)
          lvar=lenstr(vname(2,indxZ))
          ierr=nf_put_att_text (ncid, stnZ, 'long_name',
     &                         lvar, vname(2,indxZ)(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncid, stnZ, 'units',
     &                         lvar, vname(3,indxZ)(1:lvar))
        endif

        if (wrtsta(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var (ncid, vname(1,indxUb)(1:lvar),
     &                               NF_FOUT, 2, pgrd, stnUb)
          lvar=lenstr(vname(2,indxUb))
          ierr=nf_put_att_text (ncid, stnUb, 'long_name',
     &                         lvar, vname(2,indxUb)(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text (ncid, stnUb, 'units',
     &                         lvar, vname(3,indxUb)(1:lvar))
        endif

        if (wrtsta(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var (ncid, vname(1,indxVb)(1:lvar),
     &                               NF_FOUT, 2, pgrd, stnVb)
          lvar=lenstr(vname(2,indxVb))
          ierr=nf_put_att_text (ncid, stnVb, 'long_name',
     &                         lvar, vname(2,indxVb)(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text (ncid, stnVb, 'units',
     &                         lvar, vname(3,indxVb)(1:lvar))
        endif
# ifdef SOLVE3D

! 3D momenta in XI- and ETA-directions.

        if (wrtsta(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var (ncid, vname(1,indxU)(1:lvar),
     &                               NF_FOUT, 3, rgrd, stnU)
          lvar=lenstr(vname(2,indxU))
          ierr=nf_put_att_text (ncid, stnU, 'long_name',
     &                         lvar, vname(2,indxU)(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncid, stnU, 'units',
     &                         lvar, vname(3,indxU)(1:lvar))
        endif

        if (wrtsta(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var (ncid, vname(1,indxV)(1:lvar),
     &                               NF_FOUT, 3,rgrd, stnV)
          lvar=lenstr(vname(2,indxV))
          ierr=nf_put_att_text (ncid, stnV, 'long_name',
     &                         lvar, vname(2,indxV)(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncid, stnV, 'units',
     &                         lvar, vname(3,indxV)(1:lvar))
        endif

! 3D momentum component in vertical direction.

        if (wrtsta(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncid, vname(1,indxW)(1:lvar),
     &                               NF_FOUT, 3, rgrd, stnW)
          lvar=lenstr(vname(2,indxW))
          ierr=nf_put_att_text (ncid,stnW,'long_name',
     &                         lvar, vname(2,indxW)(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncid, stnW, 'units',
     &                         lvar, vname(3,indxW)(1:lvar))
        endif

! Tracer variables.

        do itrc=1,NT
          if (wrtsta(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_def_var(ncid,vname(1,indxT+itrc-1)(1:lvar),
     &                                 NF_FOUT, 3, rgrd,stnT(itrc))
            lvar=lenstr(vname(2,indxT+itrc-1))
            ierr=nf_put_att_text(ncid,stnT(itrc), 'long_name',
     &                         lvar, vname(2,indxT+itrc-1)(1:lvar))
            lvar=lenstr(vname(3,indxT+itrc-1))
            ierr=nf_put_att_text(ncid,stnT(itrc), 'units', lvar,
     &                               vname(3,indxT+itrc-1)(1:lvar))
          endif
        enddo

! Density anomaly.

        if (wrtsta(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var(ncid, vname(1,indxR)(1:lvar),
     &                                 NF_FOUT, 3, rgrd, stnR)
          lvar=lenstr(vname(2,indxR))
          ierr=nf_put_att_text(ncid, stnR, 'long_name', lvar,
     &                                   vname(2,indxR)(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text(ncid, stnR, 'units', lvar,
     &                                   vname(3,indxR)(1:lvar))
        endif

! Vertical viscosity coefficient.

        if (wrtsta(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncid, vname(1,indxAkv)(1:lvar),
     &                               NF_FOUT, 3, wgrd, stnAkv)
          lvar=lenstr(vname(2,indxAkv))
          ierr=nf_put_att_text (ncid, stnAkv, 'long_name',
     &                         lvar, vname(2,indxAkv)(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncid, stnAkv, 'units',
     &                         lvar, vname(3,indxAkv)(1:lvar))
        endif

! Vertical diffusion coefficient for potential temperature.

        if (wrtsta(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncid, vname(1,indxAkt)(1:lvar),
     &                               NF_FOUT, 3, wgrd, stnAkt)
          lvar=lenstr(vname(2,indxAkt))
          ierr=nf_put_att_text (ncid, stnAkt, 'long_name',
     &                         lvar, vname(2,indxAkt)(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncid, stnAkt, 'units',
     &                         lvar, vname(3,indxAkt)(1:lvar))
        endif
#  ifdef SALINITY
!  Vertical diffusion coefficient for salinity.

        if (wrtsta(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncid, vname(1,indxAks)(1:lvar),
     &                               NF_FOUT, 3, wgrd, stnAks)
          lvar=lenstr(vname(2,indxAks))
          ierr=nf_put_att_text (ncid, stnAks, 'long_name',
     &                         lvar, vname(2,indxAks)(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncid, stnAks, 'units',
     &                         lvar, vname(3,indxAks)(1:lvar))
        endif
#  endif
#  ifdef LMD_KPP
!  Depth of surface planetary boundary layer.

        if (wrtsta(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_def_var (ncid, vname(1,indxHbls)(1:lvar),
     &                               NF_FOUT, 2, pgrd, stnHbls)
          lvar=lenstr(vname(2,indxHbls))
          ierr=nf_put_att_text (ncid, stnHbls, 'long_name',
     &                         lvar, vname(2,indxHbls)(1:lvar))
          lvar=lenstr(vname(3,indxHbls))
          ierr=nf_put_att_text (ncid, stnHbls, 'units',
     &                         lvar, vname(3,indxHbls)(1:lvar))
        endif
#  endif
#  ifdef LMD_BKPP
!  Thickness of bottom boundary layer.

        if (wrtsta(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_def_var (ncid, vname(1,indxHbbl)(1:lvar),
     &                              NF_FOUT, 2, pgrd, stnHbbl)
          lvar=lenstr(vname(2,indxHbls))
          ierr=nf_put_att_text (ncid, stnHbbl, 'long_name',
     &                         lvar, vname(2,indxHbbl)(1:lvar))
          lvar=lenstr(vname(3,indxHbbl))
          ierr=nf_put_att_text (ncid, stnHbbl, 'units',
     &                         lvar, vname(3,indxHbbl)(1:lvar))
        endif
#  endif
# endif /* SOLVE3D */

!  Leave definition mode.
!---------------------------

        ierr=nf_enddef(ncid)
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
         write(*,'(6x,5A,I4)')  'def_statn: Created new netCDF ''',
     &                  'file ''', staname(1:lncn), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif



c--XXXXXXX


!  Write out running parameters.
!-------------------------------
!
!  Time stepping parameters.

!  Write out grid variables.
!---------------------------------------------------
!


! Write time-invariant environmental variables:
!------ -------------- ------------- ----------

!  Grid type switch.

        ierr=nf_inq_varid (ncid, 'spherical', varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_var1_text (ncid, varid, 1,
#ifdef SPHERICAL
     &                                            'T')
#else
     &                                            'F')
#endif
          if (ierr /= nf_noerr) then
            write(*,2) 'spherical', staname(1:lncn)
          endif
        else
          write(*,1) 'spherical', staname(1:lncn)
        endif
        if (ierr /= nf_noerr)  goto 99                  !--> ERROR

#ifndef SPHERICAL
! Physical Dimensions of Model Domain, xl,el (Cartesian grid only).

        ierr=nf_inq_varid (ncid, 'xl', varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_var1_FTYPE (ncid, varid, 1, xl)
          if (ierr /= nf_noerr) then
            write(*,2) 'xl', staname(1:lncn)
          endif
        else
          write(*,1) 'xl', staname(1:lncn)
        endif
        if (ierr /= nf_noerr)  goto 99                  !--> ERROR

        ierr=nf_inq_varid (ncid, 'el', varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_var1_FTYPE (ncid, varid, 1, el)
          if (ierr /= nf_noerr) then
            write(*,2) 'el', staname(1:lncn)
          endif
        else
          write(*,1) 'el', staname(1:lncn)
        endif
        if (ierr /= nf_noerr)  goto 99                  !--> ERROR
#endif

! i,j-indices of station positions

        ierr=nf_inq_varid(ncid, 'ipos', varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_vara_int(ncid,varid,1,nstation,ispos)
          if (ierr /= nf_noerr) then
            write(*,2) 'ipos', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'ipos', staname(1:lncn), nf_strerror(ierr)
        endif

        ierr=nf_inq_varid(ncid,'jpos',varid)
        if (ierr == nf_noerr) then
          ierr=nf_put_vara_int(ncid,varid,1,nstation,jspos)
          if (ierr /= nf_noerr) then
            write(*,2) 'jpos', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'jpos', staname(1:lncn), nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr)  goto 99                  !--> ERROR

# ifdef SPHERICAL
! Longitude-latitude or cartezian coordinates of RHO-points

        ierr=nf_inq_varid (ncid, 'lon_rho', varid)
        if (ierr == nf_noerr) then
          do k=1,nstation
            wrk(k)=lonr(ispos(k),jspos(k))
          enddo
          ierr=nf_put_vara_FTYPE(ncid,varid,1,nstation,wrk)
          if (ierr /= nf_noerr) then
            write(*,2) 'lon_rho', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'lon_rho', staname(1:lncn), nf_strerror(ierr)
        endif

        ierr=nf_inq_varid(ncid, 'lat_rho', varid)
        if (ierr == nf_noerr) then
          do k=1,nstation
            wrk(k)=latr(ispos(k),jspos(k))
          enddo
          ierr=nf_put_vara_FTYPE(ncid,varid,1,nstation,wrk)
          if (ierr /= nf_noerr) then
            write(*,2) 'lat_rho', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'lat_rho', staname(1:lncn), nf_strerror(ierr)
        endif
# else
        ierr=nf_inq_varid (ncid, 'x_rho', varid)
        if (ierr == nf_noerr) then
          do k=1,nstation
            wrk(k)=xr(ispos(k),jspos(k))
          enddo
          ierr=nf_put_vara_FTYPE(ncid,varid,1,nstation,wrk)
          if (ierr /= nf_noerr) then
            write(*,2) 'x_rho', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'x_rho', staname(1:lncn), nf_strerror(ierr)
        endif

        ierr=nf_inq_varid (ncid, 'y_rho', varid)
        if (ierr == nf_noerr) then
          do k=1,nstation
            wrk(k)=yr(ispos(k),jspos(k))
          enddo
          ierr=nf_put_vara_FTYPE(ncid,varid,1,nstation,wrk)
          if (ierr /= nf_noerr) then
            write(*,2) 'y_rho', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'y_rho', staname(1:lncn), nf_strerror(ierr)
        endif
# endif
        if (ierr /= nf_noerr)  goto 99                  !--> ERROR
# ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points, if appropriate.

        ierr=nf_inq_varid (ncid, 'angle', varid)
        if (ierr == nf_noerr) then
          do k=1,nstation
            wrk(k)=angler(ispos(k),jspos(k))
          enddo
          ierr=nf_put_vara_FTYPE (ncid, varid, 1, nstation, wrk)
          if (ierr /= nf_noerr) then
            write(*,2) 'angle', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'angle', staname(1:lncn), nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr)  goto 99                  !--> ERROR
# endif
! Bathymetry.

        ierr=nf_inq_varid (ncid, 'h', varid)
        if (ierr == nf_noerr) then
          do k=1,nstation
            wrk(k)=h(ispos(k),jspos(k))
          enddo
          ierr=nf_put_vara_FTYPE (ncid, varid, 1, nstation, wrk)
          if (ierr /= nf_noerr) then
            write(*,2) 'h', staname(1:lncn), nf_strerror(ierr)
          endif
        else
          write(*,1) 'h', staname(1:lncn), nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr)  goto 99                  !--> ERROR

! Open existing file and prepare for appending data:   Inquire about
! ==== ======== ==== === ======= === ========= ====    dimensions and
! variables. Check for consistency with model dimensions. In the case
! when file is rejected (whether it cannot be opened, or something is
! wrong with its dimensions) create a new file.  After that verify
! that all necessary variables are already defined, find and save
! their netCDF IDs.

      elseif (ncid==-1) then  !--> alt. to "create_new_file"
        ierr=nf_open (staname(1:lncn), nf_write, ncid)
        if (ierr == nf_noerr) then
#ifdef MPI_SILENT_MODE
          if (mynode == 0) then
#endif
            write(*,'(1x,5A,I4)')  'def_statn: Opened existing ',
     &           'netCDF file ''', staname(1:lncn),  '''.' MYID
#ifdef MPI_SILENT_MODE
          endif
#endif
        else
#if defined MPI && !defined PARALLEL_FILES
          if (mynode == 0) then
            create_new_file=.true.
            goto 10
          else
            write(*,'(/1x,4A/12x,A/)')     '### ERROR: def_statn :: ',
     &                  'Cannot open netCDF file ''', staname(1:lncn),
     &                         ''' for writing.',   nf_strerror(ierr)
            goto 99                                         !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif


!  Scan variable list from input NetCDF and activate switches for
!  stations variables. Get variable IDs.

        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),stnTime)
        if (ierr /= nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
          goto 99                                         !--> ERROR
        endif

! Free-surface and 2D barotropic velocity components

        if (wrtsta(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), stnZ)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxZ)(1:lvar),  staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrtsta(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), stnUb)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxUb)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrtsta(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), stnVb)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxVb)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

# ifdef SOLVE3D

! 3D momenta in XI- and ETA-directions.

        if (wrtsta(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), stnU)
          if (ierr /= nf_noerr) then
            write(*,1)  vname(1,indxU)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrtsta(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), stnV)
          if (ierr /= nf_noerr) then
            write(*,1)  vname(1,indxV)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

! Tracer variables.

        do itrc=1,NT
          if (wrtsta(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_inq_varid (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  stnT(itrc))
            if (ierr /= nf_noerr) then
              write(*,1)  vname(1,indxT+itrc-1)(1:lvar),
     &               staname(1:lncn), nf_strerror(ierr)
              goto 99                                     !--> ERROR
            endif
          endif
        enddo

! Density anomaly.

        if (wrtsta(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid (ncid, vname(1,indxR)(1:lvar), stnR)
          if (ierr /= nf_noerr) then
            write(*,1)  vname(1,indxR)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

! S-coordinate "omega" and "true" W- vertical velocities.

        if (wrtsta(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid (ncid, vname(1,indxO)(1:lvar), stnO)
          if (ierr /= nf_noerr) then
            write(*,1)  vname(1,indxO)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrtsta(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid (ncid, vname(1,indxW)(1:lvar), stnW)
          if (ierr /= nf_noerr) then
            write(*,1)  vname(1,indxW)(1:lvar), staname(1:lncn),
     &                                        nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

! Vertical viscosity, thermal conductivity, and salinity diffusion.

        if (wrtsta(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncid, vname(1,indxAkv)(1:lvar), stnAkv)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxAkv)(1:lvar), staname(1:lncn),
     &                                         nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrtsta(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncid,vname(1,indxAkt)(1:lvar), stnAkt)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxAkt)(1:lvar), staname(1:lncn),
     &                                         nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif

# ifdef SALINITY
        if (wrtsta(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncid,vname(1,indxAks)(1:lvar), stnAks)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxAks)(1:lvar), staname(1:lncn),
     &                                         nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif
# endif
# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        if (wrtsta(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_inq_varid (ncid,vname(1,indxHbls)(1:lvar),stnHbls)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxHbls)(1:lvar), staname(1:lncn),
     &                                          nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        if (wrtsta(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_inq_varid (ncid,vname(1,indxHbbl)(1:lvar),stnHbbl)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxHbbl)(1:lvar), staname(1:lncn),
     &                                          nf_strerror(ierr)
            goto 99                                       !--> ERROR
          endif
        endif
# endif
#endif /* SOLVE3D */
# ifdef VERBOSE
      else
        write(*,*) 'def_statn :: no action is needs to be taken'
# endif
      endif

  1   format(/1x,'### ERROR: def_statn :: Cannot find variable ''',
     &      A, ''' within netCDF file' /12x, '''', A, ''':', 1x, A/)
  2   format(/1x,'### ERROR: def_statn :: Cannot write variable ''',
     &        A, ''' into netCDF file' /12x, '''', A, ''':', 1x, A/)
 99   return
      end
#else
      subroutine def_stn_empty
      end
#endif /* STATIONS */
