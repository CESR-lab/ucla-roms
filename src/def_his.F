#ifndef AVRH
# include "cppdefs.h"

      subroutine def_his(total_rec, ierr)
#else
      subroutine def_avg(total_rec, ierr)
#endif

! Create/open history/averages netCDF file so the subsequent calls
! to wrt_his/avg can write model variables accessing them by netCDF
! IDs and proper record number.  This routine can work in two modes:
! either create a new file and save all netCDF variable IDs, and set
! record index to 1 (starting a new file), or open an existing file,
! verify that all dimensions are consistent with the running model
! and all required variables are present, find and save their variable
! IDs and set time record to start appending data.

! The difference between def_his and def_avg is as follows: they have
! different netCDF file name (hisname/avgname); netCDF file ID (passed
! as argument); time record index (hisindx/avgindx); array of switches
! to select which variables to write (wrthis/wrtavg); different sets
! of netCDF variable IDs (hisTime...hisHbl/avgTime...avgHbl); and the
! first attribute of each variable, long_name, has prefix 'averaged'.
! Because most of the code is identical for both routines, the second
! one is generated from the first entirely by CPP.

c--#define nf_def_var my_nf_def_var

#ifdef WEC
      use wec_frc, only: def_his_def_var_wec, def_avg_def_var_wec
#endif
      use tracers ! def_his_trc, def_avg_trc
      use read_write, only: output_root_name
      use bgc ! def_his_bgc_diag, def_avg_bgc_diag
      use surf_flux ! set_avg_surf_flux

      implicit none
      integer total_rec, ierr
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"
      logical create_new_file
      integer lenstr, checkdims, my_nf_def_dim
      integer rec, lfnm, lvar, timedim, prev_fill_mode,
     &        r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)
#ifdef SOLVE3D
     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc
#endif
      character(len=64) text
#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real*8, parameter :: spv_set=1.D+33
# else
      real*4, parameter :: spv_set=1.E+33
# endif
#endif

#ifndef AVRH
# define ncid nchis
# define fname hisname
# define rec_per_file nrpfhis
# define wrt wrthis
# define vidTime hisTime
# define vidZ hisZ
# define vidUb hisUb
# define vidVb hisVb
# define vidU hisU
# define vidV hisV
# define vidT hisT
# define vidR hisR
# define vidO hisO
# define vidW hisW
# define vidAkv hisAkv
# define vidAkt hisAkt
# define vidAks hisAks
# define vidHbls hisHbls
# define vidHbbl hisHbbl
# ifdef WKB_WWAVE
#  define vidWKB hisWKB
# endif
#else
# define ncid ncavg
# define fname avgname
# define rec_per_file nrpfavg
# define wrt wrtavg
# define vidTime avgTime
# define vidZ avgZ
# define vidUb avgUb
# define vidVb avgVb
# define vidU avgU
# define vidV avgV
# define vidT avgT
# define vidR avgR
# define vidO avgO
# define vidW avgW
# define vidAkv avgAkv
# define vidAkt avgAkt
# define vidAks avgAks
# define vidHbls avgHbls
# define vidHbbl avgHbbl
# ifdef WKB_WWAVE
#  define vidWKB avgWKB
# endif
#endif

! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'his_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.

      ierr=0
#ifndef AVRH
      fname=trim(output_root_name) / / '_his.nc'
#else
      fname=trim(output_root_name) / / '_avg.nc'
#endif
      lfnm=lenstr(fname)
#if defined MPI && defined PARALLEL_FILES
        ! Insert MPI node numbers to file name
        call insert_node(fname, lfnm, mynode, NNODES, ierr)
#endif
!#endif
      lfnm=lenstr(fname)
      if (rec_per_file > 0) then
        lvar=total_rec-(1+mod(total_rec-1, rec_per_file))
        call insert_time_index(fname, lfnm, lvar,  ierr)
        if (ierr /= 0) goto 99
      endif

! Decide whether to create a new file, or open an existing one.
! Overall the code below is globally organized as a 3-way switch,
!
! 10  if (create_new_file) then
!        ... create new file, save netCDF ids for all variables;
!     elseif (ncid == -1) then
!        ... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be already
!        .... known (applies only for single file output from MPI code)
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

! Create new history/averages file:    Put global attributes and
! ====== === ======= ======== =====    define all the variables.

  10  if (create_new_file) then

        ierr=nf_create(fname(1:lfnm), nf_clobber+nf_64bit_offset,
     &                                                     ncid)
        if (ierr == nf_noerr) then

          ! Set fill value handling
          ierr=nf_set_fill(ncid, nf_nofill, prev_fill_mode)

          ! Error handling
          if (ierr /= nf_noerr) then
# ifdef AVRH
            write(*,'(1x,4A,I4)') '### ERROR: def_avg :: Cannot ',
# else
            write(*,'(1x,4A,I4)') '### ERROR: def_his :: Cannot ',
# endif
     &        'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &         MYID
          endif

        ! Error: cannot create file
        else
#ifdef AVRH
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_avg :: ',
#else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_his :: ',
#endif
     &                'Cannot create ''', fname(1:lfnm), '''.',
     &                                       nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr) goto 99
        if (rec_per_file == 0) total_rec=0

! Put global attributes:
!---- ------ -----------

        call put_global_atts(ncid, ierr)

! Define dimensions for field placement on staggered grids:
! ------ ---------- --- ----- --------- -- --------- -------

        ierr=my_nf_def_dim(ncid, 'xi_rho',  xi_rho,  r2dgrd(1))
        ierr=my_nf_def_dim(ncid, 'xi_u',    xi_u,    u2dgrd(1))
        ierr=my_nf_def_dim(ncid, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=my_nf_def_dim(ncid, 'eta_v',   eta_v,   v2dgrd(2))
#ifdef SOLVE3D
        ierr=my_nf_def_dim(ncid, 's_rho',   N,       r3dgrd(3))
        ierr=my_nf_def_dim(ncid, 's_w',     N+1,     w3dgrd(3))
#endif
        ierr=my_nf_def_dim(ncid, 'time', nf_unlimited, timedim)
        ierr=my_nf_def_dim(ncid, 'auxil',   iaux,     auxil(1))

      ! Transfer all newly defined dimensions to grid types that use
      ! that dimension - e.g. u-point = ( xi_u, eta_rho, timedim )
      ! for eta_rho it can use u2dgrd(2)=r2dgrd(2) rather than nf call
      ! --------------------------------------------------------------

        auxil(2)=timedim

        r2dgrd(3)=timedim                          ! Free surface

        u2dgrd(2)=r2dgrd(2) ; u2dgrd(3)=timedim    ! 2D UBAR-type

        v2dgrd(1)=r2dgrd(1) ; v2dgrd(3)=timedim    ! 2D VBAR-type

#ifdef SOLVE3D
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2) ; r3dgrd(4)=timedim    ! 3D RHO-type

        u3dgrd(1)=u2dgrd(1) ; u3dgrd(3)=r3dgrd(3)  ! 3D U-type
        u3dgrd(2)=r2dgrd(2) ; u3dgrd(4)=timedim    ! 

        v3dgrd(1)=r2dgrd(1) ; v3dgrd(3)=r3dgrd(3)  ! 3D V-type
        v3dgrd(2)=v2dgrd(2) ; v3dgrd(4)=timedim    !

        w3dgrd(1)=r2dgrd(1) ; w3dgrd(4)=timedim    ! 3D W-type
        w3dgrd(2)=r2dgrd(2)
#endif
#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)

! Define grid variables.
! ------ ---- ----------

        if (total_rec <= 1) call def_grid(ncid, r2dgrd)
#endif

! Define evolving model variables.
! ------ -------- ----- ----------

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
#ifndef AVRH
        ierr=nf_def_var(ncid, 'time_step', nf_int, 2, auxil, hisTstep)
        ierr=nf_put_att_text (ncid, hisTstep, 'long_name', 48,
     &             'time step and record numbers from initialization')
#endif

! Time.
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var(ncid, vname(1,indxTime)(1:lvar), NF_FTYPE,
     &                                        1, timedim, vidTime)
#ifdef AVRH
        text='averaged '/ /vname(2,indxTime)
#else
        text=vname(2,indxTime)
#endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                            text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'units',  lvar,
     &                               vname(3,indxTime)(1:lvar))

! Free-surface and 2D barotropic velocity components

        if (wrt(indxZ) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var(ncid, vname(1,indxZ)(1:lvar), NF_FOUT,
     &                                         3, r2dgrd, vidZ)
#ifdef AVRH
          text='averaged '/ /vname(2,indxZ)
#else
          text=vname(2,indxZ)
#endif
          lvar=lenstr(text)
          ierr=nf_put_att_text(ncid, vidZ, 'long_name', lvar,
     &                                          text(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text(ncid, vidZ, 'units',     lvar,
     &                                 vname(3,indxZ)(1:lvar))
#ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT(ncid, vidZ, '_FillValue', NF_FOUT,
     &                                               1, spv_set)
          if (ierr /= nf_noerr) then
           write(*,*)  'nf_put_att_XXX:', nf_strerror(ierr)
          endif
#endif
        endif

        if (wrt(indxUb) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var(ncid, vname(1,indxUb)(1:lvar), NF_FOUT,
     &                                         3, u2dgrd, vidUb)
#ifdef AVRH
          text='averaged '/ /vname(2,indxUb)
#else
          text=vname(2,indxUb)
#endif
          lvar=lenstr(text)
          ierr=nf_put_att_text(ncid, vidUb, 'long_name', lvar,
     &                                          text(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text(ncid, vidUb, 'units',     lvar,
     &                               vname(3,indxUb)(1:lvar))
#ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT(ncid, vidUb, '_FillValue', NF_FOUT,
     &                                                1, spv_set)
#endif
        endif

        if (wrt(indxVb) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var(ncid, vname(1,indxVb)(1:lvar), NF_FOUT,
     &                                         3, v2dgrd, vidVb)
#ifdef AVRH
          text='averaged '/ /vname(2,indxVb)
#else
          text=vname(2,indxVb)
#endif
          lvar=lenstr(text)
          ierr=nf_put_att_text(ncid, vidVb, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text(ncid, vidVb, 'units',     lvar,
     &                                 vname(3,indxVb)(1:lvar))
#ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT(ncid, vidVb, '_FillValue', NF_FOUT,
     &                                                1, spv_set)
#endif
        endif
#ifdef SOLVE3D

! 3D momenta in XI- and ETA-directions.

        if (wrt(indxU) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var(ncid, vname(1,indxU)(1:lvar), NF_FOUT,
     &                                         4, u3dgrd, vidU)
# ifdef AVRH
          text='averaged '/ /vname(2,indxU)
# else
          text=vname(2,indxU)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidU, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncid, vidU, 'units',     lvar,
     &                                  vname(3,indxU)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidU, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

        if (wrt(indxV) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var(ncid, vname(1,indxV)(1:lvar), NF_FOUT,
     &                                           4, v3dgrd, vidV)
# ifdef AVRH
          text='averaged '/ /vname(2,indxV)
# else
          text=vname(2,indxV)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidV, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncid, vidV, 'units',     lvar,
     &                                  vname(3,indxV)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidV, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

      ! Tracer variables
      ! ----------------
# ifndef AVRH /* NOT AVRH */
        ! Define tracer his variables to output
        call def_his_trc(r3dgrd)
# else /* AVRH */
        ! Define tracers AVERAGED variables to output
        call def_avg_trc(r3dgrd)
# endif /* END AVRH */

! Density anomaly.

        if (wrt(indxR) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var (ncid, vname(1,indxR)(1:lvar), NF_FOUT,
     &                                           4, r3dgrd, vidR)
# ifdef AVRH
          text='averaged '/ /vname(2,indxR)
# else
          text=vname(2,indxR)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidR, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text (ncid, vidR, 'units',     lvar,
     &                                  vname(3,indxR)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidR, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

! S-coordinate "omega" vertical velocity.

        if (wrt(indxO) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_def_var (ncid, vname(1,indxO)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidO)
# ifdef AVRH
          text='averaged '/ /vname(2,indxO)
# else
          text=vname(2,indxO)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidO, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxO))
          ierr=nf_put_att_text (ncid, vidO, 'units',     lvar,
     &                                  vname(3,indxO)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidO, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

! True W-vertical velocity.

        if (wrt(indxW) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncid, vname(1,indxW)(1:lvar), NF_FOUT,
     &                                           4, r3dgrd, vidW)
# ifdef AVRH
          text='averaged '/ /vname(2,indxW)
# else
          text=vname(2,indxW)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidW, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncid, vidW, 'units',     lvar,
     &                                  vname(3,indxW)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidW, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

! Vertical viscosity coefficient.

        if (wrt(indxAkv) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncid, vname(1,indxAkv)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAkv)
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkv)
# else
          text=vname(2,indxAkv)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAkv, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncid, vidAkv, 'units',     lvar,
     &                                  vname(3,indxAkv)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidAkv, '_FillValue', NF_FOUT,
     &                                                   1, spv_set)
# endif
        endif

! Vertical thermal conductivity coefficient.

        if (wrt(indxAkt) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncid, vname(1,indxAkt)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAkt)
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkt)
# else
          text=vname(2,indxAkt)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAkt, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncid, vidAkt, 'units',     lvar,
     &                                  vname(3,indxAkt)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidAkt, '_FillValue', NF_FOUT,
     &                                                   1, spv_set)
# endif
        endif
# ifdef SALINITY
! Vertical diffusion coefficient for salinity.

        if (wrt(indxAks) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncid, vname(1,indxAks)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAks)
#  ifdef AVRH
          text='averaged '/ /vname(2,indxAks)
#  else
          text=vname(2,indxAks)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAks, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncid, vidAks, 'units',     lvar,
     &                                  vname(3,indxAks)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidAks, '_FillValue', NF_FOUT,
     &                                                   1, spv_set)
#  endif
        endif
# endif /* SALINITY */
# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        if (wrt(indxHbls) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_def_var (ncid, vname(1,indxHbls)(1:lvar),
     &                           NF_FOUT, 3, r2dgrd, vidHbls)
#  ifdef AVRH
          text='averaged '/ /vname(2,indxHbls)
#  else
          text=vname(2,indxHbls)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidHbls, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbls))
          ierr=nf_put_att_text (ncid, vidHbls, 'units',     lvar,
     &                                 vname(3,indxHbls)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidHbls, '_FillValue',
     &                                       NF_FOUT, 1, spv_set)
#  endif
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        if (wrt(indxHbbl) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_def_var (ncid, vname(1,indxHbbl)(1:lvar),
     &                            NF_FOUT, 3, r2dgrd, vidHbbl)
#  ifdef AVRH
          text='averaged '/ /vname(2,indxHbbl)
#  else
          text=vname(2,indxHbbl)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidHbbl, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbbl))
          ierr=nf_put_att_text (ncid, vidHbbl, 'units',     lvar,
     &                                 vname(3,indxHbbl)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidHbbl, '_FillValue',
     &                                       NF_FOUT, 1, spv_set)
#  endif
        endif
# endif
#endif /* SOLVE3D */


      ! Define WEC related variable attributes
      ! --------------------------------------
#ifdef WEC
# ifndef AVRH /* NOT AVRH */
        ! Define WEC variables to be output
        call def_his_def_var_wec(ierr, r2dgrd, u2dgrd, v2dgrd
#  ifdef SOLVE3D
     &                              , r3dgrd, u3dgrd, v3dgrd, w3dgrd
#  endif
     &                                                             )
# else /* AVRH */
        ! Define WEC AVERAGED variables to be output
        call def_avg_def_var_wec(ierr, r2dgrd, u2dgrd, v2dgrd
#  ifdef SOLVE3D
     &                              , r3dgrd, u3dgrd, v3dgrd, w3dgrd
#  endif
     &                                                             )
# endif /* END AVRH */
#endif   /* WEC */

      ! Define BGC diag variables
      ! -------------------------
#if defined BIOLOGY_BEC2 && defined BEC2_DIAG
# ifndef AVRH /* NOT AVRH */
        ! Define bgc BEC2_DIAG his variables to output
        call def_his_bgc_diag(r2dgrd
#  ifdef SOLVE3D
     &                      , r3dgrd
#  endif
     &                               )
# else /* AVRH */
        ! Define bgc BEC2_DIAG AVERAGED variables to output
        call def_avg_bgc_diag(r2dgrd
#  ifdef SOLVE3D
     &                      , r3dgrd
#  endif
     &                               )
# endif /* END AVRH */
#endif   /* BEC2_DIAG */

      ! Define surface flux variables (sustr,svstr,etc) averages only
#ifdef AVRH
      call def_avg_surf_flux(r2dgrd, u2dgrd, v2dgrd)
#endif


! Leave definition mode.
! ----- ---------- -----

        ierr=nf_enddef(ncid)
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
#ifdef AVRH
         write(*,'(7x,4A,I4)') 'def_avg :: created new netCDF file ''',
#else
         write(*,'(7x,4A,I4)') 'def_his :: created new netCDF file ''',
#endif
     &                          fname(1:lfnm), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

! Open existing file and prepare for appending data:   Inquire about
! ==== ======== ==== === ======= === ========= ====    its dimensions
! and variables; check for consistency with model dimensions -- if the
! file is rejected (it cannot be opened, or its dimensions do not match
! the ones in the mode) create a new file.  After that verify
! that all necessary variables are already defined, find and save
! their netCDF IDs.

      elseif (ncid == -1) then  !--> alt. to "create_new_file"
        ierr=nf_open(fname(1:lfnm), nf_write, ncid)
        if (ierr == nf_noerr) then
          ierr=nf_set_fill(ncid, nf_nofill, prev_fill_mode)
          if (ierr == nf_noerr) then
#ifdef MPI_SILENT_MODE
            if (mynode == 0) then
#endif
#ifdef AVRH
              write(*,'(7x,4A,I4)') 'def_avg :: open existing file ''',
#else
              write(*,'(7x,4A,I4)') 'def_his :: open existing file ''',
#endif
     &                        fname(1:lfnm),  '''.' MYID
#ifdef MPI_SILENT_MODE
            endif
#endif
          else
# ifdef AVRH
            write(*,'(1x,4A,I4)') '### ERROR: def_avg :: Cannot ',
# else
            write(*,'(1x,4A,I4)') '### ERROR: def_his :: Cannot ',
# endif
     &         'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &          MYID
            goto 99
          endif

          ierr=checkdims(ncid, fname, rec)
          if (ierr == nf_noerr) then
            if (rec_per_file == 0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, rec_per_file))
            endif
            if (ierr > 0) then
              mpi_master_only write(*,
     &                '(/1x,A,I5,1x,3A/21x,2(A,I5),1x,A/21x,A/)')
#ifdef AVRH
     &            'WARNING: def_avg :: The actual number of records',
#else
     &            'WARNING: def_his :: The actual number of records',
#endif
     &             rec, 'present in file ''',  fname(1:lfnm),  '''',
     &            'exceeds record', rec+1-ierr, '/', total_rec,
     &            'specified by restart initial conditions.',
     &        'All records beyond this number will be overwritten.'
              rec=rec-ierr
            elseif (rec_per_file == 0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI && !defined PARALLEL_FILES
              if (mynode > 0) total_rec=total_rec-1
#endif
            endif
            ierr=nf_noerr
          endif
        endif

        if (ierr /= nf_noerr) then
#if defined MPI && !defined PARALLEL_FILES
          if (mynode == 0) then
            create_new_file=.true. ; goto 10 !--> reject and create new
          else
# ifdef AVRH
            write(*,'(/1x,5A,I4/)')  '### ERROR: def_avg :: Cannot ',
# else
            write(*,'(/1x,5A,I4/)')  '### ERROR: def_his :: Cannot ',
# endif
     &        'open ''', fname(1:lfnm), '''.' MYID
            goto 99                                       !--> ERROR
          endif
#else
          create_new_file=.true. ; goto 10   !--> reject and create new
#endif
        endif

! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------

! Time-step indices, ocean time:

#ifndef AVRH
        ierr=nf_inq_varid(ncid, 'time_step', hisTstep)
        if (ierr /= nf_noerr) write(*,1) 'time_step', fname(1:lfnm)
#endif
        if (ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxTime))
          ierr=nf_inq_varid(ncid,vname(1,indxTime)(1:lvar), vidTime)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxTime)(1:lvar),
     &                                                 fname(1:lfnm)
        endif

! Free-surface and 2D barotropic velocity components

        if (wrt(indxZ) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid(ncid, vname(1,indxZ)(1:lvar), vidZ)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxZ)(1:lvar),
     &                                              fname(1:lfnm)
        endif
        if (wrt(indxUb) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid(ncid, vname(1,indxUb)(1:lvar), vidUb)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxUb)(1:lvar),
     &                                               fname(1:lfnm)
        endif
        if (wrt(indxVb) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid(ncid, vname(1,indxVb)(1:lvar), vidVb)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxVb)(1:lvar),
     &                                               fname(1:lfnm)
        endif

#ifdef SOLVE3D

! 3D momenta in XI- and ETA-directions.

        if (wrt(indxU) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid(ncid, vname(1,indxU)(1:lvar), vidU)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxU)(1:lvar),
     &                                              fname(1:lfnm)
        endif
        if (wrt(indxV)  .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid(ncid, vname(1,indxV)(1:lvar), vidV)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxV)(1:lvar),
     &                                              fname(1:lfnm)
        endif

! Tracer variables.

! DevinD: no longer need to set var id as handled by read_write.F module
!        do itrc=1,NT
!          if (wrt_t(itrc) .and. ierr == nf_noerr) then
!            lvar=lenstr(t_vname(itrc))
!            ierr=nf_inq_varid(ncid, t_vname(itrc)(1:lvar),
!     &                                                 vidT(itrc))
!            if (ierr /= nf_noerr) write(*,1)
!     &               t_vname(itrc)(1:lvar), fname(1:lfnm)
!          endif
!        enddo

! Density anomaly.

        if (wrt(indxR) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid(ncid, vname(1,indxR)(1:lvar), vidR)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxR)(1:lvar),
     &                                              fname(1:lfnm)
        endif

! S-coordinate "omega" and "true" W- vertical velocities.

        if (wrt(indxO) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid(ncid, vname(1,indxO)(1:lvar), vidO)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxO)(1:lvar),
     &                                              fname(1:lfnm)
        endif

        if (wrt(indxW) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid(ncid, vname(1,indxW)(1:lvar), vidW)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxW)(1:lvar),
     &                                              fname(1:lfnm)
        endif

! Vertical viscosity, thermal conductivity, and salinity diffusion.

        if (wrt(indxAkv) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid(ncid, vname(1,indxAkv)(1:lvar), vidAkv)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxAkv)(1:lvar),
     &                                                fname(1:lfnm)
        endif

        if (wrt(indxAkt) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid(ncid,vname(1,indxAkt)(1:lvar), vidAkt)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxAkt)(1:lvar),
     &                                                fname(1:lfnm)
        endif

# ifdef SALINITY
        if (wrt(indxAks) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid(ncid,vname(1,indxAks)(1:lvar), vidAks)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxAks)(1:lvar),
     &                                                fname(1:lfnm)
        endif
# endif
# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        if (wrt(indxHbls) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_inq_varid(ncid,vname(1,indxHbls)(1:lvar), vidHbls)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxHbls)(1:lvar),
     &                                                 fname(1:lfnm)
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        if (wrt(indxHbbl) .and. ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_inq_varid(ncid,vname(1,indxHbbl)(1:lvar), vidHbbl)
          if (ierr /= nf_noerr) write(*,1) vname(1,indxHbbl)(1:lvar),
     &                                                 fname(1:lfnm)
        endif
# endif
#endif /* SOLVE3D */



! ---------------------
! Leave definition mode
! ---------------------

        if (ierr == nf_noerr) then
#ifdef MPI_SILENT_MODE
          if (mynode == 0) then
#endif
#ifdef AVRH
            write(*,'(7x,4A,I4,2x,A,I4)')   'def_avg :: open ',
#else
            write(*,'(7x,4A,I4,2x,A,I4)')   'def_his :: open ',
#endif
     &     'existing file ''',fname(1:lfnm),''' from rec =',rec
     &      MYID
#ifdef MPI_SILENT_MODE
          endif
#endif
        endif

#if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf_open(fname(1:lfnm), nf_write, ncid)
        if (ierr == nf_noerr) then
          ierr=nf_set_fill(ncid, nf_nofill, prev_fill_mode)
          if (ierr == nf_noerr) then
# ifndef MPI_SILENT_MODE
#  ifdef AVRH
            write(*,'(7x,5A,I4)')    'def_avg :: reopen existing ',
#  else
            write(*,'(7x,5A,I4)')    'def_his :: reopen existing ',
#  endif
     &              'file ''', fname(1:lfnm), '''.' MYID
# endif
          else
# ifdef AVRH
            write(*,'(1x,4A,I4)')  '### ERROR: def_avg :: Cannot ',
# else
            write(*,'(1x,4A,I4)')  '### ERROR: def_his :: Cannot ',
# endif
     &         'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &          MYID
          endif
        else
# ifdef AVRH
          write(*,'(/1x,6A,I4/)')  '### ERROR: def_avg :: Cannot ',
# else
          write(*,'(/1x,6A,I4/)')  '### ERROR: def_his :: Cannot ',
# endif
     &                      'reopen file ''', fname(1:lfnm), '''.',
     &                       nf_strerror(ierr) MYID
        endif
#endif
      endif             !<--  create_new_file
      if (ierr /= nf_noerr) goto 99

#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)

! Write grid variables.

        if (total_rec <= 1) call wrt_grid (ncid, fname, lfnm)
#endif

  99  return

#ifdef AVRH
   1  format(/1x,'### ERROR: def_avg :: Cannot find variable ''',
#else
   1  format(/1x,'### ERROR: def_his :: Cannot find variable ''',
#endif
     &                            A, ''' in file ''', A, '''.'/)
      end

#undef ncid
#undef fname
#undef rec_per_file
#undef wrt
#undef vidTime
#undef vidZ
#undef vidUb
#undef vidVb
#undef vidU
#undef vidV
#undef vidT
#undef vidR
#undef vidO
#undef vidW
#undef vidAkv
#undef vidAkt
#undef vidAks
#undef vidHbls
#undef vidHbbl
#ifdef WKB_WWAVE
# undef vidWKB
#endif

#ifndef AVRH
# define AVRH
# ifdef AVERAGES
#  include "def_his.F"
# endif
#endif
