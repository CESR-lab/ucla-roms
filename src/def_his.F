#ifndef AVRH
# include "cppdefs.opt"

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

      use param
      use wec_frc
      use tracers                                          ! def_his_trc, def_avg_trc
      use roms_read_write
      use nc_read_write
      use bgc                                              ! def_his_bgc_diag, def_avg_bgc_diag
      use netcdf
      use ncvars
      use scalars
      use netcdf
      use grid                                             ! def_grid

      implicit none

      ! input
      integer :: total_rec, ierr
      ! local
      logical :: create_new_file
      integer :: lenstr, checkdims, my_nf_def_dim, timedimid(1), varid, itrc
      integer :: rec, lfnm, lvar, timedim, prev_fill_mode
      character(len=64) :: text
#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real(kind=8), parameter :: spv_set=1.D+33
# else
      real(kind=4), parameter :: spv_set=1.E+33
# endif
#endif

! DEVIND GET RID OF THIS
#ifndef AVRH
# define ncid nchis
# define fname hisname
# define rec_per_file nrpfhis
# define wrt wrthis
#else
# define ncid ncavg
# define fname avgname
# define rec_per_file nrpfavg
# define wrt wrtavg
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
! DevinD this is done in ncdf_create_file:
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

        ierr=nf90_create(fname(1:lfnm), NF90_NETCDF4, ncid)
        if (ierr == nf90_noerr) then

          ! Set fill value handling
          ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

          ! Error handling
          if (ierr /= nf90_noerr) then
# ifdef AVRH
            write(*,'(1x,4A,I4)') '### ERROR: def_avg :: Cannot ',
# else
            write(*,'(1x,4A,I4)') '### ERROR: def_his :: Cannot ',
# endif
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
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
     &                                       nf90_strerror(ierr)
        endif
        if (ierr /= nf90_noerr) goto 99
        if (rec_per_file == 0) total_rec=0

        call put_global_atts(ncid, ierr)

#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)
        if (total_rec <= 1) call def_grid(ncid)
#endif

! Define evolving model variables.
! ------ -------- ----- ----------

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
#ifndef AVRH
        varid = nccreate(ncid,'time_step',(/'auxil','time'/),(/iaux,0/),nf90_int)
        ierr=nf90_put_att (ncid, varid, 'long_name',
     &             'time step and record numbers from initialization')
#endif

! Time.
#ifdef AVRH
        text='averaged '/ /vname(2,indxTime)
#else
        text=vname(2,indxTime)
#endif
        ierr=nf90_inq_dimid(ncid, 'time', timedimid(1))  ! DevinD remove later?
        varid = nccreate(ncid,vname(1,indxTime),(/'time'/),(/0/),nf90_double)  ! store ocean_time as double
        ierr = nf90_put_att(ncid,varid,'long_name',text)                       ! as need extra significant figs
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxTime))

! Free-surface and 2D barotropic velocity components

        if (wrt(indxZ) .and. ierr == nf90_noerr) then
#ifdef AVRH
          text='averaged '/ /vname(2,indxZ)
#else
          text=vname(2,indxZ)
#endif
          varid = nccreate(ncid,vname(1,indxZ),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxZ))
        endif

        if (wrt(indxUb) .and. ierr == nf90_noerr) then
#ifdef AVRH
          text='averaged '/ /vname(2,indxUb)
#else
          text=vname(2,indxUb)
#endif
          varid = nccreate(ncid,vname(1,indxUb),(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxUb))
        endif

        if (wrt(indxVb) .and. ierr == nf90_noerr) then
#ifdef AVRH
          text='averaged '/ /vname(2,indxVb)
#else
          text=vname(2,indxVb)
#endif
          varid = nccreate(ncid,vname(1,indxVb),(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxVb))
        endif
#ifdef SOLVE3D

! 3D momenta in XI- and ETA-directions.

        if (wrt(indxU) .and. ierr == nf90_noerr) then
# ifdef AVRH
          text='averaged '/ /vname(2,indxU)
# else
          text=vname(2,indxU)
# endif
          varid = nccreate(ncid,vname(1,indxU),(/'xi_u','eta_rho','s_rho','time'/),(/xi_u,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxU))
        endif

        if (wrt(indxV) .and. ierr == nf90_noerr) then
# ifdef AVRH
          text='averaged '/ /vname(2,indxV)
# else
          text=vname(2,indxV)
# endif
          varid = nccreate(ncid,vname(1,indxV),(/'xi_rho','eta_v','s_rho','time'/),(/xi_rho,eta_v,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxV))
        endif

      ! Tracer variables
      ! ----------------
# ifndef AVRH /* NOT AVRH */
        ! Define tracer his variables to output
        call def_his_trc
# else /* AVRH */
        ! Define tracers AVERAGED variables to output
        call def_avg_trc
# endif /* END AVRH */

! Density anomaly.

        if (wrt(indxR) .and. ierr == nf90_noerr) then
# ifdef AVRH
          text='averaged '/ /vname(2,indxR)
# else
          text=vname(2,indxR)
# endif
          varid = nccreate(ncid,vname(1,indxR),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxR))
        endif

! S-coordinate "omega" vertical velocity.

        if (wrt(indxO) .and. ierr == nf90_noerr) then
# ifdef AVRH
          text='averaged '/ /vname(2,indxO)
# else
          text=vname(2,indxO)
# endif
          varid = nccreate(ncid,vname(1,indxO),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxO))
        endif

! True W-vertical velocity.

        if (wrt(indxW) .and. ierr == nf90_noerr) then
# ifdef AVRH
          text='averaged '/ /vname(2,indxW)
# else
          text=vname(2,indxW)
# endif
          varid = nccreate(ncid,vname(1,indxW),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxW))
        endif

! Vertical viscosity coefficient.

        if (wrt(indxAkv) .and. ierr == nf90_noerr) then
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkv)
# else
          text=vname(2,indxAkv)
# endif
          varid = nccreate(ncid,vname(1,indxAkv),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkv))
        endif

! Vertical thermal conductivity coefficient.

        if (wrt(indxAkt) .and. ierr == nf90_noerr) then
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkt)
# else
          text=vname(2,indxAkt)
# endif
          varid = nccreate(ncid,vname(1,indxAkt),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkt))
        endif
# ifdef SALINITY
! Vertical diffusion coefficient for salinity.

        if (wrt(indxAks) .and. ierr == nf90_noerr) then
#  ifdef AVRH
          text='averaged '/ /vname(2,indxAks)
#  else
          text=vname(2,indxAks)
#  endif
          varid = nccreate(ncid,vname(1,indxAks),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAks))
        endif
# endif /* SALINITY */
# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        if (wrt(indxHbls) .and. ierr == nf90_noerr) then
#  ifdef AVRH
          text='averaged '/ /vname(2,indxHbls)
#  else
          text=vname(2,indxHbls)
#  endif
          varid = nccreate(ncid,vname(1,indxHbls),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbls))
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        if (wrt(indxHbbl) .and. ierr == nf90_noerr) then
#  ifdef AVRH
          text='averaged '/ /vname(2,indxHbbl)
#  else
          text=vname(2,indxHbbl)
#  endif
          varid = nccreate(ncid,vname(1,indxHbbl),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text)
          ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbbl))
        endif
# endif
#endif /* SOLVE3D */


#ifdef WEC
# ifndef AVRH /* NOT AVRH */
!        call def_vars_wec( nchis, .false. )                ! define WEC variables
# else /* AVRH */
        call def_vars_wec( ncavg, .true.  )                ! define WEC averaged variables                                                           )
# endif /* END AVRH */
#endif   /* WEC */

#if defined BIOLOGY_BEC2 && defined BEC2_DIAG
# ifndef AVRH /* NOT AVRH */

        call def_his_bgc_diag                              ! define bgc BEC2_DIAG his variables
# else /* AVRH */
        call def_avg_bgc_diag                              ! define bgc BEC2_DIAG averaged variables
# endif /* END AVRH */
#endif   /* BEC2_DIAG */


! Leave definition mode.
! ----- ---------- -----

        ierr=nf90_enddef(ncid)
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
        ierr=nf90_open(fname(1:lfnm), nf90_write, ncid)
        if (ierr == nf90_noerr) then
          ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)
          if (ierr == nf90_noerr) then
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
     &         'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &          MYID
            goto 99
          endif

          ierr=checkdims(ncid, fname, rec)
          if (ierr == nf90_noerr) then
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
            ierr=nf90_noerr
          endif
        endif

        if (ierr /= nf90_noerr) then
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

! DevinD: no longer need to set var id as handled by read_write.F module


! ---------------------
! Leave definition mode
! ---------------------

        if (ierr == nf90_noerr) then
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
        ierr=nf90_open(fname(1:lfnm), nf_write, ncid)
        if (ierr == nf90_noerr) then
          ierr=nf_set_fill(ncid, nf_nofill, prev_fill_mode)
          if (ierr == nf90_noerr) then
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
     &         'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &          MYID
          endif
        else
# ifdef AVRH
          write(*,'(/1x,6A,I4/)')  '### ERROR: def_avg :: Cannot ',
# else
          write(*,'(/1x,6A,I4/)')  '### ERROR: def_his :: Cannot ',
# endif
     &                      'reopen file ''', fname(1:lfnm), '''.',
     &                       nf90_strerror(ierr) MYID
        endif
#endif
      endif             !<--  create_new_file
      if (ierr /= nf90_noerr) goto 99

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

#ifndef AVRH
# define AVRH
# ifdef AVERAGES
#  include "def_his.F"
# endif
#endif
