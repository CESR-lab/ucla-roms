#include "cppdefs.opt"

      subroutine def_rst(total_rec, ierr)  ! Create or open  ![
                                           ! restart netCDF file.

      use param
      use tracers,         only: t_vname, t_lname, t_units, wrt_t
      use roms_read_write, only: output_root_name
      use ncvars
      use scalars
      use netcdf
      use grid                                                       ! def_grid
      use nc_read_write

      implicit none

      ! input
      integer :: total_rec, ierr
      ! local
      logical :: create_new_file
      integer :: lenstr, checkdims, my_nf_def_dim, varid
      integer :: rec, lrst, lvar, timedim, prev_fill_mode, itrc
#ifdef MASK_LAND_DATA
# ifdef DBLEPREC
      real(kind=8),parameter :: spv_set=1.D+33
# else
      real(kind=4),parameter :: spv_set=1.E+33
# endif
#endif

! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'rst_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.

      ierr=0
      rstname=trim(output_root_name) / / '_rst.nc'
      lrst=lenstr(rstname)
#if defined MPI && defined PARALLEL_FILES
        ! Insert MPI node numbers to file name
        call insert_node(rstname, lrst, mynode, NNODES, ierr)  ! previously done in ncdf_create_file:
#endif
      lrst=lenstr(rstname)
      if (nrpfrst > 0) then
        lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
        call insert_time_index (rstname, lrst, lvar, ierr)
        if (ierr /= 0) goto 99
      endif

! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncrst==-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just reopen, no checking, all ids are assumed
!        .... to be already known (this branch is only for
!        .... output into single file from MPI code).
!     endif
!
! which is designed to implement a flexible opening policy: if
!  (i) ldefhis=.true., it forces creation of a new file [if the file
!      already exists, it will be overwritten]; on the other hand, if
! (ii) ldefhis=.false., it is assumed that the file already exists
!      and an attempt to open it is made first; if the attempt is
!      successful, the file is prepared for appending hew data;
!      if it fails, a new file is created.

      create_new_file=ldefhis
      if (ncrst /= -1) create_new_file=.false.
#if defined MPI && !defined PARALLEL_FILES
      if (mynode > 0) create_new_file=.false.
#endif

  10  if (create_new_file) then
        ierr=nf90_create(rstname(1:lrst), nf90_netcdf4, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf90_set_fill(ncrst, nf90_nofill, prev_fill_mode)
          if (ierr /= nf90_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_rst :: Cannot ',
     &               'create restart NetCDF file ''', rstname(1:lrst),
     &                                       '''',  nf90_strerror(ierr)
          goto 99              !--> ERROR
        endif

        call put_global_atts (ncrst, ierr)

#ifdef PUT_GRID_INTO_RESTART
        if (total_rec == 1) call def_grid(ncrst)
#endif

! Define evolving model variables:
! ------ -------- ----- ----------
! Time step number and time record numbers:
        varid = nccreate(ncrst,'time_step',(/'auxil','time'/),(/iaux,0/),nf90_int)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &             'time step and record numbers from initialization')

! Time.
        varid = nccreate(ncrst,vname(1,indxTime),(/'time'/),(/0/),nf90_double)
        ierr=nf90_put_att(ncrst, varid, 'long_name', vname(2,indxTime))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxTime))

! Free-surface.
        varid = nccreate(ncrst,vname(1,indxZ),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name', vname(2,indxZ))
        ierr=nf90_put_att (ncrst, varid, 'units', vname(3,indxZ))
        ierr = nf90_put_att(ncrst,varid,'_FillValue',spv_set)

! 2D momenta in XI- and ETA-components
        varid = nccreate(ncrst,vname(1,indxUb),(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                                   vname(2,indxUb))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxUb))
        ierr = nf90_put_att(ncrst,varid,'_FillValue',spv_set)
        varid = nccreate(ncrst,vname(1,indxVb),(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                                   vname(2,indxVb))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxVb))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        varid = nccreate(ncrst,'DU_avg2',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &           '<<fast-time averaged ubar(:,:,n+1/2)>>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DV_avg2',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &           '<<fast-time-averaged vbar(:,:,n+1/2)>>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DU_avg_bak',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged ubar(:,:,n-1)>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DV_avg_bak',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged vbar(:,:,n-1)>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
#  elif defined PRED_COUPLED_MODE
        varid = nccreate(ncrst,'rufrc_bak',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                        '3D-to-2D forcing, XI-component')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
        varid = nccreate(ncrst,'rvfrc_bak',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                       '3D-to-2D forcing, ETA-component')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
#  endif
# endif  /* EXACT_RESTART */

! 3D momenta in XI- and ETA-directions.
        varid = nccreate(ncrst, vname(1,indxU),(/'xi_u','eta_rho','s_rho','time'/),(/xi_u,eta_rho,N,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                     vname(2,indxU))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                     vname(3,indxU))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,vname(1,indxV),(/'xi_rho','eta_v','s_rho','time'/),(/xi_rho,eta_v,N,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                     vname(2,indxV))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                     vname(3,indxV))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

! Tracer variables.
        do itrc=1,NT
          varid = nccreate(ncrst,t_vname(itrc),(/'xi_rho','eta_rho','s_rho','time'/),
     &                                         (/xi_rho,eta_rho,N,0/), NF_FTYPE)
          ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                       t_lname(itrc))
          ierr=nf90_put_att (ncrst, varid, 'units', t_units(itrc))
          ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.
!        ierr=nf_def_var (ncrst, vname(1,indxHbls)(1:lvar),
!     &                            NF_FTYPE, 3, r2dgrd, rstHbls)
        varid = nccreate(ncrst,vname(1,indxHbls),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                 vname(2,indxHbls))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                 vname(3,indxHbls))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.
!        ierr=nf_def_var (ncrst, vname(1,indxHbbl),
!     &                            NF_FTYPE, 3, r2dgrd, rstHbbl)
        varid = nccreate(ncrst,vname(1,indxHbbl),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                 vname(2,indxHbbl))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                 vname(3,indxHbbl))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
# endif
#endif /* SOLVE3D */

        ierr=nf90_enddef(ncrst)
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
          write(*,'(7x,4A,1x,A,i4)') 'def_rst :: Created new ',
     &          'netCDF file ''', rstname(1:lrst), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Check consistency of the dimensions of fields from the file with
! model dimensions. Determine the current size of unlimited dimension
! and set initial record [in the case of MPI serialized output, at
! this moment the last time record is assumed to be **partially**
! written by MPI processes with lower rank. Thus the next write is
! expected to be into the same record rather than next one (except
! MPI-master, who initializes the record).
!
! In the case when file is rejected (whether it cannot be opened, or
! something is wrong with its dimensions, create new file.

      elseif (ncrst == -1) then
        ierr=nf90_open(rstname, nf90_write, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf90_set_fill(ncrst, nf90_nofill, prev_fill_mode)
          if (ierr /= nf90_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif

          ierr=checkdims(ncrst, rstname, rec)
          if (ierr == nf90_noerr) then
            if (nrpfrst==0) then
              ierr=rec+1 - nrecrst
            else
              ierr=rec+1 - (1+mod(nrecrst-1, abs(nrpfrst)))
            endif
            if (ierr > 0) then
              mpi_master_only write(*,
     &           '(/1x,A,I5,1x,3A/21x,A,I5,1x,A,I5,1x,A/21x,A/)')
     &           'WARNING: def_rst :: The actual number of records',
     &            rec, 'in netCDF file ''',  rstname(1:lrst), '''',
     &           'exceeds record', rec+1-ierr, '/', total_rec,
     &           'for restart initial conditions.',
     &           'All records beyond this will be overwriding.'
              rec=rec-ierr
            elseif (nrpfrst==0) then
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
          if (mynode==0) then
            create_new_file=.true.
            goto 10
          else
            write(*,'(/1x,4A, 1x,A,I4/)')  '### ERROR: def_rst :: ',
     &        'Cannot open restart netCDF file ''', rstname(1:lrst),
     &        '''.' MYID
            goto 99                                     !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif

! Find netCDF variable IDs of evolving model variables:
! ---- ------ -------- --- -- -------- ----- ----------

! Time step indices and time:

        ierr=nf90_inq_varid(ncrst, 'time_step', rstTstep)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxTime))
          ierr=nf90_inq_varid(ncrst, vname(1,indxTime)(1:lvar), rstTime)
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxTime)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) 'time_step', rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

! Free-surface and 2D momenta in XI- and ETA-directions

        lvar=lenstr(vname(1,indxZ))
        ierr=nf90_inq_varid(ncrst, vname(1,indxZ)(1:lvar), rstZ)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf90_inq_varid(ncrst, vname(1,indxUb)(1:lvar), rstUb)
          if (ierr == nf90_noerr) then
            lvar=lenstr(vname(1,indxVb))
            ierr=nf90_inq_varid(ncrst, vname(1,indxVb)(1:lvar), rstVb)
            if (ierr /= nf90_noerr) then
              write(*,1) vname(1,indxVb)(1:lvar), rstname(1:lrst)
            endif
          else
            write(*,1) vname(1,indxUb)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxZ)(1:lvar),  rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        ierr=nf90_inq_varid(ncrst,   'DU_avg2',   rst_DU_avg2)
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncrst,  'DV_avg2',     rst_DV_avg2)
          if (ierr == nf90_noerr) then
            ierr=nf90_inq_varid(ncrst, 'DU_avg_bak', rst_DU_avg_bak)
            if (ierr == nf90_noerr) then
              ierr=nf90_inq_varid(ncrst, 'DV_avg_bak', rst_DV_avg_bak)
              if (ierr /= nf90_noerr) then
                write(*,1) 'DV_avg_bak', rstname(1:lrst)
              endif
            else
              write(*,1) 'DU_avg_bak', rstname(1:lrst)
            endif
          else
            write(*,1) 'DV_avg2', rstname(1:lrst)
          endif
        else
          write(*,1) 'DU_avg2', rstname(1:lrst)
        endif
#  elif defined PRED_COUPLED_MODE
        ierr=nf90_inq_varid(ncrst, 'rufrc_bak', rst_rufrc)
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncrst, 'rvfrc_bak', rst_rvfrc)
          if (ierr /= nf90_noerr) then
            write(*,1) 'rvfrc_bak', rstname(1:lrst)
          endif
        else
          write(*,1) 'rufrc_bak', rstname(1:lrst)
        endif
#  endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR
# endif

! 3D momenta n XI- and ETA-directions.

        lvar=lenstr(vname(1,indxU))
        ierr=nf90_inq_varid (ncrst, vname(1,indxU)(1:lvar), rstU)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf90_inq_varid (ncrst, vname(1,indxV)(1:lvar), rstV)
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxV)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxU)(1:lvar), rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

! Tracer variables.

        do itrc=1,NT
          lvar=lenstr(vname(1,indxT+itrc-1))
          ierr=nf90_inq_varid(ncrst, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  rstT(itrc))
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxT+itrc-1)(1:lvar), rstname(1:lrst)
            goto 99                                       !--> ERROR
          endif
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        lvar=lenstr(vname(1,indxHbls))
        ierr=nf90_inq_varid(ncrst, vname(1,indxHbls)(1:lvar), rstHbls)
        if (ierr /= nf90_noerr) then
          write(*,1) vname(1,indxHbls)(1:lvar), rstname(1:lrst)
          goto 99                                        !--> ERROR
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.
        ierr=nf90_inq_varid(ncrst, vname(1,indxHbbl), rstHbbl)
        if (ierr /= nf90_noerr) then
        write(*,1) vname(1,indxHbbl), rstname
          goto 99                                       !--> ERROR
        endif
# endif
#endif
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(*,'(7x,4A,I4,2x,A,I4)') 'def_rst :: open existing ',
     &            'file ''', rstname(1:lrst), ''' from rec =', rec
     &             MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

#if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf90_open(rstname(1:lrst), nf_write, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf_set_fill(ncrst, nf_nofill, prev_fill_mode)
          if (ierr == nf90_noerr) then
# ifndef MPI_SILENT_MODE
            write(*,'(7x,4A,I4)') 'def_rst :: reopen existing file ''',
     &              rstname(1:lrst), '''.' MYID
# endif
          else
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A,1x,A,I4/)') '### ERROR: def_rst :: Cannot ',
     &            'reopen file ''', rstname(1:lrst),'''.'
     &             MYID
          goto 99                                         !--> ERROR
        endif
#endif
      endif              !<-- create_new_file

#ifdef PUT_GRID_INTO_RESTART

! Write grid variables.
! ----- ---- ----------

      if (total_rec == 1) call wrt_grid(ncrst, rstname, lrst)
#endif
  99  return                                              !--> ERROR
   1  format(/1x,'### ERROR: def_rst :: Cannot find variable ''', A,
     &                            ''' in netCDF file ''', A, '''.'/)

      end subroutine def_rst  !]
