#include "cppdefs.h"

      subroutine def_rst(total_rec, ierr)  ! Create or open
                                           ! restart netCDF file.
      use bgc ! def_rst_bgc_diag
      use tracers, only: t_vname, t_lname, t_units, wrt_t

      implicit none
      integer total_rec, ierr
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"
      logical create_new_file
      integer lenstr, checkdims, my_nf_def_dim
      integer rec, lrst, lvar, timedim, prev_fill_mode,
     &        r2dgrd(3), u2dgrd(3), v2dgrd(3),  auxil(2)
#ifdef SOLVE3D
     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
#endif
#ifdef MASK_LAND_DATA
# ifdef DBLEPREC
      real*8, parameter :: spv_set=1.D+33
# else
      real*4, parameter :: spv_set=1.E+33
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

      ierr=0 ; lrst=lenstr(rstname)
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

! Create new restart file:    Put global attributes
!======= === ======= =====    and define all variables.

  10  if (create_new_file) then
        ierr=nf_create(rstname(1:lrst), nf_clobber+nf_64bit_offset,
     &                                                        ncrst)
        if (ierr == nf_noerr) then
          ierr=nf_set_fill(ncrst, nf_nofill, prev_fill_mode)
          if (ierr /= nf_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_rst :: Cannot ',
     &               'create restart NetCDF file ''', rstname(1:lrst),
     &                                       '''',  nf_strerror(ierr)
          goto 99              !--> ERROR
        endif

! Put global attributes:
! --- ------ -----------

        call put_global_atts (ncrst, ierr)

! Define dimensions for field placement on staggered grids:
! ------ ---------- --- ----- --------- -- --------- -------

        ierr=my_nf_def_dim(ncrst, 'xi_rho',  xi_rho,  r2dgrd(1))
        ierr=my_nf_def_dim(ncrst, 'xi_u',    xi_u,    u2dgrd(1))
        ierr=my_nf_def_dim(ncrst, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=my_nf_def_dim(ncrst, 'eta_v',   eta_v,   v2dgrd(2))
#ifdef SOLVE3D
        ierr=my_nf_def_dim(ncrst, 's_rho',   N,       r3dgrd(3))
        ierr=my_nf_def_dim(ncrst, 's_w',     N+1,     w3dgrd(3))
#endif
        ierr=my_nf_def_dim(ncrst, 'time', nf_unlimited, timedim)
        ierr=my_nf_def_dim(ncrst, 'auxil',   iaux,     auxil(1))

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


#ifdef PUT_GRID_INTO_RESTART

! Define grid variables.
! ------ ---- ----------
        if (total_rec == 1) call def_grid(ncrst, r2dgrd)
#endif

! Define evolving model variables:
! ------ -------- ----- ----------
! Time step number and time record numbers:

        ierr=nf_def_var(ncrst, 'time_step', nf_int, 2, auxil,
     &                                                 rstTstep)
        ierr=nf_put_att_text (ncrst, rstTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')

! Time.

        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var(ncrst, vname(1,indxTime)(1:lvar), NF_FTYPE,
     &                                        1, timedim, rstTime)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncrst, rstTime, 'long_name', lvar,
     &                                   vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncrst, rstTime, 'units',     lvar,
     &                                   vname(3,indxTime)(1:lvar))

! Free-surface.

        lvar=lenstr(vname(1,indxZ))
        ierr=nf_def_var(ncrst, vname(1,indxZ)(1:lvar), NF_FTYPE,
     &                                             3, r2dgrd,rstZ)
        lvar=lenstr(vname(2,indxZ))
        ierr=nf_put_att_text (ncrst, rstZ, 'long_name', lvar,
     &                                      vname(2,indxZ)(1:lvar))
        lvar=lenstr(vname(3,indxZ))
        ierr=nf_put_att_text (ncrst, rstZ, 'units', lvar,
     &                                      vname(3,indxZ)(1:lvar))
#ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE (ncrst, rstZ, '_FillValue', NF_FTYPE,
     &                                                 1, spv_set)
#endif

! 2D momenta in XI- and ETA-components

        lvar=lenstr(vname(1,indxUb))
        ierr=nf_def_var(ncrst, vname(1,indxUb)(1:lvar), NF_FTYPE,
     &                                          3, u2dgrd, rstUb)
        lvar=lenstr(vname(2,indxUb))
        ierr=nf_put_att_text(ncrst, rstUb, 'long_name', lvar,
     &                                   vname(2,indxUb)(1:lvar))
        lvar=lenstr(vname(3,indxUb))
        ierr=nf_put_att_text(ncrst, rstUb, 'units', lvar,
     &                                   vname(3,indxUb)(1:lvar))
#ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rstUb, '_FillValue',
     &                                     NF_FTYPE, 1, spv_set)
#endif

        lvar=lenstr(vname(1,indxVb))
        ierr=nf_def_var(ncrst, vname(1,indxVb)(1:lvar), NF_FTYPE,
     &                                          3, v2dgrd, rstVb)
        lvar=lenstr(vname(2,indxVb))
        ierr=nf_put_att_text(ncrst, rstVb, 'long_name', lvar,
     &                                   vname(2,indxVb)(1:lvar))
        lvar=lenstr(vname(3,indxVb))
        ierr=nf_put_att_text(ncrst, rstVb, 'units',     lvar,
     &                                   vname(3,indxVb)(1:lvar))
#ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rstVb, '_FillValue',
     &                                     NF_FTYPE, 1, spv_set)
#endif

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        ierr=nf_def_var(ncrst, 'DU_avg2', NF_FTYPE, 3, u2dgrd,
     &                                             rst_DU_avg2)
        ierr=nf_put_att_text(ncrst, rst_DU_avg2, 'long_name',
     &            38, '<<fast-time averaged ubar(:,:,n+1/2)>>')
#   ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rst_DU_avg2, '_FillValue',
     &                                    NF_FTYPE, 1, spv_set)
#   endif

        ierr=nf_def_var(ncrst, 'DV_avg2', NF_FTYPE, 3, v2dgrd,
     &                                             rst_DV_avg2)
        ierr=nf_put_att_text(ncrst, rst_DV_avg2, 'long_name',
     &            38, '<<fast-time-averaged vbar(:,:,n+1/2)>>')
#   ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rst_DV_avg2, '_FillValue',
     &                                  NF_FTYPE,  1, spv_set)
#   endif

        ierr=nf_def_var(ncrst, 'DU_avg_bak', NF_FTYPE, 3, u2dgrd,
     &                                             rst_DU_avg_bak)
        ierr=nf_put_att_text(ncrst, rst_DU_avg_bak, 'long_name',
     &   50, '<back-step mixed fast-time-averaged ubar(:,:,n-1)>')
#   ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rst_DU_avg_bak, '_FillValue',
     &                                     NF_FTYPE,  1, spv_set)
#   endif

        ierr=nf_def_var(ncrst, 'DV_avg_bak', NF_FTYPE, 3, v2dgrd,
     &                                             rst_DV_avg_bak)
        ierr=nf_put_att_text(ncrst, rst_DV_avg_bak, 'long_name',
     &   50, '<back-step mixed fast-time-averaged vbar(:,:,n-1)>')
#   ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rst_DV_avg_bak, '_FillValue',
     &                                     NF_FTYPE, 1,  spv_set)
#   endif
#  elif defined PRED_COUPLED_MODE
        ierr=nf_def_var(ncrst, 'rufrc_bak', NF_FTYPE, 3, u2dgrd,
     &                                                rst_rufrc)
        ierr=nf_put_att_text(ncrst, rst_rufrc, 'long_name', 30,
     &                        '3D-to-2D forcing, XI-component')
#   ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rst_rufrc, '_FillValue',
     &                                NF_FTYPE, 1,  spv_set)
#   endif
        ierr=nf_def_var(ncrst, 'rvfrc_bak', NF_FTYPE, 3, v2dgrd,
     &                                                rst_rvfrc)
        ierr=nf_put_att_text(ncrst, rst_rvfrc, 'long_name', 31,
     &                       '3D-to-2D forcing, ETA-component')
#   ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE(ncrst, rst_rvfrc, '_FillValue',
     &                                NF_FTYPE, 1,  spv_set)
#   endif
#  endif
# endif  /* EXACT_RESTART */

! 3D momenta in XI- and ETA-directions.

        lvar=lenstr(vname(1,indxU))
        ierr=nf_def_var (ncrst, vname(1,indxU)(1:lvar), NF_FTYPE,
     &                                           4, u3dgrd, rstU)
        lvar=lenstr(vname(2,indxU))
        ierr=nf_put_att_text (ncrst, rstU, 'long_name', lvar,
     &                                     vname(2,indxU)(1:lvar))
        lvar=lenstr(vname(3,indxU))
        ierr=nf_put_att_text (ncrst, rstU, 'units',     lvar,
     &                                     vname(3,indxU)(1:lvar))
# ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE (ncrst, rstU, '_FillValue', NF_FTYPE,
     &                                                 1, spv_set)
# endif

        lvar=lenstr(vname(1,indxV))
        ierr=nf_def_var (ncrst, vname(1,indxV)(1:lvar), NF_FTYPE,
     &                                           4, v3dgrd, rstV)
        lvar=lenstr(vname(2,indxV))
        ierr=nf_put_att_text (ncrst, rstV, 'long_name', lvar,
     &                                     vname(2,indxV)(1:lvar))
        lvar=lenstr(vname(3,indxV))
        ierr=nf_put_att_text (ncrst, rstV, 'units',     lvar,
     &                                     vname(3,indxV)(1:lvar))
# ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE (ncrst, rstV, '_FillValue', NF_FTYPE,
     &                                                 1, spv_set)
# endif
! Tracer variables.

        do itrc=1,NT
          lvar=lenstr(t_vname(itrc))
          ierr=nf_def_var (ncrst, t_vname(itrc)(1:lvar),
     &                          NF_FTYPE, 4, r3dgrd, rstT(itrc))
          lvar=lenstr(t_lname(itrc))
          ierr=nf_put_att_text (ncrst, rstT(itrc), 'long_name',
     &                     lvar,   t_lname(itrc)(1:lvar))
          lvar=lenstr(vname(3,indxT+itrc-1))
          ierr=nf_put_att_text (ncrst, rstT(itrc), 'units', lvar,
     &                             t_units(itrc)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FTYPE (ncrst, rstT(itrc), '_FillValue',
     &                                     NF_FTYPE, 1, spv_set)
# endif
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        lvar=lenstr(vname(1,indxHbls))
        ierr=nf_def_var (ncrst, vname(1,indxHbls)(1:lvar),
     &                            NF_FTYPE, 3, r2dgrd, rstHbls)
        lvar=lenstr(vname(2,indxHbls))
        ierr=nf_put_att_text (ncrst, rstHbls, 'long_name', lvar,
     &                                 vname(2,indxHbls)(1:lvar))
        lvar=lenstr(vname(3,indxHbls))
        ierr=nf_put_att_text (ncrst, rstHbls, 'units', lvar,
     &                                 vname(3,indxHbls)(1:lvar))
#  ifdef MASK_LAND_DATA
        ierr=nf_put_att_FTYPE (ncrst, rstHbls, '_FillValue',
     &                                  NF_FTYPE, 1, spv_set)
#  endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        lvar=lenstr(vname(1,indxHbbl))
        ierr=nf_def_var (ncrst, vname(1,indxHbbl)(1:lvar),
     &                            NF_FTYPE, 3, r2dgrd, rstHbbl)
        lvar=lenstr(vname(2,indxHbbl))
        ierr=nf_put_att_text (ncrst, rstHbbl, 'long_name', lvar,
     &                                 vname(2,indxHbbl)(1:lvar))
        lvar=lenstr(vname(3,indxHbbl))
        ierr=nf_put_att_text (ncrst, rstHbbl, 'units',     lvar,
     &                                 vname(3,indxHbbl)(1:lvar))
#  ifdef MASK_LAND_DATA
        ierr=nf_put_att_FOUT (ncrst, rstHbbl, '_FillValue',
     &                                 NF_FTYPE, 1, spv_set)
#  endif
# endif
#endif /* SOLVE3D */

      ! Define BGC diag variables
      ! -------------------------
#if defined BIOLOGY_BEC2 && defined BEC2_DIAG
        call def_rst_bgc_diag(r2dgrd
# ifdef SOLVE3D
     &                      , r3dgrd
# endif
     &                               )
#endif

! Leave definition mode.                  Also initialize record
! ----- ---------- -----                  dimension size to zero.

        ierr=nf_enddef(ncrst)
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
        ierr=nf_open(rstname(1:lrst), nf_write, ncrst)
        if (ierr == nf_noerr) then
          ierr=nf_set_fill(ncrst, nf_nofill, prev_fill_mode)
          if (ierr /= nf_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &         MYID
            goto 99
          endif

          ierr=checkdims(ncrst, rstname, rec)
          if (ierr == nf_noerr) then
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
            ierr=nf_noerr
          endif
        endif

        if (ierr /= nf_noerr) then
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

        ierr=nf_inq_varid(ncrst, 'time_step', rstTstep)
        if (ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxTime))
          ierr=nf_inq_varid(ncrst, vname(1,indxTime)(1:lvar), rstTime)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxTime)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) 'time_step', rstname(1:lrst)
        endif
        if (ierr /= nf_noerr) goto 99                     !--> ERROR

! Free-surface and 2D momenta in XI- and ETA-directions

        lvar=lenstr(vname(1,indxZ))
        ierr=nf_inq_varid(ncrst, vname(1,indxZ)(1:lvar), rstZ)
        if (ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid(ncrst, vname(1,indxUb)(1:lvar), rstUb)
          if (ierr == nf_noerr) then
            lvar=lenstr(vname(1,indxVb))
            ierr=nf_inq_varid(ncrst, vname(1,indxVb)(1:lvar), rstVb)
            if (ierr /= nf_noerr) then
              write(*,1) vname(1,indxVb)(1:lvar), rstname(1:lrst)
            endif
          else
            write(*,1) vname(1,indxUb)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxZ)(1:lvar),  rstname(1:lrst)
        endif
        if (ierr /= nf_noerr) goto 99                     !--> ERROR

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        ierr=nf_inq_varid(ncrst,   'DU_avg2',   rst_DU_avg2)
        if (ierr == nf_noerr) then
          ierr=nf_inq_varid(ncrst,  'DV_avg2',     rst_DV_avg2)
          if (ierr == nf_noerr) then
            ierr=nf_inq_varid(ncrst, 'DU_avg_bak', rst_DU_avg_bak)
            if (ierr == nf_noerr) then
              ierr=nf_inq_varid(ncrst, 'DV_avg_bak', rst_DV_avg_bak)
              if (ierr /= nf_noerr) then
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
        ierr=nf_inq_varid(ncrst, 'rufrc_bak', rst_rufrc)
        if (ierr == nf_noerr) then
          ierr=nf_inq_varid(ncrst, 'rvfrc_bak', rst_rvfrc)
          if (ierr /= nf_noerr) then
            write(*,1) 'rvfrc_bak', rstname(1:lrst)
          endif
        else
          write(*,1) 'rufrc_bak', rstname(1:lrst)
        endif
#  endif
        if (ierr /= nf_noerr) goto 99                     !--> ERROR
# endif

! 3D momenta n XI- and ETA-directions.

        lvar=lenstr(vname(1,indxU))
        ierr=nf_inq_varid (ncrst, vname(1,indxU)(1:lvar), rstU)
        if (ierr == nf_noerr) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncrst, vname(1,indxV)(1:lvar), rstV)
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxV)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxU)(1:lvar), rstname(1:lrst)
        endif
        if (ierr /= nf_noerr) goto 99                     !--> ERROR

! Tracer variables.

        do itrc=1,NT
          lvar=lenstr(vname(1,indxT+itrc-1))
          ierr=nf_inq_varid(ncrst, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  rstT(itrc))
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxT+itrc-1)(1:lvar), rstname(1:lrst)
            goto 99                                       !--> ERROR
          endif
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        lvar=lenstr(vname(1,indxHbls))
        ierr=nf_inq_varid(ncrst, vname(1,indxHbls)(1:lvar), rstHbls)
        if (ierr /= nf_noerr) then
          write(*,1) vname(1,indxHbls)(1:lvar), rstname(1:lrst)
          goto 99                                        !--> ERROR
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        lvar=lenstr(vname(1,indxHbbl))
        ierr=nf_inq_varid(ncrst, vname(1,indxHbbl)(1:lvar), rstHbbl)
        if (ierr /= nf_noerr) then
        write(*,1) vname(1,indxHbbl)(1:lvar), rstname(1:lrst)
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
        ierr=nf_open(rstname(1:lrst), nf_write, ncrst)
        if (ierr == nf_noerr) then
          ierr=nf_set_fill(ncrst, nf_nofill, prev_fill_mode)
          if (ierr == nf_noerr) then
# ifndef MPI_SILENT_MODE
            write(*,'(7x,4A,I4)') 'def_rst :: reopen existing file ''',
     &              rstname(1:lrst), '''.' MYID
# endif
          else
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
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
      end
