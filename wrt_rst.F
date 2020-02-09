#include "cppdefs.h"
                              ! Write model prognostic
      subroutine wrt_rst      ! variables into restart
                              ! netCDF file.
      implicit none
#include "param.h"
#include "ncvars.h"
      integer record, start(2), count(2), ibuff(iaux),
     &        i, ierr, lstr, lvar, lenstr, ncdf_write
#ifdef SOLVE3D
     &      , itrc
#endif
      character(len=18) tstring
#include "scalars.h"
#include "ocean2d.h"
#ifdef EXACT_RESTART
# include "coupling.h"
#endif
#include "ocean3d.h"
#if defined LMD_KPP || defined LMD_BKPP
# include "mixing.h"
#endif

#include "netcdf.inc"

#if defined MPI && !defined PARALLEL_FILES
# include "mpif.h"
      integer status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 1, ocean_grid_comm, status, ierr)
#endif

! Create/open restart and set record index within the file.

      call def_rst(nrecrst, ierr)
      if (ierr /= nf_noerr) goto 99
      lstr=lenstr(rstname)
      nrecrst=max(nrecrst,1)                  !!! WARNING: Here it
      if (nrpfrst == 0) then                  !!! is assumed that
        record=nrecrst                        !!! global restart
      else                                    !!! record index
        record=1+mod(nrecrst-1, abs(nrpfrst)) !!! nrecrst is already
      endif                                   !!! advanced by main.

! Write out evolving model variables:
! ----- --- -------- ----- ----------

      ibuff(1)=iic-1  ; ibuff(2)=nrecrst        ! Time step and
      ibuff(4:iaux)=0 ; ibuff(3)=nrechis        ! record numbers.
#ifdef AVERAGES
      ibuff(4)=nrecavg
#endif
      start(1)=1      ; count(1)=iaux
      start(2)=record ; count(2)=1

      ierr=nf_put_vara_int (ncrst, rstTstep, start, count, ibuff)
      if (ierr /= nf_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_rst :: ',
     &       'Cannot write variable ''time_step'' into restart file, ',
     &       'rec =', record, nf_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif

! Time: Note that (1) time is always double precision; and (2) which
! one to write into the file, time or tdays, is decided by attribute
! "units" stored as vname(3,indxTime).  When computng time is always
! in seconds, however it is more convenient to have it in days in all
! the files; (3) this code is designed for politeness - it reads time
! units from the initial condition file (see "get_init.F") and uses
! it throughout the entire run.

      if (vname(3,indxTime)(1:6) == 'second') then
        ierr=nf_put_var1_double(ncrst, rstTime, record, time)
      else
        ierr=nf_put_var1_double(ncrst, rstTime, record, tdays)
      endif
      if (ierr /= nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,'(/1x,4A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_rst :: ',
     &        'Cannot write variable ''',    vname(1,indxTime)(1:lvar),
     &        ''' into restart file, rec =', record, nf_strerror(ierr)
     &            MYID
        goto 99                                           !--> ERROR
      endif

! Free-surface and 2D momentum components in XI- and ETA-directions.

      ierr=ncdf_write(ncrst,  rstZ,  record,  r_var,
     &                zeta(START_2D_ARRAY,knew), 1)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncrst,  rstUb,  record, u_var,
     &                  ubar(START_2D_ARRAY,knew), 1)
        if (ierr == nf_noerr) then
          ierr=ncdf_write(ncrst,  rstVb,  record, v_var,
     &                    vbar(START_2D_ARRAY,knew), 1)
          if (ierr /= nf_noerr) then
            lvar=lenstr(vname(1,indxVb))
            write(*,1) vname(1,indxVb)(1:lvar), record MYID
          endif
        else
          lvar=lenstr(vname(1,indxUb))
          write(*,1) vname(1,indxUb)(1:lvar), record MYID
        endif
      else
        lvar=lenstr(vname(1,indxZ))
        write(*,1) vname(1,indxZ)(1:lvar), record MYID
      endif
      if (ierr /= nf_noerr) goto 99                !--> ERROR

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
      ierr=ncdf_write(ncrst, rst_DU_avg2, record, u_var, DU_avg2, 1)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncrst, rst_DV_avg2, record, v_var, DV_avg2, 1)
        if (ierr == nf_noerr) then
          ierr=ncdf_write(ncrst, rst_DU_avg_bak, record, u_var,
     &                                                  DU_avg_bak, 1)
          if (ierr == nf_noerr) then
            ierr=ncdf_write(ncrst, rst_DV_avg_bak, record, v_var,
     &                                                  DV_avg_bak, 1)
            if (ierr /= nf_noerr) then
              write(*,1) 'DV_avg_bak', record MYID
            endif
          else
            write(*,1) 'DU_avg_bak', record MYID
          endif
        else
          write(*,1) 'DV_avg2', record MYID
        endif
      else
        write(*,1) 'DU_avg2', record MYID
      endif
#  elif defined PRED_COUPLED_MODE
      ierr=ncdf_write(ncrst,   rst_rufrc,  record,  u_var,
     &                  rufrc_bak(START_2D_ARRAY,nstp), 1)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncrst, rst_rvfrc,  record,  v_var,
     &                  rvfrc_bak(START_2D_ARRAY,nstp), 1)
        if (ierr /= nf_noerr) then
          write(*,1) 'rvfrc_bak', record MYID
        endif
      else
        write(*,1) 'rufrc_bak', record MYID
      endif
#  endif
      if (ierr /= nf_noerr) goto 99                !--> ERROR
# endif

! 3D momentum components in XI- and ETA-directions.

      ierr=ncdf_write(ncrst, rstU,  record,  u_var,
     &                 u(START_2D_ARRAY,1,nstp), N)
      if (ierr == nf_noerr) then
        ierr=ncdf_write(ncrst, rstV, record, v_var,
     &                 v(START_2D_ARRAY,1,nstp), N)
        if (ierr /= nf_noerr) then
          lvar=lenstr(vname(1,indxV))
          write(*,1) vname(1,indxV)(1:lvar), record MYID
        endif
      else
        lvar=lenstr(vname(1,indxU))
        write(*,1) vname(1,indxU)(1:lvar), record MYID
      endif
      if (ierr /= nf_noerr) goto 99                !--> ERROR

! Tracer variables.

      do itrc=1,NT
        ierr=ncdf_write(ncrst,  rstT(itrc), record, r_var,
     &                  t(START_2D_ARRAY,1,nstp,itrc), N)
        if (ierr /= nf_noerr) then
          lvar=lenstr(vname(1,indxT+itrc-1))
          write(*,1) vname(1,indxT+itrc-1)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      enddo

# ifdef LMD_KPP
! Depth of planetary boundary layer [m].

      ierr=ncdf_write(ncrst, rstHbls, record, r_var,
     &                hbls(START_2D_ARRAY,nstp), 1)
      if (ierr /= nf_noerr) then
        lvar=lenstr(vname(1,indxHbls))
        write(*,1) vname(1,indxHbls)(1:lvar), record MYID
        goto 99                                         !--> ERROR
      endif
# endif
# ifdef LMD_BKPP
! Thickness of KPP bottom boundary layer [m].

      ierr=ncdf_write(ncrst, rstHbbl, record, r_var,
     &                hbbl(START_2D_ARRAY,nstp), 1)
      if (ierr /= nf_noerr) then
        lvar=lenstr(vname(1,indxHbbl))
        write(*,1) vname(1,indxHbbl)(1:lvar), record MYID
        goto 99                                         !--> ERROR
      endif
# endif
#endif /* SOLVE3D */


  1   format(/1x, '### ERROR: wrt_rst :: Cannot write variable ''',
     &             A, ''' into restart file, rec =', i6, 3x,A,i4/)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to disk to the allow other
! processes to access data immediately after it is written.

#if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close (ncrst)
      if (nrpfrst > 0 .and. record >= nrpfrst) ncrst=-1
#else
      if (nrpfrst > 0 .and. record >= nrpfrst) then
        ierr=nf_close (ncrst)
        ncrst=-1
      else
        ierr=nf_sync(ncrst)
      endif
#endif
      if (ierr == nf_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(tstring,'(F18.8)') tdays
          i=1
          do while (i<18 .and. tstring(i:i)==' ')
            i=i+1
          enddo
          write(*,'(7x,A,1x,A,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &      'wrt_rst :: wrote restart, tdays =', tstring(i:i+8),
     &      'step =', ibuff(1),  'rec =', record, '/',  nrechis
     &       MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(/1x,2A/)')    '### ERROR: wrt_rst :: Cannot ',
     &                          'synchronize/close restart file.'
        if (may_day_flag == 0) may_day_flag=3
      endif

#if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                        mynode+1, 1, ocean_grid_comm, ierr)
#endif
      end
