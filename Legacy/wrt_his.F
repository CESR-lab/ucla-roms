#include "cppdefs.opt"

      subroutine wrt_his  ![
      ! writes requested model fields at requested levels into history netCDF file.
      !
      ! NOTE: wrt_his is called in middle of roms_step in main.F because
      !       even though u/ubar, etc, were computed at the end of the
      !       previous timestep, other variables like 'rho', 'Akv', etc
      !       have only just been computed.

      use param
      use wec_frc    !, only: wrt_his_wec
      use tracers    !, only: wrt_his_trc
      use bgc        ! wrt_his_bgc_diag
      use read_write ! nc_write_var
      use netcdf     ! for nf90_routines
      use eos_vars
      use grid
      use mixing
      use ncvars
      use ocean_vars
      use scalars
      use work_mod
      use mpi
      use netcdf
      use nc_read_write
      use roms_read_write

      implicit none

      integer :: ierr,  record, i, lstr, lvar, lenstr, var_id_tmp
     &        , start1d(2), count(2), ibuff(iaux), ncdf_write
      integer :: start(4)
      character(len=18) :: tstring
#ifdef SOLVE3D
      integer :: tile, itrc, j,k
#endif
#if defined MPI && !defined PARALLEL_FILES
      integer :: status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif

      call def_his(nrechis, ierr)         ! create/open history and set record index within the file.
      if (ierr /= nf90_noerr) goto 99
      lstr=lenstr(hisname)
      nrechis=max(nrechis,1)              !!! WARNING: Once time
      if (nrpfhis == 0) then              !!! stepping has been
        record=nrechis                    !!! started, it is assumed
      else                                !!! that the global history
        record=1+mod(nrechis-1, nrpfhis)  !!! record index "nrechis"
      endif                               !!! is advanced by main.

      ! write out evolving model variables:
      ibuff(1)=iic-1  ; ibuff(2)=nrecrst    ! Time step and
      ibuff(4:iaux)=0 ; ibuff(3)=nrechis    ! record numbers.
#ifdef AVERAGES
      ibuff(4)=nrecavg
#endif
      start1d(1)=1      ; count(1)=iaux
      start1d(2)=record ; count(2)=1
      start=1

      ierr=nf90_inq_varid(nchis, 'time_step', var_id_tmp)
      ierr=nf90_put_var(nchis, var_id_tmp, ibuff, start1d, count)
      if (ierr /= nf90_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_his :: ',
     &       'Cannot write variable ''time_step'' into history file, ',
     &       'rec =', record, nf90_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif


! Time: Note that (1) time is always double precision; and (2) which
! one to write into the file, time or tdays, is decided by attribute
! "units" stored as vname(3,indxTime).  When computng time is always
! in seconds, however it is more convenient to have it in days in all
! the files.

      ierr=nf90_inq_varid(nchis, 'ocean_time', var_id_tmp)
      ierr=nf90_put_var(nchis, var_id_tmp, time,(/record/))

! Barotropic mode variables: free-surface and 2D momentum
! components in XI-,ETA-directions.
! 'knew' index is normally 'n+1' timestep, but at this point
! its value hasn't evolved from 'n+1' of previous timestep,
! which makes it 'n' of this timestep still.

      start(3)=record
      if (wrthis(indxZ))  call ncwrite(nchis,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start)
      if (wrthis(indxUb)) call ncwrite(nchis,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start)
      if (wrthis(indxVb)) call ncwrite(nchis,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start)

#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_his called at the middle of next timestep as some variables only calculated there for t=n.

      start(3)=1; start(4)=record
      if (wrthis(indxU)) call ncwrite(nchis,vname(1,indxU),u( 1:i1,j0:j1,:,nstp),start)
      if (wrthis(indxV)) call ncwrite(nchis,vname(1,indxV),v(i0:i1, 1:j1,:,nstp),start)

      call wrt_his_trc(start)                              ! tracer variables

      if (wrthis(indxR)) then
# ifdef SPLIT_EOS
        call ncwrite(nchis,vname(1,indxR),rho1(i0:i1,j0:j1,:),start)
# else
        call ncwrite(nchis,vname(1,indxR), rho(i0:i1,j0:j1,:),start)
# endif
      endif
      if (wrthis(indxO)) then                              ! s-coordinate omega vertical velocity (m/s).
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
            enddo
          enddo
        enddo
        call ncwrite(nchis,vname(1,indxO),work(i0:i1,j0:j1,:),start)
      endif
      if (wrthis(indxW)) then                              ! true vertical velocity (m/s).
        do tile=0,NSUB_X*NSUB_E-1
          call wvlcty (tile, work)
        enddo
        ! need 1:N as work is 1:N+1 in size
        ! put as 0:N-1 because z-dimension array is declared in correctly in wvlcty_tile.
        call ncwrite(nchis,vname(1,indxW),work(i0:i1,j0:j1,0:N-1),start)
      endif
      if (wrthis(indxAkv)) call ncwrite(nchis,vname(1,indxAkv),Akv(i0:i1,j0:j1,:),start)
      if (wrthis(indxAkt)) call ncwrite(nchis,vname(1,indxAkt),Akt(i0:i1,j0:j1,:,itemp),start)
# ifdef SALINITY
      if (wrthis(indxAks)) call ncwrite(nchis,vname(1,indxAks),Akt(i0:i1,j0:j1,:,isalt),start)
# endif
      start(3)=record                                      ! back to 2D vars
# ifdef LMD_KPP
      if (wrthis(indxHbls)) call ncwrite(nchis,vname(1,indxHbls),hbls(i0:i1,j0:j1,nstp),start)
# endif
# ifdef LMD_BKPP
      if (wrthis(indxHbbl)) call ncwrite(nchis,vname(1,indxHbbl),hbbl(i0:i1,j0:j1,nstp),start)
# endif
#endif
      start(3)=1; start(4)=record                          ! back to 3D vars
#ifdef WEC
      call wrt_his_wec(record, start)
#endif

#if defined BIOLOGY_BEC2 && defined BEC2_DIAG
      call wrt_his_bgc_diag(record, ierr)
#endif

  1   format(/1x, '### ERROR: wrt_his :: Cannot write variable ''',
     &              A, ''' into history file, rec =', i6, 3x,A,i4)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to the disk to allow other
! processes to access data immediately after it is written.

#if defined MPI && !defined PARALLEL_FILES
      ierr=nf90_close (nchis)
      if (nrpfhis > 0 .and. record >= nrpfhis) nchis=-1
#else
      if (nrpfhis > 0 .and. record >= nrpfhis) then
        ierr=nf90_close (nchis)
        nchis=-1
      else
        ierr=nf90_sync(nchis)
      endif
#endif
      if (ierr == nf90_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(tstring,'(F18.8)') tdays
          i=1
          do while (i < 18 .and. tstring(i:i) == ' ')
            i=i+1
          enddo
          write(*,'(7x,A,1x,A,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &      'wrt_his :: wrote history, tdays =', tstring(i:i+8),
     &      'step =', ibuff(1),  'rec =', record, '/',  nrechis
     &       MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(/1x,2A/)')      '### ERROR: wrt_his :: Cannot ',
     &                       'synchronize or close history file.'
        if (may_day_flag == 0) may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                           mynode+1, 2, ocean_grid_comm, ierr)
#endif
      end subroutine wrt_his  !]
