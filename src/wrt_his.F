#include "cppdefs.h"
                               ! Writes requested model
      subroutine wrt_his       ! fields at requested levels
                               ! into history netCDF file.

      ! NOTE: wrt_his is called in middle of roms_step in main.F because
      !       even though u/ubar, etc, were computed at the end of the
      !       previous timestep, other variables like 'rho', 'Akv', etc
      !       have only just been computed.

      use param
#ifdef WEC
      use wec_frc, only: wrt_his_wec
#endif

      use tracers, only: wrt_his_trc
      use bgc ! wrt_his_bgc_diag
      use read_write ! nc_write_var
      use netcdf ! for nf90_routines
      use eos_vars
      use grid
      use mixing
      use ncvars

      implicit none

#include "scalars.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "netcdf.inc"

      integer ierr,  record, i, lstr, lvar, lenstr, var_id_tmp
     &   , start(2), count(2), ibuff(iaux), ncdf_write
      character(len=18) tstring
#ifdef SOLVE3D
      integer tile, itrc, j,k
# include "work.h"
#endif
#if defined MPI && !defined PARALLEL_FILES
# include "mpif.h"
      integer status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif


! Create/open history and set record index within the file.

      call def_his(nrechis, ierr)
      if (ierr /= nf_noerr) goto 99
      lstr=lenstr(hisname)
      nrechis=max(nrechis,1)              !!! WARNING: Once time
      if (nrpfhis == 0) then              !!! stepping has been
        record=nrechis                    !!! started, it is assumed
      else                                !!! that the global history
        record=1+mod(nrechis-1, nrpfhis)  !!! record index "nrechis"
      endif                               !!! is advanced by main.


c-#define CR
CR      write(*,*) 'wrt_his: Entry  ' MYID

! Write out evolving model variables:
! ----- --- -------- ----- ----------

      ibuff(1)=iic-1  ; ibuff(2)=nrecrst    ! Time step and
      ibuff(4:iaux)=0 ; ibuff(3)=nrechis    ! record numbers.
#ifdef AVERAGES
      ibuff(4)=nrecavg
#endif
      start(1)=1      ; count(1)=iaux
      start(2)=record ; count(2)=1

      ! get time step var_id & record time step info
      ! --------------------------------------------
      ierr=nf90_inq_varid(nchis, 'time_step', var_id_tmp)
      ierr=nf90_put_var(nchis, var_id_tmp, ibuff, start, count)
      if (ierr /= nf_noerr) then
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

      ! write ocean time
      ! ----------------
      ierr=nf90_inq_varid(nchis, 'ocean_time', var_id_tmp)
      ierr=nf90_put_var(nchis, var_id_tmp, time,(/record/))

! Barotropic mode variables: free-surface and 2D momentum
! components in XI-,ETA-directions.
! 'knew' index is normally 'n+1' timestep, but at this point
! its value hasn't evolved from 'n+1' of previous timestep,
! which makes it 'n' of this timestep still.

      if (wrthis(indxZ)) then
        call nc_write_var(nchis, zeta(START_2D_ARRAY,knew), 1,
     &                    vname(1,indxZ), rp_var, record, ierr)
      endif
      if (wrthis(indxUb)) then
        call nc_write_var(nchis, ubar(START_2D_ARRAY,knew), 1,
     &                    vname(1,indxUb), up_var, record, ierr)
      endif
      if (wrthis(indxVb)) then
        call nc_write_var(nchis, vbar(START_2D_ARRAY,knew), 1,
     &                    vname(1,indxVb), vp_var, record, ierr)
      endif

#ifdef SOLVE3D

! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! Currently not sure why Sasha chose to do wrt_his in middle of timestep,
! rather than at the beginning.

      if (wrthis(indxU)) then
        call nc_write_var(nchis, u(START_2D_ARRAY,1,nstp), N,
     &                    vname(1,indxU), up_var, record, ierr)
      endif
      if (wrthis(indxV)) then
        call nc_write_var(nchis, v(START_2D_ARRAY,1,nstp), N,
     &                    vname(1,indxV), vp_var, record, ierr)
      endif

! Tracer variables.

      call wrt_his_trc(record, ierr)

! Density anomaly.

      if (wrthis(indxR)) then
# ifdef SPLIT_EOS
        call nc_write_var(nchis, rho1, N, vname(1,indxR), rp_var, record, ierr)
# else
        call nc_write_var(nchis, rho, N, vname(1,indxR), rp_var, record, ierr)
# endif
      endif

! S-coordinate omega vertical velocity (m/s).

      if (wrthis(indxO)) then
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
            enddo
          enddo
        enddo
        call nc_write_var(nchis, work, N+1, vname(1,indxO), rp_var, record, ierr)
      endif

! True vertical velocity (m/s).

      if (wrthis(indxW)) then
        do tile=0,NSUB_X*NSUB_E-1
          call wvlcty (tile, work)
        enddo
        call nc_write_var(nchis, work, N, vname(1,indxW), rp_var, record, ierr)
      endif

! Vertical viscosity coefficient.

      if (wrthis(indxAkv)) then
        call nc_write_var(nchis, Akv, N+1, vname(1,indxAkv), rp_var, record, ierr)
      endif

! Vertical diffusion coefficient for potential temperature.

      if (wrthis(indxAkt)) then
        call nc_write_var(nchis, Akt(START_2D_ARRAY,0,itemp), N+1,
     &                    vname(1,indxAkt), rp_var, record, ierr)
      endif
# ifdef SALINITY

! Vertical diffusion coefficient for salinity.

      if (wrthis(indxAks)) then
        call nc_write_var(nchis, Akt(START_2D_ARRAY,0,isalt), N+1,
     &                    vname(1,indxAks), rp_var, record, ierr)
      endif
# endif
# ifdef LMD_KPP
! Depth of planetary boundary layer [m].

      if (wrthis(indxHbls)) then
        call nc_write_var(nchis, hbls(START_2D_ARRAY,nstp), 1,
     &                    vname(1,indxHbls), rp_var, record, ierr)
      endif
# endif
# ifdef LMD_BKPP
! Thickness of KPP bottom boundary layer [m].

      if (wrthis(indxHbbl)) then
        call nc_write_var(nchis, hbbl(START_2D_ARRAY,nstp), 1,
     &                    vname(1,indxHbbl), rp_var, record, ierr)
      endif
# endif
#endif

      ! Wave Effect on Current (WEC) variables
#ifdef WEC
      call wrt_his_wec(record, ierr, may_day_flag)
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
      ierr=nf_close (nchis)
      if (nrpfhis > 0 .and. record >= nrpfhis) nchis=-1
#else
      if (nrpfhis > 0 .and. record >= nrpfhis) then
        ierr=nf_close (nchis)
        nchis=-1
      else
        ierr=nf_sync(nchis)
      endif
#endif
      if (ierr == nf_noerr) then
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
      end
