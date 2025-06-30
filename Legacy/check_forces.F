#include "cppdefs.opt"

! Check diurnal cycle (short-wave radiation flux) as it is read from
! forcing files, interpolated in time, and modulated by an analytical
! formula.  This program is nothing, but a stripped down version of
! ROMS main.F, which goes through initialization exactly as the actual
! model, but then, instead of doing actual time stepping, it only calls
! "get_forces" and "set_forces" routines, and saves the outcome into
! netCDF file at the same interval as history file. The content of this
! file is as follows:
!
!    program check_forces    <-- exactly the same as roms main;
!    subroutine roms_thread  <-- exactly the same as roms_thread;
!    subroutine roms_init    <-- simplified roms init (only sets
!                                   up grid and init conditions);
!    subroutine roms_step    <-- reads forces;
!    subroutine wrt_forces   <-- substitute for wrt_his
!    subroutine def_frc      <-- substitute for def_his
!
! Compilation: compile this file into .o and use it as a substitute for
! main.o of the actual ROMS code (all other files must be the same).
!
! Running: exactly the same as running the ROMS code: it reads the same
! "roms.in", grid, initial, and forcing files.



!!    program check_forces       ! Open MP version of ROMS driver
      implicit none              ! with single parallel region using
      use param
      use mpi
      integer ierr               ! explicit barrier synchronization.
#ifdef MPI
      real*8 tstart, tend
C$    integer level,req_lev

c**   call system('uname -nmr')
      ierr=1
C$    req_lev=MPI_THREAD_MULTIPLE
C$    call MPI_Init_thread(req_lev, level, ierr)
C$ !!   write(*,*) 'MPI thread support levels =', req_lev,level
C$    ierr=0
      if (ierr == 1) call MPI_Init(ierr)

      call mpi_setup(ierr)
      tstart=MPI_Wtime()
      if (ierr == 0) then
#endif
        call init_scalars(ierr)          ! Initialize global scalars,
        if (ierr == 0) then              ! model tunable parameters,
C$        call omp_set_dynamic(.false.)
C$OMP PARALLEL                           ! fast-time averaging weights
          call roms_thread               ! for barotropic mode, and
C$OMP END PARALLEL                       ! launch the model in OpenMP
        endif                            ! parallel regime.
#ifdef MPI
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      tend=MPI_Wtime()
      mpi_master_only write(*,*) 'MPI_run_time =', tend-tstart
      call MPI_Finalize (ierr)
#endif
      stop
      end


      subroutine roms_thread
      use param
      use scalars
      implicit none
                                         ! Note: Because there is
      call start_timers()                ! a possibility of I/O error
      call roms_init                     ! occurring on some MPI nodes,
      if (may_day_flag /= 0) goto 99     ! but not simultaneously on
      do iic=ntstart,ntstart+ntimes      ! all, exiting is deferred
        diag_sync=.false.                ! until "may_day_flag" is
        call roms_step
        if (diag_sync .and. may_day_flag /= 0) goto 99
      enddo
  99  call stop_timers()                 ! summarized among all nodes
C$OMP BARRIER                            ! and broadcasted by "diag"
C$OMP MASTER                             ! [which sets diag_sync=T to
      call closecdf                      ! signal that "may_day_flag"
C$OMP END MASTER                         ! is in uniform state]. This
      return                             ! is THE ONLY route to exit
      end                                ! if something goes wrong.


      subroutine roms_init

      use param
      use scalars

      implicit none
      integer trd, tile, my_first, my_last, range
C$    integer omp_get_thread_num, omp_get_num_threads

      numthreads=1 ; trd=0
C$    numthreads=omp_get_num_threads()
C$    trd=omp_get_thread_num()
      proc(2)=trd

      if (mod(NSUB_X*NSUB_E,numthreads) /= 0) then
C$OMP MASTER
        mpi_master_only write(*,'(/3(1x,A,I3),A/)')
     &    '### ERROR: Wrong choice of numthreads =', numthreads,
     &         'while NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E,'.'
        may_day_flag=8
C$OMP END MASTER
C$OMP BARRIER
        goto 99 !-->  EXIT
      endif
                                       ! NOTE: This code is written
      iic=0  ; kstp=1 ; knew=1         ! under an assumption that all
#ifdef SOLVE3D
      iif=1  ; nstp=1                  ! the scalar variables assigned
      nnew=1 ; nrhs=1                  ! on the left here are placed
#endif
      synchro_flag=.true.              ! into a THREADPRIVATE common
      diag_sync=.false.                ! block so each thread must
      priv_count=0                     ! assign them to same values.

      range=(NSUB_X*NSUB_E+numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
#define my_tile_range my_first,my_last

      do tile=my_tile_range            ! Initialize or "first touch"
        call init_arrays(tile)         ! global arrays (most of them
      enddo                            ! are just set to to zero).
C$OMP BARRIER
C$OMP MASTER
      call get_grid
C$OMP END MASTER
C$OMP BARRIER
      if (may_day_flag /= 0) goto 99 !-->  EXIT

      do tile=my_tile_range
        call setup_grid1(tile)
      enddo
C$OMP BARRIER
      do tile=my_tile_range
        call setup_grid2(tile)
      enddo
C$OMP BARRIER

#ifdef SOLVE3D
C$OMP MASTER                           ! Setup vertical stretching
      call set_scoord                  ! functions for S-coordinate
C$OMP END MASTER                       ! system: needed here only
C$OMP BARRIER                          ! because get_init below
      if (may_day_flag /= 0) goto 99   ! checks for matching.
#endif

      if (nrrec > 0) then
C$OMP MASTER
        call get_init(nrrec, 1)
C$OMP END MASTER
      else
C$OMP MASTER
        write(*,*) '### ERROR: Record in initial file not specified.'
        may_day_flag=8
C$OMP END MASTER
C$OMP BARRIER
        if (may_day_flag /= 0) goto 99 !-->  EXIT
      endif    !<-- nrrec > 0

      time=start_time
      tdays=time*sec2day
  99  return
      end


!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     ***
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***


      subroutine roms_step

      use param

      use scalars

      implicit none
      integer trd, tile, my_first, my_last, range

      trd=proc(2)
      range=(NSUB_X*NSUB_E+numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      time=start_time+dt*float(iic-ntstart) !<-- corresp. to "nstp"
      tdays=time*sec2day
#ifdef SOLVE3D
      nstp=1+mod(iic-ntstart,2)
      nrhs=nstp
      nnew=3
#endif
C$OMP MASTER
      mpi_master_only write(*,'(1x,I7,F12.5)') iic-1, tdays
      call get_forces
C$OMP END MASTER
C$OMP BARRIER

#ifdef SOLVE3D
      do tile=my_tile_range                 ! interpolate forcing
        call set_forces(tile)               ! data to model time and
      enddo
C$OMP BARRIER

! Output block: write files.
!======= ====== ===== =============== ======
      if (mod(iic-ntstart,nwrt) == 0) then
C$OMP MASTER
        nrechis=nrechis+1
        call wrt_forces
C$OMP END MASTER
C$OMP BARRIER
        if (iic-ntstart > ntimes) goto 99   !-->  DONE
      endif
  99  return
      end



# define vid_srfx hisZ

      subroutine wrt_forces

      use param

      use scalars
      use mpi
      use netcdf

      implicit none

#include "forces.h"
      integer ierr,  record, i, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(iaux), ncdf_write
      character(len=18) tstring

#if defined MPI && !defined PARALLEL_FILES
      integer status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif

! Create/open history file; write grid arrays, if so needed.

      call def_frc(nrechis, ierr)
      if (ierr /= nf90_noerr) goto 99
      lstr=lenstr(hisname)

                                        !!! WARNING: Once time
! Set record within the file.           !!! stepping has been
                                        !!! started, it is assumed
      nrechis=max(nrechis,1)            !!! that the global history
      if (nrpfhis == 0) then            !!! record index "nrechis"
        record=nrechis                  !!! is advanced by main.
      else
        record=1+mod(nrechis-1, nrpfhis)
      endif

      ibuff(1)=iic-1  ; ibuff(2)=nrecrst      ! Time step and
      ibuff(4:iaux)=0 ; ibuff(3)=nrechis      ! record numbers.
      start(1)=1      ; count(1)=iaux
      start(2)=record ; count(2)=1

      ierr=nf90_put_var(nchis, hisTstep, ibuff, start, count)
      if (ierr /= nf90_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_frc :: ',
     &       'Cannot write variable ''time_step'' into history file, ',
     &       'rec =', record, nf90_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif

! Time

      ierr=nf90_put_var (nchis, hisTime, time, (/record/) )
      if (ierr /= nf90_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,'(/1x,4A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_frc :: ',
     &        'Cannot write variable ''',    vname(1,indxTime)(1:lvar),
     &        ''' into history file, rec =', record, nf90_strerror(ierr)
     &            MYID
        goto 99                                           !--> ERROR
      endif

! Short wave radiation flux

      ierr=ncdf_write(nchis, vid_srfx,  record,  r_var, srflx, 1)
      if (ierr /= nf90_noerr) then
        write(*,1) 'srflx', record MYID
        goto 99                                         !--> ERROR
      endif


  1   format(/1x, '### ERROR: wrt_frc :: Cannot write variable ''',
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
     &      'wrt_frc :: wrote history, tdays =', tstring(i:i+8),
     &      'step =', ibuff(1),  'rec =', record, '/',  nrechis
     &       MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(/1x,2A/)')      '### ERROR: wrt_frc :: Cannot ',
     &                       'synchronize or close history file.'
        if (may_day_flag == 0) may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                           mynode+1, 2, ocean_grid_comm, ierr)
#endif
      end


