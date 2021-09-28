#include "cppdefs.h"

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
      use ncvars
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
      use ncvars
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
      use ncvars
      use scalars
      use mpi

      implicit none

#include "forces.h"
      integer ierr,  record, i, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(iaux), ncdf_write
      character(len=18) tstring

#include "netcdf.inc"

#if defined MPI && !defined PARALLEL_FILES
      integer status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif

! Create/open history file; write grid arrays, if so needed.

      call def_frc(nrechis, ierr)
      if (ierr /= nf_noerr) goto 99
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

      ierr=nf_put_vara_int (nchis, hisTstep, start, count, ibuff)
      if (ierr /= nf_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_frc :: ',
     &       'Cannot write variable ''time_step'' into history file, ',
     &       'rec =', record, nf_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif

! Time

      ierr=nf_put_var1_FTYPE (nchis, hisTime, record, time)
      if (ierr /= nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,'(/1x,4A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_frc :: ',
     &        'Cannot write variable ''',    vname(1,indxTime)(1:lvar),
     &        ''' into history file, rec =', record, nf_strerror(ierr)
     &            MYID
        goto 99                                           !--> ERROR
      endif

! Short wave radiation flux

      ierr=ncdf_write(nchis, vid_srfx,  record,  r_var, srflx, 1)
      if (ierr /= nf_noerr) then
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




      subroutine def_frc(total_rec, ierr)

      use param
      use ncvars
      use scalars

      implicit none
      logical create_new_file
      integer total_rec, ierr, rec,timedim, prev_fill_mode, lfnm, lvar
     &      , r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)

      integer my_nf_def_dim     !<-- instrumented netCDF function.
      integer lenstr, checkdims
      character(len=64) text

#include "netcdf.inc"

#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real*8, parameter :: spv_set=1.D+33
# else
      real*4, parameter :: spv_set=1.E+33
# endif
#endif

      hisname='swrflux.nc' ; lfnm=lenstr(hisname)
      create_new_file=ldefhis
      if (nchis /= -1)  create_new_file=.false.
# ifdef MPI
#  ifdef PARALLEL_FILES
      ierr=0
      call insert_node(hisname, lfnm, mynode, NNODES, ierr)
#  else
      if (mynode > 0) create_new_file=.false.
#  endif
# endif
      if (nrpfhis > 0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfhis))
        call insert_time_index(hisname, lfnm, lvar, ierr)
      endif
      if (ierr /= 0) goto 99

  10  if (create_new_file) then
        ierr=nf_create(hisname(1:lfnm), nf_clobber+nf_64bit_offset,
     &                                                     nchis)
        if (ierr == nf_noerr) then
          ierr=nf_set_fill(nchis, nf_nofill, prev_fill_mode)
          if (ierr /= nf_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_frc :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_frc :: ',
     &          'Cannot create netCDF file ''', hisname(1:lfnm),
     &                              '''.',   nf_strerror(ierr)
          goto 99
        endif
        if (nrpfhis == 0) total_rec=0

        call put_global_atts(nchis, ierr)

! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------

        ierr=my_nf_def_dim(nchis, 'xi_rho',  xi_rho,  r2dgrd(1))
        ierr=my_nf_def_dim(nchis, 'xi_u',    xi_u,    u2dgrd(1))
        ierr=my_nf_def_dim(nchis, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=my_nf_def_dim(nchis, 'eta_v',   eta_v,   v2dgrd(2))

        ierr=my_nf_def_dim(nchis, 'time', nf_unlimited, timedim)
        ierr=my_nf_def_dim(nchis, 'auxil',   iaux,     auxil(1))

        auxil(2)=timedim

        r2dgrd(3)=timedim           ! 2D ZETA-type

        u2dgrd(2)=r2dgrd(2)         ! 2D UBAR-type
        u2dgrd(3)=timedim

        v2dgrd(1)=r2dgrd(1)         ! 2D VBAR-type
        v2dgrd(3)=timedim

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);

        ierr=nf_def_var (nchis, 'time_step', nf_int, 2, auxil, hisTstep)
        ierr=nf_put_att_text(nchis, hisTstep, 'long_name', 48,
     &              'time step and record numbers from initialization')

! Time.
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (nchis, vname(1,indxTime)(1:lvar), NF_FTYPE,
     &                                         1, timedim, hisTime)
        text=vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (nchis, hisTime, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (nchis, hisTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))

! Short wave radiation flux

        ierr=nf_def_var(nchis, 'srflx',  NF_FOUT, 3, r2dgrd, vid_srfx)
        text='short wave radiation flux' ; lvar=lenstr(text)
        ierr=nf_put_att_text (nchis, vid_srfx, 'long_name', lvar,
     &                                            text(1:lvar))
        text='kinematic units' ;  lvar=lenstr(text)
        ierr=nf_put_att_text (nchis, vid_srfx, 'units',     lvar,
     &                                            text(1:lvar))
#ifdef MASK_LAND_DATA
        ierr=nf_put_att_FOUT (nchis, vid_srfx, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
        if (ierr /= nf_noerr) then
           write(*,*)  'nf_put_att_XXX:', nf_strerror(ierr)
        endif
#endif

! Leave definition mode.

        ierr=nf_enddef(nchis)
        mpi_master_only write(*,'(7x,5A,I4)') 'def_frc :: ',
     &       'created new netCDF file ''', hisname(1:lfnm),
     &       '''.' MYID

! Open existing file and prepare for appending data:
! ==== ======== ==== === ======= === ========= ====

      elseif (nchis == -1) then  !--> alt. to "create_new_file"
        ierr=nf_open(hisname(1:lfnm), nf_write, nchis)
        if (ierr == nf_noerr) then
          ierr=nf_set_fill(nchis, nf_nofill, prev_fill_mode)
          if (ierr == nf_noerr) then
            mpi_master_only write(*,'(7x,5A,I4)') 'def_frc :: ',
     &                 'open existing file ''', hisname(1:lfnm),
     &              '''.' MYID
          else
            write(*,'(1x,4A,I4)') '### ERROR: def_frc :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &         MYID
            goto 99
          endif

          ierr=checkdims(nchis, hisname, rec)
          if (ierr == nf_noerr) then
            if (nrpfhis == 0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfhis))
            endif
            if (ierr > 0) then
              mpi_master_only write(*,
     &                '(/1x,A,I5,1x,3A/21x,2(A,I5),1x,A/21x,A/)')
     &          'WARNING: def_frc :: The actual number of records',
     &             rec, 'present in file ''',  hisname(1:lfnm),  '''',
     &            'exceeds record', rec+1-ierr, '/', total_rec,
     &            'specified by restart initial conditions.',
     &        'All records beyond this number will be overwritten.'
              rec=rec-ierr
            elseif (nrpfhis == 0) then
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
            create_new_file=.true.
            goto 10
          else
            write(*,'(/1x,5A,I4/)')    '### ERROR: def_frc :: ',
     &         'Cannot open ''', hisname(1:lfnm), '''.' MYID
            goto 99
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif

! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
        ierr=nf_inq_varid (nchis, 'time_step', hisTstep)
        if (ierr /= nf_noerr) then
          write(*,1) 'time_step', hisname(1:lfnm)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (nchis,vname(1,indxTime)(1:lvar),hisTime)
        if (ierr /= nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), hisname(1:lfnm)
          goto 99                                         !--> ERROR
        endif

        ierr=nf_inq_varid(nchis, 'srflx', vid_srfx)
        if (ierr /= nf_noerr) then
          write(*,1) vname(1,indxZ)(1:lvar), hisname(1:lfnm)
          goto 99                                       !--> ERROR
        endif

        mpi_master_only write(*,'(7x,4A,I4,2x,A,I4)') 'def_frc :: ',
     &                    'open existing file ''', hisname(1:lfnm),
     &                    ''' from rec =', rec MYID

#if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf_open(hisname(1:lfnm), nf_write, nchis)
        if (ierr == nf_noerr) then
# ifndef MPI_SILENT_MODE
          write(*,'(7x,5A,I4)')   'def_frc :: reopen existing ',
     &            'file ''', hisname(1:lfnm), '''.' MYID
# endif
        else
          write(*,'(/1x,6A,I4/)')  '### ERROR: def_frc :: Cannot ',
     &                       'open file ''', hisname(1:lfnm), '''.',
     &                        nf_strerror(ierr) MYID
          goto 99                                   !--> ERROR
        endif
#endif
      endif             !<--  create_new_file
  99  return

  1   format(/1x,'### ERROR: def_frc :: Cannot find variable ''',
     &                           A, ''' in file ''', A, '''.'/)
      end
