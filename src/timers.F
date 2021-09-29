#include "cppdefs.opt"
#ifndef NO_COMPILER_SUPPORT_FOR_TIMING
! Most modern compilers recognize "cpu_time" and OpenMP "omp_get_wtime"
! functions as intrinsics, so the above switch should never be defined.
! Otherwise, stub functions for start/stop_timers are provided for
! compatibility at the end of this file.

! Machine-dependent issues: on some (older, now obsolete) platforms
! Open MP threads are KERNEL-LEVEL threads, which means that they have
! distinct process IDs and their CPU consumption can be measured
! separately. Summation of CPU times is performed in this case; on
! others they are USER-LEVEL threads (as mandated by POSIX standard),
! hence it is no longer possible to distinguish PIDs and CPU times
! consumed by individual threads. Report the maximum instead of sum.
! Additionally, Open MP (standard v. 2.0) function "omp_get_wtime" may
! not be implemented on some platforms, so its use is conditionally
! avoided.

# undef KERNEL_THREADS

# if defined SGI || defined PGI
#  undef OMP_GET_WTIME
# else
#  define OMP_GET_WTIME
# endif

      subroutine start_timers

      use param
      use comm_vars
      use scalars

      implicit none
      integer getpid, trd

# ifdef OMP_GET_WTIME
C$    real*8 omp_get_wtime
# endif
C$    integer omp_get_thread_num, omp_get_num_threads

      numthreads=1     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      proc(1)=getpid() ; proc(2)=trd

# ifdef OMP_GET_WTIME
C$    WallClock=omp_get_wtime()
# endif
      call cpu_time(cpu_init)

C$OMP CRITICAL (start_timers_cr_rgn)
      if (trd_count == 0) then
# ifdef MPI
        if (mynode == 0) then
          write(*,'(/1x,2(A,I4,A,I2,A,I3),2(A,I4),A,I3)')
     &      'NUMBER OF NODES:', NNODES, '(', NP_XI, ' x', NP_ETA,
     &      ') THREADS:',  numthreads,  ' TILING:',
     &       NSUB_X,' x', NSUB_E, ' GRID:',  LLm,' x',MMm,' x',N
        endif
# else
        write(*,'(/3(1x,A,I3),4x,2(A,I4),A,I3)')
     &       'NUMBER OF THREADS:',     numthreads,     'TILING:',
     &        NSUB_X,'x',NSUB_E, 'GRID SIZE:', Lm,' x',Mm,' x',N
# endif
      endif
      trd_count=trd_count+1
# ifdef MPI
#  ifndef MPI_SILENT_MODE
      write(*,'(4x,A,I4,1x,A,I3,1x,A,I10,A)') 'Process', mynode,
     &        'thread', proc(2), '(pid=', proc(1), ') is active.'
#  endif
# else
      write(*,'(8x,A,I3,1x,A,i10,A)') 'Thread #', proc(2),
     &                           '(pid=', proc(1), ') is active.'
# endif
      if (trd_count == numthreads) then
        trd_count=0
C$      mpi_master_only write(*,'(1x,2A/)') 'This code was ',
C$   &               'built using Open MP enabled compiler.'
      endif
C$OMP END CRITICAL (start_timers_cr_rgn)
      end

      subroutine stop_timers()            ! Finalize timing

      use param
      use comm_vars
      use scalars

      implicit none                       ! for all threads.

# ifdef OMP_GET_WTIME
C$    real*8 omp_get_wtime
# endif
      if (proc(1) /= 0) then
        proc(1)=0
# ifdef OMP_GET_WTIME
C$      WallClock=omp_get_wtime()-WallClock
# endif
        call cpu_time(cpu_net) ; cpu_net=cpu_net-cpu_init

C$OMP CRITICAL (stop_timers_cr_rgn)
# ifdef MPI
#  ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#  endif
          write(*,'(1x,A,I5,2x,A,I3,2x,A,F12.2,1x,A)') 'Process',
     &    mynode, 'thread', proc(2), 'cpu time =', cpu_net, 'sec'
#  ifdef MPI_SILENT_MODE
        endif
#  endif
# else
        if (trd_count == 0) write(*,*)
        write(*,'(13x,A,I3,2x,A,F12.2,1x,A)') 'thread #', proc(2),
     &                                'cpu time =', cpu_net, 'sec'
# endif
# ifdef KERNEL_THREADS
        cpu_all(1)=cpu_all(1)      +cpu_net
# else
        cpu_all(1)=max(cpu_all(1), cpu_net)
# endif
        trd_count=trd_count+1
        if (trd_count == numthreads) then
          trd_count=0
# ifdef MPI_SILENT_MODE
          if (mynode == 0) then
# endif
# ifdef KERNEL_THREADS
C$          write(*,'(29x,A,F14.2)')   'total', cpu_all(1)
# else
C$          write(*,'(27x,A,F14.2)') 'maximum', cpu_all(1)
# endif
# ifdef OMP_GET_WTIME
C$          write(*,'(11x,A,F12.2,1x,A,F6.2,A)')
C$   &           'Wall Clock elapsed time =', WallClock, 'sec (',
C$   &  100.D0*cpu_all(1)/(WallClock*dble(numthreads)), '% CPUs)'
# endif
# ifdef MPI_SILENT_MODE
          endif
# endif
        endif
C$OMP END CRITICAL (stop_timers_cr_rgn)
      endif
      end

! The following routine is to catch loss of synchronization in Open
! MP mode. Calls to "sync_trap" are not hardcoded into the model, but
! are inserted by "mpc" if directed to do so.  The algorithm for
! trapping works as follows: every thread advances its own private
! counter "priv_count", and then global counter "barr_count" (inside
! critical region) and compares the value of global counter with its
! private. Since each thread increments the global counter by 1, it
! grows numthreads-times faster than the private, hence within each
! synchronization region the global counter, after incremented by 1
! by a thread must have values from "previous"+1 to last
! "previous"+numthreads inclussive (here "previous" means the final
! value after the previous synchronization event.  As the result,
! "itest" computed below must always match "priv_count".

      subroutine sync_trap(ibarr)
      use param
      use scalars
      implicit none
      integer ibarr, indx, itest

      indx=1+mod(ibarr-1,16)
      priv_count(indx)=priv_count(indx)+1
C$OMP CRITICAL(trap_cr_rgn)
      barr_count(indx)=barr_count(indx)+1
      itest=1+(barr_count(indx)-1)/numthreads
      if (itest /= priv_count(indx)) then
        write(*,'(A,3I10)') 'sync error', ibarr,
     &        priv_count(indx),  barr_count(indx)
      elseif (mod(priv_count(indx),4001) == 0) then
        write(*,'(A,I12,2(2x,A,I3))') 'barrier count =',
     &              priv_count(indx), 'barr# =', ibarr,
     &                                'trd =',  proc(2)
      endif
C$OMP END CRITICAL(trap_cr_rgn)
      end

#else
                                       ! These are stub-routines for
      subroutine start_timers          ! compatibility with compilers
      use param
      use scalars
      implicit none                    ! without OpenMP support.

      mpi_master_only write(*,'(/2(1x,A,I3),4x,2(1x,A,I4)/)')
     &                        'BLOCKING:', NSUB_X, 'x', NSUB_E,
     &                        'HORIZ. GRID SIZE:', Lm, 'x', Mm
      end

      subroutine stop_timers()
      end
#endif
