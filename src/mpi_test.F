#include "cppdefs.opt"

!!    program main               ! Open MP version of ROMS driver
      use param
      use mpi
      implicit none              ! with single parallel region using
      integer ierr               ! explicit barrier synchronization.
#ifdef MPI
      real*8 tstart, tend
C$    integer level,req_lev

      call system('uname -nmr')
      ierr=1
C$    req_lev=MPI_THREAD_MULTIPLE
C$    call MPI_Init_thread (req_lev, level, ierr)
C$    write(*,*) 'MPI thread support levels=', req_lev,level
C$    ierr=0
      if (ierr==1) call MPI_Init (ierr)

      call mpi_setup (ierr)
      tstart=MPI_Wtime()
      if (ierr==0) then
#endif
        if (ierr==0) then
C$        call omp_set_dynamic(.false.)
C$OMP PARALLEL
          call mpi_test_thread
C$OMP END PARALLEL
        endif
#ifdef MPI
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      tend=MPI_Wtime()
      mpi_master_only write(*,*) 'run_time =', tend-tstart
      call MPI_Finalize (ierr)
#endif
      stop
      end



      subroutine mpi_test_thread
      use param
      use scalars
      implicit none
      integer trd, tile, my_first, my_last, range
C$    integer omp_get_thread_num, omp_get_num_threads

      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      proc(2)=trd
      if (mod(NSUB_X*NSUB_E,numthreads)/=0) then
C$OMP MASTER
        mpi_master_only write(*,'(/3(1x,A,I3),A/)')
     &    '### ERROR: wrong choice of numthreads =', numthreads,
     &         'while NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E,'.'
        may_day_flag=8
C$OMP END MASTER
C$OMP BARRIER
        goto 99 !-->  EXIT
      endif
      range=(NSUB_X*NSUB_E+numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
#define my_tile_range my_first,my_last

      do tile=my_tile_range
        call mpi_test_init (tile)
      enddo
C$OMP BARRIER
C$OMP MASTER
      mpi_master_only write(*,'(/1x,A/)') 'Printing differences...'
C$OMP END MASTER
C$OMP BARRIER
      do tile=my_tile_range
        call mpi_test_check (tile)
      enddo
C$OMP BARRIER
  99  return
      end

      subroutine mpi_test_init (tile)
      use param
      use hidden_mpi_vars
      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call mpi_test_init_tile (istr,iend,jstr,jend)
      call mpi_test_copy_tile (istr,iend,jstr,jend)
      end

      subroutine mpi_test_init_tile (istr,iend,jstr,jend)

      use param
      use hidden_mpi_vars
      use mpi_test_arr
      implicit none
      integer istr,iend,jstr,jend, i,j,k, iabs,jabs
      real dx,dy, x,y

#include "compute_extended_bounds.h"

      call init_arrays_mpi_test_arr  ! DevinD added for allocatable array changes

#ifdef MPI
c      write(*,'(2(1x,A,2I5),1x,A,4I5)') 'i,j node =', inode,jnode,
c     &                         'i,j SW_corn =', iSW_corn,jSW_corn,
c     &               'istr,iend,jstr,jend =', istr,iend,jstr,jend
#endif
      dx=1./dble(LLm)         ! Note the "extended" bounds above:
      dy=1./dble(MMm)         ! arrays A2D,C2D,A3D,C3D assigned here
      do j=jstrR,jendR        ! are initialized everywhere, including
        do i=istrR,iendR      ! periodic and computational margins.
#ifdef MPI
          iabs=i+iSW_corn     ! absolute indices on the physical grid
          jabs=j+jSW_corn     ! corresponding to relative indices i,j
#else
          iabs=i
          jabs=j
#endif
#ifdef EW_PERIODIC
          if (iabs>LLm) then
            iabs=iabs-LLm
          elseif (iabs<1) then
            iabs=iabs+LLm
          endif
#endif
#ifdef NS_PERIODIC
          if (jabs>MMm) then
            jabs=jabs-MMm
          elseif (jabs<1) then
            jabs=jabs+MMm
          endif
#endif
          x=dx*dble(iabs)
          y=dy*dble(jabs)
          A2D(i,j)=x +x*x  -0.5*x*y +0.1*y
          C2D(i,j)=y+x -0.33333333333*x*y +0.00001*y*y*y*x
          A3D(i,j,0)= x+1.2*y-3.*x*y+0.1*x*x +0.2*y*y*y
          C3D(i,j,0)= -x+0.3*x*y +2*x*x*x -0.01*x*x*y*y
        enddo
      enddo

      do k=1,N
        x=dble(k)/dble(N)
        do j=jstrR,jendR
          do i=istrR,iendR
            x=dble(k)
            A3D(i,j,k)= A3D(i,j,0) + x+ 0.5*x*x  -0.025*A3D(i,j,0)*x
            C3D(i,j,k)= 0.01*C3D(i,j,0)*(x-0.0001*C3D(i,j,0))
     &                                                +0.00001*x*x*x
          enddo
        enddo
      enddo
      end

      subroutine mpi_test_copy_tile (istr,iend,jstr,jend)
      use param
      use mpi_test_arr
      implicit none
      integer istr,iend,jstr,jend, i,j,k

# include "compute_auxiliary_bounds.h"

      do j=jstrR,jendR        ! Note the "auxiliary" bounds above:
        do i=istrR,iendR      ! this copies the content of the array
          B2D(i,j)=A2D(i,j)   ! in the interior, but does not copy
          D2D(i,j)=C2D(i,j)   ! periodic and computational margins:
        enddo                 ! these must be filled by the exchange
      enddo                   ! calls below.
      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR
            B3D(i,j,k)=A3D(i,j,k)
            D3D(i,j,k)=C3D(i,j,k)
          enddo
        enddo
      enddo

c      call exchange2d_tile (istr,iend,jstr,jend, B2D)
c      call exchange_tile (istr,iend,jstr,jend, B3D,N+1)

c      call exchange2d_2_tile (istr,iend,jstr,jend, B2D, D2D)
c      call exchange_2_tile (istr,iend,jstr,jend, B3D,N+1, D3D,N+1)

c      call exchange2d_4_tile (istr,iend,jstr,jend, A2D,B2D,C2D,D2D)

c      call exchange_4_tile (istr,iend,jstr,jend, B3D,N+1, D3D,N+1,
c     &                                               B2D,1, D2D,1)


      call exchange2d_4_tile (istr,iend,jstr,jend, B3D,D3D,B2D,D2D)

      call exchange_2_tile (istr,iend,jstr,jend, B3D,N+1, D3D,N+1)

      end




      subroutine mpi_test_check (tile)
      use param
      use hidden_mpi_vars
      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call mpi_test_check_tile (istr,iend,jstr,jend)
      end

      subroutine mpi_test_check_tile (istr,iend,jstr,jend)
      use param
      use mpi_test_arr
      implicit none
      integer istr,iend,jstr,jend, i,j,k

# include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          if (B2D(i,j) /= A2D(i,j)) then
            write(*,'(4x,A,2I4,2F16.8)') '### ERR B2D: i,j =', i,j,
     &                                            A2D(i,j),B2D(i,j)
          endif
        enddo
      enddo


      do j=jstrR,jendR
        do i=istrR,iendR
          if (D2D(i,j) /= C2D(i,j)) then
            write(*,'(4x,A,2I4,2F16.8)') '### ERR D2D: i,j =', i,j,
     &                                            C2D(i,j),D2D(i,j)
          endif
        enddo
      enddo

      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR
            if (B3D(i,j,k) /= A3D(i,j,k)) then
              write(*,'(4x,A,2I4,I3,2F16.8)') '### ERR B3D: i,j,k =',
     &                                  i,j,k, A3D(i,j,k), B3D(i,j,k)
            endif
          enddo
        enddo
      enddo

      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR
            if (D3D(i,j,k) /= C3D(i,j,k)) then
              write(*,'(4x,A,2I4,I3,2F16.8)') '### ERR D3D: i,j,k =',
     &                                  i,j,k, C3D(i,j,k), D3D(i,j,k)
            endif
          enddo
        enddo
      enddo
      end

