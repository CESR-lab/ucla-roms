#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine grid_stiffness (tile)
      use param
      use hidden_mpi_vars
      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call grid_stiffness_tile (istr,iend,jstr,jend)
      return
      end

      subroutine grid_stiffness_tile (istr,iend,jstr,jend)

! Survey three-dimensional grid in order to determine maximum
! grid stiffness ratio:
!
!             z(i,j,k)-z(i-1,j,k)+z(i,j,k-1)-z(i-1,j,k-1)
!      r_x = ---------------------------------------------
!             z(i,j,k)+z(i-1,j,k)-z(i,j,k-1)-z(i-1,j,k-1)
!
! This is done for purely diagnostic purposes does not affect
! computations.

      use param
      use comm_vars
      use grid
      use ocean_vars
      use mpi

      implicit none
      integer istr,iend,jstr,jend, i,j,k, NSUB
      real my_rx0, my_rx1
# ifdef MPI
      real*8 buff(2)
      integer size, step, status(MPI_STATUS_SIZE), ierr
# endif
!
# include "compute_auxiliary_bounds.h"
!
      my_rx0=0. ; my_rx1=0.

      do j=jstr,jend
        do i=istrU,iend
# ifdef MASKING
          if (umask(i,j)>0.5) then
# endif
            my_rx0=max(my_rx0, abs( (z_w(i,j,0)-z_w(i-1,j,0))
     &                                   /(z_w(i,j,0)+z_w(i-1,j,0))
     &                                                            ))
            do k=1,N
              my_rx1=max(my_rx1, abs(
     &         (z_w(i,j,k)-z_w(i-1,j,k)+z_w(i,j,k-1)-z_w(i-1,j,k-1))
     &        /(z_w(i,j,k)+z_w(i-1,j,k)-z_w(i,j,k-1)-z_w(i-1,j,k-1))
     &                                                            ))
            enddo
# ifdef MASKING
          endif
# endif
        enddo
      enddo

      do j=jstrV,jend
        do i=istr,iend
# ifdef MASKING
          if (vmask(i,j)>0.5) then
# endif
            my_rx0=max(my_rx0, abs( (z_w(i,j,0)-z_w(i,j-1,0))
     &                                   /(z_w(i,j,0)+z_w(i,j-1,0))
     &                                                            ))
            do k=1,N
              my_rx1=max(my_rx1, abs(
     &         (z_w(i,j,k)-z_w(i,j-1,k)+z_w(i,j,k-1)-z_w(i,j-1,k-1))
     &        /(z_w(i,j,k)+z_w(i,j-1,k)-z_w(i,j,k-1)-z_w(i,j-1,k-1))
     &                                                            ))
            enddo
# ifdef MASKING
          endif
# endif
        enddo
      enddo

      if (SINGLE_TILE_MODE) then
        NSUB=1
      else
        NSUB=NSUB_X*NSUB_E
      endif

C$OMP CRITICAL (grd_stff_cr_rgn)
      if (tile_count==0) then
        rx0=my_rx0
        rx1=my_rx1
      else
        rx0=max(rx0, my_rx0)
        rx1=max(rx1, my_rx1)
      endif
      tile_count=tile_count+1
      if (tile_count==NSUB) then
        tile_count=0
# ifdef MPI
        size=NNODES
        do while(size > 1)
         step=(size+1)/2
          if (mynode>=step .and. mynode<size) then
            buff(1)=rx0
            buff(2)=rx1
            call MPI_Send (buff, 2, MPI_DOUBLE_PRECISION,
     &           mynode-step, 17, ocean_grid_comm,      ierr)
          elseif (mynode < size-step) then
            call MPI_Recv (buff, 2, MPI_DOUBLE_PRECISION,
     &        mynode+step, 17, ocean_grid_comm, status, ierr)
            rx0=max(rx0, buff(1))
            rx1=max(rx1, buff(2))
          endif
         size=step
        enddo

        buff(1)=rx0 ; buff(2)=rx1

        call MPI_Bcast(buff, 2, MPI_DOUBLE_PRECISION,
     &                      0, ocean_grid_comm, ierr)

        rx0=buff(1) ; rx1=buff(2)
# endif
        mpi_master_only write(*,'(/1x,A,F12.9,2x,A,F14.9/)')
     &   'Maximum grid stiffness ratios:   rx0 =',rx0, 'rx1 =',rx1
      endif
C$OMP END CRITICAL(grd_stff_cr_rgn)
      end
#else
      subroutine grid_stiffness_empty
      end
#endif  /* SOLVE3D */
