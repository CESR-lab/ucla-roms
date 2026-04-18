      subroutine grid_stiffness

!     Survey three-dimensional grid in order to determine maximum
!     grid stiffness ratio:
!
!             z(i,j,k)-z(i-1,j,k)+z(i,j,k-1)-z(i-1,j,k-1)
!      rx1 = ---------------------------------------------
!             z(i,j,k)+z(i-1,j,k)-z(i,j,k-1)-z(i-1,j,k-1)
!
!      rx0 =  (z_w(i,0)-z_w(i-1,j,0) )/ (z_w(i,0)+z_w(i-1,j,0) )/ (
!
!     This is done for purely diagnostic purposes does not affect
!     computations.

      use dimensions
      use grid
      use ocean_vars
      use mpi
      use roms_mpi

      implicit none

!     local
      integer i,j,k
      real              :: rx0, rx1
      real,dimension(2) :: local,global
      integer ierr

      rx0=0. ; rx1=0.

      do j=1,ny
        do i=1,nx
          if (umask(i,j)>0.5) then
            rx0=max(rx0, abs( (z_w(i,j,0)-z_w(i-1,j,0))
     &                                   /(z_w(i,j,0)+z_w(i-1,j,0))
     &                                                            ))
            do k=1,nz
              rx1=max(rx1, abs(
     &         (z_w(i,j,k)-z_w(i-1,j,k)+z_w(i,j,k-1)-z_w(i-1,j,k-1))
     &        /(z_w(i,j,k)+z_w(i-1,j,k)-z_w(i,j,k-1)-z_w(i-1,j,k-1))
     &                                                            ))
            enddo
          endif
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if (vmask(i,j)>0.5) then
            rx0=max(rx0, abs( (z_w(i,j,0)-z_w(i,j-1,0))
     &                                   /(z_w(i,j,0)+z_w(i,j-1,0))
     &                                                            ))
            do k=1,nz
              rx1=max(rx1, abs(
     &         (z_w(i,j,k)-z_w(i,j-1,k)+z_w(i,j,k-1)-z_w(i,j-1,k-1))
     &        /(z_w(i,j,k)+z_w(i,j-1,k)-z_w(i,j,k-1)-z_w(i,j-1,k-1))
     &                                                            ))
            enddo
          endif
        enddo
      enddo

      local(1) = rx0
      local(2) = rx1
      call MPI_Allreduce(local, global, 2, MPI_DOUBLE_PRECISION,
     &                   MPI_MAX, ocean_grid_comm, ierr)
      rx0 = global(1)
      rx1 = global(2)

      if (mynode.eq.0) write(*,'(/1x,A,F12.9,2x,A,F14.9/)')
     &   'Maximum grid stiffness ratios:   rx0 =',rx0, 'rx1 =',rx1
      end
