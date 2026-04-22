#include "cppdefs.opt"

      subroutine set_sponge

!     Set sponge values near the open boundaries to reduce wave
!     trapping and provide a cleaner transition to the boundary data. 
!     Currently a linear function is used, but an alternative shape is:
!       mult = 0.5 + 0.5*cos(pi*float(i)/float(isp))
!     Extended ranges (filling the full buffer) are used to avoid mpi exchanges

      use dimensions
      use mixing
      use scalars

      implicit none
      integer :: i,j,isp,dist
      real,dimension(:,:),allocatable  :: mult

      if (.not.bry_sponge) return

      ! TODO move sponge_size to an option file. 
      sponge_size = 15               ! Number of interior points in sponge layer (excludes buffer)

      isp = sponge_size+1

      if (mynode.eq.0) write(*,'(/1x,A,I3,A)')
     &  'Sponge layer = ',sponge_size,' interior grid points.'

      allocate(mult(1-bf:nx+bf,1-bf:ny+bf))

      ! find minimum distance to physical boundary
      do j=1-bf,ny+bf
        do i=1-bf,nx+bf
          dist = isp
          if (obc_west)  dist = min(dist,i+iSW_corn)
          if (obc_east)  dist = min(dist,gnx+1-i-iSW_corn)
          if (obc_south) dist = min(dist,j+jSW_corn)
          if (obc_north) dist = min(dist,gny+1-j-jSW_corn)
          mult(i,j) = float(isp-dist)/float(isp)
        enddo
      enddo

      do j=1-bf,ny+bf
        do i=1-bf,nx+bf
# ifdef UV_VIS2
          visc2_r(i,j)=visc2_r(i,j)+v_sponge*mult(i,j)
# endif
# ifdef TS_DIF2
          diff2(i,j,:)=diff2(i,j,:)+v_sponge*mult(i,j)
# endif
        enddo
      enddo
      deallocate(mult)

# ifdef UV_VIS2
      ! viscosity at psi points
      do j=1,ny+1
        do i=1,nx+1
          visc2_p(i,j)= 0.25*(visc2_r(i,j)+visc2_r(i-1,j)
     &                      + visc2_r(i,j-1)+visc2_r(i-1,j-1))
        enddo
      enddo
# endif

      end subroutine set_sponge
