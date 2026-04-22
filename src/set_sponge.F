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
      integer :: i,j,isp
      real    :: mult

      if (.not.bry_sponge) return

      ! TODO move sponge_size to an option file. 
      sponge_size = 15               ! Number of interior points in sponge layer (excludes buffer)

      isp = sponge_size+1

      if (mynode.eq.0) write(*,'(/1x,A,I3,A)')
     &  'Sponge layer = ',sponge_size,' interior grid points.'

      ! Since this code only considers sponge values in subdomains adjacent to a
      !  boundary, we're in trouble if the sponge region is wider than a subdomain.

      if ((isp>nx).or.(isp>ny)) then
        if (mynode.eq.0)
     &   print *, 'WARNING: sub-domains might be too small for the sponge layer'
      endif


      do i= 0,isp
        mult = float(isp-i)/float(isp)
        if (inode.eq.0.and.obc_west) then
# ifdef UV_VIS2
          visc2_r(i,1-bf:ny+bf)=visc2_r(i,1-bf:ny+bf)+v_sponge*mult
# endif
# ifdef TS_DIF2
          diff2(i,1-bf:ny+bf,:)=diff2(i,1-bf:ny+bf,:)+v_sponge*mult
# endif
        endif
        if (inode.eq.npx-1.and.obc_east) then
# ifdef UV_VIS2
          visc2_r(nx+1-i,1-bf:ny+bf)=visc2_r(nx+1-i,1-bf:ny+bf)+v_sponge*mult
# endif
# ifdef TS_DIF2
          diff2(nx+1-i,1-bf:ny+bf,:)=diff2(nx+1-i,1-bf:ny+bf,:)+v_sponge*mult
# endif
        endif
        if (jnode.eq.0.and.obc_south) then
# ifdef UV_VIS2
          visc2_r(1-bf:nx+bf,i)=visc2_r(1-bf:nx+bf,i)+v_sponge*mult
# endif
# ifdef TS_DIF2
          diff2(1-bf:nx+bf,i,:)=diff2(1-bf:nx+bf,i,:)+v_sponge*mult
# endif
        endif
        if (jnode.eq.npy-1.and.obc_north) then
# ifdef UV_VIS2
          visc2_r(1-bf:nx+bf,ny+1-i)=visc2_r(1-bf:nx+bf,ny+1-i)+v_sponge*mult
# endif
# ifdef TS_DIF2
          diff2(1-bf:nx+bf,ny+1-i,:)=diff2(1-bf:nx+bf,ny+1-i,:)+v_sponge*mult
# endif
        endif
      enddo

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
