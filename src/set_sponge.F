#include "cppdefs.opt"

      subroutine set_sponge

      use dimensions
      use mixing
      use scalars

      implicit none
      integer i,j, isp,ibnd
      real :: mult

      if (.not.bry_sponge) return

!     fill the full buffer to avoid mpi exchanges

      sponge_size = 15               ! Number of interior points in sponge layer (excludes buffer)
      isp = sponge_size+1

      if (mynode.eq.0) write(*,'(/1x,A,I3,A)')
     &  'Sponge layer = ',sponge_size,' interior grid points.'

      if ((isp>nx).or.(isp>ny)) then
        if (mynode.eq.0)
     &    print *, 'sub-domains might be too small for the sponge layer'
      endif


      if (inode.eq.0.and.obc_west) then
        do i= 0,isp
           mult = float(isp-i)/float(isp)
!          mult = 0.5 + 0.5*cos(pi*float(i)/float(isp))
#   ifdef UV_VIS2
          visc2_r(i,1-bf:ny+bf)=visc2_r(i,1-bf:ny+bf)+v_sponge*mult
#   endif
#   ifdef TS_DIF2
          diff2(i,1-bf:ny+bf,:)=diff2(i,1-bf:ny+bf,:)+v_sponge*mult
#   endif
        enddo
      endif
      if (inode.eq.npx-1.and.obc_east) then
        do i= 0,isp
           mult = float(isp-i)/float(isp)
#   ifdef UV_VIS2
           visc2_r(nx+1-i,1-bf:ny+bf)=visc2_r(nx+1-i,1-bf:ny+bf)+v_sponge*mult
#   endif
#   ifdef TS_DIF2
           diff2(nx+1-i,1-bf:ny+bf,:)=diff2(nx+1-i,1-bf:ny+bf,:)+v_sponge*mult
#   endif
        enddo
      endif
      if (jnode.eq.0.and.obc_south) then
        do j= 0,isp
           mult = float(isp-j)/float(isp)
#   ifdef UV_VIS2
           visc2_r(1-bf:nx+bf,j)=visc2_r(1-bf:ny+bf,j)+v_sponge*mult
#   endif
#   ifdef TS_DIF2
           diff2(1-bf:nx+bf,j,:)=diff2(1-bf:ny+bf,j,:)+v_sponge*mult
#   endif
        enddo
      endif
      if (jnode.eq.npy-1.and.obc_north) then
        do j= 0,isp
           mult = float(isp-j)/float(isp)
#   ifdef UV_VIS2
           visc2_r(1-bf:nx+bf,ny+1-j)=visc2_r(1-bf:ny+bf,ny+1-j)+v_sponge*mult
#   endif
#   ifdef TS_DIF2
           diff2(1-bf:nx+bf,ny+1-j,:)=diff2(1-bf:ny+bf,ny+1-j,:)+v_sponge*mult
#   endif
        enddo
      endif

#   ifdef UV_VIS2
      ! viscosity at psi points
      do j=1,ny+1
        do i=1,nx+1
          visc2_p(i,j)= 0.25*(visc2_r(i,j)+visc2_r(i-1,j)
     &                      + visc2_r(i,j-1)+visc2_r(i-1,j-1))
        enddo
      enddo
#   endif

      end subroutine set_sponge
