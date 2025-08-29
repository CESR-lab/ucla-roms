#include "cppdefs.opt"
! Copy the content of a private array into a shared array over
! extended range of indices.

      subroutine ext_copy_prv2shr_2d_tile (A,B)
      use param
      use dimensions
      implicit none
      integer :: i,j
      real A(PRIVATE_2D_SCRATCH_ARRAY)
      real B(GLOBAL_2D_ARRAY)

      ! GS: set from 0-ny/nx per sole usage of this routine for vars rd and r_D
      do j=0,ny
        do i=0,nx
          B(i,j)=A(i,j)
        enddo
      enddo
      end

#ifdef SOLVE3D
      subroutine ext_copy_prv2shr_tile (A,B,nmax)
      use param
      use dimensions
      implicit none
      integer :: nmax, i,j,k
      real A(PRIVATE_2D_SCRATCH_ARRAY,nmax)
      real B(GLOBAL_2D_ARRAY,nmax)

      ! GS NOTE: This appears to no longer be in use
      do k=1,nmax
        do j=0,ny
          do i=0,nx
            B(i,j,k)=A(i,j,k)
          enddo
        enddo
      enddo
      end

# define XI_ONLY
      subroutine ext_copy_prv2shr_1Dslb_tile (A,B,nmax)
      use param
      use dimensions
      implicit none
      integer :: nmax, i,j,k
      real A(PRIVATE_1D_SCRATCH_ARRAY,nmax)
      real B(GLOBAL_2D_ARRAY,nmax)

      ! GS NOTE: This appears to no longer be in use
      do k=1,nmax
        do i=0,nx
            B(i,j,k)=A(i,k)
        enddo
      enddo
      end
# undef XI_ONLY
#endif
