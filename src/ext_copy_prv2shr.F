#include "cppdefs.opt"
! Copy the content of a private array into a shared array over
! extended range of indices.

      subroutine ext_copy_prv2shr_2d_tile (istr,iend,jstr,jend, A,B)
      use param
      implicit none
      integer istr,iend,jstr,jend, i,j
      real A(PRIVATE_2D_SCRATCH_ARRAY)
      real B(GLOBAL_2D_ARRAY)

#include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          B(i,j)=A(i,j)
        enddo
      enddo
      end

#ifdef SOLVE3D
      subroutine ext_copy_prv2shr_tile (istr,iend,jstr,jend, A,B,nmax)
      use param
      implicit none
      integer istr,iend,jstr,jend, nmax, i,j,k
      real A(PRIVATE_2D_SCRATCH_ARRAY,nmax)
      real B(GLOBAL_2D_ARRAY,nmax)

# include "compute_extended_bounds.h"

      do k=1,nmax
        do j=jstrR,jendR
          do i=istrR,iendR
            B(i,j,k)=A(i,j,k)
          enddo
        enddo
      enddo
      end

# define XI_ONLY
      subroutine ext_copy_prv2shr_1Dslb_tile (istr,iend,j, A,B,nmax)
      use param
      implicit none
      integer istr,iend, nmax, i,j,k
      real A(PRIVATE_1D_SCRATCH_ARRAY,nmax)
      real B(GLOBAL_2D_ARRAY,nmax)

# include "compute_extended_bounds.h"

      do k=1,nmax
        do i=istrR,iendR
            B(i,j,k)=A(i,k)
        enddo
      enddo
      end
# undef XI_ONLY
#endif
