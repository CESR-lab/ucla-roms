#include "cppdefs.h"
#ifdef SOLVE3D
                                  ! This routine just initialised rhs
      subroutine prsgrd (tile)    ! arrays for 3D momentum equations to
      implicit none               ! zero. It is only useful for testing
      integer tile                ! or debugging purposes.
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd00_tile (istr,iend,jstr,jend, A3d(1,1),A3d(1,2))
      end

      subroutine prsgrd00_tile (istr,iend,jstr,jend, ru,rv)
      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
      do k=1,N
        do j=jstr,jend
          do i=istr,iend
            ru(i,j,k)=0.
            rv(i,j,k)=0.
          enddo
        enddo
      enddo
      end
#else
      subroutine prsgrd_empty
      end
#endif

