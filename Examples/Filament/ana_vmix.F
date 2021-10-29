#include "cppdefs.opt"
#if defined ANA_VMIX && defined SOLVE3D

      subroutine ana_vmix (tile)

      use param
      use hidden_mpi_vars

      implicit none
      integer tile
# include "compute_tile_bounds.h"
      call ana_vmix_tile (istr,iend,jstr,jend)
      end

      subroutine ana_vmix_tile (istr,iend,jstr,jend)

! Set vertical kinematic vertical viscosity "Akv" [m^2/s] for momentum
! and vertical mixing c efficients "Akt" [m^2/s] for tracer equations
! using analytical expressions.

      use param
      use grid
      use mixing
      use ocean3d
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real cff

# include "compute_auxiliary_bounds.h"

      if (FIRST_TIME_STEP) then
        do k=0,N
          do j=jstr,jend
            do i=istr,iend
             Akv(i,j,k)=0.0

             Akt(i,j,k,itemp)=0.
#  ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
            enddo
          enddo
        enddo
      endif


# ifdef EXCHANGE
      call exchange_2_tile (istr,iend,jstr,jend,  Akv, N+1,
     &                    Akt(START_2D_ARRAY,0,itemp), N+1)
#  ifdef SALINITY
      call exchange_tile (istr,iend,jstr,jend,
     &                    Akt(START_2D_ARRAY,0,isalt), N+1)
#  endif
# endif
      end
#else
      subroutine ana_vmix_empty
      end
#endif /* ANA_VMIX */
