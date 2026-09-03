      module comm_vars

#include "cppdefs.opt"

      use param
      implicit none

!   aim to remove this module

      ! computed in setup_grid2
      real area, volume !, bc_crss
      real hmin,hmax, grdmin,grdmax, Cg_min,Cg_max, Cu_Cor

#if defined BIOLOGY || defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC \
                    || defined BIOLOGY_BEC2 || defined MARBL
      real*QUAD global_sum(0:2*NT+1), global_srf_sum(0:NT)
#endif 
      real*4 cpu_all(2)
      integer trd_count, tile_count, bc_count, mcheck, first_time

      end module comm_vars
