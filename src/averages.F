      module averages
      ! time-averaged fields for output purposes

#include "cppdefs.h"

      use param
      implicit none

#ifdef AVERAGES
      real zeta_avg(GLOBAL_2D_ARRAY)
      real ubar_avg(GLOBAL_2D_ARRAY)
      real vbar_avg(GLOBAL_2D_ARRAY)
# ifdef SOLVE3D
      real u_avg(GLOBAL_2D_ARRAY,N)
      real v_avg(GLOBAL_2D_ARRAY,N)
      real rho_avg(GLOBAL_2D_ARRAY,N)
      real w_avg(GLOBAL_2D_ARRAY,0:N)
      real wvl_avg(GLOBAL_2D_ARRAY,0:N)
      real akv_avg(GLOBAL_2D_ARRAY,0:N)
      real akt_avg(GLOBAL_2D_ARRAY,0:N)
#  ifdef SALINITY
      real aks_avg(GLOBAL_2D_ARRAY,0:N)
#  endif
#  ifdef LMD_KPP
      real hbl_avg(GLOBAL_2D_ARRAY)
#  endif
#  ifdef LMD_BKPP
      real hbbl_avg(GLOBAL_2D_ARRAY)
#  endif
# endif
#endif /* AVERAGES */
 
      end module averages
