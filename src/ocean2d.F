      module ocean2d

#include "cppdefs.h"

      use param
      implicit none

! This is include file "ocean2d".
!---------------------------------------------------------
! zeta       Free surface elevation [m] and barotropic
! ubar,vbar  velocity components in XI- and ETA-directions

      real zeta(GLOBAL_2D_ARRAY,4)
      real ubar(GLOBAL_2D_ARRAY,4)
      real vbar(GLOBAL_2D_ARRAY,4)

#ifdef NHMG
      real nh_ubar(GLOBAL_2D_ARRAY)
      real nh_vbar(GLOBAL_2D_ARRAY)
      real nh_wcor(GLOBAL_2D_ARRAY)
#endif

      end module ocean2d
 
