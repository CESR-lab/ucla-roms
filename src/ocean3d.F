      module ocean3d

#include "cppdefs.h"

      use param
      implicit none

! This is include file "ocean3d".
!----- -- ------- ---- ------------
#ifdef SOLVE3D
      real u(GLOBAL_2D_ARRAY,N,3)
      real v(GLOBAL_2D_ARRAY,N,3)
# if defined NHMG 
      real w(GLOBAL_2D_ARRAY,0:N,3)
      real nhdu(GLOBAL_2D_ARRAY,1:N,2)
      real nhdv(GLOBAL_2D_ARRAY,1:N,2)
      real nhdw(GLOBAL_2D_ARRAY,0:N,2)
# endif

      real FlxU(GLOBAL_2D_ARRAY,N)
      real FlxV(GLOBAL_2D_ARRAY,N)
      real We(GLOBAL_2D_ARRAY,0:N)  ! explicit
      real Wi(GLOBAL_2D_ARRAY,0:N)  ! implicit
      
      real Hz(GLOBAL_2D_ARRAY,N)    ! height of rho-cell
      real z_r(GLOBAL_2D_ARRAY,N)   ! depth at rho-points
      real z_w(GLOBAL_2D_ARRAY,0:N) ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      real dzdxi(GLOBAL_2D_ARRAY,1:N)
      real dzdeta(GLOBAL_2D_ARRAY,1:N)
# endif

#endif  /* SOLVE3D */

      end module ocean3d
