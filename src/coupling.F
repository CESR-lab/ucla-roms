      module coupling

#include "cppdefs.opt"

      use param
      implicit none

! File "coupling":  Declare 2D arrays associated with coupling
!----- -------------  between barotropic and baroclinic modes.
! They are divided into two groups, the upper one (above !> sign)
! is 3D --> 2D forcing terms (including both direct and parametric
! forcing). They are computed by the 3D part and used as input by
! step2D. The lower group is what barotropic mode returns to 3D:
! these are fast-time-averaged barotropic variables. 
!
#ifdef SOLVE3D
      real weight(2,288)

      real rufrc(GLOBAL_2D_ARRAY)
      real rvfrc(GLOBAL_2D_ARRAY)

# ifdef VAR_RHO_2D
      real rhoA(GLOBAL_2D_ARRAY)
      real rhoS(GLOBAL_2D_ARRAY)
# endif

      real r_D(GLOBAL_2D_ARRAY)

      real Zt_avg1(GLOBAL_2D_ARRAY)
      real DU_avg1(GLOBAL_2D_ARRAY)
      real DV_avg1(GLOBAL_2D_ARRAY)
      real DU_avg2(GLOBAL_2D_ARRAY)
      real DV_avg2(GLOBAL_2D_ARRAY)

# ifdef EXTRAP_BAR_FLUXES
      real DU_avg_bak(GLOBAL_2D_ARRAY)
      real DV_avg_bak(GLOBAL_2D_ARRAY)
# endif
#endif /* SOLVE3D */ 

      end module coupling
