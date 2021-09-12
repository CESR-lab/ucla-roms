      module eos_vars

      ! Move to rho_eos when complete?

#include "cppdefs.h"

      use param
      implicit none

! Tcoef, T0   Coefficients for linear Equation of State 
! Scoef, S0     rho = Tcoef*(T-T0) + Scoef*(S-S0)
!
#ifdef SOLVE3D
# ifndef NONLIN_EOS
      real Tcoef, T0
#  ifdef SALINITY
      real Scoef, S0 
#  endif
# endif


# ifdef SPLIT_EOS
      real rho1(GLOBAL_2D_ARRAY,N)

      real qp1(GLOBAL_2D_ARRAY,N)
      real, parameter :: qp2=0.0000172
# else
      real rho(GLOBAL_2D_ARRAY,N)
# endif
# ifdef ADV_ISONEUTRAL
      real dRdx(GLOBAL_2D_ARRAY,N)
      real dRde(GLOBAL_2D_ARRAY,N)
      real idRz(GLOBAL_2D_ARRAY,0:N)
# endif
#endif

      end module eos_vars

