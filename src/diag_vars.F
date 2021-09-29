      module diag_vars

#include "cppdefs.opt"

      implicit none

      ! This should end up in a diag.F module.

! These are global summation variables used exclussively by
! "diag.F" to compute and report various running diagnostics.

      real*QUAD  avzeta, avke,   prev_ke, avpe,
     &                   avke2b, avke3bc, avkesrf

      real v2d_max
#ifdef SOLVE3D
!# ifdef MAX_ADV_CFL
      real Cu_Adv,  Cu_W
      integer i_cx_max, j_cx_max, k_cx_max
!# else
      real v3d_max, v3bc_max
!# endif
#endif

      end module diag_vars
