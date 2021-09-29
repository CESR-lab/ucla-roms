      module work_mod  ! was work .h but changed to work_mod for name clash with work array

#include "cppdefs.opt"

      use param
      implicit none

! This is "work": declaration of utility work array.
!
      real work(GLOBAL_2D_ARRAY,0:N)

      end module work_mod
 
