      module work_mod  ! was work .h but changed to work_mod for name clash with work array

#include "cppdefs.opt"

      use param
      implicit none

! This is "work": declaration of utility work array.
!
      real,allocatable,dimension(:,:,:) :: work

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_work_mod  ![
      use scalars
      implicit none

      allocate( work(GLOBAL_2D_ARRAY,0:N) ); work=init

      end subroutine init_arrays_work_mod  !]
      !----------------------------------------------------------------

      end module work_mod

