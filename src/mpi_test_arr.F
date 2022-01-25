      module mpi_test_arr

#include "cppdefs.opt"

      use param
      implicit none

      real,allocatable,dimension(:,:)   :: A2D, B2D, C2D, D2D
      real,allocatable,dimension(:,:,:) :: A3D, B3D, C3D, D3D

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_mpi_test_arr  ![
      implicit none

      allocate( A2D(GLOBAL_2D_ARRAY),     B2D(GLOBAL_2D_ARRAY),
     &          C2D(GLOBAL_2D_ARRAY),     D2D(GLOBAL_2D_ARRAY)  )

      allocate( A3D(GLOBAL_2D_ARRAY,0:N), B3D(GLOBAL_2D_ARRAY,0:N),
     &          C3D(GLOBAL_2D_ARRAY,0:N), D3D(GLOBAL_2D_ARRAY,0:N)  )

      end subroutine init_arrays_mpi_test_arr  !]

!----------------------------------------------------------------------

      end module mpi_test_arr
