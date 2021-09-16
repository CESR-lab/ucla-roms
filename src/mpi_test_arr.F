      module mpi_test_arr

#include "cppdefs.h"

      use param
      implicit none

      real, dimension(GLOBAL_2D_ARRAY) :: A2D, B2D, C2D, D2D
      real, dimension(GLOBAL_2D_ARRAY,0:N) :: A3D, B3D, C3D, D3D

      end module mpi_test_arr
