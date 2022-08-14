      module dimensions

      ! ========================================
      ! Contains global and subdomain dimensions
      ! ========================================

      ! coded by Jeroen Molemaker
      use param

      implicit none
      private

#include "cppdefs.opt"

      integer, parameter, public :: gnx = LLm
      integer, parameter, public :: gny = MMm
      integer, parameter, public ::  nz = N

      integer, public            :: nx   ! Can vary between subdomains so not compile time parameters
      integer, public            :: ny

      public :: init_dimensions

      contains

! ----------------------------------------------------------------------
      subroutine init_dimensions  ![
      implicit none
      ! need to do this here because Lm, Mm change in mpi_setup and are not known at compile time.

      nx = Lm
      ny = Mm

      end subroutine init_dimensions  !]

! ----------------------------------------------------------------------

      end module dimensions
