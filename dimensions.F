      module dimensions

      ! ========================================
      ! Contains global and subdomain dimensions
      ! ========================================

      ! coded by Jeroen Molemaker

      implicit none
!      private

#include "cppdefs.h"
#include "param.h"
#include "grid.h"
#include "ncvars.h"

!      integer, public, parameter :: gnx = LLm
!      integer, public, parameter :: gny = MMm
!      integer, public, parameter ::  nx = Lm
!      integer, public, parameter ::  ny = Mm
!      integer, public, parameter ::  nz = N

      integer, parameter :: gnx = LLm
      integer, parameter :: gny = MMm
      integer, parameter ::  nx = Lm
      integer, parameter ::  ny = Mm
      integer, parameter ::  nz = N

!      public nt,frcfile,max_frc
!      public grdname
!
!# ifdef MASKING
!      public rmask
!# endif
!
!      public mynode

      end module dimensions
