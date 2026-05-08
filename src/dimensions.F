      module dimensions

      ! ========================================
      ! Contains global and subdomain dimensions
      ! ========================================

      ! coded by Jeroen Molemaker
      use param

      implicit none
      private

#include "cppdefs.opt"

      integer,parameter, public :: gnx = LLm
      integer,parameter, public :: gny = MMm
      integer,parameter, public ::  nz = N

      integer,parameter, public :: npx = np_xi
      integer,parameter, public :: npy = np_eta

      integer,parameter, public ::  bf = 2 ! allocation buffer size

      integer,public            :: nx   ! Can vary between subdomains so not compile time parameters
      integer,public            :: ny

      integer,public :: i0,i1,j0,j1     ! accounts for physical boundary buffers
      integer,public :: inode,jnode     ! location of subdomain in the grid

      integer, public ::  xi_rho, xi_u, eta_rho, eta_v

      logical,public :: nhmg = .false.  ! engage the mg based non-hydrostatic computations

# ifdef ANA_GRID
      logical,parameter,public :: analytical_grid = .true.
# else
      logical,parameter,public :: analytical_grid = .false.
# endif
      character(len=256), public :: grdname




      public :: init_dimensions

      contains

! ----------------------------------------------------------------------
      subroutine init_dimensions  ![
      implicit none
      ! need to do this here because Lm, Mm change in mpi_setup and are not known at compile time.

!     nx = Lm
!     ny = Mm
      Lm = nx
      Mm = ny

      i0=1; i1=nx; j0=1; j1=ny

      if (inode==0)        i0=0
      if (inode==NP_XI-1)  i1=nx+1
      if (jnode==0)        j0=0
      if (jnode==NP_ETA-1) j1=ny+1

      end subroutine init_dimensions  !]

! ----------------------------------------------------------------------

      end module dimensions
