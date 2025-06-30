      module eos_vars

      ! Move to rho_eos when complete?

#include "cppdefs.opt"

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
      real,allocatable,dimension(:,:,:) :: rho1
      real,allocatable,dimension(:,:,:) :: qp1
      real, parameter :: qp2=0.0000172
# else
      real,allocatable,dimension(:,:,:) :: rho
# endif
# ifdef ADV_ISONEUTRAL
      real,allocatable,dimension(:,:,:) :: dRdx
      real,allocatable,dimension(:,:,:) :: dRde
      real,allocatable,dimension(:,:,:) :: idRz
# endif
#endif

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_eos_vars  ![
      use scalars
      implicit none

#ifdef SOLVE3D
# ifdef SPLIT_EOS
      allocate(  rho1(GLOBAL_2D_ARRAY,N) ); rho1=0.
      allocate( qp1(GLOBAL_2D_ARRAY,N) ); qp1=0.
# else
      allocate( rho(GLOBAL_2D_ARRAY,N) ); rho=0.
# endif
# ifdef ADV_ISONEUTRAL
      allocate( dRdx(GLOBAL_2D_ARRAY,N) ); dRdx=init
      allocate( dRde(GLOBAL_2D_ARRAY,N) ); dRde=init
      allocate( idRz(GLOBAL_2D_ARRAY,0:N) ); idRz=0.       ! -> loop ranges need fixing before init will work
# endif
#endif

      ! averaging variables allocated in ocean_vars to prevent circular reference from wrt_* logicals

      end subroutine  !]

!----------------------------------------------------------------------

      end module eos_vars

