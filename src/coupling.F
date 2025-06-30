      module coupling

#include "cppdefs.opt"

      use param
      implicit none

! File "coupling":  Declare 2D arrays associated with coupling
!----- -------------  between barotropic and baroclinic modes.
! They are divided into two groups, the upper one (above !> sign)
! is 3D --> 2D forcing terms (including both direct and parametric
! forcing). They are computed by the 3D part and used as input by
! step2D. The lower group is what barotropic mode returns to 3D:
! these are fast-time-averaged barotropic variables.
!
#ifdef SOLVE3D
      real weight(2,288)

      real,allocatable,dimension(:,:) :: rufrc
      real,allocatable,dimension(:,:) :: rvfrc

# ifdef VAR_RHO_2D
      real,allocatable,dimension(:,:) :: rhoA
      real,allocatable,dimension(:,:) :: rhoS
# endif

      real,allocatable,dimension(:,:) :: r_D

      real,allocatable,dimension(:,:) :: Zt_avg1
      real,allocatable,dimension(:,:) :: DU_avg1
      real,allocatable,dimension(:,:) :: DV_avg1
      real,allocatable,dimension(:,:) :: DU_avg2
      real,allocatable,dimension(:,:) :: DV_avg2

# ifdef EXTRAP_BAR_FLUXES
      real,allocatable,dimension(:,:) :: DU_avg_bak
      real,allocatable,dimension(:,:) :: DV_avg_bak
# endif
#endif /* SOLVE3D */

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_coupling  ![
      use scalars
      implicit none

#ifdef SOLVE3D
      allocate( rufrc(GLOBAL_2D_ARRAY) ) ; rufrc=init
      allocate( rvfrc(GLOBAL_2D_ARRAY) ) ; rvfrc=init
# ifdef VAR_RHO_2D
      allocate( rhoA(GLOBAL_2D_ARRAY) ) ; rhoA=init
      allocate( rhoS(GLOBAL_2D_ARRAY) ) ; rhoS=init
# endif
      allocate( r_D(GLOBAL_2D_ARRAY) ) ; r_D=init
      allocate( Zt_avg1(GLOBAL_2D_ARRAY) ) ; Zt_avg1=init
      allocate( DU_avg1(GLOBAL_2D_ARRAY) ) ; DU_avg1=init
      allocate( DV_avg1(GLOBAL_2D_ARRAY) ) ; DV_avg1=init
      allocate( DU_avg2(GLOBAL_2D_ARRAY) ) ; DU_avg2=0.    ! avg2 need to be zero as used on RHS of set_huv1
      allocate( DV_avg2(GLOBAL_2D_ARRAY) ) ; DV_avg2=0.    ! before being set in step2d.
# ifdef EXTRAP_BAR_FLUXES
      allocate( DU_avg_bak(GLOBAL_2D_ARRAY) ) ; DU_avg_bak=init
      allocate( DV_avg_bak(GLOBAL_2D_ARRAY) ) ; DV_avg_bak=init
# endif
#endif /* SOLVE3D */

      end subroutine init_arrays_coupling  !]

      !----------------------------------------------------------------

      end module coupling
