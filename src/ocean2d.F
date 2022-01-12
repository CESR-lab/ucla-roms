      module ocean2d

#include "cppdefs.opt"

      use param
      implicit none

! This is include file "ocean2d".
!---------------------------------------------------------
! zeta       Free surface elevation [m] and barotropic
! ubar,vbar  velocity components in XI- and ETA-directions

      real,allocatable,dimension(:,:,:) :: zeta
      real,allocatable,dimension(:,:,:) :: ubar
      real,allocatable,dimension(:,:,:) :: vbar

#ifdef NHMG
      real,allocatable,dimension(:,:) :: nh_ubar
      real,allocatable,dimension(:,:) :: nh_vbar
      real,allocatable,dimension(:,:) :: nh_wcor
#endif

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_ocean2d  ![
      implicit none

      allocate( zeta(GLOBAL_2D_ARRAY,4) )
      allocate( ubar(GLOBAL_2D_ARRAY,4) )
      allocate( vbar(GLOBAL_2D_ARRAY,4) )

#ifdef NHMG
      allocate( nh_ubar(GLOBAL_2D_ARRAY) )
      allocate( nh_vbar(GLOBAL_2D_ARRAY) )
      allocate( nh_wcor(GLOBAL_2D_ARRAY) )
#endif

      zeta=0.
      ubar=0.
      vbar=0.
#ifdef AVERAGES
      zeta_avg=0.
      ubar_avg=0.
      vbar_avg=v
#endif

      end subroutine init_arrays_ocean2d  !]

      !----------------------------------------------------------------

      end module ocean2d
 
