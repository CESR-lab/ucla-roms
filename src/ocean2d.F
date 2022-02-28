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

      real,allocatable,dimension(:,:) :: zeta_avg
      real,allocatable,dimension(:,:) :: ubar_avg
      real,allocatable,dimension(:,:) :: vbar_avg

#ifdef NHMG
      real,allocatable,dimension(:,:) :: nh_ubar
      real,allocatable,dimension(:,:) :: nh_vbar
      real,allocatable,dimension(:,:) :: nh_wcor
#endif

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_ocean2d  ![
      use ncvars
      use scalars

      implicit none

      allocate( zeta(GLOBAL_2D_ARRAY,4) ); zeta=0.         ! zeta(:,:,knew) needs to be =0. for set_depth_tile
      allocate( ubar(GLOBAL_2D_ARRAY,4) ); ubar=0.         ! since knew can change if exact restart or not
      allocate( vbar(GLOBAL_2D_ARRAY,4) ); vbar=0.         ! set all zeta = 0.

      if (wrtavg(indxZ)) then
        allocate( zeta_avg(GLOBAL_2D_ARRAY) )
        zeta_avg=init
      endif
      if (wrtavg(indxUb)) then
        allocate( ubar_avg(GLOBAL_2D_ARRAY) )
        ubar_avg=init
      endif
      if (wrtavg(indxVb)) then
        allocate( vbar_avg(GLOBAL_2D_ARRAY) )
        vbar_avg=init
      endif

#ifdef NHMG
      allocate( nh_ubar(GLOBAL_2D_ARRAY) )
      allocate( nh_vbar(GLOBAL_2D_ARRAY) )
      allocate( nh_wcor(GLOBAL_2D_ARRAY) )
#endif

      end subroutine init_arrays_ocean2d  !]

      !----------------------------------------------------------------

      end module ocean2d
 
