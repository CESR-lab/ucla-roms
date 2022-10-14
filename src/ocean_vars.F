
      module ocean_vars
      ! Ocean variables

#include "cppdefs.opt"

      use param
      use scalars
      implicit none
      private

      ! module preamble:  ![
      ! 2D - taken from what was ocean2d.F
      real,public,allocatable,dimension(:,:,:) :: zeta            ! free surface elevation [m] and barotropic
      real,public,allocatable,dimension(:,:,:) :: ubar            ! velocity components in XI-directions
      real,public,allocatable,dimension(:,:,:) :: vbar
      real,public,allocatable,dimension(:,:) :: zeta_avg
      real,public,allocatable,dimension(:,:) :: ubar_avg
      real,public,allocatable,dimension(:,:) :: vbar_avg
#ifdef NHMG
      real,public,allocatable,dimension(:,:) :: nh_ubar
      real,public,allocatable,dimension(:,:) :: nh_vbar
      real,public,allocatable,dimension(:,:) :: nh_wcor
#endif

      ! 3D - taken from what was ocean2d.F
#ifdef SOLVE3D
      real,public,allocatable,dimension(:,:,:,:) :: u
      real,public,allocatable,dimension(:,:,:,:) :: v
      real,public,allocatable,dimension(:,:,:)   :: u_avg
      real,public,allocatable,dimension(:,:,:)   :: v_avg
# if defined NHMG
      real,public,allocatable,dimension(:,:,:,:) :: w
      real,public,allocatable,dimension(:,:,:,:) :: nhdu
      real,public,allocatable,dimension(:,:,:,:) :: nhdv
      real,public,allocatable,dimension(:,:,:,:) :: nhdw
# endif

      real,public,allocatable,dimension(:,:,:) :: FlxU
      real,public,allocatable,dimension(:,:,:) :: FlxV
      real,public,allocatable,dimension(:,:,:) :: We      ! explicit
      real,public,allocatable,dimension(:,:,:) :: Wi      ! implicit

      real,public,allocatable,dimension(:,:,:) :: w_avg
      real,public,allocatable,dimension(:,:,:) :: wvl_avg

      real,public,allocatable,dimension(:,:,:) :: Hz      ! height of rho-cell
      real,public,allocatable,dimension(:,:,:) :: z_r     ! depth at rho-points
      real,public,allocatable,dimension(:,:,:) :: z_w     ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      real,public,allocatable,dimension(:,:,:) :: dzdxi
      real,public,allocatable,dimension(:,:,:) :: dzdeta
# endif
#endif  /* SOLVE3D */


      public :: init_arrays_ocean

      contains  !]
!----------------------------------------------------------------------
      subroutine init_arrays_ocean ![

      implicit none

      allocate( zeta(GLOBAL_2D_ARRAY,4) ); zeta=0.         ! zeta(:,:,knew) needs to be =0. for set_depth_tile
      allocate( ubar(GLOBAL_2D_ARRAY,4) ); ubar=0.         ! since knew can change if exact restart or not
      allocate( vbar(GLOBAL_2D_ARRAY,4) ); vbar=0.         ! set all zeta = 0.

#ifdef NHMG
      allocate( nh_ubar(GLOBAL_2D_ARRAY) )
      allocate( nh_vbar(GLOBAL_2D_ARRAY) )
      allocate( nh_wcor(GLOBAL_2D_ARRAY) )
#endif


#ifdef SOLVE3D
      allocate( u(GLOBAL_2D_ARRAY,N,3) )
      u(:,:,:,1)=init
      u(:,:,:,2)=0.0                             ! index 2 used on rhs u(indx) in pre_step for DC.
      u(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.
      allocate( v(GLOBAL_2D_ARRAY,N,3) )
      v(:,:,:,1)=init
      v(:,:,:,2)=0.0                             ! index 2 used on rhs v(indx) in pre_step for DC.
      v(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.

# if defined NHMG
      allocate( w(GLOBAL_2D_ARRAY,0:N,3) )
      allocate( nhdu(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdv(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdw(GLOBAL_2D_ARRAY,0:N,2) )
# endif

      allocate( FlxU(GLOBAL_2D_ARRAY,N) ) ; FlxU=init
      allocate( FlxV(GLOBAL_2D_ARRAY,N) ) ; FlxV=init
      allocate( We(GLOBAL_2D_ARRAY,0:N) ) ; We=init        ! explicit
      allocate( Wi(GLOBAL_2D_ARRAY,0:N) ) ; Wi=init        ! implicit

      allocate( Hz(GLOBAL_2D_ARRAY,N) )    ; Hz=init       ! height of rho-cell
      allocate( z_r(GLOBAL_2D_ARRAY,N) )   ; z_r=init      ! depth at rho-points
      allocate( z_w(GLOBAL_2D_ARRAY,0:N) ) ; z_w=init      ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      allocate( dzdxi(GLOBAL_2D_ARRAY,1:N)  )
      allocate( dzdeta(GLOBAL_2D_ARRAY,1:N) )
# endif
#endif  /* SOLVE3D */


      end subroutine init_arrays_ocean  !]
! ----------------------------------------------------------------------

      end module ocean_vars
