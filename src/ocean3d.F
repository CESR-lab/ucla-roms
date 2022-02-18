      module ocean3d

#include "cppdefs.opt"

      use param
      implicit none

! This is include file "ocean3d".
!----- -- ------- ---- ------------
#ifdef SOLVE3D
      real,allocatable,dimension(:,:,:,:) :: u
      real,allocatable,dimension(:,:,:,:) :: v
      real,allocatable,dimension(:,:,:)   :: u_avg
      real,allocatable,dimension(:,:,:)   :: v_avg
# if defined NHMG 
      real,allocatable,dimension(:,:,:,:) :: w
      real,allocatable,dimension(:,:,:,:) :: nhdu
      real,allocatable,dimension(:,:,:,:) :: nhdv
      real,allocatable,dimension(:,:,:,:) :: nhdw
# endif

      real,allocatable,dimension(:,:,:) :: FlxU
      real,allocatable,dimension(:,:,:) :: FlxV
      real,allocatable,dimension(:,:,:) :: We      ! explicit
      real,allocatable,dimension(:,:,:) :: Wi      ! implicit
      
      real,allocatable,dimension(:,:,:) :: w_avg
      real,allocatable,dimension(:,:,:) :: wvl_avg

      real,allocatable,dimension(:,:,:) :: Hz      ! height of rho-cell
      real,allocatable,dimension(:,:,:) :: z_r     ! depth at rho-points
      real,allocatable,dimension(:,:,:) :: z_w     ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      real,allocatable,dimension(:,:,:) :: dzdxi
      real,allocatable,dimension(:,:,:) :: dzdeta
# endif
#endif  /* SOLVE3D */

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_ocean3d  ![
      use ncvars
      implicit none

      ! DevinD - again lots of extra first touched variables that weren't in init_arrays originally?

#ifdef SOLVE3D
      allocate( u(GLOBAL_2D_ARRAY,N,3) ) ; u=0.
      allocate( v(GLOBAL_2D_ARRAY,N,3) ) ; v=0.

      if (wrtavg(indxU)) then
        allocate( u_avg(GLOBAL_2D_ARRAY,N) )
        u_avg =0.
      endif
      if (wrtavg(indxV)) then 
        allocate( v_avg(GLOBAL_2D_ARRAY,N) )
        v_avg =0.
      endif

# if defined NHMG
      allocate( w(GLOBAL_2D_ARRAY,0:N,3) )
      allocate( nhdu(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdv(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdw(GLOBAL_2D_ARRAY,0:N,2) )
# endif

      allocate( FlxU(GLOBAL_2D_ARRAY,N) ) ; FlxU=0.
      allocate( FlxV(GLOBAL_2D_ARRAY,N) ) ; FlxV =0.
      allocate( We(GLOBAL_2D_ARRAY,0:N) ) ; We =0.         ! explicit
      allocate( Wi(GLOBAL_2D_ARRAY,0:N) ) ; Wi=0.          ! implicit

      if (wrtavg(indxO)) then
        allocate( w_avg(GLOBAL_2D_ARRAY,0:N) )
        w_avg = 0.
      endif
      if (wrtavg(indxW)) then
        allocate( wvl_avg(GLOBAL_2D_ARRAY,0:N) )
        wvl_avg = 0.
      endif


      allocate( Hz(GLOBAL_2D_ARRAY,N) )    ; Hz =0.        ! height of rho-cell
      allocate( z_r(GLOBAL_2D_ARRAY,N) )   ; z_r=0.        ! depth at rho-points
      allocate( z_w(GLOBAL_2D_ARRAY,0:N) ) ; z_w=0.        ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      allocate( dzdxi(GLOBAL_2D_ARRAY,1:N)  )
      allocate( dzdeta(GLOBAL_2D_ARRAY,1:N) )
# endif
#endif  /* SOLVE3D */

#ifdef SOLVE3D
      if (wrtavg(indxU)) u_avg=0.
      if (wrtavg(indxV)) v_avg=0.
#endif /* SOLVE3D */

      end subroutine init_arrays_ocean3d  !]

      !----------------------------------------------------------------

      end module ocean3d
