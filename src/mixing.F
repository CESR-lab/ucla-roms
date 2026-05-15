      module mixing

#include "cppdefs.opt"

      use param
      use dimensions
      implicit none

! This is include file "mixing"
!------ --- ----------------------
#ifdef UV_VIS2
      real,allocatable,dimension(:,:) :: visc2_r
      real,allocatable,dimension(:,:) :: visc2_p
#endif

      real,allocatable,dimension(:,:,:) :: diff2
      real,allocatable,dimension(:,:,:) :: diff4

# ifdef ADV_ISONEUTRAL
      real,allocatable,dimension(:,:,:) :: diff3u
      real,allocatable,dimension(:,:,:) :: diff3v
# endif

      real,allocatable,dimension(:,:,:) :: Akv
      real,allocatable,dimension(:,:,:) :: Akt
      real,allocatable,dimension(:,:,:) :: bvf

! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

      real,allocatable,dimension(:,:) :: hbls
      real,allocatable,dimension(:,:) :: hbbl
      real,allocatable,dimension(:,:,:) :: swr_frac
      real,allocatable,dimension(:,:,:) :: ghat

      integer :: sponge_size  ! number of grid point of sponge layer

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_mixing  ![
      use scalars
      use tracers
      implicit none

      integer :: itrc

#ifdef UV_VIS2
      allocate( visc2_r(GLOBAL_2D_ARRAY) ); visc2_r=init
      allocate( visc2_p(GLOBAL_2D_ARRAY) ); visc2_p=init
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      allocate( diff2(GLOBAL_2D_ARRAY,nt) ); diff2=init
# endif

# ifdef ADV_ISONEUTRAL
      allocate( diff3u(GLOBAL_2D_ARRAY,nz) ); diff3u=init
      allocate( diff3v(GLOBAL_2D_ARRAY,nz) ); diff3v=init
# endif

      allocate( Akv(GLOBAL_2D_ARRAY,0:nz) ); Akv=init
      allocate( Akt(GLOBAL_2D_ARRAY,0:nz) ); Akt=init
# if defined LMD_KPP || defined LMD_BKPP
      allocate( bvf(GLOBAL_2D_ARRAY,0:nz) ); bvf=0
# endif


! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

# ifdef LMD_KPP
      allocate( hbls(GLOBAL_2D_ARRAY) ); hbls=0.

      allocate( swr_frac(GLOBAL_2D_ARRAY,0:N) ); swr_frac=init
#  ifdef LMD_NONLOCAL
      allocate( ghat(GLOBAL_2D_ARRAY,0:N) ); ghat=0.
#  endif
# endif /* LMD_KPP */
# ifdef LMD_BKPP
      allocate( hbbl(GLOBAL_2D_ARRAY) ); hbbl=0.
# endif /* LMD_BKPP */
#endif /* SOLVE3D */

#ifdef UV_VIS2
      visc2_r=visc2
      visc2_p=visc2
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      do itrc=1,nt
        diff2(:,:,itrc)=tnu2(itrc)
      enddo
# endif

! Initialize vertical mixing coefficients (see "mixing") to their
! background values. If no vertical closure scheme is selected, the
! vertical mixing coefficients are those specified by the background
! values.

# if !defined LMD_MIXING
      Akv=Akv_bak
# else
      Akv=0.
# endif

# if !defined LMD_MIXING 
      Akt = Akt_bak
# else
      Akt = 0.
# endif

      ! averaging variables allocated in ocean_vars to prevent circular reference from wrt_* logicals

#endif /* SOLVE3D */

      end subroutine init_arrays_mixing  !]

      !----------------------------------------------------------------

      end module mixing
