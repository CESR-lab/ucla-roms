      module mixing

#include "cppdefs.opt"

      use param
      use dimensions
      implicit none

      ! horizontal momentum mixing
      real,allocatable,dimension(:,:) :: visc2_r
      real,allocatable,dimension(:,:) :: visc2_p

      ! iso_neutral mixing vars ?
      real,allocatable,dimension(:,:,:) :: diff2
      real,allocatable,dimension(:,:,:) :: diff4
      real,allocatable,dimension(:,:,:) :: diff3u
      real,allocatable,dimension(:,:,:) :: diff3v

      ! vertical mixing vars
      real,allocatable,dimension(:,:,:) :: Akv
      real,allocatable,dimension(:,:,:) :: Akt

      real,allocatable,dimension(:,:,:) :: bvf

      ! currently only used in combination with kpp
      real,allocatable,dimension(:,:,:) :: swr_frac

      ! kpp related vars
      real,allocatable,dimension(:,:) :: hbls
      real,allocatable,dimension(:,:) :: hbbl
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

# ifdef TS_DIF2
      allocate( diff2(GLOBAL_2D_ARRAY,NT) ); diff2=init
# endif
# ifdef TS_DIF4
      allocate( diff4(GLOBAL_2D_ARRAY,NT) ); diff4=init
# endif
# ifdef ADV_ISONEUTRAL
      allocate( diff3u(GLOBAL_2D_ARRAY,N) ); diff3u=init
      allocate( diff3v(GLOBAL_2D_ARRAY,N) ); diff3v=init
# endif

      allocate( Akv(GLOBAL_2D_ARRAY,0:nz) ); Akv=init
      allocate( Akt(GLOBAL_2D_ARRAY,0:nz) ); Akt=init

# if defined LMD_MIXING
      allocate( bvf(GLOBAL_2D_ARRAY,0:nz) ); bvf=init
# endif


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

#ifdef UV_VIS2
      visc2_r=visc2
      visc2_p=visc2
#endif

# ifdef TS_DIF2
      do itrc=1,nt
        diff2(:,:,itrc)=tnu2(itrc)
      enddo
# endif

! Initialize vertical mixing coefficients (see "mixing") to their
! background values. If no vertical mixing closure scheme is selected, the
! vertical mixing coefficients are those specified by the background
! values.

# if defined LMD_MIXING
      Akv=0.
      bvf=0.
# else
      Akv=Akv_bak
      Akt=Akt_bak(1)
# endif

      end subroutine init_arrays_mixing  !]

      !----------------------------------------------------------------

      end module mixing
