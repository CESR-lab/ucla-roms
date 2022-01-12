      module mixing

#include "cppdefs.opt"

      use param
      implicit none

! This is include file "mixing"
!------ --- ----------------------
#ifdef UV_VIS2
      real,allocatable,dimension(:,:) :: visc2_r
      real,allocatable,dimension(:,:) :: visc2_p
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      real,allocatable,dimension(:,:,:) :: diff2
# endif
# ifdef TS_DIF4
      real,allocatable,dimension(:,:,:) :: diff4
# endif
# ifdef ADV_ISONEUTRAL
      real,allocatable,dimension(:,:,:) :: diff3u
      real,allocatable,dimension(:,:,:) :: diff3v
# endif

      real,allocatable,dimension(:,:,:) :: Akv
# ifdef SALINITY
      real,allocatable,dimension(:,:,:,:) :: Akt
# else
      real,allocatable,dimension(:,:,:,:) :: Akt
# endif
# if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
  || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING\
  || defined LMD_BKPP
      real,allocatable,dimension(:,:,:) :: bvf
# endif


# ifdef MY25_MIXING

! Mellor-Yamada (1982) Level 2.5 vertical mixing variables
! Akq     Vertical mixing coefficient [m^2/s] for TKE
! Lscale  Turbulent length scale (m).
! q2      Turbulent kinetic energy [m^2/s^2] at horizontal RHO-
!                                       and vertical W-points.
! q2l     TKE times turbulent length scale[m^3/s^2] at horizontal
!                                   RHO- and vertical W-points.

      real,allocatable,dimension(:,:,:) :: Akq
      real,allocatable,dimension(:,:,:) :: Lscale
      real,allocatable,dimension(:,:,:,:) :: q2
      real,allocatable,dimension(:,:,:,:) :: q2l
#endif /* MY25_MIXING */

! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

# ifdef LMD_KPP
      real,allocatable,dimension(:,:,:) :: hbls

      real,allocatable,dimension(:,:,:) :: swr_frac
#  ifdef LMD_NONLOCAL
      real,allocatable,dimension(:,:,:) :: ghat
#  endif
# endif /* LMD_KPP */
# ifdef LMD_BKPP
      real,allocatable,dimension(:,:,:) :: hbbl
# endif /* LMD_BKPP */
#endif /* SOLVE3D */

      integer :: sponge_size  ! number of grid point of sponge layer

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_mixing  ![
      use scalars
      use tracers
      implicit none

      integer :: itrc

#ifdef UV_VIS2
      allocate( visc2_r(GLOBAL_2D_ARRAY) )
      allocate( visc2_p(GLOBAL_2D_ARRAY) )
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      allocate( diff2(GLOBAL_2D_ARRAY,NT) )
# endif
# ifdef TS_DIF4
      allocate( diff4(GLOBAL_2D_ARRAY,NT) )
# endif
# ifdef ADV_ISONEUTRAL
      allocate( diff3u(GLOBAL_2D_ARRAY,N) )
      allocate( diff3v(GLOBAL_2D_ARRAY,N) )
# endif

      allocate( Akv(GLOBAL_2D_ARRAY,0:N) )
# ifdef SALINITY
      allocate( Akt(GLOBAL_2D_ARRAY,0:N,isalt) )
# else
      allocate( Akt(GLOBAL_2D_ARRAY,0:N,itemp) )
# endif
# if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
  || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING\
  || defined LMD_BKPP
      allocate( bvf(GLOBAL_2D_ARRAY,0:N) )
# endif


# ifdef MY25_MIXING

! Mellor-Yamada (1982) Level 2.5 vertical mixing variables
! Akq     Vertical mixing coefficient [m^2/s] for TKE
! Lscale  Turbulent length scale (m).
! q2      Turbulent kinetic energy [m^2/s^2] at horizontal RHO-
!                                       and vertical W-points.
! q2l     TKE times turbulent length scale[m^3/s^2] at horizontal
!                                   RHO- and vertical W-points.

      allocate( Akq(GLOBAL_2D_ARRAY,0:N) )
      allocate( Lscale(GLOBAL_2D_ARRAY,N) )
      allocate( q2(GLOBAL_2D_ARRAY,0:N,2) )
      allocate( q2l(GLOBAL_2D_ARRAY,0:N,2) )
#endif /* MY25_MIXING */

! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

# ifdef LMD_KPP
      allocate( hbls(GLOBAL_2D_ARRAY,2) )

      allocate( swr_frac(GLOBAL_2D_ARRAY,0:N) )
#  ifdef LMD_NONLOCAL
      allocate( ghat(GLOBAL_2D_ARRAY,0:N) )
#  endif
# endif /* LMD_KPP */
# ifdef LMD_BKPP
      allocate( hbbl(GLOBAL_2D_ARRAY,2) )
# endif /* LMD_BKPP */
#endif /* SOLVE3D */

#ifdef UV_VIS2
      visc2_r=visc2
      visc2_p=visc2
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      do itrc=1,NT
        diff2(:,:,itrc)=tnu2(itrc)
      enddo
# endif

! Initialize vertical mixing coefficients (see "mixing") to their
! background values. If no vertical closure scheme is selected, the
! vertical mixing coefficients are those specified by the background
! values.

# if !defined LMD_MIXING && !defined BVF_MIXING && !defined PP_MIXING\
                         && !defined MY2_MIXING && !defined MY25_MIXING

      Akv=Akv_bak
# else
      Akv=0.
# endif
# ifdef AVERAGES
      akv_avg=0.
# endif
# if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
  || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING \
  || defined LMD_BKPP
      bvf=0.
# endif

      do itrc=1,iTandS
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined PP_MIXING\
                         && !defined MY2_MIXING && !defined MY25_MIXING
        Akt(:,:,:,itrc)=Akt_bak(itrc)
# else
        Akt(:,:,:,itrc)=0.
# endif
      enddo
# ifdef AVERAGES
      akt_avg=0.
#  ifdef SALINITY
      aks_avg=0.
#  endif
# endif

# if defined LMD_KPP && defined LMD_NONLOCAL
      ghat=0.
# endif
# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

      hbls=0.
#  ifdef AVERAGES
      hbl_avg=0.
#  endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

      hbbl=0.
#  ifdef AVERAGES
      hbbl_avg=0.
#  endif
# endif
#endif /* SOLVE3D */

      end subroutine init_arrays_mixing  !]

      !----------------------------------------------------------------

      end module mixing
