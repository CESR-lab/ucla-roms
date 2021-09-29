      module mixing

#include "cppdefs.opt"

      use param
      implicit none

! This is include file "mixing"
!------ --- ----------------------
#ifdef UV_VIS2
      real visc2_r(GLOBAL_2D_ARRAY)
      real visc2_p(GLOBAL_2D_ARRAY)
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      real diff2(GLOBAL_2D_ARRAY,NT)
# endif
# ifdef TS_DIF4
      real diff4(GLOBAL_2D_ARRAY,NT)
# endif
# ifdef ADV_ISONEUTRAL
      real diff3u(GLOBAL_2D_ARRAY,N)
      real diff3v(GLOBAL_2D_ARRAY,N)
# endif

      real Akv(GLOBAL_2D_ARRAY,0:N)
# ifdef SALINITY
      real Akt(GLOBAL_2D_ARRAY,0:N,isalt)
# else
      real Akt(GLOBAL_2D_ARRAY,0:N,itemp)
# endif
# if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
  || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING\
  || defined LMD_BKPP
      real bvf(GLOBAL_2D_ARRAY,0:N)
# endif


# ifdef MY25_MIXING

! Mellor-Yamada (1982) Level 2.5 vertical mixing variables
! Akq     Vertical mixing coefficient [m^2/s] for TKE
! Lscale  Turbulent length scale (m).
! q2      Turbulent kinetic energy [m^2/s^2] at horizontal RHO-
!                                       and vertical W-points.
! q2l     TKE times turbulent length scale[m^3/s^2] at horizontal
!                                   RHO- and vertical W-points.

      real Akq(GLOBAL_2D_ARRAY,0:N)
      real Lscale(GLOBAL_2D_ARRAY,N)
      real q2(GLOBAL_2D_ARRAY,0:N,2)
      real q2l(GLOBAL_2D_ARRAY,0:N,2)
#endif /* MY25_MIXING */

! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

# ifdef LMD_KPP
      real hbls(GLOBAL_2D_ARRAY,2)

      real swr_frac(GLOBAL_2D_ARRAY,0:N)
#  ifdef LMD_NONLOCAL
      real ghat(GLOBAL_2D_ARRAY,0:N)
#  endif
# endif /* LMD_KPP */
# ifdef LMD_BKPP
      real hbbl(GLOBAL_2D_ARRAY,2)
# endif /* LMD_BKPP */
#endif /* SOLVE3D */

      end module mixing
