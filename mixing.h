! This is include file "mixing.h"
!------ --- ----------------------
#ifdef UV_VIS2
      real visc2_r(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE visc2_r(BLOCK_PATTERN) BLOCK_CLAUSE
      real visc2_p(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE visc2_p(BLOCK_PATTERN) BLOCK_CLAUSE
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      real diff2(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE diff2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /mixing_diff2/diff2
# endif
# ifdef TS_DIF4
      real diff4(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE diff4(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /mixing_diff4/diff4
# endif
# ifdef ADV_ISONEUTRAL
      real diff3u(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE diff3u(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real diff3v(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE diff3v(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_diffu/diff3u /ocean_diffv/diff3v
# endif

      real Akv(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE Akv(BLOCK_PATTERN,*) BLOCK_CLAUSE
# ifdef SALINITY
      real Akt(GLOBAL_2D_ARRAY,0:N,isalt)
# else
      real Akt(GLOBAL_2D_ARRAY,0:N,itemp)
# endif
CSDISTRIBUTE_RESHAPE Akt(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      common /mixing_Akv/Akv /mixing_Akt/Akt
# if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
  || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING\
  || defined LMD_BKPP
      real bvf(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE bvf(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /mixing_bvf/ bvf
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
CSDISTRIBUTE_RESHAPE Akq(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Lscale(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE Lscale(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real q2(GLOBAL_2D_ARRAY,0:N,2)
CSDISTRIBUTE_RESHAPE q2(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real q2l(GLOBAL_2D_ARRAY,0:N,2)
CSDISTRIBUTE_RESHAPE q2l(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      common /my25_mix_Akq/Akq /my25_mix_Lscale/Lscale
     &       /my25_mix_q2/q2   /my25_mix_q2l/q2l
#endif /* MY25_MIXING */

! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

# ifdef LMD_KPP
      real hbls(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE hbls(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /kpp_hbl/hbls

      real swr_frac(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE swr_frac(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /kpp_swr_frac/swr_frac
#  ifdef LMD_NONLOCAL
      real ghat(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE ghat(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /kpp_ghat/ghat
#  endif
# endif /* LMD_KPP */
# ifdef LMD_BKPP
      real hbbl(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE hbbl(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /kpp_hbbl/hbbl
# endif /* LMD_BKPP */
#endif /* SOLVE3D */

