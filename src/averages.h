/* This is include file "averages.h": time-averaged fields
 for output purposes:
*/
#ifdef AVERAGES
      real zeta_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE zeta_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real ubar_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  ubar_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      real vbar_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  vbar_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      common /avg_zeta/zeta_avg /avg_ubar/ubar_avg
     &                          /avg_vbar/vbar_avg
# ifdef SOLVE3D
      real u_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE u_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real v_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE v_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real t_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE t_avg(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real rho_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE rho_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real w_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE w_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_u/u_avg /avg_v/v_avg /avg_t/t_avg
     &                /avg_rho/rho_avg /avg_w/w_avg
      real wvl_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE wvl_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_wvl/wvl_avg
      real akv_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE akv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real akt_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE akt_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_akv/akv_avg /avg_akt/akt_avg
#  ifdef SALINITY
      real aks_avg(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE aks_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_aks/aks_avg
#  endif
#  ifdef LMD_KPP
      real hbl_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE hbl_avg(BLOCK_PATTERN) BLOCK_CLAUSE
      common /avg_hbl/hbl_avg
#  endif
#  ifdef LMD_BKPP
      real hbbl_avg(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE hbbl_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /avg_hbbl/hbbl_avg
#  endif
# endif
#endif /* AVERAGES */
 
