! This is include file "ocean2d.h".
!---------------------------------------------------------
! zeta       Free surface elevation [m] and barotropic
! ubar,vbar  velocity components in XI- and ETA-directions

      real zeta(GLOBAL_2D_ARRAY,4)
CSDISTRIBUTE_RESHAPE zeta(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ubar(GLOBAL_2D_ARRAY,4)
CSDISTRIBUTE_RESHAPE ubar(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vbar(GLOBAL_2D_ARRAY,4)
CSDISTRIBUTE_RESHAPE vbar(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_zeta/zeta /ocean_ubar/ubar /ocean_vbar/vbar

#ifdef NHMG
      real nh_ubar(GLOBAL_2D_ARRAY)
      real nh_vbar(GLOBAL_2D_ARRAY)
      real nh_wcor(GLOBAL_2D_ARRAY)
      common /nh_wcor/nh_ubar,nh_vbar,nh_wcor
#endif


 
