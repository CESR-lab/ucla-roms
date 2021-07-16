! This is include file "ocean3d.h".
!----- -- ------- ---- ------------
#ifdef SOLVE3D
      real u(GLOBAL_2D_ARRAY,N,3)
CSDISTRIBUTE_RESHAPE u(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real v(GLOBAL_2D_ARRAY,N,3)
CSDISTRIBUTE_RESHAPE v(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
! 't' has been moved to tracers.F module since it is a tracer
!      real t(GLOBAL_2D_ARRAY,N,3,NT)
!CSDISTRIBUTE_RESHAPE t(BLOCK_PATTERN,*,*,*) BLOCK_CLAUSE
      common /ocean_u/u /ocean_v/v ! /ocean_t/t
# if defined NHMG 
      real w(GLOBAL_2D_ARRAY,0:N,3)
CSDISTRIBUTE_RESHAPE w(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real nhdu(GLOBAL_2D_ARRAY,1:N,2)
      real nhdv(GLOBAL_2D_ARRAY,1:N,2)
      real nhdw(GLOBAL_2D_ARRAY,0:N,2)
      common /ocean_w/w
      common /ocean_nhdu/nhdu
      common /ocean_nhdv/nhdv
      common /ocean_nhdw/nhdw
# endif

      real FlxU(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE  FlxU(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real FlxV(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE  FlxV(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real We(GLOBAL_2D_ARRAY,0:N)  ! explicit
CSDISTRIBUTE_RESHAPE We(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Wi(GLOBAL_2D_ARRAY,0:N)  ! implicit
CSDISTRIBUTE_RESHAPE Wi(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /flx_FU/FlxU /flx_FV/FlxV /flx_We/We /flx_Wi/Wi
      
      real Hz(GLOBAL_2D_ARRAY,N)    ! height of rho-cell
CSDISTRIBUTE_RESHAPE Hz(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real z_r(GLOBAL_2D_ARRAY,N)   ! depth at rho-points
CSDISTRIBUTE_RESHAPE z_r(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real z_w(GLOBAL_2D_ARRAY,0:N) ! depth at   w-points
CSDISTRIBUTE_RESHAPE z_w(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /grid_zw/z_w /grid_zr/z_r /grid_Hz/Hz
# if defined NHMG || defined NONTRAD_COR
      real dzdxi(GLOBAL_2D_ARRAY,1:N)
      real dzdeta(GLOBAL_2D_ARRAY,1:N)
      common /ocean_dzdxi/dzdxi
      common /ocean_dzdeta/dzdeta
# endif


#endif  /* SOLVE3D */
