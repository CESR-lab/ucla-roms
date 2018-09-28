! These are global summation variables used exclussively by
! "diag.F" to compute and report various running diagnostics.

      real*QUAD  avzeta, avke,   prev_ke, avpe,
     &                   avke2b, avke3bc, avkesrf
      common /diag_vars/ avzeta, avke,    prev_ke,
     &             avpe, avke2b, avke3bc, avkesrf

#ifdef BIOLOGY
      real*QUAD global_sum(0:15)
      common /diag_vars/ global_sum
#endif

      real v2d_max
      common /diag_vars/ v2d_max
#ifdef SOLVE3D
# ifdef MAX_ADV_CFL
      real Cu_Adv,  Cu_W
      integer i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
# else
      real v3d_max, v3bc_max
      common /diag_vars/ v3d_max, v3bc_max
# endif
#endif

      
