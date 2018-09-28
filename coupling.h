! File "coupling.h":  Declare 2D arrays associated with coupling
!----- -------------  between barotropic and baroclinic modes.
! They are divided into two groups, the upper one (above !> sign)
! is 3D --> 2D forcing terms (including both direct and parametric
! forcing). They are computed by the 3D part and used as input by
! step2D. The lower group is what barotropic mode returns to 3D:
! these are fast-time-averaged barotropic variables. 
!
#ifdef SOLVE3D
      real weight(2,288)
      common /coup_weight/ weight

      real rufrc(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE rufrc(BLOCK_PATTERN) BLOCK_CLAUSE
      real rvfrc(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE rvfrc(BLOCK_PATTERN) BLOCK_CLAUSE
      common /coup_rufrc/rufrc /coup_rvfrc/rvfrc

# ifdef PRED_COUPLED_MODE
      real rufrc_bak(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE rufrc_bak(BLOCK_PATTERN,2) BLOCK_CLAUSE
      real rvfrc_bak(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE rvfrc_bak(BLOCK_PATTERN,2) BLOCK_CLAUSE
      common /coup_rufrc_bak/rufrc_bak /coup_rvfrc_bak/rvfrc_bak
# endif

# ifdef VAR_RHO_2D
      real rhoA(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE rhoA(BLOCK_PATTERN) BLOCK_CLAUSE
      real rhoS(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE rhoS(BLOCK_PATTERN) BLOCK_CLAUSE
      common /coup_rhoA/rhoA /coup_rhoS/rhoS
# endif

      real r_D(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE Zt_avg1(BLOCK_PATTERN) BLOCK_CLAUSE
      common /coup_r_D/r_D

      real Zt_avg1(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE Zt_avg1(BLOCK_PATTERN) BLOCK_CLAUSE
      real DU_avg1(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE DU_avg1(BLOCK_PATTERN) BLOCK_CLAUSE
      real DV_avg1(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE DV_avg1(BLOCK_PATTERN) BLOCK_CLAUSE
      real DU_avg2(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE DU_avg2(BLOCK_PATTERN) BLOCK_CLAUSE
      real DV_avg2(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE DV_avg2(BLOCK_PATTERN) BLOCK_CLAUSE
      common /coup_Zt_avg1/Zt_avg1
     &       /coup_DU_avg1/DU_avg1 /coup_DV_avg1/DV_avg1
     &       /coup_DU_avg2/DU_avg2 /coup_DV_avg2/DV_avg2

# ifdef EXTRAP_BAR_FLUXES
      real DU_avg_bak(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE DU_avg_bak(BLOCK_PATTERN) BLOCK_CLAUSE
      real DV_avg_bak(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE DV_avg_bak(BLOCK_PATTERN) BLOCK_CLAUSE
      common /coup_DU_avg_bak/DU_avg_bak
     &       /coup_DV_avg_bak/DV_avg_bak
# endif
#endif /* SOLVE3D */ 
