/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

#define SMC_BORE_SHELF      /* Submesocale + Bore Shelf Idealization */
/*
   Main switch starts here: model configuration choice.
*/
#if defined SMC_BORE_SHELF   
! Grid configuration
# define MASKING
c# define CURVGRID
c# define SPHERICAL
# undef ANA_GRID

! Momentum/tracer eqns
# define SOLVE3D
# define UV_ADV
# define UV_COR
c# define SALINITY

! Lateral viscosity/mixing
# define UV_VIS2
# define TS_DIF2

! Mixing Scheme
# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
# define LMD_CONVEC
# define LMD_NONLOCAL
# define LMD_RIMIX
# define MERGE_OVERLAP
c# define ANA_VMIX

!Atmospheric Forcing
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX

! Boundary Conditions
# define NS_PERIODIC
# define OBC_WEST
# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
# define T_FRC_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
c# define OBC_M2ORLANSKI
c# define NS_PERIODIC
c# define SPONGE
#endif

#include "set_global_definitions.h"

