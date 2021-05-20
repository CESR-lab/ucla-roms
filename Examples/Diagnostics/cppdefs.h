/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

/*
 * CHOOSE ONLY ONE PRIMARY FLAG FOR SWITCH LIST BELOW
 */

#define WECAGE /* Wave Effect on Current model */


/*
    Embedded (nested) grid configuration segment
*/

c--#ifndef MAX_GRID_LEVEL
c--# define MAX_GRID_LEVEL 2
c--# include "redefs.X"
c--#endif

/*
   Main switch starts here: model configuration choice.
*/

#if defined WECAGE

/* WEC */

/* DevinD - start DH's non-analytical WEC */

# define USWC_WEC /* DEFINE FOR NON-ANALYTICAL WEC - US West Coast with WEC */

/* -------------------------------------------------------  */

# if defined USWC_WEC

#  define DIAGNOSTICS
#  define DIAG3D_DEBUG
c--#  define DIAG2D_DEBUG

c--#  define BRY_EXTRACT_CHILD  /* required for online child bry interp */
#  undef GRID_ANG_DEG
#  define SOLVE3D
#  define UV_COR
#  define UV_ADV

#  define CURVGRID
#  define SPHERICAL
#  define MASKING

#  define SALINITY
#  define NONLIN_EOS
#  define SPLIT_EOS

#  define AVERAGES
c# define EXACT_RESTART
#  define NEW_S_COORD
c# define IMPLICIT_BOTTOM_DRAG

#  define UV_VIS2
#  define MIX_GP_UV
#  define TS_DIF2
#  define MIX_GP_TS

#  define LMD_MIXING
#  define LMD_KPP
#  define LMD_BKPP
#  define LMD_CONVEC
#  define LMD_NONLOCAL

#  define OBC_WEST
#  define OBC_EAST
#  define OBC_NORTH
#  define OBC_SOUTH

#  define OBC_M2FLATHER
#  define OBC_M3ORLANSKI
#  define OBC_TORLANSKI

#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
#  define SPONGE

#  undef WEC

#  ifdef WEC
#    define BRK0
#    define SURFACE_BREAK
#    undef SPEC_DD
#    undef LOG_BDRAG
#    undef WKB_WWAVE
#    undef BBL
#    define WAVE_OFFLINE
#    define WAVE_FRICTION
#    define BBL_S95
c#    define BBL_F00
#    define SUP_OFF
#    define WAVE_DIFF_FROM_LM
#  else
#   define LOG_BDRAG
#  endif
#  define BULK_FRC
#  ifdef BULK_FRC
#   define SURF_FLUX_OUTPUT_HIS /* DevinD output flux variables */
#   define SURF_FLUX_OUTPUT_AVG /* DevinD output flux variables as averages */
#  endif

# endif

/* DevinD - end DH's non-analytical WEC */

#endif


#include "set_global_definitions.h"

