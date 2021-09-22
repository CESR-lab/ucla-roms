/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

/*
 * CHOOSE ONLY ONE PRIMARY FLAG FOR SWITCH LIST BELOW
 */

#define TEST_TIDES

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

#if defined TEST_TIDES

# define TIDES
c# define RIVER_SOURCE
c# define PIPE_SOURCE
c# define ANA_RIVER_FRC

# ifdef ANA_RIVER_FRC
#  define ANA_GRID
#  define ANA_INITIAL
#  define ANA_SMFLUX
#  define ANA_SRFLUX
#  define ANA_STFLUX
#  define ANA_SSFLUX
# endif /* ANA_RIVER_FRC */

        /* Basics */
# define SOLVE3D
# define UV_ADV
# define UV_COR
        /* Equation of State */
# define NONLIN_EOS
# define SPLIT_EOS
# define SALINITY
        /* Mixing */
        /*        - lateral */
# define UV_VIS2
# define TS_DIF2
        /*        - vertical */
# define LMD_MIXING
# define LMD_KPP
# define LMD_NONLOCAL
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_BKPP
        /* Grid Configuration */
c-dd# define CURVGRID
c-dd# define SPHERICAL
# define MASKING
        /* Output Options */
# define MASK_LAND_DATA

#elif defined DUMMY_CASE

#  define AVERAGES

#endif


#include "set_global_definitions.h"

