/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

#define CALDEIRA /* Caldeira NH test */


/*
   Main switch starts here: model configuration choice.
*/


# define MPI

# define NHMG
c---# undef NHMG
# ifdef NHMG
c---#  define NONTRAD_COR
#  define NHMG_WBRY_INIT
#  define NHMG_WBRY_COUPLING
#  define NHMG_WBRY_COPY
#  undef NHMG_WBRY_ZERO
# endif

# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH

# define ANA_BRY
# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
#  ifdef NHMG
#   define W_FRC_BRY
#  endif
# endif

# define SOLVE3D
# define UV_ADV
# define UV_COR
# define SALINITY
# define NONLIN_EOS
# define SPLIT_EOS

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
                     /* OBCs algo */
# undef  OBC_M2FLATHER
# undef  OBC_M2CHARACT
# undef  OBC_M2ORLANSKI
# undef  OBC_M3ORLANSKI
# undef  OBC_TORLANSKI
# define OBC_M2SPECIFIED
# define OBC_M3SPECIFIED
# define OBC_TSPECIFIED

                      /* Sponge */
# undef SPONGE

                      /* Semi-implicit Vertical Tracer/Mom Advection */
# define  VADV_ADAPT_IMP

                      /* Vertical Mixing */
# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL


#include "set_global_definitions.h"

