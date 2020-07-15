/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

#define CALDEIRA /* Caldeira NH test */


/*
   Main switch starts here: model configuration choice.
*/

# ifdef CALDEIRA  /* Caldeira test */

# define MPI

# define NHMG
c---#  define NONTRAD_COR

#define EW_PERIODIC
c--# define OBC_EAST
c--# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH

# define ANA_BRY
# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
# endif

# define SOLVE3D
# define UV_ADV
# define UV_COR
# define SALINITY
# define NONLIN_EOS

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
                     /* OBCs algo */
c--# undef  OBC_M2ORLANSKI
c--# define OBC_M2SPECIFIED
# define OBC_M2FLATHER

# define OBC_M3SPECIFIED
c--# define  OBC_M3ORLANSKI

c--# define OBC_TSPECIFIED
# define  OBC_TORLANSKI

                      /* Sponge */
# undef SPONGE

                      /* Semi-implicit Vertical Tracer/Mom Advection */
c---# define  VADV_ADAPT_IMP

                      /* Vertical Mixing */
c---# define LINEAR_DRAG_ONLY
# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL

#endif /* END OF CONFIGURATION CHOICE */


#include "set_global_definitions.h"

