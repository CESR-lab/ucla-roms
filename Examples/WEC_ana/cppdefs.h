/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

/* #define SPLASH - DevinD commented out */

/* DevinD added flags */
#define WEC
#define ANA_WEC_FRC

/* DevinD - from DH's wave_packet cppdefs.h */

# undef VERBOSE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# undef ANA_SST

# undef SALINITY
# define ANA_SSFLUX
# define SOLVE3D
# define UV_ADV
# define NONLIN_EOS
# define SPLIT_EOS

# define EW_PERIODIC
# define NS_PERIODIC

# define MRL_WCI
# ifdef MRL_WCI
#  define ANA_WWAVE
#  define BRK0
# endif


/* DevinD - end from DH's wave_packet cppdefs.h */




#include "set_global_definitions.h"

