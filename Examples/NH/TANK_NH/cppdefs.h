/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

#define TANK
c---#define TANKINT

/*
   Main switch starts here: model configuration choice.
*/

#if defined TANK /* External NH tank test */
/*
!                       Tank Example
!                       ======= =======
!
! Chen, X.J., 2003. A fully hydrodynamic model for three-dimensional,
! free-surface flows.
! Int. J. Numer. Methods Fluids 42, 929–952.
*/
# define MPI

# define NHMG
c---# undef NHMG
# ifdef NHMG
#  define NHMG_WBRY_INIT
#  define NHMG_WBRY_COUPLING
#  define NHMG_WBRY_COPY
#  undef NHMG_WBRY_ZERO
#  undef NHMG_CHECKDIV
#  undef NHMG_DIAG
#  undef NHMG_2D_DAMPING
# endif
# define SOLVE3D
c---# define UV_ADV
# define LINEAR_DRAG_ONLY
# define NONLIN_EOS
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX

#elif defined TANKINT /* Internal NH tank test */
/*
!
!                       ===== ====== ===== =======
!
*/
# define MPI
# define NHMG
# ifdef NHMG
#  define NHMG_WBRY_INIT
#  define NHMG_WBRY_COUPLING
#  define NHMG_WBRY_COPY
#  undef NHMG_WBRY_ZERO
#  undef NHMG_CHECKDIV
#  undef NHMG_DIAG
#  undef NHMG_2D_DAMPING
# endif
# define SOLVE3D
c---# define UV_ADV
# define NONLIN_EOS
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX

# undef LMD_MIXING
# ifdef LMD_MIXING
#  undef LMD_SKPP
#  undef LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  undef  LMD_DDMIX
#  undef LMD_NONLOCAL
# endif

#endif /* END OF CONFIGURATION CHOICE */

#include "set_global_definitions.h"

