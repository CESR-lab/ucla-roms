/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

#undef BALTIC          /* Baltic Sea */
#undef DOUBLE_GYRE           /* Big Bad Basin */
#undef CANBAS2
#undef CANYON_A        /* Canyon_A */
#undef CANYON_B        /* Canyon_B */
#undef COLD_FILAMENT   /* Submesoscale cold filament */
#undef DAMEE_B         /* North Atlantic DAMEE, Big Domain */
#undef EKMAN_SPIRAL
#undef GRAV_ADJ        /* Graviational Adjustment */
#undef ISWAKE          /* Island Wake Problem */
#undef OVERFLOW        /* Graviational/Overflow */
c--#define BLACK_SEA
#undef ATLANTIC      /* Atlantic, Gulf of Mexico */
#undef PACIFIC         /* North-Equatorial Pacific */
#undef PACIFIC_2D      /* Pacific Tsunami model */
#undef RIVER           /* River runoff test problem */
#undef SEAMOUNT100        /* Seamount */
#undef SOLITON         /* Equatorial Rossby Wave */
#undef UPWELLING      /* Upwelling */
#undef USWC            /* US West Coast Yusuke 2010 */
#undef USWEST          /* US West Coast */
#undef WAVE_RAD        /* Test for wave radiation boundaries */
#define SPLASH

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

#if defined DOUBLE_GYRE    /* Big Bad Basin Configuration */
# define SOLVE3D

# define UV_ADV
# define UV_COR

# undef SALINITY
# undef NONLIN_EOS

# define UV_VIS2
# undef  TS_DIF2
# undef  TS_DIF4

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_VMIX

#elif defined CANYON_A      /* Canyon A Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define TS_DIF2
# define EW_PERIODIC
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX

#elif defined CANYON_B      /* Canyon B Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define TS_DIF2
# define EW_PERIODIC
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_VMIX

#elif defined DAMEE_B               /* North Atlantic DAMEE configuration */
# define SOLVE3D
c--# define AVERAGES
# define UV_COR
# define UV_ADV
!                       Equation of State
# define SALINITY
# define NONLIN_EOS
!                       Forcing and Climatology
# define TCLIMATOLOGY
# define TNUDGING
# define QCORRECTION
# define SFLX_CORR
!                       Lateral Mixing
# define VIS_GRID
# define DIF_GRID
!                       Vertical Mixing
# define LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
c--#  define LMD_DDMIX

c--#  define LMD_KPP
c--#  define LMD_NONLOCAL

!       Grid Configuration and Boundaries

# define CURVGRID
# define SPHERICAL
# define MASKING
# define EASTERN_WALL
# define WESTERN_WALL
# define SOUTHERN_WALL
# define NORTHERN_WALL


c--# define REST_STATE_TEST    /* Rest-state unforced problem */
# ifdef REST_STATE_TEST     /* (pressure gradient error test) */
#  define ANA_INITIAL

#  define SALINITY
#  define NONLIN_EOS

#  define ANA_SMFLUX
#  define ANA_SSFLUX
#  define ANA_STFLUX
#  define ANA_SST
#  define ANA_SRFLUX
#  undef TCLIMATOLOGY
#  undef TNUDGING
#  undef QCORRECTION
#  undef LMD_MIXING
#  define UV_VIS2
# endif


#elif defined EKMAN_SPIRAL     /* Ekman Spiral Test Problem */
# define SOLVE3D
# define UV_ADV
# define UV_COR

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX

# define EW_PERIODIC
# define NS_PERIODIC

#define STATIONS

#elif defined GRAV_ADJ     /* Gravitational Adjustment */
# define SOLVE3D

# define UV_ADV
# undef UV_COR
# define UV_VIS2

# define TS_DIF2

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX

c--# define OBC_WEST
c--# define OBC_EAST
c--# define OBC_M2ORLANSKI
c--# define OBC_M3ORLANSKI
c--# define OBC_TORLANSKI

#elif defined COLD_FILAMENT   /* Submesoscale cold filament */
# define SOLVE3D

c--# undef UV_ADV
c--# define CONST_TRACERS

# define UV_COR
# define UV_VIS2

# define TS_DIF2

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX

# define NS_PERIODIC
# define OBC_WEST
# define OBC_EAST
# define OBC_M2FLATHER
c--# define OBC_M2ORLANSKI
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI

# define ANA_BRY
# define T_FRC_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY


c--# define LMD_MIXING
c--# define LMD_KPP

#elif defined BLACK_SEA  /* Black Sea, an all-around closed basin */
# define SOLVE3D
# define UV_COR
c---# define NON_TRADITIONAL
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING
                        ! Equation of State
# define SALINITY
# define NONLIN_EOS
                       ! Lateral Mixing
# define UV_VIS2
# undef VIS_GRID

c--# define ADV_ISONEUTRAL
# define TS_DIF2
# undef DIF_GRID
                        ! Forcing
# define BULK_FLUX
# define BULK_SM_UPDATE
# define WIND_AT_RHO_POINTS

c--# define QCORRECTION
c--# define SFLX_CORR      /* salinity restoring at surface */
c--# define SEA_ICE_NOFLUX

# undef WIND_MAGN
                       ! Vertical Mixing
# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
c--# define DIURNAL_SRFLUX
# undef LMD_DDMIX


#elif defined ATLANTIC  /* Atlantic, Gulf of Mexico */
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define SALINITY
# define NONLIN_EOS

# define CURVGRID
# define SPHERICAL
# define MASKING

c--# define AVERAGES

c--# define ADV_ISONEUTRAL
# define UV_VIS2
# define TS_DIF2

# define WIND_AT_RHO_POINTS

# define DIURNAL_SRFLUX

c?? define ROBUST_DIURNAL_SRFLUX
# define QCORRECTION

# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL

# define OBC_WEST
# define OBC_EAST
# define OBC_SOUTH
c-# define OBC_NORTH

# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# define SPONGE

#elif defined SEAMOUNT100   /* North-Equatorial Pacific Configuration */
# define SOLVE3D
# define UV_COR
c---# define NON_TRADITIONAL
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING

c--# define AVERAGES
                        ! Equation of State
# define SALINITY
c--# define NONLIN_EOS
c--# define SPLIT_EOS
# define  TS_HADV_UP3
                       ! Lateral Mixing
# define UV_VIS2
# undef VIS_GRID

c--# define ADV_ISONEUTRAL
# define TS_DIF2
# undef DIF_GRID
                        ! Forcing 
c# define BULK_FLUX
c# define BULK_SM_UPDATE
c-# define WIND_AT_RHO_POINTS

c--# define QCORRECTION
c--# define SFLX_CORR     /* salinity restoring at surface */
c--# define SEA_ICE_NOFLUX

# undef WIND_MAGN
                       ! Vertical Mixing
# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
c--# define DIURNAL_SRFLUX

# undef LMD_DDMIX
                       ! Open Boundary Conditions
# define OBC_WEST
# define OBC_EAST
# define OBC_SOUTH
# define OBC_NORTH

# define OBC_M2FLATHER
c--# define OBC_M2ORLANSKI
c--# define OBC_M3ORLANSKI
c--# define OBC_TORLANSKI
# define  OBC_M3SPECIFIED
# define  OBC_TSPECIFIED
c--# define M2NUDGING
c--# define M3NUDGING

c>>># define TNUDGING
c>>># define TCLIMATOLOGY
c>>># define UCLIMATOLOGY

# define FRC_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY

# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define SPONGE


c--# define SSH_TIDES
c--# define UV_TIDES

#elif defined PACIFIC   /* North-Equatorial Pacific Configuration */
# define SOLVE3D
# define UV_COR
c---# define NON_TRADITIONAL
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING

c--# define AVERAGES
                        ! Equation of State
# define SALINITY
# define NONLIN_EOS
                       ! Lateral Mixing
# define UV_VIS2
# undef VIS_GRID

c--# define ADV_ISONEUTRAL
# define TS_DIF2
# undef DIF_GRID
                        ! Forcing 
c# define BULK_FLUX
c# define BULK_SM_UPDATE
c-# define WIND_AT_RHO_POINTS

# define QCORRECTION
# define SFLX_CORR     /* salinity restoring at surface */
# define SEA_ICE_NOFLUX

# undef WIND_MAGN
                       ! Vertical Mixing
# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
# define DIURNAL_SRFLUX

# undef LMD_DDMIX
                       ! Open Boundary Conditions
# define OBC_WEST
# define OBC_EAST
# define OBC_SOUTH
# define OBC_NORTH

# define OBC_M2FLATHER
c--# define OBC_M2ORLANSKI
c--# define OBC_M3ORLANSKI
c--# define OBC_TORLANSKI
# define  OBC_M3SPECIFIED
# define  OBC_TSPECIFIED
# define M2NUDGING
# define M3NUDGING

c>>># define TNUDGING
c>>># define TCLIMATOLOGY
c>>># define UCLIMATOLOGY

# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY

# define SPONGE


c--# define SSH_TIDES
c--# define UV_TIDES


#elif defined PACIFIC_2D   /* Pacific Tsunami model */
# define UV_COR
# define UV_ADV
# define CURVGRID
# define SPHERICAL
# define MASKING
# undef VIS_GRID
# define UV_VIS2
# define ANA_SMFLUX
# define ANA_INITIAL

# define OBC_WEST
# define OBC_SOUTH
# define OBC_M2FLATHER
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
c--# define OBC_M2ORLANSKI
c--# define OBC_VOLCONS
# define SPONGE

#elif defined OVERFLOW      /* Gravitational Overflow */
# define SOLVE3D

# define UV_ADV
# define UV_VIS2

# define TS_DIF2

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX

#elif defined SEAMOUNT     /* Seamount Configuration */
# define SOLVE3D

# define UV_ADV
# define UV_COR
c--# define UV_VIS2

# undef TS_DIF2
# undef  TS_DIF4

c--# define SALINITY
c--# define NONLIN_EOS
c--# define SPLIT_EOS

c--# define EW_PERIODIC
c--# define NS_PERIODIC
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX


#elif defined ISWAKE     /* Island Wake Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR

# define UV_VIS2
c---# define SPONGE
# define LINEAR_DRAG_ONLY

c--# define OBC_NORTH

c--# define OBC_M2ORLANSKI
c--# define OBC_M3ORLANSKI
c--# define OBC_TORLANSKI
c--# define OBC_M2SPECIFIED
# define OBC_M3SPECIFIED
# define OBC_TSPECIFIED

#define OBC_M2FLATHER

# define ANA_BRY
# define T_FRC_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY

c--# define AVERAGES

# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX

! Vertical mixing: nothing is defined here => use
! externally supplied constant background value.

#elif defined SOLITON    /* Equatorial Rossby Soliton */
# undef  SOLVE3D

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX

# define UV_COR
# define UV_ADV
# undef UV_VIS2
# define EW_PERIODIC
c--# define NS_PERIODIC

c--# define OBC_WEST
c--# define OBC_EAST
c--# define OBC_NORTH
c--# define OBC_SOUTH
c--# define OBC_M2ORLANSKI

#elif defined RIVER     /* River run-off test problem */
# define SOLVE3D
# define UV_ADV
# define UV_COR

# define SALINITY
# define NONLIN_EOS

# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX

# define PSOURCE
# define ANA_PSOURCE

# define EASTERN_WALL
# define WESTERN_WALL
# define NORTHERN_WALL
# define OBC_SOUTH
# define OBC_TORLANSKI
# define OBC_M2ORLANSKI
# define OBC_M3ORLANSKI

#elif defined UPWELLING     /* Upwelling Configuration */
# define SOLVE3D
# define UV_ADV
# define UV_COR

c--# define ADV_ISONEUTRAL

# define SALINITY
# undef NONLIN_EOS

# define EW_PERIODIC

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX

# undef ANA_VMIX


# ifndef ANA_VMIX
#  define LMD_MIXING
#  define LMD_KPP
#  define LMD_BKPP
#  define LMD_RIMIX
c--# define LMD_CONVEC
c--# define LMD_NONLOCAL
# endif





#elif defined CANBAS2   /* Canary Basin model */

# define SOLVE3D
# define UV_ADV
# define UV_COR

# define SALINITY
# define NONLIN_EOS

# define AVERAGES

# define UV_VIS2
# define TS_DIF2

# define QCORRECTION

# define LMD_MIXING
# define LMD_KPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL

# define CURVGRID
# define SPHERICAL
# define MASKING

# define OBC_WEST
# define OBC_EAST
# define OBC_NORTH
# define OBC_SOUTH

# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI

# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY

#elif defined BALTIC  /* Baltic Sea model */

# define SOLVE3D
# define UV_ADV
# define UV_COR

# define SALINITY
# define NONLIN_EOS

# undef AVERAGES

# define UV_VIS2
# define TS_DIF2


# define LMD_MIXING
# define LMD_KPP
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL

# define CURVGRID
# define SPHERICAL
# define MASKING

# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX

#elif defined USWC   /* US West Coast 2010 */
# define SOLVE3D
# define UV_ADV
# define UV_COR

# define SALINITY
# define NONLIN_EOS

# define AVERAGES

# define UV_VIS2
# define MIX_GP_UV
# define TS_DIF2

c------# define DIURNAL_SRFLUX
# undef TIDES
# ifdef TIDES
#  define SSH_TIDES
#  define UV_TIDES
c#  define TIDERAMP
# endif

# define QCORRECTION

# define LMD_MIXING
# define LMD_KPP
# define LMD_BKPP
c# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL

# define CURVGRID
# define SPHERICAL
# define MASKING

# define OBC_WEST
# define OBC_EAST
# undef OBC_NORTH
# define OBC_SOUTH

# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI

# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# define SPONGE

#elif defined SPLASH   /* US West Coast Configuration */
# define SOLVE3D
# define UV_COR
# define UV_ADV
                       /* Equation of State */
# define SALINITY
# define NONLIN_EOS
                       /* Forcing and Climatology */
# define QCORRECTION
# define SPONGE
# define SFLX_CORR
                      /* Lateral viscosity/mixing  */
# define UV_VIS2
# undef VIS_GRID
# define TS_DIF2
# undef DIF_GRID
                      /* Vertical Mixing */
# define LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# undef LMD_DDMIX
# define LMD_KPP
# define LMD_NONLOCAL
                      /* Grid Configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
                      /* Open Boundary Conditions */
# define OBC_EAST
# define OBC_WEST
c---# define OBC_NORTH
# define OBC_SOUTH
c--> # define OBC_VOLCONS
c--> # define OBC_FLUX_CORR

# undef OBC_M2ORLANSKI
# define OBC_M2FLATHER
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI

#define BRY
#ifdef BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# undef TNUDGING
# undef M3NUDGING
# undef M2NUDGING
# undef CLIMAT_TS_MIXH
#else
# define TNUDGING
# define M3NUDGING
# define M2NUDGING
#endif

# undef OBC_TSPECIFIED
# undef OBC_M2SPECIFIED
# undef OBC_M3SPECIFIED

# define AVERAGES

#elif defined WAVE_RAD
# undef  SOLVE3D
# define UV_COR
# define UV_ADV
# undef UV_VIS2

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX

# define OBC_WEST
# define OBC_EAST
# define OBC_NORTH
# define OBC_SOUTH
# define OBC_M2FLATHER
c--# define OBC_M2ORLANSKI
c--# define OBC_Z2ORLANSKI
# define LEGACY_ORLANSKI

# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY


#elif defined SLOSH
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
#endif

#include "set_global_definitions.h"

