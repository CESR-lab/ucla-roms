/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

/*
 * CHOOSE ONLY ONE PRIMARY FLAG FOR SWITCH LIST BELOW
 */

c-dd#define WEC /* Wave Effect on Current model */
c-dd#define PACIFIC_PD /* PierreD's pacific coast model with tau-correction */
#define USWC_sample

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


#if defined PACIFIC_PD || defined USWC_sample /* PierreD's pacific coast model with tau-correction */

# define FLUX_FRC /* DevinD for new flux_frc.F module - should remove flag later */

c-dd# define WEC
# ifdef WEC
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
# endif /* WEC */

/* Include standard CPP switches for UP ETH Zurich */
c-dd#include "cppdefs_UP.h"

/*
   Standard UP ETH Zurich Settings for Regional and Basin Configurations
   #### PLEASE do not change but use undef in cppdefs.h to undefine
   #### what as needed e.g., put #undef SPONGE in your personal cppdefs.h
*/
                     /* Basics */
#define SOLVE3D
#define UV_ADV
#define UV_COR
                     /*  Equation of State */
#define NONLIN_EOS
#define SPLIT_EOS
#define SALINITY
c-dd#define SFLX_CORR    /* DevinD turned off as new 25km input wrong for sss */
                     /*  Forcing */
                     /*         - surface */
#define DIURNAL_SRFLUX /* Note this is 'undef'ed' below */
c-dd#define QCORRECTION /* DevinD no longer used for bulk force */
                     /*         - lateral */
#define T_FRC_BRY
#define Z_FRC_BRY
#define M3_FRC_BRY
#define M2_FRC_BRY
#define SPONGE
#define SPONGE_WIDTH /* # of sponge points is input parameter */
                     /* Mixing */
                     /*        - lateral */
#define UV_VIS2
#define TS_DIF2
                     /*        - vertical */
#define LMD_MIXING
#define LMD_KPP
#define LMD_NONLOCAL
#define LMD_RIMIX
#define LMD_CONVEC
#define LMD_BKPP

                      /* Grid Configuration */
#define CURVGRID
#define SPHERICAL
#define MASKING

                      /* Output Options */
#define MASK_LAND_DATA

                      /* Restart */
!--> #define EXACT_RESTART

                      /* Open Boundary Conditions */
#define OBC_M2FLATHER  /* Barotop. BC: OBC_M2FLATHER, OBC_M2ORLANSKI, OBC_M2SPECIFIED */
#define OBC_M3ORLANSKI /* Baroclin. BC: OBC_M3ORLANSKI, OBC_M3SPECIFIED */
#define OBC_TORLANSKI  /* Tracer BC: OBC_TORLANSKI, OBC_TSPECIFIED */

                      /* Biology Settings */
#ifdef BIOLOGY_BEC2
# define BIOLOGY
# define DAILYPAR_BEC
#endif
#ifdef BIOLOGY_NPZDOC
# define BIOLOGY
# define DAILYPAR_PHOTOINHIBITION
#endif


/* End of UP ETH Standard Settings */

/* Open Boundaries */
#define OBC_WEST /* Open boundary in the west (in order: SO out, SO in, Ind. throughflow) */
#define OBC_NORTH  /* Open boundary North (Arctic) */
#define OBC_EAST
#define OBC_SOUTH
#undef SPONGE /* DevinD - defined in cppdefs_UP.h */


/* Switches required for Flux correction */
c-dd#define SFLX_CORR ! Already defined in cppdefs_UP.h & DEVIND IN NEW CODE
#undef VFLX_CORR
#undef QCORRECTION
c-dd#define DQDT_DATA ! DevinD not entirely sure but dont think I need it
c-dd#define TAU_CORRECTION
#undef DIURNAL_SRFLUX

     /* Output */
#define AVERAGES
#undef SLICE_AVG
/* DPD CALENDER is in def_his.F of Pierre's code but not mine  */
c-dd#define CALENDAR '365_day'     /* netCDF CF-convention CALENDAR attribute default: '360_day' */
c-dd#define STARTDATE '1980-01-01' /* Ana's Hindcast - DPD: only in init_scalars.F seemingly for netcdf stamp */

#define ADV_ISONEUTRAL

     /* Biology */
c-dd#define BIOLOGY_BEC2

c-DDDD#define BULK_FLUX
c-dd#define BULK_SM_UPDATE ! DEVIND - REMOVED AS ALWAYS NEEDED
c-dd#define WIND_AT_RHO_POINTS  ! DEVIND - DEPRECATED IN NEW CODE
#define BULK_FLUX_OUTPUT /* DevinD added this for sustr and svstr outputs in new code */

    /* Flux Analysis */
c-dd# define WRITE_DEPTHS /* For Budget Analysis Closure */

    /* Tides */
c-dd# define TIDES
# ifdef TIDES
#  define POT_TIDES
#  define SSH_TIDES
#  define UV_TIDES
c-dd-gone#  define TIDERAMP ! No longer using tideramp
# endif


#elif defined WEC

/* WEC */

# undef ANA_WEC_FRC /* DEFINE FOR ANALYTICAL WEC */

/* DevinD - from DH's analytical wave_packet cppdefs.h */

# ifdef ANA_WEC_FRC
#  undef VERBOSE
#  define ANA_GRID
#  define ANA_INITIAL
#  define ANA_SMFLUX
#  define ANA_SRFLUX
#  define ANA_STFLUX
#  undef ANA_SST

#  undef SALINITY
#  define ANA_SSFLUX
#  define SOLVE3D
#  define UV_ADV
#  define NONLIN_EOS
#  define SPLIT_EOS

#  define EW_PERIODIC
#  define NS_PERIODIC

#  ifdef ANA_WEC_FRC
#   define ANA_WWAVE
#  endif
#  ifdef WEC
#   define ANA_WWAVE
#   define BRK0
#  endif

# endif /* ANA_WEC_FRC */

/* DevinD - end from DH's analytical wave_packet */

/* DevinD - start DH's non-analytical WEC */

# define USWC_WEC /* DEFINE FOR NON-ANALYTICAL WEC - US West Coast with WEC */

/* -------------------------------------------------------  */

# if defined USWC_WEC

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

#  define AVERAGES /* DevinD uncommented */
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
#  define BULK_FLUX
#  ifdef BULK_FLUX
c---#   define COUPLED_SURF_CURR /* not used in new code */
c---#   define WND_AT_RHO_POINTS /* Not needed in new code as wind converted to u/v */
#   define BULK_FLUX_OUTPUT /* DevinD output flux variables to ncdf */
#   define TAU_CORRECTION /* PierreD used to correct bulk flux towards measured data */
#  endif

# endif

/* DevinD - end DH's non-analytical WEC */

#endif


#include "set_global_definitions.h"

