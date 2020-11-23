/* This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ==== */

/*
 * CHOOSE ONLY ONE PRIMARY FLAG FOR SWITCH LIST BELOW
 */

#define PACIFIC_PD /* PierreD's pacific coast model with tau-correction */

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


#if defined PACIFIC_PD /* PierreD's pacific coast model with tau-correction */

/* Include standard CPP switches for UP ETH Zurich */
!#include "cppdefs_UP.h"

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
#define SFLX_CORR    /* DevinD turned off for bulk_frc.F to compile */
                     /*  Forcing */
                     /*         - surface */
#define DIURNAL_SRFLUX
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
#define TAU_CORRECTION
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

#define BULK_FLUX
c-dd#define BULK_SM_UPDATE ! DEVIND - REMOVED AS ALWAYS NEEDED
c-dd#define WIND_AT_RHO_POINTS  ! DEVIND - DEPRECATED IN NEW CODE
#define BULK_FLUX_OUTPUT /* DevinD added this for sustr and svstr outputs in new code */

    /* Flux Analysis */
c-dd# define WRITE_DEPTHS /* For Budget Analysis Closure */

    /* Tides */
# define TIDES
# ifdef TIDES
#  define POT_TIDES
#  define SSH_TIDES
#  define UV_TIDES
c-dd-gone#  define TIDERAMP ! No longer using tideramp
# endif

#endif


#include "set_global_definitions.h"

