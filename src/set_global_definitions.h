/* This is "global_definitions.h":  A set of predetermined CPP macros
 and switches inserted into other files of the code during compilation
 process. General user is strongly discouraged from attempts to modify
 anything below this line.
 ------------------------------------------------------------------ */

/* Turn ON/OFF MPI parallelization. If not activated,the code becomes
 shared memory parallel code. The second switch makes each MPI process
 create its own output file (this switch has no effect if MPI is not
 defined).  */

#define MPI
#define PARALLEL_FILES

/* Turn OFF printout (other than error messages) from MPI nodes with
 rank > 0.   This does not affect the model results other than making
 less cluttered standard output.  No effect in non-MPI. */

#ifdef MPI
# define MPI_SILENT_MODE
#endif


/* Switch to allow mixed [tiled + single-block] execution in OpenMP.
 This switch enables special logical branch in "compute_tile_bounds.h"
 to recognize tile=NSUB_X*NSUB_E as covering the whole model grid, and
 increases sizes of scratch arrays declared in "private_scratch" to
 accommodate enough workspace accordingly.  Useful for debugging
 purposes only. Normally should be kept undefined. */

c--#define ALLOW_SINGLE_BLOCK_MODE
#ifdef ALLOW_SINGLE_BLOCK_MODE
# define SINGLE NSUB_X*NSUB_E,NSUB_X*NSUB_E !!!
#endif

/* Select the stage of main time-stepping algorithm at which to run
 the barotropic mode: either during predictor or corrector.  Mutually
 exclusive, only one of the two options must be defined. */

#undef PRED_COUPLED_MODE
#define CORR_COUPLED_MODE

/* In the case of CORR_COUPLED_MODE there are two viable options to
 deal with vertically-integrated advection/Coriolis terms in barotropic
 momentum equations: EITHER (1) compute them at time step "n" and
 recompute their barotropic contributions into "rubar", "rvbar" at
 every fast step; OR compute at "n+1/2" using AB-like extrapolated
 velocities and keep them constant (do not recompute) during
 fast-time stepping. */

#ifdef CORR_COUPLED_MODE
# define EXTRAP_BAR_FLUXES
#endif

/* Switch to recompute Coriolis terms in barotropic momentum equations
 at every barotropic time step.  Applicable to both predictor- and
 corrector-coupled (including the option with forward extrapolation of
 barotropic fluxes) modes. */

c--#define KEEP_CORIOLIS

/* Activation of the first switch of the following two makes computing
 vertical viscous terms as well as implicit part of vertical advection
 take place before the barotropic mode (i.e., moves tri-diagonal solver
 from step3d_uv2 to step3d_uv1) which makes it possible to compute the
 bottom drag term within the tri-diagonal solver, and yet have it
 explicitly available as forcing to the barotropic mode.  The second
 switch makes implicit no-slip b.c. at bottom be an integral part of
 the tri-diagonal solver as opposite to computing it explicitly from
 whatever latest velocity values available. */

c-# define IMPLICIT_BOTTOM_DRAG
# define IMPLCT_NO_SLIP_BTTM_BC


/* Take into account nonuniformity of density field in computation of
 barotropic pressure-gradient terms.  If not activated, then constant
 density Shallow Water Equation SWE-like term is used. No effect for
 a pure 2D configuration. */

#define VAR_RHO_2D


/* Normally, whenever nonlinear EOS is selected, use stiffened and
 split form,
                rho(S,T,z) = rho1(S,T) + qp1(S,T)*[z+qp2*z^2]

 where all pressure dependency is explicit via "z", while "rho1","qp1"
 do not depend on pressure (i.e., "z");  "qp2" is just a constant.
 The above form facilitates computation of adiabatic differences in
 pressure-gradient scheme, prsgrd32AC1.F, which is critical for its
 accuracy and robustness.  Overall this follows Shchepetkin and
 McWilliams (2003, 2011), which in its turn is an adaptation of
 Dukowicz (2001) and earlier work of Sun, Bleck, Rooth, Dukowicz,
 Chassignet, and Killworth (1999) with respect to this matter.   */

#ifdef NONLIN_EOS
# define SPLIT_EOS
#endif




/* Enable mathematically seamless restarting procedure from previously
 saved restart file.   This is done by saving two consecutive steps of
 model prognostic variables along with few additional fields (coupling
 terms, hbls/hbbl, etc).  This applies only to 3D configurations and
 has no effect in pure 2D mode. */

#ifdef SOLVE3D
# define EXACT_RESTART
#endif

/* Assign points within the land to special value rather than zero.
 For RHO-point this means all points where mask_rho == 0.  For U- and
 V-points this applies only to points which are fully inside, that both
 RHO-points adjacent to the respective velocity component in upstream
 and downstream directions are land (if only one land, the other is
 water, then the velocity point in on the coast line, so its value is
 still set to zero according to no-normal flow boundary condition).
*/

#ifdef MASKING
# define MASK_LAND_DATA
#endif


/* Switch ON/OFF double precision for real type variables (since this
 is mostly controlled by mpc and/or compiler options, this CPP-switch
 affects only on the correct choice of netCDF functions, see below)
 and the use QUAD precision for global summation variables, which is
 always desirable, but some compilers do not support it.  */

#define DBLEPREC

/* Define standard dimensions for the model arrays (vertical dimensions
 are inserted explicitly in the code, when needed).  Periodic and
 non-periodic versions may differ by the number of ghost points on each
 edge (2 or 1 respectively). This distinction is present only in the
 purely SHARED MEMORY code. In the case of message passing, when array
 dimensions correspond to a portion of the physical domain (as opposite
 to the whole domain), so two ghost zones are always provided on the
 each side. These data for these two ghost zones is then exchanged by
 message passing.  */

#ifdef MPI
# define GLOBAL_2D_ARRAY -1:Lm+2,-1:Mm+2
# define START_2D_ARRAY -1,-1
# define GLOBAL_1DX_ARRAY -1:Lm+2
# define GLOBAL_1DY_ARRAY -1:Mm+2
!# define START_1D -1
!# define END_1DX Lm+2
!# define END_1DY Mm+2
#else
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define GLOBAL_2D_ARRAY -1:Lm+2,-1:Mm+2
#   define START_2D_ARRAY -1,-1
#  else
#   define GLOBAL_2D_ARRAY -1:Lm+2,0:Mm+1
#   define START_2D_ARRAY -1,0
#  endif
# else
#  ifdef NS_PERIODIC
#   define GLOBAL_2D_ARRAY 0:Lm+1,-1:Mm+2
#   define START_2D_ARRAY 0,-1
#  else
#   define GLOBAL_2D_ARRAY 0:Lm+1,0:Mm+1
#   define START_2D_ARRAY 0,0
#  endif
# endif
#endif

#define PRIVATE_1D_SCRATCH_ARRAY istr-2:iend+2
#define PRIVATE_2D_SCRATCH_ARRAY istr-2:iend+2,jstr-2:jend+2

/* The following macros contain logical expressions which answer
 the question: ''Am I a thread working on subdomain (tile) which is
 adjacent to WESTERN[EASTERN, SOUTHERN, NORTHERN] edge (i.e. physical
 boundary) of the model domain?''  Note that ghost points associated
 with periodicity are NOT considered as physical boundary points by
 these macros. In the case of periodicity and/or MPI-partitioning in
 either direction these macros deliberately left undefined
 (non-existing): if they were, they would always be .false. in this
 case because the corresponding WEST[EAST, ...]_EXCH are .true. and
 periodicity is handled the same way as MPI messages. They are left
 undefined to prevent their mistaken use as flags for proximity to
 the edge of the grid.  */

#ifndef EW_PERIODIC
# ifdef MPI
#  define WESTERN_EDGE istr==iwest .and. .not.west_exchng
#  define EASTERN_EDGE iend==ieast .and. .not.east_exchng
# else
#  define WESTERN_EDGE istr==1
#  define EASTERN_EDGE iend==Lm
# endif
#endif

#ifndef NS_PERIODIC
# ifdef MPI
#  define SOUTHERN_EDGE jstr==jsouth .and. .not.south_exchng
#  define NORTHERN_EDGE jend==jnorth .and. .not.north_exchng
# else
#  define SOUTHERN_EDGE jstr==1
#  define NORTHERN_EDGE jend==Mm
# endif
#endif

/* The following macros are to signal the proximity of computational
 margins -- rows of ghost points points associated with either MPI
 exchange zones and/or periodicity in either direction.   Note that
 macros with _MS_ and logical variables with _ms_ are intended
 STRICTLY to control MPI message passing, while the ones _MS_(_ms_)
 are to signal proximity to EITHER ghost points to be exchanged by
 messages OR due to periodicity WITHOUT distinguishing between these
 two causes.  */

#ifdef MPI
# define WEST_MSG_EXCH  west_msg_exch.and.istr==iwest
# define EAST_MSG_EXCH  east_msg_exch.and.iend==ieast
# define SOUTH_MSG_EXCH south_msg_exch.and.jstr==jsouth
# define NORTH_MSG_EXCH north_msg_exch.and.jend==jnorth

# ifdef EW_PERIODIC
#  define WEST_EXCHNG  istr==iwest
#  define EAST_EXCHNG  iend==ieast
# else
#  define WEST_EXCHNG  west_exchng.and.istr==iwest
#  define EAST_EXCHNG  east_exchng.and.iend==ieast
# endif
# ifdef NS_PERIODIC
#  define SOUTH_EXCHNG jstr==jsouth
#  define NORTH_EXCHNG jend==jnorth
# else
#  define SOUTH_EXCHNG south_exchng.and.jstr==jsouth
#  define NORTH_EXCHNG north_exchng.and.jend==jnorth
# endif
#else
# ifdef EW_PERIODIC
#  define WEST_EXCHNG istr==1
#  define EAST_EXCHNG iend==Lm
# endif
# ifdef NS_PERIODIC
#  define SOUTH_EXCHNG jstr==1
#  define NORTH_EXCHNG jend==Mm
# endif
#endif

/* The following four macros identify position of an MPI-node
 relatively to the edge of the physical grid. They are similar to
 the above, except that: (i) that they apply only to MPI subdomain
 decomposition, hence do not refer to tiling bonds istr,...,jend;
 and, in addition to  that, (ii) their state does not depend on
 periodicity.  */

#ifdef MPI
# define WESTERN_MPI_EDGE iwest+iSW_corn==1
# define EASTERN_MPI_EDGE ieast+iSW_corn==LLm
# define SOUTHERN_MPI_EDGE jsouth+jSW_corn==1
# define NORTHERN_MPI_EDGE jnorth+jSW_corn==MMm
#endif


/* Shorthand CPP macro to flag the need to exchange periodic
 and computational margins. */

#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
# define EXCHANGE
#endif
#if defined EW_PERIODIC || defined MPI
# define EXCH_EAST_WEST
#endif
#if defined NS_PERIODIC || defined MPI
# define EXCH_NORTH_SOUTH
#endif

/* Sometimes an operation needs to be restricted to single MPI process
 (master process). Typically this happens when it is desirable to avoid
 redundant write of the same message by all MPI processes into stdout.
 There are two levels of verbosity avoidance: "mpi_master_only" which
 is always hard coded, and less restrictive "mpi_nonexit_warn" which is
 intended to be user-selectable.  Also occasionally it is needed to
 print MPI-node number into printed message. To do it conditionally
 (MPI code only) add MYID (without preceding comma) into the end of
 the message to be printed.   */

#ifdef MPI
# define mpi_master_only if (mynode==0)
# ifdef MPI_SILENT_MODE
#  define mpi_nonexit_warn if (mynode==0)
#  define MYID !
# else
#  define mpi_nonexit_warn
#  define MYID ,' node =', mynode
# endif
#else
# define mpi_master_only
# define mpi_nonexit_warn
# define MYID !
#endif

/* Similarly, if operation needed to be done by one thread only, e.g.,
 copy a redundantly computed private scalar into shared scalar, or to
 write an error message in situation where it is guaranteed that the
 error condition is discovered redundantly by every thread (and is
 the same for all) and only one needs to complain.  ZEROTH_TILE is
 intended to restrict the operation only to thread which is working
 on south-western tile.

 Occasionally a subroutine designed to process a tile may be called
 to process the whole domain. If it is necessary to distinguish
 whether it is being called for the whole domain (SINGLE_TILE_MODE)
 or a tile.

 All these switches are the same for MPI/nonMPI code.  */

#ifdef MPI
# define ZEROTH_TILE (istr==iwest .and. jstr==jsouth)
# define SINGLE_TILE_MODE (iend-istr==ieast-iwest .and. \
 jend-jstr==jnorth-jsouth)
#else
# define ZEROTH_TILE (istr==1 .and. jstr==1)
# define SINGLE_TILE_MODE (iend-istr==Lm-1 .and.+jend-jstr==Mm-1)
#endif


/* Normally the initial condition exists only as a single time record
 at given time.  This requires the use of a two-time-level scheme
 "forw_start" to start time stepping (in our case a RK2 -- forward
 Euler + trapezoidal correction is used for the initial step).  If
 the run is interrupted and restarted from a single record, the use
 of forward step causes differences between the results obtained by
 a continuous run.  Macro EXACT_RESTART activates the option of
 saving two consecutive time steps into restart file allowing exact
 restart.   */

#ifdef EXACT_RESTART
# define FIRST_TIME_STEP iic==forw_start
#else
# define FIRST_TIME_STEP iic==ntstart
#endif
#ifdef SOLVE3D
# define FIRST_2D_STEP iif==1
#else
# define FIRST_2D_STEP iic==ntstart
#endif

/* Macros to signal that the code is in predictor or corrector stage
 of 3D time stepping, and to signal the stage at which coupling with
 barotropic mode takes place.   The latter typically appears in the
 context of computing 3D --> 2D influences and is done immediately
 before the barotropic stepping.  This macro must be set differently
 for the predictor- and corrector-coupled versions.*/

#define PRED_STAGE nnew==3
#define CORR_STAGE nrhs==3

#ifdef CORR_COUPLED_MODE
# define COUPLED_STAGE_ONLY nrhs==3
#else
# define COUPLED_STAGE_ONLY nnew==3
#endif


/* Computation of global sums uses quad precision real numbers to
 avoid dependency from order of summation by different CPUs due to
 roundoff errors. However, not all compilers support quad precision,
 so it can be conditionally switched on-and-off. */

c---#if defined DBLEPREC && !defined GCC && !defined __IFC \
c---                     && !defined PGI && !defined CRAY
/*-- for now, we keep doing the double precision for global sums
   instead of the quad stuff */
#if defined JEROEN_TEST
# define QUAD 16
# define QuadZero 0.Q0
/* #  define QuadZero 0.0_16 */
#else
# define QUAD 8
# define QuadZero 0.D0
#endif

/* Turn ON/OFF double precision for real type variables, associated
 intrinsic functions and netCDF library functions. It should be noted
 that because ROMS relies on compiler options and "mpc" program (see
 mpc.F) to generate double precision executable from default precision
 source code, this switch actually does NOT affect the size of real
 data and precision of the computation. Its main effect is to select
 the correct netCDF function (nf_xxxx_double/nf_xxxx_float, see below)
 to work properly, so it must be set consistently with mpc settings
 and compiler flags (if any) according to the intended accuracy.
 Additionally, activate the use QUAD precision for global summation
 variables, which is always desirable, but some compilers do not
 support it.        */

#ifdef DBLEPREC
# define float dble
c-# define sqrt dsqrt
c-# define SQRT dsqrt
c-# define exp dexp
c-# define EXP dexp
c-# define dtanh dtanh
c-# define TANH dtanh
#endif

/* Model netCDF input/output control: decide whether to put grid data
 into output files (separate choice for each output file) and select
 appropriate double/single precision types for netCDF input (controlled
 by NF_FTYPE) and netCDF output (NF_FOUT) functions.

 Note: Even if the whole code is compiled with double precision
 accuracy, it is still possible to save history and averages netCDF
 files in single precision in order to save disk space. This happens
 if HIS_DOUBLE switch is undefined. Status of HIS_DOUBLE switch does
 not precision of restart file, which is always kept consistent with
 precision of the code.  */

/* #define HIS_DOUBLE */
#undef PUT_GRID_INTO_RESTART
#define PUT_GRID_INTO_HISTORY
#define PUT_GRID_INTO_AVERAGES

#ifdef DBLEPREC
# define NF_FTYPE nf90_double
# define nf_get_var_FTYPE nf_get_var_double
# define nf_get_var1_FTYPE nf_get_var1_double
# define nf_put_var1_FTYPE nf_put_var1_double
# define nf_get_vara_FTYPE nf_get_vara_double
# define nf_put_vara_FTYPE nf_put_vara_double
# define nf_get_att_FTYPE nf_get_att_double
# define nf_put_att_FTYPE nf_put_att_double
# ifdef HIS_DOUBLE
#  define NF_FOUT nf90_double
#  define nf_put_att_FOUT nf_put_att_double
# else
#  define NF_FOUT nf90_float
#  define nf_put_att_FOUT nf_put_att_real
# endif
#else
# define NF_FTYPE nf90_float
# define nf_get_var_FTYPE nf_get_var_real
# define nf_get_var1_FTYPE nf_get_var1_real
# define nf_put_var1_FTYPE nf_put_var1_real
# define nf_get_vara_FTYPE nf_get_vara_real
# define nf_put_vara_FTYPE nf_put_vara_real
# define nf_get_att_FTYPE nf_get_att_real
# define nf_put_att_FTYPE nf_put_att_real
# define NF_FOUT nf90_float
#endif

/* Specific to IBM XLF operating system. */

#ifdef XLF
# define flush flush_
# define etime etime_
#endif

/* 2021/07: Create single flag for closed boundary example. */

#if !defined OBC_WEST  && !defined OBC_EAST \
 && !defined OBC_SOUTH && !defined OBC_NORTH
# define OBC_NONE
#endif

/* 2021/09: moved from scoord. h */
#define VERT_COORD_TYPE_SM09

/* 2023/05:  Poor man's sea ice model */
#define SEA_ICE_NOFLUX
