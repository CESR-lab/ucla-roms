/* This is "global_definitions.h": It contains a set of predetermined
  macro definitions which are inserted into the individual files by
  C-preprocessor.  General user is strongly discouraged from attempts
  to modify anything below this line.
------------------------------------------------------------------ */

/* Switch to allow mixed [tiled + single-block] execution in Open
  MP mode. Activation of this switch enables special logical branch
  in "compute_tile_bounds" which recognizes tile=NSUB_X*NSUB_E as
  covering the whole model grid, and it increases sizes of arrays
  declared in "private_scratch" to accomodate enough workspace
  accordingly. This switch is used for debugging purposes only and
  normally should be undefined.  */

c--#define ALLOW_SINGLE_BLOCK_MODE
#ifdef ALLOW_SINGLE_BLOCK_MODE
# define SINGLE NSUB_X*NSUB_E,NSUB_X*NSUB_E !!!
#endif


/* Switch to control at what stage time-stepping of barotropic takes
  place: either during predictor or corrector phase of main time step.
  Note: only one of the two options must be defined. */

#undef PRED_COUPLED_MODE
#define CORR_COUPLED_MODE

/* In the case of CORR_COUPLED_MODE there are two viable options to
  deal with deal with vertically-integrated advection/Coriolis terms
  in barotropic momentum equations: either to compute them at time
  step "n" AND recompute their barotropic contributions at every fast
  step; or compute at "n+1/2" using AB-extrapolated velocities and
  keep them contstant (do not recompute) during fast-time stepping. */

#ifdef CORR_COUPLED_MODE
# define EXTRAP_BAR_FLUXES
#endif


/* Macro to signal that we are in predictor sub-step of 3D mode.
  This does not have alternative settings. */

#define PRE_STEP nnew.eq.3


/* Take into account nonuniform density in barotropic mode pressure-
  gradient terms. If not activated, then Shallow Water Equation (SWE)
  term is used.  This switch has no effect for a 2D problem.   */

#define VAR_RHO_2D


/* Switch ON/OFF double precision for real type variables (since this
  is mostly controlled by mpc and/or compuler options, this CPP-switch
  affects only on the correct choice of netCDF functions, see below)
  and the use QUAD precision for global summation variables, which is
  always desirable, but some compilers do not support it.  */

#define DBLEPREC


/* Turn ON/OFF MPI parallelization. If not activated,the code becomes
  shared memory parallel code. The second switch makes each MPI process
  create its own output file (this switch has no effect if MPI is not
  defined).  */

#define MPI
#define PARALLEL_FILES

/* Define standard dimensions for the model arrays (vertical
  dimensions are inserted explicitly in the code, when needed).
  Periodic and nonperiodic versions may differ by the number of
  ghost points on each edge (2 or 1 respectively). This distinction
  is present only in the purely SHARED MEMORY code. In the case of
  message passing, when array dimensions correspond to a portion of
  the physical domain (as opposite to the whole domain), so two
  ghost zones are always provided on the each side. These data for
  these two ghost zones is then exchanged by message passing.  */

#ifdef MPI
# define GLOBAL_2D_ARRAY -1:Lm+2+padd_X,-1:Mm+2+padd_E
# define START_2D_ARRAY -1,-1
#else
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define GLOBAL_2D_ARRAY -1:Lm+2+padd_X,-1:Mm+2+padd_E
#   define START_2D_ARRAY -1,-1
#  else
#   define GLOBAL_2D_ARRAY -1:Lm+2+padd_X,0:Mm+1+padd_E
#   define START_2D_ARRAY -1,0
#  endif
# else
#  ifdef NS_PERIODIC
#   define GLOBAL_2D_ARRAY 0:Lm+1+padd_X,-1:Mm+2+padd_E
#   define START_2D_ARRAY 0,-1
#  else
#   define GLOBAL_2D_ARRAY 0:Lm+1+padd_X,0:Mm+1+padd_E
#   define START_2D_ARRAY 0,0
#  endif
# endif
#endif

#define PRIVATE_1D_SCRATCH_ARRAY istr-2:iend+2
#define PRIVATE_2D_SCRATCH_ARRAY istr-2:iend+2,jstr-2:jend+2

/* The following definitions contain logical expressions which answer
  the question: ''Am I a thread working on subdomain (tile) which is
  adjacent to WESTERN[EASTERN,SOUTHERN,NORTHERN] edge (i.e. physical
  boundary) of the model domain?'' Note that ghost points associated
  with periodicity are NOT considered as physical boundary points by
  these macros.  In the case of periodicity and MPI-partitioningi in
  either direction these macros are always .false., because WEST[EAST,
  ...]_INTER are .true., and periodicty is handled by MPI messages,
  but they are also should be prevented from being used.  */

#ifdef MPI
# define WESTERN_EDGE istr.eq.iwest .and. .not.west_inter
# define EASTERN_EDGE iend.eq.ieast .and. .not.east_inter
# define SOUTHERN_EDGE jstr.eq.jsouth .and. .not.south_inter
# define NORTHERN_EDGE jend.eq.jnorth .and. .not.north_inter
#else
# define WESTERN_EDGE istr.eq.1
# define EASTERN_EDGE iend.eq.Lm
# define SOUTHERN_EDGE jstr.eq.1
# define NORTHERN_EDGE jend.eq.Mm
#endif

/* Suppress the above if periodicity applies in either direction to
  prevent their mistaken use as flags for proximity to the edge of
  the grid.  */

#ifdef EW_PERIODIC
# undef WESTERN_EDGE
# undef EASTERN_EDGE
#endif
#ifdef NS_PERIODIC
# undef SOUTHERN_EDGE
# undef NORTHERN_EDGE
#endif

/* The following four macros are to signal the proximity of
  communication boundaries: rows of ghost points associated
  with either MPI exchange zones and/or periodicity.  */

#ifdef MPI
# define WEST_INTER west_inter.and.istr.eq.iwest
# define EAST_INTER  east_inter.and.iend.eq.ieast
# define SOUTH_INTER south_inter.and.jstr.eq.jsouth
# define NORTH_INTER north_inter.and.jend.eq.jnorth
#else
# ifdef EW_PERIODIC
#  define WEST_INTER istr.eq.1
#  define EAST_INTER iend.eq.Lm
# endif
# ifdef NS_PERIODIC
#  define SOUTH_INTER jstr.eq.1
#  define NORTH_INTER jend.eq.Mm
# endif
#endif

/* The following four macros identify position of an MPI-node
  relatively to the edge of the physical grid.  They are similar to
  the above, except that: (i) that they apply only to MPI subdomain
  decomposition, hence do not refer to tiling bonds istr,...,jend;
  and, in addition to  that, (ii) their state does not depend on
  periodicity.  */

#ifdef MPI
# define WESTERN_MPI_EDGE iwest+iSW_corn.eq.1
# define EASTERN_MPI_EDGE ieast+iSW_corn.eq.LLm
# define SOUTHERN_MPI_EDGE jsouth+jSW_corn.eq.1
# define NORTHERN_MPI_EDGE jnorth+jSW_corn.eq.MMm
#endif




/* Sometimes an operation needs to be restricted to one MPI process,
  the master process. Typically this occurs when it is desirable to
  avoid redundant write of the same message by all MPI processes into
  stdout.  Also occasionally it is needed to include MPI-node number
  into printed message. To do it conditionally (MPI code only) add
  MYID (without preceeding comma) into the end of the message to be
  printed.   */

#ifdef MPI
# define MPI_master_only if (mynode.eq.0)
# define MYID ,' node =', mynode
#else
# define MPI_master_only
# define MYID !
#endif

/* Similarly, if operation needed to be done by one thread only, e.g.,
  copy a redundantly computed private scalar into shared scalar, or
  write an error message in situation where it is guaranteed that the
  error condition is discovered redundantly by every thread (and is
  the same for all) and only one needs to complain.  ZEROTH_TILE is
  intended to restrict the operation only to thread which is working
  on south-western tile.

  Occasinally a subroutine designed to process a tile may be called
  to process the whole domain. If it is necessary to dustinguish
  whether it is being called for the whole domain (SINGLE_TILE_MODE)
  or a tile.

  All these switches are the same for MPI/nonMPI code.  */

#ifdef MPI
# define ZEROTH_TILE (istr.eq.iwest .and. jstr.eq.jsouth)
# define SINGLE_TILE_MODE (iend-istr.eq.ieast-iwest .and. \
 jend-jstr.eq.jnorth-jsouth)
#else
# define ZEROTH_TILE (istr.eq.1 .and. jstr.eq.1)
# define SINGLE_TILE_MODE (iend-istr.eq.Lm-1 .and.+jend-jstr.eq.Mm-1)
#endif


/* Normally initial condition exists only as a single time record
  at given time.  This requires the use of a two-time-level scheme
  "forw_start" to start time stepping (in our case a RK2 --- forward
  Euler + trapezoidal correction is used for the initial step).  If
  the run is interrupted and restarted from a single record, the use
  of forward step causes differences between the results obtained by
  a continuous run.  Macro EXACT_RESTART activates the option of
  saving two consecutive time steps into restart file allowing exact
  restart.   */

#ifdef EXACT_RESTART
# define FIRST_TIME_STEP iic.eq.forw_start
#else
# define FIRST_TIME_STEP iic.eq.ntstart
#endif
#ifdef SOLVE3D
# define FIRST_2D_STEP iif.eq.1
#else
# define FIRST_2D_STEP iic.eq.ntstart
#endif


/* Turn ON/OFF double precision for real type variables, associated
  intrinsic functions and netCDF library functions. It should be noted
  that because ROMS relies on compiler options and "mpc" program (see
  mpc.F) to generate double precision executable from default
  precision source code, this switch actually does NOT affect the size
  of real data and precision of the computation. Its main effect is to
  select the correct netCDF function (nf_xxxx_double/nf_xxxx_float,
  see below) to work properly, so it must be set consistently with mpc
  settings and compuler flags (if any) according to the intended
  accuracy. Additionally, activate the use QUAD precision for global
  summation variables, which is always desirable, but some compilers
  do not support it.        */

#if defined DBLEPREC && !defined Linux && !defined PGI && !defined __IFC
# define QUAD 16
#  define QuadZero 0.Q0
/* #  define QuadZero 0.0_16 */
#else
# define QUAD 8
# define QuadZero 0.D0
#endif

c-#ifdef DBLEPREC
c-# define float dfloat
c-# define FLoaT dfloat
c-# define FLOAT dfloat
c-# define sqrt dsqrt
c-# define SQRT dsqrt
c-# define exp dexp
c-# define EXP dexp
c-# define dtanh dtanh
c-# define TANH dtanh
c-#endif

/* Model netCDF input/output control: decide whether to put grid data
  into output files (separate choice for each output file) and select
  appropriate double/single precision types for netCDF input
  (controlled by NF_FTYPE) and netCDF output (NF_FOUT) functions.

  NOTE: Even if the whole code is compiled with double precision
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
# define NF_FTYPE nf_double
# define nf_get_var_FTYPE nf_get_var_double
# define nf_get_att_FTYPE nf_get_att_double
# define nf_put_att_FTYPE nf_put_att_double
# define nf_get_var1_FTYPE nf_get_var1_double
# define nf_put_var1_FTYPE nf_put_var1_double
# define nf_get_vara_FTYPE nf_get_vara_double
# define nf_put_vara_FTYPE nf_put_vara_double
# ifdef HIS_DOUBLE
#  define NF_FOUT nf_double
#  define nf_put_att_FOUT nf_put_att_double
# else
#  define NF_FOUT nf_float
#  define nf_put_att_FOUT nf_put_att_real
# endif
#else
# define NF_FTYPE nf_float
# define nf_get_var_FTYPE nf_get_var_real
# define nf_get_att_FTYPE nf_get_att_real
# define nf_put_att_FTYPE nf_put_att_real
# define nf_get_var1_FTYPE nf_get_var1_real
# define nf_put_var1_FTYPE nf_put_var1_real
# define nf_get_vara_FTYPE nf_get_vara_real
# define nf_put_vara_FTYPE nf_put_vara_real
# define NF_FOUT nf_float
#endif

/* The following definitions are machine dependent macros, compiler
  directives, etc. A proper set of definitions is activated by a
  proper choice C-preprocessor flag, i.e. -DSGI for an SGI computer
  or -DCRAY for a Cray shared memory architecture (Y-MP, C-90, J-90).
  Definitions for other shared memory platforms may be appended here.
*/
#if defined SGI || defined sgi
# define CSDISTRIBUTE_RESHAPE !!
/* # define CSDISTRIBUTE_RESHAPE c$sgi distribute         */
/* # define CSDISTRIBUTE_RESHAPE c$sgi distribute_reshape */
# define BLOCK_PATTERN block,block
# define BLOCK_CLAUSE !! onto(2,*)
# define CVECTOR CDIR$ IVDEP
# define CSDOACROSS C$DOACROSS
# define CAND C$&
#elif defined CRAY || defined cray
# define CVECTOR CDIR$ IVDEP
# define CSDOACROSS CMIC$ DO ALL
# define CAND CMIC$&
# define SHARE SHARED
# define LOCAL PRIVATE
#endif




#ifdef AIX
# define flush flush_
# define etime etime_
#endif
