#include "cppdefs.h"
#if !defined GRID_LEVEL || GRID_LEVEL == 1

!!    program main               ! Open MP version of ROMS driver
      implicit none              ! with single parallel region using
      integer ierr               ! explicit barrier synchronization.
# include "param.h"
# ifdef MPI
      real*8 tstart, tend
C$    integer level,req_lev
#  include "mpif.h"
c**   call system('uname -nmr')
      ierr=1
C$    req_lev=MPI_THREAD_MULTIPLE
C$    call MPI_Init_thread(req_lev, level, ierr)
C$ !!   write(*,*) 'MPI thread support levels =', req_lev,level
C$    ierr=0
      if (ierr == 1) call MPI_Init(ierr)

      call mpi_setup(ierr)
      tstart=MPI_Wtime()
      if (ierr == 0) then
# endif
        call init_scalars(ierr)          ! Initialize global scalars,
        if (ierr == 0) then              ! model tunable parameters,
C$        call omp_set_dynamic(.false.)
C$OMP PARALLEL                           ! fast-time averaging weights
          call roms_thread               ! for barotropic mode, and
C$OMP END PARALLEL                       ! launch the model in OpenMP
        endif                            ! parallel regime.
# ifdef MPI
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      tend=MPI_Wtime()
      mpi_master_only write(*,*) 'MPI_run_time =', tend-tstart
      call MPI_Finalize (ierr)
# endif
      stop
      end


      subroutine roms_thread
      implicit none
# include "param.h"
# include "scalars.h"
                                         ! Note: Because there is
      call start_timers()                ! a possibility of I/O error
      call roms_init                     ! occurring on some MPI nodes,
      if (may_day_flag /= 0) goto 99     ! but not simultaneously on
      do iic=ntstart,ntstart+ntimes      ! all, exiting is deferred
        diag_sync=.false.                ! until "may_day_flag" is
        call roms_step
        if (diag_sync .and. may_day_flag /= 0) goto 99
      enddo
  99  call stop_timers()                 ! summarized among all nodes
C$OMP BARRIER                            ! and broadcasted by "diag"
C$OMP MASTER                             ! [which sets diag_sync=T to
      call closecdf                      ! signal that "may_day_flag"
C$OMP END MASTER                         ! is in uniform state]. This
      return                             ! is THE ONLY route to exit
      end                                ! if something goes wrong.
#endif

      subroutine roms_init
      implicit none
      integer trd, tile, my_first, my_last, range
C$    integer omp_get_thread_num, omp_get_num_threads
#include "param.h"
#include "scalars.h"

#include "ncvars.h"
#ifdef FLOATS
! grid.h is needed so that lonr and latr are readily available
# include "grid.h"
# include "floats.h"
# include "ncvars_floats.h"
#endif
      numthreads=1 ; trd=0
C$    numthreads=omp_get_num_threads()
C$    trd=omp_get_thread_num()
      proc(2)=trd

      if (mod(NSUB_X*NSUB_E,numthreads) /= 0) then
C$OMP MASTER
        mpi_master_only write(*,'(/3(1x,A,I3),A/)')
     &    '### ERROR: Wrong choice of numthreads =', numthreads,
     &         'while NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E,'.'
        may_day_flag=8
C$OMP END MASTER
C$OMP BARRIER
        goto 99 !-->  EXIT
      endif
                                       ! NOTE: This code is written
      iic=0  ; kstp=1 ; knew=1         ! under an assumption that all
#ifdef SOLVE3D
      iif=1  ; nstp=1                  ! the scalar variables assigned
      nnew=1 ; nrhs=1                  ! on the left here are placed
#endif
      synchro_flag=.true.              ! into a THREADPRIVATE common
      diag_sync=.false.                ! block so each thread must
      priv_count=0                     ! assign them to same values.

      range=(NSUB_X*NSUB_E+numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
#define my_tile_range my_first,my_last

      do tile=my_tile_range            ! Initialize or "first touch"
        call init_arrays(tile)         ! global arrays (most of them
      enddo                            ! are just set to to zero).
C$OMP BARRIER

c--#define CR
CR      write(*,*) '-11' MYID


#ifdef ANA_GRID
      do tile=my_tile_range            ! Set horizontal curvilinear
        call ana_grid(tile)            ! grid and model bathymetry
      enddo                            ! (either analytically or read
C$OMP BARRIER                          ! from grid netCDF file).
#else
C$OMP MASTER
      call get_grid
C$OMP END MASTER
C$OMP BARRIER
      if (may_day_flag /= 0) goto 99 !-->  EXIT
#endif
      do tile=my_tile_range            ! Compute various metric terms
        call setup_grid1(tile)         ! and their combinations.
      enddo
C$OMP BARRIER
CR      write(*,*) '-10' MYID
      do tile=my_tile_range
        call setup_grid2(tile)
      enddo
C$OMP BARRIER
CR      write(*,*) '-9' MYID

#ifdef SOLVE3D
C$OMP MASTER                           ! Setup vertical stretching
      call set_scoord                  ! functions for S-coordinate
C$OMP END MASTER                       ! system
C$OMP BARRIER
      if (may_day_flag /= 0) goto 99
#endif
CR      write(*,*) ' -8' MYID

#if (defined UV_VIS2 && defined VIS_GRID) ||\
    (defined TS_DIF2 && defined DIF_GRID)
      do tile=my_tile_range            ! Rescale horizontal mixing
        call visc_rescale(tile)        ! coefficients according to
      enddo                            ! local grid size.
C$OMP BARRIER
CR      write(*,*) ' -7' MYID
#endif

#ifdef SOLVE3D
      do tile=my_tile_range            ! Create three-dimensional
        call set_depth(tile)           ! S-coordinate system, which
# ifdef LMD_KPP
        call swr_frac(tile)            ! may be needed by ana_init.
# endif
      enddo
C$OMP BARRIER                          ! Here it is assumed that free
      do tile=my_tile_range            ! surface zeta is at rest state,
        call grid_stiffness(tile)      ! zeta=0). Also find and report
      enddo                            ! extremal values of topographic
C$OMP BARRIER                          ! slope parameters "rx0", "rx1".
CR      write(*,*) ' -6' MYID
#endif

#ifdef ANA_INITIAL
      do tile=my_tile_range            ! Set initial conditions for
        call set_forces(tile)          ! model prognostic variables,
      enddo                            ! either analytically or read
C$OMP BARRIER                          ! from netCDF file.  Note that
      do tile=my_tile_range            ! analytical initialization
        call ana_init(tile)            ! may require surface forcing
      enddo                            ! fields (e.g, wind stress),
C$OMP BARRIER                          ! so call forces first.
      if (nrrec > 0) then              ! Also note that ana_init may
#endif
#ifdef EXACT_RESTART
C$OMP MASTER                           ! setup environmental variables
        call get_init(nrrec-1,2)       ! (e.g. analytical boundary
C$OMP END MASTER                       ! forcing), call it first, even
C$OMP BARRIER                          ! in the case of restart run.
# ifdef SOLVE3D
        do tile=my_tile_range
          call set_depth(tile)
        enddo
C$OMP BARRIER
# endif
#endif
C$OMP MASTER
        call get_init(nrrec, 1)
C$OMP END MASTER
#ifdef ANA_INITIAL
      endif    !<-- nrrec>0
#endif
C$OMP BARRIER
      if (may_day_flag /= 0) goto 99      !--> ERROR
CR      write(*,*) ' -5' MYID
                                  ! Set initial model clock: at this
      time=start_time             ! moment "start_time" (global scalar)
      tdays=time*sec2day          ! is set by get_init or analytically
                                  ! copy it into threadprivate "time"
#ifdef SOLVE3D
      do tile=my_tile_range       ! recompute three-dimensional S-
        call set_depth(tile)      ! coordinate system: at this moment
      enddo                       ! free surface has non-zero status
C$OMP BARRIER

CR      write(*,*)  ' -4' MYID
      do tile=my_tile_range
        call set_HUV(tile)
      enddo
C$OMP BARRIER
CR      write(*,*)  ' -3' MYID

      do tile=my_tile_range
        call omega(tile)
        call rho_eos(tile)
      enddo
C$OMP BARRIER
CR      write(*,*)  ' -2' MYID
#endif

! Set up climatological environment: Set nudging coefficient for
!==== == ============== ============ sea-surface hight and tracer
! climatology; create analytical tracer and sea-surface hight
! climatology fields (if applicable); set bottom sediment grain
! size [m] and density [kg/m^3] used for bottom boundary layer
! formulation;

#if defined SPONGE || defined TCLIMATOLOGY \
  || (defined SG_BBL96 && defined ANA_BSEDIM)\
  || (defined TCLIMATOLOGY && defined ANA_TCLIMA)\
  || defined ANA_SSH

      do tile=my_tile_range
# if defined SPONGE || defined TCLIMATOLOGY
        call set_nudgcof(tile)
# endif
# if defined TCLIMATOLOGY && defined ANA_TCLIMA && defined SOLVE3D
        call ana_tclima(tile)
# endif
# ifdef ANA_SSH
        call ana_ssh(tile)
# endif
# if defined SG_BBL96 && defined ANA_BSEDIM
        call ana_bsedim(tile)
# endif
      enddo
C$OMP BARRIER
#endif
CR      write(*,*) ' -1' MYID

! Read initial input data for forcing fields; tracer and sea surface
! climatology; bottom sediment grain size and density (if applicable)
! from input netCDF files.  Note that CPP-logic here is mutually
! exclusive with respect to calls ana_tclima, ana_ssh, and ana_bsedim
! just above.

C$OMP MASTER
#ifdef ANA_GRID
        call wrt_ana_grid
#endif
        if (ldefhis .and. wrthis(indxTime)) call wrt_his
#ifdef FLOATS
! Initialization for Lagrangian floats
!-------------------------------------------------------
      nrecflt=0      ! initialization done here and not in
      ncidflt=-1     ! init_scalars since it must be done only
                     ! once (whether child levels exist or not)
      spval=1.E15    ! spval is the nodata flag for float variables

      deltac2p=2.3   ! distance from the boundary at which a float
                     ! is transferred from child to parent
      deltap2c=2.5   ! same for transfer from parent to child

      call init_arrays_floats
      call init_floats
# ifdef SPHERICAL
      call interp_r2d_type_ini (lonr(START_2D_ARRAY), iflon)

      call interp_r2d_type_ini (latr(START_2D_ARRAY), iflat)
# else
      call interp_r2d_type_ini (  xr(START_2D_ARRAY), iflon)
      call interp_r2d_type_ini (  yr(START_2D_ARRAY), iflat)
# endif
# ifdef SOLVE3D
      call fill_ini   ! fills in trackaux for ixgrd,iygrd,izgrd
                      ! and ifld (either izgrd or ifld is meaningful)
# endif
      if (ldefflt) call wrt_floats
#endif /* FLOATS */
C$OMP END MASTER
C$OMP BARRIER
CR      write(*,*) '  0' MYID
      if (may_day_flag /=0 ) goto 99     !-->  EXIT

C$OMP MASTER
        mpi_master_only write(*,'(/1x,A/)')
     &     'main :: initialization complete, started time-steping.'
C$OMP END MASTER

! Initialize child-grid model, if any. Note that this call introduces
! recursive sequence, since child-grid model may also have child on its
! own.

#ifdef GRID_LEVEL
# if GRID_LEVEL < MAX_GRID_LEVEL
#  if GRID_LEVEL == 1
        call roms_init_2
#  elif GRID_LEVEL == 2
        call roms_init_3
#  elif GRID_LEVEL == 3
        call roms_init_4
#  endif
# endif
#endif
  99  return
      end


!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     ***
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***


      subroutine roms_step
      implicit none
      integer trd, tile, my_first, my_last, range
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#ifdef FLOATS
# include "ncvars_floats.h"
# include "floats.h"
      integer chunk_size_flt, Lstr,Lend, flt_str
      common /floats_step/ flt_str
#endif
#ifdef GRID_LEVEL
      integer iter
#endif

      trd=proc(2)
      range=(NSUB_X*NSUB_E+numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

! Increment time-step index and set model clock. Note that "time" set
! below corresponds to step "n" (denoted here as "nstp"), while counter
! "iic" corresponds to "n+1", so normally, assuming that time is
! counted from zero, the following relation holds: time=dt*(iic-1).
!  Also note that the output history/restart/averages routines write
! time and all the fields at step "n" (not n+1), while the first
! element of structure "time_index" written into the files is actually
! iic-1, hence normally time=time_index*dt there.  Same rule applies
! to the diagnostic routine "diag" which prints time and time step
! (actually iic-1) on the screen.

      time=start_time+dt*float(iic-ntstart) !<-- corresp. to "nstp"
      tdays=time*sec2day
#ifdef SOLVE3D
      nstp=1+mod(iic-ntstart,2)
      nrhs=nstp ; nnew=3
#endif

#ifdef FLOATS
      nfp1=mod(nfp1+1,NFT+1)           ! rotate time
      nf  =mod(nf  +1,NFT+1)           ! indices for
      nfm1=mod(nfm1+1,NFT+1)           ! floats
      nfm2=mod(nfm2+1,NFT+1)
      nfm3=mod(nfm3+1,NFT+1)
C$OMP MASTER
      flt_str=0
C$OMP END MASTER
#endif
                                       ! Read forcing and climatology
      if (synchro_flag) then           ! data. This operation may raise
        synchro_flag=.false.           ! "may_day_flag" in the case of
C$OMP MASTER                           ! I/O errors which may occur on
        call get_forces                ! some nodes, but not on all at 
C$OMP END MASTER                       ! the same time.  Therefore to
C$OMP BARRIER                          ! avoid MPI deadlock the exit is
      endif                            ! deferred until after broadcast
                                       ! of "may_day_flag" by diag.
#ifdef SOLVE3D
      do tile=my_tile_range                 ! interpolate forcing
        call set_forces(tile)               ! data to model time and
# if defined SSH_TIDES || defined UV_TIDES
        call set_tides(tile)                ! compute surface fluxes.
# endif
        call    rho_eos(tile)
        call    set_HUV(tile)
        call       diag(tile)
# ifdef BIOLOGY
        call   bio_diag(tile)
# endif
      enddo
C$OMP BARRIER

      do tile=my_tile_range
        call omega(tile)
# if defined ANA_VMIX
        call ana_vmix(tile)
# elif defined LMD_MIXING
        call lmd_vmix(tile)
# elif defined BVF_MIXING
        call bvf_mix(tile)
# endif
      enddo
C$OMP BARRIER

      do tile=my_tile_range
        call     prsgrd(tile)
        call pre_step3d(tile)
# ifdef PRED_COUPLED_MODE
#  ifdef UV_VIS2
        call     visc3d(tile)
#  endif
# endif
# ifdef AVERAGES
        call    set_avg(tile)
# endif
      enddo
C$OMP BARRIER

# ifdef CORR_COUPLED_MODE
      do tile=my_tile_range              ! This code segment is for
        call set_HUV1(tile)              ! corrector-coupled version
      enddo                              ! of time stepping only.
C$OMP BARRIER

      nrhs=3 ; nnew=3-nstp   !!! WARNING

      do tile=my_tile_range
        call omega(tile)
        call rho_eos(tile)
      enddo
C$OMP BARRIER

      do tile=my_tile_range
        call     prsgrd(tile)
        call step3d_uv1(tile)
#  ifdef UV_VIS2
        call     visc3d(tile)
#  endif
      enddo
C$OMP BARRIER
# endif
#endif  /* SOLVE3D */

! Output block: write restart/history files.
!======= ====== ===== =============== ======

      if ( iic > ntstart .and. ( mod(iic-ntstart,nrst) == 0
#ifdef EXACT_RESTART
     &                         .or. mod(iic-ntstart+1,nrst) == 0
#endif
     &   .or. (mod(iic-ntstart,nwrt) == 0 .and. wrthis(indxTime))
#ifdef AVERAGES
     &   .or. (mod(iic-ntsavg,navg) == 0  .and. wrtavg(indxTime))
#endif
#ifdef STATIONS
     &   .or. (mod(iic-ntstart,nsta) == 0 .and. nstation > 0)
#endif
#ifdef FLOATS
     &   .or. (mod(iic-ntstart,nflt) == 0 .and. nfloats > 0)
#endif
     &                                                  )) then
C$OMP MASTER
        if (mod(iic-ntstart,nrst) == 0
#ifdef EXACT_RESTART
     &         .or. mod(iic-ntstart+1,nrst) == 0
#endif
     &                                ) nrecrst=nrecrst+1
        if (mod(iic-ntstart,nwrt) == 0) nrechis=nrechis+1
#ifdef AVERAGES
        if (mod(iic-ntsavg,navg) == 0) nrecavg=nrecavg+1
#endif
#ifdef STATIONS
        if (mod(iic-ntstart,nsta) == 0) nrecstn=nrecstn+1
#endif
#ifdef FLOATS
        if (mod(iic-ntstart,nflt) == 0) nrecflt=nrecflt+1
#endif
        if (mod(iic-ntstart,nrst) == 0
#ifdef EXACT_RESTART
     &         .or. mod(iic-ntstart+1,nrst) == 0
#endif
     &                                 ) call wrt_rst
        if (mod(iic-ntstart,nwrt) == 0 .and. wrthis(indxTime)) then
          call wrt_his
c          if (iic>60) nwrt=1 !<-- useful for debugging
        endif

#ifdef AVERAGES
        if (mod(iic-ntsavg,navg) ==0 .and. wrtavg(indxTime))
     &      call wrt_avg
#endif
#ifdef STATIONS
        if (mod(iic-ntstart,nsta)==0 .and. nstation>0)
     &      call wrt_statn
#endif
#ifdef FLOATS
        if (mod(iic-ntstart,nflt)==0 .and. nfloats>0)
     &      call wrt_floats
        diagfloats=.false.
#endif
C$OMP END MASTER
C$OMP BARRIER
        if (iic-ntstart > ntimes) goto 99   !-->  DONE
      endif

#ifdef FLOATS
! flag for diagnostic computation (for writing at next time step)
      if (mod(iic-ntstart,nflt)==0) then
        diagfloats=.true.
      endif
#endif

! Solve the 2D equations for the barotropic mode.
!------ --- -- --------- --- --- ---------- -----
#ifdef SOLVE3D
      do iif=1,nfast
#endif
#define FORW_BAK
#ifdef FORW_BAK
        kstp=knew                      ! This might look a bit silly,
        knew=kstp+1                    ! because both branches of this
        if (knew > 4) knew=1           ! "if" statement are identical.
        if (mod(knew,2) == 0) then     ! Nevertheless, it makes sense,
          do tile=my_tile_range        ! since mpc will reverse one of
# ifndef SOLVE3D
            call set_forces(tile)      ! these loops to make zig-zag
# endif
            call     step2d(tile)      ! tile-processing sequence.
          enddo
C$OMP BARRIER
        else
          do tile=my_tile_range
# ifndef SOLVE3D
            call set_forces(tile)
# endif
            call     step2d(tile)
          enddo
C$OMP BARRIER
        endif
#else
        kstp=knew
        knew=3
        do tile=my_tile_range
# ifndef SOLVE3D
          call set_forces(tile)
# endif
          call     step2d(tile)
        enddo
C$OMP BARRIER
        knew=3-kstp
        do tile=my_tile_range
           call step2d(tile)
        enddo
C$OMP BARRIER
#endif
#ifdef SOLVE3D
      enddo    ! <-- iif

# ifdef PRED_COUPLED_MODE
      do tile=my_tile_range              ! This code segment is for
        call set_HUV1(tile)              ! predictor-coupled version
      enddo                              ! of time stepping only. 
C$OMP BARRIER

      nrhs=3 ; nnew=3-nstp

      do tile=my_tile_range
        call omega(tile)
        call rho_eos(tile)
      enddo
C$OMP BARRIER
      do tile=my_tile_range
        call     prsgrd(tile)
        call step3d_uv1(tile)
      enddo
C$OMP BARRIER
# endif

      do tile=my_tile_range              ! Continue solution of
        call step3d_uv2(tile)            ! the three-dimensional
      enddo                              ! equations: finalize
C$OMP BARRIER                            ! time step for momenta
      do tile=my_tile_range              ! and tracers.
        call omega(tile)
        call step3d_t (tile)
# if defined TS_DIF2 || defined TS_DIF4
        call t3dmix(tile)
# endif
      enddo
C$OMP BARRIER
#endif /* SOLVE3D */

#ifdef FLOATS
      chunk_size_flt=32
      do while (flt_str < nfloats)
C$OMP CRITICAL
        Lstr=flt_str+1
        flt_str=Lstr+chunk_size_flt-1
C$OMP END CRITICAL
        Lend=min(Lstr+chunk_size_flt-1,nfloats)
        call step_floats(Lstr,Lend)
      enddo
c**    call step_floats(1,nfloats)    ! serial version for debugging
#endif

#ifdef GRID_LEVEL
# if GRID_LEVEL < MAX_GRID_LEVEL
      do iter=1,3
#  if GRID_LEVEL == 1
        call roms_step_2
#  elif GRID_LEVEL == 2
        call roms_step_3
#  elif GRID_LEVEL == 3
        call roms_step_4
#  endif
      enddo
# endif
#endif
  99  return
      end
