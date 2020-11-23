#include "cppdefs.h"
#if !defined GRID_LEVEL || GRID_LEVEL == 1

!!    program main               ! Open MP version of ROMS driver

#ifdef NHMG
      use nhmg
#endif
      use tracers, only: init_tracers

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
                                         ! model tunable parameters,
        call init_tracers                ! New: initiate tracer variables
        if (ierr == 0) then
C$        call omp_set_dynamic(.false.)
C$OMP PARALLEL                           ! fast-time averaging weights
          call roms_thread               ! for barotropic mode, and
C$OMP END PARALLEL                       ! launch the model in OpenMP
        endif                            ! parallel regime.
#ifdef NHMG  
      ! write nhmg statistics to fort.10
      call nhmg_clean()
#endif
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

#ifdef NHMG
      use nhmg
#endif
      implicit none

      integer trd, tile, my_first, my_last, range
C$    integer omp_get_thread_num, omp_get_num_threads

#include "param.h"
#include "scalars.h"
#include "ncvars.h"

#ifdef NHMG
#include "ocean3d.h"
#include "grid.h"
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

#ifdef NHMG
      call nhmg_init(Lm,Mm,N,NP_XI,NP_ETA)
#endif

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

!----------------------------------------------------------------------
!  Set NHMG horizontal and vertical grids 
!  then set matrices coefficients for the elliptic problem
!----------------------------------------------------------------------
#ifdef NHMG
      call nhmg_matrices(Lm,Mm,N,halo,padd_X,padd_E,dzdxi,dzdeta,Hz,
     &                   dm_r(0:Lm+1,0:Mm+1),
     &                   dn_r(0:Lm+1,0:Mm+1)    )
#endif

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
C$OMP END MASTER
C$OMP BARRIER
CR      write(*,*) '  0' MYID
      if (may_day_flag /=0 ) goto 99     !-->  EXIT

C$OMP MASTER
        mpi_master_only write(*,'(/1x,A/)')
     &     'main :: initialization complete, started time-steping.'
C$OMP END MASTER

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

#ifdef WEC
      use wec_frc, only: wec_forces
#endif
#if defined SSH_TIDES || defined UV_TIDES
      use tides, only: set_tides
#endif

      implicit none
      integer trd, tile, my_first, my_last, range
#include "param.h"
#include "scalars.h"
#include "ncvars.h"


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
        call set_tides(tile)            ! compute surface fluxes.
# endif
        call    rho_eos(tile)
        call    set_HUV(tile)
        call       diag(tile)
# ifdef BIOLOGY
        call   bio_diag(tile)
# endif
      enddo
C$OMP BARRIER

# ifdef WEC
      do tile=my_tile_range
        call  wec_forces(tile)
      enddo
C$OMP BARRIER
# endif

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
        call pre_step3d(tile)  ! u is m/s here
# ifdef AVERAGES
        call    set_avg(tile)
# endif
      enddo
C$OMP BARRIER

      do tile=my_tile_range
        call set_HUV1(tile)
      enddo
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
C$OMP END MASTER
C$OMP BARRIER
        if (iic-ntstart > ntimes) goto 99   !-->  DONE
      endif


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


  99  return
      end