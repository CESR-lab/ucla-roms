#include "cppdefs.opt"
#if !defined GRID_LEVEL || GRID_LEVEL == 1

!!    program main               ! Open MP version of ROMS driver

      use param
#ifdef NHMG
      use nhmg
#endif
      use tracers, only: init_tracers
      use mpi
      use sponge_tune


      implicit none                                        ! with single parallel region using
      integer :: ierr                                      ! explicit barrier synchronization.
# ifdef MPI
      real*8 tstart, tend
C$    integer level,req_lev

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
        ! Devin added this
        call init_tracers                ! New: initiate tracer variables (must be before init_scalars)
        call init_scalars(ierr)          ! Initialize global scalars,
                                         ! model tunable parameters,

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
      use param
      use scalars
      use ocean_vars
      use basic_output

      implicit none

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
      if (may_day_flag /=0) 
     &  call wrt_his_ocean_vars(.true.)  ! Emergency history file output
C$OMP BARRIER                            ! and broadcasted by "diag"
C$OMP MASTER                             ! [which sets diag_sync=T to
      call closecdf                      ! signal that "may_day_flag"
C$OMP END MASTER                         ! is in uniform state]. This
      return                             ! is THE ONLY route to exit
      end                                ! if something goes wrong.
#endif

      subroutine roms_init
      use param
#ifdef NHMG
      use nhmg
#endif
      use dimensions
      use grid
      use ocean_vars
      use basic_output
      use scalars
      use analytical
      use sponge_tune
      use diagnostics
      use boundary

      implicit none

      integer trd, tile, my_first, my_last, range
C$    integer omp_get_thread_num, omp_get_num_threads
      integer ierr

      numthreads=1 ; trd=0
C$    numthreads=omp_get_num_threads()
C$    trd=omp_get_thread_num()
      proc(2)=trd

      if (mod(NSUB_X*NSUB_E,numthreads) /= 0) then
C$OMP MASTER
        mpi_master_only write(*,'(/3(1x,A,I3),A/)')
     &    '### ERROR: Wrong choice of numthreads =', numthreads,
c        endif
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

      do tile=my_tile_range              ! Initialize or "first touch"
        call init_arrays(tile,my_first)  ! global arrays (most of them
      enddo                              ! are just set to to zero).
C$OMP BARRIER

#ifdef DIAGNOSTICS
      call init_diagnostics  ! call before nnew and nstp have swapped indices (as if previous timestep)
#endif

#ifdef NHMG
      call nhmg_init(Lm,Mm,N,NP_XI,NP_ETA)
#endif
#ifdef SPONGE_TUNE
      if (ub_tune) call init_orlanski_tune
#endif

c--#define CR
CR      write(*,*) '-11' MYID


      !! move the ana_grid call to within get_grid
#ifdef ANA_GRID
      do tile=my_tile_range       
        call ana_grid(tile)              ! grid and model bathymetry
      enddo                   
#else
      call get_grid
#endif
      if (may_day_flag /= 0) goto 99 !-->  EXIT

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

# ifdef BIOLOGY_BEC2
      ! pH is initialized to 0 in ecosys_init
      call ecosys2_init
# endif

#ifdef ANA_INITIAL
                                    ! Set initial conditions for
      call set_forces               ! model prognostic variables,
      call ana_init(0)              ! may require surface forcing
                                    ! either analytically or read
      if (nrrec > 0) then 
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
        call set_HUV
      enddo
C$OMP BARRIER
CR      write(*,*)  ' -3' MYID

      call omega
      call rho_eos(nrhs)

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

!     call set_global_ij_SW_NE_corners
      call wrt_his_ocean_vars


C$OMP MASTER
#ifdef ANA_GRID
      call wrt_ana_grid
#endif
!      if (ldefhis .and. wrthis(indxTime)) call wrt_his
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

      use param
      use wec_frc
      use tides
      use extract_data
      use diagnostics
      use scalars
      use analytical
      use sponge_tune
      use particles
      use random_output
      use ocean_vars
      use basic_output
      use dimensions
      use boundary

      implicit none
!     integer trd, tile, my_first, my_last, range

!     trd=proc(2)
!     range=(NSUB_X*NSUB_E+numthreads-1)/numthreads
!     my_first=trd*range
!     my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      !! Examine all outputs and try to move these to the end of the timestep, while collecting 
      !! Diagnostic terms along the way during the corrector step.
      !! 
#ifdef DIAGNOSTICS
      call do_diagnostics  ! call before nnew and nstp have swapped indices (as if previous timestep)
#endif
      call wrt_random

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

      nstp=1+mod(iic-ntstart,2)
      nrhs=nstp ; nnew=3

      ! Interpolate forcing data to model time n+1/2
      ! get_frc (both versions work with tdays, adding 0.5*dt)
      call set_forces
      call set_bry_all

      ! currrently, rho_eos works with t(nrhs), nrhs==nstp => time n
      call rho_eos(nrhs)

      ! Computes horizontal fluxes, using nrhs => time n
      ! Hz is computed by set_depth
      call set_HUV

      call diag

      ! Computes vertical flux, using Flxu,Flxv computed in set_HUV
      ! Here, omega corresponds to time n
      ! Should be already available from previous time-step
      call omega

# if defined ANA_VMIX
      call ana_vmix(0)
# elif defined LMD_MIXING
      call lmd_vmix(nstp)
# elif defined BVF_MIXING
      call bvf_mix(0)
# endif

      ! nrhs points to u/v at time level n
      call prsgrd
      call pre_step3d(0)  ! u is m/s here
      !! nnew ==n+1/2  


      call calc_avg_ocean_vars ! call here not wrt_his else
                               ! rho, omega and w will change
      !! set_HUV1 is like set_HUV, but first removes the barotr/barocl mismatch
      !! set_HUV1 works with nnew, which here is time n+1/2
      ! It's using the Hz from time n though
      ! Look into this to understand what it does
      call set_HUV1(0)

      nrhs=3 ; nnew=3-nstp   !!! WARNING

      !! nnew = n+1
      !! nrhs = n+1/2
      call omega
      call rho_eos(nrhs)  !!! use the right time index

      ! Corrector step
      call prsgrd
      call step3d_uv1(0)  ! update u/v(nnew) with dt*ru/rv

#  ifdef UV_VIS2
!     if (split_bc) then
!        call visc3d(u_slow,v_slow,visc_slow)
!        call visc3d(u_fast,v_fast,visc_fast)
!     else
        !call visc3d(u,v,visc3d)
        call visc3d(0)
!     endif
#  endif

! Output block: write restart/history files.
!======= ====== ===== =============== ======

#ifdef PARTICLES
      call do_particles
#endif
#ifdef EXTRACT_DATA
      call do_extract_data
#endif

      ! wrt_***_ocean_vars are called here, in the middle of a roms_step, because
      ! 'rho', 'Akv', etc have only just been computed for step 'n'.
      ! whereas, u & ubar, etc, remain unchanged from the end of the previous step
!     call wrt_his_ocean_vars
      call wrt_avg_ocean_vars
      call wrt_rst_ocean_vars

! Solve the 2D equations for the barotropic mode.
!------ --- -- --------- --- --- ---------- -----
      do iif=1,nfast

        kstp=knew
        knew=kstp+1
        if (knew > 4) knew=1

        call step2d(0) 

      enddo    ! <-- iif

      ! step3d_uv2 re-computes FlxU,FlV but still for n+1/2
      call step3d_uv2(0)   

      call omega
      call step3d_t(0)
# if defined TS_DIF2 || defined TS_DIF4
!       if (split_bc) then
!         call t3dmix(t_slow,diff2_slow)
!         call t3dmix(t_fast,diff2_fast)
!       else
!         call t3dmix(t,diff2)
      call t3dmix
!       endif
# endif


      ! at this point, u,v,w and tracers are up to n+1
      ! rho, omega, and fluxes are at n+1/2

      call rho_eos(nnew)
!     call lmd_vmix(nnew)

      time=start_time+dt*float(iic) ! n+1
      tdays=time*sec2day
      call wrt_his_ocean_vars


#ifdef SPONGE_TUNE
      if (ub_tune) call adjust_orlanski
#endif


  99  return
      end
