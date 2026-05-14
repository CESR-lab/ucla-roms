#include "cppdefs.opt"
#if !defined GRID_LEVEL || GRID_LEVEL == 1

!!    program main               ! Open MP version of ROMS driver

      use param
      use dimensions
#ifdef NHMG
      use nhmg
#endif
      use tracers, only: init_tracers
      use mpi
      use roms_mpi
      use obc_tune


      implicit none
      integer :: ierr
      real*8 tstart, tend

      ierr=1
      if (ierr==1) call MPI_Init(ierr)

      ocean_grid_comm=MPI_COMM_WORLD
      call MPI_Comm_size(ocean_grid_comm, nsize,  ierr)
      call MPI_Comm_rank(ocean_grid_comm, mynode, ierr)

      ierr=0
      call read_inp(ierr)

      if (ierr==0) call mpi_setup(ierr)

      tstart=MPI_Wtime()
      if (ierr==0) then

#ifdef BIOLOGY_BEC2
!     When using BEC, BGC tracer indices are set in init_tracers
!     then used in init_scalars, so this order must be maintained
         call init_tracers
         call init_scalars(ierr)
#else
        call init_scalars(ierr)          ! Initialize global scalars,
                                         ! model tunable parameters,
        call init_tracers
#endif
        if (ierr == 0) then
          call roms_thread
        endif
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
      do iic=ntstart,ntstart+ntimes-1    ! all, exiting is deferred
        diag_sync=.false.                ! until "may_day_flag" is
        call roms_step
        if (diag_sync .and. may_day_flag /= 0) goto 99
      enddo
  99  call stop_timers()                 ! summarized among all nodes
      if (may_day_flag /=0)
     &  call wrt_his_ocean_vars(.true.)  ! Emergency history file output
      call closecdf                      ! and broadcasted by "diag"
      return                             ! [which sets diag_sync=T to
      end                                ! signal uniform state].
#endif

      subroutine roms_init
      use param
#ifdef NHMG
      use nhmg
#endif
      use dimensions
      use scoord
      use grid
      use ocean_vars
      use basic_output
      use scalars
      use analytical
      use obc_tune
      use diagnostics
      use boundary
      use particles
      use cstar_output
      use random_output
      use cdr_frc
# if defined(BIOLOGY_BEC2) || defined(MARBL)
      use bec2_params
      use bgc_io
# endif

      implicit none

      integer ierr

      iic=0  ; kstp=1 ; knew=1  
      iif=1  ; nstp=1          
      nnew=1 ; nrhs=1         

      synchro_flag=.true.    
      diag_sync=.false.     
      priv_count=0         

      call init_arrays

      call set_scoord                  ! Setup vertical stretching
                                       ! functions for S-coordinate
      if (may_day_flag /= 0) goto 99

      call get_grid                    ! Can be analytic or read from netcdf

#ifdef NHMG
      call nhmg_init(nx,ny,nz,npx,npy)
#endif
      if (obc_tuning) call init_orlanski_tune

      if (may_day_flag /= 0) goto 99 !-->  EXIT

      call setup_grid1
      call setup_grid2

      call set_depth
      call lmd_swr_frac                ! may be needed by ana_init.
      call grid_stiffness              ! Report topographic slope
                                       ! parameters rx0, rx1

# if defined(BIOLOGY_BEC2)
!     pH is initialized to 0 in ecosys_init
!     FIXME: if using MARBL with diagnostics then this call is skipped and diagnostics arrays are not initialised to 0
      call ecosys2_init
# endif /*BIOLOGY_BEC2*/

      ! we should have a grid and z_r by now
      if (cdr_source) call init_cdr_frc

#ifdef ANA_INITIAL
                                    ! Set initial conditions for
      call set_forces               ! model prognostic variables,
      call ana_init                 ! may require surface forcing
                                    ! either analytically or read
      if (nrrec > 0) then
#endif


#ifdef EXACT_RESTART
        call get_init(nrrec-1,2)
        call set_depth
#endif
        call get_init(nrrec, 1)

#ifdef ANA_INITIAL
      endif    !<-- nrrec>0
#endif

#ifdef ANA_VMIX
      call ana_vmix
#endif
      if (may_day_flag /= 0) goto 99      !--> ERROR

                                  ! Set initial model clock: at this
      time=start_time             ! moment "start_time" (global scalar)
      tdays=time*sec2day          ! is set by get_init or analytically
                                  ! copy it into threadprivate "time"
      call set_forces
      call set_depth

!----------------------------------------------------------------------
!  Set NHMG horizontal and vertical grids
!  then set matrices coefficients for the elliptic problem
!----------------------------------------------------------------------
#ifdef NHMG
      call nhmg_matrices(nx,ny,N,halo,padd_X,padd_E,dzdxi,dzdeta,Hz,
     &                   dm_r(0:nx+1,0:ny+1),
     &                   dn_r(0:nx+1,0:ny+1)    )
#endif

      call set_HUV

      call omega
      call rho_eos(nrhs)

CR      write(*,*)  ' -2' MYID

      if (floats) then
        call init_particles
      endif

      call set_sponge

      ! The optional argument refers to
      ! blowup if true, and initial if false.
      ! Come up with a better name!
      call wrt_his_ocean_vars(.false.)
#if defined(BIOLOGY_BEC2) || defined(MARBL)
      call wrt_bgc(.true.)
#endif
       call init_diagnostics
      if (do_random) call init_random
      if (do_cstar)  call init_cstar

      if (may_day_flag /=0 ) goto 99     !-->  EXIT

        mpi_master_only write(*,'(/1x,A/)')
     &     'main :: initialization complete, started time-stepping.'

      call diag  ! log file output of global norms

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
      use ocean_vars
      use extract_data
      use diagnostics
      use scalars
      use obc_tune
      use particles
      use surf_flux
      use random_output
      use cstar_output
      use zslice_output
      use basic_output
      use boundary
      use analytical
      use bgc_io
      use tides
      use mixing

      implicit none

      integer :: ierr

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
      ! set_tides on the other hand work wit time, not adding 0.5*dt
      ! use time at n for surface forcing values

      tdays = (time-0.5*dt)*sec2day
      call set_forces
      ! use time at n+1/2 for boundary values
      ! we use these values at the end of pre_step, and during step2d

      tdays = (time-0.0*dt)*sec2day
      call set_bry_all

      if (bry_tides.or.pot_tides) call set_tides

      ! currrently, rho_eos works with t(nrhs), nrhs==nstp => time n
      call rho_eos(nrhs)

      ! Computes horizontal fluxes, using nrhs => time n
      ! Hz is computed by set_depth
      call set_HUV

      ! Computes vertical flux omega, using Flxu,Flxv computed in set_HUV
      ! Here, omega corresponds to time n
      ! Should be already available from previous time-step
      call omega

# if defined LMD_MIXING
      call lmd_vmix(nstp)
# endif

      ! nrhs points to u/v at time level n
      call prsgrd
      call pre_step3d   ! u is m/s here

      !! nnew ==n+1/2


      !! set_HUV1 is like set_HUV, but first removes the barotr/barocl mismatch
      !! set_HUV1 works with nnew, which here is time n+1/2
      ! It's using the Hz from time n though
      ! Look into this to understand what it does
      call set_HUV1

      nrhs=3 ; nnew=3-nstp   !!! WARNING

      !! nnew = n+1
      !! nrhs = n+1/2
      call omega
      call rho_eos(nrhs)  !!! use the right time index

      tdays = (time-0.0*dt)*sec2day
      call set_forces   ! get for time n+1/2
# if defined LMD_MIXING
      call lmd_vmix(nrhs)
# endif

      tdays = (time+0.5*dt)*sec2day
      call set_bry_all  ! get for time n+1
      if (bry_tides.or.pot_tides) call set_tides

      ! All of this seems to happen with Hz's from time n
      ! Corrector step
      call prsgrd
      call step3d_uv1  ! update u/v(nnew) with dt*ru/rv

#  ifdef UV_VIS2
      call visc3d
#  endif



! Solve the 2D equations for the barotropic mode.
!------ --- -- --------- --- --- ---------- -----
      do iif=1,nfast

        kstp=knew
        knew=kstp+1
        if (knew > 4) knew=1

        call step2d

      enddo

      ! step3d_uv2 re-computes FlxU,FlV but still for n+1/2
      call step3d_uv2

      call omega
      call step3d_t
# if defined TS_DIF2 || defined TS_DIF4
      call t3dmix
# endif


      ! at this point, u,v,w and tracers are up to n+1
      ! rho, omega, and fluxes are at n+1/2

      call rho_eos(nnew)
!     call lmd_vmix(nnew)


      time=start_time+dt*float(iic-ntstart) + dt !n+1
      tdays=time*sec2day

      call calc_avg_ocean_vars

!     Done with stepping, time for outputs

      call diag  ! log file output of global norms

      call wrt_his_ocean_vars
      call wrt_rst_ocean_vars
      call wrt_avg_ocean_vars

      if (do_extract) call do_extract_data

      if (obc_tuning)    call adjust_orlanski

      if (floats)     call do_particles

      if (wrt_smflx.or.wrt_stflx) call wrt_sflux

      if (calc_diag) call do_diagnostics

      if (do_random) call wrt_random
#if defined MARBL && defined MARBL_DIAGS
      if (do_cstar)  call wrt_cstar
#endif

      if (do_zslice) call wrt_zslice

#if defined(BIOLOGY_BEC2) || defined(MARBL)
      call wrt_bgc
#endif

      end
