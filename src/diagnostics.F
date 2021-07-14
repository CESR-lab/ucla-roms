      module diagnostics

      ! To compute diagnostics use DIAGNOSTICS flag in cppdefs.h.
      ! Also, set configuration at 'User input' section below.

      ! JM: We are not properly accounting for the Udzeta/dt,Vdzeta/dt terms
      ! These should only be a signicant error for coastal tidal flows.

      ! INFO (Optional): ![
      ! Online diagnostics for ROMS simulation.
      ! Initial coding DevinD and JeroenM - 2021/05

      ! HIS vs AVG:
      ! To chose between diagnostics at history timesteps or averages,
      ! set diag_avg to .false. for history, or .true. for averaging.
      ! Averaging is much less efficient as diagnostics must be collected
      ! at everytime step to maintain a true average.

      ! MAIN SUBROUTINE:
      ! The main subroutine is wrt_diagnostics_tile, which is called at the very
      ! beginning of main.F in order to record the diagnostics of 'nnew' index,
      ! which was t=n+1 for the previous timestep, and is now t=n since roms_step
      ! has moved onto the next timestep. I.e. 'nnew' will become 'nstp' in next timestep.
      ! But wrt_diagnostics_tile is called before nnew and nstp indicies change in main.F.

      ! If you are needing averages then the diagnotics need to be computed at everytimestep.
      ! If you just need snapshots then the diagnostics only need to be calculated at
      ! your output frequency.
      ! Both of these are controlled by the logical 'calc_diag', which is set in
      ! wrt_diagnostics_tile at the beginning of a timestep, and then that logical is
      ! used in the equation system to trigger a calculation of diagnostics.

      ! TRACERS:
      ! Since there are 5 diagnostic terms for each tracer, it can become very memory
      ! intensive to keep diag arrays for every tracer. Instead, the user selects
      ! only tracers of interest by setting T/F in wrt_t_diags in tracers_defs.h.
      ! Then the number of selected tracers is counted in init_diags_trc and an array
      ! of the correct size is allocated there. Some book keeping variables are needed
      ! to convert indices between t(itrc) the tracer variable in roms for all tracers,
      ! and Tdiag(td) in diagnostics which only has indices for selected tracers.

      ! UNITS:
      ! We need the rate of change of u or v or tracer, so the diagnostic is
      ! divided by 'dt' (the timestep) to get a diagnostics rate per second.
      ! The roms history file give u/v in m/s, and 't' tracers in their
      ! concentrations. However, since the depth of the water column is also evolving
      ! during a timestep, it is easier to work in dz*u/v/w/t / dt, since the unit
      ! of u/v/w/t is e.g. m^2/s (u*dz) when the depth is changed in the timestep as follows:
      ! step3d_uv1.F (depth not changed just to show position amongst code)
      ! step2d_tile -> set_depth_tile -> (depth (Hz) changes to t=n+1)
      ! step3d_uv2.F (depth not changed just to show position amongst code)
      !
      ! See Documentation/readme-specifics/README.diagnostics_u_v_value_changes
      ! for more details on how the units of 'u' change during a timestep.

      ! U/V evolution:
      ! During a roms_step in main.F, the values of u/v at t=n+1/2 is calculated
      ! in pre_step3d4s.F. Since we are interested in the change from t=n to t=n+1,
      ! we do not track any diagnostics in pre_step3d.
      ! The order of evolution of u/v (showns as u) and its units is as follows:
      ! 1) dz*u - Hydrostatic pressure gradient - prsgrd.F
      ! 2) dz*u - Coriolis & curvilinear coords - compute_horiz_rhs_uv_term.h (2nd time round i.e. not during pre_step3d)
      ! 3) dz*u - Advection                     - step3d_uv1.F & compute_horiz/vert_rhs_uv_term.h
      ! 4) dz*u - Viscous terms (horiz. mixing) - visc3d_GP.F
      !           only does sponge layer (zero in interior beyond sponge)
      ! -         Now depth (Hz) at t=n+1 in calculated in step2d_FB -> set_depth_tile.
      ! 5)    u - Vertical mix & 2D/3D coupling - step3d_uv2.F part 1
      ! 6)    u - Boundary changes              - u3dbc_im.F
      ! 7)    u - 2D-3D correction              - step3d_uv2.F part 2

      ! Clearly, care needs to be taken when comparing 'u' at each diagnostic step,
      ! since the units after a step can change.

      ! NOTE - steps 5)-7) stored in same array:
      ! Since the change at the boundary is (just that) only at the boundary, it
      ! would be wasteful to store it as its own 3D diagnostic term.
      ! Since all of the other changes are only over the interior (ignoring the boundary),
      ! we chose to use the boundary grid points of vertical mixing to store the boudnary
      ! changes. So the vertical mixing diagnostic term is actually the boundary diagnostic
      ! along its boundaries.
      ! Similarly, the 2D-3D correction is actually only a 2D term (applied consistently
      ! vertically), and it is really just a numerically adjustment (non-physical), we
      ! lump it in with vertical mixing, too. This again saves memory, and it's value
      ! is very small so vertical mixing is the dominant signal in that shared array.
      ! The vertical mixing term [5)] was chosen because, as with steps 6) and 7), it
      ! is calculated after Hz is updated. Therefore, it is more likely to be consistent
      ! than is steps 6) & 7) were addded to arrays of steps 1)-4) rather.

      ! TRACER evolution:
      ! 1) dz*C - h_adv - Horizontal advective fluxes                               - compute_horiz_tracer_fluxes.h
      ! 2) dz*C - v_adv - Vertical advective fluxes                                 - compute_vert_tracer_fluxes.h
      ! -         Surface flux is done here
      ! 3) dz*C - v_mix - Nonlocal transport flux for unstable (convective) forcing - step3d_t_ISO.F (lmd_kpp)
      !                   currently only temp & salt. Chose to share v_mix array
      ! 4)    C - v_mix - Implicit time step for vertical diffusion                 - step3d_t_ISO.F
      ! 5)    C - bc    - Boundary changes                                          - t3dbc_im.F
      ! 6)    C - h_mix - Horizontal Laplacian diffusion terms                      - t3dmix_GP.F

      ! As with u/v, the boundary change of tracers (step 5) is placed in step 2) for
      ! memory efficiency.

      ! The budget closing on tracers looks unusual because the error in the sponge
      ! layer is larger than in the interior beyond the sponge layer. E.g. error
      ! in sponge layer might be e-13, but beyond the sponge layer it might be
      ! between exactly 0 and e-17.

      ! W evolution:
      ! 1)    ? -                                       - prsgrd.F
      ! 2)    ? -                                       - compute_horiz_rhs_w_terms.h
      ! 3)    ? -                                       - compute_vert_rhs_w_terms.h'
      ! 4)    ? -                                       - step3d_uv2.F
      ! 5)    ? -                                       - w3dbc_im.F

      ! Currently 'w' is untested and unverified and is likely to contain errors.

      ! LOGICAL PARAMETERS:
      ! It is assumed that by using logicals that are set as known parameters at compile
      ! time, the 'if branch' will be removed during compiling, since the outcome
      ! of the branch is always known and unchanging. If this is true it will allow
      ! for if statements to be placed within loops, which is usually not advised
      ! due to the cost of the if statment at every index iteration. This has still
      ! mostly been avoided anyway.

      ! INTERIOR VS BOUNDARY GRID POINTS:
      ! Most diagnostic terms only act on the interior of the domain and exclude the
      ! boundary. The boundary term obviously only updates the boundary, whereas the
      ! 2D/3D correction updates the interior and boundary. This means that indices
      ! for loops and array statements are important. Here is a quick refresher of loop
      ! ranges (u has one less point in xi, v has one less point in eta):
      !
      !               |       WESTERN_EDGE            |        EASTERN_EDGE          |
      !     variable  | west-bry     | interior-start |  interior-end |    east-bry  |
      !  -------------|--------------|----------------|---------------|--------------|
      !        u      | istr         | istr+1=istrU   | iend          | iend+1=iendR |
      !        v      | istr-1=istrR | istr           | iend          | iend+1=iendR |
      !       rho     | istr-1=istrR | istr           | iend          | iend+1=iendR |
      !
      !
      !     variable  |      u       |      v         |     rho       |
      !  -------------------------------------------------------------|
      !  north-bry    | jend+1=jendR | jend+1=jendR   | jend+1=jendR  | <- NORTHERN_EDGE
      !  interior-end | jend         | jend           | jend          |_______________
      !  interior-str | jstr         | jstr+1=jstrV   | jstr          |
      !  south-bry    | jstr-1=jstrR | jstr           | jstr-1=jstrR  | <- SOUTHERN_EDGE

      !]

#include "cppdefs.h"
#ifdef DIAGNOSTICS

      use dimensions
      use read_write
      use netcdf
      use tracers

      implicit none
      private

#include "ocean2d.h"
#include "ocean3d.h"
! scalars.h for timestep indices 'nnew', 'knew', etc
#include "scalars.h"

      ! ***************************************************************
      ! User inputs:
      !   Budget closing check:
      !     when using for new simulation (especially if different cppflags) make sure
      !     to do a very short simulation (i.e. 10 steps) with check_diag (and check_avg
      !     if using averages) as .true. to ensure budget is closing correctly, e.g.
      !     u_budget -> 0! Then set them as false for your production run as they are a
      !     waste of resources.
      !   Notes:
      !     1) need DIAGNOSTICS flag in cppdefs.h
      !     2) units are du/dt*dz - m^2/s^2 (or v/w/tracer) - i.e. vertically integerated over cell
      !     3) if tracer diags required (diag_trc==T), select which tracers you want in wrt_t_diag array in tracers_defs.h
      !     4) for an example try Examples/Diagnostics/
      !     5) averaging frequency = rec_rate * time_step
      logical, parameter         :: diag_avg      = .true.  ! compute history (F) or averages (T)
      integer                    :: rec_rate      = 2       ! output record rate, move this to .in file rather?
      integer                    :: recs_per_file = 20      ! total recs per file

      logical, parameter, public :: diag_uv       = .true.  ! calc+write diagnostics of u/v.
      logical, parameter, public :: diag_trc      = .true.  ! calc+write diagnostics of selected tracers. Select tracers in wrt_t_diag in tracers_defs.h
      logical, parameter         :: tot_change    = .true.  ! write total rate of change in variable during timestep.
      ! End of user inputs                                  ! may not want tot_change (wasteful) since it's just the sum of diagnostic terms.
      ! ***************************************************************

      ! Developer inputs - for budget close checks:
      ! - need to have on diag_uv to use these flags: (not sure that flag still needed for this)
      logical, parameter         :: check_diag    = .true.  ! confirm sum of diagnostic terms    = gross change in u/v per timestep. Useful to ensure diagnostics not broken by different cppdefs.h config.
      logical, parameter         :: check_avg     = .true.  ! confirm sum of diagnostic avg*time = gross change in u/v over entire simulation at time of writing

      ! Preamble:  ![

      logical, public            :: calc_diag               ! flag to calculate diagnostics in equation system.
      logical, public            :: set_dz_old              ! set dz*u/v/w/t at previous timestep, needed for diag calculation when skipping steps if not averages.

      ! U & V MOMENTUM:
      ! ---------------

      integer, public, parameter                          :: nmd = 5    ! number of momentum diag terms

      real, public, allocatable, dimension(:, :, :, :, :) :: Mdiag      ! 2 for u & v
      real, public, allocatable, dimension(:, :, :, :)    :: uvdz_old   ! dz*u or v from previous timestep, need to either store this or Hz to get old u*dz.
      real, public, allocatable, dimension(:, :, :, :)    :: uv_prev    ! store incremental u/v/w (nnew)

      ! if diag_avg - for averaging of diagnostics
      real, allocatable, dimension(:, :, :, :, :)         :: Md_avg     ! averages of diagnostics. Allocatabled (in subroutine) hence memory only if diag_avg is true.

      ! if check_diag - check that diagnostics match gross change in variable
      real, allocatable, dimension(:, :, :, :)            :: uv_dif     ! gross u/v difference between current vs previous timestep u/v
      real, allocatable, dimension(:, :, :)               :: Mdsum      ! sum of Mdiag

      ! if check_avg - could share for u & v for less memory?
      real, allocatable, dimension(:, :, :, :)            :: uvdz_t0    ! original dz*u/v to compare total change in u/v over simulation.
      real, allocatable, dimension(:, :, :, :)            :: uvdz_ever  ! total change in dz*u/v over simulation.
      real, allocatable, dimension(:, :, :, :)            :: Md_a_sum   ! sum of diag averages at timestep
      real, allocatable, dimension(:, :, :, :)            :: Md_a_ever  ! total of all diag averages for entire simulation. Should add up to change in u/v/t over simulation...

      real, public                                        :: Mdtmp      ! temporary scalar for calculations of Mdiag

      integer, parameter, public ::
     &  iprsgr = 1, ! Hydrostatic pressure gradient    - internal loop ranges only (i.e. excluding domain boundary)
     &  icori  = 2, ! Coriolis & curvilinear coords    - internal loop ranges only
     &  iadv   = 3, ! Advection                        - internal loop ranges only
     &  ihmix  = 4, ! Viscous terms (horiz mixing)     - internal loop ranges only
     &  ivmix  = 5  ! Vertical mixing & 2D/3D coupling - internal loop ranges only
      ! ivmix  <--  ! Boundary changes                 - bry only - should store 3dbc changes in 2D arrays not 3D, as just 0 in interior. Goes into ivmix array for memory efficiency.
      ! ivmix  <--  ! 2D-3D correction                 - internal + bry (full domain) loop ranges. Goes into ivmix array for memory efficiency.

      character (len=20),  allocatable, dimension(:) :: vnames
      character (len=120), allocatable, dimension(:) :: lname           ! len=70 allows for 8 characters of 'Average ' if diag_avg

      character (len=26) :: units = 'm^2/s^2    (e.g. du/dt*dz)'        ! Don't go beyond 72nd character else space in units
      character (len=1),   dimension(3) :: uvw = (/ 'u', 'v', 'w' /)    ! prefix for u/v/w outputs

      character (len=120), dimension(:), allocatable :: check_lnm       ! long name for outputs checking diagnostics
      character (len=20),  dimension(:), allocatable :: check_vnm       ! short name ...
      character (len=160), dimension(:), allocatable :: check_alnm      ! long name for outputs checking diagnostics averages
      character (len=30),  dimension(:), allocatable :: check_avnm      ! short name ...

      logical, public :: diag_icori  ! flag to avoid computing 'icori' diagnostic twice in compute_horiz_uv_terms.h

      ! TRACERS:
      ! --------

      ! Since there could be many tracers, it is preferable not to carry any 't_old'
      ! or temporary arrays for each tracer to avoid the memory burden. This differs
      ! from the u,v,w implementation here currently.
      ! Need to consider pipe and river sources. Budget should still close but should label which term receives river/pipe.

      integer, public, parameter                          :: ndt = 5   ! number of tracer diag terms
      real, public, allocatable, dimension(:, :, :, :, :) :: Tdiag     ! tracer diagnostic terms
      real,         allocatable, dimension(:, : )         :: FX4, FE4  ! 4th order advective fluxes to split UP3 = ADV + DISP. Called in k loop so only need 2D array for each k
      real,         allocatable, dimension(:, :, :)       :: Td_tmp    ! temporary array to calculate change in tracer diags
      real,         allocatable, dimension(:, :, :)       :: Hz_old    ! store previous Hz to calculate old dz*t. better than storing dz*t for each tracer. Better to recalculate it?

      ! if diag_avg
      real,         allocatable, dimension(:, :, :, :, :) :: Td_avg    ! averages of tracer diagnostic terms

      ! if check_diag
      real,         allocatable, dimension(:, :, :)       :: T_dsum    ! sum of t diagnostics to compare with gross change

      ! if check_avg
      real,         allocatable, dimension(:, :, :, :)    :: Tdz_t0    ! original dz*t to compare total change in t over simulation.
      real,         allocatable, dimension(:, :, :)       :: Td_a_sum  ! sum of diag averages at timestep
      real,         allocatable, dimension(:, :, :, :)    :: Td_a_ever ! total of all diag averages for entire simulation. Should add up to change in t over simulation...

      integer                            :: nd_trc           ! number of tracers selected tracers of interest (from wrt_t_diag array)
      integer, dimension(:), allocatable :: diag_td_2_NT     ! index of tracer of selected tracers only
      integer, dimension(NT), public     :: diag_NT_2_td=-1  ! convert itrc=1,NT to index of tracer in Tdiag array. Set to bad value -1

      integer, parameter, public ::
     &  it_h_adv   = 1,  ! Horizontal advective fluxes - compute_horiz_tracer_fluxes.h
     &  it_h_disp  = 2,   ! Dissipative part of UP3 advection
     &  it_v_adv   = 3,  ! Vertical advective fluxes   - compute_vert_tracer_fluxes.h
                         ! do we need surface forcing flux as a vertical? Shouldn't be a 3D array anyway
     &  it_v_mix   = 4,  ! Implicit time step for vertical diffusion - step3d_t_ISO.F
     &  it_h_mix   = 5   ! Horizontal (along geopotential surfaces) Laplacian diffusion terms - t3dmix_GP.F
      ! it_v_mix   <--   ! Nonlocal transport flux for unstable (convective) forcing - step3d_t_ISO.F. not used if no LMD_KPP flag, hence at the end. ! step3d_t_ISO.F. This only change T & S so shouldn't be 'NT' in size, a waste.
      ! it_v_mix   <--   ! Boundary changes - t3dbc_im.F. Only at boundary so shared with v_mix (could be any other term but being consistent with u/v choice of v_mix.

      character (len=30),  allocatable, dimension(:) :: tdname
      character (len=120), allocatable, dimension(:) :: tlname  ! extra 10 characters for 'Averaged...'

      character (len=20) :: tunits = ' * m/s   (dC/dt*dz)'

      public init_diags_trc
      public set_diags_t_bc
      public set_diags_t_h_mix
      public set_diags_t_4th_adv

      ! Netcdf outputting:
      ! ------------------
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      integer            :: tn                          ! roms timestep
      character (len=60) :: file_name = '_on_diags.nc'  ! suffix to append to simulation output rootname
      integer            :: ncid=-1, prev_fill_mode


      ! W MOMENTUM:
      ! -----------

      ! NOTE: u & v are 1:N in z, but w is 0:N hence different variables needed
# ifdef NHMG
      integer, public, parameter                         :: nwd = 5    ! number of 'w' momentum diag terms (for NHMG runs)

      real, public, dimension(GLOBAL_2D_ARRAY, 0:N, nwd) :: Wdiag=0    ! 3 for u,v,w
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: wdz_old=0  ! u from previous timestep, need to either store this or Hz to get old u*dz.
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: w_dif=0    ! gross u/v difference between current vs previous timestep u/v
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: Wdsum=0    ! sum of Mdiag
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: wdif=0     ! = w_dif - Mdsum
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: w_prev=0   ! store incremental u/v/w (nnew)

      integer, parameter, public ::
     &  iwprsgr   = 1,
     &  iwhoriadv = 2,
     &  iwvertadv = 3,
     &  iwuv2     = 4,
     &  iwbc      = 5

      character (len=20), dimension(nwd) :: wvname = (/
     &  'w_prsgrd',
     &  'w_horiz_flx',
     &  'w_vert_adv',
     &  'w_2D_3D_coupling',
     &  'w_3dbc'
     &                                               /)

      character (len=60), dimension(nwd) :: wlname = (/
     &  'prsgrd.F',
     &  'compute_horiz_rhs_w_terms.h',
     &  'compute_vert_rhs_w_terms.h',
     &  'step3d_uv2.F',
     &  'w3dbc_im.F'
     &                                               /)

      character (len=20) :: wunits = 'm^2/s (dz*w)'  ! not sure on units here.

      ! Public subroutines:
      public set_diags_w_at_uv1
      public set_diags_w_at_uv2
      public set_diags_w_at_uv2_end
      public set_diags_w_at_bc
# endif /* NHMG */

      public set_diags_uv_coriolis
      public set_diags_uv_advection
      public set_diags_uv_h_mix
      public set_diags_uv_v_mix
      public set_diags_u_bc
      public set_diags_v_bc
      public set_diags_u_2d3d_cor
      public set_diags_v_2d3d_cor
      public set_diags_uv_dz_old
      public wrt_diagnostics !]

      contains

! ----------------------------------------------------------------------
      subroutine init_diags(istr, jstr, iend, jend)  ![
      ! Initilialize allocatable diagnostic arrays.
      ! Set original dz*u/v/w/t for initial timestep, hence only called once.

      implicit none

      ! inputs
      integer, intent(in) :: istr, jstr, iend, jend

      ! local
      integer :: i, j, k, diag

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      if(check_avg .and. .not. diag_avg) then                           ! catch bad config
        write(*,*) 'ERROR - diagnostics.F: check_avg must be false ',
     &  'if diag_avg is false!'
        error stop 'ERROR - diagnostics.F: check_avg!=F if diag_avg==F'
      endif

      if (diag_uv) then

        allocate(    Mdiag(GLOBAL_2D_ARRAY, N, 2, nmd) )
        allocate( uvdz_old(GLOBAL_2D_ARRAY, N, 2)      )
        allocate(  uv_prev(GLOBAL_2D_ARRAY, N, 2)      )
        Mdiag    (istrR:iendR,jstrR:jendR,:,:,:) = 0
        uvdz_old (istrR:iendR,jstrR:jendR,:,:)   = 0
        uv_prev  (istrR:iendR,jstrR:jendR,:,:)   = 0

        allocate( vnames(nmd), lname(nmd) )
        vnames(iprsgr) = '_prsgrd'
        lname (iprsgr) = 'Hydrostatic pressure gradient - prsgrd.F'
        vnames(icori)  = '_coriolis'
        lname (icori)  =
     &    'Coriolis & curvilinear coords - compute_horiz_rhs_uv_term.h'  ! if goes beyond 72 characters the text has big space, hence new line
        vnames(iadv)   = '_advection'
        lname (iadv)   = 'Advection - step3d_uv1.F'
        vnames(ihmix)  = '_h_mix'
        lname (ihmix)  =
     &    'Horizontal mixing (zero beyond sponge layer) - visc3d_GP.F'
        vnames(ivmix) = '_v_mix_and_BC'
        lname (ivmix) = 'Vertical mixing (interior), boundary ' / /
     &    'changes (bry), and 2D/3D coupling (int+bry) - '      / /
     &    'step3d_uv2.F part 1 & u/v3dbc_im.F'

        if (diag_avg) then                                       ! only if diag averaging and only done once
          allocate( Md_avg(GLOBAL_2D_ARRAY, N, 2, nmd) )         ! allocate memory for diag averages
          Md_avg(istrR:iendR,jstrR:jendR,:,:,:) = 0
          do diag=1,nmd                                          ! loop through diag terms
            lname(diag)  = 'Average ' / / trim(lname(diag))      ! and prepend 'Average'
          enddo
        endif

        if (check_diag .or. tot_change) then
          allocate( uv_dif(GLOBAL_2D_ARRAY, N, 2) )
          uv_dif(istrR:iendR,jstrR:jendR,:,:) = 0
        endif

        if (check_diag) then
          allocate(  Mdsum(GLOBAL_2D_ARRAY, N)    )
          Mdsum (istrR:iendR,jstrR:jendR,:)   = 0
        endif

        if (set_dz_old) then
          do k=1,N
            do j=jstrR,jendR   ! Full loop ranges since u-change over every point including bry.
              do i=istr,iendR  ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
                uvdz_old(i,j,k, 1) = u(i,j,k,nstp) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
              enddo
            enddo
            do j=jstr,jendR    ! Full loop ranges since u-change over every point including bry.
              do i=istrR,iendR ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
                uvdz_old(i,j,k, 2) = v(i,j,k,nstp) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
              enddo
            enddo
          enddo
        endif

        if (check_avg) then                                    ! needs to come after if (set_dz_old) as uses uvdz_old
          allocate(  Md_a_sum(GLOBAL_2D_ARRAY, N, 2) )         ! sum of diag averages at timestep
          allocate( Md_a_ever(GLOBAL_2D_ARRAY, N, 2) )         ! continuing total of sum of averages Md_a_ever += Md_a_sum
          allocate(   uvdz_t0(GLOBAL_2D_ARRAY, N, 2) )
          allocate( uvdz_ever(GLOBAL_2D_ARRAY, N, 2) )
          Md_a_ever(istrR:iendR,jstrR:jendR,:,:) = 0           ! Md_a_sum set to 0 elsewhere, so not needed here
          uvdz_t0 = uvdz_old                                   ! store the initial condition for budget confirmation at simulation end
        endif

      endif  ! <-- diag_uv

# ifdef NHMG
      if (set_dz_old) then
        do k=1,N             ! 0 is always 0?
          do j=jstrR,jendR   ! Full loop ranges since u-change over every point including bry.
            do i=istrR,iendR
              if(k==N) then
                wdz_old(i,j,k) = w(i,j,k,nstp) * Hz(i,j,k)
              else
                wdz_old(i,j,k) = w(i,j,k,nstp) * (0.5*( Hz(i,j,k+1)+Hz(i,j,k) ))
              endif
            enddo
          enddo
        enddo
      endif  ! <-- set_dz_old
# endif /* NHMG */

      if (check_diag .or. tot_change) then                           ! used for both uv and tracer
        allocate( check_lnm(3), check_vnm(3)    )
        check_vnm(1) = '_actual'
        check_lnm(1) = 'Actual rate of change in variable during ' / /
     &                 'timestep - e.g. (u_n*dz_n - u_n-1*dz_n-1) / dt'
        check_vnm(2) = '_diag_sum'
        check_lnm(2) = 'Sum of rate of diagnostics terms during '  / /
     &                 'timestep.'
        check_vnm(3) = '_budget'
        check_lnm(3) = 'Diff. of rate of change in variable vs '   / /
     &                 'sum of diags during timestep (budget '     / /
     &                 'close) - e.g. (u_actual-u_diag_sum)'
      endif

      if (check_avg) then                                             ! used for both uv and tracer
        allocate( check_alnm(3), check_avnm(3)     )
        check_avnm(1) = '_ever'
        check_alnm(1) = 'Total change in variable since start of ' / /
     &                  'run - e.g. (u_n*dz_n - u_0*dz_0)'
        check_avnm(2) = '_diag_avg_ever'
        check_alnm(2) = 'Sum of diagnostic averages multiplied '   / /
     &    'by averaging freq. for total time since start of run.'  / /
     &    ' Should match total change in var since model start.'
        check_avnm(3) = '_budget_ever'
        check_alnm(3) = 'Complete budget close using averages '    / /
     &    'since model start until current timestep  - e.g. '      / /
     &    '(u_ever - u_diag_avg_ever)'
      endif

      if (diag_trc) call init_diags_trc(istrR, jstrR, iendR, jendR)

      end subroutine init_diags  !]

! ----------------------------------------------------------------------
      subroutine init_diags_trc(istrR, jstrR, iendR, jendR)  ![
      ! setup tracer arrays only for desired tracers,
      ! since can be memory intensive if all tracers.

      implicit none

      ! inputs
      integer, intent(in)    :: istrR, jstrR, iendR, jendR

      ! local
      integer                :: itrc, cnt=0, diag, td, i, j, k
      integer, dimension(NT) :: temp_t_indx


      allocate( tdname(ndt), tlname(ndt) )
      tdname(it_h_adv)    = '_h_adv'
      tlname(it_h_adv)    =
     &  'Horizontal UP3 purely advective fluxes'                  / /
     &  ' - compute_horiz_tracer_fluxes.h'
      tdname(it_h_disp)   = '_h_disp'
      tlname(it_h_disp)   =
     &  'Horizontal UP3 advection dissipative fluxes'             / /
     &  ' - compute_horiz_tracer_fluxes.h'
      tdname(it_v_adv)    = '_v_adv'
      tlname(it_v_adv)    =
     &  'Vertical advective fluxes - compute_vert_tracer_fluxes.h'
      tdname(it_v_mix)    = '_v_mix_and_BC'
      tlname(it_v_mix)    =
     &  'Implicit vertical diffusion, nonlocal transport flux, '  / /
     &  'and boundary change - step3d_t_ISO.F & t3dbc_im.F'
      tdname(it_h_mix)    = '_h_mix'
      tlname(it_h_mix)    =
     &  'Horizontal (along geopotential surfaces) Laplacian '     / /
     &  'diffusion terms (zero beyond sponge layer) - t3dmix_GP.F'


      do itrc=1,NT
        if (wrt_t_diag(itrc)) then
          cnt=cnt+1                                          ! count tracers we want to calc+write
          diag_NT_2_td(itrc)=cnt                             ! t(itrc) = Tdiag(cnt) - to convert between index scheme
        endif
      enddo

      nd_trc = cnt                                           ! number of tracers of interest
      allocate(  Tdiag( GLOBAL_2D_ARRAY, N, nd_trc, ndt ) )  ! only for tracers we are interested in
      allocate( Td_tmp( GLOBAL_2D_ARRAY, N ) )               ! don't need to set to zero as written over in code
      allocate( Hz_old( GLOBAL_2D_ARRAY, N ) )               ! same for Hz_old
      allocate(    FX4( GLOBAL_2D_ARRAY    ) )               ! called within k loop so only need 2D slice per k
      allocate(    FE4( GLOBAL_2D_ARRAY    ) )
      Tdiag (istrR:iendR, jstrR:jendR, :, :, :) = 0

      allocate( diag_td_2_NT( nd_trc ) )                     ! to convert indices between 't' and 'Tdiag'
      diag_td_2_NT = -1                                      ! set to bad number (<1) as safeguard

      cnt=0
      do itrc=1,NT
        if (wrt_t_diag(itrc)) then
          cnt=cnt+1
          diag_td_2_NT(cnt)=itrc                             ! store the actual tracer index of 't' array
        endif                                                ! since t(NT) but Tdiag(nd_trc). NT >= nd_trc...
      enddo

      if (check_diag) then
        allocate( T_dsum( GLOBAL_2D_ARRAY, N ) )             ! don't need to set to zero as done later in code
      endif

      Hz_old(istrR:iendR,jstrR:jendR,:) = Hz(istrR:iendR,jstrR:jendR,:)  ! this should have a set_dz_old flag

      if (diag_avg) then
        allocate(  Td_avg( GLOBAL_2D_ARRAY, N, nd_trc, ndt ) )
        Td_avg (istrR:iendR, jstrR:jendR, :, :, :) = 0
        do diag=1,ndt
          tlname(diag)  = 'Average ' / / trim(tlname(diag))              ! and prepend 'Average'
        enddo
      endif

      if (check_avg) then
        allocate( Td_a_ever( GLOBAL_2D_ARRAY, N, nd_trc) )
        allocate(    Tdz_t0( GLOBAL_2D_ARRAY, N, nd_trc) )
        allocate(  Td_a_sum( GLOBAL_2D_ARRAY, N        ) )

        do td=1,nd_trc
          itrc=diag_td_2_NT(td)
          do k=1,N; do j=jstrR,jendR; do i=istrR,iendR
            Tdz_t0(i,j,k,td) = t(i,j,k,nstp,itrc) * Hz(i,j,k)  ! Set t*dz at t=0 for overall budget check of averages
          enddo; enddo; enddo
        enddo
      endif

!      if(mynode==0) print *, 'nd_trc=',nd_trc ! debug
!      if(mynode==0) print *, 'diag_td_2_NT(1)=',diag_td_2_NT(1)

      end subroutine init_diags_trc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uv_coriolis( istr, iend, jstr, jend, istrU, jstrV, k, ru, rv )   ![
      ! compute coriolis term from compute_horiz_rhs_uv_terms.h
      ! diag_icori is the logical that prevents this calculation in pre_step3d,
      ! since only want during step3d_uv1.

      implicit none

      ! inputs
      integer,                                       intent(in) :: istr, iend, jstr, jend, istrU, jstrV, k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N),   intent(in) :: ru, rv

      ! local
      integer :: i, j

      do j=jstr,jend                                           ! chose to have 2 loops rather than an if i>=istrU and j> as possible more efficient?
        do i=istrU,iend
          Mdiag(i,j,k,1,icori)=ru(i,j,k)-Mdiag(i,j,k,1,iprsgr) ! correolis & curvilinear coords
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          Mdiag(i,j,k,2,icori)=rv(i,j,k)-Mdiag(i,j,k,2,iprsgr) ! correolis & curvilinear coords
        enddo
      enddo

      end subroutine set_diags_uv_coriolis  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uv_advection( istr, iend, jstr, jend, istrA, j, dim, DC, ruv, suvstr, uv )  ![
      ! Set the diagnostic terms for prsgrd, coriolis, and advection.
      ! Works for both u and v.
      ! Prsgrd and coriolis need to also be convert to dz * u/v units
      ! Called from step3d_uv1.F

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrA, j, dim         ! 'istrA' for any as depends on inputs, which 'dim' dimension.
                                                                           ! istr, jstr, etc needed for GLOBAL_2D...
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N), intent(in) :: DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N),   intent(in) :: ruv     ! ru or rv
      real, dimension(GLOBAL_2D_ARRAY),              intent(in) :: suvstr  ! sustr or svstr
      real, dimension(GLOBAL_2D_ARRAY,N,3),          intent(in) :: uv      ! u or v

      ! local
      integer :: i, k

      do k=1,N
        do i=istrA,iend
          Mdiag(i,j,k,dim,iprsgr) = Mdiag(i,j,k,dim,iprsgr) * DC(i,0) / dt
          Mdiag(i,j,k,dim,icori)  = Mdiag(i,j,k,dim,icori)  * DC(i,0) / dt
          uv_prev(i,j,k,dim)      = uv(i,j,k,nnew)
          Mdiag(i,j,k,dim,iadv)   = DC(i,0)*ruv(i,j,k) / dt                ! all advection - includes 2nd part of compute_horiz_uv & compute_vert_uv.
     &                            - ( Mdiag(i,j,k,dim,iprsgr) + Mdiag(i,j,k,dim,icori) )
        enddo
      enddo

      ! Add surface stress to vertical mixing rather:
      do i=istrA,iend
        Mdiag(i,j,N,dim,ivmix) = suvstr(i,j)  ! * dt / dt cancels...
      enddo

      end subroutine set_diags_uv_advection  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uv_h_mix( istr, iend, jstr, jend, isA, jsA, dim, uv )  ![
      ! set diagnostics for u or v in visc3d_GP.F
      ! interior points only (i.e. excluding bry)
      ! Note this diagnostic is zero beyond the sponge layer (in the interior),
      ! but can be non-zero everywhere if model has non-zero background diffusion,
      ! hence loop ranges over full interior.

      implicit none

      ! inputs
      integer,                              intent(in) :: istr, iend, jstr, jend, isA, jsA, dim  ! isA jsA are starts depending on u or v variable
      real, dimension(GLOBAL_2D_ARRAY,N,3), intent(in) :: uv  ! u or v

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jsA,jend
          do i=isA,iend
            Mdiag(  i,j,k,dim,ihmix) = (uv(i,j,k,nnew) - uv_prev(i,j,k,dim)) / dt  ! units already in dz*u
            uv_prev(i,j,k,dim)       =  uv(i,j,k,nnew)
          enddo
        enddo
      enddo

      end subroutine set_diags_uv_h_mix  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uv_v_mix( istr, iend, jstr, jend, istrU, jstrV )  ![
      ! set diagnostics for u or v in step3d_uv2.F
      ! interior points only (i.e. excluding bry)
      ! the calc in step3d_uv2.F includes variable 'DC', hence name set_diags_uv_v_mix

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrU, jstrV

      ! local
      integer :: i, j, k

      ! Since this term includes 2D3D & boundary, need to set boundary to zero to
      ! allow for Mdiag( i,j,k,1,ivmix) = Mdiag( i,j,k,1,ivmix) + .... later on.
      ! entire column of west  or east  boundary would use jstr-1:jend+1
      ! entire row    of south or north boundary would use istr-1:iend+1

      if(WESTERN_EDGE)  Mdiag(istr-1       ,jstr-1:jend+1,:,2,ivmix) = 0 ! west  (istr-1 for v bry) Need split else surface stress over-written.
      if(WESTERN_EDGE)  Mdiag(istr         ,jstr-1:jend+1,:,1,ivmix) = 0 ! west  (istr   for u bry)
      if(EASTERN_EDGE)  Mdiag(iend+1       ,jstr-1:jend+1,:,:,ivmix) = 0 ! east  (iend+1=iendR=bry)

      if(SOUTHERN_EDGE) Mdiag(istr-1:iend+1,jstr-1       ,:,1,ivmix) = 0 ! south (jstr-1 for u bry)
      if(SOUTHERN_EDGE) Mdiag(istr-1:iend+1,jstr         ,:,2,ivmix) = 0 ! south (jstr   for v bry)
      if(NORTHERN_EDGE) Mdiag(istr-1:iend+1,jend+1       ,:,:,ivmix) = 0 ! north (jend+1=jendR=bry)

      do k=1,N-1           ! N-1 as treat k=N differently below
        do j=jstr,jend     ! interior of subdomains, exludes global domain boundaries
          do i=istrU,iend
            Mdtmp                 = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) ! convert to dz*u
            Mdiag( i,j,k,1,ivmix) = (Mdtmp - uv_prev(i,j,k,1)) / dt                        ! avoids recomputing above
            uv_prev(i,j,k,1)      =  Mdtmp
          enddo
        enddo
        do j=jstrV,jend    ! interior of subdomains, exludes global domain boundaries
          do i=istr,iend
            Mdtmp                 = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
            Mdiag( i,j,k,2,ivmix) = (Mdtmp - uv_prev(i,j,k,2)) / dt
            uv_prev(i,j,k,2)      =  Mdtmp
          enddo
        enddo
      enddo

      ! Since ivmix also contains surface stress from set_diags_uv_advection we
      ! need Mdiag( i,j,k,1,ivmix) = Mdiag( i,j,k,1,ivmix) + ...
      do j=jstr,jend     ! interior of subdomains, exludes global domain boundaries
        do i=istrU,iend
          Mdtmp                 = u(i,j,N,nnew) * (0.5*( Hz(i-1,j,N)+Hz(i,j,N) )) ! convert to dz*u
          Mdiag( i,j,N,1,ivmix) = Mdiag( i,j,N,1,ivmix)
     &                          + (Mdtmp - uv_prev(i,j,N,1)) / dt
          uv_prev(i,j,N,1)      =  Mdtmp
        enddo
      enddo
      do j=jstrV,jend    ! interior of subdomains, exludes global domain boundaries
        do i=istr,iend
          Mdtmp                 = v(i,j,N,nnew) * (0.5*( Hz(i,j-1,N)+Hz(i,j,N) ))
          Mdiag( i,j,N,2,ivmix) = Mdiag( i,j,N,2,ivmix)
     &                          + (Mdtmp - uv_prev(i,j,N,2)) / dt
          uv_prev(i,j,N,2)      =  Mdtmp
        enddo
      enddo

!      if(mynode==0) print *, 'uv2DC!!!!!!!!!'  ! debug only

      end subroutine set_diags_uv_v_mix  !]

! ----------------------------------------------------------------------
      subroutine set_diags_u_bc( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in u
      ! called from step3d_uv2.F right after call to u3dbc_im.F

      ! boundary change is stored in ivmix diagnostic to save memory.

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i

      ! currently this is not perfectly setup for openMP since we have ':' instead
      ! of istr:iend, etc. However, since it is just a 2D operation at bry, it is ok? for now.
      ! Also, the 4 boundary corners of the global domain are computed twice, i.e. bottom
      ! of west is the same as left of south. But since we have Mdiag(istr  ,:,:,1,ivmix) = ...
      ! and not Mdiag(istr  ,:,:,1,ivmix) = Mdiag(istr  ,:,:,1,ivmix) + ... it doesn't matter.
      ! Consider using jstrR instead of jstr-1 for cleanliness.

      if (WESTERN_EDGE)  uv_prev(istr  ,:,:,1) = u(istr  ,:,:,nnew) * 0.5*(Hz(istr-1,:,:)+Hz(istr  ,:,:))
      if (EASTERN_EDGE)  uv_prev(iend+1,:,:,1) = u(iend+1,:,:,nnew) * 0.5*(Hz(iend  ,:,:)+Hz(iend+1,:,:))

      if (WESTERN_EDGE)  Mdiag(istr  ,:,:,1,ivmix) = ( uv_prev(istr  ,:,:,1) - uvdz_old(istr  ,:,:, 1) )/dt ! uv_prev is already u_current, can do this because have uvdz_old already
      if (EASTERN_EDGE)  Mdiag(iend+1,:,:,1,ivmix) = ( uv_prev(iend+1,:,:,1) - uvdz_old(iend+1,:,:, 1) )/dt

      do i=istr,iend+1 ! use loop here else indexing would be even messier
        if (SOUTHERN_EDGE) uv_prev(i,jstr-1,:,1)=u(i,jstr-1,:,nnew)*0.5*(Hz(i-1,jstr-1,:)+Hz(i,jstr-1,:))
        if (NORTHERN_EDGE) uv_prev(i,jend+1,:,1)=u(i,jend+1,:,nnew)*0.5*(Hz(i-1,jend+1,:)+Hz(i,jend+1,:))
      enddo
      if (SOUTHERN_EDGE) Mdiag(:,jstr-1,:,1,ivmix) = ( uv_prev(:,jstr-1,:,1) - uvdz_old(:,jstr-1,:, 1) )/dt
      if (NORTHERN_EDGE) Mdiag(:,jend+1,:,1,ivmix) = ( uv_prev(:,jend+1,:,1) - uvdz_old(:,jend+1,:, 1) )/dt

      end subroutine set_diags_u_bc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_v_bc( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in v
      ! called from step3d_uv2.F right after call to v3dbc_im.F

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: j

      ! currently this is not perfectly setup for openMP since we have ':' instead
      ! of istr:iend, etc. However, since it is just a 2D operation at bry, it is ok? for now.

      do j=jstr,jend+1 ! use loop here else indexing would be even messier
        if (WESTERN_EDGE)  uv_prev(istr-1,j,:,2) = u(istr-1,j,:,nnew) * 0.5*(Hz(istr-1,j-1,:)+Hz(istr-1,j,:))
        if (EASTERN_EDGE)  uv_prev(iend+1,j,:,2) = u(iend+1,j,:,nnew) * 0.5*(Hz(iend+1,j-1,:)+Hz(iend+1,j,:))
      enddo
      if (WESTERN_EDGE)  Mdiag(istr-1,:,:,2,ivmix) = ( uv_prev(istr-1,:,:,2)-uvdz_old(istr-1,:,:,2) )/dt ! uv_prev is already u_current, can do this because have uvdz_old already
      if (EASTERN_EDGE)  Mdiag(iend+1,:,:,2,ivmix) = ( uv_prev(iend+1,:,:,2)-uvdz_old(iend+1,:,:,2) )/dt

      if (SOUTHERN_EDGE) uv_prev(:,jstr  ,:,2)=u(:,jstr  ,:,nnew)*0.5*(Hz(:,jstr-1,:)+Hz(:,jstr  ,:))
      if (NORTHERN_EDGE) uv_prev(:,jend+1,:,2)=u(:,jend+1,:,nnew)*0.5*(Hz(:,jend  ,:)+Hz(:,jend+1,:))
      if (SOUTHERN_EDGE) Mdiag(:,jstr  ,:,2,ivmix) = ( uv_prev(:,jstr  ,:,2)-uvdz_old(:,jstr  ,:,2) )/dt
      if (NORTHERN_EDGE) Mdiag(:,jend+1,:,2,ivmix) = ( uv_prev(:,jend+1,:,2)-uvdz_old(:,jend+1,:,2) )/dt

      end subroutine set_diags_v_bc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_u_2d3d_cor( istr, iend, jstr, jend, istrU, iendR, jstrR, jendR )  ![
      ! Final 2D-3D correction, part 2 of step3d_uv2.F.

      ! NOTE: Already done 2D-3D correction in v_mix step, then we updated the boundary
      !       and then do another 2D-3D correction. Values over the interior only change
      !       by e-17, so therefore exclude interior for efficiency.
      !       The budget closing is not quite as accurate but we are talking about
      !       e-16 difference.

      ! We add 2D3D correction to v_mix diag term for memory efficiency

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrU, iendR, jstrR, jendR  ! iend, jstr, jend needed for _EDGE's...

      ! local
      integer :: i, j, k

      if (WESTERN_EDGE) then  ! only do boundaries for efficiency
        i=istr
        Mdiag(i,:,:,1,ivmix) = Mdiag(i,:,:,1,ivmix) +
     &    ( u(i,:,:,nnew) * (0.5*( Hz(i-1,:,:)+Hz(i,:,:) )) - uv_prev(i,:,:,1) )/dt
      endif
      if (EASTERN_EDGE) then
        i=iendR
        Mdiag(i,:,:,1,ivmix) = Mdiag(i,:,:,1,ivmix) +
     &    ( u(i,:,:,nnew) * (0.5*( Hz(i-1,:,:)+Hz(i,:,:) )) - uv_prev(i,:,:,1) )/dt
      endif
      if (SOUTHERN_EDGE) then
        j=jstrR
        do i=istrU,iend       ! will cover full range, but avoid double counting SW and SE corner as already done for EAST or WEST above
        Mdiag(i,j,:,1,ivmix) = Mdiag(i,j,:,1,ivmix) +
     &    ( u(i,j,:,nnew) * (0.5*( Hz(i-1,j,:)+Hz(i,j,:) )) - uv_prev(i,j,:,1) )/dt
        enddo
      endif
      if (NORTHERN_EDGE) then
        j=jendR
        do i=istrU,iend       ! will cover full range, but avoid double counting NW and NE corner as already done for EAST or WEST above
        Mdiag(i,j,:,1,ivmix) = Mdiag(i,j,:,1,ivmix) +
     &    ( u(i,j,:,nnew) * (0.5*( Hz(i-1,j,:)+Hz(i,j,:) )) - uv_prev(i,j,:,1) )/dt
        enddo
      endif

      if (check_diag .or. tot_change) then
        do k=1,N              ! Full loop ranges since u-change over every point including bry.
          do j=jstrR,jendR    ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            do i=istr,iendR   ! can't just do (u(nnew)-u(nstp)*Hz here because actual need the old Hz for u(nstp) as Hz has evolved.
              uv_dif(i,j,k, 1) = ( u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - uvdz_old(i,j,k,1) )/dt ! do this here as units are correct
            enddo
          enddo
        enddo
      endif

      end subroutine set_diags_u_2d3d_cor  !]

! ----------------------------------------------------------------------
      subroutine set_diags_v_2d3d_cor( istr, iend, jstr, jend, istrR, iendR, jstrV, jendR )  ![
      ! Final 2D-3D correction, part 2 of step3d_uv2.F.

      ! NOTE: Already done 2D-3D correction in v_mix step, then we updated the boundary
      !       and then do another 2D-3D correction. Values over the interior only change
      !       by e-17, so therefore exclude interior for efficiency.
      !       The budget closing is not quite as accurate but we are talking about
      !       e-16 difference.

      ! We add 2D3D correction to v_mix diag term for memory efficiency

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrR, iendR, jstrV, jendR

      ! local
      integer :: i, j, k

      if (WESTERN_EDGE) then  ! only do boundaries for efficiency
        i=istrR
        do j=jstrV,jend       ! will cover full range, but avoid double counting SW and NW corner as done for SOUTH/NORTH below
        Mdiag(i,j,:,2,ivmix) = Mdiag(i,j,:,2,ivmix) +
     &  ( v(i,j,:,nnew)*(0.5*( Hz(i,j-1,:)+Hz(i,j,:) ))-uv_prev(i,j,:,2) )/dt
        enddo
      endif
      if (EASTERN_EDGE) then
        i=iendR
        do j=jstrV,jend       ! will cover full range, but avoid double counting SE and NE corner as done for SOUTH/NORTH below
        Mdiag(i,j,:,2,ivmix) = Mdiag(i,j,:,2,ivmix) +
     &  ( v(i,j,:,nnew)*(0.5*( Hz(i,j-1,:)+Hz(i,j,:) ))-uv_prev(i,j,:,2) )/dt
        enddo
      endif
      if (SOUTHERN_EDGE) then
        j=jstr
        Mdiag(:,j,:,2,ivmix) = Mdiag(:,j,:,2,ivmix) +
     &  ( v(:,j,:,nnew)*(0.5*( Hz(:,j-1,:)+Hz(:,j,:) ))-uv_prev(:,j,:,2) )/dt
      endif
      if (NORTHERN_EDGE) then
        j=jendR
        Mdiag(:,j,:,2,ivmix) = Mdiag(:,j,:,2,ivmix) +
     &  ( v(:,j,:,nnew)*(0.5*( Hz(:,j-1,:)+Hz(:,j,:) ))-uv_prev(:,j,:,2) )/dt
      endif

      if (check_diag .or. tot_change) then
        do k=1,N
          do j=jstr,jendR
            do i=istrR,iendR
              uv_dif(i,j,k, 2) = ( v(i,j,k,nnew)*(0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))-uvdz_old(i,j,k,2) )/dt ! do this here as units are correct
            enddo
          enddo
        enddo
      endif

      end subroutine set_diags_v_2d3d_cor  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uv_dz_old(istr, jstr, istrR, iendR, jstrR, jendR)  ![
      ! Set dz*u or v at previous timestep to calculate gross change in
      ! u or v overtime step, and to calculate change at boundary.
      ! Called after everything else in step3d_uv2.F.

      implicit none

      ! inputs
      integer, intent(in) :: istr, jstr, istrR, iendR, jstrR, jendR

      ! local
      integer :: i, j, k

      ! This is currently used for bry calculation so is also needed for averaging.
      ! Also for 'u_diff', but that's not needed for averaging.
      ! Should put an if statement here for diags_avg and only compute bry if averages...

      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            uvdz_old(i,j,k, 1) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
          enddo
        enddo
      enddo

      do k=1,N
        do j=jstr,jendR
          do i=istrR,iendR
            uvdz_old(i,j,k, 2) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
          enddo
        enddo
      enddo

      ! Don't forget 'w' here.

!      if(mynode==0) print *, 'DZ_OLD  $$$$$$$$$$'  ! debug

      end subroutine set_diags_uv_dz_old  !]

! ----------------------------------------------------------------------
# ifdef NHMG
      subroutine set_diags_w_at_uv1( istr, iend, jstr, jend, DC, rw )  ![
      ! Set the diagnostic terms for prsgrd, horiz & vert advection.
      ! Prsgrd and horiz adv still need to also be convert to dz * w

      implicit none

      ! inputs
      integer,                                       intent(in) :: istr, iend, jstr, jend
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N), intent(in) :: DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N), intent(in) :: rw

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstr,jend
          do i=istr,iend
            Wdiag(i,j,k,iwprsgr)   = Wdiag(i,j,k,iwprsgr)   * DC(i,0)  ! need to convert to dz*w now that we have DC
            Wdiag(i,j,k,iwhoriadv) = Wdiag(i,j,k,iwhoriadv) * DC(i,0)
            Wdiag(i,j,k,iwvertadv) = DC(i,0)*rw(i,j,k)
     &                              - ( Wdiag(i,j,k,iwprsgr) + Wdiag(i,j,k,iwhoriadv) ) ! subtract previous 2 terms included in total
            w_prev(i,j,k)          = w(i,j,k,nnew)                     ! store for next change
          enddo
        enddo
      enddo

      end subroutine set_diags_w_at_uv1  !]

! ----------------------------------------------------------------------
      subroutine set_diags_w_at_uv2( istr, iend, jstr, jend )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k

      ! This could well be wrong. Need to confirm with Jeroen.
      do k=1,N-1
        do j=jstr,jend  ! I changed istr to jstr...
          do i=Istr,Iend
            Wdiag(i,j,k,iwuv2) = w(i,j,k,nnew) * 0.5*( Hz(i,j,k+1)+Hz(i,j,k) ) - w_prev(i,j,k)
          enddo
        enddo
      enddo

      ! Need to do N seperately since can't do 0.5 * Hz(i,j,k+1)
      ! Do seperate loop rather than having if statement in loop? Not sure what is more efficient
      ! ifort compiler report could tell me if it fuses these loops...
      do j=jstr,jend  ! I changed istr to jstr...
        do i=Istr,Iend
          Wdiag(i,j,N,iwuv2) = w(i,j,N,nnew) * Hz(i,j,N) - w_prev(i,j,N)
        enddo
      enddo

      end subroutine set_diags_w_at_uv2  !]

! ----------------------------------------------------------------------
      subroutine set_diags_w_at_bc( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in w
      ! called from step3d_uv2.F right after call to w3dbc_im.F

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k

      if (WESTERN_EDGE) then
        i=istr-1
        do k=1,N-1            ! N-1 because of k+1
          do j=jstr-1,jend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (EASTERN_EDGE) then
        i=iend+1
        do k=1,N-1
          do j=jstr-1,jend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (SOUTHERN_EDGE) then
        j=jstr-1
        do k=1,N-1
          do i=istr-1,iend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (NORTHERN_EDGE) then
        j=jend+1
        do k=1,N-1
          do i=istr-1,iend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      end subroutine set_diags_w_at_bc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_w_at_uv2_end( istrR, iendR, jstrR, jendR)  ![
      implicit none

      ! inputs
      integer, intent(in) :: istrR, iendR, jstrR, jendR ! changed the name to use for u or v

      ! local
      integer :: i, j, k, tmp

      do k=1,N-1
        do j=jstrR,jendR
          do i=istrR,iendR
            ! Full loop ranges since u-change over every point including bry.
            tmp = w(i,j,k,nnew) * 0.5 * ( Hz(i,j,k+1)+Hz(i,j,k) )
            w_dif(i,j,k)   = tmp - wdz_old(i,j,k)
            wdz_old(i,j,k) = tmp
          enddo
        enddo
      enddo

      do j=jstrR,jendR
        do i=istrR,iendR
          tmp = w(i,j,N,nnew) * Hz(i,j,N)
          w_dif(i,j,N)   = tmp - wdz_old(i,j,N)
          wdz_old(i,j,N) = tmp
        enddo
      enddo

      end subroutine set_diags_w_at_uv2_end  !]
# endif /* NHMG */

! ----------------------------------------------------------------------
      subroutine set_diags_t_bc(istr, iend, jstr, jend, itrc, td)  ![
      ! Set tracer boundary changes

      use tracers, only: t
      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend, itrc, td

      ! local
      integer :: i, j

      ! We use the boundary points of the it_v_mix array to be memory efficient, as
      ! vertical mixing diag doesn't use the boundary points.

      ! currently this is not perfectly setup for openMP since we have ':' instead
      ! of istr:iend, etc. However, since it is just a 2D operation at bry, it is ok? for now.

!      if (mynode==0) print *, 'TRACER BC /////////' ! debug

      if (WESTERN_EDGE) then
        i=istr-1
        Tdiag(i,:,:,td,it_v_mix) = ( t(i,:,:,nnew,itrc) *     Hz(i,:,:)
     &                             - t(i,:,:,nstp,itrc) * Hz_old(i,:,:) ) /dt  ! divide by dt to get rate /sec
      endif
      if (EASTERN_EDGE) then
        i=iend+1
        Tdiag(i,:,:,td,it_v_mix) = ( t(i,:,:,nnew,itrc) *     Hz(i,:,:)
     &                             - t(i,:,:,nstp,itrc) * Hz_old(i,:,:) ) /dt
      endif
      if (SOUTHERN_EDGE) then
        j=jstr-1
        Tdiag(:,j,:,td,it_v_mix) = ( t(:,j,:,nnew,itrc) *     Hz(:,j,:)
     &                             - t(:,j,:,nstp,itrc) * Hz_old(:,j,:) ) /dt
      endif
      if (NORTHERN_EDGE) then
        j=jend+1
        Tdiag(:,j,:,td,it_v_mix) = ( t(:,j,:,nnew,itrc) *     Hz(:,j,:)
     &                             - t(:,j,:,nstp,itrc) * Hz_old(:,j,:) ) /dt
      endif

      end subroutine set_diags_t_bc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_t_h_mix(istr,iend,jstr,jend, itrc)  ![
      ! Horizontal mixing from t3dmix_GP.F
      ! zero in interior beyond sponge layer if zero background diffusion!

      use tracers, only: t
      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend, itrc

      ! local
      integer :: i, j, k, td

      td=diag_NT_2_td(itrc)
      do k=1,N
        do j=jstr,jend
          do i=istr,iend
            Tdiag(i,j,k,td,it_h_mix) = t(i,j,k,nnew,itrc)
          enddo
        enddo
      enddo

      end subroutine set_diags_t_h_mix  !]

! ----------------------------------------------------------------------
      subroutine set_diags_t_4th_adv( istr, iend, jstr, jend, k, itrc )  ![
      ! compute the 4th order advection.
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th


      ![ | T(i-2) | T(i-1) | T(i)   | T(i+1) |
      !                    ^        ^
      !                   UFX(i)   UFX(i+1)
      !                   TFX(i)   TFX(i+1)
      !
      ! where UFX and TFX are at u-points and T is at rho-points.
      !
      ! 4th order (centered) fluxes:
      !
      !   TFx(i) = UFx(i) * ( - T(i-2) + 7T(i-1) + 7T(i) - T(i+1) ) / 12
      !
      ! Hence divergence of the fluxes:
      !
      !   T_adv_4th = dt*pm*pn*( TFX(i+1) - TFX(i)
      !                         +TFE(i+1) - TFE(i) )
      !]

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, k, itrc

      ! local
      integer :: i, j, td

      td=diag_NT_2_td(itrc)
      do j=jstr,jend+1   ! +1 since u2rho averaging for Tdiag(it_h_disp) below
        do i=istr,iend+1

          FX4(i,j) = 0.0833333333333 * FlxU(i,j,k)  ! hoping it's fine to reuse FX4 for every tracer
     &               * (-t(i-2,j  ,k,nrhs,itrc) + 7*t(i-1,j  ,k,nrhs,itrc)
     &                  -t(i+1,j  ,k,nrhs,itrc) + 7*t(i  ,j  ,k,nrhs,itrc) )

          FE4(i,j) = 0.0833333333333 * FlxV(i,j,k)  ! hoping it's fine to reuse FX4 for every tracer
     &               * (-t(i  ,j-2,k,nrhs,itrc)+7*t(i  ,j-1,k,nrhs,itrc)
     &                  -t(i  ,j+1,k,nrhs,itrc)+7*t(i  ,j  ,k,nrhs,itrc) )

        enddo
      enddo

      do j=jstr,jend
        do i=istr,iend
          Tdiag(i,j,k,td,it_h_adv) = -pm(i,j)*pn(i,j)*( FX4(i+1,j  )-FX4(i,j)   ! no *dt since rate of flux divergence
     &                                                 +FE4(i  ,j+1)-FE4(i,j) )
        enddo
      enddo

      end subroutine set_diags_t_4th_adv  !]

! ----------------------------------------------------------------------
      subroutine set_diags_avg_tile( istr, iend, jstr, jend )  ![
      ! set the averages of diagnostic variables

      ! INFO: ![
      ! This method is arguably more efficient than set_avg.F, however, it requires more
      ! lines of code. Since diagnostic averaging will be heavy, more efficient code is argueably
      ! more important.

      ! Averaging rate is dictated by rec_rate.
      ! 'Md_avg' stores a running total of diag terms 'Mdiag' until the rec_rate
      ! is reached, at which point Md_avg is divided by the rec_rate to give the
      ! the average of the diagnostics terms.
      ! Md_avg is then refreshed for the next averaging period with Md_avg = Mdiag

      ! 'check_avg' is used to do a complete budget close since model start until current
      ! timestep using the averages of the diag terms.
      ! 'Md_a_sum' is the sum of the diag averages during that rec_rate interval.
      ! 'Md_a_ever' contains the continuous tally of diag averages since model start
      ! multiplied by rec_rate * dt, this equal the total change in the variable
      ! since model start and hence serve to check the overall budget.
      ! Since Md_a_sum is the average rate per second of the diagnostics, multiplying by
      ! rec_rate (i.e. timesteps in averaging period) * dt, gives the total change in the
      ! variable over the averaging period, added to Md_a_ever maintains the change since
      ! model start.

      ! do iMd=1,ivmix                  ! all diagnostic terms with internal loop ranges
      ! for greatest efficiency should only do the loop ranges (internal or full) of the diag terms, but doing full
      ! ranges for now just to get working.
      !]

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k, iMd, dim

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

!      if(mynode==0) print *, '--- SETTING DIAGS AVERAGES' ! debug only

      do dim=1,2

        if (mod(tn,rec_rate)==0) then      ! timestep / rec_rate

          if (check_avg) Md_a_sum(istrR:iendR,jstrR:jendR,:,dim) = 0

          do iMd=1,nmd                     ! all diagnostic terms with internal loop ranges
            do k=1,N                       ! here we divide the average by the averaging rate
              do j=jstrR,jendR
                do i=istrR,iendR
                  Md_avg(i,j,k,dim,iMd) = ( Md_avg(i,j,k,dim,iMd)+Mdiag(i,j,k,dim,iMd) ) / rec_rate
                  if (check_avg)
     &              Md_a_sum(i,j,k,dim) = Md_a_sum(i,j,k,dim) + Md_avg(i,j,k,dim,iMd)  ! if with compile time parameter
                enddo
              enddo
            enddo
          enddo

          if (check_avg) then                                                   ! see top of subroutine for more info
            do k=1,N; do j=jstrR,jendR; do i=istrR,iendR                        ! Md_a_ever, could write it now then wouldn't need
              Md_a_ever(i,j,k,dim) = Md_a_ever(i,j,k,dim)                       ! for both u and v? Same for Md_a_sum
     &                             +  Md_a_sum(i,j,k,dim) * rec_rate * dt       ! keep an endless running total for budget closing
            enddo; enddo; enddo                                                 ! * rec_rate since tot = avg * t
            call set_diags_uvdz_ever( istr, jstr, istrR, iendR, jstrR, jendR )
          endif

!          if(mynode==0 .and. dim==1) print *, '----> DIVIDE FOR AVERAGE'

        else

          if (mod(tn-1,rec_rate)==0) then   ! this is next step after we have averaged.

            ! or could have Md_avg = 0 at the top of routine but this is more efficient but more code?
            do iMd=1,nmd                    ! all diagnostic terms with internal loop ranges
              do k=1,N                      ! here we divide the average by the averaging rate
                do j=jstrR,jendR
                  do i=istrR,iendR
                    Md_avg(i,j,k,dim,iMD) = Mdiag(i,j,k,dim,iMD)
                  enddo
                enddo
              enddo
            enddo

          else                              ! all other timesteps

            do iMd=1,nmd                    ! all diagnostic terms with internal loop ranges
              do k=1,N                      ! here we divide the average by the averaging rate
                do j=jstrR,jendR
                  do i=istrR,iendR
                    Md_avg(i,j,k,dim,iMD) = Md_avg(i,j,k,dim,iMD)+Mdiag(i,j,k,dim,iMD)
                  enddo
                enddo
              enddo
            enddo

          endif  ! <-- mod(tn-1,rec_rate)==0

        endif    ! <-- mod(tn  ,rec_rate)==0

      enddo    ! <-- dim

      end subroutine set_diags_avg_tile  !]

! ----------------------------------------------------------------------
      subroutine set_diags_t_avg_tile(diagT, istr, iend, jstr, jend, td, iTd, itrc)  ![
      ! set the averages of tracer diagnostic variables

      ! NOTE: this won't work for rec_rate=1 because average will never refresh,
      ! not that you would want rec_rate=1 for averages anyway.

      implicit none

      ! inputs
      real, dimension(GLOBAL_2D_ARRAY,N), intent(in) :: diagT      ! 3D part of Tdiag for individual tracer

      integer, intent(in) :: istr, iend, jstr, jend, td, iTd, itrc ! tn  = timestep
                                                                   ! td  = tracer (local diag index)
      ! local                                                      ! iTd = diagnostic term
      integer           :: i, j, k, ierr                           ! itrc = tracer index as per 't' array
      character(len=30) :: text

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

!      if(mynode==0) print *, '^^^ SETTING TRC DIAGS AVERAGES' ! debug only

      if (mod(tn,rec_rate)==0) then                ! timestep / rec_rate

        if (check_avg .and. iTd==it_h_disp)
     &    Td_a_sum(istrR:iendR,jstrR:jendR,:) = 0  ! reset for new tracer: iTd==it_h_disp as
                                                   ! must be first diag from wrt_tracer_diags_tile
        do k=1,N
          do j=jstrR,jendR                         ! here we divide the average
            do i=istrR,iendR                       ! by the averaging rate
              Td_avg(i,j,k,td,iTd) = ( Td_avg(i,j,k,td,iTd)+diagT(i,j,k) ) / rec_rate
              if (check_avg)
     &          Td_a_sum(i,j,k) = Td_a_sum(i,j,k) + Td_avg(i,j,k,td,iTd)     ! if with parameter
            enddo
          enddo
        enddo

        text=trim(t_vname(itrc)) / / tdname(iTd)
        call nc_write_var(ncid, Td_avg(:,:,:,td,iTd), N, text, rp_var, output_rec, ierr)

        if (check_avg .and. iTd==it_v_mix) then                                       ! budget check done after final diag term from wrt_tracer_diags_tile (iTd==it_v_mix)
          do k=1,N; do j=jstrR,jendR; do i=istrR,iendR
            Td_a_ever(i,j,k,td) = Td_a_ever(i,j,k,td) + Td_a_sum(i,j,k) *rec_rate*dt  ! keep an endless running total for budget closing
          enddo; enddo; enddo                                                         ! * rec_rate * dt since total change = avg * time (where time = rec_rate*dt)
          text=trim(t_vname(itrc)) / / check_avnm(2)
          call nc_write_var(ncid, Td_a_ever(:,:,:,td), N, text, rp_var, output_rec, ierr)

          call set_diags_tdz_ever( istrR, iendR, jstrR, jendR, td, itrc )             ! calc total change in tracer since model start.
          text=trim(t_vname(itrc)) / / check_avnm(1)
          call nc_write_var(ncid, Td_a_sum(:,:,:), N, text, rp_var, output_rec, ierr)

          do k=1,N; do j=jstrR,jendR; do i=istrR,iendR                                ! check budget still closing since model start
            Td_a_sum(i,j,k) = Td_a_sum(i,j,k) - Td_a_ever(i,j,k,td)                   ! Td_a_sum on RHS is tdz_ever array as a difference array
          enddo; enddo; enddo                                                         ! Td_a_sum on LHS is re-used as a difference array for budget close.
          text=trim(t_vname(itrc)) / / check_avnm(3)
          call nc_write_var(ncid, Td_a_sum(:,:,:), N, text, rp_var, output_rec, ierr)
        endif

!        if(mynode==0) print *, '----> TRACER DIVIDE FOR AVERAGE'

      else

        if (mod(tn-1,rec_rate)==0) then   ! this is next step after we have averaged.

          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                Td_avg(i,j,k,td,iTd) = diagT(i,j,k)
              enddo
            enddo
          enddo

        else                              ! all other timesteps

          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                Td_avg(i,j,k,td,iTd) = Td_avg(i,j,k,td,iTd) + diagT(i,j,k)
              enddo
            enddo
          enddo

        endif  ! <-- mod(tn-1,rec_rate)==0

      endif    ! <-- mod(tn  ,rec_rate)==0

      end subroutine set_diags_t_avg_tile  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uvdz_ever( istr, jstr, istrR, iendR, jstrR, jendR )  ![
      ! need current dz*u/v to compare gross changes in u with sum of averages
      ! used only with 'check_avg' logical

      implicit none

      ! inputs
      integer, intent(in) :: istr, jstr, istrR, iendR, jstrR, jendR

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            uvdz_ever(i,j,k, 1) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
     &                          - uvdz_t0(i,j,k,1)
          enddo
        enddo
      enddo

      do k=1,N
        do j=jstr,jendR
          do i=istrR,iendR
            uvdz_ever(i,j,k, 2) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
     &                          - uvdz_t0(i,j,k,2)
          enddo
        enddo
      enddo

      end subroutine set_diags_uvdz_ever  !]

! ----------------------------------------------------------------------
      subroutine set_diags_tdz_ever( istrR, iendR, jstrR, jendR, td, itrc )  ![
      ! need current dz*t to compare gross changes in t with sum of averages
      ! used only with 'check_avg' logical

      implicit none

      ! inputs
      integer, intent(in) :: istrR, iendR, jstrR, jendR, td, itrc

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            Td_a_sum(i,j,k) = t(i,j,k,nnew,itrc) * Hz(i,j,k)  ! re-use Td_a_sum instead of having Tdz_ever array
     &                       - tdz_t0(i,j,k,td)
          enddo
        enddo
      enddo

      end subroutine set_diags_tdz_ever  !]

! ----------------------------------------------------------------------
      subroutine wrt_tracer_diags_tile(istr,iend,jstr,jend)  ![
      ! since the value of the tracer has been stored as dz*t, we still need to
      ! compute the differences for the contribution of each diagnostic term

      ! Unlike surf_flux.F -> wrt_surface_flux_his, here we still need to calculate
      ! averages, and so the 'if (mod(tn,rec_rate)' is split in two parts.
      ! Also, since the value of 't' is kept in Tdiag and not the change in 't',
      ! we use Td_tmp to store the change, and for memory management it is only
      ! one array, and thus needs to be written to file immediately, or added to averges.
      ! In wrt_diagnostics_tile the if(mod(tn,rec_rate)) is split where it is so that the file can
      ! be opened and 'ocean_time' written before writing averages if applicable.

      ! Since the units of 't' change from 'Conc.' to 'dz*Conc.', the diagnostics
      ! can't be calculated in one loop, but must be done individually to adjust for
      ! the correct units of 't'.

      ! We divide by / dt to get diagnostic rate per second.

      use tracers, only: t
      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer             :: itrc, td, i, j, k, ierr
      character(len=30)   :: text

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      do td=1,nd_trc                                          ! only selected tracers from wrt_t_diag array
        itrc=diag_td_2_NT(td)                                 ! set the actual tracer index in 't'

        if (check_diag) T_dsum(istrR:iendR,jstrR:jendR,:)=0   ! could just do boundaries as the rest can be overwritten in first below
        Td_tmp(istrR:iendR,jstrR:jendR,:)=0                   ! do more efficient reset of boundaries here. Above changes bry values but below doesn't use them but can be seen in output if not set to 0 again.

        ! h_disp:
        do k=1,N                                              ! first diagnostic calculated using previous step's value
          do j=jstr,jend
            do i=istr,iend
              ! DISP = UP3 - ADV_4th:
              Td_tmp(i,j,k)=( Tdiag(i,j,k,td,it_h_disp) - t(i,j,k,nstp,itrc)*Hz_old(i,j,k) )/dt
     &                      - Tdiag(i,j,k,td,it_h_adv)
              if (check_diag)
     &        T_dsum(i,j,k)=Td_tmp(i,j,k)+Tdiag(i,j,k,td,it_h_adv)  ! add to sum for budget check
            enddo                                                   ! it_h_adv already calculated so add directly (share loop)
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(Td_tmp, istr, iend, jstr, jend, td, it_h_disp, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_h_disp)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! h_adv: (already calculated in set_diags_t_4th_adv)
        if (diag_avg) then
          call set_diags_t_avg_tile(Tdiag(:,:,:,td,it_h_adv), istr, iend, jstr, jend, td, it_h_adv, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_h_adv)
          call nc_write_var(ncid, Tdiag(:,:,:,td,it_h_adv), N, text, rp_var, output_rec, ierr)
        endif

        ! v_adv:
        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              Td_tmp(i,j,k)= ( Tdiag(i,j,k,td,it_v_adv) - Tdiag(i,j,k,td,it_h_disp) )/dt  ! UP3 stored in it_h_disp
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(Td_tmp, istr, iend, jstr, jend, td, it_v_adv, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_v_adv)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif


        ! From here onwards 't' was no longer in dz*t units:

        ! h_mix:
        do k=1,N
          do j=jstr,jend
            do i=istr,iend  ! both in dz*t units and Hz is at the same timestep
              Td_tmp(i,j,k)= ( Tdiag(i,j,k,td,it_h_mix)-Tdiag(i,j,k,td,it_v_mix) )*Hz(i,j,k) /dt
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(Td_tmp, istr, iend, jstr, jend, td, it_h_mix, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_h_mix)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! v_mix and boundary: (swapped order (vmix before hmix in step3d_t) so don't need to reset boundaries of Td_tmp)
        do k=1,N
          do j=jstr,jend
            do i=istr,iend  ! it_v_mix in 't' units, but it_prev in dz*t units
              Td_tmp(i,j,k)= ( Tdiag(i,j,k,td,it_v_mix)*Hz(i,j,k) - Tdiag(i,j,k,td,it_v_adv) )/dt
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (WESTERN_EDGE) then                      ! v_mix array also includes boundary update
          i=istr-1                                  ! for memory savings. The boundary term is
          Td_tmp(i,:,:) = Tdiag(i,:,:,td,it_v_mix)  ! already the difference, so can add it
          T_dsum(i,:,:) = Tdiag(i,:,:,td,it_v_mix)  ! directly to Td_tmp. Only boundary term
        endif                                       ! contributes to sum at the boundary.
        if (EASTERN_EDGE) then
          i=iend+1
          Td_tmp(i,:,:) = Tdiag(i,:,:,td,it_v_mix)
          T_dsum(i,:,:) = Tdiag(i,:,:,td,it_v_mix)
        endif
        if (SOUTHERN_EDGE) then
          j=jstr-1
          Td_tmp(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
          T_dsum(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
        endif
        if (NORTHERN_EDGE) then
          j=jend+1
          Td_tmp(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
          T_dsum(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
        endif
        if (diag_avg) then
          call set_diags_t_avg_tile(Td_tmp, istr, iend, jstr, jend, td, it_v_mix, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_v_mix)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! need rec_rate check here as not an average
        if ( (check_diag .or. tot_change) .and. mod(tn,rec_rate) == 0) then
          do k=1,N                    ! write gross difference: don't need to set Td_sum to zero as covering full loop ranges.
            do j=jstrR,jendR          ! first diagnostic calculated using previous step's value
              do i=istrR,iendR        ! Full domain loop range incl. bry.
                Td_tmp(i,j,k)= ( t(i,j,k,nnew,itrc)*Hz(i,j,k) - t(i,j,k,nstp,itrc)*Hz_old(i,j,k) )/dt
              enddo
            enddo
          enddo
          text=trim(t_vname(itrc)) / / check_vnm(1)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        if (check_diag .and. mod(tn,rec_rate) == 0) then                       ! need rec_rate check here as not an average
          text=trim(t_vname(itrc)) / / check_vnm(2)                            ! following to compare gross change in u/v and sum of diags.
          call nc_write_var(ncid, T_dsum, N, text, rp_var, output_rec, ierr)

          Td_tmp(istrR:iendR,jstrR:jendR,:)=Td_tmp(istrR:iendR,jstrR:jendR,:)
     &                                     -T_dsum(istrR:iendR,jstrR:jendR,:)  ! difference between gross change and sum of diags.
          text=trim(t_vname(itrc)) / / check_vnm(3)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif  ! <-- check_diag
      enddo    ! <-- itrc=1,NT

      Hz_old(istrR:iendR,jstrR:jendR,:) = Hz(istrR:iendR,jstrR:jendR,:)        ! update Hz_old

      end subroutine wrt_tracer_diags_tile  !]

! ----------------------------------------------------------------------
      subroutine wrt_diagnostics( tile )  ![
      implicit none

      ! input
      integer, intent(in) :: tile

# include "compute_tile_bounds.h"

      call wrt_diagnostics_tile(istr, iend, jstr, jend)

      end subroutine wrt_diagnostics  !]

! ----------------------------------------------------------------------
      subroutine wrt_diagnostics_tile(istr, iend, jstr, jend)  ![
      ! Main subroutine to write online diagnostics
      ! Unlike wrt_surface_flux_his, no need for write_file since
      ! if you want diagnostics you will want a results file to see them.

      implicit none

      ! inputs
      integer, intent(in)     :: istr, iend, jstr, jend

      ! local
      integer                 :: r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)   ! 2D grid types
#ifdef SOLVE3D
     &                         , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)  ! 3D grid types
#endif
      integer                 :: ierr = 0, dim, diag, i, j, k, tile, idt, itrc, td

      character(len=64)       :: fname
      character(len=30)       :: text, utext
      integer, dimension(4,2) :: grdtyp3d                      ! u and v type to use in loop
      integer, dimension(2)   :: gt_wrt = (/ up_var, vp_var /) ! grid type for writing

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      tn = iic-ntstart                                ! model step, as timestep t=n. (icc = t(n+1), and ntstart=1 if not from restart). We write the previous step results during the start of the next step

                                                      ! Trigger calculation of diags if required:
      if (diag_avg .or. mod(tn+1,rec_rate) == 0) then ! +1 as calculate in the step before writing. so if tn=5 then calc_diag at tn=4, hence need +1 to trigger mod(tn...)
        calc_diag = .true.                            ! This triggers diags to be calculated in equation system. Always on for averages
      else                                            ! set here once rather than everywhere in equation system redo 'mod'
        calc_diag = .false.                           ! since this is called at the beginning of main.
      endif
      if (diag_avg .or. mod(tn+2,rec_rate) == 0) then ! +2 as need to set the previous values before calc_diag step above
        set_dz_old=.true.                             ! So if tn=5 then calc_diag at tn=4, and need set_dz_old tn=3 hence +2 to trigger mod(tn...)
      else
        set_dz_old=.false.
      endif

      if (FIRST_TIME_STEP) call init_diags(istr, jstr, iend, jend)  ! called here as need set_dz_old set above

      diag_icori = .false.  ! flag to avoid computing 'icori' diagnostic twice in compute_horiz_uv_terms.h. Won't compute during pre_step3d but will for step3d_uv1

      if (iic>ntstart) then ! catch to prevent output of zeros at inital step since mod(0,anything) = 0...

! Defining and writing netcdf variables only done by master thread:
C$OMP MASTER

        if ( mod(tn,rec_rate) == 0 ) then

          if(total_output_rec == 0) then                            ! Add root output name to file name (only done once)
            file_name = trim(output_root_name) / / trim(file_name)  ! Notice below space between '/ /', this is needed to avoid
          endif                                                     ! cpp preprocessor seeing // as a c++ comment and deleting it.


          if(output_rec == 0) then                                       ! if file doesn't yet exist or need new file!

            call ncdf_create_file(file_name, ncid, prev_fill_mode        ! lose the prev_fill_mode? r2dgrd, ...., auxil are id's to dimensions, needed to create vars in the file
     &                    ,recs_per_file, total_output_rec, auxil
#ifdef SOLVE3D
!     &                    ,r3dgrd, u3dgrd, v3dgrd, w3dgrd               ! 3D grid dimensions (usually - not diagnostics)
     &                    ,r3dgrd, grdtyp3d(:,1), grdtyp3d(:,2), w3dgrd  ! 3D grid dimensions
#endif
     &                    ,r2dgrd, u2dgrd, v2dgrd )                      ! 2D grid dimensions


            ! Define remaining variables:

            if (diag_uv) then
              do dim=1,2                                  ! u,v dims

                do diag=1,nmd                             ! diag terms
                  text  = uvw(dim) / / vnames(diag)       ! diag short name
                  call nc_define_var(ncid, text, lname(diag), units, grdtyp3d(:,dim), ierr )
                enddo

                if (check_diag .or. tot_change) then      ! following to compare gross change in u/v and sum of diags.
                  text  = uvw(dim) / / check_vnm(1)       ! gross change in dz*u or dz*v over 1 timestep
                  call nc_define_var(ncid, text, check_lnm(1), units, grdtyp3d(:,dim), ierr )
                endif

                if (check_diag) then                      ! following to compare gross change in u/v and sum of diags.
                  text  = uvw(dim) / / check_vnm(2)       ! sum of all diag terms over 1 timestep
                  call nc_define_var(ncid, text, check_lnm(2), units, grdtyp3d(:,dim), ierr )
                  text  = uvw(dim) / / check_vnm(3)       ! compare gross change against sum of diags. Should be approx zero.
                  call nc_define_var(ncid, text, check_lnm(3), units, grdtyp3d(:,dim), ierr )
                endif

                if (check_avg) then
                  text  = uvw(dim) / / check_avnm(1)      ! gross change in u/v up until now since t0
                  call nc_define_var(ncid, text, check_alnm(1), units, grdtyp3d(:,dim), ierr )
                  text  = uvw(dim) / / check_avnm(2)      ! total diag avg in u/v up until now since t0
                  call nc_define_var(ncid, text, check_alnm(2), units, grdtyp3d(:,dim), ierr )
                  text  = uvw(dim) / / check_avnm(3)      ! total budget close until now since t0
                  call nc_define_var(ncid, text, check_alnm(3), units, grdtyp3d(:,dim), ierr )
                endif

              enddo ! <-- dim
            endif   ! <-- diag_uv

            if (diag_trc) then
              do td=1,nd_trc                                    ! only tracers of interest
                itrc=diag_td_2_NT(td)                           ! set the actual tracer index in 't'
                utext = trim(t_units(itrc)) / / tunits
                do idt=1,ndt                                    ! tracer diag terms
                  text  = trim(t_vname(itrc)) / / tdname(idt)
                  call nc_define_var(ncid, text, tlname(idt), utext, r3dgrd, ierr )
                enddo

                if (check_diag .or. tot_change) then
                  text  = trim(t_vname(itrc)) / / check_vnm(1)  ! gross change in dz*t over 1 timestep
                  call nc_define_var(ncid, text,  check_lnm(1), utext, r3dgrd, ierr )
                endif

                if (check_diag) then
                  text  = trim(t_vname(itrc)) / / check_vnm(2)  ! sum of all diag terms over 1 timestep
                  call nc_define_var(ncid, text,  check_lnm(2), utext, r3dgrd, ierr )
                  text  = trim(t_vname(itrc)) / / check_vnm(3)  ! compare gross change against sum of diags. Should be approx zero.
                  call nc_define_var(ncid, text,  check_lnm(3), utext, r3dgrd, ierr )
                endif

                if (check_avg) then
                  text  = trim(t_vname(itrc)) / / check_avnm(1) ! gross change in t up until now since t0
                  call nc_define_var(ncid, text, check_alnm(1), utext, r3dgrd, ierr )
                  text  = trim(t_vname(itrc)) / / check_avnm(2) ! total diag avg in t up until now since t0
                  call nc_define_var(ncid, text, check_alnm(2), utext, r3dgrd, ierr )
                  text  = trim(t_vname(itrc)) / / check_avnm(3) ! total budget close until now since t0
                  call nc_define_var(ncid, text, check_alnm(3), utext, r3dgrd, ierr )
                endif

              enddo  ! <-- td
            endif    ! <-- diag_trc

            ierr=nf90_enddef(ncid) ! end definition to write changes to disk
            if (mynode == 0) then
              write(*,'(6x,2A)')
     &          'diagnostics :: created new netCDF file ', file_name
            endif

          endif ! <-- (output_rec == 0; needed to create a file)

          ! WRITE VARIABLES TO FILE
          ! =======================

          ierr=nf90_open(file_name, nf90_write, ncid)

          ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)   ! set fill value - nf90_nofill for optimized writing

                                                                  ! advance output record number
          output_rec = output_rec + 1                             ! current file record
          total_output_rec = total_output_rec + 1                 ! total number of recs written


          call nc_write_time(ncid, output_rec, total_output_rec)  ! write ocean time

        endif   ! <-- mod(tn,rec_rate) == 0

C$OMP END MASTER
C$OMP BARRIER

        if (diag_avg) then                                                 ! see wrt_tracer_diags_tile for description of why if(mod(tn,rec_rate)) is split here.
          if (diag_uv)  call set_diags_avg_tile   (istr, iend, jstr, jend) ! avoid very first step because diags are only all set by the end of the
          if (diag_trc) call wrt_tracer_diags_tile(istr, iend, jstr, jend) ! first step, and then added to averages and/or output at the beginning of the next step.
        endif                                                              ! Need this after defined file as wrt_tracer will output averages so needs file created.

! Defining and writing netcdf variables only done by master thread:
C$OMP MASTER

        if ( mod(tn,rec_rate) == 0 ) then

          if (diag_uv) then
            do dim=1,2                            ! over dims u/v

              do diag=1,nmd                       ! over diag terms
                text = uvw(dim) / / vnames(diag)
                if (diag_avg) then
                  call nc_write_var(ncid, Md_avg(:,:,:,dim,diag), N, text, gt_wrt(dim), output_rec, ierr)
                else
                  call nc_write_var(ncid,  Mdiag(:,:,:,dim,diag), N, text, gt_wrt(dim), output_rec, ierr)
                endif
              enddo

              if (check_diag .or. tot_change) then
                ! do uv_dif calculation at end of step3d_uv2.F as have istr, iend there already.
                ! but use 'nnew' there as still in previous timestep!
                text = uvw(dim) / / check_vnm(1)                           ! gross change in dz*u or dz*v over 1 timestep
                call nc_write_var(ncid, uv_dif(:,:,:,dim), N, text, gt_wrt(dim), output_rec, ierr)
              endif

              if (check_diag) then

                Mdsum(istrR:iendR,jstrR:jendR,:) = 0                       ! refresh for this step
                do diag=1,nmd
                  do j=jstrR,jendR
                    do i=istrR,iendR
                      Mdsum(i,j,:) = Mdsum(i,j,:) + Mdiag(i,j,:,dim,diag)  ! should make the ranges more efficient (e.g. bry term don't need interior points)
                    enddo
                  enddo
                enddo
                text = uvw(dim) / / check_vnm(2)                           ! sum of diag terms over 1 timestep
                call nc_write_var(ncid, Mdsum, N, text, gt_wrt(dim), output_rec, ierr)

                do j=jstrR,jendR
                    do i=istrR,iendR
                      Mdsum(i,j,:) = uv_dif(i,j,:,dim)-Mdsum(i,j,:)        ! difference between gross dz*u (or v) change and sum of diags of dz*u
                  enddo
                enddo
                text = uvw(dim) / / check_vnm(3)
                call nc_write_var(ncid, Mdsum, N, text, gt_wrt(dim), output_rec, ierr)
              endif

              if (check_avg) then
                text  = uvw(dim) / / check_avnm(1)                         ! gross change in u/v up until now since t0
                call nc_write_var(ncid, uvdz_ever(:,:,:,dim), N, text, gt_wrt(dim), output_rec, ierr)
                text  = uvw(dim) / / check_avnm(2)                         ! total diag avg in u/v up until now since t0
                call nc_write_var(ncid, Md_a_ever(:,:,:,dim), N, text, gt_wrt(dim), output_rec, ierr)
                do j=jstrR,jendR; do i=istrR,iendR
                  Mdsum(i,j,:) = uvdz_ever(i,j,:,dim)-Md_a_ever(i,j,:,dim) ! repurpose Mdsum array
                enddo; enddo
                text  = uvw(dim) / / check_avnm(3)
                call nc_write_var(ncid,     Mdsum(:,:,:),     N, text, gt_wrt(dim), output_rec, ierr)
              endif

            enddo  ! <-- dim
          endif    ! <-- diag_uv

          if (diag_trc .and. .not. diag_avg)                             ! for averages this has been called already
     &      call wrt_tracer_diags_tile(istr, iend, jstr, jend)           ! we call it here non-averages to get the correct writing frequency

          ierr=nf90_close (ncid)                                         ! close netcdf file

          if (mynode == 0) then
            write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &       'diagnostics :: wrote output, tdays =', tdays,
     &       'step =', iic-1, 'rec =', output_rec, '/', total_output_rec
     &        MYID
          endif

          if (output_rec == recs_per_file) then                          ! reset output_rec to create new file
            output_rec = 0
          endif

        endif ! <-- mod(tn,rec_rate)

C$OMP END MASTER
C$OMP BARRIER

      endif   ! <-- icc>ntstart

      end subroutine wrt_diagnostics_tile  !]

! ----------------------------------------------------------------------

#endif /* DIAGNOSTICS */

      end module diagnostics
