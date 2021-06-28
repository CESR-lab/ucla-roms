      module diagnostics

      ! INFO: ![
      ! Online diagnostics for ROMS simulation.
      ! Initial coding DevinD and JeroenM - 2021/05

      ! To compute diagnostics use DIAGNOSTICS flag.

      ! HIS vs AVG:
      ! To chose between history slices of diagnostics or averages,
      ! set diag_avg to .false. for history, or .true. for averaging.
      ! Averaging is much less efficient as diagnostics must be collected
      ! at everytime step to maintain a true average.

      ! The main subroutine is wrt_diagnostics_tile, which is called at the very
      ! beginning of main.F in order to record the diagnostics of 'nnew' index,
      ! which was t=n+1 for the previous timestep, and is now t=n since roms_step
      ! has moved onto the next timestep. I.e. 'nnew' will become 'nstp' in next timestep.
      ! But wrt_online_diags is called before nnew and nstp indicies change in main.F.

      ! If you are needing averages then the diagnotics need to be computed at everytimestep.
      ! If you just need snapshots then the diagnostics only need to be calculated at
      ! your output frequency.
      ! Both of these are controlled by the logical 'calc_diag', which is set in
      ! wrt_diagnostics_tile at the beginning of a timestep, and then that logical is
      ! used in the equation system to trigger a calculation of diagnostic.

      ! TRACERS:
      ! Since there are 6 diagnostic terms for each tracer, it becomes very memory
      ! intensive to keep diag arrays for every tracer. Instead, the user selects
      ! only tracers of interest by setting T/F in wrt_t_diags in tracers_defs.h.
      ! Then the number of selected tracers is counted in init_diags_trc and an array
      ! of the correct size is allocated there. Some book keeping variables are needed
      ! to convert indices between t(itrc) the tracer variable in roms for all tracers,
      ! and Tdiag(td) in diagnostics which only has indices for selected tracers.

      ! UNITS:
      ! The roms history file give u/v in m/s, and 't' tracers in their
      ! concentrations. However, since the depth of the water column is also evolving
      ! during a timestep, it is easier to work in dz*u/v/w/t, since that is the unit
      ! of u/v/w/t when the depth is changed in the timestep as follows:
      ! step3d_uv1.F (depth not changed just to show position amongst code)
      ! step2d_tile -> set_depth_tile -> (depth (Hz) changes)
      ! step3d_uv2.F (depth not changed just to show position amongst code)
      !
      ! See Documentation/readme-specifics/README.diagnostics_u_v_value_changes
      ! for more details on how the units of 'u' change during a timestep.

      ! LOGICAL PARAMETERS:
      ! It is assumed that by using logicals that are set as parameters
      ! the 'if branch' will be removed during compiling, since the outcome
      ! of the branch is always known and unchanging. If this is true it will allow
      ! for if statements to be placed within loops, which is usually not advised
      ! due to the cost of the if statment at every index iteration.

      ! We should move some parameters into .in file, as it is way more convenient,
      ! but a template of what should be in the .in file when using this module should be
      ! left at the top here, so that it can be copied into the .in, you don't need
      ! to find an old diagnostics .in file...

      ! U/V evolution:
      ! During a roms_step in main.F, the values of u/v at t=n+1/2 is calculated
      ! in pre_step3d4s.F. Since we are interested in the change from t=n to t=n+1,
      ! we do not track and diagnostics in pre_step3d.
      ! The order of evolution of u/v (showns as u) and its units is as follows:
      ! 1) dz*u - Hydrostatic pressure gradient - prsgrd.F
      ! 2) dz*u - Coriolis & curvilinear coords - compute_horiz_rhs_uv_term.h
      ! 3) dz*u - Advection                     - step3d_uv1.F & compute_horiz/vert_rhs_uv_term.h
      ! 4) dz*u - Viscous terms                 - visc3d_GP.F
      ! -         Now depth (Hz) at t=n+1 in calculated in step2d_FB -> set_depth_tile.
      ! 5)    u - Vertical mix & 2D/3D coupling - step3d_uv2.F part 1
      ! 6)    u - Boundary changes              - u3dbc_im.F
      ! 7)    u - 2D-3D correction              - step3d_uv2.F part 2

      ! Clearly, care needs to be taken when comparing 'u' at each step, since the units
      ! after a step can change.

      ! TRACER evolution:
      ! 1) dz*C - Horizontal advective fluxes                               - compute_horiz_tracer_fluxes.h
      ! 2) dz*C - Vertical advective fluxes                                 - compute_vert_tracer_fluxes.h
      ! -         Surface flux is done here
      ! 3) dz*C - Nonlocal transport flux for unstable (convective) forcing - step3d_t_ISO.F (lmd_kpp)
      ! 4)    C - Implicit time step for vertical diffusion                 - step3d_t_ISO.F
      ! 5)    C - Boundary changes                                          - t3dbc_im.F
      ! 6)    C - Horizontal Laplacian diffusion terms                      - t3dmix_GP.F

      ! W evolution:
      ! 1)    ? -                                       - prsgrd.F
      ! 2)    ? -                                       - compute_horiz_rhs_w_terms.h
      ! 3)    ? -                                       - compute_vert_rhs_w_terms.h'
      ! 4)    ? -                                       - step3d_uv2.F
      ! 5)    ? -                                       - w3dbc_im.F

      ! QUESTIONS:
      ! - leave in functionality to confirm budgets.
      ! - how to store bry term. and FC terms (does little in middle)
      ! - might we not want uv diagnostics and just tracers? I.e. a logical for that?

      ! MISC:
      ! It is my feeling that once we use averages, all diagnostic calculations should
      ! not be calculated in the step3d... loop, but in it's own subroutine to avoid
      ! having if statments in loops. This will regardless result in the cleanest code.

      ! Mdiag(2): hydrostatic pressure gradient  ! OLD JM NOTES
      ! Mdiag(3): cor(u)
      ! Mdiag(4): adv(u) (non-dissipative part)
      ! Mdiag(5): horizontal mixing (including numerical dissipation)
      ! Mdiag(6): vertical mixing   (including numerical dissipation)
      ! Mdiag(7): non-hydrostatic pressure gradient

      ! Tdiag(1): dt/dt  ! OLD JM NOTES
      ! Tdiag(2): adv(t) (non-dissipative part)
      ! Tdiag(3): adv(t) (dissipative part)
      ! Tdiag(4): vertical mixing
      ! Tdiag(5): sources/sinks
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
      ! User inputs
      logical, parameter         :: write_file    = .true.  ! whether to write a file or not (probably remove this as diagnostics not helpful without file)
      logical, parameter         :: diag_avg      = .true.  ! compute history (F) or averages (T)
      integer                    :: rec_rate      = 2       ! output record rate, move this to .in file rather?
      integer                    :: recs_per_file = 20      ! total recs per file

      logical, parameter, public :: diag_uv       = .true.  ! calc+write diagnostics of u/v. Replaces cppflag.
      logical, parameter, public :: diag_trc      = .true.  ! calc+write diagnostics of selected tracers. Replaces CPP key
      logical, parameter         :: tot_change    = .true.  ! Gross (total) change in variable over timestep. May not want this on since it's just the sum of diagnostic terms.
      ! End of user inputs
      ! ***************************************************************

      ! Developer inputs:
      ! - need to have on diag_uv to use these flags:
      logical, parameter         :: check_diag    = .true.  ! confirm sum of diagnostic terms    = gross change in u/v per timestep. Useful to ensure diagnostics not broken by different cppdefs.h config.
      logical, parameter         :: check_avg     = .true.  ! confirm sum of diagnostic avg*time = gross change in u/v over entire simulation at time of writing

      ! Preamble:  ![

      logical, public :: calc_diag  ! flag to calculate diagnostics in equation system.
      logical, public :: set_dz_old ! set dz*u/v/w/t at previous timestep, needed for diag calculation when skipping steps if not averages.

      ! U & V MOMENTUM:
      ! ---------------

      ! NOTE: u & v are 1:N in z, but w is 0:N hence different variables needed

      integer, public, parameter                          :: nmd = 7    ! number of momentum diag terms

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

      real, public, dimension(GLOBAL_1DY_ARRAY, N, 2)     :: uv_wbry    ! bry change for u and v

      integer, parameter, public ::
     &  iprsgr = 1, ! internal loop ranges only (i.e. excluding domain boundary)
     &  icori  = 2, ! internal loop ranges only
     &  iadv   = 3, ! internal loop ranges only
     &  ivisc  = 4, ! internal loop ranges only
     &  iuv2DC = 5, ! internal loop ranges only
     &  iu3d   = 6, ! bry only
     &  iuv2FC = 7  ! internal + bry (full domain) loop ranges

      ! move this to init_diags
      character (len=20), dimension(nmd) :: vnames = (/
     &  '_prsgrd',
     &  '_coriolis',
     &  '_advection',
     &  '_visc3d',
     &  '_vert_mix_2D3D',
     &  '_3dbc',       ! should store 3dbc changes in 2D arrays not 3D, as just 0 in interior
     &  '_2D3D_cor'
     &                                                /)

      character (len=70), dimension(nmd) :: lname =   (/               ! len=70 allows for 8 characters of 'Average ' if diag_avg
     &  'Hydrostatic pressure gradient - prsgrd.F',
     &  'Coriolis & curvilinear coords - compute_horiz_rhs_uv_term.h',
     &  'Advection - step3d_uv1.F',
     &  'Viscous terms - visc3d_GP.F',
     &  'Vertical mixing & 2D/3D coupling - step3d_uv2.F part 1',
     &  'Boundary changes - u/v3dbc_im.F',
     &  '2D-3D correction - step3d_uv2.F part 2'
     &                                                 /)

      character (len=20)               :: units = 'm^2/s (dz*(u or v))'
      character (len=1),  dimension(3) :: uvw = (/ 'u', 'v', 'w' /)      ! prefix for u/v/w outputs

      character (len=80), dimension(:), allocatable :: check_lnm         ! long name for outputs checking diagnostics
      character (len=20), dimension(:), allocatable :: check_vnm         ! short name ...
      character (len=80), dimension(:), allocatable :: check_alnm        ! long name for outputs checking diagnostics averages
      character (len=20), dimension(:), allocatable :: check_avnm        ! short name ...

      logical, public :: diag_icori  ! flag to avoid computing 'icori' diagnostic twice in compute_horiz_uv_terms.h

      ! TRACERS:
      ! --------
# ifdef DIAGS_TRC
      ! Since there could be many tracers, it is preferable not to carry any 't_old'
      ! or temporary arrays for each tracer to avoid the memory burden. This differs
      ! from the u,v,w implementation here currently.
      ! Need to consider pipe and river sources.

      integer, public, parameter                          :: ndt = 5   ! number of tracer diag terms
#  ifdef LMD_KPP
     &                                                           + 1
#  endif
      real, public, allocatable, dimension(:, :, :, :, :) :: Tdiag     ! tracer diagnostic terms
      real,         allocatable, dimension(:, :, :)       :: Td_tmp    ! temporary array to calculate change in tracer diags
      real,         allocatable, dimension(:, :, :)       :: Hz_old    ! store previous Hz to calculate old dz*t. better than storing dz*t for each tracer. Better to recalculate it?

      ! if diag_avg
      real,         allocatable, dimension(:, :, :, :, :) :: Td_avg    ! averages of tracer diagnostic terms

      ! if check_diag
      real,         allocatable, dimension(:, :, :)       :: T_dsum    ! sum of t diagnostics to compare with gross change

      ! if check_avg
      real, allocatable, dimension(:, :, :, :)            :: Tdz_t0    ! original dz*t to compare total change in t over simulation.
      real, allocatable, dimension(:, :, :)               :: Td_a_sum  ! sum of diag averages at timestep
      real, allocatable, dimension(:, :, :, :)            :: Td_a_ever ! total of all diag averages for entire simulation. Should add up to change in t over simulation...

      integer                            :: nd_trc           ! number of tracers selected tracers of interest (from wrt_t_diag array)
      integer, dimension(:), allocatable :: diag_td_2_NT     ! index of tracer of selected tracers only
      integer, dimension(NT), public     :: diag_NT_2_td=-1  ! convert itrc=1,NT to index of tracer in Tdiag array. Set to bad value -1

      integer, parameter, public ::
     &  it_h_adv    = 1,  ! compute_horiz_tracer_fluxes.h
     &  it_v_adv    = 2,  ! compute_vert_tracer_fluxes.h
      ! Do we need surface forcing flux as a vertical? Shouldn't be a 3D array anyway
     &  it_imp_vdif = 3,  ! step3d_t_ISO.F
     &  it_t3dmix   = 4,
     &  it_lmdkpp   = 5,  ! not used if no LMD_KPP flag, hence at the end. ! step3d_t_ISO.F. This only change T & S so shouldn't be 'NT' in size, a waste.
     &  it_bc       = 6   ! BC's shouldn't have 3D array

      ! move this to init_diags_trc
      character (len=30), dimension(ndt) :: tdname = (/
     &  '_hori_adv',
     &  '_vert_adv',
     &  '_imp_vdif',
     &  '_t3dmix',
     &  '_lmdkpp',
     &  '_t3dbc'
     &                                               /)

      character (len=90), dimension(ndt) :: tlname = (/
     &  'Horizontal advective fluxes - compute_horiz_tracer_fluxes.h',
     &  'Vertical advective fluxes - compute_vert_tracer_fluxes.h',
     &  'Implicit time step for vertical diffusion - step3d_t_ISO.F',
     &  'Horizontal (along geopotential surfaces) Laplacian diffusion terms - t3dmix_GP.F',
     &  'Nonlocal transport flux for unstable (convective) forcing - step3d_t_ISO.F',
     &  'Boundary changes - t3dbc_im.F'
     &                                               /)

      character (len=20) :: tunits = ') * m   (dC*dz)'  ! not sure on units here.

      public init_diags_trc
      public set_diags_t_at_bc
      public set_diags_t_at_t3dmix_GP
# endif /* DIAGS_TRC */

      ! Netcdf outputting:
      ! ------------------
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      integer            :: tn                          ! roms timestep
      character (len=60) :: file_name = '_on_diags.nc'  ! suffix to append to simulation output rootname
      integer            :: ncid=-1, prev_fill_mode


      ! W MOMENTUM:
      ! -----------
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

      public set_diags_at_visc3d
      public set_diags_at_uv1
      public set_diags_at_uv2_DC
      public set_diags_at_bc_u
      public set_diags_at_bc_v
      public set_diags_at_uv2_FCu
      public set_diags_at_uv2_FCv
      public set_diags_dz_old
      public wrt_diagnostics !]

      contains

! ----------------------------------------------------------------------
      subroutine init_diags(istr, jstr, iend, jend)  ![
      ! Set original dz*u/v/w/t for initial timestep, hence only called once.
      ! Called from prestep before Hz has changed in the timestep.

      implicit none

      ! inputs
      integer, intent(in) :: istr, jstr, iend, jend

      ! local
      integer :: i, j, k, diag

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      if (diag_uv) then
        allocate(    Mdiag(GLOBAL_2D_ARRAY, N, 2, nmd) )
        allocate( uvdz_old(GLOBAL_2D_ARRAY, N, 2)      )
        allocate(  uv_prev(GLOBAL_2D_ARRAY, N, 2)      )
        Mdiag    (istrR:iendR,jstrR:jendR,:,:,:) = 0
        uvdz_old (istrR:iendR,jstrR:jendR,:,:)   = 0
        uv_prev  (istrR:iendR,jstrR:jendR,:,:)   = 0
      endif

      if (set_dz_old) then

        if (diag_uv) then
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

        if(mynode==0) print *, 'INIT HZ_OLD **************'  ! debug

# ifdef NHMG
        do k=1,N ! 0 is always 0?
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
# endif /* NHMG */

      endif  ! <-- set_dz_old

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

        allocate( check_lnm(3), check_vnm(3)    )
        check_vnm(1) = '_1step'
        check_lnm(1) = 'Total change in variable over timestep.'
        check_vnm(2) = '_d_sum'
        check_lnm(2) = 'Sum of diagnostics terms over timestep.'
        check_vnm(3) = '_check'
        check_lnm(3) = 'Diff. of total change in variable ' / /
     &                 'vs sum of diags (must match) over timestep.'
      endif

      if (check_diag) then
        allocate(  Mdsum(GLOBAL_2D_ARRAY, N)    )
        Mdsum (istrR:iendR,jstrR:jendR,:)   = 0
      endif

      if (check_avg) then
        allocate(  Md_a_sum(GLOBAL_2D_ARRAY, N, 2) )         ! sum of diag averages at timestep
        allocate( Md_a_ever(GLOBAL_2D_ARRAY, N, 2) )         ! continuing total of sum of averages Md_a_ever += Md_a_sum
        allocate(   uvdz_t0(GLOBAL_2D_ARRAY, N, 2) )
        allocate( uvdz_ever(GLOBAL_2D_ARRAY, N, 2) )
        Md_a_ever(istrR:iendR,jstrR:jendR,:,:) = 0           ! Md_a_sum set to 0 elsewhere, so not needed here
        uvdz_t0 = uvdz_old                                   ! store the initial condition for budget confirmation at simulation end

        allocate( check_alnm(2), check_avnm(2)     )
        check_avnm(1) = '_dz_ever'
        check_alnm(1) = 'Total change in variable since start of run.'
        check_avnm(2) = '_avg_ever'
        check_alnm(2) = 'Sum of diagnostic averages multiplied by time'
     &                  / / ' since start of run.'
      endif

      if (diag_trc) call init_diags_trc(istrR, jstrR, iendR, jendR)

      end subroutine init_diags  !]

! ----------------------------------------------------------------------
      subroutine init_diags_trc(istrR, jstrR, iendR, jendR)  ![
      ! setup tracer arrays only for desired tracers,
      ! since can be memory intensive if all tracers.

      implicit none

      ! inputs
      integer, intent(in) :: istrR, jstrR, iendR, jendR

      ! local
      integer                :: itrc, cnt=0, diag, td, i, j, k
      integer, dimension(NT) :: temp_t_indx



      do itrc=1,NT
        if (wrt_t_diag(itrc)) then
          cnt=cnt+1                 ! count tracers we want to calc+write
          diag_NT_2_td(itrc)=cnt    ! t(itrc) = Tdiag(cnt) - to convert between index scheme
        endif
      enddo

      nd_trc = cnt                                           ! number of tracers of interest
      allocate(  Tdiag( GLOBAL_2D_ARRAY, N, nd_trc, ndt ) )  ! only for tracers we are interested in
      allocate( Td_tmp( GLOBAL_2D_ARRAY, N ) )               ! don't need to set to zero as written over in code
      allocate( Hz_old( GLOBAL_2D_ARRAY, N ) )               ! same for Hz_old
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
          tlname(diag)  = 'Average ' / / trim(tlname(diag))      ! and prepend 'Average'
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

      if(mynode==0) print *, 'nd_trc=',nd_trc ! debug
      if(mynode==0) print *, 'diag_td_2_NT(1)=',diag_td_2_NT(1)

      end subroutine init_diags_trc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv1( istr, iend, jstr, jend, istrA, j, dim, DC, ruv, suvstr, uv )  ![
      ! Set the diagnostic terms for prsgrd, coriolis, and advection.
      ! Works for both u and v.
      ! Prsgrd and coriolis need to also be convert to dz * u/v units

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrA, j, dim         ! 'istrA' for any as depends on inputs, which 'dim' dimension.

      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N), intent(in) :: DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N),   intent(in) :: ruv     ! ru or rv
      real, dimension(GLOBAL_2D_ARRAY),              intent(in) :: suvstr  ! sustr or svstr
      real, dimension(GLOBAL_2D_ARRAY,N,3),          intent(in) :: uv      ! u or v

      ! local
      integer :: i, k

      ! NOTE: currently missing WEC term

      do i=istrA,iend
        Mdiag(i,j,N,dim,iprsgr) = Mdiag(i,j,N,dim,iprsgr) * DC(i,0)  ! need to convert to dz*u now that we have DC
        Mdiag(i,j,1,dim,iprsgr) = Mdiag(i,j,1,dim,iprsgr) * DC(i,0)  ! consider taking this out of loop and doing seperately
        Mdiag(i,j,N,dim,icori)  = Mdiag(i,j,N,dim,icori)  * DC(i,0)  ! perhaps even in own subroutine to avoid clutter
        Mdiag(i,j,1,dim,icori)  = Mdiag(i,j,1,dim,icori)  * DC(i,0)
        Mdiag(i,j,N,dim,iadv)   = DC(i,0)*ruv(i,j,N) + dt*suvstr(i,j)                 ! all advection - includes 2nd part of compute_horiz_uv & compute_vert_uv.
     &                           - ( Mdiag(i,j,N,dim,iprsgr) + Mdiag(i,j,N,dim,icori) ) ! subtract previous 2 terms included in total
        uv_prev(i,j,N,dim)      = uv(i,j,N,nnew)                     ! store for next diag visc3d
        Mdiag(i,j,1,dim,iadv)   = DC(i,0)*ruv(i,j,1)
     &                           - ( Mdiag(i,j,1,dim,iprsgr) + Mdiag(i,j,1,dim,icori) )
        uv_prev(i,j,1,dim)      = uv(i,j,1,nnew)
      enddo

      do k=2,N-1
        do i=istrA,iend
          Mdiag(i,j,k,dim,iprsgr) = Mdiag(i,j,k,dim,iprsgr) * DC(i,0)
          Mdiag(i,j,k,dim,icori)  = Mdiag(i,j,k,dim,icori)  * DC(i,0)
          uv_prev(i,j,k,dim)      = uv(i,j,k,nnew)
          Mdiag(i,j,k,dim,iadv)   = DC(i,0)*ruv(i,j,k)                 ! all advection - includes 2nd part of compute_horiz_uv & compute_vert_uv.
     &                            - ( Mdiag(i,j,k,dim,iprsgr) + Mdiag(i,j,k,dim,icori) )
        enddo
      enddo

      end subroutine set_diags_at_uv1  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_visc3d( istr, iend, jstr, jend, isA, jsA, dim, uv )  ![
      ! set diagnostics for u or v in visc3d_GP.F
      ! interior points only (i.e. excluding bry)

      implicit none

      ! inputs
      integer,                              intent(in) :: istr, iend, jstr, jend, isA, jsA, dim  ! isA jsA are starts depending on u or v variable
      real, dimension(GLOBAL_2D_ARRAY,N,3), intent(in) :: uv        ! u or v

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jsA,jend    ! interior of subdomains, exludes global domain boundaries
          do i=isA,iend
            Mdiag( i,j,k,dim,ivisc) = uv(i,j,k,nnew) - uv_prev(i,j,k,dim)  ! units already in dz*u
            uv_prev(i,j,k,dim)      = uv(i,j,k,nnew)
          enddo
        enddo
      enddo

      end subroutine set_diags_at_visc3d  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_DC( istr, iend, jstr, jend, istrU, jstrV )  ![
      ! set diagnostics for u or v in step3d_uv2.F
      ! interior points only (i.e. excluding bry)

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrU, jstrV

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstr,jend    ! interior of subdomains, exludes global domain boundaries
          do i=istrU,iend
            Mdtmp                  = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) ! convert to dz*u
            Mdiag( i,j,k,1,iuv2DC) = Mdtmp - uv_prev(i,j,k,1)                        ! avoids recomputing above
            uv_prev(i,j,k,1)       = Mdtmp
          enddo
        enddo
        do j=jstrV,jend    ! interior of subdomains, exludes global domain boundaries
          do i=istr,iend
            Mdtmp                  = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
            Mdiag( i,j,k,2,iuv2DC) = Mdtmp - uv_prev(i,j,k,2)
            uv_prev(i,j,k,2)       = Mdtmp
          enddo
        enddo
      enddo

      if(mynode==0) print *, 'uv2DC!!!!!!!!!'  ! debug only

      end subroutine set_diags_at_uv2_DC  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_bc_u( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in u
      ! called from step3d_uv2.F right after call to u3dbc_im.F

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i

      ! currently this is not perfectly setup for openMP since we have ':' instead
      ! of istr:iend, etc. However, since it is just a 2D operation at bry, it is ok? for now.

      if (WESTERN_EDGE)  uv_prev(istr  ,:,:,1) = u(istr  ,:,:,nnew) * 0.5*(Hz(istr-1,:,:)+Hz(istr  ,:,:))
      if (EASTERN_EDGE)  uv_prev(iend+1,:,:,1) = u(iend+1,:,:,nnew) * 0.5*(Hz(iend  ,:,:)+Hz(iend+1,:,:))

      if (WESTERN_EDGE)  uv_wbry(:,:,1)           = uv_prev(istr  ,:,:,1) - uvdz_old(istr  ,:,:, 1)

      if (WESTERN_EDGE)  Mdiag(istr  ,:,:,1,iu3d) = uv_prev(istr  ,:,:,1) - uvdz_old(istr  ,:,:, 1) ! uv_prev is already u_current, can do this because have uvdz_old already
      if (EASTERN_EDGE)  Mdiag(iend+1,:,:,1,iu3d) = uv_prev(iend+1,:,:,1) - uvdz_old(iend+1,:,:, 1)

      do i=istr,iend+1 ! use loop here else indexing would be even messier
        if (SOUTHERN_EDGE) uv_prev(i,jstr-1,:,1)=u(i,jstr-1,:,nnew)*0.5*(Hz(i-1,jstr-1,:)+Hz(i,jstr-1,:))
        if (NORTHERN_EDGE) uv_prev(i,jend+1,:,1)=u(i,jend+1,:,nnew)*0.5*(Hz(i-1,jend+1,:)+Hz(i,jend+1,:))
      enddo
      if (SOUTHERN_EDGE) Mdiag(:,jstr-1,:,1,iu3d) = uv_prev(:,jstr-1,:,1) - uvdz_old(:,jstr-1,:, 1)
      if (NORTHERN_EDGE) Mdiag(:,jend+1,:,1,iu3d) = uv_prev(:,jend+1,:,1) - uvdz_old(:,jend+1,:, 1)

      end subroutine set_diags_at_bc_u  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_bc_v( istr, iend, jstr, jend )  ![
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
      if (WESTERN_EDGE)  Mdiag(istr-1,:,:,2,iu3d) = uv_prev(istr-1,:,:,2) - uvdz_old(istr-1,:,:,2) ! uv_prev is already u_current, can do this because have uvdz_old already
      if (EASTERN_EDGE)  Mdiag(iend+1,:,:,2,iu3d) = uv_prev(iend+1,:,:,2) - uvdz_old(iend+1,:,:,2)

      if (SOUTHERN_EDGE) uv_prev(:,jstr  ,:,2)=u(:,jstr  ,:,nnew)*0.5*(Hz(:,jstr-1,:)+Hz(:,jstr  ,:))
      if (NORTHERN_EDGE) uv_prev(:,jend+1,:,2)=u(:,jend+1,:,nnew)*0.5*(Hz(:,jend  ,:)+Hz(:,jend+1,:))
      if (SOUTHERN_EDGE) Mdiag(:,jstr  ,:,2,iu3d) = uv_prev(:,jstr  ,:,2) - uvdz_old(:,jstr  ,:,2)
      if (NORTHERN_EDGE) Mdiag(:,jend+1,:,2,iu3d) = uv_prev(:,jend+1,:,2) - uvdz_old(:,jend+1,:,2)

      end subroutine set_diags_at_bc_v  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_FCu( istr, iendR, jstrR, jendR )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istr, iendR, jstrR, jendR

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstrR,jendR  ! since this section uv2FC uses IU_RANGE & J_RANGE which includes boundary
          do i=istr,iendR
            Mdiag(i,j,k,1,iuv2FC) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - uv_prev(i,j,k,1)

            ! Full loop ranges since u-change over every point including bry.
            ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            ! can't just do (u(nnew)-u(nstp)*Hz here because actual need the old Hz for u(nstp) as Hz has evolved.
            if (check_diag .or. tot_change)  ! if with parameters
     &        uv_dif(i,j,k, 1)   = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - uvdz_old(i,j,k, 1) ! do this here as units are correct
          enddo
        enddo
      enddo

      end subroutine set_diags_at_uv2_FCu  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_FCv( istrR, iendR, jstr, jendR )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istrR, iendR, jstr, jendR

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstr,jendR  ! since this section uv2FC uses IU_RANGE & J_RANGE which includes boundary
          do i=istrR,iendR
            Mdiag(i,j,k,2,iuv2FC) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) )) - uv_prev(i,j,k,2)

            if (check_diag .or. tot_change)  ! if with parameters
     &        uv_dif(i,j,k, 2)   = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) )) - uvdz_old(i,j,k, 2) ! do this here as units are correct
          enddo
        enddo
      enddo

      end subroutine set_diags_at_uv2_FCv  !]

! ----------------------------------------------------------------------
      subroutine set_diags_dz_old(istr, jstr, istrR, iendR, jstrR, jendR)  ![
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

      if(mynode==0) print *, 'DZ_OLD  $$$$$$$$$$'  ! debug

      end subroutine set_diags_dz_old  !]

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
# ifdef DIAGS_TRC
      subroutine set_diags_t_at_bc(istr, iend, jstr, jend, itrc, td)  ![

      use tracers, only: t
      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend, itrc, td

      ! local
      integer :: i, j

      ! PROBABLY DON'T NEED THIS CAN JUST USE T(NNEW) AT THE BOUNDARY SINCE ONLY UPDATE FROM HERE.

      ! currently this is not perfectly setup for openMP since we have ':' instead
      ! of istr:iend, etc. However, since it is just a 2D operation at bry, it is ok? for now.

      if (mynode==0) print *, 'TRACER BC /////////' ! debug

      if (WESTERN_EDGE) then
        i=istr-1
        Tdiag(i,:,:,td,it_bc) = t(i,:,:,nnew,itrc) *     Hz(i,:,:)
     &                        - t(i,:,:,nstp,itrc) * Hz_old(i,:,:)
      endif
      if (EASTERN_EDGE) then
        i=iend+1
        Tdiag(i,:,:,td,it_bc) = t(i,:,:,nnew,itrc) *     Hz(i,:,:)
     &                        - t(i,:,:,nstp,itrc) * Hz_old(i,:,:)
      endif
      if (SOUTHERN_EDGE) then
        j=jstr-1
        Tdiag(:,j,:,td,it_bc) = t(:,j,:,nnew,itrc) *     Hz(:,j,:)
     &                        - t(:,j,:,nstp,itrc) * Hz_old(:,j,:)
      endif
      if (NORTHERN_EDGE) then
        j=jend+1
        Tdiag(:,j,:,td,it_bc) = t(:,j,:,nnew,itrc) *     Hz(:,j,:)
     &                        - t(:,j,:,nstp,itrc) * Hz_old(:,j,:)
      endif

      end subroutine set_diags_t_at_bc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_t_at_t3dmix_GP(istr,iend,jstr,jend, itrc)  ![
      ! from t3dmix_GP.F

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
            Tdiag(i,j,k,td,it_t3dmix) = t(i,j,k,nnew,itrc)
          enddo
        enddo
      enddo

      end subroutine set_diags_t_at_t3dmix_GP  !]

! ----------------------------------------------------------------------
      subroutine set_diags_avg_tile(istr,iend,jstr,jend)  ![
      ! set the averages of diagnostic variables

      implicit none

      ! inputs
      integer, intent(in) :: istr,iend,jstr,jend

      ! local
      integer :: i, j, k, iMd, dim

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      if(mynode==0) print *, '--- SETTING DIAGS AVERAGES' ! debug only

      ! This method is arguably more efficient than set_avg.F, however, it requires more
      ! lines of code. Since diagnostic averaging will be heavy, more efficient code is argueably
      ! more important.

      do dim=1,2

        if (mod(tn,rec_rate)==0) then      ! timestep / rec_rate

  !        do iMd=1,iuv2DC                 ! all diagnostic terms with internal loop ranges
  !       ! for greatest efficiency should only do the loop ranges (internal or full) of the diag terms, but doing full
  !       ! ranges for now just to get working.

          if (check_avg) Md_a_sum(istrR:iendR,jstrR:jendR,:,dim) = 0

          do iMd=1,nmd                     ! all diagnostic terms with internal loop ranges
            do k=1,N                       ! here we divide the average by the averaging rate
              do j=jstrR,jendR
                do i=istrR,iendR
                  Md_avg(i,j,k,dim,iMd) = ( Md_avg(i,j,k,dim,iMd)+Mdiag(i,j,k,dim,iMd) ) / rec_rate
                  if (check_avg)           ! if with parameter
     &              Md_a_sum(i,j,k,dim) = Md_a_sum(i,j,k,dim) + Md_avg(i,j,k,dim,iMd)
                enddo
              enddo
            enddo
          enddo

          if (check_avg) then
            do k=1,N; do j=jstrR,jendR; do i=istrR,iendR ! Md_a_ever, could write it now then wouldn't need for both u and v? Same for Md_a_sum
              Md_a_ever(i,j,k,dim) = Md_a_ever(i,j,k,dim) + Md_a_sum(i,j,k,dim) * rec_rate                        ! keep an endless running total for budget closing
            enddo; enddo; enddo
            call set_diags_uvdz_ever( istr, jstr, istrR, iendR, jstrR, jendR ) ! * rec_rate since tot = avg * t
          endif

          if(mynode==0 .and. dim==1) print *, '----> DIVIDE FOR AVERAGE'

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
      subroutine set_diags_t_avg_tile(istr, iend, jstr, jend, td, iTd, itrc)  ![
      ! set the averages of tracer diagnostic variables

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, td, iTd, itrc ! tn  = timestep
                                                                   ! td  = tracer (local diag index)
      ! local                                                      ! iTd = diagnostic term
      integer           :: i, j, k, ierr                           ! itrc = tracer index as per 't' array
      character(len=30) :: text

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      if(mynode==0) print *, '^^^ SETTING TRC DIAGS AVERAGES' ! debug only

      if (mod(tn,rec_rate)==0) then      ! timestep / rec_rate

        if (check_avg .and. iTd==1) Td_a_sum(istrR:iendR,jstrR:jendR,:) = 0  ! reset for new tracer

        do k=1,N                       ! here we divide the average by the averaging rate
          do j=jstrR,jendR
            do i=istrR,iendR
              Td_avg(i,j,k,td,iTd) = ( Td_avg(i,j,k,td,iTd)+Td_tmp(i,j,k) ) / rec_rate
              if (check_avg)           ! if with parameter
     &          Td_a_sum(i,j,k) = Td_a_sum(i,j,k) + Td_avg(i,j,k,td,iTd)
            enddo
          enddo
        enddo

        text=trim(t_vname(itrc)) / / tdname(iTd)
        call nc_write_var(ncid, Td_avg(:,:,:,td,iTd), N, text, rp_var, output_rec, ierr)

        if (check_avg .and. iTd==ndt) then
          do k=1,N; do j=jstrR,jendR; do i=istrR,iendR
            Td_a_ever(i,j,k,td) = Td_a_ever(i,j,k,td) + Td_a_sum(i,j,k) * rec_rate  ! keep an endless running total for budget closing
          enddo; enddo; enddo                                                       ! * rec_rate since tot = avg * t
          text=trim(t_vname(itrc)) / / check_avnm(2)
          call nc_write_var(ncid, Td_a_ever(:,:,:,td), N, text, rp_var, output_rec, ierr)

          call set_diags_tdz_ever( istrR, iendR, jstrR, jendR, td, itrc )
          text=trim(t_vname(itrc)) / / check_avnm(1)
          call nc_write_var(ncid, Td_a_sum(:,:,:), N, text, rp_var, output_rec, ierr)
        endif

        if(mynode==0) print *, '----> TRACER DIVIDE FOR AVERAGE'

      else

        if (mod(tn-1,rec_rate)==0) then   ! this is next step after we have averaged.

          do k=1,N                      ! here we divide the average by the averaging rate
            do j=jstrR,jendR
              do i=istrR,iendR
                Td_avg(i,j,k,td,iTd) = Td_tmp(i,j,k)
              enddo
            enddo
          enddo

        else                              ! all other timesteps

          do k=1,N                      ! here we divide the average by the averaging rate
            do j=jstrR,jendR
              do i=istrR,iendR
                Td_avg(i,j,k,td,iTd) = Td_avg(i,j,k,td,iTd) + Td_tmp(i,j,k)
              enddo
            enddo
          enddo

        endif  ! <-- mod(tn-1,rec_rate)==0

      endif    ! <-- mod(tn  ,rec_rate)==0

      end subroutine set_diags_t_avg_tile  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uvdz_ever( istr, jstr, istrR, iendR, jstrR, jendR )  ![
      ! need current dz*u/v to compare gross changes in u with sum of averages

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
      ! since the value of the tracer has been stored as dz*t, we need to find the
      ! differences for the contribution of each diagnostic term

      ! Unlike surf_flux.F -> wrt_surface_flux_his, here we still need to calculate
      ! averages, and so the if (write_file) is divided.
      ! Also, since the value of 't' is kept in Tdiag and not the change in 't',
      ! we use Td_tmp to store the change, and for memory management it is only
      ! one array, and thus needs to be written to file immediately, or added to averges.
      ! In wrt_diagnostics_tile the if(write_file) is split where it is so that the file can
      ! be opened and time written before writing averages if applicable.

      use tracers, only: t
      implicit none

      ! inputs
      integer, intent(in)                 :: istr, iend, jstr, jend

      ! local
      integer                             :: itrc, td, i, j, k, ierr, it_prev
      character(len=30)                   :: text

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      do td=1,nd_trc                                          ! only selected tracers from wrt_t_diag array
        itrc=diag_td_2_NT(td)                                 ! set the actual tracer index in 't'

        if (check_diag) T_dsum(istrR:iendR,jstrR:jendR,:)=0  ! could just do boundaries as the rest can be overwritten in first below
        Td_tmp(istrR:iendR,jstrR:jendR,:)=0                   ! do more efficient reset of boundaries here. Above changes bry values but below doesn't use them but can be seen in output if not set to 0 again.

        do k=1,N                                              ! first diagnostic calculated using previous step's value
          do j=jstr,jend
            do i=istr,iend
              Td_tmp(i,j,k)=Tdiag(i,j,k,td,it_h_adv) - t(i,j,k,nstp,itrc)*Hz_old(i,j,k)
              if (check_diag)
     &        T_dsum(i,j,k)=Td_tmp(i,j,k)                     ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(istr, iend, jstr, jend, td, it_h_adv, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_h_adv)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              Td_tmp(i,j,k)=Tdiag(i,j,k,td,it_v_adv) - Tdiag(i,j,k,td,it_h_adv)
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(istr, iend, jstr, jend, td, it_v_adv, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_v_adv)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif
        it_prev=it_v_adv ! store index of previous diagnostic

# ifdef LMD_KPP
        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              Td_tmp(i,j,k)=Tdiag(i,j,k,td,it_lmdkpp) - Tdiag(i,j,k,td,it_v_adv)
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(istr, iend, jstr, jend, td, it_lmdkpp, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_lmdkpp)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif
        it_prev=it_lmdkpp ! store index of previous diagnostic if lmd_kpp
# endif

        ! From here onwards 't' was no longer in dz*t units:
        do k=1,N
          do j=jstr,jend
            do i=istr,iend  ! it_imp_vdif in 't' units, but it_prev in dz*t units
              Td_tmp(i,j,k)=Tdiag(i,j,k,td,it_imp_vdif)*Hz(i,j,k) - Tdiag(i,j,k,td,it_prev)
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(istr, iend, jstr, jend, td, it_imp_vdif, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_imp_vdif)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        do k=1,N
          do j=jstr,jend
            do i=istr,iend  ! both in dz*t units
              Td_tmp(i,j,k)=(Tdiag(i,j,k,td,it_t3dmix)-Tdiag(i,j,k,td,it_imp_vdif))*Hz(i,j,k)
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg_tile(istr, iend, jstr, jend, td, it_t3dmix, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_t3dmix)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! Do boundary changes - already calculated due to need for Hz and Hz_old to differ (they don't when wrt_online_diags is called)
        if (diag_avg) then
          Td_tmp(istrR:iendR,jstrR:jendR,:)=Tdiag(istrR:iendR,jstrR:jendR,:,td,it_bc)
          call set_diags_t_avg_tile(istr, iend, jstr, jend, td, it_bc, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_bc)
          call nc_write_var(ncid, Tdiag(START_2D_ARRAY,1,td,it_bc), N, text, rp_var, output_rec, ierr)
        endif
        if (check_diag)
     &    T_dsum(istrR:iendR,jstrR:jendR,:)=T_dsum(istrR:iendR,jstrR:jendR,:)
     &                                      +Tdiag(istrR:iendR,jstrR:jendR,:,td,it_bc)  ! re-do this just for boundaries not full domain!

        ! need rec_rate check here as not an average
        if ( (check_diag .or. tot_change) .and. mod(tn,rec_rate) == 0) then
          do k=1,N                    ! write gross difference: don't need to set Td_sum to zero as covering full loop ranges.
            do j=jstrR,jendR          ! first diagnostic calculated using previous step's value
              do i=istrR,iendR        ! Full domain loop range incl. bry.
                Td_tmp(i,j,k)=t(i,j,k,nnew,itrc)*Hz(i,j,k) - t(i,j,k,nstp,itrc)*Hz_old(i,j,k)
              enddo
            enddo
          enddo
          text=trim(t_vname(itrc)) / / check_vnm(1)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! need rec_rate check here as not an average
        if (check_diag .and. mod(tn,rec_rate) == 0) then          ! following to compare gross change in u/v and sum of diags.
          text=trim(t_vname(itrc)) / / check_vnm(2)
          call nc_write_var(ncid, T_dsum, N, text, rp_var, output_rec, ierr)

          Td_tmp(istrR:iendR,jstrR:jendR,:)=Td_tmp(istrR:iendR,jstrR:jendR,:)
     &                                     -T_dsum(istrR:iendR,jstrR:jendR,:)  ! difference between gross change and sum of diags.
          text=trim(t_vname(itrc)) / / check_vnm(3)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif  ! <-- check_diag

      enddo    ! <-- itrc=1,NT

      Hz_old(istrR:iendR,jstrR:jendR,:) = Hz(istrR:iendR,jstrR:jendR,:)  ! update Hz_old

!      if(WESTERN_EDGE) then  ! just here for jstrR, shouldn't be in tracers
!        call wrt_1D_or_2D_array(ncid, output_rec, 'u_wbry', jstrR, jendR, 1, N, uv_wbry(jstrR:jendR,1:N,1) )
!      endif

      end subroutine wrt_tracer_diags_tile  !]
# endif /* DIAGS_TRC */

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
      ! main subroutine to write online diagnostics

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
      integer, dimension(3)   :: rwebrygrd(3)                  ! experimental for if using 2D bry outputs
      integer, dimension(4,2) :: grdtyp3d                      ! u and v type to using in loop
      integer, dimension(2)   :: gt_wrt = (/ up_var, vp_var /) ! grid type for writing

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      tn = iic-ntstart      ! model step, as timestep t=n. (icc = t(n+1), and ntstart=1 if not from restart). We write the previous step results during the start of the next step

      ! might be able to use diag_avg .or. calc_diag and compiler sees or is parameter!
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

        if ( write_file==.true.     .and.
     &       mod(tn,rec_rate) == 0        ) then

! Defining and writing netcdf variables only done by master thread:
C$OMP MASTER

          if(total_output_rec == 0) then                            ! Add root output name to file name (only done once)
            file_name = trim(output_root_name) / / trim(file_name)  ! Notice below space between '/ /', this is needed to avoid
          endif                                                     ! cpp preprocessor seeing // as a c++ comment and deleting it.


          if(output_rec == 0) then                                  ! if file doesn't yet exist or need new file!

            call ncdf_create_file(file_name, ncid, prev_fill_mode   ! lose the prev_fill_mode? r2dgrd, ...., auxil are id's to dimensions, needed to create vars in the file
     &                    ,recs_per_file, total_output_rec, auxil
#ifdef SOLVE3D
!     &                    ,r3dgrd, u3dgrd, v3dgrd, w3dgrd           ! 3D grid dimensions
     &                    ,r3dgrd, grdtyp3d(:,1), grdtyp3d(:,2), w3dgrd           ! 3D grid dimensions
#endif
     &                    ,r2dgrd, u2dgrd, v2dgrd )                 ! 2D grid dimensions


            ! Define remaining variables:

            do dim=1,2                                  ! u,v dims

              if (diag_uv) then
                do diag=1,nmd                           ! diag terms
                  text  = uvw(dim) / / vnames(diag)     ! diag short name
                  call nc_define_var(ncid, text, lname(diag), units, grdtyp3d(:,dim), ierr )
                enddo
              endif

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
              endif

            enddo ! <-- dim

# ifdef DIAGS_TRC
            if (diag_trc) then
              do td=1,nd_trc                                   ! only tracers of interest
                itrc=diag_td_2_NT(td)                          ! set the actual tracer index in 't'
                do idt=1,ndt                                   ! tracer diag terms
                  text  = trim(t_vname(itrc)) / / tdname(idt)
                  utext = 'd (' / / trim(t_units(itrc)) / / tunits
                  call nc_define_var(ncid, text, tlname(idt), utext, r3dgrd, ierr )
                enddo

                if (check_diag .or. tot_change) then
                  text  = trim(t_vname(itrc)) / / check_vnm(1)  ! gross change in dz*t over 1 timestep
                  call nc_define_var(ncid, text,  check_lnm(1), tunits, r3dgrd, ierr )
                endif

                if (check_diag) then
                  text  = trim(t_vname(itrc)) / / check_vnm(2)  ! sum of all diag terms over 1 timestep
                  call nc_define_var(ncid, text,  check_lnm(2), tunits, r3dgrd, ierr )
                  text  = trim(t_vname(itrc)) / / check_vnm(3)  ! compare gross change against sum of diags. Should be approx zero.
                  call nc_define_var(ncid, text,  check_lnm(3), tunits, r3dgrd, ierr )
                endif

                if (check_avg) then
                  text  = trim(t_vname(itrc)) / / check_avnm(1)  ! gross change in t up until now since t0
                  call nc_define_var(ncid, text, check_alnm(1), tunits, r3dgrd, ierr )
                  text  = trim(t_vname(itrc)) / / check_avnm(2) ! total diag avg in t up until now since t0
                  call nc_define_var(ncid, text, check_alnm(2), tunits, r3dgrd, ierr )
                endif

              enddo  ! <-- td
            endif    ! <-- diag_trc
# endif

!            ! define bry 2D vars (experimental):
!            rwebrygrd(1)=r3dgrd(2) ! set eta_rho dimension
!            rwebrygrd(2)=r3dgrd(3) ! set s_rho dimension
!            rwebrygrd(3)=r3dgrd(4) ! set time    dimension
!!          if(WESTERN_EDGE) then
!            call nc_define_var(ncid, 'u_wbry', 'u_wbry', 'u_wbry', rwebrygrd, ierr )
!!          endif

            ierr=nf90_enddef(ncid) ! end definition to write changes to disk
            if (mynode == 0) then
              write(*,'(6x,2A)') 'diagnostics :: created new netCDF file '
     &                            ,file_name
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

        endif   ! <-- write_file==.true. .and. mod(tn,rec_rate) == 0

        if (diag_avg) then                                                 ! see wrt_diagnostics_tile for description of why if(write_file) is split here.
          if (diag_uv)  call set_diags_avg_tile   (istr, iend, jstr, jend) ! avoid very first step because diags are only all set by the end of the
          if (diag_trc) call wrt_tracer_diags_tile(istr, iend, jstr, jend) ! first step, and then added to averages and/or output at the beginning of the next step.
        endif                                                              ! Need this after defined file as wrt_tracer will output averages so needs file created.

        if ( write_file==.true.     .and.
     &       mod(tn,rec_rate) == 0        ) then


          do dim=1,2                            ! over dims

            if (diag_uv) then
              do diag=1,nmd                       ! over diag terms
                text = uvw(dim) / / vnames(diag)
                if (diag_avg) then
                  call nc_write_var(ncid, Md_avg(:,:,:,dim,diag), N, text, gt_wrt(dim), output_rec, ierr)
                else
                  call nc_write_var(ncid,  Mdiag(:,:,:,dim,diag), N, text, gt_wrt(dim), output_rec, ierr)
                endif
              enddo
            endif

            if (check_diag .or. tot_change) then
              ! do uv_dif calculation at end of step3d_uv2.F as have istr, iend there already.
              ! but use 'nnew' there as still in previou timestep!
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
              call nc_write_var(ncid, uvdz_ever(:,:,:,dim),   N, text, gt_wrt(dim), output_rec, ierr)
              text  = uvw(dim) / / check_avnm(2)                         ! total diag avg in u/v up until now since t0
              call nc_write_var(ncid, Md_a_ever(:,:,:,dim),   N, text, gt_wrt(dim), output_rec, ierr)
            endif

          enddo  ! <-- dim

# ifdef DIAGS_TRC
          if (diag_trc .and. .not. diag_avg)                   ! for averages this has been called already
     &      call wrt_tracer_diags_tile(istr, iend, jstr, jend) ! we call it here non-averages to get the correct writing frequency
# endif

          ierr=nf90_close (ncid)  ! close netcdf file

          if (mynode == 0) then
            write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &       'diagnostics :: wrote output, tdays =', tdays,
     &       'step =', iic-1, 'rec =', output_rec, '/', total_output_rec
     &        MYID
          endif

          if (output_rec == recs_per_file) then  ! reset output_rec to create new file
            output_rec = 0
          endif

C$OMP END MASTER
C$OMP BARRIER

        endif ! <-- write_file==.true.

      endif   ! <-- icc>ntstart

      end subroutine wrt_diagnostics_tile  !]

! ----------------------------------------------------------------------

#endif /* DIAGNOSTICS */

      end module diagnostics
