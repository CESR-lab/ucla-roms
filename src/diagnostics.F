      module diagnostics

      ! To compute diagnostics use DIAGNOSTICS flag in cppdefs.opt.
      ! Also, set user defiend configuration in diagnostics.opt.

#include "cppdefs.opt"
#ifdef DIAGNOSTICS

      use dimensions
      use read_write
      use netcdf
      use tracers
      use ocean2d
      use ocean3d
      use scalars
      use mixing  ! visc2_r
      use hidden_mpi_vars

      implicit none
      private

#include "diagnostics.opt"

      ! INFO (Optional): ![
      ! Online diagnostics for ROMS simulation.
      ! Initial coding DevinD and JeroenM - 2021/05

      ! JM: We are not properly accounting for the Udzeta/dt,Vdzeta/dt terms
      ! These should only be a signicant error for coastal tidal flows.

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
      ! But wrt_diagnostics is called before nnew and nstp indicies change in main.F.

      ! If you are needing averages then the diagnotics need to be computed at everytimestep.
      ! If you just need snapshots then the diagnostics only need to be calculated at
      ! your output frequency.
      ! Both of these are controlled by the logical 'calc_diag', which is set in
      ! wrt_diagnostics at the beginning of a timestep, and then that logical is
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
      ! The order of evolution of u/v (shown as u) and its units is as follows:
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

      !]

      ![ Potential issues:
      !  - Since we don't have pressure in the boundary buffer, we can't get u'p' at the bry.
      !    We therefore just use the first internal pressure point instead, which introduces
      !    a dx/2 error.
      !]

      ! Preamble:  ![

      logical, public            :: calc_diag               ! flag to calculate diagnostics in equation system.
      logical, public            :: set_dz_old              ! set dz*u/v/w/t at previous timestep, needed for diag calculation when skipping steps if not averages.
      logical                    :: init_done=.false.       ! flag to init diagnostics only once
      logical                    :: p_init_done=.false.     ! flag to init pressure filter value only once

      ! U & V MOMENTUM:
      ! ---------------
      real,dimension(:,:),pointer :: wrk_xy1,wrk_xy2

      integer, public, parameter                          :: nmd = 6    ! number of momentum diag terms

      real, public, allocatable, dimension(:,:,:,:) :: Udiag 
      real, public, allocatable, dimension(:,:,:,:) :: Vdiag
      real, public, allocatable, dimension(:,:,:)   :: udz_old ! dz*u
      real, public, allocatable, dimension(:,:,:)   :: vdz_old ! dz*v
      real, public, allocatable, dimension(:,:,:)   :: u_prev
      real, public, allocatable, dimension(:,:,:)   :: v_prev
      real,         allocatable, dimension(:, :)  :: FX4, FY4   ! 4th order advective fluxes to split UP3 = ADV + DISSIP. Called in k loop so only need 2D array for each k
      real,         allocatable, dimension(:, :)  :: cen2       ! 2nd order centred scheme used to handle boundary and mask.
      real,         allocatable, dimension(:, :)  :: FZ4        ! vertical advective fluxes. u/v only. FX4/FY4 both u/v and tracers.

      ! if diag_avg - for averaging of diagnostics
      real, allocatable, dimension(:,:,:,:) :: Udiag_avg ! Averages of diagnostics
      real, allocatable, dimension(:,:,:,:) :: Vdiag_avg ! Averages of diagnostics

      real, public                                        :: Mdtmp      ! temporary scalar for calculations of Mdiag

      integer, parameter, public ::
     &  iprsgr = 1, ! Hydrostatic pressure gradient    - internal loop ranges only (i.e. excluding domain boundary)
     &  icori  = 2, ! Coriolis & curvilinear coords    - internal loop ranges only
     &  iadv   = 3, ! purely advective fluxes          - internal loop ranges only
     &  idiss  = 4, ! dissipation from advection       - internal loop ranges only
     &  ihmix  = 5, ! Viscous terms (horiz mixing)     - internal loop ranges only
     &  ivmix  = 6  ! Vertical mixing & 2D/3D coupling - internal loop ranges only
      ! ivmix  <--  ! Boundary changes                 - bry only - should store 3dbc changes in 2D arrays not 3D, as just 0 in interior. Goes into ivmix array for memory efficiency.
      ! ivmix  <--  ! 2D-3D correction                 - internal + bry (full domain) loop ranges. Goes into ivmix array for memory efficiency.

      character (len=20),  allocatable, dimension(:) :: vnames
      character (len=120), allocatable, dimension(:) :: lname           ! len=70 allows for 8 characters of 'Average ' if diag_avg

      character (len=26) :: units = 'm^2/s^2    (i.e. du/dt*dz)'        ! Don't go beyond 72nd character else space in units
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

# ifdef OBC_CHECK
      real sum_dbg ! debug sum check
# endif

      integer, public, parameter                          :: ndt = 6   ! number of tracer diag terms
      real, public, allocatable, dimension(:, :, :, :, :) :: Tdiag     ! tracer diagnostic terms
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
     &  it_h_adv   = 1,  ! Horizontal UP3 purely advective fluxes - compute_horiz_tracer_fluxes.h
     &  it_h_diss  = 2,  ! Dissipative part of UP3 hori. advection
     &  it_v_adv   = 3,  ! Vertical   UP3 purely advective fluxes - compute_vert_tracer_fluxes.h
     &  it_v_diss  = 4,  ! Dissipative part of UP3 vert. advection
                         ! do we need surface forcing flux as a vertical? Shouldn't be a 3D array anyway
     &  it_v_mix   = 5,  ! Implicit time step for vertical diffusion - step3d_t_ISO.F
     &  it_h_mix   = 6   ! Horizontal (along geopotential surfaces) Laplacian diffusion terms - t3dmix_GP.F
      ! it_v_mix   <--   ! Nonlocal transport flux for unstable (convective) forcing - step3d_t_ISO.F. not used if no LMD_KPP flag, hence at the end. ! step3d_t_ISO.F. This only change T & S so shouldn't be 'NT' in size, a waste.
      ! it_v_mix   <--   ! Boundary changes - t3dbc_im.F. Only at boundary so shared with v_mix (could be any other term but being consistent with u/v choice of v_mix.

      character (len=30),  allocatable, dimension(:) :: tdname
      character (len=120), allocatable, dimension(:) :: tlname  ! extra 10 characters for 'Averaged...'

      character (len=20) :: tunits = ' * m/s   (i.e. dC/dt*dz)'

      public init_diags_trc
      public set_diags_t_bc
      public set_diags_t_h_mix
      public set_diags_t_4th_h_adv
      public set_diags_t_4th_v_adv

      ! Netcdf outputting:
      ! ------------------
      integer            :: output_rec = 0            ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0      ! total flux output recordings so far
      integer            :: tn                        ! roms timestep
      character (len=99) :: file_name = '_dia.nc'     ! suffix to append to simulation output rootname
      integer            :: ncid=-1, prev_fill_mode


      ! PRESSURE FLUX:
      ! --------------
      real, allocatable, dimension(:, :, :), public :: up_pp, vp_pp ! u'p' (at upoint) and v'p' (at v_point) - makes it easier to use in child bry.
      real, allocatable, dimension(:, :, :) :: up, vp               ! perturbed u/v/
      real, public,allocatable, dimension(:, :, :) :: uf, vf, pf           ! filtered  u/v/pressure
      real, allocatable, dimension(:, :, :) :: pp                   ! perturbed p
      real, allocatable, dimension(:, :, :) :: u_bcl, p_bcl         ! u and p baroclinic
      real, allocatable, dimension(:, :, :) :: aup, avp, app        ! averaged perturbed u/v/p
      real, allocatable, dimension(:, :, :) :: aup_pp, avp_pp       ! averaged u'p' and v'p'

      ! -- divergence terms
      real, public,allocatable, dimension(:,:) :: dxdyi_u, dxdyi_v,dxdyi ! 1/ surface of cells
      real, public,allocatable, dimension(:,:,:) :: pxf, pxp       ! dP/dx filtered and prime
      real, public,allocatable, dimension(:,:,:) :: pyf, pyp       ! dP/dy filtered and prime
      real, allocatable, dimension(:,:) :: bt2bc    ! Barotropic to baroclinic conversion
      real, allocatable, dimension(:,:) :: divup    ! Divergence of up flux
      real, public,allocatable, dimension(:,:) :: udh_f, udh_p           ! ubar.grad(h), filtered and prime
      real, public,allocatable, dimension(:,:) :: dzt_f, dzt_p           ! dzeta/dt, filtered and prime
      real, allocatable, dimension(:,:,:) :: div_uvp, adiv_uvp    ! divergence u'p' + v'p' (and average)

      ! -- sponge terms:
      real, allocatable, dimension(:)       :: mup_ppzw,  mup_ppze  ! mean of up_ppz along orthogonal sponge from west/est boundary
      real, allocatable, dimension(:)       :: mvp_ppzs,  mvp_ppzn  ! mean of up_ppz along orthogonal sponge from south/north boundary
      real, allocatable, dimension(:)       :: box_mup_ppzw,  box_mup_ppze  ! box average
      real, allocatable, dimension(:)       :: box_mvp_ppzs,  box_mvp_ppzn
      real, allocatable, dimension(:)       :: fbox_mup_ppzw, fbox_mup_ppze ! filtered box average
      real, allocatable, dimension(:)       :: fbox_mvp_ppzs, fbox_mvp_ppzn
      real, allocatable, dimension(:)       :: sponge_fuw,    sponge_fue    ! time filtered mean of up_ppz along orthogonal sponge from west/est boundary
      real, allocatable, dimension(:)       :: sponge_fvs,    sponge_fvn
      real, allocatable, dimension(:, :)    :: visc2_r0,  visc2_p0  ! Initial sponges
      real, allocatable, dimension(:, :, :) :: diff20               ! Initial tracer sponge

      integer :: nprm = 3                      ! number of pressure flux diagnostic terms
      real, public :: alpha, beta              ! filtering constants (needed for child bry)
      integer :: sponge_loop                   ! for looping through sponge layer
      integer :: hbox_size = (box_size-1) / 2  ! half the box size. -1 since averaging point in middle, loop range either side

      logical, public :: pg_cor = .false.      ! calculate pressure on corrector step of prsgrd. false initially.

      character (len=20),  allocatable, dimension(:) :: pnames
      character (len=120), allocatable, dimension(:) :: plname           ! len=70 allows for 8 characters of 'Average ' if diag_avg
      character (len=20),  allocatable, dimension(:), public :: ppnames
      character (len=120), allocatable, dimension(:), public :: pplname

      character (len=3), public :: punits = 'W/m  (i.e. up*dz)'
      character (len=3) :: pchunits = 'TBC'

      public set_diags_uv_prime
!     public set_diags_prs_prime

      ! W MOMENTUM: (CURRENTLY NOT WORKING)
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

      public set_diags_u_4th_adv
      public set_diags_v_4th_adv
      public set_diags_uv_dz_old
      public wrt_diagnostics !]

      contains

! ----------------------------------------------------------------------
      subroutine init_diagnostics ![
      ! Allocatable and initialize diagnostic arrays.
      implicit none


      ! local
      integer :: i,j,k,diag

      ! Do some consistency checking
      if(check_avg .and. .not. diag_avg) then                           ! catch bad config
        write(*,*) 'ERROR - diagnostics.F: check_avg must be false ',
     &  'if diag_avg is false!'
        error stop 'ERROR - diagnostics.F: check_avg!=F if diag_avg==F'
      endif

      if(diag_avg .and. rec_rate==1) then
        write(*,*) 'ERROR - diagnostics.F: need rec_rate > 1 if ',
     &  'using averages!'
        error stop 'ERROR - diagnostics.F: diag_avg==T then rec_rate>1'
      endif

      if( diag_trc == .true.  .and.               ! requesting tracer diags
     &   .not. any( wrt_t_diag == .true.) ) then  ! but no tracers selected

        write(*,*) 'ERROR - diagnostics.F: diag_trc==T but ',
     &  'no tracers selected in wrt_t_diag'
        error stop 'ERROR - diagnostics.F: diag_avg==T but wrt_t_diag=F'

      endif

      if (prs_flx) then
        allocate( udh_f(nx,ny) )
        allocate( dzt_f(nx,ny) )

        allocate( pxf(nx,ny,nz) )
        allocate( pyf(nx,ny,nz) )
        
        !!! JM this is a test for reusable space
        allocate(wrk_xy1(nx,ny))
        allocate(wrk_xy2(nx,ny))

        allocate( dxdyi_u(nx,ny) )
        allocate( dxdyi_v(nx,ny) )
        allocate( dxdyi(nx,ny) )
        dxdyi_u = 0.25*( pn(0:nx-1,1:ny)+pn(1:nx,1:ny) )
     &                *( pm(0:nx-1,1:ny)+pm(1:nx,1:ny) )
        dxdyi_v = 0.25*( pn(1:nx,0:ny-1)+pn(1:nx,1:ny) )
     &                *( pm(1:nx,0:ny-1)+pm(1:nx,1:ny) )
        dxdyi = pn(1:nx,1:ny)*pm(1:nx,1:ny)
      endif

      if (diag_uv) then

        allocate(Udiag(nx,ny,nz,nmd) )
        allocate(Vdiag(nx,ny,nz,nmd) )
        allocate(udz_old(nx,ny,nz))
        allocate(vdz_old(nx,ny,nz))
        allocate(u_prev(nx,ny,nz) )
        allocate(v_prev(nx,ny,nz) )

        allocate( vnames(nmd), lname(nmd) )
        vnames(1) = '_pgr'
        lname (1) = 'Hydrostatic pressure gradient'
        vnames(2) = '_cor'
        lname (2) = 'Coriolis & curvilinear'
        vnames(3) = '_adv'
        lname (3) = 'Advection (horiz. & vert.)'
        vnames(4) = '_dis'
        lname (4) = 'Dissipation from advection (horiz. & vert.)'
        vnames(5) = '_h_mix'
        lname (5) = 'Explicit Horizontal mixing'
        vnames(6) = '_v_mix_and_BC'
        lname (6) = 'Vertical mixing plus 2D/3D coupling'

        if (diag_avg) then   
          allocate(Udiag_avg(nx,ny,nz,nmd))
          allocate(Vdiag_avg(nx,ny,nz,nmd))
          do diag=1,nmd                                   
            lname(diag)  = 'Ave ' / / trim(lname(diag))   ! prepend 'Ave'
          enddo
        endif

        if (set_dz_old) then
          do k=1,nz
            do j=1,ny  
              do i=1,nx
                udz_old(i,j,k) = u(i,j,k,nstp)*0.5*( Hz(i-1,j,k)+Hz(i,j,k) )
                vdz_old(i,j,k) = v(i,j,k,nstp)*0.5*( Hz(i,j-1,k)+Hz(i,j,k) )
              enddo
            enddo
          enddo
        endif

      endif  ! <-- diag_uv

# ifdef NHMG
      if (set_dz_old) then
        do k=1,nz
          do j=1,ny      ! Full loop ranges since u-change over every point including bry.
            do i=1,nx
              if(k==nz) then
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

      allocate( FX4( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( FY4( GLOBAL_2D_ARRAY ) )       ! called within k loop so only need 2D slice per k
      allocate( FZ4( GLOBAL_1DX_ARRAY, 0:N ) ) ! both u/v/tracer. vertical flux at z_w level.
      FZ4(:,0) = 0  ! top and bottom vertical fluxes are always zero
      FZ4(:,nz)= 0
      allocate(cen2( GLOBAL_2D_ARRAY ) )

      if (diag_trc) call init_diags_trc  
      if (prs_flx ) call init_diags_prs_flx

      init_done = .true.  ! ensure init is not triggered again

      end subroutine init_diagnostics  !]

! ----------------------------------------------------------------------
      subroutine init_diags_trc ![
      ! setup tracer arrays only for desired tracers,
      ! since can be memory intensive if all tracers.

      implicit none

      ! local
      integer                :: itrc, cnt=0, diag, td, i, j, k
      integer, dimension(NT) :: temp_t_indx


      allocate( tdname(ndt), tlname(ndt) )
      tdname(it_h_adv)    = '_h_adv'
      tlname(it_h_adv)    =
     &  'Horizontal UP3 purely advective fluxes'                  / /
     &  ' - compute_horiz_tracer_fluxes.h'
      tdname(it_h_diss)   = '_h_diss'
      tlname(it_h_diss)   =
     &  'Horizontal UP3 advection dissipative fluxes'             / /
     &  ' - compute_horiz_tracer_fluxes.h'
      tdname(it_v_adv)    = '_v_adv'
      tlname(it_v_adv)    =
     &  'Vertical purely advective fluxes'                        / / ! spline in vertical
     &  ' - compute_vert_tracer_fluxes.h'
      tdname(it_v_diss)   = '_v_diss'
      tlname(it_v_diss)   =
     &  'Vertical advection dissipative fluxes'                   / /
     &  ' - compute_vert_tracer_fluxes.h'
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


      if (diag_avg) then
        allocate(  Td_avg( GLOBAL_2D_ARRAY, N, nd_trc, ndt ) )
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
            Tdz_t0(:,:,:,td) = t(:,:,:,nstp,itrc) * Hz(:,:,:)  ! Set t*dz at t=0 for overall budget check of averages
        enddo
      endif

!      if(mynode==0) print *, 'nd_trc=',nd_trc ! debug
!      if(mynode==0) print *, 'diag_td_2_NT(1)=',diag_td_2_NT(1)

      end subroutine init_diags_trc  !]

! ----------------------------------------------------------------------
      subroutine init_diags_prs_flx
      implicit none

      ! local
      integer :: i,j,k,var

      alpha = dt / timescale
      beta  = 1 - alpha
      if(mynode==0) write(*,'(1x,A,F9.6,A,F9.6/)')
     &  'Diagnostics: - prs_flx filtering coefficients alpha=',
     &   alpha,', beta=',beta

      allocate( uf(nx+1,ny,nz) )
      allocate( vf(nx,ny+1,nz) )
      allocate( pf(0:nx,0:ny,nz) ) !! Need buffer values for rho2u/v

      do k=1,nz
        do j=1,ny
          do i=1,nx
            uf(i,j,k) = u(i,j,k,nstp) - ubar(i,j,kstp)  ! set to initial condition, to avoid ramp up
            vf(i,j,k) = v(i,j,k,nstp) - vbar(i,j,kstp)  ! although model should ramp up still to get better average.
          enddo
        enddo
      enddo

      allocate( up(GLOBAL_2D_ARRAY, N) )  ! Don't necessaryily need up, vp and pp,
      allocate( vp(GLOBAL_2D_ARRAY, N) )  ! or can perhaps reuse array, or 2D arrays.
      allocate( pp(GLOBAL_2D_ARRAY, N) )
      allocate( up_pp(GLOBAL_2D_ARRAY, N) )
      allocate( vp_pp(GLOBAL_2D_ARRAY, N) )
      allocate( div_uvp(GLOBAL_2D_ARRAY, N) )
      up=0; vp=0; pp=0
      up_pp=0; vp_pp=0

      if (check_prs) then
        allocate( u_bcl(GLOBAL_2D_ARRAY, N) )
        allocate( p_bcl(GLOBAL_2D_ARRAY, N) )
        u_bcl=0; p_bcl=0
      endif

      if (diag_avg) then
        allocate( aup_pp(GLOBAL_2D_ARRAY, N) )
        allocate( avp_pp(GLOBAL_2D_ARRAY, N) )
        allocate( adiv_uvp(GLOBAL_2D_ARRAY, N) )
        aup_pp=0; avp_pp=0; adiv_uvp=0

        if (check_prs) then
          allocate( aup(GLOBAL_2D_ARRAY, N) )
          allocate( avp(GLOBAL_2D_ARRAY, N) )
          allocate( app(GLOBAL_2D_ARRAY, N) )
          aup=0; avp=0; app=0
        endif
      endif

      allocate( ppnames(3), pplname(3) )
      ppnames(1)    = 'up'  ! change the name here!!!
      pplname(1)    = "high frequency terms from u dot P"
      ppnames(2)    = "vp"
      pplname(2)    = "high frequency terms from v dot P"
      ppnames(3)    = "div_uvp"
      pplname(3)    = "high frequency div(up)"
      if (diag_avg) then
        do var=1,3                                            ! loop through diag terms
          pplname(var)  = 'Average ' / / trim(pplname(var))   ! and prepend 'Average'
        enddo
      endif

      if (p_sponge) then  ! define only if contains domain boundary sponge layer

        ! At some point up_pp could be integrated directly into mup_ppzw, skipping the need for up_ppwz
        if (.not.west_msg_exch)  then
          allocate( mup_ppzw (-hbox_size GLOBAL_1DY_ARRAY + hbox_size ) )  ! space on either side for MPI_exchange
          allocate( sponge_fuw(GLOBAL_1DY_ARRAY    ) )                     ! to boundary neighbours
          allocate( box_mup_ppzw(GLOBAL_1DY_ARRAY  ) )   ! reset to zero each step
          mup_ppzw=0; sponge_fuw=init_sponge             ! set time filtered initial to mid point?
        endif
        if (.not.east_msg_exch)  then
          allocate( mup_ppze (-hbox_size GLOBAL_1DY_ARRAY + hbox_size ) )
          allocate( sponge_fue(GLOBAL_1DY_ARRAY     ) )
          allocate( box_mup_ppze(GLOBAL_1DY_ARRAY  ) )
!          allocate( fbox_mup_ppze(GLOBAL_1DY_ARRAY ) )
          mup_ppze=0; sponge_fue=init_sponge
        endif

        if (.not.south_msg_exch)  then
          allocate( mvp_ppzs (-hbox_size GLOBAL_1DX_ARRAY + hbox_size ) )
          allocate( sponge_fvs(GLOBAL_1DX_ARRAY     ) )
          allocate( box_mvp_ppzs(GLOBAL_1DX_ARRAY  ) )
          mvp_ppzs=0; sponge_fvs=init_sponge
        endif
        if (.not.north_msg_exch)  then
          allocate( mvp_ppzn (-hbox_size GLOBAL_1DX_ARRAY + hbox_size ) )
          allocate( sponge_fvn(GLOBAL_1DX_ARRAY     ) )
          allocate( box_mvp_ppzn(GLOBAL_1DX_ARRAY  ) )
          mvp_ppzn=0; sponge_fvn=init_sponge
        endif

        allocate( visc2_r0(GLOBAL_2D_ARRAY    ) )
        visc2_r0 = visc2_r                               ! store initial sponge values
        visc2_r  = visc2_r * init_sponge                 ! since sponge will kick in at 0.5*sponge0 after step 1, set initial to this
        allocate( visc2_p0(GLOBAL_2D_ARRAY    ) )        ! don't need this unless wanting to output initial value
        visc2_p0 = visc2_p                               ! store initial sponge values
        visc2_p  = visc2_p * init_sponge
        allocate( diff20  (GLOBAL_2D_ARRAY, NT) )
        diff20   = diff2                                 ! consider dropping NT since all the same, though this may change in future
        diff2    = diff2   * init_sponge

        sponge_loop=sponge_size-1                        ! sponge loop size, -1 since e.g. 10 sponge = 1 + 9 iterations

        ! check box_size is odd number, since including current point in middle:

        if(mynode==0) write(*,'(14x,A,I2/)')
     &    '- dynamic sponge spatial averaging box_size =',box_size

        if(mod(box_size,2)==0)
     &    error stop 'diags: box_size must be odd number!'


      endif  ! <-- p_sponge

      if (check_prs) then
        allocate( pnames(nprm), plname(nprm) )
        pnames(1)    = 'uprm'
        plname(1)    = "perturbation of u"
        pnames(2)    = "vprm"
        plname(2)    = "perturbation of v"
        pnames(3)    = "pprm"
        plname(3)    = "perturbation of pressure"
        if (diag_avg) then
          do var=1,3                                          ! loop through diag terms
            plname(var)  = 'Average ' / / trim(plname(var))   ! and prepend 'Average'
          enddo
        endif
      endif



      end subroutine init_diags_prs_flx  !]
! ----------------------------------------------------------------------
      subroutine set_diags_u_4th_adv  ![
      ! Computes 4th order advection of u
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th


      !
      ! See readme: Documentation/readme-specifics/README.grid.fluxes.2021.07
      !             section 'U-POINT FLUXES'
      !
      ! Also, see set_diags_t_4th_h_adv for explanation on why we treat the boundaries
      ! differently. Note that masking is not needed since u is always zero in the mask.
      !
      ! HORIZONTAL ADV:
      !
      !  UFX4(xi_rho, eta_rho) is at r-points
      !  UFY4(xi_u,   eta_v  ) is between v-points!
      !
      !               UFy(i,j+1)                         <---- eta_v
      !                 ^
      !  | u(i-1,j) | u(i,j)   | u(i+1,j) | u(i+2,j) |   <---- eta_rho
      !             ^          ^
      !            UFX(i-1)   UFX(i)
      !
      ! where UFX is at rho-points (faces of u-cell), and u & FlxU are at u-points.
      !
      ! 4th order (centered) fluxes:
      !
      ! (A) UFx(i) = 1/2*(FlxU(i)+FlxU(i+1)) * (-u(i-1)+7u(i)+7u(i+1)-u(i+2)) / 12
      !
      ! where FlxU is at u-points and hence need u2rho averaging.
      ! To get from a volume flux to a flux per horizontal m^2, we divide by the
      ! sub-domain length in x and y, which is the same as multiplying by pm & pn.
      ! Since they are at rho-points, we need to average pm and pn to u-points:
      !   1/2*(pm(i,j)+pm(i-1,j)) * 1/2*(pn(i,j)+pn(i-1,j))
      !
      ! Hence divergence of the fluxes: (rho2u averaging)
      !
      !   u_h_adv_4th = dt*pm*pn*( UFX(i,j) - UFX(i-1,j  )
      !                           +UFY(i,j) - UFY(i  ,j-1) )
      !
      ! VERTICAL ADV:
      !
      ! vertical flux is at w points and thus needs to be averaged vertically to get
      ! to rho points:
      !
      !   w   x   w  -> k=1    --> depth = z_w(i,j,1)
      !
      !   r   uv  r  -> k=1    --> depth = z_r(i,j,1)
      !
      !   w   x   w  -> k=0    --> depth = z_w(i,j,0)
      !   ^   ^   ^
      !  i=1 i=2 i=2
      !
      !  Hence FlxW at x(i) = 1/2 (w(i)+w(i+1)) and where w(i)=We(i)+Wi(i) in roms
      !
      !   UFZ(i,k) = 1/2*(FlxW(i-1)+FlxW(i)) * (-u(k-1)+7u(k)+7u(k+1)-u(k+2)) / 12
      !
      !     u(k+2) -z_r
      !     u(k+1) -z_r
      !   FZ4(k  ) -z_w
      !     u(k  ) -z_r
      !     u(k-1) -z_r
      !
      !]

      implicit none

      ! local
      integer :: i,j,k,inv24

      inv24 = 1./24.  ! 0.5*1/12
      do k=1,nz

        ! The Uu flux is to the east of u(i,j)
        do j=1,ny
          do i=0,nx
            FX4(i,j) = inv24*( FlxU(i,j,k) + FlxU(i+1,j,k) ) 
     &                * (-u(i-1,j,k,nrhs) + 7*u(i  ,j,k,nrhs)
     &                   -u(i+2,j,k,nrhs) + 7*u(i+1,j,k,nrhs) )
          enddo
        enddo
        ! The Vu flux is to the south of u(i,j)
        do j=1,ny+1
          do i=1,nx
            FY4(i,j) = inv24*( FlxV(i-1,j,k) + FlxV(i,j,k) ) 
     &                * (-u(i  ,j-2,k,nrhs) + 7*u(i  ,j-1,k,nrhs)
     &                   -u(i  ,j+1,k,nrhs) + 7*u(i  ,j  ,k,nrhs) )
          enddo
        enddo

        ! The first/last interior flux can only be 2nd order 
        if (.not.west_msg_exch) then 
          i=1
          do j=1,ny
            FX4(i,j) = 0.25*( FlxU(i,j,k) + FlxU(i+1,j,k) )
     &                     *( u(i,j,k,nrhs) + u(i+1,j,k,nrhs) ) 
          enddo
        endif
        if(.not.east_msg_exch) then
          i=nx
          do j=1,ny
            FX4(i,j)= 0.25*( FlxU(i,j,k) + FlxU(i+1,j,k) )
     &                    *( u(i,j,k,nrhs) + u(i+1,j,k,nrhs) )
          enddo
        endif
        if(.not.south_msg_exch) then
          j=1
          do i=1,nx
            FY4(i,j) = 0.25*( FlxV(i-1,j,k) + FlxV(i,j,k) )    
     &                     *( u(i,j,k,nrhs) + u(i,j-1,k,nrhs) )
          enddo
        endif
        if(.not.north_msg_exch) then
          j=ny+1
          do i=1,nx
            FY4(i,j) = 0.25*( FlxV(i-1,j,k) + FlxV(i,j,k) )
     &                     *( u(i,j,k,nrhs) + u(i,j-1,k,nrhs) )
          enddo
        endif

        do j=1,ny
          do i=1,nx
            Udiag(i,j,k,iadv) = -dxdyi_u(i,j)
     &                    *( FX4(i,j)-FX4(i-1,j)+FY4(i,j+1)-FY4(i,j) )
          enddo
        enddo
      enddo  ! <-- end k-loop

      ! Vertical 4th order advection
      ! The Wu flux is above of u(i,j,k), Wu(0) and Wu(nz) are always zero

      do j=1,ny

        do k=2,nz-2
          do i=1,nx
            FZ4(i,k) = inv24*( We(i-1,j,k)+We(i,j,k)+Wi(i-1,j,k)+Wi(i,j,k) )
     &                   * (-u(i,j,k-1,nrhs) + 7*u(i,j,k  ,nrhs)  
     &                      -u(i,j,k+2,nrhs) + 7*u(i,j,k+1,nrhs) )
          enddo
        enddo

        ! Only 2nd order for first/last interior flux
        do k=1,nz-1,nz-2
          do i=1,nx
            FZ4(i,k)= 0.25*( We(i-1,j,k)+We(i,j,k)+Wi(i-1,j,k)+Wi(i,j,k) )
     &                     *( u(i,j,k,nrhs)+u(i,j,k+1,nrhs) ) 
          enddo
        enddo

        do k=1,nz
          do i=1,nx
            Udiag(i,j,k,iadv) = Udiag(i,j,k,iadv) 
     &                        - dxdyi_u(i,j)*( FZ4(i,k)-FZ4(i,k-1) )
          enddo
        enddo
      enddo  ! <--end j-loop

      end subroutine set_diags_u_4th_adv  !]

! ----------------------------------------------------------------------
      subroutine set_diags_v_4th_adv  ![
      ! Compute the 4th order advection terms for v
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! local
      integer :: i,j,k,inv24
      real   ,dimension(:,:),pointer :: Uv,Vv

      Uv => wrk_xy1
      Vv => wrk_xy2

      inv24 = 1./24.  ! 0.5*1/12
      do k=1,nz

        ! the Uv flux is to the west of v(i,j)
        do j=1,ny
          do i=1,nx+1 
            FX4(i,j) = inv24*( FlxU(i,j-1,k)+FlxU(i,j,k) )
     &                * (-v(i-2,j,k,nrhs) + 7*v(i-1,j,k,nrhs)
     &                   -v(i+1,j,k,nrhs) + 7*v(i  ,j,k,nrhs) )
          enddo
        enddo
        ! the Vv flux is to the north of v(i,j)
        do j=0,ny 
          do i=1,nx
            FY4(i,j) = inv24*( FlxV(i,j,k) + FlxV(i,j+1,k) )
     &                * (-v(i,j-1,k,nrhs) + 7*v(i,j  ,k,nrhs)
     &                   -v(i,j+2,k,nrhs) + 7*v(i,j+1,k,nrhs) )
          enddo
        enddo

        ! The first/last interior flux can only be 2nd order 
        if (.not.west_msg_exch) then 
          do j=1,ny
            i=1
            FX4(i,j) = 0.25*( FlxU(i,j-1,k)+FlxU(i,j,k) )
     &                     *( v(i-1,j,k,nrhs)+v(i,j,k,nrhs) )
          enddo
        endif
        if(.not.east_msg_exch) then
          i=nx+1
          do j=1,ny
            FX4(i,j) = 0.25*( FlxU(i,j-1,k)+FlxU(i,j,k) )
     &                     *( v(i-1,j,k,nrhs)+ v(i,j,k,nrhs) )
          enddo
        endif

        if(.not.south_msg_exch) then
          j=1
          do i=1,nx
            FY4(i,j) = 0.25*( FlxV(i,j,k) + FlxV(i,j+1,k) )
     &                     *( v(i,j,k,nrhs) + v(i,j+1,k,nrhs) )
          enddo
        endif
        if(.not.north_msg_exch) then
          j=ny
          do i=1,nx
            FY4(i,j) = 0.25*( FlxV(i,j,k) + FlxV(i,j+1,k) ) 
     &                     *( v(i,j,k,nrhs) + v(i,j+1,k,nrhs) )
          enddo
        endif

        do j=1,ny
          do i=1,nx
            Vdiag(i,j,k,iadv) = -dxdyi_v(i,j)
     &                    *( FX4(i+1,j)-FX4(i,j)+FY4(i,j)-FY4(i,j-1) ) 
          enddo
        enddo
      enddo  ! <-- k

      ! Vertical 4th order advection
      ! The Wu flux is above of u(i,j,k), Wu(0) and Wu(nz) are always zero

      do j=1,ny

        do k=2,nz-2
          do i=1,nx
            FZ4(i,k) = inv24*( We(i,j-1,k)+Wi(i,j-1,k)+We(i,j,k)+Wi(i,j,k) )
     &                   * (-v(i,j,k-1,nrhs) + 7*v(i,j,k  ,nrhs)    
     &                      -v(i,j,k+2,nrhs) + 7*v(i,j,k+1,nrhs) )  
          enddo
        enddo

        ! Only 2nd order for first/last interior flux
        do k=1,nz-1,nz-2 ! do k=1 and k=nz-1
          do i=1,nx
            FZ4(i,k) = 0.25*( We(i,j-1,k)+Wi(i,j-1,k)+We(i,j,k)+Wi(i,j,k) )
     &                      *( v(i,j,k,nrhs)+v(i,j,k+1,nrhs) )

          enddo
        enddo

        do k=1,nz
          do i=1,nx
            Vdiag(i,j,k,iadv) = Vdiag(i,j,k,iadv)         
     &                        - dxdyi(i,j)*( FZ4(i,k)-FZ4(i,k-1) )
          enddo
        enddo

      enddo      ! <-- j

      end subroutine set_diags_v_4th_adv  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uv_dz_old  ![
      ! Set dz*u or v at previous timestep to calculate gross change in
      ! u or v over time step, and to calculate change at boundary.
      ! Called after everything else in step3d_uv2.F.

      implicit none

      ! local
      integer :: i,j,k

      ! This is currently used for bry calculation so is also needed for averaging.
      ! Also for 'u_diff', but that's not needed for averaging.
      ! Should put an if statement here for diags_avg and only compute bry if averages...

      do k=1,nz
        do j=1,ny
          do i=1,nx
            udz_old(i,j,k) = u(i,j,k,nnew)*0.5*( Hz(i-1,j,k)+Hz(i,j,k) )
            vdz_old(i,j,k) = v(i,j,k,nnew)*0.5*( Hz(i,j-1,k)+Hz(i,j,k) )
          enddo
        enddo
      enddo

      ! Don't forget 'w' here.

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
      subroutine set_diags_t_4th_h_adv( istr, iend, jstr, jend, k, itrc )  ![
      ! compute the horizontal 4th order advection.
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th


      ![ | T(i-2) | T(i-1) | T(i)   | T(i+1) |
      !                    ^        ^
      !                  FlxU(i)  FlxU(i+1)
      !                   TFX(i)   TFX(i+1)
      !
      ! where UFX and TFX are at u-points and T is at rho-points.
      !
      ! 4th order (centered) fluxes:
      !
      !   TFx(i) = FlxU(i) * ( - T(i-2) + 7T(i-1) + 7T(i) - T(i+1) ) / 12
      !
      ! Hence divergence of the fluxes: (u2rho)
      !
      !   T_adv_4th(i,j) = -pm*pn*( TFX(i+1,j  ) - TFX(i,j)
      !                            +TFY(i  ,j+1) - TFY(i,j) )
      !
      ! However, need to catch boundary and interior masks. We use:
      !
      ! 2nd order centred (cen2(i,j))   = 1/2 ( T(i-1,j) + T(i,j) )
      !
      ! 4th order corrector (cor4(i,j)) = - T(i-2,j) + T(i-1,j) + T(i,j) - T(i+1,j)
      !
      ! thus:
      !
      ! 4th order approx.   = ( 12 cen2 + cor4 * rmask(i-2) * rmask(i+1) ) / 12
      !
      ! (region of no mask) = ( - T(i-2,j) + 7T(i-1,j) + 7T(i,j) - T(i+1,j) ) /12
      !   4th order
      !
      ! (boundary or 1 cell = ( T(i-1,j) + T(i,j)) / 2
      !  away from mask)
      !   2nd order
      !
      !     i-2 i-1  i  i+1
      !      ^   ^   ^   ^          ### = mask
      ! a) |###|   |   |   |
      !                             a) or b) then rmask(i-2) or rmask(i+1) = 0
      ! b) |   |   |   |###|                 so 2nd order centred.
      !
      ! c) |   |###|   |   |
      !                             c) or d) then FlxU(i,j) = 0, so TFx = 0
      !
      ! This relies on rmask set to 0 beyond the physical domain, in init_arrays.
      !
      ! Note: JM and I agreed 2nd order centered scheme at boundary and near mask is
      ! adequate, even though it's possibly not quite what the code does.
      !]

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, k, itrc

      ! local
      integer :: i, j, td

      td=diag_NT_2_td(itrc)
      do j=jstr,jend+1   ! +1 since u2rho averaging for Tdiag(it_h_diss) below
        do i=istr,iend+1

          ! FX4:
          cen2(i,j) = t(i-1,j  ,k,nrhs,itrc) + t(i  ,j  ,k,nrhs,itrc)        ! 1/2 carried in FX4 calc

          FX4(i,j) = -t(i-2,j  ,k,nrhs,itrc) + t(i-1,j  ,k,nrhs,itrc)        ! temp. use FX4 as cor4 array.
     &               -t(i+1,j  ,k,nrhs,itrc) + t(i  ,j  ,k,nrhs,itrc)

          FX4(i,j) = 0.0833333333333 * FlxU(i,j,k)
     &             * ( 6 * cen2(i,j) + FX4(i,j)*rmask(i-2,j)*rmask(i+1,j) )  ! FX4 on RHS is cor4

          ! FY4:
          cen2(i,j) = t(i  ,j-1,k,nrhs,itrc) + t(i  ,j  ,k,nrhs,itrc)

          FY4(i,j) = -t(i  ,j-2,k,nrhs,itrc) + t(i  ,j-1,k,nrhs,itrc)        ! temp. use FY4 as cor4 array.
     &               -t(i  ,j+1,k,nrhs,itrc) + t(i  ,j  ,k,nrhs,itrc)

          FY4(i,j) = 0.0833333333333 * FlxV(i,j,k)
     &             * ( 6 * cen2(i,j) + FY4(i,j)*rmask(i,j-2)*rmask(i,j+1) )  ! FY4 on RHS is cor4

        enddo
      enddo

#ifdef OBC_CHECK
      sum_dbg = 0 ! refresh sum
#endif
      do j=jstr,jend
!        sum_dbg = 0 ! refresh sum
        do i=istr,iend
#if !defined OBC_CHECK
          Tdiag(i,j,k,td,it_h_adv) = -pm(i,j)*pn(i,j)*( FX4(i+1,j  )-FX4(i,j)   ! no *dt since rate of flux divergence
     &                                                 +FY4(i  ,j+1)-FY4(i,j) )
#else
          Tdiag(i,j,k,td,it_h_adv) = - ( FX4(i+1,j  )-FX4(i,j) ! pm(i,j)*pn(i,j)*  ! no *dt since rate of flux divergence
     &                                  +FY4(i  ,j+1)-FY4(i,j) )

          sum_dbg=sum_dbg+Tdiag(i,j,k,td,it_h_adv)
#endif
        enddo
!        if(k==10) print *, 'k=',k,'j=',j,'h_adv sum_dbg = ',sum_dbg
!        if(k==10) print *, 'FX4(istr,j)-(iend+1,j) = ',FX4(istr,j)-FX4(iend+1,j)
      enddo

#ifdef OBC_CHECK
!      if(k>0 .and. k<=11) then
        print *, 'k=',k,'h_adv sum_dbg = ',sum_dbg
!      endif
      if(k==10) then
        print *, 'k==10'
        print *, 'FlxU(istr  ,10,10)=', FlxU(istr  ,10,10)
        print *, 'FlxU(iend+1,10,10)=', FlxU(iend+1,10,10)
        print *, 'FX4(istr  ,10    )=', FX4(istr  ,10)
        print *, 'FX4(istr+1,10    )=', FX4(istr+1,10)
        print *, 'FX4(iend+1,10    )=', FX4(iend+1,10)
        print *, 'FY4(10    ,jstr  )=', FY4(10,jstr  )
        print *, 'FY4(10    ,jend+1)=', FY4(10,jend+1)
      endif
#endif

      end subroutine set_diags_t_4th_h_adv  !]

! ----------------------------------------------------------------------
      subroutine set_diags_t_4th_v_adv( istr, iend, jstr, jend, itrc )  ![
      ! compute the vertical 4th order advection.
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th

      ![
      !  FlxW = We + Wi (explicit & implicit parts of 'w')
      !
      !  FZ4 is at w-points
      !  Tdiag(it_v_adv) is at rho-point (T)
      !
      !  ---- w(i,j,N)   --- (surface)   here w=0
      !
      !       T(i,j,N)
      !
      !  ---- w(i,j,N-1)
      !
      !       .......
      !
      !  ---- w(i,j,1)
      !
      !       T(i,j,1)
      !
      !  ---- w(i,j,0)   --- (bottom)    here w=0
      !
      !]

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, itrc

      ! local
      integer :: i, j, k, td

      td=diag_NT_2_td(itrc)

      do j=jstr,jend   ! FlxW on same vertical axis as tracer, so same loop range.

        ! 2nd order: (due to surface/bottom contraints)
        do k=1,N-1,N-2     ! i.e. do k=1 and k=N-1
          do i=istr,iend
                                                                               ! r-point z=N
            FZ4(i,k) = 0.5 * ( We(i,j,k) + Wi(i,j,k) )                         ! w-point z=N-1  <- FZ4(i,j,N-1)
     &                 * ( t(i,j,k  ,nrhs,itrc) + t(i,j,k+1,nrhs,itrc) )       ! r-point z=N-1

          enddo
        enddo

        ! 4th order:
        do k=2,N-2
          do i=istr,iend
                                                                               ! r-point z=k+2
            FZ4(i,k) = 0.0833333333333 * ( We(i,j,k) + Wi(i,j,k) )             ! r-point z=k+1
     &                 * (-t(i,j,k-1,nrhs,itrc) + 7*t(i,j,k  ,nrhs,itrc)       ! w-point z=k  <- FZ4(i,j,k)
     &                    -t(i,j,k+2,nrhs,itrc) + 7*t(i,j,k+1,nrhs,itrc) )     ! r-point z=k
                                                                               ! r-point z=k-1
          enddo
        enddo

        do k=1,N
          do i=istr,iend
            Tdiag(i,j,k,td,it_v_adv) = - pm(i,j)*pn(i,j)* ( FZ4(i,k) - FZ4(i,k-1) )  ! FZ4 is volume integrated, hence pm*pn to get t*dz?
          enddo
        enddo

      enddo     ! <-- j

      end subroutine set_diags_t_4th_v_adv  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uv_prime  ![

      ![ uf(n+1) = (1-alpha) * uf(n) + alpha * u(n+1)
      !
      !  and set beta = 1-alpha
      !
      !  must be called every timestep for filtering of u & v.
      !  In the notation below u', uf, p' and pf actually have the barotropic mode removed.
      !  E.g. uf + u' = u - ubar
      !
      !  We want high-frequency terms, this should be called up_hf, but currently variable is up_pp.
      !
      !    up_hf = u'p' + uf p' + u' pf  ( where uf/pf is low frequency (filtered) value )
      !
      !  However, pf (p low frequency) is a difficult to get. It's large and affected by breathing grid.
      !  Thus we just use:
      !
      !  up_hf = u'p' + uf p' = ( u' + uf ) p'
      !
      !  since u = u' + uf we can simply get term u'p' = ( u - ubar) p'
      !
      !]

      implicit none

      ! local
      integer :: i, j, k

      do k=1,nz
        do j=1,ny
          do i=1,nx
            uf(i,j,k) = beta * uf(i,j,k) + alpha * ( u(i,j,k,nnew) - ubar(i,j,knew) )
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            vf(i,j,k) = beta * vf(i,j,k) + alpha * ( v(i,j,k,nnew) - vbar(i,j,knew) )
          enddo
        enddo
      enddo

      if( diag_avg .or. mod( iic-ntstart+1 ,rec_rate) == 0        ! +1 since done in step before writing
     &    .or. p_sponge                                    ) then ! also done everytimestep if averaging
                                                                  ! p_sponge needs up_pp at everystep
        do k=1,nz

          do j=1,ny
            do i=1,nx
              up(i,j,k)      = ( u(i,j,k,nnew) - ubar(i,j,knew) ) - uf(i,j,k)  ! need up for div_uvp term
              if (check_prs)  ! if with parameter
     &          u_bcl(i,j,k) = ( u(i,j,k,nnew) - ubar(i,j,knew) )
            enddo
          enddo
          do j=1,ny
            do i=1,nx
              vp(i,j,k) = ( v(i,j,k,nnew) - vbar(i,j,knew) ) - vf(i,j,k)
            enddo
          enddo

          do j=1,ny
            do i=1,nx
              up_pp(i,j,k) =          ( u(i,j,k,nnew) - ubar(i,j,knew) )
     &                       * 0.25 * ( pp(i-1,j,k) + pp(i,j,k) )               ! p and Hz from rho to u
     &                              * ( Hz(i-1,j,k) + Hz(i,j,k) )

              vp_pp(i,j,k) =          ( v(i,j,k,nnew) - vbar(i,j,knew) )
     &                       * 0.25 * ( pp(i,j-1,k) + pp(i,j,k) )
     &                              * ( Hz(i,j-1,k) + Hz(i,j,k) )               ! p and Hz from rho to v

            enddo
          enddo

        enddo     ! <-- k

!       call set_diags_div_uvp

        if(p_sponge)
     &    call set_diags_p_sponge

      endif       ! <-- diag_avg or mod()


      end subroutine set_diags_uv_prime  !]

! ----------------------------------------------------------------------
!     subroutine set_diags_prs_prime(P,ru,rv)  ![
!     ! set p' (the pressure perturbation)
!     !
!     ! Note, pressure (P) is calculated in prsgrd.F. The loop ranges there
!     ! are i=istrU-1,iend and j=jstrV-1,iend
!     ! This allows for averaging from rho to u or v.
!     ! in a sub-domain, u and v points are the first (starting) points
!     ! but rho are the end points, hence the -1 only for the start index.
!     !
!     ! ru = dP/dx * Hz
!     ! rv = dP/dy * Hz

!     use coupling
!     implicit none

!     ! inputs
!     !!!! There is a problem here since istr,iend arent' known
!     real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N), intent(in) :: P, ru, rv

!     ! local
!     integer :: i, j, k

!     if (.not. p_init_done) then

!     !!! JM; the filtered vars need to have permanence
!     !!! JM; the primes should be able to use scratch space

!       !!! figure out how to use work space for these
!       allocate( pxp(GLOBAL_2D_ARRAY, N) )
!       allocate( pyp(GLOBAL_2D_ARRAY, N) )

!       allocate( div_uvp(GLOBAL_2D_ARRAY, N) )



!     if (.not.p_init_done) then
!       p_init_done = .true.
!     endif   

!     !! This needs to run every time-step for proper exponential filtering
!     do k=1,nz
!       do j=0,ny
!         do i=0,nx
!           pf(i,j,k) = beta * pf(i,j,k) + alpha * ( P(i,j,k)  - g*z_w(i,j,N) )
!         enddo
!       enddo
!       do j=1,ny
!         do i=1,nx
!           pxf(i,j,k) = beta*pxf(i,j,k) + alpha*ru(i,j,k)*dxdyi_u(i,j)
!           pyf(i,j,k) = beta*pyf(i,j,k) + alpha*rv(i,j,k)*dxdyi_u(i,j)
!         enddo
!       enddo
!     enddo
!     do j=1,ny
!       do i=1,nx
!         udh_f(i,j) =  beta*udh_f(i,j) + alpha*(
!    &       0.5*pm(i,j)*( ubar(i  ,j,kstp)*(h(i,j)-h(i-1,j))
!    &                   + ubar(i+1,j,kstp)*(h(i+1,j)-h(i,j)) )
!    &      +0.5*pn(i,j)*( vbar(i,j  ,kstp)*(h(i,j)-h(i,j-1))
!    &                   + vbar(i,j+1,kstp)*(h(i,j+1)-h(i,j)) ))
!         dzt_f(i,j) = beta*dzt_f(i,j)  + alpha*pm(i,j)*pn(i,j)*(
!    &                dm_u(i+1,j)*ubar(i+1,j,kstp)*Du_avg1(i+1,j)
!    &              - dm_u(i  ,j)*ubar(i  ,j,kstp)*Du_avg1(i  ,j)
!    &              + dn_v(i,j+1)*vbar(i,j+1,kstp)*Dv_avg1(i,j+1)
!    &              - dn_v(i,j  )*vbar(i,j  ,kstp)*Dv_avg1(i,j  ) )
!       enddo
!     enddo

!     !! Can we replace the mod statement with a 'time_for_his' logical?
!     if( diag_avg .or. mod( iic-ntstart+1 ,rec_rate) == 0            ! +1 since done in step before writing
!    &    .or. p_sponge                                    ) then     ! also done everytimestep if averaging
!                                                                     ! p_sponge needs up_pp at everystep

!       do k=1,nz
!         do j=0,ny ! Need to start at 0 to allow rho2u averaging 
!           do i=0,nx
!             pp(i,j,k) = ( P(i,j,k)  - g*z_w(i,j,N) ) - pf(i,j,k)
!             if (check_prs)                                      
!    &          p_bcl(i,j,k) = ( P(i,j,k)  - g*z_w(i,j,N) )
!           enddo
!         enddo

!         do j=1,ny 
!           do i=1,nx
!             pxp(i,j,k) = ru(i,j,k) - pxf(i,j,k)
!             pyp(i,j,k) = rv(i,j,k) - pyf(i,j,k)
!           enddo
!         enddo
!       enddo

!       !! JM: Pretty sure that this is not neccessary 
!       call set_diags_copy_to_buffer( istr, iend, jstr, jend, pp )

!     endif

!     end subroutine set_diags_prs_prime  !]

! ----------------------------------------------------------------------
      subroutine set_diags_copy_to_buffer( istr, iend, jstr, jend, var )  ![

      implicit none

      integer,                            intent(in)    :: istr, iend, jstr, jend
      real, dimension(GLOBAL_2D_ARRAY,N), intent(inout) :: var

      if(WESTERN_EDGE)  var(istr-1,:,:)=var(istr,:,:)                 ! copy into buffer for averaging u'p' at boundary
      if(EASTERN_EDGE)  var(iend+1,:,:)=var(iend,:,:)                 ! TILING broken
      if(SOUTHERN_EDGE) var(:,jstr-1,:)=var(:,jstr,:)
      if(NORTHERN_EDGE) var(:,jend+1,:)=var(:,jend,:)

      end subroutine set_diags_copy_to_buffer  !]

! ----------------------------------------------------------------------
      subroutine E_baroclinic  ![
      implicit none

      ! Local
      integer :: i,j,k
      real    :: Wbt_f,Wbt_p

      ! Since the exponential smoother is not a great orthogonal
      ! decomposition, we need 3 terms for the hf components:
      !  hf*hf + lf*hf + hf*lf  
      ! lf-> _f, hf -> _p

      ! Barotropic to baroclinic energy flux: Wbt*b
      ! Wbt is the vertical velocity of the barotropic mode
      ! alp = ubar*hx + vbar*hy
      ! bet = z_t + ubar*zx + vbar*zy
      ! Wbt = alp + (zr + h)* (bet - alp)/(zeta + h)
      ! bt2bc = Wbt*b
      if (.not.allocated(bt2bc)) allocate(bt2bc(nx,ny))
      if (.not.allocated(divup)) allocate(divup(nx,ny))

      bt2bc = 0
      do k=1,nz
        do j=1,nx
          do i=1,ny
!           Wbt_f(i,j) = udh_f(i,j) +
!    &       (z_r(i,j,k)+h(i,j))*(dzt_f(i,j)-udh_f(i,j))/(zeta(i,j,knew) + h(i,j))
            Wbt_p = udh_p(i,j) +
     &       (z_r(i,j,k)+h(i,j))*(dzt_p(i,j)-udh_p(i,j))/(zeta(i,j,knew) + h(i,j))
!           bt2bc(i,j) = bt2bc(i,j) + Wbt_p*rho1(i,j,k)*Hz(i,j,k)
          enddo
        enddo
      enddo
      ! Divergence of pressure flux work: u*px + v*py
      divup = 0
      do k=1,nz
        do j=1,ny   
          do i=1,nx
            divup(i,j) = divup(i,j) + 0.5*(
     &              up(i,j,k)*pxp(i,j,k)+up(i+1,j,k)*pxp(i+1,j,k) 
     &            + up(i,j,k)*pxf(i,j,k)+up(i+1,j,k)*pxf(i+1,j,k) 
     &            + uf(i,j,k)*pxp(i,j,k)+uf(i+1,j,k)*pxp(i+1,j,k) 
     &            + vp(i,j,k)*pyp(i,j,k)+vp(i,j+1,k)*pyp(i,j+1,k) 
     &            + vp(i,j,k)*pyf(i,j,k)+vp(i,j+1,k)*pyf(i,j+1,k) 
     &            + vf(i,j,k)*pyp(i,j,k)+vf(i,j+1,k)*pyp(i,j+1,k))
          enddo
        enddo
      enddo

      end subroutine E_baroclinic !]

! ----------------------------------------------------------------------
      subroutine set_diags_p_sponge  ![
      ! set the sponge layer dynamically at each timestep based on
      ! pressure flux in or out of domain.
      ! only needed for s/w/n/e domain boundaries.

      ![The 3 sponge variables are visc2_r (rho-point), diff2 (rho-point) and visc2_p (psi-point).
      ! Since we already chose to put u'p' (up_pp) and v'p' at rho-points (averaged the u/v),
      ! we directly have the sponge switch at the same point as visc2_r and diff2.
      !
      ! visc2_r includes the buffer point, but we ignore the buffer when filtering u'p' since it's not relevant
      ! and we can't average u' into the buffer, so start at first interior point.

      ! Example: visc2_r has (set_nudgcof.F) isp = 4 (0-3) (incl buffer), sponge_size = 3 (1-3).
      !
      !     |
      !  r  u  r  u  r  u  r
      !     |
      !  0  1  1  2  2  3  3
      !
      ! And since up_pp and vp_pp are at rho-point, all our loops are relative to rho-points, starting
      ! from the first interior rho-point (istr) with a count of sponge_size-1=sponge_loop,
      ! need the -1 since the istr is included in the count.
      !
      ! For now, ignore masking as a constant value. This will affact the magnitude of fbox_mup_ppzw
      ! so if we move from binary on_off to something more nuansted based on magnitude and sign, then
      ! will probably need to start handling the mask.
      !
      ! Method:
      ! 1) For each sponge width along boundary:               West Boundary:
      !    Vertically integrate u'p' (up_pp)                    up_pp(i,j,k)       (i.e. u'p')
      !    and get mean across sponge width                         |
      !    (orthongonal to the boundary)                        sum in i and k     (u'p' already contains Hz)
      !    (note up_pp is already multiplied                        |
      !     by Hz in set_diags_uv_prime).                       mup_ppzw(j)        (m=mean, z=vertical integral, w=west)
      !
      ! 2) Do spatial smoothing with (box) average with             |
      !    adjacent points along the boundary.                      |
      !    box_size controls the size of the box.               box_mup_ppzw(j)    (box average)
      !    If box_size=1 then no spatial smoothing.
      !
      ! 3) Time filter result to smooth over random                 |
      !    fluctuation to try capture prevailing flux           fbox_mup_ppzw(j)   (time filter)
      !    direction rather than instantaneous direction.
      !
      ! 4) Use fbox_mup_ppzw to determine if net flux is
      !    in or out of domain, and use binary in or out            |
      !    value (on_off scalar) to be either 1 (out) or        on_off             (binary in or out)
      !    0 (in). Note, this step probably puts a hard
      !    jump in the signal.
      !
      ! 5) Time filter the in=1 or out=0 into the final             |              (sponge multiplier)
      !    sponge multiplier (sponge_fuw).                      sponge_fuw         (f=filtered,u=u,w=west)
      !    e.g. fbox_mup_ppzw -> on_off? -> sponge_fuw
      !
      ! 6) Update sponge value with the time filtered
      !    sponge multiplier sponge_fuw. Use orginal sponge
      !    values from set_nudgcof.F, which are stored here     visc2_r=visc2_r0*sponge_fuw
      !    (e.g. visc2_r(t=0)=visc2_r0).
      !    Need to use original sponge otherwise could never
      !    increase again after decreasing.
      !
      ! 7) Fix corner regions using average of u (sponge_fuw)
      !    and v (sponge_fvs) filtered values since sponge
      !    impacts both u and v signal. Current choice is to
      !    average the mean flux of u and v in the corners.
      !]

      implicit none

      ! local
      integer :: i, j, k, b                   ! b = loop index along size of box averaging
      real    :: avg                          ! temporary average coefficient
      integer :: box_low, box_high            ! upper and lower loop ranges for box (to handle edges where box is smaller)
      real    :: on_off                       ! is flux in or out and hence sponge on or off

      if (.not.west_exchng) then
        ! vertical integral and mean along sponge into domain:
        mup_ppzw=0.0                                         ! reset mean orthongal to boundary
        do j=1,ny
          do k=1,nz                                           ! loop order not contiguous but makes code easier
            do i=1, 1+sponge_loop                       ! i range already not continguous too
              mup_ppzw(j) = mup_ppzw(j) + up_pp(i,j,k)       ! up_pp already multiplied by Hz in set_diags_uv_prime
            enddo
          enddo
          mup_ppzw(j) = mup_ppzw(j) / sponge_size            ! mean value, but really not necessary since sign won't change as per below
        enddo ! <-- j                                        ! finish entire boundary first for box cart average next

        ! Box cart average and time filter:
        call set_diags_p_bry_mpi_exchange_n2s(1,ny, mup_ppzw)     ! exchange between partitions to allow for across bry avg
        do j=1,ny                                         ! set low and high end of loop for box average
          call set_diags_p_box_range_corners_j( j, 1,ny, box_low, box_high)
          box_mup_ppzw(j) = 0                                ! reset box mean (don't need array here could use scalar?)
          do b=box_low, box_high
            box_mup_ppzw(j) = box_mup_ppzw(j) + mup_ppzw(b)  ! calculate box cart average of vertically integrated and orthongal mean
          enddo
          box_mup_ppzw(j) = box_mup_ppzw(j) /
     &                     (box_high-box_low+1)              ! divide for mean
          if (.not. allocated(fbox_mup_ppzw)) then           ! set filter value to box cart average at first step,
            allocate( fbox_mup_ppzw(GLOBAL_1DY_ARRAY) )      ! as can't predict this value. No, since only uses sign below just start at zero.
            fbox_mup_ppzw = 0.0                              ! Move this to _init. Fix for tiling later
          endif
          fbox_mup_ppzw(j) =   beta * fbox_mup_ppzw(j)
     &                      + alpha *  box_mup_ppzw(j)
        enddo   ! <-- j

        ! Use box average filter value to determine sponge on or off, then time filter:
        ! Signal is altered here as net value is turned into binary 1 or 0.
        do j=1,ny
          if(fbox_mup_ppzw(j) < 0.0) then                    ! over-write mean with binary in or out of domain
            on_off = 1.0                                     ! flux out (W) - sponge
          else
            on_off = 0.0                                     ! flux in  (E)
          endif
          sponge_fuw(j) = beta*sponge_fuw(j) + alpha*on_off  ! time filter
        enddo ! <-- j
        sponge_fuw(0)   =sponge_fuw(1 )                  ! extend values into buffer for sponges
        sponge_fuw(ny+1)=sponge_fuw(ny)                  ! to update next.
        ! update sponges:
        do j=1,ny
          do i=1,1+(sponge_loop+1)                  ! here include buffers (i/j/strR) since visc2_r includes buffer
            visc2_r(i,j) = visc2_r0(i,j) * sponge_fuw(j)      ! Doesn't matter that sponge_fuw was calculated without buffer or bry.
            diff2(i,j,:) = diff20(i,j,:) * sponge_fuw(j)
          enddo                                              ! up_pp is at rho-points, so this works
        enddo                                                ! sponge_loop+1 for extra buffer point
      endif  ! < -- western_edge

      if (.not.east_exchng) then
        ! vertical integral and mean along sponge into domain:
        mup_ppze=0.0                                         ! reset mean orthongal to boundary
        do j=1,ny
          do k=1,nz
            do i=nx-sponge_loop,nx
              mup_ppze(j) = mup_ppze(j) + up_pp(i,j,k)
            enddo
          enddo
          mup_ppze(j) = mup_ppze(j) / sponge_size            ! mean value, but really not necessary since sign won't change as per below
        enddo ! <-- j

        ! Box cart average and time filter:
        call set_diags_p_bry_mpi_exchange_n2s(1,ny, mup_ppze)     ! exchange between partitions to allow for across bry avg
        do j=1,ny
          call set_diags_p_box_range_corners_j(j, 1,ny, box_low, box_high)
          box_mup_ppze(j) = 0                                ! reset box mean (don't need array here could use scalar?)
          do b=box_low, box_high
            box_mup_ppze(j) = box_mup_ppze(j) + mup_ppze(b)  ! calculate box cart average of vertically integrated and orthongal mean
          enddo
          box_mup_ppze(j) = box_mup_ppze(j) /
     &                     (box_high-box_low+1)
          if (.not. allocated(fbox_mup_ppze)) then           ! set filter value to box cart average at first step,
            allocate( fbox_mup_ppze(GLOBAL_1DY_ARRAY) )      ! as can't predict this value. No, since only uses sign below just start at zero.
            fbox_mup_ppze = 0.0                              ! Move this to _init. Fix for tiling later
          endif
          fbox_mup_ppze(j) =   beta * fbox_mup_ppze(j)
     &                      + alpha *  box_mup_ppze(j)
        enddo   ! <-- j

        ! Use box average filter value to determine sponge on or off, then time filter:
        ! Signal is altered here as net value is turned into binary 1 or 0.
        do j=1,ny
          if(fbox_mup_ppze(j) > 0.0) then                    ! over-write mean with binary in or out of domain
            on_off = 1.0                                     ! flux out (E) - sponge
          else
            on_off = 0.0                                     ! flux in  (W)
          endif
          sponge_fue(j) = beta*sponge_fue(j) + alpha*on_off  ! time filter
        enddo ! <-- j
        sponge_fue(1   )=sponge_fue(1 )                  ! extend values into buffer for sponge
        sponge_fue(ny+1)=sponge_fue(ny)                  ! update next.
        ! update sponges:
        do j=1,ny
          do i=nx-(sponge_loop+1), nx                  ! Here include buffer (istrR) since visc2_r includes buffer
            visc2_r(i,j) = visc2_r0(i,j) * sponge_fue(j)     ! up_pp is at rho-points, so this works
            diff2(i,j,:) = diff20(i,j,:) * sponge_fue(j)
          enddo                                              ! sponge_loop+1 for extra buffer point
        enddo
      endif  ! < -- eastern_edge

      if (.not.south_exchng) then
        ! vertical integral and mean along sponge into domain:
        mvp_ppzs=0.0                                         ! reset mean orthongal to boundary
        do i=1,nx                                            ! not contiguous loop order but works for algo.
          do k=1,nz
            do j=1,1+sponge_loop
              mvp_ppzs(i) = mvp_ppzs(i) + vp_pp(i,j,k)
            enddo
          enddo
          mvp_ppzs(i) = mvp_ppzs(i) / sponge_size            ! mean value, but really not necessary since sign won't change as per below
        enddo ! <-- i

        ! Box cart average and time filter:
        call set_diags_p_bry_mpi_exchange_e2w(1,nx, mvp_ppzs)
        do i=1,nx
          call set_diags_p_box_range_corners_i(i, 1,nx, box_low, box_high)
          box_mvp_ppzs(i) = 0                                ! reset box mean (don't need array here could use scalar?)
          do b=box_low, box_high
            box_mvp_ppzs(i) = box_mvp_ppzs(i) + mvp_ppzs(b)  ! calculate box cart average of vertically integrated and orthongal mean
          enddo
          box_mvp_ppzs(i) = box_mvp_ppzs(i) /
     &                     (box_high-box_low+1)
          if (.not. allocated(fbox_mvp_ppzs)) then           ! set filter value to box cart average at first step,
            allocate( fbox_mvp_ppzs(GLOBAL_1DX_ARRAY) )      ! as can't predict this value. No, since only uses sign below just start at zero.
            fbox_mvp_ppzs = 0.0                              ! Move this to _init. Fix for tiling later
          endif
          fbox_mvp_ppzs(i) =   beta * fbox_mvp_ppzs(i)
     &                      + alpha *  box_mvp_ppzs(i)
        enddo   ! <-- i

        ! Use box average filter value to determine sponge on or off, then time filter:
        ! Signal is altered here as net value is turned into binary 1 or 0.
        do i=1,nx
          if(fbox_mvp_ppzs(i) < 0.0) then                    ! over-write mean with binary in or out of domain
            on_off = 1.0                                     ! flux out (S) - sponge
          else
            on_off = 0.0                                     ! flux in  (N)
          endif
          sponge_fvs(i) = beta*sponge_fvs(i) + alpha*on_off  ! time filter
        enddo ! <-- i
        sponge_fvs(0) =sponge_fvs(1)                  ! extend values into buffer for sponge
        sponge_fvs(nx+1)=sponge_fvs(nx)                  ! update next.
        ! update sponges:
        do j=1,1+(sponge_loop+1)
          do i=1,nx
            visc2_r(i,j) = visc2_r0(i,j) * sponge_fvs(i)     ! vp_pp is at rho-points, so this works
            diff2(i,j,:) = diff20(i,j,:) * sponge_fvs(i)
          enddo                                              ! sponge_loop+1 for extra buffer point
        enddo
      endif  ! < -- southern_edge

      if (.not.north_exchng) then
        ! vertical integral and mean along sponge into domain:
        mvp_ppzn=0.0                                         ! reset mean orthongal to boundary
        do i=1,nx
          do k=1,nz
            do j=ny-sponge_loop,ny
              mvp_ppzn(i) = mvp_ppzn(i) + vp_pp(i,j,k)
            enddo
          enddo
          mvp_ppzn(i) = mvp_ppzn(i) / sponge_size            ! mean value, but really not necessary since sign won't change as per below
        enddo  !<-- i

        ! Box cart average and time filter:
        call set_diags_p_bry_mpi_exchange_e2w(1,nx, mvp_ppzn)
        do i=1,nx
          call set_diags_p_box_range_corners_i(i, 1,nx, box_low, box_high)
          box_mvp_ppzn(i) = 0                                ! reset box mean (don't need array here could use scalar?)
          do b=box_low, box_high
            box_mvp_ppzn(i) = box_mvp_ppzn(i) + mvp_ppzn(b)  ! calculate box cart average of vertically integrated and orthongal mean
          enddo
          box_mvp_ppzn(i) = box_mvp_ppzn(i) /
     &                     (box_high-box_low+1)
          if (.not. allocated(fbox_mvp_ppzn)) then           ! set filter value to box cart average at first step,
            allocate( fbox_mvp_ppzn(GLOBAL_1DX_ARRAY) )      ! as can't predict this value. No, since only uses sign below just start at zero.
            fbox_mvp_ppzn = 0.0                              ! Move this to _init. Fix for tiling later
          endif
          fbox_mvp_ppzn(i) =   beta * fbox_mvp_ppzn(i)
     &                      + alpha *  box_mvp_ppzn(i)
        enddo   ! <-- i

        ! Use box average filter value to determine sponge on or off, then time filter:
        ! Signal is altered here as net value is turned into binary 1 or 0.
        do i=1,nx
          if(fbox_mvp_ppzn(i) > 0.0) then                    ! over-write mean with binary in or out of domain
            on_off = 1.0                                     ! flux out (N) - sponge
          else
            on_off = 0.0                                     ! flux in  (S)
          endif
          sponge_fvn(i) = beta*sponge_fvn(i) + alpha*on_off  ! time filter
        enddo
        sponge_fvn(0)=sponge_fvn(1 )                  ! extend values into buffer for sponge
        sponge_fvn(nx+1)=sponge_fvn(nx )                  ! update next.
        ! update sponges:
        do j=ny-(sponge_loop+1),ny+1
          do i=1,nx
            visc2_r(i,j) = visc2_r0(i,j) * sponge_fvn(i)     ! vp_pp is at rho-points, so this works
            diff2(i,j,:) = diff20(i,j,:) * sponge_fvn(i)
          enddo                                              ! sponge_loop+1 for extra buffer point
        enddo
      endif  ! < -- northern_edge

      ! update visc2_p interpolating from visc2_r:
      visc2_p=0.0         ! reset for sum
      do j=1,ny           ! viscosity at psi points
        do i=1,nx
          visc2_p(i,j)=visc2_p(i,j) +  0.25 *                  ! equation taken from set_nudgcof.F
     &                 (  visc2_r(i,j)   + visc2_r(i-1,j)      ! where v_ponge*wrk(i,j)=visc2_r(i,j)
     &                  + visc2_r(i,j-1) + visc2_r(i-1,j-1) )  ! simple rho to psi conversion.
        enddo
      enddo

      end subroutine set_diags_p_sponge  !]

! ----------------------------------------------------------------------
      subroutine set_diags_p_box_range_corners_j( j, jstr, jend, box_low, box_high)  ![
      ! set the box averaging along sponge boundary, to handle shortened ends at domain corners

      implicit none

      integer, intent(in)  :: j, jstr, jend
      integer, intent(out) :: box_low, box_high   ! loop range below and above averaging point

      box_low  = j-hbox_size
      box_high = j+hbox_size
      if (SOUTHERN_EDGE) then
        if (j-jstr < hbox_size) box_low  = jstr   ! handle corners with smaller box averaging
      endif
      if (NORTHERN_EDGE) then
        if (jend-j < hbox_size) box_high = jend
      endif
!      if(mynode==0) print *, 'j=',j,' low=',box_low,' high=',box_high

      end subroutine set_diags_p_box_range_corners_j  !]

! ----------------------------------------------------------------------
      subroutine set_diags_p_box_range_corners_i( i, istr, iend, box_low, box_high)  ![
      ! set the box averaging along sponge boundary, to handle shortened ends at domain corners

      implicit none

      integer, intent(in)  :: i, istr, iend
      integer, intent(out) :: box_low, box_high   ! loop range below and above averaging point

      box_low  = i-hbox_size
      box_high = i+hbox_size
      if (WESTERN_EDGE) then
        if (i-istr < hbox_size) box_low  = istr   ! handle corners with smaller box averaging
      endif
      if (EASTERN_EDGE) then
        if (iend-i < hbox_size) box_high = iend
      endif
!      if(mynode==2) print *, 'i=',i,' low=',box_low,' high=',box_high,mynode

      end subroutine set_diags_p_box_range_corners_i  !]

! ----------------------------------------------------------------------
      subroutine set_diags_p_bry_mpi_exchange_n2s(jstr, jend, mup_ppzew)  ![
      ! exchange u'p' integrated and averaged across sponge
      ! to neighbouring MPI sub-domains north and south on domain E and W boundaries only.

      ![use non-blocking communication.
      ! best practice is to start receiving before sending.
      !
      ! can't use MPI_Waitall because some do left right, corners do both.
      ! from internet found MPI 3.0 standard has some grid halo exchange functions
      ! already. This would be good for ROMS but a pain to change over.
      ! didn't find function for boundaries only (i.e. not all neighbours),
      ! so using same scheme Sacha wrote in MPI_exchange8TA.F, which uses MPI_test
      ! and 'comm' variable to confirm all non-blocking communications are complete.
      !
      ! sending index starts from jstr/jend, receiving index starts from jstr-1/jend+1
      ! e.g. send east to west:
      !
      !    jend-1  jend |    x       x      (west)     (before)
      !      x       x  |  jstr    jstr+1   (east)
      !                      |       |
      !    jend-1  jend |  jend+1  jend+2   (west)     (after)
      !
      !]

      use mpi
      use hidden_mpi_vars
      implicit none

      ! inputs
      integer,                                                    intent(in)    :: jstr, jend
      real, dimension( -hbox_size GLOBAL_1DY_ARRAY + hbox_size ), intent(inout) :: mup_ppzew   ! west or east

      ! local
      integer, parameter            :: comm_size=4                                             ! no. of potential communcations
      integer                       :: ierr, status(MPI_STATUS_SIZE), mess_count, c
      integer, dimension(comm_size) :: req                                                     ! communication request handle
      integer, dimension(comm_size) :: comm
      logical                       :: flag

                        ! tags for receive      for send       tags match send to receive (e.g. 1 --> 1)
      do c=1,comm_size  !         x  1  x        x  2  x       up to 4 communications
        comm(c)=0       !         x     x        x     x       send up+down, and receive up+down
      enddo             !         x  2  x        x  1  x       corners only send 1 receive 1 per direction

      if (NORTH_MSG_EXCH) then
        !print *, 'South_recvb4 mup_ppzew(jend+hb) ',mup_ppzew(jend+hbox_size), mynode
        call MPI_Irecv(mup_ppzew(jend+1)        , hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_N, 1, ocean_grid_comm, req(1), ierr)                        ! --> 1 (receive from S)
        comm(1)=1
      endif
      if (SOUTH_MSG_EXCH) then
        !print *, 'North_recvb4 mup_ppzew(jstr-hb)',mup_ppzew(jstr-hbox_size), mynode
        call MPI_Irecv(mup_ppzew(jstr-hbox_size), hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_S, 2, ocean_grid_comm, req(2), ierr)                        ! --> 2 (receive from N)
        comm(2)=2
      endif
      if (SOUTH_MSG_EXCH) then
        !print *, 'South_send   mup_ppzew(jstr+hb-1)',mup_ppzew(jstr+(hbox_size-1)),mynode
        call MPI_Isend(mup_ppzew(jstr), hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_S, 1, ocean_grid_comm, req(3), ierr)                        ! <-- 1 (send from S)
        comm(3)=3
      endif
      if (NORTH_MSG_EXCH) then
        !print *, 'North_send   mup_ppzew(jend-hb)  ',mup_ppzew(jend-(hbox_size-1)),mynode
        call MPI_Isend(mup_ppzew(jend-(hbox_size-1)), hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_N, 2, ocean_grid_comm, req(4), ierr)                        ! <-- 2 (send from N)
        comm(4)=4
      endif

      mess_count=0
      do c=1,comm_size
        if (comm(c) > 0) mess_count=mess_count+1         ! since only exchanging 1D array, exchange was already done
      enddo                                              ! into final array, so don't need to unpack receives like in
                                                         ! mpi_exchange8TA.F. Just use while loop to confirm all complete.
      do while (mess_count>0)
        do c=1,comm_size
          if (comm(c) > 0) then
            call MPI_Test (req(c), flag, status, ierr)   ! make sure all sent messages have been received elsewhere
            if (flag) then                               ! before moving on
               mess_count=mess_count-1 ; comm(c)=0
!               print *, 'done comm & node ',c,mynode
!               if (NORTH_MSG_EXCH .and. c==1) then
!                 print *, 'South_recvA mup_ppzew(jend+hb)',mup_ppzew(jend+hbox_size),mynode  ! debug
!               endif
!               if (SOUTH_MSG_EXCH .and. c==2) then
!                 print *, 'North_recvA mup_ppzew(jstr-hb)',mup_ppzew(jstr-hbox_size),mynode  ! debug
!               endif
            endif
          endif
        enddo
      enddo  ! <-- while

!      if(mynode==2) print*, mynode, 'v(jend-hb-1:jend+hb)    ',mup_ppzew(jend-(hbox_size-1):jend+hbox_size)
!      if(mynode==5) print*, mynode, 'v(jstr-(hb):jstr+(hb-1))',mup_ppzew(jstr-hbox_size:jstr+(hbox_size-1))

      end subroutine set_diags_p_bry_mpi_exchange_n2s  !]

! ----------------------------------------------------------------------
      subroutine set_diags_p_bry_mpi_exchange_e2w(istr, iend, mvp_ppzns)  ![
      ! exchange v'p' integrated and averaged across sponge
      ! to neighbouring MPI sub-domains to the east and west on N and S domain boundaries only.
      ! see description of _n2s version for details on implementation.

      use mpi
      use hidden_mpi_vars
      implicit none

      ! inputs
      integer,                                                    intent(in)    :: istr, iend
      real, dimension( -hbox_size GLOBAL_1DX_ARRAY + hbox_size ), intent(inout) :: mvp_ppzns   ! west or east

      ! local
      integer, parameter            :: comm_size=4                                             ! no. of potential communcations
      integer                       :: ierr, status(MPI_STATUS_SIZE), mess_count, c
      integer, dimension(comm_size) :: req                                                     ! communication request handle
      integer, dimension(comm_size) :: comm
      logical                       :: flag

                        ! tags for receive      for send       tags match send to receive (e.g. 1 --> 1)
      do c=1,comm_size  !         x  x  x        x  x  x       up to 4 communications
        comm(c)=0       !         2     1        1     2       send up+down, and receive up+down
      enddo             !         x  x  x        x  x  x       corners only send 1 receive 1 per direction

      if (EAST_MSG_EXCH) then
        !print *, 'west_recvb4 mvp_ppzns(iend+hb) ',mvp_ppzns(iend+hbox_size), mynode
        call MPI_Irecv(mvp_ppzns(iend+1)        , hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_E, 1, ocean_grid_comm, req(1), ierr)                        ! --> 1 (receive from S)
        comm(1)=1
      endif
      if (WEST_MSG_EXCH) then
        !print *, 'east_recvb4 mvp_ppzns(istr-hb)',mvp_ppzns(istr-hbox_size), mynode
        call MPI_Irecv(mvp_ppzns(istr-hbox_size), hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_W, 2, ocean_grid_comm, req(2), ierr)                        ! --> 2 (receive from N)
        comm(2)=2
      endif
      if (WEST_MSG_EXCH) then
        !print *, 'west_send   mvp_ppzns(istr+hb-1)',mvp_ppzns(istr+(hbox_size-1)),mynode
        call MPI_Isend(mvp_ppzns(istr), hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_W, 1, ocean_grid_comm, req(3), ierr)                        ! <-- 1 (send from S)
        comm(3)=3
      endif
      if (EAST_MSG_EXCH) then
        !print *, 'east_send   mvp_ppzns(iend-hb)  ',mvp_ppzns(iend-(hbox_size-1)),mynode
        call MPI_Isend(mvp_ppzns(iend-(hbox_size-1)), hbox_size, MPI_DOUBLE_PRECISION,
     &                 p_E, 2, ocean_grid_comm, req(4), ierr)                        ! <-- 2 (send from N)
        comm(4)=4
      endif

      mess_count=0
      do c=1,comm_size
        if (comm(c) > 0) mess_count=mess_count+1         ! since only exchanging 1D array, exchange was already done
      enddo                                              ! into final array, so don't need to unpack receives like in
                                                         ! mpi_exchange8TA.F. Just use while loop to confirm all complete.
      do while (mess_count>0)
        do c=1,comm_size
          if (comm(c) > 0) then
            call MPI_Test (req(c), flag, status, ierr)   ! make sure all sent messages have been received elsewhere
            if (flag) then                               ! before moving on
               mess_count=mess_count-1 ; comm(c)=0
!               print *, 'done comm & node ',c,mynode
!               if (EAST_MSG_EXCH .and. c==1) then
!                 print *, 'west_recvA mvp_ppzns(iend+hb)',mvp_ppzns(iend+hbox_size),mynode  ! debug
!               endif
!               if (WEST_MSG_EXCH .and. c==2) then
!                 print *, 'east_recvA mvp_ppzns(istr-hb)',mvp_ppzns(istr-hbox_size),mynode  ! debug
!               endif
            endif
          endif
        enddo
      enddo  ! <-- while

!      if(mynode==1) print*, mynode, 'mvp_ppzns(iend-hb-1:iend+hb)    ',mvp_ppzns(iend-(hbox_size-1):iend+hbox_size)
!      if(mynode==2) print*, mynode, 'mvp_ppzns(istr-(hb):istr+(hb-1))',mvp_ppzns(istr-hbox_size:istr+(hbox_size-1))

      end subroutine set_diags_p_bry_mpi_exchange_e2w  !]

! ----------------------------------------------------------------------
      subroutine set_diags_p_conversion  ![
      ! C(z) = Div( u'p' ) + Dissipation

      implicit none


      end subroutine set_diags_p_conversion  !]

! ----------------------------------------------------------------------
      subroutine set_diags_avg_prs_flx ![
      implicit none

      if (mynode==0) print *, 'avg_prs_flx: ',mynode
      call set_diags_avg_var(up_pp, aup_pp )  ! average up_pp into aup_pp
      call set_diags_avg_var(vp_pp, avp_pp )
      call set_diags_avg_var(div_uvp, adiv_uvp )

      if (check_prs) then
        call set_diags_avg_var(up, aup )
        call set_diags_avg_var(vp, avp )
        call set_diags_avg_var(pp, app )
      endif


      end subroutine set_diags_avg_prs_flx  !]

! ----------------------------------------------------------------------
      subroutine set_diags_avg_var(var, Avar )  ![
      implicit none

      ! inputs
      real, dimension(:,:,:), intent(in)    :: var
      real, dimension(:,:,:), intent(inout) :: Avar


      if (mod(tn,rec_rate)==0) then      ! timestep / rec_rate

        Avar = (Avar + var ) / rec_rate
         
      elseif (mod(tn-1,rec_rate)==0) then   ! this is next step after we have averaged.
               
        Avar = var
         
      else                              ! all other timesteps
               
        Avar = Avar + var
         
      endif   

      end subroutine set_diags_avg_var  !]

! ----------------------------------------------------------------------
      subroutine set_diags_avg ![
      ! Update diagnostics averages
      implicit none

      if (mod(tn-1,rec_rate)==0) then
      ! This is the first step in an averaging cycle
        Udiag_avg = Udiag
        Vdiag_avg = Vdiag
      else
        Udiag_avg = Udiag_avg + Udiag
        Vdiag_avg = Vdiag_avg + Vdiag
      endif

      if (mod(tn,rec_rate)==0) then   ! timestep / rec_rate
        Udiag = Udiag/rec_rate
        Vdiag = Vdiag/rec_rate
      endif

      end subroutine set_diags_avg !]

! ----------------------------------------------------------------------
      subroutine set_diags_t_avg(diagT, td, iTd, itrc)  ![
      ! set the averages of tracer diagnostic variables

      ! NOTE: this won't work for rec_rate=1 because average will never refresh,
      ! not that you would want rec_rate=1 for averages anyway.

      implicit none

      ! inputs
      real, dimension(GLOBAL_2D_ARRAY,N), intent(in) :: diagT      ! 3D part of Tdiag for individual tracer

      integer, intent(in) :: td, iTd, itrc                         ! tn  = timestep
                                                                   ! td  = tracer (local diag index)
      ! local                                                      ! iTd = diagnostic term
      integer           :: i, j, k, ierr                           ! itrc = tracer index as per 't' array
      character(len=30) :: text


      if (mod(tn,rec_rate)==0) then                ! timestep / rec_rate

        do k=1,nz
          do j=1,ny                                ! here we divide the average
            do i=1,nx                              ! by the averaging rate
              Td_avg(i,j,k,td,iTd) = ( Td_avg(i,j,k,td,iTd)+diagT(i,j,k) ) / rec_rate
            enddo
          enddo
        enddo

        text=trim(t_vname(itrc)) / / tdname(iTd)
        call nc_write_var(ncid, Td_avg(:,:,:,td,iTd), nz, text, rp_var, output_rec, ierr)


      else

        if (mod(tn-1,rec_rate)==0) then   ! this is next step after we have averaged.

          do k=1,nz
            do j=1,ny
              do i=1,nx
                Td_avg(i,j,k,td,iTd) = diagT(i,j,k)
              enddo
            enddo
          enddo

        else                              ! all other timesteps

          do k=1,N
            do j=1,ny
              do i=1,nx
                Td_avg(i,j,k,td,iTd) = Td_avg(i,j,k,td,iTd) + diagT(i,j,k)
              enddo
            enddo
          enddo

        endif  ! <-- mod(tn-1,rec_rate)==0

      endif    ! <-- mod(tn  ,rec_rate)==0

      end subroutine set_diags_t_avg  !]

! ----------------------------------------------------------------------
      subroutine wrt_tracer_diags  ![
      ! since the value of the tracer has been stored as dz*t, we still need to
      ! compute the differences for the contribution of each diagnostic term

      ! Unlike surf_flux.F -> wrt_surface_flux_his, here we still need to calculate
      ! averages, and so the 'if (mod(tn,rec_rate)' is split in two parts.
      ! Also, since the value of 't' is kept in Tdiag and not the change in 't',
      ! we use Td_tmp to store the change, and for memory management it is only
      ! one array, and thus needs to be written to file immediately, or added to averges.
      ! In wrt_diagnostics the if(mod(tn,rec_rate)) is split where it is so that the file can
      ! be opened and 'ocean_time' written before writing averages if applicable.

      ! Since the units of 't' change from 'Conc.' to 'dz*Conc.', the diagnostics
      ! can't be calculated in one loop, but must be done individually to adjust for
      ! the correct units of 't'.

      ! We divide by / dt to get diagnostic rate per second.

      use tracers, only: t
      implicit none

      ! local
      integer             :: itrc, td, i, j, k, ierr
      character(len=30)   :: text

      do td=1,nd_trc                                          ! only selected tracers from wrt_t_diag array
        itrc=diag_td_2_NT(td)                                 ! set the actual tracer index in 't'

        Td_tmp = 0 

        ! h_diss:
        do k=1,nz                                             ! first diagnostic calculated using previous step's value
          do j=1,ny
            do i=1,nx
              ! DISP = UP3 - ADV_4th:
              Td_tmp(i,j,k)=( Tdiag(i,j,k,td,it_h_diss) - t(i,j,k,nstp,itrc)*Hz_old(i,j,k) )/dt
     &                      - Tdiag(i,j,k,td,it_h_adv)
              if (check_diag)
     &        T_dsum(i,j,k)=Td_tmp(i,j,k)+Tdiag(i,j,k,td,it_h_adv)  ! add to sum for budget check
            enddo                                                   ! it_h_adv already calculated so add directly (share loop)
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg(Td_tmp, td, it_h_diss, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_h_diss)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! h_adv: (already calculated in set_diags_t_h_4th_adv)
        if (diag_avg) then
          call set_diags_t_avg(Tdiag(:,:,:,td,it_h_adv), td, it_h_adv, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_h_adv)
          call nc_write_var(ncid, Tdiag(:,:,:,td,it_h_adv), N, text, rp_var, output_rec, ierr)
        endif

        ! v_diss:
        do k=1,N
          do j=1,ny
            do i=1,nx
              Td_tmp(i,j,k)= ( Tdiag(i,j,k,td,it_v_diss) - Tdiag(i,j,k,td,it_h_diss) )/dt  ! UP3 stored in it_h_diss
     &                      -  Tdiag(i,j,k,td,it_v_adv)
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k)+Tdiag(i,j,k,td,it_v_adv) ! it_v_adv already calculated so add directly (share loop)
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg(Td_tmp, td, it_v_diss, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_v_diss)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! v_adv: (already calculated in set_diags_t_h_4th_adv)
        if (diag_avg) then
          call set_diags_t_avg(Tdiag(:,:,:,td,it_v_adv),td, it_v_adv, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_v_adv)
          call nc_write_var(ncid, Tdiag(:,:,:,td,it_v_adv), N, text, rp_var, output_rec, ierr)
        endif

        ! From here onwards 't' was no longer in dz*t units:

        ! h_mix:
        do k=1,nz
          do j=1,ny
            do i=1,nx  ! both in dz*t units and Hz is at the same timestep
              Td_tmp(i,j,k)= ( Tdiag(i,j,k,td,it_h_mix)-Tdiag(i,j,k,td,it_v_mix) )*Hz(i,j,k) /dt
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (diag_avg) then
          call set_diags_t_avg(Td_tmp,td, it_h_mix, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_h_mix)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

        ! v_mix and boundary: (swapped order (vmix before hmix in step3d_t) so don't need to reset boundaries of Td_tmp)
        do k=1,nz
          do j=1,ny
            do i=1,nx  ! it_v_mix in 't' units, but it_prev in dz*t units
              Td_tmp(i,j,k)= ( Tdiag(i,j,k,td,it_v_mix)*Hz(i,j,k) - Tdiag(i,j,k,td,it_v_diss) )/dt
              if (check_diag)
     &        T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        if (.not.west_exchng) then                      ! v_mix array also includes boundary update
          i=0
          Td_tmp(i,:,:) = Tdiag(i,:,:,td,it_v_mix)  ! already the difference, so can add it
          T_dsum(i,:,:) = Tdiag(i,:,:,td,it_v_mix)  ! directly to Td_tmp. Only boundary term
        endif                                       ! contributes to sum at the boundary.
        if (.not.east_exchng) then   
          i=nx+1
          Td_tmp(i,:,:) = Tdiag(i,:,:,td,it_v_mix)
          T_dsum(i,:,:) = Tdiag(i,:,:,td,it_v_mix)
        endif
        if (.not.south_exchng) then   
          j=0
          Td_tmp(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
          T_dsum(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
        endif
        if (.not.north_exchng) then   
          j=ny+1
          Td_tmp(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
          T_dsum(:,j,:) = Tdiag(:,j,:,td,it_v_mix)
        endif
        if (diag_avg) then
          call set_diags_t_avg(Td_tmp, td, it_v_mix, itrc)
        else
          text=trim(t_vname(itrc)) / / tdname(it_v_mix)
          call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        endif

      enddo    ! <-- itrc=1,NT

      Hz_old = Hz

      end subroutine wrt_tracer_diags  !]

! ----------------------------------------------------------------------
      subroutine wrt_diagnostics  ![
      ! Write online diagnostics

      use nc_read_write
      implicit none

      ! local
      integer                 :: r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)   ! 2D grid types
#ifdef SOLVE3D
     &                         , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)  ! 3D grid types
#endif
     &                         , p2dgrd(3)
      integer                 :: ierr = 0, imd, i, j, k, tile, idt, itrc, td, np

      character(len=30)       :: text, utext
      integer, dimension(4,2) :: grdtyp3d                      ! u and v type to use in loop
      integer, dimension(2)   :: gt_wrt = (/ up_var, vp_var /) ! grid type for writing
      integer :: dim,diag


      if (mynode==0) print *,'wrt_diag 1',mynode
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

      if (.not. init_done) call init_diagnostics  ! called here as need set_dz_old set above

!     if (mynode==0) print *,'wrt_diag 2',mynode

      diag_icori = .false.  ! flag to avoid computing 'icori' diagnostic twice in compute_horiz_uv_terms.h. Won't compute during pre_step3d but will for step3d_uv1

      if (iic>ntstart) then ! catch to prevent output of zeros at inital step since mod(0,anything) = 0...

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

            ! create psi dimensions:
            p2dgrd(1)=u2dgrd(1) ! xi_u
            p2dgrd(2)=v2dgrd(2) ! eta_v
            p2dgrd(3)=v2dgrd(3) ! time dim

            ! Define remaining variables:

            if (diag_uv) then
              do dim=1,2                                  ! u,v dims

                do diag=1,nmd                             ! diag terms
                  text  = uvw(dim) / / vnames(diag)       ! diag short name
                  call nc_define_var(ncid, text, lname(diag), units, grdtyp3d(:,dim), ierr, prcsn )
                enddo

                if (check_diag .or. tot_change) then      ! following to compare gross change in u/v and sum of diags.
                  text  = uvw(dim) / / check_vnm(1)       ! gross change in dz*u or dz*v over 1 timestep
                  call nc_define_var(ncid, text, check_lnm(1), units, grdtyp3d(:,dim), ierr, prcsn )
                endif

                if (check_diag) then                      ! following to compare gross change in u/v and sum of diags.
                  text  = uvw(dim) / / check_vnm(2)       ! sum of all diag terms over 1 timestep
                  call nc_define_var(ncid, text, check_lnm(2), units, grdtyp3d(:,dim), ierr, prcsn )
                  text  = uvw(dim) / / check_vnm(3)       ! compare gross change against sum of diags. Should be approx zero.
                  call nc_define_var(ncid, text, check_lnm(3), units, grdtyp3d(:,dim), ierr, prcsn )
                endif

                if (check_avg) then
                  text  = uvw(dim) / / check_avnm(1)      ! gross change in u/v up until now since t0
                  call nc_define_var(ncid, text, check_alnm(1), units, grdtyp3d(:,dim), ierr, prcsn )
                  text  = uvw(dim) / / check_avnm(2)      ! total diag avg in u/v up until now since t0
                  call nc_define_var(ncid, text, check_alnm(2), units, grdtyp3d(:,dim), ierr, prcsn )
                  text  = uvw(dim) / / check_avnm(3)      ! total budget close until now since t0
                  call nc_define_var(ncid, text, check_alnm(3), units, grdtyp3d(:,dim), ierr, prcsn )
                endif

              enddo ! <-- dim
            endif   ! <-- diag_uv

            if (diag_trc) then
              do td=1,nd_trc                                    ! only tracers of interest
                itrc=diag_td_2_NT(td)                           ! set the actual tracer index in 't'
                utext = trim(t_units(itrc)) / / tunits
                do idt=1,ndt                                    ! tracer diag terms
                  text  = trim(t_vname(itrc)) / / tdname(idt)
                  call nc_define_var(ncid, text, tlname(idt), utext, r3dgrd, ierr, prcsn )
                enddo

                if (check_diag .or. tot_change) then
                  text  = trim(t_vname(itrc)) / / check_vnm(1)  ! gross change in dz*t over 1 timestep
                  call nc_define_var(ncid, text,  check_lnm(1), utext, r3dgrd, ierr, prcsn )
                endif

                if (check_diag) then
                  text  = trim(t_vname(itrc)) / / check_vnm(2)  ! sum of all diag terms over 1 timestep
                  call nc_define_var(ncid, text,  check_lnm(2), utext, r3dgrd, ierr, prcsn )
                  text  = trim(t_vname(itrc)) / / check_vnm(3)  ! compare gross change against sum of diags. Should be approx zero.
                  call nc_define_var(ncid, text,  check_lnm(3), utext, r3dgrd, ierr, prcsn )
                endif

                if (check_avg) then
                  text  = trim(t_vname(itrc)) / / check_avnm(1) ! gross change in t up until now since t0
                  call nc_define_var(ncid, text, check_alnm(1), utext, r3dgrd, ierr, prcsn )
                  text  = trim(t_vname(itrc)) / / check_avnm(2) ! total diag avg in t up until now since t0
                  call nc_define_var(ncid, text, check_alnm(2), utext, r3dgrd, ierr, prcsn )
                  text  = trim(t_vname(itrc)) / / check_avnm(3) ! total budget close until now since t0
                  call nc_define_var(ncid, text, check_alnm(3), utext, r3dgrd, ierr, prcsn )
                endif

              enddo  ! <-- td
            endif    ! <-- diag_trc

            if(prs_flx) then

              call nc_define_var(ncid, ppnames(1), pplname(1), punits, grdtyp3d(:,1), ierr, prcsn )
              call nc_define_var(ncid, ppnames(2), pplname(2), punits, grdtyp3d(:,2), ierr, prcsn )
              call nc_define_var(ncid, ppnames(3), pplname(3), punits, r3dgrd,        ierr, prcsn )

              if (check_prs) then
                call nc_define_var(ncid, pnames(1), plname(1), pchunits, grdtyp3d(:,1), ierr, prcsn )
                call nc_define_var(ncid, pnames(2), plname(2), pchunits, grdtyp3d(:,2), ierr, prcsn )
                call nc_define_var(ncid, pnames(3), plname(3), pchunits, r3dgrd, ierr, prcsn )
                ! these vars aren't averaged currently
!                text='u baroclin';      if(diag_avg) text= 'Averaged ' / / trim(text) ! catch averaging into long name
                call nc_define_var(ncid, 'u_bcl', 'u baroclin', pchunits, grdtyp3d(:,1), ierr, prcsn )
!                text='p baroclin';      if(diag_avg) text= 'Averaged ' / / trim(text)
                call nc_define_var(ncid, 'p_bcl', 'p baroclin', pchunits, r3dgrd, ierr, prcsn )
!                text='u filtered barc'; if(diag_avg) text= 'Averaged ' / / trim(text)
                call nc_define_var(ncid, 'uf', 'u filtered barc', pchunits, grdtyp3d(:,1), ierr, prcsn )
!                text='p filtered barc'; if(diag_avg) text= 'Averaged ' / / trim(text)
                call nc_define_var(ncid, 'pf', 'p filtered barc', pchunits, r3dgrd, ierr, prcsn )
              endif

              if (p_sponge) then
                if (check_sponge) then
                  call nc_define_var(ncid, 'visc2_r0', 'visc2_r0', 'visc2_r0', r2dgrd, ierr, prcsn )
                  call nc_define_var(ncid, 'visc2_r',  'visc2_r' , 'visc2_r',  r2dgrd, ierr, prcsn )
                  call nc_define_var(ncid, 'diff2',    'diff2' ,   'diff2',    r2dgrd, ierr, prcsn )
                  call nc_define_var(ncid, 'visc2_p',  'visc2_p' , 'visc2_p',  p2dgrd, ierr, prcsn )
!                  call nc_define_var(ncid, 'visc2_p0', 'visc2_p0', 'visc2_p0', p2dgrd, ierr, prcsn )
                endif
              endif  ! <-- p_sponge

            endif    ! <-- prs_flx

            ierr=nf90_enddef(ncid) ! end definition to write changes to disk
            if (mynode == 0) then
              write(*,'(6x,2A)')
     &          'diagnostics :: created new netCDF file ', file_name
            endif

          endif ! <-- (output_rec == 0; needed to create a file)

          print *,'hoi: ',mynode

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

        if (diag_avg) then                                                 ! see wrt_tracer_diags for description of why if(mod(tn,rec_rate)) is split here.
          if (diag_uv)  call set_diags_avg
          if (diag_trc) call wrt_tracer_diags
          if (prs_flx)  call set_diags_avg_prs_flx
        endif                                                              ! Need this after defined file as wrt_tracer will output averages so needs file created.

        if ( mod(tn,rec_rate) == 0 ) then !! Time for an output

          if (diag_uv) then

            do imd=1,nmd        
              if (diag_avg) then
                text = uvw(1) / / vnames(imd)
                call ncwrite(ncid,text,Udiag_avg(:,:,:,imd),output_rec)
                text = uvw(2) / / vnames(imd)
                call ncwrite(ncid,text,Vdiag_avg(:,:,:,imd),output_rec)

!               call nc_write_var(ncid,Udiag_avg(:,:,:,imd),nz,text, gt_wrt(dim), output_rec, ierr)
!               call nc_write_var(ncid,Vdiag_avg(:,:,:,imd),nz,text, gt_wrt(dim), output_rec, ierr)
              else
                text = 'u' // vnames(imd)
                call ncwrite(ncid,text,Udiag(:,:,:,imd),output_rec)
                text = 'v' / / vnames(imd)
                call ncwrite(ncid,text,Vdiag(:,:,:,imd),output_rec)

!               call nc_write_var(ncid,Udiag(:,:,:,imd),nz, text, gt_wrt(dim), output_rec, ierr)
!               call nc_write_var(ncid,Vdiag(:,:,:,imd),nz, text, gt_wrt(dim), output_rec, ierr)
              endif
            enddo

          endif    ! <-- diag_uv

          if(prs_flx) then

            if (diag_avg) then
              call nc_write_var(ncid, aup_pp,   N, ppnames(1), up_var, output_rec, ierr)
              call nc_write_var(ncid, avp_pp,   N, ppnames(2), vp_var, output_rec, ierr)
              call nc_write_var(ncid, adiv_uvp, N, ppnames(3), rp_var, output_rec, ierr)
            else
              call nc_write_var(ncid, up_pp,    N, ppnames(1), up_var, output_rec, ierr)
              call nc_write_var(ncid, vp_pp,    N, ppnames(2), vp_var, output_rec, ierr)
              call nc_write_var(ncid, div_uvp,  N, ppnames(3), rp_var, output_rec, ierr)
            endif

            if (check_prs) then

              if (diag_avg) then
                call nc_write_var(ncid, aup, N, pnames(1), up_var, output_rec, ierr)
                call nc_write_var(ncid, avp, N, pnames(2), vp_var, output_rec, ierr)
                call nc_write_var(ncid, app, N, pnames(3), rp_var, output_rec, ierr)
              else
                call nc_write_var(ncid,  up, N, pnames(1), up_var, output_rec, ierr)
                call nc_write_var(ncid,  vp, N, pnames(2), vp_var, output_rec, ierr)
                call nc_write_var(ncid,  pp, N, pnames(3), rp_var, output_rec, ierr)
              endif

              call nc_write_var(ncid, u_bcl, N, 'u_bcl', up_var, output_rec, ierr)
              call nc_write_var(ncid, p_bcl, N, 'p_bcl', rp_var, output_rec, ierr)
              call nc_write_var(ncid, uf, N, 'uf', up_var, output_rec, ierr)
              call nc_write_var(ncid, pf, N, 'pf', rp_var, output_rec, ierr)
            endif

            if (p_sponge) then
              if (check_sponge) then
                call nc_write_var( ncid, visc2_r0,     1, 'visc2_r0', rp_var, output_rec, ierr )
                call nc_write_var( ncid, visc2_r,      1, 'visc2_r',  rp_var, output_rec, ierr )
                call nc_write_var( ncid, diff2(:,:,1), 1, 'diff2',    rp_var, output_rec, ierr )
                call nc_write_var( ncid, visc2_p,      1, 'visc2_p',  qp_var, output_rec, ierr )
!                call nc_write_var( ncid, visc2_p0,     1, 'visc2_p0', qp_var, output_rec, ierr )
              endif
            endif  ! <-- p_sponge

          endif    ! <-- prs_flx

          if (diag_trc .and. .not. diag_avg)                             ! for averages this has been called already
     &      call wrt_tracer_diags                     ! we call it here non-averages to get the correct writing frequency

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

      endif   ! <-- icc>ntstart

      end subroutine wrt_diagnostics  !]

! ----------------------------------------------------------------------

#endif /* DIAGNOSTICS */

      end module diagnostics
