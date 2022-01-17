      module diagnostics

      ! Set DIAGNOSTICS flag in cppdefs.opt to activate diagnostics.
      ! Define options in diagnostics.opt

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

      ![ Potential issues:
      !  - Since we don't have pressure in the boundary buffer, we can't get u'p' at the bry.
      !    We therefore just use the first internal pressure point instead, which introduces
      !    a dx/2 error.
      !]

      ! Preamble:  ![

      logical, public            :: calc_diag=.true.        ! flag to calculate diagnostics in equation system.
      logical                    :: init_done=.false.       ! flag to init diagnostics only once
      logical                    :: p_init_done=.false.     ! flag to init pressure filter value only once
      integer                    :: navg_diag = 0           ! number of samples in average


      real,dimension(:,:),pointer :: wrku_xy1,wrkv_xy1
      real                        :: diag_time_avg  ! Averages of time

      ! u and v momentum:
      ! ---------------
      integer,parameter :: nmd = 6    ! number of momentum diag terms
      integer,parameter, public ::
     &  ipgr   = 1, ! Hydrostatic pressure gradient    
     &  icori  = 2, ! Coriolis & curvilinear coords  
     &  iadv   = 3, ! purely advective fluxes        
     &  idiss  = 4, ! dissipation from advection      
     &  ihmix  = 5, ! Viscous terms (horiz mixing)    
     &  ivmix  = 6  ! Vertical mixing & 2D/3D coupling 

      real,allocatable,dimension(:,:,:,:),public :: Udiag 
      real,allocatable,dimension(:,:,:,:),public :: Vdiag
      real,allocatable,dimension(:,:,:,:)        :: Udiag_avg ! Averages of diagnostics
      real,allocatable,dimension(:,:,:,:)        :: Vdiag_avg ! Averages of diagnostics
      real, public, allocatable, dimension(:,:)  :: u_prev
      real, public, allocatable, dimension(:,:)  :: v_prev

      real,         allocatable, dimension(:,:)  :: FX4, FY4   ! 4th order advective fluxes to split UP3 = ADV + DISSIP. Called in k loop so only need 2D array for each k
      real,         allocatable, dimension(:,:)  :: cen2       ! 2nd order centred scheme used to handle boundary and mask.
      real,         allocatable, dimension(:,:)  :: FZ4        ! vertical advective fluxes. u/v only. FX4/FY4 both u/v and tracers.

      ! Tracers: 
      ! --------

      integer, public, parameter                          :: ndt = 6   ! number of tracer diag terms
      real, public, allocatable, dimension(:, :, :, :, :) :: Tdiag     ! tracer diagnostic terms
      real,         allocatable, dimension(:, :, :)       :: Td_tmp    ! temporary array to calculate change in tracer diags

      ! if diag_avg
      real,         allocatable, dimension(:, :, :, :, :) :: Tdiag_avg ! averages of tracer diagnostic terms

      integer                            :: nd_trc           ! number of tracers selected tracers of interest (from wrt_t_diag array)
      integer, dimension(:), allocatable :: diag_td_2_NT     ! index of tracer of selected tracers only
      integer, dimension(NT), public     :: diag_NT_2_td=-1  ! convert itrc=1,NT to index of tracer in Tdiag array. Set to bad value -1
      logical, allocatable,dimension(:)  :: idia             ! index to tracer for which to do diagnostics
      integer                            :: ntdiag           ! number of tracers for which to do diagnostics

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
      public set_diags_t_h_mix
      public set_diags_t_4th_h_adv
      public set_diags_t_4th_v_adv

      ! Netcdf outputting:
      ! ------------------
      integer            :: output_rec = 0            ! record number of output. 0 indicates we need new file
      integer            :: total_rec_diag = 0      ! total flux output recordings so far
      integer            :: nrpf_diag = 20   ! number of records per file
      integer            :: tn                        ! roms timestep
      character (len=99) :: file_name = '_dia.nc'     ! suffix to append to simulation output rootname
      integer            :: ncid=-1, prev_fill_mode

      real, public,allocatable, dimension(:,:) :: dxdyi_u, dxdyi_v,dxdyi ! 1/ surface of cells

      ! Pressure flux terms
      ! --------------
      real,public,allocatable,dimension(:,:,:) :: up,vp ! u'p' (at upoint) and v'p' (at v_point) - makes it easier to use in child bry.
      real,       allocatable,dimension(:,:,:) :: p_slow ! filtered pressure 
      real,       allocatable,dimension(:,:)   :: p_fast ! perturbation pressure

      real, allocatable, dimension(:, :, :) :: up_avg,vp_avg  ! averaged u'p' and v'p'
      logical                               :: p_slow_init = .true.

      ! Baroclinic energy term
      ! --------------
      real, public,allocatable, dimension(:,:,:) :: u_slow,v_slow
      real, public,allocatable, dimension(:,:,:) :: u_fast,v_fast
      real, public,allocatable, dimension(:,:,:) :: px_slow,py_slow
      real, public,allocatable, dimension(:,:,:) :: px_fast,py_fast

      real, allocatable, dimension(:,:) :: bt2bc    ! Barotropic to baroclinic conversion
      real, allocatable, dimension(:,:) :: divup    ! Divergence of up flux
      real, public,allocatable, dimension(:,:) :: udh_f, udh_p           ! ubar.grad(h), filtered and prime
      real, public,allocatable, dimension(:,:) :: dzt_f, dzt_p           ! dzeta/dt, filtered and prime
      real, allocatable, dimension(:,:,:) :: div_uvp, adiv_uvp    ! divergence u'p' + v'p' (and average)


      integer :: nprm = 3                      ! number of pressure flux diagnostic terms
      real, public :: alpha, beta              ! filtering constants (needed for child bry)


      character (len=20),  allocatable, dimension(:) :: pnames
      character (len=120), allocatable, dimension(:) :: plname           ! len=70 allows for 8 characters of 'Average ' if diag_avg
      character (len=20),  allocatable, dimension(:), public :: ppnames
      character (len=120), allocatable, dimension(:), public :: pplname

      character (len=3), public :: punits = 'W/m  (i.e. up*dz)'
      character (len=3) :: pchunits = 'TBC'


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

      public do_diagnostics
      public set_diags_u_4th_adv
      public set_diags_v_4th_adv
      public calc_diag_pflux
           !]

      contains

! ----------------------------------------------------------------------
      subroutine init_diagnostics ![
      ! Allocate and initialize diagnostic arrays.
      implicit none

      ! Coefficients for exponential smoothing
      alpha = dt / timescale
      beta  = 1 - alpha

      if (diag_uv) then
        allocate(Udiag(nx,ny,nz,nmd) )
        allocate(Vdiag(nx,ny,nz,nmd) )
        allocate(u_prev(nx,nz) )
        allocate(v_prev(nx,nz) )
        if (diag_avg) then   
          allocate(Udiag_avg(nx,ny,nz,nmd))
          allocate(Vdiag_avg(nx,ny,nz,nmd))
        endif
      endif

      if (diag_pflx) then
        allocate( p_slow(0:nx,0:ny,nz) )
        allocate( p_fast(0:nx,0:ny) )
        allocate( up(nx,ny,nz) )
        allocate( vp(nx,ny,nz) )
        ! we will always use averaged outputs for the prs fluxes
        allocate( up_avg(nx,ny,nz) )
        allocate( vp_avg(nx,ny,nz) )
      endif

      if (diag_Ebcln) then
        allocate( u_slow(nx,ny,nz) )
        allocate( v_slow(nx,ny,nz) )
        allocate( px_slow(nx,ny,nz) )
        allocate( py_slow(nx,ny,nz) )
        allocate( div_uvp(nx,ny,nz) )

        allocate( udh_f(nx,ny) )
        allocate( dzt_f(nx,ny) )
        if (.not.diag_uv) then 
          ! We use them for px/py
          allocate(Udiag(nx,ny,nz,1) )
          allocate(Vdiag(nx,ny,nz,1) )
        endif
      endif

      !!! JM this is a test for reusable space
      allocate(wrku_xy1(nx+1,ny))
      allocate(wrkv_xy1(nx,ny+1))


      !JM: Needed for both diag_uv and diag_pflx

      allocate( dxdyi_u(nx,ny) )
      allocate( dxdyi_v(nx,ny) )
      allocate( dxdyi(nx,ny) )
      dxdyi_u = 0.25*( pn(0:nx-1,1:ny)+pn(1:nx,1:ny) )
     &              *( pm(0:nx-1,1:ny)+pm(1:nx,1:ny) )
      dxdyi_v = 0.25*( pn(1:nx,0:ny-1)+pn(1:nx,1:ny) )
     &              *( pm(1:nx,0:ny-1)+pm(1:nx,1:ny) )
      dxdyi = pn(1:nx,1:ny)*pm(1:nx,1:ny)


      allocate( FX4( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( FY4( GLOBAL_2D_ARRAY ) )       ! called within k loop so only need 2D slice per k
      allocate( FZ4( GLOBAL_1DX_ARRAY, 0:N ) ) ! both u/v/tracer. vertical flux at z_w level.
      FZ4(:,0) = 0  ! top and bottom vertical fluxes are always zero
      FZ4(:,nz)= 0
      allocate(cen2( GLOBAL_2D_ARRAY ) )

      if (diag_trc) call init_diags_trc  

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

      allocate( diag_td_2_NT( nd_trc ) )                     ! to convert indices between 't' and 'Tdiag'
      diag_td_2_NT = -1                                      ! set to bad number (<1) as safeguard

      cnt=0
      do itrc=1,NT
        if (wrt_t_diag(itrc)) then
          cnt=cnt+1
          diag_td_2_NT(cnt)=itrc                             ! store the actual tracer index of 't' array
        endif                                                ! since t(NT) but Tdiag(nd_trc). NT >= nd_trc...
      enddo


      if (diag_avg) then
        allocate(  Tdiag_avg( GLOBAL_2D_ARRAY, N, nd_trc, ndt ) )
        do diag=1,ndt
          tlname(diag)  = 'Average ' / / trim(tlname(diag))              ! and prepend 'Average'
        enddo
      endif

      end subroutine init_diags_trc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_u_4th_adv  ![
      ! Compute the 4th order advection terms for v
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! local
      integer :: i,j,k
      real    :: inv24

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
        ! JM for now, the first and last domain may start/end
        ! JM at something differend than 1 or nx/ny.
        if (.not.west_msg_exch) then 
          i=iwest
          do j=1,ny
            FX4(i,j) = 0.25*( FlxU(i,j,k) + FlxU(i+1,j,k) )
     &                     *( u(i,j,k,nrhs) + u(i+1,j,k,nrhs) ) 
          enddo
        endif
        if(.not.east_msg_exch) then
          i=ieast
          do j=1,ny
            FX4(i,j)= 0.25*( FlxU(i,j,k) + FlxU(i+1,j,k) )
     &                    *( u(i,j,k,nrhs) + u(i+1,j,k,nrhs) )
          enddo
        endif
        if(.not.south_msg_exch) then
          j=jsouth
          do i=1,nx
            FY4(i,j) = 0.25*( FlxV(i-1,j,k) + FlxV(i,j,k) )    
     &                     *( u(i,j,k,nrhs) + u(i,j-1,k,nrhs) )
          enddo
        endif
        if(.not.north_msg_exch) then
          j=jnorth+1
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
      integer :: i,j,k
      real    :: inv24
!     real   ,dimension(:,:),pointer :: Uv,Vv

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
            i=iwest
            FX4(i,j) = 0.25*( FlxU(i,j-1,k)+FlxU(i,j,k) )
     &                     *( v(i-1,j,k,nrhs)+v(i,j,k,nrhs) )
          enddo
        endif
        if(.not.east_msg_exch) then
          i=ieast+1
          do j=1,ny
            FX4(i,j) = 0.25*( FlxU(i,j-1,k)+FlxU(i,j,k) )
     &                     *( v(i-1,j,k,nrhs)+ v(i,j,k,nrhs) )
          enddo
        endif

        if(.not.south_msg_exch) then
          j=jsouth
          do i=1,nx
            FY4(i,j) = 0.25*( FlxV(i,j,k) + FlxV(i,j+1,k) )
     &                     *( v(i,j,k,nrhs) + v(i,j+1,k,nrhs) )
          enddo
        endif
        if(.not.north_msg_exch) then
          j=jnorth
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
     &                        - dxdyi_v(i,j)*( FZ4(i,k)-FZ4(i,k-1) )
          enddo
        enddo

      enddo      ! <-- j

      end subroutine set_diags_v_4th_adv  !]

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
      subroutine calc_diag_pflux(p)
      use coupling
      use param
      use ocean3d
      implicit none

      !import/export
      real,dimension(iwest-2:ieast+2,jsouth-2:jnorth+2,nz) :: p

      !local
      integer :: i,j,k
      real   ,dimension(:,:),pointer :: u_fast,v_fast

      u_fast => wrk_xy1
      v_fast => wrk_xy1

      ! the domain buffers don't have real pressure values
      if (p_slow_init) then
        do k=1,nz
          p_slow(:,:,k) = p(0:nx,0:ny,k) - g*z_w(0:nx,0:ny,nz) 
        enddo
!       p_slow = 0.0
        p_slow_init = .false.
      endif

      if (diag_pflx) then
        do k=1,nz
          p_slow(:,:,k) = beta*p_slow(:,:,k)
     &        + alpha*( p(0:nx,0:ny,k)-g*z_w(0:nx,0:ny,nz) )
        enddo

        ! prs fluxes should always be averaged, so always compute
        ! We think that during the corrector step prsgrd, nrhs points to
        ! the velocities at n+1/2, and they're not multiplied by dz.

        ! There are computations of ubar,vbar, FlxU and FlxV in set_HUV1 plus a
        ! coupling with the barotropic mode. Use those things!!
        do k=1,nz
          p_fast = p(0:nx,0:ny,k)-g*z_w(0:nx,0:ny,nz) - p_slow(:,:,k)
          do j=1,ny
            do i=1,nx
            up(i,j,k) = ( u(i,j,k,nrhs) - ubar(i,j,knew) )
     &           * 0.25*( p_fast(i-1,j)+p_fast(i,j) )*( Hz(i-1,j,k)+Hz(i,j,k) )
            vp(i,j,k) = ( v(i,j,k,nrhs) - vbar(i,j,knew) )
     &           * 0.25*( p_fast(i,j-1)+p_fast(i,j) )*( Hz(i,j-1,k)+Hz(i,j,k) )
            enddo
          enddo
        enddo
      endif
      if (mynode==0) print *,'calc up:  ',up(24,5:7,1)


      if (diag_Ebcln) then
        
        do k=1,nz
          u_slow(:,:,k)  = beta*u_slow(:,:,k) + 
     &                alpha*( u(1:nx+1,1:ny,k,nrhs)-ubar(1:nx+1,1:ny,knew) )
          v_slow(:,:,k)  = beta*v_slow(:,:,k) +
     &                alpha*( v(1:nx,1:ny+1,k,nrhs)-vbar(1:nx,1:ny+1,knew) )
        enddo
        px_slow = beta*px_slow + alpha*Udiag(:,:,:,ipgr)
        py_slow = beta*py_slow + alpha*Vdiag(:,:,:,ipgr)

        ! Divergence of pressure flux work: u*px + v*py
        if (calc_diag) then
          divup_u = 0
          divup_v = 0
          do k=1,nz
            u_fast = u(1:nx+1,1:ny,k,nnew)-ubar(1:nx+1,1:ny)-u_slow(:,:,k)
            v_fast = v(1:nx,1:ny+1,k,nnew)-vbar(1:nx,1:ny+1)-v_slow(:,:,k)
            px_fast = Udiag(:,:,k,ipgs) - px_slow(:,:,k)
            py_fast = Vdiag(:,:,k,ipgs) - py_slow(:,:,k)

            divup_u = divup_u + u_fast*px_fast +
     &                u_fast*px_slow(:,:,k)+u_slow(:,:,k)*px_fast
            divup_v = divup_v + v_fast*py_fast +
     &                v_fast*py_slow(:,:,k)+v_slow(:,:,k)*py_fast
          enddo
          divup = 0.5*( divup_u(1:nx,:) + divup_u(2:nx+1,:)
     &                  divup_v(:,1:ny) + divup_v(:,2:ny+1) )
        endif

        do j=1,ny
          do i=1,nx
!           udh_f(i,j) =  beta*udh_f(i,j) + alpha*(
!    &       0.5*pm(i,j)*( ubar(i  ,j,kstp)*(h(i,j)-h(i-1,j))
!    &                     + ubar(i+1,j,kstp)*(h(i+1,j)-h(i,j)) )
!    &        +0.5*pn(i,j)*( vbar(i,j  ,kstp)*(h(i,j)-h(i,j-1))
!    &                   + vbar(i,j+1,kstp)*(h(i,j+1)-h(i,j)) ))
!             dzt_f(i,j) = beta*dzt_f(i,j)  + alpha*dxdyi(i,j)*(
!    &                dm_u(i+1,j)*ubar(i+1,j,kstp)*Du_avg1(i+1,j)
!    &              - dm_u(i  ,j)*ubar(i  ,j,kstp)*Du_avg1(i  ,j)
!    &              + dn_v(i,j+1)*vbar(i,j+1,kstp)*Dv_avg1(i,j+1)
!    &              - dn_v(i,j  )*vbar(i,j  ,kstp)*Dv_avg1(i,j  ) )
          enddo
        enddo
      endif

      end subroutine calc_diag_pflux  !]

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
     &        u_fast(i,j,k)*px_fast(i,j,k)+u_fast(i+1,j,k)*px_fast(i+1,j,k) 
     &      + u_fast(i,j,k)*px_slow(i,j,k)+u_fast(i+1,j,k)*px_slow(i+1,j,k) 
     &      + u_slow(i,j,k)*px_fast(i,j,k)+u_slow(i+1,j,k)*px_fast(i+1,j,k) 
     &      + v_fast(i,j,k)*py_fast(i,j,k)+v_fast(i,j+1,k)*py_fast(i,j+1,k) 
     &      + v_fast(i,j,k)*py_slow(i,j,k)+v_fast(i,j+1,k)*py_slow(i,j+1,k) 
     &      + v_slow(i,j,k)*py_fast(i,j,k)+v_slow(i,j+1,k)*py_fast(i,j+1,k) )
          enddo
        enddo
      enddo

      end subroutine E_baroclinic !]


! ----------------------------------------------------------------------
      subroutine set_diags_p_conversion  ![
      ! C(z) = Div( u'p' ) + Dissipation

      implicit none


      end subroutine set_diags_p_conversion  !]

! ----------------------------------------------------------------------
      subroutine calc_diag_avg ![
      ! Update diagnostics averages
      ! The average is always scaled properly throughout
      ! reset navg_diag=0 after an output of the average
      implicit none

      ! local
      real :: coef

      navg_diag = navg_diag +1

      coef = 1./navg_diag

      if (mynode==0) print *, 'set_avg :: ',navg_diag,coef

      diag_time_avg = diag_time_avg*(1-coef) + time*coef

      if (diag_avg) then
        if (diag_uv) then
          Udiag_avg = Udiag_avg*(1-coef) + Udiag*coef
          Vdiag_avg = Vdiag_avg*(1-coef) + Vdiag*coef
        endif
        if (diag_trc) then
          Tdiag_avg = Tdiag_avg*(1-coef) + Tdiag*coef
        endif
      endif
      ! always average the prs fluxes
      if (diag_pflx) then
        up_avg = up_avg*(1-coef) + up*coef
        vp_avg = vp_avg*(1-coef) + vp*coef
      endif

      end subroutine calc_diag_avg !]

! ----------------------------------------------------------------------
      subroutine do_diagnostics  ![
      ! Write diagnostics to file
      use param
      use nc_read_write
      implicit none

      ! local
      integer                 :: ierr = 0, imd, i, j, k, tile, idt, itrc, td, np

      integer :: dim,diag, prev_fill_mode
      character(len=99)  :: fname
      integer,save  :: record


      tn = iic-ntstart                                ! model step, as timestep t=n. (icc = t(n+1), and ntstart=1 if not from restart). We write the previous step results during the start of the next step

      if (.not. init_done) call init_diagnostics 


      if (iic>ntstart) then ! catch to prevent output of zeros at inital step since mod(0,anything) = 0...

        call calc_diag_avg

        if ( mod(tn,rec_rate) == 0 ) then !! Time for an output

          if (mod(total_rec_diag,nrpf_diag)==0) then
            call create_diagnostics_file(total_rec_diag,fname)
            record = 0
          endif 
          total_rec_diag  = total_rec_diag +1
          record = record + 1

          ierr=nf90_open(fname,nf90_write,ncid)
          ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)  

!         call ncwrite(ncid,'ocean_time',(/diag_time_avg/),(/record/))
          call ncwrite(ncid,'ocean_time',(/time/),(/record/))
          if (diag_avg) then
            if (diag_uv)  call write_diag_uv(ncid,record,Udiag_avg,Vdiag_avg)
            if (diag_trc) call write_diags_trc(ncid,record,Tdiag_avg)
          else
            if (diag_uv)  call write_diag_uv(ncid,record,Udiag,Vdiag)
            if (diag_trc) call write_diags_trc(ncid,record,Tdiag)
          endif

          if (diag_pflx) call wrt_diag_pflx(ncid,record,up_avg,vp_avg)

          navg_diag = 0

          ierr=nf90_close(ncid)

          if (mynode == 0) then
            write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &       'diagnostics :: wrote output, tdays =', tdays,
     &       'step =', iic-1, 'rec =', record, '/', total_rec_diag
     &        MYID
          endif

        endif ! time for an output

      endif   ! <-- icc>ntstart

      end subroutine do_diagnostics  !]
! ----------------------------------------------------------------------
      subroutine create_diagvars(ncid)  ![
      ! Add diagnostics variables to an opened netcdf file

      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid

      !local
      character(len=20)              :: vname
      character(len=10),dimension(4) :: dimnames ! dimension names
      integer,          dimension(4) :: dimsizes ! dim lengths
      integer                        :: varid,ierr
      integer                        :: it

      if (diag_uv) then
        dimnames = (/'xi_u','eta_rho','s_rho', 'time'/)
        dimsizes = (/ xi_u,  eta_rho,  nz,      0/)

        varid = nccreate(ncid,'u_pgr',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Hydrostatic Presssure Gradient')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_cor',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Coriolis and Curvilinear')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_adv',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Advection (Hor. and Vert.')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_dis',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Numerical dissipation')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_hmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Horizontal mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_vmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Vertical mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        dimnames = (/'xi_rho','eta_v','s_rho', 'time'/)
        dimsizes = (/ xi_rho,  eta_v,  nz,    0/)

        varid = nccreate(ncid,'v_pgr',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Hydrostatic Presssure Gradient')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_cor',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Coriolis and Curvilinear')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_adv',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Advection (Hor. and Vert.')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_dis',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Numerical dissipation')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_hmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Horizontal mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_vmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Vertical mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

      endif

      if (diag_trc) then
        dimnames = (/'xi_rho','eta_rho','s_rho', 'time'/)
        dimsizes = (/ xi_rho,  eta_rho,  nz,    0/)
        do it=1,ntdiag
          vname = t_vname(idia(it)) / / '_advx'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','x-advection')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(idia(it)) / / '_advy'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','y-advection')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(idia(it)) / / '_advz'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','z-advection')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(idia(it)) / / '_mixx'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','x-mixing')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(idia(it)) / / '_mixy'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','x-mixing')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
            
        enddo
      endif

      if (diag_pflx) then
        dimnames = (/'xi_u','eta_rho','s_rho', 'time'/)
        dimsizes = (/ xi_u,  eta_rho,  nz,    0/)
        varid = nccreate(ncid,'up',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','prsflux in x-dir')
        ierr = nf90_put_att(ncid,varid,'units','W/m^2')

        dimnames = (/'xi_rho','eta_v','s_rho', 'time'/)
        dimsizes = (/ xi_rho,  eta_v,  nz,    0/)
        varid = nccreate(ncid,'vp',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','prsflux in y-dir')
        ierr = nf90_put_att(ncid,varid,'units','W/m^2')
      endif


      end subroutine create_diagvars  !]
! ----------------------------------------------------------------------
      subroutine write_diag_uv(ncid,record,ud,vd)  ![
      ! Write uv diagnostic vars to file

      use param
      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      real, dimension(:,:,:,:),intent(inout) :: ud,vd

      !local
      integer :: i0,i1,j0,j1,bfx,bfy
      integer,   dimension(4) :: start    ! start vector for writing

      !! JM showing the new netcdf functionality
      !! JM first and last domains might be smaller
      !! JM i0,i1 etc deals with that for now. 
      !! JM bfx,bfy deals with the fact that only the 
      !! JM boundary have buffers
      i0 = iwest
      i1 = ieast
      j0 = jsouth
      j1 = jnorth
      !! JM The diagnostics exclude the buffers
      bfx = 1
      bfy = 1
      if (.not.west_exchng)  bfx = 2
      if (.not.south_exchng) bfy = 2
 
      if (mynode==2) print *,'wrt_diag',i1
      if (mynode==2) print *,'pgr:   ',Ud(i1,10,10,1),Vd(i1,10,10,1)
      if (mynode==2) print *,'cor:   ',Ud(i1,10,10,2),Vd(i1,10,10,2)
      if (mynode==2) print *,'adv:   ',Ud(i1,10,10,3),Vd(i1,10,10,3)
      if (mynode==2) print *,'dis:   ',Ud(i1,10,10,4),Vd(i1,10,10,4)
      if (mynode==2) print *,'hmx:   ',Ud(i1,10,10,5),Vd(i1,10,10,5)
      if (mynode==2) print *,'vmx:   ',Ud(i1,10,10,6),Vd(i1,10,10,6)

        if (.not.west_exchng) ud(i0,:,:,:) = 0
        start = (/1, bfy, 1, record/)

        call ncwrite(ncid,'u_pgr',Ud(i0:i1,j0:j1,:,1),start)
        call ncwrite(ncid,'u_cor',Ud(i0:i1,j0:j1,:,2),start)
        call ncwrite(ncid,'u_adv',Ud(i0:i1,j0:j1,:,3),start)
        call ncwrite(ncid,'u_dis',Ud(i0:i1,j0:j1,:,4),start)
        call ncwrite(ncid,'u_hmx',Ud(i0:i1,j0:j1,:,5),start)
        call ncwrite(ncid,'u_vmx',Ud(i0:i1,j0:j1,:,6),start)
          
        if (.not.south_exchng) Vdiag_avg(:,j0,:,:) = 0
        start = (/bfx, 1, 1, record/)

        call ncwrite(ncid,'v_pgr',Vd(i0:i1,j0:j1,:,1),start)
        call ncwrite(ncid,'v_cor',Vd(i0:i1,j0:j1,:,2),start)
        call ncwrite(ncid,'v_adv',Vd(i0:i1,j0:j1,:,3),start)
        call ncwrite(ncid,'v_dis',Vd(i0:i1,j0:j1,:,4),start)
        call ncwrite(ncid,'v_hmx',Vd(i0:i1,j0:j1,:,5),start)
        call ncwrite(ncid,'v_vmx',Vd(i0:i1,j0:j1,:,6),start)


      end subroutine write_diag_uv  !]
! ----------------------------------------------------------------------
      subroutine wrt_diag_pflx(ncid,record,up,vp)  ![
      ! Write pressure fluxes

      use param
      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      real, dimension(:,:,:),intent(inout) :: up,vp

      !local
      integer :: i0,i1,j0,j1,bfx,bfy
      integer,   dimension(4) :: start    ! start vector for writing
      integer :: ubnd,vbnd

      i0 = iwest
      i1 = ieast
      j0 = jsouth
      j1 = jnorth
      !! JM The diagnostics exclude the buffers
      bfx = 1
      bfy = 1
      if (.not.west_exchng)  bfx = 2
      if (.not.south_exchng) bfy = 2

      !! write the boundary data for u(v) in the east(north)
      ubnd = 0
      vbnd = 0
      if (.not.east_exchng)  ubnd = 1
      if (.not.north_exchng) vbnd = 1

      start = (/1, bfy, 1, record/)
      call ncwrite(ncid,'up',up(i0:i1+ubnd,j0:j1,:),start)
          
      start = (/bfx, 1, 1, record/)
      call ncwrite(ncid,'vp',vp(i0:i1,j0:j1+vbnd,:),start)

      end subroutine wrt_diag_pflx  !]
! ----------------------------------------------------------------------
      subroutine create_diagnostics_file(total_rec,fname) ![
      use roms_read_write
      use nc_read_write
      implicit none

      !input/output
      integer          ,intent(in)  :: total_rec  ! total netcdf records so far for these variables
      character(len=99),intent(out) :: fname

      !local
      integer :: ncid,ierr,varid
      integer :: prec ! precision of output
      character(len=10),dimension(3) :: dimnames ! dimension names
      integer,          dimension(3) :: dimsizes 

      fname=trim(output_root_name) / / '_dia.nc'
      call create_filename(fname,total_rec)

      ierr=nf90_create(fname,nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'create_diagnostics_file: ')

      varid = nccreate(ncid,'ocean_time',(/'time'/),(/0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','Time since origin')
      ierr = nf90_put_att(ncid,varid,'units','second' )

      call put_global_atts(ncid, ierr)  ! put global attributes in file

      if (diag_uv) then
        call create_diagvars(ncid)
      endif

      ierr = nf90_enddef(ncid)


      if (mynode==0) print *,'diagnostics :: created new netCDF file: ',fname

      end subroutine create_diagnostics_file !]

#endif /* DIAGNOSTICS */

      end module diagnostics
