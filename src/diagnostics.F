      module diagnostics

      ! Set DIAGNOSTICS flag in cppdefs.opt to activate diagnostics.
      ! Define options in diagnostics.opt

#include "cppdefs.opt"

      use dimensions
      use netcdf
      use tracers
      use ocean_vars
      use scalars
      use hidden_mpi_vars
      use nc_read_write
      use grid

      implicit none
      private

#include "diagnostics.opt"

      ![ Potential issues:
      !  - Since we don't have pressure in the boundary buffer, we can't get u'p' at the bry.
      !    We therefore just use the first internal pressure point instead, which introduces
      !    a dx/2 error.
      !]

      ! Preamble:  ![

      logical, public            :: calc_diag=.true.      ! flag to calculate diagnostics in equation system.
      logical                    :: init_done=.false.     ! flag to init diagnostics only once
      integer                    :: navg_diag = 0         ! number of samples in average

      real,dimension(:,:,:),pointer :: wrk_xy

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

      integer,parameter :: ntd = 6   ! number of tracer diag terms
      integer,parameter, public ::
     &  tadvx = 1, ! Advection in x-direction
     &  tadvy = 2, ! Advection in y-direction
     &  tadvz = 3, ! Advection in z-direction
     &  tmixx = 4, ! Mixing in x-direction
     &  tmixy = 5, ! Mixing in y-direction
     &  tmixz = 6  ! Mixing in z-direction

      real,public,allocatable,dimension(:,:,:,:,:) :: Tdiag     ! tracer diagnostic terms
      real,allocatable,dimension(:,:,:,:,:) :: Tdiag_avg ! averages of tracer diagnostic terms
      integer                               :: ntdia     ! number of tracers for which to do diagnostics
      integer, dimension(nt)  :: diag_tidx = -1  ! index of tracer of selected tracers only
      integer,public,dimension(nt)  :: itrc2idia = -1  ! inverse translation from itrc to idia


      character (len=30),  allocatable, dimension(:) :: tdname
      character (len=120), allocatable, dimension(:) :: tlname  ! extra 10 characters for 'Averaged...'

      character (len=20) :: tunits = ' * m/s   (i.e. dC/dt*dz)'

      public set_diags_t_h_mix
      public set_diags_t_4th_h_adv
      public set_diags_t_4th_v_adv

      ! Netcdf outputting:
      ! -----------------
      integer            :: output_time = 0         ! record number of output. 0 indicates we need new file
      integer            :: record = nrpf           ! to trigger the first file creation
      integer            :: ncid=-1, prev_fill_mode

      integer,public     :: varid_uslow, varid_vslow, varid_p_slow

      real, public,allocatable, dimension(:,:) :: dxdyi_u, dxdyi_v,dxdyi ! 1/ surface of cells

      ! Pressure flux terms
      ! --------------
      real,allocatable,dimension(:,:),public :: up,vp         ! u'p' (at upoint) and v'p' (at v_point) 
      real,allocatable,dimension(:,:)    :: up_avg,vp_avg ! averaged u'p' and v'p'
      real,public,allocatable,dimension(:,:,:)  :: p_slow        ! filtered pressure


      real, allocatable, dimension(:,:)  :: p_fast,u_fast,v_fast

      ! Baroclinic energy term
      ! --------------
      real, public,allocatable, dimension(:,:,:) :: u_slow,v_slow
      real, public,allocatable, dimension(:,:,:) :: px_slow,py_slow                       ! will be removed once outputting brought to new system
      real, public,allocatable, dimension(:,:)   :: px_fast,py_fast

      real, allocatable,dimension(:,:) :: upx      ! Divergence of up flux 
      real, allocatable,dimension(:,:) :: vpy      ! Divergence of vp flux 
      real, allocatable,dimension(:,:) :: upx_avg  ! Divergence of up flux
      real, allocatable,dimension(:,:) :: vpy_avg  ! Divergence of up flux
      real, allocatable,dimension(:,:) :: eku      ! Baroclinic Ek
      real, allocatable,dimension(:,:) :: ekv      ! Baroclinic Ek
      real, allocatable,dimension(:,:) :: eku_avg  ! Baroclinic Ek
      real, allocatable,dimension(:,:) :: ekv_avg  ! Baroclinic Ek

      real, allocatable,dimension(:,:) :: bt2bc    ! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bt2bc_avg! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bt2bc_ss    ! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bt2bc_ss_avg! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bt2bc_ff    ! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bt2bc_ff_avg! Barotropic to baroclinic conversion

      real,public, allocatable,dimension(:,:,:) :: rho_d  ! rho for diagnostics
      real, allocatable,dimension(:,:,:) :: Wvl   ! Vertical velocity
      real, allocatable,dimension(:,:,:) :: Wvl_slow   ! Vertical velocity
      real, allocatable,dimension(:,:) :: Wvl_fast  ! Vertical velocity
      real, allocatable,dimension(:,:) :: bc2bc    ! Baroclinic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bc2bc_avg! Baroclinic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bc2bc_ss    ! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bc2bc_ss_avg! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bc2bc_ff    ! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bc2bc_ff_avg! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: ape      ! Available potential energy
      real, allocatable,dimension(:,:) :: ape_avg  ! Available potential energy
      real, allocatable,dimension(:,:) :: drhodz   ! density gradient

      real, public,allocatable, dimension(:,:) :: udh,udh_slow,udh_fast ! ubar.grad(h), filtered and prime
      real, public,allocatable, dimension(:,:) :: dzt,dzt_slow,dzt_fast ! dzeta/dt, filtered and prime
      real, public,allocatable, dimension(:,:) :: Uflx,Vflx
      real, public,allocatable, dimension(:,:) :: uhx, vhy
      real, allocatable,dimension(:,:,:) :: rho_slow
      real, allocatable,dimension(:,:) :: rho_fast

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
      public get_init_slow
      public wrt_rst_diag_pflx
           !]

      contains

! ----------------------------------------------------------------------
      subroutine init_diagnostics ![
      ! Allocate and initialize diagnostic arrays.
      implicit none

      ! local
      integer :: itrc

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
          Udiag_avg = 0
          Vdiag_avg = 0
        endif
      endif

      if (diag_trc) then
        do itrc=1,nt
          ntdia = 0
          if (trc_diag(itrc)) then
            ntdia = ntdia+1           
            diag_tidx(ntdia)= itrc
            itrc2idia(itrc) = ntdia
          endif
        enddo
        allocate(Tdiag(nx,ny,nz,ntd,ntdia) )
        if (diag_avg) then   
          allocate(Tdiag_avg(nx,ny,nz,ntd,ntdia) )
        endif
      endif

      if (diag_pflx) then
        if (mynode==0) print *,'allocating pflux vars'
        if (.not. allocated( p_slow )) then         ! in case already read in with get_init
          allocate(p_slow(GLOBAL_2D_ARRAY,nz) )     ! global 2D needed for exchange in get_init.F
          p_slow = 0
        endif
        allocate(p_fast(0:nx,0:ny) )

!       allocate( up(nx,ny,nz) )
!       allocate( vp(nx,ny,nz) )
        allocate( up(nx,ny) )
        allocate( vp(nx,ny) )
        ! we will always use averaged outputs for the prs fluxes
!       allocate( up_avg(nx,ny,nz) )
!       allocate( vp_avg(nx,ny,nz) )
        allocate( up_avg(nx,ny) )
        allocate( vp_avg(nx,ny) )
        up_avg = 0
        vp_avg = 0
      endif

      if (diag_pflx.or.diag_Ebcln) then
        if (.not. allocated( u_slow )) then         ! in case already read in with get_init
          allocate(u_slow( 1:i1,j0:j1,nz) )
          u_slow = 0
        endif
        if (.not. allocated( v_slow )) then
          allocate(v_slow(i0:i1, 1:j1,nz) )
          v_slow = 0
        endif
        allocate(u_fast(nx,ny) )
        allocate(v_fast(nx,ny) )
      endif

      if (diag_Ebcln) then

        allocate( px_slow(nx,ny,nz) )
        allocate( py_slow(nx,ny,nz) )
        px_slow = 0
        py_slow = 0
        allocate( px_fast(nx,ny) )
        allocate( py_fast(nx,ny) )
        allocate( upx(nx,ny) )
        allocate( vpy(nx,ny) )
        allocate( upx_avg(nx,ny) )
        allocate( vpy_avg(nx,ny) )
        upx_avg = 0
        vpy_avg = 0

        allocate( uhx(nx+1,ny) )
        allocate( vhy(nx,ny+1) )
        allocate( udh(nx,ny) )
        allocate( udh_slow(nx,ny) )
        udh_slow=0
        allocate( udh_fast(nx,ny) )
        allocate( dzt(nx,ny) )
        allocate( dzt_slow(nx,ny) )
        dzt_slow = 0
        allocate( dzt_fast(nx,ny) )
        allocate( Uflx(nx+1,ny) )
        allocate( Vflx(nx,ny+1) )

        allocate( bt2bc(nx,ny) )
        allocate( bt2bc_avg(nx,ny) )
        allocate( bt2bc_ss(nx,ny) )
        allocate( bt2bc_ss_avg(nx,ny) )
        allocate( bt2bc_ff(nx,ny) )
        allocate( bt2bc_ff_avg(nx,ny) )
        bt2bc_avg = 0
        bt2bc_ss_avg = 0
        bt2bc_ff_avg = 0

        allocate( rho_d(nx,ny,nz))
        allocate( rho_slow(nx,ny,nz))
        rho_slow = 0
        allocate( rho_fast(nx,ny))
        allocate( Wvl(GLOBAL_2D_ARRAY,nz))
        allocate( Wvl_slow(nx,ny,nz))
        Wvl_slow = 0;
        allocate( Wvl_fast(nx,ny))

        allocate( bc2bc(nx,ny) )
        allocate( bc2bc_avg(nx,ny) )
        allocate( bc2bc_ss(nx,ny) )
        allocate( bc2bc_ss_avg(nx,ny) )
        allocate( bc2bc_ff(nx,ny) )
        allocate( bc2bc_ff_avg(nx,ny) )
        bc2bc_avg = 0
        bc2bc_ss_avg = 0
        bc2bc_ff_avg = 0

        allocate( ape(nx,ny) )
        allocate( ape_avg(nx,ny) )
        ape_avg = 0

        allocate( drhodz(nx,ny) )

        if (.not.diag_uv) then 
          ! We use them for px/py
          allocate(Udiag(nx,ny,nz,1) )
          allocate(Vdiag(nx,ny,nz,1) )
        endif
      endif

      !!! JM this is a test for reusable space
      allocate(wrk_xy(nx,ny,4))


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


      if (mynode==0) print *,'init diagnostics', init_done
      init_done = .true.  ! ensure init is not triggered again

      end subroutine init_diagnostics  !]

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
      subroutine set_diags_t_h_mix(itrc)  ![
      ! Horizontal mixing from t3dmix_GP.F
      ! zero in interior beyond sponge layer if zero background diffusion!

      use tracers, only: t
      implicit none

      ! inputs/outputs
      integer, intent(in) :: itrc

      ! local
      integer :: i,j,k,td

      td=itrc2idia(itrc)
      do k=1,N
        do j=1,ny
          do i=1,nx
            Tdiag(i,j,k,tmixx,td) = t(i,j,k,nnew,itrc)
          enddo
        enddo
      enddo

      end subroutine set_diags_t_h_mix  !]
! ----------------------------------------------------------------------
      subroutine set_diags_t_4th_h_adv(k,itrc )  ![
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
      integer, intent(in) :: k, itrc

      ! local
      integer :: i, j, td

      td=itrc2idia(itrc)
      do j= 1,ny+1   ! +1 since u2rho averaging for Tdiag(it_h_diss) below
        do i= 1,nx+1

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

      do j= 1,ny
        do i= 1,nx
          Tdiag(i,j,k,tadvx,td) = - ( FX4(i+1,j  )-FX4(i,j) ! pm(i,j)*pn(i,j)*  ! no *dt since rate of flux divergence
     &                                  +FY4(i  ,j+1)-FY4(i,j) )

        enddo
      enddo

      end subroutine set_diags_t_4th_h_adv  !]
! ----------------------------------------------------------------------
      subroutine set_diags_t_4th_v_adv(itrc )  ![
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
      integer, intent(in) :: itrc

      ! local
      integer :: i, j, k, idia

      idia = itrc2idia(itrc)

      do j= 1,ny   ! FlxW on same vertical axis as tracer, so same loop range.

        ! 2nd order: (due to surface/bottom contraints)
        do k=1,N-1,N-2     ! i.e. do k=1 and k=N-1
          do i= 1,nx
                                                                               ! r-point z=N
            FZ4(i,k) = 0.5 * ( We(i,j,k) + Wi(i,j,k) )                         ! w-point z=N-1  <- FZ4(i,j,N-1)
     &                 * ( t(i,j,k  ,nrhs,itrc) + t(i,j,k+1,nrhs,itrc) )       ! r-point z=N-1

          enddo
        enddo

        ! 4th order:
        do k=2,N-2
          do i= 1,nx
                                                                               ! r-point z=k+2
            FZ4(i,k) = 0.0833333333333 * ( We(i,j,k) + Wi(i,j,k) )             ! r-point z=k+1
     &                 * (-t(i,j,k-1,nrhs,itrc) + 7*t(i,j,k  ,nrhs,itrc)       ! w-point z=k  <- FZ4(i,j,k)
     &                    -t(i,j,k+2,nrhs,itrc) + 7*t(i,j,k+1,nrhs,itrc) )     ! r-point z=k
                                                                               ! r-point z=k-1
          enddo
        enddo

        do k=1,nz
          do i= 1,nx
            Tdiag(i,j,k,tadvz,idia) = - pm(i,j)*pn(i,j)* ( FZ4(i,k) - FZ4(i,k-1) )  ! FZ4 is volume integrated, hence pm*pn to get t*dz?
          enddo
        enddo

      enddo     ! <-- j

      end subroutine set_diags_t_4th_v_adv  !]
! ----------------------------------------------------------------------
      subroutine calc_diag_pflux(p) ![
      use coupling
      use param
      use eos_vars
      implicit none
      !! W => kg m2 s-3
      !!  p dz => m3 s-2, u*p*dz => m4 s-3, rho0*u*p*dz=> kg m s-3 = W/m
      !! rho u p dz

      !import/export
      real,dimension(iwest-2:ieast+2,jsouth-2:jnorth+2,nz) :: p

      !local
      integer :: i,j,k
!     real   ,dimension(:,:),pointer :: u_fast,v_fast,p_fast
!     real   ,dimension(:,:),pointer :: px_fast,py_fast

# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      integer :: kbak

      kbak=knew-1
      if (kbak < 1) kbak=4
# endif



      if (diag_pflx) then

        !! JM: maybe a problem because p_fast is 0:nx,0:ny
!       p_fast => wrk_xy(:,:,1)
!       u_fast => wrk_xy(:,:,2)
!       v_fast => wrk_xy(:,:,3)

        up = 0
        vp = 0
        do k=1,nz
          do j=0,ny
            do i=0,nx
              p_slow(i,j,k) = beta*p_slow(i,j,k)
     &            + alpha*( p(i,j,k)-g*z_w(i,j,nz) )
              p_fast(i,j) = p(i,j,k)-g*z_w(i,j,nz)-p_slow(i,j,k)
            enddo
          enddo
          do j=1,ny
            do i=1,nx
              u_slow(i,j,k) = beta*u_slow(i,j,k)
     &            + alpha*( u(i,j,k,nrhs)-ubar(i,j,kbak) )
              v_slow(i,j,k) = beta*v_slow(i,j,k)
     &            + alpha*( v(i,j,k,nrhs)-vbar(i,j,kbak) )
              u_fast(i,j) = u(i,j,k,nrhs)-ubar(i,j,kbak)- u_slow(i,j,k)
              v_fast(i,j) = v(i,j,k,nrhs)-vbar(i,j,kbak)- v_slow(i,j,k)
            enddo
          enddo

          do j=1,ny
            do i=1,nx
            up(i,j) = up(i,j) + u_fast(i,j)
     &           * 0.25*( p_fast(i-1,j)+p_fast(i,j) )*( Hz(i-1,j,k)+Hz(i,j,k) )
            vp(i,j) = vp(i,j) + v_fast(i,j)
     &           * 0.25*( p_fast(i,j-1)+p_fast(i,j) )*( Hz(i,j-1,k)+Hz(i,j,k) )
            enddo
          enddo
        enddo
      endif

      if (diag_Ebcln) then

        ! the U/Vdiag term contains dz*px/py
        px_slow = beta*px_slow + alpha*Udiag(:,:,:,ipgr)
        py_slow = beta*py_slow + alpha*Vdiag(:,:,:,ipgr)

        ! Divergence of pressure flux work: u*px + v*py
        upx = 0
        vpy = 0
        do k=1,nz
          u_fast = u(1:nx,1:ny,k,nrhs)-ubar(1:nx,1:ny,kbak)-u_slow(1:nx,1:ny,k)
          v_fast = v(1:nx,1:ny,k,nrhs)-vbar(1:nx,1:ny,kbak)-v_slow(1:nx,1:ny,k)
          px_fast = Udiag(:,:,k,ipgr) - px_slow(:,:,k)
          py_fast = Vdiag(:,:,k,ipgr) - py_slow(:,:,k)

          upx = upx + u_fast*px_fast
          vpy = vpy + v_fast*py_fast

        enddo

        uhx = ubar(1:nx+1,1:ny,kbak)*(h(1:nx+1,1:ny)-h(0:nx,1:ny))*
     &        0.5*(pm(1:nx+1,1:ny)+pm(0:nx,1:ny))
        vhy = vbar(1:nx,1:ny+1,kbak)*(h(1:nx,1:ny+1)-h(1:nx,0:ny))*
     &        0.5*(pn(1:nx,1:ny+1)+pn(1:nx,0:ny))
        udh = 0.5*(uhx(1:nx,:)+uhx(2:nx+1,:) +
     &             vhy(:,1:ny)+vhy(:,2:ny+1) )

        Uflx =  sum(Flxu(1:nx+1,1:ny,:),dim=3)
        Vflx =  sum(Flxv(1:nx,1:ny+1,:),dim=3)
        dzt  = -dxdyi*( Uflx(2:nx+1,:)-Uflx(1:nx,:) +
     &                  Vflx(:,2:ny+1)-Vflx(:,1:ny) )

        dzt_slow = beta*dzt_slow + alpha*dzt
        udh_slow = beta*udh_slow + alpha*udh
        rho_slow = beta*rho_slow + alpha*rho_d(1:nx,1:ny,:)
        call wvlcty (0,Wvl)
        Wvl_slow = beta*Wvl_slow + alpha*Wvl(1:nx,1:ny,:)

        call E_baroclinic
      endif

      end subroutine calc_diag_pflux  !]
! ----------------------------------------------------------------------
      subroutine E_baroclinic  ![
      use param
      use eos_vars
      implicit none

      ! Local
      integer :: i,j,k,k1,k2
      real    :: Wbt,Wbc,Wbt_f,Wbc_f,Wbt_s,Wbc_s
      integer :: kbak

      kbak=knew-1
      if (kbak < 1) kbak=4

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


      bt2bc = 0
      bc2bc = 0
      bt2bc_ss = 0
      bc2bc_ss = 0
      bt2bc_ff = 0
      bc2bc_ff = 0
!     ape   = 0

      udh_fast = udh - udh_slow
      dzt_fast = dzt - dzt_slow
      do k=1,nz
        rho_fast = rho_d(1:nx,1:ny,k) - rho_slow(:,:,k)
        Wvl_fast = Wvl(1:nx,1:ny,k) - Wvl_slow(:,:,k)
        k2 = min(k+1,nz)
        k1 = max(k-1,1)
        drhodz = -0.5*(rho_slow(1:nx,1:ny,k2)-rho_slow(1:nx,1:ny,k1))/Hz(1:nx,1:ny,k)
        drhodz = drhodz
        do j=1,ny
          do i=1,nx
            if (drhodz(i,j)<1e-10) drhodz(i,j) = 1e-10
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            Wbt = udh(i,j) +
     &       (z_r(i,j,k)+h(i,j))*(dzt(i,j)-udh(i,j))/(zeta(i,j,kbak) + h(i,j))
            Wbt_s = udh_slow(i,j) +
     &       (z_r(i,j,k)+h(i,j))*(dzt_slow(i,j)-udh_slow(i,j))/(zeta(i,j,kbak) + h(i,j))
            Wbt_f = udh_fast(i,j) +
     &       (z_r(i,j,k)+h(i,j))*(dzt_fast(i,j)-udh_fast(i,j))/(zeta(i,j,kbak) + h(i,j))

            Wbc   = Wvl(i,j,k)    - Wbt
            Wbc_s = Wvl_slow(i,j,k) - Wbt_s
            Wbc_f = Wvl_fast(i,j) - Wbt_f

            bt2bc(i,j)    = bt2bc(i,j)    + g*Wbt  *rho_d(i,j,k)    *Hz(i,j,k)
            bt2bc_ss(i,j) = bt2bc_ss(i,j) + g*Wbt_s*rho_slow(i,j,k)*Hz(i,j,k)
            bt2bc_ff(i,j) = bt2bc_ff(i,j) + g*Wbt_f*rho_fast(i,j)  *Hz(i,j,k)
            bc2bc(i,j)    = bc2bc(i,j)    + g*Wbc  *rho_d(i,j,k)    *Hz(i,j,k)
            bc2bc_ss(i,j) = bc2bc_ss(i,j) + g*Wbc_s*rho_slow(i,j,k)*Hz(i,j,k)
            bc2bc_ff(i,j) = bc2bc_ff(i,j) + g*Wbc_f*rho_fast(i,j)  *Hz(i,j,k)
          enddo
        enddo
      enddo

      end subroutine E_baroclinic !]
! ----------------------------------------------------------------------
      subroutine calc_diag_avg ![
      ! Update diagnostics averages
      ! The average is always scaled properly throughout
      ! reset navg_diag=0 after an output of the average
      use param
      implicit none

      ! local
      real :: coef

      navg_diag = navg_diag +1

      coef = 1./navg_diag

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

      if (diag_Ebcln) then
        upx_avg = upx_avg*(1-coef) + upx*coef
        vpy_avg = vpy_avg*(1-coef) + vpy*coef
        bt2bc_avg = bt2bc_avg*(1-coef) + bt2bc*coef
        bt2bc_ss_avg = bt2bc_ss_avg*(1-coef) + bt2bc_ss*coef
        bt2bc_ff_avg = bt2bc_ff_avg*(1-coef) + bt2bc_ff*coef
        bc2bc_avg = bc2bc_avg*(1-coef) + bc2bc*coef
        bc2bc_ss_avg = bc2bc_ss_avg*(1-coef) + bc2bc_ss*coef
        bc2bc_ff_avg = bc2bc_ff_avg*(1-coef) + bc2bc_ff*coef
!       ape_avg = ape_avg*(1-coef) + ape*coef
!       eku_avg = eku_avg*(1-coef) + eku*coef
!       ekv_avg = ekv_avg*(1-coef) + ekv*coef
      endif

      end subroutine calc_diag_avg !]
! ----------------------------------------------------------------------
      subroutine do_diagnostics  ![
      ! as of now, this subroutine is called at the very beginning of step,
      ! before anything has been computed. Therefore, skip the first time
      ! Write diagnostics to file
      use param
      use nc_read_write
      implicit none

      ! local
      integer                 :: ierr = 0, imd, i, j, k, tile, idt, itrc, td, np

      integer :: dim,diag, prev_fill_mode
      character(len=99),save  :: fname

      if (.not. init_done) then
        call init_diagnostics 
        return
      endif

      call calc_diag_avg

      output_time = output_time + dt    

      if (output_time>=output_period) then  ! time for an output
        output_time= 0

        if (record==nrpf) then
          call create_diagnostics_file(fname)
          record = 0
        endif 
        record = record + 1

        ierr=nf90_open(fname,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)  

        call ncwrite(ncid,'ocean_time',(/time/),(/record/))
        if (diag_avg) then
          if (diag_uv)  call wrt_diag_uv(ncid,record,Udiag_avg,Vdiag_avg)
!         if (diag_trc) call wrt_diags_trc(ncid,record,Tdiag_avg)
        else
          if (diag_uv)  call wrt_diag_uv(ncid,record,Udiag,Vdiag)
!         if (diag_trc) call wrt_diags_trc(ncid,record,Tdiag)
        endif

        if (diag_pflx)  call wrt_diag_pflx(ncid,record,up_avg,vp_avg)
        if (diag_Ebcln) call wrt_diag_Ebcln(ncid,record) !,upx_avg,vpy_avg,bt2bc_avg)

        navg_diag = 0

        ierr=nf90_close(ncid)


        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &     'diagnostics :: wrote output, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

      endif ! time for an output


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
      character(len=10),dimension(3) :: dimnames3 ! dimension names
      integer,          dimension(3) :: dimsizes3 ! dim lengths
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
        do it=1,ntdia
          vname = t_vname(diag_tidx(it)) / / '_advx'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','x-advection')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(diag_tidx(it)) / / '_advy'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','y-advection')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(diag_tidx(it)) / / '_advz'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','z-advection')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(diag_tidx(it)) / / '_mixx'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','x-mixing')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(diag_tidx(it)) / / '_mixy'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','y-mixing')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
          vname = t_vname(diag_tidx(it)) / / '_mixz'
          varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
          ierr = nf90_put_att(ncid,varid,'long_name','z-mixing')
          ierr = nf90_put_att(ncid,varid,'units','C m^2/s')
        enddo
      endif

      if (diag_pflx) then
!       dimnames = (/'xi_u','eta_rho','s_rho', 'time'/)
!       dimsizes = (/ xi_u,  eta_rho,  nz,    0/)
!       varid = nccreate(ncid,'up',dimnames,dimsizes,diag_prec)
!       ierr = nf90_put_att(ncid,varid,'long_name','prsflux in x-dir')
!       ierr = nf90_put_att(ncid,varid,'units','W/m^2')

!       dimnames = (/'xi_rho','eta_v','s_rho', 'time'/)
!       dimsizes = (/ xi_rho,  eta_v,  nz,    0/)
!       varid = nccreate(ncid,'vp',dimnames,dimsizes,diag_prec)
!       ierr = nf90_put_att(ncid,varid,'long_name','prsflux in y-dir')
!       ierr = nf90_put_att(ncid,varid,'units','W/m^2')
        dimnames3 = (/'xi_u','eta_rho', 'time'/)
        dimsizes3 = (/ xi_u,  eta_rho,   0/)
        varid = nccreate(ncid,'up',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','prsflux xdir')
        ierr = nf90_put_att(ncid,varid,'units','W/m^2')

        dimnames3 = (/'xi_rho','eta_v', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_v,  0/)
        varid = nccreate(ncid,'vp',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','prsflux ydir')
        ierr = nf90_put_att(ncid,varid,'units','W/m^2')
      endif

      if (diag_Ebcln) then
        dimnames3 = (/'xi_u','eta_rho', 'time'/)
        dimsizes3 = (/ xi_u,  eta_rho,   0/)
        varid = nccreate(ncid,'upx',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','u times dpdx')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

        dimnames3 = (/'xi_rho','eta_v', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_v,  0/)
        varid = nccreate(ncid,'vpy',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','v times dpdy')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

        dimnames3 = (/'xi_rho','eta_rho', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_rho,  0/)
        varid = nccreate(ncid,'bt2bc_al',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','btrop to bclin all')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

        dimnames3 = (/'xi_rho','eta_rho', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_rho,  0/)
        varid = nccreate(ncid,'bt2bc_ss',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','btrop to bclin ss')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

        dimnames3 = (/'xi_rho','eta_rho', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_rho,  0/)
        varid = nccreate(ncid,'bt2bc_ff',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','btrop to bclin ff')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

        dimnames3 = (/'xi_rho','eta_rho', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_rho,  0/)
        varid = nccreate(ncid,'bc2bc_al',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','bclin to bclin all')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

        dimnames3 = (/'xi_rho','eta_rho', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_rho,  0/)
        varid = nccreate(ncid,'bc2bc_ss',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','bclin to bclin ss')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

        dimnames3 = (/'xi_rho','eta_rho', 'time'/)
        dimsizes3 = (/ xi_rho,  eta_rho,  0/)
        varid = nccreate(ncid,'bc2bc_ff',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','bclin to bclin ff')
        ierr = nf90_put_att(ncid,varid,'units','W/m^3')

      endif


      end subroutine create_diagvars  !]
! ----------------------------------------------------------------------
      subroutine wrt_diag_uv(ncid,record,ud,vd)  ![
      ! Write uv diagnostic vars to file

      use param
      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      real, dimension(:,:,:,:),intent(inout) :: ud,vd

      !local
      integer :: bfx,bfy
      integer,   dimension(4) :: start    ! start vector for writing

      ! The diagnostics exclude the domain buffers
      ! Since the diagnostic arrays are (nx,ny), the
      ! boundary domains don't write their first row/column
      ! Technically, we shouldn't write ud(1,:,:.., and vd(:,1,:....
      ! either)
      bfx = 1
      bfy = 1
      if (.not.west_exchng)  bfx = 2
      if (.not.south_exchng) bfy = 2
 
        if (.not.west_exchng) ud(1,:,:,:) = 0
        start = (/1, bfy, 1, record/)

        call ncwrite(ncid,'u_pgr',ud(:,:,:,1),start)
        call ncwrite(ncid,'u_cor',ud(:,:,:,2),start)
        call ncwrite(ncid,'u_adv',ud(:,:,:,3),start)
        call ncwrite(ncid,'u_dis',ud(:,:,:,4),start)
        call ncwrite(ncid,'u_hmx',ud(:,:,:,5),start)
        call ncwrite(ncid,'u_vmx',ud(:,:,:,6),start)
          
        if (.not.south_exchng) vd(:,1,:,:) = 0
        start = (/bfx, 1, 1, record/)

        call ncwrite(ncid,'v_pgr',vd(:,:,:,1),start)
        call ncwrite(ncid,'v_cor',vd(:,:,:,2),start)
        call ncwrite(ncid,'v_adv',vd(:,:,:,3),start)
        call ncwrite(ncid,'v_dis',vd(:,:,:,4),start)
        call ncwrite(ncid,'v_hmx',vd(:,:,:,5),start)
        call ncwrite(ncid,'v_vmx',vd(:,:,:,6),start)

      end subroutine wrt_diag_uv  !]
! ----------------------------------------------------------------------
      subroutine wrt_diag_pflx(ncid,record,up,vp)  ![
      ! Write pressure fluxes

      use param
      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      real, dimension(:,:),intent(inout) :: up,vp
!     real, dimension(:,:,:),intent(inout) :: up,vp

      !local
!     integer,   dimension(4) :: start    ! start vector for writing
      integer,   dimension(3) :: start    ! start vector for writing
      integer :: bfx,bfy

      !! JM The diagnostics exclude the buffers
      bfx = 1
      bfy = 1
      if (.not.west_exchng)  bfx = 2
      if (.not.south_exchng) bfy = 2

!     for 3D up data
!     start = (/bfx, bfy, 1, record/)
!     call ncwrite(ncid,'up',up(bfx:nx,1:ny,:),start)
      start = (/bfx, bfy, record/)
!     print *, 'up max: ',mynode, maxval(up(bfx:nx,1:ny)),maxloc(up(bfx:nx,1:ny))
      call ncwrite(ncid,'up',up(bfx:nx,1:ny),start)
!     stop
          
!     for 3D vp data
!     start = (/bfx, bfy, 1, record/)
!     call ncwrite(ncid,'vp',vp(1:nx,bfy:ny,:),start)
      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'vp',vp(1:nx,bfy:ny),start)


      end subroutine wrt_diag_pflx  !]
! ----------------------------------------------------------------------
      subroutine wrt_diag_Ebcln(ncid,record) !,upx,vpy,bt2bc)  ![
      ! Write pressure fluxes

      use param
      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record

      !local
      integer :: bfx,bfy
      integer,   dimension(3) :: start    ! start vector for writing

      !! JM The diagnostics exclude the buffers
      !! but the boundary domains include them in 
      !! the output so skip one point in the west/south
      !! for the domains that border the boundary
      bfx = 1
      bfy = 1
      if (.not.west_exchng)  bfx = 2
      if (.not.south_exchng) bfy = 2

!     start = (/bfx, bfy, record/)
      start = (/1, bfy, record/)
      call ncwrite(ncid,'upx',upx_avg,start)
          
      start = (/bfx, 1, record/)
      call ncwrite(ncid,'vpy',vpy_avg,start)

      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'bt2bc_al',bt2bc_avg(1:nx,1:ny),start)

      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'bt2bc_ss',bt2bc_ss_avg(1:nx,1:ny),start)

      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'bt2bc_ff',bt2bc_ff_avg(1:nx,1:ny),start)

      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'bc2bc_al',bc2bc_avg(1:nx,1:ny),start)

      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'bc2bc_ss',bc2bc_ss_avg(1:nx,1:ny),start)

      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'bc2bc_ff',bc2bc_ff_avg(1:nx,1:ny),start)


      end subroutine wrt_diag_Ebcln !]
! ----------------------------------------------------------------------
      subroutine wrt_rst_diag_pflx  ![
      ! write the slow components of u, v and P to restart file

      ! variables defined in ocean_vars.F for now. This should be fixed once ocean_vars
      ! variables are moved into seperate module to avoid circular reference.
      ! this is a crude hack to get diag_plfx diagnostics into restart file with current
      ! issue of circular reference. For now the variables are still defined even if not
      ! using diag_pflx with diagnostics. They won't be written though, that is done from
      ! main and only if diag_pflx is turned on.
      implicit none

      ! local
      integer :: ncid, ierr
      integer,dimension(4) :: start

      if (wrt_pflx) then                                    ! this ensures correct timestep to write

        ierr=nf90_open(fname_rst,nf90_write,ncid)

        start=1; start(4)=rec_rst
        call ncwrite(ncid,'u_slow',u_slow( 1:i1,j0:j1,:),start)
        call ncwrite(ncid,'v_slow',v_slow(i0:i1, 1:j1,:),start)
        call ncwrite(ncid,'p_slow',p_slow(i0:i1,j0:j1,:),start)

        ierr=nf90_close(ncid)

      endif

      end subroutine wrt_rst_diag_pflx  !]
! ----------------------------------------------------------------------
      subroutine create_diagnostics_file(fname) ![
      use roms_read_write
      implicit none

      !input/output
      character(len=*),intent(out) :: fname

      !local
      integer :: ncid,ierr

      call create_file('_dia',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      ! Diagnostics options
      ierr = nf90_put_att(ncid,nf90_global,'filtering_time_scale',timescale )

      call create_diagvars(ncid)

      ierr = nf90_enddef(ncid)

      end subroutine create_diagnostics_file !]
! ----------------------------------------------------------------------
      subroutine get_init_slow( ncid, record )  ![
      ! get initial slow u,v and p
      implicit none

      ! input
      integer,intent(in) :: ncid, record

      ! local
      integer :: i,j
      integer :: ierr, varid
      integer,dimension(4) :: start

      start=1; start(4)=record                                       ! 3D vars


      if (mynode==0) print *,'getting u,v,p slow'
      if (.not. allocated(u_slow)) then              ! do_diagnostics not called yet so need to init u_slow,p_slow
        allocate( u_slow(  1:i1,j0:j1, nz ))
        u_slow = 0
      endif
      if (.not. allocated(v_slow)) then
        allocate( v_slow( i0:i1, 1:j1, nz ))
        v_slow = 0
      endif
      if (.not. allocated(p_slow)) then
        allocate( p_slow( GLOBAL_2D_ARRAY, nz ))
        p_slow = 0
      endif

      ierr=nf90_inq_varid (ncid, 'u_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,'u_slow', u_slow,start)
        do j = j0,j1
          do i = 1,i1
            if (umask(i,j)<1) u_slow(i,j,:) = 0
          enddo
        enddo
      else
        u_slow = 0.0
        mpi_master_only write(*,*) ' --- WARNING: '
     &         , 'u_slow'
     &         , ' not in initial file.  Initialized to 0.0'
        ierr=nf90_noerr
      endif

      ierr=nf90_inq_varid (ncid, 'v_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,'v_slow', v_slow,start)
        do j = 1,j1
          do i = i0,i1
            if (vmask(i,j)<1) v_slow(i,j,:) = 0
          enddo
        enddo
      else
        v_slow = 0.0
        mpi_master_only write(*,*) ' --- WARNING: '
     &         , 'v_slow'
     &         , ' not in initial file.  Initialized to 0.0'
        ierr=nf90_noerr
      endif

      ierr=nf90_inq_varid (ncid, 'p_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid, 'p_slow', p_slow(i0:i1,j0:j1,:), start)
        call exchange_tile(1,nx,1,ny, p_slow(START_2D_ARRAY,1),N)  ! need exchange since r2u to get u'p'
        do j = 0,ny
          do i = 0,nx
            if (rmask(i,j)<1) p_slow(i,j,:) = 0
          enddo
        enddo
      else
        p_slow = 0.0
        mpi_master_only write(*,*) ' --- WARNING: '
     &         , 'p_slow'
     &         , ' not in initial file.  Initialized to 0.0'
        ierr=nf90_noerr
      endif

      end subroutine get_init_slow  !]
! ----------------------------------------------------------------------
!      subroutine def_rst_diag_pflx(ncid)  ![
!      implicit none
!
!      ! input
!      integer,intent(in) :: ncid
!
!      integer :: varid, ierr
!
!      varid = nccreate(ncid,'u_slow',
!     &  (/'xi_u','eta_rho','s_rho','time'/),(/xi_u,eta_rho,N,0/))
!      ierr = nf90_put_att(ncid,varid,'long_name','u_slow - fitlered u')
!      ierr = nf90_put_att(ncid,varid,'units','m/s')
!
!      end subroutine  !]
! ----------------------------------------------------------------------

      end module diagnostics
