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
      use roms_read_write
      use grid
      use tides

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
      integer,parameter,public :: nmd = 7    ! number of momentum diag terms
      integer,parameter, public ::
     &  ipgr   = 1, ! Hydrostatic pressure gradient
     &  icori  = 2, ! Coriolis & curvilinear coords
     &  iadv   = 3, ! purely advective fluxes
     &  idiss  = 4, ! dissipation from advection
     &  ihmix  = 5, ! Viscous terms (horiz mixing)
     &  ivmix  = 6, ! Vertical mixing & 2D/3D coupling
     &  icoup  = 7  ! 2D/3D coupling

      real,allocatable,dimension(:,:,:,:),public :: Udiag
      real,allocatable,dimension(:,:,:,:),public :: Vdiag
      real,allocatable,dimension(:,:,:,:)        :: Udiag_avg ! Averages of diagnostics
      real,allocatable,dimension(:,:,:,:)        :: Vdiag_avg ! Averages of diagnostics
      real, public, allocatable, dimension(:,:)  :: u_prev
      real, public, allocatable, dimension(:,:)  :: v_prev


      real,         allocatable, dimension(:,:)  :: FX4, FY4   ! 4th order advective fluxes to split UP3 = ADV + DISSIP. Called in k loop so only need 2D array for each k
      real,         allocatable, dimension(:,:)  :: FZ4        ! vertical advective fluxes. u/v only. FX4/FY4 both u/v and tracers.

      real,allocatable, dimension(:,:) :: grd,flx ! work arrays for 4th order advection

      ! Tracers:
      ! --------

      integer,parameter :: ntd = 7   ! number of tracer diag terms
      integer,parameter, public ::
     &  tadvx = 1, ! Advection in x-direction
     &  tadvy = 2, ! Advection in y-direction
     &  tadvz = 3, ! Advection in z-direction
     &  tmixx = 4, ! Mixing in x-direction
     &  tmixy = 5, ! Mixing in y-direction
     &  tmixz = 6  ! Mixing in z-direction

      real,public,allocatable,dimension(:,:) :: VFlxD     ! Divergence of implicit mixing related vertical fluxes
      real,public,allocatable,dimension(:,:) :: ZFlx      ! Implicit mixing related vertical fluxes
      real,public,allocatable,dimension(:,:,:,:,:) :: Tdiag     ! tracer diagnostic terms
      real,allocatable,dimension(:,:,:,:,:) :: Tdiag_avg ! averages of tracer diagnostic terms
      integer                               :: ntdia     ! number of tracers for which to do diagnostics
      integer,public,dimension(nt)          :: td = 0   ! inverse translation from itrc to idia


      character (len=30),  allocatable, dimension(:) :: tdname
      character (len=120), allocatable, dimension(:) :: tlname  ! extra 10 characters for 'Averaged...'

      character (len=20) :: tunits = ' * m/s   (i.e. dC/dt*dz)'

      public set_diags_t_h_mix
      public diag_t_adv_hc4
      public diag_t_adv_vc4

      ! Netcdf outputting:
      ! -----------------
      real    :: output_time = 0         ! record number of output. 0 indicates we need new file
      integer :: record = nrpf           ! to trigger the first file creation
      integer :: ncid=-1, prev_fill_mode

      real, public,allocatable, dimension(:,:) :: dxdyi_u, dxdyi_v,dxdyi ! 1/ surface of cells

      ! Pressure flux terms
      ! --------------
      real,public,allocatable,dimension(:,:,:) :: p_slow        ! filtered pressure
      real,public,allocatable,dimension(:,:,:) :: u_slow,v_slow

      real,public,allocatable,dimension(:,:) :: p_fast,u_fast,v_fast

      real,allocatable,dimension(:,:),public :: up,vp         ! u'p' (at upoint) and v'p' (at v_point)
      real,allocatable,dimension(:,:)        :: up_avg,vp_avg ! averaged u'p' and v'p'

      real,public,allocatable,dimension(:,:)  :: dzt ! dzeta/dt
      real,public,allocatable,dimension(:,:)  :: Uflx,Vflx
      real,public,allocatable,dimension(:,:,:):: rho_d  ! rho for diagnostics
      real,public,allocatable,dimension(:,:)  :: p_int  ! for barotropic/baroclinic diagnostics

      real, public :: alpha, beta              ! filtering constants (needed for child bry)

      character (len=3), public :: punits = 'W/m  (i.e. up*dz)'
      character (len=3) :: pchunits = 'TBC'


      ! W MOMENTUM: (CURRENTLY NOT WORKING)
      ! -----------

      ! NOTE: u & v are 1:N in z, but w is 0:N hence different variables needed
# ifdef NHMGDIAG
      integer, public, parameter                         :: nwd = 5    ! number of 'w' momentum diag terms (for NHMG runs)
      real,allocatable,dimension(:,:,:,:),public :: Wdiag

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
# endif /* NHMGDIAG */

      public init_diagnostics
      public do_diagnostics
      public set_diags_u_4th_adv
      public set_diags_v_4th_adv
      public calc_diag_pflux
      public get_init_slow
      public wrt_rst_diag_slow
           !]

      contains

! ----------------------------------------------------------------------
      subroutine init_diagnostics ![
      ! Allocate and initialize diagnostic arrays.
      implicit none

      ! local
      integer :: itrc,idx

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
        idx = 0
        do itrc = 1,nt
          if (wrt_t_dia(itrc)) then
            idx = idx+1
            td(itrc) = idx
          endif
        enddo
        ntdia = idx
        allocate(Tdiag(0:nx+1,0:ny+1,nz,ntd,ntdia))
        Tdiag = 0
        allocate(VFlxD(0:nx+1,nz))
        VFlxD = 0
        allocate(ZFlx(0:nx+1,nz))
        ZFlx = 0
        if (diag_avg) then
          allocate(Tdiag_avg(0:nx+1,0:ny+1,nz,ntd,ntdia))
          Tdiag_avg = 0
        endif
      endif

      if (diag_pflx) then
        if (mynode==0) print *,'allocating pflux vars'
        if (.not. allocated( p_slow )) then         ! in case already read in with get_init
          allocate(p_slow(GLOBAL_2D_ARRAY,nz) )     ! global 2D needed for exchange in get_init.F
          p_slow = 0
        endif
        allocate(p_fast(0:nx,0:ny) )

        allocate( up(nx,ny) )
        allocate( vp(nx,ny) )
        allocate( up_avg(nx,ny) )
        allocate( vp_avg(nx,ny) )
        up_avg = 0
        vp_avg = 0

        if (.not. allocated( u_slow )) then         ! in case already read in with get_init
          allocate(u_slow(1:nx+1,0:ny+1,nz) )
          u_slow = 0
        endif
        if (.not. allocated( v_slow )) then
          allocate(v_slow(0:nx+1,1:ny+1,nz) )
          v_slow = 0
        endif
        allocate(u_fast(nx+1,ny) )
        allocate(v_fast(nx,ny+1) )
      endif

      if (diag_pflx) then
        allocate( Uflx(nx+1,ny) )
        allocate( Vflx(nx,ny+1) )

        ! strictly speaking for barotropic/baroclinic balance only
        allocate( dzt(nx,ny) )
      endif

      !!! JM this is a test for reusable space
      ! allocate(wrk_xy(nx,ny,4))


      !JM: Needed for both diag_uv and diag_pflx

      allocate( dxdyi_u(nx,ny) )
      allocate( dxdyi_v(nx,ny) )
      allocate( dxdyi(nx,ny) )
      dxdyi_u = 0.25*( pn(0:nx-1,1:ny)+pn(1:nx,1:ny) )
     &              *( pm(0:nx-1,1:ny)+pm(1:nx,1:ny) )
      dxdyi_v = 0.25*( pn(1:nx,0:ny-1)+pn(1:nx,1:ny) )
     &              *( pm(1:nx,0:ny-1)+pm(1:nx,1:ny) )
      dxdyi = pn(1:nx,1:ny)*pm(1:nx,1:ny)


      allocate( grd( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( flx( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( FX4( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( FY4( GLOBAL_2D_ARRAY ) )       ! called within k loop so only need 2D slice per k
      allocate( FZ4( GLOBAL_1DX_ARRAY, 0:N ) ) ! both u/v/tracer. vertical flux at z_w level.
      FZ4(:,0) = 0  ! top and bottom vertical fluxes are always zero
      FZ4(:,nz)= 0


      if (mynode==0) print *,'init diagnostics', init_done
      init_done = .true.  ! ensure init is not triggered again

      ! Print diagnostics opts (diagnostics.opt) to netcdf attributes
      diagnostic_opt = ''
      call store_string_att(diagnostic_opt, 'Chosen Diagnostics = ')
      if (diag_uv) call store_string_att(diagnostic_opt, 'Momentum')
      if (diag_trc) call store_string_att(diagnostic_opt, ', Tracers')
      if (diag_pflx) call store_string_att(diagnostic_opt,
     & ', Baroclinic Pressure Fluxes')

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
        ! JM at something different than 1 or nx/ny.
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
# ifdef NHMGDIAG
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
# endif /* NHMGDIAG */

! ----------------------------------------------------------------------
      subroutine set_diags_t_h_mix(itrc)  ![
      ! Horizontal mixing from t3dmix_GP.F
      ! zero in interior beyond sponge layer if zero background diffusion!

      use tracers, only: t
      implicit none

      ! inputs/outputs
      integer, intent(in) :: itrc

      ! local
      integer :: i,j,k

      do k=1,N
        do j=1,ny
          do i=1,nx
            Tdiag(i,j,k,tmixx,td) = t(i,j,k,nnew,itrc)
          enddo
        enddo
      enddo

      end subroutine set_diags_t_h_mix  !]
! ----------------------------------------------------------------------
      subroutine diag_t_adv_hc4(k,itrc )  ![
      ! compute the horizontal 4th order advection.
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! inputs
      integer, intent(in) :: k, itrc

      ! local
      integer :: i,j

      do j= 1,ny
        do i= 0,nx+2
          grd(i,j) = (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))*umask(i,j)
        enddo
      enddo
      ! For non-periodic domains, we extrapolate the gradients at the
      ! domain boundaries because the buffer is only 1 point wide.
      if (.not.west_msg_exch) grd(0   ,1:ny) = grd(1   ,1:nx)
      if (.not.east_msg_exch) grd(nx+2,1:ny) = grd(nx+1,1:nx)
      do j= 1,ny
        do i= 0,nx+1
          grd(i,j) = grd(i+1,j)+grd(i,j)
        enddo
        do i= 1,nx+1
          Flx(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc)
     &            -0.1666666666666667*(grd(i,j)-grd(i-1,j))
     &                                              )*FlxU(i,j,k)
          Tdiag(i,j,k,tadvx,td(itrc)) = Flx(i,j)
        enddo
      enddo
      do j= 0,ny+2
        do i= 1,nx
          grd(i,j) = (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))*vmask(i,j)
        enddo
      enddo
      if (.not.north_msg_exch) grd(1:nx,0   ) = grd(1:nx,1   )
      if (.not.south_msg_exch) grd(1:nx,ny+2) = grd(1:nx,ny+1)
      do j= 0,ny+1
        do i= 1,nx
          grd(i,j) = grd(i,j+1)+grd(i,j)
        enddo
      enddo
      do j= 1,ny+1
        do i= 1,nx
          Flx(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc)
     &            -0.1666666666666667*(grd(i,j)-grd(i,j-1))
     &                                              )*FlxV(i,j,k)
          Tdiag(i,j,k,tadvy,td(itrc)) = Flx(i,j)
        enddo
      enddo

      end subroutine diag_t_adv_hc4  !]
! ----------------------------------------------------------------------
      subroutine diag_t_adv_vc4(itrc)  ![
      ! compute the vertical 4th order advection.
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! inputs
      integer,intent(in) :: itrc
      ! local
      integer :: i,j,k

      do j= 1,ny   ! FlxW on same vertical axis as tracer, so same loop range.

        ! 2nd order: (due to surface/bottom contraints)
        do k=1,nz-1,nz-2     ! i.e. do k=1 and k=nz-1
          do i= 1,nx

            FZ4(i,k) = 0.5 * ( We(i,j,k) + Wi(i,j,k) )
     &                 * ( t(i,j,k  ,nrhs,itrc) + t(i,j,k+1,nrhs,itrc) )

          enddo
        enddo

        ! 4th order:
        do k=2,nz-2
          do i= 1,nx

            FZ4(i,k) = 0.0833333333333 * ( We(i,j,k) + Wi(i,j,k) )
     &               * (-t(i,j,k-1,nrhs,itrc) + 7*t(i,j,k  ,nrhs,itrc)
     &                  -t(i,j,k+2,nrhs,itrc) + 7*t(i,j,k+1,nrhs,itrc) )

          enddo
        enddo

        ! fluxes at k=0 and k=nz are always zero
        ! Initialize Tdiag to zero!
        do k=1,nz-1
          do i= 1,nx
            Tdiag(i,j,k,tadvz,td(itrc)) = FZ4(i,k)
          enddo
        enddo

      enddo     ! <-- j

      end subroutine diag_t_adv_vc4  !]
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
      real,dimension(GLOBAL_2D_ARRAY,nz) :: p

      !local
      integer :: i,j,k
      real  :: uint,vint
      real  :: huint,hvint
      real  :: gHzx,gHzy,gHTx,gHty
!     real   ,dimension(:,:),pointer :: u_fast,v_fast,p_fast
!     real   ,dimension(:,:),pointer :: px_fast,py_fast
      real  :: tst

# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      integer :: kbak

      kbak=knew-1
      if (kbak < 1) kbak=4
# endif

      ! There's an issue with WEC stokes drift
      ! It's added to flxU,flxV, but not to ubar

      if (diag_pflx) then

        !! JM: maybe a problem because p_fast is 0:nx,0:ny
!       p_fast => wrk_xy(:,:,1)
!       u_fast => wrk_xy(:,:,2)
!       v_fast => wrk_xy(:,:,3)

        ! WEC terms can modify the match between
        ! integrated flux and ubar,vbar

        Uflx =  sum(Flxu(1:nx+1,1:ny,:),dim=3)
        Vflx =  sum(Flxv(1:nx,1:ny+1,:),dim=3)

        ! Strictly speaking only neccesary for barotropic/baroclinic
        ! balance
        dzt  = -dxdyi*( Uflx(2:nx+1,:)-Uflx(1:nx,:) +
     &                  Vflx(:,2:ny+1)-Vflx(:,1:ny) )

        ! Divide the barotropic fluxes into velocities [m/s]
        do j=1,ny
          do i=1,nx+1
            Uflx(i,j) =  Uflx(i,j)*(pn(i,j)+pn(i-1,j))/
     &               (z_w(i,j,nz)-z_w(i,j,0)+z_w(i-1,j,nz)-z_w(i-1,j,0))
          enddo
        enddo
        do j=1,ny+1
          do i=1,nx
            Vflx(i,j) =  Vflx(i,j)*(pm(i,j)+pm(i,j-1))/
     &               (z_w(i,j,nz)-z_w(i,j,0)+z_w(i,j-1,nz)-z_w(i,j-1,0))
          enddo
        enddo

        up = 0
        vp = 0
        tst = 0
        do k=1,nz
          do j=0,ny
            do i=0,nx
              p_slow(i,j,k) = beta*p_slow(i,j,k) + alpha*p(i,j,k)
              p_fast(i,j)   = p(i,j,k) - p_slow(i,j,k)
            enddo
          enddo
          ! u_slow, and u_fast are baroclinic velocities
          do j=1,ny
            do i=1,nx+1
              u_slow(i,j,k) = beta*u_slow(i,j,k)
     &            + alpha*( u(i,j,k,nrhs)-Uflx(i,j))
              u_fast(i,j) = u(i,j,k,nrhs)-Uflx(i,j)- u_slow(i,j,k)
            enddo
          enddo
          ! v_slow, and v_fast are baroclinic velocities
          do j=1,ny+1
            do i=1,nx
              v_slow(i,j,k) = beta*v_slow(i,j,k)
     &            + alpha*( v(i,j,k,nrhs)-Vflx(i,j) )
              v_fast(i,j) = v(i,j,k,nrhs)-Vflx(i,j)- v_slow(i,j,k)
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
      endif  ! if pflux

      end subroutine calc_diag_pflux  !]
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
      integer :: ierr = 0
      integer :: dim,diag, prev_fill_mode
      character(len=99),save  :: fname
      character(len=100) :: output_time_string
      call calc_diag_avg

      output_time = output_time + dt

      if (output_time>=output_period) then  ! time for an output

        ! store averaging timescale in diagnostic output file
        if (diag_avg) then
          diag_avg_output=''
          write(output_time_string,'(F12.1, A)') output_time, ' seconds'
          call store_string_att(diag_avg_output,output_time_string)
        endif

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
          if (diag_trc) call wrt_diag_trc(ncid,record,Tdiag_avg)
        else
          if (diag_uv)  call wrt_diag_uv(ncid,record,Udiag,Vdiag)
          if (diag_trc) call wrt_diag_trc(ncid,record,Tdiag)
        endif

        if (diag_pflx)  call wrt_diag_pflx(ncid,record,up_avg,vp_avg)

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
        dimnames = (/dn_xu,dn_yr,dn_zr, dn_tm/)
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

        varid = nccreate(ncid,'u_cpl',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'2D/3D coupling')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        dimnames = (/dn_xr,dn_yv,dn_zr, dn_tm/)
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

        varid = nccreate(ncid,'v_cpl',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'2D/3D coupling')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

      endif

      if (diag_trc) then
        dimnames = (/dn_xr,dn_yr,dn_zr, dn_tm/)
        dimsizes = (/ xi_rho,  eta_rho,  nz,    0/)
        do it=1,nt
          if (wrt_t_dia(it)) then
           vname = trim(t_vname(it)) // '_advx'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','x-advective flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_advy'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','y-advective flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_advz'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','z-advective flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_mixx'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','mixing x-flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_mixy'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','mixing y-flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_mixz'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','mixing z-flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
          endif
        enddo
      endif

      if (diag_pflx) then
        dimnames3 = (/dn_xu,dn_yr, dn_tm/)
        dimsizes3 = (/ xi_u,  eta_rho,   0/)
        varid = nccreate(ncid,'up',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','prsflux xdir')
        ierr = nf90_put_att(ncid,varid,'units','W/m^2')

        dimnames3 = (/dn_xr,dn_yv, dn_tm/)
        dimsizes3 = (/ xi_rho,  eta_v,  0/)
        varid = nccreate(ncid,'vp',dimnames3,dimsizes3,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name','prsflux ydir')
        ierr = nf90_put_att(ncid,varid,'units','W/m^2')
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
      integer,   dimension(4) :: start    ! start vector for writing

        if (.not.west_exchng) ud(1,:,:,:) = 0
        start = (/1, bfy, 1, record/)

        call ncwrite(ncid,'u_pgr',ud(:,:,:,1),start)
        call ncwrite(ncid,'u_cor',ud(:,:,:,2),start)
        call ncwrite(ncid,'u_adv',ud(:,:,:,3),start)
        call ncwrite(ncid,'u_dis',ud(:,:,:,4),start)
        call ncwrite(ncid,'u_hmx',ud(:,:,:,5),start)
        call ncwrite(ncid,'u_vmx',ud(:,:,:,6),start)
        call ncwrite(ncid,'u_cpl',ud(:,:,:,7),start)

        if (.not.south_exchng) vd(:,1,:,:) = 0
        start = (/bfx, 1, 1, record/)

        call ncwrite(ncid,'v_pgr',vd(:,:,:,1),start)
        call ncwrite(ncid,'v_cor',vd(:,:,:,2),start)
        call ncwrite(ncid,'v_adv',vd(:,:,:,3),start)
        call ncwrite(ncid,'v_dis',vd(:,:,:,4),start)
        call ncwrite(ncid,'v_hmx',vd(:,:,:,5),start)
        call ncwrite(ncid,'v_vmx',vd(:,:,:,6),start)
        call ncwrite(ncid,'v_cpl',vd(:,:,:,7),start)

      end subroutine wrt_diag_uv  !]
! ----------------------------------------------------------------------
      subroutine wrt_diag_trc(ncid,record,trcd)  ![
      ! Write uv diagnostic vars to file

      use param
      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      real, dimension(0:nx+1,0:ny+1,nz,ntd,ntdia),intent(in) :: trcd

      !local
      integer :: it
      integer,   dimension(4) :: start    ! start vector for writing

      do it = 1,nt
        if (wrt_t_dia(it)) then
          start = (/1, 1, 1, record/)

          call ncwrite(ncid,trim(t_vname(it))//'_advx' ,trcd(i0:i1,j0:j1,:,1,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_advy' ,trcd(i0:i1,j0:j1,:,2,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_advz' ,trcd(i0:i1,j0:j1,:,3,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_mixx' ,trcd(i0:i1,j0:j1,:,4,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_mixy' ,trcd(i0:i1,j0:j1,:,5,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_mixz' ,trcd(i0:i1,j0:j1,:,6,td(it)),start)
        endif
      enddo

      end subroutine wrt_diag_trc  !]
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

!     for 3D up data
!     start = (/bfx, bfy, 1, record/)
!     call ncwrite(ncid,'up',up(bfx:nx,1:ny,:),start)
      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'up',up(bfx:nx,1:ny),start)

!     for 3D vp data
!     start = (/bfx, bfy, 1, record/)
!     call ncwrite(ncid,'vp',vp(1:nx,bfy:ny,:),start)
      start = (/bfx, bfy, record/)
      call ncwrite(ncid,'vp',vp(1:nx,bfy:ny),start)


      end subroutine wrt_diag_pflx  !]
! ----------------------------------------------------------------------
      subroutine wrt_rst_diag_slow(ncid,record)  ![
      ! write the slow components of u, v and P to restart file

      implicit none

      ! import/exprot
      integer,intent(in) :: ncid
      integer,intent(in) :: record

      call ncwrite(ncid,'u_slow',u_slow( 1:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'v_slow',v_slow(i0:i1, 1:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'p_slow',p_slow(i0:i1,j0:j1,:),(/1,1,1,record/))


      end subroutine wrt_rst_diag_slow  !]
! ----------------------------------------------------------------------
      subroutine create_diagnostics_file(fname) ![
      implicit none

      !input/output
      character(len=*),intent(out) :: fname

      !local
      integer :: ncid,ierr

      call create_file('_dia',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      ! Diagnostics options
      ierr=nf90_put_att (ncid, nf90_global, 'diagnostic_options', adjustl(diagnostic_opt))
      if (diag_avg) then
        ierr=nf90_put_att (ncid, nf90_global, 'diag_avg_time', adjustl(diag_avg_output))
        ierr=nf90_put_att(ncid,nf90_global,'filtering_time_scale',timescale )
      endif

      call create_diagvars(ncid)

      ierr = nf90_close(ncid)

      end subroutine create_diagnostics_file !]
! ----------------------------------------------------------------------
      subroutine get_init_slow(ncid,record,tindx )  ![
      ! get initial slow u,v and p
      implicit none

      ! input
      integer,intent(in) :: ncid, record,tindx

      ! local
      integer :: i,j
      integer :: ierr, varid
      integer,dimension(4) :: start

      start=1; start(4)=record                                       ! 3D vars


      if (mynode==0) print *,'getting u,v,p slow'
      if (.not. allocated(u_slow)) then              ! do_diagnostics not called yet so need to init u_slow,p_slow
        allocate(u_slow(nx+1,0:ny+1,nz ))
        u_slow = 0
      endif
      if (.not. allocated(v_slow)) then
        allocate( v_slow(0:nx+1,ny+1,nz ))
        v_slow = 0
      endif
      if (.not. allocated(p_slow)) then
        allocate( p_slow( GLOBAL_2D_ARRAY, nz ))
        p_slow = 0
      endif

      ierr=nf90_inq_varid (ncid, 'u_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,'u_slow', u_slow(1:i1,j0:j1,:),start)
#ifdef MASKING
        do j = j0,j1
          do i = 1,i1
            if (umask(i,j)<1) u_slow(i,j,:) = 0
          enddo
        enddo
#endif
      else
!       u_slow = 0.0
        u_slow(1:i1,j0:j1,:) = u(1:i1,j0:j1,:,tindx)
        mpi_master_only write(*,*) ' --- WARNING: '
     &    ,'u_slow not in initial file.'
     &    ,'  Initialized to u'
        ierr=nf90_noerr
      endif

      ierr=nf90_inq_varid (ncid, 'v_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,'v_slow', v_slow(i0:i1,1:j1,:),start)
#ifdef MASKING
        do j = 1,j1
          do i = i0,i1
            if (vmask(i,j)<1) v_slow(i,j,:) = 0
          enddo
        enddo
#endif
      else
!       v_slow = 0.0
        v_slow(i0:i1,1:j1,:) = v(i0:i1,1:j1,:,tindx)
        mpi_master_only write(*,*) ' --- WARNING: '
     &    ,'v_slow not in initial file.'
     &    ,'  Initialized to v'
        ierr=nf90_noerr
      endif

      ierr=nf90_inq_varid (ncid, 'p_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid, 'p_slow', p_slow(i0:i1,j0:j1,:), start)
        call exchange_xxx(p_slow)  ! need exchange since r2u to get u'p'
#ifdef MASKING
        do j = 0,ny
          do i = 0,nx
            if (rmask(i,j)<1) p_slow(i,j,:) = 0
          enddo
        enddo
#endif
      else
        p_slow = 0.0
        mpi_master_only write(*,*) ' --- WARNING: '
     &         , 'p_slow'
     &         , ' not in initial file.  Initialized to 0.0'
        ierr=nf90_noerr
      endif

      end subroutine get_init_slow  !]
! ----------------------------------------------------------------------

      end module diagnostics
