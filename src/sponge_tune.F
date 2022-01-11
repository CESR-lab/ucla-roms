      module sponge_tune

      ! Tunes the sponge values near the boundaries to attempt
      ! to match the incoming baroclinic wave flux from the parent

#include "cppdefs.opt"
#ifdef SPONGE_TUNE

      use scalars
      use dimensions
      use diagnostics
      use roms_read_write

      implicit none
      private

      logical           :: tune_init = .true.
      integer,parameter :: isp = 15  !! width of the sponge region in grid points
      real              :: sp_timscale = 24*3600
      real              :: pflux_norm = 1e5 ! typical pflux scale

      real,allocatable,dimension(:) :: cflx_west,cflx_east,cflx_north,cflx_south
      real,allocatable,dimension(:) :: pflx_west,pflx_east,pflx_north,pflx_south
      real,allocatable,dimension(:) :: sp_west,sp_east,sp_north,sp_south
      real,allocatable,dimension(:) :: sp_mult
      real,allocatable,dimension(:,:):: wrk

      logical           :: obc_west =.false.
      logical           :: obc_east =.false.
      logical           :: obc_north=.false.
      logical           :: obc_south=.false.

      public adjust_sponge

#include "sponge_tune.opt"
      contains

! ----------------------------------------------------------------------
      subroutine init_sponge_tune
      ! Initialize sponge tuning
      implicit none

      ! local
      integer :: ierr,i

      allocate(wrk(nx,ny))
      allocate(sp_mult(isp))
      do i=1,isp
       sp_mult(i) = float(i-0.5)/float(isp)
      enddo
# ifdef OBC_WEST
      obc_west =.true.
# endif
# ifdef OBC_EAST
      obc_east =.true.
# endif
# ifdef OBC_NORTH
      obc_north=.true.
# endif
# ifdef OBC_SOUTH
      obc_south=.true.
# endif


      if (obc_west.and.(.not.west_exchng)) then
        if (nx<isp) ierr = 1
        allocate(pflx_west(ny))
        allocate(cflx_west(ny))
        allocate(sp_west(ny))
        sp_west = v_sponge
        cflx_west = 0
      endif
      if (obc_east.and.(.not.east_exchng)) then
        if (nx<isp) ierr = 1
        allocate(pflx_east(ny))
        allocate(cflx_east(ny))
        allocate(sp_east(ny))
        sp_east = v_sponge
        cflx_east = 0
      endif
      if (obc_north.and.(.not.north_exchng)) then
        if (ny<isp) ierr = 1
        allocate(pflx_north(nx))
        allocate(cflx_north(nx))
        allocate(sp_north(nx))
        sp_north = v_sponge
        cflx_north = 0
      endif
      if (obc_south.and.(.not.south_exchng)) then
        if (ny<isp) ierr = 1
        allocate(pflx_south(nx))
        allocate(cflx_south(nx))
        allocate(sp_south(nx))
        sp_south = v_sponge
        cflx_south = 0
      endif

      tune_init = .false.

      ! Check that the subdomains are larger than the sponge_size
      if (ierr>0) then
        print *,'Subdomain size is less than sponge_size'
        stop 'fatal sponge tune error'
      endif

      end subroutine init_sponge_tune
  
! ----------------------------------------------------------------------
      subroutine set_pflx
      ! Read parent grid baroclinic pressure fluxes from bry file
      ! and interpolate to the correct time

      if (obc_west.and.(.not.west_exchng)) then
        set_bry_data(nc_pflx_w,pflx_west)
      endif
!     ncread

      end subroutine read_pflx
      
! ----------------------------------------------------------------------
      subroutine comp_cflx
      ! Compute integrated, averaged and filtered baroclinic pressure fluxes

      ! wrk should point to some allocated space of size(nx,ny)

      if (obc_west.and.(.not.west_exchng)) then
        wrk(1:isp,:) = sum(up_pp(1:isp,1:ny,:),dim=3)
        cflx_west = cflx_west*beta + 
     &             alpha*sum(wrk(1:isp,:),dim=1)/isp
      endif
      if (obc_east.and.(.not.east_exchng)) then
        wrk(1:isp,:) = sum(up_pp(nx-isp+1:nx,1:ny,:),dim=3)
        cflx_east = cflx_east*beta + 
     &             alpha*sum(wrk(1:isp,:),dim=1)/isp
      endif
      
      end subroutine comp_cflx
      
! ----------------------------------------------------------------------
      subroutine adjust_sponge
      ! Adjust the value of the sponge up or down in order
      ! to match the Incoming parent baroclinic pressure flux
      ! Parent and child fluxes are time filtered

      use mixing
      implicit none

      ! local 
      real :: i,j,alpha

      alpha = dt/sp_timscale/pflux_norm

      if (tune_init) call init_sponge_tune

      call set_pflx
      call comp_cflx

      if (obc_west.and.(.not.west_exchng)) then
        sp_west = sp_west - alpha*(pflx_west - cflx_west)
        sp_west = max(0.,sp_west)
        do j = 1,ny
          visc2_r(1:isp,j) = sp_mult*sp_west(j)
        enddo
      endif

      if (obc_east.and.(.not.east_exchng)) then
        sp_east = sp_east - alpha*(pflx_east - cflx_east)
        sp_east = max(0.,sp_east)
        do j = 1,ny
          visc2_r(nx:-1:nx-isp+1,j) = sp_mult*sp_west(j)
        enddo
      endif

      if (obc_south.and.(.not.south_exchng)) then
        sp_south = sp_south + alpha*(pflx_south - cflx_south)
        sp_south = max(0.,sp_south)
        do i = 1,nx
          visc2_r(i,1:ny) = sp_mult*sp_south(i)
        enddo
      endif

      ! north boundary
      if (obc_north.and.(.not.north_exchng)) then
        sp_north = sp_north+ alpha*(pflx_north- cflx_north)
        sp_north = max(0.,sp_north)
        do i = 1,nx
          visc2_r(i,ny:-1:ny-isp+1) = sp_mult*sp_south(i)
        enddo
      endif

# ifdef EXCHANGE
      call exchange2d_tile(1,nx,1,ny,visc2_r)
# endif

      visc2_p(1:nx+1,1:ny+1) = 0.25*(
     &           visc2_r(1:nx+1,1:ny+1) + visc2_r(0:nx,1:ny+1)    
     &         + visc2_r(1:nx+1,0:ny  ) + visc2_r(0:nx,0:ny  ) )


      end subroutine adjust_sponge

#endif /* SPONGE_TUNE */

      end module sponge_tune
