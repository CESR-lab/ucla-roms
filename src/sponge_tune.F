      module sponge_tune

      ! Tunes the sponge values near the boundaries to attempt
      ! to match the incoming baroclinic wave flux from the parent

#include "cppdefs.opt"
#ifdef SPONGE_TUNE

      use netcdf
      use nc_read_write
      use scalars
      use dimensions
      use diagnostics
      use roms_read_write
      use hidden_mpi_vars
      use mixing
      use read_write, only: output_root_name,handle_ierr

      implicit none
      private

#include "sponge_tune.opt"

      integer,parameter :: isp   = 15  !! width of the sponge region in grid points
      logical,parameter :: onoff = .true.
      real,   parameter :: sp_timscale = 24*3600
      real,   parameter :: pflux_norm = 1e-2 ! typical pflux scale
      integer,parameter :: nspt = 10   ! number of records per file
      real   ,parameter :: output_period = 20
      logical           :: tune_init = .true.

      integer           :: total_rec=0       ! in the future, read in from restart file
      real              :: last_time = -99999! in the future, read in from restart file

      real              :: sp_mx             ! max value of sponge


      real,allocatable,dimension(:) :: pflx_west,pflx_east,pflx_north,pflx_south
      real,allocatable,dimension(:) :: cflx_west,cflx_east,cflx_north,cflx_south
      real,allocatable,dimension(:) :: sp_west,sp_east,sp_north,sp_south
      real,allocatable,dimension(:) :: sp_mult
      real,allocatable,dimension(:,:):: wrk

      logical           :: obc_west =.false.
      logical           :: obc_east =.false.
      logical           :: obc_north=.false.
      logical           :: obc_south=.false.

      public adjust_sponge

      contains

! ----------------------------------------------------------------------
      subroutine init_sponge_tune
      ! Initialize sponge tuning
      implicit none

      ! local
      integer :: ierr,i

      if (.not.diag_pflx) then
        if(mynode==0) print *,'For sponge tuning,'
        if(mynode==0) print *,'pflx diagnostics must be turned on!'
        stop
      endif
      allocate(wrk(nx,ny))
      allocate(sp_mult(isp))
      do i=1,isp
       sp_mult(i) = float(isp-i+0.5)/float(isp)
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

      !! max value based of cfl condition
      sp_mx = 0.1/(minval(pm)*minval(pn)*dt)

      !! reset viscosity
      visc2_r = 0

      if (obc_west.and.(.not.west_exchng)) then
        if (nx<isp) ierr = 1
        if (.not.onoff) then
          allocate(nc_pflx_w%vdata(ny,1,2))
          allocate(pflx_west(ny))
        endif
        allocate(cflx_west(ny))
        allocate(sp_west(ny))
        sp_west = v_sponge
        cflx_west = 0
      endif
      if (obc_east.and.(.not.east_exchng)) then
        if (nx<isp) ierr = 1
        if (.not.onoff) then
          allocate(nc_pflx_e%vdata(ny,1,2))
          allocate(pflx_east(ny))
        endif
        allocate(cflx_east(ny))
        allocate(sp_east(ny))
        sp_east = v_sponge
        cflx_east = 0
      endif
      if (obc_north.and.(.not.north_exchng)) then
        if (ny<isp) ierr = 1
        if (.not.onoff) then
          allocate(nc_pflx_n%vdata(nx,1,2))
          allocate(pflx_north(nx))
        endif
        allocate(cflx_north(nx))
        allocate(sp_north(nx))
        sp_north = v_sponge
        cflx_north = 0
      endif
      if (obc_south.and.(.not.south_exchng)) then
        if (ny<isp) ierr = 1
        if (.not.onoff) then
          allocate(nc_pflx_s%vdata(nx,1,2))
          allocate(pflx_south(nx))
        endif
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
        call set_frc_data(nc_pflx_w,pflx_west)
      endif
      if (obc_east.and.(.not.east_exchng)) then
        call set_frc_data(nc_pflx_e,pflx_east)
      endif
      if (obc_south.and.(.not.south_exchng)) then
        call set_frc_data(nc_pflx_s,pflx_south)
      endif
      if (obc_north.and.(.not.north_exchng)) then
        call set_frc_data(nc_pflx_n,pflx_north)
      endif

      end subroutine set_pflx
      
! ----------------------------------------------------------------------
      subroutine comp_cflx
      ! Compute integrated, averaged and filtered baroclinic pressure fluxes
      ! The fluxes are all into the domain!!

      if (obc_west.and.(.not.west_exchng)) then
!       wrk(1:isp,:) = sum(up_pp(1:isp,1:ny,:),dim=3)
        wrk(1:isp,:) = up(1:isp,1:ny)  !! up is now already integrated
        cflx_west = cflx_west*beta + 
     &             alpha*sum(wrk(1:isp,:),dim=1)/isp
      endif
      if (obc_east.and.(.not.east_exchng)) then
!       wrk(1:isp,:) = sum(up_pp(nx-isp+1:nx,1:ny,:),dim=3)
        wrk(1:isp,:) = up(nx-isp+1:nx,1:ny)
        cflx_east = cflx_east*beta -
     &             alpha*sum(wrk(1:isp,:),dim=1)/isp
      endif

      if (obc_south.and.(.not.south_exchng)) then
!       wrk(:,1:isp) = sum(vp_pp(1:nx,1:isp,:),dim=3)
        wrk(:,1:isp) = vp(1:nx,1:isp)
        cflx_south = cflx_south*beta + 
     &             alpha*sum(wrk(:,1:isp),dim=2)/isp
!       print *,'comp_cflx: ', sum(cflx_south)/nx
      endif
      if (obc_north.and.(.not.north_exchng)) then
!       wrk(:,1:isp) = sum(vp_pp(1:nx,ny+1-isp:ny,:),dim=3)
        wrk(:,1:isp) = vp(1:nx,ny+1-isp:ny)
        cflx_south = cflx_south*beta - 
     &             alpha*sum(wrk(:,1:isp),dim=2)/isp
      endif
      
      end subroutine comp_cflx
      
! ----------------------------------------------------------------------
      subroutine adjust_sponge ![
      ! Adjust the value of the boundary sponge 
      ! based on the incoming pressure flux.
      ! Onoff: Based on the sign of the child flux
      ! Else:  Based on the difference between parent and child flux
      implicit none

      ! local 
      real    :: alpha,beta,eps
      integer :: i,j,idx,itrc

      eps = 1e-6
      if (onoff) then
        alpha = dt/sp_timscale
        beta  = 1-alpha
      else
        alpha = dt/sp_timscale/pflux_norm
      endif

      if (tune_init) call init_sponge_tune

      if (.not.onoff) call set_pflx

      call comp_cflx

      if (obc_south.and.(.not.south_exchng)) then
!       print *, 'adjusting: ',sum(sp_south)/nx,mynode
!       print *, 'adjusting: ',sum(cflx_south)/nx,mynode
!       print *,sum(0.5*(1-cflx_south/abs(cflx_south+eps)))/nx
        if (onoff) then
          sp_south = sp_south*beta + 
     &       alpha*v_sponge*0.5*(1-cflx_south/abs(cflx_south+eps))
        else
          sp_south = sp_south + alpha*(pflx_south - cflx_south)
          sp_south = max(0.,sp_south)
          sp_south = min(sp_mx,sp_south)
        endif
!       print *, 'adjusting: ',sum(sp_south)/nx,mynode
        do i = 1,nx
          visc2_r(i,1:isp) = sp_mult*sp_south(i)
        enddo
      endif

      if (obc_north.and.(.not.north_exchng)) then
        if (onoff) then
          sp_north = sp_north*beta + 
     &       alpha*v_sponge*0.5*(1-cflx_north/abs(cflx_north+eps))
        else
          sp_north = sp_north + alpha*(pflx_north- cflx_north)
          sp_north = max(0.,sp_north)
          sp_north = min(sp_mx,sp_north)
        endif
        do i = 1,nx
          visc2_r(i,ny:ny-isp+1:-1) = sp_mult*sp_north(i)
        enddo
      endif
      
      idx = isp
      if (obc_west.and.(.not.west_exchng)) then
        if (onoff) then
          sp_west = sp_west*beta +
     &       alpha*v_sponge*0.5*(1-cflx_west/abs(cflx_west+eps))
        else
          sp_west = sp_west - alpha*(pflx_west - cflx_west)
          sp_west = max(0.,sp_west)
          sp_west = min(sp_mx,sp_west)
        endif
        do j = 1,ny
          if (.not.south_exchng) idx = min(j,isp)
          if (.not.north_exchng) idx = min(ny-j+1,isp)
          visc2_r(1:idx,j) = sp_mult(1:idx)*sp_west(j)
        enddo
      endif

      idx = isp
      if (obc_east.and.(.not.east_exchng)) then
        if (onoff) then
          sp_east = sp_east*beta +
     &       alpha*v_sponge*0.5*(1-cflx_east/abs(cflx_east+eps))
        else
          sp_east = sp_east - alpha*(pflx_east - cflx_east)
          sp_east = max(0.,sp_east)
          sp_east = min(sp_mx,sp_east)
        endif
        do j = 1,ny
          if (.not.south_exchng) idx = min(j,isp)
          if (.not.north_exchng) idx = min(ny-j+1,isp)
          visc2_r(nx:nx-idx+1:-1,j) = sp_mult(1:idx)*sp_east(j)
        enddo
      endif

# ifdef EXCHANGE
      call exchange2d_tile(1,nx,1,ny,visc2_r)
# endif

      do itrc = 1,nt
        diff2(:,:,itrc) = visc2_r
      enddo

      visc2_p(1:nx+1,1:ny+1) = 0.25*(
     &           visc2_r(1:nx+1,1:ny+1) + visc2_r(0:nx,1:ny+1)    
     &         + visc2_r(1:nx+1,0:ny  ) + visc2_r(0:nx,0:ny  ) )

      if (time>=last_time+output_period) then
        call write_sp_tune
        last_time = time
      endif

      end subroutine adjust_sponge !]

! ----------------------------------------------------------------------
      subroutine write_sp_tune ![
      implicit none

      !local
      integer            :: ncid,ierr
      integer            :: record
      integer            :: i0,i1,j0,j1
      integer,dimension(3) :: start    
      character(len=99)  :: fname
      save fname,record
      
      if (mod(total_rec,nspt)==0) then
        call create_sp_tune_file(total_rec,fname)
        record = 0
      endif
      total_rec  = total_rec +1
      record = record + 1

      ierr=nf90_open(fname, nf90_write, ncid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'write sponge_tune: ',fname)

      !! JM Simplify when submains can be unequal in size
      !! JM Already, we could move this to mpi_setup
      i0 = iwest
      i1 = ieast
      j0 = jsouth
      j1 = jnorth
      if (inode==0) i0=i0-1
      if (inode==NP_XI-1) i1=i1+1
      if (jnode==0) j0=j0-1
      if (jnode==NP_ETA-1) j1=j1+1
      
      start = (/1, 1, record/)
      call ncwrite(ncid,'visc2d',visc2_r(i0:i1,j0:j1),start)
      
      ierr=nf90_close(ncid)    

      end subroutine write_sp_tune !]

! ----------------------------------------------------------------------
      subroutine create_sp_tune_file(total_rec,fname) ![
      implicit none

      !input/output
      integer          ,intent(in)  :: total_rec  ! total netcdf records so far for these variables
      character(len=99),intent(out) :: fname

      !local
      integer :: ncid,ierr,varid
      integer :: prec ! precision of output
      character(len=10),dimension(3) :: dimnames ! dimension names
      integer,          dimension(3) :: dimsizes 

      fname=trim(output_root_name) / / '_spn.nc'
      call create_filename(fname,total_rec)
      if (mynode==0) print *,'sponge tune :: created new netCDF file: ',fname

      ierr=nf90_create(fname,nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'create_sp_tune_file: ')

      call put_global_atts(ncid, ierr)  ! put global attributes in file

      dimnames = (/'xi_rho','eta_rho','time'/)
      dimsizes = (/xi_rho,eta_rho,0/)

      prec = nf90_float
      varid = nccreate(ncid,'visc2d',dimnames,dimsizes,prec)
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'horizontal mixing rate')
      ierr = nf90_put_att(ncid,varid,'units','m^2/s' )

      ierr = nf90_enddef(ncid)

      end subroutine create_sp_tune_file !]



#endif /* SPONGE_TUNE */

      end module sponge_tune
