      module sponge_tune

      ! Tunes the sponge values near the boundaries to attempt
      ! to match the incoming baroclinic wave flux from the parent

#include "cppdefs.opt"
#ifdef SPONGE_TUNE

      use netcdf
      use nc_read_write
      use scalars
      use dimensions
      use ocean_vars
      use diagnostics
      use roms_read_write
      use hidden_mpi_vars
      use mixing

      implicit none
      private

#include "sponge_tune.opt"

      logical           :: tune_init = .true.

      real              :: last_time = -99999! in the future, read in from restart file
      integer           :: record=nrpf

      real              :: sp_mx             ! max value of sponge

      real,allocatable,dimension(:) :: pflx_west,pflx_east,pflx_north,pflx_south
      real,allocatable,dimension(:) :: cflx_west,cflx_east,cflx_north,cflx_south
      real,allocatable,dimension(:) :: sp_west,sp_east,sp_north,sp_south
      real,allocatable,dimension(:) :: sp_mult
      real,allocatable,dimension(:,:):: wrk


      public adjust_sponge

      contains

! ----------------------------------------------------------------------
      subroutine init_sponge_tune ![
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

      !! max value based of cfl condition
!     sp_mx = 0.06/(maxval(pm)*maxval(pn)*dt)
      sp_mx = v_sponge

      !! reset viscosity
      visc2_r = 0

      ierr = 0
      if (obc_west.and.(.not.west_exchng)) then
        if (nx<isp) ierr = 1
        allocate(nc_pflx_w%vdata(ny,1,2))
        allocate(pflx_west(ny))
        pflx_west = 0
        allocate(cflx_west(ny))
        allocate(sp_west(ny))
        sp_west = v_sponge
        cflx_west = 0
      endif
      if (obc_east.and.(.not.east_exchng)) then
        if (nx<isp) ierr = 1
        allocate(nc_pflx_e%vdata(ny,1,2))
        allocate(pflx_east(ny))
        pflx_east = 0
        allocate(cflx_east(ny))
        allocate(sp_east(ny))
        sp_east = v_sponge
        cflx_east = 0
      endif
      if (obc_north.and.(.not.north_exchng)) then
        if (ny<isp) ierr = 1
        allocate(nc_pflx_n%vdata(nx,1,2))
        allocate(pflx_north(nx))
        pflx_north = 0
        allocate(cflx_north(nx))
        allocate(sp_north(nx))
        sp_north = v_sponge
        cflx_north = 0
      endif
      if (obc_south.and.(.not.south_exchng)) then
        if (ny<isp) ierr = 1
        allocate(nc_pflx_s%vdata(nx,1,2))
        allocate(pflx_south(nx))
        pflx_south = 0
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

      end subroutine init_sponge_tune !]
! ----------------------------------------------------------------------
      subroutine set_pflx ![
      ! Read parent grid baroclinic pressure fluxes from bry file
      ! and interpolate to the correct time

      use grid

      !local
      integer:: i,j,ierr

      if (obc_west.and.(.not.west_exchng)) then
        call set_frc_data(nc_pflx_w,pflx_west)
        do j=1,ny
          if (abs(pflx_west(j))>100) pflx_west(j) = 0.
        enddo
      endif
      if (obc_east.and.(.not.east_exchng)) then
        call set_frc_data(nc_pflx_e,pflx_east)
        do j=1,ny
          if (abs(pflx_east(j))>100) pflx_east(j) = 0.
        enddo
      endif
      if (obc_south.and.(.not.south_exchng)) then
        call set_frc_data(nc_pflx_s,pflx_south)
        do i=1,nx
          if (abs(pflx_south(i))>100) pflx_south(i) = 0.
        enddo
      endif
      if (obc_north.and.(.not.north_exchng)) then
        call set_frc_data(nc_pflx_n,pflx_north)
        do i=1,nx
          if (rmask(i,ny+1)<1) pflx_north(i) = 0.
          if (abs(pflx_north(i))>100) pflx_north(i) = 0.
        enddo
      endif

      end subroutine set_pflx !]
! ----------------------------------------------------------------------
      subroutine comp_cflx ![
      ! Compute integrated, averaged and filtered baroclinic pressure fluxes
      ! up/vp fluxes are positive to the east/north
      ! Instead of taking the mean flux over the sponge region, we're 
      ! taking the flux on the 'inside' boundary of the sponge region.

      if (obc_west.and.(.not.west_exchng)) then
        cflx_west = up(isp,1:ny) 
      endif
      if (obc_east.and.(.not.east_exchng)) then
        cflx_east = up(nx-isp+1,1:ny)
      endif
      if (obc_south.and.(.not.south_exchng)) then
        cflx_south = vp(1:nx,isp)
      endif
      if (obc_north.and.(.not.north_exchng)) then
        cflx_north = vp(1:nx,ny-isp+1)
      endif
      
      end subroutine comp_cflx !]
! ----------------------------------------------------------------------
      subroutine adjust_sponge ![
      ! Adjust the value of the boundary sponge 
      ! based on the incoming pressure flux.
      ! Onoff: Based on the sign of the child flux
      ! Else:  Based on the difference between parent and child flux
      implicit none

      ! local 
      real    :: alpha,beta,eps,sp_mn
      integer :: i,j,k,idx,itrc
      real    :: damp,u_corr,v_corr

      eps = 1e-6
      sp_mn = 0.0*sp_mx

      alpha = sp_mx*2*dt/sp_timscale

      if (tune_init) call init_sponge_tune

      call set_pflx
      call comp_cflx

      if (obc_south.and.(.not.south_exchng)) then
!       sp_south = sp_south*(1-alpha* 
!    &  (pflx_south-cflx_south)/(abs(pflx_south)+abs(cflx_south)+eps))
        sp_south = sp_south-alpha* 
     &  (pflx_south-cflx_south)/(abs(pflx_south)+abs(cflx_south)+eps)
        sp_south = max(sp_mn,sp_south)
        sp_south = min(sp_mx,sp_south)
        if (mynode==10) write(10,*),pflx_south(nx/2),cflx_south(nx/2),sp_south(nx/2)
        do i = 1,nx
!         visc2_r(i,1:isp) = sp_mult*sp_south(i)
!         visc2_r(i,0) = visc2_r(i,1)
          visc2_r(i,1:isp) = sp_mult*0.2*sp_mx
          visc2_r(i,0) = visc2_r(i,1)
        enddo
        print *,'max/min: ',maxval(sp_south)/sp_mx
        do k = 1,nz
         do j = 1,isp
          do i = 1,nx
            damp = sp_south(i)/sp_mx
            u_corr = -damp*sp_mult(j)*(u(i,j,k,nnew)-u_slow(i,j,k))
            u(i,j,k,nnew) = u(i,j,k,nnew) + u_corr
            v_corr = -damp*sp_mult(j)*(v(i,j,k,nnew)-v_slow(i,j,k))
            v(i,j,k,nnew) = u(i,j,k,nnew) + v_corr
          enddo
         enddo
        enddo
      endif

      if (obc_north.and.(.not.north_exchng)) then
!       sp_north = sp_north*(1+alpha*
!    &   (pflx_north-cflx_north)/(abs(pflx_north)+abs(cflx_north)+eps))
        sp_north = sp_north+alpha*
     &   (pflx_north-cflx_north)/(abs(pflx_north)+abs(cflx_north)+eps)
        sp_north = max(sp_mn,sp_north)
        sp_north = min(sp_mx,sp_north)
        do i = 1,nx
          visc2_r(i,ny:ny-isp+1:-1) = sp_mult*sp_north(i)
          visc2_r(i,ny+1) = visc2_r(i,ny)
        enddo
      endif
      
      idx = isp
      if (obc_west.and.(.not.west_exchng)) then
        sp_west = sp_west-alpha*
     &  (pflx_west-cflx_west)/(abs(pflx_west)+abs(cflx_west)+eps)
        sp_west = max(sp_mn,sp_west)
        sp_west = min(sp_mx,sp_west)
        do j = 1,ny ! to get a 'diagonal' in the corners
          if (.not.south_exchng) idx = min(j,isp)
          if (.not.north_exchng) idx = min(ny-j+1,isp)
          visc2_r(1:idx,j) = sp_mult(1:idx)*sp_west(j)
          visc2_r(0,j) = visc2_r(1,j)
        enddo
      endif

      idx = isp
      if (obc_east.and.(.not.east_exchng)) then
        sp_east = sp_east+alpha* 
     &  (pflx_east-cflx_east)/(abs(pflx_east)+abs(cflx_east)+eps)
        sp_east = max(sp_mn,sp_east)
        sp_east = min(sp_mx,sp_east)
        if (mynode==10) write(11,*) pflx_east(ny/2),cflx_east(ny/2),sp_east(ny/2)
        do j = 1,ny
          if (.not.south_exchng) idx = min(j,isp)
          if (.not.north_exchng) idx = min(ny-j+1,isp)
          visc2_r(nx:nx-idx+1:-1,j) = sp_mult(1:idx)*sp_east(j)
          visc2_r(nx+1,j) = visc2_r(nx,j)
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

      if (time>=last_time+output_period .and. wrt_sponge) then
        call write_sp_tune
        last_time = time
      endif

      end subroutine adjust_sponge !]
! ----------------------------------------------------------------------
      subroutine write_sp_tune ![
      implicit none

      !local
      integer            :: ncid,ierr
      integer            :: i0,i1,j0,j1
      integer,dimension(3) :: start    
      character(len=99)  :: fname
      save fname
      
      if (record==nrpf) then
        call create_sp_tune_file(fname)
        record = 0
      endif
      record = record + 1

      ierr=nf90_open(fname, nf90_write, ncid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'write sponge_tune: ',fname)

      call ncwrite(ncid,'ocean_time',(/time/),(/record/))

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
      subroutine create_sp_tune_file(fname) ![
      implicit none

      !input/output
      character(len=*),intent(out) :: fname

      !local
      integer :: ncid,ierr,varid
      character(len=10),dimension(3) :: dimnames ! dimension names
      integer,          dimension(3) :: dimsizes 

      call create_file('_spn',fname)

      dimnames = (/'xi_rho','eta_rho','time'/)
      dimsizes = (/xi_rho,eta_rho,0/)

      ierr=nf90_open(fname,nf90_write,ncid)

      varid = nccreate(ncid,'visc2d',dimnames,dimsizes,nf90_float)
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'horizontal mixing rate')
      ierr = nf90_put_att(ncid,varid,'units','m^2/s' )

      ierr = nf90_close(ncid)

      end subroutine create_sp_tune_file !]



#endif /* SPONGE_TUNE */

      end module sponge_tune
