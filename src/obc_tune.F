      module obc_tune

      ! Tunes the boundary conditions at open boundaries to attempt
      ! to match the incoming baroclinic wave flux from the parent

#include "cppdefs.opt"

      use netcdf
      use nc_read_write
      use scalars
      use dimensions
      use diagnostics
      use roms_read_write
      use roms_mpi
      use mixing

      implicit none
      private

#include "obc_tune.opt"

      logical   :: tune_init = .true.

      real      :: output_time = 0
      real      :: navg = 0
      integer   :: record=nrpf     ! triggers creation of initial file


      real,allocatable,dimension(:) :: pflx_west,pflx_east,pflx_north,pflx_south
      real,allocatable,dimension(:) :: cflx_west,cflx_east,cflx_north,cflx_south
      real,allocatable,public,dimension(:) :: ub_west,ub_east,ub_north,ub_south

      real,allocatable,dimension(:) :: cflx_south_avg
      real,allocatable,dimension(:) :: pflx_south_avg
      real,allocatable,dimension(:) ::   ub_south_avg

      real,allocatable,dimension(:) :: cflx_north_avg
      real,allocatable,dimension(:) :: pflx_north_avg
      real,allocatable,dimension(:) ::   ub_north_avg

      real,allocatable,dimension(:) :: cflx_west_avg
      real,allocatable,dimension(:) :: pflx_west_avg
      real,allocatable,dimension(:) ::   ub_west_avg

      real,allocatable,dimension(:) :: cflx_east_avg
      real,allocatable,dimension(:) :: pflx_east_avg
      real,allocatable,dimension(:) ::   ub_east_avg

      public adjust_orlanski
      public init_orlanski_tune
      public wrt_rst_tune
      public get_init_tune

      contains

! ----------------------------------------------------------------------
      subroutine init_orlanski_tune ![
      ! Initialize obc tuning
      implicit none

      ! local
      integer :: i

      if (.not.diag_pflx) then
        if (mynode==0) then
          print *,'For Orlanski tuning,'
          print *,'pflx diagnostics must be turned on!'
          print *,'(diag_pflx=.true. in diagnostics.opt)'
        endif
        stop
      endif

      call store_string_att(bc_options,'obc_tuning, ')


      !! max value based of cfl condition
!     sp_mx = 0.06/(maxval(pm)*maxval(pn)*dt)

      if (obc_west.and.inode.eq.0) then
        allocate(nc_pflx_w%vdata(j0:j1,1,2))
        allocate(pflx_west(0:ny+1))
        pflx_west = 0
        allocate(pflx_west_avg(0:ny+1))
        pflx_west_avg = 0
        allocate(cflx_west(0:ny+1))
        cflx_west = 0
        allocate(cflx_west_avg(0:ny+1))
        cflx_west_avg = 0
        if (.not.allocated(ub_west)) then
          allocate(ub_west(0:ny+1))
          ub_west = 1.
        endif
        allocate(ub_west_avg(0:ny+1))
        ub_west_avg = 0.
      endif
      if (obc_east.and.inode.eq.npx-1) then
        allocate(nc_pflx_e%vdata(j0:j1,1,2))
        allocate(pflx_east(0:ny+1))
        pflx_east = 0
        allocate(pflx_east_avg(0:ny+1))
        pflx_east_avg = 0
        allocate(cflx_east(0:ny+1))
        cflx_east = 0
        allocate(cflx_east_avg(0:ny+1))
        cflx_east_avg = 0
        if (.not.allocated(ub_east)) then
          allocate(ub_east(0:ny+1))
          ub_east = 1.
        endif
        allocate(ub_east_avg(0:ny+1))
        ub_east_avg = 0.
      endif
      if (obc_north.and.jnode.eq.npy-1) then
        allocate(nc_pflx_n%vdata(i0:i1,1,2))
        allocate(pflx_north(0:nx+1))
        pflx_north = 0
        allocate(pflx_north_avg(0:nx+1))
        pflx_north_avg = 0
        allocate(cflx_north(0:nx+1))
        cflx_north = 0
        allocate(cflx_north_avg(0:nx+1))
        cflx_north_avg = 0
        if (.not.allocated(ub_north)) then
          allocate(ub_north(0:nx+1))
          ub_north = 1.
        endif
        allocate(ub_north_avg(0:nx+1))
        ub_north_avg = 0
      endif
      if (obc_south.and.jnode.eq.0) then
        allocate(nc_pflx_s%vdata(i0:i1,1,2))
        allocate(pflx_south(0:nx+1))
        pflx_south = 0
        allocate(pflx_south_avg(0:nx+1))
        pflx_south_avg = 0
        allocate(cflx_south(0:nx+1))
        cflx_south = 0
        allocate(cflx_south_avg(0:nx+1))
        cflx_south_avg = 0
        if (.not.allocated(ub_south)) then
          allocate(ub_south(0:nx+1))
          ub_south = 1.
        endif
        allocate(ub_south_avg(0:nx+1))
        ub_south_avg = 0
      endif

      tune_init = .false.

      end subroutine init_orlanski_tune !]
! ----------------------------------------------------------------------
      subroutine set_pflx ![
      ! Read parent grid baroclinic pressure fluxes from bry file
      ! and interpolate to the correct time

      use grid

      !local
      integer:: i,j,ierr

      if (obc_west.and.inode.eq.0) then
        call set_frc_data(nc_pflx_w,pflx_west(j0:j1))
        do j=1,ny
          if (abs(pflx_west(j))>100) pflx_west(j) = 0.
        enddo
      endif
      if (obc_east.and.inode.eq.npx-1) then
        call set_frc_data(nc_pflx_e,pflx_east(j0:j1))
        do j=1,ny
          if (abs(pflx_east(j))>100) pflx_east(j) = 0.
        enddo
      endif
      if (obc_south.and.jnode.eq.0) then
        call set_frc_data(nc_pflx_s,pflx_south(i0:i1))
        do i=1,nx
          if (abs(pflx_south(i))>100) pflx_south(i) = 0.
        enddo
      endif
      if (obc_north.and.jnode.eq.npy-1) then
        call set_frc_data(nc_pflx_n,pflx_north(i0:i1))
        do i=1,nx
!         if (rmask(i,ny+1)<1) pflx_north(i) = 0.
          if (abs(pflx_north(i))>100) pflx_north(i) = 0.
        enddo
      endif

      end subroutine set_pflx !]
! ----------------------------------------------------------------------
      subroutine comp_cflx ![
      ! Compute  baroclinic pressure fluxes
      ! up/vp fluxes are positive to the east/north
      ! Instead of taking the mean flux near the boundary, we're 
      ! taking the flux just inside the domain.

      ! local
      integer :: pos

      pos = 2
      if (obc_west.and.inode.eq.0) then
        cflx_west(1:ny) = up(pos,1:ny) 
      endif
      if (obc_east.and.inode.eq.npx-1) then
        cflx_east(1:ny) = up(nx-pos+1,1:ny)
      endif
      if (obc_south.and.jnode.eq.0) then
        cflx_south(1:nx) = vp(1:nx,pos)
      endif
      if (obc_north.and.jnode.eq.npy-1) then
        cflx_north(1:nx) = vp(1:nx,ny-pos+1)
      endif
      
      end subroutine comp_cflx !]
! ----------------------------------------------------------------------
      subroutine adjust_orlanski ![
      ! Adjust the restoring rate at the boundary
      ! by modifying ubind based on the difference 
      ! between parent and child pressure flux
      implicit none

      ! local 
      real    :: alpha,beta,eps
      real    :: ub_mn,ub_mx
      integer :: i,j

      eps = 1e-6
      ub_mx = 2.0 
      ub_mn =-1.0

      alpha = 0.5*ub_mx*dt/tune_timscale

      if (tune_init) call init_orlanski_tune

      call set_pflx
      call comp_cflx

      if (obc_south.and.jnode.eq.0) then
        do i=1,nx
          ub_south(i) = ub_south(i)+alpha*(pflx_south(i)-cflx_south(i))
          ub_south(i) = max(ub_mn,ub_south(i))
          ub_south(i) = min(ub_mx,ub_south(i))
        enddo
      endif

      if (obc_north.and.jnode.eq.npy-1) then
        do i=1,nx
          ub_north(i) = ub_north(i)-alpha*(pflx_north(i)-cflx_north(i))
          ub_north(i) = max(ub_mn,ub_north(i))
          ub_north(i) = min(ub_mx,ub_north(i))
        enddo
      endif
      
      if (obc_west.and.inode.eq.0) then
        do j=1,ny
          ub_west(j) = ub_west(j)+alpha*(pflx_west(j)-cflx_west(j))
          ub_west(j) = max(ub_mn,ub_west(j))
          ub_west(j) = min(ub_mx,ub_west(j))
        enddo
      endif

      if (obc_east.and.inode.eq.npx-1) then
        do j=1,ny
          ub_east(j) = ub_east(j)-alpha*(pflx_east(j)-cflx_east(j))
          ub_east(j) = max(ub_mn,ub_east(j))
          ub_east(j) = min(ub_mx,ub_east(j))
        enddo
      endif

      if (wrt_obc_tune) then
        output_time = output_time + dt
        call calc_tune_avg
        if (output_time>=output_period) then
          call write_obc_tune
          output_time = 0
        endif
      endif

      end subroutine adjust_orlanski !]
! ----------------------------------------------------------------------
      subroutine calc_tune_avg ![
      ! Update diagnostics averages
      ! The average is always scaled properly throughout
      ! reset navg_diag=0 after an output of the average
      use param
      implicit none

      ! local
      real :: coef

      navg = navg +1

      coef = 1./navg

      if (tune_avg) then
        if (obc_south.and.jnode.eq.0) then
          cflx_south_avg = cflx_south_avg*(1-coef) + cflx_south*coef
          pflx_south_avg = pflx_south_avg*(1-coef) + pflx_south*coef
          ub_south_avg   =   ub_south_avg*(1-coef) +   ub_south*coef
        endif

        if (obc_north.and.jnode.eq.npy-1) then
          cflx_north_avg = cflx_north_avg*(1-coef) + cflx_north*coef
          pflx_north_avg = pflx_north_avg*(1-coef) + pflx_north*coef
          ub_north_avg   =   ub_north_avg*(1-coef) +   ub_north*coef
        endif

        if (obc_west.and.inode.eq.0) then
          cflx_west_avg = cflx_west_avg*(1-coef) + cflx_west*coef
          pflx_west_avg = pflx_west_avg*(1-coef) + pflx_west*coef
          ub_west_avg   =   ub_west_avg*(1-coef) +   ub_west*coef
        endif

        if (obc_east.and.inode.eq.npx-1) then
          cflx_east_avg = cflx_east_avg*(1-coef) + cflx_east*coef
          pflx_east_avg = pflx_east_avg*(1-coef) + pflx_east*coef
          ub_east_avg   =   ub_east_avg*(1-coef) +   ub_east*coef
        endif
      endif

      end subroutine calc_tune_avg !]
! ----------------------------------------------------------------------
      subroutine write_obc_tune ![
      implicit none

      !local
      integer            :: ncid,ierr
      character(len=99)  :: fname
      save fname
      
      if (record==nrpf) then
        call create_obc_tune_file(fname)
        record = 0
      endif
      record = record + 1

      ierr=nf90_open(fname, nf90_write, ncid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'write obc_tune: ',fname)

      call ncwrite(ncid,'ocean_time',(/time/),(/record/))

      ! fluxes and ub coefficients are defined as nx, ny sized arrays
      ! so use method 2 for output (see roms_read_write)
      if (obc_south.and.(jnode==0)) then
        call ncwrite(ncid,'cf_south',cflx_south_avg(i0:i1),(/1,record/))
        call ncwrite(ncid,'pf_south',pflx_south_avg(i0:i1),(/1,record/))
        call ncwrite(ncid,'ub_south',  ub_south_avg(i0:i1),(/1,record/))
      endif
      if (obc_north.and.(jnode==npy-1)) then
        call ncwrite(ncid,'cf_north',cflx_north_avg(i0:i1),(/1,record/))
        call ncwrite(ncid,'pf_north',pflx_north_avg(i0:i1),(/1,record/))
        call ncwrite(ncid,'ub_north',  ub_north_avg(i0:i1),(/1,record/))
      endif
      if (obc_east.and.(inode==npx-1)) then
        call ncwrite(ncid,'cf_east',cflx_east_avg(j0:j1),(/1,record/))
        call ncwrite(ncid,'pf_east',pflx_east_avg(j0:j1),(/1,record/))
        call ncwrite(ncid,'ub_east',  ub_east_avg(j0:j1),(/1,record /))
      endif
      if (obc_west.and.(inode==0)) then
        call ncwrite(ncid,'cf_west',cflx_west_avg(j0:j1),(/1,record/))
        call ncwrite(ncid,'pf_west',pflx_west_avg(j0:j1),(/1,record/))
        call ncwrite(ncid,'ub_west',  ub_west_avg(j0:j1),(/1,record/))
      endif
      
      ierr=nf90_close(ncid)    

      navg = 0

      if (mynode == 0) then
         write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &    'orlanski_tune :: wrote output, tdays =', tdays,
     &    'step =', iic-1, 'rec =', record
        endif


      end subroutine write_obc_tune !]
! ----------------------------------------------------------------------
      subroutine create_obc_tune_file(fname) ![
      implicit none

      !input/output
      character(len=*),intent(out) :: fname

      !local
      integer :: ncid,ierr,varid

      call create_file('_tun',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      varid = nccreate(ncid,'cf_south',(/dn_xr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'South boundary child flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'pf_south',(/dn_xr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'South boundary parent flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'ub_south',(/dn_xr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'South boundary binding velocity')
      ierr = nf90_put_att(ncid,varid,'units','m/s' )

      varid = nccreate(ncid,'cf_north',(/dn_xr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'North boundary child flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'pf_north',(/dn_xr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'North boundary parent flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'ub_north',(/dn_xr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'North boundary binding velocity')
      ierr = nf90_put_att(ncid,varid,'units','m/s' )

      varid = nccreate(ncid,'cf_west',(/dn_yr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'West boundary child flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'pf_west',(/dn_yr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'West boundary parent flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'ub_west',(/dn_yr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'West boundary binding velocity')
      ierr = nf90_put_att(ncid,varid,'units','m/s' )

      varid = nccreate(ncid,'cf_east',(/dn_yr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'East boundary child flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'pf_east',(/dn_yr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'East boundary parent flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m' )

      varid = nccreate(ncid,'ub_east',(/dn_yr,dn_tm/),(/xi_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                            ,'East boundary binding velocity')
      ierr = nf90_put_att(ncid,varid,'units','m/s' )

      ierr = nf90_close(ncid)

      end subroutine create_obc_tune_file !]
! ----------------------------------------------------------------------
      subroutine wrt_rst_tune(ncid,record)  ![
      ! Write the bc tuning coefficients to the restart file

      implicit none
      ! import/export
      integer,intent(in) :: ncid
      integer,intent(in) :: record

      if (tune_init) call init_orlanski_tune

      if (obc_south.and.(jnode==0)) then
        call ncwrite(ncid,'ub_south',ub_south(1:i1),(/1,record/))
      endif
      if (obc_north.and.(jnode==npy-1)) then
        call ncwrite(ncid,'ub_north',ub_north(1:i1),(/1,record/))
      endif
      if (obc_east.and.(inode==npx-1)) then
        call ncwrite(ncid,'ub_east' ,ub_east(1:j1) ,(/1,record/))
      endif
      if (obc_west.and.(inode==0)) then
        call ncwrite(ncid,'ub_west' ,ub_west(1:j1) ,(/1,record/))
      endif

      end subroutine wrt_rst_tune !]
! ----------------------------------------------------------------------
      subroutine get_init_tune(ncid,record)  ![
      ! get initial ub bry coupling coefficients
      implicit none

      ! input
      integer,intent(in) :: ncid,record

      ! local
      integer :: i,j
      integer :: ierr, varid

      ! If you want to use obc tuning, you'll need to compute
      ! the diagnostic pressure fluxes, or: diag_pflx=.true.
      ! in diagnostics.opt

      if (.not.diag_pflx) then
        if (mynode.eq.0) then
          print *, 'For obc tuning, ',
     &     'you must have diag_pflx=.true.'
        endif
      endif

      if (mynode==0) print *,'getting ub coefficients'

      if (obc_south.and.(jnode==0)) then
        if (.not. allocated(ub_south)) then
           allocate(ub_south(0:nx+1))
        endif
        ierr=nf90_inq_varid (ncid, 'ub_south', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_south', ub_south(i0:i1),(/1,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_south'
     &      ,' not in initial file, setting to one.'
          ub_south = 1
        endif
      endif
      if (obc_north.and.(jnode==npy-1)) then
        if (.not. allocated(ub_north)) then
           allocate(ub_north(0:nx+1))
        endif
        ierr=nf90_inq_varid (ncid, 'ub_north', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_north', ub_north(i0:i1),(/1,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_north'
     &      ,' not in initial file, setting to one.'
          ub_north = 1
        endif
      endif

      if (obc_east.and.(inode==npx-1)) then
        if (.not. allocated(ub_east)) then
           allocate(ub_east(0:ny+1))
        endif
        ierr=nf90_inq_varid (ncid, 'ub_east', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_east', ub_east(i0:j1),(/1,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_east'
     &      ,' not in initial file, setting to one.'
          ub_east = 1
        endif
      endif

      if (obc_west.and.(inode==0)) then
        if (.not. allocated(ub_west)) then
           allocate(ub_west(0:ny+1))
        endif
        ierr=nf90_inq_varid (ncid, 'ub_west', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_west', ub_west(j0:j1),(/1,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_west'
     &      ,' not in initial file, setting to one.'
          ub_west = 1
        endif

      endif
      
      end subroutine get_init_tune  !]

      end module obc_tune
