      module sponge_tune

      ! Tunes the sponge values near the boundaries to attempt
      ! to match the incoming baroclinic wave flux from the parent

#include "cppdefs.opt"

      use netcdf
      use nc_read_write
      use scalars
      use dimensions
      use diagnostics
      use roms_read_write
      use hidden_mpi_vars
      use mixing

      implicit none
      private

#include "sponge_tune.opt"

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

      real,allocatable,dimension(:) :: cflx_west_avg
      real,allocatable,dimension(:) :: pflx_west_avg
      real,allocatable,dimension(:) ::   ub_west_avg


      public adjust_orlanski
      public init_orlanski_tune
      public wrt_rst_ub
      public get_init_ub

      contains

! ----------------------------------------------------------------------
      subroutine init_orlanski_tune ![
      ! Initialize sponge tuning
      implicit none

      ! local
      integer :: i

      if (.not.diag_pflx.and.mynode==0) then
        print *,'For Orlanski tuning,'
        print *,'pflx diagnostics must be turned on!'
        stop
      endif

      call store_string_att(bc_options,'UB_TUNING, ')


      !! max value based of cfl condition
!     sp_mx = 0.06/(maxval(pm)*maxval(pn)*dt)

      if (obc_west.and.(.not.west_exchng)) then
        allocate(nc_pflx_w%vdata(ny,1,2))
        allocate(pflx_west(ny))
        pflx_west = 0
        allocate(pflx_west_avg(ny))
        pflx_west_avg = 0
        allocate(cflx_west(ny))
        cflx_west = 0
        allocate(cflx_west_avg(ny))
        cflx_west_avg = 0
        if (.not.allocated(ub_west)) then
          allocate(ub_west(ny))
          ub_west = 0.
        endif
        allocate(ub_west_avg(ny))
        ub_west_avg = 0.
      endif
      if (obc_east.and.(.not.east_exchng)) then
        allocate(nc_pflx_e%vdata(ny,1,2))
        allocate(pflx_east(ny))
        pflx_east = 0
        allocate(cflx_east(ny))
        cflx_east = 0
        if (.not.allocated(ub_east)) then
          allocate(ub_east(ny))
          ub_east = 0.
        endif
      endif
      if (obc_north.and.(.not.north_exchng)) then
        allocate(nc_pflx_n%vdata(nx,1,2))
        allocate(pflx_north(nx))
        pflx_north = 0
        allocate(cflx_north(nx))
        cflx_north = 0
        if (.not.allocated(ub_north)) then
          allocate(ub_north(nx))
          ub_north = 0.
        endif
      endif
      if (obc_south.and.(.not.south_exchng)) then
        allocate(nc_pflx_s%vdata(nx,1,2))
        allocate(pflx_south(nx))
        pflx_south = 0
        allocate(pflx_south_avg(nx))
        pflx_south_avg = 0
        allocate(cflx_south(nx))
        cflx_south = 0
        allocate(cflx_south_avg(nx))
        cflx_south_avg = 0
        if (.not.allocated(ub_south)) then
          allocate(ub_south(nx))
          ub_south = 0.
        endif
        allocate(ub_south_avg(nx))
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
      ! Compute  baroclinic pressure fluxes
      ! up/vp fluxes are positive to the east/north
      ! Instead of taking the mean flux over the sponge region, we're
      ! taking the flux on the 'inside' boundary of the sponge region.

      ! local
      integer :: pos

      pos = 2
      if (obc_west.and.(.not.west_exchng)) then
        cflx_west = up(pos,1:ny)
      endif
      if (obc_east.and.(.not.east_exchng)) then
        cflx_east = up(nx-pos+1,1:ny)
      endif
      if (obc_south.and.(.not.south_exchng)) then
        cflx_south = vp(1:nx,pos)
      endif
      if (obc_north.and.(.not.north_exchng)) then
        cflx_north = vp(1:nx,ny-pos+1)
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

      eps = 1e-6
      ub_mx = 2.0
      ub_mn =-1.0

      alpha = 0.5*ub_mx*dt/sp_timscale

      if (tune_init) call init_orlanski_tune

      call set_pflx
      call comp_cflx

      if (obc_south.and.(.not.south_exchng)) then
        ub_south = ub_south+alpha*(pflx_south-cflx_south)
        ub_south = max(ub_mn,ub_south)
        ub_south = min(ub_mx,ub_south)
!       if (mynode==6) write(10,*),pflx_south(50),cflx_south(50),ub_south(50)
      endif

      if (obc_north.and.(.not.north_exchng)) then
        ub_north = ub_north-alpha*(pflx_north-cflx_north)
        ub_north = max(ub_mn,ub_north)
        ub_north = min(ub_mx,ub_north)
      endif

      if (obc_west.and.(.not.west_exchng)) then
        ub_west = ub_west+alpha*(pflx_west-cflx_west)
        ub_west = max(ub_mn,ub_west)
        ub_west = min(ub_mx,ub_west)
      endif

      if (obc_east.and.(.not.east_exchng)) then
        ub_east = ub_east-alpha*(pflx_east-cflx_east)
        ub_east = max(ub_mn,ub_east)
        ub_east = min(ub_mx,ub_east)
      endif

      output_time = output_time + dt
      call calc_spn_avg
      if (output_time>=output_period .and. wrt_sponge) then
        call write_sp_tune
        output_time = 0
      endif

      end subroutine adjust_orlanski !]
! ----------------------------------------------------------------------
      subroutine calc_spn_avg ![
      ! Update diagnostics averages
      ! The average is always scaled properly throughout
      ! reset navg_diag=0 after an output of the average
      use param
      implicit none

      ! local
      real :: coef

      navg = navg +1

      coef = 1./navg

      if (spn_avg) then
        if (obc_south.and.(.not.south_exchng)) then
          cflx_south_avg = cflx_south_avg*(1-coef) + cflx_south*coef
          pflx_south_avg = pflx_south_avg*(1-coef) + pflx_south*coef
          ub_south_avg   =   ub_south_avg*(1-coef) +   ub_south*coef
        endif

        if (obc_west.and.(.not.west_exchng)) then
          cflx_west_avg = cflx_west_avg*(1-coef) + cflx_west*coef
          pflx_west_avg = pflx_west_avg*(1-coef) + pflx_west*coef
          ub_west_avg   =   ub_west_avg*(1-coef) +   ub_west*coef
        endif
      endif

      end subroutine calc_spn_avg !]
! ----------------------------------------------------------------------
      subroutine write_sp_tune ![
      implicit none

      !local
      integer            :: ncid,ierr
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

      ! fluxes and ub coefficients are defined as nx, ny sized arrays
      ! so use method 2 for output (see roms_read_write)
      if (obc_south.and.(jnode==0)) then
        call ncwrite(ncid,'cf_south',cflx_south_avg,(/bfx,record /))
        call ncwrite(ncid,'pf_south',pflx_south_avg,(/bfx,record /))
        call ncwrite(ncid,'ub_south',  ub_south_avg,(/bfx,record /))
      endif
      if (obc_east.and.(inode==np_xi-1)) then
!       call ncwrite(ncid,'ub_east',ub_east(j0:j1),(/1,record /))
      endif
      if (obc_west.and.(inode==0)) then
        call ncwrite(ncid,'cf_west',cflx_west_avg,(/ bfy,record /))
        call ncwrite(ncid,'pf_west',pflx_west_avg,(/ bfy,record /))
        call ncwrite(ncid,'ub_west',  ub_west_avg,(/ bfy,record /))
      endif

      ierr=nf90_close(ncid)

      navg = 0

      if (mynode == 0) then
         write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &    'orlanski_tune :: wrote output, tdays =', tdays,
     &    'step =', iic-1, 'rec =', record
        endif


      end subroutine write_sp_tune !]
! ----------------------------------------------------------------------
      subroutine create_sp_tune_file(fname) ![
      implicit none

      !input/output
      character(len=*),intent(out) :: fname

      !local
      integer :: ncid,ierr,varid

      call create_file('_spn',fname)

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

      ierr = nf90_close(ncid)

      end subroutine create_sp_tune_file !]
! ----------------------------------------------------------------------
      subroutine wrt_rst_ub(ncid,record)  ![
      ! Write the bc tuning coefficients to the restart file

      implicit none
      ! import/export
      integer,intent(in) :: ncid
      integer,intent(in) :: record

      if (mynode==0) print *,'writing ub in restart file'
      if (tune_init) call init_orlanski_tune

      if (obc_south.and.(jnode==0)) then
        call ncwrite(ncid,'ub_south',ub_south,(/bfx,record/))
      endif
      if (obc_north.and.(jnode==np_eta-1)) then
        call ncwrite(ncid,'ub_north',ub_north,(/bfx,record/))
      endif
      if (obc_east.and.(inode==np_xi-1)) then
        call ncwrite(ncid,'ub_east' ,ub_east ,(/bfy,record/))
      endif
      if (obc_west.and.(inode==0)) then
        call ncwrite(ncid,'ub_west' ,ub_west ,(/bfy,record/))
      endif

      end subroutine wrt_rst_ub !]
! ----------------------------------------------------------------------
      subroutine get_init_ub(ncid,record)  ![
      ! get initial ub bry coupling coefficients
      implicit none

      ! input
      integer,intent(in) :: ncid,record

      ! local
      integer :: i,j
      integer :: ierr, varid

      if (mynode==0) print *,'getting ub coefficients'
      if (obc_south.and.(jnode==0)) then
        if (.not. allocated(ub_south)) allocate(ub_south(nx))
        ierr=nf90_inq_varid (ncid, 'ub_south', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_south', ub_south,(/bfx,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_south'
     &      ,' not in initial file, setting to zero.'
          ub_south = 0
        endif
      endif
      if (obc_north.and.(jnode==np_eta-1)) then
        if (.not. allocated(ub_north)) allocate(ub_north(nx))
        ierr=nf90_inq_varid (ncid, 'ub_north', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_north', ub_north,(/bfx,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_north'
     &      ,' not in initial file, setting to zero.'
          ub_north = 0
        endif
      endif

      if (obc_east.and.(inode==np_xi-1)) then
        if (.not. allocated(ub_east)) allocate(ub_east(ny))
        ierr=nf90_inq_varid (ncid, 'ub_east', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_east', ub_east,(/bfy,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_east'
     &      ,' not in initial file, setting to zero.'
          ub_east = 0
        endif
      endif

      if (obc_west.and.(inode==0)) then
        if (.not. allocated(ub_west)) allocate(ub_west(ny))
        ierr=nf90_inq_varid (ncid, 'ub_west', varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid,'ub_west', ub_west,(/bfy,record/))
        else
          if (mynode==0) print *,'--WARNING: ub_west'
     &      ,' not in initial file, setting to zero.'
          ub_west = 0
        endif

      endif

      end subroutine get_init_ub  !]

      end module sponge_tune
