      module random_output
      ! Collection of random variables for output

#include "cppdefs.opt"

      use param
      use dimensions
      use roms_read_write
      use nc_read_write
      use netcdf
      use scalars
      use read_write
      use ocean_vars
      use diagnostics
      implicit none

      private

#include "random_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation

      ! Public functions
      public wrt_random

      contains


!----------------------------------------------------------------------
      subroutine def_vars_random(ncid)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer                        :: ierr, varid

      varid = nccreate(ncid,'zeta',(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','sea surface heigth')
      ierr = nf90_put_att(ncid,varid,'units','m')

      varid = nccreate(ncid,'u_surf',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','surface x velocity')
      ierr = nf90_put_att(ncid,varid,'units','m/s')

      varid = nccreate(ncid,'v_surf',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','surface y velocity')
      ierr = nf90_put_att(ncid,varid,'units','m/s')

      end subroutine def_vars_random  !]
!----------------------------------------------------------------------
      subroutine wrt_random  ![
      ! Call wrt_random after completion of the time-step
      ! (After step3d_uv2)
      implicit none

      ! local
      character(len=99),save :: fname
      integer,dimension(3)   :: start
      integer                :: ncid,ierr
      logical,save           :: first_step=.true.
      
      ! do nothing if do_random is false
      if (.not.do_random) return

      if (first_step) then
        first_step=.false.
      else
        output_time = output_time + dt                           ! only start count after first timestep
      endif

      if (output_time>=output_period) then

        if (record==nrpf) then
          call create_file('_rnd',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_random(ncid)
          ierr = nf90_enddef(ncid)
          record = 0
        endif
        record = record+1

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'wrt_random, opening: ',fname)

!       ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)     ! set fill value - nf90_nofill for optimized writing

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))


        start=1; start(3)=record
        call ncwrite(ncid,'zeta'  ,zeta(i0:i1,j0:j1,knew),start)
        call ncwrite(ncid,'u_surf',u( 1:i1,j0:j1,nz,nnew),start)
        call ncwrite(ncid,'v_surf',v(i0:i1, 1:j1,nz,nnew),start)

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_random :: wrote random, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        output_time=0
      endif

      end subroutine wrt_random !]
!----------------------------------------------------------------------
      subroutine create_file_random(fname)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      ! local
      integer :: ierr,varid,did,ncid
      character(len=10),dimension(4) :: dimnames           ! dimension names
      integer,          dimension(4) :: dimsizes

      fname=trim(adjustl(output_root_name)) / / '_rnd'
      call append_date_node(fname)

      ierr=nf90_create(trim(fname),nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) 
     &   call handle_ierr(ierr,'create_random_file: ')

      varid = nccreate(ncid,'ocean_time',(/'time'/),(/0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','time since 2000/1/1')
      ierr = nf90_put_att(ncid,varid,'units','second' )

      call put_global_atts(ncid, ierr)                     ! put global attributes in file
      if (ierr/=nf90_noerr) call handle_ierr(ierr,'create_random: ')

      call def_vars_random(ncid)

      ierr = nf90_enddef(ncid)
      if (mynode == 0) then
        write(*,'(7x,2A)')
     &    'wrt_random :: created new netCDF file ', trim(fname)
      endif

      end subroutine create_file_random !]

!-----------------------------------------------------------------------

      end module random_output
