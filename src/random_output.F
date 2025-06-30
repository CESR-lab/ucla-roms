      module random_output
      ! Collection of random variables for output

#include "cppdefs.opt"

      use param
      use dimensions
      use roms_read_write
      use nc_read_write
      use netcdf
      use scalars
      use ocean_vars
      use diagnostics
      implicit none

      private

#include "random_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation

      ! Public functions
      public wrt_random,init_random

      contains

!----------------------------------------------------------------------
      subroutine init_random ![
      ! Allocate and initialize arrays.
      implicit none

      ! local
      logical,save :: done=.false.

      if (done) then
        return
      else
        done = .true.
      endif

      ! put the relevant part of your code here


      if (mynode==0) print *,'init random'

      end subroutine init_random  !]
!----------------------------------------------------------------------
      subroutine def_vars_random(ncid)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer                        :: ierr, varid

      varid = nccreate(ncid,'zeta',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','sea surface heigth')
      ierr = nf90_put_att(ncid,varid,'units','m')

      varid = nccreate(ncid,'u_surf',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','surface x velocity')
      ierr = nf90_put_att(ncid,varid,'units','m/s')

      varid = nccreate(ncid,'v_surf',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/))
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

      output_time = output_time + dt

      if (output_time>=output_period) then

        if (record==nrpf) then
          call create_file('_rnd',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_random(ncid)
          ierr = nf90_close(ncid)
          record = 0
        endif

        record = record+1

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr)
     &      call handle_ierr(ierr,'wrt_random, opening: ',fname)

!       ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)     ! set fill value - nf90_nofill for optimized writing

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        call ncwrite(ncid,'zeta'  ,zeta(i0:i1,j0:j1,knew),(/1,1,record/))
        call ncwrite(ncid,'u_surf',u( 1:i1,j0:j1,nz,nnew),(/1,1,record/))
        call ncwrite(ncid,'v_surf',v(i0:i1, 1:j1,nz,nnew),(/1,1,record/))

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

      end module random_output
