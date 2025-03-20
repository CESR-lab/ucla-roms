      module cstar_output
      ! Collection of variables dedicated to the C-star system
      ! Jeroen Molemaker Feb 2025

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

#include "cstar_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation
      integer :: navg = 0   

      ! Needed for averaging
      real,allocatable,dimension(:,:) :: zeta_avg
      real,allocatable,dimension(:,:) :: us_avg
      real,allocatable,dimension(:,:) :: vs_avg

      ! Public functions
      public wrt_cstar,init_cstar

      contains
!----------------------------------------------------------------------
      subroutine init_cstar ![
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

      if (do_avg) then
        allocate(zeta_avg(GLOBAL_2D_ARRAY) )
        allocate(us_avg(GLOBAL_2D_ARRAY) )
        allocate(vs_avg(GLOBAL_2D_ARRAY) )
      endif

      if (mynode==0) print *,'init cstar'

      end subroutine init_cstar  !]
!----------------------------------------------------------------------
      subroutine calc_average ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      use param
      implicit none

      ! local
      real :: coef

      navg = navg+1

      coef = 1./navg

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'cstar :: started averaging. ',
     &   'output_period (s) =', output_period
      endif

      zeta_avg = zeta_avg*(1-coef) + zeta(:,:,knew)*coef
      us_avg = us_avg*(1-coef) + u(:,:,nz,nnew)*coef
      vs_avg = vs_avg*(1-coef) + v(:,:,nz,nnew)*coef

      end subroutine calc_average !]
!----------------------------------------------------------------------
      subroutine def_vars_cstar(ncid)  ![
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

      if (do_avg) then
        ierr=nf90_put_att(ncid,nf90_global,'averaging',
     &      'All variables are averaged in time')   
      endif

      end subroutine def_vars_cstar  !]
!----------------------------------------------------------------------
      subroutine wrt_cstar  ![
      ! Call wrt_cstar after completion of the time-step
      implicit none

      ! local
      character(len=99),save :: fname
      integer,dimension(3)   :: start
      integer                :: ncid,ierr

      if (do_avg) call calc_average

      output_time = output_time + dt

      if (output_time>=output_period) then

        if (record==nrpf) then
          call create_file('_cstar',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_cstar(ncid)
          ierr = nf90_close(ncid)
          record = 0
        endif

        record = record+1

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'wrt_cstar, opening: ',fname)

!       ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)     ! set fill value - nf90_nofill for optimized writing

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        if (do_avg) then
          call ncwrite(ncid,'zeta'  ,zeta_avg(i0:i1,j0:j1),(/1,1,record/))
          call ncwrite(ncid,'u_surf',us_avg( 1:i1,j0:j1),(/1,1,record/))
          call ncwrite(ncid,'v_surf',vs_avg(i0:i1, 1:j1),(/1,1,record/))
        else
          call ncwrite(ncid,'zeta'  ,zeta(i0:i1,j0:j1,knew),(/1,1,record/))
          call ncwrite(ncid,'u_surf',u( 1:i1,j0:j1,nz,nnew),(/1,1,record/))
          call ncwrite(ncid,'v_surf',v(i0:i1, 1:j1,nz,nnew),(/1,1,record/))
        endif

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_cstar :: wrote cstar, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        output_time=0
        navg = 0
      endif
!
      end subroutine wrt_cstar !]
!----------------------------------------------------------------------

      end module cstar_output
