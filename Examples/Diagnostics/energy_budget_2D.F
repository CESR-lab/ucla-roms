      module energy_budget
      ! Compute kinetic energy budget
      ! Vertically integrated version

#include "cppdefs.opt"

      use param
      use dimensions
      use roms_read_write
      use nc_read_write
      use netcdf
      use scalars
      use surf_flux
      use tracers
      use diagnostics
      use grid
      use ocean_vars
      implicit none

      private

#include "energy_budget.opt"

      real    :: output_time  = 0
      integer :: record  = nrpf ! to trigger the first file creation
      integer :: navg = 0      ! number of samples in average



      real,allocatable, dimension(:,:,:) :: u2_bdg,v2_bdg,w2_bdg
      real,allocatable, dimension(:,:,:) :: ek_bdg,ek_bdg_avg 

      real,allocatable, dimension(:,:) :: u_hlf,v_hlf,uf_hlf,vf_hlf,om
!     real,allocatable, dimension(:,:) :: ustr,vstr


      ! Public functions
      public wrt_energy_budget,init_energy_budget

      contains
!----------------------------------------------------------------------
      subroutine init_energy_budget ![
      ! Allocate and initialize arrays.
      implicit none

      ! local
      logical,save :: done=.false.

      if (.not.diag_uv) then
        stop 'diag_uv is needed for energy budgets'
      endif

      if (done) then ! already done
        if (mynode.eq.0) print *, 'init energy_budget already done'
        return
      else
        done = .true.
      endif

      allocate(u_hlf(nx+1,ny+1))
      allocate(v_hlf(nx+1,ny+1))
      allocate(uf_hlf(nx+1,ny+1))
      allocate(vf_hlf(nx+1,ny+1))
      allocate(om(nx,ny))

!     allocate(ustr(nx+1,ny+1))
!     allocate(vstr(nx+1,ny+1))

      allocate(u2_bdg(GLOBAL_2D_ARRAY,nmd))
      allocate(v2_bdg(GLOBAL_2D_ARRAY,nmd))
      allocate(w2_bdg(GLOBAL_2D_ARRAY,1))
      u2_bdg = 0
      v2_bdg = 0
      w2_bdg = 0

      allocate(ek_bdg(GLOBAL_2D_ARRAY,nmd))
      ek_bdg = 0
      if (output_avg) then
        allocate(ek_bdg_avg(GLOBAL_2D_ARRAY,nmd))
        ek_bdg_avg = 0
      endif

      if (mynode.eq.0) print *, 'initialized energy_budget'

      end subroutine init_energy_budget  !]
!----------------------------------------------------------------------
      subroutine calc_energy_budget ![
      ! Calculate variables for energy budget

      ! Special emphasis is directed at the conversion between kinetic and 
      ! potential energy. Traditionally this is diagnosed as w*b, where w is
      ! the 'true' vertical velocity w= omega + u*z_x + v*z_y (see wvlcty.F) 
      ! We now recognize that wb should be diagnosed as:
      ! wb = om*dz*b + u*dz_x*b, where the buoyancy terms are taken from the
      ! pressure gradient subroutine to be consistent with the model computations
      ! the terms dz*b and z_x*b are stored in Wdiag, and U/Vdiag

      ! Barotropic/baroclinic omega:
      ! om_bt = z_t*(zr+h)/(zeta+h) ! dz/dt is computed in diagnostic.F
      ! om_bc = om - om_bt !! 
      ! the grid corrected omega is in fact the baroclinic omega

      implicit none

      ! local
      integer :: i,j,k,itrm

      ! Uflx, Vflx,dzdt are computed in diagnostics as part of the pflux
      ! u(nrhs) is about at n+1/2, but 0.5*(u(n)+u(n+1)) is more accurate
      ! nstp -> n
      ! nnew -> n+1

      ! Hack!! to extract wind work from the vmx term 
      ! Move windstress to the hmx term (5) (which is mostly zero)

      ! No need to remove dz factor from sustr/svstr since we work here
      ! with dz*diag terms (non-normalized)

      Udiag(1:nx,1:ny,nz,6) = Udiag(1:nx,1:ny,nz,6) - sustr(1:nx,1:ny)
      Vdiag(1:nx,1:ny,nz,6) = Vdiag(1:nx,1:ny,nz,6) - svstr(1:nx,1:ny)
      Udiag(1:nx,1:ny,nz,5) = Udiag(1:nx,1:ny,nz,5) + sustr(1:nx,1:ny)
      Vdiag(1:nx,1:ny,nz,5) = Vdiag(1:nx,1:ny,nz,5) + svstr(1:nx,1:ny)

      ! KE budget
      om = 0.0
      u2_bdg = 0.0
      v2_bdg = 0.0
      w2_bdg = 0.0
      do k = 1,nz
        do j = 1,ny+1
          do i = 1,nx+1
            u_hlf(i,j) = 0.5*( u(i,j,k,nstp)+u(i,j,k,nnew) )
            v_hlf(i,j) = 0.5*( v(i,j,k,nstp)+v(i,j,k,nnew) )
            uf_hlf(i,j) = 0.5*u_hlf(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))*dn_u(i,j)
            vf_hlf(i,j) = 0.5*v_hlf(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))*dm_v(i,j)
          enddo
        enddo
        ! vertical flux in m/s
        do j = 1,ny
          do i = 1,nx
            om(i,j) = -(uf_hlf(i+1,j)-uf_hlf(i,j)+vf_hlf(i,j+1)-vf_hlf(i,j))
     &               *pm(i,j)*pn(i,j) + om(i,j)
          enddo
        enddo
        do itrm = 1,nmd
          do j = 1,ny
            do i = 1,nx
              u2_bdg(i,j,itrm) = u2_bdg(i,j,itrm) + u_hlf(i,j)*Udiag(i,j,k,itrm)
              v2_bdg(i,j,itrm) = v2_bdg(i,j,itrm) + v_hlf(i,j)*Vdiag(i,j,k,itrm)
            enddo
          enddo
        enddo

        ! omega == 0 at k=0 so ok to do om*wdiag from 1:nz
        ! Wdiag(:,:,nz) is already half the amplitude to we're good to add 
        ! it in to the vertical integral
        do itrm = 1,1
          do j = 1,ny
            do i = 1,nx
              w2_bdg(i,j,itrm) = w2_bdg(i,j,itrm) + om(i,j)*Wdiag(i,j,k,itrm)
            enddo
          enddo
        enddo

      enddo ! k-loop

      if (.not.ew_periodic.and.inode==0) u2_bdg(1,:,:) = 0
      if (.not.ns_periodic.and.jnode==0) v2_bdg(:,1,:) = 0

!     exchange_xxx assumes that arrays are global_2d_array in size
      do itrm=1,nmd
        call exchange_xxx(u2_bdg(:,:,itrm),v2_bdg(:,:,itrm))
      enddo

      do j = 1,ny
        do i = 1,nx
          ek_bdg(i,j,:) = 0.5*(u2_bdg(i,j,:)+u2_bdg(i+1,j,:))
     &                  + 0.5*(v2_bdg(i,j,:)+v2_bdg(i,j+1,:))
          ek_bdg(i,j,9) = ek_bdg(i,j,9) - w2_bdg(i,j,1)
        enddo
      enddo

      end subroutine calc_energy_budget  !]
!----------------------------------------------------------------------
      subroutine calc_average ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg=0 after an output of the average
      use param
      implicit none

      ! local
      real :: coef

      navg = navg+1

      coef = 1./navg

      ek_bdg_avg = ek_bdg_avg*(1-coef) + ek_bdg*coef

      end subroutine calc_average !]
!----------------------------------------------------------------------
      subroutine def_vars_energy(ncid)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer                        :: ierr, varid

      ierr = nf90_redef(ncid)

       varid = nccreate(ncid,'ek_pgr',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; pressure gradient')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_cor',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; coriolis term')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_adv',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; advection')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_dis',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; dissipation')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_hmx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; horizontal mixing')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_vmx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; vertical mixing')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_cpl',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; coupling')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_dzt',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; time varying dz sizes')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_omb',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name',
     &                 'ek balance; WB')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

      end subroutine def_vars_energy  !]
!----------------------------------------------------------------------
      subroutine wrt_energy_budget  ![
      ! Call wrt_energy_budg after completion of the time-step
      implicit none

      ! local
      character(len=99),save :: fname
      integer,dimension(3)   :: start
      integer                :: ncid,ierr
      
      ! for 2d energy budgets, we don't want the 
      ! normalization in diagnostics
      if (diag_norm) then
        stop 'turn off normalization in diagnostics'
      endif
      call calc_energy_budget
      if (output_avg) call calc_average

      output_time = output_time  + dt

      if (output_time>=output_period) then
        if (record==nrpf) then
          call create_file('_ekb',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_energy(ncid)
          ierr = nf90_close(ncid)
          record = 0
        endif
        record = record+1

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'wrt_energy_budg, opening: ',fname)

!       ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)     ! set fill value - nf90_nofill for optimized writing

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        if (output_avg) then
          ek_bdg = ek_bdg_avg
          navg = 0
        endif

!       if (mynode.eq.0) then
!         print *,'ed write: ',ek_bdg(12,12,20,1)
!       endif

        call ncwrite(ncid,'ek_pgr',ek_bdg(i0:i1,j0:j1,1),(/1,1,record/))
        call ncwrite(ncid,'ek_cor',ek_bdg(i0:i1,j0:j1,2),(/1,1,record/))
        call ncwrite(ncid,'ek_adv',ek_bdg(i0:i1,j0:j1,3),(/1,1,record/))
        call ncwrite(ncid,'ek_dis',ek_bdg(i0:i1,j0:j1,4),(/1,1,record/))
        call ncwrite(ncid,'ek_hmx',ek_bdg(i0:i1,j0:j1,5),(/1,1,record/))
        call ncwrite(ncid,'ek_vmx',ek_bdg(i0:i1,j0:j1,6),(/1,1,record/))
        call ncwrite(ncid,'ek_cpl',ek_bdg(i0:i1,j0:j1,7),(/1,1,record/))
        call ncwrite(ncid,'ek_dzt',ek_bdg(i0:i1,j0:j1,8),(/1,1,record/))
        call ncwrite(ncid,'ek_omb',ek_bdg(i0:i1,j0:j1,9),(/1,1,record/))

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_energy_budget :: wrote budget, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        output_time=0
      endif

      end subroutine wrt_energy_budget !]
!----------------------------------------------------------------------

      end module energy_budget
