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
      use surf_flux
      use roms_mpi
      implicit none

      private

#include "random_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation
      integer :: navg = 0      ! number of samples in average


      real,allocatable, dimension(:,:,:,:) :: u2_bdg,v2_bdg,w2_bdg,ek_bdg  ! total ek
      real,allocatable, dimension(:,:,:,:) :: ek_bdg_avg
      real, allocatable,dimension(:,:,:)   :: Wvl        ! Vertical velocity
      real, allocatable,dimension(:,:,:)   :: wb,wb_avg  ! Release of pot energy


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
      allocate(ek_bdg(0:nx+1,0:ny+1,nz,nmd+1))
      allocate(ek_bdg_avg(0:nx+1,0:ny+1,nz,nmd+1))
      ek_bdg_avg = 0

      allocate(u2_bdg(GLOBAL_2D_ARRAY,nz,nmd))
      allocate(v2_bdg(GLOBAL_2D_ARRAY,nz,nmd))
      allocate(w2_bdg(GLOBAL_2D_ARRAY,0:nz,1))

      allocate(Wvl(GLOBAL_2D_ARRAY,0:nz))
      allocate(rho_d(nx,ny,nz))
      rho_d = 0

      allocate(wb(GLOBAL_2D_ARRAY,1:nz))
      allocate(wb_avg(GLOBAL_2D_ARRAY,1:nz))
      wb_avg = 0

      if (mynode==0) print *,'init random'

      end subroutine init_random  !]
!----------------------------------------------------------------------
      subroutine calc_random ![
      ! Calculate variables for random output
      implicit none

      ! local
      integer      :: i,j,k,itrm
   
      ! Hack!! to extract wind work from the vmx term
      ! Move windstress to the hmx term (5) (which is mostly zero)
      Udiag(1:nx,1:ny,nz,6) = Udiag(1:nx,1:ny,nz,6) - sustr(1:nx,1:ny)
      Vdiag(1:nx,1:ny,nz,6) = Vdiag(1:nx,1:ny,nz,6) - svstr(1:nx,1:ny)
      Udiag(1:nx,1:ny,nz,5) = Udiag(1:nx,1:ny,nz,5) + sustr(1:nx,1:ny)
      Vdiag(1:nx,1:ny,nz,5) = Vdiag(1:nx,1:ny,nz,5) + svstr(1:nx,1:ny)

!     Full ke budget
      u2_bdg = 0
      v2_bdg = 0

      do k = 1,nz
        do itrm = 1,nmd
          do j = 1,ny
            do i = 1,nx
!             u2_bdg(i,j,k,itrm) = 
!    &              u(i,j,k,nrhs)*Udiag(i,j,k,itrm)
!             v2_bdg(i,j,k,itrm) = 
!    &              v(i,j,k,nrhs)*Vdiag(i,j,k,itrm)
              u2_bdg(i,j,k,itrm) = Udiag(i,j,k,itrm)*
     &            0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))
              v2_bdg(i,j,k,itrm) = Vdiag(i,j,k,itrm)*
     &            0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))
            enddo
          enddo
        enddo
      enddo
!     if (mynode.eq.0) then
!       print *,'u local: ',u(11,16,16,nrhs),Udiag(11,16,16,2)
!       print *,'v local: ',v(10,16,16,nrhs),Vdiag(10,16,16,2)
!     endif
!     stop

      w2_bdg = 0
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            w2_bdg(i,j,k,1) =
     &      (Wi(i,j,k)+We(i,j,k))*Wdiag(i,j,k,1)
          enddo
        enddo
      enddo

      do itrm=1,nmd
        call exchange_xxx(u2_bdg(:,:,:,itrm),v2_bdg(:,:,:,itrm))
      enddo

!     call exchange_xxx(w2_bdg(:,:,:,1))

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ek_bdg(i,j,k,1:8) = 0.5*(u2_bdg(i,j,k,:)+u2_bdg(i+1,j,k,:))
     &                        + 0.5*(v2_bdg(i,j,k,:)+v2_bdg(i,j+1,k,:))
          enddo
        enddo
      enddo
!     if (mynode.eq.0) then
!       print *,'ek local: ',ek_bdg(10,16,16,1),ek_bdg(10,16,16,8)
!       print *,'ek local: ',v2_bdg(10,16,16,1),v2_bdg(10,16,16,8)
!     endif

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ek_bdg(i,j,k,9) = 0.5*(w2_bdg(i,j,k-1,1)+w2_bdg(i,j,k,1))
          enddo
        enddo
      enddo

      call wvlcty (0,Wvl)
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            wb(i,j,k) = g*Wvl(i,j,k)*rho_d(i,j,k)*Hz(i,j,k)
          enddo
        enddo
      enddo

      

      end subroutine calc_random !]
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

      ek_bdg_avg = ek_bdg_avg*(1-coef) + ek_bdg*coef
      wb_avg = wb_avg*(1-coef) + wb*coef

      end subroutine calc_average !]
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

      varid = nccreate(ncid,'Ek_prg',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_cor',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_adv',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_dis',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_hmx',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_vmx',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_cpl',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_zxb',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_dzb',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

      varid = nccreate(ncid,'Ek_wb',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &              (/xi_rho,eta_rho,nz,nmd,0/))
      ierr = nf90_put_att(ncid,varid,'long_name','Kinetic energy budget terms')
      ierr = nf90_put_att(ncid,varid,'units','m3/s3')

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
      

      call calc_random
      call calc_average
      
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

        call ncwrite(ncid,'Ek_prg',ek_bdg_avg(i0:i1,j0:j1,:,1),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_cor',ek_bdg_avg(i0:i1,j0:j1,:,2),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_adv',ek_bdg_avg(i0:i1,j0:j1,:,3),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_dis',ek_bdg_avg(i0:i1,j0:j1,:,4),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_hmx',ek_bdg_avg(i0:i1,j0:j1,:,5),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_vmx',ek_bdg_avg(i0:i1,j0:j1,:,6),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_cpl',ek_bdg_avg(i0:i1,j0:j1,:,7),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_zxb',ek_bdg_avg(i0:i1,j0:j1,:,8),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_dzb',ek_bdg_avg(i0:i1,j0:j1,:,9),(/1,1,1,record/))
        call ncwrite(ncid,'Ek_wb' ,wb_avg(i0:i1,j0:j1,:),(/1,1,1,record/))

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_random :: wrote random, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        output_time=0
        navg =0
      endif

      end subroutine wrt_random !]
!----------------------------------------------------------------------

      end module random_output
