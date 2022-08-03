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
      use surf_flux
      use grid
      implicit none

      private

#include "random_output.opt"

      real    :: output_time  = 0
      real    :: output_time2 = 0
      integer :: record  = nrpf ! to trigger the first file creation
      integer :: record2 = nrpf2
      integer :: navg = 0      ! number of samples in average


      ! Barotropic/Baroclinic energy terms
      ! --------------
      real,allocatable, dimension(:,:,:)  :: Ud2d,Vd2d 
      real,allocatable, dimension(:,:)    :: ub_slow,vb_slow
      real,allocatable, dimension(:,:)    :: ub_fast,vb_fast

      real,allocatable, dimension(:,:)    :: px_int,py_int

      real,allocatable, dimension(:,:,:,:):: Udiag_slow,Vdiag_slow 
      real,allocatable, dimension(:,:)    :: Udiag_fast,Vdiag_fast 

      real,allocatable, dimension(:,:) :: u2,v2,ek  ! baroclinic ek
      real,allocatable, dimension(:,:) :: ek_avg 

      real,allocatable, dimension(:,:,:):: u2_bt,v2_bt,ek_bt  !  barotropic ek balance
      real,allocatable, dimension(:,:,:):: u2_bc,v2_bc,ek_bc  !  baroclinic ek balance

      real,allocatable, dimension(:,:,:):: u2_bta,v2_bta,ek_bta !  barotropic ek balance
      real,allocatable, dimension(:,:,:):: u2_bca,v2_bca,ek_bca !  baroclinic ek balance
      real,allocatable, dimension(:,:,:):: u2_bts,v2_bts,ek_bts !  barotropic ek balance
      real,allocatable, dimension(:,:,:):: u2_bcs,v2_bcs,ek_bcs !  baroclinic ek balance

      real,allocatable, dimension(:,:,:):: ek_bt_avg,ek_bc_avg 
      real,allocatable, dimension(:,:,:):: ek_bta_avg,ek_bca_avg 
      real,allocatable, dimension(:,:,:):: ek_bts_avg,ek_bcs_avg 

      real,allocatable, dimension(:,:) :: dzt_slow,dzt_fast ! dzeta/dt, filtered and prime
      real,allocatable, dimension(:,:) :: udh,udh_slow,udh_fast ! ubar.grad(h), filtered and prime
      real,allocatable, dimension(:,:) :: uhx, vhy

      real, allocatable,dimension(:,:,:):: rho_slow
      real, allocatable,dimension(:,:)  :: rho_fast

      real, allocatable,dimension(:,:,:) :: Wvl   ! Vertical velocity
      real, allocatable,dimension(:,:,:) :: Wvl_slow   ! Vertical velocity
      real, allocatable,dimension(:,:) :: Wvl_fast  ! Vertical velocity
      real, allocatable,dimension(:,:) :: bt2bc    ! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bt2bc_avg! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bc2bc    ! Barotropic to baroclinic conversion
      real, allocatable,dimension(:,:) :: bc2bc_avg! Barotropic to baroclinic conversion

      real, allocatable,dimension(:,:) :: ape      ! Available potential energy
      real, allocatable,dimension(:,:) :: ape_avg  ! Available potential energy
      real, allocatable,dimension(:,:) :: drhodz   ! density gradient

      ! Public functions
      public wrt_random

      contains
!----------------------------------------------------------------------
      subroutine init_random ![
      ! Allocate and initialize arrays.
      implicit none

      ! local
      integer :: itrc


      allocate( p_int(nx,ny) )
      allocate( px_int(nx,ny) )
      allocate( py_int(nx,ny) )
      p_int = 0

      allocate( uhx(nx+1,ny) )
      allocate( vhy(nx,ny+1) )
      allocate( udh(nx,ny) )
      allocate( udh_slow(nx,ny) )
      udh_slow=0
      allocate( udh_fast(nx,ny) )

      allocate( dzt_slow(nx,ny) )
      dzt_slow = 0
      allocate( dzt_fast(nx,ny) )

      allocate( bt2bc(nx,ny) )
      allocate( bt2bc_avg(nx,ny) )
      bt2bc_avg = 0

      allocate( bc2bc(nx,ny) )
      allocate( bc2bc_avg(nx,ny) )
      bc2bc_avg = 0

      allocate( rho_d(nx,ny,nz))
      rho_d = 0
      allocate( rho_slow(nx,ny,nz))
      rho_slow = 0
      allocate( rho_fast(nx,ny))

      allocate( Wvl(GLOBAL_2D_ARRAY,nz))
      allocate( Wvl_slow(nx,ny,nz))
      Wvl_slow = 0;
      allocate( Wvl_fast(nx,ny))

      allocate(ub_slow(nx,ny))
      allocate(vb_slow(nx,ny))
      ub_slow = 0
      vb_slow = 0
      allocate(ub_fast(nx,ny))
      allocate(vb_fast(nx,ny))

      allocate(Ud2d(nx,ny,nmd))
      allocate(Vd2d(nx,ny,nmd))

      allocate(Udiag_slow(nx,ny,nz,nmd))
      allocate(Vdiag_slow(nx,ny,nz,nmd))
      Udiag_slow = 0
      Vdiag_slow = 0
      allocate(Udiag_fast(nx,ny))
      allocate(Vdiag_fast(nx,ny))

      allocate(u2(GLOBAL_2D_ARRAY))
      allocate(v2(GLOBAL_2D_ARRAY))
      allocate(u2_bc(GLOBAL_2D_ARRAY,nmd))
      allocate(v2_bc(GLOBAL_2D_ARRAY,nmd))
      allocate(u2_bt(GLOBAL_2D_ARRAY,nmd))
      allocate(v2_bt(GLOBAL_2D_ARRAY,nmd))

      allocate(ek_bt(nx,ny,nmd))
      allocate(ek_bc(nx,ny,nmd))
      allocate(ek(nx,ny))

      allocate(u2_bca(GLOBAL_2D_ARRAY,nmd))
      allocate(v2_bca(GLOBAL_2D_ARRAY,nmd))
      allocate(u2_bta(GLOBAL_2D_ARRAY,nmd))
      allocate(v2_bta(GLOBAL_2D_ARRAY,nmd))
      allocate(u2_bcs(GLOBAL_2D_ARRAY,nmd))
      allocate(v2_bcs(GLOBAL_2D_ARRAY,nmd))
      allocate(u2_bts(GLOBAL_2D_ARRAY,nmd))
      allocate(v2_bts(GLOBAL_2D_ARRAY,nmd))
      allocate(ek_bta(nx,ny,nmd))
      allocate(ek_bca(nx,ny,nmd))
      allocate(ek_bts(nx,ny,nmd))
      allocate(ek_bcs(nx,ny,nmd))
      allocate(ek_bta_avg(nx,ny,nmd))
      allocate(ek_bca_avg(nx,ny,nmd))
      allocate(ek_bts_avg(nx,ny,nmd))
      allocate(ek_bcs_avg(nx,ny,nmd))
      ek_bta_avg = 0
      ek_bca_avg = 0
      ek_bts_avg = 0
      ek_bcs_avg = 0

      allocate(ek_bt_avg(nx,ny,nmd))
      allocate(ek_bc_avg(nx,ny,nmd))
      allocate(ek_avg(nx,ny))
      ek_bt_avg = 0
      ek_bc_avg = 0
      ek_avg = 0

      allocate( ape(nx,ny) )
      allocate( ape_avg(nx,ny) )
      ape_avg = 0

      allocate( drhodz(nx,ny) )

      if (mynode==0) print *,'init random'

      end subroutine init_random  !]
!----------------------------------------------------------------------
      subroutine calc_random ![
      ! Calculate variables for random output
      implicit none

      ! local
      integer :: i,j,k,k1,k2,itrm
      real    :: Wbt,Wbc,Wbt_f,Wbc_f,Wbt_s,Wbc_s

      ! barotropic and baroclinic kinetic energy;
      ! high pass filtered and vertically averaged.
      ! Also averaged in time

      do j=1,ny
        do i=1,nx
          ub_slow(i,j) = beta*ub_slow(i,j) + alpha*Uflx(i,j)
          ub_fast(i,j) = Uflx(i,j) - ub_slow(i,j)
          vb_slow(i,j) = beta*vb_slow(i,j) + alpha*Vflx(i,j)
          vb_fast(i,j) = Vflx(i,j) - vb_slow(i,j)
        enddo
      enddo

      ! Hack!! to extract wind work from the vmx term 
      ! Move windstress to the hmx term (5) (which is mostly zero)
      Udiag(:,:,nz,6) = Udiag(:,:,nz,6) - sustr(1:nx,1:ny)
      Vdiag(:,:,nz,6) = Vdiag(:,:,nz,6) - svstr(1:nx,1:ny)
      Udiag(:,:,nz,5) = Udiag(:,:,nz,5) + sustr(1:nx,1:ny)
      Vdiag(:,:,nz,5) = Vdiag(:,:,nz,5) + svstr(1:nx,1:ny)

      Udiag_slow = beta*Udiag_slow + alpha*Udiag
      Vdiag_slow = beta*Vdiag_slow + alpha*Vdiag

!     u2,v2 are sized as global_2D_Arrays

      u2 = 0
      v2 = 0
      u2_bc = 0
      v2_bc = 0
      u2_bt = 0
      v2_bt = 0
      ! Uflx/Vflx come from diagnostics and must be in velocity units
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            u_fast(i,j) = u(i,j,k,nrhs)-Uflx(i,j)- u_slow(i,j,k)
            v_fast(i,j) = v(i,j,k,nrhs)-Vflx(i,j)- v_slow(i,j,k)
          enddo
        enddo
        do j = 1,ny
          do i = 1,nx
            u2(i,j) = u2(i,j) + 
     &       0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u_fast(i,j)*u_fast(i,j)
            v2(i,j) = v2(i,j) + 
     &       0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v_fast(i,j)*v_fast(i,j)
          enddo
        enddo
        do itrm = 1,nmd
          Udiag_fast = Udiag(:,:,k,itrm) - Udiag_slow(:,:,k,itrm)
          Vdiag_fast = Vdiag(:,:,k,itrm) - Vdiag_slow(:,:,k,itrm)
          do j = 1,ny
            do i = 1,nx
              u2_bc(i,j,itrm)= u2_bc(i,j,itrm)+
     &              u_fast(i,j)*Udiag_fast(i,j)
              v2_bc(i,j,itrm)= v2_bc(i,j,itrm)+ 
     &              v_fast(i,j)*Vdiag_fast(i,j)
              u2_bt(i,j,itrm)= u2_bt(i,j,itrm)+
     &              ub_fast(i,j)*Udiag_fast(i,j)
              v2_bt(i,j,itrm)= v2_bt(i,j,itrm)+ 
     &              vb_fast(i,j)*Vdiag_fast(i,j)
            enddo
          enddo
        enddo
      enddo ! do k=1,nz

      u2_bca = 0
      v2_bca = 0
      u2_bta = 0
      v2_bta = 0
      do k = 1,nz
        do itrm = 1,nmd
          do j = 1,ny
            do i = 1,nx
!             u2_bca(i,j,itrm)= u2_bca(i,j,itrm)+
!    &              (u(i,j,k,nrhs)-Uflx(i,j))*Udiag(i,j,k,itrm)
!             v2_bca(i,j,itrm)= v2_bca(i,j,itrm)+ 
!    &              (v(i,j,k,nrhs)-Vflx(i,j))*Vdiag(i,j,k,itrm)
              u2_bca(i,j,itrm)= u2_bca(i,j,itrm)+
     &              u(i,j,k,nrhs)*Udiag(i,j,k,itrm)
              v2_bca(i,j,itrm)= v2_bca(i,j,itrm)+ 
     &              v(i,j,k,nrhs)*Vdiag(i,j,k,itrm)
              u2_bta(i,j,itrm)= u2_bta(i,j,itrm)+
     &              Uflx(i,j)*Udiag(i,j,k,itrm)
              v2_bta(i,j,itrm)= v2_bta(i,j,itrm)+ 
     &              Vflx(i,j)*Vdiag(i,j,k,itrm)
            enddo
          enddo
        enddo
      enddo

      u2_bcs = 0
      v2_bcs = 0
      u2_bts = 0
      v2_bts = 0
      do k = 1,nz
        do itrm = 1,nmd
          do j = 1,ny
            do i = 1,nx
              u2_bcs(i,j,itrm)= u2_bcs(i,j,itrm)+
     &              u_slow(i,j,k)*Udiag_slow(i,j,k,itrm)
              v2_bcs(i,j,itrm)= v2_bcs(i,j,itrm)+ 
     &              v_slow(i,j,k)*Vdiag_slow(i,j,k,itrm)
              u2_bts(i,j,itrm)= u2_bts(i,j,itrm)+
     &              ub_slow(i,j)*Udiag_slow(i,j,k,itrm)
              v2_bts(i,j,itrm)= v2_bts(i,j,itrm)+ 
     &              vb_slow(i,j)*Vdiag_slow(i,j,k,itrm)
            enddo
          enddo
        enddo
      enddo


!     In order to use mpi_exchange, they all need to be
!     global_2d_array in size
!     call mpi_exchange(u2_bt,v2_bt,u2_bc,v2_bc,u2,v2)
      call exchange2d_2_tile(1,nx,1,ny,u2,v2)
      do itrm=1,nmd
        call exchange2d_4_tile(1,nx,1,ny,u2_bt(:,:,itrm),v2_bt(:,:,itrm)
     &                                  ,u2_bc(:,:,itrm),v2_bc(:,:,itrm))
        call exchange2d_4_tile(1,nx,1,ny,u2_bta(:,:,itrm),v2_bta(:,:,itrm)
     &                                  ,u2_bca(:,:,itrm),v2_bca(:,:,itrm))
        call exchange2d_4_tile(1,nx,1,ny,u2_bts(:,:,itrm),v2_bts(:,:,itrm)
     &                                  ,u2_bcs(:,:,itrm),v2_bcs(:,:,itrm))
      enddo


      do j = 1,ny
        do i = 1,nx
          ek_bt(i,j,:) = 0.5*(u2_bt(i,j,:)+u2_bt(i+1,j,:))
     &                 + 0.5*(v2_bt(i,j,:)+v2_bt(i,j+1,:))
          ek_bc(i,j,:) = 0.5*(u2_bc(i,j,:)+u2_bc(i+1,j,:))
     &                 + 0.5*(v2_bc(i,j,:)+v2_bc(i,j+1,:))
          ek(i,j)      = 0.5*(u2(i,j)+u2(i+1,j))
     &                 + 0.5*(v2(i,j)+v2(i,j+1))
          ek_bta(i,j,:) = 0.5*(u2_bta(i,j,:)+u2_bta(i+1,j,:))
     &                  + 0.5*(v2_bta(i,j,:)+v2_bta(i,j+1,:))
          ek_bca(i,j,:) = 0.5*(u2_bca(i,j,:)+u2_bca(i+1,j,:))
     &                  + 0.5*(v2_bca(i,j,:)+v2_bca(i,j+1,:))
          ek_bts(i,j,:) = 0.5*(u2_bts(i,j,:)+u2_bts(i+1,j,:))
     &                  + 0.5*(v2_bts(i,j,:)+v2_bts(i,j+1,:))
          ek_bcs(i,j,:) = 0.5*(u2_bcs(i,j,:)+u2_bcs(i+1,j,:))
     &                  + 0.5*(v2_bcs(i,j,:)+v2_bcs(i,j+1,:))
        enddo
      enddo

      ! 'Missing' terms from the KE balance: PE2KE.
      ! Barotropic to baroclinic energy flux: Wbt*b
      ! Wbt is the vertical velocity of the barotropic mode
      ! alp = ubar*hx + vbar*hy
      ! bet = z_t + ubar*zx + vbar*zy
      ! Wbt = alp + (zr + h)* (bet - alp)/(zeta + h)
      ! bt2bc = Wbt*b
      ! dz/dt is computed in diagnostic.F

      ! ubar.grad(h)
        uhx = Uflx(1:nx+1,1:ny)*(h(1:nx+1,1:ny)-h(0:nx,1:ny))*
     &        0.5*(pm(1:nx+1,1:ny)+pm(0:nx,1:ny))
        vhy = Vflx(1:nx,1:ny+1)*(h(1:nx,1:ny+1)-h(1:nx,0:ny))*
     &        0.5*(pn(1:nx,1:ny+1)+pn(1:nx,0:ny))
        udh = 0.5*(uhx(1:nx,:)+uhx(2:nx+1,:) +
     &             vhy(:,1:ny)+vhy(:,2:ny+1) )

      call wvlcty (0,Wvl)
      Wvl_slow = beta*Wvl_slow + alpha*Wvl(1:nx,1:ny,:)
      rho_slow = beta*rho_slow + alpha*rho_d(1:nx,1:ny,:)
      dzt_slow = beta*dzt_slow + alpha*dzt
      udh_slow = beta*udh_slow + alpha*udh

      udh_fast = udh - udh_slow
      dzt_fast = dzt - dzt_slow

      bt2bc = 0
      bc2bc = 0
      ape   = 0
      do k=1,nz
        rho_fast = rho_d(1:nx,1:ny,k) - rho_slow(:,:,k)
        Wvl_fast = Wvl(1:nx,1:ny,k)   - Wvl_slow(:,:,k)
        k2 = min(k+1,nz)
        k1 = max(k-1,1)
        drhodz = -0.5*(rho_slow(1:nx,1:ny,k2)-rho_slow(1:nx,1:ny,k1))/Hz(1:nx,1:ny,k)
        do j=1,ny
          do i=1,nx
            Wbt_f = udh_fast(i,j) + (dzt_fast(i,j)-udh_fast(i,j))*
     &         (z_r(i,j,k)+h(i,j))/(z_w(i,j,nz)-z_w(i,j,0))

            Wbc_f = Wvl_fast(i,j) - Wbt_f

            bt2bc(i,j) = bt2bc(i,j) + g*Wbt_f*rho_fast(i,j)  *Hz(i,j,k)
            bc2bc(i,j) = bc2bc(i,j) + g*Wbc_f*rho_fast(i,j)  *Hz(i,j,k)

            if (drhodz(i,j)<1e-10) drhodz(i,j) = 1e-10
            ape(i,j) = ape(i,j)
     &        + 0.5*rho_fast(i,j)*rho_fast(i,j)/drhodz(i,j)*Hz(i,j,k)/h(i,j)
          enddo
        enddo
      enddo
      

      ! Don't split in fast/slow, but write at fast interval
      ! (p_int, px_int,py_int,ubar,vbar,zeta,u_surf,v_surf)
      do j=1,ny
        do i=1,nx
          px_int(i,j) = sum(Udiag(i,j,:,ipgr))
          py_int(i,j) = sum(Vdiag(i,j,:,ipgr))
        enddo
      enddo

      end subroutine calc_random  !]
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

      ek_bt_avg = ek_bt_avg*(1-coef) + ek_bt*coef
      ek_bc_avg = ek_bc_avg*(1-coef) + ek_bc*coef
      ek_avg = ek_avg*(1-coef) + ek*coef

      bt2bc_avg = bt2bc_avg*(1-coef) + bt2bc*coef
      bc2bc_avg = bc2bc_avg*(1-coef) + bc2bc*coef

      ek_bta_avg = ek_bta_avg*(1-coef) + ek_bta*coef
      ek_bca_avg = ek_bca_avg*(1-coef) + ek_bca*coef
      ek_bts_avg = ek_bts_avg*(1-coef) + ek_bts*coef
      ek_bcs_avg = ek_bcs_avg*(1-coef) + ek_bcs*coef

      end subroutine calc_average !]
!----------------------------------------------------------------------
      subroutine def_vars_random(ncid,vers)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      integer,intent(in) :: vers
      ! local
      integer                        :: ierr, varid

      ierr = nf90_redef(ncid)
      if (vers==1) then
       varid = nccreate(ncid,'zeta',(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','sea surface heigth')
       ierr = nf90_put_att(ncid,varid,'units','m')
 
       varid = nccreate(ncid,'ubar',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','barotropic velocity')
       ierr = nf90_put_att(ncid,varid,'units','m/s')
 
       varid = nccreate(ncid,'vbar',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','barotropic velocity')
       ierr = nf90_put_att(ncid,varid,'units','m/s')
 
       varid = nccreate(ncid,'u_surf',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','surface x velocity')
       ierr = nf90_put_att(ncid,varid,'units','m/s')
 
       varid = nccreate(ncid,'v_surf',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','surface y velocity')
       ierr = nf90_put_att(ncid,varid,'units','m/s')
 
       varid = nccreate(ncid,'p_int',(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','vertical int of p')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
 
       varid = nccreate(ncid,'px_int',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','vertical int of px')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
 
       varid = nccreate(ncid,'py_int',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','vertical int of py')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
      endif

      if (vers==2) then
       varid = nccreate(ncid,'ek_bt',(/'xi_rho','eta_rho','nmd','time'/),(/xi_rho,eta_rho,nmd,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF barotrop balance')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_bc',(/'xi_rho','eta_rho','nmd','time'/),(/xi_rho,eta_rho,nmd,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF baroclin balance')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek',(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF baroclinic ek')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'bt2bc',(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF btrop to bclin')
       ierr = nf90_put_att(ncid,varid,'units','W/m^3')

       varid = nccreate(ncid,'bc2bc',(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF bclin to bclin')
       ierr = nf90_put_att(ncid,varid,'units','W/m^3')

       varid = nccreate(ncid,'ek_bta',(/'xi_rho','eta_rho','nmd','time'/),(/xi_rho,eta_rho,nmd,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF barotrop balance')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_bca',(/'xi_rho','eta_rho','nmd','time'/),(/xi_rho,eta_rho,nmd,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF baroclin balance')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_bts',(/'xi_rho','eta_rho','nmd','time'/),(/xi_rho,eta_rho,nmd,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF barotrop balance')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

       varid = nccreate(ncid,'ek_bcs',(/'xi_rho','eta_rho','nmd','time'/),(/xi_rho,eta_rho,nmd,0/))
       ierr = nf90_put_att(ncid,varid,'long_name','HF baroclin balance')
       ierr = nf90_put_att(ncid,varid,'units','m^2/s^3')

      endif

      end subroutine def_vars_random  !]
!----------------------------------------------------------------------
      subroutine wrt_random  ![
      ! Call wrt_random after completion of the time-step
      ! (After step3d_uv2)
      implicit none

      ! local
      character(len=99),save :: fname,fname2
      integer,dimension(3)   :: start
      integer                :: ncid,ierr
      logical,save           :: first_step=.true.
      
      ! do nothing if do_random is false
      if (.not.do_random) return

      if (first_step) then
        first_step=.false.
        call init_random
        output_time = output_time  - dt ! Hack, but needed because of the timing of the call within the time step
!       output_time2= output_time2 - dt                        
        return
      else
        output_time = output_time  + dt       ! only start count after first timestep
        output_time2= output_time2 + dt       ! only start count after first timestep
      endif

      call calc_random

      call calc_average


      if (output_time>=output_period) then
        if (record==nrpf) then
          call create_file('_rnd',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_random(ncid,1)
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
        call ncwrite(ncid,'ubar'  ,ubar(1:i1,j0:j1,knew),(/1,1,record/))
        call ncwrite(ncid,'vbar'  ,vbar(i0:i1,1:j1,knew),(/1,1,record/))
        call ncwrite(ncid,'u_surf',u( 1:i1,j0:j1,nz,nnew),(/1,1,record/))
        call ncwrite(ncid,'v_surf',v(i0:i1, 1:j1,nz,nnew),(/1,1,record/))

        call ncwrite(ncid,'p_int'  ,p_int,(/bfx,bfy,record/))

        call ncwrite(ncid,'px_int',px_int,(/1,bfy,record/))
        call ncwrite(ncid,'py_int',py_int,(/bfx,1,record/))

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_random :: wrote random, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        output_time=0
      endif

      if (output_time2>=output_period2) then
        if (record2==nrpf2) then
          call create_file('_rnd2',fname2)
          ierr=nf90_open(fname2,nf90_write,ncid)
          call def_vars_random(ncid,2)
          ierr = nf90_close(ncid)
          record2= 0
        endif
        record2 = record2+1

        ierr=nf90_open(fname2,nf90_write,ncid)

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record2/))

        call ncwrite(ncid,'ek',ek_avg(1:nx,1:ny),(/bfx,bfy,record2/))
        call ncwrite(ncid,'ek_bt',ek_bt_avg(1:nx,1:ny,:),(/bfx,bfy,1,record2/))
        call ncwrite(ncid,'ek_bc',ek_bc_avg(1:nx,1:ny,:),(/bfx,bfy,1,record2/))
        call ncwrite(ncid,'bt2bc',bt2bc_avg(1:nx,1:ny),(/bfx,bfy,record2/))
        call ncwrite(ncid,'bc2bc',bc2bc_avg(1:nx,1:ny),(/bfx,bfy,record2/))
        call ncwrite(ncid,'ek_bta',ek_bta_avg(1:nx,1:ny,:),(/bfx,bfy,1,record2/))
        call ncwrite(ncid,'ek_bca',ek_bca_avg(1:nx,1:ny,:),(/bfx,bfy,1,record2/))
        call ncwrite(ncid,'ek_bts',ek_bts_avg(1:nx,1:ny,:),(/bfx,bfy,1,record2/))
        call ncwrite(ncid,'ek_bcs',ek_bcs_avg(1:nx,1:ny,:),(/bfx,bfy,1,record2/))
        ierr=nf90_close (ncid)

        output_time2=0
        navg = 0
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_random :: wrote random, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif
      endif

      end subroutine wrt_random !]
!----------------------------------------------------------------------

      end module random_output
