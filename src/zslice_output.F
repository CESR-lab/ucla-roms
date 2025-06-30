      module zslice_output
      ! Collection of zslice variables for output

#include "cppdefs.opt"

      use param
      use grid
      use dimensions
      use roms_read_write
      use nc_read_write
      use netcdf
      use scalars
      use ocean_vars
      use diagnostics
      use tracers
      use surf_flux

      implicit none

      private

#include "zslice_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf    ! to trigger the first file creation
      integer :: navg = 0         ! number of samples in average
      real    :: FillValue=1.E+33 !

      ! Variables to z-slice
      real,allocatable,dimension(:,:,:,:) :: Tz, Tz_avg
      real,allocatable,dimension(:,:,:) :: Uz, Uz_avg
      real,allocatable,dimension(:,:,:) :: Vz, Vz_avg

      ! Needed for z-slicing
      real,allocatable,dimension(:,:) :: zz
      real,allocatable,dimension(:,:) :: var
      real,allocatable,dimension(:,:) :: var_zlv

      ! Public functions
      public wrt_zslice

      contains

!----------------------------------------------------------------------
      subroutine init_zslice ![
      ! Allocate and initialize arrays.
      implicit none

      integer m

      if (wrt_T_zsl)  then
          allocate( Tz(GLOBAL_2D_ARRAY,ndep,nt_z) )
          Tz=0.
          if (do_avg)  then
             allocate( Tz_avg(GLOBAL_2D_ARRAY,ndep,nt_z) )
             Tz_avg=0.
          endif
      endif

      if (wrt_U_zsl)  then
          allocate( Uz(GLOBAL_2D_ARRAY,ndep) )
          Uz=0.
          if (do_avg)  then
             allocate( Uz_avg(GLOBAL_2D_ARRAY,ndep) )
             Uz_avg=0.
          endif
      endif

      if (wrt_V_zsl)  then
          allocate( Vz(GLOBAL_2D_ARRAY,ndep) )
          Vz=0.
          if (do_avg)  then
             allocate( Vz_avg(GLOBAL_2D_ARRAY,ndep) )
             Vz_avg=0.
          endif
      endif

      allocate(zz(1:nx,0:N+1))
      allocate(var(1:nx,1:nz))
      allocate(var_zlv(nx,ndep))

      if (mynode==0) write(*,'(7x,A,I2,A,I2,A)')
     &    'init zslice : nb depth =',ndep , ' to be done on ',
     &    nt_z , ' tracers '
      if (mynode==0) then
         do m=1,ndep
         write(*,'(7x,A,I2,A,F7.2)') 'zslice : dep(',m,')=',vecdep(m)
         enddo
         do m=1,nt_z
         write(*,'(7x,A,I2,A,A)') 'zslice : tracers(',trc2zsc(m),')=',t_vname(trc2zsc(m))
         enddo
      endif

      end subroutine init_zslice  !]
!----------------------------------------------------------------------
      subroutine calc_zslice ![
      ! Calculate variables for zslice output
      ! inspired from zslice in Tools-Roms
      ! only the vertical RHO-point (kmin==1) are implemented
      !

      implicit none

      ! local
      integer i, j, k, m, trcz
      real zlev, dpth
      integer km(1:nx)

      if (wrt_T_zsl) then
      do j=1,ny
         do k=1,N
         do i=1,nx
            zz(i,k)=z_r(i,j,k)-z_w(i,j,N)
         enddo
         enddo
         do i=1,nx
            zz(i,N+1)=z_w(i,j,N)-z_w(i,j,N)
            zz(i,  0)=z_w(i,j,0)-z_w(i,j,N)
         enddo
         do trcz=1,nt_z
!         if (mynode==0 .and. j==10) print *, j, zz(10,0), zz(10,N), zz(10,N+1)
!         if (mynode==0 .and. j==10) print *, t(10,j,nz,nnew,trc2zsc(trcz))
            call sigma_to_z(var_zlv,t(1:nx,j,1:nz,nnew,trc2zsc(trcz)),zz,j)
            do i=1,nx
            do m=1,ndep
            Tz(i,j,m,trcz)=var_zlv(i,m)
            enddo
            enddo
          enddo
      enddo
      endif

      if (wrt_U_zsl) then
      do j=1,ny
         do k=1,N
         do i=1,nx
            zz(i,k)=0.5*(z_r(i,j,k)+z_r(i-1,j,k))-0.5*(z_w(i-1,j,N)+z_w(i,j,N))
         enddo
         enddo
         do i=1,nx
            zz(i,N+1)=0.5*(z_w(i-1,j,N)+z_w(i,j,N))-0.5*(z_w(i-1,j,N)+z_w(i,j,N))
            zz(i,  0)=0.5*(z_w(i-1,j,0)+z_w(i,j,0))-0.5*(z_w(i-1,j,N)+z_w(i,j,N))
         enddo
         call sigma_to_z(var_zlv,u(1:nx,j,1:nz,nnew),zz,j)
         do i=1,nx
         do m=1,ndep
         Uz(i,j,m)=var_zlv(i,m)
         enddo
         enddo
      enddo
      call exchange_xxx(Uz)

      endif

      if (wrt_V_zsl) then
      do j=1,ny
         do k=1,N
         do i=1,nx
            zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))-0.5*(z_w(i,j,N)+z_w(i,j-1,N))
         enddo
         enddo
         do i=1,nx
            zz(i,N+1)=0.5*(z_w(i,j,N)+z_w(i,j-1,N))-0.5*(z_w(i,j,N)+z_w(i,j-1,N))
            zz(i,  0)=0.5*(z_w(i,j,0)+z_w(i,j-1,0))-0.5*(z_w(i,j,N)+z_w(i,j-1,N))
         enddo
         call sigma_to_z(var_zlv,v(1:nx,j,1:nz,nnew),zz,j)
         do i=1,nx
         do m=1,ndep
         Vz(i,j,m)=var_zlv(i,m)
         enddo
         enddo
      enddo
      call exchange_xxx(Vz)
      endif

      end subroutine calc_zslice  !]
!----------------------------------------------------------------------
      subroutine sigma_to_z(var_zlv,var,zz,j) ![

      implicit none

      !import/export
      integer, intent(in) :: j
      real, dimension(1:nx,0:N+1), intent(in) :: zz
      real, dimension(:,:), intent(in) :: var
      real, dimension(:,:), intent(out) :: var_zlv

      ! local
      integer i, k, m
      real zlev, dpth
      integer km(1:nx)

!         if (mynode==0 .and. j==10) print *, j, zz(10,0), zz(10,N), zz(10,N+1)
!         if (mynode==0 .and. j==10) print *, var(10,nz)
         var_zlv=0

         do m=1,ndep

          zlev=vecdep(m)

          do i=1,nx
              dpth=zz(i,N+1)-zz(i,0)
              if (rmask(i,j) < 0.5) then
                km(i)=-3          !--> masked out
              elseif (dpth*(zlev-zz(i,N+1)) > 0.) then
                km(i)=N+2         !<-- above surface
              elseif (dpth*(zlev-zz(i,N)) > 0.) then
                km(i)=N           !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-zlev) > 0.) then
                km(i)=-2          !<-- below bottom
              elseif (dpth*(zz(i,1)-zlev) > 0.) then
                km(i)=0           !<-- above bottom, but below z_r(1)
              else
                km(i)=-1          !--> to search
              endif
          enddo

          do k=N-1,1,-1
            do i=1,nx
              if (km(i) == -1) then
                if ((zz(i,k+1)-zlev)*(zlev-zz(i,k)) >= 0.) then
                  km(i)=k
                endif
              endif
            enddo
          enddo

          do i=1,nx
            if (km(i) == -3) then
              var_zlv(i,m)=0.             !<-- masked out
            elseif (km(i) == -2) then
              var_zlv(i,m)=0.             !<-- below bottom
            elseif (km(i) == N+2) then
              var_zlv(i,m)=0.             !<-- above surface
            elseif (km(i) == N) then
              var_zlv(i,m)=var(i,N)       !-> R-point, above z_r(N)
     &              +(zlev-zz(i,N))*(var(i,N)-var(i,N-1))
     &                                   /(zz(i,N)-zz(i,N-1))
            elseif (km(i) == 0) then   !-> R-point below z_r(1),
              var_zlv(i,m)=var(i,1)  !     but above bottom
     &          -(zz(i,1)-zlev)*(var(i,2)-var(i,1))
     &                                  /(zz(i,2)-zz(i,1))
            else
              k=km(i)
              var_zlv(i,m)=( var(i,k)*(zz(i,k+1)-zlev)
     &                        +var(i,k+1)*(zlev-zz(i,k))
     &                               )/(zz(i,k+1)-zz(i,k))
            endif
          enddo

        enddo

      end subroutine sigma_to_z  !]
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
     &   'zslice :: started averaging. ',
     &   'output_period (s) =', output_period
      endif

      if (wrt_T_zsl)  Tz_avg = Tz_avg*(1-coef) + Tz*coef
      if (wrt_U_zsl)  Uz_avg = Uz_avg*(1-coef) + Uz*coef
      if (wrt_V_zsl)  Vz_avg = Vz_avg*(1-coef) + Vz*coef

      end subroutine calc_average !]
!----------------------------------------------------------------------
      subroutine def_vars_zslice(ncid)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer                        :: ierr, varid, n

      varid = nccreate(ncid,'depth',(/'depth'/),(/ndep/))

      if (wrt_T_zsl) then
      do n=1,nt_z
      varid = nccreate(ncid,t_vname(trc2zsc(n)),
     &         (/dn_xr,dn_yr,'depth',dn_tm/),(/xi_rho,eta_rho,ndep,0/))
      if (do_avg) then
        ierr = nf90_put_att(ncid,varid,'long_name','averaged'//t_lname(trc2zsc(n)))
      else
        ierr = nf90_put_att(ncid,varid,'long_name',t_lname(trc2zsc(n)))
      endif
      ierr = nf90_put_att(ncid,varid,'units',t_units(trc2zsc(n)))
      enddo
      endif

      if (wrt_U_zsl) then
      varid = nccreate(ncid,'u',(/dn_xu,dn_yr,'depth',dn_tm/),(/xi_u,eta_rho,ndep,0/))
      if (do_avg) then
        ierr = nf90_put_att(ncid,varid,'long_name','averaged u-momentum component')
      else
        ierr = nf90_put_att(ncid,varid,'long_name','u-momentum component')
      endif
      ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif

      if (wrt_V_zsl) then
      varid = nccreate(ncid,'v',(/dn_xr,dn_yv,'depth',dn_tm/),(/xi_rho,eta_v,ndep,0/))
      if (do_avg) then
         ierr = nf90_put_att(ncid,varid,'long_name','averaged v-momentum component')
      else
         ierr = nf90_put_att(ncid,varid,'long_name','v-momentum component')
      endif
      ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif

      end subroutine def_vars_zslice  !]
!----------------------------------------------------------------------
      subroutine wrt_zslice  ![
      ! Call wrt_zslice after completion of the time-step
      ! (After step3d_uv2)
      implicit none

      ! local
      character(len=99),save :: fname
      integer,dimension(3)   :: start
      integer                :: ncid,ierr
      integer                :: i,j,k,n
      logical,save           :: first_step=.true.

      if (first_step) then
        first_step=.false.
        call init_zslice
      endif

      if (do_avg) call calc_zslice
      if (do_avg) call calc_average

      output_time = output_time + dt

      if (output_time>=output_period) then

        if (record==nrpf) then
          call create_file('_zsl',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_zslice(ncid)
          call ncwrite(ncid,'depth',vecdep)
          ierr = nf90_close(ncid)
          record = 0
        endif

        if (.not. do_avg) call calc_zslice
        record = record+1

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr)
     &      call handle_ierr(ierr,'wrt_zslice, opening: ',fname)
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))
        if (wrt_T_zsl) then
          do n=1,nt_z
          if (do_avg) then
          call ncwrite(ncid,trim(t_vname(trc2zsc(n))),Tz_avg(1:nx,1:ny,:,n),
     &        (/bfx,bfy,1,record/))
          else
          call ncwrite(ncid,trim(t_vname(trc2zsc(n))),Tz(1:nx,1:ny,:,n),(/bfx,bfy,1,record/))
          endif
          enddo
        endif
        if (wrt_U_zsl) then
          if (do_avg) then
          call ncwrite(ncid,'u',Uz_avg(1:nx,1:ny,:),(/1,bfy,1,record/))
          else
          call ncwrite(ncid,'u',Uz(1:nx,1:ny,:),(/1,bfy,1,record/))
          endif
        endif
        if (wrt_V_zsl) then
          if (do_avg) then
          call ncwrite(ncid,'v',Vz_avg(1:nx,1:ny,:),(/bfx,1,1,record/))
          else
          call ncwrite(ncid,'v',Vz(1:nx,1:ny,:),(/bfx,1,1,record/))
          endif
        endif

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_zslice :: wrote zslice, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        output_time=0
        navg = 0

      endif

      end subroutine wrt_zslice !]
!----------------------------------------------------------------------

      end module zslice_output
