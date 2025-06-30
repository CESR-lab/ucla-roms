      module particles
      ! Define and advect particles

      ! Particles live in index space: [0:nx]x[0:ny][0:nz]

#include "cppdefs.opt"

      use dimensions    ! has dimensions, a list of forcing files, and rmask, glob_istr, etc
      use netcdf
      use ocean_vars
      use scalars
      use hidden_mpi_vars
      use grid
      use nc_read_write

      implicit none
      private


#include "particles.opt"

      logical :: terminal_vel= .false. ! terminal velocities on particles (not yet implemented)

      ! We choose not to define a derived particle type that has all the
      ! information for a particles because the memory storage would be
      ! less efficient. Looping through all the x-positions of all the
      ! particles would be slower than going through an array that is
      ! stored in contiguous memory.

      character(len=130) :: particle_string
      character(len=130) :: particle_num_string
      character(len=30) :: particle_conc_string
      character(len=30) :: buffer_string


      ! indices of various quantities in particle array
      integer,parameter :: itag = 1   ! identification tag
      integer,parameter :: ipx  = 2   ! x-position (index space)
      integer,parameter :: ipy  = 3   ! y-position (index space)
      integer,parameter :: ipz  = 4   ! z-position (index space)
      integer,parameter :: ipu  = 5   ! x-velocity (m/s)
      integer,parameter :: ipv  = 6   ! y-velocity (m/s)
      integer,parameter :: ipw  = 7   ! z-velocity (m/s)
      integer,parameter :: iprx = 8   ! dx-velocity (index space)
      integer,parameter :: ipry = 9   ! dy-velocity (index space)
      integer,parameter :: iprz = 10  ! dz-velocity (index space)
      integer,parameter :: iprxm= 11  ! previous dx-velocity (index space)
      integer,parameter :: iprym= 12  ! previous dy-velocity (index space)
      integer,parameter :: iprzm= 13  ! previous dz-velocity (index space)
      integer           :: npv  = 13  ! total particle vars

      integer     :: ntag                  ! current max tag (lifetime particles)
      logical :: do_init = .true.
      logical :: part_first = .true.             ! First time step for particles

      real,dimension(:,:),allocatable,target :: part    ! particles
      real   ,dimension(:),allocatable :: val ! needed for sorting
      integer,dimension(:),allocatable :: idx ! needed for sorting
      real,dimension(:),allocatable :: ptmp  ! needed for sorting
      real,dimension(:,:,:),allocatable :: Wp    ! total vertical flux
      real,dimension(:),allocatable :: zw
      real,dimension(:) ,allocatable :: nbp_east,nbp_west ! for reseeding near open boundaries
      real,dimension(:) ,allocatable :: nbp_north,nbp_south

      integer :: szx,szy,szc   !! size of mpi buffers for particle exchanges

      real                  :: output_time = 500 ! produces output at first timestep
      integer               :: total_rec = 0
      integer ::  npmx ! allocated space for particle array
      integer :: big_number= 10000000      ! added to particles tags
                                            ! to get unique numbers for all subdomains
      integer :: n_sur=0,n_bot=0

      interface interp
        module procedure  interp_2D, interp_3D
      end interface

      real,allocatable,dimension(:) :: sn_nw, sn_n, sn_ne, rv_nw, rv_n, rv_ne,
     &                                 sn_w, sn_e, rv_w, rv_e, sn_sw, sn_s,
     &                                 sn_se, rv_sw, rv_s, rv_SE
      integer :: is_s, is_w, is_n, is_e
      integer :: is_sw,is_se,is_nw,is_ne
      integer :: rn_s, rn_w, rn_n, rn_e                       ! received msg buffer size
      integer :: rn_sw,rn_se,rn_nw,rn_ne

      public do_particles,init_particles

      contains
!-------------------------------------------------
      subroutine do_particles![
      implicit none

      !local
      integer :: ierr,ip

      if (part_first) then
!       call init_particles
!       call wrt_particles
!       output_time = 0
      endif

      output_time = output_time + dt

      ! possibly not every eulerian time-step
      ! or... with sub step (for when working offline)

!      call idealized_flow

      call advance_particles

      call exchange_particles ! marks particles as out of domain

      if (full_seed) call seed_bry

      call sort_particles

      if (output_time > output_period) then
        call wrt_particles
        output_time = 0
      endif

      end subroutine do_particles !]
!-------------------------------------------------
      subroutine init_particles![
      use mpi
      implicit none
      !local
      integer :: ip,i,j,lnp
      real    :: h_mn,dx_mn,dy_mn,est_np
      real    :: hvol


      allocate(Wp(0:nx+1,0:ny+1,0:nz)) ! total vertical flux in m/s
      allocate(zw(0:nz));

      if (full_seed) then ! seed for constant particle density
        h_mn = sum(h(1:nx,1:ny))/(nx*ny)
        dx_mn= 1./(sum(pm(1:nx,1:ny))/(nx*ny))
        dy_mn= 1./(sum(pn(1:nx,1:ny))/(nx*ny))
        np = h_mn*dx_mn*dy_mn*ppm3*nx*ny
        npmx = np*extra_space_fac
        npmx = max(npmx,pmin)
        allocate(part(npmx,npv))
!       print *,'Seeding ',np,' particles ',mynode

        np = 0;
        do j=1,ny
         do i=1,nx
           zw = z_w(i,j,:)
           hvol = h(i,j)/(pm(i,j)*pn(i,j))
           call seed_particles(hvol,zw,i,j)
         enddo
        enddo
      else
        np = 100
        npmx = np*extra_space_fac
        allocate(part(npmx,npv))

        ! this is 'homogenous' in index space, not in physical space
        call random_number(part(1:np,2))  ! x-position
        call random_number(part(1:np,3))  ! y-position
        call random_number(part(1:np,4))  ! z-position
        part(1:np,2) = part(1:np,2)*nx
        part(1:np,3) = part(1:np,3)*ny
        part(1:np,4) = part(1:np,4)*nz
      endif

      ! Tag the particles and set velocities to zero
      do ip = 1,np
          part(ip,itag) = ip + big_number*mynode
          part(ip,ipu) = 0.
          part(ip,ipv) = 0.
          part(ip,ipw) = 0.
      enddo
      ntag = np

      call init_arrays_mpi_roms

      call wrt_particles
      output_time = 0

      ! Print user options (particles.opt) to netcdf attributes
      particle_opt = ''
      if (full_seed) then
        call store_string_att(particle_opt, 'full seed')
        call store_string_att(particle_opt, ', constant density')
        write(particle_conc_string, '(A, E10.2)') ', concentration(m^-3):', ppm3
        call store_string_att(particle_opt, particle_conc_string)
        write(particle_num_string, '(A, I2)') ', # of partciles:', np

        call store_string_att(particle_opt, particle_num_string)
      else
        call store_string_att(particle_opt, 'seeded randomly')
      endif
      write(buffer_string, '(A, F5.3)') ', buffer:', extra_space_fac

      call store_string_att(particle_opt, buffer_string)
      particle_opt = trim(adjustl(particle_opt))

      end subroutine init_particles !]
!-------------------------------------------------
      subroutine advance_particles![
      ! using AB2 time stepping
      implicit none
      ! we can reduce permanent particle storage
      ! For AB2, we only need px,py,pz, and pdxm, etc...
      ! For RK like stepping, only positions are permanent

      ! local
      real,dimension(:),pointer:: px,py,pz
      real,dimension(:),pointer:: dpx,dpy,dpz
      real,dimension(:),pointer:: dpxm,dpym,dpzm
      integer :: ip

!     if (mynode==5) then
!       print *,'before rhs, px: ',px(30)
!     endif
      call rhs_particles

      px   => part(1:np,ipx)
      py   => part(1:np,ipy)
      pz   => part(1:np,ipz)
      dpx  => part(1:np,iprx)
      dpy  => part(1:np,ipry)
      dpz  => part(1:np,iprz)
      dpxm => part(1:np,iprxm)
      dpym => part(1:np,iprym)
      dpzm => part(1:np,iprzm)

      if (part_first) then !!! but some particles are birthed later
        part_first = .false.
        dpxm = dpx
        dpym = dpy
        dpzm = dpz
       endif

       px =  px + 1.5*dpx - 0.5*dpxm
       py =  py + 1.5*dpy - 0.5*dpym
       pz =  pz + 1.5*dpz - 0.5*dpzm

       do ip = 1,np
         if (pz(ip)<0) then
           pz(ip) = 0.02
           n_bot = n_bot+1
         endif
         if (pz(ip)>nz) then
           pz(ip) = nz-0.02
           n_sur = n_sur+1
         endif
       enddo

       dpxm = dpx
       dpym = dpy
       dpzm = dpz

      end subroutine advance_particles !]
!-------------------------------------------------
      subroutine seed_bry![
      ! put particles in the first row next to open boundaries
      ! in order to maintain full particle coverage
      implicit none

      !local
      integer :: i,j,ip,np_old,ntag_old
      real    :: hvol
      real,dimension(:),pointer:: px,py

      ! Check the particle density between 0.0 and 1.0
      ! and add particles if needed

      px => part(:,ipx)
      py => part(:,ipy)

      np_old = np

      if (obc_west.and.(.not.west_msg_exch)) then
        if (.not.allocated(nbp_west)) allocate(nbp_west(0:ny+1))
        nbp_west = 0
        do ip = 1,np
          if (floor(px(ip))==0) then
            nbp_west(floor(py(ip))+1) = nbp_west(floor(py(ip))+1) +1
          endif
        enddo
        do j=1,ny
          zw = z_w(1,j,:)
          hvol = h(1,j)/(pm(1,j)*pn(1,j))
          call seed_particles(hvol,zw,1,j,nbp_west(j))
        enddo
      endif

      if (obc_east.and.(.not.east_msg_exch)) then
        if (.not.allocated(nbp_east)) allocate(nbp_east(0:ny+1))
        nbp_east = 0
        do ip = 1,np
          if (floor(px(ip))==nx-1) then
            nbp_east(floor(py(ip))+1) = nbp_east(floor(py(ip))+1) +1
          endif
        enddo
        do j=1,ny
          zw = z_w(nx,j,:)
          hvol = h(nx,j)/(pm(nx,j)*pn(nx,j))
          call seed_particles(hvol,zw,nx,j,nbp_east(j))
        enddo
      endif

      if (obc_south.and.(.not.south_msg_exch)) then
        if (.not.allocated(nbp_south)) allocate(nbp_south(0:nx+1))
        nbp_south = 0
        do ip = 1,np
          if (floor(py(ip))==0) then
            nbp_south(floor(px(ip))+1) = nbp_south(floor(px(ip))+1) +1
          endif
        enddo
        do i=1,nx
          zw = z_w(i,1,:)
          hvol = h(i,1)/(pm(i,1)*pn(i,1))
          call seed_particles(hvol,zw,i,1,nbp_south(i))
        enddo
      endif

      if (obc_north.and.(.not.north_msg_exch)) then
        if (.not.allocated(nbp_north)) allocate(nbp_north(0:nx+1))
        nbp_north = 0
        do ip = 1,np
          if (floor(py(ip))==ny-1) then
            nbp_north(floor(px(ip))+1) = nbp_north(floor(px(ip))+1) +1
          endif
        enddo
        do i=1,nx
          zw = z_w(i,ny,:)
          hvol = h(i,ny)/(pm(i,ny)*pn(i,ny))
          call seed_particles(hvol,zw,i,ny,nbp_north(i))
        enddo
      endif

      do ip = np_old+1,np
        part(ip,itag) = ip + ntag-np_old + big_number*mynode
        ! first time step with these particles will be off by 1.5 for AB2
        part(ip,iprxm) = 0
        part(ip,iprym) = 0
        part(ip,iprzm) = 0
      enddo
      ntag = ntag + np-np_old

      end subroutine seed_bry !]
!-------------------------------------------------
      subroutine seed_particles(hvol,zw,i,j,np_cur)![
      implicit none
      ! import/export
      real                ,intent(in) :: hvol ! volume of the water column
      real,dimension(0:nz),intent(in) :: zw   ! w-point z-level
      integer             ,intent(in) :: i,j  ! local index position
      real,optional       ,intent(in) :: np_cur ! particles already here

      ! local
      integer :: lnp
      real,dimension(:),allocatable :: z

      lnp = floor(hvol*ppm3)
      if (present(np_cur)) then
        lnp = lnp - np_cur
      endif
      if (lnp>0) then
        allocate(z(lnp))
        call random_number(z);
        z = z*(zw(nz)-zw(0)) + zw(0)
        call zlev2kidx(z,part(np+1:np+lnp,ipz),zw(:))
        call random_number(part(np+1:np+lnp,ipx))
        call random_number(part(np+1:np+lnp,ipy))
        part(np+1:np+lnp,ipx)=part(np+1:np+lnp,ipx)+i-1
        part(np+1:np+lnp,ipy)=part(np+1:np+lnp,ipy)+j-1
        np = np + lnp
        deallocate(z)
      endif

      end subroutine seed_particles !]
!-------------------------------------------------
      subroutine wrt_particles ![
      use nc_read_write
      implicit none

      ! local
      integer :: ierr,diag, prev_fill_mode,ncid
      character(len=99),save  :: fname
      integer,save            :: record=nrpf
      integer,dimension(3) :: start

      if (record==nrpf) then
        call create_particles_file(fname)
        record = 0
      endif
      record = record + 1

      ierr=nf90_open(fname,nf90_write,ncid)
      ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

      call ncwrite(ncid,'ocean_time',(/time/),(/record/))


      if (np>0) then
        start = (/1, 1, record/)
!       print *,'maxmin: ',maxval(part(1:np,1:7)),minval(part(1:np,1:7)),mynode

        ! translate to global indices (local 0-nx):
        part(1:np,ipx) = part(1:np,ipx) + (iSW_corn+0.5)   ! start of sub-domain is actually 0.5 rho
        part(1:np,ipy) = part(1:np,ipy) + (jSW_corn+0.5)   ! but locally we defined it as 0.

        call ncwrite(ncid,'particles',part(1:np,1:7),start)

        ! translate back to local indices (0-nx):
        part(1:np,ipx) = part(1:np,ipx) - (iSW_corn+0.5)
        part(1:np,ipy) = part(1:np,ipy) - (jSW_corn+0.5)

      endif

      ierr=nf90_close(ncid)

      if (mynode == 0) then
       write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &       'particles :: wrote output, tdays =', tdays,
     &       'step =', iic-1, 'rec =', record, '/', total_rec
     &        MYID
      endif
!     print *, 'n_sur,n_bot: ',n_sur,n_bot,mynode

      end subroutine wrt_particles !]
!-------------------------------------------------
      subroutine create_particles_file(fname)![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      ! local
      integer :: ncid,ierr,varid
      character(len=7),dimension(3) :: dimnames ! dimension names
      character(len=7)  :: dn_seven = 'seven'
      character(len=7)  :: dn_npart = 'npart'
      integer,          dimension(3) :: dimsizes

      call create_file('_part',fname)

      dimnames = (/dn_npart,dn_seven, dn_tm/)
      dimsizes = (/ npmx,  7,      0/)

      ierr=nf90_open(fname,nf90_write,ncid)

      varid = nccreate(ncid,'particles',dimnames,dimsizes,nf90_double)
      if (ierr/=nf90_noerr)
     &    call handle_ierr(ierr,'create_particles_file: ')
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Tag/x/y/z/u/v/w')

      ierr=nf90_put_att( ncid, nf90_global, 'big_number', big_number )

      ierr = nf90_close(ncid)

      end subroutine create_particles_file !]
!-------------------------------------------------
      subroutine rhs_particles![
      implicit none

      ! local
      integer :: i,j,k,iu,jv,kw,ip
      real    :: x,y,z,xu,yv,zw
      real    :: pdxi,pdyi,pdz
      real,dimension(:),pointer:: px,py,pz
      real,dimension(:),pointer:: pu,pv,pw
      real,dimension(:),pointer:: prx,pry,prz

      real, dimension(2,2)   :: tmp2d   ! DPD temp to remove debug warnings
      real, dimension(2,2,2) :: tmp3d

      px  => part(1:np,ipx)
      py  => part(1:np,ipy)
      pz  => part(1:np,ipz)
      pu  => part(1:np,ipu)
      pv  => part(1:np,ipv)
      pw  => part(1:np,ipw)
      prx => part(1:np,iprx)
      pry => part(1:np,ipry)
      prz => part(1:np,iprz)

!     do k=0,nz
!       Wp(:,:,k) = (We(0:nx+1,0:ny+1,k)+Wi(0:nx+1,0:ny+1,k))
!    &                *pm(0:nx+1,0:ny+1)*pn(0:nx+1,0:ny+1)
!     enddo
      ! first data at first rho-point, assume no gradient below that
      do ip = 1,np
        if (pz(ip)< 2*nz) then

          i = floor(px(ip)+0.5)
          j = floor(py(ip)+0.5)
          k = floor(pz(ip)+0.5)
          k = max(k,1)
          k = min(k,nz-1)
          iu= floor(px(ip)+1.0)
          jv= floor(py(ip)+1.0)
          kw= floor(pz(ip))

          x =  px(ip)-i +0.5
          y =  py(ip)-j +0.5
          z =  pz(ip)-k +0.5
          xu=  px(ip)-iu+1.0
          yv=  py(ip)-jv+1.0
          zw=  pz(ip)-kw


          if (px(ip)<-0.5.or.px(ip)>nx+0.5) then
            print *,'px: ', ip,px(ip),mynode
            stop
          endif
          if (py(ip)<-0.5.or.py(ip)>ny+0.5) then
            print *,'py: ', ip,py(ip),mynode
            stop
          endif
          if (pz(ip)< 0.0.or.pz(ip)>nz) then
            print *,'pz: ', ip,pz(ip),mynode
            stop
          endif

!          call interp(pu(ip),u(iu:iu+1,j:j+1,k:k+1,nnew),xu,y,z)
!          call interp(pv(ip),v(i:i+1,jv:jv+1,k:k+1,nnew),x,yv,z)
!          call interp(pw(ip),Wp(i:i+1,j:j+1,kw:kw+1),x,y,zw)
          tmp3d=u(iu:iu+1,j:j+1,k:k+1,nnew)
          call interp(pu(ip),tmp3d,xu,y,z)
          tmp3d=v(i:i+1,jv:jv+1,k:k+1,nnew)
          call interp(pv(ip),tmp3d,x,yv,z)
          tmp3d=Wp(i:i+1,j:j+1,kw:kw+1)
          call interp(pw(ip),tmp3d,x,y,zw)

          ! local grid distances
!          call interp(pdxi,pm(i:i+1,j:j+1),x,y)
!          call interp(pdyi,pn(i:i+1,j:j+1),x,y)
!          call interp(pdz,Hz(i:i+1,j:j+1,k:k+1),x,y,z)
          tmp2d=pm(i:i+1,j:j+1)
          call interp(pdxi,tmp2d,x,y)
          tmp2d=pn(i:i+1,j:j+1)
          call interp(pdyi,tmp2d,x,y)
          tmp3d=Hz(i:i+1,j:j+1,k:k+1)
          call interp(pdz,tmp3d,x,y,z)

          prx(ip) = dt*pu(ip)*pdxi
          pry(ip) = dt*pv(ip)*pdyi
          prz(ip) = dt*pw(ip)/pdz
!         if (terminal_vel) then
!           prz(ip) = prz(ip) + dt*pwt(ip)
!         endif
        else
          prx(ip) = 0.
          pry(ip) = 0.
          prz(ip) = 0.
        endif
      enddo

!     if (mynode==1) print *,'pw(409860): ',part(409860,ipw)

      end subroutine rhs_particles !]
!-------------------------------------------------
      subroutine idealized_flow ![
      ! Idealized flow for particle testing
      implicit none
!     import/export
!     local
      integer :: i,j,k
      integer :: z

      if (mynode<3) then
        do k = 1,nz
         do j = -1,ny+2
          do i = -1,nx+2
            z = z_r(i,j,k)
            u(i,j,k,:) = 0.1*(tanh(0.5*(z+30))+1)
            v(i,j,k,:) = 0.
          enddo
         enddo
        enddo
      else
        u = 0;
        v = 0;
      endif

      do k = 1,nz
       do j = -1,ny+2
        do i = -1,nx+2
         FlxU(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*u(i,j,k,1)
         FlxV(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*v(i,j,k,1)
        enddo
       enddo
      enddo
      Wp(:,:,0) = 0
      do k = 1,nz
        do j = 0,ny+1
         do i = 0,nx+1
           Wp(i,j,k) = Wp(i,j,k-1)-FlxU(i+1,j,k)+FlxU(i,j,k)
     &                            -FlxV(i,j+1,k)+FlxV(i,j,k)
         enddo
        enddo
      enddo
      do k = 1,nz
         Wp(:,:,k) = Wp(:,:,k)*pm(0:nx+1,0:ny+1)*pn(0:nx+1,0:ny+1)
      enddo

      end subroutine idealized_flow !]
!-------------------------------------------------
       subroutine interp_2D(fi,f,x,y) ![
       ! 2D Linear interpolation
       implicit none
!      import/export
       real(kind=8)               ,intent(out):: fi
       real(kind=8),dimension(2,2),intent(in) :: f
       real(kind=8)               ,intent(in) :: x,y
!      local
       real(kind=8),dimension(2,2) :: wt

       wt(1,1) = (1-x)*(1-y);
       wt(2,1) =    x *(1-y);
       wt(1,2) = (1-x)*   y ;
       wt(2,2) =    x *   y ;
!      print *,'int_lin_2d: ', sum(wt)
       fi   = sum(f*wt)
       end subroutine interp_2D !]
!-------------------------------------------------
       subroutine interp_3D(fi,f,x,y,z) ![
       !! 3D Linear interpolation
       implicit none
!      import/export
       real(kind=8)                 ,intent(out):: fi
       real(kind=8),dimension(2,2,2),intent(in) :: f
       real(kind=8)                 ,intent(in) :: x,y,z
!      local
       real(kind=8),dimension(2,2,2) :: wt

       wt(1,1,1) = (1-x)*(1-y)*(1-z);
       wt(2,1,1) =    x *(1-y)*(1-z);
       wt(1,2,1) = (1-x)*   y *(1-z);
       wt(2,2,1) =    x *   y *(1-z);
       wt(1,1,2) = (1-x)*(1-y)*   z ;
       wt(2,1,2) =    x *(1-y)*   z ;
       wt(1,2,2) = (1-x)*   y *   z ;
       wt(2,2,2) =    x *   y *   z ;
!      print *,'int_lin_3d: ', sum(wt)
       fi   = sum(f*wt)
       end subroutine interp_3D !]
!------------------------------------------
      subroutine exchange_particles![
      ! Deal with particles entering/leaving the domain.
      ! Deal with particles switching sub-domains
      use hidden_mpi_vars
      implicit none

      ! local
      integer :: i,ip,ierr

      is_e = 0;is_w = 0;is_n= 0;is_s= 0
      is_ne= 0;is_nw= 0;is_sw= 0;is_se = 0

      ! remove nx,ny when sending to e,n
      do ip = 1,np

       if (is_n>szx-13.or.is_s>szx-13) then
         print *, 'is_n/s: ',is_n,is_s,szx
         stop
       endif
       if (is_e>szy-13.or.is_w>szy-13) then
         print *, 'is_e/w: ',is_e,is_w,szy
         stop
       endif

       if (part(ip,ipz)<nz) then
        if (part(ip,ipx)>nx) then
          if (east_msg_exch) then
            if (part(ip,ipy)>ny.and.north_msg_exch) then
              part(ip,ipx) = part(ip,ipx) - nx
              part(ip,ipy) = part(ip,ipy) - ny
              sn_ne(is_ne+1:is_ne+npv) = part(ip,:)
              is_ne = is_ne+npv
            elseif (part(ip,ipy)<0.and.south_msg_exch) then
              part(ip,ipx) = part(ip,ipx) - nx
              sn_se(is_se+1:is_se+npv) = part(ip,:)
              is_se = is_se+npv
            else
              part(ip,ipx) = part(ip,ipx) - nx
              sn_e(is_e+1:is_e+npv) = part(ip,:)
              is_e = is_e+npv
            endif
          endif
          part(ip,ipz) = 10*nz                   ! indicate particle no longer in domain
        elseif (part(ip,ipx)< 0.0) then
          if (west_msg_exch) then
            if (part(ip,ipy)>ny) then
              part(ip,ipy) = part(ip,ipy) - ny
              sn_nw(is_nw+1:is_nw+npv) = part(ip,:)
              is_nw = is_nw+npv
            elseif (part(ip,ipy)<0) then
              sn_sw(is_sw+1:is_sw+npv) = part(ip,:)
              is_sw = is_sw+npv
            else
              sn_w(is_w+1:is_w+npv) = part(ip,:)
              is_w = is_w+npv
            endif
          endif
          part(ip,ipz) = 10*nz
        elseif (part(ip,ipy)>ny) then              ! corners already accounted for so only if north
          if (north_msg_exch) then
            part(ip,ipy) = part(ip,ipy) - ny
            sn_n(is_n+1:is_n+npv) = part(ip,:)
            is_n = is_n+npv
          endif
          part(ip,ipz) = 10*nz
        elseif (part(ip,ipy)< 0.0) then
          if (south_msg_exch) then
            sn_s(is_s+1:is_s+npv) = part(ip,:)
            is_s = is_s+npv
          endif
          part(ip,ipz) = 10*nz
        endif
       endif ! (part(ip,ipz)<nz)
      enddo

      call exchange_data

      if (mynode==1) then
        do ip = 1,np
          if (part(ip,itag)==200000177) then
            print *,'trck 1: ', part(ip,itag),part(ip,ipx),mynode
          endif
        enddo
      endif

      ! add nx,ny when receiving from the e,n
      ! We order this in sw,s,se,w,e,nw,n,ne to match the desired sorting of the particle array
      if (rn_sw>0) then         ! divide buffer size by npv to get particle count
        rn_sw = rn_sw/npv
        do ip = 1,rn_sw
          do i = 1,npv
            part(np+ip,i) = rv_sw(i+(ip-1)*npv)
          enddo
        enddo
        np = np+rn_sw
      endif

      if (rn_s>0) then
        rn_s = rn_s/npv
        do ip = 1,rn_s
          do i = 1,npv
            part(np+ip,i) = rv_s(i+(ip-1)*npv)
          enddo
        enddo
        np = np+rn_s
!       if (mynode==1) print *,'rn_s: ',rn_s
      endif

      if (rn_se>0) then
        rn_se = rn_se/npv
        do ip = 1,rn_se
          do i = 1,npv
            part(np+ip,i) = rv_se(i+(ip-1)*npv)
          enddo
          part(np+ip,ipx) = part(np+ip,ipx)+nx
        enddo
        np = np+rn_se
!       if (mynode==1) print *,'rn_se: ',rn_se
      endif

      if (rn_w>0) then
        rn_w = rn_w/npv
        do ip = 1,rn_w
          do i = 1,npv
            part(np+ip,i) = rv_w(i+(ip-1)*npv)
          enddo
        enddo
        np = np+rn_w
!       if (mynode==1) print *,'rn_w: ',rn_w,np
!       if (mynode==1) print *,'rv_w: ',rv_w(1:rn_w*npv)
      endif

      if (rn_e>0) then
        rn_e = rn_e/npv
        do ip = 1,rn_e
          do i = 1,npv
            part(np+ip,i) = rv_e(i+(ip-1)*npv)
          enddo
          part(np+ip,ipx) = part(np+ip,ipx)+nx
        enddo
        np = np+rn_e
!       if (mynode==1) print *,'rn_e: ',rn_e
      endif

      if (rn_nw>0) then
        rn_nw = rn_nw/npv
        do ip = 1,rn_nw
          do i = 1,npv
            part(np+ip,i) = rv_nw(i+(ip-1)*npv)
          enddo
          part(np+ip,ipy) = part(np+ip,ipy)+ny
        enddo
        np = np+rn_nw
!       if (mynode==1) print *,'rn_nw: ',rn_nw
      endif

      if (rn_n>0) then
        rn_n = rn_n/npv
        do ip = 1,rn_n
          do i = 1,npv
            part(np+ip,i) = rv_n(i+(ip-1)*npv)
          enddo
          part(np+ip,ipy) = part(np+ip,ipy)+ny
        enddo
        np = np+rn_n
!       if (mynode==1) print *,'rn_n: ',rn_n
      endif

      if (rn_ne>0) then
        rn_ne = rn_ne/npv
        do ip = 1,rn_ne
          do i = 1,npv
            part(np+ip,i) = rv_ne(i+(ip-1)*npv)
          enddo
          part(np+ip,ipx) = part(np+ip,ipx)+nx
          part(np+ip,ipy) = part(np+ip,ipy)+ny
        enddo
        np = np+rn_ne
      endif

      end subroutine exchange_particles !]
!------------------------------------------
      subroutine extend_size_particles![
      implicit none

      ! local
!     real,dimension(:,:),allocatable :: tmp_part

!     allocate(tmp_part(np,npv))

!     tmp_part = part

!     deallocate(part)
!     allocate(part(np+extra,npv))
!     part(1:np,:) = tmp_part
!     deallocate(tmp_part)

      print *,'extended particles array'

      end subroutine extend_size_particles !]
!------------------------------------------
      subroutine remsort(val,idx,val_mx) ![
      implicit none
      !import/export
      integer                      ,intent(in)    :: val_mx
      real   ,dimension(:)         ,intent(inout) :: val
      integer,dimension(:),optional,intent(inout) :: idx
      !local
      integer :: i,i2,ierr,i3

      ! remove big entries
!     print *, mynode,maxval(val(1:np)),np

      i2 = 0
      do i=1,np
        if (val(i)>val_mx) then
          i2=i
          exit
        endif
      enddo

      if (i2>0) then ! only do this if large val was found
!       i3 = i2
!       print *,mynode,i2,val(i2)
        do i=i2+1,np
          if (val(i)<val_mx) then
            val(i2) = val(i)
            if (present(idx)) idx(i2) = idx(i)
            i2 = i2+1
          endif
        enddo
        np = i2-1
!       print *,mynode,i3,val(i3)
      endif

!     print *, mynode,maxval(val(1:np)),np
!     call MPI_barrier (ocean_grid_comm, ierr)
!     stop


      end subroutine remsort !]
!------------------------------------------
      subroutine sort_particles  ![
      implicit none

      ! local
      integer :: ip,iv,val_mx,ierr,np_old

      if (.not.allocated(val)) allocate(val(npmx))
      if (.not.allocated(idx)) allocate(idx(npmx))
      if (.not.allocated(ptmp)) allocate(ptmp(npmx))

      np_old = np
      do ip = 1,np
        idx(ip) = ip
        val(ip) = floor(part(ip,ipx)) + floor(part(ip,ipy))*nx
     &          + floor(part(ip,ipz))*ny*nx
      enddo

      val_mx = 2*nx*ny*nz
      call remsort(val,idx,val_mx)

      call qsort(val,idx,1,np)

      do iv = 1,npv                     ! loop over components of particle
        ptmp(1:np_old) = part(1:np_old,iv)
        do ip = 1,np
          part(ip,iv) = ptmp(idx(ip))
        enddo
      enddo

      end subroutine sort_particles !]
! ---------------------------------------------------------------------
      subroutine exchange_data  ![
      ! orginally to exchange particles

      use hidden_mpi_vars
      use mpi

      implicit none

      ! local
      integer :: i, j, ierr, mess_count, comm(16), req(16), status(MPI_STATUS_SIZE)
!# ifdef ALT_SEND_ORD
!      integer :: ipass
!# endif
      logical :: flag

      ! SEND EVERYTHING:
                     ! tags for receive      for send
      do i=1,16      !         3  5  1        4  6  2
        comm(i)=0    !         8     7        7     8
      enddo          !         2  6  4        1  5  3

      ! must send and receive even if there are no particles to send,
      ! as the receiver does not know this.

      if (west_msg_exch) then
        call MPI_Irecv (rv_w, szy, MPI_DOUBLE_PRECISION,   ! rv_w is array to store incoming data and
     &          p_W, 8, ocean_grid_comm, req(1), ierr)     !  szy is its size
        comm(1)=1                                          ! p_W   is mpi_rank of sub-domain to the west
      endif
                                                           ! ocean_grid_comm = MPI_COMM_WORLD is the group of all mpi procs in ROMS run
      if (east_msg_exch) then                              ! req() stores the mpi request handle to later use
        call MPI_Irecv (rv_e, szy, MPI_DOUBLE_PRECISION,   ! mpi_test to see if it has been completed.
     &          p_E, 7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif

      if (south_msg_exch) then
        call MPI_Irecv (rv_s, szx, MPI_DOUBLE_PRECISION,
     &          p_S, 6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif

      if (north_msg_exch) then
        call MPI_Irecv (rv_n, szx, MPI_DOUBLE_PRECISION,
     &          p_N, 5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif

! corners ...

      if (west_msg_exch .and. south_msg_exch) then
        call MPI_Irecv (rv_sw, szc, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif

      if (east_msg_exch .and. north_msg_exch) then
        call MPI_Irecv (rv_ne, szc, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif

      if (east_msg_exch .and. south_msg_exch) then
        call MPI_Irecv (rv_se, szc, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif

      if (west_msg_exch .and. north_msg_exch) then
        call MPI_Irecv (rv_nw, szc, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif



! Send everything, sides...
!----- --------------------

!      do ipass=0,1
!        if (mod(inode+ipass,2)==0) then
          if (west_msg_exch) then
            call MPI_Isend (sn_w(1:is_w), is_w, MPI_DOUBLE_PRECISION,
     &              p_W, 7, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
!        else
          if (east_msg_exch) then
            call MPI_Isend (sn_e(1:is_e), is_e, MPI_DOUBLE_PRECISION,
     &             p_E, 8, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
!        endif
!      enddo

!      do ipass=0,1
!        if (mod(jnode+ipass,2)==0) then
          if (south_msg_exch) then
            call MPI_Isend (sn_s(1:is_s), is_s, MPI_DOUBLE_PRECISION,
     &             p_S, 5, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
!        else
          if (north_msg_exch) then
            call MPI_Isend (sn_n(1:is_n), is_n, MPI_DOUBLE_PRECISION,
     &             p_N, 6, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
!        endif
!      enddo




! ...and corners:

        if (west_msg_exch .and. south_msg_exch) then
          call MPI_Isend (sn_sw(1:is_sw), is_sw, MPI_DOUBLE_PRECISION,
     &             p_SW, 1, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif

        if (east_msg_exch .and. north_msg_exch) then
          call MPI_Isend (sn_ne(1:is_ne), is_ne, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif

        if (east_msg_exch .and. south_msg_exch) then
          call MPI_Isend (sn_se(1:is_se), is_se, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif

        if (west_msg_exch .and. north_msg_exch) then
          call MPI_Isend (sn_nw(1:is_nw), is_nw, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif

! Test messages for completion and unpack them in whenever they
! are ready.


      mess_count=0
      do i=1,16
        if (comm(i) > 0) mess_count=mess_count+1
      enddo

!     if (mess_count>0) print *, mynode, 'mess_count=',mess_count

      do while (mess_count>0)
        if (comm(1)>0) then
          call MPI_Test (req(1), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_w,ierr)
            mess_count=mess_count-1 ; comm(1)=0
          endif
        endif

        if (comm(2)>0) then
          call MPI_Test (req(2), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_e,ierr)
            mess_count=mess_count-1 ; comm(2)=0
          endif
        endif

        if (comm(3)>0) then
          call MPI_Test (req(3), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_s,ierr)
            mess_count=mess_count-1 ; comm(3)=0
          endif
        endif

        if (comm(4)>0) then
          call MPI_Test (req(4), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_n,ierr)
            mess_count=mess_count-1 ; comm(4)=0
          endif
        endif

        if (comm(5)>0) then
          call MPI_Test (req(5), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_sw,ierr)
            mess_count=mess_count-1 ; comm(5)=0
          endif
        endif

        if (comm(6)>0) then
          call MPI_Test (req(6), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_ne,ierr)
            mess_count=mess_count-1 ; comm(6)=0
          endif
        endif

        if (comm(7)>0) then
          call MPI_Test (req(7), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_se,ierr)
            mess_count=mess_count-1 ; comm(7)=0
          endif
        endif

        if (comm(8)>0) then
          call MPI_Test (req(8), flag, status, ierr)
          if (flag) then
            call MPI_get_count(status,MPI_DOUBLE_PRECISION,rn_nw,ierr)
            mess_count=mess_count-1 ; comm(8)=0
          endif
        endif

        do j=9,16
          if (comm(j) > 0) then
            call MPI_Test (req(j), flag, status, ierr)         ! make sure all sent messages have been received elsewhere
            if (flag) then                                     ! before moving on
               mess_count=mess_count-1 ; comm(j)=0
            endif
          endif
        enddo

      enddo !<--  while (mess_count>0)


      end subroutine exchange_data  !]
! ---------------------------------------------------------------------
      subroutine init_arrays_mpi_roms  ![
      implicit none

      ! Assume max particle transfer is proportional to length of boundary
      ! and number of particles in the sub domain

      szx = exchange_facx*npmx*npv  ! N/S border   - as a % of number of particles
      szy = exchange_facy*npmx*npv  ! E/W border
      szc = exchange_facc*npmx*npv  ! corners


      allocate(sn_nw(szc),sn_n(szx),sn_ne(szc),
     &         rv_nw(szc),rv_n(szx),rv_ne(szc),
     &         sn_sw(szc),sn_s(szx),sn_se(szc),
     &         rv_sw(szc),rv_s(szx),rv_se(szc),
     &         sn_w(szy),sn_e(szy),
     &         rv_w(szy),rv_e(szy)    )

      end subroutine init_arrays_mpi_roms  !]
! ---------------------------------------------------------------------
      recursive subroutine qsort(a, d, first, last)  ![
      implicit none

      ! inputs
      real   ,dimension(:)         ,intent(inout) :: a
      integer,dimension(:),optional,intent(inout) :: d
      integer :: first, last

      ! local
      real :: x, t, s
      integer :: i, j

      x = a( (first+last) / 2 )
      i = first
      j = last
      do
         do while (a(i) < x)
            i=i+1
         end do
         do while (x < a(j))
            j=j-1
         end do
         if (i >= j) exit
         t = a(i);  a(i) = a(j);  a(j) = t
         if (present(d)) then
           t = d(i);  d(i) = d(j);  d(j) = t
         endif
         i=i+1
         j=j-1
      end do
      if (present(d)) then
        if (first < i-1) call qsort(a, d, first, i-1)
        if (j+1 < last)  call qsort(a, d, j+1, last)
      else
        if (first < i-1) call qsort(a,first=first,last=i-1)
        if (j+1 < last)  call qsort(a,first=j+1,last=last)
      endif

      end subroutine qsort  !]
!------------------------------------------
      subroutine zlev2kidx(zlev,kidx,zw)![
      implicit none
      ! import/export
      real,dimension(:)   ,intent(inout):: zlev
      real,dimension(0:nz),intent(in)   :: zw
      real,dimension(:)   ,intent(out)  :: kidx
      ! local
      integer :: i,k,npz

      npz = size(zlev)

      call qsort(zlev,first=1,last=npz)

      k = 0;
      do i = 1,npz
        do while (zw(k+1)<zlev(i))
          k = k+1;
        enddo
        kidx(i) = (zlev(i)-zw(k))/(zw(k+1)-zw(k))+k;
      enddo

      end subroutine zlev2kidx !]

      end module particles
