      module particles
      ! Define and advect particles

      ! Particles live in index space: [0:nx]x[0:ny][0:nz]

#include "cppdefs.opt"
#ifdef PARTICLES

      use dimensions    ! has dimensions, a list of forcing files, and rmask, glob_istr, etc
      use read_write
      use netcdf
      use ocean2d
      use ocean3d
      use scalars

      implicit none
      private


      ! user inputs:
      integer :: np=50                     ! local number of particles
      real    :: extra_space_fac = 1.5     ! buffer space to receive extra exchanged particles 
      real    :: exchange_facx = 0.1       ! maximum number of particles expected to transfer in NS
      real    :: exchange_facy = 0.1       ! maximum number of particles expected to transfer in EW
      real    :: exchange_facc = 0.01      ! maximum number of particles expected to transfer in corners    

      
      real    :: output_period = 400
      integer :: nrpf = 100
      real    :: ppm3 = 1e-6     ! target particles per cubic meter

      ! end user inputs

      ! We choose not to define a derived particle type that has all the
      ! information for a particles because the memory storage would be 
      ! less efficient. Looping through all the x-positions of all the
      ! particles would be slower than going through an array that is
      ! stored in contiguous memory.
      
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

      integer           :: ntag                  ! current max tag (lifetime particles)
      logical :: do_init = .true.
      logical :: part_first = .true.             ! First time step for particles
      logical :: terminal_vel= .true.            ! terminal velocities on particles
      real,dimension(:,:),allocatable,target :: part    ! particles
      integer,dimension(:),allocatable :: val,idx ! needed for sorting
      real,dimension(:),allocatable :: ptmp  ! needed for sorting
      real,dimension(:,:,:),allocatable :: Wp    ! total vertical flux
      real,dimension(:) ,allocatable :: nbp_east,nbp_west ! for reseeding near open boundaries
      real,dimension(:) ,allocatable :: nbp_north,nbp_south


      real                  :: output_time = 0
      integer               :: total_rec = 0
      logical :: full_seed = .true.            ! seed with constant density
      integer :: big_number= 100000000      ! added to particles tags 
                                            ! to get unique numbers for all subdomains

      interface interp
        module procedure  interp_2D, interp_3D
      end interface


      real,allocatable,dimension(:,:) :: sn_nw, sn_n, sn_ne, rv_nw, rv_n, rv_ne,
     &                                   sn_w, sn_e, rv_w, rv_e, sn_sw, sn_s,
     &                                   sn_se, rv_sw, rv_s, rv_SE
      integer :: is_e,is_w,is_n,is_s,is_ne,is_nw,is_sw,is_se     ! exchange count
      integer :: rn_s,  rn_w,  rn_n,  rn_e                                     ! received msg buffer size
      integer :: rn_sw, rn_se, rn_nw, rn_ne

      public do_particles

      contains
!-------------------------------------------------
      subroutine do_particles![
      implicit none

      output_time = output_time + dt

      if (part_first) 
     &  call init_particles

      ! possibly not every eulerian time-step
      ! or... with sub step (for when working offline)
      call advance_particles
      call exchange_particles

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

      if (full_seed) then ! seed for constant particle density
        h_mn = sum(h(1:nx,1:ny))/(nx*ny)
        dx_mn= 1./(sum(pm(1:nx,1:ny))/(nx*ny))
        dy_mn= 1./(sum(pn(1:nx,1:ny))/(nx*ny))
        np = h_mn*dx_mn*dy_mn*ppm3
        npmx = np*extra_space_fac

        allocate(part(npmx,npv)
        if (mynode==0) print *,'Seeding ',np,' particles'


        np = 0;
        do j=1,ny
         do i=1,nx
           hvol = h(i,j)/(pm(i,j)*pn(i,j))
           call seed_particles(hvol,z_w(i,j,:),i,j)
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

      ! Tag the particles
      do ip = 1,np 
          part(ip,itag) = ip + big_number*mynode
      enddo
      ntag = np
      
      call init_arrays_mpi_roms

      end subroutine init_particles !]
!-------------------------------------------------
      subroutine advance_particles![
      ! using AB2 time stepping
      implicit none

      ! local
      real,dimension(:),pointer:: px,py,pz
      real,dimension(:),pointer:: dpx,dpy,dpz
      real,dimension(:),pointer:: dpxm,dpym,dpzm

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

      if (part_first) then
        part_first = .false.
        dpxm = dpx
        dpym = dpy
        dpzm = dpz
       endif

       px =  px + 1.5*dpx - 0.5*dpxm
       py =  py + 1.5*dpy - 0.5*dpym
       pz =  pz + 1.5*dpz - 0.5*dpzm

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
      integer :: i,j,ip
      real    :: hvol
      real,dimension(:),pointer:: px,py

      ! Check the particle density between 0.0 and 1.0
      ! and add particles if needed

      px => part(:,ipx)
      py => part(:,ipy)

      if (obc_west.and.(.not.west_exchng)) then
        if (.not.allocated(nbp_west)) allocate(nbp_west(0:nx+1))
        nbp_west = 0
        do ip = 1,np
          if (floor(px(ip))==0) then
            nbp_west(floor(py(ip))+1) = nbp_west(floor(py(ip))+1) +1
          endif
        enddo
        do j=1,ny
          hvol = h(1,j)/(pm(1,j)*pn(1,j))
          call seed_particles(hvol,z_w(1,j,:),1,j,nbp_west(j))
        enddo
      endif

      if (obc_east.and.(.not.east_exchng)) then
        if (.not.allocated(nbp_east)) allocate(nbp_east(0:nx+1))
        nbp_east = 0
        do ip = 1,np
          if (floor(px(ip))==nx-1) then
            nbp_east(floor(py(ip))+1) = nbp_east(floor(py(ip))+1) +1
          endif
        enddo
        do j=1,ny
          hvol = h(nx,j)/(pm(nx,j)*pn(nx,j))
          call seed_particles(hvol,z_w(nx,j,:),nx,j,nbp_east(j))
        enddo
      endif

      if (obc_south.and.(.not.south_exchng)) then
        if (.not.allocated(nbp_south)) allocate(nbp_south(0:ny+1))
        nbp_south = 0
        do ip = 1,np
          if (floor(py(ip))==0) then
            nbp_south(floor(px(ip))+1) = nbp_south(floor(px(ip))+1) +1
          endif
        enddo
        do i=1,nx
          hvol = h(i,1)/(pm(i,1)*pn(i,1))
          call seed_particles(hvol,z_w(i,1,:),i,1,nbp_south(i))
        enddo
      endif

      if (obc_north.and.(.not.north_exchng)) then
        if (.not.allocated(nbp_north)) allocate(nbp_north(0:ny+1))
        nbp_north = 0
        do ip = 1,np
          if (floor(py(ip))==ny-1) then
            nbp_north(floor(px(ip))+1) = nbp_north(floor(px(ip))+1) +1
          endif
        enddo
        do i=1,nx
          hvol = h(i,ny)/(pm(i,ny)*pn(i,ny))
          call seed_particles(hvol,z_w(i,ny,:),i,ny,nbp_north(i))
        enddo
      endif

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
        call random_number(z); z = z*(zw(nz)-zw(0))
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
      subroutine wrt_particles![
      use nc_read_write
      implicit none

      ! local
      integer :: ierr,diag, prev_fill_mode,ncid
      character(len=99),save  :: fname
      integer,save            :: record
      integer,dimension(3) :: start

      if (mod(total_rec,nrpf)==0) then
        call create_particles_file(total_rec,fname)
        record = 0
      endif 
      total_rec  = total_rec +1
      record = record + 1

      ierr=nf90_open(fname,nf90_write,ncid)
      ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)  

      call ncwrite(ncid,'ocean_time',(/time/),(/record/))

      start = (/1, 1, record/)
      call ncwrite(ncid,'particles',part(1:np,1:7),start)

      ierr=nf90_close(ncid)

      if (mynode == 0) then
       write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)') 
     &       'particles :: wrote output, tdays =', tdays,
     &       'step =', iic-1, 'rec =', record, '/', total_rec
     &        MYID
      endif

      end subroutine wrt_particles !]
      
!-------------------------------------------------
      subroutine create_particles_file(total_rec,fname)![
      use roms_read_write
      use nc_read_write
      implicit none

      !input/output
      integer          ,intent(in)  :: total_rec  ! total netcdf records so far for these variables
      character(len=99),intent(out) :: fname

      ! local
      integer :: ncid,ierr,varid
      character(len=10),dimension(3) :: dimnames ! dimension names
      integer,          dimension(3) :: dimsizes

      fname=trim(output_root_name) / / '_part.nc'
      call create_filename(fname,total_rec)

      ierr=nf90_create(fname,nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'create_particles_file: ')

      varid = nccreate(ncid,'ocean_time',(/'time'/),(/0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','Time since origin')
      ierr = nf90_put_att(ncid,varid,'units','second' )

      call put_global_atts(ncid, ierr)  ! put global attributes in file

      dimnames = (/'npart','seven', 'time'/)
      dimsizes = (/ np,  7,      0/)

      varid = nccreate(ncid,'particles',dimnames,dimsizes,nf90_float)
      ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Tag/x/y/z/u/v/w')
      ierr = nf90_enddef(ncid)

      if (mynode==0) print *,'particles :: created new netCDF file: ',fname

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

      px   => part(1:np,ipx)
      py   => part(1:np,ipy)
      pz   => part(1:np,ipz)
      pu   => part(1:np,ipu)
      pv   => part(1:np,ipv)
      pw   => part(1:np,ipw)
      prx   => part(1:np,iprx)
      pry   => part(1:np,ipry)
      prz   => part(1:np,iprz)

      do k=0,nz
        Wp(:,:,k) = (We(0:nx+1,0:ny+1,k)+Wi(0:nx+1,0:ny+1,k))
     &                *pm(0:nx+1,0:ny+1)*pn(0:nx+1,0:ny+1)
      enddo
      ! first data at first rho-point, assume no gradient below that
      do ip = 1,np
        if (pz(ip)< 2*nz) then

          i = floor(px(ip)+0.5)
          j = floor(py(ip)+0.5)
          k = floor(pz(ip)+0.5)
          k = max(k,1)
          iu= floor(px(ip)+1.0)
          jv= floor(py(ip)+1.0)
          kw= floor(pz(ip))
          kw= min(kw,nz-1)

          x =  px(ip)-i +0.5
          y =  py(ip)-j +0.5
          z =  pz(ip)-k +0.5
          xu=  px(ip)-iu+1.0
          yv=  py(ip)-jv+1.0
          zw=  pz(ip)-kw

          if (px(ip)<0..or.px(ip)>nx) then
            pz(ip) = 10*nz
            cycle
!           print *,'x: ', px(ip), mynode
!           stop
          endif
          if (py(ip)<0..or.py(ip)>ny) then
            pz(ip) = 10*nz
            cycle
!           print *,'y: ', py(ip), mynode
!           stop
          endif
          if (pz(ip)<0..or.pz(ip)>nz) then
            pz(ip) = 10*nz
            cycle
!           print *,'z: ', pz(ip), mynode
!           stop
          endif

          if (k<1..or.k>nz-1) then
            pz(ip) = 10*nz
            print *,'k: ', k, mynode
            stop
          endif
          if (kw<0..or.kw>nz-1) then
            print *,'kw: ', kw, mynode
            stop
          endif

          if (pz(ip)<nz) then
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
          endif
        endif
      enddo

      end subroutine rhs_particles !]
      
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

       wt(1,1,1) = (1-z)*(1-y)*(1-x);
       wt(2,1,1) =    z *(1-y)*(1-x);
       wt(1,2,1) = (1-z)*   y *(1-x);
       wt(2,2,1) =    z *   y *(1-x);
       wt(1,1,2) = (1-z)*(1-y)*   x ;
       wt(2,1,2) =    z *(1-y)*   x ;
       wt(1,2,2) = (1-z)*   y *   x ;
       wt(2,2,2) =    z *   y *   x ;
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
      integer :: ip

      is_e = 0;is_w = 0;is_n= 0;is_s= 0
      is_ne= 0;is_nw= 0;is_sw= 0;is_se = 0

      ! remove nx,ny when sending to e,n
      do ip = 1,np
        if (part(ip,ipx)>nx) then
          if (part(ip,ipy)>ny) then
            part(ip,ipx) = part(ip,ipx) - nx
            part(ip,ipy) = part(ip,ipy) - ny
            is_ne = is_ne+1
            sn_ne(is_ne,:) = part(ip,:)
          elseif (part(ip,ipy)<0) then
            part(ip,ipy) = part(ip,ipy) - ny
            is_se = is_se+1
            sn_se(is_se,:) = part(ip,:)
          else
            part(ip,ipy) = part(ip,ipy) - ny
            is_e = is_e+1
            sn_e(is_e,:) = part(ip,:)
          endif
          part(ip,ipz) = 10*nz                   ! indicate particle no longer in domain
        elseif (part(ip,ipx)< 0.0) then
          if (part(ip,ipy)>ny) then
            part(ip,ipx) = part(ip,ipx) - nx
            is_nw = is_nw+1
            sn_nw(is_nw,:) = part(ip,:)
          elseif (part(ip,ipy)<0) then
            is_sw = is_sw+1
            sn_sw(is_sw,:) = part(ip,:)
          else
            is_w = is_w+1
            sn_w(is_w,:) = part(ip,:)
          endif
          part(ip,ipz) = 10*nz
        elseif (part(ip,ipy)>ny) then            ! corners already accounted for so only if north
          part(ip,ipx) = part(ip,ipx) - nx
          is_n = is_n+1
          sn_n(is_n,:) = part(ip,:)
          part(ip,ipz) = 10*nz
        elseif (part(ip,ipy)< 0.0) then
          is_s = is_s+1
          sn_S(is_s,:) = part(ip,:)
          part(ip,ipz) = 10*nz
        endif
      enddo
      
      if (inode<2) then
        if (is_e > 0) print *, 'EXCH E!',mynode
      endif
      if (inode>0) then
        if (is_w > 0) print *, 'EXCH W!',mynode
      endif
      if (jnode>0) then
        if (is_s > 0) print *, 'EXCH S!',mynode 
      endif
      if (jnode<1) then
        if (is_n > 0) print *, 'EXCH N!',mynode
      endif
      if (is_sw > 0) print *, 'EXCH SW!',mynode
      if (is_se > 0) print *, 'EXCH SE!',mynode
      if (is_nw > 0) print *, 'EXCH NW!',mynode
      if (is_ne > 0) print *, 'EXCH NE!',mynode

      call exchange_data

      ! add nx,ny when receiving from the e,n 
      ! We order this in sw,s,se,w,e,nw,n,ne to match the desired sorting of the particle array
      if (rn_sw>0) then         ! divide buffer size by npv to get particle count  
        rn_sw = rn_sw/npv     
        part(np+1:np+rn_sw,:) = rv_sw(1:rn_sw,:)
        np = np+rn_sw
      endif

      if (rn_s>0) then
        rn_s = rn_s/npv     
        part(np+1:np+rn_s,:) = rv_s(1:rn_s,:)
        np = np+rn_s
      endif

      if (rn_se>0) then
        rn_se = rn_se/npv     
        part(np+1:np+rn_se,:) = rv_SE(1:rn_se,:)
        part(np+1:np+rn_se,ipx) = part(np+1:np+rn_se,ipx) + nx
        np = np+rn_se
      endif

      if (rn_w>0) then
        rn_w = rn_w/npv 
        part(np+1:np+rn_w,:) = rv_w(1:rn_w,:)
        np = np+rn_w
      endif

      if (rn_e>0) then
        rn_e = rn_e/npv     
        part(np+1:np+rn_e,:) = rv_e(1:rn_e,:)
        part(np+1:np+rn_e,ipx) = part(np+1:np+rn_e,ipx) + nx
        np = np+rn_e
      endif

      if (rn_nw>0) then
        rn_nw = rn_nw/npv     
        part(np+1:np+rn_nw,:) = rv_nw(1:rn_nw,:)
        part(np+1:np+rn_nw,ipy) = part(np+1:np+rn_nw,ipy) + ny
        np = np+rn_nw
      endif

      if (rn_n>0) then
        rn_n = rn_n/npv     
        part(np+1:np+rn_n,:) = rv_n(1:rn_n,:)
        part(np+1:np+rn_n,ipy) = part(np+1:np+rn_n,ipy) + ny
        np = np+rn_n
      endif

      if (rn_ne>0) then
        rn_ne = rn_ne/npv     
        part(np+1:np+rn_ne,:) = rv_ne(1:rn_ne,:)
        part(np+1:np+rn_ne,ipx) = part(np+1:np+rn_ne,ipx) + nx
        part(np+1:np+rn_ne,ipy) = part(np+1:np+rn_ne,ipy) + ny
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
      integer,dimension(:)         ,intent(inout) :: val
      integer,dimension(:),optional,intent(inout) :: idx
      !local
      integer :: nval,i,i2

      ! remove big entries
      do i=1,np
        if (val(i)>val_mx) then
          i2=i
          exit
        endif
      enddo

      do i=i2+1,nval
        if (val(i)<val_mx) then
          val(i2) = val(i)
          if (present(idx)) idx(i2) = idx(i)
          i2 = i2+1
        endif
      enddo
      np = i2-1

      end subroutine remsort !]
!------------------------------------------
      subroutine sort_particles  ![
      implicit none

      ! local
      integer :: ip,iv,val_mx

      if (.not.allocated(val)) allocate(val(npmx))
      if (.not.allocated(idx)) allocate(idx(npmx))
      if (.not.allocated(ptmp)) allocate(ptmp(npmx))

      do ip = 1,np
        idx(ip) = ip
        val(ip) = floor(part(ip,ipx)) + floor(part(ip,ipy)-1)*nx 
     &          + floor(part(ip,ipz)-1)*ny*nx
      enddo

      val_mx = 2*nx*ny*nz
      call remsort(val,idx,val_mx)

      call qsort(val,idx,1,np)

      do iv = 1,npv                     ! loop over components of particle
        ptmp = part(1:np,iv)            
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
      integer :: ssize, wsize, nsize, esize, swsize, sesize, nwsize, nesize    ! size of buffer to send

      ssize=is_s*npv; wsize=is_w*npv      ! size of buffers to send
      nsize=is_n*npv; esize=is_e*npv      ! inline this later below so only done if necessary
      swsize=is_sw*npv;   sesize=is_se*npv
      nwsize=is_nw*npv;   nesize=is_ne*npv

      ! SEND EVERYTHING:
                     ! tags for receive      for send
      do i=1,16      !         3  5  1        4  6  2
        comm(i)=0    !         8     7        7     8
      enddo          !         2  6  4        1  5  3
      
      ! must send and receive even if there are no particles to send, 
      ! as the receiver does not know this.

      if (west_msg_exch) then
        call MPI_Irecv (rv_w, rn_w, MPI_DOUBLE_PRECISION,   ! rv_w is array to store incoming data
     &          p_W, 8, ocean_grid_comm, req(1), ierr)       ! rn_w is amount of data to receive of (up to 4 arrays - compute_message_bounds.h)
        comm(1)=1                                            ! p_W   is mpi_rank of sub-domain to the west
      endif                                                  ! jtg   is used to define tile location. No tiling jtg=0
                                                             ! ocean_grid_comm = MPI_COMM_WORLD is the group of all mpi procs in ROMS run
      if (east_msg_exch) then                                ! req() stores the mpi request handle to later use
        call MPI_Irecv (rv_e, rn_e, MPI_DOUBLE_PRECISION,  !       mpi_test to see if it has been completed.
     &          p_E, 7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif

      if (south_msg_exch) then
        call MPI_Irecv (rv_s, rn_s, MPI_DOUBLE_PRECISION,
     &          p_S, 6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif

      if (north_msg_exch) then
        call MPI_Irecv (rv_n, rn_n, MPI_DOUBLE_PRECISION,
     &          p_N, 5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif

! corners ...

      if (west_msg_exch .and. south_msg_exch) then
        call MPI_Irecv (rv_sw, rn_sw, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif

      if (east_msg_exch .and. north_msg_exch) then
        call MPI_Irecv (rv_ne, rn_ne, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif

      if (east_msg_exch .and. south_msg_exch) then
        call MPI_Irecv (rv_SE, rn_se, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif

      if (west_msg_exch .and. north_msg_exch) then
        call MPI_Irecv (rv_nw, rn_nw, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif



! Send everything, sides...
!----- --------------------

!      do ipass=0,1
!        if (mod(inode+ipass,2)==0) then
          if (west_msg_exch) then
            call MPI_Isend (sn_w, wsize, MPI_DOUBLE_PRECISION,
     &              p_W, 7, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
!        else
          if (east_msg_exch) then
            call MPI_Isend (sn_e, esize, MPI_DOUBLE_PRECISION,
     &             p_E, 8, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
!        endif
!      enddo

!      do ipass=0,1
!        if (mod(jnode+ipass,2)==0) then
          if (south_msg_exch) then
            call MPI_Isend (sn_s, ssize, MPI_DOUBLE_PRECISION,
     &             p_S, 5, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
!        else
          if (north_msg_exch) then
            call MPI_Isend (sn_n, nsize, MPI_DOUBLE_PRECISION,
     &             p_N, 6, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
!        endif
!      enddo




! ...and corners:

        if (west_msg_exch .and. south_msg_exch) then
          call MPI_Isend (sn_sw, swsize, MPI_DOUBLE_PRECISION,
     &             p_SW, 1, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif

        if (east_msg_exch .and. north_msg_exch) then
          call MPI_Isend (sn_ne, nesize, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif

        if (east_msg_exch .and. south_msg_exch) then
          call MPI_Isend (sn_se, sesize, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif

        if (west_msg_exch .and. north_msg_exch) then
          call MPI_Isend (sn_nw, nwsize, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif

! Test messages for completion and unpack them in whenever they
! are ready.


      mess_count=0
      do i=1,16
        if (comm(i) > 0) mess_count=mess_count+1
      enddo

      do while (mess_count>0)
        if (comm(1)>0) then
          call MPI_Test (req(1), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(1)=0
          endif
        endif

        if (comm(2)>0) then
          call MPI_Test (req(2), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(2)=0
          endif
        endif

        if (comm(3)>0) then
          call MPI_Test (req(3), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(3)=0
          endif
        endif

        if (comm(4)>0) then
          call MPI_Test (req(4), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(4)=0
          endif
        endif

        if (comm(5)>0) then
          call MPI_Test (req(5), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(5)=0
          endif
        endif

        if (comm(6)>0) then
          call MPI_Test (req(6), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(6)=0
          endif
        endif

        if (comm(7)>0) then
          call MPI_Test (req(7), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(7)=0
          endif
        endif

        if (comm(8)>0) then
          call MPI_Test (req(8), flag, status, ierr)
          if (flag) then
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

      !local
      integer :: szx,szy,szc

      ! assume max particle transfer is proportional to length of boundary:

      szx = exchange_facx*np  ! N/S border   - as a % of number of particles
      szy = exchange_facy*np  ! E/W border
      szc = exchange_facc*np  ! corners
      

      allocate( sn_nw(szc,npv),   sn_n(szx,npv),   sn_ne(szc,npv),   ! size on outside for
     &          rv_nw(szc,npv),   rv_n(szx,npv),   rv_ne(szc,npv),   ! contiguous memory
     &          sn_w(szy,npv),                        sn_e(szy,npv),
     &          rv_w(szy,npv),                        rv_e(szy,npv),
     &          sn_sw(szc,npv),   sn_s(szx,npv),   sn_se(szc,npv),
     &          rv_sw(szc,npv),   rv_s(szx,npv),   rv_SE(szc,npv)  )

      sn_nw=init ;   sn_n=init ;   sn_ne=init
      rv_nw=init ;   rv_n=init ;   rv_ne=init
      sn_w=init ;                  sn_e=init
      rv_w=init ;                  rv_e=init
      sn_sw=init ;   sn_s=init ;   sn_se=init
      rv_sw=init ;   rv_s=init ;   rv_SE=init
      
      print *, 'INIT P MPI BUFFERS'

      end subroutine init_arrays_mpi_roms  !]

! ---------------------------------------------------------------------
      subroutine test_dummy_exchange  ![

      implicit none

      ! inputs
!      integer :: istr,iend,jstr,jend

      ! local
      real :: wipe = -8.0

      ! set the send values to something to confirm exchange works well

      sn_nw=mynode ;   sn_n=mynode ;   sn_ne=mynode  ! slight error now at beginning as these vars not init yet.
      sn_w=mynode ;                    sn_e=mynode
      sn_sw=mynode ;   sn_s=mynode ;   sn_se=mynode

      call exchange_data

      ! confirm correct receive
      if (mynode==iic-1) then
        print *, 'mynode=',mynode
        print *, 'rv_nw(4,3)=',rv_nw(4,3)
        print *, 'rv_n(1,1)=',rv_n(1,1)
        print *, 'rv_ne(1,1)=',rv_ne(1,1)
        print *, 'rv_w(1,1)=',rv_w(1,1)
        print *, 'rv_e(1,1)=',rv_e(1,1)
        print *, 'rv_sw(1,1)=',rv_sw(1,1)
        print *, 'rv_s(1,1)=',rv_s(1,1)
        print *, 'rv_SE(1,1)=',rv_SE(1,1)
      endif

      ! reset receive arrays to confirm works each timestep

      rv_nw=wipe ;   rv_n=wipe ;   rv_ne=wipe
      rv_w=wipe ;                  rv_e=wipe
      rv_sw=wipe ;   rv_s=wipe ;   rv_SE=wipe

      end subroutine test_dummy_exchange !]

! ---------------------------------------------------------------------
      recursive subroutine qsort(a, d, first, last)  ![
      implicit none

      ! inputs
      real ::  a(:), d(:)

      ! local
      real :: x, t, s
      integer :: first, last
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
         t = d(i);  d(i) = d(j);  d(j) = t
         i=i+1
         j=j-1
      end do
      if (first < i-1) call qsort(a, d, first, i-1)
      if (j+1 < last)  call qsort(a, d, j+1, last)
      end subroutine qsort  !]

!------------------------------------------
      subroutine zlev2kidx(zlev,kidx,zw)![
      implicit none
      ! import/export
      real,dimension(:),intent(inout):: zlev
      real,dimension(:),intent(in)   :: zw
      real,dimension(:),intent(out)  :: kidx
      ! local
      integer :: i,k,nt

      call qsort(zlev)
      nt = size(zlev)
      k = 0;
      do i = 1,nt
        do while (zw(k+1)<zlev(i))
          k = k+1;
        enddo
        kidx(i) = (zlev(i)-zw(k))/(zw(k+1)-zw(k));
      enddo

      end subroutine zlev2kidx !]

! ---------------------------------------------------------------------
#endif /* PARTICLES */
      end module particles
