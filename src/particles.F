      module particles
      ! Define and advect particles

      ! STEPS:
      !
      ! 0) initialize particles:
      !    - analytical
      !    - from netcdf file
      !    - 'full'
      !
      ! 1) Move particles:
      !    - use index space for positions
      !    - ipx(np),ipy(np),ipz(np)
      !
      ! 2) Exchange with other subdomains
      !    - Can be done every few times steps when using computational margin
      !    - Add/remove particles at open boundaries
      !
      ! 3) Data structure
      !    - Simple arrays
      !    - Sort to maintain locality in memory space
      !
      ! 4) Write particles when needed
      !

#include "cppdefs.opt"
#ifdef PARTICLES

      use dimensions                                   ! has dimensions, a list of forcing files, and rmask, glob_istr, etc
      use read_write
      use netcdf
      use ocean2d
      use ocean3d
      use scalars

      implicit none
      private

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

      logical :: do_init = .true.
      logical :: part_first = .true.             ! First time step for particles
      logical :: terminal_vel= .true.            ! terminal velocities on particles
      integer           :: np                    ! local number of particles
      real,dimension(:,:),pointer    :: part    ! particle tag (1), position (2-4), velocity (5:7)
      real,dimension(:,:),allocatable :: val,idx ! needed for sorting
      real,dimension(:,:,:),allocatable :: Wp 

      real                  :: output_time = 0
      integer               :: total_rec = 0
      real                  :: output_period = 200
      integer               :: nrpf = 10

      interface interp
        module procedure  interp_2D, interp_3D
      end interface

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
      real    :: ppm

      ! we'll need z to index translation if we want to 
      ! have a distribution homogenous in space

!     ppm = 0.1 !! particle per vertical meter
!     np = 0;
!     do j=1,ny
!      do i=1,nx
!        np = np +  round(h(i,j)*ppm)
!      enddo
!     enddo

!     allocate(part(np,npv))
!     ip = 1;
!     do j=1,ny
!      do i=1,nx
!        lnp = round(h(i,j)*ppm)
!        call random_number(z(1:lnp))*h(i,j)
!        call zlev2kidx(z(1:lnp),part(ip:ip+lnp-1,ipz),z_w(i,j,:))
!        call random_number(part(ip:ip+lnp-1,ipx))
!        call random_number(part(ip:ip+lnp-1,ipy))
!        part(ip:ip+lnp-1,ipx) = part(ip:ip+lnp-1,ipx) + i-0.5
!        part(ip:ip+lnp-1,ipy) = part(ip:ip+lnp-1,ipy) + j-0.5
!        ip = ip + lnp
!      enddo
!     enddo
      allocate(Wp(0:nx+1,0:ny+1,0:nz)) ! total vertical flux in m/s


      np = 100
      allocate(part(np,npv))

      !! tag the particles
      !! possibly add more than nynode*np if np is not identical per
      !! subdomain
      do ip = 1,np 
        part(ip,itag) = ip + 10*mynode*np
      enddo

      ! A domain goes from i = [0.5,nx+0.5] (rho-indices)
      ! this is 'homogenous' in index space, not in physical space
      call random_number(part(:,2))  ! x-position
      call random_number(part(:,3))  ! y-position
      call random_number(part(:,4))  ! z-position
      part(:,2) = part(:,2)*nx
      part(:,3) = part(:,3)*ny
      part(:,4) = part(:,4)*nz


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

      px   => part(:,ipx)
      py   => part(:,ipy)
      pz   => part(:,ipz)
      dpx  => part(:,iprx)
      dpy  => part(:,ipry)
      dpz  => part(:,iprz)
      dpxm => part(:,iprxm)
      dpym => part(:,iprym)
      dpzm => part(:,iprzm)

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

      px   => part(:,ipx)
      py   => part(:,ipy)
      pz   => part(:,ipz)
      pu   => part(:,ipu)
      pv   => part(:,ipv)
      pw   => part(:,ipw)
      prx   => part(:,iprx)
      pry   => part(:,ipry)
      prz   => part(:,iprz)

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
            continue
!           print *,'x: ', px(ip), mynode
!           stop
          endif
          if (py(ip)<0..or.py(ip)>ny) then
            pz(ip) = 10*nz
            continue
!           print *,'y: ', py(ip), mynode
!           stop
          endif
          if (pz(ip)<0..or.pz(ip)>nz) then
            pz(ip) = 10*nz
            continue
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
          call interp(pu(ip),u(iu:iu+1,j:j+1,k:k+1,nnew),xu,y,z)
          call interp(pv(ip),v(i:i+1,jv:jv+1,k:k+1,nnew),x,yv,z)
          call interp(pw(ip),Wp(i:i+1,j:j+1,kw:kw+1),x,y,zw)

          ! local grid distances
          call interp(pdxi,pm(i:i+1,j:j+1),x,y)
          call interp(pdyi,pn(i:i+1,j:j+1),x,y)
          call interp(pdz,Hz(i:i+1,j:j+1,k:k+1),x,y,z)

          prx(ip) = dt*pu(ip)*pdxi ! del i
          pry(ip) = dt*pv(ip)*pdyi ! del j
!         if (terminal_vel) then
!           prz(ip) = (pw(ip)+ pwt(ip))/pdz
!         else
            prz(ip) = dt*pw(ip)/pdz
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
      implicit none

      ! local

!     ib_east = 0
!     ib_west = 0
!     ib_north= 0
!     ib_south= 0
!     ib_ne   = 0
!     ib_nw   = 0
!     ib_sw   = 0
!     ib_se   = 0
!     do ip = 1,np
!       if (part(ip,ipx)>nx) then
!         if (part(ip,ipy)>ny) then
!           ib_ne = ib_ne+1
!           sn_NE(ib_east,:) = part(ip,:)
!         elseif (part(ip,ipy)<0) then
!           ib_se = ib_se+1
!           sn_SE(ib_se,:) = part(ip,:)
!         elseif
!         ib_east = ib_east+1
!         sendE(ib_east,:) = part(ip,:)
!         part(ip,ipz) = 10*nz
!         continue
!       elseif (part(ip,ipx)< 0.0) then
!         if (part(ip,ipy)>ny) then
!           ib_nw = ib_nw+1
!           sn_NW(ib_nw,:) = part(ip,:)
!         elseif (part(ip,ipy)<0) then
!           ib_sw = ib_sw+1
!           sn_SW(ib_sw,:) = part(ip,:)
!         elseif
!         sendW(ib_west,:) = part(ip,:)
!         ib_west = ib_west+1
!         part(ip,ipz) = 10*nz
!         continue
!       elseif (part(ip,ipy)>ny) then
!         sendN(ib_north,:) = part(ip,:)
!         ib_north = ib_north+1
!         part(ip,ipz) = 10*nz
!       elseif (part(ip,ipy)< 0.0) then
!         sendS(ib_south,:) = part(ip,:)
!         ib_south = ib_south+1
!         part(ip,ipz) = 10*nz
!     enddo

!     call exchange_data

!     part(np:np+rWcnt,:) = recW(1:rWcnt,:)
!     np = np+rWcnt
!     part(np:np+rEcnt,:) = recE(1:rEcnt,:)
!     np = np+rEcnt
!     part(np:np+rScnt,:) = recS(1:rScnt,:)
!     np = np+rScnt
!     part(np:np+rNcnt,:) = recN(1:rNcnt,:)
!     np = np+rNcnt

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
      subroutine sort_particles![
      implicit none

      ! local
      real,dimension(np) :: val

!     val = ip + (jp-1)*nx + (zp-1)*ny*nx
!     do ip = 1,np
!       idx(ip) = ip
!     enddo

!     call qsort(val,idx)

!     ptmp = px
!     do ip = 1,np
!       px(ip) = ptmp(idx(ip))
!     enddo
!     ptmp = py
!     do ip = 1,np
!       py(ip) = ptmp(idx(ip))
!     enddo
!     ptmp = pz
!     do ip = 1,np
!       pz(ip) = ptmp(idx(ip))
!       part(ip,:) = ptmp(idx(ip),:)
!     enddo

      end subroutine sort_particles !]
!------------------------------------------
      subroutine zlev2kidx(zlev,kidx,zw)![
      implicit none
      ! import/export
      real,dimension(:),intent(inout):: zlev
      real,dimension(:),intent(in)   :: zw
      real,dimension(:),intent(out)  :: kidx
      ! local
      integer :: i,k,nt

!     call qsort(zlev)
!     nt = size(zlev)
!     k = 0;
!     do i = 1:nt
!       while (zw(k+1)<zlev(i))
!         k = k+1;
!       end while
!       kidx(i) = (zlev(i)-zw(k))/(zw(k+1)-zw(k));
!     enddo
      kidx = 1

      end subroutine zlev2kidx !]


#endif /* PARTICLES */
      end module particles
