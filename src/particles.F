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
      use ocean_vars
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
      real,dimension(:,:),allocatable :: part    ! particle tag (1), position (2-4), velocity (5:7)
      real,dimension(:,:),allocatable :: val,idx ! needed for sorting

      contains

!------------------------------------------
      subroutine init_particles![
      use mpi
      implicit none
      !local
      integer :: ip,i,j,lnp
      real    :: ppm

      ! we'll need z to index translation if we want to 
      ! have a distribution homogenous in space

      ppm = 0.1 !! particle per vertical meter
      np = 0;
      do j=1,ny
       do i=1,nx
         np = np +  round(h(i,j)*ppm)
       enddo
      enddo

      allocate(part(np,npv))
      ip = 1;
      do j=1,ny
       do i=1,nx
         lnp = round(h(i,j)*ppm)
         call random_number(z(1:lnp))*h(i,j)
         call zlev2kidx(z(1:lnp),part(ip:ip+lnp-1,ipz),z_w(i,j,:))
         call random_number(part(ip:ip+lnp-1,ipx))
         call random_number(part(ip:ip+lnp-1,ipy))
         part(ip:ip+lnp-1,ipx) = part(ip:ip+lnp-1,ipx) + i-0.5
         part(ip:ip+lnp-1,ipy) = part(ip:ip+lnp-1,ipy) + j-0.5
         ip = ip + lnp
       enddo
      enddo


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
      part(:,2) = part(:,2)*nx + 0.5
      part(:,3) = part(:,3)*ny + 0.5
      part(:,4) = part(:,4)*nz + 0.5


      end subroutine init_particles !]
!------------------------------------------
      subroutine particles![
      implicit none

      if (part_first) 
     &  call init_particles

      ! possibly not every eulerian time-step
      ! or... with sub step (for when working offline)
      call advance_particles
      call exchange_particles

      call sort_particles

      if (output_time > output_period) 
     &  call output_particles


      end subroutine particles !]
!------------------------------------------
      subroutine advance_particles![
      ! using AB2 time stepping
      implicit none

      ! local
      real,dimension(:)pointer:: dpx,dpy,dpz
      real,dimension(:)pointer:: dpxm,dpym,dpzm

      call rhs_particles

      dpx  => part(:,iprx)
      dpy  => part(:,ipry)
      dpz  => part(:,iprz)
      dpxm => part(:,idpxm)
      dpym => part(:,idpym)
      dpzm => part(:,idpzm)

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
!------------------------------------------
      subroutine rhs_particles![
      implicit none

      ! local
      integer :: i,j,k,iu,jv,kw
      real    :: x,y,z,xu,yv,zw
      real    :: pdxi,pdyi,pdz

      ! first data at first rho-point, assume no gradient below that
      do ip = 1,np
        if ((pz(ip)< 10*nz) then
          i = floor(px(ip))
          j = floor(py(ip))
          k = max(floor(pz(ip)),1)
          iu= floor(px(ip)+0.5)
          jv= floor(py(ip)+0.5)
          kw= floor(pz(ip)-0.5)

          x =  px(ip)-i
          y =  py(ip)-j
          z =  max(pz(ip)-k,0)
          xu=  px(ip)-iu+0.5
          yv=  py(ip)-jv+0.5
          zw=  pz(ip)-kw-0.5

          call interp(pu(ip),u(iu:iu+1,j:j+1,k:k+1,nnew),xu,y,z)
          call interp(pv(ip),v(i:i+1,jv:jv+1,k:k+1,nnew),x,yv,z)
          call interp(pw(ip),w(i:i+1,j:j+1,kw:kw+1,nnew),x,y,zw)

          ! local grid distances
          call interp(pdxi,dxi(i:i+1,j:j+1),x,y)
          call interp(pdyi,dyi(i:i+1,j:j+1),x,y)
          call interp(pdz,Hz(i:i+1,j:j+1,k:k+1),x,y,z)

          prx(ip) = pu(ip)*pdxi ! del i
          pry(ip) = pv(ip)*pdyi ! del j
          if (terminal_vel) then
            prz(ip) = (pw(ip)+ pwt(ip))/pdz
          else
            prz(ip) = pw(ip)/pdz
          endif
        endif
      enddo

      end subroutine rhs_particles !]
!------------------------------------------
      subroutine exchange_particles![
      ! Deal with particles entering/leaving the domain.
      ! Deal with particles switching sub-domains
      implicit none

      ! local

      do ip = 1,np
        if (part(ip,ipx)>nx+0.5) then
          buf_east(ib_east,:) = part(ip,:)
          ib_east = ib_east+1
          part(ip,itag) = -999
        elseif (part(ip,ipx)< 0.5) then
          buf_west(ib_west,:) = part(ip,:)
          ib_west = ib_west+1
          part(ip,itag) = -999
        elseif (part(ip,ipy)>ny+0.5) then
          buf_north(ib_north,:) = part(ip,:)
          ib_north = ib_north+1
          part(ip,itag) = -999
        elseif (part(ip,ipy)< 0.5) then
          buf_south(ib_south,:) = part(ip,:)
          ib_south = ib_south+1
          part(ip,itag) = -999
      enddo


     


      end subroutine exchange_particles !]
!------------------------------------------
      subroutine sort_particles![
      implicit none

      ! local
      real,dimension(np) :: val

      val = (ip + (jp-1)*nx + (zp-1)*ny
      do ip = 1,np
        idx(ip) = ip
      enddo

      call qsort(val,idx)

      ptmp = px
      do ip = 1,np
        px(ip) = ptmp(idx(ip))
      enddo
      ptmp = py
      do ip = 1,np
        py(ip) = ptmp(idx(ip))
      enddo
      ptmp = pz
      do ip = 1,np
        pz(ip) = ptmp(idx(ip))
        part(ip,:) = ptmp(idx(ip),:))
      enddo

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

      call qsort(zlev)
      nt = size(zlev)
      k = 0;
      do i = 1:nt
        while (zw(k+1)<zlev(i))
          k = k+1;
        end while
        kidx(i) = (zlev(i)-zw(k))/(zw(k+1)-zw(k));
      enddo

      end subroutine zlev2kidx !]


#endif /* PARTICLES */
      end module particles
