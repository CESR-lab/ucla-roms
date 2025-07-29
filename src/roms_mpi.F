      module roms_mpi

      ! This module provides mpi based parrallelization to roms

!------------------------------------
!      Pack data in buffers (depends on bf,2d/3d/, and nvars)
!      Use mpi to send and receive buffers
!      Unpack data in buffers (depends on bf,2d/3d/, and nvars)
 
!      We can have bf size and number of mpi-exchanged rows be
!      different (for instance, for the nhmg solver, we need only 1 row
!      of u, and v. Also, in that case, we only need to transfer to the
!      left and south

!      possible options for mpi_buffer_exchange: full 8, 4 (no corners),
!      2 (left, down only), send NE, recv SW

      use param
      use dimensions
      use mpi
      use netcdf 
      use nc_read_write

      implicit none
   
      private
 
      integer,public :: p_W, p_SW, p_S, p_SE, p_E, p_NE, p_N,p_NW

      real,dimension(:),allocatable :: sendW,sendE
      real,dimension(:),allocatable :: sendN,sendS
      real,dimension(:),allocatable :: sendNE,sendSE
      real,dimension(:),allocatable :: sendNW,sendSW
      real,dimension(:),allocatable :: recvW,recvE
      real,dimension(:),allocatable :: recvN,recvS
      real,dimension(:),allocatable :: recvNE,recvSE
      real,dimension(:),allocatable :: recvNW,recvSW
      integer :: szEW,szNS,szCr
      integer :: szW,szE,szN,szS,szNE,szNW,szSE,szSW
      logical :: init_buffer_done = .false.

      integer :: hl         ! Exchange halo can be less than the buffer size
      logical :: do_corners ! Argument to modify the some of the exchanges
      integer :: bl         ! possible argument to make the buffer size of
                            ! the arrays flexible

      ! Default values are halo=bf, do_corners = .true., and bl = bf

      interface exchange_xxx
        module procedure 
     &        exchange_2,exchange_22,exchange_222,exchange_2222
     &       ,exchange_3,exchange_33,exchange_333,exchange_3333
     &                  ,exchange_32,exchange_332
        ! etc, etc
      end interface

      public :: exchange_xxx
      public :: mpi_setup

      contains

!--------------------------------------------------------
      subroutine mpi_setup(ierr)  ![
      implicit none
      integer :: ierr,ncid
      integer,dimension(8) :: neighbors
      integer,dimension(2) :: subdompos
      integer,dimension(4) :: partition 
      character(len=128) :: fname
      logical :: new_part
      integer :: offx,offy


      if (analytical_grid) then
        nnodes = np_xi*np_eta
        new_part = .false.
      else
        ! Check in gridfile to check if mpi masking is present
        ierr=nf90_open(grdname,nf90_nowrite,ncid)
        if (ierr/=nf90_noerr) call handle_ierr(ierr)
        ierr=nf90_get_att(ncid,nf90_global,"neighbors",neighbors)
        ! new style partitioning has the neighbors global attribute
        if (ierr.eq.nf90_noerr) then
          if (mynode.eq.0) print * ,'neighbors att found',grdname
          new_part = .true.
          ierr=nf90_get_att(ncid,nf90_global,"partition",partition)
          nnodes=partition(2)
          ierr=nf90_get_att(ncid,nf90_global,"neighbors",neighbors)
          ierr=nf90_get_att(ncid,nf90_global,"subdompos",subdompos)
        else
          if (mynode.eq.0) print * ,'neighbors att not found',grdname
          new_part = .false. 
          nnodes = np_xi*np_eta
        endif
        ierr=nf90_close(ncid)
      endif

      ! Check whether the number of expected nodes 
      ! matches the mpirun or mpiexec argument

      ! Also check npx,npy somewhere

      if (nsize/=nnodes) then
        if (mynode.eq.0) write(*,'(/1x,A,I4,1x,A,I3,A/)')
     &   '### ERROR: mpi_setup: number of MPI-nodes should be',
     &               nnodes, 'instead of ', nsize, '.'
        ierr=99
        return
      endif

!     Determine position of subdomain in logical grid. inodes, jnodes 
!     will be needed to determine whether the subdomain is at a physical 
!     boundary.
!     Afterwards, determine the mpi ranks of the neighboring subdomains.
!     These are needed to designate MPI ranks for sources of
!     incoming and targets for outgoing messages. Negative numbers are invalid
!     and indicate that there is not a valid subdomain in that direction

      if (new_part) then
        inode=subdompos(1)
        jnode=subdompos(2)
        p_N  = neighbors(1)
        p_NE = neighbors(2)
        p_E  = neighbors(3)
        p_SE = neighbors(4)
        p_S  = neighbors(5)
        p_SW = neighbors(6)
        p_W  = neighbors(7)
        p_NW = neighbors(8)
      else
        jnode=mynode/npx
        inode=mynode-jnode*npx
        p_N = ijnode2rank(inode  ,jnode+1,npx,npy)
        p_NE= ijnode2rank(inode+1,jnode+1,npx,npy)
        p_E = ijnode2rank(inode+1,jnode  ,npx,npy)
        p_SE= ijnode2rank(inode+1,jnode-1,npx,npy)
        p_S = ijnode2rank(inode  ,jnode-1,npx,npy)
        p_SW= ijnode2rank(inode-1,jnode-1,npx,npy)
        p_W = ijnode2rank(inode-1,jnode  ,npx,npy)
        p_NW= ijnode2rank(inode-1,jnode+1,npx,npy)
      endif

!     Determine the sizes of subdomains. 
!     The first and last subdomain are possibly smaller to add
!     to the exact size of the domain. Half of the surplus is 
!     removed from those. If the surplus (offx,offy) is odd, the
!     last subdomain will be the smallest.

      nx = ceiling(1.0*gnx/npx)
      ny = ceiling(1.0*gny/npy)

        ! Check against original weird integer devision
        if (nx.ne.(gnx+npx-1)/npx ) 
     &       stop 'integer devision error in mpi_setup'
        if (ny.ne.(gny+npy-1)/npy ) 
     &       stop 'integer devision error in mpi_setup'


        offx = npx*nx-gnx
        if (inode.eq.0) then
          iSW_corn = 0
        else
          iSW_corn = inode*nx - offx/2
        endif

        if (inode.eq.0) then
          nx = nx - offx/2
        endif
        if (inode.eq.npx-1) then
          nx = nx -(offx+1)/2
        endif

        offy= npy*ny-gny
        if (jnode == 0) then
          jSW_corn = 0
        else
          jSW_corn = jnode*ny - offy/2
        endif
        if (jnode == 0) then
          ny = ny - offy/2
        endif
        if (jnode == NP_ETA-1) then
          ny = ny -(offy+1)/2
        endif

        iwest = 1;  ieast = nx
        jsouth = 1; jnorth = ny

        ! Set sizes for output files. The first and last partial files
        ! include the buffers (for the rho-dimension), the others are
        ! all sized: [1:nx],[1:ny]

        i0 = 1; i1 = nx
        j0 = 1; j1 = ny

        xi_rho = nx
        xi_u   = nx
        if (inode.eq.0) then
          i0 = 0
          xi_rho = nx+1
        endif
        if (inode.eq.npx-1) then
          i1 = nx+1
          xi_rho = nx+1
          xi_u   = nx+1
        endif

        eta_rho = ny
        eta_v   = ny
        if (jnode.eq.0) then
          j0=0
          eta_rho = ny+1
        endif
        if (jnode.eq.npy-1) then
          j1=ny+1
          eta_rho = ny+1
          eta_v   = ny+1
        endif

      end subroutine mpi_setup !]
!--------------------------------------------------------
      function ijnode2rank(inode,jnode,npx,npy) result(rank) ![

! Convert processor-grid indices inode,jnode into MPI rank. Check whether
! the indices are within the bounds of processor grid. Return them back
! into periodicity domain, if there is periodicity in either direction,
! but only if the number of subdomains is greater than one (hence
! periodicity must go through message passing).

      implicit none

      ! import/export
      integer,intent(in) :: inode,jnode
      integer,intent(in) :: npx,npy

      integer :: rank

      ! local
      integer :: i,j

      i=inode ; j=jnode
      if (ew_periodic) then
        if (i < 0) then
          i=i+npx
        elseif (i > npx-1) then
          i=i-np_xi
        endif
      endif

      if (ns_periodic) then
        if (j < 0) then
          j=j+np_eta
        elseif (j > np_eta-1) then
          j=j-np_eta
        endif
      endif

      if (0<=i .and. i<npx .and. 0<=j .and. j<npy) then
        rank = i + j*npx
      else
        rank = -1  ! meaning that the rank does not exist
      endif
      end  function ijnode2rank !]
!--------------------------------------------------------
      subroutine exchange_2(A,exch_bf,array_bf,do_corn) ![
      ! 2D variable, 2 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:) :: A
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      integer,        optional :: array_bf ! non-default allocated buffer for array
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      ! Provide the ability to work with arrays with non-standard 
      ! allocated buffer sizes
      bl = bf 
      if (present(array_bf)) then
        bl = array_bf
      endif

      ! Provide the ability to exchange less than the allocated buffer
      hl = bl
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_2dvar(A)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_2dvar(A)

      end subroutine exchange_2 !]
!--------------------------------------------------------
      subroutine exchange_22(A,B,exch_bf,do_corn) ![
      ! 2D variable, 2 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:) :: A
      real,intent(inout),dimension(:,:) :: B
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_2dvar(A)
      call pack_2dvar(B)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_2dvar(A)
      call unpack_2dvar(B)

      end subroutine exchange_22 !]
!--------------------------------------------------------
      subroutine exchange_222(A,B,C,exch_bf,do_corn) ![
      ! 2D variable, 2 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:) :: A
      real,intent(inout),dimension(:,:) :: B
      real,intent(inout),dimension(:,:) :: C
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_2dvar(A)
      call pack_2dvar(B)
      call pack_2dvar(C)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_2dvar(A)
      call unpack_2dvar(B)
      call unpack_2dvar(C)

      end subroutine exchange_222 !]
!--------------------------------------------------------
      subroutine exchange_2222(A,B,C,D,exch_bf,do_corn) ![
      ! 2D variable, 2 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:) :: A
      real,intent(inout),dimension(:,:) :: B
      real,intent(inout),dimension(:,:) :: C
      real,intent(inout),dimension(:,:) :: D
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_2dvar(A)
      call pack_2dvar(B)
      call pack_2dvar(C)
      call pack_2dvar(D)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_2dvar(A)
      call unpack_2dvar(B)
      call unpack_2dvar(C)
      call unpack_2dvar(D)

      end subroutine exchange_2222 !]
!--------------------------------------------------------
      subroutine exchange_3(A,exch_bf,do_corn) ![
      ! 3D variable, 1 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_3dvar(A)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_3dvar(A) 

      end subroutine exchange_3 !]
!--------------------------------------------------------
      subroutine exchange_33(A,B,exch_bf,do_corn) ![
      ! 3D variable, 1 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A
      real,intent(inout),dimension(:,:,:) :: B
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_3dvar(A)
      call pack_3dvar(B)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_3dvar(A)
      call unpack_3dvar(B)

      end subroutine exchange_33 !]
!--------------------------------------------------------
      subroutine exchange_333(A,B,C,exch_bf,do_corn) ![
      ! 3D variable, 1 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A
      real,intent(inout),dimension(:,:,:) :: B
      real,intent(inout),dimension(:,:,:) :: C
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_3dvar(A)
      call pack_3dvar(B)
      call pack_3dvar(C)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_3dvar(A)
      call unpack_3dvar(B)
      call unpack_3dvar(C)

      end subroutine exchange_333 !]
!--------------------------------------------------------
      subroutine exchange_3333(A,B,C,D,exch_bf,do_corn) ![
      ! 3D variable, 1 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A
      real,intent(inout),dimension(:,:,:) :: B
      real,intent(inout),dimension(:,:,:) :: C
      real,intent(inout),dimension(:,:,:) :: D
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_3dvar(A)
      call pack_3dvar(B)
      call pack_3dvar(C)
      call pack_3dvar(D)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_3dvar(A)
      call unpack_3dvar(B)
      call unpack_3dvar(C)
      call unpack_3dvar(D)

      end subroutine exchange_3333 !]
!--------------------------------------------------------
      subroutine exchange_32(A,B,exch_bf,do_corn) ![
      ! 3D variable, 1 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A
      real,intent(inout),dimension(:,:)   :: B
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_3dvar(A)
      call pack_2dvar(B)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_3dvar(A)
      call unpack_2dvar(B)

      end subroutine exchange_32 !]
!--------------------------------------------------------
      subroutine exchange_332(A,B,C,exch_bf,do_corn) ![
      ! 3D variable, 1 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A
      real,intent(inout),dimension(:,:,:) :: B
      real,intent(inout),dimension(:,:)   :: C
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      if (.not.init_buffer_done) call init_mpi_buffers

      bl = bf
      hl = bf
      if (present(exch_bf)) hl = exch_bf

      do_corners = .true.
      if (present(do_corn)) do_corners = do_corn

      call reset_buffer_sizes
      call pack_3dvar(A)
      call pack_3dvar(B)
      call pack_2dvar(C)

      call mpi_buffer_exchange

      call reset_buffer_sizes
      call unpack_3dvar(A)
      call unpack_3dvar(B)
      call unpack_2dvar(C)

      end subroutine exchange_332 !]
!--------------------------------------------------------
      subroutine reset_buffer_sizes ![
      ! reset the current size of the mpi buffers
      implicit none

      szW = 0
      szE = 0
      szS = 0
      szN = 0
      szSW = 0
      szSE = 0
      szNW = 0
      szNE = 0

      end subroutine reset_buffer_sizes !]
!--------------------------------------------------------
      subroutine pack_2dvar(A) ![
      ! add a single 2D variable to the mpi exchange buffers
      implicit none

      !input/output
      real,intent(inout),dimension(:,:) :: A

      !local
      integer(kind=4),dimension(2) :: dims
      integer(kind=4) :: nxl,nyl

      dims = shape(A)
      nxl = dims(1) -2*bl
      nyl = dims(2) -2*bl
      call pack_buffers(A,nxl,nyl,1)
        
      end subroutine pack_2dvar !]
!--------------------------------------------------------
      subroutine pack_3dvar(A) ![
      ! add a single 2D variable to the mpi exchange buffers
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A

      !local
      integer(kind=4),dimension(3) :: dims
      integer(kind=4) :: nxl,nyl,nzl

      dims = shape(A)
      nxl = dims(1) -2*bl
      nyl = dims(2) -2*bl
      nzl = dims(3)
      call pack_buffers(A,nxl,nyl,nzl)
        
      end subroutine pack_3dvar !]
!--------------------------------------------------------
      subroutine unpack_2dvar(A) ![
      ! add a single 2D variable to the mpi exchange buffers
      implicit none

      !input/output
      real,intent(inout),dimension(:,:) :: A

      !local
      integer(kind=4),dimension(2) :: dims
      integer(kind=4) :: nxl,nyl

      dims = shape(A)
      nxl = dims(1) -2*bl
      nyl = dims(2) -2*bl
      call unpack_buffers(A,nxl,nyl,1)
        
      end subroutine unpack_2dvar !]
!--------------------------------------------------------
      subroutine unpack_3dvar(A) ![
      ! add a single 2D variable to the mpi exchange buffers
      implicit none

      !input/output
      real,intent(inout),dimension(:,:,:) :: A

      !local
      integer(kind=4),dimension(3) :: dims
      integer(kind=4) :: nxl,nyl,nzl

      dims = shape(A)
      nxl = dims(1) -2*bl
      nyl = dims(2) -2*bl
      nzl = dims(3)
      call unpack_buffers(A,nxl,nyl,nzl)
        
      end subroutine unpack_3dvar !]
!--------------------------------------------------------
      subroutine init_mpi_buffers ![
      ! Allocate space for mpi exchange buffers
      implicit none

      !local
      integer(kind=4) :: mxEW,mxNS,mxCr

      if (mynode.eq.0) print *, 'init buffers'
      init_buffer_done = .true.

      mxEW = 4*bf*(j1+1-j0)*(nz+1)  ! maximum of 4 3D arrays at the time
      mxNS = 4*bf*(i1+1-i0)*(nz+1)  ! maximum of 4 3D arrays at the time
      mxCr = 4*bf*bf*(nz+1)         ! maximum of 4 3D arrays at the time

      allocate(sendW(mxEW))
      allocate(sendE(mxEW))
      allocate(sendN(mxNS))
      allocate(sendS(mxNS))
      allocate(sendNE(mxCr))
      allocate(sendSE(mxCr))
      allocate(sendNW(mxCr))
      allocate(sendSW(mxCr))

      allocate(recvW(mxEW))
      allocate(recvE(mxEW))
      allocate(recvN(mxNS))
      allocate(recvS(mxNS))
      allocate(recvNE(mxCr))
      allocate(recvSE(mxCr))
      allocate(recvNW(mxCr))
      allocate(recvSW(mxCr))

      end subroutine init_mpi_buffers !]
!--------------------------------------------------------
      subroutine pack_buffers(A,nxl,nyl,nzl) ![
      implicit none

      !input/output
      integer,intent(in) :: nxl,nyl,nzl
      real,intent(in),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl) :: A

      !local
      integer :: il0,il1,jl0,jl1

      il0=1; il1=nxl; jl0=1; jl1=nyl

      if (inode==0)        il0=0
      if (inode==NP_XI-1)  il1=nxl+1
      if (jnode==0)        jl0=0
      if (jnode==NP_ETA-1) jl1=nyl+1

      szEW = hl*(jl1+1-jl0)*nzl
      szNS = hl*(il1+1-il0)*nzl
      szCr = hl*hl*nzl

! Append buffers with data from the array to be send
!----------------------------------------------
      if (p_w>=0) then
        sendW(szW+1:szW+szEW) = reshape(A(1:hl,jl0:jl1,:),(/szEW/))
        szW = szW+szEW
      endif

      if (p_e>=0) then
        sendE(szE+1:szE+szEW) = reshape(A(nxl+1-hl:nxl,jl0:jl1,:),(/szEW/))
        szE = szE+szEW
      endif

      if (p_s>=0) then
        sendS(szS+1:szS+szNS) = reshape(A(il0:il1,1:hl,:),(/szNS/))
        szS = szS+szNS
      endif

      if (p_n>=0) then
        sendN(szN+1:szN+szNS) = reshape(A(il0:il1,nyl+1-hl:nyl,:),(/szNS/))
        szN = szN+szNS
      endif

      if (p_sw>=0) then
        sendSW(szSW+1:szSW+szCr) = reshape(A(1:hl,1:hl,:),(/szCr/))
        szSW = szSW+szCr
      endif

      if (p_se>=0) then
        sendSE(szSE+1:szSE+szCr) = reshape(A(nxl+1-hl:nxl,1:hl,:),(/szCr/))
        szSE = szSE+szCr
      endif

      if (p_ne>=0) then
        sendNE(szNE+1:szNE+szCr) = reshape(A(nxl+1-hl:nxl,nyl+1-hl:nyl,:),(/szCr/))
        szNE = szNE+szCr
      endif

      if (p_nw>=0) then
        sendNW(szNW+1:szNW+szCr) = reshape(A(1:hl,nyl+1-hl:nyl,:),(/szCr/))
        szNW = szNW+szCr
      endif

      end subroutine pack_buffers !]
!--------------------------------------------------------
      subroutine unpack_buffers(A,nxl,nyl,nzl) ![
      implicit none

      !input/output
      integer,intent(in) :: nxl,nyl,nzl
      real,intent(inout),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl) :: A

      !local
      integer,dimension(3) :: shEW,shNS,shCr
      integer :: il0,il1,jl0,jl1

      il0=1; il1=nxl; jl0=1; jl1=nyl

      if (inode==0)        il0=0
      if (inode==NP_XI-1)  il1=nxl+1
      if (jnode==0)        jl0=0
      if (jnode==NP_ETA-1) jl1=nyl+1

      shEW = shape(A(1:hl,jl0:jl1,:))
      shNS = shape(A(il0:il1,1:hl,:))
      shCr = shape(A(1:hl,1:hl,:))
      szEW = hl*(jl1+1-jl0)*nzl
      szNS = hl*(il1+1-il0)*nzl
      szCr = hl*hl*nzl

! Unpack data from exchanged buffers into array
!----------------------------------------------
      if (p_w>=0) then
        A(1-hl:0,jl0:jl1,:) = reshape(recvW(szW+1:szW+szEW),shEW)
        szW = szW+szEW
      elseif (inode>0) then
        A(1-hl:0,jl0:jl1,:) = A(1:hl,jl0:jl1,:)
      endif
 
      if (p_e>=0) then
        A(nxl+1:nxl+hl,jl0:jl1,:) = reshape(recvE(szE+1:szE+szEW),shEW)
        szE = szE+szEW
      elseif (inode<npx-1) then
        A(nxl+1:nxl+hl,jl0:jl1,:) =  A(nxl+1-hl:nxl,jl0:jl1,:)
      endif
 
      if (p_s>=0) then
        A(il0:il1,1-hl:0,:) = reshape(recvS(szS+1:szS+szNS),shNS)
        szS = szS+szNS
      elseif (jnode>0) then
        A(il0:il1,1-hl:0,:) = A(il0:il1,1:hl,:)
      endif

      if (p_n>=0) then
        A(il0:il1,nyl+1:nyl+hl,:) = reshape(recvN(szN+1:szN+szNS),shNS)
        szN = szN+szNS
      elseif (jnode<npy-1) then
        A(il0:il1,nyl+1:nyl+hl,:) = A(il0:il1,nyl+1-hl:nyl,:)
      endif

      if (p_sw>=0.and.do_corners) then
        A(1-hl:0,1-hl:0,:) = reshape(recvSW(szSW+1:szSW+szCr),shCr)
        szSW = szSW+szCr
      elseif (inode>0.and.jnode>0) then
        A(1-hl:0,1-hl:0,:) = A(1:hl,1:hl,:)
      endif

      if (p_se>=0.and.do_corners) then
        A(nxl+1:nxl+hl,1-hl:0,:) = reshape(recvSE(szSE+1:szSE+szCr),shCr)
        szSE = szSE+szCr
      elseif (inode<npx-1.and.jnode>0) then
        A(nxl+1:nxl+hl,1-hl:0,:) = A(nxl+1-hl:nxl,1:hl,:)
      endif

      if (p_nw>=0.and.do_corners) then
        A(1-hl:0,nyl+1:nyl+hl,:) = reshape(recvNW(szNW+1:szNW+szCr),shCr)
        szNW = szNW+szCr
      elseif (inode>0.and.jnode<npy-1) then
        A(1-hl:0,nyl+1:nyl+hl,:) = A(1:hl,nyl+1-hl:nyl,:)
      endif
 
      if (p_ne>=0.and.do_corners) then
        A(nxl+1:nxl+hl,nyl+1:nyl+hl,:) = reshape(recvNE(szNE+1:szNE+szCr),shCr)
        szNE = szNE+szCr
      elseif (inode<npx-1.and.jnode<npy-1) then
        A(nxl+1:nxl+hl,nyl+1:nyl+hl,:) = A(nxl+1-hl:nxl,nyl+1-hl:nyl,:)
      endif

      end subroutine unpack_buffers !]
!--------------------------------------------------------
      subroutine mpi_buffer_exchange ![
      implicit none

      ! local
      integer(kind=4) mess_count, comm(16), req(16),
     &                                  status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      integer(kind=4) i,ierr
      logical flag

! Permutation array comm(1:16) keeps track which messages are actually
! being received -- hence comm(indx)=0  means that no messages are
! expected from the direction labelled "indx", while for active messages
! "comm" keeps index of the corresponding request handle "req".
! This is needed because later in this code array "req" is subject to
! rearrangement in order to ignore directions from which no message is
! expected, as well as to ignore requests from which messages are
      
                     ! tags for receive      for send        ! each sub-domain can receive up to 8 exchanges.
      do i=1,16      !         3  5  1        4  6  2        ! the 'tag' value indicates which side/corner
        comm(i)=0    !         8     7        7     8        ! the message came from.
      enddo          !         2  6  4        1  5  3

! Prepare to receive:  
      if (p_w>=0) then
        call MPI_Irecv (recvW, szW, MPI_DOUBLE_PRECISION,
     &          p_W, 8, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (p_e>=0) then
        call MPI_Irecv (recvE, szE, MPI_DOUBLE_PRECISION,
     &          p_E, 7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (p_s>=0) then
        call MPI_Irecv (recvS, szS, MPI_DOUBLE_PRECISION,
     &          p_S, 6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (p_n>=0) then
        call MPI_Irecv (recvN, szN, MPI_DOUBLE_PRECISION,
     &          p_N, 5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (p_sw>=0) then
        call MPI_Irecv (recvSW, szSW, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (p_ne>=0) then
        call MPI_Irecv (recvNE, szNE, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (p_se>=0) then
        call MPI_Irecv (recvSE, szSE, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (p_nw>=0) then
        call MPI_Irecv (recvNW, szNW, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif

! Send everything
!----------------------------------------------
      if (p_w>=0) then
        call MPI_Isend (sendW, szW, MPI_DOUBLE_PRECISION,
     &            p_W,7,ocean_grid_comm, req(9), ierr)
        comm(9)=9
      endif
      if (p_e>=0) then
        call MPI_Isend (sendE, szE, MPI_DOUBLE_PRECISION,
     &             p_E,8,ocean_grid_comm, req(10), ierr)
        comm(10)=10
      endif
      if (p_s>=0) then
        call MPI_Isend (sendS, szS, MPI_DOUBLE_PRECISION,
     &         p_S, 5, ocean_grid_comm, req(11), ierr)
        comm(11)=11
      endif
      if (p_n>=0) then
        call MPI_Isend (sendN, szN, MPI_DOUBLE_PRECISION,
     &         p_N, 6, ocean_grid_comm, req(12), ierr)
        comm(12)=12
      endif
      if (p_sw>=0) then
        call MPI_Isend (sendSW, szSW, MPI_DOUBLE_PRECISION,
     &           p_SW, 1, ocean_grid_comm, req(13), ierr)
        comm(13)=13
      endif
      if (p_ne>=0) then
        call MPI_Isend (sendNE, szNE, MPI_DOUBLE_PRECISION,
     &           p_NE, 2, ocean_grid_comm, req(14), ierr)
        comm(14)=14
      endif
      if (p_se>=0) then
        call MPI_Isend (sendSE, szSE, MPI_DOUBLE_PRECISION,
     &           p_SE, 3, ocean_grid_comm, req(15), ierr)
        comm(15)=15
      endif
      if (p_nw>=0) then
        call MPI_Isend (sendNW, szNW, MPI_DOUBLE_PRECISION,
     &           p_NW, 4, ocean_grid_comm, req(16), ierr)
        comm(16)=16
      endif

! Verify that everything has been succesfully transferred
!----------------------------------------------
      ! 1 each for send and receive
      mess_count=0
      do i=1,16
        if (comm(i) > 0) mess_count=mess_count+1
      enddo

!  Stay in this loop untill every message has been received
      do while (mess_count>0)
        do i=1,16
          if (comm(i) > 0) then
            call MPI_Test (req(i), flag, status, ierr)
            if (flag) then
               mess_count=mess_count-1 ; comm(i)=0
            endif
          endif
        enddo
      enddo

      end subroutine mpi_buffer_exchange !]
! ----------------------------------------------------------------------

      end module roms_mpi
 
