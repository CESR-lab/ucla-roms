      module mpi_exchanges
      ! overloaded mpi exchange subroutines

      use param
      use dimensions
      use hidden_mpi_vars
      use roms_read_write
      use mpi
      implicit none

      private

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

!     logical,dimension(:,:,:),allocatable :: maskNS,maskEW,maskCr
      
      ! optional arguments are halo size and corner exchanges
      ! Default values are halo=bf, and do_corners = .true.
      interface exchange_xxx
        module procedure exchange2D_1var,exchange2D_2var
!       module procedure exchange3D_1var,exchange3D_2var
!       module procedure exchange2D_1var_coarse
      end interface

      public :: exchange_xxx

      contains

! new imagined structure

!     pack data in buffers (depends on bf,2d/3d/, and nvars)
!     use mpi to send and receive buffers
!     uppack data in buffers (depends on bf,2d/3d/, and nvars)
!     we can have bf size and number of mpi-exchanged rows be
!     different (for instance, for the nhmg solver, we need only 1 row
!     of u, and v. Also, in that case, we only need to transfer to the
!     left and south

!     possible options for mpi_buffer_exchange: full 8, 4 (no corners),
!     2 (left, down only)

!--------------------------------------------------------
      subroutine exchange2D_1var(A,exch_bf,do_corn) ![
      ! 2D variable, 1 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(1-bf:nx+bf,1-bf:ny+bf) :: A
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      !local
      integer(kind=4) :: halo
      logical  :: do_corners

      if (.not.init_buffer_done) call init_mpi_buffers

      if (present(exch_bf)) then
        halo= min(exch_bf,bf)
      else
        halo = bf
      endif
      if (present(do_corn)) then
        do_corners = do_corn
      else
        do_corners = .true.
      endif

      call pack_buffers(A,nx,ny,1,bf,halo)

      call mpi_buffer_exchange

      call unpack_buffers(A,nx,ny,1,bf,halo,do_corners)
        
      end subroutine exchange2D_1var !]
!--------------------------------------------------------
      subroutine exchange2D_2var(A,B,exch_bf,do_corn) ![
      ! 2D variable, 2 at the time
      ! Can exchange less data than the allocated buffer space
      implicit none

      !input/output
      real,intent(inout),dimension(1-bf:nx+bf,1-bf:ny+bf) :: A
      real,intent(inout),dimension(1-bf:nx+bf,1-bf:ny+bf) :: B
      integer(kind=4),optional :: exch_bf  ! amount of data to exchange: exch_bf <= bf
      logical,optional :: do_corn

      !local
      integer(kind=4) :: halo
      logical :: do_corners

      if (.not.init_buffer_done) call init_mpi_buffers

      if (present(exch_bf)) then
        halo= min(exch_bf,bf) 
      else
        halo = bf
      endif
      if (present(do_corn)) then
        do_corners = do_corn
      else
        do_corners = .true.
      endif

      call pack_buffers(A,nx,ny,1,bf,halo,B)

      call mpi_buffer_exchange

      call unpack_buffers(A,nx,ny,1,bf,halo,do_corners,B)
        
      end subroutine exchange2D_2var !]
!--------------------------------------------------------
      subroutine init_mpi_buffers ![
      ! Allocate space for mpi exchange buffers
      implicit none

      !local
      integer(kind=4) :: mxEW,mxNS,mxCr

      if (mynode.eq.0) print *, 'init buffers'
      init_buffer_done = .true.

      mxEW = 4*bf*(j1+1-j0)*nz  ! maximum of 4 3D arrays at the time
      mxNS = 4*bf*(i1+1-i0)*nz  ! maximum of 4 3D arrays at the time
      mxCr = 4*bf*bf*nz      ! maximum of 4 3D arrays at the time

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
      subroutine pack_buffers(A,nxl,nyl,nzl,bl,hl,B,C,D) ![
      implicit none

      !input/output
      integer,intent(in) :: nxl,nyl,nzl,bl,hl
      real,intent(in),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl) :: A
      real,intent(in),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl),optional :: B
      real,intent(in),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl),optional :: C
      real,intent(in),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl),optional :: D

      szEW = hl*(j1+1-j0)
      szNS = hl*(i1+1-i0)
      szCr = hl*hl

! Pack all data in buffers to be send
!----------------------------------------------
      if (west_msg_exch) then
        sendW(1:szEW) = reshape(A(1:hl,j0:j1,:),(/szEW/))
        szW = szEW
        if ((present(B))) then
          sendW(szEW+1:2*szEW) = reshape(B(1:hl,j0:j1,:),(/szEW/))
          szW = 2*szEW
        endif
        if ((present(C))) then
          sendW(2*szEW+1:3*szEW) = reshape(C(1:hl,j0:j1,:),(/szEW/))
          szW = 3*szEW
        endif
        if ((present(D))) then
          sendW(3*szEW+1:4*szEW) = reshape(D(1:hl,j0:j1,:),(/szEW/))
          szW = 4*szEW
        endif
      endif

      if (east_msg_exch) then
        sendE(1:szEW) = reshape(A(nx+1-hl:nx,j0:j1,:),(/szEW/))
        szE = szEW
        if (present(B)) then
          sendE(szEW+1:2*szEW) = reshape(B(1:hl,j0:j1,:),(/szEW/))
          szE = 2*szEW
        endif
        if (present(C)) then
          sendE(2*szEW+1:3*szEW) = reshape(C(1:hl,j0:j1,:),(/szEW/))
          szE = 3*szEW
        endif
        if (present(D)) then
          sendE(3*szEW+1:4*szEW) = reshape(D(1:hl,j0:j1,:),(/szEW/))
          szE = 4*szEW
        endif
      endif

      if (south_msg_exch) then
        sendS(1:szNS) = reshape(A(i0:i1,1:hl,:),(/szNS/))
        szS = szNS
        if (present(B)) then
          sendS(szNS+1:2*szNS) = reshape(B(i0:i1,1:hl,:),(/szNS/))
          szS = 2*szNS
        endif
        if (present(C)) then
          sendS(2*szNS+1:3*szNS) = reshape(C(i0:i1,1:hl,:),(/szNS/))
          szS = 3*szNS
        endif
        if (present(D)) then
          sendS(3*szNS+1:4*szNS) = reshape(D(i0:i1,1:hl,:),(/szNS/))
          szS = 4*szNS
        endif
      endif

      if (north_msg_exch) then
        sendN(1:szNS) = reshape(A(i0:i1,ny+1-hl:ny,:),(/szNS/))
        szN = szNS
        if (present(B)) then
          sendE(szNS+1:2*szNS) = reshape(B(i0:i1,ny+1-hl:ny,:),(/szNS/))
          szN = 2*szNS
        endif
        if (present(C)) then
          sendE(2*szNS+1:3*szNS) = reshape(C(i0:i1,ny+1-hl:ny,:),(/szNS/))
          szN = 3*szNS
        endif
        if (present(D)) then
          sendE(3*szNS+1:4*szNS) = reshape(D(i0:i1,ny+1-hl:ny,:),(/szNS/))
          szN = 4*szNS
        endif
      endif

      if (south_msg_exch.and.west_msg_exch) then
        sendSW(1:szCr) = reshape(A(1:hl,1:hl,:),(/szCr/))
        szSW = szCr
        if (present(B)) then
          sendSW(szCr+1:2:szCr) = reshape(B(1:hl,1:hl,:),(/szCr/))
          szSW = 2*szCr
        endif
        if (present(C)) then
          sendSW(2*szCr+1:3*szCr) = reshape(D(1:hl,1:hl,:),(/szCr/))
          szSW = 3*szCr
        endif
        if (present(D)) then
          sendSW(3*szCr+1:4*szCr) = reshape(D(1:hl,1:hl,:),(/szCr/))
          szSW = 4*szCr
        endif
      endif

      if (south_msg_exch.and.east_msg_exch) then
        sendSE(1:szCr) = reshape(A(nx+1-hl:nx,1:hl,:),(/szCr/))
        szSE = szCr
        if (present(B)) then
          sendSE(szCr+1:2:szCr) = reshape(B(nx+1-hl:nx,1:hl,:),(/szCr/))
          szSE = 2*szCr
        endif
        if (present(C)) then
         sendSE(2*szCr+1:3:szCr) = reshape(C(nx+1-hl:nx,1:hl,:),(/szCr/))
          szSE = 3*szCr
        endif
        if (present(D)) then
         sendSE(3*szCr+1:4:szCr) = reshape(D(nx+1-hl:nx,1:hl,:),(/szCr/))
          szSE = 4*szCr
        endif
      endif

      if (north_msg_exch.and.east_msg_exch) then
        sendNE(1:szCr) = reshape(A(nx+1-hl:nx,ny+1-hl:ny,:),(/szCr/))
        szNE = szCr
        if (present(B)) then
          sendNE(szCr+1:2:szCr) = reshape(B(nx+1-hl:nx,ny+1-hl:ny,:),(/szCr/))
          szNE = 2*szCr
        endif
        if (present(C)) then
          sendNE(2*szCr+1:3:szCr) = reshape(C(nx+1-hl:nx,ny+1-hl:ny,:),(/szCr/))
          szNE = 3*szCr
        endif
        if (present(D)) then
          sendNE(3*szCr+1:4:szCr) = reshape(D(nx+1-hl:nx,ny+1-hl:ny,:),(/szCr/))
          szNE = 4*szCr
        endif
      endif

      if (north_msg_exch.and.west_msg_exch) then
        sendNW(1:szCr) = reshape(A(1:hl,ny+1-hl:ny,:),(/szCr/))
        szNW = szCr
        if (present(B)) then
          sendNW(szCr+1:2:szCr) = reshape(B(1:hl,ny+1-hl:ny,:),(/szCr/))
          szNW = 2*szCr
        endif
        if (present(C)) then
          sendNW(2*szCr+1:3:szCr) = reshape(C(1:hl,ny+1-hl:ny,:),(/szCr/))
          szNW = 3*szCr
        endif
        if (present(D)) then
          sendNW(3*szCr+1:4:szCr) = reshape(D(1:hl,ny+1-hl:ny,:),(/szCr/))
          szNW = 4*szCr
        endif
      endif

      end subroutine pack_buffers !]
!--------------------------------------------------------
      subroutine unpack_buffers(A,nxl,nyl,nzl,bl,hl,do_corners,B,C,D) ![
      implicit none

      !input/output
      integer,intent(in) :: nxl,nyl,nzl,bl,hl
      logical,intent(in) :: do_corners ! Corners (.true.), no Corners (.false.)
      real,intent(inout),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl) :: A
      real,intent(inout),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl),optional :: B
      real,intent(inout),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl),optional :: C
      real,intent(inout),dimension(1-bl:nxl+bl,1-bl:nyl+bl,nzl),optional :: D

      !local
      integer,dimension(3) :: shEW,shNS,shCr

      shEW = shape(A(1:hl,j0:j1,:))
      shNS = shape(A(i0:i1,1:hl,:))
      shCr = shape(A(1:hl,1:hl,:))


! Unpack all data from exchanged buffers into arrays
!----------------------------------------------
      if (west_msg_exch) then
!       A(1-hl:0,j0:j1,:) = unpack(recvW(1:szEW),maskEW,field=0.)
        A(1-hl:0,j0:j1,:) = reshape(recvW(1:szEW),shEW)
        if (present(B)) then
          B(1-hl:0,j0:j1,:)=reshape(recvW(szEW+1:2*szEW),shEW)
        endif
        if (present(C)) then
          C(1-hl:0,j0:j1,:)=reshape(recvW(2*szEW+1:3*szEW),shEW)
        endif
        if (present(D)) then
          D(1-hl:0,j0:j1,:)=reshape(recvW(3*szEW+1:4*szEW),shEW)
        endif
      endif

      if (east_msg_exch) then
        A(nx+1:nx+hl,j0:j1,:) = reshape(recvE(1:szEW),shEW)
        if (present(B)) then
          B(nx+1:nx+hl,j0:j1,:) = reshape(recvE(szEW+1:2*szEW),shEW)
        endif
        if (present(C)) then
          C(nx+1:nx+hl,j0:j1,:) = reshape(recvE(2*szEW+1:3*szEW),shEW)
        endif
        if (present(D)) then
          D(nx+1:nx+hl,j0:j1,:) = reshape(recvE(3*szEW+1:4*szEW),shEW)
        endif
      endif

      if (south_msg_exch) then
        A(i0:i1,1-hl:0,:) = reshape(recvS(1:szNS),shNS)
        if (present(B)) then
          B(i0:i1,1-hl:0,:) = reshape(recvS(szNS+1:2*szNS),shNS)
        endif
        if (present(C)) then
          C(i0:i1,1-hl:0,:) = reshape(recvS(2*szNS+1:3*szNS),shNS)
        endif
        if (present(D)) then
          D(i0:i1,1-hl:0,:) = reshape(recvS(3*szNS+1:4*szNS),shNS)
        endif
      endif

      if (north_msg_exch) then
        A(i0:i1,ny+1:ny+hl,:) = reshape(recvN(1:szNS),shNS)
        if (present(B)) then
          B(i0:i1,ny+1:ny+hl,:) = reshape(recvN(szNS+1:2*szNS),shNS)
        endif
        if (present(C)) then
          C(i0:i1,ny+1:ny+hl,:) = reshape(recvN(2*szNS+1:3*szNS),shNS)
        endif
        if (present(D)) then
          D(i0:i1,ny+1:ny+hl,:) = reshape(recvN(3*szNS+1:4*szNS),shNS)
        endif
      endif

      if (south_msg_exch.and.west_msg_exch.and.do_corners) then
        A(1-hl:0,1-hl:0,:) = reshape(recvSW(1:szSW),shCr)
        if (present(B)) then
          B(1-hl:0,1-hl:0,:) = reshape(recvSW(szSW+1:2*szSW),shCr)
        endif
        if (present(C)) then
          C(1-hl:0,1-hl:0,:) = reshape(recvSW(2*szSW+1:3*szSW),shCr)
        endif
        if (present(D)) then
          D(1-hl:0,1-hl:0,:) = reshape(recvSW(3*szSW+1:4*szSW),shCr)
        endif
      endif

      if (south_msg_exch.and.east_msg_exch.and.do_corners) then
        A(nx+1:nx+hl,1-hl:0,:) = reshape(recvSE(1:szSE),shCr)
        if (present(B)) then
          B(nx+1:nx+hl,1-hl:0,:) = reshape(recvSE(szSE+1:2*szSE),shCr)
        endif
        if (present(C)) then
          C(nx+1:nx+hl,1-hl:0,:) = reshape(recvSE(2*szSE+1:3*szSE),shCr)
        endif
        if (present(D)) then
          D(nx+1:nx+hl,1-hl:0,:) = reshape(recvSE(3*szSE+1:4*szSE),shCr)
        endif
      endif

      if (north_msg_exch.and.west_msg_exch.and.do_corners) then
        A(1-hl:0,ny+1:ny+hl,:) = reshape(recvNW(1:szNW),shCr)
        if (present(B)) then
          B(1-hl:0,ny+1:ny+hl,:) = reshape(recvNW(szNW+1:2*szNW),shCr)
        endif
        if (present(C)) then
          C(1-hl:0,ny+1:ny+hl,:) = reshape(recvNW(2*szNW+1:3*szNW),shCr)
        endif
        if (present(D)) then
          D(1-hl:0,ny+1:ny+hl,:) = reshape(recvNW(3*szNW+1:4*szNW),shCr)
        endif
      endif

      if (north_msg_exch.and.east_msg_exch.and.do_corners) then
        A(nx+1:nx+hl,ny+1:ny+hl,:) = reshape(recvNE(1:szNE),shCr)
        if (present(B)) then
          B(nx+1:nx+hl,ny+1:ny+hl,:) = reshape(recvNE(szNE+1:2*szNE),shCr)
        endif
        if (present(C)) then
          C(nx+1:nx+hl,ny+1:ny+hl,:) = reshape(recvNE(2*szNE+1:3*szNE),shCr)
        endif
        if (present(D)) then
          D(nx+1:nx+hl,ny+1:ny+hl,:) = reshape(recvNE(3*szNE+1:4*szNE),shCr)
        endif
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
      if (mynode.eq.5) then
        print *,size(recvS),szS
        print *,p_S
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)

     
      
                     ! tags for receive      for send        ! each sub-domain can receive up to 8 exchanges.
      do i=1,16      !         3  5  1        4  6  2        ! the 'tag' value indicates which side/corner
        comm(i)=0    !         8     7        7     8        ! the message came from.
      enddo          !         2  6  4        1  5  3

! Prepare to receive:  
      if (west_msg_exch) then
        call MPI_Irecv (recvW, szW, MPI_DOUBLE_PRECISION,
     &          p_W, 8, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch) then
        call MPI_Irecv (recvE, szE, MPI_DOUBLE_PRECISION,
     &          p_E, 7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch) then
        call MPI_Irecv (recvS, szS, MPI_DOUBLE_PRECISION,
     &          p_S, 6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch) then
        call MPI_Irecv (recvN, szN, MPI_DOUBLE_PRECISION,
     &          p_N, 5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.south_msg_exch) then
        call MPI_Irecv (recvSW, szSW, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.north_msg_exch) then
        call MPI_Irecv (recvNE, szNE, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.south_msg_exch) then
        call MPI_Irecv (recvSE, szSE, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.north_msg_exch) then
        call MPI_Irecv (recvNW, szNW, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif

! Send everything
!----------------------------------------------
      if (west_msg_exch) then
        call MPI_Isend (sendW, szW, MPI_DOUBLE_PRECISION,
     &            p_W,7,ocean_grid_comm, req(9), ierr)
        comm(9)=9
      endif
      if (east_msg_exch) then
        call MPI_Isend (sendE, szE, MPI_DOUBLE_PRECISION,
     &             p_E,8,ocean_grid_comm, req(10), ierr)
        comm(10)=10
      endif
      if (south_msg_exch) then
        call MPI_Isend (sendS, szS, MPI_DOUBLE_PRECISION,
     &         p_S, 5, ocean_grid_comm, req(11), ierr)
        comm(11)=11
      endif
      if (north_msg_exch) then
        call MPI_Isend (sendN, szN, MPI_DOUBLE_PRECISION,
     &         p_N, 6, ocean_grid_comm, req(12), ierr)
        comm(12)=12
      endif
      if (west_msg_exch.and.south_msg_exch) then
        call MPI_Isend (sendSW, szSW, MPI_DOUBLE_PRECISION,
     &           p_SW, 1, ocean_grid_comm, req(13), ierr)
        comm(13)=13
      endif
      if (east_msg_exch.and.north_msg_exch) then
        call MPI_Isend (sendNE, szNE, MPI_DOUBLE_PRECISION,
     &           p_NE, 2, ocean_grid_comm, req(14), ierr)
        comm(14)=14
      endif
      if (east_msg_exch.and.south_msg_exch) then
        call MPI_Isend (sendSE, szSE, MPI_DOUBLE_PRECISION,
     &           p_SE, 3, ocean_grid_comm, req(15), ierr)
        comm(15)=15
      endif
      if (west_msg_exch.and.north_msg_exch) then
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

      end module mpi_exchanges
