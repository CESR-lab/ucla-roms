module mg_mpi_exchange

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids

  implicit none

  interface fill_halo
     module procedure   &
          fill_halo_2D, &
          fill_halo_3D, &
          fill_halo_3D_relax, &
          fill_halo_4D
  end interface fill_halo

  interface set_rurvbc2zero
     module procedure   &
          set_rurvbc2zero_3D
  end interface set_rurvbc2zero

contains

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_2D(lev,a2D,lbc_null)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(:,:), pointer, intent(inout)::a2D
    character(len=1), optional, intent(in) :: lbc_null

    integer(kind=ip) :: nx, ny
    integer(kind=ip) :: nh
    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    logical :: lbc
    character(len=1) :: cuv

    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    logical :: flag_sw_s, flag_sw_w
    logical :: flag_se_s, flag_se_e
    logical :: flag_ne_n, flag_ne_e
    logical :: flag_nw_n, flag_nw_w

    real(kind=rp), dimension(:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    if (present(lbc_null)) then
       cuv=trim(lbc_null)
       lbc=.true.
    else
       lbc=.false.
    endif

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nh = 1

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    !- 
    !-

    sendS => gbuffers(lev)%sendS2D1
    recvS => gbuffers(lev)%recvS2D1
    sendN => gbuffers(lev)%sendN2D1
    recvN => gbuffers(lev)%recvN2D1

    sendE => gbuffers(lev)%sendE2D1
    recvE => gbuffers(lev)%recvE2D1
    sendW => gbuffers(lev)%sendW2D1
    recvW => gbuffers(lev)%recvW2D1

    sendSW => gbuffers(lev)%sendSW2D1
    sendSE => gbuffers(lev)%sendSE2D1
    sendNW => gbuffers(lev)%sendNW2D1
    sendNE => gbuffers(lev)%sendNE2D1

    recvSW => gbuffers(lev)%recvSW2D1
    recvSE => gbuffers(lev)%recvSE2D1
    recvNW => gbuffers(lev)%recvNW2D1
    recvNE => gbuffers(lev)%recvNE2D1

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103
    swnetag = 104
    senwtag = 105
    neswtag = 106
    nwsetag = 107

    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvS,nx*nh,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    elseif ((lbc).and.(trim(cuv)=='v')) then
!       a2d(1,:) = zero
    else !!Homogenous Neumann  
       a2D(1-nh:0,1:nx) = a2D(nh:1:-1,1:nx)
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvE,ny*nh,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    elseif ((lbc).and.(trim(cuv)=='u')) then
!       a2d(:,nx+1) = zero
    else !!Homogenous Neumann 
       a2D(1:ny,nx+1) = a2D(1:ny,nx)
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvN,nx*nh,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    elseif ((lbc).and.(trim(cuv)=='v')) then 
!       a2d(ny+1,:) = zero
    else !!Homogenous Neumann  
       a2D(ny+1:ny+nh,1:nx) = a2D(ny:ny-nh+1:-1,1:nx)
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvW,ny*nh,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    elseif ((lbc).and.(trim(cuv)=='u')) then
!       a2d(:,1) = zero
    else !!Homogenous Neumann
       a2D(1:ny,1-nh:0) = a2D(1:ny,nh:1:-1)
    endif

    flag_sw_s = .false.
    flag_sw_w = .false.
    if (southwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvSW,nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            neswtag,MPI_COMM_WORLD,req(5),ierr)
       comm(5)=5
    elseif ((lbc).and.(west.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')) then
!       a2d(1-nh:0,1-nh:0) = zero
    elseif (south.ne.MPI_PROC_NULL) then
       flag_sw_s = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_sw_w = .true.
    else !!Homogenous Neumann
       a2D(1-nh:0,1-nh:0) = a2D(nh:1:-1,nh:1:-1)
    endif

    flag_se_s = .false.
    flag_se_e = .false.
    if (southeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvSE,nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            nwsetag,MPI_COMM_WORLD,req(6),ierr)
       comm(6)=6
    elseif ((lbc).and.(east.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')) then
!       a2d(1-nh:0,nx+1:nx+nh) = zero
    elseif (south.ne.MPI_PROC_NULL) then
       flag_se_s = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_se_e = .true.
    else !!Homogenous Neumann
       a2D(1-nh:0,nx+1:nx+nh) = a2D(nh:1:-1,nx:nx-nh+1:-1)
    endif

    flag_ne_n = .false.
    flag_ne_e = .false.
    if (northeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvNE,nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            swnetag,MPI_COMM_WORLD,req(7),ierr)
       comm(7)=7
    elseif (((lbc).and.(east.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')).or. &
         ((lbc).and.(trim(cuv)=='v')) ) then
!       a2d(ny+1:ny+nh,nx+1:nx+nh) = zero
    elseif (north.ne.MPI_PROC_NULL) then
       flag_ne_n = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_ne_e = .true.
    else !!Homogenous Neumann
       a2D(ny+1:ny+nh,nx+1:nx+nh) = a2D(ny:ny-nh+1:-1,nx:nx-nh+1:-1)
    endif

    flag_nw_n = .false.
    flag_nw_w = .false.
    if (northwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvNW,nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            senwtag,MPI_COMM_WORLD,req(8),ierr)
       comm(8)=8
    elseif ( ((lbc).and.(west.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')).or. &
         ((lbc).and.(trim(cuv)=='v')) ) then
!       a2d(ny+1:ny+nh,1-nh:0) = zero
    elseif (north.ne.MPI_PROC_NULL) then
       flag_nw_n = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_nw_w = .true.
    else !!Homogenous Neumann  
       a2D(ny+1:ny+nh,1-nh:0) = a2D(ny:ny-nh+1:-1,nh:1:-1)
    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS(:,:) = a2D(1:nh,1:nx)  
       call MPI_ISend(                              &
            sendS,nx*nh,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE(:,:) = a2D(1:ny,nx-nh+1:nx) 
       call MPI_ISend(                             &
            sendE,ny*nh,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN(:,:) = a2D(ny-nh+1:ny,1:nx)
       call MPI_ISend(                              &
            sendN,nx*nh,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW(:,:) = a2D(1:ny,1:nh)  
       call MPI_ISend(                             &
            sendW,ny*nh,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif


    if (southwest.ne.MPI_PROC_NULL) then
       sendSW(:,:) = a2D(1:nh,1:nh)  
       call MPI_ISend(                                   &
            sendSW,nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            swnetag,MPI_COMM_WORLD,req(13),ierr)
       comm(13)=13
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE(:,:) = a2D(1:nh,nx-nh+1:nx)  
       call MPI_ISend(                                   &
            sendSE,nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            senwtag,MPI_COMM_WORLD,req(14),ierr)
       comm(14)=14
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE(:,:) = a2D(ny-nh+1:ny,nx-nh+1:nx) 
       call MPI_ISend(                                   &
            sendNE,nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            neswtag,MPI_COMM_WORLD,req(15),ierr)
       comm(15)=15
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       sendNW(:,:) = a2D(ny-nh+1:ny,1:nh)
       call MPI_ISend(                                   &
            sendNW,nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            nwsetag,MPI_COMM_WORLD,req(16),ierr)
       comm(16)=16
    endif


    !- Wait for completion of receive and fill ghost points

    icount=0                       ! Compress arrays "comm" and
    do i=1,16                      ! "req" to disregard directions
       if (comm(i).gt.0) then      ! in which no message was sent
          icount=icount+1          ! or is expected from.  At the
          if (icount.lt.i) then    ! end of this segment icount
             comm(icount)=comm(i)  ! is equal to the actual number
             req(icount)=req(i)    ! of messages sent and received, 
          endif                    ! arrays comm,req(1:icount)
       endif                       ! store directional indices
    enddo

    do while (icount > 0)

       call MPI_Waitany(icount, req, j, status, ierr)

       indx=comm(j)           ! Save directional index for
       icount=icount-1        ! message received and ready to
       do i=j,icount          ! unpack, then erase its "req"
          req(i)=req(i+1)     ! and "comm" and "req" by 
          comm(i)=comm(i+1)   ! by compressing the arrays, so
       enddo                  ! that the same message will 

       ! be unpacked only once.
       if (indx.eq.1) then ! south
          a2D(1-nh:0,1:nx)  = recvS(:,:)

       elseif (indx.eq.2) then ! east
          a2D(1:ny,nx+1:nx+nh) = recvE(:,:)

       elseif (indx.eq.3) then ! north
          a2D(ny+1:ny+nh,1:nx)  = recvN (:,:)

       elseif (indx.eq.4) then ! west
          a2D(1:ny,1-nh:0) = recvW(:,:)

       elseif (indx.eq.5) then ! southwest
          a2D(1-nh:0,1-nh:0) = recvSW(:,:)

       elseif (indx.eq.6) then ! southeast
          a2D(1-nh:0,nx+1:nx+nh) = recvSE(:,:)

       elseif (indx.eq.7) then ! northeast
          a2D(ny+1:ny+nh,nx+1:nx+nh) = recvNE(:,:)

       elseif (indx.eq.8) then ! northwest
          a2D(ny+1:ny+nh,1-nh:0) = recvNW(:,:)

       endif

    enddo      !<-- while

    !- corners on physical boundarie if a flag is true-!
    if (flag_sw_s) then
       a2D(1-nh:0,1-nh:0) = a2D(1-nh:0,nh:1:-1)
    elseif (flag_sw_w) then
       a2D(1-nh:0,1-nh:0) = a2D(nh:1:-1,1-nh:0)
    endif

    if (flag_se_s) then
       a2D(1-nh:0,nx+1:nx+nh) = a2D(1-nh:0,nx:nx-nh+1:-1)
    elseif (flag_se_e) then
       a2D(1-nh:0,nx+1:nx+nh) = a2D(nh:1:-1,nx+1:nx+nh)
    endif

    if (flag_ne_n) then
       a2D(ny+1:ny+nh,nx+1:nx+nh) = a2D(ny+1:ny+nh,nx:nx-nh+1:-1)
    elseif (flag_ne_e) then
       a2D(ny+1:ny+nh,nx+1:nx+nh) = a2D(ny:ny-nh+1:-1,nx+1:nx+nh)
    endif

    if (flag_nw_n) then
       a2D(ny+1:ny+nh,1-nh:0) = a2D(ny+1:ny+nh,nh:1:-1)
    elseif (flag_nw_w) then
       a2D(ny+1:ny+nh,1-nh:0) = a2D(ny:ny-nh+1:-1,1-nh:0)
    endif


  end subroutine fill_halo_2D

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_3D_relax(lev,p,nx,ny,nz)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(nz,0:ny+1,0:nx+1), intent(inout)::p

    integer(kind=ip) :: nx, ny, nz
    integer(kind=ip) :: nh
    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    logical :: flag_sw_s, flag_sw_w
    logical :: flag_se_s, flag_se_e
    logical :: flag_ne_n, flag_ne_e
    logical :: flag_nw_n, flag_nw_w

    real(kind=rp), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = size(p,dim=1)

    nh = 1

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)


    if (nz == grid(lev)%nz) then

       sendS => gbuffers(lev)%sendS
       recvS => gbuffers(lev)%recvS
       sendN => gbuffers(lev)%sendN
       recvN => gbuffers(lev)%recvN
       sendE => gbuffers(lev)%sendE
       recvE => gbuffers(lev)%recvE
       sendW => gbuffers(lev)%sendW
       recvW => gbuffers(lev)%recvW

       sendSW => gbuffers(lev)%sendSW
       sendSE => gbuffers(lev)%sendSE
       sendNW => gbuffers(lev)%sendNW
       sendNE => gbuffers(lev)%sendNE
       recvSW => gbuffers(lev)%recvSW
       recvSE => gbuffers(lev)%recvSE
       recvNW => gbuffers(lev)%recvNW
       recvNE => gbuffers(lev)%recvNE

    elseif (nz == (grid(lev)%nz+1)) then

       sendS => gbuffers(lev)%sendSp
       recvS => gbuffers(lev)%recvSp
       sendN => gbuffers(lev)%sendNp
       recvN => gbuffers(lev)%recvNp
       sendE => gbuffers(lev)%sendEp
       recvE => gbuffers(lev)%recvEp
       sendW => gbuffers(lev)%sendWp
       recvW => gbuffers(lev)%recvWp

       sendSW => gbuffers(lev)%sendSWp
       sendSE => gbuffers(lev)%sendSEp
       sendNW => gbuffers(lev)%sendNWp
       sendNE => gbuffers(lev)%sendNEp
       recvSW => gbuffers(lev)%recvSWp
       recvSE => gbuffers(lev)%recvSEp
       recvNW => gbuffers(lev)%recvNWp
       recvNE => gbuffers(lev)%recvNEp

    else
       STOP
    endif

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103
    swnetag = 104
    senwtag = 105
    neswtag = 106
    nwsetag = 107

    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvS,nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    else !!Homogenous Neumann  
       p(:,1-nh:0,1:nx) = p(:,nh:1:-1,1:nx)
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvE,nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    else  !!Homogenous Neumann
       p(:,1:ny,nx+1:nx+nh) = p(:,1:ny,nx:nx-nh+1:-1)
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvN,nz*nx*nh,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    else  !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1:nx) = p(:,ny:ny-nh+1:-1,1:nx)
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvW,nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    else   !!Homogenous Neumann
       p(:,1:ny,1-nh:0) = p(:,1:ny,nh:1:-1)
    endif

    flag_sw_s = .false.
    flag_sw_w = .false.
    if (southwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            neswtag,MPI_COMM_WORLD,req(5),ierr)
       comm(5)=5
    elseif (south.ne.MPI_PROC_NULL) then
       flag_sw_s = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_sw_w = .true.
    else !!Homogenous Neumann  
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,nh:1:-1)
    endif

    flag_se_s = .false.
    flag_se_e = .false.
    if (southeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            nwsetag,MPI_COMM_WORLD,req(6),ierr)
       comm(6)=6
    elseif (south.ne.MPI_PROC_NULL) then
       flag_se_s = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_se_e = .true.
    else!!Homogenous Neumann  
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx:nx-nh+1:-1)
    endif

    flag_ne_n = .false.
    flag_ne_e = .false.
    if (northeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            swnetag,MPI_COMM_WORLD,req(7),ierr)
       comm(7)=7
    elseif (north.ne.MPI_PROC_NULL) then
       flag_ne_n = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_ne_e = .true.
    else!!Homogenous Neumann  
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx:nx-nh+1:-1)
    endif

    flag_nw_n = .false.
    flag_nw_w = .false.
    if (northwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            senwtag,MPI_COMM_WORLD,req(8),ierr)
       comm(8)=8
    elseif (north.ne.MPI_PROC_NULL) then
       flag_nw_n = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_nw_w = .true.
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,nh:1:-1)
    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS = p(:,1:nh,1:nx)  
       call MPI_ISend(                              &
            sendS,nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = p(:,1:ny,nx-nh+1:nx) 
       call MPI_ISend(                             &
            sendE,nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN = p(:,ny-nh+1:ny,1:nx)
       call MPI_ISend(                              &
            sendN,nz*nx,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW = p(:,1:ny,1:nh)  
       call MPI_ISend(                             &
            sendW,nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif

    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = p(:,1:nh,1:nh)  
       call MPI_ISend(                                &
            sendSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            swnetag,MPI_COMM_WORLD,req(13),ierr)
       comm(13)=13
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = p(:,1:nh,nx-nh+1:nx)  
       call MPI_ISend(                                &
            sendSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            senwtag,MPI_COMM_WORLD,req(14),ierr)
       comm(14)=14
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = p(:,ny-nh+1:ny,nx-nh+1:nx) 
       call MPI_ISend(                                &
            sendNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            neswtag,MPI_COMM_WORLD,req(15),ierr)
       comm(15)=15
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = p(:,ny-nh+1:ny,1:nh)
       call MPI_ISend(                                &
            sendNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            nwsetag,MPI_COMM_WORLD,req(16),ierr)
       comm(16)=16
    endif

    !- Wait for completion of receive and fill ghost points

    icount=0                       ! Compress arrays "comm" and
    do i=1,16                      ! "req" to disregard directions
       if (comm(i).gt.0) then      ! in which no message was sent
          icount=icount+1          ! or is expected from.  At the
          if (icount.lt.i) then    ! end of this segment icount
             comm(icount)=comm(i)  ! is equal to the actual number
             req(icount)=req(i)    ! of messages sent and received, 
          endif                    ! arrays comm,req(1:icount)
       endif                       ! store directional indices
    enddo

    do while (icount > 0)

       call MPI_Waitany(icount, req, j, status, ierr)

       indx=comm(j)           ! Save directional index for
       icount=icount-1        ! message received and ready to
       do i=j,icount          ! unpack, then erase its "req"
          req(i)=req(i+1)     ! and "comm" and "req" by 
          comm(i)=comm(i+1)   ! by compressing the arrays, so
       enddo                  ! that the same message will 

       ! be unpacked only once.
       if (indx.eq.1) then ! south
          p(:,1-nh:0,1:nx)  = recvS

       elseif (indx.eq.2) then ! east
          p(:,1:ny,nx+1:nx+nh) = recvE

       elseif (indx.eq.3) then ! north
          p(:,ny+1:ny+nh,1:nx)  = recvN 

       elseif (indx.eq.4) then ! west
          p(:,1:ny,1-nh:0) = recvW

       elseif (indx.eq.5) then ! southwest
          p(:,1-nh:0,1-nh:0) = recvSW

       elseif (indx.eq.6) then ! southeast
          p(:,1-nh:0,nx+1:nx+nh) = recvSE

       elseif (indx.eq.7) then ! northeast
          p(:,ny+1:ny+nh,nx+1:nx+nh) = recvNE

       elseif (indx.eq.8) then ! northwest
          p(:,ny+1:ny+nh,1-nh:0) = recvNW
       endif

    enddo      !<-- while

    !- corners on physical boundarie if a flag is true-!
    if (flag_sw_s) then
       p(:,1-nh:0,1-nh:0) = p(:,1-nh:0,nh:1:-1)
    elseif (flag_sw_w) then
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,1-nh:0)
    endif

    if (flag_se_s) then
       p(:,1-nh:0,nx+1:nx+nh) = p(:,1-nh:0,nx:nx-nh+1:-1)
    elseif (flag_se_e) then
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx+1:nx+nh)
    endif

    if (flag_ne_n) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny+1:ny+nh,nx:nx-nh+1:-1)
    elseif (flag_ne_e) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx+1:nx+nh)
    endif

    if (flag_nw_n) then
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny+1:ny+nh,nh:1:-1)
    elseif (flag_nw_w) then
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,1-nh:0)
    endif

  end subroutine fill_halo_3D_relax

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_3D(lev,p,lbc_null)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)::p
    character(len=1), optional, intent(in) :: lbc_null

    integer(kind=ip) :: nx, ny, nz
    integer(kind=ip) :: nh, ih
    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    logical :: lbc
    character(len=1) :: cuv

    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    logical :: flag_sw_s, flag_sw_w
    logical :: flag_se_s, flag_se_e
    logical :: flag_ne_n, flag_ne_e
    logical :: flag_nw_n, flag_nw_w

    real(kind=rp), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    if (present(lbc_null)) then
       cuv=trim(lbc_null)
       lbc=.true.
    else
       lbc=.false.
    endif

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = size(p,dim=1)

    nh = (size(p,dim=2) - ny ) / 2

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    if (nh == 1) then
       if (nz == grid(lev)%nz) then

          sendS => gbuffers(lev)%sendS
          recvS => gbuffers(lev)%recvS
          sendN => gbuffers(lev)%sendN
          recvN => gbuffers(lev)%recvN
          sendE => gbuffers(lev)%sendE
          recvE => gbuffers(lev)%recvE
          sendW => gbuffers(lev)%sendW
          recvW => gbuffers(lev)%recvW

          sendSW => gbuffers(lev)%sendSW
          sendSE => gbuffers(lev)%sendSE
          sendNW => gbuffers(lev)%sendNW
          sendNE => gbuffers(lev)%sendNE
          recvSW => gbuffers(lev)%recvSW
          recvSE => gbuffers(lev)%recvSE
          recvNW => gbuffers(lev)%recvNW
          recvNE => gbuffers(lev)%recvNE

       elseif (nz == (grid(lev)%nz+1)) then

          sendS => gbuffers(lev)%sendSp
          recvS => gbuffers(lev)%recvSp
          sendN => gbuffers(lev)%sendNp
          recvN => gbuffers(lev)%recvNp
          sendE => gbuffers(lev)%sendEp
          recvE => gbuffers(lev)%recvEp
          sendW => gbuffers(lev)%sendWp
          recvW => gbuffers(lev)%recvWp

          sendSW => gbuffers(lev)%sendSWp
          sendSE => gbuffers(lev)%sendSEp
          sendNW => gbuffers(lev)%sendNWp
          sendNE => gbuffers(lev)%sendNEp
          recvSW => gbuffers(lev)%recvSWp
          recvSE => gbuffers(lev)%recvSEp
          recvNW => gbuffers(lev)%recvNWp
          recvNE => gbuffers(lev)%recvNEp

       else
          write(*,*) 'mg_mpi_exchange: fill_halo_3D: pb 1 !'
          STOP
       endif

    elseif (nh == 2) then

       if (nz == grid(lev)%nz) then

          sendS => gbuffers(lev)%sendS3D2
          recvS => gbuffers(lev)%recvS3D2
          sendN => gbuffers(lev)%sendN3D2
          recvN => gbuffers(lev)%recvN3D2
          sendE => gbuffers(lev)%sendE3D2
          recvE => gbuffers(lev)%recvE3D2
          sendW => gbuffers(lev)%sendW3D2
          recvW => gbuffers(lev)%recvW3D2

          sendSW => gbuffers(lev)%sendSW3D2
          sendSE => gbuffers(lev)%sendSE3D2
          sendNW => gbuffers(lev)%sendNW3D2
          sendNE => gbuffers(lev)%sendNE3D2
          recvSW => gbuffers(lev)%recvSW3D2
          recvSE => gbuffers(lev)%recvSE3D2
          recvNW => gbuffers(lev)%recvNW3D2
          recvNE => gbuffers(lev)%recvNE3D2

       elseif (nz == (grid(lev)%nz+1)) then

          sendS => gbuffers(lev)%sendS3D2p
          recvS => gbuffers(lev)%recvS3D2p
          sendN => gbuffers(lev)%sendN3D2p
          recvN => gbuffers(lev)%recvN3D2p
          sendE => gbuffers(lev)%sendE3D2p
          recvE => gbuffers(lev)%recvE3D2p
          sendW => gbuffers(lev)%sendW3D2p
          recvW => gbuffers(lev)%recvW3D2p

          sendSW => gbuffers(lev)%sendSW3D2p
          sendSE => gbuffers(lev)%sendSE3D2p
          sendNW => gbuffers(lev)%sendNW3D2p
          sendNE => gbuffers(lev)%sendNE3D2p
          recvSW => gbuffers(lev)%recvSW3D2p
          recvSE => gbuffers(lev)%recvSE3D2p
          recvNW => gbuffers(lev)%recvNW3D2p
          recvNE => gbuffers(lev)%recvNE3D2p

       else
          write(*,*) 'mg_mpi_exchange: fill_halo_3D: pb 2 !'
          STOP
       endif

    else
       write(*,*) 'mg_mpi_exchange: fill_halo_3D: pb 3 !', nh,nx,ny,nz
       stop
    endif

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103
    swnetag = 104
    senwtag = 105
    neswtag = 106
    nwsetag = 107

    !-
    ! Fill_halo nh=2 for ZR and ZW 
    !- If you don't have a neighbour
    !- FOR Homogenous Neumann CONDITIONS
    !-
    ! With a halo of 1
    ! ... | D | C   | B      | A      | *           <= initial state
    ! ... | D | C   | B      | A      | A           <= fill_halo 1
    ! ... |   | B­D | A­C    | A­B    | A           <= calculation x +1  – x ­1  (1 to n)
    ! ... |   | B­D | A­C    | A­B    | A­B         <= fill_halo 1
    ! ... |   |     | A­2B+D | ­B+C   | A­B         <= calculation x +1  – x ­1  (1 to n)
    !
    ! With a halo of 2
    ! ... | D | C   | B      | A      | *     | *   <= initial state
    ! ... | D | C   | B      | A      | A     | ?   <= fill_halo 2
    ! ... |   | B­D | A­C    | A­B    | ?­A   | ?   <= calculation x +1  – x ­1  (0 to n+1)
    ! ... |   |     | A­2B+D | ?­2A+C | ?­A+B | ?   <= calculation x +1  – x ­1  (0 to n+1)
    !
    ! What is the value of ?: 
    ! ? ­ 2A + C = ­B + C
    ! ? = +2A – C – B + C
    ! ? = 2A – B
    !
    ! Verification: (halo of 2)
    ! ... | D | C   | B      | A      | *     | *      <= initial state
    ! ... | D | C   | B      | A      | A     | 2A­B   <= fill_halo 2
    ! ... |   | B­D | A­C    | A­B    | A­B   | 2A­B   <= calculation x +1  – x ­1  (0 to n+1) in define_matrix
    ! ... |   |     | A­2B+D | ­B+C   | A     | 2A­B   <= calculation x +1  – x ­1  (0 to n+1) in define_matrix
    !-

    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvS,nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    elseif ((lbc).and.(trim(cuv)=='v')) then
       p(:,1,:) = zero
    else !!Homogenous Neumann
!!$       p(:,1-nh:0,1:nx) = p(:,nh:1:-1,1:nx)
       do ih = 1, nh
          if (ih == 1) then
             p(:,0,1:nx) = p(:,1,1:nx)
          elseif(ih == 2) then
             p(:,-1,1:nx) = two * p(:,1,1:nx) - p(:,2,1:nx)
          else
             stop
          endif
       enddo
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvE,nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    elseif ((lbc).and.(trim(cuv)=='u')) then
       p(:,:,nx+1) = zero
    else  !!Homogenous Neumann
!!$       p(:,1:ny,nx+1:nx+nh) = p(:,1:ny,nx:nx-nh+1:-1)
       do ih = 1, nh
          if (ih == 1) then
             p(:,1:ny,nx+1) = p(:,1:ny,nx)
          elseif(ih == 2) then
             p(:,1:ny,nx+2) = two * p(:,1:ny,nx) - p(:,1:ny,nx-1)
          else
             stop
          endif
       enddo
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvN,nz*nx*nh,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    elseif ((lbc).and.(trim(cuv)=='v')) then 
       p(:,ny+1,:) = zero
    else  !!Homogenous Neumann  
!!$       p(:,ny+1:ny+nh,1:nx) = p(:,ny:ny-nh+1:-1,1:nx)
       do ih = 1, nh
          if (ih == 1) then
             p(:,ny+1,1:nx) = p(:,ny,1:nx)
          elseif(ih == 2) then
             p(:,ny+2,1:nx) = two * p(:,ny,1:nx) - p(:,ny-1,1:nx)
          else
             stop
          endif
       enddo
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvW,nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    elseif ((lbc).and.(trim(cuv)=='u')) then
       p(:,:,1) = zero
    else   !!Homogenous Neumann
!!$       p(:,1:ny,1-nh:0) = p(:,1:ny,nh:1:-1)
       do ih = 1, nh
          if (ih == 1) then
             p(:,1:ny,0) = p(:,1:ny,1)
          elseif(ih == 2) then
             p(:,1:ny,-1) = two * p(:,1:ny,1) - p(:,1:ny,2)
          else
             stop
          endif
       enddo
    endif

    flag_sw_s = .false.
    flag_sw_w = .false.
    if (southwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            neswtag,MPI_COMM_WORLD,req(5),ierr)
       comm(5)=5
    elseif ((lbc).and.(west.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')) then
!       p(:,1-nh:0,1-nh:0) = zero
    elseif (south.ne.MPI_PROC_NULL) then
       flag_sw_s = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_sw_w = .true.
    else !!Homogenous Neumann  
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,nh:1:-1)
    endif

    flag_se_s = .false.
    flag_se_e = .false.
    if (southeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            nwsetag,MPI_COMM_WORLD,req(6),ierr)
       comm(6)=6
    elseif ((lbc).and.(east.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')) then
!       p(:,1-nh:0,nx+1:nx+nh) = zero
    elseif (south.ne.MPI_PROC_NULL) then
       flag_se_s = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_se_e = .true.
    else!!Homogenous Neumann  
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx:nx-nh+1:-1)
    endif

    flag_ne_n = .false.
    flag_ne_e = .false.
    if (northeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            swnetag,MPI_COMM_WORLD,req(7),ierr)
       comm(7)=7
    elseif (((lbc).and.(east.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')).or. &
         ((lbc).and.(trim(cuv)=='v')) ) then
!       p(:,ny+1:ny+nh,nx+1:nx+nh) = zero
    elseif (north.ne.MPI_PROC_NULL) then
       flag_ne_n = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_ne_e = .true.
    else!!Homogenous Neumann  
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx:nx-nh+1:-1)
    endif

    flag_nw_n = .false.
    flag_nw_w = .false.
    if (northwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            senwtag,MPI_COMM_WORLD,req(8),ierr)
       comm(8)=8
    elseif ( ((lbc).and.(west.eq.MPI_PROC_NULL).and.(trim(cuv)=='u')).or. &
         ((lbc).and.(trim(cuv)=='v')) ) then
!       p(:,ny+1:ny+nh,1-nh:0) = zero
    elseif (north.ne.MPI_PROC_NULL) then
       flag_nw_n = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_nw_w = .true.
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,nh:1:-1)
    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS = p(:,1:nh,1:nx)  
       call MPI_ISend(                              &
            sendS,nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = p(:,1:ny,nx-nh+1:nx) 
       call MPI_ISend(                             &
            sendE,nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN = p(:,ny-nh+1:ny,1:nx)
       call MPI_ISend(                              &
            sendN,nz*nx*nh,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW = p(:,1:ny,1:nh)  
       call MPI_ISend(                             &
            sendW,nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif

    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = p(:,1:nh,1:nh)  
       call MPI_ISend(                                &
            sendSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            swnetag,MPI_COMM_WORLD,req(13),ierr)
       comm(13)=13
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = p(:,1:nh,nx-nh+1:nx)  
       call MPI_ISend(                                &
            sendSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            senwtag,MPI_COMM_WORLD,req(14),ierr)
       comm(14)=14
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = p(:,ny-nh+1:ny,nx-nh+1:nx) 
       call MPI_ISend(                                &
            sendNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            neswtag,MPI_COMM_WORLD,req(15),ierr)
       comm(15)=15
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = p(:,ny-nh+1:ny,1:nh)
       call MPI_ISend(                                &
            sendNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            nwsetag,MPI_COMM_WORLD,req(16),ierr)
       comm(16)=16
    endif

    !- Wait for completion of receive and fill ghost points

    icount=0                       ! Compress arrays "comm" and
    do i=1,16                      ! "req" to disregard directions
       if (comm(i).gt.0) then      ! in which no message was sent
          icount=icount+1          ! or is expected from.  At the
          if (icount.lt.i) then    ! end of this segment icount
             comm(icount)=comm(i)  ! is equal to the actual number
             req(icount)=req(i)    ! of messages sent and received, 
          endif                    ! arrays comm,req(1:icount)
       endif                       ! store directional indices
    enddo

    do while (icount > 0)

       call MPI_Waitany(icount, req, j, status, ierr)

       indx=comm(j)           ! Save directional index for
       icount=icount-1        ! message received and ready to
       do i=j,icount          ! unpack, then erase its "req"
          req(i)=req(i+1)     ! and "comm" and "req" by 
          comm(i)=comm(i+1)   ! by compressing the arrays, so
       enddo                  ! that the same message will 

       ! be unpacked only once.
       if (indx.eq.1) then ! south
          p(:,1-nh:0,1:nx)  = recvS

       elseif (indx.eq.2) then ! east
          p(:,1:ny,nx+1:nx+nh) = recvE

       elseif (indx.eq.3) then ! north
          p(:,ny+1:ny+nh,1:nx)  = recvN 

       elseif (indx.eq.4) then ! west
          p(:,1:ny,1-nh:0) = recvW

       elseif (indx.eq.5) then ! southwest
          p(:,1-nh:0,1-nh:0) = recvSW

       elseif (indx.eq.6) then ! southeast
          p(:,1-nh:0,nx+1:nx+nh) = recvSE

       elseif (indx.eq.7) then ! northeast
          p(:,ny+1:ny+nh,nx+1:nx+nh) = recvNE

       elseif (indx.eq.8) then ! northwest
          p(:,ny+1:ny+nh,1-nh:0) = recvNW
       endif

    enddo      !<-- while

    !- corners on physical boundarie if a flag is true-!
    if (flag_sw_s) then
       p(:,1-nh:0,1-nh:0) = p(:,1-nh:0,nh:1:-1)
    elseif (flag_sw_w) then
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,1-nh:0)
    endif

    if (flag_se_s) then
       p(:,1-nh:0,nx+1:nx+nh) = p(:,1-nh:0,nx:nx-nh+1:-1)
    elseif (flag_se_e) then
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx+1:nx+nh)
    endif

    if (flag_ne_n) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny+1:ny+nh,nx:nx-nh+1:-1)
    elseif (flag_ne_e) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx+1:nx+nh)
    endif

    if (flag_nw_n) then
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny+1:ny+nh,nh:1:-1)
    elseif (flag_nw_w) then
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,1-nh:0)
    endif

  end subroutine fill_halo_3D

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_4D(lev,cA)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(:,:,:,:), pointer, intent(inout)::cA

    integer(kind=ip) :: nx, ny, nz, nd
    integer(kind=ip) :: nh
    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    real(kind=rp), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:,:)  , pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:,:)  , pointer :: sendNW,recvNW,sendNE,recvNE

    call tic(lev,'fill_halo_4D')

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = 1

    nd = size(cA(:,:,:,:),dim=1)

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    sendS => gbuffers(lev)%sendS4D
    recvS => gbuffers(lev)%recvS4D
    sendN => gbuffers(lev)%sendN4D
    recvN => gbuffers(lev)%recvN4D

    sendE => gbuffers(lev)%sendE4D
    recvE => gbuffers(lev)%recvE4D
    sendW => gbuffers(lev)%sendW4D
    recvW => gbuffers(lev)%recvW4D

    sendSW => gbuffers(lev)%sendSW4D
    sendSE => gbuffers(lev)%sendSE4D
    sendNW => gbuffers(lev)%sendNW4D
    sendNE => gbuffers(lev)%sendNE4D

    recvSW => gbuffers(lev)%recvSW4D
    recvSE => gbuffers(lev)%recvSE4D
    recvNW => gbuffers(lev)%recvNW4D
    recvNE => gbuffers(lev)%recvNE4D

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103
    swnetag = 104
    senwtag = 105
    neswtag = 106
    nwsetag = 107

    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                 &
            recvS,nd*nz*nx,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    else !!Homogenous Neumann  
!       cA(:,:,0,1:nx) = zero 
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvE,nd*nz*ny,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    else !!Homogenous Neumann
!       cA(:,:,1:ny,nx+1) = zero
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                 &
            recvN,nd*nz*nx,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    else !!Homogenous Neumann  
!       cA(:,:,ny+1,1:nx) = zero
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvW,nd*nz*ny,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    else !!Homogenous Neumann
!       cA(:,:,1:ny,0) = zero
    endif

    if (southwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvSW,nd*nz,MPI_DOUBLE_PRECISION,southwest, &
            neswtag,MPI_COMM_WORLD,req(5),ierr)
       comm(5)=5
    else !!Homogenous Neumann  
!       cA(:,:,0,0) = zero
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvSE,nd*nz,MPI_DOUBLE_PRECISION,southeast, &
            nwsetag,MPI_COMM_WORLD,req(6),ierr)
       comm(6)=6
    else !!Homogenous Neumann  
!       cA(:,:,0,nx+1) = zero
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvNE,nd*nz,MPI_DOUBLE_PRECISION,northeast, &
            swnetag,MPI_COMM_WORLD,req(7),ierr)
       comm(7)=7
    else !!Homogenous Neumann  
!       cA(:,:,ny+1,nx+1) = zero
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvNW,nd*nz,MPI_DOUBLE_PRECISION,northwest, &
            senwtag,MPI_COMM_WORLD,req(8),ierr)
       comm(8)=8
    else !!Homogenous Neumann  
!       cA(:,:,ny+1,0) = zero
    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS = cA(:,:,1,1:nx)  
       call MPI_ISend(                                 &
            sendS,nd*nz*nx,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = cA(:,:,1:ny,nx) 
       call MPI_ISend(                                &
            sendE,nd*nz*ny,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN = cA(:,:,ny,1:nx)
       call MPI_ISend(                                 &
            sendN,nd*nz*nx,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW = cA(:,:,1:ny,1)  
       call MPI_ISend(                                &
            sendW,nd*nz*ny,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif

    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = cA(:,:,1,1)  
       call MPI_ISend(                                   &
            sendSW,nd*nz,MPI_DOUBLE_PRECISION,southwest, &
            swnetag,MPI_COMM_WORLD,req(13),ierr)
       comm(13)=13
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = cA(:,:,1,nx)  
       call MPI_ISend(                                   &
            sendSE,nd*nz,MPI_DOUBLE_PRECISION,southeast, &
            senwtag,MPI_COMM_WORLD,req(14),ierr)
       comm(14)=14
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = cA(:,:,ny,nx) 
       call MPI_ISend(                                   &
            sendNE,nd*nz,MPI_DOUBLE_PRECISION,northeast, &
            neswtag,MPI_COMM_WORLD,req(15),ierr)
       comm(15)=15
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = cA(:,:,ny,1)
       call MPI_ISend(                                   &
            sendNW,nd*nz,MPI_DOUBLE_PRECISION,northwest, &
            nwsetag,MPI_COMM_WORLD,req(16),ierr)
       comm(16)=16
    endif

    !- Wait for completion of receive and fill ghost points

    icount=0                       ! Compress arrays "comm" and
    do i=1,16                      ! "req" to disregard directions
       if (comm(i).gt.0) then      ! in which no message was sent
          icount=icount+1          ! or is expected from.  At the
          if (icount.lt.i) then    ! end of this segment icount
             comm(icount)=comm(i)  ! is equal to the actual number
             req(icount)=req(i)    ! of messages sent and received, 
          endif                    ! arrays comm,req(1:icount)
       endif                       ! store directional indices
    enddo

    do while (icount > 0)

       call MPI_Waitany(icount, req, j, status, ierr)

       indx=comm(j)           ! Save directional index for
       icount=icount-1        ! message received and ready to
       do i=j,icount          ! unpack, then erase its "req"
          req(i)=req(i+1)     ! and "comm" and "req" by 
          comm(i)=comm(i+1)   ! by compressing the arrays, so
       enddo                  ! that the same message will 

       ! be unpacked only once.
       if (indx.eq.1) then ! south
          cA(:,:,0,1:nx)  = recvS

       elseif (indx.eq.2) then ! east
          cA(:,:,1:ny,nx+1) = recvE

       elseif (indx.eq.3) then ! north
          cA(:,:,ny+1,1:nx)  = recvN 

       elseif (indx.eq.4) then ! west
          cA(:,:,1:ny,0) = recvW

       elseif (indx.eq.5) then ! southwest
          cA(:,:,0,0) = recvSW

       elseif (indx.eq.6) then ! southeast
          cA(:,:,0,nx+1) = recvSE

       elseif (indx.eq.7) then ! northeast
          cA(:,:,ny+1,nx+1) = recvNE

       elseif (indx.eq.8) then ! northwest
          cA(:,:,ny+1,0) = recvNW
       endif

    enddo      !<-- while  

    call toc(lev,'fill_halo_4D')

  end subroutine fill_halo_4D

  !----------------------------------------
  subroutine set_rurvbc2zero_3D(a3D,gt)

    real(kind=rp), dimension(:,:,:), pointer, intent(inout)::a3D
    character(len=1) :: gt

    integer(kind=ip) :: nx, ny
    integer(kind=ip) :: south, east, north, west

    nx = grid(1)%nx
    ny = grid(1)%ny

    south     = grid(1)%neighb(1)
    east      = grid(1)%neighb(2)
    north     = grid(1)%neighb(3)
    west      = grid(1)%neighb(4)

    if (trim(gt) == 'u') then

       if (east == MPI_PROC_NULL) then
          a3D(:,:,nx+1) = zero
       endif

       if (west == MPI_PROC_NULL) then
          a3D(:,:,1) = zero
       endif

    elseif (trim(gt) == 'v') then

       if (south == MPI_PROC_NULL) then
          a3D(:,1,:) = zero
       endif

       if (north == MPI_PROC_NULL) then
          a3D(:,ny+1,:) = zero
       endif

    else
       stop
    endif

  end subroutine set_rurvbc2zero_3D

  !----------------------------------------
  subroutine global_max(maxloc)
    ! return the global max: maxglo
    ! using the local max on each subdomain
    real(kind=rp),intent(inout) :: maxloc

    real(kind=rp)    :: maxglo
    integer(kind=ip) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    ! this is not optimal, but not wrong
    call MPI_ALLREDUCE(maxloc,maxglo,1,MPI_DOUBLE_PRECISION,MPI_max,MPI_COMM_WORLD,ierr)   

    maxloc = maxglo

  end subroutine global_max

  !----------------------------------------
  subroutine global_sum(lev,sumloc,sumglo)
    ! return the global sum: sumglo
    ! using the local sum on each subdomain
    integer(kind=ip),intent(in) :: lev
    real(kind=rp),intent(in) :: sumloc
    real(kind=rp),intent(out) :: sumglo

    integer(kind=ip) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    call MPI_ALLREDUCE(sumloc,sumglo,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)   

    ! therefore we need to rescale the global sum
    sumglo = sumglo * (grid(lev)%npx*grid(lev)%npy)/(grid(1)%npx*grid(1)%npy)

  end subroutine global_sum

end module mg_mpi_exchange
