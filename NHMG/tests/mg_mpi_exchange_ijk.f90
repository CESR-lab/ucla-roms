module mg_mpi_exchange_ijk

  use mg_mpi
  use mg_grids

  implicit none

  interface fill_halo_ijk
     module procedure        &
          fill_halo_2D     , &
          fill_halo_3D_ijk
  end interface fill_halo_ijk

contains

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_2D(nx,ny,npx,npy, a2D)

    integer(kind=ip), intent(in):: nx,ny
    integer(kind=ip), intent(in):: npx,npy
    real(kind=rp), dimension(:,:), pointer, intent(inout)::a2D

    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: sntag, ewtag, nstag, wetag

    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    real(kind=rp), dimension(:), allocatable :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:), allocatable :: sendE,recvE,sendW,recvW

    integer(kind=ip) :: pi, pj

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    pj = myrank/npx
    pi = mod(myrank,npx)

    if (pj >= 1) then ! south
       south = (pj-1)*npx+pi
    else
       south = MPI_PROC_NULL
    endif

    if (pi < npx-1) then ! east
       east = pj*npx+pi+1
    else
       east = MPI_PROC_NULL
    endif

    if (pj < npy-1) then ! north
       north = (pj+1)*npx+pi
    else
       north = MPI_PROC_NULL
    endif

    if (pi >= 1) then ! west
       west = pj*npx+pi-1
    else
       west = MPI_PROC_NULL
    endif

    allocate(sendN(nx))
    allocate(recvN(nx))
    allocate(sendS(nx))
    allocate(recvS(nx))

    allocate(sendE(ny))
    allocate(recvE(ny))
    allocate(sendW(ny))
    allocate(recvW(ny))

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103


    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                           &
            recvS,nx,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    else !!Homogenous Neumann  
       a2D(0,1:nx) = a2D(1,1:nx)
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                          &
            recvE,ny,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    else !!Homogenous Neumann
       a2D(1:ny,nx+1) = a2D(1:ny,nx)
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                           &
            recvN,nx,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    else !!Homogenous Neumann  
       a2D(ny+1,1:nx) = a2D(ny,1:nx)
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                          &
            recvW,ny,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    else !!Homogenous Neumann
       a2D(1:ny,0) = a2D(1:ny,1)
    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS(:) = a2D(1,1:nx)  
       call MPI_ISend(                           &
            sendS,nx,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE(:) = a2D(1:ny,nx) 
       call MPI_ISend(                          &
            sendE,ny,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN(:) = a2D(ny,1:nx)
       call MPI_ISend(                           &
            sendN,nx,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW(:) = a2D(1:ny,1)  
       call MPI_ISend(                          &
            sendW,ny,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
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
          a2D(0,1:nx)  = recvS(:)

       elseif (indx.eq.2) then ! east
          a2D(1:ny,nx+1) = recvE(:)

       elseif (indx.eq.3) then ! north
          a2D(ny+1,1:nx)  = recvN (:)

       elseif (indx.eq.4) then ! west
          a2D(1:ny,0) = recvW(:)
       endif

    enddo      !<-- while

  end subroutine fill_halo_2D


  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_3D_ijk(nx,ny,p,type)

   integer(kind=ip),intent(in) :: nx,ny
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)::p
    character(len=1) :: type

    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    integer(kind=ip) :: lev
    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    integer(kind=ip) :: nz

    logical :: flag_sw_s, flag_sw_w
    logical :: flag_se_s, flag_se_e
    logical :: flag_ne_n, flag_ne_e
    logical :: flag_nw_n, flag_nw_w

    real(kind=rp), dimension(:,:), allocatable :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:), allocatable :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:)  , allocatable :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:)  , allocatable :: sendNW,recvNW,sendNE,recvNE

    lev = 1

    nz = size(p,dim=3)

    flag_sw_s = .false.
    flag_sw_w = .false.
    flag_se_s = .false.
    flag_se_e = .false.
    flag_ne_n = .false.
    flag_ne_e = .false.
    flag_nw_n = .false. 
    flag_nw_w = .false.

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

!!$    if (trim(type)=='u') then
!!$
!!$       allocate(sendN(nx,nz))
!!$       allocate(recvN(nx,nz))
!!$       allocate(sendS(nx,nz))
!!$       allocate(recvS(nx,nz))
!!$
!!$       allocate(sendW(ny,nz))
!!$       allocate(recvE(ny,nz))
!!$
!!$       allocate(sendNW(nz))
!!$       allocate(recvNE(nz))
!!$       allocate(sendSW(nz))
!!$       allocate(recvSE(nz))
!!$
!!$    elseif (trim(type)=='v') then
!!$
!!$       allocate(sendS(nx,nz))
!!$       allocate(recvN(nx,nz))
!!$
!!$       allocate(sendE(ny,nz))
!!$       allocate(recvE(ny,nz))
!!$       allocate(sendW(ny,nz))
!!$       allocate(recvW(ny,nz))
!!$
!!$       allocate(sendSE(nz))
!!$       allocate(recvNE(nz))
!!$       allocate(sendSW(nz))
!!$       allocate(recvNW(nz))
!!$
!!$    elseif (trim(type)=='w') then

       allocate(sendN(nx,nz))
       allocate(recvN(nx,nz))
       allocate(sendS(nx,nz))
       allocate(recvS(nx,nz))

       allocate(sendE(ny,nz))
       allocate(recvE(ny,nz))
       allocate(sendW(ny,nz))
       allocate(recvW(ny,nz))

       allocate(sendNE(nz))
       allocate(recvNE(nz))
       allocate(sendSE(nz))
       allocate(recvSE(nz))

       allocate(sendNW(nz))
       allocate(recvNW(nz))
       allocate(sendSW(nz))
       allocate(recvSW(nz))

!!$    else
!!$       stop
!!$    endif


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

!!$    if ((trim(type)=='u').or.(trim(type)=='w')) then

       if (south.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                              &
               recvS,nz*nx,MPI_DOUBLE_PRECISION,south, &
               nstag,MPI_COMM_WORLD,req(1),ierr)
          comm(1)=1
       else !!Homogenous Neumann  
          p(1:nx,0,:) = p(1:nx,1,:)
       endif

!!$    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvE,nz*ny,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    else !!Homogenous Neumann
       p(nx+1,1:ny,:) = p(nx,1:ny,:)
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvN,nz*nx,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    else !!Homogenous Neumann  
       p(1:nx,ny+1,:) = p(1:nx,ny,:)
    endif

!!$    if ((trim(type)=='v').or.(trim(type)=='w')) then

       if (west.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                             &
               recvW,nz*ny,MPI_DOUBLE_PRECISION,west, &
               ewtag,MPI_COMM_WORLD,req(4),ierr)
          comm(4)=4
       else !!Homogenous Neumann
          p(0,1:ny,:) = p(1,1:ny,:)
       endif

!!$    endif


!!$    if (trim(type)=='w') then

       flag_sw_s = .false.
       flag_sw_w = .false.
       if (southwest.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                                &
               recvSW,nz,MPI_DOUBLE_PRECISION,southwest, &
               neswtag,MPI_COMM_WORLD,req(5),ierr)
          comm(5)=5
       elseif (south.ne.MPI_PROC_NULL) then
          flag_sw_s = .true.
       elseif (west.ne.MPI_PROC_NULL) then
          flag_sw_w = .true.
       else !!Homogenous Neumann  
          p(0,0,:) = p(1,1,:)
       endif

!!$    endif

!!$    if ((trim(type)=='u').or.(trim(type)=='w')) then

       flag_se_s = .false.
       flag_se_e = .false.
       if (southeast.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                                &
               recvSE,nz,MPI_DOUBLE_PRECISION,southeast, &
               nwsetag,MPI_COMM_WORLD,req(6),ierr)
          comm(6)=6
       elseif (south.ne.MPI_PROC_NULL) then
          flag_se_s = .true.
       elseif (east.ne.MPI_PROC_NULL) then
          flag_se_e = .true.
       else !!Homogenous Neumann  
          p(nx+1,0,:) = p(nx,1,:)
       endif

!!$    endif

    flag_ne_n = .false.
    flag_ne_e = .false.
    if (northeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvNE,nz,MPI_DOUBLE_PRECISION,northeast, &
            swnetag,MPI_COMM_WORLD,req(7),ierr)
       comm(7)=7
    elseif (north.ne.MPI_PROC_NULL) then
       flag_ne_n = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_ne_e = .true.
    else !!Homogenous Neumann  
       p(nx+1,ny+1,:) = p(nx,ny,:)
    endif


!!$    if ((trim(type)=='v').or.(trim(type)=='w')) then

       flag_nw_n = .false.
       flag_nw_w = .false.
       if (northwest.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                                &
               recvNW,nz,MPI_DOUBLE_PRECISION,northwest, &
               senwtag,MPI_COMM_WORLD,req(8),ierr)
          comm(8)=8
       elseif (north.ne.MPI_PROC_NULL) then
          flag_nw_n = .true.
       elseif (west.ne.MPI_PROC_NULL) then
          flag_nw_w = .true.
       else !!Homogenous Neumann  
          p(0,ny+1,:) = p(1,ny,:)
       endif

!!$    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS = p(1:nx,1,:)
       call MPI_ISend(                              &
            sendS,nz*nx,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

!!$    if ((trim(type)=='v').or.(trim(type)=='w')) then

       if (east.ne.MPI_PROC_NULL) then
          sendE = p(nx,1:ny,:) 
          call MPI_ISend(                             &
               sendE,nz*ny,MPI_DOUBLE_PRECISION,east, &
               ewtag,MPI_COMM_WORLD,req(10),ierr)
          comm(10)=10
       endif

!!$    endif

!!$    if ((trim(type)=='u').or.(trim(type)=='w')) then

       if (north.ne.MPI_PROC_NULL) then
          sendN = p(1:nx,ny,:)
          call MPI_ISend(                              &
               sendN,nz*nx,MPI_DOUBLE_PRECISION,north, &
               nstag,MPI_COMM_WORLD,req(11),ierr)
          comm(11)=11
       endif

!!$    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW = p(1,1:ny,:)  
       call MPI_ISend(                             &
            sendW,nz*ny,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif


    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = p(1,1,:)  
       call MPI_ISend(                                &
            sendSW,nz,MPI_DOUBLE_PRECISION,southwest, &
            swnetag,MPI_COMM_WORLD,req(13),ierr)
       comm(13)=13
    endif


!!$    if ((trim(type)=='v').or.(trim(type)=='w')) then

       if (southeast.ne.MPI_PROC_NULL) then
          sendSE = p(nx,1,:)  
          call MPI_ISend(                                &
               sendSE,nz,MPI_DOUBLE_PRECISION,southeast, &
               senwtag,MPI_COMM_WORLD,req(14),ierr)
          comm(14)=14
       endif

!!$    endif

!!$    if (trim(type)=='w') then

       if (northeast.ne.MPI_PROC_NULL) then
          sendNE = p(nx,ny,:) 
          call MPI_ISend(                                &
               sendNE,nz,MPI_DOUBLE_PRECISION,northeast, &
               neswtag,MPI_COMM_WORLD,req(15),ierr)
          comm(15)=15
       endif

!!$    endif

 !!$   if ((trim(type)=='u').or.(trim(type)=='w')) then

       if (northwest.ne.MPI_PROC_NULL) then
          sendNW = p(1,ny,:)
          call MPI_ISend(                                &
               sendNW,nz,MPI_DOUBLE_PRECISION,northwest, &
               nwsetag,MPI_COMM_WORLD,req(16),ierr)
          comm(16)=16
       endif

!!$    endif

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
          p(1:nx,0,:)  = recvS

       elseif (indx.eq.2) then ! east
          p(nx+1,1:ny,:) = recvE

       elseif (indx.eq.3) then ! north
          p(1:nx,ny+1,:)  = recvN 

       elseif (indx.eq.4) then ! west
          p(0,1:ny,:) = recvW

       elseif (indx.eq.5) then ! southwest
          p(0,0,:) = recvSW

       elseif (indx.eq.6) then ! southeast
          p(nx+1,0,:) = recvSE

       elseif (indx.eq.7) then ! northeast
          p(nx+1,ny+1,:) = recvNE

       elseif (indx.eq.8) then ! northwest
          p(0,ny+1,:) = recvNW
       endif

    enddo      !<-- while

    !- corners on physicical boundarie if a flag is true-!
    if (flag_sw_s) then
       p(0,0,:) = p(1,0,:)
    elseif (flag_sw_w) then
       p(0,0,:) = p(0,1,:)
    endif

    if (flag_se_s) then
       p(nx+1,0,:) = p(nx,0,:)
    elseif (flag_se_e) then
       p(nx+1,0,:) = p(nx+1,1,:)
    endif

    if (flag_ne_n) then
       p(nx+1,ny+1,:) = p(nx,ny+1,:)
    elseif (flag_ne_e) then
       p(nx+1,ny+1,:) = p(nx+1,ny,:)
    endif

    if (flag_nw_n) then
       p(0,ny+1,:) = p(1,ny+1,:)
    elseif (flag_nw_w) then
       p(0,ny+1,:) = p(0,ny,:)
    endif

    if (allocated(sendN)) deallocate(sendN)
    if (allocated(recvN)) deallocate(recvN)

    if (allocated(sendS)) deallocate(sendS)
    if (allocated(recvS)) deallocate(recvS)

    if (allocated(sendE)) deallocate(sendE)
    if (allocated(recvE)) deallocate(recvE)

    if (allocated(sendW)) deallocate(sendW)
    if (allocated(recvW)) deallocate(recvW)

    if (allocated(sendNE)) deallocate(sendNE)
    if (allocated(recvNE)) deallocate(recvNE)

    if (allocated(sendSE)) deallocate(sendSE)
    if (allocated(recvSE)) deallocate(recvSE)

    if (allocated(sendNW)) deallocate(sendNW)
    if (allocated(recvNW)) deallocate(recvNW)

    if (allocated(sendSW)) deallocate(sendSW)
    if (allocated(recvSW)) deallocate(recvSW)


  end subroutine fill_halo_3D_ijk

end module mg_mpi_exchange_ijk
