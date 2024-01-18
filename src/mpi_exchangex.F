      module mpi_exchanges
      ! overloaded mpi exchange subroutines

      use param
      use hidden_mpi_vars
      use mess_buffers
      use mpi

      implicit none
      private
      
      interface mpi_exchange_new
        module procedure mpi_exchange1,mpi_exchange2,mpi_exchange3,mpi_exchange4
      end interface

      public :: mpi_exchange_new

      contains

!--------------------------------------------------------
      subroutine mpi_exchange1(A,nx,ny,nz)
      implicit none

      ! input/output
!     integer(kind=4) 
      real(kind=8),dimension(1-bf:nx+bf,1-bf:ny+bf,nz) :: A
      integer(kind=4) i,j,k, kshft, ierr, mess_count, comm(16), req(16),
     &                                  status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      logical flag
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (west .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend==ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr==jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend==jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=4*(istr+iend-2*iwest)/(iend-istr+1)
      jtg=4*(jstr+jend-2*jsouth)/(jend-jstr+1)
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+8, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend==ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif

! Send 4 sides and 4 corners
!----------------------------------------------
      EWsize = nz*ny*bf
      NSsize = nz*nx*bf

      if (west_msg_exch) then
        do k=1,nz
          kshft=bf*ny*(k-1) +1
          do j=1,ny
            jshft=bf*(j-1) +1
            do i=1,bf
              sendW(i+jshft+kshft)=A(i,j,k)
            enddo
          enddo
        enddo
        call MPI_Isend (sendW, EWsize, MPI_DOUBLE_PRECISION,
     &            p_W,7,ocean_grid_comm, req(9), ierr)
        comm(9)=9
      endif

      if (east_msg_exch) then
        do k=1,nz
          kshft=bf*ny*(k-1) +1
          do j=1,ny
            jshft=bf*(j-1) +1
            do i=1,bf
              sendE(i+jshft+kshft)=A(nx+1-i,j,k)
            enddo
          enddo
        enddo
        call MPI_Isend (sendE, EWsize, MPI_DOUBLE_PRECISION,
     &             p_E,8,ocean_grid_comm, req(10), ierr)
        comm(10)=10
      endif

      if (south_msg_exch) then
        do k=1,nz
          kshft=bf*nx*(k-1) +1
          do j=1,bf
            jshft=nx*(j-1) +1
            do i=1,nx
              sendS(i+jshft+kshft)=A(i,j,k)
            enddo
        enddo
        call MPI_Isend (sendS, NSsize, MPI_DOUBLE_PRECISION,
     &         p_S, 5, ocean_grid_comm, req(11), ierr)
        comm(11)=11
      endif

      if (north_msg_exch) then
        do k=1,nz
          kshft=bf*nx*(k-1) +1
          do j=1,bf
            jshft=nx*(j-1) +1
            do i=1,nx
              sendN(i+jshft+kshft)=A(i,ny+1-j,k)
            enddo
          enddo
        enddo
        call MPI_Isend (sendN, NSsize, MPI_DOUBLE_PRECISION,
     &         p_N, 6, ocean_grid_comm, req(12), ierr)
        comm(12)=12
      endif

!  JM 01/11/2024
!  Communicate corners. It seems that often these are not neccesary
!  In addition, for now, we're going to only exchange a single point,
!  since I can only think to the averaging for Coriolis that needs
!  corner averaging (and exchanges for coarsening in surface forcing)

      if (west_msg_exch.and.south_msg_exch) then
        do k=1,nz
          sn_SW(ishft+k)=A(i,1:i,k)
        enddo
      enddo
        call MPI_Isend (sn_SW, Csize, MPI_DOUBLE_PRECISION,
     &           p_SW, 1, ocean_grid_comm, req(13), ierr)
        comm(13)=13
      endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i) > 0) mess_count=mess_count+1
      enddo
      do while (mess_count>0)
        if (comm(1)>0) then
          call MPI_Test (req(1), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(1)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(iwest-2,j,k)=recvW(j-jmin       +kshft)
                A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(2)>0) then
          call MPI_Test (req(2), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(2)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(ieast+1,j,k)=recvE(j-jmin       +kshft)
                A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(3)>0) then
          call MPI_Test (req(3), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(3)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jsouth-2,k)=recvS(i-imin       +kshft)
                A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(4)>0) then
          call MPI_Test (req(4), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(4)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jnorth+1,k)=recvN(i-imin       +kshft)
                A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(5)>0) then
          call MPI_Test (req(5), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(5)=0
            do k=1,nmaxA
              A(iwest-2,jsouth-2,k)=rv_SW(k        )
              A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
              A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
              A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
            enddo
          endif
        endif
        if (comm(6)>0) then
          call MPI_Test (req(6), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(6)=0
            do k=1,nmaxA
              A(ieast+1,jnorth+1,k)=rv_NE(k        )
              A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
              A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
              A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
            enddo
          endif
        endif
        if (comm(7)>0) then
          call MPI_Test (req(7), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(7)=0
            do k=1,nmaxA
              A(ieast+1,jsouth-2,k)=rv_SE(k        )
              A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
              A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
              A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
            enddo
          endif
        endif
        if (comm(8)>0) then
          call MPI_Test (req(8), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(8)=0
            do k=1,nmaxA
              A(iwest-2,jnorth+1,k)=rv_NW(k        )
              A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
              A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
              A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
            enddo
          endif
        endif
        do j=9,16
          if (comm(j) > 0) then
            call MPI_Test (req(j), flag, status, ierr)
            if (flag) then
               mess_count=mess_count-1 ; comm(j)=0
            endif
          endif
        enddo
      enddo
      end














      subroutine mpi_exchange8_2_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                                      B, nmaxB)
      use param
      use hidden_mpi_vars
      use mess_buffers
      use mpi
      implicit none
      integer(kind=4) istr,iend,jstr,jend, nmaxA
      real(kind=8) A(-1:Lm+2,-1:Mm+2,nmaxA)
      integer(kind=4) nmaxB, offset
      real(kind=8) B(-1:Lm+2,-1:Mm+2,nmaxB)
      integer(kind=4) i,j,k, kshft, ierr, mess_count, comm(16), req(16),
     &                                  status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      logical flag
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (istr==iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend==ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr==jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend==jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      ksize=ksize+nmaxB
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=4*(istr+iend-2*iwest)/(iend-istr+1)
      jtg=4*(jstr+jend-2*jsouth)/(jend-jstr+1)
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch.and.istr==iwest) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+8, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend==ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif
      do ipass=0,1
        if (mod(inode+ipass,2)==0) then
          if (west_msg_exch.and.istr==iwest) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=B(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=B(iwest+1,j,k)
              enddo
            enddo
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+7, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
        else
          if (east_msg_exch.and.iend==ieast) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=B(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=B(ieast  ,j,k)
              enddo
            enddo
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+8, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
        endif
      enddo
      do ipass=0,1
        if (mod(jnode+ipass,2)==0) then
          if (south_msg_exch.and.jstr==jsouth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=B(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=B(i,jsouth+1,k)
              enddo
            enddo
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+5, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
        else
          if (north_msg_exch.and.jend==jnorth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=B(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=B(i,jnorth  ,k)
              enddo
            enddo
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+6, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
        endif
      enddo
        if (west_msg_exch.and.istr==iwest .and. 
     &               south_msg_exch.and.jstr==jsouth) then
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SW(k        +offset)=B(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxB+offset)=B(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxB+offset)=B(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxB+offset)=B(iwest+1,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 1, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NE(k        +offset)=B(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxB+offset)=B(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxB+offset)=B(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxB+offset)=B(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SE(k        +offset)=B(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxB+offset)=B(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxB+offset)=B(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxB+offset)=B(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NW(k        +offset)=B(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxB+offset)=B(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxB+offset)=B(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxB+offset)=B(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i) > 0) mess_count=mess_count+1
      enddo
      do while (mess_count>0)
        if (comm(1)>0) then
          call MPI_Test (req(1), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(1)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(iwest-2,j,k)=recvW(j-jmin       +kshft)
                A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(iwest-2,j,k)=recvW(j-jmin       +kshft)
                B(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(2)>0) then
          call MPI_Test (req(2), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(2)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(ieast+1,j,k)=recvE(j-jmin       +kshft)
                A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(ieast+1,j,k)=recvE(j-jmin       +kshft)
                B(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(3)>0) then
          call MPI_Test (req(3), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(3)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jsouth-2,k)=recvS(i-imin       +kshft)
                A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jsouth-2,k)=recvS(i-imin       +kshft)
                B(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(4)>0) then
          call MPI_Test (req(4), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(4)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jnorth+1,k)=recvN(i-imin       +kshft)
                A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jnorth+1,k)=recvN(i-imin       +kshft)
                B(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(5)>0) then
          call MPI_Test (req(5), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(5)=0
            do k=1,nmaxA
              A(iwest-2,jsouth-2,k)=rv_SW(k        )
              A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
              A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
              A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              B(iwest-1,jsouth-2,k)=rv_SW(k + nmaxB+offset)
              B(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxB+offset)
              B(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxB+offset)
            enddo
          endif
        endif
        if (comm(6)>0) then
          call MPI_Test (req(6), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(6)=0
            do k=1,nmaxA
              A(ieast+1,jnorth+1,k)=rv_NE(k        )
              A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
              A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
              A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              B(ieast+2,jnorth+1,k)=rv_NE(k + nmaxB+offset)
              B(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxB+offset)
              B(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxB+offset)
            enddo
          endif
        endif
        if (comm(7)>0) then
          call MPI_Test (req(7), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(7)=0
            do k=1,nmaxA
              A(ieast+1,jsouth-2,k)=rv_SE(k        )
              A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
              A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
              A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              B(ieast+2,jsouth-2,k)=rv_SE(k + nmaxB+offset)
              B(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxB+offset)
              B(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxB+offset)
            enddo
          endif
        endif
        if (comm(8)>0) then
          call MPI_Test (req(8), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(8)=0
            do k=1,nmaxA
              A(iwest-2,jnorth+1,k)=rv_NW(k        )
              A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
              A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
              A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              B(iwest-1,jnorth+1,k)=rv_NW(k + nmaxB+offset)
              B(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxB+offset)
              B(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxB+offset)
            enddo
          endif
        endif
        do j=9,16
          if (comm(j) > 0) then
            call MPI_Test (req(j), flag, status, ierr)
            if (flag) then
               mess_count=mess_count-1 ; comm(j)=0
            endif
          endif
        enddo
      enddo
      end
      subroutine mpi_exchange8_3_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                            B, nmaxB, C, nmaxC)
      use param
      use hidden_mpi_vars
      use mess_buffers
      use mpi
      implicit none
      integer(kind=4) istr,iend,jstr,jend, nmaxA
      real(kind=8) A(-1:Lm+2,-1:Mm+2,nmaxA)
      integer(kind=4) nmaxB, offset
      real(kind=8) B(-1:Lm+2,-1:Mm+2,nmaxB)
      integer(kind=4) nmaxC
      real(kind=8) C(-1:Lm+2,-1:Mm+2,nmaxC)
      integer(kind=4) i,j,k, kshft, ierr, mess_count, comm(16), req(16),
     &                                  status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      logical flag
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (istr==iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend==ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr==jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend==jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      ksize=ksize+nmaxB
      ksize=ksize+nmaxC
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=4*(istr+iend-2*iwest)/(iend-istr+1)
      jtg=4*(jstr+jend-2*jsouth)/(jend-jstr+1)
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch.and.istr==iwest) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+8, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend==ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif
      do ipass=0,1
        if (mod(inode+ipass,2)==0) then
          if (west_msg_exch.and.istr==iwest) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=B(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=B(iwest+1,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=C(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=C(iwest+1,j,k)
              enddo
            enddo
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+7, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
        else
          if (east_msg_exch.and.iend==ieast) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=B(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=B(ieast  ,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=C(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=C(ieast  ,j,k)
              enddo
            enddo
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+8, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
        endif
      enddo
      do ipass=0,1
        if (mod(jnode+ipass,2)==0) then
          if (south_msg_exch.and.jstr==jsouth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=B(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=B(i,jsouth+1,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=C(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=C(i,jsouth+1,k)
              enddo
            enddo
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+5, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
        else
          if (north_msg_exch.and.jend==jnorth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=B(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=B(i,jnorth  ,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=C(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=C(i,jnorth  ,k)
              enddo
            enddo
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+6, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
        endif
      enddo
        if (west_msg_exch.and.istr==iwest .and. 
     &               south_msg_exch.and.jstr==jsouth) then
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SW(k        +offset)=B(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxB+offset)=B(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxB+offset)=B(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxB+offset)=B(iwest+1,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SW(k        +offset)=C(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxC+offset)=C(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxC+offset)=C(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxC+offset)=C(iwest+1,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 1, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NE(k        +offset)=B(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxB+offset)=B(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxB+offset)=B(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxB+offset)=B(ieast  ,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NE(k        +offset)=C(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxC+offset)=C(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxC+offset)=C(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxC+offset)=C(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SE(k        +offset)=B(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxB+offset)=B(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxB+offset)=B(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxB+offset)=B(ieast  ,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SE(k        +offset)=C(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxC+offset)=C(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxC+offset)=C(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxC+offset)=C(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NW(k        +offset)=B(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxB+offset)=B(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxB+offset)=B(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxB+offset)=B(iwest+1,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NW(k        +offset)=C(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxC+offset)=C(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxC+offset)=C(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxC+offset)=C(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i) > 0) mess_count=mess_count+1
      enddo
      do while (mess_count>0)
        if (comm(1)>0) then
          call MPI_Test (req(1), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(1)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(iwest-2,j,k)=recvW(j-jmin       +kshft)
                A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(iwest-2,j,k)=recvW(j-jmin       +kshft)
                B(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                C(iwest-2,j,k)=recvW(j-jmin       +kshft)
                C(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(2)>0) then
          call MPI_Test (req(2), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(2)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(ieast+1,j,k)=recvE(j-jmin       +kshft)
                A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(ieast+1,j,k)=recvE(j-jmin       +kshft)
                B(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                C(ieast+1,j,k)=recvE(j-jmin       +kshft)
                C(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(3)>0) then
          call MPI_Test (req(3), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(3)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jsouth-2,k)=recvS(i-imin       +kshft)
                A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jsouth-2,k)=recvS(i-imin       +kshft)
                B(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                C(i,jsouth-2,k)=recvS(i-imin       +kshft)
                C(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(4)>0) then
          call MPI_Test (req(4), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(4)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jnorth+1,k)=recvN(i-imin       +kshft)
                A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jnorth+1,k)=recvN(i-imin       +kshft)
                B(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                C(i,jnorth+1,k)=recvN(i-imin       +kshft)
                C(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(5)>0) then
          call MPI_Test (req(5), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(5)=0
            do k=1,nmaxA
              A(iwest-2,jsouth-2,k)=rv_SW(k        )
              A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
              A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
              A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              B(iwest-1,jsouth-2,k)=rv_SW(k + nmaxB+offset)
              B(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxB+offset)
              B(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              C(iwest-1,jsouth-2,k)=rv_SW(k + nmaxC+offset)
              C(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxC+offset)
              C(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxC+offset)
            enddo
          endif
        endif
        if (comm(6)>0) then
          call MPI_Test (req(6), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(6)=0
            do k=1,nmaxA
              A(ieast+1,jnorth+1,k)=rv_NE(k        )
              A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
              A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
              A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              B(ieast+2,jnorth+1,k)=rv_NE(k + nmaxB+offset)
              B(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxB+offset)
              B(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              C(ieast+2,jnorth+1,k)=rv_NE(k + nmaxC+offset)
              C(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxC+offset)
              C(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxC+offset)
            enddo
          endif
        endif
        if (comm(7)>0) then
          call MPI_Test (req(7), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(7)=0
            do k=1,nmaxA
              A(ieast+1,jsouth-2,k)=rv_SE(k        )
              A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
              A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
              A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              B(ieast+2,jsouth-2,k)=rv_SE(k + nmaxB+offset)
              B(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxB+offset)
              B(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              C(ieast+2,jsouth-2,k)=rv_SE(k + nmaxC+offset)
              C(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxC+offset)
              C(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxC+offset)
            enddo
          endif
        endif
        if (comm(8)>0) then
          call MPI_Test (req(8), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(8)=0
            do k=1,nmaxA
              A(iwest-2,jnorth+1,k)=rv_NW(k        )
              A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
              A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
              A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              B(iwest-1,jnorth+1,k)=rv_NW(k + nmaxB+offset)
              B(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxB+offset)
              B(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              C(iwest-1,jnorth+1,k)=rv_NW(k + nmaxC+offset)
              C(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxC+offset)
              C(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxC+offset)
            enddo
          endif
        endif
        do j=9,16
          if (comm(j) > 0) then
            call MPI_Test (req(j), flag, status, ierr)
            if (flag) then
               mess_count=mess_count-1 ; comm(j)=0
            endif
          endif
        enddo
      enddo
      end
      subroutine mpi_exchange8_4_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                  B, nmaxB, C, nmaxC, D, nmaxD)
      use param
      use hidden_mpi_vars
      use mess_buffers
      use mpi
      implicit none
      integer(kind=4) istr,iend,jstr,jend, nmaxA
      real(kind=8) A(-1:Lm+2,-1:Mm+2,nmaxA)
      integer(kind=4) nmaxB, offset
      real(kind=8) B(-1:Lm+2,-1:Mm+2,nmaxB)
      integer(kind=4) nmaxC
      real(kind=8) C(-1:Lm+2,-1:Mm+2,nmaxC)
      integer(kind=4) nmaxD
      real(kind=8) D(-1:Lm+2,-1:Mm+2,nmaxD)
      integer(kind=4) i,j,k, kshft, ierr, mess_count, comm(16), req(16),
     &                                  status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      logical flag
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (istr==iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend==ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr==jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend==jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      ksize=ksize+nmaxB
      ksize=ksize+nmaxC
      ksize=ksize+nmaxD
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=4*(istr+iend-2*iwest)/(iend-istr+1)
      jtg=4*(jstr+jend-2*jsouth)/(jend-jstr+1)
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch.and.istr==iwest) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+8, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend==ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+7, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif
      do ipass=0,1
        if (mod(inode+ipass,2)==0) then
          if (west_msg_exch.and.istr==iwest) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=B(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=B(iwest+1,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=C(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=C(iwest+1,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=D(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=D(iwest+1,j,k)
              enddo
            enddo
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+7, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
        else
          if (east_msg_exch.and.iend==ieast) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=B(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=B(ieast  ,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=C(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=C(ieast  ,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=D(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=D(ieast  ,j,k)
              enddo
            enddo
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+8, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
        endif
      enddo
      do ipass=0,1
        if (mod(jnode+ipass,2)==0) then
          if (south_msg_exch.and.jstr==jsouth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=B(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=B(i,jsouth+1,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=C(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=C(i,jsouth+1,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=D(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=D(i,jsouth+1,k)
              enddo
            enddo
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+5, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
        else
          if (north_msg_exch.and.jend==jnorth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=B(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=B(i,jnorth  ,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=C(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=C(i,jnorth  ,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=D(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=D(i,jnorth  ,k)
              enddo
            enddo
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+6, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
        endif
      enddo
        if (west_msg_exch.and.istr==iwest .and. 
     &               south_msg_exch.and.jstr==jsouth) then
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SW(k        +offset)=B(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxB+offset)=B(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxB+offset)=B(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxB+offset)=B(iwest+1,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SW(k        +offset)=C(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxC+offset)=C(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxC+offset)=C(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxC+offset)=C(iwest+1,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_SW(k        +offset)=D(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxD+offset)=D(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxD+offset)=D(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxD+offset)=D(iwest+1,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 1, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NE(k        +offset)=B(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxB+offset)=B(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxB+offset)=B(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxB+offset)=B(ieast  ,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NE(k        +offset)=C(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxC+offset)=C(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxC+offset)=C(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxC+offset)=C(ieast  ,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_NE(k        +offset)=D(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxD+offset)=D(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxD+offset)=D(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxD+offset)=D(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend==ieast .and. 
     &               south_msg_exch.and.jstr==jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SE(k        +offset)=B(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxB+offset)=B(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxB+offset)=B(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxB+offset)=B(ieast  ,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SE(k        +offset)=C(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxC+offset)=C(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxC+offset)=C(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxC+offset)=C(ieast  ,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_SE(k        +offset)=D(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxD+offset)=D(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxD+offset)=D(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxD+offset)=D(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr==iwest .and. 
     &               north_msg_exch.and.jend==jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NW(k        +offset)=B(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxB+offset)=B(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxB+offset)=B(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxB+offset)=B(iwest+1,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NW(k        +offset)=C(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxC+offset)=C(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxC+offset)=C(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxC+offset)=C(iwest+1,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_NW(k        +offset)=D(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxD+offset)=D(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxD+offset)=D(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxD+offset)=D(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i) > 0) mess_count=mess_count+1
      enddo
      do while (mess_count>0)
        if (comm(1)>0) then
          call MPI_Test (req(1), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(1)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(iwest-2,j,k)=recvW(j-jmin       +kshft)
                A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(iwest-2,j,k)=recvW(j-jmin       +kshft)
                B(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                C(iwest-2,j,k)=recvW(j-jmin       +kshft)
                C(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                D(iwest-2,j,k)=recvW(j-jmin       +kshft)
                D(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(2)>0) then
          call MPI_Test (req(2), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(2)=0
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(ieast+1,j,k)=recvE(j-jmin       +kshft)
                A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(ieast+1,j,k)=recvE(j-jmin       +kshft)
                B(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                C(ieast+1,j,k)=recvE(j-jmin       +kshft)
                C(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                D(ieast+1,j,k)=recvE(j-jmin       +kshft)
                D(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(3)>0) then
          call MPI_Test (req(3), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(3)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jsouth-2,k)=recvS(i-imin       +kshft)
                A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jsouth-2,k)=recvS(i-imin       +kshft)
                B(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                C(i,jsouth-2,k)=recvS(i-imin       +kshft)
                C(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                D(i,jsouth-2,k)=recvS(i-imin       +kshft)
                D(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(4)>0) then
          call MPI_Test (req(4), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(4)=0
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                A(i,jnorth+1,k)=recvN(i-imin       +kshft)
                A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jnorth+1,k)=recvN(i-imin       +kshft)
                B(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                C(i,jnorth+1,k)=recvN(i-imin       +kshft)
                C(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                D(i,jnorth+1,k)=recvN(i-imin       +kshft)
                D(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
        if (comm(5)>0) then
          call MPI_Test (req(5), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(5)=0
            do k=1,nmaxA
              A(iwest-2,jsouth-2,k)=rv_SW(k        )
              A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
              A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
              A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              B(iwest-1,jsouth-2,k)=rv_SW(k + nmaxB+offset)
              B(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxB+offset)
              B(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              C(iwest-1,jsouth-2,k)=rv_SW(k + nmaxC+offset)
              C(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxC+offset)
              C(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxC+offset)
            enddo
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              D(iwest-1,jsouth-2,k)=rv_SW(k + nmaxD+offset)
              D(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxD+offset)
              D(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxD+offset)
            enddo
          endif
        endif
        if (comm(6)>0) then
          call MPI_Test (req(6), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(6)=0
            do k=1,nmaxA
              A(ieast+1,jnorth+1,k)=rv_NE(k        )
              A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
              A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
              A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              B(ieast+2,jnorth+1,k)=rv_NE(k + nmaxB+offset)
              B(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxB+offset)
              B(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              C(ieast+2,jnorth+1,k)=rv_NE(k + nmaxC+offset)
              C(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxC+offset)
              C(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxC+offset)
            enddo
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              D(ieast+2,jnorth+1,k)=rv_NE(k + nmaxD+offset)
              D(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxD+offset)
              D(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxD+offset)
            enddo
          endif
        endif
        if (comm(7)>0) then
          call MPI_Test (req(7), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(7)=0
            do k=1,nmaxA
              A(ieast+1,jsouth-2,k)=rv_SE(k        )
              A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
              A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
              A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              B(ieast+2,jsouth-2,k)=rv_SE(k + nmaxB+offset)
              B(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxB+offset)
              B(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              C(ieast+2,jsouth-2,k)=rv_SE(k + nmaxC+offset)
              C(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxC+offset)
              C(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxC+offset)
            enddo
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              D(ieast+2,jsouth-2,k)=rv_SE(k + nmaxD+offset)
              D(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxD+offset)
              D(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxD+offset)
            enddo
          endif
        endif
        if (comm(8)>0) then
          call MPI_Test (req(8), flag, status, ierr)
          if (flag) then
            mess_count=mess_count-1 ; comm(8)=0
            do k=1,nmaxA
              A(iwest-2,jnorth+1,k)=rv_NW(k        )
              A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
              A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
              A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
            enddo
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              B(iwest-1,jnorth+1,k)=rv_NW(k + nmaxB+offset)
              B(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxB+offset)
              B(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxB+offset)
            enddo
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              C(iwest-1,jnorth+1,k)=rv_NW(k + nmaxC+offset)
              C(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxC+offset)
              C(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxC+offset)
            enddo
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              D(iwest-1,jnorth+1,k)=rv_NW(k + nmaxD+offset)
              D(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxD+offset)
              D(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxD+offset)
            enddo
          endif
        endif
        do j=9,16
          if (comm(j) > 0) then
            call MPI_Test (req(j), flag, status, ierr)
            if (flag) then
               mess_count=mess_count-1 ; comm(j)=0
            endif
          endif
        enddo
      enddo
      end
