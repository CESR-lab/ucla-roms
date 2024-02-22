/* This package contains MPI halo exchange subroutines designed to
 communicate 2-point-wide computational margins of from 1 to up to 4
 arrays at once. All four versions are generated from the same code via
 self-expansion controlled by multiply defined CPP-macro NARGS.

 The two other CPP-switches control the following:

 NONBLK_SEND introduces the use of non-blocking MPI_Isend -- with
 subsequent Wait (combined into Waitany for Irecvs) instead of the
 generic MPI_Send;

 ALT_SEND_ORD introduces east-west vs. west-east alternation for
 sending for even and odd nodes in XI-direction, and, similar
 alternation in ETA-direction.

 Settings of these two switches does not affect the result, but may
 affect performance.  */

#ifndef NARGS
# include "cppdefs.opt"
# define NONBLK_SEND
# define ALT_SEND_ORD
#endif

#ifdef MPI
# ifndef NARGS
      subroutine mpi_exchange8_tile (istr,iend,jstr,jend, A, nmaxA)
# elif NARGS == 2
      subroutine mpi_exchange8_2_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                                      B, nmaxB)
# elif NARGS == 3
      subroutine mpi_exchange8_3_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                            B, nmaxB, C, nmaxC)
# elif NARGS == 4
      subroutine mpi_exchange8_4_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                  B, nmaxB, C, nmaxC, D, nmaxD)
# endif

      use param
      use hidden_mpi_vars
      use mess_buffers
      use mpi

      implicit none

      integer istr,iend,jstr,jend, nmaxA
      real A(GLOBAL_2D_ARRAY,nmaxA)
# if NARGS > 1
      integer nmaxB, offset
      real B(GLOBAL_2D_ARRAY,nmaxB)
# endif
# if NARGS > 2
      integer nmaxC
      real C(GLOBAL_2D_ARRAY,nmaxC)
# endif
# if NARGS > 3
      integer nmaxD
      real D(GLOBAL_2D_ARRAY,nmaxD)
# endif
      integer i,j,k, kshft, ierr, mess_count, comm(16), req(16),
     &                                  status(MPI_STATUS_SIZE)
# ifdef ALT_SEND_ORD
      integer ipass
# endif
      logical flag

# include "compute_message_bounds.h"

      ! DevinD (2021/10):
      ! can't use MPI_Waitall for non-blocking send receives because
      ! domain edges don't send/receive in all directions.
      ! Therefore Sacha using his own method, which uses MPI_test
      ! and 'comm' variable to confirm all non-blocking communications are complete.
      !
      ! From internet found MPI 3.0 standard has some grid halo exchange functions
      ! already, whereby it can create a Cartesian mpi sub-domain grid automatically.
      ! This would be good for ROMS but a pain to change over.
      ! Also, didn't find function for boundaries only (i.e. not all neighbours),
      ! but sure it exists.


! Prepare to receive everything:  Permutation array comm(1:16) keeps
!-------- -- ------- -----------  track which messages are actually
! being received -- hence comm(indx)=0  means that no messages are
! expected the the direction labelled "indx", while for active messages
! "comm" keeps index of the corresponding requist handle "req".
!  This is needed because later in this code array "req" is subject to
! rearrangement in order to ignore directions from which no message is
! expected, as well as to ignore requests from which messages are
! already received and unpacked

c--#define CR
                     ! tags for receive      for send        ! each sub-domain can receive up to 8 exchanges.
      do i=1,16      !         3  5  1        4  6  2        ! the 'tag' value indicates which side/corner
        comm(i)=0    !         8     7        7     8        ! the message came from.
      enddo          !         2  6  4        1  5  3

      if (WEST_MSG_EXCH) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,  ! recvW is array to store incoming data
     &          p_W, jtg+8, ocean_grid_comm, req(1), ierr)   ! jsize is amount of data to receive of (up to 4 arrays - compute_message_bounds.h)
        comm(1)=1                                            ! p_W   is mpi_rank of sub-domain to the west
      endif                                                  ! jtg   is used to define tile location. No tiling jtg=0
                                                             !       +8 is the tag which indicates side/corrner message came from.
      if (EAST_MSG_EXCH) then                                ! ocean_grid_comm = MPI_COMM_WORLD is the group of all mpi procs in ROMS run
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,  ! req() stores the mpi request handle to later use
     &          p_E, jtg+7, ocean_grid_comm, req(2), ierr)   !       mpi_test to see if it has been completed.
        comm(2)=2
      endif

      if (SOUTH_MSG_EXCH) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+6, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif

      if (NORTH_MSG_EXCH) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+5, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif

! corners ...

      if (WEST_MSG_EXCH .and. SOUTH_MSG_EXCH) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &            p_SW, 2, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif

      if (EAST_MSG_EXCH .and. NORTH_MSG_EXCH) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &            p_NE, 1, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif

      if (EAST_MSG_EXCH .and. SOUTH_MSG_EXCH) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &            p_SE, 4, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif

      if (WEST_MSG_EXCH .and. NORTH_MSG_EXCH) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &            p_NW, 3, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif


! Send everything, sides...
!----- --------------------

# ifdef ALT_SEND_ORD
      do ipass=0,1
        if (mod(inode+ipass,2)==0) then
# endif
          if (WEST_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange: 1.1', mynode
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
# if NARGS > 1
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=B(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=B(iwest+1,j,k)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=C(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=C(iwest+1,j,k)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=D(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=D(iwest+1,j,k)
              enddo
            enddo
# endif
# ifdef NONBLK_SEND
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+7, ocean_grid_comm, req(9), ierr)

            comm(9)=9
# else
            call MPI_Send  (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+7, ocean_grid_comm,         ierr)
# endif
          endif
# ifdef ALT_SEND_ORD
        else
# endif

          if (EAST_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange: 1.2', mynode
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
# if NARGS > 1
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=B(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=B(ieast  ,j,k)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=C(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=C(ieast  ,j,k)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=D(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=D(ieast  ,j,k)
              enddo
            enddo
# endif
# ifdef NONBLK_SEND
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+8, ocean_grid_comm, req(10), ierr)

            comm(10)=10
# else
            call MPI_Send  (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+8, ocean_grid_comm,          ierr)
# endif
          endif
# ifdef ALT_SEND_ORD
        endif
      enddo

      do ipass=0,1
        if (mod(jnode+ipass,2)==0) then
# endif
          if (SOUTH_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange: 1.3', mynode
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
# if NARGS > 1
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=B(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=B(i,jsouth+1,k)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=C(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=C(i,jsouth+1,k)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=D(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=D(i,jsouth+1,k)
              enddo
            enddo
# endif
# ifdef NONBLK_SEND
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+5, ocean_grid_comm, req(11), ierr)

            comm(11)=11
# else
            call MPI_Send  (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+5, ocean_grid_comm,          ierr)
# endif
          endif

# ifdef ALT_SEND_ORD
        else
# endif
          if (NORTH_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange: 1.4', mynode
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
# if NARGS > 1
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=B(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=B(i,jnorth  ,k)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=C(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=C(i,jnorth  ,k)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=D(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=D(i,jnorth  ,k)
              enddo
            enddo
# endif
# ifdef NONBLK_SEND
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+6, ocean_grid_comm, req(12), ierr)

            comm(12)=12
# else
            call MPI_Send  (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+6, ocean_grid_comm,          ierr)
# endif
          endif
# ifdef ALT_SEND_ORD
        endif
      enddo
# endif

! ...and corners:

        if (WEST_MSG_EXCH .and. SOUTH_MSG_EXCH) then
CR        write(*,*) 'mpi_exchange: 1.5', mynode
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
# if NARGS > 1
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SW(k        +offset)=B(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxB+offset)=B(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxB+offset)=B(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxB+offset)=B(iwest+1,jsouth+1,k)
          enddo
# endif
# if NARGS > 2
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SW(k        +offset)=C(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxC+offset)=C(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxC+offset)=C(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxC+offset)=C(iwest+1,jsouth+1,k)
          enddo
# endif
# if NARGS > 3
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_SW(k        +offset)=D(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxD+offset)=D(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxD+offset)=D(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxD+offset)=D(iwest+1,jsouth+1,k)
          enddo
# endif
# ifdef NONBLK_SEND
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 1, ocean_grid_comm, req(13), ierr)

          comm(13)=13
# else
          call MPI_Send  (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 1, ocean_grid_comm,          ierr)
# endif
        endif

        if (EAST_MSG_EXCH .and. NORTH_MSG_EXCH) then
CR        write(*,*) 'mpi_exchange: 1.6', mynode
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
# if NARGS > 1
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NE(k        +offset)=B(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxB+offset)=B(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxB+offset)=B(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxB+offset)=B(ieast  ,jnorth  ,k)
          enddo
# endif
# if NARGS > 2
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NE(k        +offset)=C(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxC+offset)=C(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxC+offset)=C(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxC+offset)=C(ieast  ,jnorth  ,k)
          enddo
# endif
# if NARGS > 3
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_NE(k        +offset)=D(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxD+offset)=D(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxD+offset)=D(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxD+offset)=D(ieast  ,jnorth  ,k)
          enddo
# endif
# ifdef NONBLK_SEND
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm, req(14), ierr)

          comm(14)=14
# else
          call MPI_Send  (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 2, ocean_grid_comm,          ierr)
# endif
        endif

        if (EAST_MSG_EXCH .and. SOUTH_MSG_EXCH) then
CR        write(*,*) 'mpi_exchange: 1.7', mynode
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
# if NARGS > 1
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SE(k        +offset)=B(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxB+offset)=B(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxB+offset)=B(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxB+offset)=B(ieast  ,jsouth+1,k)
          enddo
# endif
# if NARGS > 2
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SE(k        +offset)=C(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxC+offset)=C(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxC+offset)=C(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxC+offset)=C(ieast  ,jsouth+1,k)
          enddo
# endif
# if NARGS > 3
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_SE(k        +offset)=D(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxD+offset)=D(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxD+offset)=D(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxD+offset)=D(ieast  ,jsouth+1,k)
          enddo
# endif
# ifdef NONBLK_SEND
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm, req(15), ierr)

          comm(15)=15
# else
          call MPI_Send  (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 3, ocean_grid_comm,          ierr)
# endif
        endif

        if (WEST_MSG_EXCH .and. NORTH_MSG_EXCH) then
CR        write(*,*) 'mpi_exchange: 1.8', mynode
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
# if NARGS > 1
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NW(k        +offset)=B(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxB+offset)=B(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxB+offset)=B(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxB+offset)=B(iwest+1,jnorth  ,k)
          enddo
# endif
# if NARGS > 2
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NW(k        +offset)=C(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxC+offset)=C(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxC+offset)=C(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxC+offset)=C(iwest+1,jnorth  ,k)
          enddo
# endif
# if NARGS > 3
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_NW(k        +offset)=D(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxD+offset)=D(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxD+offset)=D(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxD+offset)=D(iwest+1,jnorth  ,k)
          enddo
# endif
# ifdef NONBLK_SEND
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm, req(16), ierr)

          comm(16)=16
# else
          call MPI_Send  (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 4, ocean_grid_comm,          ierr)
# endif
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

            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                A(iwest-2,j,k)=recvW(j-jmin       +kshft)
                A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
# if NARGS > 1
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(iwest-2,j,k)=recvW(j-jmin       +kshft)
                B(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                C(iwest-2,j,k)=recvW(j-jmin       +kshft)
                C(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                D(iwest-2,j,k)=recvW(j-jmin       +kshft)
                D(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
              enddo
            enddo
# endif
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
# if NARGS > 1
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                B(ieast+1,j,k)=recvE(j-jmin       +kshft)
                B(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                C(ieast+1,j,k)=recvE(j-jmin       +kshft)
                C(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                D(ieast+1,j,k)=recvE(j-jmin       +kshft)
                D(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
# endif
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
# if NARGS > 1
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jsouth-2,k)=recvS(i-imin       +kshft)
                B(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                C(i,jsouth-2,k)=recvS(i-imin       +kshft)
                C(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                D(i,jsouth-2,k)=recvS(i-imin       +kshft)
                D(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
# endif
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
# if NARGS > 1
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                B(i,jnorth+1,k)=recvN(i-imin       +kshft)
                B(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
# endif
# if NARGS > 2
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                C(i,jnorth+1,k)=recvN(i-imin       +kshft)
                C(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
# endif
# if NARGS > 3
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                D(i,jnorth+1,k)=recvN(i-imin       +kshft)
                D(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
# endif
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
# if NARGS > 1
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              B(iwest-1,jsouth-2,k)=rv_SW(k + nmaxB+offset)
              B(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxB+offset)
              B(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxB+offset)
            enddo
# endif
# if NARGS > 2
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              C(iwest-1,jsouth-2,k)=rv_SW(k + nmaxC+offset)
              C(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxC+offset)
              C(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxC+offset)
            enddo
# endif
# if NARGS > 3
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
              D(iwest-1,jsouth-2,k)=rv_SW(k + nmaxD+offset)
              D(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxD+offset)
              D(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxD+offset)
            enddo
# endif
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
# if NARGS > 1
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              B(ieast+2,jnorth+1,k)=rv_NE(k + nmaxB+offset)
              B(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxB+offset)
              B(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxB+offset)
            enddo
# endif
# if NARGS > 2
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              C(ieast+2,jnorth+1,k)=rv_NE(k + nmaxC+offset)
              C(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxC+offset)
              C(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxC+offset)
            enddo
# endif
# if NARGS > 3
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
              D(ieast+2,jnorth+1,k)=rv_NE(k + nmaxD+offset)
              D(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxD+offset)
              D(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxD+offset)
            enddo
# endif
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
# if NARGS > 1
            offset=4*nmaxA
            do k=1,nmaxB
              B(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              B(ieast+2,jsouth-2,k)=rv_SE(k + nmaxB+offset)
              B(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxB+offset)
              B(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxB+offset)
            enddo
# endif
# if NARGS > 2
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              C(ieast+2,jsouth-2,k)=rv_SE(k + nmaxC+offset)
              C(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxC+offset)
              C(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxC+offset)
            enddo
# endif
# if NARGS > 3
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
              D(ieast+2,jsouth-2,k)=rv_SE(k + nmaxD+offset)
              D(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxD+offset)
              D(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxD+offset)
            enddo
# endif
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
# if NARGS > 1
            offset=4*nmaxA
            do k=1,nmaxB
              B(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              B(iwest-1,jnorth+1,k)=rv_NW(k + nmaxB+offset)
              B(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxB+offset)
              B(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxB+offset)
            enddo
# endif
# if NARGS > 2
            offset=offset + 4*nmaxB
            do k=1,nmaxC
              C(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              C(iwest-1,jnorth+1,k)=rv_NW(k + nmaxC+offset)
              C(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxC+offset)
              C(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxC+offset)
            enddo
# endif
# if NARGS > 3
            offset=offset + 4*nmaxC
            do k=1,nmaxD
              D(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
              D(iwest-1,jnorth+1,k)=rv_NW(k + nmaxD+offset)
              D(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxD+offset)
              D(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxD+offset)
            enddo
# endif
          endif
        endif

# ifdef NONBLK_SEND
        do j=9,16
          if (comm(j) > 0) then
            call MPI_Test (req(j), flag, status, ierr)         ! make sure all sent messages have been received elsewhere
            if (flag) then                                     ! before moving on
               mess_count=mess_count-1 ; comm(j)=0
            endif
          endif
        enddo
# endif
      enddo !<--  while (mess_count>0)
      end

# ifndef NARGS
#  define NARGS 2
# elif NARGS == 2
#  undef NARGS
#  define NARGS 3
# elif NARGS == 3
#  undef NARGS
#  define NARGS 4
# elif NARGS == 4
#  undef NARGS
# endif
# ifdef NARGS
#  include "mpi_exchange8TA.F"
# endif
#else
      subroutine mpi_exchange8_empty
      end
#endif  /* MPI */
