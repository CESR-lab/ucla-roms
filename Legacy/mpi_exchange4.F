/* This package contains MPI halo exchange subroutines designed to
 communicate 2-point-wide computational margins of from 1 to up to 4
 arrays at once. All four versions are generated from common code via
 self-expansion controlled by multiply defined CPP-macro NARGS.
 Macro ALT_SEND_ORD switches between using strictly paired
 Send-Recv -- Recv-Send vs. using MPI_SendRecv.
*/

#ifndef NARGS
# include "cppdefs.h"
# define ALT_SEND_ORD
#endif
#ifdef MPI
# ifndef NARGS
      subroutine mpi_exchange4_tile (istr,iend,jstr,jend, A, nmaxA)
# elif NARGS == 2
      subroutine mpi_exchange4_2_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                                      B, nmaxB)
# elif NARGS == 3
      subroutine mpi_exchange4_3_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                            B, nmaxB, C, nmaxC)
# elif NARGS == 4
      subroutine mpi_exchange4_4_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                  B, nmaxB, C, nmaxC, D, nmaxD)
# endif
      implicit none
# include "mpif.h"
# include "param.h"
# include "mess_buffers.h"
# include "hidden_mpi_vars.h"
      integer istr,iend,jstr,jend, nmaxA
      real A(GLOBAL_2D_ARRAY,nmaxA)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if NARGS > 1
      integer nmaxB, offset
      real B(GLOBAL_2D_ARRAY,nmaxB)
CSDISTRIBUTE_RESHAPE B(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if NARGS > 2
      integer nmaxC
      real C(GLOBAL_2D_ARRAY,nmaxC)
CSDISTRIBUTE_RESHAPE C(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if NARGS > 3
      integer nmaxD
      real D(GLOBAL_2D_ARRAY,nmaxD)
CSDISTRIBUTE_RESHAPE D(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
      integer i,j,k, kshft, md_XI,md_ETA, ipass, ierr,
     &        stts1(MPI_STATUS_SIZE),  stts2(MPI_STATUS_SIZE),
     &        stts3(MPI_STATUS_SIZE),  stts4(MPI_STATUS_SIZE)

# define FOUR_MESSAGES
# include "compute_message_bounds.h"

! This version of mpi_exchange uses two-stage algorithm which first
! updates north and south computational margins of the subdomain by
! sending messages which contain only internal points and physical
! boundaries computed by the sending node.  After completion of this
! stage, it sends east- and west-bound messages containing not only
! internal and physical boundary points, but also the two ghost
! points on each end which belong to north and south computational
! margins filled-in during the first stage.  As the result, there is
! no need to send separate messages in the diagonal direction to fill
! corners, however, corner data travels in two stages, so that the
! second stage must start only after the first is complete.

! During the both north-south and east-west stages the messages are
! paired in Send-Recv and Recv-Send arrangements to avoid head-on
! collisions (a principle known as "MPI deadlock safety" --- in
! principle sends and receives can be logically synchronous and
! buffer-less in this code).

! In the code below loop "ipass" is merely to arrange pairing of
! Send-and-Recv in such a way that if one subdomain sends message to
! his neighbor (say on the on the south), the neighbor is receiving
! this message first (i.e. message coming from his north), rather than
! than send his south-bound message, as it is done by first subdomain.
! Note that even-numbered jnodes use Send-Recv sequence, while the
! order is reversed for odd-numbered jnodes.    As the result, at any
! north-south contact between two MPI nodes, the south-bound message
! is transmitted and received first; the north-bound follows.
! Similar pairing takes place on east-west sides.

! In addition to that, both north-south and east-west comminications
! are directionally alternated to reduce interference on hardware nodes
! with multiple CPUs. This is done via inner control logic "md_XI==0"
! or "md_ETA==0" with the rationale is as follows: suppose a pair of
! MPI processes is working on subdomains adjacent in ETA-direction
! belong to the same dual-CPU hardware node and share the same
! communication device.  Then messsages are arranged in such a way,
! that when one MPI process from the pair sends, the other one receives
! messages to/from members of a similar pair residing on neigboring
! hardware node.  The polarity reverses during the second stage
! "ipass", once "md_ETA" changes its parity.  This resilts in balanced
! full-duplex communication between the hardware nodes during both
! "ipass" stages.  Similar arrangement takes place for north-south
! communication if subdomains residing on the same hardware node are
! adjacent in XI-direction.  This optimization has neutral effect
! (no advantage, no penalty) for single-CPU nodes.

! Use synchroneous (Ssend) version of MPI Send, it makes it faster.

c--#define CR
#define MPI_Send MPI_Ssend

      do ipass=0,1
        md_XI=mod(inode+ipass,2)
        md_ETA=mod(jnode+ipass,2)
        if (md_ETA==0) then
          if (SOUTH_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange4: 1.3', mynode
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

# ifdef ALT_SEND_ORD
            if (md_XI==0) then
              call MPI_Send (sendS, isize, MPI_DOUBLE_PRECISION,
     &                 p_S, itg+5, ocean_grid_comm,        ierr)

              call MPI_Recv (recvS, isize, MPI_DOUBLE_PRECISION,
     &                 p_S, itg+6, ocean_grid_comm, stts2, ierr)
            else
              call MPI_Recv (recvS, isize, MPI_DOUBLE_PRECISION,
     &                 p_S, itg+6, ocean_grid_comm, stts2, ierr)

              call MPI_Send (sendS, isize, MPI_DOUBLE_PRECISION,
     &                 p_S, itg+5, ocean_grid_comm,        ierr)
            endif
# else
            call MPI_SendRecv (sendS, isize, MPI_DOUBLE_PRECISION,
     &                 p_S, itg+5,
     &                         recvS, isize, MPI_DOUBLE_PRECISION,
     &                 p_S, itg+6,
     &                              ocean_grid_comm, stts2, ierr)
# endif
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
        else
          if (NORTH_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange4: 1.4', mynode
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

# ifdef ALT_SEND_ORD
            if (md_XI==0) then
              call MPI_Recv (recvN, isize, MPI_DOUBLE_PRECISION,
     &                 p_N, itg+5, ocean_grid_comm, stts1, ierr)

              call MPI_Send (sendN, isize, MPI_DOUBLE_PRECISION,
     &                 p_N, itg+6, ocean_grid_comm,        ierr)
            else
              call MPI_Send (sendN, isize, MPI_DOUBLE_PRECISION,
     &                 p_N, itg+6, ocean_grid_comm,        ierr)

              call MPI_Recv (recvN, isize, MPI_DOUBLE_PRECISION,
     &                 p_N, itg+5, ocean_grid_comm, stts1, ierr)
            endif
# else
            call MPI_SendRecv (sendN, isize, MPI_DOUBLE_PRECISION,
     &                 p_N, itg+6,
     &                         recvN, isize, MPI_DOUBLE_PRECISION,
     &                 p_N, itg+5,
     &                              ocean_grid_comm, stts1, ierr)
# endif

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
      enddo  !<-- ipass


      do ipass=0,1
        md_XI=mod(inode+ipass,2)
        md_ETA=mod(jnode+ipass,2)
        if (md_XI==0) then
          if (WEST_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange4: 1.1', mynode
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

# ifdef ALT_SEND_ORD
            if (md_ETA==0) then
              call MPI_Send (sendW, jsize, MPI_DOUBLE_PRECISION,
     &                 p_W, jtg+7, ocean_grid_comm,        ierr)

              call MPI_Recv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &                 p_W, jtg+8, ocean_grid_comm,stts4,ierr)
            else
              call MPI_Recv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &                 p_W, jtg+8, ocean_grid_comm,stts4,ierr)

              call MPI_Send (sendW, jsize, MPI_DOUBLE_PRECISION,
     &                 p_W, jtg+7, ocean_grid_comm,        ierr)
            endif
# else
            call MPI_SendRecv (sendW, jsize, MPI_DOUBLE_PRECISION,
     &                 p_W, jtg+7,
     &                         recvW, jsize, MPI_DOUBLE_PRECISION,
     &                 p_W, jtg+8,
     &                              ocean_grid_comm, stts4, ierr)
# endif

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
        else
          if (EAST_MSG_EXCH) then
CR          write(*,*) 'mpi_exchange4: 1.2', mynode
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

# ifdef ALT_SEND_ORD
            if (md_ETA==0) then
              call MPI_Recv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &                 p_E, jtg+7, ocean_grid_comm, stts3, ierr)

              call MPI_Send (sendE, jsize, MPI_DOUBLE_PRECISION,
     &                 p_E, jtg+8, ocean_grid_comm,        ierr)
            else
              call MPI_Send (sendE, jsize, MPI_DOUBLE_PRECISION,
     &                 p_E, jtg+8, ocean_grid_comm,        ierr)

              call MPI_Recv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &                 p_E, jtg+7, ocean_grid_comm, stts3, ierr)
            endif
# else
            call MPI_SendRecv (sendE, jsize, MPI_DOUBLE_PRECISION,
     &                 p_E, jtg+8,
     &                         recvE, jsize, MPI_DOUBLE_PRECISION,
     &                 p_E, jtg+7,
     &                              ocean_grid_comm, stts3, ierr)
# endif

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
      enddo

CR    write(*,*) 'leaving mpi_exchange4 ', mynode
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
#  include "mpi_exchange4.F"
# endif
#else
      subroutine mpi_exchange_empty
      end
#endif
