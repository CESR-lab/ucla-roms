      module mess_buffers

#include "cppdefs.opt"

      use param
      implicit none

! Note: these buffers are designed to accommodate MPI computational
! margins for up to 4 full 3D arrays to be exchanged simultaneously.

      integer :: size_Z, size_X, size_E

      real,allocatable,dimension(:) :: sn_NW, sendN, sn_NE, rv_NW, recvN, rv_NE,
     &                                 sendW, sendE, recvW, recvE, sn_SW, sendS,
     &                                 sn_SE, rv_SW, recvS, rv_SE

C$OMP THREADPRIVATE(sn_NW,rv_NW, sendN,recvN, sn_NE,rv_NE)
C$OMP THREADPRIVATE(sendW,recvW, sendE,recvE)
C$OMP THREADPRIVATE(sn_SW,rv_SW, sendS,recvS, sn_SE,rv_SE)

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_mess_buffers  ![
      use scalars
      implicit none

      size_Z=16*(N+1); size_X=8*(N+1)*(Lm+4)   ! DevinD 8=2x4= 2 for margin, up to 4 different arrays
      size_E=8*(N+1)*(Mm+4)                    ! size_Z = 4 in corner, up to 4 different arrays

      allocate( sn_NW(size_Z),   sendN(size_X),   sn_NE(size_Z),
     &          rv_NW(size_Z),   recvN(size_X),   rv_NE(size_Z),
     &          sendW(size_E),                    sendE(size_E),
     &          recvW(size_E),                    recvE(size_E),
     &          sn_SW(size_Z),   sendS(size_X),   sn_SE(size_Z),
     &          rv_SW(size_Z),   recvS(size_X),   rv_SE(size_Z)  )

      sn_NW=init ;   sendN=init ;   sn_NE=init
      rv_NW=init ;   recvN=init ;   rv_NE=init
      sendW=init ;                  sendE=init
      recvW=init ;                  recvE=init
      sn_SW=init ;   sendS=init ;   sn_SE=init
      rv_SW=init ;   recvS=init ;   rv_SE=init

      end subroutine init_arrays_mess_buffers  !]

!----------------------------------------------------------------------

      end module mess_buffers
