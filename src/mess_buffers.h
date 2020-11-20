! Note: these buffers are designed to accommodate MPI computational
! margins for up to 4 full 3D arrays to be exchanged simultaneously.

      integer, parameter :: size_Z=16*(N+1), size_X=8*(N+1)*(Lm+4),
     &                                       size_E=8*(N+1)*(Mm+4)

      real sn_NW(size_Z),   sendN(size_X),   sn_NE(size_Z),
     &     rv_NW(size_Z),   recvN(size_X),   rv_NE(size_Z),
     &     sendW(size_E),                    sendE(size_E),
     &     recvW(size_E),                    recvE(size_E),
     &     sn_SW(size_Z),   sendS(size_X),   sn_SE(size_Z),
     &     rv_SW(size_Z),   recvS(size_X),   rv_SE(size_Z)

      common /mess_buffers/ sn_NW,rv_NW, sendN,recvN, sn_NE,rv_NE,
     &                    sendW,recvW,                sendE,recvE,
     &                     sn_SW,rv_SW, sendS,recvS, sn_SE,rv_SE
C$OMP THREADPRIVATE(/mess_buffers/)
