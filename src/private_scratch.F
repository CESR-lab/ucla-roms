      module private_scratch

#include "cppdefs.opt"

      use param
      implicit none

! Auxiliary module "private_scratch":   A set of tile-size arrays
!---------- ------ --------------------   to provide workspace for
! intermediate computations individually for each thread.  The arrays 
! passed as arguments to physical routines called by their drivers
! and are used internally there.  In most cases the result is no
! longer needed upon completion of physical routines; occasionally
! there arrays are used to transmit data between physical routines
! working on the same tile within the same parallel region.
!
! Note that "sse-ssz" below are to make N2d=size_XI*max(size_ETA,N+1)
! without using "max" function inside parameter statement.
  
#ifdef ALLOW_SINGLE_BLOCK_MODE
      integer, parameter :: size_XI=6+Lm, size_ETA=6+Mm,
#else
      integer, parameter :: size_XI=7+(Lm+NSUB_X-1)/NSUB_X,
     &                      size_ETA=7+(Mm+NSUB_E-1)/NSUB_E,
#endif
     &         sse=size_ETA/(N+1),  ssz=(N+1)/size_ETA,
     &         N2d=size_XI*(sse*size_ETA+ssz*(N+1))/(sse+ssz),
     &         N3d=size_XI*size_ETA*(N+1)

      real A2d(N2d,32)
#ifdef SOLVE3D
      real A3d(N3d,6)
      integer iA2d(N2d,2)
#endif
C$OMP THREADPRIVATE( A3d, A2d, iA2d, A2d )

      end module private_scratch
