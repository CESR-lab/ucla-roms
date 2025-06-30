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
      integer :: size_XI, size_ETA
#else
      integer :: size_XI,
     &           size_ETA,
#endif
     &         sse, ssz,
     &         N2d,
     &         N3d


      real,allocatable,dimension(:,:) :: A2d
#ifdef SOLVE3D
      real,allocatable,dimension(:,:) :: A3d
      integer,allocatable,dimension(:,:) :: iA2d
#endif
C$OMP THREADPRIVATE( A3d, A2d, iA2d, A2d )

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_private_scratch  ![
      use scalars
      implicit none

#ifdef ALLOW_SINGLE_BLOCK_MODE
      size_XI=6+Lm; size_ETA=6+Mm
#else
      size_XI=7+(Lm+NSUB_X-1)/NSUB_X
      size_ETA=7+(Mm+NSUB_E-1)/NSUB_E
#endif
      sse=size_ETA/(N+1);  ssz=(N+1)/size_ETA
      N2d=size_XI*(sse*size_ETA+ssz*(N+1))/(sse+ssz)
      N3d=size_XI*size_ETA*(N+1)

      ! Initialize (first touch) private
      ! scratch arrays in parallel by each thread.
      allocate( A2d(N2d,32) ); A2d=0.
#ifdef SOLVE3D
      allocate( A3d(N3d,6)  ); A3d=0.            ! if private scratch not set to zero result changes
      allocate( iA2d(N2d,2) ); iA2d=0.           ! this hasn't been followed to see why
#endif

      end subroutine init_arrays_private_scratch  !]

      !----------------------------------------------------------------

      end module private_scratch
