      module work_mod  ! was work .h but changed to work_mod for name clash with work array

      use dimensions
      implicit none

! This is "work": declaration of utility work array.
!
      real,allocatable,dimension(:,:,:) :: work

      real,allocatable,dimension(:),  target :: wrk_x
      real,allocatable,dimension(:,:),target :: wrk_xy
      real,allocatable,dimension(:,:),target :: wrk_xz

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_work_mod  ![
      use scalars
      implicit none

      allocate( work(1-bf:nx+bf,1-bf:ny+bf,0:nz) );

      allocate( wrk_x(1-bf:nx+bf) ); 
      allocate( wrk_xy(1-bf:nx+bf,1-bf:ny+bf) );
      allocate( wrk_xz(1-bf:nx+bf,0:nz) );

      end subroutine init_arrays_work_mod  !]
      !----------------------------------------------------------------

      end module work_mod
 
