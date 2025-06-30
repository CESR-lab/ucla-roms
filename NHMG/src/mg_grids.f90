module mg_grids

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist

  implicit none

  !-
  !- Derived type are created to place all grid information in one structure.
  !- The array memory allocation is called only one time for each grid level
  !- at the begining of the program. This avoid memory allocation and deallocation
  !- during the execution of the program.
  !- Pointers are used to simplify the readability of the code.
  !- In subroutines which use grids arrays, we use this tip to improve readability:
  !-     real(kind=rp),dimension(:,:,:), pointer:: p
  !-     p  => grid(lev)%p
  !-

  !-
  !- GRID_TYPE is a derived type which has the function
  !- to reserve memory space for all the main variables of the multigrid solver.
  !-
  type grid_type

     integer(kind=ip) :: nx,ny,nz                     ! Subdomain dimensions (whithout halo points)
     integer(kind=ip) :: npx, npy                     ! Number of processes
     integer(kind=ip) :: incx, incy                   ! ??
     integer(kind=ip) :: Ng2D, Ng ,Ngp                ! Size for 2D and 3D arrays (included halo=1)
     integer(kind=ip) :: gather                       ! Gathering activated (1) or not (0)
     integer(kind=ip) :: ngx, ngy                     ! Gathering in x and y (is 1 or 2 (and generally 2))
     integer(kind=ip) :: key                          ! Gathering key
     integer(kind=ip) :: localcomm                    ! Gathering MPI comm (should be integer (output of MPI_SPLIT))
     integer(kind=ip),dimension(8)::neighb            ! MPI: neighbours

     real(kind=rp),dimension(:,:,:,:),pointer :: cA => null() ! Matrix of coefficients

     real(kind=rp),dimension(:,:,:),pointer :: p => null()    ! Pressure
     real(kind=rp),dimension(:,:,:),pointer :: b => null()    ! Right-hand side
     real(kind=rp),dimension(:,:,:),pointer :: r => null()    ! Residual

     real(kind=rp),dimension(:,:,:),pointer :: px => null()    ! Pressure
     real(kind=rp),dimension(:,:,:),pointer :: py => null()    ! Right-hand side
     real(kind=rp),dimension(:,:,:),pointer :: pz => null()    ! Residual

     real(kind=rp),dimension(:,:)  ,pointer :: dx => null()   ! Mesh in x  (1 halo point)
     real(kind=rp),dimension(:,:)  ,pointer :: dy => null()   ! Mesh in y  (1 halo point)
     real(kind=rp),dimension(:,:)  ,pointer :: dxu => null()  ! Mesh in x  (1 halo point)
     real(kind=rp),dimension(:,:)  ,pointer :: dyv => null()  ! Mesh in y  (1 halo point)
     real(kind=rp),dimension(:,:,:),pointer :: dz => null()   ! Mesh in z at w point   (nz  , 2 halo points)

     ! All these variables are dependent on dx, dy, and dz
     real(kind=rp),dimension(:,:,:),pointer :: dzw   => null() ! Cell height at w-points
     real(kind=rp),dimension(:,:,:),pointer :: Arx   => null() ! Cell face surface at u-points
     real(kind=rp),dimension(:,:,:),pointer :: Ary   => null() ! Cell face surface at v-points
     real(kind=rp),dimension(:,:)  ,pointer :: Arz   => null() !
     real(kind=rp),dimension(:,:)  ,pointer :: beta  => null() !
     real(kind=rp),dimension(:,:)  ,pointer :: gamu  => null() !
     real(kind=rp),dimension(:,:)  ,pointer :: gamv  => null() !
     real(kind=rp),dimension(:,:,:),pointer :: zxdy  => null() ! Slopes in x-direction defined at rho-points * dy
     real(kind=rp),dimension(:,:,:),pointer :: zydx  => null() ! Slopes in y-direction defined at rho-points * dx
     real(kind=rp),dimension(:,:,:),pointer :: alpha => null() ! All levels

     ! Dummy array to calculate uf, vf and wf
     ! Remark: We can gain memory place using p instead of dummy3Dnz
     real(kind=rp),dimension(:,:,:),pointer :: dummy3Dnz  => null()
     real(kind=rp),dimension(:,:,:),pointer :: dummy3Dnzp  => null()

     ! Gathering buffers, allocated only if gathering is activated (coarsen grids!)
     real(kind=rp),dimension(:,:,:)    ,pointer :: dummy3 => null()         ! A dummy 3D array for gathering
     real(kind=rp),dimension(:,:,:,:)  ,pointer :: gatherbuffer2D => null() ! 2D
     real(kind=rp),dimension(:,:,:,:,:),pointer :: gatherbuffer   => null() ! 3D nz

     ! croco exchange arrays: resized and reshaped (i,j,k)=>(k,j,i)
     real(kind=rp),dimension(:,:,:)  ,pointer :: du => null()
     real(kind=rp),dimension(:,:,:)  ,pointer :: dv => null()
     real(kind=rp),dimension(:,:,:)  ,pointer :: dw => null()

  end type grid_type

  !-
  !- MPI_BUFFERS is a derived type which has the function
  !- to reserve memory for MPI exchanges at the boundaries
  !-
  type mpi_buffers

     ! MPI: 2D array buffers (halo=1 and halo=2 => h, dz)
     real(kind=rp),dimension(:,:),pointer :: sendN2D1,recvN2D1,sendS2D1,recvS2D1
     real(kind=rp),dimension(:,:),pointer :: sendE2D1,recvE2D1,sendW2D1,recvW2D1
     real(kind=rp),dimension(:,:),pointer :: sendSW2D1,recvSW2D1,sendSE2D1,recvSE2D1
     real(kind=rp),dimension(:,:),pointer :: sendNW2D1,recvNW2D1,sendNE2D1,recvNE2D1

     ! MPI: 3D array buffers (nz and nz+1) (halo=1)
     real(kind=rp),dimension(:,:,:),pointer :: sendN,recvN,sendS,recvS
     real(kind=rp),dimension(:,:,:),pointer :: sendE,recvE,sendW,recvW
     real(kind=rp),dimension(:,:,:),pointer :: sendSW,recvSW,sendSE,recvSE
     real(kind=rp),dimension(:,:,:),pointer :: sendNW,recvNW,sendNE,recvNE
     real(kind=rp),dimension(:,:,:),pointer :: sendNp,recvNp,sendSp,recvSp
     real(kind=rp),dimension(:,:,:),pointer :: sendEp,recvEp,sendWp,recvWp
     real(kind=rp),dimension(:,:,:),pointer :: sendSWp,recvSWp,sendSEp,recvSEp
     real(kind=rp),dimension(:,:,:),pointer :: sendNWp,recvNWp,sendNEp,recvNEp

     ! MPI: 3D array buffers (nz and nz+1) (halo=2)
     real(kind=rp),dimension(:,:,:),pointer :: sendN3D2,recvN3D2,sendS3D2,recvS3D2
     real(kind=rp),dimension(:,:,:),pointer :: sendE3D2,recvE3D2,sendW3D2,recvW3D2
     real(kind=rp),dimension(:,:,:),pointer :: sendSW3D2,recvSW3D2,sendSE3D2,recvSE3D2
     real(kind=rp),dimension(:,:,:),pointer :: sendNW3D2,recvNW3D2,sendNE3D2,recvNE3D2
     real(kind=rp),dimension(:,:,:),pointer :: sendN3D2p,recvN3D2p,sendS3D2p,recvS3D2p
     real(kind=rp),dimension(:,:,:),pointer :: sendE3D2p,recvE3D2p,sendW3D2p,recvW3D2p
     real(kind=rp),dimension(:,:,:),pointer :: sendSW3D2p,recvSW3D2p,sendSE3D2p,recvSE3D2p
     real(kind=rp),dimension(:,:,:),pointer :: sendNW3D2p,recvNW3D2p,sendNE3D2p,recvNE3D2p

     ! MPI: 4D matrix coefficient buffers (halo=1)
     real(kind=rp),dimension(:,:,:),pointer :: sendN4D,recvN4D,sendS4D,recvS4D
     real(kind=rp),dimension(:,:,:),pointer :: sendE4D,recvE4D,sendW4D,recvW4D
     real(kind=rp),dimension(:,:)  ,pointer :: sendSW4D,recvSW4D,sendSE4D,recvSE4D
     real(kind=rp),dimension(:,:)  ,pointer :: sendNW4D,recvNW4D,sendNE4D,recvNE4D

  end type mpi_buffers


  !-
  !- GRID is an vector of grid_type (derived type) which stores the informations of
  !- all the different grid levels of the multigrid solver.
  !- It is the main array of this application.
  !- It is the most important to :-)
  !-
  type(grid_type)  , dimension(:), pointer :: grid
  type(mpi_buffers), dimension(:), pointer :: gbuffers

  ! Number of grid levels of the multigrid solver
  integer(kind=ip):: nlevs ! index of the coarsest level (1 is the finest)q

  real(kind=rp) :: nhhc, nhtheta_b, nhtheta_s

contains

  !----------------------------------------
  subroutine define_grids(npxg, npyg, nxl, nyl, nzl)

    integer(kind=ip), intent(in) :: npxg,npyg     ! Number of processes in x and y
    integer(kind=ip), intent(in) :: nxl, nyl, nzl ! Subdomain dimensions


    ! Number of matrix coefficient
    ! The matrix coefficient is symmetric and compact
    ! 2D : nd = 5
    ! 3D : nd = 8 (equivalent to 15 points = 7 + 7' + 1 where 7 and 7' are symmetric)
    integer(kind=ip) :: nd

    integer(kind=ip) :: nx, ny, nz

    integer(kind=ip) :: lev ! loop's variable on levels

    if (myrank==0) write(*,*)'  - define grids:'
    if (myrank==0) write(*,*)'    - define grid levels'

    ! Calculate the number of levels of the multigrid solver method
    call find_grid_levels(npxg, npyg, nxl, nyl, nzl)

    ! Allocate the memory for the nlevs grids
    allocate(grid(nlevs))
    allocate(gbuffers(nlevs))

    ! At lev = 1
    grid(1)%nx = nxl
    grid(1)%ny = nyl
    grid(1)%nz = nzl

    grid(1)%npx = npxg
    grid(1)%npy = npyg

    grid(1)%incx=1
    grid(1)%incy=1

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    allocate(grid(1)%px(nz,0:ny+1,0:nx+1))
    allocate(grid(1)%py(nz,0:ny+1,0:nx+1))
    allocate(grid(1)%pz(nz+1,0:ny+1,0:nx+1))

    if (myrank==0) write(*,*)'    - define grid dims'
    ! define grid dimensions at each level
    call define_grid_dims()

    !-
    !- Allocate memory for grids at each level
    !- We use a lot of loops to reduce the compilation time which can be very long
    !- if we use only one loop for these memory allocations.
    !-

    do lev=1,nlevs ! set_horiz_grids
       nx = grid(lev)%nx
       ny = grid(lev)%ny
       allocate(grid(lev)%dx(0:ny+1,0:nx+1))
       allocate(grid(lev)%dy(0:ny+1,0:nx+1))
       allocate(grid(lev)%dxu(0:ny+1,0:nx+1))
       allocate(grid(lev)%dyv(0:ny+1,0:nx+1))
    enddo

    do lev=1,nlevs ! set_vert_grids
       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       allocate(grid(lev)%dz(nz,-1:ny+2,-1:nx+2)) ! 2 extra points
    enddo

    do lev=1,nlevs ! set_vert_grids
       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       allocate(grid(lev)%dzw( nz+1,0:ny+1,0:nx+1)) ! at w point
       allocate(grid(lev)%Arx( nz  ,0:ny+1,1:nx+1)) ! at u point
       allocate(grid(lev)%Ary( nz  ,1:ny+1,0:nx+1)) ! at v point
       allocate(grid(lev)%Arz(      0:ny+1,0:nx+1)) ! at w point
       allocate(grid(lev)%beta(     0:ny+1,0:nx+1)) !
       allocate(grid(lev)%gamu(     0:ny+1,0:nx+1)) !
       allocate(grid(lev)%gamv(     0:ny+1,0:nx+1)) !
       allocate(grid(lev)%zxdy(nz  ,0:ny+1,0:nx+1)) ! at rho point
       allocate(grid(lev)%zydx(nz  ,0:ny+1,0:nx+1)) ! at rho point
       allocate(grid(lev)%alpha(nz ,0:ny+1,0:nx+1)) ! at rho point
    enddo

    lev = 1 ! Some intermediate arrays for define matrices and compute rhs
    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    allocate(grid(lev)%dummy3Dnz( nz  ,0:ny+1,0:nx+1))
    allocate(grid(lev)%dummy3Dnzp(nz+1,0:ny+1,0:nx+1))

    allocate(grid(lev)%du(1:nz  ,0:ny+1,0:nx+1))
    allocate(grid(lev)%dv(1:nz  ,0:ny+1,0:nx+1))
    allocate(grid(lev)%dw(1:nz+1,0:ny+1,0:nx+1))

    do lev=1,nlevs !
       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       allocate(grid(lev)%p(nz,0:ny+1,0:nx+1))
       allocate(grid(lev)%b(nz,0:ny+1,0:nx+1))
       allocate(grid(lev)%r(nz,0:ny+1,0:nx+1))
    enddo

    do lev=1,nlevs !
       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       if (nz == 1) then
          nd = 5
       else
          nd = 8
       endif
       allocate(grid(lev)%cA(nd,nz,0:ny+1,0:nx+1))
    enddo

    ! MPI exhanges for 2D arrays
    !-Halo 1-!
    do lev=1,nlevs
       nx = grid(lev)%nx
       allocate(gbuffers(lev)%sendS2D1(1,nx))
       allocate(gbuffers(lev)%recvS2D1(1,nx))
       allocate(gbuffers(lev)%sendN2D1(1,nx))
       allocate(gbuffers(lev)%recvN2D1(1,nx))
    enddo

    do lev=1,nlevs
       ny = grid(lev)%ny
       allocate(gbuffers(lev)%sendE2D1(ny,1))
       allocate(gbuffers(lev)%recvE2D1(ny,1))
       allocate(gbuffers(lev)%sendW2D1(ny,1))
       allocate(gbuffers(lev)%recvW2D1(ny,1))
    enddo

    do lev=1,nlevs
       allocate(gbuffers(lev)%sendSW2D1(1,1))
       allocate(gbuffers(lev)%sendSE2D1(1,1))
       allocate(gbuffers(lev)%sendNW2D1(1,1))
       allocate(gbuffers(lev)%sendNE2D1(1,1))
    enddo

    do lev=1,nlevs
       allocate(gbuffers(lev)%recvSW2D1(1,1))
       allocate(gbuffers(lev)%recvSE2D1(1,1))
       allocate(gbuffers(lev)%recvNW2D1(1,1))
       allocate(gbuffers(lev)%recvNE2D1(1,1))
    enddo

    ! MPI exhanges for 3D arrays (halo=1)
    do lev=1,nlevs
       nx = grid(lev)%nx
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendS(nz,1,nx))
       allocate(gbuffers(lev)%recvS(nz,1,nx))
       allocate(gbuffers(lev)%sendN(nz,1,nx))
       allocate(gbuffers(lev)%recvN(nz,1,nx))
    enddo

    do lev=1,nlevs
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendE(nz,ny,1))
       allocate(gbuffers(lev)%recvE(nz,ny,1))
       allocate(gbuffers(lev)%sendW(nz,ny,1))
       allocate(gbuffers(lev)%recvW(nz,ny,1))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendSW(nz,1,1))
       allocate(gbuffers(lev)%sendSE(nz,1,1))
       allocate(gbuffers(lev)%sendNW(nz,1,1))
       allocate(gbuffers(lev)%sendNE(nz,1,1))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%recvSW(nz,1,1))
       allocate(gbuffers(lev)%recvSE(nz,1,1))
       allocate(gbuffers(lev)%recvNW(nz,1,1))
       allocate(gbuffers(lev)%recvNE(nz,1,1))
    enddo

    do lev=1,nlevs
       nx = grid(lev)%nx
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendSp(nz+1,1,nx))
       allocate(gbuffers(lev)%recvSp(nz+1,1,nx))
       allocate(gbuffers(lev)%sendNp(nz+1,1,nx))
       allocate(gbuffers(lev)%recvNp(nz+1,1,nx))
    enddo

    do lev=1,nlevs
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendEp(nz+1,ny,1))
       allocate(gbuffers(lev)%recvEp(nz+1,ny,1))
       allocate(gbuffers(lev)%sendWp(nz+1,ny,1))
       allocate(gbuffers(lev)%recvWp(nz+1,ny,1))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendSWp(nz+1,1,1))
       allocate(gbuffers(lev)%sendSEp(nz+1,1,1))
       allocate(gbuffers(lev)%sendNWp(nz+1,1,1))
       allocate(gbuffers(lev)%sendNEp(nz+1,1,1))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%recvSWp(nz+1,1,1))
       allocate(gbuffers(lev)%recvSEp(nz+1,1,1))
       allocate(gbuffers(lev)%recvNWp(nz+1,1,1))
       allocate(gbuffers(lev)%recvNEp(nz+1,1,1))
    enddo

    ! MPI exhanges for 3D arrays (halo=2 ) and dz
    do lev=1,nlevs
       nx = grid(lev)%nx
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendS3D2(nz,2,nx))
       allocate(gbuffers(lev)%recvS3D2(nz,2,nx))
       allocate(gbuffers(lev)%sendN3D2(nz,2,nx))
       allocate(gbuffers(lev)%recvN3D2(nz,2,nx))
    enddo

    do lev=1,nlevs
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendE3D2(nz,ny,2))
       allocate(gbuffers(lev)%recvE3D2(nz,ny,2))
       allocate(gbuffers(lev)%sendW3D2(nz,ny,2))
       allocate(gbuffers(lev)%recvW3D2(nz,ny,2))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendSW3D2(nz,2,2))
       allocate(gbuffers(lev)%sendSE3D2(nz,2,2))
       allocate(gbuffers(lev)%sendNW3D2(nz,2,2))
       allocate(gbuffers(lev)%sendNE3D2(nz,2,2))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%recvSW3D2(nz,2,2))
       allocate(gbuffers(lev)%recvSE3D2(nz,2,2))
       allocate(gbuffers(lev)%recvNW3D2(nz,2,2))
       allocate(gbuffers(lev)%recvNE3D2(nz,2,2))
    enddo

    do lev=1,nlevs
       nx = grid(lev)%nx
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendS3D2p(nz+1,2,nx))
       allocate(gbuffers(lev)%recvS3D2p(nz+1,2,nx))
       allocate(gbuffers(lev)%sendN3D2p(nz+1,2,nx))
       allocate(gbuffers(lev)%recvN3D2p(nz+1,2,nx))
    enddo

    do lev=1,nlevs
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendE3D2p(nz+1,ny,2))
       allocate(gbuffers(lev)%recvE3D2p(nz+1,ny,2))
       allocate(gbuffers(lev)%sendW3D2p(nz+1,ny,2))
       allocate(gbuffers(lev)%recvW3D2p(nz+1,ny,2))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%sendSW3D2p(nz+1,2,2))
       allocate(gbuffers(lev)%sendSE3D2p(nz+1,2,2))
       allocate(gbuffers(lev)%sendNW3D2p(nz+1,2,2))
       allocate(gbuffers(lev)%sendNE3D2p(nz+1,2,2))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       allocate(gbuffers(lev)%recvSW3D2p(nz+1,2,2))
       allocate(gbuffers(lev)%recvSE3D2p(nz+1,2,2))
       allocate(gbuffers(lev)%recvNW3D2p(nz+1,2,2))
       allocate(gbuffers(lev)%recvNE3D2p(nz+1,2,2))
    enddo

    ! MPI exhanges for 4D CA array
    do lev=1,nlevs
       nx = grid(lev)%nx
       nz = grid(lev)%nz
       if (nz == 1) then
          nd = 5
       else
          nd = 8
       endif
       allocate(gbuffers(lev)%sendS4D(nd,nz,nx))
       allocate(gbuffers(lev)%recvS4D(nd,nz,nx))
       allocate(gbuffers(lev)%sendN4D(nd,nz,nx))
       allocate(gbuffers(lev)%recvN4D(nd,nz,nx))
    enddo

    do lev=1,nlevs
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       if (nz == 1) then
          nd = 5
       else
          nd = 8
       endif
       allocate(gbuffers(lev)%sendE4D(nd,nz,ny))
       allocate(gbuffers(lev)%recvE4D(nd,nz,ny))
       allocate(gbuffers(lev)%sendW4D(nd,nz,ny))
       allocate(gbuffers(lev)%recvW4D(nd,nz,ny))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       if (nz == 1) then
          nd = 5
       else
          nd = 8
       endif
       allocate(gbuffers(lev)%sendSW4D(nd,nz))
       allocate(gbuffers(lev)%sendSE4D(nd,nz))
       allocate(gbuffers(lev)%sendNW4D(nd,nz))
       allocate(gbuffers(lev)%sendNE4D(nd,nz))
    enddo

    do lev=1,nlevs
       nz = grid(lev)%nz
       if (nz == 1) then
          nd = 5
       else
          nd = 8
       endif
       allocate(gbuffers(lev)%recvSW4D(nd,nz))
       allocate(gbuffers(lev)%recvSE4D(nd,nz))
       allocate(gbuffers(lev)%recvNW4D(nd,nz))
       allocate(gbuffers(lev)%recvNE4D(nd,nz))
    enddo

    ! Initialize pressure to zero
    grid(1)%p(:,:,:) = zero

    if (myrank==0) write(*,*)'    - define gather informations'
    ! Call the routine which calculate gather informations
    call define_gather_informations()

  end subroutine define_grids

  !------------------------------------------------------------!
  subroutine find_grid_levels(npxg, npyg, nx,ny,nz)

    integer(kind=ip), intent(in) :: npxg, npyg ! Number of processes in x and y
    integer(kind=ip), intent(in) :: nx, ny, nz ! Subdomain dimensions (no halo points)

    integer(kind=ip) :: nxg, nyg, nzg

    integer(kind=ip) :: nhmin,nhg,nzmin, nl1,nl2

    nxg = npxg * nx ! Global domain dimension in x
    nyg = npyg * ny ! Global domain dimension in y
    nzg = nz        ! Global domain dimension in z

    nhmin = 6 ! Smallest horizontal dimension of the coarsest grid

    nzmin = 3 ! Smallest vertical dimension of the coarsest grid

    nhg = min(nxg,nyg) ! smallest horizontal dimension of the finest grid

    ! nhg = nhmin * 2^(nl1-1) and therefore:

    nl1 = 1+floor( log( nhg * one / nhmin * one) / log(two) )

    nl2 = 1+floor( log( nzg * one / nzmin * one) / log(two) )

    nlevs = min(nl1,nl2) ! The number of grid levels in the MG.

  end subroutine find_grid_levels

  !----------------------------------------
  subroutine define_grid_dims()

    integer(kind=ip) :: nx, ny, nz
    integer(kind=ip) :: npx, npy
    integer(kind=ip) :: lev, incx, incy

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    npx = grid(1)%npx
    npy = grid(1)%npy

    incx = 1
    incy = 1

    lev=1
    grid(lev)%gather=0
    grid(lev)%ngx = 1
    grid(lev)%ngy = 1

    do lev = 2, nlevs

       if (nz.eq.1) then ! 2D coarsening
          nx = nx/2
          ny = ny/2
       else              ! 3D coarsening
          nx = nx/2
          ny = ny/2
          nz = nz/2
       endif

       ! determine if gathering is needed
       !- assumes squarish nxg nyg dimensions !

       grid(lev)%gather=0
       grid(lev)%ngx = 1
       grid(lev)%ngy = 1

       if((min(nx,ny)<nsmall).and.(npx*npy>1))then
          grid(lev)%gather = 1
          if (npx > 1)then
             npx  = npx/2
             nx   = nx*2
             grid(lev)%ngx = 2
          endif
          if (npy > 1)then
             npy  = npy/2
             ny   = ny*2
             grid(lev)%ngy = 2
          endif
          incx=incx*2
          incy=incy*2

       endif

       grid(lev)%nx   = nx
       grid(lev)%ny   = ny
       grid(lev)%nz   = nz
       grid(lev)%npx  = npx
       grid(lev)%npy  = npy
       grid(lev)%incx = incx
       grid(lev)%incy = incy

    enddo

  end subroutine define_grid_dims

  !----------------------------------------
  subroutine define_neighbours(neighb)
    integer(kind=ip), dimension(4), optional, intent(in) :: neighb ! S, E, N, W

    integer(kind=ip) :: lev
    integer(kind=ip) :: npx, npy
    integer(kind=ip) :: incx, incy
    integer(kind=ip) :: pi, pj

    if (myrank==0) write(*,*)'    - define neighbours'

    npx = grid(1)%npx
    npy = grid(1)%npy

    pj = myrank/npx
    pi = mod(myrank,npx)

    ! Neighbours
    do lev=1,nlevs
       ! incx is the distance to my neighbours in x (1, 2, 4, ...)
       incx = grid(lev)%incx
       incy = grid(lev)%incy

       if (pj >= incy) then ! south
          grid(lev)%neighb(1) = (pj-incy)*npx+pi
       elseif (north_south_perio) then
          grid(lev)%neighb(1) =  mod(pj-incy+npy,npy)*npx+pi
       else
          grid(lev)%neighb(1) = MPI_PROC_NULL
       endif

       if (pi < npx-incx) then ! east
          grid(lev)%neighb(2) = pj*npx+pi+incx
       elseif (east_west_perio) then
          grid(lev)%neighb(2) =  pj*npx + mod(pi+incx,npx)
       else
          grid(lev)%neighb(2) = MPI_PROC_NULL
       endif

       if (pj < npy-incy) then ! north
          grid(lev)%neighb(3) = (pj+incy)*npx+pi
       elseif (north_south_perio) then
          grid(lev)%neighb(3) =  mod(pj+incy,npy)*npx+pi
       else
          grid(lev)%neighb(3) = MPI_PROC_NULL
       endif

       if (pi >= incx) then ! west
          grid(lev)%neighb(4) = pj*npx+pi-incx
       elseif (east_west_perio) then
          grid(lev)%neighb(4) =  pj*npx + mod(npx+pi-incx,npx)
       else
          grid(lev)%neighb(4) = MPI_PROC_NULL
       endif

       if ((pj >= incy).and.(pi >= incx)) then ! south west
          grid(lev)%neighb(5) = (pj-incy)*npx+ pi-incx
       elseif (east_west_perio.and.north_south_perio) then
          grid(lev)%neighb(5) = mod(pj-incy+npy,npy)*npx+ mod(pi-incx+npx,npx)
       else
          grid(lev)%neighb(5) = MPI_PROC_NULL
       endif

       if ((pj >= incy).and.(pi < npx-incx)) then ! south east
          grid(lev)%neighb(6) = (pj-incy)*npx+ pi+incx
       elseif (east_west_perio.and.north_south_perio) then
          grid(lev)%neighb(6) = mod(pj-incy+npy,npy)*npx+ mod(pi+incx,npx)
       else
          grid(lev)%neighb(6) = MPI_PROC_NULL
       endif

       if ((pj < npy-incy).and.(pi < npx-incx)) then ! north east
          grid(lev)%neighb(7) = (pj+incy)*npx + pi+incx
       elseif (east_west_perio.and.north_south_perio) then
          grid(lev)%neighb(7) = mod(pj+incy,npy)*npx+ mod(pi+incx,npx)
       else
          grid(lev)%neighb(7) = MPI_PROC_NULL
       endif
       if ((pj < npy-incy).and.(pi >= incx)) then ! north west
          grid(lev)%neighb(8) = (pj+incy)*npx + pi-incx
       elseif (east_west_perio.and.north_south_perio) then
          grid(lev)%neighb(8) = mod(pj+incy,npy)*npx+ mod(pi-incx+npx,npx)
       else
          grid(lev)%neighb(8) = MPI_PROC_NULL
       endif
    enddo

    if (present(neighb)) then
       ! Test for grid level 1 the consistancy with the ocean model
       If  ((grid(1)%neighb(1) /= neighb(1)).or. &
            (grid(1)%neighb(2) /= neighb(2)).or. &
            (grid(1)%neighb(3) /= neighb(3)).or. &
            (grid(1)%neighb(4) /= neighb(4))) then
          write(*,*)'Error: neighbour definition problem !'
          stop -1
       endif
    endif

  end subroutine define_neighbours

  !---------------------------------------------------------------------
  subroutine define_gather_informations()

    integer(kind=ip) :: nx, ny, nz, nd
    integer(kind=ip) :: npx, npy
    integer(kind=ip) :: incx, incy

    integer(kind=ip) :: pi, pj
    integer(kind=ip) :: lev

    ! for the gathering
    integer(kind=ip) :: ngx, ngy
    integer(kind=ip) :: N, family, nextfamily, color, key, localcomm, ierr

    ! Global indexing is always used to locate each core
    ! A core that has coordinates (2,3) on the finest decomposition
    ! will remain at this location after gathering
    npx = grid(1)%npx ! grid(1) is not a bug!
    npy = grid(1)%npy
    pj = myrank/npx
    pi = mod(myrank,npx)

    do lev=2,nlevs

       if(grid(lev)%gather.eq.1)then

          nx = grid(lev)%nx
          ny = grid(lev)%ny
          nz = grid(lev)%nz
          nd = size(grid(lev)%cA,1)
          incx=grid(lev)%incx / 2
          incy=grid(lev)%incy / 2
          ngx=grid(lev)%ngx
          ngy=grid(lev)%ngy

          ! Gather cores by quadruplets (and marginally by pair, for the coarsest grid)

          ! cores having the same family index share the same subdomain
          ! All cores in a family do the same work

          ! Integer division; for instance, cores with pi=0,1 have the same family index
          family=(pi/incx)*incx*incy + (npx)*incy*(pj/incy)

          nextfamily = (pi/(2*incx))*incx*incy*4 + (npx)*2*incy*(pj/(incy*2))

          ! Assign a color to each core: make a cycling ramp index
          ! through 2 or 4 close families
          ! - cores having the same color should be a pair or a quadruplet
          ! - colors are all distinct *within* a family
          color=nextfamily + mod(pi,incx)+mod(pj,incy)*incx

          N=incx*npx;
          key = mod(mod(family,N)/(incx*incy),2)+2*mod( (family/N),2)

          grid(lev)%key=key

          call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, localcomm, ierr)
          grid(lev)%localcomm = localcomm

          ! this dummy 3D array is to store the restriction from lev-1, before the gathering
          ! its size can be deduced from the size after the gathering

          nx = nx/ngx ! ngx is 1 or 2 (and generally 2)
          ny = ny/ngy ! ngy is 1 or 2 (and generally 2)
          allocate(grid(lev)%dummy3(nz,0:ny+1,0:nx+1))

          allocate(grid(lev)%gatherbuffer2D(0:ny+1,0:nx+1,0:ngx-1,0:ngy-1))
          allocate(grid(lev)%gatherbuffer (1:nz  ,0:ny+1,0:nx+1,0:ngx-1,0:ngy-1))

          ! number of elements of dummy3
          grid(lev)%Ng2D=(nx+2)*(ny+2)
          grid(lev)%Ng  =(nx+2)*(ny+2)*nz
          grid(lev)%Ngp  =(nx+2)*(ny+2)*(nz+1)

       endif

    enddo

  end subroutine define_gather_informations

  !---------------------------------------------------------------------
  subroutine print_grids()
    integer(kind=ip) :: lev,ierr

    if (myrank==0) write(*,*)'    - print grid information:'

    call MPI_Barrier( MPI_COMM_WORLD ,ierr)
    if (myrank.eq.0)then
       do lev=1,nlevs
          if (grid(lev)%gather.eq.0)then
             write(*,100)"    lev=",lev,": ", &
                  grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                  " on ",grid(lev)%npx,' x',grid(lev)%npy," procs"
          else
             write(*,100)"    lev=",lev,": ", &
                  grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                  " on ",grid(lev)%npx,' x',grid(lev)%npy," procs / gather"
          endif
       enddo
    endif
100 format (A,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)

  end subroutine print_grids

  !---------------------------------------------------------------------
  subroutine grids_dealloc()

    integer(kind=ip) :: lev

    ! Check if we deallocate everything

    if (associated(grid)) then

       do lev = 1, nlevs

          if (associated(grid(lev)%cA))             deallocate(grid(lev)%cA)
          if (associated(grid(lev)%p))              deallocate(grid(lev)%p)
          if (associated(grid(lev)%b))              deallocate(grid(lev)%b)
          if (associated(grid(lev)%r))              deallocate(grid(lev)%r)
          if (associated(grid(lev)%dx))             deallocate(grid(lev)%dx)
          if (associated(grid(lev)%dy))             deallocate(grid(lev)%dy)
          if (associated(grid(lev)%dz))             deallocate(grid(lev)%dz)
          if (associated(grid(lev)%dzw))            deallocate(grid(lev)%dzw)
          if (associated(grid(lev)%zxdy))           deallocate(grid(lev)%zxdy)
          if (associated(grid(lev)%zydx))           deallocate(grid(lev)%zydx)
          if (associated(grid(lev)%alpha))          deallocate(grid(lev)%alpha)
          if (associated(grid(lev)%dummy3Dnz))      deallocate(grid(lev)%dummy3Dnz)
          if (associated(grid(lev)%dummy3Dnzp))     deallocate(grid(lev)%dummy3Dnzp)
          if (associated(grid(lev)%dummy3))         deallocate(grid(lev)%dummy3)
          if (associated(grid(lev)%gatherbuffer2D)) deallocate(grid(lev)%gatherbuffer2D)
          if (associated(grid(lev)%gatherbuffer))   deallocate(grid(lev)%gatherbuffer)

       enddo

       deallocate(grid)

    endif

    if (associated(gbuffers)) then

       do lev = 1, nlevs
          if (associated(gbuffers(lev)%sendN2D1)) deallocate(gbuffers(lev)%sendN2D1)
          if (associated(gbuffers(lev)%recvN2D1)) deallocate(gbuffers(lev)%recvN2D1)
          if (associated(gbuffers(lev)%sendS2D1)) deallocate(gbuffers(lev)%sendS2D1)
          if (associated(gbuffers(lev)%recvS2D1)) deallocate(gbuffers(lev)%recvS2D1)
          if (associated(gbuffers(lev)%sendE2D1)) deallocate(gbuffers(lev)%sendE2D1)
          if (associated(gbuffers(lev)%recvE2D1)) deallocate(gbuffers(lev)%recvE2D1)
          if (associated(gbuffers(lev)%sendW2D1)) deallocate(gbuffers(lev)%sendW2D1)
          if (associated(gbuffers(lev)%recvW2D1)) deallocate(gbuffers(lev)%recvW2D1)

          if (associated(gbuffers(lev)%sendN)) deallocate(gbuffers(lev)%sendN)
          if (associated(gbuffers(lev)%recvN)) deallocate(gbuffers(lev)%recvN)
          if (associated(gbuffers(lev)%sendS)) deallocate(gbuffers(lev)%sendS)
          if (associated(gbuffers(lev)%recvS)) deallocate(gbuffers(lev)%recvS)
          if (associated(gbuffers(lev)%sendE)) deallocate(gbuffers(lev)%sendE)
          if (associated(gbuffers(lev)%recvE)) deallocate(gbuffers(lev)%recvE)
          if (associated(gbuffers(lev)%sendW)) deallocate(gbuffers(lev)%sendW)
          if (associated(gbuffers(lev)%recvW)) deallocate(gbuffers(lev)%recvW)
          if (associated(gbuffers(lev)%sendSW)) deallocate(gbuffers(lev)%sendSW)
          if (associated(gbuffers(lev)%recvSW)) deallocate(gbuffers(lev)%recvSW)
          if (associated(gbuffers(lev)%sendSE)) deallocate(gbuffers(lev)%sendSE)
          if (associated(gbuffers(lev)%recvSE)) deallocate(gbuffers(lev)%recvSE)
          if (associated(gbuffers(lev)%sendNW)) deallocate(gbuffers(lev)%sendNW)
          if (associated(gbuffers(lev)%recvNW)) deallocate(gbuffers(lev)%recvNW)
          if (associated(gbuffers(lev)%sendNE)) deallocate(gbuffers(lev)%sendNE)
          if (associated(gbuffers(lev)%recvNE)) deallocate(gbuffers(lev)%recvNE)
          if (associated(gbuffers(lev)%sendNp)) deallocate(gbuffers(lev)%sendNp)
          if (associated(gbuffers(lev)%recvNp)) deallocate(gbuffers(lev)%recvNp)
          if (associated(gbuffers(lev)%sendSp)) deallocate(gbuffers(lev)%sendSp)
          if (associated(gbuffers(lev)%recvSp)) deallocate(gbuffers(lev)%recvSp)
          if (associated(gbuffers(lev)%sendEp)) deallocate(gbuffers(lev)%sendEp)
          if (associated(gbuffers(lev)%recvEp)) deallocate(gbuffers(lev)%recvEp)
          if (associated(gbuffers(lev)%sendWp)) deallocate(gbuffers(lev)%sendWp)
          if (associated(gbuffers(lev)%recvWp)) deallocate(gbuffers(lev)%recvWp)
          if (associated(gbuffers(lev)%sendSWp)) deallocate(gbuffers(lev)%sendSWp)
          if (associated(gbuffers(lev)%recvSWp)) deallocate(gbuffers(lev)%recvSWp)
          if (associated(gbuffers(lev)%sendSEp)) deallocate(gbuffers(lev)%sendSEp)
          if (associated(gbuffers(lev)%recvSEp)) deallocate(gbuffers(lev)%recvSEp)
          if (associated(gbuffers(lev)%sendNWp)) deallocate(gbuffers(lev)%sendNWp)
          if (associated(gbuffers(lev)%recvNWp)) deallocate(gbuffers(lev)%recvNWp)
          if (associated(gbuffers(lev)%sendNEp)) deallocate(gbuffers(lev)%sendNEp)
          if (associated(gbuffers(lev)%recvNEp)) deallocate(gbuffers(lev)%recvNEp)

          if (associated(gbuffers(lev)%sendN3D2)) deallocate(gbuffers(lev)%sendN3D2)
          if (associated(gbuffers(lev)%recvN3D2)) deallocate(gbuffers(lev)%recvN3D2)
          if (associated(gbuffers(lev)%sendS3D2)) deallocate(gbuffers(lev)%sendS3D2)
          if (associated(gbuffers(lev)%recvS3D2)) deallocate(gbuffers(lev)%recvS3D2)
          if (associated(gbuffers(lev)%sendE3D2)) deallocate(gbuffers(lev)%sendE3D2)
          if (associated(gbuffers(lev)%recvE3D2)) deallocate(gbuffers(lev)%recvE3D2)
          if (associated(gbuffers(lev)%sendW3D2)) deallocate(gbuffers(lev)%sendW3D2)
          if (associated(gbuffers(lev)%recvW3D2)) deallocate(gbuffers(lev)%recvW3D2)
          if (associated(gbuffers(lev)%sendSW3D2)) deallocate(gbuffers(lev)%sendSW3D2)
          if (associated(gbuffers(lev)%recvSW3D2)) deallocate(gbuffers(lev)%recvSW3D2)
          if (associated(gbuffers(lev)%sendSE3D2)) deallocate(gbuffers(lev)%sendSE3D2)
          if (associated(gbuffers(lev)%recvSE3D2)) deallocate(gbuffers(lev)%recvSE3D2)
          if (associated(gbuffers(lev)%sendNW3D2)) deallocate(gbuffers(lev)%sendNW3D2)
          if (associated(gbuffers(lev)%recvNW3D2)) deallocate(gbuffers(lev)%recvNW3D2)
          if (associated(gbuffers(lev)%sendNE3D2)) deallocate(gbuffers(lev)%sendNE3D2)
          if (associated(gbuffers(lev)%recvNE3D2)) deallocate(gbuffers(lev)%recvNE3D2)
          if (associated(gbuffers(lev)%sendN3D2p)) deallocate(gbuffers(lev)%sendN3D2p)
          if (associated(gbuffers(lev)%recvN3D2p)) deallocate(gbuffers(lev)%recvN3D2p)
          if (associated(gbuffers(lev)%sendS3D2p)) deallocate(gbuffers(lev)%sendS3D2p)
          if (associated(gbuffers(lev)%recvS3D2p)) deallocate(gbuffers(lev)%recvS3D2p)
          if (associated(gbuffers(lev)%sendE3D2p)) deallocate(gbuffers(lev)%sendE3D2p)
          if (associated(gbuffers(lev)%recvE3D2p)) deallocate(gbuffers(lev)%recvE3D2p)
          if (associated(gbuffers(lev)%sendW3D2p)) deallocate(gbuffers(lev)%sendW3D2p)
          if (associated(gbuffers(lev)%recvW3D2p)) deallocate(gbuffers(lev)%recvW3D2p)
          if (associated(gbuffers(lev)%sendSW3D2p)) deallocate(gbuffers(lev)%sendSW3D2p)
          if (associated(gbuffers(lev)%recvSW3D2p)) deallocate(gbuffers(lev)%recvSW3D2p)
          if (associated(gbuffers(lev)%sendSE3D2p)) deallocate(gbuffers(lev)%sendSE3D2p)
          if (associated(gbuffers(lev)%recvSE3D2p)) deallocate(gbuffers(lev)%recvSE3D2p)
          if (associated(gbuffers(lev)%sendNW3D2p)) deallocate(gbuffers(lev)%sendNW3D2p)
          if (associated(gbuffers(lev)%recvNW3D2p)) deallocate(gbuffers(lev)%recvNW3D2p)
          if (associated(gbuffers(lev)%sendNE3D2p)) deallocate(gbuffers(lev)%sendNE3D2p)
          if (associated(gbuffers(lev)%recvNE3D2p)) deallocate(gbuffers(lev)%recvNE3D2p)

          if (associated(gbuffers(lev)%sendN4D)) deallocate(gbuffers(lev)%sendN4D)
          if (associated(gbuffers(lev)%recvN4D)) deallocate(gbuffers(lev)%recvN4D)

          if (associated(gbuffers(lev)%sendS4D)) deallocate(gbuffers(lev)%sendS4D)
          if (associated(gbuffers(lev)%recvS4D)) deallocate(gbuffers(lev)%recvS4D)
          if (associated(gbuffers(lev)%sendE4D)) deallocate(gbuffers(lev)%sendE4D)
          if (associated(gbuffers(lev)%recvE4D)) deallocate(gbuffers(lev)%recvE4D)
          if (associated(gbuffers(lev)%sendW4D)) deallocate(gbuffers(lev)%sendW4D)
          if (associated(gbuffers(lev)%recvW4D)) deallocate(gbuffers(lev)%recvW4D)

          if (associated(gbuffers(lev)%sendSW4D)) deallocate(gbuffers(lev)%sendSW4D)
          if (associated(gbuffers(lev)%recvSW4D)) deallocate(gbuffers(lev)%recvSW4D)
          if (associated(gbuffers(lev)%sendSE4D)) deallocate(gbuffers(lev)%sendSE4D)
          if (associated(gbuffers(lev)%recvSE4D)) deallocate(gbuffers(lev)%recvSE4D)

          if (associated(gbuffers(lev)%sendNW4D)) deallocate(gbuffers(lev)%sendNW4D)
          if (associated(gbuffers(lev)%recvNW4D)) deallocate(gbuffers(lev)%recvNW4D)
          if (associated(gbuffers(lev)%sendNE4D)) deallocate(gbuffers(lev)%sendNE4D)
          if (associated(gbuffers(lev)%recvNE4D)) deallocate(gbuffers(lev)%recvNE4D)

       enddo

       deallocate(gbuffers)

    endif

  end subroutine grids_dealloc

end module mg_grids
