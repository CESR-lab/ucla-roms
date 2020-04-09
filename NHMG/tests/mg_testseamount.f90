program mg_testseamount

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_zr_hz
  use mg_mpi_exchange_ijk
  use mg_setup_tests
  use nhmg

  implicit none

  integer(kind=ip) :: nx, ny, nz  ! local dimensions
  integer(kind=ip) :: it          ! iteration
  integer(kind=ip) :: np, ierr, rank

  real(kind=rp), dimension(:,:), pointer :: dx, dy
  real(kind=rp), dimension(:,:), pointer :: dxu, dyv
  real(kind=rp), dimension(:,:), pointer :: zeta, h
  real(kind=rp), dimension(:,:,:), pointer :: z_r,z_w
  real(kind=rp), dimension(:,:,:), pointer :: Hz

  real(kind=rp), dimension(:,:,:), allocatable, target :: u,v,w
  real(kind=rp), dimension(:,:,:), pointer :: up,vp,wp

  integer(kind=ip) :: nit=1
  integer(kind=ip):: nxg  = 64       ! global x dimension
  integer(kind=ip):: nyg  = 64       ! global y dimension
  integer(kind=ip):: nzg  = 64       ! z dimension
  integer(kind=ip):: npxg  = 1       ! number of processes in x
  integer(kind=ip):: npyg  = 1       ! number of processes in y
  real(kind=rp) :: Lx      = 1.e4_rp ! Domain length in x [meter]
  real(kind=rp) :: Ly      = 1.e4_rp ! Domain length in y [meter]
  real(kind=rp) :: Htot    = 4.e3_rp ! Depth [meter]
  real(kind=rp) :: hc      = 4.e3_rp !
  real(kind=rp) :: theta_b = 0._rp   !
  real(kind=rp) :: theta_s = 0._rp   !

  integer(kind=ip)  :: lun_nml = 4
  logical :: nml_exist=.false.

  integer(kind=ip)  :: is

  namelist/tsparam/ &
       nit        , &
       nxg        , &
       nyg        , &
       nzg        , &
       npxg       , &
       npyg       , &
       Lx         , &
       Ly         , &
       Htot       , &
       hc         , &
       theta_b    , &
       theta_s

  call tic(1,'test_seamount')

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  !---------------------!
  !- Namelist (or not) -!
  !---------------------!

  !- Check if a ts_namelist file exist
  inquire(file='ts_namelist', exist=nml_exist)

  !- Read namelist file if it is present, else use default values
  if (nml_exist) then
     if (rank == 0) write(*,*)' Reading ts_namelist file'
     open(unit=lun_nml, File='ts_namelist', ACTION='READ')
     rewind(unit=lun_nml)
     read(unit=lun_nml, nml=tsparam)
  endif

  if (rank == 0) then
     write(*,*)'  test seamount parameters:'
     write(*,*)'  - nit     : ', nit 
     write(*,*)'  - nxg     : ', nxg
     write(*,*)'  - nyg     : ', nyg
     write(*,*)'  - nzg     : ', nzg
     write(*,*)'  - npxg    : ', npxg 
     write(*,*)'  - npyg    : ', npyg
     write(*,*)'  - Lx      : ', Lx
     write(*,*)'  - Ly      : ', Ly
     write(*,*)'  - Htot    : ', Htot
     write(*,*)'  - hc      : ', hc
     write(*,*)'  - theta_b : ', theta_b
     write(*,*)'  - theta_s : ', theta_s
     write(*,*)'  '
  endif

  !---------------------!
  !- Global/local dim  -!
  !---------------------!

  if (np /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !", np, npxg, npyg
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  !-------------------!
  !- Initialise nhmg -!
  !-------------------!
  call nhmg_init(nx,ny,nz,npxg,npyg)

  !---------------------!
  !- Setup seamount    -!
  !---------------------!
  if (rank == 0) write(*,*)' Initialise seamount test'

  if (rank==0) then
     write(*,'(A,3(x,F12.2))')'    - Lx, Ly, Htot        :', Lx, Ly, Htot
     write(*,'(A,3(x,F12.2))')'    - hc, theta_b, theta_s:', hc, theta_b, theta_s
  endif

  is= 2

  allocate(    h(1-is:nx+is,1-is:ny+is))
  allocate( zeta(1-is:nx+is,1-is:ny+is))

  allocate(  dx (1-is:nx+is,1-is:ny+is))
  allocate(  dy (1-is:nx+is,1-is:ny+is))
  allocate(  dxu(1-is:nx+is,1-is:ny+is))
  allocate(  dyv(1-is:nx+is,1-is:ny+is))

  allocate(   z_r(1-is:nx+is,1-is:ny+is,1:nz))
  allocate(   z_w(1-is:nx+is,1-is:ny+is,1:nz+1))
  allocate(   Hz (1-is:nx+is,1-is:ny+is,1:nz))

  call setup_seamount(  &
       nx,ny,npxg,npyg, &
       Lx,Ly,Htot,      &
       dx,dy,           &
       zeta,h        )

  !- stretching vertical grid -!   
  call setup_zr_zw_hz(hc,theta_b,theta_s,zeta,h,z_r,z_w,Hz,'new_s_coord')
  !- linear vertical grid -!
  !  call setup_zr_zw_hz(h,z_r,z_w,Hz)

  call nhmg_matrices(nx,ny,nz,z_r(0:nx+1,0:ny+1,1:nz),Hz(0:nx+1,0:ny+1,1:nz),dx,dy)

  !-------------------------------------!
  !- U,V,W initialisation (model vars) -!
  !-------------------------------------!
  if (rank == 0) write(*,*)' Allocate u, v, w'

  allocate(u(1-is:nx+is,1-is:ny+is,1:nz))
  allocate(v(1-is:nx+is,1-is:ny+is,1:nz))
  allocate(w(1-is:nx+is,1-is:ny+is,0:nz))

  do it = 1, nit

     if (rank==0) then
        write(*,*)' U=0, V=0 and W=-1 except at bottom'
     endif
     u(:,:,:)    =  0._8
     v(:,:,:)    =  0._8
     w(:,:,0)    =  0._8
     w(:,:,1:nz) = -1._8

     if (netcdf_output) then
        call write_netcdf(u,vname='u',netcdf_file_name='u.nc',rank=rank,iter=it)
        call write_netcdf(v,vname='v',netcdf_file_name='v.nc',rank=rank,iter=it)
        call write_netcdf(w,vname='w',netcdf_file_name='w.nc',rank=rank,iter=it)
     endif

     grid(1)%p=0._rp

     !--------------------!
     !- Call nhmg solver -!
     !--------------------!
     if (rank == 0) write(*,*)'Call nhmg solver',it
     call nhmg_solve(u,v,w,z_w,Hz,.true.)

     if (netcdf_output) then
        call write_netcdf(u,vname='uc',netcdf_file_name='uc.nc',rank=rank,iter=it)
        call write_netcdf(v,vname='vc',netcdf_file_name='vc.nc',rank=rank,iter=it)
        call write_netcdf(w,vname='wc',netcdf_file_name='wc.nc',rank=rank,iter=it)
     endif

  enddo

  !---------------------!
  !- Deallocate memory -!
  !---------------------!
  if (rank == 0) write(*,*)' Clean memory before to finish the program.'
  call nhmg_clean()

  !----------------------!
  !- End Bench-seamount -!
  !----------------------!
  call mpi_finalize(ierr)

  call toc(1,'test_seamount')
  if(rank == 0) call print_tictoc(rank)

end program mg_testseamount

