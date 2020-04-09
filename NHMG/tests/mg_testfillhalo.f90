program mg_testfillhalo

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_mpi_exchange
  use nhmg

  implicit none

  integer(kind=ip) :: nx, ny, nz  ! local dimensions
  integer(kind=ip) :: np, ierr, rank

  real(kind=rp), dimension(:,:)  , pointer :: a2D
  real(kind=rp), dimension(:,:,:), pointer :: a3D
  real(kind=rp), dimension(:,:,:), pointer :: a3Dp

  integer(kind=ip)  :: nbv =9
  real(kind=rp), dimension(:)  , pointer :: myres
  real(kind=rp), dimension(:,:), pointer :: res

  integer(kind=ip):: nxg  = 64       ! global x dimension
  integer(kind=ip):: nyg  = 64       ! global y dimension
  integer(kind=ip):: nzg  = 64       ! z dimension
  integer(kind=ip):: npxg  = 1       ! number of processes in x
  integer(kind=ip):: npyg  = 1       ! number of processes in y

  integer(kind=ip)  :: lun_nml = 4
  integer(kind=ip)  :: ii
  logical :: nml_exist=.false.

  namelist/fhparam/ &
       nxg        , &
       nyg        , &
       nzg        , &
       npxg       , &
       npyg 

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  !---------------------!
  !- Namelist (or not) -!
  !---------------------!
  !- Check if a fh_namelist file exist
  inquire(file='fh_namelist', exist=nml_exist)

  !- Read namelist file if it is present, else use default values
  if (nml_exist) then
     if (rank == 0) write(*,*)' Reading fh_namelist file'
     open(unit=lun_nml, File='fh_namelist', ACTION='READ')
     rewind(unit=lun_nml)
     read(unit=lun_nml, nml=fhparam)
  endif

  if (rank == 0) then
     write(*,*)'  test fill halo parameters:'
     write(*,*)'  - nxg     : ', nxg
     write(*,*)'  - nyg     : ', nyg
     write(*,*)'  - nzg     : ', nzg
     write(*,*)'  - npxg    : ', npxg 
     write(*,*)'  - npyg    : ', npyg
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
  !- Setup fill halo   -!
  !---------------------!
  if (rank == 0) then
     open(unit=20,file="results.txt",form="formatted")
     write(20,*)''
     write(20,*)'Process topology:'
     write(20,*)'-----------------'
     allocate(res(npyg,npxg))
     res = reshape([(ii,ii=0,np-1)],[npyg,npxg])
     do ii=npyg,1,-1
        write(20,'(16(x,I3))')int(res(:,ii))
     enddo
     deallocate(res)
  endif

  if (rank == 0) write(*,*)' Allocate fill halo test arrays'

  ! S E N W SW SE NE NW '
  allocate( myres(nbv))
  allocate( res(nbv,0:np-1))

  allocate( a2D (       0:ny+1,0:ny+1))
  allocate( a3D (1:nz  ,0:ny+1,0:nx+1))
  allocate( a3Dp(1:nz+1,0:ny+1,0:nx+1))

  if (rank == 0) write(*,*)' Initialize fill halo test arrays'

  a2D (:,:)   = -one
  a3D (:,:,:) = -one
  a3Dp(:,:,:) = -one
  a2D (  1:ny,1:nx) = rank
  a3D (:,1:ny,1:nx) = rank
  a3Dp(:,1:ny,1:nx) = rank

  if (netcdf_output) then
     call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2D.nc' ,rank=rank)
     call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3D.nc' ,rank=rank)
     call write_netcdf(a3Dp,vname='a3Dp',netcdf_file_name='a3Dp.nc',rank=rank)
  endif

  !-----------------!
  !- Fill halo a2D -!
  !-----------------!
  call fill_halo(1,a2D)

  myres(1) = sum(a2D(   0,1:nx))/nx ! S
  myres(2) = sum(a2D(1:ny,nx+1))/ny ! E
  myres(3) = sum(a2D(ny+1,1:nx))/nx ! N
  myres(4) = sum(a2D(1:ny,   0))/ny ! W
  myres(5) =     a2D(   0,   0)     ! SW
  myres(6) =     a2D(   0,nx+1)     ! SE
  myres(7) =     a2D(ny+1,nx+1)     ! NE
  myres(8) =     a2D(ny+1,   0)     ! NW
  myres(9) = sum(a2D(1:ny,1:nx))/(ny*nx) ! domain
  call MPI_GATHER(myres,nbv,MPI_DOUBLE_PRECISION,res,nbv, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (rank == 0 ) then
     write(20,*)''
     write(20,*)'a2D results:'
     write(20,*)'------------'
     write(20,*)'r    S      E      N      W      SW     SE     NE     NW     =rank'
     do ii=0,np-1
        write(20,'(I2,a,9(x,F6.2))') ii,':',res(:,ii)
     enddo
  endif

  !-----------------!
  !- Fill halo a3D -!
  !-----------------!
  call fill_halo(1,a3D)

  myres(1) = sum(a3D(:,   0,1:nx))/(nx*nz) ! S
  myres(2) = sum(a3D(:,1:ny,nx+1))/(ny*nz) ! E
  myres(3) = sum(a3D(:,ny+1,1:nx))/(nx*nz) ! N
  myres(4) = sum(a3D(:,1:ny,   0))/(ny*nz) ! W
  myres(5) = sum(a3D(:,   0,   0))/ nz     ! SW
  myres(6) = sum(a3D(:,   0,nx+1))/ nz     ! SE
  myres(7) = sum(a3D(:,ny+1,nx+1))/ nz     ! NE
  myres(8) = sum(a3D(:,ny+1,   0))/ nz     ! NW
  myres(9) = sum(a3D(:,1:ny,1:nx))/(ny*nx*nz) ! domain
  call MPI_GATHER(myres,nbv,MPI_DOUBLE_PRECISION,res,nbv, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (rank == 0 ) then
     write(20,*)''
     write(20,*)'a3D(nz,:,:) results:'
     write(20,*)'--------------------'
     write(20,*)'r    S      E      N      W      SW     SE     NE     NW     =rank'
     do ii=0,np-1
        write(20,'(I2,a,9(x,F6.2))') ii,':',res(:,ii)
     enddo
  endif

  !------------------!
  !- Fill halo a3Dp -!
  !------------------!
  call fill_halo(1,a3Dp)

  myres(1) = sum(a3Dp(:,   0,1:nx))/(nx*(nz+1)) ! S
  myres(2) = sum(a3Dp(:,1:ny,nx+1))/(ny*(nz+1)) ! E
  myres(3) = sum(a3Dp(:,ny+1,1:nx))/(nx*(nz+1)) ! N
  myres(4) = sum(a3Dp(:,1:ny,   0))/(ny*(nz+1)) ! W
  myres(5) = sum(a3Dp(:,   0,   0))/ (nz+1)     ! SW
  myres(6) = sum(a3Dp(:,   0,nx+1))/ (nz+1)     ! SE
  myres(7) = sum(a3Dp(:,ny+1,nx+1))/ (nz+1)     ! NE
  myres(8) = sum(a3Dp(:,ny+1,   0))/ (nz+1)     ! NW
  myres(9) = sum(a3Dp(:,1:ny,1:nx))/(ny*nx*(nz+1)) ! domain
  call MPI_GATHER(myres,nbv,MPI_DOUBLE_PRECISION,res,nbv, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (rank == 0 ) then
     write(20,*)''
     write(20,*)'a3D(nz+1,:,:) results:'
     write(20,*)'----------------------'
     write(20,*)'r     S      E      N      W      SW     SE     NE     NW    =rank'
     do ii=0,np-1
        write(20,'(I2,a,9(x,F6.2))') ii,':',res(:,ii)
     enddo
  endif

  if (netcdf_output) then
     if (rank == 0) write(*,*)' Writting arrays after fill_halo'
     call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2Dfh.nc' ,rank=rank)
     call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3Dfh.nc' ,rank=rank)
     call write_netcdf(a3Dp,vname='a3Dp',netcdf_file_name='a3Dpfh.nc',rank=rank)
  endif

  if (rank == 0) write(*,*)' Initialize fill halo test arrays (u)'

  a2D (:,:)   = -one
  a3D (:,:,:) = -one
  a2D (  1:ny,1:nx) = rank
  a3D (:,1:ny,1:nx) = rank


  !-------------------!
  !- Fill halo a2D u -!
  !-------------------!
  call fill_halo(1,a2D,lbc_null='u')

  myres(1) = sum(a2D(   0,1:nx))/nx ! S
  myres(2) = sum(a2D(1:ny,nx+1))/ny ! E
  myres(3) = sum(a2D(ny+1,1:nx))/nx ! N
  myres(4) = sum(a2D(1:ny,   0))/ny ! W
  myres(5) =     a2D(   0,   0)     ! SW
  myres(6) =     a2D(   0,nx+1)     ! SE
  myres(7) =     a2D(ny+1,nx+1)     ! NE
  myres(8) =     a2D(ny+1,   0)     ! NW
  myres(9) = sum(a2D(1:ny,1:nx))/(ny*nx) ! domain
  call MPI_GATHER(myres,nbv,MPI_DOUBLE_PRECISION,res,nbv, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (rank == 0 ) then
     write(20,*)''
     write(20,*)'a2D u results:'
     write(20,*)'--------------'
     write(20,*)'r    S      E      N      W      SW     SE     NE     NW    =rank'
     do ii=0,np-1
        write(20,'(I2,a,9(x,F6.2))') ii,':',res(:,ii)
     enddo
  endif

  !-------------------!
  !- Fill halo a3D u -!
  !-------------------!
  call fill_halo(1,a3D,lbc_null='u')

  myres(1) = sum(a3D(:,   0,1:nx))/(nx*nz) ! S
  myres(2) = sum(a3D(:,1:ny,nx+1))/(ny*nz) ! E
  myres(3) = sum(a3D(:,ny+1,1:nx))/(nx*nz) ! N
  myres(4) = sum(a3D(:,1:ny,   0))/(ny*nz) ! W
  myres(5) = sum(a3D(:,   0,   0))/ nz     ! SW
  myres(6) = sum(a3D(:,   0,nx+1))/ nz     ! SE
  myres(7) = sum(a3D(:,ny+1,nx+1))/ nz     ! NE
  myres(8) = sum(a3D(:,ny+1,   0))/ nz     ! NW
  myres(9) = sum(a3D(:,1:ny,1:nx))/(ny*nx*nz) ! domain
  call MPI_GATHER(myres,nbv,MPI_DOUBLE_PRECISION,res,nbv, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (rank == 0 ) then
     write(20,*)''
     write(20,*)'a3D(nz,:,:) u results:'
     write(20,*)'----------------------'
     write(20,*)'r    S      E      N      W      SW     SE     NE     NW    =rank'
     do ii=0,np-1
        write(20,'(I2,a,9(x,F6.2))') ii,':',res(:,ii)
     enddo
  endif

  if (netcdf_output) then
     if (rank == 0) write(*,*)' Writting arrays (u) after fill_halo'
     call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2Dufh.nc' ,rank=rank)
     call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3Dufh.nc' ,rank=rank)
  endif

  if (rank == 0) write(*,*)' Initialize fill halo test arrays (v)'

  a2D (:,:)   = -one
  a3D (:,:,:) = -one
  a2D (  1:ny,1:nx) = rank
  a3D (:,1:ny,1:nx) = rank

  !-------------------!
  !- Fill halo a2D v -!
  !-------------------!
  call fill_halo(1,a2D,lbc_null='v')

  myres(1) = sum(a2D(   0,1:nx))/nx ! S
  myres(2) = sum(a2D(1:ny,nx+1))/ny ! E
  myres(3) = sum(a2D(ny+1,1:nx))/nx ! N
  myres(4) = sum(a2D(1:ny,   0))/ny ! W
  myres(5) =     a2D(   0,   0)     ! SW
  myres(6) =     a2D(   0,nx+1)     ! SE
  myres(7) =     a2D(ny+1,nx+1)     ! NE
  myres(8) =     a2D(ny+1,   0)     ! NW
  myres(9) = sum(a2D(1:ny,1:nx))/(ny*nx) ! domain
  call MPI_GATHER(myres,nbv,MPI_DOUBLE_PRECISION,res,nbv, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (rank == 0 ) then
     write(20,*)''
     write(20,*)'a2D v results:'
     write(20,*)'--------------'
     write(20,*)'r    S      E      N      W      SW     SE     NE     NW    =rank'
     do ii=0,np-1
        write(20,'(I2,a,9(x,F6.2))') ii,':',res(:,ii)
     enddo
  endif

  !-------------------!
  !- Fill halo a3D v -!
  !-------------------!
  call fill_halo(1,a3D,lbc_null='v')

  myres(1) = sum(a3D(:,   0,1:nx))/(nx*nz) ! S
  myres(2) = sum(a3D(:,1:ny,nx+1))/(ny*nz) ! E
  myres(3) = sum(a3D(:,ny+1,1:nx))/(nx*nz) ! N
  myres(4) = sum(a3D(:,1:ny,   0))/(ny*nz) ! W
  myres(5) = sum(a3D(:,   0,   0))/ nz     ! SW
  myres(6) = sum(a3D(:,   0,nx+1))/ nz     ! SE
  myres(7) = sum(a3D(:,ny+1,nx+1))/ nz     ! NE
  myres(8) = sum(a3D(:,ny+1,   0))/ nz     ! NW
  myres(9) = sum(a3D(:,1:ny,1:nx))/(ny*nx*nz) ! domain
  call MPI_GATHER(myres,nbv,MPI_DOUBLE_PRECISION,res,nbv, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (rank == 0 ) then
     write(20,*)''
     write(20,*)'a3D(nz,:,:) v results:'
     write(20,*)'----------------------'
     write(20,*)'r    S      E      N      W      SW     SE     NE     NW    =rank'
     do ii=0,np-1
        write(20,'(I2,a,9(x,F6.2))') ii,':',res(:,ii)
     enddo
  endif

  if (netcdf_output) then
     if (rank == 0) write(*,*)' Writting arrays (v) after fill_halo'
     call write_netcdf(a2D ,vname='a2D' ,netcdf_file_name='a2Dvfh.nc' ,rank=rank)
     call write_netcdf(a3D ,vname='a3D' ,netcdf_file_name='a3Dvfh.nc' ,rank=rank)
  endif

  !---------------------!
  !- Deallocate memory -!
  !---------------------!
  if (rank == 0) write(*,*)' Clean memory before to finish the program.'
  call nhmg_clean()

  !----------------------!
  !- End Bench-fillhalo -!
  !----------------------!
  call mpi_finalize(ierr)

end program mg_testfillhalo

