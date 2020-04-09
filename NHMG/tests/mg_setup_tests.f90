module mg_setup_tests

  use mg_mpi 
  use mg_tictoc
  use mg_mpi_exchange_ijk
  use mg_netcdf_out

  implicit none

contains

  !-------------------------------------------------------------------------     
  subroutine setup_cuc(    &
       nx, ny, npxg, npyg, &
       dx, dy,             &
       zeta, h         )

    integer(kind=ip), intent(in) :: nx,ny  ! local dimsg
    integer(kind=ip), intent(in) :: npxg,npyg ! nb procs
    real(kind=rp), dimension(:,:), pointer, intent(out) :: dx, dy, zeta, h

    integer(kind=ip):: is_err,nc_id,varid
    integer(kind=ip):: i,j,i0,j0,pi,pj
    real(kind=rp), dimension(:,:), allocatable   :: dummy2d
    integer(kind=ip), dimension(2)   :: starts,counts

    character*80 :: file,varname
    integer(kind=ip) :: ierr

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    pj = myrank/npxg
    pi = mod(myrank,npxg)

    ! dummy 2D to read from netcdf
    allocate(dummy2d(1:nx,1:ny))

    i0=2+pi*nx
    j0=2+pj*ny

    starts(1)=i0
    starts(2)=j0
    counts(1)=nx
    counts(2)=ny

    !!------------------------------------------------!!
    !! Here enter the directory where cuc_nhgrd.nc is !!
    !!------------------------------------------------!!
    file='../DATA/cuc_nhgrd.nc'
    !file='/workgpfs/rech/dgw/rdgw004/MGROMS/DATA' ! IDRIS ADA $workdir

    is_err = nf90_open( trim(file), NF90_NOWRITE ,nc_id  )

    if (is_err /= nf90_noerr) then
       write(*,*)'Error: problem to open file: ', trim(file)
       stop
    end if

    ! --- read h ---
    varname='h'
    is_err = nf90_inq_varid(nc_id,trim(varname),varid)
    is_err = nf90_get_var(nc_id,varid,dummy2d,starts,counts)

    do i=1,nx
       do j=1,ny
          h(i,j)=dummy2d((i-1)+1,(j-1)+1)
       enddo
    enddo

    ! --- read pn ---
    varname='pn'
    is_err = nf90_inq_varid(nc_id,trim(varname),varid)
    is_err = nf90_get_var(nc_id,varid,dummy2d,starts,counts)

    do j=1,ny
       do i=1,nx
          dy(i,j)=1._rp/dummy2d((i-1)+1,(j-1)+1)
       enddo
    enddo

    ! --- read pm ---
    varname='pm'
    is_err = nf90_inq_varid(nc_id,trim(varname),varid)
    is_err = nf90_get_var(nc_id,varid,dummy2d,starts,counts)

    do j=1,ny
       do i=1,nx
          dx(i,j)=1._rp/dummy2d((i-1)+1,(j-1)+1)
       enddo
    enddo

    is_err = nf90_close(nc_id)

    call fill_halo_ijk(nx,ny,npxg,npyg,dx)
    call fill_halo_ijk(nx,ny,npxg,npyg,dy)
    call fill_halo_ijk(nx,ny,npxg,npyg, h)

    do j= 1,ny
       do i= 1,nx
          dx(i,j)=max(1._rp,dx(i,j))
          dy(i,j)=max(1._rp,dy(i,j))
       enddo
    enddo

    zeta(:,:) = 0._rp

  end subroutine setup_cuc

  !-------------------------------------------------------------------------     
  subroutine setup_seamount( &
       nx, ny, npxg, npyg,   &
       Lx, Ly, Htot,         &
       dx, dy,               &
       zeta, h             )

    integer(kind=ip), intent(in) :: nx,ny  ! local dims
    integer(kind=ip), intent(in) :: npxg,npyg ! nb procs
    real(kind=rp)   , intent(in) :: Lx, Ly, Htot
    real(kind=rp), dimension(:,:), pointer, intent(out) :: dx, dy, zeta, h

    integer(kind=ip), parameter :: ip=4, rp=8
    integer(kind=ip):: nxg, nyg  ! global dims
    integer(kind=ip):: pi, pj
    integer(kind=ip):: i,j

    real(kind=rp) :: x, y
    real(kind=rp) :: x0, y0
    integer(kind=ip) :: ierr

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    nxg = npxg * nx
    nyg = npyg * ny

    ! grid definition

    dx(:,:) = Lx/real(nxg,kind=rp)
    dy(:,:) = Ly/real(nyg,kind=rp)

    pj = myrank/npxg
    pi = myrank-pj*npxg

    x0 = Lx * 0.5_rp
    y0 = Ly * 0.5_rp
    do i = 0,nx+1 !!!  I need to know my global index range
       do j = 0,ny+1
          x = (real(i+(pi*nx),kind=rp)-0.5_rp) * dx(i,j)
          y = (real(j+(pj*ny),kind=rp)-0.5_rp) * dy(i,j)
          zeta(i,j) = 0._rp
          h(i,j) = Htot * (1._rp - 0.5_rp * exp(-(x-x0)**2._rp/(Lx/5._rp)**2._rp -(y-y0)**2._rp/(Ly/5._rp)**2._rp))
          ! h(i,j) = Htot
       enddo
    enddo

    if (netcdf_output) then
       call write_netcdf(dx,vname='dx',netcdf_file_name='dx.nc',rank=myrank)
       call write_netcdf(dy,vname='dy',netcdf_file_name='dy.nc',rank=myrank)
       call write_netcdf(zeta, vname= 'zeta',netcdf_file_name= 'zeta.nc',rank=myrank)
       call write_netcdf(h, vname= 'h',netcdf_file_name= 'h.nc',rank=myrank)
    endif

  end subroutine setup_seamount

  !-------------------------------------------------------------------------     
  subroutine setup_rndtopo( &
       nx, ny, npxg, npyg,   &
       Lx, Ly, Htot,         &
       dx, dy,               &
       zeta, h             )

    integer(kind=ip), intent(in) :: nx,ny  ! local dims
    integer(kind=ip), intent(in) :: npxg,npyg ! nb procs
    real(kind=rp)   , intent(in) :: Lx, Ly, Htot
    real(kind=rp), dimension(:,:), pointer, intent(out) :: dx, dy, zeta, h

    integer(kind=ip), parameter :: ip=4, rp=8
    integer(kind=ip):: nxg, nyg  ! global dims
    integer(kind=ip):: pi, pj
    integer(kind=ip):: i,j

    real(kind=rp) :: x, y
    real(kind=rp) :: x0, y0
    integer(kind=ip) :: ierr

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    nxg = npxg * nx
    nyg = npyg * ny

    ! grid definition

    dx(:,:) = Lx/real(nxg,kind=rp)
    dy(:,:) = Ly/real(nyg,kind=rp)

    pj = myrank/npxg
    pi = myrank-pj*npxg

    call random_number(h(:,:)) ! between 0 and 1

    x0 = Lx * 0.5_rp
    y0 = Ly * 0.5_rp
    do i = 0,nx+1 !!!  I need to know my global index range
       do j = 0,ny+1
          x = (real(i+(pi*nx),kind=rp)-0.5_rp) * dx(i,j)
          y = (real(j+(pj*ny),kind=rp)-0.5_rp) * dy(i,j)
          zeta(i,j) = 0._rp
          h(i,j) = Htot * (20._rp/100._rp) * abs(h(i,j)) ! between 0% and 20% of Htot
!          h(i,j) = Htot * (1._rp - 0.5_rp * exp(-(x-x0)**2._rp/(Lx/5._rp)**2._rp -(y-y0)**2._rp/(Ly/5._rp)**2._rp))
       enddo
    enddo

    if (netcdf_output) then
       call write_netcdf(dx,vname='dx',netcdf_file_name='dx.nc',rank=myrank)
       call write_netcdf(dy,vname='dy',netcdf_file_name='dy.nc',rank=myrank)
       call write_netcdf(zeta, vname= 'zeta',netcdf_file_name= 'zeta.nc',rank=myrank)
       call write_netcdf(h, vname= 'h',netcdf_file_name= 'h.nc',rank=myrank)
    endif

  end subroutine setup_rndtopo

end module mg_setup_tests
