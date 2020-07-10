module nhmg

  use mg_cst
  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_tictoc
  use mg_mpi_exchange
  use mg_autotune
  use mg_horiz_grids
  use mg_vert_grids
  use mg_projection
  use mg_solvers
  use mg_netcdf_out

  implicit none

  private

    integer(kind=ip) :: tscount = 1
    integer(kind=ip),public,parameter :: halo = 2
    integer(kind=ip),public :: iprec1 = 1,iprec2 = 2

    public nhmg_matrices
    public nhmg_init
    public nhmg_clean
    public nhmg_solve

contains

  !--------------------------------------------------------------
  subroutine nhmg_init(nx,ny,nz,npxg,npyg)
      
    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg

    call tic(1,'nhmg_init')

    call mg_mpi_init()

    if (myrank==0) write(*,*)' nhmg_init:'

    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg,npyg,nx,ny,nz)

    call define_neighbours()

    call print_grids()

    call toc(1,'nhmg_init')

  end subroutine nhmg_init

  !--------------------------------------------------------------
  subroutine nhmg_matrices(nx,ny,nz,hl,pdx,pdy,zxa,zya,Hza,dxa,dya)

    integer(kind=ip), intent(in) :: nx,ny,nz
    integer(kind=ip), intent(in) :: hl,pdx,pdy

    real(kind=rp), dimension(1-hl:nx+hl+pdx,1-hl:ny+hl+pdy,1:nz),intent(in) :: zxa,zya
    real(kind=rp), dimension(1-hl:nx+hl+pdx,1-hl:ny+hl+pdy,1:nz),intent(in) :: Hza
    real(kind=rp), dimension(0:nx+1,0:ny+1)     , optional, intent(in) :: dxa
    real(kind=rp), dimension(0:nx+1,0:ny+1)     , optional, intent(in) :: dya

    integer(kind=ip) :: i,j,k

    call tic(1,'nhmg_matrices')

    if (present(dxa) .and. present(dya)) then

       do i = 0,nx+1
          do j = 0,ny+1
             grid(1)%dx(j,i) = dxa(i,j)
             grid(1)%dy(j,i) = dya(i,j)
          enddo
       enddo

       call set_horiz_grids()

    end if

    do i = -1,nx+2
       do j = -1,ny+2
          do k = 1,nz
             grid(1)%dz(k,j,i) = Hza(i,j,k)
          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1, nz
             grid(1)%zxdy(k,j,i) = zxa(i,j,k) * grid(1)%dy(j,i)
             grid(1)%zydx(k,j,i) = zya(i,j,k) * grid(1)%dx(j,i)
          enddo
       enddo
    enddo

    call fill_outer_halos(grid(1)%dz,nx,ny,nz)

    call set_vert_grids()

    call set_matrices()

    call toc(1,'nhmg_matrices')

  end subroutine nhmg_matrices

  !--------------------------------------------------------------
  subroutine fill_outer_halos(z,nx,ny,nz)
    real(kind=rp), dimension(1:nz,-1:ny+2,-1:nx+2), intent(inout) :: z

    integer :: i,j,k,nx,ny,nz

    if(grid(1)%neighb(1).eq.MPI_PROC_NULL) then
       ! south
       do i=-1,nx+2 ! <=extended range to fill the corners (not a bug)
          do k=1,nz
             z(k,-1,i)=z(k,1,i)
             z(k, 0,i)=z(k,1,i)
          enddo
       enddo
    endif
    if(grid(1)%neighb(3).eq.MPI_PROC_NULL) then
       ! north
       do i=-1,nx+2
          do k=1,nz
             z(k,ny+1,i)=z(k,ny,i)
             z(k,ny+2,i)=z(k,ny,i)
          enddo
       enddo
    endif
    if(grid(1)%neighb(2).eq.MPI_PROC_NULL) then
       ! east
       do j=-1,ny+2
          do k=1,nz
             z(k,j,nx+1)=z(k,j,nx)
             z(k,j,nx+2)=z(k,j,nx)
          enddo
       enddo
    endif
    if(grid(1)%neighb(4).eq.MPI_PROC_NULL) then
       ! west
       do j=-1,ny+2
          do k=1,nz
             z(k,j,-1)=z(k,j,1)
             z(k,j, 0)=z(k,j,1)
          enddo
       enddo
    endif

  end subroutine fill_outer_halos
  !--------------------------------------------------------------
  subroutine nhmg_solve(nx,ny,nz,hl,pdx,pdy,ua,va,wa)

    integer(kind=ip), intent(in) :: nx,ny,nz
    integer(kind=ip), intent(in) :: hl,pdx,pdy

    real(kind=rp), dimension(1-hl:nx+hl+pdx,1-hl:ny+hl+pdy,1:nz  ),intent(in) :: ua
    real(kind=rp), dimension(1-hl:nx+hl+pdx,1-hl:ny+hl+pdy,1:nz  ),intent(in) :: va
    real(kind=rp), dimension(1-hl:nx+hl+pdx,1-hl:ny+hl+pdy,1:nz+1),intent(in) :: wa

    integer(kind=ip) :: i,j,k

    call tic(1,'nhmg_solve')

    !  Fill the rhs for the poisson equation
    do i = 1,nx
       do j = 1,ny 
          do k = 1,nz
             grid(1)%b(k,j,i) = ua(i+1,j,k) - ua(i,j,k) &
                              + va(i,j+1,k) - va(i,j,k) &
                              + wa(i,j,k+1) - wa(i,j,k)
          enddo
       enddo
    enddo

    !- auto tuning tests if autotune = .true.
    if ((tscount == autotune_ts).and.(autotune)) then
       call sb_autotune()  !- test of autotuning
    end if

    call solve_p()   

    tscount = tscount + 1

    call correction_uvw()

!   do i = 1,nx
!      do j = 1,ny 
!         do k = 1,nz
!            grid(1)%b(k,j,i) = ( ua(i+1,j,k) + grid(1)%du(k,j,i+1) )  &
!                             - ( ua(i  ,j,k) + grid(1)%du(k,j,i  ) )  &
!                             + ( va(i,j+1,k) + grid(1)%dv(k,j+1,i) )  &
!                             - ( va(i,j  ,k) + grid(1)%dv(k,j  ,i) )  &
!                             + ( wa(i,j,k+1) + grid(1)%dw(k+1,j,i) )  &
!                             - ( wa(i,j,k  ) + grid(1)%dw(k  ,j,i) )
!         enddo
!      enddo
!   enddo

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='so.nc',rank=myrank,iter=1)
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='so.nc',rank=myrank,iter=1)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='so.nc',rank=myrank,iter=1)
       call write_netcdf(grid(1)%du,vname='du',netcdf_file_name='so.nc',rank=myrank,iter=1)
       call write_netcdf(grid(1)%dv,vname='dv',netcdf_file_name='so.nc',rank=myrank,iter=1)
       call write_netcdf(grid(1)%dw,vname='dw',netcdf_file_name='so.nc',rank=myrank,iter=1)
    endif

    call toc(1,'nhmg_solve')

  end subroutine nhmg_solve

!!$  !--------------------------------------------------------------
!!$  subroutine nhmg_checkdivergence(ua,va,wa,Hza)
!!$!! Move this to the Croco side
!!$!! It's here to provide a template
!!$
!!$    real(kind=rp), dimension(:,:,:), intent(in) :: ua
!!$    real(kind=rp), dimension(:,:,:), intent(in) :: va
!!$    real(kind=rp), dimension(:,:,:), intent(in) :: wa
!!$    real(kind=rp), dimension(:,:,:), intent(in) :: Hza
!!$
!!$    real(kind=rp), dimension(:,:),   pointer :: dx,dy
!!$    real(kind=rp), dimension(:,:,:), pointer :: u,v,w,dz
!!$
!!$    integer(kind=ip) :: i,j,k,is,js,ishift
!!$    integer(kind=ip) :: nx,ny,nz
!!$
!!$    integer(kind=ip), save :: iter_checkdiv=0
!!$    iter_checkdiv = iter_checkdiv + 1
!!$
!!$
!!$    nx = grid(1)%nx
!!$    ny = grid(1)%ny
!!$    nz = grid(1)%nz
!!$
!!$    dx => grid(1)%dx
!!$    dy => grid(1)%dy
!!$
!!$    ishift=2
!!$
!!$    ! need to update dz because define_matrices may not be called every time step
!!$    dz => grid(1)%dz
!!$    do k=1,nz
!!$       do j=0,ny+1
!!$          js=j+ishift
!!$          do i=0,nx+1
!!$             is=i+ishift
!!$             dz(k,j,i) = Hza(is,js,k)
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    ! set fluxes
!!$    u  => grid(1)%u
!!$    v  => grid(1)%v
!!$    w  => grid(1)%w
!!$
!!$    do k=1,nz
!!$       do j=1,ny
!!$          js=j+ishift
!!$          do i=1,nx+1
!!$             is=i+ishift
!!$             u(k,j,i) = ua(is,js,k) * &
!!$                  qrt * (dz(k,j,i) + dz(k,j,i-1)) * (dy(j,i)+dy(j,i-1))
!!$          enddo
!!$       enddo
!!$       do j=1,ny+1
!!$          js=j+ishift
!!$          do i=1,nx
!!$             is=i+ishift
!!$             v(k,j,i) = va(is,js,k) * &
!!$                  qrt * (dz(k,j,i) + dz(k,j-1,i)) * (dx(j,i)+dx(j-1,i))
!!$          enddo
!!$       enddo
!!$       do j=1,ny
!!$          js=j+ishift
!!$          do i=1,nx
!!$             is=i+ishift
!!$             w(k+1,j,i) = wa(is,js,k+1) * dx(j,i) * dy(j,i)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    w(1,:,:) = zero
!!$
!!$!   call set_rhs()
!!$
!!$    write(*,*) 'check div',maxval(abs(grid(1)%b))
!!$
!!$    call write_netcdf(grid(1)%b,vname='rhs',netcdf_file_name='chk.nc',rank=myrank,iter=iter_checkdiv)
!!$
!!$  end subroutine nhmg_checkdivergence

  !--------------------------------------------------------------

  subroutine nhmg_clean()

    !   call grids_dealloc()

    if (myrank==0) then
       call print_tictoc(myrank)
    endif

  end subroutine nhmg_clean

end module nhmg

