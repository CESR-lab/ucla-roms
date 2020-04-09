module mg_gather

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids

  implicit none

  interface gather
     module procedure &
          gather_2D,  &
          gather_3D
  end interface gather

contains
  !----------------------------------------
  subroutine gather_2D(lev,x,y)

    ! lev is the level of y, after the gathering
    ! the grid information for x is not available because x does not belong to a grid
    ! x is a temporary array (see mg_intergrids.f90)

    integer(kind=ip),intent(in) :: lev
    real(kind=rp),dimension(:,:),pointer,intent(in) :: x
    real(kind=rp),dimension(:,:),pointer,intent(out) :: y

    integer(kind=ip):: nx,ny
    integer(kind=ip):: ngx,ngy,Ng
    integer(kind=ip):: i,j,l,m,ii,jj
    integer(kind=ip):: i0,i1
    integer(kind=ip):: ierr
    real(kind=rp),dimension(:,:,:,:),pointer :: buffer

    ! number of cores per direction involved in this gathering (1 or 2)
    ngx = grid(lev)%ngx
    ngy = grid(lev)%ngy

    nx = grid(lev)%nx
    ny = grid(lev)%ny

    ! numel(x)
    Ng = grid(lev)%Ng2D

    buffer => grid(lev)%gatherbuffer2D

    call MPI_ALLGATHER( x, Ng, MPI_DOUBLE_PRECISION, buffer, Ng, MPI_DOUBLE_PRECISION, grid(lev)%localcomm,ierr)

    nx = nx / ngx
    ny = ny / ngy
    do m=0,ngy-1
       ! copy only the inner points of x into y because 
       ! the halo of x is corrupted
       ! Indeed, x comes from the coarsening
       ! after which we didn't update the halo
       ! because fill_halo is not available for the intermediate grid...
       !
       do l=0,ngx-1
          !
          if(l==0)then
             i0=0
          else
             i0=1
          endif
          if(l==ngx-1)then
             i1=nx+1
          else
             i1=nx
          endif
          !
          ii = 1+l*nx
          do i=1,nx
             jj = 1+m*ny
             do j=1,ny
                y(jj,ii) = buffer(j,i,l,m)
                jj=jj+1
             enddo
             ii=ii+1
          enddo
       enddo
    enddo

  end subroutine gather_2D
  !----------------------------------------
  subroutine gather_3D(lev,x,y)

    ! lev is the level of y, after the gathering

    ! the grid information for x is not available because x does not belong to a grid
    ! x is a temporary array (see mg_intergrids.f90)

    integer(kind=ip),intent(in) :: lev
    real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
    real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y

    integer(kind=ip):: nx,ny,nz
    integer(kind=ip):: ngx,ngy,Ng
    integer(kind=ip):: i,j,k,l,m,ii,jj
    integer(kind=ip):: i0,i1
    integer(kind=ip):: ierr
    real(kind=rp),dimension(:,:,:,:,:),pointer :: buffer

    !  ngx = 1,2 the number of subdomains to be gathered in x
    !  ngy = 1,2 the number of subdomains to be gathered in y

    ngx = grid(lev)%ngx
    ngy = grid(lev)%ngy

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz

    Ng = grid(lev)%Ng
    buffer => grid(lev)%gatherbuffer

    call MPI_ALLGATHER( x, Ng, MPI_DOUBLE_PRECISION, buffer, Ng, MPI_DOUBLE_PRECISION, grid(lev)%localcomm,ierr)


    nx = nx / ngx
    ny = ny / ngy

    do m=0,ngy-1
       ! copy only the inner points of x into y because 
       ! the halo of x is corrupted
       ! Indeed, x comes from the coarsening
       ! after which we didn't update the halo
       ! because fill_halo is not available for the intermediate grid...
       !
       do l=0,ngx-1
          !
          if(l==0)then
             i0=0
          else
             i0=1
          endif
          if(l==ngx-1)then
             i1=nx+1
          else
             i1=nx
          endif
          !
          ii = 1+l*nx
          do i=1,nx
             jj = 1+m*ny
             do j=1,ny
                do k=1,nz
                   y(k,jj,ii) = buffer(k,j,i,l,m)
                enddo
                jj=jj+1
             enddo
             ii=ii+1
          enddo
       enddo
    enddo

  end subroutine gather_3D

  !----------------------------------------
  subroutine split(lev,x,y)

    ! lev is the level of x, where it has to be split
    ! y is the dummy 3D intermediate array, before interpolation

    integer(kind=ip),intent(in) :: lev
    real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
    real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y

    integer(kind=ip):: nx,ny,nz
    integer(kind=ip):: ngx,ngy
    integer(kind=ip):: i,j,k,l,m,ii,jj,key

    !  ngx = 1,2 the number of subdomains to be gathered in x
    !  ngy = 1,2 the number of subdomains to be gathered in y
    ngx = grid(lev)%ngx
    ngy = grid(lev)%ngy

    nx = grid(lev)%nx / ngx
    ny = grid(lev)%ny / ngy
    nz = grid(lev)%nz

!! TODO this requires a comment, but Guillaume does not remember what it should be
    key = grid(lev)%key

    l = mod(key,2)
    m = key/2

    ii = l*nx

    do i= 0,nx+1

       jj = 0+m*ny
       do j= 0,ny+1

          do k=1,nz
             y(k,j,i) = x(k,jj,ii)
          enddo
          jj=jj+1

       enddo
       ii=ii+1

    enddo

  end subroutine split

end module mg_gather
