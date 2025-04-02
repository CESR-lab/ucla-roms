module mg_intergrids

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather

  implicit none

contains

  !------------------------------------------------------------
  subroutine fine2coarse(lev)

    ! coarsen grid(lev)%r to grid(lev+1)%b and set grid(lev+1)%p=0

    integer(kind=ip), intent(in) :: lev

    real(kind=rp),dimension(:,:,:),pointer :: r
    real(kind=rp),dimension(:,:,:),pointer :: b

    integer(kind=ip) :: nx, ny, nz

    nx = grid(lev+1)%nx
    ny = grid(lev+1)%ny
    nz = grid(lev+1)%nz

    r => grid(lev)%r

    if (grid(lev+1)%gather == 1) then
       !- b here is sub level rhs, we will gather at the end
       !- in the level+1 rhs
       b => grid(lev+1)%dummy3
       !- ngx and ngy are the number of processes being gathered in x and y
       !- values are 1 or 2
       nx = grid(lev+1)%nx / grid(lev+1)%ngx
       ny = grid(lev+1)%ny / grid(lev+1)%ngy
    else
       b => grid(lev+1)%b
    endif

    if (grid(lev)%nz == 1) then
       call fine2coarse_2D(r,b,nx,ny)

    else
       call tic(lev,'fine2coarse_3D')

       call fine2coarse_3D(r,b,nx,ny,nz)

       call toc(lev,'fine2coarse_3D')

    end if

    if (grid(lev+1)%gather == 1) then
       r => grid(lev+1)%dummy3
       b => grid(lev+1)%b
       call gather(lev+1,r,b)
    endif

    call fill_halo(lev+1,grid(lev+1)%b)

    grid(lev+1)%p = zero

  end subroutine fine2coarse

  !---------------------------------------------------------
  subroutine fine2coarse_2D(x,y,nx,ny)
!  be aware that the 2D mg is not tested yet

    real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
!   real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y
    real(kind=rp),dimension(:,:,:),pointer :: y
    integer(kind=ip), intent(in) :: nx, ny

    integer(kind=ip) :: i,j,i2,j2
    integer(kind=ip) :: d

    d = size(x,1) ! vertical dimension of the fine level, can 2 or 1

   if(d==1)then
       ! x was already 2D
       do i2=1,nx
          i=2*i2-1
          do j2=1,ny
             j=2*j2-1      
             y(1,j2,i2) = (x(1,j,i)+x(1,j,i+1)+x(1,j+1,i)+x(1,j+1,i+1))
          enddo
       enddo
    else
       ! x was 3D
       do i2=1,nx
          i=2*i2-1
          do j2=1,ny
             j=2*j2-1     
             y(1,j2,i2) = (x(1,j,i)+x(1,j,i+1)+x(1,j+1,i)+x(1,j+1,i+1)&
                          +x(2,j,i)+x(2,j,i+1)+x(2,j+1,i)+x(2,j+1,i+1))
          enddo
       enddo
    endif

  end subroutine fine2coarse_2D

  !------------------------------------------------------------
  subroutine fine2coarse_3D(x,y,nx,ny,nz)
    !
    ! Fine2coarse 'x' from fine level l1 to 'y' on coarse level l2=l1+1
    real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
!   real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y
    real(kind=rp),dimension(:,:,:),pointer :: y
    integer(kind=ip), intent(in) :: nx, ny, nz
    ! local
    integer(kind=ip) :: i,j,k,i2,j2,k2
    real(kind=rp):: z

    !! No division by 8 because we work in volume integrated form
    !! See Brink et al., 'Tutorial on Multi Grid'
    do i2=1,nx
       i=2*i2-1
       do j2=1,ny
          j=2*j2-1
          do k2=1,nz
             k=2*k2-1
             z = x(k  ,j,i)+x(k  ,j,i+1)+x(k  ,j+1,i)+x(k  ,j+1,i+1) &
               + x(k+1,j,i)+x(k+1,j,i+1)+x(k+1,j+1,i)+x(k+1,j+1,i+1)
             y(k2,j2,i2) = z
          enddo
       enddo
    enddo

  end subroutine fine2coarse_3D

  !------------------------------------------------------------
  subroutine coarse2fine(lev)

    ! interpolate grid(lev+1)%p to grid(lev)%r and add it to grid(lev)%p

    integer(kind=ip), intent(in) :: lev

    real(kind=rp),dimension(:,:,:),pointer :: rf
    real(kind=rp),dimension(:,:,:),pointer :: pc

    integer(kind=ip) :: nxc, nyc, nzc

    nxc = grid(lev+1)%nx
    nyc = grid(lev+1)%ny
    nzc = grid(lev+1)%nz

    pc => grid(lev+1)%p

    if (grid(lev+1)%gather == 1) then
       rf => grid(lev+1)%dummy3
       call split(lev+1,pc,rf)
       pc => grid(lev+1)%dummy3
       nxc = grid(lev+1)%nx / grid(lev+1)%ngx
       nyc = grid(lev+1)%ny / grid(lev+1)%ngy
       nzc = grid(lev+1)%nz       
    endif

    rf => grid(lev)%r

    if (grid(lev)%nz == 1) then

       call coarse2fine_2D_linear(rf,pc,nxc,nyc)

    else

       call coarse2fine_3D_linear(rf,pc,nxc,nyc,nzc)

    endif

    call fill_halo(lev,grid(lev)%r)

    grid(lev)%p = grid(lev)%p + grid(lev)%r

  end subroutine coarse2fine

  !------------------------------------------------------------
  subroutine coarse2fine_2D_linear(xf,xc,nx,ny)

!   real(kind=rp),dimension(:,:,:),pointer,intent(out) :: xf
    real(kind=rp),dimension(:,:,:),pointer :: xf
    real(kind=rp),dimension(:,:,:),pointer,intent(in)  :: xc
    integer(kind=ip),intent(in) :: nx, ny

    ! local
  integer(kind=ip) :: i,j,i2,j2,k,k2
  real(kind=rp) :: a,b,c
  !
  ! weights for bilinear in (i,j)
  a = 9._8 / 16._8
  b = 3._8 / 16._8
  c = 1._8 / 16._8

  k  = 1
  k2 = 1

  do i2=1,nx
     i=2*i2-1
     do j2=1,ny
        j=2*j2-1
        xf(k  ,j  ,i  ) =  &
             + a * xc(k2,j2  ,i2) + c * xc(k2,j2-1,i2-1) &
             + b * xc(k2,j2-1,i2) + b * xc(k2,j2  ,i2-1)
        xf(k  ,j+1,i  ) =  &
             + a * xc(k2,j2  ,i2) + c * xc(k2,j2+1,i2-1) &
             + b * xc(k2,j2+1,i2) + b * xc(k2,j2  ,i2-1)
        xf(k  ,j  ,i+1) =  &
             + a * xc(k2,j2  ,i2) + c * xc(k2,j2-1,i2+1) &
             + b * xc(k2,j2-1,i2) + b * xc(k2,j2  ,i2+1)
        xf(k  ,j+1,i+1) =  &
             + a * xc(k2,j2  ,i2) + c * xc(k2,j2+1,i2+1) &
             + b * xc(k2,j2+1,i2) + b * xc(k2,j2  ,i2+1)
     enddo
  enddo


  end subroutine coarse2fine_2D_linear

  !------------------------------------------------------------
  subroutine coarse2fine_3D_linear(xf,xc,nx,ny,nz)

!   real(kind=rp),dimension(:,:,:),pointer,intent(out) :: xf
    real(kind=rp),dimension(:,:,:) :: xf
    real(kind=rp),dimension(:,:,:),pointer,intent(in)  :: xc
    integer(kind=ip),intent(in) :: nx, ny, nz

    ! local
    integer(kind=ip) :: i,j,k,i2,j2,k2,kp
    real(kind=rp) :: a,b,c,d,e,f,g

    integer(kind=rp) :: dirichlet_flag

    if (surface_neumann) then
       dirichlet_flag = 0
    else
       dirichlet_flag = 1
    endif

    !
    ! weights for bilinear in (i,j), nearest in k
    a = 9._8 / 16._8
    b = 3._8 / 16._8
    c = 1._8 / 16._8
    !
    ! weights for trilinear in (i,j,k)
    d = 27._8 / 64._8
    e =  9._8 / 64._8
    f =  3._8 / 64._8
    g =  1._8 / 64._8
    ! 
    do i2=1,nx
       i=2*i2-1
       do j2=1,ny
          j=2*j2-1

          ! bottom level
          k  = 1
          k2 = 1
          xf(k  ,j  ,i  ) =  &
               + a * xc(k2,j2  ,i2) + c * xc(k2,j2-1,i2-1) &
               + b * xc(k2,j2-1,i2) + b * xc(k2,j2  ,i2-1)
          xf(k  ,j+1,i  ) =  &
               + a * xc(k2,j2  ,i2) + c * xc(k2,j2+1,i2-1) &
               + b * xc(k2,j2+1,i2) + b * xc(k2,j2  ,i2-1)
          xf(k  ,j  ,i+1) =  &
               + a * xc(k2,j2  ,i2) + c * xc(k2,j2-1,i2+1) &
               + b * xc(k2,j2-1,i2) + b * xc(k2,j2  ,i2+1)
          xf(k  ,j+1,i+1) =  &
               + a * xc(k2,j2  ,i2) + c * xc(k2,j2+1,i2+1) &
               + b * xc(k2,j2+1,i2) + b * xc(k2,j2  ,i2+1)

          ! interior level
          do k=2,nz*2-1
             k2 = ((k+1)/2)
             kp = k2-(mod(k,2)*2-1)
             xf(k  ,j  ,i  ) =  &
                  + d * xc(k2,j2  ,i2) + f * xc(k2,j2-1,i2-1) &
                  + e * xc(k2,j2-1,i2) + e * xc(k2,j2  ,i2-1) &
                  + e * xc(kp,j2  ,i2) + g * xc(kp,j2-1,i2-1) &
                  + f * xc(kp,j2-1,i2) + f * xc(kp,j2  ,i2-1)
             xf(k  ,j+1,i  ) =  &
                  + d * xc(k2,j2  ,i2) + f * xc(k2,j2+1,i2-1) &
                  + e * xc(k2,j2+1,i2) + e * xc(k2,j2  ,i2-1) &
                  + e * xc(kp,j2  ,i2) + g * xc(kp,j2+1,i2-1) &
                  + f * xc(kp,j2+1,i2) + f * xc(kp,j2  ,i2-1)
             xf(k  ,j  ,i+1) = &
                  + d * xc(k2,j2  ,i2) + f * xc(k2,j2-1,i2+1) &
                  + e * xc(k2,j2-1,i2) + e * xc(k2,j2  ,i2+1) &
                  + e * xc(kp,j2  ,i2) + g * xc(kp,j2-1,i2+1) &
                  + f * xc(kp,j2-1,i2) + f * xc(kp,j2  ,i2+1)
             xf(k  ,j+1,i+1) = &
                  + d * xc(k2,j2  ,i2) + f * xc(k2,j2+1,i2+1) &
                  + e * xc(k2,j2+1,i2) + e * xc(k2,j2  ,i2+1) &
                  + e * xc(kp,j2  ,i2) + g * xc(kp,j2+1,i2+1) &
                  + f * xc(kp,j2+1,i2) + f * xc(kp,j2  ,i2+1)
          enddo

          ! top level
          k = nz*2
          xf(k  ,j  ,i  ) =  (1-hlf*dirichlet_flag)  * (     &
               a * xc(k2,j2  ,i2) + c * xc(k2,j2-1,i2-1) +   &
               b * xc(k2,j2-1,i2) + b * xc(k2,j2  ,i2-1)   )
          xf(k  ,j+1,i  ) =  (1-hlf*dirichlet_flag) * (      &
               a * xc(k2,j2  ,i2) + c * xc(k2,j2+1,i2-1) +   &
               b * xc(k2,j2+1,i2) + b * xc(k2,j2  ,i2-1)   )
          xf(k  ,j  ,i+1) =  (1-hlf*dirichlet_flag)  * (     &
               a * xc(k2,j2  ,i2) + c * xc(k2,j2-1,i2+1) +   &
               b * xc(k2,j2-1,i2) + b * xc(k2,j2  ,i2+1)   ) 
          xf(k  ,j+1,i+1) =  (1-hlf*dirichlet_flag)  * (     &
               a * xc(k2,j2  ,i2) + c * xc(k2,j2+1,i2+1) +   &
               b * xc(k2,j2+1,i2) + b * xc(k2,j2  ,i2+1)   )
       enddo
    enddo

  end subroutine coarse2fine_3D_linear

end module mg_intergrids
