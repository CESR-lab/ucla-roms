module mg_diagnostics

  use mg_cst
  use mg_mpi
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine diag_momentum
    
    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:)  , pointer :: Arz
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer :: um,vm,wm

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
   
    dx    => grid(1)%dx
    dy    => grid(1)%dy
    dxu   => grid(1)%dxu
    dyv   => grid(1)%dyv
    dzw   => grid(1)%dzw
    Arx   => grid(1)%Arx
    Ary   => grid(1)%Ary
    Arz   => grid(1)%Arz
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx

    u => grid(1)%u
    v => grid(1)%v
    w => grid(1)%w
    um => grid(1)%um
    vm => grid(1)%vm
    wm => grid(1)%wm

    ! wm

!   do i = 1,nx
!      do j = 1,ny
!         k = 1 
!         wm(k,j,i) = 1/Arz(j,i) * (w(k,j,i) + hlf*( &
!              +zx(k  ,j,i)*u(k  ,j,i  ) &
!              +zx(k  ,j,i)*u(k  ,j,i+1) &
!              +zy(k  ,j,i)*v(k  ,j  ,i) &
!              +zy(k  ,j,i)*v(k  ,j+1,i) ))
!         do k = 2,nz 
!            wm(k,j,i) = 1/Arz(j,i) * (w(k,j,i) + qrt*( &
!                  +zx(k-1,j,i)*u(k-1,j,i  ) &
!                  +zx(k-1,j,i)*u(k-1,j,i+1) &
!                  +zx(k  ,j,i)*u(k  ,j,i  ) &
!                  +zx(k  ,j,i)*u(k  ,j,i+1) &
!                  +zy(k-1,j,i)*v(k-1,j  ,i) &
!                  +zy(k-1,j,i)*v(k-1,j+1,i) &
!                  +zy(k  ,j,i)*v(k  ,j  ,i) &
!                  +zy(k  ,j,i)*v(k  ,j+1,i) ))
!         enddo
!         k = nz+1
!         wm(k,j,i) = 1/Arz(j,i) * (w(k,j,i) + hlf*( &
!              +zx(k-1,j,i)*u(k-1,j,i  ) &
!              +zx(k-1,j,i)*u(k-1,j,i+1) &
!              +zy(k-1,j,i)*v(k-1,j  ,i) &
!              +zy(k-1,j,i)*v(k-1,j+1,i) )) 
!      enddo
!   enddo

    call fill_halo(1,wm)

    ! um and vm 

    do i = 1,nx+1  
       do j = 0,ny+1
          k = 1
          um(k,j,i) = 1/Arx(k,j,i) * (u(k,j,i) + qrt*( &
               + zxdy(k,j,i  )* two *dzw(k  ,j,i  )*wm(k  ,j,i  ) & 
               + zxdy(k,j,i  )      *dzw(k+1,j,i  )*wm(k+1,j,i  ) &
               + zxdy(k,j,i-1)* two *dzw(k  ,j,i-1)*wm(k  ,j,i-1) & 
               + zxdy(k,j,i-1)      *dzw(k+1,j,i-1)*wm(k+1,j,i-1) ))
          do k = 2,nz-1
             um(k,j,i) = 1/Arx(k,j,i) * (u(k,j,i) + qrt*( &
                  + zxdy(k,j,i  ) *dzw(k  ,j,i  )*wm(k  ,j,i  ) & 
                  + zxdy(k,j,i  ) *dzw(k+1,j,i  )*wm(k+1,j,i  ) &
                  + zxdy(k,j,i-1) *dzw(k  ,j,i-1)*wm(k  ,j,i-1) & 
                  + zxdy(k,j,i-1) *dzw(k+1,j,i-1)*wm(k+1,j,i-1) ))
          enddo
          k = nz
          um(k,j,i) = 1/Arx(k,j,i) * (u(k,j,i) + qrt*( &
               + zxdy(k,j,i  )      *dzw(k  ,j,i  )*wm(k  ,j,i  ) & 
               + zxdy(k,j,i  )* two *dzw(k+1,j,i  )*wm(k+1,j,i  ) &
               + zxdy(k,j,i-1)      *dzw(k  ,j,i-1)*wm(k  ,j,i-1) & 
               + zxdy(k,j,i-1)* two *dzw(k+1,j,i-1)*wm(k+1,j,i-1) ))
       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1
          k = 1
          vm(k,j,i) = 1/Ary(k,j,i) * (v(k,j,i) + qrt*( &
               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i)*wm(k  ,j  ,i) & 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i)*wm(k+1,j  ,i) & 
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i)*wm(k  ,j-1,i) & 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i)*wm(k+1,j-1,i) ))
          do k = 2,nz-1
             vm(k,j,i) = 1/Ary(k,j,i) * (v(k,j,i) + qrt*( &
                  + zydx(k,j  ,i)* dzw(k  ,j  ,i)*wm(k  ,j  ,i) & 
                  + zydx(k,j  ,i)* dzw(k+1,j  ,i)*wm(k+1,j  ,i) &  
                  + zydx(k,j-1,i)* dzw(k  ,j-1,i)*wm(k  ,j-1,i) & 
                  + zydx(k,j-1,i)* dzw(k+1,j-1,i)*wm(k+1,j-1,i) ))
          enddo
          k = nz
          vm(k,j,i) = 1/Ary(k,j,i) * (v(k,j,i) + qrt*( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i)*wm(k  ,j  ,i) & 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*wm(k+1,j  ,i) & 
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*wm(k  ,j-1,i) & 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*wm(k+1,j-1,i) )) 
       enddo
    enddo
    
  end subroutine diag_momentum

  !-----------------------------------------------------------------------------------
  subroutine diag_kin_energy

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer :: um,vm,wm
    real(kind=rp), dimension(:,:,:), pointer :: ke

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
   
    u => grid(1)%u
    v => grid(1)%v
    w => grid(1)%w
    um => grid(1)%um
    vm => grid(1)%vm
    wm => grid(1)%wm

    ke => grid(1)%ke

    do i = 1,nx
       do j = 1,ny
          do k = 1,nz 
             ke(k,j,i) = qrt *( & 
                +  (u(k,j,i  )*um(k,j,i  ) + u(k,j,i+1)*um(k,j,i+1)) &
                +  (v(k,j  ,i)*vm(k,j  ,i) + v(k,j+1,i)*vm(k,j+1,i)) &
                +  (w(k  ,j,i)*wm(k  ,j,i) + w(k+1,j,i)*wm(k+1,j,i)) )
          enddo 
       enddo
    enddo

  end subroutine diag_kin_energy

end module mg_diagnostics

