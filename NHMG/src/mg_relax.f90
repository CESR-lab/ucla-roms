module mg_relax

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange

  implicit none

contains

  !----------------------------------------
  subroutine relax(lev,nsweeps)
    integer(kind=ip), intent(in):: lev
    integer(kind=ip), intent(in):: nsweeps

    real(kind=rp),dimension(:,:,:), pointer:: p
    real(kind=rp),dimension(:,:,:), pointer:: b
    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nd

    p  => grid(lev)%p
    b  => grid(lev)%b
    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then
!      Be aware that the 2D mg is not tested yet
!      The definition of the matrix in 2D has not been implemented yet

       call relax_2D_5(lev,p,b,cA,nsweeps,nx,ny)

    else

       select case(trim(relax_method))

       case('Gauss-Seidel','GS')
          call relax_3D_8_GS(lev,p,b,cA,nsweeps,nx,ny,nz)

       case('Red-Black','RB')
          call relax_3D_8_RB(lev,p,b,cA,nsweeps,nx,ny,nz)

       case('Four-Color','FC')
          call relax_3D_8_FC(lev,p,b,cA,nsweeps,nx,ny,nz)

       end select

    end if

  end subroutine relax

  !----------------------------------------
  subroutine relax_2D_5(lev,p,b,cA,nsweeps,nx,ny)

    integer(kind=ip)                         , intent(in)   :: lev
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip)           :: i,j,k, it,rb
    integer(kind=ip)            :: ib,ie,jb,je,rbb,rbe,rbi
    real(kind=rp) :: z,gamma,g1,g2
    logical :: red_black=.true.

    gamma = one
    g1 = gamma
    g2 = one - gamma

    k=1

    do it = 1,nsweeps
       if (mod(it,1) == 0) then
          ib = 1
          ie = nx
          jb = 1
          je = ny
       else
          ib = 0
          ie = nx+1
          jb = 0
          je = ny+1
       endif

       if (red_black) then
          rbb = 1
          rbe = 2
          rbi = 2
       else
          rbb = 0
          rbe = 0
          rbi = 1
       endif

       do rb = rbb,rbe
          do i = ib,ie
             do j = jb+mod(i+rb,rbi),je,rbi

                z =    b(k,j,i)                                           &
                     - cA(2,k,j,i)*p(k,j-1,i  ) - cA(2,k,j+1,i  )*p(k,j+1,  i)&
                     - cA(3,k,j,i)*p(k,j  ,i-1) - cA(3,k,j  ,i+1)*p(k,j  ,i+1)&
                     - cA(4,k,j,i)*p(k,j-1,i-1) - cA(4,k,j+1,i+1)*p(k,j+1,i+1)&
                     - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1)

                p(k,j,i) = z / cA(1,k,j,i)

             enddo
          enddo
       enddo

       call fill_halo(lev,p)

    enddo

  end subroutine relax_2D_5

 !----------------------------------------
  subroutine relax_3D_8_GS(lev,p,b,cA,nsweeps,nx,ny,nz)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip)            :: i,j,it

    call tic(lev,'relax_3D_8_GS')

    ! add a loop on smoothing
    do it = 1,nsweeps

       do i = 1, nx
          do j = 1, ny

             call relax_3D_8_heart(p,b,cA,i,j,nz)

          enddo ! j
       enddo    ! i

       call fill_halo(lev,p)

    enddo  !it

    call toc(lev,'relax_3D_8_GS')

  end subroutine relax_3D_8_GS

 !----------------------------------------
  subroutine relax_3D_8_RB(lev,p,b,cA,nsweeps,nx,ny,nz)
! This is the recommended relaxation method for production

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip) :: i,j,it
    integer(kind=ip) :: rb

    call tic(lev,'relax_3D_8_RB')

    do it = 1,nsweeps
       do rb = 1, 2 ! Red black loop
          do i = 1, nx
             do j = 1+mod(i+rb,2),ny,2

                call relax_3D_8_heart(p,b,cA,i,j,nz)

             enddo
          enddo
          call fill_halo(lev,p)
       enddo
    enddo

    call toc(lev,'relax_3D_8_RB')

  end subroutine relax_3D_8_RB

 !----------------------------------------
  subroutine relax_3D_8_FC(lev,p,b,cA,nsweeps,nx,ny,nz)
! This is the recommended relaxation method for MPI debugging with topography

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip)            :: i,j,it
    integer(kind=ip)            :: fc1,fc2

    call tic(lev,'relax_3D_8_FC')

    do it = 1,nsweeps
       do fc1 = 1, 2
          do fc2 = 1, 2
             do i = 1 + mod(fc1-1,2), nx, 2
                do j = 1 + mod(fc2-1,2), ny, 2

                   call relax_3D_8_heart(p,b,cA,i,j,nz)

                enddo
             enddo
             call fill_halo(lev,p)
          enddo
       enddo
    enddo

    call toc(lev,'relax_3D_8_FC')

  end subroutine relax_3D_8_FC

  !----------------------------------------
  subroutine relax_3D_8_heart(p,b,cA,i,j,nz)

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip)                          , intent(in)  :: i
    integer(kind=ip)                          , intent(in)  :: j
    integer(kind=ip)                          , intent(in)  :: nz

    !- Local -!
    integer(kind=ip) :: k
    real(kind=rp), dimension(nz) :: rhs, d, ud

    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    k=1 !lower level
    rhs(k) = b(k,j,i)                                              &
         - cA(3,k,j,i)*p(k+1,j-1,i)                                &
         - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
         - cA(5,k+1,j+1,i)*p(k+1,j+1,i)                            &
         - cA(6,k,j,i)*p(k+1,j,i-1)                                &
         - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
         - cA(8,k+1,j,i+1)*p(k+1,j,i+1)                            &
!!  Special cross terms
         - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1) &
         - cA(8,k,j,i)*p(k,j-1,i-1) - cA(8,k,j+1,i+1)*p(k,j+1,i+1)

    d(k)   = cA(1,k,j,i)
    ud(k)  = cA(2,k+1,j,i)

    do k = 2,nz-1 !interior levels
       rhs(k) = b(k,j,i) &
            - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i) &
            - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
            - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i) &
            - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1) &
            - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
            - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1)
       d(k)   = cA(1,k,j,i)
       ud(k)  = cA(2,k+1,j,i)
    enddo

    k=nz !upper level
    rhs(k) = b(k,j,i)                                              &
         - cA(3,k-1,j+1,i)*p(k-1,j+1,i) &
         - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
         - cA(5,k,j,i)*p(k-1,j-1,i)                                &
         - cA(6,k-1,j,i+1)*p(k-1,j,i+1) &
         - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
         - cA(8,k,j,i)*p(k-1,j,i-1)
    d(k)   = cA(1,k,j,i)

    call tridiag(nz,d,ud,rhs,p(:,j,i)) ! solve 1D tridiagonal system

  end subroutine relax_3D_8_heart

  !----------------------------------------
  subroutine tridiag(l,d,dd,b,xc)
    !     Axc = b
    !     Solve tridiagonal system

    integer                   ,intent(in)  :: l
    real(kind=rp),dimension(l),intent(in)  :: d,b
    real(kind=rp),dimension(l),intent(in)  :: dd
    real(kind=rp),dimension(l),intent(out) :: xc

    integer                  :: k
    real(kind=rp),dimension(l):: gam
    real(kind=rp)             :: bet

    bet   = one / d(1)
    xc(1) = b(1)*bet
    do k=2,l
       gam(k)= dd(k-1)*bet
       bet     = one /(d(k)-dd(k-1)*gam(k))
       xc(k) = (b(k)-dd(k-1)*xc(k-1))*bet
    enddo
    do k=l-1,1,-1
       xc(k) = xc(k)-gam(k+1)*xc(k+1)
    enddo

  end subroutine tridiag

  !----------------------------------------
  subroutine compute_residual(lev,res)
    integer(kind=ip), intent(in) :: lev
    real(kind=rp)   , intent(out):: res

    real(kind=rp),dimension(:,:,:)  , pointer:: p
    real(kind=rp),dimension(:,:,:)  , pointer:: b
    real(kind=rp),dimension(:,:,:)  , pointer:: r
    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nd
    real(kind=rp) ::resloc

    p  => grid(lev)%p
    b  => grid(lev)%b
    r  => grid(lev)%r
    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then

       call compute_residual_2D_5(res,p,b,r,cA,nx,ny)

    else

       call compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)

    end if

    call fill_halo(lev,r)

    if (lev >-1) then
       resloc=res
       call global_sum(lev,resloc,res)
       res = sqrt(res)
    else
       res = -999._rp
    endif

  end subroutine compute_residual

  !----------------------------------------
  subroutine compute_residual_2D_5(res,p,b,r,cA,nx,ny)

    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip) :: i,j,k
    real(kind=rp)  :: z

    res = zero

    k=1

    do i = 1,nx
       do j = 1,ny

          z = b(k,j,i) - cA(1,k,j,i)*p(k,j,i)                           &
               - cA(2,k,j,i)*p(k,j-1,i  ) - cA(2,k,j+1,i  )*p(k,j+1,  i)&
               - cA(3,k,j,i)*p(k,j  ,i-1) - cA(3,k,j  ,i+1)*p(k,j  ,i+1)&
               - cA(4,k,j,i)*p(k,j-1,i-1) - cA(4,k,j+1,i+1)*p(k,j+1,i+1)&
               - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1)

          r(k,j,i) = z
          !          res = max(res,abs(r(k,j,i)))
          res = res+z*z

       enddo
    enddo

  end subroutine compute_residual_2D_5

  !----------------------------------------
  subroutine compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)

    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip)           :: i,j,k

    res = zero

    do i = 1,nx
       do j = 1,ny

          k=1 !lower level
          r(k,j,i) = b(k,j,i)                                           &
               - cA(1,k,j,i)*p(k,j,i)                                   &
               - cA(2,k+1,j,i)*p(k+1,j,i)                               &
               - cA(3,k,j,i)*p(k+1,j-1,i)                               &
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(5,k+1,j+1,i)*p(k+1,j+1,i)                           &
               - cA(6,k,j,i)*p(k+1,j,i-1)                               &
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
               - cA(8,k+1,j,i+1)*p(k+1,j,i+1)                           &
!!  Special cross terms
               - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1)&
               - cA(8,k,j,i)*p(k,j-1,i-1) - cA(8,k,j+1,i+1)*p(k,j+1,i+1)

          res = res+r(k,j,i)*r(k,j,i)

          do k = 2,nz-1 !interior levels
             r(k,j,i) = b(k,j,i)                                           &
                  - cA(1,k,j,i)*p(k,j,i)                                   &
                  - cA(2,k,j,i)*p(k-1,j,i)   - cA(2,k+1,j,i)*p(k+1,j,i)    &
                  - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i)&
                  - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
                  - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i)&
                  - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1)&
                  - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
                  - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1)

             res = res+r(k,j,i)*r(k,j,i)
          enddo

          k=nz !upper level
          r(k,j,i) = b(k,j,i)                                           &
               - cA(1,k,j,i)*p(k,j,i)                                   &
               - cA(2,k,j,i)*p(k-1,j,i)                                 &
               - cA(3,k-1,j+1,i)*p(k-1,j+1,i)                           &
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(5,k,j,i)*p(k-1,j-1,i)                               &
               - cA(6,k-1,j,i+1)*p(k-1,j,i+1)                           &
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
               - cA(8,k,j,i)*p(k-1,j,i-1)

          res = res+r(k,j,i)*r(k,j,i)

       enddo
    enddo

  end subroutine compute_residual_3D_8

end module mg_relax
