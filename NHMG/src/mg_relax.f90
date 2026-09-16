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

             if (ghost_bottom) then

                call relax_3D_8_heart_g(lev,p,b,cA,i,j,nz)

             else

                call relax_3D_8_heart(p,b,cA,i,j,nz)

             endif

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

                if (ghost_bottom) then

                   call relax_3D_8_heart_g(lev,p,b,cA,i,j,nz)

                else

                   call relax_3D_8_heart(p,b,cA,i,j,nz)

                endif

             enddo
          enddo
          call fill_halo(lev,p)
          if (ghost_bottom) call fill_halo(lev,grid(lev)%pg)
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

                   if (ghost_bottom) then

                      call relax_3D_8_heart_g(lev,p,b,cA,i,j,nz)

                   else

                      call relax_3D_8_heart(p,b,cA,i,j,nz)

                   endif

                enddo
             enddo
             call fill_halo(lev,p)
             if (ghost_bottom) call fill_halo(lev,grid(lev)%pg)
          if (ghost_bottom) call fill_halo(lev,grid(lev)%pg)
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
  subroutine relax_3D_8_heart_g(lev,p,b,cA,i,j,nz)
    ! ghost-bottom twin of relax_3D_8_heart: the column holds the ghost
    ! row (index 0) below the nz cells.  The ghost couples to its own
    ! cell 1 through cA(2,1) and to the four neighbouring cells 1 through
    ! cAg(3..6); the cell-1 row couples to the four neighbouring ghosts
    ! through cA(5,1), cA(8,1) and the neighbours' cAg(3), cAg(5).  No
    ! bottom xy pairs.
    integer(kind=ip), intent(in) :: lev
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip), intent(in) :: i,j,nz
    integer(kind=ip) :: k, n
    real(kind=rp), dimension(0:nz) :: rhs, d, ud, x
    real(kind=rp), dimension(:,:), pointer :: pg,bg
    real(kind=rp), dimension(:,:,:), pointer :: cAg

    pg => grid(lev)%pg;  bg => grid(lev)%bg;  cAg => grid(lev)%cAg
    n = nz+1

    ! ghost row (bottom face)
    rhs(0) = bg(j,i) - cAg(3,j,i)*p(1,j-1,i) - cAg(4,j,i)*p(1,j+1,i) &
                     - cAg(5,j,i)*p(1,j,i-1) - cAg(6,j,i)*p(1,j,i+1)
    d(0)  = cAg(1,j,i)
    ud(0) = cA(2,1,j,i)

    k=1 !lower level
    rhs(k) = b(k,j,i)                                              &
         - cA(3,k,j,i)*p(k+1,j-1,i)                                &
         - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
         - cA(5,k+1,j+1,i)*p(k+1,j+1,i)                            &
         - cA(6,k,j,i)*p(k+1,j,i-1)                                &
         - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
         - cA(8,k+1,j,i+1)*p(k+1,j,i+1)                            &
         - cA(5,k,j,i)*pg(j-1,i) - cAg(3,j+1,i)*pg(j+1,i)          &
         - cA(8,k,j,i)*pg(j,i-1) - cAg(5,j,i+1)*pg(j,i+1)
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
    ud(k)  = zero

    call tridiag(n,d,ud,rhs,x)
    pg(j,i)     = x(0)
    p(1:nz,j,i) = x(1:nz)

  end subroutine relax_3D_8_heart_g
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
  subroutine compute_residual(lev,res,after_rb)
    integer(kind=ip), intent(in) :: lev
    real(kind=rp)   , intent(out):: res
    logical, optional, intent(in) :: after_rb   ! p just relaxed by the RB line sweep

    real(kind=rp),dimension(:,:,:)  , pointer:: p
    real(kind=rp),dimension(:,:,:)  , pointer:: b
    real(kind=rp),dimension(:,:,:)  , pointer:: r
    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nd
    real(kind=rp) ::resloc
    logical :: half

    half = .false.
    if (present(after_rb)) half = after_rb .and. &
         (trim(relax_method)=='RB' .or. trim(relax_method)=='Red-Black')

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

       if (ghost_bottom) then
          call compute_residual_3D_8_g(lev,res,p,b,r,cA,nx,ny,nz,half,   &
               grid(lev)%neighb(4)==MPI_PROC_NULL, grid(lev)%neighb(2)==MPI_PROC_NULL, &
               grid(lev)%neighb(1)==MPI_PROC_NULL, grid(lev)%neighb(3)==MPI_PROC_NULL)
       else
       call compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz,half,          &
            grid(lev)%neighb(4)==MPI_PROC_NULL, grid(lev)%neighb(2)==MPI_PROC_NULL, &
            grid(lev)%neighb(1)==MPI_PROC_NULL, grid(lev)%neighb(3)==MPI_PROC_NULL)
       endif
    end if

    ! no halo fill of r: the restriction reads the interior only and the
    ! gather path refreshes the coarse halos itself

    if (lev < 2) then
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
  subroutine compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz,half,west,east,south,north)

    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz
    logical                                  , intent(in)   :: half
    logical                                  , intent(in)   :: west,east,south,north
    ! half: p was just relaxed by relax_3D_8_RB, whose second colour
    ! (i+j odd) solves each column exactly with the first colour fixed:
    ! those columns have zero residual except at k = 1, where the bottom
    ! cross terms couple same-colour columns updated in the same half
    ! sweep, and except next to a physical wall (west..north: this rank's
    ! walls), where the mirrored halo (coarse levels) changes with the
    ! column itself.  Elsewhere only the k = 1 row is evaluated there
    ! (about half a pass).

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

          if (half .and. mod(i+j,2) == 1 .and. .not. (         &
               (west .and. i==1) .or. (east .and. i==nx) .or.  &
               (south .and. j==1) .or. (north .and. j==ny) )) then
             r(2:nz,j,i) = zero
             cycle
          endif

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

  !----------------------------------------
  subroutine compute_residual_3D_8_g(lev,res,p,b,r,cA,nx,ny,nz,half,west,east,south,north)
    ! ghost-bottom twin of compute_residual_3D_8: the k=1 rows couple to
    ! the ghost rows instead of the bottom xy pairs, and the ghost rows
    ! (grid%rg) are part of the residual.  After the RB line sweep the
    ! second-colour columns are exact in every row, ghost included (no
    ! same-colour couplings remain), except next to a physical wall.
    integer(kind=ip), intent(in) :: lev
    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz
    logical                                  , intent(in)   :: half
    logical                                  , intent(in)   :: west,east,south,north
    real(kind=rp), dimension(:,:), pointer :: pg,bg,rg
    real(kind=rp), dimension(:,:,:), pointer :: cAg
    integer(kind=ip)           :: i,j,k

    pg => grid(lev)%pg;  bg => grid(lev)%bg;  rg => grid(lev)%rg;  cAg => grid(lev)%cAg
    res = zero

    do i = 1,nx
       do j = 1,ny

          if (half .and. mod(i+j,2) == 1 .and. .not. (         &
               (west .and. i==1) .or. (east .and. i==nx) .or.  &
               (south .and. j==1) .or. (north .and. j==ny) )) then
             r(1:nz,j,i) = zero;  rg(j,i) = zero
             cycle
          endif

          rg(j,i) = bg(j,i) - cAg(1,j,i)*pg(j,i) - cA(2,1,j,i)*p(1,j,i)  &
               - cAg(3,j,i)*p(1,j-1,i) - cAg(4,j,i)*p(1,j+1,i)           &
               - cAg(5,j,i)*p(1,j,i-1) - cAg(6,j,i)*p(1,j,i+1)
          res = res + rg(j,i)*rg(j,i)

          k=1 !lower level
          r(k,j,i) = b(k,j,i)                                           &
               - cA(1,k,j,i)*p(k,j,i)                                   &
               - cA(2,k,j,i)*pg(j,i)                                    &
               - cA(2,k+1,j,i)*p(k+1,j,i)                               &
               - cA(3,k,j,i)*p(k+1,j-1,i)                               &
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(5,k+1,j+1,i)*p(k+1,j+1,i)                           &
               - cA(6,k,j,i)*p(k+1,j,i-1)                               &
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
               - cA(8,k+1,j,i+1)*p(k+1,j,i+1)                           &
               - cA(5,k,j,i)*pg(j-1,i) - cAg(3,j+1,i)*pg(j+1,i)         &
               - cA(8,k,j,i)*pg(j,i-1) - cAg(5,j,i+1)*pg(j,i+1)
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

  end subroutine compute_residual_3D_8_g

end module mg_relax
