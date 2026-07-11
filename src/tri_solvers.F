      module tri_solvers

      ! uses a generic tri-diagonal solved algorithm to solve
      ! for vertical mixing and implicit vertical advection

      use dimensions
      use ocean_vars
      use tracers
      use grid
      use scalars

      implicit none
 
      private

      real,allocatable,dimension(:,:) :: a,b,c,d
      real,allocatable,dimension(:,:) :: cf
      real,allocatable,dimension(:,:) :: gam
      real,allocatable,dimension(:)   :: bet_inv

      public init_tri_solve
      public pre_t_solve

      contains

      
!----------------------------------------------------------------------
      subroutine init_tri_solve ![

      allocate(cf(nx,0:nz))
      cf = 0
      allocate(gam(nx,nz))
      allocate(bet_inv(nx))
      allocate(a(nx,nz))
      allocate(b(nx,nz))
      allocate(c(nx,nz))
      allocate(d(nx,nz))

      end subroutine init_tri_solve !]
!----------------------------------------------------------------------
      subroutine pre_t_solve(t,dt,dz,Av,om,j,drag)  ![
      implicit none

      ! Solve:
      !  dz*T(k) - dt*( 2*Ak(k  )*(T(k+1)-T(k  ))/(dz(k+1)+dz(k  ))
      !          + dt*( 2*Ak(k-1)*(T(k  )-T(k-1))/(dz(k  )+dz(k-1))
      !          + om(k  )*T(k+1) (om(k  )<0) / T(k)   (om(k)>0)
      !          - om(k-1)*T(k  ) (om(k-1)<0) / T(k-1) (om(k-1)>0)
      !    = dz*T(k)
      !  or: a(k)*T(k-1) + b(k)*T(k) + c(k)*T(k+1) = d(k)

      ! import/export
      real,intent(in) :: dt
      real,dimension(1-bf:nx+bf,1-bf:ny+bf,nz),intent(inout) :: t
      real,dimension(nx,nz),intent(in) :: dz
      real,dimension(nx,0:nz),intent(in) :: om,Av
      real,dimension(nx),intent(in),optional :: drag
      integer,intent(in) :: j

      ! local
      integer :: i,k

      do k=1,nz-1
        do i=1,nx
          cf(i,k)=2.*dt*Av(i,k)/(dz(i,k+1)+dz(i,k))
        enddo
      enddo

      do k=1,nz
        do i=1,nx
          a(i,k) = -cf(i,k-1)  - max(om(i,k-1),0.)
          b(i,k) =  dz(i,k) +cf(i,k)+cf(i,k-1)
     &                      + max(om(i,k),0.) - min(om(i,k-1),0.)
          c(i,k) = -cf(i,k) + min(om(i,k),0.)
          d(i,k) = t(i,j,k);
        enddo
      enddo

      if (present(drag)) then
        b(:,1) = b(:,1) + drag
      endif

      call tri_solve(a,b,c,d)

      do k=1,nz
        do i=1,nx
          t(i,j,k) = d(i,k)
        enddo
      enddo


      end subroutine pre_t_solve  !]
!----------------------------------------------------------------------
      subroutine tri_solve(a,b,c,d)  ![
      implicit none

      ! import/export
      real,    intent(in)    :: a(nx,nz), b(nx,nz), c(nx,nz)
      real,    intent(inout) :: d(nx,nz)

      ! local
      integer :: i,k

      do i=1,nx
        bet_inv(i) = 1./b(i,1)
        d(i,1)     = d(i,1)*bet_inv(i)                 !! x(1) = r(1)/bet
      enddo

      do k=2,nz
        do i=1,nx
          gam(i,k)   = c(i,k-1)*bet_inv(i)             !! c(k-1)/bet
          bet_inv(i) = 1./( b(i,k) - a(i,k)*gam(i,k) ) !! 1./(b(k) - a(k)*gam(k))
          d(i,k)     = ( d(i,k) - a(i,k)*d(i,k-1) )*bet_inv(i)
        enddo
      enddo

      do k=nz-1,1,-1
        do i=1,nx
          d(i,k) = d(i,k) - gam(i,k+1)*d(i,k+1)        !! x(k) = x(k) - gam(k+1)*x(k+1)
        enddo
      enddo

      end subroutine tri_solve  !]
!----------------------------------------------------------------------

      end module tri_solvers
