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
      public vmix_solve
      public vmix_solve_w

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
      subroutine vmix_solve(t,dt,dz,Av,om,j,drag)  ![
!DIR$ ATTRIBUTES FORCEINLINE :: vmix_solve
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
      cf(:,0) = 0

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
        b(:,1) = b(:,1) + dt*drag
      endif

      call tri_solve(a,b,c,d)

      do k=1,nz
        do i=1,nx
          t(i,j,k) = d(i,k)
        enddo
      enddo

      end subroutine vmix_solve  !]
!----------------------------------------------------------------------
      subroutine vmix_solve_w(w,dt,dz,Av,om,j)  ![
!DIR$ ATTRIBUTES FORCEINLINE :: vmix_solve_w
      implicit none

      ! Solve:
      !  0.5*(dz(k)+dz(k+1))*w(k) - dt*( 0.5*(Av(k)+Av(k+1))*(w(k+1)-w(k  ))/dz(k+1)
      !                + dt*( 0.5*(Av(k)+Av(k-1))*(w(k  )-w(k-1))/dz(k  )
      !                + 0.5*(om(k)+om(k+1))*w(k+1) (om<0) / w(k)  (om>0)
      !                - 0.5*(om(k)+om(k-1))*w(k  ) (om<0) / w(k-1)(om>0)
      !    = 0.5*(dz(k)+dz(k+1))*w(k)
      !  w(0) = 0; This means we can just ignore a(1) and leave b(1) and c(1) as is.
      !  0.5*dz(nz)*w(nz) -
      !          + dt*( 0.5*(Av(nz)+Av(nz-1))*(w(nz)-w(nz-1))/dz(nz )
      !          - 0.5*(om(nz-1))*w(nz  ) (om<0) / w(nz-1)(om>0)
      !  = 0.5+dz(nz)*w(nz)
      !
      !  or: a(k)*w(k-1) + b(k)*w(k) + c(k)*w(k+1) = w(k)
      ! with:
      !    a(k) = -dt*0.5*(Av(k)+Av(k-1)) /dz(k  ) - max(0.5*(om(k)+om(k-1),0)
      !    b(k) = -dt*0.5*(Av(k)+Av(k-1)) /dz(k  ) - min(0.5*(om(k)+om(k-1),0)
      !           +dt*0.5*(Av(k)+Av(k+1)) /dz(k+1) - max(0.5*(om(k)+om(k+1),0)
      !    c(k) = +dt*0.5*(Av(k)+Av(k+1)) /dz(k+1) - min(0.5*(om(k)+om(k+1),0)

      ! import/export
      real,intent(in) :: dt
      real,dimension(1-bf:nx+bf,1-bf:ny+bf,0:nz),intent(inout) :: w
      real,dimension(nx,nz),intent(in) :: dz
      real,dimension(nx,0:nz),intent(in) :: om,Av
      integer,intent(in) :: j

      ! local
      integer :: i,k

      do k=0,nz-1
        do i=1,nx
          cf(i,k)=dt*0.5*(Av(i,k)+Av(i,k+1))/dz(i,k+1)
        enddo
      enddo

      do k=1,nz-1
        do i=1,nx
          a(i,k) = -cf(i,k-1)- max(0.5*(om(i,k)+om(i,k-1)),0.)
          b(i,k) = 0.5*(dz(i,k)+dz(i,k+1))
     &                   +cf(i,k)+cf(i,k-1)
     &                       + max(0.5*(om(i,k)+om(i,k+1)),0.)
     &                       - min(0.5*(om(i,k)+om(i,k-1)),0.)
          c(i,k) = -cf(i,k  )+ min(0.5*(om(i,k)+om(i,k+1)),0.)
          d(i,k) = w(i,j,k);
        enddo
      enddo
      do i=1,nx
        a(i,nz) = -cf(i,nz-1)- max(0.5*(om(i,nz)+om(i,nz-1)),0.)
        b(i,nz) = 0.5*dz(i,nz)
     &                +cf(i,nz)+cf(i,nz-1)
     &                       - min(0.5*(om(i,nz)+om(i,nz-1)),0.)
        d(i,nz) = w(i,j,nz);
      enddo

      call tri_solve(a,b,c,d)

      do k=1,nz
        do i=1,nx
          w(i,j,k) = d(i,k)
        enddo
      enddo


      end subroutine vmix_solve_w  !]
!----------------------------------------------------------------------
      subroutine tri_solve(a,b,c,d)  ![
!DIR$ ATTRIBUTES FORCEINLINE :: tri_solve
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
