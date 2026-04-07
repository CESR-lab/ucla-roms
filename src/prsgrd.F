      module pressure_grad

      use dimensions
      use diagnostics
      use tides
      use grid
      use ocean_vars
      use eos_vars

      use coupling
      use scalars
      use roms_mpi

      implicit none

      private

      ! These are the coefficients for the 6th order line integrals

      real,parameter :: C1 = 11./42.
      real,parameter :: C2 = 13./110.
      real,parameter :: C3 =  1./42.
      real,parameter :: C4 =  3./20.
      real,parameter :: C5 =  1./630.
      real,parameter :: C6 =  1./8.

      ! These are the coefficients for the 4th order compact schemes
      real,parameter :: alp1 = 0.25
      real,parameter :: alp2 = 0.10
!     real,dimension(:),allocatable :: cp1,cp2
      real,dimension(:,:),allocatable :: bdr,bdz
      real,dimension(:,:),allocatable :: bddr,bddz

      logical,parameter :: split_eos = .false.

      real,dimension(:,:,:),allocatable :: P
      real,dimension(:,:,:),allocatable :: dR,dZ,ddR,ddZ

      public prsgrd
      
      contains

!---------------------------------------------------------------------------------------
      subroutine ini_prsgrd  ![

      logical,save  :: initialized = .false.
      integer :: n

      if (initialized) return

      allocate(  P(1-bf:nx+bf,1-bf:ny+bf,nz) )
      allocate( dR(1-bf:nx+bf,1-bf:ny+bf,nz) )
      allocate( dZ(1-bf:nx+bf,1-bf:ny+bf,nz) )
      allocate(ddR(1-bf:nx+bf,1-bf:ny+bf,nz) )
      allocate(ddZ(1-bf:nx+bf,1-bf:ny+bf,nz) )
      dR = 0
      dZ = 0
      ddR = 0
      ddZ = 0

      n = max(ny,nz)
      allocate( bdr(0:nx,0:n) )
      allocate( bdz(0:nx,0:n) )
      allocate( bddr(0:nx,0:n) )
      allocate( bddz(0:nx,0:n) )

      initialized = .true.

      end subroutine ini_prsgrd  !]
!---------------------------------------------------------------------------------------
      subroutine prsgrd  ![

      implicit none

      integer  :: i,j,k
      real :: dpth
      real :: delz,delr
      real :: delzp,delrp
      real :: sumzp,sumrp
      real :: delzpp,delrpp
      real :: sumzpp,sumrpp

      real, parameter :: i24 = 1./24.

      real :: grho,HalfGRho


      call ini_prsgrd

      grho=g/rho0
      HalfGRho=0.5*grho

      ru = 0
      rv = 0

      call compact_gradients_z

      do j=0,ny
        do i= 0,nx
          P(i,j,nz)=g*z_w(i,j,nz) + grho*( rho(i,j,nz)
     &       +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &          /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
          if (pot_tides) then
            P(i,j,N) = P(i,j,N) - g*ptide(i,j)
          endif
        enddo
        do k=nz-1,1,-1
          do i=0,nx
            delz = z_r(i,j,k+1)-z_r(i,j,k) 
            delr = rho(i,j,k+1)-rho(i,j,k)
            delzp= dZ(i,j,k+1)-dZ(i,j,k)
            delrp= dR(i,j,k+1)-dR(i,j,k)

            P(i,j,k)=P(i,j,k+1)+HalfGRho*( 
     &            (rho(i,j,k+1)+rho(i,j,k))*delz
     &     -C1*( ( dR(i,j,k+1)- dR(i,j,k))*( delz - C2*( dZ(i,j,k+1)+ dZ(i,j,k)) )
     &          -( dZ(i,j,k+1)- dZ(i,j,k))*( delr - C2*( dR(i,j,k+1)+ dR(i,j,k)) ) )
     &     +C3*( (ddR(i,j,k+1)+ddR(i,j,k))*( delz - C4*( dZ(i,j,k+1)+ dZ(i,j,k)) )
     &          -(ddZ(i,j,k+1)+ddZ(i,j,k))*( delr - C4*( dR(i,j,k+1)+ dR(i,j,k)) ) )
     &     +C5*( (ddR(i,j,k+1)-ddR(i,j,k))*( delzp- C6*(ddZ(i,j,k+1)+ddZ(i,j,k)) )
     &          -(ddZ(i,j,k+1)-ddZ(i,j,k))*( delrp- C6*(ddR(i,j,k+1)+ddR(i,j,k)) ) )
     &                                                                      )
          enddo
        enddo
      enddo

      call compact_gradients_x

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            delz = z_r(i,j,k)-z_r(i-1,j,k) 
            delr = rho(i,j,k)-rho(i-1,j,k)
            delzp= dZ(i,j,k)-dZ(i-1,j,k)
            delrp= dR(i,j,k)-dR(i-1,j,k)
            sumzp= dZ(i,j,k)+dZ(i-1,j,k)
            sumrp= dR(i,j,k)+dR(i-1,j,k)
            delzpp = ddZ(i,j,k)-ddZ(i-1,j,k)
            delrpp = ddR(i,j,k)-ddR(i-1,j,k)
            sumzpp = ddZ(i,j,k)+ddZ(i-1,j,k)
            sumrpp = ddR(i,j,k)+ddR(i-1,j,k)

            ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &                              P(i-1,j,k)-P(i,j,k)-HalfGRho*(
     &            (rho(i,j,k)+rho(i-1,j,k))*delz
     &     -C1*( delrp *( delz - C2*sumzp ) - delzp *( delr - C2*sumrp ) )
     &     +C3*( sumrpp*( delz - C4*sumzp ) - sumzpp*( delr - C4*sumrp ) )
     &     +C5*( delrpp*( delzp- C6*sumzpp) - delzpp*( delrp- C6*sumrpp) )))
          enddo
        enddo
      enddo

      call compact_gradients_y

!     Udiag(1:nx,1:ny,:,ipgr) = dR(1:nx,1:ny,:)
!     Udiag(1:nx,1:ny,:,izxb) = dZ(1:nx,1:ny,:)
!     Vdiag(1:nx,1:ny,:,ipgr) = ddR(1:nx,1:ny,:)
!     Vdiag(1:nx,1:ny,:,izxb) = ddZ(1:nx,1:ny,:)

      do k = 1,nz
        do j= 1,ny
          do i=1,nx
            delz = z_r(i,j,k)-z_r(i,j-1,k) 
            delr = rho(i,j,k)-rho(i,j-1,k)
            delzp= dZ(i,j,k)-dZ(i,j-1,k)
            delrp= dR(i,j,k)-dR(i,j-1,k)
            sumzp= dZ(i,j,k)+dZ(i,j-1,k)
            sumrp= dR(i,j,k)+dR(i,j-1,k)
            delzpp = ddZ(i,j,k)-ddZ(i,j-1,k)
            delrpp = ddR(i,j,k)-ddR(i,j-1,k)
            sumzpp = ddZ(i,j,k)+ddZ(i,j-1,k)
            sumrpp = ddR(i,j,k)+ddR(i,j-1,k)

            rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &                              P(i,j-1,k)-P(i,j,k)-HalfGRho*(
     &            (rho(i,j,k)+rho(i,j-1,k))*delz
     &     -C1*( delrp *( delz - C2*sumzp ) - delzp *( delr - C2*sumrp ) )
     &     +C3*( sumrpp*( delz - C4*sumzp ) - sumzpp*( delr - C4*sumrp ) )
     &     +C5*( delrpp*( delzp- C6*sumzpp) - delzpp*( delrp- C6*sumrpp) )))
            enddo
        enddo 
      enddo


      if (nrhs==3) then ! corrector stage
        if ((calc_diag.and.diag_uv)) then
          do k=1,nz
            Udiag(1:nx,1:ny,k,ipgr) = ru(1:nx,1:ny,k)*dxdyi_u*umask(1:nx,1:ny)
            Vdiag(1:nx,1:ny,k,ipgr) = rv(1:nx,1:ny,k)*dxdyi_v*vmask(1:nx,1:ny)

            Udiag(1:nx,1:ny,k,izxb) = Udiag(1:nx,1:ny,k,ipgr) - 
     &        0.5*(Hz(1:nx,1:ny,k)+Hz(0:nx-1,1:ny,k))*dn_u(1:nx,1:ny)
     &          *(P(0:nx-1,1:ny,k)-P(1:nx,1:ny,k))*dxdyi_u(1:nx,1:ny)
            Vdiag(1:nx,1:ny,k,izxb) = Vdiag(1:nx,1:ny,k,ipgr) - 
     &        0.5*(Hz(1:nx,1:ny,k)+Hz(1:nx,0:ny-1,k))*dm_v(1:nx,1:ny)
     &          *(P(1:nx,0:ny-1,k)-P(1:nx,1:ny,k))*dxdyi_v(1:nx,1:ny)
          enddo
          do k=1,nz-1
            Wdiag(1:nx,1:ny,k,idzb) = (P(1:nx,1:ny,k)-P(1:nx,1:ny,k+1))*dxdyi
          enddo
        endif
        if (diag_pflx) then
          call calc_diag_pflux(p)
        endif
      endif


# ifdef NHMG
!
!--------------------------------------------------------------------
! hydrostatic pressure gradient for the vertical flux
!       dpdz = -zx*dpdx - zy*dpdy (See Molemaker et al., 2018)
!--------------------------------------------------------------------
!
!     We are re-using dR as a work array 

      do j=Jstr,Jend
         do k = 1,N
            do i=Istr,Iend
              dR(i,k) = dzdxi(i,j,k) *(                                
     $              ru(i  ,j,k  )/(Hz(i,j,k)+Hz(i-1,j,k)) 
     $            + ru(i+1,j,k  )/(Hz(i,j,k)+Hz(i+1,j,k)))
     $            +     dzdeta(i,j,k) *(                               
     $              rv(i,j  ,k  )/(Hz(i,j,k)+Hz(i,j-1,k)) 
     $            + rv(i,j+1,k  )/(Hz(i,j,k)+Hz(i,j+1,k)))
            enddo
         enddo

         !! rw is a V*dwdt units object (m4/s2), just like ru, rv
         do k = 1,N-1 
            do i=Istr,Iend
              rw(i,j,k)=-0.25*(dR(i,k+1)+dR(i,k))*(Hz(i,j,k+1)+Hz(i,j,k))
            enddo
         enddo
         do i=Istr,Iend
            rw(i,j,N) = -0.5*dR(i,N)*Hz(i,j,N)
         enddo
      enddo
      if (diag_w.and.calc_diag) then
        Wdiag(1:nx,1:ny,:,iwprsgr) =rw(1:nx,1:ny,:)*dxdyi
      endif
# endif /* NHMG */

      end subroutine prsgrd !]
!-----------------------------------------------------------------------
      subroutine compact_gradients_z  ![

      implicit none

      ! Computes vertical 1st and 2nd derivates using compact schemes

      integer :: i,j,k

      ! First derivative:
      !  0.25 f'(k-1) + f'(k) + 0.25 f'(k+1) = 0.75(f(k+1) - f(k-1))
      ! Second derivative:
      !  0.1 f''(k-1) + f''(k) + 0.1 f''(k+1) = 1.2( f(k+1)-2f(k))+f(k-1))
      ! Assume f' and f'' zero in the buffer of the physical domain


      ! Compute 2nd order gradients for k=1,k=nz and use as 
      ! boundary conditions for the tri-diagonal system
      do j=0,ny
        do i=0,nx
          if (split_eos) then
# ifdef SPLIT_EOS
            dpth= -0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dR(i,j,1) = 0.5*(-3*rho1(i,j, 1)+4*rho1(i,j,2)-rho1(i,j,3))
     &                + 0.5*(-3* qp1(i,j, 1)+4* qp1(i,j,2)- qp1(i,j,3))
     &                     *dpth*(1.-qp2*dpth) 
            dR(i,j,nz)= 0.5*(-3*rho1(i,j,nz)+4*rho1(i,j,nz-1)+rho1(i,j,nz-2))
     &                + 0.5*(-3*qp1(i,j,nz) +4*qp1(i,j,nz-1) +qp1(i,j,nz-2))
     &                     *dpth*(1.-qp2*dpth) 
# endif
          else
            dR(i,j, 1) = 0.5*(-3*rho(i,j, 1)+4*rho(i,j,   2)-rho(i,j,   3))
            dR(i,j,nz) = 0.5*( 3*rho(i,j,nz)-4*rho(i,j,nz-1)+rho(i,j,nz-2))
          endif
          dZ(i,j, 1) = 0.5*(-3*z_r(i,j, 1)+4*z_r(i,j,   2)-z_r(i,j,   3))
          dZ(i,j,nz) = 0.5*( 3*z_r(i,j,nz)-4*z_r(i,j,nz-1)+z_r(i,j,nz-2))
        enddo

        !  0.25 f'(k-1) + f'(k) + 0.25 f'(k+1) = 0.75(f(k+1) - f(k-1))
        do k=2,nz-1
          do i=0,nx
            if (split_eos) then
            else
              bdr(i,k) = 0.75*(rho(i,j,k+1)-rho(i,j,k-1))
              bdz(i,k) = 0.75*(z_r(i,j,k+1)-z_r(i,j,k-1))
            endif
          enddo
        enddo
        do i=0,nx
          bdr(i, 2)   = bdr(i, 2)   - 0.25*dR(i,j, 1)
          bdr(i,nz-1) = bdr(i,nz-1) - 0.25*dR(i,j,nz)
          bdz(i, 2)   = bdz(i, 2)   - 0.25*dZ(i,j, 1)
          bdz(i,nz-1) = bdz(i,nz-1) - 0.25*dZ(i,j,nz)

          call thomas_tds(nz-2,alp1,dR(i,j,2:nz-1),bdr(i,2:nz-1))
          call thomas_tds(nz-2,alp1,dZ(i,j,2:nz-1),bdz(i,2:nz-1))
        enddo

        !  0.1 f''(k-1) + f''(k) + 0.1 f''(k+1) = 1.2( f(k+1)-2f(k))+f(k-1))
        do i=0,nx
          if (split_eos) then
            print *,'not ready'
            stop
          else
            ddZ(i,j, 1) = 2*z_r(i,j, 1)-5*z_r(i,j,   2)+4*z_r(i,j,   3)-z_r(i,j,   4)
            ddZ(i,j,nz) = 2*z_r(i,j,nz)-5*z_r(i,j,nz-1)+4*z_r(i,j,nz-2)-z_r(i,j,nz-3)
            ddR(i,j, 1) = 2*rho(i,j, 1)-5*rho(i,j,2)+4*rho(i,j,3)-rho(i,j,4)
            ddR(i,j,nz) = 2*rho(i,j,nz)-5*rho(i,j,nz-1)+4*rho(i,j,nz-2)-rho(i,j,nz-3)
          endif
        enddo
        do k=2,nz-1
          do i=0,nx
            if (split_eos) then
            else
              bddr(i,k) = 1.2*(rho(i,j,k+1)-2*rho(i,j,k)+rho(i,j,k-1))
              bddz(i,k) = 1.2*(z_r(i,j,k+1)-2*z_r(i,j,k)+z_r(i,j,k-1))
            endif
          enddo
        enddo
        do i=0,nx
          bddr(i, 2)   = bddr(i, 2)   - 0.1*ddR(i,j, 1)
          bddr(i,nz-1) = bddr(i,nz-1) - 0.1*ddR(i,j,nz)
          bddz(i, 2)   = bddz(i, 2)   - 0.1*ddZ(i,j, 1)
          bddz(i,nz-1) = bddz(i,nz-1) - 0.1*ddZ(i,j,nz)

          call thomas_tds(nz-2,alp2,ddR(i,j,2:nz-1),bddr(i,2:nz-1))
          call thomas_tds(nz-2,alp2,ddZ(i,j,2:nz-1),bddz(i,2:nz-1))
        enddo
      enddo

      end subroutine compact_gradients_z  !]
!-----------------------------------------------------------------------
      subroutine compact_gradients_x  ![

      ! Computes horizontal 1st and 2nd derivates using compact schemes
      ! Works with distributed memory parallelization.
      ! cp1 and cp2 can pre-computed in prsgrd_ini

      ! First derivative:
      !  0.25 f'(i-1) + f'(i) + 0.25 f'(i+1) = 0.75(f(i+1) - f(i-1))
      ! Second derivative:
      !  0.1 f''(k-1) + f''(k) + 0.1 f''(k+1) = 1.2( f(k+1)-2f(k))+f(k-1))
      ! Assume f' and f'' zero in the buffer of the physical domain

      integer :: i,j,k
      real :: alp1i,alp2i,det1,det2
      real :: cf1,cf2,cor0
      real,dimension(nx) :: cp1,cp2
      real,dimension(nx) :: cr1_lft,cr2_lft
      real,dimension(nx) :: cr1_rgt,cr2_rgt
      real,dimension(ny,nz) :: cr_dR,cr_dZ
      real,dimension(ny,nz) :: cr_ddR,cr_ddZ

      alp1i = 1./alp1
      alp2i = 1./alp2

      ! compute cp1, cp2
      cp1(1) = alp1 
      cp2(1) = alp2 
      do i=2,nx
        cp1(i) = alp1 /(1.0 - alp1*cp1(i-1) )
        cp2(i) = alp2 /(1.0 - alp2*cp2(i-1) )
      enddo

      ! compute left and right correction functions
      cr1_lft(1) = 1
      cr2_lft(1) = 1
      do i = 2,nx
        cr1_lft(i) = - cr1_lft(i-1)*cp1(i);
        cr2_lft(i) = - cr2_lft(i-1)*cp2(i);
      enddo
      do i = nx-1,1,-1
        cr1_lft(i) = cr1_lft(i) - cp1(i) * cr1_lft(i+1);
        cr2_lft(i) = cr2_lft(i) - cp2(i) * cr2_lft(i+1);
      enddo

      ! the right hand correction functions are just mirrors
      do i = 1,nx
        cr1_rgt(nx+1-i) = cr1_lft(i)
        cr2_rgt(nx+1-i) = cr2_lft(i)
      enddo

      do k=1,nz
        do j=1,ny
          do i=1,nx
            if (split_eos) then
               stop 'split eos not ready'
            else
              bdr(i,j)  = 0.75*(rho(i+1,j,k)-rho(i-1,j,k));
              bddr(i,j) = 1.20*(rho(i+1,j,k)-2*rho(i,j,k)+rho(i-1,j,k))
            endif
              bdz(i,j)  = 0.75*(z_r(i+1,j,k)-z_r(i-1,j,k));
              bddz(i,j) = 1.20*(z_r(i+1,j,k)-2*z_r(i,j,k)+z_r(i-1,j,k))
          enddo
    
          !Solve locally at each sub-domain using a forward-backward sweep

          dR(1,j,k)  = bdr(1,j)
          dZ(1,j,k)  = bdz(1,j)
          ddR(1,j,k) = bddr(1,j)
          ddZ(1,j,k) = bddz(1,j)
          do i = 2,nx
            dR(i,j,k) = (bdr(i,j)*alp1i - dR(i-1,j,k))*cp1(i);
            dZ(i,j,k) = (bdz(i,j)*alp1i - dZ(i-1,j,k))*cp1(i);
            ddR(i,j,k) = (bddr(i,j)*alp2i - ddR(i-1,j,k))*cp2(i);
            ddZ(i,j,k) = (bddz(i,j)*alp2i - ddZ(i-1,j,k))*cp2(i);
          enddo
          do i = nx-1,1,-1
            dR(i,j,k) = dR(i,j,k)   - cp1(i)*dR(i+1,j,k);
            dZ(i,j,k) = dZ(i,j,k)   - cp1(i)*dZ(i+1,j,k);
            ddR(i,j,k) = ddR(i,j,k) - cp2(i)*ddR(i+1,j,k);
            ddZ(i,j,k) = ddZ(i,j,k) - cp2(i)*ddZ(i+1,j,k);
          enddo

        enddo
      enddo

      call exchange_xxx(dR,dZ,ddR,ddZ)

      cf1 = alp1*cr1_lft(1)
      det1 = 1 - cf1*cf1
      cf2 = alp2*cr2_lft(1)
      det2 = 1 - cf2*cf2

      if (inode>0) then
        do k=1,nz
          do j=1,ny
            cr_dR(j,k) = ( dR(0,j,k) - cf1* dR(1,j,k) )/det1
            cr_dZ(j,k) = ( dZ(0,j,k) - cf1* dZ(1,j,k) )/det1
            cr_ddR(j,k)= (ddR(0,j,k) - cf2*ddR(1,j,k) )/det2
            cr_ddZ(j,k)= (ddZ(0,j,k) - cf2*ddZ(1,j,k) )/det2
            ! We need corrected gradients at i=0 as well and have all the information
            cor0 = ( dR(1,j,k) - cf1* dR(0,j,k) )/det1
            dR(0,j,k) =  dR(0,j,k) - alp1* cor0*cr1_rgt(nx)
            cor0 = ( dZ(1,j,k) - cf1* dZ(0,j,k) )/det1
            dZ(0,j,k) =  dZ(0,j,k) - alp1* cor0*cr1_rgt(nx)
            cor0 = ( ddR(1,j,k) - cf2*ddR(0,j,k) )/det2
            ddR(0,j,k)=  ddR(0,j,k) - alp2* cor0*cr2_rgt(nx)
            cor0 = ( ddZ(1,j,k) - cf2* ddZ(0,j,k) )/det2
            ddZ(0,j,k)=  ddZ(0,j,k) - alp2* cor0*cr2_rgt(nx)
            do i=1,nx
              dR(i,j,k) =  dR(i,j,k) - alp1* cr_dR(j,k)*cr1_lft(i)
              dZ(i,j,k) =  dZ(i,j,k) - alp1* cr_dZ(j,k)*cr1_lft(i)
              ddR(i,j,k)= ddR(i,j,k) - alp2*cr_ddR(j,k)*cr2_lft(i)
              ddZ(i,j,k)= ddZ(i,j,k) - alp2*cr_ddZ(j,k)*cr2_lft(i)
            enddo
          enddo
        enddo
      endif

      if (inode<npx-1) then
        do k=1,nz
          do j=1,ny
            cr_dR(j,k) = ( dR(nx+1,j,k) - cf1* dR(nx,j,k) )/det1
            cr_dZ(j,k) = ( dZ(nx+1,j,k) - cf1* dZ(nx,j,k) )/det1
            cr_ddR(j,k)= (ddR(nx+1,j,k) - cf2*ddR(nx,j,k) )/det2
            cr_ddZ(j,k)= (ddZ(nx+1,j,k) - cf2*ddZ(nx,j,k) )/det2
            do i=1,nx
              dR(i,j,k) =  dR(i,j,k) - alp1* cr_dR(j,k)*cr1_rgt(i)
              dZ(i,j,k) =  dZ(i,j,k) - alp1* cr_dZ(j,k)*cr1_rgt(i)
              ddR(i,j,k)= ddR(i,j,k) - alp2*cr_ddR(j,k)*cr2_rgt(i)
              ddZ(i,j,k)= ddZ(i,j,k) - alp2*cr_ddZ(j,k)*cr2_rgt(i)
            enddo
          enddo
        enddo
      endif

      end subroutine compact_gradients_x !]
!------------------------------------------
      subroutine compact_gradients_y  ![

      ! Computes horizontal 1st and 2nd derivates using compact schemes
      ! Works with distributed memory parallelization.
      ! cp1 and cp2 can pre-computed in prsgrd_ini

      ! First derivative:
      !  0.25 f'(i-1) + f'(i) + 0.25 f'(i+1) = 0.75(f(i+1) - f(i-1))
      ! Second derivative:
      !  0.1 f''(k-1) + f''(k) + 0.1 f''(k+1) = 1.2( f(k+1)-2f(k))+f(k-1))
      ! Assume f' and f'' zero in the buffer of the physical domain

      integer :: i,j,k
      real :: alp1i,alp2i,det1,det2
      real :: cf1,cf2,cor0
      real,dimension(ny) :: cp1,cp2
      real,dimension(ny) :: cr1_lft,cr2_lft
      real,dimension(ny) :: cr1_rgt,cr2_rgt
      real,dimension(nx,nz) :: cr_dR,cr_dZ
      real,dimension(nx,nz) :: cr_ddR,cr_ddZ

      alp1i = 1./alp1
      alp2i = 1./alp2

      ! compute cp1, cp2
      cp1(1) = alp1 
      cp2(1) = alp2 
      do j=2,ny
        cp1(j) = alp1 /(1.0 - alp1*cp1(j-1) )
        cp2(j) = alp2 /(1.0 - alp2*cp2(j-1) )
      enddo

      ! compute left and right correction functions
      cr1_lft(1) = 1
      cr2_lft(1) = 1
      do j = 2,ny
        cr1_lft(j) = - cr1_lft(j-1)*cp1(j);
        cr2_lft(j) = - cr2_lft(j-1)*cp2(j);
      enddo
      do j = ny-1,1,-1
        cr1_lft(j) = cr1_lft(j) - cp1(j) * cr1_lft(j+1);
        cr2_lft(j) = cr2_lft(j) - cp2(j) * cr2_lft(j+1);
      enddo

      ! the right hand correction functions are just mirrors
      do j = 1,ny
        cr1_rgt(ny+1-j) = cr1_lft(j)
        cr2_rgt(ny+1-j) = cr2_lft(j)
      enddo

      do k=1,nz
        do i=1,nx
          do j=1,ny
            if (split_eos) then
               stop 'split eos not ready'
            else
              bdr(i,j)  = 0.75*(rho(i,j+1,k)-rho(i,j-1,k));
              bddr(i,j) = 1.20*(rho(i,j+1,k)-2*rho(i,j,k)+rho(i,j-1,k))
            endif
            bdz(i,j)  = 0.75*(z_r(i,j+1,k)-z_r(i,j-1,k));
            bddz(i,j) = 1.20*(z_r(i,j+1,k)-2*z_r(i,j,k)+z_r(i,j-1,k))
          enddo
    
          !Solve locally at each sub-domain using a forward-backward sweep

          dR(i,1,k)  = bdr(i,1)
          dZ(i,1,k)  = bdz(i,1)
          ddR(i,1,k) = bddr(i,1)
          ddZ(i,1,k) = bddz(i,1)
          do j = 2,ny
            dR(i,j,k) = (bdr(i,j)*alp1i - dR(i,j-1,k))*cp1(j);
            dZ(i,j,k) = (bdz(i,j)*alp1i - dZ(i,j-1,k))*cp1(j);
            ddR(i,j,k) = (bddr(i,j)*alp2i - ddR(i,j-1,k))*cp2(j);
            ddZ(i,j,k) = (bddz(i,j)*alp2i - ddZ(i,j-1,k))*cp2(j);
          enddo
          do j = ny-1,1,-1
            dR(i,j,k) = dR(i,j,k)   - cp1(j)*dR(i,j+1,k);
            dZ(i,j,k) = dZ(i,j,k)   - cp1(j)*dZ(i,j+1,k);
            ddR(i,j,k) = ddR(i,j,k) - cp2(j)*ddR(i,j+1,k);
            ddZ(i,j,k) = ddZ(i,j,k) - cp2(j)*ddZ(i,j+1,k);
          enddo

        enddo
      enddo

      call exchange_xxx(dR,dZ,ddR,ddZ)

      cf1 = alp1*cr1_lft(1)
      det1 = 1 - cf1*cf1
      cf2 = alp2*cr2_lft(1)
      det2 = 1 - cf2*cf2

      if (jnode>0) then
        do k=1,nz
          do i=1,nx
            cr_dR(i,k) = ( dR(i,0,k) - cf1* dR(i,1,k) )/det1
            cr_dZ(i,k) = ( dZ(i,0,k) - cf1* dZ(i,1,k) )/det1
            cr_ddR(i,k)= (ddR(i,0,k) - cf2*ddR(i,1,k) )/det2
            cr_ddZ(i,k)= (ddZ(i,0,k) - cf2*ddZ(i,1,k) )/det2
            ! We need corrected gradients at j=0 as well and have all the information
            cor0 = ( dR(i,1,k) - cf1* dR(i,0,k) )/det1
            dR(i,0,k) =  dR(i,0,k) - alp1* cor0*cr1_rgt(ny)
            cor0 = ( dZ(i,1,k) - cf1* dZ(i,0,k) )/det1
            dZ(i,0,k) =  dZ(i,0,k) - alp1* cor0*cr1_rgt(ny)
            cor0 = (ddR(i,1,k) - cf2*ddR(i,0,k) )/det2
            ddR(i,0,k)= ddR(i,0,k) - alp2* cor0*cr2_rgt(ny)
            cor0 = (ddZ(i,1,k) - cf2*ddZ(i,0,k) )/det2
            ddZ(i,0,k)= ddZ(i,0,k) - alp2* cor0*cr2_rgt(ny)
            do j=1,ny
              dR(i,j,k) =  dR(i,j,k) - alp1* cr_dR(i,k)*cr1_lft(j)
              dZ(i,j,k) =  dZ(i,j,k) - alp1* cr_dZ(i,k)*cr1_lft(j)
              ddR(i,j,k)= ddR(i,j,k) - alp2*cr_ddR(i,k)*cr2_lft(j)
              ddZ(i,j,k)= ddZ(i,j,k) - alp2*cr_ddZ(i,k)*cr2_lft(j)
            enddo
          enddo
        enddo
      endif

      if (jnode<npy-1) then
        do k=1,nz
          do i=1,nx
            cr_dR(i,k) = ( dR(i,ny+1,k) - cf1* dR(i,ny,k) )/det1
            cr_dZ(i,k) = ( dZ(i,ny+1,k) - cf1* dZ(i,ny,k) )/det1
            cr_ddR(i,k)= (ddR(i,ny+1,k) - cf2*ddR(i,ny,k) )/det2
            cr_ddZ(i,k)= (ddZ(i,ny+1,k) - cf2*ddZ(i,ny,k) )/det2
            do j=1,ny
              dR(i,j,k) =  dR(i,j,k) - alp1* cr_dR(i,k)*cr1_rgt(j)
              dZ(i,j,k) =  dZ(i,j,k) - alp1* cr_dZ(i,k)*cr1_rgt(j)
              ddR(i,j,k)= ddR(i,j,k) - alp2*cr_ddR(i,k)*cr2_rgt(j)
              ddZ(i,j,k)= ddZ(i,j,k) - alp2*cr_ddZ(i,k)*cr2_rgt(j)
            enddo
          enddo
        enddo
      endif

      end subroutine compact_gradients_y !]
!------------------------------------------
      subroutine thomas_tds(n,alpha,x,b) ![
      ! Solves for the tri-diagonal, symmetric Toeplitz matrix with alpha at the
      ! upper and lower diagonal and ones on the diagonal.
      ! cp can be pre-computed for efficciency

      implicit none

      ! import/export
      integer           :: n     ! size of the problem
      real              :: alpha ! values on off-diagonal
      real,dimension(n) :: x     ! solution
      real,dimension(n) :: b     ! rhs 

      ! local
      integer           :: i
      real,dimension(n) :: cp ! (pre-computed) coefficients

      ! compute cp
      cp(1) = alpha
      do i=2,n
        cp(i) = alpha /(1.0 - alpha * cp(i-1))
      enddo

      ! --- Forward Sweep ---
      x(1) = b(1)
      do i = 2,n
        x(i) = (b(i)/alpha - x(i-1))*cp(i);
      enddo
    
      !--- Backward Substitution ---
      do i = n-1,1,-1
        x(i) = x(i) - cp(i) * x(i+1);
      enddo

      end subroutine thomas_tds  !]

      end module pressure_grad
