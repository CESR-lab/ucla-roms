
      ! Everything after the implicit none

      ! This initial condition assumes linera eof
      ! with Tcoef=2.0e-4

      integer :: istr,iend,jstr,jend, i,j,k, itrc, ierr
      real :: b,h_sbl,H,g_alpha,dbdx,bf_int,dzdx,alpha


      real :: b0          = 5.0e-2
      real :: B_cff       = 0.025    ! gamma in JCM
      real :: lambda_inv  =  8.0     ! m. scale of transition from surface to pycnocline
      real :: Nb          = 1.0e-7   ! surface s^-2
      real :: N0          = 3.0e-5   ! background
      real :: h0          = 60.0     ! m. boundary layer beyond filament
      real :: dh0         = 15.      ! m
      real :: L           = 2000.0   ! m


      ! necessary for FIRST_TIME_STEP flag to work, as per set_global_definitions.h
      ! this had no value for analytical examples previously.
#ifdef EXACT_RESTART
      forw_start=ntstart
#endif

      
      ! Calculate temperature from buoyancy:
      ! Make sure that you are using linear equation of state

      g_alpha = g*Tcoef/rho0
      alpha = Tcoef/rho0


      b=0.0
      do j=-1,ny+2
        do i=-1,nx+2

!         Filament
          h_sbl = h0 + dh0*exp(-(( xr(i,j) /L )**2))

!         Double front
!         h_sbl = h0 + dh0*( tanh( 2e-3*(xr(i,j) + 3.15e3)) - tanh( 2e-3*(xr(i,j) - 3.15e3)))

          do k=1,nz
            b  = b0 + Nb * (z_r(i,j,k) + HD) + 
     &           0.5*N0*( (1+ B_cff) * z_r(i,j,k) - ( 1- B_cff)
     &          *( h_sbl + lambda_inv * 
     &           log(cosh((1./lambda_inv) *(z_r(i,j,k) + h_sbl)))))

            t(i,j,k,1,itemp) = b/(g*alpha)
          enddo
        enddo
      enddo
     
      bf_int = 0
      do k=1,nz
        bf_int = bf_int + Hz(1,1,k)*(
     &           b0 + Nb * (z_r(1,1,k) + HD) + 
     &           0.5*N0*( (1+ B_cff) * z_r(1,1,k) - ( 1- B_cff)
     &          *( h0 + lambda_inv * 
     &           log(cosh((1./lambda_inv) *(z_r(1,1,k) + h0)))))
     &               )/g
      enddo

!     print *,'bf_int: ',bf_int


      do k=1,nz
        do j=-1,ny+2
          do i=-1,nx+2
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
# ifdef SALINITY
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
# endif
          enddo
        enddo
      enddo

      ! Assumes that zeta = 0 for the computation of vertical integrals
      ! This introduces a small error, but short of an iterative approach,
      ! I'm not sure how to fix that

      zeta = 0
      do k=1,nz
        do j=-1,ny+2
          do i=-1,nx+2
	     zeta(i,j,1)= zeta(i,j,1) + t(i,j,k,1,1)*alpha*Hz(i,j,k) 
          enddo
        enddo
      enddo
      zeta(:,:,1) = (zeta(:,:,1) - bf_int)

      ! f*v = g*zx
      do j=0,ny+1
        do i=0,nx+1
	   dzdx =  0.5*pm(i,j)*(zeta(i+1,j,1)-zeta(i-1,j,1) )
	   v(i,j,nz,1) =  g*dzdx/f(i,j)
        enddo
      enddo
      ! Set geostrophic flow in v-direction based on temperature as proxy for buoyancy:
      do k=nz-1,1,0-1
        do j=0,ny+1
          do i=0,nx+1
            dbdx = 0.25*pm(i,j)*g*alpha*(
     &          t(i+1,j,k  ,1,1)-t(i-1,j,k  ,1,1) +
     &          t(i+1,j,k+1,1,1)-t(i-1,j,k+1,1,1) )
            v(i,j,k,1)= v(i,j,k+1,1) -
     &           dbdx*( z_r(i,j,k+1)-z_r(i,j,k) )/f(i,j)
          enddo
        enddo
      enddo

      vbar = 0.
      do k=nz-1,1,0-1
        do j=0,ny+1
          do i=0,nx+1
            vbar(i,j,1)= vbar(i,j,1) + v(i,j,k,1)*Hz(i,j,k)/HD
          enddo
        enddo
      enddo

!     v(:,:,1,1)=0.0
!     do k=2,nz
!       do j=0,ny+1
!         do i=0,nx+1
!           dbdx = 0.25*pm(i,j)*g_alpha*(
!    &          t(i+1,j,k  ,1,1)-t(i-1,j,k  ,1,1) +
!    &          t(i+1,j,k-1,1,1)-t(i-1,j,k-1,1,1) )
!           v(i,j,k,1)= v(i,j,k-1,1) +
!    &           dbdx*( z_r(i,j,k)-z_r(i,j,k-1) )/f(i,j)
!         enddo
!       enddo
!     enddo




      do j=1,ny
        do i=1,nx
          ubar(i,j,1)=0.
          vbar(i,j,2)=vbar(i,j,1)
          zeta(i,j,2)=zeta(i,j,1)
          do k=1,nz
            u(i,j,k,1)=0.
            u(i,j,k,2)=u(i,j,k,1)
            v(i,j,k,2)=v(i,j,k,1)
          enddo
        enddo
      enddo


