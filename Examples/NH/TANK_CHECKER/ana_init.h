 
      ! Everything after the implicit none
      ! Replace with something less trivial when needed

      integer :: i,j,k,ig,jg
      real :: x
      real :: drho,thickn,cff
      real :: amp,Lt
      
      drho   = 10
      thickn = 2
      cff=2.*(log((0.99+1.)/(-0.99+1.))/2.)/thickn
      amp = 0.05
      Lt = 10.
      do k=1,nz
        do j=1-bf,ny+bf
          do i=1-bf,nx+bf
       
# ifdef TANKINT
            t(i,j,k,1,itemp)=rho0-0.5*drho*
     &      tanh(cff*(z_r(i,j,k)+h(i,j)/2-amp*cos(pi*xr(i,j)/Lt)))
     &      - 1000.
# else
            t(i,j,k,1,itemp)= 20
# endif
            t(i,j,k,2,itemp)= t(i,j,k,1,itemp)

            t(i,j,k,1,isalt)= 0.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo

      do j=1-bf,ny+bf
        do i=1-bf,nx+bf
          ubar(i,j,1)=0.0
          ubar(i,j,2)=ubar(i,j,1)
          vbar(i,j,1)=0
          vbar(i,j,2)=vbar(i,j,1)
# ifdef TANKINT
          zeta(i,j,1)=0.
# else
          ! 2dx checkerboard of 1 mm from global index parity (x_mid = SizeX/2 = 5 m)
          ig = nint((xr(i,j)+5.0)*pm(i,j)+0.5)
          jg = nint((yr(i,j)+5.0)*pn(i,j)+0.5)
          zeta(i,j,1)=1.e-3*(1-2*mod(ig+jg,2))
# endif
          zeta(i,j,2)=zeta(i,j,1)
          do k=1,nz
            u(i,j,k,1)=0.0
            u(i,j,k,2)=u(i,j,k,1)
            v(i,j,k,1)=0.
            v(i,j,k,2)=v(i,j,k,1)
          enddo
        enddo
      enddo

