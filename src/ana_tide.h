
      real    :: f0,D,c,Ld,amp,kh
      complex :: cI
      integer :: i,j

      cI = cmplx(0.0,1.0)
      f0 = 9e-5
      f0 = 1e-11
      D = 4000
      c = sqrt(g*D)   ! = om/k
      Ld = c/f0
      amp = 0.01

      ! Kelvin wave solution
      ! u = uc*exp(-yr/Ld)*exp(cI(kx - om*t))
      ! z = zc*exp(-yr/Ld)*exp(cI(kx - om*t))
      ! just a single tidal component

      ftide(1) = 0.0001405189     ! [s^-1] M2 frequency
      kh  = ftide(1)/c

      if (bry_tides) then
       do j=-1,ny+2          ! Extended ranges for x,y arrays
        do i=-1,nx+2
	  ztide_r(i,j,1) = amp*real(exp(-yr(i,j)/Ld)*exp(cI*kh*xr(i,j)))
	  ztide_i(i,j,1) =-amp*aimag(exp(-yr(i,j)/Ld)*exp(cI*kh*xr(i,j)))
        enddo
       enddo
       do j= 0,ny+1          ! Extended ranges for x,y arrays
        do i= 0,nx+1
	  utide_r(i,j,1) = 0.5*(ztide_r(i,j,1)+ztide_r(i-1,j,1))*g/c
	  utide_i(i,j,1) = 0.5*(ztide_i(i,j,1)+ztide_i(i-1,j,1))*g/c
	  vtide_r(i,j,1) = 0.
	  vtide_i(i,j,1) = 0.
        enddo
       enddo
      endif
      kh = 2*pi/(2e6*1.4)
      if (pot_tides) then
       do j=-1,ny+2          ! Extended ranges for x,y arrays
        do i=-1,nx+2
	  ptide_r(i,j,1) = amp*real(exp(cI*kh*xr(i,j)))
	  ptide_i(i,j,1) =-amp*aimag(exp(cI*kh*xr(i,j)))
        enddo
       enddo
      endif
