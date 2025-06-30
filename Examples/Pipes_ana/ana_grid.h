
! ----------------------------------------------------------------------
      ! Basic variables that are needed
      real    :: SizeX, SizeY
      real    :: f0,beta
      real    :: x0,y0,dx,dy
      integer :: i,j

      ! Additional vars; depends on your case

      real dh, shelf, slope, land, coast
      real :: depth,max_depth
      real riv_west, riv_east,riv_cells
      real psz,px,py,pipe_cells

      SizeX = 30.0e3  !! Domain size in x-direction [m]
      SizeY = 30.0e3  !! Domain size in y-direction [m]

      f0 = 1.0e-4      !! Coriolis, [1/s]
      beta = 0

      dx = SizeX/gnx   ! Grid spacing
      dy = SizeY/gny

# ifdef MPI
      x0=dx*dble(iSW_corn)             ! Coordinates of south-west
      y0=dy*dble(jSW_corn)             ! corner of MPI subdomain
# else
      x0=0. ; y0=0.
# endif

      do j=-1,ny+2                      ! Extended ranges for x,y arrays
        do i=-1,nx+2                    !
          xr(i,j)=x0+dx*(dble(i)-0.5D0) !
          yr(i,j)=y0+dy*(dble(j)-0.5D0) !

          pm(i,j)=1./dx
          pn(i,j)=1./dy
        enddo
      enddo

! Set Coriolis parameter [1/s] at RHO-points.

      x0=SizeX/2.
      y0=SizeY/2.
      do j=-1,ny+2
        do i=-1,nx+2
          f(i,j)=f0+beta*( yr(i,j)-y0 )
# if defined NONTRAD_COR
!         feta(i,j) = f0*cos(pi/4)
!         fxi(i,j)  = f0*sin(pi/4)
# endif
        enddo
      enddo

      depth =  10
      max_depth = 100;
      shelf=sizeY/5 ! shelf location in meters from south
      slope=(max_depth-depth)/(sizeY*4/5) ! Similar triangles o/a=dh/pm=(max_depth-depth)/(MMm*4/5)
      do j=-1,ny+2
        do i=-1,nx+2

          if(yr(i,j)<shelf) then
            ! Constant shallow region 20% of domain in south.
            h(i,j)=depth
          else
            ! Uniform gradient from south (shallow) to north (deep).
            dh=(yr(i,j)-shelf)*slope
            h(i,j)=depth+dh
          endif

        enddo
      enddo

      ! Set up land masking for river channel
      land  = sizeY*0.1  ! Land extends 10% of domain from south
      coast = sizeY*0.02 ! Coast is not as far
      riv_west=sizeX*0.4 ! River west bank at 40% from west
      riv_east=sizeX*0.6 ! River west bank at 60% from west

      do j=-1,ny+2
        do i=-1,nx+2
          ! default is water
          rmask(i,j) = 1

          if(yr(i,j)<land) then
            if (xr(i,j)<riv_west .or. xr(i,j)>riv_east) then
              rmask(i,j)=0.0
            endif
          endif
          if(yr(i,j)<coast) then !! All land in the far south
            rmask(i,j) = 0.0
          endif
        enddo
      enddo

      if (pipe_source) then
        psz = sizeX*0.02 ! Width of the pipe
        px  = sizeX*.5  ! x location pipe
        py  = sizeY*.5  ! y location pipe
        pipe_cells = nint(psz/dx)**2 !number of cells in this pipe
        do j=-1,ny+2
          do i=-1,nx+2
            pipe_fraction(i,j) = 0.0
            pipe_idx(i,j) = 0
            if (xr(i,j)> px-0.5*psz .and. xr(i,j)<px+0.5*psz) then
              if (yr(i,j)> py-0.5*psz .and. yr(i,j)<py+0.5*psz) then
                 pipe_fraction(i,j) = 1.0/pipe_cells
                 pipe_idx(i,j) = 1
              endif
            endif
          enddo
        enddo
      endif
