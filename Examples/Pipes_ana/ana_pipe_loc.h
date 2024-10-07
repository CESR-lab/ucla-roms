
      real psz, px, py, pipe_cells
      real sizeX, sizeY
      real dx, dy
      integer :: i, j


      SizeX = 30.0e3  !! Domain size in x-direction [m]
      SizeY = 30.0e3  !! Domain size in y-direction [m]
      dx = SizeX/gnx
      dy = SizeY/gny

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
