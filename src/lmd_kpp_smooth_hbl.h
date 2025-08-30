! Apply horizontal smoothing operator to a private array "wrk" (which
! is actually a CPP-alias for "hbl" or "hbbl").  The array is computed
! within a tile with one row of extra points all-around except lateral
! physical boundaries. The smoothed field is computed over the internal
! range of indices within the tile and is placed back to the same array
! "wrk".  The smoothing is performed while avoiding values under land
! mask which is accomplished by expressing everything via elementary
! differences subject to masking by U- and V-rules (hence applying
! Neumann condition at the coastline).  Without masking the stencil of
! the smoothing operator has the following weights, depending on
! coefficient settings in the code segment below:
!
!   cff=0, cff1=1/8        cff=1/12, cff1=3/16        cff=1/8, cff1=1/4
!
!         1/8                1/32  1/8  1/32           1/16  1/8  1/16
!
!    1/8  1/2  1/8           1/8   3/8  1/8            1/8   1/4  1/8
!
!         1/8                1/32  1/8  1/32           1/16  1/8  1/16
!
!       5-point                 isotropic                 2D 1-2-1
!      Laplacian                Laplacian                 Hanning
!       smoother                smoother                  window 
!
! All three smoothing operators suppress the checkerboard mode in just
! after a single iteration, however, only the last one eliminates flat-
! front 2dx-modes in one iteration; the first and the second attenuate
! the 2dx-mode by factors of 1/2 and 1/4 per iteration.

      cff=1.D0/12.D0 ; cff1=3.D0/16.D0

#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=J_EXT_RANGE
          wrk(0,j)=wrk(1,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=J_EXT_RANGE
          wrk(nx+1,j)=wrk(nx,j)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=I_EXT_RANGE
          wrk(i,0)=wrk(i,1)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=I_EXT_RANGE
          wrk(i,ny+1)=wrk(i,ny)
        enddo
      endif
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        wrk(0,0)=wrk(1,1)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        wrk(0,ny+1)=wrk(1,ny)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        wrk(nx+1,0)=wrk(nx,1)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        wrk(nx+1,ny+1)=wrk(nx,ny)
      endif
#   endif
#  endif

      do j=0,ny+1                 ! The smoothing isotropy is
        do i=1,nx+1                 ! achieved by computing masked
          FX(i,j)=(wrk(i,j)-wrk(i-1,j))  ! elementary differences in
#  ifdef MASKING
     &                      *umask(i,j)  ! each direction first, then
#  endif
        enddo                            ! adding transversal terms
      enddo                              ! expressed via the very same
      do j=1,ny+1                   ! differences.
        do i=0,nx+1
          FE1(i,j)=(wrk(i,j)-wrk(i,j-1))
#  ifdef MASKING
     &                      *vmask(i,j)
#  endif
        enddo
        do i=1,nx
          FE(i,j)=FE1(i,j) + cff*( FX(i+1,j)+FX(i  ,j-1)
     &                            -FX(i  ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=1,ny
        do i=1,nx+1
          FX(i,j)=FX(i,j) + cff*( FE1(i,j+1)+FE1(i-1,j  )
     &                           -FE1(i,j  )-FE1(i-1,j+1))
        enddo
        do i=1,nx
          wrk(i,j)=wrk(i,j) + cff1*( FX(i+1,j)-FX(i,j)
     &                              +FE(i,j+1)-FE(i,j))
#  ifdef MASKING
          wrk(i,j)=wrk(i,j)*rmask(i,j)
#  endif
        enddo
      enddo              !--> discard FX,FE,FE1
