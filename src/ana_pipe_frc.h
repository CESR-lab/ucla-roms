! ONLY FOR ANALYTICAL PIPES:
! Edit/expand routine as you require.
! Placed here to avoid repo update clash with pipe_frc.F

      implicit none

      ! local
      integer :: i,j,ip

      ! Set river volume and tracer values for each time step
      pipe_vol(1)   = 5e2  ! Volume flux in m3/s
      pipe_trc(1,1) = 24.0 ! Temperature in  Degrees C
      pipe_trc(1,2) =  1.0  ! Salinity in PSU

      pipe_prf(1,:) = 0.0  ! Dispersion profile
      pipe_prf(1,1) = 0.5  ! Dispersion profile
      pipe_prf(1,2) = 0.5  ! Dispersion profile

      pipe_flx = pipe_fraction*pipe_vol(1)
!     do j = 1,ny
!       do i = 1,nx
!         pipe_flx(i,j) = 0.0
!         if (pipe_idx(i,j) > 0) then
!           pidx = pipe_idx(i,j)
!           pipe_flx(i,j) = pipe_frc(i,j)*pipe_vol(pidx)
!         endif
!       enddo
!     enddo

