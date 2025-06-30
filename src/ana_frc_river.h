! ONLY FOR ANALYTICAL RIVERS:
! Edit/expand routine as you require.
! Placed here to avoid repo update clash with river_frc.F

      implicit none

      ! Set river volume and tracer values for each time step
      riv_vol(1)   = 5e2  ! Volume flux in m3/s
      riv_trc(1,1) = 24.0 ! Temperature in  Degrees C
#ifdef SALINITY
      riv_trc(1,2) =  1.0  ! Salinity in PSU
#endif
