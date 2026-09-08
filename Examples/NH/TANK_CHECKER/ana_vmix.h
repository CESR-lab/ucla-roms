      ! modify to your liking
      ! Non-trivial example in the Filament example

      integer :: i,j,k

      if (FIRST_TIME_STEP) then
        do k=0,nz
          do j=0,ny+1
            do i=0,nx+1
              
!             Akv(i,j,k) = Akv_bak + 1e-3*exp(-(z_w(i,j,k)-z_w(i,j,1) )/5)
              Akv(i,j,k) = Akv_bak

!             Akt(i,j,k,itemp)= Akt_bak(itemp)
!             Akt(i,j,k,itemp)= Akt_bak(itemp)*exp( z_w(i,j,k)/10)
              Akt(i,j,k)= Akt_bak*exp( z_w(i,j,k)/10)
#  ifdef SALINITY
!             Akt(i,j,k)= Akt_bak(isalt)
#  endif
            enddo
          enddo
        enddo
      endif

