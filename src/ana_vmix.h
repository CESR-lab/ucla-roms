      ! modify to your liking
      ! Non-trivial example in the Filament example

      integer :: i,j,k

      if (FIRST_TIME_STEP) then
        do k=0,nz
          do j=0,ny+1
            do i=0,nx+1
              
              Akv(i,j,k) = 1e-5

              Akt(i,j,k,itemp)= 1e-5
#  ifdef SALINITY
              Akt(i,j,k,isalt)= 1e-5
#  endif
            enddo
          enddo
        enddo
      endif

