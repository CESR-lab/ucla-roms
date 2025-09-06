      ! modify to your liking
      ! Non-trivial example in the Filament example

      integer :: i,j,k

      if (FIRST_TIME_STEP) then
        do k=0,nz
          do j=0,ny+1
            do i=0,nx+1
              Akv(i,j,k) = Akv_bak
              Akt(i,j,k) = Akt_bak(1)
            enddo
          enddo
        enddo
      endif

