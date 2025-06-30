      ! modify to your liking
      ! Non-trivial example in the Filament example

      integer :: i,j,k

      if (FIRST_TIME_STEP) then
        do k=0,nz
          do j=0,ny+1
            do i=0,nx+1

              Akv(i,j,k) = Akv_bak

              Akt(i,j,k,itemp)= Akt_bak(itemp)
#  ifdef SALINITY
              Akt(i,j,k,isalt)= Akt_bak(isalt)
#  endif
            enddo
          enddo
        enddo
      endif

