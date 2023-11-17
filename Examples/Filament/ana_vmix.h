      integer :: i,j,k
      real    :: cff
      real    :: hbl,ssgm

      hbl = 70;  ! match the initial mixed layer depth

      if (FIRST_TIME_STEP) then
        do k=0,nz
          do j=0,ny+1
            do i=0,nx+1
              ssgm=(z_w(i,j,nz)-z_w(i,j,k))/hbl

              if (ssgm < 1.) then
                if (ssgm<0.07D0) then
                  cff=0.5*(ssgm-0.07D0)**2/0.07D0
                else
                  cff=0.D0
                endif
                cff=cff + ssgm*(1.-ssgm)**2
              else
                cff=0
              endif

              ! attempt to get a max value of 0.005 m^2/s
              Akv(i,j,k) = cff*0.005/0.15

              Akt(i,j,k,itemp)=cff*0.005/0.15
#  ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
            enddo
          enddo
        enddo
      endif

