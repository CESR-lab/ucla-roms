!-------------------------------------------------------------------------
      subroutine set_wec_frc_tile (time,istr,iend,jstr,jend)
      implicit none
      ! import/export
      integer,intent(in) :: istr,iend,jstr,jend
      integer,intent(out):: ierr
      ! local
      integer            :: i,j,k,tmp
      real               :: cff1,cff2
      real               :: tday
     
      if (analytical) then
        call ana_wec_frc(time)
      else 

        tday = time/(24*3600) ! convert from time seconds to days
        if (wec_tim(it2) < tday) then  !! We need to refresh data
C$OMP MASTER 
          if (wec_tim(it1) < 0) then 
            !! First time, we need to read an extra slice of data
            !! Read in the last record with wec_time < tday
            call read_wec_frc(time)
          else ! switch time indices
            tmp = it1
            it1 = it2
            it2 = tmp
          endif
          !! Read in the first record with wec_time > tday
          call read_wec_frc(time)
C$OMP END MASTER 
C$OMP BARRIER
        endif

        ! We now have data with wec_time(it1)<tday and wec_time(it2)<tday
        cf1 = ( wec_tim(it2) - tday )/(wec_tim(it2) - wec_tim(it1) )
        cf2 = ( tday - wec_tim(it1) )/(wec_tim(it2) - wec_tim(it1) )

        ! linear interpolation in time
        awv = cf1*awv(:,:,it1) + cf2*awv(:,:,it2)
        pwv = cf1*pwv(:,:,it1) + cf2*pwv(:,:,it2)
        lwv
        uob
        vob
        sup
        usa
        vsa
        ust
        vst



      endif

      end subroutine set_wec_frc
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
      subroutine read_wec_frc(time)
      ! This version is intended to exclusively read WEC variables that are
      ! part of 'Leonels method'
      implicit none

      character(len=5) :: wec_awv_name = 'Awave'   !! Amplitude
      character(len=5) :: wec_pwv_name = 'Pwave'
      character(len=3) :: wec_lwv_name = 'lmw'
      character(len=4) :: wec_uob_name = 'uorb'
      character(len=4) :: wec_vob_name = 'vorb'
      character(len=3) :: wec_sup_name = 'sup'
      character(len=5) :: wec_usa_name = 'ust2d'    !! Vertical average of u_stokes
      character(len=5) :: wec_vsa_name = 'vst2d'    !! Vertical average of v_stokes
      character(len=4) :: wec_ust_name = 'ust0'     !! Surface u_stokes
      character(len=4) :: wec_vst_name = 'vst0'     !! Surface v_stokes

      character(len=8) :: wec_tim_name = 'wwv_time' !! Wave forcing time

      integer,parameter :: r_var = 0 !! grid type of var (u_var=1,v_var=2,q_var=3)


        ! ifile and irec should be in the module pre-amble, initial as 0
        ! Look for a file that has the requested variable and search
        ! for the last record for which wec_tim < time.

        call find_rec_new(time,time,wec_awv_name,wec_tim_name,ifile,irec)

        ! open the file (maybe check if the file is open already)
        ierr=nf_open (frcfile(ifile), nf_nowrite, ncid)

        ! read the time
        ierr=nf_inq_varid (ncid,wec_tim_name,varid)
        ierr=nf_get_var1_double(ncid, varid, irec, wec_time(it1) )
!       call ncread(frcfile(ifile),wec_tim_name,wec_tim(it1),irec,1)

        ! read the vars
        ierr=nf_inq_varid (ncid,wec_awv_name,varid)
        ierr=ncdf_read(ncid,varid,irec,r_var,awv(START_2D_ARRAY,it1),1)
!       call ncread(frcfile(ifile),wec_tim_name,awv(START_2D_ARRAY,it1),irec,1)
                                      .
                                      .
                                      .
                                      .
                                      .
        ierr=nf_inq_varid (ncid,wec_vst_name,varid)
        ierr=ncdf_read(ncid,varid,irec,r_var,vst(START_2D_ARRAY,it1),1)
!       call ncread(frcfile(ifile),wec_vst_name,vst(START_2D_ARRAY,it1),irec,1)

      end
