#include "cppdefs.h"
#ifdef FLOATS
                                       ! Read initial conditions
      subroutine get_init_flts         ! for the floats variables
                                       ! from netCDF file.
      implicit none
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"
#include "floats.h"
      real bufftrack(1:6,NFT+1,Mfloats)
      integer i, j, ncid, indx, varid,  ierr, lstr, lvar, lenstr,
     &    start(2), count(2),  ibuff(Mfloats),  checkdims, temp(4),
     &    startTinfo(3),countTinfo(3), starttrack(4),counttrack(4)


! Open initial conditions netCDF file for floats reading.
! set record from which the data will be read.
!
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of
! value of nrrec supplied from the parameter file; (2) if the
! file has multiple records and nrrec is positive, then nrrec is
! used, provided that nrrec is within the available records;
! (3) if the file has multiple records and nrrec<0, then THE LAST
! available record is used.

      if (may_day_flag/=0) return      !-->  EXIT
      lstr=lenstr(ininame)
      ierr=nf_open(ininame(1:lstr), nf_nowrite, ncid)
      if (ierr == nf_noerr) write (*,*) 'DAUBE'
        ierr=checkdims (ncid, ininame, indx)
        if (ierr /= nf_noerr) write (*,*) 'DAUBE2'
          if (indx==0) then
            indx=1
          elseif (indx>0 .and. nrrec>0 .and. nrrec<=indx) then
            indx=nrrec
          endif

      else
        write(*,'(/1x,4A/12x,A)')     '### ERROR: get_init_flts :: ',
                'Cannot open netCDF file ''', ininame(1:lstr) ,'''.',
                 nf_strerror(ierr)
       goto 99                                           !--> ERROR
      endif

! ---- -- -------- ----- ----------
! nfloats: find netCDF id, read value

      ierr=nf_inq_varid (ncid, 'nfloats', varid)
      ierr=nf_get_var1_int (ncid, varid, indx, nfloats)
      if (ierr /= nf_noerr) then
        write (*,'(/1x,2A/8x,2A/)') '### ERROR: get_init_flts ',
     &        ':: Cannot read nfloats in netCDF file ''',
     &             ininame(1:lstr),'''.'
              goto 99                                    !--> ERROR
      endif

!------------------------------------
! Tinfo: find netCDF id and read values

      ierr=nf_inq_varid (ncid, 'Tinfo', varid)
      if (ierr == nf_noerr) then
          startTinfo(1)=1
          startTinfo(2)=1
          startTinfo(3)=indx
          countTinfo(1)=5
          countTinfo(2)=nfloats
          countTinfo(3)=1
          ierr=nf_get_vara_FTYPE (ncid, varid, startTinfo,
     &                          countTinfo, Tinfo)
          if (ierr /= nf_noerr) then
            write (*,'(/1x,4A/)') '### ERROR: get_init_flts :: ',
     &       'Cannot read Tinfo from file ''', ininame(1:lstr), '''.'
              goto 99                                    !--> ERROR
          endif
      endif

!------------------------------------
! grid level :  find netCDF id and read values

      ierr=nf_inq_varid (ncid, 'grid_level', varid)
      if (ierr == nf_noerr) then
          start(1)=1
          start(2)=indx
          count(1)=nfloats
          count(2)=1
          ierr=nf_get_vara_int (ncid, varid, start, count,fltgrd)
          if (ierr /= nf_noerr) then
            write (*,'(/1x,4A/)') '### ERROR: get_init_flts :: ',
     &    'Cannot read grid_level from file ''',ininame(1:lstr),'''.'
              goto 99                                    !--> ERROR
          endif
      endif

!------------------------------------
! track: find netCDF id and read values

      ierr=nf_inq_varid (ncid, 'track', varid)
      if (ierr == nf_noerr) then
          starttrack(1)=1
          starttrack(2)=1
          starttrack(3)=1
          starttrack(4)=indx
          counttrack(1)=6
          counttrack(2)=NFT+1
          counttrack(3)=nfloats
          counttrack(4)=1

          ierr=nf_get_vara_FTYPE (ncid, varid, starttrack,
     &                          counttrack, bufftrack)

          if (ierr == nf_noerr) then

            temp(1)=nfp1  ! WARNING: This takes into
            temp(2)=nf    ! account that these time step
            temp(3)=nfm1  ! indices will be rotated once
            temp(4)=nfm2  ! before the first computation.

              do j=1,nfloats
                do i=1,NFT+1
                  track(ixgrd,temp(i),j)=bufftrack(1,i,j)
                  track(iygrd,temp(i),j)=bufftrack(2,i,j)
                  track(izgrd,temp(i),j)=bufftrack(3,i,j)
                  track(ixrhs,temp(i),j)=bufftrack(4,i,j)
                  track(iyrhs,temp(i),j)=bufftrack(5,i,j)
                  track(izrhs,temp(i),j)=bufftrack(6,i,j)
                enddo
              enddo
          else
            write (*,'(/1x,4A/)') '### ERROR: get_init_flts :: ',
     &         'Cannot read track from file ''',ininame(1:lstr),'''.'
              goto 99                                    !--> ERROR
          endif
      endif

!  Close input NetCDF file.

      ierr=nf_close(ncid)
#ifdef MPI_SILENT_MODE
      if (mynode==0) then
#endif
        write(*,'(1x,4A)')  'get_init_flts :: initialized floats ',
     &               'from netCDF file ''', ininame(1:lstr), '''.'
#ifdef MPI_SILENT_MODE
      endif
#endif
      return
  99  may_day_flag=2
      return
      end
#else
      subroutine get_init_flts_empty
      end
#endif  /* FLOATS */
