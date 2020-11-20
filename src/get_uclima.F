#include "cppdefs.h"
#if defined UCLIMATOLOGY && !defined ANA_UCLIMA

                                      ! Read momentum climatology
      subroutine get_uclima(ierr)     ! fields from netCDF file
      implicit none
      integer ierr i, lstr,lvar,lenstr, ncdf_read
      real(kind=8) dtime
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "netcdf.inc"
# include "climat.h"

! Initialization:  Check, whether climatological data file is already
!----------------  opened, and if not, open it.  Find and save netCDF
! IDs for relevant variables.  Determine whether there is cycling to
! reuse the input data and find cycling period "uclm_cycle", set
! initial cycling index "uclm_ncycle" and record index "uclm_rec",
! time index "iuclm" and both time record bounds to large negative
! artificial values, so that it will trigger the logic in reading
! part below.

      ierr=nf_noerr
      lstr=lenstr(clm_file)
      if (iic==ntstart) then
        ituclm=1
        uclm_time(1)=-1.E+20
        uclm_time(2)=-1.E+20

        if (ncidclm==-1) then
          ierr=nf_open (clm_file(1:lstr), nf_nowrite, ncidclm)
          if (ierr /= nf_noerr) write(*,'(/1x,4A/)') '### ERROR: ',
     &                    'get_uclima :: Cannot open netCDF file ''',
     &                                     clm_file(1:lstr),   '''.'
        endif
        if (ierr==nf_noerr) then
          ierr=nf_inq_varid (ncidclm, 'uclm_time', uclm_tid)
          if (ierr==nf_noerr) then
            lvar=lenstr(vname(1,indxUb))
            ierr=nf_inq_varid (ncidclm, vname(1,indxUb)(1:lvar),
     &                                                 ubclm_id)
            if (ierr==nf_noerr) then
              lvar=lenstr(vname(1,indxVb))
              ierr=nf_inq_varid (ncidclm, vname(1,indxVb)(1:lvar),
     &                                                  vbclm_id)
              if (ierr==nf_noerr) then
# ifdef SOLVE3D
                lvar=lenstr(vname(1,indxU))
                ierr=nf_inq_varid (ncidclm, vname(1,indxU)(1:lvar),
     &                                                     uclm_id)
                if (ierr==nf_noerr) then
                  lvar=lenstr(vname(1,indxV))
                  ierr=nf_inq_varid (ncidclm, vname(1,indxV)(1:lvar),
     &                                                      vclm_id)
                  if (ierr==nf_noerr) then
# endif
                    call find_record(ncidclm,uclm_tid, uclm_rec,ntuclm,
     &                             uclm_ncycle,uclm_cycle,  1,1, ierr)
# ifdef SOLVE3D
                  else
                    write(*,1) vname(1,indxV)(1:lvar),
     &                                     clm_file(1:lstr)
                  endif
                else
                  write(*,1) vname(1,indxU)(1:lvar),
     &                                     clm_file(1:lstr)
                endif
# endif
              else
                write(*,1) vname(1,indxVb)(1:lvar),
     &                                    clm_file(1:lstr)
              endif
            else
              write(*,1) vname(1,indxUb)(1:lvar),
     &                                    clm_file(1:lstr)
            endif
          else
            write(*,1) 'uclm_time', clm_file(1:lstr)
          endif
        endif
      endif
  1   format(1x,'### ERROR: get_uclima :: Cannot find variable ''',
     &                               A, ''' in file ''', A, '''.')

! Read data from the file: Check if model time is bounded by the past
!----- ---- ---- --- ----- and the future data times: if not, increment
! record and cycling indices, flip time index and read a new portion of
! the data. Repeat until model time falls between the two data times.

      do while(uclm_time(ituclm) < time+0.5*dt .and. ierr == nf_noerr)
        call advance_record(uclm_rec,ntuclm, uclm_ncycle,uclm_cycle,
     &                                                    1,1, ierr)

        if (ierr == nf_noerr) then
          ierr=nf_get_var1_double(ncidclm, uclm_tid, uclm_rec, cff)
          if (ierr == nf_noerr) then
            ituclm=min(3-ituclm,ntuclm)
            uclm_time(ituclm)=cff*day2sec + uclm_cycle*uclm_ncycle
            ierr=ncdf_read(ncidclm, ubclm_id, uclm_rec, u_var,
     &                     ubclima(START_2D_ARRAY,ituclm), 1)
            if (ierr == nf_noerr) then
              ierr=ncdf_read(ncidclm, vbclm_id, uclm_rec, v_var,
     &                       vbclima(START_2D_ARRAY,ituclm), 1)
              if (ierr == nf_noerr) then
# ifdef SOLVE3D
                ierr=ncdf_read(ncidclm, uclm_id,  uclm_rec,  u_var,
     &                         uclima(START_2D_ARRAY,1,ituclm), N)
                if (ierr == nf_noerr) then
                  ierr=ncdf_read(ncidclm, vclm_id,  uclm_rec,  v_var,
     &                           vclima(START_2D_ARRAY,1,ituclm), N)
                  if (ierr == nf_noerr) then
# endif
#ifdef MPI_SILENT_MODE
                    if (mynode == 0) then
#endif
                      write(*,'(4x,A,5x,A,G12.4,2(1x,A,I4))')
     &               'get_uclima :: read momentum climatology',
     &               'uclm_time =',cff, 'rec =',uclm_rec MYID
# ifdef MPI_SILENT_MODE
                    endif
# endif
                    if (ntuclm == 1) return
# ifdef SOLVE3D
                  else
                    lvar=lenstr(vname(1,indxV))
                    write(*,2) vname(1,indxV)(1:lvar), uclm_rec
                  endif
                else
                  lvar=lenstr(vname(1,indxU))
                  write(*,2) vname(1,indxU)(1:lvar), uclm_rec
                endif
# endif
              else
                lvar=lenstr(vname(1,indxVb))
                write(*,2) vname(1,indxVb)(1:lvar), uclm_rec
              endif
            else
              lvar=lenstr(vname(1,indxUb))
              write(*,2) vname(1,indxUb)(1:lvar), uclm_rec
            endif
          else
            write(*,2) 'clm_time', uclm_rec
          endif
        else
           write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &    '### ERROR: get_uclima :: requested time record', uclm_rec,
     &    'exeeds the last record',  ntuclm,  'available in netCDF ',
     &    'file ''',    clm_file(1:lstr),  '''',   'tdays = ', tdays,
     &    '  but the last available  uclm_time =',
     &                                   uclm_time(ituclm)*sec2day
        endif
      enddo
  2   format(/1x,'ERROR in get_uclima: cannot read variable ''',
     &                                 A, ''' for record ', I4/)
      return
      end

      subroutine set_uclima_tile (istr,iend,jstr,jend, ierr)

!  Set-up momentum climatology for current tile.

      implicit none
      integer istr,iend,jstr,jend, ierr, i,j,k, it1,it2
      real cff, cff1, cff2
# include "param.h"
# include "climat.h"
# include "scalars.h"

# include "compute_extended_bounds.h"

! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.

      if (ntuclm > 1) then                  !--> interpolate velocity
        tmid=time+0.5*dt                       ! climatolog in time
        it2=ituclm ; cff1=uclm_time(it2)-tmid
        it1=3-it2  ; cff2=tmid-uclm_time(it1)
        if (cff1 < dt) synchro_flag=.true.
        cff=1./(cff1+cff2)
        cff1=cff1*cff ; cff2=cff2*cff

        do j=jstrR,jendR
          do i=istrR,iendR
            ubclm(i,j)=cff1*ubclima(i,j,it1)+cff2*ubclima(i,j,it2)
            vbclm(i,j)=cff1*vbclima(i,j,it1)+cff2*vbclima(i,j,it2)
          enddo
        enddo
# ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
            do i=istrR,iendR
              uclm(i,j,k)=cff1*uclima(i,j,k,it1)+cff2*uclima(i,j,k,it2)
              vclm(i,j,k)=cff1*vclima(i,j,k,it1)+cff2*vclima(i,j,k,it2)
            enddo
          enddo
        enddo
# endif
        if (cff1 < 0. .or. cff2 < 0.) then
          if (ZEROTH_TILE) then
            write(*,'(/1x,3A/3(1x,A,F16.10)/)')      '### WARNING: ',
     &   'set_uclima_tile :: Current model time is out of bounds of',
     &   '''uclm_time''.',  'uclm_tstart =',  uclm_time(it1)*sec2day,
     &   'tdays =', tdays,  'uclm_tend =',    uclm_time(it2)*sec2day
          endif
          ierr=ierr+1
        endif
      elseif (iic == ntstart ) then         ! load time-invariant 
        do j=jstrR,jendR                    ! momentum climatology
          do i=istrR,iendR
            ubclm(i,j)=ubclima(i,j,ituclm)
            vbclm(i,j)=vbclima(i,j,ituclm)
          enddo
        enddo
# ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
            do i=istrR,iendR
              uclm(i,j,k)=uclima(i,j,k,ituclm)
              vclm(i,j,k)=vclima(i,j,k,ituclm)
            enddo
          enddo
        enddo
# endif
      endif
      end
#else
      subroutine get_uclima_empty
      end
#endif /* UCLIMATOLOGY && !ANA_UCLIMA */
