#include "cppdefs.h"
#if defined SOLVE3D && ( defined TCLIMATOLOGY || ( defined TNUDGING &&\
                         defined T_FRC_BRY  )) && !defined ANA_TCLIMA

                                       ! Read climatology fields for
      subroutine get_tclima(ierr)      ! tracers from climatological

      use param
      use climat
      use ncvars

      implicit none                    ! file at appropriate time.
# define TCLIMA_DATA
# include "scalars.h"
# include "netcdf.inc"
      real(kind=8) dtime
      integer itrc, ierr, lfnm,lvar, lenstr, ncdf_read

! Initialization: Check, whether climatological data file is already
!---------------- opened, and if not, open it.  Find and save netCDF
! IDs for relevant variables.  Determine whether there is cycling to
! reuse the input data and find cycling period "tclm_cycle", set
! initial cycling index "tclm_ncycle" and record index "tclm_rec",
! time index "itclm" and both time record bounds to large negative
! artificial values, so that it will trigger the logic in reading
! part below.

      if (iic == ntstart) then
        do itrc=1,NT
          ittclm(itrc)=1
          tclm_id(itrc)=-1  ;  tclm_time(1,itrc)=-1.E+33
          tclm_rec(itrc)=-1 ;  tclm_time(2,itrc)=-1.E+33
        enddo
      endif

      ierr=nf_noerr
      lfnm=lenstr(clm_file)

      if (iic == ntstart) then
        if (ncidclm == -1) then
          ierr=nf_open(clm_file(1:lfnm), nf_nowrite, ncidclm)
          if (ierr /= nf_noerr) write(*,'(/1x,4A/)')  '### ERROR: ',
     &                         'get_tclima :: can not open file ''',
     &                                    clm_file(1:lfnm),   '''.'
        endif

        itrc=0
        do while(itrc < NT .and. ierr == nf_noerr)
          itrc=itrc+1
          lvar=lenstr(vname(1,indxT+itrc-1))
          ierr=nf_inq_varid(ncidclm, vname(1,indxT+itrc-1)(1:lvar),
     &                                               tclm_id(itrc))
          if (ierr == nf_noerr) then
            if (itrc == itemp) then
              ierr=nf_inq_varid(ncidclm, 'tclm_time', tclm_tid(itrc))
              if (ierr /= nf_noerr) write(*,1) 'tclm_time',
     &                                         itrc, clm_file(1:lfnm)
# ifdef SALINITY
            elseif (itrc == isalt) then
              ierr=nf_inq_varid(ncidclm, 'sclm_time', tclm_tid(itrc))
              if (ierr /= nf_noerr) write(*,1) 'sclm_time',
     &                                         itrc, clm_file(1:lfnm)
# endif
            endif
            if (ierr == nf_noerr) then
              call find_record(      ncidclm, tclm_tid(itrc),
     &                        tclm_rec(itrc), nttclm(itrc),
     &                     tclm_ncycle(itrc), tclm_cycle(itrc),
     &                                     1, 1,         ierr)

              if (ierr == nf_noerr) then
# ifdef MPI_SILENT_MODE
                if (mynode == 0) then
# endif
                  write(*,'(1x,5A,I4)')       'Found ''',
     &                     vname(1,indxT+itrc-1)(1:lvar),
     &         ''' and its time in ''', clm_file(1:lfnm),
     &                        ''' rec =', tclm_rec(itrc)
# ifdef MPI_SILENT_MODE
                endif
# endif
              endif
            endif
# ifdef BIOLOGY
          elseif (itrc >= iNO3) then
            write(*,2) vname(1,indxT+itrc-1)(1:lvar),
     &                        itrc, clm_file(1:lfnm)
# endif 
          else
            write(*,1) vname(1,indxT+itrc-1)(1:lvar),
     &                        itrc, clm_file(1:lfnm)
          endif
        enddo  !<-- while itrc < NT

# ifdef BIOLOGY

! Note that for tracers beyond NO3 it is possible to use analytical
! b.c. input instead of instead of reading data from file, so cancel
! "not-found" error status, if it occurs: the corresponding netCDF
! variable ID retains its -1 value set above which will be used as
! logical condition to do analytical input instead of time
! interpolation.

        if (ierr /= nf_noerr .and. itrc >= iNO3_) ierr=nf_noerr
# endif
      endif !<-- iic == ntstart

  1   format(1x,'### ERROR: get_tclima :: Cannot find variable ''',
     &             A, ''', itrc =', I3,1x, 'in file ''', A, '''.')
# ifdef BIOLOGY
  2   format(1x,  'WARNING: get_tclima :: Cannot find variable ''',
     &          A, ''', itrc =', I3, 1x, 'in file ''', A,'''.'/12x,
     &        'Setting analytical climatology for this variable.')
# endif

! Read data from the file:  Check if model time is bounded by past and
!----- ---- ---- --- -----  future data times: if not, increment record
! and cycling indices, flip time index and read a new portion of data.
! Repeat until model time falls between the two data times.

      do itrc=1,NT
        if (tclm_id(itrc) > 0) then
          do while (tclm_time(ittclm(itrc),itrc) < time+0.5*dt
     &                                  .and. ierr == nf_noerr)

            call advance_record( tclm_rec(itrc), nttclm(itrc),
     &                        tclm_ncycle(itrc), tclm_cycle(itrc),
     &                                        1, 1,  ierr )

            if (ierr == nf_noerr) then
              ierr=nf_get_var1_double(ncidclm, tclm_tid(itrc),
     &                                 tclm_rec(itrc), dtime)
              if (ierr == nf_noerr) then
                ittclm(itrc)=min(3-ittclm(itrc),nttclm(itrc))
                tclm_time(ittclm(itrc),itrc)=dtime*day2sec
     &                 + tclm_cycle(itrc)*tclm_ncycle(itrc)

                ierr=ncdf_read(ncidclm, tclm_id(itrc), tclm_rec(itrc),
     &           r_var, tclima(START_2D_ARRAY,1,ittclm(itrc),itrc), N)

                if (ierr == nf_noerr) then
# ifdef MPI_SILENT_MODE
                  if (mynode == 0) then
# endif
                    write(*,'(3x,A,I2,1x,A,F12.4,2(1x,A,I4))')
     &             'get_tclima :: read climatology for tracer', itrc,
     &             'tclm_time =', dtime,  'rec =', tclm_rec(itrc) MYID
# ifdef MPI_SILENT_MODE
                  endif
# endif
                  if (nttclm(itrc) == 1) goto 3
                else
                  lvar=lenstr(vname(1,indxT+itrc-1))
                  write(*,4) vname(1,indxT+itrc-1)(1:lvar),
     &                                     itrc, tclm_rec(itrc)
                endif
              else
                write(*,4) 'Xclm_time', itrc, tclm_rec(itrc)
              endif
            else
              write(*,'(/1x,A,I3,1x,A,I3/7x,3A/7x,2(A,G12.4)/)')
     &       '### ERROR: get_tclima :: requested time record ',
     &        tclm_rec(itrc), 'exeeds the last record', nttclm(itrc),
     &       'available in netCDF file ''',  clm_file(1:lfnm),  '''',
     &       'tdays = ', tdays,' but the last available tclm_time =',
     &                          tclm_time(ittclm(itrc),itrc)*sec2day
            endif
          enddo  !<-- while()
        endif
  3    continue
      enddo  !<-- itrc
  4   format(1x,'### ERROR: get_tclima :: Cannot read variable ''',
     &                     A, ''', itrc =', I3, ' for record', I3)
      end


      subroutine set_tclima_tile(istr,iend,jstr,jend, ierr)

! Set-up tracer climatology for current tile.

      use param
      use climat
      use ocean3d

      implicit none
      integer istr,iend,jstr,jend, ierr, itrc,i,j,k, it1,it2
      real tmid, cff, cff1, cff2
# define TCLIMA_DATA
# include "scalars.h"

# include "compute_extended_bounds.h"

! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.

      do itrc=1,NT
        if (tclm_id(itrc) > 0) then    !--> available in input file
          if (nttclm(itrc) > 1) then     !--> interpolate in time
            tmid=time+0.5*dt
            it2=ittclm(itrc) ; cff1=tclm_time(it2,itrc)-tmid
            it1=3-it2        ; cff2=tmid-tclm_time(it1,itrc)
            if (cff1 < dt) synchro_flag=.true.
            cff=1./(cff1+cff2)
            cff1=cff1*cff ; cff2=cff2*cff

            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  tclm(i,j,k,itrc)=cff1*tclima(i,j,k,it1,itrc)
     &                            +cff2*tclima(i,j,k,it2,itrc)
                enddo
              enddo
            enddo

            if (cff1 < 0. .or. cff2 < 0.) then
              if (ZEROTH_TILE) then
                write(*,'(/1x,3A/3(1x,A,F16.10)/)')   '### WARNING: ',
     &            'set_tclima_tile :: Current model time is outside ',
     &            'bounds of ''tclm_time''.', 'tclm_tstart =',
     &             tclm_time(it1,itrc)*sec2day, 'tdays =',      tdays,
     &            'tclm_tend =', tclm_time(it2,itrc)*sec2day
              endif
              ierr=ierr+1
            endif
          elseif (iic == ntstart) then        !--> load time-invariant
            do k=1,N                          !    tracer climatology
              do j=jstrR,jendR
                do i=istrR,iendR
                  tclm(i,j,k,itrc)=tclima(i,j,k,ittclm(itrc),itrc)
                enddo
              enddo
            enddo
          endif

# ifdef BIOLOGY
#  define temp cff
#  define SiO4 cff1
        elseif (itrc == iNO3_) then     ! <-- analytical setup for
          do k=1,N                      ! climatology for biological
            do j=jstrR,jendR            ! constituents.
              do i=istrR,iendR
                temp=tclm(i,j,k,itemp)
                if (temp < 8.) then
                  SiO4=30.
                elseif (temp > =8. .and. temp <= 11.) then
                  SiO4=30.-((temp-8.)*(20./3.))
                elseif (temp > 11. .and. temp <= 13.) then
                  SiO4=10.-((temp-11.)*(8./2.))
                elseif (temp > 13. .and. temp <= 16.) then
                  SiO4=2.-((temp-13.)*(2./3.))
                elseif (temp > 16.) then
                  SiO4=0.
                endif

                tclm(i,j,k,iNO3_)=1.67 + SiO4*( 0.5873 + SiO4*(
     &                                0.0144 + SiO4*0.0003099 ))

! Cut off surface NO3

                tclm(i,j,k,iNO3_)=tclm(i,j,k,iNO3_)
     &               *(1-0.5*(tanh((z_r(i,j,k)+100.)/8.)+1.))
              enddo
            enddo
          enddo
        elseif (iic == ntstart .and. itrc == iNH4_) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                if (tclm_tid(iNH4_)<0) tclm(i,j,k,iNH4_)=0.10
                if (tclm_tid(iChla)<0) tclm(i,j,k,iChla)=0.08
                if (tclm_tid(iPhyt)<0) tclm(i,j,k,iPhyt)=0.06
                if (tclm_tid(iZoo_)<0) tclm(i,j,k,iZoo_)=0.04
                if (tclm_tid(iSDet)<0) tclm(i,j,k,iSDet)=0.02
                if (tclm_tid(iLDet)<0) tclm(i,j,k,iLDet)=0.02
              enddo
            enddo
          enddo
#  undef SiO4
#  undef temp
# endif
        endif   ! <-- tclm_tid(itrc)>0
      enddo   ! <-- itrc
      end
#else
      subroutine get_tclima_empty
      end
#endif
