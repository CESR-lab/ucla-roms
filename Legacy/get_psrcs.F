#include "cppdefs.h"
#if defined PSOURCE && !defined ANA_PSOURCE

      subroutine get_psrcs (ierr)

! Read in point sources of flow and tracer at the appropriate
! time from psource NetCDF files.  Written by Eileen Deng-Idica
! (05-03-2007) & modified by Yusuke Uchiyama (06-20-2011)

      implicit none
# include "param.h"
# include "forces.h"
# define PSOURCE_DATA
# include "sources.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncvars.h"
# ifdef MPI
#  include "hidden_mpi_vars.h"
# endif
      real cff, csum
      integer ierr, lstr,lvar,lenstr
      integer is, itrc, ncid, varid, i,k
#ifdef MPI
     &                        , i1,i2,j1,j2
#endif
      integer, parameter :: pdims=2
      integer            :: start(pdims), count(pdims)

c      real,    dimension(:,:), allocatable :: D
      character(len=7)  snsrc !, sk
      character(len=64) fmt

c# define CR

! Initialization: Inquire about the contents of psource NetCDF file:
!================  variables and dimensions. Check for consistency.

      if (may_day_flag /= 0) return      !-->  EXIT

      lstr=lenstr(psrc_file)
      ncid=ncidpsrcs
      ierr=nf_noerr

      if (iic == ntstart) then
        itpsrc=1
        psrc_time(1)=-1.E+20
        psrc_time(2)=-1.E+20


! If not opened yet, open psource NetCDF file for reading.
! Find and save IDs for relevant variables.

        if (ncid == -1) then
          ierr=nf_open(psrc_file(1:lstr), nf_nowrite, ncid)
          if (ierr == nf_noerr) then
            ncidpsrcs=ncid
          else
            write(*,'(/1x,4A,1x,2A,I4 /)') '### ERROR: get_psrcs',
     &         ' :: Cannot open netCDF file ''', psrc_file(1:lstr),
     &         '''.',nf_strerror(ierr) MYID
            goto 99                                      !--> ERROR
          endif
        endif
! Nsrc: number of point sources
        ierr=nf_inq_dimid(ncid, 'Nsrc', varid)
        if (ierr == nf_noerr) then
          ierr=nf_inq_dimlen(ncid, varid, Nsrc)
          if (ierr /= nf_noerr) then
            write(*,1) 'Nsrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          else
            if (Nsrc > Msrc) then
              write(*,*) '### ERROR: get_psrcs :: Found Nsrc > Msrc'
              goto 99                                    !--> ERROR
            endif
          endif
        else
          write(*,1) 'Nsrc', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

        write(snsrc,'(i7)') Nsrc
c        if (.not.allocated(D)) allocate (D(Nsrc,N))
! time
        ierr=nf_inq_varid(ncid, 'psrc_time', psrc_tid)
        if (ierr == nf_noerr) then
          call find_record(ncid, psrc_tid,  psrc_rec, ntpsrc,
     &                   psrc_ncycle, psrc_cycle,  1,1, ierr)

        else
          write(*,1) 'psrc_time', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

CR        write(*,*) 'ntpsrc=',ntpsrc
CR        write(*,*) 'psrc_rec=',psrc_rec
CR        write(*,*) 'psrc_cycle=',psrc_cycle
CR        write(*,*) 'psrc_ncycle=',psrc_ncycle
CR        write(*,*) 'psrc_time=',psrc_time(1)

! Isrc & Jsrc (Nsrc): point source locations relative to the global grid

        ierr=nf_inq_varid (ncid, 'Isrc', varid)
        if (ierr == nf_noerr) then
          ierr=nf_get_var_int (ncid, varid, Isrc)
          if (ierr /= nf_noerr) then
            write(*,1) 'Isrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Isrc id', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

        ierr=nf_inq_varid (ncid, 'Jsrc', varid)
        if (ierr == nf_noerr) then
          ierr=nf_get_var_int (ncid, varid, Jsrc)
          if (ierr /= nf_noerr) then
            write(*,1) 'Jsrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Jsrc id', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

# ifdef MPI

! Initially indices "Isrc,Jsrc" are assumed to be defined relative
! to the absolute physical grid. Translate them into would be local
! indices within the current MPI-subdomain and check whether they
! fall within. If not, set unrealistic negative value to ignore
! these particular point sources.

        do is=1,Nsrc
          Isrc(is)=Isrc(is)-iSW_corn   ! Transl
          Jsrc(is)=Jsrc(is)-jSW_corn
          if ( Isrc(is)<iwest  .or. Isrc(is)>ieast  .or.
     &         Jsrc(is)<jsouth .or. Jsrc(is)>jnorth ) then
            Isrc(is)=-100
            Jsrc(is)=-100
          endif
        enddo
# endif

! Dsrc (Nsrc): direction of mass flux (0:xi-, 1:eta-, 2:s-direction)

        ierr=nf_inq_varid (ncid, 'Dsrc', varid)
        if (ierr == nf_noerr) then
          ierr=nf_get_var_int (ncid, varid, Dsrc)
          if (ierr /= nf_noerr) then
            write(*,1) 'Dsrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Dsrc id', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

CR        fmt = '(6x,' / / snsrc / / 'i4)'
CR        write(*,*) 'get_psrcs: Isrc'
CR        write(*,fmt) (Isrc(i),i=1,Nsrc)
CR        write(*,*) 'get_psrcs: Jsrc'
CR        write(*,fmt) (Jsrc(i),i=1,Nsrc)
CR        write(*,*) 'get_psrcs: Dsrc'
CR        write(*,fmt) (Dsrc(i),i=1,Nsrc)

! Qshape (Nsrc,N)

        ierr=nf_inq_varid (ncid, 'Qshape', varid)
        if (ierr == nf_noerr) then
          ierr=nf_get_var_FTYPE (ncid, varid, D)
          if (ierr == nf_noerr) then
            Qshape=0.D0
            do is=1,Nsrc
              cff=0.D0; csum=0.D0
              do k=1,N
                Qshape(is,k)=D(is,k)
                cff=cff+Qshape(is,k)
              enddo
              if (cff/=0.D0) then
                csum=1./cff
              else
                csum=0.D0
              endif
              if (abs(cff-1.0)>1.D-6) then
                write(*,*) 'WARNING : get_psrcs :: Qshape is not ',
     &                     'normalized.  Correction has been made.'
                do k=1,N            ! normalize Qshape
                  Qshape(is,k)=Qshape(is,k)*csum
                enddo
              endif
            enddo
          else
            write(*,1) 'Qshape', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Qshape', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif
CR
CR        write(*,*) 'get_psrcs: Qshape'
CR        write(sk,'(i7)') N
CR        fmt = '(6x,i4,' / / sk / / 'f6.3)'
CR        do i=1,Nsrc
CR          write(*,fmt) i,(Qshape(i,k),k=1,N)
CR        enddo
CR
! Qbar id
        ierr=nf_inq_varid (ncid, 'Qbar', Qbar_id)
        if (ierr /= nf_noerr) then
          write(*,1) 'Qbar', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif
! Tsrc id
        do itrc=1,NT
          lvar=lenstr(vname(1,indxT+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                                            Tsrc_id(itrc))
          if (ierr /= nf_noerr) then
            write(*,1) vname(1,indxT+itrc-1)(1:lvar),
     &            psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        enddo

!  Print diagnosis

# define CD
CD        write(*,*)'[ get_psrcs diagnosis ]' MYID
CD        write(*,'(6x,a,i4)')'number of point sources, Nsrc = ',Nsrc
CD        fmt = '(6x,a,' / / snsrc / / 'i2)'
CD        write(*,fmt)'direction : ', (Dsrc(i),i=1,Nsrc)
CD        fmt = '(6x,3a,' / / snsrc / / 'i2)'
CD        do itrc=1,NT
CD          lvar=lenstr(vname(1,indxT+itrc-1))
CD          write(*,fmt) 'tracer ', vname(1,indxT+itrc-1)(1:lvar)
CD        enddo


      endif             ! iic==ntstart
  1   format(1x,'WARNING: get_psrcs :: Cannot find variable ''',
     &                          A, ''' in file ''', A, '''.')

! Read data from the file: Check if model time is bounded by the past
!----- ---- ---- --- ----- and the future data times: if not, increment
! record and cycling indices, flip time index and read a new portion of
! the data. Repeat until model time falls between the two data times.
CR
CR        write(*,*) 'get_psrcs: ntpsrc=',ntpsrc
CR
      do while (psrc_time(itpsrc) < time+0.5*dt .and. ierr == nf_noerr)
        call advance_record(psrc_rec, ntpsrc,  psrc_ncycle, psrc_cycle,
     &                                                     1,1,  ierr)

        if (ierr == nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, psrc_tid, psrc_rec, cff)
          if (ierr == nf_noerr) then
            itpsrc=3-itpsrc
            psrc_time(itpsrc)=psrc_cycle*psrc_ncycle + cff*day2sec

            start(1)=1        ; count(1)=Nsrc
            start(2)=psrc_rec ; count(2)=1
! Qbar
            ierr=nf_get_vara_FTYPE (ncid, Qbar_id, start, count,
     &                                           QbarG(1,itpsrc))
            if (ierr /= nf_noerr) then
              write(*,2) 'Qbar', psrc_rec, nf_strerror(ierr)
              goto 99                                    !--> ERROR
            endif

CR        write(snsrc,'(i7)') Nsrc
CR        write(*,*) 'get_psrcs: QbarG'
CR        fmt = '(6x,' / / snsrc / / 'f6.2)'
CR        write(*,fmt) (QbarG(i,itpsrc),i=1,Nsrc)

! Tsrc
            do itrc=1,NT
              lvar=lenstr(vname(1,indxT+itrc-1))
              ierr=nf_get_vara_FTYPE (ncid, Tsrc_id(itrc), start,
     &                                count,  TsrcG(1,itpsrc,itrc))
              if (ierr /= nf_noerr) then
                write(*,2) vname(1,indxT+itrc-1)(1:lvar), psrc_rec,
     &                                            nf_strerror(ierr)
                goto 99                                  !--> ERROR
              endif
CR
CR        write(*,*) 'get_psrcs: ',vname(1,indxT+itrc-1)
CR        fmt = '(6x,' / / snsrc / / 'f6.2)'
CR        write(*,fmt) (TsrcG(i,itpsrc,itrc),i=1,Nsrc)
CR
            enddo
              write(*,'(6x,2A,F12.4,2(1x,A,I4))') 'get_psrcs',
     &                '-- read point source file    for psrc_time =',
     &                   cff, 'rec =', psrc_rec MYID
          else
            write(*,2) 'psrc_time', psrc_rec,  nf_strerror(ierr)
          endif
        else
          write(*,'(/2(1x,A,I4)/12x,3A/12x,A,G12.4,A,G12.4/)')
     &          '### ERROR: get_psrcs :: Requested time record ',
     &           psrc_rec,    'exeeds the last record',   ntpsrc,
     &          'available in netCDF file ''',  psrc_file(1:lstr),
     &          '''',  'tdays = ',   tdays,  '   while the last ',
     &         'available psrc_time =', psrc_time(itpsrc)*sec2day
        endif
      enddo


  2   format(/1x,'### ERROR: get_psrcs :: Cannot read variable ''',
     &                                A, ''' for rec =', I4/12x,A/)
      return

  99  may_day_flag=2
      return
      end


      subroutine set_psource_tile(istr,iend,jstr,jend, ierr)

! Set-up psource data for current tile.

      implicit none
      integer istr,iend,jstr,jend, ierr, it1,it2
      integer is, itrc, k
      real cff, cff1,cff2
# include "param.h"
# include "forces.h"
# define PSOURCE_DATA
# include "sources.h"
# include "scalars.h"
# include "ocean3d.h"

# include "compute_extended_bounds.h"

      it1=3-itpsrc ; it2=itpsrc
      cff=time+0.5*dt
      cff1=psrc_time(it2)-cff
      cff2=cff-psrc_time(it1)

!  Load time invariant data.

      if (psrc_cycle < 0.) then
        if (FIRST_TIME_STEP) then
          do is=1,Nsrc
            Qbar(is)=QbarG(is,itpsrc)
            do k=1,N
              Qsrc(is,k)=Qbar(is)*Qshape(is,k)
            enddo
            do itrc=1,NT
              do k=1,N
                Tsrc(is,k,itrc)=TsrcG(is,itpsrc,itrc)
              enddo
            enddo
          enddo
        endif

! Time-interpolate Qbar and Tsrc from point data.  Check that for
! the next time step [when time=time+dt] time+dt is still between
! psrc_time(it1) and psrc_time(it2); if not set set synchro_flag to
! signal that the new forcing data should be  read from the netCDF
! input file (master thread only).

      elseif (cff1 >= 0. .and. cff2 >= 0.) then
        if (cff1 < dt) synchro_flag=.true.
c*        if (ZEROTH_TILE .and. cff1 < dt) synchro_flag=.TRUE.

        cff=1./(cff1+cff2)
        cff1=cff1*cff ; cff2=cff2*cff
        do is=1,Nsrc
          Qbar(is)=cff1*QbarG(is,it1)+cff2*QbarG(is,it2)
c*          if (QbarG(is,it1) == 0. .or. QbarG(is,it2) == 0) then
c*              Qbar(is) = 0.    ! makes sure set_psource does not
c*          endif                ! interpolate more flow than you want
          do k=1,N
            Qsrc(is,k)=Qbar(is)*Qshape(is,k)
          enddo
          do itrc=1,NT
            cff=cff1*TsrcG(is,it1,itrc)+cff2*TsrcG(is,it2,itrc)
            do k=1,N
               Tsrc(is,k,itrc)=cff
            enddo
          enddo
        enddo

!  Unable to set-up Qbar and Tsrc:
!  Complain about the error and signal to quit.

      else
        if (ZEROTH_TILE) then
          write(*,'(/1x,3A/3(1x,A,F16.10)/)')         '### WARNING: ',
     &    'set_psource_tile :: Current model time is outside bounds ',
     &    'of ''psrc_time''.','psrc_tstart =',psrc_time(it1)*sec2day,
     &                        'tdays =',      tdays,
     &                        'psrc_tend =',  psrc_time(it2)*sec2day
          may_day_flag=2
        endif
        ierr=ierr+1
      endif
      end
#else
      subroutine get_psource_empty
      end
#endif /* PSOURCE && !ANA_PSOURCE */
