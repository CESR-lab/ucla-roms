#include "cppdefs.h"
! DevinD added BULK_FLUX as get sst & sst accounted for in flux modules
! when BULK_FLUX is turned on
#if defined SOLVE3D && (!defined ANA_STFLUX || !defined ANA_SSFLUX) && !defined BULK_FLUX

                                            ! Read surface flux for
      subroutine get_stflux(itrc, ierr_ext) ! tracer itrc from forcing
      implicit none                         ! netCDF file.
      integer itrc, ierr_ext
# define STFLUX_DATA
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "ncvars.h"
# include "netcdf.inc"
      integer ncid, ifrc, ierr, lfrc,lvar,ltvar,lenstr, ncdf_read
      character(len=max_name_size) frcname
      character(len=8) stf_var, stf_time_var
      real(kind=8) dtime
      character(len=64) fmt

# ifdef VERBOSE
      write(*,*) 'enter get_stflux' MYID
# endif
                                      ! Specify netCDF variable name
      if (itrc == itemp) then         ! of the surface forcing flux for
        stf_var='shflux'              ! each tracer and name for the
        stf_time_var='shf_time'       ! corresponding time variable.
# ifdef SALINITY
      elseif (itrc == isalt) then     ! Currently it is done for T and
        stf_var='swflux'              ! S only, but in principle it can
        stf_time_var='swf_time'       ! be expanded to any number of
# endif
      endif                           ! tracers right here without
      lvar=lenstr(stf_var)            ! changing anything throughout
      ltvar=lenstr(stf_time_var)      ! any other parts of the code.


! Initialization:  Check, whether forcing netCDF file is already
!----------------  opened, an if not, open it. Find and save netCDF
! IDs for relevant variables. Determine whether there is cycling to
! reuse the input data and what is cycling period "stf_cycle", find
! initial cycling index "stf_ncycle", and record index "stf_rec".


! Note that these three must be initialized here no matter what: not
! doing it causes segmentation fault and uninitialized variable errors
! in do while(stf_time(itstf(itsmth,...)...) loop  below in the case
! when something goes wrong with reading data, even thought the
! .and.ierr==nf_noerr  logic should prevent execution of this while.
! Set both bounding times to unrealistacally large negative values to
! trigger the while-loop.

      if (iic == ntstart) then
        itstf(itrc)=1    ; stf_time(1,itrc)=-1.E+33
        stf_rec(itrc)=-1 ; stf_time(2,itrc)=-1.E+33
      endif

      ifrc=max(stf_file_id(itrc),1)
 100  ncid=ncfrc(ifrc)
      frcname=frcfile(ifrc)
      lfrc=lenstr(frcname)
      ierr=nf_noerr

      if (stf_rec(itrc) < 0) then  !<-- initialization or file shift
        if (ncid == -1) then
          ierr=nf_open (frcname(1:lfrc), nf_nowrite, ncid)
          if (ierr == nf_noerr) then
            ncfrc(ifrc)=ncid
          else
            write(*,'(/1x,4A,1x,2A,I4 /)') '### ERROR: get_stflux :: ',
     &            'Cannot open netCDF file ''', frcname(1:lfrc), '''.',
     &             nf_strerror(ierr) MYID
          endif
        endif

        if (ierr == nf_noerr) then
          ierr=nf_inq_varid(ncid, stf_var(1:lvar), stf_id(itrc))
          if (ierr == nf_noerr) then
            stf_file_id(itrc)=ifrc
            mpi_nonexit_warn write(*,'(1x,5A,I3)')    'Found ''',
     &              stf_var(1:lvar), ''' in ''', frcname(1:lfrc),
     &                                       ''', ifrc =', ifrc

            ierr=nf_inq_varid (ncid, stf_time_var(1:ltvar),
     &                                       stf_tid(itrc))
            if (ierr == nf_noerr) then
              call find_record(         ncid, stf_tid(itrc),
     &                         stf_rec(itrc), ntstf(itrc),
     &                      stf_ncycle(itrc), stf_cycle(itrc),
     &                                  ifrc, max_frc,  ierr)

              if (stf_rec(itrc) < 0) then
c**             ierr=nf_close(ncid) ; ncfrc(ifrc)=-1
                ifrc=ifrc+1 ;  goto 100  !--> proceed to the next file
              endif

            else
              write(*,2) stf_time_var(1:ltvar), itrc, frcname(1:lfrc)
            endif
          elseif (ifrc < max_frc) then
# ifdef VERBOSE
            mpi_nonexit_warn write(*,1) stf_var(1:lvar), itrc,
     &                                         frcname(1:lfrc)
# endif
c**         ierr=nf_close(ncid) ; ncfrc(ifrc)=-1 ; stf_rec(itrc)=-1
            ifrc=ifrc+1 ;  goto 100   !--> proceed to the next file
          else
            write(*,2)  stf_var(1:lvar), itrc, frcname(1:lfrc)
          endif
        endif
      endif

# ifdef VERBOSE
  1   format(1x,  'WARNING: get_stflux :: Cannot find variable ''',
     &             A, ''', itrc =', I3, 1x, 'in file ''', A, '''.')
# endif
  2   format(1x, '### ERROR: get_stflux :: Cannot find variable ''',
     &             A, ''', itrc =', I3, 1x, 'in file ''', A, '''.')


! Read data from the file:   Check if model time is bounded by past
!----- ---- ---- --- -----   and future data times: if not, increment
! record and cycling indices, flip time index and read a new portion
! of data. Repeat until model time falls between the two data times.

      do while(stf_time(itstf(itrc),itrc) < time+0.5*dt .and.
     &                                      ierr == nf_noerr)
        call advance_record(   stf_rec(itrc), ntstf(itrc),
     &                      stf_ncycle(itrc), stf_cycle(itrc),
     &                                  ifrc, max_frc,  ierr)
     &
        if (stf_rec(itrc) < 0) then
c**       ierr=nf_close(ncid) ; ncfrc(ifrc)=-1
          ifrc=ifrc+1 ; goto 100   !--> proceed to the next file

        elseif (ierr == nf_noerr) then
          ierr=nf_get_var1_double(ncid, stf_tid(itrc), stf_rec(itrc),
     &                                                        dtime)
          if (ierr == nf_noerr) then
            itstf(itrc)=min(3-itstf(itrc),ntstf(itrc))
            stf_time(itstf(itrc),itrc)=dtime*day2sec
     &                     + stf_cycle(itrc)*stf_ncycle(itrc)

            ierr=ncdf_read(ncid, stf_id(itrc),  stf_rec(itrc),  r_var,
     &                     stflxg(START_2D_ARRAY,itstf(itrc),itrc), 1)
            if (ierr == nf_noerr) then
# ifdef MPI_SILENT_MODE
              if (mynode == 0) then
# endif
                write(fmt,'(A,I2,A)') '(10x,3A,I3,T', 60-ltvar,
     &                                 ',2A,F12.4,2(1x,A,I5))'

                write(*,fmt) 'get_stflux :: read ', stf_var(1:lvar),
     &              ' flux for tracer', itrc, stf_time_var(1:ltvar),
     &              ' =', dtime, 'rec =', stf_rec(itrc) MYID
# ifdef MPI_SILENT_MODE
              endif
# endif
              if (ntstf(itrc) == 1) return
            else
              write(*,3) 'stflux', itrc, stf_rec, nf_strerror(ierr)
            endif
          else
            write(*,3) 'stf_time', itrc, stf_rec(itstf(itrc)),
     &                                      nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,2A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')  '### ',
     &   'ERROR: get_stflux :: requested time record ', stf_rec(itrc),
     &   'exceeds the last record',   ntstf(itrc),    'available in ',
     &   'netCDF file ''',    frcname(1:lfrc),    '''',    'tdays = ',
     &    tdays, '  but the last available  stf_time =',
     &                           stf_time(itstf(itrc),itrc)*sec2day
        endif
      enddo
      if (ierr_ext == nf_noerr) ierr_ext=ierr

  3   format(/1x,'### ERROR: get_stflux :: cannot read variable ''',
     &                    A, ''', itrc =', I2, ' rec =', I4, 12x,A/)
# ifdef VERBOSE
      write(*,*) '...leaving get_stflux' MYID
# endif
      end
#else
      subroutine get_stflux_empty
      end
#endif /* SOLVE3D  && !ANA_STFLUX || (SALINITY && !ANA_SSFLUX) */
