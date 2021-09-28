#include "cppdefs.h"
c--#define VERBOSE

      subroutine find_record(ncid, varid,  irec, ntime,  icycle,
     &                     cycle_length, ifile, max_files, ierr)

! Find time record number "irec" and cycling index "icycle" which
! enclose the current model time as
!
!            data_time(irec)  <=  time  <  data_time(irec+1)
!
! where, in the case of cycling data_time my be incremented by
!
!            data_time = data_time + icycle*cycle_length
!
! Having/not having  cycling of input data which is identified by
! having positive value of "cycle_length" stored as an attribute to
! timing variable known to this routine by its netCDF variable ID
! "varid".

! Arguments:  Note that arguments are logically arranged in pairs
!------------                     like ..., index,its_max_value, ...
! input:  ncid   netCDF file ID
!        varid   netCDF variable ID for time of forcing field
!        ifile   file index in sequence of files
!    max_files   total number of files in the sequence.

! output: irec   record number stepped back by one (two in the case
!                   of initialization, so the calling get_smth routine
!                   increments it forward by 1 and reads the data at
!                   proper time (does it twice during initialization),
!                   or irec=-1 to signal calling routine that all
!                   records in current file are in the past relatively
!                   to model time, so proceed examining next file;
!        ntime   total number of records in the file (may be
!                    restricted by periodicity range);
!       icycle   cycling index
! cycle_length   duration of time cycling period [seconds]
!         ierr   error status

! Error signalling:  Three outcomes are possible upon execution:
!------ -----------  (i) success: an appropriate record found,
!                                        irec > 0, ierr = 0 ;
!                    (ii) signal to renew search: record not found
!                         within the present file, by list of files
!                         is not exhausted yet:   irec=-1, ierr=0
!                    (iii) error: either a netCDF error within this
!                          routine (return ierr=netCDF status), or
!                          time bounding condition is not achieved
!                          and there are no more files to try next:
!                          ierr=-999 (irec=1, but it does not matter)
      use param
      use scalars
      use netcdf
      implicit none
      integer ncid, varid, irec, ntime, icycle, ifile, max_files, ierr
      real(kind=8) cycle_length

#include "netcdf.inc"

      integer nvdims, nvatts, vartype, vdims(8),
     &        size, i, it1,it2, lvar,ldim,lenstr
      real(kind=8) tstart, tend, test, t1,t2
      character(len=16) varname, dimname

      ntime=1         ! Initialize: if not found later otherwise, it
      cycle_length=0. ! is assumed that the file contains only a single
      icycle=0        ! record for the requested variable, and no
      irec=1          ! temporal interpolation needed/can be performed.

! Make general inquiry into for variable with netCDF ID "varid" and
! find whether it has only is a single dimension.   If there is more
! than one record available check whether the variable has cycling
! attribute, and if so, retrieve cycling period (leaving it set to
! zero means that no recycling will be performed).

      ierr=nf_inq_var(ncid, varid, varname, vartype,
     &                      nvdims,  vdims,  nvatts)
      if (ierr == nf_noerr) then
        lvar=lenstr(varname)
        if (nvdims == 1) then
          ierr=nf_inq_dim(ncid, vdims, dimname, size)
          if (ierr == nf_noerr) then
            ldim=lenstr(dimname)
            ntime=size
            if (ntime > 1) then
              ierr=nf_get_att_double(ncid, varid, 'cycle_length',
     &                                             cycle_length)

                  !  No need for error checking here because
                  ! "cycle_length" may be legitimately absent.

! Read the last record of time variable perform one of the two checks:

              ierr=nf_get_var1_double(ncid, varid, ntime, tend)
              if (ierr == nf_noerr) then

! either, (1) in the case of cycling, check whether all time records
! lie within the declared cycling interval, i.e.,

!     time(rec=1)  <  time(rec)  <=  time(rec=1)+cycle_length
!          for all   1  <  rec  <=  ntime

! and if not restrict the last record to be read "ntime" to those which
! are within (this means that all the remaining records in the file
! will be ignored)...


                if (cycle_length > 0.) then
                  i=ntime
                  do while(tend > cycle_length .and. i > 1 .and.
     &                                          ierr == nf_noerr)
                    i=i-1
                    ierr=nf_get_var1_double(ncid, varid, i, tend)
                    if (ierr /= nf_noerr) write(*,1)
     &                      varname(1:lvar), i, nf_strerror(ierr)
                  enddo

                  if (i < ntime) then
#ifdef MPI_SILENT_MODE
                    if (mynode == 0) then
#endif
                      write(*,'(8x,A,2(A,I4,1x),A)')       'WARNING: ',
     &               'find_record :: Adjusted number of records from ',
     &                ntime, 'to', i,  'restricted by "cycle_length".'
#ifdef MPI_SILENT_MODE
                    endif
#endif
                    ntime=i
                  endif

! ...or, (2), conversely, in the case of non-cycling, check whether
! the last record is in future relative to the current model time, and
! if not, which is equivalent to say that the file does not have usable
! records, get out, while signalling to proceed to the next file.

                elseif (time > tend*day2sec) then
                  if (ifile < max_files) then
                    irec=-1 ;  return   !--> proceed to next file
                  else
                    write(*,'(/1x,2A/)') '### ERROR: find_record :: ',
     &             'run out of time records in non-recycling regime.'
                    ierr=-9999
                  endif
                endif
              else
                write(*,1) varname(1:lvar), ntime, nf_strerror(ierr)
              endif
            endif  !<-- ntime > 1
#ifdef VERBOSE
# ifdef MPI_SILENT_MODE
            if (mynode == 0) then
# endif
              write(*,'(1x,3A,I5,2x,A,ES16.9)')   'find_record: ',
     &              'varname = ''', varname(1:lvar), ''' ntime =',
     &                     ntime, 'cycle_length =', cycle_length
# ifdef MPI_SILENT_MODE
            endif
# endif
#endif
          else
            write(*,'(/1x,4A/12x,A/)')    '### ERROR: find_record :: ',
     &           'Cannot make inquiry about dimension for variable ''',
     &                      varname(1:lvar), '''.', nf_strerror(ierr)
          endif
        else
          write(*,'(1x,4A,I3,A)')    '### ERROR: find_record :: ',
     &                 'variable ''',  varname(1:lvar), ''' has ',
     &                 'wrong number of dimensions,', nvdims, '.'
          ierr=-9999
        endif
      else
        write(*,'(/1x,2A,I4,1x,A/12x,A/)')   '### ERROR: find_record ',
     &               ':: Cannot make general inquiry about variable #',
     &                     varid, 'in netCDF file.', nf_strerror(ierr)
      endif
      if (ntime == 1 .or. ierr /= nf_noerr) return


!> Search for the minimum possible record number which has its time
!>======= in the future relatively to model time: this boils down to
! 2x2 scenarios [cycling vs. not] x [bounded on both sides vs. on
! future side only], which in the code below reduced to just three
! logical branches by merging [cycling vs. not] in both-side case.


! Both sides vs. future only:  In non-cycling case future-only bounding
!----- ----- -- ------- -----  is possible only during roll-over from
! one file to the next in time sequence. In this case the "past" record
! comes from the previous file, while the "future" is the first record
! of the new. In this case only one record needs to be read, and it is
! OK to not have the past bound.


! Cycling: First look at the first record and, in the case of cycling,
!--------- compare its value "tstart" with current model time to set
! initial cycling index "icycle" in such a way, that current model
! time is guaranteed to be within the bounds:

!  tstart+icycle*cycle_length <= time < tstart+(icycle+1)*cycle_length

! because "tend" is restricted to be within one "cycle_length" from
! "tstart",

!        tstart+(icycle+1)*cycle_length > tend + icycle*cycle_length

! this leaves two possibilities:
! either
!        tstart+icycle*cycle_length <= time < tend+icycle*cycle_length
! or
!     tend+icycle*cycle_length <= time < tstart+(icycle+1)*cycle_length

! so the first case model "time" is bounded by two consecutive records
! in the file, while in the second in can be bounded by periodicity
! cycle increment.

! Merged:    If model time is bounded by two consecutive time records,
!--------    it does not matter whether there is cycling or not, since
! the same number of whole cycles is added to both "tstart" and "tend".
! In this case, to accelerate the search, use division in halves
! instead linear scanning of all the records in the file. Here it
! should be kept in mind that intervals between time records may be
! uneven, and some records may be even repeated, so a more efficient
! secant-type algorithm may actually fail.


! Because no actual data has been read thus far, back step "irec" so
!-------- -- ------ ---- --- ---- ---- the caller (e.g., get_smth) will
! advance it forward to what it should be (recall that the first thing
! get_smth does is incrementing record index).  NOTE: it is OK when
! irec=0 as outcome of the back step, because "advance_record" call
! will make it 1 again.

! Limitations:  Two conditions must be met for this algorithm to be
!-------------  applicable:
!               (1) In the case of periodic cycling the entire time
!                   series data must be within a SINGLE file;
!               (2) In the case of using multiple files in time series
!                   model start (including restart) is possible only if
!                   BOTH records bounding time of the initial condition
!                   belong to THE SAME file.

                                        ! convert "tend" (already read
      tend=tend*day2sec                 ! above) and cycle_length from
      cycle_length=cycle_length*day2sec ! days to seconds.

      ierr=nf_get_var1_double(ncid, varid, irec, tstart)
      if (ierr == nf_noerr) then
        tstart=tstart*day2sec
        if (cycle_length > 0.) then
          icycle=int((time-tstart)/cycle_length)
          if (time < tstart) icycle=icycle-1
          tstart=tstart + icycle*cycle_length
          tend  =tend   + icycle*cycle_length
        endif

        if (time<tstart .and. iic>ntstart) then   !--> file roll-over
          irec=0 !<-- this takes into account that advance_record will
                 !  turn it into irec=1 before reading the actual data

        elseif (tstart <= time .and. time < tend) then  !--> 2-sided
          it1=1 ; t1=tstart ; it2=ntime ; t2=tend       !    bounding
          do while(it1+1 < it2 .and. ierr == nf_noerr)
            i=(it1+it2)/2                                 ! Search
            ierr=nf_get_var1_double(ncid, varid, i, test) ! algorithm
            if (ierr == nf_noerr) then                    ! by division
              test=test*day2sec+icycle*cycle_length       ! in halves.
              if (time < test) then
                it2=i ; t2=test
              else
                it1=i ; t1=test
              endif
            else
              write(*,1) varname(1:lvar), i, nf_strerror(ierr)
            endif
          enddo
          irec=it2  !<-- initially set to future record, then step back

          if (irec > 1) then
            irec=irec-1                  !--> back step record normally
          elseif (cycle_length > 0.) then
            irec=ntime ; icycle=icycle-1 !--> rewind time within file
c**       elseif (ifile > 1) then
c**         irec=-2                      !--> proceed to previous file
          else
            write(*,'(/1x,2A/)')  '### ERROR: find_record :: run out ',
     &                      'of time records in non-recycling regime.'
            ierr=-9999
          endif

        elseif (tend <= time .and.cycle_length > 0.) then !--> bounded
          irec=ntime          ! from two sides via periodicity; Records
                              ! ntime and 1 (with incremented "icycle")
                              ! will be read  during initialization.
        else
          write(*,*) '### ERROR: find_record :: Algorithm limitation.'
          ierr=-9999
        endif

! Furthermore, during the initialization not just one, but two time
! records should be read, so back step "irec" one more time, however
! without performing the check for periodicity cycling: it is OK when
! irec=0 as the outcome of the second back step because the subsequent
! call to advance_record will make it 1 again.

        if (iic == ntstart) irec=irec-1

      else
        write(*,1) varname(1:lvar), irec, nf_strerror(ierr)
      endif

  1   format(/1x, '### ERROR: find_record :: Cannot read variable ''',
     &                               A, ''' for record ', I4/12x,A/)
      end

! Note that the following routine EITHER, (1) advance record number
! "irec" if possible, while NOT touching error status "ierr" at all
! OR, (2) set "ierr", while leaving "irec" not touched.

      subroutine advance_record(irec, ntime,  icycle, cycle_length,
     &                                     ifile, max_files,  ierr)
      implicit none
      integer irec, ntime, icycle, ifile, max_files, ierr
      real cycle_length
      if (irec < ntime) then
        irec=irec+1                !--> advance record normally
      elseif (cycle_length > 0.) then
        irec=1 ; icycle=icycle+1   !--> rewind time within the file
      elseif (ifile < max_files) then
        irec=-1                    !--> signal to proceed to next file
      else
        write(*,'(/1x,2A/)')  '### ERROR: advance_record :: run out ',
     &                     'of time records in non-recycling regime.'
        ierr=-9999
      endif
      end
