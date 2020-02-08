      subroutine get_2D_fld(ierr_ext)

! A generic prototype routine (literally stands for "get something")
! to read a named 2D forcing field from a netCDF file.  For a specific
! forcing data the following variable names must be redefined by CPP to
! to match the variable names within the designated part of "forces.h"
! accessed via macro SMTH_DATA, as well as  variable names within the
! netCDF file to be read.   This file is intended to be used as an
! include file and CANNOT BE COMPILED independently as a stand-alone
! file because all the to-be-redefined variables rely on declarations
! inside "forces.h" and will remain undefined unless SMTH_DATA is
! replaced with a meaningful CPP-macro matching a pre-existing macro
! in "forces.h".

! Variables to-be-redefined are:
!---------- --------------- -----
!    get_smth   name of the routine;
!     varname   name of the variable (field) inside netCDF file;
!  [ varname2,3,4 ] optional, name(s) of the second, third, and fourth
!                      variable (fields) sharing timing with the first;
!  [ timename ] name of the timing variable inside netCDF file,
!               optional: if left undefined, a code is generated to
!               (1) assume that it is 'varname' appended by suffix
!               '_time' and check whether such variable exists; if
!               not, (2) make an assumption that netCDF file follows
!               CF-conventions, so the name of the timing variable
!               should be the same as the name of time dimension which
!               in its turn is the last dimension of the variable
!               'varname'.

!   SMTH_DATA   a macro to access the specific portion of "forces.h"
!               containing relevant variables and indices; All of the
!               following macros must be defined to match specific
!               variable names from that portion of "forces.h";

!         var   name of array from "forces.h" to hold the data;
!  [ var2,3,4 ] optional, name(s) of arrays from "forces.h" to hold
!                     data for the second, third, and fourth field(s);
!     timevar   name of the associated timing variable in "forces.h";
!     file_id   index within array of netCDF file IDs of forcing files
!                   which containing data for varname, [varname2,3.4];
!       varid   name of the variable to hold netCDF variable ID for
!                               'varname' (must exist in "forces.h");
! [ varid2,3,4 ] optional, variable IDs for the second, third, and
!                                      fourth fields;

!   timevarid   variable ID for the associated timing variable;
!      record   record index for data to be read;
!
!      ntsmth   number of time records for the variable in netCDF file;
!      itsmth   time index index of array var(:,:,itsmth);
!      icycle   cycling index (number of whole cycles to add to data
!               time to make model time be within data time bounds);
! cycle_length  cycle length;

! [ VAR_GRD_TYPE,2,3,4  ]  optionally redefine the staggered grid
!               placement for the field. If not defined RHO-type is
!               assumed;  meaningful settings are r_var, u_var,
!               v_var, and q_var;


! Note: if varname2 is defined, var2,varid2 must be defined as well;
! same applies to varname3 and 4.

! Note that macro SMTH_DATA is undefined within "forces.h", while
! all others are undefined at the end of this file.  This makes it
! possible to safely include this file multiple times to generate
! routines for reading different forcing.


! Algorithm:    Scan all the forcing files one-by-one in order to
!------------   find the first file with matching variable "varname".
! Once found look for the timing variable "timename" within the same
! file; if found determine its total number of records, cycling period
! (if exists), and read the "timename" records one-by-one in order to
! find two records which bound current model time
!
!         data_time(rec-1)  <=  model_time  <  data_time(rec)
!
! (this may involve adding an integer number of cycling periods to data
! time if cycle length is positive). Once the bounding time records are
! found proceed with reading the actual data fields.  In the event when
! both "varname" and "timename" are found, but time bounding condition
! cannot be achieved, proceed with examining the next file in the
! sequence, and so on until the variables are found and the condition
! is satisfied (except when the file is the last in the sequence, so
! there is nothing left to examine).  If more than one file in the
! sequence contains the variable along with data times satisfying
! bounding condition, only the first find will be processed, while
! the second (and so) will never be looked at, so THE ORDER OF FILES
! in the sequence DOES MATTER.  In addition to that, if files in the
! sequence are alike (containing both "varname" and "timename"), but
! for the different times, this algorithm can seamlessly read the
! consecutive time records across the sequence by re-initializing
! search for time record once model time goes out of bounds by
! exceeding the time of the last record in the previous file.

! De facto the inputs are:
!-- ------ --- ------ ----
!      frcfile(:)  -- a sequence of forcing file names;
!      varname, timename -- netCDF variable names for the data
!                                     and its timing variable;
!      time -- current model time;

! outputs:

!      itsmth, timevar -- time index (1 or 2) and data time such that
!                         current model time is bounded as
!                         timevar(3-itsmth) < time <= timevar(itsmth)
!      var -- two time slices of 2D-data corresponding to these times;
!
!      file_id         -- are essentially for internal book-keeping
!      varid,             needed by this routine itself, but not to
!      timevarid          be accessed from outside.  Basically when
!      record             this routine is called again and again it
!      ntsmth             needs to do just a single-record read of
!      icycle             time and the 2D field to keep them current
!      cycle_length       with the model time.


! Error signalling:   The incoming argument (IN/OUT) "ierr_ext"
!------ -----------   is external error status which will remain
!                     unchanged and NOT TOUCHED by this routine
!                     unless
!                            (i) its incoming status is no-error;
!                     AND
!                            (ii) an error occurs inside this
!                                 routine itself.
      implicit none
      integer ierr_ext
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#ifdef SSH_DATA
# include "climat.h"
#else
# include "forces.h"
#endif
#include "netcdf.inc"
      character(len=max_name_size) frcname
      integer ncid, ifrc, ierr, lfrc, lenstr, check_vardims, ncdf_read
      real(kind=8) dtime
      character(len=40) fmt
#ifndef timename
      character(len=16) tname
      integer ndims, dimids(8), ltvr

! Note that if the name of time variable is not known, it will be
! determined from netCDF file and reported as the actual name during
! the initialization.  Thereafter the variable will be accessed by its
! netCDF variable ID and its name will never be looked at again during
! the entire period of the file use.  Nor it will be retained between
! the calls to this subroutine.  Because of this, it is initialized
! here to a trivial name rather than left uninitialized, so the
! subsequent after initialization calls print it as 'time =' rather
! than risk an attempt to access an uninitialized variable.

       tname='time' ; ltvr=4
#endif

! Initialization:  Check, whether forcing netCDF file is already in
!----------------  open state, open it if not.   Find and save netCDF
! IDs the variable itself and for the associated timing variable;
! Determine whether there is cycling to reuse the input data and what
! is cycling period "cycle_length", set initial cycling index "icycle",
! record index "record", time index "itsmth".

#ifdef VERBOSE
      write(*,*) 'entering get_smth' MYID
#endif

! Note that these three must be initialized here no matter what: not
! doing it causes segmentation fault and uninitialized variable errors
! in do while(timevar(itsmth)...)...) loop  below in the case when
! something goes wrong with reading data, even thought the
! (... .and. ierr == nf_noerr)  logic should prevent execution of this
! while.  Set both bounding times to unrealistically large negative
! values to trigger the while-loop.

      if (iic == ntstart) then
        itsmth=1  ; timevar(1)=-1.E+33
        record=-1 ; timevar(2)=-1.E+33
      endif
                                ! Note: negative record is used here
      ifrc=max(file_id,1)       ! to signal that time records bounding
 100  ncid=ncfrc(ifrc)          ! model time are not found yet.  It can
      frcname=frcfile(ifrc)     ! occur for two reasons: initialization
      lfrc=lenstr(frcname)      ! and switch-over of consecutive files
                                ! in time sequence.
      ierr=nf_noerr
      if (record < 0) then      !<-- initialization or file switch
        if (ncid == -1) then
          ierr=nf_open(frcname(1:lfrc), nf_nowrite, ncid)
          if (ierr == nf_noerr) then
            ncfrc(ifrc)=ncid
          else
            write(*,'(/1x,4A/12x,A/)')  '### ERROR: get_smth :: ',
     &                     'Cannot open file ''', frcname(1:lfrc),
     &                     ''' for reading.',   nf_strerror(ierr)
          endif
        endif
        if (ierr == nf_noerr) then
          ierr=nf_inq_varid(ncid, varname, varid)
          if (ierr == nf_noerr) then
            file_id=ifrc
#ifdef varname2
            ierr=nf_inq_varid(ncid, varname2, varid2)
            if (ierr == nf_noerr) then
# ifdef varname3
              ierr=nf_inq_varid(ncid, varname3, varid3)
              if (ierr == nf_noerr) then
#  ifdef varname4
                ierr=nf_inq_varid(ncid, varname4, varid4)
                if (ierr == nf_noerr) then
#  endif
# endif
#endif

! Determine netCDF variable ID for timing variable: if its name is
! defined by CPP-macro just hardwire it;  otherwise attempt to guess
! it by appending primary variable name with suffix '_time' and check
! whether the matching variable exists.   It this fails, make a second
! attempt by assuming that the file is "CF-compliant", so the name of
! timing variable should match the name of time dimension (as it should
! be a "coordinate variable" in CF terms), while the time dimension
! itself is expected to be the last dimension of the primary variable.

#ifdef timename
                  ierr=nf_inq_varid(ncid, timename, timevarid)
#else
                  tname=varname/ /'_time' ; ltvr=lenstr(tname)
                  ierr=nf_inq_varid(ncid, tname(1:ltvr), timevarid)
                  if (ierr /= nf_noerr) then
                    ierr=nf_inq_varndims(ncid, varid, ndims)
                    if (ierr == nf_noerr) then
                      ierr=nf_inq_vardimid(ncid, varid, dimids)
                      if (ierr == nf_noerr) then
                        ierr=nf_inq_dimname(ncid, dimids(ndims), tname)
                        if (ierr == nf_noerr) then
                          ltvr=lenstr(tname)
                          ierr=nf_inq_varid(ncid, tname(1:ltvr),
     &                                               timevarid)
# ifdef VERBOSE
                          if (ierr == nf_noerr) then
                            write(*,'(6x,6A)')   'get_smth :: ',
     &                           'identified ''', tname(1:ltvr),
     &                            ''' as time for ''', varname,
     &                                ''' using CF convention'
                          endif
# endif
                        endif
                      endif
                    endif
                  endif
#endif
                  if (ierr == nf_noerr) then
                    call find_record(ncid, timevarid,  record, ntsmth,
     &                   icycle, cycle_length,  ifrc, max_frc,  ierr)
                    if (record < 0) then
c**                   ierr=nf_close(ncid) ; ncfrc(ifrc)=-1
                      ifrc=ifrc+1 ; goto 100  !--> to the next file
                    endif
#ifndef VAR_GRD_TYPE
# define VAR_GRD_TYPE r_var
#endif
                    if (ierr == nf_noerr) then
                      ierr=check_vardims(ncid, frcname, varid,
     &                                         VAR_GRD_TYPE, 0)
#ifdef varname2
# ifndef VAR_GRD_TYPE2
#  define VAR_GRD_TYPE2 r_var
# endif
                      if (ierr == nf_noerr) then
                        ierr=check_vardims(ncid, frcname, varid2, 
     &                                          VAR_GRD_TYPE2, 0)
# ifdef varname3
#  ifndef VAR_GRD_TYPE3
#   define VAR_GRD_TYPE3 r_var
#  endif
                        if (ierr == nf_noerr) then
                          ierr=check_vardims(ncid, frcname, varid3,
     &                                            VAR_GRD_TYPE3, 0)
#  ifdef varname4
#   ifndef VAR_GRD_TYPE4
#    define VAR_GRD_TYPE4 r_var
#   endif
                          if (ierr == nf_noerr) then
                            ierr=check_vardims(ncid, frcname, varid4,
     &                                             VAR_GRD_TYPE4, 0)
#  endif
# endif
#endif
                            if (ierr == nf_noerr) then
#ifdef MPI_SILENT_MODE
                              if (mynode == 0) then
#endif
                    write(*,'(1x,7A,I4,A,I5)') 'Found ''', varname
#ifdef varname2
     &                                            / /','/ /varname2
# ifdef varname3
     &                                            / /','/ /varname3
#  ifdef varname4
     &                                            / /','/ /varname4
#  endif
# endif
#endif
     &                                    ,''' and its time ''',
#ifdef timename
     &                                  timename,
#else
     &                                  tname(1:ltvr),
#endif
     &          ''' in ''', frcname(1:lfrc), ''' ifrc =', ifrc,
     &                                     '  rec =', record+1
#ifdef MPI_SILENT_MODE
                              endif
#endif
                            endif
#ifdef varname2
# ifdef varname3
#  ifdef varname4
                          endif
#  endif
                        endif
# endif
                      endif
#endif
                    endif

                  else
#ifdef timename
                    write(*,3) timename,      varname, frcname(1:lfrc)
#else
                    write(*,3) tname(1:ltvr), varname, frcname(1:lfrc)
#endif
                  endif
#ifdef varname2
# ifdef varname3
#  ifdef varname4
                else
                  write(*,2) varname4, frcname(1:lfrc)
                endif
#  endif
              else
                write(*,2) varname3, frcname(1:lfrc)
              endif
# endif
            else
              write(*,2) varname2, frcname(1:lfrc)
            endif
#endif
          elseif (ifrc < max_frc) then
#ifdef VERBOSE
# ifdef MPI_SILENT_MODE
            if (mynode == 0) then
# endif
              write(*,1) varname, frcname(1:lfrc)
# ifdef MPI_SILENT_MODE
            endif
# endif
#endif
c**         ierr=nf_close(ncid) ; ncfrc(ifrc)=-1
            record=-1 ; ifrc=ifrc+1 ; goto 100   !--> next file
          else
            write(*,2) varname, frcname(1:lfrc)
          endif
        endif
      endif  !<-- iic == ntstart

#ifdef VERBOSE
  1   format(1x,   'WARNING: get_smth :: Cannot find variable ''',
     &                             A, ''' in file ''', A, '''.')
#endif
  2   format(1x, '### ERROR: get_smth :: Cannot find variable ''',
     &                             A, ''' in file ''', A, '''.')

  3   format(1x, '### ERROR: get_smth :: Cannot identify ''', A,
     &                        ''' as timing variable for ''', A,
     &                               ''' in file ''', A, '''.')


! Read data from the file: Check if model time is bounded by the past
!----- ---- ---- --- ----- and the future data times: if not, increment
! record and cycling indices, flip time index and read a new portion of
! the data. Repeat until model time falls between the two data times.

      do while( timevar(itsmth) < time+0.5*dt .and. ierr == nf_noerr )
        call advance_record(record, ntsmth,  icycle, cycle_length,
     &                                       ifrc, max_frc,  ierr)
        if (record < 0) then
c**       ierr=nf_close(ncid) ; ncfrc(ifrc)=-1
          ifrc=ifrc+1         ; goto 100   !-->  next file

        else if (ierr == nf_noerr) then
          ierr=nf_get_var1_double(ncid, timevarid, record, dtime)
          if (ierr == nf_noerr) then
            itsmth=3-itsmth
            timevar(itsmth)=dtime*day2sec + cycle_length*icycle

            ierr=ncdf_read(ncid, varid, record, VAR_GRD_TYPE,
     &                           var(START_2D_ARRAY,itsmth), 1)
            if (ierr == nf_noerr) then
#ifdef varname2
              ierr=ncdf_read(ncid, varid2, record, VAR_GRD_TYPE2,
     &                           var2(START_2D_ARRAY,itsmth), 1)
              if (ierr == nf_noerr) then
# ifdef varname3
                ierr=ncdf_read(ncid, varid3, record, VAR_GRD_TYPE3,
     &                             var3(START_2D_ARRAY,itsmth), 1)
                if (ierr == nf_noerr) then
#  ifdef varname4
                  ierr=ncdf_read(ncid, varid4, record, VAR_GRD_TYPE4,
     &                               var4(START_2D_ARRAY,itsmth), 1)
                  if (ierr == nf_noerr) then
#  endif
# endif
#endif
#ifdef MPI_SILENT_MODE
                    if (mynode == 0) then
#endif
                      write(fmt,'(A,I2,A)') '(12x,2A,T',
#ifdef timename
     &                                        60-len(timename),
#else
     &                                        60-ltvr,
#endif
     &                                  ',2A,F12.4,2(1x,A,I5))'

                      write(*,fmt) 'get_smth :: read ', varname
#ifdef varname2
     &                                         / /','/ /varname2
# ifdef varname3
     &                                         / /','/ /varname3
#  ifdef varname4
     &                                         / /','/ /varname4
#  endif
# endif
#endif
#ifdef timename
     &                     , timename,
#else
     &                     , tname(1:ltvr),
#endif
     &                      ' =', dtime, 'rec =', record MYID
#ifdef MPI_SILENT_MODE
                    endif
#endif
#ifdef varname2
# ifdef varname3
#  ifdef varname4
                  else
                    write(*,4) varname4, record, nf_strerror(ierr)
                  endif
#  endif
                else
                  write(*,4) varname3, record, nf_strerror(ierr)
                endif
# endif
              else
                write(*,4) varname2, record, nf_strerror(ierr)
              endif
#endif
            else
              write(*,4) varname, record, nf_strerror(ierr)
            endif
          else
#ifdef timename
            write(*,4) timename,      record, nf_strerror(ierr)
#else
            write(*,4) tname(1:ltvr), record, nf_strerror(ierr)
#endif
          endif
        else
          write(*,'(/2(1x,A,I5)/8x,4A/7x,2(A,G12.4)/)')
     &    '### ERROR: get_smth :: requested time record ',  record,
     &    'exceeds the last record', ntsmth, 'available in forcing ',
     &    'netCDF file ''',   frcname(1:lfrc),   '''',   'tdays = ',
     &     tdays,        '  but the last available  time =',
     &                                      timevar(itsmth)*sec2day
        endif
      enddo
      if (ierr_ext == nf_noerr) ierr_ext=ierr

  4   format(/1x,'### ERROR: get_smth :: Cannot read variable ''',
     &                            A, ''' for record ', I4/12x,A/)
#ifdef VERBOSE
      write(*,*) '...leaving get_smth' MYID
#endif
      end

#undef get_smth
#undef file_id
#ifdef timename
# undef timename
#endif
#undef timevar
#undef timevarid
#undef record
#undef ntsmth
#undef itsmth
#undef icycle
#undef cycle_length
#undef varname
#undef varid
#undef var
#undef VAR_GRD_TYPE
#ifdef varname2
# undef varname2
# undef varid2
# undef var2
# undef VAR_GRD_TYPE2
#endif
#ifdef varname3
# undef varname3
# undef varid3
# undef var3
# undef VAR_GRD_TYPE3
#endif
#ifdef varname4
# undef varname4
# undef varid4
# undef var4
# undef VAR_GRD_TYPE4
#endif
