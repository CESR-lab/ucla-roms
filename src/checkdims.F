#include "cppdefs.opt"

      integer function checkdims(ncid, fname, max_rec)  ![

! Verify that all spatial dimensions in file with netCDF ID=ncid
! and file name=fname are consistent with the corresponding model
! dimensions and find the current size of unlimited dimension.
!
! input: ncid      netCDF ID and file name of input netCDF file
!        fname     [both arguments are assumed to be defined,
!                  however fname is used only to write error
!                  messages something goes wrong.]
!
! output: checkdims (returned value) error status [=nf90_noerr=0,
!                                            if no error occurs]
!         max_rec   current size of unlimited dimension
!                   [=0, if there is no unlimited dimension].

      use param
      use netcdf
      use roms_read_write
      use nc_read_write

      implicit none
      integer ncid, max_rec, ierr, icount, ndims, recdim,
     &                       dsize, id, ldim, lname, lenstr
      character(len=*) fname
      character(len=16) dname

#ifdef SOLVE3D
      integer check_scoord, chk_Cs_wr
      chk_Cs_wr=0
#endif
      lname=lenstr(fname)

      icount=0      !--> initialize: error counter (must remain 0);
      recdim=0      ! netCDF ID of unlimited dimension (if exists);
      max_rec=0     ! current size unlimited dimension (if exists).

      ierr=nf90_inquire(ncid, ndimensions = ndims)
      if (ierr == nf90_noerr) then
        ierr=nf90_inquire(ncid,unlimiteddimid = recdim)
        if (ierr/=0) call handle_ierr(ierr,'checkdim:: find unldimid ')
        do id=1,ndims
          ierr=nf90_inquire_dimension (ncid, id, dname, dsize)
          if (ierr == nf90_noerr) then
            ldim=lenstr(dname)
            if ( (ldim == 6 .and. dname(1:ldim) == 'xi_rho') .or.
     &           (ldim == 4 .and. dname(1:ldim) == 'xi_v' ) ) then
              if (dsize/=xi_rho) then
                write(*,1) dname(1:ldim), dsize,fname(1:lname),xi_rho
                icount=icount+1
              endif
            elseif ( (ldim == 4 .and. dname(1:ldim) == 'xi_u'  ) .or.
     &               (ldim == 6 .and. dname(1:ldim) == 'xi_psi') ) then
              if (dsize /= xi_u) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), xi_u
                icount=icount+1
              endif
            elseif ( (ldim == 7 .and. dname(1:ldim) == 'eta_rho') .or.
     &               (ldim == 5 .and. dname(1:ldim) == 'eta_u' ) ) then
              if (dsize/=eta_rho) then
                write(*,1) dname(1:ldim), dsize,fname(1:lname),eta_rho
                icount=icount+1
              endif
            elseif ( (ldim == 5 .and. dname(1:ldim) == 'eta_v' ) .or.
     &               (ldim == 7 .and. dname(1:ldim) =='eta_psi') ) then
              if (dsize /= eta_v) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), eta_v
                icount=icount+1
              endif
#ifdef SOLVE3D
            elseif (ldim == 3 .and. dname(1:ldim) == 's_w') then
              chk_Cs_wr=chk_Cs_wr+1
              if (dsize /= N+1) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), N+1
                icount=icount+1
              endif
            elseif (ldim == 5 .and. dname(1:ldim) == 's_rho') then
              chk_Cs_wr=chk_Cs_wr+2
              if (dsize /= N) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), N
                icount=icount+1
              endif
#endif
            elseif (id==recdim) then
              max_rec=dsize
            endif
          else
            write(*,'(/1x,2A,I3/12x,3A/12x,A/)')      '### ERROR: ',
     &            'checkdims :: Cannot get size of dimension #', id,
     &            'from netCDF file ''',    fname(1:lname),   '''.',
     &                                           nf90_strerror(ierr)
            icount=icount+1
          endif
        enddo
      else
        write(*,'(/1x,4A/12x,A/)')       '### ERROR: checkdims :: ',
     &          'Cannot get number of dimensions in netCDF file ''',
     &                     fname(1:lname), '''.', nf90_strerror(ierr)
        icount=icount+1
      endif
#ifdef SOLVE3D
      if (chk_Cs_wr>0) then
        icount=icount+check_scoord(ncid, fname, chk_Cs_wr)
      endif
#endif
  1   format(/' ### ERROR: checkdims :: wrong size of dimension ''',
     &                  A, ''' =', i5 / 12x, 'in netCDF file ''', A,
     &                            ''': must be', i5,1x, 'instead.'/)
      checkdims=icount

      end  !]

      integer function check_vardims(ncid,fname, varid, horiz_type,n3)  ![

! Verify that spatial dimensions of a specified variable in netCDF
! file match that expected by the running mode.
!
! input: ncid       netCDF ID and file name of input netCDF file
!        fname      [both arguments are assumed to be defined, however
!                   fname is used only for writing error messages if
!                   something goes wrong.]
!        varid      variable ID
!        horiz_type horizontal type defining placement of the variable
!                   in staggered grid [=0,1,2,3 matching the definition
!                       of one of the named parameters in "ncvars"]
!        n3         expected size of vertical dimension [=0 means do
!                   not check this dimension]

! returned value    error status [=nf90_noerr when no mismatches]

      use param
      use netcdf

      implicit none
      integer ncid, varid, horiz_type, n3
      character(len=*) fname
      character(len=32) vname
      integer ierr, vtype, ndims, natts, dimid(8),
     &                size, i, lfnm, lvar, lenstr

#define NO_RECORD_CHECK
#include "compute_starts_counts.h"
#if defined MPI && !defined PARALLEL_FILES
      count(1)=LLm+2-imin  !<-- dimensions of the variable
      count(2)=MMm+2-jmin  !<-- on the entire physical grid
#endif

      lfnm=lenstr(fname)
      ierr=nf90_inquire_variable(ncid, varid, vname, vtype, ndims, dimid, natts)
      if (ierr == nf90_noerr) then
        lvar=lenstr(vname)
        do i=1,2
          if (ierr == nf90_noerr) then
            ierr=nf90_inquire_dimension(ncid, dimid(i), len=size)
            if (ierr == nf90_noerr) then
              if (count(i) /= size) then
                write(*,'(/2A,I3,1x,5A,2(I5,1x,A)/)')   '### ERROR: ',
     &                    'checkvardims :: mismatch in dimension #', i,
     &                    'for variable ''', vname(1:lvar), ''' in ''',
     &                     fname(1:lfnm),  ''': expecting',   count(i),
     &                              'instead of', size, 'in the file.'
                ierr=nf90_noerr-1
              endif
            else
              write(*,'(/1x,2A,I3,1x,6A/)') '### ERROR: checkvardims ',
     &                    ':: Cannot determine size of dimension #', i,
     &                    'for variable ''', vname(1:lvar), ''' in ''',
     &                       fname(1:lfnm), ''', ', nf90_strerror(ierr)
            endif
          endif
        enddo

        if (n3 > 0 .and. ierr == nf90_noerr) then
!          ierr=nf_inq_dimlen(ncid, dimid(3), size)
          ierr=nf90_inquire_dimension(ncid, dimid(3), len=size)
          if (ierr == nf90_noerr) then
            if (n3 /= size) then
              write(*,'(/6A,2(I4,1x,A)/)')  '### ERROR: checkvardims ',
     &                    ':: mismatch in dimension #3 of variable ''',
     &                        vname(1:lvar), ''' in ''', fname(1:lfnm),
     &         ''': expecting', n3, 'instead of', size, 'in the file.'
              ierr=nf90_noerr-1
            endif
          else
            write(*,'(/1x,7A/)')  '### ERROR: checkvardims :: Cannot ',
     &                'determine size of dimension #3 for variable ''',
     &                vname(1:lvar), ''' in ''', fname(1:lfnm), ''', ',
     &                                               nf90_strerror(ierr)
          endif
        endif ! <-- n3 > 0
      else
        write(*,'(/1x,2A,I3,1x,4A/)')    '### ERROR: checkvardims :: ',
     &      'Cannot make inquiry for variable varid =', varid, 'in ''',
     &                       fname(1:lfnm), ''', ', nf90_strerror(ierr)
      endif
      check_vardims=ierr

      end  !]


#ifdef SOLVE3D

! Check whether S-coordinate parameter hc, and Cs_w- Cs_r-curves found
! in netCDF file "fname" match their counterparts in the running model.
! Print error messages and return non-zero status, in the case of
! mismatch.

      integer function check_scoord(ncid, fname, what_to_check)  ![

      use param
      use scoord

      implicit none
      integer ncid, what_to_check
      character(len=*) fname
      real, parameter :: epsil=1.D-7
      real tst_val(N+1)
      integer ierr, k, icount, lfile, lenstr, read_nc1dat
      logical chk_Cs_w, chk_Cs_r

      if (mod(what_to_check,2) == 1) then
        chk_Cs_w=.true.                    ! Decode incoming argument
      else                                 ! "what_to_check" into two
        chk_Cs_w=.false.                   ! logical switched to check
      endif                                ! "Cs_w", "Cs_r", or both.
      if (what_to_check > 1) then          ! The policy here is to
        chk_Cs_r=.true.                    ! check only curve(s) which
      else                                 ! has the corresponding
        chk_Cs_r=.false.                   ! vertical dimension
      endif

      icount=0
      lfile=lenstr(fname)
      ierr=read_nc1dat(ncid, fname, 'hc', 1, tst_val)
      if (ierr == 0) then
        if (abs((hc-tst_val(1))/hc)>epsil) then
                print *, hc,tst_val(1)
          write(*,'(/1x,2A,F12.5,1x,A,F12.5,1x,3A/)')  '### ERROR: ',
     &                  'Mismatch in ''hc'': should be', hc, 'found',
     &               tst_val(1), 'in file ''', fname(1:lfile), '''.'
          icount=icount+1
# ifdef VERBOSE
        else
          write(*,*) 'check_scoord :: checking hc:  OK'
# endif
        endif
      endif

      if (chk_Cs_w) then
        ierr=read_nc1dat(ncid, fname, 'Cs_w', N+1, tst_val)
        if (ierr == 0) then
          do k=N,0,-1
            if (abs(Cs_w(k)-tst_val(k+1)) > epsil) then
              ierr=ierr+1
            endif
          enddo
          if (ierr /= 0) then
            write(*,'(/1x,4A/)')  '### ERROR: Mismatch in ''Cs_w'' ',
     &                    'values in file ''', fname(1:lfile), '''.'
            icount=icount+1
# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking Cs_w: OK'
# endif
          endif
        endif
      endif

      if (chk_Cs_r) then
        ierr=read_nc1dat(ncid, fname, 'Cs_r', N, tst_val)
        if (ierr == 0) then
          do k=N,1,-1
            if (abs(Cs_r(k)-tst_val(k)) > epsil) then
              ierr=ierr+1
            endif
          enddo
          if (ierr /= 0) then
            write(*,'(/1x,4A/)')  '### ERROR: Mismatch in ''Cs_r'' ',
     &                    'values in file ''', fname(1:lfile), '''.'
            icount=icount+1
# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking Cs_r: OK'
# endif
          endif
        endif
      endif

      if (icount /= 0) then
        ierr=read_nc1dat(ncid, fname, 'theta_s', 1,  tst_val)
        if (ierr == 0) then
          if (abs(theta_s-tst_val(1)) > epsil) then
            write(*,'(/1x,A,2(A,F9.4,1x),3A/)')       '### ERROR: ',
     &              'Mismatch in ''theta_s'': should be',   theta_s,
     &      'found', tst_val(1), 'in file ''', fname(1:lfile), '''.'
# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking theta_s: OK'
# endif
          endif
        endif
        ierr=read_nc1dat(ncid, fname, 'theta_b', 1,  tst_val)
        if (ierr == 0) then
          if (abs(theta_b-tst_val(1)) > epsil) then
            write(*,'(/1x,A,2(A,F9.4,1x),3A/)')        '### ERROR: ',
     &      'Mismatch in ''theta_b'': should be',   theta_b, 'found',
     &               tst_val(1), 'in file ''', fname(1:lfile), '''.'
# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking theta_b: OK'
# endif
          endif
        endif
      endif

      check_scoord=icount

      end  !]
#endif



      function read_nc1dat(ncid, fname, vname, nlen, value)  ![

! A basic function to retrieve a simple variable [stored as either
! a scalar or fixed-size array], which can be stored in netCDF file
! as a global attribute or a variable.
!
! input: ncid     netCDF ID and name of input file [both arguments
!        fname    are assumed to be defined at entry (hence the file
!                 in already open state) however fname is used only
!                 to write error messages something goes wrong.]
!
!        vname    name of the variable to be read and its expected
!        nlen     length, which will be checked against the actual
!                 size of the global attribute (or variable) in
!                 the file.
!
! output: value   value(s, if array) of the variables. Note that
!                 the value is read as a whole [there is no
!                 provision to read a subset of values, a record.]
!
! returned value  error status =0 means no error.

      use param
      use netcdf

      implicit none
      integer read_nc1dat, ncid, nlen
      character(len=*) fname, vname
      real value(nlen)

      integer ndims, dimids(8), varid, vtype, size,
     &               i,j, ierr,  lvar, lfile, lenstr

      lvar=lenstr(vname) ; lfile=lenstr(fname)

      ierr=nf90_inquire_attribute(ncid, nf90_global, vname, vtype, size)
      if (ierr == nf90_noerr) then
        if (size == nlen) then
          ierr=nf90_get_att(ncid, nf90_global, vname, value)
          if (ierr /= nf90_noerr) then
            write(*,'(/1x,6A/12x,A/)')    '### ERROR: read_nc1dat :: ',
     &          'Cannot retrieve global attribute ''',   vname(1:lvar),
     &          ''' from ''', fname(1:lfile), '''.', nf90_strerror(ierr)
          endif
        endif
      else
        size=1
        ierr=nf90_inq_varid(ncid, vname, varid)
        if (ierr == nf90_noerr) then
!          ierr=nf_inq_varndims(ncid, varid, ndims)
          ierr=nf90_inquire_variable(ncid, varid, ndims=ndims)
          if (ierr == nf90_noerr) then
            if (ndims > 0) then
!              ierr=nf_inq_vardimid(ncid, varid, dimids)
              ierr=nf90_inquire_variable(ncid, varid, dimids=dimids)
              if (ierr == nf90_noerr) then
                i=0
                do while (i<ndims .and. ierr==nf90_noerr)
                  i=i+1
!                  ierr=nf_inq_dimlen(ncid, dimids(i), j)
                  ierr=nf90_inquire_dimension(ncid, dimids(i), len=j)
                  if (ierr == nf90_noerr) then
                    size=size*j
                  endif
                enddo
                if (ierr /= nf90_noerr) then
                  write(*,'(/1x,2A,I3,1x,3A/12x,A/)')    '### ERROR: ',
     &                'read_nc1dat :: Cannot find size of dimension #',
     &                 dimids(i),  'in file ''', fname(1:lfile), '''.',
     &                                             nf90_strerror(ierr)
                endif
              else
                write(*,'(/1x,6A/12x,A/)')   '### ERROR: read_nc1dat ',
     &       ':: Cannot determine dimension IDs for ''', vname(1:lvar),
     &       ''' in file ''', fname(1:lfile), '''.', nf90_strerror(ierr)
              endif
            endif
          else
            write(*,'(/1x,6A/12x,A/)')    '### ERROR: read_nc1dat :: ',
     &   'Cannot determine number of dimensions for ''', vname(1:lvar),
     &      ''' in file ''',  fname(1:lfile), '''.', nf90_strerror(ierr)
          endif
          if (size == nlen .and. ierr == nf90_noerr) then
            ierr=nf90_get_var(ncid, varid, value)
            if (ierr/=nf90_noerr) then
              write(*,'(/1x,6A/12x,A/)')  '### ERROR: read_nc1dat :: ',
     &                       'Cannot read variable ''',  vname(1:lvar),
     &     ''' from file ''', fname(1:lfile), '''.', nf90_strerror(ierr)
            endif
          endif
        else
          mpi_master_only write(*,'(/1x,6A/14x,A/)')
     &             '### WARNING: read_nc1dat :: ',
     &            'Cannot determine netCDF ID for ''',  vname(1:lvar),
     &      ''' in file ''', fname(1:lfile), '''.', nf90_strerror(ierr)
        endif
      endif
      if (size /= nlen .and. ierr == nf90_noerr) then
        write(*,'(/1x,6A/12x,A,I7,1x,A,I7/)')    '### ERROR: ',
     &   'read_nc1dat :: Unexpected size of ''', vname(1:lvar),
     &   ''' in file ''',  fname(1:lfile),  ''':', 'should be',
     &                                 nlen, 'but found', size
        ierr=+1
      endif
      read_nc1dat=ierr

      end  !]


