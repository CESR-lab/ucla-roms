! This package contains a complete set of basic operators for reading
! and writing ROMS-style netCDF data files. The "standard" implies that
! horizontal dimensions are named as "xi_","eta_", vertical "s_" with
! corresponding suffix, "rho", "u", and "v" for horizontal dimensions;
! "rho" and "w" for vertical consistently with grid staggering rules
! within ROMS code. Time dimension (whether or not is "unlimited" from
! netCDF point of view) has its name ending with "time".

! Other than the spatial grid staggering rules, all other aspects
! related to netCDF file structure are expected to follow so called
! "CF conventions" as closely as possible.   As the result, all get_*
! and put_* routines from this package are known to work for files
! other than ROMS-standard (leaving only 4 of them "read_roms_grid",
! "write_roms_grid", "roms_find_dims", and "roms_check_dims" be
! strictly ROMS specific).

! It should be noted that somewhat similar functionality for reading
! and writing netCDF files can be found in "nc_read_write.F", however
! the distinction between there and the routines in this package is
! that all the ones below having argument "ncid" are expected to have
! netCDF file in open state, while access to a specific variable is
! done by name, hence it is "file-by-ID -- var-by-name" semantics,
! while "nc_read_write.F" uses "file-by-name -- var-by-name". For this
! reason argument "ncid" is always placed before the filename, and the
! latter is used only to write error messages, but has no effect other
! than that.

! Another note is that FORTRAN 2003 standard mandates that the rank
! of argument (scalar vs. array) should be the same for both calling
! routine and the callee, even in the trivial case where the array
! consists of just a single element.   Thus, it is formally illegal
! (thought works correctly in practice) to pass a scalar as an
! argument to a routine expecting an array of size 1.  The fact
! that size is equal to 1 in known only during runtime, but not at
! compiling, so the compiler instrumented to verify F2003 compliance
! issues an error message and quits.  It is for this and only this
! reason routines containing _sclr_ in their names in the list below
! were introduced, even thought their functionality may seem to be
! redundant (below "value" is scalar, while "var" is array).

! The content is:

!    init_time(ncid, fname, tname, nrecs, init_year, ierr)

!    read_roms_grid  (fname, Lm,Mm)
!    write_roms_grid (fname, Lm,Mm)
!    roms_find_dims  (ncid, fname, Lm,Mm,N)
!    roms_check_dims (ncid, fname, Lm,Mm,N)

!    put_sclr_by_name_real (ncid, vname, value)
!    get_sclr_by_name_real (ncid, vname, value)
!    put_sclr_by_name_double (ncid, vname, value)
!    get_sclr_by_name_double (ncid, vname, value)

!    put_var_by_name_real (ncid, vname, var)
!    get_var_by_name_real (ncid, vname, var)
!    put_var_by_name_double (ncid, vname, var)
!    get_var_by_name_double (ncid, vname, var)

!    put_sclr_rec_by_name_real   (ncid, fname, vname, rec, value)
!    get_sclr_rec_by_name_real   (ncid, fname, vname, rec, value)
!    put_sclr_rec_by_name_double (ncid, fname, vname, rec, value)
!    get_sclr_rec_by_name_double (ncid, fname, vname, rec, value)

!    put_rec_by_name_real   (ncid, fname, vname, n1,n2,n3,rec, var)
!    get_rec_by_name_real   (ncid, fname, vname, n1,n2,n3,rec, var)
!    put_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)
!    get_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)

!    put_patch_by_name_real   (ncid, fname, vname, iwest,jsouth,
!                                            n1,n2,n3,rec, var)
!    get_patch_by_name_real   (ncid, fname, vname, iwest,jsouth,
!                                            n1,n2,n3,rec, var)
!    put_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
!                                            n1,n2,n3,rec, var)
!    get_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
!                                            n1,n2,n3,rec, var)

! All are subroutines designed to provide sufficient diagnostic
! messages and terminate the execution if something goes wrong
! rather than functions returning non-zero status.

! With the exception of "read_roms_grid" which opens the named file,
! creates its netCDF file ID as an internal variable, reads all the
! relevant data, and closes it after that, all the above procedures
! imply that the file is in opened state, hence input argument "ncid"
! has meaningful value at entry, while argument "fname" is used only
! for error messages id something goes wrong.

#ifndef WRITER
# ifndef DBPRCS
#  define TIMING
#  undef VERBOSE

! The following module is designed to be completely initialized by
! "read_roms_grid", which includes both allocation arrays with proper
! dimensions matching the actual grid file and filling them with data.
! Note that "angle" is no longer part of the module because all what
! is needed in most cases when grid file is read is cos and sin of
! angle to rotate vector componets, but not angle itself.

      module roms_grid_vars
       real(kind=8), allocatable, dimension(:,:) :: lon_r,lat_r,
     &            lon_p,lat_p, pm,pn, f, rmask, csA,snA, orterr
       logical curv_grid
      end module roms_grid_vars

! The following module contains parameter for grid generation and is
! shared by only by "read_grid_params.F" (where they are initialized),
! "grid.F", and "read_roms_grid" (later in this file, used only as
! input to leave signature in netCDF file).

      module roms_grid_params
       real(kind=8) lon_min,lon_max, lat_min,lat_max,
     &    size_x,size_y, cent_lat, psi0,theta0,alpha
       integer nx,ny, flip_xy
      end module

! This module is to monitor the total amount of data read and written.

      module mod_io_size_acct
       integer(kind=8), save :: read_clk=0,  sz_read_acc=0,
     &                          write_clk=0, sz_write_acc=0
      end module mod_io_size_acct



      subroutine init_time(ncid, fname, tname, nrecs, units,
     &                                     init_year, ierr)

! Takes netCDF ID "ncid" of a file in open state and name of timing
! variable "tname" which is expected to have a single dimension (not
! necessarily having the same named as the variable itself (as in CF-
! compliant case) and finds the number of records "nrecs" available
! in the file, time units (days, hours, seconds, etc), and the initial
! year from which the time is counted. Basically it expects the timing
! variable to have attribute "units" which looks like
!
!          swf_time:units = "days since 2009-01-01 00:00:00" ;
!
! and just reads it.

      implicit none
      character(len=*) fname, tname, units
      integer ncid, nrecs, init_year, ierr,  varid, vartype, vardims,
     &   vardimids(8), varatts, i,is,ie, lfnm, ltnm,lunt,lstr, lenstr
      character(len=64) str, varname, dimname

#  include "netcdf.inc"

      nrecs=-1     ; lfnm=lenstr(fname)
      init_year=-1 ; ltnm=lenstr(tname)

      ierr=nf_inq_varid(ncid, tname(1:ltnm), varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid,  varname,  vartype,
     &                      vardims, vardimids, varatts)
        if (ierr == nf_noerr) then
          if (vardims == 1) then
            ierr=nf_inq_dim(ncid, vardimids(vardims), dimname, nrecs)
            if (ierr == nf_noerr) then
              ierr=nf_get_att_text(ncid, varid, 'units', str)
              if (ierr == nf_noerr) then

! Decode 'units' attribute: the first word is expected to be units
! as such: days, hours, seconds, the initial year is expected to be
! the first set of four digits.

                lstr=lenstr(str)
                i=1
                do while(str(i:i) /= ' ' .and. i < lstr)
                  i=i+1
                enddo
                if (str(i:i) == ' ') i=i-1
                units=str(1:i) !<-- the first word in the string
                lunt=i

                do while(i<lstr .and. (str(i:i)<'0'.or.'9'<str(i:i)))
                  i=i+1
                enddo
                if ('0' <= str(i:i) .and. str(i:i) <= '9') then
                  is=i ; ie=0
                  do while('0'<=str(i:i).and.str(i:i)<='9'.and.i<lstr)
                    i=i+1
                  enddo
                  ie=i-1
                  if (ie == is+3) then
                    init_year=0
                    do i=is,ie
                      init_year = 10*init_year +ichar(str(i:i))-48
                    enddo
                  endif
                endif

                write(*,'(1x,A,2(I6,1x,3A))')  'init_time :: found',
     &                    nrecs, 'records, units = ''', units(1:lunt),
     &                    ''' starting from year', init_year, 'in ''',
     &                                           fname(1:lfnm), '''.'
                if (init_year < 0) then
                  write(*,'(/1x,6A/)')    '### ERROR: init_time :: ',
     &            'Cannot find initial year segment within attribute ',
     &                    '''units'' for variable ''', tname(1:ltnm),
     &                            ''' in file ''', fname(1:lfnm), '''.'

                  write(*,*) 'units =''', str(1:lstr), ''''
                  ierr=ierr-999
                endif
              else
                write(*,'(/1x,6A/12x,A/)')   '### ERROR: init_time ',
     &            ':: Cannot find attribute ''units'' for variable ''',
     &                tname(1:ltnm),  ''' in file ''',  fname(1:lfnm),
     &                                      '''.',  nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,2A,I3,1x,3A/12x,A/)')        '### ERROR: ',
     &             'init_time :: Cannot determine name and size of ',
     &             'dimension #',  vardimids(vardims), ''' in file ''',
     &                         fname(1:lfnm), '''.', nf_strerror(ierr)
            endif
          elseif (vardims > 1) then
            write(*,'(/1x,6A,I3/)')       '### ERROR: init_time :: ',
     &        'Ambiguous: variable  ''',   tname(1:ltnm),   ''' in ''',
     &         fname(1:lfnm),''' has more than one dimension:',vardims
            ierr=ierr+1
          else
            write(*,'(/1x,6A/)')       '### ERROR: init_time :: No ',
     &                'dimension found for variable ''', tname(1:ltnm),
     &                                ''' in ''', fname(1:lfnm), '''.'
            ierr=ierr+1
          endif
        else
         write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: init_time ',
     &          ':: Cannot make general inquiry for variable #', varid,
     &                'in ''', fname(1:lfnm), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A/)') '### ERROR: init_time :: Cannot ',
     &              'find variable ''', tname(1:ltnm), ''' in file ''',
     &                         fname(1:lfnm), '''.', nf_strerror(ierr)
      endif
      end


      subroutine read_roms_grid(fname, Lm,Mm)

! Open ROMS grid netCDF file (argument "fname", input), read its
! dimensions (arguments Lm,Mm used in output mode), allocate all the
! arrays defined in the module above, read their values, and compute
! csA,snA, which are cos and sin of the angle between the geographical
! east direction and XI-direction of the curvilinear grid.
! For universality all the variables are allocated and read regardless
! whether they are actually needed or not.

      use roms_grid_vars
      use mod_io_size_acct
      implicit none
      character(len=*) fname
      integer Lm,Mm, ncgrd, i,j, ierr, lfnm,lenstr

#  include "netcdf.inc"

      if (.not.allocated(lon_r)) then
        lfnm=lenstr(fname)
        ierr=nf_open(fname, nf_nowrite, ncgrd)
        if (ierr == nf_noerr) then
          write(*,'(/1x,4A)') 'Reading ROMS grid file ''',
     &                             fname(1:lfnm), '''...'
          call roms_find_dims(ncgrd, fname, Lm,Mm, i)

          allocate( lon_r(0:Lm+1,0:Mm+1), lon_p(1:Lm+1,1:Mm+1),
     &              lat_r(0:Lm+1,0:Mm+1), lat_p(1:Lm+1,1:Mm+1),
     &                 pm(0:Lm+1,0:Mm+1),    pn(0:Lm+1,0:Mm+1),
     &                csA(0:Lm+1,0:Mm+1),   snA(0:Lm+1,0:Mm+1),
     &                  f(0:Lm+1,0:Mm+1), rmask(0:Lm+1,0:Mm+1),
     &                                             stat=ierr )
          if (ierr /= 0) then
            write(*,'(/1x,2A/)') '### ERROR: read_roms_grid ',
     &                           ':: Cannot allocate memory.'
            stop
          endif

          call get_var_by_name_double(ncgrd,  'lon_rho', lon_r)
          call get_var_by_name_double(ncgrd,  'lat_rho',  lat_r)
          call get_var_by_name_double(ncgrd,       'pm',     pm)
          call get_var_by_name_double(ncgrd,       'pn',     pn)
          call get_var_by_name_double(ncgrd,        'f',      f)
          call get_var_by_name_double(ncgrd, 'mask_rho',  rmask)
          sz_read_acc=sz_read_acc +  8*6*(Lm+2)*(Mm+2)

          ierr=nf_inq_varid(ncgrd, 'lon_psi', i)
          if (ierr == nf_noerr) then
            call get_var_by_name_double(ncgrd, 'lon_psi', lon_p)
            call get_var_by_name_double(ncgrd, 'lat_psi', lat_p)
            sz_read_acc=sz_read_acc + 8*2*(Lm+1)*(Mm+1)
          else
            do j=1,Mm+1
              do i=1,Lm+1
                lon_p(i,j)=0.25D0*( lon_r(i,j)+lon_r(i-1,j)
     &                        +lon_r(i,j-1) +lon_r(i-1,j-1))
                lat_p(i,j)=0.25D0*( lat_r(i,j)+lat_r(i-1,j)
     &                        +lat_r(i,j-1) +lat_r(i-1,j-1))
              enddo
            enddo
          endif

          call read_angle(ncgrd,fname, 1,1,Lm+2,Mm+2, csA,snA)
          curv_grid=.false.
          do j=0,Mm+1
            do i=0,Lm+1
              if (abs(snA(i,j)) > 1.D-12) curv_grid=.true.
            enddo
          enddo
          if (curv_grid) then
            write(*,'(1x,2A/)') 'Curvilinear grid discovered. ',
     &                     'Vector components will be rotated.'
          endif
          ierr=nf_close (ncgrd)
          return !---> successful return
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: read_roms_grid ',
     &             ':: Cannot open netCDF file ''', fname(1:lfnm),
     &                      ''' for reading,',  nf_strerror(ierr)
          stop
        endif
      else
        write(*,*) 'WARNING: Grid is already initialized.'
      endif !<-- .not.allocated(lon_r)
      end


      subroutine read_angle(ncid,fname, iwest,jsouth, isize,jsize,
     &                                                    csA,snA)

! Read a rectangular portion (subdomain) of netCDF variable "angle"
! from ROMS grid file, which the angle between local XI-coordinate of
! ROMS grid and true EAST direction and compute its sine and cosine,
! csA=cos(angle) and snA=sin(angle).  Arguments iwest,jsouth specify
! the western and southern edges, while isize,jsize -- the sizes of
! the subdomain within the i- and j-dimensions (with this respect it
! is conformal to "get_patch_by_name_double" routine defined below
! within in this file). The angle may be in either radians or degrees,
! which is determined by checking attribute "units".  NetCDF file is
! expected to be in open state ("ncid" is a valid netCDF file ID;
! argument "fname" is merely to write error messages on the screen).

      implicit none
      character(len=*) fname
      integer iwest,jsouth, isize,jsize, ncid,varid, ierr,
     &           start(2),count(2), i,j, lfnm,lstr,lenstr
      real(kind=8) csA(isize,jsize),snA(isize,jsize), cff
      character(len=16) str

#include "phys_const.h"
#include "netcdf.inc"

      lfnm=lenstr(fname)
      ierr=nf_inq_varid(ncid, 'angle', varid)
      if (ierr == nf_noerr) then
        start(1)=iwest  ; count(1)=isize
        start(2)=jsouth ; count(2)=jsize
        ierr=nf_get_vara_double(ncid, varid, start,count, csA)
        if (ierr == nf_noerr) then
          ierr=nf_get_att_text(ncid, varid, 'units', str)
          if (ierr == nf_noerr) then
            lstr=lenstr(str)
            write(*,'(2x,5A)')  'retrieved east angle from ''',
     &        fname(1:lfnm), ''', units = ''',str(1:lstr),'''.'
            if (str(1:6) == 'degree') then
C$OMP PARALLEL SHARED(isize,jsize, csA,snA) PRIVATE(i,j, cff)
C$OMP  DO
              do j=1,jsize
                do i=1,isize
                  cff=deg2rad*csA(i,j)
                  csA(i,j)=cos(cff)
                  snA(i,j)=sin(cff)
                enddo
              enddo
C$OMP  END DO
C$OMP END PARALLEL
            elseif (str(1:6) == 'radian') then
C$OMP PARALLEL SHARED(isize,jsize, csA,snA) PRIVATE(i,j, cff)
C$OMP  DO
              do j=1,jsize
                do i=1,isize
                  cff=csA(i,j)
                  csA(i,j)=cos(cff)
                  snA(i,j)=sin(cff)
                enddo
              enddo
C$OMP  END DO
C$OMP END PARALLEL
            else
              write(*,'(/1x,4A/)')    '### ERROR: Unknown units for ',
     &                'variable ''angle'' in ''', fname(1:lfnm), '''.'
              stop
            endif
          else
            write(*,'(/1x,4A/)')  '### ERROR: Cannot read attribute ',
     &    '''units'' for variable ''angle'' in ''',fname(1:lfnm),'''.'
          endif
        else
          write(*,'(/1x,4A/)')   '### ERROR: Cannot read variable ',
     &                    '''angle'' from ''',  fname(1:lfnm), '''.'
        endif
      else
        write(*,'(/1x,4A/)')  '### ERROR: Cannot find variable ',
     &                    '''angle'' in ''', fname(1:lfnm), '''.'
      endif
      if (ierr /= nf_noerr) stop
      end


      subroutine write_roms_grid(fname, Lm,Mm)

! Create ROMS grid file and write all the variables associated with
! horizontal curvilinear coordinates (these are stored inside module
! "roms_grid" where all the arrays are expected to be allocated with
! their dimensions specified by Lm,Mm and assigned meaningful values).
! This routine also creates netcdf variables for model topography and
! mask, "hraw", "h", and "mask_rho", however they are left
! uninitialized to be filled in later.

      use roms_grid_vars
      use roms_grid_params
      implicit none
      character(len=*) fname
      integer Lm,Mm, xi_rho,eta_rho, xi_u,eta_v, ncgrd, varid, ierr,
     &                r2dgrd(2), p2dgrd(2), lfnm, lstt, lstr, lenstr
      character(len=32) str
      character(len=256) settings

#  include "netcdf.inc"

      xi_rho=Lm+2   ; eta_rho=Mm+2
      xi_u=xi_rho-1 ; eta_v=eta_rho-1

! Create netcdf file.
!------- ------ -----

      lfnm=lenstr(fname)
      ierr=nf_create(fname(1:lfnm), nf_64bit_offset, ncgrd)
      if (ierr == nf_noerr) then
        ierr=nf_def_dim(ncgrd, 'xi_rho',  xi_rho,    r2dgrd(1))
        ierr=nf_def_dim(ncgrd, 'xi_u',    xi_u,      p2dgrd(1))
        ierr=nf_def_dim(ncgrd, 'eta_rho', eta_rho,   r2dgrd(2))
        ierr=nf_def_dim(ncgrd, 'eta_v',   eta_v,     p2dgrd(2))

! Grid type switch: Spherical or Cartesian.

        ierr=nf_def_var(ncgrd, 'spherical', nf_char, 0, 0, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',24,
     &                                   'grid type logical switch')
        ierr=nf_put_att_text(ncgrd,varid, 'option_T', 9, 'spherical')
        ierr=nf_put_att_text(ncgrd,varid, 'option_F', 9, 'cartesian')

! Longitude/latitude at RHO-points.

        ierr=nf_def_var(ncgrd, 'lon_rho', nf_double, 2, r2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 23,
     &                                     'longitude of RHO-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 11, 'degree_east')


        ierr=nf_def_var(ncgrd, 'lat_rho', nf_double, 2, r2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid,'long_name',22,
     &                                      'latitude of RHO-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 12,'degree_north')

! Longitude/latitude at PSI-points.

        ierr=nf_def_var(ncgrd, 'lon_psi', nf_double, 2, p2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 23,
     &                                    'longitude of PSI-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 11, 'degree_east')


        ierr=nf_def_var(ncgrd, 'lat_psi', nf_double, 2, p2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid,'long_name',22,
     &                                     'latitude of PSI-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 12,'degree_north')

! Curvilinear coordinate metric coefficients pm,pn.

        ierr=nf_def_var(ncgrd, 'pm', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 35,
     &                      'curvilinear coordinate metric in XI')
        ierr=nf_put_att_text(ncgrd,varid, 'units',  7,  'meter-1')

        ierr=nf_def_var(ncgrd, 'pn', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 36,
     &                     'curvilinear coordinate metric in ETA')
        ierr=nf_put_att_text(ncgrd,varid, 'units',  7,  'meter-1')

! Angle between direction to the EAST and XI-axis, at RHO-points

        ierr=nf_def_var(ncgrd, 'angle', nf_double, 2, r2dgrd,varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',30,
     &                           'angle between EAST and XI-axis')
        ierr=nf_put_att_text(ncgrd, varid, 'units', 7, 'degrees')

! Coriolis Parameter

        ierr=nf_def_var(ncgrd, 'f', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',   32,
     &                       'Coriolis parameter at RHO-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 8, 'second-1')


! Land-Sea mask at RHO-points.

        ierr=nf_def_var(ncgrd,'mask_rho',nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',18,
     &                                         'mask on RHO-points')
        ierr=nf_put_att_text(ncgrd, varid, 'option_0',   4,  'land' )
        ierr=nf_put_att_text(ncgrd, varid, 'option_1',   5,  'water')

! Raw and smoothed bathymetry.

        ierr=nf_def_var(ncgrd, 'hraw', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',  28,
     &                             'raw bathymetry at RHO-points')
        ierr=nf_put_att_text(ncgrd, varid, 'units', 5, 'meter')

        ierr=nf_def_var(ncgrd, 'h', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',  24,
     &                              'bathymetry at RHO-points')
        ierr=nf_put_att_text(ncgrd, varid, 'units', 5, 'meter')

        if (allocated(orterr)) then
          ierr=nf_def_var(ncgrd,'ort_error', nf_double, 2,r2dgrd,varid)
          ierr=nf_put_att_text(ncgrd, varid, 'long_name',  19,
     &                                           'orthogonality error')
          ierr=nf_put_att_text(ncgrd, varid, 'units', 7, 'degrees')
        endif

! Create signature containing parameters used for generating grid and
! save it a global attribute so the grid can be reproduced if needed.

        write(str,*) nx ; lstr=lenstr(str)
        settings='nx='/ /str(1:lstr); lstt=lenstr(settings)
        write(str,*) ny ; lstr=lenstr(str)
        settings=settings(1:lstt)/ /' ny='/ /str(1:lstr)
        lstt=lenstr(settings)

        if (lat_max > lat_min .or. lon_max > lon_min) then

          write(str,*) lat_min ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' lat_min='/ /str(1:lstr)
          lstt=lenstr(settings)

          write(str,*) lat_max ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' lat_max='/ /str(1:lstr)
          lstt=lenstr(settings)

          write(str,*) lon_min ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' lon_min='/ /str(1:lstr)
          lstt=lenstr(settings)

          write(str,*) lon_max ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' lon_max='/ /str(1:lstr)
          lstt=lenstr(settings)

        elseif (size_x>0.D0 .or. size_y>0.D0) then !--> convert to km

          write(str,*) size_x * 1.0D-3 ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' size_x='/ /str(1:lstr)
          lstt=lenstr(settings)

          write(str,*) size_y * 1.0D-3 ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' size_y='/ /str(1:lstr)
          lstt=lenstr(settings)
        endif

        if (cent_lat /= 0.D0) then
          write(str,*) cent_lat ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' cent_lat='/ /str(1:lstr)
          lstt=lenstr(settings)
        endif

        if (psi0 /= 0.D0) then
          write(str,*)  psi0 ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' Lon='/ /str(1:lstr)
          lstt=lenstr(settings)
        endif

        if (theta0 /= 0.D0) then
          write(str,*)  theta0 ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' Lat='/ /str(1:lstr)
          lstt=lenstr(settings)
        endif

        if (alpha /= 0.D0) then
          write(str,*)  alpha ; lstr=lenstr(str)
          settings=settings(1:lstt)/ /' rotate='/ /str(1:lstr)
          lstt=lenstr(settings)
        endif

        write(str,*) flip_xy; lstr=lenstr(str)
        settings=settings(1:lstt)/ /' flip_xy='/ /str(1:lstr)
        lstt=lenstr(settings)

        ierr=nf_put_att_text(ncgrd, nf_global, 'Settings',
     &                                    lstt, settings)

! Leave definition mode.
! ----- ---------- -----

        ierr=nf_enddef(ncgrd)

! Grid type switch: ALWAYS SPHERICAL

        ierr=nf_inq_varid (ncgrd, 'spherical', varid)
        ierr=nf_put_var1_text (ncgrd, varid, 1, 'T')

! Longitude/latitude at RHO- and PSI-points.

        call put_var_by_name_double(ncgrd, 'lon_rho', lon_r)
        call put_var_by_name_double(ncgrd, 'lat_rho', lat_r)
        call put_var_by_name_double(ncgrd, 'lon_psi', lon_p)
        call put_var_by_name_double(ncgrd, 'lat_psi', lat_p)

! Curvilinear coordinate metric coefficients pm,pn.

        call put_var_by_name_double(ncgrd, 'pm', pm)
        call put_var_by_name_double(ncgrd, 'pn', pn)

! Angle between XI-axis and EAST at RHO-points

c*       call put_var_by_name_double(ncgrd, 'angle', angle)

! Coriolis Parameter.

        call put_var_by_name_double(ncgrd, 'f', f)

        if (allocated(orterr)) then
          call put_var_by_name_double(ncgrd, 'ort_error', orterr)
        endif

! Close netCDF file
! ----- ------ ----

        ierr=nf_close(ncgrd)
      else
        write(*,'(/1x,4A/12x,A)')  '### ERROR: Cannot create netCDF ',
     &          'file ''',  fname(1:lfnm), '''.',  nf_strerror(ierr)
      endif
      end


! The following two routines are either to find dimensions from already
! opened netCDF file, or to check whether dimensions in the file match
! the supplied values.   The second variant is needed merely as a check
! point to detect the mismatch and to stop the execution if it happens.
! In both cases the presence of both horizontal dimensions is mandatory
! while vertical is optional (i.e., "find_dims" does not touch its
! argument N, if no dimension, so N retains its previous value or
! remains uninitialized, if not assigned by the caller; "check_dims"
! ignores vertical mismatch, if input value of "N_ck" is zero, or if
! vertical dimension cannot be determined from the file.  Both routines
! are designed to terminate execution if something goes wrong.


      subroutine roms_find_dims(ncgrd,  fname, Lm,    Mm,    N)
# else
      subroutine roms_check_dims(ncgrd, fname, Lm_ck, Mm_ck, N_ck)

#  define CHECK_DIMS
# endif /* ! DBPRCS */
      implicit none
      character(len=*) fname
      integer ncgrd, Lm,Mm,N,  xi_rho,xi_u, eta_rho,eta_v, s_rho,s_w,
     &               ndims, size, id, i,is, ierr, lvar, lfnm, lenstr
# ifdef CHECK_DIMS
     &        , Lm_ck, Mm_ck, N_ck
# endif
      character(len=16) dname
      character(len=128) string

# include "netcdf.inc"

      xi_rho=0  ; xi_u=0  ; s_rho=0
      eta_rho=0 ; eta_v=0 ;   s_w=0

      lfnm=lenstr(fname)
      ierr=nf_inq_ndims(ncgrd, ndims)
      if (ierr == nf_noerr) then
        do id=1,ndims
          dname='                '
          ierr=nf_inq_dim (ncgrd, id, dname, size)
          if (ierr == nf_noerr) then
            lvar=lenstr(dname)
            if (lvar == 6 .and. dname(1:lvar) == 'xi_rho') then
              xi_rho=size
            elseif (lvar == 4 .and. dname(1:lvar) == 'xi_u') then
              xi_u=size

            elseif (lvar == 7 .and. dname(1:lvar) == 'eta_rho') then
              eta_rho=size
            elseif (lvar == 5 .and. dname(1:lvar) == 'eta_v') then
              eta_v=size

            elseif (lvar == 5 .and. dname(1:lvar) == 's_rho') then
              s_rho=size
            elseif (lvar == 3 .and. dname(1:lvar) == 's_w') then
              s_w=size

            elseif (lvar == 5 .and. dname(1:lvar) == 'depth') then
              s_rho=size
            elseif (lvar == 7 .and. dname(1:lvar) == 'rho_ntr') then
              s_rho=size
            endif
          else
            write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     &              'determine name and size of dimension #',     id,
     &              'in ''', fname(1:lfnm), '''.', nf_strerror(ierr)
          endif
        enddo

        write(string,'(A,6(1x,A,I4))')
# ifndef CHECK_DIMS
     &       ' roms_find_dims ::',
# else
     &       'roms_check_dims ::',
# endif
     &       'xi_rho=',xi_rho,  'xi_u=',xi_u,  'eta_rho=',eta_rho,
     &        'eta_v=',eta_v,  's_rho=',s_rho,     's_w=',s_w
        lvar=lenstr(string)
        i=0                                  ! Write dimensions into
        do while(i < lvar)                   ! character string first,
          i=i+1                              ! and then suppress blank
          if (string(i:i) == '=') then       ! characters after =sign.
            i=i+1                            ! This is merely to make
            if (string(i:i) == ' ') then
              is=1
              do while(string(i+is:i+is) == ' ' .and. i+is < lvar)
                is=is+1
              enddo
              string(i:lvar-is)=string(i+is:lvar) ; lvar=lvar-is
            endif
          endif
        enddo                                ! a narrower printout
        write(*,'(2x,A)') string(1:lvar)     ! on the screen.

        ierr=0
        if (xi_rho > 0) then
          Lm=xi_rho-2
        elseif (xi_u > 0) then
          Lm=xi_u-1
        else
          write(*,'(/1x,4A/)')  '### ERROR: Cannot determine size ',
     &               'of horizontal XI-dimension in netCDF file ''',
     &                                         fname(1:lfnm), '''.'
          ierr=ierr+1
        endif
        if (eta_rho > 0) then
          Mm=eta_rho-2
        elseif (eta_v > 0) then
          Mm=eta_rho-1
        else
          write(*,'(/1x,4A/)')  '### ERROR: Cannot determine size ',
     &              'of horizontal ETA-dimension in netCDF file ''',
     &                                         fname(1:lfnm), '''.'
          ierr=ierr+1
        endif                      ! The policy here is that vertical
# ifdef CHECK_DIMS
        N=0 !<-- initialize        ! dimension is optional, therefore
# endif
        if (s_rho > 0) then        ! it is filled up only if found, and
          N=s_rho                  ! "not touched" otherwise (hence it
        elseif (s_w > 0) then      ! is possible to call this function
          N=s_w-1                  ! while passing a constant, if no
        endif                      ! vertical dimension is expected to
                                   ! exist in the file.
# ifdef CHECK_DIMS
        if (Lm /= Lm_ck) then
          write(*,'(/1x,2A,I4,1x,3A,I4/)')  '### ERROR: Size of XI-',
     &           'dimension Lm =', Lm, 'from file ''', fname(1:lfnm),
     &                 ''' does not match the previous size ', Lm_ck
          ierr=ierr+1
        endif
        if (Mm /= Mm_ck) then
          write(*,'(/1x,2A,I4,1x,3A,I4/)') '### ERROR: Size of ETA-',
     &           'dimension Mm =', Mm, 'from file ''', fname(1:lfnm),
     &                 ''' does not match the previous size ', Mm_ck
          ierr=ierr+1
        endif
        if (N > 0  .and.  N_ck > 0  .and.  N /= N_ck) then
          write(*,'(/1x,2A,I4,1x,3A,I4/)')     '### ERROR: Size of ',
     &    'vertical dimension N =', N, 'from file ''', fname(1:lfnm),
     &                  ''' does not match the previous size ', N_ck
          ierr=ierr+1
        endif
# endif
        if (ierr /= 0) stop !--> ERROR
      else
        write(*,'(/1x,4A/12x,A/)')   '### ERROR: Cannot determine ',
     &                     'number of dimensions in netCDF file ''',
     &                      fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      end
#endif  /* !WRITER */

! The rest are the reading-writting subroutins generated by CPP
! from the same source code package  (in fact, quadrupled: X2 due
! to read/write functionality and another X2 due to single/double
! precision version).

#undef KIND_TYPE
#undef nf_get_var_TYPE
#undef nf_put_var_TYPE
#undef nf_get_vara_TYPE
#undef nf_put_vara_TYPE
#undef nf_get_att_TYPE
#undef nf_put_att_TYPE
#undef get_sclr_by_name_TYPE
#undef put_sclr_by_name_TYPE
#undef get_var_by_name_TYPE
#undef put_var_by_name_TYPE
#undef get_sclr_rec_by_name_TYPE
#undef put_sclr_rec_by_name_TYPE
#undef get_rec_by_name_TYPE
#undef put_rec_by_name_TYPE
#undef get_patch_by_name_TYPE
#undef put_patch_by_name_TYPE

#ifdef DBPRCS
# define KIND_TYPE 8
# define nf_get_var_TYPE nf_get_var_double
# define nf_put_var_TYPE nf_put_var_double
# define nf_get_vara_TYPE nf_get_vara_double
# define nf_put_vara_TYPE nf_put_vara_double
# define nf_get_att_TYPE nf_get_att_double
# define nf_put_att_TYPE nf_put_att_double
# define get_sclr_by_name_TYPE get_sclr_by_name_double
# define put_sclr_by_name_TYPE put_sclr_by_name_double
# define get_var_by_name_TYPE get_var_by_name_double
# define put_var_by_name_TYPE put_var_by_name_double
# define get_sclr_rec_by_name_TYPE get_sclr_rec_by_name_double
# define put_sclr_rec_by_name_TYPE put_sclr_rec_by_name_double
# define get_rec_by_name_TYPE get_rec_by_name_double
# define put_rec_by_name_TYPE put_rec_by_name_double
# define get_patch_by_name_TYPE get_patch_by_name_double
# define put_patch_by_name_TYPE put_patch_by_name_double
#else
# define KIND_TYPE 4
# define nf_get_var_TYPE nf_get_var_real
# define nf_put_var_TYPE nf_put_var_real
# define nf_get_vara_TYPE nf_get_vara_real
# define nf_put_vara_TYPE nf_put_vara_real
# define nf_get_att_TYPE nf_get_att_real
# define nf_put_att_TYPE nf_put_att_real
# define get_sclr_by_name_TYPE get_sclr_by_name_real
# define put_sclr_by_name_TYPE put_sclr_by_name_real
# define get_var_by_name_TYPE get_var_by_name_real
# define put_var_by_name_TYPE put_var_by_name_real
# define get_sclr_rec_by_name_TYPE get_sclr_rec_by_name_real
# define put_sclr_rec_by_name_TYPE put_sclr_rec_by_name_real
# define get_rec_by_name_TYPE get_rec_by_name_real
# define put_rec_by_name_TYPE put_rec_by_name_real
# define get_patch_by_name_TYPE get_patch_by_name_real
# define put_patch_by_name_TYPE put_patch_by_name_real
#endif

! The following eight routines are just instrumented wrappers around
! the standard sequence of netCDF calls which (1) inquire variable ID
! and (2) put/get the ENTIRE variable into/from netCDF file.  These
! wrapper is needed solely to write error messages if something goes
! wrong. These are
!
!        put/get_sclr/var_by_name_TYPE (ncid, vname, value/var)
!
! where get/put and TYPE=real/double occur in all permutations (hence
! it adds up to a total of eight).  Because of semantically identical
! code real/double is implemented by CPP-redefinition of basic netCDF
! functions using the same source code.

#ifdef WRITER
      subroutine put_sclr_by_name_TYPE(ncid, vname, value)
#else
      subroutine get_sclr_by_name_TYPE(ncid, vname, value)
#endif

! These two routines are for reading or writing just a single number
! which may exist in netCDF file either as a variable or a global
! attribute containing just a single number of the proper type.
! Selection between variable or attribute is by the file, while no
! attempt to change is format is made here. For this reason is either
! variable of attribute with given name must pre-exist in order for
! these operations to succeed.

      implicit none
      integer ncid, varid, type, size, ierr, lvar, lenstr
      character(len=*) vname
      real(kind=KIND_TYPE) value
#include "netcdf.inc"

      lvar=lenstr(vname)
      ierr=nf_inq_varid(ncid, vname(1:lvar), varid)
      if (ierr == nf_noerr) then
#ifdef WRITER
        ierr=nf_put_var_TYPE(ncid, varid, value)
        if (ierr == nf_noerr) then
          write(*,'(8x,3A)') 'wrote ''', vname(1:lvar), ''''
        else
          write(*,'(/1x,4A,1x,A/)') '### ERROR: Cannot write ',
     &     'netCDF variable ''', vname, ''':', nf_strerror(ierr)
        endif
#else
        ierr=nf_get_var_TYPE(ncid, varid, value)
        if (ierr == nf_noerr) then
          write(*,'(9x,3A)') 'read ''', vname(1:lvar), ''''
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read netCDF ',
     &         'variable ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
#endif
      else
        ierr=nf_inq_att(ncid, nf_global, vname(1:lvar), type, size)
        if (ierr == nf_noerr) then

#ifdef DBPRCS
          if (size == 1 .and. type == nf_double) then
#else
          if (size == 1 .and. type == nf_real) then
#endif
#ifdef WRITER
            ierr=nf_redef(ncid)
            if (ierr == nf_noerr) then
              ierr=nf_put_att_TYPE(ncid, nf_global, vname(1:lvar),
     &                                                     value)
              if (ierr == nf_noerr) then
                ierr=nf_enddef(ncid)
                if (ierr == nf_noerr) then
                  write(*,'(8x,3A)') 'wrote ''', vname(1:lvar),
     &                                ''' as global attribute'
                else
                  write(*,'(/1x,2A/12x,A/)')     '### ERROR: Cannot ',
     &                     'close redefinition mode for netCDF file.',
     &                                              nf_strerror(ierr)
                endif
              else
                write(*,'(/1x,4A/12x,A/)')    '### ERROR: Cannot put ',
     &                            'global attribute ''', vname(1:lvar),
     &                       ''' into netCDF file.', nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,2A/12x,A/)')   '### ERROR: Cannot switch ',
     &        'netCDF file into redefinition mode.', nf_strerror(ierr)
            endif
#else
            ierr=nf_get_att_TYPE(ncid, nf_global, vname(1:lvar), value)
            if (ierr == nf_noerr) then
              write(*,'(9x,4A)') 'read ''', vname(1:lvar),
     &                           ''' as global attribute'
            else
              write(*,'(/1x,4A/12x,A/)')    '### ERROR: Cannot read ',
     &                           'global attribute ''', vname(1:lvar),
     &                      ''' from netCDF file.', nf_strerror(ierr)
            endif
#endif
          else
            if (size /= 1)  then
              write(*,'(/1x,5A,I4,1x,A/)')        '### ERROR: Global ',
     &                'attribute ''', vname(1:lvar), ''' is present, ',
     &                    'but has wrong size', size, 'instead of 1.'
              ierr=ierr-1
            endif
#ifdef DBPRCS
            if (type /= nf_double) then
#else
            if (type /= nf_real) then
#endif
              write(*,'(/1x,5A,I4,1x,I4,1x,A/)')  '### ERROR: Global ',
     &                'attribute ''', vname(1:lvar), ''' is present, ',
     &                        'but has wrong type', type, 'instead of',
#ifdef DBPRCS
     &                         nf_double, 'which is double precision.'
#else
     &                         nf_real, 'which is real.'
#endif
              ierr=ierr-1
            endif
          endif
        else
          write(*,'(/1x,4A/)')  '### ERROR: Neither variable, nor ',
     &                           'global attribute named ''', vname,
     &                              ''' is present in netCDF file.'
        endif
      endif
      if (ierr /= nf_noerr) stop
      end





#ifdef WRITER
      subroutine put_var_by_name_TYPE(ncid, vname, var)
#else
      subroutine get_var_by_name_TYPE(ncid, vname, var)
#endif
      implicit none
      integer ncid, ierr, varid, lvar,lenstr
      character(len=*) vname
      real(kind=KIND_TYPE) var(*)
#include "netcdf.inc"
      lvar=lenstr(vname)
      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
#ifdef WRITER
        ierr=nf_put_var_TYPE(ncid, varid, var)
        if (ierr == nf_noerr) then
          write(*,'(8x,3A)') 'wrote ''', vname(1:lvar), ''''
#else
        ierr=nf_get_var_TYPE(ncid, varid, var)
        if (ierr == nf_noerr) then
          write(*,'(9x,3A)')  'read ''', vname(1:lvar), ''''
#endif
        else
#ifdef WRITER
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot write ',
#else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot read ',
#endif
     &               'netCDF variable ''', vname(1:lvar), '''.',
     &                                       nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)')    '### ERROR: Cannot find ',
     &        'netCDF variable ID for ''', vname(1:lvar), '''.',
     &                                        nf_strerror(ierr)
      endif
      if (ierr /= nf_noerr) stop
      end


! The following routines are to put/get just a single number into
! a specified record of one-dimensional netCDF variable.  The netCDF
! file "fname" is expected to be open (hence input argument "ncid" is
! a valid file ID, while the name of the file is needed only to write
! error messages if something goes wrong, but otherwise is not used);
! "vname" is name of the variable (to be translated inside into netCDF
! ID with error message if not found), and "rec" is record number;
! "var" is input for put_ (output for get_) is just a scalar (single
! number). The corresponding netCDF variable is expected to be either
! a one-dimensional array (having more dimensions results in error
! message) or a scalar.  If array it puts/gets the value at location
! "rec" with performing necessary error checking; if scalar it takes
! the only value, while argument "rec" is not used.  Again, get/put
! and TYPE=real/double occur in all four permutations.


#ifdef WRITER
      subroutine put_sclr_rec_by_name_TYPE(ncid, fname, vname,
     &                                             rec, value)
#else
      subroutine get_sclr_rec_by_name_TYPE(ncid, fname, vname,
     &                                             rec, value)
#endif
      implicit none
      integer ncid, rec
      character(len=*) fname, vname
      real(kind=KIND_TYPE) value

      character(len=16) name
      integer varid, vtype, ndims, natts, dimid(8), size,
     &        start(4), count(4), ierr, lfnm, lvar, lenstr
#include "netcdf.inc"

      lfnm=lenstr(fname) ; lvar=lenstr(vname)

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then
          if (ndims == 1) then
            ierr=nf_inq_dimlen(ncid, dimid(1), size)
            if (ierr == nf_noerr) then
              start(1)=rec ; count(1)=1
#ifdef WRITER
              ierr=nf_put_vara_TYPE(ncid,varid, start,count, value)
              if (ierr == nf_noerr) then
                write(*,'(5x,A,I5,1x,5A)')     'wrote rec', rec,
     &                  'of scalar ''', vname(1:lvar), ''' into ''',
     &                                         fname(1:lfnm), '''.'
                return  !---> successful return
              else
                write(*,'(/1x,5A/12x,A)')  '### ERROR: Cannot write ',
     &                     'variable ''', vname(1:lvar), ''' into ''',
     &                         fname(1:lfnm),''':', nf_strerror(ierr)
              endif
#else
              if (0 < rec .and. rec <= size) then
                ierr=nf_get_vara_TYPE(ncid,varid, start,count, value)
                if (ierr == nf_noerr) then
                  write(*,'(6x,A,I5,1x,5A)') 'read rec', rec,
     &                    'of scalar ''', vname(1:lvar), ''' from ''',
     &                                           fname(1:lfnm), '''.'
                  return  !---> successful return
                else
                  write(*,'(/1x,7A/12x,A)') '### ERROR: Cannot read ',
     &              'variable ''', vname(1:lvar), ''' from netCDF ',
     &              'file ''',fname(1:lfnm), ''':', nf_strerror(ierr)
                endif
              else
                write(*,'(/1x,2A,I4,1x,6A,I4/)')       '### ERROR: ',
     &            'Requested record number ',   rec,  'for scalar ',
     &            'variable ''',   vname(1:lvar),  ''' in file ''',
     &             fname(1:lfnm), ''' exceeds dimension bound', size
              endif
#endif
            else
              write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     &              'Cannot determine size of dimension #',   dimid(1),
     &          'in file ''', fname(1:lfnm),  '''.', nf_strerror(ierr)
            endif
          elseif (ndims == 0) then

#ifdef WRITER
            ierr=nf_put_var_TYPE(ncid, varid, value)
            if (ierr == nf_noerr) then
              write(*,'(6x,5A)')  'wrote scalar ''', vname(1:lvar),
     &                           ''' into ''', fname(1:lfnm), '''.'
              return  !---> successful return
            else
              write(*,'(/1x,5A/12x,A)')    '### ERROR: Cannot write ',
     &            'scalar variable ''', vname(1:lvar), ''' into ''',
     &                         fname(1:lfnm),''':', nf_strerror(ierr)
            endif
#else
            ierr=nf_get_var_TYPE(ncid, varid, value)
            if (ierr == nf_noerr) then
              write(*,'(7x,5A)') 'read scalar ''', vname(1:lvar),
     &                             ''' from ''', fname(1:lfnm), '''.'
              return  !---> successful return
            else
              write(*,'(/1x,5A/12x,A)')     '### ERROR: Cannot read ',
     &              'scalar variable ''', vname(1:lvar), ''' from ''',
     &                         fname(1:lfnm),''':', nf_strerror(ierr)
            endif
#endif
          else
            write(*,'(/1x,5A,I4/)')       '### ERROR: Variable ''',
     &           vname(1:lvar), ''' from file ''', fname(1:lfnm),
     &          ''' has more then one dimension, ndims =', ndims
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)')    '### ERROR: Cannot ',
     &           'make general inquiry for variable ID =',     varid,
     &           'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     &         'for variable ''', vname(1:lvar), ''' from file ''',
     &                      fname(1:lfnm),  ''':', nf_strerror(ierr)
      endif
      stop
      end


! A more sophisticated functions, which read/write a portion of data
! into/from an existing netCDF variable accessing it by name. In doing
! so, it checks whether the first one, two, or three dimensions for
! that variable as defined in the netCDF file are consistent with the
! shape of array "var" specified as n1,n2,n3.  If not, it complains
! about the error and quits.   The variable "var" may be one, two, or
! three-dimensional, and, in addition to that may have record dimension
! (not necessarily unlimited).  If the variable is two dimensional,
! then the calling routine must be specify n3=0 to avoid checking of
! the non-existing dimension.  Similar policy applies for n1 and n2
! (e.g., n1,n2,n3=0,0,0 means that the variable is a scalar). However,
! if the variable in netCDF file has record dimension (identified as
! either unlimited netCDF dimension, or as the last dimension of the
! variable with its name ending with "...time" AND the actual number
! of spatial dimensions of netCDF variable in the file is LESS that
! the number of non-zero arguments THEN the excess spatial dimensions
! (n3, or both n3 and n2) will be ignored, while "rec" will be
! interpreted normally.

! Note: in the code below arguments n1,n2,n3 are used EXCLUSIVELY for
! checking, while the starts and counts are computed from the the
! dimensions of the variable in netCDF file.

#ifdef WRITER
      subroutine put_rec_by_name_TYPE(ncid, fname, vname,
     &                                n1,n2,n3, rec, var)
#else
      subroutine get_rec_by_name_TYPE(ncid, fname, vname,
     &                                n1,n2,n3, rec, var)
#endif
      use mod_io_size_acct
      implicit none
      integer ncid, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=KIND_TYPE), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts,  dimid(8),
     &        rec_dimid, rec_size, size, start(4),count(4),
     &                  id, ierr, lfnm, lvar, ldim, lenstr
      logical matched_dims
#include "netcdf.inc"
#ifdef TIMING
      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)
#endif
      lfnm=lenstr(fname) ;  lvar=lenstr(vname)
      do id=1,4
        start(id)=0         ! Determine the number of spatial
        count(id)=0         ! dimensions "nspc" for the variable.
      enddo
      nspc=0                ! Note that "nspc" found here should be
      if (n1 > 0) then      ! either
        nspc=nspc+1         !      less by 1 than the actual number of
        count(nspc)=n1      !      dimensions "ndims" of the variable
      endif                 !      stored in the file -- in this case
      if (n2 > 0) then      !      the extra dimension will be treated
        nspc=nspc+1         !      as record dimension;
        count(nspc)=n2      ! or
      endif                 !      be the same as "ndims" -- in this
      if (n3 > 0) then      !      case there is no record dimension,
        nspc=nspc+1         !      and argument "rec" is ignored.
        count(nspc)=n3
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then

#ifdef VERBOSE
        write(*,'(6(2x,A,I5))') 'n1 =',n1,  'n2 =', n2, 'n3 =',n3,
     &                 'rec =',rec, 'nspc =',nspc, 'ndims =',ndims
#endif

          if (nspc == ndims  .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
!      (1) unlimited dimension, or
!      (2) dimension with name ending as "...time", or
!   if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
#ifdef VERBOSE
            if (ierr == nf_noerr) then
              if (dimid(ndims) == rec_dimid) then
                write(*,'(1x,3A,I3,A,I3)')          'Variable ''',
     &            vname(1:lvar), ''' has unlimited dimension #',
     &            ndims,  ', its netCDF dimension #',  rec_dimid
              endif
            endif
#endif
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              ldim=lenstr(name)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  start(id)=1   ; rec_size=rec_size*size

                  if (count(id) /= size) then
                    ldim=lenstr(name)
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     &                                   '### ERROR: ',
#ifdef WRITER
# ifdef DBPRCS
     &                                   'put_rec_by_name_double',
# else
     &                                   'put_rec_by_name_real',
# endif
#else
# ifdef DBPRCS
     &                                   'get_rec_by_name_double',
# else
     &                                   'get_rec_by_name_real',
# endif
#endif
     &                ' :: Mismatch of dimension #', id, 'named ''',
     &                 name(1:ldim),     '''',    'for variable ''',
     &                 vname(1:lvar),  ''': attempted', count(id),
     &                'instead of',   size,   'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)')    '### ERROR: ',
     &                'Cannot get name and size of dimension #',  id,
     &                'for variable ''',    vname(1:lvar),   '''.',
     &                                             nf_strerror(ierr)
                stop
              endif
            enddo  !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
#ifdef VERBOSE
              elseif (start(ndims) == rec .and. count(ndims) == 1) then
                do id=ndims, nspc
                  write(*,'(3(1x,A),I3,1x,3A)') 'WARNING:',
# ifdef WRITER
#  ifdef DBPRCS
     &                                   'put_rec_by_name_double',
#  else
     &                                   'put_rec_by_name_real',
#  endif
# else
#  ifdef DBPRCS
     &                                   'get_rec_by_name_double',
#  else
     &                                   'get_rec_by_name_real',
#  endif
# endif
     &                        ':: Ignore spatial dimension #', id,
     &                               'for ''', vname(1:lvar), '''.'
                enddo
              else
                write(*,'(3(1x,A),2A)') 'WARNING:',
# ifdef WRITER
#  ifdef DBPRCS
     &                                   'put_rec_by_name_double',
#  else
     &                                   'put_rec_by_name_real',
#  endif
# else
#  ifdef DBPRCS
     &                                   'get_rec_by_name_double',
#  else
     &                                   'get_rec_by_name_real',
#  endif
# endif
     &              ':: Insufficient number of dimensions for ''',
     &             vname(1:lvar), '''. Record number ignored.'
#endif
              endif
#ifdef WRITER
              ierr=nf_put_vara_TYPE(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                write(*,'(5x,A,I5,1x,5A)') 'wrote rec', rec, 'of ''',
     &             vname(1:lvar), ''' into ''', fname(1:lfnm), '''.'
                sz_write_acc = sz_write_acc + rec_size * KIND_TYPE
#else
              ierr=nf_get_vara_TYPE(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                write(*,'(6x,A,I5,1x,5A)')  'read rec', rec, 'of ''',
     &             vname(1:lvar), ''' from ''', fname(1:lfnm), '''.'
                sz_read_acc = sz_read_acc + rec_size * KIND_TYPE
#endif
#ifdef TIMING
                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max
# ifdef WRITER
                write_clk = write_clk + inc_clk
# else
                read_clk = read_clk + inc_clk
# endif
#endif
                return  !---> successful return
              else
#ifdef WRITER
                write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot write ',
     &          'variable ''',vname(1:lvar),''' into netCDF file ''',
     &                       fname(1:lfnm), '''.', nf_strerror(ierr)
#else
                write(*,'(/1x,7A/12x,A)')  '### ERROR: Cannot read ',
     &          'variable ''',vname(1:lvar),''' from netCDF file ''',
     &                      fname(1:lfnm),  '''.', nf_strerror(ierr)
#endif
              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     &                 'of dimensions for variable ''', vname(1:lvar),
     &          ''': requested ', nspc, '[+1], but found in file is',
     &                                                    ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)')    '### ERROR: Cannot ',
     &           'make general inquiry for variable ID =',     varid,
     &           'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     &           'for variable ''', vname(1:lvar), ''' from file ''',
     &                      fname(1:lfnm),  ''':', nf_strerror(ierr)
      endif
      stop
      end


! The following pair is essentially an instrumented versions of
! put_vara and get_vara to read a single record of a 2D-subdomain
!         [iwest : iwest+n1-1] x [jsouth : jsouth+n2-1]
! within the (i,j)-index-space of netCDF array. The third dimension
! (if exists) is treated as the whole: n3 is checked against the
! actual dimension of netCDF variable and the mismatch is treated
! as an error; n1,n2 are checked only for upper-bound overrun, e.g.,
! iwest+n1-1 exceeds the actual size of the first dimension
! resulting in error.

#ifdef WRITER
      subroutine put_patch_by_name_TYPE(ncid, fname, vname,
     &                   iwest,jsouth, n1,n2,n3, rec, var)
#else
      subroutine get_patch_by_name_TYPE(ncid, fname, vname,
     &                   iwest,jsouth, n1,n2,n3, rec, var)
#endif
c--#define VERBOSE
      use mod_io_size_acct
      implicit none
      integer ncid, iwest,jsouth, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=KIND_TYPE), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts,  dimid(8),
     &        rec_dimid, rec_size, size, start(4),count(4),
     &                  id, ierr, lfnm, lvar, ldim, lenstr
      logical matched_dims

#include "netcdf.inc"
#ifdef TIMING
      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)
#endif
      lfnm=lenstr(fname) ;  lvar=lenstr(vname)
      do id=1,4
        start(id)=0
        count(id)=0
      enddo                 ! Determine the number of spatial
      nspc=0                ! dimensions "nspc" for the variable.
      if (n1 > 0) then
        nspc=nspc+1         ! Note that "nspc" found here should be
        start(nspc)=iwest   !
        count(nspc)=n1      ! either
      endif                 !      less by 1 than the actual number of
      if (n2 > 0) then      !      dimensions "ndims" of the variable
        nspc=nspc+1         !      stored in the file -- in this case
        start(nspc)=jsouth  !      the extra dimension will be treated
        count(nspc)=n2      !      as record dimension,
      endif                 ! or
      if (n3 > 0) then      !      be the same as "ndims" -- in this
        nspc=nspc+1         !      case there is no record dimension,
        count(nspc)=n3      !      and argument "rec" is ignored.
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then

#ifdef VERBOSE
        write(*,'(6(2x,A,I5))') 'n1 =',n1,  'n2 =', n2, 'n3 =',n3,
     &                 'rec =',rec, 'nspc =',nspc, 'ndims =',ndims
#endif
          if (nspc == ndims  .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
!      (1) unlimited dimension, or
!      (2) dimension with name ending as "...time", or
!   if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
#ifdef VERBOSE
            if (ierr == nf_noerr) then
              if (dimid(ndims) == rec_dimid) then
                write(*,'(1x,3A,I3,A,I3)')          'Variable ''',
     &            vname(1:lvar), ''' has unlimited dimension #',
     &            ndims,  ', its netCDF dimension #',  rec_dimid
              endif
            endif
#endif
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              ldim=lenstr(name)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  rec_size=rec_size*count(id)
                  if (start(id) == 0) then
                    start(id)=1
                    if (count(id) /= size) then
                      write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     &                                   '### ERROR: ',
#ifdef WRITER
# ifdef DBPRCS
     &                                   'put_patch_by_name_double',
# else
     &                                   'put_patch_by_name_real',
# endif
#else
# ifdef DBPRCS
     &                                   'get_patch_by_name_double',
# else
     &                                   'get_patch_by_name_real',
# endif
#endif
     &                ' :: Mismatch of dimension #', id, 'named ''',
     &                 name(1:ldim),     '''',    'for variable ''',
     &                 vname(1:lvar),  ''': attempted', count(id),
     &                'instead of',   size,   'in the netCDF file.'
                      stop
                    endif
                  elseif (start(id)+count(id)-1 > size) then
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,3(I5,1x,A)/)')
     &                                   '### ERROR: ',
#ifdef WRITER
# ifdef DBPRCS
     &                                   'put_patch_by_name_double',
# else
     &                                   'put_patch_by_name_real',
# endif
#else
# ifdef DBPRCS
     &                                   'get_patch_by_name_double',
# else
     &                                   'get_patch_by_name_real',
# endif
#endif
     &              ' :: Overrun dimension bound #', id, 'named ''',
     &               name(1:ldim),      '''',     'for variable ''',
     &               vname(1:lvar), ''': attempted', start(id), '+',
     &               count(id), '-1 >', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)')    '### ERROR: ',
     &                'Cannot get name and size of dimension #',  id,
     &                'for variable ''',    vname(1:lvar),   '''.',
     &                                             nf_strerror(ierr)
                stop
              endif
            enddo  !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
#ifdef VERBOSE
              elseif (start(ndims) == rec .and. count(ndims) == 1) then
                do id=ndims, nspc
                  write(*,'(3(1x,A),I3,1x,3A)') 'WARNING:',
# ifdef WRITER
#  ifdef DBPRCS
     &                                   'put_patch_by_name_double',
#  else
     &                                   'put_patch_by_name_real',
#  endif
# else
#  ifdef DBPRCS
     &                                   'get_patch_by_name_double',
#  else
     &                                   'get_patch_by_name_real',
#  endif
# endif
     &                        ':: Ignore spatial dimension #', id,
     &                               'for ''', vname(1:lvar), '''.'
                enddo
              elseif (rec > 0) then
                write(*,'(3(1x,A),2A)') 'WARNING:',
# ifdef WRITER
#  ifdef DBPRCS
     &                                   'put_patch_by_name_double',
#  else
     &                                   'put_patch_by_name_real',
#  endif
# else
#  ifdef DBPRCS
     &                                   'get_patch_by_name_double',
#  else
     &                                   'get_patch_by_name_real',
#  endif
# endif
     &              ':: Insufficient number of dimensions for ''',
     &             vname(1:lvar), '''. Record number ignored.'
#endif
              endif
#ifdef WRITER
              ierr=nf_put_vara_TYPE(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                if (start(ndims) == rec .and. count(ndims) == 1) then
                  write(*,'(5x,A,I5,1x,5A)') 'wrote rec',rec,'of ''',
     &              vname(1:lvar), ''' into ''', fname(1:lfnm), ''''
                else
                  write(*,'(5x,5A)') 'wrote''', vname(1:lvar),
     &                            ''' into ''', fname(1:lfnm), ''''
                endif
                sz_write_acc = sz_write_acc + rec_size * KIND_TYPE
#else
              ierr=nf_get_vara_TYPE(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                if (start(ndims) == rec .and. count(ndims) == 1) then
                  write(*,'(6x,A,I5,1x,5A)')  'read rec', rec,'of ''',
     &             vname(1:lvar), ''' from ''', fname(1:lfnm), ''''
                else
                  write(*,'(6x,5A)')  'read ''', vname(1:lvar),
     &                             ''' from ''', fname(1:lfnm), ''''
                endif
                sz_read_acc = sz_read_acc + rec_size * KIND_TYPE
#endif
#ifdef TIMING
                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max
# ifdef WRITER
                write_clk = write_clk + inc_clk
# else
                read_clk = read_clk + inc_clk
# endif
#endif
                return  !---> successful return
              else
#ifdef WRITER
                write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot write ',
     &          'variable ''',vname(1:lvar),''' into netCDF file ''',
     &                       fname(1:lfnm), '''.', nf_strerror(ierr)
#else
                write(*,'(/1x,7A/12x,A)')  '### ERROR: Cannot read ',
     &          'variable ''',vname(1:lvar),''' from netCDF file ''',
     &                      fname(1:lfnm),  '''.', nf_strerror(ierr)
#endif
              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     &                 'of dimensions for variable ''', vname(1:lvar),
     &          ''': requested ', nspc, '[+1], but found in file is',
     &                                                    ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)')    '### ERROR: Cannot ',
     &           'make general inquiry for variable ID =',     varid,
     &           'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     &           'for variable ''', vname(1:lvar), ''' from file ''',
     &                      fname(1:lfnm),  ''':', nf_strerror(ierr)
      endif
      stop
      end


#if !defined WRITER && !defined DBPRCS
# define WRITER
# include "roms_read_write.F"
# undef WRITER
# define DBPRCS
# include "roms_read_write.F"
# define WRITER
# include "roms_read_write.F"
#endif
