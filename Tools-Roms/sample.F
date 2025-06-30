#ifndef VAR_TYPE
c--# define VERBOSE

      program sample

! A generic "compile-once --> use-forever" command-line tool to extract
! a desired subdomain from a ROMS-style netCDF file or sequence of such
! files. Usage: "sample" takes at least 6 arguments, where the first 4
! are integer numbers specifying subdomain; arguments from 5 to the
! one-before-last are the name(s) of input files; and the last argument
! is always the name of output file, for example

!       sample imin imax jmin jmax source.nc target.nc
!       sample imin imax jmin jmax sr1.nc sr2.nc [...srN.nc] target.nc
!       sample imin imax jmin jmax src_*.nc target.nc

! Bounds imin, imax, jmin, jmax are defined for RHO-points using "C"
! or "ncview"- indexing convention, that is, acceptable values start
! from 0 and run to size_of_dimension-1, hence

!              0  <=  imin  <=  imax  <=  xi_rho-1
!              0  <=  jmin  <=  jmax  <=  eta_rho-1,

! so basically they can be selected by pointing at the desired pixel
! and reading its X,Y-indices when using "ncview" to visualize a RHO-
! type field ("ncview" uses C-style index range 0:dim_size-1).
! The dimensions of the resultant subdomain netCDF file will be set
! according to ROMS rules, that is one point less for U-velocity
! in XI-direction, and one point less for V- in ETA-direction, except
! in the cases when imin == imax or jmin == jmax, which is extraction
! of vertical section.  For vertical section a single row of RHO-points
! will always be extracted, however for U-(V-)points there are two
! possibilities:
!   (i) if the crossection is in the middle, two rows of U-(V-)points
!       flanking the RHO-rows are extracted;
!   (ii) if either imin == 0 or imin == xi_rho-1, only one row of U-
!        points -- the one which on the inside side from RHO-point is
!        extracted;
! This is shown on the diagram (uppercase U,R are boundary points):

!   imin=0                0  < imin < xi_rho-1            imin=xi_rho-1
!-----------                 ---------------                ---------
!   R  -  U  -  r  -  u  ...  u  -  r  -  u ...  u  -  r  -  U  -  R
!   0     0     1     1     i-1     i     i           xi-2  xi-2  xi-1
! [extracted]              [ extracted rows ]             [ extracted ]

! The structure, dimensions, and composition of variables in the output
! file will be exactly the same as for the source [taken from the first
! file of the sequence in the case of multiple files], except that
! dimensions xi_rho, xi_u, eta_rho, and eta_v are changed as determined
! by imin,imax,jmin,jmax.  Variables which do not have these dimensions
! are copied as a whole, while variables having only one of these
! dimensions (like lateral boundary data) will be copied within the
! imin:imax or jmin:jmax range of indices, whichever applicable.

! If multiple source files have different sets of variables, then only
! the first file is allowed have more variables than others, and the
! extra variables are not allowed to have record dimension: e.g., it is
! a widely used practice to put grid variables into starting history
! file of output sequence, while not into the subsequent. In this case
! if the whole sequence is given to "sample", then it puts grid
! variables into output as well.

! This code does not perform any data type conversions, it just copies
! the data exactly as is in the source file(s).

! Created by Alexander Shchepetkin: old_galaxy@yahoo.com

      implicit none
      logical var_with_recs
      integer nargs, iarg, imin, imax, jmin, jmax, ncsrc, nctarg,
     &        ndims, nvars, ngatts, vtype, varatts, vdims, newvdims,
     &        varid, newvarid, nrecords, irec, newrec, last_rec_out,
     &        unlimdimid, ierr, i, ipass, id, size, old_fill_mode,
     &        var_grid_type, lsrc, ltarg, latt, lvar, ldim,
     &        lenstr, read_str_int

      integer, allocatable, dimension(:) :: dimid, vdimids,newvdimids,
     &                                 start,count, newstart,newcount

      integer targ_size,      max_bffrc_size, max_targc_size,
     &        max_ibfr2_size, max_itrg2_size, max_ibffr_size,
     &        max_itarg_size, max_bffr4_size, max_targ4_size,
     &                        max_bffr8_size, max_targ8_size

      character(len=1), allocatable, dimension(:) :: bffrc, targc
      integer(kind=2),  allocatable, dimension(:) :: ibfr2, itrg2
      integer(kind=4),  allocatable, dimension(:) :: ibffr, itarg
      real(kind=4),     allocatable, dimension(:) :: bffr4, targ4
      real(kind=8),     allocatable, dimension(:) :: bffr8, targ8

      character(len=128) src_file, trg_fname
      character(len=64) dname, attname, vname, string
      integer, parameter :: not_found=-32768

# include "netcdf.inc"

      max_bffrc_size=0 ; max_targc_size=0    ! max_TYPE_size are to
      max_ibfr2_size=0 ; max_itrg2_size=0    ! leep tack of the actual
      max_ibffr_size=0 ; max_itarg_size=0    ! size of allocated buffer
      max_bffr4_size=0 ; max_targ4_size=0    ! arrays for each type
      max_bffr8_size=0 ; max_targ8_size=0

      imin=not_found ; imax=not_found
      jmin=not_found ; jmax=not_found

      nargs=0                           ! At first, sort out arguments
      nargs=iargc()                     ! in the command line and check
      if (nargs >= 6) then              ! whether "sample" is invoked
        do iarg=1,4                     ! correctly.  Complain about
          call getarg(iarg,string)      ! error and quit, if something
          ierr=read_str_int(string,i)   ! goes wrong.
          if (ierr == 0) then
            if (iarg == 1) then
              imin=i
            elseif (iarg == 2) then
              imax=i
            elseif (iarg == 3) then
              jmin=i
            elseif (iarg == 4) then
              jmax=i
            endif
          else
            lvar=lenstr(string)
            write(*,'(1x,A,I3,1x,3A/12x,A)') '### ERROR: Argument #',
     &               iarg, 'has illegal value of ''', string(1:lvar),
     &                         '''.', 'It must be a number instead.'
          endif
        enddo
      else
        write(*,'(/1x,2A/1x,A/ /8x,A/1x,A/8x,2A/1x,A/8x,A/)')
     &     'Command-line operator "sample" extracts subdomain from ',
     &     'ROMS-style netCDF file.', 'Its usage should be:',
     &     'sample imin imax jmin jmax source_file.nc target_file.nc',
     &     'or',     'sample imin imax jmin jmax file1.nc ',
     &                                  '[file2.nc ...] target.nc',
     &     'or',     'sample imin imax jmin jmax file*.nc target.nc'

        write(*,'(5(1x,2A/))')    'where imin,imax,jmin,jmax are ',
     &     'bounding indices defining rectangular subdomain',
     &     'using netCDF C-convention: the range of indices ',
     &     'should be from 0 to the number', 'of RHO-points-1 in the ',
     &     'corresponding horizontal direction (c.f., the indices',
     &     'are defined exactly as they reported "ncview" when ',
     &     'pointing mouse cursor to', 'the specific grid point of ',
     &                                          'a RHO-type field.)'

       write(*,'(3(1x,2A/))')   'Any number of input netCDF files is ',
     & 'acceptable, including the use of wild', 'cards, however only ',
     & 'a single output file containing all the temporal records',
     &                  'for the entire sequence ', 'will be created.'
        stop
      endif

      ierr=0
      if (imin == not_found) ierr=ierr+1
      if (imax == not_found) ierr=ierr+1
      if (jmin == not_found) ierr=ierr+1
      if (jmax == not_found) ierr=ierr+1
      if (ierr == 0) then
        if (imax < imin) then
          write(*,'(/1x,2A,2I5/)')   '### ERROR: wrong XI-bounds, ',
     &                                    'imin,imax =', imin,imax
          ierr=ierr+1
        endif
        if (jmax < jmin) then
          write(*,'(/1x,2A,2I5/)')   '### ERROR: wrong ETA-bounds, ',
     &                                    'jmin,jmax =', jmin,jmax
          ierr=ierr+1
        endif
      endif
      if (ierr /= 0) stop

      write(*,'(/1x,2A,4I4/)')     'Extracting subdomain ',
     &        'imin,imax,jmin,jmax=', imin,imax,jmin,jmax

      iarg=5
      call getarg(iarg,  src_file)    !<-- first file to process
      lsrc=lenstr(src_file)
      call getarg(nargs, trg_fname)   !<-- output file
      ltarg=lenstr(trg_fname)
      last_rec_out=0

      ierr=nf_open(src_file(1:lsrc), nf_nowrite, ncsrc)
      if (ierr == nf_noerr) then
        ierr=nf_inq (ncsrc, ndims, nvars, ngatts, unlimdimid)
        if (ierr == nf_noerr) then

          allocate( dimid(1:ndims), vdimids(1:ndims), start(1:ndims),
     &                           count(1:ndims), newvdimids(1:ndims),
     &                        newstart(1:ndims), newcount(1:ndims) )

          ierr=nf_create(trg_fname(1:ltarg),nf_64bit_offset,nctarg)
          if (ierr == nf_noerr) then
            ierr=nf_set_fill(nctarg, nf_nofill, old_fill_mode)
            if (ierr == nf_noerr) then
               write(*,'(/4x,3A/)')      'Created netCDF file ''',
     &                                    trg_fname(1:ltarg),'''.'
            else
              write(*,'(/1x,A/)') '### ERROR: Cannot set nofill mode.'
              stop
            endif

! Sort out dimensions: for each dimension find its name, size and
! create the corresponding dimension in the target file, which has
! the same order in the sequence, the same name, and the same
! length, with the exception of "xi_*" and "eta_*" dimensions which
! have their length determined from to the desired subdomain bounds.

            do id=1,ndims
              ierr=nf_inq_dim(ncsrc, id, dname, size)
              if (ierr == nf_noerr) then
                ldim=lenstr(dname)
                if (id == unlimdimid) then
                  size=nf_unlimited
                else
                  if (dname(1:3)=='xi_') then
                    if ( (dname(4:6)=='rho' .or. dname(4:4)=='v')
     &                                       .and. imax < size ) then
                      size=imax-imin+1
                    elseif ( (dname(4:6)=='psi' .or. dname(4:4)=='u')
     &                                     .and. imax < size+1 ) then
                      if (imin == imax) then
                        if (imin == 0 .or. imax == size) then
                          size=1      !<-- extract U-boundary row
                        else
                          size=2      !<-- U - R - U in the middle
                        endif
                      else
                        size=imax-imin
                      endif
                    else
                      write(string,*) size-1
                      write(*,'(/1x,2A,I6,2(1x,3A)/)')   '### ERROR: ',
     &                 'Specified imax =',   imax,    'is outside of ',
     &                 'allowed range 0:',    string(1:lenstr(string)),
     &                 'for dimension ''',dname(1:lenstr(dname)),'''.'
                      stop
                    endif
                  elseif (dname(1:4)=='eta_') then
                    if ( (dname(5:7)=='rho' .or. dname(5:5)=='u')
     &                                       .and. jmax <size ) then
                        size=jmax-jmin+1
                    elseif ( (dname(5:7)=='psi' .or. dname(5:5)=='v')
     &                                    .and. jmax < size+1 ) then
                      if (jmin == jmax) then
                        if (jmin == 0 .or. jmax == size) then
                          size=1      !<-- extract V-boundary row
                        else
                          size=2      !<-- V - R - V in the middle
                        endif
                      else
                        size=jmax-jmin
                      endif
                    else
                      write(string,*) size-1
                      write(*,'(/1x,2A,I6,2(1x,3A)/)')   '### ERROR: ',
     &                 'Specified jmax =',   jmax,    'is outside of ',
     &                 'allowed range 0:',    string(1:lenstr(string)),
     &                 'for dimension ''',dname(1:lenstr(dname)),'''.'
                      stop
                    endif
                  endif
                endif
                ierr=nf_def_dim(nctarg, dname(1:ldim), size, dimid(id))
                if (ierr == nf_noerr) then
                  write(*,'(6x,3A,I4)')     'created dimension ''',
     &                                 dname(1:ldim), ''' =', size
                else
                  write(*,'(/1x,6A/12x,A/)')     '### ERROR: Cannot ',
     &                           'create dimension ''', dname(1:ldim),
     &                 ''' within netCDF file ''', trg_fname(1:ltarg),
     &                                       ''':', nf_strerror(ierr)
                  stop
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A/)')    '### ERROR: ',
     &            'Cannot determine name and size of dimension #', id,
     &                 'in netCDF file ''',  src_file(1:lsrc),  '''.',
     &                                              nf_strerror(ierr)
                stop
              endif
            enddo  ! <-- id

! Copy global attributes

            do id=1,ngatts
              ierr=nf_inq_attname(ncsrc, nf_global, id, attname)
              if (ierr == nf_noerr) then
                latt=lenstr(attname)
                if (attname(1:latt) == 'nc_format') then
                  write(*,'(8x,3A)') 'skipping global attribute ''',
     &                                       attname(1:latt), '''.'
                else
                  ierr=nf_copy_att(ncsrc, nf_global, attname(1:latt),
     &                                            nctarg, nf_global)
                  if (ierr == nf_noerr) then
                    write(*,'(8x,3A)')  'copied global attribute ''',
     &                                        attname(1:latt), '''.'
                  else
                    write(*,'(/1x,4A/12x,3A/12x,A/)')   '### ERROR: ',
     &              'Cannot copy global attribute ''',attname(1:latt),
     &              '''', 'into netCDF file ''',   trg_fname(1:ltarg),
     &                                       '''.', nf_strerror(ierr)
                  endif
                endif
              else
                write(*,'(/1x,3A,I3/12x,3A/12x,A/)')   '### ERROR: ',
     &             'Cannot determine name of global attribute with ',
     &             'no', id, 'from netCDF file ''', src_file(1:lsrc),
     &                                      '''.', nf_strerror(ierr)
              endif
            enddo

! The rest of this program is organized as a triple nested structure,

!           do ipass=0,....  <-- loop over files (except ipass == 0)
!             do while (varid=...)  <-- loop over variables
!               --> for each variable id: examine dimensions, find
!               --> whether it is (i) a scalar (has no xi-, eta-
!               --> dimensions), and (ii) does have a time dimension
!               --> (ether unlimited dimension regardless of name,
!               --> or named as [smth_]time); calculate buffer size
!               --> needed to hold one record of it (the whole, if
!               --> it does not have time dimension);
!               do while (irec=... )<-- loop over records within file
!                 if (ipass == 0) then
!                   --> define variable in target netCDF file,
!                   --> but do not write anything into it yet.
!                 else
!                   --> read variable or record from source file
!                   --> write variable or record
!                 endif
!               enddo
!             enddo
!             if (ipass == 0) then
!               --> switch target file from definition to writing
!             else
!               --> change source file for ipass=1,2,..., nagrs-1
!             endif
!           enddo
!
! Scan through variable list and process them one-by-one. Depending
! on its dimensions, each variable is classified by two categories:
! (i)  var_grid_type=0,1,2,3 where 0 that it does not have spatial
!      dimensions xi_ and eta_ (so there is no need to cut a subdomain
!      out of this variable); 1 have xi-only; 2 is for eta-only,
!      3 for both, and
! (ii) "var_with_recs=T/F" --- which affect the manner of how the
!      variable needs to be subsequently processed, either read/write
!      it record-by-record, or the whole variable at once.


!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   *   ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     ***
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***


            do ipass=0,nargs-5
              varid=0
              do while(varid < nvars .and. ierr == nf_noerr)
                varid=varid+1
                start=0 ; var_grid_type=0
                count=0 ; var_with_recs=.false.

                ierr=nf_inq_var(ncsrc, varid, vname,  vtype,
     &                              vdims, vdimids, varatts)
                if (ierr == nf_noerr) then
                  lvar=lenstr(vname)
                  do id=1,vdims
                    ierr=nf_inq_dim(ncsrc, vdimids(id), dname,
     &                                              count(id))
                    if (ierr == nf_noerr) then
                      ldim=lenstr(dname)
                      start(id)=1
                      if (dname(1:3) == 'xi_') then
                        var_grid_type=var_grid_type+1
                      elseif (dname(1:4) == 'eta_') then
                        var_grid_type=var_grid_type+2
                      elseif (vdimids(id) == unlimdimid) then
                        var_with_recs=.true.
                      elseif (ldim > 3) then
                        if (dname(ldim-3:ldim) == 'time') then
                          var_with_recs=.true.
                        endif
                      endif
                    else
                      write(*,'(/1x,2A,I3/12x,A/)')     '### ERROR: ',
     &                'Cannot get name  and size of dimension with #',
     &                                vdimids(id), nf_strerror(ierr)
                    endif
                  enddo

                  if (ipass == 0) then
                    ierr=nf_def_var(nctarg, vname, vtype, vdims,
     &                                         vdimids, newvarid)
                    if (ierr == nf_noerr) then
                        write(*,'(6x,3A)') 'created variable ''',
     &                                   vname(1:lvar), '''.'
                      do id=1,varatts
                        ierr=nf_inq_attname(ncsrc,varid,id,attname)
                        if (ierr == nf_noerr) then
                          latt=lenstr(attname)
                          ierr=nf_copy_att(ncsrc, varid, attname,
     &                                          nctarg, newvarid)
                          if (ierr == nf_noerr) then
                            write(*,'(8x,4A)') 'copied attribute ''',
     &                                        attname(1:latt), '''.'
                          else
                            write(*,'(/1x,5A/12x,3A/12x,A/)')
     &                         '### ERROR: Cannot copy attribute ''',
     &                          attname(1:latt),''' for variable ''',
     &                          vname(1:lvar), '''',  'in file ''',
     &                          trg_fname(1:ltarg), '''.',
     &                                             nf_strerror(ierr)
                          endif
                        else
                                write(*,'(/1x,2A,I3,1x,3A/12x,A/)')
     &                              '### ERROR: Cannot get name of ',
     &                         'attribute no', id, 'for variable ''',
     &                          vname(1:lvar), '''.',
     &                                             nf_strerror(ierr)
                        endif
                      enddo
                    else
                      write(*,'(/1x,4A/12x,3A/12x,A/)')  '### ERROR: ',
     &                   'Cannot create variable ''',    vname(1:lvar),
     &                   ''' in netCDF file', '''', trg_fname(1:ltarg),
     &                                        '''.', nf_strerror(ierr)
                    endif

                  else !--> ipass > 0

                    ierr=nf_inq_varid(nctarg, vname, newvarid)
                    if (ierr == nf_noerr) then
                      ierr=nf_inq_varndims(nctarg, newvarid, newvdims)
                      if (ierr == nf_noerr) then
                        ierr=nf_inq_vardimid(nctarg, newvarid,
     &                                               newvdimids)
                        if (ierr == nf_noerr) then
                          do id=1,newvdims
                            ierr=nf_inq_dimlen(nctarg, vdimids(id),
     &                                                newcount(id))
                            if (ierr == nf_noerr) then
                              newstart(id)=1
                            else
                              newstart(id)=0
                              write(*,'(/1x,2A,I3,1x,3A/12x,A/)')
     &                        '### ERROR: Cannot determine size of ',
     &                        'dimension #',  id,  'for variable ''',
     &                         vname(1:lvar),'''.',nf_strerror(ierr)
                            endif
                          enddo
                          do id=1,newvdims
                            if (newstart(id) == 0) ierr=1 !<-- ERROR
                          enddo
                        else
                          write(*,'(/1x,4A/12x,A/)')   '### ERROR: ',
     &                            'varndimids for ''', vname(1:lvar),
     &                                      '''?', nf_strerror(ierr)
                        endif
                      else
                        write(*,'(/1x,4A/12x,A/)')     '### ERROR: ',
     &                              'varndims for ''', vname(1:lvar),
     &                                      '''?', nf_strerror(ierr)
                      endif
                    else
                      write(*,'(/1x,4A/12x,A/)') '### ERROR: varid ',
     &                                'for ''', vname(1:lvar), '''?',
     &                                             nf_strerror(ierr)
                    endif

! if everything goes normally so far, both old-variable start:count,
! and new-variable newstart:newcount are determined at this moment.
! Set bounds for loop over record (if any) for this variable
! depending in the status of its "var_with_recs". If the variable is
! to be processed record-by-record, change the counts of the last
! dimension to one.

                    if (var_with_recs) then
                      irec=0                   !<-- input record counter
                      nrecords=count(vdims)    !<-- and its upper bound
                      count(vdims)=1
                      newcount(vdims)=1
                    else
                      irec=-1           !<-- these make while loop
                      nrecords=0        !<-- to be executed only once
                    endif

! Check whether the sizes of allocated arrays are sufficient: if not
! increase them.

                    size=1 ; targ_size=1
                    do id=1,vdims
                      size=size*count(id)
                      targ_size=targ_size*newcount(id)
                    enddo

                    if (vtype == nf_char) then
                     if (size > max_bffrc_size .or.
     &                   targ_size > max_targc_size) then
                        max_bffrc_size=size ; max_targc_size=targ_size
                        if (allocated(targc)) deallocate(targc)
                        if (allocated(bffrc)) deallocate(bffrc)
                        allocate( bffrc(max_bffrc_size),
     &                            targc(max_targc_size) )
                        write(*,'(/1x,2A,8I5)') 'Allocated array ',
     &                   'bffrc: size =', (count(id),id=1,vdims)
                        write(*,'(21x,A,8I5/)') 'targ_size =',
     &                                 (newcount(id),id=1,vdims)
                      endif
                    elseif (vtype == nf_int2) then
                      if (size > max_ibfr2_size .or.
     &                    targ_size > max_itrg2_size) then
                        max_ibfr2_size=size ; max_itrg2_size=targ_size
                        if (allocated(itrg2)) deallocate(itrg2)
                        if (allocated(ibfr2)) deallocate(ibfr2)
                        allocate( ibfr2(max_ibfr2_size),
     &                            itrg2(max_itrg2_size) )
                        write(*,'(/1x,2A,8I5)') 'Allocated array ',
     &                   'ibfr2: size =', (count(id),id=1,vdims)
                        write(*,'(21x,A,8I5/)') 'targ_size =',
     &                                 (newcount(id),id=1,vdims)
                      endif
                    elseif (vtype == nf_int) then
                      if (size > max_ibffr_size .or.
     &                    targ_size > max_itarg_size) then
                        max_ibffr_size=size ; max_itarg_size=targ_size
                        if (allocated(itarg)) deallocate(itarg)
                        if (allocated(ibffr)) deallocate(ibffr)
                        allocate( ibffr(max_ibffr_size),
     &                            itarg(max_itarg_size) )
                        write(*,'(/1x,2A,8I5)') 'Allocated array ',
     &                   'ibffr: size =', (count(id),id=1,vdims)
                        write(*,'(21x,A,8I5/)') 'targ_size =',
     &                                 (newcount(id),id=1,vdims)
                      endif
                    elseif (vtype == nf_real) then
                      if (size > max_bffr4_size .or.
     &                    targ_size > max_targ4_size) then
                        max_bffr4_size=size ; max_targ4_size=targ_size
                        if (allocated(targ4)) deallocate(targ4)
                        if (allocated(bffr4)) deallocate(bffr4)
                        allocate( bffr4(max_bffr4_size),
     &                            targ4(max_targ4_size) )
                        write(*,'(/1x,2A,8I5)') 'Allocated array ',
     &                   'bffr4: size =', (count(id),id=1,vdims)
                        write(*,'(21x,A,8I5/)') 'targ_size =',
     &                                 (newcount(id),id=1,vdims)
                      endif
                    elseif (vtype == nf_double) then
                     if (size > max_bffr8_size .or.
     &                   targ_size > max_targ8_size) then
                        max_bffr8_size=size ; max_targ8_size=targ_size
                        if (allocated(targ8)) deallocate(targ8)
                        if (allocated(bffr8)) deallocate(bffr8)
                        allocate( bffr8(max_bffr8_size),
     &                            targ8(max_targ8_size) )
                        write(*,'(/1x,2A,8I5)') 'Allocated array ',
     &                   'bffr8: size =', (count(id),id=1,vdims)
                        write(*,'(21x,A,8I5/)') 'targ_size =',
     &                                 (newcount(id),id=1,vdims)
                      endif
                    endif

#ifdef VERBOSE
                    write(*,'(1x,3A,2I3,2x,L1)')       'vname = ''',
     &                             vname(1:lvar), ' vtype =', vtype,
     &                                 var_grid_type, var_with_recs
                    do id=1,vdims
                      write(*,'(2I7,1x,A,2I6 )') start(id),count(id),
     &                             '-->', newstart(id), newcount(id)
                    enddo
#endif
                    newrec=last_rec_out !<- record counter for output

                    do while(irec < nrecords .and. ierr == nf_noerr)
                      irec=irec+1         ; newrec=newrec+1
                      if (var_with_recs) then
                        start(vdims)=irec ; newstart(newvdims)=newrec
                      endif
! Reading... -->

                      if (vtype == nf_char) then
                        ierr=nf_get_vara_text (ncsrc, varid, start,
     &                                                count, bffrc)
                      elseif (vtype == nf_int2) then
                        ierr=nf_get_vara_int2 (ncsrc, varid, start,
     &                                                count, ibfr2)
                      elseif (vtype == nf_int) then
                        ierr=nf_get_vara_int (ncsrc, varid, start,
     &                                                count, ibffr)
                      elseif (vtype == nf_real) then
                        ierr=nf_get_vara_real (ncsrc, varid, start,
     &                                                count, bffr4)
                      elseif (vtype == nf_double) then
                        ierr=nf_get_vara_double(ncsrc, varid, start,
     &                                                count, bffr8)
                      else
                        write(*,'(/1x,3A/12x,A/)')    '### ERROR: ',
     &                        '"sample" does not support variable ',
     &                                           'types other than',
     &                          'char, int2, int, real, or double.'
                        ierr=1000
                      endif
                      if (ierr /= nf_noerr) then
                        write(*,'(/1x,6A/12x,A)')      '### ERROR: ',
     &                         'Cannot read ''',     vname(1:lvar),
     &                         ''' from ''', src_file(1:lsrc), '''.',
     &                                             nf_strerror(ierr)
                        write(*,*) 'start =', start
                        write(*,*) 'count =', count
#ifdef VERBOSE
                      elseif (var_with_recs) then
                        write(*,'(5x,A,I5,1x,3A)')  'read rec',
     &                      irec, 'of ''', vname(1:lvar), ''''
                      else
                        write(*,'(1x,4A)')    'read complete ',
     &                      'variable ''', vname(1:lvar), ''''
#endif
                      endif
                      if (ierr == nf_noerr) then

! Writing without extraction...

                        if (var_grid_type == 0) then
                          if (vtype == nf_char) then
                            ierr=nf_put_vara_text (nctarg, newvarid,
     &                                     newstart,newcount, bffrc)
                          elseif (vtype == nf_int2) then
                            ierr=nf_put_vara_int2  (nctarg, newvarid,
     &                                      newstart,newcount, ibfr2)
                          elseif (vtype == nf_int) then
                            ierr=nf_put_vara_int   (nctarg, newvarid,
     &                                      newstart,newcount, ibffr)
                          elseif (vtype == nf_real) then
                            ierr=nf_put_vara_real  (nctarg, newvarid,
     &                                      newstart,newcount, bffr4)
                          elseif (vtype == nf_double) then
                            ierr=nf_put_vara_double(nctarg, newvarid,
     &                                     newstart,newcount, bffr8)
                          endif
                        else

! Extract subdomain and write...

                          if (vtype == nf_char) then
                            call extr_sub_char (imin,imax,jmin,jmax,
     &                                         var_grid_type, vdims,
     &                                  count,newcount, bffrc,targc)
                            ierr=nf_put_vara_text (nctarg, newvarid,
     &                                     newstart,newcount, targc)
                          elseif (vtype == nf_int2) then
                            call extr_sub_int2 (imin,imax,jmin,jmax,
     &                                         var_grid_type, vdims,
     &                                  count,newcount, ibfr2,itrg2)
                            ierr=nf_put_vara_int2 (nctarg, newvarid,
     &                                     newstart,newcount, itrg2)
                          elseif (vtype == nf_int) then
                            call extr_sub_int  (imin,imax,jmin,jmax,
     &                                         var_grid_type, vdims,
     &                                  count,newcount, ibffr,itarg)
                            ierr=nf_put_vara_int   (nctarg, newvarid,
     &                                     newstart,newcount, itarg)

                          elseif (vtype == nf_real) then
                            call extr_sub_real  (imin,imax,jmin,jmax,
     &                                          var_grid_type, vdims,
     &                                  count,newcount, bffr4,targ4)
                            ierr=nf_put_vara_real  (nctarg, newvarid,
     &                                     newstart,newcount, targ4)

                          elseif (vtype == nf_double) then
                            call extr_sub_double(imin,imax,jmin,jmax,
     &                                          var_grid_type, vdims,
     &                                   count,newcount, bffr8,targ8)
                            ierr=nf_put_vara_double(nctarg, newvarid,
     &                                      newstart,newcount, targ8)
                          endif
                        endif
                        if (ierr == nf_noerr) then
                          if (var_with_recs) then
                            write(*,'(4x,A,I5,1x,3A)')   'wrote rec',
     &                          newrec, 'of ''', vname(1:lvar), ''''
                          else
                            write(*,'(1x,4A)')     'wrote complete ',
     &                            'variable ''', vname(1:lvar), ''''
                          endif
                        else
                          write(*,'(/1x,2A/12x,A/)')   '### ERROR: ',
     &                              'Cannot write ''', vname(1:lvar),
     &                            ''' into ''',   trg_fname(1:ltarg),
     &                                      '''.', nf_strerror(ierr)
                        endif
                      endif   !<-- ierr
                    enddo   ! <-- while irec
                  endif   !<-- ipass == 0 .vs. ipass > 0
                else
                  write(*,'(/1x,3A,I3/12x,3A/12x,A/)') '### ERROR: ',
     &              'Cannot make general inquiry for variable with ',
     &              'netCDF ID =', varid, 'in ''',  src_file(1:lsrc),
     &                                     '''.',  nf_strerror(ierr)
                endif
              enddo  !<-- varid
              if (ipass == 0) then
                ierr=nf_enddef(nctarg)
                ierr=nf_sync(nctarg)
                write(*,'(/1x,3A/)') 'Leaving definition mode for ''',
     &                                      trg_fname(1:ltarg), '''.'
              else
                last_rec_out=newrec         ! sync output only when
                ierr=nf_sync(nctarg)        ! going to the next input
                ierr=nf_close(ncsrc)        ! file: doing so too often
                iarg=iarg+1                 ! slows down the execution
                if (iarg < nargs) then
                  call getarg(iarg,  src_file)
                  lsrc=lenstr(src_file)
                  ierr=nf_open(src_file(1:lsrc), nf_nowrite, ncsrc)
                  if (ierr == nf_noerr) then
                    ierr=nf_inq(ncsrc, ndims,nvars,ngatts, unlimdimid)
                    if (ierr == nf_noerr) then
                      write(*,'(/1x,4A,I5,1x,A,I5/)')  'Processing ',
     &                        'file ''',  src_file(1:lsrc),  ''', #',
     &                                   iarg-4,' out of',  nargs-5
                    else
                      write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot',
     &                   ' make general inquiry into netCDF file ''',
     &                    src_file(1:lsrc), '''.', nf_strerror(ierr)
                    endif
                  else
                    write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot ',
     &                       'open netCDF file ''', src_file(1:lsrc),
     &                       ''' for reading.',    nf_strerror(ierr)
                  endif
                endif
              endif
            enddo  !<-- ipass

            ierr=nf_close(nctarg)
            write(*,'(/1x,3A/)')  'Closed netCDF file ''',
     &                          trg_fname(1:ltarg), '''.'
          else
            write(*,'(/1x,4A/12x,A/)')   '### ERROR: Cannot create ',
     &                    'netCDF file ''', trg_fname(1:lsrc), '''.',
     &                                             nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,4A/12x,A/)')       '### ERROR: Cannot make ',
     &       'general inquiry into netCDF file ''', src_file(1:lsrc),
     &                                     '''.',  nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot open netCDF ',
     &                'file ''', src_file(1:lsrc), ''' for reading.',
     &                                             nf_strerror(ierr)
      endif
      end
#endif   /* !VAR_TYPE */

#ifdef VAR_TYPE
# undef VAR_TYPE
# undef extr_sub_TYPE
#endif
#ifndef TYPE
# define extr_sub_TYPE extr_sub_char
# define VAR_TYPE character(len=1)
#elif TYPE == 2
# define extr_sub_TYPE extr_sub_int2
# define VAR_TYPE integer(kind=2)
#elif TYPE == 3
# define extr_sub_TYPE extr_sub_int
# define VAR_TYPE integer(kind=4)
#elif TYPE == 4
# define extr_sub_TYPE extr_sub_real
# define VAR_TYPE real(kind=4)
#elif TYPE == 8
# define extr_sub_TYPE extr_sub_double
# define VAR_TYPE real(kind=8)
#endif

      subroutine extr_sub_TYPE(imin,imax,jmin,jmax, var_grid_type,
     &                          vdims,count, newcount, bffr, targ)
      implicit none
      VAR_TYPE bffr(*), targ(*)
      integer imin,imax,jmin,jmax, var_grid_type, vdims, count(vdims),
     &     newcount(vdims), i,j,k, istr,jstr, itr,isrc, ktr,ksrc,kmax


! Here var_grid_type=1 is for a variable with XI-dimension only;
!                    2 for ETA-only (e.g., in the case of BRY arrays);
!                    3 for having both XI- and ETA-dimensions;

      if ( (var_grid_type == 1 .or. var_grid_type == 3)  .and.
     &    (imin == imax .and. imin > 0 .and. newcount(1) == 2) ) then
        istr=imin-1
      else
        istr=imin
      endif
      if ( (var_grid_type == 2 .or. var_grid_type == 3) .and.
     &    (jmin == jmax .and. jmin > 0 .and. newcount(2) == 2) ) then
        jstr=jmin-1
      else
        jstr=jmin
      endif
      if (var_grid_type == 3) then
        isrc=3
      else
        isrc=2                    ! Fold all dimensions above
      endif                       ! XI and/or ETA into single
      kmax=1                      ! dimension of length "kmax"
      do i=isrc,vdims
        kmax=kmax*count(i)
      enddo

      if (var_grid_type == 3) then
        do k=0,kmax-1
          ktr=k*newcount(1)*newcount(2) ; ksrc=k*count(1)*count(2)
          do j=0,newcount(2)-1
            itr=j*newcount(1)+ktr ; isrc=istr+(j+jstr)*count(1)+ksrc
            do i=1,newcount(1)
              targ(i+itr)=bffr(i+isrc)
            enddo
          enddo
        enddo
      else      !--> vertical slice in either XI- or ETA-directions;
        do k=0,kmax-1
          itr=k*newcount(1)
          if (var_grid_type == 1) then
            isrc=istr +k*count(1)
          elseif (var_grid_type == 2) then
            isrc=jstr +k*count(1)
          endif
          do i=1,newcount(1)
            targ(i+itr)=bffr(i+isrc)
          enddo
        enddo
      endif
      end

#if TYPE < 8
# ifndef TYPE
#  define TYPE 2
# elif TYPE == 2
#  undef TYPE
#  define TYPE 3
# elif TYPE == 3
#  undef TYPE
#  define TYPE 4
# elif TYPE == 4
#  undef TYPE
#  define TYPE 8
# endif
# include "sample.F"
#endif
