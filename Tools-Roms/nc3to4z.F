#define ZLIB_COMPRESSION
#ifdef ZLIB_COMPRESSION
      program nc3to4z
#else
      program nc4to3unz
#endif

! A generic program to convert existing netCDF-3 files into netCDF-4
! with zlib compression (CPP-switch ZLIB_COMPRESSION is activated) or,
! alternatively, convert netCDF-4 into netCDF-3 (ZLIB_COMPRESSION is
! undefined) with inflation of data if compressed.   The latter may be
! useful for compatibility with legacy software tools.  Note that the
! only differences between nc3to4z and nc4to3unz are (i) the argument
! in nf_create() call (netCDF-4  vs. -3); (ii) set deflation for all
! variables which have size of one record (one index of unlimited
! dimension) greater than 16;  and (iii) write label (as a global
! attribute) to mark the file as compressed netCDF-4 format
! (conversely eliminate such attribute in netCDF-4 to 3 conversion).

! Method: for each given file create new file containing the same
! dimensions, variables, and attributes -- everything the same except
! that deflation flag is now set for all arrays larger than a certain
! threshold (very small in practice).  Identify unlimited dimension,
! if exist, and copy data variable-by-variable, record-by-record while
! doing compression of each individual record for each individual
! variable.  Once completed, delete the original file and rename the
! new file into the original.

! CPP-switch ZLIB_COMPRESSION is activated externally and should not
! be defined inside within file.

! CPP-switch EXEC_REPLACE_COMMANDS is for "safety": if activated,
! the original netCDF file will be overwritten. Otherwise kept.

! Relevant man pages about chunking and compression can be found at

!   http://fossies.org/dox/netcdf-4.3.2/netcdf_perf_chunking.html

!   https://www.unidata.ucar.edu/software/netcdf/workshops/ -->
!                               --> 2012/nc4chunking/index.html

!   http://www.hdfgroup.org/pubs/papers/ -->
!                      --> 2008-06_netcdf4_perf_report.pdf

!   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/
!   http://www.iges.org/grads/gadoc/compression.html

#define TIMING
#define EXEC_REPLACE_COMMANDS
c--#define VERBOSE

      implicit none
      character(len=256) src_fname, trg_fname, string
      character(len=64) dimname, attname, varname
      logical set_to_write, exit_no_err
      integer nargs, iarg, ip, getpid, ncsrc, nctarg, ndims,
     &        nvars, ngatts, unlimdimid, varid, vartype, vardims,
     &        varatts, nrecs, rec, recs_in_var, size,
     &        max_cbfr, max_ibfr1, max_ibfr2, max_ibffr, max_bffr4,
     &        max_bffr8, i,j, ierr, lsrc,ltrg, lvar,latt, lstr

      integer lenstr !<-- function

      integer, allocatable, dimension(:) :: dimid, vardimids, dsize,
     &                                       vardsizes, start, count
#ifdef ZLIB_COMPRESSION
     &                                , num_chunks, nblocks, blksize
#endif
      character(len=1), allocatable, dimension(:) :: cbfr
      integer(kind=1),  allocatable, dimension(:) :: ibfr1
      integer(kind=2),  allocatable, dimension(:) :: ibfr2
      integer(kind=4),  allocatable, dimension(:) :: ibffr
      real(kind=4),     allocatable, dimension(:) :: bffr4
      real(kind=8),     allocatable, dimension(:) :: bffr8

      integer :: alloc_dims=0,  alloc_cbfr=0,  alloc_ibfr1=0,
     &           alloc_ibfr2=0, alloc_ibffr=0, alloc_bffr4=0,
     &                                         alloc_bffr8=0

      character(len=9), parameter :: signame='nc_format'
#ifdef ZLIB_COMPRESSION
      integer :: num_dims=0,  alloc_num_chunks=0
      ! DevinD: Sasha had deflate level at 6, but find minimal difference with 1,
      !         and 1 should be faster.
      !         Found 10% or so savings with shuffle on (1).
      integer, parameter :: shuffle=1, deflate=1, deflate_level=1
      integer is,ie, read_str_int
      character(len=25), parameter :: signature =
     &                     'netCDF-4, zlib-compressed'
#endif
#ifdef EXEC_REPLACE_COMMANDS
      character(len=320) command
#endif
#ifdef TIMING
      logical print_time
      real(kind=4) tstart, run_time
      real(kind=8) ReadSize, ReadTime, WrtTime, GrayTime
      integer(kind=4) iclk(2), nclk, clk_rate, clk_max
      integer(kind=8) :: inc_clk, net_read_size=0, net_read_clk=0,
     &                 net_clk=0, net_wrt_clk=0,   net_gray_clk=0
#endif
#include "netcdf.inc"

      nargs=iargc()
      if (nargs == 0) then
        write(*,'(2(/1x,2A)/1x,A/ /8x,2A/1x,A/8x,A/8(/1x,2A)/)')
     &    'Command-line operators "nc3to4z" and "nc4to3unz" perform ',
     &  'netCDF-3 to netCDF-4', 'format conversion with simultaneous ',
     &                    'zlib-compression and the reverse operation',
     &      'respectively. Usage:', 'nc3to4z [--nchunks=m1,m2,m3,..] ',
     &                        'file.nc [file2.nc file3.nc ...]', 'and',
     &                     'nc4to3unz file.nc [file2.nc file3.nc ...]',
     &      'where m1,m2,m3,... are the desired numbers of chunks in ',
     &        'each dimension of the',  'variable (i.e., counting by ',
     &            'array dimension number - not netCDF file dimension',
     &            'IDs).  Setting mN=1 makes chunk size equal to the ',
     &       'whole dimension, while mN=0', 'causes chunk size 1 and ',
     &        'the number of chunks be the same as the dimension.  If',
     &    'the variable has more dimensions than the number of "mN"s ',
     &   'specified, the chunk', 'sizes for the remaining dimensions ',
     &                      'are defaulted to 1.  Unlimited dimension',
     &       'always gets chunk size 1 -- the corresponding "mN" (if ',
     &           'present) is ignored.',  'If argument --nchunks= is ',
     &          'omitted chunking relies on netCDF internal defaults.'

        stop
      endif

#ifdef TIMING
      call cpu_time(tstart)
      nclk=1
      call system_clock(iclk(nclk), clk_rate, clk_max)
# ifdef VERBOSE
      write(*,*) 'clk_rate=', clk_rate, ' clk_max =', clk_max
# endif
#endif


! Create unique name for temporal netCDF "target" file using PID.
! This is useful to avoid collusions if this program needs to be
! run concurrently from a multi-threaded script.

      ip=getpid()
      write(string,'(I10)') ip
      lsrc=lenstr(string)
      trg_fname='tmp'/ /string(1:lsrc)/ /'nc3to4z.nc'
      ltrg=lenstr(trg_fname)
      exit_no_err=.true.
      nrecs=0

! Note that algorithm below allows per-file individual specification
! of chunking parameters, e.g.,
!
!   nc3to4z --nchunks=L,M,N file1.nc --nchunks=mX,mY file2.nc file3.nc
!
! where L,M,N applies to  file1.nc, while mX,mY applies to file2.nc
! and file3.nc.

      iarg=0
   1   iarg=iarg+1
        call getarg(iarg,string) ; lstr=lenstr(string)
#ifdef ZLIB_COMPRESSION
        if (string(1:1) == '-') then
          if (string(1:10) == '--nchunks=') then
            is=0 ; num_dims=1 !<-- at least one; each comma adds one more.
            do i=1,lstr
              if (string(i:i) == '=') then
                is=i+1
              elseif (string(i:i) == ',') then
                num_dims=num_dims+1
              endif
            enddo
            if (is > 0) then
              if (num_dims > alloc_num_chunks) then
                alloc_num_chunks=num_dims
                if (allocated(num_chunks)) deallocate(num_chunks)
                allocate(num_chunks(num_dims))
              endif
              num_chunks=0  !<-- reset
              j=0
              do while(is < lstr+1)
                ie=is
                do while(string(ie:ie) /= ',' .and. ie < lstr)
                  ie=ie+1
                enddo
                if (string(ie:ie) == ',') ie=ie-1
                j=j+1 ; ierr=read_str_int(string(is:ie), num_chunks(j))
                if (ierr /= 0 .or. num_chunks(j) < 0) then
                  write(*,'(/1x,5A/)') '### ERROR: Illegal entry in ''',
     &                          string(1:lstr), ''': ''', string(is:ie),
     &                       ''' is not a non-negative integer number.'
                  stop
                endif
                is=ie+2
              enddo
              write(*,'(/1x,A,8I6)') 'num_chunks =', (num_chunks(j),
     &                                               j=1,num_dims)
              write(*,*)
            else
              write(*,*) '### ERROR: Improper argument ',string(1:lstr)
            endif
          else
            write(*,*) '### ERROR: Unknown parameter input.'; stop
          endif
          if (iarg < nargs) then
             goto 1 !--> proceed to the next iarg
           else
             write(*,*) '### ERROR: No file name is given.'; stop
           endif
        endif    !<-- string(1:1) == '-'
#endif

        src_fname=string(1:lstr) ; lsrc=lenstr(src_fname)
        ierr=nf_open(src_fname(1:lsrc), nf_nowrite, ncsrc)
        if (ierr == nf_noerr) then
#ifdef ZLIB_COMPRESSION
          lvar=len(string)   ! Blank out string before reading it from
          do i=1,lvar        ! netCDF file: the lenght of the attribute
            string(i:i)=' '  ! in may be too short to overwrite it all.
          enddo
          ierr=nf_get_att_text(ncsrc, nf_global, signame, string)
          if (ierr == nf_noerr) then
            lvar=lenstr(string) ; latt=lenstr(signature)
c**         write(*,*) lvar, '''', string(1:lvar), ''''
c**         write(*,*) latt, '''', signature(1:latt), ''''
            if (string(1:lvar) == signature(1:latt)) then
              write(*,'(/1x,3A,1x,A/)')  'File ''', src_fname(1:lsrc),
     &                 ''' is already in compressed netCDF-4 format.',
     &                                                '==> Skipping.'

              ierr=nf_close(ncsrc)
              goto 1 !--> proceed to the next iarg
            endif
          endif
#endif
          write(*,'(/1x,3A)') 'Processing ''',src_fname(1:lsrc),'''...'
          exit_no_err=.true.
          ierr=nf_create(trg_fname(1:ltrg),
#ifdef ZLIB_COMPRESSION
     &                          nf_netcdf4,
#else
     &          nf_clobber+nf_64bit_offset,
#endif
     &                             nctarg)

          if (ierr == nf_noerr) then
            write(*,'(4x,3A)') 'Created netCDF file ''',
     &                             trg_fname(1:ltrg), '''.'
            ierr=nf_inq(ncsrc, ndims, nvars, ngatts, unlimdimid)
            if (ierr == nf_noerr) then
#ifdef VERBOSE
              write(*,'(4(1x,A,I4))') 'ncsrc: ndims =',ndims,'nvars =',
     &            nvars, 'ngatts =',ngatts, 'unlimdimid =',unlimdimid
#endif
              if (ndims > alloc_dims) then
                alloc_dims=ndims
                if (allocated(dimid)) then
                  deallocate(count)      ;  deallocate(start)
#ifdef ZLIB_COMPRESSION
                  deallocate(nblocks)    ;  deallocate(blksize)
#endif
                  deallocate(vardsizes)  ;  deallocate(dsize)
                  deallocate(vardimids)  ;  deallocate(dimid)
                endif
                allocate(dimid(ndims))   ;  allocate(vardimids(ndims))
                allocate(dsize(ndims))   ;  allocate(vardsizes(ndims))
#ifdef ZLIB_COMPRESSION
                allocate(nblocks(ndims)) ;  allocate(blksize(ndims))
                               nblocks=0 ;  blksize=0
#endif
                allocate(start(ndims))   ;  allocate(count(ndims))
                               start=0   ;  count=0
              endif

! Sort out dimensions:

              nrecs=0
              do i=1,ndims
                ierr=nf_inq_dim(ncsrc, i, dimname, dsize(i))
                if (ierr == nf_noerr) then
                  lvar=lenstr(dimname)
                  if (i == unlimdimid) then
                    nrecs=dsize(i) ; dsize(i)=nf_unlimited
                  endif
                  ierr=nf_def_dim(nctarg, dimname(1:lvar), dsize(i),
     &                                                     dimid(i))
                  if (ierr == nf_noerr) then
                    write(*,'(6x,3A,I5)') 'Created dimension ''',
     &                              dimname(1:lvar), '''=', dsize(i)
                  else
                    write(*,'(/1x,6A/12x,A/)')   '### ERROR: Cannot ',
     &                         'define dimension ''', dimname(1:lvar),
     &                      ''' in netCDF file ''', trg_fname(1:ltrg),
     &                                       '''.', nf_strerror(ierr)
                    exit_no_err=.false.
                  endif
                else
                  write(*,'(/1x,2A,I3/12x,3A/12x,A/)')  '### ERROR: ',
     &                'Cannot determine name and size of dimension #',
     &                  i, 'in netCDF file ''',     src_fname(1:lsrc),
     &                                       '''.', nf_strerror(ierr)
                endif
              enddo

! Copy global attributes, however do not copy netCDF format version
! signature, since it is going to be changed in the new file.

              do i=1,ngatts
                ierr=nf_inq_attname(ncsrc, nf_global, i, attname)
                if (ierr == nf_noerr) then
                  latt=lenstr(attname) ; lvar=lenstr(signame)
                  if (attname(1:latt) /= signame(1:lvar)) then
                    ierr=nf_copy_att(ncsrc, nf_global, attname(1:latt),
     &                                              nctarg, nf_global)
                    if (ierr == nf_noerr) then
                      write(*,'(8x,3A)')  'Copied global attribute ''',
     &                                          attname(1:latt), '''.'
                    else
                      write(*,'(/1x,4A/12x,3A/12x,A/)')  '### ERROR: ',
     &               'Cannot copy global attribute ''',attname(1:latt),
     &                '''',  'into netCDF file ''',  trg_fname(1:ltrg),
     &                                        '''.', nf_strerror(ierr)
                      exit_no_err=.false.
                    endif
                  endif
                else
                  write(*,'(/1x,2A,I3,1x,3A/12x,A/)')    '### ERROR: ',
     &                'Cannot determine name of global attribute #', i,
     &            'in ''', src_fname(1:lsrc), '''.', nf_strerror(ierr)
                  exit_no_err=.false.
                endif
              enddo

#ifdef ZLIB_COMPRESSION                             /* Put label so  */
              latt=lenstr(signature)                 ! the file can be
              ierr=nf_put_att_text(nctarg,nf_global, ! later identified
     &                    signame, latt, signature)  ! as compressed
#endif                                              /*  netCDF-4.   */

! Define variables and their attributes.

              max_cbfr=0  ; max_ibfr1=0 ; max_ibfr2=0
              max_ibffr=0 ; max_bffr4=0 ; max_bffr8=0

              do i=1,nvars
                ierr=nf_inq_var(ncsrc, i, varname,  vartype,
     &                          vardims, vardimids, varatts)
                if (ierr == nf_noerr) then
                  lvar=lenstr(varname)
                  ierr=nf_def_var(nctarg, varname(1:lvar),  vartype,
     &                                    vardims, vardimids, varid)
                  if (ierr == nf_noerr) then
                    size=1
                    do j=1,vardims
                      vardsizes(j)=dsize(vardimids(j))
                      if (vardimids(j) /= unlimdimid) then
                        size=size*vardsizes(j)
                      endif
                    enddo
                    if (vartype == nf_char) then
                      if (size > max_cbfr) max_cbfr=size
                    elseif (vartype == nf_byte) then
                      if (size > max_ibfr1) max_ibfr1=size
                    elseif (vartype == nf_short) then
                      if (size > max_ibfr2) max_ibfr2=size
                    elseif (vartype == nf_int) then
                      if (size > max_ibffr) max_ibffr=size
                    elseif (vartype == nf_float) then
                      if (size > max_bffr4) max_bffr4=size
                    elseif (vartype == nf_double) then
                      if (size > max_bffr8) max_bffr8=size
                    endif
                    write(*,'(6x,3A,I10)') 'Created variable ''',
     &                     varname(1:lvar), ''', size =', size
#ifdef ZLIB_COMPRESSION
                    if (size > 16) then
                      do j=1,vardims
                        if (vardimids(j) == unlimdimid) then
                          nblocks(j)=vardsizes(j) ; blksize(j)=1

                        elseif (num_dims > 0) then
                          if (j <= num_dims) then
                            if (num_chunks(j) > 0) then
                              nblocks(j)=num_chunks(j)
                              blksize(j)=( vardsizes(j)+nblocks(j)-1
     &                                                )/nblocks(j)

                            else
                              nblocks(j)=vardsizes(j) ; blksize(j)=1
                            endif
                          else
                            nblocks(j)=vardsizes(j) ; blksize(j)=1
                          endif
# ifdef CUSTOM_DEFAULT
                        else
                          if (j == 1 .or. j == 2) then
                            nblocks(j)=max(1, vardsizes(j)/200)
                            blksize(j)=( vardsizes(j)+nblocks(j)-1
     &                                              )/nblocks(j)
                          else
                            blksize(j)=1 ; nblocks(j)=vardsizes(j)
                          endif
# endif
                       endif
# ifdef VERBOSE
                       write(*,'(1x,A,I3,4(1x,A,I5))') 'dim',j,'size=',
     &                        vardsizes(j), 'chunk size =', blksize(j),
     &                             'nblocks =', nblocks(j), 'excess =',
     &                              blksize(j)*nblocks(j)-vardsizes(j)
# endif
                      enddo
# ifndef CUSTOM_DEFAULT
                      if (num_dims > 0) then
# endif
                        ierr=nf_def_var_chunking(nctarg, varid,
     &                                       nf_chunked, blksize)
                        if (ierr /= nf_noerr) then
                          write(*,'(1x,4A/12x,A)')       '### ERROR: ',
     &                        'Cannot define chunking for variable ''',
     &                       varname(1:lvar), '''.', nf_strerror(ierr)
                        endif
# ifndef CUSTOM_DEFAULT
                      endif
# endif
                      if (ierr == nf_noerr) then
                        ierr=nf_def_var_deflate(nctarg, varid, shuffle,
     &                                          deflate, deflate_level)
                        if (ierr == nf_noerr) then
                          write(*,'(7x,A,3(1x,A,I3))')
     &                     'nf_def_var_deflate:', 'shuffle =', shuffle,
     &                     'deflate =',deflate, 'level =',deflate_level
                        else
                          write(*,'(/1x,4A/12x,A/)')     '### ERROR: ',
     &                          'Cannot set deflation for variable ''',
     &                       varname(1:lvar), '''.' , nf_strerror(ierr)
                        endif
                      endif
                    endif  !<-- size
#endif
                    do j=1,varatts
                      ierr=nf_inq_attname(ncsrc, i, j, attname)
                      latt=lenstr(attname)
                      ierr=nf_copy_att(ncsrc, i,  attname(1:latt),
     &                                             nctarg, varid)
                      if (ierr == nf_noerr) then
                        write(*,'(8x,3A)')    'Copied attribute ''',
     &                                         attname(1:latt), '''.'
                      else
                        write(*,'(/1x,8A/12x,A/)')       '### ERROR: ',
     &                     'Cannot copy attribute ''', attname(1:latt),
     &                           ''' for variable ''', varname(1:lvar),
     &                                 ''' into ''', trg_fname(1:ltrg),
     &                                       '''.',  nf_strerror(ierr)
                      endif
                    enddo
                  else
                    write(*,'(/1x,4A/12x,3A/12x,A/)')    '### ERROR: ',
     &                    'Cannot create variable ''', varname(1:lvar),
     &                    '''', 'in netCDF file ''', trg_fname(1:ltrg),
     &                                        '''.', nf_strerror(ierr)
                  endif
                else
                  write(*,'(/1x,2A,I3,1x,3A/12x,A/)')   '### ERROR: ',
     &               'Cannot make general inquiry for variable # ', i,
     &            'in ''', src_fname(1:lsrc), '''.', nf_strerror(ierr)
                endif
                if (ierr /= nf_noerr) exit_no_err=.false.
              enddo

              ierr=nf_enddef(nctarg)
              write(*,'(6x,3A)')  'Leaving definition mode for ''',
     &                                    trg_fname(1:ltrg), '''.'
            else
               write(*,'(/1x,4A/12x,A/)')    '### ERROR: Cannot make ',
     &        'general inquiry into netCDF file ''', src_fname(1:lsrc),
     &                                        '''.', nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,4A/12x,A/)')     '### ERROR: Cannot create ',
     &     'netCDF file ''',trg_fname(1:ltrg),'''.', nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot open netCDF ',
     &                  'file ''', src_fname(1:lsrc),''' for reading.',
     &                                               nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr) exit_no_err=.false.

        if (max_cbfr > alloc_cbfr) then
          if (allocated(cbfr)) deallocate(cbfr)
          allocate(cbfr(max_cbfr)) ;   alloc_cbfr=max_cbfr
          write(*,*) 'allocated nf_char workspace', max_cbfr
        endif

        if (max_ibfr1 > alloc_ibfr1) then
          if (allocated(ibfr1)) deallocate(ibfr1)
          allocate(ibfr1(max_ibfr1)) ; alloc_ibfr1=max_ibfr1
          write(*,*) 'allocated nf_byte workspace',  max_ibfr1
        endif

        if (max_ibfr2 > alloc_ibfr2) then
          if (allocated(ibfr2)) deallocate(ibfr2)
          allocate(ibfr2(max_ibfr2)) ; alloc_ibfr2=max_ibfr2
          write(*,*) 'allocated nf_int2 workspace',  max_ibfr2
        endif

        if (max_ibffr > alloc_ibffr) then
          if (allocated(ibffr)) deallocate(ibffr)
          allocate(ibffr(max_ibffr)) ; alloc_ibffr=max_ibffr
          write(*,*) 'allocated nf_int workspace',   max_ibffr
        endif

        if (max_bffr4 > alloc_bffr4) then
          if (allocated(bffr4)) deallocate(bffr4)
          allocate(bffr4(max_bffr4)) ; alloc_bffr4=max_bffr4
          write(*,*) 'allocated nf_real workspace',  max_bffr4
        endif

        if (max_bffr8 > alloc_bffr8) then
          if (allocated(bffr8)) deallocate(bffr8)
          allocate(bffr8(max_bffr8)) ; alloc_bffr8=max_bffr8
          write(*,*) 'allocated nf_double workspace', max_bffr8
        endif


!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   *   ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     ***
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***


        do rec=1,max(nrecs,1)
#ifdef TIMING
          print_time=.true.
#endif
          do i=1,nvars
            ierr=nf_inq_var(ncsrc, i, varname,  vartype,
     &                      vardims, vardimids, varatts)
            if (ierr == nf_noerr) then
              if (rec == 1) then            ! During the first record
                set_to_write=.true.         ! write every the variable
              else                          ! regardless of whether it
                set_to_write=.false.        ! has or does not have
              endif                         ! unlimited dimension.
              size=1 ; recs_in_var=1
              do j=1,vardims
                if (vardimids(j) /= unlimdimid) then
                  start(j)=1   ; count(j)=dsize(vardimids(j))
                  size=size*count(j)
                else
                  start(j)=rec ; count(j)=1 ; recs_in_var=nrecs
                  set_to_write=.true.
                endif
              enddo
              if (set_to_write) then
#ifdef TIMING
! Note: there is a possibility that "iclk" goes beyond "clk_max" and
! is folded back into the range of [0 ... clk_max]. The logic needs to
! account for this. Recal that all the arguments of system_clock() are
! 4-byte integers, while all integers computed below system_clock()
! call are 8-byte integers.

                nclk=3-nclk
                call system_clock(iclk(nclk), clk_rate,clk_max)
                inc_clk=iclk(nclk)-iclk(3-nclk)
                if (inc_clk < 0) inc_clk=inc_clk+clk_max
                net_clk=net_clk+inc_clk
                net_gray_clk=net_gray_clk+inc_clk
#endif
                lvar=lenstr(varname)
                if (vartype == nf_char) then
                  ierr=nf_get_vara_text  (ncsrc, i, start,count, cbfr)
                elseif (vartype == nf_byte) then
                  ierr=nf_get_vara_int1  (ncsrc, i, start,count, ibfr1)
                elseif (vartype == nf_short) then
                  ierr=nf_get_vara_int2  (ncsrc, i, start,count, ibfr2)
                elseif (vartype == nf_int) then
                  ierr=nf_get_vara_int   (ncsrc, i, start,count, ibffr)
                elseif (vartype == nf_float) then
                  ierr=nf_get_vara_real  (ncsrc, i, start,count, bffr4)
                elseif (vartype == nf_double) then
                  ierr=nf_get_vara_double(ncsrc, i, start,count, bffr8)
                else
                  write(*,'(/8x,4A/)') '### ERROR: Unrecognized type ',
     &                        'of variable ''', varname(1:lvar), '''.'
                  ierr=nf_noerr-1
                endif
#ifdef TIMING
                nclk=3-nclk
                call system_clock(iclk(nclk), clk_rate, clk_max)
                inc_clk=iclk(nclk)-iclk(3-nclk)
                if (inc_clk < 0) inc_clk=inc_clk+clk_max
                net_clk=net_clk+inc_clk
                net_read_clk=net_read_clk+inc_clk

                if (vartype == nf_char) then
                  net_read_size=net_read_size+size
                elseif (vartype == nf_byte) then
                  net_read_size=net_read_size+size
                elseif (vartype == nf_short) then
                  net_read_size=net_read_size+2*size
                elseif (vartype == nf_int) then
                  net_read_size=net_read_size+4*size
                elseif (vartype == nf_float) then
                  net_read_size=net_read_size+4*size
                elseif (vartype == nf_double) then
                  net_read_size=net_read_size+8*size
                endif
#endif
                if (ierr == nf_noerr) then
                  if (vartype == nf_char) then
                    ierr=nf_put_vara_text(nctarg,i, start,count, cbfr)
                  elseif (vartype == nf_byte) then
                    ierr=nf_put_vara_int1(nctarg,i, start,count, ibfr1)
                  elseif (vartype == nf_short) then
                    ierr=nf_put_vara_int2(nctarg,i, start,count, ibfr2)
                  elseif (vartype == nf_int) then
                    ierr=nf_put_vara_int (nctarg,i, start,count, ibffr)
                  elseif (vartype == nf_float) then
                    ierr=nf_put_vara_real(nctarg,i, start,count, bffr4)
                  elseif (vartype == nf_double) then
                    ierr=nf_put_vara_double(nctarg,
     &                                          i, start,count, bffr8)
                  endif
#ifdef TIMING
                  nclk=3-nclk
                  call system_clock(iclk(nclk), clk_rate, clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  if (inc_clk < 0) inc_clk=inc_clk+clk_max
                  net_clk=net_clk+inc_clk
                  net_wrt_clk=net_wrt_clk+inc_clk
#endif
                  if (ierr == nf_noerr) then
#ifdef TIMING
                    if (print_time .and. clk_rate > 0) then
                      print_time=.false.
                      call system_clock(iclk(nclk), clk_rate, clk_max)
                      write(*,'(F12.2,2x,A,2(I5,1x,A),2A)')
     &                                   dble(net_clk)/dble(clk_rate),
     &                             'wrote rec', rec, '/', recs_in_var,
     &                                 'of ''', varname(1:lvar), ''''
                    else
#endif
                      write(*,'(14x,A,2(I5,1x,A),2A)') 'wrote rec',rec,
     &                '/', recs_in_var, 'of ''', varname(1:lvar), ''''
#ifdef TIMING
                    endif
#endif
                  else
                    write(*,'(/1x,6A/12x,A/)')    '### ERROR: Cannot ',
     &              'write variable ''', varname(1:lvar), ''' into ''',
     &                    trg_fname(1:ltrg),  '''.', nf_strerror(ierr)
                  endif
                else
                  write(*,'(/1x,6A/12x,A/)') '### ERROR: Cannot read ',
     &                    'variable ''', varname(1:lvar), ''' from ''',
     &                     src_fname(1:lsrc), '''.', nf_strerror(ierr)
                endif
              endif  !<-- set_to_write
            else
              write(*,'(/1x,2A,I3/12x,4A/12x,A/)')       '### ERROR: ',
     &                   'Cannot make general inquiry for variable # ',
     &                   i, 'in ', 'netCDF file ''', src_fname(1:lsrc),
     &                                        '''.', nf_strerror(ierr)
            endif
            if (ierr /= nf_noerr) exit_no_err=.false.
          enddo   !<-- i=1,nvars
        enddo   !<-- rec

        ierr=nf_close(nctarg)
        if (ierr == nf_noerr) then
          write(*,'(4x,3A)')  'Closed netCDF file ''',
     &                         trg_fname(1:ltrg), '''.'
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot close netCDF ',
     &         'file ''', trg_fname(1:ltrg), '''.', nf_strerror(ierr)
          exit_no_err=.false.
        endif
        ierr=nf_close(ncsrc)
        if (ierr == nf_noerr) then
          write(*,'(4x,3A)')  'Closed netCDF file ''',
     &                         src_fname(1:lsrc), '''.'
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot close netCDF ',
     &         'file ''', src_fname(1:lsrc), '''.', nf_strerror(ierr)
          exit_no_err=.false.
        endif

#ifdef EXEC_REPLACE_COMMANDS
        if (exit_no_err) then
          write(*,'(/1x,A/)') 'Executing commands:'
          command='/bin/rm -f '/ /src_fname(1:lsrc)
          lvar=lenstr(command)
          write(*,'(8x,3A)')  '''', command(1:lvar), ''''
          call system(command(1:lvar))
          command='/bin/mv -f '/ /trg_fname(1:ltrg)
     &                   / /' '/ /src_fname(1:lsrc)
          lvar=lenstr(command)
          write(*,'(8x,3A/)')  '''', command(1:lvar), ''''
          call system(command(1:lvar))
        endif
#endif
      if (iarg < nargs) goto 1 !--> proceed to the next iarg

#ifdef TIMING
      call cpu_time(run_time) ; run_time=run_time-tstart
      write(*,'(/2x,A,F11.2,1x,A/)') 'cpu time =', run_time, 'sec'

      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate, clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      net_gray_clk=net_gray_clk+inc_clk

      if (clk_rate > 0) then
        ReadSize=1.0D-6*net_read_size
        GrayTime=net_gray_clk/dble(clk_rate)
        ReadTime=net_read_clk/dble(clk_rate)
        WrtTime = net_wrt_clk/dble(clk_rate)

        write(*,'(6x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &         'Total data read :', ReadSize, 'MBytes in',  ReadTime,
     &                           'sec (', ReadSize/ReadTime, 'MB/sec)'
# ifdef ZLIB_COMPRESSION
        write(*,'(1x,A,22x,F12.2,1x,A,F8.2,1x,A)')
     &                                        'Compressing, writing :',
# else
        write(*,'(14x,A,22x,F12.2,1x,A,F8.2,1x,A)')        'Writing :',

# endif
     &                    WrtTime, 'sec (', ReadSize/WrtTime, 'MB/sec)'
        write(*,'(5x,A,22x,F12.2,1x,A)') 'Flush operations :',
     &                                                 GrayTime, 'sec'

        write(*,'(44x,A/10x,A,12x,F12.2,1x,A/)')  '------------------',
     &  'Elapsed wall-clock time', dble(net_clk)/dble(clk_rate), 'sec'
      endif
#endif
      stop
      end
