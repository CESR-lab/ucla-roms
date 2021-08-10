      program partit

! ROMS netCDF file partitioning tool: reads netCDF files containing
! data (usually forcing, initial fields, and boundary conditions) for
! the whole physical grid and prepares multiple files containing only
! data for individual MPI-subdomains. These files are be then read in
! parallel by different MPI processes of running model.

! Usage:  partit NP_XI NP_ETA ncname1 ... ncnameN
! ------  where NP_XI  number of subdomains along XI-direction
!               NP_ETA number of subdomains along ETA-direction
!               ncname1 ... ncnameN  names of netCDF files

! The partitioning algorithm works as follows:  The partitionable
! dimensions ('xi_rho', 'xi_u', 'eta_rho' and 'eta_v') are identified
! by their names, then then their values are read and compared in pairs
! to detect if any of the directions have periodicity.  It is assumed
! that ghost points corresponding to physical boundaries are stored in
! the file, but computational margins (including periodic margins) are
! not.  Consequently, if 'xi_rho' and 'xi_u' are equal to each other
! then XI-direction is periodic, and if they differ by one, it is not.
! ETA-direction is treated similarly.  Once periodicity type is
! determined, the internal number of internal points in each direction
! (i.e. excluding ghost points corresponding to physical boundaries)
! is divided by the number of subdomains in that direction and then
! physical boundary points are attached to subdomains which are
! adjacent to the boundaries. This results in slightly different
! dimension sizes of netCDF files corresponding to different
! subdomains.

! Once all dimensions are sorted out, data corresponding to subdomains
! is extracted from the source file and copied into partial files.

! Non-partitionable objects of netCDF files, such as scalar variables
! and attributes (both global and attributes to variables) are copied
! redundantly into the partitioned files, while partitionable array
! data is subdivided into subdomains and distributed among the
! partitioned files in such a manner that all files contain individual
! data without any overlap or redundantly stored data.

! CPP-switch "FAST_MODE" activates mode where each field (a time record
! there of, if it has record dimension) for the entire physical grid is
! read at once, then subdomain data is extracted from the array and is
! written into partial file for each subdomain.  This is expected to
! speed up reading in comparison with where a rectangular block of data
! corresponding to subdomain is read every time via nf_get_vara_TYPE
! routines.

! To comply with Fortran 2003 standard this code has been revised to
! use individual buffer arrays for each type of data: character, byte,
! int2, integer, single and double precision floating point, so each
! call to nf_get/put_vara_TYPE and extract_TYPE matches exactly the
! internally-declared type of argument, even though it is necessary
! to just match length of each data type in bytes because the buffer
! arrays are just passed from one subroutine to another without
! performing any mathematical operations.

c--#define TRANSPOSED_MPI_NODE_ORDER

#define FAST_MODE
c--#define AUTORENICE
#define TIMING
c--#define VERBOSE


      implicit none
#include "netcdf.inc"
      character(len=128) src_name, string
      integer nargs, nnodes, NP_XI,  xi_rho,  id_xi_rho,  id_xi_psi,
     &        arg,   node,   NP_ETA, eta_rho, id_eta_rho, id_eta_psi,
     &        ncsrc, ndims,  ngatts, xi_u,    id_xi_u,    id_xi_v,
     &        nvars, varatts,        eta_v,   id_eta_v,   id_eta_u,
     &        i,j,k, ierr,   tsize,  rec,     unlimdimid, ierr_all,
     &        lstr,  lsrc,   lfnm,   lvar,    lenstr,     chk_alloc

      character(len=128), allocatable, dimension(:) :: fname
      character(len=32),  allocatable, dimension(:) :: dimname, vname

      integer, allocatable, dimension(:,:) :: vid, dimids
      integer, allocatable, dimension(:) :: ncid,  dimid,  dimsize,
     &                     xi_start, xi_size,   eta_start, eta_size,
     &              vtype, vdims,    part_type, start, count, start1
#ifdef FAST_MODE
     &                                                      , count1
#endif
      logical, allocatable, dimension(:) ::  series, western_edge,
     &                 eastern_edge, southern_edge, northern_edge

      integer :: maxdims=0, maxvars=0 !<-- sizes of allocated arrays


! Buffer arrays for reading and writing, and their allocated sizes

      integer size, size_cbfr,  size_ibfr1, size_ibfr2,
     &              size_ibffr, size_bffr4, size_bffr8

      integer :: alloc_cbfr=0,  alloc_ibfr1=0, alloc_ibfr2=0,
     &           alloc_ibffr=0, alloc_bffr4=0, alloc_bffr8=0

      character(len=1), allocatable, dimension(:) :: cbfr
      integer(kind=1),  allocatable, dimension(:) :: ibfr1
      integer(kind=2),  allocatable, dimension(:) :: ibfr2
      integer(kind=4),  allocatable, dimension(:) :: ibffr
      real(kind=4), allocatable, dimension(:) :: bffr4
      real(kind=8), allocatable, dimension(:) :: bffr8

#ifdef FAST_MODE
      integer size1, size_cbf_in,  size_ibf1_in, size_ibf2_in,
     &               size_ibff_in, size_bff4_in, size_bff8_in

      integer :: allc_cbf_in=0,  allc_ibf1_in=0, allc_ibf2_in=0,
     &           allc_ibff_in=0, allc_bff4_in=0, allc_bff8_in=0

      character(len=1), allocatable, dimension(:) :: cbf_in
      integer(kind=1),  allocatable, dimension(:) :: ibf1_in
      integer(kind=2),  allocatable, dimension(:) :: ibf2_in
      integer(kind=4),  allocatable, dimension(:) :: ibff_in
      real(kind=4), allocatable, dimension(:) :: bff4_in
      real(kind=8), allocatable, dimension(:) :: bff8_in
#endif
#ifdef TIMING
      real(kind=4) tstart, run_time
      real(kind=8) InitTime, ReadSize, ReadTime, WrtSize,  WrtTime,
     &                                           SyncTime, GrayTime
      integer(kind=4) iclk(2), nclk, clk_rate, clk_max, vtypsize
      integer(kind=8) inc_clk
      integer(kind=8) :: net_clk=0, net_read_size=0, net_wrt_size=0,
     &              net_init_clk=0,  net_read_clk=0,  net_wrt_clk=0,
     &                               net_sync_clk=0, net_gray_clk=0
# ifdef FAST_MODE
      real(kind=8)  AssmTime
      integer(kind=8) :: net_assm_clk=0
# endif
#endif              /* Function "iargc" is intrinsic in most modern */
#ifdef INT_IARGC    /* compilers and does not need to be declared.  */
      integer iargc /* IBM xlf95 is a notable exclusion. So do 7.x  */
#endif              /* and earlier versions of Intel compiler, but  */
                    /* it is intrinsic in 9.x and later versions.   */
#ifdef AUTORENICE
      integer getpid, pid            ! Sometimes it makes sense to
      character(len=32) cmd          ! run partit in the background
      pid=getpid()                   ! mode with lowered priority to
      write(cmd,'(I8)') pid          ! not interfere with running MPI
      lstr=lenstr(cmd)               ! job. This code segment catches
      cmd(11:lstr+10)=cmd(1:lstr)    ! pid of its own process and
      lstr=lstr+10                   ! executes re-nice command.
      cmd(1:10)='renice 19 '
      write(*,'(/3A/)') 'Autorenice: executing ''',cmd(1:lstr),'''.'
      call system(cmd(1:lstr))
      write(*,*)
#endif
#ifdef TIMING
      call cpu_time(tstart)
      nclk=1
      call system_clock(iclk(nclk), clk_rate, clk_max)
#endif

! Check how many arguments are given, complain about the error,
! if too few, otherwise extract NP_X and NP_E from the first two
! arguments.

      nargs=iargc()
      if (nargs < 3) then
        write(*,'(/1x,A/30x,A/)')      'Usage of partit should be:',
     &                     'partit NP_X NP_E file1.nc file2.nc ...'
        stop
      endif

      call getarg(1,string) ; lstr=lenstr(string) ; NP_XI=0
      do i=1,lstr
        j=ichar(string(i:i))-48
        if (j < 0 .or. j > 9) then
          write(*,'(/8x,3A/)')     '### ERROR: First argument ',
     &           string(1:lstr), ', must be an integer number.'
          stop
        else
          NP_XI=10*NP_XI+j
        endif
      enddo

      call getarg(2,string) ; lstr=lenstr(string) ; NP_ETA=0
      do i=1,lstr
        j=ichar(string(i:i))-48
        if (j < 0 .or. j > 9) then
          write(*,'(/8x,3A/)')    '### ERROR: Second argument ',
     &           string(1:lstr), ', must be an integer number.'
          stop
        else
          NP_ETA=10*NP_ETA+j
        endif
      enddo
      nnodes=NP_XI*NP_ETA
      write(*,'(/1x,2(4x,A,I3)/)') 'NP_XI =',NP_XI, 'NP_ETA =',NP_ETA

      allocate( ncid(0:nnodes-1),          fname(0:nnodes-1),
     &          xi_size(0:nnodes-1),       eta_size(0:nnodes-1),
     &          xi_start(0:nnodes-1),      eta_start(0:nnodes-1),
     &          western_edge(0:nnodes-1),  eastern_edge(0:nnodes-1),
     &          southern_edge(0:nnodes-1), northern_edge(0:nnodes-1) )

! Process netCDF files: open, check whether it is already a partitioned
! file -- print a warning message and skip the file if it is, otherwise
! proceed with general inquiry, check sizes of book-keeping arrays, and
! if needed allocate them for proper numbers of dimensions and
! variables actually found in netCDF file.

      do arg=3,nargs
        ierr_all=0
        call getarg(arg,src_name)
        lsrc=lenstr(src_name)
        ierr=nf_open(src_name(1:lsrc), nf_nowrite, ncsrc)
        if (ierr == nf_noerr) then
          ierr=nf_inq_att (ncsrc, nf_global, 'partition', i,j)
          if (ierr == nf_noerr) then
            write(*,'(/1x,4A/14x,2A/)')'### WARNING: netCDF file ''',
     &        src_name(1:lsrc),  ''' is already ', 'a partial file ',
     &       'and cannot be partitioned any further ==> ignoring it.'
            goto 97     !--> next file
          endif
        else
          write(*,'(/1x,4A/14x,A)')      '### WARNING: Cannot open ',
     &    'netCDF file ''',src_name(1:lsrc),'''.', nf_strerror(ierr)
          goto 97     !--> next file
        endif

        write(*,'(1x,3A)') 'Processing netCDF file ''',
     &                        src_name(1:lsrc), '''...'

        ierr=nf_inq(ncsrc, ndims, nvars, ngatts, unlimdimid)
        if (ierr == nf_noerr) then
          if (ndims > maxdims .or. nvars > maxvars) then
            if (allocated(dimname)) then
              deallocate( dimname, dimid, dimsize, start, count,
     &                                                       start1,
#ifdef FAST_MODE
     &                                                       count1,
#endif
     &        vname, vtype, vdims, part_type, series, dimids, vid )
            endif
            maxdims=max(ndims,4) ; maxvars=nvars

            allocate( dimname(maxdims), dimid(maxdims),
     &                dimsize(maxdims), start(maxdims), count(maxdims),
     &                                 start1(maxdims)
#ifdef FAST_MODE
     &                                               , count1(maxdims)
#endif
     &                                                               )
            allocate(vname(maxvars),  vtype(maxvars),  vdims(maxvars),
     &                            part_type(maxvars), series(maxvars))

            allocate(dimids(maxdims,maxvars), vid(maxvars,0:nnodes-1))

            write(*,'(1x,2A,I3,1x,A,I3)')   'allocated book-keeping ',
     &           'arrays, maxdims =', maxdims, 'maxvars =', maxvars
          endif
        else
          write(*,'(/1x,4A/12x,A/)')       '### ERROR: Cannot make ',
     &       'general inquiry into netCDF file ''', src_name(1:lsrc),
     &                                     '''.', nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr) stop

! Sort out dimensions: For each dimension find and save its name and
! size. Then check whether all partitionable dimensions (identified
! by names 'xi_rho', 'xi_u', 'eta_rho' and 'eta_v')  are present and
! save their IDs and sizes.

        tsize=-1      ! <-- default value.
        do j=1,ndims
          ierr=nf_inq_dim(ncsrc, j, dimname(j), dimsize(j))
          if (ierr == nf_noerr) then
            if (j == unlimdimid) then
              tsize=dimsize(j)
              dimsize(j)=nf_unlimited
            endif
          else
            write(*,'(/1x,2A,I4/12x,3A/12x,A/)')       '### ERROR: ',
     &           'Cannot determine name and size for dimension #', j,
     &                  'in netCDF file ''', src_name(1:lsrc), '''.',
     &                                             nf_strerror(ierr)
             goto 97     !--> next file
          endif
        enddo

! Determine IDs and sizes of partitionable dimensions, 'xi_rho',
! 'xi_u', 'eta_rho' and 'eta_v'. Also save IDs of obsolete dimensions
! 'xi_psi', 'xi_v', 'eta_psi' and and 'eta_u'. These are used to
! readress obsolete (redundant) dimensions according to the rules:

#ifdef VERBOSE
        write(*,'(4x,A)') 'Identifying dimensions:'
#endif
                                   ! Mapping redundant dimensions:
        xi_rho=0   ; id_xi_rho=0   !
        xi_u=0     ; id_xi_u=0     !         xi_psi  --> xi_u
        eta_rho=0  ; id_eta_rho=0  !
        eta_v=0    ; id_eta_v=0    !         xi_v    --> xi_rho
                                   !
        id_xi_v=0  ; id_xi_psi=0   !         eta_psi --> eta_v
        id_eta_u=0 ; id_eta_psi=0  !
                                   !         eta_u   --> eta_rho
        do j=1,ndims
          lvar=lenstr(dimname(j))
          if (lvar==6 .and. dimname(j)(1:lvar)=='xi_rho') then
            id_xi_rho=j  ; xi_rho=dimsize(j)
          elseif (lvar==4 .and. dimname(j)(1:lvar)=='xi_u') then
            id_xi_u=j    ; xi_u=dimsize(j)
          elseif (lvar==7 .and. dimname(j)(1:lvar)=='eta_rho') then
            id_eta_rho=j ; eta_rho=dimsize(j)
          elseif (lvar==5 .and. dimname(j)(1:lvar)=='eta_v') then
            id_eta_v=j   ; eta_v=dimsize(j)
          elseif (lvar==6 .and. dimname(j)(1:lvar)=='xi_psi') then
            id_xi_psi=j
          elseif (lvar==4 .and. dimname(j)(1:lvar)=='xi_v') then
            id_xi_v=j
          elseif (lvar==7 .and. dimname(j)(1:lvar)=='eta_psi') then
            id_eta_psi=j
          elseif (lvar==5 .and. dimname(j)(1:lvar)=='eta_u') then
            id_eta_u=j
          endif
#ifdef VERBOSE
          write(*,'(8x,I3,1x,A,T24,I3)') j, dimname(j)(1:lvar),
     &                                              dimsize(j)
#endif
        enddo

        if (id_xi_rho /= 0 .and.  id_xi_u == 0) then
          xi_u=xi_rho-1
        elseif (id_xi_rho == 0 .and.  id_xi_u /= 0) then
          xi_rho=xi_u+1
        endif
        if (id_eta_rho /= 0 .and. id_eta_v == 0) then
          eta_v=eta_rho-1
        elseif (id_eta_rho == 0 .and. id_eta_v == 0) then
          eta_rho=eta_v+1
        endif
        if (xi_rho == 0  .or. xi_u == 0  .or.
     &      eta_rho == 0 .or. eta_v == 0) then
          write(*,'(/8x,2A/15x,3A/)')  '### ERROR: not all ',
     &       'partitionable dimensions are found in netCDF ',
     &                    'file ''', src_name(1:lsrc), '''.'
          goto 97     !--> next file
        endif

! Set horizontal dimensions for each subdomain
!---- ---------- ---------- --- ---- ---------

        call mpi_setup (NP_XI,NP_ETA, xi_rho,eta_rho,
     &           xi_start,xi_size, eta_start,eta_size,
     &                    western_edge, eastern_edge,
     &                   southern_edge, northern_edge)


! Create partitioned files:
!======= =========== ======              ! Check whether the source
                                         ! file name contains leading
        j=1                              ! directories, and if it does,
        do i=1,lsrc                      ! take them out:  "partit" is
          if (src_name(i:i)=='/') j=i+1  ! is designed to always place
        enddo                            ! its output files into
        fname(nnodes-1)=src_name(j:lsrc) ! current working directory.

        do node=0,nnodes-1
          fname(node)=fname(nnodes-1) ; lfnm=lenstr(fname(node))
          call insert_node(fname(node), lfnm, node, nnodes, ierr)

          if (ierr == 0) then
            ierr=nf_create(fname(node)(1:lfnm),  nf_clobber +
     &                           nf_64bit_offset, ncid(node))
            if (ierr == nf_noerr) then
              if ( node < 16 .or. ( nnodes > 16 .and.
     &                              node == nnodes-1 )) then
                write(*,'(4x,3A)') 'Created partitioned file ''',
     &                                fname(node)(1:lfnm), '''.'

              elseif (nnodes > 16 .and. node < 18) then
                write(*,'(16x,A)') '.................................'
              endif

              ierr=nf_set_fill(ncid(node), nf_nofill, i)
              if (ierr /= nf_noerr) then
                 write(*,'(/8x,A,1x,3A/8x,A)')   '### ERROR: cannot ',
     &                      'set "nf_nofill" mode for netCDF file ''',
     &                  fname(node)(1:lfnm), '''.', nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,4A/8x,A)' )   '### ERROR: cannot create ',
     &                          'netCDF file ''', fname(node)(1:lfnm),
     &                                      '''.',  nf_strerror(ierr)
            endif
          endif
          if (ierr /=  nf_noerr) stop

! Define dimensions for the partitioned files:
!------- ---------- --- --- ----------- -----

          do j=1,ndims
            if (j == id_xi_rho) then
              size=xi_size(node)
              if (western_edge(node)) size=size+1
              if (eastern_edge(node)) size=size+1

            elseif (j == id_xi_u) then
              size=xi_size(node)
              if (eastern_edge(node)) size=size+1

            elseif (j == id_eta_rho) then
              size=eta_size(node)
              if (southern_edge(node)) size=size+1
              if (northern_edge(node)) size=size+1

            elseif (j == id_eta_v) then
              size=eta_size(node)
              if (northern_edge(node)) size=size+1
            else
              size=dimsize(j)
            endif

            dimid(j)=0 ; lvar=lenstr(dimname(j))

            if (j /= unlimdimid .and. size == 0) then
              if (node == 0) write(*,'(4x,4A)')   'Suppressing ',
     &       'zero-size dimension ''', dimname(j)(1:lvar), '''.'
            elseif  (j == unlimdimid .and. tsize == 0) then
              if (node == 0) write(*,'(4x,4A)')   'Suppressing ',
     &                        'zero-size unlimited dimension ''',
     &                                 dimname(j)(1:lvar), '''.'
            elseif (j == id_xi_psi) then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                          'obsolete dimension ''xi_psi''.'
            elseif (j == id_xi_v) then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                            'obsolete dimension ''xi_v''.'
            elseif (j == id_eta_psi) then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                         'obsolete dimension ''eta_psi''.'
            elseif (j == id_eta_u)  then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                            'obsolete dimension ''xi_u''.'
            else
              ierr=nf_def_dim(ncid(node), dimname(j)(1:lvar),
     &                                        size, dimid(j))
              if (ierr /= nf_noerr) then
                ierr_all=ierr_all+1
                write(*,'(/1x,4A,I8,A,I4/12x,A)')    '### ERROR: ',
     &            'Cannot define dimension ''', dimname(j)(1:lvar),
     &                       ''' of size =', size,',  node=', node,
     &                                           nf_strerror(ierr)
#ifdef VERBOSE
              else
                write(*,'(8x,2A,I2,2x,2A,I4,T50,A,I4)') 'Defined ',
     &                  'dimension #',dimid(i), dimname(j)(1:lvar),
     &                                 ' =', size, 'node =', node
#endif
              endif
            endif
          enddo  !<-- loop over dimensions

! After this moment array dimid(1:ndims) contains the set of NEW
! dimension IDs.  If the original file contains any of the four
! obsolete dimensions, 'xi_psi', 'eta_psi', 'xi_v', and 'eta_u'
! which have been eliminated, then dimid(i) does not correspond to
! the set of dimension IDs of the original file [which would be
! just dimid(i)=i]. Consequently, array dimid(1:ndims) will later
! be used later to remap old dimension IDs into new ones.
!
! Put global attributes:    The new attribute 'partition' identifies
!---- ------ -----------    identifies position of each individual
!                           subdomain within the processor grid

          start(1)=node    ; start(3)=xi_start(node)
          start(2)=nnodes  ; start(4)=eta_start(node)
          ierr=nf_put_att_int(ncid(node), nf_global, 'partition',
     &                                       nf_int,   4, start)
        enddo  !<-- node=0,nnodes-1
        if (ierr_all /= 0) stop

! Copy global attributes

        do i=1,ngatts
          ierr=nf_inq_attname (ncsrc, nf_global, i, string)
          if (ierr == nf_noerr) then
            lvar=lenstr(string)
            do node=0,nnodes-1
              ierr=nf_copy_att (ncsrc, nf_global, string(1:lvar),
     &                                     ncid(node), nf_global)
              if (ierr. ne. nf_noerr) then
                ierr_all=ierr_all+1
                write(*,'(/1x,7A/12x,A)')  '### ERROR: Cannot copy ',
     &             'global attribute ''', string(1:lvar), ''' into ',
     &             'netCDF file ''',   fname(node)(1:lfnm),   '''.',
     &                                             nf_strerror(ierr)
                goto 97
              endif
            enddo
          else
            write(*,'(/1x,2A,I3,1x,4A/12x,A/)') '### ERROR: Cannot ',
     &       'determine name of global attribute #', i, 'in netCDF ',
     &       'file ''', src_name(1:lsrc), '''.',  nf_strerror(ierr)
            goto 97
          endif
        enddo

! Define variables and their attributes in the partitioned files.
!------- --------- --- ----- ---------- -- --- ----------- ------

        do i=1,nvars
          ierr=nf_inq_var(ncsrc,  i,    vname(i),  vtype(i),
     &                    vdims(i), dimids(1,i),   varatts)
          if (ierr /= nf_noerr) then
             write(*,'(/1x,2A,I3/12x,3A/12x,A/)')     '### ERROR: ',
     &            'Cannot make general inquiry about variable ID =',
     &             i, 'in netCDF file ''',  src_name(1:lsrc), ''',',
     &                                            nf_strerror(ierr)

            goto 97
          endif

! Readress obsolete dimensions, if any:

          do j=1,vdims(i)
            if (dimids(j,i) == id_xi_psi) then
              dimids(j,i)=id_xi_u
            elseif (dimids(j,i) == id_xi_v) then
              dimids(j,i)=id_xi_rho
            elseif (dimids(j,i) == id_eta_psi) then
              dimids(j,i)=id_eta_v
            elseif (dimids(j,i) == id_eta_u) then
              dimids(j,i)=id_eta_rho
            endif
          enddo

! Determine whether partitionable dimensions or unlimited dimension
! are present for this variable: the convention adopted here is:
!         part_type = 0 -- non-partitionable array;
!                   = 1 -- has partitionable XI-dimension only;
!                   = 2 -- has partitionable ETA-dimension only;
!                   = 3 -- partitionable in both XI and ETA.

          series(i)=.false.
          part_type(i)=0
          do j=1,vdims(i)
            if ( dimids(j,i) == id_xi_rho .or.
     &           dimids(j,i) == id_xi_u  ) then
              part_type(i)=part_type(i)+1
            elseif ( dimids(j,i) == id_eta_rho .or.
     &               dimids(j,i) == id_eta_v ) then
              part_type(i)=part_type(i)+2
            elseif (dimids(j,i)==unlimdimid) then
              series(i)=.true.
            endif
          enddo

          if (tsize > 0 .or. .not.series(i)) then

! WARNING: Since dimids(1:vdims(i),i) contains dimension IDs
! corresponding to the set of IDs of the ORIGINAL file, and since
! some of the original dimensions were eliminated (merged), the set
! of dimension IDs in the NEW definitions is obtained by the inverse
! mapping of dimids(j,i) onto start(j) using dimid(k) as a mapping
! array.

            do j=1,vdims(i)
              do k=1,ndims
                if (dimids(j,i) == k) start(j)=dimid(k)
              enddo
            enddo
c**         write(*,*) 'old_dimids:', (dimids(j,i),j=1,vdims(i))
c**         write(*,*) 'new_dimids:',    (start(j),j=1,vdims(i))

            lvar=lenstr(vname(i))
            do node=0,nnodes-1
              vid(i,node)=0
              if (lvar > 5) then
                if (vname(i)(lvar-4:lvar) == '_west' .and.
     &              .not.western_edge(node)) vid(i,node)=-1000

                if (vname(i)(lvar-4:lvar) == '_east'  .and.
     &              .not.eastern_edge(node)) vid(i,node)=-1000
              endif
              if (lvar > 6) then
                if (vname(i)(lvar-5:lvar) == '_south' .and.
     &              .not.southern_edge(node)) vid(i,node)=-1000

                if (vname(i)(lvar-5:lvar) == '_north'  .and.
     &              .not.northern_edge(node)) vid(i,node)=-1000
              endif
              if (vid(i,node) == 0) then
                ierr=nf_def_var(ncid(node), vname(i)(1:lvar), vtype(i),
     &                                   vdims(i), start, vid(i,node))
                if (ierr /= nf_noerr) then
                  write(*,'(/1x,4A,I4/12x,A/)') '### ERROR: Cannot ',
     &                    'create variable ''',   vname(i)(1:lvar),
     &                  ''', node =',   node,    nf_strerror(ierr)
#ifdef VERBOSE
                else
                 write(*,'(I3,3A,T36,2(A,I3,1x),2(1x,A,I2),1x,A,I4)')
     &             i, ' Created variable ''',  vname(i)(1:lvar), '''',
     &                'varid =',  vid(i,node),   'vdims =',  vdims(i),
     &               'part_type =',part_type(i), 'series =',series(i),
     &                                                'node =', node
#endif
                endif

#ifdef DUMMY_BRY_VAR
              else
                ierr=nf_def_var(ncid(node), vname(i)(1:lvar), vtype(i),
     &                                          0, start, vid(i,node))
# ifdef VERBOSE
                write(*,'(I3,4A,T84,A,I4)') i,' Created dummy scalar ',
     &           'variable ''', vname(i)(1:lvar), '''', 'node =', node
# endif
#else
# ifdef VERBOSE
              else
                write(*,'(I3,4A,T84,A,I4)')          i, ' Suppressed ',
     &           'variable ''', vname(i)(1:lvar), '''', 'node =', node
# endif
#endif
              endif
            enddo     !--> discard start
            do j=1,varatts
              ierr=nf_inq_attname(ncsrc, i, j, string)
              lvar=lenstr(string)
              do node=0,nnodes-1
                if (vid(i,node)>0) then
                  ierr=nf_copy_att(ncsrc, i, string(1:lvar),
     &                              ncid(node), vid(i,node))
                  if (ierr /= nf_noerr) then
                    write(*,'(/1x,4A,I4/12x,A/)') '### ERROR: ',
     &               'Cannot copy attribute ''', string(1:lvar),
     &               ''' for variable ''',   vname(i)(1:lvar),
     &               ''', node =',   node,    nf_strerror(ierr)

                  endif
                endif
              enddo
            enddo
          endif    !<-- tsize > 0 .or. .not.series(i)
        enddo    ! <-- i, loop over variables

! Leave definition mode:
!------ ---------- -----

        do node=0,nnodes-1
          ierr=nf_enddef(ncid(node))
          if (ierr /= nf_noerr) then
            write(*,'(/1x,5A/)')  '### ERROR: Cannot switch partial ',
     &                          'netCDF file ''', fname(node)(1:lfnm),
     &                      ''' into input mode, ', nf_strerror(ierr)
            ierr_all=ierr_all+1
          endif
        enddo
        if (ierr_all /= 0) stop
#ifdef VERBOSE
        write(*,'(/4x,A/)') 'Left definition mode.'
#endif


! Allocate buffer arrays: Scan all the variables to find the largest
!--------- ------ ------- memory requirement within each variable type.

  1     size_cbfr=0     ; size_ibfr1=0    ; size_ibfr2=0
        size_ibffr=0    ; size_bffr4=0    ; size_bffr8=0
#ifdef FAST_MODE
        size_cbf_in=0   ; size_ibf1_in=0  ; size_ibf2_in=0
        size_ibff_in=0  ; size_bff4_in=0  ; size_bff8_in=0
#endif

        do i=1,nvars
          do node=0,nnodes-1
            if (vid(i,node) > 0) then
              size=1
              do j=1,vdims(i)
                if ( dimids(j,i) == id_xi_rho .or.
     &               dimids(j,i) == id_xi_u ) then
                  count(j)=xi_size(node)
                  if ( western_edge(node) .and.
     &                 dimids(j,i) == id_xi_rho ) then
                    count(j)=count(j)+1
                  endif
                  if (eastern_edge(node)) then
                    count(j)=count(j)+1
                  endif
                elseif ( dimids(j,i) == id_eta_rho .or.
     &                   dimids(j,i) == id_eta_v ) then
                  count(j)=eta_size(node)
                  if ( southern_edge(node) .and.
     &                 dimids(j,i) == id_eta_rho) then
                    count(j)=count(j)+1
                  endif
                  if (northern_edge(node)) then
                    count(j)=count(j)+1
                  endif
                elseif (dimids(j,i) == unlimdimid) then
                  count(j)=1
                else
                  count(j)=dimsize(dimids(j,i))
                endif
                size=size*count(j)
              enddo

              if (vtype(i) == nf_char) then
                size_cbfr=max(size_cbfr,size)
              elseif (vtype(i) == nf_byte) then
                size_ibfr1=max(size_ibfr1,size)
              elseif (vtype(i) == nf_short) then
                size_ibfr2=max(size_ibfr2,size)
              elseif (vtype(i) == nf_int) then
                size_ibffr=max(size_ibffr,size)
              elseif (vtype(i) == nf_float) then
                size_bffr4=max(size_bffr4,size)
              elseif (vtype(i) == nf_double) then
                size_bffr8=max(size_bffr8,size)
              else
                write(*,'(/8x,3A/)')  '### ERROR: Variable ''',
     &              vname(i)(1:lvar), ''' is of unknown type.'
                stop
             endif
            endif
          enddo  !<--  node=0,nnodes-1

#ifdef FAST_MODE
          if (part_type(i) > 0) then
            size1=1
            do j=1,vdims(i)
              if (dimids(j,i) /= unlimdimid) then
                size1=size1*dimsize(dimids(j,i))
              endif
            enddo

            if (vtype(i) == nf_char) then
              size_cbf_in=max(size_cbf_in,size1)
            elseif (vtype(i) == nf_byte) then
              size_ibf1_in=max(size_ibf1_in,size1)
            elseif (vtype(i) == nf_short) then
              size_ibf2_in=max(size_ibf2_in,size1)
            elseif (vtype(i) == nf_int) then
              size_ibff_in=max(size_ibff_in,size1)
            elseif (vtype(i) == nf_float) then
              size_bff4_in=max(size_bff4_in,size1)
            elseif (vtype(i) == nf_double) then
              size_bff8_in=max(size_bff8_in,size1)
            endif
          endif
#endif
        enddo   !<-- i=1,nvars

        if (size_cbfr > alloc_cbfr) then
          if (allocated(cbfr)) deallocate(cbfr)
          allocate(cbfr(size_cbfr)) ;  alloc_cbfr=size_cbfr
          write(*,*) 'allocated nf_char workspace', size_cbfr
        endif
        if (size_ibfr1 > alloc_ibfr1) then
          if (allocated(ibfr1)) deallocate(ibfr1)
          allocate(ibfr1(size_ibfr1)) ; alloc_ibfr1=size_ibfr1
          write(*,*) 'allocated nf_byte workspace', size_ibfr1
        endif
        if (size_ibfr2 > alloc_ibfr2) then
          if (allocated(ibfr2)) deallocate(ibfr2)
          allocate(ibfr2(size_ibfr2)) ; alloc_ibfr2=size_ibfr2
          write(*,*) 'allocated nf_int2 workspace', size_ibfr2
        endif
        if (size_ibffr > alloc_ibffr) then
          if (allocated(ibffr)) deallocate(ibffr)
          allocate(ibffr(size_ibffr)) ; alloc_ibffr=size_ibffr
          write(*,*)  'allocated nf_int workspace', size_ibffr
        endif
        if (size_bffr4 > alloc_bffr4) then
          if (allocated(bffr4)) deallocate(bffr4)
          allocate(bffr4(size_bffr4)) ; alloc_bffr4=size_bffr4
          write(*,*) 'allocated nf_real workspace', size_bffr4
        endif
        if (size_bffr8 > alloc_bffr8) then
          if (allocated(bffr8)) deallocate(bffr8)
          allocate(bffr8(size_bffr8)) ; alloc_bffr8=size_bffr8
          write(*,*)'allocated nf_double workspace',size_bffr8
        endif

#ifdef FAST_MODE
        if (size_cbf_in > allc_cbf_in) then
          if (allocated(cbf_in)) deallocate(cbf_in)
          allocate(cbf_in(size_cbf_in)) ; allc_cbf_in=size_cbf_in
          write(*,*) 'allocated nf_char workspace', size_cbf_in
        endif
        if (size_ibf1_in > allc_ibf1_in) then
          if (allocated(ibf1_in)) deallocate(ibf1_in)
          allocate(ibf1_in(size_ibf1_in)) ; allc_ibf1_in=size_ibf1_in
          write(*,*) 'allocated nf_byte workspace',size_ibf1_in
        endif
        if (size_ibf2_in > allc_ibf2_in) then
          if (allocated(ibf2_in)) deallocate(ibf2_in)
          allocate(ibf2_in(size_ibf2_in)) ; allc_ibf2_in=size_ibf2_in
          write(*,*) 'allocated nf_int2 workspace', size_ibf2_in
        endif
        if (size_ibff_in > allc_ibff_in) then
          if (allocated(ibff_in)) deallocate(ibff_in)
          allocate(ibff_in(size_ibff_in)) ; allc_ibff_in=size_ibff_in
          write(*,*) 'allocated nf_int workspace',size_ibff_in
        endif
        if (size_bff4_in > allc_bff4_in) then
          if (allocated(bff4_in)) deallocate(bff4_in)
          allocate(bff4_in(size_bff4_in)) ; allc_bff4_in=size_bff4_in
          write(*,*) 'allocated nf_real workspace',size_bff4_in
        endif
        if (size_bff8_in > allc_bff8_in) then
          if (allocated(bff8_in)) deallocate(bff8_in)
          allocate(bff8_in(size_bff8_in)) ; allc_bff8_in=size_bff8_in
            write(*,*) 'allocated nf_double workspace', size_bff8_in
        endif
#endif

#ifdef TIMING
! Note: there is a possibility that "iclk" goes beyond "clk_max" and
! is folded back into the range of [0 ... clk_max]. The logic needs to
! account for this. Recall that all the arguments of system_clock() are
! 4-byte integers, while all integers computed below system_clock()
! call are 8-byte integers.

        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_init_clk=net_init_clk+inc_clk
#endif


!     **     *     ***  *******    ***  *********  ********
!      *    ***   ***   ***   ***  ***  *  ***  *  ***    *
!       *   ***   ***   ***   ***  ***     ***     ***
!       *  *** * ***    ***   **   ***     ***     ******
!        * **  * **     ******     ***     ***     ***
!        ***    ***     ***  **    ***     ***     ***    *
!         *     **      ***   ***  ***     ***     ********



! Transfer variables into partial files.
!========= ========= ==== ======= =====

        do rec=1,max(tsize,1)
          if (tsize > 1) then
#ifdef TIMING
            nclk=3-nclk
            call system_clock(iclk(nclk), clk_rate,clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            if (inc_clk < 0) inc_clk=inc_clk+clk_max
            net_clk=net_clk+inc_clk
            net_gray_clk=net_gray_clk+inc_clk
            write(*,'(F8.1,1x,2(A,I8,1x),A)')
     &                         dble(net_clk)/dble(clk_rate),
#else
            write(*,'(4x,2(A,I8,1x),A)')
#endif
     &       'rec', rec, 'out of', tsize,  '...'
          endif

          do i=1,nvars
            lvar=lenstr(vname(i))
            if (rec==1 .or. series(i)) then
#ifdef TIMING
              if (vtype(i) == nf_char .or.        ! Determine the size
     &            vtype(i) == nf_byte) then       ! of one element of
                vtypsize=1                        ! the variable in
              elseif (vtype(i) == nf_short) then  ! bytes (needed only
                vtypsize=2                        ! to report the total
              elseif (vtype(i) == nf_int .or.     ! amount of data in
     &                vtype(i) == nf_float) then  ! MBytes).
                vtypsize=4
              elseif (vtype(i) == nf_double) then
                vtypsize=8
              endif
#endif
              if (part_type(i) == 0 .and. .not.series(i)) then

! Scalar variable (actually it may be an array, but having neither
! partitionable spatial dimention, nor record dimension, i.e., it is
! treated as a whole and only during rec=1):

                write(*,'(21x,A,I4,1x,A,I4,1x,3A)') 'var', i, 'out of',
     &           nvars,'non-part.  non-rec. ''',vname(i)(1:lvar),'''.'

                size=1
                do j=1,vdims(i)
                  size=size * dimsize(dimids(j,i))
                enddo

                if (vtype(i) == nf_char) then
                  ierr=chk_alloc(size, alloc_cbfr, 'nf_char')
                elseif (vtype(i) == nf_byte) then
                  ierr=chk_alloc(size, alloc_ibfr1, 'nf_byte')
                elseif (vtype(i) == nf_short) then
                  ierr=chk_alloc(size, alloc_ibfr2, 'nf_short')
                elseif (vtype(i) == nf_int) then
                  ierr=chk_alloc(size, alloc_ibffr,   'nf_int')
                elseif (vtype(i) == nf_float) then
                  ierr=chk_alloc(size, alloc_bffr4, 'nf_float')
                elseif (vtype(i) == nf_double) then
                  ierr=chk_alloc(size,alloc_bffr8, 'nf_double')
                else
                  write(*,'(/1x,4A/)')  '### ERROR: Unknown type ',
     &                 'of variable ''', vname(i)(1:lvar), '''.'
                  stop
                endif
                if (ierr /= 0) stop

                if (vtype(i) == nf_char) then
                  ierr=nf_get_var_text (ncsrc, i, cbfr)
                elseif (vtype(i) == nf_byte) then
                  ierr=nf_get_var_int1   (ncsrc, i, ibfr1)
                elseif (vtype(i) == nf_short) then
                  ierr=nf_get_var_int2   (ncsrc, i, ibfr2)
                elseif (vtype(i) == nf_int) then
                  ierr=nf_get_var_int    (ncsrc, i, ibffr)
                elseif (vtype(i) == nf_float) then
                  ierr=nf_get_var_real   (ncsrc, i, bffr4)
                elseif (vtype(i) == nf_double) then
                  ierr=nf_get_var_double (ncsrc, i, bffr8)
                endif

                if (ierr /= nf_noerr) then
                  write(*,'(/1x,6A/12x,A/)')      '### ERROR: ',
     &              'Cannot read variable ''', vname(i)(1:lvar),
     &              ''' from netCDF file ''',  src_name(1:lsrc),
     &                                 '''.', nf_strerror(ierr)
                  goto 97
#ifdef TIMING
                else
                  nclk=3-nclk
                  call system_clock(iclk(nclk),clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  if (inc_clk < 0) inc_clk=inc_clk+clk_max
                  net_clk=net_clk+inc_clk
                  net_read_clk=net_read_clk+inc_clk
                  net_read_size=net_read_size+size*vtypsize
#endif
                endif

                do node=0,nnodes-1
                  if (vid(i,node) > 0) then
                    if (vtype(i) == nf_char) then
                      ierr=nf_put_var_text (ncid(node), vid(i,node),
     &                                                        cbfr)
                    elseif (vtype(i) == nf_byte) then
                      ierr=nf_put_var_int1 (ncid(node), vid(i,node),
     &                                                       ibfr1)
                    elseif (vtype(i) == nf_short) then
                      ierr=nf_put_var_int2 (ncid(node), vid(i,node),
     &                                                       ibfr2)
                    elseif (vtype(i) == nf_int) then
                      ierr=nf_put_var_int  (ncid(node), vid(i,node),
     &                                                       ibffr)
                    elseif (vtype(i) == nf_float) then
                      ierr=nf_put_var_real (ncid(node), vid(i,node),
     &                                                       bffr4)
                    elseif (vtype(i) == nf_double) then
                      ierr=nf_put_var_double(ncid(node),vid(i,node),
     &                                                       bffr8)
                    endif

                    if (ierr /= nf_noerr) then
                      write(*,'(/1x,6A/12x,A/)')       '### ERROR:',
     &                'Cannot write variable ''',  vname(i)(1:lvar),
     &                ''' into netCDF file ''', fname(node)(1:lfnm),
     &                                     '''.', nf_strerror(ierr)
                      goto 97
#ifdef TIMING
                    else
                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      if (inc_clk < 0) inc_clk=inc_clk+clk_max
                      net_clk=net_clk+inc_clk
                      net_wrt_clk=net_wrt_clk+inc_clk
                      net_wrt_size=net_wrt_size+size*vtypsize
#endif
                    endif
                  endif
                enddo
              elseif (part_type(i) == 0) then

! Non-partitionable array (but with records):

#ifndef VERBOSE
                if (rec < 5) then
#endif
                  write(*,'(21x,A,I4,1x,A,I4,1x,3A)') 'var',i,'out of',
     &                 nvars, 'non-part.  ''', vname(i)(1:lvar), '''.'
#ifndef VERBOSE
                elseif (rec < 6) then
                   write(*,'(24x,A)') '.........................'
                endif
#endif
                size=1
                do j=1,vdims(i)
                  if (dimids(j,i)==unlimdimid) then
                    start(j)=rec ; count(j)=1
                  else
                    start(j)=1 ; count(j)=dimsize(dimids(j,i))
                  endif
                  size=size*count(j)
                enddo

                if (vtype(i) == nf_char) then
                  ierr=chk_alloc(size, alloc_cbfr,  'nf_char')
                elseif (vtype(i) == nf_byte) then
                  ierr=chk_alloc(size, alloc_ibfr1,  'nf_byte')
                elseif (vtype(i) == nf_short) then
                  ierr=chk_alloc(size, alloc_ibfr2,  'nf_short')
                elseif (vtype(i) == nf_int) then
                  ierr=chk_alloc(size, alloc_ibffr,    'nf_int')
                elseif (vtype(i) == nf_float) then
                  ierr=chk_alloc(size, alloc_bffr4,  'nf_float')
                elseif (vtype(i) == nf_double) then
                  ierr=chk_alloc(size, alloc_bffr8, 'nf_double')
                else
                  write(*,'(/1x,4A/)')   '### ERROR: Unknown ',
     &            'type of variable ''',vname(i)(1:lvar),'''.'
                  stop
                endif
                if (ierr /= 0) stop

                if (vtype(i) == nf_char) then
                  ierr=nf_get_vara_text (ncsrc, i, start,count, cbfr)
                elseif (vtype(i) == nf_byte) then
                  ierr=nf_get_vara_int1  (ncsrc, i, start,count, ibfr1)
                elseif (vtype(i) == nf_short) then
                  ierr=nf_get_vara_int2  (ncsrc, i, start,count, ibfr2)
                elseif (vtype(i) == nf_int) then
                  ierr=nf_get_vara_int   (ncsrc, i, start,count, ibffr)
                elseif (vtype(i) == nf_float) then
                  ierr=nf_get_vara_real  (ncsrc, i, start,count, bffr4)
                elseif (vtype(i) == nf_double) then
                  ierr=nf_get_vara_double(ncsrc, i, start,count, bffr8)
                endif

                if (ierr /= nf_noerr) then
                  write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A/)')
     &              '### ERROR: Cannot read time record =',    rec,
     &              'for nonpartitionable array ''', vname(i)(1:lvar),
     &              ''' from netCDF',   'file ''',   src_name(1:lsrc),
     &                                       '''.', nf_strerror(ierr)
                  goto 97
#ifdef TIMING
                else
                  nclk=3-nclk
                  call system_clock(iclk(nclk),clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  if (inc_clk < 0) inc_clk=inc_clk+clk_max
                  net_clk=net_clk+inc_clk
                  net_read_clk=net_read_clk+inc_clk
                  net_read_size=net_read_size+size*vtypsize
#endif
                endif

                do node=0,nnodes-1
                  if (vid(i,node) > 0) then
                    if (vtype(i) == nf_char) then
                      ierr=nf_put_vara_text (ncid(node), vid(i,node),
     &                                             start,count, cbfr)
                    elseif (vtype(i) == nf_byte) then
                      ierr=nf_put_vara_int1  (ncid(node), vid(i,node),
     &                                            start,count, ibfr1)
                    elseif (vtype(i) == nf_short) then
                      ierr=nf_put_vara_int2  (ncid(node), vid(i,node),
     &                                            start,count, ibfr2)
                    elseif (vtype(i) == nf_int) then
                      ierr=nf_put_vara_int   (ncid(node), vid(i,node),
     &                                            start,count, ibffr)
                    elseif (vtype(i) == nf_float) then
                      ierr=nf_put_vara_real  (ncid(node), vid(i,node),
     &                                            start,count, bffr4)
                    elseif (vtype(i) == nf_double) then
                      ierr=nf_put_vara_double(ncid(node), vid(i,node),
     &                                            start,count, bffr8)
                  endif

                    if (ierr /= nf_noerr) then
                      write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A/)')
     &               '### ERROR: Cannot write time record =',     rec,
     &               'for nonpartitionable array ''',vname(i)(1:lvar),
     &               '''into netCDF',  'file ''', fname(node)(1:lfnm),
     &                                       '''.', nf_strerror(ierr)
                      goto 97
#ifdef TIMING
                    else
                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      if (inc_clk < 0) inc_clk=inc_clk+clk_max
                      net_clk=net_clk+inc_clk
                      net_wrt_clk=net_wrt_clk+inc_clk
                      net_wrt_size=net_wrt_size+size*vtypsize
#endif
                    endif
                  endif
                enddo
              elseif (part_type(i) > 0) then

! Partitioned array:

#ifndef VERBOSE
                if (rec < 5) then
#endif
                  write(*,'(21x,A,I4,1x,A,I4,1x,3A)') 'var',i,'out of',
     &              nvars,'partitioned array ''',vname(i)(1:lvar),''''
#ifndef VERBOSE
                elseif (rec < 6) then
                   write(*,'(24x,A)') '.........................'
                endif
#endif
                                                    ! This applies to
                                                    ! both FAST and not
                do j=1,vdims(i)                     ! FAST_MODE because
                  if (dimids(j,i)==unlimdimid) then ! "start1" is used
                    start1(j)=rec                   ! for reading the
                  else                              ! whole array (fast
                    start1(j)=1                     ! mode only) and
                  endif                             ! for writing
                enddo                               ! partial arrays
                                                    ! in both modes.
#ifdef FAST_MODE
                size1=1
                do j=1,vdims(i)
                  if (dimids(j,i)==unlimdimid) then
                    count1(j)=1
                  else
                    count1(j)=dimsize(dimids(j,i))
                  endif
                  size1=size1*count1(j)
                enddo

                if (vtype(i) == nf_char) then
                  ierr=chk_alloc(size1, allc_cbf_in, 'nf_char')
                elseif (vtype(i) == nf_byte) then
                  ierr=chk_alloc(size1, allc_ibf1_in,  'nf_byte')
                elseif (vtype(i) == nf_short) then
                  ierr=chk_alloc(size1, allc_ibf2_in,  'nf_short')
                elseif (vtype(i) == nf_int) then
                  ierr=chk_alloc(size1, allc_ibff_in,    'nf_int')
                elseif (vtype(i) == nf_float) then
                  ierr=chk_alloc(size1, allc_bff4_in,  'nf_float')
                elseif (vtype(i) == nf_double) then
                  ierr=chk_alloc(size1, allc_bff8_in, 'nf_double')
                else
                  write(*,'(/1x,4A/)')   '### ERROR: Unknown ',
     &            'type of variable ''',vname(i)(1:lvar),'''.'
                  stop
                endif
                if (ierr /= 0) stop

                if (vtype(i) == nf_char) then
                  ierr=nf_get_vara_text (ncsrc, i, start1,count1,
     &                                                     cbf_in)
                elseif (vtype(i) == nf_byte) then
                  ierr=nf_get_vara_int1   (ncsrc, i, start1,count1,
     &                                                    ibf1_in)
                elseif (vtype(i) == nf_short) then
                  ierr=nf_get_vara_int2   (ncsrc, i, start1,count1,
     &                                                    ibf2_in)
                elseif (vtype(i) == nf_int) then
                  ierr=nf_get_vara_int    (ncsrc, i, start1,count1,
     &                                                    ibff_in)
                elseif (vtype(i) == nf_float) then
                  ierr=nf_get_vara_real   (ncsrc, i, start1,count1,
     &                                                    bff4_in)
                elseif (vtype(i) == nf_double) then
                  ierr=nf_get_vara_double (ncsrc, i, start1,count1,
     &                                                    bff8_in)
                endif

                if (ierr /= nf_noerr) then
                  write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A)')
     &                '### ERROR: Cannot read time record =',  rec,
     &               'of partitioned array ''',   vname(i)(1:lvar),
     &               '''', 'from netCDF file ''', src_name(1:lsrc),
     &                                   '''.',  nf_strerror(ierr)
                  write(*,*) 'start1 =', (start1(j), j=1,vdims(i))
                  write(*,*) 'count1 =', (count1(j), j=1,vdims(i))
                  goto 97
# ifdef TIMING
                else
                  nclk=3-nclk
                  call system_clock(iclk(nclk),clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  if (inc_clk < 0) inc_clk=inc_clk+clk_max
                  net_clk=net_clk+inc_clk
                  net_read_clk=net_read_clk+inc_clk
                  net_read_size=net_read_size+size1*vtypsize
# endif
                endif
#endif   /* FAST_MODE */


                do node=0,nnodes-1
                  if (vid(i,node) > 0) then
#ifdef VERBOSE
                    write(*,'(24x,A,I4,4x,3A,I3)')  'part var id =',
     &                             i, 'name = ''', vname(i)(1:lvar),
     &                                           ''' node =', node
#endif
                    size=1
                    do j=1,vdims(i)
                      if (dimids(j,i) == id_xi_rho) then
                        start(j)=xi_start(node)
                        count(j)=xi_size(node)
                        if (western_edge(node)) then
                          count(j)=count(j)+1
                        endif
                        if (eastern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == id_xi_u) then
                        start(j)=xi_start(node)
                        count(j)=xi_size(node)
                        if (.not.western_edge(node)) then
                          start(j)=start(j)-1
                        endif
                        if (eastern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == id_eta_rho) then
                        start(j)=eta_start(node)
                        count(j)=eta_size(node)
                        if (southern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                        if (northern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == id_eta_v) then
                        start(j)=eta_start(node)
                        count(j)=eta_size(node)
                        if (.not.southern_edge(node)) then
                          start(j)=start(j)-1
                        endif
                        if (northern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == unlimdimid) then
                        start(j)=rec  ; count(j)=1
                      else
                        start(j)=1 ; count(j)=dimsize(dimids(j,i))
                      endif
                      size=size*count(j)
                    enddo

                    if (vtype(i) == nf_char) then
                      ierr=chk_alloc(size, alloc_cbfr,  'nf_char')
                    elseif (vtype(i) == nf_byte) then
                      ierr=chk_alloc(size, alloc_ibfr1,  'nf_byte')
                    elseif (vtype(i) == nf_short) then
                      ierr=chk_alloc(size, alloc_ibfr2,  'nf_short')
                    elseif (vtype(i) == nf_int) then
                      ierr=chk_alloc(size, alloc_ibffr,    'nf_int')
                    elseif (vtype(i) == nf_float) then
                      ierr=chk_alloc(size, alloc_bffr4,  'nf_float')
                    elseif (vtype(i) == nf_double) then
                      ierr=chk_alloc(size, alloc_bffr8, 'nf_double')
                    else
                      write(*,'(/1x,4A/)')   '### ERROR: Unknown ',
     &                'type of variable ''',vname(i)(1:lvar),'''.'
                      stop
                    endif
                    if (ierr /= 0) stop

#ifdef FAST_MODE
                    if (vtype(i) == nf_char) then
                      call extract_text (cbf_in, cbfr, start, count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_byte) then
                      call extract_int1  (ibf1_in, ibfr1, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_short) then
                      call extract_int2  (ibf2_in, ibfr2, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_int) then
                      call extract_int   (ibff_in, ibffr, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_float) then
                      call extract_real  (bff4_in, bffr4, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_double) then
                      call extract_double(bff8_in, bffr8, start,count,
     &                                              count1, vdims(i))
                    endif
# ifdef TIMING
                    nclk=3-nclk
                    call system_clock(iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    if (inc_clk < 0) inc_clk=inc_clk+clk_max
                    net_clk=net_clk+inc_clk
                    net_assm_clk=net_assm_clk+inc_clk
# endif
#else
                    if (vtype(i) == nf_char) then
                      ierr=nf_get_vara_text (ncsrc, i, start,count,
     &                                                         cbfr)
                    elseif (vtype(i) == nf_byte) then
                      ierr=nf_get_vara_int1   (ncsrc, i, start,count,
     &                                                         ibfr1)
                    elseif (vtype(i) == nf_short) then
                      ierr=nf_get_vara_int2   (ncsrc, i, start,count,
     &                                                         ibfr2)
                    elseif (vtype(i) == nf_int) then
                      ierr=nf_get_vara_int    (ncsrc, i, start,count,
     &                                                         ibffr)
                    elseif (vtype(i) == nf_float) then
                      ierr=nf_get_vara_real   (ncsrc, i, start,count,
     &                                                         bffr4)
                    elseif (vtype(i) == nf_double) then
                      ierr=nf_get_vara_double (ncsrc, i, start,count,
     &                                                         bffr8)
                    endif

                    if (ierr /= nf_noerr) then
                      write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A)')
     &                 '### ERROR: Cannot read time record =',  rec,
     &                 'of partitioned array ''',   vname(i)(1:lvar),
     &                 '''', 'from netCDF file ''', src_name(1:lsrc),
     &                                     '''.',  nf_strerror(ierr)
                      write(*,*) 'start =', (start(j), j=1,vdims(i))
                      write(*,*) 'count =', (count(j), j=1,vdims(i))
                      goto 97
# ifdef TIMING
                    else
                      nclk=3-nclk
                      call system_clock(iclk(nclk),clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      if (inc_clk < 0) inc_clk=inc_clk+clk_max
                      net_clk=net_clk+inc_clk
                      net_read_clk=net_read_clk+inc_clk
                      net_read_size=net_read_size+size*vtypsize
# endif
                    endif
#endif /* FAST_MODE */

                    if (vtype(i) == nf_char) then
                      ierr=nf_put_vara_text (ncid(node), vid(i,node),
     &                                           start1, count, cbfr)
                    elseif (vtype(i) == nf_byte) then
                      ierr=nf_put_vara_int1  (ncid(node), vid(i,node),
     &                                         start1,  count, ibfr1)
                    elseif (vtype(i) == nf_short) then
                      ierr=nf_put_vara_int2  (ncid(node), vid(i,node),
     &                                          start1, count, ibfr2)
                    elseif (vtype(i) == nf_int) then
                      ierr=nf_put_vara_int   (ncid(node), vid(i,node),
     &                                        start1, count, ibffr)
                    elseif (vtype(i) == nf_float) then
                      ierr=nf_put_vara_real  (ncid(node), vid(i,node),
     &                                          start1, count, bffr4)
                    elseif (vtype(i) == nf_double) then
                      ierr=nf_put_vara_double(ncid(node), vid(i,node),
     &                                          start1, count, bffr8)
                    endif

                    if (ierr /= nf_noerr) then
                      write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A)')
     &                   '### ERROR: Cannot write time record =',  rec,
     &                     'of partitioned array ''', vname(i)(1:lvar),
     &                '''', 'into netCDF file ''', fname(node)(1:lfnm),
     &                                        '''.', nf_strerror(ierr)
                      write(*,*) 'start1 =', (start1(j), j=1,vdims(i))
                      write(*,*) ' count =',  (count(j), j=1,vdims(i))
                      goto 97
#ifdef TIMING
                    else
                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      if (inc_clk < 0) inc_clk=inc_clk+clk_max
                      net_clk=net_clk+inc_clk
                      net_wrt_clk=net_wrt_clk+inc_clk
                      net_wrt_size=net_wrt_size+size*vtypsize
#endif
                    endif
                  endif       !<--  vid(i,node) > 0
                enddo       ! <-- node=0,nnodes-1
              endif       ! <-- part_type, series, 3-way
            endif       ! <--series(i) .or. rec==1
          enddo       ! <-- i=1,nvars
        enddo       ! <-- rec=1,max(tsize,1)

! Close all netCDF files

  97    write(*,*) 'closing files...'
        ierr=nf_close(ncsrc)
#if defined TIMING
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_gray_clk=net_gray_clk+inc_clk
#endif
        write(*,*) '...........input'

        do node=0,nnodes-1
         ierr=nf_close(ncid(node))
        enddo
#if defined TIMING
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_sync_clk=net_sync_clk+inc_clk
#endif
        write(*,*) '..........output'
      enddo       ! <-- arg=3,nargs

#ifdef TIMING
      call cpu_time(run_time) ; run_time=run_time-tstart
      write(*,'(/1x,A,F11.2,1x,A/)') 'cpu time =', run_time, 'sec'

      if (clk_rate > 0) then
        ReadSize=1.0D-6*dble(net_read_size)
        WrtSize=1.0D-6*dble(net_wrt_size)
        ReadTime=net_read_clk/dble(clk_rate)
        WrtTime = net_wrt_clk/dble(clk_rate)
        InitTime=net_init_clk/dble(clk_rate)
        SyncTime=net_sync_clk/dble(clk_rate)

        write(*,'(/1x,A,22x,F12.2,1x,A)') 'Analysis/file creation :',
     &                                               InitTime, 'sec'
        write(*,'(8x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &         'Total data read :', ReadSize, 'MBytes in',  ReadTime,
     &                          'sec (', ReadSize/ReadTime, 'MB/sec)'
        write(*,'(5x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &      'Total data written :', WrtSize,  'MBytes in',   WrtTime,
     &                          'sec (',  WrtSize/WrtTime,  'MB/sec)'
# ifdef FAST_MODE
        AssmTime=net_assm_clk/dble(clk_rate)
        write(*,'(3x,A,22x,F12.2,1x,A)')    'Array splitting time :',
     &                                               AssmTime, 'sec'
# endif
        write(*,'(2x,A,22x,F12.2,1x,A)')   'Output file sync time :',
     &                                               SyncTime, 'sec'

        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_gray_clk=net_gray_clk+inc_clk
        GrayTime=dble(net_gray_clk)/dble(clk_rate)
        write(*,'(14x,A,22x,F12.2,1x,A)') 'Gray time :',GrayTime,'sec'
        write(*,'(47x,A/12x,A,11x,F12.2,1x,A/)') '------------------',
     &   'Elapsed wall-clock time:',dble(net_clk)/dble(clk_rate),'sec'
      endif
#endif
      stop
      end

! Setup horizontal dimensions and associated variables for each
! subdomain. The following code is extracted into a separate entity,
! which is mathematically consistent with the actual mpi_setup, and
! compute_starts_and_counts.h.  In essence, this subroutine receives
! four integer numbers: dimension of the grid (including boundary
! rows on the side) -- xi_rho,eta_rho; and number of partitions in
! each direction, NP_XI,NP_ETA.  These four are translated into
! dimensions of subdomains, Lm,Mm, bounds of used portions of arrays
! iwest,ieast,jsouth,jnorth and global-to-relative index shift
! translations, SW_corn,jSW_corn (exactly the same way as in
! mpi_setup.F), which are then further translated into
! xi_start(node),xi_size(node), and eta_start(node),eta_size(node),
! which have meaning of starting netCDF indices for RHO-point sub-
! array in netCDF file belonging to each individual subdomain, and
! the sizes of subarrays.  The other four variables, western_,
! eastern_, southern_, and northern_edge are logical flags to
! identify the proximity of side boundary on each side for each
! subdomain.

! Note that the code above - the main part of "partit" - is written
! in such a way that it makes no assumption about the structure of the
! "processor grid" (the arrangement of subdomains corresponding to MPI
! nodes relatively to the physical grid), but rather relies
! exclusively on the eight variables defined in the code below.


      subroutine mpi_setup (NP_XI,NP_ETA, xi_rho,eta_rho,
     &               xi_start,xi_size, eta_start,eta_size,
     &                        western_edge, eastern_edge,
     &                       southern_edge, northern_edge)
      implicit none
      integer, intent(in) :: NP_XI,NP_ETA, xi_rho,eta_rho
! out
      integer, intent(out), dimension(0:NP_XI*NP_ETA-1) ::
     &              xi_start, xi_size, eta_start, eta_size
      logical, intent(out), dimension(0:NP_XI*NP_ETA-1) ::
     & western_edge, eastern_edge, southern_edge, northern_edge
! internal
      integer LLm,Lm, MMm,Mm, nnodes,node, inode,jnode,
     &        iwest,ieast,jsouth,jnorth, iSW_corn,jSW_corn,
     &                                   off_XI,off_ETA

      nnodes=NP_XI*NP_ETA
      LLm=xi_rho-2           ; MMm=eta_rho-2
      Lm=(LLm+NP_XI-1)/NP_XI ; Mm=(MMm+NP_ETA-1)/NP_ETA
#ifdef VERBOSE
        write(*,'(4x,2A,2I5,4x,A,2I5)') 'mpi_setup: found grid ',
     &          'sizes: LLm,MMm =',  LLm,MMm,  'Lm,Mm =', Lm,Mm
#endif


      do node=0,nnodes-1
#ifdef TRANSPOSED_MPI_NODE_ORDER
        inode=node/NP_ETA        ! The following code segment
        jnode=node-inode*NP_ETA  ! maps exactly onto logic inside
#else
        jnode=node/NP_XI         ! "mpi_setup" in the actual ROMS
        inode=node-jnode*NP_XI   ! code.
#endif

        off_XI=NP_XI*Lm-LLm
        iSW_corn=inode*Lm-off_XI/2
        if (inode == 0) then
          iwest=1+off_XI/2
        else
          iwest=1
        endif
        if (inode < NP_XI-1) then
          ieast=Lm
        else
          ieast=Lm -(off_XI+1)/2
        endif

        off_ETA=NP_ETA*Mm-MMm
        jSW_corn=jnode*Mm-off_ETA/2
        if (jnode == 0) then
          jsouth=1+off_ETA/2
        else
          jsouth=1
        endif
        if (jnode < NP_ETA-1) then
          jnorth=Mm
        else
          jnorth=Mm -(off_ETA+1)/2
        endif

                                           ! These two correspond to
        xi_size(node)=ieast-iwest+1        ! "xi_rho" and "eta_rho" in
        eta_size(node)=jnorth-jsouth+1     ! "mpi_setup" inside
                                           ! #ifdef PARALLEL_FILES
        if (inode == 0) then               ! CPP-branch
          xi_start(node)=iSW_corn+iwest
        else
          xi_start(node)=iSW_corn+iwest+1
        endif

        if (jnode == 0) then
          eta_start(node)=jSW_corn+jsouth
        else
          eta_start(node)=jSW_corn+jsouth+1
        endif

        if (inode == 0) then
          western_edge(node)=.true.
        else
          western_edge(node)=.false.
        endif
        if (inode < NP_XI-1) then
          eastern_edge(node)=.false.
        else
          eastern_edge(node)=.true.
        endif
        if (jnode == 0) then
          southern_edge(node)=.true.
        else
          southern_edge(node)=.false.
        endif
        if (jnode < NP_ETA-1) then
          northern_edge(node)=.false.
        else
          northern_edge(node)=.true.
        endif

#ifdef VERBOSE
        write(*,'(A,3I5,2x,A,4(1x,L1))') 'node =', node, inode,jnode,
     &                    'western,eastern,southern,northern_edge =',
     &                      western_edge(node),   eastern_edge(node),
     &                      southern_edge(node), northern_edge(node)
#endif
      enddo   !--> discard iwest,ieast,jsouth,jnorth
      end


      function chk_alloc(size, alloc_size, type_str)
      implicit none
      integer chk_alloc, size, alloc_size
      character(len=*) type_str
      if (size > alloc_size) then
        write(*,'(/1x,4A/12x,A,I10,2x,A,I10/)')          '### ERROR: ',
     & 'Insufficient size of allocated workspace of type ''', type_str,
     &     ''':', 'needed', size, 'but actually allocated', alloc_size
        chk_alloc=1
      else
       chk_alloc=0
      endif
      end

! The following six routines extract rectangular block of data "buff"
! from the whole-grid array "bfr_in"; all the routines have common core
! and differ only by type declarations.

#ifdef FAST_MODE
# define PARTIT

      subroutine extract_text (bfr_in, bffr, start,count,
     &                                      count1,vdims)
      implicit none
      character(len=1) bfr_in(*), bffr(*)
# include "copy_subarray.h"
      end

      subroutine extract_int1 (bfr_in, bffr, start,count,
     &                                     count1, vdims)
      implicit none
      integer(kind=1) ::  bfr_in(*), bffr(*)
# include "copy_subarray.h"
      end

      subroutine extract_int2 (bfr_in, bffr, start,count,
     &                                     count1, vdims)
      implicit none
      integer(kind=2) :: bfr_in(*), bffr(*)
# include "copy_subarray.h"
      end

      subroutine extract_int   (bfr_in, bffr, start,count,
     &                                      count1, vdims)
      implicit none
      integer(kind=4) :: bfr_in(*), bffr(*)
# include "copy_subarray.h"
      end

      subroutine extract_real   (bfr_in, bffr, start,count,
     &                                      count1, vdims)
      implicit none
      real(kind=4) :: bfr_in(*), bffr(*)
# include "copy_subarray.h"
      end

      subroutine extract_double (bfr_in, bffr, start,count,
     &                                      count1, vdims)
      implicit none
      real(kind=8) :: bfr_in(*), bffr(*)
# include "copy_subarray.h"
      end
#endif
