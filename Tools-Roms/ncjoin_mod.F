      module ncjoin_mod  ! contains variables and generic routines to be used in both
                         ! ncjoin and ncjoin_mpi to avoid duplications and for clarity.

      use mpi
      use netcdf  ! for nf90_ calls

      implicit none

![ CPP DEFS:
! Delete partial files after joining:
#define DEL_PART_FILES
! Document program execution time:
#define TIMING
! Verbose terminal output:
c--#define VERBOSE
! Verbose terminal output for MPI specifics:
c--#define MPIVERBOSE
!]

      ! INCLUDES:
      ! remove netcdf.inc once fully converted to nf90_ routines
#include "netcdf.inc"

      ![ DECLARATION OF VARIABLES:

      integer :: n_x=1, n_y=1, n_z=1 ! number of gpoints in x, y & z (temporarily 1)

      ! MPI vars:
      integer :: n_procs, my_rank    ! number of processors, rank of each processor
      integer :: nprocs_x            ! number of mpi processes in x
      integer :: nprocs_y            ! number of mpi processes in y
      integer :: iproc, jproc        ! MPI proc sub-domain positions in x and y
      integer :: p_indx               ! loop index over processes

      ! Chunking vars:
      integer, dimension(:), allocatable  :: chunksize_x, chunksize_y ! Chunk sizes in x & y dimensions per variable
      integer, dimension(:), allocatable  :: n_chunks_x, n_chunks_y   ! no. of chunks to cover each dim of full domain per var
      integer, dimension(:), allocatable  :: p_chunks_x, p_chunks_y   ! no. of chunks to cover each dim of sub-domain (mpi proc)
      integer                             :: p_points_x, p_points_y   ! no. of grid points for each dim of proc sub-domain
      integer, dimension(4)               :: set_req_chunk_sizes      ! Set required chunk sizes in each dim
      integer, dimension(4)               :: chunk_sizes              ! Confirm chunk sizes in newly defined variable

      real, allocatable, dimension(:,:,:) :: data                     ! store 2D/3D data

      ! R1) Read in: Index ranges as per each partitioned input file.
      integer, dimension(4) :: start, count                      ! netcdf vars for non-partitioned variables
      integer, dimension(:,:), allocatable :: start_in, count_in ! reading start & count ranges for each dim of 2D/3D arrays
      integer :: x_start_in   ! temp var to adjust for xi_rho  vs xi_u  points of read partial file
      integer :: y_start_in   ! temp var to adjust for eta_rho vs eta_v points of read partial file
      integer :: x_in_dimsize ! temp var to store x dimension size of read partial file
      integer :: y_in_dimsize ! temp var to store y dimension size of read partial file

      ! R2) Index ranges to store read in data in proc's temp 'data' array
      integer, dimension(:), allocatable :: data_x_in_str, data_x_in_end ! proc stores data in array 'data' with local (not global)
      integer, dimension(:), allocatable :: data_y_in_str, data_y_in_end ! indices, hence set inbound indices for proc's data array.

      ! W1) Index ranges from data in temp array to write to joined output file
      integer :: p_str_out_x, p_str_out_y ! mpi (sub-domain) proc's temporary 'data' array start x and y grid points
      integer :: p_end_out_x, p_end_out_y ! mpi (sub-domain) proc's temporary 'data' array end   x and y grid points

      ! W2) Write out: Index ranges as per full model domain.
      integer, dimension(4) :: start_out, count_out ! proc's data indices relative to global output file

      ! A)  Assign read in partitions to procs:
      ! Each mpi process will loop through the nodes (partitions) it is responsible for.
      integer, dimension(:), allocatable :: lnodes_in_proc ! list of nodes (partial files) read into proc
      integer :: nnodes_in_proc  ! number of nodes into proc for looping.
      integer :: p_in_nodes      ! counter for nodes linked to proc. nnodes_in_proc=p_in_nodes
      integer :: in2p            ! looping index for nodes read in by proc
      logical :: sw_corner, se_corner, nw_corner, ne_corner ! if node's corner is contained within proc range
      integer :: x_nc_sw,x_nc_se,x_nc_nw,x_nc_ne            ! x index of corners of input node
      integer :: y_nc_sw,y_nc_se,y_nc_nw,y_nc_ne            ! y index of corners of input node

      ! Misc vars:
      character(len=4)  :: char_arg       ! Store program input argument
      integer           :: d              ! loop index over dimensions
      integer           :: v_dimsize      ! Size of current dimension length
      character(len=1)  :: char_in        ! Text read in. I.e. scalar 'spherical' variable
      character(len=15) :: var_name_debug ! debug. should delete when confirmed collective/independent I/O
      integer           :: scalar_int     ! variable for joining scalar integers
      real              :: scalar_real    ! variable for joining scalar reals (float or double)
      logical           :: mayday         ! flag to break out of subroutine with error
      integer           :: guess_nnodes   ! guessed number of input nodes (files) linked to mpi proc for storage array size
      integer           :: deflate_level  ! compression level for joined file
      logical           :: shuffle        ! shuffle on for extra compression (nc3to4z had this).

      ! ------ CHILD_BRY VARS -------------
      integer :: global_south_bry_size    ! size of join dimension in child boundary file
      integer :: global_west_bry_size
      integer, dimension(4) :: vid_out    ! variable ID's in joined output file

      ! ------ ORIGINAL NCJOIN VARS -------

      integer, parameter :: maxdims=32  ! max allowed dimensions in file (change as required)
      integer, parameter :: maxvars=220 ! max allowed vars. BGC has 201 vars in restart file

      logical complete, clean_set,  digit, var_mask,    lnewvar
      integer nargs, nnodes,  size_XI,  XI_rho, id_xi_rho,  id_xi_u,
     &        arg,   node,    size_ETA, ETA_rho,id_eta_rho, id_eta_v,
     &        ierr,  maxnodes,size_S,   tsize,  unlimdimid, rec,
     &        ntest, nctarg,  ndims,    size,   code_size,  lvar,
     &        nvars, ngatts,  varatts,  size1,  code_size_bak,
     &        i,j,k, is,ie,   lncn,     ltrg,   lstr, lbak ! , lenstr 0325

      integer, external :: lenstr    ! 0325 addition

      character(len=8) sffx, sffx_bak
      character(len=32) vname(maxvars), dimname(maxdims) ! Variable and dimension names
      character(len=64) nctestname, nctargname, root, root_bak, string
      character(len=64), dimension(:), allocatable :: ncname

      integer, dimension(:),   allocatable :: ncid, xi_start, eta_start
      integer, dimension(:,:), allocatable :: vid, dimsize
      logical, dimension(:),   allocatable :: western_edge,  eastern_edge,
     &                                        southern_edge, northern_edge

      logical series(maxvars)
      integer, dimension(maxvars) :: varid, vnode, vdims, vartype, part_type
      integer, dimension(maxdims) :: dimid, ldim,  ibuff, start1  ! DevinD I doubt start1 is still needed
      integer, dimension(maxdims,maxvars) :: dimids=1 ! DevinD set to value for /=tsize to work
      integer max_buff_size, alloc_buff_size
      real*8, allocatable, dimension(:) :: buff       ! DevinD now just using for 1D arrays

#ifdef DEL_PART_FILES
      logical del_part_files
      character(len=128) rmcmd
#endif
#ifdef TIMING
      real*4 tstart, RUN_time, CPU_time(2)
      integer iclk(2), nclk, clk_rate, clk_max, iclk_init
      integer*8 net_read_size, net_wrt_size, net_fcrt_clk,
     &          net_read_clk,  net_wrt_clk,  net_assm_clk,
     &          net_sync_clk,  net_gray_clk, inc_clk
      real*8 ReadSize, ReadTime, WrtSize,  WrtTime,
     &       FcrtTime, AssmTime, SyncTime, GrayTime

# ifdef DEL_PART_FILES
      integer*8 net_rmcmd_clk
# endif
#endif

      !]

      contains

! ----------------------------------------------------------------------
      subroutine init_timing_and_vars  ![
      implicit none

#ifdef TIMING
      call etime(CPU_time, tstart)
      nclk=1
      call system_clock (iclk(nclk), clk_rate, clk_max)
      iclk_init=iclk(nclk)
      net_read_clk=0       ! Set all timing counters
      net_read_size=0
      net_wrt_size=0
      net_wrt_clk =0
      net_sync_clk=0
      net_assm_clk=0
      net_gray_clk=0
      net_fcrt_clk=0
#endif
#ifdef DEL_PART_FILES
      del_part_files=.false.
# ifdef TIMING
      net_rmcmd_clk=0
# endif
#endif
      ntest=-1                         ! initialize sizes of buffer
      maxnodes=-1                      ! arrays to be allocated. Here
      max_buff_size=0                  ! "max_*" means the needed size,
      alloc_buff_size=0                ! and "alloc_*" is the size of
                                       ! the actually allocated array.
      nargs=iargc()  ! Intrinsic function to get total number of arguements
      arg=2          ! start at 2+1 since first 2 arguements are mpi sub-domain values np_x & np_y

      end subroutine init_timing_and_vars  !]

! ----------------------------------------------------------------------
      subroutine setup_mpi  ![
                            !  setup MPI system and
      implicit none         !  error checking of user MPI config

      call MPI_Init(ierr)                               ! initialize MPI
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr) ! get proc's rank (rank starts from 0)
      call MPI_Comm_size(MPI_COMM_WORLD, n_procs, ierr) ! get total number of processors

      if(nargs < 4) then   ! Incorrect program arguments: give instructions
        if(my_rank==0) then
          write(*,'(/1x,A//10x,2A)') 'Correct usage of ncjoin_mpi:',
     &      'mpiexec  -n  np  ncjoin_mpi  ',
     &      'np_x  np_y  his.0000.*.nc'
          write(*,'(1x,A/10x,2A/)') 'or',
     &      'mpiexec  -n  np  ncjoin_mpi  ',
     &      'np_x  np_y  -d  his.0000.*.nc'
          write (*,'(/1x,2A/1x,2A//1x,2A/1x,A//1x,A//1x,A//)')
     &      'np_x and np_y are the ',
     &      'number of sub-domains you choose in x and y...',
     &      '...for ncjoin_mpi not your input partition files, ',
     &      'they can be different!',
     &      'For efficiency, try to keep keep the ratio of np_x to ',
     &      'np_y similar ','to the sub-domaining of your partitions.',
     &      '-d deletes the partitioned files after joining.',
     &      'ncjoin_mpi can only be used on expanse!'
        endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      endif

      call getarg(1,char_arg)    ! get proc sub-domains from arguements
      read(char_arg,*) nprocs_x  ! e.g. ncjoin_mpi 8 6 ....
      call getarg(2,char_arg)    ! read program argument as text
      read(char_arg,*) nprocs_y  ! convert argument to integer

      if (n_procs /= nprocs_x*nprocs_y) then  ! Abort if wrong mpi sub-domain tiling
        write(*,'(/1x,2A/)') 'Error: number of mpi ',
     &             'processes /= nprocs_x * nprocs_y (sub-domains)'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      endif

      if(my_rank==0) write(*,'(/1x,A,/1x,A,I3,2(/1x,A,I2))')
     &                'MPI procs:', 'n_procs= ', n_procs,
     &                'nprocs_x= ', nprocs_x, 'nprocs_y= ', nprocs_y

      end subroutine setup_mpi !]

! ----------------------------------------------------------------------
      subroutine check_partial_file_set  ![
                                         !  determine if there is a complete
      implicit none                      !  set of partial files to join

  1     nnodes=-1
        root_bak(1:1)=' '  ! Reset variables which
        sffx_bak(1:1)=' '  ! define the file set.
        code_size_bak=-1

  2     arg=arg+1 ! cycle through argument list of names of partitioned files by
                  ! coming back here with 'goto 2' after each filename.
                  ! Struggled to remove all goto's here. Can atleast use 'do while'
                  ! with 'cycle'.

        call getarg(arg,nctestname) ! get arguements from calling ncjoin arg1 ...
        lncn=lenstr(nctestname)     ! e.g. ncjoin his.0000.*.nc will actually send in his.0000.00.nc his.0000.01.nc ...
#ifdef DEL_PART_FILES
        if (arg.eq.3 .and. (lncn.eq.2 .and. nctestname(1:2).eq.'-d'
     &     .or. lncn.eq.8 .and. nctestname(1:8).eq.'--delete') ) then
           write(*,'(/1x,2A/)') '>>>> Flag to delete partial files ', ! Already in master
     &                                                  'is raised.'  ! only section
           del_part_files=.true.
          goto 2
        endif
#endif
        if (ntest.ne.-1) then
          ierr=nf_close(ntest)
          ntest=-1
        endif
        ierr=nf_open (nctestname, nf_nowrite, ntest)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_att (ntest, nf_global, 'partition', i, lvar) ! e.g. :partition = 0, 6, 1, 1 ;
                                                                   ! i=attribute type, lvar= number of values stored in attribute.
          if (ierr .eq. nf_noerr) then
            if (i.eq.nf_int .and. lvar.eq.4) then
              ierr=nf_get_att_int (ntest,nf_global, 'partition', ibuff) ! ibuff array stores the 4 integers of 'partition'.
              if (ierr .eq. nf_noerr) then
                if (nnodes.eq.-1) then
                  nnodes=ibuff(2)
                  if (nnodes.gt.maxnodes) then
                    maxnodes=nnodes
                    if (allocated(ncid)) then
                      deallocate(dimsize)
                      deallocate(vid)
                      deallocate(northern_edge)
                      deallocate(southern_edge)
                      deallocate(eastern_edge)
                      deallocate(western_edge)
                      deallocate(eta_start)
                      deallocate(xi_start)
                      deallocate(ncname)
                      deallocate(ncid)
                    endif
                    allocate (ncid(0:nnodes-1))
                    allocate (ncname(0:nnodes-1))
                    allocate (xi_start(0:nnodes-1))
                    allocate (eta_start(0:nnodes-1))
                    allocate (western_edge(0:nnodes-1))
                    allocate (eastern_edge(0:nnodes-1))
                    allocate (southern_edge(0:nnodes-1))
                    allocate (northern_edge(0:nnodes-1))
                    allocate (vid(maxvars,0:nnodes-1))
                    allocate (dimsize(maxdims,0:nnodes))
                  endif
                                          ! Reset variables defining
                  complete=.false.        ! a complete set of partitial
                  do node=0,nnodes-1      ! files. These variables will
                    ncid(node)=-1         ! receive meaningful values
                    xi_start(node)=-1     ! from data read from netCDF
                    eta_start(node)=-1    ! file headers, subsequently
                  enddo                   ! be used to verify the set
                                          ! completeness.

                elseif (nnodes.ne.ibuff(2)) then
                  write(*,'(/1x,2A,I4/14x,3A/14x,A,I4,4x,A/)')
     &                 '### WARNING: Number of MPI nodes in global ',
     &                 'attribute ''partition'', nnodes =', ibuff(2),
     &                 'in netCDF file ''',       nctestname(1:lncn),
     &                 ''' contradicts that from the initial',
     &                 'file in the sequence, nnodes =',      nnodes,
     &                                   ' ==> The file is ignored.'
                  arg=arg-1
                  goto 5
                endif

                node=ibuff(1)
                if (ncid(node).ne.-1) then
                  write(*,'(/1x,2A,I4,1x,A)') '### ERROR: netCDF ID ',
     &                         'for file corresponding to MPI-node =',
     &                                   node,  'is already in use.'
                  stop
                endif


                if (ncid(node).eq.-1 .and. xi_start(node).eq.-1
     &                         .and.  eta_start(node).eq.-1) then
                  ncid(node)=ntest
                  ncname(node)=nctestname
                  xi_start(node)=ibuff(3)  ! DevinD: i SW of partition WRT global numbering
                  eta_start(node)=ibuff(4) ! DevinD: j SW of partition WRT global numbering

#define ntest illegal


! Lexical analysis of the file name: It is assumed that name of the
! file consists of root name (eg.: "history"); integer number which
! contains MPI node number (eg.: 03) and; suffix (eg.: ".nc").
! Files which belong to the same set normally have the same (1) root
! and (2) suffix names; the same (3) number of digits in the MPI node
! number segment in the filename and; (4) MPI node number from the
! file name should match the number determined from global attribute
! "partition".
                                               ! Determine positions
                  digit=.false.                ! of starting and ending
                  is=0                         ! characters of MPI node
                  ie=0                         ! segment (is:ie)
                  i=lncn+1
                  do while (is.eq.0 .and. i.gt.1)
                    i=i-1
                    if (nctestname(i:i).ge.'0' .and.
     &                  nctestname(i:i).le.'9') then
                      if (.not.digit) then
                        if (i.lt.lncn) then
                          if (nctestname(i+1:i+1).eq.'.') then
                            ie=i
                            digit=.true.         ! check that node
                          endif                  ! segment and suffix
                        else                     ! are separated by '.'
                          ie=i                   ! no-suffix case
                          digit=.true.
                        endif
                      endif
                    elseif (digit .and. nctestname(i:i).eq.'.') then
                      digit=.false.
                      is=i+1
                    endif
                  enddo

                  if (is.gt.0 .and. ie.ge.is) then
                    root=nctestname(1:is-1)
                    if (ie.lt.lncn) then         ! Extract common
                      sffx=nctestname(ie+1:lncn) ! part of file names,
                    else                         ! MPI node number
                      sffx(1:1)=' '              ! and suffix (if any)
                    endif
                    k=0
                    do i=is,ie
                      k=10*k + ichar(nctestname(i:i))-48
                    enddo
                    code_size=ie-is+1
                  else
                    write(*,'(/1x,3A/)')        '### ERROR: Cannot ',
     &                  'determine MPI node number from filename ''',
     &                                     nctestname(1:lncn), '''.'
                  endif
# ifdef VERBOSE
                  if(my_rank==0)
     &              write(*,'(1x,3A,I3,1x,A,I4,1x,A,I4,3x,A,2I4)')
     &             'fname = ''', nctestname(1:lncn),  ''' code_size =',
     &              code_size,  'code =', k, 'node =', node, 'i,jSW =',
     &              xi_start(node), eta_start(node)
# endif

! Checking consistency of root name with previously found.

                  ierr=nf_noerr
                  if (root_bak(1:1).eq.' ') then
                    root_bak=root
                  else
                    lvar=lenstr(root)
                    lbak=lenstr(root_bak)
                    if (lvar.ne.lbak .or. root.ne.root_bak) then
                      ierr=ierr+1
                      write(*,'(/8x,6A/17x,3A/)') 'WARNING: file ''',
     &                     nctestname(1:lncn),   ''' has different ',
     &                    'root name ''',  root(1:lvar),   ''' than',
     &                    'previously found root name ''',
     &                     root_bak(1:lbak), ''' from the same set.'
                    endif
                  endif

! Checking consistency of suffix with previously found..

                  if (sffx_bak(1:1).eq.' ') then
                    sffx_bak=sffx
                  else
                    lvar=lenstr(sffx)
                    lbak=lenstr(sffx_bak)
                    if (lvar.ne.lbak .or. sffx.ne.sffx_bak) then
                      ierr=ierr+1
                      write(*,'(/8x,7A/17x,3A/)')       'WARNING: ',
     &                  'file ''',  nctestname(1:lncn),   ''' has ',
     &                  'different suffix name ''',    sffx(1:lvar),
     &                  ''' than','previously found suffix name ''',
     &                  sffx_bak(1:lbak),   ''' from the same set.'
                    endif
                  endif

! Checking consistency of length of node number segment

                  if (code_size_bak.eq.-1) then
                    code_size_bak=code_size
                  elseif (code_size .ne. code_size_bak) then
                    ierr=ierr+1
                    write(*,'(/8x,A,I2,1x,A/17x,3A,I2,A/)')
     &              'WARNING: number of digits in MPI node segment',
     &               code_size, 'in filename', '''',
     &               nctestname(1:lncn),
     &                ''' is different than previously determined',
     &                                     code_size_bak, '.'
                  endif

! Checking consistency of node number with the file name.

                  if (k.ne.node) then
                    ierr=ierr+1
                    write(*,'(/8x,3A,I3/17x,2A/17x,A,I3,A/)')
     &                   'WARNING: file ''', nctestname(1:lncn),
     &                   ''' belongs to different MPI node',   node,
     &                   '(as determined from its global attribute',
     &                   '''partition'')', 'than node', k,
     &                   ' determined from to the file name.'
                  endif

! Stop, if something is wrong.

                  if (ierr.ne.nf_noerr) then
                    mayday=.true.; return ! goto 97 can't have go to
                  endif
                else
                  arg=arg-1
                  goto 5
                endif
              else
                write(*,'(/1x,2A/14x,3A/)')   '### WARNING: Cannot ',
     &           'aquire global attribute ''partition'' from netCDF',
     &                                 'file ''', nctestname(1:lncn),
     &                               '''. ==> This file is ignored.'
              endif
            else
              write(*,'(/1x,2A/14x,3A/)')  '### WARNING: Wrong type ',
     &                'or size of global attribute ''partition'' in ',
     &                           'netCDF file ''', nctestname(1:lncn),
     &                                '''. ==> This file is ignored.'
            endif
          else
            write(*,'(/1x,3A/)') '### WARNING: ''', nctestname(1:lncn),
     &      ''' is not a partial netCDF file: ==> The file is ignored.'
            if (arg<nargs) then  ! Not last file
              goto 2             ! -> next file
            else
              write(*,'(/1x,A/)') 'Final file. Terminating ncjoin_mpi.'
              complete=.false.   ! Probably not needed anymore... Use this to skip parallel read/write section of code.
              mayday=.true.; return ! goto 23            ! goto end of mater only pre-processing.
            endif
          endif
        else
          write(*,'(/1x,4A/14x,A/)')   '### WARNING: Cannot open ''',
     &                   nctestname(1:lncn), ''' as a netCDF file: ',
     &                nf_strerror(ierr), ' ==> The file is ignored.'
        endif

#define nctestname illegal

! Verify, whether ncname(0:nnodes-1) and ncid(0:nnodes-1) > 0 (i.e.,
! successfully opened for reading) comprise a complete set of partial
! files.  Keep searching, if not.

   5    continue
        if (nnodes.gt.0) then
          complete=.true.
          do node=0,nnodes-1
            if (ncid(node).lt.0) complete=.false.
          enddo
        endif

      if (.not.complete .and. arg.lt.nargs) goto 2  !--> next file


! Once a complete set is identified, print the finenames.

      if (complete) then
        lncn=lenstr(ncname(0))
        if(my_rank==0) write(*,'(2(1x,A,I4),1x,A,2x,A,2I5)')   'Processing set of ',
     &                         nnodes, 'files', 0, ncname(0)(1:lncn),
     &                          'i,jSW =', xi_start(0), eta_start(0)
        do node=1,nnodes-1
          if (node.lt.16 .or. (nnodes.gt.16 .and.
     &                         node.eq.nnodes-1 )) then
            if(my_rank==0) write(*,'(29x,I4,1x,A,2x,A,2I5)') node,
     &                  ncname(node)(1:lncn), 'i,jSW =',
     &                  xi_start(node), eta_start(node)
          elseif (nnodes.gt.16 .and. node.lt.18) then
            if(my_rank==0) write(*,'(24x,A)') '.................................'
          endif
        enddo

#undef ntest
        if (ntest.ne.-1) then       ! Thus far netCDF file id array
          ierr=nf_close(ntest)      ! "ncid(0:nnodes-1)" was used just
          ntest=-1                  ! to signal that a complete set of
        endif                       ! partitioned files has been
        do node=0,nnodes-1          ! identified, but all the files are
          ncid(node)=-1             ! actually closed at this moment.
        enddo                       ! Reset the ids accordingly.

      elseif (arg.lt.nargs) then
        goto 1
      else
        write(*,*) 'stop at 466'
        stop
      endif

      end subroutine check_partial_file_set  !]

! ----------------------------------------------------------------------
      subroutine check_ndims_nvars_natts( node )  ![
                                                  !  check consistent number of dimensions,
      implicit none                               !  variables and attributes.

      ! Verify that ndims, ngatts, unlimdimid are the same for all nodes,
      ! however, note that different files may store different composition
      ! of variables, and netCDF variable IDs for the same variable (with
      ! the same name) may be different across the set of files.

      integer, intent(in) :: node

            ierr=nf_inq (ncid(node), ibuff(1), ibuff(2),  ! returns: 1-ndims,  2-nvars
     &                               ibuff(3), ibuff(4))  !          2-ngatts, 4-unlimid
            if (ierr .ne. nf_noerr) then
              write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot make ',
     &                        'general inquiry into netCDF file ''',
     &               ncname(node)(1:lncn), '''.', nf_strerror(ierr)
              mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
            elseif (ibuff(1) .gt. maxdims) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of dimensions', ibuff(1), 'in netCDF file ''',
     &          ncname(node)(1:lncn),    '''',   'exceeds limit.  ',
     &             'Increase parameter maxdims in file "ncjoin.F".'
              mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
            elseif (ibuff(2) .gt. maxvars) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of variables',  ibuff(2), 'in netCDF file ''',
     &          ncname(node)(1:lncn),  '''',      'exceeds limit. ',
     &             'Increase parameter maxvars in file "ncjoin.F".'
              mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
            elseif (node.eq.0) then
              ndims=ibuff(1)
c**           nvars=ibuff(2)
              ngatts=ibuff(3)
              unlimdimid=ibuff(4)
            else
              if (ibuff(1) .ne. ndims) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                    'file ''', ncname(node)(1:lncn), ''' has ',
     &                    'different number of dimensions than ''',
     &                                       ncname(0)(1:lstr), '''.'
                ierr=ierr+1
              endif
c**           if (ibuff(2) .ne. nvars) then
c**             write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
c**  &                    'file ''', ncname(node)(1:lncn), ''' has ',
c**  &                       'different number of variables than ''',
c**  &                                      ncname(0)(1:lstr), '''.'
c**             ierr=ierr+1
c**           endif
              if (ibuff(3) .ne. ngatts) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                   'file ''',  ncname(node)(1:lncn), ''' has ',
     &               'different number of global attributes than ''',
     &                                       ncname(0)(1:lstr),'''.'
                ierr=ierr+1
              endif
              if (ibuff(4) .ne. unlimdimid) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                    'file ''', ncname(node)(1:lncn), ''' has ',
     &                'different ID for unlimited dimension than ''',
     &                                      ncname(0)(1:lstr), '''.'
                ierr=ierr+1
              endif
              if (ierr .ne. nf_noerr) then
                mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
              endif
            endif

! Verify that the sequence of dimension names is consistent
! throughout the entire set of variables.

#define i ilegal
            do j=1,ibuff(1)
              ierr=nf_inq_dim (ncid(node),j,string,dimsize(j,node))
              if (ierr .eq. nf_noerr) then
                lstr=lenstr(string)
                if (node.eq.0) then
                  ldim(j)=lstr
                  dimname(j)=string(1:lstr)
                elseif (lstr.ne.ldim(j) .or. string(1:lstr).ne.
     &                                 dimname(j)(1:ldim(j)) )then
                  write(*,'(/1x,2A,I3,3A/12x,6A/12x,3A/)')    '### ',
     &                 'ERROR: Name of dimension #', j, ', named ''',
     &                  string(1:lstr),   ''' in netCDF file',  '''',
     &                  ncname(node)(1:lncn),   ''' does not match ',
     &                 'name ''',  dimname(j)(1:ldim(j)), ''' with ',
     &                 'the corresponding name from netCDF file ''',
     &                                  ncname(0)(1:lncn),    '''.'
                  mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
                endif                   ! error not based on ierr so using mayday.
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &            'Cannot determine name and size of dimension #', j,
     &            'in netCDF file ''',   ncname(node)(1:lncn), '''.',
     &                                             nf_strerror(ierr)
                mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
              endif
            enddo
#undef i


      end subroutine check_ndims_nvars_natts  !]

! ----------------------------------------------------------------------
      subroutine create_catalog_of_var_names_IDs_ranks( node )  ![

      ![Create catalog of variable names, IDs, and ranks (number of
      ! dimensions) throughout the entire set. The meaning of the arrays
      ! defined here is as follows:
      !
      !     nvars    -- total number of variables discovered;
      !     vname(i), where i=1:nvars -- variable name;
      !     vid(i,node) -- netCDF ID for that variable in netCDF file for
      !                      MPI node "node", where node=0:nnodes-1;
      !     vdims(i) -- "rank", i.e. namber of dimensions of that variable;
      !     vnode(i) -- the file index for the lowest MPI node where the
      !                 variable has been found, or has achieved first time
      !                 its full rank (this is needed to differentiate
      !                 between the true variable and its proxy "dummy
      !                 scalar", if partit creates one, e.g., in the case
      !                 of boundary forcing variable in inner MPI node,
      !]                 where it is not needed).

      implicit none

      integer, intent(in) :: node

            if (node.eq.0) nvars=0
            do i=1,ibuff(2)
              ierr=nf_inq_varname (ncid(node), i, string)
              if (ierr .eq. nf_noerr) then
                lstr=lenstr(string)
                ierr=nf_inq_varndims (ncid(node), i, k)
                if (ierr .eq. nf_noerr) then
                  lnewvar=.true.
                  do j=1,nvars
                    lvar=lenstr(vname(j))
                    if (lstr.eq.lvar .and. string(1:lstr)
     &                               .eq.vname(j)(1:lvar)) then
                      lnewvar=.false.
                      vid(j,node)=i
                      if (k.gt.vdims(j)) then
                        vdims(j)=k
                        vnode(j)=node
                      endif
                    endif
                  enddo
                  if (lnewvar) then
                    nvars=nvars+1
                    vname(nvars)=string(1:lstr)
                    vid(nvars,node)=i
                    vnode(nvars)=node
                    vdims(nvars)=k
                  endif
                else
                  write(*,'(/1x,3A,I3/12x,5A/12x,A)')  '### ERROR: ',
     &                      'Cannot determine number of dimensions ',
     &                      'for variable with id =', i,  'named ''',
     &                       string(1:lstr),  ''' in netCDF file ''',
     &                       ncname(node)(1:lncn), '''.',
     &                                             nf_strerror(ierr)
                  mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &                 'Cannot determine name of variable with id =',
     &                  i, 'in netCDF file ''', ncname(node)(1:lncn),
     &                                      '''.', nf_strerror(ierr)
                mayday=.true.; return ! goto 97 ! can't have goto out of subroutine
              endif
            enddo

      end subroutine create_catalog_of_var_names_IDs_ranks  !]

! ----------------------------------------------------------------------
      subroutine determine_joined_file_dim_sizes  ![
      implicit none

![ Determine sizes of dimensions for combined file: For partitionable
! dimensions 'xi_rho', 'xi_u', 'eta_rho' and 'eta_v' determine the
! extent of the physical grid in each direction as the maximum over all
! subdomains of the dimension of each partial file combined with its
! starting index, "xi_start" or "eta_start".   This is straghtforward
! for RHO-points, but for U- and V-dimensions it requires to take into
! account the fact that the subdomains adjacent to eastern or southern
! edge have one point less than the corresponding RHO-dimension.
! Consequently, all subsequent subdomains receive one-point shift.
! For all other dimensions, verify that the sizes are the same for
! all nodes.  Also find size of unlimited dimension, if it is present.
! Note that variable "tsize" is set to its default value 1 (meaning
! one record), and may or may not be overwritten by the actual size of
! unlimited dimension (if it exists). If the unlimited dimension does
! not exist, it retains its value of 1 so that the loop over records
! is still executed, but only once. !]

        id_xi_rho=0
        id_eta_rho=0
        id_xi_u =0
        id_eta_v=0

        XI_rho=0
        ETA_rho=0

        size_XI=1
        size_ETA=1
        size_S=1
        tsize=1

        do i=1,ndims
          dimsize(i,nnodes)=0
          lvar=lenstr(dimname(i))
          if (lvar.eq.6 .and. dimname(i)(1:lvar).eq.'xi_rho') then
            id_xi_rho=i
            do node=0,nnodes-1
              dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &               dimsize(i,node) +xi_start(node)-1 )
              size_XI=max(size_XI,dimsize(i,node))
              XI_rho=max(XI_rho, dimsize(i,nnodes))
            enddo

          elseif (lvar.eq.4 .and.dimname(i)(1:lvar).eq.'xi_u') then
            id_xi_u=i
            do node=0,nnodes-1
              if (xi_start(node).gt.1) then
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                 dimsize(i,node) +xi_start(node)-2 )
              else
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                                  dimsize(i,node) )
              endif
              size_XI=max(size_XI,dimsize(i,node))
              XI_rho=max(XI_rho, dimsize(i,nnodes)+1)
            enddo

          elseif (lvar.eq.7.and. dimname(i)(1:lvar).eq.'eta_rho') then
            id_eta_rho=i
            do node=0,nnodes-1
              dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &              dimsize(i,node) +eta_start(node)-1 )
              size_ETA=max(size_ETA,dimsize(i,node))
              ETA_rho=max(ETA_rho, dimsize(i,nnodes))
            enddo

          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_v') then
            id_eta_v=i
            do node=0,nnodes-1
              if (eta_start(node).gt.1) then
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                 dimsize(i,node) +eta_start(node)-2 )
              else
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                                   dimsize(i,node))
              endif
              size_ETA=max(size_ETA,dimsize(i,node))
              ETA_rho=max(ETA_rho, dimsize(i,nnodes)+1)
            enddo

          else
            dimsize(i,nnodes)=dimsize(i,0)
            do node=1,nnodes-1
              if (dimsize(i,0).ne.dimsize(i,node)) then
                lncn=lenstr(ncname(node))
                write(*,'(/1x,A,I3,3A,I4,1x,A/12x,4A/12x,3A,I4,A/)')
     &                 '### ERROR: Nonpartitionable dimension #',  i,
     &                 ' named ''', dimname(i)(1:lvar), ''', size =',
     &                  dimsize(i,node),   'in netCDF',    'file ''',
     &                  ncname(node)(1:lncn),    ''' has different ',
     &                 'size than the corresponding',
     &                 'dimension from file ''',   ncname(0)(1:lncn),
     &                     ''', which has size =', dimsize(i,0), '.'
                mayday=.true.; return ! goto 97
              endif
            enddo
            if (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'s_rho') then
              size_S=max(size_S, dimsize(i,0))
            elseif (lvar.eq.3.and.dimname(i)(1:lvar).eq.'s_w') then
              size_S=max(size_S, dimsize(i,0))
            endif
          endif
          if (i.eq. unlimdimid) then
            tsize=dimsize(i,nnodes)
            dimsize(i,nnodes)=nf_unlimited
          endif
        enddo ! <-- i loop over dimensions

      end subroutine determine_joined_file_dim_sizes  !]

! ----------------------------------------------------------------------
      subroutine indentify_boundary_edges  ![
      implicit none

#ifdef VERBOSE
        write(*,'(1x,A)') 'Identifying presense of boundary edges:'
#endif
        do node=0,nnodes-1
          western_edge(node)=.true.
          eastern_edge(node)=.true.
          southern_edge(node)=.true.
          northern_edge(node)=.true.

          if (xi_start(node).gt.1) then
            western_edge(node)=.false.
          endif
          if (id_xi_rho.gt.0) then
            if ( xi_start(node)+dimsize(id_xi_rho,node)
     &          .lt.XI_rho ) eastern_edge(node)=.false.
          endif
          if (id_xi_u.gt.0) then
            if ( xi_start(node)+dimsize(id_xi_u,node)
     &          .lt.XI_rho ) eastern_edge(node)=.false.
          endif
          if (eta_start(node).gt.1) then
            southern_edge(node)=.false.
          endif
          if (id_eta_rho.gt.0) then
            if ( eta_start(node)+dimsize(id_eta_rho,node)
     &           .lt.ETA_rho ) northern_edge(node)=.false.
          endif
          if (id_eta_v.gt.0) then
            if ( eta_start(node)+dimsize(id_eta_v,node)
     &          .lt.ETA_rho ) northern_edge(node)=.false.
          endif

#ifdef VERBOSE
          if(my_rank==0) then
            if (node.eq.0) then
              write(*,'(8x,A,I4,4(2x,A,L1))') 'node =', node,
     &       'WST=', western_edge(node),   'EST=', eastern_edge(node),
     &       'SOU=', southern_edge(node),  'NOR=', northern_edge(node)
            else
                write(*,'(14x,I4,4(6x,L1))') node,
     &                       western_edge(node),   eastern_edge(node),
     &                       southern_edge(node),  northern_edge(node)
            endif
          endif
#endif
        enddo ! <- node=0,nnodes-1

      end subroutine indentify_boundary_edges  !]

! ----------------------------------------------------------------------
      subroutine create_joined_chunked_file  ![
      implicit none                          !  create joined file with chunking for compression

      ! Create combined netCDF file:   Once the completeness of the set of
      !------- -------- ------ -----   partial files have been established
      ! and dimensions survive consistency check, create the combined file,
      ! define its dimensions and copy global attributes.

      call create_joined_empty_file
      if(mayday) return                ! return if error occured in create_joined_empty_file

! Define dimensions: also compute the size of buffer needed to
! accommodate the largest array.

# ifdef VERBOSE
        write(*,'(/1x,A,3x,A,1x,A,1x,A)')
     &                 'Dimensions:', 'id', 'size', 'name'
# endif
        size_XI=1
        size_ETA=1
        size_S=1

        do i=1,ndims
          lvar=lenstr(dimname(i))
          ! Set time dimension to actual length, not unlimited.
          ! This is necessary for netcdf independent parallel I/O (can't have unlimited dim).
          if(dimsize(i,nnodes) /= 0) then ! Not time dimension
            ierr=nf_def_dim (nctarg, dimname(i)(1:lvar),
     &                       dimsize(i,nnodes), dimid(i))
          else                            ! is  time dimension
            ierr=nf_def_dim (nctarg, dimname(i)(1:lvar),
     &                       tsize, dimid(i))
          endif
          if (ierr .eq. nf_noerr) then
            if (dimid(i) .eq. i) then
# ifdef VERBOSE
              write(*,'(14x,I3,I5,1x,3A)') dimid(i),
     &                               dimsize(i,nnodes), '''',
     &                               dimname(i)(1:lvar), ''''
# endif
              if (dimname(i)(1:3) .eq. 'xi_') then
                size_XI=max(size_XI, dimsize(i,nnodes))
              elseif (dimname(i)(1:4) .eq. 'eta_') then
                size_ETA=max(size_ETA, dimsize(i,nnodes))
              elseif (dimname(i)(1:5) .eq. 's_rho' .or.
     &                dimname(i)(1:3) .eq. 's_w') then
                size_S=max(size_S, dimsize(i,nnodes))
              endif
            else
              write(*,'(/1x,2A,I3,1x,5A/12x,2A,I3,A/)')  '### ERROR: ',
     &        'id =', dimid(i), 'for dimension ''', dimname(i)(1:lvar),
     &        ''' from netCDF file ''',    nctargname(1:ltrg),    '''',
     &                    'differs from ', 'the original id =', i, '.'
              mayday=.true.; return ! goto 97
            endif
          else
           write(*,'(/1x,4A/12x,A/)')      '### ERROR: Cannot define ',
     &    'dimension ''', dimname(i)(1:lvar), '''.', nf_strerror(ierr)
              mayday=.true.; return ! goto 97
          endif
        enddo ! <-- i loop over dimensions

        call copy_global_attributes
        if(mayday) return

! Define variables and copy their attributes.

#ifdef VERBOSE
        if(my_rank==0) write(*,'(1x,2A)') 'Variables, their dimensions and ',
     &                                           'attributes:'
#endif

        if(my_rank==0) write(*,'(/1x,A)')
     &                  'Calculating chunk sizes for each variable.'

        ! Need to assign read partitions to MPI procs:
        ! Use roms subdomains method    |  ...  |
        ! to divide full domain:        | 3 4 5 |
        !                               | 0 1 2 |
        !
        ! Note this division rounds down but i/j node uses zero counting.
        ! This is also for all vars!
        jproc=my_rank/nprocs_x            ! indices iproc,jproc identifying
        iproc=my_rank-jproc*nprocs_x      ! the location of current subdomain

        if (allocated(chunksize_x)) then  ! If several sets of partial files then
          deallocate(chunksize_x)         ! need to deallocate for next set else seg fault.
          deallocate(chunksize_y)
          deallocate(n_chunks_x)
          deallocate(n_chunks_y)
          deallocate(p_chunks_x)
          deallocate(p_chunks_y)
        endif
        allocate(chunksize_x(nvars))      ! Allocate chunking variables per variable.
        allocate(chunksize_y(nvars))      ! Note this should probably be done per variable
        allocate(n_chunks_x(nvars))       ! type (r/u/v) and not for every variable.
        allocate(n_chunks_y(nvars))
        allocate(p_chunks_x(nvars))
        allocate(p_chunks_y(nvars))

        do i=1,nvars ! LOOP THROUGH VARS

#ifdef MPIVERBOSE
        if(my_rank==0) write(*,'(/1x,2A)') 'Defining var:', vname(i)
#endif

          node=vnode(i)
          lncn=lenstr(ncname(node))
          if (ncid(node).eq.-1) ierr=nf_open  (ncname(node),
     &                               nf_nowrite, ncid(node))
          if (ierr .eq. nf_noerr) then
            ierr=nf_inq_var (ncid(node), vid(i,node), vname(i),
     &                vartype(i), vdims(i), dimids(1,i),  varatts)
            ! DevinD: set variable as a collective I/O variable: ! Now done by master only but leave here for record of independent mode!
!            ierr=nf90_var_par_access(ncid(node), vid(i,node), NF90_COLLECTIVE)
!        print *, 'read collective variable set, ierr = ', nf90_strerror(ierr)

            if (ierr .eq. nf_noerr) then
              lvar=lenstr(vname(i))

              if (.not. vdims(i)<=1 .and.                   ! Not scalar, 1D vars, so calc chunk sizing
     &            dimsize(dimids(2,i),node) /= tsize) then  ! Not time_step var (rank 2 array)

                do d=1,vdims(i) ! Loop through var dimensions to adjust for rho/u/v points
                                ! u/v points are 1 shorter than rho points.

                  k=dimids(d,i) ! Value of variables current dimension ID
                                ! Use dimid to identify x/y/z rho/u/v point
                  if     (k.eq.id_xi_rho)  then      ! xi_rho   points in x of full domain.
                    n_x = dimsize(dimids(d,i),nnodes) ! For dimsize array x is last dimension hence vdims(i)
                  elseif (k.eq.id_xi_u)    then      ! xi_u    points in x of full domain.
                    n_x = dimsize(dimids(d,i),nnodes)
                  elseif (k.eq.id_eta_rho) then      ! eta_rho points in y of full domain.
                    n_y = dimsize(dimids(d,i),nnodes)
                  elseif (k.eq.id_eta_v)   then      ! eta_v   points in y of full domain.
                    n_y = dimsize(dimids(d,i),nnodes)
                  endif

                  if (d==4) then ! 4D variable (x,y,z,t), thus has z dim.
                    n_z = dimsize(dimids(3,i),nnodes) ! Num points in z of full domain.
                  else
                    n_z = 1 ! 2D or 1D variable
                  endif

                end do ! <- d=1,vdims(v)

                ! Chunksizing
                ! - currently for 1 chunk per process. This potentially needs further work.
                chunksize_x(i) = (n_x+nprocs_x-1) / nprocs_x        ! rounding up of n_x / nprocs_x
                chunksize_y(i) = (n_y+nprocs_y-1) / nprocs_y

                ! Number of chunks in x and y
                n_chunks_x(i) = (n_x+chunksize_x(i)-1) / chunksize_x(i) ! rounding up of n_x / chunksize_x
                n_chunks_y(i) = (n_y+chunksize_y(i)-1) / chunksize_y(i)

                ! Number of chunks for each proc in x and y
                p_chunks_x = (n_chunks_x(i)+nprocs_x-1) / nprocs_x ! rounding up of n_chunks_x / nprocs_x
                p_chunks_y = (n_chunks_y(i)+nprocs_y-1) / nprocs_y

#ifdef MPIVERBOSE
                print *,'Full domain size of var:', vname(i)
                write(*,'(3(A,i5))')' n_x= ',n_x,' n_y= ',n_y,
     &                                             ' n_z= ',n_z
                print *,'chunk sizing of var:', vname(i)
                print *,'chunksize_x= ',chunksize_x(i),
     &                 ' chunksize_y= ',chunksize_y(i)
                print *,'n_chunks_x= ',n_chunks_x(i),
     &                                ' n_chunks_y= ',n_chunks_y(i)
                print *,'p_chunks_x= ',p_chunks_x(i),
     &                 ' p_chunks_y= ',p_chunks_y(i)
#endif

                ! Since ncjoin_mpi writes at each timestep, the chunk size in time
                ! dim needs to be 1, else the chunk would need to be decompressed before
                ! adding another timestep to that chunk. This is inefficient.
                ! Hence set time dim to 1. We don't bother with dividing the z dim.
                set_req_chunk_sizes=0 ! reset chunk sizes for new variable
                do d=1,vdims(i)       ! Loop through number of variable's dimensions

                  if(d == vdims(i) .and.
     &               dimsize(dimids(d,i),nnodes)==tsize) then ! Should be time (tiny risk of bug here if dim randomly = tsize...)
                    set_req_chunk_sizes(d)=1 ! time dim
                  elseif(d == 1) then        ! x dim, (not time dimension caught above)
                    set_req_chunk_sizes(d)=chunksize_x(i)
                  elseif(d == 2) then        ! y dim
                    set_req_chunk_sizes(d)=chunksize_y(i)
                  elseif(d == 3) then        ! z dim
                    set_req_chunk_sizes(d)=n_z
                  endif
                end do

                ! OLD F77 method:
                ! ierr=nf_def_var (nctarg, vname(i)(1:lvar),vartype(i), vdims(i), dimids(1,i), varid(i))

                ! NEW F90 method:
                ! - nf90 routine no longer requires number of dimensions for
                !   the variable (hence vdims(i) not required).
                ! - However, dimids needs to be passed as an array.
                ierr=nf90_def_var (nctarg, vname(i)(1:lvar),vartype(i),
     &                             dimids(1:vdims(i),i), varid(i),
     &                             chunksizes=set_req_chunk_sizes(1:vdims(i)),
     &                             deflate_level=deflate_level, shuffle=shuffle)

#ifdef MPIVERBOSE
                chunk_sizes=0 ! Refresh debug values
                ierr=nf90_inquire_variable(nctarg,varid(i),chunksizes=chunk_sizes)
                print *,'var: ',vname(i)(1:lvar),
     &                 ' set chunksizes: ',chunk_sizes
#endif

              else  ! is 1D data

                ierr=nf90_def_var (nctarg, vname(i)(1:lvar),vartype(i),
     &                             dimids(1:vdims(i),i), varid(i))      ! Doesn't need compression?

              end if ! <- if (vdims(i) /= 1)

              ! Non-master pre-processing left here for reference until collective/independent story resolved.
              ! Set variables for parllel I/O 'collective' mode, this apparently
              ! allows for besst optimization by underlying MPI_IO library.
!              ierr=nf90_var_par_access(nctarg, varid(i), NF90_COLLECTIVE)
!        print *, 'write collective variable, ierr = ', nf90_strerror(ierr)

              if (ierr .eq. nf_noerr) then
#ifdef VERBOSE
                write(*,'(8x,3A,8I3)')   '''', vname(i)(1:lvar),
     &              ''',  dimids =', (dimids(j,i), j=1,vdims(i))
#endif
                do j=1,varatts
                  ierr=nf_inq_attname (ncid(node), vid(i,node),
     &                                              j, string)
                  if (ierr .eq. nf_noerr) then
                    lstr=lenstr(string)
                    ierr=nf_copy_att (ncid(node), vid(i,node),
     &                       string(1:lstr), nctarg, varid(i))
                    if (ierr. ne. nf_noerr) then
                      write(*,'(/1x,2A,I3,3A/12x,4A)')   '### ERROR: ',
     &                 'Cannot copy attribute #', j,' for variable ''',
     &                  vname(i)(1:lvar),  ''' into netCDF', 'file ''',
     &                  nctargname(1:ltrg), '''.  ', nf_strerror(ierr)
                      mayday=.true.; return ! goto 97
                    endif
#ifdef VERBOSE
                    write(*,'(16x,3A)') '''', string(1:lstr), ''''
#endif
                  else
                    write(*,'(/1x,2A,I3/12x,3A/12x,A/)') '### ERROR: ',
     &                             'Cannot get name of attribute #', j,
     &                      'for variable ''', vname(i)(1:lvar), '''.',
     &                                               nf90_strerror(ierr)
                    mayday=.true.; return ! goto 97
                  endif
                enddo
              else
                write(*,'(/8x,5A/)') 'ERROR: Cannot define ',
     &                  'variable ''', vname(i)(1:lvar), '''.',
     &                   nf90_strerror(ierr)
                mayday=.true.; return ! goto 97
              endif
            else

            write(*,'(/8x,2A/15x,A,I3,1x,3A/)')  '### ERROR: Cannot ',
     &        'determine name, type and attributes for variable #', i,
     &            'from netCDF file ''', ncname(node)(1:lncn),  '''.'
              mayday=.true.; return ! goto 97
            endif

! Determine whether partitionable dimensions or unlimited dimension
! are present for this variable: the convention adopted here is:
!           part_type = 0 -- non-partitionable array;
!                     = 1 -- has partitionable XI-dimension only;
!                     = 2 -- has partitionable ETA-dimension only;
!                     = 3 -- partitionable in both XI and ETA.

            series(i)=.false.
            part_type(i)=0
            do j=1,vdims(i)
              if (dimids(j,i).eq.id_xi_rho .or.
     &            dimids(j,i).eq.id_xi_u) then
                part_type(i)=part_type(i)+1
              elseif (dimids(j,i).eq.id_eta_rho .or.
     &                dimids(j,i).eq.id_eta_v) then
                part_type(i)=part_type(i)+2
              elseif (dimids(j,i).eq.unlimdimid) then
                series(i)=.true.
              endif
            enddo


            if (node.gt.0) then
              ierr=nf90_close(ncid(node)) ! DevinD: Most vars in node==0, hence keep
              ncid(node)=-1               ! it open whilst looping through vars, close
            endif                         ! any node>0 as might not be used again.

          else
            write(*,'(/1x,A,1x,3A/12x,A)')  '### ERROR: Cannot open ',
     &                  'netCDF file ''', ncname(node)(1:lncn), '''.',
     &                                           nf90_strerror(ierr)
            mayday=.true.; return ! goto 97
          endif

        enddo  ! <-- i=1,nvars, variable IDs.

        ierr=nf90_close(ncid(0)) ! DevinD: completed master only reads, close remaining
        ncid(0)=-1               ! open read file.

! Leave definition mode

        ierr=nf_enddef (nctarg)
        ierr=nf90_close(nctarg)  ! DevinD: close output file for parallel open later.

# ifdef VERBOSE
        if(my_rank==0) write(*,'(/1x,A)') 'Leaving definition mode.'
# endif

      end subroutine create_joined_chunked_file  !]

! ----------------------------------------------------------------------
      subroutine create_joined_empty_file  ![
      implicit none

        i=lenstr(root_bak)
        if (sffx_bak(1:1).ne.' ') then
          j=lenstr(sffx_bak)
          if (root_bak(i:i).eq.'.' .and. sffx_bak(1:1).eq.'.') then
            nctargname=root_bak(1:i)/ /sffx_bak(2:j)
          else
            nctargname=root_bak(1:i)/ /sffx_bak(1:j)
          endif
        else
          nctargname=root_bak(1:i)
        endif
        ltrg=lenstr(nctargname)

        j=0
        do i=1,ltrg
          if (nctargname(i:i).eq.'/') j=i+1
        enddo
        if (j.gt.0) then
          nctargname=nctargname(j:ltrg)
          ltrg=ltrg-j+1
        endif

!--> create...

!        ierr=nf90_create (nctargname(1:ltrg), IOR(NF90_NETCDF4, NF90_MPIIO),       ! for not master only version that worked with independent mode!
!     &                    nctarg, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
        ierr=nf90_create (nctargname(1:ltrg), NF90_NETCDF4, nctarg)                 ! Master only works!
!        ierr=nf90_create (nctargname(1:ltrg), IOR(NF90_NETCDF4, NF90_MPIIO), nctarg) ! Master only didn't work (hdf error at write)
        if (ierr .eq. nf_noerr) then
          if(my_rank==0) write(*,'(/1x,3A)')  'Created netCDF file ''',
     &                        nctargname(1:ltrg), '''.'
          if(my_rank==0) write(*,'(1x,A,I1)')
     &      '-> compression level (deflate_level) = ', deflate_level
          if(my_rank==0) write(*,*)
     &      '-> compression shuffle = ', shuffle
        else
          write(*,'(/1x,4A/12x,A/)')     '### ERROR: Cannot create ',
     &                          'netCDF file ''', nctargname(1:ltrg),
     &                                      '''.', nf_strerror(ierr)
          mayday=.true.; return ! goto 97
        endif

      end subroutine create_joined_empty_file

! ----------------------------------------------------------------------
      subroutine copy_global_attributes  ![
      implicit none

! Copy all global attributes, except 'partition'.

# ifdef VERBOSE
        write(*,'(1x,A)') 'Copying global attributes:'
# endif
        lncn=lenstr(ncname(0))
        if (ncid(0).eq.-1) ierr=nf_open (ncname(0), nf_nowrite,
     &                                                 ncid(0))
        if (ierr .eq. nf_noerr) then

          do i=1,ngatts
            ierr=nf_inq_attname (ncid(0), nf_global, i, string)
            if (ierr. eq. nf_noerr) then
              lvar=lenstr(string)
              if (string(1:lvar) .ne. 'partition') then
                ierr=nf_copy_att (ncid(0), nf_global, string(1:lvar),
     &                                             nctarg, nf_global)
                if (ierr .ne. nf_noerr) then
                  write(*,'(/1x,4A/12x,3A/12x,A)')     '### ERROR: ',
     &             'Cannot copy global attribute ''', string(1:lvar),
     &             ''' into netCDF',  'file ''',  nctargname(1:ltrg),
     &                                     '''.',  nf_strerror(ierr)
                  mayday=.true.; return ! goto 97
                endif
# ifdef VERBOSE
                if(my_rank==0) write(*,'(20x,3A)') '''', string(1:lvar), ''''
# endif
              endif
            else
            write(*,'(/1x,2A,I3/12x,3A/12x,A/)') '### ERROR: Cannot',
     &                    ' determine name of global attribute #', i,
     &                     'from netCDF file ''',  ncname(0)(1:lncn),
     &                                     '''.',  nf_strerror(ierr)
              mayday=.true.; return ! goto 97
            endif
          enddo
        else
          write(*,'(/1x,A,1x,3A/14x,A)')     '### ERROR: Cannot open ',
     &   'netCDF file ''', ncname(0)(1:lncn), '''.', nf_strerror(ierr)
          mayday=.true.; return ! goto 97
        endif

      end subroutine copy_global_attributes  !]

! ----------------------------------------------------------------------
      subroutine exchange_mpi_master_pre_processing_to_all  ![
                                                            !  exchange all info about file pre-processing
      implicit none                                         !  to all other mpi procs

      call MPI_BCAST(nvars,  1, MPI_INT, 0, MPI_COMM_WORLD, ierr) ! Once working group scalars together and send in
      call MPI_BCAST(nnodes, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr) ! 1 long array to spare latency of repeated calls.
      if(my_rank /= 0) then                                       ! Even with 120 cores didn't make a difference...
        jproc=my_rank/nprocs_x        ! indices iproc,jproc identifying
        iproc=my_rank-jproc*nprocs_x  ! the location of current mpi proc subdomain
        if (allocated(chunksize_x)) then
          deallocate(chunksize_x)
          deallocate(chunksize_y)
          deallocate(p_chunks_x)
          deallocate(p_chunks_y)
        endif
        allocate(chunksize_x(nvars))
        allocate(chunksize_y(nvars))  ! We previously only allocated all these arrays
        allocate(p_chunks_x(nvars))   ! for the master proc. Hence need to do again
        allocate(p_chunks_y(nvars))   ! for all other procs.

        if (allocated(ncid)) then
          deallocate(dimsize)
          deallocate(vid)
          deallocate(northern_edge)
          deallocate(southern_edge)
          deallocate(eastern_edge)
          deallocate(western_edge)
          deallocate(eta_start)
          deallocate(xi_start)
          deallocate(ncname)
          deallocate(ncid)
        endif
        allocate (ncid(0:nnodes-1))
        allocate (ncname(0:nnodes-1))
        allocate (xi_start(0:nnodes-1))
        allocate (eta_start(0:nnodes-1))
        allocate (western_edge(0:nnodes-1))
        allocate (eastern_edge(0:nnodes-1))
        allocate (southern_edge(0:nnodes-1))
        allocate (northern_edge(0:nnodes-1))
        allocate (vid(maxvars,0:nnodes-1))
        allocate (dimsize(maxdims,0:nnodes))
      endif  ! <- my_rank /= 0

      call MPI_BCAST(tsize,           1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(part_type, maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(chunksize_x,     nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(chunksize_y,     nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(p_chunks_x,  nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(p_chunks_y,  nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'

      call MPI_BCAST(dimsize, maxdims*(nnodes+1), MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(dimids,     maxdims*maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(varid,              maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(vid,         maxvars*nnodes, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(vdims,              maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'

      call MPI_BCAST(vname, maxvars*32, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(nctargname,    64, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(ncname, nnodes*64, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'

      call MPI_BCAST(western_edge,  nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(southern_edge, nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(eastern_edge,  nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(northern_edge, nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(series,       maxvars, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'

      call MPI_BCAST(vartype,  maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(xi_start,  nnodes, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(eta_start, nnodes, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(id_xi_u,        1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(id_xi_rho,      1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(id_eta_v,       1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(id_eta_rho,     1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'
      call MPI_BCAST(unlimdimid,     1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=0) print *, 'POOR MPI EXCHANGE'

#ifdef TIMING
      nclk=3-nclk
      call system_clock (iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      net_fcrt_clk=net_fcrt_clk+inc_clk  ! Timing for file creation
#endif

      end subroutine  !]

! ----------------------------------------------------------------------
      subroutine collectively_open_input_output_files  ![
      implicit none

        ierr=nf90_open (nctargname,                               ! MPI I/O: open output file
!     &                  IOR(NF90_NETCDF4, NF90_MPIIO), nctarg,
     &                  IOR(NF90_WRITE, NF90_MPIIO), nctarg,
     &                  comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
        if(ierr/=0) print *, 'BAD OPEN OF OUTFILE'
        ! Set variables for parllel I/O 'collective' mode, this apparently
        ! allows for best optimization by underlying MPI_IO library.
        ! Would like to try in default (independent) mode, but couldn't get working
        ! in master only pre-processing version of ncjoin_mpi.
        if(my_rank==0) write(*,'(1x,A)')
     &                'Set ''collective'' write of variables (MPI-I/O).'
        do i=1,nvars
          if (part_type(i)/=0) then  ! catch to enable master only of non-partitioned arrays
          ierr=nf90_var_par_access(nctarg, varid(i), NF90_COLLECTIVE)  ! Otherwise hang for master only waiting for other procs
#ifdef MPIVERBOSE
            if(my_rank==0)
     &        print *, 'nf90_var_par_access = collective: ', vname(i) ! debug output
#endif
          endif
!          ierr=nf90_var_par_access(nctarg, varid(i), NF90_INDEPENDENT) ! Doesn't work?
!          ierr=nf90_inquire_variable(nctarg, varid(i), var_name_debug)
!          if(ierr/=0) print *, 'BAD INQUIRE VAR'
!          if(my_rank==3) print *, 'var name from varid: ', var_name_debug
        enddo

        do node=0,nnodes-1 ! MPI I/O: open partition input files
          ierr=nf90_open (ncname(node),
!     &                    IOR(NF90_NETCDF4, NF90_MPIIO), ncid(node), ! not sure if I need nf_nowrite in here somehow
     &                    IOR(NF90_NOWRITE, NF90_MPIIO), ncid(node),  ! nowrite is faster for read only.
     &                    comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
          if(ierr/=0) print *, 'BAD OPEN OF INFILE'
!          ierr=nf90_inquire_variable(ncid(node), vid(1,node), var_name_debug)
!          if(ierr/=0) print *, 'BAD INQUIRE VAR'
!          if(my_rank==3) print *, 'var name from vid: ', var_name_debug

          ! Note, all read input variables currently in 'independent' mode (default).
          ! This seemed to give good scaling of the read portion of code.
          ! To set variables to 'collective' mode, would have to loop through variables
          ! here.
        enddo

      end subroutine collectively_open_input_output_files  !]

! ----------------------------------------------------------------------
      subroutine assign_read_in_partial_files_to_mpi_proc  ![
                                                           !  Assign the input files to mpi proc
      implicit none                                        !  if it fall within proc's range

          if (tsize.gt.1) then
#ifdef TIMING
            nclk=3-nclk
            call system_clock (iclk(nclk), clk_rate,clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            net_gray_clk=net_gray_clk+inc_clk

            if(my_rank==0) write(*,'(/F8.1,1x,A,I8,1x,A,I8,1x,A)')
     &         dble(iclk(nclk)-iclk_init)/dble(clk_rate),
#else
            if(my_rank==0) write(*,'(4x,A,I8,1x,A,I8,1x,A)')
#endif
     &     'Processing variable', i, 'out of', nvars,  '...'
          endif

          if (part_type(i).gt.0) then           ! Assign nodes to proc once per variable
                                                ! Only necessary for partitioned variable (2D/3D vars)

            start_out=1                             ! init array, z dimension always starts at 1 as not subdivided.
            p_points_x=p_chunks_x(i)*chunksize_x(i) ! sub-domain covers (chunk_size_x X processor_chunks_in_x)
            p_points_y=p_chunks_y(i)*chunksize_y(i)
            start_out(1)=iproc*p_points_x+1     ! Start of sub-domain in x
            p_str_out_x=start_out(1)
                                                ! is not full domain, it will differ.
            p_end_out_x=(iproc+1)*p_points_x
            start_out(2)=jproc*p_points_y+1     ! Start of sub-domain in y
            p_str_out_y=start_out(2)
            p_end_out_y=(jproc+1)*p_points_y

            if (vdims(i)==4) then               ! 4D variable (x,y,z,t), thus has z dim.
              n_z = dimsize(dimids(3,i),nnodes) ! Num points in z of full domain.
            else
              n_z = 1                           ! 2D or 1D variable
            endif

            count_out=(/p_points_x,p_points_y,n_z,1/) ! If 2D var then n_z is set to 1 so this works

            ! catch N & E boundaries with partially full chunks:
            if(iproc==(nprocs_x-1)) then        ! -1 as nodes 0 counting
              n_x=dimsize(dimids(1,i),nnodes)
              p_end_out_x=n_x
              p_points_x=n_x-(iproc*p_points_x) ! Total x points - start of eastern sub-domain
              count_out=(/p_points_x,p_points_y,n_z,1/)
            end if
            if(jproc==(nprocs_y-1)) then        ! -1 as nodes 0 counting
              n_y=dimsize(dimids(2,i),nnodes)
              p_end_out_y=n_y
              p_points_y=n_y-(jproc*p_points_y) ! Total y points - start of northern sub-domain
              count_out=(/p_points_x,p_points_y,n_z,1/)
            end if

#ifdef MPIVERBOSE
            if(my_rank==0) print *,'  ' ![
            if(my_rank==0) print *,'p_points_x= ',p_points_x,' p_points_y= ',p_points_y
            do p_indx=0,n_procs-1
              if(my_rank==p_indx) then
                print *, '------ my_rank= ',my_rank
                print *, 'iproc= ', iproc,' jproc= ', jproc
                print *, 'p_str_out_x= ',p_str_out_x
                print *, 'p_end_out_x= ',p_end_out_x
                print *, 'count_out_x=  ',count_out(1)
                print *, 'p_str_out_y= ',p_str_out_y
                print *, 'p_end_out_y= ',p_end_out_y
                print *, 'count_out_y=  ',count_out(2)
                print *, 'count_out_z=  ',count_out(3)
              end if
              call MPI_Barrier(MPI_COMM_WORLD, ierr)
            end do ! <- p_indx
            if(my_rank==0) print *,'  '
            call MPI_Barrier(MPI_COMM_WORLD, ierr) !]
#endif

            if (allocated(data)) deallocate(data)
            allocate(data(p_points_x,p_points_y,n_z)) ! Match proc ranges
#ifdef MPIVERBOSE
            if(my_rank==0) print *,'proc data var size x= ', p_points_x
            if(my_rank==0) print *,'proc data var size y= ', p_points_y
            if(my_rank==0) print *,'proc data var size z= ', n_z
            if(my_rank==0) print *,'  '
#endif

            ![A) ASSIGN READ IN PARTITIONS TO PROCS:
            ! Allocate list of nodes read into proc (lnodes_in_proc)
            ! Also store number of nodes into proc (nnodes_in_proc) for looping.
            ! (i-index fastest varying (to be consistent with Sasha's vars)
            ! If using only 1 proc, then that proc will read in all nodes (nnodes).
            ! However, each mpi proc has its own copy of all arrays, so if we are using many
            ! procs, we don't want each proc carrying arrays sizes of all nodes (nnodes),
            ! since they will not be reading in from every node.
            ! This will need a smarter algorithm to chose a reasonable array size,
            ! for now do some crude hacks.
            ! If using few cores, then ok to use nnodes. If not then use smaller value. !]
            if (n_procs<=8) then
              guess_nnodes=nnodes
            else if (n_procs<=20) then
              guess_nnodes=nnodes/4    ! hybrid number > guess_nnodes but < nnodes
            endif
#ifdef MPIVERBOSE
            if(my_rank==0) print *,'No. of local MPI arrays copies= ', n_procs
            if(my_rank==0) print *,'Estimated read nodes per proc= ', guess_nnodes
            if(my_rank==0) print *,'  '
#endif

            if (.not. allocated(lnodes_in_proc))
     &        allocate(lnodes_in_proc(guess_nnodes)) ! All nodes if using single proc, as all partitions relevant to 1 proc...
            lnodes_in_proc=-1                        ! Set to safe value
            if (.not. allocated(start_in))
     &        allocate(start_in(4,guess_nnodes))     ! In theory for a new set of partial files
            if (.not. allocated(count_in))           ! these arrays should be deallocated and
     &        allocate(count_in(4,guess_nnodes))     ! reallocated. But since the size if based on
            start_in(1:2,:)=-1                       ! fixed 'guess_nnodes', should be fine to leave.
            count_in(1:2,:)=-1  ! set x & y to safe (wrong) value
            start_in(3,:)=1     ! set z
            count_in(3,:)=n_z   ! if 2D variable n_z is set to 1
            start_in(4,:)=1     ! set t
            count_in(4,:)=1

            if (.not. allocated(data_x_in_str))
     &        allocate(data_x_in_str(guess_nnodes))  ! All nodes if using single proc, as all partitions relevant to 1 proc...
            if (.not. allocated(data_x_in_end))
     &        allocate(data_x_in_end(guess_nnodes))
            if (.not. allocated(data_y_in_str))
     &        allocate(data_y_in_str(guess_nnodes))
            if (.not. allocated(data_y_in_end))
     &        allocate(data_y_in_end(guess_nnodes))
            data_x_in_str=-1    ! Set to safe (wrong) value
            data_x_in_end=-1
            data_y_in_str=-1
            data_y_in_end=-1

            p_in_nodes=0        ! restart count
            do node=0,nnodes-1  ! check all nodes to see if within proc's range
              sw_corner=.false. ! reset corner flags
              se_corner=.false.
              nw_corner=.false.
              ne_corner=.false.

              x_in_dimsize=dimsize(dimids(1,i),node)  ! 2D/3D var 1 is x, 2 is y.
              y_in_dimsize=dimsize(dimids(2,i),node)

              if     (dimids(1,i).eq.id_xi_u)   then  ! Adjust node's x start position. (dimids indices: (maxdims,maxvars))
                x_start_in=max(xi_start(node)-1,1)    ! u-point
              elseif (dimids(1,i).eq.id_xi_rho) then
                x_start_in=xi_start(node)             ! rho-point
              end if
              if     (dimids(2,i).eq.id_eta_v)   then ! Adjust node's y start position.
                y_start_in=max(eta_start(node)-1,1)   ! v-point
              elseif (dimids(2,i).eq.id_eta_rho) then
                y_start_in=eta_start(node)            ! rho-point
              end if

              x_nc_sw=x_start_in                 ! input node's i SW corner WRT global indexing
              y_nc_sw=y_start_in                 ! input node's j SW corner WRT global indexing
              x_nc_se=x_start_in+x_in_dimsize-1  ! SE corner
              y_nc_se=y_start_in
              x_nc_nw=x_start_in                 ! NW corner
              y_nc_nw=y_start_in+y_in_dimsize-1
              x_nc_ne=x_start_in+x_in_dimsize-1  ! NE corner
              y_nc_ne=y_start_in+y_in_dimsize-1

              if( ! SW corner
     &           x_nc_sw >= p_str_out_x .and.  ! Check partition's SW corner
     &           x_nc_sw <= p_end_out_x .and.  ! falls within node.
     &           y_nc_sw >= p_str_out_y .and.
     &           y_nc_sw <= p_end_out_y) then
                sw_corner=.true.
              end if
              if( ! SE corner
     &           x_nc_se >= p_str_out_x .and.  ! Check partition's SE corner
     &           x_nc_se <= p_end_out_x .and.  ! falls within node.
     &           y_nc_se >= p_str_out_y .and.
     &           y_nc_se <= p_end_out_y) then
                se_corner=.true.
              end if
              if( ! NW corner
     &           x_nc_nw >= p_str_out_x .and.  ! Check partition's NW corner
     &           x_nc_nw <= p_end_out_x .and.  ! falls within node.
     &           y_nc_nw >= p_str_out_y .and.
     &           y_nc_nw <= p_end_out_y) then
                nw_corner=.true.
              end if
              if( ! NE corner
     &           x_nc_ne >= p_str_out_x .and.  ! Check partition's NE corner
     &           x_nc_ne <= p_end_out_x .and.  ! falls within node.
     &           y_nc_ne >= p_str_out_y .and.
     &           y_nc_ne <= p_end_out_y) then
                ne_corner=.true.
              end if

              ! Count node if any corners contained in proc range
              if(sw_corner .or. se_corner .or. nw_corner .or. ne_corner ) then
                p_in_nodes=p_in_nodes+1         ! Count nodes linked to proc
                lnodes_in_proc(p_in_nodes)=node ! Stored linked node number to proc

                ! Set start_in and count_in for node relative to proc:

                ! ALL WEST NODE FACES (x)
                ! Subtract west (start) corners of node x from proc x.
                ! If negative then read in from beginning of input file.
                start_in(1,p_in_nodes)=p_str_out_x-x_nc_sw+1
                count_in(1,p_in_nodes)=x_in_dimsize-(p_str_out_x-x_nc_sw) ! Less the indented start
                if(start_in(1,p_in_nodes)<1) then ! x_nc_sw is east of p_str_out_x
                  start_in(1,p_in_nodes)=1
                  count_in(1,p_in_nodes)=x_in_dimsize ! x size entire dim
                end if

                ! ALL EAST NODE FACES (x)
                ! catch nodes that go beyond proc range on east side
                if(p_end_out_x-x_nc_se<0) then
                  count_in(1,p_in_nodes)=x_in_dimsize-(x_nc_se-p_end_out_x)
                end if

                ! ALL SOUTH NODE FACES (y)
                ! Subtract south (start) corners of node y from proc y.
                ! If negative then read in from beginning of input file.
                start_in(2,p_in_nodes)=p_str_out_y-y_nc_sw+1
                count_in(2,p_in_nodes)=y_in_dimsize-(p_str_out_y-y_nc_sw) ! Less the indented start
                if(start_in(2,p_in_nodes)<1) then ! y_nc_sw is north of p_str_out_y
                  start_in(2,p_in_nodes)=1
                  count_in(2,p_in_nodes)=y_in_dimsize ! x size entire dim
                end if

                ! ALL NORTH NODE FACES (y)
                ! catch nodes that go beyond proc range on north side
                if(p_end_out_y-y_nc_ne<0) then
                  count_in(2,p_in_nodes)=y_in_dimsize-(y_nc_ne-p_end_out_y)
                end if

                ! Also, proc stores data in array 'data' with local (not global)
                ! indices, hence set inbound indices for proc's data array.
                ! Set data_x_in_str, data_x_in_end, & for y
                data_x_in_str(p_in_nodes)=max(x_nc_sw-p_str_out_x+1,1)          ! if node sw corner is west of proc sw corner, then use 1.
                data_x_in_end(p_in_nodes)=min(x_nc_se-p_str_out_x+1,p_points_x) ! if node se corner is east of proc se corner, then use end of data array in x. +1 for typical subtraction count issue.
                data_y_in_str(p_in_nodes)=max(y_nc_sw-p_str_out_y+1,1)          ! if node sw corner is west of proc sw corner, then use 1.
                data_y_in_end(p_in_nodes)=min(y_nc_ne-p_str_out_y+1,p_points_y) ! if node se corner is east of proc se corner, then use end of data array in x

              end if ! <- sw_corner .or. ... (any corner)

            end do   ! <- nodes
            nnodes_in_proc=p_in_nodes
            if (nnodes_in_proc > guess_nnodes) then                   ! ERROR: must recomplie ncjoin_mpi
              write(*,'(/1x,A,I3,2A,I3,A/)') 'ERROR:: guess_nnodes =',  !        with increased guess_nnodes
     &          guess_nnodes, ' is too low. Must be just >= ',
     &          'nnodes_in_proc = ', nnodes_in_proc,
     &          '. Increase guess_nnodes just enough, and recompile!'
              error stop 'ERROR: guess_nnodes<nnodes_in_proc. See log.'
            endif

#ifdef MPIVERBOSE
            do p_indx=0,n_procs-1 ![
              if(my_rank==p_indx) then
                if(my_rank==0) print *, 'ASSIGN READ NODES TO PROC:'
                print *, '      '
                print *, '------ my_rank= ',my_rank
                print *, 'nnodes_in_proc= ', nnodes_in_proc
                print *, 'lnodes_in_proc= ', lnodes_in_proc(1:4)
                do node =1,nnodes_in_proc
                  print *, '______'
                  print *, 'read node: ', lnodes_in_proc(node)
                  print *, 'start_in: ',start_in(:,node)
                  print *, 'count_in: ',count_in(:,node)
                  print *, 'data_x_in_str: ',data_x_in_str(node)
                  print *, 'data_x_in_end: ',data_x_in_end(node)
                  print *, 'data_y_in_str: ',data_y_in_str(node)
                  print *, 'data_y_in_end: ',data_y_in_end(node)
                end do
              end if
              call MPI_Barrier(MPI_COMM_WORLD, ierr)
            end do !] <- p_indx
#endif

          endif ! <- end assign partitions to mpi proc

      end subroutine assign_read_in_partial_files_to_mpi_proc  !]

! ----------------------------------------------------------------------
      subroutine read_write_scalar_var  ![
      implicit none

                lvar=lenstr(vname(i))
                if(my_rank==0) write(*,'(16x,3A)')
     &            'Copy scalar variable: ''',vname(i)(1:lvar),'''...'

                if (vartype(i) .eq. nf_char)       then
                  ierr=nf90_get_var(ncid(0), vid(i,0), char_in)
                elseif (vartype(i) .eq. nf_int)    then
                  ierr=nf90_get_var(ncid(0), vid(i,0), scalar_int)
                elseif (vartype(i) .eq. nf_float)  then
                  ierr=nf90_get_var(ncid(0), vid(i,0), scalar_real) ! Hoping double precision variable works fine here.
                elseif (vartype(i) .eq. nf_double) then             ! Compiler default for real is double.
                  ierr=nf90_get_var(ncid(0), vid(i,0), scalar_real)
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A/)') '### ERROR: scalar variable ',
     &              '''', vname(i)(1:lvar), ''' has unknown type.'
                  mayday=.true.; return ! goto 97
                endif
                if (ierr .eq. nf_noerr) then
                  if (vartype(i) .eq. nf_char)       then
                    ierr=nf90_put_var(nctarg ,varid(i),char_in)
                  elseif (vartype(i) .eq. nf_int)    then
                    ierr=nf90_put_var(nctarg, varid(i), scalar_int)
                  elseif (vartype(i) .eq. nf_float)  then
                    ierr=nf90_put_var(nctarg ,varid(i), scalar_real)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf90_put_var(nctarg, varid(i), scalar_real)
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/1x,4A/12x,3A/12x,A)')  '### ERROR: ',
     &                            'Cannot write scalar variable ''',
     &                           vname(i)(1:lvar), ''' into netCDF',
     &                                'file ''', nctargname(1:ltrg),
     &                                    '''.',  nf_strerror(ierr)
                    mayday=.true.; return ! goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot ',
     &                 'read scalar variable ''', vname(i)(1:lvar),
     &                                    '''.', nf_strerror(ierr)
                  mayday=.true.; return ! goto 97
                endif

      end subroutine read_write_scalar_var  !]

! ----------------------------------------------------------------------
      subroutine read_write_non_partitioned_var  ![
      implicit none

                lvar=lenstr(vname(i))
                if(my_rank==0 .and. rec==1) write(*,'(16x,3A)')
     &            'Copy non-partitioned array: ''', vname(i)(1:lvar), '''...'
                size=1
                do j=1,vdims(i)
                  if (dimids(j,i).eq.unlimdimid) then
                    start(j)=rec
                    count(j)=1
                  else
                    start(j)=1
                    count(j)=dimsize(dimids(j,i),0)
                  endif
                  size=size*count(j)
                enddo
                if (vartype(i) .eq. nf_char .or.
     &              vartype(i) .eq. nf_byte) then
                  size=size*1
                elseif (vartype(i) .eq. nf_short) then
                  size=size*2
                elseif (vartype(i) .eq. nf_int .or.
     &                  vartype(i) .eq. nf_float) then
                  size=size*4
                elseif (vartype(i) .eq. nf_double) then
                  size=size*8
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,3A/)')  '### ERROR: variable ''',
     &                  vname(i)(1:lvar), ''' has unknown type.'
                  mayday=.true.; return ! goto 97
                endif
!                if (size .gt. 8*max_buff_size) then ! DevinD uncommented quick fix
                  if (allocated(buff)) deallocate(buff)
                  max_buff_size=(size+7)/8
                  allocate(buff(max_buff_size))
#ifdef MPIVERBOSE
                  if(my_rank==0) write(*,*)
     &              '1D var: allocated "buff" with max_buff_size =', max_buff_size
#endif
!                endif

                if (vartype(i) .eq. nf_char) then

                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_vara_int    (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_float) then
                  ierr=nf_get_vara_real   (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_vara_double (ncid(0), vid(i,0),
     &                                    start,count, buff)
                endif
                if (ierr .eq. nf_noerr) then

                  if (vartype(i) .eq. nf_char) then

                  elseif (vartype(i) .eq. nf_int) then
                    ierr=nf_put_vara_int   (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_float) then
#ifdef MPIVERBOSE
                    print *, 'WRITING FLOAT non-par: ', 'my_rank', my_rank
#endif
                    ierr=nf_put_vara_real  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf_put_vara_double(nctarg, varid(i),
     &                                    start,count, buff)
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/8x,4A,I3/15x,3A/)')    '### ERROR: ',
     &                 'Cannot write variable ''', vname(i)(1:lvar),
     &              ''' for time record',rec, 'into netCDF file ''',
     &                  nctargname(1:ltrg),'''.', nf_strerror(ierr)
                    mayday=.true.; return ! goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A,I3,A/15x,A/)')     '### ERROR: ',
     &              'Cannot read variable ''',     vname(i)(1:lvar),
     &              ''' for time record',rec,'.', nf_strerror(ierr)
                  mayday=.true.; return ! goto 97
                endif

      end subroutine read_write_non_partitioned_var  !]

! ----------------------------------------------------------------------
      subroutine read_write_partitioned_var  ![
      implicit none

                lvar=lenstr(vname(i))
                if(my_rank==0 .and. rec==1) write(*,'(16x,2A,I3,1x,3A)')
     &           'Assembly partitioned ','array type', part_type(i),  'name = ''',
     &                                    vname(i)(1:lvar), ''''

!                do node=0,nnodes-1 ! Old ncjoin
                do in2p=1,nnodes_in_proc ! loop through only relevant nodes read into proc
                  node=lnodes_in_proc(in2p)

                  var_mask=.false.
                  if (part_type(i).eq.1 .and. lvar.gt.6) then
                    if (vname(i)(lvar-5:lvar).eq.'_south'
     &                        .and. southern_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')  'node =', node,
     &               'identified XI-partitioned southern ',
     &               'boundary array ''', vname(i)(1:lvar), ''''
#endif
                    elseif (vname(i)(lvar-5:lvar).eq.'_north'
     &                           .and. northern_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')   'node =', node,
     &               'identified XI-partitioned northern ',
     &               'boundary array ''', vname(i)(1:lvar), ''''
#endif
                    endif
                  elseif (part_type(i).eq.2 .and. lvar.gt.5) then
                    if (vname(i)(lvar-4:lvar).eq.'_west'
     &                                .and. western_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')     'node =', node,
     &               'identified ETA-partitioned western boundary ',
     &               'array ''', vname(i)(1:lvar), ''''
#endif
                    elseif (vname(i)(lvar-4:lvar).eq.'_east'
     &                          .and. eastern_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')     'node =', node,
     &               'identified ETA-partitioned eastern boundary ',
     &               'array ''', vname(i)(1:lvar), ''''
#endif
                    endif
                  elseif (part_type(i).eq.3) then
                    var_mask=.true.
#ifdef VERBOSE
                    write(*,'(3x,A,I4,1x,4A)')   'node =',  node,
     &                      'identified 2D-partitioned array ''',
     &                       vname(i)(1:lvar), ''''
#endif
                  endif
                  ! DevinD can remove this as start,count not used below this anymore.
                  if (var_mask) then

                    do j=1,vdims(i)
                      k=dimids(j,i)
                      if (k.eq.unlimdimid) then
                        start_in(j,in2p)=rec ! DevinD: set input  start to current timestep
                        start_out(j)=rec     ! DevinD: set output start to current timestep
                      endif
                    enddo

#ifdef VERBOSE
                      write(*,'(3x,A,I4,2x,A,I7,2x,A,I4,1x,A,I4,1x,A)')
     &                     'node =', node, 'ncid=',ncid(node),
     &                     'xi_start =', xi_start(node),
     &                     'eta_start =', eta_start(node), 'reading...'
#endif

# ifdef TIMING
                      nclk=3-nclk
                      call system_clock (iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      net_gray_clk=net_gray_clk+inc_clk
# endif

                      if (vartype(i) .eq. nf_char) then

                      elseif (vartype(i) .eq. nf_float) then

                        ierr=nf90_get_var(ncid(node), vid(i,node),
     &                     data(data_x_in_str(in2p):data_x_in_end(in2p),
     &                          data_y_in_str(in2p):data_y_in_end(in2p),:),
     &                          start_in(:,in2p), count_in(:,in2p) )
                        if(ierr/=0) print *, 'BAD READ', nf90_strerror(ierr),
     &                                        'my_rank', my_rank

                      elseif (vartype(i) .eq. nf_double) then
!                        ierr=nf_get_vara_double(ncid(node),vid(i,node),
!     &                                             start, count, buff)
                      endif
# ifdef TIMING
                      if (ierr.eq.nf_noerr) then
                        net_read_size=net_read_size+size
                        nclk=3-nclk
                        call system_clock (iclk(nclk),clk_rate,clk_max)
                        inc_clk=iclk(nclk)-iclk(3-nclk)
                        net_read_clk=net_read_clk+inc_clk
                      else
# else
                      if (ierr.ne.nf_noerr) then
# endif
                        lvar=lenstr(vname(i))
                        lncn=lenstr(ncname(node))
                        write(*,'(/1x,4A,I3/15x,3A/15x,A/)')  '### ',
     &                              'ERROR: Cannot read variable ''',
     &                   vname(i)(1:lvar), ''' for time record', rec,
     &                  'from netCDF file ''',  ncname(node)(1:lncn),
     &                                      '''.', nf_strerror(ierr)
!dd                        goto 97
                      endif

                  endif  ! <-- var_mask

                enddo   !<-- in2p=1,nnodes_in_proc

#ifdef TIMING
                nclk=3-nclk  ! DevinD moved here
                call system_clock (iclk(nclk), clk_rate, clk_max)
                inc_clk=iclk(nclk)-iclk(3-nclk)
                net_gray_clk=net_gray_clk+inc_clk
#endif

#ifdef MPIVERBOSE
                print *, 'WRITING FLOAT: ', 'my_rank', my_rank, 'tstep= ', rec
#endif

                ierr=nf90_put_var(nctarg, varid(i),
     &                            data, start_out, count_out)
                if(ierr/=nf90_noerr) print *, 'BAD WRITE!',
     &            nf90_strerror(ierr),'my_rank', my_rank


#ifdef TIMING
                if (ierr.eq.nf_noerr) then       ! DevinD moved here
                  net_wrt_size=net_wrt_size+size

                  nclk=3-nclk
                  call system_clock(iclk(nclk), clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  net_wrt_clk=net_wrt_clk+inc_clk
                endif
#endif

      end subroutine read_write_partitioned_var  !]

! ----------------------------------------------------------------------
      subroutine display_timing_summary  ![
      implicit none

#ifdef TIMING
      call etime(CPU_time, RUN_time)
      RUN_time=RUN_time-tstart

      ! DevinD commented this out as etime now takes array not 2 arguements. Might fix in future.
!      if(my_rank==0) write(*,'(/3(1x,A,F11.2,1x))') 'CPU_time:  run =', RUN_time,
!     &                   'usr =', CPU_time(1),  'sys =', CPU_time(2)

      if (clk_rate.gt.0) then
        ReadSize=1.0D-6*net_read_size
        WrtSize=1.0D-6*net_wrt_size
        ReadTime=net_read_clk/dble(clk_rate)
        AssmTime=net_assm_clk/dble(clk_rate)
        WrtTime = net_wrt_clk/dble(clk_rate)
        SyncTime=net_sync_clk/dble(clk_rate)
        FcrtTime=net_fcrt_clk/dble(clk_rate)

        if(my_rank==0) write(*,'(/1x,A,22x,F12.2,1x,A)')
     &      'Analysis/file creation :', FcrtTime, 'sec'
        if(my_rank==0) write(*,'(7x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &         'Master data read :', ReadSize, 'MBytes in',  ReadTime,
     &                          'sec (', ReadSize/ReadTime, 'MB/sec)'
        if(my_rank==0) write(*,'(8x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)') ! DevinD
     &   'Total data read :', ReadSize*n_procs, 'MBytes in', ReadTime,
     &   'sec (', ReadSize*n_procs/ReadTime, 'MB/sec)'
        if(my_rank==0) write(*,'(4x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &      'Master data written :', WrtSize,  'MBytes in',   WrtTime,
     &                          'sec (',  WrtSize/WrtTime,  'MB/sec)'
        if(my_rank==0) write(*,'(5x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &      'Total data written :', WrtSize*n_procs,  'MBytes in', WrtTime,
     &      'sec (',  WrtSize*n_procs/WrtTime,  'MB/sec)'
        if(my_rank==0) write(*,'(2x,A,22x,F12.2,1x,A)')
     &       'Output file sync time :',  SyncTime, 'sec'

# ifdef DEL_PART_FILES
        if(my_rank==0) then  ! Master only
          if (del_part_files) then
            write(*,'(1x,A,22x,F12.2,1x,A)')
     &      'Removing partial files :',net_rmcmd_clk/dble(clk_rate),'sec'
          endif
        endif
# endif
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
        GrayTime=net_gray_clk/dble(clk_rate)
        if(my_rank==0) write(*,'(14x,A,22x,F12.2,1x,A)')
     &                 'Gray time :',GrayTime,'sec'
        inc_clk=iclk(nclk)-iclk_init
        if(my_rank==0) write(*,'(47x,A/12x,A,11x,F12.2,1x,A/)')
     &                 '------------------',
     &   'Elapsed wall-clock time:', inc_clk/dble(clk_rate), 'sec'
      endif
#endif

      end subroutine display_timing_summary  !]

! ----------------------------------------------------------------------
      subroutine close_file_set  ![
      implicit none

! Close all files                        ! At this moment open/closed
                                         ! status of partial files
        if(my_rank==0) write(*,'(/1x,A)') 'closing files...'    ! depends on the state of CPP ! was a 98 here
        do node=0,nnodes-1               ! switch KEEP_CLOSED.   If the
          if (ncid(node).ne.-1) then     ! switch is defined, then only
            ierr=nf_close(ncid(node))    ! node=0 file is expected to
            ncid(node)=-1                ! be opened here.   Otherwise
          endif                          ! the entire set is opened and
        enddo                            ! needs to be closed.   Either
        if(my_rank==0) write(*,*) '...........input'    ! way, as ncid(node).eq/ne.-1
#if defined TIMING
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
#endif
        ierr=nf_close (nctarg)           ! is used as flag indicating
#if defined TIMING
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_sync_clk=net_sync_clk+inc_clk
#endif
        if(my_rank==0) write(*,*) '...........output'   ! status of each file.

#ifdef DEL_PART_FILES
        if(my_rank==0) then  ! Master only to delete partial files
          if (del_part_files) then
            if (clean_set) then
# ifdef TIMING
              nclk=3-nclk
              call system_clock (iclk(nclk), clk_rate, clk_max)
              inc_clk=iclk(nclk)-iclk(3-nclk)
              net_gray_clk=net_gray_clk+inc_clk
# endif
              write(*,'(/1x,A)') 'Deleting partial files...'
              do node=0,nnodes-1
                rmcmd='/bin/rm -f '/ /ncname(node)
                lstr=lenstr(rmcmd)
                if (node.lt.16 .or. (nnodes.gt.16 .and.
     &                           node.eq.nnodes-1 )) then
                  write(*,'(27x,3A)')  '''', rmcmd(1:lstr), ''''

                elseif (nnodes.gt.16 .and. node.lt.18) then
                  write(*,'(24x,A)') '.................................'
                endif
                call system (rmcmd(1:lstr))
              enddo
              write(*,*)
# ifdef TIMING
              nclk=3-nclk
              call system_clock (iclk(nclk), clk_rate, clk_max)
              inc_clk=iclk(nclk)-iclk(3-nclk)
              net_rmcmd_clk=net_rmcmd_clk+inc_clk
# endif
            else
              write(*,'(/1x,2A/)')  '### ERROR: Not removing ',
     &                       'partial files because of errors.'
            endif
          endif ! <- del_part_files
        endif   ! <- my_rank==0
#endif

      end subroutine close_file_set  !]

! ----------------------------------------------------------------------
! JOIN_CHILD_BRY ROUTINES:
! ----------------------------------------------------------------------
      subroutine determine_global_child_bry_dims      ![
      implicit none

      ! Determine the size of global child dimensions
      ! in order to create joined child file

      ierr=nf90_get_att(ncid(0), NF90_GLOBAL, 'global_south_size', global_south_bry_size)
      ierr=nf90_get_att(ncid(0), NF90_GLOBAL, 'global_west_size', global_west_bry_size)

      end subroutine determine_global_child_bry_dims  !]

! ----------------------------------------------------------------------
      subroutine create_joined_child_file  ![
      implicit none

      call create_joined_empty_file
      if(mayday) return

      call copy_global_attributes
      if(mayday) return

      ! Define dimensions

      ierr=nf_def_dim (nctarg, 'xi_rho',                 ! dimname(i)(1:lvar),
     &                 global_south_bry_size, dimid(1))  ! hardcoded 1 and dimsize(i,nnodes)
      ierr=nf_def_dim (nctarg, 'eta_rho',
     &                 global_west_bry_size,  dimid(2))
      ierr=nf_def_dim (nctarg, 'rec_time', 0, dimid(3))  ! 0 for unlimited

      ! Define variables

      ierr=nf90_def_var (nctarg, 'i_south', NF90_FLOAT, ! vartype(i),
     &                  (/dimid(1),dimid(3)/), vid_out(1)) ! dimids(1:vdims(i),i), varid(i)
      ierr=nf90_def_var (nctarg, 'i_west', NF90_FLOAT,
     &                  (/dimid(2),dimid(3)/), vid_out(2))
      ierr=nf90_def_var (nctarg, 'i_north', NF90_FLOAT,
     &                  (/dimid(1),dimid(3)/), vid_out(3))
      ierr=nf90_def_var (nctarg, 'i_east', NF90_FLOAT,
     &                  (/dimid(2),dimid(3)/), vid_out(4))

      ! Leave definition mode

!      ierr=nf90_close(ncid(0)) ! close last input file. completed master only reads
!      ncid(0)=-1

      ierr=nf_enddef (nctarg)
!      ierr=nf90_close(nctarg)  ! no need to close since serial

# ifdef VERBOSE
      if(my_rank==0) write(*,'(/1x,A)') 'Leaving definition mode.'
# endif

      end subroutine create_joined_child_file  !]

! ----------------------------------------------------------------------
      subroutine read_write_child_bry  ![
      implicit none

      ! local
      integer :: in_varid, len_dim
      integer, dimension(2)           :: child_start, child_count=1, dimids
      real, dimension(:), allocatable :: bry_data
      character(len=12) :: dimname

      ! get start and count attributes
      if (ncid(node).eq.-1)
     &  ierr=nf90_open (ncname(node), nf90_nowrite, ncid(node))

      ierr = nf90_get_att(ncid(node), vid(i,node), 'start', child_start(1))
      ierr = nf90_get_att(ncid(node), vid(i,node), 'count', child_count(1))

      if (child_count(1) /= 0) then  ! only write if parent has child boundary (0 if not)

        ! read variable
        if (allocated(bry_data)) deallocate (bry_data)
        allocate (bry_data(child_count(1)))
        ierr = nf90_get_var( ncid(node), vid(i,node), bry_data )
        if(ierr/=0) print *, 'read ierr=', nf90_strerror(ierr), 'node=', node

        ! write variable
        child_start(2)=1 ! rec  ! set recorded timestep
!        ierr = nf90_inq_varid( nctarg, 'i_south', vid_out )

!        ierr = nf90_inquire_variable(nctarg, vid_out, dimids=dimids)
!        ierr = nf90_inquire_dimension(nctarg, dimids(1), dimname, len_dim)
!        print *, 'dim1=',dimname, 'len=',len_dim, 'node=',node
!        ierr = nf90_inquire_dimension(nctarg, dimids(2), dimname, len_dim)
!        print *, 'dim1=',dimname, 'len=',len_dim, 'node=',node
!        print *, 'child_start=',child_start, 'node=',node
!        print *, 'child_count=',child_count, 'node=',node

        ierr = nf90_put_var( nctarg, vid_out(i), bry_data,
     &                       start=child_start, count=child_count )

!        print *, 'wrt ierr=', nf90_strerror(ierr), 'node=', node

      endif

      ! Need to close and sync for time being.
      ierr=nf90_sync (nctarg) ! debug only

      end subroutine read_write_child_bry  !]

! ----------------------------------------------------------------------

      end module ncjoin_mod
