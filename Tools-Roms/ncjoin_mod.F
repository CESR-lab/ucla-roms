      module ncjoin_mod  ! contains variables and generic routines to be used in both
                         ! ncjoin and ncjoin_mpi to avoid duplications and for clarity.

      use mpi

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


      contains

! ----------------------------------------------------------------------
      subroutine setup_mpi  ! setup MPI system and
                            ! error checking of user MPI config
      implicit none

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

      if(n_procs /= nprocs_x*nprocs_y) then     ! Abort if wrong mpi sub-domain tiling
        print *, 'Error: number of mpi processes /= nprocs_x * nprocs_y (sub-domains)'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      endif

      if(my_rank==0) write(*,'(/1x,A,/1x,A,I3,2(/1x,A,I2)/)')
     &                'MPI procs:', 'n_procs= ', n_procs,
     &                'nprocs_x= ', nprocs_x, 'nprocs_y= ', nprocs_y

      end subroutine setup_mpi

! ----------------------------------------------------------------------
      subroutine check_partial_file_set  ! determine is there is a complete
                                         ! set of partial files to join

  1     nnodes=-1
        root_bak(1:1)=' '  ! Reset variables which
        sffx_bak(1:1)=' '  ! define the file set.
        code_size_bak=-1

  2     arg=arg+1 ! cycle through argument list of names of partitioned files by
                  ! coming back here with 'goto 2' after each filename.

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

                  if (ierr.ne.nf_noerr) return ! goto 97 can't have go to
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
              complete=.false.   ! Use this to skip parallel read/write section of code.
              return ! goto 23            ! goto end of mater only pre-processing.
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

#ifdef VERBOSE
      if(my_rank==0) write(*,*) ' line 433, complete =', complete,
     &                          ' nnodes =', nnodes
#endif

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

      end subroutine check_partial_file_set
! ----------------------------------------------------------------------

      end module ncjoin_mod
