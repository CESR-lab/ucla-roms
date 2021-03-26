
      program ncjoin_mpi  ! Join partial result files into one full-domain file

      ![ NCJOIN INFO:

      ! Creator: Sasha (Alex S.)      (Serial version)
      !          Devin.D & Jeroen.M   (  MPI  version)

      ! Edits:
      ! - 2021-Jan: DevinD changed nf_create and nf_def_var to nf90
      !             routines to create netcdf-4 files to allow for
      !             file compression.
      !             Change compression level with 'deflate_level='.
      ! - 2021-Feb: DevinD added parallel I/O functionality using MPI.

! Generic netCDF assembly tool: reads ROMS partial netCDF files and
! assembles them into a file for the whole physical grid.
! It is the inverse operation to partit.
! MPI Usage:
!
!       mpiexec -n np ncjoin_mpi np_x np_y files.???.nc
! or
!       mpiexec -n np ncjoin_mpi np_x np_y -d files.???.nc
!
! where files.???.nc matches a complete set of partial files (or
! several complete sets) and "-d" or "--delete" forces deletion of
! partial files upon successful joining. 

! np_x and np_y are the number of sub-domains you choose in x and y
! for ncjoin_mpi not your input partition files, they can be different!
! For efficiency, try to keep keep the ratio of np_x to np_y similar to
! the sub-domaining of your partitions.

! CPP-switch "DEL_PART_FILES" enables the user-activated option to
! delete partial files just after they are joined.  If compiled with
! this switch is defined, ncjoin becomes sensitive to whether its first
! argument is "-d" or "--delete".  The deletion occurs one-set-at-at-
! a-time and ONLY if the joining is successful, i.e., if an error of
! any kind occurs, so the set cannot be joined for whatever reason (a
! missing file, or one file has fewer records than others, or reading
! error; etc), the files stay intact.  If ncjoin is interrupted, the
! partial files stay intact.

!] --------------------------------------------------------------------

! PROGRAM SECTIONS:
!           (parts 1-3 -> master proc only)
! - PART 1: Confirm all input files are correct and related to each other
! - PART 2: Veryify all variables and dimensions in input files
! - PART 3: Create output file
!           (parts 4+  -> all procs)
! - PART 4: Exchange all partial file info from master to all procs
! - PART 5: Collectively (MPI-I/O) open all input and output files
! - PART 6: Assign read in partial files to neighbouring MPI procs
! - PART 7: Read/write scalar variables
! - PART 8: Read/write non-partioned variables (e.g. time_step)
! - PART 9: Read/write partitioned variables   (bulk of work done here)
! - PART10: Close all files and display program timing summary

![NCJOIN_MPI INFO:

![LIMITATIONS:
! Currently ncjoin_mpi only works on Xsede's Expanse platform, as it has up-to-date
! netCDF libraries that can do parallel I/O with compression.
! Xsede's comet can handle parallel I/O but not with compression.
! Maya currently does not support parallel I/O at all.
! Scaling limitations: based on initial testing during implementation, the reading
! of input files was scaling well with increasing numbers of procs. However, the
! writing was not scaling very well, with diminishing return beyond a certain point. !]

![USEFUL LINKS:
! https://www.archer.ac.uk/documentation/white-papers/fortanIO_netCDF/fortranIO_netCDF.pdf
! http://www.spscicomp.org/ScicomP14/talks/Latham.pdf
! https://www.researchgate.net/publication/332190037_Best_Practice_Guide_-_Parallel_IO
!]

! MPI SUB-DOMAINS OF NCJOIN_MPI:
! Unlike the old version of ncjoin that ran on one core, ncjoin_mpi uses many cores (procs).
! The work to join files is shared among procs by dividing the full domain into
! sub-domains within x & y. This is similar (but not the same) to how roms uses sub-domains.
! However, the sub-domain configuration of ncjoin_mpi is completely independent from
! the sub-domain configuration of the roms partial files to be joined!
! For roms we use inode and jnode to locate a sub-domain in the full domain.
! For ncjoin_mpi we use iproc and jproc to locate a sub-domain in the full domain.
! Note again, iproc/jproc does not need to match inode/jnode.
! In that sense, you should avoided conflating roms sub-domains (referred to henceforth
! as nodes, partitions or partial files) with that ncjoin_mpi sub-domains (referred to
! as mpi processes or procs).
! For example, a set of partial result files might come from a roms simulation with
! 30x20 sub-domains (and hence files from 000-599). However, we may chose to run
! ncjoin_mpi to join those result files with 10x5 sub-domains using 50 cores.
! The potentially uncoupled nature of roms partitioning to ncjoin_mpi proc sub-domains
! is illustrated in EXAMPLE A and DIAGRAM A below.

! NEED FOR CHUNKING:
! In order to compress variables with netCDF (and the underlying HDF5 library), a
! variable must first be sub-divided into uniform chunks along each of its dimensions.
! The reason for this is when you read a specific bit of compressed data, the entire
! chunk that contains it needs to be decompressed to get that specific point. So if
! there was no chunking, all of the variable's data would need to be decompressed (costly!)
! each time you want a small section of its data. Hence, the need for chunking the data.
! For example: if we have temperature 'T' with dimensions (1000x500x100), and with
! 10 timesteps. We could then divide 'T' into chunks of 100x50x100x1, where we have
! divided x & y into 10x10 chunks, we have left z as 100, we have split each timestep
! into 1.
! Note, chunksize selection affects speed of data access depending on how you plan to use
! the data. If you want to look at the full domain at a single timestep, then the above
! chunking of 'T' is useful as we have one chunk per time step. If you wanted to get a
! timeseries at a single point, this would not be so efficient as you'd have to decompress
! many chunks to access everytime step at a point. Currently we have prioritized domain
! access (for ncview) rather than time series access.

! NEED UNIFORM MPI PROC SUB-DOMAINS FOR EFFICIENT WRITES:
! The principle behind this routine is that in order for efficient writing
! to be done,and hence the compression (which is the cpu expensive process),
! it is necessary to have each mpi process (proc) write into whole chunks where possible.
! Hence, an mpi proc range should overlap exactly with the prescribed chunk sizes of
! the variable.

! ROMS HAS NON-UNIFORM SUB-DOMAINS:
! In roms, the sub-domains are created such that the internal sub-domains are
! all the same size, but the boundary sub-domains are modified such that
! the total number of grid points in x & y match the domain. I.e. the E and W
! boundary sub-domains might be a different size in x to the internal sub-domains.
! The N and S boundary sub-domains might be a different size in y to the internal
! sub-domains.

! NETCDF NEEDS UNIFORM CHUNKSIZING:
! However, with the netCDF library (and thus its underlying HDF5 library), it
! is not possible to have non-uniform chunksizes. Consequently, it is not
! possible to efficiently write the roms sub-domains directly into the joined file
! with compression. (Note, it is possible without compression as no chunking is
! required when uncompressed, but obviously we want compression to save space).

! RANGE OF NCJOIN_MPI PROCS TO FOLLOW CHUNKSIZE:
! The writing range to the joined output file of each MPI process must therefore follow
! the prescribed chunksizing. As such, starting from the SW corner of the domain,
! each proc will write complete chunks, until the eastern boundary and northern
! boundary of the full domain where it is possible 1 column/row of partially filled
! chunks will be written. (See DIAGRAM A below - procs C and D do not fill thier
! northern most chunks, and procs B and D do not fill their eastern most chunks).

! The need for this is that if procs had to share overlapping chunks, there is
! contention between procs for that memory, and when the second proc goes to
! add its data to the chunk, it will have to fully uncompress the existing partial
! data in the chunk, and recompress the whole chunk. This is inefficient...

! ASSIGN READ IN PARTITIONS TO PROCS:
! Since the MPI proc write range does not match the roms partitions, we must
! decide which partitions are to be read in by each proc. Unfortunately here
! there will be some repeated work, as ingoing partitions that overlap MPI proc
!  ranges will have to be read (and decompressed if roms output is compressed) twice.
! See INFO A) below for more details.

! ROUTINE VARIABLE INDICES:
! There are two sets of variable indices for data read in, and two sets for data
! written out. Match labels R1/2 & W1/2 with variable declarations in the
! code below. Note, variable's with 'p_' belong to each mpi proc of ncjoin_mpi.
! '_in_'  refers to variables involved in reading in data.
! '_out_' refers to variables involved in writing outt data.

  ! READ VARIABLES:

    ! R1) Index ranges as per each partitioned input file.
    !     Read in from file: these are ranges of the data in the input
    !     partition file.

    ! R2) Index ranges as per MPI proc temporary array.
    !     These are the ranges with which to populate the temporary array 'data'
    !     to store the read in data. Indices local to 'data' array.
    !     Since the MPI write ranges don't match the roms partition sub-domain sizes,
    !     these ranges differ from R1).

    !       partition_infile( R1 ) -> data( R2 )

  ! WRITE VARIABLES:

    ! W1) Index ranges local to MPI proc temporary array.
    !     These are the ranges to take from the temporary array
    !     of each MPI proc to write out. The indices are local the temporary
    !     array 'data' of the MPI proc.

    ! W2) Index ranges as per full model domain.
    !     Write to joined file: these are ranges of the full domain to write
    !     to the joined file from a proc's 'data' variable.

    !       data( W1 ) -> joined_outfile ( W2 )


! EXTRA INFO:
! -----------

! INFO A) ASSIGN READ IN PARTITIONS TO PROCS (cont.):

! Once we've calculated the MPI proc's write ranges to the joined full
! domain outfile, the SW corner and NW corner of those ranges is determined.
! - p_str_out_x, p_end_out_x
! - p_str_out_y, p_end_out_y

! From there we need to find all input partitions that contribute to this
! range. To do so we calculate x & y of the 4 corners of input partition:
! x_nc_sw, y_nc_sw, ...
! Then we see if those corners fall within the proc's corners.
! If so, then proc must read in some or all data from that partitions.

! EXAMPLE A:
! Say we use ncjoin_mpi with 4 procs. Let's choose a 2x2 split in x & y sub-domains
! of the proc's, labelled A-D with double lines in the diagram below.
! Let's say the partial result files had a partitioning of 3x2, labelled 0-5
! with single lines.
! If total number of grid points in x was not exactly divisible by nprocs_x (i.e. 2),
! the ncjoin_mpi procs on the eastern edge would have partial incomplete ranges.
! Similarly if the num. grid points in y was not exactly divisible by nprocs_y (i.e. 2),
! the ncjoin_mpi procs on the northern edge would have partial incomplete ranges.

! We can assign read in partial files to procs as follows:
! 0 & 1 -> A,    1 & 2 -> B,    3 & 4 -> C,    4 & 5 -> D

!    ═══════════════════════════════════
!   ║─────────────────║──────────────   ║
!   ║          │      ║     │        │  ║
!   ║        C │      ║     │  D     │  ║
!   ║     3    │     4║     │    5   │  ║        (DIAGRAM A)
!   ║          │      ║     │        │  ║
!   ║          │      ║     │        │  ║     - nodes (partitions 0-5)
!   ║─═─═─═─═─═─═─═─══─═─═─═─═─══─═─════║
!   ║          │      ║     │        │  ║     = mpi procs of ncjoin_mpi (A-D)
!   ║          │      ║     │        │  ║
!   ║        A │      ║     │  B     │  ║   -=- overlap of proc range and partition range
!   ║     0    │     1║     │    2   │  ║
!   ║          │      ║     │        │  ║
!   ║          │      ║     │        │  ║
!    ═══════════════════════════════════

! Unfortunately partial files 1 and 4 are each read in by 2 procs.
! E.g. the left part of 1 is read by proc A and the right part of 1 by proc B.
! This is not ideal but deemed unavoiable with our method.

! Once we have assigned partial files (nodes) to proc's, then each proc just loops
! over those specific nodes to fill its 'data' variable.

![ISSUES TO ONE DAY RESOLVE:
!
! Note, initially it was attempted to parallelize ncjoin using openMP. It turns out
! netCDF does not have good support for openMP, so this approach was abandoned.
!
! Parallel I/O with netcdf has been around since atleast 2010. However, parallel I/O
! with compression is very recent, only really since 2020. NetCDF uses hdf5 for netcdf 4
! files, and only recently did hdf5 allow for parellel compression. Subsequently netCDF
! followed suit. It is my suspicion that the netCDF wrapper around hdf5 might improve
! with subsequent updates, so the performance from ncjoin_mpi might improve passively
! with netCDF library upgrades.
!
! Collective/Independent access to netcdf files:
! This is something that needs to be fine tuned as I found conflicting information.
! The official netcdf fortran 4.5.3 webpage says:
! "Set collective access on this variable. This will cause all
! reads/writes to happen together on every processor.
! call handle_err(nf90_var_par_access(ncid, varid, nf90_collective))"
! My intuition suggests that `independent' writes would be the most efficient with many
! cores, since files can be accessed independently. However, many sources said collective
! setting should be used. When I had the code written without only master doing the
! pre-processing, I was able to use both collective or independent modes, and initial
! testing gave me similar timings.
!
! However, when I wrote ncjoin_mpi such that only the master did the pre-processing,
! it required closing the files and re-opening the newly created file. For some reason I
! could only get it working in collective mode and not independent, but this test should
! be looked into further at some point. This is likely the way I am calling the parallel
! I/O commands, but could also be bugs in the netCDF wrappers around the underlying hdf5
! library.

! Chunk sizing:
! Currently the chunk sizes are split such there is one chunk per MPI proc. If you are
! using many cores with ncjoin-mpi this is probably efficient, but if you have a huge
! domain and use very few cores, these chunks might be too big. An investigation into the
! effect of chunk size would be nice. Prevailing wisdom online seems to be 'bigger the
! better' provided several chunks still fit into your chunk cache.
!
! Legacy vars and code:
! Haven't fully cleaned up old ncjoin vars and code that is not actually used anymore.
!
! Misc:
! Improve scaling of writing.
! Create module and move sections into subroutines.
! guess_nnodes variable needs a better choice.
! Memory demand is higher on expanse than it was for ncjoin. All the auxiliary variables
! are using up too much memory, I need to streamline this.
! change all netcdf calls to nf90.
! non-partioned sections it doesn't catch all types of variables as it did in ncjoin.
! Could have used global incides on proc's data array. I.e. for a 9x9 ncjoin_mpi run
! the proc range with the NE corner of the global domain for 300x300 domain could be
! data(200:300,200:300,z). So really the local array is 100x100 in x and y, but the
! indices stick with global indices, and therefore don't have to carry extra index
! conversion.
! Could print out each proc's write time to see more about bad write scaling.
! Currently I use nf90_put_var and nf90_get_var using multi-dimensional array instead
! of 1D buffer. E.g nf90_put_var(ncid, varid, data(istr:iend,jstr:jend,z:..). But it's
! possible that a local copy is made into a buffer by netCDF, this would be inefficient.
!
!]

!]    ===================================================================

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

      use ncjoin_mod
      use netcdf  ! for nf90_create & nf90_def_var

      implicit none

      ! remove netcdf.inc once fully converted to nf90_ routines
!#include "netcdf.inc"  ! added to ncjoin_mod
      ! Change this to 'use mpi' for module rather... as per 'simple_xy_par_wr2.f90'
!#include "mpif.h" ! added use mpi to ncjoin_mod

      ! ***************************************************************
      ! START USER INPUT          ! mostly you will use the default values!
      integer :: deflate_level=1  ! compression level for joined file. For no compression 'deflate_level = 0'.
      logical :: shuffle=.true.   ! shuffle on for extra conmpression (nc3to4z had this).
                                  ! in tests shuffle gave an extra ~10% compression over just d1 compression.
      integer :: guess_nnodes=16  ! guessed number of input nodes (files) linked to mpi proc for storage array size
                                  ! Must have guess_nnodes > nnodes_in_proc! But try to keep this number low for good memory management.
      !   END USER INPUT          ! ncjoin_mpi will tell you at run time if number is too low (i.e. guess_nnodes < nnodes_in_proc).
      ! ***************************************************************

      ![ DECLARATION OF VARIABLES:

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
      nargs=iargc()
      arg=2          ! start at 2+1 since first 2 arguements are mpi sub-domain values np_x & np_y

      !] ---------------- PROGRAM STARTS ---------------------


      call setup_mpi

  11  nnodes=-1 ! DevinD repeated to get 'goto' outside of master only region, as called by all procs.

      if(my_rank==0) then  ! MPI MASTER ONLY: PRE-PROCESSING. i.e. check input files, creating output file, etc.
                           ! Extract a set of files which cover the whole physical grid.

        call check_partial_file_set    ! PART 1: CHECK ALL INPUT PARTIAL FILES ARE CORRECT
        if     (ierr.ne.nf_noerr) then
          goto 97  ! Taken from within check_partial_file_set as can't have goto
        elseif (complete==.false.) then
          goto 23  ! Taken from within check_partial_file_set as can't have goto
        endif

      ![ PART 2: VERIFY ALL VARIABLES & DIMENSIONS:

! At this moment a set of files recorded as ncname(0:nnodes-1),
! xi_start(0:nnodes-1), eta_start(0:nnodes-1) comprise a complete
! set, but the files are actually in "closed" state and all netCDF
! IDs are reset to -1.

! Verify that ndims, ngatts, unlimdimid are the same for all nodes,
! however, note that different files may store different composition
! of variables, and netCDF variable IDs for the same variable (with
! the same name) may be different across the set of files.

#ifdef TIMING
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
#endif
        do node=0,nnodes-1
          lncn=lenstr(ncname(node))

!          if (ncid(node).eq.-1) ierr=nf90_open (ncname(node), ! Original MPI-IO version
!!     &         IOR(NF90_NETCDF4, NF90_MPIIO), ncid(node), ! not sure if I need nf_nowrite in here somehow
!     &         IOR(NF90_NOWRITE, NF90_MPIIO), ncid(node), ! DevinD NOTE SURE nowrite or netcdf4. netcdf4 seemed to be slower on laptop test.
!     &         comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)

          if (ncid(node).eq.-1) ierr=nf90_open(ncname(node),nf90_nowrite, ncid(node))

          if (ierr .eq. nf_noerr) then

            ierr=nf_inq (ncid(node), ibuff(1), ibuff(2),
     &                               ibuff(3), ibuff(4))
            if (ierr .ne. nf_noerr) then
              write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot make ',
     &                        'general inquiry into netCDF file ''',
     &               ncname(node)(1:lncn), '''.', nf_strerror(ierr)
              goto 97
            elseif (ibuff(1) .gt. maxdims) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of dimensions', ibuff(1), 'in netCDF file ''',
     &          ncname(node)(1:lncn),    '''',   'exceeds limit.  ',
     &             'Increase parameter maxdims in file "ncjoin.F".'
              goto 97
            elseif (ibuff(2) .gt. maxvars) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of variables',  ibuff(2), 'in netCDF file ''',
     &          ncname(node)(1:lncn),  '''',      'exceeds limit. ',
     &             'Increase parameter maxvars in file "ncjoin.F".'
              goto 97
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
              if (ierr .ne. nf_noerr) goto 97 
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
                  goto 97
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &            'Cannot determine name and size of dimension #', j,
     &            'in netCDF file ''',   ncname(node)(1:lncn), '''.',
     &                                             nf_strerror(ierr)
                goto 97
              endif
            enddo
#undef i

! Create catalog of variable names, IDs, and ranks (number of
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
!                 where it is not needed).   

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
                  goto 97
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &                 'Cannot determine name of variable with id =',
     &                  i, 'in netCDF file ''', ncname(node)(1:lncn),
     &                                      '''.', nf_strerror(ierr)
                goto 97
              endif
            enddo

            if (node.gt.0) then           ! close all the files, except for node=0.
              ierr=nf_close(ncid(node))   ! since master only, need to close files to open collectively later.
              ncid(node)=-1               ! keep node=0 open as still using below.
            endif

          else
            write(*,'(/1x,A,1x,3A/14x,A)')    '### ERROR: Cannot ',
     &                 'open netCDF file ''', ncname(node)(1:lncn),
     &                                    '''.', nf_strerror(ierr)
            goto 97
          endif
        enddo  !<-- node=0,nnodes-1 

#ifdef VERBOSE
        if(my_rank==0) then
          write(*,'(/1x,A,I3)') 'Inventory of variables: nvars =',nvars
          do i=1,nvars
            lvar=lenstr(vname(i))
            write(*,'(4I4,2x,3A)') i, vid(i,vnode(i)), vnode(i),
     &                   vdims(i), '''', vname(i)(1:lvar), ''''
          enddo
          write(*,*) '...............................'
        endif
#endif

! Determine sizes of dimensions for combined file: For partitionable
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
! is still executed, but only once.

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
                goto 97
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


#ifdef VERBOSE
      if(my_rank==0) write(*,'(1x,A)') 'Identifying presense of boundary edges:'
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


! set buffer size for needed to accommodate data for the largest
! possible MPI-subdomain (i.e., partial file).
! 
        max_buff_size=size_XI*size_ETA*size_S

      !] END PART 2: VERIFY ALL VARIABLES & DIMENSIONS

      ![ PART 3: CREATE OUTPUT FILE:

! Create combined netCDF file:   Once the completeness of the set of
!------- -------- ------ -----   partial files have been established
! and dimensions survive consistency check, create the combined file,
! define its dimensions and copy global attributes.

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
          goto 97
        endif

! Define dimensions: also compute the size of buffer needed to
! accommodate the largest array.

# ifdef VERBOSE
        if(my_rank==0) write(*,'(/1x,A,3x,A,1x,A,1x,A)')
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
              if(my_rank==0) write(*,'(14x,I3,I5,1x,3A)') dimid(i),
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
              goto 97
            endif
          else
           write(*,'(/1x,4A/12x,A/)')      '### ERROR: Cannot define ',
     &    'dimension ''', dimname(i)(1:lvar), '''.', nf_strerror(ierr)
              goto 97
          endif 
        enddo ! <-- i loop over dimensions

! Copy all global attributes, except 'partition'.

# ifdef VERBOSE
        if(my_rank==0) write(*,'(1x,A)') 'Copying global attributes:'
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
                  goto 97
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
              goto 97
            endif
          enddo
        else
          write(*,'(/1x,A,1x,3A/14x,A)')     '### ERROR: Cannot open ',
     &   'netCDF file ''', ncname(0)(1:lncn), '''.', nf_strerror(ierr)
          goto 97
        endif

! Define variables and copy their attributes.

#ifdef VERBOSE
        if(my_rank==0) write(*,'(1x,2A)') 'Variables, their dimensions and ',
     &                                           'attributes:'
#endif

        ! Need to assign read partitions to MPI procs:
        ! Use roms subdomains method    |  ...  |
        ! to divide full domain:        | 3 4 5 |
        !                               | 0 1 2 |

        if(my_rank==0) write(*,'(/1x,A)')
     &                  'Calculating chunk sizes for each variable.'

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
                if(my_rank==0) print *,'Full domain size of var:', vname(i)
                if(my_rank==0) write(*,'(3(A,i5))')' n_x= ',n_x,' n_y= ',n_y,
     &                                             ' n_z= ',n_z
                if(my_rank==0) print *,'chunk sizing of var:', vname(i)
                if(my_rank==0) print *,'chunksize_x= ',chunksize_x(i),' chunksize_y= ',chunksize_y(i)
                if(my_rank==0) print *,'n_chunks_x= ',n_chunks_x(i),
     &                                ' n_chunks_y= ',n_chunks_y(i)
                if(my_rank==0) print *,'p_chunks_x= ',p_chunks_x(i),
     &                                ' p_chunks_y= ',p_chunks_y(i)
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
                if(my_rank==0) print *,'var: ',vname(i)(1:lvar),
     &                                ' set chunksizes: ',chunk_sizes
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
                if(my_rank==0) write(*,'(8x,3A,8I3)')   '''', vname(i)(1:lvar),
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
                      goto 97
                    endif
#ifdef VERBOSE
                    if(my_rank==0) write(*,'(16x,3A)') '''', string(1:lstr), ''''
#endif
                  else
                    write(*,'(/1x,2A,I3/12x,3A/12x,A/)') '### ERROR: ',
     &                             'Cannot get name of attribute #', j,
     &                      'for variable ''', vname(i)(1:lvar), '''.',
     &                                               nf90_strerror(ierr)
                    goto 97
                  endif
                enddo
              else
                write(*,'(/8x,5A/)') 'ERROR: Cannot define ',
     &                  'variable ''', vname(i)(1:lvar), '''.',
     &                   nf90_strerror(ierr)
                goto 97
              endif
            else

            write(*,'(/8x,2A/15x,A,I3,1x,3A/)')  '### ERROR: Cannot ',
     &        'determine name, type and attributes for variable #', i,
     &            'from netCDF file ''', ncname(node)(1:lncn),  '''.'
              goto 97
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
            goto 97
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

  23    if(my_rank==0) write(*,'(/1x,A/)')
     &                          'End of master proc pre-processing.'

      endif ! <- if(my_rank==0) END OF MASTER ONLY PRE-PROCESSING

      !] END PART 3: CREATE OUTPUT FILE

      ![ PART 4: EXCHANGE ALL PARTIAL FILE INFO FROM MASTER TO ALL PROCS:

      call MPI_Barrier(MPI_COMM_WORLD, ierr) ! Might not be necessary.

      call MPI_BCAST(complete,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(arg,   1, MPI_INT, 0, MPI_COMM_WORLD, ierr)          ! Need for catch at end of program
      call MPI_BCAST(nargs, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)

      if (complete) then  ! Covers entire read/write section of code
                          ! Needed to skip last file is not partial file (i.e. his.0000.nc)

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
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(part_type, maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(chunksize_x,     nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(chunksize_y,     nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(p_chunks_x,  nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(p_chunks_y,  nvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'

        call MPI_BCAST(dimsize, maxdims*(nnodes+1), MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(dimids,     maxdims*maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(varid,              maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(vid,         maxvars*nnodes, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(vdims,              maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'

        call MPI_BCAST(vname, maxvars*32, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(nctargname,    64, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(ncname, nnodes*64, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'

        call MPI_BCAST(western_edge,  nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(southern_edge, nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(eastern_edge,  nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(northern_edge, nnodes, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(series,       maxvars, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'

        call MPI_BCAST(vartype,  maxvars, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(xi_start,  nnodes, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(eta_start, nnodes, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(id_xi_u,        1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(id_xi_rho,      1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(id_eta_v,       1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(id_eta_rho,     1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'
        call MPI_BCAST(unlimdimid,     1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
        if(ierr/=0) print *, 'POOR EXCHANGE'

      !] END PART 4: EXCHANGE ALL PARTIAL FILE INFO FROM MASTER TO ALL PROCS

      ![ PART 5: COLLECTIVELY (MPI) OPEN ALL INPUT & OUTPUT FILES:

!     **     *     ***  *******    ***  *********  ********
!      *    ***   ***   ***   ***  ***  *  ***  *  ***    *
!       *   ***   ***   ***   ***  ***     ***     ***
!       *  *** * ***    ***   **   ***     ***     ******
!        * **  * **     ******     ***     ***     ***
!        ***    ***     ***  **    ***     ***     ***    *
!         *     **      ***   ***  ***     ***     ********
!

! Allocate necessary buffer arrays: here "max_*_size" is needed size
! of buffer array, as determined by processing of current file, while
! "alloc_*_size" is size of array already allocated.  Both are
! initialized to zero at the beginning.   Basically this anticipates
! the possibility of gradual increase of needed size of buffer array,
! so that if it is allocated at this moment, but insufficient, it is
! deallocated first.

        if (max_buff_size .gt. alloc_buff_size) then
!          if (allocated(buff)) deallocate(buff)
!          allocate(buff(max_buff_size))
          alloc_buff_size=max_buff_size
!          if(my_rank==0) write(*,*) 'allocated "buff" with  max_buff_size =',
!     &                                         max_buff_size
        endif
#ifdef TIMING
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_fcrt_clk=net_fcrt_clk+inc_clk
#endif

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
          ierr=nf90_var_par_access(nctarg, varid(i), NF90_COLLECTIVE)
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

        !] END PART 5: COLLECTIVELY (MPI) OPEN ALL INPUT & OUTPUT FILES:

        do i=1,nvars ! LOOP THROUGH VARIABLES & WRITE EACH ONE

          ![ PART 6: ASSIGN READ IN PARTIAL FILES TO NEIGHBOURING MPI PROC:

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

          !] END PART 6: ASSIGN READ IN PARTIAL FILES TO NEIGHBOURING MPI PROC:

          do rec=1,tsize ! LOOP THROUGH TIMESTEPS OF CURRENT VAR

            if(my_rank==0) write(*,'(4x,A,I8,1x,A,I8,1x,A)')
     &       'Processing record', rec, 'out of', tsize,  '...'

            ![ PART 7: READ/WRITE Scalar (zero-dimensional) variables:
            !  =======================================================

            if (rec.eq.1 .or. series(i)) then ! series is a 1D time_series variable
              if (part_type(i).eq.0 .and. .not.series(i)) then

                ! Scalar (zero-dimensional) variables:
                ! DevinD: this should also be done by one proc.

                lvar=lenstr(vname(i))
                if(my_rank==0) write(*,'(16x,3A)')
     &            'Copy scalar variable: ''',vname(i)(1:lvar),'''...'

                if (vartype(i) .eq. nf_char) then
                  ierr=nf90_get_var(ncid(0), vid(i,0),char_in)
!                elseif (vartype(i) .eq. nf_float) then
!                  ierr=nf_get_var_real   (ncid(0), vid(i,0), buff)
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A/)') '### ERROR: scalar variable ',
     &              '''', vname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (ierr .eq. nf_noerr) then
                  if (vartype(i) .eq. nf_char) then
                    ierr=nf90_put_var(nctarg ,varid(i),char_in)
!                  elseif (vartype(i) .eq. nf_float) then
!                    ierr=nf_put_var_real   (nctarg ,varid(i), buff)
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/1x,4A/12x,3A/12x,A)')  '### ERROR: ',
     &                            'Cannot write scalar variable ''',
     &                           vname(i)(1:lvar), ''' into netCDF',
     &                                'file ''', nctargname(1:ltrg),
     &                                    '''.',  nf_strerror(ierr)
                    goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot ',
     &                 'read scalar variable ''', vname(i)(1:lvar),
     &                                    '''.', nf_strerror(ierr)
                  goto 97
                endif
              elseif (part_type(i).eq.0) then !] END PART 7: READ/WRITE Scalar variables

                ![ PART 8: READ/WRITE Non-partitionable array:
                !  ===========================================

                ! This should only really be done by one proc.

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
                  goto 97
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
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/8x,4A,I3/15x,3A/)')    '### ERROR: ',
     &                 'Cannot write variable ''', vname(i)(1:lvar),
     &              ''' for time record',rec, 'into netCDF file ''',
     &                  nctargname(1:ltrg),'''.', nf_strerror(ierr)
                    goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A,I3,A/15x,A/)')     '### ERROR: ',
     &              'Cannot read variable ''',     vname(i)(1:lvar),
     &              ''' for time record',rec,'.', nf_strerror(ierr)
                  goto 97
                endif
              elseif (part_type(i).gt.0) then !] END PART 8: READ/WRITE Non-partitionable array:

                ![ PART 9: READ/WRITE Partitioned array:
                !  =====================================

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
                    size=1
                    size1=1
                    do j=1,vdims(i)
                      k=dimids(j,i)
                      if (k.eq.id_xi_rho  .or. k.eq.id_xi_u  .or.
     &                    k.eq.id_eta_rho .or. k.eq.id_eta_v) then
                        start(j)=1
                        count(j)=dimsize(k,node)

                        if (k.eq.id_xi_rho) then
                          start1(j)=xi_start(node)
                        elseif (k.eq.id_xi_u) then
                          start1(j)=max(xi_start(node)-1,1)
                        elseif (k.eq.id_eta_rho) then
                          start1(j)=eta_start(node)
                        elseif (k.eq.id_eta_v) then
                          start1(j)=max(eta_start(node)-1,1)
                        endif

                      elseif (k.eq.unlimdimid) then
                        start(j)=rec
                        count(j)=1
                        start1(j)=rec
                        start_in(j,in2p)=rec ! DevinD: set input  start to current timestep
                        start_out(j)=rec     ! DevinD: set output start to current timestep

                      else
                        start(j)=1
                        count(j)=dimsize(k,nnodes)
                        start1(j)=1
                      endif
                      size=size*count(j)
                      size1=size*count(j)
                    enddo

! Convert sizese to Bytes ! DevinD can remove this, not used.

                    if (vartype(i) .eq. nf_char .or.
     &                  vartype(i) .eq. nf_byte) then
                      size=size*1
                      size1=size1*1
                    elseif (vartype(i) .eq. nf_short) then
                      size=size*2
                      size1=size1*2
                    elseif (vartype(i) .eq. nf_int .or.
     &                      vartype(i) .eq. nf_float) then
                      size=size*4
                      size1=size1*4
                    elseif (vartype(i) .eq. nf_double) then
                      size=size*8
                      size1=size1*8
                    else
                      lvar=lenstr(vname(i))
                      write(*,'(/8x,4A/)') '### ERROR: variable ''',
     &                     vname(i)(1:lvar), ''' has unknown type.'
!dd                      goto 97
                    endif

                    if (size .gt. 8*max_buff_size) then
!                      if (allocated(buff)) deallocate(buff)
                      max_buff_size=(size+7)/8
!                      allocate(buff(max_buff_size))
!                      write(*,*) 'allocated "buff" with ',
!     &                   'max_buff_size =', max_buff_size
                    endif

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

# ifdef VERBOSE
                    write(*,'(1x,A)')  'copying ...'
# endif

! Note that all values of "count1(k)" [except for "k" corresponding to
! unlimited dimension] are equal to the actual dimensions of the
! variable in netCDF file, while all all "start(k)" [except unlimited
! dimension] are equal to 1.


!# ifdef TIMING
!                    nclk=3-nclk
!                    call system_clock (iclk(nclk), clk_rate, clk_max)
!                    inc_clk=iclk(nclk)-iclk(3-nclk)
!                    net_gray_clk=net_gray_clk+inc_clk
!# endif
                    ! original writing here
                    if (vartype(i) .eq. nf_char) then
!                    elseif (vartype(i) .eq. nf_float) then
!                       ierr=nf_put_vara_real  (nctarg, varid(i),
!     &                                     start1, count, buff)
                    endif
!# ifdef TIMING
!                    if (ierr.eq.nf_noerr) then
!                      net_wrt_size=net_wrt_size+size
!
!                      nclk=3-nclk
!                      call system_clock(iclk(nclk), clk_rate,clk_max)
!                      inc_clk=iclk(nclk)-iclk(3-nclk)
!                      net_wrt_clk=net_wrt_clk+inc_clk
!                    else
!# else
                    if (ierr.ne.nf_noerr) then
!# endif
                      lvar=lenstr(vname(i))
                      lncn=lenstr(vname(i))
                      write(*,'(/1x,3A,I3/12x,3A/12x,A/)')
     &                         '### ERROR: Cannot write variable ''',
     &                     vname(i)(1:lvar),''' for time record',rec,
     &                     'into netCDF file ''', nctargname(1:ltrg),
     &                                     '''.',  nf_strerror(ierr)
!dd                      goto 97
                    endif

                  endif  ! <-- var_mask 

                enddo   !<-- node=0,nnodes-1

! DevinD moved here:
#ifdef TIMING
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_gray_clk=net_gray_clk+inc_clk
#endif

#ifdef MPIVERBOSE
                print *, 'WRITING FLOAT: ', 'my_rank', my_rank, 'tstep= ', rec
#endif

                ierr=nf90_put_var(nctarg, varid(i),
     &                            data, start_out, count_out)
                if(ierr/=nf90_noerr) print *, 'BAD WRITE!', nf90_strerror(ierr),
     &                                           'my_rank', my_rank

! DevinD moved here:
#ifdef TIMING
                if (ierr.eq.nf_noerr) then
                  net_wrt_size=net_wrt_size+size

                  nclk=3-nclk
                  call system_clock(iclk(nclk), clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  net_wrt_clk=net_wrt_clk+inc_clk
                endif
#endif


              endif    !<-- part_type .eq./.gt. 0 switch
            endif    !<-- rec.eq.1 .or. series(i) switch

            !] END PART 9: READ/WRITE Partitioned array

          enddo    !<--- rec, loop over records

![
! Use intermediate one-per-record nf_sync of the target file only in
! verbose mode because it slows down the execution speed.   Note that
! time spent by nf_sync is counted as writing time, and so does time
! spent to close the target file. 

c--#ifdef VERBOSE
# if defined TIMING
          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate, clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          net_gray_clk=net_gray_clk+inc_clk
# endif
!          ierr=nf_sync (nctarg)

c          write(*,*) 'close' 
c          ierr=nf_close (nctarg) 
c          write(*,*) 'reopen'
c          ierr=nf_open (nctargname(1:ltrg), nf_write, nctarg)
# if defined TIMING
          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate, clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          net_sync_clk=net_sync_clk+inc_clk
# endif
c--#endif !]

        enddo    !<-- i=1,nvars, loop over variables

      endif ! <- complete (all the way at start of read/write parallel section)

        ![ PART10: CLOSE ALL FILES & DISPLAY TIMING SUMMARY:

        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false. 

! Close all files                        ! At this moment open/closed
                                         ! status of partial files
  98    if(my_rank==0) write(*,'(/1x,A)') 'closing files...'    ! depends on the state of CPP
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

      if (arg .lt. nargs)  goto 11 !--> next set of partial files.
                                   ! DevinD changed from 1 to 11 to avoid goto into master only.

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

#ifdef DEL_PART_FILES
        if(my_rank==0) then  ! Master only
          if (del_part_files) then
            write(*,'(1x,A,22x,F12.2,1x,A)')
     &      'Removing partial files :',net_rmcmd_clk/dble(clk_rate),'sec'
          endif
        endif
#endif
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

      !] END PART10: CLOSE ALL FILES & DISPLAY TIMING SUMMARY:

      call MPI_Finalize (ierr) ! cleanup MPI

      end


