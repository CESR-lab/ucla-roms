
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
! - PART10: Close all files
! - PART11: Ddisplay program timing summary

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

      implicit none

      ! remove netcdf.inc once fully converted to nf90_ routines
!#include "netcdf.inc"  ! added to ncjoin_mod
      ! Change this to 'use mpi' for module rather... as per 'simple_xy_par_wr2.f90'
!#include "mpif.h" ! added use mpi to ncjoin_mod

#include "ncjoin_mpi.opt"


      call init_timing_and_vars

      call setup_mpi

      do while (arg .lt. nargs)

      nnodes=-1       ! DevinD repeated to get 'goto' outside of master only region, as called by all procs. ! used to be 11 marker here
      mayday=.false.  ! reset mayday flag

      if(my_rank==0) then  ! MPI MASTER ONLY: PRE-PROCESSING. i.e. check input files, creating output file, etc.
                           ! Extract a set of files which cover the whole physical grid.

        write(*,'(/1x,A/)') 'Pre-processing input files...'

        call check_partial_file_set  ! PART 1: CHECK ALL INPUT PARTIAL FILES ARE CORRECT
        if (mayday) goto 23          ! Using goto the break from if(my_rank==0)
                                     ! Only other idea I can think of is using
                                     ! select case (my_rank)
                                     !   case (0) ... as I can use exit with this but not for if.

      ![ PART 2: VERIFY ALL VARIABLES & DIMENSIONS:

! At this moment a set of files recorded as ncname(0:nnodes-1),
! xi_start(0:nnodes-1), eta_start(0:nnodes-1) comprise a complete
! set, but the files are actually in "closed" state and all netCDF
! IDs are reset to -1.

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

            call check_ndims_nvars_natts(node)
            if (mayday) goto 23

            call create_catalog_of_var_names_IDs_ranks(node)
            if (mayday) goto 23

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

        call determine_joined_file_dim_sizes
        if (mayday) goto 97

        call indentify_boundary_edges  ! not sure this is needed anymore

      !] END PART 2: VERIFY ALL VARIABLES & DIMENSIONS

        call create_joined_chunked_file                                ! PART 3:
        if (mayday) goto 97

  23    if(my_rank==0) write(*,'(/1x,A/)')
     &                          'End of master proc pre-processing.'

      endif ! <- if(my_rank==0) END OF MASTER ONLY PRE-PROCESSING

      ![ PART 4: EXCHANGE ALL PARTIAL FILE INFO FROM MASTER TO ALL PROCS:

      call MPI_Barrier(MPI_COMM_WORLD, ierr)                             ! exchange vars needed
      call MPI_BCAST(mayday,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)  ! after if (mayday==.true.) goto 97
      call MPI_BCAST(complete, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(arg,      1, MPI_INT,     0, MPI_COMM_WORLD, ierr)  ! need for catch at end of program
      call MPI_BCAST(nargs,    1, MPI_INT,     0, MPI_COMM_WORLD, ierr)
      if (mayday) goto 97 !20231013 got rid of ==.true. due to compiler issues

      call exchange_mpi_master_pre_processing_to_all                   !]PART 4:

      call collectively_open_input_output_files                        ! PART 5:

        do i=1,nvars ! LOOP THROUGH VARIABLES & WRITE EACH ONE

          call assign_read_in_partial_files_to_mpi_proc                ! PART 6:

          do rec=1,tsize ! LOOP THROUGH TIMESTEPS OF CURRENT VAR

            if(my_rank==0) write(*,'(4x,A,I8,1x,A,I8,1x,A)')
     &       'Processing record', rec, 'out of', tsize,  '...'

            if (rec.eq.1 .or. series(i)) then

              if (part_type(i).eq.0 .and. .not.series(i)) then
                                                                       ! PART 7: ONLY DONE BY MASTER
                if (my_rank==0) call read_write_1D_var                 ! all procs needs to know if 'mayday' activated
                call MPI_BCAST(mayday, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
                if (mayday) goto 97


              elseif (part_type(i).eq.0) then
                                                                       ! PART 8:
                if (my_rank==0) call read_write_non_partitioned_var    ! all procs needs to know if 'mayday' activated
                call MPI_BCAST(mayday, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
                if (mayday) goto 97


              elseif (part_type(i).gt.0) then

                call read_write_partitioned_var                        ! PART 9:


              endif    !<-- part_type .eq./.gt. 0 switch

            endif    !<-- rec.eq.1 .or. series(i) switch

          enddo    !<--- rec, loop over records

        enddo    !<-- i=1,nvars, loop over variables


        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false.  ! most errors arrive here

  98    call close_file_set                                            ! PART 10

      enddo !<- do while (arg .lt. nargs) replaces ! if (arg .lt. nargs)  goto 11 !--> next set of partial files.
                                                   ! DevinD changed from 1 to 11 to avoid goto into master only.

      call display_timing_summary                                      ! PART 11

      call MPI_Finalize (ierr) ! cleanup MPI

      end


