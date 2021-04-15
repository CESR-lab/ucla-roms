
      program join_child_bry

      ![ INFO:
      !
      !  Tool to join output files of parent grid containing child bry files for
      !  input into child simulation.
      !
      !  Chose not to add to ncjoin_mpi because it differs too much and would pollute
      !  ncjoin_mpi with too many catches for joining child bry. Can still use several
      !  of the same routines from ncjoin_mod module.
      !
      !]

      ! PROGRAM SECTIONS:
      ! - PART 1: Confirm all input files are correct and related to each other
      ! - PART 2: Veryify all variables and dimensions in input files
      ! - PART 3: Create output file
      ! - PART X: Collectively (MPI-I/O) open all input and output files
      ! - PART X: Read/write partitioned variables   (bulk of work done here)
      ! - PART X: Close all files
      ! - PART X: Ddisplay program timing summary


![ CPP DEFS:
! Delete partial files after joining:
c#define DEL_PART_FILES
! Document program execution time:
#define TIMING
! Verbose terminal output:
c--#define VERBOSE
!]

      use ncjoin_mod

      implicit none


      deflate_level=0             ! leave as is - no compression and no shuffle.
      shuffle=.false.             ! just needed for compatibility with shared ncjoin_mod routines


      call init_timing_and_vars
      arg=0                       ! set to 0 for serial version as arguments of program start at 1, unlike ncjoin_mpi.

      my_rank=0                   ! allows easy use of ncjoin_mpi routines

      do while (arg .lt. nargs)

      nnodes=-1       ! DevinD repeated to get 'goto' outside of master only region, as called by all procs. ! used to be 11 marker here
      mayday=.false.  ! reset mayday flag

      if(my_rank==0) then  ! leave here incase make mpi version MPI MASTER ONLY: PRE-PROCESSING. i.e. check input files, creating output file, etc.
                           ! Extract a set of files which cover the whole physical grid.

        write(*,'(/1x,A/)') 'Pre-processing input files...'

        call check_partial_file_set  ! PART 1: CHECK ALL INPUT PARTIAL FILES ARE CORRECT
        if (mayday) goto 23          ! Using goto the break from if(my_rank==0)
                                     ! Only other idea I can think of is using select case (my_rank) case (0) ... as I can use exit with this but not for if.

      ![ PART 2: VERIFY ALL VARIABLES & DIMENSIONS:

#ifdef TIMING
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
#endif
        do node=0,nnodes-1
          lncn=lenstr(ncname(node))

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

        call determine_global_child_bry_dims

!        call determine_joined_file_dim_sizes
!        if (mayday) goto 97

!        call indentify_boundary_edges  ! not sure this is needed anymore

      !] END PART 2: VERIFY ALL VARIABLES & DIMENSIONS

        call create_joined_child_file                                  ! PART 3:
        if (mayday) goto 23

  23    if(my_rank==0) write(*,'(/1x,A/)')
     &                          'End of master proc pre-processing.'

      endif ! <- if(my_rank==0) END OF MASTER ONLY PRE-PROCESSING


      if (mayday==.true.) goto 97

!      call exchange_mpi_master_pre_processing_to_all                   ! PART 4:

!      call collectively_open_input_output_files                        ! PART 5:

        do i=1,nvars ! LOOP THROUGH VARIABLES & WRITE EACH ONE

          lvar=lenstr(vname(i))
          if(my_rank==0) write(*,'(16x,3A)')
     &      'Assembly of array name = ''', vname(i)(1:lvar), ''''

          do rec=1, tsize ! LOOP THROUGH TIMESTEPS OF CURRENT VAR

            if(my_rank==0) write(*,'(4x,A,I8,1x,A,I8,1x,A)')
     &       'Processing record', rec, 'out of', tsize,  '...'

            do node=0,nnodes-1
              call read_write_child_bry                                 ! PART X:
            enddo  !<-- nodes

          enddo    !<-- rec, loop over records

        enddo      !<-- i=1,nvars, loop over variables

        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false.  ! most errors arrive here

  98    call close_file_set                                            ! PART 10

      enddo !<- do while (arg .lt. nargs)


      call display_timing_summary


      end program join_child_bry


