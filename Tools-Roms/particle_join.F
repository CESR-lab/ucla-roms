
      program particle_join

      ![ INFO:
      !
      !  Tool to join output files of particle tracking module.
      !
      !  Long term fixes:
      !  - ....
      !
      !]

      ! PROGRAM SECTIONS:
      ! - PART 1: Confirm all input files are correct and related to each other
      ! - PART 2: Veryify all variables and dimensions in input files
      ! - PART 3:
      ! - PART X:


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


      deflate_level=1                ! leave as is - no compression and no shuffle.
      shuffle=.true.                 ! just needed for compatibility with shared ncjoin_mod routines


      call init_timing_and_vars
      arg=0                          ! set to 0 for serial version as arguments of program start at 1, unlike ncjoin_mpi.

      my_rank=0                      ! allows easy use of ncjoin_mpi routines, even though serial program

      call particle_check_args       ! NEW

      do while (arg .lt. nargs)

      nnodes=-1       ! DevinD repeated to get 'goto' outside of master only region, as called by all procs. ! used to be 11 marker here
      mayday=.false.  ! reset mayday flag

      if(my_rank==0) then  ! leave here incase make mpi version MPI MASTER ONLY: PRE-PROCESSING. i.e. check input files, creating output file, etc.
                           ! Extract a set of files which cover the whole physical grid.

        write(*,'(/1x,A/)') 'Pre-processing input files...'

        ! NEW: believe this is legal. Partition still exists in particle file. Some redundant allocation perhaps.
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

            ! NEW: not sure if should bother using this. Have so few dims in file...
            call check_ndims_nvars_natts(node)
            if (mayday) goto 23

!            call cbry_find_relevant_parent_nodes

            ! ONLY DO THIS CATALOGUE for relevant nodes!
            call create_catalog_of_var_names_IDs_ranks(node)
            if (mayday) goto 23

            ! NEW: don't need to do this anymore, might want to do it though for efficiency later.
            !      might not like loads of files open...
!            if (node.gt.0) then           ! close all the files, except for node=0.
!              ierr=nf_close(ncid(node))   ! since master only, need to close files to open collectively later.
!              ncid(node)=-1               ! keep node=0 open as still using below.
!            endif

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

!        call cbry_determine_global_dims

!        call determine_joined_file_dim_sizes
!        if (mayday) goto 97  ! NEW: do i still need this here?

!        call indentify_boundary_edges  ! not sure this is needed anymore

      !] END PART 2: VERIFY ALL VARIABLES & DIMENSIONS

        call particle_create_joined_file                                ! PART 3:
        if (mayday) goto 23

  23    if(my_rank==0) write(*,'(/1x,A/)')
     &                          'End of master proc pre-processing.'

      endif ! <- if(my_rank==0) END OF MASTER ONLY PRE-PROCESSING

      if (mayday==.true.) goto 97

!      call exchange_mpi_master_pre_processing_to_all                   ! PART 4:

!      call collectively_open_input_output_files                        ! PART 5:

        ! WRITE OCEAN_TIME:
        if(my_rank==0) write(*,'(16x,A)')
     &    'Assembly of ocean_time'
        node=0; i=1                                                     ! set variable to ocean_time (should be 1)
        vartype(i)=nf90_double  ! hack because vartype wasn't done in pre-proc for some reason
        do rec=1, tsize
          call read_write_non_partitioned_var                                         ! write ocean_time from 1 node only
        enddo

        ! WRITE PARTICLES:
        ! assumption is that individual writes to file are far worse than extra parsing of data.
        ! - do 1 sub-domain (node) at a time
        ! - read in entire particle array for sub-domain
        ! - parse through all tags of every step, and compile list of contained tags
        ! - populate one particle array by searching each time record
        ! - write that one particle, writing continuous time-records where possible.
        ! - one all tags complete move on to next sub-domain.

        ! allocate output array - can be reused
        ! do one particle at a time because otherwise output array is just too big
        if (allocated(p1_out)) deallocate(p1_out)          ! might not be big enough for next node, reallocate rather
        allocate( p1_out(6,tsize) )                        ! 1 particle

        i=2  ! NEW: hard-code particle var id
        do node=0,nnodes-1                                 ! loop over all input sub-domains

          if(my_rank==0) write(*,'(16x,A,I4)')
     &      'Assembly of particles in node =', node

          npart = dimsize(2,node)                          ! get number of particles in read in array
          print *, 'npart = ',npart

          if (allocated(p_in)) deallocate(p_in)            ! array for read in particles. might not be big enough for next node, reallocate rather
          allocate( p_in(npart,7,tsize) )

          ! read in full array (may want to do per t-step for memory efficiency later:
!          ierr=nf90_inq_varid( ncid(node), 'particles',vid(i,node) )
!          if(ierr/=0) print *, 'BAD ID', nf90_strerror(ierr), 'node', node
          ierr=nf90_get_var(ncid(node), vid(i,node), p_in )
          if(ierr/=0) print *, 'BAD READ', nf90_strerror(ierr), 'node', node
!          print *, 'p_in(1,:,1)=',p_in(1,:,2)
!          print *, 'p_in(npart-10,:,1)=',p_in(npart-10,:,2)  ! hoping for nan or fill, yes = E+36!
!          stop 'STOPPED'

          ! find ptag's contained in node
          if (allocated(ptags_in)) deallocate(ptags_in)       ! array of ptags contained in node
          allocate( ptags_in(npart) )
          cnt_ptag=0                                          ! set ptag counter
          do rec=1,tsize
            do p=1,npart
!              print *, 'p_in(p,1,rec)=',p_in(p,1,rec)
              if (p_in(p,1,rec) < 1.0e30) then                ! valid tag (indx 1) - less than fill value
                found_tag=.false.                             ! only reset if valid tag
                if (cnt_ptag > 0) then                        ! first tag so no checking for conflict required                                         ! check if haven't found tag before (unique)
                  do it=1,cnt_ptag
                    if (ptags_in(it) == p_in(p,1,rec)) then   ! found it - go on to next particle
!                      print *, 'found particle - cycle to next'
                      found_tag = .true.                       ! need flag to cycle p loop
                      cycle
                    endif
                  enddo
                  if (found_tag) cycle
                endif
                cnt_ptag=cnt_ptag+1                           ! haven't found - we have a new particle
!                print *, 'counting new! cnt_ptag=',cnt_ptag
                ptags_in(cnt_ptag) = p_in(p,1,rec)            ! store new particle tag in list of known particles
              else
                exit                                          ! extent of particles for this record,
              endif                                           ! exit particle loop to next record?
            enddo  ! nparts
            print *, 'rec=',rec,'cnt_ptag=',cnt_ptag
          enddo    ! recs

!          stop 'stopped'


          ! output ptag records (doing continguous times in one go for efficiency of write calls (assumed bottle neck) )
          do tag=1,cnt_ptag                 ! all unique tags
            p1_out = -1e20                  ! catch non-continuous times
            do rec=1,tsize
              do p=1,npart
                if(p_in(p,1,rec) < 1.0e20) then               ! valid tag - less than fill value

                  if (p_in(p,1,rec) == ptags_in(tag)) then
                    p1_out(:,rec) = p_in(p,2:7,rec)           ! store 6 details of particle
                    exit                                      ! next time record as partilce found
                  endif
                else
                  exit                                        ! extent of particles for this record, particle not present in this timestep
                endif                                         ! exit particle loop to next record?

              enddo  ! nparts
            enddo    ! recs

            ! need to write tag - loop through concurrent times to reduce writes:
            cnt_rec=0
            do rec=1,tsize
              if (p1_out(1,rec) /= -1e20) then
                if (cnt_rec==0) rec_str=rec
                cnt_rec=cnt_rec+1
              else
                if (cnt_rec>0) then                      ! catch incase nothing found yet
!                  start_out =  ! set this up
!                  count_out =  ! set this up
                  ierr=nf90_put_var(nctarg, varid(i), p1_out(:,rec_str:rec-1), start_out, count_out)
                  if(ierr/=nf90_noerr) print *, 'BAD WRITE!', nf90_strerror(ierr),'my_rank', my_rank
                  cnt_rec=0
                endif
              endif
            enddo
          enddo      ! tags

        enddo  !<-- nodes

        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false.  ! most errors arrive here

  98    call close_file_set                                            ! PART 10

      enddo !<- do while (arg .lt. nargs)


      call display_timing_summary


      end program particle_join
