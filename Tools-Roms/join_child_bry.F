      ![ INFO:
      !
      !  Tool to join output files of parent grid containing child bry files for
      !  input into child simulation.
      !
      !  Chose not to add to ncjoin_mpi because it differs too much and would pollute
      !  ncjoin_mpi with too many catches for joining child bry. Can still use several
      !  of the same routines from ncjoin_mod module.
      !
      !  Added cbry_ to beginning of routines that are specifically for child boundary,
      !  not ncjoin_mpi.
      !
      !  Long term fixes:
      !  - really online bry shouldn't write empty files, and therefore this shouldn't
      !    need to deal with empty files.
      !    Online files should contain list of node numbers to join, and use that rather
      !    than pre-process.
      !  - Could parallelize in future.
      !
      !]

      ! PROGRAM SECTIONS:
      ! - PART 1: Confirm all input files are correct and related to each other
      ! - PART 2: Veryify all variables and dimensions in input files
      ! - PART 3: Create output file
      ! - PART X: Collectively (MPI-I/O) open all input and output files
      ! - PART X: Read/write partitioned variables   (bulk of work done here)
      ! - PART X: Close all files
      ! - PART X: Display program timing summary


![ CPP DEFS:
! Delete partial files after joining:
c#define DEL_PART_FILES
! Document program execution time:
#define TIMING
! Verbose terminal output:
c--#define VERBOSE
!]

      module join_child_bry_mod

      use ncjoin_mod

      implicit none

      ! ------ CHILD_BRY VARS -------------
      !  Added cbry_ to beginning of routines that are specifically for child boundary,
      !  not ncjoin_mpi.
      integer :: global_south_bry_size    ! size of join dimension in child boundary file
      integer :: global_west_bry_size
      integer :: tot_pwithc               ! total parent domains with child bry
      integer :: max_cvars_p              ! max number of child variables contained in sub-domain
      integer :: nd, vr                   ! loop index
      logical :: new_var = .true.         ! trigger reading start and count for new variable
      integer, dimension(4) :: vid_out    ! variable ID's in joined output file
      integer, dimension(3) :: child_start, child_count, startin, countin  ! start+counts for read/write
      integer, dimension(:),   allocatable :: tot_cvars_p_all, tot_cvars_p ! child vars in parent (all nodes and actual)
      integer, dimension(:),   allocatable :: pwithc_nodes                 ! list of whether parent contains child
      integer, dimension(:,:), allocatable :: cvars_p_all, cvars_p         ! each child var number in parent
      real,    dimension(:,:), allocatable :: bry_data

      contains

!----------------------------------------------------------------------
      subroutine cbry_check_args  ![

      ! Ensure reasonable number of arguements used,
      ! if not give usage instructions.

      implicit none

      if(nargs < 4) then                                               ! too few program arguments
        write(*,'(/1x,A//10x,A)') 'Correct usage of join_child_bry:',
     &    'join_child_bry      child_bry.0000.*.nc'
        write(*,'(1x,A/10x,A/)') 'or',
     &    'join_child_bry  -d  child_bry.0000.*.nc'
        write (*,'(/1x,A//1x,A/1x,A/)')
     &    'where -d deletes the partial files.',
     &    'This tool joins partial files of online interpolated child',
     &    'boundary files created by a ROMS parent simulation.'


        error stop
      endif

      end subroutine cbry_check_args  !]

! ---------------------------------------------------------------------
      subroutine cbry_determine_global_dims      ![

      ! Determine the size of global child dimensions
      ! in order to create joined child file

      implicit none

      ! local
      integer :: dimid

      ierr=nf90_get_att(ncid(0), NF90_GLOBAL, 'global_south_size', global_south_bry_size)
      ierr=nf90_get_att(ncid(0), NF90_GLOBAL, 'global_west_size', global_west_bry_size)

      ierr=nf90_inq_dimid( ncid(0), 's_rho', dimid)  ! get s_rho veritcal levels
      ierr=nf90_inquire_dimension(ncid(0), dimid, len=s_rho)

      end subroutine cbry_determine_global_dims  !]

! ----------------------------------------------------------------------
      subroutine cbry_create_joined_file  ![
      implicit none  ! create joined child bry file

      call create_joined_empty_file
      if(mayday) return

      call copy_global_attributes
      if(mayday) return

      ! Define dimensions
      ! input file dimnames are different from output names. E.g. xi_rho_south in input
      ! and just xi_rho in output file. Therefore cannot use ncjoin method of dimname(i)(1:lvar)

      ierr=nf_def_dim (nctarg, 'xi_rho',
     &                 global_south_bry_size,   dimid(1))
      ierr=nf_def_dim (nctarg, 'eta_rho',
     &                 global_west_bry_size,    dimid(2))
      ierr=nf_def_dim (nctarg, 's_rho',
     &                 s_rho,  dimid(3))
      ierr=nf_def_dim (nctarg, 'rec_time', 0,   dimid(4))  ! 0 for unlimited
      ierr=nf_def_dim (nctarg, 'xi_u',
     &                 global_south_bry_size-1, dimid(5))
      ierr=nf_def_dim (nctarg, 'eta_v',
     &                 global_west_bry_size -1, dimid(6))

      ! Define variables & atts

      call cbry_copy_vars_atts

      ! Leave definition mode

!      ierr=nf90_close(ncid(0)) ! close last input file. completed master only reads
!      ncid(0)=-1

      ierr=nf_enddef (nctarg)
!      ierr=nf90_close(nctarg)  ! no need to close since serial

# ifdef VERBOSE
      if(my_rank==0) write(*,'(/1x,A)') 'Leaving definition mode.'
# endif

#ifdef TIMING
      nclk=3-nclk
      call system_clock (iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      net_fcrt_clk=net_fcrt_clk+inc_clk  ! Timing for file creation
#endif

      end subroutine cbry_create_joined_file  !]

! ----------------------------------------------------------------------
      subroutine cbry_read_write  ![
      ! due to cbry_find_relevant_parent routine, we know this
      ! variable is contained in parent and msut be written
      implicit none

# ifdef TIMING
      nclk=3-nclk
      call system_clock (iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      net_gray_clk=net_gray_clk+inc_clk
# endif


      if (new_var) then                                                          ! once per variable
        startin=1; countin=1; child_start=1; child_count=1

        if (ncid(node).eq.-1)
     &    ierr=nf90_open (ncname(node), nf90_nowrite, ncid(node))

        if (vdims(i) >= 2) then  ! exclude time (1D vars)
          ierr = nf90_get_att(ncid(node), vid(i,node), 'start', child_start(1))  ! start position in global child bry
          ierr = nf90_get_att(ncid(node), vid(i,node), 'count', child_count(1))  ! count of sub-domain entry in global child bry
        endif
        countin(1:2)=(/ child_count(1), nmax(i) /)

!        print *, 'ALLLLLOC BRY_DATA'
        if (allocated(bry_data)) deallocate (bry_data)                           ! variable for reading data
        allocate ( bry_data( child_count(1), nmax(i) ) )
        new_var = .false.
!        print *, 'DONE ALLLLLOC'
!        print *, 'child_count(1)',child_count(1)
!        print *, 'nmax(i)',nmax(i)
      endif  !<-- new_var

      if     (vdims(i) == 1) then                 ! set new record
        startin(1)     = rec
        child_start(1) = rec
        countin(1)=1
      elseif (vdims(i) == 2) then                 ! handle 2D or 3D variable dims
        startin(2)     = rec                      ! set recorded timestep for read in partial bry
        child_start(2) = rec                      ! set recorded timestep for written join bry
      elseif (vdims(i) == 3) then
        startin(3)     = rec                      ! set recorded timestep for read in partial bry
        child_start(3) = rec                      ! set recorded timestep for written join bry
        child_count(2) = s_rho                    ! should be nmax rather?
      else
        error stop 'ERROR: dimensions issue when writing!'
      endif

      ! read variable:

      ierr = nf90_get_var( ncid(node), vid(i,node), bry_data,
     &                     start=startin(1:vdims(i)), count=countin(1:vdims(i)) )
      if(ierr/=0) then
        print *, 'read ierr=', nf90_strerror(ierr), 'node=', node
        error stop 'READING ERROR'
      endif

# ifdef TIMING
      net_read_size=net_read_size+size
      nclk=3-nclk
      call system_clock (iclk(nclk),clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      net_read_clk=net_read_clk+inc_clk

!        nclk=3-nclk  ! DevinD moved here
!        call system_clock (iclk(nclk), clk_rate, clk_max)
!        inc_clk=iclk(nclk)-iclk(3-nclk)
!        net_gray_clk=net_gray_clk+inc_clk
# endif

      ! write variable:

      ierr = nf90_put_var( nctarg, varid(i), bry_data,
     &                     start=child_start(1:vdims(i)), count=child_count(1:vdims(i)) )
      if(ierr/=0) then
        print *, 'wrt ierr=', nf90_strerror(ierr), 'node=', node
        error stop 'WRITING ERROR'
      endif

#ifdef TIMING
      net_wrt_size=net_wrt_size+size
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      net_wrt_clk=net_wrt_clk+inc_clk
#endif

      ! Need to close and sync for time being.
!      ierr=nf90_sync (nctarg) ! debug only

      end subroutine cbry_read_write  !]

! ----------------------------------------------------------------------
      subroutine cbry_copy_vars_atts  ![
      implicit none                   !  copy the variables and attributes for join_child_bry.
                                      !  can't use ncjoin_mpi routine because of chunking and dimension issues,
                                      !  since child bry sizes are different on east and west as orientation different from parent.

#ifdef VERBOSE
        if(my_rank==0) write(*,'(1x,2A)') 'Variables, their dimensions and ',
     &                                           'attributes:'
#endif

        do i=1,nvars ! LOOP THROUGH VARS

          node=vnode(i)
          lncn=lenstr(ncname(node))
          if (ncid(node).eq.-1) ierr=nf_open  (ncname(node),
     &                               nf_nowrite, ncid(node))
          if (ierr .eq. nf_noerr) then
            ierr=nf_inq_var (ncid(node), vid(i,node), vname(i),
     &                vartype(i), vdims(i), dimids(1,i),  varatts)

            if (ierr .eq. nf_noerr) then

              if     (vdims(i) == 1) then
                dimids(1,i) = dimid(4)
                nmax(i)=1
              elseif (vdims(i) == 2) then                         ! handle 2D or 3D variable dimids
                dimids(2:vdims(i),i) = (/ dimid(4) /)             ! same for r/u/v vars
                nmax(i)=1                                         ! record vertical levels of var
              elseif (vdims(i) == 3) then
                dimids(2:vdims(i),i) = (/ dimid(3), dimid(4) /)
                nmax(i)=s_rho
              else
                error stop 'ERROR: dimensions when creating file!'
              endif

              if (vdims(i) >= 2) then                       ! not 1D time variable
                lvar=lenstr(vname(i))
                if (vname(i)(lvar-4:lvar)=='south'  .or.
     &              vname(i)(lvar-4:lvar)=='north') then

                  dimids(1,i) = dimid(1)                      ! xi_rho for rho and v points
                  if (vname(i)(1:2)=='u_'     .or.
     &                vname(i)(1:5)=='ubar_') then
                    dimids(1,i) = dimid(5)                    ! xi_u for u points
                  endif

                else                                          ! else west or east bry

                  dimids(1,i) = dimid(2)                      ! eta_rho for rho and u points
                  if (vname(i)(1:2)=='v_'     .or.
     &                vname(i)(1:5)=='vbar_') then
                    dimids(1,i) = dimid(6)                    ! eta_v for v points
                  endif

                endif
              endif    ! <-- vdims >= 2

              ierr=nf90_def_var (nctarg, vname(i)(1:lvar),vartype(i),  ! add compression here if desired
     &                           dimids(1:vdims(i),i), varid(i))

              if (ierr .eq. nf_noerr) then
!#ifdef VERBOSE
!                write(*,'(8x,3A,8I3)')   '''', vname(i)(1:lvar),
!     &              ''',  dimids =', (dimids(j,i), j=1,vdims(i))
!#endif
                do j=1,varatts                                         ! copy attributes
                  ierr=nf_inq_attname (ncid(node), vid(i,node),
     &                                              j, string)
                  if (ierr .eq. nf_noerr) then
                    lstr=lenstr(string)

                    if (string /= 'start' .and.                ! ignore 'start' and 'count' attributes, not relevant to joined file.
     &                  string /= 'count'      ) then

                      ierr=nf_copy_att (ncid(node), vid(i,node),
     &                         string(1:lstr), nctarg, varid(i))
                      if (ierr. ne. nf_noerr) then
                        write(*,'(/1x,2A,I3,3A/12x,4A)')   '### ERROR: ',
     &                   'Cannot copy attribute #', j,' for variable ''',
     &                    vname(i)(1:lvar),  ''' into netCDF', 'file ''',
     &                    nctargname(1:ltrg), '''.  ', nf_strerror(ierr)
                        mayday=.true.; return ! goto 97
                      endif
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

!            series(i)=.false.
!            part_type(i)=0
!            do j=1,vdims(i)
!              if (dimids(j,i).eq.id_xi_rho .or.
!     &            dimids(j,i).eq.id_xi_u) then
!                part_type(i)=part_type(i)+1
!              elseif (dimids(j,i).eq.id_eta_rho .or.
!     &                dimids(j,i).eq.id_eta_v) then
!                part_type(i)=part_type(i)+2
!              elseif (dimids(j,i).eq.unlimdimid) then
!                series(i)=.true.
!              endif
!            enddo

!            if (node.gt.0) then
!              ierr=nf90_close(ncid(node)) ! DevinD: Most vars in node==0, hence keep
!              ncid(node)=-1               ! it open whilst looping through vars, close
!            endif                         ! any node>0 as might not be used again.

          else
            write(*,'(/1x,A,1x,3A/12x,A)')  '### ERROR: Cannot open ',
     &                  'netCDF file ''', ncname(node)(1:lncn), '''.',
     &                                           nf90_strerror(ierr)
            mayday=.true.; return ! goto 97
          endif

        enddo  ! <-- i=1,nvars, variable IDs.


        do i=1,ndims                  ! set time record size
          if (i.eq. unlimdimid) then
            tsize=dimsize(i,0)
          endif
        enddo

      end subroutine cbry_copy_vars_atts  !]

! ----------------------------------------------------------------------
      subroutine cbry_find_relevant_parent_nodes  ![
      ! find all parent nodes that actually contain child boundaries
      ! to limit loops to only nodes containing child data

      ! If count>=1 for any variable, then exit loop, bool true, add node number to array.
      !

      implicit none

      integer :: c, nd, vrs
      logical :: vfirst

      if (node==0) then
        allocate( tot_cvars_p_all( 0:nnodes-1) )           ! number of child variables in parent
        allocate(     cvars_p_all( 0:nnodes-1, ibuff(2) )) ! each child var number in parent, ibuff(2) is total vars
        tot_cvars_p_all=0
        cvars_p_all=0
      endif

      vfirst = .true.
      vrs=0                                                ! reset
      do i=1,ibuff(2)                                      ! nvars, which also equates to all the variable ID's (clearly range 1:nvars)
        c=0
        ierr = nf90_get_att(ncid(node), i, 'count', c)     ! count of sub-domain entry in global child bry
        if (ierr /= 0) c=0                                 ! if didn't read correctly, had issues so reset c
        if (c>0) then                                      ! var contains child bry
          if(vfirst) then
            tot_pwithc = tot_pwithc+1                      ! count total number of parent with child
            vfirst = .false.
          endif
          vrs = vrs+1
          cvars_p_all(node,i) = i                          ! store child var number to pnode
        endif
      enddo  !<-- nvars

      tot_cvars_p_all(node) = vrs
      if (vrs > max_cvars_p) max_cvars_p = vrs             ! need to establish parent with most child vars for array

      if (node==nnodes-1) then                             ! consolodate array list of parents containing child
        allocate( pwithc_nodes( tot_pwithc ))
        allocate(  tot_cvars_p( tot_pwithc ))
        allocate( cvars_p( tot_pwithc, max_cvars_p ))      ! has to be as long as the parent with most
        c=0                                                ! child vars to share array
        do nd=0,nnodes-1                                   ! all nodes to add to final node list, only know size now
          if (tot_cvars_p_all( nd )>0) then
            c=c+1
            pwithc_nodes(c) = nd
            tot_cvars_p(c) = tot_cvars_p_all(nd)
            vrs=0
            do i=1,ibuff(2)
              if (cvars_p_all(nd,i) > 0) then              ! parent contains child var
                vrs=vrs+1
                cvars_p( c, vrs ) = i                      ! store variable number
              endif
            enddo
          endif
        enddo    !<-- nodes

!        print *, 'pwithc_nodes=', pwithc_nodes
!        print *, 'tot_cvars_p=', tot_cvars_p
!        do nd=1,tot_pwithc
!          print *, 'cvars_p(nd)=', cvars_p(nd,:)
!        enddo

        deallocate( tot_cvars_p_all )
        deallocate( cvars_p_all     )

      endif

      end subroutine cbry_find_relevant_parent_nodes  !]

! ----------------------------------------------------------------------
!      subroutine cbry_read_write_ocean_time  ![
!      ! write ocean_time outside node loop, nice all nodes contain ocean_time and
!      ! easier to have own routine than put catches in the main cbry_read_write routine
!
!      implicit none
!
!
!      ! check first variable is called ocean_time
!
!      end subroutine cbry_read_write_ocean_time  !]

!----------------------------------------------------------------------

      end module join_child_bry_mod

!----------------------------------------------------------------------

      program join_child_bry

      use ncjoin_mod
      use join_child_bry_mod

      implicit none


      deflate_level=0                                      ! leave as is - no compression and no shuffle.
      shuffle=.false.                                      ! just needed for compatibility with shared ncjoin_mod routines


      call init_timing_and_vars
      arg=0                                                ! set to 0 for serial version as arguments of program start at 1, unlike ncjoin_mpi.

      my_rank=0                                            ! allows easy use of ncjoin_mpi routines

      call cbry_check_args

      do while (arg .lt. nargs)

      nnodes=-1                                            ! repeated to get 'goto' outside of master only region, as called by all procs. ! used to be 11 marker here
      mayday=.false.                                       ! reset mayday flag

      if(my_rank==0) then                                  ! leave here incase make mpi version MPI MASTER ONLY: PRE-PROCESSING. i.e. check input files, creating output file, etc.
                                                           ! Extract a set of files which cover the whole physical grid.

        write(*,'(/1x,A/)') 'Pre-processing input files...'

        call check_partial_file_set                        ! PART 1: CHECK ALL INPUT PARTIAL FILES ARE CORRECT
        if (mayday) goto 23                                ! Using goto the break from if(my_rank==0)
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

            ! use check_ndims... above to get number of variables etc.
            ! make a list of nodes that actually contain child bry
            ! this speeds up program to avoid reading empty nodes every time
            ! once roms is corrected to only produce files at specific partitions
            ! then remove this routine. The buff(2) partition attribute will need
            ! to be fixed to match only number of parents with children.
            call cbry_find_relevant_parent_nodes

            ! ONLY DO THIS CATALOGUE for relevant nodes!
            call create_catalog_of_var_names_IDs_ranks(node)
            if (mayday) goto 23

            if (node.gt.0) then                            ! close all the files, except for node=0.
              ierr=nf_close(ncid(node))                    ! since master only, need to close files to open collectively later.
              ncid(node)=-1                                ! keep node=0 open as still using below.
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

        call cbry_determine_global_dims

!        call determine_joined_file_dim_sizes
!        if (mayday) goto 97

!        call indentify_boundary_edges  ! not sure this is needed anymore

      !] END PART 2: VERIFY ALL VARIABLES & DIMENSIONS

        call cbry_create_joined_file                                  ! PART 3:
        if (mayday) goto 23

  23    if(my_rank==0) write(*,'(/1x,A/)')
     &                          'End of master proc pre-processing.'

      endif ! <- if(my_rank==0) END OF MASTER ONLY PRE-PROCESSING


      if (mayday==.true.) goto 97

!      call exchange_mpi_master_pre_processing_to_all      ! PART 4:

!      call collectively_open_input_output_files           ! PART 5:

        if(my_rank==0) write(*,'(16x,A)')
     &    'Assembly of ocean_time'
        node=0; i=1                                        ! set variable to ocean_time (should be 1)
        do rec=1, tsize
          call cbry_read_write                             ! write ocean_time from 1 node only
        enddo

        do nd=1,tot_pwithc                                 ! loop total parents with child brys
          node = pwithc_nodes(nd)                          ! set parent node number from stored list.

          if(my_rank==0) write(*,'(16x,A,I4)')
     &      'Assembly of bry vars in node =', node
          if(my_rank==0) write(*,'(18x,A,I4)')
     &      'bry vars to copy =', tot_cvars_p(nd)

          do vr=1,tot_cvars_p(nd)                          ! LOOP THROUGH VARIABLES & WRITE EACH ONE (only the variables within parent)
            i = cvars_p(nd,vr)                             ! set actual variable number in terms of all variables
            new_var = .true.

            do rec=1, tsize                                ! LOOP THROUGH TIMESTEPS OF CURRENT VAR

              call cbry_read_write                         ! PART X:

            enddo  !<-- rec, loop over records
          enddo    !<-- i=1,nvars, loop over variables

        enddo      !<-- nodes

        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false.                                  ! most errors arrive here

  98    call close_file_set                                ! PART 10

      enddo        !<- do while (arg .lt. nargs)


      call display_timing_summary


      end program join_child_bry


