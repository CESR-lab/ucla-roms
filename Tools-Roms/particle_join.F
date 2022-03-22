
      ! INFO:  ![
      !
      !  Tool to join output files of particle tracking module.
      !
      !  Inital tool developement - Devin Dollery & Jereon Molemaker - 2022 March
      !
      !  The program is called 'particle_join'. It gets most of its variables
      !  and subroutines from the particle_join_mod module. It uses some
      !  utility joining routines from the ncjoin_mod module.
      !
      !  METHOD:
      !  - PART 1: Confirm all input files are correct and related to each other (same as ncjoin)
      !  - PART 2: Veryify all variables and dimensions in input files (same as ncjoin)
      !  - PART 3: Create list of unique particles:
      !              for A) each sub-domain, and B) entire global domain
      !              by reading and parsing all input partitioned files.
      !              (Due to exchanges between sub-domains there will be some repeated particles in
      !               adjacent sub-domains. Don't double count in global unique list.)
      !  - PART 4: Create joined output file with known number of global unique particles.
      !              Use a fixed (not unlimited) time dimension, since the time records can be written
      !              out of order. E.g. Particle moving westward, recs 1 & 2 in node1, recs 3 & 4 in
      !              node0. But node0 will be written first.
      !              Also, fixed time-dimension is better for parallel writes, if we do that in future.
      !  - PART 5: Write particles one sub-domain at a time:
      !              Read in sub-domain data.
      !              Transfer all (non-ordered) particle input data into particle sorted output array.
      !              Write one particle at a time for continuous timesteps (since particle may not
      !              exist in sub-domain for all timesteps due to movement).
      !
      !  - PART 6: Compress file?
      !
      !  CONSIDERATIONS:
      !
      !    Joining particles requires many small writes to the output file.
      !    Compression for small writes is very inefficient, since the chunk
      !    of data needs to be uncompressed each time it is written to.
      !    Hence, don't use compression, rather compress all at the end if desired.
      !
      !    Assumption is that individual writes to file are far worse than extra parsing of data.
      !    We realized writing with compression is way too slow.
      !    This means we can't have a max_number multiple sized file and rely on compression.
      !    Rather, we will first parse all sub-domains and create a list of unique particle numbers
      !    sort this list as you go for quicker searching if contained already.
      !    Then can create output file with exact total number of particles.
      !    Write all particles.
      !    Compress all at very end if even desired.
      !
      ! MISC:
      !
      !    large number of particle tags represented as real's requires double precision
      !    real's to capture enough signficant figures. This is done by the 'mpc' program,
      !    probably better to do this by a compiler directive in future.
      !
      ! Long term fixes:
      !  - ....
      !]

![ CPP DEFS:
! Note: these cppdefs work for both module particle_join_mod and program
!       particle_join, because they are in the same file. Also, make
!       sure these are consistent with the cppdefs set in ncjoin_mod.
! Document program execution time:
#define TIMING
! Verbose terminal output:
c--#define VERBOSE
!]

      module particle_join_mod

      use ncjoin_mod
      use mpi
      use netcdf
      implicit none

      ! PARTICLES:
      integer :: p, it, glob_n_ptags, tag, cnt_rec, rec_str
      integer :: npart                                               ! number of particles in read in node
      logical :: found_tag
      real    :: no_rec_val = -1.0e-20

      integer,allocatable,dimension(:)  :: loc_n_ptags               ! local number of unique particles
      real,allocatable,dimension(:,:,:) :: p_in                      ! particles read in from 1 node. May change to one timestep
      real,allocatable,dimension(:,:,:) :: p_out                     ! output one particle at a time, but all records, for fewer writes
      real,allocatable,dimension(:)     :: ptags_in                  ! list of ptag contained in node
      real,allocatable,dimension(:)     :: glob_ptags                ! list of all unique ptags globally

      ! use a derived type for sub-domain arrays with list of unique particle tags.
      ! this way you don't have to allocate a 2D array big enough to capture sub-domain
      ! list of longest array rather only allocate array large enough for each sub-domain.
      type arrays_1d
        real,allocatable,dimension(:) :: arr
      end type arrays_1d
      type(arrays_1d), allocatable :: loc_ptags(:)

      contains

! ---------------------------------------------------------------------
      subroutine particle_check_args  ![
      ! ensure reasonable number of arguements used,
      ! if not give usage instructions.
      implicit none

      if(nargs < 4) then                                               ! too few program arguments
        write(*,'(/1x,A//10x,A)') 'Correct usage of particle_join:',
     &    'particle_join      sample_part.00000.*.nc'
        write(*,'(1x,A/10x,A/)') 'or',
     &    'particle_join  -d  sample_part.00000.*.nc'
        write (*,'(/1x,A//1x,A/1x,A/)')
     &    'where -d deletes the partial files.',
     &    'This tool joins partial files of particle tracking',
     &    'result files created by a ROMS simulation.'

        error stop
      endif

      end subroutine particle_check_args  !]

! ----------------------------------------------------------------------
      subroutine particle_list_unique_particles  ![
      ! create sub-domain list and global list of unique particles,
      ! can then create joined file with this.
      !
      ![assumption is that individual writes to file are far worse than extra parsing of data.
      ! read in entire particle array for sub-domain.
      ! parse through all tags of every step, and compile list of unique tags in sub-domain.
      ! sort the list by tag value, and store local list for sub-domain for use again later.
      ! use local list to compile global list of unique particles (also sorted).
      ! can create joined file with known number of unique particles.
      !
      ! global list of unique particles, the final size is only known at the end.
      ! to prevent allocations and copies of tmp_array to global_array for each
      ! sub-domain, we rather allocate a much larger tmp_array. Once it starts
      ! to get full, we reallocate more space.
      !]
      implicit none

      ! local
      integer                       :: cnt_p, ndup, low_dupl, buff_size
      logical                       :: duplicate
      real,dimension(:),allocatable :: tmp_ptags

      tsize=dimsize(unlimdimid,0)                          ! set here as set to 0 somewhere else for unlim size, but we don't want that
      print *, 'tsize = ', tsize

      allocate(   loc_ptags(0:nnodes-1) )                  ! list of local (sub-domain) unique particles
      allocate( loc_n_ptags(0:nnodes-1) )                  ! local number of unique particles

      glob_n_ptags=0                                       ! set ptag counter for global unique particles

      i=2                                                  ! hard-code particle var id for time-being
      do node=0,nnodes-1                                   ! loop over all input sub-domains

        if(my_rank==0) write(*,'(16x,A,I4)')
     &    'Pre-processing of particles in node =', node

        npart = dimsize(2,node)                            ! get size of particle array in read in array
!        print *, 'npart = ',npart                         ! not actual number of particles as it includes
                                                           ! extra space potential exchanged/seeded particles
        if (allocated(tmp_ptags)) deallocate( tmp_ptags )  ! tmp array of local (sub-domain) unique particles
        allocate( tmp_ptags( npart ) )
        if (allocated(p_in)) deallocate(p_in)              ! array for read in particles. might not be big enough
        allocate( p_in(npart,7,tsize) )                    ! for next node, need to reallocate
        cnt_p=0

        ierr=nf90_get_var(ncid(node), vid(i,node), p_in )  ! read in full array (may want to do per t-step for memory efficiency in future):
        if(ierr/=0) print *, 'BAD READ unique list', nf90_strerror(ierr), 'node', node

        ! FIND UNIQUE PTAG'S CONTAINED IN NODE:
        do rec=1,tsize
          do p=1,npart

            if (p_in(p,1,rec) < 1.0e30) then               ! valid tag (tag indx=1) < less than fill value

              duplicate=.false.
              if (cnt_p > 0) then                          ! no checking needed for first tag for pre-existing required
                call lookup_duplicate( tmp_ptags, p_in(p,1,rec), cnt_p, duplicate )  ! will return duplicate=true if duplicate
              endif
              if (.not. duplicate) then
                cnt_p=cnt_p+1                              ! new unique particle
                tmp_ptags(cnt_p) = p_in(p,1,rec)           ! store new particle tag in list of known particles
                call insertion_sort_1p( tmp_ptags, cnt_p ) ! sort new particle in list of sub-domain unique particles
              endif

            else
              exit                                         ! extent of particles for this record
            endif
          enddo    ! nparts
          print *, 'rec=',rec,'cnt_p=',cnt_p
        enddo      ! recs

        loc_n_ptags(node)=cnt_p                            ! store number of unique tags for sub-domain

        if (cnt_p>0) then                                  ! cnt_p>0 only if particles found in sub-domain

          allocate( loc_ptags(node)%arr( cnt_p ) )
          loc_ptags(node)%arr = tmp_ptags(1:cnt_p)         ! store list of unique tags for sub-domain

          ! ADD NEW PARTICLES TO GLOBAL UNIQUE TAG LIST:
          if (glob_n_ptags==0) then                        ! first sub-domain containing particles
            glob_n_ptags=cnt_p                             ! global count of unique particles
            buff_size=cnt_p * (nnodes/2)                   ! chosen to multiply unique particles by no. sub-domains/2
            allocate(glob_ptags(buff_size))                ! allocate with plenty extra space to avoid reallocation.
            glob_ptags(1:cnt_p)=loc_ptags(node)%arr        ! copy first unique sub-domain list of particles (already sorted)
          else                                             ! all subsequent sub-domains containing particles
            call check_buff_space_and_expand( glob_n_ptags+cnt_p, buff_size, glob_ptags )
                                                           ! need buff space as don't know how many particles duplicated, thus
            ndup=0                                         ! final global array size unknown.
            low_dupl=1                                     ! set low limit lookup search index
            do it=1,cnt_p                                  ! look in list of tags for all previous sub-domains for duplicate

              ! now since sub-domain list and global list are both sorted, don't have to search
              ! entire global list, once you know index of previous tag in global list, that
              ! becomes the lower limit. If duplicate, low is updated in lookup_duplicate.
              call lookup_duplicate( glob_ptags, loc_ptags(node)%arr(it), glob_n_ptags, duplicate, low_dupl )
              if (duplicate) then
                ndup=ndup+1                                ! count duplicate
              else
                glob_ptags(glob_n_ptags+it-ndup) = loc_ptags(node)%arr(it)  ! add non-duplicate to tmp global array
                ! starts from end so not sure how to improve search range? Ask JM?
                ! set new low search bound as well. Since copy array, need it's own low_sort
                call insertion_sort_1p( glob_ptags(1:glob_n_ptags+it-ndup), glob_n_ptags+it-ndup )
              endif
            enddo
            glob_n_ptags=glob_n_ptags+cnt_p-ndup           ! global count of unique particles (plus sub-domain count less duplicates)
          endif

          print *, 'glob_n_ptags=',glob_n_ptags            ! debug
!          print *, 'glob_ptags=',glob_ptags;print *, ' ';print *, ' '

        endif      ! <- cnt_p>0
      enddo        ! <- nodes

      deallocate(tmp_ptags)                                ! make glob_ptags array exact size
      allocate(tmp_ptags( glob_n_ptags ))
      tmp_ptags=glob_ptags(1:glob_n_ptags)
      deallocate(glob_ptags)
      allocate(glob_ptags( glob_n_ptags ))
      glob_ptags=tmp_ptags
      deallocate(tmp_ptags)

      end subroutine particle_list_unique_particles  !]

! ----------------------------------------------------------------------
      subroutine particle_create_joined_file  ![
      ! create joined particle file
      ! input files the particle outputs not in order,
      ! output file dimensions therefore different,
      ! hence hard-coding of output file dimensions here
      implicit none

      call create_joined_empty_file
      if(mayday) return

      call copy_global_attributes
      if(mayday) return

      ierr=nf90_def_dim (nctarg, 'six',  6,           dimid(1))      ! define dimensions:
      ierr=nf90_def_dim (nctarg, 'tag', glob_n_ptags, dimid(2))
      ierr=nf90_def_dim (nctarg, 'time', tsize,       dimid(3))      ! normally time is unlimited dimension (size=0 here),
                                                                     ! but we will write out of sequence so need known time dimension size
                                                                     ! dimsize(unlimdimid,0) = time recs
      ! Define variables & atts:
      ierr=nf90_def_var (nctarg, 'ocean_time', NF90_DOUBLE, dimid(3),   varid(1),
     &                   deflate_level=deflate_level, shuffle=shuffle)

      ierr=nf90_def_var (nctarg, 'ptag',       NF90_DOUBLE, dimid(2),   varid(3),
     &                   deflate_level=deflate_level, shuffle=shuffle)

      ierr=nf90_def_var (nctarg, 'particles',  NF90_DOUBLE, dimid(1:3), varid(2),
     &                   deflate_level=deflate_level, shuffle=shuffle)

      ! still need to copy attributes?

      ierr=nf_enddef (nctarg)                                        ! leave definition mode

# ifdef VERBOSE
      if(my_rank==0) write(*,'(/1x,A)') 'Leaving definition mode.'
# endif

#ifdef TIMING
      nclk=3-nclk
      call system_clock (iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      net_fcrt_clk=net_fcrt_clk+inc_clk                              ! timing for file creation
#endif

      end subroutine particle_create_joined_file  !]

! ----------------------------------------------------------------------
      subroutine particle_write  ![
      ! write all particles from sub-domain into global output file
      !
      ![per sub-domain, populate p_out array which is sorted by all unique particles in sub-domain.
      ! next, write one particle by writing continuous time-records where possible (as particle
      ! can move in and out of domains).
      ! once all tags complete, move on to next sub-domain. !]
      implicit none

      ! local
      integer :: indx, low_indx
                                                           ! write all  particle tags in one go:
      ierr=nf90_put_var(nctarg, varid(3), glob_ptags )     ! hard-coded varid in index (3)
      if(ierr/=nf90_noerr) print *, 'BAD WRITE - ptags', nf90_strerror(ierr),'my_rank', my_rank

      i=2                                                  ! hard-code particle var id
      do node=0,nnodes-1                                   ! loop over all input sub-domains

# ifdef TIMING
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
# endif

        if(my_rank==0) write(*,'(16x,A,I4)')
     &    'Assembly of particles in node =', node

        npart = dimsize(2,node)                            ! buffer size of particle array
        if (allocated(p_out)) deallocate(p_out)
        allocate( p_out(6,loc_n_ptags(node),tsize) )       ! don't need tag each time, hence 1-6 not 1-7 dim size
        p_out = no_rec_val                                 ! catch value for non-continuous times (could improve this)
        if (allocated(p_in))  deallocate(p_in)             ! array for read in particles. might not be big enough
        allocate( p_in(npart,7,tsize) )                    ! for next node, reallocate rather

        ! POPULATE OUTPUT ARRAY for all read in particles:
        ! might want to do for one tag at a time so that not continually looking up indx in array?
        ierr=nf90_get_var(ncid(node), vid(i,node), p_in )
        if(ierr/=0) print *, 'BAD READ - particle_write', nf90_strerror(ierr), 'node', node

        do rec=1,tsize
          do p=1,npart
            if (p_in(p,1,rec) < 1.0e30) then               ! valid tag (indx 1) < less than fill value +E36
                                                           ! need to find local index to populate local array p_out array
              ! here p_in tags not sorted so must serch from low_bound=1
              ! motivation to sort p_in online in roms rather!
              call lookup_index( loc_ptags(node)%arr, p_in(p,1,rec), loc_n_ptags(node), indx )
              p_out(:,indx,rec) = p_in(p,2:7,rec)          ! store new particle tag in list of known particles
            else
              exit                                         ! extent of particles for this record
            endif
          enddo    ! nparts
        enddo      ! recs
!        print *, 'p_out=',p_out

# ifdef TIMING
        net_read_size=net_read_size+size
        nclk=3-nclk
        call system_clock (iclk(nclk),clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_read_clk=net_read_clk+inc_clk
# endif

        start_out=1                                        ! set default of 1
        count_out=1                                        ! set default of 1
        count_out(1)=6                                     ! always 1:6 in size

        ! WRITE ONE PARTICLE AT A TIME:
        low_indx=1                                         ! set low bound of index lookup
        do p=1,loc_n_ptags(node)                           ! find particles index in global array
          call lookup_index( glob_ptags, loc_ptags(node)%arr(p), glob_n_ptags, indx, low_indx )
          start_out(2)=indx

          do rec=1,tsize                                   ! write continuous timsteps (possible not all timesteps not contained in sub-domain)
            if (p_out(1,p,rec) /= no_rec_val) then         ! no record value set everywhere in array originally
                if (cnt_rec==0) rec_str=rec
                cnt_rec=cnt_rec+1
            endif

            ! write previous records when found a no_rec_val, or if reached end of array
            if (cnt_rec>0 .and. (p_out(1,p,rec) == no_rec_val .or. rec==tsize) ) then
              start_out(3) = rec_str
              count_out(3) = cnt_rec                       ! continous timestep count
              ierr=nf90_put_var(nctarg, varid(i), p_out(:,p,rec_str:rec_str+cnt_rec-1),  ! consider array order to prevent local copy?
     &                          start_out(1:3), count_out(1:3))
              if(ierr/=nf90_noerr) then
                print *, 'BAD WRITE!', nf90_strerror(ierr),'my_rank', my_rank
                stop 'BAD WRITE'
              endif
              cnt_rec=0                    ! reset count
            endif
          enddo    ! <-- rec

        enddo      ! <-- p

#ifdef TIMING
        net_wrt_size=net_wrt_size+size
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_wrt_clk=net_wrt_clk+inc_clk
#endif

      enddo  !<-- nodes

      end subroutine particle_write  !]

! ----------------------------------------------------------------------
      subroutine insertion_sort_1p( A, n )  ![
      ! improve this to a divisional search later
      ! sort 1 value (at end of array) into correct position to already sorted list
      ! can't check for duplicate too, as algo shifts values in array before it knows the outcome.
      ! checks from top so lower bound not important
      implicit none

      ! inputs
      real,dimension(:) :: A
      integer           :: n
!      integer,optional  :: low

      ! local
      integer :: l
      real    :: val

      if (n>1) then                      ! only sort if array larger than 1

        val = A(n)                       ! store value to sort
        l = n-1                          ! start comparing 1 left of value

        do while ( val<A(l) )
          A(l+1)=A(l)                    ! store left value in position of A(l)
          l=l-1
          if (l==0) exit                 ! reached end of array (couldn't have in do while since A(l=0) illegal
        enddo

        A(l+1)=val                       ! set value in correct position (+1 as already incremented down)
      endif

      ! for when glob_ptags is same as tmp_g_ptags
!      if (present(low)) low=l+2          ! since located index is l+1, low is l+1+1 (next index up)

      end subroutine insertion_sort_1p  !]

! ---------------------------------------------------------------------
      subroutine lookup_index( A, val, asize, indx, low_indx )  ![
      ! find corresponding index in ordered list
      ! bisectional search, looks in middle of known range each iteration.
      implicit none

      ! inputs
      real, dimension(:), intent(in)    :: A
      real              , intent(in)    :: val
      integer           , intent(in)    :: asize
      integer           , intent(out)   :: indx
      integer,optional  , intent(inout) :: low_indx

      ! local
      integer :: low, imid, high
      logical :: found_ip

      if(present(low_indx)) then                           ! looking up from order list to ordered list then
          low=low_indx                                     ! previously found index is clearly the new lower bound
      else
          low=1
      endif

      found_ip=.false.
      low=1
      high=asize

      if(asize==1) then                                    ! catch if only 1 value
        if (val==A(asize)) then
          found_ip=.true.
          indx=asize
        endif
      else
        do while (.not. found_ip)
          imid=(low+high)/2                                  ! integer division - will round down
          if     (val> A(imid)) then
            low  = imid+1                                    ! index must be larger
            cycle
          elseif (val< A(imid)) then
            high = imid-1
            cycle
          elseif (val==A(imid)) then
            found_ip=.true.
            indx=imid
            cycle                                            ! cycle will exit since found_ip=.true.
          elseif (low > high)   then                         ! this stop condition is possibly wrong
            stop 'lookup_index val not found as low > high'
          else
            stop 'lookup_index algo broken'
          endif
        enddo
      endif

      if(present(low_indx)) low_indx=imid+1

      end subroutine lookup_index  !]

! ---------------------------------------------------------------------
      subroutine lookup_duplicate( A, val, asize, duplicate, low_in )  ![
      ! find if value exists already in ordered list with bisection
      ! could combine this with lookup_index into one routine
      implicit none

      ! inputs
      real, dimension(:), intent(in)    :: A
      real              , intent(in)    :: val
      integer           , intent(in)    :: asize
      logical           , intent(out)   :: duplicate
      integer, optional,  intent(inout) :: low_in

      ! local
      integer :: low, imid, high

      ! optional arguement additional if statements.
      ! Maybe should just have seperate routine for the optional arguement to avoid 'if's'?

      if(present(low_in))then
          low=low_in
      else
          low=1
      endif

      duplicate=.false.
      high=asize

      if(asize==1) then                                    ! catch if only 1 value
        if (val==A(asize)) then
          duplicate=.true.
        endif
      else
        do while (.not. duplicate)
          imid=(low+high)/2                                ! integer division - will round down
          if     (low> high   ) then
            exit                                           ! not found - exit
          elseif (val> A(imid)) then
            low  = imid+1                                  ! index must be larger
            cycle
          elseif (val< A(imid)) then
            high = imid-1
            cycle
          elseif (val==A(imid)) then
            duplicate=.true.
            exit                                           ! cycle will exit since found_ip=.true.
          else
            stop 'lookup_index algo broken'
          endif

        enddo
      endif

      if(present(low_in) .and. duplicate) low_in=imid+1

      end subroutine lookup_duplicate  !]

! ---------------------------------------------------------------------
      recursive subroutine check_buff_space_and_expand( req_size, buff_size, A ) ![
      ! check array has enough buffer space to receive more particles
      ! expand array if not.
      ! this algo needs more consideration, could get major over-allocation
      ! for very irregular particle distributions.
      implicit none

      ! inputs
      integer,                      intent(in)    :: req_size
      integer,                      intent(inout) :: buff_size
      real,dimension(:),allocatable,intent(inout) :: A

      ! local
      real :: ratio_rem                                    ! ratio of remaining nodes to estimate remaining particles
      real,dimension(:),allocatable :: A_tmp
      integer :: buff_size_old

      if (req_size > buff_size) then

        print *, 'reallocting larger global array size'
        print *, '  old buff_size=', buff_size

        buff_size_old = buff_size
        ratio_rem = 1 + real(nnodes-node+1) / real(node+1) ! +1 since node starts from 0. node is denominator as required increase.
        buff_size = buff_size * ratio_rem
        if (node==nnodes-1) buff_size = buff_size * 1.01   ! as ratio_rem = 1 for last node

        allocate(A_tmp(buff_size_old))
        A_tmp = A
        deallocate(A)
        allocate(A(buff_size))
        A(1:buff_size_old) = A_tmp
        deallocate(A_tmp)

        if (req_size > buff_size) then                      ! incase still not big enough
          call check_buff_space_and_expand( req_size, buff_size, A )
        endif

        print *, '  new buff_size=', buff_size

      endif

      end subroutine check_buff_space_and_expand  !]
! ---------------------------------------------------------------------

      end module particle_join_mod

! ---------------------------------------------------------------------

      program particle_join

      use ncjoin_mod
      use particle_join_mod

      implicit none

      deflate_level=0                  ! leave as is - no compression and no shuffle.
      shuffle=.false.                  ! just needed for compatibility with shared ncjoin_mod routines.

      call init_timing_and_vars
      arg=0                            ! set to 0 for serial version as arguments of program start at 1, unlike ncjoin_mpi.

      my_rank=0                        ! allows easy use of ncjoin_mpi routines, even though serial program

      call particle_check_args         ! NEW

      do while (arg .lt. nargs)

      nnodes=-1                        ! DevinD repeated to get 'goto' outside of master only region, as called by all procs. ! used to be 11 marker here
      mayday=.false.                   ! reset mayday flag

      if(my_rank==0) then              ! leave here incase make mpi version MPI MASTER ONLY: PRE-PROCESSING. i.e. check input files, creating output file, etc.
                                       ! Extract a set of files which cover the whole physical grid.

        write(*,'(/1x,A/)') 'Pre-processing input files...'

        ! NEW: believe this is legal. Partition still exists in particle file. Some redundant allocation perhaps.
        call check_partial_file_set    ! PART 1: CHECK ALL INPUT PARTIAL FILES ARE CORRECT
        if (mayday) goto 23            ! Using goto the break from if(my_rank==0)
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

            call check_ndims_nvars_natts(node)                       ! NEW: not sure if should bother using this. Have so few dims in file...
            if (mayday) goto 23

            call create_catalog_of_var_names_IDs_ranks(node)
            if (mayday) goto 23

                                                                     ! NEW: don't need to do this anymore, might want to do it though for efficiency later.
                                                                     !      might not like loads of files open...
!            if (node.gt.0) then                                     ! close all the files, except for node=0.
!              ierr=nf_close(ncid(node))                             ! since master only, need to close files to open collectively later.
!              ncid(node)=-1                                         ! keep node=0 open as still using below.
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

      !] END PART 2 (verify all dimensions and variables)

        call particle_list_unique_particles                          ! PART 3:

        ! NEW: need to move to after parsing of all files for number of particles
        call particle_create_joined_file                             ! PART 4:

        if (mayday) goto 23

  23    if(my_rank==0) write(*,'(/1x,A/)')
     &                          'End of master proc pre-processing.'

      endif ! <- if(my_rank==0) END OF MASTER ONLY PRE-PROCESSING

      if (mayday) goto 97

        if(my_rank==0) write(*,'(16x,A)')                             ! WRITE OCEAN_TIME:
     &    'Assembly of ocean_time'
        node=0; i=1                                                   ! set variable to ocean_time (should be 1)
        vartype(i)=nf90_double                                        ! hack because vartype wasn't done in pre-proc for some reason
        do rec=1, tsize
          call read_write_non_partitioned_var                         ! write ocean_time from 1 node only
        enddo

        call particle_write                                           ! PART 5:

        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false.  ! most errors arrive here

  98    call close_file_set

      enddo !<- do while (arg .lt. nargs)

      call display_timing_summary


      end program particle_join
