      ![ INFO:
      !
      !  Tool to join output files containing output data from extract_data.F.
      !  Typically, this data is for child boundary conditions or mooring data.
      !
      !  Chose not to add to ncjoin_mpi because it differs too much and would pollute
      !  ncjoin_mpi with too many catches for joining child bry. Can still use several
      !  of the same routines from ncjoin_mod module.
      !
      !  Long term fixes:
      !  - currently only handles s_rho sized vertical dimension.
      !  - really online extract_data shouldn't write empty files, and therefore this shouldn't
      !    need to deal with empty files.
      !    Online files should contain list of node numbers to join, and use that rather
      !    than pre-process.
      !  - Could parallelize in future.
      !
      !  Methodology:
      !
      !  - inputs: variables share dimensions where possible. Dimensions are simply numbered sequentially.
      !            the actual dimension name is stored as an attribute and used for the joined file.
      !            dims:
      !            time_*(digits) - needs time since different sets can have different output rates
      !            np_*(digits)   - number of points in sub-domain (goes with matching time dimension)
      !            vars:
      !            *(object)_time - one time variable for object group that share dimension and time
      !            atts:
      !            dname          - actual dimension name for the joined file
      !            dsize          - global joined dimension size
      !
      !  - join output:
      !          - assume best to do all vars and records in file and then move on, rather than
      !            continually looping through files (appeared to slow in earlier versions of this tool)
      !          - efficient to store variable ID so don't have to look it up for each file?
      !            not sure because unlike ncjoin this typically only has a few relevant sub-domains.
      !]

      ! PROGRAM SECTIONS:
!      ! - PART 1: Confirm all input files are correct and related to each other
!      ! - PART 2: Veryify all variables and dimensions in input files
!      ! - PART 3: Create output file
!      ! - PART X: Collectively (MPI-I/O) open all input and output files
!      ! - PART X: Read/write partitioned variables   (bulk of work done here)
!      ! - PART X: Close all files
!      ! - PART X: Display program timing summary


![ CPP DEFS:
! Delete partial files after joining:
c#define DEL_PART_FILES
! Document program execution time:
#define TIMING
! Verbose terminal output:
c--#define VERBOSE
!]

      module extract_data_join_mod

      use ncjoin_mod

      implicit none

      ! ------ CHILD_BRY VARS -------------
      !  Added xdata_ to beginning of routines that are specifically for child boundary,
      !  not ncjoin_mpi.
      integer :: global_south_bry_size    ! size of join dimension in child boundary file
      integer :: global_west_bry_size
      integer :: tot_pwithc               ! total parent domains with child bry
      integer :: max_cvars_p              ! max number of child variables contained in sub-domain
      integer :: nd, vr                   ! loop index
      logical :: new_var                  ! trigger reading start and count for new variable
      integer, dimension(4) :: vid_out    ! variable ID's in joined output file NOT SURE I NEED THIS???
      integer, dimension(3) :: child_start, child_count, startin, countin  ! start+counts for read/write
      integer, dimension(:),   allocatable :: tot_cvars_p_all, tot_cvars_p ! child vars in parent (all nodes and actual)
      integer, dimension(:),   allocatable :: pwithc_nodes                 ! list of whether parent contains child
      integer, dimension(:,:), allocatable :: cvars_p_all, cvars_p         ! each child var number in parent
      real,    dimension(:,:), allocatable :: bry_data

      ! XDATA VARS:
      integer,parameter :: max_sets=20                               ! max number of expected sets - need check for not exceeded!!!
      integer,parameter :: max_evars=36                              ! max vars per set - z,S,T,ub,u,up,vb,v,vp - need check (9x4=36)
      integer :: nsets=0                                             ! count of sets (root name)
      logical :: found
      integer :: iv, g, id
      integer,dimension(4) :: vdimids, vdimids_out
      integer,dimension(:),allocatable :: sub_cnt

      type extract_set                                               ! all variables associated with e.g. child grid, mooring, etc
        integer,pointer :: nvars                                     ! number of vars in set (actual want target but can't do it in derived type)
        integer :: s_rho=1                                           ! vertical dimension size
        logical :: time_done=.false.
        integer :: tsize                                             ! time size
        integer :: nstr                                              ! index after set name in full name
        character(len=30) :: tname = ' '
        character(len=30) :: sname = ' '                             ! set name (root name)    - set as blank in order
        character(len=:),dimension(:),allocatable :: vname           ! variable  name (nvars)  - to compare to read in name
        character(len=:),dimension(:),allocatable :: dname           ! dimension name
        character(len=:),dimension(:),allocatable :: lname           ! long name
        character(len=:),dimension(:),allocatable :: uname           ! long name
        integer,dimension(:),allocatable      :: vndims              ! dimension of object (nvars)
        integer,dimension(:),allocatable      :: dsize               ! dimension of object (nvars)
        integer,dimension(:),allocatable      :: nmax                ! number of vertical levels
        integer,dimension(:),allocatable      :: nsubs               ! number of sub-domains with var
        integer,dimension(:,:),allocatable    :: lsubs               ! list of sub-domains with var (nnodes,nvars)
        integer,dimension(:),allocatable      :: vid_out             ! varid in output file - NOT sure if needed as do at once...
!        integer,dimension(:),allocatable      :: vid_out
        ! don't need start and count as can just get it when I need it since writing one variable
        ! and one sub-domain at a time
      end type extract_set

      type(extract_set),dimension(max_sets) :: s

      contains

!----------------------------------------------------------------------
      subroutine xdata_check_args  ![
      ! ensure reasonable number of arguements used,
      ! if not give usage instructions.
      implicit none

      if(nargs < 4) then                                                       ! too few program arguments
        write(*,'(/1x,A//10x,A)') 'Correct usage of extract_data_join:',
     &    'extract_data_join      sample_extract.00000.*.nc'
        write(*,'(1x,A/10x,A/)') 'or',
     &    'extract_data_join  -d  sample_extract.00000.*.nc'
        write (*,'(/1x,A//1x,A/1x,A/)')
     &    'where -d deletes the partial files.',
     &    'This tool joins partial files of online interpolated',
     &    'variables created by a ROMS simulation.'

        error stop
      endif

      end subroutine xdata_check_args  !]

! ----------------------------------------------------------------------
      subroutine xdata_create_joined_file  ![
      implicit none  ! create joined child bry file

      call xdata_create_joined_empty_file
      if(mayday) return

      call copy_global_attributes
      if(mayday) return

      ierr=nf90_def_dim(nctarg, 'time', s(g)%dsize(1), vdimids_out(3) )
      ierr=nf90_def_dim(nctarg, 's_rho', s(g)%s_rho, vdimids_out(2) )


      do iv=1,s(g)%nvars
        ierr=nf90_def_dim(nctarg, s(g)%dname(iv), s(g)%dsize(iv), vdimids(1))  ! define dimensions: (repeated but doesn't matter)
      enddo

      ierr=nf90_def_var(nctarg, 'ocean_time',nf90_double,
     &                      vdimids_out(3), s(g)%vid_out(1))
      do iv=2,s(g)%nvars
        if (s(g)%vndims(iv)==1) then
!         ! ignore redundant time dimensions
        else
          ierr=nf90_inq_dimid(nctarg, s(g)%dname(iv), vdimids_out(1))          ! only dimension that can varys
          lncn=lenstr(s(g)%vname(iv))
          if (s(g)%vndims(iv)==2) then
            ierr=nf90_def_var(nctarg, s(g)%vname(iv)( s(g)%nstr:lncn ),nf90_double,
     &                        (/ vdimids_out(1), vdimids_out(3) /), s(g)%vid_out(iv))
          else
            ierr=nf90_def_var(nctarg, s(g)%vname(iv)( s(g)%nstr:lncn ),nf90_double,              ! add compression here if desired
     &                      vdimids_out(1:3), s(g)%vid_out(iv))
          endif
        endif
        ierr=nf90_put_att(nctarg,s(g)%vid_out(iv),'long_name',s(g)%lname(iv))
        ierr=nf90_put_att(nctarg,s(g)%vid_out(iv),'units',    s(g)%uname(iv))
      enddo

      ! Leave definition mode

!      ierr=nf90_close(ncid(0)) ! close last input file. completed master only reads
!      ncid(0)=-1

      ierr=nf90_enddef (nctarg)
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

      end subroutine xdata_create_joined_file  !]

! ----------------------------------------------------------------------
      subroutine xdata_read_write  ![
      ! due to xdata_find_relevant_parent routine, we know this
      ! variable is contained in parent and msut be written
      implicit none

      integer,save :: vid_in

# ifdef TIMING
      nclk=3-nclk
      call system_clock (iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      net_gray_clk=net_gray_clk+inc_clk
# endif


      if (new_var) then                                              ! once per variable (don't repeat per record)
        startin=1; countin=1; child_start=1; child_count=1

        if (ncid(node).eq.-1)
     &    ierr=nf90_open (ncname(node), nf90_nowrite, ncid(node))

        if (s(g)%vndims(iv) >= 2) then  ! exclude time (1D vars)
!          print *, 'xdata_read_write read_write - s(g)%vname(iv)=',s(g)%vname(iv)
          ierr = nf90_inq_varid(ncid(node), s(g)%vname(iv), vid_in)
          if (ierr/=0) call handle_ierr(ierr,'xdata_read_write nf90_inq_varid var=',s(g)%vname(iv))
          ierr = nf90_get_att(ncid(node), vid_in, 'start', child_start(1))  ! start position in global child bry
          ierr = nf90_get_att(ncid(node), vid_in, 'count', child_count(1))  ! count of sub-domain entry in global child bry                            ! first var is time in input files
        else
!          print *, 'xdata_read_write read_write - s(g)%tname=',s(g)%tname
          ierr = nf90_inq_varid(ncid(node), s(g)%tname, vid_in)
          if (ierr/=0) call handle_ierr(ierr,'xdata_read_write nf90_inq_varid var=',s(g)%tname)
        endif
        countin(1:2)=(/ child_count(1), s(g)%nmax(iv) /)

        if (allocated(bry_data)) deallocate (bry_data)                           ! variable for reading data
        allocate ( bry_data( child_count(1), s(g)%nmax(iv) ) )
        new_var = .false.
!        print *, 'DONE ALLLLLOC'
!        print *, 'child_count(1)',child_count(1)
!        print *, 's(g)%nmax(iv)',s(g)%nmax(iv)
      endif  !<-- new_var

      if     (s(g)%vndims(iv) == 1) then                   ! set new record
        startin(1)     = rec
        child_start(1) = rec
        countin(1)=1
      elseif (s(g)%vndims(iv) == 2) then                   ! handle 2D or 3D variable dims
        startin(2)     = rec                               ! set recorded timestep for read in partial bry
        child_start(2) = rec                               ! set recorded timestep for written join bry
      elseif (s(g)%vndims(iv) == 3) then
        startin(3)     = rec                               ! set recorded timestep for read in partial bry
        child_start(3) = rec                               ! set recorded timestep for written join bry
        child_count(2) = s(g)%s_rho                        ! should be nmax rather?
      else
        error stop 'ERROR: dimensions issue when writing!'
      endif

      ierr = nf90_get_var( ncid(node), vid_in, bry_data,   ! read variable
     &                     start=startin(1:s(g)%vndims(iv)), count=countin(1:s(g)%vndims(iv)) )
      if (ierr/=0) then
        print *, 'node=', node
        call handle_ierr(ierr,'xdata_read_write read var=',s(g)%vname(iv))
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

      ierr = nf90_put_var( nctarg, iv, bry_data,           ! write variable
     &                     start=child_start(1:s(g)%vndims(iv)), count=child_count(1:s(g)%vndims(iv)) )
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

      end subroutine xdata_read_write  !]

! ----------------------------------------------------------------------
      subroutine xdata_copy_vars_atts  ![
      implicit none                   !  copy the variables and attributes for extract_data_join.
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

      end subroutine xdata_copy_vars_atts  !]

! ----------------------------------------------------------------------
      logical function findstr(all_string,want_string,istart) ![
      ! this utility function should be moved somewhere more generic...
      ! copied from extract_data.F.
      ! should be in roms_read_write and then that module should be available to
      ! Tools-Roms/ as well as src/
      ! didn't do it now because roms_read_write depends on other modules so would
      ! require a bit of work (do once ncvars is removed)
      implicit none

      !input/output
      character(len=*),intent(in)  :: all_string           ! entire string
      character(len=*),intent(in)  :: want_string          ! desired string to find within all_string
      integer,optional,intent(out) :: istart

      !local
      integer :: nl,nlv,i

      nl  = len(trim(all_string))
      nlv = len(want_string)

      findstr = .false.
      do i = 1,nl-nlv+1
         if (all_string(i:i+nlv-1) == want_string) then
          findstr = .true.
          if (present(istart)) then
            istart=i                                       ! return string starting index
          endif
          exit
         endif
      enddo

      end function findstr !]

!----------------------------------------------------------------------
      subroutine xdata_preprocess_lists  ![
      ! create list of sub-domains contains particles
      ! input dimensions are meaningless
      ! each input variable contains the output variable dimension
      ! current method relies on sets being grouped together in input files
      implicit none

      ! local
      integer :: fnvars, vndims, istart, varid, sdimid
      integer,pointer :: s_nv
      character(len=60) :: name

      print *, 'Preprocessing node = ',node

      ierr=nf90_inquire(ncid(node), nVariables=fnvars)
      if (ierr/=0) call handle_ierr(ierr,'xdata_preprocess_lists nf90_inquire')

      if (fnvars>0) then                                   ! only preprocess if sub-domain contains data
        do varid=1,fnvars                                  ! varids are simply 1:fnvars
          ierr=nf90_inquire_variable(ncid(node), varid, name, ndims=vndims, dimids=vdimids)
          if (ierr/=0) call handle_ierr(ierr,'xdata_preprocess_lists nf90_inquire_var')

          found=findstr(name,'_',istart)                   ! find set (root) name before first underscore (istart = _)
          found=.false.
          if (nsets>0) then                                ! see if set found before
            do g=1,nsets
              if (trim( s(g)%sname )==name(1:istart-1)) then
                found=.true.
                exit
              endif
            enddo
          endif
          if (.not. found) then                            ! new set, allocate new set variables
            nsets=nsets+1
            s(nsets)%sname=name(1:istart-1)
            print *, 'found set =',s(nsets)%sname
            print *, 'nsets =',nsets
            s(nsets)%nstr=istart+1
            allocate ( s(nsets)%nvars )
            s(nsets)%nvars=0
            allocate(character(50) :: s(nsets)%vname( max_evars) )
            allocate(character(50) :: s(nsets)%dname( max_evars) )
            allocate(character(90) :: s(nsets)%lname( max_evars) )
            allocate(character(30) :: s(nsets)%uname( max_evars) )
            allocate( s(nsets)%vndims( max_evars ))
            allocate( s(nsets)%dsize(  max_evars ))
            allocate( s(nsets)%nmax(   max_evars ))
            allocate( s(nsets)%nsubs(  max_evars ))
            s(nsets)%nsubs=0
            allocate( s(nsets)%lsubs( nnodes, max_evars )) ! space for all sub-domains (could reduce this?)
            allocate( s(nsets)%vid_out( max_evars ))

            ierr=nf90_inq_dimid(ncid(node), 's_rho', sdimid)
            ierr=nf90_inquire_dimension(ncid(node), sdimid, len=s(nsets)%s_rho)  ! store s_rho dimension

            g=nsets
          endif

          s_nv => s(g)%nvars
          found=.false.                                    ! search set variables for this variable
          if (s_nv>0 .and. vndims>1) then
            do iv=1,s_nv
              if (trim( s(g)%vname(iv) ) == trim(name)) then
                found=.true.
!                print *, found, 'iv=',iv
                exit
              endif
            enddo
          endif
          if (s_nv>0 .and. vndims==1) found=.true.         ! ignore additional time variables as redundant
          if (.not. found) then                            ! new variable, store details
            s_nv=s_nv+1
            iv=s_nv                                        ! need iv later as used if (found) above
            s(g)%vname(s_nv)=name
            s(g)%vndims(s_nv)=vndims
            print *, 'found var =',s(g)%vname(s_nv)
            print *, 'nvars =',s_nv

            ! since will need to find variable ID in file anyway, it I can also get num dims and dim name
            ! so might not need to store this but doing it for now.
            ! definitely don't need start and count size doing entire var at once for sub-domain
            ! so may as well get it once then and not store for all vars..
            if (vndims==1) then                            ! must be time if 1 dim only
              s(g)%tname=name
              print *, 's(g)%tname',s(g)%tname
              s(g)%vname(s_nv)='ocean_time'
              s(g)%dname(s_nv)='time'
              ierr=nf90_inquire_dimension(ncid(node), vdimids(1), len=s(g)%dsize(s_nv))
              s(g)%tsize = s(g)%dsize(s_nv)                ! save the size fo timedimension
            else
              ierr=nf90_get_att(ncid(node), varid, 'dname', s(g)%dname(s_nv))
              if (ierr/=0) call handle_ierr(ierr,'xdata_preprocess_lists attribute dname missing')
              ierr=nf90_get_att(ncid(node), varid, 'dsize', s(g)%dsize(s_nv))
              if (ierr/=0) call handle_ierr(ierr,'xdata_preprocess_lists attribute dsize missing')
            endif
!            print *, 'joined dname =',s(g)%dname(s_nv)  ! store name only, don't have to group can just inquire dimid from dim name...
!            print *, 'joined dsize =',s(g)%dsize(s_nv)  ! store name only, don't have to group can just inquire dimid from dim name...
!            print *, ' '

            if (vndims==3) then
              s(g)%nmax(s_nv)=s(g)%s_rho
            else
              s(g)%nmax(s_nv)=1
            endif

            ierr=nf90_get_att(ncid(node), varid, 'long_name', s(g)%lname(s_nv))
            if (ierr/=0) call handle_ierr(ierr,'xdata_preprocess_lists attribute longname missing')
            ierr=nf90_get_att(ncid(node), varid, 'units', s(g)%uname(s_nv))
            if (ierr/=0) call handle_ierr(ierr,'xdata_preprocess_lists attribute units missing')
          endif

          if (s_nv==1 .or. vndims>1) then                ! ignore additional redundant time variables
            s(g)%nsubs(iv)=s(g)%nsubs(iv)+1                  ! count of sub-domains that contain var
            s(g)%lsubs( s(g)%nsubs(iv), iv ) = node        ! store sub-domain number for writing loop
          endif

        enddo      !<- varid

      endif

      if (node==nnodes-1) then
        do g=1,nsets
          print *, ' '
          print *, 'sname=',s(g)%sname
          do iv=1,s(g)%nvars
            print *, 'vname=',s(g)%vname(iv)
            print *, 'nsubs=',s(g)%nsubs(iv)
            print *, 'lsubs=',s(g)%lsubs(1:s(g)%nsubs(iv),iv)
            print *, ' '
          enddo
        enddo
      endif

      end subroutine xdata_preprocess_lists  !]

! ----------------------------------------------------------------------
      subroutine xdata_create_joined_empty_file  ![
      implicit none

      integer :: idot

        i=lenstr(root_bak)
        found=findstr(root_bak,'.',idot)
        if (idot==1) stop 'cannot have . at start of input file names'
        j=lenstr(sffx_bak)
        if (root_bak(i:i).eq.'.' .and. sffx_bak(1:1).eq.'.') then
          nctargname=trim( s(g)%sname )/ /root_bak(idot:i)/ /sffx_bak(2:j)
        else
          nctargname=trim( s(g)%sname )/ /root_bak(idot:i)/ /sffx_bak(1:j)
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

!        ierr=nf90_create (nctargname(1:ltrg), IOR(NF90_NETCDF4, NF90_MPIIO),       ! for not master only version that worked with independent mode!
!     &                    nctarg, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
        ierr=nf90_create (nctargname(1:ltrg), NF90_NETCDF4, nctarg)                 ! Master only works!
!        ierr=nf90_create (nctargname(1:ltrg), IOR(NF90_NETCDF4, NF90_MPIIO), nctarg) ! Master only didn't work (hdf error at write)
        if (ierr .eq. nf90_noerr) then
          if(my_rank==0) write(*,'(/1x,3A)')  'Created netCDF file ''',
     &                        nctargname(1:ltrg), '''.'
          if(my_rank==0) write(*,'(1x,A,I1)')
     &      '-> compression level (deflate_level) = ', deflate_level
          if(my_rank==0) write(*,*)
     &      '-> compression shuffle = ', shuffle
        else
          write(*,'(/1x,4A/12x,A/)')     '### ERROR: Cannot create ',
     &                          'netCDF file ''', nctargname(1:ltrg),
     &                                      '''.', nf90_strerror(ierr)
          mayday=.true.; return ! goto 97
        endif

      end subroutine xdata_create_joined_empty_file  !]

!----------------------------------------------------------------------

      end module extract_data_join_mod

!----------------------------------------------------------------------

      program extract_data_join  ![
      use ncjoin_mod
      use extract_data_join_mod

      implicit none


      deflate_level=0                                      ! leave as is - no compression and no shuffle.
      shuffle=.false.                                      ! just needed for compatibility with shared ncjoin_mod routines

      call init_timing_and_vars
      arg=0                                                ! set to 0 for serial version as arguments of program start at 1, unlike ncjoin_mpi.
      my_rank=0                                            ! allows easy use of ncjoin_mpi routines

      call xdata_check_args

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
            lncn=lenstr(ncname(node))  !? might be able to delete

            if (ncid(node).eq.-1) ierr=nf90_open(ncname(node),nf90_nowrite, ncid(node))
            if (ierr .eq. nf_noerr) then

              call xdata_preprocess_lists

              if (node.gt.0) then                            ! close all the files, except for node=0.
                ierr=nf90_close(ncid(node))                  ! since master only, need to close files to open collectively later.
                ncid(node)=-1                                ! keep node=0 open as still using below.
              endif

            else
              write(*,'(/1x,A,1x,3A/14x,A)')    '### ERROR: Cannot ',
     &                 'open netCDF file ''', ncname(node)(1:lncn),
     &                                    '''.', nf_strerror(ierr)
              goto 97
            endif
          enddo  !<-- node=0,nnodes-1

  23      if(my_rank==0) write(*,'(/1x,A/)')
     &                          'End of master proc pre-processing.'

      !] END PART 2: VERIFY ALL VARIABLES & DIMENSIONS

        endif ! <- if(my_rank==0) END OF MASTER ONLY PRE-PROCESSING

        if (mayday==.true.) goto 97

        do g=1,nsets

          call xdata_create_joined_file                                ! PART 3:

          if (allocated(sub_cnt)) deallocate(sub_cnt)
          allocate( sub_cnt( s(g)%nvars ))
          sub_cnt=1

          ! prob don't need if statement and time_done because not looping through other times
          if (.not. s(g)%time_done) then                               ! time variable comes first from online inputs
            if(my_rank==0) write(*,'(16x,A)')
     &            'Assembly of ocean_time'
            node=s(g)%lsubs(1,1)                                       ! first node containing time
            iv=1                                                       ! time variable varid=1
            new_var=.true.
            do rec=1, s(g)%tsize
              call xdata_read_write                                    ! write ocean_time from 1 node only
            enddo
            s(g)%time_done=.true.                                      ! don't need this!
          endif

          do node=0,nnodes-1                                           ! need all nodes since don't know which are contained for set
            do iv=2,s(g)%nvars                                         ! =2 since time already done above
              if (s(g)%lsubs( sub_cnt(iv),iv )==node) then             ! var contained in sub-domain (node)

                if(my_rank==0) write(*,'(14x,A,I4,2A)')
     &            'joining node ',node,' var=',s(g)%vname(iv)

                new_var=.true.
                do rec=1, s(g)%tsize
                  call xdata_read_write
                enddo

                sub_cnt(iv)=sub_cnt(iv)+1
              endif  !<-- lsubs
            enddo    !<-- vars
          enddo      !<-- nodes
        enddo        !<-- nsets

        if (ierr.eq.nf90_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false.                                  ! most errors arrive here

  98    call close_file_set                                ! PART 10

      enddo          !<- do while (arg .lt. nargs)

      call display_timing_summary

      end program extract_data_join  !]


