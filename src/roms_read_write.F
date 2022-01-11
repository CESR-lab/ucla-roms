      module roms_read_write
      ! ROMS specific read write functionality

      use netcdf
      use ncvars
      use scalars
!     use param

      implicit none
      private

      ! Type ncbry contains all required netcdf variables and input data array for frc variable.
      ! Set the name and time_name of variable to match what is in the input forcing files.
      ! Use same name as bulk_frc variable name, e.g. uwnd has nc_uwnd for netcdf vars.
      ! type (ncvs) :: nc_uwnd  = ncvs( name='uwnd',       time_name='wnd_time'  )
      type, public :: ncbry
        character(len=20)              :: vname  ! name of variable in input file
        character(len=20)              :: tname  ! time variable name for variable
        real,dimension(:,:,:),pointer  :: vdata  ! need=0 to initialize all vars other than name and time_name. Slightly inefficient to set=0 though. stores raw input data for 2 times: time_A < mode
        character(len=max_name_size),dimension(:),pointer :: flist  ! list of file names to look into for the variable
        integer              :: coarse = 0       ! Flag to indicate interpolation of the input data
        integer              :: ifile=0          ! Index to file in list of filenames
        integer              :: irec =0          ! Record number in file
        integer :: it1 = 1, it2 = 2              ! used to cycle between correct entries of 'data' above
        real,   dimension(2)  :: times = -99 ! [-99,-99]       ! stores 2 times that go with 'data' above
      end type ncbry

      interface set_frc_data
        module procedure  set_frc_data_2D
      end interface

      public :: set_frc_data

      contains

! ----------------------------------------------------------------------
      subroutine set_frc_data_2D(nc_var,var)
      implicit none

      ! input/outputs
      type(ncbry)        ,intent(inout):: nc_var ! derived type containing all neccesary supporting data
      real,dimension(:,:),intent(out)  :: var    ! time interpolated forcing data

      ! local
      character(len=20) :: vname  ! name of variable in input file
      character(len=20) :: tname  ! time variable name for variable
      real,dimension(2) :: vtimes ! times of the var data slices
      integer           :: it1,it2
      integer           :: ifile,irec
      integer           :: coarse
      real,dimension(:,:,:),pointer :: vdata
      character(len=max_name_size),dimension(:),pointer :: flist
      integer :: tmp, i, j,ierr,ncid
      real    :: cff1,cff2  ! for time interpolations
      real    :: modtime ! model time

      vname  = nc_var%vname
      tname  = nc_var%tname
      vtimes = nc_var%times
      it1    = nc_var%it1
      it2    = nc_var%it2
      ifile  = nc_var%ifile
      irec   = nc_var%irec
      coarse = nc_var%coarse

      vdata => nc_var%vdata
      flist => nc_var%flist

      modtime = tdays + 0.5*dt*sec2day    ! input data times in days. Interpolated time is based on 1/2 step time

      if (vtimes(it2) < modtime) then  ! need to refresh data

        if (ifile == 0) then   ! the first time, we need to fill both slots

          call find_new_record(vname,tname,modtime,ifile,irec,
     &                           vtimes(it1),flist)

          !! read forcing from file

          ierr=nf90_open(flist(ifile),nf90_nowrite, ncid)
          if (coarse) then
            print *,'coarse to fine not implemented here'
            stop
!           cdata => nc_var%cdata
!           ncread(ncid,vname,cdata,(/1,1,irec/))
!           call coarse2fine(cdata,vdata(:,:,it1))
          else
            !! Allocate vdata to be the right size for the
            !! subdomain (including whether we need a buffer or not
            call ncread(ncid,vname,vdata(:,:,it1),(/1,1,irec/))
          endif
          ierr = nf90_close(ncid)

          if (mynode == 0) call display_read_time_to_log(vname,vtimes(it1),irec )

        else
          tmp = it1
          it1 = it2 ! if it1 = 1, it now equals 2, and vice-versa
          it2 = tmp
        endif

        call find_new_record(vname,tname,modtime,ifile,irec,
     &                       vtimes(it2),flist)

        ierr=nf90_open(flist(ifile),nf90_nowrite, ncid)
        call ncread(ncid,vname,vdata(:,:,it1),(/1,1,irec/))
        ierr = nf90_close(ncid)

        if (mynode == 0) call display_read_time_to_log(vname,vtimes(it2),irec )

      endif ! 

      nc_var%times = vtimes
      nc_var%it1   = it1
      nc_var%it2   = it2
      nc_var%ifile = ifile
      nc_var%irec  = irec

      ! Temporal interpolation

      cff1=( vtimes(it2)-modtime )/( vtimes(it2)-vtimes(it1) )
      cff2=( modtime-vtimes(it1) )/( vtimes(it2)-vtimes(it1) )

      if (cff1.ge.0. .and. cff2.ge.0.) then

        var = cff1*vdata(:,:,it1) + cff2*vdata(:,:,it2)

      else
        write(*,'(/1x,4A/3(1x,A,F16.10)/)') 'ERROR: set_forcing',
     &  ':: Model time outside bounds of variable (times):',
     &     vname, tname,
     &                    'start =',    vtimes(it1),
     &                    'tmid_days=', modtime,
     &                    'end =',      vtimes(it2)
        error stop 'ERROR: set_forcing :: time interpolation error'
      endif

      end subroutine set_frc_data_2D 
      
! ----------------------------------------------------------------------
      subroutine find_new_record(vname,tname,time,ifile,irec,vtime,flist) ![

      ![================================================================= 
      ! Finds the index for the filename and record number of the first record
      ! with a time that is larger than the model time. It will stores
      ! the corresponding time of the record in vtime.

      ! When called for the first time, it will find the last record with a 
      ! time smaller than the model time.
      !
      ! Find: if var contained in file
      !
      ! Check: see if model time is bounded by model
      !
      !   -> yes: save record (irec), save ifile and save input time (vtime)
      !
      !   ->  no: try another file
      !
      ! For the very first timestep we need the record just before
      ! vtime(1)<model_time. Therefore irec = irec_old.
      ! E.g. if record 3 gives  vtime(1)>time_model
      ! then record 2 will give vtime(1)<time_model
      !
      ! Label (A):
      ! Handling first data read if model time is between 2 files.
      ! The record 1 of the later file will give vtime(1)>time_model,
      ! we use the last record for the previous file for that
      ! variable. This does assume there are no missing periods between
      ! the data, and the data is chronological.
      !
      ! However, if the model time is between two input files for the
      ! very first timestep: the record 1 of the first file will give
      ! vtime(1)>time_model, but we can't simply use irec = irec_old,
      ! as record 0 doesn't exist.
      ! Thus if irec=irec_old=0, this is not possible and is caught by the code.
      !
      ! Limitations
      ! -----------
      ! 1) Files for a variable must be in chronological order in input
      !    file list of forcing files
      ! 2) Files for the same variable do not need to be grouped together
      !    but it is slightly more efficient if they are.
      ! 3) Algorithm does not detect a missing year of data.
      !    E.g. if modeling 2005-2007, and annual data for 2006 is missing,
      !    then times will just be interpolated between end of 2005 & start 2007.
      !
      !]================================================================= 

      implicit none

      ! input/output
      character(len=*)              ,intent(in)   :: vname   ! name of the variable
      character(len=*)              ,intent(in)   :: tname   ! name of the corresponding time
      real                          ,intent(in)   :: time    ! target time
      integer                       ,intent(inout):: ifile   ! index for file name in list
      integer                       ,intent(inout):: irec    ! variable record entry
      real,dimension(1)             ,intent(out)  :: vtime   ! read in variable time
      character(len=*) ,dimension(:),intent(in)   :: flist   ! list of file names to search in

      ! local
      integer, dimension(1) :: dimids     ! time dimension ID
      real    :: time_old
      logical :: found_var                ! is variable in file
      logical :: found_rec                ! is correct record found of var
      logical :: found_var_ever           ! if variable was found in any previous files
      integer :: ncid, ierr
      integer :: nfiles                   ! total number of forcing files
      integer :: irec_old, ifile_old
      integer :: max_recs                 ! total variable records in file
      integer :: varid                    ! variable ID
      integer :: timeid                   ! variable time id
      logical :: first                    ! first time call for this variable


      if (ifile==0) then  ! First time for this variable
        first = .true.
        ifile = 1
        irec_old  = 0
        ifile_old = 0
      else
        first = .false.
      endif
      
      found_rec      = .false.
      found_var_ever = .false.

      do while (.not.found_rec .and. ifile <= max_frc) ! While record not found & still more files to check

        found_var = .false.  ! reset for new file
        ierr = nf90_open(frcfile(ifile), nf90_nowrite, ncid)
        if (ierr/=0) call handle_ierr(ierr,
     &   'Find_new_record: frc_file not found:  ',frcfile(ifile))

        ierr = nf90_inq_varid(ncid, vname, varid)   ! Check if this file contains the variable

        if (ierr==0) then ! Variable found, now check if file times contain model time

          found_var      = .true.
          found_var_ever = .true.
           
          ierr = nf90_inq_varid(ncid,tname,timeid)  ! Get time ID

          ! Number of records in file
          ierr = nf90_inquire_variable(ncid, timeid, dimids = dimids)
          ierr = nf90_inquire_dimension(ncid, dimids(1), len=max_recs)

          ! check if the time has a cycle length
          ierr=nf90_get_att(ncid, timeid, 'cycle_length', cycle_length) 
          if (ierr==nf90_noerr) then
            call find_new_record_cycle(ncid,timeid,max_recs,cycle_length)
            found_rec = .true.
          endif

          do while (irec<max_recs .and. .not. found_rec)                ! Search through the records until correct time is found

            irec=irec+1
            ierr=nf90_get_var(ncid,timeid,vtime, (/irec/), (/1/))       ! Get vtime (nf90 needs array hence syntax)


            if (vtime(1)>time) then   ! Correct time
              found_rec = .true.
            else                      ! Not yet
              irec_old = irec
              ifile_old = ifile     ! Only needed initially. See label (A)
              time_old = vtime(1)
            endif

          enddo  ! while not found and irec<nrecs

        endif    ! found var in file

        if (.not.found_rec) then  ! correct record not found, try next file
          ifile = ifile+1 
          irec = 0                ! Reset irec for new file
        endif

        ierr=nf90_close (ncid)

      enddo  ! search through files: while .not. found_rec .and. ifile<=nfiles.

      if (.not. found_var_ever) then
        write(*,'(/1x,4A,I3/)') 'ERROR: find_new_record:: ',
     &              'Could not find var: ', vname, ' mynode=', mynode
        error stop 'find_new_record: var not found in forcing files'
      endif
      if (.not. found_rec) then
        write(*,'(/1x,3A,I3/)')
     &  'ERROR: find_new_record: Found variable, but ran out of time records for ',
     &                                      vname, ' mynode=', mynode
        error stop 'find_new_record: found var, but ran out of time records'
      endif

      if (first) then  ! Handle initital timestep. Need vtime < model_time. As per (A) above.

        ! find_new_rec_cycle deals internally with this
        if(cycle_length==0.) then  ! Not cyclical data
          if (irec==1) then        ! Model time between input files.
            ifile = ifile_old
          endif
          irec = irec_old
          vtime(1) = time_old
        endif

      endif ! Handle initital timestep.

      if (irec==0) then
        write(*,'(/1x,4A,I3/)') 'ERROR: find_new_record:: ',
     &  'First available forcing record is past current time for var: ',
     &     vname, ' mynode=', mynode
        error stop 'First available forcing record is past current time'
      endif

      end subroutine find_new_record  !]

! ----------------------------------------------------------------------
      subroutine find_new_record_cycle(ncid,tid,irec,vtime)
      implicit none

      ! input/output
      integer                       ,intent(in) :: ncid ! Handle of open netcdf file
      integer                       ,intent(in) :: tid  ! ID of time variable
      integer                       ,intent(out):: irec ! correct record number
      real,dimension(1)             ,intent(out):: vtime ! correct record time

      ! local
      real,dimension(1) :: tim0,tim1

      
      ierr=nf90_get_var(ncid,tid,tim0, (/1/), (/1/))        ! first time in file
      ierr=nf90_get_var(ncid,tid,tim1, (/max_rec/), (/1/))  ! last time in file

      icycle = floor(mtime -tim1)

      
      end subroutine find_new_record_cycle

      end module roms_read_write
