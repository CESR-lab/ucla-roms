      module roms_read_write
      ! ROMS specific read write functionality

      use netcdf
      use ncvars
      use scalars
      use nc_read_write
!     use param

      implicit none

#include "cppdefs.opt"

      private

      ! Type ncbry contains all required netcdf variables and input data array for frc variable.
      ! Set the name and time_name of variable to match what is in the input forcing files.
      ! Use same name as bulk_frc variable name, e.g. uwnd has nc_uwnd for netcdf vars.
      ! type (ncvs) :: nc_uwnd  = ncvs( name='uwnd',       time_name='wnd_time'  )
      type,public  :: ncbry
        character(len=20)                 :: vname     ! name of variable in input file
        character(len=20)                 :: tname     ! time variable name for variable
        real,dimension(:,:,:),allocatable :: vdata     ! need=0 to initialize all vars other than name and time_name. Slightly inefficient to set=0 though. stores raw input data for 2 times: time_A < mode
        integer                           :: coarse=0  ! Flag to indicate interpolation of the input data
        real,dimension(:,:),allocatable   :: cdata     ! Array to read in coarse surface forcing data
        integer                           :: ifile=0   ! Index to file in list of filenames
        integer                           :: irec =0   ! Record number in file
        integer         :: it1 = 1, it2 = 2              ! used to cycle between correct entries of 'data' above
        real,   dimension(2)              :: times = -99 ! [-99,-99]       ! stores 2 times that go with 'data' above
      end type ncbry

      integer, target :: test_var

!     interface set_frc_data
!       module procedure  set_frc_data_2D
!     end interface

      public :: set_frc_data,create_filename

      contains

! ----------------------------------------------------------------------
      subroutine set_frc_data(nc,var1d,var2d)
      implicit none

      ! input/outputs
      type(ncbry)        ,intent(inout) :: nc    ! derived type containing all neccesary supporting data
      real,dimension(:,:),optional,intent(out)   :: var2d ! time interpolated forcing data
      real,dimension(:)  ,optional,intent(out)   :: var1d  ! time interpolated forcing data

      ! local
      real    :: cff1,cff2  ! for time interpolations
      real    :: modtime    ! model time
      logical :: d1

      ! unpack ncbry derived type
      integer           :: it1,it2
      real,dimension(2) :: vtimes

      if (present(var1d)) then
        d1 = .true.
      else
        d1 = .false.
      endif

      modtime = tdays + 0.5*dt*sec2day    ! input data times in days. Interpolated time is based on 1/2 step time

      it1 = nc%it1
      it2 = nc%it2

      if (nc%times(it2) < modtime) then  ! need to refresh data

        if (nc%ifile == 0) then   ! the first time, we need to fill both slots
          call fill_frc_slice(nc,modtime,it1,d1)
        else ! cycle it1,it2 between 1 and 2
          it1 = mod(it1,2) + 1
          it2 = mod(it2,2) + 1
        endif

        call fill_frc_slice(nc,modtime,it2,d1)
      endif 

      nc%it1 = it1
      nc%it2 = it2

      ! Temporal interpolation

      vtimes  = nc%times
      cff1=( vtimes(it2)-modtime )/( vtimes(it2)-vtimes(it1) )
      cff2=( modtime-vtimes(it1) )/( vtimes(it2)-vtimes(it1) )

      if (d1) then
        var1d = cff1*nc%vdata(:,1,it1) + cff2*nc%vdata(:,1,it2)
      else
        var2d = cff1*nc%vdata(:,:,it1) + cff2*nc%vdata(:,:,it2)
      endif


      if (vtimes(it1)>modtime.or.vtimes(it2)<modtime) then
        error stop 'ERROR: set_frc_data :: time interpolation error'
      endif

      end subroutine set_frc_data
      
! ----------------------------------------------------------------------
      subroutine find_new_record(vname,tname,time,ifile,irec,vtime) ![
      ! Finds the index for the filename and record number of the first record
      ! with a time that is larger than the model time. It will stores
      ! the corresponding time of the record in vtime.

      ! When called for the first time, it will find the last record with a 
      ! time smaller than the model time.
      implicit none

      ! input/output
      character(len=*)              ,intent(in)   :: vname   ! name of the variable
      character(len=*)              ,intent(in)   :: tname   ! name of the corresponding time
      real                          ,intent(in)   :: time    ! target time
      integer                       ,intent(inout):: ifile   ! index for file name in list
      integer                       ,intent(inout):: irec    ! variable record entry
      real                          ,intent(out)  :: vtime   ! read in variable time

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
      integer :: vid,tid                  ! ID's of variable and its time
      real    :: cycle_length             ! cyclical data
      logical :: first                    ! first time call for this variable

      real,dimension(1)    :: ltime

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

!     if (mynode==0) print *, 'find_new_rec, trying: ',frcfile(ifile)
        found_var = .false.  ! reset for new file
        ierr = nf90_open(frcfile(ifile), nf90_nowrite, ncid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'Find_new_record, opening: ',frcfile(ifile))

        ierr = nf90_inq_varid(ncid, vname,vid)    ! Check if this file contains the variable

        if (ierr==nf90_noerr) then ! Variable found, now check if file times contain model time

          found_var      = .true.
          found_var_ever = .true.
           
          ierr = nf90_inq_varid(ncid,tname,tid)  ! Get time ID
          if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'Find_new_record, searching: ',tname)

          ! Check if the time has a cycle length
          ierr=nf90_get_att(ncid, tid, 'cycle_length', cycle_length) 
          if (ierr==nf90_noerr) then
            call find_new_record_cycle(ncid,tid,time,irec,vtime,cycle_length,first)
            found_rec = .true.
          endif

          ! Number of records in file
          ierr = nf90_inquire_variable(ncid,tid, dimids = dimids)
          ierr = nf90_inquire_dimension(ncid, dimids(1), len=max_recs)

          do while (irec<max_recs .and. .not. found_rec)                ! Search through the records until correct time is found

            irec=irec+1
            ierr=nf90_get_var(ncid,tid,ltime,(/irec/), (/1/))       ! Get vtime (nf90 needs array hence syntax)
            vtime = ltime(1)

            if (vtime>time) then   ! Correct time
              found_rec = .true.
            else                      ! Not yet
              irec_old = irec
              ifile_old = ifile     ! Only needed initially. See label (A)
              time_old = vtime
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
          vtime = time_old
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
      subroutine find_new_record_cycle(ncid,tid,time,irec,vtime,cycle_length,first)
      ! Find the first record that has vtime>time, or if it's
      ! the first time calling, the record preceeding that one.
      implicit none

      ! input/output
      integer                       ,intent(in)  :: ncid ! Handle of open netcdf file
      integer                       ,intent(in)  :: tid  ! ID of time variable
      real                          ,intent(in)  :: time ! target time
      integer                       ,intent(out) :: irec ! correct record number
      real                          ,intent(out) :: vtime! correct record time
      real                          ,intent(in)  :: cycle_length
      logical                       ,intent(in)  :: first

      ! local
      integer                       :: icycle,max_rec,ierr
      integer, dimension(1) :: dimids     ! time dimension ID
      real,allocatable,dimension(:) :: tim

      ! Number of records in file
      ierr = nf90_inquire_variable(ncid,tid, dimids = dimids)
      ierr = nf90_inquire_dimension(ncid, dimids(1), len=max_rec)
      
      allocate(tim(max_rec))
      ierr=nf90_get_var(ncid,tid,tim)

      icycle = floor((time - tim(max_rec))/cycle_length) + 1

      do irec = 1,max_rec
         vtime = tim(irec) + icycle*cycle_length
         if (vtime>time) then
           exit
         endif
      enddo

      if (first) then
        irec = irec-1;
        if (irec==0) then
          icycle = icycle-1
          irec = max_rec
        endif
        vtime = tim(irec) + icycle*cycle_length
      endif

      deallocate(tim)

      end subroutine find_new_record_cycle

! ----------------------------------------------------------------------
      subroutine fill_frc_slice(nc,modtime,it,d1)
      ! Fill a time slice of forcing data
      implicit none

      ! input/outputs
      type(ncbry)        ,intent(inout) :: nc    ! derived type containing all neccesary supporting data
      real               ,intent(in)    :: modtime
      integer            ,intent(in)    :: it
      logical            ,intent(in)    :: d1  ! 1 or 2 dimensional forcing data

      ! local
      integer :: ierr,ncid

      ! unpacking derived type ncbry
      character(len=20)  :: vname  ! name of variable in input file
      character(len=20)  :: tname  ! time variable name for variable
      integer            :: irec,ifile

      irec  = nc%irec
      ifile = nc%ifile
      vname = nc%vname
      tname = nc%tname

      call find_new_record(vname,tname,modtime,ifile,irec,
     &                         nc%times(it) )

      ierr=nf90_open(frcfile(ifile),nf90_nowrite, ncid)
      if (nc%coarse) then
        print *,'coarse to fine not yet implemented here'
        stop
!       cdata => nc_var%cdata
!       ncread(ncid,vname,cdata,(/1,1,irec/))
!       call exchange_2d
!       call coarse2fine(cdata,vdata(:,:,it))
      else
        !! Allocate vdata to be the right size for the
        !! subdomain (including whether we need a buffer or not
        if (d1) then
          call ncread(ncid,vname,nc%vdata(:,1,it),(/1,irec/))
        else
          call ncread(ncid,vname,nc%vdata(:,:,it),(/1,1,irec/))
        endif
      endif
      ierr = nf90_close(ncid)

      if (mynode == 0) 
     &      call display_read_time_to_log(vname,nc%times(it),irec )

      nc%irec  = irec
      nc%ifile = ifile

      end subroutine fill_frc_slice

! ----------------------------------------------------------------------
      subroutine handle_ierr(ierr, err_msg1, err_msg2)  ![
      ! Handle fatal ierr with error messages.
      ! err_msg1 & 2 are optional messages
      implicit none

      ! Inputs
      integer ierr
      character(len=*), optional :: err_msg1, err_msg2
      ! Local
      character(:), allocatable :: f_msg, f_msg1, f_msg2

      ! Check if error message 1 is available
      if(present(err_msg1))then
          f_msg1=err_msg1
      else
          f_msg1=''
      endif
      ! Check if error message 2 is available
      if(present(err_msg2))then
          f_msg2=err_msg2
      else
          f_msg2=''
      endif
      ! Combine optional messages:
      f_msg = trim(f_msg1) / / ' ' / / trim(f_msg2)

      write(*,'(/3x,2A/5x,A/12x,A/)') 'ERROR: read_write :: ',
     &       'netcdf ierr != 0 ', trim(f_msg), nf90_strerror(ierr)

      error stop ! Stop simulation assumed fatal error

      end subroutine handle_ierr  !]
! ----------------------------------------------------------------------
      subroutine display_read_time_to_log(var_name, var_time, var_irec)  ![
      implicit none   
        
      character(len=*), intent(in) :: var_name
      integer,          intent(in) :: var_irec
      real,             intent(in) :: var_time
          
      write(*,'(9x,A,A12,5x,A,G14.6,1x,A,I4)') 
     &  'set_frc :: ',         var_name,
     &  'input time (days) =', var_time,         ! G is number for time 
     &  'rec =',               var_irec

      end subroutine display_read_time_to_log  !]
! ----------------------------------------------------------------------
      subroutine create_filename(fname,total_rec)  ![
      ! takes a filename for output and inserts the mpi node number 
      ! (if appropriate), and the label that matches the current record
      ! number
      implicit none

      ! input/output
      character(len=*),  intent(inout) :: fname              ! desired netcdf file name
      integer,           intent(in)    :: total_rec          ! total netcdf records so far for these variables

      ! local
      integer           :: lfnm, lvar, lenstr, lstr, ierr=0
      character(len=64) :: fname_tmp                         ! taken from read_inp.F

      fname_tmp=fname
      lstr=lenstr(fname_tmp)
#if defined MPI && defined PARALLEL_FILES
      call insert_node(fname_tmp, lstr, mynode, NNODES, ierr)  ! insert MPI node numbers to file name
#endif
      fname = fname_tmp(1:lstr)
      call insert_time_index(fname,lstr,total_rec,ierr)     ! insert time index to file name

      end subroutine create_filename !]
! ----------------------------------------------------------------------
      subroutine coarse2fine(cdata,fdata)
      ! interpolates from (possibly partitioned) coarse to fine surface forcing data
      implicit none


      ! input/output
      real,dimension(:,:),intent(in) :: cdata   ! coarse data
      real,dimension(:,:),intent(out):: fdata   ! fine data

      ! local
      integer :: nx,ny
      integer :: i,j
      integer :: ic,jc
      real  :: a,b

      a = 0.75
      b = 0.25


!     call exchange_2d(cdata)

      do j=1,ny,2
        do i=1,nx,2
          fdata(i  ,j  ) = a*(a*cdata(ic,jc  )+b*cdata(ic+1,jc  ))
     &                   + b*(a*cdata(ic,jc+1)+b*cdata(ic+1,jc+1))
          fdata(i+1,j  ) = a*(b*cdata(ic,jc  )+a*cdata(ic+1,jc  ))
     &                   + b*(b*cdata(ic,jc+1)+a*cdata(ic+1,jc+1))
          fdata(i  ,j+1) = b*(a*cdata(ic,jc  )+b*cdata(ic+1,jc  ))
     &                   + a*(a*cdata(ic,jc+1)+b*cdata(ic+1,jc+1))
          fdata(i+1,j+1) = b*(b*cdata(ic,jc  )+a*cdata(ic+1,jc  ))
     &                   + a*(b*cdata(ic,jc+1)+a*cdata(ic+1,jc+1))
        enddo
      enddo


      end subroutine coarse2fine

      end module roms_read_write
