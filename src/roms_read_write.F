      module roms_read_write
      ! ROMS specific read write functionality
      !
      ! output_root_name is used as prefix of all output files
      ! Length is limited to 18 characters, because output_file_names
      ! are appended with frame and node number and are limited to 32

      use netcdf
      use scalars
      use nc_read_write

      implicit none

#include "cppdefs.opt"

      private


      ! Time is in seconds relative to a reference date
      
      ! ---- user input begins
      integer,parameter,dimension(3) :: reference_date = (/2000,1,1/)  ! year, month, day
      integer :: dt_format = 0
      ! ---- user input ends

      character(len=21) :: refdatestr
      integer :: offset = 0


      ! Type ncforce contains all required netcdf variables and input data array for frc variable.
      ! Set the name and time_name of variable to match what is in the input forcing files.
      ! Use same name as bulk_frc variable name, e.g. uwnd has nc_uwnd for netcdf vars.
      ! type (ncvs) :: nc_uwnd  = ncvs( name='uwnd',       time_name='wnd_time'  )
      type,public  :: ncforce
        character(len=20)                 :: vname     ! name of variable in input file
        character(len=20)                 :: tname     ! time variable name for variable
        real,dimension(:,:,:),allocatable :: vdata     ! currently must be exact size of sub-domain (no buffer)
        real,dimension(:,:),allocatable   :: cdata     ! Array to read in coarse surface forcing data
        integer                           :: coarse=0  ! Flag to indicate interpolation of the input data
        integer                           :: ifile=0   ! Index to file in list of filenames
        integer                           :: irec =0   ! Record number in file
        integer                           :: it1 = 1, it2 = 2              ! used to cycle between correct entries of 'data' above
        real,dimension(2)                 :: times = -99 ! [-99,-99]       ! stores 2 times that go with 'data' above
      end type ncforce

      integer, target :: test_var
      character(len=99), public  :: output_root_name
      character(len=41), public  :: git_hash

      ! index ranges for netcdf output (See init_output_indices)
      integer,public :: i0,i1,j0,j1
      integer,public :: bfx,bfy

      integer,parameter :: max_options_string=2000
      character(len=max_options_string),public::surf_forcing_strings=' '
      character(len=max_options_string),public::bc_options = ' '

      ! from old ncvars **************:
      integer, parameter, public :: max_frc_files=360
      integer, public :: max_frc, ncfrc(max_frc_files), nrrec
      ! Horizontal Grid Type Codes =  0,1,2,3 for RHO-, U-, V-, PSI-points
      integer, parameter, public :: r_var=0, u_var=1, v_var=2, q_var=3
      !            Horizontal array dimensions in netCDF files. In the case
      ! xi_rho     of MPI code with PARALLEL_FILES activated these dimensions
      ! xi_u       depend on the corresponding sizes of individual subdomains
      ! eta_rho    rather than the whole physical grid, therefore they become
      ! eta_v      live variables placed in common block and assigned values
      !            in "mpi_setup" rather tnan be parameters defined here.

#if defined MPI && defined PARALLEL_FILES
      integer, public ::  xi_rho, xi_u, eta_rho, eta_v
#else
      integer, parameter, public :: xi_rho=LLm+2, eta_rho=MMm+2,
     &                              xi_u=xi_rho-1, eta_v=eta_rho-1
#endif

      integer, parameter, public :: max_name_size=256
      character, public :: date_str*44, title*80
      character(len=max_name_size), public :: ininame, grdname, frcfile(max_frc_files)


      public :: findstr,append_date_node
      public :: set_frc_data
      public :: read_output_root_name
      public :: init_output_indices
      public :: find_new_record
      public :: nc_check_units
      public :: put_global_atts
      public :: store_string_att
      public :: create_file

      contains

! ----------------------------------------------------------------------
      subroutine init_output_indices  ![
      ! these indices also relate to reading input fields!
      !
      ! There are 2 approaches to outputting to netcdf files
      ! in the interior, we alway output (1:nx,1:ny), the boundary
      ! subdomains netcdf files are larger because they include buffers
      ! Method 1:
      ! For arrays that have buffers, we use i0:i1,j0:j1.
      ! These indices are expanded for boundary adjacent 
      ! subdomains because there, we output domain buffers

      ! Method 2:
      ! Some arrays are size (nx,ny) only. For those, we just
      ! skip the buffers for the boundary subdomain outputs.

      use dimensions
      use hidden_mpi_vars
      implicit none

      ! Method 1
      i0=1; i1=nx; j0=1; j1=ny

      if (inode==0)        i0=0
      if (inode==NP_XI-1)  i1=nx+1
      if (jnode==0)        j0=0
      if (jnode==NP_ETA-1) j1=ny+1

      ! Method 2
      bfx = 1
      bfy = 1
      if (inode==0) bfx = 2  ! Start writing at 2 to skip buffer
      if (jnode==0) bfy = 2  ! Start writing at 2 to skip buffer

      end subroutine init_output_indices  !]
! ----------------------------------------------------------------------
      subroutine set_frc_data(nc,var1d,var2d,obry)  ![
      implicit none

      ! input/outputs
      type(ncforce),               intent(inout) :: nc    ! derived type containing all neccesary supporting data
      real,dimension(:)  ,optional,intent(out)   :: var1d ! time interpolated forcing data
      real,dimension(:,:),optional,intent(out)   :: var2d ! time interpolated forcing data
      integer,optional,            intent(in)    :: obry  ! bry=1, suppress terminal info, bry=2 display boundary

      ! local
      real    :: cff1,cff2  ! for time interpolations
      real    :: modtime    ! model time
      logical :: d1
      integer :: bry

      integer           :: it1,it2
      real,dimension(2) :: vtimes

      if (present(var1d)) then
        d1 = .true.
      else
        d1 = .false.
      endif

      if (present(obry)) then
        bry = obry
      else
        bry = 0
      endif

      modtime = tdays + 0.5*dt*sec2day    ! input data times in days. Interpolated time is based on 1/2 step time

      it1 = nc%it1
      it2 = nc%it2

      if (nc%times(it2) < modtime) then  ! need to refresh data

        if (nc%ifile == 0) then   ! the first time, we need to fill both slots
          call fill_frc_slice(nc,modtime,it1,d1,bry)
        else ! cycle it1,it2 between 1 and 2
          it1 = mod(it1,2) + 1
          it2 = mod(it2,2) + 1
        endif

        call fill_frc_slice(nc,modtime,it2,d1,bry)
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
        print * ,'set_frc_data :: ', nc%vname,frcfile(nc%ifile)
        print * ,'set_frc_data :: ', vtimes(it1), vtimes(it2), modtime
        error stop 'ERROR: set_frc_data :: time interpolation error'
      endif

      end subroutine set_frc_data  !]
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

      if (ifile==0) then                                             ! first time for this variable
        first = .true.
        ifile = 1
        irec_old  = 0
        ifile_old = 0
      else
        first = .false.
      endif
      
      found_rec      = .false.
      found_var_ever = .false.

      do while (.not.found_rec .and. ifile <= max_frc)               ! while record not found & still more files to check

        found_var = .false.                                          ! reset for new file
        ierr = nf90_open(frcfile(ifile), nf90_nowrite, ncid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'Find_new_record, opening: ',frcfile(ifile))

        ierr = nf90_inq_varid(ncid, vname,vid)                       ! check if this file contains the variable

        if (ierr==nf90_noerr) then                                   ! variable found, now check if file times contain model time

          found_var      = .true.
          found_var_ever = .true.
           
          ierr = nf90_inq_varid(ncid,tname,tid)                      ! get time ID
          if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'Find_new_record, searching: ',tname)

          ! Check if the time has a cycle length
          ierr=nf90_get_att(ncid, tid, 'cycle_length', cycle_length)
          if (ierr==nf90_noerr) then
            call find_new_record_cycle(ncid,tid,time,irec,vtime,cycle_length,first)
            found_rec = .true.
          elseif (ierr==nf90_enotatt) then                           ! non-cyclical var won't have attribute
            cycle_length = -1                                        ! nf90_enotatt = attribute not found
          else                                                       ! any other netcdf error is a problem.
            call handle_ierr(ierr,'find_new_record: cycle length att')
          endif

          ! Number of records in file
          ierr = nf90_inquire_variable(ncid,tid, dimids = dimids)
          ierr = nf90_inquire_dimension(ncid, dimids(1), len=max_recs)

          do while (irec<max_recs .and. .not. found_rec)             ! Search through the records until correct time is found

            irec=irec+1                                              ! why not use ncread????
            ierr=nf90_get_var(ncid,tid,ltime,(/irec/), (/1/))        ! Get vtime (nf90 needs array hence syntax)
            vtime = ltime(1)

            if (vtime>time) then                                     ! correct time
              found_rec = .true.
            else                                                     ! not yet
              irec_old = irec
              ifile_old = ifile                                      ! only needed initially. See label (A)
              time_old = vtime
            endif

          enddo                                                      ! while not found and irec<nrecs

        endif                                                        ! found var in file

        if (.not.found_rec) then                                     ! correct record not found, try next file
          ifile = ifile+1 
          irec = 0                                                   ! Reset irec for new file
        endif

        ierr=nf90_close (ncid)

      enddo                            ! search through files: while .not. found_rec .and. ifile<=nfiles.


      if (.not. found_var_ever) then
        write(*,'(/1x,4A,I3/)') 'ERROR: find_new_record:: ',
     &              'Could not find var: ', vname, ' mynode=', mynode
        error stop 'find_new_record: var not found in forcing files'
      endif
      if (.not. found_rec) then
        write(*,'(/1x,3A,I3/)')
     &  'ERROR: find_new_record: Ran out of time records for ',
     &                               vname, ' mynode=', mynode
        error stop 
      endif

      if (first) then  ! Handle initital timestep. Need vtime < model_time. As per (A) above.

        ! find_new_rec_cycle deals internally with this
        if(cycle_length<0) then     ! not cyclical data
          if (irec==1) then         ! model time between input files.
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
      subroutine find_new_record_cycle(ncid,tid,time,irec,vtime,cycle_length,first) ![
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

      end subroutine find_new_record_cycle !]
! ----------------------------------------------------------------------
      subroutine fill_frc_slice(nc,modtime,it,d1,bry) ![
      ! Fill a time slice of forcing data
      implicit none

      ! input/outputs
      type(ncforce)      ,intent(inout) :: nc      ! derived type containing all neccesary supporting data
      real               ,intent(in)    :: modtime ! current time of the model
      integer            ,intent(in)    :: it      ! index (1 or 2) of the forcing data
      logical            ,intent(in)    :: d1      ! 1 or 2 dimensional forcing data
      integer            ,intent(in)    :: bry     ! suppress terminal output of boundary

      ! local
      integer :: ierr,ncid

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
        !! Allocate vdata to be the right size for the (DevinD - this should have been put somewhere clearer...
        !! subdomain (including whether we need a buffer or not
        if (d1) then
          call ncread(ncid,vname,nc%vdata(:,1,it),(/1,irec/))
        else
          call ncread(ncid,vname,nc%vdata(:,:,it),(/1,1,irec/))
        endif
      endif
      ierr = nf90_close(ncid)

      if (bry==0) then
        if (mynode==0)
     &    call display_read_time_to_log(vname,nc%times(it),irec )
      elseif (bry==1) then
        ! do nothing (suppress terminal output)
      elseif (bry==2) then
        if (mynode==0 .or. mynode==nnodes-1)                         ! 2 nodes catch e/w/s/n boundaries
     &    call display_read_time_to_log(vname,nc%times(it),irec )
      else
        stop 'fill_frc_slice - bry error'
      endif


      nc%irec  = irec
      nc%ifile = ifile

      end subroutine fill_frc_slice !]
! ----------------------------------------------------------------------
      subroutine display_read_time_to_log(var_name, var_time, var_irec)  ![
      implicit none   
        
      character(len=*), intent(in) :: var_name
      integer,          intent(in) :: var_irec
      real,             intent(in) :: var_time
          
      write(*,'(9x,A,A12,5x,A,G14.6,1x,A,I4)') 
     &  'set_frc :: ',         var_name,
     &  'input time (days) =', var_time,  
     &  'rec =',               var_irec

      end subroutine display_read_time_to_log  !]
! ----------------------------------------------------------------------
      subroutine read_output_root_name(keyword, kwlen)  ![
      ! read from .in file the root-name
      ! for all output netcdf result files
      ! this is called from read_inp.F
      ! The following needs to be in the roms.in input file to use this
      ! subroutine (file name needs to be 5 spaces from left margin):
      !output_root_name:
      !     rootname (<- insert desired rootname)

      implicit none

      ! inputs
      character(len=*) :: keyword
      integer          :: kwlen                            ! keyword length
      ! local
      integer           :: ierr = nf90_noerr
      integer,parameter :: input=15

      call cancel_kwd (keyword(1:kwlen), ierr)             ! remove keyword from keyword list for read_inp.F
                                                           ! to know it's accounted for
      read(input,'(A)',err=95) output_root_name            ! read & save filename root

      goto 100
                                                           ! error for read(input,'(A)',err=95)
  95  write(*,'(/1x,4A/)') '### ERROR: roms_read_write :: Cannot ',
     &                    'read entry ''', keyword(1:kwlen), '''.'
      error stop
 100  continue


      end subroutine read_output_root_name  !]
! ----------------------------------------------------------------------
      subroutine create_file(ftype,fname,nodate)  ![
      ! Creates a file using the root_name plus the extension
      ! It appends a date/time string and, if neccesary, a node number
      ! It puts the global attributes and creates an ocean_time variable
      implicit none

      ! input/output
      character(len=*), intent(in) :: ftype     ! desired netcdf file extension
      character(len=*), intent(out):: fname     ! desired netcdf file name
      logical,optional, intent(in) :: nodate    ! optional argument to skip date label and time variable

      ! local
      integer :: ierr,ncid,varid

      fname=trim(adjustl(output_root_name)) / / trim(ftype)
      if (present(nodate)) then
        call append_date_node(fname,nodate)
      else
        call append_date_node(fname)
      endif

!     if (mynode==0) then
!       print *,'creating: ',fname
!       stop
!     endif

      ierr=nf90_create(trim(fname),nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) 
     &   call handle_ierr(ierr,'create_file: '/ /trim(fname) )

      if (.not.present(nodate)) then
        varid = nccreate(ncid,'ocean_time',(/'time'/),(/0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name', refdatestr)
        ierr = nf90_put_att(ncid,varid,'units','second' )
      endif

      call put_global_atts(ncid, ierr)                     ! put global attributes in file
      if (ierr/=nf90_noerr) 
     %   call handle_ierr(ierr,'create_file: '/ / trim(fname) )

      ierr = nf90_close(ncid)
      if (mynode == 0) then
        write(*,'(7x,2A)')
     &    'created new netcdf file ', trim(fname)
      endif

      end subroutine create_file !]
! ----------------------------------------------------------------------
      subroutine coarse2fine(cdata,fdata) ![
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

      end subroutine coarse2fine !]
! ----------------------------------------------------------------------
      logical function findstr(string,pattern,istart) ![
      ! Should be in roms_read_write and then that module should be available to
      ! Tools-Roms/ as well as src/
      ! didn't do it now because roms_read_write depends on other modules so would
      ! require a bit of work (do once ncvars is removed)
      implicit none

      !input/output
      character(len=*),intent(in)  :: string        ! string
      character(len=*),intent(in)  :: pattern       ! desired pattern to find within string
      integer,optional,intent(out) :: istart

      !local
      integer :: nl,nlv,i

      nl  = len(trim(string))
      nlv = len(pattern)

      findstr = .false.
      do i = 1,nl-nlv+1
         if (string(i:i+nlv-1) == pattern) then
          findstr = .true.
          exit
         endif
      enddo

      if (present(istart)) then
        if (findstr) then
          istart=i                                  ! return string starting index
        else
          istart=0
        endif
      endif

      end function findstr !]
! ----------------------------------------------------------------------
      subroutine sec2date(time,date) ![
      ! input  time in seconds since 1970.0
      ! output date  Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
      ! Adapted from: Clive Page, Leicester University, UK.   1995-may-2
      implicit none

      ! import/export
      real,                intent(in)  :: time
      integer,dimension(6),intent(out) :: date

      ! Local
      integer :: mjday, nsecs,utime
      real    :: day

      !! number of seconds from 1970,1,1 0:00:00 to 2000,1,1 0:00:00
      if (offset==0) then
        call init_refdate
      endif

      utime = time + offset

      ! Note the MJD algorithm only works from years 1901 to 2099.
      mjday   = int(utime/86400 + 40587)
      date(1) = 1858 + int( (mjday + 321.51) / 365.25)
      day     = aint( mod(mjday + 262.25, 365.25) ) + 0.5
      date(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
      date(3) = 1 + int(mod(day,30.6))
      nsecs   = mod(utime, 86400)
      date(6) = mod(nsecs, 60)
      nsecs   = nsecs / 60
      date(5) = mod(nsecs, 60)
      date(4) = nsecs / 60

      end subroutine sec2date !]
! ----------------------------------------------------------------------
      subroutine init_refdate ![
      ! input: date  Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
      ! output  time in seconds since 1970.1,1 0:00:00
      implicit none

      ! Local
      integer :: jul0,jul,y,m,d
      real    :: day

      y = 1970
      m = 1
      d = 1

      jul0 = d-32075 + 1461*(y + 4800 + (m-14)/12)/4 +
     &       367*(m-2-((m-14)/12)*12)/12 -
     &       3*((y+4900+(m-14)/12)/100)/4

!     jul0 = floor(jul0)

      y = reference_date(1)
      m = reference_date(2)
      d = reference_date(3)

      jul = d-32075 + 1461*(y + 4800 + (m-14)/12)/4 +
     &       367*(m-2-((m-14)/12)*12)/12 -
     &       3*((y+4900+(m-14)/12)/100)/4

      offset = (jul-jul0)*3600*24

      write(refdatestr,'(A,I4,A,I0.2,A,I0.2)') 
     &    'Time since ',reference_date(1),'/'
     &,       reference_date(2),'/',reference_date(3)

      if (mynode==0) then
        print *,'Reference time used for simulation is: ',
     &           refdatestr(12:21)
      endif

             


      end subroutine init_refdate !]
! ----------------------------------------------------------------------
      subroutine append_date_node(fname,nodate) ![
      ! Insert date and node number string into filename
      ! If the nodate argument is present, only add the node number
      implicit none

      ! import/export
      character(len=*),intent(inout):: fname 
      logical,optional,intent(in)   :: nodate    ! leave out date extention

      ! Local
      integer,dimension(6) :: date
      character(len=15)  :: datestr
      character(len=6)   :: indxstr
      integer :: label
      integer,parameter :: period=900

      if (.not.present(nodate)) then
        call sec2date(time,date)
        select case (dt_format)
          case(1)  !! omit the year
            write(datestr,'(A,5I0.2)') '.',date(2:6)
          case(2)  !! omit the year and the month
            write(datestr,'(A,4I0.2)') '.',date(3:6)
          case(3)  !! omit the seconds
            write(datestr,'(A,I4,4I0.2)') '.',date(1:5)
          case(4)  !! Get old style number, based on output period
            label = int(time/period)
            write(datestr,'(A,I0.5)') '.',label
          case default
            write(datestr,'(A,I4,5I0.2)') '.',date
        end select
        fname = trim(fname) //datestr
      endif

#if defined MPI && defined PARALLEL_FILES
      if (nnodes<10) then
        write(indxstr,'(A,I0.1)') '.',mynode
      elseif (nnodes<100) then
        write(indxstr,'(A,I0.2)') '.',mynode
      elseif (nnodes<1000) then
        write(indxstr,'(A,I0.3)') '.',mynode
      elseif (nnodes<10000) then
        write(indxstr,'(A,I0.4)') '.',mynode
      elseif (nnodes<100000) then
        write(indxstr,'(A,I0.5)') '.',mynode
      elseif (nnodes>100000) then
         print *, 'ROMS is not ready for this future'
         stop
      endif
      fname = trim(fname) //trim(indxstr)
#endif

      fname = trim(fname) // '.nc'

      end subroutine append_date_node !]
!-----------------------------------------------------------------------
      subroutine nc_check_units(file_name,var_name, var_units )  ![
!     subroutine nc_check_units(ifile,var_name, var_units )
      implicit none

      ! input/output
!     integer,          intent(in) :: ifile      ! index
      character(len=*), intent(in) :: file_name
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: var_units

      ! local
      integer           :: ierr,ncid
      integer           :: var_id
      character(len=20) :: file_var_units ! variable units in forcing file

!     ierr = nf90_open(frcfile(ifile), nf90_nowrite, ncid)
      ierr = nf90_open(file_name,nf90_nowrite,ncid)
      if (ierr/=nf90_noerr) 
!    &    call handle_ierr(ierr,'check_units, opening: ',frcfile(ifile))
     &    call handle_ierr(ierr,'check_units, opening: ',file_name)

      ierr = nf90_inq_varid(ncid, var_name, var_id)                    
      if (ierr/=0) call handle_ierr(ierr,'check_units:: cant find var')
      ierr = nf90_get_att(ncid, var_id, 'units', file_var_units)        
      if (ierr/=0) call handle_ierr(ierr,'check_units:: cant find units')

      if (trim(var_units) /= trim(file_var_units)) then
        write(*,'(/3x,3A/5x,3A/)') 'ERROR:: check_units: ',
     &   trim(file_var_units),' for var: ', 
     &   trim(var_name),',  Units must be: ', trim(var_units)
        error stop
      end if
      ierr = nf90_close(ncid)

      end subroutine nc_check_units  !]
! ---------------------------------------------------------------------
      subroutine put_global_atts (ncid, ierr)  ![

      use param
      use eos_vars
      use scalars
      use scoord
      use strings
      use netcdf

      implicit none
      integer :: ncid, ierr, dimid


#if defined MPI && defined PARALLEL_FILES
      integer :: ibuff(4), npart

      npart=4                                    ! Put global attribute 'partition' which identifies subdomain
      ibuff(1)=mynode                            ! within the processor grid individually for each file.
      ibuff(2)=nnodes
      if (WESTERN_MPI_EDGE) then
        ibuff(3)=iSW_corn+iwest
      else
        ibuff(3)=iSW_corn+iwest+1
      endif
      if (SOUTHERN_MPI_EDGE) then
        ibuff(4)=jSW_corn+jsouth
      else
        ibuff(4)=jSW_corn+jsouth+1
      endif

      ! Add 4 main x and y dimensions so that ncjoin will work,
      ierr=nf90_def_dim(ncid,'xi_rho', xi_rho, dimid)
      ierr=nf90_def_dim(ncid,'xi_u',   xi_u,   dimid)
      ierr=nf90_def_dim(ncid,'eta_rho',eta_rho,dimid)
      ierr=nf90_def_dim(ncid,'eta_v',  eta_v,  dimid)

      ierr=nf90_put_att(ncid, nf90_global, 'partition', ibuff)
#endif

      ierr=nf90_put_att(ncid, nf90_global, 'title',title)         ! Title of the run

#ifndef ANA_GRID
      ierr=nf90_put_att(ncid, nf90_global, 'grid_file',grdname)
#endif
#ifndef ANA_INITIAL
      ierr=nf90_put_att(ncid, nf90_global, 'init_file',ininame)
#endif

      ierr=nf90_put_att(ncid, nf90_global, 'ntimes', ntimes)         ! Time stepping parameters
      ierr=nf90_put_att(ncid, nf90_global, 'ndtfast',ndtfast)
      ierr=nf90_put_att(ncid, nf90_global, 'dt', dt)
      ierr=nf90_put_att(ncid, nf90_global, 'dtfast', dtfast)

! Vertical coordinate type, control parameters "theta_s", "theta_b",
! "hc", and stretching curves "Cs_w", "Cs_r" at W- and RHO-points.

#ifdef SOLVE3D

      ierr=nf90_put_att(ncid, nf90_global, 'theta_s', theta_s)
      ierr=nf90_put_att(ncid, nf90_global, 'theta_b', theta_b)
      ierr=nf90_put_att(ncid, nf90_global, 'hc', hc)
      ierr=nf90_put_att(ncid, nf90_global, 'Cs_w', Cs_w)
      ierr=nf90_put_att(ncid, nf90_global, 'Cs_r', Cs_r)

! Boussinesq reference density and EOS coefficients for linear EOS.

      ierr=nf90_put_att(ncid, nf90_global, 'rho0', rho0)
      ierr=nf90_put_att (ncid, nf90_global, 'rho0_units', 'kg/m^3')
# ifndef NONLIN_EOS
      ierr=nf90_put_att(ncid,nf90_global,'Tcoef', Tcoef)
      ierr=nf90_put_att(ncid,nf90_global,'Tcoef_units','deg Celsius-1')

      ierr=nf90_put_att(ncid, nf90_global, 'T0', T0)
      ierr=nf90_put_att (ncid, nf90_global, 'T0_units', 'deg Celsius')

#  ifdef SALINITY
      ierr=nf90_put_att(ncid, nf90_global, 'Scoef', Scoef)
      ierr=nf90_put_att (ncid, nf90_global, 'Scoef_units', 'PSU-1')

      ierr=nf90_put_att(ncid, nf90_global, 'S0', S0)
      ierr=nf90_put_att (ncid, nf90_global, 'S0_units', 'PSU')
#  endif
# endif
#endif /* SOLVE3D */

#ifdef UV_VIS2
      ierr=nf90_put_att(ncid,nf90_global, 'visc2', visc2)            ! Horizontal viscosity and mixing coefficients.
      ierr=nf90_put_att (ncid,nf90_global, 'visc2_units', 'm^2/s')

      ierr=nf90_put_att(ncid, nf90_global, 'gamma2', gamma2)         ! Lateral slipperiness parameter
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      ierr=nf90_put_att(ncid,nf90_global, 'tnu2', tnu2)
      ierr=nf90_put_att (ncid,nf90_global, 'tnu2_units', 'm^2/s')
# endif
# ifdef MY25_MIXING
#  ifdef Q2_DIF2
      ierr=nf90_put_att(ncid,nf90_global, 'q2nu2', q2nu2)
      ierr=nf90_put_att (ncid,nf90_global, 'q2nu2_units', 'm^2/s')
#  endif
#  ifdef Q2_DIF4
      ierr=nf90_put_att(ncid,nf90_global,'q2nu4', q2nu4)
      ierr=nf90_put_att (ncid,nf90_global,'q2nu4_units', 'm^4/s')
#  endif
# endif
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined PP_MIXING\
                         && !defined MY2_MIXING && !defined MY25_MIXING

      ierr=nf90_put_att(ncid, nf90_global, 'Akv_bak', Akv_bak)       ! Background vertical viscosity tracer mixing
      ierr=nf90_put_att (ncid, nf90_global, 'Akv_bak_units', 'm^2/s')

      ierr=nf90_put_att(ncid, nf90_global, 'Akt_bak', Akt_bak)
      ierr=nf90_put_att (ncid, nf90_global, 'Akt_bak_units', 'm^2/s')
# endif
# ifdef MY25_MIXING
      ierr=nf90_put_att(ncid, nf90_global, 'Akq_bak',Akq_bak)
      ierr=nf90_put_att (ncid, nf90_global, 'Akq_bak_units', 'm^2/s')
# endif
#endif

#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined M3_FRC_BRY
      ierr=nf90_put_att(ncid, nf90_global, 'ubind', ubind) ! Binding velocity for external data forcing at open boundaries
      ierr=nf90_put_att (ncid, nf90_global, 'ubind_units', 'm/s')
#endif

#ifdef SPONGE
      ierr=nf90_put_att(ncid, nf90_global, 'v_sponge', v_sponge)     ! Maximum viscosity in sponge layer
      ierr=nf90_put_att (ncid, nf90_global, 'v_sponge_units','m^2/s')
#endif

      ierr=nf90_put_att(ncid,nf90_global, 'rdrg', rdrg)              ! Bottom drag parameters.
      ierr=nf90_put_att (ncid,nf90_global, 'rdrg_units', 'm/s')

      ierr=nf90_put_att(ncid,nf90_global, 'rdrg2',rdrg2)
      ierr=nf90_put_att (ncid,nf90_global, 'rdrg2_units',
     &                                      'nondimensional')
#ifdef SOLVE3D
      ierr=nf90_put_att(ncid,nf90_global, 'Zob', Zob)
      ierr=nf90_put_att (ncid,nf90_global, 'Zob_units', 'm')
#endif

      ierr=nf90_put_att (ncid, nf90_global, 'SRCS', srcs)  ! List of source-code file names and activated CPP-switches
      ierr=nf90_put_att (ncid, nf90_global, 'CPPS', cpps)

      ! string set added to get values from modules without struggling with circular references
      ! as those modules use roms_read_write so this module can't use them.
      ! rather the string is filled in the module's init and can be written here to all files.
      ierr=nf90_put_att (ncid, nf90_global, 'surf_forcing_strings', surf_forcing_strings)
      ierr=nf90_put_att (ncid, nf90_global, 'bc_options', bc_options)

      call add_git_hash(ncid)

      end subroutine put_global_atts  !]
! ---------------------------------------------------------------------
      subroutine store_string_att( char_string, string_add )  ![
      ! Add settings that affect the solution to a string
      ! which will be added to global attributes of all nc files
      implicit none

      ! input
      character(len=*),intent(inout) :: char_string
      character(len=*),intent(in)    :: string_add

      ! local
      integer :: is, ie

      is=len_trim( char_string )
      is=is+2
      ie=len_trim( string_add )
      ie=ie+is-1

      if (ie > max_options_string) then
        print *, 'ERROR roms_read_write:: ',
     &    'store_string_att character full=', trim(string_add)
        stop 'ROMS stopped'
      endif

      char_string(is:ie) = string_add

      end subroutine store_string_att  !]
! ---------------------------------------------------------------------
      function  log2str(arg) ![
      ! Translate a logical into a string
      implicit none

      ! input
      logical,intent(in) :: arg

      character(len=5) :: log2str

      if (arg) then
        log2str = 'True '
      else
        log2str = 'False'
      endif


      end function log2str  !]
! ---------------------------------------------------------------------

      end module roms_read_write
