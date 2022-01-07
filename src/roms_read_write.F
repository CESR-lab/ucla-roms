      module roms_read_write
      ! ROMS specific read write functionality

      use netcdf
!     use param

      implicit none
      private

      interface ncwrite
        module procedure  ncwrite_1D, ncwrite_2D, ncwrite_3D
      end interface

      public :: ncwrite
      public :: nccreate

      contains

!----------------------------------------------------------
      subroutine ncwrite_1D(ncid,vname,dat,record)
        implicit none

        ! Input
        integer          :: ncid
        character(len=*) :: vname
        integer,optional :: record

        real,dimension(:) :: dat

        ! local
        integer              :: varid,ierr
        integer,dimension(1) :: dims
        integer,dimension(2) :: start,sizes

        dims = shape(dat)
        ierr = nf90_inq_varid(ncid,vname,varid)

        if (present(record)) then
          start = (/1,record/)
          sizes = (/dims, 1/)
          ierr = nf90_put_var(ncid,varid,dat,start,sizes) 
        else
          ierr = nf90_put_var(ncid,varid,dat)
        endif

      end subroutine ncwrite_1D

!----------------------------------------------------------
      subroutine ncwrite_2D(ncid,vname,dat,start)
        implicit none

        ! Input
        integer          :: ncid
        character(len=*) :: vname
        integer,dimension(:),optional :: start

        real,dimension(:,:) :: dat


        ! local
        integer              :: varid,ierr
        integer,dimension(2) :: dims
        integer,dimension(3) :: sizes

        dims = shape(dat)
        ierr = nf90_inq_varid(ncid,vname,varid)
        if (ierr.ne.nf90_noerr) then
          print *, trim(nf90_strerror(ierr))
          stop
        endif

        if (present(start)) then
          sizes = (/dims, 1/)
          ierr = nf90_put_var(ncid,varid,dat,start,sizes) 
        else
          ierr = nf90_put_var(ncid,varid,dat)
        endif
        if (ierr.ne.nf90_noerr) then
          print *, trim(nf90_strerror(ierr))
          stop
        endif

      end subroutine ncwrite_2D

!----------------------------------------------------------
      subroutine ncwrite_3D(ncid,vname,dat,start)
        implicit none

        ! Input
        integer          :: ncid
        character(len=*) :: vname
        integer,dimension(:),optional :: start

        real,dimension(:,:,:) :: dat

        ! local
        integer              :: varid,ierr
        integer,dimension(3) :: dims
        integer,dimension(4) :: sizes

        dims = shape(dat)
        ierr = nf90_inq_varid(ncid,vname,varid)
        if (ierr.ne.nf90_noerr) then
          print *, trim(nf90_strerror(ierr))
          stop
        endif

        if (present(start)) then
          sizes = (/dims, 1/)
          ierr = nf90_put_var(ncid,varid,dat,start,sizes) 
        else
          ierr = nf90_put_var(ncid,varid,dat)
        endif
        if (ierr.ne.nf90_noerr) then
          print *, trim(nf90_strerror(ierr))
          stop
        endif

      end subroutine ncwrite_3D

! ----------------------------------------------------------------------
      subroutine set_forcing_2D ( 
     &                var_name,      var_time_name, !  Names of var and its time in the netcdf file
     &                var_data,      var,           !  Arrays for the variable
     &                var_time,                     !  Input times
     &                var_file_indx, var_irec,      !  File indx & Current input record
     &                it1, it2 )                    !  Time index placeholders

      implicit none

      ! input/outputs
      character(len=*), intent(in) :: var_name         ! Name of the desired variable in the file
      character(len=*), intent(in) :: var_time_name    ! Name of the corresponding time variable

      real,dimension(:,:,2),intent(inout) :: var_data  ! 2 time levels of data from file to interpolate
      real,dimension(2)    ,intent(inout) :: var_time  ! Corresponding times of data from file to interpolate
      real,dimension(:,:)  ,intent(out)   :: var       ! Forcing data, interpolated to the right time
      integer              ,intent(inout) :: ifile     ! index in the list of filenames
      integer              ,intent(inout) :: irec      ! record in the file
      integer              ,intent(inout) :: it1, it2  ! Indices to the 2 instances of data and times

      ! local
      integer :: tmp, i, j
      real    :: cff1,cff2  ! for time interpolations
      real    :: tmid_days

      tmid_days = tdays + 0.5*dt*sec2day    ! input data times in days. Interpolated time is based on 1/2 step time

      if (var_time(it2) < tmid_days) then  ! need to refresh data

        if (var_file_indx == 0) then        ! initially only

          call find_new_record( var_name, var_time_name, tmid_days,
     &                          ifile,irec, var_time(it1) )

          call read_forcing( ifile,irec,var_name,var_time(it1),
     &                            var_data(:,:,it1) )

          if (mynode == 0) call display_read_time_to_log( var_name, var_time(it1), irec )

        else
          tmp = it1
          it1 = it2 ! if it1 = 1, it now equals 2, and vice-versa
          it2 = tmp
        endif

        call find_new_record( var_name, var_time_name, tmid_days,
     &                        ifile, irec, var_time(it2) )

        call read_forcing( ifile,irec,var_name,var_time(it2),
     &                          var_data(:,:,it2) )

        if (mynode == 0) call display_read_time_to_log( var_name, var_time(it2), irec )

      endif ! -> END OF UPDATING VARIABLE INPUT DATA

      ! Temporal interpolation

      cff1=( var_time(it2)-tmid_days ) / ( var_time(it2)-var_time(it1) )
      cff2=( tmid_days-var_time(it1) ) / ( var_time(it2)-var_time(it1) )

      if (cff1.ge.0. .and. cff2.ge.0.) then

        var = cff1*var_data(:,:,it1) + cff2*var_data(:,:,it2)

      else
        write(*,'(/1x,4A/3(1x,A,F16.10)/)') 'ERROR: set_small_array',
     &  ':: Model time outside bounds of variable (times):',
     &     var_name, var_time_name,
     &                    'start =',    var_times(it1),
     &                    'tmid_days=', tmid_days,
     &                    'end =',      var_times(it2)
        error stop 'ERROR: set_forcing :: time interpolation error'
      endif

      end subroutine set_forcing_2D 

      end module roms_read_write
