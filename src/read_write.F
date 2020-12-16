      module read_write

      ! =================================================
      ! Contains functions and subroutine for interacting
      ! with netcdf file input and output
      ! =================================

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020)
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use netcdf
      implicit none

      ! Includes
      ! ========

#include "cppdefs.h"
      ! param.h needed for GLOBAL_2D_ARRAY
!#include "param.h" - Aug 7th get from 'dimensions'
      ! scalars.h - needed for may_day_flag in error handling
#include "scalars.h"

      private

      ! Root name for prefix of all output files, chosen as 26 characters
      ! as less than 32 characters of output_file_names as nodes & time
      ! numbers need to be appended to name, so root must be shorter.
      character(len=26), public :: output_root_name

      ! DevinD - this clashes with dimensions which contains ncvars.h
      ! where these variables live. But soon ncvars should be moved from
      ! dimensions and put into this module or another plan be made.
      ! Horizontal Grid Type Codes =  0,1,2,3 (RHO-, U-, V-, PSI-points)
      ! for reading and writing to correct grid type
      integer, parameter, public :: rp_var=0, up_var=1, vp_var=2, qp_var=3
      real :: cycle_length ! Deal with cyclical forcing files

      ! Public subroutines
      public find_rec_new_mod
      public nc_define_var
      public nc_write_var
      public nc_read_var
!      public nc_check_units
      public set_frc_var_tile
      public set_small_arrays
      public ncdf_create_file
      public nc_write_time
      public read_output_root_name
      public ncdf_read_mod
      public handle_ierr

      public ncdf_read_coarser_grid

      contains

! ----------------------------------------------------------------------

      subroutine find_rec_new_mod(model_time,v_name,vt_name,ifile,irec)

      ! =================================================
      ! Find netcdf file (ifile) and needed record (irec)
      ! =================================================

      use netcdf

! Find the first time record number "irec" and its forcing file
! "ifile" for variable "var", with accompanying time "var_time". 
! for which var_time > ocean_time.
!
! When called for the first time, it will the last record for which
! var_time < ocean_time
!
! On entry, ifile and irec should contain results from the last call to
! find_record for this variable, or zeros if called for the first time

!------------    

!         frcfile    list of forcing file names (from module dimensions)
!         time       model time (days/secs to match input data time)
!         var_name   name of variable in file
!         var_time   name of corresponding time variable in file

! in/out: ifile      index of correct file in frcfile
!         irec       record number in file

      implicit none

      ! input/output
      ! ------------
      ! model_time must already be converted if necessary
      ! to match data format of seconds or days
      real            ,intent(in)   :: model_time
      character(len=*),intent(in)   :: v_name,vt_name
      integer         ,intent(inout):: ifile,irec

      ! local
      real,dimension(1):: vtime ! Read in variable time
      integer :: found_var,found_rec
      integer :: ierr,v_id,ncid,vt_id
      integer :: nfiles,nrecords
      integer :: ifile_old,irec_old
      integer :: dimid
      logical :: first

      ! When called for the first time (ifile==0) for this variable, this 
      ! routine will return the last record with vtime < model_time.

      if (ifile==0) then 
        first    =.true.
        ifile    = 1
        irec_old = 0
      else
        first    =.false.
      endif

      found_var = 0
      found_rec = 0

      nfiles = max_frc !! max_frc is the number of files in the list frcfile

      do while ((found_rec==0).and.(ifile<=nfiles))

        if (mynode==0) then ! ''' 3 of these ' writes in terminal 'v_name'
          write(*,'(1x,4A)') 'Searching for variable ''', ! 1x makes 1 space
     &                v_name, ''' in file ',frcfile(ifile) ! 4A says join 4 strings
        endif

        ierr = nf90_open(frcfile(ifile), nf90_nowrite, ncid)
        if (ierr.ne.0) then ! This is not supposed to happen
          if (mynode==0) then
            print *,ifile,frcfile(ifile)
            stop 'FNR: Forcing file not found'
          endif
        endif

        ! Check if this file contains the variable
        ierr = nf90_inq_varid(ncid, v_name, v_id)
        if (ierr.eq.0) then ! Variable found, now check times
          if (mynode==0) then
            write(*,'(1x,4A)') 'Found variable ''',v_name,
     &         ''' in file ',frcfile(ifile)
          endif
          found_var = 1
          irec = 1

          ! Find the size of a dimension in netcdf file
          ierr=nf90_inq_dimid(ncid, vt_name, dimid)
          ierr=nf90_inquire_dimension(ncid, dimid, len=nrecords)

          do while (irec<=nrecords.and.found_rec==0)
           

            ierr = nf90_inq_varid(ncid, vt_name, vt_id) ! Get time ID (vt_id)
            ! nf90_get_var needs an arry for last 2 terms hence syntax
            ierr = nf90_get_var(ncid, vt_id, vtime, (/irec/), (/1/))


            ! hard code conversion of input time in days to seconds
!             vtime = vtime * 86400
            if (vtime(1)>model_time) then
              ! Causes heisenbug - print vtime
!              if(mynode==0) print *,'vtime=',vtime ! Debug
              found_rec = 1
            else
              irec_old = irec
              ifile_old = ifile
              irec = irec+1
            endif
          enddo ! while not found and irec<nrecs

        else ! Var not found in file
          if (mynode==0) then
            write(*,'(1x,4A)') 'NOT found variable ''',v_name,
     &                          ''' in file:',frcfile(ifile)
          endif
        endif

        if (found_rec == 0) then
          !! go to next file if ifile<nfiles
          ifile = ifile+1
        endif

      enddo !while ifile<=nfiles

      if (first) then
        ifile = ifile_old
        irec  = irec_old
      endif

      if (found_var==0) then
!          write(*,'(/1x,2A/)')  '### ERROR: find_record: Could not find variable ', var_name
        stop 'variable not found in forcing files'
      endif
      if (found_rec==0) then
!          write(*,'(/1x,2A/)')  '### ERROR: find_record: Ran out of time records for ', var_name
        stop 'ran out of time records'
      endif
      if (irec==0) then
        stop 'First available forcing record is past current time'
      endif

      end subroutine find_rec_new_mod


!-----------------------------------------------------------------------

      subroutine set_frc_var_tile (
     &                istr,iend,     jstr,jend,     ! Tile bounds indices
     &                var_name,      var_time_name, ! Text names
     &                var_data,      var,           ! Variable arrays
     &                var_times,     var_grd_type,  ! Input times & Grid type
     &                var_file_indx, var_irec,      ! File indx & Current input record
     &                it1, it2,      coarse )       ! Time index placeholders & interp flag

      ! ===========================================
      ! WRAPPER TO READ & SET SINGLE FORCE VARIABLE
      ! VARIABLE FROM NETCDF INPUT FILE
      ! ===============================

      ! Bulk forcing variables only interpolated in time, not converted
      ! in any other way currently, so can all use the same, generic
      ! set function

      ! INPUT DATA - TIME MUST BE IN DAYS

      implicit none

      ! input/outputs
      integer,intent(in) :: istr, iend, jstr, jend
      character(len=*)   :: var_name
      character(len=*)   :: var_time_name

      real :: var_data(GLOBAL_2D_ARRAY,2) ! Read in raw data for 2 times to interpolate
      real :: var(GLOBAL_2D_ARRAY)
      real, dimension(2) :: var_times ! Store input file record times

      integer,intent(in) :: var_grd_type
      integer :: var_file_indx, var_irec ! File indx & Current input record
      integer :: it1, it2 ! Placeholders for 2 read in times
      integer, intent(in) :: coarse ! Whether function is interpolated from coarse data
      ! Try make coarse set to default value if not included using
      ! OPTIONAL argument

      ! local
      integer :: tmp,i,j
      real cff1,cff2 ! for time interpolations
      real tmid_days

! Needed for istrR,iendR, etc
# include "compute_extended_bounds.h"


      ! note: var_times(it2) > var_times(it1), but values of it1
      ! and it2 change between each reading of values so that only one
      ! value is read in.

      ! Interpolated time is based on 1/2 step time
      tmid_days = tdays + 0.5*dt*sec2day ! Input data times in days

      ! UPDATE LATEST READ IN VARIABLE?
      ! ===============================

      ! Need var_times(it1) < time < var_times(it2), thus:
      if (var_times(it2) < tmid_days) then  ! We need to refresh data

        ! FIRST TIMESTEP
        ! --------------
        ! var_times(it1) & var_times(it2) declared < 0 at top their module:
!        if (var_times(it1) < 0) then ! Only first timestep
        if (FIRST_TIME_STEP) then ! Only first timestep

          ! Reading is only done by master hence C$OMP MASTER
C$OMP MASTER

          ! ifile and irec are in their module pre-amble, initially as 0
          ! Find the force file (file_indx) and the record (irec) for
          ! the last record with var_times < model_time
          ! First time, we need to read an extra slice of data
          ! so this routine provides irec = irec - 1 for first time.
          call find_new_record( var_name, var_time_name, tmid_days,
     &                          var_file_indx, var_irec )

          ! For the first time irec = irec - 1, as above, such that
          ! closest record to var_times(it1) < model_time.
          ! read in raw input variables and
          ! read and set input time in var_times(it1) to days (units)
          call read_var_frc( var_file_indx, var_irec,
     &                       var_data(GLOBAL_2D_ARRAY,it1),
     &                       var_name, var_grd_type, var_times(it1),
     &                       var_time_name, coarse )

C$OMP END MASTER
C$OMP BARRIER
          ! Barrier needed because master processor must have read in
          ! all data before child processors try to access read in data

          ! Confirm to terminal success reading of variable
          ! show time of read in data
          if (mynode == 0) then
            write(*,'(4x,2A,5x,A,G14.6,1x,A,I4)')
     &               'set_frc_var_tile :: read var: ', var_name,
     &               'var_input_time (days) =',var_times(it1),
     &               'rec =',var_irec MYID
          endif

!          ! Convert read inputs into useable, rotated & rmask'ed variables
!          ! for the first, earlier read in time variables(time) < time.
!          call set_var_inp2vars (istr,iend,jstr,jend,it1)

        ! ALL OTHER TIMESTEPS
        ! -------------------
        ! Every time step switch indices such that
        ! var_data(i,j,1) & var_times(1) -> earlier && var_data(i,j,2) & var_times(2) -> later: step t=x
        ! var_data(i,j,1) & var_times(1) -> later   && var_data(i,j,2) & var_times(2) -> earlier: step t=x+1
        ! and so on. The benifit is only need to change one variable each time.
        else

          tmp = it1
          it1 = it2 ! If it1 = 1, it now equals 2, and vice-versa
          it2 = tmp ! If it2 = 2, it now equals 1, and vice-versa

        endif ! End of first step

        ! Reading is only done by master hence C$OMP MASTER
C$OMP MASTER

        ! Find the force file (file_indx) and the record (irec) for
        ! the first record with var_time > model_time
        call find_new_record( var_name, var_time_name, tmid_days,
     &                        var_file_indx, var_irec )

        ! Read in the first record with var_time > model_time
        ! Remember, it2 alternates between index 1 & 2 for each reading
        ! read in raw input variables and
        ! read and set input time in var_times(it2) to days (units)
        call read_var_frc( var_file_indx, var_irec,
     &                     var_data(GLOBAL_2D_ARRAY,it2),
     &                     var_name, var_grd_type, var_times(it2),
     &                     var_time_name, coarse )

C$OMP END MASTER
C$OMP BARRIER

        ! Confirm to terminal reading of variables and
        ! time for read in data
        if (mynode == 0) then                  ! Text formatting:
          write(*,'(4x,2A,5x,A,G14.6,1x,A,I4)') ! 4x is spaces, A is string. / at end adds blank line
     &               'set_frc_var_tile :: read var: ', var_name, ! 5x is 5 spaces, A is string
     &               'var_input_time (days) =',var_times(it2), ! G is number for time
     &               'rec =', var_irec MYID ! I is integer
        endif

!        ! Convert read inputs into useable, rotated & rmask'ed variables
!        ! for the latter read in time variables(time) > time.
!        call set_wec_inp2vars (istr,iend,jstr,jend,it2)

      endif ! -> END OF UPDATING VARIABLE INPUT DATA

      ! Temporal interpolation
      ! ======================

      ! Set coefficients for interpolation.
      ! Time difference in seconds. Not yet ratios of 1.
      ! cff1 with it2 because if time=it1 then difference=0, but want cff1 = 1.
      cff1=( var_times(it2)-tmid_days ) / ( var_times(it2)-var_times(it1) )
      cff2=( tmid_days-var_times(it1) ) / ( var_times(it2)-var_times(it1) )

      ! Do r-masking? Jeroen didn't necessarily want to do that.
      ! Could be here or at set_inp2vars but might be unessary extra loop
      ! when just * rmask less of an overhead?

      ! Check the interpolation coefficients such that modeltime is
      ! bounded by read in variable times.
      if (cff1.ge.0. .and. cff2.ge.0.) then

        do j=jstrR,jendR
          do i=istrR,iendR

            var(i,j)=cff1*var_data(i,j,it1)+cff2*var_data(i,j,it2)

          enddo
        enddo

      ! Throw warning if interpolation times do not bound model time.
      elseif (ZEROTH_TILE) then
        write(*,'(/1x,4A/3(1x,A,F16.10)/)') 'WARNING: set_frc_var_tile',
     &  ':: Model time outside bounds of variable (times):',
     &     var_name, var_time_name,
     &                    'start =', var_times(it1),
     &                    'tmid_days =',     tmid_days,
     &                    'end =',   var_times(it2)
!        ierr=ierr+1
        error stop 'set_frc_var_tile :: time interpolation error '
      endif


      end subroutine set_frc_var_tile
! ----------------------------------------------------------------------

      subroutine find_new_record( var_name, var_time_name, time_model,
     &                            var_file_indx, irec )

      ! ====================================
      ! FIND NETCDF FILE CONTAINING VARIABLE
      ! STORE VAR_FILE_INDX (File record)
      ! STORE irec (needed time record)
      ! ===============================

      ! Input: var_name

      ! Find: if var contained in file

      ! Check: see if model time is bounded by model

        ! -> yes: save record (irec) and save var_file_indx

        ! ->  no: try another file

      ! Limitations
      ! -----------
      ! 1) Files for a variable must be in chronological order in input
      !    file list of forcing files
      ! 2) Files for the same variable do not need to be grouped together
      !    but it is slightly more efficient if they are.
      ! 3) Algorithm does not detect a missing year of data.
      !    E.g. if modeling 2005-2007, and annual data for 2006 is missing,
      !    then times will just be interpolated between end of 2005 & start 2007.

      implicit none


      ! input/output
      ! ------------
      character(len=*),intent(in)   :: var_name, var_time_name
      real            ,intent(in)   :: time_model ! Model time same units as data
      integer         ,intent(inout):: var_file_indx
      integer         ,intent(inout):: irec ! Variable record entry

      ! local
      ! -----
      real,dimension(1) :: vtime ! Read in variable time
      integer :: found_var, found_rec, found_var_log
      integer :: ncid
      integer :: ierr
      integer :: nfiles ! Number of forcing files
      integer :: irec_old, file_indx_old
      integer :: var_max_recs ! Total variable records in file
      integer :: v_id, dimid  ! variable ID & dimension ID
      integer :: var_tid      ! Variable time id
      logical :: first        ! Check if first time finding an entry
!      real    :: cycle_length ! Deal with cyclical forcing files

      ! Routine
      ! =======

      ! First timestep (var_file_indx==0) for each variable
      if (var_file_indx==0) then
        first    = .true. ! used to set irec = irec-1 for first read
        irec_old = 0
        file_indx_old = 0
        var_file_indx = 1
      else
        first = .false.
      endif

      found_rec = 0
      found_var_log = 0 ! Track if variable was found in previous files

      ! max_frc = number of files in list 'frcfile' of force file names
      nfiles = max_frc

      ! Continue searching while record not found & still more files to check
      do while ((found_rec==0).and.(var_file_indx <= nfiles))

        found_var = 0 ! Reset variable to not found

!        if (mynode==0) then ! ''' 3 of these ' writes in terminal 'var_name'
!          write(*,'(1x,4A)') 'Searching for variable ''', ! 1x makes 1 space
!     &     var_name, ''' in file ',frcfile(var_file_indx) ! 4A says join 4 strings
!        endif

        ! Try to open force file
        ierr = nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)
        if (ierr.ne.nf90_noerr) then ! Force file not found. Stop.
          print *,frcfile(var_file_indx),'mynode=',mynode
          error stop 'Find_new_record: Forcing file not found' ! Any node not just 0
        endif

        ! Check if this file contains the variable
        ierr = nf90_inq_varid(ncid, var_name, v_id)

        ! Variable found, now check if file times contain model time
        if (ierr.eq.nf90_noerr) then

!          if (mynode==0) then
!            write(*,'(3x,4A)') 'Found variable ''',var_name,
!     &           ''' in file ',frcfile(var_file_indx)
!          endif
          found_var = 1
          ! Needed to deal with found_var = 0 further above, so that
          ! If variable is found in a file but time is wrong, and it
          ! gets to last file without finding variable, that error msg correct
          found_var_log = found_var_log + found_var

          ! Get time dimension ID
          ierr=nf90_inq_dimid(ncid, var_time_name, dimid)
          ! Find the size of dimensions (var_max_recs)
          ! I.e. total number of records per variable, using time dimension
          ierr=nf90_inquire_dimension(ncid, dimid, len=var_max_recs)

          ! Find variable ID for time variable
          ierr = nf90_inq_varid(ncid, var_time_name, var_tid) ! Get time ID (var_tid)

          ! Reset and get cycle length (if available for cyclical forcing) for new variable
          cycle_length=0.
          ierr=nf90_get_att(ncid, var_tid, 'cycle_length', cycle_length) ! Will continue if not available
          ! Set irec back to 0 for cyclical loading if irec==var_max_recs
          ! else won't be able to enter do while loop below, and must stay on same file.
          if(cycle_length/=0. .and. irec==var_max_recs) irec=0

          ! Search through all time records until correct time is found
          do while (irec<var_max_recs.and.found_rec==0)
          ! Had to change the do while because irec still needs +1
          ! therefore irec < var_max_recs not <=

            ! Advance irec. If at start then irec=0 -> irec=1 here.
            ! Also, reset irec = 0 below at label (B) if record not found in file
            ! otherwise irec = irec+1 so doesn't start from 1, but
            ! continues from previous irec.
            ! E.g. previous record was irec=3, now check time of irec=4
            irec = irec + 1

            ! Get vtime - variable time for irec
            ! nf90_get_var needs an array for last 2 terms hence syntax
            ierr = nf90_get_var(ncid, var_tid, vtime, (/irec/), (/1/))
!            if(mynode==0) print*,'vtime(1)=',vtime(1),' time_model=',time_model ! debug

            ! Fix time for cyclical forcing file if applicable
            if(cycle_length/=0.) then
              call cycle_length_handling(vtime(1), ncid, var_tid,
     &                                var_time_name, irec )
            endif

            if (vtime(1)>time_model) then ! Correct time
              found_rec = 1

            else ! Wrong time
              irec_old = irec
              ! file_indx_old is only needed for very first timestep
              file_indx_old = var_file_indx ! see label (A)

            endif

          enddo ! while not found and irec<nrecs

        ! Variable not found in file
        else
!          if (mynode==0) then
!            write(*,'(3x,4A)') 'NOT found variable ''',var_name,
!     &                            ''' in file ',frcfile(var_file_indx)
!          endif
        endif

        ! Variable times don't contain model time, try another file
        if (found_rec == 0) then

          ! If the variable was found but the record not then time not bounded
!          if (mynode==0 .and. found_var==1) then
!            write(*,'(3x,5A)') 'times outside of model time for var ''',
!     &            var_name, ''' in file ',frcfile(var_file_indx)
!          endif

          ! Go to next file if ifile<nfiles
          var_file_indx = var_file_indx+1 ! Go to next file
          irec = 0 ! Reset irec for search in new file. See label (B)

        endif

        ! Close netcdf file - necessary otherwise roms crashes if many variables.
        ierr=nf90_close (ncid)

      enddo ! while ifile<=nfiles.

      ! End of variable & record search
      ! -------------------------------

      if (first) then ! For first time step only

        ! For the very first timestep we need the record just before
        ! vtime(1)<time_model. Therefore irec = irec_old.
        ! E.g. if record 3 gives  vtime(1)>time_model
        ! then record 2 will give vtime(1)<time_model

        ! (A) Handling first data read if model time is between 2 files.
        ! However, if the model time is between two input files for the
        ! very first timestep. The record 1 of the later file will give
        ! vtime(1)>time_model, but we can't simply use irec = irec_old,
        ! as record 0 doesn't exist.
        ! Instead, we use the last record for the previous file for that
        ! variable. This does assume there are no missing periods between
        ! the data, and the data is chronological.

        if(cycle_length==0.) then ! Not cyclical data
          if (irec==1) then ! Model time between input files.
            ! set file index to previous file so that read in time is
            ! from the last record of the previous file containing
            ! this variable.
            var_file_indx = file_indx_old
            ! Explain we are using previous file's last record
            if (mynode==0 .and. var_file_indx>0) then ! var_file_indx>0 so not 1st file
              write(*,'(3x,3A)') 'Model time between input files: ',
     &        'using last entry in previous file ',frcfile(var_file_indx)
            endif
          endif
          irec = irec_old

        else ! Cyclical data
          if (irec==1) then ! Model time between input files.
            irec=var_max_recs ! Use last entry of cyclical file t_data(end)<model_time<t_data(1)
          else
            irec=irec_old ! Use previous record where t(irec_old) < t < t(irec)
          endif
        endif

      endif

      if (found_var_log==0) then
        write(*,'(/1x,4A,I3/)') 'ERROR: find_new_record:: ',
     &              'Could not find var: ', var_name, ' mynode=', mynode
        error stop 'find_new_record: var not found in forcing files'
      endif
      if (found_rec==0) then
        if(mynode==0) write(*,'(/1x,2A/)')
     &  'ERROR: find_new_record: Ran out of time records for ', var_name
        error stop 'ran out of time records'
      endif

      ! For the first timestep we set irec=irec-1, if the first file's
      ! first record (irec=1) is past current time, then irec=1-1=0
      if (irec==0) then
        error stop 'First available forcing record is past current time'
      endif


      end subroutine find_new_record
!-----------------------------------------------------------------------

      subroutine cycle_length_handling(vtime, ncid, var_tid,
     &                                 var_time_name, irec )
      ! Find correct timestep for cyclical forcing files
      ! e.g. 1 year repeating data
      ! NOTE: assumes model_time is mid-point of time-step: time+0.5*dt

      implicit none

      ! inputs/outputs
      real, dimension(1) :: vtime ! Cyclical timestep time to convert to model time if applicable
      integer, intent(in) :: ncid, var_tid
      character(len=*), intent(in) :: var_time_name
      integer :: irec ! Variable record entry

      ! local
!      real    :: cycle_length ! cycle_length of forcing file
      real    :: time_model
      integer :: icycle=0
      integer :: var_max_recs ! Total variable records in file
      integer :: ierr
      integer :: dimid  ! variable ID & dimension ID
      real,dimension(1) :: t_max, t_min


      time_model=tdays+0.5*dt*sec2day ! Calculate model time (assumed mid-point!)

      ! Should never go to the next var_file_indx because this method ensures
      ! time is within cycle_length. So if you hit last cycle record,
      ! then for the next record: icycle = icycle_old+1 and hence
      ! still within cycle_length.

      ! Need index of last record (var_max_recs):
      ierr=nf90_inq_dimid(ncid, var_time_name, dimid) ! Get time dimension ID
      ! Find the size of dimension (var_max_recs)
      ierr=nf90_inquire_dimension(ncid, dimid, len=var_max_recs)
      ! Find the maximum time in cyclical data (t_max)
      ierr=nf90_get_var(ncid, var_tid, t_max, (/var_max_recs/), (/1/))
      ! Find the minimum time in cyclical data (t_min)
      ierr=nf90_get_var(ncid, var_tid, t_min, (/1/), (/1/))

      ! Get time to nearest model time using 'icycle' using 2 steps:
      ! 1) icycle=floor( (time_model/cycle_length) ) ! Floor for integer division
      ! 2) vtime=vtime+icycle*cycle_length
      ! However, need to handle irec=1 and irec=var_max_recs

      if(irec==1) then

        ! If irec==1 is nearest for t_data > t_model, it is possible that
        ! model_time can sit in different icycle's.
        ! E.g. cycle_length=360, t_data=[15;45;...;345]
        ! then if t=346, icycle=0, but if t=374, icycle=1.
        ! To avoid this we add the difference: (cycle_length-t_data(end)) to model time.
        ! icycle=floor(( (time_model+cycle_length-t_max(1)) /cycle_length))
        ! icycle=floor((346+360-345)/360)) = 361/360 = 1
        ! icycle=floor((374+360-345)/360)) = 389/360 = 1
        ! This ensures icycle=1 for irec=1, to get correct relative time.

        ! Find number of cycle_length's to get within model time
        ! Use floor() instead of just integer divsion so that always uses
        ! lower number even when times are negative.
        icycle=floor(( (time_model+cycle_length-t_max(1)) /cycle_length)) ! Force lower integer division
        ! Set input data time record as per model time
        vtime=vtime+icycle*cycle_length

      elseif(irec==var_max_recs) then

        ! From the above example, if model startup time=370s then:
        ! 345 < tmodel < 15+360
        ! In this case for the 'first' step in find_new_record, irec is
        ! set to var_max_recs, as that time < model time.
        ! For the icycle calc to still work we must subtract t_min(1):
        ! icycle=floor(( (time_model-t_min(1)) /cycle_length) )
        ! icycle=floor((370-15)/360)) = 355/360 = 0
        ! This will ensure correct 'icycle' even for 'first' switch in find_record

        icycle=floor(( (time_model-t_min(1)) /cycle_length) ) ! Force integer division
        ! Set input data time record as per model time
        vtime=vtime+icycle*cycle_length

      else

        ! Find number of cycle_length's to get within model time
        icycle=floor( (time_model/cycle_length) ) ! Force integer division
        ! Set input data time record as per model time
        vtime=vtime+icycle*cycle_length

      endif

      ! If cycle_length==0 then not cyclical forcing, vtime remains unchanged.

!      if(mynode==0) print*,'vtimeA=',vtime,' time_model=',time_model ! debug

      end subroutine cycle_length_handling
!-----------------------------------------------------------------------

      subroutine read_var_frc( var_file_indx, irec, var_array,
     &    var_name, var_type, var_time, var_time_name, coarse )

      ! ========================================
      ! READ ONE VARIABLE FROM INPUT NETCDF FILE
      ! ========================================

      implicit none

      ! input/output
      ! ------------
      integer, intent(in) :: var_file_indx
      integer, intent(in) :: irec

      ! param.h needed for GLOBAL_2D_ARRAY
!#include "param.h"
      real, intent(in) :: var_array(GLOBAL_2D_ARRAY,1) ! variables array to record

      ! Similar to vname(1,inx) from old ncvars.h
      character(len=*), intent(in) :: var_name ! Variable short name
      integer, intent(in)          :: var_type ! variable type: u-, v- or rho-point

      ! Set var_time for read in variables
      real,dimension(1)   :: var_time
      character(len=*)    :: var_time_name

      integer, intent(in) :: coarse ! Whether function is interpolated from coarse data

      ! local
      ! -----
      integer ierr,ncid
      integer :: var_tid ! Variable time_id

      ! open the file (maybe check if the file is open already)
      ! =============
      ierr=nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)
      ! Need error handling

      ! read the time
      ! =============
      ierr = nf90_inq_varid(ncid, var_time_name, var_tid) ! Get time ID (varid)
      ierr = nf90_get_var(ncid, var_tid, var_time, (/irec/), (/1/)) ! Get var_time

      ! Fix time for cyclical forcing file if applicable
      if(cycle_length/=0.) then
        call cycle_length_handling(var_time, ncid, var_tid,
     &                             var_time_name, irec )
      endif
!      if(mynode==0) then ! debug
!      print*,'read_var_frc: var_time=',var_time,' time_model=',tdays+0.5*dt*sec2day ! debug
!      endif

      ! read the var
      ! ============
      call nc_read_var(ncid,var_array,1, ! 1 is nmax (vertical layers)
     &                       var_name,var_type,irec,ierr, coarse)

      ! Close netcdf file (necessary otherwise roms crashes if many variables.)
      ! =================
      ierr=nf90_close (ncid)

      ! NOTE about masking here. C$OMP beware

      ! Error handling
      ! --------------

!  2   if (ierr .ne. nf90_noerr) then
!        write(*,1) var_name, record MYID
!        goto 99                                         !--> ERROR
!      endif
!
!      ! text format for '1' in write(*,1) above
!  1   format(/1x, '### ERROR: nc_read_var :: Cannot read variable ''',
!     &              A, ''' from file, rec =', i6, 3x,A,i4)
!      goto 100 ! Skip 99
!  99  if (may_day_flag == 0) may_day_flag=3
! 100  continue



      end subroutine read_var_frc
! ----------------------------------------------------------------------

      subroutine set_small_arrays (
     &                var_name,      var_time_name, ! Text names
     &                var_data,      var,           ! Variable arrays
     &                var_times,                    ! Input times
     &                dim1,          dim2,          ! small array dimensions (dim2=1 if 1D)
     &                var_file_indx, var_irec,      ! File indx & Current input record
     &                it1, it2 )                    ! Time index placeholders

      ! Wrapper to set small, odd sized arrays 1D or 2D (but not 2D of total grid)
      ! e.g. river/pipe sources
      ! Input data - time must be in days

      implicit none

      ! input/outputs
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: var_time_name
      integer, intent(in) :: dim1, dim2 ! Dimensions of variable (dim2=1 if 1D)

      real :: var_data(dim1,dim2,2) ! Read in raw data for 2 times to interpolate
      real :: var(dim1,dim2)
      real, dimension(2) :: var_times ! Store input file record times
      integer :: var_file_indx, var_irec ! File indx & Current input record
      integer :: it1, it2 ! Placeholders for 2 read in times

      ! local
      integer :: tmp, i, j
      real cff1,cff2 ! for time interpolations
      real tmid_days

      ! Interpolated time is based on 1/2 step time
      tmid_days = tdays + 0.5*dt*sec2day ! Input data times in days

      ! UPDATE LATEST READ IN VARIABLE?
      ! Need var_times(it1) < time < var_times(it2), thus:
      if (var_times(it2) < tmid_days) then  ! We need to refresh data

        if (FIRST_TIME_STEP) then ! Only first timestep

C$OMP MASTER
          call find_new_record( var_name, var_time_name, tmid_days,
     &                          var_file_indx, var_irec )

          call read_small_arrays( var_file_indx, var_irec,
     &                       var_name, var_times(it1), var_time_name,
     &                       dim1, dim2,var_data(:,:,it1) )
C$OMP END MASTER
C$OMP BARRIER

          if (mynode == 0) then
            write(*,'(4x,2A,5x,A,G14.6,1x,A,I4)')
     &               'set_small_array :: read var: ', var_name,
     &               'var_input_time (days) =',var_times(it1),
     &               'rec =',var_irec MYID
          endif

        ! ALL OTHER TIMESTEPS
        ! Every time step switch indices such that
        ! var_data(i,j,1) & var_times(1) -> earlier && var_data(i,j,2) & var_times(2) -> later: step t=x
        ! var_data(i,j,1) & var_times(1) -> later   && var_data(i,j,2) & var_times(2) -> earlier: step t=x+1
        ! and so on. The benifit is only need to change one variable each time.
        else
          tmp = it1
          it1 = it2 ! If it1 = 1, it now equals 2, and vice-versa
          it2 = tmp ! If it2 = 2, it now equals 1, and vice-versa
        endif ! End of first step

C$OMP MASTER
        call find_new_record( var_name, var_time_name, tmid_days,
     &                        var_file_indx, var_irec )

        call read_small_arrays( var_file_indx, var_irec,
     &                     var_name, var_times(it2), var_time_name,
     &                     dim1, dim2, var_data(:,:,it2) )
C$OMP END MASTER
C$OMP BARRIER

        if (mynode == 0) then
          write(*,'(4x,2A,5x,A,G14.6,1x,A,I4)')
     &               'set_small_array :: read var: ', var_name,
     &               'var_input_time (days) =',var_times(it2),
     &               'rec =', var_irec MYID
        endif

      endif ! -> END OF UPDATING VARIABLE INPUT DATA

      ! Temporal interpolation
      ! ======================

      ! Set coefficients for interpolation.
      cff1=( var_times(it2)-tmid_days ) / ( var_times(it2)-var_times(it1) )
      cff2=( tmid_days-var_times(it1) ) / ( var_times(it2)-var_times(it1) )

      ! Check the interpolation coefficients such that modeltime is
      ! bounded by read in variable times.
      if (cff1.ge.0. .and. cff2.ge.0.) then
        do j=1,dim2
          do i=1,dim1
            var(i,j)=cff1*var_data(i,j,it1)+cff2*var_data(i,j,it2)
          enddo
        enddo
      ! Throw warning if interpolation times do not bound model time.
      else
        write(*,'(/1x,4A/3(1x,A,F16.10)/)') 'ERROR: set_small_array',
     &  ':: Model time outside bounds of variable (times):',
     &     var_name, var_time_name,
     &                    'start =', var_times(it1),
     &                    'tmid_days =',     tmid_days,
     &                    'end =',   var_times(it2)
        error stop 'set_small_array :: time interpolation error '
      endif

      end subroutine set_small_arrays
! ----------------------------------------------------------------------

      subroutine read_small_arrays( var_file_indx, irec,
     &    var_name, var_time, var_time_name, dim1, dim2, var_array )

      ! Wrapper to read small array 1D or 2D from netcdf input file

      implicit none

      ! input/output
      integer, intent(in) :: var_file_indx
      integer, intent(in) :: irec
      character(len=*), intent(in) :: var_name ! Variable short name
      real, dimension(1) :: var_time
      character(len=*), intent(in) :: var_time_name
      integer, intent(in) :: dim1, dim2 ! Dimensions of variable (dim2=1 if 1D)
      real :: var_array(dim1, dim2,1) ! variables array to record

      ! local
      integer ierr,ncid, start(4), count(4), i, j
      integer :: var_tid, var_id ! Variable time_id
      real    :: buff(dim1*dim2) ! 1D version of variable array

      ! read the time
      ierr=nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)  ! Open file
      ierr=nf90_inq_varid(ncid, var_time_name, var_tid)           ! Get time ID (varid)
      ierr=nf90_get_var(ncid, var_tid, var_time, (/irec/), (/1/)) ! Get var_time

      ! Fix time for cyclical forcing file if applicable
      if(cycle_length/=0.) then
        call cycle_length_handling(var_time, ncid, var_tid,
     &                             var_time_name, irec )
      endif
!      if(mynode==0) then ! debug
!      print*,'read_var_frc: var_time=',var_time,' time_model=',tdays+0.5*dt*sec2day ! debug
!      endif

      ! Setup netcdf variables
      do i=1,4          ! start,count(1:2) correspond to dim1 and dim2
        start(i)=1      ! dimensions, while start(3) is either
        count(i)=1      ! time record for (2D-fields)
      enddo             ! or start(2) is time record for 1D field
      count(1)=dim1     ! Size of dim1 to read
      if (dim2 > 1) then! If 2D variable
        count(2)=dim2   ! Size of dim2 to read
        start(3)=irec
      else              ! If 1D variable
        start(2)=irec
      endif

      ! read variable
      ierr=nf90_inq_varid(ncid, var_name, var_id) ! get variable ID var_id
      ierr=nf90_get_var(ncid, var_id, buff, start,count) ! buff is vector of coarse grid array
      if(ierr/=0) call handle_ierr(ierr,var_name)
      ierr=nf90_close(ncid)

!      if(mynode==0) print *, var_name,' buff= ',buff
      ! Convert buffer array into roms variable
      do j=1,dim2
        do i=1,dim1
          ! buff stored as rows starting i=1, j=1
          ! Hence var(2,1) = buff(1+dim2*1)
          var_array(i,j,1) = buff(i + dim1*(j-1))
          if(mynode==0) then
!            print*,'var_array(i,j)=,i=,j=',var_array(i,j,1),i,j
          endif
        enddo
      enddo

      end subroutine read_small_arrays
! ----------------------------------------------------------------------

      subroutine read_output_root_name(keyword, kwlen)

      ! ==============================
      ! Read from input file root-name
      ! for output netcdf result files
      ! ==============================

      ! This is called from read_inp.F

      ! 1) Save desired prefix name (root name) for
      !    all additional output files

      ! The following needs to be in the roms.in input file to use this
      ! subroutine (file name needs to be 5 spaces from left margin):
      !output_root_name:
      !     rootname (<- insert desired rootname)

      implicit none

      ! Inputs
      character(len=32) :: keyword
      integer kwlen ! Keyword length

      ! Local
      integer :: ierr = nf90_noerr
      integer, parameter :: input=15


      ! Remove keyword from keyword list to know it's accounted for.
      call cancel_kwd (keyword(1:kwlen), ierr)

      ! Read & save filename root
      read(input,'(A)',err=95) output_root_name

      ! Error handling
      ! --------------
      goto 100
      ! Error for read(input,'(A)',err=95)
  95  write(*,'(/1x,4A/)') '### ERROR: read_write :: Cannot read ',
     &                       'entry ''', keyword(1:kwlen), '''.'
      error stop
 100  continue


      end subroutine read_output_root_name


!-----------------------------------------------------------------------

      subroutine ncdf_create_file( fname, ncid, prev_fill_mode,
     &                   rec_per_file, total_rec, auxil,
#ifdef SOLVE3D
     &                   r3dgrd, u3dgrd, v3dgrd, w3dgrd, ! 3D grid dimensions
#endif
     &                   r2dgrd, u2dgrd, v2dgrd )        ! 2D grid dimensions

      ! =============================================
      ! CREATE NEW NETCDF FILE with global attributes
      ! dimensions and time only
      ! ========================

      implicit none

      ! Needed for 'N'
!#include "param.h"

      ! Inputs
      ! ------
      integer, intent(out) :: ncid ! netcdf file ID
      character(len=32), intent(inout) :: fname ! desired netcdf file name
      integer, intent(inout) :: prev_fill_mode  ! Needed for nf90_set_fill
      integer, intent(in)    :: rec_per_file    ! records per file
      integer, intent(in)    :: total_rec       ! Total netcdf records so far for these variables
      ! Associated dimension arrays (e.g. r2dgrd = xi_r, eta_r, time)
      integer, intent(out) :: r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)
#ifdef SOLVE3D
     &                      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif

      ! Local
      ! -----
      integer var_id_tmp, lfnm, lvar, lenstr, lstr
      character(len=64) :: fname_tmp ! Taken from read_inp.F
      integer :: ierr = 0 ! Set local no_error
      integer timedim(1) ! Made it (1). Can remain local as time defined here

      ! External function needed
      integer, external :: my_nf_def_dim

      ! Common
      ! Need ncvars.h for xi_rho, xi_u, eta_rho, eta_v and iaux
#include "ncvars.h"


      ! Add node number to file name (Functionality normaly in read_inp.F)
      ! ----------------------------

      ! Only done once - commented out if read_inp_net_flux is used
      if(total_rec == 0) then

        fname_tmp=fname
        lstr=lenstr(fname_tmp) ! Function in lenstr.F
#if defined MPI && defined PARALLEL_FILES
        ! Insert MPI node numbers to file name
        call insert_node(fname_tmp, lstr, mynode, NNODES, ierr)
#endif
        fname=fname_tmp(1:lstr)

      endif ! (nrechis_flux == 0)


      ! Insert time index to file name
      ! ------------------------------

      lfnm=lenstr(fname) ! Function in lenstr.F
        ! DevinD: a fix to deal with file ncrechis +2 at the beginning
        ! for wrt_his, because rec=rec+1 in different order for new ncdf stuff
      if (total_rec==0) then
        lvar=total_rec-(1+mod(total_rec-1, rec_per_file))
      else
        lvar=total_rec
      endif
      call insert_time_index(fname, lfnm, lvar,  ierr)
      if (ierr /= 0) call handle_ierr(ierr) ! Error termination

      ! Create new file
      ! ---------------

      ierr=nf90_create(fname, nf90_clobber+nf90_64bit_offset, ncid)

      if (ierr == nf90_noerr) then

        ! Set fill value - nf90_nofill produces optimized writes:
        ! "Use "no fill" mode, omitting the initialization of variable values
        ! with fill values. Creation of large files much faster, but
        ! eliminates possibility of detecting reading of values not yet written"
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        ! Error handling
        if (ierr /= nf90_noerr) then
          write(*,'(1x,4A,I4)') '### ERROR: read_write::Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
        endif

      ! Error: cannot create file
      else

        if(mynode==0) then
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: read_write :: ',
     &       'Cannot create ''', fname, '''.', nf90_strerror(ierr)
        endif
        error stop ! Stop simulation
      endif

      ! Put global attributes:
      !---- ------ -----------

      call put_global_atts(ncid, ierr)

      ! Define dimensions for field placement on staggered grids:
      ! ------ ---------- --- ----- --------- -- --------- -------

      ierr=my_nf_def_dim(ncid, 'xi_rho',  xi_rho,  r2dgrd(1))
      ierr=my_nf_def_dim(ncid, 'xi_u',    xi_u,    u2dgrd(1))
      ierr=my_nf_def_dim(ncid, 'eta_rho', eta_rho, r2dgrd(2))
      ierr=my_nf_def_dim(ncid, 'eta_v',   eta_v,   v2dgrd(2))
#ifdef SOLVE3D
      ierr=my_nf_def_dim(ncid, 's_rho',   N,       r3dgrd(3))
      ierr=my_nf_def_dim(ncid, 's_w',     N+1,     w3dgrd(3))
#endif
      ierr=my_nf_def_dim(ncid, 'time', nf90_unlimited, timedim(1))
      ierr=my_nf_def_dim(ncid, 'auxil',   iaux,     auxil(1))

      ! Transfer prescribed dimensions to rest of dimension arrays

      auxil(2)=timedim(1)

      r2dgrd(3)=timedim(1)                          ! Free surface

      u2dgrd(2)=r2dgrd(2) ; u2dgrd(3)=timedim(1)    ! 2D UBAR-type

      v2dgrd(1)=r2dgrd(1) ; v2dgrd(3)=timedim(1)    ! 2D VBAR-type

#ifdef SOLVE3D
      r3dgrd(1)=r2dgrd(1)
      r3dgrd(2)=r2dgrd(2) ; r3dgrd(4)=timedim(1)    ! 3D RHO-type

      u3dgrd(1)=u2dgrd(1) ; u3dgrd(3)=r3dgrd(3)     ! 3D U-type
      u3dgrd(2)=r2dgrd(2) ; u3dgrd(4)=timedim(1)    !

      v3dgrd(1)=r2dgrd(1) ; v3dgrd(3)=r3dgrd(3)     ! 3D V-type
      v3dgrd(2)=v2dgrd(2) ; v3dgrd(4)=timedim(1)    !

      w3dgrd(1)=r2dgrd(1) ; w3dgrd(4)=timedim(1)    ! 3D W-type
      w3dgrd(2)=r2dgrd(2)
#endif


!#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
! || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)
!
!! Define grid variables.
!! ------ ---- ----------
!
!      if (total_rec <= 1) call def_grid(ncid, r2dgrd)
!#endif

      ! Insert time variables
      ! ---------------------

      ! Time-step number and time-record indices: (history file only, this
      ! may be needed in the event when a history record is used to restart
      ! the current model run);
      ierr=nf90_def_var (ncid, 'time_step', nf90_int,
     &                          auxil, var_id_tmp )

      ! Set variable long name
      ierr=nf90_put_att (ncid, var_id_tmp, 'long_name',
     &             'time step and record numbers from initialization')

      ! Define ocean_time variable
      call nc_define_var(ncid, 'ocean_time',
     &            'Time since initialization', 'second', timedim,ierr)


      end subroutine ncdf_create_file


! ----------------------------------------------------------------------

      subroutine nc_define_var(ncid,var_name,var_long_name,var_units,
     &                              var_grd,ierr)

      ! =======================================
      ! Define variables for output netcdf file
      ! =======================================

      ! Flow of function calls taken from old def_his.F

      implicit none

      ! Input variables
      ! ---------------

      integer, intent(in) :: ncid ! netcdf file ID
      ! Variable attributes - similar to vname() from old roms ncvars.h
      character(len=*), intent(in) :: var_name, var_long_name, var_units
      integer, dimension(:), intent(in)  :: var_grd ! Could be 2d or 3d grid
      integer ierr ! Track netcdf errors

      ! Local variables
      ! ---------------

      integer var_id ! variable ID

      ! Taken from def_his.F
#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real*8, parameter :: spv_set=1.D+33
# else
      real*4, parameter :: spv_set=1.E+33
# endif
#endif

      ! Define variables
      ! ================

      ! Get variable ID
      ierr = nf90_inq_varid (ncid,var_name,var_id)
!      if(ierr .ne. nf90_noerr) goto 2 ! Don't know why this creates an error???

      ! Set variable short name
        ! Note: nf90 functionality no longer requires number of dimensions
        ! assume it can work it out from size of array var_grd.
      ierr=nf90_def_var (ncid, var_name, nf90_float,
     &                         var_grd, var_id)
      if(ierr .ne. nf90_noerr) goto 2
      ! Old code:
      !      ierr=nf_def_var (ncid, var_name(1:lvar), nf_float,
      !     &                       var_dims, var_grd, his_var)

      ! Set variable long name
      ierr=nf90_put_att (ncid, var_id, 'long_name', var_long_name)
      if(ierr .ne. nf90_noerr) goto 2

      ! Set variable units
      ierr=nf90_put_att (ncid, var_id, 'units', var_units)
      if(ierr .ne. nf90_noerr) goto 2

      ! Set variable masking fill value
#ifdef MASK_LAND_DATA
      ierr=nf90_put_att (ncid, var_id, '_FillValue', spv_set)
      if(ierr .ne. nf90_noerr) goto 2
#endif

      ! Error handling
      ! --------------

  2   if (ierr .ne. nf90_noerr) then
        write(*,1) var_name MYID
        goto 99                                         !--> ERROR
      endif

      ! text format for '1' in write(*,1) above
  1   format(/1x, '### ERROR: nc_define_var :: Cannot def. variable ''',
     &              A, ''' into history file, rec =', i6, 3x,A,i4)
      goto 100 ! Skip 99
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue


      end subroutine nc_define_var


! ----------------------------------------------------------------------

      subroutine nc_write_var(ncid,var_array,nmax,var_name,
     &                             var_type,record,ierr)

      ! =====================================
      ! Write variables to output netcdf file
      ! =====================================

      ! Taken from old wrt_his.F

      implicit none

      ! Input variables
      ! ---------------

      integer, intent(in) :: ncid ! netcdf file ID
      ! nmax - number of vertical indices (=1 if 2D variable)
      ! this process needed so ncdf_write can handle both
      ! 2D or 3D variable arrays
      integer, intent(in) :: nmax
      ! param.h & cppdefs.h needed for var_array for GLOBAL_2D_ARRAY
!#include "param.h"
      real, intent(in)    :: var_array(GLOBAL_2D_ARRAY,nmax) ! variables array to record
      integer, intent(in) :: record ! Timestep to record
      ! Similar to vname(1,inx) from old ncvars.h
      character(len=*), intent(in) :: var_name ! Variable short name
      integer, intent(in) :: var_type ! variable type: u-, v- or rho-point
      integer ierr ! Track netcdf errors

      ! Local variables
      ! ---------------
      integer var_id ! variable ID

      ! Subroutine
      ! ==========

      ! Get variable ID
      ierr = nf90_inq_varid (ncid,var_name,var_id)
      if(ierr .ne. nf90_noerr) goto 2

      ! Write variable
      ierr = ncdf_write_mod(ncid, var_id, record, var_type, var_array, nmax)

      ! Error handling
      ! --------------

  2   if (ierr .ne. nf90_noerr) then
        write(*,1) var_name, record MYID
        goto 99                                         !--> ERROR
      endif

      ! text format for '1' in write(*,1) above
  1   format(/1x, '### ERROR: nc_write_var :: Cannot write variable ''',
     &              A, ''' into history file, rec =', i6, 3x,A,i4)
      goto 100 ! Skip 99
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue


      end subroutine nc_write_var


! ----------------------------------------------------------------------

      subroutine nc_write_time( ncid, record, total_recs )

      ! =============================
      ! Write timestep and ocean time
      ! to output netcdf file
      ! =====================

      implicit none

      ! Inputs
      ! ------
      integer :: ncid
      integer :: record     ! Current file record number
      integer :: total_recs ! Total records for variable

      ! Local
      ! -----
      integer ibuff(6), start(2), count(2) ! iaux = 6 from wrt_his.F
      integer :: var_id_tmp, ierr


      ibuff(1)=iic-1 ; ibuff(2)=999999     ! Time step and nrecrst=999999 DevinD hardcode
      ibuff(4:6)=0   ; ibuff(3)=total_recs ! record numbers. iaux = 6 in ncvars.h
!#ifdef AVERAGES
!      ibuff(4)=nrecavg
!#endif
      start(1)=1      ; count(1)=6 ! iaux = 6 in ncvars.h
      start(2)=record ; count(2)=1

      ! get time step var_id & record time step info
      ! --------------------------------------------
      ierr=nf90_inq_varid(ncid, 'time_step', var_id_tmp)
      ierr=nf90_put_var(ncid, var_id_tmp, ibuff, start, count)
      if (ierr /= nf90_noerr) then
        if (mynode==0) then
          write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') 'ERROR: read_write :: ',
     &       'Cannot write variable ''time_step'' into history file, ',
     &       'rec =', record, nf90_strerror(ierr) MYID
        endif
        error stop
      endif


      ! get ocean_time var_id & record ocean_time
      ! -----------------------------------------

      ! Which one to write into the file, time or tdays, is decided by attribute
      ! "units" stored as vname(3,indxTime).  When computng time is always
      ! in seconds, however it is more convenient to have it in days in all
      ! the files.

!        if (vname(3,indxTime)(1:6) == 'second') then
      ierr=nf90_inq_varid(ncid, 'ocean_time', var_id_tmp)
      ierr=nf90_put_var(ncid, var_id_tmp, time,(/record/))
!        else
!          ierr=nf_put_var1_double(nchis, hisTime, record, tdays)
!        endif

      if (ierr /= nf90_noerr) then
        write(*,'(/1x,4A,i6/11x,A,3x,A,i4/)') 'ERROR: read_write :: ',
     &        'Cannot write variable ''', 'ocean_time',
     &        ''' into history file, rec =', record, nf90_strerror(ierr)
     &            MYID
      endif

      end subroutine nc_write_time


! ----------------------------------------------------------------------

      subroutine nc_read_var(ncid,var_array,nmax,var_name,
     &                             var_type,record,ierr, coarse)

      ! =====================================
      ! Read variables from input netcdf file
      ! =====================================

      ! Similar to old get_forces.F

      ! This is based on ncid being common to all variables.
      ! Would be better if takes in file_id, and then computes ncid
      ! WEC used this because ncid already calculated in read_wec_var

      implicit none

      ! Input variables
      ! ---------------

      integer, intent(in) :: ncid ! netcdf file ID
      ! nmax - number of vertical indices (=1 if 2D variable)
      ! this process needed so ncdf_read can handle both
      ! 2D or 3D variable arrays
      integer, intent(in) :: nmax
      ! param.h & cppdefs.h needed for var_array for GLOBAL_2D_ARRAY
!#include "param.h"
      real, intent(in) :: var_array(GLOBAL_2D_ARRAY,nmax) ! variables array to record
      integer, intent(in) :: record ! Timestep to record
      ! Similar to vname(1,inx) from old ncvars.h
      character(len=*), intent(in) :: var_name ! Variable short name
      integer, intent(in)  :: var_type ! variable type: u-, v- or rho-point
      integer ierr ! Track netcdf errors

      integer, optional :: coarse ! Whether function is interpolated from coarse data

      ! Local variables
      ! ---------------
      integer var_id ! variable ID

      ! Subroutine
      ! ==========

      ! Check if 'coarse' is available as optional arguement, if not set to 0
      if(.not. present(coarse))then
        coarse=0
      endif

      ! Get variable ID
      ierr = nf90_inq_varid (ncid,var_name,var_id)
      if(ierr .ne. nf90_noerr) goto 2

      ! Read variable

      if(coarse==1) then ! INTERPOLATED CASE: read coarse grid data

        ierr = ncdf_read_coarser_grid(ncid, var_id, record,
     &                           var_type, var_array, nmax) ! grid type, array, nmax

      else ! NORMAL CASE: Reading forcing from same sized grid

        ierr = ncdf_read_mod(ncid, var_id, record,
     &                           var_type, var_array, nmax)

      endif

      ! Error handling
      ! --------------

  2   if (ierr .ne. nf90_noerr) then
        write(*,1) var_name, record MYID
        goto 99                                         !--> ERROR
      endif

      ! text format for '1' in write(*,1) above
  1   format(/1x, '### ERROR: nc_read_var :: Cannot read variable ''',
     &              A, ''' from file, rec =', i6, 3x,A,i4)
      goto 100 ! Skip 99
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue


      end subroutine nc_read_var
!-----------------------------------------------------------------------

!      subroutine nc_check_units( var_file_indx, var_name, var_units )
!
!      ! ===========================================
!      ! Ensure units in file are same as roms units
!      ! ===========================================
!
!      ! DevinD - not tested yet might not work. Parked for now.
!
!      implicit none
!
!      ! input/output
!      ! ------------
!      integer, intent(in) :: var_file_indx
!      character(len=*), intent(in) :: var_name ! Variable short name
!      character(len=*), intent(in) :: var_units
!
!      ! local
!      ! -----
!      integer ierr,ncid
!      integer var_id
!      character(len=20) :: file_var_units ! Variable units of forcing file
!
!
!      ! open the file (maybe check if the file is open already)
!      ! =============
!      ierr=nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)
!      call handle_ierr(ierr)
!
!      ! Get units
!      ! =========
!      ierr = nf90_inq_varid(ncid, var_name, var_id) ! Get variable ID (varid)
!      ierr = nf90_get_att(ncid, var_id, 'units', file_var_units)
!
!      if (mynode==0) print *, 'file_var_units',file_var_units
!
!      if (var_name /= file_var_units) then
!        if(mynode==0) then
!          write(*,'(/3x,3A/)') 'ERROR: read_write.F :: ',
!     &       'nc_check_units - units not equal for var: ', var_name
!        endif
!        error stop ! Stop simulation, fatal error
!      end if
!
!
!      end subroutine nc_check_units
! ----------------------------------------------------------------------

      function ncdf_write_mod(ncid, varid, record, horiz_type, A, nmax)

      ! =============================
      ! Write variable to output file (low-level)
      ! =============================

      ! Routine is tailored to account for MPI tile size differences
      ! and boundary nodes.

      ! Routine is an exact copy of old code's ncdf_write function,
      ! which came from ncdf_read_write.F. However, netcdf calls
      ! changed from e.g. nf_def_var to nf90_def_var.

      ! NOTE: ncdf_write was combined with ncdf_read in ncdf_read_write
      ! to ensure consistency, thus any changes here should probably be
      ! made to read equivalent in this module!

      ! ---------------------------------------------

      ! Write a floating point array into an output netCDF file.

      ! Arguments:
      !            A       real array of standard horizontal dimensions
      !                                  which is to be read or written.
      !            ncid    netCDF ID of in the file.
      !            varid   variable ID of that variable in netCDF file.
      !            record  record number.
      !            type    type of the grid (RHO-, U, V, W, PSI etc.)

      implicit none


      ! Output
      ! ------
      integer ncdf_write_mod

      ! Inputs
      ! ------
      integer ncid, varid, record, horiz_type, nmax
!#include "param.h"
      ! param.h & cppdefs.h needed for A for GLOBAL_2D_ARRAY
      real A(GLOBAL_2D_ARRAY,nmax)
      ! commented: CSDISTRIBUTE_RESHAPE  A...

      ! Local
      ! -----

      logical mask_land_data
      integer vid, i,j,k, shft, ierr
      integer datatype, ndims, natts, dimid(8)
      character(len=16) vname
#include "buffer.h"
#ifdef MASK_LAND_DATA
# include "grid.h"
      real*8, parameter :: spv_set=1.D+33
#endif
#include "compute_starts_counts.h"

      if (varid > 0) then          ! Normally netCDF IDs are positive.
        vid=varid                  ! Negative "varid" is used here as
        mask_land_data=.true.      ! flag to signal that land masking
      else                         ! does not need to be applied for
        vid=-varid                 ! this variable (typically this is
        mask_land_data=.false.     ! reserved for grid variables and
      endif                        ! topography).

      ! Write array from the disk.
      ! ===== ===== ==== === =====

      ! Note that expression for "shft" is exactly the same in all five
      ! cases below, while application of land mask is different for the
      ! variables of different grid staggering; also note effectively .or.
      ! rather than .and. logic in setting velocity values to infinity:
      ! velocity components at the boundary (normal to it) are set to 0,
      ! while the ones fully inside (between two land points) to spv.

#ifdef MASK_LAND_DATA
      if (mask_land_data) then
        if (horiz_type == 0) then
          do k=1,nmax
            do j=jmin,jmax
              shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
              do i=imin,imax
                if (rmask(i,j) > 0.5) then

                  buff(i+shft)=A(i,j,k)

                else

                  buff(i+shft)=spv_set

                endif
              enddo
            enddo
          enddo
        elseif (horiz_type == 1) then
          do k=1,nmax
            do j=jmin,jmax
              shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
              do i=imin,imax
                if (rmask(i,j)+rmask(i-1,j) > 0.5) then

                  buff(i+shft)=A(i,j,k)

                else

                  buff(i+shft)=spv_set

                endif
              enddo
            enddo
          enddo
        elseif (horiz_type == 2) then
          do k=1,nmax
            do j=jmin,jmax
              shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
              do i=imin,imax
                if (rmask(i,j)+rmask(i,j-1) > 0.5) then

                  buff(i+shft)=A(i,j,k)

                else

                  buff(i+shft)=spv_set

                endif
              enddo
            enddo
          enddo
        elseif (horiz_type == 3) then
          do k=1,nmax
            do j=jmin,jmax
              shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
              do i=imin,imax
                if ( rmask(i,j)+rmask(i-1,j)+rmask(i,j-1)
     &                           +rmask(i-1,j-1) > 0.5 ) then

                  buff(i+shft)=A(i,j,k)

                else

                  buff(i+shft)=spv_set

                endif
              enddo
            enddo
          enddo
        endif  !<-- horiz_type == 0,1,2,3
      else  !<-- mask_land_data
#endif
        do k=1,nmax
          do j=jmin,jmax
            shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
            do i=imin,imax

              buff(i+shft)=A(i,j,k)

            enddo
          enddo
        enddo
#ifdef MASK_LAND_DATA
      endif  !<-- mask_land_data
#endif

      ! Put variable in netcdf file
      ! ---------------------------

      ierr=nf90_put_var(ncid, vid, buff, start, count)

      ! Error handling
      ! --------------

      if (ierr /= nf90_noerr) then
        write(*,'(/1x,2A,3x,A,I4/)')  '### ERROR: ncdf_write_mod :: ',
     &             nf90_strerror(ierr) MYID
        write(*,'(12x,A,I7,3x,A,I7/12x,A,I3,7I6)') 'ncid =', ncid,
     &   'varid =', vid, 'start,count =', (start(i),count(i), i=1,4)

        ! get netcdf variable information for error message
        i=nf90_inquire_variable(ncid, vid,vname, datatype,
     &                                  ndims,dimid,natts)
        if (i == nf90_noerr) then
          write(*,'(1x,2A,1x,A,I2,2x,A,I3,2x,A,8I3)') 'vname = ',
     &             vname, 'datatype =', datatype, 'ndims =', ndims,
     &                             'dimid =', (dimid(i), i=1,ndims)

          ! Cycle through bad variable's dimensions to see which
          ! has the error
          do i=1,ndims
            ! Return dimension length (k) for error message
            j=nf90_inquire_dimension(ncid, dimid(i), vname, k)
            if (j == nf90_noerr) write(*,'(29x,2A,I5)') vname,' =',k
          enddo

        endif
      endif

      ncdf_write_mod=ierr


      end function ncdf_write_mod


! ----------------------------------------------------------------------

      function ncdf_read_mod(ncid, varid, record, horiz_type, A, nmax)

      ! =============================
      ! Read variable from input file (low-level)
      ! =============================

      ! Routine is tailored to account for MPI tile size differences
      ! and boundary nodes.

      ! Routine is an exact copy of old code's ncdf_read function,
      ! which came from ncdf_read_write.F. However, netcdf calls
      ! changed from e.g. nf_def_var to nf90_def_var.

      ! NOTE: ncdf_write was combined with ncdf_write in ncdf_read_write.F
      ! to ensure consistency, thus any changes here should probably be
      ! made to write equivalent in this module!

      ! ---------------------------------------------

      ! Read a floating point array from an input netCDF file.

      ! Arguments:
      !            A       real array of standard horizontal dimensions
      !                                  which is to be read or written.
      !            ncid    netCDF ID of in the file.
      !            varid   variable ID of that variable in netCDF file.
      !            record  record number.
      !            type    type of the grid (RHO-, U, V, W, PSI etc.)

      implicit none


      ! Output
      ! ------
      integer ncdf_read_mod

      ! Inputs
      ! ------
      integer ncid, varid, record, horiz_type, nmax
!#include "param.h"
      ! param.h & cppdefs.h needed for A for GLOBAL_2D_ARRAY
      real A(GLOBAL_2D_ARRAY,nmax)
      ! commented: CSDISTRIBUTE_RESHAPE  A...

      ! Local
      ! -----

      logical mask_land_data
      integer vid, i,j,k, shft, ierr
#include "buffer.h"
#ifdef MASK_LAND_DATA
# include "grid.h"
      real*8, parameter :: spv_set=1.D+33
#endif
#include "compute_starts_counts.h"

      if (varid > 0) then          ! Normally netCDF IDs are positive.
        vid=varid                  ! Negative "varid" is used here as
        mask_land_data=.true.      ! flag to signal that land masking
      else                         ! does not need to be applied for
        vid=-varid                 ! this variable (typically this is
        mask_land_data=.false.     ! reserved for grid variables and
      endif                        ! topography).

      ! Read array from the disk.
      !===== ===== ==== === =====

      ierr=nf90_get_var(ncid, vid, buff, start,count)
      if (ierr /= nf90_noerr) then
        write(*,'(/1x,2A,3x,A,I4/)') '### ERROR: ncdf_read_mod :: ',
     &             nf90_strerror(ierr) MYID
      else

      ! Note that expression for "shft" is exactly the same in all five
      ! cases below, while application of land mask is different for the
      ! variables of different grid staggering; also note effectively .or.
      ! rather than .and. logic in setting velocity values to infinity:
      ! velocity components at the boundary (normal to it) are set to 0,
      ! while the ones fully inside (between two land points) to spv.

#ifdef MASK_LAND_DATA
        if (mask_land_data) then
          if (horiz_type == 0) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if (rmask(i,j) > 0.5) then

                    A(i,j,k)=buff(i+shft)

                  else

                    A(i,j,k)=0.D0

                  endif
                enddo
              enddo
            enddo
          elseif (horiz_type == 1) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if (rmask(i,j)+rmask(i-1,j) > 0.5) then

                    A(i,j,k)=buff(i+shft)

                  else

                    A(i,j,k)=0.D0

                  endif
                enddo
              enddo
            enddo
          elseif (horiz_type == 2) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if (rmask(i,j)+rmask(i,j-1) > 0.5) then

                    A(i,j,k)=buff(i+shft)

                  else

                    A(i,j,k)=0.D0

                  endif
                enddo
              enddo
            enddo
          elseif (horiz_type == 3) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax
                  if ( rmask(i,j)+rmask(i-1,j)+rmask(i,j-1)
     &                           +rmask(i-1,j-1) > 0.5 ) then

                    A(i,j,k)=buff(i+shft)

                  else

                    A(i,j,k)=0.D0

                  endif
                enddo
              enddo
            enddo
          ! No masking for river locations as on land cells only
          elseif (horiz_type == 100) then
            do k=1,nmax
              do j=jmin,jmax
                shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
                do i=imin,imax

                    A(i,j,k)=buff(i+shft)

                enddo
              enddo
            enddo
          else
            error stop 'ncdf_read_mod: horiz_type not supported'
          endif  !<-- horiz_type == 0,1,2,3,100
        else  !<-- mask_land_data
#endif
          do k=1,nmax
            do j=jmin,jmax
              shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
              do i=imin,imax

                A(i,j,k)=buff(i+shft)

              enddo
            enddo
          enddo
#ifdef MASK_LAND_DATA
        endif  !<-- mask_land_data
#endif
      endif

      ! Exchange periodic and computational margins (reader only).

#ifdef EXCHANGE
# ifdef MPI
#  define EXCH_ARR_RANGE iwest,ieast,jsouth,jnorth
# else
#  define EXCH_ARR_RANGE 1,Lm,1,Mm
# endif
# ifdef SOLVE3D
      call exchange_tile(EXCH_ARR_RANGE, A,nmax)
# else
      call exchange2d_tile(EXCH_ARR_RANGE, A)
# endif
#endif

      ncdf_read_mod=ierr


      end function ncdf_read_mod


! ----------------------------------------------------------------------
      subroutine handle_ierr(ierr, err_msg)
      ! Handle ierr with error message. Assumed fatal error.
      ! err_msg is optional message to include to know more about
      ! location of error

      implicit none

      ! Inputs
      integer ierr
      character(len=*), optional :: err_msg
      ! Local
      character(len=100) :: f_msg

      ! Check if error message is available
      if(present(err_msg))then
          f_msg=err_msg
      else
          f_msg=''
      endif

      if(mynode==0) then
        write(*,'(/3x,2A/5x,A/12x,A/)') 'ERROR: read_write :: ',
     &       'netcdf ierr != 0 ', trim(f_msg), nf90_strerror(ierr)
      endif

      error stop ! Stop simulation assumed fatal error

      end subroutine handle_ierr


! ----------------------------------------------------------------------

      function ncdf_read_coarser_grid(ncid, varid, record,
     &                                horiz_type, A, nmax)

      ! =========================
      ! Read coarse grid variable
      ! from input file
      ! =========================

      ! Routine is tailored to account for MPI tile size differences
      ! and boundary nodes.
      ! Routine is tailored from ncdf_read_mod.

      ! Initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)

      ! LIMITATIONS
      ! -----------
      ! It can read in forcing data from files with a grid exactly half the
      ! number of grid points as the simulation grid points. It has the
      ! following limitations (note the 'c' is used to signify the coarser
      ! grid):
      ! A) LLm = 2 x LLmc    and    MMm = 2 x MMmc
      !    Note: for the forcing file xi_rho = LLm + 2 (boundary node on each side)
      !    so forcing files won't appear exactly half the size, but LLm=2xLLmc!
      ! B) Boundary MPI tiles need to be used completely (off_xi=0 &
      !    off_eta =0). I.e. the number of grid points in a direction
      !    over the domain are perfectly divisble by the
      !    number of MPI subdomains in that direction. This is necessary
      !    so that the MPI subdomains cover the same physical grid regions
      !    for the refined and coarse grids.
      !    (If LLm/NP_XI divides perfectly, but LLmc/NP_XI does not, then
      !     the width of the subdomains in the xi direction will be a different
      !     size, which the current algorithm cannot handle as interpolation
      !     coefficients are hard-coded and constant.)
      ! C) Input variables & resulting interpolated variable 'A' are both
      !    rho-point variables.
      ! D) Only for 2D arrays.

      ! ----------------------------------------------------------------

      ! Read a floating point array from input netCDF file with coarse grid.

      ! Arguments:
      !            A       real array of standard horizontal dimensions
      !                    which is to be populated after reading in coarse
      !                    grid data to 'buff' variable.
      !            ncid    netCDF ID of in the file.
      !            varid   variable ID of that variable in netCDF file.
      !            record  record number.
      !            type    type of the grid (RHO-, U, V, W, PSI etc.)

      ! Interpolation process
      ! =====================

      ! Interpolation scheme: bilinear interpolation.

      ! Tile relations
      ! --------------

      ! Relation between refined (grid) points to coarse points:
      ! Refined point is closer to tile border, hence cannot interpolate
      ! refined border points until after an MPI exchange.
      ! For interior points,we find between coarse points AA and BB
      ! we interpolate refined points 22,32,23,33 as per diagram below.
      ! (Letters for numbers AA=1,1 & BB=2,2 for coarse points)
      ! (numbers for refined points)
      !
      !       (Internal tile)
      !     ___________________
      !    |                   |
      !    | 14    24 34    44 |
      !    |    AB       BB    |
      !    | 13    23 33    43 |
      !    | 12    22 32    42 |
      !    |    AA       BA    |
      !    | 11    21 31    41 |
      !    |___________________|

      ! Part 1 - Interior points
      ! ------------------------
      ! 1a) Interpolate interior points (excluding refined grid tile
      !     boundary).In this case only 22,32,23,33 are surrounded by
      !     coarse data and hence can be interpolated fully.
      ! 1b) Adjust refined grid indices to interpolate tile boundaries
      !     for tiles that sit on the domain boundary. [See 2d) below]

      ! Part 2 - Tile boundary
      ! ----------------------
      ! 2a) Interpolate north and south boundary in x-direction (exclude corner points)
      !     E.g.: x-interpolation of 24 and 34 using coarse nodes AB & BB
      ! 2b) Interpolate west and east boundary in y-direction (exclude corner points)
      !     E.g.: y-interpolation of 42 and 43 using coarse nodes BA & BB
      ! 2c) Populate the corner nodes with raw value of coarse corner points
      !     (for MPI exchange), since cannot interpolate from 1 point.
      !     E.g.: 11=AA, 41=BA, 14=AB, 44=BB
      ! 2d) Deal with boundary edges of boundary tiles since the domain
      !     has extra rho-point beyond each grid boundary, it is possible
      !     to immediately interpolate these boundaries.
      !     In the 'boundary tile' diagram below we can see that for the
      !     western edge refined nodes 00-03 and 10-13 can be directly
      !     interpolated from coarse nodes ZZ-ZB and AZ-AB.
      !     Similar applies for the southern edge.
      !     Note, a 1D interpolation in x-direction for refined nodes
      !     04 & 14 is needed from coarse nodes ZB & AB. These can only
      !     be completed after an MPI exchange.
      !     Same applies for refined nodes 40 & 41 in y-direction.
      !
      !          (Boundary tile)
      !       |___________________|_
      !       |                   |
      !    04 | 14    24 34    44 |
      ! ZB    |    AB       BB    |
      !    03 | 13    23 33    43 |
      !    02 | 12    22 32    42 |
      ! ZA    |    AA       BA    |
      !    01 | 11    21 31    41 |
      !       |___________________|___
      !
      !    00   10    20 30    40
      ! ZZ         AZ       BZ

      ! Part 3 - MPI exchange
      ! ---------------------
      ! In order to complete the interpolation of the tile boundary, an
      ! MPI exchange is needed to get the coarse point values in the
      ! computational margin beyond the tile boundary from neighbouring tiles.

      ! Part 4 - Interpolate boundary
      ! -----------------------------
      ! 4a) Complete outstanding 1D tile boundary interpolations from
      !     part 2a) & 2b) in remaining perpendicular direction using
      !     values received from MPI exchange.
      ! 4b) Complete tile corner interpolations of raw corner value
      !     now available from MPI exchange.

      ! Note: There is no need to do another MPI exchange as computational
      ! margins of tile are complete and interpolated.

      ! After the MPI exchange we find the following, let's look at the
      ! computational margin of an interior tile at its SW corner:
      !
      ! Current tile computational margins filled by MPI exchange
      !
      !  .       |
      !  .  from |     current
      !  .    W  |      tile
      !  .       |
      !  .       |
      !  . . . . |_______________           ! Solid line is tile boundary
      !  .       .
      !  .  from .     from                 ! 2 point computational margin
      !  .   SW  .       S
      !  . . . . . . . . . . . . .          ! Dotted lines to show tile margin
      !
      !
      !  (Remaining interpolations)
      !
      !  . o   x | x   o   o   o
      !  .       |
      !  . o   x | x   o   o   o
      !  .       |
      !  . y   c | c   y   y   y
      !  . . . . |_______________           ! Solid line is tile boundary
      !  . y   c . c   y   y   y
      !  .       .                          ! 2 point computational margin
      !  . o   x . x   o   o   0
      !  . . . . . . . . . . . . .          ! Dotted lines to show tile margin
      !
      !
      ! After the MPI exchange:
      ! - points marked 'o' are complete already from 'part 1'
      ! - points marked 'x' have had 1D interpolation in y already in 'part 2b'
      !   and all that remains is remaining 1D interpolation in x to complete bilinear interp.
      ! - points marked 'y' have had 1D interpolation in x already in 'part 2a'
      !   and all that remains is remaining 1D interpolation in y to complete bilinear interp.
      ! - points marked 'c' are raw corner values for coarse grid tiles,
      !   and still need full bilinear interpolation.
      !
      ! For boundary tiles:
      !   We only have the 'c' corner exchange for interior corners of the
      !   boundary tiles. Corners on the domain boundary have already
      !   been calculated in part 1.

      ! *Part 1b (adjustment of refined indices)
      ! -------- (Additional description)
      ! For internal MPI tiles, it is not possible to calculate the single
      ! row/column of points directly against the tile border, since the
      ! border coarse points sit closer to the tile centre than the border
      ! refined points (as seen above in 'internal tile' diagram).
      ! This means interpolating the refined points is
      ! not immediately possible without and MPI exchange to get coarse
      ! point values in the tile's computational margin. As per refined
      ! points 11-41 in 'Internal tile' diagram above.
      ! However, since there is an extra rho-point beyond the domain
      ! boundary for all tiles on the domain boundary (see 'Boundary tile'
      ! diagram points 00-30), it is possible to interpolate all the
      ! refined points on the domain boundary immediately.
      ! The most efficient way to do the domain boundary is to use the same
      ! routine to populate the interior points, but just use a mechanism to
      ! shift the refined point indices to correspond to the domain
      ! boundary tiles.
      ! Here is an example of how coarse point indices ic/jc relate to
      ! refined grid indices for interior/boundary tiles.
      ! For interior tiles the loop populates as follows:
      ! ic=1,jc=1 -> ir=2,jr=2 (refined nodes adjacent to the NE direction)
      ! E.g. AA -> 22 in the 'internal tile' diagram
      !
      ! However, for interior and boundary tiles:
      ! TILE       ic jc   ir jr
      ! interior   1  1    2  2
      ! SW corner  0  0    0  0
      ! southern   1  0    2  0
      ! western    0  1    0  2

      ! This are only necessary on the minimum side (west) as the maximum
      ! side (east) is satisfied already by the algorithm.
      ! The goal is thus to get the Cimin & Cjmin to correspond to the
      ! correct ir and jr points using an adjustment index (adj_ir/adj_jr).
      ! It is necessary to look at how ir/jr are calculated relative
      ! to ic/jc in order to get the correct adj_ir/adj_jr.

      ! ----------------------------------------------------------------

      implicit none

      ! Output
      ! ------
      integer ncdf_read_coarser_grid

      ! Inputs
      ! ------
      integer ncid, varid, record, horiz_type, nmax
      ! param.h & cppdefs.h needed for A for GLOBAL_2D_ARRAY
      ! 'A' represents variable used within simulation (refined grid)
      real A(GLOBAL_2D_ARRAY,nmax)

      ! Local
      ! -----
      logical mask_land_data
      integer vid, i,j,k, shft, ierr
#include "buffer.h"
#ifdef MASK_LAND_DATA
# include "grid.h"
      real*8, parameter :: spv_set=1.D+33
#endif

      ! New variables
      integer Cimin,Cimax,Cjmin,Cjmax, Cstart(4),Ccount(4)
      integer ic,jc,ir,jr,Lmc,Mmc,c,adj_ir,adj_jr ! indices of coarse and refined variables
      integer corn_imin,corn_imax,corn_jmin,corn_jmax ! tile corner indices
      real    shftC
      real    temp_int ! Temporary internal rho-point for calculations
      real    temp_ext ! Temporary external rho-point for calculations
      ! Raw coarse corner values and indices from MPI exchange:
      real    SE_raw,SW_raw,NE_raw,NW_raw ! temporary variables
      integer, dimension(4) :: x_SW,y_SW,x_SE,y_SE,x_NW,y_NW,x_NE,y_NE ! temperary variables
      logical corner_bool
      real, parameter :: c9_16 = 9./16. ! interpolation coefficents for efficiency
      real, parameter :: c3_16 = 3./16. ! Decimal required to prevent integer division
      real, parameter :: c1_16 = 1./16.

#include "compute_starts_counts.h"


      ! Compute coarse grid size & indices
      ! ----------------------------------

      ! Currently algorithm is for perfectly
      ! divisible grid points vs mpi tiling only:
      ! For perfect tile size scaling, number of internal nodes in each
      ! direction is half of refined grid for the coarser grid.
      Lmc = Lm/2 ! Lmc = coarse grid number of rho-points in xi-direction of MPI tile
      Mmc = Mm/2

      ! Error handling: if refined grid is not divisble by 2 this routine
      ! cannot be used:
      if(FIRST_TIME_STEP) then ! Only need to check this once for efficiency
        if(modulo(Lm,2) /= 0 .or. modulo(Mm,2) /= 0) then
          write(*,*) 'Error: ncdf_read_coarser_grid: refined subdomain',
     &    ' dimensions (Lm and Mm) not exactly divisible by 2.'
          error stop
        end if
        if(modulo(LLm,NP_XI) /= 0 .or. modulo(MMm,NP_ETA) /= 0) then
          write(*,*) 'Error: ncdf_read_coarser_grid: Domain dimensions',
     &    ' (LLm and MMm) not exactly divisible by MPI tiling (NP_XI ',
     &    'and NP_ETA)'
          error stop
        end if
        if(horiz_type /= 0) then
          write(*,*) 'Error: ncdf_read_coarser_grid: Only supports ',
     &    'rho-point variable interpolation'
          error stop
        end if
        if(nmax /= 1) then
          write(*,*) 'Error: ncdf_read_coarser_grid: Only supports ',
     &    '2D variable interpolation (nmax=1)'
          error stop
        end if
      endif

      ! Coarse grid indices:
      Cimin = imin        ! Interior imin=1 Cimin=1, west bound: imin=0 Cimin=0
      if (EASTERN_MPI_EDGE) then ! Could switch this round (!EAST) for faster code?
        Cimax=Lmc+1       ! east bound: imax=Lm+1 Cimax=Lmc+1
      else
        Cimax=Lmc         ! Interior imax=Lm Cimax=Lmc
      endif
      Cjmin = jmin        ! Interior jmin=1 Cjmin=1, west bound: jmin=0 Cjmin=0
      if (NORTHERN_MPI_EDGE) then
        Cjmax=Mmc+1
      else
        Cjmax=Mmc         ! Interior jmax=Mm Cjmax=Mmc
      endif

      do i=1,4          ! start,count(1:2) correspond to XI- and ETA-
        Cstart(i)=1     ! dimensions, while 3 is either for vertical
        Ccount(i)=1     ! dimension (if any) or for time record
      enddo             ! (2D-fields); 4 is for time record only (3D var)

      Ccount(1)=Cimax-Cimin+1 ! No. points in x. Usually Lmc, E/W boundary = Lmc+1
      Ccount(2)=Cjmax-Cjmin+1 ! No. points in y. Usually Mmc, N/S boundary = Mmc+1
      Cstart(3)=record


      ! Read coarse array from the disk.
      !===== ====== ===== ==== === =====

      ierr=nf90_noerr ! Set ierr to no error
      ierr=nf90_get_var(ncid, varid, buff, Cstart,Ccount) ! buff is vector of coarse grid array
      if (ierr /= nf90_noerr) then ! Error termination
        print *,'In ncdf_read_coarser_grid:'; call handle_ierr(ierr)
      endif
      if(mynode==0) print*, 'Interp-var' ! Confirm interpolating to terminal output

      ! Note: unlike ncdf_read_mod, masking has been removed as should
      !       be done in force calculation

      ! Part 1 (Interpolate interior of tile)
      ! ------

      ! See part 1b) additional description for adj_jr & adj_ir info
      if (SOUTHERN_MPI_EDGE) then
        adj_jr = -1  ! If Cjmin=0, need jr=0, but jr= 2*jc-Cjmin+1= 2*0-0+1= 1
      else           ! hence need adj_jr=-1
        adj_jr = 0
      endif
      if (WESTERN_MPI_EDGE) then
        adj_ir = -1  ! If Cimin=0, need ir=0, but ir= 2*ic-Cimin+1= 2*0-0+1= 1
      else           ! hence need adj_ir=-1
        adj_ir = 0
      endif

      ! Loop through coarse point indices:
      do jc=Cjmin,Cjmax-1
        ! -1 above because 2nd last course point row covers remaining
        ! internal refined points.

        jr=2*jc-Cjmin + 1 + adj_jr
        ! +1 above because refined point is closer to border so
        ! 1ic1jc -> 2ir,2jr are the required neighbours.

        shftC=1-Cimin+Ccount(1)*(jc-Cjmin) ! shft equivalent for coarse variable
        ! ncdf_read_mod version of shft:
        ! shft=1-imin+count(1)*(j-jmin+(k-1)*count(2))
        ! shft=1... because netcdf indices start from 1 not 0.
        ! For rho points imin=jmin=1 for interior nodes.
        ! Ccount(1) is number of nodes in xi-direction for coarse
        ! tile. This is needed because coarse variable is read into
        ! an array (1D vector) called 'buff' which is stored based on rows
        ! in xi-direction.
        ! For coarse point indices in 'buff' for internal tiles we find:
        ! A(2,1) = buff(1+Lmc) | A(2,2) = buff(2+Lmc) | ... | A(Lmc,2) = buff(Lmc*2)
        ! A(1,1) = buff(1)     | A(2,1) = buff(2)     | ... | A(Lmc,1) = buff(Lmc)

        do ic=Cimin,Cimax-1
          ! +1 below because refined node is closer to border so
          ! 1ic1jc -> 2ir,2jr are neighbours.
          ir=2*ic-Cimin + 1 + adj_ir

          ! 1N (north),1E (east) of ic,jc
          ! (see 'internal tile' diagram above for AA->22)
          A(ir  ,jr  ,1)= c9_16*buff(ic+shftC)             ! AA -> 22
     &                  + c3_16*buff(ic+shftC+1)           ! BA -> 22
     &                  + c3_16*buff(ic+shftC+Ccount(1))   ! AB -> 22
     &                  + c1_16*buff(ic+shftC+Ccount(1)+1) ! BB -> 22

          ! 1N,2E of ic,jc
          A(ir+1,jr  ,1)= c3_16*buff(ic+shftC)             ! AA -> 32
     &                  + c9_16*buff(ic+shftC+1)           ! BA -> 32
     &                  + c1_16*buff(ic+shftC+Ccount(1))   ! AB -> 32
     &                  + c3_16*buff(ic+shftC+Ccount(1)+1) ! BB -> 32

          ! 2N,1E of ic,jc
          A(ir  ,jr+1,1)= c3_16*buff(ic+shftC)             ! AA -> 23
     &                  + c1_16*buff(ic+shftC+1)           ! BA -> 23
     &                  + c9_16*buff(ic+shftC+Ccount(1))   ! AB -> 23
     &                  + c3_16*buff(ic+shftC+Ccount(1)+1) ! BB -> 23

          ! 2N,2E of ic,jc
          A(ir+1,jr+1,1)= c1_16*buff(ic+shftC)             ! AA -> 33
     &                  + c3_16*buff(ic+shftC+1)           ! BA -> 33
     &                  + c3_16*buff(ic+shftC+Ccount(1))   ! AB -> 33
     &                  + c9_16*buff(ic+shftC+Ccount(1)+1) ! BB -> 33


          ! Part 2a (Tile boundary - xi-direction)
          ! -------
          ! South tile border - 1D xi-interpolation only on
          ! southern row of refined points.
          ! Must actually be exchanging with an MPI neighbour to the
          ! south, else tile is on southern domain boundary and excluded.
          ! Same goes for east/north/west tile borders to follow.
          ! This also includes domain boundary tiles, e.g. points 04-14.
          if(south_exchng .and. jc==Cjmin) then ! If first course jc index
            A(ir  ,jmin,1)= 0.75*buff(ic+shftC)      ! AA -> 21
     &                    + 0.25*buff(ic+shftC+1)    ! BA -> 21

            A(ir+1,jmin,1)= 0.25*buff(ic+shftC)      ! AA -> 31
     &                    + 0.75*buff(ic+shftC+1)    ! BA -> 31
          endif
          ! North tile border - 1D xi-interpolation only on
          ! northern row of refined points
          if(north_exchng .and. jc==Cjmax-1) then ! If last course jc index
            A(ir  ,jmax,1)= 0.75*buff(ic+shftC+Ccount(1))  ! AB -> 24
     &                    + 0.25*buff(ic+shftC+Ccount(1)+1)! BB -> 24

            A(ir+1,jmax,1)= 0.25*buff(ic+shftC+Ccount(1))  ! AB -> 34
     &                    + 0.75*buff(ic+shftC+Ccount(1)+1)! BB -> 34
          endif

          ! Part 2b (Tile boundary - eta-direction)
          ! -------
          ! West boundary - 1D eta-interpolation only on
          ! western column of course points
          if(west_exchng .and. ic==Cimin) then ! If first course ic index
            A(imin,jr  ,1)= 0.75*buff(ic+shftC)            ! AA -> 12
     &                    + 0.25*buff(ic+shftC+Ccount(1))  ! AB -> 12

            A(imin,jr+1,1)= 0.25*buff(ic+shftC)            ! AA -> 13
     &                    + 0.75*buff(ic+shftC+Ccount(1))  ! AB -> 13
          endif
          ! East boundary - 1D eta-interpolation only on
          ! eastern column of course points
          if(east_exchng .and. ic==Cimax-1) then ! If last course ic index
            A(imax,jr  ,1)= 0.75*buff(ic+shftC+1)          ! BA -> 42
     &                    + 0.25*buff(ic+shftC+Ccount(1)+1)! BB -> 42

            A(imax,jr+1,1)= 0.25*buff(ic+shftC+1)          ! BA -> 43
     &                    + 0.75*buff(ic+shftC+Ccount(1)+1)! BB -> 43
          endif

        enddo ! i
      enddo ! j

      ! Part 2c (Corner nodes)
      ! -------
      ! Store raw value of coarse tile corners into refined tile corners
      ! to transfer in MPI exchange later for interpolation of
      ! corners.
      ! This should only be for interior corners not tile corners
      ! that sit on the domain boundary (use exchng to determine this).
      ! E.g. Tile with NE corner that is considered interior will
      ! exchange both north and east. if north_exchng .and. east_exchang
      ! Note, have not used 'else if' as tile could exchng in all directions.
      if (south_exchng) then
        ! Point 11 -> SW corner
        if (west_exchng) A(imin,jmin,1)=buff(1)
        ! Point 41 -> SE corner
        if (east_exchng) A(imax,jmin,1)=buff(Ccount(1))
      endif
      if (north_exchng) then
      ! Point 41 -> NW corner
        if (west_exchng) A(imin,jmax,1)=buff(1+Ccount(1)*(Cjmax-Cjmin))
      ! Point 44 -> NE corner
        ! need: Cjmax-Cjmin+1 as bottom boundary tiles start from
        ! Cjmin=0 thus need +1 to account for that
        if (east_exchng) A(imax,jmax,1)=buff(Ccount(1)*(Cjmax-Cjmin+1))
      endif


      ! Part 3 (MPI exchange)
      ! ------
#ifdef EXCHANGE
# ifdef MPI
#  define EXCH_ARR_RANGE iwest,ieast,jsouth,jnorth
# else
#  define EXCH_ARR_RANGE 1,Lm,1,Mm
# endif
# ifdef SOLVE3D
      call exchange_tile(EXCH_ARR_RANGE, A,1) ! Swapped nmax for 1 for now
# else
      call exchange2d_tile(EXCH_ARR_RANGE, A)
# endif
#endif


      ! Part 4a) (Complete tile boundary interpolation)
      ! --------
      ! We introduce indices markers corn_imin/corn_jmin/corn_imax/
      ! corn_jmax that act as imin/jmin/imax/jmax for interior corner points
      ! 'c' (see 'remaining interpolations' diagram), in order to prevent
      ! 1D interpolation of 'c' points that still need 2D interp in part 4b).
      ! For domain boundary tiles we only have the corner exchange for
      ! points 'c' for interior corners of the boundary tile. Corners on
      ! the domain boundary have already been calculated in part 1.
      if (SOUTHERN_MPI_EDGE) then ! using if (not SOUTHERN) might be more efficient order.
        ! Arbitrary negative number (-100) beyond loop range to prevent
        ! corner interpolation of 'c' point, since it is not an internal tile corner.
        corn_jmin = -100
      else
        corn_jmin = jmin
      endif
      if (NORTHERN_MPI_EDGE) then
        corn_jmax = -100
      else
        corn_jmax = jmax
      endif
      if (WESTERN_MPI_EDGE) then
        corn_imin = -100
      else
        corn_imin = imin
      endif
      if (EASTERN_MPI_EDGE) then
        corn_imax = -100
      else
        corn_imax = imax
      endif

      ! Points 'x' (1D interpolation in xi-direction)
      ! ---------- (west + east sides of tile)
      ! Loop through entire tile range including margin
      do j=imin-2,jmax+2

        ! Avoid overridding 'c' points. (Note, overhead of if statement)
        if(j/=corn_jmin-1 .and. j/=corn_jmin .and.
     &     j/=corn_jmax   .and. j/=corn_jmax+1) then

          ! Tile's west boundary:
          if (west_exchng) then ! Not for western domain boundary tiles
                                ! as they are accounted for in part 1
            temp_int      = A(imin,j,1) ! Need to store one value temporarily
            A(imin  ,j,1) = 0.75*temp_int + 0.25*A(imin-1,j,1)
            A(imin-1,j,1) = 0.25*temp_int + 0.75*A(imin-1,j,1)
          endif
          ! Tile's east boundary:
          if (east_exchng) then ! Not for eastern domain boundary tiles
            temp_int      = A(imax,j,1)
            A(imax  ,j,1) = 0.75*temp_int + 0.25*A(imax+1,j,1)
            A(imax+1,j,1) = 0.25*temp_int + 0.75*A(imax+1,j,1)
          endif

        endif
      end do

      ! Points 'y' (1D interpolation in eta-direction)
      ! ---------- (north + south sides of tile)
      do i=imin-2,imax+2

        ! Avoid overridding 'c' points.
        if(i/=corn_imin-1 .and. i/=corn_imin .and.
     &     i/=corn_imax   .and. i/=corn_imax+1) then

          ! Tile's south boundary:
          if(south_exchng) then
            temp_int      = A(i,jmin,1)
            A(i,jmin  ,1) = 0.75*temp_int + 0.25*A(i,jmin-1,1)
            A(i,jmin-1,1) = 0.25*temp_int + 0.75*A(i,jmin-1,1)
          endif
          ! Tile's north boundary:
          if(north_exchng) then
            temp_int      = A(i,jmax,1)
            A(i,jmax  ,1) = 0.75*temp_int + 0.25*A(i,jmax+1,1)
            A(i,jmax+1,1) = 0.25*temp_int + 0.75*A(i,jmax+1,1)
          endif

        endif
      end do

      ! 4b) Corner points 'c' (bilinear interpolation)
      ! ---------------------
      ! I)  Can either loop through 4 corners but requires array,
      ! II) or write it out explicitly for each corner,
      ! Chose option I) for less code.
      ! Store x/y indices in arrays to allow for looping
      ! Indices: 1-SW corner of tile, 2-SE corner, 3-NW corner, 4-NE corner
      ! And within each corner there are 4 'c' points again split
      ! by SW/SE/NW/NE.
      x_SW = [ imin-1, imax  , imin-1, imax   ] ! x indices of SW point at each tile corner 'c'
      y_SW = [ jmin-1, jmin-1, jmax  , jmax   ] ! y indices of SW point at each tile corner 'c'
      x_SE = [ imin  , imax+1, imin  , imax+1 ]
      y_SE = [ jmin-1, jmin-1, jmax  , jmax   ]
      x_NW = [ imin-1, imax  , imin-1, imax   ]
      y_NW = [ jmin  , jmin  , jmax+1, jmax+1 ]
      x_NE = [ imin  , imax+1, imin  , imax+1 ]
      y_NE = [ jmin  , jmin  , jmax+1, jmax+1 ]

      ! Loop through 4 corners
      do c=1,4

        ! Only internal corners remaining, hence to be an internal
        ! corner it needs to be exchanging in both directions
        ! i.e. internal SW corner will exchange south and exchange west
        corner_bool = .false. ! Refresh corner bool
        if      (c==1) then
          if (south_exchng .and. west_exchng) then
            corner_bool = .true.
          endif
        else if (c==2) then
          if (south_exchng .and. east_exchng) then
            corner_bool = .true.
          endif
        else if (c==3) then
          if (north_exchng .and. west_exchng) then
            corner_bool = .true.
          endif
        else if (c==4) then
          if (north_exchng .and. east_exchng) then
            corner_bool = .true.
          endif
        endif

        ! If it is an interior tile corner, interpolate corner
        if (corner_bool==.true.) then

        ! store values so not over-written
        SW_raw = A(x_SW(c),y_SW(c),1)
        SE_raw = A(x_SE(c),y_SE(c),1)
        NW_raw = A(x_NW(c),y_NW(c),1)
        NE_raw = A(x_NE(c),y_NE(c),1) ! Could ignore last one for efficiency as never over-written


        A(x_SW(c),y_SW(c),1)= c3_16*NW_raw + c1_16*NE_raw
     &                      + c9_16*SW_raw + c3_16*SE_raw

        A(x_SE(c),y_SE(c),1)= c1_16*NW_raw + c3_16*NE_raw
     &                      + c3_16*SW_raw + c9_16*SE_raw

        A(x_NW(c),y_NW(c),1)= c9_16*NW_raw + c3_16*NE_raw
     &                      + c3_16*SW_raw + c1_16*SE_raw

        A(x_NE(c),y_NE(c),1)= c3_16*NW_raw + c9_16*NE_raw
     &                      + c1_16*SW_raw + c3_16*SE_raw

        end if ! corner_bool
      end do ! c

      ncdf_read_coarser_grid=ierr ! Return function value as ierr


      end function ncdf_read_coarser_grid
! ----------------------------------------------------------------------

! Copy of subroutine in case anyone wants to use it in future to have
! control of netcdf output file in roms.in file.

!      subroutine read_inp_net_flux(keyword, ierr, kwlen)
!
!      ! =============================================
!      ! Read from input file which net-flux variables
!      ! to output to own netcdf results file
!      ! ====================================
!
!      ! This will serve two functions from read_inp.F
!
!      ! 1) Save output file specifications:
!      !    output intervals, max number records per file, etc
!
!      ! 2) Choose which variables to turn on/off for outputting:
!      !    E.g.: sustr (True or False)
!
!      ! 3) Save desired name of output file
!
!      ! The following needs to be in the roms.in input file to use this
!      ! subroutine (file name needs to be 5 spaces from left margin):
!      !net_flux_fields: write_file, rec_rate, recs/file | sustr  / filename
!      !       T         1    4       T
!      !     net_flux_results.nc
!
!      implicit none
!
!      ! Inputs
!      character(len=32) :: keyword
!      integer ierr, kwlen
!      integer, parameter :: input=15
!
!      ! Local
!      character(len=64) :: fname ! Taken from read_inp.F
!      integer lstr, lenstr ! length of string, function to get lstr
!
!
!      ! Just have one keyword to trigger both 1), 2) & 3) functions
!      ! Remove keyword from keyword list to know it's accounted for.
!      call cancel_kwd (keyword(1:kwlen), ierr)
!
!      read(input,*,err=95) ! Read variables below
!
!      ! Part 1) Specifications of output file
!      ! =====================================
!
!     &  write_file, rec_rate, recs_per_file,
!
!      ! Part 2) Select variables to write
!      ! =================================
!
!      ! Store input file output information (yes or no)
!      ! in wrt_wec_vars which is an array of bools (logical)
!     &  wrt_net_flux_vars(indxSUSTR)
!
!      ! Part 3) store filename
!      ! ======================
!
!      read(input,'(A)',err=95) fname
!      lstr=lenstr(fname) ! Function in lenstr.F
!#if defined MPI && defined PARALLEL_FILES
!      ! Insert MPI node numbers to file name
!      call insert_node(fname, lstr, mynode, NNODES, ierr)
!#endif
!      file_name_flux_net=fname(1:lstr)
!
!      ! Write to terminal output in simulation pre-amble text which
!      ! result variables are being stored
!      ! ---------------------------------
!
!      mpi_master_only write(*,'(/1x,A,L1,2x,A,I5,2x,A,I4,2x,2A)')
!     &     'surf_flux: save results net_flux = ', write_file,
!     &     'nwrt =', rec_rate, 'recs/file =', recs_per_file,
!     &     'file = ', file_name_flux_net(1:lstr)
!
!      mpi_master_only write(*,'(/1x,A,3(/8x,A,T16,L1,T20,A))')
!     &               'fields to be saved in 2D Flux history: (T/F)'
!     &     ,'sustr',  wrt_net_flux_vars(indxSUSTR),   sustr_long_name
!
!! Template for any 3D variables in future
!!#  ifdef SOLVE3D
!!      mpi_master_only write(*,'(/1x,A,8(/8x,A,T16,L1,T20,A))')
!!     &               'fields to be saved in WEC 3D history: (T/F)'
!!     &     ,  'UST',    wrt_wec_vars(indxUST),   ust_long_name
!!#  endif
!
!      ! Error handling
!      ! --------------
!
!      goto 100
!
!      ! Text format for error message: read(input,*,err=95) above
!  95  write(*,'(/1x,4A/)') '### ERROR: surf_flux :: Cannot read ',
!     &                       'entry ''', keyword(1:kwlen), '''.'
!      ierr=ierr+1
!
! 100  continue
!
!
!
!      end subroutine read_inp_net_flux


!-----------------------------------------------------------------------


! ----------------------------------------------------------------------

      ! =======================================================
      ! COPY OF ORIGINAL ROUTINE THAT WORKED FOR DEF_HIS.F, TOO
      ! =======================================================

!      subroutine ncdf_create_new_file( fname,ncid, prev_fill_mode,
!     &                   rec_per_file, total_rec,
!     &                   r2dgrd, u2dgrd, v2dgrd,
!#ifdef SOLVE3D
!     &                   r3dgrd, u3dgrd, v3dgrd, w3dgrd,
!#endif
!     &                   timedim, auxil, N, ierr ) ! ,create_err
!
!      ! =============================================
!      ! CREATE NEW NETCDF FILE with global attributes
!      ! and dimensions only
!      ! ===================
!
!      implicit none
!
!      ! Inputs
!      ! ------
!      integer:: ncid ! netcdf file ID
!      character(len=32), intent(in) :: fname ! desired netcdf file name
!      integer prev_fill_mode ! Needed for nf90_set_fill
!      integer rec_per_file ! records per file
!      integer total_rec
!      ! dimension numbers
!      integer r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)
!#ifdef SOLVE3D
!     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
!#endif
!!      integer r2dgrd1,u2dgrd1,v2dgrd2
!!      integer r3dgrd3,w3dgrd3,timedim,auxil1
!      integer timedim
!      integer N ! Vertical number of rho-points
!      integer ierr ! Track netcdf errors
!
!
!      integer create_err ! DevinD created this to return error from def_his.F
!
!      ! External function needed
!      integer, external :: my_nf_def_dim
!
!      ! Common
!      ! Need ncvars.h for xi_rho, xi_u, eta_rho, eta_v and iaux
!#include "ncvars.h"
!
!
!      ! Create new file
!      ! ---------------
!
!      ierr=nf90_create(fname, nf90_clobber+nf90_64bit_offset, ncid)
!
!      if (ierr == nf90_noerr) then
!
!        ! Set fill value handling
!        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)
!
!        ! Error handling
!        if (ierr /= nf90_noerr) then
!!# ifdef AVRH
!!          write(*,'(1x,4A,I4)') '### ERROR: read_write::Cannot ',
!!# else
!          write(*,'(1x,4A,I4)') '### ERROR: read_write::Cannot ',
!!# endif
!     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
!     &         MYID
!        endif
!
!      ! Error: cannot create file
!      else
!
!!#ifdef AVRH
!!        write(*,'(/1x,4A/12x,A/)')  '### ERROR: read_write :: ',
!!#else
!        write(*,'(/1x,4A/12x,A/)')  '### ERROR: read_write :: ',
!!#endif
!     &       'Cannot create ''', fname, '''.', nf90_strerror(ierr)
!
!        ! Set return error flag
!        create_err = ierr
!
!      endif
!
!      if (rec_per_file == 0) total_rec=0
!
!      ! Put global attributes:
!      !---- ------ -----------
!
!      call put_global_atts(ncid, ierr)
!
!      ! Define dimensions for field placement on staggered grids:
!      ! ------ ---------- --- ----- --------- -- --------- -------
!
!      ierr=my_nf_def_dim(ncid, 'xi_rho',  xi_rho,  r2dgrd(1))
!      ierr=my_nf_def_dim(ncid, 'xi_u',    xi_u,    u2dgrd(1))
!      ierr=my_nf_def_dim(ncid, 'eta_rho', eta_rho, r2dgrd(2))
!      ierr=my_nf_def_dim(ncid, 'eta_v',   eta_v,   v2dgrd(2))
!#ifdef SOLVE3D
!      ierr=my_nf_def_dim(ncid, 's_rho',   N,       r3dgrd(3))
!      ierr=my_nf_def_dim(ncid, 's_w',     N+1,     w3dgrd(3))
!#endif
!      ierr=my_nf_def_dim(ncid, 'time', nf90_unlimited, timedim)
!      ierr=my_nf_def_dim(ncid, 'auxil',   iaux,     auxil(1))
!
!      ! Transfer prescribed dimensions to rest of dimension arrays
!
!      auxil(2)=timedim
!
!      r2dgrd(3)=timedim                          ! Free surface
!
!      u2dgrd(2)=r2dgrd(2) ; u2dgrd(3)=timedim    ! 2D UBAR-type
!
!      v2dgrd(1)=r2dgrd(1) ; v2dgrd(3)=timedim    ! 2D VBAR-type
!
!#ifdef SOLVE3D
!      r3dgrd(1)=r2dgrd(1)
!      r3dgrd(2)=r2dgrd(2) ; r3dgrd(4)=timedim    ! 3D RHO-type
!
!      u3dgrd(1)=u2dgrd(1) ; u3dgrd(3)=r3dgrd(3)  ! 3D U-type
!      u3dgrd(2)=r2dgrd(2) ; u3dgrd(4)=timedim    !
!
!      v3dgrd(1)=r2dgrd(1) ; v3dgrd(3)=r3dgrd(3)  ! 3D V-type
!      v3dgrd(2)=v2dgrd(2) ; v3dgrd(4)=timedim    !
!
!      w3dgrd(1)=r2dgrd(1) ; w3dgrd(4)=timedim    ! 3D W-type
!      w3dgrd(2)=r2dgrd(2)
!#endif
!
!
!#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
! || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)
!
!! Define grid variables.
!! ------ ---- ----------
!
!        if (total_rec <= 1) call def_grid(ncid, r2dgrd)
!#endif
!
!
!      end subroutine ncdf_create_new_file


! ----------------------------------------------------------------------


      end module read_write

