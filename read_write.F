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

      ! Public subroutines
      public find_rec_new_mod
      public nc_define_var
      public nc_write_var
      public nc_read_var
!      public nc_check_units
      public set_frc_var_tile
      public ncdf_create_file
      public nc_write_time
      public read_output_root_name

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
     &                it1,           it2 )          ! Time index placeholders

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
     &                       var_time_name )

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
     &                     var_time_name )

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
          if (mynode==0) then
            print *,var_file_indx,frcfile(var_file_indx)
            error stop 'FNR: Forcing file not found'
          endif
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

          ! Get & save time dimension ID
          ierr=nf90_inq_dimid(ncid, var_time_name, dimid)
          ! Find the size of dimensions (var_max_recs)
          ! I.e. total number of records per variable, using time dimension
          ierr=nf90_inquire_dimension(ncid, dimid, len=var_max_recs)

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

            ! Find variable ID for time variable
            ierr = nf90_inq_varid(ncid, var_time_name, var_tid) ! Get time ID (vt_id)
            ! nf90_get_var needs an array for last 2 terms hence syntax
            ierr = nf90_get_var(ncid, var_tid, vtime, (/irec/), (/1/))

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

        if (irec==1) then ! Model time between input files.
          ! set file index to previous file so that read in time is
          ! from the last record of the previous file containing
          ! this variable.
          var_file_indx = file_indx_old
          ! Explain we are using previous file's last record
          if (mynode==0 .and. var_file_indx>0) then ! var_file_indx>0 so not 1st file
            write(*,'(3x,3A)') 'Model time between input files - using',
     &           ' last entry in previous file ',frcfile(var_file_indx)
          endif
        endif

        irec = irec_old

      endif

      if (found_var_log==0) then
        write(*,'(/1x,3A/)') 'ERROR: find_new_record:: ',
     &                       'Could not find var:', var_name
        error stop 'find_new_record: var not found in forcing files'
      endif
      if (found_rec==0) then
        mpi_master_only write(*,'(/1x,2A/)')
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

      subroutine read_var_frc( var_file_indx, irec, var_array,
     &               var_name, var_type, var_time,var_time_name )

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
      integer, intent(in)  :: var_type ! variable type: u-, v- or rho-point

      ! Set var_time for read in variables
      real,dimension(1) :: var_time
      character(len=*)  :: var_time_name

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
      ierr = nf90_get_var(ncid, var_tid, var_time, (/irec/), (/1/))

      ! read the var
      ! ============
      call nc_read_var(ncid,var_array,1,
     &                       var_name,var_type,irec,ierr)

!     call nf90_close(ncid)

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
     &                             var_type,record,ierr)

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

      ! Local variables
      ! ---------------
      integer var_id ! variable ID

      ! Subroutine
      ! ==========

      ! Get variable ID
      ierr = nf90_inq_varid (ncid,var_name,var_id)
      if(ierr .ne. nf90_noerr) goto 2

      ! Read variable
      ierr = ncdf_read_mod(ncid, var_id, record,
     &                           var_type, var_array, nmax)

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

      ! ============================
      ! Write variable to input file (low-level)
      ! ============================

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

      ! ===========================
      ! Read variable to input file (low-level)
      ! ===========================

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
        write(*,'(/1x,2A,3x,A,I4/)') '### ERROR: ncdf_read :: ',
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
          endif  !<-- horiz_type == 0,1,2,3
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
      subroutine handle_ierr(ierr)
      ! Handle ierr with error message. Assumed fatal error.

      implicit none

      integer ierr

      if(mynode==0) then
        write(*,'(/3x,2A/12x,A/)') 'ERROR: read_write :: ',
     &       'netcdf ierr != 0 ', nf90_strerror(ierr)
      endif

      error stop ! Stop simulation assumed fatal error

      end subroutine handle_ierr

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

