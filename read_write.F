      module read_write
      ! =================================================
      ! Contains functions and subroutine for interacting
      ! with netcdf file input and output
      ! =================================
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use netcdf
      implicit none

#include "cppdefs.h"

      ! scalars.h - needed for may_day_flag in error handling
#include "scalars.h"

      private

      ! Needed in order to use the ncdf_write & ncdf_read subroutines
      integer, external :: ncdf_write, ncdf_read

      public find_rec_new_mod
      public nc_define_var
      public nc_write_var
      public nc_read_var

      contains

! --------------------------------------------------------------------
      subroutine find_rec_new_mod(time,v_name,vt_name,ifile,irec)

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
!         time       model time
!         var_name   name of variable in file
!         var_time   name of corresponding time variable in file

! in/out: ifile      index of correct file in frcfile
!         irec       record number in file

      implicit none


      ! input/output
      real(kind=8)    ,intent(in)   :: time
      character(len=*),intent(in)   :: v_name,vt_name
      integer         ,intent(inout):: ifile,irec
      ! local
      real,dimension(1):: vtime
      integer     :: found_var,found_rec
      integer     :: ierr,v_id,ncid,vt_id
      integer     :: nfiles,nrecords
      integer     :: ifile_old,irec_old
      integer     :: dimid ! DevinD added from ncio-nc_size()

      integer,dimension(1) :: ir1,ic1
    
      logical :: first

    
      ! When called for the first time (ifile==0) for this variable, this 
      ! routine will return the last record with vtime<time. 
      
      if (ifile==0) then 
        first    =.true.
        ifile    = 1
        irec_old = 0
      endif

      found_var = 0
      found_rec = 0

      nfiles = max_frc !! max_frc is the number of files in the list frcfile

         do while ((found_rec==0).and.(ifile<=nfiles))

           if (mynode==0) then
!            print *, 'ifile, nfiles: ',ifile,nfiles
             print *, 'Searching for ',v_name,' in file ',frcfile(ifile)
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
               print *, 'Variable ',v_name,' found in file ',frcfile(ifile)
             endif
             found_var = 1
             irec = 1

             ! Find the size of a dimension in netcdf file
             ierr=nf90_inq_dimid(ncid, vt_name, dimid) ! from ncio-nc_size
             ierr=nf90_inquire_dimension(ncid, dimid, len=nrecords)
!            if (mynode==0) write(*,*) 'nrecords=',nrecords ! debug

             do while (irec<=nrecords.and.found_rec==0)

               ierr = nf90_inq_varid(ncid, vt_name, vt_id) ! Get time ID (vt_id_
               ! nf90_get_var needs an arry for last 2 terms hence syntax
               ierr = nf90_get_var(ncid, vt_id, vtime, (/irec/), (/1/))

               ! hard code conversion of input time in days to seconds
               vtime = vtime * 86400
               if (vtime(1)>time) then
                 found_rec = 1
               else
                 irec_old = irec
                 ifile_old = ifile
                 irec = irec+1
               endif
             enddo ! while not found and irec<nrecs
           else ! Var not found in file
             if (mynode==0) then
               print *, 'Variable ',v_name,' NOT found in file:',frcfile(ifile)
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
      character(len=64) text

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

      ! add if variable averaged or not to variable long name
#ifdef AVRH
      text='averaged '/ /var_long_name
#else
      text=var_long_name
#endif

      ! Set variable long name
      ierr=nf90_put_att (ncid, var_id, 'long_name', text)
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
#include "param.h"
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
      ierr = ncdf_write(ncid, var_id, record, var_type, var_array, nmax)

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

      subroutine nc_read_var(ncid,var_array,nmax,var_name,
     &                             var_type,record,ierr)

      ! =====================================
      ! Read variables from input netcdf file
      ! =====================================

      ! Similar to old get_forces.F

      implicit none

      ! Input variables
      ! ---------------

      integer, intent(in) :: ncid ! netcdf file ID
      ! nmax - number of vertical indices (=1 if 2D variable)
      ! this process needed so ncdf_read can handle both
      ! 2D or 3D variable arrays
      integer, intent(in) :: nmax
      ! param.h & cppdefs.h needed for var_array for GLOBAL_2D_ARRAY
#include "param.h"
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
      ierr = ncdf_read(ncid, var_id, record, var_type, var_array, nmax)

      ! Error handling
      ! --------------

  2   if (ierr .ne. nf90_noerr) then
        write(*,1) var_name, record MYID
        goto 99                                         !--> ERROR
      endif

      ! text format for '1' in write(*,1) above
  1   format(/1x, '### ERROR: nc_read_var :: Cannot write variable ''',
     &              A, ''' into history file, rec =', i6, 3x,A,i4)
      goto 100 ! Skip 99
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue


      end subroutine nc_read_var


      end module read_write

