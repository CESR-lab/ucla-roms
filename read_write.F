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
      ierr = ncdf_read_mod(ncid, var_id, record,
     &                           var_type, var_array, nmax)

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
#include "param.h"
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
#include "param.h"
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
!
!      subroutine ncdf_create_new_file(fname,ncid,prev_fill_mode,
!     &                          rec_per_file,total_rec,
!     &                   r2dgrd1,u2dgrd1,r2dgrd2,v2dgrd2,
!     &                   r3dgrd3,w3dgrd3,timedim,auxil1,N,
!     &                   ierr,create_err)
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
!      integer r2dgrd1,u2dgrd1,r2dgrd2,v2dgrd2
!      integer r3dgrd3,w3dgrd3,timedim,auxil1
!      integer N ! Vertical number of rho-points
!      integer ierr ! Track netcdf errors
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
!# ifdef AVRH
!          write(*,'(1x,4A,I4)') '### ERROR:def_avg read_write::Cannot ',
!# else
!          write(*,'(1x,4A,I4)') '### ERROR:def_his read_write::Cannot ',
!# endif
!     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
!     &         MYID
!        endif
!
!      ! Error: cannot create file
!      else
!
!#ifdef AVRH
!        write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_avg read_write :: ',
!#else
!        write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_his read_write :: ',
!#endif
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
!      ierr=my_nf_def_dim(ncid, 'xi_rho',  xi_rho,  r2dgrd1)
!      ierr=my_nf_def_dim(ncid, 'xi_u',    xi_u,    u2dgrd1)
!      ierr=my_nf_def_dim(ncid, 'eta_rho', eta_rho, r2dgrd2)
!      ierr=my_nf_def_dim(ncid, 'eta_v',   eta_v,   v2dgrd2)
!#ifdef SOLVE3D
!      ierr=my_nf_def_dim(ncid, 's_rho',   N,       r3dgrd3)
!      ierr=my_nf_def_dim(ncid, 's_w',     N+1,     w3dgrd3)
!#endif
!      ierr=my_nf_def_dim(ncid, 'time', nf90_unlimited, timedim)
!      ierr=my_nf_def_dim(ncid, 'auxil',   iaux,     auxil1)
!
!
!      end subroutine ncdf_create_new_file


! ----------------------------------------------------------------------


      end module read_write

