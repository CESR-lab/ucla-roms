      module read_write
      ! =================================================
      ! Contains functions and subroutine for interacting
      ! with netcdf file input and output
      ! =================================
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use netcdf
      implicit none

      private

      public find_rec_new_mod

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


      end module read_write

