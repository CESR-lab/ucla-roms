      module river_frc
      ! Contains all the neccesary components to produce the riv_uflx,riv_vflx
      ! arrays which have the the horizontal river volume flux in m2/s at the current time
      ! riv_uflx, riv_vflx should be on faces between a masked and unmasked cell,
      ! in the direction of the unmasked cell.

      ! call set_river_source_tile(istr,iend,jstr,jend, ierr)
  
       use dimensions !! has dimensions, a list of forcing files, and rmask
       implicit none
       ! gnx full domain, nx  subdomain
       private

#include "cppdefs.h"
#include "param.h"
#include "grid.h"
#include "ncvars.h"

      ! These are usable by code that uses this module
       real,public, dimension(nx+1,ny) :: riv_uflx
       real,public, dimension(nx,ny+1) :: riv_vflx

       ! These are private to this module
       character(len=10) :: var_name = 'river_flux'
       character(len=15) :: var_time = 'river_flux_time'
       character(len=15) :: var_time = 'river_flux_time'
       character(len=15) :: var_time = 'river_flux_time'
       character(len=15) :: var_time = 'river_flux_time'

       integer    :: ifile = 0
       integer    :: irec  = 0
  
       integer                  :: it1 = 1
       integer                  :: it2 = 2
       real, dimension(2)       :: frc_tim = [-99,-99]
       real, dimension(nx,ny,3) :: rflx 

       public nx
       public set_river_frc

      contains
       ! --------------------------------------------------------------------
        subroutine ana_river_frc
         ! Analytical river forcing data
        end subroutine ana_river_frc

       ! --------------------------------------------------------------------
        subroutine ini_river_frc
         use ncio
         implicit none
         ! Initialize river forcing
        end subroutine ini_river_frc


       ! --------------------------------------------------------------------
       subroutine read_river_frc(time)
         ! Reads river forcing data from input file
         ! The first time, it reads two instances
         use ncio
         implicit none

         real(kind=8), intent(in) :: time
         integer  :: tmp
         integer  :: is,js,ierr

         !! We need to deal with the fact that currently, the first subdomain
         !! includes the west (or south) buffer, which we'll want to skip.
         !! Find_rec_new checks for this and will change is or js from 1 to 2 if needed


         if (frc_tim(it1) < 0) then !! first time
           ! calling find_rec_new with irec==0 returns the last record with var_time<time
           call init_river_frc
           call find_rec_new(time,var_name,var_time,ifile,irec)
           if (mynode==0) then
!             print *, '----------------------'
!             print *, 'READ_RIVER'
!             print *, ifile,irec,is,js
!             print *, frcfile(ifile)
!             print *, var_name,' ',var_time
!             print *, 'it1: ', it1
!             print *, 'nx, ny: ', nx,ny
!             print *, 'irec: ', irec
              print *, '-------'
           endif
           
           call nc_read(frcfile(ifile),var_time,frc_tim(it1),[irec],[1])
           frc_tim(it1) = frc_tim(it1) + 360  !! place holder for proper cycle_length treatment

           call nc_read(frcfile(ifile),var_name,rflx(:,:,it1),[is,js,irec],[nx,ny,1])
!          rsum = sum(rflx(:,:,it1))
!          if (rsum>0) then
!            print *, 'rsum: ',mynode,rsum
!          endif
!          call mpi_barrier(ocean_grid_comm,ierr)
!          stop 'returning from 2D nc_read'
         else
           tmp =it1
           it1 =it2
           it2 =tmp
         endif

!        if (mynode==0) then
!           print *,'read_river_frc, it1: ', frc_tim(it1)
!        endif

         call find_rec_new(time,var_name,var_time,ifile,irec,is,js)
!        if (mynode==0) then
!          print *,'back from FRN ', ifile,irec
!        endif
         call nc_read(frcfile(ifile),var_time,frc_tim(it2),[irec],[1])
         frc_tim(it2) = frc_tim(it2) + 360  !! place holder for proper cycle_length'
!        if (mynode==0) then
!          print *,'back from nc_read ', frc_tim(it2)
!        endif
!        call mpi_barrier(ocean_grid_comm,ierr)

         call nc_read(frcfile(ifile),var_name,rflx(:,:,it2),[is,js,irec],[nx,ny,1])

!        if (mynode==0) then
!           print *,'read_river_frc: ', frc_tim(it1),frc_tim(it2)
!        endif

       end subroutine read_river_frc

       ! --------------------------------------------------------------------
       subroutine set_river_frc(time)
         ! Use an analytical forcing or interpolate two instances 
         ! of the forcing that are read from input files
         implicit none
         ! import/export
         real(kind=8), intent(in) :: time
         ! local
         integer :: i,j
         integer :: faces
         logical :: analytical 
         real    :: cf1,cf2
         real    :: tday

         tday = time/(24*3600)
!        if (mynode==0) then
!          print *,'SRF  ', mynode
!          print *,'SRF  ', tday
!          print *,'SRF  ', var_name
!          print *,'SRF  ', var_time
!          stop 'testing set_river_frc'
!        endif
         analytical = .false.
         if (analytical) then
            !! cal ana_river_frc
            stop 'no analytical river data'
         else
           ! check whether we need to update forcing data
           if (frc_tim(it2) < tday) then
             call read_river_frc(tday)
           endif

           ! linear interpolation in time
           cf1 = ( frc_tim(it2) - tday )/(frc_tim(it2) - frc_tim(it1) ) 
           cf2 = ( tday - frc_tim(it1) )/(frc_tim(it2) - frc_tim(it1) ) 
           rflx(:,:,3) = cf1*rflx(:,:,it1) + cf2 * rflx(:,:,it2)

           ! compute uflx,vflx from rflx
           do i = 1,nx
             do j = 1,ny
               if (rflx(i,j,3) > 0) then ! distribute mass flux to all available unmasked cells
                 faces =  rmask(i-1,j)+rmask(i+1,j)+rmask(i,j-1)+rmask(i,j+1) !! amount of unmasked cells around
!                faces = 3
                 if ( faces == 0 .or. rmask(i,j)>0  ) then
                   stop 'river forcing error'
                 endif
                 if (rmask(i-1,j)>0 ) riv_uflx(i  ,j) =-rflx(i,j,3)/faces
                 if (rmask(i+1,j)>0 ) riv_uflx(i+1,j) = rflx(i,j,3)/faces
                 if (rmask(i,j-1)>0 ) riv_vflx(i,j  ) =-rflx(i,j,3)/faces
                 if (rmask(i,j+1)>0 ) riv_vflx(i,j+1) = rflx(i,j,3)/faces
!                print *,'river flux: ',rflx(i,j,3),rflx(i,j,1),faces
!                print *,'river flux: ',cf1, cf2
!                stop
               endif
             enddo
           enddo
     
         endif
       end subroutine set_river_frc
       ! --------------------------------------------------------------------
 
       subroutine find_rec_new(time,v_name,vt_name,ifile,irec,is,js)

         use netcdf
         use ncio

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
! INPUT/OUTPUT
         real(kind=8)    ,intent(in)    :: time
         character(len=*),intent(in)  :: v_name,vt_name
         integer     ,intent(inout) :: ifile,irec
         integer     ,intent(inout) :: is,js
! LOCAL
         real        :: vtime
         integer     :: found_var,found_rec

         integer     :: ierr,v_id,ncid
         integer     :: nfiles,nrecords
         integer     :: ifile_old,irec_old
         integer     :: partition_info(4)

         logical :: first
      
      ! The first time, this routine will return the record just before
      ! 'time'
!        if (mynode==0) then
!          print *, 'ifile, nfiles: ',ifile,nfiles
!          print *, 'Entering FNR: ',ifile,irec
!        endif

         if (ifile==0) then 
           first=.true.
           ifile = 1
           irec_old = 0
         endif

         found_var = 0
         found_rec = 0


         nfiles = max_frc

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
             nrecords= nc_size(frcfile(ifile),vt_name)
             do while (irec<=nrecords.and.found_rec==0)

               call nc_read(frcfile(ifile),vt_name,vtime,[irec],[1])
               vtime = vtime + 360
!              if (mynode==0) then
!                print *,'vtime: ',vtime
!              endif
               if (vtime>time) then
                 found_rec = 1
               else
                 irec_old = irec
                 ifile_old = ifile
                 irec = irec+1
               endif
             enddo ! while not found and irec<nrecs
           else
             if (mynode==0) then
               print *, 'Variable ',v_name,' NOT found in file ',frcfile(ifile)
             endif
           endif

           if (found_rec == 0) then
             !! go to next file if ifile<nfiles
             ifile = ifile+1
           endif

         enddo !while ifile<=nfiles

         call mpi_barrier(ocean_grid_comm,ierr)
         if (found_rec==1) then
           if (mynode==0) then
             print *,'SUCCES: vtime,ifile,irec: ',vtime,ifile,irec
           endif
         else
           if (mynode==0) then
             print *,'FAIL: vtime,ifile,irec: ',vtime,ifile,irec
           endif
         endif

         ! Hack to deal with the fact that currently roms partitions
         ! the buffer in some of the partitioned files
         ierr=nf90_get_att(ncid,nf90_global,'partition',partition_info)
         if (partition_info(3)==1) then
           is = 2;
         else
           is = 1;
         endif
         if (partition_info(4)==1) then
           js = 2;
         else
           js = 1;
         endif
          
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

         if (mynode==0) then
           print *,'FNR: ',found_rec,irec,ifile
         endif

        end subroutine find_rec_new

      end module river_frc
