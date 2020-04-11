      module river_frc
      ! Contains all the neccesary components to produce the riv_uflx,riv_vflx
      ! arrays which have the the horizontal river volume flux in m2/s at the current time
      ! riv_uflx, riv_vflx should be on faces between a masked and unmasked cell,
      ! in the direction of the unmasked cell.

      ! TODO
      ! * Deal with rivers in the buffer that impact boundary uflx,vflx
      ! * we could do the faces thingy in the distribution field
      ! already (better for doing the boundary)
      

      ! call set_river_source_tile(istr,iend,jstr,jend, ierr)
  
       use dimensions !! has dimensions, a list of forcing files, and rmask
       implicit none

       private

      ! These are usable by code that uses this module
       integer,public,parameter :: nriv = 10 !! max number of rivers
       real,public, dimension(nx+1,ny) :: riv_uflx
       real,public, dimension(nx,ny+1) :: riv_vflx
       real,public, dimension(nriv)    :: riv_vol
       real,public, dimension(nriv,nt) :: riv_trc

       integer,public :: iriver
       real,   public :: depth
       real,   public :: riv_uvel,riv_vvel
       real,   public :: river_flux
       real,   public, dimension(nriv,2)   :: rvol
       real,   public, dimension(nriv,nt,2):: rtrc

       ! These are private to this module
       character(len=10) :: riv_flx_name = 'river_flux'  !! stored in the grid file
       character(len=12) :: riv_vol_name = 'river_volume'!! stored in a forcing file
       character(len=12) :: riv_trc_name = 'river_tracer'!! stored in a forcing file
       character(len=10) :: riv_tim_name = 'river_time'  !! stored in a forcing file

       character(len=6) :: nriv_dim_name = 'nriver'      !! dimension name for number of rivers in file
       character(len=8) :: ntrc_dim_name = 'ntracers'    !! dimension name for number of tracers in file

       logical :: analytical = .false.                   !! analytical or forcing files

       integer    :: ifile = 0
       integer    :: irec  = 0
  
       integer                          :: it1 = 1
       integer                          :: it2 = 2
       integer                          :: nriv_frc ! number of rivers in forcing file
       integer                          :: ntrc_frc ! number of tracers in forcing file
       real, dimension(2)               :: rtim = [-99,-99]
       real, dimension(0:nx+1,0:ny+1)   :: rflx  ! include 1 width buffer
       real, parameter :: pi=3.14159265358979323

       public set_river_frc

      contains
       ! --------------------------------------------------------------------
        subroutine ana_river_frc
         ! Analytical river forcing data
        end subroutine ana_river_frc

       ! --------------------------------------------------------------------
        subroutine ini_river_frc
         use ncio
         use netcdf
         implicit none
         ! Initialize river forcing

         integer :: is,js,ierr,ncid
         integer :: partition_info(4)
         integer :: i,j
         integer :: faces

         !! We need to deal with the fact that currently, the first subdomain
         !! includes the west (or south) buffer, which we'll want to skip.

         !! Hack to deal with the fact that currently roms partitions
         !! the buffer in some of the partitioned files

         ierr = nf90_open(grdname, nf90_nowrite, ncid)
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

         call nc_read(grdname,riv_flx_name,rflx(1:nx,1:ny),[is,js],[nx,ny])

         ! compute uflx,vflx from rflx 
         do i = 1,nx
           do j = 1,ny
             if (rflx(i,j) > 0) then ! distribute mass flux to all available unmasked cells
               iriver = floor(rflx(i,j)-1e-5)

               !----DevinD - won't compile if MASKING not defined
               ! Need to discuss with JeroenM how to deal with this.

# ifdef MASKING

               faces =  rmask(i-1,j)+rmask(i+1,j)+rmask(i,j-1)+rmask(i,j+1) !! amount of unmasked cells around
               if ( faces == 0 .or. rmask(i,j)>0  ) then
                  stop 'river forcing error'
               endif
               if (rmask(i-1,j)>0 ) riv_uflx(i  ,j) =-(rflx(i,j)-iriver)/faces + 10*iriver
               if (rmask(i+1,j)>0 ) riv_uflx(i+1,j) = (rflx(i,j)-iriver)/faces + 10*iriver
               if (rmask(i,j-1)>0 ) riv_vflx(i,j  ) =-(rflx(i,j)-iriver)/faces + 10*iriver
               if (rmask(i,j+1)>0 ) riv_vflx(i,j+1) = (rflx(i,j)-iriver)/faces + 10*iriver


# endif

               !----DevinD - end MASKING

               endif

             enddo
           enddo
          
        end subroutine ini_river_frc


       ! --------------------------------------------------------------------
       subroutine read_river_frc(time)
         ! Reads river forcing data from input file
         ! The first time, it reads two instances
         use ncio
         implicit none

         real(kind=8), intent(in) :: time
         integer  :: tmp


         if (rtim(it1) < 0) then !! first time
           call ini_river_frc

           ! calling find_rec_new with irec==0 returns the last record with var_time<time
           call find_rec_new(time,riv_vol_name,riv_tim_name,ifile,irec)
           call nc_read(frcfile(ifile),riv_tim_name,rtim(it1),[irec],[1])
           rtim(it1) = rtim(it1) + 360  !! place holder for proper cycle_length treatment

           nriv_frc = nc_size(frcfile(ifile),nriv_dim_name)
           ntrc_frc = nc_size(frcfile(ifile),ntrc_dim_name)
           if (ntrc_frc>nt) then
             ntrc_frc= nt;
             print *, 'Ignoring some tracers in river input'
           endif

           if (nriv_frc>nriv) then
             if (mynode==0) then
               print *, 'Too many rivers in forcing file,'
               print *, 'Increase nriv in river_frc.F'
             endif
             stop
           endif

           call nc_read(frcfile(ifile),riv_vol_name,rvol(1:nriv_frc,it1),[1,irec],[nriv_frc,1])
           call nc_read(frcfile(ifile),riv_trc_name,rtrc(1:nriv_frc,1:ntrc_frc,it1),[1,1,irec],[nriv_frc,ntrc_frc,1])

         else
           tmp =it1
           it1 =it2
           it2 =tmp
         endif

         call find_rec_new(time,riv_vol_name,riv_tim_name,ifile,irec)

         call nc_read(frcfile(ifile),riv_tim_name,rtim(it2),[irec],[1])
         rtim(it2) = rtim(it2) + 360  !! place holder for proper cycle_length'

         call nc_read(frcfile(ifile),riv_vol_name,rvol(1:nriv_frc,it2),[1,irec],[nriv_frc,1])
         call nc_read(frcfile(ifile),riv_trc_name,rtrc(1:nriv_frc,1:ntrc_frc,it2),[1,1,irec],[nriv_frc,ntrc_frc,1])

         if (mynode==0) then
            print *,'read_river_frc: ', rtim(it1),rtim(it2)
         endif

       end subroutine read_river_frc

       ! --------------------------------------------------------------------
       subroutine set_river_frc(time)
         ! Use an analytical forcing or interpolate two instances 
         ! of the forcing that are read from input files
         implicit none
         ! import/export
         real(kind=8), intent(in) :: time
         ! local
         real    :: cf1,cf2
         real    :: tday

         tday = time/(24*3600) ! convert from time seconds to days

         if (analytical) then
            !! cal ana_river_frc
            !! fill riv_uflx, riv_vflx (river distribution)
            !! fill riv_vol, riv_trc (volume and tracer value data)
            stop 'no analytical river data'
         else
           ! check whether we need to update forcing data
           if (rtim(it2) < tday) then
             call read_river_frc(tday)
           endif

           ! linear interpolation in time
           cf1 = ( rtim(it2) - tday )/(rtim(it2) - rtim(it1) ) 
           cf2 = ( tday - rtim(it1) )/(rtim(it2) - rtim(it1) ) 
           riv_vol = cf1*rvol(:,it1) + cf2 * rvol(:,it2)
           riv_trc = cf1*rtrc(:,:,it1) + cf2 * rtrc(:,:,it2)
!          if (mynode==0) then
!            print *, 'volume: ',riv_vol(1:nriv_frc)
!          endif
           ! tidal modulation of river volume
           riv_vol = (cos(2*pi*tday) + 1)*riv_vol

     
         endif
       end subroutine set_river_frc
       ! --------------------------------------------------------------------
 
       subroutine find_rec_new(time,v_name,vt_name,ifile,irec) 

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
! LOCAL
         real        :: vtime
         integer     :: found_var,found_rec

         integer     :: ierr,v_id,ncid
         integer     :: nfiles,nrecords
         integer     :: ifile_old,irec_old

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
               print *, 'Variable ',v_name,' NOT found in file:',frcfile(ifile)
             endif
           endif

           if (found_rec == 0) then
             !! go to next file if ifile<nfiles
             ifile = ifile+1
           endif

         enddo !while ifile<=nfiles

!        call mpi_barrier(ocean_grid_comm,ierr)
         if (found_rec==1) then
           if (mynode==0) then
!            print *,'SUCCES: vtime,ifile,irec: ',vtime,ifile,irec
           endif
         else
           if (mynode==0) then
!            print *,'FAIL: vtime,ifile,irec: ',vtime,ifile,irec
           endif
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
           print *,'FNR: ',frcfile(ifile),irec
         endif

        end subroutine find_rec_new

      end module river_frc
