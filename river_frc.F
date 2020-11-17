      module river_frc

      ! ====================
      ! RIVER FORCING MODULE
      ! ====================

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)

#include "cppdefs.h"
#ifdef RIVER_SOURCE
      ! Contains all the neccesary components to produce the riv_uflx,riv_vflx
      ! arrays which have the the horizontal river volume flux in m2/s at the current time
      ! riv_uflx, riv_vflx should be on faces between a masked and unmasked cell,
      ! in the direction of the unmasked cell.
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use read_write

      implicit none

      private

      ! ****************************************************************
      ! user inputs
      integer,public,parameter :: nriv = 1 ! number of rivers
      logical :: analytical = .false.      ! analytical or forcing files
      ! end user inputs
      ! ****************************************************************

      ! Variables used for equation system calculations
      real,public, dimension(GLOBAL_2D_ARRAY) :: riv_uflx
      real,public, dimension(GLOBAL_2D_ARRAY) :: riv_vflx
      real,public, dimension(nriv)     :: riv_vol  ! River volume
      real,public, dimension(nriv,nt)  :: riv_trc  ! River tracer conc.
      real, dimension(GLOBAL_2D_ARRAY) :: rflx  ! Raw read in river locations

      integer,public :: iriver
      real,   public :: depth
      real,   public :: riv_uvel,riv_vvel
      real,   public :: river_flux

      ! Netcdf names
      character(len=10) :: riv_flx_name = 'river_flux'  !! stored in the grid file
      character(len=12) :: riv_vol_name = 'river_volume'!! stored in a forcing file
      character(len=12) :: riv_trc_name = 'river_tracer'!! stored in a forcing file
      character(len=10) :: riv_tim_name = 'river_time'  !! stored in a forcing file
      character(len=6) :: nriv_dim_name = 'nriver'      !! dimension name for number of rivers in file
      character(len=8) :: ntrc_dim_name = 'ntracers'    !! dimension name for number of tracers in file
      ! Netcdf variables used for both riv_vol & riv_trc
      real, dimension(nriv,2) :: riv_vol_data ! rvol - river data at input times
      real, dimension(nriv,nt,2) :: riv_trc_data !rtrc - river tracer data at input times
      integer, dimension(2) :: r_file_indx = 0
      integer, dimension(2) :: r_irec      = 0
      integer, dimension(2) :: r_it1 = 1, r_it2 = 2
      real, dimension(2,2)  :: r_times = -99 ! [-99,-99]

      public set_river_frc_tile

      contains
! ----------------------------------------------------------------------

      subroutine ana_river_frc
         ! Analytical river forcing data
      end subroutine ana_river_frc
! ----------------------------------------------------------------------

      subroutine init_river_frc(istr,iend,jstr,jend)

      ! Initialize river forcing:
      ! Read in a grid file with locations of river mouths.
      ! Done only once as river mouth positon does not change.
      ! river_flux = iriver + fraction of river's flux through grid point.
      ! E.g. River 3 is over 2 grid points (half flux through each point),
      ! hence river_flux = 3 + 0.5 = 3.5.

      use netcdf

      implicit none

      ! inputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer :: ierr,ncid
      integer :: i,j, v_id
      integer :: faces

      ! Read 'river_flux' (River mouth locations)
      ! Should put this in wrapper in read_write.F as 'nc_read_once'
      ierr = nf90_open(grdname, nf90_nowrite, ncid)
      ierr = nf90_inq_varid(ncid, riv_flx_name, v_id)
      ierr = ncdf_read_mod(ncid, v_id, 1, rp_var, rflx(START_2D_ARRAY), 1)
      if(ierr/=0) call handle_ierr(ierr,'init_river:: not in gridfile!')

      ! compute uflx,vflx from rflx
      do j = jstr,jend
        do i = istr,iend
          if (rflx(i,j) > 0) then ! distribute mass flux to all available unmasked cells
            ! subtract 1e-5 in case only 1 grid point for river, so that floor still
            ! produces correct iriver number.
            iriver = floor(rflx(i,j)-1e-5)
#ifdef MASKING
            faces =  rmask(i-1,j)+rmask(i+1,j)+rmask(i,j-1)+rmask(i,j+1) !! amount of unmasked cells around
            if ( faces == 0 .or. rmask(i,j)>0  ) then
              error stop 'river grid position error'
            endif
            ! 10*iriver needed because uflx/vflx can be positive or negative around
            ! the iriver number, and hence nearest integer is safest done with 10*.
            if (rmask(i-1,j)>0 ) riv_uflx(i  ,j) =-(rflx(i,j)-iriver)/faces + 10*iriver
            if (rmask(i+1,j)>0 ) riv_uflx(i+1,j) = (rflx(i,j)-iriver)/faces + 10*iriver
            if (rmask(i,j-1)>0 ) riv_vflx(i,j  ) =-(rflx(i,j)-iriver)/faces + 10*iriver
            if (rmask(i,j+1)>0 ) riv_vflx(i,j+1) = (rflx(i,j)-iriver)/faces + 10*iriver
#endif
          endif
        enddo
      enddo

      if(mynode==0) print*, 'rivers: read in river locations'
          
      end subroutine init_river_frc
! ----------------------------------------------------------------------

      subroutine set_river_frc_tile(istr,iend,jstr,jend)

      ! SET RIVER FORCES (REALISTIC OR ANALYTICAL FORCING):
      ! - read and interpolation all river forcing.
      ! - All river variables need time interpolation only
      !   here so can use same generic routine.
      ! - Input data in days!

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local (scalars for FIRST_TIME_STEP to work)
#include "scalars.h"

      ! Read in river flux locations (once off)
      if(FIRST_TIME_STEP) call init_river_frc(istr,iend,jstr,jend)

      ! Set river volume flux for all rivers at current time
      call set_small_arrays(
     &          riv_vol_name,   riv_tim_name, ! Text names
     &          riv_vol_data,   riv_vol,      ! Variable arrays
     &          r_times(:,1),                 ! Input times
     &          nriv,           1,            ! dim1 and dim2=1 if 1D
     &          r_file_indx(1), r_irec(1),    ! File indx & Current input record
     &          r_it1(1),       r_it2(1) )    ! Time index placeholders                   ! Online spatial interpolation flag

!      if(mynode==0) print*, 'riv_vol= ', riv_vol ! debug

      ! Set river tracers flux for all rivers at current time
      call set_small_arrays(
     &          riv_trc_name,   riv_tim_name, ! Text names
     &          riv_trc_data,   riv_trc,      ! Variable arrays
     &          r_times(:,2),                 ! Input times
     &          nriv,           nt,           ! dim1 and dim2
     &          r_file_indx(2), r_irec(2),    ! File indx & Current input record
     &          r_it1(2),       r_it2(2) )    ! Time index placeholders

!      if(mynode==0) print*, 'riv_trc(1,1)= ', riv_trc(1,1) ! debug
!      if(mynode==0) print*, 'riv_trc(1,2)= ', riv_trc(1,2) ! debug

      end subroutine set_river_frc_tile
! ----------------------------------------------------------------------

!      subroutine set_river_frc_tile_old_DPD(istr,iend,jstr,jend,ierr,time)
!
!      ! ================
!      ! SET RIVER FORCES
!      ! (REALISTIC OR ANALYTICAL FORCING)
!      ! =================================
!
!      ! Use an analytical forcing or interpolate two instances
!      ! of the forcing that are read from input files
!      implicit none
!
!      ! input/outputs
!      integer istr,iend,jstr,jend ! ,intent(in) ::
!      integer ierr ! ,intent(out)::
!      real(kind=8), intent(in) :: time ! Model time in seconds
!      ! local
!      integer :: tmp,i,j
!      real    :: cf1,cf2
!      real    :: tday
!
!      ! convert model time from seconds to days
!      ! as input data is currently in days
!      tday = time/(24*3600)
!
!      ! IF ANALYTICAL FORCING
!      ! =====================
!# ifdef ANA_RIVER_SOURCE
!!      if (analytical) then
!      !! cal ana_river_frc
!      !! fill riv_uflx, riv_vflx (river distribution)
!      !! fill riv_vol, riv_trc (volume and tracer value data)
!      stop 'no analytical river data'
!
!
!# else /* REALISTIC MODEL */
!
!      ! IF REALISTIC FORCING
!      ! ====================
!
!      ! note: rtim(it2) > rtim(it1), but values of it1
!      ! and it2 change between each reading of values so that only one
!      ! value is read in.
!
!      ! UPDATE LATEST READ IN VARIABLE?
!      ! -------------------------------
!      ! Need rtim(it1) < time < rtim(it2), thus:
!      if (rtim(it2) < time) then  ! We need to refresh data
!
!        ! FIRST TIMESTEP
!        ! --------------
!        ! rtim(it1) & rtim(its) declared < 0 at top of module:
!        if (rtim(it1) < 0) then ! First timestep
!
!          ! ifile and irec should be in the module pre-amble, initial as 0
!          ! Find the right force file (ifile) and the record (irec) for
!          ! the last record with rtime < tday
!          ! First time, we need to read an extra slice of data
!          ! so this routine provides irec = irec - 1 for first time.
!          ! Must make sure model time format (secs/days) matches data time
!          ! Hence tday as input data times are in days
!          call find_rec_new_mod(tday,riv_flx_name,riv_tim_name,
!     &                               ifile,irec)
!
!          ! Reading is only done by master processor hence C$OMP MASTER
!C$OMP MASTER
!          ! For the first time irec = irec - 1, as above, such that
!          ! the closest to rtim(it1) < time.
!          ! read in raw input variables and
!          ! read and set input time in rtim(it2) in seconds
!          call read_river_frc(frcfile(ifile),irec,it1,rtim(it1))
!C$OMP END MASTER
!C$OMP BARRIER
!          ! Barrier needed because master processor must have read in all
!          ! WEC data before child processors try to access read in data!
!
!          ! Confirm to terminal reading of variables and
!          ! time for read in data
!          if (mynode == 0) then
!            write(*,'(4x,A,5x,A,G14.6,1x,A,I4)')
!     &               'set_river_frc :: read RIVER inputs',
!     &               'river_input_time (s) =',rtim(it1),
!     &               'rec =',irec MYID
!          endif
!
!          ! Convert read inputs into useable, masked variables
!          ! for the first, earlier read in time variables(time) < time.
!          call set_river_inp2vars (istr,iend,jstr,jend,it1)
!
!        ! ALL OTHER TIMESTEPS
!        ! -------------------
!        ! Every time step switch indices such that
!        ! var(i,j,1) & rtim(1) -> earlier && var(i,j,2) & rtim(2) -> later: step t=x
!        ! var(i,j,1) & rtim(1) -> later   && var(i,j,2) & rtim(2) ->  earlier: step t=x+1
!        ! and so on. The benifit is only need to change one variable each time.
!        else
!          tmp = it1
!          it1 = it2 ! If it1 = 1, it now equals 2, and vice-versa
!          it2 = tmp ! If it2 = 2, it now equals 1, and vice-versa
!
!        endif ! End of first step
!
!        ! Find the force file (ifile) and the record (irec) for
!        ! the first record with input time > model time
!        ! Must make sure model time format (secs/days) matches data time
!        ! Hence tday as input data times are in days
!        call find_rec_new_mod(tday,riv_flx_name,riv_tim_name,
!     &                             ifile,irec)
!
!          ! Causes heisenbug - WEC but test with rivers? Might be read issue
!          ! if (mynode==0) print *,'irec=',irec ! Debug
!
!        ! Reading is only done by master processor hence C$OMP MASTER
!C$OMP MASTER
!        ! Read in the first record with wec_time > tday
!        ! Remember, it2 alternates between index 1 & 2 for each reading
!        ! read in raw input variables and
!        ! read and set input time in rtim(it2) to second
!        call read_river_frc(frcfile(ifile),irec,it2,rtim(it2))
!C$OMP END MASTER
!C$OMP BARRIER
!
!        ! Confirm to terminal reading of variables and
!        ! time for read in data
!        if (mynode == 0) then                  ! Text formatting:
!          write(*,'(4x,A,5x,A,G14.6,1x,A,I4)') ! 4x is spaces, A is string
!     &               'set_river_frc :: read RIVER inputs', ! 5x is 5 spaces, A is string
!     &               'rvier_input_time (s) =',rtim(it2), ! ! G is number for time
!     &               'rec =',irec MYID ! I is integer
!        endif
!
!        ! Convert read inputs into useable, rotated & rmask'ed variables
!        ! for the latter read in time variables(time) > time.
!        call set_river_inp2vars (istr,iend,jstr,jend,it2)
!
!      endif
!
!      ! Temporal interpolation
!      ! ======================
!
!
!      ! linear interpolation in time
!      cf1 = ( rtim(it2) - tday )/(rtim(it2) - rtim(it1) )
!      cf2 = ( tday - rtim(it1) )/(rtim(it2) - rtim(it1) )
!      riv_vol = cf1*rvol(:,it1) + cf2 * rvol(:,it2)
!      riv_trc = cf1*rtrc(:,:,it1) + cf2 * rtrc(:,:,it2)
!!          if (mynode==0) then
!!            print *, 'volume: ',riv_vol(1:nriv_frc)
!!          endif
!      ! tidal modulation of river volume
!      riv_vol = (cos(2*pi*tday) + 1)*riv_vol
!
!
!!      endif ! for if (analytical) then
!
!# endif /* ANA_RIVER_SOURCE or REALISTIC MODEL */
!
!      end subroutine set_river_frc_tile_old_DPD


! ----------------------------------------------------------------------
!       subroutine read_river_frc_old(time)
!         ! Reads river forcing data from input file
!         ! The first time, it reads two instances
!         use ncio
!         implicit none
!
!         real(kind=8), intent(in) :: time
!         integer  :: tmp
!
!
!         if (rtim(it1) < 0) then !! first time
!           call ini_river_frc
!
!           ! calling find_rec_new with irec==0 returns the last record with var_time<time
!           call find_rec_new(time,riv_vol_name,riv_tim_name,ifile,irec)
!           call nc_read(frcfile(ifile),riv_tim_name,rtim(it1),[irec],[1])
!           rtim(it1) = rtim(it1) + 360  !! place holder for proper cycle_length treatment
!
!           nriv_frc = nc_size(frcfile(ifile),nriv_dim_name)
!           ntrc_frc = nc_size(frcfile(ifile),ntrc_dim_name)
!           if (ntrc_frc>nt) then
!             ntrc_frc= nt;
!             print *, 'Ignoring some tracers in river input'
!           endif
!
!           if (nriv_frc>nriv) then
!             if (mynode==0) then
!               print *, 'Too many rivers in forcing file,'
!               print *, 'Increase nriv in river_frc.F'
!             endif
!             stop
!           endif
!
!           call nc_read(frcfile(ifile),riv_vol_name,rvol(1:nriv_frc,it1),[1,irec],[nriv_frc,1])
!           call nc_read(frcfile(ifile),riv_trc_name,rtrc(1:nriv_frc,1:ntrc_frc,it1),[1,1,irec],[nriv_frc,ntrc_frc,1])
!
!         else
!           tmp =it1
!           it1 =it2
!           it2 =tmp
!         endif
!
!         call find_rec_new(time,riv_vol_name,riv_tim_name,ifile,irec)
!
!         call nc_read(frcfile(ifile),riv_tim_name,rtim(it2),[irec],[1])
!         rtim(it2) = rtim(it2) + 360  !! place holder for proper cycle_length'
!
!         call nc_read(frcfile(ifile),riv_vol_name,rvol(1:nriv_frc,it2),[1,irec],[nriv_frc,1])
!         call nc_read(frcfile(ifile),riv_trc_name,rtrc(1:nriv_frc,1:ntrc_frc,it2),[1,1,irec],[nriv_frc,ntrc_frc,1])
!
!         if (mynode==0) then
!            print *,'read_river_frc: ', rtim(it1),rtim(it2)
!         endif
!
!       end subroutine read_river_frc_old

! ----------------------------------------------------------------------
 
!       subroutine find_rec_new(time,v_name,vt_name,ifile,irec)
!
!         use netcdf
!         use ncio
!
!! Find the first time record number "irec" and its forcing file
!! "ifile" for variable "var", with accompanying time "var_time".
!! for which var_time > ocean_time.
!!
!! When called for the first time, it will the last record for which
!! var_time < ocean_time
!!
!! On entry, ifile and irec should contain results from the last call to
!! find_record for this variable, or zeros if called for the first time
!
!!------------
!
!!         frcfile    list of forcing file names (from module dimensions)
!!         time       model time
!!         var_name   name of variable in file
!!         var_time   name of corresponding time variable in file
!
!! in/out: ifile      index of correct file in frcfile
!!         irec       record number in file
!
!         implicit none
!! INPUT/OUTPUT
!         real(kind=8)    ,intent(in)    :: time
!         character(len=*),intent(in)  :: v_name,vt_name
!         integer     ,intent(inout) :: ifile,irec
!! LOCAL
!         real        :: vtime
!         integer     :: found_var,found_rec
!
!         integer     :: ierr,v_id,ncid
!         integer     :: nfiles,nrecords
!         integer     :: ifile_old,irec_old
!
!         logical :: first
!
!      ! The first time, this routine will return the record just before
!      ! 'time'
!!        if (mynode==0) then
!!          print *, 'ifile, nfiles: ',ifile,nfiles
!!          print *, 'Entering FNR: ',ifile,irec
!!        endif
!
!         if (ifile==0) then
!           first=.true.
!           ifile = 1
!           irec_old = 0
!         endif
!
!         found_var = 0
!         found_rec = 0
!
!         nfiles = max_frc
!
!         do while ((found_rec==0).and.(ifile<=nfiles))
!
!           if (mynode==0) then
!!            print *, 'ifile, nfiles: ',ifile,nfiles
!             print *, 'Searching for ',v_name,' in file ',frcfile(ifile)
!           endif
!
!           ierr = nf90_open(frcfile(ifile), nf90_nowrite, ncid)
!           if (ierr.ne.0) then ! This is not supposed to happen
!             if (mynode==0) then
!               print *,ifile,frcfile(ifile)
!               stop 'FNR: Forcing file not found'
!             endif
!           endif
!
!           ! Check if this file contains the variable
!           ierr = nf90_inq_varid(ncid, v_name, v_id)
!           if (ierr.eq.0) then ! Variable found, now check times
!             if (mynode==0) then
!               print *, 'Variable ',v_name,' found in file ',frcfile(ifile)
!             endif
!             found_var = 1
!             irec = 1
!             nrecords= nc_size(frcfile(ifile),vt_name)
!             do while (irec<=nrecords.and.found_rec==0)
!
!               call nc_read(frcfile(ifile),vt_name,vtime,[irec],[1])
!               vtime = vtime + 360
!!              if (mynode==0) then
!!                print *,'vtime: ',vtime
!!              endif
!               if (vtime>time) then
!                 found_rec = 1
!               else
!                 irec_old = irec
!                 ifile_old = ifile
!                 irec = irec+1
!               endif
!             enddo ! while not found and irec<nrecs
!           else
!             if (mynode==0) then
!               print *, 'Variable ',v_name,' NOT found in file:',frcfile(ifile)
!             endif
!           endif
!
!           if (found_rec == 0) then
!             !! go to next file if ifile<nfiles
!             ifile = ifile+1
!           endif
!
!         enddo !while ifile<=nfiles
!
!!        call mpi_barrier(ocean_grid_comm,ierr)
!         if (found_rec==1) then
!           if (mynode==0) then
!!            print *,'SUCCES: vtime,ifile,irec: ',vtime,ifile,irec
!           endif
!         else
!           if (mynode==0) then
!!            print *,'FAIL: vtime,ifile,irec: ',vtime,ifile,irec
!           endif
!         endif
!
!         if (first) then
!           ifile = ifile_old
!           irec  = irec_old
!         endif
!
!         if (found_var==0) then
!!          write(*,'(/1x,2A/)')  '### ERROR: find_record: Could not find variable ', var_name
!           stop 'variable not found in forcing files'
!         endif
!         if (found_rec==0) then
!!          write(*,'(/1x,2A/)')  '### ERROR: find_record: Ran out of time records for ', var_name
!           stop 'ran out of time records'
!         endif
!         if (irec==0) then
!           stop 'First available forcing record is past current time'
!         endif
!
!         if (mynode==0) then
!           print *,'FNR: ',frcfile(ifile),irec
!         endif
!
!        end subroutine find_rec_new

# endif
      end module river_frc
