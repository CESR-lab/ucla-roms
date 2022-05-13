      module river_frc

      ! RIVER FORCING MODULE
      ! --------------------

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)

#include "cppdefs.opt"
#if defined RIVER_SOURCE

      ! Contains all the necessary components to produce the riv_uflx,riv_vflx
      ! arrays which have the the horizontal river volume flux in m2/s at the current time
      ! riv_uflx, riv_vflx should be on faces between a masked and unmasked cell,
      ! in the direction of the unmasked cell.
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use read_write
      use scalars    ! for FIRST_STEP

      implicit none

      private

#include "river_frc.opt"

      ! Variables used for equation system calculations
      real,public,allocatable,dimension(:,:) :: riv_uflx
      real,public,allocatable,dimension(:,:) :: riv_vflx
      real,public,allocatable,dimension(:,:) :: rflx         ! River locations
      real,public,dimension(nriv)    :: riv_vol              ! River volume
      real,public,dimension(nriv,nt) :: riv_trc              ! River tracer conc.

      integer,public :: iriver            ! River index for looping through rivers
      real,   public :: riv_depth
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
      real, dimension(nriv,2)    :: riv_vol_data ! rvol - river data at input times
      real, dimension(nriv,nt,2) :: riv_trc_data ! rtrc - river tracer data at input times
      integer, dimension(2) :: r_file_indx = 0
      integer, dimension(2) :: r_irec      = 0
      integer, dimension(2) :: r_it1 = 1, r_it2 = 2
      real, dimension(2,2)  :: r_times = -99 ! [-99,-99]

      ! Misc:
      logical, public :: init_riv_done = .false.  ! if river variables have been initialized yet

      public set_river_frc
      public init_river_frc

      contains

! ----------------------------------------------------------------------
      subroutine set_river_frc

      ! SET RIVER FORCES (REALISTIC OR ANALYTICAL FORCING):
      ! - read and interpolation all river forcing.
      ! - All river variables need time interpolation only
      !   here so can use same generic routine.
      ! - Input data in days!

      implicit none

      if(.not. init_riv_done) call init_river_frc ! initialize once river flux locations & arrays

      ! Set river flux volumes and tracer data
      if(analytical) then
!       if(mynode==0) print*, 'rivers: going to set_ana'

        call set_ana_river_frc ! cppflags needed else won't link without the analytical.F

      else

        ! Set river volume flux for all rivers at current time
        call set_small_arrays(
     &          riv_vol_name,   riv_tim_name, ! Text names
     &          riv_vol_data,   riv_vol,      ! Variable arrays
     &          r_times(:,1),                 ! Input times
     &          nriv,           1,            ! dim1 and dim2=1 if 1D
     &          r_file_indx(1), r_irec(1),    ! File indx & Current input record
     &          r_it1(1),       r_it2(1) )    ! Time index placeholders                   ! Online spatial interpolation flag

!        if(mynode==0) print*, 'riv_vol= ', riv_vol ! debug

        ! Set river tracers flux for all rivers at current time
        call set_small_arrays(
     &          riv_trc_name,   riv_tim_name, ! Text names
     &          riv_trc_data,   riv_trc,      ! Variable arrays
     &          r_times(:,2),                 ! Input times
     &          nriv,           nt,           ! dim1 and dim2
     &          r_file_indx(2), r_irec(2),    ! File indx & Current input record
     &          r_it1(2),       r_it2(2) )    ! Time index placeholders


      endif

      end subroutine set_river_frc

! ----------------------------------------------------------------------
      subroutine init_river_frc

      ! Initialize river forcing:
      ! Read in a grid file with locations of river mouths and flux contribution per cell.
      ! Done only once as river mouth position does not change.

      use netcdf
      implicit none

      ! local
      integer :: ierr,ncid,v_id

      allocate( riv_uflx(GLOBAL_2D_ARRAY) )
      allocate( riv_vflx(GLOBAL_2D_ARRAY) )
      if(.not. allocated(rflx))                  ! for analytical it is already allocated.
     &  allocate( rflx(GLOBAL_2D_ARRAY) )
      riv_uflx = 0.; riv_vflx = 0.
!      rflx = 0.

      ! Should put this in wrapper in read_write.F as 'nc_read_once'
      if (analytical) then

        ! rflx is defined in ana_grid

      else ! Read 'river_flux' from grid file (River mouth locations)
        ierr=nf90_open(grdname, nf90_nowrite, ncid)
        ierr=nf90_inq_varid(ncid, riv_flx_name, v_id)
        v_id = -v_id  ! set v_id to negative so that ncdf_read_mod knows not to mask!
        ierr=ncdf_read_mod(ncid, v_id, 0, 0, rflx(START_2D_ARRAY), 1)
        if(ierr/=0) call handle_ierr(ierr,'init_riv:: riv not in file!')

      endif

      ! Compute uflx,vflx from rflx
      call calc_river_flux

      init_riv_done = .true.

      if(mynode==0) print*, 'rivers: init river locations'
          
      end subroutine init_river_frc

! ----------------------------------------------------------------------
      subroutine calc_river_flux
      ! Calculate the river flux contributions to each cell.
      ! river_flux = iriver + fraction of river's flux through grid point.
      ! E.g. River 3 is over 2 grid points (half flux through each point),
      ! hence river_flux = 3 + 0.5 = 3.5
      implicit none

      ! local
      integer :: i,j,faces

      ! compute uflx,vflx from rflx
      do j = 0,ny+1   ! Loop over -1 and +1 because rflx cell only flows into
        do i = 0,nx+1 ! neighbour, hence cell next to boundary could flow into cell.
          if (rflx(i,j) > 0) then ! distribute mass flux to all available unmasked cells
            ! subtract 1e-5 in case only 1 grid point for river, so that floor still
            ! produces correct iriver number.
!            write(*,*) 'mynode=',mynode,'i,j',i,j,rflx(i,j),'rflx(i,j)'
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

      end subroutine calc_river_flux

! ----------------------------------------------------------------------
      subroutine set_ana_river_frc
      ! Analytical river forcing volume and tracer data

#include "river_frc_ana.opt"

      end subroutine set_ana_river_frc

! ----------------------------------------------------------------------

# endif /* ifdef RIVER_SOURCE */
      end module river_frc
