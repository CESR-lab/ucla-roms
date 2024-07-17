      module river_frc

      ! river forcing module
      ! initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)
      ! contains all the necessary components to produce the riv_uflx,riv_vflx
      ! arrays which have the the horizontal river volume flux in m2/s at the current time
      ! riv_uflx, riv_vflx should be on faces between a masked and unmasked cell,
      ! in the direction of the unmasked cell.

#include "cppdefs.opt"
  
      use dimensions
      use roms_read_write
      use nc_read_write
      use scalars                                                    ! for FIRST_STEP
      use grid

      implicit none

      private

#include "river_frc.opt"

      ! Variables used for equation system calculations
      real,public,allocatable,dimension(:,:) :: riv_uflx
      real,public,allocatable,dimension(:,:) :: riv_vflx
      real,public,allocatable,dimension(:,:) :: rflx                 ! river locations
      real,public,dimension(nriv)    :: riv_vol                      ! river volume
      real,public,dimension(nriv,nt) :: riv_trc                      ! river tracer conc.

      integer,public :: iriver                                       ! river index for looping through rivers
      real,   public :: riv_depth
      real,   public :: riv_uvel,riv_vvel
      real,   public :: river_flux

      ! Netcdf names
      character(len=10) :: riv_flx_name = 'river_flux'               ! stored in the grid file
      character(len=12) :: riv_vol_name = 'river_volume'             ! stored in a forcing file
      character(len=12) :: riv_trc_name = 'river_tracer'             ! stored in a forcing file
      character(len=10) :: riv_tim_name = 'river_time'               ! stored in a forcing file
      character(len=6) :: nriv_dim_name = 'nriver'                   ! dimension name for number of rivers in file
      character(len=8) :: ntrc_dim_name = 'ntracers'                 ! dimension name for number of tracers in file

      ! Misc:
      logical, public :: init_riv_done = .false.                     ! if river variables have been initialized yet

      public set_river_frc
      public init_river_frc

      contains

! ----------------------------------------------------------------------
      subroutine set_river_frc  ![
      ! SET RIVER FORCES (REALISTIC OR ANALYTICAL FORCING):
      ! - read and interpolation all river forcing.
      ! - All river variables need time interpolation only
      !   here so can use same generic routine.
      ! - Input data in days!

      implicit none

      if(.not. init_riv_done) call init_river_frc ! initialize once river flux locations & arrays

      ! set river flux volumes and tracer data:
      if(analytical) then

        call set_ana_river_frc ! cppflags needed else won't link without the analytical.F

      else

        call set_frc_data(nc_rvol,riv_vol)                 ! set river volume flux for all rivers at current time

        call set_frc_data(nc_rtrc,var2d=riv_trc)           ! set river tracers flux for all rivers at current time

      endif

      end subroutine set_river_frc  !]
! ----------------------------------------------------------------------
      subroutine init_river_frc  ![
      ! Initialize river forcing:
      ! Read in a grid file with locations of river mouths and flux contribution per cell.
      ! done only once as river mouth position does not change.
      ! Or: .... if analytical, define the river fluxes in this
      ! subroutine
      use netcdf
      implicit none

      ! local
      integer :: i,j
      integer :: ierr,ncid,varid
      real  :: riv_cells,riv_east,riv_west

      allocate( riv_uflx(GLOBAL_2D_ARRAY) ); riv_uflx = 0.
      allocate( riv_vflx(GLOBAL_2D_ARRAY) ); riv_vflx = 0.
      allocate( rflx(GLOBAL_2D_ARRAY) )    ; rflx = 0.


      if (analytical) then
        riv_west=xl*0.4 ! River west bank at 40% from west
        riv_east=xl*0.6 ! River west bank at 60% from west
        ! pm is constant for this case
        riv_cells = nint( (riv_east - riv_west)*pm(1,1)) !number of cells in this river
        do j=0,ny+1
          do i=0,nx+1
            if (xr(i,j)>riv_west .and. xr(i,j)<riv_east) then
              ! find 'coastline' masked cells
# ifdef MASKING
              if (rmask(i,j)==0 .and. rmask(i,j+1)==1) then
                rflx(i,j) = 1.0+1.0/riv_cells
              endif
# endif
            endif


          enddo
        enddo

        ierr=nf90_open(ana_grdname,nf90_write,ncid)
        varid = nccreate(ncid,'river_flux',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
        ierr=nf90_put_att(ncid, varid,'long_name','River volume flux') 
!       ierr=nf90_close(ncid)
!       print *,'added river_flux',mynode
!       ierr=nf90_open(ana_grdname,nf90_write,ncid)
        call ncwrite(ncid,'river_flux', rflx(i0:i1,j0:j1))
        ierr=nf90_close(ncid)

      else ! Read 'river_flux' from grid file (River mouth locations)

        allocate(nc_rvol%vdata(nriv,1 ,2))
        allocate(nc_rtrc%vdata(nriv,nt,2))
        ierr=nf90_open(grdname, nf90_nowrite, ncid)
        call ncread(ncid, riv_flx_name, rflx(i0:i1,j0:j1) )
        if(ierr/=0) call handle_ierr(ierr,'init_riv:: riv not in file!')
      endif

      call calc_river_flux                                           ! compute uflx,vflx from rflx

      init_riv_done = .true.

      if(mynode==0) write(*,'(/7x,A/)')
     &  'river_frc: init river locations'

      end subroutine init_river_frc  !]
! ----------------------------------------------------------------------
      subroutine calc_river_flux  ![
      ! calculate the river flux contributions to each cell.
      ! river_flux = iriver + fraction of river's flux through grid point.
      ! e.g. River 3 is over 2 grid points (half flux through each point),
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

      end subroutine calc_river_flux  !]
! ----------------------------------------------------------------------
      subroutine set_ana_river_frc  ![
      ! Analytical river forcing volume and tracer data

#include "ana_frc_river.h"

      end subroutine set_ana_river_frc  !]

! ----------------------------------------------------------------------

      end module river_frc
