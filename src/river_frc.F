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
      real   ,public,allocatable,dimension(:,:) :: riv_uflx
      real   ,public,allocatable,dimension(:,:) :: riv_vflx
      real   ,public,allocatable,dimension(:,:) :: rflx ! river locations
      real   ,public,allocatable,dimension(:,:) :: rfrc ! River fraction
      real   ,public,allocatable,dimension(:,:) :: ridx_real ! River indices (read as real by ncread)
      integer,public,allocatable,dimension(:,:) :: ridx ! River indices (stored as int by ROMS)
      
      real,public,dimension(nriv)    :: riv_vol                      ! river volume
      real,public,dimension(nriv,nt) :: riv_trc                      ! river tracer conc.

!integer,public :: iriver  ! river index for looping through rivers
      integer,public :: iriver                                       ! river index for looping through rivers      
      real,   public :: riv_depth
      real,   public :: riv_uvel,riv_vvel
      real,   public :: river_flux

      ! Netcdf names
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

      if (.not. init_riv_done) then
         allocate(nc_rvol%vdata(nriv,1 ,2))
         allocate(nc_rtrc%vdata(nriv,nt,2))
      end if

      ! set river flux volumes and tracer data:
      if(analytical) then

        call set_ana_river_frc ! cppflags needed else won't link without the analytical.F

      else
         call set_frc_data(nc_rvol,riv_vol) ! set river volume flux for all rivers at current time
        call set_frc_data(nc_rtrc,var2d=riv_trc)           ! set river tracers flux for all rivers at current time

      endif

      if(.not. init_riv_done) call init_river_frc ! initialize once river flux locations & arrays
      
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
      allocate( ridx(GLOBAL_2D_ARRAY) )    ; ridx = 0
      allocate( ridx_real(GLOBAL_2D_ARRAY) )    ; ridx_real = 0.
      allocate( rfrc(GLOBAL_2D_ARRAY) )    ; rfrc = 0.      

      
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

      else ! Not analytical, read from file

!     First, check forcing file for separate index and fraction variables
         ierr=nf90_open(frcfile(nc_rvol%ifile), nf90_nowrite, ncid) ! open river forcing file
         ierr = nf90_inq_varid(ncid, "river_index", varid) ! check river forcing file for index...
         ierr = ierr * nf90_inq_varid(ncid, "river_fraction", varid) ! ... and fraction variables
         
         if (ierr == nf90_noerr) then ! Found the variables in the forcing file
            call ncread(ncid, "river_index", ridx_real(i0:i1,j0:j1))
            call ncread(ncid, "river_fraction", rfrc(i0:i1,j0:j1))
            ierr = nf90_close(ncid)            
!     Check for non-integer values
            ridx(i0:i1, j0:j1) = int(ridx_real(i0:i1, j0:j1))
            if (any(abs(ridx_real(i0:i1,j0:j1) - ridx(i0:i1,j0:j1)) > 1.0e-6)) then
               write(*,*) 'ERROR: river_index contains non-integers!'
               stop
            endif
!     Convert to integer
         else !     Fallback, check grid for combined index and fraction variable (deprecated):
            ierr = nf90_close(ncid) ! Close the forcing file

            ierr=nf90_open(grdname, nf90_nowrite, ncid) ! open grid file
            ierr = nf90_inq_varid(ncid, "river_flux", varid) ! check grid file for variable
         
            if (ierr /= nf90_noerr) then ! True if required vars in neither file
               call handle_ierr(ierr, 'init_riv:: '//
     &              'unable to find a single variable (`river_flux`)'//
     &              ' or separate variables '//
     &              '(`river_index`,`river_fraction`) in either the'//
     &              ' grid (' // trim(grdname) //') '//
     &              ' or river forcing (' // trim(frcfile(nc_rvol%ifile)) //
     &              ') files. ROMS is unable to determine '//
     &              'how to distribute river inputs')
            else                !It's in one of the two files (current ncid)
               if (mynode == 0) then
                  write(*,'(/7x,A/)') 'DEPRECATION WARNING, init_riv: '//
     &                 'ROMS has located a single variable '//
     &                 '(river_flux) containing both river index '//
     &                 ' and fraction. This behavior is deprecated '//
     &                 'as of March 2025 and will be removed in a '//
     &                 'future version. Please provide separate ' //
     &                 'variables: river_index and river_fraction.'
               end if

!     Read from whichever of the two files contains the variable
               call ncread(ncid, "river_flux", rflx(i0:i1,j0:j1))
               ierr = nf90_close(ncid)
!     Split the variable up here for later use
!     ridx(i0:i1,j0:j1) = floor(rflx(i0:i1,j0:j1)-1e-5)
               where (rflx(i0:i1, j0:j1) > 0)
                  ridx(i0:i1, j0:j1) = floor(rflx(i0:i1, j0:j1) - 1e-5)
                  rfrc(i0:i1,j0:j1) = rflx(i0:i1,j0:j1) - ridx(i0:i1,j0:j1)                  
               elsewhere
                  ridx(i0:i1, j0:j1) = 0
                  rfrc(i0:i1, j0:j1) = 0
               end where
            end if
         end if                 ! separate river_index, river_fraction vars (T) single river_flux var (F)
      endif                     ! analytical

      call calc_river_flux                                           ! compute uflx,vflx from rflx

      init_riv_done = .true.

      if(mynode==0) write(*,'(/7x,A/)')
     &  'river_frc: init river locations'

      end subroutine init_river_frc  !]
! ----------------------------------------------------------------------
      subroutine calc_river_flux  ![
      ! calculate the river flux contributions to each cell.
      ! river_flux = index + fraction of river's flux through grid point.
      ! e.g. River 3 is over 2 grid points (half flux through each point),
      ! hence river_flux = 3 + 0.5 = 3.5
      implicit none

      ! local
      integer :: i,j,faces

! compute uflx,vflx from rflx
      do j = 0,ny+1   ! Loop over -1 and +1 because rflx cell only flows into
        do i = 0,nx+1 ! neighbour, hence cell next to boundary could flow into cell.
          if (rfrc(i,j) > 0) then ! distribute mass flux to all available unmasked cells
            ! subtract 1e-5 in case only 1 grid point for river, so that floor still
            ! produces correct iriver number.
!            write(*,*) 'mynode=',mynode,'i,j',i,j,rflx(i,j),'rflx(i,j)'

#ifdef MASKING
            faces =  rmask(i-1,j)+rmask(i+1,j)+rmask(i,j-1)+rmask(i,j+1) !! amount of unmasked cells around
            if ( faces == 0 .or. rmask(i,j)>0  ) then
              error stop 'river grid position error'
            endif
            ! 10*ridx needed because uflx/vflx can be positive or negative around
!           the ridx number, and hence nearest integer is safest done with 10*.
            if (rmask(i-1,j)>0 ) riv_uflx(i  ,j) =-rfrc(i,j)/faces + 10*ridx(i,j)
            if (rmask(i+1,j)>0 ) riv_uflx(i+1,j) = rfrc(i,j)/faces + 10*ridx(i,j)
            if (rmask(i,j-1)>0 ) riv_vflx(i,j  ) =-rfrc(i,j)/faces + 10*ridx(i,j)
            if (rmask(i,j+1)>0 ) riv_vflx(i,j+1) = rfrc(i,j)/faces + 10*ridx(i,j)
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
