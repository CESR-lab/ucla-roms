      module river_frc

      ! RIVER FORCING MODULE
      ! --------------------

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)

#include "cppdefs.h"
#if defined RIVER_SOURCE

      ! Contains all the necessary components to produce the riv_uflx,riv_vflx
      ! arrays which have the the horizontal river volume flux in m2/s at the current time
      ! riv_uflx, riv_vflx should be on faces between a masked and unmasked cell,
      ! in the direction of the unmasked cell.
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use read_write

      implicit none

      private

      ! scalars.h needed for FIRST_TIME_STEP
#include "scalars.h"
#include "river_frc.opt"

      ! Variables used for equation system calculations
      real,public, dimension(GLOBAL_2D_ARRAY) :: riv_uflx
      real,public, dimension(GLOBAL_2D_ARRAY) :: riv_vflx
      real,public, dimension(nriv)    :: riv_vol      ! River volume
      real,public, dimension(nriv,nt) :: riv_trc      ! River tracer conc.
      real,public, dimension(GLOBAL_2D_ARRAY) :: rflx ! River locations

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

      public set_river_frc_tile

      contains
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

      ! Initialize river flux locations (One time only)
      if(FIRST_TIME_STEP) call init_river_frc(istr,iend,jstr,jend)

      ! Set river flux volumes and tracer data
      if(analytical) then
!       if(mynode==0) print*, 'rivers: going to set_ana'
#ifdef ANA_RIVER_FRC
        call set_ana_river_frc ! cppflags needed else won't link without the analytical.F
#endif

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

!        if(mynode==0) print*, 'riv_trc(1,1)= ', riv_trc(1,1) ! debug
!        if(mynode==0) print*, 'riv_trc(1,2)= ', riv_trc(1,2) ! debug

      endif

      end subroutine set_river_frc_tile
! ----------------------------------------------------------------------

      subroutine init_river_frc(istr,iend,jstr,jend)

      ! Initialize river forcing:
      ! Read in a grid file with locations of river mouths and flux contribution per cell.
      ! Done only once as river mouth position does not change.

      use netcdf
      implicit none

      ! inputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer :: ierr,ncid,v_id

      ! Should put this in wrapper in read_write.F as 'nc_read_once'
      if (analytical) then

        ! rflx is defined in ana_grid

      else ! Read 'river_flux' from grid file (River mouth locations)
        ierr=nf90_open(grdname, nf90_nowrite, ncid)
        ierr=nf90_inq_varid(ncid, riv_flx_name, v_id)
        ! 0 for only one record, 100 to prevent land masking in ncdf_read_mod as rivers
        ! defined on land cells only.
        ierr=ncdf_read_mod(ncid, v_id, 0, 100, rflx(START_2D_ARRAY), 1)
        if(ierr/=0) call handle_ierr(ierr,'init_riv:: riv not in file!')

      endif

      ! Compute uflx,vflx from rflx
      call calc_river_flux(istr,iend,jstr,jend)

      if(mynode==0) print*, 'rivers: init river locations'
          
      end subroutine init_river_frc

! ----------------------------------------------------------------------
      subroutine calc_river_flux(istr,iend,jstr,jend)
      ! Calculate the river flux contributions to each cell.
      ! river_flux = iriver + fraction of river's flux through grid point.
      ! E.g. River 3 is over 2 grid points (half flux through each point),
      ! hence river_flux = 3 + 0.5 = 3.5
      implicit none

      ! inputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! local
      integer :: i,j,faces

      ! compute uflx,vflx from rflx
      do j = jstr-1,jend+1   ! Loop over -1 and +1 because rflx cell only flows into
        do i = istr-1,iend+1 ! neighbour, hence cell next to boundary could flow into cell.
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

# endif /* ifdef RIVER_SOURCE */
      end module river_frc
