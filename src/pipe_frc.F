      module pipe_frc

      ! Pipes forcing module
      ! --------------------

      ! for pipes:  Qbar(nx,y), the total pipe flux for that grid point ! (m3/s) function of time
      !             Qshape(npipes), the vertical shape of the pipe distribution ! fractions sum(Qshape) = 1

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)

#include "cppdefs.h"
#if defined PIPE_SOURCE
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use read_write

      implicit none

      private

      ! scalars.h needed for FIRST_TIME_STEP
#include "scalars.h"

      ! ****************************************************************
      ! user inputs
      integer,public,parameter :: npip = 1 ! number of pipes
#ifdef ANA_PIPE_FRC
      logical :: analytical = .true.      ! analytical forcing 
#else
      logical :: analytical = .false.     ! forcing netcdf files
#endif
      ! end user inputs
      ! ****************************************************************

      ! Variables used in the evolution equations
      integer,public  :: pidx
      integer,public, dimension(GLOBAL_2D_ARRAY) :: pipe_idx  ! pipe index
      real   ,public, dimension(GLOBAL_2D_ARRAY) :: pipe_ratio! pipe fractional flux

      real   ,public, dimension(GLOBAL_2D_ARRAY) :: pipe_flx  ! pipe flux
      real   ,public, dimension(npip,N)          :: pipe_prf  ! Pipe vertical profile
      real   ,public, dimension(npip,nt)         :: pipe_trc  ! Pipe tracer conc.

      real   ,public, dimension(npip)            :: pipe_vol  ! Pipe volume


      public set_pipe_frc_tile

      contains
! ----------------------------------------------------------------------

      subroutine set_pipe_frc_tile(istr,iend,jstr,jend)

      ! SET RIVER FORCES (REALISTIC OR ANALYTICAL FORCING):
      ! - read and interpolation all river forcing.
      ! - All river variables need time interpolation only
      !   here so can use same generic routine.
      ! - Input data in days!

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! Initialize pipe flux locations (One time only)
      if(FIRST_TIME_STEP) call init_pipe_frc(istr,iend,jstr,jend)

      ! Set river flux volumes and tracer data
      if(analytical) then
!       if(mynode==0) print*, 'rivers: going to set_ana_pip'
        call set_ana_pipe_frc
        
      else



      endif

      end subroutine set_pipe_frc_tile
! ----------------------------------------------------------------------

      subroutine init_pipe_frc(istr,iend,jstr,jend)

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

      else ! Read 'pipe_flux' from grid file (River mouth locations)
!        ierr=nf90_open(grdname, nf90_nowrite, ncid)
!        ierr=nf90_inq_varid(ncid, riv_flx_name, v_id)
!        ! 0 for only one record, 100 to prevent land masking in ncdf_read_mod as rivers
!        ! defined on land cells only.
!        ierr=ncdf_read_mod(ncid, v_id, 0, 100, rflx(START_2D_ARRAY), 1)
!        if(ierr/=0) call handle_ierr(ierr,'init_riv:: riv not in file!')

      endif

      ! Compute uflx,vflx from rflx
      call calc_pipe_flux(istr,iend,jstr,jend)

      if(mynode==0) print*, 'pipes: init pipe locations'
          
      end subroutine init_pipe_frc

! ----------------------------------------------------------------------
      subroutine calc_pipe_flux(istr,iend,jstr,jend)
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
!      do j = jstr,jend
!        do i = istr,iend
!          if (rflx(i,j) > 0) then ! distribute mass flux to all available unmasked cells
!            ! subtract 1e-5 in case only 1 grid point for river, so that floor still
!            ! produces correct iriver number.
!!            write(*,*) 'mynode=',mynode,'i,j',i,j,rflx(i,j),'rflx(i,j)'
!            iriver = floor(rflx(i,j)-1e-5)
!#ifdef MASKING
!            faces =  rmask(i-1,j)+rmask(i+1,j)+rmask(i,j-1)+rmask(i,j+1) !! amount of unmasked cells around
!            if ( faces == 0 .or. rmask(i,j)>0  ) then
!              error stop 'river grid position error'
!            endif
!            ! 10*iriver needed because uflx/vflx can be positive or negative around
!            ! the iriver number, and hence nearest integer is safest done with 10*.
!            if (rmask(i-1,j)>0 ) riv_uflx(i  ,j) =-(rflx(i,j)-iriver)/faces + 10*iriver
!            if (rmask(i+1,j)>0 ) riv_uflx(i+1,j) = (rflx(i,j)-iriver)/faces + 10*iriver
!            if (rmask(i,j-1)>0 ) riv_vflx(i,j  ) =-(rflx(i,j)-iriver)/faces + 10*iriver
!            if (rmask(i,j+1)>0 ) riv_vflx(i,j+1) = (rflx(i,j)-iriver)/faces + 10*iriver
!#endif
!          endif
!        enddo
!      enddo

      end subroutine calc_pipe_flux
! ----------------------------------------------------------------------

# endif /* ifdef PIPE_SOURCE */
      end module pipe_frc
