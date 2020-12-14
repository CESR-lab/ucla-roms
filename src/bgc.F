      module bgc

      ! BGC - bio-geo-chemical module
      ! -----------------------------

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Dec)
      ! Refactoring of ETH Zurich roms code with BEC2, code which has been used
      ! by Pierre Damien to run full pacific models.

#include "cppdefs.h"
#ifdef BIOLOGY_BEC2

      implicit none

      private ! Make all variable private to this module unless public specified

      ! Includes:
      ! param.h needed for GLOBAL_2D_array to work. NT = number tracer from param.h
#include "param.h"

      ! ****************************************************************
      ! User inputs

      integer, parameter :: nbgc_flx = 2 ! Number of surface bgc flux forcings
      ! Interpolate forcing from coarser input grid (=1) or not (=0)
      integer :: interp_frc = 0 ! Factor 2 only for now
      ! dust forcing
      character(len=4) :: dust_vname = 'dust'
      character(len=9) :: dust_tname = 'dust_time'
      ! iron forcing
      character(len=4) :: iron_vname = 'iron'
      character(len=9) :: iron_tname = 'iron_time'

      ! End of user inputs
      ! *************************************************************

      ! Netcdf variables for surface forcing
      real                         :: flx_data(GLOBAL_2D_ARRAY,2,nbgc_flx)
      integer, dimension(nbgc_flx) :: flx_grd_type  = 0
      integer, dimension(nbgc_flx) :: flx_file_indx = 0
      integer, dimension(nbgc_flx) :: flx_irec      = 0
      integer, dimension(nbgc_flx) :: flx_it1 = 1, flx_it2 = 2
      real, dimension(2,nbgc_flx)  :: flx_times = -99 ! [-99,-99]


      public set_bgc_surf_frc

      contains
! ----------------------------------------------------------------------

      subroutine set_bgc_surf_frc(istr,iend,jstr,jend)
      ! read in bgc surface flux and interpolate to model time.
      ! Taken from get_smth.F & set_smth.F of ETH code.

      ! Since bgc is set up such that iron and dust have flux at surface but
      ! are not directly bgc tracers, I have separated them into set_bgc_frc

      use read_write, only: set_frc_var_tile

      implicit none

! Needed for iron and dust variables, not yet ported to this module. Should do.
#include "bgc_forces.h"

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! dust
      call set_frc_var_tile(
     &        istr,iend,          jstr,jend,        ! Tile bounds indices
     &        dust_vname,         dust_tname,       ! Text name of array and time
     &        flx_data(:,:,:,1),  dust,             ! Variable arrays
     &        flx_times(:,1),     flx_grd_type(1),  ! Input times & Grid type
     &        flx_file_indx(1),   flx_irec(1),      ! File indx & Current input record
     &        flx_it1(1),         flx_it2(1),       ! Time index placeholders
     &        interp_frc )                          ! Online spatial interpolation flag

      ! iron
      call set_frc_var_tile(
     &        istr,iend,          jstr,jend,        ! Tile bounds indices
     &        iron_vname,         iron_tname,       ! Text name of array and time
     &        flx_data(:,:,:,2),  iron,             ! Variable arrays
     &        flx_times(:,2),     flx_grd_type(2),  ! Input times & Grid type
     &        flx_file_indx(2),   flx_irec(2),      ! File indx & Current input record
     &        flx_it1(2),         flx_it2(2),       ! Time index placeholders
     &        interp_frc )                          ! Online spatial interpolation flag

!      if(mynode==0) print*, 'dust(20,20) ',dust(20,20)
!      if(mynode==0) print*, 'iron(20,20) ',iron(20,20)

      end subroutine set_bgc_surf_frc
! ----------------------------------------------------------------------

#endif /* BIOLOGY_BEC2 - entire module */
      end module bgc
