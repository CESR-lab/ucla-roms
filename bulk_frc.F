      module bulk_frc

      ! ===================
      ! BULK FORCING MODULE
      ! ===================

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.h"
! CPP flag BULK_FRC_MOD covering entire module so that the module is empty if
! BULK_FRC_MOD is unused, and thus module will not take up unnecessary size in
! the executable.
#if defined BULK_FRC_DPD /* for whole module! */

      ! Modules required
      ! ----------------

      use read_write, only: set_frc_var_tile

      implicit none

      ! Make all variable private to this module unless public specified
      private

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_array to work
#include "param.h"

      ! Variables
      ! =========

      ! wind
      ! ----

       ! u-wind
       ! ------

        ! stores raw data u-wind for 2 times: time_A < model_time < time_B
      real uwnd_data(GLOBAL_2D_ARRAY,2)
      real uwnd(GLOBAL_2D_ARRAY) ! Stores time interpolated u-wind
      ! grid type are needed for ncdf_read in read_write ncdf module
      integer :: uwnd_grd_type = 0 ! should be = r_var but need to include ncvars.h
      character(len=4) :: uwnd_name = 'uwnd'
      ! timestep variables for wind bulk force
      ! Notice time here is called something different wnd not uwnd.
      ! get_smth.F has a method to deal with different time names, might
      ! be worth using. Hard coded now for ease.
      character(len=8) :: uwnd_time_name = 'wnd_time'
      integer :: uwnd_file_indx = 0 ! Set file index in force file list as zero to begin with
      integer :: uwnd_irec      = 0 ! Set netcdf record number to zero to start
      ! These were only one variable (itwnd) in old code, but maybe it is clearer
      ! with 2. Although if it's one it's forcing you to realise that
      ! it jumps between values.
      integer :: uwnd_it1 = 1, uwnd_it2 = 2  ! Placeholder for read in time
      real, dimension(2) :: uwnd_times = [-99,-99] ! Store input file times

       ! v-wind
       ! ------
      real    :: vwnd_data(GLOBAL_2D_ARRAY,2)
      real    :: vwnd(GLOBAL_2D_ARRAY)
      integer :: vwnd_grd_type  = 0
      integer :: vwnd_file_indx = 0
      integer :: vwnd_irec      = 0
      integer :: vwnd_it1 = 1, vwnd_it2 = 2
      real, dimension(2) :: vwnd_times = [-99,-99]
      character(len=4)   :: vwnd_name = 'vwnd'
      character(len=8)   :: vwnd_time_name = 'wnd_time'

      ! Atmospheric tracers
      ! -------------------


      public set_bulk_frc_tile

      contains


! ----------------------------------------------------------------------

      subroutine set_bulk_frc_tile(istr,iend,jstr,jend)

      ! ===================
      ! SET ALL BULK FORCES
      ! ===================

      ! Landing subroutine from set_forces
      ! where all bulk forces are set.

      ! All bulk variables need time interpolation only
      ! here so can use same generic routine

      implicit none

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend

      ! Wind

        ! u-wind (time in days)

          ! 1) Find & save force variable file_indx & correct record (irec)
          ! 2) Read in forcing data
          ! 3) Interpolate input forcing to correct model time

      ! Could use arrays for bulk forces only then don't need to pass
      ! through so many variables.
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          uwnd_name,      uwnd_time_name, ! Text names
     &          uwnd_data,      uwnd,           ! Variable arrays
     &          uwnd_times,     uwnd_grd_type,  ! Input times & Grid type
     &          uwnd_file_indx, uwnd_irec,      ! File indx & Current input record
     &          uwnd_it1,       uwnd_it2 )      ! Time index placeholders

!      mpi_master_only print*, 'uwnd(1,1)',uwnd(1,1)

        ! v-wind

      ! Atmospheric tracers  - TIME INTERP ONLY
        ! get_bulk_tra in old code

        ! tair (air temp)

        ! qair (specific humidity) rhum in new code!

      ! prec -> precipitations - TIME INTERP ONLY (scale but =1)

        ! rain

      ! rad -> short/long radiation - TIME INTERP ONLY (scale but = 1)

        ! lwrad - Long wave radiation

        ! swrad - short wave radiot

      end subroutine set_bulk_frc_tile


!-----------------------------------------------------------------------

#endif /* BULK_FRC_MOD for whole module */

      end module bulk_frc
