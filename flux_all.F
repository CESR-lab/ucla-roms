      module flux_all

      ! ===================
      ! FLUX FORCING MODULE
      ! ===================

      ! For generic flux functionality only.
      ! See bulk_frc.F, flux & analytical modules for specific fluxes

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.h"
! CPP flag XXXX covering entire module so that the module is empty if
! BULK_FRC_MOD is unused, and thus module will not take up unnecessary size in
! the executable.
! Currently only set-up for BULK_FLUX but eventually should be used for all
! fluxes once forces.h is fully ported into flux modules.
#if defined BULK_FLUX /* for whole module! */

      ! Modules required
      ! ================

      use read_write, only: set_frc_var_tile

      implicit none

      ! Make all variable private to this module unless public specified
      private

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_array to work. NT = number tracer from param.h
#include "param.h"

      ! Variables
      ! =========

      ! Total Flux variables (used by flux, bulk & analytical forcing)
      ! --------------------

      ! The following is taken from forces.h from the new code:

      ! Surface momentum flux (wind stress):
      !======== ======== ==== ===== ========
      ! sustr   XI- and ETA-components of kinematic surface momentum
      ! svstr   flux (wind stresses) at defined horizontal U- and
      !         V-points,  dimensioned as [m^2/s^2].

      ! uwind  two-time level gridded data for XI- anf ETA-componets
      ! vwind  of wind stess (normally assumed to be in [Newton/m^2].

      real, public :: sustr(GLOBAL_2D_ARRAY)
      real, public :: svstr(GLOBAL_2D_ARRAY)

      ! uwind & vwind of new code is currently only in bulk_frc.F
      ! renamed as uwind_data & vwind_data
!#ifdef WIND_MAGN
!      real wndmag(GLOBAL_2D_ARRAY)
!#endif

      ! Currently haven't incorporated dQdSST from forces.h:
      ! Data for net heat flux sensitivity to sea-surface temperature,
      ! dQdSST = d NetHeatFlux/d SST [Watts/m^2/Celsius, always negative]
      ! See Eq.(6) from B. Barnier, L. Siefridt, and P. Marchesiello, 1995.

      ! ----------------------------------------------------------------

      ! The following is taken from forces.h from the new code:

      ! Surface tracer fluxes (used by flux, bulk & analytical forcing)
      ! ---------------------
      !  stflx   kinematic surface fluxes for tracer-type variables at
      !          horizontal RHO-points. Physical dimensions [degC m/s]
      !          - temperature; [PSU m/s] - salinity.
      real, public :: stflx(GLOBAL_2D_ARRAY,NT) ! NT = number tracer from param.h

        ! sea-surface temperature (SST) data
        ! -----------------------
      real, public :: sst_data(GLOBAL_2D_ARRAY,2)
      ! The actual variable lives in stflx(i,j,itemp)
      integer :: sst_grd_type  = 0
      integer :: sst_file_indx = 0
      integer :: sst_irec      = 0
      integer :: sst_it1 = 1, sst_it2 = 2
      real, dimension(2) :: sst_times = [-99,-99]
      character(len=5)   :: sst_name = 'shflux' ! Name taken from get_stflux.F
      character(len=8)   :: sst_time_name = 'shf_time'

        ! Sea-surface salinity (SSS) data
        ! --------------------
      real, public :: sss_data(GLOBAL_2D_ARRAY,2)
      ! The actual variable lives in stflx(i,j,isalt)
      integer :: sss_grd_type  = 0
      integer :: sss_file_indx = 0
      integer :: sss_irec      = 0
      integer :: sss_it1 = 1, sss_it2 = 2
      real, dimension(2) :: sss_times = [-99,-99]
      character(len=5)   :: sss_name = 'swflux' ! Name taken from get_stflux.F
      character(len=8)   :: sss_time_name = 'swf_time'


      contains


! ----------------------------------------------------------------------

      subroutine set_flux_correctors(istr,iend,jstr,jend)

      ! ===============================
      ! READ IN SURFACE FLUX CORRECTORS
      ! FOR TEMPERATURE & SALINITY
      ! ==========================

      implicit none

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend

#if !defined ANA_STFLUX
      ! Sea-surface temperature (SST) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sst_name,       sst_time_name,  ! Text names
     &          sst_data,                       ! Variable input data arrays
     &          stflx(GLOBAL_2D_ARRAY,itemp),   ! Variable arrays
     &          sst_times,      sst_grd_type,   ! Input times & Grid type
     &          sst_file_indx,  sst_irec,       ! File indx & Current input record
     &          sst_it1,        sst_it2 )       ! Time index placeholders
#endif

#if !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sss_name,       sss_time_name,  ! Text names
     &          sss_data,                       ! Variable input data arrays
     &          stflx(GLOBAL_2D_ARRAY,isalt),   ! Variable arrays
     &          sss_times,      sss_grd_type,   ! Input times & Grid type
     &          sss_file_indx,  sss_irec,       ! File indx & Current input record
     &          sss_it1,        sss_it2 )       ! Time index placeholders
#endif

      end subroutine set_flux_correctors


!-----------------------------------------------------------------------

#endif /* BULK_FLUX for whole module */

      end module flux_all
