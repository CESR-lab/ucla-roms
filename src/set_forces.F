#include "cppdefs.h"

      subroutine set_forces (tile)

! Using either data read from netCDF files or created analytically,
! prepare surface and bottom boundary fluxes as well as relevant
! climatological fields, so they can be applied to the model. This
! procedure essentially interpolates the available data to current
! model time and converts units to make all fluxes be kinematic
! fluxes, i.e,
!                      input data       conversion      kinematic
!                      units            factor          flux units
!
!  wind stress         [Newton/m^2]      1/rho0          [m^2/s^2]
!
!  heat, SWR fluxes    [Watts/m^2]       1/(rho*Cp)     [deg C * m/s]
!
!  fresh water flux     [cm/day]     S_surf*0.01/86400  [PSU *  m/s]
!
!     dQdSST       [Watts/(m^2*deg C)]   1/(rho*Cp)        [m/s]
!
! where S_surf is current model salinity near surface (i.e., fresh
! water precipitation/evaporation flux is converted into equivalent
! "salinity" flux. Units of climatological data: ssh[m], sst[deg C],
! tclima[deg C], and uclima[deg C] remain unchanged and only temporal
! interpolation has to be performed for them.
!
! Bottom drag is computed using either Styles and Glenn(1996) bottom
! boundary layer formulation, or linear/quadratic friction law..
!
! NOTE: Since some of the surface flux formulation may require
! climatological data (e.g., salinity at surface) update tclima
! first.

# if defined BULK_FRC
      use bulk_frc,  only: set_bulk_frc_tile
# else
      use flux_frc,  only: set_flux_frc
# endif
      use surf_flux, only: set_surf_field_corr
      use tracers,   only: set_surf_tracer_flx, iTandS

# ifdef RIVER_SOURCE
      use river_frc
# endif
# ifdef WEC
      use wec_frc
# endif

      implicit none


      real  tmid
      integer tile, ierr
#include "param.h"
#include "scalars.h"

#include "compute_tile_bounds.h"

      tmid=time+0.5*dt

      ierr=0
! External data to supply at open boundaries. Note that there are
! two mutually exclusive mechanisms for each variable: either _BRY
! or climatological data. 

#if defined T_FRC_BRY  || defined M2_FRC_BRY || \
    defined M3_FRC_BRY || defined Z_FRC_BRY
# ifdef ANA_BRY
c***              no code here
# else
      call set_bry_all_tile(istr,iend,jstr,jend, ierr)
# endif
#endif

!--> Climatological data

#if defined M2NUDGING && !defined M2_FRC_BRY
# ifdef ANA_SSH
      call ana_ssh_tile(istr,iend,jstr,jend)
# else
      call set_ssh_tile(istr,iend,jstr,jend, ierr)
# endif
#endif
#ifdef UCLIMATOLOGY
# ifdef ANA_UCLIMA
      call ana_uclima_tile(istr,iend,jstr,jend)
# else
      call set_uclima_tile (istr,iend,jstr,jend, ierr)
# endif
#endif
#if defined SOLVE3D && defined TCLIMATOLOGY
# ifdef ANA_TCLIMA
      call ana_tclima_tile(istr,iend,jstr,jend)
# else
      call set_tclima_tile(istr,iend,jstr,jend, ierr)
# endif
#endif

                                        !--> Surface fluxes


#ifdef ANA_SMFLUX
      call ana_smflux_tile(istr,iend,jstr,jend)
#endif

! Thermodynamic forcing: Note that BULK_FLUX requires computing the
! short-wave radiation flux first because bulk flux routine performs
! the final assembly of everything. Conversely if model is forced by
! precomputed total flux (which includes daily averaged short-wave 
! radiation interpolated in time), then to introduce DIURNAL CYCLE
! modulation set_srflux routine must interpolate short-wave flux in
! time first, then subtract it from total, then modulate short-wave,
! and, finally, add it back to total -- hence it must be called after.


#ifdef SOLVE3D

# if defined QCORRECTION || defined SFLX_CORR
      call set_surf_field_corr(istr,iend,jstr,jend)
# endif

# ifdef BULK_FRC
#  ifdef LMD_KPP
#   ifdef ANA_SRFLUX
      call ana_srflux_tile(istr,iend,jstr,jend) ! Should move this to bulk module
#   endif
#  endif
      call set_bulk_frc_tile(istr,iend,jstr,jend) ! bulk flux module

# else

      ! DevinD not sure what flag to use here to avoid set_flux for purely analytical
#  ifndef ANA_SMFLUX
      call set_flux_frc(istr,iend,jstr,jend)
#  endif

#  ifdef ANA_STFLUX
      call ana_stflux_tile(istr,iend,jstr,jend, itemp) ! DevinD Should move to module
#  endif
#  ifdef LMD_KPP
#   ifdef ANA_SRFLUX
      call ana_srflux_tile(istr,iend,jstr,jend) ! DevinD Should move to module
#   endif
#  endif
#  ifdef SALINITY
#   ifdef ANA_SSFLUX
      call ana_stflux_tile(istr,iend,jstr,jend, isalt) ! DevinD Should move to module
#   endif
#  endif
# endif

      ! DevinD, this would probably be faster if using CPPFLAG instead,
      ! e.g. PASSIVE_TRACERS for now just use number of tracers to decide:
      if(NT>iTandS) call set_surf_tracer_flx(istr,iend,jstr,jend)

# ifdef RIVER_SOURCE
      call set_river_frc_tile(istr,iend,jstr,jend)
#  endif

# ifdef WEC
      call set_wec_frc_tile(istr,iend,jstr,jend,ierr)
# endif

#endif  /* SOLVE3D */

                                 
!--> Bottom boundary fluxes [Styles and Glenn (1996) bottom
!    boundary layer formulation.  Not implemented in this code]

#if defined ANA_BMFLUX
      call ana_bmflux ILLEGAL
#elif defined SG_BBL96
# ifdef ANA_WWAVE
      call ana_wwave ILLEGAL
# else
      call set_wwave_tile(istr,iend,jstr,jend)
# endif
      call sg_bbl96 ILLEGAL
#endif

c>>>      if (ierr /= 0) may_day_flag=2
#ifdef ANA_PSOURCE
      if (ZEROTH_TILE) call ana_psource
#endif
      end
