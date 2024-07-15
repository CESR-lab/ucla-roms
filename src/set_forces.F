#include "cppdefs.opt"

      subroutine set_forces

! Using either data read from netCDF files or created analytically,
! prepare surface and bottom boundary fluxes, so they can be 
! applied to the model. This
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
      use param
      use bulk_frc
      use flux_frc
      use surf_flux, only: set_surf_field_corr
      use tracers,   only: set_surf_tracer_flx, iTandS
      use river_frc
      use pipe_frc
      use tides
      use wec_frc
      use bgc ! , only: set_bgc_surf_frc
      use hidden_mpi_vars
      use scalars
      use analytical
      use boundary

      implicit none

      
      integer ierr
      integer,save:: tile=0

#include "compute_tile_bounds.h"

      ierr=0
! External data to supply at open boundaries. Note that there are
! two mutually exclusive mechanisms for each variable: either _BRY

#if  !defined OBC_NONE &&                         \
     (defined T_FRC_BRY  || defined M2_FRC_BRY || \
      defined M3_FRC_BRY || defined Z_FRC_BRY )
# ifdef ANA_BRY
c***              no code here
# else
!     call set_bry_all_tile(istr,iend,jstr,jend, ierr)
!     call set_bry_all
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
      call set_surf_field_corr
# endif

# ifdef BULK_FRC
#  ifdef LMD_KPP
#   ifdef ANA_SRFLUX
      call ana_srflux_tile(istr,iend,jstr,jend) ! Should move this to bulk module
#   endif
#  endif
      call set_bulk_frc(istr,iend,jstr,jend)
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

      if(NT>iTandS) call set_surf_tracer_flx(istr,iend,jstr,jend)

      ! BGC surface flux:
# if defined(BIOLOGY_BEC2) || defined(MARBL)
      call set_bgc_surf_frc(istr,iend,jstr,jend)
# endif


      if (river_source) call set_river_frc

      if (pipe_source) call set_pipe_frc

# ifdef WEC
      call set_wec_frc_tile(istr,iend,jstr,jend,ierr)
      call wec_forces
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
