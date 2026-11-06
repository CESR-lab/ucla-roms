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
      use param
      use bulk_frc
      use flux_frc
      use surf_flux, only: set_surf_field_corr
      use tracers,   only: set_surf_tracer_flx
      use river_frc
      use pipe_frc
      use cdr_frc
      use wec_frc
      use bgc_io ! , only: set_bgc_surf_frc
      use roms_mpi
      use scalars
      use analytical
      use boundary

      implicit none

      ! Surface fluxes

# if defined QCORRECTION || defined SFLX_CORR
      call set_surf_field_corr
# endif

# ifdef BULK_FRC
      call set_bulk_frc
# endif
# ifdef ANA_SURFLUX
      call ana_smflux
      call ana_stflux
      call ana_srflux
# endif
# if !defined BULK_FRC && !defined ANA_SURFLUX
      call set_flux_frc
# endif

      if (nt>2) call set_surf_tracer_flx

      ! BGC surface flux:
# if defined(BIOLOGY_BEC2) || defined(MARBL)
      call set_bgc_surf_frc
# endif

      if (river_source) call set_river_frc

      if (pipe_source)  call set_pipe_frc

      if (cdr_source)   call set_cdr_frc

# ifdef WEC
      call set_wec_frc_tile(1,nx,1,ny)
      call wec_forces
# endif



      end
