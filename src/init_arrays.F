#include "cppdefs.opt"

      subroutine init_arrays

      ! 'First touch' principle is automatically satisfied for a pure MPI code

      use private_scratch,   only: init_arrays_private_scratch
      use ocean_vars,        only: init_arrays_ocean
      use basic_output,      only: init_avg_arrays
      use coupling,          only: init_arrays_coupling
      use grid,              only: init_arrays_grid
      use eos_vars,          only: init_arrays_eos_vars
      use surf_flux,         only: init_arrays_surf_flx
      use mixing,            only: init_arrays_mixing
      use work_mod,          only: init_arrays_work_mod
      use boundary,          only: init_arrays_boundary
      use pipe_frc,          only: init_arrays_pipes
      use advection,         only: init_advection
#ifdef BULK_FRC
      use bulk_frc,          only: init_arrays_bulk_frc
#else
      use flux_frc,          only: init_arrays_flux_frc
#endif
#ifdef WEC
      use wec_frc,           only: init_arrays_wec_tile
      use dimensions,        only: nx, ny
#endif
#if defined(BIOLOGY_BEC2) || defined(MARBL)
      use bgc_shared_vars,   only: init_arrays_bgc_shared_vars
      use bgc_forces,        only: init_arrays_bgc_forces
      use bgc_io,            only: init_arrays_bgc_frc
#endif

      implicit none

      ! This whole thing could use some cleaning

      call init_arrays_private_scratch

      call init_arrays_ocean

      call init_advection

      call init_avg_arrays

      call init_arrays_coupling

      call init_arrays_grid

      call init_arrays_eos_vars

      call init_arrays_surf_flx

#ifdef BULK_FRC
      call init_arrays_bulk_frc
#else
      call init_arrays_flux_frc
#endif

#ifdef WEC
      ! Change this at a later stage
      call init_arrays_wec_tile(1,nx,1,ny)
#endif

! Set variable horizontal viscosities and tracer diffusion
! coefficients (see "mixing") to their background values.

      call init_arrays_mixing
      call init_arrays_work_mod   ! DevinD added this here to allocate
      call init_arrays_boundary

      call init_arrays_pipes

#if defined(BIOLOGY_BEC2) || defined(MARBL)
      call init_arrays_bgc_shared_vars
      call init_arrays_bgc_forces
      call init_arrays_bgc_frc
#endif

      end
