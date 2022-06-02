#include "cppdefs.opt"

      subroutine get_forces

! Read various forcing and climatological fields from their netCDF
! files and save them as globally accessible arrays in common blocks
! "forces.h" and "climat.F". NOTE that NO EFFORT has been made to
! convert units of these fluxes at this stage (i.e., convert them
! into kinematic fluxes): typically data is assumed to have dynamical
! physical units, thus wind stress is in [Newton/m^2]; surface heat
! and short wave radiation fluxes are in [Watts/m^2]; SST sensitivity
! to heat flux dQdSST is in [Watts/(deg C * m^2)]; surface fresh
! water flux is in [mm/day].

      use param
      use scalars
      implicit none
      integer ierr

c--#define VERBOSE
#ifdef VERBOSE
      write(*,*) 'entering get_forces' MYID
#endif
      ierr=0
#if defined M2NUDGING && !defined M2_FRC_BRY && !defined ANA_SSH
      call get_ssh(ierr)
#endif
#if defined UCLIMATOLOGY && !defined ANA_UCLIMA
      call get_uclima(ierr)
#endif
#ifdef SOLVE3D
# if ( defined TCLIMATOLOGY ||( defined TNUDGING && \
       defined T_FRC_BRY )) && !defined ANA_TCLIMA
      call get_tclima(ierr)
# endif
#endif

      if (ierr /= 0) may_day_flag=2
#ifdef VERBOSE
      write(*,*) 'leaving get_forces' MYID
#endif
      end

