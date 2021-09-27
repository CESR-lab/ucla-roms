#include "cppdefs.h"

! ANALYTICAL PACKAGE:   Set up various analytical forcing fields
!=========== ========   for the model.
!
! ana_bsedim_tile   Analytical bottom sediment grain size
!                          and density.
! ana_meanRHO_tile  Analytical mean density anomaly.
! ana_smflux_tile   Analytical kinematic surface momentum flux
!                          (wind stress).
! ana_srflux_tile   Analytical kinematic surface shortwave
!                          radiation.
! ana_ssh_tile      Analytical sea surface height climatology.
! ana_sst_tile      Analytical sea surface temperature and dQdSST
!                         which are used during heat flux correction.
! ana_stflux_tile   Analytical kinematic surface flux of tracer type
!                          variables.
! ana_tclima_tile   Analytical tracer climatology fields.
! ana_uclima_tile   Analytical tracer climatology fields.
! ana_wwave_tile    Analytical wind induced wave amplitude,
!                         direction and period.


#ifdef ANA_SMFLUX
      subroutine ana_smflux_tile(istr,iend,jstr,jend)

      use param
      use surf_flux, only: sustr,svstr
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
      real TauX, TauY, cff

# include "compute_extended_bounds.h"

      if (FIRST_TIME_STEP) then  ! Set wind stress (kinematic surface
        do j=jstrR,jendR         ! momentum flux [m^2/s^2]) components
          do i=istrR,iendR       ! in XI- and ETA-directions [m^2/s^2])
            sustr(i,j)=0.        ! at horizontal U- and V-points.
            svstr(i,j)=0.        ! For code compactness, both are set
          enddo                  ! to zero at first (default values),
        enddo                    ! then either one of them, or both may
      endif                      ! be overwritten by nontrivial values.

      end
#endif /* ANA_SMFLUX */

#if defined ANA_SRFLUX && defined SOLVE3D

      subroutine ana_srflux_tile(istr,iend,jstr,jend)

! Set kinematic surface solar shortwave radiation flux "srflx"
! [degC m/s] using an analytical expression.
      use param
      use surf_flux, only: srflx
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j

# include "compute_extended_bounds.h"

      do j=jstrR,jendR
        do i=istrR,iendR
          srflx(i,j)=0.
        enddo
      enddo
      end
#endif /* ANA_SRFLUX */


#if defined ANA_SSH
      subroutine ana_ssh(tile)            ! Set analytical sea-surface
      use param
      use hidden_mpi_vars
      implicit none                       ! height climatology [m]
      integer tile
#include "compute_tile_bounds.h"
      call ana_ssh_tile(istr,iend,jstr,jend)
      end

      subroutine ana_ssh_tile(istr,iend,jstr,jend)

      use param
      use climat
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
#ifdef JET
      real TL0, TC0, cff1, cff2, cff3, cff4
#endif

# include "compute_extended_bounds.h"

      end
#endif



#ifdef SOLVE3D
# if defined ANA_STFLUX || defined ANA_SSFLUX
      subroutine ana_stflux_tile(istr,iend,jstr,jend, itrc)

! Set kinematic surface flux of tracer type variables "stflx" (tracer
! units m/s) using analytical expressions.
      use param
      use surf_flux, only: stflx
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, itrc, i,j

# include "compute_extended_bounds.h"

      if (itrc==itemp) then

! Set kinematic surface heat flux [degC m/s] at horizontal
! RHO-points.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,itemp)=0.
          enddo
        enddo

#  ifdef SALINITY
      elseif (itrc == isalt) then

! Set kinematic surface freshwater flux (m/s) at horizontal
! RHO-points, scaling by surface salinity is done in STEP3D.

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,isalt)=0.
          enddo
        enddo
#  endif
      else

!  Set kinematic surface flux of additional tracers, if any.

      endif
      end
# endif /* ANA_STFLUX || ANA_SSFLUX */

#endif /* SOLVE3D */

      subroutine ana_empty
      end

! ----------------------------------------------------------------------
# if defined RIVER_SOURCE
      subroutine set_ana_river_frc
      ! Analytical river forcing volume and tracer data

      use river_frc
      implicit none

      ! Set river volume and tracer values for each time step
      riv_vol(1)   = 5e2  ! Volume flux in m3/s
      riv_trc(1,1) = 24.0 ! Temperature in  Degrees C
      riv_trc(1,2) =  1.0  ! Salinity in PSU

      end subroutine set_ana_river_frc
#endif /* RIVER_SOURCE */


