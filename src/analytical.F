      module analytical

      ! all analytical routines contained in this module

#include "cppdefs.opt"

      use param
      use ocean_vars
      use grid
      use scalars
      use tracers
      use mixing
      use scalars

      implicit none
      private

! ANALYTICAL PACKAGE:   Set up various analytical forcing fields
!=========== ========   for the model.
!
! ana_meanRHO_tile  Analytical mean density anomaly.
! ana_smflux_tile   Analytical kinematic surface momentum flux (wind stress).
! ana_srflux_tile   Analytical kinematic surface shortwave radiation.
! ana_sst_tile      Analytical sea surface temperature and dQdSST
!                         which are used during heat flux correction.
! ana_stflux_tile   Analytical kinematic surface flux of tracer type variables.
! ana_wwave_tile    Analytical wind induced wave amplitude, direction and period.

      public ana_init
      public ana_vmix
      public ana_smflux_tile
      public ana_srflux_tile
      public ana_stflux_tile
      public ana_init_generic

      contains

! ----------------------------------------------------------------------
      subroutine ana_init ![

      use boundary
      use eos_vars

      implicit none

#include "ana_init.h"

      call ana_init_generic(1,nx,1,ny)

      end subroutine ana_init !]
! ----------------------------------------------------------------------
      subroutine ana_init_generic(istr,iend,jstr,jend) ![

      use param
      use boundary
      implicit none

      ! input
      integer, intent(in) :: istr,iend,jstr,jend
      ! local
      integer :: i,j,k, itrc

! The rest of the code is independent of the specific configuration.

! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.


#  ifdef ANA_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=1,ny                       ! Save perimeter values
          zeta_west(j)=zeta(0,j,1)      ! of the initial fields
          ubar_west(j)=ubar(1,j,1)      ! to be used as external
          vbar_west(j)=vbar(0,j,1)      ! forcing data for side
        enddo                           ! boundaries.
#    endif
#    ifdef SOLVE3D
        do k=1,nz
          do j=1,ny
#     ifdef M3_FRC_BRY
            u_west(j,k)=u(1,j,k,1)
            v_west(j,k)=v(0,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_west(j,k,:)=t(0,j,k,1,:)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_EAST
      if (EASTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=1,ny
          zeta_east(j)=zeta(nx+1,j,1)
          ubar_east(j)=ubar(nx+1,j,1)
          vbar_east(j)=vbar(nx+1,j,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,nz
          do j=1,ny
#     ifdef M3_FRC_BRY
            u_east(j,k)=u(nx+1,j,k,1)
            v_east(j,k)=v(nx+1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_east(j,k,:)=t(nx+1,j,k,1,:)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif


#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=1,nx
          zeta_south(i)=zeta(i,0,1)
          ubar_south(i)=ubar(i,0,1)
          vbar_south(i)=vbar(i,1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,nz
          do i=1,nx
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,0,k,1)
            v_south(i,k)=v(i,1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,:)=t(i,0,k,1,:)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=1,nx
          zeta_north(i)=zeta(i,ny+1,1)
          ubar_north(i)=ubar(i,ny+1,1)
          vbar_north(i)=vbar(i,ny+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,nz
          do i=1,nx
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,ny+1,k,1)
            v_north(i,k)=v(i,ny+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,:)=t(i,ny+1,k,1,:)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  endif /* ANA_BRY */


# ifdef EXCHANGE
      call exchange_xxx(zeta(:,:,1),ubar(:,:,1),vbar(:,:,1))
#  ifdef SOLVE3D
      call exchange_xxx(u(:,:,:,1),v(:,:,:,1))
      do itrc=1,nt
        call exchange_xxx(t(:,:,:,1,itrc) )
      enddo
#  endif
# endif  /* EXCHANGE */


      end subroutine ana_init_generic !]
! ---------------------------------------------------------------------
      subroutine ana_vmix ![

! Set vertical kinematic vertical viscosity "Akv" [m^2/s] for momentum
! and vertical mixing c efficients "Akt" [m^2/s] for tracer equations
! using analytical expressions.

      implicit none

#include "ana_vmix.h"

# ifdef EXCHANGE
      call exchange_xxx(Akv,Akt(:,:,:,itemp))

#  ifdef SALINITY
      call exchange_xxx(Akt(:,:,:,isalt))
#  endif

# endif  /* EXCHANGE */


      end subroutine ana_vmix !]
! ----------------------------------------------------------------------
      subroutine ana_smflux_tile(istr,iend,jstr,jend) ![

      use param
      use surf_flux !, only: sustr,svstr
      use grid
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
      real TauX, TauY, cff

# include "compute_extended_bounds.h"

      if (FIRST_TIME_STEP) then  ! Set wind stress (kinematic surface
        do j=jstrR,jendR         ! momentum flux [m^2/s^2]) components
          do i=istrR,iendR       ! in XI- and ETA-directions [m^2/s^2])
            sustr(i,j)=0.0       ! at horizontal U- and V-points.
            svstr(i,j)=0.0       ! For code compactness, both are set
          enddo                  ! to zero at first (default values),
        enddo                    ! then either one of them, or both may
      endif                      ! be overwritten by nontrivial values.

      end subroutine ana_smflux_tile !]
! ----------------------------------------------------------------------
      subroutine ana_srflux_tile(istr,iend,jstr,jend)  ![

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
      end subroutine ana_srflux_tile  !]
! ----------------------------------------------------------------------
      subroutine ana_stflux_tile(istr,iend,jstr,jend, itrc)  ![

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
      end subroutine ana_stflux_tile  !]
! ---------------------------------------------------------------------

      end module analytical
