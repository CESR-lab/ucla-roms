      module analytical

      ! all analytical routines contained in this module

#include "cppdefs.opt"

      use param
      use dimensions
      use ocean_vars
      use grid
      use scalars
      use tracers
      use mixing
      use scalars
      use surf_flux

      implicit none
      private

! ANALYTICAL PACKAGE:   Set up various analytical forcing fields
!=========== ========   for the model.
!
! ana_smflux        Analytical kinematic surface momentum flux (wind stress).
! ana_srflux        Analytical kinematic surface shortwave radiation.
! ana_sst           Analytical sea surface temperature and dQdSST
!                       which are used during heat flux correction.
! ana_stflux        Analytical kinematic surface flux of tracer type variables.

      public ana_init
      public ana_vmix
      public ana_smflux
      public ana_srflux
      public ana_stflux
      public ana_init_generic

      contains

! ----------------------------------------------------------------------
      subroutine ana_init ![

      use boundary
      use eos_vars

      implicit none

#include "ana_init.h"

      call ana_init_generic

      end subroutine ana_init !]
! ----------------------------------------------------------------------
      subroutine ana_init_generic ![

      use param
      use boundary
      implicit none

      ! Local
      integer :: i,j,k, itrc

! The rest of the code is independent of the specific configuration.
! It fills analytical boundary values and exchanges u,v and tracers to 
! fill halos

! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.


#  ifdef ANA_BRY
      if (obc_west.and.west_bnd) then
#    ifdef M2_FRC_BRY
        do j=1,ny                       ! Save perimeter values
          zeta_west(j)=zeta(0,j,1)      ! of the initial fields
          ubar_west(j)=ubar(1,j,1)      ! to be used as external
          vbar_west(j)=vbar(0,j,1)      ! forcing data for side
        enddo                           ! boundaries.
#    endif
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
      endif

      if (obc_east.and.east_bnd) then
#    ifdef M2_FRC_BRY
        do j=1,ny
          zeta_east(j)=zeta(nx+1,j,1)
          ubar_east(j)=ubar(nx+1,j,1)
          vbar_east(j)=vbar(nx+1,j,1)
        enddo
#    endif
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
      endif


      if (obc_south.and.south_bnd) then
#    ifdef M2_FRC_BRY
        do i=1,nx
          zeta_south(i)=zeta(i,0,1)
          ubar_south(i)=ubar(i,0,1)
          vbar_south(i)=vbar(i,1,1)
        enddo
#    endif
        do k=1,nz
          do i=1,nx
#    ifdef M3_FRC_BRY
            u_south(i,k)=u(i,0,k,1)
            v_south(i,k)=v(i,1,k,1)
#    endif
#    ifdef T_FRC_BRY
            t_south(i,k,:)=t(i,0,k,1,:)
#    endif
          enddo
        enddo
      endif

      if (obc_north.and.north_bnd) then
#    ifdef M2_FRC_BRY
        do i=1,nx
          zeta_north(i)=zeta(i,ny+1,1)
          ubar_north(i)=ubar(i,ny+1,1)
          vbar_north(i)=vbar(i,ny+1,1)
        enddo
#    endif
        do k=1,nz
          do i=1,nx
#    ifdef M3_FRC_BRY
            u_north(i,k)=u(i,ny+1,k,1)
            v_north(i,k)=v(i,ny+1,k,1)
#    endif
#    ifdef T_FRC_BRY
            t_north(i,k,:)=t(i,ny+1,k,1,:)
#    endif
          enddo
        enddo
      endif
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
      call exchange_xxx(Akv,Akt)
# endif  /* EXCHANGE */


      end subroutine ana_vmix !]
! ----------------------------------------------------------------------
      subroutine ana_smflux ![

      ! Analytical surface momentum stress
      ! provide ana_smflux.h for non-trivial forcing


      implicit none
      integer i,j

      if (iic==ntstart) then     ! Set wind stress (kinematic surface
        do j=1-bf,ny+bf          ! momentum flux [m^2/s^2]) components
          do i=1-bf,nx+bf        ! in x- and y-directions [m^2/s^2])
            sustr(i,j)=0.0       ! at horizontal U- and V-points.
            svstr(i,j)=0.0       ! For code compactness, both are set
          enddo                  ! to zero at first (default values),
        enddo                    ! then either one of them, or both may
      endif                      ! be overwritten by nontrivial values.

      end subroutine ana_smflux !]
! ----------------------------------------------------------------------
      subroutine ana_srflux ![

! Set kinematic surface solar shortwave radiation flux "srflx"
! [degC m/s] using an analytical expression.

      implicit none
      integer :: i,j

      do j=1-bf,ny+bf
        do i=1-bf,nx+bf
          srflx(i,j)=0.
        enddo
      enddo
      end subroutine ana_srflux !]
! ----------------------------------------------------------------------
      subroutine ana_stflux  ![

! Set kinematic surface flux of tracer type variables "stflx" (tracer
! units m/s) using analytical expressions.

      implicit none
      integer i,j

! Set kinematic surface heat flux [degC m/s] at horizontal
! RHO-points.

        do j=1-bf,ny+bf
          do i=1-bf,nx+bf
            stflx(i,j,itemp)=0.
            stflx(i,j,isalt)=0.
          enddo
        enddo

      end subroutine ana_stflux !]
! ---------------------------------------------------------------------

      end module analytical
