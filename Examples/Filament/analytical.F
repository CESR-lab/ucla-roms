      module analytical

      ! all analytical routines contained in this module

#include "cppdefs.opt"

      use param
      use ocean_vars
      use grid
      use scalars
      use tracers
      use mixing

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

      ! user input:
      logical, parameter :: restore = .false.

! REMOVE ALL FLAGS LATER
#ifdef ANA_INITIAL
      public ana_init
#endif
#ifdef ANA_VMIX
      public ana_vmix
#endif

#ifdef ANA_SMFLUX
      public ana_smflux_tile
#endif
#if defined ANA_SRFLUX && defined SOLVE3D
      public ana_srflux_tile
#endif
#if defined ANA_SSH
      public ana_ssh_tile
#endif
#ifdef SOLVE3D
# if defined ANA_STFLUX || defined ANA_SSFLUX
      public ana_stflux_tile
# endif
#endif

      contains

! ----------------------------------------------------------------------
#ifdef ANA_INITIAL
      subroutine ana_init ![

      use param
      use scalars
      use boundary
      use tracers
      use eos_vars
      use mixing

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

# include "compute_auxiliary_bounds.h"


! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.


# if defined OBC_WEST  || defined OBC_EAST ||\
     defined OBC_SOUTH || defined OBC_NORTH
#  ifdef ANA_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR                     ! Save perimeter values
          zeta_west(j)=zeta(istr-1,j,1)      ! of the initial fields
          ubar_west(j)=ubar(istr  ,j,1)      ! to be used as external
          vbar_west(j)=vbar(istr-1,j,1)      ! forcing data for side
        enddo                                ! boundaries.
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_EAST
      if (EASTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR
          zeta_east(j)=zeta(iend+1,j,1)
          ubar_east(j)=ubar(iend+1,j,1)
          vbar_east(j)=vbar(iend+1,j,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_east(j,k)=u(iend+1,j,k,1)
            v_east(j,k)=v(iend+1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_east(j,k,itemp)=t(iend+1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif


#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_south(i)=zeta(i,jstr-1,1)
          ubar_south(i)=ubar(i,jstr-1,1)
          vbar_south(i)=vbar(i,jstr  ,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,jstr-1,k,1)
            v_south(i,k)=v(i,jstr  ,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,itemp)=t(i,jstr-1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_north(i)=zeta(i,jend+1,1)
          ubar_north(i)=ubar(i,jend+1,1)
          vbar_north(i)=vbar(i,jend+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,jend+1,k,1)
            v_north(i,k)=v(i,jend+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,itemp)=t(i,jend+1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  else /* alternative to ANA_BRY */
#   ifdef UCLIMATOLOGY
      do j=jstrR,jendR                   ! Save initial data into
        do i=istrR,iendR                 ! climatology for the
          ubclm(i,j)=ubar(i,j,1)         ! subsequent use as inflow
          vbclm(i,j)=vbar(i,j,1)         ! boundary conditions
        enddo
      enddo
#   endif
#   ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#    ifdef UCLIMATOLOGY
            uclm(i,j,k)=u(i,j,k,1)
            vclm(i,j,k)=v(i,j,k,1)
#    endif
#    ifdef TCLIMATOLOGY
            tclm(i,j,k,itemp)=t(i,j,k,1,itemp)
#    endif
          enddo
        enddo
      enddo
#   endif /* SOLVE3D */
#  endif /* ANA_BRY */
# endif /* at least one of OBC_XXXXX */


# ifdef EXCHANGE
      call exchange2d_3_tile(istr,iend,jstr,jend,
     &                       zeta(START_2D_ARRAY,1),
     &                       ubar(START_2D_ARRAY,1),
     &                       vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
      call exchange_2_tile(istr,iend,jstr,jend,
     &                     u(START_2D_ARRAY,1,1), N,
     &                     v(START_2D_ARRAY,1,1), N)
      do itrc=1,nt
        call exchange_tile(istr,iend,jstr,jend,
     &                t(START_2D_ARRAY,1,1,itrc), N)
      enddo
#  endif
#  ifdef ISWAKE
#   ifdef UCLIMATOLOGY
      call exchange2d_2_tile(istr,iend,jstr,jend, ubclm,vbclm)
#   endif
#   ifdef SOLVE3D
#    ifdef UCLIMATOLOGY
      call exchange_2_tile(istr,iend,jstr,jend, uclm,N, vclm,N)
#    endif
#    ifdef TCLIMATOLOGY
      call exchange_tile(istr,iend,jstr,jend,
     &                    tclm(START_2D_ARRAY,1,itemp), N)
#    endif
#   endif
#  elif defined COLD_FILAMENT
      call exchange_tile(istr,iend,jstr,jend, Akv,N+1)

#  endif
# endif  /* EXCHANGE */


# if defined OBC_WEST && (defined NS_PERIODIC || defined MPI)
      if (WESTERN_EDGE) then
#  ifdef M2_FRC_BRY
!       call exch_bry_EW_tile (jstr,jend, zeta_west,1)
!       call exch_bry_EW_tile (jstr,jend, ubar_west,1)
!       call exch_bry_EW_tile (jstr,jend, vbar_west,1)
#  endif
#  ifdef M3_FRC_BRY
!       call exch_bry_EW_tile (jstr,jend, u_west,N)
!       call exch_bry_EW_tile (jstr,jend, v_west,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
!         call exch_bry_EW_tile (jstr,jend, t_west(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_EAST && (defined NS_PERIODIC || defined MPI)
      if (EASTERN_EDGE) then
#  ifdef M2_FRC_BRY
!       call exch_bry_EW_tile (jstr,jend, zeta_east,1)
!       call exch_bry_EW_tile (jstr,jend, ubar_east,1)
!       call exch_bry_EW_tile (jstr,jend, vbar_east,1)
#  endif
#  ifdef M3_FRC_BRY
!       call exch_bry_EW_tile (jstr,jend, u_east,N)
!       call exch_bry_EW_tile (jstr,jend, v_east,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
!         call exch_bry_EW_tile (jstr,jend, t_east(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_SOUTH && (defined EW_PERIODIC || defined MPI)
      if (SOUTHERN_EDGE) then
#  ifdef M2_FRC_BRY
!       call exch_bry_NS_tile (istr,iend, zeta_south,1)
!       call exch_bry_NS_tile (istr,iend, ubar_south,1)
!       call exch_bry_NS_tile (istr,iend, vbar_south,1)
#  endif
#  ifdef M3_FRC_BRY
!       call exch_bry_NS_tile (istr,iend, u_south,N)
!       call exch_bry_NS_tile (istr,iend, v_south,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
!         call exch_bry_NS_tile (istr,iend, t_south(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_NORTH && (defined EW_PERIODIC || defined MPI)
      if (NORTHERN_EDGE) then
#  ifdef M2_FRC_BRY
!       call exch_bry_NS_tile (istr,iend, zeta_north,1)
!       call exch_bry_NS_tile (istr,iend, ubar_north,1)
!       call exch_bry_NS_tile (istr,iend, vbar_north,1)
#  endif
#  ifdef M3_FRC_BRY
!       call exch_bry_NS_tile (istr,iend, u_north,N)
!       call exch_bry_NS_tile (istr,iend, v_north,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
!         call exch_bry_NS_tile (istr,iend, t_north(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

      end subroutine ana_init_generic !]
#endif /* ANA_INITIAL */
! ----------------------------------------------------------------------
#if defined ANA_VMIX && defined SOLVE3D
      subroutine ana_vmix (tile) ![

      implicit none
      integer tile
      call ana_vmix_tile 
      end subroutine ana_vmix !]
! ----------------------------------------------------------------------
      subroutine ana_vmix_tile ![

! Set vertical kinematic vertical viscosity "Akv" [m^2/s] for momentum
! and vertical mixing c efficients "Akt" [m^2/s] for tracer equations
! using analytical expressions.


      implicit none
      integer :: i,j,k
      real    :: cff
      real    :: hbl,ssgm

      hbl = 70;  ! match the initial mixed layer depth

      if (FIRST_TIME_STEP) then
        do k=0,nz
          do j=0,ny+1
            do i=0,nx+1
              ssgm=(z_w(i,j,nz)-z_w(i,j,k))/hbl

              if (ssgm < 1.) then
                if (ssgm<0.07D0) then
                  cff=0.5*(ssgm-0.07D0)**2/0.07D0
                else
                  cff=0.D0
                endif
                cff=cff + ssgm*(1.-ssgm)**2
              else
                cff=0
              endif

              ! attempt to get a max value of 0.005 m^2/s
              Akv(i,j,k) = cff*0.005/0.15

              Akt(i,j,k,itemp)=cff*0.005/0.15
#  ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
            enddo
          enddo
        enddo
      endif


# ifdef EXCHANGE
!     call exchange_2_tile (istr,iend,jstr,jend,  Akv, N+1,
!    &                    Akt(START_2D_ARRAY,0,itemp), N+1)
!  ifdef SALINITY
!     call exchange_tile (istr,iend,jstr,jend,
!    &                    Akt(START_2D_ARRAY,0,isalt), N+1)
!  endif
# endif

      end subroutine ana_vmix_tile !]
#endif /* ANA_VMIX */

! ----------------------------------------------------------------------

! ONLY THE _TILE ROUTINES ARE USED in set_forces! GET RID OF THE NON-TILE VERSIONS

! ----------------------------------------------------------------------
#ifdef ANA_SMFLUX
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
#endif /* ANA_SMFLUX */

! ----------------------------------------------------------------------
#if defined ANA_SRFLUX && defined SOLVE3D
      subroutine ana_srflux_tile(istr,iend,jstr,jend) ![

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
          srflx(i,j) = 0.
        enddo
      enddo

      end subroutine ana_srflux_tile !]
#endif /* ANA_SRFLUX */

! ----------------------------------------------------------------------
# if defined ANA_STFLUX || defined ANA_SSFLUX
      subroutine ana_stflux_tile(istr,iend,jstr,jend, itrc) ![

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
      end subroutine ana_stflux_tile !]
# endif /* ANA_STFLUX || ANA_SSFLUX */

! ----------------------------------------------------------------------

      end module analytical
