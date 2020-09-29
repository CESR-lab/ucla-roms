      module tides

#include "cppdefs.h"
! Empty module is TIDES not defined.
#ifdef TIDES

      implicit none


      ! Make all variable private to this module unless public specified
      private

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_array to work.
#include "param.h"
      ! tides.h for variable Ntides
#include "tides.h"

      ! ****************************************************************

#ifdef POT_TIDES
      real, public :: POT_Tamp(GLOBAL_2D_ARRAY,Ntides) ! Public?
      real, public :: POT_Tphase(GLOBAL_2D_ARRAY,Ntides) ! Public?
      real, public ::  PTide(GLOBAL_2D_ARRAY)
#endif

      public set_tides_mod

      contains
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
      subroutine set_tides_mod(tile)

      ! DevinD created seperate step as can't have compute_tile_bounds
      ! and compute_auxilliary_bounds in the same subroutine as they both
      ! declare and compute. Fortan doesn't allow.

      implicit none

      integer tile

#include "compute_tile_bounds.h"



      call set_tides_mod_tile(istr,iend,jstr,jend)

      end subroutine set_tides_mod
!-----------------------------------------------------------------------

! ----------------------------------------------------------------------
      subroutine set_tides_mod_tile(istr,iend,jstr,jend)

      implicit none

! scalars.h for iic and ntstart
#include "scalars.h"
! grid.h for rmask
# include "grid.h"

      integer tile, itide, i, j, istr,iend,jstr,jend
      real ramp, omega

#include "compute_auxiliary_bounds.h"

! CURRENTLY READING STILL DONE IN GET_FORCES, NEED TO BE PORTED TO MODULE.

      ! Tides only read in once, hence if
!      if (iic==ntstart) then



!#ifdef POT_TIDES
!      call set_pot_tides(istr,iend,jstr,jend)
!#endif

!      endif

      ramp = 1.

# if defined POT_TIDES
      do j=jstrR-1,jendR
        do i=istrR-1,iendR
          Ptide(i,j)=0.
        enddo
      enddo
# endif /* POT_TIDES */

      do itide=1,Ntides
        if (Tperiod(itide) > 0.) then

          omega=2.*pi*time/Tperiod(itide)

# if defined POT_TIDES
!-----------------------------------------------------------------------
!  Compute tidal potential (m)
!  -- > to be applied as a pressure-gradient force
!       in prsgrd (if SOLVE3D) or step2d (if not SOLVE3D)
!-----------------------------------------------------------------------
!
          do j=jstrR-1,jendR
            do i=istrR-1,iendR
              Ptide(i,j)=Ptide(i,j)+
     &                   ramp*POT_Tamp(i,j,itide)*
     &                   cos(omega-POT_Tphase(i,j,itide))
#  ifdef MASKING
              Ptide(i,j)=Ptide(i,j)*rmask(i,j) ! DevinD change to *rmask only
#  endif
            enddo
          enddo
# endif /* POT_TIDES */

          ! DevinD debug
          if(mynode==0) then
            print *, 'itide',itide,'omega',omega
            print *, 'POT_Tamp(20,20,itide)',POT_Tamp(20,20,itide)
            print *, 'POT_Tphase(20,20,itide)',POT_Tphase(20,20,itide)
            print *, 'mynode==0: Ptide(20,20)', Ptide(20,20)
          end if

        endif  !<--- period > 0
      enddo  !<-- itide

!      stop 'Stop in set_tides_mod'

      end subroutine set_tides_mod_tile
!-----------------------------------------------------------------------

! ----------------------------------------------------------------------
!      subroutine set_pot_tides(istr,iend,jstr,jend)
!
!      implicit none
!
!      ! input/outputs
!      integer,intent(in)  :: istr,iend,jstr,jend
!!      integer itide
!
!
!
!      end subroutine set_pot_tides
!-----------------------------------------------------------------------

#endif /* TIDES for whole module */

      end module tides
