#define GRID_LEVEL 2
#include "cppdefs.h"
#ifdef GRID_LEVEL
# if GRID_LEVEL > 1

      subroutine u2dbc_tile (istr,iend,jstr,jend, zeta_new,buff)

      use param
      use grid

      implicit none
      integer istr,iend,jstr,jend, i,j, jmin,jmax
      real zeta_new(PRIVATE_2D_SCRATCH_ARRAY),
     &         buff(PRIVATE_2D_SCRATCH_ARRAY)

#  include "scalars.h"
#  include "ocean2d.h"
                                ! Preload broker array "buff" with
      if (SOUTHERN_EDGE) then   ! masked cross-sections of the fine-
        jmin=jstr-1             ! grid elements, accross which normal
      else                      ! velocity component exports/imports
        jmin=jstr-2             ! mass in/out of the domain. These
      endif                     ! These crossections are then used
      if (NORTHERN_EDGE) then   ! inside the interpolation routine
        jmax=jend+1             ! to enforce flux-consistency between
      else                      ! the fine and coarse grids within
        jmax=jend+2             ! each fine-grid triad. Note that the
      endif                     ! interpolation routines return mass
                                ! fluxes, which are then need to be
      if (WESTERN_EDGE) then    ! divided by the crossections.
        do j=jmin,jmax
          buff(istr-1,j)=0.5*( zeta_new(istr-1,j)+h(istr-1,j)
     &                        +zeta_new(istr  ,j)+h(istr  ,j))
     &                                           *dm_u(istr,j)
          buff(istr,j)=buff(istr-1,j)
#  ifdef MASKING
     &                 *umask(istr,j)
#  endif
        enddo
#  if GRID_LEVEL == 2
        call int_u2d_west_1 (istr,iend,jstr,jend, buff)
#  elif GRID_LEVEL == 3
        call int_u2d_west_2 (istr,iend,jstr,jend, buff)
#  elif GRID_LEVEL == 4
        call int_u2d_west_3 (istr,iend,jstr,jend, buff)
#  endif
        do j=jstr,jend
          ubar(istr,j,knew)=buff(istr,j)/buff(istr-1,j)
        enddo
      endif

      if (EASTERN_EDGE) then
        do j=jmin,jmax
          buff(iend+2,j)=0.5*( zeta_new(iend+1,j)+h(iend+1,j)
     &                        +zeta_new(iend  ,j)+h(iend  ,j))
     &                                         *dm_u(iend+1,j)
          buff(iend+1,j)=buff(iend+2,j)
#  ifdef MASKING
     &                 *umask(iend+1,j)
#  endif
        enddo
#  if GRID_LEVEL == 2
        call int_u2d_east_1 (istr,iend,jstr,jend, buff)
#  elif GRID_LEVEL == 3
        call int_u2d_east_2 (istr,iend,jstr,jend, buff)
#  elif GRID_LEVEL == 4
        call int_u2d_east_3 (istr,iend,jstr,jend, buff)
#  endif
        do j=jstr,jend
          ubar(iend+1,j,knew)=buff(iend+1,j)/buff(iend+2,j)
        enddo
      endif





      if (SOUTHERN_EDGE) then
      endif
      if (NORTHERN_EDGE) then
      endif




      return
      end
# endif
# if GRID_LEVEL < MAX_GRID_LEVEL

      subroutine int_u2d_west (istr,iend,jstr,jend, buff)

      use param
      use coupling

      implicit none
      integer istr,iend,jstr,jend, jmin,jmax, j,jf
      real buff(PRIVATE_2D_SCRATCH_ARRAY)
      real T1,T4,T5, crs,crsR,crsL, uctr,dltR,dltL, cff,cffR,cffL
      parameter (T1=1./27., T4=4./27., T5=5./27.)
#  include "scalars.h"
#  include "ocean2d.h"

      jmin=jmin_child+(jstr-1)/3
      jmax=jmin_child+(jend-1)/3

      do j=jmin,jmax                  ! At entry "buff" is preloaded
        jf=2+3*(j-jmin_child)         ! with masked crossections from
        crs=buff(istr,jf  )           ! the fine grid.  Save them to
        crsR=buff(istr,jf+1)          ! restore integral flux in each
        crsL=buff(istr,jf-1)          ! triade.

        uctr=ubar(imin_child,j,knew)
        dltR=ubar(imin_child,j+1,knew)-uctr
        dltL=uctr-ubar(imin_child,j-1,knew)

        cffL=4.*dltL
        cffR=4.*dltR
        if (dltL*dltR <0.) then
          dltL=0.
          dltR=0.
        elseif (abs(dltL) > abs(cffR)) then
          dltL=cffR
        elseif (abs(dltR) > abs(cffL)) then
          dltR=cffL
        endif

        buff(istr,jf  )=crs *( uctr - T1*(dltR-dltL) )
        buff(istr,jf+1)=crsR*( uctr +T5*dltR +T4*dltL )
        buff(istr,jf-1)=crsL*( uctr -T4*dltR -T5*dltL )

        cff=crs+crsR+crsL
        if (cff>0.) cff=( DU_avg2(imin_child,j) -buff(istr,jf)
     &                     -buff(istr,jf+1)-buff(istr,jf-1) )/cff

        buff(istr,jf  )=buff(istr,jf  ) +cff*crs
        buff(istr,jf+1)=buff(istr,jf+1) +cff*crsR
        buff(istr,jf-1)=buff(istr,jf-1) +cff*crsL
      enddo
      return
      end

      subroutine int_u2d_east (istr,iend,jstr,jend, buff)

      use param
      use coupling

      implicit none
      integer istr,iend,jstr,jend, jmin,jmax, j,jf
      real buff(PRIVATE_2D_SCRATCH_ARRAY)
      real T1,T4,T5, crs,crsR,crsL, uctr,dltR,dltL, cff,cffR,cffL
      parameter (T1=1./27., T4=4./27., T5=5./27.)
#  include "scalars.h"
#  include "ocean2d.h"

      jmin=jmin_child+(jstr-1)/3
      jmax=jmin_child+(jend-1)/3

      do j=jmin,jmax
        jf=2+3*(j-jmin_child)
        crs=buff(iend+1,jf  )
        crsR=buff(iend+1,jf+1)
        crsL=buff(iend+1,jf-1)

        uctr=ubar(imax_child,j,knew)
        dltR=ubar(imax_child,j+1,knew)-uctr
        dltL=uctr-ubar(imax_child,j-1,knew)

        cffL=4.*dltL
        cffR=4.*dltR
        if (dltL*dltR <0.) then
          dltL=0.
          dltR=0.
        elseif (abs(dltL) > abs(cffR)) then
          dltL=cffR
        elseif (abs(dltR) > abs(cffL)) then
          dltR=cffL
        endif

        buff(iend+1,jf  )=crs *( uctr- T1*(dltR-dltL) )
        buff(iend+1,jf+1)=crsR*( uctr+T5*dltR +T4*dltL )
        buff(iend+1,jf-1)=crsL*( uctr-T4*dltR -T5*dltL )

        cff=crs+crsR+crsL
        if (cff>0.) cff=( DU_avg2(imin_child,j) -buff(iend+1,jf)
     &                   -buff(iend+1,jf+1)-buff(iend+1,jf-1) )/cff

        buff(iend+1,jf  )=buff(iend+1,jf  ) +cff*crs
        buff(iend+1,jf+1)=buff(iend+1,jf+1) +cff*crsR
        buff(iend+1,jf-1)=buff(iend+1,jf-1) +cff*crsL
      enddo
      return
      end
# endif
#else
      subroutine int_2d_empty
      end
#endif

