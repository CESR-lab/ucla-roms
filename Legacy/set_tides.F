#include "cppdefs.h"
#if defined SSH_TIDES || defined UV_TIDES

      subroutine set_tides (tile)

!================================================== Robert Hetland ===
! Copyright (c) 2000 Rutgers/UCLA
!================================================ Hernan G. Arango ===
!
! This routine adds tidal elevation (m) and tidal currents (m/s) to
! sea surface height and 2D momentum climatologies, respectively.
!
!=====================================================================

      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_tides_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                                          A2d(1,3), A2d(1,4))
      end

      subroutine set_tides_tile (istr,iend,jstr,jend, Cangle, Sangle,
     &                                                Cphase, Sphase)

      use tides

      implicit none
# include "param.h"
# include "scalars.h"
# include "climat.h"
# include "grid.h"
!# include "tides.h"
# include "boundary.h"

      integer istr,iend,jstr,jend,    itide, i,j
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Cangle, Cphase,
     &                                             Sangle, Sphase
      real angle, phase, omega, ramp

#  include "compute_auxiliary_bounds.h"

!# ifdef TIDES_DEBUG
!      if (iic==ntstart) then
!        call read_tides(istr,iend,jstr,jend)
!      endif
!# endif

      ramp=1.
      ! JM confirmed no longer want tide ramp
!# ifdef TIDERAMP
!      ramp=TANH(dt*sec2day*float(iic-ntstart))
!# endif

      do itide=1,Ntides
        if (Tperiod(itide) > 0.) then
          omega=2.*pi*time/Tperiod(itide)

# if defined SSH_TIDES
#  if defined OBC_WEST
          if (WESTERN_EDGE) then
            do j=jstrR,jendR
#   ifdef Z_FRC_BRY
              zeta_west(j) =zeta_west(j)
#   else
              ssh(istr-1,j)=ssh(istr-1,j)
#   endif
     &              +ramp*SSH_Tamp(istr-1,j,itide)
     &       *cos(omega-SSH_Tphase(istr-1,j,itide))
            enddo
          endif
#  endif
#  if defined OBC_EAST
          if (EASTERN_EDGE) then
            do j=jstrR,jendR
#   ifdef Z_FRC_BRY
              zeta_east(j) =zeta_east(j)
#   else
              ssh(iend+1,j)=ssh(iend+1,j)
#   endif
     &                   +ramp*SSH_Tamp(iend+1,j,itide)
     &            *cos(omega-SSH_Tphase(iend+1,j,itide))
            enddo
          endif
#  endif
#  if defined OBC_SOUTH
          if (SOUTHERN_EDGE) then
            do i=istrR,iendR
#   ifdef Z_FRC_BRY
              zeta_south(i)=zeta_south(i)
#   else
              ssh(i,jstr-1)=ssh(i,jstr-1)
#   endif
     &                   +ramp*SSH_Tamp(i,jstr-1,itide)
     &            *cos(omega-SSH_Tphase(i,jstr-1,itide))
            enddo
          endif
#  endif
#  if defined OBC_NORTH
          if (NORTHERN_EDGE) then
            do i=istrR,iendR
#   ifdef Z_FRC_BRY
              zeta_north(i)=zeta_north(i)
#   else
              ssh(i,jend+1)=ssh(i,jend+1)
#   endif
     &                   +ramp*SSH_Tamp(i,jend+1,itide)
     &            *cos(omega-SSH_Tphase(i,jend+1,itide))
            enddo
          endif
#  endif
# endif /* SSH_TIDES */
# if defined UV_TIDES
#  ifdef OBC_WEST
          if (WESTERN_EDGE) THEN
            do i=istr-1,istr
              do j=jstr-1,jendR
                angle=UV_Tangle(i,j,itide)-angler(i,j)
                phase=omega-UV_Tphase(i,j,itide)
                Cangle(i,j)=cos(angle)
                Cphase(i,j)=cos(phase)
                Sangle(i,j)=sin(angle)
                Sphase(i,j)=sin(phase)
              enddo
            enddo
            do j=jstrR,jendR
              i=istrU-1
#   ifdef M2_FRC_BRY
              ubar_west(j)=ubar_west(j)
#   else
              ubclm(i,j)=ubclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i-1,j,itide)+UV_Tmajor(i  ,j,itide))
     &       *(Cangle(i-1,j)+Cangle(i,j))*(Cphase(i-1,j)+Cphase(i,j))

     &               -(UV_Tminor(i-1,j,itide)+UV_Tminor(i  ,j,itide))
     &       *(Sangle(i-1,j)+Sangle(i,j))*(Sphase(i-1,j)+Sphase(i,j))
     &                                                              )
            enddo
            do j=jstr,jendR
              i=istr-1
#   ifdef M2_FRC_BRY
              vbar_west(j)=vbar_west(j)
#   else
              vbclm(i,j)=vbclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i,j-1,itide)+UV_Tmajor(i,j  ,itide))
     &       *(Sangle(i,j-1)+Sangle(i,j))*(Cphase(i,j-1)+Cphase(i,j))

     &               +(UV_Tminor(i,j-1,itide)+UV_Tminor(i,j  ,itide))
     &       *(Cangle(i,j-1)+Cangle(i,j))*(Sphase(i,j-1)+Sphase(i,j))
     &                                                              )
            enddo
          endif
#  endif

#  ifdef OBC_EAST
          if (EASTERN_EDGE) THEN
            do i=iend,iend+1
              do j=jstr-1,jendR
                angle=UV_Tangle(i,j,itide)-angler(i,j)
                phase=omega-UV_Tphase(i,j,itide)
                Cangle(i,j)=cos(angle)
                Cphase(i,j)=cos(phase)
                Sangle(i,j)=sin(angle)
                Sphase(i,j)=sin(phase)
              enddo
            enddo
            do j=jstrR,jendR
              i=iend+1
#   ifdef M2_FRC_BRY
              ubar_east(j)=ubar_east(j)
#   else
              ubclm(i,j)=ubclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i-1,j,itide)+UV_Tmajor(i  ,j,itide))
     &                                   *(Cangle(i-1,j)+Cangle(i,j))
     &                                   *(Cphase(i-1,j)+Cphase(i,j))

     &               -(UV_Tminor(i-1,j,itide)+UV_Tminor(i  ,j,itide))
     &                                   *(Sangle(i-1,j)+Sangle(i,j))
     &                                   *(Sphase(i-1,j)+Sphase(i,j))
     &                                                              )
            enddo
            do j=jstr,jendR
              i=iend+1
#   ifdef M2_FRC_BRY
              vbar_east(j)=vbar_east(j)
#   else
              vbclm(i,j)=vbclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i,j-1,itide)+UV_Tmajor(i,j  ,itide))
     &       *(Sangle(i,j-1)+Sangle(i,j))*(Cphase(i,j-1)+Cphase(i,j))

     &               +(UV_Tminor(i,j-1,itide)+UV_Tminor(i,j  ,itide))
     &       *(Cangle(i,j-1)+Cangle(i,j))*(Sphase(i,j-1)+Sphase(i,j))
     &                                                              )
            enddo
          endif
#  endif

#  ifdef OBC_SOUTH
          if (SOUTHERN_EDGE) THEN
            do j=jstr-1,jstr
              do i=istr-1,iendR
                angle=UV_Tangle(i,j,itide)-angler(i,j)
                phase=omega-UV_Tphase(i,j,itide)
                Cangle(i,j)=cos(angle)
                Cphase(i,j)=cos(phase)
                Sangle(i,j)=sin(angle)
                Sphase(i,j)=sin(phase)
              enddo
            enddo
            do i=istr,iendR
              j=jstr-1
#   ifdef M2_FRC_BRY
              ubar_south(i)=ubar_south(i)
#   else
              ubclm(i,j)=ubclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i-1,j,itide)+UV_Tmajor(i  ,j,itide))
     &       *(Cangle(i-1,j)+Cangle(i,j))*(Cphase(i-1,j)+Cphase(i,j))

     &               -(UV_Tminor(i-1,j,itide)+UV_Tminor(i  ,j,itide))
     &       *(Sangle(i-1,j)+Sangle(i,j))*(Sphase(i-1,j)+Sphase(i,j))
     &                                                              )
            enddo
            do i=istrR,iendR
              j=jstrV-1
#   ifdef M2_FRC_BRY
              vbar_south(i)=vbar_south(i)
#   else
              vbclm(i,j)=vbclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i,j-1,itide)+UV_Tmajor(i,j  ,itide))
     &       *(Sangle(i,j-1)+Sangle(i,j))*(Cphase(i,j-1)+Cphase(i,j))

     &               +(UV_Tminor(i,j-1,itide)+UV_Tminor(i,j  ,itide))
     &       *(Cangle(i,j-1)+Cangle(i,j))*(Sphase(i,j-1)+Sphase(i,j))
     &                                                              )
            enddo
          endif
#  endif
#  ifdef OBC_NORTH
          if (NORTHERN_EDGE) THEN
            do j=jend,jend+1
              do i=istr-1,iendR
                angle=UV_Tangle(i,j,itide)-angler(i,j)
                phase=omega-UV_Tphase(i,j,itide)
                Cangle(i,j)=cos(angle)
                Cphase(i,j)=cos(phase)
                Sangle(i,j)=sin(angle)
                Sphase(i,j)=sin(phase)
              enddo
            enddo
            do i=istr,iendR
              j=jend+1
#   ifdef M2_FRC_BRY
              ubar_north(i)=ubar_north(i)
#   else
              ubclm(i,j)=ubclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i-1,j,itide)+UV_Tmajor(i  ,j,itide))
     &       *(Cangle(i-1,j)+Cangle(i,j))*(Cphase(i-1,j)+Cphase(i,j))

     &               -(UV_Tminor(i-1,j,itide)+UV_Tminor(i  ,j,itide))
     &       *(Sangle(i-1,j)+Sangle(i,j))*(Sphase(i-1,j)+Sphase(i,j))
     &                                                              )
            enddo

            do i=istrR,iendR
              j=jend+1
#   ifdef M2_FRC_BRY
              vbar_north(i)=vbar_north(i)
#   else
              vbclm(i,j)=vbclm(i,j)
#   endif
     &  +ramp*0.125*( (UV_Tmajor(i,j-1,itide)+UV_Tmajor(i,j  ,itide))
     &       *(Sangle(i,j-1)+Sangle(i,j))*(Cphase(i,j-1)+Cphase(i,j))

     &               +(UV_Tminor(i,j-1,itide)+UV_Tminor(i,j  ,itide))
     &       *(Cangle(i,j-1)+Cangle(i,j))*(Sphase(i,j-1)+Sphase(i,j))
     &                                                              )
            enddo
          endif
#  endif
# endif /* UV_TIDES */

        endif  !<--- period > 0
      enddo  !<-- itide

      ! DevinD debug
      if(mynode==0) then
        print*,'OLD VERSION'
        print*,'zeta_north(20)',zeta_north(20)
        print*,'zeta_west(20)',zeta_west(20)
        print*,'zeta_south(20)',zeta_south(20)
        print*,'zeta_east(20)',zeta_east(20)
        !
        print*,'ubar_north(20)',ubar_north(20),'vbar_north(20)',vbar_north(20)
        print*,'ubar_west(20)',ubar_west(20),'vbar_west(20)',vbar_west(20)
        print*,'ubar_south(20)',ubar_south(20),'vbar_south(20)',vbar_south(20)
        print*,'ubar_east(20)',ubar_east(20),'vbar_east(20)',vbar_east(20)
      endif

      end
#else
      subroutine set_tides_empty
      end
#endif /* SSH_TIDES || UV_TIDES */
