      module tides

#include "cppdefs.h"
! Empty module if TIDES not defined.
#ifdef TIDES

      implicit none


      ! Make all variable private to this module unless public specified
      private

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_array to work.
#include "param.h"
      ! For iic, ntstart, deg2rad
#include "scalars.h"

      ! ****************************************************************
      ! user inputs

!** Ntides       Number of tidal components to consider.
#if defined SSH_TIDES || defined UV_TIDES
      integer, parameter, public :: Ntides=15
#endif

      ! end user inputs
      ! ****************************************************************


!** Tperiod      Tidal period (s).
!** SSH_Tamp     Tidal elevation amplitude (m) at RHO-points.
!** SSH_Tphase   Tidal elevation phase (degrees/360) at RHO-points.
!** UV_Tangle    Tidal current angle (radians; counterclockwise
!**                from EAST and rotated to curvilinear grid) at
!**                RHO-points.
!** UV_Tmajor    Maximum tidal current: tidal ellipse major axis
!**                (m/s) at RHO-points.
!** UV_Tminor    Minimum tidal current: tidal ellipse minor axis
!**                (m/s) at RHO-points.
!** UV_Tphase    Tidal current phase (degrees/360) at RHO-points.

      ! Get rid of some flags. Also decide what needs to be public if any.
#if defined SSH_TIDES || defined UV_TIDES
!      integer, parameter, public :: Ntides=15
      real, public :: Tperiod(Ntides)
#endif
#if defined SSH_TIDES
      real, public :: SSH_Tamp(GLOBAL_2D_ARRAY,Ntides)
      real, public :: SSH_Tphase(GLOBAL_2D_ARRAY,Ntides)
#endif
#if defined UV_TIDES
      real, public :: UV_Tangle(GLOBAL_2D_ARRAY,Ntides)
      real, public :: UV_Tmajor(GLOBAL_2D_ARRAY,Ntides)
      real, public :: UV_Tminor(GLOBAL_2D_ARRAY,Ntides)
      real, public :: UV_Tphase(GLOBAL_2D_ARRAY,Ntides)
#endif

#ifdef POT_TIDES
      real, public :: POT_Tamp(GLOBAL_2D_ARRAY,Ntides) ! Public?
      real, public :: POT_Tphase(GLOBAL_2D_ARRAY,Ntides) ! Public?
      real, public :: PTide(GLOBAL_2D_ARRAY)
#endif

      public set_tides

      contains
! ----------------------------------------------------------------------

      subroutine set_tides(tile)

      ! DevinD created seperate step as can't have compute_tile_bounds
      ! and compute_auxilliary_bounds in the same subroutine as they both
      ! declare and compute. Fortan doesn't allow.

      implicit none

      integer tile
#include "private_scratch.h"
#include "compute_tile_bounds.h"

      call set_tides_tile(istr,iend,jstr,jend,A2d(1,1),A2d(1,2),
     &                                            A2d(1,3),A2d(1,4))

      end subroutine set_tides
!-----------------------------------------------------------------------

      subroutine set_tides_tile(istr,iend,jstr,jend,Cangle,Sangle,
     &                                                  Cphase,Sphase)

      implicit none

! scalars.h for iic and ntstart
! grid.h for rmask
# include "grid.h"
# include "boundary.h"

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Cangle, Cphase,
     &                                             Sangle, Sphase
      ! local
      integer tile, itide, i, j
      real angle, phase, ramp, omega

#include "compute_auxiliary_bounds.h"

      ! Read tides (only once)
      if (iic==ntstart) then
        call read_tides(istr,iend,jstr,jend)
      endif

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

# if defined POT_TIDES
!-----------------------------------------------------------------------
!  Compute tidal potential (m)
!  -- > to be applied as a pressure-gradient force
!       in prsgrd (if SOLVE3D) or step2d (if not SOLVE3D)
!-----------------------------------------------------------------------
!
          do j=jstrR-1,jendR ! DevinD does jstrR-1 not go beyond computational margin?
            do i=istrR-1,iendR
              Ptide(i,j)=( Ptide(i,j)+
     &                   ramp*POT_Tamp(i,j,itide)*
     &                   cos(omega-POT_Tphase(i,j,itide)) )
#  ifdef MASKING
     &                   *rmask(i,j)
#  endif
            enddo
          enddo
# endif /* POT_TIDES */

        endif  !<--- period > 0
      enddo  !<-- itide

      end subroutine set_tides_tile
!-----------------------------------------------------------------------

      subroutine read_tides(istr,iend,jstr,jend)
      ! Read tides from forcing file. Only read once at first timestep.

      use read_write
      use netcdf

      implicit none

      ! For max_frc
#include "ncvars.h"
      ! scalars.h for deg2rad

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer ierr, var_file_indx, ncid, var_id, itide, found_rec,ierr2,i,j

      found_rec = 0     ! Set found variable flag to not found
      var_file_indx = 1 ! Forcing file index number
      ! Since tides are only read once, don't need to store netcdf vars

      ! Continue searching while tides not found & still more files to check
      do while ((found_rec==0).and.(var_file_indx <= max_frc))

        ! Try to open force file
 100    ierr=nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)
        do itide=1,Ntides
          ierr=nf90_inq_varid(ncid, 'tide_period', var_id)
          if(ierr==nf90_noerr) then

            ! Tperiod
            ierr=nf90_get_var(ncid,var_id,Tperiod(itide),(/itide/))
            if(ierr/=0) call handle_ierr(ierr,'tides::tide_period')
#ifdef SSH_TIDES
            ! tidal elevation amplitude and phase.
            ierr =nf90_inq_varid(ncid,'tide_Eamp',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    SSH_Tamp(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Eamp')
            ! tide_Ephase
            ierr =nf90_inq_varid(ncid,'tide_Ephase',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    SSH_Tphase(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Ephase')
#endif
#ifdef UV_TIDES
            ! Read tidal currents angle, phase, major and minor ellipse axis.
            ierr =nf90_inq_varid(ncid,'tide_Cangle',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    UV_Tangle(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Cangle')
            ! tide_Cphase
            ierr =nf90_inq_varid(ncid,'tide_Cphase',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    UV_Tphase(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Cphase')
            ! tide_Cmax
            ierr =nf90_inq_varid(ncid,'tide_Cmax',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    UV_Tmajor(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Cmax')
            ! tide_Cmin
            ierr =nf90_inq_varid(ncid,'tide_Cmin',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    UV_Tminor(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Cmin')

#endif
#ifdef POT_TIDES
            ! tide_Pamp
            ierr =nf90_inq_varid(ncid,'tide_Pamp',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    POT_Tamp(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Pamp')
            ! tide_Pphase
            ierr =nf90_inq_varid(ncid,'tide_Pphase',var_id)
            ierr2=ncdf_read_mod(ncid, var_id, itide, r_var,
     &                    POT_Tphase(START_2D_ARRAY,itide), 1)
            if(ierr/=0 .or. ierr2/=0) call handle_ierr(ierr,'tides::tide_Pphase')
#endif

            ! Convert tidal phase and angle from degrees to radians.
            Tperiod(itide)=Tperiod(itide)*3600.
            do j=0,Mm+1
              do i=0,Lm+1
# ifdef SSH_TIDES
                SSH_Tphase(i,j,itide)=SSH_Tphase(i,j,itide)*deg2rad
# endif
# ifdef UV_TIDES
                UV_Tangle(i,j,itide)=UV_Tangle(i,j,itide)*deg2rad
                UV_Tphase(i,j,itide)=UV_Tphase(i,j,itide)*deg2rad
# endif
# ifdef POT_TIDES
                POT_Tphase(i,j,itide)=POT_Tphase(i,j,itide)*deg2rad
# endif
                enddo
              enddo

            found_rec=1 ! Set record to found to exit do while loop
          else
            var_file_indx=var_file_indx+1 ! Set forcing file index to next file
            go to 100 ! Try next file
          endif

        end do ! itide=1,Ntides
      end do ! found_rec==0 && var_file_indx <= max_frc

      if (found_rec==0 .and. var_file_indx==max_frc) then
        error stop 'tides.F: tide input data not found'
      endif

      if(mynode==0) write(*,'(5x,A)')'tides:: read all tides'

      end subroutine read_tides
!-----------------------------------------------------------------------

#endif /* TIDES for whole module */

      end module tides
