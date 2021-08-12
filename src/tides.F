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
#include "tides.opt"


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

      ! Computes current forcing for boundaries and surface tidal
      ! potential

      ! DevinD created seperate step as can't have compute_tile_bounds
      ! and compute_auxilliary_bounds in the same subroutine as they both
      ! declare and compute. Fortan doesn't allow.

      implicit none

      integer tile
#include "private_scratch.h"
#include "compute_tile_bounds.h"

      call set_tides_tile(istr,iend,jstr,jend)

      end subroutine set_tides
!-----------------------------------------------------------------------

      subroutine set_tides_tile(istr,iend,jstr,jend)

      implicit none

! scalars.h for iic and ntstart
! grid.h for rmask
# include "grid.h"
# include "boundary.h"

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer tile, itide, i, j
      real om

#include "compute_auxiliary_bounds.h"

      ! Read tides (only once)
      if (iic==ntstart) then
        call read_tides(istr,iend,jstr,jend)
      endif

      ramp = 1.

# if defined POT_TIDES
# endif /* POT_TIDES */

      do itide=1,Ntides

        omega=2.*pi*time/Tperiod(itide)

        if (bry_tides) then

#  ifdef OBC_WEST
          if (WESTERN_EDGE) THEN
            i=istr-1
            do j=jstrR,jendR
              zeta_west(j) =zeta_west(j)
     &        + ztide_r(i,j,itide)*cos(om) + ztide_i(i,j,itide)*sin(om)
            enddo
            i=istrU-1
            do j=jstrR,jendR
              ubar_west(j)=ubar_west(j)
     &        + utide_r(i,j,itide)*cos(om) + utide_i(i,j,itide)*sin(om)
            enddo
            i=istr-1
            do j=jstr,jendR
              vbar_west(j)=vbar_west(j)
     &        + vtide_r(i,j,itide)*cos(om) + vtide_i(i,j,itide)*sin(om)
            enddo
          endif
#  endif

#  ifdef OBC_EAST
          if (EASTERN_EDGE) THEN
            i = iend+1
            do j=jstrR,jendR
              zeta_east(j) =zeta_east(j)
     &        + ztide_r(i,j,itide)*cos(om) + ztide_i(i,j,itide)*sin(om)
            enddo
            do j=jstrR,jendR
              ubar_east(j)=ubar_east(j)
     &        + utide_r(i,j,itide)*cos(om) + utide_i(i,j,itide)*sin(om)
            enddo
            do j=jstr,jendR
              vbar_east(j)=vbar_east(j)
     &        + vtide_r(i,j,itide)*cos(om) + vtide_i(i,j,itide)*sin(om)
            enddo
          endif
#  endif

#  ifdef OBC_SOUTH
          if (SOUTHERN_EDGE) THEN
            j=jstr-1
            do i=istrR,iendR
              zeta_south(i)=zeta_south(i)
     &        + ztide_r(i,j,itide)*cos(om) + ztide_i(i,j,itide)*sin(om)
            enddo
            j=jstr-1
            do i=istr,iendR
              ubar_south(i)=ubar_south(i)
     &        + utide_r(i,j,itide)*cos(om) + utide_i(i,j,itide)*sin(om)
            enddo
            j=jstrV-1
            do i=istrR,iendR
              vbar_south(i)=vbar_south(i)
     &        + vtide_r(i,j,itide)*cos(om) + vtide_i(i,j,itide)*sin(om)
            enddo
          endif
#  endif
#  ifdef OBC_NORTH
          if (NORTHERN_EDGE) THEN
            j=jend+1
            do i=istrR,iendR
              zeta_north(i)=zeta_north(i)
     &        + ztide_r(i,j,itide)*cos(om) + ztide_i(i,j,itide)*sin(om)
            enddo
            do i=istr,iendR
              ubar_north(i)=ubar_north(i)
     &        + utide_r(i,j,itide)*cos(om) + utide_i(i,j,itide)*sin(om)
            enddo
            do i=istrR,iendR
              vbar_north(i)=vbar_north(i)
     &        + vtide_r(i,j,itide)*cos(om) + vtide_i(i,j,itide)*sin(om)
            enddo
          endif
#  endif
        endif  ! bry_tides

        if (pot_tides) then
!-----------------------------------------------------------------------
!  Compute tidal potential (m)
!  -- > to be applied as a pressure-gradient force
!       in prsgrd (if SOLVE3D) or step2d (if not SOLVE3D)
!-----------------------------------------------------------------------
          if (itide==1) then ! initialize at the first tidal component
            do j=jstrR-1,jendR
              do i=istrR-1,iendR
                Ptide(i,j)=0.
     &        + ptide_r(i,j,itide)*cos(om) + ptide_i(i,j,itide)*sin(om)
              enddo
            enddo
          else
            do j=jstrR-1,jendR ! DevinD does jstrR-1 not go beyond computational margin?
              do i=istrR-1,iendR
                Ptide(i,j)=( Ptide(i,j)+
     &          + ptide_r(i,j,itide)*cos(om) + ptide_i(i,j,itide)*sin(om)
              enddo
            enddo
          endif

        endif  ! bry_tides

      enddo  !  itide=1,ntides

      end subroutine set_tides_tile
!-----------------------------------------------------------------------

      subroutine read_tides(istr,iend,jstr,jend)
      ! JM  shouldn't this be master only and so 
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
