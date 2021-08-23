      module tides

#include "cppdefs.h"
! Empty module if TIDES not defined.
#ifdef TIDES

      implicit none


      ! Make all variable private to this module unless public specified
      private

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_Array to work.
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


      real, public,dimension(:)    ,allocatable :: omega
      real, public,dimension(:,:,:),allocatable :: ztide_r,ztide_i
      real, public,dimension(:,:,:),allocatable :: utide_r,utide_i
      real, public,dimension(:,:,:),allocatable :: vtide_r,vtide_i
      real, public,dimension(:,:,:),allocatable :: ptide_r,ptide_i
      real, public,dimension(:,:),  allocatable :: ptide

      public set_tides

      contains
!-----------------------------------------------------------------------
      subroutine init_tides
      ! Allocate tidal variables and read tidal data

      implicit none

      allocate(omega(ntides))

      if (bry_tides) then
        allocate(ztide_r(GLOBAL_2D_ARRAY,ntides))
        allocate(ztide_i(GLOBAL_2D_ARRAY,ntides))
        allocate(utide_r(GLOBAL_2D_ARRAY,ntides))
        allocate(utide_i(GLOBAL_2D_ARRAY,ntides))
        allocate(vtide_r(GLOBAL_2D_ARRAY,ntides))
        allocate(vtide_i(GLOBAL_2D_ARRAY,ntides))
      endif

      if (pot_tides) then
        allocate(ptide_r(GLOBAL_2D_ARRAY,ntides))
        allocate(ptide_i(GLOBAL_2D_ARRAY,ntides))
        allocate(ptide(GLOBAL_2D_ARRAY))
      endif

      call read_tides

      end subroutine init_tides

!-----------------------------------------------------------------------
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
      real omT

#include "compute_auxiliary_bounds.h"

      if (.not.allocated(omega)) then 
        call init_tides
      endif

      do itide=1,ntides

        omT = omega(itide)*time

        if (bry_tides) then

#  ifdef OBC_WEST
          if (WESTERN_EDGE) THEN
            i=istr-1
            do j=jstrR,jendR
              zeta_west(j) =zeta_west(j)
     &        + ztide_r(i,j,itide)*cos(omT)- ztide_i(i,j,itide)*sin(omT)
            enddo
            i=istrU-1
            do j=jstrR,jendR
              ubar_west(j)=ubar_west(j)
     &        + utide_r(i,j,itide)*cos(omT)- utide_i(i,j,itide)*sin(omT)
            enddo
            i=istr-1
            do j=jstr,jendR
              vbar_west(j)=vbar_west(j)
     &        + vtide_r(i,j,itide)*cos(omT)- vtide_i(i,j,itide)*sin(omT)
            enddo
          endif
#  endif

#  ifdef OBC_EAST
          if (EASTERN_EDGE) THEN
            i = iend+1
            do j=jstrR,jendR
              zeta_east(j) =zeta_east(j)
     &        + ztide_r(i,j,itide)*cos(omT)- ztide_i(i,j,itide)*sin(omT)
            enddo
            do j=jstrR,jendR
              ubar_east(j)=ubar_east(j)
     &        + utide_r(i,j,itide)*cos(omT)- utide_i(i,j,itide)*sin(omT)
            enddo
            do j=jstr,jendR
              vbar_east(j)=vbar_east(j)
     &        + vtide_r(i,j,itide)*cos(omT)- vtide_i(i,j,itide)*sin(omT)
            enddo
          endif
#  endif

#  ifdef OBC_SOUTH
          if (SOUTHERN_EDGE) THEN
            j=jstr-1
            do i=istrR,iendR
              zeta_south(i)=zeta_south(i)
     &        + ztide_r(i,j,itide)*cos(omT)- ztide_i(i,j,itide)*sin(omT)
            enddo
            j=jstr-1
            do i=istr,iendR
              ubar_south(i)=ubar_south(i)
     &        + utide_r(i,j,itide)*cos(omT)- utide_i(i,j,itide)*sin(omT)
            enddo
            j=jstrV-1
            do i=istrR,iendR
              vbar_south(i)=vbar_south(i)
     &        + vtide_r(i,j,itide)*cos(omT)- vtide_i(i,j,itide)*sin(omT)
            enddo
          endif
#  endif
#  ifdef OBC_NORTH
          if (NORTHERN_EDGE) THEN
            j=jend+1
            do i=istrR,iendR
              zeta_north(i)=zeta_north(i)
     &        + ztide_r(i,j,itide)*cos(omT)- ztide_i(i,j,itide)*sin(omT)
            enddo
            do i=istr,iendR
              ubar_north(i)=ubar_north(i)
     &        + utide_r(i,j,itide)*cos(omT)- utide_i(i,j,itide)*sin(omT)
            enddo
            do i=istrR,iendR
              vbar_north(i)=vbar_north(i)
     &        + vtide_r(i,j,itide)*cos(omT)- vtide_i(i,j,itide)*sin(omT)
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
                ptide(i,j)=
     &        + ptide_r(i,j,itide)*cos(omT)- ptide_i(i,j,itide)*sin(omT)
              enddo
            enddo
          else
            do j=jstrR-1,jendR ! DevinD does jstrR-1 not go beyond computational margin?
              do i=istrR-1,iendR
                ptide(i,j)= ptide(i,j)+
     &          + ptide_r(i,j,itide)*cos(omT)-ptide_i(i,j,itide)*sin(omT)
              enddo
            enddo
          endif

        endif  ! bry_tides

      enddo  !  itide=1,ntides

      end subroutine set_tides_tile
!-----------------------------------------------------------------------

      subroutine read_tides
      ! Read tides from forcing file. Only read once at first timestep.

      use read_write
      use netcdf

      implicit none

      ! For max_frc
#include "ncvars.h"
      ! scalars.h for deg2rad

      ! local
      integer i,j
      integer ierr, var_file_indx, ncid, var_id, itide, found_rec


      found_rec = 0     ! Set found variable flag to not found
      var_file_indx = 1 ! Forcing file index number
      ! Since tides are only read once, don't need to store netcdf vars

      ! Continue searching while tides not found & still more files to check

      !! search for forcing file with tidal information
      do while ((found_rec==0).and.(var_file_indx <= max_frc))

        ierr=nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)
        if(ierr/=0) call handle_ierr(ierr,'tides::opening forcing file')

        ierr=nf90_inq_varid(ncid, 'omega', var_id)
        if(ierr/=0) then 
          var_file_indx = var_file_indx +1
        else
          found_rec=1 ! Set to 1 to exit do while loop
        endif
      enddo 

      if (found_rec==0) then
        error stop 'tides.F: tide input data not found'
      endif

      do itide = 1,ntides

        ierr=nf90_inq_varid(ncid, 'omega', var_id)
        ierr=nf90_get_var(ncid,var_id,omega(itide),(/itide/))
        if(ierr/=0) call handle_ierr(ierr,'tides::omega')

        if (bry_tides) then

          ! read tidal elevations
          call nc_read_var(ncid,ztide_r(:,:,itide),1,'ssh_Re',     
     &                 r_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::ssh_Re')

          call nc_read_var(ncid,ztide_i(:,:,itide),1,'ssh_Im',     
     &                 r_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::ssh_Im')


          ! Read tidal barotropic velocities. These currents must be on u,v-points
          ! and rotated to align with the grid

          call nc_read_var(ncid,utide_r(:,:,itide),1,'u_Re',     
     &                 u_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::u_Re')
          call nc_read_var(ncid,utide_i(:,:,itide),1,'u_Im',     
     &                 u_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::u_Im')

          call nc_read_var(ncid,vtide_r(:,:,itide),1,'v_Re',     
     &                 v_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::v_Re')
          call nc_read_var(ncid,vtide_i(:,:,itide),1,'v_Im',     
     &                 v_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::v_Im')


        endif !bry_tides

        if (pot_tides) then

          call nc_read_var(ncid,ptide_r(:,:,itide),1,'pot_Re',     
     &                 r_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::pot_Re')
          call nc_read_var(ncid,ptide_i(:,:,itide),1,'pot_Im',     
     &                 r_var, itide, ierr)
          if(ierr/=0) call handle_ierr(ierr,'tides::pot_Im')

        endif !pot_tides

      end do ! itide=1,ntides

      if(mynode==0) write(*,'(5x,A)')'tides:: read all tides'

      end subroutine read_tides
!-----------------------------------------------------------------------

#endif /* TIDES for whole module */

      end module tides
