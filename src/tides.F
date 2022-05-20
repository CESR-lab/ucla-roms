      module tides

#include "cppdefs.opt"
! Empty module if TIDES not defined.
! JM The code is currently broken for pot_tides=.false.
! JM Broken in step2d and prsgrd
#ifdef TIDES

      use param
      use scalars                                          ! For iic, ntstart, deg2rad
      use nc_read_write
      use roms_read_write

      implicit none

      private                                              ! default variables private to module

#include "tides.opt"

      real, dimension(:)    ,allocatable :: ftide
      real, dimension(:,:,:),allocatable :: ztide_r,ztide_i
      real, dimension(:,:,:),allocatable :: utide_r,utide_i
      real, dimension(:,:,:),allocatable :: vtide_r,vtide_i
      real, dimension(:,:,:),allocatable :: ptide_r,ptide_i
      real, public,dimension(:,:),  allocatable :: ptide

      public set_tides

      contains

!-----------------------------------------------------------------------
      subroutine init_tides  ![
      ! allocate tidal variables and read tidal data

      implicit none

      allocate(ftide(ntides))

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

      end subroutine init_tides  !]

!-----------------------------------------------------------------------
      subroutine set_tides(tile)  ![
      ! computes current forcing for boundaries and surface tidal
      ! potential
      ! DevinD created seperate step as can't have compute_tile_bounds
      ! and compute_auxilliary_bounds in the same subroutine as they both
      ! declare and compute. Fortan doesn't allow.

      use hidden_mpi_vars
      use private_scratch  ! this isn't needed?

      implicit none

      integer :: tile

#include "compute_tile_bounds.h"

      call set_tides_tile(istr,iend,jstr,jend)

      end subroutine set_tides  !]

!-----------------------------------------------------------------------
      subroutine set_tides_tile(istr,iend,jstr,jend)  ![

      use boundary
      use grid
      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer :: tile, itide, i, j
      real :: omT

#include "compute_auxiliary_bounds.h"

      if (.not.allocated(ftide)) then 
        call init_tides
      endif

      do itide=1,ntides

        omT = ftide(itide)*time

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

      end subroutine set_tides_tile  !]

!-----------------------------------------------------------------------
      subroutine read_tides  ![
      ! read tides from forcing file. Only read once at first timestep.

      use netcdf
      implicit none

      ! local
      integer :: i,j
      integer :: ierr, var_file_indx, ncid, var_id, itide, found_rec
      integer,dimension(3) :: start


      found_rec = 0                                        ! set found variable flag to not found
      var_file_indx = 1                                    ! forcing file index number

      ! since tides are only read once, don't need to store netcdf vars
      ! continue searching while tides not found & still more files to check
      ! search for forcing file with tidal information
      do while ((found_rec==0).and.(var_file_indx <= max_frc))

        ierr=nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)
        if(ierr/=0) call handle_ierr(ierr,'tides :: opening forcing file')

        ierr=nf90_inq_varid(ncid, 'omega', var_id)
        if(ierr/=0) then 
          var_file_indx = var_file_indx +1
        else
          found_rec=1 ! Set to 1 to exit do while loop
        endif
      enddo 

      if (found_rec==0) then
        error stop 'tides :: tide input data not found'
      endif

      start = 1                                                      ! netcdf start indices

      do itide = 1,ntides

        ierr=nf90_inq_varid(ncid, 'omega', var_id)
        ierr=nf90_get_var(ncid,var_id,ftide(itide),(/itide/))
        if(ierr/=0) call handle_ierr(ierr,'tides :: read omega')

        start(3) = itide

        if (bry_tides) then

          call ncread(ncid,'ssh_Re',ztide_r(i0:i1,j0:j1,itide),start)          ! read tidal elevations
          call ncread(ncid,'ssh_Im',ztide_i(i0:i1,j0:j1,itide),start)

          call ncread(ncid,'u_Re',  utide_r( 1:i1,j0:j1,itide),start)          ! read tidal barotropic velocities.
          call ncread(ncid,'u_Im',  utide_i( 1:i1,j0:j1,itide),start)          ! these currents must be on u,v-points
          call ncread(ncid,'v_Re',  vtide_r(i0:i1, 1:j1,itide),start)          ! and rotated to align with the grid
          call ncread(ncid,'v_Im',  vtide_i(i0:i1, 1:j1,itide),start)

        endif

        if (pot_tides) then

          call ncread(ncid,'pot_Re',ptide_r(i0:i1,j0:j1,itide),start)
          call ncread(ncid,'pot_Im',ptide_i(i0:i1,j0:j1,itide),start)

        endif

      end do ! <-- itide=1,ntides

      if(mynode==0) write(*,'(7x,A)')'tides :: read all tides'

      end subroutine read_tides  !]

!-----------------------------------------------------------------------

#endif /* TIDES for whole module */

      end module tides
