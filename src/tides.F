      module tides

#include "cppdefs.opt"

      use param
      use scalars                                          ! For iic, ntstart, deg2rad
      use nc_read_write
      use roms_read_write
      use mpi_exchanges
      use grid

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

      character(len=30) :: string 

      if (.not.pot_tides.and..not.bry_tides) then
         return
      endif
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

      if (pot_tides.or.bry_tides) then
        if (ana_tides) then
          call make_tides
        else
          call read_tides
        endif
      endif

      ! Print user options (tides.opt) to netcdf attributes
      tide_opt = ''
      write(string, '(A, I2)') '# of constituents:', ntides
      call store_string_att(tide_opt, string)

      if (ana_tides) then
        call store_string_att(tide_opt, ', analytical')
      else
        call store_string_att(tide_opt, ', real')
      endif

      if (pot_tides) call store_string_att(tide_opt, ', tidalpotential')
      if (bry_tides) call store_string_att(tide_opt, 
     &                                     ', barotropictide@boundary') 

      tide_opt = trim(adjustl(tide_opt))


      end subroutine init_tides  !]
!-----------------------------------------------------------------------
      subroutine set_tides  ![
      ! computes current forcing for boundaries and surface tidal
      ! potential

      use hidden_mpi_vars
      use boundary
      use grid
      implicit none

      ! local
      integer :: tile, itide, i, j
      real :: omT

      if (.not.allocated(ftide)) then
        call init_tides
      endif

      do itide=1,ntides

        omT = ftide(itide)*(time+0.5*dt)

        if (bry_tides) then
          if (inode==0.and.obc_west) then                  ! Western BC
            if (ana_bry) then
              zeta_west(j0:j1) = 0
              ubar_west(j0:j1) = 0
              vbar_west(1:j1)  = 0
            endif
            do j=j0,j1
              zeta_west(j)=zeta_west(j)
     &        + ztide_r(0,j,itide)*cos(omT)- ztide_i(0,j,itide)*sin(omT)
            enddo
            do j=j0,j1
              ubar_west(j)=ubar_west(j)
     &        + utide_r(1,j,itide)*cos(omT)- utide_i(1,j,itide)*sin(omT)
            enddo
            do j=1,j1
              vbar_west(j)=vbar_west(j)
     &        + vtide_r(0,j,itide)*cos(omT)- vtide_i(0,j,itide)*sin(omT)
            enddo
          endif

          if (inode==np_xi-1.and.obc_east) then           ! Eastern BC
            if (ana_bry) then
              zeta_east(j0:j1) = 0
              ubar_east(j0:j1) = 0
              vbar_east(1:j1)  = 0
            endif
            do j=j0,j1
              zeta_east(j) =zeta_east(j)
     &        + ztide_r(nx+1,j,itide)*cos(omT)- ztide_i(nx+1,j,itide)*sin(omT)
            enddo
            do j=j0,j1
              ubar_east(j)=ubar_east(j)
     &        + utide_r(nx+1,j,itide)*cos(omT)- utide_i(nx+1,j,itide)*sin(omT)
            enddo
            do j=1,j1
              vbar_east(j)=vbar_east(j)
     &        + vtide_r(nx+1,j,itide)*cos(omT)- vtide_i(nx+1,j,itide)*sin(omT)
            enddo
          endif

          if (jnode==0.and.obc_south) then                 ! Southern BC
            if (ana_bry) then
              zeta_south(i0:i1) = 0
              ubar_south(1:i1)  = 0
              vbar_south(i0:i1) = 0
            endif
            do i=i0,i1
              zeta_south(i)=zeta_south(i)
     &        + ztide_r(i,0,itide)*cos(omT)- ztide_i(i,0,itide)*sin(omT)
            enddo
            do i=1,i1
              ubar_south(i)=ubar_south(i)
     &        + utide_r(i,0,itide)*cos(omT)- utide_i(i,0,itide)*sin(omT)
            enddo
            do i=i0,i1
              vbar_south(i)=vbar_south(i)
     &        + vtide_r(i,1,itide)*cos(omT)- vtide_i(i,1,itide)*sin(omT)
            enddo
          endif

          if (jnode==np_eta-1.and.obc_north) then          ! Northern BC
            if (ana_bry) then
              zeta_north(i0:i1) = 0
              ubar_north(i0:i1) = 0
              vbar_north(i0:i1) = 0
            endif
            do i=i0,i1
              zeta_north(i)=zeta_north(i)
     &        + ztide_r(i,ny+1,itide)*cos(omT)- ztide_i(i,ny+1,itide)*sin(omT)
            enddo
            do i=1,i1
              ubar_north(i)=ubar_north(i)
     &        + utide_r(i,ny+1,itide)*cos(omT)- utide_i(i,ny+1,itide)*sin(omT)
            enddo
            do i=i0,i1
              vbar_north(i)=vbar_north(i)
     &        + vtide_r(i,ny+1,itide)*cos(omT)- vtide_i(i,ny+1,itide)*sin(omT)
            enddo
          endif
        endif                                              ! <-- end bry_tides

        if (pot_tides) then
!-----------------------------------------------------------------------
!  Compute tidal potential (m)
!  -- > to be applied as a pressure-gradient force
!       in prsgrd (if SOLVE3D) or step2d (if not SOLVE3D)
!-----------------------------------------------------------------------
          if (itide==1) then                               ! initialize at the first tidal component
            do j=0,ny
              do i=0,nx
                ptide(i,j)=
     &        + ptide_r(i,j,itide)*cos(omT)- ptide_i(i,j,itide)*sin(omT)
              enddo
            enddo
          else
            do j=0,ny
              do i=0,nx
                ptide(i,j)= ptide(i,j)+
     &          + ptide_r(i,j,itide)*cos(omT)-ptide_i(i,j,itide)*sin(omT)
              enddo
            enddo
          endif
        endif                                              ! <-- end pot_tides

      enddo                                                ! <-- end itide=1,ntides

      end subroutine set_tides  !]
!-----------------------------------------------------------------------
      subroutine read_tides  ![
      ! read tides from forcing file. Only read once at first timestep.

      use netcdf
      use dimensions
      implicit none

      ! local
      integer :: i,j
      integer :: ierr, var_file_indx, ncid, var_id, itide, found_rec
      integer,dimension(3) :: start
      character(len=100) :: info
      character(len=300) :: frcstr


      found_rec = 0                                        ! set found variable flag to not found
      var_file_indx = 1                                    ! forcing file index number

      ! since tides are only read once, don't need to store netcdf vars
      ! continue searching while tides not found & still more files to check
      ! search for forcing file with tidal information
      do while ((found_rec==0).and.(var_file_indx <= max_frc))

        ierr=nf90_open(frcfile(var_file_indx), nf90_nowrite, ncid)
        if(ierr/=0) 
     &     call handle_ierr(ierr,'tides :: opening forcing file')

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
        if (ierr/=0) print *, 'no omega in file'
        ierr=nf90_get_var(ncid,var_id,ftide(itide),(/itide/))
        if(ierr/=0) call handle_ierr(ierr,'tides :: read omega')

        start(3) = itide

        if (bry_tides) then
          call ncread(ncid,'ssh_Re',ztide_r(i0:i1,j0:j1,itide),start)   ! read tidal elevations
          call ncread(ncid,'ssh_Im',ztide_i(i0:i1,j0:j1,itide),start)
!         call exchange_xxx(ztide_r(:,:,itide),ztide_i(:,:,itide))

          call ncread(ncid,'u_Re',utide_r(1:i1,j0:j1,itide),start)      ! read tidal barotropic velocities.
          call ncread(ncid,'u_Im',utide_i(1:i1,j0:j1,itide),start)      ! these currents must be on u,v-points
!         call exchange_xxx(utide_r(:,:,itide),utide_i(:,:,itide))      ! and rotated to align with the grid

          call ncread(ncid,'v_Re',  vtide_r(i0:i1,1:j1,itide),start)    
          call ncread(ncid,'v_Im',  vtide_i(i0:i1,1:j1,itide),start)
!         call exchange_xxx(vtide_r(:,:,itide),vtide_i(:,:,itide))
        endif

        if (pot_tides) then
          call ncread(ncid,'pot_Re',ptide_r(i0:i1,j0:j1,itide),start)
          call ncread(ncid,'pot_Im',ptide_i(i0:i1,j0:j1,itide),start)
!         call exchange_xxx(ptide_r(:,:,itide),ptide_i(:,:,itide))
        endif

      end do ! <-- itide=1,ntides

      ! add tide info to forcing attributes - embedded in here due to tides init
      if (mynode == 0) then
        ierr = nf90_get_att(ncid,nf90_global,'data_source',info)
        if (len(trim(info)) /= 0) then
          if (index(force_info,info) == 0) then
            frctype = 'TIDES:'
            write(frcstr,'(2A)') trim(adjustl(frctype)),
     &            trim(adjustl(info))
            call store_string_att(force_info,frcstr)
          endif
        endif
      endif

      if (bry_tides) then
        call exchange_xxx(ztide_r,ztide_i)
        call exchange_xxx(utide_r,utide_i)
        call exchange_xxx(vtide_r,vtide_i)
      endif
      if (pot_tides) then
        call exchange_xxx(ptide_r,ptide_i)
      endif

      if(mynode==0) write(*,'(7x,A)')'tides :: read all tides'

      end subroutine read_tides  !]
!-----------------------------------------------------------------------
      subroutine make_tides  ![
      ! Construct an analytical tidal forcing

      use dimensions
      implicit none

      ! local

#include "ana_tide.h"

      if(mynode==0) write(*,*) 'tides :: constructed analytical tides'

      end subroutine make_tides  !]
!-----------------------------------------------------------------------

      end module tides
