      module cdr_frc

      ! CDR forcing module
      ! --------------------

      ! for cdrs:  Qbar(nx,y), the total cdr flux for that grid point ! (m3/s) function of time
      !            Qshape(ncdrs), the vertical shape of the cdr distribution ! fractions sum(Qshape) = 1

      ! done: add MPI_reduce_scatter
      ! force all tracers or single tracer forcing

      ! Coding by Jeroen Molemaker (March 2025)

#include "cppdefs.opt"

      use dimensions
      use roms_read_write
      use nc_read_write
      use scalars
      use grid
      use ocean_vars
      use hidden_mpi_vars

      implicit none

      private

#include "cdr_frc.opt"

      ! Variables used in the evolution equations
      integer,public                            :: cdr_nprf     ! number of locations with a release
      real   ,public,allocatable,dimension(:,:) :: cdr_prf      ! cdr fractions
      integer,public,allocatable,dimension(:)   :: cdr_iloc     ! cdr local i location
      integer,public,allocatable,dimension(:)   :: cdr_jloc     ! cdr local j location
      real   ,public,allocatable,dimension(:)   :: cdr_icdr     ! cdr release number

      integer,public  :: icdr,cidx !   index for looping through cdr locations

      real,public,dimension(ncdr)           :: cdr_vol      ! cdr volume
      real,public,dimension(ncdr,nt)        :: cdr_flx      ! cdr tracer flux

      real,dimension(ncdr,nt)        :: cdr_trc      ! cdr tracer concentration
      real,dimension(ncdr)           :: cdr_lon,cdr_lat         ! cdr location
      real,dimension(ncdr)           :: cdr_hsc,cdr_vsc,cdr_dep ! cdr scales and depth
      integer,dimension(0:ncdr)   :: cdr_nloc=0  ! number of active grid points for per cdr

     ! Netcdf names
      character(len=9)  :: cdr_flx_name = 'cdr_flux'  !! stored in the grid file
      character(len=11) :: cdr_vol_name = 'cdr_volume'!! stored in a forcing file
      character(len=11) :: cdr_trc_name = 'cdr_tracer'!! stored in a forcing file
      character(len=9)  :: cdr_tim_name = 'cdr_time'  !! stored in a forcing file
!     character(len=5)  :: ncdr_dim_name = 'ncdr'      !! dimension name for number of cdrs in file
!     character(len=8)  :: ntrc_dim_name = 'ntracers'   !! dimension name for number of tracers in file


      logical :: init_cdr_done = .false.

      public set_cdr_frc
      public init_cdr_frc

      contains

! ----------------------------------------------------------------------
      subroutine set_cdr_frc  ![
      ! set cdr forces (realistic and analytical)
      ! - read and interpolation all cdr forcing.
      ! - all cdr variables need time interpolation only
      !   here so can use same generic routine.
      ! - input data in days!
      implicit none

      ! local
      integer :: i,j

      if (.not. init_cdr_done) call init_cdr_frc

      if (cdr_analytical) then                                ! Set cdr flux volumes and tracer data

!       call set_ana_cdr_frc

      else

        if (cdr_volume) then
          call set_frc_data(nc_cdrvol,cdr_vol)                ! set cdr volume for all cdrs at current time
          call set_frc_data(nc_cdrtrc,var2d=cdr_trc)          ! set cdr tracers conc. for all cdrs at current time

          do icdr=1,ncdr
            cdr_flx(icdr,:) = cdr_vol(icdr)*cdr_trc(icdr,:)
          enddo

        else
          call set_frc_data(nc_cdrflx,var2d=cdr_flx)          ! set cdr tracer flux for all cdrs at current time

        endif

      endif

      end subroutine set_cdr_frc !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc  ![

      ! Initialize cdr forcing:
      ! Read in a grid file with locations of cdrs and flux contribution per cell.
      ! Done only once as cdr position does not change.
      ! Realistic case - stored as one value in NetCDF file where
      ! cdr grid point value = pidx + cdr_fraction

      use netcdf
      use mpi
      implicit none

      ! local
      integer :: ierr,ncid,i,j,k
      real :: norm,vint,arg
      real :: local_int,global_int
      real, dimension(nz) :: arg_nz
      real :: local_min_val
      integer, dimension(2) :: lmi
      real, dimension(2) :: local_minloc
      real, dimension(2) :: global_minloc
      real,dimension(:,:)  ,allocatable :: dist
      real,dimension(:,:,:),allocatable :: frac

      call init_arrays_cdr

      if (cdr_analytical) then

        ! cdr_flx is defined in ana_cdr_frc.h

      else

        ierr=nf90_open(cdr_file, nf90_nowrite, ncid)
        if(ierr/=0)
     &      call handle_ierr(ierr,'init_cdr_frc:: Cant open cdr forcing file')

        call ncread(ncid,cdr_loc_lon,cdr_lon)
        call ncread(ncid,cdr_loc_lat,cdr_lat)
        call ncread(ncid,cdr_loc_dep,cdr_dep)
        call ncread(ncid,cdr_scl_hor,cdr_hsc)
        call ncread(ncid,cdr_scl_vrt,cdr_vsc)
        ierr = nf90_close(ncid)

        ! find fractional cdr fluxes for this subdomain
        allocate(dist(GLOBAL_2D_ARRAY)); dist = 1e6 ! Some large number so we don't get single-point releases in the halos
        allocate(frac(GLOBAL_2D_ARRAY,ncdr))

        cidx = 0
        do icdr= 1,ncdr

          if (cdr_hsc(icdr) == 0) then
            call ll2dist(lonr,latr,cdr_lon(icdr),cdr_lat(icdr),dist);

            local_min_val = MINVAL(dist)
            ! Find local index that is closest to release location
            lmi = MINLOC(dist)
            ! Adjust index to account for the fact that GLOBAL_ARRAY does not start at index 1
            lmi(:) = lmi(:) + lbound(dist) - 1

            ! Pack local minimum and MPI rank into a vector
            local_minloc(1) = local_min_val
            local_minloc(2) = mynode

            ! Find global minimum and rank where the minimum occurs
            call MPI_Reduce( local_minloc, global_minloc, 2, mpi_2double_precision,
     &      mpi_minloc, 0, ocean_grid_comm, ierr)

            ! Broadcast rank where minimum occurs to all processes
            call MPI_Bcast(global_minloc, 2,mpi_double_precision,0,ocean_grid_comm,ierr)

            if (mynode == global_minloc(2)) then
              print *, 'The minimum distance to single-point Release', icdr, 'is', global_minloc(1)
              print *, 'This is on rank:', mynode
              print *, 'at point', lmi
              print *, 'The intended release location was Lon:', cdr_lon(icdr), 'Lat:', cdr_lat(icdr)
              print *, 'The release will take place at Lon:', lonr(lmi(1),lmi(2)), 'Lat:', latr(lmi(1),lmi(2))
              if (rmask(lmi(1),lmi(2))==0) then
                error stop 'ERROR: single-point CDR release requested on a land point'
              else
                frac(lmi(1),lmi(2),icdr) = 1
                cidx = cidx+1
              endif
            endif

          else

!           close  = max(3*cdr_hsc(icdr),dx)
            norm = pi*cdr_hsc(icdr)*cdr_hsc(icdr) ! 2d integral of exp(-(dist/hscl)^2 )
            call ll2dist(lonr,latr,cdr_lon(icdr),cdr_lat(icdr),dist);
            frac(:,:,icdr) = exp(-(dist/cdr_hsc(icdr))**2)
            frac(:,:,icdr) = frac(:,:,icdr)/pm/pn/norm
            do j=1,ny
              do i=1,nx
!               if (dist(i,j)<close.and.rmask(i,j)>0) then
                if (frac(i,j,icdr) >1e-3.and.rmask(i,j)>0) then
                  cidx = cidx+1
                endif
              enddo
            enddo
          endif ! hscl = 0
        enddo ! icdr

        cdr_nprf = cidx

        ! allocate total number of profiles for all releases
        allocate(cdr_prf(cdr_nprf,nz))
        cdr_prf(:,:) = 0
        
        allocate(cdr_icdr(cdr_nprf))
        allocate(cdr_iloc(cdr_nprf))
        allocate(cdr_jloc(cdr_nprf))

        ! populate the profiles with fractional releases
        cidx = 0
        do icdr= 1,ncdr
          do j=1,ny
            do i=1,nx
              if (frac(i,j,icdr) >1e-3.and.rmask(i,j)>0) then
                 cidx = cidx+1

                 ! Location and specific release for this profile
                 cdr_icdr(cidx) = icdr
                 cdr_iloc(cidx) = i
                 cdr_jloc(cidx) = j

                 ! figure vertical profile, depths are relative to surface
                 ! prf = exp( - ((z-d)/vsc )^2)
                 if (cdr_vsc(icdr) == 0) then
                   do k=1,nz
                     arg_nz(k) = abs(z_r(i,j,k) + cdr_dep(icdr))
                   enddo
                     arg = minloc(arg_nz,1)
                     cdr_prf(cidx,arg) = frac(i,j,icdr)
                 else
                   vint = 0
                   do k=1,nz
                     arg = ( (z_r(i,j,k) + cdr_dep(icdr) )/cdr_vsc(icdr) )**2
                     cdr_prf(cidx,k) = exp(-arg)*Hz(i,j,k)
                     vint = vint + cdr_prf(cidx,k)
                   enddo
                   cdr_prf(cidx,:) = frac(i,j,icdr)*cdr_prf(cidx,:)/vint
                 endif
               endif
            enddo
          enddo
          cdr_nloc(icdr) = cidx

          local_int = sum(cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),:))
          call MPI_Reduce(local_int,global_int,1,
     %           mpi_double_precision,mpi_sum,0,ocean_grid_comm,ierr)
          call MPI_Bcast(global_int,1,mpi_double_precision,0,ocean_grid_comm,ierr)

          cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),:) =
     %        cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),:)/global_int
          local_int = sum(cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),:))

        enddo

        deallocate(dist)
        deallocate(frac)

      endif

      init_cdr_done = .true.
      if(mynode==0) write(*,'(/7x,A/)') 'cdr_frc: init cdr locations'

      end subroutine init_cdr_frc  !]
! ----------------------------------------------------------------------
      subroutine init_arrays_cdr  ![
      implicit none

      character(len=30) :: string

      if (.not. cdr_analytical) then
        if (cdr_volume) then
          allocate(nc_cdrvol%vdata(ncdr,1 ,2))
          allocate(nc_cdrtrc%vdata(ncdr,nt,2))
        else
          allocate(nc_cdrflx%vdata(ncdr,nt,2))
        endif
      endif

      ! Print user options (cdr_frc.opt) to netcdf attributes
      cdr_frc_opt = ''
      write(string, '(A,I3)') 'ncdr =', ncdr
      call store_string_att(cdr_frc_opt, string)
      if (cdr_analytical) then
        call store_string_att(cdr_frc_opt, ', Analytical')
      endif

      end subroutine init_arrays_cdr  !]
! ----------------------------------------------------------------------
!     subroutine set_ana_cdr_frc  ![
      ! Analytical cdr forcing volume and tracer data
      ! Put here to avoid circular reference if in analytical.F

!!!!!#include "ana_cdr_frc.h"

!     end subroutine set_ana_cdr_frc  !]
! ----------------------------------------------------------------------
      subroutine ll2dist(lon,lat,lon0,lat0,dist)  ![
      ! Great circle distance from lon,lat fields to point (lon0,lat0)
      ! Using haversine function
      implicit none

      ! import/export
      integer :: i,j
      real,dimension(GLOBAL_2D_ARRAY),intent(in) :: lon,lat
      real,                           intent(in) :: lon0,lat0
      real,dimension(GLOBAL_2D_ARRAY),intent(out):: dist

      ! local
      real,parameter ::radius=6.371e6; ! radius earth in m
      real,dimension(:,:),allocatable :: a,dLon2,dLat2

      allocate(dLon2(GLOBAL_2D_ARRAY))
      allocate(dLat2(GLOBAL_2D_ARRAY))
      dLon2 = 0.5*(lon-lon0)*deg2rad
      dLat2 = 0.5*(lat-lat0)*deg2rad

      ! Need to loop over only the points in the computational domain (i.e. not the halos)
      ! This will prevent single-point releases from occurring in the halos, in which case
      ! they aren't actually included.
      do j=1,ny
        do i=1,nx
          dist(i,j) = sin(dLat2(i,j))**2 + cos(lat(i,j)*deg2rad)*cos(lat0*deg2rad) * sin(dLon2(i,j))**2
          dist(i,j)=2*atan2(sqrt(dist(i,j)),sqrt(1-dist(i,j)));
          dist(i,j) = radius*dist(i,j);    !Haversine distance
        enddo
      enddo
      deallocate(dLon2)
      deallocate(dLat2)

      end subroutine ll2dist  !]
! ----------------------------------------------------------------------

      end module cdr_frc
