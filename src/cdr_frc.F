      module cdr_frc

      ! CDR forcing module
      ! --------------------

      ! for cdrs:  Qbar(nx,y), the total cdr flux for that grid point ! (m3/s) function of time
      !            Qshape(ncdrs), the vertical shape of the cdr distribution ! fractions sum(Qshape) = 1

      ! Coding by Jeroen Molemaker (March 2025)

#include "cppdefs.opt"
  
      use dimensions                                       ! has dimensions, a list of forcing files, and rmask
      use roms_read_write
      use nc_read_write
      use scalars
      use grid
      use ocean_vars

      implicit none

      private

#include "cdr_frc.opt"

      ! Variables used in the evolution equations
      integer,public  :: icdr,cidx ! Pipe index for looping through cdrs
      integer,public,allocatable,dimension(:,:) :: cdr_idx      ! cdr indices at grid points
      real   ,public,allocatable,dimension(:,:) :: cdr_prf      ! cdr fractions
      real   ,public,allocatable,dimension(:)   :: cdr_icdr     ! cdr release number

      real,public,dimension(ncdr)           :: cdr_vol      ! cdr volume
      real,public,dimension(ncdr,nt)        :: cdr_flx      ! cdr tracer flux

      real,dimension(ncdr,nt)        :: cdr_trc      ! cdr tracer concentration
      real,dimension(ncdr)           :: cdr_lon,cdr_lat         ! cdr location
      real,dimension(ncdr)           :: cdr_hsc,cdr_vsc,cdr_dep ! cdr scales and depth

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

        call set_ana_cdr_frc
        
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
      implicit none

      ! local
      integer :: ierr,ncid,i,j,k
      real :: norm,vint,arg
      real,dimension(:,:),allocatable :: dist,frac

      call init_arrays_cdr

      if (cdr_analytical) then

        ! cdr_flx is defined in ana_cdr_frc.h

      else

        ierr=nf90_open(cdr_file, nf90_nowrite, ncid)
        if(ierr/=0) 
     &      call handle_ierr(ierr,'init_cdr_frc:: Cant open cdr forcing file')

        call ncread(ncid,cdr_loc_lon,cdr_lon)
        call ncread(ncid,cdr_loc_lon,cdr_lat)
        call ncread(ncid,cdr_loc_dep,cdr_dep)
        call ncread(ncid,cdr_scl_hor,cdr_hsc)
        call ncread(ncid,cdr_scl_vrt,cdr_vsc)
        ierr = nf90_close(ncid)

        ! find fractional cdr fluxes for this subdomain
        allocate(dist(GLOBAL_2D_ARRAY))
        allocate(frac(GLOBAL_2D_ARRAY))
        do icdr= 1,ncdr
          norm = pi*cdr_hsc(icdr) ! 2d integral of exp(-(dist/hscl)^2 )
          call ll2dist(lonr,latr,cdr_lon(icdr),cdr_lat(icdr),dist);
          frac = exp(-(dist/cdr_hsc(icdr))**2)
          frac = frac/pm/pn/norm

          allocate(cdr_idx(GLOBAL_2D_ARRAY))
          do j=1,ny
            do i=1,nx
              if (frac(i,j) >1e-3) then
                cidx = cidx+1
                cdr_idx(i,j) = cidx
              endif
            enddo
          enddo
        enddo
        ! This presents a problem if 2 releases overlap in space

        allocate(cdr_prf(cidx,nz))
        allocate(cdr_icdr(cidx))
        do icdr= 1,ncdr
          do j=1,ny
            do i=1,nx
               if (cdr_idx(i,j)>0) then
                 cidx = cdr_idx(i,j)

                 ! figure vertical profile, depths are relative to surface
                 ! prf = exp( - ((z-d)/vsc )^2)
                 vint = 0
                 do k=1,nz
                   arg = ( (z_r(i,j,k) + cdr_dep(icdr) )/cdr_vsc(icdr) )**2
                   cdr_prf(cidx,k) = exp(-arg)*Hz(i,j,k)
                   vint = vint + cdr_prf(cidx,k)
                 enddo
                 cdr_prf(cidx,:) = frac(i,j)*cdr_prf(cidx,:)/vint
                 cdr_icdr(cidx) = icdr  ! identifies a specific release with this profile
               endif
            enddo
          enddo

          ! Todo, use mpi gather to obtain discretely unit integrals of fractions across sub-domains

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
      allocate( cdr_idx(GLOBAL_2D_ARRAY) );      cdr_idx=0.

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
!!       subroutine set_ana_cdr_frc  ![
      ! Analytical cdr forcing volume and tracer data
      ! Put here to avoid circular reference if in analytical.F

!! #include "ana_cdr_frc.h"

!     end subroutine set_ana_cdr_frc  !]
! ----------------------------------------------------------------------
      subroutine ll2dist(lon,lat,lon0,lat0,dist)  ![
      ! Great circle distance from lon,lat fields to point (lon0,lat0)
      ! Using haversine function
      implicit none

      ! import/export
      real,dimension(GLOBAL_2D_ARRAY),intent(in) :: lon,lat
      real,                           intent(in) :: lon0,lat0
      real,dimension(GLOBAL_2D_ARRAY),intent(out):: dist

      ! local
      real,parameter ::radius=6.371e6; ! radius earth in m
      real,dimension(:,:),allocatable :: a,dLon2,dLat2

      allocate(dLon2(GLOBAL_2D_ARRAY))
      allocate(dLat2(GLOBAL_2D_ARRAY))
      dLon2 = 0.5*(lon-lon0)
      dLat2 = 0.5*(lat-lat0)
      dist = sin(dLat2)**2 + cos(lat)*cos(lat0) * sin(dLon2)**2
      dist=2*atan2(sqrt(dist),sqrt(1-dist));

      dist = radius*dist;    !Haversine distance
      deallocate(dLon2)
      deallocate(dLat2)

      end subroutine ll2dist  !]
! ----------------------------------------------------------------------

      end module cdr_frc
