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

      implicit none

      private

#include "cdr_frc.opt"

      ! Variables used in the evolution equations
      integer,public  :: pidx ! Pipe index for looping through cdrs
      integer,public,allocatable,dimension(:,:) :: cdr_idx      ! cdr indices at grid points
      real   ,public,allocatable,dimension(:,:) :: cdr_fraction ! cdr fractional flux at grid points
      real   ,public,allocatable,dimension(:,:) :: cdr_flx      ! cdr flux
      real   ,public, dimension(ncdr,N)         :: cdr_prf      ! Pipe vertical profile
      real,public,allocatable,dimension(:,:)    :: cdr_trc      ! Pipe tracer conc.
      real   ,public, dimension(ncdr)           :: cdr_vol      ! Pipe volume

     ! Netcdf names
      character(len=9)  :: cdr_flx_name = 'cdr_flux'  !! stored in the grid file
      character(len=11) :: cdr_vol_name = 'cdr_volume'!! stored in a forcing file
      character(len=11) :: cdr_trc_name = 'cdr_tracer'!! stored in a forcing file
      character(len=9)  :: cdr_tim_name = 'cdr_time'  !! stored in a forcing file
      character(len=5)  ::ncdr_dim_name = 'ncdr'      !! dimension name for number of cdrs in file
      character(len=8)  :: ntrc_dim_name = 'ntracers'   !! dimension name for number of tracers in file


      logical :: init_cdr_done = .false.

      public set_cdr_frc
      public init_arrays_cdrs

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

      if (p_analytical) then                                ! Set cdr flux volumes and tracer data

        call set_ana_cdr_frc
        
      else

        call set_frc_data(nc_pvol,cdr_vol)                ! set cdr volume for all cdrs at current time

        call set_frc_data(nc_ptrc,var2d=cdr_trc)          ! set cdr tracers conc. for all cdrs at current time

        call set_cdr_vert_prf                             ! set cdr vertical profiles

        do j = 1,ny                                        ! calculate fluxes with time interpolated volumes and conc.:
          do i = 1,nx
            if (cdr_idx(i,j) > 0.) then
              cdr_flx(i,j)=cdr_fraction(i,j)*cdr_vol( cdr_idx(i,j) )
            endif
          enddo
        enddo

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
      integer :: ierr,ncid,v_id,i,j


      if (p_analytical) then

        ! cdr_flx is defined in ana_cdr_frc.h

      else

        ! Read 'cdr_lon,cdr_lat' from crd file.
        allocate(cdr_lon(ncdr))
        allocate(cdr_lat(ncdr))
        allocate(cdr_hsc(ncdr))
        allocate(cdr_vsc(ncdr))
        allocate(cdr_dep(ncdr))

        ierr=nf90_open(cdrname, nf90_nowrite, ncid)
        if(ierr/=0) 
     &      call handle_ierr(ierr,'init_cdr_frc:: Cant open cdr forcing file')

        call ncread(ncid,cdr_loc_lon,cdr_lon)
        call ncread(ncid,cdr_loc_lon,cdr_lat)
        call ncread(ncid,cdr_loc_hsc,cdr_hsc)
        call ncread(ncid,cdr_loc_vsc,cdr_vsc)
        call ncread(ncid,cdr_loc_dep,cdr_dep)
        ierr = nf90_close(cdrname)

        ! find fractional cdr fluxes for this subdomain
        allocate(dist(GLOBAL_2D_ARRAY))
        allocate(frac(GLOBAL_2D_ARRAY))
        do icdr= 1,ncdr
          norm = pi*hsc*cdr_hsc(icdr) ! 2d integral of exp(-(dist/hscl)^2 )
          call ll2dist(lonr,latr,cdr_lon(icdr),cdr_lat(icdr),dist);
          frac = exp(-(dist/cdr_hsc)**2)
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

        allocate(cdr_prf(cidx,nz))
        do icdr= 1,ncdr
          do j=1,ny
            do i=1,nx
               if (cdr_idx(i,j)>0) then
                 cidx = cdr_idx(i,j)

                 ! figure vertical profile, depths are relative to surface
                 ! prf = exp( - ((z-d)/vsc )^2)
                 vint = 0
                 do k=1,nz
                   arg = ( (z_r(i,j,k) + cdr_dep(icdr) )/cdr_vsc)**2
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
      subroutine init_arrays_cdrs  ![
      implicit none

      character(len=30) :: string  
      allocate( cdr_idx(GLOBAL_2D_ARRAY) );      cdr_idx=0.
      allocate( cdr_fraction(GLOBAL_2D_ARRAY) ); cdr_fraction=0.
      allocate( cdr_flx(GLOBAL_2D_ARRAY) );      cdr_flx=0.
      allocate( cdr_trc(ncdr,nt) ) ;             cdr_trc=0.

      if (.not. cdr_analytical) then
        allocate(nc_cdrvol%vdata(ncdr,1 ,2))
        allocate(nc_cdrtrc%vdata(ncdr,nt,2))
      endif

      ! Print user options (cdr_frc.opt) to netcdf attributes 
      cdr_frc_opt = ''
      write(string, '(A,I3)') 'ncdr =', ncdr
      call store_string_att(cdr_frc_opt, string)
      if (p_analytical) then
        call store_string_att(cdr_frc_opt, ', Analytical')
      endif

      end subroutine init_arrays_cdrs  !]
! ----------------------------------------------------------------------
      subroutine set_cdr_vert_prf  ![
      ! set vertical discharge profile of cdr releases

      ! This is a time-independent equation for now, but
      ! it is a placeholder for a more sophisticated time-evolving
      ! profile in future as per the requirement.

      implicit none

      ! local
      integer :: i,j,icdr

      ! Loop through all cdr releases end set the profile:
      do icdr=1,ncdr
        cdr_prf(icdr,:)=0     ! Set all values to zero
        cdr_prf(icdr,1)= 0.5  ! Dispersion profile bottom cell
        cdr_prf(icdr,2)= 0.5  ! Dispersion profile 2nd from bottom cell
        do k=1,nz
          cdr_prf(icdr,k)= exp( z_r(i,j,k)
        enddo
     
      enddo

      end subroutine set_cdr_vert_prf  !]
! ----------------------------------------------------------------------
      subroutine set_ana_cdr_frc  ![
      ! Analytical cdr forcing volume and tracer data
      ! Put here to avoid circular reference if in analytical.F

#include "ana_cdr_frc.h"

      end subroutine set_ana_cdr_frc  !]
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
      dist=sin((dLat2)**2 + cos(lat)*cos(lat0) * sin(dLon2)**2
      dist=2*atan2(sqrt(dist),sqrt(1-dist));

      dist = radius*dist;    !Haversine distance
      deallocate(dLon2)
      deallocate(dLat2)

      end subroutine ll2dist  !]
! ----------------------------------------------------------------------

      end module cdr_frc
