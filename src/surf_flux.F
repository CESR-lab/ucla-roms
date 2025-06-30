      module surf_flux
      ! Declaration of surface flux variables
      ! Output of surface fluxes
      ! Sets fields for flux correction

      ! Devin Dollery & Jeroen Molemaker (2020-2022)

#include "cppdefs.opt"

      use param
      use dimensions
      use roms_read_write
      use nc_read_write
      use netcdf
      use scalars

      implicit none

      private

#include "surf_flux.opt"

#if defined SALINITY
      logical,parameter :: salinity=.true.
#else
      logical,parameter :: salinity=.false.
#endif

      logical :: init_done=.false.        ! flag to init surf_flux
      integer :: navg_sflx = 0            ! number of samples in average

      ! Surface momemtum flux [m^2/s^2] units as per Eq.Sys. m^2/s^2 not N/m^2
      ! possibly make sustr_r,svstr_r private to bulk_frc, and add here ustar instead
      real,public,allocatable,dimension(:,:)   :: sustr    ! stress u-point: used in Eq.System
      real,public,allocatable,dimension(:,:)   :: sustr_r  ! rho-point, only used in lmd_kpp to get ustar
      real,public,allocatable,dimension(:,:)   :: svstr    ! v-point: used in Eq.System
      real,public,allocatable,dimension(:,:)   :: svstr_r  ! rho-point, only used in lmd_kpp to get ustar
      real,public,allocatable,dimension(:,:,:) :: stflx    ! Surface fluxes of tracer type variables (rho-points)
      real,public,allocatable,dimension(:,:)   :: srflx    ! Short-wave radiation surface flux

      real,public,allocatable,dimension(:,:) ::  uwnd             ! time interpolated u-wind
      real,public,allocatable,dimension(:,:) ::  vwnd             ! time interpolated v-wind

      character(len=5) :: sustr_name = 'sustr'
      character(len=5) :: svstr_name = 'svstr'
      character(len=5) :: shflx_name = 'shflx'
      character(len=5) :: ssflx_name = 'ssflx'

      ! averages of surface fluxes
      real,allocatable,dimension(:,:)  :: sustr_avg     ! surf stress u-point average
      real,allocatable,dimension(:,:)  :: svstr_avg     ! surf stress v-point average
      real,allocatable,dimension(:,:,:):: stflx_avg     ! surf tracer flux average


      ! Sea-surface temperature (SST) and salinity (SSS) data for restoring
      real,public,allocatable,dimension(:,:) :: sst
      real,public :: dSSTdt
      real,public,allocatable,dimension(:,:) :: sss
      real,public :: dSSSdt                                ! input units (cm/day)

      ! Netcdf outputting:
      real    :: output_time = 0
      integer :: record = nrpf     ! to trigger the first file creation
      integer :: ncid = -1, prev_fill_mode

      ! Public functions
      public init_arrays_surf_flx
      public set_surf_field_corr
      public wrt_sflux

      contains

! ----------------------------------------------------------------------
      subroutine init_arrays_surf_flx  ![
      use scalars
      implicit none

      ! local
      allocate(  uwnd  (GLOBAL_2D_ARRAY)    );
      allocate(  vwnd  (GLOBAL_2D_ARRAY)    );
      allocate( sustr  (GLOBAL_2D_ARRAY)    ); sustr=init
      allocate( sustr_r(GLOBAL_2D_ARRAY)    ); sustr_r=init
      allocate( svstr  (GLOBAL_2D_ARRAY)    ); svstr=init
      allocate( svstr_r(GLOBAL_2D_ARRAY)    ); svstr_r=init
      allocate( stflx  (GLOBAL_2D_ARRAY,nt) ); stflx=init
      allocate( srflx  (GLOBAL_2D_ARRAY)    ); srflx=init

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      allocate( sss(GLOBAL_2D_ARRAY)        ); sss=init
      allocate(nc_sss%vdata(GLOBAL_2D_ARRAY,2) )

#endif

      if (sflx_avg) then
        allocate(sustr_avg(1:i1,j0:j1))
        allocate(svstr_avg(i0:i1,1:j1))
        allocate(stflx_avg(i0:i1,j0:j1,nt))
      endif

      end subroutine init_arrays_surf_flx  !]
! ----------------------------------------------------------------------
      subroutine set_surf_field_corr ![
      ! Set surface fields that will be restored towards

      implicit none

#if defined QCORRECTION && !defined ANA_SST
      ! Sea-surface temperature (SST) data
      print *,'Allocate space for sst first!! '
      call set_frc_data(nc_sst,sst,'r')
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_data(nc_sss,sss,'r')
#endif

      end subroutine set_surf_field_corr  !]
! ----------------------------------------------------------------------
      subroutine calc_sflx_avg  ![
      implicit none

      ! local
      real :: coef

      navg_sflx = navg_sflx +1
      coef = 1./navg_sflx

      if (wrt_smflx) then  ! surface momentum fluxes
        sustr_avg = sustr_avg*(1-coef)+sustr( 1:i1,j0:j1)*coef
        svstr_avg = svstr_avg*(1-coef)+svstr(i0:i1, 1:j1)*coef
      endif
      if (wrt_stflx) then  ! surface tracer fluxes
        stflx_avg = stflx_avg*(1-coef)+stflx(i0:i1,j0:j1,:)*coef
      endif

      end subroutine calc_sflx_avg !]
!----------------------------------------------------------------------
      subroutine create_sflx_vars(ncid)  ![
      ! Add sflux  variables to an opened netcdf file
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer           :: ierr, varid

      ! output surface flux as per Eq.Sys. units m^2/s^2 not N/m^2
      if (wrt_smflx) then
        varid = nccreate(ncid,'sustr',(/dn_xu,dn_yr,dn_tm/),
     &                   (/xi_u,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',
     &         'wind stress in x-direction')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
        varid = nccreate(ncid,'svstr',(/dn_xr,dn_yv,dn_tm/),
     &                   (/xi_rho,eta_v,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',
     &         'wind stress in y-direction')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
      endif

      if (wrt_stflx) then
        varid = nccreate(ncid,shflx_name,(/dn_xr,dn_yr,dn_tm/),
     &                 (/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',
     &                            'Surface heat flux')
        ierr = nf90_put_att(ncid,varid,'units','degC m/s')
        if (salinity) then
          varid = nccreate(ncid,'ssflx',(/dn_xr,dn_yr,dn_tm/),
     &                 (/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',
     &                            'Surface Salinity flux')
          ierr = nf90_put_att(ncid,varid,'units','PSU m/s')
        endif
      endif

      end subroutine create_sflx_vars  !]
! ----------------------------------------------------------------------
      subroutine wrt_sflux  ![
      ! write surface flux variables to netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      implicit none

      ! local
      integer,dimension(4)   :: start
      character(len=99),save :: fname
      integer                :: ierr

      output_time = output_time + dt

      if (sflx_avg) call calc_sflx_avg

      if (output_time>=output_period) then  ! time for an output
        output_time = 0
        navg_sflx   = 0

        if (record==nrpf) then
          call create_sflx_file(fname)
          record = 0
        endif
        record = record + 1

        ierr=nf90_open(fname,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        start=1; start(3)=record
        if (sflx_avg) then
          if (wrt_smflx) then
            call ncwrite(ncid,'sustr',sustr_avg,start)
            call ncwrite(ncid,'svstr',svstr_avg,start)
          endif
          if (wrt_stflx) then
            call ncwrite(ncid,'shflx',stflx_avg(:,:,1),start)
            if (salinity) then
             call ncwrite(ncid,'ssflx',stflx_avg(:,:,2),start)
            endif
          endif
        else  ! snapshots
          if (wrt_smflx) then
            call ncwrite(ncid,'sustr',sustr(1:i1,j0:j1),start)
            call ncwrite(ncid,'svstr',svstr(i0:i1,1:j1),start)
          endif
          if (wrt_stflx) then
            call ncwrite(ncid,'shflx',stflx(i0:i1,j0:j1,1),start)
            if (salinity) then
             call ncwrite(ncid,'ssflx',stflx(i0:i1,j0:j1,2),start)
            endif
          endif
        endif
        ierr=nf90_close(ncid)

        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'surf_flux :: wrote surface flux, tdays =', tdays,
     &        'step =', iic, 'rec =', record
        endif
      endif  ! time for an output

      end subroutine wrt_sflux  !]
!----------------------------------------------------------------------
      subroutine create_sflx_file(fname)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      ! local
      integer :: ierr,varid
      character(len=10),dimension(4) :: dimnames           ! dimension names
      integer,          dimension(4) :: dimsizes

      if (sflx_avg) then
        call create_file('_flx_avg',fname)
      else
        call create_file('_flx_his',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      call create_sflx_vars(ncid)

      if (sflx_avg) then
        ierr=nf90_put_att(ncid,nf90_global,'type','surface flux average')
      else
        ierr=nf90_put_att(ncid,nf90_global,'type','surface flux history')
      endif

      ierr = nf90_close(ncid)

      end subroutine create_sflx_file !]
!-----------------------------------------------------------------------

      end module surf_flux
