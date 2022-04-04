      module surf_flux
      ! Declaration of surface flux variables
      ! Output of surface fluxes
      ! Sets fields for flux correction 

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.opt"

      use param
      use read_write
      use roms_read_write
      use ncvars
      use scalars

      implicit none

      private

#include "surf_flux.opt"

      ! Surface momemtum flux [m^2/s^2] units as per Eq.Sys. m^2/s^2 not N/m^2
      ! possibly make sustr_r,svstr_r private to bulk_frc, and add here ustar instead
      real,public,allocatable,dimension(:,:)   :: sustr    ! stress u-point: used in Eq.System
      real,public,allocatable,dimension(:,:)   :: sustr_r  ! rho-point, only used in lmd_kpp to get ustar
      real,public,allocatable,dimension(:,:)   :: svstr    ! v-point: used in Eq.System
      real,public,allocatable,dimension(:,:)   :: svstr_r  ! rho-point, only used in lmd_kpp to get ustar
      real,public,allocatable,dimension(:,:,:) :: stflx    ! Surface fluxes of tracer type variables (rho-points)
      real,public,allocatable,dimension(:,:)   :: srflx    ! Short-wave radiation surface flux
      character(len=5) :: sustr_name = 'sustr'
      character(len=5) :: svstr_name = 'svstr'
      character(len=5) :: shflx_name = 'shflx'
      character(len=5) :: ssflx_name = 'ssflx'

      ! Output surface flux averages
      real,allocatable,dimension(:,:) :: sustr_avg         ! surf stress u-point average
      real,allocatable,dimension(:,:) :: svstr_avg         ! surf stress v-point average
      real,allocatable,dimension(:,:) :: shflx_avg         ! surf heat flux average
      real,allocatable,dimension(:,:) :: ssflx_avg         ! surf salinity flux average
      character(len=9) :: sustr_avg_name = 'sustr_avg'
      character(len=9) :: svstr_avg_name = 'svstr_avg'
      character(len=9) :: shflx_avg_name = 'shflx_avg'
      character(len=9) :: ssflx_avg_name = 'ssflx_avg'


      ! Currently haven't incorporated dQdSST

      ! SST and SSS are needed to do flux-correction which is basically
      ! a restoring terms towards a climatological SST/SSS field.
      ! ----------------------------------------------------------------
#if defined QCORRECTION && !defined ANA_SST
      ! sea-surface temperature (SST) data
      real,public,allocatable,dimension(:,:) :: sst
      real,public :: dSSTdt
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      real,public,allocatable,dimension(:,:) :: sss
      real,public :: dSSSdt                                ! input units (cm/day)
#endif


      ! Netcdf outputting
      ! -----------------
      integer            :: ncid = -1
      integer            :: output_rec = 0               ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0         ! total flux output recordings so far
      character (len=60) :: file_name = '_surf_flux.nc'  ! suffix to append to simulation output rootname

      ! Public functions
      public init_arrays_surf_flx
      public set_surf_field_corr
      public wrt_surface_flux_his
#if defined AVERAGES
      public set_avg_surf_flux
      public def_avg_surf_flux
      public wrt_avg_surf_flux
#endif

      contains

! ----------------------------------------------------------------------
      subroutine init_arrays_surf_flx  ![
      use scalars
      implicit none

      allocate( sustr  (GLOBAL_2D_ARRAY)    ); sustr=init
      allocate( sustr_r(GLOBAL_2D_ARRAY)    ); sustr_r=init
      allocate( svstr  (GLOBAL_2D_ARRAY)    ); svstr=init
      allocate( svstr_r(GLOBAL_2D_ARRAY)    ); svstr_r=init
      allocate( stflx  (GLOBAL_2D_ARRAY,NT) ); stflx=init
      allocate( srflx  (GLOBAL_2D_ARRAY)    ); srflx=init

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      allocate( sss(GLOBAL_2D_ARRAY)        ); sss=init
      allocate(  nc_sss%data( GLOBAL_2D_ARRAY,2) )
#endif

      if (wrt_surf_flx_avg) then
        allocate( sustr_avg(GLOBAL_2D_ARRAY) )
        allocate( svstr_avg(GLOBAL_2D_ARRAY) )   ! surf stress v-point average
        allocate( shflx_avg(GLOBAL_2D_ARRAY) )   ! surf heat flux average
        allocate( ssflx_avg(GLOBAL_2D_ARRAY) )
      endif

      end subroutine init_arrays_surf_flx  !]

! ----------------------------------------------------------------------
      subroutine set_surf_field_corr(istr,iend,jstr,jend)  ![
      ! Set surface fields that will be restored towards

      implicit none

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend


#if defined QCORRECTION && !defined ANA_SST
      ! Sea-surface temperature (SST) data
      call set_frc_var_tile( istr, iend, jstr, jend, nc_sst, sst, interp_frc )
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_var_tile( istr, iend, jstr, jend, nc_sss, sss, interp_frc )
#endif


      end subroutine set_surf_field_corr  !]

!-----------------------------------------------------------------------
      subroutine set_avg_surf_flux  ![
      ! set surface flux averages

      implicit none

      ! local
      real :: cff, cff1

      cff=1./dble(rec_rate)    !<-- scaling
      if (mod(iic-ntstart-1,rec_rate) == 0) then ! iic-ntstart gives timestep. -1 give next averaging cycle.
        if (mynode==0) print *, 'Refresh averaging cycle - surf flux'
        cff1=0.                !--> initialize
      else
        cff1=1.                !--> keep adding
      endif

      if (wrt_sustr_avg) sustr_avg = cff1*sustr_avg + cff*sustr
      if (wrt_svstr_avg) svstr_avg = cff1*svstr_avg + cff*svstr
      if (wrt_shflx_avg) shflx_avg = cff1*shflx_avg + cff*stflx(:,:,itemp)
# ifdef SALINITY
      if (wrt_ssflx_avg) ssflx_avg = cff1*ssflx_avg + cff*stflx(:,:,isalt)
# endif

      if (mynode==0) print *, 'sustr_avg(1,1)=', sustr_avg(1,1), 'sustr(1,1)=', sustr(1,1)

      end subroutine set_avg_surf_flux  !]

!-----------------------------------------------------------------------
      subroutine def_avg_surf_flux(r2dgrd, u2dgrd, v2dgrd)  ![
      ! Define surface flux averaged variables in average output file.

      implicit none

      ! Inputs
      integer r2dgrd(3), u2dgrd(3), v2dgrd(3)
      ! Local
      integer ierr

      if (wrt_sustr) then
        call nc_define_var(ncid, sustr_avg_name,
     &   'Averaged wind stress in x-direction','m^2/s^2', u2dgrd, ierr )
      endif
      if (wrt_svstr) then
        call nc_define_var(ncid, svstr_avg_name,
     &   'Averaged wind stress in y-direction', 'm^2/s^2', v2dgrd, ierr )
      endif
      if (wrt_shflx) then
        call nc_define_var(ncid, shflx_avg_name,
     &   'Averaged Surface heat flux', 'degC m/s', r2dgrd, ierr )
      endif
      if (wrt_ssflx) then
        call nc_define_var(ncid, ssflx_avg_name,
     &   'Averaged Surface salinity flux','PSU m/s', r2dgrd, ierr )
      endif


      end subroutine def_avg_surf_flux  !]

!-----------------------------------------------------------------------
      subroutine wrt_avg_surf_flux  ![
      ! rite averaged surface flux variables to averaged output file
      implicit none

      ! local
      integer :: ierr = 0

      if (wrt_sustr_avg) then
        call nc_write_var(ncid, sustr_avg, 1, sustr_avg_name, up_var, output_rec, ierr)
      endif
      if (wrt_svstr_avg) then
        call nc_write_var(ncid, svstr_avg, 1, svstr_avg_name, vp_var, output_rec, ierr)
      endif
      if (wrt_shflx_avg) then
        call nc_write_var(ncid, shflx_avg, 1, shflx_avg_name, rp_var, output_rec, ierr)
      endif
      if (wrt_ssflx_avg) then
        call nc_write_var(ncid, ssflx_avg, 1, ssflx_avg_name, rp_var, output_rec, ierr)
      endif

      end subroutine wrt_avg_surf_flux  !]

!-----------------------------------------------------------------------
      subroutine wrt_surface_flux_his  ![
      ! Write surface fluxes to output netcdf file

      use netcdf
      use read_write

      implicit none

      ! local
      integer :: lenstr, my_nf_def_dim
      integer :: lfnm, lvar, prev_fill_mode,
     &           r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)
#ifdef SOLVE3D
     &         , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif

      integer :: ierr = 0        ! set to no-error
      character(len=64) :: fname ! taken from read_inp.F

      if (wrt_surf_flx_avg) call set_avg_surf_flux                              ! averaging needed every step

      if (wrt_surf_flx_file==.true. .and. mod(iic-ntstart,rec_rate) == 0) then  ! when timestep divisble by record rate

        if(total_output_rec == 0) then                            ! Add root output name to file name (only done once)
          file_name = trim(output_root_name) / / trim(file_name)  ! Notice below space between '/ /', this is needed to avoid
        endif                                                     ! cpp preprocessor seeing // as a c++ comment and deleting it.


        if(output_rec == 0) then                                  ! if file doesn't yet exist or need new file!

          call ncdf_create_file(file_name, ncid, prev_fill_mode   ! lose the prev_fill_mode? r2dgrd, ...., auxil are id's to dimensions, needed to create vars in the file
     &                  ,recs_per_file, total_output_rec, auxil
#ifdef SOLVE3D
     &                  ,r3dgrd, u3dgrd, v3dgrd, w3dgrd           ! 3D grid dimensions
#endif
     &                  ,r2dgrd, u2dgrd, v2dgrd )                 ! 2D grid dimensions


          ! Define remaining variables:
          if (wrt_sustr) then
            ! Output surface flux as per Eq.Sys. units m^2/s^2 not N/m^2
            call nc_define_var(ncid, sustr_name,
     &            'wind stress in x-direction','m^2/s^2', u2dgrd, ierr )
          endif
          if (wrt_svstr) then
            ! Output surface flux as per Eq.Sys. units m^2/s^2 not N/m^2
            call nc_define_var(ncid, svstr_name,
     &           'wind stress in y-direction', 'm^2/s^2', v2dgrd, ierr )
          endif
          if (wrt_shflx) then
            call nc_define_var(ncid, shflx_name,
     &                   'Surface heat flux', 'degC m/s', r2dgrd, ierr )
          endif
          if (wrt_ssflx) then
            call nc_define_var(ncid, ssflx_name,
     &                 'Surface salinity flux','PSU m/s', r2dgrd, ierr )
          endif
!          ! Output of sustr_rho and svstr_rho
!            call nc_define_var(ncid, 'sustr_r',
!     &            'rho wind stress in x-direction','m^2/s^2', r2dgrd, ierr )
!            call nc_define_var(ncid, 'svstr_r',
!     &            'rho wind stress in y-direction','m^2/s^2', r2dgrd, ierr )

          if (wrt_surf_flx_avg) call def_avg_surf_flux(r2dgrd, u2dgrd, v2dgrd)


          ierr=nf90_enddef(ncid) ! end definition to write changes to disk
          if (mynode == 0) then
            write(*,'(6x,2A)') 'surf_flux :: created new netCDF file '
     &                          ,file_name
          endif

        endif ! (output_rec == 0; needed to create a file)


        ! Info to terminal output
        ! =======================

        if(total_output_rec == 0) then                                 ! Only done once - commented out if read_inp_net_flux is used

          mpi_master_only write(*,'(/1x,A,L1,2x,A,I5,2x,A,I4,2x,2A)')  ! write to terminal output in simulation pre-amble text which
     &     'surf_flux: save results net_flux = ', wrt_surf_flx_file,          ! result variables are being stored
     &     'rec_rate =', rec_rate, 'recs/file =', recs_per_file,
     &     'file = ', file_name

          mpi_master_only write(*,'(/1x,A,4(/8x,A,T16,L1),/)')         ! Change to number of variables 4(
     &               '2D history fields to be saved: (T/F)'
     &     ,sustr_name,  wrt_sustr
     &     ,svstr_name,  wrt_svstr
     &     ,shflx_name,  wrt_shflx
     &     ,ssflx_name,  wrt_ssflx

          mpi_master_only write(*,'(/1x,A,4(/8x,A,T16,L1),/)')         ! Change to number of variables 4(
     &               '2D AVERAGED fields to be saved: (T/F)'
     &     ,sustr_name,  wrt_sustr_avg
     &     ,svstr_name,  wrt_svstr_avg
     &     ,shflx_name,  wrt_shflx_avg
     &     ,ssflx_name,  wrt_ssflx_avg

        endif ! (total_output_rec == 0)


        ! WRITE VARIABLES TO FILE
        ! =======================

        ierr=nf90_open(file_name, nf90_write, ncid)

        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)  ! set fill value - nf90_nofill for optimized writing

                                                ! advance output record number
        output_rec = output_rec + 1             ! current file record
        total_output_rec = total_output_rec + 1 ! total number of recs written


        call nc_write_time(ncid, output_rec, total_output_rec)  ! write ocean time


        if (wrt_sustr) then  ! write variables
          call nc_write_var(ncid, sustr, 1, sustr_name, up_var, output_rec, ierr)
        endif
        if (wrt_svstr) then
          call nc_write_var(ncid, svstr, 1, svstr_name, vp_var, output_rec, ierr)
        endif
        if (wrt_shflx) then
          call nc_write_var(ncid, stflx(GLOBAL_2D_ARRAY,itemp), 1,
     &                      shflx_name, rp_var, output_rec, ierr)
        endif
#ifdef SALINITY
        if (wrt_ssflx) then
          call nc_write_var(ncid, stflx(GLOBAL_2D_ARRAY,isalt), 1,
     &                      ssflx_name, rp_var, output_rec, ierr)
        endif
#endif
!          ! Output of sustr_rho and svstr_rho
!          call nc_write_var(ncid, sustr_r, 1, 'sustr_r', rp_var, output_rec, ierr)
!          call nc_write_var(ncid, svstr_r, 1, 'svstr_r', rp_var, output_rec, ierr)

        if (wrt_surf_flx_avg) call wrt_avg_surf_flux

        ierr=nf90_close (ncid)  ! close netcdf file

        if (mynode == 0) then
          write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &     'surf_flux :: wrote file, tdays =', tdays,
     &     'step =', iic, 'rec =', output_rec, '/', total_output_rec
     &      MYID
        endif

        if (output_rec == recs_per_file) then  ! reset output_rec to create new file
          output_rec = 0 
        endif

      endif ! (wrt_surf_flx_file==.true.)

      end subroutine wrt_surface_flux_his  !]

!-----------------------------------------------------------------------

      end module surf_flux
