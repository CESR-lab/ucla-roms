      module surf_flux
      ! Declaration of surface flux variables
      ! Output of surface fluxes
      ! Sets fields for flux correction 

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.h"
! CPP flag XXXX covering entire module so that the module is empty if
! BULK_FLUX is unused, and thus module will not take up unnecessary size in
! the executable.
! Currently only set-up for BULK_FLUX but eventually should be used for all
! fluxes once forces.h is fully ported into flux modules.
#if defined BULK_FLUX /* for whole module! */  

      implicit none

      ! Make all variable private to this module unless public specified
      private

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_array to work. NT = number tracer from param.h
#include "param.h"

      ! ****************************************************************
      ! User inputs

      logical :: write_file = .true. ! Whether to write a flux file or not
      integer :: rec_rate = 1        ! Output record rate,
      integer :: recs_per_file = 100 ! total recs per file

      logical :: wrt_sustr = .True. ! Output surface u-stress
      logical :: wrt_svstr = .True. ! Output surface v-stress
      logical :: wrt_shflx = .True. ! Output surface heat flux
      logical :: wrt_ssflx = .True. ! Output surface salinity flux
      !logical :: wrt_scflx = .true. ! for instance all surface chemical fluxes

      ! End of user inputs
      ! *************************************************************

      ! Surface momemtum flux [m^2/s^2] units as per Eq.Sys. m^2/s^2 not N/m^2
      ! possibly make sustr_r,svstr_r private to bulk_frc, and add here ustar instead
      real, public :: sustr(GLOBAL_2D_ARRAY)    ! stress u-point: used in Eq.System
      real, public :: sustr_r(GLOBAL_2D_ARRAY)  ! rho-point, only used in lmd_kpp to get ustar
      real, public :: svstr(GLOBAL_2D_ARRAY)    ! v-point: used in Eq.System
      real, public :: svstr_r(GLOBAL_2D_ARRAY)  ! rho-point, only used in lmd_kpp to get ustar
      real, public :: stflx(GLOBAL_2D_ARRAY,NT) ! Surface fluxes of tracer type variables (rho-points)
      character(len=5) :: sustr_name        = 'sustr'
      character(len=5) :: svstr_name        = 'svstr'
      character(len=5) :: shflx_name        = 'shflx'
      character(len=5) :: ssflx_name        = 'ssflx'


      ! Currently haven't incorporated dQdSST

      ! SST and SSS are needed to do flux-correction which is basically
      ! a restoring terms towards a climatological SST/SSS field.
      ! ----------------------------------------------------------------
      ! sea-surface temperature (SST) data
#if defined QCORRECTION && !defined ANA_SST
      real, public :: sst(GLOBAL_2D_ARRAY)
      real         :: sst_data(GLOBAL_2D_ARRAY,2)
      integer :: sst_grd_type  = 0
      integer :: sst_file_indx = 0
      integer :: sst_irec      = 0
      integer :: sst_it1 = 1, sst_it2 = 2
      real, dimension(2) :: sst_times = [-99,-99]
      character(len=5)   :: sst_name = 'sst'
      character(len=8)   :: sst_time_name = 'sst_time'
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      real, public :: sss(GLOBAL_2D_ARRAY)
      real         :: sss_data(GLOBAL_2D_ARRAY,2)
      integer :: sss_grd_type  = 0
      integer :: sss_file_indx = 0
      integer :: sss_irec      = 0
      integer :: sss_it1 = 1, sss_it2 = 2
      real, dimension(2) :: sss_times = [-99,-99]
      character(len=5)   :: sss_name = 'SSS'
      character(len=8)   :: sss_time_name = 'sss_time'

      real, public :: dSSSdt ! Input units (cm/day)
#endif


      ! Netcdf outputting
      ! -----------------
      integer :: output_rec = 0 ! record number of output. 0 indicates we need new file
      integer :: total_output_rec = 0 ! Total flux output recordings so far
      ! Suffix to append to simulation output rootname
      character (len=32) :: file_name = '_surf_flux.nc'


      public set_surf_field_corr
      public wrt_surface_flux_his

! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin
      ! Debug test read in coarse grid.
      ! -------------------------------
      real, public :: h_interp(GLOBAL_2D_ARRAY)
      real, public :: lat_interp(GLOBAL_2D_ARRAY)
#endif

! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin_frc
      ! Test read in coarse forcing.
      ! -------------------------------
      real, public :: uwnd_int(GLOBAL_2D_ARRAY)
      real         :: uwnd_int_data(GLOBAL_2D_ARRAY,2)
      integer :: uwnd_int_grd_type  = 0
      integer :: uwnd_int_file_indx = 0
      integer :: uwnd_int_irec      = 0
      integer :: uwnd_int_it1 = 1, uwnd_int_it2 = 2
      real, dimension(2) :: uwnd_int_times = [-99,-99]
      ! Note: need to rename uwnd in coarse input file else variable name
      ! will be confused with the refined grid input file containing 'uwnd'
      character(len=8)   :: uwnd_int_name = 'uwnd-int'
      character(len=8)   :: uwnd_int_time_name = 'wnd_time'

      ! Used to copy uwnd from bulk force as can't refer directly due to
      ! circular reference
      real, public :: uwnd_bf(GLOBAL_2D_ARRAY)
#endif



      contains
! ----------------------------------------------------------------------

      subroutine set_surf_field_corr(istr,iend,jstr,jend)
      ! Set surface fields that will be restored towards

      use read_write, only: set_frc_var_tile
     &  , ncdf_read_coarser_grid ! coarser only for ddevin debug

      implicit none

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend

      integer :: ierr = 0 ! debug coarse temporary delete later.

! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin
      integer :: tmp1=0 ! debug delete later
      integer :: tmp2=0
      integer :: tmp3=0
      integer :: tmp4=0
#endif


#if defined QCORRECTION && !defined ANA_SST
      ! Sea-surface temperature (SST) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sst_name,       sst_time_name,  ! Netcdf names for var and time
     &          sst_data,       sst,            ! Input data and var arrays
     &          sst_times,      sst_grd_type,   ! Input times & Grid type
     &          sst_file_indx,  sst_irec,       ! File indx & current input record
     &          sst_it1,        sst_it2 )       ! Indices to input data and time
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sss_name,       sss_time_name,  ! Netcdf names for var and time
     &          sss_data,       sss,            ! Input data and var arrays
     &          sss_times,      sss_grd_type,   ! Input times & Grid type
     &          sss_file_indx,  sss_irec,       ! File indx & current input record
     &          sss_it1,        sss_it2 )       ! Indices to input data and time
#endif

! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin
      ! Debug test read in coarse grid.
      ! -------------------------------
!      ierr = ncdf_read_coarser_grid(1,1,0,0,lat_interp,1) ! ncid,varid,record,horiz_type,A,nmax)
      ierr = ncdf_read_coarser_grid(tmp1,tmp2,tmp3,tmp4,lat_interp,1)
#endif

! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin_frc
      ! Debug test read in coarse forcing
      call set_frc_var_tile(
     &          istr,iend,           jstr,jend,           ! Tile bounds indices
     &          uwnd_int_name,       uwnd_int_time_name,  ! Netcdf names for var and time
     &          uwnd_int_data,       uwnd_int,            ! Input data and var arrays
     &          uwnd_int_times,      uwnd_int_grd_type,   ! Input times & Grid type
     &          uwnd_int_file_indx,  uwnd_int_irec,       ! File indx & current input record
     &          uwnd_int_it1,        uwnd_int_it2, 1 )       ! Indices to input data and time
#endif

      end subroutine set_surf_field_corr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine wrt_surface_flux_his
      ! Write surface fluxes to output netcdf file

      use netcdf ! for nf90_routines
      use read_write, only: ncdf_create_file, nc_define_var, nc_write_var
     &   , nc_write_time, output_root_name, rp_var, up_var, vp_var

      implicit none

      ! Local
      ! -----
#include "scalars.h"
      integer lenstr, my_nf_def_dim
      integer lfnm, lvar, prev_fill_mode,
     &        r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)
#ifdef SOLVE3D
     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif

      integer var_id_tmp
      integer lstr ! length of string
      character(len=64) :: fname ! Taken from read_inp.F
      integer :: ierr = 0 ! This is local so need to set to no-error
      integer :: ncid = -1 ! File ID number - will be set below


      ! Only plot if write_file is true from roms.in or top of module
      ! and timestep is divisble by record rate
      ! =======================================


      if (write_file==.true. .and. mod(iic-1,rec_rate) == 0) then ! 1 = ntstart

        ! Create file
        ! ===========

        ! Add root output name to file name (only done once)
        if(total_output_rec == 0) then
          ! Notice below space between '/ /', this is needed to avoid
          ! cpp preprocessor seeing // as a c++ comment and deleting it.
          file_name = trim(output_root_name) / / trim(file_name)
        endif

        ! If file doesn't yet exist or need new file!
        if(output_rec == 0) then

          !! lose the prev_fill_mode?
          ! r2dgrd, ...., auxil are id's to dimensions, needed to
          ! create vars in the file
          call ncdf_create_file( file_name, ncid, prev_fill_mode,
     &                   recs_per_file, total_output_rec, auxil,
#ifdef SOLVE3D
     &                   r3dgrd, u3dgrd, v3dgrd, w3dgrd, ! 3D grid dimensions
#endif
     &                   r2dgrd, u2dgrd, v2dgrd )        ! 2D grid dimensions


          ! Define remaining variables
          ! --------------------------

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


! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin
          ! debug only: output of hraw interpolation from coarse grid
          call nc_define_var(ncid, 'lat_rho',
     &      'lat_rho interpolated from coarse grid','m', r2dgrd, ierr )
#endif

#ifdef debug_ddevin_frc
          ! debug only:
            call nc_define_var(ncid, 'uwnd-int', ! Note: rho-point
     &            'uwnd interpolated from coarse grid','m', r2dgrd, ierr )
            call nc_define_var(ncid, 'uwnd', ! Note: rho-point
     &            'uwnd bulk_frc','m', r2dgrd, ierr )
#endif

          ! Leave definition mode.
          ! ----------------------
          ierr=nf90_enddef(ncid) ! End definition to write changes to disk
          if (mynode == 0) then
            write(*,'(6x,2A)') 'surf_flux :: created new netCDF file '
     &                          ,file_name
          endif

        endif ! (output_rec == 0; needed to create a file)


        ! Info to terminal output
        ! =======================

        ! Only done once - commented out if read_inp_net_flux is used
        if(total_output_rec == 0) then

          ! Write to terminal output in simulation pre-amble text which
          ! result variables are being stored
          ! ---------------------------------
          mpi_master_only write(*,'(/1x,A,L1,2x,A,I5,2x,A,I4,2x,2A)')
     &     'surf_flux: save results net_flux = ', write_file,
     &     'rec_rate =', rec_rate, 'recs/file =', recs_per_file,
     &     'file = ', file_name

          mpi_master_only write(*,'(/1x,A,4(/8x,A,T16,L1),/)') ! Change to number of variables 4(
     &               'fields to be saved in 2D Flux history: (T/F)'
     &     ,sustr_name,  wrt_sustr
     &     ,svstr_name,  wrt_svstr
     &     ,shflx_name,  wrt_shflx
     &     ,ssflx_name,  wrt_ssflx

        endif ! (total_output_rec == 0)


        ! WRITE VARIABLES TO FILE
        ! =======================

        ! Open netcdf file
        ierr=nf90_open(file_name, nf90_write, ncid)

        ! Set fill value - nf90_nofill for optimized writing
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        ! Advance output record number
        output_rec = output_rec + 1 ! Current file record
        total_output_rec = total_output_rec + 1 ! total number of recs written

        ! write ocean time
        ! ----------------
        call nc_write_time(ncid, output_rec, total_output_rec)

        ! write variables
        ! ---------------
        if (wrt_sustr) then
          call nc_write_var(ncid, sustr, 1, sustr_name, up_var, output_rec, ierr)
        endif
        if (wrt_svstr) then
          call nc_write_var(ncid, svstr, 1, svstr_name, vp_var, output_rec, ierr)
        endif
        if (wrt_shflx) then
          call nc_write_var(ncid, stflx(GLOBAL_2D_ARRAY,itemp), 1,
     &                      shflx_name, rp_var, output_rec, ierr)
        endif
        if (wrt_ssflx) then
          call nc_write_var(ncid, stflx(GLOBAL_2D_ARRAY,isalt), 1,
     &                      ssflx_name, rp_var, output_rec, ierr)
        endif

!          ! Output of sustr_rho and svstr_rho
!          call nc_write_var(ncid, sustr_r, 1, 'sustr_r', rp_var, output_rec, ierr)
!          call nc_write_var(ncid, svstr_r, 1, 'svstr_r', rp_var, output_rec, ierr)

! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin
        ! debug only: output of hraw interpolation from coarse grid
        call nc_write_var(ncid,lat_interp,1,'lat_rho',rp_var,output_rec,ierr)
#endif

#ifdef debug_ddevin_frc
        ! debug only:
          call nc_write_var(ncid, uwnd_int, 1, 'uwnd-int', rp_var, output_rec, ierr)
          call nc_write_var(ncid, uwnd_bf, 1, 'uwnd', rp_var, output_rec, ierr)
#endif

        ! Close netcdf file
        ! -----------------
        ierr=nf90_close (ncid)

        ! Confirm work completed
        ! ----------------------
        if (mynode == 0) then
          write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'flux_his :: wrote history, tdays =', tdays,
     &     'step =', iic-1, 'rec =', output_rec, '/', total_output_rec
     &      MYID
        endif

        ! Reset output_rec to create new file
        if (output_rec == recs_per_file) then
          output_rec = 0 
        endif


      endif ! (write_file==.true.)

! Leaving debug functionality in until fully testsed
#ifdef debug_ddevin
      ! Debug test read in coarse grid.
      ! -------------------------------
      stop 'stop wrt_surf_flux_his debug'
#endif

      end subroutine wrt_surface_flux_his


!-----------------------------------------------------------------------

#endif /* BULK_FLUX for whole module */

      end module surf_flux
