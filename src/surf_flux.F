      module surf_flux
      ! Declaration of surface flux variables
      ! Output of surface fluxes
      ! Sets fields for flux correction 

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.h"

      use read_write ! rp_var, up_var, vp_var

      implicit none

      ! Make all variable private to this module unless public specified
      private

      ! Includes:
      ! param.h needed for GLOBAL_2D_array to work. NT = number tracer from param.h
#include "param.h"
      ! ncvars.h for ncavg
#include "ncvars.h"

      ! ****************************************************************
      ! User inputs

      ! Interpolate forcing from coarser input grid (=1) or not (=0)
      integer :: interp_frc = 0 ! For SST and SSS correction

      logical :: write_file    = .false. ! Whether to write a flux file or not
      integer :: rec_rate      = 2       ! Output record rate,
      integer :: recs_per_file = 3     ! total recs per file

      logical :: wrt_sustr = .True. ! Output surface u-stress
      logical :: wrt_svstr = .True. ! Output surface v-stress
      logical :: wrt_shflx = .True. ! Output surface heat flux
      logical :: wrt_ssflx = .True. ! Output surface salinity flux
      !logical :: wrt_scflx = .true.! for instance all surface chemical fluxes

      ! Writing averages to the main _avg.F file. Output averaging controlled
      ! by roms.in, but T/F for printing variables is controlled here.
      logical :: wrt_sustr_avg = .True. ! Output average surface u-stress
      logical :: wrt_svstr_avg = .True. ! Output average surface v-stress
      logical :: wrt_shflx_avg = .True. ! Output average surface heat flux
      logical :: wrt_ssflx_avg = .True. ! Output average surface salinity flux

      ! End of user inputs
      ! *************************************************************

      ! Surface momemtum flux [m^2/s^2] units as per Eq.Sys. m^2/s^2 not N/m^2
      ! possibly make sustr_r,svstr_r private to bulk_frc, and add here ustar instead
      real, public :: sustr(GLOBAL_2D_ARRAY)    ! stress u-point: used in Eq.System
      real, public :: sustr_r(GLOBAL_2D_ARRAY)  ! rho-point, only used in lmd_kpp to get ustar
      real, public :: svstr(GLOBAL_2D_ARRAY)    ! v-point: used in Eq.System
      real, public :: svstr_r(GLOBAL_2D_ARRAY)  ! rho-point, only used in lmd_kpp to get ustar
      real, public :: stflx(GLOBAL_2D_ARRAY,NT) ! Surface fluxes of tracer type variables (rho-points)
      real, public :: srflx(GLOBAL_2D_ARRAY)    ! Short-wave radiation surface flux
      character(len=5) :: sustr_name = 'sustr'
      character(len=5) :: svstr_name = 'svstr'
      character(len=5) :: shflx_name = 'shflx'
      character(len=5) :: ssflx_name = 'ssflx'

      ! Output surface flux averages
#if defined AVERAGES && defined SURF_FLUX_OUTPUT_AVG
      real :: sustr_avg(GLOBAL_2D_ARRAY)    ! surf stress u-point average
      real :: svstr_avg(GLOBAL_2D_ARRAY)    ! surf stress v-point average
      real :: shflx_avg(GLOBAL_2D_ARRAY)    ! surf heat flux average
      real :: ssflx_avg(GLOBAL_2D_ARRAY)    ! surf salinity flux average
      character(len=9) :: sustr_avg_name = 'sustr_avg'
      character(len=9) :: svstr_avg_name = 'svstr_avg'
      character(len=9) :: shflx_avg_name = 'shflx_avg'
      character(len=9) :: ssflx_avg_name = 'ssflx_avg'
#endif


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
      real, public :: dSSTdt
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
      character (len=60) :: file_name = '_surf_flux.nc'

      ! Public functions
      public set_surf_field_corr
      public wrt_surface_flux_his
#if defined AVERAGES && defined SURF_FLUX_OUTPUT_AVG
      public set_avg_surf_flux
      public def_avg_surf_flux
      public wrt_avg_surf_flux
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


#if defined QCORRECTION && !defined ANA_SST
      ! Sea-surface temperature (SST) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sst_name,       sst_time_name,  ! Netcdf names for var and time
     &          sst_data,       sst,            ! Input data and var arrays
     &          sst_times,      sst_grd_type,   ! Input times & Grid type
     &          sst_file_indx,  sst_irec,       ! File indx & current input record
     &          sst_it1,        sst_it2,        ! Time index placeholders
     &          interp_frc )                    ! Online spatial interpolation flag
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sss_name,       sss_time_name,  ! Netcdf names for var and time
     &          sss_data,       sss,            ! Input data and var arrays
     &          sss_times,      sss_grd_type,   ! Input times & Grid type
     &          sss_file_indx,  sss_irec,       ! File indx & current input record
     &          sss_it1,        sss_it2,        ! Time index placeholders
     &          interp_frc )                    ! Online spatial interpolation flag
#endif


      end subroutine set_surf_field_corr
!-----------------------------------------------------------------------

#if defined AVERAGES && defined SURF_FLUX_OUTPUT_AVG
      subroutine set_avg_surf_flux(istrR,iendR,jstrR,jendR,cff,cff1)
      ! Set surface flux averages

      implicit none

      ! inputs
      integer istrR,iendR,jstrR,jendR
      real cff,cff1
      ! local
      integer i,j,k,itrc

      do j=jstrR,jendR
        do i=istrR,iendR
          if (wrt_sustr_avg) then
              sustr_avg(i,j)=cff1*sustr_avg(i,j)
     &                      + cff*sustr(i,j)
          endif
          if (wrt_svstr_avg) then
              svstr_avg(i,j)=cff1*svstr_avg(i,j)
     &                      + cff*svstr(i,j)
          endif
          if (wrt_shflx_avg) then
              shflx_avg(i,j)=cff1*shflx_avg(i,j)
     &                      + cff*stflx(i,j,itemp)
          endif
          if (wrt_ssflx_avg) then
              ssflx_avg(i,j)=cff1*ssflx_avg(i,j)
     &                      + cff*stflx(i,j,isalt)
          endif
        enddo
      enddo

      end subroutine set_avg_surf_flux
#endif /* AVERAGES && SURF_FLUX_OUTPUT_AVG*/
!-----------------------------------------------------------------------

#if defined AVERAGES && defined SURF_FLUX_OUTPUT_AVG
      subroutine def_avg_surf_flux(r2dgrd, u2dgrd, v2dgrd)
      ! Define surface flux averaged variables in average output file.

      implicit none

      ! Inputs
      integer r2dgrd(3), u2dgrd(3), v2dgrd(3)
      ! Local
      integer ierr

      if (wrt_sustr) then
        call nc_define_var(ncavg, sustr_avg_name,
     &   'Averaged wind stress in x-direction','m^2/s^2', u2dgrd, ierr )
      endif
      if (wrt_svstr) then
        call nc_define_var(ncavg, svstr_avg_name,
     &   'Averaged wind stress in y-direction', 'm^2/s^2', v2dgrd, ierr )
      endif
      if (wrt_shflx) then
        call nc_define_var(ncavg, shflx_avg_name,
     &   'Averaged Surface heat flux', 'degC m/s', r2dgrd, ierr )
      endif
      if (wrt_ssflx) then
        call nc_define_var(ncavg, ssflx_avg_name,
     &   'Averaged Surface salinity flux','PSU m/s', r2dgrd, ierr )
      endif


      end subroutine def_avg_surf_flux
#endif /* AVERAGES && SURF_FLUX_OUTPUT_AVG*/
!-----------------------------------------------------------------------

#if defined AVERAGES && defined SURF_FLUX_OUTPUT_AVG
      subroutine wrt_avg_surf_flux(record, ierr)
      ! Write averaged surface flux variables to averaged output file

      implicit none

      ! inputs
      integer, intent(in) :: record ! The timestep to be recorded
      integer :: ierr

      if (wrt_sustr_avg) then
        call nc_write_var(ncavg, sustr_avg, 1, sustr_avg_name, up_var,
     &                    record, ierr)
      endif
      if (wrt_svstr_avg) then
        call nc_write_var(ncavg, svstr_avg, 1, svstr_avg_name, vp_var,
     &                    record, ierr)
      endif
      if (wrt_shflx_avg) then
        call nc_write_var(ncavg, shflx_avg, 1,
     &                    shflx_avg_name, rp_var, record, ierr)
      endif
      if (wrt_ssflx_avg) then
        call nc_write_var(ncavg, ssflx_avg, 1,
     &                    ssflx_avg_name, rp_var, record, ierr)
      endif


      end subroutine wrt_avg_surf_flux
#endif /* AVERAGES && SURF_FLUX_OUTPUT_AVG*/
!-----------------------------------------------------------------------

      subroutine wrt_surface_flux_his
      ! Write surface fluxes to output netcdf file

      use netcdf ! for nf90_routines
      use read_write !, only: ncdf_create_file, nc_define_var, nc_write_var
                     !&   , nc_write_time, output_root_name, rp_var, up_var, vp_var

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
          call ncdf_create_file(file_name, ncid, prev_fill_mode
     &                  ,recs_per_file, total_output_rec, auxil
#ifdef SOLVE3D
     &                  ,r3dgrd, u3dgrd, v3dgrd, w3dgrd ! 3D grid dimensions
#endif
     &                  ,r2dgrd, u2dgrd, v2dgrd )        ! 2D grid dimensions


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

      end subroutine wrt_surface_flux_his
!-----------------------------------------------------------------------


      end module surf_flux
