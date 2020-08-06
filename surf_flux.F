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
      integer :: rec_rate = 2        ! Output record rate,
      integer :: recs_per_file = 3   ! total recs per file

      logical :: wrt_sustr = .True.  ! Output surface u-stress
      logical :: wrt_svstr = .True. ! Output surface v-stress
      logical :: wrt_shflx = .True. ! Output surface heat flux
      logical :: wrt_ssflx = .True. ! Output surface salinity flux
      !logical :: wrt_scflx = .true. ! for instance all surface chemical fluxes

      ! End of user inputs
      ! *************************************************************

      real, public :: sustr(GLOBAL_2D_ARRAY) ! Surface momemtum flux [N/m^2]
      real, public :: svstr(GLOBAL_2D_ARRAY) ! Surface momemtum flux [N/m^2]
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
      real         :: sst(GLOBAL_2D_ARRAY)
      real, public :: sst_data(GLOBAL_2D_ARRAY,2)
      integer :: sst_grd_type  = 0
      integer :: sst_file_indx = 0
      integer :: sst_irec      = 0
      integer :: sst_it1 = 1, sst_it2 = 2
      real, dimension(2) :: sst_times = [-99,-99]
      character(len=5)   :: sst_name = 'sst'
      character(len=8)   :: sst_time_name = 'sst_time'

      ! Sea-surface salinity (SSS) data
      real        :: sss(GLOBAL_2D_ARRAY)
      real,public :: sss_data(GLOBAL_2D_ARRAY,2)
      integer :: sss_grd_type  = 0
      integer :: sss_file_indx = 0
      integer :: sss_irec      = 0
      integer :: sss_it1 = 1, sss_it2 = 2
      real, dimension(2) :: sss_times = [-99,-99]
      character(len=5)   :: sss_name = 'sss'
      character(len=8)   :: sss_time_name = 'sss_time'


      ! Netcdf outputting
      ! -----------------
      integer :: output_rec = 0 ! record number of output. 0 indicates we need new file
      integer :: total_output_rec = 0 ! Total flux output recordings so far
      ! Suffix to append to simulation output rootname
      character (len=32) :: file_name = '_surf_flux.nc'


      public wrt_surface_flux_his

      contains
! ----------------------------------------------------------------------

      subroutine set_field_corr(istr,iend,jstr,jend)
      ! Set surface fields that will be restored towards

      use read_write, only: set_frc_var_tile

      implicit none

      ! input/outputs
      integer,intent(in)  :: istr,iend,jstr,jend

!#if !defined ANA_STFLUX
      ! Sea-surface temperature (SST) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sst_name,       sst_time_name,  ! Netcdf names for var and time
     &          sst_data,       sst,            ! Input data and var arrays
     &          sst_times,      sst_grd_type,   ! Input times & Grid type
     &          sst_file_indx,  sst_irec,       ! File indx & current input record
     &          sst_it1,        sst_it2 )       ! Indices to input data and time

!#if !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_var_tile(
     &          istr,iend,      jstr,jend,      ! Tile bounds indices
     &          sss_name,       sss_time_name,  ! Netcdf names for var and time
     &          sss_data,       sss,            ! Input data and var arrays
     &          sss_times,      sss_grd_type,   ! Input times & Grid type
     &          sss_file_indx,  sss_irec,       ! File indx & current input record
     &          sss_it1,        sss_it2 )       ! Indices to input data and time

      end subroutine set_field_corr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine wrt_surface_flux_his
      ! Write surface fluxes to output netcdf file

      use netcdf ! for nf90_routines
      use read_write, only: ncdf_create_file, nc_define_var, nc_write_var
     &   , nc_write_time, output_root_name, r_var, u_var, v_var

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
          file_name = trim(output_root_name) // trim(file_name)
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
            call nc_define_var(ncid, sustr_name,
     &              'wind stress in x-direction','N/m^2', u2dgrd, ierr )
          endif
          if (wrt_svstr) then
            call nc_define_var(ncid, svstr_name,
     &             'wind stress in y-direction', 'N/m^2', v2dgrd, ierr )
          endif
          if (wrt_shflx) then
            call nc_define_var(ncid, shflx_name,
     &                   'Surface heat flux', 'degC m/s', r2dgrd, ierr )
          endif
          if (wrt_ssflx) then
            call nc_define_var(ncid, ssflx_name,
     &                 'Surface salinity flux','PSU m/s', r2dgrd, ierr )
          endif

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
          call nc_write_var(ncid, sustr, 1, sustr_name, u_var, output_rec, ierr)
        endif
        if (wrt_svstr) then
          call nc_write_var(ncid, svstr, 1, svstr_name, v_var, output_rec, ierr)
        endif
        if (wrt_shflx) then
          call nc_write_var(ncid, stflx(GLOBAL_2D_ARRAY,itemp), 1,
     &                      shflx_name, r_var, output_rec, ierr)
        endif
        if (wrt_ssflx) then
          call nc_write_var(ncid, stflx(GLOBAL_2D_ARRAY,isalt), 1,
     &                      ssflx_name, r_var, output_rec, ierr)
        endif

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

#endif /* BULK_FLUX for whole module */

      end module surf_flux
