      module surf_flux
      ! Declaration of surface flux variables
      ! Output of surface fluxes
      ! Sets fields for flux correction 

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.opt"

      use param
      use read_write
      use ncvars

      implicit none

      private

#include "surf_flux.opt"

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
#if defined QCORRECTION && !defined ANA_SST
      ! sea-surface temperature (SST) data
      real, public :: sst(GLOBAL_2D_ARRAY)
      real, public :: dSSTdt
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      real, public :: sss(GLOBAL_2D_ARRAY)
      real, public :: dSSSdt                  ! input units (cm/day)
#endif


      ! Netcdf outputting
      ! -----------------
      integer            :: output_rec = 0               ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0         ! total flux output recordings so far
      character (len=60) :: file_name = '_surf_flux.nc'  ! suffix to append to simulation output rootname

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
      call set_frc_var_tile( istr, iend, jstr, jend, nc_sst, sst, interp_frc )
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_var_tile( istr, iend, jstr, jend, nc_sss, sss, interp_frc )
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
      use scalars

      implicit none

      ! local
      integer :: lenstr, my_nf_def_dim
      integer :: lfnm, lvar, prev_fill_mode,
     &           r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)
#ifdef SOLVE3D
     &         , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif

      integer :: ierr = 0        ! This is local so need to set to no-error
      integer :: ncid = -1       ! File ID number - will be set below
      character(len=64) :: fname ! Taken from read_inp.F

      ! Only plot if write_file is true from roms.in or top of module
      ! and timestep is divisble by record rate

      if (write_file==.true. .and. mod(iic-ntstart,rec_rate) == 0) then

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
     &     'surf_flux: save results net_flux = ', write_file,          ! result variables are being stored
     &     'rec_rate =', rec_rate, 'recs/file =', recs_per_file,
     &     'file = ', file_name

          mpi_master_only write(*,'(/1x,A,4(/8x,A,T16,L1),/)')         ! Change to number of variables 4(
     &               'fields to be saved in 2D Flux history: (T/F)'
     &     ,sustr_name,  wrt_sustr
     &     ,svstr_name,  wrt_svstr
     &     ,shflx_name,  wrt_shflx
     &     ,ssflx_name,  wrt_ssflx

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
        if (wrt_ssflx) then
          call nc_write_var(ncid, stflx(GLOBAL_2D_ARRAY,isalt), 1,
     &                      ssflx_name, rp_var, output_rec, ierr)
        endif
!          ! Output of sustr_rho and svstr_rho
!          call nc_write_var(ncid, sustr_r, 1, 'sustr_r', rp_var, output_rec, ierr)
!          call nc_write_var(ncid, svstr_r, 1, 'svstr_r', rp_var, output_rec, ierr)


        ierr=nf90_close (ncid)  ! close netcdf file

        if (mynode == 0) then
          write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &     'flux_his :: wrote history, tdays =', tdays,
     &     'step =', iic-1, 'rec =', output_rec, '/', total_output_rec
     &      MYID
        endif

        if (output_rec == recs_per_file) then  ! reset output_rec to create new file
          output_rec = 0 
        endif

      endif ! (write_file==.true.)

      end subroutine wrt_surface_flux_his
!-----------------------------------------------------------------------


      end module surf_flux
