      module surf_flux
      ! Declaration of surface flux variables
      ! Output of surface fluxes
      ! Sets fields for flux correction 

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020)

#include "cppdefs.opt"

      use param
      use read_write, only: set_frc_var_tile, ncvs
      use roms_read_write
      use nc_read_write
      use netcdf
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


      ! netcdf outputting:
      integer :: ncid = -1, prev_fill_mode
      real    :: t_avg_surf=0
      integer :: navg_surf = 0                              ! number of samples in average


      ! Public functions
      public init_arrays_surf_flx
      public set_surf_field_corr
      public wrt_his_surf_flux
      public wrt_avg_surf_flux

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

      if (wrt_file_avg_surf_flux) then
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

! ----------------------------------------------------------------------
      subroutine calc_avg_surf_flux  ![
      ! set surf_flux averaged variables
      ! taken from set_avg.F in old code
      implicit none

      ! local
      real :: coef

      navg_surf = navg_surf +1
      coef = 1./navg_surf

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'surf_flux :: started averaging. ',
     &   'output_period_avg (s) =', output_period_avg
      endif

      t_avg_surf = t_avg_surf*(1-coef) + time*coef

      ! need i0:i1 indices because arrays still GLOBAL_2D therefore wasted margin
      if (wrt_sustr) sustr_avg( 1:i1,j0:j1) = sustr_avg( 1:i1,j0:j1)*(1-coef) + sustr( 1:i1,j0:j1)*coef
      if (wrt_svstr) svstr_avg(i0:i1, 1:j1) = svstr_avg(i0:i1, 1:j1)*(1-coef) + svstr(i0:i1, 1:j1)*coef
      if (wrt_shflx) shflx_avg(i0:i1,j0:j1) = shflx_avg(i0:i1,j0:j1)*(1-coef) + stflx(i0:i1,j0:j1,itemp)*coef
      if (wrt_ssflx) ssflx_avg(i0:i1,j0:j1) = ssflx_avg(i0:i1,j0:j1)*(1-coef) + stflx(i0:i1,j0:j1,isalt)*coef

      end subroutine calc_avg_surf_flux  !]

!----------------------------------------------------------------------
      subroutine def_vars_surf_flux( ncid, avg )  ![
      ! define output variable & attributes in netcdf results file
      ! for instantaneous or averaged variables
      implicit none

      ! input
      integer,intent(in) :: ncid
      logical,intent(in) :: avg
      ! local
      integer           :: ierr, varid
      character(len=64) :: text_lname

      ! output surface flux as per Eq.Sys. units m^2/s^2 not N/m^2
      if (wrt_sustr) then
        varid = nccreate(ncid,sustr_name,(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/))
        if (.not. avg) text_lname='wind stress in x-direction'
        if (      avg) text_lname='averaged '/ /'wind stress in x-direction'
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
      endif
      if (wrt_svstr) then
        varid = nccreate(ncid,svstr_name,(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/))
        if (.not. avg) text_lname='wind stress in y-direction'
        if (      avg) text_lname='averaged '/ /'wind stress in y-direction'
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
      endif
      if (wrt_shflx) then
        varid = nccreate(ncid,shflx_name,(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname='Surface heat flux'
        if (      avg) text_lname='averaged '/ /'Surface heat flux'
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units','degC m/s')
      endif
      if (wrt_ssflx) then
        varid = nccreate(ncid,ssflx_name,(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
        if (.not. avg) text_lname='Surface salinity flux'
        if (      avg) text_lname='averaged '/ /'Surface salinity flux'
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units','PSU m/s')
      endif

      end subroutine def_vars_surf_flux  !]

!----------------------------------------------------------------------
      subroutine wrt_his_surf_flux  ![
      ! write surface fluxes to output netcdf file
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_his                              ! current file output record
      integer,save           :: total_rec_his=0                      ! total his output records so far
      real,save              :: output_time_his=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_his
      integer                :: tile, ierr

      if (.not. first_step) then
        output_time_his = output_time_his + dt                       ! only start count after first timestep
      else
        first_step=.false.                                           ! as first step the surf_flux values are for t=0
      endif

      if (output_time_his>=output_period_his .or.
     &    output_time_his==0                      ) then             ! catch to get t=0 output

        if (mod(total_rec_his,nrpf_his)==0) then
          call create_file_surf_flux(total_rec_his,fname_his,.false.)
          rec_his = 0
          if (output_time_his==0) call display_output_settings_to_terminal( fname_his )
        endif
        total_rec_his = total_rec_his +1
        rec_his = rec_his + 1

        ierr=nf90_open(fname_his,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)     ! set fill value - nf90_nofill for optimized writing

        call ncwrite(ncid,'ocean_time',(/time/),(/rec_his/))

        start=1; start(3)=rec_his
        if (wrt_sustr) call ncwrite(ncid,sustr_name, sustr( 1:i1,j0:j1),start)
        if (wrt_svstr) call ncwrite(ncid,svstr_name, svstr(i0:i1, 1:j1),start)
        if (wrt_shflx) call ncwrite(ncid,shflx_name, stflx(i0:i1,j0:j1,itemp),start)
        if (wrt_ssflx) call ncwrite(ncid,ssflx_name, stflx(i0:i1,j0:j1,isalt),start)

        ierr=nf90_close (ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'surf_flux :: wrote history, tdays =', tdays,
     &     'step =', iic-1, 'rec =', rec_his, '/', total_rec_his     ! -1 since iic=t(n+1)
        endif

        output_time_his=0
      endif

      end subroutine wrt_his_surf_flux  !]

! ----------------------------------------------------------------------
      subroutine wrt_avg_surf_flux  ![
      ! write averaged variables to output netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_avg                              ! current file output record
      integer,save           :: total_rec_avg=0                      ! total avg output records so far
      real,save              :: output_time_avg=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_avg
      integer :: tile, tn, ierr

      if (.not. first_step) then
        call calc_avg_surf_flux                                      ! don't include t=0 in averaging
        output_time_avg = output_time_avg + dt                       ! only start count after first timestep
      endif
      first_step=.false.                                             ! as first step the surf_flux values are for t=0

      if (output_time_avg>=output_period_avg) then                   ! catch to get t=0 output

        if (mod(total_rec_avg,nrpf_avg)==0) then
          call create_file_surf_flux(total_rec_avg,fname_avg,.true.)
          rec_avg = 0
        endif
        total_rec_avg = total_rec_avg +1
        rec_avg = rec_avg + 1

        ierr=nf90_open(fname_avg,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/t_avg_surf/),(/rec_avg/))

        start=1; start(3)=rec_avg                                    ! back to 2D vars
        if (wrt_sustr) call ncwrite(ncid,sustr_name, sustr_avg( 1:i1,j0:j1),start)
        if (wrt_svstr) call ncwrite(ncid,svstr_name, svstr_avg(i0:i1, 1:j1),start)
        if (wrt_shflx) call ncwrite(ncid,shflx_name, shflx_avg(i0:i1,j0:j1),start)
        if (wrt_ssflx) call ncwrite(ncid,ssflx_name, ssflx_avg(i0:i1,j0:j1),start)

        navg_surf=0
        output_time_avg=0

        ierr=nf90_close(ncid)
        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'surf_flux :: wrote averages, tdays =', tdays,
     &        'step =', iic-1, 'rec =', rec_avg, '/', total_rec_avg  ! -1 since iic=t(n+1)
     &         MYID
        endif
      endif

      end subroutine wrt_avg_surf_flux  !]

!----------------------------------------------------------------------
      subroutine create_file_surf_flux(total_rec,fname,avg)  ![
      implicit none

      !input/output
      integer          ,intent(in)  :: total_rec           ! total netcdf records so far for these variables
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ierr,varid
      character(len=10),dimension(4) :: dimnames           ! dimension names
      integer,          dimension(4) :: dimsizes

      if (      avg) fname=trim(output_root_name) / / '_surf_flux_avg.nc'
      if (.not. avg) fname=trim(output_root_name) / / '_surf_flux.nc'
      call create_filename(fname,total_rec)

      ierr=nf90_create(fname,nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) call handle_ierr(ierr,'create_file_surf_flux: ')

      ! ADD 4 DEFAULT DIMENSIONS SO NCJOIN WILL WORK!!!

      varid = nccreate(ncid,'ocean_time',(/'time'/),(/0/),nf90_double)
      if (      avg) ierr = nf90_put_att(ncid,varid,'long_name','Average time since origin')
      if (.not. avg) ierr = nf90_put_att(ncid,varid,'long_name','Time since origin')
      ierr = nf90_put_att(ncid,varid,'units','second' )

      call put_global_atts(ncid, ierr)                     ! put global attributes in file
      if (ierr/=nf90_noerr)
     &  call handle_ierr(ierr,'create_file_surf_flux: ')

      if (      avg) call def_vars_surf_flux( ncid, .true.  )
      if (.not. avg) call def_vars_surf_flux( ncid, .false. )

      ierr = nf90_enddef(ncid)
      if (mynode == 0) then
        write(*,'(7x,2A)')
     &    'surf_flux_frc :: created new netCDF file ', trim(fname)
      endif

      end subroutine create_file_surf_flux !]

!-----------------------------------------------------------------------
      subroutine display_output_settings_to_terminal( fname )  ![
      implicit none

      character(len=99) :: fname

      if (mynode==0) write(*,'(/7x,2A,F6.1,2x,A,I4,2x,2A)')          ! write to terminal output in simulation pre-amble text which
     &  'surf_flux:: net flux history file.',                        ! result variables are being stored
     &  'ouput_period =', output_period_his,
     &  'recs/file =', nrpf_his,
     &  'file = ', trim(fname)

      if (mynode==0) write(*,'(9x,A,4(/11x,A,T20,L1),/)')            ! 4(....) repeats formatting 4 times.
     &               '2D his/avg fields to be saved: (T/F)'          ! T20 moves to the 20th character on line
     &     ,sustr_name,  wrt_sustr
     &     ,svstr_name,  wrt_svstr
     &     ,shflx_name,  wrt_shflx
     &     ,ssflx_name,  wrt_ssflx

      end subroutine display_output_settings_to_terminal  !]

!-----------------------------------------------------------------------

      end module surf_flux
