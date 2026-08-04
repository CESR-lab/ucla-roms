
module basic_output
  ! Basic output routines for ocean variables
  ! his, avg, and rst files

#include "cppdefs.opt"
  use namelist_open_mod, only: open_namelist_file
  use netcdf, only:&
  &nf90_double, nf90_int, nf90_global, nf90_write,&
  &nf90_inq_varid, nf90_put_var, nf90_put_att, nf90_open, nf90_close,&
  &nf90_nofill, nf90_set_fill, nf90_redef, nf90_enddef,&
  &nf90_create, nf90_clobber, nf90_64bit_data, nf90_noerr
  use nc_read_write, only: nccreate, ncwrite
  use dimensions, only:  eta_rho, eta_v, xi_rho, xi_u&
  &, ds_xr, ds_yr, ds_zr, ds_xu, ds_yv, ds_zw
  use roms_read_write, only:&
  &max_name_size, dn_tm, dn_xr, dn_xu,&
  &dn_yr, dn_yv, dn_zr, dn_zw,&
  &avg_output, create_file, sec2date, output_root_name,&
  &refdatestr, put_global_atts, append_date_node
  use sponge_tune, only:  wrt_rst_ub, ub_tune
  use calc_pflx_mod, only: wrt_rst_diag_slow, calc_pflx
  use scalars, only:&
  &NT, iic, knew, nstp, time,&
  &nz, dt, start_time, nnew, tdays
  use ocean_vars, only:&
  &zeta_avg, ubar_avg, vbar_avg, u_avg, v_avg,&
  &w_avg, wvl_avg, zeta, ubar, vbar, u, v, z_r, We, Wi
  use param, only: Lm, Mm, isalt, itemp, mynode, nt_passive, ocean_grid_comm, nt_cdr_oae, nt_cdr_dor
  use error_handling_mod, only: error_log
  use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type
  use pio, only : PIO_openfile, PIO_closefile, PIO_write
#endif
  use mpi_f08, only: MPI_CHARACTER, mpi_bcast

  implicit none
  private

  real(kind=8), public :: output_period_rst = 100        ! output period in seconds
  real(kind=8)    :: nrpf_rst          = 10

  real(kind=8)    :: output_period_his = 100         ! output period in seconds
  integer(kind=4) :: nrpf_his          = 20          ! total recs per file

  real(kind=8)    :: output_period_avg = 100         ! output averaging period in seconds
  integer(kind=4) :: nrpf_avg          = 20          ! total recs per file

  logical, public ::&
  &wrt_file_rst, wrt_file_avg, wrt_file_his, monthly_restarts,&
  &wrt_Z, wrt_Ub, wrt_Vb, wrt_U, wrt_V, wrt_R, wrt_O, wrt_W,&
  &wrt_Akv, wrt_Akt, wrt_aks, wrt_Hbls, wrt_Hbbl,&
  &wrt_avg_Z, wrt_avg_Ub, wrt_avg_Vb, wrt_avg_U, wrt_avg_V,&
  &wrt_avg_R, wrt_avg_O, wrt_avg_W, wrt_avg_Akv, wrt_avg_Akt,&
  &wrt_avg_Aks, wrt_avg_Hbls, wrt_avg_Hbbl

  namelist /BASIC_OUTPUT_SETTINGS/&
  &output_period_rst, nrpf_rst,&
  &output_period_his, nrpf_his,&
  &output_period_avg, nrpf_avg,&
  &wrt_file_rst, wrt_file_avg, wrt_file_his, monthly_restarts,&
  &wrt_Z, wrt_Ub, wrt_Vb, wrt_U, wrt_V, wrt_R, wrt_O, wrt_W,&
  &wrt_Akv, wrt_Akt, wrt_aks, wrt_Hbls, wrt_Hbbl,&
  &wrt_avg_Z, wrt_avg_Ub, wrt_avg_Vb, wrt_avg_U, wrt_avg_V,&
  &wrt_avg_R, wrt_avg_O, wrt_avg_W, wrt_avg_Akv, wrt_avg_Akt,&
  &wrt_avg_Aks, wrt_avg_Hbls, wrt_avg_Hbbl


  character(len=12) :: module_name = "basic_output"
# ifdef SOLVE3D
  real(kind=8),allocatable,dimension(:,:,:) :: akv_avg
  real(kind=8),allocatable,dimension(:,:,:) :: akt_avg
#  ifdef SALINITY
  real(kind=8),allocatable,dimension(:,:,:) :: aks_avg
#  endif
#  ifdef LMD_KPP
  real(kind=8),allocatable,dimension(:,:) :: hbl_avg
#  endif
#  ifdef LMD_BKPP
  real(kind=8),allocatable,dimension(:,:) :: hbbl_avg
#  endif
  real(kind=8),allocatable,dimension(:,:,:) :: rho_avg
# endif

  ![   Taken from old ncvars (should be improved eventually):
  integer(kind=4), parameter, public :: &
#ifdef SOLVE3D
  & indxU=5, indxV=6, indxO=7, indxW=8&
  &, indxR=9, indxT=10,&
# ifdef SALINITY
  &indxS=indxT+1, &
# endif
  & indxTime=1, indxZ=2, indxUb=3, indxVb=4

  integer(kind=4), public :: &
# ifdef SALINITY
  & indxAks, &
# endif
# ifdef LMD_KPP
#  ifdef SALINITY
  & indxHbls, &
#  else
  & indxHbls, &
#  endif
# ifdef LMD_BKPP
  & indxHbbl, &
# endif
# else
# ifdef LMD_BKPP
  & indxHbbl, &
# endif
# endif
  & indxAkt, indxAkv
#endif
  !]

  character(len=42), public ::  vname(3, 180) ! DevinD 180 copied from old code only
  integer(kind=4), parameter, public :: iaux=6           ! Length of netCDF variable "time_step"

  ! end old ncvars

  ! netcdf outputting:
  integer(kind=4) :: ncid=-1, prev_fill_mode
  real(kind=8)    :: t_avg_ovars=0
  integer(kind=4) :: navg_ovars = 0                             ! number of samples in average
  character(len=max_name_size),public :: fname_rst                 ! restart filename needed by diagnostics
  integer(kind=4),public         :: rec_rst=2          ! current file output record needed by diags
  integer(kind=4)                :: month_at_prev_timestep
  logical                :: write_another_step = .false.

  public :: wrt_his_ocean_vars
  public :: wrt_avg_ocean_vars
  public :: calc_avg_ocean_vars
  public :: wrt_rst_ocean_vars
  public :: init_avg_arrays
  public :: init_restarts
  public :: read_nml_basic_output

contains                  !]

  subroutine read_nml_basic_output

!     Read the "BASIC_OUTPUT_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_basic_output"
    character(len=512) :: msg = ""

    ! initialize indxAkv
    indxAkv = indxT + NT

#ifdef SOLVE3D
    indxAkv = indxT + NT
    indxAkt = indxAkv + 1
# ifdef SALINITY
    indxAks = indxAkt + 1
# endif
# ifdef LMD_KPP
#  ifdef SALINITY
    indxHbls = indxAks + 1
#  else
    indxHbls = indxAkt + 1
#  endif
# ifdef LMD_BKPP
    indxHbbl = indxHbls + 1
# endif
# else
# ifdef LMD_BKPP
    indxHbbl = indxAks + 1
# endif
# endif
#endif

    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=BASIC_OUTPUT_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read BASIC_OUTPUT_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)

  end subroutine read_nml_basic_output

!----------------------------------------------------------------------
  subroutine init_restarts  ![

    implicit none

    integer(kind=4),dimension(6)   :: date

    if (wrt_file_rst .and. monthly_restarts) then
      ! Look one timestep ahead to initialize the "previous month", so that
      ! when we restart we don't immediately write out a file.
      call sec2date(start_time+dt,date)
      month_at_prev_timestep = date(2)
    endif

  end subroutine init_restarts !]
!----------------------------------------------------------------------
  subroutine init_avg_arrays  ![
    implicit none
    character(len=15) :: sr_name = "init_avg_arrays"
    character(len=1024) :: error_info


    if (wrt_file_his .and. mod(output_period_his, dt) /= 0) then
      write(error_info, *)&
      &"dt (", dt, ") is not a factor of ",&
      &"output_period_his (",&
      &output_period_his, ")"

      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif

    if (wrt_file_avg .and. mod(output_period_avg,dt) /= 0) then
      write(error_info, *)&
      &"dt (", dt, ") is not a factor of ",&
      &"output_period_avg (",&
      &output_period_avg, ")"
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif
    call error_log%abort_check()

    if (wrt_file_avg) then
      if (wrt_avg_Z)  then
        allocate( zeta_avg(GLOBAL_2D_ARRAY) )
        zeta_avg=0._8                            ! avg needs to be 0. because Nan x 0. = Nan in set_avg.F
      endif
      if (wrt_avg_Ub) then
        allocate( ubar_avg(GLOBAL_2D_ARRAY) )
        ubar_avg=0._8
      endif
      if (wrt_avg_Vb) then
        allocate( vbar_avg(GLOBAL_2D_ARRAY) )
        vbar_avg=0._8
      endif

      if (wrt_avg_U) then
        allocate( u_avg(GLOBAL_2D_ARRAY,nz) )
        u_avg=0._8
      endif
      if (wrt_avg_V) then
        allocate( v_avg(GLOBAL_2D_ARRAY,nz) )
        v_avg=0._8
      endif
      if (wrt_avg_O) then
        allocate( w_avg(GLOBAL_2D_ARRAY,0:nz) )
        w_avg=0
      endif
      if (wrt_avg_W) then
        allocate( wvl_avg(GLOBAL_2D_ARRAY,0:nz) )
        wvl_avg=0._8
      endif

      ! the following variables do not 'live' in this module, but in order to prevent
      ! circular reference cause by e.g. wrt_R logical, they are allocated here.
      if (wrt_avg_R) then
        allocate( rho_avg(GLOBAL_2D_ARRAY,nz) )
        rho_avg=0._8
      endif
      if (wrt_avg_Akv)  allocate( akv_avg(GLOBAL_2D_ARRAY,0:nz) )
      if (wrt_avg_Akt)  allocate( akt_avg(GLOBAL_2D_ARRAY,0:nz) )
# ifdef SALINITY
      if (wrt_avg_Aks)  allocate( aks_avg(GLOBAL_2D_ARRAY,0:nz) )
# endif
# ifdef LMD_KPP
      if (wrt_avg_Hbls) allocate( hbl_avg(GLOBAL_2D_ARRAY) )
# endif
# ifdef LMD_BKPP
      if (wrt_avg_Hbbl) allocate( hbbl_avg(GLOBAL_2D_ARRAY) )
# endif

    endif  ! <-- wrt_file_avg

! Names of variables in NetCDF output files. The first element
! is the name of the variable; the other two are are attributes.

    vname(1,indxTime)='ocean_time'
    vname(2,indxTime)='Time since initialization'
    vname(3,indxTime)='second'

    vname(1,indxZ)='zeta'
    vname(2,indxZ)='free-surface elevation'
    vname(3,indxZ)='meter'

    vname(1,indxUb)='ubar'
    vname(2,indxUb)='vertically averaged u-momentum component'
    vname(3,indxUb)='meter second-1'

    vname(1,indxVb)='vbar'
    vname(2,indxVb)='vertically averaged v-momentum component'
    vname(3,indxVb)='meter second-1'

#ifdef SOLVE3D
    vname(1,indxU)='u'
    vname(2,indxU)='u-momentum component'
    vname(3,indxU)='meter second-1'

    vname(1,indxV)='v'
    vname(2,indxV)='v-momentum component'
    vname(3,indxV)='meter second-1'

    vname(1,indxO)='omega'
    vname(2,indxO)='S-coordinate vertical velocity'
    vname(3,indxO)='meter second-1'

    vname(1,indxW)='w'
    vname(2,indxW)='vertical velocity'
    vname(3,indxW)='meter second-1'

    vname(1,indxR)='rho'
    vname(2,indxR)='density anomaly'
    vname(3,indxR)='kilogram meter-3'

    vname(1,indxT)='temp'
    vname(2,indxT)='potential temperature'
    vname(3,indxT)='Celsius'

# ifdef SALINITY
    vname(1,indxS)='salt'
    vname(2,indxS)='salinity'
    vname(3,indxS)='PSU'
# endif

    vname(1,indxAkv)='Akv'
    vname(2,indxAkv)='vertical viscosity coefficient'
    vname(3,indxAkv)='meter2 second-1'

    vname(1,indxAkt)='Akt'
    vname(2,indxAkt)='vertical thermal conductivity coefficient'
    vname(3,indxAkt)='meter2 second-1'
# ifdef SALINITY
    vname(1,indxAks)='AKs'
    vname(2,indxAks)='salinity vertical diffusion coefficient'
    vname(3,indxAks)='meter2 second-1'
# endif
# ifdef LMD_KPP
    vname(1,indxHbls)='hbls'
    vname(2,indxHbls)='Thickness of KPP surface boundary layer'
    vname(3,indxHbls)='meter'
# endif
# ifdef LMD_BKPP
    vname(1,indxHbbl)='hbbl'
    vname(2,indxHbbl)='Thickness of KPP bottom boundary layer'
    vname(3,indxHbbl)='meter'
# endif
#endif

    ! NOTE: call this here to prevent output in between timesteps which would break code
! checking script.
    if (wrt_file_rst) call display_output_settings_to_terminal_rst
    if (wrt_file_his) call display_output_settings_to_terminal_his
    if (wrt_file_avg) call display_output_settings_to_terminal_avg

  end subroutine init_avg_arrays !]
! ----------------------------------------------------------------------
  subroutine wrt_his_ocean_vars(special)  ![
    ! write ocean_vars variables to output netcdf file
    ! ocean_vars variables are calculated for t=n in timestep t=n
    ! (unlike u/v/etc which are calculated for t=n+1 in timestep t=n)

    use mixing, only:&
#ifdef LMD_KPP
    &hbls,&
#endif
#ifdef LMD_BKPP
    & hbbl, &
#endif
    & Akv, Akt
#ifdef SPLIT_EOS
    use eos_vars, only: rho1
#else
    use eos_vars, only: rho
#endif
    use work_mod, only: work
    use grid, only: pn, pm
    use dimensions, only: i0,i1,j0,j1
    use tracers, only: wrt_his_trc
    use wvlcty_mod, only: wvlcty

    implicit none

    ! import/export
    logical,optional       :: special !! for blowup and initial history file
    ! local
    integer(kind=4),dimension(4)   :: start
    integer(kind=4),save           :: rec_his                     ! current file output record
    integer(kind=4)                :: total_rec_his=0                      ! total his output records so far
    real(kind=8),save              :: output_time_his=0
    logical,save           :: first_step=.true.
    character(len=max_name_size),save :: fname_his
    integer(kind=4)                :: tile, ierr, i, j, k
    logical,save           :: supress=.false.
    logical                :: blowup,initial

    if (first_step) then
      output_time_his=output_period_his+10
      rec_his = nrpf_his
      first_step = .false.
    end if

    blowup = .false.
    initial= .false.
    if (present(special)) then
      if (special) then
        blowup = .true.
      else
        initial = .true.
      endif
    endif


    if (wrt_file_his.or.blowup) then

      if (.not.present(special)) then
        output_time_his = output_time_his + dt
      endif


      if (output_time_his>=output_period_his.or.present(special) ) then

        if (initial) then
          if (mynode==0) print *,'Writing initial history file'
        elseif (blowup) then
          if (mynode==0) print *,'Writing emergency history file'
        endif

        if (rec_his==nrpf_his) then
#ifdef PARALLEL_IO
          rec_his = 0
          if (mynode == 0) then
            call create_file_ocean_vars(fname_his,.false.)     !     call create_file('_rnd',fname, nonode=.true.)
!          ierr=nf90_open(fname,nf90_write,ncid)
!          call def_vars_random(ncid)
!          ierr = nf90_close(ncid)
          endif
          call MPI_Bcast(fname_his,max_name_size,MPI_CHARACTER,0,ocean_grid_comm,ierr)
          call MPI_Barrier(ocean_grid_comm, ierr)
        endif

        rec_his = rec_his + 1

        if (mynode == 0) then
          ierr=nf90_open(fname_his,nf90_write,ncid)
          ! always add time
          call ncwrite(ncid,'ocean_time',(/time/),(/rec_his/))
          call write_time_step( ncid, rec_his, total_rec_his)
          ierr=nf90_close (ncid)
        endif
        ! abort_check uses MPI collectives; must not run only on rank 0
        call error_log%abort_check()
        call MPI_Barrier(ocean_grid_comm, ierr)

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname_his), PIO_write)

        start=1; start(3)=rec_his                                    ! back to 2D vars
        pio_gtype = '2Drw'
        if (wrt_Z)  call ncwrite(ncid,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start, .true.)
        pio_gtype = '2Duw'
        if (wrt_Ub) call ncwrite(ncid,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start, .true.)
        pio_gtype = '2Dvw'
        if (wrt_Vb) call ncwrite(ncid,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start, .true.)

#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_his called at the middle of next timestep as some variables only calculated there for t=n.

        start(3)=1; start(4)=rec_his
        pio_gtype = '3Duw'
        if (wrt_U) call ncwrite(ncid,vname(1,indxU),u( 1:i1,j0:j1,:,nnew),start, .true.)
        pio_gtype = '3Dvw'
        if (wrt_V) call ncwrite(ncid,vname(1,indxV),v(i0:i1, 1:j1,:,nnew),start, .true.)

        pio_gtype = '3Drw'
        call wrt_his_trc(ncid,start)                              ! tracer variables

        if (wrt_R) then
# ifdef SPLIT_EOS
          pio_gtype = '3Drw'
          call ncwrite(ncid,vname(1,indxR),rho1(i0:i1,j0:j1,:),start, .true.)
# else
          pio_gtype = '3Drw'
          call ncwrite(ncid,vname(1,indxR), rho(i0:i1,j0:j1,:),start, .true.)
# endif
        endif
        if (wrt_O) then                              ! s-coordinate omega vertical velocity (m/s).
          do k=0,nz
            do j=0,Mm+1
              do i=0,Lm+1
                work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
                ! Implicit omega only
!               work(i,j,k)=pm(i,j)*pn(i,j)*Wi(i,j,k)
              enddo
            enddo
          enddo
          pio_gtype = '3Dww'
          call ncwrite(ncid,vname(1,indxO),work(i0:i1,j0:j1,:),start, .true.)
        endif
        if (wrt_W) then                              ! true vertical velocity (m/s).
!         do tile=0,NSUB_X*NSUB_E-1
!           call wvlcty (tile, work)
!         enddo
          call wvlcty (0, work)
          pio_gtype = '3Drw'
          call ncwrite(ncid,vname(1,indxW),work(i0:i1,j0:j1,1:nz),start, .true.)
        endif
        pio_gtype = '3Dww'
        if (wrt_Akv) call ncwrite(ncid,vname(1,indxAkv),Akv(i0:i1,j0:j1,:),start, .true.)
        pio_gtype = '3Dww'
        if (wrt_Akt) call ncwrite(ncid,vname(1,indxAkt),Akt(i0:i1,j0:j1,:,itemp),start, .true.)
# ifdef SALINITY
        pio_gtype = '3Dww'
        if (wrt_Aks) call ncwrite(ncid,vname(1,indxAks),Akt(i0:i1,j0:j1,:,isalt),start, .true.)
# endif
        start(3)=rec_his                                      ! back to 2D vars
# ifdef LMD_KPP
        pio_gtype = '2Drw'
        if (wrt_Hbls) call ncwrite(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1),start, .true.)
# endif
# ifdef LMD_BKPP
        pio_gtype = '2Drw'
        if (wrt_Hbbl) call ncwrite(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1),start, .true.)
# endif
#endif
        output_time_his=0
        call PIO_closefile(pio_FileDesc)
      endif

#else   ! PARALLEL_IO

          call create_file_ocean_vars(fname_his,.false.)
          rec_his = 0
        endif
        rec_his = rec_his + 1

        ierr=nf90_open(fname_his,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/rec_his/))

        call write_time_step( ncid, rec_his, total_rec_his )

        start=1; start(3)=rec_his                                    ! back to 2D vars
        if (wrt_Z)  call ncwrite(ncid,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start)
        if (wrt_Ub) call ncwrite(ncid,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start)
        if (wrt_Vb) call ncwrite(ncid,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start)

#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_his called at the middle of next timestep as some variables only calculated there for t=n.

        start(3)=1; start(4)=rec_his
        if (wrt_U) call ncwrite(ncid,vname(1,indxU),u( 1:i1,j0:j1,:,nnew),start)
        if (wrt_V) call ncwrite(ncid,vname(1,indxV),v(i0:i1, 1:j1,:,nnew),start)

        call wrt_his_trc(ncid,start)                              ! tracer variables

        if (wrt_R) then
# ifdef SPLIT_EOS
          call ncwrite(ncid,vname(1,indxR),rho1(i0:i1,j0:j1,:),start)
# else
          call ncwrite(ncid,vname(1,indxR), rho(i0:i1,j0:j1,:),start)
# endif
        endif
        if (wrt_O) then                              ! s-coordinate omega vertical velocity (m/s).
          do k=0,nz
            do j=0,Mm+1
              do i=0,Lm+1
                work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
                ! Implicit omega only
!               work(i,j,k)=pm(i,j)*pn(i,j)*Wi(i,j,k)
              enddo
            enddo
          enddo
          call ncwrite(ncid,vname(1,indxO),work(i0:i1,j0:j1,:),start)
        endif
        if (wrt_W) then                              ! true vertical velocity (m/s).
          call wvlcty (0, work)
          call ncwrite(ncid,vname(1,indxW),work(i0:i1,j0:j1,1:nz),start)
        endif
        if (wrt_Akv) call ncwrite(ncid,vname(1,indxAkv),Akv(i0:i1,j0:j1,:),start)
        if (wrt_Akt) call ncwrite(ncid,vname(1,indxAkt),Akt(i0:i1,j0:j1,:,itemp),start)
# ifdef SALINITY
        if (wrt_Aks) call ncwrite(ncid,vname(1,indxAks),Akt(i0:i1,j0:j1,:,isalt),start)
# endif
        start(3)=rec_his                                      ! back to 2D vars
# ifdef LMD_KPP
        if (wrt_Hbls) call ncwrite(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1),start)
# endif
# ifdef LMD_BKPP
        if (wrt_Hbbl) call ncwrite(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1),start)
# endif
#endif
        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&
          &'ocean_vars :: wrote history, tdays =', tdays,&
          &'step =', iic, 'rec =', rec_his
        endif

        output_time_his=0
      endif

#endif ! PARALLEL_IO

    endif  ! <-- wrt_file_his

  end subroutine wrt_his_ocean_vars  !]
! ----------------------------------------------------------------------
  subroutine wrt_avg_ocean_vars  ![
    ! write averaged variables to output netcdf file
    ! don't include t=0 in averaging. This create 0.5dt error in averaging,
    ! but this 0.5dt error has always been in ROMS.
    ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5_8*t2, but we've never done that.

    use grid, only: pm, pn
    use dimensions, only: i0, i1, j0, j1
    use tracers, only: wrt_avg_trc
    use roms_read_write, only: store_string_att
    implicit none

    ! local
    integer(kind=4),dimension(4)   :: start
    integer(kind=4)                :: rec_avg
    integer(kind=4),save           :: total_rec_avg=0                      ! total avg output records so far
    real(kind=8),save              :: output_time_avg=0                    ! time since last output
    character(len=max_name_size),save :: fname_avg
    character(len=99)      :: output_time_string
    character(len=99)      :: formatted_string
    integer(kind=4) :: tile, tn, ierr, k
    logical, save          :: first_step

    if (first_step) then
      rec_avg = nrpf_avg
      first_step = .false.
    end if

    if (wrt_file_avg) then

      output_time_avg = output_time_avg + dt

      if (output_time_avg>=output_period_avg) then

        ! save info to netcdf
        avg_output = ''
        write(output_time_string,'(F12.1)') output_time_avg
        write(formatted_string,'(A,A)') trim(adjustl(output_time_string)),&
        &' seconds'
        call store_string_att(avg_output,formatted_string)

        if (rec_avg == nrpf_avg) then
#ifdef PARALLEL_IO
          rec_avg = 0
          if (mynode == 0) then
            pio_gtype = '3Drw' ! So it initializes the variables with the correct size
            call create_file_ocean_vars(fname_avg,.true.)
          endif
          call MPI_Bcast(fname_avg,max_name_size,MPI_CHARACTER,0,ocean_grid_comm,ierr)
          call MPI_Barrier(ocean_grid_comm, ierr)
        endif
        rec_avg = rec_avg + 1

        if (mynode == 0) then
          ierr=nf90_open(fname_avg,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/t_avg_ovars/),(/rec_avg/))
          ierr=nf90_close (ncid)
        endif

        call MPI_Barrier(ocean_grid_comm, ierr)

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname_avg), PIO_write)

        start=1; start(3)=rec_avg                                    ! back to 2D vars
        pio_gtype = '2Drw'
        if (wrt_avg_Z)   call ncwrite(ncid, vname(1,indxZ),  zeta_avg(i0:i1,j0:j1), start,.true.)
        pio_gtype = '2Duw'
        if (wrt_avg_Ub)  call ncwrite(ncid, vname(1,indxUb), ubar_avg( 1:i1,j0:j1), start,.true.)
        pio_gtype = '2Dvw'
        if (wrt_avg_Vb)  call ncwrite(ncid, vname(1,indxVb), vbar_avg(i0:i1, 1:j1), start,.true.)
# ifdef SOLVE3D
        start(3)=1; start(4)=rec_avg
        pio_gtype = '3Duw'
        if (wrt_avg_U)   call ncwrite(ncid, vname(1,indxU),     u_avg( 1:i1,j0:j1,:), start,.true.)
        pio_gtype = '3Dvw'
        if (wrt_avg_V)   call ncwrite(ncid, vname(1,indxV),     v_avg(i0:i1, 1:j1,:), start,.true.)

        pio_gtype = '3Drw'
        call wrt_avg_trc(ncid,start)     ! STILL NEED TO UPDATE                                   ! tracer variables
        pio_gtype = '3Drw'
        if (wrt_avg_R)   call ncwrite(ncid, vname(1,indxR),   rho_avg(i0:i1,j0:j1,:), start,.true.)
        if (wrt_avg_O) then
          do k=0,nz
            w_avg(i0:i1,j0:j1,k) = w_avg(i0:i1,j0:j1,k)*pm(i0:i1,j0:j1)*pn(i0:i1,j0:j1)  ! convert before write.
          enddo
          pio_gtype = '3Dww'
          call ncwrite(ncid, vname(1,indxO),     w_avg(i0:i1,j0:j1,:), start,.true.)            ! here rather than calc_
        endif
        pio_gtype = '3Dww'                                                                 ! for efficiency
        if (wrt_avg_W)   call ncwrite(ncid, vname(1,indxW),   wvl_avg(i0:i1,j0:j1,1:nz), start,.true.)
        pio_gtype = '3Dww'
        if (wrt_avg_Akv) call ncwrite(ncid, vname(1,indxAkv), akv_avg(i0:i1,j0:j1,:), start,.true.)
        pio_gtype = '3Dww'
        if (wrt_avg_Akt) call ncwrite(ncid, vname(1,indxAkt), akt_avg(i0:i1,j0:j1,:), start,.true.)
#  ifdef SALINITY
        pio_gtype = '3Dww'
        if (wrt_avg_Aks) call ncwrite(ncid, vname(1,indxAks), aks_avg(i0:i1,j0:j1,:), start,.true.)
#  endif
        start(3)=rec_avg                                                ! back to 2D vars
#  ifdef LMD_KPP
        pio_gtype = '2Drw'
        if (wrt_avg_Hbls) call ncwrite(ncid, vname(1,indxHbls),  hbl_avg(i0:i1,j0:j1), start,.true.)
#  endif
#  ifdef LMD_BKPP
        pio_gtype = '2Drw'
        if (wrt_avg_Hbbl) call ncwrite(ncid, vname(1,indxHbbl), hbbl_avg(i0:i1,j0:j1), start,.true.)
#  endif
# endif /* SOLVE3D */

        navg_ovars=0
        output_time_avg=0

        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&  ! confirm work completed
          &'ocean_vars :: wrote averages, tdays =', tdays,&
          &'step =', iic, 'rec =', rec_avg, '/', total_rec_avg
        endif
        call PIO_closefile(pio_FileDesc)
      endif
#else ! PARALLEL_IO
          call create_file_ocean_vars(fname_avg,.true.)
          rec_avg = 0
        endif
        rec_avg = rec_avg + 1

        ierr=nf90_open(fname_avg,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/t_avg_ovars/),(/rec_avg/))

        start=1; start(3)=rec_avg                                    ! back to 2D vars
        if (wrt_avg_Z)   call ncwrite(ncid, vname(1,indxZ),  zeta_avg(i0:i1,j0:j1), start)
        if (wrt_avg_Ub)  call ncwrite(ncid, vname(1,indxUb), ubar_avg( 1:i1,j0:j1), start)
        if (wrt_avg_Vb)  call ncwrite(ncid, vname(1,indxVb), vbar_avg(i0:i1, 1:j1), start)
# ifdef SOLVE3D
        start(3)=1; start(4)=rec_avg
        if (wrt_avg_U)   call ncwrite(ncid, vname(1,indxU),     u_avg( 1:i1,j0:j1,:), start)
        if (wrt_avg_V)   call ncwrite(ncid, vname(1,indxV),     v_avg(i0:i1, 1:j1,:), start)

        call wrt_avg_trc(ncid,start)     ! STILL NEED TO UPDATE                                   ! tracer variables

        if (wrt_avg_R)   call ncwrite(ncid, vname(1,indxR),   rho_avg(i0:i1,j0:j1,:), start)
        if (wrt_avg_O) then
          do k=0,nz
            w_avg(i0:i1,j0:j1,k) = w_avg(i0:i1,j0:j1,k)*pm(i0:i1,j0:j1)*pn(i0:i1,j0:j1)  ! convert before write.
          enddo
          call ncwrite(ncid, vname(1,indxO),     w_avg(i0:i1,j0:j1,:), start)            ! here rather than calc_avg
        endif                                                                            ! for efficiency
        if (wrt_avg_W)   call ncwrite(ncid, vname(1,indxW),   wvl_avg(i0:i1,j0:j1,1:nz), start)
        if (wrt_avg_Akv) call ncwrite(ncid, vname(1,indxAkv), akv_avg(i0:i1,j0:j1,:), start)
        if (wrt_avg_Akt) call ncwrite(ncid, vname(1,indxAkt), akt_avg(i0:i1,j0:j1,:), start)
#  ifdef SALINITY
        if (wrt_avg_Aks) call ncwrite(ncid, vname(1,indxAks), aks_avg(i0:i1,j0:j1,:), start)
#  endif
        start(3)=rec_avg                                                ! back to 2D vars
#  ifdef LMD_KPP
        if (wrt_avg_Hbls) call ncwrite(ncid, vname(1,indxHbls),  hbl_avg(i0:i1,j0:j1), start)
#  endif
#  ifdef LMD_BKPP
        if (wrt_avg_Hbbl) call ncwrite(ncid, vname(1,indxHbbl), hbbl_avg(i0:i1,j0:j1), start)
#  endif
# endif /* SOLVE3D */

        navg_ovars=0
        output_time_avg=0

        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&  ! confirm work completed
          &'ocean_vars :: wrote averages, tdays =', tdays,&
          &'step =', iic, 'rec =', rec_avg, '/', total_rec_avg
        endif
      endif

#endif ! PARALLEL_IO

    endif  ! <-- wrt_file_avg

  end subroutine wrt_avg_ocean_vars  !]
!----------------------------------------------------------------------
  subroutine wrt_rst_ocean_vars  ![
    ! write ocean_vars variables to output netcdf restart file
    ! ocean_vars variables are calculated for t=n in timestep t=n
    ! (unlike u/v/etc which are calculated for t=n+1 in timestep t=n)
    implicit none

    ! local
    real(kind=8),save              :: output_time_rst=0                    ! time since last output
    integer(kind=4),dimension(6)   :: date
    character(len=15)      :: datestr

    if (wrt_file_rst) then

      if (monthly_restarts) then
        ! We need to save the last two timesteps of the month, so that when the new simulation
        ! loads the restart file the first timestep it takes will be the first timestep of the
        ! new month.
        ! time + 2*dt here so that it looks two timesteps ahead to see if we're changing months.
        call sec2date(time+2*dt,date)

        if (write_another_step) then
          call wrt_restart_file
          write_another_step = .false.
        endif

        if ((date(2) - month_at_prev_timestep) /= 0) then
          call wrt_restart_file
          write_another_step = .true.
        endif

        month_at_prev_timestep = date(2)
      else
!     if (.not. first_step) output_time_rst = output_time_rst + dt   ! only start count after first timestep
!     first_step=.false.                                             ! as first step the ocean_vars values are for t=0
        output_time_rst = output_time_rst + dt

        if (output_time_rst    >=output_period_rst&
#ifdef EXACT_RESTART
        &.or.&
        &output_time_rst+dt >=output_period_rst&         ! timestep before
#endif
        &) then
          call wrt_restart_file
          if (output_time_rst>=output_period_rst) then
            output_time_rst=0                             ! catch for exact restart to not
          endif                                             ! reset the time
        endif ! it output_time > output_period
      endif ! monthly_restarts

    endif ! wrt_file_rst

  end subroutine wrt_rst_ocean_vars  !]
!----------------------------------------------------------------------
  subroutine wrt_restart_file  ![
    ! write ocean_vars variables to output netcdf restart file
    ! ocean_vars variables are calculated for t=n in timestep t=n
    ! (unlike u/v/etc which are calculated for t=n+1 in timestep t=n)
#if defined LMD_KPP
    use mixing, only: hbls, hbbl
#endif
    use dimensions, only: i0, i1, j0,j1
    use tracers, only: wrt_rst_trc
    use coupling, only:&
    &du_avg1, dv_avg1, du_avg2, dv_avg2,&
    &du_avg_bak, dv_avg_bak
    use grid, only: riv_umask, riv_vmask
#ifdef MARBL
!     add MARBL saved state to restart to restart file
    use marbl_driver, only: marbldrv_write_ss_vars_to_rst
#endif

    implicit none

    ! local
    integer(kind=4),dimension(4)   :: start
    integer(kind=4),save           :: total_rec_rst=0                      ! total his output records so far
    real(kind=8),save              :: output_time_rst=0                    ! time since last output
    logical,save           :: first_step=.true.
    integer(kind=4)                :: tile, ierr, i, j, k

#ifdef PARALLEL_IO
    if (rec_rst==2) then
#ifdef EXACT_RESTART
      time = time + dt
#endif
      call create_file_rst_ocean_vars(fname_rst)
#ifdef EXACT_RESTART
      time = time - dt
#endif
      rec_rst = 0
    endif
    total_rec_rst = total_rec_rst +1
    rec_rst = rec_rst + 1

    if (mynode == 0) then
      ierr=nf90_open(fname_rst,nf90_write,ncid)

      call ncwrite(ncid,'ocean_time',(/time/),(/rec_rst/))

      call write_time_step( ncid, rec_rst, total_rec_rst)
      ierr=nf90_close (ncid)
    endif
    call error_log%abort_check()
    call MPI_Barrier(ocean_grid_comm, ierr)

    ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname_rst), PIO_write)

    start=1; start(3)=rec_rst                                    ! back to 2D vars
    pio_gtype = '2Drw'
    call ncwrite(ncid,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start,.true.)
    pio_gtype = '2Duw'
    call ncwrite(ncid,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start,.true.)
    pio_gtype = '2Dvw'
    call ncwrite(ncid,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start,.true.)
#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_rst called at the middle of next timestep as some variables only calculated there for t=n.

    start(3)=1; start(4)=rec_rst
    pio_gtype = '3Duw'
    call ncwrite(ncid,vname(1,indxU),u( 1:i1,j0:j1,:,nnew),start,.true.)
    pio_gtype = '3Dvw'
    call ncwrite(ncid,vname(1,indxV),v(i0:i1, 1:j1,:,nnew),start,.true.)
    pio_gtype = '3Drw'
    call wrt_rst_trc(ncid,start)                              ! tracer variables

    start(3)=rec_rst                                      ! back to 2D vars

# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
    pio_gtype = '2Duw'
    call ncwrite(ncid,'DU_avg1',    DU_avg1( 1:i1,j0:j1),start,.true.)
    pio_gtype = '2Dvw'
    call ncwrite(ncid,'DV_avg1',    DV_avg1(i0:i1, 1:j1),start,.true.)
    pio_gtype = '2Duw'
    call ncwrite(ncid,'DU_avg2',    DU_avg2( 1:i1,j0:j1),start,.true.)
    pio_gtype = '2Dvw'
    call ncwrite(ncid,'DV_avg2',    DV_avg2(i0:i1, 1:j1),start,.true.)
    pio_gtype = '2Duw'
    call ncwrite(ncid,'DU_avg_bak', DU_avg_bak( 1:i1,j0:j1),start,.true.)
    pio_gtype = '2Dvw'
    call ncwrite(ncid,'DV_avg_bak', DV_avg_bak(i0:i1, 1:j1),start,.true.)

#  elif defined PRED_COUPLED_MODE
    pio_gtype = '2Duw'
    call ncwrite(ncid,'rufrc',rufrc_bak( 1:i1,j0:j1,nstp),start,.true.)
    pio_gtype = '2Dvw'
    call ncwrite(ncid,'rvfrc',rvfrc_bak(i0:i1, 1:j1,nstp),start,.true.)
#  endif
# endif


# ifdef LMD_KPP
    pio_gtype = '2Drw'
    call ncwrite(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1),start,.true.)
# endif
# ifdef LMD_BKPP
    pio_gtype = '2Drw'
    call ncwrite(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1),start,.true.)
# endif
#endif
    ! Added for perfect restart reasons
!       start(3)=1; start(4)=rec_rst
!       call ncwrite(ncid,'Akv',Akv(i0:i1,j0:j1,:),start)
!       call ncwrite(ncid,'Akt',Akt(i0:i1,j0:j1,:,itemp),start)
!       call ncwrite(ncid,'Aks',Akt(i0:i1,j0:j1,:,isalt),start)

#ifdef SPONGE_TUNE
    if (ub_tune)   call wrt_rst_ub(ncid,rec_rst)
#endif
    if (calc_pflx) call wrt_rst_diag_slow(ncid,rec_rst)

    pio_gtype = '2Duw'
    call ncwrite(ncid,'riv_umask',    riv_umask( 1:i1,j0:j1),start,.true.)
    pio_gtype = '2Dvw'
    call ncwrite(ncid,'riv_vmask',    riv_vmask(i0:i1, 1:j1),start,.true.)

#ifdef MARBL
    pio_gtype = '3Drw'
    call marbldrv_write_ss_vars_to_rst(ncid,rec_rst)
#endif

    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&
      &'ocean_vars :: wrote restart, tdays =', tdays,&
      &'step =', iic, 'rec =', rec_rst, '/', total_rec_rst
    endif

    call PIO_closefile(pio_FileDesc)

#else ! PARALLEL_IO

    if (rec_rst==2) then
#ifdef EXACT_RESTART
      time = time + dt
#endif
      call create_file_rst_ocean_vars(fname_rst)
#ifdef EXACT_RESTART
      time = time - dt
#endif
      rec_rst = 0
    endif
    total_rec_rst = total_rec_rst +1
    rec_rst = rec_rst + 1

    ierr=nf90_open(fname_rst,nf90_write,ncid)
    ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

    call ncwrite(ncid,'ocean_time',(/time/),(/rec_rst/))

    call write_time_step( ncid, rec_rst, total_rec_rst )

    start=1; start(3)=rec_rst                                    ! back to 2D vars
    call ncwrite(ncid,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start)
    call ncwrite(ncid,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start)
    call ncwrite(ncid,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start)
#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_rst called at the middle of next timestep as some variables only calculated there for t=n.

    start(3)=1; start(4)=rec_rst
    call ncwrite(ncid,vname(1,indxU),u( 1:i1,j0:j1,:,nnew),start)
    call ncwrite(ncid,vname(1,indxV),v(i0:i1, 1:j1,:,nnew),start)

    call wrt_rst_trc(ncid,start)                              ! tracer variables

    start(3)=rec_rst                                      ! back to 2D vars

# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
    call ncwrite(ncid,'DU_avg1',    DU_avg1( 1:i1,j0:j1),start)
    call ncwrite(ncid,'DV_avg1',    DV_avg1(i0:i1, 1:j1),start)
    call ncwrite(ncid,'DU_avg2',    DU_avg2( 1:i1,j0:j1),start)
    call ncwrite(ncid,'DV_avg2',    DV_avg2(i0:i1, 1:j1),start)
    call ncwrite(ncid,'DU_avg_bak', DU_avg_bak( 1:i1,j0:j1),start)
    call ncwrite(ncid,'DV_avg_bak', DV_avg_bak(i0:i1, 1:j1),start)

#  elif defined PRED_COUPLED_MODE
    call ncwrite(ncid,'rufrc',rufrc_bak( 1:i1,j0:j1,nstp),start)
    call ncwrite(ncid,'rvfrc',rvfrc_bak(i0:i1, 1:j1,nstp),start)
#  endif
# endif


# ifdef LMD_KPP
    call ncwrite(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1),start)
# endif
# ifdef LMD_BKPP
    call ncwrite(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1),start)
# endif
#endif
    ! Added for perfect restart reasons
!       start(3)=1; start(4)=rec_rst
!       call ncwrite(ncid,'Akv',Akv(i0:i1,j0:j1,:),start)
!       call ncwrite(ncid,'Akt',Akt(i0:i1,j0:j1,:,itemp),start)
!       call ncwrite(ncid,'Aks',Akt(i0:i1,j0:j1,:,isalt),start)

#ifdef SPONGE_TUNE
    if (ub_tune)   call wrt_rst_ub(ncid,rec_rst)
#endif
    if (calc_pflx) call wrt_rst_diag_slow(ncid,rec_rst)

    call ncwrite(ncid,'riv_umask',    riv_umask( 1:i1,j0:j1),start)
    call ncwrite(ncid,'riv_vmask',    riv_vmask(i0:i1, 1:j1),start)

#ifdef MARBL
    call marbldrv_write_ss_vars_to_rst(ncid,rec_rst)
#endif

    ierr=nf90_close(ncid)
    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&
      &'ocean_vars :: wrote restart, tdays =', tdays,&
      &'step =', iic, 'rec =', rec_rst, '/', total_rec_rst
    endif

#endif ! PARALLEL_IO

  end subroutine wrt_restart_file  !]
!----------------------------------------------------------------------
  subroutine create_file_ocean_vars(fname,avg)  ![
    implicit none

    !input/output
    character(len=max_name_size),intent(out) :: fname
    logical,          intent(in)  :: avg                 ! his or average file

    ! local
    integer(kind=4) :: ierr,varid

#ifdef PARALLEL_IO
    if (avg) then
      fname=trim(adjustl(output_root_name)) // '_avg'
    else
      fname=trim(adjustl(output_root_name)) // '_his'
    endif
    call append_date_node(fname, nonode=.true.)
    ierr = nf90_create(trim(fname), ior(nf90_clobber, nf90_64bit_data), ncid)
    if (ierr/=nf90_noerr) then
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/create_file_ocean_vars",&
      &info="unable to create file "//trim(fname))
    end if
    call put_global_atts(ncid, ierr)
    varid = nccreate(ncid,'ocean_time',(/dn_tm/),(/0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name', refdatestr)
    ierr = nf90_put_att(ncid,varid,'units','second' )
    if (avg) then
      call def_vars_ocean_vars(  .true.  )
      ierr=nf90_put_att(ncid,nf90_global,'type','ROMS averages file')
      ierr=nf90_put_att(ncid,nf90_global,'averaging_timescale',adjustl(avg_output))
    else
      call def_vars_ocean_vars(  .false. )
      ierr=nf90_put_att(ncid,nf90_global,'type','ROMS history file')
    endif
    ierr=nf90_enddef(ncid)
    if (ierr/=nf90_noerr) then
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/create_file_ocean_vars",&
      &info="enddef for file "//trim(fname))
    end if
    ierr=nf90_close(ncid)
    if (ierr/=nf90_noerr) then
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/create_file_ocean_vars",&
      &info="close for file "//trim(fname))
    end if
    if (mynode == 0) then
      write(*,'(7x,2A)') 'created new netcdf file ', trim(fname)
    endif
#else
    if (avg) then
      call create_file('_avg',fname)
    else
      call create_file('_his',fname)
    endif

    ierr=nf90_open(fname,nf90_write,ncid)
    ierr=nf90_redef(ncid)

    if (avg) then
      call def_vars_ocean_vars(  .true.  )
      ierr=nf90_put_att(ncid,nf90_global,'type','ROMS averages file')
      ierr=nf90_put_att(ncid,nf90_global,'averaging_timescale',adjustl(avg_output))
    else
      call def_vars_ocean_vars(  .false. )
      ierr=nf90_put_att(ncid,nf90_global,'type','ROMS history file')
    endif

    ierr=nf90_enddef(ncid)
    ierr=nf90_close(ncid)
#endif

  end subroutine create_file_ocean_vars !]
!----------------------------------------------------------------------
  subroutine create_file_rst_ocean_vars(fname)  ![
    implicit none

    !input/output
    character(len=max_name_size),intent(out) :: fname

    ! local
    integer(kind=4) :: ierr,varid

#ifdef PARALLEL_IO
    if (mynode == 0) then
      fname=trim(adjustl(output_root_name)) // '_rst'
      call append_date_node(fname, nonode=.true.)
      ierr = nf90_create(trim(fname), ior(nf90_clobber, nf90_64bit_data), ncid)
      if (ierr/=nf90_noerr) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//"/create_file_rst_ocean_vars",&
        &info="unable to create file "//trim(fname))
      end if
      call put_global_atts(ncid, ierr)
      ierr=nf90_put_att(ncid,nf90_global,'type','ROMS restart file')
      call def_vars_rst_ocean_vars
      ierr=nf90_enddef(ncid)
      if (ierr/=nf90_noerr) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//"/create_file_rst_ocean_vars",&
        &info="enddef for file "//trim(fname))
      end if
      ierr=nf90_close(ncid)
      if (ierr/=nf90_noerr) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//"/create_file_rst_ocean_vars",&
        &info="close for file "//trim(fname))
      end if
      write(*,'(7x,2A)') 'created new netcdf file ', trim(fname)
    endif
    call error_log%abort_check()
    call MPI_Bcast(fname,max_name_size,MPI_CHARACTER,0,ocean_grid_comm,ierr)
    call MPI_Barrier(ocean_grid_comm, ierr)
#else
    call create_file('_rst',fname)

    ierr=nf90_open(fname,nf90_write,ncid)
    ierr=nf90_redef(ncid)

    ierr=nf90_put_att(ncid,nf90_global,'type','ROMS restart file')

    call def_vars_rst_ocean_vars

    ierr=nf90_enddef(ncid)
    ierr=nf90_close(ncid)
#endif

  end subroutine create_file_rst_ocean_vars !]
! ----------------------------------------------------------------------
  subroutine def_vars_ocean_vars( avg )  ![
    ! define output variable & attributes in netcdf results file
    ! for instantaneous or averaged variables

    use tracers, only: def_avg_trc, def_his_trc
    implicit none

    ! input
    logical,intent(in) :: avg
    ! local
    integer(kind=4)           :: ierr, varid
    character(len=64) :: text_lname
    character(len=7)  :: dn_aux = 'auxil'

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
    varid = nccreate(ncid,'time_step',(/dn_aux,dn_tm/),(/iaux,0/),nf90_int)
    ierr=nf90_put_att (ncid, varid, 'long_name',&
    &'time step and record numbers from initialization')

    if ( (wrt_Z .and. .not. avg) .or. (wrt_avg_Z .and. avg) ) then ! .or. needed for his or avg output selection
      if (.not. avg) text_lname=vname(2,indxZ)
      if (      avg) text_lname='averaged '//vname(2,indxZ)
      varid = nccreate(ncid,vname(1,indxZ),(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxZ))
    endif
    if ( (wrt_Ub .and. .not. avg) .or. (wrt_avg_Ub .and. avg) ) then
      if (.not. avg) text_lname=vname(2,indxUb)
      if (      avg) text_lname='averaged '//vname(2,indxUb)
      varid = nccreate(ncid,vname(1,indxUb),(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxUb))
    endif
    if ( (wrt_Vb .and. .not. avg) .or. (wrt_avg_Vb .and. avg) ) then
      if (.not. avg) text_lname=vname(2,indxVb)
      if (      avg) text_lname='averaged '//vname(2,indxVb)
      varid = nccreate(ncid,vname(1,indxVb),(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxVb))
    endif
#ifdef SOLVE3D
    if ( (wrt_U .and. .not. avg) .or. (wrt_avg_U .and. avg) ) then
      if (.not. avg) text_lname=vname(2,indxU)
      if (      avg) text_lname='averaged '//vname(2,indxU)
      varid = nccreate(ncid,vname(1,indxU),(/dn_xu,dn_yr,dn_zr,dn_tm/),(/ds_xu,ds_yr,ds_zr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxU))
    endif
    if ( (wrt_V .and. .not. avg) .or. (wrt_avg_V .and. avg) ) then
      if (.not. avg) text_lname=vname(2,indxV)
      if (      avg) text_lname='averaged '//vname(2,indxV)
      varid = nccreate(ncid,vname(1,indxV),(/dn_xr,dn_yv,dn_zr,dn_tm/),(/ds_xr,ds_yv,ds_zr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxV))
    endif

    if (avg) then
      call def_avg_trc( ncid )
    else
      call def_his_trc( ncid )
    endif

    if ( (wrt_R .and. .not. avg) .or. (wrt_avg_R .and. avg) ) then
      if (.not. avg) text_lname=vname(2,indxR)
      if (      avg) text_lname='averaged '//vname(2,indxR)
      varid =&
      &nccreate(ncid,vname(1,indxR),(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxR))
    endif
    if ( (wrt_O .and. .not. avg) .or. (wrt_avg_O .and. avg) ) then    ! s-coordinate "omega" vertical velocity.
      if (.not. avg) text_lname=vname(2,indxO)
      if (      avg) text_lname='averaged '//vname(2,indxO)
      varid =&
      &nccreate(ncid,vname(1,indxO),(/dn_xr,dn_yr,dn_zw,dn_tm/),(/ds_xr,ds_yr,ds_zw,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxO))
    endif
    if ( (wrt_W .and. .not. avg) .or. (wrt_avg_W .and. avg) ) then   ! true W-vertical velocity.
      if (.not. avg) text_lname=vname(2,indxW)
      if (      avg) text_lname='averaged '//vname(2,indxW)
      varid =&
      &nccreate(ncid,vname(1,indxW),(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxW))
    endif
    if ( (wrt_Akv .and. .not. avg) .or. (wrt_avg_Akv .and. avg)) then ! vertical viscosity coefficient.
      if (.not. avg) text_lname=vname(2,indxAkv)
      if (      avg) text_lname='averaged '//vname(2,indxAkv)
      varid = nccreate(ncid,vname(1,indxAkv),&
      &(/dn_xr,dn_yr,dn_zw,dn_tm/),(/ds_xr,ds_yr,ds_zw,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkv))
    endif
    if ( (wrt_Akt .and. .not. avg) .or. (wrt_avg_Akt .and. avg) ) then  ! vertical thermal conductivity coefficient.
      if (.not. avg) text_lname=vname(2,indxAkt)
      if (      avg) text_lname='averaged '//vname(2,indxAkt)
      varid = nccreate(ncid,vname(1,indxAkt),&
      &(/dn_xr,dn_yr,dn_zw,dn_tm/),(/ds_xr,ds_yr,ds_zw,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkt))
    endif
# ifdef SALINITY
    if ( (wrt_Aks .and. .not. avg) .or. (wrt_avg_Aks .and. avg) ) then! vertical diffusion coefficient for salinity.
      if (.not. avg) text_lname=vname(2,indxAks)
      if (      avg) text_lname='averaged '//vname(2,indxAks)
      varid =&
      &nccreate(ncid,vname(1,indxAks),(/dn_xr,dn_yr,dn_zw,dn_tm/),(/ds_xr,ds_yr,ds_zw,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAks))
    endif
# endif /* SALINITY */
# ifdef LMD_KPP
    if ( (wrt_Hbls .and. .not. avg) .or. (wrt_avg_Hbls .and. avg) ) then                                   ! depth of
      if (.not. avg) text_lname=vname(2,indxHbls)
      if (      avg) text_lname='averaged '//vname(2,indxHbls)
      varid = nccreate(ncid,vname(1,indxHbls),(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbls))
    endif
# endif
# ifdef LMD_BKPP
    if ( (wrt_Hbbl .and. .not. avg) .or. (wrt_avg_Hbbl .and. avg) ) then       ! thickness of bottom boundary layer.
      if (.not. avg) text_lname=vname(2,indxHbbl)
      if (      avg) text_lname='averaged '//vname(2,indxHbbl)
      varid = nccreate(ncid,vname(1,indxHbbl),&
      &(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbbl))
    endif
# endif
#endif /* SOLVE3D */

  end subroutine def_vars_ocean_vars  !]
! ----------------------------------------------------------------------
  subroutine def_vars_rst_ocean_vars  ![
    ! define output variable & attributes in netcdf restart file
    use tracers, only: t_vname, t_lname, t_units
#ifdef MARBL
!     add MARBL saved state to restart to restart file
    use marbl_driver, only: marbldrv_create_ss_vars_in_rst
#endif

    implicit none

    ! local
    integer(kind=4) :: ierr, varid, itrc
    character(len=7)  :: dn_aux = 'auxil'

#ifdef PARALLEL_IO
    varid = nccreate(ncid,'ocean_time',(/dn_tm/),(/0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name', refdatestr)
    ierr = nf90_put_att(ncid,varid,'units','second' )
#endif

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
    varid = nccreate(ncid,'time_step',(/dn_aux,dn_tm/),(/iaux,0/),nf90_int)

    ierr=nf90_put_att (ncid, varid, 'long_name',&
    &'time step and record numbers from initialization')

    varid = nccreate(ncid,vname(1,indxZ),(/dn_xr,dn_yr,dn_tm/),&
    &(/ds_xr,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxZ))
    ierr = nf90_put_att(ncid,varid,'units',vname(3,indxZ))

    varid = nccreate(ncid,vname(1,indxUb),(/dn_xu,dn_yr,dn_tm/),&
    &(/ds_xu,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxUb))
    ierr = nf90_put_att(ncid,varid,'units',vname(3,indxUb))

    varid = nccreate(ncid,vname(1,indxVb),(/dn_xr,dn_yv,dn_tm/),&
    &(/ds_xr,ds_yv,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxVb))
    ierr = nf90_put_att(ncid,varid,'units',vname(3,indxVb))
#ifdef MARBL
    call marbldrv_create_ss_vars_in_rst(ncid)
#endif
#ifdef SOLVE3D
    varid = nccreate(ncid,vname(1,indxU),(/dn_xu,dn_yr,dn_zr,dn_tm/),&
    &(/ds_xu,ds_yr,nz,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxU))
    ierr = nf90_put_att(ncid,varid,'units',vname(3,indxU))

    varid = nccreate(ncid,vname(1,indxV),(/dn_xr,dn_yv,dn_zr,dn_tm/),&
    &(/ds_xr,ds_yv,nz,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxV))
    ierr = nf90_put_att(ncid,varid,'units',vname(3,indxV))

    do itrc=1,NT
      varid = nccreate(ncid,t_vname(itrc),(/dn_xr,dn_yr,dn_zr,dn_tm/),&
      &(/ds_xr,ds_yr,nz,0/), nf90_double)
      ierr=nf90_put_att (ncid, varid, 'long_name', t_lname(itrc))
      ierr=nf90_put_att (ncid, varid, 'units', t_units(itrc))
    enddo

# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
    varid = nccreate(ncid,'DU_avg1',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'<<fast-time averaged uflx>>')

    varid = nccreate(ncid,'DV_avg1',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'<<fast-time-averaged vflx>>')

    varid = nccreate(ncid,'DU_avg2',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'<<fast-time averaged ubar(:,:,n+1/2)>>')

    varid = nccreate(ncid,'DV_avg2',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'<<fast-time-averaged vbar(:,:,n+1/2)>>')

    varid = nccreate(ncid,'DU_avg_bak',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'<back-step mixed fast-time-averaged ubar(:,:,n-1)>')

    varid = nccreate(ncid,'DV_avg_bak',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'<back-step mixed fast-time-averaged vbar(:,:,n-1)>')
#  elif defined PRED_COUPLED_MODE
    varid = nccreate(ncid,'rufrc_bak',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'3D-to-2D forcing, XI-component')
    varid = nccreate(ncid,'rvfrc_bak',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'3D-to-2D forcing, ETA-component')
#  endif
# endif  /* EXACT_RESTART */

# ifdef LMD_KPP
    varid = nccreate(ncid,vname(1,indxHbls),&
    &(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxHbls))
    ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbls))
# endif
# ifdef LMD_BKPP
    varid = nccreate(ncid,vname(1,indxHbbl),(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxHbbl))
    ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbbl))
# endif

    if (calc_pflx) then
      varid = nccreate(ncid,'u_slow',(/dn_xu,dn_yr,dn_zr,dn_tm/),&
      &(/ds_xu,ds_yr,nz,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','time filtered u')
      ierr = nf90_put_att(ncid,varid,'units','m/s')

      varid = nccreate(ncid,'v_slow',(/dn_xr,dn_yv,dn_zr,dn_tm/),&
      &(/ds_xr,ds_yv,nz,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','time filtered v')
      ierr = nf90_put_att(ncid,varid,'units','m/s')

      varid = nccreate(ncid,'p_slow',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/),&
      &(/ds_xr,ds_yr,nz,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',&
      &'time filtered pressure')
      ierr = nf90_put_att(ncid,varid,'units','Pa??')
    endif

#ifdef SPONGE_TUNE
    if (ub_tune) then
      varid = nccreate(ncid,'ub_west',&
      &(/dn_yr,dn_tm/),&
      &(/ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',&
      &'west bc tuning coeff')
      varid = nccreate(ncid,'ub_east',&
      &(/dn_yr,dn_tm/),&
      &(/ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',&
      &'east bc tuning coeff')
      varid = nccreate(ncid,'ub_north',&
      &(/dn_xr,dn_tm/),&
      &(/ds_xr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',&
      &'north bc tuning coeff')
      varid = nccreate(ncid,'ub_south',&
      &(/dn_xr,dn_tm/),&
      &(/ds_xr,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',&
      &'south bc tuning coeff')
    endif
#endif
    varid = nccreate(ncid,'riv_umask',(/dn_xu,dn_yr,dn_tm/),&
    &(/ds_xu,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'river mask at u points')
    ierr = nf90_put_att(ncid,varid,'units','nondim')

    varid = nccreate(ncid,'riv_vmask',(/dn_xr,dn_yv,dn_tm/),&
    &(/ds_xr,ds_yv,0/),nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'river mask at v points')
    ierr = nf90_put_att(ncid,varid,'units','nondim')

#endif /* SOLVE3D */

  end subroutine def_vars_rst_ocean_vars  !]
! ----------------------------------------------------------------------
  subroutine calc_avg_ocean_vars  ![
! set averaged variables
#ifdef SPLIT_EOS
    use eos_vars, only: rho1, qp1
#else
    use eos_vars, only: rho
#endif
    use work_mod, only: work
    use mixing, only: &
#ifdef LMD_KPP
    & hbls, hbbl, &
#endif
    & akv, akt
    use private_scratch, only: a2d
    use dimensions, only: i0, i1, j0, j1, nx, ny
    use tracers, only: set_avg_trc
    use wvlcty_mod, only: wvlcty_tile
    implicit none

    ! local
    real(kind=8)         :: coef
!     logical,save :: first_step=.true.

    if (wrt_file_avg) then

!     if (.not. first_step) then                           ! don't include t=0 in average

      navg_ovars = navg_ovars +1
      coef = 1._8/navg_ovars

      if (coef==1) then                                  ! this refreshes average (1-coef)=0
        if (mynode==0) write(*,'(7x,2A,F9.1)')&
        &'ocean_vars :: started averaging. ',&
        &'output_period_avg (s) =', output_period_avg
      endif

      t_avg_ovars = t_avg_ovars*(1-coef) + time*coef

      ! need i0:i1 indices because arrays still GLOBAL_2D therefore wasted margin
      if (wrt_avg_Z)  zeta_avg(i0:i1,j0:j1) = zeta_avg(i0:i1,j0:j1)*(1-coef) + zeta(i0:i1,j0:j1,knew)*coef
      if (wrt_avg_Ub) ubar_avg( 1:i1,j0:j1) = ubar_avg( 1:i1,j0:j1)*(1-coef) + ubar( 1:i1,j0:j1,knew)*coef
      if (wrt_avg_Vb) vbar_avg(i0:i1, 1:j1) = vbar_avg(i0:i1, 1:j1)*(1-coef) + vbar(i0:i1, 1:j1,knew)*coef
# ifdef SOLVE3D
      if (wrt_avg_U)   u_avg( 1:i1,j0:j1,:) = u_avg( 1:i1,j0:j1,:)*(1-coef) + u( 1:i1,j0:j1,:,nstp)*coef
      if (wrt_avg_V)   v_avg(i0:i1, 1:j1,:) = v_avg(i0:i1, 1:j1,:)*(1-coef) + v(i0:i1, 1:j1,:,nstp)*coef

      call set_avg_trc(coef)                                        ! tracer variables

      if (wrt_avg_R) then
#  ifdef SPLIT_EOS
        rho_avg(i0:i1,j0:j1,:) = rho_avg(i0:i1,j0:j1,:)*(1-coef) +&
        &( rho1(i0:i1,j0:j1,:) - qp1(i0:i1,j0:j1,:) * z_r(i0:i1,j0:j1,:) )*coef
#  else
        rho_avg(i0:i1,j0:j1,:) = rho_avg(i0:i1,j0:j1,:)*(1-coef) + rho(i0:i1,j0:j1,:)*coef
#  endif
      endif
      if (wrt_avg_O) w_avg(i0:i1,j0:j1,:) = w_avg(i0:i1,j0:j1,:)*(1-coef)&
      &+ ( We(i0:i1,j0:j1,:)+Wi(i0:i1,j0:j1,:) ) *coef
      if (wrt_avg_W) then
        ! the w averaging, calculation and output should be checked at some point.
        call wvlcty_tile(1,nx,1,ny, work, A2d(1,1), A2d(1,2), A2d(1,3))  ! as per set_avg

        wvl_avg(i0:i1,j0:j1,1:nz) = wvl_avg(i0:i1,j0:j1,1:nz)*(1-coef)&
        &+ work(i0:i1,j0:j1,1:nz)*coef
      endif
      if (wrt_avg_Akv) akv_avg(i0:i1,j0:j1,:) = akv_avg(i0:i1,j0:j1,:)*(1-coef) + akv(i0:i1,j0:j1,:)*coef
      if (wrt_avg_Akt) akt_avg(i0:i1,j0:j1,:) = akt_avg(i0:i1,j0:j1,:)*(1-coef) + akt(i0:i1,j0:j1,:,itemp)*coef
#  ifdef SALINITY
      if (wrt_avg_Aks) aks_avg(i0:i1,j0:j1,:) = aks_avg(i0:i1,j0:j1,:)*(1-coef) + akt(i0:i1,j0:j1,:,isalt)*coef
#  endif
#  ifdef LMD_KPP
      if (wrt_avg_Hbls) hbl_avg(i0:i1,j0:j1)  =  hbl_avg(i0:i1,j0:j1)*(1-coef) + hbls(i0:i1,j0:j1)*coef
#  endif
#  ifdef LMD_BKPP
      if (wrt_avg_Hbbl) hbbl_avg(i0:i1,j0:j1) = hbbl_avg(i0:i1,j0:j1)*(1-coef) + hbbl(i0:i1,j0:j1)*coef
#  endif
# endif /* SOLVE3D */

!     endif  ! <-- .not. first_step

!     first_step=.false.

    endif ! <-- wrt_file_avg

  end subroutine calc_avg_ocean_vars  !]
! ----------------------------------------------------------------------
  subroutine write_time_step( ncid, record, total_recs)  ![
    implicit none
    character(len=16) :: sr_name = "write_time_step"
    ! inputs
    integer(kind=4), intent(in) :: ncid
    integer(kind=4), intent(in) :: record                        ! current file record number
    integer(kind=4), intent(in) :: total_recs                    ! total records for variable

    ! local
    integer(kind=4) :: ibuff(iaux), start(2), count(2)           ! iaux = 6 from wrt_his.F
    integer(kind=4) :: var_id_tmp, ierr

    ibuff(1)=iic    ; ibuff(2)=record                    ! time step and
    ibuff(4:iaux)=0 ; ibuff(3)=total_recs                ! rechis numbers.

    start(1)=1      ; count(1)=6                         ! iaux = 6 in ncvars
    start(2)=record ; count(2)=1

    ierr=nf90_inq_varid(ncid, 'time_step', var_id_tmp)
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &context=module_name//"/"//sr_name,&
    &info="variable: time_step")

    ierr=nf90_put_var(ncid, var_id_tmp, ibuff, start, count) ! & record time step info
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &context=module_name//"/"//sr_name,&
    &info="variable: time_step")
    ! can't use ncwrite() as ibuff is an integer array not real

  end subroutine write_time_step  !]
!-----------------------------------------------------------------------
  subroutine display_output_settings_to_terminal_his  ![
    use tracers, only: itands, wrt_t, t_vname
    implicit none

    ! local
    integer(kind=4) :: itrc
    if (mynode==0) then
      write(*,'(/7x,2A,F6.1,2x,A,I4)')&
      &'ocean_vars :: history file ',&
      &'output_period =', output_period_his,&
      &'recs/file =', nrpf_his ! write to terminal output in simulation pre-amble text which result variables are being stored

      write(*,'(9x,A)')&
      &'his fields to be saved: (T/F)' ! T20 moves to the 20th character on line
      write(*,'(9x,A)')  repeat('-',62)
      write(*, '(11x,A,T20,A,T36,A)')&
      &"Name","Write (T/F)","Long name"
      write(*,'(9x,A)')  repeat('-',62)
      ! 13(....) repeats formatting 11 times:
      write(*,'(13(/11x,A,T30,L1,T36,A))')&
      &'zeta',   wrt_Z,    vname(2,indxZ)&
      &, 'ubar',   wrt_Ub,   vname(2,indxUb)&
      &, 'vbar',   wrt_Vb,   vname(2,indxVb)&
#ifdef SOLVE3D
      &, 'u',      wrt_U,    vname(2,indxU)&
      &, 'v',      wrt_V,    vname(2,indxV)&

      &, 'rho',    wrt_R,    vname(2,indxR)&
      &, 'Omega',  wrt_O,    vname(2,indxO)&
      &, 'W',      wrt_W,    vname(2,indxW)&

      &, 'Akv',    wrt_Akv,  vname(2,indxAkv)&
      &, 'Akt',    wrt_Akt,  vname(2,indxAkt)&
# ifdef SALINITY
      &, 'Aks',    wrt_Aks,  vname(2,indxAks)&
# endif
# ifdef LMD_KPP
      &, 'hbls',   wrt_Hbls, vname(2,indxHbls)&
# endif
# ifdef LMD_BKPP
      &, 'hbbl',   wrt_Hbbl, vname(2,indxHbbl)&
# endif
      &, ''
      do itrc=1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
      write(*,'(11x,A,I2,A,T30,L1,T36,A)') 't(',&
      &itrc, ')', wrt_t(itrc), t_vname(itrc)
      enddo
      write(*,'(9x,A)')  repeat('-',62)
      end if
#endif /* SOLVE3D */

    end subroutine display_output_settings_to_terminal_his  !]

    subroutine display_output_settings_to_terminal_avg  ![
      use tracers, only: itands, wrt_t_avg, t_vname
      implicit none

      ! local
      integer(kind=4) :: itrc
      if (mynode==0) then
         ! write to terminal output in simulation pre-amble text which result variables are being stored
        write(*,'(/7x,2A,F6.1,2x,A,I4)')&
        &'ocean_vars :: average file ',&
        &'output_period =', output_period_avg,&
        &'recs/file =', nrpf_avg

        write(*,'(9x,A)')&
        &'avg fields to be saved: (T/F)' ! T20 moves to the 20th character on line
        write(*,'(9x,A)')  repeat('-',62)
        write(*, '(11x,A,T20,A,T36,A)')&
        &"Name","Write (T/F)","Long name"
        write(*,'(9x,A)')  repeat('-',62)
        ! 13(....) repeats formatting 11 times.
        write(*,'(13(/11x,A,T30,L1,T36,A))')&
        &'zeta',   wrt_avg_Z,    vname(2,indxZ)&
        &, 'ubar',   wrt_avg_Ub,   vname(2,indxUb)&
        &, 'vbar',   wrt_avg_Vb,   vname(2,indxVb)&
#ifdef SOLVE3D
        &, 'u',      wrt_avg_U,    vname(2,indxU)&
        &, 'v',      wrt_avg_V,    vname(2,indxV)&

        &, 'rho',    wrt_avg_R,    vname(2,indxR)&
        &, 'Omega',  wrt_avg_O,    vname(2,indxO)&
        &, 'W',      wrt_avg_W,    vname(2,indxW)&

        &, 'Akv',    wrt_avg_Akv,  vname(2,indxAkv)&
        &, 'Akt',    wrt_avg_Akt,  vname(2,indxAkt)&
# ifdef SALINITY
        &, 'Aks',    wrt_avg_Aks,  vname(2,indxAks)&
# endif
# ifdef LMD_KPP
        &, 'hbls',   wrt_avg_Hbls, vname(2,indxHbls)&
# endif
# ifdef LMD_BKPP
        &, 'hbbl',   wrt_avg_Hbbl, vname(2,indxHbbl)&
# endif
        &,''
        do itrc=1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
        write(*,'(11x,A,I2,A,T30,L1,T36,A)') 't(',&
        &itrc, ')', wrt_t_avg(itrc), t_vname(itrc)
        enddo
        write(*,'(9x,A)')  repeat('-',62)
        end if ! (mynode==0)
#endif /* SOLVE3D */

      end subroutine display_output_settings_to_terminal_avg  !]

      subroutine display_output_settings_to_terminal_rst
        character(len=120) :: stdout_str
        character(len=120) :: period_info
        if (mynode==0) then
          write(stdout_str,'(7x,A,2x,A,I4)')&
          &'ocean_vars :: restart file ',&
          &'recs/file = 2'
          if (monthly_restarts) then
            write(period_info,'(A,1L)')&
            &'monthly_restarts = ', monthly_restarts
          else
            write(period_info,'(7x,A,F6.1)')&
            &'output_period_rst =', output_period_rst
          end if
          write(*,'(A,2x,A)') trim(stdout_str),trim(period_info)
        end if
      end subroutine display_output_settings_to_terminal_rst

    end module basic_output
