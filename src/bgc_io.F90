module bgc_io

  ! BGC - bio-geo-chemical module
  ! -----------------------------

  ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Dec)
  ! Refactoring of ETH Zurich roms code with BEC2, code which has been used
  ! by Pierre Damien to run full pacific models.
#include "cppdefs.opt"

#if defined BIOLOGY_BEC2 || defined MARBL

  ! param needed for GLOBAL_2D_array to work. NT = number tracer from param
  use param, only: ieast, itemp, iwest, jnorth, jsouth, ocean_grid_comm, nt_cdr_oae, nt_cdr_dor
  use tracers, only: t_avg, wrt_t_avg, nt_2_t_avg, t_units
#ifdef MARBL_DIAGS
  use marbl_driver, only: marbldrv_compute_init_diagnostics
#endif

! imports from bgc_shared_vars
  use bgc_shared_vars, only: wrt_bgc_his, wrt_bgc_avg,&
       & wrt_bgc_dia_his, wrt_bgc_dia_avg, interp_bgc_frc,&
       & nrpf_bgc_avg, nrpf_bgc_his, nrpf_bgc_avg_dia,&
       & nrpf_bgc_his_dia, t, nc_dust, nc_iron,&
#ifdef PCO2AIR_FORCING
       & nc_xco2air,&
#ifdef MARBL
       & nc_xco2air_alt,&
#endif
#endif
       & bgc_idx, t_idx, lm, mm, itands,&
       & mynode, nt_passive, output_period_bgc_avg, output_period_bgc_his,&
       & t_vname, wrt_t, t_lname
#ifdef NOX_FORCING
  use bgc_shared_vars, only: nc_nox
#endif
#ifdef NHY_FORCING
  use bgc_shared_vars, only: nc_nhy
#endif

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
  use bgc_shared_vars, only: nc_swrad_avg
#endif
#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
  use bgc_shared_vars, only: bgc_diag_2d_avg, bgc_diag_3d_avg,&
       &bgc_diag_2d, bgc_diag_3d, nr_bgc_diag_2d, nr_bgc_diag_3d,&
       &nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d, output_period_bgc_avg_dia,&
       &output_period_bgc_his_dia,&
       &wrt_bgc_diag_2d, idx_bgc_diag_2d, vname_bgc_diag_2d,&
       &wrt_bgc_diag_3d, idx_bgc_diag_3d, vname_bgc_diag_3d
#endif
! /bgc_shared_vars
  use nc_read_write, only: nccreate, ncwrite
  use roms_read_write, only:&
  &dn_tm, dn_xr, dn_yr, dn_zr,&
  &create_file, set_frc_data, max_name_size

  use netcdf, only:&
  &nf90_write, nf90_open, nf90_close, nf90_redef, nf90_enddef,&
  &nf90_double, nf90_put_att
  use scalars, only: nstp, nt, time, nz, dt, iic, nnew, tdays
  use dimensions, only: i0, i1, j0, j1, xi_rho, eta_rho&
  &, ds_xr, ds_yr, ds_zr, ds_xu, ds_yv, ds_zw
  use error_handling_mod, only: error_log
#ifdef PARALLEL_IO
  use pio_roms, only: pio_file_is_open, pio_FileDesc, pio_IoSystem,&
  &pio_type, pio_gtype
  use pio, only : PIO_openfile, PIO_closefile, PIO_write
#endif
  use mpi_f08, only: MPI_CHARACTER, mpi_bcast

  implicit none

  private ! Make all variable private to this module unless public specified

  ! Includes:

! netcdf outputting:
  character(len=6) :: module_name = "bgc_io"
  integer(kind=4) :: ncid=-1, prev_fill_mode
  real(kind=8)    :: t_avg_bgc=0
  real(kind=8)    :: t_avg_dia_bgc=0
  integer(kind=4),save :: navg_bgc = 0               ! number of samples in average
  integer(kind=4),save :: navg_dia_bgc = 0           ! number of samples in average
  integer(kind=4) :: record_avg = 0           ! Triggers making of initial file
  integer(kind=4) :: record_his = 0
  integer(kind=4) :: record_dia_avg = 0   ! Triggers making of initial file
  integer(kind=4) :: record_dia_his = 0
  real(kind=8),save :: output_time_his = 0
  real(kind=8),save :: output_time_avg = 0
  real(kind=8),save :: output_time_dia_his = 0
  real(kind=8),save :: output_time_dia_avg = 0

  public bgc_idx, t_idx
  public set_bgc_surf_frc
  public init_arrays_bgc_frc

  public wrt_bgc_diags
  public wrt_bgc

contains

! ----------------------------------------------------------------------
  subroutine set_bgc_surf_frc(istr,iend,jstr,jend)  ![

    ! read in bgc surface flux and interpolate to model time.
    ! Taken from get_smth.F & set_smth.F of ETH code.

    ! Since bgc is set up such that iron and dust have flux at surface but
    ! are not directly bgc tracers, I have separated them into set_bgc_frc

    ! this routine is called by # if defined BIOLOGY_BEC2 in set_forces.F of ETH code

    use scalars, only: cp, nrhs, rho0
    use dimensions, only: inode, jnode
    use bgc_forces, only: &
#ifdef PCO2AIR_FORCING
    & xco2air, &
#ifdef MARBL
    & xco2air_alt, &
#endif
#endif
#ifdef NOX_FORCING
    & nox, &
#endif
#ifdef NHY_FORCING
    & nhy,&
#endif
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
    & swrad_avg, &
#endif
    & dust, iron
    use dimensions, only: npx, npy
    implicit none

! Needed for iron and dust variables, not yet ported to this module. Should do.
!#include "bgc_forces.h"

    ! input/outputs
    integer(kind=4),intent(in) :: istr,iend,jstr,jend

    ! local
    integer(kind=4) :: i, j

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
# include "compute_extended_bounds.h"
#endif

#ifdef PARALLEL_IO
    pio_file_is_open = 0
#endif
    call set_frc_data(nc_dust,dust,'r')
    call set_frc_data(nc_iron,iron,'r')


#ifdef PCO2AIR_FORCING
    call set_frc_data(nc_xco2air,xco2air,'r')
#ifdef MARBL
    call set_frc_data(nc_xco2air_alt,xco2air_alt,'r')
#endif
#endif /* PCO2AIR_FORCING */

#ifdef NOX_FORCING
    call set_frc_data(nc_nox,nox,'r')
#endif /* NOX_FORCING */
#ifdef NHY_FORCING
    call set_frc_data(nc_nhy,nhy,'r')
#endif /* NHY_FORCING */

    call error_log%abort_check()
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
    ! Daily average Short wave radiation
    call set_frc_data(nc_swrad_avg,swrad_avg,'r')
    call error_log%abort_check()
    do j=jstrR,jendR
      do i=istrR,iendR
        swrad_avg(i,j)= swrad_avg(i,j)/(rho0*Cp)

#ifdef SEA_ICE_NOFLUX
        if( t(i,j,nz,nrhs,itemp) .le. -1.8_8 ) then    ! Restrict stflx to prevent surface temperature to go below -2 d
#   if defined LMD_KPP
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
          swrad_avg(i,j)=0._8
#endif
#   endif
        endif
#endif
      enddo !<- istrR
    enddo   !<- jstrR
#endif /* DAILYPAR_PHOTOINHIBITION */

#ifdef PARALLEL_IO
    if (pio_file_is_open == 1) then
      call PIO_closefile(pio_FileDesc)
    endif
    pio_file_is_open = 0
#endif

  end subroutine set_bgc_surf_frc  !]
! ----------------------------------------------------------------------
  subroutine wrt_bgc_tracers(avg)  ![
    ! write variables to output netcdf file
    ! bgc variables are in fact tracers, but we choose to output
    ! them to their own file.
    ! this file also include bgc diagnostic variables.
    ! selection of variables to output is still done in tracers.opt
    ! otherwise would cause circular reference as needed in both directions.
    !
    ! TIME: since bgc is calculated in step3d_t_ISO.F and uses 'nnew'
    ! for the tracers (i.e. t=n+1). The bgc tracers are for t=n+1._8
    ! since we now write them out immediately, we need to set the time to
    ! t=n+1, i.e. time+dt.
    implicit none

    ! import/export
    logical, intent(in) :: avg

! local
    character(len=15) :: sr_name = "wrt_bgc_tracers"
    character(len=max_name_size),save :: fname_his,fname_avg
    character(len=1024) :: error_info
    integer(kind=4)                :: itrc, itavg, ierr
    logical,save           :: first_step_avg=.true.
    logical,save           :: first_step_his=.true.

    if (wrt_bgc_his .and. mod(output_period_bgc_his, dt) /= 0) then
      write(error_info, *)&
      &"dt (", dt, ") is not a factor of ",&
      &"output_period_bgc_his (",&
      &output_period_bgc_his, ")"

      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif

    if (wrt_bgc_avg .and. mod(output_period_bgc_avg,dt) /= 0) then
      write(error_info, *)&
      &"dt (", dt, ") is not a factor of ",&
      &"output_period_bgc_avg (",&
      &output_period_bgc_avg, ")"

      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif

    if (first_step_his) then
      record_his = nrpf_bgc_his
      record_dia_his = nrpf_bgc_his_dia
    end if
    if (first_step_avg) then
      record_avg = nrpf_bgc_avg
      record_dia_avg = nrpf_bgc_avg_dia
    end if

    call error_log%abort_check()

#ifdef PARALLEL_IO
    if (avg) then

!        if (.not. first_step_avg) then
      call calc_avg_bgc                          ! don't include t=0 in averaging
      output_time_avg = output_time_avg + dt     ! only start count after first timestep
!        endif
!        first_step_avg=.false.


!        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,F11.4,1x,A,I4)')
!     &       'wrt_bgc_tracers :: average, time =', output_time_avg,
!     &       'period =', output_period_bgc_avg, 'rec =', record_avg
!        endif

      if (output_time_avg>=output_period_bgc_avg) then  ! time for an output
        if (mod(record_avg,nrpf_bgc_avg)==0) then
          if (mynode == 0) then
            call create_bgc_file(fname_avg,avg)
          endif
          call MPI_Bcast(fname_avg,max_name_size,MPI_CHARACTER,0,ocean_grid_comm,ierr)
          call MPI_Barrier(ocean_grid_comm, ierr)
          record_avg = 0
        endif
        record_avg = record_avg + 1
        output_time_avg = 0
        navg_bgc = 0
        if (mynode == 0) then
          ierr=nf90_open(fname_avg,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_avg/))
          ierr=nf90_close(ncid)
        endif

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname_avg), PIO_write)

        pio_gtype = '3Drw'
        do itrc=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1, nt
          if (wrt_t_avg(itrc)) then
            if (mynode == 0) then
              write(*,'(7x,A,1x,A)')&
              &'wrt_bgc_tracers :: average, trc = ', t_vname(itrc)
            endif
            itavg = NT_2_t_avg(itrc)                         ! get respective index for t_avg(itavg) -> t(itrc)
            call ncwrite(ncid,t_vname(itrc),t_avg(i0:i1,j0:j1,:,itavg),(/1,1,1,record_avg/),.true.)
          endif
        enddo
        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc :: wrote average, tdays =', tdays,&
          &'step =', iic, 'rec =', record_avg
        endif
        ! close inside the output block: the file is only open when writing
        call PIO_closefile(pio_FileDesc)
      endif

    else

      if (.not. first_step_his) then
        output_time_his = output_time_his + dt
      endif
      first_step_his=.false.

!        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,F11.4,1x,A,I4)')
!     &       'wrt_bgc_tracers :: history, time =', output_time_his,
!     &       'period =', output_period_bgc_his, 'rec =', record_his
!        endif

      if (output_time_his>=output_period_bgc_his .or.&
      &output_time_his==0                 ) then  ! time for an output
        if (mod(record_his,nrpf_bgc_his)==0) then
          if (mynode == 0) then
            call create_bgc_file(fname_his,avg)
          endif
          call MPI_Bcast(fname_his,max_name_size,MPI_CHARACTER,0,ocean_grid_comm,ierr)
          call MPI_Barrier(ocean_grid_comm, ierr)
          record_his = 0
        endif
        record_his = record_his + 1
        output_time_his=0
        if (mynode == 0) then
          ierr=nf90_open(fname_his,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_his/))
          ierr=nf90_close(ncid)
        endif
        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname_his), PIO_write)
        pio_gtype = '3Drw'
        do itrc=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1, nt
          if (wrt_t(itrc)) then
            if (mynode == 0) then
              write(*,'(7x,A,1x,A)')&
              &'wrt_bgc_tracers :: history , trc = ', t_vname(itrc)
            endif
!            call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nstp,itrc),(/1,1,1,record_his/))
            call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record_his/),.true.)
          endif
        enddo
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc :: wrote history, tdays =', tdays,&
          &'step =', iic-1, 'rec =', record_his
        endif
        ! close inside the output block: the file is only open when writing
        call PIO_closefile(pio_FileDesc)
      endif  ! <-- wrt_file_his

    endif

#else ! PARALLEL_IO

    if (avg) then

!        if (.not. first_step_avg) then
      call calc_avg_bgc                          ! don't include t=0 in averaging
      output_time_avg = output_time_avg + dt     ! only start count after first timestep
!        endif
!        first_step_avg=.false.


!        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,F11.4,1x,A,I4)')
!     &       'wrt_bgc_tracers :: average, time =', output_time_avg,
!     &       'period =', output_period_bgc_avg, 'rec =', record_avg
!        endif

      if (output_time_avg>=output_period_bgc_avg) then  ! time for an output
        if (mod(record_avg,nrpf_bgc_avg)==0) then
          call create_bgc_file(fname_avg,avg)
          record_avg = 0
        endif
        record_avg = record_avg + 1
        output_time_avg = 0
        navg_bgc = 0
        ierr=nf90_open(fname_avg,nf90_write,ncid)
        call ncwrite(ncid,'ocean_time',(/time/),(/record_avg/))
        do itrc=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1, nt
          if (wrt_t_avg(itrc)) then
            if (mynode == 0) then
              write(*,'(7x,A,1x,A)')&
              &'wrt_bgc_tracers :: average, trc = ', t_vname(itrc)
            endif
            itavg = NT_2_t_avg(itrc)                         ! get respective index for t_avg(itavg) -> t(itrc)
            call ncwrite(ncid,t_vname(itrc),t_avg(i0:i1,j0:j1,:,itavg),(/1,1,1,record_avg/))
          endif
        enddo
        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc :: wrote average, tdays =', tdays,&
          &'step =', iic, 'rec =', record_avg
        endif
      endif

    else

      if (.not. first_step_his) then
        output_time_his = output_time_his + dt
      endif
      first_step_his=.false.

!        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,F11.4,1x,A,I4)')
!     &       'wrt_bgc_tracers :: history, time =', output_time_his,
!     &       'period =', output_period_bgc_his, 'rec =', record_his
!        endif

      if (output_time_his>=output_period_bgc_his .or.&
      &output_time_his==0                 ) then  ! time for an output
        if (mod(record_his,nrpf_bgc_his)==0) then
          call create_bgc_file(fname_his,avg)
          record_his = 0
        endif
        record_his = record_his + 1
        output_time_his=0
        ierr=nf90_open(fname_his,nf90_write,ncid)
        call ncwrite(ncid,'ocean_time',(/time/),(/record_his/))
        do itrc=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1, nt
          if (wrt_t(itrc)) then
            if (mynode == 0) then
              write(*,'(7x,A,1x,A)')&
              &'wrt_bgc_tracers :: history , trc = ', t_vname(itrc)
            endif
!            call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nstp,itrc),(/1,1,1,record_his/))
            call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record_his/))
          endif
        enddo
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc :: wrote history, tdays =', tdays,&
          &'step =', iic-1, 'rec =', record_his
        endif
      endif

    endif

    ierr=nf90_close(ncid)
#endif
  end subroutine wrt_bgc_tracers  !]
! ----------------------------------------------------------------------
  subroutine wrt_bgc_diags ![
    implicit none
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
! local
    character(len=13) :: sr_name = "wrt_bgc_diags"
    character(len=max_name_size),save :: fname_his,fname_avg
    character(len=1024) :: error_info
    integer(kind=4)                :: itrc, ierr, idiag
    logical,save           :: first_step_dia_avg=.true.
    logical,save           :: first_step_dia_his=.true.

    if (wrt_bgc_dia_his .and. mod(output_period_bgc_his_dia, dt) /= 0) then
      write(error_info, *)&
      &"dt (", dt, ") is not a factor of ",&
      &"output_period_bgc_his_dia (",&
      &output_period_bgc_his_dia, ")"

      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif

    if (wrt_bgc_dia_avg .and. mod(output_period_bgc_avg_dia,dt) /= 0) then
      write(error_info, *)&
      &"dt (", dt, ") is not a factor of ",&
      &"output_period_bgc_avg_dia (",&
      &output_period_bgc_avg_dia, ")"

      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif
    call error_log%abort_check()

    if (wrt_bgc_dia_avg) then

#ifdef PARALLEL_IO
!         if (.not. first_step_dia_avg) then
      call calc_avg_dia_bgc                              ! don't include t=0 in averaging
      output_time_dia_avg = output_time_dia_avg + dt     ! only start count after first timestep
!         endif
!         first_step_dia_avg=.false.

      if (output_time_dia_avg>=output_period_bgc_avg_dia) then  ! time for an output
        if (mod(record_dia_avg,nrpf_bgc_avg_dia)==0) then
          if (mynode == 0) then
            call create_bgc_dia_file(fname_avg,.true.)
          endif
          call MPI_Bcast(fname_avg,max_name_size,MPI_CHARACTER,0,ocean_grid_comm,ierr)
          call MPI_Barrier(ocean_grid_comm, ierr)
          record_dia_avg = 0
        endif
        record_dia_avg = record_dia_avg + 1
        output_time_dia_avg = 0
        navg_dia_bgc = 0
        if (mynode == 0) then
          ierr=nf90_open(fname_avg,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_dia_avg/))
          ierr=nf90_close(ncid)
        endif

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname_avg), PIO_write)

        pio_gtype = '2Drw'
        do itrc=1,nr_bgc_diag_2d
          if (wrt_bgc_diag_2d(itrc)) then
            idiag = idx_bgc_diag_2d(itrc)
            call ncwrite(ncid,vname_bgc_diag_2d(1,itrc),&
            &bgc_diag_2d_avg(i0:i1,j0:j1,idiag),(/1,1,record_dia_avg/),.true.)
          endif
        enddo
        pio_gtype = '3Drw'
        do itrc=1,nr_bgc_diag_3d
          if (wrt_bgc_diag_3d(itrc)) then
            idiag = idx_bgc_diag_3d(itrc)
            call ncwrite(ncid,vname_bgc_diag_3d(1,itrc),&
            &bgc_diag_3d_avg(i0:i1,j0:j1,:,idiag),(/1,1,1,record_dia_avg/),.true.)
          endif
        enddo

        call PIO_closefile(pio_FileDesc)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc diag :: wrote average, tdays =', tdays,&
          &'step =', iic, 'rec =', record_dia_avg
        endif
      endif
#else ! PARALLEL_IO
!         if (.not. first_step_dia_avg) then
      call calc_avg_dia_bgc                              ! don't include t=0 in averaging
      output_time_dia_avg = output_time_dia_avg + dt     ! only start count after first timestep
!         endif
!         first_step_dia_avg=.false.

      if (output_time_dia_avg>=output_period_bgc_avg_dia) then  ! time for an output
        if (mod(record_dia_avg,nrpf_bgc_avg_dia)==0) then
          call create_bgc_dia_file(fname_avg,.true.)
          record_dia_avg = 0
        endif
        record_dia_avg = record_dia_avg + 1
        output_time_dia_avg = 0
        navg_dia_bgc = 0
        ierr=nf90_open(fname_avg,nf90_write,ncid)
        call ncwrite(ncid,'ocean_time',(/time/),(/record_dia_avg/))
        do itrc=1,nr_bgc_diag_2d
          if (wrt_bgc_diag_2d(itrc)) then
            idiag = idx_bgc_diag_2d(itrc)
            call ncwrite(ncid,vname_bgc_diag_2d(1,itrc),&
            &bgc_diag_2d_avg(i0:i1,j0:j1,idiag),(/1,1,record_dia_avg/))
          endif
        enddo
        do itrc=1,nr_bgc_diag_3d
          if (wrt_bgc_diag_3d(itrc)) then
            idiag = idx_bgc_diag_3d(itrc)
            call ncwrite(ncid,vname_bgc_diag_3d(1,itrc),&
            &bgc_diag_3d_avg(i0:i1,j0:j1,:,idiag),(/1,1,1,record_dia_avg/))
          endif
        enddo
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc diag :: wrote average, tdays =', tdays,&
          &'step =', iic, 'rec =', record_dia_avg
        endif
      endif
      ierr=nf90_close(ncid)
#endif ! PARALLEL_IO
    endif

    if (wrt_bgc_dia_his) then

#ifdef PARALLEL_IO
      if (.not. first_step_dia_his) then
        output_time_dia_his = output_time_dia_his + dt
      endif
      first_step_dia_his=.false.

      if (output_time_dia_his>=output_period_bgc_his_dia .or.&
      &output_time_dia_his==0                     ) then  ! time for an output
        if (mod(record_dia_his,nrpf_bgc_his_dia)==0) then
          if (mynode == 0) then
            call create_bgc_dia_file(fname_his,.false.)
          endif
          call MPI_Bcast(fname_his,max_name_size,MPI_CHARACTER,0,ocean_grid_comm,ierr)
          call MPI_Barrier(ocean_grid_comm, ierr)
          record_dia_his = 0
        endif
        record_dia_his = record_dia_his + 1
        output_time_dia_his = 0
        if (mynode == 0) then
          ierr=nf90_open(fname_his,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_dia_his/))
          ierr=nf90_close(ncid)
        endif

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname_his), PIO_write)

        pio_gtype = '2Drw'
        do itrc=1,nr_bgc_diag_2d
          if (wrt_bgc_diag_2d(itrc)) then
            idiag = idx_bgc_diag_2d(itrc)
            call ncwrite(ncid,vname_bgc_diag_2d(1,itrc),&
            &bgc_diag_2d(i0:i1,j0:j1,idiag),(/1,1,record_dia_his/),.true.)
          endif
        enddo
        pio_gtype = '3Drw'
        do itrc=1,nr_bgc_diag_3d
          if (wrt_bgc_diag_3d(itrc)) then

            idiag = idx_bgc_diag_3d(itrc)
            call ncwrite(ncid,vname_bgc_diag_3d(1,itrc),&
            &bgc_diag_3d(i0:i1,j0:j1,:,idiag),(/1,1,1,record_dia_his/),.true.)
          endif
        enddo

        call PIO_closefile(pio_FileDesc)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc diag :: wrote history, tdays =', tdays,&
          &'step =', iic, 'rec =', record_dia_his
        endif
      endif  ! <-- wrt_file_his

#else ! PARALLEL_IO

      if (.not. first_step_dia_his) then
        output_time_dia_his = output_time_dia_his + dt
      endif
      first_step_dia_his=.false.

      if (output_time_dia_his>=output_period_bgc_his_dia .or.&
      &output_time_dia_his==0                     ) then  ! time for an output
        if (mod(record_dia_his,nrpf_bgc_his_dia)==0) then
          call create_bgc_dia_file(fname_his,.false.)
          record_dia_his = 0
        endif
        record_dia_his = record_dia_his + 1
        output_time_dia_his = 0
        ierr=nf90_open(fname_his,nf90_write,ncid)
        call ncwrite(ncid,'ocean_time',(/time/),(/record_dia_his/))
        do itrc=1,nr_bgc_diag_2d
          if (wrt_bgc_diag_2d(itrc)) then
            idiag = idx_bgc_diag_2d(itrc)
            call ncwrite(ncid,vname_bgc_diag_2d(1,itrc),&
            &bgc_diag_2d(i0:i1,j0:j1,idiag),(/1,1,record_dia_his/))
          endif
        enddo
        do itrc=1,nr_bgc_diag_3d
          if (wrt_bgc_diag_3d(itrc)) then

            idiag = idx_bgc_diag_3d(itrc)
            call ncwrite(ncid,vname_bgc_diag_3d(1,itrc),&
            &bgc_diag_3d(i0:i1,j0:j1,:,idiag),(/1,1,1,record_dia_his/))
          endif
        enddo
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
          &'bgc diag :: wrote history, tdays =', tdays,&
          &'step =', iic, 'rec =', record_dia_his
        endif

        ierr=nf90_close(ncid)

      endif  ! <-- wrt_file_his

#endif ! PARALLEL_IO

    endif

#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */

  end subroutine wrt_bgc_diags  !]
!----------------------------------------------------------------------
  subroutine create_bgc_file(fname,avg)  ![
    implicit none

    !input/output
    character(len=max_name_size),intent(out) :: fname
    logical,          intent(in)  :: avg                 ! his or average file

    ! local
    integer(kind=4) :: ncid,ierr,varid

#ifdef PARALLEL_IO
    if (avg) then
      call create_file('_bgc_avg',fname, nonode=.true.)
    else
      call create_file('_bgc',fname, nonode=.true.)
    endif
#else
    if (avg) then
      call create_file('_bgc_avg',fname)
    else
      call create_file('_bgc',fname)
    endif
#endif
    ierr=nf90_open(fname,nf90_write,ncid)
    ierr=nf90_redef(ncid)

    call def_vars_bgc( ncid,avg )

    ierr = nf90_enddef(ncid)
  end subroutine create_bgc_file !]
!----------------------------------------------------------------------
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
  subroutine create_bgc_dia_file(fname,avg)  ![
    implicit none

    !input/output
    character(len=max_name_size),intent(out) :: fname
    logical,          intent(in)  :: avg                 ! his or average file

    ! local
    integer(kind=4) :: ncid,ierr,varid

#ifdef PARALLEL_IO
    if (avg) then
      call create_file('_bgc_dia_avg',fname,nonode=.true.)
    else
      call create_file('_bgc_dia',fname,nonode=.true.)
    endif
#else
    if (avg) then
      call create_file('_bgc_dia_avg',fname)
    else
      call create_file('_bgc_dia',fname)
    endif
#endif
    ierr=nf90_open(fname,nf90_write,ncid)
    ierr=nf90_redef(ncid)

    call def_bgc_diag(ncid,avg)

    ierr = nf90_enddef(ncid)
  end subroutine create_bgc_dia_file !]
#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */
! ----------------------------------------------------------------------
  subroutine def_vars_bgc( ncid, avg )  ![
    ! define output variable & attributes in netcdf results file
    ! for actual or averaged variables
    implicit none

    ! input
    integer(kind=4),intent(in) :: ncid
    logical,intent(in) :: avg
    ! local
    integer(kind=4)           :: ierr, varid, itrc
    character(len=64) :: text_lname

    if (avg) then
      do itrc=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1,NT
        if (wrt_t_avg(itrc)) then
          text_lname='avg_'//t_lname(itrc)
          varid = nccreate(ncid,t_vname(itrc),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/), nf90_double)

          ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
          ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
        endif
      enddo
    else
      do itrc=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1,NT
        if (wrt_t(itrc)) then
          text_lname=t_lname(itrc)
          varid = nccreate(ncid,t_vname(itrc),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/), nf90_double)
          ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
          ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
        endif
      enddo
    endif

  end subroutine def_vars_bgc  !]
! ----------------------------------------------------------------------
  subroutine calc_avg_bgc  ![
    implicit none

    ! local
    real(kind=8)    :: coef
    integer(kind=4) :: itrc, itavg, k

    navg_bgc = navg_bgc +1
    coef = 1._8/navg_bgc

    if (coef==1) then                                    ! this refreshes average (1-coef)=0
      if (mynode==0) write(*,'(7x,2A,F9.1)')&
      &'bgc :: started averaging. ',&
      &'output_period_bgc_avg (s) =', output_period_bgc_avg
    endif

    t_avg_bgc = t_avg_bgc*(1-coef) + time*coef

    ! need i0:i1 indices because arrays still GLOBAL_2D therefore wasted margin
    do itrc=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1, NT
      if (wrt_t_avg(itrc)) then
        itavg = NT_2_t_avg(itrc)                         ! get respective index for t_avg(itavg) -> t(itrc)
        t_avg(i0:i1,j0:j1,:,itavg) = t_avg(i0:i1,j0:j1,:,itavg)    *(1-coef)&
        &+     t(i0:i1,j0:j1,:,nstp,itrc)*coef     ! CONFIRM NSTP OR NNEW!!!
      endif
    enddo

  end subroutine calc_avg_bgc  !]
! ----------------------------------------------------------------------
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
  subroutine calc_avg_dia_bgc  ![
    implicit none

    ! local
    real(kind=8)    :: coef
    integer(kind=4) :: itrc, itavg, k

    navg_dia_bgc = navg_dia_bgc +1
    coef = 1._8/navg_dia_bgc

    if (coef==1) then                                    ! this refreshes average (1-coef)=0
      if (mynode==0) write(*,'(7x,2A,F9.1)')&
      &'bgc diag :: started averaging. ',&
      &'output_period_bgc_avg_dia (s) =', output_period_bgc_avg_dia
    endif

    t_avg_dia_bgc = t_avg_dia_bgc*(1-coef) + time*coef

    do itrc=1,nr_bgc_wrdiag_2d
      bgc_diag_2d_avg(i0:i1,j0:j1,itrc) = bgc_diag_2d_avg(i0:i1,j0:j1,itrc) *(1-coef)&
      &+ bgc_diag_2d(i0:i1,j0:j1,itrc)     *coef
    end do
    do itrc=1,nr_bgc_wrdiag_3d
      do k=1,nz
        bgc_diag_3d_avg(i0:i1,j0:j1,k,itrc) = bgc_diag_3d_avg(i0:i1,j0:j1,k,itrc) *(1-coef)&
        &+ bgc_diag_3d(i0:i1,j0:j1,k,itrc)     *coef
      end do
    end do

  end subroutine calc_avg_dia_bgc  !]
#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */
! ----------------------------------------------------------------------
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
  subroutine def_bgc_diag(ncid,avg)  ![

    ! Define history/avg file variables for BEC2_DIAG
    ! Taken from def_his.F of ETH code.
    ! Called from def_his.F

    implicit none

    ! input
    integer(kind=4),intent(in) :: ncid
    logical,intent(in) :: avg
    ! local
    integer(kind=4) :: idiag, ierr, varid

    ! 2d diagnostics:
    do idiag=1,nr_bgc_diag_2d
      if (wrt_bgc_diag_2d(idiag)) then
        varid = nccreate(ncid,vname_bgc_diag_2d(1,idiag),&
        &(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/), nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',vname_bgc_diag_2d(2,idiag))
        ierr = nf90_put_att(ncid,varid,'units',vname_bgc_diag_2d(3,idiag))
      endif
    enddo

    ! 3d diagnostics:
    do idiag=1,nr_bgc_diag_3d
      if (wrt_bgc_diag_3d(idiag)) then
        varid = nccreate(ncid,vname_bgc_diag_3d(1,idiag),&
        &(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,nz,0/), nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',vname_bgc_diag_3d(2,idiag))
        ierr = nf90_put_att(ncid,varid,'units',vname_bgc_diag_3d(3,idiag))
      endif
    enddo

  end subroutine def_bgc_diag  !]
#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */
! ----------------------------------------------------------------------
  subroutine init_arrays_bgc_frc  ![
    implicit none

    allocate( nc_dust%vdata( GLOBAL_2D_ARRAY,2) )
    if (interp_bgc_frc) nc_dust%coarse = 1
    allocate( nc_iron%vdata( GLOBAL_2D_ARRAY,2) )
    if (interp_bgc_frc) nc_iron%coarse = 1

#ifdef PCO2AIR_FORCING
    allocate( nc_xco2air%vdata( GLOBAL_2D_ARRAY,2) )
    if (interp_bgc_frc) nc_xco2air%coarse = 1
#ifdef MARBL
    allocate( nc_xco2air_alt%vdata( GLOBAL_2D_ARRAY,2) )
    if (interp_bgc_frc) nc_xco2air_alt%coarse = 1
#endif
#endif

#ifdef NOX_FORCING
    allocate( nc_nox%vdata( GLOBAL_2D_ARRAY,2) )
    if (interp_bgc_frc) nc_nox%coarse = 1
#endif
#ifdef NHY_FORCING
    allocate( nc_nhy%vdata( GLOBAL_2D_ARRAY,2) )
    if (interp_bgc_frc) nc_nhy%coarse = 1
#endif

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
    allocate( nc_swrad_avg%vdata( GLOBAL_2D_ARRAY,2) )
    if (interp_bgc_frc) nc_swrad_avg%coarse = 1
#endif

  end subroutine init_arrays_bgc_frc !]
! ----------------------------------------------------------------------
  subroutine wrt_bgc(init)  ![
#ifdef MARBL_DIAGS
    use dimensions, only: inode, jnode
    use param, only: nsub_e, nsub_x
#endif
    implicit none

!     import/export
    logical,optional,intent(in)       :: init
    integer(kind=4) :: tile = 0
#ifdef MARBL_DIAGS
#include "compute_tile_bounds.h"
#endif
    ! initialisation writting
    if (wrt_bgc_his.and.present(init)) then
      call wrt_bgc_tracers(.false.) ! write his data
    else
      ! local
      if (wrt_bgc_his)     call wrt_bgc_tracers(.false.) ! write his data
      if (wrt_bgc_avg)     call wrt_bgc_tracers(.true.) ! write avg data
    endif

#if defined BEC2_DIAG || defined MARBL_DIAGS
#if defined MARBL
    if (present(init)) then
      call marbldrv_compute_init_diagnostics(istr,iend,jstr,jend,t)
    end if
#endif
    call wrt_bgc_diags
#endif


  end subroutine wrt_bgc  !]
! ----------------------------------------------------------------------
#endif /*(BIOLOGY_BEC2 || MARBL)*/

end module bgc_io
