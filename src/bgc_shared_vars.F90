module bgc_shared_vars
  ! formerly bgc_ecosys_bec2.h

! # TODO marbldrve_configure_diagnostics for pointers example

#include "cppdefs.opt"
#if defined(BIOLOGY_BEC2) || defined(MARBL)
  use namelist_open_mod, only: open_namelist_file
  use error_handling_mod, only: error_log
  use param, only: mynode, nt_passive, nt_bgc, lm, mm, nz
  use roms_read_write, only: ncforce
  use bgc_forces, only: xco2air_default   ! declared in bgc_forces; read here via BGC_SETTINGS
  use tracers, only: itands, t_vname, wrt_t, t_lname, t
#if defined(BIOLOGY_BEC2)
  use tracers, only: ialk, idic, ife, io2, ipo4, isio3, idoc,&
       &idofe, idon, idonr, idop, idopr, inh4, ino3, izooc, idiatc,&
       &idiatchl, idiatfe, idiatsi, idiazc, idiazchl, idiazfe, ispc,&
       &ispcaco3, ispchl, ispfe
#endif


#if defined(BIOLOGY_BEC2)
  use bec2_vars, only:&
#if defined(BEC2_DIAG)
       &nr_bec2_diag_2d, nr_bec2_diag_3d,&
       &vname_bec2_diag_2d, wrt_bec2_diag_2d, vname_bec2_diag_3d,&
       &wrt_bec2_diag_3d, bec2_diag_2d, bec2_diag_3d,&
       &ammox_idx_t, anammox_idx_t, atmpress_idx_t, caco3fluxin_idx_t,&
       &caco3fluxtosed_idx_t, caco3prod_idx_t, caco3remin_idx_t,&
       &caco3sedloss_idx_t, co2star_idx_t, dco2star_idx_t,&
       &denitrif1_idx_t, denitrif2_idx_t, denitrif3_idx_t,&
       &denitrif_idx_t, diatagg_idx_t, diatczero_idx_t,&
       &diatfeuptake_idx_t, diatlightlim_idx_t, diatloss_idx_t,&
       &diatnh4uptake_idx_t, diatnlim_idx_t, diatno2uptake_idx_t,&
       &diatno3uptake_idx_t, diatphotoacc_idx_t, diatplim_idx_t,&
       &diatpo4uptake_idx_t, diatsio3uptake_idx_t, diatsiuptake_idx_t,&
       &diazagg_idx_t, diazczero_idx_t, diazfeuptake_idx_t,&
       &diazlightlim_idx_t, diazloss_idx_t, diaznfix_idx_t,&
       &diaznh4uptake_idx_t, diazno2uptake_idx_t, diazno3uptake_idx_t,&
       &diazphotoacc_idx_t, diazplim_idx_t, diazpo4uptake_idx_t,&
       &docprod_idx_t, docremin_idx_t, doczero_idx_t, dofeprod_idx_t,&
       &doferemin_idx_t, donprod_idx_t, donremin_idx_t, donrremin_idx_t,&
       &dopprod_idx_t, dopremin_idx_t, dustfluxin_idx_t,&
       &dustfluxtosed_idx_t, dustremin_idx_t, fescavenge_idx_t,&
       &fescavengerate_idx_t, fesedflux_idx_t, fgco2_idx_t, fgn2_idx_t,&
       &fgn2o_idx_t, fgo2_idx_t, fluxtosed_idx_t, grazediat_idx_t,&
       &grazediatzoo_idx_t, grazediaz_idx_t, grazediazzoo_idx_t,&
       &grazedicdiat_idx_t, grazedicdiaz_idx_t, grazedicsp_idx_t,&
       &grazesp_idx_t, grazespzoo_idx_t, idx_bec2_diag_2d,&
       &idx_bec2_diag_3d, ironflux_idx_t, ironuptakediat_idx_t,&
       &ironuptakediaz_idx_t, ironuptakesp_idx_t, lossdicdiat_idx_t,&
       &lossdicdiaz_idx_t, lossdicsp_idx_t, n2oammox_idx_t,&
       &n2osat_idx_t, n2sat_idx_t, nitrif_idx_t, nitrox_idx_t,&
       &o2cons_idx_t, o2prod_idx_t, o2sat_idx_t, otherremin_idx_t,&
       &par_idx_t, parinc_idx_t, pcaco3prod_idx_t, xco2_idx_t,&
       &xco2air_idx_t, xco2oc_idx_t, ph_idx_t, photocdiat_idx_t,&
       &photocdiaz_idx_t, photocsp_idx_t, pironfluxin_idx_t,&
       &pironfluxtosed_idx_t, pironprod_idx_t, pironremin_idx_t,&
       &pocfluxin_idx_t, pocprod_idx_t, pocremin_idx_t,&
       &pocsedloss_idx_t, pvco2_idx_t, pvn2_idx_t, pvn2o_idx_t,&
       &pvo2_idx_t, schmidt_n2_idx_t, schmidt_n2o_idx_t,&
       &schmidtco2_idx_t, schmidto2_idx_t, seddenitrif_idx_t,&
       &sio2fluxin_idx_t, sio2fluxtosed_idx_t, sio2prod_idx_t,&
       &sio2remin_idx_t, sio2sedloss_idx_t, spagg_idx_t,&
       &spcaco3zero_idx_t, spczero_idx_t, spfeuptake_idx_t,&
       &splightlim_idx_t, sploss_idx_t, spnh4uptake_idx_t, spnlim_idx_t,&
       &spno2uptake_idx_t, spno3uptake_idx_t, spphotoacc_idx_t,&
       &spplim_idx_t, sppo4uptake_idx_t, spqcaco3_idx_t, totchl_idx_t,&
       &totphytoc_idx_t, totprod_idx_t, ws10m_idx_t, xkw_idx_t,&
       &zooczero_idx_t, zooloss_idx_t, zoolossdic_idx_t, &
#endif
       &init_arrays_bec2_vars
#elif defined(MARBL)
  use marbl_driver, only: &
# ifdef MARBL_DIAGS
       & marbldrv_configure_diagnostics&
       &,nr_marbl_diag_2d,nr_marbl_diag_3d&
       &,marbl_diag_2d,marbl_diag_3d&
       &,idx_marbl_diag_2d,idx_marbl_diag_3d&
       &,wrt_marbl_diag_2d,wrt_marbl_diag_3d, &
# endif
       marbldrv_configure_saved_state
#endif

  implicit none
  character(len=18) :: module_name = "bgc_shared_vars.F"
! BGC forcing file settings
  type (ncforce) :: nc_dust = ncforce(&
  &vname='dust', tname='dust_time' ) ! dust forcing
  type (ncforce) :: nc_iron = ncforce(&
  &vname='iron', tname='iron_time' ) ! iron forcing
#ifdef PCO2AIR_FORCING
  type (ncforce) :: nc_xco2air = ncforce(&
  &vname='xco2_air',   tname='xco2_time' )
#ifdef MARBL
  type (ncforce) :: nc_xco2air_alt = ncforce(&
  &vname='xco2_air_alt',   tname='xco2_time' )
#endif
#endif
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
  type (ncforce) :: nc_swrad_avg = ncforce(&
  &vname='swrad_davg', tname='srf_davg_time' )
#endif
#ifdef NOX_FORCING
  type (ncforce) :: nc_nox = ncforce(&
  &vname='nox', tname='nox_time' )
#endif
#ifdef NHY_FORCING
  type (ncforce) :: nc_nhy = ncforce(&
  &vname='nhy', tname='nhy_time' )
#endif

  ! note: choice of bgc tracers to output is still selected in tracers.opt.

  real(kind=8)    :: output_period_bgc_his = 0._8        ! output period in seconds
  integer(kind=4) :: nrpf_bgc_his          = 0                  ! total recs per file

  real(kind=8)    :: output_period_bgc_avg = 0._8        ! output averaging period in seconds
  integer(kind=4) :: nrpf_bgc_avg          = 0                  ! total recs per file

  real(kind=8)    :: output_period_bgc_his_dia = 0._8  ! output period in seconds
  integer(kind=4) :: nrpf_bgc_his_dia          = 0             ! total recs per file

  real(kind=8)    :: output_period_bgc_avg_dia = 0._8    ! output period in seconds
  integer(kind=4) :: nrpf_bgc_avg_dia          = 0              ! total recs per file

  integer(kind=4) :: nbgc_flx = 0                   ! number of surface bgc flux forcings

  logical, public :: wrt_bgc_his, wrt_bgc_avg,&
  &wrt_bgc_dia_his, wrt_bgc_dia_avg,&
  &interp_bgc_frc

  namelist /BGC_SETTINGS/&
  &wrt_bgc_his, wrt_bgc_avg, wrt_bgc_dia_his, wrt_bgc_dia_avg,&
  &interp_bgc_frc,&
  &output_period_bgc_his, nrpf_bgc_his,&
  &output_period_bgc_avg, nrpf_bgc_avg,&
  &output_period_bgc_his_dia, nrpf_bgc_his_dia,&
  &output_period_bgc_avg_dia, nrpf_bgc_avg_dia,&
  &nbgc_flx,&
  &xco2air_default
#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
!
! Diagnostic variables appearing in average and history files:
!
  integer(kind=4) :: nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d, nr_bgc_wrdiag
#if defined(BEC2_DIAG)
  integer(kind=4), parameter :: nr_bgc_diag_2d = nr_bec2_diag_2d
  integer(kind=4), parameter :: nr_bgc_diag_3d = nr_bec2_diag_3d
#elif defined(MARBL_DIAGS)
  ! Just hardcoding this for now but would be nice to obtain it dynamically
  integer(kind=4), parameter ::  nr_bgc_diag_3d = nr_marbl_diag_3d
  integer(kind=4), parameter ::  nr_bgc_diag_2d = nr_marbl_diag_2d
#endif
  integer(kind=4), parameter :: nr_bgc_diag = nr_bgc_diag_2d + nr_bgc_diag_3d

  real(kind=8),allocatable,dimension(:,:,:,:),target   :: bgc_diag_3d
  real(kind=8),allocatable,dimension(:,:,:)  ,target   :: bgc_diag_2d

! Control BEC2_DIAG Output vars
  logical, target :: wrt_bgc_diag_2d(nr_bgc_diag_2d)
  logical, target :: wrt_bgc_diag_3d(nr_bgc_diag_3d)
  integer(kind=4), target :: idx_bgc_diag_2d(nr_bgc_diag_2d)
  integer(kind=4), target :: idx_bgc_diag_3d(nr_bgc_diag_3d)
  ! Arrays storing information (name, unit, longname, fill value) about each diagnostic variable:
  character(len=72), target ::  vname_bgc_diag_2d(4,nr_bgc_diag_2d)
  character(len=72), target ::  vname_bgc_diag_3d(4,nr_bgc_diag_3d)

  real(kind=8),allocatable,dimension(:,:,:,:),target :: bgc_diag_3d_avg
  real(kind=8),allocatable,dimension(:,:,:)  ,target :: bgc_diag_2d_avg

#endif /* BEC2_DIAG  || MARBL_DIAGS */

  public bgc_precheck, read_nml_bgc

contains

  subroutine bgc_precheck

    implicit none

    character(len=12) :: sr_name = "bgc_precheck"
    character(len=1024) :: error_info

#if defined MARBL && defined BIOLOGY_BEC2
    write(error_info,*) "Both MARBL and BEC2 have been enabled"&
    &" in cppdefs.opt, but only one (or none) may be used in any"&
    &" given run."
    call error_log%raise_global(&
    &context=module_name//"/"//sr_name,&
    &info=error_info)
    call error_log%abort_check()
#endif

  end subroutine bgc_precheck


  integer(kind=4) function bgc_idx(t_idx)
    ! Takes a tracer array tracer index t(:,:,:,:,t_idx) and returns
    ! the corresponding bgc tracer index.
    ! Were previously harcoded individually as, e.g.
    ! po4_ind_t     = iPO4-iTandS-nt_passive
    integer(kind=4), intent(in) :: t_idx
    bgc_idx = t_idx-iTandS-nt_passive

  end function bgc_idx

  real(kind=8) function t_idx(b_idx)
    ! Takes a bgc tracer index (e.g. dtracer_module(:,:,:,b_idx)) and
    ! returns the corresponding tracer array index t(:,:,:,:,t_idx)
    integer(kind=4), intent(in) ::b_idx
    t_idx = b_idx+iTandS+nt_passive
  end function t_idx


#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
  subroutine find_write_bgc_diag
    implicit none

    integer(kind=4) idiag
#if defined(MARBL_DIAGS)

    wrt_marbl_diag_2d => wrt_bgc_diag_2d
    wrt_marbl_diag_3d => wrt_bgc_diag_3d
    idx_marbl_diag_2d => idx_bgc_diag_2d
    idx_marbl_diag_3d => idx_bgc_diag_3d


    call marbldrv_configure_diagnostics(&
    &vname_bgc_diag_2d,vname_bgc_diag_3d&
    &)

#elif defined(BEC2_DIAG)

    wrt_bec2_diag_2d => wrt_bgc_diag_2d
    wrt_bec2_diag_3d => wrt_bgc_diag_3d
    idx_bec2_diag_2d => idx_bgc_diag_2d
    idx_bec2_diag_3d => idx_bgc_diag_3d
    vname_bec2_diag_2d => vname_bgc_diag_2d
    vname_bec2_diag_3d => vname_bgc_diag_3d

#include "BEC_2Ddiagnostics.h"
#include "BEC_3Ddiagnostics.h"

#endif /* MARBL_DIAGS */


    nr_bgc_wrdiag_2d=0
    do idiag=1,nr_bgc_diag_2d
      if (wrt_bgc_diag_2d(idiag)) then
        nr_bgc_wrdiag_2d=nr_bgc_wrdiag_2d+1
        idx_bgc_diag_2d(idiag)=nr_bgc_wrdiag_2d
      else
        idx_bgc_diag_2d(idiag)=0
      endif
    enddo

    nr_bgc_wrdiag_3d=0
    do idiag=1,nr_bgc_diag_3d
      if (wrt_bgc_diag_3d(idiag)) then
        nr_bgc_wrdiag_3d=nr_bgc_wrdiag_3d+1
        idx_bgc_diag_3d(idiag)=nr_bgc_wrdiag_3d
      else
        idx_bgc_diag_3d(idiag)=0
      endif
    enddo


  end subroutine find_write_bgc_diag  !]
#endif /* BEC2_DIAG || MARBL_DIAGS */

  subroutine read_nml_bgc

!     Read the "BGC_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=13) :: sr_name = "read_nml_bgc"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=BGC_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read BGC_SETTINGS section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    nbgc_flx = 0+ &
#ifdef PCO2AIR_FORCING
    &+1&
#endif
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
    &+1 &
#endif
    & + nbgc_flx
    close(namelist_unit)

  end subroutine read_nml_bgc


  subroutine init_arrays_bgc_shared_vars
    implicit none

#if defined(BEC2_DIAG) || defined (MARBL_DIAGS)
    call find_write_bgc_diag
    call display_bgc_output_settings_to_terminal

    allocate( bgc_diag_3d(GLOBAL_2D_ARRAY,nz,nr_bgc_wrdiag_3d) )
    allocate( bgc_diag_2d(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )

    bgc_diag_2d=0.0_8
    bgc_diag_3d=0.0_8
#if defined(MARBL_DIAGS)
    marbl_diag_2d => bgc_diag_2d
    marbl_diag_3d => bgc_diag_3d
#elif defined(BEC2_DIAG)
    bec2_diag_2d => bgc_diag_2d
    bec2_diag_3d => bgc_diag_3d
#endif

    if (wrt_bgc_dia_avg) then
      allocate( bgc_diag_3d_avg(GLOBAL_2D_ARRAY,nz,nr_bgc_wrdiag_3d) )
      allocate( bgc_diag_2d_avg(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )
      bgc_diag_2d_avg=0.0_8
      bgc_diag_3d_avg=0.0_8
    endif

    if (mynode == 0) then
      write(*,'(7x,A,I3,I3)') 'BGC diags allocation :: ',&
      &nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d
    endif

#endif /* BEC2_DIAG || MARBL_DIAGS */

#if defined(MARBL)
!     Allocate arrays for MARBL saved state variables
    call marbldrv_configure_saved_state
#elif defined(BIOLOGY_BEC2)
    call init_arrays_bec2_vars ! bec2 version
#endif

  end subroutine init_arrays_bgc_shared_vars

  subroutine display_bgc_output_settings_to_terminal()
    implicit none
    integer(kind=4) :: idx, tidx, lineidx, n_bgc_fields, n_bgc_wrt_fields
    character(len=200), allocatable :: stdoutlines(:)
    character(len=200) :: stdoutline

    n_bgc_fields = nt_bgc
    n_bgc_wrt_fields = 0
#if defined(MARBL_DIAGS) || defined(BEC2_DIAG)
    n_bgc_fields=n_bgc_fields+nr_bgc_diag_2d+nr_bgc_diag_3d
#endif /* BEC2_DIAG || MARBL_DIAGS */
    allocate(stdoutlines(n_bgc_fields))
    stdoutlines=''
    lineidx=1

    ! BGC tracers
    do idx = 1, nt_bgc
      tidx = t_idx(idx)      ! Tracer index corresponding to this BGC tracer
      write(stdoutline, '(11X, A, T50, L1, 4X, A)')&
      &trim(t_vname(tidx)),& ! vname
      &wrt_t(tidx),&      !Write (T/F)
      &trim(t_lname(tidx)) ! Long name
      if (wrt_t(tidx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
      stdoutlines(lineidx) = stdoutline
      lineidx=lineidx+1
    end do

#if defined(MARBL_DIAGS) || defined(BEC2_DIAG)
    ! BGC 2D diagnostics
    do idx = 1, nr_bgc_diag_2d
      write(stdoutline, '(11X, A, T50, L1, 4X, A)')&
      &trim(vname_bgc_diag_2d(1,idx)),& ! Short name
      &wrt_bgc_diag_2d(idx),& !Write (T/F)
      &trim(vname_bgc_diag_2d(2,idx)) ! Long name
      if (wrt_bgc_diag_2d(idx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
      stdoutlines(lineidx) = stdoutline
      lineidx=lineidx+1
    end do

    ! BGC 3D diagnostics
    do idx = 1, nr_bgc_diag_3d
      write(stdoutline, '(11X, A, T50, L1, 4X, A)')&
      &trim(vname_bgc_diag_3d(1,idx)),& ! Short name
      &wrt_bgc_diag_3d(idx),& ! Write (T/F)
      &trim(vname_bgc_diag_3d(2,idx)) ! Long name
      if (wrt_bgc_diag_3d(idx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
      stdoutlines(lineidx) = stdoutline
      lineidx=lineidx+1
    end do
#endif /* BEC2_DIAG || MARBL_DIAGS */

    if (mynode==0) then
      if (wrt_bgc_his) then
        write(*,'(/7x,2A,F6.1,2x,A,I4)')& ! write to terminal output in simulation pre-amble text which
        &'bgc :: history file ',& ! result variables are being stored
        &'output_period =', output_period_bgc_his,&
        &'recs/file =', nrpf_bgc_his
      end if
      if (wrt_bgc_avg) then
        write(*,'(/7x,2A,F6.1,2x,A,I4)')& ! write to terminal output in simulation pre-amble text which
        &'bgc :: average file ',& ! result variables are being stored
        &'output_period =', output_period_bgc_avg,&
        &'recs/file =', nrpf_bgc_avg
      end if
      write(*,'(9x,A)') 'bgc fields to be saved (T/F)'
      write(*,'(11x,A)') repeat('-',62)
      write(*, '(11x,A,T40,A,T64,A)')&
      &"Name","Write (T/F)","Long name"
      write(*,'(11x,A/)') repeat('-',62)
      do idx = 1,lineidx-1
        write (*, '(A)') trim(stdoutlines(idx))
      end do
      write(*,'(11x,A)') repeat('-',62)
      write(*,'(11x,A,I3)') 'TOTAL NO. OF BGC FIELDS: ', n_bgc_fields
      write(*,'(11x,A,I3)') 'FIELDS TO WRITE: ', n_bgc_wrt_fields
      write(*,'(11x,A)') repeat('-',62)
    end if

    deallocate(stdoutlines)

  end subroutine display_bgc_output_settings_to_terminal
#endif  /* BIOLOGY_BEC2 || MARBL  */
end module bgc_shared_vars
