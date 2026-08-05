module tracers

  ! initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)
  !
  ! INFO:  ![
  ! 1) in order to specify tracer variable details, user must
  ! set them in tracers.opt (NT still set in param for now).
  !
  ! 2) tracer name in netcdf file for _bry.nc and _init.nc
  ! must have same name as trace variable. ptrace1 -> ptrace1
  ! However, in the surface flux forcing file, the name of the
  ! variable must be appended with _flx. E.g. ptrace1 -> ptrace1_flx
  !
  ! 3) Make sure tracer surface flux units are correct!
  !
  ! BGC: bgc tracers are stored in the tracer array, but the outputting
  ! of bgc tracers is done to a seperate bgc output file in bgc.F
  !]

#include "cppdefs.opt"
  use param, only: isalt, itemp, lm, mm, mynode, nt_passive, nt_cdr_oae, nt_cdr_dor
  use dimensions, only: i0, i1, j0, j1, nx, ny, eta_rho, xi_rho&
  &, ds_xr, ds_yr, ds_zr
  use surf_flux, only: stflx                          ! surface tracer flux should possibly live in this module rath
#ifdef CDR_TRACER
  use surf_flux, only: ddic_dco2, ddic_dalk, k_gas, uwnd, vwnd
  use ocean_vars, only: Hz
#endif
  use scalars, only: nz, nstp, nnew, forw_start, iic, nt
! for 'FIRST_TIME_STEP' and nstp, only:
  use nc_read_write, only: nccreate, ncwrite
  use roms_read_write, only:&
  &ncforce, dn_tm, dn_xr, dn_yr, dn_zr,&
  &set_frc_data
  use netcdf, only: nf90_put_att, nf90_double
#ifdef MARBL
  use marbl_driver, only: marbldrv_configure_tracers
#endif
#ifdef PARALLEL_IO
  use pio_roms, only: pio_file_is_open, pio_FileDesc
  use pio, only : PIO_closefile
#endif

  implicit none
  private
  character(len=7) :: module_name = "tracers"
  character(len=1024) :: error_info = ""
  real(kind=8),public,allocatable,dimension(:,:,:,:,:) :: t              ! array of tracers

  ! tracer netcdf variables (user input required in i
  character(len=42), dimension(:), public, allocatable :: t_vname                 ! short name
  character(len=60), dimension(:), public, allocatable :: t_lname                 ! long name
  character(len=42), dimension(:), public, allocatable :: t_units                 ! tracer units
  character(len=47), dimension(:), public, allocatable :: t_tname                 ! tracer input time variable name
  logical, dimension(:), public, allocatable      :: wrt_t                   ! t/f output tracer.
  logical, dimension(:), public, allocatable      :: wrt_t_avg               ! t/f output avg tracer.
  logical, dimension(:), public, allocatable      :: wrt_t_dia      ! t/f diagnostics tracer.

  integer(kind=4), dimension(:), allocatable              :: t_ana_frc               ! whether surface flux is read

  integer, dimension(:), public, allocatable :: itrc_alk_pair ! Pairing logic

  !-- Tracer netcdf variables as arrays/matrices of 'NT' length:
  ! Final tracer concentrations live in 't' in ocean3d
  ! Surface tracer flux lives in 'stflx' in surf_flux.F module.
  type (ncforce), allocatable    :: nc_t( : )                                   ! array of ncvs for each tracer
  integer(kind=4),public :: iTandS                                       ! combined index of temperature and salinit
  ! public as used in set_forces.F. Need to be parame
!     integer,public :: nt_passive=0                                 ! total number of passive tracers
  integer(kind=4)        :: itot=0                                       ! index counter to count total number of tr
  ! this term prevents the need to hardcode tracer in
  integer(kind=4)        :: interp_t                                     ! interpolate forcing from coarser input gr

  real(kind=8), public, allocatable, dimension(:,:,:,:) :: t_avg         ! tracer averages. Memory only allocated fo
  integer(kind=4),      allocatable, dimension(:)       :: t_avg_2_NT    ! convert index from t_avg(itavg) to t(itrc

  integer(kind=4),dimension(:),public, allocatable :: NT_2_t_avg                     ! convert index from t(itrc) to
  integer(kind=4),public        :: n_t_avg                               ! number of tracer averages to output <= NT

  logical, public :: wrt_temp, wrt_salt, wrt_temp_dia, wrt_salt_dia
  namelist /TS_OUTPUT_SETTINGS/ wrt_temp, wrt_salt, wrt_temp_dia, wrt_salt_dia
#if defined(BIOLOGY_BEC2) && !defined(MARBL)
# include "BEC_tracers_indx.h"
#endif

  public set_surf_tracer_flx
  public init_tracers
  public def_his_trc
  public wrt_his_trc
  public wrt_rst_trc
  public def_avg_trc
  public wrt_avg_trc
  public set_avg_trc
  public read_nml_tracers
contains

!     ----------------------------------------------------------------------
  subroutine read_nml_tracers
    use error_handling_mod, only: error_log
    use namelist_open_mod, only: open_namelist_file
!     Read the "TS_OUTPUT_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_tracers"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=TS_OUTPUT_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read TS_OUTPUT_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)

  end subroutine read_nml_tracers

  subroutine set_surf_tracer_flx ![
    ! set tracer flux at surface
    use error_handling_mod, only: error_log
    implicit none
    character(len=20) :: sr_name = "set_surf_tracer_flx"
    ! local
    integer(kind=4)           :: itrc       ! tracer number for loop index
    character(len=46) :: t_flx_name ! Tracer time name

#ifdef CDR_TRACER
      call exchange_xxx(t(:,:,nz,nrhs,itemp) )
      call set_gas_transfer_velocity(istr,iend,jstr,jend)
#endif

#ifdef PARALLEL_IO
    pio_file_is_open = 0
#endif

    do itrc=iTandS+1,nt

      if (t_ana_frc(itrc)==0) then ! Read in forcing data (not analytical)

        call set_frc_data(nc_t(itrc), stflx(:,:,itrc),'r' )

      elseif(t_ana_frc(itrc)==1) then ! Analytical forcing

        call set_ana_surf_tracer_flx(itrc)

#ifdef CDR_TRACER
      elseif (t_ana_frc(itrc)==2) then
            call set_frc_data(nc_t(itrc), stflx(:,:,itrc), 'r' )
            call exchange_xxx(stflx(:,:,itrc))

            call exchange_xxx(t(:,:,nz,nrhs,itrc) )
            call subtract_gas_exchange_from_tracer_flx(itrc,istr,iend,jstr,jend)
            call exchange_xxx(stflx(:,:,itrc))
#endif

      else

        write(error_info, *)&
        &'Forcing type not supported: t_ana_frc(itrc)= ', t_ana_frc(itrc),&
        &', for tracer: ', nc_t(itrc)%vname
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info=error_info)

      endif ! if(t_ana_frc(itrc)==0)
!       call exchange_xxx(stflx)

    enddo
#ifdef PARALLEL_IO
    if (pio_file_is_open==1) then
      call PIO_closefile(pio_FileDesc)
    endif
    pio_file_is_open = 0
#endif

    call error_log%abort_check()
  end subroutine set_surf_tracer_flx  !]
! ----------------------------------------------------------------------
  subroutine set_ana_surf_tracer_flx(itrc)  ![
    ! Set analytical surface tracer flux
    implicit none

    ! input/outputs
    integer(kind=4) itrc ! Current tracer index number

    ! local
    integer(kind=4) i,j
    real(kind=8)    u_pist

    ! Currently set up for zero surface tracer flux
    ! Also, time invariant so only set once for efficiency (doesn't change each timestep)
    if (FIRST_TIME_STEP) then ! Only first timestep

      do j=0,ny+1
        do i=0,nx+1
          stflx(i,j,itrc)= 0._8
        enddo
      enddo

    endif

!     Could have various switches here for the different tracers, e.g.:
!     This should be an include file instead: ana_frc_trc.h
!     u_pist = 5.55D-5 ! 0.20 cm/hour
!     u_pist = u_pist/20._8  ! linearization of carbon chemistry
!     if(itrc==3) then
!        do j=0,ny+1
!          do i=0,nx+1
!            stflx(i,j,itrc)= -u_pist*t(i,j,nz,nrhs,itrc)
!          enddo
!        enddo
!     endif

  end subroutine set_ana_surf_tracer_flx  !]
! ----------------------------------------------------------------------
      subroutine set_gas_transfer_velocity(istr,iend,jstr,jend) ![

      implicit none
      integer istr, iend, jstr, jend

#if defined CDR_TRACER
      ! coefficients to compute Schmidt number
      real, parameter :: a = 2116.8
      real, parameter :: b = -136.25
      real, parameter :: c = 4.7353
      real, parameter :: d = -0.092307
      real, parameter :: e = 0.0007555
      real, parameter :: xkw_coef = 6.97e-07
      real schmidt_nr, sst
      integer i, j

      do j=jstr,jend
        do i=istr,iend+1

          sst = t(i,j,nz,nrhs,itemp)
          schmidt_nr = a + sst * (b + sst * (c + sst * (d + e * sst)))
          k_gas(i,j) = xkw_coef * (uwnd(i,j)**2 + vwnd(i,j)**2) * sqrt(660.0 / schmidt_nr)
#ifdef SEA_ICE_NOFLUX
          if( sst .le. -1.8 ) then
              k_gas(i,j)=0.                    ! If SST colder than -1.8 C, assume zero piston velocity due to sea ice
          endif
#endif
        enddo
      enddo

#endif /* defined CDR_TRACER */

      end subroutine set_gas_transfer_velocity  !]
! ----------------------------------------------------------------------
      subroutine subtract_gas_exchange_from_tracer_flx(itrc,istr,iend,jstr,jend)  ![
      ! Modify surface tracer flux by subtracting the air-sea gas exchange term

      implicit none

      integer itrc ! Current tracer index number
      integer istr, iend, jstr, jend

#ifdef CDR_TRACER

      ! local
      integer i, j, iALK
      real beta, eta, cALK

      iALK = itrc_alk_pair(itrc) ! Identify if this tracer has a pair

      do j=jstr,jend
        do i=istr,iend+1
          beta = max(ddic_dco2(i,j), 1.0e-5) ! avoid division by zero
          eta = ddic_dalk(i,j)
          if (iALK > 0) then
             cALK = t(i,j,nz,nrhs,iALK)
          else
             cALK = 0.0
          endif

          ! Implementation: Flux = Flux - (k / beta) * (C_dic - ddic_dalk * C_alk)
          stflx(i,j,itrc) = stflx(i,j,itrc) -
     &      ( k_gas(i,j) / beta ) * ( t(i,j,nz,nrhs,itrc) - eta * cALK )

          ! Here, stflx is stored as (stflx * Hz) to maintain consistency
          ! with t, which is also in (t * Hz) form. This ensures the units match
          ! when stflx is added to t in step3d_t_ISO.
          ! After the implicit vertical mixing step in step3d_t_ISO, (t * Hz)
          ! is divided by Hz to yield the updated tracer t.
        enddo
      enddo
#endif
      end subroutine subtract_gas_exchange_from_tracer_flx  !]
! ----------------------------------------------------------------------
  subroutine def_his_trc( ncid )  ![
    ! Define history file variables in def_his.F

    implicit none

    ! input
    integer(kind=4),intent(in) :: ncid
    ! local
    integer(kind=4) itrc, ierr, varid

    do itrc=1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
      if (wrt_t(itrc)) then
        varid = nccreate(ncid,t_vname(itrc),&
        &(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',t_lname(itrc))
        ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
      endif
    enddo

  end subroutine def_his_trc  !]
! ----------------------------------------------------------------------
  subroutine def_avg_trc( ncid )  ![
    implicit none

    ! input
    integer(kind=4),intent(in) :: ncid
    ! local
    integer(kind=4) :: itrc, ierr, varid
    character(len=64) :: long_name

    do itrc=1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
      if (wrt_t_avg(itrc)) then
        long_name='averaged '//t_lname(itrc) ! Add averaged to long name
        varid = nccreate(ncid,t_vname(itrc),&
        &(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',long_name)
        ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
      endif
    enddo

  end subroutine def_avg_trc  !]
! ----------------------------------------------------------------------
  subroutine wrt_his_trc (ncid,start)  ![
    ! write tracers to history file
    ! temp+salinity and passive tracers (no bgc)
    implicit none

    ! inputs
    integer(kind=4),             intent(in) :: ncid
    integer(kind=4),dimension(:),intent(in) :: start
    ! local
    integer(kind=4) :: itrc

    do itrc=1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
      if (wrt_t(itrc)) call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nnew,itrc),start, .true.)
    enddo

  end subroutine wrt_his_trc  !]
! ----------------------------------------------------------------------
  subroutine wrt_avg_trc (ncid,start)  ![
    ! Write tracers to avg file
    implicit none

    ! inputs
    integer(kind=4),             intent(in) :: ncid
    integer(kind=4),dimension(:),intent(in) :: start

    ! local
    integer(kind=4) :: itrc, itavg

    do itrc=1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
      if (wrt_t_avg(itrc)) then
        itavg = NT_2_t_avg(itrc)                                   ! get respective index for t_avg(itavg) -> t(itrc)
        call ncwrite(ncid,t_vname(itrc),t_avg(i0:i1,j0:j1,:,itavg),start, .true.)
      endif
    enddo

  end subroutine wrt_avg_trc  !]
! ----------------------------------------------------------------------
  subroutine wrt_rst_trc (ncid, start)  ![
    ! write tracers to restart file
    ! restart includes bgc tracers for simplicity
    implicit none

    ! inputs
    integer(kind=4), intent(in) :: ncid
    integer(kind=4),dimension(:),intent(in) :: start
    ! local
    integer(kind=4) :: itrc

    do itrc=1,NT
      call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nnew,itrc),start,.true.)
    enddo

  end subroutine wrt_rst_trc  !]
! ----------------------------------------------------------------------
  subroutine set_avg_trc(coef)  ![

    implicit none

    ! inputs
    real(kind=8)    :: coef
    ! local
    integer(kind=4) :: itrc, itavg

    do itrc=1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
      if (wrt_t_avg(itrc)) then
        itavg = NT_2_t_avg(itrc)                         ! get respective index for t_avg(itavg) -> t(itrc)
        t_avg(i0:i1,j0:j1,:,itavg) = t_avg(i0:i1,j0:j1,:,itavg)     *(1-coef) +&
        &t(i0:i1,j0:j1,:,nstp,itrc) * coef
      endif
    enddo

  end subroutine set_avg_trc  !]
! ----------------------------------------------------------------------
  subroutine init_tracers  ![

    ! Set all tracer variable values. This was placed at bottom of
    ! the module for ease of reading the rest of the module's code.

    ! It is necessary to keep track of tracer index number in order to
    ! correctly locate and calculate upon the tracer in 't' the array
    ! of all the tracers.

    ! Attempted to atleast here keep all variables together.
    ! Yes, still need to declare them in module preamble, however,
    ! they take no value there so order/value is not a worry.
    ! Old code: index is set in param, and variable values in
    ! init_scalars.F.
    ! New method, both set here.

    ! save
    ! SHOULD DECLARE TRACER INDICES IN TRACERS_DEFS.H & USE THE SAVE COMMAND HERE???

    implicit none

    ! local
    integer(kind=4) :: cnt=0, itrc
    character(len=8) :: passive_tracer_num
    allocate(t_vname(nt))
    allocate(t_lname(nt))
    allocate(t_units(nt))
    allocate(t_tname(nt))
    allocate(wrt_t(nt))
    allocate(wrt_t_avg(nt))
    allocate(wrt_t_dia(nt))
    allocate(t_ana_frc(nt))
    allocate(nc_t(nt))
    allocate(NT_2_t_avg(nt))
    allocate(itrc_alk_pair(nt))
    wrt_t(:) = .false.
    wrt_t_avg(:) = .false.
    wrt_t_dia(:) = .false.

    itrc_alk_pair(:) = 0

    ! Core tracers - temp and salt:
    t_vname(itemp)='temp';        t_units(itemp)='Celsius'
    t_lname(itemp)='potential temperature'
    if(wrt_temp) wrt_t(itemp) = .true.
    if(wrt_temp_dia) wrt_t_dia(itemp) = .true.
    iTandS = 1                ! if only temp, no salt.
#ifdef SALINITY
    t_vname(isalt)='salt';        t_units(isalt)='PSU'
    t_lname(isalt)='salinity'
    if(wrt_salt) wrt_t(isalt) = .true.
    if(wrt_salt_dia) wrt_t_dia(isalt) = .true.
    iTandS = 2         ! if both temp and salt.
#endif
    itot=iTandS        ! set up counting for additional tracers

    do itrc=iTandS+1,iTandS+nt_passive
      write(passive_tracer_num, '(I0)') (itrc-iTandS)
      t_vname(itrc)='passive_tracer' // TRIM(passive_tracer_num)
      t_units(itrc)='mMol m-3'
      t_lname(itrc)='passive tracer' // TRIM(passive_tracer_num)
      wrt_t(itrc) = .false.
      wrt_t_dia(itrc) = .false.
      t_ana_frc(itrc)=1
      itot = itot+1
    enddo

    do itrc=iTandS+nt_passive+1,iTandS+nt_passive+nt_cdr_oae,2
      write(passive_tracer_num, '(I0)') (itrc-iTandS-nt_passive)
      t_vname(itrc)='CDR_OAE_ALK' // TRIM(passive_tracer_num)
      t_units(itrc)='mMol m-3'
      t_lname(itrc)='CDR OAE ALK tracer' // TRIM(passive_tracer_num)
      wrt_t(itrc) = .false.
      wrt_t_dia(itrc) = .false.
      t_ana_frc(itrc)=2
      itot = itot+1

      write(passive_tracer_num, '(I0)') (itrc+1-iTandS-nt_passive)
      t_vname(itrc+1)='CDR_OAE_DIC' // TRIM(passive_tracer_num)
      t_units(itrc+1)='mMol m-3'
      t_lname(itrc+1)='CDR OAE DIC tracer' // TRIM(passive_tracer_num)
      wrt_t(itrc+1) = .false.
      wrt_t_dia(itrc+1) = .false.
      t_ana_frc(itrc+1)=2
      ! Identify the corresponding ALK tracer
      itrc_alk_pair(itrc+1) = itrc
      itot = itot+1
    enddo

    do itrc=iTandS+nt_passive+2*nt_cdr_oae+1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
      write(passive_tracer_num, '(I0)') (itrc-iTandS-nt_passive-2*nt_cdr_oae)
      t_vname(itrc)='CDR_DOR_DIC' // TRIM(passive_tracer_num)
      t_units(itrc)='mMol m-3'
      t_lname(itrc)='CDR DOR DIC tracer' // TRIM(passive_tracer_num)
      wrt_t(itrc) = .false.
      wrt_t_dia(itrc) = .false.
      t_ana_frc(itrc)=2
      itot = itot+1
    enddo

    ! Additional passive tracers:
#ifdef BIOLOGY_BEC2
#include "BEC_tracers.h"
#endif

#ifdef MARBL
    call marbldrv_configure_tracers(&
    &itot,t_vname,t_lname,t_units,t_tname,wrt_t,wrt_t_avg,t_ana_frc)
#endif

    if (mynode==0) then
      print *, 'metadata about ',NT, ' tracers:'
      do itrc=1,NT
        print *, '-----------'
        print *, 'TRACER NO.: ', itrc
        print *,'SHORT NAME: ', t_vname(itrc)
        print *,'LONG NAME: ', t_lname(itrc)
        print *,'UNITS: ', t_units(itrc)
        print *, '-----------'
      end do
    end if

    allocate( t(GLOBAL_2D_ARRAY,nz,3,NT) )
    t=0._8

    ! remove averages flag above but do wrt_file_avg flag over this to avoid any allocation
    ! Allocate memory for only the tracer averages required for output
    do itrc=1,NT
      if (wrt_t_avg(itrc)) then
        cnt=cnt+1                                    ! count tracer averagesto calc+write
        NT_2_t_avg(itrc)=cnt                         ! t(itrc) = t_avg(cnt) - to convert between index scheme
      endif
    enddo

    n_t_avg = cnt                                    ! number of tracers of interest
    allocate( t_avg( GLOBAL_2D_ARRAY, nz, n_t_avg) )  ! only for tracers we are interested in
    t_avg=0._8

    allocate( t_avg_2_NT( n_t_avg ) )                ! to convert indices between t_avg(itavg) and t(itrc)
    t_avg_2_NT = -1                                  ! set to bad number (<1) as safeguard

    cnt=0
    do itrc=1,NT
      if (wrt_t_avg(itrc)) then
        cnt=cnt+1
        t_avg_2_NT(cnt)=itrc                         ! store the actual tracer index of 't' array
      endif                                          ! since t(NT) but t_avg(n_t_avg). NT >= n_t_avg...
    enddo

    ! initialize read in forcing data arrays
    do itrc=iTandS+1,NT
      if (t_ana_frc(itrc)==0) then

        allocate( nc_t(itrc)%vdata( GLOBAL_2D_ARRAY,2) )

        ! set nc_v%vname and nc_t%tname only once: currently set in t_vname & t_tname,
        ! left it like this so people don't need to change their tracers.opt files.

        nc_t(itrc)%vname = trim(t_vname(itrc)) // '_flx' ! Forcing file flux name
        nc_t(itrc)%tname = t_tname(itrc)

      endif
    enddo

  end subroutine init_tracers  !]
! ----------------------------------------------------------------------
end module tracers

