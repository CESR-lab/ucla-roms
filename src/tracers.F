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

      use param
      use dimensions
      use surf_flux,  only: stflx                          ! surface tracer flux should possibly live in this module rather?
      use scalars                                          ! for 'FIRST_TIME_STEP' and nstp
      use nc_read_write
      use roms_read_write
      use netcdf
#ifdef MARBL
      use marbl_driver
#endif
      use mpi_exchanges
      implicit none
      private

      real,public,allocatable,dimension(:,:,:,:,:) :: t              ! array of tracers

                                                                     ! tracer netcdf variables (user input required in init_tracers)
      character*42, dimension(nt), public :: t_vname                 ! short name
      character*60, dimension(nt), public :: t_lname                 ! long name
      character*42, dimension(nt), public :: t_units                 ! tracer units
      character*47, dimension(nt), public :: t_tname                 ! tracer input time variable name
      logical, dimension(nt), public      :: wrt_t                   ! t/f output tracer.
      logical, dimension(nt), public      :: wrt_t_avg               ! t/f output avg tracer.
      logical, dimension(nt), public      :: wrt_t_dia = .false.     ! t/f diagnostics tracer.

      integer, dimension(nt)              :: t_ana_frc               ! whether surface flux is read in (0), or analytical (1)

      !-- Tracer netcdf variables as arrays/matrices of 'NT' length:
      ! Final tracer concentrations live in 't' in ocean3d
      ! Surface tracer flux lives in 'stflx' in surf_flux.F module.
      type (ncforce)    :: nc_t( NT )                                   ! array of ncvs for each tracer
      integer,public :: iTandS                                       ! combined index of temperature and salinity
                                                                     ! public as used in set_forces.F. Need to be parameter
!     integer,public :: nt_passive=0                                 ! total number of passive tracers
      integer        :: itot=0                                       ! index counter to count total number of tracers.
                                                                     ! this term prevents the need to hardcode tracer indices.
      integer        :: interp_t                                     ! interpolate forcing from coarser input grid (=1) or not (=0)

      real, public, allocatable, dimension(:,:,:,:) :: t_avg         ! tracer averages. Memory only allocated for tracers of interest.
      integer,      allocatable, dimension(:)       :: t_avg_2_NT    ! convert index from t_avg(itavg) to t(itrc) since averages array may be smaller

      integer,dimension(NT),public :: NT_2_t_avg                     ! convert index from t(itrc) to t_avg(itavg) since averages array may be smaller
      integer,public        :: n_t_avg                               ! number of tracer averages to output <= NT

#if defined(BIOLOGY_BEC2) && !defined(MARBL)
# include "bgc_tracers_indx.opt"
#endif

      public set_surf_tracer_flx
      public init_tracers
      public def_his_trc
      public wrt_his_trc
      public wrt_rst_trc
      public def_avg_trc
      public wrt_avg_trc
      public set_avg_trc

      contains

! ----------------------------------------------------------------------
      subroutine set_surf_tracer_flx ![
      ! set tracer flux at surface

      implicit none

      ! local
      integer           :: itrc       ! tracer number for loop index
      character(len=46) :: t_flx_name ! Tracer time name

      do itrc=iTandS+1,nt

        if (t_ana_frc(itrc)==0) then ! Read in forcing data (not analytical)

          call set_frc_data(nc_t(itrc), stflx(:,:,itrc),'r' )

        elseif(t_ana_frc(itrc)==1) then ! Analytical forcing

          call set_ana_surf_tracer_flx(itrc)

        else

          if(mynode==0) then
            write(*,'(/1x,2A,I3,2A/)') 'ERROR: set_surf_tracer_flx:: ',
     &      'Forcing type not supported: t_ana_frc(itrc)= ', t_ana_frc(itrc),
     &      ', for tracer: ', nc_t(itrc)%vname
          endif
          error stop

        endif ! if(t_ana_frc(itrc)==0)

!       call exchange_xxx(stflx)

      enddo

      end subroutine set_surf_tracer_flx  !]
! ----------------------------------------------------------------------
      subroutine set_ana_surf_tracer_flx(itrc)  ![
      ! Set analytical surface tracer flux
      implicit none

      ! input/outputs
      integer itrc ! Current tracer index number

      ! local
      integer i,j
      real    u_pist

      ! Currently set up for zero surface tracer flux
      ! Also, time invariant so only set once for efficiency (doesn't change each timestep)
      if (FIRST_TIME_STEP) then ! Only first timestep

        do j=0,ny+1
          do i=0,nx+1
            stflx(i,j,itrc)= 0.
          enddo
        enddo

      endif

!     Could have various switches here for the different tracers, e.g.:
!     This should be an include file instead: ana_frc_trc.h
!     u_pist = 5.55e-5 ! 0.20 cm/hour
!     u_pist = u_pist/20.  ! linearization of carbon chemistry
!     if(itrc==3) then
!        do j=0,ny+1
!          do i=0,nx+1
!            stflx(i,j,itrc)= -u_pist*t(i,j,nz,nrhs,itrc)
!          enddo
!        enddo
!     endif

      end subroutine set_ana_surf_tracer_flx  !]
! ----------------------------------------------------------------------
      subroutine def_his_trc( ncid )  ![
      ! Define history file variables in def_his.F

      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer itrc, ierr, varid

      do itrc=1,iTandS+nt_passive
        if (wrt_t(itrc)) then
          varid = nccreate(ncid,t_vname(itrc),
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',t_lname(itrc))
          ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
        endif
      enddo

      end subroutine def_his_trc  !]
! ----------------------------------------------------------------------
      subroutine def_avg_trc( ncid )  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer :: itrc, ierr, varid
      character(len=64) :: long_name

      do itrc=1,iTandS+nt_passive
        if (wrt_t_avg(itrc)) then
          long_name='averaged '/ /t_lname(itrc) ! Add averaged to long name
          varid = nccreate(ncid,t_vname(itrc),
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
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
      integer,             intent(in) :: ncid
      integer,dimension(:),intent(in) :: start
      ! local
      integer :: itrc

      do itrc=1,iTandS+nt_passive
        if (wrt_t(itrc)) call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nnew,itrc),start)
      enddo

      end subroutine wrt_his_trc  !]
! ----------------------------------------------------------------------
      subroutine wrt_avg_trc (ncid,start)  ![
      ! Write tracers to avg file
      implicit none

      ! inputs
      integer,             intent(in) :: ncid
      integer,dimension(:),intent(in) :: start

      ! local
      integer :: itrc, itavg

      do itrc=1,iTandS+nt_passive
        if (wrt_t_avg(itrc)) then
          itavg = NT_2_t_avg(itrc)                                   ! get respective index for t_avg(itavg) -> t(itrc)
          call ncwrite(ncid,t_vname(itrc),t_avg(i0:i1,j0:j1,:,itavg),start)
        endif
      enddo

      end subroutine wrt_avg_trc  !]
! ----------------------------------------------------------------------
      subroutine wrt_rst_trc (ncid, start)  ![
      ! write tracers to restart file
      ! restart includes bgc tracers for simplicity
      implicit none

      ! inputs
      integer, intent(in) :: ncid
      integer,dimension(:),intent(in) :: start
      ! local
      integer :: itrc
      do itrc=1,NT
        call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nnew,itrc),start)
      enddo

      end subroutine wrt_rst_trc  !]
! ----------------------------------------------------------------------
      subroutine set_avg_trc(coef)  ![

      implicit none

      ! inputs
      real    :: coef
      ! local
      integer :: itrc, itavg

      do itrc=1,iTandS+nt_passive
        if (wrt_t_avg(itrc)) then
          itavg = NT_2_t_avg(itrc)                         ! get respective index for t_avg(itavg) -> t(itrc)
          t_avg(i0:i1,j0:j1,:,itavg) = t_avg(i0:i1,j0:j1,:,itavg)     *(1-coef) +
     &                                     t(i0:i1,j0:j1,:,nstp,itrc) * coef
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
      integer :: cnt=0, itrc

      ! Core tracers - temp and salt:
      t_vname(itemp)='temp';        t_units(itemp)='Celsius'
      t_lname(itemp)='potential temperature'
      iTandS = 1         ! if only temp, no salt.
#ifdef SALINITY
      t_vname(isalt)='salt';        t_units(isalt)='PSU'
      t_lname(isalt)='salinity'
      iTandS = 2         ! if both temp and salt.
#endif
      itot=iTandS        ! set up counting for additional tracers

      ! Additional passive tracers:
#include "tracers.opt"


#ifdef MARBL
      call marbldrv_configure_tracers(
     &     itot,t_vname,t_lname,t_units,t_tname,wrt_t,wrt_t_avg,t_ana_frc)
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

      allocate( t(GLOBAL_2D_ARRAY,N,3,NT) )
      t=0.

      ! remove averages flag above but do wrt_file_avg flag over this to avoid any allocation
      ! Allocate memory for only the tracer averages required for output
      do itrc=1,NT
        if (wrt_t_avg(itrc)) then
          cnt=cnt+1                                    ! count tracer averagesto calc+write
          NT_2_t_avg(itrc)=cnt                         ! t(itrc) = t_avg(cnt) - to convert between index scheme
        endif
      enddo

      n_t_avg = cnt                                    ! number of tracers of interest
      allocate( t_avg( GLOBAL_2D_ARRAY, N, n_t_avg) )  ! only for tracers we are interested in
      t_avg=0.

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

          nc_t(itrc)%vname = trim(t_vname(itrc)) / / '_flx' ! Forcing file flux name
          nc_t(itrc)%tname = t_tname(itrc)

        endif
      enddo

      end subroutine init_tracers  !]

! ----------------------------------------------------------------------

      end module tracers

