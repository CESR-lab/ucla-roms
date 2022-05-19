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
      use read_write, only: set_frc_var_tile, ncvs
      use surf_flux,  only: stflx                          ! surface tracer flux should possibly live in this module rather?
      use ncvars
      use scalars                                          ! for 'FIRST_TIME_STEP' and nstp
      use nc_read_write
      use roms_read_write
      use netcdf
      implicit none
      private

      !-- Tracer indices specific to each tracer
      ! Set to zero as possible means to avoid errors later? Might not be
      ! Necessary. itmp -> itemp, had to spell differently for now as
      ! itemp lives in param still (avoid clash). Same for islt -> isalt
      integer :: itmp=0, islt=0

      real,public,allocatable,dimension(:,:,:,:,:) :: t              ! array of tracers

                                                                     ! tracer netcdf variables (user input required in init_tracers)
      character*42, dimension(nt), public :: t_vname                 ! short name
      character*60, dimension(nt), public :: t_lname                 ! long name
      character*42, dimension(nt), public :: t_units                 ! tracer units
      character*47, dimension(nt), public :: t_tname                 ! tracer input time variable name
      logical, dimension(nt), public      :: wrt_t                   ! t/f output tracer.
      logical, dimension(nt), public      :: wrt_t_avg               ! t/f output avg tracer.
      integer, dimension(nt)              :: t_ana_frc               ! whether surface flux is read in (0), or analytical (1)

      logical, dimension(nt), public      :: wrt_t_diag = .false.    ! only with DIAGNOSTICS. T/F output tracer diagnostics (from diagnostics.F). Declare .false. so that doesn't need to be filled in tracers_defs.h if not used.

      !-- Tracer netcdf variables as arrays/matrices of 'NT' length:
      ! Final tracer concentrations live in 't' in ocean3d
      ! Surface tracer flux lives in 'stflx' in surf_flux.F module.
      type (ncvs)    :: nc_t( NT )                                   ! array of ncvs for each tracer
      integer,public :: iTandS                                       ! combined index of temperature and salinity
                                                                     ! public as used in set_forces.F. Need to be parameter
      integer,public :: nt_passive=0                                 ! total number of passive tracers
      integer        :: itot=0                                       ! index counter to count total number of tracers.
                                                                     ! this term prevents the need to hardcode tracer indices.
      integer        :: interp_t                                     ! interpolate forcing from coarser input grid (=1) or not (=0)

      ! 3 below currently in param as per old code
!      integer :: itrc_bio=0 ! Starting tracer index for bgc tracers. Set in tracers_defs.h
      ! Number of 'base' bgc tracers, i.e. minimum number of bgc tracers without cocco
      ! and nitrogen cycle additions
!      integer :: ntrc_bio_base
!      integer :: ntrc_bio !  - from (param - ntrc_bio=ntrc_bio_base+ntrc_bio_cocco+ntrc_bio_ncycle)

#ifdef AVERAGES
      real, public, allocatable, dimension(:,:,:,:) :: t_avg         ! tracer averages. Memory only allocated for tracers of interest.
      integer,      allocatable, dimension(:)       :: t_avg_2_NT    ! convert index from t_avg(itavg) to t(itrc) since averages array may be smaller

      integer,dimension(NT),public :: NT_2_t_avg                     ! convert index from t(itrc) to t_avg(itavg) since averages array may be smaller
      integer,public        :: n_t_avg                               ! number of tracer averages to output <= NT
#endif

#ifdef BIOLOGY_BEC2
! bgc needs public indexes as used in other files
# include "bgc_tracers_indx.opt"
#endif

      public set_surf_tracer_flx
      public init_tracers
      public set_bry_tracer_defaults
      public def_his_trc
      public wrt_his_trc
      public wrt_rst_trc
#ifdef AVERAGES
      public def_avg_trc
      public wrt_avg_trc
      public set_avg_trc
#endif
      contains

! ----------------------------------------------------------------------
      subroutine set_surf_tracer_flx(istr,iend,jstr,jend) ![ set tracer flux at surface

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! local
      integer           :: itrc       ! tracer number for loop index
      character(len=46) :: t_flx_name ! Tracer time name

      do itrc=iTandS+1,NT

        if (t_ana_frc(itrc)==0) then ! Read in forcing data (not analytical)

          call set_frc_var_tile( istr, iend, jstr, jend, nc_t(itrc), stflx(:,:,itrc), interp_t )

        elseif(t_ana_frc(itrc)==1) then ! Analytical forcing

          call set_ana_surf_tracer_flx(istr,iend,jstr,jend,itrc)

        else

          if(mynode==0) then
            write(*,'(/1x,2A,I3,2A/)') 'ERROR: set_surf_tracer_flx:: ',
     &      'Forcing type not supported: t_ana_frc(itrc)= ', t_ana_frc(itrc),
     &      ', for passive tracer: ', t_vname(itrc)
          endif
          error stop

        endif ! if(t_ana_frc(itrc)==0)

      enddo

      end subroutine set_surf_tracer_flx  !]

! ----------------------------------------------------------------------
      subroutine set_ana_surf_tracer_flx(istr,iend,jstr,jend,itrc)  ![
      ! Set analytical surface tracer flux
      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      integer itrc ! Current tracer index number

      ! local
      integer i,j

#include "compute_extended_bounds.h"

      ! Currently set up for zero surface tracer flux
      ! Also, time invariant so only set once for efficiency (doesn't change each timestep)
      if (FIRST_TIME_STEP) then ! Only first timestep

        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,itrc)= 0.
          enddo
        enddo

      endif

!     Could have various switches here for the different tracers, e.g.:
!     if(itrc==iptrace1) then
!        do j=jstrR,jendR
!          do i=istrR,iendR
!            stflx(i,j,itrc)=sin(time)*100. ! Made up example
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
          varid = nccreate(ncid,t_vname(itrc),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
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
          varid = nccreate(ncid,t_vname(itrc),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
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
        if (wrt_t(itrc)) call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nstp,itrc),start)
      enddo

      end subroutine wrt_his_trc  !]

! ----------------------------------------------------------------------
#ifdef AVERAGES
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
#endif /* AVERAGES */

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
        call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nstp,itrc),start)
      enddo

      end subroutine wrt_rst_trc  !]

! ----------------------------------------------------------------------
#ifdef AVERAGES
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
#endif /* AVERAGES */

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
      itmp=1;
      t_vname(itmp)='temp';        t_units(itmp)='Celsius'
      t_lname(itmp)='potential temperature'
      iTandS = 1         ! if only temp, no salt.
#ifdef SALINITY
      islt=2;
      t_vname(islt)='salt';        t_units(islt)='PSU'
      t_lname(islt)='salinity'
      iTandS = 2         ! if both temp and salt.
#endif
      itot=iTandS        ! set up counting for additional tracers

      ! Additional passive tracers:
#include "tracers.opt"

      allocate( t(GLOBAL_2D_ARRAY,N,3,NT) )
      t=0.

#if defined AVERAGES
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

#endif  /* AVERAGES */

      ! initialize read in forcing data arrays
      do itrc=iTandS+1,NT
        if (t_ana_frc(itrc)==0) then

          allocate( nc_t(itrc)%data( GLOBAL_2D_ARRAY,2) )

          ! set nc_v%name and nc_t%time_name only once: currently set in t_vname & t_tname,
          ! left it like this so people don't need to change their tracers.opt files.

          nc_t(itrc)%name      = trim(t_vname(itrc)) / / '_flx' ! Forcing file flux name
          nc_t(itrc)%time_name = t_tname(itrc)

        endif
      enddo

      end subroutine init_tracers  !]

! ----------------------------------------------------------------------
      subroutine set_bry_tracer_defaults(t_bry_dt, itrc, itbry, LMmpp)  ![
      ! Sets default values for tracers at the boundary if missing
      ! in the _bry.nc file.
      ! In the long run this should be put into an include file to avoid
      ! conflicts when module is updated.
      ! Taken from Pierre's code get_bry_all.F

      use param
      implicit none
      real t_bry_dt(0:LMmpp,N,2,NT), cff,cff1
      integer, intent(in):: LMmpp, itbry,itrc
      integer :: ij, k
#  define temp cff
#  define SiO4 cff1

      if (itrc == itemp) then
        t_bry_dt(:,:,itbry,itrc) = 5.0 ! default temperature (only used if temp not in boundary file
#  ifdef SALINITY
      else if (itrc == isalt) then
        t_bry_dt(:,:,itbry,itrc) = 35.0 ! default salinity
#  endif
#  if defined BIOLOGY_BEC2
! provide defaults for selected biogeochemical variables
      elseif (itrc .eq. iNH4 .or. itrc .eq. iDOC .or.
     &    itrc .eq. iDON .or. itrc. eq. iDOFE .or.
     &    itrc .eq. iDOP) then
        t_bry_dt(:,:,itbry,itrc) = 0.0
      else if (itrc .eq. iPO4) then
        do k=1,N
          do ij=0,LMmpp
            temp=t_bry_dt(ij,k,itbry,itemp)
            if (temp.lt.8.) then
              SiO4=30.
            elseif (temp.ge.8. .and. temp.le.11.) then
              SiO4=30.-((temp-8.)*(20./3.))
            elseif (temp.gt.11. .and. temp.le.13.) then
              SiO4=10.-((temp-11.)*(8./2.))
            elseif (temp.gt.13. .and. temp.le.16.) then
              SiO4=2.-((temp-13.)*(2./3.))
            elseif (temp.gt.16.) then
              SiO4=0.
            endif
          t_bry_dt(ij,k,itbry,itrc)=(1.67+0.5873*SiO4
     &         +0.0144*SiO4**2
     &         +0.0003099*SiO4**3)/16
#   ifdef NO3_SFC_CUTOFF
          t_bry_dt(ij,k,itbry,itrc)=t_bry_dt(ij,k,itbry,itrc)*
     &                             (1-.5*(tanh((z_r(iwest,ij,k)+100.)/8.)+1))
#   endif /* NO3_SFC_CUTOFF */
          enddo
        enddo
      ! DevinD - would be better to use switch statement here rather than if's
      else if (itrc .eq. iNO3) then
        t_bry_dt(:,:,itbry,itrc) = 16*t_bry_dt(:,:,itbry,iPO4)
      else if (itrc .eq. iSiO3) then
        t_bry_dt(:,:,itbry,itrc) = t_bry_dt(:,:,itbry,iNO3)
      else if (itrc .eq. iFE) then
        t_bry_dt(:,:,itbry,itrc) = 5.395e-4  ! mmol/m3
      else if (itrc .eq. iO2) then
        t_bry_dt(:,:,itbry,itrc) = 177.9  ! mmol/m3
      else if (itrc .eq. iDIC) then
        t_bry_dt(:,:,itbry,itrc) = 2255  ! mmol/m3
      else if (itrc .eq. iALK) then
        t_bry_dt(:,:,itbry,itrc) = 2363  ! mmol/m3
      else if (itrc .eq. iSPC) then
        t_bry_dt(:,:,itbry,itrc) = 0.10
      else if (itrc .eq. iSPCHL) then
        t_bry_dt(:,:,itbry,itrc) = 0.1 *
     &     t_bry_dt(:,:,itbry,iSPC)
      else if (itrc .eq. iSPFE) then
        t_bry_dt(:,:,itbry,itrc) = 4.d-5 *
     &     t_bry_dt(:,:,itbry,iSPC)
      else if (itrc .eq. iSPCACO3) then
        t_bry_dt(:,:,itbry,itrc) = 0.03 *
     &     t_bry_dt(:,:,itbry,iSPC)
      else if (itrc .eq. iDIATC) then
        t_bry_dt(:,:,itbry,itrc) = 0.10
      else if (itrc .eq. iDIATCHL) then
        t_bry_dt(:,:,itbry,itrc) = 0.1 *
     &     t_bry_dt(:,:,itbry,iDIATC)
      else if (itrc .eq. iDIATSI) then
        t_bry_dt(:,:,itbry,itrc) = 0.2 *
     &     t_bry_dt(:,:,itbry,iDIATC)
      else if (itrc .eq. iDIATFE) then
        t_bry_dt(:,:,itbry,itrc) = 3.d-5 *
     &     t_bry_dt(:,:,itbry,iDIATC)
      else if (itrc .eq. iDIAZC) then
        t_bry_dt(:,:,itbry,itrc) = 0.01
      else if (itrc .eq. iDIAZCHL) then
        t_bry_dt(:,:,itbry,itrc) = 0.1 *
     &     t_bry_dt(:,:,itbry,iDIAZC)
      else if (itrc .eq. iDIAZFE) then
        t_bry_dt(:,:,itbry,itrc) = 3.d-5 *
     &     t_bry_dt(:,:,itbry,iDIAZC)
      else if (itrc .eq. iZOOC) then
        t_bry_dt(:,:,itbry,itrc) = 0.05
      else if (itrc >= itrc_bio) then
        t_bry_dt(:,:,itbry,itrc) = 1e-30
#    endif /* BIOLOGY_BEC2 */
      else
        t_bry_dt(:,:,itbry,itrc) = 0.0
      endif ! tracer cases


      end subroutine set_bry_tracer_defaults  !]

! ----------------------------------------------------------------------

      end module tracers

