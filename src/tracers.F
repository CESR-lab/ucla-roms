      module tracers

      ! initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)

#include "cppdefs.h"

      use read_write, only: set_frc_var_tile ! Modules required
      use surf_flux,  only: stflx ! Surface tracer flux should possibly
                                  ! live in this module rather?

      implicit none

      private ! All variables private to module unless public specified

      ! Includes
      ! ========

      ! param.h needed for GLOBAL_2D_array to work. NT = number tracer from param.h
#include "param.h"

      ! ****************************************************************
      ! User inputs

        ! 1) in order to specify tracer variable details, user must
        ! set them in tracers_defs.h. (NT still set in param.h for now)

        ! 2) tracer name in netcdf file for _bry.nc and _init.nc
        ! must have same name as trace variable. ptrace1 -> ptrace1
        ! However, in the surface flux forcing file, the name of the
        ! variable must be appended with _flx. E.g. ptrace1 -> ptrace1_flx

        ! 3) Make sure tracer surface flux units are correct!

      !-- Tracer indices specific to each tracer
      ! THIS SHOULD GO INTO TRACERS_DEFS.H AND JUST USE SAVE ATTRIBUTE SO
      ! THAT THEY ARE RETAINED AFTER INITIALIZING
      ! Set to zero as possible means to avoid errors later? Might not be
      ! Necessary. itmp -> itemp, had to spell differently for now as
      ! itemp lives in param.h still (avoid clash). Same for islt -> isalt
      integer :: itmp=0, islt=0
      ! Remaining passive/bgc tracer indices declared in tracers_defs_idx.h
#include "tracers_defs_idx.h"

      ! Interpolate forcing from coarser input grid (=1) or not (=0)
      integer :: interp_t = 0 ! Factor 2 scaling only for now

      ! End of user inputs
      ! ****************************************************************

      !-- Array of tracers
      real, public :: t(GLOBAL_2D_ARRAY,N,3,NT)

      !-- Tracer variables (user input required in init_tracers)
      character*42, dimension(nt), public :: t_vname ! Short name
      character*60, dimension(nt), public :: t_lname ! long name
      character*42, dimension(nt), public :: t_units ! Tracer units
      character*47, dimension(nt), public :: t_tname ! Tracer input time variable name
      logical, dimension(nt), public      :: wrt_t   ! T/F output tracer.
      integer, dimension(nt) :: t_ana_frc ! Whether surf flux is read in (0), or analytical (1)

      !-- Tracer netcdf variables as arrays/matrices of 'NT' length
      ! Abstracted from user.
      ! Final tracer concentrations live in 't' in ocean3d.h
      ! Surface tracer flux lives in 'stflx' in surf_flux.F module.
      real                   :: t_data(GLOBAL_2D_ARRAY,2,NT) ! Could set to zero here?
      integer, dimension(NT) :: t_grd_type  = 0
      integer, dimension(NT) :: t_file_indx = 0
      integer, dimension(NT) :: t_irec      = 0
      integer, dimension(NT) :: t_it1 = 1, t_it2 = 2
      real, dimension(2,NT)  :: t_times = -99 ! [-99,-99]

      integer, public :: iTandS ! Combined index of temperature and salinity
                                ! public as used in set_forces.F. Need to be parameter
      integer :: pt_tot ! total number of passive tracers
      integer :: itot=0 ! index counter to count total number of tracers.
      ! 3 below currently in param.h as per old code
!      integer :: itrc_bio=0 ! Starting tracer index for bgc tracers. Set in tracers_defs.h
      ! Number of 'base' bgc tracers, i.e. minimum number of bgc tracers without cocco
      ! and nitrogen cycle additions
!      integer :: ntrc_bio_base
!      integer :: ntrc_bio !  - from (param.h - ntrc_bio=ntrc_bio_base+ntrc_bio_cocco+ntrc_bio_ncycle)


      public set_surf_tracer_flx
      public init_tracers
      public set_bry_tracer_defaults

      contains
! ----------------------------------------------------------------------

      subroutine set_surf_tracer_flx(istr,iend,jstr,jend) ! set tracer flux at surface

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend

      ! local
      integer itrc ! Tracer number for loop index
      character(len=46) :: t_flx_name ! Tracer time name

      ! Need ocean3d.h for 't' tracer variables:
#include "ocean3d.h"

      do itrc=iTandS+1,NT ! Instead of 3 use iTandS variable DevinD todo

        if(t_ana_frc(itrc)==0) then ! Read in forcing data (not analytical)

          ! For now set up variable _flx name here but might be better to
          ! create in init_tracers function where this is done once and
          ! storedin an array for speed.

          t_flx_name = trim(t_vname(itrc)) / / '_flx' ! Forcing file flux name

          call set_frc_var_tile(
     &        istr,iend,          jstr,jend,        ! Tile bounds indices
     &        t_flx_name,         t_tname(itrc),    ! Text name of array and time
     &        t_data(:,:,:,itrc), stflx(:,:,itrc),  ! Variable arrays
     &        t_times(:,itrc),    t_grd_type(itrc), ! Input times & Grid type
     &        t_file_indx(itrc),  t_irec(itrc),     ! File indx & Current input record
     &        t_it1(itrc),        t_it2(itrc),      ! Time index placeholders
     &        interp_t )                            ! Online spatial interpolation flag

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

      end subroutine set_surf_tracer_flx
! ----------------------------------------------------------------------

      subroutine set_ana_surf_tracer_flx(istr,iend,jstr,jend,itrc)
      ! Set analytical surface tracer flux
      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      integer itrc ! Current tracer index number

      ! local
      integer i,j

! Needed for 'FIRST_TIME_STEP'
#include "scalars.h"

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

      end subroutine set_ana_surf_tracer_flx
! ----------------------------------------------------------------------

      subroutine init_tracers()

      ! Set all tracer variable values. This was placed at bottom of
      ! the module for ease of reading the rest of the module's code.

      implicit none

      ! save
      ! SHOULD DECLARE TRACER INDICES IN TRACERS_DEFS.H & USE THE SAVE COMMAND HERE

      ! It is necessary to keep track of tracer index number in order to
      ! correctly locate and calculate upon the tracer in 't' the array
      ! of all the tracers.

      ! Attempted to atleast here keep all variables together.
      ! Yes, still need to declare them in module preamble, however,
      ! they take no value there so order/value is not a worry.
      ! Old code: index is set in param.h, and variable values in
      ! init_scalars.F.
      ! New method, both set here.

! Following is just to test using indxPO4, should delete after confirming not needed.
#ifdef BIOLOGY_BEC2
# include "ncvars.h"
#endif

      ! Core tracers - temp and salt:
      itmp=1;                      wrt_t(itmp) =.True.
      t_vname(itmp)='temp';        t_units(itmp)='Celsius'
      t_lname(itmp)='potential temperature'
      iTandS = 1 ! If only temp, no salt.
#ifdef SALINITY
      islt=2;                      wrt_t(islt) =.True.
      t_vname(islt)='salt';        t_units(islt)='PSU'
      t_lname(islt)='salinity'
      iTandS = 2 ! If both temp and salt.
#endif
      itot=iTandS ! Set total number of tracers so far as iTandS

      ! Additional passive tracers:
#include "tracers_defs.h"

      end subroutine init_tracers
! ----------------------------------------------------------------------

      subroutine set_bry_tracer_defaults(t_bry_dt, itrc, itbry, LMmpp)
      ! Sets default values for tracers at the boundary if missing
      ! in the _bry.nc file.
      ! In the long run this should be put into an include file to avoid
      ! conflicts when module is updated.
      ! Taken from Pierre's code get_bry_all.F

      implicit none
#  include "param.h"
      real t_bry_dt(0:LMmpp,N,2,NT), cff,cff1
      integer, intent(in):: LMmpp, itbry,itrc
      integer :: ij, k
#  define temp cff
#  define SiO4 cff1

      if (itrc == itemp) then
        t_bry_dt(:,:,itbry,itrc) = 5.0 ! default temperature
#  ifdef SALINITY
      else if (itrc == isalt) then
        t_bry_dt(:,:,itbry,itrc) = 35.0 ! default salinity
#  endif
#  ifdef BIOLOGY_NPZDOC
      elseif (itrc == iNO3_) then
        cff=20./3.
        cff1=2./3.
        do k=1,N
          do ij=0,LMmpp
            temp=t_bry_dt(j,k,itbry,itemp)
              if (temp.lt.8.) then
                SiO4=30.
              elseif (temp.le.11.) then
                SiO4=30.-((temp-8.)*cff)
              elseif (temp.le.13.) then
                SiO4=10.-((temp-11.)*4.)
              elseif (temp.le.16.) then
                SiO4=2.-((temp-13.)*cff1)
              else
                SiO4=0.
              endif
              t_bry_dt(j,k,itbry,itrc)=1.67+0.5873*SiO4
     &             +0.0144*SiO4**2
     &             +0.0003099*SiO4**3
#   ifdef NO3_SFC_CUTOFF
              t_bry_dt(ij,k,itbry,itrc)=t_bry_dt(ij,k,itbry,itrc)*
     &             (1-.5*(tanh((z_r(iwest,ij,k)+100.)/8.)+1))
#   endif /* NO3_SFC_CUTOFF */
          enddo
        enddo
      elseif (itrc .eq. iNH4_) then ! DevinD why is there an underscore?? Flag not used anyway
        t_bry_dt(:,:,itbry,itrc) = 0.0
      elseif (itrc .eq. iChla) then
        t_bry_dt(:,:,itbry,itrc) = 0.10
      elseif (itrc .eq. iPhyt) then
        t_bry_dt(:,:,itbry,itrc) = 0.05
      elseif (itrc .eq. iZoo_) then ! DevinD why is there an underscore?? Flag not used anyway
        t_bry_dt(:,:,itbry,itrc) = 0.02
      elseif (itrc .eq. iSDet) then
        t_bry_dt(:,:,itbry,itrc) = 0.0
      elseif (itrc .eq. iLDet) then
        t_bry_dt(:,:,itbry,itrc) = 0.0
#   ifdef OXYGEN
      elseif (itrc .eq. iO2) then
        t_bry_dt(:,:,itbry,itrc) = 177.9  ! mmol/m3
#    ifdef CARBON
      elseif (itrc .eq. iDIC) then
        t_bry_dt(:,:,itbry,itrc) = 2200.
      elseif (itrc .eq. iTALK) then
        t_bry_dt(:,:,itbry,itrc) = 2300.
      elseif (itrc .eq. iSDetC) then
        t_bry_dt(:,:,itbry,itrc) = 0.0
      elseif (itrc .eq. iLDetC) then
        t_bry_dt(:,:,itbry,itrc) = 0.0
      elseif (itrc .eq. iCaCO3) then
        t_bry_dt(:,:,itbry,itrc) = 0.0
#    endif
#   endif /* OXYGEN */
#  endif /* BIOLOGY_NPZDOC */

#    if defined BIOLOGY_BEC2
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


      end subroutine set_bry_tracer_defaults
! ----------------------------------------------------------------------

      end module tracers

