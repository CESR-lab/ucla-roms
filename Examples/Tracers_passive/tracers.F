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
      ! Set to zero as possible means to avoid errors later? Might not be
      ! Necessary. itmp -> itemp, had to spell differently for now as
      ! itemp lives in param.h still (avodi clash). Same for islt -> isalt
      integer :: itmp=0, islt=0, itrace1=0, isalt2=0 ! , iptrace2=0

      ! Interpolate forcing from coarser input grid (=1) or not (=0)
      integer :: interp_t = 0 ! Factor 2 scaling only for now

      ! End of user inputs
      ! ****************************************************************

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
      real                   :: t_data(GLOBAL_2D_ARRAY,2,NT)
      integer, dimension(NT) :: t_grd_type  = 0
      integer, dimension(NT) :: t_file_indx = 0
      integer, dimension(NT) :: t_irec      = 0
      integer, dimension(NT) :: t_it1 = 1, t_it2 = 2
      real, dimension(2,NT)  :: t_times = -99 ! [-99,-99]

      integer, public :: iTandS ! Combined index of temperature and salinity
                                ! public as used in set_forces.F


      public set_surf_tracer_flx
      public init_tracers

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
            stflx(i,j,itrc)= 0.5 ! 0.
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

      ! It is necessary to keep track of tracer index number in order to
      ! correctly locate and calculate upon the tracer in 't' the array
      ! of all the tracers.

      ! Attempted to atleast here keep all variables together.
      ! Yes, still need to declare them in module preamble, however,
      ! they take no value there so order/value is not a worry.
      ! Old code: index is set in param.h, and variable values in
      ! init_scalars.F.
      ! New method, both set here.

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

      ! Additional passive tracers:
#include "tracers_defs.h"

      end subroutine init_tracers
! ----------------------------------------------------------------------

      end module tracers

