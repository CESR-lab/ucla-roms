#include "cppdefs.opt"

      subroutine init_scalars (ierr)

! Set initial values for globally accessable (stored in common blocks)
! scalar variables.  Typically these are associated with time indices,
! global summation counters, model utility switches, i/o contlol, etc.
! Note that NONE of the model physical parameters is initialized here.

      use param
      use bgc ! init_scalars_bec2
      use boundary
      use comm_vars
      use mixing
      use ncvars
      use scalars
      use roms_read_write
      use dimensions

      use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
      use, intrinsic :: iso_fortran_env, only: real32

      implicit none
      integer ierr, i

      real :: nan  !real(real32) :: nan
      nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
      init =  nan ! 0. ! nan  ! for checking no bad inits

#define ALL_DATA

c      integer size, kmp_get_stacksize_s
c      size=kmp_get_stacksize_s()
c      write(*,*) 'kmp_stacksize =', size

      start_time=0. ; iic=0          ! Set time and evolving time
      nfast=1       ; kstp=1         ! stepping indices to their
      ntstart=1     ; knew=1         ! initial values needed for

#ifdef SOLVE3D
      iif=1   ; nstp=1               ! proper logical initialization
      nrhs=1  ; nnew=1               ! of model fields.
#endif
                                     ! Reset intra-thread
      may_day_flag=0                 ! communication signal
      trd_count=0                    ! variables and counters.
      cpu_all=0.   !<-- array
      barr_count=0 !<-- array

      mcheck=0
      first_time=0
      tile_count=0                   ! Reset global sumation variables
      bc_count=0                     ! and counters needed for volume
      volume=0.                      ! and energy diagnostics, and for
c      avke=0.                       ! volume conservation constraint.
c      avpe=0.
c      prev_ke=0.
c      avke2b=0.
c      avke3bc=0.
c      v2d_max=0.
c      v3d_max=0.
c      v3bc_max=0.

      hmin=+1.E+20   ; hmax=-1.E+20     ! Set extreme values to large
      grdmin=+1.E+20 ; grdmax=-1.E+20   ! unrealistic numbers, positive
      Cg_min=+1.E+20 ; Cg_max=-1.E+20   ! or negative, depending on the
#ifdef SOLVE3D
      rx0=-1.E+20    ; rx1=-1.E+20      ! context.
#endif
      bc_crss=QuadZero
                             ! Set slipperness parameter to "free slip"
      gamma2=1.              ! along the straight boundaries: this may
                             ! be overwritten later by "read_inp".

      do i=1,max_frc_files
        ncfrc(i)=-1
      enddo
#ifndef ANA_BRY
       bry_id=-1
#endif

      call get_date (date_str)

      call init_dimensions                                 ! needs to come before init_output_indices
      call init_output_indices

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

#ifdef BIOLOGY_BEC2
        call init_scalars_bec2
#endif

      ierr=0
      call read_inp (ierr)
      if (ierr==0) then
#ifdef SOLVE3D
        call set_weights
#endif
#ifdef GRID_LEVEL
# if GRID_LEVEL < MAX_GRID_LEVEL
#  if GRID_LEVEL == 1
        call init_scalars_2 (ierr)
#  elif GRID_LEVEL == 2
        call init_scalars_3 (ierr)
#  elif GRID_LEVEL == 3
        call init_scalars_4 (ierr)
#  endif
# endif
#endif
      endif
      end
