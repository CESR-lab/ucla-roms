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

      implicit none
      integer ierr, i
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

      ncrst=-1   ; nrecrst=0       ! Set netCDF file IDs to negative
      nchis=-1   ; nrechis=0       ! values to indicate their initially
#ifdef AVERAGES
      ncavg=-1   ; nrecavg=0       ! closed status. Also reset record
#endif
      do i=1,max_frc_files
        ncfrc(i)=-1
      enddo
      ncidclm=-1
#ifndef ANA_BRY
       bry_id=-1
#endif

      call get_date (date_str)

! Names of variables in NetCDF output files. The first element
! is the name of the variable; the other two are are attributes.

c--#define CR
CR      write(*,'(1x,A,I3)') 'indxTime =',indxTime
      vname(1,indxTime)='ocean_time'
      vname(2,indxTime)='Time since initialization'
      vname(3,indxTime)='second'

CR      write(*,'(1x,A,I3)') 'indxZ =',indxZ
      vname(1,indxZ)='zeta'
      vname(2,indxZ)='free-surface elevation'
      vname(3,indxZ)='meter'

CR      write(*,'(1x,A,I3)') 'indxUb =',indxUb
      vname(1,indxUb)='ubar'
      vname(2,indxUb)='vertically averaged u-momentum component'
      vname(3,indxUb)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxVb = ',indxVb
      vname(1,indxVb)='vbar'
      vname(2,indxVb)='vertically averaged v-momentum component'
      vname(3,indxVb)='meter second-1'

#ifdef SOLVE3D
CR      write(*,'(1x,A,I3)') 'indxU = ', indxU
      vname(1,indxU)='u'
      vname(2,indxU)='u-momentum component'
      vname(3,indxU)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxV = ', indxV
      vname(1,indxV)='v'
      vname(2,indxV)='v-momentum component'
      vname(3,indxV)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxO =',indxO
      vname(1,indxO)='omega'
      vname(2,indxO)='S-coordinate vertical velocity'
      vname(3,indxO)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxW =', indxW
      vname(1,indxW)='w'
      vname(2,indxW)='vertical velocity'
      vname(3,indxW)='meter second-1'

CR      write(*,'(1x,A,I3)') 'indxR =', indxR
      vname(1,indxR)='rho'
      vname(2,indxR)='density anomaly'
      vname(3,indxR)='kilogram meter-3'

CR      write(*,'(1x,A,I3)') 'indxT = ', indxT
      vname(1,indxT)='temp'
      vname(2,indxT)='potential temperature'
      vname(3,indxT)='Celsius'

# ifdef SALINITY
CR      write(*,'(1x,A,I3)') 'indxS = ', indxS
      vname(1,indxS)='salt'
      vname(2,indxS)='salinity'
      vname(3,indxS)='PSU'
# endif
!# ifdef BIOLOGY
!CR      write(*,'(1x,A,I3)') 'indxNO3 =', indxNO3
!      vname(1,indxNO3)='NO3'
!      vname(2,indxNO3)='NO3 Nutrient'
!      vname(3,indxNO3)='mMol N m-3'
!
!CR      write(*,'(1x,A,I3)') 'indxNH4 =', indxNH4
!      vname(1,indxNH4='NH4'
!      vname(2,indxNH4)='NH4 Nutrient'
!      vname(3,indxNH4)='mMol N m-3'
!
!CR      write(*,'(1x,A,I3)') 'indxDet =', indxDet
!      vname(1,indxDet)='DET'
!      vname(2,indxDet)='Detritus'
!      vname(3,indxDet)='mMol N m-3'
!
!CR      write(*,'(1x,A,I3)') 'indxPhyt =', indxPhyt
!      vname(1,indxPhyt)='PHYTO'
!      vname(2,indxPhyt)='Phytoplankton'
!      vname(3,indxPhyt)='mMol N m-3'
!
!CR      write(*,'(1x,A,I3)') 'indxZoo =', indxZoo
!      vname(1,indxZoo)='ZOO'
!      vname(2,indxZoo)='Zooplankton'
!      vname(3,indxZoo)='mMol N m-3'
!# endif

CR      write(*,'(1x,A,I3)') 'indxAkv =', indxAkv
      vname(1,indxAkv)='Akv'
      vname(2,indxAkv)='vertical viscosity coefficient'
      vname(3,indxAkv)='meter2 second-1'

CR      write(*,'(1x,A,I3)') 'indxAkt =',indxAkt
      vname(1,indxAkt)='Akt'
      vname(2,indxAkt)='vertical thermal conductivity coefficient'
      vname(3,indxAkt)='meter2 second-1'
# ifdef SALINITY
CR      write(*,'(1x,A,I3)') 'indxAkt =',indxAks
      vname(1,indxAks)='AKs'
      vname(2,indxAks)='salinity vertical diffusion coefficient'
      vname(3,indxAks)='meter2 second-1'
# endif
# ifdef LMD_KPP
CR      write(*,'(1x,A,I3)') 'indxHbls =',indxHbls
      vname(1,indxHbls)='hbls'
      vname(2,indxHbls)='Thickness of KPP surface boundary layer'
      vname(3,indxHbls)='meter'
# endif
# ifdef LMD_BKPP
CR      write(*,'(1x,A,I3)') 'indxHbbl =',indxHbbl
      vname(1,indxHbbl)='hbbl'
      vname(2,indxHbbl)='Thickness of KPP bottom boundary layer'
      vname(3,indxHbbl)='meter'
# endif
#endif

#ifdef SOLVE3D
# ifdef SG_BBL96
#  ifndef ANA_WWAVE
      vname(1,indxWWA)='Awave'
      vname(2,indxWWA)='wind induced wave amplitude'
      vname(3,indxWWA)='meter'

      vname(1,indxWWD)='Dwave'
      vname(2,indxWWD)='wind induced wave direction'
      vname(3,indxWWD)='degrees'

      vname(1,indxWWP)='Pwave'
      vname(2,indxWWP)='wind induced wave Period'
      vname(3,indxWWP)='second'
#  endif
#  ifndef ANA_BSEDIM
      vname(1,indxBSS)='Ssize'
      vname(2,indxBSS)='bottom sediment grain diameter size'
      vname(3,indxBSS)='meter'

      vname(1,indxBSD)='Sdens'
      vname(2,indxBSD)='bottom sediment grain density'
      vname(3,indxBSD)='kilogram meter-3'
#  endif
# endif
#endif
#ifdef ICE
      vname(1,indxAi)='aice'
      vname(2,indxAi)='fraction of cell covered by ice'
      vname(3,indxAi)='nondimensional'

      vname(1,indxUi)='uice'
      vname(2,indxUi)='u-component of ice velocity'
      vname(3,indxUi)='meter sec-1'

      vname(1,indxVi)='vice'
      vname(2,indxVi)='v-component of ice velocity'
      vname(3,indxVi)='meter sec-1'

      vname(1,indxHi)='hice'
      vname(2,indxHi)='depth of ice cover'
      vname(3,indxHi)='meter'

      vname(1,indxHS)='snow_thick'
      vname(2,indxHS)='depth of snow cover'
      vname(3,indxHS)='meter'

      vname(1,indxTIsrf)='tisrf'
      vname(2,indxTIsrf)='temperature of ice surface'
      vname(3,indxTIsrf)='Kelvin'
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
