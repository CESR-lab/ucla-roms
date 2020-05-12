#include "cppdefs.h"


      subroutine check_switches1(ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new CPP-switches  !!!!!!
!!!!!! were introduced into "cppdefs.h".  NO ACTION IS NEEDED  !!!!!!
!!!!!! if changes in "cppdefs.h" are limited to activation or  !!!!!!
!!!!!! deactivation of previously known switches.              !!!!!!
!!!!!! To refresh this file compile and execute "cppcheck.F"   !!!!!!
!!!!!! as an independent program, or use commands              !!!!!!
!!!!!! "make checkdefs" or "make depend".                      !!!!!!
!!!!!! Number of Configuration Choices:                      1 !!!!!!
!!!!!! Total number of CPP-switches:                        99 !!!!!!

      implicit none
      integer ierr, is,ie, iexample
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      do is=1,max_opt_size
        cpps(is:is)=' '
      enddo
      iexample=0
!
! cppdefs.h                                                       
!
      is=1 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='<cppdefs.h>'
      cpps(ie+1:ie+1)=' '

#ifdef WEC
      iexample=iexample+1
      is=ie+2 ; ie=is+2
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='WEC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DOUBLE_GYRE
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DOUBLE_GYRE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SOLVE3D
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SOLVE3D'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_ADV
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='UV_ADV'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_COR
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='UV_COR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SALINITY
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SALINITY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef NONLIN_EOS
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='NONLIN_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_VIS2
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='UV_VIS2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef TS_DIF2
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='TS_DIF2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef TS_DIF4
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='TS_DIF4'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_GRID
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_GRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_INITIAL
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_INITIAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SMFLUX
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_SMFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_STFLUX
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_STFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_VMIX
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_VMIX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CANYON_A
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='CANYON_A'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef EW_PERIODIC
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='EW_PERIODIC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CANYON_B
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='CANYON_B'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DAMEE_B
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DAMEE_B'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef TCLIMATOLOGY
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='TCLIMATOLOGY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef TNUDGING
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='TNUDGING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef QCORRECTION
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='QCORRECTION'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SFLX_CORR
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SFLX_CORR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef VIS_GRID
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='VIS_GRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DIF_GRID
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DIF_GRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_MIXING
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_MIXING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_RIMIX
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_RIMIX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_CONVEC
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_CONVEC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CURVGRID
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='CURVGRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPHERICAL
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPHERICAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef MASKING
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='MASKING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef EASTERN_WALL
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='EASTERN_WALL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef WESTERN_WALL
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='WESTERN_WALL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SOUTHERN_WALL
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SOUTHERN_WALL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef NORTHERN_WALL
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='NORTHERN_WALL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef REST_STATE_TEST
      is=ie+2 ; ie=is+14
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='REST_STATE_TEST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SSFLUX
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_SSFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SST
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_SST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SRFLUX
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_SRFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef EKMAN_SPIRAL
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='EKMAN_SPIRAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef NS_PERIODIC
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='NS_PERIODIC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef GRAV_ADJ
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='GRAV_ADJ'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef COLD_FILAMENT
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='COLD_FILAMENT'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_WEST
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_WEST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_EAST
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_EAST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M2FLATHER
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M2FLATHER'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M3ORLANSKI
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M3ORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_TORLANSKI
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_TORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_BRY
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef T_FRC_BRY
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='T_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef Z_FRC_BRY
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='Z_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef M2_FRC_BRY
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='M2_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef M3_FRC_BRY
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='M3_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BLACK_SEA
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BLACK_SEA'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BULK_FLUX
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BULK_FLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BULK_SM_UPDATE
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BULK_SM_UPDATE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef WIND_AT_RHO_POINTS
      is=ie+2 ; ie=is+17
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='WIND_AT_RHO_POINTS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef WIND_MAGN
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='WIND_MAGN'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_KPP
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_KPP'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_BKPP
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_BKPP'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_NONLOCAL
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_NONLOCAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_DDMIX
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_DDMIX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ATLANTIC
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ATLANTIC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DIURNAL_SRFLUX
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DIURNAL_SRFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_SOUTH
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_SOUTH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPONGE
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPONGE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SEAMOUNT100
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SEAMOUNT100'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef TS_HADV_UP3
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='TS_HADV_UP3'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_NORTH
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_NORTH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M3SPECIFIED
      is=ie+2 ; ie=is+14
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M3SPECIFIED'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_TSPECIFIED
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_TSPECIFIED'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef FRC_BRY
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef PACIFIC
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='PACIFIC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SEA_ICE_NOFLUX
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SEA_ICE_NOFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef M2NUDGING
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='M2NUDGING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef M3NUDGING
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='M3NUDGING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef PACIFIC_2D
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='PACIFIC_2D'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OVERFLOW
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OVERFLOW'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SEAMOUNT
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SEAMOUNT'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ISWAKE
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ISWAKE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LINEAR_DRAG_ONLY
      is=ie+2 ; ie=is+15
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LINEAR_DRAG_ONLY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SOLITON
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SOLITON'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef RIVER
      is=ie+2 ; ie=is+4
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='RIVER'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef PSOURCE
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='PSOURCE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_PSOURCE
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_PSOURCE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M2ORLANSKI
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M2ORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UPWELLING
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='UPWELLING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CANBAS2
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='CANBAS2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef AVERAGES
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='AVERAGES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPLASH
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPLASH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef RIVER_SOURCE
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='RIVER_SOURCE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BRY
      is=ie+2 ; ie=is+2
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M2SPECIFIED
      is=ie+2 ; ie=is+14
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M2SPECIFIED'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_WEC_FRC
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_WEC_FRC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef VERBOSE
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='VERBOSE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPLIT_EOS
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPLIT_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef MRL_WCI
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='MRL_WCI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_WWAVE
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_WWAVE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BRK0
      is=ie+2 ; ie=is+3
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BRK0'
      cpps(ie+1:ie+1)=' '
#endif
      if (iexample == 0) then
        mpi_master_only write(*,'(/1x,A)')
     & '### ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample > 1) then
        mpi_master_only write(*,'(/1x,2A/)')
     &  '### ERROR: More than one main configuration ',
     &             'switch is defined in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99  mpi_master_only write(*,'(/1x,2A/14x,A)')
     &  '### ERROR: checkdefs :: Insufficient size of string "cpps" ',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
