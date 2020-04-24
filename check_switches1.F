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
!!!!!! Total number of CPP-switches:                        20 !!!!!!

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
#ifdef ANA_SRFLUX
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_SRFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_STFLUX
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_STFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SST
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_SST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SALINITY
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SALINITY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef ANA_SSFLUX
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='ANA_SSFLUX'
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
#ifdef NONLIN_EOS
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='NONLIN_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPLIT_EOS
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPLIT_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef EW_PERIODIC
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='EW_PERIODIC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef NS_PERIODIC
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='NS_PERIODIC'
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
