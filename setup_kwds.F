#include "cppdefs.h"

      subroutine setup_kwds (ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new keywords were !!!!!!
!!!!!! introduced into "read_inp.F". To create or refresh this !!!!!!
!!!!!! file use compile and execute "checkkwds.F" as an        !!!!!!
!!!!!! independent program, or use commands "make checkkwds"   !!!!!!
!!!!!! or "make depend".                                       !!!!!!

      implicit none
      integer ierr, is,ie
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      do is=1,max_opt_size
        kwds(is:is)=' '
      enddo
      is=1
#ifdef WEC
#endif
#ifndef NONLIN_EOS
#endif
#ifdef MPI
#endif
#ifdef SOLVE3D
#endif
#ifdef INT_IARGC
#endif
#ifdef SOLVE3D
# ifdef LMD_KPP
# endif
#endif
#if defined SOLITON
#elif defined SEAMOUNT
#elif defined WAVE_RAD
#elif defined ISWAKE
#elif defined UPWELLING
#elif defined GRAV_ADJ
#elif defined IBERIA
#elif defined RIVER
#elif defined PACIFIC
#elif defined PACIFIC_2D
#elif defined CANBAS2
#elif defined SPIRAL
#elif defined PACBIG
#elif defined NWPAC
#elif defined SMSC
#elif defined SMSC_SUB
#elif defined USWEST
#elif defined SCB_TEST
#elif defined DAMEE_B
#else
#endif
#ifdef MPI
#endif
#ifdef MPI
#endif
#ifdef AVERAGES
#endif
      ie=is + 5
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='title'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +13
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='time_stepping'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifndef SOLVE3D
#endif
#ifdef SOLVE3D
      ie=is + 7
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='S-coord'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='rho0'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifndef NONLIN_EOS
      ie=is +11
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='lin_rho_eos'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  ifdef SALINITY
#  endif
#  ifdef SALINITY
#  endif
# endif
#endif /* SOLVE3D */
#ifdef UV_VIS2
      ie=is +12
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='lateral_visc'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if !defined EW_PERIODIC || !defined NS_PERIODIC
      ie=is + 6
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='gamma2'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      ie=is +12
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='tracer_diff2'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif
      ie=is +11
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='bottom_drag'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined PP_MIXING\
                         && !defined MY2_MIXING && !defined MY25_MIXING
      ie=is +15
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='vertical_mixing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif
# ifdef MY25_MIXING
      ie=is +13
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='MY_bak_mixing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  ifdef Q_DIF2
#  endif
#  ifdef Q_DIF4
#  endif
# endif
# if defined SFLX_CORR && defined SALINITY
      ie=is +14
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='SSS_correction'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif /* SOLVE3D */
#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined TNUDGING \
  || defined Z_FRC_BRY || defined M3_FRC_BRY || defined M2NUDGING \
                       || defined M3NUDGING  || defined WKB_FRC_BRY
      ie=is + 5
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='ubind'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
# endif
#endif
#ifdef SPONGE
      ie=is + 8
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='v_sponge'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifndef ANA_GRID
      ie=is + 4
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='grid'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif
      ie=is + 7
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='initial'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifdef ANA_INITIAL
#endif
#if defined MPI && defined PARALLEL_FILES
#endif
#ifdef ANA_INITIAL
#endif
#if !defined ANA_SMFLUX  ||  defined SOLVE3D  &&\
   ( !defined ANA_STFLUX  ||\
     (defined SG_BBL96    && !defined ANA_BSEDIM) ||\
     (defined SG_BBL96    && !defined ANA_WWAVE)  ||\
     (defined QCORRECTION && !defined ANA_SST)    ||\
     (defined SALINITY    && !defined ANA_SSFLUX) ||\
     (defined LMD_KPP     && !defined ANA_SRFLUX))
      ie=is + 7
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='forcing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef MPI_SILENT_MODE
# endif
# ifdef MPI_SILENT_MODE
# endif
#endif
#if defined T_FRC_BRY || defined M2_FRC_BRY || defined M3_FRC_BRY \
                      || defined Z_FRC_BRY  || defined WKB_FRC_BRY
# ifndef ANA_BRY
      ie=is + 8
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='boundary'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
#  ifdef MPI_SILENT_MODE
#  endif
#  ifdef MPI_SILENT_MODE
#  endif
# endif
#endif
#if defined SOLVE3D &&( defined TCLIMATOLOGY ||( defined TNUDGING && \
                        defined T_FRC_BRY )) && !defined ANA_TCLIMA
      ie=is +11
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='climatology'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif
      ie=is + 7
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='restart'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#if defined MPI && defined PARALLEL_FILES
#endif
#ifndef SOLVE3D
#endif
      ie=is + 7
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='history'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#if defined MPI && defined PARALLEL_FILES
#endif
#ifndef SOLVE3D
#endif
      ie=is +22
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='primary_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
      ie=is +24
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='auxiliary_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef SALINITY
# endif
# ifdef LMD_KPP
# endif
# ifdef LMD_BKPP
# endif
# ifdef SALINITY
# endif
# ifdef LMD_KPP
# endif
# ifdef LMD_BKPP
# endif
# ifdef SALINITY
# endif
# ifdef LMD_KPP
# endif
# ifdef LMD_BKPP
# endif
#endif /* SOLVE3D */
#ifdef AVERAGES
      ie=is + 8
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifndef SOLVE3D
# endif
      ie=is +16
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='primary_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
      ie=is +18
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='auxiliary_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef SALINITY
#  endif
#  ifdef LMD_KPP
#  endif
#  ifdef LMD_BKPP
#  endif
# endif /* SOLVE3D */
#endif /* AVERAGES */
#ifdef WEC
      ie=is +18
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='wci_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +18
      if (ie>=max_opt_size) goto 99
      kwds(is:ie)='wci_average_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
# endif /* AVERAGES */
#endif /* WEC */
#ifdef MPI
#endif
      return
  99  mpi_master_only write(*,'(/1x,2A/25x,A/)')
     &  '### ERROR: setup_kwds :: Insufficient size of str "kwds" ',
     &  'in file "strings.h".',  'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
