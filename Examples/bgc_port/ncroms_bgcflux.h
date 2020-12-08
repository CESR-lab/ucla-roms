! This header file contains all variables and parameters for the 
! netcdf output of biogeochemical fluxes.
!
! NOTE: This file must always be included AFTER bgcflux.h as it uses
! some parameters that are defined in that file.

#if defined SOLVE3D && defined BGC_FLUX_ANALYSIS

      logical new_bgc_flux_his
      integer n_bgc_flux_his, nrpf_bgc_flux_his
     &     , ncid_bgc_flux_his, nrec_bgc_flux_his
      common /scalars_bgc/
     &     new_bgc_flux_his
     &     , n_bgc_flux_his, nrpf_bgc_flux_his
     &     , ncid_bgc_flux_his, nrec_bgc_flux_his

#ifdef AVERAGES
      integer ncid_bgc_flux_avg, nrec_bgc_flux_avg
     &     , nrpf_bgc_flux_avg
      common /scalars_bgc_avg/ ncid_bgc_flux_avg, nrec_bgc_flux_avg
     &     , nrpf_bgc_flux_avg

      character*80 bgc_flux_avg_name
      common /c_bgcflux_avg/ bgc_flux_avg_name

#endif

# ifdef BIOLOGY_NPZDOC
      integer, parameter :: num_bgcflux = 20
     &     + NumFluxTerms + NumVSinkTerms
!DL: gas exchange fluxes are always computed and output:
!#   ifdef OXYGEN
!     &     + NumGasExcTerms 
!#   endif /* OXYGEN */
#   ifdef SEDIMENT_BIOLOGY
     &     + NumSedFluxTerms
#   endif /* SEDIMENT_BIOLOGY */
# elif defined BIOLOGY_BEC
      integer, parameter :: num_bgcflux_2d = 12
      integer, parameter :: num_bgcflux = 74 
      integer, dimension(num_bgcflux) :: vid_bec_flux_his
      common /c_bgcflux_bec/ vid_bec_flux_his
# endif

      character*80 bgc_flux_his_name,
     &     vname_bgcflux(4, num_bgcflux)
      common /c_bgcflux/ bgc_flux_his_name, vname_bgcflux

!!!!!!!!!!!!!!!!!!!!!!!! NPZD model !!!!!!!!!!!!!!!!!!!
# ifdef BIOLOGY_NPZDOC
! indices for the netcdf output
# ifdef OXYGEN
      integer, parameter :: indxU10 = 1
      integer, parameter :: indxKvO2 = 2
      integer, parameter :: indxO2sat = 3
#   ifdef CARBON
      integer, parameter :: indxKvCO2 = 4
      integer, parameter :: indxCO2sol = 5
      integer, parameter :: indxPCO2 = 6
      integer, parameter :: indxPCO2air = 7
      integer, parameter :: indxPH = 8
#   endif /* CARBON */
# endif /* OXYGEN */

! if OXYGEN and/or CARBON are not defined, there will be a gap in indices
      integer, parameter :: indxPAR = 9
      integer, parameter :: indxPARinc = 10

! deliberately leave space for additional non-flux variables
      integer, parameter :: indxFlux = 20 
! first vertical sinking flux
      integer, parameter :: indxVSinkFlux = indxFlux + NumFluxTerms 
#  ifdef OXYGEN
! first gas exchange flux
!      integer, parameter :: indxGasExcFlux = 
!     &     indxFlux + NumFluxTerms + NumVSinkTerms 
#  endif
#  ifdef SEDIMENT_BIOLOGY
! first sediment-related flux
      integer, parameter :: indxSedFlux = indxFlux + NumFluxTerms + 
     &     NumVSinkTerms
#  endif /* SEDIMENT_BIOLOGY */

      integer hisPAR_flux, hisPARinc_flux
     &     , hisFlux(NumFluxTerms)
     &     , hisVSinkFlux(NumVSinkTerms)
#  ifdef OXYGEN
     &     , hisONNO3, hisONNH4
     &     , hisU10, hisKvO2, hisO2sat
#  ifdef OXYLIM
     &     , hisNCDET
#  endif /* OXYLIM */
#   ifdef CARBON
     &     , hisCNP, hisCNZ, hisrCaCO3orgC
     &     , hisKvCO2, hisCO2sol
#   endif /* CARBON */
#  endif /* OXYGEN */
!#   ifdef OXYGEN
!     &     , hisGasExcFlux(NumGasExcTerms) 
!#   endif /* OXYGEN */
#   ifdef SEDIMENT_BIOLOGY
     &     , hisSedFlux(NumSedFluxTerms)
#   endif /* SEDIMENT_BIOLOGY */
     &     , bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ

      common /ncids_bgc_flux/ hisPAR_flux, hisPARinc_flux
     &     , hisFlux
     &     , hisVSinkFlux
#  ifdef OXYGEN
     &     , hisONNO3, hisONNH4
     &     , hisU10, hisKvO2, hisO2sat
#  ifdef OXYLIM
     &     , hisNCDET
#  endif /* OXYLIM */
#   ifdef CARBON
     &     , hisCNP, hisCNZ, hisrCaCO3orgC
     &     , hisKvCO2, hisCO2sol
#   endif /* CARBON */
#  endif /* OXYGEN */
!#   ifdef OXYGEN
!     &     , hisGasExcFlux
!#   endif /* OXYGEN */
#   ifdef SEDIMENT_BIOLOGY
     &     , hisSedFlux
#   endif /* SEDIMENT_BIOLOGY */
     &     , bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ
#endif /* BIOLOGY_NPZDOC */      
#ifdef BIOLOGY_BEC
      integer bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ
      common /ncids_bgc_flux/ bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ      
#endif /* BIOLOGY_BEC */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ifdef AVERAGES

      logical new_bgc_flux_avg
      integer nts_bgc_flux_avg, n_bgc_flux_avg
      common /scalars_bgc_avg/
     &     new_bgc_flux_avg,
     &     nts_bgc_flux_avg, n_bgc_flux_avg
      real time_bgc_flux_avg
      common /scalars_bgc_avg_real/ time_bgc_flux_avg

# ifdef BIOLOGY_NPZDOC
      integer avgPAR_flux, avgPARinc_flux
     &     , avgFlux(NumFluxTerms)
     &     , avgVSinkFlux(NumVSinkTerms)
#  ifdef OXYGEN
     &     , avgONNO3, avgONNH4
     &     , avgU10, avgKvO2, avgO2sat
#  ifdef OXYLIM
     &     , avgNCDET
#  endif /* OXYLIM */
#   ifdef CARBON
     &     , avgCNP, avgCNZ, avgrCaCO3orgC
     &     , avgKvCO2, avgCO2sol
#   endif /* CARBON */
#  endif /* OXYGEN */
!#   ifdef OXYGEN
!     &     , avgGasExcFlux(NumGasExcTerms) 
!#   endif /* OXYGEN */
#   ifdef SEDIMENT_BIOLOGY
     &     , avgSedFlux(NumSedFluxTerms)
#   endif /* SEDIMENT_BIOLOGY */
     &     , bgc_flux_avgTime, bgc_flux_avgTstep
     &     , bgc_flux_avgZ

      common /ncids_bgc_flux_avg/ avgPAR_flux, avgPARinc_flux
     &     , avgFlux
     &     , avgVSinkFlux
#  ifdef OXYGEN
     &     , avgONNO3, avgONNH4
     &     , avgU10, avgKvO2, avgO2sat
#  ifdef OXYLIM
     &     , avgNCDET
#  endif
#   ifdef CARBON
     &     , avgCNP, avgCNZ, avgrCaCO3orgC
     &     , avgKvCO2, avgCO2sol
#   endif /* CARBON */
#  endif /* OXYGEN */
!#   ifdef OXYGEN
!     &     , avgGasExcFlux
!#   endif /* OXYGEN */
#   ifdef SEDIMENT_BIOLOGY
     &     , avgSedFlux
#   endif /* SEDIMENT_BIOLOGY */
     &     , bgc_flux_avgTime, bgc_flux_avgTstep
     &     , bgc_flux_avgZ
#endif /* BIOLOGY_NPZDOC */

#  ifdef BIOLOGY_BEC
      integer, dimension(num_bgcflux) :: vid_bec_flux_avg
      integer :: bgc_flux_avgTstep, bgc_flux_avgTime,
     &    bgc_flux_avgZ
      common /c_bgcflux_avg_bec/ vid_bec_flux_avg, bgc_flux_avgTstep,
     &    bgc_flux_avgTime, bgc_flux_avgZ
#  endif
# endif /* AVERAGES */

#endif /* SOLVE3D && BGC_FLUX_ANALYSIS */
