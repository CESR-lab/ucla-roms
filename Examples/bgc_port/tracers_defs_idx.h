! Declare your tracer indices & set to zero:
!
!   You need to declare a tracer index corresponding to the index
!   you used in tracers_defs.h.
!   Set the index equal to zero here as it will be set to final
!   value in tracers_defs.h

! Repo example - Examples/Tracers_passive:
!      integer :: itrace1=0, isalt2=0


! Passive tracers:

! BGC tracers:

      integer, public :: iPO4=0
     & , iNO3=0
     & , iSIO3=0
     & , iNH4=0
     & , iFE=0
     & , iO2=0
     & , iDIC=0
     & , iALK=0
     & , iDOC=0
     & , iDon=0
     & , iDofe=0
     & , iDop=0
     & , iDopr=0
     & , iDonr=0
     & , iZOOC=0
     & , iSPCHL=0
     & , iSPC=0
     & , iSPFE=0
     & , iSPCACO3=0
     & , iDIATCHL=0
     & , iDIATC=0
     & , iDIATFE=0
     & , iDIATSI=0
     & , iDiazchl=0
     & , iDiazc=0
     & , iDiazfe=0

#ifdef Ncycle_SY
     & , iNO2=0
     & , iN2 =0
     & , iN2O=0
# ifdef N2O_TRACER_DECOMP
     & , iN2O_AO1=0
     & , iN2O_SIDEN=0
     & , iN2O_SODEN=0
     & , iN2O_ATM=0
     & , iN2_SED=0
# endif /* N2O_TRACER_DECOMP */
# ifdef N2O_NEV
     & , iN2O_NEV=0
# endif /* N2O_NEV*/
#elif defined N2O_NEV
     & , iN2O_NEV=0
#endif /* Ncycle_SY */
