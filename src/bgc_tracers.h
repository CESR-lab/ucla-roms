! DevinD: this was tracers.h in the ETH code.

!
! Biological and other tracers
!
#if defined BIOLOGY_NPZDOC || defined  LEGACY_NPZD
      integer, parameter :: indxNO3=LAST_I+1
     &           , indxNH4 =indxNO3+1, indxChla=indxNO3+2
     &           , indxPhyt=indxNO3+3, indxZoo =indxNO3+4
     &           , indxSDet=indxNO3+5, indxLDet=indxNO3+6
# undef LAST_I
# define LAST_I indxLDet
# ifdef OXYGEN
     &           , indxO2 = LAST_I+1
#  undef LAST_I
#  define LAST_I indxO2
#  ifdef CARBON
      integer, parameter :: indxDIC=LAST_I+1
     &           , indxTALK=indxDIC+1
     &           , indxSDetC=indxTALK+1
     &           , indxLDetC=indxSDetC+1
     &           , indxCaCO3=indxLDetC+1
#   undef LAST_I
#   define LAST_I indxCaCO3
#  endif /* CARBON */
# endif /* OXYGEN */
#endif /* BIOLOGY_NPZDOC || LEGACY_NPZD*/

#ifdef BIOLOGY_BEC2
      integer, parameter :: indxPO4=indxS+1 ! indxPO4=LAST_I+1 DevinD hard-coded indxS
     &          , indxNo3=indxPO4+1, indxSiO3=indxPO4+2
     &          , indxNh4=indxPO4+3, indxFe=indxPO4+4
     &          , indxO2=indxPO4+5,  indxDic=indxPO4+6
     &          , indxAlk=indxPO4+7, indxDOC=indxPO4+8
     &          , indxDon=indxPO4+9, indxDOFe=indxPO4+10
     &          , indxDop=indxPO4+11, indxDOPr=indxPO4+12
     &          , indxDonr=indxPO4+13, indxZooC=indxPO4+14
     &          , indxSpC=indxPO4+15, indxSpchl=indxPO4+16
     &          , indxSpfe=indxPO4+17, indxSpCaCO3=indxPO4+18
     &          , indxDiatC=indxPO4+19, indxDiatchl=indxPO4+20
     &          , indxDiatfe=indxPO4+21, indxDiatSi=indxPO4+22
     &          , indxDiazC=indxPO4+23, indxDiazchl=indxPO4+24
     &          , indxDiazfe=indxPO4+25
# undef LAST_I
# define LAST_I indxDiazfe
# ifdef BEC_COCCO
     &          , indxCoccoc=indxPO4+26, indxCoccochl=indxPO4+27
     &          , indxCoccocal=indxPO4+28, indxCoccofe=indxPO4+29
     &          , indxCal=indxPO4+30
#  undef LAST_I
#  define LAST_I indxCal
# endif

# ifdef Ncycle_SY
     &          , indxno2=LAST_I+1, indxn2=LAST_I+2
     &          , indxn2o=LAST_I+3
#  undef LAST_I
#  define LAST_I indxn2o
# ifdef N2O_TRACER_DECOMP
     &          , indxn2o_ao1=LAST_I+1, indxn2o_siden=LAST_I+2
     &          , indxn2o_soden=LAST_I+3, indxn2o_atm=LAST_I+4
     &          , indxn2_sed=LAST_I+5
#  undef LAST_I
#  define LAST_I indxn2_sed
# endif
# endif
# ifdef N2O_NEV
     &          , indxn2o_nev=LAST_I+1
#  undef LAST_I
#  define LAST_I indxn2o_nev
# endif
# ifdef USE_EXPLICIT_VSINK
     &          , indxdusthard=LAST_I+1, indxpochard=indxdusthard+1
     &          , indxpcaco3hard=indxdusthard+2, indxpsio2hard=indxdusthard+3
     &          , indxpironhard=indxdusthard+4, indxdustsoft=indxdusthard+5
     &          , indxpocsoft=indxdusthard+6, indxpcaco3soft=indxdusthard+7
     &          , indxpsio2soft=indxdusthard+8, indxpironsoft=indxdusthard+9
#  undef LAST_I
#  define LAST_I indxpironsoft
# endif
#endif /* BIOLOGY_BEC2 */


