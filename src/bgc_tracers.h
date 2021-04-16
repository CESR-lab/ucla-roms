! DevinD: this was tracers.h in the ETH code.

!
! Biological and other tracers
!

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
# endif
#endif /* BIOLOGY_BEC2 */


