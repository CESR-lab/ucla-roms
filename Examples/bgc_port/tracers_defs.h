! This is include file "tracers_defs.h"
!------ --- ---------------------------

! Passive tracers are defined here.

! iTandS represents the last index of the combination
! of 'temp' and 'salt'.
! Number passive tracers sequentially from iTandS+1

! Note: user must also put in index storing integer in
! the top of the tracers.h module for each passive tracer.
! E.g. if you add ptrace1, then at the top of tracers.h
! in the user input section, you need to add: integer iptrace1=0.

! (The reason for seperation is that the index variable needs to be
!  defined at compile time, and is not stored in an array like the
!  the rest of the variables that belong to each tracer.)

! itrace1   = integer to identify tracer in tracer array 't'
! wrt_t     = True/False whether to print tracer to output file
! t_vname   = Tracer short name
! t_units   = Tracer units (for outputing tracer)
! t_lname   = Tracer long name (for outputing tracer)
! t_ana_frc = Whether surf flux is read in (0), or analytical (1)
!             Could be extended 2,3,... depending on different types of
!             analytical forcing, e.g. time variant/invariant, and further
!             specific for each variable using its 'itrace' index...


! - PASSIVE TRACERS:

!      itrace1=1+iTandS;              wrt_t(itrace1) =.True.
!      t_vname(itrace1)='trace1';     t_units(itrace1)='%/%/%'
!      t_tname(itrace1)='trace1_time';t_ana_frc(itrace1)=0
!      t_lname(itrace1)='long trace1'

!      isalt2=2+iTandS;               wrt_t(isalt2) =.True.
!      t_vname(isalt2)='salt2';       t_units(isalt2)='PSUuu'
!      t_tname(isalt2)='salt2_time';  t_ana_frc(isalt2)=0
!      t_lname(isalt2)='long salt2'

! - BGC TRACERS:

!      - Number bgc tracers from 1 and add on itrace index
!        from passive tracers above: e.g. ibgc1 = 1+isalt2
!      - If no passive tracers then add on iTandS index:
!        e.g. ibgc1=1+iTandS
!      - This section of code is modified from ETH code's
!        file init_scalars_bec2.F

!     idea: would be useful to do:
!           indxPO4=itotal+iTandS; itotal=1+itotal
!           This way if you comment out a tracer you no longer want, it doesn't
!           mess up your numbering. Not hard-coded.
!           Should then put in error checking to make sure NT = itotal at the end!

! -- default bgc tracers

      indxPO4=1+iTandS;              wrt_t(indxPO4) =.True.
      t_vname(indxPO4)='PO4';        t_units(indxPO4)='mMol P m-3'
      t_tname(indxPO4)=''; t_ana_frc(indxPO4)=1
      t_lname(indxPO4)='Phosphate'

      indxNO3=2+iTandS;              wrt_t(indxNO3) =.True.
      t_vname(indxNO3)='NO3';        t_units(indxNO3)='mMol N m-3'
      t_tname(indxNO3)='';t_ana_frc(indxNO3)=1
      t_lname(indxNO3)='Nitrate'

      indxSIO3=3+iTandS;             wrt_t(indxSIO3) =.True.
 	  t_vname(indxSIO3)='SiO3';      t_units(indxSIO3)='mMol Si m-3'
 	  t_tname(indxSIO3)='';t_ana_frc(indxSIO3)=1
 	  t_lname(indxSIO3)='Silicate'

 	  indxNH4=4+iTandS;              wrt_t(indxNH4) =.True.
 	  t_vname(indxNH4)='NH4';        t_units(indxNH4)='mMol N m-3'
 	  t_tname(indxNH4)='';t_ana_frc(indxNH4)=1
 	  t_lname(indxNH4)='Ammonium'

      indxFE=5+iTandS;               wrt_t(indxFE) =.True.
	  t_vname(indxFE)='Fe';          t_units(indxFE)='mMol Fe m-3'
	  t_tname(indxFE)='';t_ana_frc(indxFE)=1
	  t_lname(indxFE)='Iron'

      indxO2=6+iTandS;               wrt_t(indxO2) =.True.
	  t_vname(indxO2)='O2';          t_units(indxO2)='mMol O2 m-3'
	  t_tname(indxO2)='';t_ana_frc(indxO2)=1
	  t_lname(indxO2)='Oxygen'

      indxDIC=7+iTandS;              wrt_t(indxDIC) =.True.
	  t_vname(indxDIC)='DIC';        t_units(indxDIC)='mMol C m-3'
	  t_tname(indxDIC)='';t_ana_frc(indxDIC)=1
	  t_lname(indxDIC)='Dissolved inorganic carbon'

      indxALK=8+iTandS;              wrt_t(indxALK) =.True.
	  t_vname(indxALK)='Alk';        t_units(indxALK)='mMol m-3'
	  t_tname(indxALK)='';t_ana_frc(indxALK)=1
	  t_lname(indxALK)='Alkalinity'

      indxDOC=9+iTandS;              wrt_t(indxDOC) =.True.
	  t_vname(indxDOC)='DOC';        t_units(indxDOC)='mMol C m-3'
	  t_tname(indxDOC)='';t_ana_frc(indxDOC)=1
	  t_lname(indxDOC)='Dissolved organic carbon'

      indxDon=10+iTandS;             wrt_t(indxDon) =.True.
	  t_vname(indxDon)='DON';        t_units(indxDon)='mMol N m-3'
	  t_tname(indxDon)='';t_ana_frc(indxDon)=1
	  t_lname(indxDon)='Dissolved organic nitrogen'

      indxDofe=11+iTandS;            wrt_t(indxDofe) =.True.
	  t_vname(indxDofe)='DOFE';      t_units(indxDofe)='mMol Fe m-3'
	  t_tname(indxDofe)='';t_ana_frc(indxDofe)=1
	  t_lname(indxDofe)='Dissolved organic iron'

      indxDop=12+iTandS;             wrt_t(indxDop) =.True.
	  t_vname(indxDop)='DOP';        t_units(indxDop)='mMol P m-3'
	  t_tname(indxDop)='';t_ana_frc(indxDop)=1
	  t_lname(indxDop)='Dissolved organic phosphorus'

      indxDopr=13+iTandS;            wrt_t(indxDopr) =.True.
	  t_vname(indxDopr)='DOPR';      t_units(indxDopr)='mMol P m-3'
	  t_tname(indxDopr)='';t_ana_frc(indxDopr)=1
	  t_lname(indxDopr)='Refractory dissolved organic phosphorus'

      indxDonr=14+iTandS;            wrt_t(indxDonr) =.True.
	  t_vname(indxDonr)='DONR';      t_units(indxDonr)='mMol N m-3'
	  t_tname(indxDonr)='';t_ana_frc(indxDonr)=1
	  t_lname(indxDonr)='Refractory dissolved organic nitrogen'

      indxZOOC=15+iTandS;            wrt_t(indxZOOC) =.True.
	  t_vname(indxZOOC)='ZOOC';      t_units(indxZOOC)='mMol C m-3'
	  t_tname(indxZOOC)='';t_ana_frc(indxZOOC)=1
	  t_lname(indxZOOC)='Zooplankton'

      indxSPCHL=16+iTandS;           wrt_t(indxSPCHL) =.True.
	  t_vname(indxSPCHL)='SPCHL';    t_units(indxSPCHL)='mg Chl-a m-3'
	  t_tname(indxSPCHL)='';t_ana_frc(indxSPCHL)=1
	  t_lname(indxSPCHL)='Small phytoplankton chlorophyll'

      indxSPC=17+iTandS;             wrt_t(indxSPC) =.True.
	  t_vname(indxSPC)='SPC';        t_units(indxSPC)='mMol C m-3'
	  t_tname(indxSPC)='';t_ana_frc(indxSPC)=1
	  t_lname(indxSPC)='Small phytoplankton carbon'

      indxSPFE=18+iTandS;            wrt_t(indxSPFE) =.True.
	  t_vname(indxSPFE)='SPFE';      t_units(indxSPFE)='mMol Fe m-3'
	  t_tname(indxSPFE)='';t_ana_frc(indxSPFE)=1
	  t_lname(indxSPFE)='Small phytoplankton iron'

      indxSPCACO3=19+iTandS;         wrt_t(indxSPCACO3) =.True.
	  t_vname(indxSPCACO3)='SPCACO3';t_units(indxSPCACO3)='mMol CaCO3 m-3'
	  t_tname(indxSPCACO3)='';t_ana_frc(indxSPCACO3)=1
	  t_lname(indxSPCACO3)='Small phytoplankton CaCO3'

      indxDIATCHL=20+iTandS;         wrt_t(indxDIATCHL) =.True.
	  t_vname(indxDIATCHL)='DIATCHL';t_units(indxDIATCHL)='mg Chl-a m-3'
	  t_tname(indxDIATCHL)='';t_ana_frc(indxDIATCHL)=1
	  t_lname(indxDIATCHL)='Diatom chlorophyll'

      indxDIATC=21+iTandS;           wrt_t(indxDIATC) =.True.
	  t_vname(indxDIATC)='DIATC';    t_units(indxDIATC)='mMol C m-3'
	  t_tname(indxDIATC)='';t_ana_frc(indxDIATC)=1
	  t_lname(indxDIATC)='Diatom carbon'

      indxDIATFE=22+iTandS;          wrt_t(indxDIATFE) =.True.
	  t_vname(indxDIATFE)='DIATFE';  t_units(indxDIATFE)='mMol Fe m-3'
	  t_tname(indxDIATFE)='';t_ana_frc(indxDIATFE)=1
	  t_lname(indxDIATFE)='Diatom Iron'

      indxDIATSI=23+iTandS;          wrt_t(indxDIATSI) =.True.
	  t_vname(indxDIATSI)='DIATSI';  t_units(indxDIATSI)='mMol Si m-3'
	  t_tname(indxDIATSI)='';t_ana_frc(indxDIATSI)=1
	  t_lname(indxDIATSI)='Diatom silicon'

      indxDiazchl=24+iTandS;         wrt_t(indxDiazchl) =.True.
	  t_vname(indxDiazchl)='DIAZCHL';t_units(indxDiazchl)='mg Chl-a m-3'
	  t_tname(indxDiazchl)='';t_ana_frc(indxDiazchl)=1
	  t_lname(indxDiazchl)='Diazotroph chlorophyll'

      indxDiazc=25+iTandS;           wrt_t(indxDiazc) =.True.
	  t_vname(indxDiazc)='DIAZC';    t_units(indxDiazc)='mMol C m-3'
	  t_tname(indxDiazc)='';t_ana_frc(indxDiazc)=1
	  t_lname(indxDiazc)='Diazotroph carbon'

      indxDiazfe=26+iTandS;          wrt_t(indxDiazfe) =.True.
	  t_vname(indxDiazfe)='DIAZFE';  t_units(indxDiazfe)='mMol Fe m-3'
	  t_tname(indxDiazfe)='';t_ana_frc(indxDiazfe)=1
	  t_lname(indxDiazfe)='Diazotroph iron'

! # ifdef BEC_COCCO

!	  =XX+iTandS;              wrt_t() =.True.
!	  t_vname()='';       t_units()=''
!	  t_tname()='';t_ana_frc()=1
!	  t_lname()=''











