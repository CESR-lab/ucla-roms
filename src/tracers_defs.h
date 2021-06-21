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
! wrt_t     = True/False whether to write tracer to output file
! wrt_t_avg = True/False whether to write tracer avg to output file
! t_vname   = Tracer short name
! t_units   = Tracer units (for outputing tracer)
! t_lname   = Tracer long name (for outputing tracer)
! t_ana_frc = Whether surf flux is read in (0), or analytical (1)
!             Could be extended 2,3,... depending on different types of
!             analytical forcing, e.g. time variant/invariant, and further
!             specific for each variable using its 'itrace' index...

! - TEMP & SALT:

!     Defined in tracers.F, but outputting controller here:
      wrt_t(itmp) =.True.;           wrt_t_avg(itmp) =.True.
      wrt_t_diag(itmp) =.False.
#ifdef SALINITY
      wrt_t(islt) =.True.;           wrt_t_avg(islt) =.True.
      wrt_t_diag(islt) =.True.
#endif


! - PASSIVE TRACERS:

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

!    For code-dev only, don't bgc for physics only.
#ifdef BIOLOGY_BEC2

      ! itrc_bio=1+iTandS ! Starting tracer index for bgc tracers. Set in tracers_defs.h

      iPO4=1+iTandS;
      wrt_t(iPO4) =.True.;        wrt_t_avg(iPO4) =.True.
      t_vname(iPO4)='PO4';        t_units(iPO4)='mMol P m-3'
      t_tname(iPO4)=''; t_ana_frc(iPO4)=1
      t_lname(iPO4)='Phosphate'

      iNO3=2+iTandS;
      wrt_t(iNO3) =.True.;        wrt_t_avg(iNO3) =.True.
      t_vname(iNO3)='NO3';        t_units(iNO3)='mMol N m-3'
      t_tname(iNO3)='';t_ana_frc(iNO3)=1
      t_lname(iNO3)='Nitrate'

      iSIO3=3+iTandS;
      wrt_t(iSIO3) =.True.;       wrt_t_avg(iSIO3) =.True.
 	  t_vname(iSIO3)='SiO3';      t_units(iSIO3)='mMol Si m-3'
 	  t_tname(iSIO3)='';t_ana_frc(iSIO3)=1
 	  t_lname(iSIO3)='Silicate'

 	  iNH4=4+iTandS;
 	  wrt_t(iNH4) =.True.;        wrt_t_avg(iNH4) =.True.
 	  t_vname(iNH4)='NH4';        t_units(iNH4)='mMol N m-3'
 	  t_tname(iNH4)='';t_ana_frc(iNH4)=1
 	  t_lname(iNH4)='Ammonium'

      iFE=5+iTandS;
 	  wrt_t(iFE) =.True.;         wrt_t_avg(iFE) =.True.
	  t_vname(iFE)='Fe';          t_units(iFE)='mMol Fe m-3'
	  t_tname(iFE)='';t_ana_frc(iFE)=1
	  t_lname(iFE)='Iron'

      iO2=6+iTandS;
	  wrt_t(iO2) =.True.;         wrt_t_avg(iO2) =.True.
	  t_vname(iO2)='O2';          t_units(iO2)='mMol O2 m-3'
	  t_tname(iO2)='';t_ana_frc(iO2)=1
	  t_lname(iO2)='Oxygen'

      iDIC=7+iTandS;
	  wrt_t(iDIC) =.True.;        wrt_t_avg(iDIC) =.True.
	  t_vname(iDIC)='DIC';        t_units(iDIC)='mMol C m-3'
	  t_tname(iDIC)='';t_ana_frc(iDIC)=1
	  t_lname(iDIC)='Dissolved inorganic carbon'

      iALK=8+iTandS;
	  wrt_t(iALK) =.True.;        wrt_t_avg(iALK) =.True.
	  t_vname(iALK)='Alk';        t_units(iALK)='mMol m-3'
	  t_tname(iALK)='';t_ana_frc(iALK)=1
	  t_lname(iALK)='Alkalinity'

      iDOC=9+iTandS;
	  wrt_t(iDOC) =.True.;        wrt_t_avg(iDOC) =.True.
	  t_vname(iDOC)='DOC';        t_units(iDOC)='mMol C m-3'
	  t_tname(iDOC)='';t_ana_frc(iDOC)=1
	  t_lname(iDOC)='Dissolved organic carbon'

      iDon=10+iTandS;
	  wrt_t(iDon) =.True.;        wrt_t_avg(iDon) =.True.
	  t_vname(iDon)='DON';        t_units(iDon)='mMol N m-3'
	  t_tname(iDon)='';t_ana_frc(iDon)=1
	  t_lname(iDon)='Dissolved organic nitrogen'

      iDofe=11+iTandS;
	  wrt_t(iDofe) =.True.;       wrt_t_avg(iDofe) =.True.
	  t_vname(iDofe)='DOFE';      t_units(iDofe)='mMol Fe m-3'
	  t_tname(iDofe)='';t_ana_frc(iDofe)=1
	  t_lname(iDofe)='Dissolved organic iron'

      iDop=12+iTandS;
	  wrt_t(iDop) =.True.;        wrt_t_avg(iDop) =.True.
	  t_vname(iDop)='DOP';        t_units(iDop)='mMol P m-3'
	  t_tname(iDop)='';t_ana_frc(iDop)=1
	  t_lname(iDop)='Dissolved organic phosphorus'

      iDopr=13+iTandS;
	  wrt_t(iDopr) =.True.;       wrt_t_avg(iDopr) =.True.
	  t_vname(iDopr)='DOPR';      t_units(iDopr)='mMol P m-3'
	  t_tname(iDopr)='';t_ana_frc(iDopr)=1
	  t_lname(iDopr)='Refractory dissolved organic phosphorus'

      iDonr=14+iTandS;
	  wrt_t(iDonr) =.True.;       wrt_t_avg(iDonr) =.True.
	  t_vname(iDonr)='DONR';      t_units(iDonr)='mMol N m-3'
	  t_tname(iDonr)='';t_ana_frc(iDonr)=1
	  t_lname(iDonr)='Refractory dissolved organic nitrogen'

      iZOOC=15+iTandS;
	  wrt_t(iZOOC) =.True.;       wrt_t_avg(iZOOC) =.True.
	  t_vname(iZOOC)='ZOOC';      t_units(iZOOC)='mMol C m-3'
	  t_tname(iZOOC)='';t_ana_frc(iZOOC)=1
	  t_lname(iZOOC)='Zooplankton'

	  iSPC=16+iTandS;
	  wrt_t(iSPC) =.True.;        wrt_t_avg(iSPC) =.True.
	  t_vname(iSPC)='SPC';        t_units(iSPC)='mMol C m-3'
	  t_tname(iSPC)='';t_ana_frc(iSPC)=1
	  t_lname(iSPC)='Small phytoplankton carbon'

      iSPCHL=17+iTandS;
	  wrt_t(iSPCHL) =.True.;      wrt_t_avg(iSPCHL) =.True.
	  t_vname(iSPCHL)='SPCHL';    t_units(iSPCHL)='mg Chl-a m-3'
	  t_tname(iSPCHL)='';t_ana_frc(iSPCHL)=1
	  t_lname(iSPCHL)='Small phytoplankton chlorophyll'

      iSPFE=18+iTandS;
	  wrt_t(iSPFE) =.True.;       wrt_t_avg(iSPFE) =.True.
	  t_vname(iSPFE)='SPFE';      t_units(iSPFE)='mMol Fe m-3'
	  t_tname(iSPFE)='';t_ana_frc(iSPFE)=1
	  t_lname(iSPFE)='Small phytoplankton iron'

      iSPCACO3=19+iTandS;
	  wrt_t(iSPCACO3) =.True.;    wrt_t_avg(iSPCACO3) =.True.
	  t_vname(iSPCACO3)='SPCACO3';t_units(iSPCACO3)='mMol CaCO3 m-3'
	  t_tname(iSPCACO3)='';t_ana_frc(iSPCACO3)=1
	  t_lname(iSPCACO3)='Small phytoplankton CaCO3'

	  iDIATC=20+iTandS;
	  wrt_t(iDIATC) =.True.;      wrt_t_avg(iDIATC) =.True.
	  t_vname(iDIATC)='DIATC';    t_units(iDIATC)='mMol C m-3'
	  t_tname(iDIATC)='';t_ana_frc(iDIATC)=1
	  t_lname(iDIATC)='Diatom carbon'

      iDIATCHL=21+iTandS;
	  wrt_t(iDIATCHL) =.True.;    wrt_t_avg(iDIATCHL) =.True.
	  t_vname(iDIATCHL)='DIATCHL';t_units(iDIATCHL)='mg Chl-a m-3'
	  t_tname(iDIATCHL)='';t_ana_frc(iDIATCHL)=1
	  t_lname(iDIATCHL)='Diatom chlorophyll'

      iDIATFE=22+iTandS;
	  wrt_t(iDIATFE) =.True.;     wrt_t_avg(iDIATFE) =.True.
	  t_vname(iDIATFE)='DIATFE';  t_units(iDIATFE)='mMol Fe m-3'
	  t_tname(iDIATFE)='';t_ana_frc(iDIATFE)=1
	  t_lname(iDIATFE)='Diatom Iron'

      iDIATSI=23+iTandS;
	  wrt_t(iDIATSI) =.True.;     wrt_t_avg(iDIATSI) =.True.
	  t_vname(iDIATSI)='DIATSI';  t_units(iDIATSI)='mMol Si m-3'
	  t_tname(iDIATSI)='';t_ana_frc(iDIATSI)=1
	  t_lname(iDIATSI)='Diatom silicon'

	  iDiazc=24+iTandS;
	  wrt_t(iDiazc) =.True.;      wrt_t_avg(iDiazc) =.True.
	  t_vname(iDiazc)='DIAZC';    t_units(iDiazc)='mMol C m-3'
	  t_tname(iDiazc)='';t_ana_frc(iDiazc)=1
	  t_lname(iDiazc)='Diazotroph carbon'

      iDiazchl=25+iTandS;
	  wrt_t(iDiazchl) =.True.;    wrt_t_avg(iDiazchl) =.True.
	  t_vname(iDiazchl)='DIAZCHL';t_units(iDiazchl)='mg Chl-a m-3'
	  t_tname(iDiazchl)='';t_ana_frc(iDiazchl)=1
	  t_lname(iDiazchl)='Diazotroph chlorophyll'

      iDiazfe=26+iTandS;
	  wrt_t(iDiazfe) =.True.;     wrt_t_avg(iDiazfe) =.True.
	  t_vname(iDiazfe)='DIAZFE';  t_units(iDiazfe)='mMol Fe m-3'
	  t_tname(iDiazfe)='';t_ana_frc(iDiazfe)=1
	  t_lname(iDiazfe)='Diazotroph iron'

      !ntrc_bio_base=26 ! Total number of base bgc tracers. Hard-coded for now. Use itot later.

#define LAST_I iDIAZFE

#ifdef Ncycle_SY
      iNO2=LAST_I+1;
	  wrt_t(iNO2) =.True.;        wrt_t_avg(iNO2) =.True.
	  t_vname(iNO2)='NO2';        t_units(iNO2)='mMol N m-3'
	  t_tname(iNO2)='';t_ana_frc(iNO2)=1
	  t_lname(iNO2)='Nitrite'

      iN2=LAST_I+2;
	  wrt_t(iN2) =.True.;         wrt_t_avg(iN2) =.True.
	  t_vname(iN2)='N2';          t_units(iN2)='mMol N2 m-3'
	  t_tname(iN2)='';t_ana_frc(iN2)=1
	  t_lname(iN2)='Dinitrogen'

      iN2O=LAST_I+3;
	  wrt_t(iN2O) =.True.;        wrt_t_avg(iN2O) =.True.
	  t_vname(iN2O)='N2O';        t_units(iN2O)='mMol N2O m-3'
	  t_tname(iN2O)='';t_ana_frc(iN2O)=1
	  t_lname(iN2O)='Nitrous oxide'
# undef LAST_I
# define LAST_I iN2O
#endif /* Ncycle_SY */

!	  i=               wrt_t() =.True.
!	  t_vname()='';        t_units()=''
!	  t_tname()='';t_ana_frc()=1
!	  t_lname()=''

      ! total number of bgc tracers
      !ntrc_bio=ntrc_bio_base ! +ntrc_bio_cocco +ntrc_bio_ncycle


#endif /* BIOLOGY_BEC2 */


