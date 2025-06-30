!
! - BGC DIAGNOSTICS
!
! nr_bec2_diag_3d : number of available 3d diags set in bgc_ecosys_vars
!
! MARBL USERS NOTE:
!--------------------------------------------------------------------------------
!      if running with MARBL, do not specify BGC diagnostics in this file: they are
!      provided to ROMS by MARBL. Control the MARBL config (therefore diagnostics)
!      using a "marbl_in" file in your run directory.
!      Regarding output: by default, all MARBL diagnostics are written to output.
!      To control which particular diagnostics to write out, make a
!      text file "marbl_diagnostic_output_list" in your run directory
!      and add the shortname of each required diagnostic to a new line.
!      If the text file is empty, does not exist,
!      contains only comments (!) or no recognised diagnostics,
!      ROMS will revert to default behaviour and output all MARBL diagnostics.
!      Output frequency etc. are still controlled using "bgc.opt".
!      For a list of available MARBL diagnostics, see
!      $MARBL_ROOT/tests/regression_tests/requested_diags/requested_diags.py


#ifdef BIOLOGY_BEC2
#ifdef BEC2_DIAG

!!!!!!!!!!!! bec2_diag_3d !!!!!!!!!!!!!!!!!!!

      par_idx_t=1
      vname_bec2_diag_3d(1,par_idx_t)='PAR'
      vname_bec2_diag_3d(2,par_idx_t)='Photosynthetically active radiation'
      vname_bec2_diag_3d(3,par_idx_t)='W m-2'
      vname_bec2_diag_3d(4,par_idx_t)='  '
      wrt_bec2_diag_3d(par_idx_t)=.True.

      pocfluxin_idx_t=2
      vname_bec2_diag_3d(1,pocfluxin_idx_t)='POC_FLUX_IN'
      vname_bec2_diag_3d(2,pocfluxin_idx_t)='POC flux into cell'
      vname_bec2_diag_3d(3,pocfluxin_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,pocfluxin_idx_t)=' '
      wrt_bec2_diag_3d(pocfluxin_idx_t)=.True.

      pocprod_idx_t=3
      vname_bec2_diag_3d(1,pocprod_idx_t)='POC_PROD'
      vname_bec2_diag_3d(2,pocprod_idx_t)='POC production'
      vname_bec2_diag_3d(3,pocprod_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,pocprod_idx_t)=' '
      wrt_bec2_diag_3d(pocprod_idx_t)=.false.

      pocremin_idx_t=4
      vname_bec2_diag_3d(1,pocremin_idx_t)='POC_REMIN'
      vname_bec2_diag_3d(2,pocremin_idx_t)='POC remineralization'
      vname_bec2_diag_3d(3,pocremin_idx_t)='mmol C/m2/s'
      vname_bec2_diag_3d(4,pocremin_idx_t)=' '
      wrt_bec2_diag_3d(pocremin_idx_t)=.false.

      caco3fluxin_idx_t=5
      vname_bec2_diag_3d(1,caco3fluxin_idx_t)='CACO3_FLUX_IN'
      vname_bec2_diag_3d(2,caco3fluxin_idx_t)='Incoming flux of large detrital CaCO3'
      vname_bec2_diag_3d(3,caco3fluxin_idx_t)='mmol CaCO3/m2/s '
      vname_bec2_diag_3d(4,caco3fluxin_idx_t)=' '
      wrt_bec2_diag_3d(caco3fluxin_idx_t)=.false.

      pcaco3prod_idx_t=6
      vname_bec2_diag_3d(1,pcaco3prod_idx_t)='PCACO3_PROD'
      vname_bec2_diag_3d(2,pcaco3prod_idx_t)='Production of large detrital CaCO3'
      vname_bec2_diag_3d(3,pcaco3prod_idx_t)='mmol CaCO3/m2/s '
      vname_bec2_diag_3d(4,pcaco3prod_idx_t)=' '
      wrt_bec2_diag_3d(pcaco3prod_idx_t)=.false.

      caco3remin_idx_t=7
      vname_bec2_diag_3d(1,caco3remin_idx_t)='CACO3_REMIN'
      vname_bec2_diag_3d(2,caco3remin_idx_t)='Remineralization of large detrital CaCO3'
      vname_bec2_diag_3d(3,caco3remin_idx_t)='mmol CaCO3/m2/s '
      vname_bec2_diag_3d(4,caco3remin_idx_t)=' '
      wrt_bec2_diag_3d(caco3remin_idx_t)=.false.

      sio2fluxin_idx_t=8
      vname_bec2_diag_3d(1,sio2fluxin_idx_t)='SIO2_FLUX_IN'
      vname_bec2_diag_3d(2,sio2fluxin_idx_t)='Incoming flux of large detritus SiO2'
      vname_bec2_diag_3d(3,sio2fluxin_idx_t)='mmol Si/m2/s '
      vname_bec2_diag_3d(4,sio2fluxin_idx_t)=' '
      wrt_bec2_diag_3d(sio2fluxin_idx_t)=.false.

      sio2prod_idx_t=9
      vname_bec2_diag_3d(1,sio2prod_idx_t)='SIO2_PROD'
      vname_bec2_diag_3d(2,sio2prod_idx_t)='Production of large detritus SiO2'
      vname_bec2_diag_3d(3,sio2prod_idx_t)='mmol Si/m2/s '
      vname_bec2_diag_3d(4,sio2prod_idx_t)=' '
      wrt_bec2_diag_3d(sio2prod_idx_t)=.false.

      sio2remin_idx_t=10
      vname_bec2_diag_3d(1,sio2remin_idx_t)='SIO2_REMIN'
      vname_bec2_diag_3d(2,sio2remin_idx_t)='Remineralization of large detritus SiO2'
      vname_bec2_diag_3d(3,sio2remin_idx_t)='mmol Si/m2/s'
      vname_bec2_diag_3d(4,sio2remin_idx_t)=' '
      wrt_bec2_diag_3d(sio2remin_idx_t)=.false.

      dustfluxin_idx_t=11
      vname_bec2_diag_3d(1,dustfluxin_idx_t)='DUST_FLUX_IN'
      vname_bec2_diag_3d(2,dustfluxin_idx_t)='Incoming dust flux'
      vname_bec2_diag_3d(3,dustfluxin_idx_t)='kg dust/m2/s '
      vname_bec2_diag_3d(4,dustfluxin_idx_t)=' '
      wrt_bec2_diag_3d(dustfluxin_idx_t)=.false.

      dustremin_idx_t=12
      vname_bec2_diag_3d(1,dustremin_idx_t)='DUST_REMIN'
      vname_bec2_diag_3d(2,dustremin_idx_t)='Remineralization of dust'
      vname_bec2_diag_3d(3,dustremin_idx_t)='kg dust/m2/s '
      vname_bec2_diag_3d(4,dustremin_idx_t)=' '
      wrt_bec2_diag_3d(dustremin_idx_t)=.false.

      pironfluxin_idx_t=13
      vname_bec2_diag_3d(1,pironfluxin_idx_t)='P_IRON_FLUX_IN'
      vname_bec2_diag_3d(2,pironfluxin_idx_t)='Incoming flux of large detritus iron'
      vname_bec2_diag_3d(3,pironfluxin_idx_t)='mmol Fe/m2/s '
      vname_bec2_diag_3d(4,pironfluxin_idx_t)=' '
      wrt_bec2_diag_3d(pironfluxin_idx_t)=.false.

      pironprod_idx_t=14
      vname_bec2_diag_3d(1,pironprod_idx_t)='P_IRON_PROD'
      vname_bec2_diag_3d(2,pironprod_idx_t)='Production of large detritus iron'
      vname_bec2_diag_3d(3,pironprod_idx_t)='mmol Fe/m2/s '
      vname_bec2_diag_3d(4,pironprod_idx_t)=' '
      wrt_bec2_diag_3d(pironprod_idx_t)=.false.

      pironremin_idx_t=15
      vname_bec2_diag_3d(1,pironremin_idx_t)='P_IRON_REMIN'
      vname_bec2_diag_3d(2,pironremin_idx_t)='Remineralization of large detritus iron'
      vname_bec2_diag_3d(3,pironremin_idx_t)='mmol Fe/m2/s '
      vname_bec2_diag_3d(4,pironremin_idx_t)=' '
      wrt_bec2_diag_3d(pironremin_idx_t)=.false.

      grazesp_idx_t=16
      vname_bec2_diag_3d(1,grazesp_idx_t)='GRAZE_SP'
      vname_bec2_diag_3d(2,grazesp_idx_t)='Grazing rate on small phytoplankton'
      vname_bec2_diag_3d(3,grazesp_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,grazesp_idx_t)=' '
      wrt_bec2_diag_3d(grazesp_idx_t)=.false.

      grazediat_idx_t=17
      vname_bec2_diag_3d(1,grazediat_idx_t)='GRAZE_DIAT'
      vname_bec2_diag_3d(2,grazediat_idx_t)='Grazing rate on diatoms'
      vname_bec2_diag_3d(3,grazediat_idx_t)='mmol C/m2/s'
      vname_bec2_diag_3d(4,grazediat_idx_t)=' '
      wrt_bec2_diag_3d(grazediat_idx_t)=.false.

      grazediaz_idx_t=18
      vname_bec2_diag_3d(1,grazediaz_idx_t)='GRAZE_DIAZ'
      vname_bec2_diag_3d(2,grazediaz_idx_t)='Grazing rate on diazotrophs'
      vname_bec2_diag_3d(3,grazediaz_idx_t)='mmol C/m2/s'
      vname_bec2_diag_3d(4,grazediaz_idx_t)=' '
      wrt_bec2_diag_3d(grazediaz_idx_t)=.false.

      sploss_idx_t=19
      vname_bec2_diag_3d(1,sploss_idx_t)='SP_LOSS'
      vname_bec2_diag_3d(2,sploss_idx_t)='Small phytoplankton non-grazing mortality'
      vname_bec2_diag_3d(3,sploss_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,sploss_idx_t)=' '
      wrt_bec2_diag_3d(sploss_idx_t)=.false.

      diatloss_idx_t=20
      vname_bec2_diag_3d(1,diatloss_idx_t)='DIAT_LOSS'
      vname_bec2_diag_3d(2,diatloss_idx_t)='Diatom non-grazing mortality'
      vname_bec2_diag_3d(3,diatloss_idx_t)='mmol/m2/s '
      vname_bec2_diag_3d(4,diatloss_idx_t)=' '
      wrt_bec2_diag_3d(diatloss_idx_t)=.false.

      zooloss_idx_t=21
      vname_bec2_diag_3d(1,zooloss_idx_t)='ZOO_LOSS'
      vname_bec2_diag_3d(2,zooloss_idx_t)='Mortality due to higher trophic grazing on zooplankton'
      vname_bec2_diag_3d(3,zooloss_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,zooloss_idx_t)=' '
      wrt_bec2_diag_3d(zooloss_idx_t)=.false.

      spagg_idx_t=22
      vname_bec2_diag_3d(1,spagg_idx_t)='SP_AGG'
      vname_bec2_diag_3d(2,spagg_idx_t)='Aggregation of small phytoplankton'
      vname_bec2_diag_3d(3,spagg_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,spagg_idx_t)=' '
      wrt_bec2_diag_3d(spagg_idx_t)=.false.

      diatagg_idx_t=23
      vname_bec2_diag_3d(1,diatagg_idx_t)='DIAT_AGG'
      vname_bec2_diag_3d(2,diatagg_idx_t)='Aggregation of diatoms'
      vname_bec2_diag_3d(3,diatagg_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,diatagg_idx_t)=' '
      wrt_bec2_diag_3d(diatagg_idx_t)=.false.

      photocsp_idx_t=24
      vname_bec2_diag_3d(1,photocsp_idx_t)='PHOTOC_SP'
      vname_bec2_diag_3d(2,photocsp_idx_t)='C fixation rate by small phytoplankton'
      vname_bec2_diag_3d(3,photocsp_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,photocsp_idx_t)=' '
      wrt_bec2_diag_3d(photocsp_idx_t)=.false.

      photocdiat_idx_t=25
      vname_bec2_diag_3d(1,photocdiat_idx_t)='PHOTOC_DIAT'
      vname_bec2_diag_3d(2,photocdiat_idx_t)='C fixation rate by diatoms'
      vname_bec2_diag_3d(3,photocdiat_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,photocdiat_idx_t)=' '
      wrt_bec2_diag_3d(photocdiat_idx_t)=.false.

      totprod_idx_t=26
      vname_bec2_diag_3d(1,totprod_idx_t)='TOT_PROD'
      vname_bec2_diag_3d(2,totprod_idx_t)='Total autotroph production'
      vname_bec2_diag_3d(3,totprod_idx_t)='mmol/m3/s '
      vname_bec2_diag_3d(4,totprod_idx_t)=' '
      wrt_bec2_diag_3d(totprod_idx_t)=.true.

      docprod_idx_t=27
      vname_bec2_diag_3d(1,docprod_idx_t)='DOC_PROD'
      vname_bec2_diag_3d(2,docprod_idx_t)='DOC production'
      vname_bec2_diag_3d(3,docprod_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,docprod_idx_t)=' '
      wrt_bec2_diag_3d(docprod_idx_t)=.false.

      docremin_idx_t=28
      vname_bec2_diag_3d(1,docremin_idx_t)='DOC_REMIN'
      vname_bec2_diag_3d(2,docremin_idx_t)='Remineralization rate of DOC'
      vname_bec2_diag_3d(3,docremin_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,docremin_idx_t)=' '
      wrt_bec2_diag_3d(docremin_idx_t)=.false.

      fescavenge_idx_t=29
      vname_bec2_diag_3d(1,fescavenge_idx_t)='FE_SCAVENGE'
      vname_bec2_diag_3d(2,fescavenge_idx_t)='Loss of dissolved iron, scavenging'
      vname_bec2_diag_3d(3,fescavenge_idx_t)='mmol/m2/s '
      vname_bec2_diag_3d(4,fescavenge_idx_t)=' '
      wrt_bec2_diag_3d(fescavenge_idx_t)=.false.

      spnlim_idx_t=30
      vname_bec2_diag_3d(1,spnlim_idx_t)='SP_N_LIM'
      vname_bec2_diag_3d(2,spnlim_idx_t)='N limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,spnlim_idx_t)=' '
      vname_bec2_diag_3d(4,spnlim_idx_t)=' '
      wrt_bec2_diag_3d(spnlim_idx_t)=.false.

      spfeuptake_idx_t=31
      vname_bec2_diag_3d(1,spfeuptake_idx_t)='SP_FE_LIM'
      vname_bec2_diag_3d(2,spfeuptake_idx_t)='Fe limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,spfeuptake_idx_t)=' '
      vname_bec2_diag_3d(4,spfeuptake_idx_t)=' '
      wrt_bec2_diag_3d(spfeuptake_idx_t)=.false.

      sppo4uptake_idx_t=32
      vname_bec2_diag_3d(1,sppo4uptake_idx_t)='SP_PO4_LIM'
      vname_bec2_diag_3d(2,sppo4uptake_idx_t)='PO4 limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,sppo4uptake_idx_t)=' '
      vname_bec2_diag_3d(4,sppo4uptake_idx_t)=' '
      wrt_bec2_diag_3d(sppo4uptake_idx_t)=.false.

      splightlim_idx_t=33
      vname_bec2_diag_3d(1,splightlim_idx_t)='SP_LIGHT_LIM'
      vname_bec2_diag_3d(2,splightlim_idx_t)='Small phytoplankton light limitation factor'
      vname_bec2_diag_3d(3,splightlim_idx_t)=' '
      vname_bec2_diag_3d(4,splightlim_idx_t)=' '
      wrt_bec2_diag_3d(splightlim_idx_t)=.false.

      diatnlim_idx_t=34
      vname_bec2_diag_3d(1,diatnlim_idx_t)='DIAT_N_LIM'
      vname_bec2_diag_3d(2,diatnlim_idx_t)='N limitation factor of diatoms'
      vname_bec2_diag_3d(3,diatnlim_idx_t)=' '
      vname_bec2_diag_3d(4,diatnlim_idx_t)=' '
      wrt_bec2_diag_3d(diatnlim_idx_t)=.false.

      diatfeuptake_idx_t=35
      vname_bec2_diag_3d(1,diatfeuptake_idx_t)='DIAT_FE_LIM'
      vname_bec2_diag_3d(2,diatfeuptake_idx_t)='Fe limitation factor of diatoms'
      vname_bec2_diag_3d(3,diatfeuptake_idx_t)=' '
      vname_bec2_diag_3d(4,diatfeuptake_idx_t)=' '
      wrt_bec2_diag_3d(diatfeuptake_idx_t)=.false.

      diatpo4uptake_idx_t=36
      vname_bec2_diag_3d(1,diatpo4uptake_idx_t)='DIAT_PO4_LIM'
      vname_bec2_diag_3d(2,diatpo4uptake_idx_t)='PO4 limitation factor of diatoms'
      vname_bec2_diag_3d(3,diatpo4uptake_idx_t)=' '
      vname_bec2_diag_3d(4,diatpo4uptake_idx_t)=' '
      wrt_bec2_diag_3d(diatpo4uptake_idx_t)=.false.

      diatsio3uptake_idx_t=37
      vname_bec2_diag_3d(1,diatsio3uptake_idx_t)='DIAT_SIO3_LIM'
      vname_bec2_diag_3d(2,diatsio3uptake_idx_t)='SiO3 limitation factor of diatoms'
      vname_bec2_diag_3d(3,diatsio3uptake_idx_t)=' '
      vname_bec2_diag_3d(4,diatsio3uptake_idx_t)=' '
      wrt_bec2_diag_3d(diatsio3uptake_idx_t)=.false.

      diatlightlim_idx_t=38
      vname_bec2_diag_3d(1,diatlightlim_idx_t)='DIAT_LIGHT_LIM'
      vname_bec2_diag_3d(2,diatlightlim_idx_t)='Diatom light limitation factor'
      vname_bec2_diag_3d(3,diatlightlim_idx_t)=' '
      vname_bec2_diag_3d(4,diatlightlim_idx_t)=' '
      wrt_bec2_diag_3d(diatlightlim_idx_t)=.false.

      caco3prod_idx_t=39
      vname_bec2_diag_3d(1,caco3prod_idx_t)='CACO3_PROD'
      vname_bec2_diag_3d(2,caco3prod_idx_t)='Production of CaCO3 by small phytoplankton'
      vname_bec2_diag_3d(3,caco3prod_idx_t)='mmol CaCO3/m2/s '
      vname_bec2_diag_3d(4,caco3prod_idx_t)=' '
      wrt_bec2_diag_3d(caco3prod_idx_t)=.false.

      diaznfix_idx_t=40
      vname_bec2_diag_3d(1,diaznfix_idx_t)='DIAZ_NFIX'
      vname_bec2_diag_3d(2,diaznfix_idx_t)='Total nitrogen fixation by diazotrophs'
      vname_bec2_diag_3d(3,diaznfix_idx_t)='mmol/m2/s   '
      vname_bec2_diag_3d(4,diaznfix_idx_t)=' '
      wrt_bec2_diag_3d(diaznfix_idx_t)=.false.

      diazloss_idx_t=41
      vname_bec2_diag_3d(1,diazloss_idx_t)='DIAZ_LOSS'
      vname_bec2_diag_3d(2,diazloss_idx_t)='Diazotroph non-grazing mortality'
      vname_bec2_diag_3d(3,diazloss_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,diazloss_idx_t)=' '
      wrt_bec2_diag_3d(diazloss_idx_t)=.false.

      photocdiaz_idx_t=42
      vname_bec2_diag_3d(1,photocdiaz_idx_t)='PHOTOC_DIAZ'
      vname_bec2_diag_3d(2,photocdiaz_idx_t)='Diazotroph C-fixation'
      vname_bec2_diag_3d(3,photocdiaz_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,photocdiaz_idx_t)=' '
      wrt_bec2_diag_3d(photocdiaz_idx_t)=.false.

      diazpo4uptake_idx_t=43
      vname_bec2_diag_3d(1,diazpo4uptake_idx_t)='DIAZ_PO4_LIM'
      vname_bec2_diag_3d(2,diazpo4uptake_idx_t)='PO4 limitation factor of diazotrophs'
      vname_bec2_diag_3d(3,diazpo4uptake_idx_t)=' '
      vname_bec2_diag_3d(4,diazpo4uptake_idx_t)=' '
      wrt_bec2_diag_3d(diazpo4uptake_idx_t)=.false.

      diazfeuptake_idx_t=44
      vname_bec2_diag_3d(1,diazfeuptake_idx_t)='DIAZ_FE_LIM'
      vname_bec2_diag_3d(2,diazfeuptake_idx_t)='Fe limitation factor of diazotrophs'
      vname_bec2_diag_3d(3,diazfeuptake_idx_t)=' '
      vname_bec2_diag_3d(4,diazfeuptake_idx_t)=' '
      wrt_bec2_diag_3d(diazfeuptake_idx_t)=.false.

      diazlightlim_idx_t=45
      vname_bec2_diag_3d(1,diazlightlim_idx_t)='DIAZ_LIGHT_LIM'
      vname_bec2_diag_3d(2,diazlightlim_idx_t)='Diazotroph light limitation'
      vname_bec2_diag_3d(3,diazlightlim_idx_t)=' '
      vname_bec2_diag_3d(4,diazlightlim_idx_t)=' '
      wrt_bec2_diag_3d(diazlightlim_idx_t)=.false.

      fescavengerate_idx_t=46
      vname_bec2_diag_3d(1,fescavengerate_idx_t)='FE_SCAVENGE_RATE'
      vname_bec2_diag_3d(2,fescavengerate_idx_t)='Annual scavenging rate of iron as % of ambient'
      vname_bec2_diag_3d(3,fescavengerate_idx_t)=' '
      vname_bec2_diag_3d(4,fescavengerate_idx_t)=' '
      wrt_bec2_diag_3d(fescavengerate_idx_t)=.false.

      donprod_idx_t=47
      vname_bec2_diag_3d(1,donprod_idx_t)='DON_PROD'
      vname_bec2_diag_3d(2,donprod_idx_t)='Production of dissolved organic N'
      vname_bec2_diag_3d(3,donprod_idx_t)='mmol /m2/s'
      vname_bec2_diag_3d(4,donprod_idx_t)=' '
      wrt_bec2_diag_3d(donprod_idx_t)=.false.

      donremin_idx_t=48
      vname_bec2_diag_3d(1,donremin_idx_t)='DON_REMIN'
      vname_bec2_diag_3d(2,donremin_idx_t)='Remineralization rate of DON'
      vname_bec2_diag_3d(3,donremin_idx_t)='mmol N/m2/s  '
      vname_bec2_diag_3d(4,donremin_idx_t)=' '
      wrt_bec2_diag_3d(donremin_idx_t)=.false.

      dofeprod_idx_t=49
      vname_bec2_diag_3d(1,dofeprod_idx_t)='DOFE_PROD'
      vname_bec2_diag_3d(2,dofeprod_idx_t)='Production of dissolved organic iron'
      vname_bec2_diag_3d(3,dofeprod_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,dofeprod_idx_t)=' '
      wrt_bec2_diag_3d(dofeprod_idx_t)=.false.

      doferemin_idx_t=50
      vname_bec2_diag_3d(1,doferemin_idx_t)='DOFE_REMIN'
      vname_bec2_diag_3d(2,doferemin_idx_t)='Remineralization rate of DOFE'
      vname_bec2_diag_3d(3,doferemin_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,doferemin_idx_t)=' '
      wrt_bec2_diag_3d(doferemin_idx_t)=.false.

      dopprod_idx_t=51
      vname_bec2_diag_3d(1,dopprod_idx_t)='DOP_PROD'
      vname_bec2_diag_3d(2,dopprod_idx_t)='Production of dissolved organic P'
      vname_bec2_diag_3d(3,dopprod_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,dopprod_idx_t)=' '
      wrt_bec2_diag_3d(dopprod_idx_t)=.false.

      dopremin_idx_t=52
      vname_bec2_diag_3d(1,dopremin_idx_t)='DOP_REMIN'
      vname_bec2_diag_3d(2,dopremin_idx_t)='Remineralization rate of DOP'
      vname_bec2_diag_3d(3,dopremin_idx_t)='mmol/m2/s '
      vname_bec2_diag_3d(4,dopremin_idx_t)=' '
      wrt_bec2_diag_3d(dopremin_idx_t)=.false.

      diatsiuptake_idx_t=53
      vname_bec2_diag_3d(1,diatsiuptake_idx_t)='DIAT_SI_UPTAKE'
      vname_bec2_diag_3d(2,diatsiuptake_idx_t)='Silicon uptake rate by diatoms'
      vname_bec2_diag_3d(3,diatsiuptake_idx_t)='mmol Si/m2/s  '
      vname_bec2_diag_3d(4,diatsiuptake_idx_t)=' '
      wrt_bec2_diag_3d(diatsiuptake_idx_t)=.false.

      ironuptakesp_idx_t=54
      vname_bec2_diag_3d(1,ironuptakesp_idx_t)='PHOTOFE_SP'
      vname_bec2_diag_3d(2,ironuptakesp_idx_t)='Iron uptake rate by small phytoplankton'
      vname_bec2_diag_3d(3,ironuptakesp_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,ironuptakesp_idx_t)=' '
      wrt_bec2_diag_3d(ironuptakesp_idx_t)=.false.

      ironuptakediat_idx_t=55
      vname_bec2_diag_3d(1,ironuptakediat_idx_t)='PHOTOFE_DIAT'
      vname_bec2_diag_3d(2,ironuptakediat_idx_t)='Iron uptake rate by diatoms'
      vname_bec2_diag_3d(3,ironuptakediat_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,ironuptakediat_idx_t)=' '
      wrt_bec2_diag_3d(ironuptakediat_idx_t)=.false.

      ironuptakediaz_idx_t=56
      vname_bec2_diag_3d(1,ironuptakediaz_idx_t)='PHOTOFE_DIAZ'
      vname_bec2_diag_3d(2,ironuptakediaz_idx_t)='Iron uptake rate by diatotrophs'
      vname_bec2_diag_3d(3,ironuptakediaz_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,ironuptakediaz_idx_t)=' '
      wrt_bec2_diag_3d(ironuptakediaz_idx_t)=.false.

      nitrif_idx_t=57
      vname_bec2_diag_3d(1,nitrif_idx_t)='NITRIF'
      vname_bec2_diag_3d(2,nitrif_idx_t)='Nitrification'
      vname_bec2_diag_3d(3,nitrif_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,nitrif_idx_t)=' '
      wrt_bec2_diag_3d(nitrif_idx_t)=.false.

      denitrif_idx_t=58
      vname_bec2_diag_3d(1,denitrif_idx_t)='DENITRIF'
      vname_bec2_diag_3d(2,denitrif_idx_t)='Denitrification'
      vname_bec2_diag_3d(3,denitrif_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,denitrif_idx_t)=' '
      wrt_bec2_diag_3d(denitrif_idx_t)=.false.

      spno3uptake_idx_t=59
      vname_bec2_diag_3d(1,spno3uptake_idx_t)='SP_NO3_UPTAKE'
      vname_bec2_diag_3d(2,spno3uptake_idx_t)='NO3 uptake rate of small phyto'
      vname_bec2_diag_3d(3,spno3uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,spno3uptake_idx_t)=' '
      wrt_bec2_diag_3d(spno3uptake_idx_t)=.true.

      diatno3uptake_idx_t=60
      vname_bec2_diag_3d(1,diatno3uptake_idx_t)='DIAT_NO3_UPTAKE'
      vname_bec2_diag_3d(2,diatno3uptake_idx_t)='NO3 uptake rate of diatoms'
      vname_bec2_diag_3d(3,diatno3uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,diatno3uptake_idx_t)=' '
      wrt_bec2_diag_3d(diatno3uptake_idx_t)=.true.

      diazno3uptake_idx_t=61
      vname_bec2_diag_3d(1,diazno3uptake_idx_t)='DIAZ_NO3_UPTAKE'
      vname_bec2_diag_3d(2,diazno3uptake_idx_t)='NO3 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,diazno3uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,diazno3uptake_idx_t)=' '
      wrt_bec2_diag_3d(diazno3uptake_idx_t)=.true.

      spnh4uptake_idx_t=62
      vname_bec2_diag_3d(1,spnh4uptake_idx_t)='SP_NH4_UPTAKE'
      vname_bec2_diag_3d(2,spnh4uptake_idx_t)='NH4 uptake rate of small phyto'
      vname_bec2_diag_3d(3,spnh4uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,spnh4uptake_idx_t)=' '
      wrt_bec2_diag_3d(spnh4uptake_idx_t)=.false.

      diatnh4uptake_idx_t=63
      vname_bec2_diag_3d(1,diatnh4uptake_idx_t)='DIAT_NH4_UPTAKE'
      vname_bec2_diag_3d(2,diatnh4uptake_idx_t)='NH4 uptake rate of diatoms'
      vname_bec2_diag_3d(3,diatnh4uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,diatnh4uptake_idx_t)=' '
      wrt_bec2_diag_3d(diatnh4uptake_idx_t)=.false.

      diaznh4uptake_idx_t=64
      vname_bec2_diag_3d(1,diaznh4uptake_idx_t)='DIAZ_NH4_UPTAKE'
      vname_bec2_diag_3d(2,diaznh4uptake_idx_t)='NH4 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,diaznh4uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,diaznh4uptake_idx_t)=' '
      wrt_bec2_diag_3d(diaznh4uptake_idx_t)=.false.

      grazedicsp_idx_t=65
      vname_bec2_diag_3d(1,grazedicsp_idx_t)='SP_GRAZE_DIC'
      vname_bec2_diag_3d(2,grazedicsp_idx_t)='Small phyto grazing rate routed to DIC'
      vname_bec2_diag_3d(3,grazedicsp_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,grazedicsp_idx_t)=' '
      wrt_bec2_diag_3d(grazedicsp_idx_t)=.false.

      grazedicdiat_idx_t=66
      vname_bec2_diag_3d(1,grazedicdiat_idx_t)='DIAT_GRAZE_DIC'
      vname_bec2_diag_3d(2,grazedicdiat_idx_t)='Diatom grazing rate routed to DIC'
      vname_bec2_diag_3d(3,grazedicdiat_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,grazedicdiat_idx_t)=' '
      wrt_bec2_diag_3d(grazedicdiat_idx_t)=.false.

      grazedicdiaz_idx_t=67
      vname_bec2_diag_3d(1,grazedicdiaz_idx_t)='DIAZ_GRAZE_DIC'
      vname_bec2_diag_3d(2,grazedicdiaz_idx_t)='Diazotroph grazing rate routed to DIC'
      vname_bec2_diag_3d(3,grazedicdiaz_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,grazedicdiaz_idx_t)=' '
      wrt_bec2_diag_3d(grazedicdiaz_idx_t)=.false.

      lossdicsp_idx_t=68
      vname_bec2_diag_3d(1,lossdicsp_idx_t)='SP_LOSS_DIC'
      vname_bec2_diag_3d(2,lossdicsp_idx_t)='Small phyto non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,lossdicsp_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,lossdicsp_idx_t)=' '
      wrt_bec2_diag_3d(lossdicsp_idx_t)=.false.

      lossdicdiat_idx_t=69
      vname_bec2_diag_3d(1,lossdicdiat_idx_t)='DIAT_LOSS_DIC'
      vname_bec2_diag_3d(2,lossdicdiat_idx_t)='Diatom non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,lossdicdiat_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,lossdicdiat_idx_t)=' '
      wrt_bec2_diag_3d(lossdicdiat_idx_t)=.false.

      lossdicdiaz_idx_t=70
      vname_bec2_diag_3d(1,lossdicdiaz_idx_t)='DIAZ_LOSS_DIC'
      vname_bec2_diag_3d(2,lossdicdiaz_idx_t)='Diazotroph non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,lossdicdiaz_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,lossdicdiaz_idx_t)=' '
      wrt_bec2_diag_3d(lossdicdiaz_idx_t)=.false.

      zoolossdic_idx_t=71
      vname_bec2_diag_3d(1,zoolossdic_idx_t)='ZOO_LOSS_DIC'
      vname_bec2_diag_3d(2,zoolossdic_idx_t)='Zooplankton loss routed to DIC'
      vname_bec2_diag_3d(3,zoolossdic_idx_t)='mmol C/m2/s  '
      vname_bec2_diag_3d(4,zoolossdic_idx_t)=' '
      wrt_bec2_diag_3d(zoolossdic_idx_t)=.false.

      diazagg_idx_t=72
      vname_bec2_diag_3d(1,diazagg_idx_t)='DIAZ_AGG'
      vname_bec2_diag_3d(2,diazagg_idx_t)='Aggregation of diatoms'
      vname_bec2_diag_3d(3,diazagg_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,diazagg_idx_t)=' '
      wrt_bec2_diag_3d(diazagg_idx_t)=.false.

      grazespzoo_idx_t=73
      vname_bec2_diag_3d(1,grazespzoo_idx_t)='GRAZE_SP_ZOO'
      vname_bec2_diag_3d(2,grazespzoo_idx_t)='Grazing rate on small phytoplankton routed to zoo'
      vname_bec2_diag_3d(3,grazespzoo_idx_t)='mmol C/m2/s '
      vname_bec2_diag_3d(4,grazespzoo_idx_t)=' '
      wrt_bec2_diag_3d(grazespzoo_idx_t)=.false.

      grazediatzoo_idx_t=74
      vname_bec2_diag_3d(1,grazediatzoo_idx_t)='GRAZE_DIAT_ZOO'
      vname_bec2_diag_3d(2,grazediatzoo_idx_t)='Grazing rate on diatoms routed to zoo'
      vname_bec2_diag_3d(3,grazediatzoo_idx_t)='mmol C/m2/s'
      vname_bec2_diag_3d(4,grazediatzoo_idx_t)=' '
      wrt_bec2_diag_3d(grazediatzoo_idx_t)=.false.

      grazediazzoo_idx_t=75
      vname_bec2_diag_3d(1,grazediazzoo_idx_t)='GRAZE_DIAZ_ZOO'
      vname_bec2_diag_3d(2,grazediazzoo_idx_t)='Grazing rate on diazotrophs routed to zoo'
      vname_bec2_diag_3d(3,grazediazzoo_idx_t)='mmol C/m2/s'
      vname_bec2_diag_3d(4,grazediazzoo_idx_t)=' '
      wrt_bec2_diag_3d(grazediazzoo_idx_t)=.false.

      spqcaco3_idx_t=76
      vname_bec2_diag_3d(1,spqcaco3_idx_t)='QCACO3'
      vname_bec2_diag_3d(2,spqcaco3_idx_t)='Small phyto CaCO3/C ratio'
      vname_bec2_diag_3d(3,spqcaco3_idx_t)='mmol CaCO3/mmol C'
      vname_bec2_diag_3d(4,spqcaco3_idx_t)=' '
      wrt_bec2_diag_3d(spqcaco3_idx_t)=.false.

      spphotoacc_idx_t=77
      vname_bec2_diag_3d(1,spphotoacc_idx_t)='PHOTOACC_SP'
      vname_bec2_diag_3d(2,spphotoacc_idx_t)='Chl synthesis in photoadaptation for small phyto'
      vname_bec2_diag_3d(3,spphotoacc_idx_t)='mg Chl/m2/s'
      vname_bec2_diag_3d(4,spphotoacc_idx_t)=' '
      wrt_bec2_diag_3d(spphotoacc_idx_t)=.false.

      diatphotoacc_idx_t=78
      vname_bec2_diag_3d(1,diatphotoacc_idx_t)='PHOTOACC_DIAT'
      vname_bec2_diag_3d(2,diatphotoacc_idx_t)='Chl synthesis in photoadaptation for diatoms'
      vname_bec2_diag_3d(3,diatphotoacc_idx_t)='mg Chl/m2/s'
      vname_bec2_diag_3d(4,diatphotoacc_idx_t)=' '
      wrt_bec2_diag_3d(diatphotoacc_idx_t)=.false.

      diazphotoacc_idx_t=79
      vname_bec2_diag_3d(1,diazphotoacc_idx_t)='PHOTOACC_DIAZ'
      vname_bec2_diag_3d(2,diazphotoacc_idx_t)='Chl synthesis in photoadaptation for diazotrophs'
      vname_bec2_diag_3d(3,diazphotoacc_idx_t)='mg Chl/m2/s'
      vname_bec2_diag_3d(4,diazphotoacc_idx_t)=' '
      wrt_bec2_diag_3d(diazphotoacc_idx_t)=.false.

      spczero_idx_t=80
      vname_bec2_diag_3d(1,spczero_idx_t)='SPC_ZERO'
      vname_bec2_diag_3d(2,spczero_idx_t)='Change caused by setting negative SPC values to zero'
      vname_bec2_diag_3d(3,spczero_idx_t)='mMol C/m3/s'
      vname_bec2_diag_3d(4,spczero_idx_t)=' '
      wrt_bec2_diag_3d(spczero_idx_t)=.false.

      diatczero_idx_t=81
      vname_bec2_diag_3d(1,diatczero_idx_t)='DIATC_ZERO'
      vname_bec2_diag_3d(2,diatczero_idx_t)='Change caused by setting negative DIATC values to zero'
      vname_bec2_diag_3d(3,diatczero_idx_t)='mMol C/m3/s'
      vname_bec2_diag_3d(4,diatczero_idx_t)=' '
      wrt_bec2_diag_3d(diatczero_idx_t)=.false.

      diazczero_idx_t=82
      vname_bec2_diag_3d(1,diazczero_idx_t)='DIAZC_ZERO'
      vname_bec2_diag_3d(2,diazczero_idx_t)='Change caused by setting negative DIAZC values to zero'
      vname_bec2_diag_3d(3,diazczero_idx_t)='mMol C/m3/s'
      vname_bec2_diag_3d(4,diazczero_idx_t)=' '
      wrt_bec2_diag_3d(diazczero_idx_t)=.false.

      doczero_idx_t=83
      vname_bec2_diag_3d(1,doczero_idx_t)='DOC_ZERO'
      vname_bec2_diag_3d(2,doczero_idx_t)='Change caused by setting negative DOC values to zero'
      vname_bec2_diag_3d(3,doczero_idx_t)='mMol C/m3/s'
      vname_bec2_diag_3d(4,doczero_idx_t)=' '
      wrt_bec2_diag_3d(doczero_idx_t)=.false.

      zooczero_idx_t=84
      vname_bec2_diag_3d(1,zooczero_idx_t)='ZOOC_ZERO'
      vname_bec2_diag_3d(2,zooczero_idx_t)='Change caused by setting negative ZOOC values to zero'
      vname_bec2_diag_3d(3,zooczero_idx_t)='mMol C/m3/s'
      vname_bec2_diag_3d(4,zooczero_idx_t)=' '
      wrt_bec2_diag_3d(zooczero_idx_t)=.false.

      spcaco3zero_idx_t=85
      vname_bec2_diag_3d(1,spcaco3zero_idx_t)='SPCACO3_ZERO'
      vname_bec2_diag_3d(2,spcaco3zero_idx_t)='Change caused by setting negative SPCACO3 values to zero'
      vname_bec2_diag_3d(3,spcaco3zero_idx_t)='mMol CaCO3/m3/s'
      vname_bec2_diag_3d(4,spcaco3zero_idx_t)=' '
      wrt_bec2_diag_3d(spcaco3zero_idx_t)=.false.

      donrremin_idx_t=86
      vname_bec2_diag_3d(1,donrremin_idx_t) = 'DONr_REMIN'
      vname_bec2_diag_3d(2,donrremin_idx_t)='Portion of refractory DON remineralized'
      vname_bec2_diag_3d(3,donrremin_idx_t)='mMol N/m2/s'
      vname_bec2_diag_3d(4,donrremin_idx_t)=' '
      wrt_bec2_diag_3d(donrremin_idx_t)=.false.

      totchl_idx_t=87
      vname_bec2_diag_3d(1,totchl_idx_t)='TOT_CHL'
      vname_bec2_diag_3d(2,totchl_idx_t)='Total Chlorophyll'
      vname_bec2_diag_3d(3,totchl_idx_t)='mg Chl/m3'
      vname_bec2_diag_3d(4,totchl_idx_t)=' '
      wrt_bec2_diag_3d(totchl_idx_t)=.True.

      spplim_idx_t=88
      vname_bec2_diag_3d(1,spplim_idx_t)='SP_P_LIM'
      vname_bec2_diag_3d(2,spplim_idx_t)='Small phytoplankton P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,spplim_idx_t)=' '
      vname_bec2_diag_3d(4,spplim_idx_t)=' '
      wrt_bec2_diag_3d(spplim_idx_t)=.false.

      diatplim_idx_t=89
      vname_bec2_diag_3d(1,diatplim_idx_t)='DIAT_P_LIM'
      vname_bec2_diag_3d(2,diatplim_idx_t)='Diatom P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,diatplim_idx_t)=' '
      vname_bec2_diag_3d(4,diatplim_idx_t)=' '
      wrt_bec2_diag_3d(diatplim_idx_t)=.false.

      diazplim_idx_t=90
      vname_bec2_diag_3d(1,diazplim_idx_t)='DIAZ_P_LIM'
      vname_bec2_diag_3d(2,diazplim_idx_t)='Diazotroph P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,diazplim_idx_t)=' '
      vname_bec2_diag_3d(4,diazplim_idx_t)=' '
      wrt_bec2_diag_3d(diazplim_idx_t)=.false.

      totphytoc_idx_t=91
      vname_bec2_diag_3d(1,totphytoc_idx_t)='TOT_PHYTOC'
      vname_bec2_diag_3d(2,totphytoc_idx_t)='Total Phytoplankton Carbon'
      vname_bec2_diag_3d(3,totphytoc_idx_t)='mMol C/m3'
      vname_bec2_diag_3d(4,totphytoc_idx_t)=' '
      wrt_bec2_diag_3d(totphytoc_idx_t)=.True.

      o2cons_idx_t=92
      vname_bec2_diag_3d(1,o2cons_idx_t)='O2_CONSUMPTION'
      vname_bec2_diag_3d(2,o2cons_idx_t)='O2 consumption rate'
      vname_bec2_diag_3d(3,o2cons_idx_t)='mmol O2/m2/s'
      vname_bec2_diag_3d(4,o2cons_idx_t)=' '
      wrt_bec2_diag_3d(o2cons_idx_t)=.True.

      o2prod_idx_t=93
      vname_bec2_diag_3d(1,o2prod_idx_t)='O2_PRODUCTION'
      vname_bec2_diag_3d(2,o2prod_idx_t)='O2 production rate'
      vname_bec2_diag_3d(3,o2prod_idx_t)='mmol O2/m2/s'
      vname_bec2_diag_3d(4,o2prod_idx_t)=' '
      wrt_bec2_diag_3d(o2prod_idx_t)=.True.

      ammox_idx_t=94                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,ammox_idx_t)='AMMOX'
      vname_bec2_diag_3d(2,ammox_idx_t)='Ammonium oxidation (NH4 -> NO2)'
      vname_bec2_diag_3d(3,ammox_idx_t)='mmol N/m2/s '
      vname_bec2_diag_3d(4,ammox_idx_t)=' '
      wrt_bec2_diag_3d(ammox_idx_t)=.false.

      nitrox_idx_t=95                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,nitrox_idx_t)='NITROX'
      vname_bec2_diag_3d(2,nitrox_idx_t)='Nitrite oxidation (NO2 -> NO3)'
      vname_bec2_diag_3d(3,nitrox_idx_t)='mmol N/m2/s'
      vname_bec2_diag_3d(4,nitrox_idx_t)=' '
      wrt_bec2_diag_3d(nitrox_idx_t)=.true.

      anammox_idx_t=96                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,anammox_idx_t)='ANAMMOX'
      vname_bec2_diag_3d(2,anammox_idx_t)='Anammox (NO2+NH4 -> N2)'
      vname_bec2_diag_3d(3,anammox_idx_t)='mmol N/m2/s'
      vname_bec2_diag_3d(4,anammox_idx_t)=' '
      wrt_bec2_diag_3d(anammox_idx_t)=.true.

      denitrif1_idx_t=97                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,denitrif1_idx_t)='DENITRIF1'
      vname_bec2_diag_3d(2,denitrif1_idx_t)='Denitrification1'
      vname_bec2_diag_3d(3,denitrif1_idx_t)='mmol/m2/s'
      vname_bec2_diag_3d(4,denitrif1_idx_t)=' '
      wrt_bec2_diag_3d(denitrif1_idx_t)=.true.

      denitrif2_idx_t=98                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,denitrif2_idx_t)='DENITRIF2'
      vname_bec2_diag_3d(2,denitrif2_idx_t)='Denitrification2'
      vname_bec2_diag_3d(3,denitrif2_idx_t)='mmol/m2/s'
      vname_bec2_diag_3d(4,denitrif2_idx_t)=' '
      wrt_bec2_diag_3d(denitrif2_idx_t)=.false.

      denitrif3_idx_t=99                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,denitrif3_idx_t)='DENITRIF3'
      vname_bec2_diag_3d(2,denitrif3_idx_t)='Denitrification3'
      vname_bec2_diag_3d(3,denitrif3_idx_t)='mmol/m2/s'
      vname_bec2_diag_3d(4,denitrif3_idx_t)=' '
      wrt_bec2_diag_3d(denitrif3_idx_t)=.false.

      spno2uptake_idx_t=100                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,spno2uptake_idx_t)='SP_NO2_UPTAKE'
      vname_bec2_diag_3d(2,spno2uptake_idx_t)='NO2 uptake rate of small phyto'
      vname_bec2_diag_3d(3,spno2uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,spno2uptake_idx_t)=' '
      wrt_bec2_diag_3d(spno2uptake_idx_t)=.false.

      diatno2uptake_idx_t=101                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,diatno2uptake_idx_t)='DIAT_NO2_UPTAKE'
      vname_bec2_diag_3d(2,diatno2uptake_idx_t)='NO3 uptake rate of diatoms'
      vname_bec2_diag_3d(3,diatno2uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,diatno2uptake_idx_t)=' '
      wrt_bec2_diag_3d(diatno2uptake_idx_t)=.false.

      diazno2uptake_idx_t=102                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,diazno2uptake_idx_t)='DIAZ_NO2_UPTAKE'
      vname_bec2_diag_3d(2,diazno2uptake_idx_t)='NO3 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,diazno2uptake_idx_t)='mmol/m2/s  '
      vname_bec2_diag_3d(4,diazno2uptake_idx_t)=' '
      wrt_bec2_diag_3d(diazno2uptake_idx_t)=.false.

      n2oammox_idx_t=103                 ! Available with NCYCLE_SY
      vname_bec2_diag_3d(1,n2oammox_idx_t)='N2OAMMOX'
      vname_bec2_diag_3d(2,n2oammox_idx_t)='N2O produced by Ammox'
      vname_bec2_diag_3d(3,n2oammox_idx_t)='mmol/m2/s'
      vname_bec2_diag_3d(4,n2oammox_idx_t)=''
      wrt_bec2_diag_3d(n2oammox_idx_t)=.false.

#endif
#endif
