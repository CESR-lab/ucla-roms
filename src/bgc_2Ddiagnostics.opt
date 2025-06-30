!
! - BGC DIAGNOSTICS
!
! nr_bec2_diag_2d : number of available 2d diags set in bgc_ecosys_vars

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

!!!!!!!!!!!! bec2_diag_2d !!!!!!!!!!!!!!

      pco2air_idx_t=1
      vname_bec2_diag_2d(1,pco2air_idx_t)='pCO2air'
      vname_bec2_diag_2d(2,pco2air_idx_t)='Surface water pCO2'
      vname_bec2_diag_2d(3,pco2air_idx_t)='ppm'
      vname_bec2_diag_2d(4,pco2air_idx_t)='  '
      wrt_bec2_diag_2d(pco2air_idx_t)=.false.     ! pCO2air

      parinc_idx_t=2
      vname_bec2_diag_2d(1,parinc_idx_t)='PARinc'
      vname_bec2_diag_2d(2,parinc_idx_t)='Incoming photosynth. active radiation '
      vname_bec2_diag_2d(3,parinc_idx_t)='W m-2'
      vname_bec2_diag_2d(4,parinc_idx_t)='  '
      wrt_bec2_diag_2d(parinc_idx_t)=.true.

      fgo2_idx_t=3
      vname_bec2_diag_2d(1,fgo2_idx_t)='FG_O2'
      vname_bec2_diag_2d(2,fgo2_idx_t)='Air-sea flux of O2'
      vname_bec2_diag_2d(3,fgo2_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,fgo2_idx_t)=' '
      wrt_bec2_diag_2d(fgo2_idx_t)=.true.

      fgco2_idx_t=4
      vname_bec2_diag_2d(1,fgco2_idx_t)='FG_CO2'
      vname_bec2_diag_2d(2,fgco2_idx_t)='Air-sea flux of CO2'
      vname_bec2_diag_2d(3,fgco2_idx_t)='mmol/m2/s '
      vname_bec2_diag_2d(4,fgco2_idx_t)=' '
      wrt_bec2_diag_2d(fgco2_idx_t)=.true.

      ws10m_idx_t=5
      vname_bec2_diag_2d(1,ws10m_idx_t)='WS10m'
      vname_bec2_diag_2d(2,ws10m_idx_t)='Wind speed at 10m'
      vname_bec2_diag_2d(3,ws10m_idx_t)='m/s '
      vname_bec2_diag_2d(4,ws10m_idx_t)=' '
      wrt_bec2_diag_2d(ws10m_idx_t)=.false.

      xkw_idx_t=6
      vname_bec2_diag_2d(1,xkw_idx_t)='XKW'
      vname_bec2_diag_2d(2,xkw_idx_t)='XKW'
      vname_bec2_diag_2d(3,xkw_idx_t)='m/s '
      vname_bec2_diag_2d(4,xkw_idx_t)=' '
      wrt_bec2_diag_2d(xkw_idx_t)=.false.

      atmpress_idx_t=7
      vname_bec2_diag_2d(1,atmpress_idx_t)='ATM_PRESS'
      vname_bec2_diag_2d(2,atmpress_idx_t)='Atmospheric pressure'
      vname_bec2_diag_2d(3,atmpress_idx_t)='atm '
      vname_bec2_diag_2d(4,atmpress_idx_t)=' '
      wrt_bec2_diag_2d(atmpress_idx_t)=.false.

      schmidto2_idx_t=8
      vname_bec2_diag_2d(1,schmidto2_idx_t)='SCHMIDT_O2'
      vname_bec2_diag_2d(2,schmidto2_idx_t)='Schmidt number for O2'
      vname_bec2_diag_2d(3,schmidto2_idx_t)=' '
      vname_bec2_diag_2d(4,schmidto2_idx_t)=' '
      wrt_bec2_diag_2d(schmidto2_idx_t)=.false.

      o2sat_idx_t=9
      vname_bec2_diag_2d(1,o2sat_idx_t)='O2SAT'
      vname_bec2_diag_2d(2,o2sat_idx_t)='O2 saturation concentration'
      vname_bec2_diag_2d(3,o2sat_idx_t)='mmol/m3 '
      vname_bec2_diag_2d(4,o2sat_idx_t)=' '
      wrt_bec2_diag_2d(o2sat_idx_t)=.false.

      schmidtco2_idx_t=10
      vname_bec2_diag_2d(1,schmidtco2_idx_t)='SCHMIDT_CO2'
      vname_bec2_diag_2d(2,schmidtco2_idx_t)='Schmidt number for CO2'
      vname_bec2_diag_2d(3,schmidtco2_idx_t)=' '
      vname_bec2_diag_2d(4,schmidtco2_idx_t)=' '
      wrt_bec2_diag_2d(schmidtco2_idx_t)=.false.

      pvo2_idx_t=11
      vname_bec2_diag_2d(1,pvo2_idx_t)='PV_O2'
      vname_bec2_diag_2d(2,pvo2_idx_t)='Piston velocity for O2'
      vname_bec2_diag_2d(3,pvo2_idx_t)='m/s '
      vname_bec2_diag_2d(4,pvo2_idx_t)=' '
      wrt_bec2_diag_2d(pvo2_idx_t)=.false.

      pvco2_idx_t=12
      vname_bec2_diag_2d(1,pvco2_idx_t)='PV_CO2'
      vname_bec2_diag_2d(2,pvco2_idx_t)='Piston velocity for CO2'
      vname_bec2_diag_2d(3,pvco2_idx_t)='m/s '
      vname_bec2_diag_2d(4,pvco2_idx_t)=' '
      wrt_bec2_diag_2d(pvco2_idx_t)=.false.

      ironflux_idx_t=13
      vname_bec2_diag_2d(1,ironflux_idx_t)='IRON_FLUX'
      vname_bec2_diag_2d(2,ironflux_idx_t)='Iron surface flux'
      vname_bec2_diag_2d(3,ironflux_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,ironflux_idx_t)=' '
      wrt_bec2_diag_2d(ironflux_idx_t)=.false.

      seddenitrif_idx_t=14
      vname_bec2_diag_2d(1,seddenitrif_idx_t)='SED_DENITRIF'
      vname_bec2_diag_2d(2,seddenitrif_idx_t)='Sediment denitrification'
      vname_bec2_diag_2d(3,seddenitrif_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,seddenitrif_idx_t)=' '
      wrt_bec2_diag_2d(seddenitrif_idx_t)=.true.

      ph_idx_t=15
      vname_bec2_diag_2d(1,ph_idx_t)='PH'
      vname_bec2_diag_2d(2,ph_idx_t)='surf pH value'
      vname_bec2_diag_2d(3,ph_idx_t)=' '
      vname_bec2_diag_2d(4,ph_idx_t)=' '
      wrt_bec2_diag_2d(ph_idx_t)=.false.

      pco2_idx_t=16
      vname_bec2_diag_2d(1,pco2_idx_t)='pCO2'
      vname_bec2_diag_2d(2,pco2_idx_t)='Surface water pCO2'
      vname_bec2_diag_2d(3,pco2_idx_t)='ppm'
      vname_bec2_diag_2d(4,pco2_idx_t)='  '
      wrt_bec2_diag_2d(pco2_idx_t)=.false.

      co2star_idx_t=17
      vname_bec2_diag_2d(1,co2star_idx_t)='CO2STAR'
      vname_bec2_diag_2d(2,co2star_idx_t)='CO2STAR'
      vname_bec2_diag_2d(3,co2star_idx_t)='mmol/m3 '
      vname_bec2_diag_2d(4,co2star_idx_t)=' '
      wrt_bec2_diag_2d(co2star_idx_t)=.false.

      pco2oc_idx_t=18
      vname_bec2_diag_2d(1,pco2oc_idx_t)='PCO2OC'
      vname_bec2_diag_2d(2,pco2oc_idx_t)='PCO2OC'
      vname_bec2_diag_2d(3,pco2oc_idx_t)='not looked up yet '
      vname_bec2_diag_2d(4,pco2oc_idx_t)=' '
      wrt_bec2_diag_2d(pco2oc_idx_t)=.false.

      dco2star_idx_t=19
      vname_bec2_diag_2d(1,dco2star_idx_t)='DCO2STAR'
      vname_bec2_diag_2d(2,dco2star_idx_t)='DCO2STAR'
      vname_bec2_diag_2d(3,dco2star_idx_t)='mmol/m3'
      vname_bec2_diag_2d(4,dco2star_idx_t)=' '
      wrt_bec2_diag_2d(dco2star_idx_t)=.false.

      fesedflux_idx_t=20
      vname_bec2_diag_2d(1,fesedflux_idx_t)='SED_FE_FLUX'
      vname_bec2_diag_2d(2,fesedflux_idx_t)='Benthic Iron Flux'
      vname_bec2_diag_2d(3,fesedflux_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,fesedflux_idx_t)=' '
      wrt_bec2_diag_2d(fesedflux_idx_t)=.false.

      fluxtosed_idx_t=21
      vname_bec2_diag_2d(1,fluxtosed_idx_t)='FLUX_TO_SED'
      vname_bec2_diag_2d(2,fluxtosed_idx_t)='POC reaching the sediments'
      vname_bec2_diag_2d(3,fluxtosed_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,fluxtosed_idx_t)=' '
      wrt_bec2_diag_2d(fluxtosed_idx_t)=.false.

      caco3fluxtosed_idx_t=22
      vname_bec2_diag_2d(1,caco3fluxtosed_idx_t)='CACO3_FLUX_TO_SED'
      vname_bec2_diag_2d(2,caco3fluxtosed_idx_t)='CaCO3 reaching the sed'
      vname_bec2_diag_2d(3,caco3fluxtosed_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,caco3fluxtosed_idx_t)=' '
      wrt_bec2_diag_2d(caco3fluxtosed_idx_t)=.false.

      sio2fluxtosed_idx_t=23
      vname_bec2_diag_2d(1,sio2fluxtosed_idx_t)='SIO2_FLUX_TO_SED'
      vname_bec2_diag_2d(2,sio2fluxtosed_idx_t)='Opal reaching the sed'
      vname_bec2_diag_2d(3,sio2fluxtosed_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,sio2fluxtosed_idx_t)=' '
      wrt_bec2_diag_2d(sio2fluxtosed_idx_t)=.false.

      pironfluxtosed_idx_t=24
      vname_bec2_diag_2d(1,pironfluxtosed_idx_t)='PIRON_FLUX_TO_SED'
      vname_bec2_diag_2d(2,pironfluxtosed_idx_t)='P iron reaching the sed'
      vname_bec2_diag_2d(3,pironfluxtosed_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,pironfluxtosed_idx_t)=' '
      wrt_bec2_diag_2d(pironfluxtosed_idx_t)=.false.

      dustfluxtosed_idx_t=25
      vname_bec2_diag_2d(1,dustfluxtosed_idx_t)='DUST_FLUX_TO_SED'
      vname_bec2_diag_2d(2,dustfluxtosed_idx_t)='Dust reaching the sed'
      vname_bec2_diag_2d(3,dustfluxtosed_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,dustfluxtosed_idx_t)=' '
      wrt_bec2_diag_2d(dustfluxtosed_idx_t)=.false.

      pocsedloss_idx_t=26
      vname_bec2_diag_2d(1,pocsedloss_idx_t)='POC_SED_LOSS'
      vname_bec2_diag_2d(2,pocsedloss_idx_t)='POC lost to sed at bottom'
      vname_bec2_diag_2d(3,pocsedloss_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,pocsedloss_idx_t)=' '
      wrt_bec2_diag_2d(pocsedloss_idx_t)=.false.

      otherremin_idx_t=27
      vname_bec2_diag_2d(1,otherremin_idx_t)='OTHER_REMIN'
      vname_bec2_diag_2d(2,otherremin_idx_t)='other remin pathways in sed'
      vname_bec2_diag_2d(3,otherremin_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,otherremin_idx_t)=' '
      wrt_bec2_diag_2d(otherremin_idx_t)=.false.

      caco3sedloss_idx_t=28
      vname_bec2_diag_2d(1,caco3sedloss_idx_t)='P_CACO3_SED_LOSS'
      vname_bec2_diag_2d(2,caco3sedloss_idx_t)='CaCO3 lost to sed'
      vname_bec2_diag_2d(3,caco3sedloss_idx_t)='mmol/m2/s'
      vname_bec2_diag_2d(4,caco3sedloss_idx_t)=' '
      wrt_bec2_diag_2d(caco3sedloss_idx_t)=.false.

      sio2sedloss_idx_t=29
      vname_bec2_diag_2d(1,sio2sedloss_idx_t)='SIO2_SED_LOSS'
      vname_bec2_diag_2d(2,sio2sedloss_idx_t)='SiO2 lost to sed'
      vname_bec2_diag_2d(3,sio2sedloss_idx_t)='mmol Si/m2/s'
      vname_bec2_diag_2d(4,sio2sedloss_idx_t)=' '
      wrt_bec2_diag_2d(sio2sedloss_idx_t)=.false.

      fgn2o_idx_t=30                 ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,fgn2o_idx_t)='FG_N2O'
      vname_bec2_diag_2d(2,fgn2o_idx_t)='Air-sea flux of N2O '
      vname_bec2_diag_2d(3,fgn2o_idx_t)='mmol/m2/s '
      vname_bec2_diag_2d(4,fgn2o_idx_t)=' '
      wrt_bec2_diag_2d(fgn2o_idx_t)=.True.

      schmidt_n2o_idx_t=31           ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,schmidt_n2o_idx_t)='SCHMIDT_N2O'
      vname_bec2_diag_2d(2,schmidt_n2o_idx_t)='Schmidt number for N2O'
      vname_bec2_diag_2d(3,schmidt_n2o_idx_t)=' '
      vname_bec2_diag_2d(4,schmidt_n2o_idx_t)=' '
      wrt_bec2_diag_2d(schmidt_n2o_idx_t)=.false.

      n2osat_idx_t=32                ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,n2osat_idx_t)='N2OSAT'
      vname_bec2_diag_2d(2,n2osat_idx_t)='N2O saturation concentration'
      vname_bec2_diag_2d(3,n2osat_idx_t)='mmol/m3 '
      vname_bec2_diag_2d(4,n2osat_idx_t)=' '
      wrt_bec2_diag_2d(n2osat_idx_t)=.false.

      pvn2o_idx_t=33                 ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,pvn2o_idx_t)='PV_N2O'
      vname_bec2_diag_2d(2,pvn2o_idx_t)='Piston velocity for N2O'
      vname_bec2_diag_2d(3,pvn2o_idx_t)='m/s '
      vname_bec2_diag_2d(4,pvn2o_idx_t)=' '
      wrt_bec2_diag_2d(pvn2o_idx_t)=.false.

      fgn2_idx_t=34                  ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,fgn2_idx_t)='FG_N2'
      vname_bec2_diag_2d(2,fgn2_idx_t)='Air-sea flux of N2 '
      vname_bec2_diag_2d(3,fgn2_idx_t)='mmol/m2/s '
      vname_bec2_diag_2d(4,fgn2_idx_t)=' '
      wrt_bec2_diag_2d(fgn2_idx_t)=.true.

      schmidt_n2_idx_t=35            ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,schmidt_n2_idx_t)='SCHMIDT_N2'
      vname_bec2_diag_2d(2,schmidt_n2_idx_t)='Schmidt number for N2'
      vname_bec2_diag_2d(3,schmidt_n2_idx_t)=' '
      vname_bec2_diag_2d(4,schmidt_n2_idx_t)=' '
      wrt_bec2_diag_2d(schmidt_n2_idx_t)=.false.

      n2sat_idx_t=36                 ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,n2sat_idx_t)='N2SAT'
      vname_bec2_diag_2d(2,n2sat_idx_t)='N2 saturation concentration'
      vname_bec2_diag_2d(3,n2sat_idx_t)='mmol/m3 '
      vname_bec2_diag_2d(4,n2sat_idx_t)=' '
      wrt_bec2_diag_2d(n2sat_idx_t)=.false.

      pvn2_idx_t=37                 ! Available with NCYCLE_SY
      vname_bec2_diag_2d(1,pvn2_idx_t)='PV_N2'
      vname_bec2_diag_2d(2,pvn2_idx_t)='Piston velocity for N2'
      vname_bec2_diag_2d(3,pvn2_idx_t)='m/s '
      vname_bec2_diag_2d(4,pvn2_idx_t)=' '
      wrt_bec2_diag_2d(pvn2_idx_t)=.false.




#endif
#endif
