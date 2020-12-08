       real,dimension(GLOBAL_2D_ARRAY):: WS_AVG, XKW_AVG,
     &     AP_AVG, SCHMIDT_O2_AVG, O2SAT_AVG, FG_O2_AVG,
     &     PH_AVG, pCO2_AVG,
     &     pCO2air_AVG,
     &     SCHMIDT_CO2_AVG, DCO2STAR_AVG,
     &     CO2STAR_AVG,
     &     FG_CO2_AVG, IRON_FLUX_AVG, PARinc_avg,
     &     PARinc_flux_avg, zeta_bgc_flux_avg

#ifdef CH_CARBON_DEPTH
       real,dimension(GLOBAL_2D_ARRAY,N):: HCO3d_AVG, CO3d_AVG, CO2STARd_AVG,
     &     PHd_AVG

       common /time_avg1/
     &    HCO3d_AVG, CO3d_AVG, CO2STARd_AVG, PHd_AVG
#endif /*CH_CARBON_DEPTH*/

        real,dimension(GLOBAL_2D_ARRAY,N)::
     &    PAR_avg,
     &    PO4_RESTORE_AVG, NO3_RESTORE_AVG,
     &    SiO3_RESTORE_AVG, PAR_flux_avg, PO4STAR_AVG,
     &    POC_FLUX_IN_AVG, POC_PROD_AVG, POC_REMIN_AVG,
     &    CaCO3_FLUX_IN_AVG, CaCO3_PROD_AVG,
     &    CaCO3_REMIN_AVG,  SiO2_FLUX_IN_AVG,
     &    SiO2_PROD_AVG, SiO2_REMIN_AVG, dust_FLUX_IN_AVG,
     &    dust_REMIN_AVG, P_iron_FLUX_IN_AVG,
     &    P_iron_PROD_AVG, P_iron_REMIN_AVG,
     &    graze_sp_AVG, graze_diat_AVG, graze_tot_AVG,
     &    sp_loss_AVG, diat_loss_AVG, zoo_loss_AVG,
     &    sp_agg_AVG, diat_agg_AVG,
     &    photoC_sp_AVG, photoC_diat_AVG, tot_prod_AVG,
     &    DOC_prod_AVG, DOC_remin_AVG, Fe_scavenge_AVG,
     &    sp_N_lim_AVG, sp_Fe_lim_AVG, sp_PO4_lim_AVG,
     &    sp_light_lim_AVG, diat_N_lim_AVG, diat_Fe_lim_AVG,
     &    diat_PO4_lim_AVG, diat_SiO3_lim_AVG
        real,dimension(GLOBAL_2D_ARRAY,N)::
     &    diat_light_lim_AVG, CaCO3_form_AVG,
     &    diaz_Nfix_AVG, graze_diaz_AVG, diaz_loss_AVG,
     &     photoC_diaz_AVG, diaz_P_lim_AVG,
     &    diaz_Fe_lim_AVG, diaz_light_lim_AVG,
     &     Fe_scavenge_rate_AVG, DON_prod_AVG,
     &    DON_remin_AVG, DOFe_prod_AVG,
     &    DOFe_remin_AVG, DOP_prod_AVG,
     &    DOP_remin_AVG, bSI_form_AVG,
     &    photoFe_diaz_AVG, photoFe_diat_AVG,
#ifdef OXYLIM_BEC
     &    photoFe_sp_AVG, nitrif_AVG, denitr_DOC_AVG,
     &    denitr_POC_AVG
#else
     &    photoFe_sp_AVG,nitrif_AVG
#endif

       common /time_avg1/WS_AVG, XKW_AVG,
     &    AP_AVG, SCHMIDT_O2_AVG, O2SAT_AVG, FG_O2_AVG,
     &    PH_AVG, pCO2_AVG, pCO2air_AVG,
     &    SCHMIDT_CO2_AVG, CO2STAR_AVG, DCO2STAR_AVG,
     &    FG_CO2_AVG, IRON_FLUX_AVG, PARinc_avg,
     &    PARinc_flux_avg, zeta_bgc_flux_avg, PAR_avg,
     &    PO4_RESTORE_AVG, NO3_RESTORE_AVG,
     &    SiO3_RESTORE_AVG, PAR_flux_avg, PO4STAR_AVG,
     &    POC_FLUX_IN_AVG, POC_PROD_AVG, POC_REMIN_AVG,
     &    CaCO3_FLUX_IN_AVG, CaCO3_PROD_AVG,
     &    CaCO3_REMIN_AVG,  SiO2_FLUX_IN_AVG,
     &    SiO2_PROD_AVG, SiO2_REMIN_AVG, dust_FLUX_IN_AVG,
     &    dust_REMIN_AVG, P_iron_FLUX_IN_AVG,
     &    P_iron_PROD_AVG, P_iron_REMIN_AVG,
     &    graze_sp_AVG, graze_diat_AVG, graze_tot_AVG,
     &    sp_loss_AVG, diat_loss_AVG, zoo_loss_AVG


       common /time_avg2/
     &    sp_agg_AVG, diat_agg_AVG,
     &    photoC_sp_AVG, photoC_diat_AVG, tot_prod_AVG,
     &    DOC_prod_AVG, DOC_remin_AVG, Fe_scavenge_AVG,
     &    sp_N_lim_AVG, sp_Fe_lim_AVG, sp_PO4_lim_AVG,
     &    sp_light_lim_AVG, diat_N_lim_AVG, diat_Fe_lim_AVG,
     &    diat_PO4_lim_AVG, diat_SiO3_lim_AVG,
     &    diat_light_lim_AVG, CaCO3_form_AVG,
     &    diaz_Nfix_AVG, graze_diaz_AVG, diaz_loss_AVG,
     &     photoC_diaz_AVG, diaz_P_lim_AVG,
     &    diaz_Fe_lim_AVG, diaz_light_lim_AVG,
     &     Fe_scavenge_rate_AVG, DON_prod_AVG,
     &    DON_remin_AVG, DOFe_prod_AVG,
     &    DOFe_remin_AVG, DOP_prod_AVG,
     &    DOP_remin_AVG, bSI_form_AVG,
     &    photoFe_diaz_AVG, photoFe_diat_AVG,
#ifdef OXYLIM_BEC
     &    photoFe_sp_AVG, nitrif_AVG, denitr_DOC_AVG,
     &    denitr_POC_AVG
#else
     &    photoFe_sp_AVG, nitrif_AVG
#endif

#ifdef SLICE_AVG
      real PAR_slavg(GLOBAL_2D_ARRAY)
      real PARinc_slavg(GLOBAL_2D_ARRAY)
      real pco2_slavg(GLOBAL_2D_ARRAY)
      real pCO2air_slavg(GLOBAL_2D_ARRAY)
      real pH_slavg(GLOBAL_2D_ARRAY)
      common /time_slavg1/ PAR_slavg, PARinc_slavg,
     &        pco2_slavg, pCO2air_slavg, pH_slavg
#endif /* SLICE_AVG */
