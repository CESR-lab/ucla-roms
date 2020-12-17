!
!  tracer: values passed to biological model by ROMS
!  DTRACER_MODULE: stores the tendencies of the biological tracers
!
      real tracer(GLOBAL_2D_ARRAY,N,ntrc_bio)
      real DTRACER_MODULE(GLOBAL_2D_ARRAY,N,ntrc_bio)
      ! DevinD changed /tracers/ as clashes with module name 'tracers'
      common /tracers_com/ tracer, DTRACER_MODULE

      integer counter_no3(288),counter_coccochl(288)      
      common /counter_neg/ counter_no3, counter_coccochl

#ifdef BEC2_DIAG
!
! Diagnostic variables appearing in average and history files:
!
      integer nr_bec2_diag_2d, nr_bec2_diag_3d, nr_bec2_diag
      integer nr_cchem_mocsy_2d, nr_cchem_mocsy_3d
# ifdef CCHEM_MOCSY
      ! Parameters nr_cchem_mocsy_2d and nr_cchem_mocsy_3d give the numbers of *additional*
      ! diagnostic variables if MOCSY code is used for carbon chemistry (relative to the 
      ! OCMIP code):
#  ifdef CCHEM_TODEPTH
      parameter( nr_cchem_mocsy_2d=-3, nr_cchem_mocsy_3d=7 )
#  else
      parameter( nr_cchem_mocsy_2d=2, nr_cchem_mocsy_3d=0 )
#  endif /* CCHEM_TODEPTH */
# else /* CCHEM_MOCSY */
      parameter( nr_cchem_mocsy_2d=0, nr_cchem_mocsy_3d=0 )
# endif /* CCHEM_MOCSY */

/* Still in BEC2 Diag now start counting for nr_bec2_diag_3d */

      parameter( nr_bec2_diag_3d=93+nr_cchem_mocsy_3d
# ifdef Ncycle_SY
     &  +14
#endif
#ifdef N2O_NEV
     &  +2 
#endif
# ifdef USE_EXPLICIT_VSINK
     &  +10
#endif
# ifdef BEC_COCCO
     &  +18
#endif
     & , nr_bec2_diag_2d=29+nr_cchem_mocsy_2d 
# ifdef Ncycle_SY
     &  +8
# ifdef N2O_TRACER_DECOMP
     &  +4 ! 0 from coccos, 0 from impl sinking
#endif /* N2O_TRACER_DECOMP*/
#endif
#ifdef N2O_NEV
     &  +1 ! 0 from coccos, 0 from impl sinking
#endif /* N2O_NEV*/
     &  )
      parameter( nr_bec2_diag=nr_bec2_diag_2d+nr_bec2_diag_3d )
# ifdef BEC2_DIAG_USER
      real, pointer, dimension(:,:,:,:) :: bec2_diag_3d
      real, pointer, dimension(:,:,:) :: bec2_diag_2d
      logical bec2_diag_3d_l(nr_bec2_diag_3d), bec2_diag_2d_l(nr_bec2_diag_2d)
      integer nr_bec2_diag_2d_user, nr_bec2_diag_3d_user
      integer bec2_diag_3d_idx(nr_bec2_diag_3d), bec2_diag_2d_idx(nr_bec2_diag_2d)
# else
      real bec2_diag_3d(GLOBAL_2D_ARRAY,N,nr_bec2_diag_3d)
      real bec2_diag_2d(GLOBAL_2D_ARRAY,nr_bec2_diag_2d)
# endif /* BEC2_DIAG_USER */

      common /bec2_diag1/ bec2_diag_2d, bec2_diag_3d
# ifdef BEC2_DIAG_USER
     &      , bec2_diag_3d_l, bec2_diag_2d_l, bec2_diag_3d_idx, bec2_diag_2d_idx
     &      , nr_bec2_diag_2d_user, nr_bec2_diag_3d_user
# endif

      ! Indices to be used in bec2_diag_2d or bec2_diag_3d: these are 3d in space if
      ! MOCSY is selected and the C chemistry is computed at depth. Otherwise they are
      ! 2d in space.
!      integer ph_idx_t, pco2oc_idx_t, co2star_idx_t
# ifdef CCHEM_MOCSY
      ! Additional variables if if MOCSY is used instead of the OCMIP code for carbon
      ! chemistry. The numbers of these variables are given by nr_cchem_mocsy_2d and
      ! nr_cchem_mocsy_3d:
     &        , co3_idx_t, hco3_idx_t
#  ifdef CCHEM_TODEPTH
     &        , omega_calc_idx_t, omega_arag_idx_t
#  endif
# endif /* CCHEM_MOCSY */

  ! MF: Control BEC2_DIAG Output vars
      logical wrtavg_bec2_diag_2d(nr_bec2_diag_2d)
      logical wrtavg_bec2_diag_3d(nr_bec2_diag_3d)

      common wrtavg_bec2_diag_2d, wrtavg_bec2_diag_3d

      ! Indices to be used in bec2_diag_3d only:
      integer, parameter :: par_idx_t=1,pocfluxin_idx_t=par_idx_t+1,
     &   pocprod_idx_t=par_idx_t+2,pocremin_idx_t=par_idx_t+3,caco3fluxin_idx_t=par_idx_t+4,
     &   pcaco3prod_idx_t=par_idx_t+5,caco3remin_idx_t=par_idx_t+6,sio2fluxin_idx_t=par_idx_t+7,
     &   sio2prod_idx_t=par_idx_t+8,sio2remin_idx_t=par_idx_t+9,dustfluxin_idx_t=par_idx_t+10,
     &   dustremin_idx_t=par_idx_t+11,pironfluxin_idx_t=par_idx_t+12,pironprod_idx_t=par_idx_t+13,
     &   pironremin_idx_t=par_idx_t+14,grazesp_idx_t=par_idx_t+15,grazediat_idx_t=par_idx_t+16,
     &   grazediaz_idx_t=par_idx_t+17,sploss_idx_t=par_idx_t+18,diatloss_idx_t=par_idx_t+19,
     &   zooloss_idx_t=par_idx_t+20,spagg_idx_t=par_idx_t+21,diatagg_idx_t=par_idx_t+22,
     &   photocsp_idx_t=par_idx_t+23,photocdiat_idx_t=par_idx_t+24,totprod_idx_t=par_idx_t+25,
     &   docprod_idx_t=par_idx_t+26,docremin_idx_t=par_idx_t+27,fescavenge_idx_t=par_idx_t+28,
     &   spnlim_idx_t=par_idx_t+29,spfeuptake_idx_t=par_idx_t+30,sppo4uptake_idx_t=par_idx_t+31,
     &   splightlim_idx_t=par_idx_t+32,diatnlim_idx_t=par_idx_t+33,diatfeuptake_idx_t=par_idx_t+34,
     &   diatpo4uptake_idx_t=par_idx_t+35,diatsio3uptake_idx_t=par_idx_t+36,diatlightlim_idx_t=par_idx_t+37,
     &   caco3prod_idx_t=par_idx_t+38,diaznfix_idx_t=par_idx_t+39,diazloss_idx_t=par_idx_t+40,
     &   photocdiaz_idx_t=par_idx_t+41,diazpo4uptake_idx_t=par_idx_t+42,diazfeuptake_idx_t=par_idx_t+43,
     &   diazlightlim_idx_t=par_idx_t+44,fescavengerate_idx_t=par_idx_t+45,donprod_idx_t=par_idx_t+46,
     &   donremin_idx_t=par_idx_t+47,dofeprod_idx_t=par_idx_t+48,doferemin_idx_t=par_idx_t+49,
     &   dopprod_idx_t=par_idx_t+50,dopremin_idx_t=par_idx_t+51,diatsiuptake_idx_t=par_idx_t+52,
     &   ironuptakesp_idx_t=par_idx_t+53,ironuptakediat_idx_t=par_idx_t+54,ironuptakediaz_idx_t=par_idx_t+55,
     &   nitrif_idx_t=par_idx_t+56,denitrif_idx_t=par_idx_t+57,spno3uptake_idx_t=par_idx_t+58,
     &   diatno3uptake_idx_t=par_idx_t+59,diazno3uptake_idx_t=par_idx_t+60,spnh4uptake_idx_t=par_idx_t+61,
     &   diatnh4uptake_idx_t=par_idx_t+62,diaznh4uptake_idx_t=par_idx_t+63,grazedicsp_idx_t=par_idx_t+64,
     &   grazedicdiat_idx_t=par_idx_t+65,grazedicdiaz_idx_t=par_idx_t+66,lossdicsp_idx_t=par_idx_t+67,
     &   lossdicdiat_idx_t=par_idx_t+68,lossdicdiaz_idx_t=par_idx_t+69,zoolossdic_idx_t=par_idx_t+70,
     &   diazagg_idx_t=par_idx_t+71,grazespzoo_idx_t=par_idx_t+72,grazediatzoo_idx_t=par_idx_t+73,
     &   grazediazzoo_idx_t=par_idx_t+74,spqcaco3_idx_t=par_idx_t+75,spphotoacc_idx_t=par_idx_t+76,
     &   diatphotoacc_idx_t=par_idx_t+77,diazphotoacc_idx_t=par_idx_t+78,spczero_idx_t=par_idx_t+79,
     &   diatczero_idx_t=par_idx_t+80,diazczero_idx_t=par_idx_t+81,doczero_idx_t=par_idx_t+82,
     &   zooczero_idx_t=par_idx_t+83,spcaco3zero_idx_t=par_idx_t+84,donrremin_idx_t=par_idx_t+85,
     &   totchl_idx_t=par_idx_t+86,
     &   spplim_idx_t=par_idx_t+87,diatplim_idx_t=par_idx_t+88,diazplim_idx_t=par_idx_t+89,
     &   totphytoc_idx_t=par_idx_t+90,o2cons_idx_t=par_idx_t+91, o2prod_idx_t=par_idx_t+92
#  undef LAST_I
#  define LAST_I o2prod_idx_t
# ifdef USE_EXPLICIT_VSINK
     &   ,pironhardremin_idx_t=LAST_I+1, caco3hardremin_idx_t=LAST_I+2, sio2hardremin_idx_t=LAST_I+3
     &   ,pochardremin_idx_t=LAST_I+4, dusthardremin_idx_t=LAST_I+5
     &   ,pironsoftremin_idx_t=LAST_I+6, caco3softremin_idx_t=LAST_I+7, sio2softremin_idx_t=LAST_I+8
     &   ,pocsoftremin_idx_t=LAST_I+9, dustsoftremin_idx_t=LAST_I+10
# undef LAST_I
# define LAST_I dustsoftremin_idx_t
# else /* USE_EXPLICIT_VSINK */
! already defined above,no additional diagnostics
# endif /* USE_EXPLICIT_VSINK */
# if defined CCHEM_MOCSY && defined CCHEM_TODEPTH
     &   ,ph_idx_t=LAST_I+1, pco2oc_idx_t=ph_idx_t+1, co3_idx_t=ph_idx_t+2
     &   ,hco3_idx_t=ph_idx_t+3, co2star_idx_t=ph_idx_t+4
     &   ,omega_calc_idx_t=ph_idx_t+5, omega_arag_idx_t=ph_idx_t+6
#  undef LAST_I
#  define LAST_I omega_arag_idx_t
# endif
# ifdef BEC_COCCO
      integer, parameter :: grazecocco_idx_t=LAST_I+1,coccoloss_idx_t=LAST_I+2,
     &   coccoagg_idx_t=LAST_I+3,photoccocco_idx_t=LAST_I+4,
     &   cocconlim_idx_t=LAST_I+5,
     &   coccopo4uptake_idx_t=LAST_I+6,coccofeuptake_idx_t=LAST_I+7,
     &   coccolightlim_idx_t=LAST_I+8,
     &   caco3prodcocco_idx_t=LAST_I+9,ironuptakecocco_idx_t=LAST_I+10,
     &   coccono3uptake_idx_t=LAST_I+11,cocconh4uptake_idx_t=LAST_I+12,
     &   coccograzedic_idx_t=LAST_I+13,
     &   coccolossdic_idx_t=LAST_I+14,grazecoccozoo_idx_t=LAST_I+15,
     &   coccoqcaco3_idx_t=LAST_I+16,
     &   coccophotoacc_idx_t=LAST_I+17,
     &   coccoplim_idx_t=LAST_I+18
#  undef LAST_I
#  define LAST_I coccoplim_idx_t
# endif        

# ifdef Ncycle_SY
      integer, parameter :: ammox_idx_t=LAST_I+1,nitrox_idx_t=LAST_I+2,
     &   anammox_idx_t=LAST_I+3,denitrif1_idx_t=LAST_I+4, denitrif2_idx_t=LAST_I+5,
     &   denitrif3_idx_t=LAST_I+6,spno2uptake_idx_t=LAST_I+7,
     &   diatno2uptake_idx_t=LAST_I+8,diazno2uptake_idx_t=LAST_I+9,
     &   n2oammox_idx_t=LAST_I+10,n2osoden_cons_idx_t=LAST_I+11,
     &   n2oao1_cons_idx_t=LAST_I+12,n2oatm_cons_idx_t=LAST_I+13,
     &   n2osiden_cons_idx_t=LAST_I+14
#  undef LAST_I
#  define LAST_I n2osiden_cons_idx_t
# endif

# ifdef N2O_NEV
      integer, parameter :: n2oprodnev_idx_t=LAST_I+1,n2oconsnev_idx_t=LAST_I+2
# endif

      ! Indices to be used in bec2_diag_2d only:
      integer, parameter :: pco2air_idx_t=1,
     &   parinc_idx_t=pco2air_idx_t+1, fgo2_idx_t=pco2air_idx_t+2, fgco2_idx_t=pco2air_idx_t+3,
     &   ws10m_idx_t=pco2air_idx_t+4, xkw_idx_t=pco2air_idx_t+5, atmpress_idx_t=pco2air_idx_t+6,
     &   schmidto2_idx_t=pco2air_idx_t+7, o2sat_idx_t=pco2air_idx_t+8, schmidtco2_idx_t=pco2air_idx_t+9,
     &   pvo2_idx_t=pco2air_idx_t+10,pvco2_idx_t=pco2air_idx_t+11,ironflux_idx_t=pco2air_idx_t+12,
     &   seddenitrif_idx_t=pco2air_idx_t+13,ph_idx_t=pco2air_idx_t+14,pco2_idx_t=pco2air_idx_t+15,
     &   co2star_idx_t=pco2air_idx_t+16,pco2oc_idx_t=pco2air_idx_t+17,dco2star_idx_t=pco2air_idx_t+18
# undef LAST_I
# define LAST_I dco2star_idx_t
#ifdef Ncycle_SY
     &   ,schmidt_n2o_idx_t=LAST_I+1, pvn2o_idx_t=LAST_I+2, n2osat_idx_t=LAST_I+3,
     &    fgn2o_idx_t=LAST_I+4, schmidt_n2_idx_t=LAST_I+5, pvn2_idx_t=LAST_I+6, 
     &    fgn2_idx_t=LAST_I+7, n2sat_idx_t=LAST_I+8
# undef LAST_I
# define LAST_I n2sat_idx_t
#ifdef N2O_TRACER_DECOMP
     &   ,fgn2o_ao1_idx_t=LAST_I+1, fgn2o_siden_idx_t=LAST_I+2,
     &   fgn2o_soden_idx_t=LAST_I+3, fgn2o_atm_idx_t=LAST_I+4
# undef LAST_I
# define LAST_I fgn2o_atm_idx_t
# endif
# endif
#ifdef N2O_NEV
     &   ,fgn2o_nev_idx_t=LAST_I+1
# undef LAST_I
# define LAST_I fgn2o_nev_idx_t
# endif
# ifdef CCHEM_MOCSY
     &   ,ph_idx_t=pco2air_idx_t+16, pco2oc_idx_t=ph_idx_t+1, co3_idx_t=ph_idx_t+2
#  ifndef CCHEM_TODEPTH
     &   ,hco3_idx_t=ph_idx_t+3, co2star_idx_t=ph_idx_t+4
# undef LAST_I
# define LAST_I co2star_idx_t
#  endif
# endif /* CCHEM_MOCSY */
      integer, parameter :: fesedflux_idx_t=LAST_I+1,
     &   fluxtosed_idx_t=LAST_I+2,caco3fluxtosed_idx_t=LAST_I+3,
     &   sio2fluxtosed_idx_t=LAST_I+4,pironfluxtosed_idx_t=LAST_I+5,dustfluxtosed_idx_t=LAST_I+6,
     &   pocsedloss_idx_t=LAST_I+7,otherremin_idx_t=LAST_I+8,caco3sedloss_idx_t=LAST_I+9,
     &   sio2sedloss_idx_t=LAST_I+10 

      ! Array for storing the Netcdf variable IDs of the diagnostics:
      ! The IDs of the 2d vars are first, the those of the 3d.
      integer hisT_bec2_diag(nr_bec2_diag), avgT_bec2_diag(nr_bec2_diag), slavgT_bec2_diag(nr_bec2_diag),
     &        rstT_bec2_diag(nr_bec2_diag)

      ! Arrays storing information (name, unit, fill value) about each diagnostic variable:
      character*72  vname_bec2_diag_2d(4,nr_bec2_diag_2d)
      character*72  vname_bec2_diag_3d(4,nr_bec2_diag_3d)

      common /bec2_diag2/ hisT_bec2_diag, avgT_bec2_diag, slavgT_bec2_diag, rstT_bec2_diag,
     &   vname_bec2_diag_2d, vname_bec2_diag_3d

# ifdef BEC2_DIAG_USER
#  ifdef AVERAGES
      real, pointer, dimension(:,:,:,:) :: bec2_diag_3d_avg
      real, pointer, dimension(:,:,:) ::  bec2_diag_2d_avg
#   ifdef SLICE_AVG
      real, pointer, dimension(:,:,:) ::   bec2_diag_3d_slavg
      real, pointer, dimension(:,:,:) ::   bec2_diag_2d_slavg
#   endif
#  endif /* AVERAGES */
# else /* BEC2_DIAG_USER */
#  ifdef AVERAGES
      real bec2_diag_3d_avg(GLOBAL_2D_ARRAY,N,nr_bec2_diag_3d)
      real bec2_diag_2d_avg(GLOBAL_2D_ARRAY,nr_bec2_diag_2d)
#   ifdef SLICE_AVG
      real bec2_diag_3d_slavg(GLOBAL_2D_ARRAY,nr_bec2_diag_3d)
      real bec2_diag_2d_slavg(GLOBAL_2D_ARRAY,nr_bec2_diag_2d)
#   endif
#  endif /* AVERAGES */
# endif /* BEC2_DIAG_USER */
# ifdef AVERAGES
      common /bec2_diag3/ bec2_diag_3d_avg, bec2_diag_2d_avg
#  ifdef SLICE_AVG
     &    , bec2_diag_3d_slavg, bec2_diag_2d_slavg
#  endif
# endif /* AVERAGES */
!#else
!mm does not compile:
!      real PH_HIST(GLOBAL_2D_ARRAY)
!      common /ph_hist_com/ PH_HIST
#endif /* BEC2_DIAG */

!     IFRAC  sea ice fraction (non-dimensional)
!     PRESS  sea level atmospheric pressure (Pascals)
      real ifrac(GLOBAL_2D_ARRAY),
     &   press(GLOBAL_2D_ARRAY)
      common /fic_ap/ ifrac,press

      logical landmask(GLOBAL_2D_ARRAY)
      common /calcation/landmask

      logical lsource_sink,lflux_gas_o2, lflux_gas_co2
#if defined Ncycle_SY || defined N2O_NEV
     &  ,lflux_gas_n2o
#endif
#ifdef Ncycle_SY
     &  ,lflux_gas_n2
# endif
     &  ,liron_flux,ldust_flux
#ifdef RIVER_LOAD_N
     &  ,lriver_load_n
#endif
#ifdef RIVER_LOAD_P
     &  ,lriver_load_p
#endif
#ifdef RIVER_LOAD_ALK_DIC_SI
     &  ,lriver_load_alk,lriver_load_dic,lriver_load_si
#endif
      common /ecoflag/lsource_sink,lflux_gas_o2,lflux_gas_co2,
     &   liron_flux,ldust_flux
#ifdef RIVER_LOAD_N
     &  ,lriver_load_n
#endif
#ifdef RIVER_LOAD_P
     &  ,lriver_load_p
#endif
#ifdef RIVER_LOAD_ALK_DIC_SI
     &  , lriver_load_alk,lriver_load_dic,lriver_load_si
#endif

!
! Relative tracer indices for prognostic variables:
!
      integer po4_ind_t, no3_ind_t, sio3_ind_t, nh4_ind_t, fe_ind_t, dic_ind_t, alk_ind_t,
     &        o2_ind_t, doc_ind_t, don_ind_t, dofe_ind_t, dop_ind_t, dopr_ind_t, donr_ind_t,
     &        zooc_ind_t, spchl_ind_t, spc_ind_t, spfe_ind_t, spcaco3_ind_t, diatchl_ind_t,
     &        diatc_ind_t, diatfe_ind_t, diatsi_ind_t, diazchl_ind_t, diazc_ind_t, diazfe_ind_t
      parameter( po4_ind_t=1, no3_ind_t=2, sio3_ind_t=3, nh4_ind_t=4, fe_ind_t=5,
     &    o2_ind_t=6, dic_ind_t=7, alk_ind_t=8, doc_ind_t=9, don_ind_t=10, dofe_ind_t=11,
     &    dop_ind_t=12, dopr_ind_t=13, donr_ind_t=14, zooc_ind_t=15, spc_ind_t=16,
     &    spchl_ind_t=17, spfe_ind_t=18, spcaco3_ind_t=19, diatc_ind_t=20, diatchl_ind_t=21,
     &    diatfe_ind_t=22, diatsi_ind_t=23, diazc_ind_t=24, diazchl_ind_t=25, diazfe_ind_t=26
# undef LAST_I
# define LAST_I diazfe_ind_t
     &)
#ifdef BEC_COCCO
      integer, parameter ::
     &     coccoc_ind_t=27, coccochl_ind_t=28, coccocal_ind_t=29, coccofe_ind_t=30, cal_ind_t=31
#  undef LAST_I
#  define LAST_I cal_ind_t
#endif
#ifdef USE_EXPLICIT_VSINK
      integer, parameter ::
     &     dusthard_ind_t=LAST_I+1, pochard_ind_t=LAST_I+2, pcaco3hard_ind_t=LAST_I+3,
     &     psio2hard_ind_t=LAST_I+4, pironhard_ind_t=LAST_I+5,
     &     dustsoft_ind_t=LAST_I+6, pocsoft_ind_t=LAST_I+7, pcaco3soft_ind_t=LAST_I+8,
     &     psio2soft_ind_t=LAST_I+9, pironsoft_ind_t=LAST_I+10
#  undef LAST_I
#  define LAST_I pironsoft_ind_t
#endif
#ifdef Ncycle_SY
      integer, parameter ::
     &     no2_ind_t=LAST_I+1, n2_ind_t=LAST_I+2,  n2o_ind_t=LAST_I+3
#  undef LAST_I
#  define LAST_I n2o_ind_t
#ifdef N2O_TRACER_DECOMP
     &     ,n2o_ao1_ind_t=LAST_I+1, n2o_siden_ind_t=LAST_I+2, 
     &     n2o_soden_ind_t=LAST_I+3, n2o_atm_ind_t=LAST_I+4,
     &     n2_sed_ind_t=LAST_I+5
#  undef LAST_I
#  define LAST_I n2_sed_ind_t
# endif
# endif

# ifdef N2O_NEV
      integer, parameter ::
     &      n2o_nev_ind_t=LAST_I+1
#  undef LAST_I
#  define LAST_I n2o_nev_ind_t
#endif

!
! Parameters related to sinking particles:
!

      real 
     &   POC_diss,       ! diss. length (m), modified by TEMP
     &   POC_mass,       ! molecular weight of POC
     &   P_CaCO3_diss,   ! diss. length (m)
     &   P_CaCO3_gamma,  ! prod frac -> hard subclass
     &   P_CaCO3_mass,   ! molecular weight of CaCO
     &   P_CaCO3_rho,    ! QA mass ratio for CaCO3
     &   P_SiO2_diss,    ! diss. length (m), modified by TEMP
     &   P_SiO2_gamma,   ! prod frac -> hard subclass
     &   P_SiO2_mass,    ! molecular weight of SiO2
     &   P_SiO2_rho,     ! QA mass ratio for SiO2
     &   dust_diss,      ! diss. length (m)
     &   dust_gamma,     ! prod frac -> hard subclass
     &   dust_mass,      ! base units are already grams
     &   dust_rho,       ! QA mass ratio for dust
     &   P_iron_gamma,    ! prod frac -> hard subclass
     &   POC_gamma       ! prod frac -> hard subclass


      common /sinking_particles1/ POC_diss,POC_mass,P_CaCO3_diss,
     &   P_CaCO3_mass,P_CaCO3_rho,P_SiO2_diss,P_SiO2_mass,P_SiO2_rho,
     &   dust_diss,dust_mass,dust_rho

      ! The gamma parameters are set init_particulate_terms and so they
      ! need to be in a common block:
      common /sinking_particles1/ POC_gamma, P_CaCO3_gamma, P_SiO2_gamma,
     &   dust_gamma, P_iron_gamma
!
! Arrays related to sinking particles:
!
!  *_hflux_in: incoming flux of hard subclass (base units/m^2/sec)
!  *_hflux_out: outgoing flux of hard subclass (base units/m^2/sec)
!  *_sflux_in: incoming flux of soft subclass (base units/m^2/sec)
!  *_sflux_out: outgoing flux of soft subclass (base units/m^2/sec)
!  *_sed_loss: loss to sediments (base units/m^2/sec)
!  *_remin: remineralization term (base units/m^3/sec)

      real, dimension(GLOBAL_2D_ARRAY) ::
     &   P_CaCO3_sflux_out, P_CaCO3_hflux_out,
     &   P_SiO2_sflux_out, P_SiO2_hflux_out,
     &   dust_sflux_out, dust_hflux_out,
     &   P_iron_sflux_out, P_iron_hflux_out,
     &   POC_sflux_out, POC_hflux_out, 
     &   P_CaCO3_sflux_in, P_CaCO3_hflux_in,
     &   P_SiO2_sflux_in, P_SiO2_hflux_in, 
     &   dust_sflux_in, dust_hflux_in,
     &   P_iron_sflux_in, P_iron_hflux_in,
     &   POC_sflux_in, POC_hflux_in,
     &   P_CaCO3_sed_loss, P_SiO2_sed_loss, 
     &   P_iron_sed_loss,POC_sed_loss,
     &   dust_sed_loss,
     &   DOP_remin, DOPr_remin

      real, dimension(GLOBAL_2D_ARRAY) ::
#ifdef USE_EXPLICIT_VSINK
     &   dusthard_remin, POChard_remin, P_CaCO3hard_remin,
     &   P_SiO2hard_remin, P_ironhard_remin,
     &   dustsoft_remin, POCsoft_remin, P_CaCO3soft_remin,
     &   P_SiO2soft_remin, P_ironsoft_remin,
#endif
     &   POC_remin, P_iron_remin, P_SiO2_remin, P_CaCO3_remin

      common /sinking_particles2/
     &   P_CaCO3_sflux_out, P_CaCO3_hflux_out,
     &   P_SiO2_sflux_out, P_SiO2_hflux_out,
     &   dust_sflux_out, dust_hflux_out,
     &   P_iron_sflux_out, P_iron_hflux_out,
     &   POC_sflux_out, POC_hflux_out,
     &   P_CaCO3_sflux_in, P_CaCO3_hflux_in,
     &   P_SiO2_sflux_in, P_SiO2_hflux_in, 
     &   dust_sflux_in, dust_hflux_in,
     &   P_iron_sflux_in, P_iron_hflux_in,
     &   POC_sflux_in, POC_hflux_in,
     &   P_CaCO3_sed_loss, P_SiO2_sed_loss, 
     &   P_iron_sed_loss,POC_sed_loss,
     &   dust_sed_loss,
     &   DOP_remin, DOPr_remin,
#ifdef USE_EXPLICIT_VSINK
     &   dusthard_remin, POChard_remin, P_CaCO3hard_remin,
     &   P_SiO2hard_remin, P_ironhard_remin,
     &   dustsoft_remin, POCsoft_remin, P_CaCO3soft_remin,
     &   P_SiO2soft_remin, P_ironsoft_remin,
#endif
     &   POC_remin, P_iron_remin, P_SiO2_remin, P_CaCO3_remin

!
! VSinkFlux array and related indices:
!
#ifdef USE_EXPLICIT_VSINK
      ! Number of sinking components:
      ! Indices to be used in VSinkFlux only (!!)
      integer, parameter :: nsink=10,iDUSTHARD_VSink=1,iPOCHARD_VSink=iDUSTHARD_VSink+1,
     &   iPCACO3HARD_VSink=iDUSTHARD_VSink+2,iPSIO2HARD_VSink=iDUSTHARD_VSink+3,  
     &   iPIRONHARD_VSink=iDUSTHARD_VSink+4,iDUSTSOFT_VSink=iDUSTHARD_VSink+5,  
     &   iPOCSOFT_VSink=iDUSTHARD_VSink+6,iPCACO3SOFT_VSink=iDUSTHARD_VSink+7,  
     &   iPSIO2SOFT_VSink=iDUSTHARD_VSink+8,iPIRONSOFT_VSink=iDUSTHARD_VSink+9  
      ! Vertical sink fluxes [mmol m-2 s-1], upward flux is positive
      real VSinkFlux(GLOBAL_2D_ARRAY,0:N,nsink)
      ! Indices of sinking variables in temprorary array tracer (used in ecosys_bec2.F):
!      integer, dimension(nsink), parameter :: 
!     &    tidx_vsink = (/ dust_ind_t, poc_ind_t, 
!     &        pcaco3_ind_t, psio2_ind_t, piron_ind_t /)
      common /bec2_vsink/ VSinkFlux
#endif

!
! Arrays related to carbon chemistry: these are in bec2_diag_2d or
! bec2_diag_3d if BEC2_DIAG is defined
!
#ifndef BEC2_DIAG
      real, dimension(GLOBAL_2D_ARRAY) ::
     &   ph_hist, pCO2sw, PARinc
      real, dimension(GLOBAL_2D_ARRAY,N) ::
     &   PAR
      common /c_chem/ ph_hist, pCO2sw, PARinc, PAR
# ifndef PCO2AIR_FORCING
!     otherwise defined in bgc_forces.h
      real pco2air
      common /c_chem/pco2air
# endif

      real, dimension(GLOBAL_2D_ARRAY) ::
     &   ph_avg, pCO2_avg, pCO2air_avg, PARinc_avg
      real, dimension(GLOBAL_2D_ARRAY,N) ::
     &   PAR_avg
      common /time_avg1/ PAR_avg, PARinc_avg, 
     &        pco2_avg, pCO2air_avg, pH_avg
# ifdef SLICE_AVG
      real PAR_slavg(GLOBAL_2D_ARRAY)
      real PARinc_slavg(GLOBAL_2D_ARRAY)
      real pco2_slavg(GLOBAL_2D_ARRAY)
      real pCO2air_slavg(GLOBAL_2D_ARRAY)
      real pH_slavg(GLOBAL_2D_ARRAY)
      common /time_slavg1/ PAR_slavg, PARinc_slavg, 
     &        pco2_slavg, pCO2air_slavg, pH_slavg
# endif /* SLICE_AVG */
#endif /* !BEC2_DIAG */

!
! Options related to MOCSY:
!
#ifdef CCHEM_MOCSY
      character*10 optcon, optt, optp, optb, optkf, optk1k2, optgas
      common /mocsy_opt/  optcon, optt, optp, optb, optkf, optk1k2, optgas
#endif

