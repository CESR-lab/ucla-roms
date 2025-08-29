      module bec2_vars
#include "cppdefs.opt"
#ifdef BIOLOGY_BEC2

      use param
      use tracers               ! for iPO4, etc, indices of bgc tracers (formerly in param.h)
      use roms_read_write


      implicit none
      ! DTRACER_MODULE: stores the tendencies of the biological tracers
      ! in the long run this should be removed, no need to copy arrays?
      ! local copies are made in ecosys_set_interior already anyway?
      real,allocatable,dimension(:,:,:,:) :: DTRACER_MODULE

C      integer, pointer :: nr_bec2_diag_2d, nr_bec2_diag_3d, nr_bec2_diag
C      integer, pointer :: nr_bec2_wrdiag_2d, nr_bec2_wrdiag_3d, nr_bec2_wrdiag

      real,dimension(:,:,:,:),pointer   :: bec2_diag_3d
      real,dimension(:,:,:)  ,pointer   :: bec2_diag_2d

! Control BEC2_DIAG Output vars
      logical, dimension(:), pointer :: wrt_bec2_diag_2d
      logical, dimension(:), pointer :: wrt_bec2_diag_3d
      integer, dimension(:), pointer :: idx_bec2_diag_2d
      integer, dimension(:), pointer :: idx_bec2_diag_3d

      character*72, pointer ::  vname_bec2_diag_2d(:,:)
      character*72, pointer ::  vname_bec2_diag_3d(:,:)

#if defined(BEC2_DIAG)
      ! Indices to be used in bec2_diag_3d only:
      integer :: par_idx_t, pocfluxin_idx_t,
     &   pocprod_idx_t, pocremin_idx_t, caco3fluxin_idx_t,
     &   pcaco3prod_idx_t, caco3remin_idx_t, sio2fluxin_idx_t,
     &   sio2prod_idx_t, sio2remin_idx_t, dustfluxin_idx_t,
     &   dustremin_idx_t, pironfluxin_idx_t, pironprod_idx_t,
     &   pironremin_idx_t, grazesp_idx_t, grazediat_idx_t,
     &   grazediaz_idx_t, sploss_idx_t, diatloss_idx_t,
     &   zooloss_idx_t, spagg_idx_t, diatagg_idx_t,
     &   photocsp_idx_t, photocdiat_idx_t, totprod_idx_t,
     &   docprod_idx_t, docremin_idx_t, fescavenge_idx_t,
     &   spnlim_idx_t, spfeuptake_idx_t, sppo4uptake_idx_t,
     &   splightlim_idx_t, diatnlim_idx_t, diatfeuptake_idx_t,
     &   diatpo4uptake_idx_t, diatsio3uptake_idx_t, diatlightlim_idx_t,
     &   caco3prod_idx_t, diaznfix_idx_t, diazloss_idx_t,
     &   photocdiaz_idx_t, diazpo4uptake_idx_t, diazfeuptake_idx_t,
     &   diazlightlim_idx_t, fescavengerate_idx_t, donprod_idx_t,
     &   donremin_idx_t, dofeprod_idx_t, doferemin_idx_t,
     &   dopprod_idx_t, dopremin_idx_t, diatsiuptake_idx_t,
     &   ironuptakesp_idx_t, ironuptakediat_idx_t, ironuptakediaz_idx_t,
     &   nitrif_idx_t, denitrif_idx_t, spno3uptake_idx_t,
     &   diatno3uptake_idx_t, diazno3uptake_idx_t, spnh4uptake_idx_t,
     &   diatnh4uptake_idx_t, diaznh4uptake_idx_t, grazedicsp_idx_t,
     &   grazedicdiat_idx_t, grazedicdiaz_idx_t, lossdicsp_idx_t,
     &   lossdicdiat_idx_t, lossdicdiaz_idx_t, zoolossdic_idx_t,
     &   diazagg_idx_t, grazespzoo_idx_t, grazediatzoo_idx_t,
     &   grazediazzoo_idx_t, spqcaco3_idx_t, spphotoacc_idx_t,
     &   diatphotoacc_idx_t, diazphotoacc_idx_t, spczero_idx_t,
     &   diatczero_idx_t, diazczero_idx_t, doczero_idx_t,
     &   zooczero_idx_t, spcaco3zero_idx_t, donrremin_idx_t,
     &   totchl_idx_t, spplim_idx_t, diatplim_idx_t, diazplim_idx_t,
     &   totphytoc_idx_t, o2cons_idx_t, o2prod_idx_t, ammox_idx_t,
     &   nitrox_idx_t, anammox_idx_t, denitrif1_idx_t, denitrif2_idx_t,
     &   denitrif3_idx_t, spno2uptake_idx_t,
     &   diatno2uptake_idx_t, diazno2uptake_idx_t, n2oammox_idx_t

      ! Indices to be used in bec2_diag_2d only:
      integer :: pco2air_idx_t, parinc_idx_t, fgo2_idx_t, fgco2_idx_t, ws10m_idx_t, xkw_idx_t,
     &   atmpress_idx_t, schmidto2_idx_t, o2sat_idx_t, schmidtco2_idx_t, pvo2_idx_t, pvco2_idx_t,
     &   ironflux_idx_t, seddenitrif_idx_t,ph_idx_t,pco2_idx_t, co2star_idx_t, pco2oc_idx_t,
     &   dco2star_idx_t, fesedflux_idx_t, fluxtosed_idx_t, caco3fluxtosed_idx_t, sio2fluxtosed_idx_t,
     &   pironfluxtosed_idx_t, dustfluxtosed_idx_t, pocsedloss_idx_t, otherremin_idx_t,
     &   caco3sedloss_idx_t, sio2sedloss_idx_t, schmidt_n2o_idx_t, pvn2o_idx_t, n2osat_idx_t,
     &   fgn2o_idx_t, schmidt_n2_idx_t, pvn2_idx_t, fgn2_idx_t, n2sat_idx_t

!      ! Array for storing the Netcdf variable IDs of the diagnostics:
!      ! The IDs of the 2d vars are first, the those of the 3d.
!      integer :: hisT_bec2_diag(nr_bec2_diag)  , avgT_bec2_diag(nr_bec2_diag),
!     &           slavgT_bec2_diag(nr_bec2_diag), rstT_bec2_diag(nr_bec2_diag)
      integer, parameter, public :: nr_bec2_diag_2d = 37
      integer, parameter, public :: nr_bec2_diag_3d = 103
C
#endif /* BEC2_DIAG */


      real,allocatable,dimension(:,:) :: ifrac, press

      logical,allocatable,dimension(:,:) :: landmask



      logical :: lsource_sink,lflux_gas_o2, lflux_gas_co2
#if defined Ncycle_SY
     &  ,lflux_gas_n2o, lflux_gas_n2
# endif
     &  ,liron_flux,ldust_flux


!
! Parameters related to sinking particles:
!
      real ::
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


      real,allocatable,dimension(:,:) ::
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

      real,allocatable,dimension(:,:) ::
     &   POC_remin, P_iron_remin, P_SiO2_remin, P_CaCO3_remin


!!
!! Arrays related to carbon chemistry:
!!
      real,allocatable,dimension(:,:) :: ph_srf

      ! misc from porting ETH code:
      character(len=40) tclm_name(NT) ! Taken from ETH code ncvars


      public init_scalars_bec2
      contains
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine init_arrays_bec2_vars  ![
      implicit none

      allocate( DTRACER_MODULE(GLOBAL_2D_ARRAY,N,ntrc_bio) )

      allocate( ifrac(GLOBAL_2D_ARRAY), press(GLOBAL_2D_ARRAY) )

      allocate( landmask(GLOBAL_2D_ARRAY) )

      allocate(P_CaCO3_sflux_out(GLOBAL_2D_ARRAY), P_CaCO3_hflux_out(GLOBAL_2D_ARRAY) )
      allocate( P_SiO2_sflux_out(GLOBAL_2D_ARRAY),  P_SiO2_hflux_out(GLOBAL_2D_ARRAY) )
      allocate(   dust_sflux_out(GLOBAL_2D_ARRAY),    dust_hflux_out(GLOBAL_2D_ARRAY) )
      allocate( P_iron_sflux_out(GLOBAL_2D_ARRAY),  P_iron_hflux_out(GLOBAL_2D_ARRAY) )
      allocate(    POC_sflux_out(GLOBAL_2D_ARRAY),     POC_hflux_out(GLOBAL_2D_ARRAY) )
      allocate( P_CaCO3_sflux_in(GLOBAL_2D_ARRAY),  P_CaCO3_hflux_in(GLOBAL_2D_ARRAY) )
      allocate(  P_SiO2_sflux_in(GLOBAL_2D_ARRAY),   P_SiO2_hflux_in(GLOBAL_2D_ARRAY) )
      allocate(    dust_sflux_in(GLOBAL_2D_ARRAY),     dust_hflux_in(GLOBAL_2D_ARRAY) )
      allocate(  P_iron_sflux_in(GLOBAL_2D_ARRAY),   P_iron_hflux_in(GLOBAL_2D_ARRAY) )
      allocate(     POC_sflux_in(GLOBAL_2D_ARRAY),      POC_hflux_in(GLOBAL_2D_ARRAY) )
      allocate( P_CaCO3_sed_loss(GLOBAL_2D_ARRAY),   P_SiO2_sed_loss(GLOBAL_2D_ARRAY) )
      allocate(  P_iron_sed_loss(GLOBAL_2D_ARRAY),      POC_sed_loss(GLOBAL_2D_ARRAY) )
      allocate(    dust_sed_loss(GLOBAL_2D_ARRAY) )
      allocate(        DOP_remin(GLOBAL_2D_ARRAY),        DOPr_remin(GLOBAL_2D_ARRAY) )
      allocate(        POC_remin(GLOBAL_2D_ARRAY),      P_iron_remin(GLOBAL_2D_ARRAY) )
      allocate(     P_SiO2_remin(GLOBAL_2D_ARRAY),     P_CaCO3_remin(GLOBAL_2D_ARRAY) )


      allocate( ph_srf(GLOBAL_2D_ARRAY) )

      end subroutine init_arrays_bec2_vars  !]

!----------------------------------------------------------------------
      subroutine init_scalars_bec2()  ![
!
! Set initial values for  globally accessible (stored in common
! blocks) scalar variables of the BEC model.
!
      use comm_vars
      use scalars

      implicit none
      integer i,j, itrc ,lvar,lenstr ! ierr, DevinD removed ierr as arguement
      integer omp_get_num_threads

#ifdef SOLVE3D

# ifdef BEC2_DIAG

# endif /* BEC2_DIAG */

      tclm_name(iPO4) = 'po4_time'
      tclm_name(iNO3) = 'no3_time'
      tclm_name(iSIO3) = 'sio3_time'
      tclm_name(iO2) = 'o2_time'
      tclm_name(iFE) = 'fe_time'
      tclm_name(iDIC) = 'dic_time'
      tclm_name(iALK) = 'alk_time'
      tclm_name(iNH4) = 'nh4_time'
      tclm_name(iDOC) = 'doc_time'
      tclm_name(iDOP) = 'dop_time'
      tclm_name(iDOP) = 'dop_time'
      tclm_name(iDOPR) = 'dop_time'
      tclm_name(iDON) = 'don_time'
      tclm_name(iDONR) = 'don_time'
      tclm_name(iDOFE) = 'dofe_time'
      tclm_name(iSPC) = 'spc_time'
      tclm_name(iSPCHL) = 'spchl_time'
      tclm_name(iSPCACO3) = 'spcaco3_time'
      tclm_name(iSPFE) = 'spfe_time'
      tclm_name(iDIATC) = 'diatc_time'
      tclm_name(iDIATCHL) = 'diatchl_time'
      tclm_name(iDIATSI) = 'diatsi_time'
      tclm_name(iDIATFE) = 'diatfe_time'
      tclm_name(iZOOC) = 'zooc_time'
      tclm_name(iDIAZC) = 'diazc_time'
      tclm_name(iDIAZCHL) = 'diazchl_time'
      tclm_name(iDIAZFE) = 'diazfe_time'
# ifdef Ncycle_SY
      tclm_name(iNO2) = 'no2_time'
      tclm_name(iN2) = 'n2_time'
      tclm_name(iN2O) = 'n2o_time'
# endif
#endif /* SOLVE3D */

      end subroutine init_scalars_bec2  !]
! ----------------------------------------------------------------------
C
#endif /*BIOLOGY_BEC2 */
      end module bec2_vars
