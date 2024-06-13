      module bgc_ecosys_vars
      ! formerly bgc_ecosys_bec2.h

#include "cppdefs.opt"
#if defined(BIOLOGY_BEC2) || defined(MARBL)

      use param
      use tracers               ! for iPO4, etc, indices of bgc tracers (formerly in param.h)
      use roms_read_write
#ifdef MARBL      
      use marbl_driver,only: marbldrv_configure_saved_state
# ifdef MARBL_DIAGS      
     &     ,marbldrv_configure_diagnostics
     &     ,nr_marbl_diag_2d,nr_marbl_diag_3d
     &     ,diag_array_2d,diag_array_3d
     &     ,idx_bgc_diag_2d,idx_bgc_diag_3d
     &     ,wrt_bgc_diag_2d,wrt_bgc_diag_3d
# endif      
#endif      

      implicit none

#include "bgc.opt"

      ! tracer: values passed to biological model by ROMS
      ! DTRACER_MODULE: stores the tendencies of the biological tracers
      ! in the long run this should be removed, no need to copy arrays?
      ! local copies are made in ecosys_set_interior already anyway?
      real,allocatable,dimension(:,:,:,:) :: tracer
      real,allocatable,dimension(:,:,:,:) :: DTRACER_MODULE

#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
!
! Diagnostic variables appearing in average and history files:
!
      integer :: nr_bec2_diag_2d, nr_bec2_diag_3d, nr_bec2_diag
      integer :: nr_bec2_wrdiag_2d, nr_bec2_wrdiag_3d, nr_bec2_wrdiag
#ifndef MARBL_DIAGS
      parameter( nr_bec2_diag_3d=103, nr_bec2_diag_2d=37)
#else
      ! Just hardcoding this for now but would be nice to obtain it dynamically
      parameter( nr_bec2_diag_3d=nr_marbl_diag_3d,
     &           nr_bec2_diag_2d=nr_marbl_diag_2d)
#endif
      parameter( nr_bec2_diag=nr_bec2_diag_2d+nr_bec2_diag_3d )
      real,allocatable,dimension(:,:,:,:),target   :: bec2_diag_3d
      real,allocatable,dimension(:,:,:)  ,target   :: bec2_diag_2d

! Control BEC2_DIAG Output vars
      logical, target :: wrt_bec2_diag_2d(nr_bec2_diag_2d)
      logical, target :: wrt_bec2_diag_3d(nr_bec2_diag_3d)
      integer, target :: idx_bec2_diag_2d(nr_bec2_diag_2d)
      integer, target :: idx_bec2_diag_3d(nr_bec2_diag_3d) 
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
#endif /* BEC2_DIAG */
      ! Arrays storing information (name, unit, fill value) about each diagnostic variable:
      character*72 ::  vname_bec2_diag_2d(4,nr_bec2_diag_2d)
      character*72 ::  vname_bec2_diag_3d(4,nr_bec2_diag_3d)

      real,allocatable,dimension(:,:,:,:) :: bec2_diag_3d_avg
      real,allocatable,dimension(:,:,:)   :: bec2_diag_2d_avg

#endif /* BEC2_DIAG  || MARBL_DIAGS */

      real,allocatable,dimension(:,:) :: ifrac, press

      logical,allocatable,dimension(:,:) :: landmask

      
#ifdef BIOLOGY_BEC2
      logical :: lsource_sink,lflux_gas_o2, lflux_gas_co2
#if defined Ncycle_SY
     &  ,lflux_gas_n2o, lflux_gas_n2
# endif
     &  ,liron_flux,ldust_flux

!
! Relative tracer indices for prognostic variables:
!
      integer :: po4_ind_t, no3_ind_t, sio3_ind_t, nh4_ind_t, fe_ind_t, dic_ind_t, alk_ind_t,
     &        o2_ind_t, doc_ind_t, don_ind_t, dofe_ind_t, dop_ind_t, dopr_ind_t, donr_ind_t,
     &        zooc_ind_t, spchl_ind_t, spc_ind_t, spfe_ind_t, spcaco3_ind_t, diatchl_ind_t,
     &        diatc_ind_t, diatfe_ind_t, diatsi_ind_t, diazchl_ind_t, diazc_ind_t, diazfe_ind_t
#ifdef Ncycle_SY
      integer :: no2_ind_t, n2_ind_t, n2o_ind_t
# endif

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
#endif /*BIOLOGY_BEC2*/
!!
!! Arrays related to carbon chemistry: these are in bec2_diag_2d or
!! bec2_diag_3d if BEC2_DIAG is defined
!!
!#ifndef BEC2_DIAG
!      real,allocatable,dimension(:,:)   :: ph_hist, pCO2sw, PARinc
!      real,allocatable,dimension(:,:,:) :: PAR
!# ifndef PCO2AIR_FORCING
!!     otherwise defined in bgc_forces.h
!      real :: pco2air
!# endif
!
!      real,allocatable,dimension(:,:) :: ph_avg, pCO2_avg, pCO2air_avg, PARinc_avg
!      real,allocatable,dimension(:,:,:) :: PAR_avg
!#endif /* !BEC2_DIAG */

      contains

!----------------------------------------------------------------------
#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
      subroutine find_write_bgc_diag
      implicit none

      integer idiag
#if defined(BEC2_DIAG)
#include "bgc_2Ddiagnostics.opt" 
#include "bgc_3Ddiagnostics.opt"
      
#elif defined(MARBL_DIAGS)

      wrt_bgc_diag_2d => wrt_bec2_diag_2d
      wrt_bgc_diag_3d => wrt_bec2_diag_3d
      idx_bgc_diag_2d => idx_bec2_diag_2d
      idx_bgc_diag_3d => idx_bec2_diag_3d
      
      
      call marbldrv_configure_diagnostics(
     &     vname_bec2_diag_2d,vname_bec2_diag_3d
     & )
#endif /* BEC2_DIAG */ 
    
      nr_bec2_wrdiag_2d=0
      do idiag=1,nr_bec2_diag_2d
        if (wrt_bec2_diag_2d(idiag)) then
           nr_bec2_wrdiag_2d=nr_bec2_wrdiag_2d+1 
           idx_bec2_diag_2d(idiag)=nr_bec2_wrdiag_2d 
           if (mynode == 0) then
              write(*,'(7x,A,I3,I3,A,A)') '2D BGC diags :: ', 
     &        nr_bec2_wrdiag_2d, idiag, ' ', vname_bec2_diag_2d(1,idiag) 
           endif 
        else
           idx_bec2_diag_2d(idiag)=0
        endif
      enddo

      nr_bec2_wrdiag_3d=0
      do idiag=1,nr_bec2_diag_3d
        if (wrt_bec2_diag_3d(idiag)) then
           nr_bec2_wrdiag_3d=nr_bec2_wrdiag_3d+1
           idx_bec2_diag_3d(idiag)=nr_bec2_wrdiag_3d  
           if (mynode == 0) then
              write(*,'(7x,A,I3,I3,A,A)') '3D BGC diags :: ', 
     &        nr_bec2_wrdiag_3d, idiag, ' ',vname_bec2_diag_3d(1,idiag)
           endif
        else
           idx_bec2_diag_3d(idiag)=0
        endif
      enddo

      end subroutine find_write_bgc_diag  !]
#endif /* BEC2_DIAG || MARBL_DIAGS */

!----------------------------------------------------------------------
      subroutine init_arrays_bgc_ecosys_vars  ![
      implicit none

      allocate( tracer(GLOBAL_2D_ARRAY,N,ntrc_bio) )
      allocate( DTRACER_MODULE(GLOBAL_2D_ARRAY,N,ntrc_bio) )

#if defined(BEC2_DIAG) || defined (MARBL_DIAGS)

      call find_write_bgc_diag 
   
      allocate( bec2_diag_3d(GLOBAL_2D_ARRAY,N,nr_bec2_wrdiag_3d) )
      allocate( bec2_diag_2d(GLOBAL_2D_ARRAY,nr_bec2_wrdiag_2d) )
      
      bec2_diag_2d=0.0
      bec2_diag_3d=0.0
#ifdef MARBL_DIAGS
      diag_array_2d => bec2_diag_2d
      diag_array_3d => bec2_diag_3d
#endif      
      
      if (wrt_avg_dia) then
        allocate( bec2_diag_3d_avg(GLOBAL_2D_ARRAY,N,nr_bec2_wrdiag_3d) )
        allocate( bec2_diag_2d_avg(GLOBAL_2D_ARRAY,nr_bec2_wrdiag_2d) )
        bec2_diag_2d_avg=0.0
        bec2_diag_3d_avg=0.0
      endif

      if (mynode == 0) then
              write(*,'(7x,A,I3,I3)') 'BGC diags allocation :: ',
     &        nr_bec2_wrdiag_2d, nr_bec2_wrdiag_3d
      endif 

#endif /* BEC2_DIAG || MARBL_DIAGS */
#ifdef MARBL      
!     Allocate arrays for MARBL saved state variables
      call marbldrv_configure_saved_state
#endif 
      allocate( ifrac(GLOBAL_2D_ARRAY), press(GLOBAL_2D_ARRAY) )

      allocate( landmask(GLOBAL_2D_ARRAY) )
#if defined(BIOLOGY_BEC2)
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

!      allocate(
!     &   P_CaCO3_sflux_out, P_CaCO3_hflux_out,
!     &   P_SiO2_sflux_out, P_SiO2_hflux_out,
!     &   dust_sflux_out, dust_hflux_out,
!     &   P_iron_sflux_out, P_iron_hflux_out,
!     &   POC_sflux_out, POC_hflux_out,
!     &   P_CaCO3_sflux_in, P_CaCO3_hflux_in,
!     &   P_SiO2_sflux_in, P_SiO2_hflux_in,
!     &   dust_sflux_in, dust_hflux_in,
!     &   P_iron_sflux_in, P_iron_hflux_in,
!     &   POC_sflux_in, POC_hflux_in,
!     &   P_CaCO3_sed_loss, P_SiO2_sed_loss,
!     &   P_iron_sed_loss,POC_sed_loss,
!     &   dust_sed_loss,
!     &   DOP_remin, DOPr_remin,
!     &   mold=ifrac )                           ! DevinD added mold statement, to avoid
!                                                ! repeated GLOBAL.. for every array
!      allocate( POC_remin, P_iron_remin, P_SiO2_remin, P_CaCO3_remin, mold=ifrac )

      allocate( ph_srf(GLOBAL_2D_ARRAY) )

!#ifndef BEC2_DIAG
!      allocate( ph_hist, pCO2sw, PARinc, mold=ifrac )
!      allocate( PAR(GLOBAL_2D_ARRAY,N) )
!
!      allocate( ph_avg, pCO2_avg, pCO2air_avg, PARinc_avg, mold=ifrac )
!      allocate( PAR_avg(GLOBAL_2D_ARRAY,N) )
!#endif /* !BEC2_DIAG */
#endif /* BIOLOGY_BEC2*/
      end subroutine init_arrays_bgc_ecosys_vars  !]

!----------------------------------------------------------------------

#endif /* BIOLOGY_BEC2 || MARBL  */

      end module bgc_ecosys_vars
