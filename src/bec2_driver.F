#include "cppdefs.opt"

      ! DevinD notes: subroutines that take in a slice of arguements, e.g. landmask in cschmidt_n20
      ! is not good because it makes a local copy everytime the routine is called.
      ! Better to send in full array. Same for O2satu routine.

#ifdef BIOLOGY_BEC2

       subroutine ecosys_bec2_tile(Istr,Iend,Jstr,Jend)
       ! Taken from PDamien's ecosys_bec2_td2.F

       use param
       use tracers, only: t, iTandS
       use surf_flux ! , only: sustr,sustr_r,svstr,svstr_r,srflx
# ifdef BULK_FRC
       use bulk_frc, only: wspd_used ! to compute gas exchange
# endif
       use bgc_io
       use grid
       use ocean_vars
       use scalars
       use bgc_shared_vars
       use bec2_vars
       use bgc_forces
       use bec2_params

       implicit none


       integer i,j,k,m
       integer Istr,Iend,Jstr,Jend
       integer su_i1, su_i2, sv_j1, sv_j2
       integer idiag

       real
     &    O2SAT_USED      ! used O2 saturation (mmol/m^3)
#if defined Ncycle_SY
     &    , N2OSAT_USED   ! used O2 saturation (mmol/m^3)
     &    , N2SAT_USED    ! used N2 saturation (mmol/m^3)
# endif

       real,dimension(istr:iend,jstr:jend)::
     &    SSTT,         ! sea surface temperature (C)
     &    SSSS,         ! sea surface salinity (psu)
     &    SHF_QSW,      ! penetrative solar heat flux (W/m^2)
     &    QA_dust_def,  ! incoming deficit in the QA(dust) POC flux
     &    PAR_out,      ! photosynthetically available radiation (W/m^2)
     &    FICE_USED,    ! used ice fraction (non-dimensional)
     &    WS_USED,      ! used wind speed (m/s) ! HF: changed from (cm/s)
     &    XKW,          ! a * WS ** 2 (m/s) ! HF: changed from (cm/s)
     &    AP_USED,      ! used atmospheric pressure (atm)
     &    XKW_ICE,      ! common portion of piston vel., a*(1-fice)*u**2 (m/s) ! HF: changed from (cm/s)
     &    SCHMIDT_USED, ! used Schmidt number
     &    PV,           ! piston velocity (m/s) ! HF: changed from (cm/s)
     &    O2SAT_1atm,   ! O2 saturation @ 1 atm (mmol/m^3)
     &    FLUX          ! tracer flux (mmol/m^2/s)
#if defined DAILYPAR_PHOTOINHIBITION
     &   ,SHF_QSW_AVG   ! daily avg penetrative solar heat flux (W/m^2)
     &   ,PAR_out_AVG   ! daily avg photosynthetically available radiation (W/m^2)
#endif
# ifndef BULK_FRC
     &   ,sustr_rho     ! wind stress in xi direction on rho grid (N/m^2)
     &   ,svstr_rho     ! wind stress in xi direction on rho grid (N/m^2)
# endif

       real STF(istr:iend,jstr:jend,ntrc_bio)

       real,dimension(istr:iend)::
     &    PHLO          ! lower bound for ph in solver
     &   ,PHHI          ! upper bound for ph in solver
       real, PARAMETER ::
     &    phlo_init = 7.0,  ! low bound for ph for no prev soln
     &    phhi_init = 9.0,  ! high bound for ph for no prev soln
     &    del_ph = 0.20     ! delta-ph for prev soln

       real, PARAMETER ::
     &    a = 8.6e-7       ! a = 0.31 cm/hr s^2/m^2 in (s/m)
       real,dimension(istr:iend)::
     &    PH_NEW,       ! computed PH from solver
     &    XCO2,         ! atmospheric co2 conc. (dry-air, 1 atm)
     &    CO2STAR_ROW,  ! CO2STAR from solver
     &    DCO2STAR_ROW, ! DCO2STAR from solver
     &    pco2oc_ROW,   ! pco2oc from solver
     &    DpCO2_ROW     ! DpCO2 from solver

# if defined Ncycle_SY
       real,dimension(istr:iend,jstr:jend)::
     &    N2OSAT_1atm   ! N2O saturation @ 1 atm (mmol/m^3)

       real,dimension(istr:iend)::
     &    XN2O         ! atmospheric n2o conc. (dry-air, 1 atm)
       real,dimension(istr:iend,jstr:jend)::
     &    N2SAT_1atm   ! N2 saturation @ 1 atm (mmol/m^3)
# endif

      if(mynode==0) write(*,'(7x,A)') 'ecosys_bec2_td:: computing bgc' ! DevinD added

# ifndef BULK_FRC
!
! Interpolate/extrapolate wind stress to rho grid:
!
       do j=jstr,jend
          if (j == jstr) then
             sv_j1 = j+1
          else
             sv_j1 = j
          endif
          if (j .lt. jend) then
             sv_j2 = j+1
          else
             sv_j2 = j-1
          end if
          do i=istr,iend
             if (i == istr) then
                su_i1 = i+1
             else
                su_i1 = i
             end if
             if (i .lt. iend) then
                su_i2 = i+1
             else
                su_i2 = i-1
             end if
             ! DevinD: should change to sustr_r from surf_flux module, have already
             sustr_rho(i,j) = 0.5*(sustr(su_i1,j)+sustr(su_i2,j))
             svstr_rho(i,j) = 0.5*(svstr(i,sv_j1)+svstr(i,sv_j2))
          end do
       end do
# endif

!
! Initialize local copies of tracer array:
!


       sstt(istr:iend,jstr:jend) = t(istr:iend,jstr:jend,N,nnew,1)
       ssss(istr:iend,jstr:jend) = t(istr:iend,jstr:jend,N,nnew,2)

!
! Get short wave radiation:
!
#ifdef DAILYPAR_BEC
       shf_qsw(istr:iend,jstr:jend) =
     &               swrad_avg(istr:iend,jstr:jend)*rho0*Cp
#elif defined DAILYPAR_PHOTOINHIBITION
       shf_qsw_avg(istr:iend,jstr:jend) =
     &               swrad_avg(istr:iend,jstr:jend)*rho0*Cp
       shf_qsw(istr:iend,jstr:jend)=srflx(istr:iend,jstr:jend)*rho0*Cp
#else
       shf_qsw(istr:iend,jstr:jend)=srflx(istr:iend,jstr:jend)*rho0*Cp
#endif

!
! Gas exchange:
!
#ifdef ICEOBS
# define IFRAC sic /* Ice fraction has different name in this case */
#endif
       FICE_USED(istr:iend,jstr:jend) = IFRAC(istr:iend,jstr:jend)

# ifdef BULK_FRC
       ws_used(istr:iend,jstr:jend) = wspd_used(istr:iend,jstr:jend)
# else
       call WS(sustr_rho(istr:iend,jstr:jend),
     &    svstr_rho(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &    ws_used(istr:iend,jstr:jend), istr,iend,jstr,jend)
# endif

       XKW(istr:iend,jstr:jend) = a * WS_USED(istr:iend,jstr:jend)*WS_USED(istr:iend,jstr:jend)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(ws10m_idx_t)) then
          idiag = idx_bec2_diag_2d(ws10m_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = WS_USED(istr:iend,jstr:jend)
       endif
# endif

     !---------------------------------------------------------------------
     !   convert PRESS from Pascals to atm,no, it is atm already
     !----------------------------------------------------------------------
       AP_USED(istr:iend,jstr:jend) = PRESS(istr:iend,jstr:jend)

     !----------------------------------------------------------------------
     !   Set bad AP values to 1. This is necessary for runs restarting off
     !   a run in which the flux coupler didnt restart on AP correctly.
     !----------------------------------------------------------------------
       WHERE (AP_USED(istr:iend,jstr:jend) > 1.5 .OR. AP_USED(istr:iend,jstr:jend) < 0.5)
          AP_USED(istr:iend,jstr:jend) = c1
       END WHERE

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(xkw_idx_t)) then
          idiag = idx_bec2_diag_2d(xkw_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = XKW(istr:iend,jstr:jend)
       endif
       if (wrt_bec2_diag_2d(atmpress_idx_t)) then
          idiag = idx_bec2_diag_2d(atmpress_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = AP_USED(istr:iend,jstr:jend)
       endif
# endif

     !------------------------------------------------------------------------
     !   Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
     !-------------------------------------------------------------------------

       IF (lflux_gas_o2 .OR. lflux_gas_co2) THEN
         XKW_ICE(istr:iend,jstr:jend) = XKW(istr:iend,jstr:jend)
         WHERE (FICE_USED(istr:iend,jstr:jend) > 0.2 .AND. FICE_USED(istr:iend,jstr:jend) < 0.9999)
            XKW_ICE(istr:iend,jstr:jend) = (c1 - FICE_USED(istr:iend,jstr:jend)) * XKW_ICE(istr:iend,jstr:jend)
         END WHERE
         WHERE (FICE_USED(istr:iend,jstr:jend) >= 0.9999)
           XKW_ICE(istr:iend,jstr:jend) = c0
         END WHERE
       END IF

#if defined Ncycle_SY
     !------------------------------------------------------------------------
     !   compute N2O flux
     !------------------------------------------------------------------------
!        lflux_gas_n2o = .true.
       IF (lflux_gas_n2o) THEN
         call CSCHMIDT_N2O(SSTT(istr:iend,jstr:jend),
     &     landmask(istr:iend,jstr:jend),SCHMIDT_USED(istr:iend,jstr:jend),
     &     istr,iend,jstr,jend)


#if defined BEC2_DIAG && defined Ncycle_SY
       if (wrt_bec2_diag_2d(schmidt_n2o_idx_t)) then
          idiag = idx_bec2_diag_2d(schmidt_n2o_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = SCHMIDT_USED(istr:iend,jstr:jend)
       endif
#endif

         call N2OSATU(SSTT(istr:iend,jstr:jend),
     &     SSSS(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &     N2OSAT_1atm(istr:iend,jstr:jend),istr, iend,jstr,jend)
         do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
        XN2O = 300.0 * 1e-9 !!! need to figure out namelists here
               N2OSAT_USED = AP_USED(i,j) * N2OSAT_1atm(i,j) * XN2O(i)
# ifdef MASKING
     &                      * rmask(i,j)
# endif

               !!! calculating gas exchange for all n2o tracer
               !!! we let n2o_atm trace the invasion of atmospheric n2o
               !!! and set the atmopspheric concentration of all other
               !!! tracer to 0 So that everything is additive.
               FLUX(i,j) = PV(i,j) * (N2OSAT_USED - t(i,j,N,nnew,iN2O))
               STF(i,j,bgc_idx(iN2O)) = FLUX(i,j)
               !!! Update tracer
               t(i,j,N,nnew,iN2O) = t(i,j,N,nnew,iN2O)+
     &                 stf(i,j,bgc_idx(iN2O))*dt/Hz(i,j,N)

# if defined BEC2_DIAG && defined Ncycle_SY
       if (wrt_bec2_diag_2d(fgn2o_idx_t)) then
          idiag = idx_bec2_diag_2d(fgn2o_idx_t)
          bec2_diag_2d(i,j,idiag) = STF(i,j,bgc_idx(iN2O))
       endif
       if (wrt_bec2_diag_2d(pvn2o_idx_t)) then
           idiag = idx_bec2_diag_2d(pvn2o_idx_t)
           bec2_diag_2d(i,j,idiag) = PV(i,j)
       endif
       if (wrt_bec2_diag_2d(n2osat_idx_t)) then
          idiag = idx_bec2_diag_2d(n2osat_idx_t)
          bec2_diag_2d(i,j,idiag) = N2OSAT_USED
       endif
#endif

            end do  ! i
         end do  ! j
       END IF
# endif

# ifdef Ncycle_SY
     !------------------------------------------------------------------------
     !   compute N2 flux
     !------------------------------------------------------------------------
!        lflux_gas_n2 = .true.
       IF (lflux_gas_n2) THEN
         call CSCHMIDT_N2(SSTT(istr:iend,jstr:jend),
     &     landmask(istr:iend,jstr:jend),SCHMIDT_USED(istr:iend,jstr:jend),
     &     istr,iend,jstr,jend)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(schmidt_n2_idx_t)) then
          idiag = idx_bec2_diag_2d(schmidt_n2_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = SCHMIDT_USED(istr:iend,jstr:jend)
       endif
# endif
         call N2SATU(SSTT(istr:iend,jstr:jend),
     &     SSSS(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &     N2SAT_1atm(istr:iend,jstr:jend),istr, iend,jstr,jend)
         do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
               N2SAT_USED = AP_USED(i,j) * N2SAT_1atm(i,j) * 0.78084 !Molar ratio of N2 in air
# ifdef MASKING
     &                      * rmask(i,j)
# endif
!!!!!---------------Tracking only N2 excess for no2---------!!!!! to track full use:
               FLUX(i,j) = PV(i,j) * (c0 - t(i,j,N,nnew,iN2))
               STF(i,j,bgc_idx(iN2)) = FLUX(i,j)
               t(i,j,N,nnew,iN2) = t(i,j,N,nnew,iN2)+
     &                 stf(i,j,bgc_idx(iN2))*dt/Hz(i,j,N)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(pvn2_idx_t)) then
          idiag = idx_bec2_diag_2d(pvn2_idx_t)
          bec2_diag_2d(i,j,idiag) = PV(i,j)
       endif
       if (wrt_bec2_diag_2d(n2sat_idx_t)) then
          idiag = idx_bec2_diag_2d(n2sat_idx_t)
          bec2_diag_2d(i,j,idiag) = N2SAT_USED
       endif
       if (wrt_bec2_diag_2d(fgn2_idx_t)) then
          idiag = idx_bec2_diag_2d(fgn2_idx_t)
          bec2_diag_2d(i,j,idiag) = FLUX(i,j)
       endif
# endif
            end do  ! i
         end do  ! j
       END IF

# endif
    !------------------------------------------------------------------------
    !   compute O2 flux
    !------------------------------------------------------------------------

       IF (lflux_gas_o2) THEN
         call CSCHMIDT_O2(SSTT(istr:iend,jstr:jend),
     &     landmask(istr:iend,jstr:jend),SCHMIDT_USED(istr:iend,jstr:jend),
     &     istr,iend,jstr,jend)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(schmidto2_idx_t)) then
          idiag = idx_bec2_diag_2d(schmidto2_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = SCHMIDT_USED(istr:iend,jstr:jend)
       endif
# endif
         call O2SATU(SSTT(istr:iend,jstr:jend),
     &     SSSS(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &     O2SAT_1atm(istr:iend,jstr:jend),istr, iend,jstr,jend)
         do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
               O2SAT_USED = AP_USED(i,j) * O2SAT_1atm(i,j)
# ifdef MASKING
     &                      * rmask(i,j)
# endif
               FLUX(i,j) = PV(i,j) * (O2SAT_USED - t(i,j,N,nnew,iO2))
               STF(i,j,bgc_idx(iO2)) = FLUX(i,j)
               t(i,j,N,nnew,iO2) = t(i,j,N,nnew,iO2)+
     &                 stf(i,j,bgc_idx(iO2))*dt/Hz(i,j,N)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(pvo2_idx_t)) then
          idiag = idx_bec2_diag_2d(pvo2_idx_t)
          bec2_diag_2d(i,j,idiag) = PV(i,j)
       endif
       if (wrt_bec2_diag_2d(o2sat_idx_t)) then
          idiag = idx_bec2_diag_2d(o2sat_idx_t)
          bec2_diag_2d(i,j,idiag) = O2SAT_USED
       endif
       if (wrt_bec2_diag_2d(fgo2_idx_t)) then
          idiag = idx_bec2_diag_2d(fgo2_idx_t)
          bec2_diag_2d(i,j,idiag) = FLUX(i,j)
       endif
# endif
            end do  ! i
         end do  ! j
       END IF

    !------------------------------------------------------------------------
    !   compute CO2 flux, computing disequilibrium one row at a time
    !------------------------------------------------------------------------
       IF (lflux_gas_co2) THEN
          call CSCHMIDT_CO2(SSTT(istr:iend,jstr:jend),
     &       landmask(istr:iend,jstr:jend),
     &       SCHMIDT_USED(istr:iend,jstr:jend),
     &       istr,iend,jstr,jend)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(schmidtco2_idx_t)) then
          idiag = idx_bec2_diag_2d(schmidtco2_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = SCHMIDT_USED(istr:iend,jstr:jend)
       endif
# endif

          do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(pvco2_idx_t)) then
          idiag = idx_bec2_diag_2d(pvco2_idx_t)
          bec2_diag_2d(i,j,idiag) = PV(i,j)
       endif
# endif
            end do
          end do

          DO j = jstr,jend

# ifdef PCO2AIR_FORCING
             XCO2(istr:iend) = pco2air(istr:iend,j)
# endif /* PCO2AIR_FORCING */

            WHERE (ph_srf(istr:iend,j) .ne. c0)
                PHLO(istr:iend) = ph_srf(istr:iend,j) - del_ph
                PHHI(istr:iend) = ph_srf(istr:iend,j) + del_ph
            ELSEWHERE
                PHLO = phlo_init
                PHHI = phhi_init
            END WHERE

             CALL co2calc_row(LANDMASK(istr:iend,j),
     &         SSTT(istr:iend,j), SSSS(istr:iend,j),
     &         t(istr:iend,j,N,nnew,iDIC),
     &         t(istr:iend,j,N,nnew,iALK),
     &         t(istr:iend,j,N,nnew,iPO4),
     &         t(istr:iend,j,N,nnew,iSIO3),
     &         PHLO(istr:iend), PHHI(istr:iend), PH_NEW(istr:iend), XCO2(istr:iend),
     &         AP_USED(istr:iend,j), CO2STAR_ROW(istr:iend),
     &         DCO2STAR_ROW(istr:iend), pco2oc_ROW(istr:iend),
     &         DpCO2_ROW(istr:iend),istr,iend)

             ph_srf(istr:iend,j) = PH_NEW(istr:iend)
#  ifdef BEC2_DIAG
        if (wrt_bec2_diag_2d(ph_idx_t)) then
           idiag = idx_bec2_diag_2d(ph_idx_t)
           bec2_diag_2d(istr:iend,j,idiag) = PH_NEW(istr:iend)
        endif
        if (wrt_bec2_diag_2d(co2star_idx_t)) then
           idiag = idx_bec2_diag_2d(co2star_idx_t)
           bec2_diag_2d(istr:iend,j,idiag)  = CO2STAR_ROW(istr:iend)
        endif
        if (wrt_bec2_diag_2d(dco2star_idx_t)) then
           idiag = idx_bec2_diag_2d(dco2star_idx_t)
           bec2_diag_2d(istr:iend,j,idiag) = DCO2STAR_ROW(istr:iend)
        endif
        if (wrt_bec2_diag_2d(pco2oc_idx_t)) then
           idiag = idx_bec2_diag_2d(pco2oc_idx_t)
           bec2_diag_2d(istr:iend,j,idiag)   = pco2oc_ROW(istr:iend)
        endif
        if (wrt_bec2_diag_2d(pco2air_idx_t)) then
           idiag = idx_bec2_diag_2d(pco2air_idx_t)
           bec2_diag_2d(istr:iend,j,idiag)  = XCO2(istr:iend)
        endif
#  endif /* BEC2_DIAG */

             FLUX(istr:iend,j) = PV(istr:iend,j) * DCO2STAR_ROW(istr:iend)
          END DO  ! j = jstr,jend

          STF(istr:iend,jstr:jend,bgc_idx(iDIC)) = FLUX(istr:iend,jstr:jend)

          t(istr:iend,jstr:jend,N,nnew,iDIC) =
     &       t(istr:iend,jstr:jend,N,nnew,iDIC) +
     &       stf(istr:iend,jstr:jend,bgc_idx(iDIC))
     &       *dt/Hz(istr:iend,jstr:jend,N)
# ifdef BEC2_DIAG
        if (wrt_bec2_diag_2d(fgco2_idx_t)) then
           idiag = idx_bec2_diag_2d(fgco2_idx_t)
           bec2_diag_2d(istr:iend,jstr:jend,idiag) = FLUX(istr:iend,jstr:jend)
        endif
# endif
       ENDIF  ! lflux_gas_co2

    !---------------------------------------------------------------------
    !   Surface fluxes:
    !----------------------------------------------------------------------

       !     iron fluxes in forcing file are in nmol Fe cm-2 s-1
       !     BEC needs them in mmol Fe m-2 s-1
       if (liron_flux) then
          FLUX(istr:iend,jstr:jend)=iron(istr:iend,jstr:jend)
       else
          FLUX = c0
       endif
       FLUX = FLUX * parm_Fe_bioavail * nmol_cm2_to_mmol_m2
       STF(istr:iend,jstr:jend,bgc_idx(iFE)) = FLUX(istr:iend,jstr:jend)
       t(istr:iend,jstr:jend,N,nnew,iFE) =
     &         t(istr:iend,jstr:jend,N,nnew,iFE) +
     &         stf(istr:iend,jstr:jend,bgc_idx(iFE))
     &         *dt/Hz(istr:iend,jstr:jend,N)

# ifdef NOX_FORCING
       ! NOx: added to nitrate (convert units from kg N to mmol N)
       t(istr:iend,jstr:jend,N,nnew,iNO3) =
     &         t(istr:iend,jstr:jend,N,nnew,iNO3) +
     &         nox(istr:iend,jstr:jend)*71394.200220751
     &         *dt/Hz(istr:iend,jstr:jend,N)
# endif
# ifdef NHY_FORCING
       ! NHy: added to ammonium (convert units from kg N to mmol N)
       t(istr:iend,jstr:jend,N,nnew, iNH4) =
     &         t(istr:iend,jstr:jend,N,nnew,iNH4) +
     &         nhy(istr:iend,jstr:jend)*71394.200220751
     &         *dt/Hz(istr:iend,jstr:jend,N)
# endif

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(ironflux_idx_t)) then
          idiag = idx_bec2_diag_2d(ironflux_idx_t)
          bec2_diag_2d(istr:iend,jstr:jend,idiag) = FLUX(istr:iend,jstr:jend)
       endif
# endif

     !---------------------------------------------------------------------
     ! Now we compute the tendencies for the biological tracers
     ! (stored in array DTRACER_MODULE).
     !
     ! NOTE: the k-loop needs to be in this order
     ! because PAR etc. is initialized for k==N
     !---------------------------------------------------------------------

       do j = jstr, jend

          do k = N,1,-1

          call ecosys_set_interior(k,t(istr:iend,j,k,nnew,1),
     &             SHF_QSW(istr:iend,j),
     &             PAR_out(istr:iend,j),
#if defined DAILYPAR_PHOTOINHIBITION
     &             SHF_QSW_AVG(istr:iend,j),
     &             PAR_out_AVG(istr:iend,j),
#endif
     &             qa_dust_def(istr:iend,j),istr,iend,j,
     &             dt,dust(istr:iend,j))
          end do
       end do

     !---------------------------------------------------------------------
     ! Compute and write back new tracer values to ROMS:
     ! the explicit Euler method is used for computing the new
     ! biological tracer values.
     !---------------------------------------------------------------------

       do k=1,n
          do m=1,ntrc_bio
             t(istr:iend,jstr:jend,k,nnew,isalt+nt_passive+m) =
     &             t(istr:iend,jstr:jend,k,nnew,isalt+nt_passive+m)
     &             + dt * DTRACER_MODULE(istr:iend,jstr:jend,k,m)
          enddo

       enddo

!       call wrt_his_bgc_tracers
!       call wrt_avg_bgc_tracers
!       call wrt_bgc

       end subroutine ecosys_bec2_tile

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE ECOSYS_SET_INTERIOR   ********************
!
!-------------------------------------------------------------------------

       subroutine ecosys_set_interior(k,temp,SHF_QSW,
     &  PAR_out,
#if defined DAILYPAR_PHOTOINHIBITION
     &  SHF_QSW_AVG,PAR_out_AVG,
#endif
     &  qa_dust_def,istr,iend,j,dt2,dust_flux)

       use param
       use grid
       use ocean_vars
       use scalars
       use bgc_io
       use bgc_shared_vars
       use bec2_params

        implicit none

       integer istr,iend,i,j,k,n_loc,auto_ind,auto_ind2
       integer idiag
       real dt2,tmp1,tmp2,tmp3,tmp4
       real temp(istr:iend), SHF_QSW(istr:iend),
     &    PAR_out(istr:iend), dust_flux(istr:iend)

#if defined DAILYPAR_PHOTOINHIBITION
      real SHF_QSW_AVG(istr:iend), PAR_out_AVG(istr:iend)
#endif

       REAL, PARAMETER ::
     &    epsC      = 1.00e-8,  ! small C concentration (mmol C/m^3)
     &    epsTinv   = 3.17e-8,  ! small inverse time scale (1/year) (1/sec)
     &    epsnondim = 1.00e-6,  ! small non-dimensional number (non-dim)
     &    cks = 9.0,         ! constant used in Fe quota modification
     &    cksi = 5.0         ! constant used in Si quota modification

       REAL, DIMENSION(istr:iend) ::
     &    PO4_loc,         ! local copy of model PO4
     &    NO3_loc,         ! local copy of model NO3
     &    SiO3_loc,        ! local copy of model SiO3
     &    NH4_loc,         ! local copy of model NH4
     &    Fe_loc,          ! local copy of model Fe
     &    O2_loc,          ! local copy of model O2
     &    DOC_loc,         ! local copy of model DOC
     &    spC_loc,         ! local copy of model spC
     &    spChl_loc,       ! local copy of model spChl
     &    spCaCO3_loc,     ! local copy of model spCaCO3
     &    diatC_loc,       ! local copy of model diatC
     &    diatChl_loc,     ! local copy of model diatChl
     &    zooC_loc,        ! local copy of model zooC
     &    spFe_loc,        ! local copy of model spFe
     &    diatSi_loc,      ! local copy of model diatSi
     &    diatFe_loc,      ! local copy of model diatFe
     &    diazC_loc,       ! local copy of model diazC
     &    diazChl_loc,     ! local copy of model diazChl
     &    diazFe_loc,      ! local copy of model diazFe
     &    DON_loc,         ! local copy of model DON
     &    DOFe_loc,        ! local copy of model DOFe
     &    DOP_loc,         ! local copy of model DOP
     &    DOPr_loc,        ! local copy of model DOPr
     &    DONr_loc,        ! local copy of model DONr
# ifdef Ncycle_SY
     &    NO2_loc,         ! local copy of model no2
     &    N2_loc,         ! local copy of model n2
     &    N2O_loc,         ! N2O_ao1_loc + N2O_soden_loc - N2O_siden_loc
     &    ammox,          ! ammonium oxidation (NH4 -> NO2) (mmol N/m^3/sec)
     &    nitrox,          ! nitrite oxidation (NO2 -> NO3) (mmolN/m^3/sec)
     &    anammox,          ! anaerobic ammonium oxidation (NH4 -> N2) (mmolN/m^3/sec)
     &    epsN,         ! small value for denitrification calculations
     &    R_oxic,        ! Pseudo oxic remineralization (NO3-> NO2, no units)
     &    R_den1,        ! Pseudo denitrification 1 (NO3-> NO2, no units)
     &    R_den2,        ! pseudo denitrification 2 (NO2-> N2O), no units)
     &    R_den3,        ! Pseudo denitrification 3 (N2O-> N2, no units)
     &    roxic,       ! portion of C remin done by denitrif that is done by denitrif
     &    rden1,       ! portion of C remin done by denitrif that is done by denitrif
     &    rden2,       ! portion of C remin done by denitrif that is done by denitrif
     &    rden3,       ! portion of C remin done by denitrif that is done by denitrif
     &    Rden,       ! portion of C remin done by denitrif
     &    DENITRIF1,       ! NO3 consumed during denitrification
     &    DENITRIF2,       ! NO2 consumed during denitrification
     &    DENITRIF3,       ! N2O consumed during denitrification
     &    n2o_nh4_yield1,   ! N2O yield during ammonium oxidation (0 --> 1)
     &    no2_nh4_yield1,   ! NO2 yield during ammonium oxidation (0 --> 1)
# endif
     &    QA_dust_def,     ! incoming deficit in the QA(dust) POC flux
     &    NITRIF,          ! nitrification (NH4 -> NO3) (mmol N/m^3/sec)
     &    DENITRIF,        ! water column nitrification (NO3 -> N2) (mmol N/m^3/sec)
     &    SED_DENITRIF,    ! sedimentary denitrification (mmol N/m^3/sec)
     &    OTHER_REMIN,     ! organic C remin not due oxic or denitrif (mmolC/m^3/sec)
     &    flux_oxidated,   ! carbon oxidation rate in sediment
     &    O2_CONSUMPTION,  ! O2 consumption
     &    O2_PRODUCTION    ! O2 production

       REAL, DIMENSION(istr:iend) ::
     &    Tfunc,           ! temp response function GD98 (non-dim)
     &    Tfunc_diat,      ! temp response function GD98 (non-dim) for diatoms
     &    PAR_in,          ! photosynthetically available radiation (W/m^2)
#if defined DAILYPAR_PHOTOINHIBITION
     &    PAR_in_AVG,      ! daily average photosynthetically available radiation (W/m^2)
# endif
     &    KPARdz,          ! PAR adsorption coefficient (non-dim)
     &    PAR_lay,         ! average PAR over layer depth (W/m^2)
     &    f_loss_thres,    ! fraction of grazing loss reduction at depth
     &    C_loss_thres,    ! bio-C threshold at which losses go to zero (mmol C/m^3)
     &    VFe,             ! C-specific Fe uptake (non-dim)
     &    VSiO3,           ! C-specific SiO3 uptake (non-dim)
     &    f_nut,           ! nut limitation factor, modifies C fixation (non-dim)
     &    PCmax,           ! max value of PCphoto at temperature TEMP (1/sec)
     &    light_lim,       ! light limitation factor
     &    PCphoto,         ! C-specific rate of photosynth. (1/sec)
     &    pChl,            ! Chl synth. regulation term (mg Chl/mmol N)
     &    VNC,             ! C-specific N uptake rate (mmol N/mmol C/sec)
     &    z_umax,          ! max. zoo growth rate at local T (1/sec)
     &    remaining_P,     ! used in routing P from autotrophs w/ Qp different from Qp_zoo_pom
     &    Zprime,          ! used to limit zoo mort at low biomass (mmol C/m^3)
     &    zoo_loss,        ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
     &    zoo_loss_doc,    ! zoo_loss routed to doc (mmol C/m^3/sec)
     &    zoo_loss_dic,    ! zoo_loss routed to dic (mmol C/m^3/sec)
     &    f_zoo_detr_loc,  ! frac of zoo losses into large detrital pool (non-dim)
     &    DOC_prod,        ! production of DOC (mmol C/m^3/sec)
     &    DON_prod,        ! production of dissolved organic N
     &    DOFe_prod,       ! production of dissolved organic Fe
     &    DOP_prod,        ! production of dissolved organic P
     &    DOC_remin,       ! remineralization of DOC (mmol C/m^3/sec)
     &    DON_remin,       ! portion of DON remineralized
     &    DOFe_remin,      ! portion of DOFe remineralized
     &    DONr_remin,      ! portion of refractory DON remineralized
     &    POC_prod,        ! POC production (mmol C/m^3/sec)
     &    P_CaCO3_prod,    ! large detritus CaCO3, base units = mmol CaCO3
     &    P_SiO2_prod,     ! prod. of large detritus SiO2
     &    Fe_scavenge,     ! loss of dissolved iron, scavenging (mmol Fe/m^3/sec)
     &    P_iron_prod,     ! P iron production (base units = mmol Fe)
     &    Fe_scavenge_rate, ! annual scavenging rate of iron as % of ambient
     &    WORK1,WORK2      ! temporary arrays

       REAL, DIMENSION(istr:iend,autotroph_cnt) ::
     &    gQfe,            ! fe/C for growth
     &    Qfe,             ! init fe/C ratio (mmolFe/mmolC)
     &    Qsi,             ! Diatom initial Si/C ratio (mmol Si/mmol C)
     &    gQsi,            ! diatom Si/C ratio for growth (new biomass)
     &    thetaC,          ! local Chl/C ratio (mg Chl/mmol C)
     &    Pprime,          ! used to limit autotroph mort at low biomass (mmol C/m^3)
     &    C_loc,           ! amount of carbon in each autotroph species
     &    Chl_loc,         ! amount of chlorophyll in each autotroph species
     &    autoFe_loc,      ! amount of iron in each autotroph species
     &    Si_loc,          ! amount of Si in each autotroph species
     &    CaCO3_loc,       ! amount of CaCO3 in each autotroph species
     &    QCaCO3,          ! CaCO3/C ratio (mmol CaCO3/mmol C)
# ifdef Ncycle_SY
     &    VNO3,            ! NO3 uptake rate (non-dim)
     &    VNO2,            ! NO2 uptake rate (non-dim)
# else
     &    VNO3,            ! NO3 uptake rate (non-dim)
# endif
     &    VNH4,            ! NH4 uptake rate (non-dim)
     &    VNtot,           ! total N uptake rate (non-dim)
     &    VPO4,            ! C-specific PO4 uptake (non-dim)
     &    VDOP,            ! C-specific DOP uptake rate (non-dim)
     &    VPtot,           ! total P uptake rate/limitation (non-dim)
     &    photoC,          ! C-fixation (mmol C/m^3/sec)
     &    photoFe,         ! iron uptake
     &    photoSi,         ! silicon uptake (mmol Si/m^3/sec)
     &    NO3_V,           ! nitrate uptake (mmol NO3/m^3/sec)
#ifdef Ncycle_SY
     &    NO2_V,           ! nitrate uptake (mmol NO3/m^3/sec)
# endif
     &    NH4_V,           ! ammonium uptake (mmol NH4/m^3/sec)
     &    PO4_V,           ! PO4 uptake (mmol PO4/m^3/sec)
     &    DOP_V,           ! DOP uptake (mmol DOP/m^3/sec)
     &    CaCO3_PROD,      ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
     &    photoacc,        ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
     &    auto_loss,       ! autotroph non-grazing mort (mmol C/m^3/sec)
     &    auto_loss_poc,   ! auto_loss routed to poc (mmol C/m^3/sec)
     &    auto_loss_doc,   ! auto_loss routed to doc (mmol C/m^3/sec)
     &    auto_loss_dic,   ! auto_loss routed to dic (mmol C/m^3/sec)
     &    auto_agg,        ! autotroph aggregation (mmol C/m^3/sec)
     &    auto_graze,      ! autotroph grazing rate (mmol C/m^3/sec)
     &    auto_graze_zoo,  ! auto_graze routed to zoo (mmol C/m^3/sec)
     &    auto_graze_poc,  ! auto_graze routed to poc (mmol C/m^3/sec)
     &    auto_graze_doc,  ! auto_graze routed to doc (mmol C/m^3/sec)
     &    auto_graze_dic,  ! auto_graze routed to dic (mmol C/m^3/sec)
     &    remaining_P_dop, ! remaining_P from mort routed to DOP pool
     &    remaining_P_dip, ! remaining_P from mort routed to remin
     &    Nfix,            ! total Nitrogen fixation (mmol N/m^3/sec)
     &    Nexcrete         ! fixed N excretion


    !---------------------------------------------------------------------------
    !   exit immediately if computations are not to be performed
    !---------------------------------------------------------------------------

       IF (.NOT. lsource_sink) RETURN

    !---------------------------------------------------------------------------
    !   create local copies of model tracers
    !   treat negative values as zero
    !---------------------------------------------------------------------------
       PO4_loc  = MAX(c0, T(istr:iend,j,k,nnew,iPO4))
       NO3_loc  = MAX(c0, T(istr:iend,j,k,nnew,iNO3))
       SiO3_loc = MAX(c0, T(istr:iend,j,k,nnew,iSIO3))
       NH4_loc  = MAX(c0, T(istr:iend,j,k,nnew,iNH4))
       Fe_loc   = MAX(c0, T(istr:iend,j,k,nnew,iFE))
       O2_loc   = MAX(c0, T(istr:iend,j,k,nnew,iO2))
       DOC_loc  = MAX(c0, T(istr:iend,j,k,nnew,iDOC))
       spC_loc  = MAX(c0, T(istr:iend,j,k,nnew,iSPC))
       spChl_loc= MAX(c0, T(istr:iend,j,k,nnew,iSPCHL))
       spCaCO3_loc  = MAX(c0,T(istr:iend,j,k,nnew,iSPCACO3))
       diatC_loc= MAX(c0, T(istr:iend,j,k,nnew,iDIATC))
       diatChl_loc  = MAX(c0, T(istr:iend,j,k,nnew,iDIATCHL))
       zooC_loc  = MAX(c0, T(istr:iend,j,k,nnew,iZOOC))
       spFe_loc  = MAX(c0, T(istr:iend,j,k,nnew,iSPFE))
       diatSi_loc = MAX(c0, T(istr:iend,j,k,nnew,iDIATSI))
       diatFe_loc  = MAX(c0, T(istr:iend,j,k,nnew,iDIATFE))
       diazC_loc    = MAX(c0, T(istr:iend,j,k,nnew,iDIAZC))
       diazChl_loc  = MAX(c0, T(istr:iend,j,k,nnew,iDIAZCHL))
       diazFe_loc   = MAX(c0, T(istr:iend,j,k,nnew,iDIAZFE))
       DON_loc  = MAX(c0, T(istr:iend,j,k,nnew,iDON))
       DOFe_loc = MAX(c0, T(istr:iend,j,k,nnew,iDOFE))
       DOP_loc  = MAX(c0, T(istr:iend,j,k,nnew,iDOP))
       DOPr_loc = MAX(c0, T(istr:iend,j,k,nnew,iDOPR))
       DONr_loc = MAX(c0, T(istr:iend,j,k,nnew,iDONR))
# ifdef Ncycle_SY
       NO2_loc  = MAX(c0, T(istr:iend,j,k,nnew,iNO2))
       N2_loc  = MAX(c0, T(istr:iend,j,k,nnew,iN2))
       N2O_loc = MAX(c0, t(istr:iend,j,k,nnew,iN2O))
# endif

#ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(totchl_idx_t)) then
          idiag = idx_bec2_diag_3d(totchl_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = t(istr:iend,j,k,nnew,iSPCHL)
     &     + t(istr:iend,j,k,nnew,iDIATCHL) + t(istr:iend,j,k,nnew,iDIAZCHL)
       endif
       if (wrt_bec2_diag_3d(totphytoc_idx_t)) then
          idiag = idx_bec2_diag_3d(totphytoc_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = t(istr:iend,j,k,nnew,iSPC)
     &     + t(istr:iend,j,k,nnew,iDIATC) + t(istr:iend,j,k,nnew,iDIAZC)
       endif
       if (wrt_bec2_diag_3d(spczero_idx_t)) then
          idiag = idx_bec2_diag_3d(spczero_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = max(c0,-T(istr:iend,j,k,nnew,iSPC))/dt
       endif
       if (wrt_bec2_diag_3d(diatczero_idx_t)) then
          idiag = idx_bec2_diag_3d(diatczero_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = max(c0,-T(istr:iend,j,k,nnew,iDIATC))/dt
       endif
       if (wrt_bec2_diag_3d(diazczero_idx_t)) then
          idiag = idx_bec2_diag_3d(diazczero_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = max(c0,-T(istr:iend,j,k,nnew,iDIAZC))/dt
       endif
       if (wrt_bec2_diag_3d(doczero_idx_t)) then
          idiag = idx_bec2_diag_3d(doczero_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = max(c0,-T(istr:iend,j,k,nnew,iDOC))/dt
       endif
       if (wrt_bec2_diag_3d(zooczero_idx_t)) then
          idiag = idx_bec2_diag_3d(zooczero_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = max(c0,-T(istr:iend,j,k,nnew,iZOOC))/dt
       endif
       if (wrt_bec2_diag_3d(spcaco3zero_idx_t)) then
          idiag = idx_bec2_diag_3d(spcaco3zero_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = max(c0,-T(istr:iend,j,k,nnew,iSPCACO3))/dt
       endif
#endif

    !---------------------------------------------------------------------------
    !   If any phyto box are zero, set others to zeros.
    !---------------------------------------------------------------------------

       WHERE (spC_loc == c0 .OR. spChl_loc == c0 .OR. spFe_loc == c0)
          spC_loc = c0
          spChl_loc = c0
          spCaCO3_loc = c0
          spFe_loc = c0
       END WHERE

       WHERE (diatC_loc == c0 .OR. diatChl_loc == c0 .OR.
     &     diatFe_loc == c0
     &    .OR. diatSi_loc == c0)
          diatC_loc = c0
          diatChl_loc = c0
          diatFe_loc = c0
          diatSi_loc = c0
       END WHERE

       WHERE (diazC_loc == c0 .OR. diazChl_loc == c0 .OR.
     &   diazFe_loc == c0)
          diazC_loc = c0
          diazChl_loc = c0
          diazFe_loc = c0
       END WHERE

     !---------------------------------------------------------------------------
     ! Initialize C_loc, Chl_loc, CaCO3_loc, autoFe_loc:
     !---------------------------------------------------------------------------

       C_loc(istr:iend,sp_ind) = spC_loc
       C_loc(istr:iend,diat_ind) = diatC_loc
       C_loc(istr:iend,diaz_ind) = diazC_loc
       Chl_loc(istr:iend,sp_ind) = spChl_loc
       Chl_loc(istr:iend,diat_ind) = diatChl_loc
       Chl_loc(istr:iend,diaz_ind) = diazChl_loc
       CaCO3_loc(istr:iend,sp_ind) = spCaCO3_loc
       CaCO3_loc(istr:iend,diat_ind) = c0
       CaCO3_loc(istr:iend,diaz_ind) = c0
       Si_loc(istr:iend,sp_ind) = c0
       Si_loc(istr:iend,diat_ind) = diatSi_loc
       Si_loc(istr:iend,diaz_ind) = c0
       autoFe_loc(istr:iend,sp_ind) = spFe_loc
       autoFe_loc(istr:iend,diat_ind) = diatFe_loc
       autoFe_loc(istr:iend,diaz_ind) = diazFe_loc

     !---------------------------------------------------------------------------
     !   set local variables, with incoming ratios
     !---------------------------------------------------------------------------

       do auto_ind = 1, autotroph_cnt
          thetaC(istr:iend,auto_ind) = Chl_loc(istr:iend,auto_ind) / (C_loc(istr:iend,auto_ind) + epsC)
          Qfe(istr:iend,auto_ind) = autoFe_loc(istr:iend,auto_ind) / (C_loc(istr:iend,auto_ind) + epsC)
          if (kSiO3(auto_ind) > c0) then
             Qsi(istr:iend,auto_ind) = min(Si_loc(istr:iend,auto_ind) / (C_loc(istr:iend,auto_ind) + epsC), gQsi_max)
          end if
       end do

     !---------------------------------------------------------------------------
     !   DETERMINE NEW ELEMENTAL RATIOS FOR GROWTH (NEW BIOMASS)
     ! Modify these initial ratios under low ambient iron conditions
     ! Modify the initial si/C ratio under low ambient Si conditions
     !---------------------------------------------------------------------------

       do auto_ind = 1, autotroph_cnt
          gQfe(istr:iend,auto_ind) = gQfe_0(auto_ind)
          where (Fe_loc < cks * kFe(auto_ind))
             gQfe(istr:iend,auto_ind) =
     &           max((gQfe(istr:iend,auto_ind) * Fe_loc / (cks * kFe(auto_ind))),
     &           gQfe_min(auto_ind))
          end where

          if (kSiO3(auto_ind) > c0) then
             gQsi(istr:iend,auto_ind) = gQsi_0
             where ((Fe_loc < cksi * kFe(auto_ind)) .and. (Fe_loc > c0) .and.
     &           (SiO3_loc > (cksi * kSiO3(auto_ind))))
                gQsi(istr:iend,auto_ind) = min((gQsi(istr:iend,auto_ind) * cksi * kFe(auto_ind) / Fe_loc), gQsi_max)
             end where

             where (Fe_loc == c0)
                gQsi(istr:iend,auto_ind) = gQsi_max
             end where

             where (SiO3_loc < (cksi * kSiO3(auto_ind)))
                gQsi(istr:iend,auto_ind) = max((gQsi(istr:iend,auto_ind) * SiO3_loc / (cksi * kSiO3(auto_ind))),
     &                                  gQsi_min)
             end where
          endif

          if (CaCO3_ind(auto_ind) > 0) then
             QCaCO3(istr:iend,auto_ind) = min(CaCO3_loc(istr:iend,auto_ind)
     &                      / (C_loc(istr:iend,auto_ind) + epsC), QCaCO3_max)
          end if
       end do  ! auto_ind = 1, autotroph_cnt
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(spqcaco3_idx_t)) then
          idiag = idx_bec2_diag_3d(spqcaco3_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = QCaCO3(istr:iend,sp_ind)
       endif
# endif

     !---------------------------------------------------------------------------
     !   various k==N initializations
     !
     !   f_qsw_par   fraction of incoming SW -> PAR (non-dim)
     !---------------------------------------------------------------------------

       if (k == N) then
#if defined DAILYPAR_PHOTOINHIBITION
          PAR_out_AVG = MAX(c0, f_qsw_par * SHF_QSW_AVG)
# ifdef MASKING
     &            * rmask(istr:iend,j)
# endif
# endif
          PAR_out = MAX(c0, f_qsw_par * SHF_QSW)
# ifdef MASKING
     &            * rmask(istr:iend,j)
# endif

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(parinc_idx_t)) then
          idiag = idx_bec2_diag_2d(parinc_idx_t)
          bec2_diag_2d(istr:iend,j,idiag) = PAR_out(istr:iend)
       endif
# endif
          CALL init_particulate_terms(QA_dust_def,istr,iend,j,dust_flux)
       end if
     !-----------------------------------------------------------------------
     !  compute PAR related quantities
     !  Morel, Maritorena, JGR, Vol 106, No. C4, pp 7163--7180, 2001
     !-----------------------------------------------------------------------

       do i = istr, iend
#if defined DAILYPAR_PHOTOINHIBITION
          PAR_in_AVG(i) = PAR_out_AVG(i)
#endif
          PAR_in(i) = PAR_out(i)
          KPARdz(i) = max(spChl_loc(i)+diatChl_loc(i)+diazChl_loc(i), 0.02)
      ! in the following computation of KPARdz, it was taken into account that Hz is in m
          if (KPARdz(i) < 0.13224) then
             KPARdz(i) = 0.0919*(KPARdz(i)**0.3536) * Hz(i,j,k)
          else
             KPARdz(i) = 0.1131*(KPARdz(i)**0.4562) * Hz(i,j,k)
          end if
#if defined DAILYPAR_PHOTOINHIBITION
          PAR_out_AVG(i) = PAR_in_AVG(i) * exp(-KPARdz(i))
#endif
          PAR_out(i) = PAR_in(i) * exp(-KPARdz(i))
          PAR_lay(i) = PAR_in(i) * (c1 - exp(-KPARdz(i))) / KPARdz(i)
# ifdef BEC2_DIAG
        if (wrt_bec2_diag_3d(par_idx_t)) then
          idiag = idx_bec2_diag_3d(par_idx_t)
          bec2_diag_3d(i,j,k,idiag) = PAR_lay(i)
        endif
# endif
       end do

     !-----------------------------------------------------------------------
     !  Tref = 30.0 reference temperature (deg. C)
     !
     !  Using q10 formulation with Q10 value of 2.0 (Doney et al., 1996).
     !  growth, mort and grazing rates scaled by Tfunc where they are computed
     !-----------------------------------------------------------------------

       Tfunc = Q_10**(((temp + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)
       Tfunc_diat  = 1.55**(((temp + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)

     !-----------------------------------------------------------------------
     !  calculate the loss threshold interpolation factor
     !-----------------------------------------------------------------------

       do i = istr, iend
          if (-z_r(i,j,k) > thres_z1) then
             if (-z_r(i,j,k) < thres_z2) then
                f_loss_thres(i) = (thres_z2 + z_r(i,j,k))/(thres_z2 - thres_z1)
             else
                f_loss_thres(i) = c0
             endif
          else
             f_loss_thres(i) = c1
          endif
       end do

     !-----------------------------------------------------------------------
     !  Compute Pprime for all autotrophs, used for loss terms
     !-----------------------------------------------------------------------

       do auto_ind = 1, autotroph_cnt
          where (temp < temp_thres(auto_ind))
             C_loss_thres = f_loss_thres * loss_thres2(auto_ind)
          elsewhere
             C_loss_thres = f_loss_thres * loss_thres(auto_ind)
          end where
          Pprime(istr:iend,auto_ind) = max(C_loc(istr:iend,auto_ind) - C_loss_thres, c0)
       end do

     !-----------------------------------------------------------------------
     !  Get relative nutrient uptake rates for autotrophs,
     !  min. relative uptake rate modifies C fixation in the manner
     !  that the min. cell quota does in GD98.
     !-----------------------------------------------------------------------

       do auto_ind = 1, autotroph_cnt

# ifdef Ncycle_SY
          VNO3(istr:iend,auto_ind) = (NO3_loc / kNO3(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NO2_loc / kNO2(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNO2(istr:iend,auto_ind) = (NO2_loc / kNO2(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NO2_loc / kNO2(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNH4(istr:iend,auto_ind) = (NH4_loc / kNH4(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NO2_loc / kNO2(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNtot(istr:iend,auto_ind) = VNO3(istr:iend,auto_ind) + VNO2(istr:iend,auto_ind) + VNH4(istr:iend,auto_ind)
          if (Nfixer(auto_ind)) VNtot(istr:iend,auto_ind) = c1
# else
          VNO3(istr:iend,auto_ind) = (NO3_loc / kNO3(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNH4(istr:iend,auto_ind) = (NH4_loc / kNH4(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNtot(istr:iend,auto_ind) = VNO3(istr:iend,auto_ind) + VNH4(istr:iend,auto_ind)
          if (Nfixer(auto_ind)) VNtot(istr:iend,auto_ind) = c1
# endif

          VFe = Fe_loc / (Fe_loc + kFe(auto_ind))
          f_nut = min(VNtot(istr:iend,auto_ind), VFe)

          VPO4(istr:iend,auto_ind) = (PO4_loc / kPO4(auto_ind)) /
     &           (c1 + (PO4_loc / kPO4(auto_ind)) + (DOP_loc / kDOP(auto_ind)))
          VDOP(istr:iend,auto_ind) = (DOP_loc / kDOP(auto_ind)) /
     &           (c1 + (PO4_loc / kPO4(auto_ind)) + (DOP_loc / kDOP(auto_ind)))
          VPtot(istr:iend,auto_ind) = VPO4(istr:iend,auto_ind) + VDOP(istr:iend,auto_ind)

          !f_nut = min(f_nut, VPO4(istr:iend,auto_ind))

          ! consider VPtot instead of VPO4 only! DOP is taken up, so
          ! it should also potentially be limiting (that's why half-sat
          ! constants are defined for DOP!)
          f_nut = min(f_nut, VPtot(istr:iend,auto_ind))

# ifdef BEC2_DIAG
          if (auto_ind == sp_ind) then
             if (wrt_bec2_diag_3d(spfeuptake_idx_t)) then
                idiag = idx_bec2_diag_3d(spfeuptake_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = VFe
             endif
             if (wrt_bec2_diag_3d(sppo4uptake_idx_t)) then
                idiag = idx_bec2_diag_3d(sppo4uptake_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = VPO4(istr:iend,sp_ind)
             endif
          else if (auto_ind == diat_ind) then
             if (wrt_bec2_diag_3d(diatfeuptake_idx_t)) then
                idiag = idx_bec2_diag_3d(diatfeuptake_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = VFe
             endif
             if (wrt_bec2_diag_3d(diatpo4uptake_idx_t)) then
                idiag = idx_bec2_diag_3d(diatpo4uptake_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) =VPO4(istr:iend,diat_ind)
             endif
          else
             if (wrt_bec2_diag_3d(diazfeuptake_idx_t)) then
                idiag = idx_bec2_diag_3d(diazfeuptake_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = VFe
             endif
             if (wrt_bec2_diag_3d(diazpo4uptake_idx_t)) then
                idiag = idx_bec2_diag_3d(diazpo4uptake_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) =VPO4(istr:iend,diaz_ind)
             endif
          end if
# endif
          if (kSiO3(auto_ind) > c0) then
             VSiO3 = SiO3_loc / (SiO3_loc + kSiO3(auto_ind))
             f_nut = min(f_nut, VSiO3)
          endif
# ifdef BEC2_DIAG
          if (auto_ind == diat_ind) then
             if (wrt_bec2_diag_3d(diatsio3uptake_idx_t)) then
                idiag = idx_bec2_diag_3d(diatsio3uptake_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = VSiO3
             endif
          end if
# endif /* BEC2_DIAG */

     !-----------------------------------------------------------------------
     !     get photosynth. rate, phyto C biomass change, photoadapt
     !-----------------------------------------------------------------------

          if (auto_ind==diat_ind) then
             PCmax = PCref(auto_ind) * f_nut * Tfunc_diat
          else
             PCmax = PCref(auto_ind) * f_nut * Tfunc
          endif

          where (temp < temp_thres(auto_ind)) PCmax = c0
          light_lim = c1 - exp((-c1 * alphaPI(auto_ind) * thetaC(istr:iend,auto_ind) * PAR_lay) /
     &                      (PCmax + epsTinv))
          PCphoto = PCmax * light_lim

          photoC(istr:iend,auto_ind) = PCphoto * C_loc(istr:iend,auto_ind)
# ifdef BEC2_DIAG
          if (auto_ind == sp_ind) then
             if (wrt_bec2_diag_3d(splightlim_idx_t)) then
                idiag = idx_bec2_diag_3d(splightlim_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = light_lim
             endif
          else if (auto_ind == diat_ind) then
             if (wrt_bec2_diag_3d(diatlightlim_idx_t)) then
                idiag = idx_bec2_diag_3d(diatlightlim_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = light_lim
             endif
          else
             ! Diazotrophs:
             if (wrt_bec2_diag_3d(diazlightlim_idx_t)) then
                idiag = idx_bec2_diag_3d(diazlightlim_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = light_lim
             endif
          end if
# endif /* BEC2_DIAG */

     !-----------------------------------------------------------------------
     !  Get nutrient uptakes by autotrophs based on calculated C fixation
     !  total N uptake VNC is used in photoadaption
     !-----------------------------------------------------------------------
# ifdef Ncycle_SY
       epsN = 1.0e-30 ! small value to avoid dividing by 0
          where (VNtot(istr:iend,auto_ind) > c0)
             where (NO3_loc .le. c0)
                NO3_V(istr:iend,auto_ind) = c0
             elsewhere
                NO3_V(istr:iend,auto_ind) = (VNO3(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                                        * photoC(istr:iend,auto_ind) * Q
             end where
             where (NO2_loc .le. c0)
                NO2_V(istr:iend,auto_ind) = c0
             elsewhere
                NO2_V(istr:iend,auto_ind) = (VNO2(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                                        * photoC(istr:iend,auto_ind) * Q
             end where
             NH4_V(istr:iend,auto_ind) = (VNH4(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                    * photoC(istr:iend,auto_ind) * Q

             VNC = PCphoto * Q
          elsewhere
             NO3_V(istr:iend,auto_ind) = c0
             NO2_V(istr:iend,auto_ind) = c0
             NH4_V(istr:iend,auto_ind) = c0
             VNC = c0
          end where
# else
          where (VNtot(istr:iend,auto_ind) > c0)
             NO3_V(istr:iend,auto_ind) = (VNO3(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                   * photoC(istr:iend,auto_ind) * Q
             NH4_V(istr:iend,auto_ind) = (VNH4(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                   * photoC(istr:iend,auto_ind) * Q
             VNC = PCphoto * Q
          elsewhere
             NO3_V(istr:iend,auto_ind) = c0
             NH4_V(istr:iend,auto_ind) = c0
             VNC = c0
          end where
# endif
          where (VPtot(istr:iend,auto_ind) > c0)
             PO4_V(istr:iend,auto_ind) = (VPO4(istr:iend,auto_ind) / VPtot(istr:iend,auto_ind))
     &                    * photoC(istr:iend,auto_ind) * Qp(auto_ind)
             DOP_V(istr:iend,auto_ind) = (VDOP(istr:iend,auto_ind) / VPtot(istr:iend,auto_ind))
     &                    * photoC(istr:iend,auto_ind) * Qp(auto_ind)
          elsewhere
             PO4_V(istr:iend,auto_ind) = c0
             DOP_V(istr:iend,auto_ind) = c0
          end where
          photoFe(istr:iend,auto_ind) = photoC(istr:iend,auto_ind) * gQfe(istr:iend,auto_ind)
# ifdef BEC2_DIAG
          if (auto_ind == sp_ind) then
             if (wrt_bec2_diag_3d(ironuptakesp_idx_t)) then
                idiag = idx_bec2_diag_3d(ironuptakesp_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = photoFe(istr:iend,auto_ind)*Hz(istr:iend,j,k)
             endif
          else if (auto_ind == diat_ind) then
             if (wrt_bec2_diag_3d(ironuptakediat_idx_t)) then
                idiag = idx_bec2_diag_3d(ironuptakediat_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = photoFe(istr:iend,auto_ind)*Hz(istr:iend,j,k)
             endif
          else if (auto_ind == diaz_ind) then
             if (wrt_bec2_diag_3d(ironuptakediaz_idx_t)) then
                idiag = idx_bec2_diag_3d(ironuptakediaz_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = photoFe(istr:iend,auto_ind)*Hz(istr:iend,j,k)
             endif
          end if

# endif

     !-----------------------------------------------------------------------
     !  Get nutrient uptake by diatoms based on C fixation
     !-----------------------------------------------------------------------

          if (kSiO3(auto_ind) > c0) then
             photoSi(istr:iend,auto_ind) = photoC(istr:iend,auto_ind) * gQsi(istr:iend,auto_ind)
# ifdef BEC2_DIAG
             if (auto_ind == diat_ind) then
                if (wrt_bec2_diag_3d(diatsiuptake_idx_t)) then
                   idiag = idx_bec2_diag_3d(diatsiuptake_idx_t)
                   bec2_diag_3d(istr:iend,j,k,idiag) = photoSi(istr:iend,auto_ind)*Hz(istr:iend,j,k)
                endif
             end if
# endif
          endif

     !-----------------------------------------------------------------------
     !  calculate pChl, (used in photoadapt., GD98)
     !  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
     !  GD 98 Chl. synth. term
     !-----------------------------------------------------------------------

          WORK1 = alphaPI(auto_ind) * thetaC(istr:iend,auto_ind) * PAR_lay

     ! changed from c0 to value slightly larger than c0 (parameter
     ! PAR_thres_pChl, BEC blew up in SO setup otherwise):
          where (WORK1 > PAR_thres_pChl)
             pChl = thetaN_max(auto_ind) * PCphoto / WORK1
             photoacc(istr:iend,auto_ind) = (pChl * VNC / thetaC(istr:iend,auto_ind))
     &                       * Chl_loc(istr:iend,auto_ind)
          elsewhere
             photoacc(istr:iend,auto_ind) = c0
          end where
# ifdef BEC2_DIAG
          if (wrt_bec2_diag_3d(spphotoacc_idx_t)) then
             idiag = idx_bec2_diag_3d(spphotoacc_idx_t)
             bec2_diag_3d(istr:iend,j,k,idiag) = photoacc(istr:iend,sp_ind)*Hz(istr:iend,j,k)
          endif
          if (wrt_bec2_diag_3d(diatphotoacc_idx_t)) then
             idiag = idx_bec2_diag_3d(diatphotoacc_idx_t)
             bec2_diag_3d(istr:iend,j,k,idiag) = photoacc(istr:iend,diat_ind)*Hz(istr:iend,j,k)
          endif
          if (wrt_bec2_diag_3d(diazphotoacc_idx_t)) then
             idiag = idx_bec2_diag_3d(diazphotoacc_idx_t)
             bec2_diag_3d(istr:iend,j,k,idiag) = photoacc(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
          endif
# endif

     !-----------------------------------------------------------------------
     !  implicit CaCO3 Production, parameterized as function of small phyto production:
     !    - decrease CaCO3 as function of nutrient limitation
     !    - decrease CaCO3 prod at low temperatures
     !    - increase CaCO3 prod under bloom conditions
     !  Maximum calcification rate is 40% of primary production.
     !
     !  explicit CaCO3 Production, parametrized as function of photoC of PFT
     !    think about: temp dependence, f_nut dependence, CaCO3:Corg
     !    ratio (right now set to one)
     !-----------------------------------------------------------------------

          if (imp_calcifier(auto_ind)) then
             CaCO3_PROD(istr:iend,auto_ind) = parm_f_prod_sp_CaCO3 * photoC(istr:iend,auto_ind)
             CaCO3_PROD(istr:iend,auto_ind) = CaCO3_PROD(istr:iend,auto_ind) * f_nut

             where (temp < CaCO3_temp_thres1)
     &          CaCO3_PROD(istr:iend,auto_ind) = CaCO3_PROD(istr:iend,auto_ind) * max((temp-CaCO3_temp_thres2), c0) /
     &                    (CaCO3_temp_thres1-CaCO3_temp_thres2)

             where (C_loc(istr:iend,auto_ind) > CaCO3_sp_thres)
     &          CaCO3_PROD(istr:iend,auto_ind) = min((CaCO3_PROD(istr:iend,auto_ind)
     &                     * C_loc(istr:iend,auto_ind) / CaCO3_sp_thres),
     &                    (f_photosp_CaCO3 * photoC(istr:iend,auto_ind)))
          elseif (exp_calcifier(auto_ind)) then
!!  calcification: use different light_lim than for photosynthesis (see Zondervan 2007)
!             QCaCO3(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind)*light_lim_calc
!             CaCO3_PROD(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind)*photoC(istr:iend,auto_ind)*(1-light_lim)
             CaCO3_PROD(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind)*photoC(istr:iend,auto_ind)
          else
             CaCO3_PROD(istr:iend,auto_ind) = c0
          endif
# ifdef BEC2_DIAG
          if (auto_ind == sp_ind) then
             if (wrt_bec2_diag_3d(caco3prod_idx_t)) then
                idiag = idx_bec2_diag_3d(caco3prod_idx_t)
                bec2_diag_3d(istr:iend,j,k,idiag) = CaCO3_PROD(istr:iend,auto_ind)*Hz(istr:iend,j,k)
             endif
          end if
# endif

     !-----------------------------------------------------------------------
     !  get autotroph loss (in C units)
     !  autotroph agg loss
     !-----------------------------------------------------------------------


       if (auto_ind == diat_ind) then
          auto_loss(istr:iend,auto_ind) =  mort(auto_ind) * Pprime(istr:iend,auto_ind) * Tfunc_diat
       else
          auto_loss(istr:iend,auto_ind) =  mort(auto_ind) * Pprime(istr:iend,auto_ind) * Tfunc
       endif

          auto_agg(istr:iend,auto_ind) = min(agg_rate_max(auto_ind) * dps * Pprime(istr:iend,auto_ind),
     &                               mort2(auto_ind) * Pprime(istr:iend,auto_ind) * Pprime(istr:iend,auto_ind))
          auto_agg(istr:iend,auto_ind) = max((agg_rate_min(auto_ind) * dps) * Pprime(istr:iend,auto_ind),
     &                               auto_agg(istr:iend,auto_ind))

     !-----------------------------------------------------------------------
     !  get grazing rate (graze_sp) on autotroph (in C units)
     !  compute sum of carbon in the grazee class including auto_ind
     !-----------------------------------------------------------------------

          WORK1 = c0

          do auto_ind2 = 1, autotroph_cnt
             if (grazee_ind(auto_ind2) == grazee_ind(auto_ind))
     &          WORK1 = WORK1 + Pprime(istr:iend,auto_ind)
          end do

          z_umax = z_umax_0(auto_ind) * Tfunc
          where (WORK1 > c0)
             auto_graze(istr:iend,auto_ind) = (Pprime(istr:iend,auto_ind) / WORK1) *
     &           z_umax * zooC_loc * (WORK1 / (WORK1 + z_grz(auto_ind)))
          elsewhere
             auto_graze(istr:iend,auto_ind) = c0
          end where

     !-----------------------------------------------------------------------
     !  Get N fixation by diazotrophs based on C fixation,
     !  Diazotrophs fix more than they need then 20% is excreted
     !-----------------------------------------------------------------------
# ifdef Ncycle_SY
          if (Nfixer(auto_ind)) then
             WORK1 = photoC(istr:iend,auto_ind) * Q
             Nfix(istr:iend,auto_ind)     = (WORK1 * r_Nfix_photo) - NO3_V(istr:iend,auto_ind) - NO2_V(istr:iend,auto_ind)
     &                 - NH4_V(istr:iend,auto_ind)
             Nexcrete(istr:iend,auto_ind) = Nfix(istr:iend,auto_ind) + NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind)
     &                 + NH4_V(istr:iend,auto_ind) - WORK1
          endif
# else
          if (Nfixer(auto_ind)) then
             WORK1 = photoC(istr:iend,auto_ind) * Q
             Nfix(istr:iend,auto_ind)     = (WORK1 * r_Nfix_photo) - NO3_V(istr:iend,auto_ind)
     &                 - NH4_V(istr:iend,auto_ind)
             Nexcrete(istr:iend,auto_ind) = Nfix(istr:iend,auto_ind) + NO3_V(istr:iend,auto_ind)
     &                 + NH4_V(istr:iend,auto_ind) - WORK1
          endif
#endif

# ifdef BEC2_DIAG
          if (auto_ind == diaz_ind) then
             if (wrt_bec2_diag_3d(diaznfix_idx_t)) then
                idiag = idx_bec2_diag_3d(diaznfix_idx_t)
                bec2_diag_3d(istr:iend,j,k,diaznfix_idx_t) = Nfix(istr:iend,auto_ind)*Hz(istr:iend,j,k)
             endif
          end if
# endif

     !-----------------------------------------------------------------------
     !  CALCULATE GRAZING AND OTHER MORT
     !-----------------------------------------------------------------------

     !-----------------------------------------------------------------------
     !  routing of grazing and loss terms
     !  all aggregation goes to POC
     !  currently assumes that 33% of grazed caco3 is remineralized
     !  if autotrophs(sp_ind)%graze_zoo ever changes, coefficients on routing grazed sp must change!
     !  min.%C routed to POC from grazing for ballast requirements = 0.4 * Qcaco3
     !  min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho
     !  NOTE: if autotrophs(diat_ind)%graze_zoo is changed, coeff.s for poc,doc and dic must change!
     !-----------------------------------------------------------------------

          auto_graze_zoo(istr:iend,auto_ind) = graze_zoo(auto_ind) * auto_graze(istr:iend,auto_ind)
          if (imp_calcifier(auto_ind)) then
             auto_graze_poc(istr:iend,auto_ind) = auto_graze(istr:iend,auto_ind)
     &       * max((caco3_poc_min * QCaCO3(istr:iend,auto_ind)),
     &             min(spc_poc_fac * max(1.0,Pprime(istr:iend,auto_ind)),
     &                     f_graze_sp_poc_lim))
          else
             auto_graze_poc(istr:iend,auto_ind) = graze_poc(auto_ind) * auto_graze(istr:iend,auto_ind)
          endif
          auto_graze_doc(istr:iend,auto_ind) = graze_doc(auto_ind) * auto_graze(istr:iend,auto_ind)
          auto_graze_dic(istr:iend,auto_ind) = auto_graze(istr:iend,auto_ind)
     &        - (auto_graze_zoo(istr:iend,auto_ind) + auto_graze_poc(istr:iend,auto_ind)
     &        + auto_graze_doc(istr:iend,auto_ind))

          if (imp_calcifier(auto_ind)) then
             auto_loss_poc(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind) * auto_loss(istr:iend,auto_ind)
          else
             auto_loss_poc(istr:iend,auto_ind) = loss_poc(auto_ind) * auto_loss(istr:iend,auto_ind)
          endif
          auto_loss_doc(istr:iend,auto_ind) = (c1 - parm_labile_ratio) * (auto_loss(istr:iend,auto_ind)
     &         - auto_loss_poc(istr:iend,auto_ind))
          auto_loss_dic(istr:iend,auto_ind) = parm_labile_ratio * (auto_loss(istr:iend,auto_ind)
     &        - auto_loss_poc(istr:iend,auto_ind))

     !-----------------------------------------------------------------------
     ! P from some autotrophs w/ Qp different from Qp_zoo_pom must be routed differently than other
     ! elements to ensure that sinking detritus and zooplankton pools get their fixed P/C ratios.
     ! The remaining P is split evenly between DOP and PO4.
     !-----------------------------------------------------------------------

          if (Qp(auto_ind) /= Qp_zoo_pom) then
             remaining_P = ((auto_graze(istr:iend,auto_ind) + auto_loss(istr:iend,auto_ind)
     &                  + auto_agg(istr:iend,auto_ind)) * Qp(auto_ind))
     &                  - ((auto_graze_zoo(istr:iend,auto_ind)) * Qp_zoo_pom)
     &                  - ((auto_graze_poc(istr:iend,auto_ind) + auto_loss_poc(istr:iend,auto_ind)
     &                  + auto_agg(istr:iend,auto_ind)) * Qp_zoo_pom)
             remaining_P_dop(istr:iend,auto_ind) = (c1 - parm_labile_ratio) * remaining_P
             remaining_P_dip(istr:iend,auto_ind) = parm_labile_ratio * remaining_P
          endif

       end do  ! auto_ind = 1, autotroph_cnt

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(grazesp_idx_t)) then
          idiag = idx_bec2_diag_3d(grazesp_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazediat_idx_t)) then
          idiag = idx_bec2_diag_3d(grazediat_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazediaz_idx_t)) then
          idiag = idx_bec2_diag_3d(grazediaz_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(sploss_idx_t)) then
          idiag = idx_bec2_diag_3d(sploss_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_loss(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diatloss_idx_t)) then
          idiag = idx_bec2_diag_3d(diatloss_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_loss(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diazloss_idx_t)) then
          idiag = idx_bec2_diag_3d(diazloss_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_loss(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(spagg_idx_t)) then
          idiag = idx_bec2_diag_3d(spagg_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_agg(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diatagg_idx_t)) then
          idiag = idx_bec2_diag_3d(diatagg_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_agg(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diazagg_idx_t)) then
          idiag = idx_bec2_diag_3d(diazagg_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_agg(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(photocsp_idx_t)) then
          idiag = idx_bec2_diag_3d(photocsp_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = photoC(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(photocdiat_idx_t)) then
          idiag = idx_bec2_diag_3d(photocdiat_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = photoC(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(photocdiaz_idx_t)) then
          idiag = idx_bec2_diag_3d(photocdiaz_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = photoC(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(totprod_idx_t)) then
          idiag = idx_bec2_diag_3d(totprod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = sum(photoC,dim=2)
       endif
       if (wrt_bec2_diag_3d(spnlim_idx_t)) then
          idiag = idx_bec2_diag_3d(spnlim_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = VNtot(istr:iend,sp_ind)
       endif
       if (wrt_bec2_diag_3d(diatnlim_idx_t)) then
          idiag = idx_bec2_diag_3d(diatnlim_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = VNtot(istr:iend,diat_ind)
       endif
       if (wrt_bec2_diag_3d(spplim_idx_t)) then
          idiag = idx_bec2_diag_3d(spplim_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = VPtot(istr:iend,sp_ind)
       endif
       if (wrt_bec2_diag_3d(diatplim_idx_t)) then
          idiag = idx_bec2_diag_3d(diatplim_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = VPtot(istr:iend,diat_ind)
       endif
       if (wrt_bec2_diag_3d(diazplim_idx_t)) then
          idiag = idx_bec2_diag_3d(diazplim_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = VPtot(istr:iend,diaz_ind)
       endif
# ifdef Ncycle_SY
       if (wrt_bec2_diag_3d(spno2uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(spno2uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NO2_V(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diatno2uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(diatno2uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NO2_V(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diazno2uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(diazno2uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NO2_V(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
# endif
       if (wrt_bec2_diag_3d(spno3uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(spno3uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NO3_V(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diatno3uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(diatno3uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NO3_V(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diazno3uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(diazno3uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NO3_V(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(spnh4uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(spnh4uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NH4_V(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diatnh4uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(diatnh4uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NH4_V(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(diaznh4uptake_idx_t)) then
          idiag = idx_bec2_diag_3d(diaznh4uptake_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NH4_V(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazedicsp_idx_t)) then
          idiag = idx_bec2_diag_3d(grazedicsp_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze_dic(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazedicdiat_idx_t)) then
          idiag = idx_bec2_diag_3d(grazedicdiat_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze_dic(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazedicdiaz_idx_t)) then
          idiag = idx_bec2_diag_3d(grazedicdiaz_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze_dic(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(lossdicsp_idx_t)) then
          idiag = idx_bec2_diag_3d(lossdicsp_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_loss_dic(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(lossdicdiat_idx_t)) then
          idiag = idx_bec2_diag_3d(lossdicdiat_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_loss_dic(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(lossdicdiaz_idx_t)) then
          idiag = idx_bec2_diag_3d(lossdicdiaz_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_loss_dic(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazespzoo_idx_t)) then
          idiag = idx_bec2_diag_3d(grazespzoo_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze_zoo(istr:iend,sp_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazediatzoo_idx_t)) then
          idiag = idx_bec2_diag_3d(grazediatzoo_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze_zoo(istr:iend,diat_ind)*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(grazediazzoo_idx_t)) then
          idiag = idx_bec2_diag_3d(grazediazzoo_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = auto_graze_zoo(istr:iend,diaz_ind)*Hz(istr:iend,j,k)
       endif
# endif /* BEC2_DIAG */


     !-----------------------------------------------------------------------
     !  get fractional factor for routing of zoo losses, based on food supply
     !  more material is routed to large detrital pool when diatoms eaten
     !-----------------------------------------------------------------------

       WORK1 = c0
       WORK2 = c0
       do auto_ind = 1, autotroph_cnt
          WORK1 = WORK1 + f_zoo_detr(auto_ind) * (auto_graze(istr:iend,auto_ind) + epsC * epsTinv)
          WORK2 = WORK2 + (auto_graze(istr:iend,auto_ind) + epsC * epsTinv)
       end do
       f_zoo_detr_loc = WORK1 / WORK2

     !-----------------------------------------------------------------------
     !  0.01 small zoo threshold C concentration (mmol C/m^3)
     !  zoo losses, scaled by Tfunc
     !-----------------------------------------------------------------------

       C_loss_thres = f_loss_thres * loss_thres_zoo

       Zprime = max(zooC_loc - C_loss_thres, c0)

       zoo_loss = (parm_z_mort2_0 * Zprime**1.5 + parm_z_mort_0 * Zprime) * Tfunc

       zoo_loss_doc = (c1 - parm_labile_ratio) * (c1 - f_zoo_detr_loc) * zoo_loss
       zoo_loss_dic = parm_labile_ratio * (c1 - f_zoo_detr_loc) * zoo_loss

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(zooloss_idx_t)) then
          idiag = idx_bec2_diag_3d(zooloss_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = zoo_loss*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(zoolossdic_idx_t)) then
          idiag = idx_bec2_diag_3d(zoolossdic_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = zoo_loss_dic*Hz(istr:iend,j,k)
       endif
# endif

     !-----------------------------------------------------------------------
     !  compute terms for DOM
     !-----------------------------------------------------------------------

       DOC_prod = zoo_loss_doc + sum(auto_loss_doc, dim=2) + sum(auto_graze_doc, dim=2)
       DON_prod = Q * DOC_prod
       DOP_prod = Qp_zoo_pom * zoo_loss_doc
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(docprod_idx_t)) then
          idiag = idx_bec2_diag_3d(docprod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DOC_prod*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(donprod_idx_t)) then
          idiag = idx_bec2_diag_3d(donprod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DON_prod*Hz(istr:iend,j,k)
       endif
# endif
       do auto_ind = 1, autotroph_cnt
          if (Qp(auto_ind) == Qp_zoo_pom) then
             DOP_prod = DOP_prod + Qp(auto_ind) * (auto_loss_doc(istr:iend,auto_ind)
     &                  + auto_graze_doc(istr:iend,auto_ind))
          else
             DOP_prod = DOP_prod + remaining_P_dop(istr:iend,auto_ind)
          endif
       end do
       DOFe_prod = Qfe_zoo * zoo_loss_doc
       do auto_ind = 1, autotroph_cnt
          DOFe_prod = DOFe_prod + Qfe(istr:iend,auto_ind) *
     &    (auto_loss_doc(istr:iend,auto_ind) + auto_graze_doc(istr:iend,auto_ind))
       end do
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(dofeprod_idx_t)) then
          idiag = idx_bec2_diag_3d(dofeprod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DOFe_prod*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(dopprod_idx_t)) then
          idiag = idx_bec2_diag_3d(dopprod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DOP_prod*Hz(istr:iend,j,k)
       endif
# endif

       DOC_remin  = DOC_loc  * DOC_reminR
       DON_remin  = DON_loc  * DON_reminR
       DOFe_remin = DOFe_loc * DOFe_reminR
       DOP_remin(istr:iend,j)  = DOP_loc  * DOP_reminR

     !-----------------------------------------------------------------------
     !  Refractory remin rate due to photochemistry
     !  below euphotic zone remin rate sharply decrease
     !-----------------------------------------------------------------------

       where (PAR_lay > 1.0)
          DONr_remin = DONr_loc * DONr_reminR
          DOPr_remin(istr:iend,j) = DOPr_loc * DOPr_reminR
       elsewhere
          DONr_remin = DONr_loc * (c1/(365.0*670.0)) * dps ! 1/670 yrs
          DOPr_remin(istr:iend,j) = DOPr_loc * (c1/(365.0*460.0)) * dps ! 1/460 yrs
          DOC_remin = DOC_remin * 0.0685
          DON_remin = DON_remin * 0.1
          DOFe_remin = DOFe_remin * 0.05
          DOP_remin(istr:iend,j) = DOP_remin(istr:iend,j) * 0.05
       end where
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(docremin_idx_t)) then
          idiag = idx_bec2_diag_3d(docremin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DOC_remin*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(donremin_idx_t)) then
          idiag = idx_bec2_diag_3d(donremin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DON_remin*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(donrremin_idx_t)) then
          idiag = idx_bec2_diag_3d(donrremin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DONr_remin*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(doferemin_idx_t)) then
          idiag = idx_bec2_diag_3d(doferemin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DOFe_remin*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(dopremin_idx_t)) then
          idiag = idx_bec2_diag_3d(dopremin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DOP_remin(istr:iend,j)*Hz(istr:iend,j,k)
       endif
# endif

     !-----------------------------------------------------------------------
     !  large detritus C
     !-----------------------------------------------------------------------

       POC_prod = f_zoo_detr_loc * zoo_loss + sum(auto_graze_poc, dim=2)
     &            + sum(auto_agg, dim=2) + sum(auto_loss_poc, dim=2)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(pocprod_idx_t)) then
          idiag = idx_bec2_diag_3d(pocprod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = POC_prod(istr:iend)*Hz(istr:iend,j,k)
       endif
# endif

     !-----------------------------------------------------------------------
     !  large detrital CaCO3
     !  33% of CaCO3 is remin when phyto are grazed
     !-----------------------------------------------------------------------

       P_CaCO3_prod = ((c1 - f_graze_CaCO3_REMIN) * auto_graze(istr:iend,sp_ind) +
     &                auto_loss(istr:iend,sp_ind) + auto_agg(istr:iend,sp_ind)) * QCaCO3(istr:iend,sp_ind)

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(pcaco3prod_idx_t)) then
          idiag = idx_bec2_diag_3d(pcaco3prod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = P_CaCO3_prod(istr:iend)*Hz(istr:iend,j,k)
       endif
# endif

     !-----------------------------------------------------------------------
     !  large detritus SiO2
     !  grazed diatom SiO2, 60% is remineralized
     !-----------------------------------------------------------------------

       P_SiO2_prod = c0
       do auto_ind = 1, autotroph_cnt
          if (kSiO3(auto_ind) > c0) then
             P_SiO2_prod = P_SiO2_prod + Qsi(istr:iend,auto_ind)
     &          * ((c1 - f_graze_si_remin) * auto_graze(istr:iend,auto_ind) + auto_agg(istr:iend,auto_ind)
     &          + loss_poc(auto_ind) * auto_loss(istr:iend,auto_ind))
          endif
       end do
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(sio2prod_idx_t)) then
          idiag = idx_bec2_diag_3d(sio2prod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = P_SiO2_prod(istr:iend)*Hz(istr:iend,j,k)
       endif
# endif

     !-----------------------------------------------------------------------
     !  Compute iron scavenging (Moore et al. 2008):
     !  1) compute in terms of loss per year per unit iron (%/year/fe)
     !  2) scale by sinking POMx10 + Dust + bSi + CaCO3 flux
     !                      POMx10 = 120 ;
     !  3) increase scavenging at higher iron (>0.6nM) (we are currently using 0.8)
     !  4) convert to net loss per second
     !-----------------------------------------------------------------------

       Fe_scavenge_rate = parm_Fe_scavenge_rate0

       Fe_scavenge_rate = Fe_scavenge_rate *
     &     ((POC_sflux_out(istr:iend,j) + POC_hflux_out(istr:iend,j)) * 120.1 +
     &      (P_CaCO3_sflux_out(istr:iend,j) + P_CaCO3_hflux_out(istr:iend,j)) * P_CaCO3_mass +
     &      (P_SiO2_sflux_out(istr:iend,j) + P_SiO2_hflux_out(istr:iend,j)) * P_SiO2_mass +
     &      (dust_sflux_out(istr:iend,j) + dust_hflux_out(istr:iend,j)) * dust_fescav_scale)

       where (Fe_loc > Fe_scavenge_thres1)
          Fe_scavenge_rate = Fe_scavenge_rate +
     &                    (Fe_loc - Fe_scavenge_thres1) * fe_max_scale2
       end where

       Fe_scavenge = yps * Fe_loc * Fe_scavenge_rate

       P_iron_prod = (zoo_loss * f_zoo_detr_loc * Qfe_zoo) + Fe_scavenge

       do auto_ind = 1, autotroph_cnt
          P_iron_prod = P_iron_prod + Qfe(istr:iend,auto_ind)
     &       * (auto_agg(istr:iend,auto_ind) + auto_graze_poc(istr:iend,auto_ind)
     &          + auto_loss_poc(istr:iend,auto_ind))
       end do
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(pironprod_idx_t)) then
          idiag = idx_bec2_diag_3d(pironprod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = P_iron_prod*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(fescavenge_idx_t)) then
          idiag = idx_bec2_diag_3d(fescavenge_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = Fe_scavenge*Hz(istr:iend,j,k)
       endif
       if (wrt_bec2_diag_3d(fescavengerate_idx_t)) then
          idiag = idx_bec2_diag_3d(fescavengerate_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = Fe_scavenge_rate
       endif
# endif

       call compute_particulate_terms(k,QA_dust_def, TEMP, O2_loc, NO3_loc,POC_prod,
     &        P_iron_prod,P_CaCO3_prod,P_SiO2_prod,SED_DENITRIF, OTHER_REMIN,istr,iend,j)

     !-----------------------------------------------------------------------
     !  nitrate & ammonium
     !  nitrification in low light
     !  use exponential decay of PAR across model level to compute taper factor
     !-----------------------------------------------------------------------

# ifdef Ncycle_SY
        ! ammonium oxidation
#if defined DAILYPAR_PHOTOINHIBITION
        where (NH4_loc .gt. c0 .and. O2_loc .gt. c0 .and. PAR_out_AVG < parm_nitrif_par_lim)
#else
        where (NH4_loc .gt. c0 .and. O2_loc .gt. c0 .and. PAR_out < parm_nitrif_par_lim)
#endif
            ammox = parm_kao * (O2_loc / (parm_ko2_ao + O2_loc))
     &      * (NH4_loc/(parm_knh4_ao + NH4_loc))
#if defined DAILYPAR_PHOTOINHIBITION
                where (PAR_in_AVG > parm_nitrif_par_lim)
                        ammox = ammox * log(PAR_out_AVG / parm_nitrif_par_lim) / (-KPARdz)
                endwhere
#else
                where (PAR_in > parm_nitrif_par_lim)
                        ammox = ammox * log(PAR_out / parm_nitrif_par_lim) / (-KPARdz)
                endwhere
#endif
        elsewhere
            ammox = c0
        end where

#if defined DAILYPAR_PHOTOINHIBITION
        where (NO2_loc .gt. c0 .and. O2_loc .gt. c0 .and. PAR_out_AVG < parm_nitrif_par_lim)
#else
        where (NO2_loc .gt. c0 .and. O2_loc .gt. c0 .and. PAR_out < parm_nitrif_par_lim)
#endif
               nitrox = parm_kno * (O2_loc / (parm_ko2_no + O2_loc))
     &         * (NO2_loc/(parm_kno2_no + NO2_loc))
#if defined DAILYPAR_PHOTOINHIBITION
               where (PAR_in_AVG > parm_nitrif_par_lim)
                        nitrox = nitrox * log(PAR_out_AVG /parm_nitrif_par_lim) / (-KPARdz)
               endwhere
#else
               where (PAR_in > parm_nitrif_par_lim)
                        nitrox = nitrox * log(PAR_out /parm_nitrif_par_lim) / (-KPARdz)
               endwhere
#endif
        elsewhere
            nitrox = c0
        end where

        where (O2_loc .gt. 1e-5)
                n2o_nh4_yield1 = c1/(2.0 + 2.0 / ((parm_n2o_ji_a / O2_loc + parm_n2o_ji_b)/100.0));
                no2_nh4_yield1 = c1/(1.0 * ((parm_n2o_ji_a / O2_loc +parm_n2o_ji_b)/100.0) + 1.0);
        elsewhere
                n2o_nh4_yield1 = 0.5
                no2_nh4_yield1 = c0
        end where

# ifdef BEC2_DIAG
      if (wrt_bec2_diag_3d(ammox_idx_t)) then
          idiag = idx_bec2_diag_3d(ammox_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = ammox*Hz(istr:iend,j,k)
      endif
      if (wrt_bec2_diag_3d(nitrox_idx_t)) then
          idiag = idx_bec2_diag_3d(nitrox_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = nitrox*Hz(istr:iend,j,k)
      endif
# endif

# else
       WORK1 = (O2_loc - parm_o2_min) / parm_o2_min_delta
       WORK1 = min(max(WORK1,c0),c1)
#if defined DAILYPAR_PHOTOINHIBITION
       where (PAR_out_AVG < parm_nitrif_par_lim)
#else
       where (PAR_out < parm_nitrif_par_lim)
#endif
          NITRIF = WORK1 * parm_kappa_nitrif * NH4_loc
#if defined DAILYPAR_PHOTOINHIBITION
          where (PAR_in_AVG > parm_nitrif_par_lim)
             NITRIF = NITRIF * log(PAR_out_AVG / parm_nitrif_par_lim) / (-KPARdz)
          end where
#else
          where (PAR_in > parm_nitrif_par_lim)
             NITRIF = NITRIF * log(PAR_out / parm_nitrif_par_lim) / (-KPARdz)
          end where
#endif
       elsewhere
          NITRIF = c0
       end where

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(nitrif_idx_t)) then
          idiag = idx_bec2_diag_3d(nitrif_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = NITRIF*Hz(istr:iend,j,k)
       endif
# endif
# endif

     !-----------------------------------------------------------------------
     !  Compute denitrification under low O2 conditions
     !-----------------------------------------------------------------------

       WORK1 = ((parm_o2_min + parm_o2_min_delta) - O2_loc) / parm_o2_min_delta
       WORK1 = min(max(WORK1,c0),c1)

       WORK1 = merge(c0, WORK1, NO3_loc == c0)

# ifdef Ncycle_SY
! split denitrification into three components i.e., one done by no3, one by no2 and one by n2o
       ! Total organic matter denitrified in C units
       R_oxic = parm_koxic * (O2_loc / (parm_ko2_oxic + O2_loc))
       R_den1 = parm_kden1 * (NO3_loc / (parm_kno3_den1 + NO3_loc)) * exp(-max(O2_loc,c0)/parm_ko2_den1) + epsN
       R_den2 = parm_kden2 * (NO2_loc / (parm_kno2_den2 + NO2_loc)) * exp(-max(O2_loc,c0)/parm_ko2_den2) + epsN
       R_den3 = parm_kden3 * (N2O_loc / (parm_kn2o_den3 + N2O_loc)) * exp(-max(O2_loc,c0)/parm_ko2_den3) + epsN
       roxic = R_oxic / (R_den1 + R_den2 + R_den3 + R_oxic)
       rden1 = R_den1 / (R_den1 + R_den2 + R_den3 + R_oxic)
       rden2 = R_den2 / (R_den1 + R_den2 + R_den3 + R_oxic)
       rden3 = R_den3 / (R_den1 + R_den2 + R_den3 + R_oxic)

       Rden =  (DOC_remin + POC_remin(istr:iend,j)
     &                - OTHER_REMIN - SED_DENITRIF * denitrif_C_N)

! Calculate denitrification rates
! --> if  consumption is such that NO3, NO2 or N2O goes negative after one
! timestep then reduce consumption such that it is completely consumed

! ANAMMOX NO2 + NH4 --> N2
       where  (NH4_loc .gt. c0 .and. NO2_loc .gt. c0)
                anammox = parm_kax * (NH4_loc / (parm_knh4_ax + NH4_loc)) *
     &                  (NO2_loc / (parm_kno2_ax + NO2_loc)) * exp(-max(O2_loc,c0)/parm_ko2_ax)
        elsewhere
                anammox = c0
        end where

! DENITRIF1 NO3 --> NO2
       where (NO3_loc + (nitrox - rden1 * Rden * denitrif_NO3_C -
     &         SED_DENITRIF - sum(NO3_V, dim=2)) * dt .ge. c0)
                DENITRIF1 = rden1 * Rden * denitrif_NO3_C
       elsewhere (NO3_loc .le. c0)
                DENITRIF1 =  c0
       elsewhere
                DENITRIF1 = NO3_loc / dt + nitrox - SED_DENITRIF -
     &          sum(NO3_V, dim=2)
       end where

! DENITRIF2 NO2 -->0.5  N2O
       where (NO2_loc  +  (ammox * no2_nh4_yield1 - nitrox +
     &        DENITRIF1 - rden2 * Rden * denitrif_NO2_C - anammox - sum(NO2_V, dim=2)) *
     &        dt .ge. c0)
                DENITRIF2 = rden2 * Rden * denitrif_NO2_C
       elsewhere (NO2_loc .le. c0)
                DENITRIF2 = c0
       elsewhere
                DENITRIF2 = NO2_loc / dt + ammox * no2_nh4_yield1 - nitrox +
     &          DENITRIF1
       end where

! DENITRIF3 N2O --> N2
       where (N2O_loc + (ammox * n2o_nh4_yield1 + 0.5 * DENITRIF2
     &        - rden3 * Rden * denitrif_N2O_C) * dt  .ge. c0)
                        DENITRIF3 = rden3 * Rden * denitrif_N2O_C
        elsewhere (N2O_loc .le. c0)
                DENITRIF3 = c0
        elsewhere
                DENITRIF3 = N2O_loc / dt + ammox * n2o_nh4_yield1 + 0.5 * DENITRIF2
        end where

# ifdef BEC2_DIAG
      if (wrt_bec2_diag_3d(anammox_idx_t)) then
          idiag = idx_bec2_diag_3d(anammox_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = anammox*Hz(istr:iend,j,k)
      endif
      if (wrt_bec2_diag_3d(denitrif1_idx_t)) then
          idiag = idx_bec2_diag_3d(denitrif1_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DENITRIF1*Hz(istr:iend,j,k)
      endif
      if (wrt_bec2_diag_3d(denitrif2_idx_t)) then
          idiag = idx_bec2_diag_3d(denitrif2_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DENITRIF2*Hz(istr:iend,j,k)
      endif
      if (wrt_bec2_diag_3d(denitrif3_idx_t)) then
          idiag = idx_bec2_diag_3d(denitrif3_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DENITRIF3*Hz(istr:iend,j,k)
      endif
# endif

# else

       DENITRIF = WORK1 * ((DOC_remin
     &                  + POC_remin(istr:iend,j)
     &                  - OTHER_REMIN) / denitrif_C_N - SED_DENITRIF)

       where (NO3_loc < parm_denitrif_NO3_limit)
            DENITRIF = 1/parm_denitrif_NO3_limit * NO3_loc * DENITRIF
       end where

# ifdef NO3_RESTORE
       if (k .eq. N)
             restore_no3 = c1/2592000 * (no3_restore(istr:iend,j,k)-NO3_loc)
       else
             restore_no3 = c0
# endif

# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(denitrif_idx_t)) then
          idiag = idx_bec2_diag_3d(denitrif_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = DENITRIF*Hz(istr:iend,j,k)
       endif
# endif
# endif

     !-----------------------------------------------------------------------
     !  Collect tendencies of biological processes in DTRACER_MODULE
     !-----------------------------------------------------------------------

     !-----------------------------------------------------------------------
     !  nitrate & ammonium
     !-----------------------------------------------------------------------
# ifdef Ncycle_SY
       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNO3)) = nitrox - DENITRIF1 -
     &                                 SED_DENITRIF - sum(NO3_V, dim=2)+
     &                                 (r_no2tonh4_ax - c1) * anammox
# ifdef NO3_RESTORE
     &                                 + restore_no3
# endif

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNO2)) = ammox * no2_nh4_yield1 - nitrox +
     &                                 DENITRIF1 - DENITRIF2 - r_no2tonh4_ax * anammox - sum(NO2_V, dim=2)

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNH4)) = -sum(NH4_V, dim=2) - ammox - anammox + DON_remin + DONr_remin
     &     + Q * (zoo_loss_dic + sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2)
     &     + POC_remin(istr:iend,j)
     &     * (c1 - DONrefract))

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iN2O)) = ammox * n2o_nh4_yield1 + 0.5 * DENITRIF2 - DENITRIF3

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iN2)) = DENITRIF3 + 0.5 * SED_DENITRIF + anammox

# else
       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNO3)) = NITRIF - DENITRIF -
     &                                 SED_DENITRIF - sum(NO3_V, dim=2)
# ifdef NO3_RESTORE
     &                                 + restore_no3
# endif

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNH4)) = -sum(NH4_V, dim=2) - NITRIF + DON_remin + DONr_remin
     &     + Q * (zoo_loss_dic + sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2)
     &     + POC_remin(istr:iend,j)
     &     * (c1 - DONrefract))
# endif

       do auto_ind = 1, autotroph_cnt
          if (Nfixer(auto_ind))
     &        DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNH4)) = DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNH4))
     &                                          + Nexcrete(istr:iend,auto_ind)
       end do

     !-----------------------------------------------------------------------
     !  dissolved iron
     !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iFE)) = P_iron_remin(istr:iend,j)
     &           + (Qfe_zoo * zoo_loss_dic) + DOFe_remin - sum(photoFe, dim=2) - Fe_scavenge

       do auto_ind = 1, autotroph_cnt
          DTRACER_MODULE(istr:iend,j,k,bgc_idx(iFE)) = DTRACER_MODULE(istr:iend,j,k,bgc_idx(iFE))
     &      + (Qfe(istr:iend,auto_ind) * (auto_loss_dic(istr:iend,auto_ind) + auto_graze_dic(istr:iend,auto_ind)))
     &      + auto_graze_zoo(istr:iend,auto_ind) * (Qfe(istr:iend,auto_ind)-Qfe_zoo)
       end do

     !-----------------------------------------------------------------------
     !  dissolved SiO3
     !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iSIO3)) = P_SiO2_remin(istr:iend,j)

       do auto_ind = 1, autotroph_cnt
          if (kSiO3(auto_ind) > 0) then
             DTRACER_MODULE(istr:iend,j,k,bgc_idx(iSIO3)) = DTRACER_MODULE(istr:iend,j,k,bgc_idx(iSIO3))
     &                 - photoSi(istr:iend,auto_ind) + Qsi(istr:iend,auto_ind)
     &                 * (f_graze_si_remin * auto_graze(istr:iend,auto_ind)
     &                 + (c1 - loss_poc(auto_ind)) * auto_loss(istr:iend,auto_ind))
          endif
       end do

     !-----------------------------------------------------------------------
     !  phosphate
     !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iPO4)) = DOP_remin(istr:iend,j) + DOPr_remin(istr:iend,j)
     &     - sum(PO4_V, dim=2)
     &     + Qp_zoo_pom * ((c1 - DOPrefract)
     &     * POC_remin(istr:iend,j)
     &     + zoo_loss_dic)

       do auto_ind = 1, autotroph_cnt
          if (Qp(auto_ind) == Qp_zoo_pom) then
             DTRACER_MODULE(istr:iend,j,k,bgc_idx(iPO4)) = DTRACER_MODULE(istr:iend,j,k,bgc_idx(iPO4))
     &       + Qp(auto_ind) * (auto_loss_dic(istr:iend,auto_ind) + auto_graze_dic(istr:iend,auto_ind))
          else
             DTRACER_MODULE(istr:iend,j,k,bgc_idx(iPO4)) = DTRACER_MODULE(istr:iend,j,k,bgc_idx(iPO4))
     &                  + remaining_P_dip(istr:iend,auto_ind)
          endif
       end do


     !-----------------------------------------------------------------------
     !  autotroph Carbon
     !  autotroph Chlorophyll
     !  autotroph Fe
     !  autotroph Si
     !  autotroph CaCO3
     !-----------------------------------------------------------------------

       do auto_ind = 1, autotroph_cnt
          WORK1 = auto_graze(istr:iend,auto_ind) + auto_loss(istr:iend,auto_ind) + auto_agg(istr:iend,auto_ind)

          n_loc = C_ind(auto_ind)
          DTRACER_MODULE(istr:iend,j,k,n_loc) = photoC(istr:iend,auto_ind) - WORK1

          n_loc = Chl_ind(auto_ind)
          DTRACER_MODULE(istr:iend,j,k,n_loc) = photoacc(istr:iend,auto_ind) - thetaC(istr:iend,auto_ind) * WORK1

          n_loc = Fe_ind(auto_ind)
          DTRACER_MODULE(istr:iend,j,k,n_loc) =  photoFe(istr:iend,auto_ind) - Qfe(istr:iend,auto_ind) * WORK1

          n_loc = Si_ind(auto_ind)
          if (n_loc > 0) then
             DTRACER_MODULE(istr:iend,j,k,n_loc) =  photoSi(istr:iend,auto_ind) - Qsi(istr:iend,auto_ind) * WORK1
          endif

          n_loc = CaCO3_ind(auto_ind)
          if (n_loc > 0) then
             DTRACER_MODULE(istr:iend,j,k,n_loc) = CaCO3_PROD(istr:iend,auto_ind)
     &              - QCaCO3(istr:iend,auto_ind) * WORK1
          endif
       end do

     !-----------------------------------------------------------------------
     !  zoo Carbon
     !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iZOOC)) = sum(auto_graze_zoo, dim=2) - zoo_loss

     !-----------------------------------------------------------------------
     !  dissolved organic matter
     !  from sinking remin small fraction to refractory pool
     !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDOC)) = DOC_prod - DOC_remin

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDON)) = (DON_prod * (c1 - DONrefract)) - DON_remin

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDONR)) = (DON_prod * DONrefract) - DONr_remin
     &         + ( POC_remin(istr:iend,j)
     &         * DONrefract * Q)

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDOP)) = (DOP_prod * (c1 - DOPrefract)) - DOP_remin(istr:iend,j)
     &         - sum(DOP_V, dim=2)

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDOPR)) = (DOP_prod * DOPrefract) - DOPr_remin(istr:iend,j)
     &         + ( POC_remin(istr:iend,j)
     &         * DOPrefract * Qp_zoo_pom)

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDOFE)) = DOFe_prod - DOFe_remin

     !-----------------------------------------------------------------------
     !   dissolved inorganic Carbon
     !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDIC)) =
     &       sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2) - sum(photoC, dim=2)
     &       + DOC_remin
     &       + POC_remin(istr:iend,j) + P_CaCO3_remin(istr:iend,j)
     &       + zoo_loss_dic

       do auto_ind = 1, autotroph_cnt
          if (CaCO3_ind(auto_ind) > 0)
     &        DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDIC)) = DTRACER_MODULE(istr:iend,j,k,bgc_idx(iDIC))
     &        + f_graze_CaCO3_REMIN * auto_graze(istr:iend,auto_ind) * QCaCO3(istr:iend,auto_ind)
     &        - CaCO3_PROD(istr:iend,auto_ind)
       end do

     !-----------------------------------------------------------------------
     !  alkalinity
     !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iALK)) = -DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNO3)) +
     &                  DTRACER_MODULE(istr:iend,j,k,bgc_idx(iNH4)) + c2
     &     * P_CaCO3_remin(istr:iend,j)

       do auto_ind = 1, autotroph_cnt
          if (CaCO3_ind(auto_ind) > 0)
     &        DTRACER_MODULE(istr:iend,j,k,bgc_idx(iALK)) = DTRACER_MODULE(istr:iend,j,k,bgc_idx(iALK))
     &        + c2 * (f_graze_CaCO3_REMIN * auto_graze(istr:iend,auto_ind) * QCaCO3(istr:iend,auto_ind)
     &        - CaCO3_PROD(istr:iend,auto_ind))
       end do

     !-----------------------------------------------------------------------
     !  oxygen
     !-----------------------------------------------------------------------

       O2_PRODUCTION = c0
# ifdef Ncycle_SY
       do auto_ind = 1, autotroph_cnt
          if (.not. Nfixer(auto_ind)) then
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &          + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &          + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &          + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind))) / parm_Remin_D_C_O2
     &          + (NO2_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &          + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind))) / parm_Red_D_C_O2_NO2V)
             end where
          else
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &             + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &           + (NO2_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind)
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2_NO2V
     &           + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind)
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Remin_D_C_O2
     &           + (Nfix(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind)
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2_diaz)
             end where
          endif
       end do
#else
       do auto_ind = 1, autotroph_cnt
          if (.not. Nfixer(auto_ind)) then
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &          + NH4_V(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &          + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &          + NH4_V(istr:iend,auto_ind))) / parm_Remin_D_C_O2)
             end where
          else
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &             + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &           + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Remin_D_C_O2
     &           + (Nfix(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind)
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2_diaz)
             end where
          endif
       end do
# endif

# ifdef Ncycle_SY
       where (O2_loc .le. c0)
                WORK1 = c0
       elsewhere
       WORK1 = roxic
       endwhere
# else
       WORK1 = (O2_loc - parm_o2_min) / parm_o2_min_delta
       WORK1 = min(max(WORK1,c0),c1)
# endif
# ifdef Ncycle_SY
       O2_CONSUMPTION = WORK1 *
     &  ((
     &   POC_remin(istr:iend,j)
     &   + DOC_remin - (SED_DENITRIF*denitrif_C_N) - OTHER_REMIN)/ parm_Remin_D_C_O2) + (zoo_loss_dic
     &   + sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2)) / parm_Remin_D_C_O2
     &   + 1.5 * ammox + 0.5 * nitrox
# else
       O2_CONSUMPTION = WORK1 *
     &  ((
     &   POC_remin(istr:iend,j)
     &   + DOC_remin - (SED_DENITRIF*denitrif_C_N) - OTHER_REMIN + zoo_loss_dic
     &   + sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2)) / parm_Remin_D_C_O2 + (c2 * NITRIF))
# endif
       DTRACER_MODULE(istr:iend,j,k,bgc_idx(iO2)) = O2_PRODUCTION - O2_CONSUMPTION

# ifdef BEC2_DIAG
      if (wrt_bec2_diag_3d(o2prod_idx_t)) then
          idiag = idx_bec2_diag_3d(o2prod_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = O2_PRODUCTION*Hz(istr:iend,j,k)
      endif
      if (wrt_bec2_diag_3d(o2cons_idx_t)) then
          idiag = idx_bec2_diag_3d(o2cons_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = O2_CONSUMPTION*Hz(istr:iend,j,k)
      endif
# endif

       return
       end subroutine ecosys_set_interior


!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE INIT_PARTICULATE_TERMS   ****************
!
!-------------------------------------------------------------------------

        SUBROUTINE init_particulate_terms(QA_dust_def,istr,
     &      iend,j,net_dust_in)

        use param
        use ocean_vars
        use scalars
        use bgc_shared_vars
        use bec2_params

        implicit none

     !---------------------------------------------------------------------------
     !   Set incoming fluxes (put into outgoing flux for first level usage).
     !   Set dissolution length, production fraction and mass terms.
     !---------------------------------------------------------------------------

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------
        integer istr, iend, j
        REAL, DIMENSION(istr:iend), INTENT(OUT) ::
     &          QA_dust_def     ! incoming deficit in the QA(dust) POC flux

     !---------------------------------------------------------------------------
     !   local variables
     !---------------------------------------------------------------------------

     ! the incoming dust flux has the same unit as in the ROMS forcing file.
     ! It is assumed to be kg/m2/s.

        REAL, DIMENSION(istr:iend) ::
     &         net_dust_in        ! net incoming dust flux

     !-----------------------------------------------------------------------
     !  parameters, from Armstrong et al. 2000
     !
     !  July 2002, length scale for excess POC and bSI modified by temperature
     !  Value given here is at Tref of 30 deg. C, JKM
     !
     !   diss       dissolution length for soft subclass
     !   gamma      fraction of production -> hard subclass
     !   mass       mass of 1e6 base units in kg    (WAS: 1e9 base units in g)
     !   rho        QA mass ratio of POC to this particle class
     !
     !   Base units:
     !     POC:        mmol C      (WAS: nmol C)
     !     P_CaCO3:    mmol CaCO3  (WAS: nmol CaCO3)
     !     P_SiO2:     mmol SiO2   (WAS: nmol SiO2)
     !     dust:       kg dust     (WAS: g dust)
     !     P_iron:     mmol Fe     (WAS: nmol Fe)
     !
     !  Units of fluxes:
     !     sflux_in:    incoming flux of soft subclass (base units/m^2/sec)
     !     hflux_in:    incoming flux of hard subclass (base units/m^2/sec)
     !     prod:        production term (base units/m^3/sec)
     !     sflux_out:   outgoing flux of soft subclass (base units/m^2/sec)
     !     hflux_out:   outgoing flux of hard subclass (base units/m^2/sec)
     !     remin:       remineralization term (base units/m^3/sec)
     !    NOTE: Area/volume units were cm^2 and cm^3!
     !-----------------------------------------------------------------------

        POC_diss      = parm_POC_diss   ! diss. length (m), modified by TEMP
        POC_mass      = 12.01           ! molecular weight of POC
        POC_gamma     = c0

        P_CaCO3_diss  = parm_CaCO3_diss ! diss. length (m)
        P_CaCO3_gamma = 0.30            ! prod frac -> hard subclass
        P_CaCO3_mass  = 100.09          ! molecular weight of CaCO3
        P_CaCO3_rho   = 0.05 * P_CaCO3_mass / POC_mass ! QA mass ratio for CaCO3

        P_SiO2_diss   = parm_SiO2_diss  ! diss. length (m), modified by TEMP
        P_SiO2_gamma  = 0.030           ! prod frac -> hard subclass
        P_SiO2_mass   = 60.08           ! molecular weight of SiO2
        P_SiO2_rho    = 0.05 * P_SiO2_mass / POC_mass ! QA mass ratio for SiO2

        dust_diss     = 200.0           ! diss. length (m) (DL: changed from cm)
        dust_gamma    = 0.97            ! prod frac -> hard subclass
        dust_mass     = 1.0e6           ! base units are already kg
        dust_rho      = 0.05 * dust_mass / POC_mass ! QA mass ratio for dust

        P_iron_gamma  = c0              ! prod frac -> hard subclass

     !-----------------------------------------------------------------------
     !  Set incoming fluxes
     !-----------------------------------------------------------------------

        P_CaCO3_sflux_out = c0
        P_CaCO3_hflux_out = c0
        P_SiO2_sflux_out  = c0
        P_SiO2_hflux_out  = c0

        dust_sflux_out(istr:iend,j) = (c1 - dust_gamma) * net_dust_in
        dust_hflux_out(istr:iend,j) = dust_gamma * net_dust_in

        P_iron_sflux_out = c0
        P_iron_hflux_out = c0

     !-----------------------------------------------------------------------
     !  Hard POC is QA flux and soft POC is excess POC.
     !
     !  Note (MF):
     !  These names are convenient given the mineral associated soft and
     !  hard components, but rather confusing when reading about the
     !  particulate/ballast models in Armstrong et al 2002, and Lima et
     !  al 2014 (there is no such thing as POC hard, POC_gamma
     !  doesnt exist here or is set to c0 in CESM BEC). Soft/Hard actually
     !  applies only to the mineral associated components. Think about
     !  POC_hflux as being the QA component of POC and POC_sflux the excess POC.
     !-----------------------------------------------------------------------

        POC_sflux_out = c0
        POC_hflux_out = c0

     !-----------------------------------------------------------------------
     !  Compute initial QA(dust) POC flux deficit.
     !      Note : dust_rho converts dust_flux from kg/m2/s to mmolC/m2/s
     !-----------------------------------------------------------------------

        QA_dust_def = dust_rho * (dust_sflux_out(istr:iend,j) + dust_hflux_out(istr:iend,j))

     !-----------------------------------------------------------------------
     !  Apply dust forcing to surface layer (only if explicitly tracking
     !  dust as a tracer)
     !-----------------------------------------------------------------------

        return
        end subroutine init_particulate_terms


!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE COMPUTE_PARTICULATE_TERMS   ****************
!
!-------------------------------------------------------------------------

       subroutine compute_particulate_terms(k,QA_dust_def, TEMP, O2_loc, NO3_loc,
     &        POC_prod,P_iron_prod,P_CaCO3_prod,P_SiO2_prod,SED_DENITRIF,
     &        OTHER_REMIN,istr,iend,j)

       use param
       use grid
       use ocean_vars
       use bgc_shared_vars
       use bec2_params

       implicit none

       integer istr, iend, i, j, k, n_loc, iB
       integer idiag
       real, dimension(istr:iend) ::
     &   temp, O2_loc, NO3_loc
     &   ,POC_prod,P_CaCO3_prod,P_SiO2_prod,P_iron_prod
     &   ,QA_dust_def   ! incoming deficit in the QA (dust) POC flux
     &   ,dust_remin    ! remineralization of dust (base units/m^3/sec)
     &   ,SED_DENITRIF  ! sedimentary denitrification (mmolN/m^3/s)
     &   ,OTHER_REMIN   ! sedimentary remin not due to oxic or denitrification
     &   ,flux_oxidated ! carbon oxidation rate in sediment
     &   ,FESEDFLUX     ! sedimentary Fe inputs
     &   ,scalelength
     &   ,TfuncS         ! temperature scaling from soft POM remin (right now just applied to ballast)
#ifdef TDEP_REMIN
     &   ,Tfunc_soft     ! temperature scaling from soft POM remin (Laufkoetter 2017)
#endif
     &   ,DECAY_Hard       ! scaling factor for dissolution of Hard Ballast
     &   ,DECAY_HardDust   ! scaling factor for dissolution of Hard dust
       real
     &   poc_diss_loc,      ! diss. length used (m)
     &   sio2_diss_loc,     ! diss. length varies spatially with O2 (m)
     &   caco3_diss_loc,
     &   dust_diss_loc,
     &   decay_POC_E,       ! scaling factor for dissolution of excess POC
     &   decay_SiO2,        ! scaling factor for dissolution of SiO2
     &   decay_CaCO3,       ! scaling factor for dissolution of CaCO3
     &   decay_dust,        ! scaling factor for dissolution of dust
     &   POC_PROD_avail,    ! POC production available for excess POC flux (mmol/m^3/s)
     &   new_QA_dust_def,   ! outgoing deficit in the QA(dust) POC flux
     &   flux, flux_alt,    ! temp variables used to update sinking flux
     &   dz_loc, dzr_loc    ! Hz and its inverse at a particular i,j,k location

       logical
     &   poc_error      ! POC error flag

     !-----------------------------------------------------------------------
     !  incoming fluxes are outgoing fluxes from previous level
     !-----------------------------------------------------------------------

       P_CaCO3_sflux_in(istr:iend,j) = P_CaCO3_sflux_out(istr:iend,j)
       P_CaCO3_hflux_in(istr:iend,j) = P_CaCO3_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(caco3fluxin_idx_t)) then
          idiag = idx_bec2_diag_3d(caco3fluxin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = P_CaCO3_sflux_in(istr:iend,j)
     &          + P_CaCO3_hflux_in(istr:iend,j)
       endif
# endif

       P_SiO2_sflux_in(istr:iend,j) = P_SiO2_sflux_out(istr:iend,j)
       P_SiO2_hflux_in(istr:iend,j) = P_SiO2_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(sio2fluxin_idx_t)) then
          idiag = idx_bec2_diag_3d(sio2fluxin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = P_SiO2_sflux_in(istr:iend,j)
     &          + P_SiO2_hflux_in(istr:iend,j)
       endif
# endif

       dust_sflux_in(istr:iend,j) = dust_sflux_out(istr:iend,j)
       dust_hflux_in(istr:iend,j) = dust_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(dustfluxin_idx_t)) then
          idiag = idx_bec2_diag_3d(dustfluxin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = dust_sflux_in(istr:iend,j)
     &          + dust_hflux_in(istr:iend,j)
       endif
# endif

       POC_sflux_in(istr:iend,j) = POC_sflux_out(istr:iend,j)
       POC_hflux_in(istr:iend,j) = POC_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(pocfluxin_idx_t)) then
          idiag = idx_bec2_diag_3d(pocfluxin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = POC_sflux_in(istr:iend,j)
     &          + POC_hflux_in(istr:iend,j)
       endif
# endif

       P_iron_sflux_in(istr:iend,j) = P_iron_sflux_out(istr:iend,j)
       P_iron_hflux_in(istr:iend,j) = P_iron_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
       if (wrt_bec2_diag_3d(pironfluxin_idx_t)) then
          idiag = idx_bec2_diag_3d(pironfluxin_idx_t)
          bec2_diag_3d(istr:iend,j,k,idiag) = P_iron_sflux_in(istr:iend,j)
     &          + P_iron_hflux_in(istr:iend,j)
       endif
# endif

     !-----------------------------------------------------------------------
     !  initialize loss to sediments = 0 and local copy of percent sed
     !-----------------------------------------------------------------------

       P_iron_sed_loss(istr:iend,j) = c0
       POC_sed_loss(istr:iend,j) = c0
       P_CaCO3_sed_loss(istr:iend,j) = c0
       P_SiO2_sed_loss(istr:iend,j) = c0
       dust_sed_loss(istr:iend,j) = c0
       SED_DENITRIF(istr:iend) = c0
       OTHER_REMIN(istr:iend) = c0
       FESEDFLUX(istr:iend) = c0

     !-----------------------------------------------------------------------
     !  compute scalelength and decay factors
     !-----------------------------------------------------------------------

       do i = istr, iend
          if (-z_w(i,j,k) < parm_scalelen_z(1)) then
             scalelength(i) = parm_scalelen_vals(1)
          else if (-z_w(i,j,k) >= parm_scalelen_z(size(parm_scalelen_z))) then
             scalelength(i) = parm_scalelen_vals(size(parm_scalelen_z))
          else
             do n_loc = 2, size(parm_scalelen_z)
                if (-z_w(i,j,k) < parm_scalelen_z(n_loc)) then
                   scalelength(i) = parm_scalelen_vals(n_loc-1)
     &                 + (parm_scalelen_vals(n_loc) - parm_scalelen_vals(n_loc-1))
     &                 * (-z_w(i,j,k) - parm_scalelen_z(n_loc-1))/
     &                 (parm_scalelen_z(n_loc) - parm_scalelen_z(n_loc-1))
                   exit
                endif
             end do
          endif
          DECAY_Hard(i)     = exp(-Hz(i,j,k) / 4.0e4)  ! DL: denominator converted to m
          DECAY_HardDust(i) = exp(-Hz(i,j,k) / 1.2e5)  ! DL: denominator converted to m

     !----------------------------------------------------------------------
     !   Tref = 30.0 reference temperature (deg. C)
     !-----------------------------------------------------------------------

       TfuncS(i) = Q_10**((temp(i) - Tref) / c10)
#ifdef TDEP_REMIN
       Tfunc_soft(i) = exp(parm_ktfunc_soft*(temp(i)-Tref))
#endif

       end do

       poc_error = .false.

       do i = istr, iend
          dz_loc = Hz(i,j,k)
          if (LANDMASK(i,j)) then
             dzr_loc = c1 / dz_loc

             poc_diss_loc = POC_diss
             sio2_diss_loc = P_SiO2_diss
             caco3_diss_loc = P_CaCO3_diss
             dust_diss_loc = dust_diss

     !-----------------------------------------------------------------------
     !  increase POC diss length scale where O2 concentrations are low
     !-----------------------------------------------------------------------

             if ((O2_loc(i) >= 5.0) .and. (O2_loc(i) < 40.0)) then
                poc_diss_loc = POC_diss*(c1+(parm_lowo2_remin_factor-c1)*(40.0 - O2_loc(i))/35.0)
             else if (O2_loc(i) < 5.0) then
                poc_diss_loc = POC_diss * parm_lowo2_remin_factor
             endif

     !-----------------------------------------------------------------------
     !  apply scalelength factor to length scales
     !-----------------------------------------------------------------------

# ifdef TDEP_REMIN
             poc_diss_loc = scalelength(i) * poc_diss_loc / Tfunc_soft(i)
             sio2_diss_loc = scalelength(i) * sio2_diss_loc / Tfunc_soft(i)
             caco3_diss_loc = scalelength(i) * caco3_diss_loc / Tfunc_soft(i)
             dust_diss_loc = scalelength(i) * dust_diss_loc / Tfunc_soft(i)
# else
             poc_diss_loc = scalelength(i) * poc_diss_loc
             sio2_diss_loc = scalelength(i) * sio2_diss_loc
             caco3_diss_loc = scalelength(i) * caco3_diss_loc
             dust_diss_loc = scalelength(i) * dust_diss_loc
# endif

     !-----------------------------------------------------------------------
     !  apply temperature dependence to sio2_diss length scale
     !-----------------------------------------------------------------------
             ! MF: Not sure whether this should still be done? In CESM
             ! BEC, they dont do this anymore, but only apply the
             ! scaling factors, which in the end also represent some
             ! kind of temperature dependence!
             sio2_diss_loc = sio2_diss_loc / TfuncS(i)
     !-----------------------------------------------------------------------
     !  decay_POC_E and decay_SiO2 set locally, modified by O2
     !-----------------------------------------------------------------------

             decay_POC_E = exp(-dz_loc / poc_diss_loc)
             decay_SiO2  = exp(-dz_loc / sio2_diss_loc)
             decay_CaCO3 = exp(-dz_loc / caco3_diss_loc)
             decay_dust  = exp(-dz_loc / dust_diss_loc)

     !-----------------------------------------------------------------------
     !  Set outgoing fluxes for non-iron pools.
     !  The outoing fluxes for ballast materials are from the
     !  solution of the coresponding continuous ODE across the model
     !  level. The ODE has a constant source term and linear decay.
     !  It is assumed that there is no sub-surface dust production.
     !-----------------------------------------------------------------------

             P_CaCO3_sflux_out(i,j) = P_CaCO3_sflux_in(i,j) * decay_CaCO3 +
     &           P_CaCO3_prod(i) * ((c1 - P_CaCO3_gamma) * (c1 - decay_CaCO3)
     &             * caco3_diss_loc)

             P_CaCO3_hflux_out(i,j) = P_CaCO3_hflux_in(i,j) * DECAY_Hard(i) +
     &           P_CaCO3_prod(i) * (P_CaCO3_gamma * dz_loc)

             P_SiO2_sflux_out(i,j) = P_SiO2_sflux_in(i,j) * decay_SiO2 +
     &           P_SiO2_prod(i) * ((c1 - P_SiO2_gamma) * (c1 - decay_SiO2)
     &             * sio2_diss_loc)

             P_SiO2_hflux_out(i,j) = P_SiO2_hflux_in(i,j) * DECAY_Hard(i) +
     &           P_SiO2_prod(i) * (P_SiO2_gamma * dz_loc)

             dust_sflux_out(i,j) = dust_sflux_in(i,j) * decay_dust

             dust_hflux_out(i,j) = dust_hflux_in(i,j) * DECAY_HardDust(i)

     !-----------------------------------------------------------------------
     !  Compute how much POC_PROD is available for deficit reduction
     !  and excess POC flux after subtracting off fraction of non-dust
     !  ballast production from net POC_PROD.
     !-----------------------------------------------------------------------

             POC_PROD_avail = POC_prod(i) -
     &          P_CaCO3_rho * P_CaCO3_prod(i) -
     &          P_SiO2_rho * P_SiO2_prod(i)

     !-----------------------------------------------------------------------
     !  Check for POC production bounds violations
     !-----------------------------------------------------------------------

             if (POC_PROD_avail < c0) then
                print *, 'subroutine compute_particulate_terms:',
     &            ' non_dust ballast production exceeds POC production'
                print*,'k',k
                print *, 'POC_prod_avail: ', POC_prod_avail
                print*,'POC_prod(i)',POC_prod(i)
                print*,'P_CaCO3_rho*P_CaCO3_prod(i)',P_CaCO3_rho * P_CaCO3_prod(i)
                print*,'P_SiO2_rho * P_SiO2_prod(i)',P_SiO2_rho * P_SiO2_prod(i)
                poc_error = .true.
             endif

     !-----------------------------------------------------------------------
     !  Compute 1st approximation to new QA_dust_def, the QA_dust
     !  deficit leaving the cell (implicit). In the case of explicit
     !  sinking, new_QA_dust_def equals dust flux leaving the layer.
     !  Ignore POC_PROD_avail at this stage.
     !-----------------------------------------------------------------------

             if (QA_dust_def(i) > 0) then
                new_QA_dust_def = QA_dust_def(i) *
     &             (dust_sflux_out(i,j) + dust_hflux_out(i,j)) /
     &             (dust_sflux_in(i,j) + dust_hflux_in(i,j))
             else
                new_QA_dust_def = c0
             endif

     !-----------------------------------------------------------------------
     !  Use POC_PROD_avail to reduce new_QA_dust_def.
     !-----------------------------------------------------------------------

             if (new_QA_dust_def > c0) then
                new_QA_dust_def = new_QA_dust_def - POC_PROD_avail * dz_loc
                if (new_QA_dust_def < c0) then
                   POC_PROD_avail = -new_QA_dust_def * dzr_loc
                   new_QA_dust_def = c0
                else
                   POC_PROD_avail = c0
                endif
             endif

             QA_dust_def(i) = new_QA_dust_def

     !-----------------------------------------------------------------------
     !  Compute outgoing POC fluxes. QA POC flux is computing using
     !  ballast fluxes and new_QA_dust_def. If no QA POC flux came in
     !  and no production occured, then no QA POC flux goes out. This
     !  shortcut is present to avoid roundoff cancellation errors from
     !  the dust_rho * dust_flux_out - QA_dust_def computation.
     !  Any POC_PROD_avail still remaining goes into excess POC flux.
     !-----------------------------------------------------------------------

             if (POC_hflux_in(i,j) == c0 .and. POC_prod(i) == c0) then
                POC_hflux_out(i,j) = c0
             else
                POC_hflux_out(i,j) = P_CaCO3_rho *
     &             (P_CaCO3_sflux_out(i,j) + P_CaCO3_hflux_out(i,j)) +
     &             P_SiO2_rho *
     &             (P_SiO2_sflux_out(i,j) + P_SiO2_hflux_out(i,j)) +
     &             dust_rho *
     &             (dust_sflux_out(i,j) + dust_hflux_out(i,j)) -
     &             new_QA_dust_def
                POC_hflux_out(i,j) = max(POC_hflux_out(i,j), c0)
             endif

             POC_sflux_out(i,j) = POC_sflux_in(i,j) * decay_POC_E +
     &          POC_PROD_avail *((c1 - decay_POC_E) *
     &          poc_diss)

     !-----------------------------------------------------------------------
     !  Compute remineralization terms. It is assumed that there is no
     !  sub-surface dust production.
     !-----------------------------------------------------------------------

             P_CaCO3_remin(i,j) = P_CaCO3_prod(i) +
     &          ((P_CaCO3_sflux_in(i,j) - P_CaCO3_sflux_out(i,j)) +
     &          (P_CaCO3_hflux_in(i,j) - P_CaCO3_hflux_out(i,j))) * dzr_loc

             P_SiO2_remin(i,j) = P_SiO2_prod(i) +
     &          ((P_SiO2_sflux_in(i,j) - P_SiO2_sflux_out(i,j)) +
     &          (P_SiO2_hflux_in(i,j) - P_SiO2_hflux_out(i,j))) * dzr_loc

             POC_remin(i,j) = POC_prod(i) +
     &          ((POC_sflux_in(i,j) - POC_sflux_out(i,j)) +
     &          (POC_hflux_in(i,j) - POC_hflux_out(i,j))) * dzr_loc

             dust_remin(i) =
     &          ((dust_sflux_in(i,j) - dust_sflux_out(i,j)) +
     &          (dust_hflux_in(i,j) - dust_hflux_out(i,j))) * dzr_loc

     !-----------------------------------------------------------------------
     !  Compute iron remineralization and flux out.
     !-----------------------------------------------------------------------

             if ((POC_sflux_in(i,j)+POC_hflux_in(i,j)) == c0) then
                P_Iron_remin(i,j) = (POC_remin(i,j) * parm_Red_Fe_C)
             else
                P_Iron_remin(i,j) = (POC_remin(i,j) *
     &             (P_iron_sflux_in(i,j) + P_iron_hflux_in(i,j)) /
     &             (POC_sflux_in(i,j) + POC_hflux_in(i,j)))
             endif
             P_Iron_remin(i,j) = P_Iron_remin(i,j)
     &             + P_iron_sflux_in(i,j) * 1.5e-5
             ! MF: What does the hard-coded 1.5e-5 stand for? Applied to
             ! a flux of mmol Fe m-2 s-1? CESM BEC says Fe desorption
             ! from particulate pool. If this is the case, shouldnt it
             ! be treate like the other mineral pools with a dissolution
             ! length and an associated decay?

             P_iron_sflux_out(i,j) = P_iron_sflux_in(i,j) + dz_loc *
     &           ((c1 - P_iron_gamma) * P_iron_prod(i) - P_Iron_remin(i,j))

             if (P_iron_sflux_out(i,j) < c0) then
                P_iron_sflux_out(i,j) = c0
                P_Iron_remin(i,j) = P_iron_sflux_in(i,j) * dzr_loc +
     &            (c1 - P_iron_gamma) * P_iron_prod(i)
             endif

     !-----------------------------------------------------------------------
     !  Compute iron release from dust remin/dissolution
     !
     !  dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
     !                      gFe     molFe     nmolFe
     !  Also add in Fe source from sediments if applicable to this cell.
     !-----------------------------------------------------------------------


             P_Iron_remin(i,j) = P_Iron_remin(i,j)
     &           + dust_remin(i) * dust_to_Fe

!        if (wrt_bec2_diag_3d(pironremin_idx_t)) then
!           idiag = idx_bec2_diag_3d(pironremin_idx_t)
!           bec2_diag_3d(istr:iend,j,k,idiag) = P_Iron_remin(istr:iend,j)*dz_loc
!        endif

             P_iron_hflux_out(i,j) = P_iron_hflux_in(i,j)

          else  ! LANDMASK(i,j)
             P_CaCO3_sflux_out(i,j) = c0
             P_CaCO3_hflux_out(i,j) = c0
             P_CaCO3_remin(i,j) = c0

             P_SiO2_sflux_out(i,j) = c0
             P_SiO2_hflux_out(i,j) = c0
             P_SiO2_remin(i,j) = c0

             dust_sflux_out(i,j) = c0
             dust_hflux_out(i,j) = c0
             dust_remin(i) = c0

             POC_sflux_out(i,j) = c0
             POC_hflux_out(i,j) = c0
             POC_remin(i,j) = c0

             P_iron_sflux_out(i,j) = c0
             P_iron_hflux_out(i,j) = c0
             P_Iron_remin(i,j) = c0
          endif  ! LANDMASK(i,j)

     !-----------------------------------------------------------------------
     !  Bottom Sediments Cell?
     !  If so compute sedimentary burial and denitrification N losses.
     !  Using empirical relations from Bohlen et al., 2012 (doi:10.1029/2011GB004198) for Sed Denitrification
     !  OTHER_REMIN estimates organic matter remineralized in the sediments
     !      by the processes other than oxic remin and denitrification (SO4 and CO2,
     !      etc..)
     !      based on Soetaert et al., 1996, varies between 10% and 50%
     !      0.4_r8 is a coefficient with units mmolC/cm2/yr sinking flux,
     !      OTHER_REMIN is 50% above this high flux value,
     !      In special case where bottom O2 has been depleted to < 1.0 uM,
     !               all sedimentary remin is due to DENITRIFICATION + OTHER_REMIN
     !  POC burial from Dunne et al. 2007 (doi:10.1029/2006GB002907), maximum of 80% burial efficiency imposed
     !  Bsi preservation in sediments based on
     !     Ragueneau et al. 2000 (doi:10.1016/S0921-8181(00)00052-7)
     !  Calcite is preserved in sediments above the lysocline, dissolves below.
     !       Here a constant depth is used for lysocline.
     !-----------------------------------------------------------------------

          if (LANDMASK(i,j) .and. (k == 1)) then

             flux = POC_sflux_out(i,j)+POC_hflux_out(i,j) ! mmol C/m^2/s

             if (flux > c0) then
               flux_alt = flux*spd ! convert to mmol C/m^2/day

               POC_sed_loss(i,j) = flux * min(0.8, parm_POMbury
     &             * (0.013 + 0.53 * flux_alt*flux_alt / (7.0 + flux_alt)**2))

            ! set denitr to zero if NO3<=0 in bottom water (to not
            ! further deplete already not existing NO3)
            ! if this works: actually, this might still take out too
            ! much NO3 if NO3_BW is just slightly above zero. Include
            ! something like "if SED_DENITRIF > NO3_loc,
            ! reduce SED_DENITRIF to available NO3"
            ! check: is some NO3 threshold implemented for water col
            ! denitrif? could also cause neg NO3...
            ! FOR NOW: NO3min=5 (see Simon's draft!! he uses 2, but said
            ! it's kind of arbitrary...), below this
            ! sed_denitr is reduced! same for water column denitrf!

            ! constrain SED_DENITRIF by what is left to
            ! be remineralized after POC_sed_loss
               if (NO3_loc(i) < parm_sed_denitrif_NO3_limit) then
               SED_DENITRIF(i) = dzr_loc * min(flux * (0.06 + 0.19 *
     &                 0.99**(O2_loc(i)-NO3_loc(i))),
     &                flux-POC_sed_loss(i,j)*denitrif_C_N)
               SED_DENITRIF(i) = 1/parm_sed_denitrif_NO3_limit*NO3_loc(i)*SED_DENITRIF(i)
               else
               SED_DENITRIF(i) = dzr_loc * min(flux * (0.06 + 0.19 *
     &                 0.99**(O2_loc(i)-NO3_loc(i))),
     &                flux-POC_sed_loss(i,j)*denitrif_C_N)
               endif

#ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(seddenitrif_idx_t)) then
          idiag = idx_bec2_diag_2d(seddenitrif_idx_t)
          bec2_diag_2d(i,j,idiag) = SED_DENITRIF(i)*dz_loc
       endif
       if (wrt_bec2_diag_2d(fluxtosed_idx_t)) then
          idiag = idx_bec2_diag_2d(fluxtosed_idx_t)
          bec2_diag_2d(i,j,idiag) = flux
       endif
       if (wrt_bec2_diag_2d(pocsedloss_idx_t)) then
          idiag = idx_bec2_diag_2d(pocsedloss_idx_t)
          bec2_diag_2d(i,j,idiag) = POC_sed_loss(i,j)
       endif
#endif

               flux_alt = flux*1e-4*spd*365.0 ! convert to mmol C/cm^2/year
               OTHER_REMIN(i) = dzr_loc
     &             * min(min(0.1 + flux_alt,0.5) * (flux - POC_sed_loss(i,j)),
     &                   (flux - POC_sed_loss(i,j) - (SED_DENITRIF(i)*dz_loc*denitrif_C_N)))

     !----------------------------------------------------------------------------------
     !              if bottom water O2 is depleted, assume all remin is denitrif + other
     !----------------------------------------------------------------------------------

               if (O2_loc(i) < c1) then
                 OTHER_REMIN(i) = dzr_loc *
     &                (flux - POC_sed_loss(i,j) - (SED_DENITRIF(i)*dz_loc*denitrif_C_N))
              endif

#ifdef BEC2_DIAG
       if (wrt_bec2_diag_2d(otherremin_idx_t)) then
          idiag = idx_bec2_diag_2d(otherremin_idx_t)
          bec2_diag_2d(i,j,idiag) = OTHER_REMIN(i)*dz_loc
       endif
#endif

            ! FESEDFLUX parametrization from Dale 2015
            ! carbon oxidation rate in mmol m-2 s-1
            flux_oxidated(i) = flux - POC_sed_loss(i,j)

            ! convert from mmol m-2 s-1 to mmol m-2 d-1
            flux_oxidated(i) = flux_oxidated(i) * 86400  ! day2sec=86400

            ! constant 170 in umol m-2 d-1, flux_oxidated in mmol m-2
            ! day-1, o2 in muM
            if (flux_oxidated(i) > 0) then
                FESEDFLUX(i) = 170.0 * tanh(flux_oxidated(i)/O2_loc(i))
            endif
            ! convert from umol m-2 d-1 to mmol m-2 s-1
!            if (latr(i,j) > 5.0 .and. latr(i,j) < 20.0) then
!                FESEDFLUX(i) = FESEDFLUX(i)*0.001/86400*1.0     !sec2day=1/86400
!            else
                FESEDFLUX(i) = FESEDFLUX(i)*0.001/86400*1.0
!            endif
# ifdef BEC2_DIAG
        if (wrt_bec2_diag_2d(fesedflux_idx_t)) then
           idiag = idx_bec2_diag_2d(fesedflux_idx_t)
           bec2_diag_2d(i,j,idiag) = FESEDFLUX(i)
        endif
# endif
            endif  ! flux > c0

            flux = P_SiO2_sflux_out(i,j)+P_SiO2_hflux_out(i,j)
            flux_alt = flux*spd ! convert to mmol/m^2/day
            ! first compute burial efficiency, then compute loss to sediments
            if (flux_alt > c1) then
               P_SiO2_sed_loss(i,j) = 0.3 * flux_alt - 0.06
            else
               P_SiO2_sed_loss(i,j) = 0.04
            endif
            P_SiO2_sed_loss(i,j) = flux * parm_BSIbury * P_SiO2_sed_loss(i,j)

# ifdef BEC2_DIAG
        if (wrt_bec2_diag_2d(sio2fluxtosed_idx_t)) then
           idiag = idx_bec2_diag_2d(sio2fluxtosed_idx_t)
           bec2_diag_2d(i,j,idiag) = flux
        endif
        if (wrt_bec2_diag_2d(sio2sedloss_idx_t)) then
           idiag = idx_bec2_diag_2d(sio2sedloss_idx_t)
           bec2_diag_2d(i,j,idiag) = P_SiO2_sed_loss(i,j)
        endif
# endif

            flux = c0  ! set flux back to zero
            if (-z_w(i,j,k) < 3300.0) then
               flux = P_CaCO3_sflux_out(i,j) + P_CaCO3_hflux_out(i,j)
               P_CaCO3_sed_loss(i,j) = flux
            endif

# ifdef BEC2_DIAG
        if (wrt_bec2_diag_2d(caco3fluxtosed_idx_t)) then
           idiag = idx_bec2_diag_2d(caco3fluxtosed_idx_t)
           bec2_diag_2d(i,j,idiag) = flux
        endif
        if (wrt_bec2_diag_2d(caco3sedloss_idx_t)) then
           idiag = idx_bec2_diag_2d(caco3sedloss_idx_t)
           bec2_diag_2d(i,j,idiag) = P_CaCO3_sed_loss(i,j)
        endif
# endif

     !----------------------------------------------------------------------------------
     !  Update sinking fluxes and remin fluxes, accounting for sediments.
     !  flux used to hold sinking fluxes before update.
     !----------------------------------------------------------------------------------

            flux = P_CaCO3_sflux_out(i,j) + P_CaCO3_hflux_out(i,j)
            if (flux > c0) then
               P_CaCO3_remin(i,j) = P_CaCO3_remin(i,j)
     &             + ((flux - P_CaCO3_sed_loss(i,j)) * dzr_loc)
            endif

            flux = P_SiO2_sflux_out(i,j) + P_SiO2_hflux_out(i,j)
            if (flux > c0) then
               P_SiO2_remin(i,j) = P_SiO2_remin(i,j)
     &             + ((flux - P_SiO2_sed_loss(i,j)) * dzr_loc)
            endif

            flux = POC_sflux_out(i,j) + POC_hflux_out(i,j)
            if (flux > c0) then
               POC_remin(i,j) = POC_remin(i,j)
     &             + ((flux - POC_sed_loss(i,j)) * dzr_loc)
            endif

     !-----------------------------------------------------------------------
     !   Remove all Piron and dust that hits bottom, sedimentary Fe source
     !        accounted for by FESEDFLUX elsewhere.
     !-----------------------------------------------------------------------


            flux = c0  ! set flux back to zero
            flux = (P_iron_sflux_out(i,j) + P_iron_hflux_out(i,j))
            if (flux > c0) then
               P_iron_sed_loss(i,j) = flux
            endif

            dust_sed_loss(i,j) = dust_sflux_out(i,j) + dust_hflux_out(i,j)

# ifdef BEC2_DIAG
        if (wrt_bec2_diag_2d(pironfluxtosed_idx_t)) then
           idiag = idx_bec2_diag_2d(pironfluxtosed_idx_t)
           bec2_diag_2d(i,j,idiag) = flux
        endif
        if (wrt_bec2_diag_2d(dustfluxtosed_idx_t)) then
           idiag = idx_bec2_diag_2d(dustfluxtosed_idx_t)
           bec2_diag_2d(i,j,idiag) = dust_sed_loss(i,j)
        endif
# endif

     !-----------------------------------------------------------------------
     !   Bottom layer: set all outgoing fluxes to 0.0
     !-----------------------------------------------------------------------

            if (k == 1) then
               P_CaCO3_sflux_out(i,j) = c0
               P_CaCO3_hflux_out(i,j) = c0

               P_SiO2_sflux_out(i,j) = c0
               P_SiO2_hflux_out(i,j) = c0

               dust_sflux_out(i,j) = c0
               dust_hflux_out(i,j) = c0

               POC_sflux_out(i,j) = c0
               POC_hflux_out(i,j) = c0

               P_iron_sflux_out(i,j) = c0
               P_iron_hflux_out(i,j) = c0
            endif

             P_Iron_remin(i,j) = P_Iron_remin(i,j)
     &           + (FESEDFLUX(i) * dzr_loc)
         end if  ! LANDMASK(i,j) .and. (k == 1)

       end do  ! i = istr, iend

# ifdef BEC2_DIAG
        if (wrt_bec2_diag_3d(pironremin_idx_t)) then
           idiag = idx_bec2_diag_3d(pironremin_idx_t)
           bec2_diag_3d(istr:iend,j,k,idiag) = P_Iron_remin(istr:iend,j)*dz_loc
        endif
        if (wrt_bec2_diag_3d(caco3remin_idx_t)) then
           idiag = idx_bec2_diag_3d(caco3remin_idx_t)
           bec2_diag_3d(istr:iend,j,k,idiag) = P_CaCO3_remin(istr:iend,j)*dz_loc
        endif
        if (wrt_bec2_diag_3d(sio2remin_idx_t)) then
           idiag = idx_bec2_diag_3d(sio2remin_idx_t)
           bec2_diag_3d(istr:iend,j,k,idiag) = P_SiO2_remin(istr:iend,j)*dz_loc
        endif
        if (wrt_bec2_diag_3d(pocremin_idx_t)) then
           idiag = idx_bec2_diag_3d(pocremin_idx_t)
           bec2_diag_3d(istr:iend,j,k,idiag) = POC_remin(istr:iend,j)*dz_loc
        endif
        if (wrt_bec2_diag_3d(dustremin_idx_t)) then
           idiag = idx_bec2_diag_3d(dustremin_idx_t)
           bec2_diag_3d(istr:iend,j,k,idiag) = dust_remin(istr:iend)*dz_loc
        endif
# endif

       return
       end subroutine compute_particulate_terms

# ifndef BULK_FRC
!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE WS            ****************************
!
!-------------------------------------------------------------------------

       subroutine WS(SMFTX, SMFTY,landmaskc,work,istr,iend,jstr,jend)
    !-------------------------------------------------------------------------
    !   result declaration
    !-------------------------------------------------------------------------

! needed for scalars
        use param
        use scalars

        implicit none

! needed for rho0
        integer::istr,iend,jstr,jend
        REAL, DIMENSION(istr:iend,jstr:jend), INTENT(OUT)  :: WORK ! 10m windspeed (m/s)

     !-------------------------------------------------------------------------
     !   Compute 10m windspeed from surface velocity fluxes where LAND_MASK is
     !   true. Give zero where LAND_MASK is false.
     !-------------------------------------------------------------------------

     !-------------------------------------------------------------------------
     !   Derivation
     !
     !   stress (N/m^2) is a function of surface velocity fluxes
     !   (taux,tauy) = (SMFTX,SMFTY) / momentum_factor
     !
     !   turbulent velocity scale (m/s) is a function of stress
     !   ustar^2 = sqrt(taux^2 + tauy^2) / rho_air
     !           = sqrt(SMFTX^2 + SMFTY^2) / (momentum_factor rho_air)
     !
     !   turbulent velocity scale (m/s) is a function of 10m windspeed (m/s)
     !   ustar^2 = coef_1 ws + coef_2 WS^2 + coef_3 WS^3
     !   ustar^2 = c_d WS^2
     !
     !   Compute ustar^2 from velocity fluxes and solve cubic equation
     !   for 10m windspeed. Initial guess for iterative solver comes
     !   from assuming a constant drag coefficient c_d.
     !
     !-----------------------------------------------------------------------
        real,parameter::
     &    rho_air   = 1.2            ! ambient air density (kg/m^3)

     !-----------------------------------------------------------------------
     !   arguments
     !-------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) ::
     &    SMFTX,  ! 'zonal' surface velocity flux (m^2/s^2)
     &    SMFTY    ! 'meridional' surface velocity flux (m^2/s^2)

     !------------------------------------------------------------------------
     !   local variables
     !-------------------------------------------------------------------------

        REAL, PARAMETER ::
     &    coef_1  = 0.0027,     ! (m/s)
     &    coef_2  = 0.000142,   ! (non-dimensional)
     &    coef_3  = 0.0000764,  ! (s/m)
     &    c_d     = 1.7e-3      ! drag coefficient for initial u

        REAL, DIMENSION(istr:iend,jstr:jend) ::
     &    ustar_squared  ! square of turbulent velocity scale (m^2/s^2)
        LOGICAL,DIMENSION(istr:iend,jstr:jend),intent(in) :: landmaskc
     !------------------------------------------------------------------------
     !   solve for WS in (m/s)
     !-------------------------------------------------------------------------

       WHERE (landmaskc)
         ustar_squared = SQRT(SMFTX**2 + SMFTY**2) * rho0 /
     &         rho_air

         WORK = SQRT(ustar_squared / c_d)

       !----------------------------------------------------------------------
       !   perform 3 Newton iterations
       !   loop constructs are not allowed inside WHERE constructs
       !----------------------------------------------------------------------

          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3))
     &         - ustar_squared) /
     &         (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))

          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3))
     &          - ustar_squared) /
     &         (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))

          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3))
     &           - ustar_squared) /
     &       (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))

        ELSEWHERE
           WORK = 0.0
        END WHERE

        RETURN
        END subroutine WS
# endif

# if defined Ncycle_SY
!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE CSCHMIDT_N2O   ****************************
!
!-------------------------------------------------------------------------

         subroutine CSCHMIDT_N2O(SSTT,landmaskc,SCHMIDT_N2O,istr,
     &          iend,jstr,jend)
        implicit none

       integer,intent(in)::istr,iend,jstr,jend
       logical,dimension(istr:iend,jstr:jend),intent(in) :: landmaskc  ! DevinD slicing the array makes a local copy.
                                                             ! should just send through full array
     !---------------------------------------------------------------------------
     !   Compute Schmidt number of N2O in seawater as function of SST
     !   where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   ref : Sarmiento and Gruber 2006 book (OBCD, table 3.3.1 page 85
     !---------------------------------------------------------------------------
     !
     ! over land, do NOT assign zero to the result, but use very large value instead
     ! in order to prevent division by zero.

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

           REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

           REAL, DIMENSION(istr:iend,jstr:jend),intent(out) :: SCHMIDT_N2O

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

          REAL, PARAMETER ::
     &      a = 2301.1,
     &      b = 151.1,
     &      c = 4.7364,
     &      d = 0.059431

         WHERE (LANDMASKC)
           SCHMIDT_N2O = a + SSTT* (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
            SCHMIDT_N2O = 1e30
          END WHERE
          return
          END subroutine CSCHMIDT_N2O

          subroutine N2OSATU(SSTT, SSSS,landmaskc, N2OSAT,istr,
     &         iend,jstr,jend)

        use bec2_params

        implicit none
! needed for c0
!#include "bgc_param_bec2.h"
       integer,intent(in)::istr,iend,jstr,jend
       LOGICAL,DIMENSION(istr:iend,jstr:jend),intent(in) :: landmaskc

     !---------------------------------------------------------------------------
     !
     !   Computes n2o saturation concentration at 1 atm n2o pressure
     !   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
     !   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   FROM WEISS AND PRICE (1980), Marine chemistry.
     !
     !---------------------------------------------------------------------------

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) ::
     &    SSTT,  ! sea surface temperature (C)
     &    SSSS    ! sea surface salinity (psu)

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

         REAL, DIMENSION(istr:iend,jstr:jend),intent(out) :: N2OSAT

     !---------------------------------------------------------------------------
     !   local variables
     !---------------------------------------------------------------------------

         REAL, DIMENSION(istr:iend,jstr:jend) :: TS, S_LOC

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

        REAL, PARAMETER ::
     &    a_1 = -165.8802,
     &    a_2 = 222.8743,
     &    a_3 = 92.0792,
     &    a_4 = -1.48425,
     &    b_1 = -0.056235,
     &    b_2 = 0.031619,
     &    b_3 = -0.0048472

        S_LOC = max(1.00e-4,SSSS)

        WHERE (LANDMASKC)
           TS = T0_Kelvin + SSTT

         N2OSAT =  exp(a_1 + a_2 * (100.0/TS) + a_3 * log(TS/100.0) +
     &          a_4 * (TS/100.0)**2 + S_LOC * (b_1 +b_2 * (TS/100.0) + b_3 *
     &          (TS/100.0)**2))
        ELSEWHERE
          N2OSAT = c0
       END WHERE
     !---------------------------------------------------------------------------
     !   Convert from mol/l to mmol/m^3
     !---------------------------------------------------------------------------

        N2OSAT = N2OSAT * 1000 * 1000

        return
        END subroutine N2OSATU

     !---------------------------------------------------------------------------
     !---------------------------------------------------------------------------
     !   CSCHMIDT_N2
     !---------------------------------------------------------------------------
     !---------------------------------------------------------------------------


         subroutine CSCHMIDT_N2(SSTT,landmaskc,SCHMIDT_N2,istr,
     &          iend,jstr,jend)
        implicit none

        integer,intent(in)::istr,iend,jstr,jend
        LOGICAL,DIMENSION(istr:iend,jstr:jend),intent(in) :: landmaskc

     !---------------------------------------------------------------------------
     !   Compute Schmidt number of N2 in seawater as function of SST
     !   where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   ref : Sarmiento and Gruber 2006 (Table 3.3.1)
     !
     !---------------------------------------------------------------------------
     !
     ! over land, do NOT assign zero to the result, but use very large value instead
     ! in order to prevent division by zero.

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend),intent(out) :: SCHMIDT_N2

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

          REAL, PARAMETER ::
     &      a = 2206.1,
     &      b = 144.86,
     &      c = 4.5413,
     &      d = 0.056988

         WHERE (LANDMASKC)
           SCHMIDT_N2 = a + SSTT* (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
            SCHMIDT_N2 = 1e30
          END WHERE
          return
          END subroutine CSCHMIDT_N2

     !*****************************************************************************


          subroutine N2SATU(SSTT, SSSS, landmaskc, N2SAT,istr,
     &         iend,jstr,jend)

        use bec2_params

        implicit none
       integer,intent(in)::istr,iend,jstr,jend
       LOGICAL,DIMENSION(istr:iend,jstr:jend),intent(in) :: landmaskc

     !---------------------------------------------------------------------------
     !
     !   Computes dinitrogen saturation concentration at 1 atm total pressure
     !   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
     !   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   FROM Sarmiento and Gruber 2006
     !   --> Solubility : table 3.2.2
     !   --> apply correction for moist air following Panel 3.2.1
     !
     !   N2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
     !   Simon Yang, Feb 2018
     !---------------------------------------------------------------------------

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) ::
     &    SSTT,  ! sea surface temperature (C)
     &    SSSS    ! sea surface salinity (psu)

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend),intent(out) :: N2SAT

     !---------------------------------------------------------------------------
     !   local variables
     !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend) :: TS, BUNSEN, PH2O_to_P,
     &     S_LOC

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

        REAL, PARAMETER ::
     &    a_1 = -59.6274,
     &    a_2 = 85.7761,
     &    a_3 = 24.3696,
     &    b_1 = -0.051580,
     &    b_2 = 0.026329,
     &    b_3 = -0.0037252,
     &    vbarn2o = 22.4136 ! volar molume of N2 at standard temperature and pressure

        S_LOC = max(1.00e-4,SSSS)

        WHERE (LANDMASKC)
           TS = T0_Kelvin + SSTT

          BUNSEN=  exp(a_1 + a_2 * (100.0/TS) + a_3 * log(TS/100.0) +
     &          S_LOC * (b_1 +b_2 * (TS/100.0) + b_3 *
     &          (TS/100.0)**2))

        PH2O_to_P = exp(24.4543-67.4509*(100.0/TS)-4.8489*log(TS/100.0)
     &              -0.000544*S_LOC)

        N2SAT = ( BUNSEN / vbarn2o ) * (c1 - PH2O_to_P)

        ELSEWHERE
          N2SAT = c0
       END WHERE
     !---------------------------------------------------------------------------
     !   Convert from ml/l to mmol/m^3
     !---------------------------------------------------------------------------

        N2SAT = N2SAT * 1e6

        return
        END subroutine N2SATU

     !*****************************************************************************
# endif
     !---------------------------------------------------------------------------
     !   Subroutine CSCHMIDT_O2
     !---------------------------------------------------------------------------

         subroutine CSCHMIDT_O2(SSTT,landmaskc,SCHMIDT_O2,istr,
     &          iend,jstr,jend)
        implicit none

        integer,intent(in)::istr,iend,jstr,jend
        LOGICAL,DIMENSION(istr:iend,jstr:jend),intent(in) :: landmaskc

     !---------------------------------------------------------------------------
     !   Compute Schmidt number of O2 in seawater as function of SST
     !   where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   ref : Keeling et al, Global Biogeochem. Cycles, Vol. 12,
     !         No. 1, pp. 141-163, March 1998
     !---------------------------------------------------------------------------
     !
     ! over land, do NOT assign zero to the result, but use very large value instead
     ! in order to prevent division by zero.

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend),intent(out) :: SCHMIDT_O2

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

          REAL, PARAMETER ::
     &      a = 1638.0,
     &      b = 81.83,
     &      c = 1.483,
     &      d = 0.008004

         WHERE (LANDMASKC)
           SCHMIDT_O2 = a + SSTT* (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
            SCHMIDT_O2 = 1e30
          END WHERE
          return
          END subroutine CSCHMIDT_O2

     !*****************************************************************************

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE O2SATU   ****************************
!
!-------------------------------------------------------------------------

          subroutine O2SATU(SSTT, SSSS,landmaskc, O2SAT,istr,
     &         iend,jstr,jend)

        use bec2_params

        implicit none
        integer,intent(in)::istr,iend,jstr,jend
        LOGICAL,DIMENSION(istr:iend,jstr:jend),intent(in) :: landmaskc

     !---------------------------------------------------------------------------
     !
     !   Computes oxygen saturation concentration at 1 atm total pressure
     !   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
     !   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
     !   THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
     !
     !   *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
     !   *** IT SHOULDN'T BE THERE.                                ***
     !
     !   O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
     !   0 permil <= S <= 42 permil
     !   CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
     !   O2SAT = 282.015 mmol/m^3
     !
     !---------------------------------------------------------------------------

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) ::
     &    SSTT,  ! sea surface temperature (C)
     &    SSSS    ! sea surface salinity (psu)

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend),intent(out) :: O2SAT

     !---------------------------------------------------------------------------
     !   local variables
     !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend) :: TS, S_LOC

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

        REAL, PARAMETER ::
     &    a_0 = 2.00907,
     &    a_1 = 3.22014,
     &    a_2 = 4.05010,
     &    a_3 = 4.94457,
     &    a_4 = -2.56847E-1,
     &    a_5 = 3.88767,
     &    b_0 = -6.24523E-3,
     &    b_1 = -7.37614E-3,
     &    b_2 = -1.03410E-2,
     &    b_3 = -8.17083E-3,
     &    c_0 = -4.88682E-7

        S_LOC = max(1.00e-4,SSSS)

        WHERE (LANDMASKC)
           TS = LOG( ((T0_Kelvin+25.0) - SSTT) / (T0_Kelvin + SSTT) )

         O2SAT = EXP(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) +
     &         S_LOC*( (b_0+TS*(b_1+TS*(b_2+TS*b_3))) + S_LOC*c_0 ))
        ELSEWHERE
          O2SAT = c0
       END WHERE
     !---------------------------------------------------------------------------
     !   Convert from ml/l to mmol/m^3
     !---------------------------------------------------------------------------

        O2SAT = O2SAT * 44.6596

        return
        END subroutine O2SATU

     !*****************************************************************************

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE CSCHMIDT_CO2   ****************************
!
!-------------------------------------------------------------------------

          subroutine CSCHMIDT_CO2(SSTT,landmaskc,SCHMIDT_CO2,istr,
     &                     iend,jstr,jend )

        use bec2_params

        implicit none
        integer,intent(in)::istr,iend,jstr,jend
        LOGICAL,DIMENSION(istr:iend,jstr:jend),intent(in) :: landmaskc

     !---------------------------------------------------------------------------
     !   Compute Schmidt number of CO2 in seawater as function of SST
     !   where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
     !   pp. 7373-7382, May 15, 1992
     !---------------------------------------------------------------------------
     ! over land, do NOT assign zero to the result, but use very large value instead
     ! in order to prevent division by zero.

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(OUT) :: SCHMIDT_CO2

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

         REAL, PARAMETER ::
     &    a = 2073.1,
     &    b = 125.62,
     &    c = 3.6276,
     &    d = 0.043219

         WHERE (LANDMASKC)
           SCHMIDT_CO2 = a + SSTT * (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
           SCHMIDT_CO2 = 1e30
         END WHERE
         return
        END subroutine  CSCHMIDT_CO2

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE CO2CALC_ROW   ***************************
!-------------------------------------------------------------------------

       SUBROUTINE co2calc_row(mask, t, s, dic_in, ta_in, pt_in,
     &    sit_in, phlo, phhi, ph, xco2_in, atmpres, co2star,
     &    dco2star, pco2oc, dpco2, istr,iend)


     !---------------------------------------------------------------------------
     !   SUBROUTINE CO2CALC
     !
     !   PURPOSE : Calculate delta co2*, etc. from total alkalinity, total CO2,
     !             temp, salinity (s), etc.
     !---------------------------------------------------------------------------
        implicit none
       real,parameter::c0=0.0, c1=1.0, c10=10.0,
     &   c1000=1000.0, T0_Kelvin=273.16,rho_sw=4.1/3.996
       integer::istr,iend
     !---------------------------------------------------------------------------
     !   input arguments
     !---------------------------------------------------------------------------
       REAL, PARAMETER :: xacc = 1e-10
       INTEGER, PARAMETER :: maxit = 100

     !-----------------------------------------------------------------------------
     !   declarations for function coefficients & species concentrations
     !-----------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend) ::
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff,
     &  bt, st, ft, dic, ta, pt, sit

         LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask
         REAL, DIMENSION(istr:iend), INTENT(IN) ::
     &    t,         ! temperature (degrees C)
     &    s,         ! salinity (PSU)
     &    dic_in,    ! total inorganic carbon (mmol/m^3)
     &    ta_in,     ! total alkalinity (meq/m^3)
     &    pt_in,     ! inorganic phosphate (mmol/m^3)
     &    sit_in,    ! inorganic silicate (mmol/m^3)
     &    phlo,      ! lower limit of pH range
     &    phhi,      ! upper limit of pH range
     &    xco2_in,   ! atmospheric mole fraction CO2 in dry air (ppmv)
     &    atmpres     ! atmospheric pressure (atmosphere)


     !---------------------------------------------------------------------------
     !   output arguments
     !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend), INTENT(OUT) ::
     &    ph,        ! computed ph values, for initial guess on next time step
     &    co2star,   ! CO2*water (mmol/m^3)
     &    dco2star,  ! delta CO2 (mmol/m^3)
     &    pco2oc,    ! oceanic pCO2 (ppmv)
     &    dpco2      ! Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)


     !---------------------------------------------------------------------------
     !   local variable declarations
     !---------------------------------------------------------------------------

       INTEGER :: i

        REAL ::
     &    mass_to_vol,   ! (mol/kg) -> (mmol/m^3)
     &    vol_to_mass,   ! (mmol/m^3) -> (mol/kg)
     &    tk,            ! temperature (K)
     &    is,            ! ionic strength
     &    scl,           ! chlorinity
     &    co2starair,    ! co2star saturation
     &    tk100, tk1002, invtk, dlogtk, is2, sqrtis,
     &    s2, sqrts, s15, htotal2

        REAL, DIMENSION(istr:iend) ::
     &    xco2,          ! atmospheric CO2 (atm)
     &    htotal,        ! free concentration of H ion
     &    x1, x2,        ! bounds on htotal for solver
     &    s_loc          ! local copy of salinity with fixed negative values

     !---------------------------------------------------------------------------
     !   check for existence of ocean points
     !---------------------------------------------------------------------------

       s_loc = max(1.00e-4,s)

       IF (COUNT(mask) == 0) THEN
         ph          = c0
         co2star     = c0
         dco2star    = c0
         pCO2oc      = c0
         dpCO2       = c0
       RETURN
      END IF


     !---------------------------------------------------------------------------
     !   set unit conversion factors
     !---------------------------------------------------------------------------

       mass_to_vol = 1e6 * rho_sw
       vol_to_mass = c1 / mass_to_vol

     !---------------------------------------------------------------------------
     !   convert tracer units to per mass  xco2 from uatm to atm
     !---------------------------------------------------------------------------

        DO i = istr,iend
          IF (mask(i)) THEN
            dic(i)  = dic_in(i)  * vol_to_mass
            ta(i)   = ta_in(i)   * vol_to_mass
            pt(i)   = pt_in(i)   * vol_to_mass
            sit(i)  = sit_in(i)  * vol_to_mass
           xco2(i) = xco2_in(i) * 1e-6

          !---------------------------------------------------------------------
          !   Calculate all constants needed to convert between various
          !   measured carbon species. References for each equation are
          !   noted in the code.  Once calculated, the constants are stored
          !   and passed in the common block "const". The original version
          !   of this code was based on the code by Dickson in Version 2 of
          !   "Handbook of Methods for the Analysis of the Various Parameters
          !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
          !   p25-26).
          !   Derive simple terms used more than once
          !---------------------------------------------------------------------

            tk       = T0_Kelvin + t(i)
            tk100    = tk * 1e-2
            tk1002   = tk100 * tk100
            invtk    = c1 / tk
            dlogtk   = LOG(tk)

            is       = 19.924 * s_loc(i) / (c1000 - 1.005 * s_loc(i))
            is2      = is * is
            sqrtis   = SQRT(is)
            sqrts    = SQRT(s_loc(i))
            s15      = s_loc(i) ** 1.5
            s2       = s_loc(i) ** 2
            scl      = s_loc(i) / 1.80655

          !---------------------------------------------------------------------
          !   f = k0(1-pH2O)*correction term for non-ideality
          !   Weiss  Price (1980, Mar. Chem., 8, 347-359;
          !                 Eq 13 with table 6 values)
          !---------------------------------------------------------------------

            ff(i) = EXP(-162.8301 + 218.2968/tk100 +
     &          90.9241*LOG(tk100) -
     &          1.47696*tk1002 + s_loc(i)*(.025695 - .025225*tk100 +
     &          0.0049867*tk1002))

          !---------------------------------------------------------------------
          !   K0 from Weiss 1974
          !---------------------------------------------------------------------

            k0(i) = EXP(93.4517/tk100 - 60.2409 + 23.3585*LOG(tk100) +
     &          s_loc(i)*(0.023517 - 0.023656 * tk100 +
     &          0.0047036 * tk1002))

          !---------------------------------------------------------------------
          !   k1 = [H][HCO3]/[H2CO3]
          !   k2 = [H][CO3]/[HCO3]
          !   Millero p.664 (1995) using Mehrbach et al. data on seawater scale
          !---------------------------------------------------------------------

          k1(i) = 10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -
     &         0.0118*s_loc(i) + 0.000116*s2))

          k2(i) = 10**(-1*(1394.7*invtk + 4.777 -
     &           0.0184*s_loc(i) + 0.000118*s2))

          !---------------------------------------------------------------------
          !   kb = [H][BO2]/[HBO2]
          !   Millero p.669 (1995) using data from Dickson (1990)
          !---------------------------------------------------------------------

          kb(i) = EXP((-8966.90 - 2890.53*sqrts - 77.942*s_loc(i) +
     &          1.728*s15 - 0.0996*s2)*invtk +
     &          (148.0248 + 137.1942*sqrts + 1.62142*s_loc(i)) +
     &          (-24.4344 - 25.085*sqrts - 0.2474*s_loc(i)) *
     &          dlogtk + 0.053105*sqrts*tk)

          !---------------------------------------------------------------------
          !   k1p = [H][H2PO4]/[H3PO4]
          !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k1p(i) = EXP(-4576.752*invtk + 115.525 - 18.453 * dlogtk +
     &         (-106.736*invtk + 0.69171) * sqrts +
     &          (-0.65643*invtk - 0.01844) * s_loc(i))

          !---------------------------------------------------------------------
          !   k2p = [H][HPO4]/[H2PO4]
          !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
          !---------------------------------------------------------------------

          k2p(i) = EXP(-8814.715*invtk + 172.0883 - 27.927 * dlogtk +
     &          (-160.340*invtk + 1.3566) * sqrts +
     &          (0.37335*invtk - 0.05778) * s_loc(i))

          !---------------------------------------------------------------------
          !   k3p = [H][PO4]/[HPO4]
          !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k3p(i) = EXP(-3070.75*invtk - 18.141 +
     &          (17.27039*invtk + 2.81197) * sqrts +
     &          (-44.99486*invtk - 0.09984) * s_loc(i))

          !---------------------------------------------------------------------
          !   ksi = [H][SiO(OH)3]/[Si(OH)4]
          !   Millero p.671 (1995) using data from Yao and Millero (1995)
          !---------------------------------------------------------------------

          ksi(i) = EXP(-8904.2*invtk + 117.385 - 19.334 * dlogtk +
     &         (-458.79*invtk + 3.5913) * sqrtis +
     &         (188.74*invtk - 1.5998) * is +
     &          (-12.1652*invtk + 0.07871) * is2 +
     &          LOG(1.0-0.001005*s_loc(i)))

          !---------------------------------------------------------------------
          !   kw = [H][OH]
          !   Millero p.670 (1995) using composite data
          !---------------------------------------------------------------------

          kw(i) = EXP(-13847.26*invtk + 148.9652 - 23.6521 * dlogtk +
     &        (118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     &         sqrts - 0.01615 * s_loc(i))

          !---------------------------------------------------------------------
          !   ks = [H][SO4]/[HSO4]
          !   Dickson (1990, J. chem. Thermodynamics 22, 113)
          !---------------------------------------------------------------------

          ks(i) = EXP(-4276.1*invtk + 141.328 - 23.093*dlogtk +
     &          (-13856*invtk + 324.57 - 47.986*dlogtk) *
     &          sqrtis +
     &          (35474*invtk - 771.54 + 114.723*dlogtk) * is -
     &          2698*invtk*is**1.5 + 1776*invtk*is2 +
     &          LOG(1.0 - 0.001005*s_loc(i)))

          !---------------------------------------------------------------------
          !   kf = [H][F]/[HF]
          !   Dickson and Riley (1979) -- change pH scale to total
          !---------------------------------------------------------------------

          kf(i) = EXP(1590.2*invtk - 12.641 + 1.525*sqrtis +
     &         LOG(1.0 - 0.001005*s_loc(i)) +
     &          LOG(1.0 + (0.1400/96.062)*(scl)/ks(i)))

          !---------------------------------------------------------------------
          !   Calculate concentrations for borate, sulfate, and fluoride
          !   bt : Uppstrom (1974)
          !   st : Morris  Riley (1966)
          !   ft : Riley (1965)
          !---------------------------------------------------------------------

          bt(i) = 0.000232 * scl/10.811
          st(i) = 0.14 * scl/96.062
          ft(i) = 0.000067 * scl/18.9984

          x1(i) = c10 ** (-phhi(i))
          x2(i) = c10 ** (-phlo(i))

        END IF ! if mask

       END DO ! i loop


     !---------------------------------------------------------------------------
     !   If DIC and TA are known then either a root finding or iterative
     !   method must be used to calculate htotal. In this case we use
     !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
     !   (function "rtsafe.f" with error trapping removed).
     !
     !   As currently set, this procedure iterates about 12 times. The
     !   x1 and x2 values set below will accomodate ANY oceanographic
     !   values. If an initial guess of the pH is known, then the
     !   number of iterations can be reduced to about 5 by narrowing
     !   the gap between x1 and x2. It is recommended that the first
     !   few time steps be run with x1 and x2 set as below. After that,
     !   set x1 and x2 to the previous value of the pH +/- ~0.5.
     !---------------------------------------------------------------------------

        CALL drtsafe_row(mask, x1, x2, xacc, htotal,istr,
     &    iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
     !---------------------------------------------------------------------------
     !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
     !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
     !---------------------------------------------------------------------------

        DO i = istr,iend
          IF (mask(i)) THEN

          htotal2 = htotal(i) ** 2
          co2star(i) = dic(i) * htotal2 /
     &          (htotal2 + k1(i)*htotal(i) + k1(i)*k2(i))
          co2starair = xco2(i) * ff(i) * atmpres(i)
          dco2star(i) = co2starair - co2star(i)
          ph(i) = -LOG10(htotal(i))

          !---------------------------------------------------------------------
          !   Add two output arguments for storing pco2oc
          !   Should we be using K0 or ff for the solubility here?
          !---------------------------------------------------------------------

          pco2oc(i) = co2star(i) / ff(i)
          dpCO2(i)    = pco2oc(i) - xco2(i) * atmpres(i)

          !---------------------------------------------------------------------
          !   Convert units of output arguments
          !   Note: pco2oc and dpCO2 are calculated in atm above.
          !---------------------------------------------------------------------

          co2star(i)  = co2star(i) * mass_to_vol
          dco2star(i) = dco2star(i) * mass_to_vol

          pco2oc(i) = pco2oc(i) * 1e6
          dpCO2(i)    = dpCO2(i) * 1e6

        ELSE ! if mask

          ph(i)       = c0
          co2star(i)  = c0
          dco2star(i) = c0
          pco2oc(i) = c0
          dpCO2(i)    = c0


        END IF ! if mask
       END DO ! i loop

      END SUBROUTINE co2calc_row

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE TALK_ROW   ****************************
!
!-------------------------------------------------------------------------

       SUBROUTINE talk_row(mask, x, fn, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)

        implicit none
       REAL, PARAMETER :: xacc = 1e-10
       INTEGER, PARAMETER :: maxit = 100

     !-----------------------------------------------------------------------------
     !   declarations for function coefficients & species concentrations
     !-----------------------------------------------------------------------------
         integer :: istr,iend

        REAL, DIMENSION(istr:iend) ::
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff,
     &  bt, st, ft, dic, ta, pt, sit

     !---------------------------------------------------------------------------
     !   This routine computes TA as a function of DIC, htotal and constants.
     !   It also calculates the derivative of this function with respect to
     !   htotal. It is used in the iterative solution for htotal. In the call
     !   "x" is the input value for htotal, "fn" is the calculated value for
     !   TA and "df" is the value for dTA/dhtotal.
     !---------------------------------------------------------------------------
         real,parameter::c1=1.0, c2=2.0, c3=3.0

     !---------------------------------------------------------------------------
     !   input arguments
     !---------------------------------------------------------------------------
        LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask
        REAL, DIMENSION(istr:iend), INTENT(IN) :: x

     !---------------------------------------------------------------------------
     !   output arguments
     !---------------------------------------------------------------------------

         REAL, DIMENSION(istr:iend), INTENT(OUT) :: fn, df

     !---------------------------------------------------------------------------
     !   local variable declarations
     !---------------------------------------------------------------------------

         INTEGER :: i

          REAL ::
     &    x1, x2, x3, k12, k12p, k123p, a, a2, da, b, b2, db, c

     !---------------------------------------------------------------------------

      DO i = istr,iend
        IF (mask(i)) THEN
          x1 = x(i)
          x2 = x1 * x1
          x3 = x2 * x1
          k12 = k1(i) * k2(i)
          k12p = k1p(i) * k2p(i)
          k123p = k12p * k3p(i)
          a = x3 + k1p(i) * x2 + k12p * x1 + k123p
          a2 = a * a
          da = c3 * x2 + c2 * k1p(i) * x1 + k12p
          b = x2 + k1(i) * x1 + k12
          b2 = b * b
          db = c2 * x1 + k1(i)
          c = c1 + st(i)/ks(i)
          !---------------------------------------------------------------------
          !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
          !---------------------------------------------------------------------

          fn(i) = k1(i) * x1 * dic(i)/b +
     &          c2 * dic(i) * k12/b +
     &          bt(i)/(c1 + x1/kb(i)) +
     &          kw(i)/x1 +
     &          pt(i) * k12p * x1/a +
     &          c2 * pt(i) * k123p/a +
     &          sit(i)/(c1 + x1/ksi(i)) -
     &          x1/c -
     &          st(i)/(c1 + ks(i)/x1/c) -
     &          ft(i)/(c1 + kf(i)/x1) -
     &          pt(i) * x3/a -
     &          ta(i)

          !---------------------------------------------------------------------
          !   df = d(fn)/dx
          !---------------------------------------------------------------------

          df(i) = ((k1(i)*dic(i)*b) - k1(i)*x1*dic(i)*db)/b2 -
     &          c2 * dic(i) * k12 * db/b2 -
     &          bt(i)/kb(i)/(c1+x1/kb(i)) ** 2 -
     &          kw(i)/x2 +
     &          (pt(i) * k12p * (a - x1 * da))/a2 -
     &          c2 * pt(i) * k123p * da/a2 -
     &          sit(i)/ksi(i)/(c1+x1/ksi(i)) ** 2 -
     &          c1/c +
     &          st(i) * (c1 + ks(i)/x1/c)**(-2) * (ks(i)/c/x2) +
     &          ft(i) * (c1 + kf(i)/x1)**(-2) * kf(i)/x2 -
     &          pt(i) * x2 * (c3 * a - x1 * da)/a2

        END IF ! if mask
        END DO ! i loop

      END SUBROUTINE talk_row

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE DRTSAFE_ROW   ****************************
!
!-------------------------------------------------------------------------

        SUBROUTINE drtsafe_row(mask_in, x1, x2, xacc, soln,istr,
     &   iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
        implicit none
       INTEGER, PARAMETER :: maxit = 100
       integer::istr,iend
     !-----------------------------------------------------------------------------
     !   declarations for function coefficients & species concentrations
     !-----------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend) ::
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff,
     &  bt, st, ft, dic, ta, pt, sit

     !---------------------------------------------------------------------------
     !   Vectorized version of drtsafe, which was a modified version of
     !      Numerical Recipes algorithm.
     !   Keith Lindsay, Oct 1999
     !
     !   Algorithm comment :
     !      Iteration from Newton s method is used unless it leaves
     !      bracketing interval or the dx is > 0.5 the previous dx.
     !      In that case, bisection method is used.
     !---------------------------------------------------------------------------
        real,parameter::c0=0.0, c2=2.0

     !---------------------------------------------------------------------------
     !   input arguments
     !---------------------------------------------------------------------------

        LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask_in
        REAL, DIMENSION(istr:iend), INTENT(IN) :: x1, x2
        REAL, INTENT(IN) :: xacc

     !---------------------------------------------------------------------------
     !   output arguments
     !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend), INTENT(OUT) :: soln

     !---------------------------------------------------------------------------
     !   local variable declarations
     !---------------------------------------------------------------------------

       LOGICAL :: leave_bracket, dx_decrease
       LOGICAL, DIMENSION(istr:iend) :: mask
       INTEGER ::  i, it
       REAL :: temp
       REAL, DIMENSION(istr:iend) :: xlo, xhi, flo, fhi, f,
     &     df, dxold, dx

     !---------------------------------------------------------------------------
     !   bracket root at each location and set up first iteration
     !---------------------------------------------------------------------------

       mask = mask_in

       CALL talk_row(mask, x1, flo, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
       CALL talk_row(mask, x2, fhi, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)

       DO i = istr,iend
         IF (mask(i)) THEN
           IF (flo(i) .LT. c0) THEN
              xlo(i) = x1(i)
              xhi(i) = x2(i)
           ELSE
             xlo(i) = x2(i)
             xhi(i) = x1(i)
             temp = flo(i)
             flo(i) = fhi(i)
             fhi(i) = temp
          END IF
          soln(i) = 0.5 * (xlo(i) + xhi(i))
          dxold(i) = ABS(xlo(i) - xhi(i))
          dx(i) = dxold(i)
         END IF
        END DO

        CALL talk_row(mask, soln, f, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)

     !---------------------------------------------------------------------------
     !   perform iterations, zeroing mask when a location has converged
     !---------------------------------------------------------------------------

       DO it = 1,maxit
        DO i = istr,iend
           IF (mask(i)) THEN
             leave_bracket = ((soln(i)-xhi(i))*df(i)-f(i)) *
     &             ((soln(i)-xlo(i))*df(i)-f(i)) .GE. 0
             dx_decrease = ABS(c2 * f(i)) .LE. ABS(dxold(i) * df(i))
             IF (leave_bracket .OR. .NOT. dx_decrease) THEN
                dxold(i) = dx(i)
                dx(i) = 0.5 * (xhi(i) - xlo(i))
                soln(i) = xlo(i) + dx(i)
                IF (xlo(i) .EQ. soln(i)) mask(i) = .FALSE.
             ELSE
                dxold(i) = dx(i)
                dx(i) = -f(i) / df(i)
                temp = soln(i)
                soln(i) = soln(i) + dx(i)
                IF (temp .EQ. soln(i)) mask(i) = .FALSE.
             END IF
             IF (ABS(dx(i)) .LT. xacc) mask(i) = .FALSE.
          END IF
        END DO

        IF (.NOT. ANY(mask)) RETURN

        CALL talk_row(mask, soln, f, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)

        DO i = istr,iend
          IF (mask(i)) THEN
             IF (f(i) .LT. c0) THEN
                xlo(i) = soln(i)
                flo(i) = f(i)
             ELSE
                xhi(i) = soln(i)
                fhi(i) = f(i)
             END IF
          END IF
        END DO

       END DO ! iteration loop

       END SUBROUTINE drtsafe_row


#else  /* BIOLOGY_BEC2 */
       subroutine ecosys_bec2_tile_empty ()

       end
#endif  /* BIOLOGY_BEC2 */
