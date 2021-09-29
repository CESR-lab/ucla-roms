#include "cppdefs.opt"

      subroutine step2d(tile)

! Advance barotropic mode variables (2D momenta and free surface) by
! one barotropic time step using Generalized Forward-Backward AB3-AM4
! stepping algorithm [Sec. 2.3 from SM2005 starting with Eq. (2.49)]
! and perform fast-time averaging to interact with baroclinic mode.

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none
      integer tile

#include "compute_tile_bounds.h"
      call step2d_FB_tile( istr,iend,jstr,jend, A2d(1,1),  A2d(1,2),
     &                    A2d(1, 3), A2d(1, 4), A2d(1, 5), A2d(1, 6),
     &                    A2d(1, 7), A2d(1, 8), A2d(1, 9), A2d(1,10),
     &                               A2d(1,11), A2d(1,12), A2d(1,13))
      end

      subroutine step2d_FB_tile( istr,iend,jstr,jend, zeta_new,Dnew,
     &                           rubar,rvbar, urhs,vrhs,  DUon,DVom,
     &                                       Drhs, UFx,UFe,VFx,VFe)

      use param
      use river_frc
      use pipe_frc
      use wec_frc
#if defined POT_TIDES && !defined SOLVE3D
      use tides, only: Ptide
#endif
      use climat
      use coupling
      use grid
      use mixing
      use ocean2d
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j, kbak, kold
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: zeta_new, Dnew,
     &                         rubar,rvbar,  urhs,vrhs,  DUon,DVom,
     &                                       Drhs, UFx,UFe,VFx,VFe
      real fwd,fwd1,fwd2, bkw_new,bkw,bkw1,bkw2,
     &            cff,cff1,cff2, DUnew,DVnew

#include "compute_auxiliary_bounds.h"

! Set coefficients for AB3-AM4 forward-backward algorithm.
!---- ------------ --- ------- ------- -------- ----------
! Because Forward Euler step is used to update "zeta" during the first
! barotropic step, the pressure-gradient term in momentum equation must
! be computed via Backward step to keep it numerically stable. However,
! this interferes with the computation of forcing terms "rufrc" and
! "rvfrc" because the free surface in pressure gradient computation in 3D
! is exactly at time corresponding to baroclinic step "nstp" (rather
! than ahead by one barotropic step after it is updated by a normal
! forward-backward step).  To resolve this conflict, the pressure
! gradient term is computed in two  stages during the first barotropic
! step: use zeta(:,:,kstp) at first to ensure exact consistency with
! 3D mode; then, after vertical integrals of 3D r.h.ss "rufrc" and
! "rvfrc" are converted into forcing terms, add correction based on
! the difference zeta_new(:,:)-zeta(:,:,kstp) to "rubar" and "rvbar"
! to make them consistent with Backward step for pressure gradient.
! gradient terms (search for label PGF_FB_CORRECTION below).

      if (FIRST_2D_STEP) then            ! Meaning of time indices
        kbak=kstp ; kold=kstp            ! ------- -- ---- -------
        fwd=1. ; fwd1=0. ; fwd2=0.       ! m-2   m-1   m     m+1
#ifdef SOLVE3D
        bkw_new=0. ; bkw=1.              ! kold  kbak  kstp  knew
#else
        bkw_new=1. ; bkw=0.              ! fwd2  fwd1  fwd
#endif
        bkw1=0.    ; bkw2=0.             ! bkw2  bkw1  bkw   bkw_new
      elseif (FIRST_2D_STEP+1) then
        kbak=kstp-1
        if (kbak < 1) kbak=4
        kold=kbak


        fwd=1.; fwd1=0. ; fwd2=0.        ! Logically AB2-AM3 forward-
        bkw_new=1.0833333333333          ! backward scheme with maximum
        bkw=   -0.1666666666666          ! stability coefficients while
        bkw1=   0.0833333333333          ! maintaining third-order
        bkw2=0.                          ! accuracy, alpha_max=1.73
      else
        kbak=kstp-1
        if (kbak < 1) kbak=4
        kold=kbak-1
        if (kold < 1) kold=4


         fwd=1.781105  ; fwd1=-1.06221 ; fwd2=0.281105
         bkw_new=0.614 ; bkw=0.285 ; bkw1=0.088 ; bkw2=0.013
      endif

! Preliminary step: compute total depth of water column and vertically
! ----------- ----- integrated fluxes needed for computing horizontal
! divergence to advance free surface and for nonlinear advection terms
! for the barotropic momentum equations. DUon == H*dy*Ubar, DVom ==
! H*dx*Vbar, at u, and v positions respectively
! In case of WEC, ust2d is added to Ubar (and vst2d added to Vbar)

      do j=jstrV-2,jend+1
        do i=istrU-2,iend+1
          Drhs(i,j)=h(i,j) +fwd*zeta(i,j,kstp) +fwd1*zeta(i,j,kbak)
     &                                          +fwd2*zeta(i,j,kold)
        enddo
      enddo
      do j=jstr-1,jend+1
        do i=istrU-1,iend+1
          urhs(i,j)=fwd*ubar(i,j,kstp) +fwd1*ubar(i,j,kbak)
     &                                         +fwd2*ubar(i,j,kold)
#  if defined WEC && defined MASKING
          urhs(i,j)=urhs(i,j)*umask(i,j)+ust2d(i,j)*(umask(i,j)-1.0)
#  endif
          DUon(i,j)=0.5*(Drhs(i,j)+Drhs(i-1,j))*dn_u(i,j)*( urhs(i,j)
#  ifdef WEC
     &                                                  + ust2d(i,j)
#  endif
     &                                                               )
        enddo
      enddo
      do j=jstrV-1,jend+1
        do i=istr-1,iend+1
          vrhs(i,j)=fwd*vbar(i,j,kstp) +fwd1*vbar(i,j,kbak)
     &                                         +fwd2*vbar(i,j,kold)
#  if defined WEC && defined MASKING
          vrhs(i,j)=vrhs(i,j)*vmask(i,j)+vst2d(i,j)*(vmask(i,j)-1.0)
#  endif
          DVom(i,j)=0.5*(Drhs(i,j)+Drhs(i,j-1))*dm_v(i,j)*(vrhs(i,j)
#  ifdef WEC
     &                                                  + vst2d(i,j)
#  endif
     &                                                               )
        enddo
      enddo


! Advance free-surface:  Compute "zeta_new" at new time step and
!-------- ---- --------  interpolate it half-step backward, "zwrk"
! for the subsequent computation of barotropic pressure gradient.

#define zwrk UFx
#define rzeta  UFe
#define rzeta2  VFe
#define rzetaSA VFx
      do j=jstrV-1,jend !! why are these from ...-1?
        do i=istrU-1,iend
          zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &            *(DUon(i,j)-DUon(i+1,j)+DVom(i,j)-DVom(i,j+1))
#ifdef PIPE_SOURCE
          zeta_new(i,j)=zeta_new(i,j) + dtfast*pm(i,j)*pn(i,j)*pipe_flx(i,j)
#endif
#ifdef MASKING
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j) ! not sure this is needed
#endif
          Dnew(i,j)=zeta_new(i,j)+h(i,j)

          zwrk(i,j)=bkw_new*zeta_new(i,j) +bkw*zeta(i,j,kstp)
     &              +bkw1*zeta(i,j,kbak) +bkw2*zeta(i,j,kold)

#if defined VAR_RHO_2D && defined SOLVE3D
          rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
#else
          rzeta(i,j)=zwrk(i,j)
#endif
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j) !! zeta^2
        enddo
      enddo


! Set boundary conditions for the free-surface, after which copy
! the newly computed free-surface zeta_new into shared array.

      call zetabc_tile (istr,iend,jstr,jend, zeta_new)

      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,knew)=zeta_new(i,j)
        enddo
      enddo

! Compute fast-time-averaged fields over all short time steps.
!-------- ---- ---- -------- ------ ---- --- ----- ----------
! Reset/initialize arrays for averaged fields during the first
! barotropic time step; Accumulate averages after that. Include
! physical boundary points, but not periodic ghost points or
! computation  MPI computational margins.

#ifdef SOLVE3D
        cff1=weight(1,iif) ; cff2=weight(2,iif)

        if (FIRST_2D_STEP) then
          do j=jstrR,jendR
            do i=istrR,iendR
# ifdef EXTRAP_BAR_FLUXES
! PAC23 verified setting with theoretical alpha_max=1.0877. 
! the number 0.1024390... is the ratio of delta/gamma=0.21/2.05
! this should match the simlarly-labelled setting of NOW MID BAK
! in set_depth.F.   All the other settings are to be considered
! as experimental.

              DU_avg_bak(i,j)=DU_avg1(i,j)
     &                       -0.1024390243902439D0*DU_avg2(i,j)
              DV_avg_bak(i,j)=DV_avg1(i,j)
     &                       -0.1024390243902439D0*DV_avg2(i,j)


# endif
              Zt_avg1(i,j)=cff1*zeta(i,j,knew)
              DU_avg1(i,j)=0. ; DV_avg1(i,j)=0.
              DU_avg2(i,j)=cff2*DUon(i,j)
              DV_avg2(i,j)=cff2*DVom(i,j)
            enddo
          enddo
        else
          do j=jstrR,jendR
            do i=istrR,iendR
              Zt_avg1(i,j)=Zt_avg1(i,j) + cff1*zeta(i,j,knew)
              DU_avg2(i,j)=DU_avg2(i,j) + cff2*DUon(i,j)
              DV_avg2(i,j)=DV_avg2(i,j) + cff2*DVom(i,j)
            enddo
          enddo
        endif
#endif

! Compute pressure-gradient terms:  NOTE that "rubar" and "rvbar" are
!-------- -------- -------- ------  computed within the same fused loop
! despite the fact that their normal index ranges are different. Fusing
! the loops causes redundant computation of one column of "rubar" on
! the western physical boundary and one row of "rvbar" on the southern,
! but, at the same time it allows to share references to array elements
! (i,j) which results in an increase of computational density by almost
! a factor of 1.5 resulting in overall more efficient code.

      cff=0.5*g
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=cff*dn_u(i,j)*( (h(i-1,j)+h(i,j))*(rzeta(i-1,j)
     &                        -rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)

#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(zwrk(i-1,j)-zwrk(i,j)) )
#endif

#if defined TIDES && !defined SOLVE3D
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                   *( Ptide(i,j)-Ptide(i-1,j) )
#endif

#ifdef WEC
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i-1,j) )
#endif
     &                                                              )
!>
          rvbar(i,j)=cff*dm_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)

#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(zwrk(i,j-1)-zwrk(i,j)) )
#endif

#if defined TIDES && !defined SOLVE3D
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                   *( Ptide(i,j)-Ptide(i,j-1) )
#endif

#ifdef WEC
     &                  + ( h(i,j-1)+h(i,j)+rzeta(i,j-1)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i,j-1) )
#endif
     &                                                              )
        enddo
      enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA

#undef rzetaSA
#undef rzeta2
#undef rzeta
#undef zwrk

! Suppress computation of momentum advection, Coriolis, and lateral
! viscosity terms in the case of 3D model configuration. Because all
! these terms are already included in the baroclinic-to-barotropic
! forcing, "rufrc" and "rvfrc", undefining some or all of them just
! below does not mean entirely omitting them, but is rather choice
! between recomputing them at every barotropic step, or keeping
! them "frozen" in fast time.

#ifdef SOLVE3D
# ifdef CORR_COUPLED_MODE
#  ifdef EXTRAP_BAR_FLUXES
#   undef UV_ADV
#  endif
# else
#  undef UV_ADV
# endif
# ifndef KEEP_CORIOLIS
#  undef UV_COR
# endif
# undef UV_VIS2
#endif

! Compute horizontal advection terms for momentum equations (2D only)
!-------- ---------- --------- ----- --- -------- --------- --- -----
! NOTE: mathematically necessary (minimal) index ranges for momentum-
! flux components are
!
!      UFx(istrU-1:iend,jstr:jend)   VFx(istr:iend+1,jstrV:jend)
!      UFe(istrU:iend,jstr:jend+1)   VFe(istr,iend,jstrV-1,jend)
!
! however, for the purpose computational efficiency, these ranges are
! unified by suppressing U,V-suffices in order to allow fusion of the
! consecutive loops. This leads to slight increase of the redundant
! computations near western and southern boundaries in non-periodic
! directions.

#ifdef UV_ADV
      do j=jstr,jend
        do i=istr-1,iend
          UFx(i,j)=0.25*(DUon(i,j)+DUon(i+1,j))
     &                     *(urhs(i,j)+urhs(i+1,j))

          VFx(i+1,j)=0.25*(DUon(i+1,j)+DUon(i+1,j-1))
     &                       *(vrhs(i+1,j)+vrhs(i,j))
# ifdef MASKING
     &                                 *pmask(i+1,j)
# endif
        enddo
      enddo
      do j=jstr-1,jend
        do i=istr,iend
          VFe(i,j)=0.25*(DVom(i,j)+DVom(i,j+1))
     &                      *(vrhs(i,j)+vrhs(i,j+1))

          UFe(i,j+1)=0.25*(DVom(i,j+1)+DVom(i-1,j+1))
     &                       *(urhs(i,j+1)+urhs(i,j))
# ifdef MASKING
     &                                 *pmask(i,j+1)
# endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)

          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo !--> discard UFx,VFe,UFe,VFx, DUon,DVom
#endif


! Compute Coriolis terms and contribution to to the advection terms
! associated with the curvilinear nature of horizontal coordinates.


#if defined UV_COR || (defined CURVGRID && defined UV_ADV)
      do j=jstrV-1,jend
        do i=istrU-1,iend
          cff=Drhs(i,j)*(
# ifdef UV_COR
     &                             fomn(i,j)
# endif
# if (defined CURVGRID && defined UV_ADV)
     &  +0.5*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &        -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
# endif
     &                                          )

#  ifdef WEC
#    if defined CURVGRID && defined UV_ADV
          cff1 = Drhs(i,j)*(
     &    0.5*( dndx(i,j)*(vst2d(i,j)+vst2d(i,j+1))
     &          -dmde(i,j)*(ust2d(i,j)+ust2d(i+1,j)) ))
#    else
          cff1 = 0.0
#    endif
          UFx(i,j)=(cff+cff1)*(vrhs(i,j)+vrhs(i,j+1))
     &                 +cff*(vst2d(i,j)+vst2d(i,j+1))
          VFe(i,j)=(cff+cff1)*(urhs(i,j)+urhs(i+1,j))
     &                 +cff*(ust2d(i,j)+ust2d(i+1,j))
#  else
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
#  endif
        enddo
      enddo
      do j=jstr,jend
        do i=istrU,iend
          rubar(i,j)=rubar(i,j)+0.25*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          rvbar(i,j)=rvbar(i,j)-0.25*(VFe(i,j)+VFe(i,j-1))
        enddo
      enddo
#endif

! Compute horizontal viscous stress terms (2D only).

#ifdef UV_VIS2
      do j=jstr-1,jend
        do i=istr-1,iend
          cff=2.*Drhs(i,j)*visc2_r(i,j)
          UFx(i,j)=cff*(ubar(i+1,j,kstp)-ubar(i,j,kstp))*pm(i,j)
     &                                                *dn_r(i,j)
          VFe(i,j)=cff*(vbar(i,j+1,kstp)-vbar(i,j,kstp))*pn(i,j)
     &                                                *dm_r(i,j)


          cff1=0.0625*visc2_p(i+1,j+1)*(  Drhs(i,j)+Drhs(i+1,j)
     &                                +Drhs(i,j+1)+Drhs(i+1,j+1)

     &             )*( (pn(i+1,j+1)+pn(i,j+1)+pn(i+1,j)+pn(i,j))
     &                    *(ubar(i+1,j+1,kstp)-ubar(i+1,j,kstp))

     &                +(pm(i+1,j+1)+pm(i,j+1)+pm(i+1,j)+pm(i,j))
     &                    *(vbar(i+1,j+1,kstp)-vbar(i,j+1,kstp))
     &                                                        )
# ifdef MASKING
     &                     *pmask(i+1,j+1)
# endif
          UFe(i+1,j+1)=cff1*dm_p(i+1,j+1)
          VFx(i+1,j+1)=cff1*dn_p(i+1,j+1)
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=rubar(i,j)+UFx(i,j)-UFx(i-1,j)
     &                         +UFe(i,j+1)-UFe(i,j)

          rvbar(i,j)=rvbar(i,j)+VFx(i+1,j)-VFx(i,j)
     &                         +VFe(i,j)-VFe(i,j-1)
        enddo
      enddo
#endif /* UV_VIS2 */

! Bottom drag: linear and/or quadratic law, 2D mode only.

#ifndef SOLVE3D
# ifndef LINEAR_DRAG_ONLY
#  define rd VFx
      if (rdrg2 > 0.) then
        do j=jstrV-1,jend
          do i=istrU-1,iend
            rd(i,j)=rdrg + rdrg2*sqrt( 0.333333333333*(
     &              ubar(i,j,kstp)**2 +ubar(i+1,j,kstp)**2
     &                    +ubar(i,j,kstp)*ubar(i+1,j,kstp)
     &              +vbar(i,j,kstp)**2+vbar(i,j+1,kstp)**2
     &                    +vbar(i,j,kstp)*vbar(i,j+1,kstp)
     &                                                 ))
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend
            rubar(i,j)=rubar(i,j) - 0.5*(rd(i,j)+rd(i-1,j))
     &                  *dm_u(i,j)*dn_u(i,j)*ubar(i,j,kstp)

            rvbar(i,j)=rvbar(i,j) - 0.5*(rd(i,j)+rd(i,j-1))
     &                  *dm_v(i,j)*dn_v(i,j)*vbar(i,j,kstp)
          enddo
        enddo
      elseif (rdrg > 0.) then
# else
      if (rdrg > 0.) then
# endif
        do j=jstr,jend
          do i=istr,iend
            rubar(i,j)=rubar(i,j) -rdrg*dm_u(i,j)*dn_u(i,j)
     &                                      *ubar(i,j,kstp)

            rvbar(i,j)=rvbar(i,j) -rdrg*dm_v(i,j)*dn_v(i,j)
     &                                      *vbar(i,j,kstp)
          enddo
        enddo
      endif
# undef rd
#endif

! Add 2DH vortex-force terms combined with advection terms
#if (defined WEC && defined UV_ADV)
      do j=jstr,jend
        do i=istrU,iend
          vstu = 0.25*( vst2d(i,j) + vst2d(i,j+1)
     &                 +vst2d(i-1,j)+vst2d(i-1,j+1) )
          dudx = 0.5*( urhs(i+1,j)-urhs(i-1,j) )
          dvdx = 0.5*( vrhs(i,j) - vrhs(i-1,j)
     &                 +vrhs(i,j+1) - vrhs(i-1,j+1) )
          rubar(i,j) = rubar(i,j) + 0.5*dn_u(i,j)*
     &                      ( Drhs(i-1,j)+Drhs(i,j) )
     &               *( ust2d(i,j)*dudx + vstu*dvdx )
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          ustv = 0.25*( ust2d(i,j) + ust2d(i+1,j)
     &                 +ust2d(i,j-1)+ust2d(i+1,j-1) )
          dude = 0.5*( urhs(i,j) - urhs(i,j-1)
     &                 +urhs(i+1,j) - urhs(i+1,j-1) )
          dvde = 0.5*( vrhs(i,j-1)-vrhs(i,j+1) )
          rvbar(i,j) = rvbar(i,j) + 0.5*dm_v(i,j)*
     &                      ( Drhs(i,j-1)+Drhs(i,j) )
     &               *( ustv*dude + vst2d(i,j)*dvde )
        enddo
      enddo
#endif




! Coupling between 2D and 3D parts.
!--------- ------- -- --- -- ------
! Before the first barotropic time step arrays "rufrc" and "rvfrc"
! contain vertically integrals of the 3D r.h.s. terms for the momentum
! equations (including surface and bottom stresses).  During the first
! barotropic time step convert them into forcing terms by subtracting
! the fast-time "rubar" and "rvbar" from them;

! In the case of predictor-coupled mode the resultant forcing terms
! "rufrc" and "rvfrc" are extrapolated forward in time, so they become
! centered effectively at time n+1/2.   This is done using optimized
! Adams-Bashforth weights.  In the code below rufrc_bak(:,:,nstp) is
! at (n-1)time step, while rufrc_bak(:,:,3-nstp) is at (n-2).  After
! its use as input, the latter is overwritten by the value at time
! step "nstp" (mathematically "n") during the next step.

! From now on, the computed forcing terms "rufrc" and "rvfrc" will
! remain constant during  the fast-time stepping and will added to
! "rubar" and "rvbar" during all subsequent barotropic steps.

#ifdef SOLVE3D
      if (FIRST_2D_STEP) then
# ifdef CORR_COUPLED_MODE
#  if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
        kbak=kstp-1
        if (kbak < 1) kbak=4              ! Compensate offset in time
        do j=jstrV-1,jend                 ! between 3D and 2D modes in
          do i=istrU-1,iend               ! computing Coriolis terms.
            cff=Drhs(i,j)*fomn(i,j)
            UFx(i,j)=cff*(vbar(i,j,kbak)+vbar(i,j+1,kbak)
     &                             -vrhs(i,j)-vrhs(i,j+1))
            VFe(i,j)=cff*(ubar(i,j,kbak)+ubar(i+1,j,kbak)
     &                             -urhs(i,j)-urhs(i+1,j))
          enddo
        enddo
#  endif
# else  /*--> PRED_COUPLED_MODE */
        if (FIRST_TIME_STEP) then           ! Set coefficients for
          fwd=1.  ; fwd1=0. ;  fwd2=0.      ! AB3-like forward-in-time
        elseif (FIRST_TIME_STEP+1) then     ! extrapolation of 3D -> 2D
          fwd=1.5 ; fwd1=-0.5 ; fwd2=0.     ! forcing terms "rufrc" and
        else                                ! and "rvfrc" in the case
          fwd2=0.281105                     ! of predictor-coupled
          fwd1=-0.5-2.*fwd2                 ! barotropic mode.
          fwd=1.5+fwd2
        endif
# endif
        do j=jstr,jend
          do i=istr,iend

! Compensate for (cancel out) bottom drag terms: at input into
! step2d "rufrc" and "rvfrc" contain bottom drag terms computed by
! 3D mode, however there are no 2D couterparts in "rubar" and "rvbar"
! because 2D bottom drag will be computed implicitly during the final
! stage of updating ubar(:,:,knew) and vbar(:,:,knew) below.  Note
! that unlike other terms, bottom drag should not be extrapolated
! forward, if "rufrc" is, so this cancelation is done now rather
! than at the bottom of this loop.

            rufrc(i,j)=rufrc(i,j) + 0.5*(r_D(i,j)+r_D(i-1,j))
     &                    *dm_u(i,j)*dn_u(i,j)*ubar(i,j,kstp) !! not sure about this

            rvfrc(i,j)=rvfrc(i,j) + 0.5*(r_D(i,j)+r_D(i,j-1))
     &                    *dm_v(i,j)*dn_v(i,j)*vbar(i,j,kstp) !! not sure about this

!! Here, we remove all the terms that are in 2d rhs
            rufrc(i,j)=rufrc(i,j)-rubar(i,j)  ! Barotropic mode runnig
#  if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
     &           -0.25*(UFx(i,j)+UFx(i-1,j))  ! during corrector stage:
#  endif
            rvfrc(i,j)=rvfrc(i,j)-rvbar(i,j)  ! there is no need for
#  if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
     &           +0.25*(VFe(i,j)+VFe(i,j-1))  ! forward extrapolation.
#  endif

          enddo
        enddo

! Add correction term to shift pressure-gradient terms from "kstp"
! to "knew": in essence, convert the fist 2D step from forward-Euler
! to forward- backward (this is PGF_FB_CORRECTION mentioned above).


# define zwrk UFx
# define rzeta  UFe
# define rzeta2  VFe
# define rzetaSA VFx

        do j=jstrV-1,jend
          do i=istrU-1,iend
            zwrk(i,j)=zeta_new(i,j)-zeta(i,j,kstp)
# if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
# endif
          enddo
        enddo

        cff=0.5*g
        do j=jstr,jend
          do i=istr,iend
            rubar(i,j)=rubar(i,j) +cff*dn_u(i,j)*( (h(i-1,j)+h(i,j))
     &          *(rzeta(i-1,j)-rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)

# if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(zwrk(i-1,j)-zwrk(i,j)) )
# endif
     &                                                              )
!>
            rvbar(i,j)=rvbar(i,j) +cff*dm_v(i,j)*( (h(i,j-1)+h(i,j))
     &          *(rzeta(i,j-1)-rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)

# if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(zwrk(i,j-1)-zwrk(i,j)) )
# endif
     &                                                              )
          enddo
        enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA

# undef rzetaSA
# undef rzeta2
# undef rzeta
# undef zwrk
      endif   !<-- FIRST_2D_STEP
#endif

! Advance 2D momentum components while simultaneously adding them to
! accumulate fast-time-averages to compute barotropic fluxes. Doing so
! "on the fly" yields a more computationally dense code and eliminates
! repeated multiplication by Dnew (since the fluxes are available first
! as volatile variables DUnew, DVnew at this moment). However the fast-
! timme averaged fluxes are needed not only in the interior but at the
! physical boundaries as well, so this approach requires as set of
! loops along the boundaries just below the b.c. call.

#define Dstp DUon

      do j=jstrV-1,jend
        do i=istrU-1,iend
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo

      cff=0.5*dtfast
c<<     cff=2.*dtfast
#ifdef SOLVE3D
      cff1=0.5*weight(1,iif)
#else
      cff2=2.*dtfast
#endif
      do j=jstr,jend
        do i=istrU,iend
          DUnew=( (Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
c<<     &     +cff*iA_u(i,j)
#ifdef SOLVE3D
     &                            *(rubar(i,j)+rufrc(i,j))
#else
     &                            *rubar(i,j)
     &                                    +cff2*sustr(i,j)

# ifdef WEC 
     &                                    +cff2*brk2dx(i,j)
#  ifdef WAVE_FRICTION
     &                                    +cff2*frc2dx(i,j)
#  endif
# endif
#endif
     &                                                   )
#ifdef MASKING
     &                                         *umask(i,j)
#endif
          ubar(i,j,knew)=DUnew/( Dnew(i,j)+Dnew(i-1,j)
#ifdef SOLVE3D
     &                       +dtfast*(r_D(i,j)+r_D(i-1,j))
     &                                                   )

          DU_avg1(i,j)=DU_avg1(i,j)+cff1*dn_u(i,j)*( DUnew
# ifdef WEC
     &                 +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
# endif
     &                                                   )

#else
     &                                                   )
#endif
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          DVnew=( (Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
c<<     &     +cff*iA_v(i,j)
#ifdef SOLVE3D
     &                            *(rvbar(i,j)+rvfrc(i,j))
#else
     &                            *rvbar(i,j)
# ifdef WIND_AT_RHO_POINTS
     &                   +dtfast*(svstr(i,j-1)+svstr(i,j))
# else
     &                                    +cff2*svstr(i,j)
# endif
# ifdef WEC
     &                                    +cff2*brk2de(i,j)
#  ifdef WAVE_FRICTION
     &                                    +cff2*frc2de(i,j)
#  endif
# endif
#endif
     &                                                   )
#ifdef MASKING
     &                                         *vmask(i,j)
#endif
          vbar(i,j,knew)=DVnew/( Dnew(i,j)+Dnew(i,j-1)
#ifdef SOLVE3D
     &                       +dtfast*(r_D(i,j)+r_D(i,j-1))
     &                                                   )
 
          DV_avg1(i,j)=DV_avg1(i,j) +cff1*dm_v(i,j) * ( DVnew
# ifdef WEC 
     &                 +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
# endif
     &                                                   )

#else
     &                                                   )
#endif
        enddo
      enddo

! Set boundary conditions and compute integral mass flux across all
! open boundaries, if any.

      call    u2dbc_tile(istr,iend,jstr,jend, UFx)
      call    v2dbc_tile(istr,iend,jstr,jend, UFx)

! Fast-time-averaged barotropic fluxes along physical boundaries.

#ifdef SOLVE3D
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr-1,jendR
          Dnew(istr-1,j)=h(istr-1,j)+zeta_new(istr-1,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr-1,jendR
          Dnew(iend+1,j)=h(iend+1,j)+zeta_new(iend+1,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr-1,iendR
          Dnew(i,jstr-1)=h(i,jstr-1)+zeta_new(i,jstr-1)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr-1,iendR
          Dnew(i,jend+1)=h(i,jend+1)+zeta_new(i,jend+1)
        enddo
      endif
# endif
      cff1=0.5*weight(1,iif)
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstrR,jendR
          DU_avg1(istrU-1,j)=DU_avg1(istrU-1,j)+cff1*(Dnew(istrU-1,j)

     &         +Dnew(istrU-2,j))*( ubar(istrU-1,j,knew)
# ifdef WEC
     &                                             +ust2d(istrU-1,j)
# endif
     &                                             )*dn_u(istrU-1,j)  
        enddo
        do j=jstrV,jend
          DV_avg1(istr-1,j)=DV_avg1(istr-1,j) +cff1*(Dnew(istr-1,j)

     &       +Dnew(istr-1,j-1) )*( vbar(istr-1,j,knew)
# ifdef WEC
     &                                              +vst2d(istr-1,j)
# endif
     &                                              )*dm_v(istr-1,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstrR,jendR
          DU_avg1(iend+1,j)=DU_avg1(iend+1,j) +cff1*( Dnew(iend+1,j)
     &            +Dnew(iend,j) )*( ubar(iend+1,j,knew)
# ifdef WEC 
     &                                              +ust2d(iend+1,j)
# endif
     &                                              )*dn_u(iend+1,j)
        enddo
        do j=jstrV,jend
          DV_avg1(iend+1,j)=DV_avg1(iend+1,j) +cff1*( Dnew(iend+1,j)
     &        +Dnew(iend+1,j-1) )*( vbar(iend+1,j,knew)
# ifdef WEC 
     &                                              +vst2d(iend+1,j)
# endif
     &                                              )*dm_v(iend+1,j)   
       enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istrU,iend
          DU_avg1(i,jstr-1)=DU_avg1(i,jstr-1) +cff1*( Dnew(i,jstr-1)
     &        +Dnew(i-1,jstr-1) )*( ubar(i,jstr-1,knew)
# ifdef WEC
     &                                              +ust2d(i,jstr-1)
# endif
     &                                              )*dn_u(i,jstr-1)
        enddo
        do i=istrR,iendR
          DV_avg1(i,jstrV-1)=DV_avg1(i,jstrV-1)+cff1*(Dnew(i,jstrV-1)
     &         +Dnew(i,jstrV-2))*( vbar(i,jstrV-1,knew)
# ifdef WEC 
     &                                              +vst2d(i,jstrV-1)
# endif
     &                                              )*dm_v(i,jstrV-1)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istrU,iend
          DU_avg1(i,jend+1)=DU_avg1(i,jend+1) +cff1*( Dnew(i,jend+1)
     &        +Dnew(i-1,jend+1) )*( ubar(i,jend+1,knew)
# ifdef WEC 
     &                                               +ust2d(i,jend+1)
# endif
     &                                               )*dn_u(i,jend+1)
        enddo
        do i=istrR,iendR
          DV_avg1(i,jend+1)=DV_avg1(i,jend+1) +cff1*( Dnew(i,jend+1)
     &            +Dnew(i,jend) )*( vbar(i,jend+1,knew)
# ifdef WEC 
     &                                               +vst2d(i,jend+1)
# endif
     &                                               )*dm_v(i,jend+1)
        enddo
      endif
# endif
#endif


#ifdef RIVER_SOURCE
      do j=jstr,jend
        do i=istrU,iend
          if (abs(riv_uflx(i,j)).gt.1e-3) then
            iriver = nint(riv_uflx(i,j)/10)
            river_flux = riv_vol(iriver)*(riv_uflx(i,j)-10*iriver)
            ubar(i,j,knew) = river_flux*
     &           2/( dn_u(i,j)*(Dnew(i-1,j)+Dnew(i,j)) )
# ifdef SOLVE3D
            DU_avg1(i,j) = river_flux
# endif
          endif
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          if (abs(riv_vflx(i,j)).gt.1e-3) then
            iriver = nint(riv_vflx(i,j)/10)
            river_flux = riv_vol(iriver)*(riv_vflx(i,j)-10*iriver)
            vbar(i,j,knew) = river_flux*
     &           2/( dm_v(i,j)*(Dnew(i,j-1)+Dnew(i,j)) )
# ifdef SOLVE3D
            DV_avg1(i,j) = river_flux
# endif
          endif
        enddo
      enddo
#endif   /* RIVER_SOURCE */

! At the end of the last 2D step replace the new free-surface field
! "zeta(:,:,knew)" with its fast-time-averaged "_avg1" version: this
! is the only version of free surface which communicates with 3D part
! via Zt_avg1 --> zeta(:,:,knew) --> z_w(:,:,N) in set_depth and it
! is also used as initial condition for "zeta" for fast-time stepping
! during the next 3D-step.

# ifdef OBC_CHECK
      if(iif == nfast) then
        print *, 'step2d:'
        print *, 'DU_avg_bak(istr,10)=',DU_avg_bak(istr,10)
        print *, 'DU_avg1(istr,10)=',DU_avg1(istr,10)
        print *, 'DU_avg2(istr,10)=',DU_avg2(istr,10)
      endif
# endif

#ifdef SOLVE3D
      if (iif == nfast) then
        do j=jstrR,jendR
          do i=istrR,iendR
            zeta(i,j,knew)=Zt_avg1(i,j)
          enddo
        enddo
        call set_depth_tile(istr,iend,jstr,jend)
      endif
#else
      j=min(iic,2*ninfo)-1          ! This introduces "logarithmic"
      i=1                           ! scaling of diag call frequency
      do while(i < j)               ! at the beginning of model run.
        i=2*i
      enddo
      if (mod(iic-1,min(i,ninfo)) == 0) then
        call diag_tile(istr,iend,jstr,jend, UFx,UFe,VFx)
      endif
#endif
#ifdef EXCHANGE
      call exchange2d_3_tile(istr,iend,jstr,jend,
     &                  zeta(START_2D_ARRAY,knew),
     &                  ubar(START_2D_ARRAY,knew),
     &                  vbar(START_2D_ARRAY,knew))
#endif

      end
