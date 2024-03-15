#include "cppdefs.opt"

      subroutine step2d

! Advance barotropic mode variables (2D momenta and free surface) by
! one barotropic time step using Generalized Forward-Backward AB3-AM4
! stepping algorithm [Sec. 2.3 from SM2005 starting with Eq. (2.49)]
! and perform fast-time averaging to interact with baroclinic mode.

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none
      integer,save :: tile=0

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
      use tides
      use coupling
      use grid
      use mixing
      use ocean_vars
      use scalars
      use mpi_exchanges

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
        bkw_new=0. ; bkw=1.              ! kold  kbak  kstp  knew
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
          if (pipe_source) then          ! If pipe_source is false, the branch will
            if (pipe_idx(i,j) > 0.) then ! be eliminated by compile optimization
              zeta_new(i,j)=zeta_new(i,j) + dtfast*pm(i,j)*pn(i,j)*pipe_flx(i,j)
            endif
          endif

#ifdef MASKING
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j) ! not sure this is needed
#endif
          Dnew(i,j)=zeta_new(i,j)+h(i,j)

          zwrk(i,j)=bkw_new*zeta_new(i,j) +bkw*zeta(i,j,kstp)
     &              +bkw1*zeta(i,j,kbak) +bkw2*zeta(i,j,kold)

#if defined VAR_RHO_2D
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

#if defined VAR_RHO_2D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(zwrk(i-1,j)-zwrk(i,j)) )
#endif

#ifdef WEC
     &                  + ( h(i-1,j)+h(i,j)+rzeta(i-1,j)+rzeta(i,j) )
     &                                       *( sup(i,j)-sup(i-1,j) )
#endif
     &                                                              )
!>
          rvbar(i,j)=cff*dm_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)

#if defined VAR_RHO_2D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(zwrk(i,j-1)-zwrk(i,j)) )
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

# ifdef CORR_COUPLED_MODE
#  ifdef EXTRAP_BAR_FLUXES
#   undef UV_ADV
#  endif
# else
#  undef UV_ADV
# endif


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

!! Here, we remove all the terms that are in 2d rhs
            rufrc(i,j)=rufrc(i,j)-rubar(i,j)  ! Barotropic mode runnig
            rvfrc(i,j)=rvfrc(i,j)-rvbar(i,j)  ! there is no need for
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
# if defined VAR_RHO_2D
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
# if defined VAR_RHO_2D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(zwrk(i-1,j)-zwrk(i,j)) )
# endif
     &                                                              )
!>
            rvbar(i,j)=rvbar(i,j) +cff*dm_v(i,j)*( (h(i,j-1)+h(i,j))
     &          *(rzeta(i,j-1)-rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)
# if defined VAR_RHO_2D
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

      cff =0.5*dtfast
      cff1=0.5*weight(1,iif)

      do j=jstr,jend
        do i=istrU,iend
          DUnew=( (Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
c<<     &     +cff*iA_u(i,j)
     &                            *(rubar(i,j)+rufrc(i,j))
     &                                                   )
#ifdef MASKING
     &                                         *umask(i,j)
#endif
          ubar(i,j,knew)=DUnew/( Dnew(i,j)+Dnew(i-1,j)
     &                                                   )

          DU_avg1(i,j)=DU_avg1(i,j)+cff1*dn_u(i,j)*( DUnew
# ifdef WEC
     &                 +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
# endif
     &                                                   )

        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          DVnew=( (Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
c<<     &     +cff*iA_v(i,j)
     &                            *(rvbar(i,j)+rvfrc(i,j))
     &                                                   )
#ifdef MASKING
     &                                         *vmask(i,j)
#endif
          vbar(i,j,knew)=DVnew/( Dnew(i,j)+Dnew(i,j-1)
     &                                                   )
 
          DV_avg1(i,j)=DV_avg1(i,j) +cff1*dm_v(i,j) * ( DVnew
# ifdef WEC 
     &                 +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
# endif
     &                                                   )
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

      if (river_source) then
        do j=jstr,jend
          do i=istrU,iend
            if (abs(riv_uflx(i,j)).gt.1e-3) then
              iriver = nint(riv_uflx(i,j)/10)
              river_flux = riv_vol(iriver)*(riv_uflx(i,j)-10*iriver)
              ubar(i,j,knew) = river_flux*
     &             2/( dn_u(i,j)*(Dnew(i-1,j)+Dnew(i,j)) )
              DU_avg1(i,j) = river_flux
            endif
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            if (abs(riv_vflx(i,j)).gt.1e-3) then
              iriver = nint(riv_vflx(i,j)/10)
              river_flux = riv_vol(iriver)*(riv_vflx(i,j)-10*iriver)
              vbar(i,j,knew) = river_flux*
     &             2/( dm_v(i,j)*(Dnew(i,j-1)+Dnew(i,j)) )
              DV_avg1(i,j) = river_flux
            endif
          enddo
        enddo
      endif ! <-- river_source

! At the end of the last 2D step replace the new free-surface field
! "zeta(:,:,knew)" with its fast-time-averaged "_avg1" version: this
! is the only version of free surface which communicates with 3D part
! via Zt_avg1 --> zeta(:,:,knew) --> z_w(:,:,N) in set_depth and it
! is also used as initial condition for "zeta" for fast-time stepping
! during the next 3D-step.

      if (iif == nfast) then
        do j=jstrR,jendR
          do i=istrR,iendR
            zeta(i,j,knew)=Zt_avg1(i,j)
          enddo
        enddo
        call set_depth_tile(istr,iend,jstr,jend)
      endif

#ifdef EXCHANGE
      call exchange_xxx(zeta(:,:,knew),ubar(:,:,knew),vbar(:,:,knew))
#endif

      end
