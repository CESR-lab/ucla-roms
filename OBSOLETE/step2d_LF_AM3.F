#include "cppdefs.h"

      subroutine step2D (tile)
      implicit none
      integer tile
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"
      call step2D_LFAM3_tile (istr,iend,jstr,jend, A2d(1,1),
     &                            A2d(1,2),  A2d(1, 3), A2d(1, 4),
     &                  A2d(1,5), A2d(1, 6), A2d(1, 7), A2d(1, 8),
     &                  A2d(1,9), A2d(1,10), A2d(1,11), A2d(1,12))
      end

      subroutine step2D_LFAM3_tile (istr,iend,jstr,jend, zeta_new,
     &                              Dnew, rubar, rvbar, DUon, DVom,
     &                              Drhs, Drhs_p, UFe,UFx,VFe,VFx)
      implicit none
#include "param.h"
      integer istr,iend,jstr,jend, i,j, krhs, kbak
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: zeta_new,Dnew,
     &      rubar,rvbar, DUon,DVom, Drhs, Drhs_p, UFx,UFe,VFe,VFx
      real cff,cff1,cff2,cff3
!
! The following stability limits are obtained empirically using 3/4
! degree Atlantic model configuration. In all these cases barotropic
! mode time step is about 180...250 seconds, which is much less than
! the inertial period. The maximum stability coefficients turned out
! to be slightly different than predicted by linear theory, although
! all theoretical tendencies agree with practice. Note the nearly 70%
! gain in stability in comparison with LF-TR for appropriate set of
! coefficients (linear theory predicts beta=0.166, epsil=0.84).
!
      real, parameter ::
c**  &         gamma=0.,     beta=0.,   epsil=0.    !--> Cu=0.818
c**  &         gamma=1./12., beta=0.,   epsil=0.    !--> Cu=0.878
c*   &         gamma=1./12., beta=0.1,  epsil=0.6   !--> Cu=1.050
     &         gamma=0.,     beta=0.14, epsil=0.74  !==> Cu=1.341

#include "grid.h"
#include "ocean2d.h"
#ifdef SOLVE3D
# include "coupling.h"
#else
# include "forces.h"
#endif
#include "mixing.h"
#include "climat.h"
#include "scalars.h"

#include "compute_auxiliary_bounds.h"

#define PREDICTOR_2D_STEP knew==3
      if (PREDICTOR_2D_STEP) then      ! NOTE: In the following code
        krhs=kstp                      ! it is assumed that variables
      else                             ! with time index "krhs" are
        krhs=3                         ! time-centered at step "n"
      endif                            ! in barotropic time during
      if (FIRST_2D_STEP) then          ! predictor sub-step and
        kbak=kstp                      ! "n+1/2" during corrector.
      else
        kbak=3-kstp                    ! "kbak" is used as "from"
      endif                            ! time index for LF time step.
!
! Preliminary step: compute total depth of the water column and
!------------ ----- vertically integrated mass fluxes, which are
! used in computation of free-surface elevation time tendency and
! in nonlinear advection terms for the barotropic momentum
! equations.
!
# define JR_RANGE jstrV-2,jend+1
# define IR_RANGE istrU-2,iend+1

      do j=JR_RANGE
        do i=IR_RANGE
          Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
        enddo
      enddo
#undef IR_RANGE
#undef JR_RANGE

      do j=jstr-1,jendR
        do i=istrU-1,iend+1
          DUon(i,j)=0.5*(Drhs(i,j)+Drhs(i-1,j))
     &                 *ubar(i,j,krhs)*dn_u(i,j)
        enddo
      enddo
      do j=jstrV-1,jend+1
        do i=istr-1,iendR
          DVom(i,j)=0.5*(Drhs(i,j)+Drhs(i,j-1))
     &                 *vbar(i,j,krhs)*dm_v(i,j)
        enddo
      enddo
#ifdef OBC_VOLCONS
      call set_DUV_bc_tile (istr,iend,jstr,jend, Drhs, DUon,DVom)
#endif

#ifdef SOLVE3D
!
! Fields averaged over all barotropic time steps.  NOTE: The index
!------- -------- ---- --- ---------- ---- ------  ranges here are
! designed to include physical boundaries only. Periodic ghost
! points and internal MPI computational margins are NOT included.
!
                                   ! Reset alltarget arrays for the
                                   ! barotropic mode averages arrays
      if (PREDICTOR_2D_STEP) then  ! during the first predictor step.
        if (FIRST_2D_STEP) then    ! At all subsequent time steps
          do j=jstrR,jendR         ! accumulate averages of the first
            do i=istrR,iendR       ! kind using DELAYED mode, i.e.,
              Zt_avg1(i,j)=0.      ! Zt_avg1 are accumulated not
              DU_avg1(i,j)=0.      ! immediately after corrector step
              DV_avg1(i,j)=0.      ! when they are computed, but
              DU_avg2(i,j)=0.      ! during the subsequent predictor
              DV_avg2(i,j)=0.      ! substep. This allows to save
            enddo                  ! operations because DUon,DVom are
          enddo                    ! computed anyway. The last time
        else                       ! step has special code to add all
          cff=weight(1,iif-1)      ! three barotropic variables after
          do j=jstrR,jendR         ! the last corrector substep.
            do i=istrR,iendR
              Zt_avg1(i,j)=Zt_avg1(i,j)+cff*zeta(i,j,krhs)
              DU_avg1(i,j)=DU_avg1(i,j)+cff*DUon(i,j)
              DV_avg1(i,j)=DV_avg1(i,j)+cff*DVom(i,j)
            enddo
          enddo
        endif
      else                   !--> CORRECTOR_2D_STEP
        cff=weight(2,iif)
        do j=jstrR,jendR
          do i=istrR,iendR
            DU_avg2(i,j)=DU_avg2(i,j)+cff*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff*DVom(i,j)
          enddo
        enddo
      endif
#else
      if (PREDICTOR_2D_STEP) call diag_tile (istr,iend,jstr,jend,
     &                                               UFx,UFe,VFx)
#endif /* SOLVE3D */
!
! Advance free-surface:   Compute zeta_new, which is at new time
!-------- ---- --------   step, and interpolate backward for the
! subsequent computation of barotropic pressure-gradient terms.
! It should be noted that during the predictor of the first 2D step
! in 3D mode, the pressure gradient terms are computed using just
! zeta(:,:,kstp), i.e. like in Forward Euler step, rather than the
! more accurate predictor of generalozed RK2.  This is to keep it
! consistent with computation of pressure gradient in 3D mode, which
! uses exactly the initial value of "zeta", rather than value
! changed by the first barotropic predictor step.  Later in this
! code, just after "rufrc,rvfrc" are finalized, a correction term
! based on the difference zeta_new(:,:)-zeta(:,:,kstp) to
! "rubar,rvbar" to make them consistent with generalozed RK2
! stepping for pressure gradient terms.
!
#define zwrk UFx
#define rzeta  UFe
#define rzeta2  VFe
#define rzetaSA VFx

      if (PREDICTOR_2D_STEP) then
        if (FIRST_2D_STEP) then     ! Modified RK2 time step (with
          cff=dtfast                ! Forward-Backward feedback with
#ifdef SOLVE3D
          cff1=0.                   !==> Forward Euler
          cff2=1.
#else
          cff1=0.333333333333       ! optimally chosen beta=1/3 and
          cff2=0.666666666667       ! epsilon=2/3, see below) is used
#endif
          cff3=0.                   ! here for the start up.
        else
          cff=2.*dtfast             ! In the code below "zwrk" is
          cff1=beta                 ! time-centered at time step "n"
          cff2=1.-2.*beta           ! in the case of LF (for all but
          cff3=beta                 ! the first time step)
        endif

        do j=jstrV-1,jend
          do i=istrU-1,iend
            zeta_new(i,j)=zeta(i,j,kbak) + cff*pm(i,j)*pn(i,j)
     &         *(DUon(i,j)-DUon(i+1,j) + DVom(i,j)-DVom(i,j+1))
#ifdef MASKING
            zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
#endif
            Dnew(i,j)=zeta_new(i,j)+h(i,j)

            zwrk(i,j)=cff1*zeta_new(i,j) + cff2*zeta(i,j,kstp)
     &                                   + cff3*zeta(i,j,kbak)
#if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
#else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
#endif
          enddo
        enddo
      else                              !--> CORRECTOR_2D_STEP
        if (FIRST_2D_STEP) then
          cff=0.333333333333            ! Modified RK2 weighting:
          cff1=0.333333333333           ! here "zwrk" is time-
          cff2=0.333333333333           ! centered at "n+1/2".
          cff3=0.
        else
          cff=1.-epsil                  ! zwrk is always time-
          cff1=(0.5-gamma)*epsil        ! centered at n+1/2
          cff2=(0.5+2.*gamma)*epsil     ! during corrector sub-
          cff3=       -gamma *epsil     ! step.
        endif

        do j=jstrV-1,jend
          do i=istrU-1,iend
            zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &            *(DUon(i,j)-DUon(i+1,j) + DVom(i,j)-DVom(i,j+1))
#ifdef MASKING
            zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
#endif
            Dnew(i,j)=zeta_new(i,j)+h(i,j)

            zwrk(i,j)=cff*zeta(i,j,krhs) + cff1*zeta_new(i,j)
     &              +cff2*zeta(i,j,kstp) + cff3*zeta(i,j,kbak)

#if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
#else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
#endif
          enddo
        enddo
      endif
!
! Set boundary conditions for the free-surface, after which copy
! the newly computed free-surface zeta_new into shared array. It
! should be noted that "zeta_new" is always centered at time step
! "m+1", while zeta(:,:,knew) should be centered either at "m+1/2"
! after predictor step and at "m+1" after corrector. Chosing it
! to be this way makes it possible avoid storing r.h.s. for zeta,
! ubar and vbar between predictor and corrector sub-steps.
!
      call zetabc_tile (istr,iend,jstr,jend, zeta_new)

      if (PREDICTOR_2D_STEP) then
        if (FIRST_2D_STEP) then
          cff1=0.5
          cff2=0.5
          cff3=0.
        else
          cff1=0.5-gamma
          cff2=0.5+2.*gamma
          cff3=      -gamma
        endif
        do j=jstrR,jendR
          do i=istrR,iendR
            zeta(i,j,knew)=cff1*zeta_new(i,j)+cff2*zeta(i,j,kstp)
     &                                       +cff3*zeta(i,j,kbak)
          enddo
        enddo
      else
        do j=jstrR,jendR
          do i=istrR,iendR
            zeta(i,j,knew)=zeta_new(i,j)
          enddo
        enddo
      endif
!
! Compute pressure-gradient terms  NOTE that "rubar" and "rvbar"
!-------- -------- -------- -----  are computed within the same
! fused loop despite the fact that their normal indexe ranges are
! different. Fusing loops causes redundant computation of one
! column of "rubar" on the western physical boundary and one row
! of "rvbar" on the southern, but, at the same time it allows to
! share references to array elements (i,j) which results in an
! increase of computational density by almost a factor of 1.5
! resulting in overall more efficient code pipelined in 26 cycles
! (61% of peak speed) on R10000 vs. 16+16 cycles of separate loop
! version for the case when both CPP switches below are defined.
!
      cff=0.5*g
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=cff*dn_u(i,j)*( (h(i-1,j)+h(i,j))*(rzeta(i-1,j)
     &                        -rzeta(i,j)) +rzeta2(i-1,j)-rzeta2(i,j)

#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i-1,j)-h(i,j))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(zwrk(i-1,j)-zwrk(i,j)))
#endif
     &                                                              )

          rvbar(i,j)=cff*dm_v(i,j)*( (h(i,j-1)+h(i,j))*(rzeta(i,j-1)
     &                        -rzeta(i,j)) +rzeta2(i,j-1)-rzeta2(i,j)

#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(i,j-1)-h(i,j))*( rzetaSA(i,j-1)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(zwrk(i,j-1)-zwrk(i,j)))
#endif
     &                                                              )
        enddo
      enddo            !--> discard  zwrk, rzeta, rzeta2, rzetaSA

#undef rzetaSA
#undef rzeta2
#undef rzeta
#undef zwrk
!
! Suppress computation of advection and viscous terms (lateral and
! bottom drag) in the case of 3D model configuration. Note that this
! does not mean omission of these terms: they are already implicitly
! present within the baroclinic-barotropic coupling terms "rufrc"
! and "rvfrc".  This merely means that advection and viscous terms
! are frozen during barotropic mode time stepping.
!
#ifdef SOLVE3D
# undef UV_ADV
# undef UV_VIS2
#endif
!
! Compute horizontal advection terms for momentum equations (2D only)
!
#ifdef UV_ADV
      do j=jstr,jend
        do i=istrU-1,iend
          UFx(i,j)=0.25*(DUon(i,j)+DUon(i+1,j))
     &        *(ubar(i,j,krhs)+ubar(i+1,j,krhs))
        enddo
      enddo
      do j=jstrV-1,jend
        do i=istr,iend
          VFe(i,j)=0.25*(DVom(i,j)+DVom(i,j+1))
     &        *(vbar(i,j,krhs)+vbar(i,j+1,krhs))
        enddo
      enddo
      do j=jstr,jend+1
        do i=istrU,iend
          UFe(i,j)=0.25*(DVom(i,j)+DVom(i-1,j))
     &        *(ubar(i,j,krhs)+ubar(i,j-1,krhs))
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend+1
          VFx(i,j)=0.25*(DUon(i,j)+DUon(i,j-1))
     &        *(vbar(i,j,krhs)+vbar(i-1,j,krhs))
        enddo
      enddo            !--> discard DUon,DVom
      do j=jstr,jend
        do i=istrU,iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo           !--> discard UFx,UFe,VFx,VFe
#endif
!
! Compute Coriolis (2D and 3D) term and advective curvilinear metric
! terms (2D only).
!
#if defined UV_COR || (defined CURVGRID && defined UV_ADV)
      do j=jstrV-1,jend
        do i=istrU-1,iend
          cff=Drhs(i,j)*(
# ifdef UV_COR
     &                    fomn(i,j)
# endif
# if (defined CURVGRID && defined UV_ADV)
     &  +0.5*( dndx(i,j)*(vbar(i,j,krhs)+vbar(i,j+1,krhs))
     &        -dmde(i,j)*(ubar(i,j,krhs)+ubar(i+1,j,krhs)))
# endif
     &                                                    )
          UFx(i,j)=cff*(vbar(i,j,krhs)+vbar(i,j+1,krhs))
          VFe(i,j)=cff*(ubar(i,j,krhs)+ubar(i+1,j,krhs))
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
      enddo                   !--> discard UFx,UFe,VFx,VFe
#endif /* UV_COR || (CURVGRID && UV_ADV) */
#ifdef UV_VIS2
!
! Compute viscous dissipation terms, if any:  Preliminary step: if
!-------- ------- ----------- ------ -- ----  necessary, fix peroodic
! boundary conditions for "Drhs". Then compute total depth at PSI-
! points "Drhs_p".
!
      do j=jstr,jend+1
        do i=istr,iend+1
          Drhs_p(i,j)=0.25*( Drhs(i,j  )+Drhs(i-1,j  )
     &                      +Drhs(i,j-1)+Drhs(i-1,j-1))
        enddo
      enddo

      do j=jstrV-1,jend
        do i=istrU-1,iend
          cff=2.*Drhs(i,j)*visc2_r(i,j)
          UFx(i,j)=cff*(ubar(i+1,j,krhs)-ubar(i,j,krhs))
     &                                 *pm(i,j)*dn_r(i,j)
          VFe(i,j)=cff*(vbar(i,j+1,krhs)-vbar(i,j,krhs))
     &                                 *pn(i,j)*dm_r(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend+1
          cff=visc2_p(i,j)*Drhs_p(i,j)*0.25*(
     &                 (pn(i,j)+pn(i-1,j)+pn(i,j-1)+pn(i-1,j-1))
     &                        *(ubar(i,j,krhs)-ubar(i,j-1,krhs))
     &                +(pm(i,j)+pm(i-1,j)+pm(i,j-1)+pm(i-1,j-1))
     &                        *(vbar(i,j,krhs)-vbar(i-1,j,krhs))
     &                                                         )
# ifdef MASKING
     &                     *pmask(i,j)
# endif
          UFe(i,j)=cff*dm_p(i,j)
          VFx(i,j)=cff*dn_p(i,j)
        enddo
      enddo
      do j=jstr,jend
        do i=istrU,iend
          rubar(i,j)=rubar(i,j)+UFx(i,j)-UFx(i-1,j)
     &                         +UFe(i,j+1)-UFe(i,j)
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          rvbar(i,j)=rvbar(i,j)+VFx(i+1,j)-VFx(i,j)
     &                         +VFe(i,j)-VFe(i,j-1)
        enddo
      enddo
#endif /* UV_VIS2 */
!
! Add in linear and/or quadratic bottom stress.
!---- -- ------ ------ --------- ------ -------
!
#ifndef SOLVE3D
      if (rdrg2>0.) then
        do j=jstr,jend
          do i=istrU,iend
            cff=0.25*( vbar(i  ,j,krhs)+vbar(i  ,j+1,krhs)
     &                +vbar(i-1,j,krhs)+vbar(i-1,j+1,krhs))
            rubar(i,j)=rubar(i,j)-ubar(i,j,krhs)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,krhs)*ubar(i,j,krhs)+cff*cff)
     &                                     )*dm_u(i,j)*dn_u(i,j)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            cff=0.25*( ubar(i,j  ,krhs)+ubar(i+1,j  ,krhs)
     &                +ubar(i,j-1,krhs)+ubar(i+1,j-1,krhs))
            rvbar(i,j)=rvbar(i,j)-vbar(i,j,krhs)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,krhs)*vbar(i,j,krhs))
     &                                     )*dm_v(i,j)*dn_v(i,j)
          enddo
        enddo
      else if (rdrg>0.0) then
        do j=jstr,jend
          do i=istrU,iend
            rubar(i,j)=rubar(i,j)-rdrg*ubar(i,j,krhs)
     &                            *dm_u(i,j)*dn_u(i,j)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rvbar(i,j)=rvbar(i,j)-rdrg*vbar(i,j,krhs)
     &                            *dm_v(i,j)*dn_v(i,j)
          enddo
        enddo
      endif
#endif
!
! Coupling 2D and 3D parts:  Before the predictor step of the first
!--------- -- --- -- ------  barotropic time step arrays "rufrc" and
! "rvfrc" contain vertically integrals of the 3D right-hand-side
! terms for the momentum equations (including surface and bottom
! stresses, if so prescribed). During the first barotropic time step
! connvert them into forcing terms by subtracting the fast-time
! "rubar" and "rvbar" from them; These forcing terms are then
! extrapolated forward in time using optimized Adams-Bashforth
! weights, so that the resultant "rufrc" and "rvfrc" are centered
! effectively at time n+1/2 in baroclinic time. From now on, these
! newly computed forcing terms remain unchanged during the fast time
! stepping and will added to "rubar" and "rvbar" during all
! subsequent barotropic time steps.
!
#ifdef SOLVE3D
      if (FIRST_2D_STEP .and. PREDICTOR_2D_STEP) then
# ifdef PRED_COUPLED_MODE
        if (FIRST_TIME_STEP) then
          cff3=0.                        ! This version is designed
          cff2=0.                        ! for coupling during 3D
          cff1=1.                        ! predictor sub-step: here
        elseif (FIRST_TIME_STEP+1) then  ! forcing term "rufrc" is
          cff3=0.                        ! computed as instantaneous
          cff2=-0.5                      ! value at 3D time step
          cff1=1.5                       ! "nstp" first, and then
        else                             ! extrapolated half-step
          cff3=0.281105                  ! forward using  AM3-like
          cff2=-0.5-2.*cff3              ! weights optimized for
          cff1=1.5+cff3                  ! maximum stability (with
        endif                            ! special care for startup)

        do j=jstr,jend
          do i=istrU,iend
            cff=rufrc(i,j)-rubar(i,j)
            rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                             + cff3*rufrc_bak(i,j,nstp)
            rufrc_bak(i,j,nstp)=cff
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            cff=rvfrc(i,j)-rvbar(i,j)
            rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                             + cff3*rvfrc_bak(i,j,nstp)
            rvfrc_bak(i,j,nstp)=cff
          enddo
        enddo
# else
        do j=jstr,jend                       ! This version is
          do i=istr,iend                     ! designed for coupling
            rufrc(i,j)=rufrc(i,j)-rubar(i,j) ! during 3D corrector
            rvfrc(i,j)=rvfrc(i,j)-rvbar(i,j) ! sub-step: no forward
          enddo                              ! extrapolation needs
        enddo                                ! to be performed.
# endif
!
! Since coupling requires that pressure gradient term is computed
! using zeta(:,:,kstp) instead of 1/3 toward zeta_new(:,:) as needed
! by generalozed RK2 scheme, apply compensation to shift pressure
! gradient terms from "kstp" to 1/3 toward "knew".
!
# define zwrk UFx
# define rzeta  UFe
# define rzeta2  VFe
# define rzetaSA VFx

        do j=jstrV-1,jend
          do i=istrU-1,iend
            zwrk(i,j)=0.333333333333*(zeta_new(i,j)-zeta(i,j,kstp))
# if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*( 0.333333333333*zeta_new(i,j)
     &                              +1.666666666666*zeta(i,j,kstp))
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*( 0.333333333333*zeta_new(i,j)
     &                              +1.666666666666*zeta(i,j,kstp))
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
      endif   !<-- FIRST_2D_STEP .and. PREDICTOR_2D_STEP
#endif
!
! Perform time step for 2D momenta:   During predictor sub-step,
!-------- ---- ---- --- -- --------   once newly computed "ubar"
! and "vbar" become available, interpolate them half-step backward
! in barotropic time (i.e. they end up time-centered at n+1/2) in
! order to use it during subsequent corrector sub-step.
!
#define Dstp DUon
      if (FIRST_2D_STEP .or. .not. PREDICTOR_2D_STEP) then
        do j=jstrV-1,jend
          do i=istrU-1,iend
            Dstp(i,j)=h(i,j)+zeta(i,j,kstp)
          enddo
        enddo
      else
        do j=jstrV-1,jend
          do i=istrU-1,iend
            Dstp(i,j)=h(i,j)+zeta(i,j,kbak)
          enddo
        enddo
      endif

      if (PREDICTOR_2D_STEP) then
        if (FIRST_2D_STEP) then
          cff=0.5*dtfast
          cff1=0.5
          cff2=0.5
          cff3=0.
        else
          cff=dtfast
          cff1=0.5-gamma
          cff2=0.5+2.*gamma
          cff3=      -gamma
        endif

        do j=jstr,jend
          do i=istrU,iend
            ubar(i,j,knew)=( ubar(i,j,kbak)*(Dstp(i,j)+Dstp(i-1,j))
     &                  +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
#ifdef SOLVE3D
     &                                      *(rubar(i,j)+rufrc(i,j))
#else
     &                               *rubar(i,j) + 4.*cff*sustr(i,j)
#endif
     &                                     )/(Dnew(i,j)+Dnew(i-1,j))
#ifdef MASKING
     &                                                  *umask(i,j)
#endif
            ubar(i,j,knew)=cff1*ubar(i,j,knew)+cff2*ubar(i,j,kstp)
     &                                        +cff3*ubar(i,j,kbak)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            vbar(i,j,knew)=( vbar(i,j,kbak)*(Dstp(i,j)+Dstp(i,j-1))
     &                  +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
#ifdef SOLVE3D
     &                                      *(rvbar(i,j)+rvfrc(i,j))
#else
     &                               *rvbar(i,j) + 4.*cff*svstr(i,j)
#endif
     &                                     )/(Dnew(i,j)+Dnew(i,j-1))
#ifdef MASKING
     &                                                  *vmask(i,j)
#endif
            vbar(i,j,knew)=cff1*vbar(i,j,knew)+cff2*vbar(i,j,kstp)
     &                                        +cff3*vbar(i,j,kbak)
          enddo
        enddo
      else                        !--> CORRECTOR_2D_STEP
        cff=0.5*dtfast
        do j=jstr,jend
          do i=istrU,iend
            ubar(i,j,knew)=( ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))
     &                  +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
#ifdef SOLVE3D
     &                                      *(rubar(i,j)+rufrc(i,j))
#else
     &                               *rubar(i,j) + 4.*cff*sustr(i,j)
#endif
     &                                     )/(Dnew(i,j)+Dnew(i-1,j))
#ifdef MASKING
     &                                                  *umask(i,j)
#endif
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            vbar(i,j,knew)=( vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))
     &                  +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
#ifdef SOLVE3D
     &                                      *(rvbar(i,j)+rvfrc(i,j))
#else
     &                               *rvbar(i,j) + 4.*cff*svstr(i,j)
#endif
     &                                     )/(Dnew(i,j)+Dnew(i,j-1))
#ifdef MASKING
     &                                                  *vmask(i,j)
#endif
          enddo
        enddo
      endif
#undef Dstp
c??
c??
c      do j=jstr,jend
c        do i=istr,iend
c          ubar(i,j,knew)=0.
c          vbar(i,j,knew)=0.
c        enddo
c      enddo
c??
c??

!
! Set boundary conditions and compute integral mass flux accross
! all open boundaries, if any.
!
      call    u2dbc_tile (istr,iend,jstr,jend, UFx)
      call    v2dbc_tile (istr,iend,jstr,jend, UFx)
#ifdef OBC_VOLCONS
      call obc_flux_tile (istr,iend,jstr,jend, UFx)
#endif
#ifdef PSOURCE
      do is=1,Nsrc           ! Apply point mass sources associated
        i=Isrc(is)           ! with river runoff simulations.
        j=Jsrc(is)
        if (istrR<=i .and. i<=iendR .and.
     &      jstrR<=j .and. j<=jendR) then
          if (Dsrc(is)==0) then
            ubar(i,j,knew)=2.*Qbar(is)/( dn_u(i,j)
     &                   *(Dnew(i-1,j)+Dnew(i,j)) )
          else
            vbar(i,j,knew)=2.*Qbar(is)/( dm_v(i,j)
     &                   *(Dnew(i,j-1)+Dnew(i,j)) )
          endif
        endif
      enddo
#endif
#ifdef SOLVE3D

! Finalize computation of barotropic mode averages. This procedure
! starts with filling in boundary rows of total depths at the new
! time step, which is needed to be done only during the last
! barotropic time step, since normally computation of averages
! occurs at the beginning of the next pedictor step (because DUon
! DVom are being computed any way). Strictly speaking, filling the
! boundaries is necessary only in the case of open boundaries,
! since otherwise the associated fluxes are all zeros.

      if (iif==nfast .and. knew<3) then
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

        cff=weight(1,iif)
        cff1=0.5*cff
        do j=jstrR,jendR
          do i=istrR,iendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff*zeta(i,j,knew)
            DU_avg1(i,j)=DU_avg1(i,j)+cff1*ubar(i,j,knew)
     &                 *(Dnew(i,j)+Dnew(i-1,j))*dn_u(i,j)
            DV_avg1(i,j)=DV_avg1(i,j)+cff1*vbar(i,j,knew)
     &                 *(Dnew(i,j)+Dnew(i,j-1))*dm_v(i,j)


! At the end of the last 2D step replace the new free-surface field
! "zeta(:,:,knew)" with its fast-time-averaged "_avg1" version: this
! is the only version of free surface which communicates with 3D part
! via Zt_avg1 --> zeta(:,:,knew) --> z_w(:,:,N) in set_depth and it
! is also used as initial condition for "zeta" for fast-time stepping
! during the next 3D-step.

            zeta(i,j,knew)=Zt_avg1(i,j)
          enddo
        enddo
        call set_depth_tile (istr,iend,jstr,jend, UFx)
      endif
#endif

#ifdef EXCHANGE
      call exchange2d_3_tile (istr,iend,jstr,jend,
     &                   zeta(START_2D_ARRAY,knew),
     &                   ubar(START_2D_ARRAY,knew),
     &                   vbar(START_2D_ARRAY,knew))
#endif
      end
