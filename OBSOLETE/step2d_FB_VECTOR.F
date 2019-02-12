#include "cppdefs.h"

c**   module step2d_mod
c**   contains

                                ! This is an experimental rewrite of
      subroutine step2d (tile)  ! the standard loop-syntax routine
      implicit none             ! into F90 vector notation syntax for
      integer tile              ! performance evaluation purposes
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"
      call step2d_FB_tile (istr,iend,jstr,jend, A2d(1,1),  A2d(1,2),
     &                    A2d(1, 3), A2d(1, 4), A2d(1, 5), A2d(1, 6),
     &                    A2d(1, 7), A2d(1, 8), A2d(1, 9), A2d(1,10),
     &                               A2d(1,11), A2d(1,12), A2d(1,13))
      end subroutine step2d

      subroutine step2d_FB_tile (istr,iend,jstr,jend, zeta_new,Dnew,
     &                           rubar,rvbar, urhs,vrhs,  DUon,DVom,
     &                                        Drhs, UFx,UFe,VFx,VFe)

! Prform one time step for barotropic mode (free-surface and baro-
! tropic 2D momentum equations) using Generalized Forward-Backward
! AB3-AM4 algorithm. Also calculate fast-time averages to interact
! with baroclinic mode.

      implicit none
#include "param.h"
      integer istr,iend,jstr,jend, i,j, kbak, kold
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: zeta_new, Dnew,
     &                         rubar,rvbar,  urhs,vrhs,  DUon,DVom,
     &                                       Drhs, UFx,UFe,VFx,VFe
      real cff, cff0,cff1,cff2,cff3,  DUnew,DVnew
#ifdef PSOURCE
      integer is
#endif
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
#include "sources.h"

#include "compute_auxiliary_bounds.h"

! Preliminary step: compute total depth (meters) of the water column
! ----------------- and vertically-integrated mass fluxes needed for
! computing divergence in rhs_zeta terms and for nonlinear advection
! terms for the barotropic momentum equations.

      if (FIRST_2D_STEP) then         ! Meaning of temporal indices
        kbak=kstp                     ! ------- -- -------- -------
        kold=kstp                     ! m-2     m-1      m      m+1
        cff1=1.                       ! kold    kbak    kstp    knew
        cff2=0.
        cff3=0.
      elseif (FIRST_2D_STEP+1) then
        kbak=kstp-1
        if (kbak<1) kbak=4
        kold=kbak

cxx     cff1= 1.5
cxx     cff2=-0.5
cxx     cff3= 0.

        cff1=1.                  ! Logically AB2-AM3 forward-backward
        cff2=0.                  ! scheme with coefficients chosen for
        cff3=0.                  ! maximum stability ... (see below)
      else
        kbak=kstp-1
        if (kbak<1) kbak=4
        kold=kstp-2
        if (kold<1) kold=4

c       cff1= 1.72
c       cff2=-0.94
c       cff3= 0.22

        cff1= 1.781105
        cff2=-1.06221
        cff3= 0.281105
      endif

          Drhs(istrU-2:iend+1,jstrV-2:jend+1)=
     &                   h(istrU-2:iend+1,jstrV-2:jend+1)
     &          +cff1*zeta(istrU-2:iend+1,jstrV-2:jend+1,kstp)
     &          +cff2*zeta(istrU-2:iend+1,jstrV-2:jend+1,kbak)
     &          +cff3*zeta(istrU-2:iend+1,jstrV-2:jend+1,kold)


          urhs(istrU-1:iend+1,jstr-1:jend+1)=
     &            cff1*ubar(istrU-1:iend+1,jstr-1:jend+1,kstp)
     &           +cff2*ubar(istrU-1:iend+1,jstr-1:jend+1,kbak)
     &           +cff3*ubar(istrU-1:iend+1,jstr-1:jend+1,kold)

          DUon(istrU-1:iend+1,jstr-1:jend+1)=
     &                0.5*( Drhs(istrU-1:iend+1,jstr-1:jend+1)
     &                     +Drhs(istrU-2:iend  ,jstr-1:jend+1)
     &                    )*dn_u(istrU-1:iend+1,jstr-1:jend+1)
     &                     *urhs(istrU-1:iend+1,jstr-1:jend+1)


          vrhs(istr-1:iend+1,jstrV-1:jend+1)=
     &            cff1*vbar(istr-1:iend+1,jstrV-1:jend+1,kstp)
     &           +cff2*vbar(istr-1:iend+1,jstrV-1:jend+1,kbak)
     &           +cff3*vbar(istr-1:iend+1,jstrV-1:jend+1,kold)

          DVom(istr-1:iend+1,jstrV-1:jend+1)=
     &                0.5*( Drhs(istr-1:iend+1,jstrV-1:jend+1)
     &                     +Drhs(istr-1:iend+1,jstrV-2:jend  )
     &                    )*dm_v(istr-1:iend+1,jstrV-1:jend+1)
     &                     *vrhs(istr-1:iend+1,jstrV-1:jend+1)



#ifdef OBC_VOLCONS
      call set_DUV_bc_tile (istr,iend,jstr,jend, Drhs, DUon,DVom)
#endif

      if (FIRST_2D_STEP) then
        cff0=1.
        cff1=0.
        cff2=0.
        cff3=0.
      elseif (FIRST_2D_STEP+1) then
        cff0=0.779
        cff1=-0.058
        cff2=0.279
        cff3=0.
      else
c       cff0=0.6133333333333
c       cff1=0.2821333333333
c       cff2=0.0957333333333
c       cff3=0.0088

        cff0=0.614
        cff1=0.285
        cff2=0.088
        cff3=0.013
      endif

#define zwrk UFx
#define rzeta  UFe
#define rzeta2  VFe
#define rzetaSA VFx


          zeta_new(istrU-1:iend,jstrV-1:jend)=
     &       zeta(istrU-1:iend,jstrV-1:jend,kstp)
     &       +dtfast*pm(istrU-1:iend,jstrV-1:jend)
     &              *pn(istrU-1:iend,jstrV-1:jend)
     &            *( DUon(istrU-1:iend,jstrV-1:jend)
     &              -DUon(istrU:iend+1,jstrV-1:jend)
     &              +DVom(istrU-1:iend,jstrV-1:jend)
     &              -DVom(istrU-1:iend,jstrV:jend+1)
     &                                             )
#ifdef MASKING
          zeta_new(istrU-1:iend,jstrV-1:jend)=
     &           zeta_new(istrU-1:iend,jstrV-1:jend)
     &             *rmask(istrU-1:iend,jstrV-1:jend)
#endif
          Dnew(istrU-1:iend,jstrV-1:jend)=
     &           zeta_new(istrU-1:iend,jstrV-1:jend)
     &                 +h(istrU-1:iend,jstrV-1:jend)

          zwrk(istrU-1:iend,jstrV-1:jend)=
     &       cff0*zeta_new(istrU-1:iend,jstrV-1:jend)
     &          +cff1*zeta(istrU-1:iend,jstrV-1:jend,kstp)
     &          +cff2*zeta(istrU-1:iend,jstrV-1:jend,kbak)
     &          +cff3*zeta(istrU-1:iend,jstrV-1:jend,kold)

#if defined VAR_RHO_2D && defined SOLVE3D
          rzeta(istrU-1:iend,jstrV-1:jend)=
     &          (1.+rhoS(istrU-1:iend,jstrV-1:jend))
     &             *zwrk(istrU-1:iend,jstrV-1:jend)

          rzeta2(istrU-1:iend,jstrV-1:jend)=
     &             rzeta(istrU-1:iend,jstrV-1:jend)
     &             *zwrk(istrU-1:iend,jstrV-1:jend)

          rzetaSA(istrU-1:iend,jstrV-1:jend)=
     &              zwrk(istrU-1:iend,jstrV-1:jend)
     &           *( rhoS(istrU-1:iend,jstrV-1:jend)
     &             -rhoA(istrU-1:iend,jstrV-1:jend))
#else
          rzeta(istrU-1:iend,jstrV-1:jend)=
     &              zwrk(istrU-1:iend,jstrV-1:jend)
          rzeta2(istrU-1:iend,jstrV-1:jend)=
     &              zwrk(istrU-1:iend,jstrV-1:jend)
     &             *zwrk(istrU-1:iend,jstrV-1:jend)
#endif





!
! Set boundary conditions for the free-surface, after which copy
! the newly computed free-surface zeta_new into shared array.
!
      call zetabc_tile (istr,iend,jstr,jend, zeta_new)

          zeta(istrR:iendR,jstrR:jendR,knew)=
     &       zeta_new(istrR:iendR,jstrR:jendR)

! Compute fast-time-averaged fields over all short timesteps.
!-------- ---- ---- -------- ------ ---- --- ----- ----------
! Reset/initialise arrays for averaged fields during the first
! barotropic time step; Accumulate averages after that. Include
! physical boundary points, but not periodic ghost points or
! computation  MPI computational margins.

#ifdef SOLVE3D
        cff1=weight(1,iif)
        cff2=weight(2,iif)

        if (FIRST_2D_STEP) then
          do j=jstrR,jendR
            do i=istrR,iendR
              Zt_avg1(i,j)=cff1*zeta(i,j,knew)
              DU_avg1(i,j)=0.
              DV_avg1(i,j)=0.
              DU_avg2(i,j)=cff2*DUon(i,j)
              DV_avg2(i,j)=cff2*DVom(i,j)
            enddo
          enddo
        else
          do j=jstrR,jendR
            do i=istrR,iendR
              Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
              DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
              DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
            enddo
          enddo
        endif
#endif
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

      rubar(istr:iend,jstr:jend)=cff*dn_u(istr:iend,jstr:jend)*(
     &       (h(istr-1:iend-1,jstr:jend)+h(istr:iend,jstr:jend))
     &   *(rzeta(istr-1:iend-1,jstr:jend)-rzeta(istr:iend,jstr:jend))
     &   +rzeta2(istr-1:iend-1,jstr:jend)-rzeta2(istr:iend,jstr:jend)

#if defined VAR_RHO_2D && defined SOLVE3D
     &              +(h(istr-1:iend-1,jstr:jend)-h(istr:iend,jstr:jend))*( rzetaSA(i-1,j)+rzetaSA(i,j)
     &                        +0.333333333333*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(zwrk(i-1,j)-zwrk(i,j)))
#endif
     &                                                              )

      rvbar(istr:iend,jstr:jend)=cff*dm_v(istr:iend,jstr:jend)*(
     &            (h(istr:iend,jstr-1:jend-1)+h(istr:iend,jstr:jend))
     &   *(rzeta(istr:iend,jstr-1:jend-1)-rzeta(istr:iend,jstr:jend))
     &   +rzeta2(istr:iend,jstr-1:jend-1)-rzeta2(istr:iend,jstr:jend)

#if defined VAR_RHO_2D && defined SOLVE3D
     &           +(h(istr:iend,jstr-1:jend-1)-h(istr:iend,jstr:jend))
     &    *( rzetaSA(istr:iend,jstr-1:jend-1)
     &             +rzetaSA(istr:iend,jstr:jend)
     &             +0.333333333333*(rhoA(istr:iend,jstr-1:jend-1)
     &               -rhoA(istr:iend,jstr:jend))
     &    *(zwrk(istr:iend,jstr-1:jend-1)-zwrk(istr:iend,jstr:jend)))
#endif
     &                                                              )
                     !--> discard  zwrk, rzeta, rzeta2, rzetaSA

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
        UFx(istr-1:iend,jstr:jend)=0.25*( DUon(istr-1:iend,jstr:jend)
     &                                   +DUon(istr:iend+1,jstr:jend)
     &                                )*( urhs(istr-1:iend,jstr:jend)
     &                                   +urhs(istr:iend+1,jstr:jend)
     &                                                              )

        VFx(istr:iend+1,jstr:jend)=0.25*( DUon(istr:iend+1,jstr:jend)
     &                               +DUon(istr:iend+1,jstr-1:jend-1)
     &                                )*( vrhs(istr-1:iend,jstr:jend)
     &                                   +vrhs(istr:iend+1,jstr:jend)
     &                                                              )

        VFe(istr:iend,jstr-1:jend)=0.25*( DVom(istr:iend,jstr-1:jend)
     &                                   +DVom(istr:iend,jstr:jend+1)
     &                                )*( vrhs(istr:iend,jstr-1:jend)
     &                                   +vrhs(istr:iend,jstr:jend+1)
     &                                                              )

        UFe(istr:iend,jstr:jend+1)=0.25*( DVom(istr:iend,jstr:jend+1)
     &                               +DVom(istr-1:iend-1,jstr:jend+1)
     &                                )*( urhs(istr:iend,jstr-1:jend)
     &                                   +urhs(istr:iend,jstr:jend+1)
     &                                                              )

        rubar(istr:iend,jstr:jend)=rubar(istr:iend,jstr:jend)
     &         -UFx(istr:iend,jstr:jend)+UFx(istr-1:iend-1,jstr:jend)
     &         -UFe(istr:iend,jstr+1:jend+1)+UFe(istr:iend,jstr:jend)

        rvbar(istr:iend,jstr:jend)=rvbar(istr:iend,jstr:jend)
     &         -VFx(istr+1:iend+1,jstr:jend)+VFx(istr:iend,jstr:jend)
     &         -VFe(istr:iend,jstr:jend)+VFe(istr:iend,jstr-1:jend-1)

                      !--> discard UFx,VFe,UFe,VFx, DUon,DVom
#endif
!
! Compute Coriolis (2D and 3D) term and advective curvilinear metric
! terms (2D only).
!
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)

#define cff UFe

        cff(istrU-1:iend,jstrV-1:jend)=
     &                Drhs(istrU-1:iend,jstrV-1:jend)*(
#  ifdef UV_COR
     &                   fomn(istrU-1:iend,jstrV-1:jend)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &      +0.5*( dndx(istrU-1:iend,jstrV-1:jend)
     &             *( vrhs(istrU-1:iend,jstrV-1:jend)
     &               +vrhs(istrU-1:iend,jstrV:jend+1))
     &            -dmde(istrU-1:iend,jstrV-1:jend)
     &             *( urhs(istrU-1:iend,jstrV-1:jend)
     &               +urhs(istrU:iend+1,jstrV-1:jend)))
#  endif
     &                                              )
        UFx(istrU-1:iend,jstrV-1:jend)=
     &                 cff(istrU-1:iend,jstrV-1:jend)*(
     &                   vrhs(istrU-1:iend,jstrV-1:jend)
     &                  +vrhs(istrU-1:iend,jstrV:jend+1)
     &                                                 )
        VFe(istrU-1:iend,jstrV-1:jend)=
     &                 cff(istrU-1:iend,jstrV-1:jend)*(
     &                   urhs(istrU-1:iend,jstrV-1:jend)
     &                  +urhs(istrU:iend+1,jstrV-1:jend)
     &                                                 )

#undef cff

        rubar(istrU:iend,jstr:jend)=rubar(istrU:iend,jstr:jend)
     &                         +0.25*( UFx(istrU:iend,jstr:jend)
     &                           +UFx(istrU-1:iend-1,jstr:jend))

        rvbar(istr:iend,jstrV:jend)=rvbar(istr:iend,jstrV:jend)
     &                         -0.25*( VFe(istr:iend,jstrV:jend)
     &                           +VFe(istr:iend,jstrV-1:jend-1))
# endif
!
! Compute horizontal viscous stress terms (2D only).
!
# ifdef UV_VIS2
      do j=jstrV-1,jend
        do i=istrU-1,iend
          cff=2.*Drhs(i,j)*visc2_r(i,j)
          UFx(i,j)=cff*(urhs(i+1,j)-urhs(i,j))*pm(i,j)*dn_r(i,j)
          VFe(i,j)=cff*(vrhs(i,j+1)-vrhs(i,j))*pn(i,j)*dm_r(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend+1
          cff=visc2_p(i,j)*Drhs_p(i,j)*0.25*(
     &                 (pn(i,j)+pn(i-1,j)+pn(i,j-1)+pn(i-1,j-1))
     &                        *(urhs(i,j)-urhs(i,j-1))
     &                +(pm(i,j)+pm(i-1,j)+pm(i,j-1)+pm(i-1,j-1))
     &                        *(vrhs(i,j)-vrhs(i-1,j))
     &                                                         )
#  ifdef MASKING
     &                     *pmask(i,j)
#  endif
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
# endif /* UV_VIS2 */
!
! Linear and/or quadratic bottom stress.
!
#ifndef SOLVE3D
      if (rdrg2>0.) then
        do j=jstr,jend
          do i=istrU,iend
            cff=0.25*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))

            rubar(i,j)=rubar(i,j)-ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                                     )*dm_u(i,j)*dn_u(i,j)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            cff=0.25*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))

            rvbar(i,j)=rvbar(i,j)-vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                                     )*dm_v(i,j)*dn_v(i,j)
          enddo
        enddo
      else if (rdrg>0.0) then
        do j=jstr,jend
          do i=istrU,iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *dm_u(i,j)*dn_u(i,j)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *dm_v(i,j)*dn_v(i,j)
          enddo
        enddo
      endif
#endif

! Coupling between 2D and 3D parts.
!--------- ------- -- --- -- ------
! Before the predictor step of the first barotropic time step
! arrays "rufrc" and "rvfrc" contain vertically integrals of the
! 3D right-hand-side terms for the momentum equations (including
! surface and bottom stresses, if so prescribed).
!
! During the first barotropic time step connvert them into forcing
! terms by subtracting the fast-time "rubar" and "rvbar" from them;
! These forcing terms are then extrapolated forward in time using
! optimized Adams-Bashforth weights, so that the resultant rufrc
! and rvfrc are centered effectively at time n+1/2. From now on,
! these newly computed forcing terms will remain constant during
! the fast time stepping and will added to "rubar" and "rvbar"
! during all subsequent barotropic time steps.

#ifdef SOLVE3D
      if (FIRST_2D_STEP) then
/*      if (FIRST_TIME_STEP) then
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
*/
        do j=jstr,jend                       ! This version is
          do i=istr,iend                     ! designed for coupling
            rufrc(i,j)=rufrc(i,j)-rubar(i,j) ! during 3D corrector
            rvfrc(i,j)=rvfrc(i,j)-rvbar(i,j) ! sub-step: no forward
          enddo                              ! extrapolation needs
        enddo                                ! to be performed.
      endif
#endif

! Perform time step for the 2D momentum equations. Also compute
! fast-time averaged barotropic mass fluxes. Doing so on the fly
! yields a more computationally dense code and eliminates repeated
! multiplication by Dnew (since mass fluxes are actually available
! as volatile variables DUnew,DVnew at this moment. However, as the
! result of this arrangement, a special code is needed to compute
! fast-time averages along the physical boundaries, which is done
! immediately below.

#define Dstp DUon
#define DUnew DVom
#define DVnew DVom

          Dstp(istrU-1:iend,jstrV-1:jend)=+h(istrU-1:iend,jstrV-1:jend)
     &                            +zeta(istrU-1:iend,jstrV-1:jend,kstp)

      cff=0.5*dtfast
#ifdef SOLVE3D
      cff1=0.5*weight(1,iif)
#else
      cff2=2.*dtfast
#endif

          DUnew(istrU:iend,jstr:jend)=( (Dstp(istrU:iend,jstr:jend)
     &                              +Dstp(istrU-1:iend-1,jstr:jend)
     &                               )*ubar(istrU:iend,jstr:jend,kstp)
     &     +cff*(pm(istrU:iend,jstr:jend)+pm(istrU-1:iend-1,jstr:jend))
     &         *(pn(istrU:iend,jstr:jend)+pn(istrU-1:iend-1,jstr:jend))
#ifdef SOLVE3D
     &        *(rubar(istrU:iend,jstr:jend)+rufrc(istrU:iend,jstr:jend))
#else
     &     *rubar(istrU:iend,jstr:jend)+cff2*sustr(istrU:iend,jstr:jend)
#endif
     &                                                    )
#ifdef MASKING
     &                                *umask(istrU:iend,jstr:jend)
#endif
           ubar(istrU:iend,jstr:jend,knew)=DUnew(istrU:iend,jstr:jend)
     &     /(Dnew(istrU:iend,jstr:jend)+Dnew(istrU-1:iend-1,jstr:jend))
#ifdef SOLVE3D
          DU_avg1(istrU:iend,jstr:jend)=DU_avg1(istrU:iend,jstr:jend)
     &    +cff1*DUnew(istrU:iend,jstr:jend)*dn_u(istrU:iend,jstr:jend)
#endif


          DVnew(istr:iend,jstrV:jend)=( ( Dstp(istr:iend,jstrV:jend)
     &                                   +Dstp(istr:iend,jstrV-1:jend-1)
     &                                 )*vbar(istr:iend,jstrV:jend,kstp)
     &      +cff*(pm(istr:iend,jstrV:jend)+pm(istr:iend,jstrV-1:jend-1))
     &          *(pn(istr:iend,jstrV:jend)+pn(istr:iend,jstrV-1:jend-1))
#ifdef SOLVE3D
     &        *(rvbar(istr:iend,jstrV:jend)+rvfrc(istr:iend,jstrV:jend))
#else
     &     *rvbar(istr:iend,jstrV:jend)+cff2*svstr(istr:iend,jstrV:jend)
#endif
     &                                                    )
#ifdef MASKING
     &                                      *vmask(istr:iend,jstrV:jend)
#endif
          vbar(istr:iend,jstrV:jend,knew)=DVnew(istr:iend,jstrV:jend)
     &      /(Dnew(istr:iend,jstrV:jend)+Dnew(istr:iend,jstrV-1:jend-1))
#ifdef SOLVE3D
          DV_avg1(istr:iend,jstrV:jend)=DV_avg1(istr:iend,jstrV:jend)
     &      +cff1*DVnew(istr:iend,jstrV:jend)*dm_v(istr:iend,jstrV:jend)
#endif

! Set boundary conditions and compute integral mass flux accross
! all open boundaries, if any.

      call    u2dbc_tile (istr,iend,jstr,jend, UFx)
      call    v2dbc_tile (istr,iend,jstr,jend, UFx)
#ifdef OBC_VOLCONS
      call obc_flux_tile (istr,iend,jstr,jend, UFx)
#endif

! Compute fast-time averaged barotropic mass fluxes along physical
! boundaries.

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
     &         +Dnew(istrU-2,j))*ubar(istrU-1,j,knew)*dn_u(istrU-1,j)
        enddo
        do j=jstrV,jend
          DV_avg1(istr-1,j)=DV_avg1(istr-1,j) +cff1*(Dnew(istr-1,j)
     &       +Dnew(istr-1,j-1) )*vbar(istr-1,j,knew)*dm_v(istr-1,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstrR,jendR
          DU_avg1(iend+1,j)=DU_avg1(iend+1,j) +cff1*( Dnew(iend+1,j)
     &            +Dnew(iend,j) )*ubar(iend+1,j,knew)*dn_u(iend+1,j)
        enddo
        do j=jstrV,jend
          DV_avg1(iend+1,j)=DV_avg1(iend+1,j) +cff1*( Dnew(iend+1,j)
     &        +Dnew(iend+1,j-1) )*vbar(iend+1,j,knew)*dm_v(iend+1,j)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istrU,iend
          DU_avg1(i,jstr-1)=DU_avg1(i,jstr-1) +cff1*( Dnew(i,jstr-1)
     &        +Dnew(i-1,jstr-1) )*ubar(i,jstr-1,knew)*dn_u(i,jstr-1)
        enddo
        do i=istrR,iendR
          DV_avg1(i,jstrV-1)=DV_avg1(i,jstrV-1)+cff1*(Dnew(i,jstrV-1)
     &         +Dnew(i,jstrV-2))*vbar(i,jstrV-1,knew)*dm_v(i,jstrV-1)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istrU,iend
          DU_avg1(i,jend+1)=DU_avg1(i,jend+1) +cff1*( Dnew(i,jend+1)
     &        +Dnew(i-1,jend+1) )*ubar(i,jend+1,knew)*dn_u(i,jend+1)
        enddo
        do i=istrR,iendR
          DV_avg1(i,jend+1)=DV_avg1(i,jend+1) +cff1*( Dnew(i,jend+1)
     &            +Dnew(i,jend) )*vbar(i,jend+1,knew)*dm_v(i,jend+1)
        enddo
      endif
# endif
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
# ifdef SOLVE3D
            DU_avg1(i,j)=Qbar(is)
# endif
          else
            vbar(i,j,knew)=2.*Qbar(is)/( dm_v(i,j)
     &                   *(Dnew(i,j-1)+Dnew(i,j)) )
# ifdef SOLVE3D
            DV_avg1(i,j)=Qbar(is)
# endif
          endif
        endif
      enddo
#endif
#ifdef SOLVE3D
      if (iif==nfast) call set_depth_tile (istr,iend,jstr,jend,UFx)
#else
      call diag_tile (istr,iend,jstr,jend, UFx,UFe,VFx)
#endif
#ifdef EXCHANGE
      call exchange2d_3_tile (istr,iend,jstr,jend,
     &                   zeta(START_2D_ARRAY,knew),
     &                   ubar(START_2D_ARRAY,knew),
     &                   vbar(START_2D_ARRAY,knew))

#endif
      end subroutine step2d_FB_tile

c**   end module step2d_mod
