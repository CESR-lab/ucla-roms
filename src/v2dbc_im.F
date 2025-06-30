#include "cppdefs.opt"

      subroutine v2dbc_tile (istr,iend,jstr,jend, grad)

! Set lateral boundary conditions for the barotropic (i.e. vertically
! integrated) ETA-component velocity vbar(:,:,knew).

      use param
#ifdef WEC
      use wec_frc
#endif
      use boundary
      use grid
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j
      real grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cext, cff,
     &                                     dft,dfx,dfy, zx,hx
      real, parameter :: eps=1.D-33
#ifdef OBC_FLUX_CORR
      real uv_crss
#endif

#include "compute_auxiliary_bounds.h"

#ifdef OBC_FLUX_CORR
      uv_crss=ubar_xs*dtfast/tau_obc
#endif

#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
# ifdef OBC_SOUTH
#  ifdef OBC_M2FLATHER
        do i=istr,iend
          cff=0.5*(h(i,jstr-1)+h(i,jstr))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5*(pn(i,jstr-1)+pn(i,jstr))

          zx=(0.5+cx)*zeta(i,jstr,kstp)+(0.5-cx)*zeta(i,jstr-1,kstp)
          if (cx > 0.292893218813452) then
            zx=zx + ( zeta(i,jstr,knew) +cx*zeta(i,jstr-1,kstp)
     &                               -(1.+cx)*zeta(i,jstr,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
          endif

          vbar(i,jstr,knew)=0.5*( (1.-cx)*vbar(i,jstr,kstp)
     &                               +cx*vbar(i,jstr+1,kstp)
#   ifdef M2_FRC_BRY
     &                 +vbar_south(i) -hx*(zx-zeta_south(i))
#   else
     &                 +vbclm(i,jstr) -hx*(zx-ssh(i,jstr-1))
#   endif
     &                                                    )
#   ifdef MASKING
     &                                      *vmask(i,jstr)
#   endif
        enddo

#  elif defined OBC_M2ORLANSKI
        do i=istr,iend+1
          grad(i,jstr  )=(vbar(i,jstr  ,kstp)-vbar(i-1,jstr  ,kstp))
#   ifdef MASKING
     &                                                *pmask(i,jstr)
#   endif
          grad(i,jstr+1)=(vbar(i,jstr+1,kstp)-vbar(i-1,jstr+1,kstp))
#   ifdef MASKING
     &                                              *pmask(i,jstr+1)
#   endif
        enddo
        do i=istr,iend
          dft=vbar(i,jstr+1,kstp)-vbar(i,jstr+1,knew)
          dfx=vbar(i,jstr+1,knew)-vbar(i,jstr+2,knew)

#   ifdef OBC_RAD_NORMAL
          cy=0. ; cff=max(dfx*dfx, eps)
#   else
          if (dft*(grad(i,jstr+1)+grad(i+1,jstr+1)) > 0.) then
            dfy=grad(i  ,jstr+1)
          else
            dfy=grad(i+1,jstr+1)
          endif
          cff=max(dfx*dfx+dfy*dfy, eps)
#    ifdef OBC_RAD_NPO
          cy=0.
#    else
          cy=min(cff,max(dft*dfy,-cff))
#    endif
#   endif

          cx=dft*dfx
          if (cx < 0.) then       !--> cancel, if inflow
            cx=0. ; cy=0.
#   if defined M2_FRC_BRY || defined M2NUDGING
            cext=dtfast*ubind * 0.5*(pn(i,jstr-1)+pn(i,jstr))
          else
            cext=0.
#   endif
          endif

          vbar(i,jstr,knew)=( cff*vbar(i,jstr,kstp)
     &                       +cx*vbar(i,jstr+1,knew)
     &                     -max(cy,0.)*grad(i  ,jstr)
     &                     -min(cy,0.)*grad(i+1,jstr)
     &                                    )/(cff+cx)
#   ifdef OBC_FLUX_CORR
     &                                    -uv_crss
#   endif
#   if defined M2_FRC_BRY || defined M2NUDGING
          vbar(i,jstr,knew)=(1.-cext)*vbar(i,jstr,knew)
#    ifdef M2_FRC_BRY
     &                                     +cext*vbar_south(i)
#    else
     &                                     +cext*vbclm(i,jstr)
#    endif
     &          -attnM2*dtfast*g*0.5*(pn(i,jstr-1)+pn(i,jstr))
#    ifdef M2_FRC_BRY
     &                      *(zeta(i,jstr,knew)-zeta_south(i))
#    else
     &                      *(zeta(i,jstr,knew)-ssh(i,jstr-1))
#    endif
#   endif
#   ifdef MASKING
          vbar(i,jstr,knew)=vbar(i,jstr,knew)*vmask(i,jstr)
#   endif
        enddo

#  elif defined OBC_M2SPECIFIED
        do i=istr,iend
#   ifdef M2_FRC_BRY
          vbar(i,jstr,knew)=vbar_south(i)
#   else
          vbar(i,jstr,knew)=vbclm(i,jstr)
#   endif
#   ifdef MASKING
     &                     *vmask(i,jstr)
#   endif
        enddo
#  else
        do i=istr,iend
          vbar(i,jstr,knew)=vbar(i,jstr+1,knew)
#   ifdef MASKING
     &                           *vmask(i,jstr)
#   endif
        enddo
#  endif
# else
                                             ! Southern edge closed
        do i=istr,iend                       ! ======== ==== ======
#  ifdef WEC
          vbar(i,jstr,knew)=-vst2d(i,jstr)   ! no Lagrangian flux
#   ifdef MASKING
     &                      *vmask(i,jstr)
#   endif
#  else
          vbar(i,jstr,knew)=0.               ! (no-flux, default)
#  endif

        enddo
# endif          /* OBC_SOUTH */
      endif   !<-- SOUTHERN_EDGE



      if (NORTHERN_EDGE) then
# ifdef OBC_NORTH
#  if defined OBC_M2FLATHER
        do i=istr,iend
          cff=0.5*(h(i,jend)+h(i,jend+1))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5*(pn(i,jend)+pn(i,jend+1))

          zx=(0.5+cx)*zeta(i,jend,kstp)+(0.5-cx)*zeta(i,jend+1,kstp)
          if (cx > 0.292893218813452) then
            zx=zx + ( zeta(i,jend,knew) +cx*zeta(i,jend+1,kstp)
     &                               -(1.+cx)*zeta(i,jend,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
          endif

          vbar(i,jend+1,knew)=0.5*( (1.-cx)*vbar(i,jend+1,kstp)
     &                                     +cx*vbar(i,jend,kstp)
#   ifdef M2_FRC_BRY
     &                     +vbar_north(i) +hx*(zx-zeta_north(i))
#   else
     &                   +vbclm(i,jend+1) +hx*(zx-ssh(i,jend+1))
#   endif
     &                                                        )
#   ifdef MASKING
     &                                        *vmask(i,jend+1)
#   endif
        enddo

#  elif defined OBC_M2ORLANSKI
        do i=istr,iend+1
          grad(i,jend  )=(vbar(i,jend  ,kstp)-vbar(i-1,jend  ,kstp))
#   ifdef MASKING
     &                                                *pmask(i,jend)
#   endif
          grad(i,jend+1)=(vbar(i,jend+1,kstp)-vbar(i-1,jend+1,kstp))
#   ifdef MASKING
     &                                              *pmask(i,jend+1)
#   endif
        enddo
        do i=istr,iend
          dft=vbar(i,jend,kstp)-vbar(i,jend  ,knew)
          dfx=vbar(i,jend,knew)-vbar(i,jend-1,knew)

#   ifdef OBC_RAD_NORMAL
          cy=0. ; cff=max(dfx*dfx, eps)
#   else
          if (dft*(grad(i,jend)+grad(i+1,jend)) > 0.) then
            dfy=grad(i  ,jend)
          else
            dfy=grad(i+1,jend)
          endif
          cff=max(dfx*dfx+dfy*dfy, eps)
#    ifdef OBC_RAD_NPO
          cy=0.
#    else
          cy=min(cff,max(dft*dfy,-cff))
#    endif
#   endif

          cx=dft*dfx
          if (cx < 0.) then       !--> cancel, if inflow
            cx=0. ; cy=0.
#   if defined M2_FRC_BRY || defined M2NUDGING
            cext=dtfast*ubind * 0.5*(pn(i,jend)+pn(i,jend+1))
          else
            cext=0.
#   endif
          endif

          vbar(i,jend+1,knew)=( cff*vbar(i,jend+1,kstp)
     &                             +cx*vbar(i,jend,knew)
     &                      -max(cy,0.)*grad(i  ,jend+1)
     &                      -min(cy,0.)*grad(i+1,jend+1)
     &                                      )/(cff+cx)
#   ifdef OBC_FLUX_CORR
     &                                      +uv_crss
#   endif
#   if defined M2_FRC_BRY || defined M2NUDGING
          vbar(i,jend+1,knew)=(1.-cext)*vbar(i,jend+1,knew)
#    ifdef M2_FRC_BRY
     &                                    +cext*vbar_north(i)
#    else
     &                                   +cext*vbclm(i,jend+1)
#    endif
     &          +attnM2*dtfast*g*0.5*(pn(i,jend)+pn(i,jend+1))
#    ifdef M2_FRC_BRY
     &                      *(zeta(i,jend,knew)-zeta_north(i))
#    else
     &                      *(zeta(i,jend,knew)-ssh(i,jend+1))
#    endif
#   endif
#   ifdef MASKING
          vbar(i,jend+1,knew)=vbar(i,jend+1,knew)*vmask(i,jend+1)
#   endif
        enddo
#  elif defined OBC_M2SPECIFIED
        do i=istr,iend
#   ifdef M2_FRC_BRY
          vbar(i,jend+1,knew)=vbar_north(i)
#   else
          vbar(i,jend+1,knew)=vbclm(i,jend+1)
#   endif
#   ifdef MASKING
     &                       *vmask(i,jend+1)
#   endif
        enddo
#  else
        do i=istr,iend
          vbar(i,jend+1,knew)=vbar(i,jend,knew)
#   ifdef MASKING
     &                         *vmask(i,jend+1)
#   endif
        enddo
#  endif
# else
        do i=istr,iend                        ! Northern edge closed
#  ifdef WEC
          vbar(i,jend+1,knew)=-vst2d(i,jend+1)  ! no Lagrangian flux
#   ifdef MASKING
     &                        *vmask(i,jend+1)
#   endif
#  else
          vbar(i,jend+1,knew)=0.                ! (no-flux, default)
#  endif
        enddo
# endif /* OBC_NORTH */
      endif          !<--  NORTHERN_EDGE
#endif                  /* !NS_PERIODIC */




c****
c****
#ifdef OBC_M2FLATHER
# undef OBC_M2FLATHER
# define OBC_M2ORLANSKI
#endif

#ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
# ifdef OBC_WEST
#  if defined OBC_M2ORLANSKI
        do j=jstrV-1,jend
          grad(istr-1,j)=vbar(istr-1,j+1,kstp)-vbar(istr-1,j,kstp)
          grad(istr  ,j)=vbar(istr  ,j+1,kstp)-vbar(istr  ,j,kstp)
        enddo
        do j=jstrV,jend
#   ifdef LEGACY_ORLANSKI
          dft=vbar(istr,j,kstp)-vbar(istr  ,j,knew)
          dfx=vbar(istr,j,knew)-vbar(istr+1,j,knew)

#    ifdef OBC_RAD_NORMAL
          cy=0 ; cff=max(dfx*dfx, eps)
#    else
          if (dft*(grad(istr,j-1)+grad(istr,j)) > 0.) then
            dfy=grad(istr,j-1)
          else
            dfy=grad(istr,j  )
          endif
          cff=max(dfx*dfx+dfy*dfy, eps)
#     ifdef OBC_RAD_NPO
          cy=0.
#     else
          cy=min(cff,max(dft*dfy,-cff))
#     endif
#    endif

          cx=dft*dfx
          if (cx < 0.) then       !--> cancel, if inflow
            cx=0. ; cy=0.
#    if defined M2_FRC_BRY || defined M2NUDGING
            cext=dtfast*ubind * 0.5*(pm(istr-1,j-1)+pm(istr-1,j))
          else
            cext=0.
#    endif
          endif

          vbar(istr-1,j,knew)=( cff*vbar(istr-1,j,kstp)
     &                             +cx*vbar(istr,j,knew)
     &                      -max(cy,0.)*grad(istr-1,j-1)
     &                      -min(cy,0.)*grad(istr-1,j  )
     &                                       )/(cff+cx)

#   else

          cx=-0.125*dtfast*(ubar(istr,j,kstp)+ubar(istr,j-1,kstp))
     &      *(pm(istr-1,j)+pm(istr-1,j-1)+pm(istr,j)+pm(istr,j-1))

          cy= 0.125*dtfast*(vbar(istr-1,j,kstp)+vbar(istr,j,kstp))
     &      *(pn(istr-1,j)+pn(istr-1,j-1)+pn(istr,j)+pn(istr,j-1))

          if (cx > 0.) then
            cext=0.
          else
            cext=-cx ; cx=0.
          endif

          vbar(istr-1,j,knew)=(1.-cx)*( vbar(istr-1,j,kstp)
     &                          -max(cy,0.)*grad(istr-1,j-1)
     &                          -min(cy,0.)*grad(istr-1,j  )
     &                                                     )
     &                       +cx*(        vbar(istr,j,kstp)
     &                            -max(cy,0.)*grad(istr,j-1)
     &                            -min(cy,0.)*grad(istr,j  )
     &                                                     )

#   endif /* LEGACY_ORLANSKI */

#   if defined M2_FRC_BRY || defined M2NUDGING
          vbar(istr-1,j,knew)=(1.-cext)*vbar(istr-1,j,knew)
#    ifdef M2_FRC_BRY
     &                                   +cext*vbar_west(j)
#    else
     &                                +cext*vbclm(istr-1,j)
#    endif
#   endif
#   ifdef MASKING
          vbar(istr-1,j,knew)=vbar(istr-1,j,knew)*vmask(istr-1,j)
#   endif
        enddo
#  else          /* alternative open: Chapman, specified, gradient */
        do j=jstrV,jend
#   if defined OBC_M2FLATHER
          cx=dtfast*sqrt(0.5*g*(h(istr-1,j-1)+h(istr-1,j)))
     &                    *0.5*(pm(istr-1,j-1)+pm(istr-1,j))
          vbar(istr-1,j,knew)=( vbar(istr-1,j,kstp)
     &               +cx*vbar(istr,j,knew) )/(1.+cx)
#   elif defined OBC_M2SPECIFIED
#    ifdef M2_FRC_BRY
          vbar(istr-1,j,knew)=vbar_west(j)
#    else
          vbar(istr-1,j,knew)=vbclm(istr-1,j)
#    endif
#   else
          vbar(istr-1,j,knew)=vbar(istr,j,knew)
#   endif
#   ifdef MASKING
          vbar(istr-1,j,knew)=vbar(istr-1,j,knew)*vmask(istr-1,j)
#   endif
        enddo
#  endif
# else
#  ifdef NS_PERIODIC
#   define J_RANGE jstrV,jend
#  else
#   define J_RANGE jstr,jendR
#  endif
                                   ! Closed BC: free-slip (gamma2=+1)
        do j=J_RANGE               ! ====== ===   no-slip (gamma2=-1)
          vbar(istr-1,j,knew)=gamma2*vbar(istr,j,knew)
#  ifdef MASKING
     &                                *vmask(istr-1,j)
#  endif
        enddo
#  undef J_RANGE
# endif              /* OBC_WEST */
      endif       !<-- WESTERN_EDGE


c--#ifdef ISWAKE_XXX
c--# define OBC_M2ORLANSKI
c--#endif



      if (EASTERN_EDGE) then
# ifdef OBC_EAST
#  if defined OBC_M2ORLANSKI
        do j=jstrV-1,jend
          grad(iend  ,j)=vbar(iend  ,j+1,kstp)-vbar(iend  ,j,kstp)
          grad(iend+1,j)=vbar(iend+1,j+1,kstp)-vbar(iend+1,j,kstp)
        enddo
        do j=jstrV,jend
#   ifdef LEGACY_ORLANSKI
          dft=vbar(iend,j,kstp)-vbar(iend  ,j,knew)
          dfx=vbar(iend,j,knew)-vbar(iend-1,j,knew)

#    ifdef OBC_RAD_NORMAL
          cy=0. ; cff=max(dfx*dfx, eps)
#    else
          if (dft*(grad(iend,j-1)+grad(iend,j)) > 0.) then
            dfy=grad(iend,j-1)
          else
            dfy=grad(iend,j  )
          endif
          cff=max(dfx*dfx+dfy*dfy, eps)
#     ifdef OBC_RAD_NPO
          cy=0.
#     else
          cy=min(cff,max(dft*dfy,-cff))
#     endif
#    endif

          cx=dft*dfx
          if (cx < 0.) then       !--> cancel, if inflow
            cx=0. ; cy=0.
#    if defined M2_FRC_BRY || defined M2NUDGING
            cext=dtfast*ubind * 0.5*(pm(iend+1,j-1)+pm(iend+1,j))
          else
            cext=0.
#    endif
          endif

          vbar(iend+1,j,knew)=( cff*vbar(iend+1,j,kstp)
     &                             +cx*vbar(iend,j,knew)
     &                       -max(cy,0.)*grad(iend+1,j-1)
     &                       -min(cy,0.)*grad(iend+1,j  )
     &                                        )/(cff+cx)

#   else

          cx=0.125*dtfast*(ubar(iend+1,j,kstp)+ubar(iend+1,j-1,kstp))
     &         *(pm(iend,j)+pm(iend,j-1)+pm(iend+1,j)+pm(iend+1,j-1))

          cy=0.125*dtfast*(vbar(iend,j,kstp)+vbar(iend+1,j,kstp))
     &         *(pn(iend,j)+pn(iend,j-1)+pn(iend+1,j)+pn(iend+1,j-1))

          if (cx > 0.) then
            cext=0.
          else
            cext=-cx ; cx=0.
          endif

          vbar(iend+1,j,knew)=(1.-cx)*( vbar(iend+1,j,kstp)
     &                          -max(cy,0.)*grad(iend+1,j-1)
     &                          -min(cy,0.)*grad(iend+1,j  )
     &                                                     )
     &                       +cx*(        vbar(iend,j,kstp)
     &                            -max(cy,0.)*grad(iend,j-1)
     &                            -min(cy,0.)*grad(iend,j  )
     &                                                     )

#   endif /* LEGACY_ORLANSKI */

#   if defined M2_FRC_BRY || defined M2NUDGING
          vbar(iend+1,j,knew)=(1.-cext)*vbar(iend+1,j,knew)
#    ifdef M2_FRC_BRY
     &                                   +cext*vbar_east(j)
#    else
     &                                +cext*vbclm(iend+1,j)
#    endif
#   endif
#   ifdef MASKING
          vbar(iend+1,j,knew)=vbar(iend+1,j,knew)*vmask(iend+1,j)
#   endif
        enddo
#  else          /* alternative open: Chapman, specified, gradient */
        do j=jstrV,jend
#   if defined OBC_M2FLATHER
          cx=dtfast*sqrt(0.5*g*(h(iend+1,j-1)+h(iend+1,j)))
     &                    *0.5*(pm(iend+1,j-1)+pm(iend+1,j))
          vbar(iend+1,j,knew)=( vbar(iend+1,j,kstp)
     &                +cx*vbar(iend,j,knew))/(1.+cx)
#   elif defined OBC_M2SPECIFIED
#    ifdef M2_FRC_BRY
          vbar(iend+1,j,knew)=vbar_east(j)
#    else
          vbar(iend+1,j,knew)=vbclm(iend+1,j)
#    endif
#   else
          vbar(iend+1,j,knew)=vbar(iend,j,knew)
#   endif
#   ifdef MASKING
          vbar(iend+1,j,knew)=vbar(iend+1,j,knew)*vmask(iend+1,j)
#   endif
        enddo
#  endif
# else
#  ifdef NS_PERIODIC
#   define J_RANGE jstrV,jend
#  else
#   define J_RANGE jstr,jendR
#  endif
                                        ! Wall: free-slip (gamma2=+1)
        do j=J_RANGE                    ! =====   no-slip (gamma2=-1)
          vbar(iend+1,j,knew)=gamma2*vbar(iend,j,knew)
#  ifdef MASKING
     &                                *vmask(iend+1,j)
#  endif
        enddo
#  undef J_RANGE
# endif            /* OBC_EAST */
      endif      !<-- EASTERN_EDGE
#endif             /* !EW_PERIODIC */

                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========

#if defined OBC_SOUTH && defined OBC_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        vbar(istr-1,jstr,knew)=0.5*( vbar(istr-1,jstr+1,knew)
     &                              +vbar(istr  ,jstr  ,knew))
      endif
#endif
#if defined OBC_SOUTH && defined OBC_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        vbar(iend+1,jstr,knew)=0.5*( vbar(iend+1,jstr+1,knew)
     &                              +vbar(iend  ,jstr  ,knew))
      endif
#endif
#if defined OBC_NORTH && defined OBC_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        vbar(istr-1,jend+1,knew)=0.5*( vbar(istr-1,jend,knew)
     &                                +vbar(istr,jend+1,knew))
      endif
#endif
#if defined OBC_NORTH && defined OBC_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        vbar(iend+1,jend+1,knew)=0.5*( vbar(iend+1,jend,knew)
     &                                +vbar(iend,jend+1,knew))
      endif
#endif

! Daniel add
!====================================================================
!       NO LAGRANGIAN FLUX CONDITION AT LAND-SEA INTERFACES
!====================================================================
#if defined WEC && defined MASKING
      do i=istrR,iendR
        do j=jstr,jendR
          vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
     &                  +vst2d(i,j)*(vmask(i,j)-1.0)
        enddo
      enddo
#endif

      end

