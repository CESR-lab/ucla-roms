#include "cppdefs.h"
#if defined SOLVE3D && defined LMD_KPP
c--#define BULK_DIFFERENCE
#define INT_AT_RHO_POINTS
c--#define INT_AT_W_POINTS
c--#define SMOOTH_HBL
c--#define LIMIT_UNSTABLE_ONLY


      subroutine lmd_kpp_tile (istr,iend,jstr,jend, Kv,Kt,Ks, Rib,
     &                         ustar, hbl,srf_dpth, Bo,Bosol, Bfsfc,
     &                         Gm1,dGm1dS,  Gt1,dGt1dS,  Gs1,dGs1dS,
     &                                           zscale,ws,wm,
     &                                Hz_inv,du,dv,bvr, kbl)
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# define tind nstp

      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k, itrc
      real Kv(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &     Kt(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &     Ks(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &    Rib(PRIVATE_2D_SCRATCH_ARRAY,0:N),

     &      ustar(PRIVATE_2D_SCRATCH_ARRAY),
     &         Bo(PRIVATE_2D_SCRATCH_ARRAY),
     &      Bosol(PRIVATE_2D_SCRATCH_ARRAY),
     &      Bfsfc(PRIVATE_2D_SCRATCH_ARRAY),

     &        hbl(PRIVATE_2D_SCRATCH_ARRAY),
     &   srf_dpth(PRIVATE_2D_SCRATCH_ARRAY),

     &        Gm1(PRIVATE_2D_SCRATCH_ARRAY),
     &     dGm1dS(PRIVATE_2D_SCRATCH_ARRAY),
     &        Gt1(PRIVATE_2D_SCRATCH_ARRAY),
     &     dGt1dS(PRIVATE_2D_SCRATCH_ARRAY),
     &        Gs1(PRIVATE_2D_SCRATCH_ARRAY),
     &     dGs1dS(PRIVATE_2D_SCRATCH_ARRAY),

     &     zscale(PRIVATE_2D_SCRATCH_ARRAY),
     &         ws(PRIVATE_2D_SCRATCH_ARRAY),
     &         wm(PRIVATE_2D_SCRATCH_ARRAY)

      integer kmo(PRIVATE_2D_SCRATCH_ARRAY)
      integer kbl(PRIVATE_2D_SCRATCH_ARRAY), km,kp
      real h_MO(PRIVATE_2D_SCRATCH_ARRAY)
      real Bfsfc_up(PRIVATE_2D_SCRATCH_ARRAY), z_up

      real Hz_inv(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &         du(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &         dv(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &        bvr(PRIVATE_1D_SCRATCH_ARRAY,0:N)

      real hekman, Vtc, dVsq,   Vtsq,   f1,
     &     cff,    z_bl,   a1,  Av_bl,  At_bl,  As_bl,
     &     cff_up, sigma,  a2,  dAv_bl, dAt_bl, dAs_bl,
     &     cff_dn, Ritop,  a3,  Gm,     Gt,     Gs



      real nubl, nu0c, cs, Cv, Ricr, betaT, epssfc,
     &             C_Ek, C_MO, Cstar, Cg, eps
      parameter (
     &   nubl=0.01,     ! Maximum allowed boundary layer

     &   nu0c=0.1,      ! convective adjustment
                        ! viscosity and diffusivity [m^2/s].

     &   cs=98.96,      ! see parameter associated with turbulent
                        !     velocity scales in lmd_wscale.F

     &   Cv=1.8,        ! Ratio of interior Brunt-Vaisala frequency
                        ! "N" to that at the entrainment depth "he".

     &   Ricr=0.45,     ! Critical bulk Richardson number (0.3).

     &   betaT=-0.2,    ! Ratio of entrainment flux to
                        !       to surface buoyancy flux.

     &   epssfc=0.1,    ! Nondimensional extent of the surface layer.

     &   C_Ek=0.7,      ! Constant for computation of Ekman depth.

     &   C_MO=1.,       ! Constant for computaion Monin-Obukhov depth.

     &   Cstar=10.,     ! Proportionality coefficient parameterizing
                        ! nonlocal transport.
     &   eps=1.E-20)

      Cg=Cstar * vonKar * (cs*vonKar*epssfc)**(1./3.)
      Vtc=Cv * sqrt(-betaT/(cs*epssfc)) / (Ricr*vonKar**2)

# ifdef SMOOTH_HBL
#  ifdef EW_PERIODIC
      imin=istr-1
      imax=iend+1
#  else
      if (WESTERN_EDGE) then
        imin=istr
      else
        imin=istr-1
      endif
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
#  endif
#  ifdef NS_PERIODIC
      jmin=jstr-1
      jmax=jend+1
#  else
      if (SOUTHERN_EDGE) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=jend
      else
        jmax=jend+1
      endif
#  endif
# else
      imin=istr
      imax=iend
      jmin=jstr
      jmax=jend
# endif
# define I_EXT_RANGE imin,imax
# define J_EXT_RANGE jmin,jmax


! Compute thermal expansion coefficient "alpha" [kg/m^3/decC] and
! saline contraction coefficient "beta" [kg/m^3/PSU] at surface, then
! compute surface turbulent buoyancy forcing "Bo" [m^2/s^3] (in doing
! so remove incoming solar shortwave radiation component and save it
! separately as "Bosol"). Also get an approximation for ssurface
! layer depth using "epssfc" and boundary layer depth from previous
! time step (this is needed to estimate turbulent velocity scale
! in computation of Vterm in Rib, before new hbl is found). Also
! compute turbulent  friction velocity "ustar" from wind stress at
! RHO-points. Finally, initialize boundary layer depth "hbl" and
! index "kbl" to the maximum (bottomed out) values.

#define alpha wm
#define beta  ws
      call alfabeta_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                             J_EXT_RANGE, alpha,beta)
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          Bo(i,j)=g*( alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
# ifdef SALINITY
     &                              -beta(i,j)*stflx(i,j,isalt)
# endif
     &                                                        )
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
          Bfsfc(i,j)=Bo(i,j)  !<-- at surface

          srf_dpth(i,j)=epssfc*hbls(i,j,tind)

          ustar(i,j)=sqrt(0.5*sqrt( (sustr(i,j)+sustr(i+1,j))**2
     &                             +(svstr(i,j)+svstr(i,j+1))**2))

          hbl(i,j)=z_w(i,j,N)-z_w(i,j,0) +eps
          kbl(i,j)=0
          Rib(i,j,N)=0.
          Rib(i,j,0)=0.
          h_MO(i,j)=C_MO*ustar(i,j)**3/vonKar
          kmo(i,j)=0
        enddo                      !--> discard alpha,beta;
      enddo                        !    keep Bo,Bosol to the very end.
#undef beta
#undef alpha

! Compute fraction of the solar shortwave flux "swdk" penetrating
! to grid level depth (at vertical RHO-points). After which compute
! total surface buoyancy flux "Bfsfc" as the sum of contributions
! from surface turbulent buoyancy forcing "Bo" and radiative flux
! down to boundary layer depth [Bosol*(1-swdk)]. After that compute
! turbulent velocity scales "wm,ws"
!
! Compute bulk Richardson number "Rib" at every grid level, and find
! index "k" when it becomes > Ricr for the first time.  Then find
! depth of the oceanic planetary boundary layer "hbl", such that
! Rib(hbl)=Ric by linear interpolation using nearest points.
!
!                  [Br - B(d)] * d
!     Rib(d) = ----------------------- ;       Rib(hbl)=Ricr
!              |Vr - V(d)|^2 + Vt(d)^2
!
! To do so, first compute numerator of bluk Richardson number,
! Ritop=(Br-B)*d, where Br is the near-surface reference buoyancy,
! B is the mean buoyancy as function of d, and d is the distance
! coordinate from the boundary.  Then compute the square of
! velocity shear relative to reference velocities, dVsq=|Vr-V|^2,
! at horizontal and vertical RHO-points.

# if defined BULK_DIFFERENCE
#  define swdk wm
      do k=N-1,1,-1
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            zscale(i,j)=z_w(i,j,N)-z_r(i,j,k)
          enddo
        enddo
        call lmd_swfrac_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                                J_EXT_RANGE, zscale,swdk)
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            Bfsfc(i,j)=Bo(i,j)+Bosol(i,j)*(1.-swdk(i,j))
#  undef swdk
#  ifdef LIMIT_UNSTABLE_ONLY
            if (Bfsfc(i,j)<0.) then
              zscale(i,j)=min(zscale(i,j), srf_dpth(i,j))
            endif
#  else
            zscale(i,j)=min(zscale(i,j), srf_dpth(i,j))
#  endif
          enddo
        enddo

        call lmd_wscale_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                J_EXT_RANGE, ustar, Bfsfc, zscale, wm,ws)
        cff=g/rho0
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            Ritop=-cff*(rho1(i,j,N)-rho1(i,j,k))
     &                   *(z_r(i,j,N)-z_r(i,j,k))

            dVsq=0.25*( ( u(i  ,j,N,tind)-u(i  ,j,k,tind)
     &                   +u(i+1,j,N,tind)-u(i+1,j,k,tind) )**2
     &                 +( v(i,j  ,N,tind)-v(i,j  ,k,tind)
     &                   +v(i,j+1,N,tind)-v(i,j+1,k,tind) )**2)

            Vtsq=Vtc*ws(i,j)*(z_r(i,j,N)-z_r(i,j,k))
     &         *sqrt(max(0., 0.5*(bvf(i,j,k)+bvf(i,j,k-1))))

            Rib(i,j,k)= Ritop/(dVsq+Vtsq+eps)

            if (kbl(i,j)==0 .and. Rib(i,j,k)>Ricr) kbl(i,j)=k
          enddo
        enddo
      enddo

      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          if (kbl(i,j)>0) then
            k=kbl(i,j)
            hbl(i,j)=z_w(i,j,N)-( z_r(i,j,k)*(Rib(i,j,k+1)-Ricr)
     &                            +z_r(i,j,k+1)*(Ricr-Rib(i,j,k))
     &                                )/(Rib(i,j,k+1)-Rib(i,j,k))
          endif
#  ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
#  endif
        enddo     !--> discard kbl
      enddo
# elif defined INT_AT_RHO_POINTS

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

      cff=1./Ricr
      do k=N-1,1,-1
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            Rib(i,j,k)=Rib(i,j,k+1)

     &               +0.5*( ( u(i  ,j,k+1,tind)-u(i  ,j,k,tind)
     &                       +u(i+1,j,k+1,tind)-u(i+1,j,k,tind) )**2
     &                     +( v(i,j  ,k+1,tind)-v(i,j  ,k,tind)
     &                       +v(i,j+1,k+1,tind)-v(i,j+1,k,tind) )**2
     &                    )/(Hz(i,j,k)+Hz(i,j,k+1))

     &             -0.5*(Hz(i,j,k)+Hz(i,j,k+1))*( cff*bvf(i,j,k)
     &                                      + 400.*f(i,j)*f(i,j) )
          enddo
        enddo
      enddo

#  define swdk wm
      do k=N-1,1,-1
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            zscale(i,j)=z_w(i,j,N)-z_r(i,j,k)
          enddo
        enddo
        call lmd_swfrac_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                                J_EXT_RANGE, zscale,swdk)
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            Bfsfc_up(i,j)=Bfsfc(i,j)
            Bfsfc(i,j)=Bo(i,j)+Bosol(i,j)*(1.-swdk(i,j))
#  undef swdk
#  ifdef LIMIT_UNSTABLE_ONLY
            if (Bfsfc(i,j)<0.) then
              zscale(i,j)=min(zscale(i,j), srf_dpth(i,j))
            endif
#  else
            zscale(i,j)=min(zscale(i,j), srf_dpth(i,j))
#  endif
            if (Bfsfc(i,j)*(z_w(i,j,N)-z_r(i,j,k)) > h_MO(i,j)
     &                                  .and. kmo(i,j)==0) then
               kmo(i,j)=k

               cff=h_MO(i,j)

               cff_dn=max(0., Bfsfc(i,j))
               cff_up=max(0., Bfsfc_up(i,j))
               if (k==N) then
                 z_up=z_w(i,j,N)
               else
                 z_up=z_r(i,j,k+1)
               endif

            h_MO(i,j)=z_w(i,j,N) +cff*(cff_up*z_up-cff_dn*z_r(i,j,k))
     &      /( cff_up*cff_dn*(z_up-z_r(i,j,k)) + cff*(cff_dn-cff_up) )

            endif
          enddo
        enddo


        call lmd_wscale_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                J_EXT_RANGE, ustar, Bfsfc, zscale, wm,ws)
        cff=g/rho0
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            Vtsq=Vtc*ws(i,j) * sqrt( max(0., 0.5*(bvf(i,j,k)
     &                                       +bvf(i,j,k-1))))
            Rib(i,j,k)=Rib(i,j,k)+Vtsq

            if (kbl(i,j)==0 .and. Rib(i,j,k)<0.) kbl(i,j)=k
          enddo
        enddo
      enddo   !<-- k

      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          if (kbl(i,j)>0) then
            k=kbl(i,j)
            hbl(i,j)=z_w(i,j,N)-( z_r(i,j,k)*Rib(i,j,k+1)
     &                            -z_r(i,j,k+1)*Rib(i,j,k)
     &                         )/(Rib(i,j,k+1)-Rib(i,j,k))
          endif

          if (kmo(i,j)>0) hbl(i,j)=min(hbl(i,j),h_MO(i,j))


#  ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
#  endif
        enddo     !--> discard kbl
      enddo

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

# elif defined INT_AT_W_POINTS

      do j=J_EXT_RANGE
        do k=1,N-1
          do i=I_EXT_RANGE
            Hz_inv(i,k)=1.D0/(Hz(i,j,k+1)+Hz(i,j,k))
            du(i,k)=Hz_inv(i,k)*(u(i  ,j,k+1,tind)-u(i  ,j,k,tind)
     &                          +u(i+1,j,k+1,tind)-u(i+1,j,k,tind))
            dv(i,k)=Hz_inv(i,k)*(v(i,j  ,k+1,tind)-v(i,j  ,k,tind)
     &                          +v(i,j+1,k+1,tind)-v(i,j+1,k,tind))
          enddo
        enddo

        do k=2,N-1
          do i=I_EXT_RANGE
            cff=(bvf(i,j,k)*Hz_inv(i,k)+bvf(i,j,k-1)*Hz_inv(i,k-1))
     &                                  /(Hz_inv(i,k)+Hz_inv(i,k-1))
            if (bvf(i,j,k)>0.D0  .and.  bvf(i,j,k-1)>0.D0) then
              bvr(i,k)=min(cff, 2.D0*min(bvf(i,j,k),bvf(i,j,k-1)))
            elseif (bvf(i,j,k)<0.D0.and.bvf(i,j,k-1)<0.D0) then
              bvr(i,k)=max(cff, 2.D0*max(bvf(i,j,k),bvf(i,j,k-1)))
            else
              bvr(i,k)=0.D0
            endif
          enddo
        enddo
        do i=I_EXT_RANGE
          bvr(i,N)=bvf(i,j,N-1)
          bvr(i,1)=bvf(i,j,  1)
           du(i,N)=du(i,N-1)
           du(i,0)=du(i,1)
           dv(i,N)=dv(i,N-1)
           dv(i,0)=dv(i,1)
        enddo

        cff=1./Ricr

        do k=N,1,-1
          kp=min(k+1,N)
          km=max(k-1,1)
          do i=I_EXT_RANGE
            Rib(i,j,k-1)=Rib(i,j,k) + Hz(i,j,k)*(  0.375D0*(
     &            du(i,k)**2+du(i,k-1)**2 +dv(i,k)**2+dv(i,k-1)**2
     &          ) + 0.25D0*( du(i,k)*du(i,k-1)+dv(i,k)*dv(i,k-1) ))

     &    -cff*( 0.5*Hz(i,j,k)*(bvf(i,j,k)+bvf(i,j,k-1)) -0.0625*(
     &                  (Hz(i,j,kp)+Hz(i,j,k))*(bvr(i,kp)-bvr(i,k))
     &                 -(Hz(i,j,k)+Hz(i,j,km))*(bvr(i,k)-bvr(i,km))
     &                                                           ))
          enddo
        enddo

c       cff=1./Ricr
c       do k=2,N-1
c         do i=I_EXT_RANGE
c            Rib(i,j,k)=Rib(i,j,k+1)
c     &          + 0.5*( ( u(i  ,j,k+1,tind)-u(i  ,j,k,tind)
c     &                   +u(i+1,j,k+1,tind)-u(i+1,j,k,tind) )**2
c     &                 +( v(i,j  ,k+1,tind)-v(i,j  ,k,tind)
c     &                   +v(i,j+1,k+1,tind)-v(i,j+1,k,tind) )**2
c     &                    )/(Hz(i,j,k)+Hz(i,j,k+1))
c     &          - cff*bvf(i,j,k) * 0.5*(Hz(i,j,k)+Hz(i,j,k+1))
c          enddo
c        enddo
c        do i=I_EXT_RANGE
c          Rib(i,j,0)=Rib(i,j,1)
c        enddo
c        do k=1,N-1,+1
c          do i=I_EXT_RANGE
c            Rib(i,j,k)=( Hz(i,j,k+1)*Rib(i,j,k)
c     &                   +Hz(i,j,k)*Rib(i,j,k+1)
c     &                 )/(Hz(i,j,k)+Hz(i,j,k+1))
c          enddo
c        enddo
c        do i=I_EXT_RANGE
c          Rib(i,j,0)=2.*Rib(i,j,0)-Rib(i,j,1)
c        enddo
       enddo    !<-- k

#  define swdk wm
      do k=N,1,-1
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            zscale(i,j)=z_w(i,j,N)-z_w(i,j,k-1)
          enddo
        enddo
        call lmd_swfrac_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                                J_EXT_RANGE, zscale,swdk)
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            Bfsfc(i,j)=Bo(i,j)+Bosol(i,j)*(1.-swdk(i,j))
#  undef swdk
#  ifdef LIMIT_UNSTABLE_ONLY
            if (Bfsfc(i,j)<0.) then
              zscale(i,j)=min(zscale(i,j), srf_dpth(i,j))
            endif
#  else
            zscale(i,j)=min(zscale(i,j), srf_dpth(i,j))
#  endif
          enddo
        enddo

        call lmd_wscale_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                J_EXT_RANGE, ustar, Bfsfc, zscale, wm,ws)
        cff=g/rho0
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            Vtsq=Vtc*ws(i,j) * sqrt( max(0., bvf(i,j,k-1) ))

            Rib(i,j,k-1)=Rib(i,j,k-1)+Vtsq

            if (kbl(i,j)==0 .and. Rib(i,j,k-1)<0.) kbl(i,j)=k
          enddo
        enddo
      enddo
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          if (kbl(i,j)>0) then
            k=kbl(i,j)
            hbl(i,j)=z_w(i,j,N)-( z_w(i,j,k-1)*Rib(i,j,k)
     &                            -z_w(i,j,k)*Rib(i,j,k-1)
     &                         )/(Rib(i,j,k)-Rib(i,j,k-1))
          endif
#  ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
#  endif
        enddo     !--> discard kbl
      enddo
# else
      ERROR: Algorithm of Rib is not specified.
# endif

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

# ifdef SMOOTH_HBL
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=J_EXT_RANGE
          hbl(istr-1,j)=hbl(istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=J_EXT_RANGE
          hbl(iend+1,j)=hbl(iend,j)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=I_EXT_RANGE
          hbl(i,jstr-1)=hbl(i,jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=I_EXT_RANGE
          hbl(i,jend+1)=hbl(i,jend)
        enddo
      endif
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        hbl(istr-1,jstr-1)=hbl(istr,jstr)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        hbl(istr-1,jend+1)=hbl(istr,jend)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        hbl(iend+1,jstr-1)=hbl(iend,jstr)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        hbl(iend+1,jend+1)=hbl(iend,jend)
      endif
#   endif               /* ! The actual smoothing via 2 passes     */
#  endif                /* ! of 4-point averaging; arrays "wm,ws"  */
      do j=jstr,jend+1     ! are used as scratch arrays here.
        do i=istr,iend+1
          wm(i,j)=rmask(i,j)+rmask(i-1,j)+rmask(i,j-1)+rmask(i-1,j-1)
          ws(i,j)=  hbl(i,j)  +hbl(i-1,j)  +hbl(i,j-1)  +hbl(i-1,j-1)
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          cff_dn=wm(i,j)+wm(i+1,j)+wm(i,j+1)+wm(i+1,j+1)
          cff_up=ws(i,j)+ws(i+1,j)+ws(i,j+1)+ws(i+1,j+1)
          if (cff_dn>0) hbl(i,j)=cff_up/cff_dn
        enddo
      enddo                     !--> discard wm,ws
# endif
# undef I_EXT_RANGE
# undef J_EXT_RANGE
# define imin illegal
# define imax illegal
# define jmin illegal
# define jmax illegal

!  Find stability and buoyancy forcing "Bfsfc" for boundary layer.
!  Insure that Bfsfc is never equal to zero

#define swdk wm
      call lmd_swfrac_tile (istr,iend,jstr,jend,
     &                      istr,iend,jstr,jend, hbl,swdk)
      do j=jstr,jend
        do i=istr,iend
          Bfsfc(i,j)=Bo(i,j)+Bosol(i,j)*(1.-swdk(i,j))
        enddo
      enddo      !--> discard swdk
#undef swdk

!  Compare "hbl" with physically limiting cases (Ekman depth
!  and Monin-Obukhov depth) and correct.

      do j=jstr,jend
        do i=istr,iend
          if (Bfsfc(i,j)>0.) then
c**            hekman=cekman*ustar(i,j)/max(abs(f(i,j)),eps)
c**            hbl(i,j)=min(hbl(i,j), hekman)
c**

c**            hbl(i,j)=min(hbl(i,j), h_MO)
          endif

# ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
# endif
          kbl(i,j)=N             !<-- initialize search
        enddo
      enddo
                                 ! Find new boundary
      do k=N-1,1,-1              ! layer index "kbl".
        do j=jstr,jend
          do i=istr,iend
            if (z_w(i,j,k) > z_w(i,j,N)-hbl(i,j)) kbl(i,j)=k
          enddo
        enddo
      enddo

! Find buoyancy forcing for final "hbl" values, and compute
! tubulent velocity scales (wm,ws) at "hbl".

#define swdk wm
      call lmd_swfrac_tile (istr,iend,jstr,jend, istr,iend,jstr,jend,
     &                                                      hbl,swdk)

      do j=jstr,jend
        do i=istr,iend
          Bfsfc(i,j)=Bo(i,j)+Bosol(i,j)*(1.-swdk(i,j))
#undef swdk
# ifdef MASKING
          Bfsfc(i,j)=Bfsfc(i,j)*rmask(i,j)
# endif
# ifdef LIMIT_UNSTABLE_ONLY
          if (Bfsfc(i,j)<0.) then
            zscale(i,j)=hbl(i,j)*epssfc
          else
            zscale(i,j)=hbl(i,j)
          endif
# else
          zscale(i,j)=hbl(i,j)*epssfc
# endif
        enddo
      enddo

      call lmd_wscale_tile (istr,iend,jstr,jend, istr,iend,jstr,jend,
     &                                   ustar, Bfsfc, zscale, wm,ws)

! Compute nondimensional shape function coefficients Gx( ) by
! matching values and vertical derivatives of  interior mixing
! coefficients at hbl (sigma=1).

      do j=jstr,jend
        do i=istr,iend
# ifdef LIMIT_UNSTABLE_ONLY
          f1=5.0 * max(0., Bfsfc(i,j)) * vonKar/(ustar(i,j)**4+eps)
# else
          f1=0.
#endif
          k=kbl(i,j)
          z_bl=z_w(i,j,N)-hbl(i,j)

          cff=1.D0/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)

c      if (cff_up<0. .or. cff_dn<0. ) then
c        write(*,*) cff_up,cff_dn, k
c      endif


          Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
          Gm1(i,j)=Av_bl/(hbl(i,j)*wm(i,j)+eps)
          dGm1dS(i,j)=min(0.D0, Av_bl*f1-dAv_bl/(wm(i,j)+eps))

          At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
          Gt1(i,j)=At_bl/(hbl(i,j)*ws(i,j)+eps)
          dGt1dS(i,j)=min(0.D0, At_bl*f1-dAt_bl/(ws(i,j)+eps))

# ifdef SALINITY
          As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          Gs1(i,j)=As_bl/(hbl(i,j)*ws(i,j)+eps)
          dGs1dS(i,j)=min(0.D0, As_bl*f1-dAs_bl/(ws(i,j)+eps))
# endif
        enddo
      enddo

! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.

      do k=1,N-1
        do j=jstr,jend
          do i=istr,iend
            zscale(i,j)=z_w(i,j,N)-z_w(i,j,k)
# ifdef LIMIT_UNSTABLE_ONLY
            if (Bfsfc(i,j)<0.) then
              zscale(i,j)=min(zscale(i,j), hbl(i,j)*epssfc)
            endif
# else
            zscale(i,j)=min(zscale(i,j), hbl(i,j)*epssfc)
# endif
          enddo
        enddo
        call lmd_wscale_tile (istr,iend,jstr,jend,
     &                        istr,iend,jstr,jend,
     &            ustar, Bfsfc, zscale, wm,ws)

! Compute vertical mixing coefficients

        do j=jstr,jend
          do i=istr,iend
            if (k>=kbl(i,j)) then
              sigma=(z_w(i,j,N)-z_w(i,j,k))/(hbl(i,j)+eps)
# ifdef MASKING
              sigma=sigma*rmask(i,j)
# endif
              a1=sigma-2.
              a2=3.-2.*sigma
              a3=sigma-1.

              Gm=a1+a2*Gm1(i,j)+a3*dGm1dS(i,j)
              Gt=a1+a2*Gt1(i,j)+a3*dGt1dS(i,j)
# ifdef SALINITY
              Gs=a1+a2*Gs1(i,j)+a3*dGs1dS(i,j)
# endif

              Kv(i,j,k)=hbl(i,j)*wm(i,j)*sigma*(1.+sigma*Gm)
              Kt(i,j,k)=hbl(i,j)*ws(i,j)*sigma*(1.+sigma*Gt)
# ifdef SALINITY
              Ks(i,j,k)=hbl(i,j)*ws(i,j)*sigma*(1.+sigma*Gs)
# endif
c**              Kv(i,j,k)=min(nubl,Kv(i,j,k))
c**              Kt(i,j,k)=min(nubl,Kt(i,j,k))
# ifdef SALINITY
c**              Ks(i,j,k)=min(nubl,Ks(i,j,k))
# endif
# ifdef LMD_NONLOCAL
              if (Bfsfc(i,j) < 0.) then
                ghat(i,j,k)=Cg/(ws(i,j)*hbl(i,j)+eps)
              endif
# endif
            else
# ifdef LMD_NONLOCAL
              ghat(i,j,k)=0.
# endif
# ifdef LMD_CONVEC
              if (bvf(i,j,k)<0.) then     ! Add convective
                Kv(i,j,k)=Kv(i,j,k) + nu0c   ! adjustment
                Kt(i,j,k)=Kt(i,j,k) + nu0c
#  ifdef SALINITY
                Ks(i,j,k)=Ks(i,j,k) + nu0c
#  endif
              endif
# endif
            endif
          enddo    /* !!!! This loop was not pipelined */
        enddo

!  Finalize: Set vertical mixing coefficients "Akv" and "Akt" using a
!  time  filter to avoid temporal instability of the implicit scheme.

        do j=jstr,jend
          do i=istr,iend
            Akv(i,j,k)=Kv(i,j,k)
            Akt(i,j,k,itemp)=Kt(i,j,k)

# ifdef SALINITY
            Akt(i,j,k,isalt)=Ks(i,j,k)
# endif

# ifdef BIOLOGY
            do itrc=3,NT
              Akt(i,j,k,itrc)=Ks(i,j,k)
            enddo
# endif
          enddo
        enddo
      enddo

      do j=jstr,jend                  ! Copy "hbl" into its shared
        do i=istr,iend                ! array and pad out alues at
          hbls(i,j,3-nstp)=hbl(i,j)   ! lateral side-boundaries.
        enddo
      enddo
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          hbls(istr-1,j,3-nstp)=hbls(istr,j,3-nstp)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          hbls(iend+1,j,3-nstp)=hbls(iend,j,3-nstp)
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          hbls(i,jstr-1,3-nstp)=hbls(i,jstr,3-nstp)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          hbls(i,jend+1,3-nstp)=hbls(i,jend,3-nstp)
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and.
     &    SOUTHERN_EDGE) then
        hbls(istr-1,jstr-1,3-nstp)=hbls(istr,jstr,3-nstp)
      endif
      if (WESTERN_EDGE .and.
     &    NORTHERN_EDGE) then
        hbls(istr-1,jend+1,3-nstp)=hbls(istr,jend,3-nstp)
      endif
      if (EASTERN_EDGE .and.
     &    SOUTHERN_EDGE) then
        hbls(iend+1,jstr-1,3-nstp)=hbls(iend,jstr,3-nstp)
      endif
      if (EASTERN_EDGE .and.
     &    NORTHERN_EDGE) then
        hbls(iend+1,jend+1,3-nstp)=hbls(iend,jend,3-nstp)
      endif
#  endif
# endif

# ifdef EXCHANGE
#  ifdef SALINITY
      call exchange_4_tile ( istr,iend,jstr,jend,
     &                 Akt(START_2D_ARRAY,0,isalt), N+1,
#  else
      call exchange_3_tile ( istr,iend,jstr,jend,
#  endif
     &       Akv, N+1,  Akt(START_2D_ARRAY,0,itemp), N+1,
     &                    hbls(START_2D_ARRAY,3-nstp), 1
     &                                                 )
# endif
      return
      end


# undef imin
# undef imax
# undef jmin
# undef jmax


      subroutine lmd_wscale_tile (istr,iend,jstr,jend,
     &                            imin,imax,jmin,jmax,
     &                     ustar, Bfsfc,zscale, wm,ws)

! Compute turbulent velocity scales for momentum and tracer
!
! input:  Bfsfc   boyancy forcing
!          zscale  boundary layer depth [m].
!
! output: wm   turbulent velocity scale [m/s] for momentum
!         ws   turbulent velocity scale [m/s] for tracer
!
! This routine was adapted from Bill Large 1995 code.

      implicit none
# include "param.h"
# include "scalars.h"
# ifdef MASKING
#  include "grid.h"
# endif

      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j
      real ustar(PRIVATE_2D_SCRATCH_ARRAY),
     &     Bfsfc(PRIVATE_2D_SCRATCH_ARRAY),
     &    zscale(PRIVATE_2D_SCRATCH_ARRAY),   ustar3,   r2,
     &        wm(PRIVATE_2D_SCRATCH_ARRAY),   zetahat,  r3,
     &        ws(PRIVATE_2D_SCRATCH_ARRAY),     r4
      parameter (r2=0.5, r3=1./3., r4=0.25)

      real zetam,am,cm,  zetas,as,cs
      parameter (
     &    zetam=-0.2, ! Maximum stability parameters "zeta"
     &    am=1.257,   ! value of the 1/3 power law regime of
     &    cm=8.360,   ! flux profile for momentum and tracers
     &    zetas=-1.0, ! and coefficients of flux profile for
     &    as=-28.86,  ! momentum and tracers in their 1/3-power
     &    cs=98.96)   ! law regime;


      do j=jmin,jmax
        do i=imin,imax
          ustar3=ustar(i,j)*ustar(i,j)*ustar(i,j)
          zetahat=vonKar*zscale(i,j)*Bfsfc(i,j)
# ifdef MASKING
          zetahat=zetahat*rmask(i,j)
# endif

! Stable regime.

          if (zetahat>=0.) then
            wm(i,j)=vonKar*ustar(i,j)*ustar3/max( ustar3
     &                               +5.*zetahat, 1.E-20)
            ws(i,j)=wm(i,j)

! Unstable regime, note that zetahat is always negative here, also
! negative are constants "zetam" and "zetas".

          else
            if (zetahat > zetam*ustar3) then
              wm(i,j)=vonKar*( ustar(i,j)*(ustar3-16.*zetahat) )**r4
            else
              wm(i,j)=vonKar*(am*ustar3-cm*zetahat)**r3
            endif
            if (zetahat > zetas*ustar3) then
              ws(i,j)=vonKar*( (ustar3-16.*zetahat)/ustar(i,j) )**r2
            else
              ws(i,j)=vonKar*(as*ustar3-cs*zetahat)**r3
            endif
          endif
        enddo       /* !!! This loop wan not pipelined. */
      enddo
      return
      end
#else
      subroutine lmd_kpp_empty
      end
#endif /* LMD_KPP */

