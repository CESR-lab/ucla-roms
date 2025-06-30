#include "cppdefs.opt"
#if defined SOLVE3D && defined LMD_KPP

# define INT_AT_RHO_POINTS
# define SMOOTH_HBL

      subroutine lmd_kpp_tile (istr,iend,jstr,jend,       Kv,Kt,Ks,
     &                         ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,
     &                         Cr,FC, wrk1,wrk2,
     &                         Gm1,dGm1dS,  Gt1,dGt1dS,  Gs1,dGs1dS,
     &                                                   kbls, kmo,
     &                                                         tind)

      use param
      use wec_frc
      use surf_flux ! bulk_frc-> only: sustr_r, svstr_r, stflx, srflx
      use grid      ! flux_frc-> only: sustr, svstr, stflx, srflx
      use mixing
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: Kv, Kt, Ks
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: ustar, Bo, Bosol
     &                                                  , hbl, bbl
# ifdef SMOOTH_HBL
     &                                              , FX, FE, FE1
# else
      real FX, FE, FE1  !<-- not in use, except as arguments
# endif
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: Cr,FC, wrk1,wrk2
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY) ::  Bfsfc_bl,
     &                               Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS
      integer, dimension(PRIVATE_1D_SCRATCH_ARRAY) :: kbls, kmo
     &                                                    , kbbl
      integer :: tind

      real, parameter ::
      ! Ricr - chosen as 0.15 by JM & PD 2020-09
     &   Ricr=0.15,     ! Critical bulk Richardson number (0.3)
     &   Ri_inv=1./Ricr,
     &   epssfc=0.1,    ! nondimensional extent of the surface layer
     &   betaT=-0.2,    ! ratio of entrainment flux to surface buoyancy
                        !      forcing flux (the "20% convectin rule")
     &   nubl=0.01,     ! maximum allowed boundary layer
     &   nu0c=0.1,      ! convective adjustment for viscosity and
                        !                      diffusivity [m^2/s]
     &   Cv=1.8,        ! ratio of interior Brunt-Vaisala frequency
                        ! "N" to that at the entrainment depth "he".
     &   C_MO=1.,       ! constant for computaion Monin-Obukhov depth.
     &   C_Ek=258.,     ! constant for computating stabilization term
                        !    due to Coriolis force (Ekman depth limit)
     &   Cstar=10.,     ! proportionality coefficient parameterizing
                        !                         nonlocal transport
     &   zeta_m=-0.2,   ! Maximum stability parameters "zeta"
     &   a_m=1.257,     ! value of the 1/3 power law regime of
     &   c_m=8.360,     ! flux profile for momentum and tracers
     &   zeta_s=-1.0,   ! and coefficients of flux profile for
     &   a_s=-28.86,    ! momentum and tracers in their 1/3-power
     &   c_s=98.96      ! law regime;

      real, parameter :: r2=0.5, r3=1./3., r4=0.25, EPS=1.E-20


      real Cg, ustar3, Bfsfc, zscale, zetahat, ws,wm, Kern, Vtc,Vtsq,
     &         ssgm, z_bl, cff,cff1, cff_up,cff_dn

      real Kv0, Kt0, Ks0, wmb,wsb, sgmb,amp

      !! When smoothing hbl, we try to not
      !! smooth over the mask and the domain
      !! boundaries
# ifdef SMOOTH_HBL
      integer imin,imax,jmin,jmax
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
#  define I_EXT_RANGE imin,imax
#  define J_EXT_RANGE jmin,jmax
# else
#  define I_EXT_RANGE istr,iend
#  define J_EXT_RANGE jstr,jend
# endif
# ifdef LMD_KPP
! Nondimensional constants for computing non-local flux and convective
! deepening of surface boundary layer. The only reason why they are not
! parameters and computed here is because of power function involved.

      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
      Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)
# endif

! Compute thermal expansion, "alpha" [kg/m^3/degC], and saline
! contraction, "beta" [kg/m^3/PSU], coefficients at surface; then
! compute surface turbulent buoyancy forcing "Bo" [m^2/s^3] (in doing
! so remove incoming solar shortwave radiation component and save it
! separately as "Bosol"). Also get an approximation for the surface
! layer depth using "epssfc" and boundary layer depth from previous
! time step (this is needed to estimate turbulent velocity scale in
! computation of "Vterm" in "Cr", before new hbl is found).  Also
! compute turbulent  friction velocity "ustar" from wind stress at
! RHO-points. Initialize boundary layer depth "hbl" and index "kbls"
! to the maximum (bottomed out) values.

# ifdef LMD_KPP
#  define alpha Bosol
#  define beta  Bo
      call alfabeta_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                             J_EXT_RANGE, alpha,beta)
# endif
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          Bo(i,j)=g*( alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
#  ifdef SALINITY
     &                              -beta(i,j)*stflx(i,j,isalt)
#  endif
     &                                                        )
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
#  undef beta
#  undef alpha

! DevinD & JM: better to use rho variables for sustr & svstr since we
! currently calculated both u/v and rho versions. If we converted u/v back to rho
! it would be a double average and dilute the answer further.
#   ifdef BULK_FRC
          ustar(i,j)=sqrt( sqrt( sustr_r(i,j)**2 +svstr_r(i,j)**2 ))
#   else
! Orginal for non-bulk-flux
          ustar(i,j)=sqrt( sqrt( 0.333333333333*(
     &      sustr(i,j)**2 + sustr(i+1,j)**2 + sustr(i,j)*sustr(i+1,j)
     &     +svstr(i,j)**2 + svstr(i,j+1)**2 + svstr(i,j)*svstr(i,j+1)
     &                                                            )))
#   endif /* BULK_FRC */

          hbl(i,j)=hbls(i,j)     ! use old-time-step values
          bbl(i,j)=hbbl(i,j)     ! as the initial guess

          kbls(i)=0
          Cr(i,N)=0.
          Cr(i,0)=0.
          FC(i,N)=0.
        enddo  !--> discard alpha,beta; keep Bo,Bosol to the very end.



# if defined BULK_DIFFERENCE

# elif defined INT_AT_W_POINTS

# elif defined INT_AT_RHO_POINTS

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

        do k=N-1,1,-1
          do i=I_EXT_RANGE
            cff_up=(z_w(i,j,N)-z_w(i,j,k))**2
            cff_dn=(z_w(i,j,k)-z_w(i,j,0))**2
            Kern=cff_up*cff_dn/( (cff_up +(epssfc*hbl(i,j))**2)
     &                          *(cff_dn +(epssfc*bbl(i,j))**2) )
            FC(i,k)=FC(i,k+1) + Kern*(
     &                0.5*( ( u(i,j,k+1,tind)+u(i+1,j,k+1,tind)
     &                       -u(i,j,k  ,tind)-u(i+1,j,k  ,tind) )**2
     &                     +( v(i,j,k+1,tind)+v(i,j+1,k+1,tind)
     &                       -v(i,j,k  ,tind)-v(i,j+1,k  ,tind) )**2
     &                      )/(Hz(i,j,k)+Hz(i,j,k+1))

     &               -0.5*(Hz(i,j,k)+Hz(i,j,k+1))*( Ri_inv*bvf(i,j,k)

     &                                            +C_Ek*f(i,j)*f(i,j)
     &                                                             ))
          enddo
        enddo


        do i=I_EXT_RANGE
          z_bl=z_w(i,j,0)+0.25*Hz(i,j,1)
          cff_up=(z_w(i,j,N)-z_bl)**2
          cff_dn=(z_bl-z_w(i,j,0))**2
          Kern=cff_up*cff_dn/( (cff_up +(epssfc*hbl(i,j))**2)
     &                        *(cff_dn +(epssfc*bbl(i,j))**2) )
            FC(i,0)=FC(i,1) + Kern*(
     &                   0.5*( (u(i,j,1,tind)+u(i+1,j,1,tind))**2
     &                        +(v(i,j,1,tind)+v(i,j+1,1,tind))**2
     &                       )/Hz(i,j,1)
     &                  -0.5*Hz(i,j,1)*( Ri_inv*bvf(i,j,1)
     &                                  +C_Ek*f(i,j)*f(i,j)
     &                                                       ))
        enddo

#   define swdk_r wrk1
        do k=N,1,-1
          do i=I_EXT_RANGE
            swdk_r(i,k)=sqrt(swr_frac(i,j,k)*swr_frac(i,j,k-1))

            zscale=z_w(i,j,N)-z_r(i,j,k)
            Bfsfc=Bo(i,j)+Bosol(i,j)*(1.-swdk_r(i,k))

#   include "lmd_wscale_ws_only.h"

c--         Vtsq=Vtc*ws*sqrt(max(0., 0.5*(bvf(i,j,k)+bvf(i,j,k-1)) ))

            Vtsq=Vtc*ws*sqrt(max(0., bvf(i,j,k-1) ))  !<-- fix by Gokhan

            Cr(i,k)=FC(i,k)+Vtsq

            if (kbls(i) == 0 .and. Cr(i,k) < 0.) kbls(i)=k
          enddo
        enddo

        do i=I_EXT_RANGE
          if (kbls(i) > 0) then
            k=kbls(i)
            if (k == N) then
              hbl(i,j)=z_w(i,j,N)-z_r(i,j,N)
            else
              hbl(i,j)=z_w(i,j,N)-( z_r(i,j,k)*Cr(i,k+1)
     &                              -z_r(i,j,k+1)*Cr(i,k)
     &                              )/(Cr(i,k+1)-Cr(i,k))
            endif
          else
            hbl(i,j)=z_w(i,j,N)-z_w(i,j,0)
          endif
#   ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
#   endif
        enddo
#   undef swdk_r

        do i=I_EXT_RANGE
          kbbl(i)=0       ! reset Cr at bottom and kbls for BKPP
          Cr(i,0)=0.D0
        enddo
        do k=1,N,+1
          do i=I_EXT_RANGE
            Cr(i,k)=FC(i,k)-FC(i,0)
            if (kbbl(i) == 0 .and. Cr(i,k) > 0.) kbbl(i)=k
          enddo
        enddo
        do i=I_EXT_RANGE
          if (kbbl(i) > 0) then
            k=kbbl(i)
            if (k == 1) then
              bbl(i,j)=z_r(i,j,1)-z_w(i,j,0)   ! no BBL
            else
              bbl(i,j)=( z_r(i,j,k-1)*Cr(i,k)-z_r(i,j,k)*Cr(i,k-1)
     &                            )/(Cr(i,k)-Cr(i,k-1)) -z_w(i,j,0)
            endif
          else
            bbl(i,j)=z_w(i,j,N)-z_w(i,j,0) ! total depth
          endif
#   ifdef MASKING
          bbl(i,j)=bbl(i,j)*rmask(i,j)
#   endif
        enddo     !--> discard FC, Cr and kbbl
# else
      ERROR: Algorithm for Bulk Richardson number is not specified.
# endif

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

# ifdef SMOOTH_HBL
      enddo  !<-- break j-loop if smoothing of hbl or bbl takes place

#  ifdef LMD_KPP
#   define wrk hbl
#   include "lmd_kpp_smooth_hbl.h"
#   undef wrk
#  endif
#  ifdef LMD_BKPP
#   define wrk bbl
#   include "lmd_kpp_smooth_hbl.h"
#   undef wrk
#  endif

      do j=jstr,jend   !--> restart j-loop till the very end
# endif /* SMOOTH_HBL */
# undef I_EXT_RANGE
# undef J_EXT_RANGE
# define imin illegal
# define imax illegal
# define jmin illegal
# define jmax illegal


# if defined LMD_KPP || defined LMD_BKPP
        if (.not. FIRST_TIME_STEP) then
# endif
          do i=istr,iend
# ifdef LMD_KPP
            hbl(i,j)=0.5*(hbl(i,j) + hbls(i,j))
# endif
# ifdef LMD_BKPP
            bbl(i,j)=0.5*(bbl(i,j) + hbbl(i,j))
# endif
          enddo
# if defined LMD_KPP || defined LMD_BKPP
        endif
# endif

# if defined MERGE_OVERLAP && defined LMD_KPP && defined LMD_BKPP
        do i=istr,iend
          cff=z_w(i,j,N)-z_w(i,j,0)          ! if surface and bottom
          if (hbl(i,j)+bbl(i,j) > cff) then  ! boundary layers overlap
            hbl(i,j) = cff                   ! set both of them to the
            bbl(i,j) = cff                   ! total depth of water
          endif                              ! column.
        enddo
# endif

# ifdef LMD_KPP
        do i=istr,iend
          kbls(i)=N          ! initialize search, then find
        enddo                ! new boundary layer index "kbls"
        do k=N-1,1,-1
          do i=istr,iend
            if (z_w(i,j,k) > z_w(i,j,N)-hbl(i,j)) kbls(i)=k
          enddo
        enddo

! Compare "hbl" with physically limiting cases (Ekman depth and
! Monin-Obukhov depth) and correct.


! Find buoyancy forcing for final "hbl" values, and compute tubulent
! velocity scales (wm,ws) at "hbl".  Then compute nondimensional shape
! function coefficients Gx( ) by matching values and vertical
! derivatives of  interior mixing coefficients at hbl (ssgm=1).

        do i=istr,iend
          k=kbls(i)
          z_bl=z_w(i,j,N)-hbl(i,j)
          zscale=hbl(i,j)

          if (swr_frac(i,j,k-1) > 0.) then
            Bfsfc=Bo(i,j) +Bosol(i,j)*( 1. -swr_frac(i,j,k-1)
     &              *swr_frac(i,j,k)*(z_w(i,j,k)-z_w(i,j,k-1))
     &               /( swr_frac(i,j,k  )*(z_w(i,j,k)   -z_bl)
     &                 +swr_frac(i,j,k-1)*(z_bl -z_w(i,j,k-1))
     &                                                     ) )
          else
            Bfsfc=Bo(i,j)+Bosol(i,j)
          endif

#  include "lmd_wscale_wm_and_ws.h"

          Bfsfc_bl(i)=Bfsfc
        enddo

! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.

        do i=istr,iend
          do k=N,0,-1
            Bfsfc=Bfsfc_bl(i)
            zscale=z_w(i,j,N)-z_w(i,j,k)

#  include "lmd_wscale_wm_and_ws.h"

! Compute vertical mixing coefficients

            ssgm=(z_w(i,j,N)-z_w(i,j,k))/max(hbl(i,j),EPS)

            if (ssgm < 1.) then
              if (ssgm<0.07D0) then
                cff=0.5*(ssgm-0.07D0)**2/0.07D0
              else
                cff=0.D0
              endif
              cff=cff + ssgm*(1.-ssgm)**2
              ! JM
              ! not entirely sure where this amp factor came from
              ! and whether we need it or not
              ! Check in with Jon on this
              !Kv(i,j,k)=Kv(i,j,k) + wm*hbl(i,j)*cff
              !Kt(i,j,k)=Kt(i,j,k) + ws*hbl(i,j)*cff
              amp = ssgm**2
              Kv(i,j,k)=sqrt((amp*Kv(i,j,k))**2 + (wm*hbl(i,j)*cff)**2 )
              Kt(i,j,k)=sqrt((amp*Kt(i,j,k))**2 + (ws*hbl(i,j)*cff)**2 )
#  ifdef SALINITY
              !Ks(i,j,k)=Ks(i,j,k) + ws*hbl(i,j)*cff
              Ks(i,j,k)=sqrt((amp*Ks(i,j,k))**2 + (ws*hbl(i,j)*cff)**2 )
#  endif

#  ifdef LMD_NONLOCAL
              if (Bfsfc < 0.) then
                ghat(i,j,k)=Cg * ssgm*(1.-ssgm)**2
              else
                ghat(i,j,k)=0.
              endif
#  endif
            else
#  ifdef LMD_NONLOCAL
              ghat(i,j,k)=0.
#  endif
            endif
          enddo
        enddo
# endif   /* LMD_KPP */

        do i=istr,iend
          kbbl(i)=N          !<-- initialize search
        enddo
        do k=N-1,1,-1       ! find new boundary layer index "kbls"
          do i=istr,iend
            if (z_w(i,j,k) > z_w(i,j,0)+bbl(i,j)) kbbl(i)=k
          enddo
        enddo

! Compute nondimensional shape function coefficients Gx( ) by matching
! values and vertical derivatives of interior mixing coefficients at
! bbl (sgmb=1).  Below wm,ws are turbulent velocity scales associated
! with bottom drag, assuming no buoyancy forcing effects.


        do i=istr,iend
          ! wmb = vonKar*ustar_bottom (rd*u)
          wmb=vonKar*vonKar*sqrt( 0.333333333333*(
     &               u(i,j,1,nstp)**2 +u(i+1,j,1,nstp)**2
     &                      +u(i,j,1,nstp)*u(i+1,j,1,nstp)
     &              +v(i,j,1,nstp)**2 +v(i,j+1,1,nstp)**2
     &                      +v(i,j,1,nstp)*v(i,j+1,1,nstp)
     &                                               ) )
     &                  / log(1.+0.5*Hz(i,j,1)/Zob)

          wsb=wmb

! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.

          do k=0,N
            sgmb=(z_w(i,j,k)-z_w(i,j,0)+Zob)/(bbl(i,j)+Zob)
            if (sgmb < 1.) then

              cff1=sgmb*(1.-sgmb)**2

              Kv(i,j,k)=sqrt( Kv(i,j,k)**2 + (wmb*bbl(i,j)*cff1)**2 )
              Kt(i,j,k)=sqrt( Kt(i,j,k)**2 + (wsb*bbl(i,j)*cff1)**2 )
#  ifdef SALINITY
              Ks(i,j,k)=sqrt( Ks(i,j,k)**2 + (wsb*bbl(i,j)*cff1)**2 )
#  endif

            endif
          enddo      !<-- k
        enddo      !<-- i


        do i=istr,iend
# ifdef MASKING
          if (rmask(i,j) > 0.5) then
# endif
            if (FIRST_TIME_STEP) then
              do k=0,N
                Akv(i,j,k)=Kv(i,j,k)
                Akt(i,j,k,itemp)=Kt(i,j,k)
# ifdef SALINITY
                Akt(i,j,k,isalt)=Ks(i,j,k)
# endif
              enddo
            else  !! Average 2 time steps of mixing coefficients
              do k=0,N
!               Akv(i,j,k)       = 0.5*Akv(i,j,k)       + 0.5*Kv(i,j,k)
!               Akt(i,j,k,itemp) = 0.5*Akt(i,j,k,itemp) + 0.5*Kt(i,j,k)
                Akv(i,j,k)       = Kv(i,j,k)
                Akt(i,j,k,itemp) = Kt(i,j,k)
# ifdef SALINITY
!               Akt(i,j,k,isalt) = 0.5*Akt(i,j,k,isalt) + 0.5*Ks(i,j,k)
                Akt(i,j,k,isalt) = Ks(i,j,k)
# endif
              enddo

# ifdef WEC
              do k=1,N      ! Apply wave effects
                Akv(i,j,k)       =Akv(i,j,k)+Akb(i,j,k)
                Akt(i,j,k,itemp) =Akt(i,j,k,itemp)+Akb(i,j,k)+Akw(i,j,k)
#  ifdef SALINITY
                Akt(i,j,k,isalt) =Akt(i,j,k,isalt)+Akb(i,j,k)+Akw(i,j,k)
#  endif
              enddo
              Akv(i,j,0)       =2.*Akv(i,j,1)-Akv(i,j,2)
              Akt(i,j,0,itemp) =2.*Akt(i,j,1,itemp)-Akt(i,j,2,itemp)
#  ifdef SALINITY
              Akt(i,j,0,isalt) =2.*Akt(i,j,1,isalt)-Akt(i,j,2,isalt)
#  endif
# endif

            endif  !<-- FIRST_TIME_STEP
# ifdef MASKING
          else
            do k=0,N
              Akv(i,j,k)=0.
              Akt(i,j,k,itemp)=0.
#  ifdef SALINITY
              Akt(i,j,k,isalt)=0.
#  endif
            enddo
          endif
# endif
        enddo
      enddo

# ifdef LMD_KPP
      do j=jstr,jend                   ! Copy "hbl" into shared array
        do i=istr,iend                 ! and padd lateral ghost points
          hbls(i,j)=hbl(i,j)
        enddo
      enddo
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          hbls(istr-1,j)=hbls(istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          hbls(iend+1,j)=hbls(iend,j)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          hbls(i,jstr-1)=hbls(i,jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          hbls(i,jend+1)=hbls(i,jend)
        enddo
      endif
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        hbls(istr-1,jstr-1)=hbls(istr,jstr)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        hbls(istr-1,jend+1)=hbls(istr,jend)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        hbls(iend+1,jstr-1)=hbls(iend,jstr)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        hbls(iend+1,jend+1)=hbls(iend,jend)
      endif
#   endif
#  endif
# endif   /* LMD_KPP */
# ifdef LMD_BKPP
      do j=jstr,jend                   ! Copy "bbl" into shared array
        do i=istr,iend                 ! and padd lateral ghost points
          hbbl(i,j)=bbl(i,j)
        enddo
      enddo
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          hbbl(istr-1,j)=hbbl(istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          hbbl(iend+1,j)=hbbl(iend,j)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          hbbl(i,jstr-1)=hbbl(i,jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          hbbl(i,jend+1)=hbbl(i,jend)
        enddo
      endif
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        hbbl(istr-1,jstr-1)=hbbl(istr,jstr)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        hbbl(istr-1,jend+1)=hbbl(istr,jend)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        hbbl(iend+1,jstr-1)=hbbl(iend,jstr)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        hbbl(iend+1,jend+1)=hbbl(iend,jend)
      endif
#   endif
#  endif
# endif   /* LMD_KPP */

# ifdef EXCHANGE
#  ifdef LMD_KPP
#   ifdef LMD_BKPP
      call exchange_xxx(Akv)
      call exchange_xxx(hbls(:,:))
      call exchange_xxx(hbbl(:,:))
#   else
      call exchange_xxx(Akv)
      call exchange_xxx(hbls(:,:))
#   endif
#  else
      call exchange_xxx(Akv)
      call exchange_xxx(hbbl(:,:))
#  endif
#  ifdef SALINITY
      call exchange_xxx(Akt(:,:,:,itemp))
      call exchange_xxx(Akt(:,:,:,isalt))
#  else
      call exchange_xxx(Akt(:,:,:,itemp))
#  endif
# endif
      end


      subroutine check_kpp_switches (ierr)

! This routine is to keep track of the status of CPP-switches defined
! locally in this file above. It must be placed here rather than stored
! in a separate file in order to be exposed to relevant CPP-settings.
! This routine does not affect model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      use param
      use strings

      implicit none
      integer ierr, is,ie, lenstr

      ie=lenstr(cpps)
      is=ie+2 ; ie=is+10
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<lmd_kpp.F>'
# ifdef INT_AT_RHO_POINTS
      is=ie+2 ; ie=is+16
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='INT_AT_RHO_POINTS'
# endif
# ifdef SMOOTH_HBL
      is=ie+2 ; ie=is+9
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='SMOOTH_HBL'
# endif
# ifdef MERGE_OVERLAP
      is=ie+2 ; ie=is+13
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='MERGE_OVERLAP'
# endif
      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
#else
      subroutine lmd_kpp_empty
      end
#endif /* LMD_KPP || LMD_BKPP */
