#include "cppdefs.opt"

# define INT_AT_RHO_POINTS
# define SMOOTH_HBL

      subroutine lmd_kpp (ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,
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
      use roms_mpi

      implicit none
      integer :: i,j,k
      real, dimension(GLOBAL_2D_ARRAY) :: ustar, Bo, Bosol
      real, dimension(GLOBAL_2D_ARRAY) :: hbl,bbl
      real, dimension(GLOBAL_2D_ARRAY) :: FX,FE,FE1
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:nz) :: Cr,FC, wrk1,wrk2
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

      real Kv0, Kt0, wmb,wsb, sgmb,amp

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
c-#  define alpha Bosol
c-#  define beta  Bo
      call alfabeta_tile (Bosol,Bo)
# endif
      do j=0,ny+1
        do i=0,nx+1
          Bo(i,j)=g*( Bosol(i,j)*(stflx(i,j,itemp)-srflx(i,j))
     &                              -bo(i,j)*stflx(i,j,isalt)
     &                                                        )
          Bosol(i,j)=g*Bosol(i,j)*srflx(i,j)
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
# ifdef LMD_BKPP          
          bbl(i,j)=hbbl(i,j)     ! as the initial guess
# endif
          kbls(i)=0
          Cr(i,nz)=0.
          Cr(i,0)=0.
          FC(i,nz)=0.
        enddo  !--> discard alpha,beta; keep Bo,Bosol to the very end.



# if defined BULK_DIFFERENCE

# elif defined INT_AT_W_POINTS

# elif defined INT_AT_RHO_POINTS

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

        do k=nz-1,1,-1
          do i=0,nx+1
            cff_up=(z_w(i,j,nz)-z_w(i,j,k))**2
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


        do i=0,nx+1
          z_bl=z_w(i,j,0)+0.25*Hz(i,j,1)
          cff_up=(z_w(i,j,nz)-z_bl)**2
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
        do k=nz,1,-1
          do i=0,nx+1
            swdk_r(i,k)=sqrt(swr_frac(i,j,k)*swr_frac(i,j,k-1))

            zscale=z_w(i,j,nz)-z_r(i,j,k)
            Bfsfc=Bo(i,j)+Bosol(i,j)*(1.-swdk_r(i,k))

#   include "lmd_wscale_ws_only.h"

c--         Vtsq=Vtc*ws*sqrt(max(0., 0.5*(bvf(i,j,k)+bvf(i,j,k-1)) ))

            Vtsq=1.8*Vtc*ws*sqrt(max(1.e-5, bvf(i,j,k-1) ))  !<-- fix by Jeroen

            Cr(i,k)=FC(i,k)+Vtsq

            if (kbls(i) == 0 .and. Cr(i,k) < 0.) kbls(i)=k
          enddo
        enddo

        do i=0,nx+1
          if (kbls(i) > 0) then
            k=kbls(i)
            if (k == nz) then
              hbl(i,j)=z_w(i,j,nz)-z_r(i,j,nz)
            else
              hbl(i,j)=z_w(i,j,nz)-( z_r(i,j,k)*Cr(i,k+1)
     &                              -z_r(i,j,k+1)*Cr(i,k)
     &                              )/(Cr(i,k+1)-Cr(i,k))
            endif
          else
            hbl(i,j)=z_w(i,j,nz)-z_w(i,j,0)
          endif
#   ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
#   endif
        enddo
#   undef swdk_r

        do i=0,nx+1
          kbbl(i)=0       ! reset Cr at bottom and kbls for BKPP
          Cr(i,0)=0.D0
        enddo
        do k=1,nz,+1
          do i=0,nx+1
            Cr(i,k)=FC(i,k)-FC(i,0)
            if (kbbl(i) == 0 .and. Cr(i,k) > 0.) kbbl(i)=k
          enddo
        enddo
        do i=0,nx+1
          if (kbbl(i) > 0) then
            k=kbbl(i)
            if (k == 1) then
              bbl(i,j)=z_r(i,j,1)-z_w(i,j,0)   ! no BBL
            else
              bbl(i,j)=( z_r(i,j,k-1)*Cr(i,k)-z_r(i,j,k)*Cr(i,k-1)
     &                            )/(Cr(i,k)-Cr(i,k-1)) -z_w(i,j,0)
            endif
          else
            bbl(i,j)=z_w(i,j,nz)-z_w(i,j,0) ! total depth
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

      do j=1,ny   !--> restart j-loop till the very end

# if defined LMD_KPP || defined LMD_BKPP
        if (.not. FIRST_TIME_STEP) then
# endif
          do i=1,nx
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
        do i=1,nx
          cff=z_w(i,j,nz)-z_w(i,j,0)          ! if surface and bottom
          if (hbl(i,j)+bbl(i,j) > cff) then  ! boundary layers overlap
            hbl(i,j) = cff                   ! set both of them to the
            bbl(i,j) = cff                   ! total depth of water
          endif                              ! column.
        enddo
# endif

# ifdef LMD_KPP
        do i=1,nx
          kbls(i)=nz         ! initialize search, then find
        enddo                ! new boundary layer index "kbls"
        do k=nz-1,1,-1
          do i=1,nx
            if (z_w(i,j,k) > z_w(i,j,nz)-hbl(i,j)) kbls(i)=k
          enddo
        enddo

! Compare "hbl" with physically limiting cases (Ekman depth and
! Monin-Obukhov depth) and correct.


! Find buoyancy forcing for final "hbl" values, and compute tubulent
! velocity scales (wm,ws) at "hbl".  Then compute nondimensional shape
! function coefficients Gx( ) by matching values and vertical
! derivatives of  interior mixing coefficients at hbl (ssgm=1).

        do i=1,nx
          k=kbls(i)
          z_bl=z_w(i,j,nz)-hbl(i,j)
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

        do i=1,nx
          do k=nz,0,-1
            Bfsfc=Bfsfc_bl(i)
            zscale=z_w(i,j,nz)-z_w(i,j,k)

#  include "lmd_wscale_wm_and_ws.h"

! Compute vertical mixing coefficients

            ssgm=(z_w(i,j,nz)-z_w(i,j,k))/max(hbl(i,j),EPS)

            if (ssgm < 1.) then
              if (ssgm<0.07D0) then
                cff=0.5*(ssgm-0.07D0)**2/0.07D0
              else
                cff=0.D0
              endif
              cff=cff + ssgm*(1.-ssgm)**2
              ! Add kpp values and previously computed values,
              ! tailed away near the surface.
              amp = ssgm**2
              AKv(i,j,k)=sqrt((amp*AKv(i,j,k))**2 + (wm*hbl(i,j)*cff)**2 )
              AKt(i,j,k)=sqrt((amp*AKt(i,j,k))**2 + (ws*hbl(i,j)*cff)**2 )

#  ifdef LMD_NONLOCAL
              if (Bfsfc < 0.) then
                ghat(i,j,k)=-Cg * ssgm*(1.-ssgm)**2
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

        do i=1,nx
          kbbl(i)=nz         !<-- initialize search
        enddo
        do k=nz-1,1,-1       ! find new boundary layer index "kbls"
          do i=1,nx
            if (z_w(i,j,k) > z_w(i,j,0)+bbl(i,j)) kbbl(i)=k
          enddo
        enddo

! Compute nondimensional shape function coefficients Gx( ) by matching
! values and vertical derivatives of interior mixing coefficients at
! bbl (sgmb=1).  Below wm,ws are turbulent velocity scales associated
! with bottom drag, assuming no buoyancy forcing effects.


        do i=1,nx
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

          do k=0,nz
            sgmb=(z_w(i,j,k)-z_w(i,j,0)+Zob)/(bbl(i,j)+Zob)
            if (sgmb < 1.) then

              cff1=sgmb*(1.-sgmb)**2

              AKv(i,j,k)=sqrt( AKv(i,j,k)**2 + (wmb*bbl(i,j)*cff1)**2 )
              AKt(i,j,k)=sqrt( AKt(i,j,k)**2 + (wsb*bbl(i,j)*cff1)**2 )

            endif
          enddo      !<-- k
        enddo      !<-- i


        do i=1,nx
# ifdef MASKING
          if (rmask(i,j) > 0.5) then
# endif
            if (FIRST_TIME_STEP) then
            else  !! Average 2 time steps of mixing coefficients
# ifdef WEC
              do k=1,nz     ! Apply wave effects
                Akv(i,j,k) = Akv(i,j,k) + Akb(i,j,k)
                Akt(i,j,k) = Akt(i,j,k) + Akb(i,j,k)+Akw(i,j,k)
              enddo
              Akv(i,j,0) = 2.*Akv(i,j,1)-Akv(i,j,2)
              Akt(i,j,0) = 2.*Akt(i,j,1)-Akt(i,j,2)
# endif

            endif  !<-- FIRST_TIME_STEP
# ifdef MASKING
          else
            do k=0,nz
              Akv(i,j,k) = 0.
              Akt(i,j,k) = 0.
            enddo
          endif
# endif
        enddo
      enddo

# ifdef LMD_KPP
      do j=1,ny                        ! Copy "hbl" into shared array
        do i=1,nx                      ! and padd lateral ghost points
          hbls(i,j)=hbl(i,j)
        enddo
      enddo
      if (west_bnd) then
        do j=1,ny
          hbls(0,j)=hbls(1,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          hbls(nx+1,j)=hbls(nx,j)
        enddo
      endif
      if (south_bnd) then
        do i=1,nx
          hbls(i,0)=hbls(i,1)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          hbls(i,ny+1)=hbls(i,ny)
        enddo
      endif
      if (west_bnd .and. south_bnd) then
        hbls(0,0)=hbls(1,1)
      endif
      if (west_bnd .and. north_bnd) then
        hbls(0,ny+1)=hbls(1,ny)
      endif
      if (east_bnd .and. south_bnd) then
        hbls(nx+1,0)=hbls(nx,1)
      endif
      if (east_bnd .and. north_bnd) then
        hbls(nx+1,ny+1)=hbls(nx,ny)
      endif
# endif   /* LMD_KPP */

# ifdef LMD_BKPP
      do j=1,ny                   ! Copy "bbl" into shared array
        do i=1,nx                 ! and padd lateral ghost points
          hbbl(i,j)=bbl(i,j)
        enddo
      enddo
      if (west_bnd) then
        do j=1,ny
          hbbl(0,j)=hbbl(1,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          hbbl(nx+1,j)=hbbl(nx,j)
        enddo
      endif
      if (south_bnd) then
        do i=1,nx
          hbbl(i,0)=hbbl(i,1)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          hbbl(i,ny+1)=hbbl(i,ny)
        enddo
      endif
      if (west_bnd .and. south_bnd) then
        hbbl(0,0)=hbbl(1,1)
      endif
      if (west_bnd .and. north_bnd) then
        hbbl(0,ny+1)=hbbl(1,ny)
      endif
      if (east_bnd .and. south_bnd) then
        hbbl(nx+1,0)=hbbl(nx,1)
      endif
      if (east_bnd .and. north_bnd) then
        hbbl(nx+1,ny+1)=hbbl(nx,ny)
      endif
# endif   /* LMD_BKPP */

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
# endif
      end
