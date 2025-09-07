#include "cppdefs.opt"
      module lmd_mixing

      ! Parameterizations of vertical mixing based on:
      ! Large, W. G., J. C. McWilliams, and S. C. Doney, 1994:
      ! Oceanic vertical mixing: a review and a model
      ! with a nonlocal boundary layer parameterization,
      ! Rev. Geophys., 32, 363-403.

      use dimensions
      use roms_mpi

      use tracers
      use grid
      use mixing
      use ocean_vars
      use scalars
      use surf_flux

#include "lmd_mixing.opt"

      ! work space for smooth_2d
      real,allocatable,dimension(:,:) :: fx,fy,fx1,fy1

      public lmd_vmix

      contains

!---------------------------------------------------------
      subroutine lmd_vmix(tind) ![
!     Compute Akv, Akt, Aks; mixing coefficients for 
!     momenta and tracers.

      implicit none

      integer,intent(in)  :: tind

      ! Local
      integer  :: i,j,k
      real :: nu_sx, cff,cff1, dudz,dvdz
      real :: dist,mult
      real,dimension(:,:),allocatable :: Rig

      integer  :: istr,iend,jstr,jend
 
      if (.not.allocated(fx)) then
        allocate(fx(0:nx+1,0:ny+1))
        allocate(fy(0:nx+1,0:ny+1))
        allocate(fx1(0:nx+1,0:ny+1))
        allocate(fy1(0:nx+1,0:ny+1))
      endif

       istr = 1
       jstr = 1
       iend = nx
       jend = ny


      if (lmd_rimix) then

        allocate(Rig(0:nx+1,0:ny+1))

!  Compute horizontal velocity shear (du/dz)^2+(dv/dz)^2 at
!  horizontal rho- and vertical w-points. Then compute gradient
!  Richardson number (already divided by its critical value.
!  Compute "interior" viscosities and diffusivities everywhere as
!  Richardson number instability due to resolved vertical shear
!  and internal wave breaking.

        do k=1,nz-1
          do j=0,ny+1
            do i=0,nx+1
              cff=0.5/(z_r(i,j,k+1)-z_r(i,j,k))
              dudz=cff*( u(i  ,j,k+1,tind)-u(i  ,j,k,tind)
     &                  +u(i+1,j,k+1,tind)-u(i+1,j,k,tind))
              dvdz=cff*( v(i,j  ,k+1,tind)-v(i,j  ,k,tind)
     &                  +v(i,j+1,k+1,tind)-v(i,j+1,k,tind))

              Rig(i,j) = bvf(i,j,k)/( Ri0*max(
     &                      dudz*dudz+dvdz*dvdz,1.e-10 ))
            enddo
          enddo

          ! We're using 4 2d arrays inside. Possibly pass work space or
          ! use pointers

          if (smooth_ri) call smooth_2d(Rig)

          do j=1,ny
            do i=1,nx

              cff=min(1., max(0., Rig(i,j)))
              nu_sx=1. - cff*cff
              nu_sx=nu_sx*nu_sx*nu_sx
              AKv(i,j,k)=nuwm + nu0m*nu_sx
              AKt(i,j,k)=nuws + nu0s*nu_sx

              if (lmd_convec.and.(Rig(i,j)<0.) ) then
                AKv(i,j,k) = AKv(i,j,k) + nu0c
                AKt(i,j,k) = AKt(i,j,k) + nu0c
              endif
            enddo
          enddo
        enddo

      else  ! internal wave breaking only

        do k=1,nz-1
          do j=1,ny
            do i=1,nx
              AKv(i,j,k) = nuwm
              AKt(i,j,k) = nuws
            enddo
          enddo
        enddo
      endif

! Supress mixing away near the bottom.
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            dist = z_w(i,j,k)-z_w(i,j,0)
            if (dist<Lturb) then
              mult = sin(0.5*pi*(z_w(i,j,k) - z_w(i,j,0))/Lturb)
              AKv(i,j,k) = AKv(i,j,k)*mult
              AKt(i,j,k) = AKt(i,j,k)*mult
            endif
          enddo
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          AKv(i,j,nz)=AKv(i,j,nz-1)
          AKt(i,j,nz)=AKt(i,j,nz-1)
          AKv(i,j,0 )=AKv(i,j,   1)
          AKt(i,j,0 )=AKt(i,j,   1)
        enddo
      enddo

      ! vertical smoothing of interior mixing
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            AKv(i,j,k)=0.5*AKv(i,j,k)+0.25*AKv(i,j,k-1)+0.25*AKv(i,j,k+1)
            AKt(i,j,k)=0.5*AKt(i,j,k)+0.25*AKt(i,j,k-1)+0.25*AKt(i,j,k+1)
          enddo
        enddo
      enddo

      ! Add surface and bottom boundary layer mixing
      if (lmd_kpp.or.lmd_bkpp) then
        call lmd_tbl(tind)
      endif

      !! Do we really need to exchange the mixing coefficients?
      call exchange_xxx(Akv,Akt)

      deallocate(Rig)

      end   subroutine lmd_vmix !]
!---------------------------------------------------------
      subroutine lmd_tbl(tind)   ![

# define INT_AT_RHO_POINTS
# define SMOOTH_HBL

      implicit none
      
      ! Local
      integer i,j,k
      real, dimension(1-bf:nx+bf,1-bf:ny+bf) :: ustar, Bo, Bosol
     &                                                  , hbl, bbl
# ifdef SMOOTH_HBL
     &                                              , FX, FE, FE1
# else
      real FX, FE, FE1  !<-- not in use, except as arguments
# endif
      real, dimension(0:nx+1,0:nz) :: Cr,FC, wrk1,wrk2
      real, dimension(0:nx+1) ::  Bfsfc_bl,
     &                               Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS
      integer, dimension(0:nx+1) :: kbls, kmo,kbbl
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

      real,parameter :: r2=0.5, r3=1./3., r4=0.25, EPS=1.E-20

      real,parameter :: Cg = Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
      real,parameter :: Vtc= Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)


      real ustar3, Bfsfc, zscale, zetahat, ws,wm, Kern, Vtsq,
     &         ssgm, z_bl, cff,cff1, cff_up,cff_dn
!     real Cg, ustar3, Bfsfc, zscale, zetahat, ws,wm, Kern, Vtc,Vtsq,
!    &         ssgm, z_bl, cff,cff1, cff_up,cff_dn

      real Kv0, Kt0, Ks0, wmb,wsb, sgmb,amp
      integer  :: istr,iend,jstr,jend



      !! When smoothing hbl, we try to not 
      !! smooth over the mask and the domain
      !! boundaries
# ifdef SMOOTH_HBL
      integer imin,imax,jmin,jmax
#  ifdef EW_PERIODIC
      imin=0
      imax=nx+1
#  else
      if (WESTERN_EDGE) then
        imin=1
      else
        imin=0
      endif
      if (EASTERN_EDGE) then
        imax=nx
      else
        imax=nx+1
      endif
#  endif
#  ifdef NS_PERIODIC
      jmin=0
      jmax=ny+1
#  else
      if (SOUTHERN_EDGE) then
        jmin=1
      else
        jmin=0
      endif
      if (NORTHERN_EDGE) then
        jmax=ny
      else
        jmax=ny+1
      endif
#  endif
       istr = 1
       jstr = 1
       iend = nx
       jend = ny
#  define I_EXT_RANGE imin,imax
#  define J_EXT_RANGE jmin,jmax
# else
#  define I_EXT_RANGE 1,nx
#  define J_EXT_RANGE 1,ny
# endif

! Nondimensional constants for computing non-local flux and convective
! deepening of surface boundary layer. The only reason why they are not
! parameters and computed here is because of power function involved.

!     Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
!     Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)

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
      call alfabeta_tile (1,nx,1,ny, I_EXT_RANGE,
     &                             J_EXT_RANGE, alpha,beta)
# endif
      ! JM if we move away from a virtual salt flux, this will have to 
      ! change: stflx(isalt) = - pme*sss
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          Bo(i,j)=g*( alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
     &                              -beta(i,j)*stflx(i,j,isalt)
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
          cff=z_w(i,j,N)-z_w(i,j,0)          ! if surface and bottom
          if (hbl(i,j)+bbl(i,j) > cff) then  ! boundary layers overlap
            hbl(i,j) = cff                   ! set both of them to the
            bbl(i,j) = cff                   ! total depth of water
          endif                              ! column.
        enddo
# endif

# ifdef LMD_KPP
        do i=1,nx
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
              AKv(i,j,k)=sqrt((amp*AKv(i,j,k))**2 + (wm*hbl(i,j)*cff)**2 )
              AKt(i,j,k)=sqrt((amp*AKt(i,j,k))**2 + (ws*hbl(i,j)*cff)**2 )
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

          do k=0,nz
            sgmb=(z_w(i,j,k)-z_w(i,j,0)+Zob)/(bbl(i,j)+Zob)
            if (sgmb < 1.) then

              cff1=sgmb*(1.-sgmb)**2

              AKv(i,j,k)=sqrt( AKv(i,j,k)**2 + (wmb*bbl(i,j)*cff1)**2 )
              AKt(i,j,k)=sqrt( AKt(i,j,k)**2 + (wsb*bbl(i,j)*cff1)**2 )

            endif
          enddo      !<-- k
        enddo      !<-- i


        do i=istr,iend
# ifdef MASKING
          if (rmask(i,j) > 0.5) then
# endif
            if (FIRST_TIME_STEP) then
              do k=0,nz
!               Akv(i,j,k)=Kv(i,j,k)
!               Akt(i,j,k,itemp)=Kt(i,j,k)
              enddo
            else  !! Average 2 time steps of mixing coefficients
              do k=0,N
!               Akv(i,j,k)       = 0.5*Akv(i,j,k)       + 0.5*Kv(i,j,k)
!               Akt(i,j,k,itemp) = 0.5*Akt(i,j,k,itemp) + 0.5*Kt(i,j,k)
!               Akv(i,j,k)       = Kv(i,j,k)
!               Akt(i,j,k,itemp) = Kt(i,j,k)
              enddo

# ifdef WEC
              do k=1,N      ! Apply wave effects
                Akv(i,j,k) = Akv(i,j,k)+Akb(i,j,k)
                Akt(i,j,k) = Akt(i,j,k)+Akb(i,j,k)+Akw(i,j,k)
              enddo
              Akv(i,j,0) = 2.*Akv(i,j,1)-Akv(i,j,2)
              Akt(i,j,0) = 2.*Akt(i,j,1)-Akt(i,j,2,itemp)
# endif

            endif  !<-- FIRST_TIME_STEP
# ifdef MASKING
          else
            do k=0,nz
              Akv(i,j,k)=0.
              Akt(i,j,k)=0.
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
      call exchange_xxx(hbls)
      call exchange_xxx(hbbl)
#   else
      call exchange_xxx(Akv)
      call exchange_xxx(hbls)
#   endif
#  else
      call exchange_xxx(Akv)
      call exchange_xxx(hbbl)
#  endif
      call exchange_xxx(Akt)
# endif
      end subroutine lmd_tbl !]
!---------------------------------------------------------
      subroutine smooth_2d(f) ![

! Smooth function in 2 dimensions ![
! Apply horizontal smoothing operator to a private array "wrk" (which
! is actually a CPP-alias for "hbl" or "hbbl").  The array is computed
! within a tile with one row of extra points all-around except lateral
! physical boundaries. The smoothed field is computed over the internal
! range of indices within the tile and is placed back to the same array
! "wrk".  The smoothing is performed while avoiding values under land
! mask which is accomplished by expressing everything via elementary
! differences subject to masking by U- and V-rules (hence applying
! Neumann condition at the coastline).  Without masking the stencil of
! the smoothing operator has the following weights, depending on
! coefficient settings in the code segment below:
!
!   cff=0, cff1=1/8        cff=1/12, cff1=3/16        cff=1/8, cff1=1/4
!
!         1/8                1/32  1/8  1/32           1/16  1/8  1/16
!
!    1/8  1/2  1/8           1/8   3/8  1/8            1/8   1/4  1/8
!
!         1/8                1/32  1/8  1/32           1/16  1/8  1/16
!
!       5-point                 isotropic                 2D 1-2-1
!      Laplacian                Laplacian                 Hanning
!
! All three smoothing operators suppress the checkerboard mode in just
! after a single iteration, however, only the last one eliminates flat-
! front 2dx-modes in one iteration; the first and the second attenuate
! the 2dx-mode by factors of 1/2 and 1/4 per iteration.   !]

      implicit none

      real,dimension(0:nx+1,0:ny+1),intent(inout) :: f
      real,parameter ::cff = 1.0/12.0 
      real,parameter ::cff1= 3.0/16.0

      ! Local
      integer :: i,j


      if (.not.ew_periodic) then
        if (inode.eq.0) then
          f(0,0:ny+1) = f(1,0:ny+1)
        endif
        if (inode.eq.npx-1) then
          f(nx+1,0:ny+1) = f(nx,0:ny+1)
        endif
      endif
      if (.not.ns_periodic) then
        if (jnode.eq.0) then
          f(0:nx+1,0) = f(0:nx+1,1)
        endif
        if (jnode.eq.npy-1) then
          f(0:nx+1,ny+1) = f(0:nx+1,ny)
        endif
      endif 

     ! Isotropic smoothing by first computing masked           
     ! elementary differences in each direction,
     ! after which expressing the transversal terms
     ! via the very same masked differences.

      do j=0,ny+1
        do i=1,nx+1
          fx1(i,j)=(f(i,j)-f(i-1,j))*umask(i,j)
        enddo
      enddo
      do j=1,ny+1
        do i=0,nx+1
          fy1(i,j)=(f(i,j)-f(i,j-1))*vmask(i,j)
        enddo
      enddo
      do j=1,ny
        do i=1,nx+1
          fx(i,j)=fx1(i,j) + cff*( fy1(i,j+1)+fy1(i-1,j  )
     &                            -fy1(i,j  )-fy1(i-1,j+1))
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx
          fy(i,j)=fy1(i,j) + cff*( fx1(i+1,j)+fx1(i  ,j-1)
     &                            -fx1(i  ,j)-fx1(i+1,j-1))
        enddo
      enddo
      do j=1,ny
        do i=1,nx
          f(i,j) = f(i,j) + cff1*( fx(i+1,j)-fx(i,j)
     &                            +fy(i,j+1)-fy(i,j))
        enddo
      enddo 

      end subroutine smooth_2d  !]
!---------------------------------------------------------

      end module lmd_mixing
