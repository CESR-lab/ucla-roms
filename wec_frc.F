      module wec_frc
      ! Wave Effect on Current (WEC) module which provides all the terms that enter
      ! in the momentum equations

       use dimensions !! has dimensions, a list of forcing files and rmask

       implicit none

! DevinD - need the following includes for GLOBAL_2D_ARRAY to work
#include "cppdefs.h"
#include "param.h"
#include "grid.h"

! DevinD - get universal scalars, g and pi
# include "scalars.h"
! DevinD - get ocean variables
# include "ocean2d.h"
# ifdef SOLVE3D
#   include "ocean3d.h"
# endif

       private

      ! Variables usable by code that uses this module
       real, public, dimension(GLOBAL_2D_ARRAY) :: brk2dx, brk2de !! xi/eta-direction 2D breaking dissipation (rho)
       real, public, dimension(GLOBAL_2D_ARRAY) :: frc2dx,frc2de  !! xi/eta-direction 2D frictional dissipation (rho)
       real, public, dimension(GLOBAL_2D_ARRAY) :: ust2d     !! xi-direction 2D Stokes transport (u-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: vst2d     !! eta-direction 2D Stokes transport (v-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: sup            !! quasi-static wave set-up (rho-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: calP           !! pressure correction term (rho-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: Kapsrf         !! Bernoulli head term at the surface (rho-point)
!      real, public, dimension(nx,ny) :: lm             !! Daniel's term: wavelength for diffusivity (not used)


       ! 3D terms
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: brk3dx, brk3de !! xi/eta-direction 3D breaking dissipation (rho)
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: ust       !! xi-direction 3D Stokes drift velocity (u-point)
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: vst       !! eta-direction 3D Stokes drift velocity (u-point)
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: wst,      !! vertical 3D Stokes drift velocity (rho-point)
     &                                                 Kappa,    !! 3D Bernoulli head term (rho-point)
     &                                                 kvf,      !! vertical vortex force term (K term, 3D, rho-point)
     &                                                 frc3dx, frc3de    !! xi/eta-direction 3D frictional dissipation (rho)


       ! terms with extra z-index
       real, public, dimension(GLOBAL_2D_ARRAY,0:nz) :: Akb,  !! breaking-wave-induced additional diffusivity (w-point
     &                                                  Akw,  !! wave-induced additional diffusivity (rho-point)
     &                                                  E_pre !! previous time-step value for Akw estimation (rho)


       real, public :: vstu,ustv,dudx,dvdx,dude,dvde


       logical :: analytical = .true.  ! analytical forcing
       
       ! Can be called from outside this module
       !public wstokes_tile

      !!!
      !!!!! DevinD - variables added as per DH's WAVE_PACKET forces.h for old MRL_WCI flag
      !!!

       ! wfrq | BBL/MRL | wind-induced wave frequency [rad/s]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wfrq
       ! wdrx | MRL     | cosine of wave direction [non dimension]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdrx
       ! wdre | MRL     | sine of   wave direction [non dimension]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdre
       ! whrm | MRL     | (RMS) wave height (twice the wave amplitude) [m]
       real, public, dimension(GLOBAL_2D_ARRAY) :: whrm
       ! wdsp | MRL     | breaking dissipation rate (\epsilon_b term) [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdsp
       ! wdrg | MRL     | frictional dissipation rate (\epsilon_d term) [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdrg
       ! wbst | MRL/BKPP| frictional dissipation stress (e_d k/sigma) [m2/s2]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wbst

       ! breaking dissipation [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wveb
       ! frictional dissipation [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wved


      !!!
      !!!!! DevinD - end variables added as per forces.h
      !!!

      ! DevinD - variables used in ana_wec_frc_tile and ana_init_wec_frc_tile

       ! For PART 1: Declare some variables from analytical.F
       real :: kkw, env, Clw, Cg



       ! DevinD - timestep variables as per river_frc.F
       integer              :: it1 = 1  ! step n. Note: should be parameter as unchanging
       integer              :: it2 = 2  ! step n+1. Note: should be parameter as unchanging
       real, dimension(2)   :: rtim = [-99,-99] ! used to store both time steps (n, n+1)


      public set_wec_frc
      public ana_wec_frc_tile
      public init_arrays_wec_tile
      public ana_init_wec_frc_tile


      contains 

       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------

        subroutine init_arrays_wec_tile (istr,iend,jstr,jend)

          ! This is called from init_arrays_tile of init_arrays.F
          ! It is needed to "first touch" arrays to optimize
          ! NUMA shared memory allocation on linux.

          implicit none

          integer, intent(in) :: istr, iend, jstr, jend
          ! Declare looping indices
          integer :: i, j

          real, parameter :: init=0.

          ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
          ! This had to declared last else won't compile
#include "compute_auxiliary_bounds.h"

          ! This section of code is taken from DH's WAVE_PACKET init_arrays.F


              whrm(i,j) =init
              wdrx(i,j) =init
              wdre(i,j) =init
              wbst(i,j) =init

          do j=jstrR,jendR
            do i=istrR,iendR
              sup(i,j)    =init
              ust2d(i,j)  =init
              vst2d(i,j)  =init
              brk2dx(i,j) =init
              brk2de(i,j) =init
              frc2dx(i,j) =init
              frc2de(i,j) =init
# ifdef SOLVE_3D
              calP(i,j)   =init
              Kapsrf(i,j) =init
              do k=1,N
#  ifndef SURFACE_BREAK
                brk3dx(i,j,k) =init
                brk3de(i,j,k) =init
#  endif
#  ifdef BODY_FRICTION
                frc3dx(i,j,k) =init
                frc3de(i,j,k) =init
#  endif
                ust(i,j,k) =init
                vst(i,j,k) =init
                wst(i,j,k) =init
                kvf(i,j,k) =init
              enddo
              do k=0,N
#  if defined LMD_KPP || defined LMD_BKPP
                Akb(i,j,k) =init
#  endif
                Akw(i,j,k) =init
                E_pre(i,j,k) =init
              enddo
# endif /* SOLVE_3D */
              wdsp(i,j) =init
#    ifdef SURFACE_ROLLER
              rdsp(i,j) =init
#    endif
            enddo
          enddo


        end subroutine init_arrays_wec_tile

       ! --------------------------------------------------------------------

        subroutine ana_init_wec_frc_tile (istr,iend,jstr,jend)

         ! Initialize wec forcing - chose to include because DH's WAVE_PACKET
         ! version is separate, and init version doesn't include time, whereas
         ! ana_wec_frc_tile does include time.
         ! Might be able to just add if statement to ana_wec_frc_tile instead?

         implicit none

         integer, intent(in) :: istr, iend, jstr, jend
         ! Declare looping indices
         integer :: i, j, k

         real :: cff1 ! This could be declared in module but not sure - here for now.

          ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
          ! This had to appear after other variables declared else wont compile
#include "compute_auxiliary_bounds.h"


         ! Taken from ana_init.F of DH's WAVE_PACKET code


          kkw=2.0*pi/1.0 ! Changed kw to kkw so as not to conflict part 2 also has kw
          env=0.0001
          Clw=sqrt(g*h(2,2))
          Cg=0.5*sqrt(g/kkw)
           do j=jstrR,jendR
            do i=istrR,iendR
             do k=1,N
              t(i,j,k,1,itemp)=18.0
              t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
             enddo
             cff1=(dble(i)+dble(iSW_corn)-(LLm/2))/pm(2,2)
             whrm(i,j)  = 0.001*2.0*exp(-env*cff1**2)       ! Hrms=2*A
             wfrq(i,j)   = sqrt(g*kkw)
             wdrx(i,j)   = 1.0                             ! cosine wave direction (xi)
             wdre(i,j)   = 0.0                             ! sine wave direction (eta)
             wdsp(i,j)   = 0.0
             sup(i,j)    = 0.0
    !         sup(i,j)    = -(0.25*whrm(i,j)**2)*kw/(2.*sinh(2*kw*h(i,j)))   ! A**2 = 0.25 Hrms**2
             ust2d(i,j)  = (0.125*whrm(i,j)**2)*sqrt(g*kkw)/h(i,j)
             ubar(i,j,1) = -((Clw**2)/((Clw**2)-(Cg**2)))*
     &                         ust2d(i,j)/h(i,j)
             zeta(i,j,1) = ubar(i,j,1)*Cg/g
             brk2dx(i,j) = 0.0
             brk2de(i,j) = 0.0
             frc2dx(i,j) = 0.0
             frc2de(i,j) = 0.0

             vbar(i,j,1)= 0.0
             vst2d(i,j) = 0.0
             do k=1,N
              u(i,j,k,1) = ubar(i,j,1)                     ! actually need
              v(i,j,k,1) = vbar(i,j,1)                     ! to init u even if baro
              ust(i,j,k) = ust2d(i,j)
              vst(i,j,k) = vst2d(i,j)
             enddo
           enddo
          enddo


        end subroutine ana_init_wec_frc_tile

       ! --------------------------------------------------------------------




        subroutine ana_wec_frc_tile (istr,iend,jstr,jend)
         ! Analytical river forcing data

         implicit none

         integer, intent(in) :: istr, iend, jstr, jend
         ! Declare looping indices
         integer :: i, j

         ! get all cppdefs flags
! # include "cppdefs.h" - included in top of module (left here to remember where they were needed)
         ! get universal scalars, g and pi
! # include "scalars.h" - included in top of module
         ! get iSW_corn variable
! # include "param.h" - included in top of module
         ! get pm coordinate transformer from grid.h
! # include "grid.h" - included in top of module
         ! note haven't included DH's forces.h as all wave parameters declared in WEC

         ! For PART 2 below: the include ocean.h files for mrl_wci part
         ! included in top of module
!# include "ocean2d.h"
!# ifdef SOLVE3D
!#   include "ocean3d.h"
!# endif

!          ! For PART 1: Declare some variables from analytical.F
!          real :: kkw, env, Cg - in top of module

          ! For PART 2: declare variables needed
          ! note - commented out SEEMINGLY unneeded variables
          real :: cff,cff1,cff2,cff3,cff4,
     &      inv_g,
     &      khd, kh,
     &      umag,
!     &      eps,
     &      wave_ramp,
!     &      wramp2,kbrk,fb,fb0,fb1,fb2,inv_fbs,
!     &      bconst,tauc,tauw,Uorb1,ka_f00,hz0,z_tide,
!     &      dtinv,inv_zb,a_brk,a_kv,c1o3,c4o3,inv_k,
     &      dd,
     &      ust_ker
!     &      intfb,
!     &      fbsrf,fn1,fn2,khmax,bz1,bz2,beta,abot,delta,
!     &      a_frc,kfrc,

          real ::  wh(PRIVATE_2D_SCRATCH_ARRAY),
     &             fr(PRIVATE_2D_SCRATCH_ARRAY),
     &             kw(PRIVATE_2D_SCRATCH_ARRAY),
!     &            brk(PRIVATE_2D_SCRATCH_ARRAY),
     &            stk(PRIVATE_2D_SCRATCH_ARRAY),
     &           Dstp(PRIVATE_2D_SCRATCH_ARRAY),
     &            act(PRIVATE_2D_SCRATCH_ARRAY),
     &             kD(PRIVATE_2D_SCRATCH_ARRAY),
     &          inv_d(PRIVATE_2D_SCRATCH_ARRAY),
     &          inv_f(PRIVATE_2D_SCRATCH_ARRAY)
!     &            frc(PRIVATE_2D_SCRATCH_ARRAY),
!     &           ebrk(PRIVATE_2D_SCRATCH_ARRAY),
!     &           erol(PRIVATE_2D_SCRATCH_ARRAY)



          real, parameter :: eps=1.e-10,
# ifdef STOKESEKMAN
     &                bconst=0.D0,
# else
     &                bconst=0.1D0, !0.03D0,               ! breaking contribution to KPP
# endif
     &                a_kv=1.2D0, !1.2D0                 ! breaking scale for eddy visc.
     &                a_brk=0.2D0,                ! breaking scale for body force
     &                a_frc=3.D0,                 ! friction scale for body force
     &                khmax=20.D0,                ! deep-water limit for k x dep
     &                c1o3=0.3333333333333333,    ! 1/3
     &                c4o3=1.3333333333333333     ! 4/3
!     &               )


# ifdef SOLVE3D
          integer  k,kk
          real     wrk1(PRIVATE_2D_SCRATCH_ARRAY,0:nz),
     &             wrk2(PRIVATE_2D_SCRATCH_ARRAY,0:nz),
!     &             wrk3(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &             wrk4(0:N),wrk5(0:N),wrk6(0:N), kvsurf
# endif


          ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
          ! This had to appear after eps=1.e-10 else wont compile
#include "compute_auxiliary_bounds.h"


          !!!!!!!!!!
          ! If first step - Cancelled as put in ana_init.F instead
          !!!!!!!!!!        This is taken from river_frc.F read

!          if (rtim(it1) < 0) then !! first time
!
!
!            call ana_init_wec_frc_tile(istr,iend,jstr,jend)
!
!
!          !!!!!!!!!!
!          ! Else timestep > 0
!          !!!!!!!!!!
!
!          else


            !*** PART 1 - from DH's WAVE_PACKET analytical.F ***!

          ! The following code is taken from analytical.F L995 from DH's
          ! WAVE_PACKET example code



          kkw=2.0*pi/1.0 ! Changed kw to kkw so as not to conflict part 2
          env=0.0001
          Cg=0.5*sqrt(g/kkw)
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              cff1=modulo(time,(gnx/Cg))        ! Changed LLm to gnx
              cff2=(dble(i)-(Cg*cff1)
     &                +dble(iSW_corn)-(gnx/2))/pm(2,2)    ! initial shift mid domain
              whrm(i,j)= 0.001*2.0*(exp(-env*cff2**2)
     &                      +exp(-env*(cff2+gnx/pm(2,2))**2))
              wfrq(i,j)= sqrt(g*2.0*pi)
              wdrx(i,j)= 1.0                             ! cosine wave direction (xi)
              wdre(i,j)= 0.0                             ! sine wave direction (eta)
              wveb(i,j)= 0.0
              wved(i,j)= 0.0
            enddo
          enddo

            !*** END PART 1 ***!

          !!!!!
          ! Taken from mrl_wci.F
          !!!!!

            !*** PART 2 - from DH's WAVE_PACKET mrl_wci.F ***!

          wave_ramp = 1.0
!          wramp2= wave_ramp**2
          inv_g = 1.0/g
!
!        ! Evaluate 2DH wave-current interaction variable.
!        ! ================================================
!        !
!        ! 2DH Stokes velocities, breaking, roller and bottom-friction
!        ! dissipation terms, defined at horizontal rho-points.

          ! Note: this would normally involve A3d and A2d scratch
          ! variables, but I have ignored them for ease now.

          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              Dstp(i,j)=max(h(i,j)+zeta(i,j,knew),eps)
              inv_d(i,j)=1.0/Dstp(i,j)
!!
              fr(i,j) =wfrq(i,j)
              wh(i,j) =wave_ramp*whrm(i,j)
              inv_f(i,j) =1.0/max(fr(i,j),eps)
              khd = Dstp(i,j)*(fr(i,j)**2)*inv_g
!#  elif defined KH_HUNT
              kh = sqrt( khd*khd + khd/(1.0 + khd*(0.6666666666
     &               +khd*(0.3555555555 + khd*(0.1608465608
     &               +khd*(0.0632098765 + khd*(0.0217540484
     &                            +khd*0.0065407983)))))) )
!#  endif
              kD(i,j) =kh
              kw(i,j) =kh*inv_d(i,j)
!!              ebrk(i,j) =wramp2*wveb(i,j)            ! epsilon_b (m3/s3)
!!              ebrk(i,j) =0            ! epsilon_b (m3/s3) kaida DD
!!
!!              cff =ebrk(i,j)
!!
!!              frc(i,j) =wramp2*wved(i,j)*kw(i,j)
!!     &                           *inv_f(i,j)     ! 2DH friction term
!!              brk(i,j) =cff*kw(i,j)*inv_f(i,j)       ! 2DH breaking term
              act(i,j) =inv_f(i,j)*(
     &                      0.125*g*(wh(i,j)**2) ! wave action density
                ! ifdef surface roller
     &                  )
!!
              stk(i,j) = 0.125*(whrm(i,j)**2)*wfrq(i,j)   ! DD
     &                      /Dstp(i,j)
!!
!!              ! DevinD - NOTE THERE is some masking stuff here I skipped
!!
            enddo
          enddo              ! <-- discard inv_d (ifndef BBL_F00)
!# ifdef BRK0
!    !      print *, 'depth-induced wave breaking deactivated'
!          do j=jstr-1,jend+1
!            do i=istr-1,iend+1
!              brk(i,j)=0.D0
!            enddo
!          enddo
!# endif

        ! brk(i,j)=0

        ! sup(i,j)=0

        ! DevinD - Dissipation terms

        !
        ! 2DH depth-averaged Stokes drift and breaking acceleration
        ! =========================================================
        !

          do j=jstrR,jendR
            do i=istrR,iendR
!              !
              ust2d(i,j)  =0.5*(stk(i-1,j)*wdrx(i-1,j)+stk(i,j)*wdrx(i,j))
!              brk2dx(i,j) =0.5*(brk(i-1,j)*wdrx(i-1,j)+brk(i,j)*wdrx(i,j))
!              frc2dx(i,j) =0.5*(frc(i-1,j)*wdrx(i-1,j)+frc(i,j)*wdrx(i,j))
!!#  ifdef MASKING
!!#   ifdef SHOREFACE
!!              ust2d(i,j)  = ust2d(i,j)*umask(i,j)
!!#   endif
!!              brk2dx(i,j) =brk2dx(i,j)*umask(i,j)
!!              frc2dx(i,j) =frc2dx(i,j)*umask(i,j)
!!#  endif
            enddo
          enddo

            ! DevinD - I don't know why vst2d is not included in loop
            ! above for efficiency, but I have kept it separate.

          do j=jstr,jendR
            do i=istrR,iendR
              vst2d(i,j)  =0.5*(stk(i,j-1)*wdre(i,j-1)+stk(i,j)*wdre(i,j))
!              brk2de(i,j) =0.5*(brk(i,j-1)*wdre(i,j-1)+brk(i,j)*wdre(i,j))
!              frc2de(i,j) =0.5*(frc(i,j-1)*wdre(i,j-1)+frc(i,j)*wdre(i,j))
!#  ifdef MASKING
!#   ifdef SHOREFACE
!              vst2d(i,j)  = vst2d(i,j)*vmask(i,j)
!#   endif
!              brk2de(i,j) =brk2de(i,j)*vmask(i,j)
!              frc2de(i,j) =frc2de(i,j)*vmask(i,j)
!#  endif
            enddo
          enddo

#ifdef SOLVE3D
!
! Conservative 3D wave-current interaction variables.
! ===================================================
!
! We estimate cell-averaged 3D Stokes drift velocity (ust_r) by
! presuming that roller Stokes drift has the same z-dependency as
! the primary waves.  Analytical volume-averaging operation is
! applied to 3D Stokes drift velocity such that:
!
!   u^st = 0.5*g*A^2/frq*k/dz/(1-exp(-4*k*D))*
!         ( exp(2*k*(z_up +h-D)) -exp(-2*k*(z_up +h+D))
!          -exp(2*k*(z_low+h-D)) +exp(-2*k*(z_low+h+D)) )
!
! in which sinh, cosh, and tanh functions in u^st are expanded to
! exponential function so as not to be singular when kD gets very
! large.  Here 0.5*g*A^2/frq is equivalent to wave action density.
! The similar rule is also adapted in computation of CalP, Kapsrf,
! and Akw terms (see the code below).
!
# define ust_r wrk1
# define inv_dz wrk2
# define inv_ex inv_f
!
! 3D depth-dependent Stokes drift velocities
! ==========================================
!

      do j=jstr-1,jend+1
        do i=istr-1,iend+1
!
          inv_ex(i,j)=1.0/(1.0-exp(-4.*kD(i,j))
     &                        -2.0*exp(-2.*kD(i,j)))
!!!!         inv_ex(i,j)=1.0/max(1.0-exp(-4.*kD(i,j))
!!!!     &                        -2.0*exp(-2.*kD(i,j)),eps)
          ust_ker=(act(i,j)*(sqrt(g*kw(i,j)*tanh(kD(i,j)))**2)/9.81)
     &                    *inv_ex(i,j)
!!!!          ust_ker=act(i,j)*kw(i,j)*inv_ex(i,j)
          cff1 = exp( 2.*kw(i,j)*(z_w(i,j,0)+h(i,j)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,0)+h(i,j)+Dstp(i,j)))
          do k=1,N,+1              ! DevinD N should be nz <-- irreversible
            inv_dz(i,j,k)=1./Hz(i,j,k)
            cff2 = exp( 2.*kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &            -exp(-2.*kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
            ust_r(i,j,k) =ust_ker*inv_dz(i,j,k)*( cff2-cff1 )
            cff1=cff2
!!
!! Horizontal interpolation of ust_r to u- & v- points.
!! *** Notice that mask should NOT be applied here.
!!
            if (i>istr-1) then
                    ust(i,j,k)=0.5*( ust_r(i-1,j,k)*wdrx(i-1,j)
     &                        +ust_r(i,j,k)*wdrx(i,j) )
            end if
            if (j>jstr-1) then
                    vst(i,j,k)=0.5*( ust_r(i,j-1,k)*wdre(i,j-1)
     &                        +ust_r(i,j,k)*wdre(i,j) )
            end if
          enddo
        enddo
      enddo



!
! calP & Kapsrf: surface pressure & Bernoulli head at rho-point
! =============================================================
!
# define kv wrk4
# define dkvdz wrk5
# define d2kv  wrk6
# define actp umag
      do j=jstr,jend        ! kv (k dot v) at rho-point
        do i=istr,iend      ! and its 1st & 2nd derivertives
          do k=1,N  ! DevinD change N to nz    ! at rho-point
            kv(k) =0.5*kw(i,j)*(
     &         wdrx(i,j)* ( u(i,j,k,nstp)+u(i+1,j,k,nstp) )
     &       + wdre(i,j)* ( v(i,j,k,nstp)+v(i,j+1,k,nstp) ) )
          enddo
          kvsurf =1.5*kv(N)-0.5*kv(N-1) ! extrapolate to surface
          do k=1,N-1
            dkvdz(k) =2.0*(kv(k+1)-kv(k))/(Hz(i,j,k+1)+Hz(i,j,k))
          enddo
          dkvdz(0) = dkvdz(1)    !2.*dkvdz(1)-dkvdz(2) ! severe!
          dkvdz(N) = dkvdz(N-1)  !2.*dkvdz(N-1)-dkvdz(N-2)
          do k=1,N
            d2kv(k) =dkvdz(k)-dkvdz(k-1) ! d^2kv/dz^2 x Hz
          enddo
          cff3 = 0.0
          do k=1,N
            dd   = z_r(i,j,k)-z_w(i,j,N)
            cff3 = cff3 + d2kv(k)*(
     &             exp( 2.*kw(i,j)*(dd-Dstp(i,j)))
     &           + exp(-2.*kw(i,j)*(dd+Dstp(i,j))) )
          enddo
          cff1 =-2.0*exp(-2.*kD(i,j))*inv_ex(i,j)*dkvdz(N)
          cff2 = dkvdz(0)/max(tanh(2.*kD(i,j)),eps)
          cff3 = cff3*inv_ex(i,j)
          cff4 =-2.0*kw(i,j)*kvsurf
!ccccc          actp =0.125*g*(wh(i,j)**2)/max(fr(i,j),eps)
!ccccc          calP(i,j) = actp*inv_g*tanh(kD(i,j))
          actp=act(i,j) ! DevinD - this is unnecessary as is define actp
          calP(i,j) = actp*tanh(kD(i,j))
     &                    *( cff1+cff2+cff3+cff4 )
!#  ifdef MASKING
!     &                                 *rmask(i,j)
!#  endif
          cff=0.0
          do k=1,N
            cff = cff +  d2kv(k)*
     &         ( exp( 2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)+Dstp(i,j))) )
          enddo
#  undef KAPSRF_SAFE
!#  ifdef KAPSRF_SAFE
!c#   define FLT 0.98
!#   define FLT 0.5
!          Kapsrf(i,j) = ( FLT*Kapsrf(i,j)  ! time-filter to suppress noises
!     &                   +(1.0-FLT)*cff*actp*inv_ex(i,j) )
!#  else
          Kapsrf(i,j) = cff*actp*inv_ex(i,j)
!#  endif
!#  ifdef MASKING
!     &                           *rmask(i,j)
!#  endif
        enddo
      enddo      ! <-- discard kv, dkvdz, d2kv, & inv_ex
# undef kv
# undef dkvdz
# undef d2kv
# undef inv_ex

# undef ust_r


#endif
! endif for SOLVE3D

!          endif ! This ends main if statement for first or other timestep

        end subroutine ana_wec_frc_tile

       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------


        subroutine read_wec_frc  !(time)
         ! Reads WEC forcing data from input file

         ! Not applicable for analytical version

        end subroutine read_wec_frc


               ! --------------------------------------------------------------------
       subroutine set_wec_frc(time)
         ! Use an analytical forcing or interpolate two instances
         ! of the forcing that are read from input files
         implicit none
         ! import/export
         real(kind=8), intent(in) :: time
         ! local
!         real    :: cf1,cf2
         real    :: tday

         tday = time/(24*3600) ! convert from time seconds to days

         if (analytical) then

!            call ana_wec_frc(tday)
        ! This is from river_frc but I don't believe it is relevant
        ! because ana_wec_frc_tile is called from set_forces.F

         else
           ! check whether we need to update forcing data
!           if (rtim(it2) < tday) then
!             call read_wec_frc(tday)
!           endif

!           ! linear interpolation in time
!           cf1 = ( rtim(it2) - tday )/(rtim(it2) - rtim(it1) )
!           cf2 = ( tday - rtim(it1) )/(rtim(it2) - rtim(it1) )
!           riv_vol = cf1*rvol(:,it1) + cf2 * rvol(:,it2)
!           riv_trc = cf1*rtrc(:,:,it1) + cf2 * rtrc(:,:,it2)
!          if (mynode==0) then
!            print *, 'volume: ',riv_vol(1:nriv_frc)
!          endif
!           ! tidal modulation of river volume
!           riv_vol = (cos(2*pi*tday) + 1)*riv_vol

           stop 'no analytical river data'

         endif
       end subroutine set_wec_frc


       ! --------------------------------------------------------------------

       subroutine find_wec_new  !(time,v_name,vt_name,ifile,irec)

        ! Get nc file model input data

       end subroutine find_wec_new

!        ! -------------------------------------
!        subroutine ana_wec_tile(istr,jstr,iend,jend)
!         !analytical wec forcing
!         ! Set everything to zero for testing
!         implicit none
!         real(kind=8), intent(in) :: istr,jstr,iend,jend
!         integer :: i, j, k
!
!          ! rho-point variables
!          do i=istr,iend
!            do j=jstr,jend
!               sup(i,j)    = 0.D0
!               brk2dx(i,j) = 0.D0
!               brk2de(i,j) = 0.D0
!               frc2dx(i,j) = 0.D0
!               frc2de(i,j) = 0.D0
!               calP(i,j)   = 0.D0
!               Kapsrf(i,j) = 0.D0
!!               lm(i,j)     = 0.D0
!               do k=1,nz
!                  brk3dx(i,j,k) = 0.D0
!                  brk3de(i,j,k) = 0.D0
!                  wst(i,j,k)    = 0.D0
!                  Kappa(i,j,k)  = 0.D0 !may not need this variable
!                  kvf(i,j,k)    = 0.D0
!                  Akw(i,j,k)    = 0.D0
!                  frc3dx(i,j,k) = 0.D0
!                  frc3de(i,j,k) = 0.D0
!                enddo
!                do k=0,nz
!                   Akb(i,j,k) = 0.D0
!                enddo
!             enddo
!          enddo
!
!          !u-point variables
!          do i=istr,iend-1
!           do j=jstr,jend
!              ust2d(i,j) = 0.D0
!              do k=1,nz
!                 ust(i,j,k) = 0.D0
!              enddo
!           enddo
!          enddo
!
!          !v-point variables
!          do i=istr,iend
!           do j=jstr,jend-1
!              vst2d(i,j) = 0.D0
!              do k=1,nz
!                 vst(i,j,k) =0.D0
!              enddo
!            enddo
!          enddo
!
!        end subroutine ana_wec_tile


        ! ----------------------------------------


!        subroutine init_wec_arrays ! DevinD added this routine.
!
!        ! Initiate arrays used by WEC (from DH's MRL init_arrays.F)
!        !      note: JM's river_frc module doesn't have an init step, it is done in main loop.
!
!            implicit none
!
!            ! beware ANA_SMFLUX and ANA_SRFLUX entry in init_arrays.F
!
!
!
!        end subroutine init_wec_arrays


        ! ----------------------------------------
                 
!        subroutine wstokes_tile(istr,jstr,iend,jend, pm, pn)
         ! Compute vertical stokes velocity         
!         implicit none
!         real(kind=8), intent(in) :: istr,jstr,iend,jend
!         integer :: i, j, k
!         real Wrk(0:nz)

!         do j=jstr,jend
!           do i=istr,iend
!              Wrk(0) = 0.D0
!              do k=1,nz,+1 ! --> recursive
!                  Wrk(k) = Wrk(k-1) -Hz(i,j,k)
!     &       *( pm(i,j)*(ust(i+1,j,k)-ust(i,j,k))
!     &         +pn(i,j)*(vst(i,j+1,k)-vst(i,j,k)) )
!             enddo 
!           wst(i,j,nz)=+0.375*Wrk(nz) +0.75*Wrk(nz-1) -0.125*Wrk(nz-2)
!           do k=nz-1,2,-1
!              wst(i,j,k)=+0.5625*(Wrk(k  )+Wrk(k-1))
!     &                 -0.0625*(Wrk(k+1)+Wrk(k-2))
!           enddo
!           wst(i,j,1)= -0.125*Wrk(2) +0.75*Wrk(1) +0.375*Wrk(0)
!          enddo
!         enddo
!# ifndef EW_PERIODIC
!         if (WESTERN_EDGE) then                       ! Set lateral
!           do k=1,nz                                   ! boundary
!             do j=jstr,jend                           ! conditions
!               wst(istr-1,j,k)=wst(istr,j,k)
!           enddo
!          enddo
!         endif
!        if (EASTERN_EDGE) then
!          do k=1,nz
!            do j=jstr,jend
!              wst(iend+1,j,k)=wst(iend,j,k)
!            enddo
!          enddo
!        endif
!# endif
!# ifndef NS_PERIODIC
!         if (SOUTHERN_EDGE) then
!           do k=1,nz
!             do i=istr,iend
!               wst(i,jstr-1,k)=wst(i,jstr,k)
!             enddo
!          enddo
!         endif
!         if (NORTHERN_EDGE) then
!           do k=1,nz
!             do i=istr,iend
!               wst(i,jend+1,k)=wst(i,jend,k)
!             enddo
!           enddo
!         endif
!#  ifndef EW_PERIODIC
!         if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
!           do k=1,nz
!             wst(istr-1,jstr-1,k)=wst(istr,jstr,k)
!           enddo
!         endif
!         if (WESTERN_EDGE .and. NORTHERN_EDGE) then
!           do k=1,nz
!             wst(istr-1,jend+1,k)=wst(istr,jend,k)
!          enddo
!         endif
!         if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
!           do k=1,nz
!             wst(iend+1, jstr-1,k)=wst(iend,jstr,k)
!           enddo
!         endif
!         if (EASTERN_EDGE .and. NORTHERN_EDGE) then
!           do k=1,nz
!             wst(iend+1,jend+1,k)=wst(iend,jend,k)
!           enddo
!         endif
!#  endif
!# endif
   
!        end subroutine wstokes_tile
        ! ----------------------------------------
        
      end module wec_frc
