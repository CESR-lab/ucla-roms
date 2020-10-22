#include "cppdefs.h"
#ifdef ANA_INITIAL
                                  ! Set initial conditions for momentum
      subroutine ana_init (tile)  ! free surface, and tracer variables
      implicit none               ! using analytical expressions.
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call ana_init_tile( istr,iend,jstr,jend,
# ifdef SOLVE3D
     &                                 A3d(1,1)
# else
     &                                 A2d(1,1)
# endif
     &                                        )
      end

      subroutine ana_init_tile (istr,iend,jstr,jend, AA)

      use tracers, only: t

      implicit none

# include "param.h"
# include "scalars.h"
      integer istr,iend,jstr,jend, i,j,k, itrc, i0,j0
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: AA

      real Ampl, fac, x,y, x0,y0, cff, cff0,cff1,cff2,cff3
# if defined COLD_FILAMENT

      real, parameter :: f0=0.8E-4,              ! Coriolis parameter
     &     h_0=60., dh=100.,  wdth=5.e+3,
     &     b_0=7.3e-3, N_0=6.3e-3, B=0.02, lam=1./50.

      real xx,zz, hbl, sgm
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: hq

#  define VISC_TERM_WIND_BALANCE
#  ifdef VISC_TERM_WIND_BALANCE
#   include "mixing.h"
      real  b11,b12,b21,b22, det, TauX,TauY
      real, dimension(0:N) :: Av,FC, c11,c12,c21,c22, d1,d2
#  endif
# elif defined ISWAKE
      integer jwake,mid_wake
      real AmpR,AmpU, cffS, Zthrm, dZthrm, xitm, alpha, csa,sna
# endif
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "eos_vars.h"
# if defined OBC_WEST  || defined OBC_EAST ||\
     defined OBC_SOUTH || defined OBC_NORTH
#  ifdef ANA_BRY
#   include "boundary.h"
#  else
#   include "climat.h"
#  endif
# endif

# include "compute_auxiliary_bounds.h"

      do j=jstrR,jendR         ! Set everything (except temperature
        do i=istrR,iendR       ! and salinity) to all-zero state, then
          zeta(i,j,1)=0.       ! modify some of the variables, if a
          ubar(i,j,1)=0.       ! non-trivial initialization required.
          vbar(i,j,1)=0.       ! Note: A code to initialize T [and S]
        enddo                  ! must always be supplied for 3D
      enddo                    ! applications.
# ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            u(i,j,k,1)=0.
            u(i,j,k,2)=0.
            v(i,j,k,1)=0.
            v(i,j,k,2)=0.
          enddo
        enddo
      enddo
# endif

# if defined SOLITON

! References for Equatorial Rossby Soliton:
!
!   Boyd, J. P., 1980:  Equatorial solitary waves. part i: Rossby
!                solitons.  J. Phys. Oceanogr., 10:1699-1717.
!
!   Boyd, J. P., 1985:  Equatorial solitary waves. part 3: Westward-
!                traveling modons.  J. Phys. Oceanogr., 15:46-54.

! For the purpose of testing the EXACT symmetry of the code there
! is an option to reverse beta-effect in "ana_grid" and, consistently,
! reverse velocities here.  This results in eastward traveling
! soliton, however it has the same energy and theoretically this
! solution should pass ETALON CHECK with exactly the same check
! points. In either case the soliton is initialized 1/3 into the
! domain in the direction of its propagation (so that it travels
! the remaining 2/3 of the domain before exiting through the
! boundary.


      if (fomn(istr,jend)-fomn(istr,jstr) < 0.) then
C$OMP MASTER
          mpi_master_only write(*,'(/1x,2A/)')  '>>> Eastward-',
     &      'traveling soliton: reversing initial velocities.'
C$OMP END MASTER
         x0=xl/3.
      else
         x0=2.*xl/3.
      endif

c**   x0=xl/2.

      y0=el/2.
      cff1=0.395
      cff2=0.771*(cff1*cff1)
      do j=jstrR,jendR
        do i=istrR,iendR
          x=xr(i,j)-x0
          y=yr(i,j)-y0
          cff3=exp(-cff1*x)
          fac=cff2*(2.*cff3/(1.+cff3*cff3))**2
          zeta(i,j,1)=0.25*fac*(6.*y*y+3.)*exp(-0.5*y*y)
        enddo
      enddo

      do j=jstrR,jendR
        do i=istr,iendR
          x=0.5*(xr(i-1,j)+xr(i,j))-x0
          y=0.5*(yr(i-1,j)+yr(i,j))-y0
          cff3=exp(-cff1*x)
          fac=cff2 * (2.*cff3/(1.+cff3*cff3))**2
          ubar(i,j,1)=0.25*fac*(6.*y*y-9.)*exp(-0.5*y*y)
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          x=0.5*(xr(i,j-1)+xr(i,j))-x0
          y=0.5*(yr(i,j-1)+yr(i,j))-y0
          cff3=exp(-cff1*x)
          fac=cff2 * (2.*cff3/(1.+cff3*cff3))**2
          vbar(i,j,1)=2.*fac*y*(-2.*cff1*tanh(cff1*x))
     &                                  *exp(-0.5*y*y)
        enddo
      enddo
      if (fomn(istr,jend)-fomn(istr,jstr) < 0.) then
        do j=jstrR,jendR               ! Reverse velocities, if the
          do i=istrR,iendR             ! sign of gradient of Coriolis
            ubar(i,j,1)=-ubar(i,j,1)   ! parameter "beta" is negative
            vbar(i,j,1)=-vbar(i,j,1)   ! (see ana_grid). In this case
          enddo                        ! soliton should travel to the
        enddo                          ! east, while maintaining all
      endif                            ! energetics exactly the same.

# elif defined WAVE_RAD

      x0=xl/2. ; y0=el/2.              ! Set initial perturbation in
      cff=64./xl                       ! in free surface and velocities
      do j=jstr,jend                   ! for isotropically expanding
        do i=istr,iend                 ! circular wave.
          x=cff*(xr(i,j)-x0)
          y=cff*(yr(i,j)-y0)
          cff1=sqrt(x*x+y*y)
          if (cff1 > 9. .and. cff1 < 11.) then
            zeta(i,j,1)=1.+cos(pi*(cff1-10.))
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          x=cff*(0.5*(xr(i-1,j)+xr(i,j))-x0)
          y=cff*(               yr(i,j) -y0)
          cff1=sqrt(x*x+y*y)
          if (cff1 > 9. .and. cff1 < 11.) then
            ubar(i,j,1)=(1.+cos(pi*(cff1-10.)))
     &      *(x/cff1)*sqrt(2.*g/(h(i-1,j)+h(i,j)))
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          x=cff*(               xr(i,j) -x0)
          y=cff*(0.5*(yr(i,j-1)+yr(i,j))-y0)
          cff1=sqrt(x*x+y*y)
          if (cff1 > 9. .and. cff1 < 11.) then
            vbar(i,j,1)=(1.+cos(pi*(cff1-10.)))
     &      *(y/cff1)*sqrt(2.*g/(h(i,j-1)+h(i,j)))
          endif
        enddo
      enddo

# elif defined DOUBLE_GYRE
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=3.8*exp(z_r(i,j,k)/800.)
     &                *(0.6-0.4*tanh(z_r(i,j,k)/800.))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined CANYON_A
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=T0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined CANYON_B
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=3.488*exp(z_r(i,j,k)/800.)
     &        *(1.-0.666666666666*tanh(z_r(i,j,k)/800.))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined DAMEE_B
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
c            t(i,j,k,1,itemp)=0.2+6.*exp(z_r(i,j,k)/2500.)
c     &              +30./(1.+exp((-200.-z_r(i,j,k))/200.))

            t(i,j,k,1,itemp)=0.2+6.*exp(z_r(i,j,k)/2500.)
     &                    +20.*exp(z_r(i,j,k)/250.)

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,isalt)=35.
#   endif
          enddo
        enddo
      enddo
#  endif
# elif defined EKMAN_SPIRAL

# elif defined GRAV_ADJ
#  ifdef SOLVE3D
      i0=LLm/2
#   ifdef MPI
      i0=i0 -iSW_corn
#   endif
      cff=Tcoef*(35.-5.)*20./(4.*rho0)   ! A provision is made to form
      do j=jstrR,jendR                   ! a 2-point transitional layer
        do i=istrR,iendR                 ! between the left- and right-
          if (i < i0) then               ! constant values in order to
            t(i,j,1,1,itemp)=5.          ! "soften" the step function.
            zeta(i,j,1)=-cff
          elseif (i == i0) then          ! Setting free surface to
            t(i,j,1,1,itemp)=5.  +3.75
            zeta(i,j,1)=-cff*0.75        !         Tcoef*DeltaT
          elseif (i == i0+1) then        ! zeta = -------------- * h
            t(i,j,1,1,itemp)=35. -3.75   !           2.*rho0
            zeta(i,j,1)=+cff*0.75
          else                           ! balances the barotropic
            t(i,j,1,1,itemp)=35.         ! mode making the solution be
            zeta(i,j,1)=+cff             ! almost free of barotropic
          endif                          ! seiche mode.
          zeta(i,j,2)=zeta(i,j,1)
          t(i,j,1,2,itemp)=t(i,j,1,1,itemp)    ! Canonically, when
        enddo                                  ! initialized to zero,
        do k=2,N                               ! free surface
          do i=istrR,iendR                     ! oscillates with its
            t(i,j,k,1,itemp)=t(i,j,1,1,itemp)  ! own frequency,
            t(i,j,k,2,itemp)=t(i,j,1,1,itemp)  ! period = 7900 sec,
          enddo                                ! amplitude=-0.06...0.06
        enddo   !<-- k                         ! provided that
      enddo   !<-- j                           ! Tcoef=0.2.
#  endif

# elif defined SLOSH
      cff=pi/xl
      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,1)=cos(cff*xr(i,j))
          zeta(i,j,2)=zeta(i,j,1)
        enddo
      enddo
   


# elif defined COLD_FILAMENT
#  ifdef VISC_TERM_WIND_BALANCE
#   define I_EXT_RANGE istrU-1,iend+1
#   define J_EXT_RANGE jstrV-1,jend+1
#  else
#   define I_EXT_RANGE istr,iend+1
#   define J_EXT_RANGE jstr,jend+1
#  endif

      x0=xl/2. !<-- sets fillament to the middle of the domain

      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          xx=0.25*(xr(i,j)+xr(i-1,j)+xr(i,j-1)+xr(i-1,j-1)) -x0

          hq(i,j)=h_0+dh/sqrt(1.+(xx/wdth)**2)    !<-- hbls
        enddo
      enddo

      cff=rho0/(g*Tcoef)  !<-- translation from boyancy to temperature
      cff0=cff*b_0
      cff1=cff*0.5*(N_0**2)*(1.+B)
      cff2=cff*0.5*(N_0**2)*(1.-B)
      cff3=cff2/lam

      do k=1,N,+1
        do j=J_EXT_RANGE
          do i=I_EXT_RANGE
            zz=0.25*( z_r(i,j,k)+z_r(i-1,j,k)+z_r(i,j-1,k)
     &                                    +z_r(i-1,j-1,k) )

            AA(i,j,k)=cff0 + cff1*zz -cff2*hq(i,j)
     &               -cff3*log(cosh(lam*(zz+hq(i,j))))
          enddo
        enddo

        do j=jstr,jend
          do i=istr,iend
            t(i,j,k,1,itemp)=0.25*( AA(i,j,k) +AA(i+1,j,k)
     &                           +AA(i,j+1,k) +AA(i+1,j+1,k)
     &                                                     )
          enddo
        enddo

        if (k == 1) then            ! The following segment converts
          do j=J_EXT_RANGE          ! "AA" from density anomaly to
            do i=I_EXT_RANGE        ! hydrostatic pressure anomaly.
              AA(i,j,0)=0.          ! <-- initialize vertical
            enddo                   !           integration
          enddo
        else
          cff=0.125*Tcoef/rho0
          do j=J_EXT_RANGE
            do i=I_EXT_RANGE
              xx=AA(i,j,0)  !<-- save (k-1)th value

              AA(i,j,0)=AA(i,j,0) +cff*(AA(i,j,k)+AA(i,j,k-1))*(
     &            z_r(i,j,k)+z_r(i-1,j,k)+z_r(i,j-1,k)+z_r(i-1,j-1,k)
     &           -z_r(i,j,k-1)-z_r(i-1,j,k-1)
     &                               -z_r(i,j-1,k-1)-z_r(i-1,j-1,k-1)
     &                                                              )
              AA(i,j,k-1)=xx
            enddo
          enddo
        endif
      enddo !<-- k

      do j=J_EXT_RANGE           ! After the above steps "AA(:,:,k)"
        do i=I_EXT_RANGE         ! becomes baroclinic pressure anomaly
          AA(i,j,N)=AA(i,j,0)    ! normalized as P/(g*rho0) and defined
        enddo                    ! at horizontal  PSI- and vertical
      enddo                      ! RHO-points.
#  undef I_EXT_RANGE
#  undef J_EXT_RANGE

#  ifdef VISC_TERM_WIND_BALANCE
      do j=jstr-1,jend           ! Compute horizontal pressure gradient
        do i=istr-1,iend         ! components and place them into d1,d2
          cff=0.5*g
          do k=1,N
            d1(k)=cff*pm(i,j)*Hz(i,j,k)*( AA(i,j  ,k)-AA(i+1,j  ,k)
     &                                   +AA(i,j+1,k)-AA(i+1,j+1,k))

            d2(k)=cff*pn(i,j)*Hz(i,j,k)*( AA(i,j  ,k)+AA(i+1,j  ,k)
     &                                   -AA(i,j+1,k)-AA(i+1,j+1,k))
          enddo

          hbl=0.25*(hq(i,j)+hq(i+1,j)+hq(i,j+1)+hq(i+1,j+1))
          cff=1./hbl

! Note that setting TauX,TauY here must match that in ana_smflux.

          TauY=0.4/rho0 ; TauX=-TauY

          cff0=hbl * 0.41*sqrt(sqrt(TauX**2 +TauY**2))

          do k=0,N
            sgm=-cff*z_w(i,j,k)
            if (sgm < 1.) then
              Av(k)=cff0*(sgm+0.1)*(1.-sgm)**2
            else
              Av(k)=0.
            endif
          enddo
          if (istrR <= i .and. i <= iendR .and.    ! Copy vertical
     &        jstrR <= j .and. j <= jendR) then    ! viscosity into
            do k=0,N                               ! shared array to
              Akv(i,j,k)=Av(k)                     ! write it out
            enddo                                  ! into history
          endif                                    ! file.

          do k=1,N-1
            FC(k)=2.*Av(k)/(Hz(i,j,k)+Hz(i,j,k+1))
          enddo
          FC(0)=0.
          FC(N)=0.

! Block-tri-diagonal problem for u,v-velocity components which are in
! simultaneous geostrophic balance with the density field initialized
! above and in viscous-Coriolis (Ekman) balance with each other.

          b11=FC(1)                     ! Free-slip bottom boundary
          b22=FC(1)                     ! condition: everything is the
          b12=-f0*Hz(i,j,1)             ! same as it would be in k-loop
          b21= f0*Hz(i,j,1)             ! for k=1 with all (k-1)-terms
                                        ! vanishing,
          det=1./(b11*b22-b12*b21)      ! FC(k-1) --> FC(0) --> 0

          c11(1)= det*b22*FC(1) ; c12(1)=-det*b12*FC(1)
          c21(1)=-det*b21*FC(1) ; c22(1)= det*b11*FC(1)

          d1(1)= det*( b22*d1(1) -b12*d2(1))
          d2(1)= det*(-b21*d1(1) +b11*d2(1))

          do k=2,N-1,+1                          !--> forward sweep
            b11=FC(k)    +FC(k-1)*(1.-c11(k-1))  ! of block-Gaussian
            b22=FC(k)    +FC(k-1)*(1.-c22(k-1))  ! elimination
            b12=-f0*Hz(i,j,k) -FC(k-1)*c12(k-1)
            b21= f0*Hz(i,j,k) -FC(k-1)*c21(k-1)

            det=1./(b11*b22-b12*b21)

            c11(k)= det*b22*FC(k) ; c12(k)=-det*b12*FC(k)
            c21(k)=-det*b21*FC(k) ; c22(k)= det*b11*FC(k)

            cff1=d1(k)+FC(k-1)*d1(k-1)
            cff2=d2(k)+FC(k-1)*d2(k-1)
            d1(k)= det*( b22*cff1 -b12*cff2)
            d2(k)= det*(-b21*cff1 +b11*cff2)
          enddo

          b11=          FC(N-1)*(1.-c11(N-1))    ! surface boundary
          b22=          FC(N-1)*(1.-c22(N-1))    ! condition: apply
          b12=-f0*Hz(i,j,N) -FC(N-1)*c12(N-1)    ! wind stress, i.e.,
          b21= f0*Hz(i,j,N) -FC(N-1)*c21(N-1)    ! replace what would
                                                 ! be upper interface
          det=1./(b11*b22-b12*b21)               ! implicit viscous
                                                 ! flux
          cff1=d1(N)+TauX +FC(N-1)*d1(N-1)       ! FC(k)*(u(k+1)-u(k))
          cff2=d2(N)+TauY +FC(N-1)*d2(N-1)       ! at k=N with TauX,y.
          d1(N)= det*( b22*cff1 -b12*cff2)
          d2(N)= det*(-b21*cff1 +b11*cff2)

          do k=N-1,1,-1                          !--> backsubstitution
            d1(k)=d1(k) +c11(k)*d1(k+1)+c12(k)*d2(k+1)
            d2(k)=d2(k) +c21(k)*d1(k+1)+c22(k)*d2(k+1)
          enddo

! After this moment d1,d2 contain horizontal velocity components u,v
! computed at common location on C-grid (at RHO-points).  Interpolate
! them to their native placements.

          if (i < iend .and. j >= jstr) then     ! Note: this code
            do k=1,N                             ! segment relies on
              u(i+1,j,k,1)=0.5*d1(k)             ! the fact that
            enddo                                ! i,j-indices are
          endif                                  ! increasing in their
          if (i >= istrU .and. j >= jstr) then   ! respective loops
            do k=1,N                             ! enclosing this part,
              u(i,j,k,1)=u(i,j,k,1)+0.5*d1(k)    ! hence for each
            enddo                                ! individual value of
          endif                                  ! "i" in computation
                                                 ! of "u" the upper
          if (i >= istr .and. j < jend) then     ! "if" statement is
            do k=1,N                             ! executed first.
              v(i,j+1,k,1)=0.5*d2(k)             ! Same applies to "j"
            enddo                                ! and "v".
          endif
          if (i >= istr .and. j >= jstrV) then
            do k=1,N
              v(i,j,k,1)=v(i,j,k,1)+0.5*d2(k)
            enddo
          endif
        enddo  !<-- i
      enddo  !<-- j

#  else
      do k=1,N                 ! Alternative to VTWB: compute velocity
       cff=0.5*g/f0            ! components from geostrophic balance at
       do j=jstr,jend          ! their native C-grid locations.
         do i=istrU,iend
           u(i,j,k,1)=-cff*(pn(i,j)+pn(i-1,j))*(AA(i,j+1,k)-AA(i,j,k))
         enddo
       enddo
       do j=jstrV,jend
         do i=istr,iend
           v(i,j,k,1)=+cff*(pm(i,j)+pm(i,j-1))*(AA(i+1,j,k)-AA(i,j,k))
         enddo
       enddo
      enddo
#  endif

      do k=1,N
        do j=jstr,jend
          do i=istrU,iend
            ubar(i,j,1)=ubar(i,j,1) +(Hz(i-1,j,k)+Hz(i,j,k))*u(i,j,k,1)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            vbar(i,j,1)=vbar(i,j,1) +(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,1)
          enddo
        enddo
      enddo  !<-- k

      do j=jstr,jend
        do i=istr,iend
          zeta(i,j,1)=0.09165 + 0.25*( AA(i,j,0)+AA(i+1,j,0)
     &                              +AA(i,j+1,0)+AA(i+1,j+1,0))
        enddo
      enddo
      do j=jstr,jend
        do i=istrU,iend
          ubar(i,j,1)=ubar(i,j,1)/( z_w(i,j,N)+z_w(i-1,j,N)
     &                             -z_w(i,j,0)-z_w(i-1,j,0))
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          vbar(i,j,1)=vbar(i,j,1)/( z_w(i,j,N)+z_w(i,j-1,N)
     &                             -z_w(i,j,0)-z_w(i,j-1,0))
        enddo
      enddo





#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          zeta(istr-1,j,1)=zeta(istr,j,1)
          ubar(istr-1,j,1)=0.
        enddo
        do k=1,N
          do j=jstr,jend
            t(istr-1,j,k,1,itemp)=t(istr,j,k,1,itemp)
            u(istrU-1,j,k,1)=0.
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          zeta(iend+1,j,1)=zeta(iend,j,1)
          ubar(iend+1,j,1)=0.
        enddo
        do k=1,N
          do j=jstr,jend
            t(iend+1,j,k,1,itemp)=t(iend,j,k,1,itemp)
            u(iend+1,j,k,1)=0.
          enddo
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          zeta(i,jstr-1,1)=zeta(i,jstr,1)
          vbar(i,jstr-1,1)=0.
        enddo
        do i=istrR,iendR
          ubar(i,jstr-1,1)=ubar(i,jstr,1)
        enddo
        do k=1,N
          do i=istr,iend
            t(i,jstr-1,k,1,itemp)=t(i,jstr,k,1,itemp)
            v(i,jstrV-1,k,1)=0.
          enddo
          do i=istrR,iendR
            u(i,jstr-1,k,1)=u(i,jstr,k,1)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          zeta(i,jend+1,1)=zeta(i,jend,1)
          vbar(i,jend+1,1)=0.
        enddo
        do i=istrR,iendR
          ubar(i,jend+1,1)=ubar(i,jend,1)
        enddo
        do k=1,N
          do i=istr,iend
            t(i,jend+1,k,1,itemp)=t(i,jend,k,1,itemp)
            v(i,jend+1,k,1)=0.
          enddo
          do i=istrR,iendR
            u(i,jend+1,k,1)=u(i,jend,k,1)
          enddo
        enddo
      endif
#  endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstrR,jendR
          vbar(istr-1,j,1)=vbar(istr,j,1)
        enddo
        do k=1,N
          do j=jstrR,jendR
            v(istr-1,j,k,1)=v(istr,j,k,1)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstrR,jendR
          vbar(iend+1,j,1)=vbar(iend,j,1)
        enddo
        do k=1,N
          do j=jstrR,jendR
            v(iend+1,j,k,1)= v(iend,j,k,1)
          enddo
        enddo
      endif
#   ifndef NS_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        zeta(istr-1,jstr-1,1)=zeta(istr,jstr,1)
        do k=1,N
          t(istr-1,jstr-1,k,1,itemp)=t(istr,jstr,k,1,itemp)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        zeta(iend+1,jstr-1,1)=zeta(iend,jstr,1)
        do k=1,N
          t(iend+1,jstr-1,k,1,itemp)=t(iend,jstr,k,1,itemp)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        zeta(istr-1,jend+1,1)=zeta(istr,jend,1)
        do k=1,N
          t(istr-1,jend+1,k,1,itemp)=t(istr,jend,k,1,itemp)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        zeta(iend+1,jend+1,1)=zeta(iend,jend,1)
        do k=1,N
          t(iend+1,jend+1,k,1,itemp)=t(iend,jend,k,1,itemp)
        enddo
      endif
#   endif
#  endif


# elif defined OVERFLOW
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
c           t(i,j,k,1,itemp)=T0-0.5*T0*(1.+tanh((yr(i,j)-25000.)
c     &                                                   /1000.))

            t(i,j,k,1,itemp)=T0*(0.5-0.5*tanh( yr(i,j)/1000.-25.))

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  endif
# elif defined SEAMOUNT
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#   ifdef NONLIN_EOS
            t(i,j,k,1,itemp)=4.0  ! Just spatially uniform T,S to
            t(i,j,k,2,isalt)=33.0 ! illustrate compressibility effect.
#   else

! There are four possible choices here: exponential profile which can
! be initialized as either a set of values at locations z_r (to test
! prsgrd31, 32, 32A); or in finite volume sense (to test prsgrd40, 41,
! 42, 44). Also there are two similar versions of linear ! profiles.

#    define AMPL 3.
#    define STRAT 500.

           t(i,j,k,1,itemp)=AMPL*exp(z_r(i,j,k)/STRAT)

c**      t(i,j,k,1,itemp)=AMPL*exp((z_r(i,j,k)-0.0003*yr(i,j))/STRAT)

c            t(i,j,k,1,itemp)=AMPL*STRAT*( exp(z_w(i,j,k)/STRAT)
c     &                                   -exp(z_w(i,j,k-1)/STRAT)
c     &                                )/(z_w(i,j,k)-z_w(i,j,k-1))

! Linear profile to check exact error cancellation.


c            t(i,j,k,1,itemp)=AMPL*(1.+z_r(i,j,k)/5000.)

c           t(i,j,k,1,itemp)=0.5*(z_w(i,j,k)+z_w(i,j,k-1))

c            t(i,j,k,1,itemp)=AMPL*(1.+0.5*(z_w(i,j,k)
c     &                            +z_w(i,j,k-1))/5000.)

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)

#    undef STRAT
#    undef AMPL
#   endif
          enddo
        enddo
      enddo
/*
#   define AREA 7.e+9
      do j=jstrR,jendR
        do i=istrR,iendR
          cff=(xr(i,j)-0.8e+5)**2 + (yr(i,j)-0.8e+5)**2
          if (cff < AREA) then
            zeta(i,j,1)=1.0*( 1.-cff/AREA )**2
          endif

          cff=(xr(i,j)-2.4e+5)**2 + (yr(i,j)-2.4e+5)**2
          if (cff < AREA) then
            zeta(i,j,1)=1.0*( 1.-cff/AREA )**2
          endif

          cff=(xr(i,j)-0.8e+5)**2 + (yr(i,j)-2.4e+5)**2
          if (cff < AREA) then
            zeta(i,j,1)=-1.0*( 1.-cff/AREA )**2
          endif

          cff=(xr(i,j)-2.4e+5)**2 + (yr(i,j)-0.8e+5)**2
          if (cff < AREA) then
            zeta(i,j,1)=-1.0*( 1.-cff/AREA )**2
          endif
        enddo
      enddo
*/

#  endif
# elif defined ISWAKE
      y0=0. !<-- here y-coordinate is centered at the channel mid axis
#  ifdef SOLVE3D

! Analytical initial conditions for 3D Island Wake problem:
!-------------------- ------- ---------- --- -- ------ ---- --------
! free-surface, velocity and temperature (actually negative density --
! it is assumed that thermal expansion coefficient is equal to 1 for
! this problem) fields are set to satisfy thermal wind balance
! equation
!                      du      g       d rho
!                 f * ---- = ------ * -------     [here rho=-T]
!                      dz     rho0       dy
!
! subject to the conditions that surface velocity is in geostrophic
! balance with free surface, and velocity vanishes at infinitely large
! negative depth.  Both velocity and temperature field are set via
! the same hyperbolic shape function
!
!                  1     1         z-Zthrm(y)
!                 --- + --- * tanh ----------
!                  2     2           dZthrm
!
! where "Zthrm" and "dZthrm" set thermocline depth and sharpness
! respectively, and y-dependency of "Zthrm" creates thermocline slope
! needed for thermal balance; free surface is determined by vertical
! integral of temperature profile and the assumption that total
! pressure (baroclinic+barotropic) gradient vanishes at infinite
! negative depth. Integration is done analytically.

      AmpR=3.0           ! magnitude of density perturbation [kg/m^3]
      AmpU=0.2           ! velocity magnitude above thermocline [m/s]
      Zthrm=-120.        ! median depth of thermocline [m]
      dZthrm=60.         ! thickness of thermocline [m]

                               ! set thermocline slope parameter to
      cffS=AmpU*rho0/(AmpR*g)  ! achieve geostrophic balance between
                               ! XI-component of baroclinic velocity
                               ! and density field.

      cff=1./(1.+exp(-2.*Zthrm/dZthrm))
      do j=jstrR,jendR
        do i=istrR,iendR
          xitm=(-Zthrm  -cffS*f(i,j)*(yr(i,j)-y0))/dZthrm
          zeta(i,j,1)=(AmpR/rho0) *0.5*dZthrm
     &                       *log( cff*(1.+exp(2.*xitm)) )
#   ifdef MASKING
     &                                 *rmask(i,j)
#   endif
          ubar(i,j,1)=0.
          vbar(i,j,1)=0.
        enddo
      enddo
      do k=N,1,-1
        do j=jstrR,jendR
          do i=istrR,iendR
            xitm=(z_r(i,j,k)-Zthrm -cffS*f(i,j)*(yr(i,j)-y0))/dZthrm
            t(i,j,k,1,itemp)= 0.5*AmpR*tanh(xitm)
#   ifdef MASKING
     &                                 *rmask(i,j)
#   endif
          enddo
        enddo
        do j=jstrR,jendR
          do i=istr,IendR
            xitm=( 0.5*(z_r(i,j,k)+z_r(i-1,j,k)) -Zthrm
     &                    -cffS*0.5*(f(i,j)+f(i-1,j))
     &                  *(0.5*(yr(i,j)+yr(i-1,j))-y0) )/dZthrm

            u(i,j,k,1)=AmpU*(0.5 + 0.5*tanh(xitm))
#   ifdef MASKING
     &                                  *umask(i,j)
#   endif
            ubar(i,j,1)=ubar(i,j,1) +(Hz(i-1,j,k)+Hz(i,j,k))
     &                                           *u(i,j,k,1)
          enddo
        enddo
      enddo
      do j=jstrR,jendR
        do i=istr,IendR
          ubar(i,j,1)=ubar(i,j,1)/( z_w(i,j,N)+z_w(i-1,j,N)
     &                             -z_w(i,j,0)-z_w(i-1,j,0))
        enddo
      enddo

#   define SIGNELL
#   ifdef SIGNELL

      AmpU=0.15           ! This setup is originally from Rich Signell,
      cff=AmpU/g          ! who chose to pose a non-stratified problem
                          ! to be solved by full 3D code.
#    ifdef MPI
      j0=MMm/2 -jSW_corn
#    else
      j0=MMm/2
#    endif

      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,1)=-cff*f(i,j)*(yr(i,j)-y0)
          ubar(i,j,1)=AmpU
          vbar(i,j,1)=0.
        enddo
      enddo
      do k=1,N
        do j=jstrR,jendR                ! Temperature field just below
          do i=istrR,iendR              ! is passive tracer initialized
            u(i,j,k,1)=ubar(i,j,1)      ! as a "smoothed" step function
            v(i,j,k,1)=vbar(i,j,1)      ! for the flow visualization
                                        ! purpose. It does not affect
            if (j < j0) then            ! the dynamics of the flow.
              t(i,j,k,1,itemp)=-1.
            elseif (j == j0) then
              t(i,j,k,1,itemp)=-0.333333333333333
            elseif (j == j0+1) then
              t(i,j,k,1,itemp)=+0.333333333333333
            else
              t(i,j,k,1,itemp)=+1.
            endif

          enddo
        enddo
      enddo
#   endif
#  else
      AmpU=0.2                ! 2D setup for flow around cylinder
      cff=AmpU/g              ! problem: set u-velocity constant and
      alpha=0.3               ! free surface from geostrophic balance.
      csa=cos(alpha)
      sna=sin(alpha)

      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,1)=-cff*f(i,j)*((yr(i,j)-y0)*csa-(xr(i,j)-y0)*sna)
          ubar(i,j,1)=AmpU*csa
          vbar(i,j,1)=AmpU*sna
        enddo
      enddo
#  endif  /* SOLVE3D */

# elif defined SHELFRONT
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=T0+2.5*tanh((yr(i,j)-50000.0)/20000.0)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,1,isalt)=S0  +  tanh((yr(i,j)-50000.0)/20000.0)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif
# elif defined UPWELLING
#  ifdef SOLVE3D
      j0=MMm/2
#   ifdef MPI
      j0=j0 -jSW_corn
#   endif

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR

c**         t(i,j,k,1,itemp)=14.+8.*exp(z_r(i,j,k)/50.)


#   define Z0 (-35.)
#   define THKNSS 6.5
#   define Z1 (-75.)
#   define STRAT 150.

#   define THKNSS 25.0
#   define Z1 (-120.)
#   define STRAT 75.

c**         t(i,j,k,1,itemp)=14.+4.*tanh((z_r(i,j,k)-Z0)/THKNSS)
c**  &                                     +(z_r(i,j,k)-Z1)/STRAT

            t(i,j,k,1,itemp)=14.+4.*THKNSS*log(
     &                             cosh((z_w(i,j,k  )-Z0)/THKNSS)
     &                            /cosh((z_w(i,j,k-1)-Z0)/THKNSS)
     &                                )/(z_w(i,j,k)-z_w(i,j,k-1))
     &                  +((z_w(i,j,k)+z_w(i,j,k-1))/2. -Z1)/STRAT

#   undef Z0
#   undef THKNSS
#   undef Z1
#   undef STRAT

c            t(i,j,k,1,itemp)=14.

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)

#   ifdef SALINITY
c*            t(i,j,k,1,isalt)=1.E-4*yr(i,j)-S0 +0.1
            t(i,j,k,1,isalt)=1.

            if (j < j0) then           ! A provision is made to form
              t(i,j,k,1,isalt)=0.      ! a 2-point transitional layer
            elseif (j == j0) then      ! between S=0 and S=+1, hence
              t(i,j,k,1,isalt)=0.2     ! "softening" the step function.
            elseif (j == j0+1) then
              t(i,j,k,1,isalt)=0.8
            else
              t(i,j,k,1,isalt)=1.
            endif

            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
          enddo
        enddo
      enddo
#  endif
# elif defined  RIVER_SOURCE
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=4.+10.*exp(z_r(i,j,k)/50.)
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif
# elif defined PACIFIC_2D
      Ampl=1.0
      cff1=800./dble(LLm)

      if (LLm==840 .and. MMm==512) then
        i0=519
        j0=283
      elseif (LLm==1680 .and. MMm==1024) then
        i0=1036
        j0=565
      else
        i0=LLm/2
        j0=MMm/2
      endif
#  ifdef MPI
      i0=i0 -iSW_corn
      j0=j0 -jSW_corn
#  endif

      cff1=800./dble(LLm)
      do j=jstrR,jendR
        do i=istrR,iendR
          cff=cff1*sqrt(dble(i-i0)**2+dble(j-j0)**2)
          if (cff<pi) then
            zeta(i,j,1)=Ampl*(0.5+0.5*cos(cff))
          endif
        enddo
      enddo
# elif defined USWEST
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
c**         t(i,j,k,1,itemp)=10.+10.*exp(z_r(i,j,k)/200.)
c**         t(i,j,k,1,itemp)=20. + z_r(i,j,k)/4500.
c*          t(i,j,k,1,itemp)=4.+13.*exp(z_r(i,j,k)/450.)


c           t(i,j,k,1,itemp)=4.+13.*450.*( exp(z_w(i,j,k)/450.)
c     &                                   -exp(z_w(i,j,k-1)/450.))
c     &                                 /(z_w(i,j,k)-z_w(i,j,k-1))


#   define Z0 (-80.)
#   define THKNSS 50.
#   define Z1 0.
#   define STRAT 1000.

c*         t(i,j,k,1,itemp)=4.+4.*tanh((z_r(i,j,k)-Z0)/THKNSS)
c*     &                                  +(z_r(i,j,k)-Z1)/STRAT

           t(i,j,k,1,itemp)=14.+5.*THKNSS*log(
     &                            cosh((z_w(i,j,k )-Z0)/THKNSS)
     &                           /cosh((z_w(i,j,k-1)-Z0)/THKNSS)
     &                               )/(z_w(i,j,k)-z_w(i,j,k-1))

     &               +((z_w(i,j,k)+z_w(i,j,k-1))/2.-Z1)/STRAT
#   undef Z0
#   undef THKNSS
#   undef Z1
#   undef STRAT

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)

#   ifdef SALINITY
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
#   ifdef BIOLOGY
            t(i,j,k,1,ino3_)=0.
            t(i,j,k,2,ino3_)=t(i,j,k,1,ino3_)
            t(i,j,k,1,inh4_)=0.
            t(i,j,k,2,inh4_)=t(i,j,k,1,inh4_)
            t(i,j,k,1,idet_)=0.
            t(i,j,k,2,idet_)=t(i,j,k,1,idet_)
            t(i,j,k,1,iphyt)=0.
            t(i,j,k,2,iphyt)=t(i,j,k,1,iphyt)
            t(i,j,k,1,izoo_)=0.
            t(i,j,k,2,izoo_)=t(i,j,k,1,izoo_)
#   endif
          enddo
        enddo
      enddo
#  endif
# elif defined BALTIC
#  ifdef SOLVE3D
      i0=116
#   ifdef MPI
      i0=i0 -iSW_corn
#   endif

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            if (i<i0) then

              t(i,j,k,1,itemp)=8.+8.*exp(z_r(i,j,k)/40.) !   Baltic sea code must go here.
              t(i,j,k,1,isalt)=16.
            else

            endif




            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif

# else
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=???
            t(i,j,k,2,itemp)=???
            t(i,j,k,1,isalt)=???
            t(i,j,k,2,isalt)=???
          enddo
        enddo
      enddo
#  endif
# endif  /* end of configuration selection switch */

! The rest of the code is independent of the specific configuration.

! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.


# if defined OBC_WEST  || defined OBC_EAST ||\
     defined OBC_SOUTH || defined OBC_NORTH
#  ifdef ANA_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR                     ! Save perimeter values
          zeta_west(j)=zeta(istr-1,j,1)      ! of the initial fields
          ubar_west(j)=ubar(istr  ,j,1)      ! to be used as external
          vbar_west(j)=vbar(istr-1,j,1)      ! forcing data for side
        enddo                                ! boundaries.
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_EAST
      if (EASTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR
          zeta_east(j)=zeta(iend+1,j,1)
          ubar_east(j)=ubar(iend+1,j,1)
          vbar_east(j)=vbar(iend+1,j,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_east(j,k)=u(iend+1,j,k,1)
            v_east(j,k)=v(iend+1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_east(j,k,itemp)=t(iend+1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif


#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_south(i)=zeta(i,jstr-1,1)
          ubar_south(i)=ubar(i,jstr-1,1)
          vbar_south(i)=vbar(i,jstr  ,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,jstr-1,k,1)
            v_south(i,k)=v(i,jstr  ,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,itemp)=t(i,jstr-1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_north(i)=zeta(i,jend+1,1)
          ubar_north(i)=ubar(i,jend+1,1)
          vbar_north(i)=vbar(i,jend+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,jend+1,k,1)
            v_north(i,k)=v(i,jend+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,itemp)=t(i,jend+1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  else /* alternative to ANA_BRY */
#   ifdef UCLIMATOLOGY
      do j=jstrR,jendR                   ! Save initial data into
        do i=istrR,iendR                 ! climatology for the
          ubclm(i,j)=ubar(i,j,1)         ! subsequent use as inflow
          vbclm(i,j)=vbar(i,j,1)         ! boundary conditions
        enddo
      enddo
#   endif
#   ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#    ifdef UCLIMATOLOGY
            uclm(i,j,k)=u(i,j,k,1)
            vclm(i,j,k)=v(i,j,k,1)
#    endif
#    ifdef TCLIMATOLOGY
            tclm(i,j,k,itemp)=t(i,j,k,1,itemp)
#    endif
          enddo
        enddo
      enddo
#   endif /* SOLVE3D */
#  endif /* ANA_BRY */
# endif /* at least one of OBC_XXXXX */





# ifdef EXCHANGE
      call exchange2d_3_tile(istr,iend,jstr,jend,
     &                       zeta(START_2D_ARRAY,1),
     &                       ubar(START_2D_ARRAY,1),
     &                       vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
      call exchange_2_tile(istr,iend,jstr,jend,
     &                     u(START_2D_ARRAY,1,1), N,
     &                     v(START_2D_ARRAY,1,1), N)
      do itrc=1,NT
        call exchange_tile(istr,iend,jstr,jend,
     &                t(START_2D_ARRAY,1,1,itrc), N)
      enddo
#  endif
#  ifdef ISWAKE
#   ifdef UCLIMATOLOGY
      call exchange2d_2_tile(istr,iend,jstr,jend, ubclm,vbclm)
#   endif
#   ifdef SOLVE3D
#    ifdef UCLIMATOLOGY
      call exchange_2_tile(istr,iend,jstr,jend, uclm,N, vclm,N)
#    endif
#    ifdef TCLIMATOLOGY
      call exchange_tile(istr,iend,jstr,jend,
     &                    tclm(START_2D_ARRAY,1,itemp), N)
#    endif
#   endif
#  elif defined COLD_FILAMENT
      call exchange_tile(istr,iend,jstr,jend, Akv,N+1)

#  endif
# endif  /* EXCHANGE */


# if defined OBC_WEST && (defined NS_PERIODIC || defined MPI)
      if (WESTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_west,1)
        call exch_bry_EW_tile (jstr,jend, ubar_west,1)
        call exch_bry_EW_tile (jstr,jend, vbar_west,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_west,N)
        call exch_bry_EW_tile (jstr,jend, v_west,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_west(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_EAST && (defined NS_PERIODIC || defined MPI)
      if (EASTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_east,1)
        call exch_bry_EW_tile (jstr,jend, ubar_east,1)
        call exch_bry_EW_tile (jstr,jend, vbar_east,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_east,N)
        call exch_bry_EW_tile (jstr,jend, v_east,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_east(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_SOUTH && (defined EW_PERIODIC || defined MPI)
      if (SOUTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_south,1)
        call exch_bry_NS_tile (istr,iend, ubar_south,1)
        call exch_bry_NS_tile (istr,iend, vbar_south,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_south,N)
        call exch_bry_NS_tile (istr,iend, v_south,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_south(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_NORTH && (defined EW_PERIODIC || defined MPI)
      if (NORTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_north,1)
        call exch_bry_NS_tile (istr,iend, ubar_north,1)
        call exch_bry_NS_tile (istr,iend, vbar_north,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_north,N)
        call exch_bry_NS_tile (istr,iend, v_north,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_north(0,1,itrc),N)
        enddo
#  endif
      endif
# endif
      end
#else
      subroutine ana_init_empty
      end
#endif /* ANA_INITIAL */
