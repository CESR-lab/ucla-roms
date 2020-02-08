#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd44_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                      A3d(1,3), A3d(1,4), A3d(1,5), A3d(1,6),
     &                                A2d(1,1), A2d(1,2), A2d(1,3),
     &                                A2d(1,4), A2d(1,5), A2d(1,6))
      end

      subroutine prsgrd44_tile (istr,iend,jstr,jend, ru,rv, P,r,d,
     &                                  FX, FC,r1, aR,dR,aL,dL)

      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv,FX
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: P, r, d
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: FC, r1,
     &                                           aR,dR, aL,dL
      real dh,dP, cff,cff1,cffR,cffL, deltaR,deltaL, Ampl,Hdd,rr,limtr
      real, parameter :: Zero=0., Half=0.5, One=1., ThreeHalfth=1.5,
     &                   Two=2.,  Three=3., Four=4., OneSixth=1./6.,
     &                   OneFifth=0.2, OneTwelfth=1./12., eps=1.E-33

# include "grid.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "scalars.h"

# include "compute_auxiliary_bounds.h"

c--#define NEUMANN

! Parabolic WENO reconstruction: The second and third loops below
!---------- ---- --------------- compute left and right side limits
! aL,aR of the field rho assuming monotonized parabolic distributions
! within each grid box. Also computed are dL,dR, which are then used
! as a measure of quadratic variation during sabsequent WENO
! reconciliation of side limits.

      do j=jstrV-1,jend
        do k=N-1,1,-1
          do i=istrU-1,iend
            FC(i,k)=One/(Hz(i,j,k+1)+Hz(i,j,k))
            r(i,j,k)=FC(i,k)*(rho(i,j,k+1)*Hz(i,j,k)
     &                        +rho(i,j,k)*Hz(i,j,k+1))
            d(i,j,k)=FC(i,k)*(rho(i,j,k+1)-rho(i,j,k))
          enddo
        enddo

        do k=2,N-1
          do i=istrU-1,iend
            deltaR=Hz(i,j,k)*d(i,j,k)
            deltaL=Hz(i,j,k)*d(i,j,k-1)

            if (deltaR*deltaL < Zero) then
              deltaR=Zero
              deltaL=Zero
            endif
            cff=Hz(i,j,k-1)+Two*Hz(i,j,k)+Hz(i,j,k+1)
            cffR=cff*d(i,j,k)
            cffL=cff*d(i,j,k-1)
            if (abs(deltaR) > abs(cffL)) deltaR=cffL
            if (abs(deltaL) > abs(cffR)) deltaL=cffR

            cff=(deltaR-deltaL)/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
            deltaR=deltaR-cff*Hz(i,j,k+1)
            deltaL=deltaL+cff*Hz(i,j,k-1)

            aR(i,k)=rho(i,j,k)+deltaR
            aL(i,k)=rho(i,j,k)-deltaL

            dR(i,k)=( Two*deltaR-deltaL )**2
            dL(i,k)=( Two*deltaL-deltaR )**2
          enddo
        enddo

        do i=istrU-1,iend
          aL(i,N)=aR(i,N-1)
          aR(i,N)=Two*rho(i,j,N)-aL(i,N)

          dR(i,N)=( Two*aR(i,N)+aL(i,N)-Three*rho(i,j,N) )**2
          dL(i,N)=( Three*rho(i,j,N)-Two*aL(i,N)-aR(i,N) )**2

          aR(i,1)=aL(i,2)
          aL(i,1)=Two*rho(i,j,1)-aR(i,1)

          dR(i,1)=( Two*aR(i,1)+aL(i,1)-Three*rho(i,j,1) )**2
          dL(i,1)=( Three*rho(i,j,1)-Two*aL(i,1)-aR(i,1) )**2
        enddo

        do k=1,N-1
          do i=istrU-1,iend
             deltaL=max(dL(i,k),eps)
             deltaR=max(dR(i,k+1),eps)
             r1(i,k)=(deltaR*aR(i,k)+deltaL*aL(i,k+1))
     &                                /(deltaR+deltaL)
          enddo
        enddo

        do i=istrU-1,iend
#ifdef NEUMANN
          r1(i,N)=ThreeHalfth*rho(i,j,N)-Half*r1(i,N-1)
          r1(i,0)=ThreeHalfth*rho(i,j,1)-Half*r1(i,1)
#else
          r1(i,N)=Two*rho(i,j,N)-r1(i,N-1)
          r1(i,0)=Two*rho(i,j,1)-r1(i,1)
#endif
        enddo

! Power-law reconciliation: It starts with computation of side
!------ --- --------------- limits dR,dL of the first derivative
! assuming parabolic distributions within each grid box. In this
! version of the code, before doing so (see "else" branch of 3-way
! switch below), in situation when interfacial deviations deltaR
! and deltaL are differ by more than a factor of two (hence
! monotonic parabolic fit becomes impossible), the parabolic
! assumption is switched to power-law function,  such that its
! derivative is zero at one end and, consequently, larger than
! taht of (would be) limited parabolic on the other end. The basic
! parabolic version of the code is commented out, but left here
! for reference.

        do k=1,N
          do i=istrU-1,iend

c**         cff=Two/Hz(i,j,k)
c**         dR(i,k)=cff*(Two*r1(i,k)+r1(i,k-1)-Three*rho(i,j,k))
c**         dL(i,k)=cff*(Three*rho(i,j,k)-Two*r1(i,k-1)-r1(i,k))
c**         cff=r(i,j,k)-r(i,j,k-1)
c**         if (cff*dR(i,k)<Zero) dR(i,k)=Zero
c**         if (cff*dL(i,k)<Zero) dL(i,k)=Zero

            deltaR=r1(i,k)-rho(i,j,k)
            deltaL=rho(i,j,k)-r1(i,k-1)
            cff=deltaR*deltaL
            if (cff>eps) then
              cff=(deltaR+deltaL)/cff
            else
              cff=Zero
            endif
            cffL=cff*deltaL
            cffR=cff*deltaR
            if (cffL>3.) then
              cffL=cffL*deltaL
              cffR=Zero
            elseif (cffR>3.) then
              cffL=Zero
              cffR=cffR*deltaR
            else
              cffL=Four*deltaL-Two*deltaR
              cffR=Four*deltaR-Two*deltaL
            endif
            cff=One/Hz(i,j,k)
            dR(i,k)=cff*cffR
            dL(i,k)=cff*cffL
          enddo
        enddo

! Compute final value of derivative at each interface by reconciling
! two side limits dR(k) and dL(k+1) coming from adjacent grid boxes.
! the difference between these two also causes change of interfacial
! value r(k) by Ampl. The commented code (left here for reference)
! computes the exact value of Ampl assuming power law reconciliation
! and solving associated quadratic equation. The code segment below
! it corresponds to Pade fit to exact solution, which avoids
! computation of sqrt for the sake of computational efficiency.

        do k=N-1,1,-1
          do i=istrU-1,iend
            d(i,j,k)=FC(i,k)*(Hz(i,j,k+1)*dL(i,k+1)+Hz(i,j,k)*dR(i,k))

            cffR=8.*(dR(i,k)+Two*dL(i,k))
            cffL=8.*(dL(i,k+1)+Two*dR(i,k+1))
            if (abs(d(i,j,k))>abs(cffR)) d(i,j,k)=cffR
            if (abs(d(i,j,k))>abs(cffL)) d(i,j,k)=cffL

            if ((dL(i,k+1)-dR(i,k))*(rho(i,j,k+1)
     &               -rho(i,j,k)) > Zero) then
              Hdd=Hz(i,j,k)*(d(i,j,k)-dR(i,k))
              rr=rho(i,j,k)-r1(i,k-1)
            else
              Hdd=Hz(i,j,k+1)*(dL(i,k+1)-d(i,j,k))
              rr=r1(i,k+1)-rho(i,j,k+1)
            endif
            rr=abs(rr)

c*             Ampl=TwoFifth*Hdd*rr
c*             Hdd=abs(Hdd)
c*             cff=rr*(rr+0.16*Hdd)
c*             if (cff>eps) Ampl=Ampl/(rr+sqrt(cff))

            Ampl=OneFifth*Hdd*rr
            Hdd=abs(Hdd)
            cff=rr*rr+0.0763636363636363636*Hdd
     &                *(rr+0.004329004329004329*Hdd)
            if (cff>eps) then
              Ampl=Ampl*(rr+0.0363636363636363636*Hdd)/cff
            else
              Ampl=Zero
            endif

            r(i,j,k)=r1(i,k)+Ampl
          enddo
        enddo
        do i=istrU-1,iend
#ifdef NEUMANN
          r(i,j,0)=ThreeHalfs*rho(i,j,1)-Half*r(i,j,1)
          r(i,j,N)=ThreeHalfs*rho(i,j,N)-Half*r(i,j,N-1)
          d(i,j,0)=Zero
          d(i,j,N)=Zero
#else
          r(i,j,0)=Two*rho(i,j,1)-r(i,j,1)
          r(i,j,N)=Two*rho(i,j,N)-r(i,j,N-1)
          d(i,j,0)=d(i,j,1)
          d(i,j,N)=d(i,j,N-1)
#endif
        enddo

! (2) Compute pressure P and lateral pressure force FX:

        do i=istrU-1,iend  !--> set pressure at free surface,
          P(i,j,N)=Zero    !    initialize vertical integration
        enddo
        do k=N,1,-1
          do i=istrU-1,iend
            P(i,j,k-1)=P(i,j,k)+Hz(i,j,k)*rho(i,j,k)

            FX(i,j,k)=Half*Hz(i,j,k)*( P(i,j,k)+P(i,j,k-1)
     &            +OneFifth*Hz(i,j,k)*( r(i,j,k)-r(i,j,k-1)
     &           -OneTwelfth*Hz(i,j,k)*(d(i,j,k)+d(i,j,k-1)
     &                                                   ))
     &                                                    )
          enddo
        enddo

! (3) Compute net pressure gradient forces, XI- and ETA-components.

        if (j>=jstr) then
          do i=istrU,iend
            FC(i,N)=Zero    !--> set pressure at free surface.
          enddo
          cff=0.5*g
          cff1=g/rho0
          do k=N,1,-1
            do i=istrU,iend
              dh=z_w(i,j,k-1)-z_w(i-1,j,k-1)

              dP=P(i-1,j,k-1)-P(i,j,k-1)
              rr=Half*dh*( r(i,j,k-1)+r(i-1,j,k-1)
     &          -OneSixth*dh*(d(i,j,k-1)-d(i-1,j,k-1)))

              limtr=Two*dP*rr
              rr=rr*rr + dP*dP
              if (limtr>eps*rr) then
                limtr=limtr/rr
              else
                limtr=0.
              endif

              FC(i,k-1)=Half*dh*( P(i,j,k-1)+P(i-1,j,k-1)
     &        +limtr*OneFifth*dh*( r(i,j,k-1)-r(i-1,j,k-1)
     &             -OneTwelfth*dh*(d(i,j,k-1)+d(i-1,j,k-1)
     &                                                  ))
     &                                                   )
              ru(i,j,k)=( cff*(Hz(i-1,j,k)+Hz(i,j,k))
     &                          *(z_w(i-1,j,N)-z_w(i,j,N))

     &                      +cff1*( FX(i-1,j,k)-FX(i,j,k)
     &                                 +FC(i,k)-FC(i,k-1))
     &                                         )*dn_u(i,j)
            enddo
          enddo
        endif

        if (j>=jstrV) then
          do i=istr,iend
            FC(i,N)=Zero    !--> set pressure at free surface.
          enddo
          cff=0.5*g
          cff1=g/rho0
          do k=N,1,-1
            do i=istr,iend
              dh=z_w(i,j,k-1)-z_w(i,j-1,k-1)

              dP=P(i,j-1,k-1)-P(i,j,k-1)
              rr=Half*dh*( r(i,j,k-1)+r(i,j-1,k-1)
     &          -OneSixth*dh*(d(i,j,k-1)-d(i,j-1,k-1)))

              limtr=Two*dP*rr
              rr=rr*rr + dP*dP
              if (limtr>eps*rr) then
                limtr=limtr/rr
              else
                limtr=0.
              endif

              FC(i,k-1)=Half*dh*( P(i,j,k-1)+P(i,j-1,k-1)
     &        +limtr*OneFifth*dh*( r(i,j,k-1)-r(i,j-1,k-1)
     &             -OneTwelfth*dh*(d(i,j,k-1)+d(i,j-1,k-1)
     &                                                  ))
     &                                                   )
              rv(i,j,k)=( cff*(Hz(i,j-1,k)+Hz(i,j,k))
     &                          *(z_w(i,j-1,N)-z_w(i,j,N))

     &                      +cff1*( FX(i,j-1,k)-FX(i,j,k)
     &                                 +FC(i,k)-FC(i,k-1))
     &                                         )*dm_v(i,j)
            enddo
          enddo
        endif
      enddo    !<-- j
      end

#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */

