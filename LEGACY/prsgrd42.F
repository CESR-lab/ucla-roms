#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd42_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                      A3d(1,3), A3d(1,4), A3d(1,5), A2d(1,1),
     &                      A2d(1,2), A2d(1,3), A2d(1,4), A2d(1,5))
      end

      subroutine prsgrd42_tile (istr,iend,jstr,jend, ru,rv, P,r,FX,
     &                                              FC,aR,dR,aL,dL)

      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: P,r,FX
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: FC,aL,aR,dL,dR

      real dh,dP,rr, cff,cff1, cffL,cffR, deltaL,deltaR
      real, parameter :: Zero=0., Half=0.5, One=1., ThreeHalfth=1.5,
     &               Two=2., Three=3., OneSixth=1./6., epsil=1.E-33

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

      do j=jstrV-2,jend+1
        do k=N-1,1,-1
          do i=istrU-2,iend+1
            FC(i,k)=(rho(i,j,k+1)-rho(i,j,k))/(Hz(i,j,k+1)+Hz(i,j,k))
          enddo
        enddo

        do k=2,N-1
          do i=istrU-2,iend+1
            deltaR=Hz(i,j,k)*FC(i,k)
            deltaL=Hz(i,j,k)*FC(i,k-1)

            if (deltaR*deltaL < Zero) then
              deltaR=Zero
              deltaL=Zero
            endif
            cff=Hz(i,j,k-1)+Two*Hz(i,j,k)+Hz(i,j,k+1)
            cffR=cff*FC(i,k)
            cffL=cff*FC(i,k-1)
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

        do i=istrU-2,iend+1
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
          do i=istrU-2,iend+1
             deltaL=max(dL(i,k),epsil)
             deltaR=max(dR(i,k+1),epsil)
             r(i,j,k)=(deltaR*aR(i,k)+deltaL*aL(i,k+1))
     &                                 /(deltaR+deltaL)
          enddo
        enddo

        do i=istrU-2,iend+1
#ifdef NEUMANN
          r(i,j,N)=ThreeHalfth*rho(i,j,N)-Half*r(i,j,N-1)
          r(i,j,0)=ThreeHalfth*rho(i,j,1)-Half*r(i,j,1)
#else
          r(i,j,N)=Two*rho(i,j,N)-r(i,j,N-1)
          r(i,j,0)=Two*rho(i,j,1)-r(i,j,1)
#endif
        enddo
!
! (2) Compute pressure P and lateral pressure force FX:
!
        do i=istrU-2,iend+1 !--> set pressure at free surface,
          P(i,j,N)=Zero     !    initialize vertical integration
        enddo
        do k=N,1,-1
          do i=istrU-2,iend+1
            P(i,j,k-1)=P(i,j,k)+Hz(i,j,k)*rho(i,j,k)

            deltaR=r(i,j,k)-rho(i,j,k)
            deltaL=rho(i,j,k)-r(i,j,k-1)
            if (deltaR*deltaL < Zero) then
              rr=Zero
            elseif (abs(deltaR) > Two*abs(deltaL)) then
              rr=Three*deltaL
            elseif (abs(deltaL) > Two*abs(deltaR)) then
              rr=Three*deltaR
            else
              rr=deltaR+deltaL
            endif

            FX(i,j,k)=Half*Hz(i,j,k)*( P(i,j,k)+P(i,j,k-1)
     &                              +OneSixth*rr*Hz(i,j,k))
          enddo
        enddo

! (3) Compute net horizontal (along geopetential) XI- and
! ETA-components of pressure gradients

        if (j>=jstr .and. j<=jend) then
          do i=istrU-1,iend+1
            FC(i,N)=Zero      !--> set pressure at free surface.
          enddo
          do k=N,1,-1
            do i=istrU-1,iend+1
              dP=P(i-1,j,k-1)-P(i,j,k-1)
              dh=z_w(i,j,k-1)-z_w(i-1,j,k-1)

              deltaR=dh*r(i,j,k-1)-dP
              deltaL=dP-dh*r(i-1,j,k-1)
              if (deltaR*deltaL < Zero) then
                rr=Zero
              elseif (abs(deltaR) > Two*abs(deltaL)) then
                rr=Three*deltaL
              elseif (abs(deltaL) > Two*abs(deltaR)) then
                rr=Three*deltaR
              else
                rr=deltaR+deltaL
              endif

              FC(i,k-1)=Half*dh*( P(i,j,k-1)+P(i-1,j,k-1)
     &                                       +OneSixth*rr)

              ru(i,j,k)=Two*( FX(i-1,j,k)-FX(i,j,k)+FC(i,k)
     &                 -FC(i,k-1) )/(Hz(i-1,j,k)+Hz(i,j,k))
# ifdef MASKING
     &                                         *umask(i,j)
# endif
            enddo
          enddo        !--> discard FC, keep ru,FX
        endif

        if (j>=jstrV-1) then
          do i=istr,iend
            FC(i,N)=Zero    !--> set pressure at free surface.
          enddo
          do k=N,1,-1
            do i=istr,iend
              dP=P(i,j-1,k-1)-P(i,j,k-1)
              dh=z_w(i,j,k-1)-z_w(i,j-1,k-1)

              deltaR=dh*r(i,j,k-1)-dP
              deltaL=dP-dh*r(i,j-1,k-1)
              if (deltaR*deltaL < Zero) then
                rr=Zero
              elseif (abs(deltaR) > Two*abs(deltaL)) then
                rr=Three*deltaL
              elseif (abs(deltaL) > Two*abs(deltaR)) then
                rr=Three*deltaR
              else
                rr=deltaR+deltaL
              endif

              FC(i,k-1)=Half*dh*( P(i,j,k-1)+P(i,j-1,k-1)
     &                                       +OneSixth*rr)

              rv(i,j,k)=Two*( FX(i,j-1,k)-FX(i,j,k)+FC(i,k)
     &                 -FC(i,k-1) )/(Hz(i,j-1,k)+Hz(i,j,k))
# ifdef MASKING
     &                                         *vmask(i,j)
# endif
            enddo
          enddo
        endif
      enddo    !<-- j        !--> discard FC,FX,  keep ru,rv

! After this moment arrays ru,rv contain baroclinic part of
! horizontal (along geopotential) pressure gradient components,
! computed with higher-order accuracy with respect to vertical
! derivatives of density, but only with second-order with respect
! to horizontal (along sigma) derivatives. The second-order
! truncation term comes from the fact that in the code above,
! in computation of upper/lower facet fluxes "FC"s partial
! derivatives of pressure along sigma levels were replaced with
! vertical derivatives of pressure (i.e., -density) projected by
! the coordinate slope. Since,
!
!           dP |                 dP |               dP     dz
!          ----|     =          ----|          +   ---- * ----
!           dx |                 dx |               dz     dx
!              | along sigma        | along z
!
!                             |<-  missing  ->|  |<- projected ->|
!                                in the code         vertical
!                                   above            derivative

! The following code segments compensate the second-order truncation
! error by adding in the missing part into "FCs" above.

      rr=g/(24.*rho0)
      cff=0.5*g
      cff1=0.5*g/rho0
      do j=jstr,jend
        do k=N-1,1,-1
          do i=istrU,iend
            dh=rr*(z_w(i,j,k)-z_w(i-1,j,k))
            FC(i,  k)=max(dh,Zero)*( ru(i,j,k+1) +ru(i+1,j,k)
     &                              -ru(i,j,k) -ru(i-1,j,k+1))
     &               +min(dh,Zero)*( ru(i+1,j,k+1) +ru(i,j,k)
     &                              -ru(i,j,k+1) -ru(i-1,j,k))
          enddo
        enddo
        do i=istrU,iend
          FC(i,N)=0.
          dh=rr*(z_w(i,j,0)-z_w(i-1,j,0))
          FC(i,0)=max(dh,Zero)*(ru(i,j,1) -ru(i-1,j,1))
     &           +min(dh,Zero)*(ru(i+1,j,1) -ru(i,j,1))
        enddo
        do k=1,N
          do i=istrU,iend
            ru(i,j,k)=( (Hz(i-1,j,k)+Hz(i,j,k))*(
     &         cff*(z_w(i-1,j,N)-z_w(i,j,N)) + cff1*ru(i,j,k)
     &                      ) + FC(i,k)-FC(i,k-1) )*dn_u(i,j)
          enddo
        enddo
      enddo

      do j=jstrV,jend
        do k=N-1,1,-1
          do i=istr,iend
            dh=rr*(z_w(i,j,k)-z_w(i,j-1,k))
            FX(i,j,k)=max(dh,Zero)*( rv(i,j,k+1) +rv(i+1,j,k)
     &                              -rv(i,j,k) -rv(i,j-1,k+1))
     &               +min(dh,Zero)*( rv(i+1,j,k+1) +rv(i,j,k)
     &                              -rv(i,j,k+1) -rv(i,j-1,k))
          enddo
        enddo
        do i=istr,iend
          FX(i,j,N)=0.
          dh=rr*(z_w(i,j,0)-z_w(i,j-1,0))
          FX(i,j,0)=max(dh,Zero)*(rv(i,j,1) -rv(i,j-1,1))
     &             +min(dh,Zero)*(rv(i+1,j,1) -rv(i,j,1))
        enddo
      enddo
      do j=jstrV,jend
        do k=1,N
          do i=istr,iend
            rv(i,j,k)=( (Hz(i,j-1,k)+Hz(i,j,k))*(
     &         cff*(z_w(i,j-1,N)-z_w(i,j,N)) + cff1*rv(i,j,k)
     &                  ) + FX(i,j,k)-FX(i,j,k-1) )*dm_v(i,j)
          enddo
        enddo
      enddo
      end

#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */
