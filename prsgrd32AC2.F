#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd32A_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                        A3d(1,3), A2d(1,1), A2d(1,2), A2d(1,3),
     &                                  A2d(1,1), A2d(1,2), A2d(1,3),
     &                                 A2d(1,4), A2d(1,5), A2d(1,6))
      end

      subroutine prsgrd32A_tile (istr,iend,jstr,jend, ru,rv, P,
     &                         dR,dZ,dq, FC,dZx,rx, dRx,qx,dqx)
      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv, P
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: dR, dZ, dq
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FC,dZx,rx,dRx,qx,dqx

      real GRho,  HalfGRho, cff,cfr,cfq
      real, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=1.E-33

# include "grid.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "scalars.h"

! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensional coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid. These
! derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.

! If CPP-switch SPLIT_EOS is defined, the Equation of State (EOS) is
! assumed to have form
!
!           rho(T,S,z) = rho1(T,S) + (zeta-z) * qp1(T,S)
!
! where rho1 is potential density at 1 atm and qp1 is compressibility
! coefficient, which does not depend on z (or weakly dependens on z).
! Matematical algorithm in the code below is based on identity
!
!                                          1
!          J(rho1+qp*z, z) = J(rho1, z) + --- J(qp1, z^2)
!                                          2
! where both Jacobians in the r.h.s. are then discretized in pseudo-
! flux form using cubic polynomial fits for z, rho1 and qp1. Doing so
! requires approximation for spatial derivatives or rho1 and qp1 at
! the same locations as rho1 and qp1 themselves. They are computed as
! harmonic averages of their respective elementary differences
! independently one from the other.

# include "compute_auxiliary_bounds.h"

! Preliminary step (same for XI- and ETA-components:
!------------ ---- ----- --- --- --- ---------------

      GRho=g/rho0
      HalfGRho=0.5*GRho

      do j=jstrV-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)
# ifdef SPLIT_EOS
               dq(i,k)=qp1(i,j,k+1)-qp1(i,j,k)
# endif
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
# ifdef SPLIT_EOS
          dq(i,N)=dq(i,N-1)
          dq(i,0)=dq(i,1)
# endif
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1               !--> irreversible
          do i=istrU-1,iend
            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr>epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0.
            endif
# ifdef SPLIT_EOS
            cfq=2.*dq(i,k)*dq(i,k-1)
            if (cfq>epsil) then
              dq(i,k)=cfq/(dq(i,k)+dq(i,k-1))
            else
              dq(i,k)=0.
            endif
# endif
            dZ(i,k)=2.*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))
          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho1(i,j,N)
     &       +0.5*(rho1(i,j,N)-rho1(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &           /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (z_r(i,j,k+1)-z_r(i,j,k))
     &                                   *( rho1(i,j,k+1)+rho1(i,j,k)
# ifdef SPLIT_EOS
     &                  +( z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k)) )
     &                                     *(qp1(i,j,k+1)+qp1(i,j,k))
# endif
     &                                                              )

     &   -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     &            -(dZ(i,k+1)-dZ(i,k))*( rho1(i,j,k+1)-rho1(i,j,k)
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
# ifdef SPLIT_EOS
     &         -(dq(i,k+1)-dq(i,k))*( 0.5*(z_r(i,j,k+1)+z_r(i,j,k))
     &                                   *(z_r(i,j,k+1)-z_r(i,j,k))
     &      -OneTwelfth*(z_r(i,j,k+1)*dZ(i,k+1)+z_r(i,j,k)*dZ(i,k)) )

     &                 +(z_r(i,j,k+1)*dZ(i,k+1)-z_r(i,j,k)*dZ(i,k))
     &     *(qp1(i,j,k+1)-qp1(i,j,k)-OneTwelfth*(dq(i,k+1)+dq(i,k)))
# endif
     &                                                             ))
          enddo
        enddo
      enddo

! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----

      do k=N,1,-1
        do j=jstr,jend
          do i=istrU-1,iend+1
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
# ifdef MASKING
     &                              *umask(i,j)
# endif
            rx(i,j)=(rho1(i,j,k)-rho1(i-1,j,k))
# ifdef MASKING
     &                              *umask(i,j)
# endif
# ifdef SPLIT_EOS
            qx(i,j)=(qp1(i,j,k)-qp1(i-1,j,k))
#  ifdef MASKING
     &                              *umask(i,j)
#  endif
# endif
          enddo
        enddo

        do j=jstr,jend
          do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff>epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif
            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr>epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            cfq=2.*qx(i,j)*qx(i+1,j)
            if (cfq>epsil) then
              dqx(i,j)=cfq/(qx(i,j)+qx(i+1,j))
            else
              dqx(i,j)=0.
            endif
# endif
          enddo
        enddo                 !--> discard FC, rx, qx

        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &      P(i-1,j,k)-P(i,j,k)-HalfGRho*( (z_r(i,j,k)-z_r(i-1,j,k))
     &                                   *( rho1(i,j,k)+rho1(i-1,j,k)
# ifdef SPLIT_EOS
     &                                 +0.5*(qp1(i,j,k)+qp1(i-1,j,k))
     &                                    *( z_w(i,j,N)+z_w(i-1,j,N)
     &                                      -z_r(i,j,k)-z_r(i-1,j,k))
# endif
     &                                                              )

     &  -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     &           -(dZx(i,j)-dZx(i-1,j))*( rho1(i,j,k)-rho1(i-1,j,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
# ifdef SPLIT_EOS
     &       -(dqx(i,j)-dqx(i-1,j))*( 0.5*(z_r(i,j,k)+z_r(i-1,j,k))
     &                                   *(z_r(i,j,k)-z_r(i-1,j,k))
     &    -OneTwelfth*(z_r(i,j,k)*dZx(i,j)+z_r(i-1,j,k)*dZx(i-1,j)) )

     &               +(z_r(i,j,k)*dZx(i,j)-z_r(i-1,j,k)*dZx(i-1,j))
     &   *(qp1(i,j,k)-qp1(i-1,j,k)-OneTwelfth*(dqx(i,j)+dqx(i-1,j)))
# endif
     &                                                            )))
          enddo
        enddo       !--> discard dZx, dRx, dqx
!
! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----
!
        do j=jstrV-1,jend+1
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
# ifdef MASKING
     &                              *vmask(i,j)
# endif
            rx(i,j)=(rho1(i,j,k)-rho1(i,j-1,k))
# ifdef MASKING
     &                              *vmask(i,j)
# endif
# ifdef SPLIT_EOS
            qx(i,j)=(qp1(i,j,k)-qp1(i,j-1,k))
#  ifdef MASKING
     &                              *vmask(i,j)
#  endif
# endif
          enddo
        enddo

        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff>epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif
            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr>epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            cfq=2.*qx(i,j)*qx(i,j+1)
            if (cfq>epsil) then
              dqx(i,j)=cfq/(qx(i,j)+qx(i,j+1))
            else
              dqx(i,j)=0.
            endif
# endif
          enddo
        enddo                 !--> discard FC, rx, qx

        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &      P(i,j-1,k)-P(i,j,k)-HalfGRho*( (z_r(i,j,k)-z_r(i,j-1,k))
     &                                   *( rho1(i,j,k)+rho1(i,j-1,k)
# ifdef SPLIT_EOS
     &                                 +0.5*(qp1(i,j,k)+qp1(i,j-1,k))
     &                                    *( z_w(i,j,N)+z_w(i,j-1,N)
     &                                      -z_r(i,j,k)-z_r(i,j-1,k))
# endif
     &                                                              )

     &  -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     &           -(dZx(i,j)-dZx(i,j-1))*( rho1(i,j,k)-rho1(i,j-1,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
# ifdef SPLIT_EOS
     &       -(dqx(i,j)-dqx(i,j-1))*( 0.5*(z_r(i,j,k)+z_r(i,j-1,k))
     &                                   *(z_r(i,j,k)-z_r(i,j-1,k))
     &    -OneTwelfth*(z_r(i,j,k)*dZx(i,j)+z_r(i,j-1,k)*dZx(i,j-1)) )

     &               +(z_r(i,j,k)*dZx(i,j)-z_r(i,j-1,k)*dZx(i,j-1))
     &   *(qp1(i,j,k)-qp1(i,j-1,k)-OneTwelfth*(dqx(i,j)+dqx(i,j-1)))
# endif
     &                                                            )))
          enddo
        enddo                    !--> discard dZx, dRx, dqx
      enddo     !<-- k
      end
#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */
