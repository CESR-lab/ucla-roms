#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd32A_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                                A3d(1,3),  A2d(1,1), A2d(1,2),
     &                       A2d(1,1), A2d(1,2), A2d(1,3), A2d(1,4))
      end

      subroutine prsgrd32A_tile (istr,iend,jstr,jend, ru,rv, P,
     &                                   dR,dZ, FC,aux,dRx,dZx)
      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv, P
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: dR,dZ
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: aux,FC,dRx,dZx

      real GRho, HalfGRho, cff,cff1
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

! NOTE: This routine is an alternative form of prsgrd32 and it
!       produces results identical to that if its prototype.

# include "compute_auxiliary_bounds.h"

! Preliminary step (same for XI- and ETA-components:
!------------ ---- ----- --- --- --- ---------------
      GRho=g/rho0
      HalfGRho=0.5*GRho

      do j=jstrV-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dZ(i,N)=dZ(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1               !--> irreversible
          do i=istrU-1,iend
            cff=2.*dR(i,k)*dR(i,k-1)
            if (cff>epsil) then
              dR(i,k)=cff/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0.
            endif
            dZ(i,k)=2.*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))
          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     &        +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &           /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                                     *(z_r(i,j,k+1)-z_r(i,j,k))

     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                             ))
          enddo
        enddo
      enddo   !<-- j

! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----

      do k=N,1,-1
        do j=jstr,jend
          do i=istrU-1,iend+1
            aux(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
#ifdef MASKING
     &                              *umask(i,j)
#endif
            FC(i,j)=(rho(i,j,k)-rho(i-1,j,k))
#ifdef MASKING
     &                              *umask(i,j)
#endif
          enddo
        enddo

        do j=jstr,jend
          do i=istrU-1,iend
            cff=2.*aux(i,j)*aux(i+1,j)
            if (cff>epsil) then
              dZx(i,j)=cff/(aux(i,j)+aux(i+1,j))
            else
              dZx(i,j)=0.
            endif
            cff1=2.*FC(i,j)*FC(i+1,j)
            if (cff1>epsil) then
              dRx(i,j)=cff1/(FC(i,j)+FC(i+1,j))
            else
              dRx(i,j)=0.
            endif
          enddo
        enddo                 !--> discard aux, FC

        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &                           P(i-1,j,k)-P(i,j,k)-HalfGRho*(
     &          (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))

     &   -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     &              -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                            )))
          enddo
        enddo

! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----

        do j=jstrV-1,jend+1
          do i=istr,iend
            aux(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
#ifdef MASKING
     &                              *vmask(i,j)
#endif
            FC(i,j)=(rho(i,j,k)-rho(i,j-1,k))
#ifdef MASKING
     &                              *vmask(i,j)
#endif
          enddo
        enddo                 !--> discard aux, FC

        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*aux(i,j)*aux(i,j+1)
            if (cff>epsil) then
              dZx(i,j)=cff/(aux(i,j)+aux(i,j+1))
            else
              dZx(i,j)=0.
            endif
            cff1=2.*FC(i,j)*FC(i,j+1)
            if (cff1>epsil) then
              dRx(i,j)=cff1/(FC(i,j)+FC(i,j+1))
            else
              dRx(i,j)=0.
            endif
          enddo
        enddo

        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &      P(i,j-1,k)-P(i,j,k)-HalfGRho*( (rho(i,j,k)+rho(i,j-1,k))
     &                                     *(z_r(i,j,k)-z_r(i,j-1,k))

     &   -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     &              -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                            )))
          enddo
        enddo
      enddo
      end
#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */
