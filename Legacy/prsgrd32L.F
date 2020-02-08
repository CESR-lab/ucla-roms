#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine prsgrd (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call prsgrd32_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                                A3d(1,3), A2d(1,1), A2d(1,2),
     &                      A2d(1,1), A2d(1,2), A2d(1,3), A2d(1,4))
      end

      subroutine prsgrd32_tile (istr,iend,jstr,jend, ru,rv, FX,
     &                                   dR,dZ, aux,dRx,dZx,FC)
      implicit none
      integer istr,iend,jstr,jend, i,j,k, k1,k2
# include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv, FX
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: dR,dZ
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: aux,dRx,dZx
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: FC

      real GRho,HalfGRho
      real, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=1.E-33

# include "grid.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "scalars.h"

! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensional coordinates xi,
! eta, and s (basically their respective fortran indices).

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
c            cff=2.*dR(i,k)*dR(i,k-1)
c            if (cff>epsil) then
c              dR(i,k)=cff/(dR(i,k)+dR(i,k-1))
c            else
c              dR(i,k)=0.
c            endif
c            dZ(i,k)=2.*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))

             dR(i,k)=0.5*(dR(i,k)+dR(i,k-1))
             dZ(i,k)=0.5*(dZ(i,k)+dZ(i,k-1))

          enddo
        enddo
        do k=1,N-1
          do i=istrU-1,iend
            FX(i,j,k)=HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                               *(z_r(i,j,k+1)-z_r(i,j,k))

     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                             ))
          enddo
        enddo
        do i=istrU-1,iend
          FX(i,j,N)=GRho*( rho(i,j,N) +0.5*(rho(i,j,N)-rho(i,j,N-1))
     &             *(z_w(i,j,N)-z_r(i,j,N))/(z_r(i,j,N)-z_r(i,j,N-1))
     &                                      )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
      enddo

! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----

      do k=N,1,-1     !--> irreversible
        k1=1+mod(k,2)
        k2=3-k1
        do j=jstr,jend
          do i=istrU-1,iend+1
            aux(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
#ifdef MASKING
     &                               *umask(i,j)
#endif
            FC(i,j,k1)=(rho(i,j,k)-rho(i-1,j,k))
#ifdef MASKING
     &                               *umask(i,j)
#endif
          enddo
        enddo

        do j=jstr,jend
          do i=istrU-1,iend
c            cff=2.*aux(i,j)*aux(i+1,j)
c            if (cff>epsil) then
c              dZx(i,j)=cff/(aux(i,j)+aux(i+1,j))
c            else
c              dZx(i,j)=0.
c            endif
c            cff1=2.*FC(i,j,k1)*FC(i+1,j,k1)
c            if (cff1>epsil) then
c              dRx(i,j)=cff1/(FC(i,j,k1)+FC(i+1,j,k1))
c            else
c              dRx(i,j)=0.
c            endif

             dZx(i,j)=0.5*(aux(i,j)+aux(i+1,j))
             dRx(i,j)=0.5*(FC(i,j,k1)+FC(i+1,j,k1))

          enddo
        enddo                 !--> discard aux, FC(:,:,k1)

        do j=jstr,jend
          do i=istrU,iend
            FC(i,j,k1)=HalfGRho*( (rho(i,j,k)+rho(i-1,j,k))
     &                                 *(z_r(i,j,k)-z_r(i-1,j,k))

     & -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                          -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     &            -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                          -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                           ))
          enddo
        enddo
        if (k==N) then
          do j=jstr,jend
            do i=istrU,iend
              ru(i,j,N)=g*(z_w(i-1,j,N)-z_w(i,j,N))
     &                 +FX(i-1,j,N)-FX(i,j,N)-FC(i,j,k1)
            enddo
          enddo
        else
          do j=jstr,jend
            do i=istrU,iend
              ru(i,j,k)=ru(i,j,k+1) +FX(i-1,j,k)-FX(i,j,k)
     &                               +FC(i,j,k2)-FC(i,j,k1)

              ru(i,j,k+1)=ru(i,j,k+1)*0.5*( Hz(i,j,k+1)+
     &                            Hz(i-1,j,k+1) )*dn_u(i,j)
            enddo
          enddo
        endif
      enddo
      do j=jstr,jend
        do i=istrU,iend
          ru(i,j,1)=ru(i,j,1)*0.5*(Hz(i,j,1)+Hz(i-1,j,1))*dn_u(i,j)
        enddo
      enddo

! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----

      do k=N,1,-1     !--> irreversible
        k1=1+mod(k,2)
        k2=3-k1
        do j=jstrV-1,jend+1
          do i=istr,iend
            aux(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
#ifdef MASKING
     &                               *vmask(i,j)
#endif
            FC(i,j,k1)=(rho(i,j,k)-rho(i,j-1,k))
#ifdef MASKING
     &                               *vmask(i,j)
#endif
          enddo
        enddo

        do j=jstrV-1,jend
          do i=istr,iend
c            cff=2.*aux(i,j)*aux(i,j+1)
c            if (cff>epsil) then
c              dZx(i,j)=cff/(aux(i,j)+aux(i,j+1))
c            else
c              dZx(i,j)=0.
c            endif
c            cff1=2.*FC(i,j,k1)*FC(i,j+1,k1)
c            if (cff1>epsil) then
c              dRx(i,j)=cff1/(FC(i,j,k1)+FC(i,j+1,k1))
c            else
c              dRx(i,j)=0.
c            endif

            dZx(i,j)=0.5*(aux(i,j)+aux(i,j+1))
            dRx(i,j)=0.5*(FC(i,j,k1)+FC(i,j+1,k1))

          enddo         !--> discard aux, FC(:,:,k1)
        enddo

        do j=jstrV,jend
          do i=istr,iend
            FC(i,j,k1)=HalfGRho*( (rho(i,j,k)+rho(i,j-1,k))
     &                                 *(z_r(i,j,k)-z_r(i,j-1,k))

     &  -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                          -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     &             -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                          -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                           ))
          enddo         !-->  discard dRx, dZx
        enddo
        if (k==N) then
          do j=jstrV,jend
            do i=istr,iend
              rv(i,j,N)=g*(z_w(i,j-1,N)-z_w(i,j,N))
     &                 +FX(i,j-1,N)-FX(i,j,N)-FC(i,j,k1)
            enddo
          enddo
        else
          do j=jstrV,jend
            do i=istr,iend
              rv(i,j,k)=rv(i,j,k+1) +FX(i,j-1,k)-FX(i,j,k)
     &                               +FC(i,j,k2)-FC(i,j,k1)

              rv(i,j,k+1)=rv(i,j,k+1)*0.5*( Hz(i,j,k+1)+
     &                            Hz(i,j-1,k+1) )*dm_v(i,j)
            enddo
          enddo
        endif
      enddo
      do j=jstrV,jend
        do i=istr,iend
          rv(i,j,1)=rv(i,j,1)*0.5*(Hz(i,j,1)+Hz(i,j-1,1))*dm_v(i,j)
        enddo
      enddo
      end
#else
      subroutine prsgrd_empty
      end
#endif /* SOLVE3D */
