#include "cppdefs.opt"
#define SEDIMENT

#ifdef SEDIMENT
      subroutine sediment (tile)

      use param
      use grid
      use hidden_mpi_vars
      use private_scratch

      implicit none
      integer tile

# include "compute_tile_bounds.h"
      call sediment_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                     A2d(1,3), A2d(1,4), A2d(1,5), A2d(1, 6),
     &                     A2d(1,7), A2d(1,8), A2d(1,9), A2d(1,10))
      return
      end

      subroutine sediment_tile (istr,iend,jstr,jend, FC, qc, qR,qL,
     &                      WR,WL, Hz_inv,Hz_inv2,Hz_inv3, ksource)

      use param
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend,   ised,indx, i,j,k,ks
      real FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),   cff,
     &       qc(PRIVATE_1D_SCRATCH_ARRAY,N),   cff1, cu,
     &       qR(PRIVATE_1D_SCRATCH_ARRAY,N),   cff2, cff3,
     &       qL(PRIVATE_1D_SCRATCH_ARRAY,N),
     &       WR(PRIVATE_1D_SCRATCH_ARRAY,N),   cffR, dltR,
     &       WL(PRIVATE_1D_SCRATCH_ARRAY,N),   cffL, dltL,
     &   Hz_inv(PRIVATE_1D_SCRATCH_ARRAY,N),
     &  Hz_inv2(PRIVATE_1D_SCRATCH_ARRAY,N),   tau_water,
     &  Hz_inv3(PRIVATE_1D_SCRATCH_ARRAY,N),   ero_flux
      integer ksource(PRIVATE_1D_SCRATCH_ARRAY,N)

# include "forces.h"
c# include "bbl.h"
c# include "sediment.h"

c??
c??
c??
      integer NST, ithck
      parameter (NST=1, ithck=1)
      integer idsed(NST)
      real bed(GLOBAL_2D_ARRAY,N,ithck)
      real Wsed(NST),Erate(NST),Srho(NST),poros(NST),tau_ce(NST)
c??
c??

      do j=jstr,jend                        ! Auxiliary step: save
        do k=1,N                            ! inverses of grid box
          do i=istr,iend                    ! heights to avoid
            Hz_inv(i,k)=1./Hz(i,j,k)        ! repeated divisions
          enddo
        enddo
        do k=1,N-1
          do i=istr,iend
            Hz_inv2(i,k)=1./(Hz(i,j,k)+Hz(i,j,k+1))
          enddo
        enddo
        do k=2,N-1
          do i=istr,iend
            Hz_inv3(i,k)=1./(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          enddo
        enddo
!
! Vertical sinking of suspended particles:   Copy concentration
!========= ======= == ========= ==========   of suspended sediment
! into scratch array "qc" (q-central) which are hereafter interpreted
! as a set of grid-box averaged values for concentration.
! Then reconstruct vertical profile of "qc" in terms of a set of
! parabolic segments within each grid box; and, finally, compute
! semi-Lagrangian flux due to sinking.
!
        do ised=1,NST
          indx=idsed(ised)
          do k=1,N
            do i=istr,iend
              qc(i,k)=t(i,j,k,nnew,indx)
            enddo
          enddo

          do k=N-1,1,-1
            do i=istr,iend
              FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
            enddo
          enddo
          do k=2,N-1
            do i=istr,iend
              dltR=Hz(i,j,k)*FC(i,k)
              dltL=Hz(i,j,k)*FC(i,k-1)
              cff=Hz(i,j,k-1)+2.*Hz(i,j,k)+Hz(i,j,k+1)
              cffR=cff*FC(i,k)
              cffL=cff*FC(i,k-1)            ! Apply PPM monotonicity
                                            ! constraint to prevent
              if (dltR*dltL <= 0.) then   ! oscillation within the
                dltR=0.                     ! grid box
                dltL=0.
              elseif (abs(dltR) > abs(cffL)) then
                dltR=cffL
              elseif (abs(dltL) > abs(cffR)) then
                dltL=cffR
              endif                         ! Compute right and left
                                            ! side values qR,qL of
              cff=(dltR-dltL)*Hz_inv3(i,k)  ! parabolic segments
              dltR=dltR-cff*Hz(i,j,k+1)     ! within grid box Hz(k)
              dltL=dltL+cff*Hz(i,j,k-1)     ! (WR,WL are measures of
              qR(i,k)=qc(i,k)+dltR          ! quadratic variations).
              qL(i,k)=qc(i,k)-dltL
              WR(i,k)=( 2.*dltR-dltL )**2   ! NOTE: Although each
              WL(i,k)=( dltR-2.*dltL )**2   ! parabolic segment is
            enddo                           ! monotone within its
          enddo          !--> discard FC    ! grid box, monotonicity
                                            ! of the whole profile is
          cff=1.0E-14                       ! not guaranteed, because
          do k=2,N-2                        ! qL(k+1)-qR(k) may still
            do i=istr,iend                  ! have different sign
              dltL=max(WL(i,k),   cff)      ! than qc(k+1)-qc(k)...
              dltR=max(WR(i,k+1), cff)
              qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
              qL(i,k+1)=qR(i,k)
            enddo                           !   ...this possibility
          enddo      !--> discard WR,WL     ! is excluded, after qL
                                            ! and qR are reconciled
          do i=istr,iend                    ! using WENO procedire.
            FC(i,N)=0.     !<-- no-flux BC
# if defined LINEAR_CONTINUATION
            qL(i,N)=qR(i,N-1)
            qR(i,N)=2.*qc(i,N)-qL(i,N)
# elif defined NEUMANN
            qL(i,N)=qR(i,N-1)
            qR(i,N)=1.5*qc(i,N)-0.5*qL(i,N)
# else
            qR(i,N)=qc(i,N)                 ! Strictly monotone
            qL(i,N)=qc(i,N)                 ! version as the default:
            qR(i,N-1)=qc(i,N)               ! distributions at top...
# endif
# if defined LINEAR_CONTINUATION
            qR(i,1)=qL(i,2)
            qL(i,1)=2.*qc(i,1)-qR(i,1)
# elif defined NEUMANN
            qR(i,1)=qL(i,2)
            qL(i,1)=1.5*qc(i,1)-0.5*qR(i,1)
# else
            qL(i,2)=qc(i,1)                 !     ...and bottom grid
            qR(i,1)=qc(i,1)                 ! boxes re assumed to be
            qL(i,1)=qc(i,1)                 ! piecewise constant.
# endif
          enddo

          do k=1,N                          ! Since the reconciled
            do i=istr,iend                  ! interfacial values may
              dltR=qR(i,k)-qc(i,k)          ! cause non-monotonic
              dltL=qc(i,k)-qL(i,k)          ! behavior of parabolic
              cffR=2.*dltR                  ! segments inside grid
              cffL=2.*dltL                  ! box apply monotonicity
                                            ! constraint again.
              if (dltR*dltL < 0.) then
                dltR=0.
                dltL=0.
              elseif (abs(dltR) > abs(cffL)) then
                dltR=cffL
              elseif (abs(dltL) > abs(cffR)) then
                dltL=cffR
              endif

              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
            enddo
          enddo     !--> discard everything, except qR,qL
!
! After this moment reconstruction is considered complete. The next
! stage is to compute vertical advective fluxes FC. It is expected
! that sinking may occurs relatively fast, the algorithm is designed
! to be free of CFL criterion, which is achieved by allowing
! integration bounds for semi-Lagrangian advective flux to use as
! many grid boxes in upstream direction as necessary.
!
          cff=dt*abs(Wsed(ised))          ! In the two code segments
          do k=1,N                        ! WL is z-coordinate of the
            do i=istr,iend                ! departure point for grid
              FC(i,k-1)=0.                ! box interface z_w with
              WL(i,k)=z_w(i,j,k-1)+cff    ! the same indices;
              WR(i,k)=Hz(i,j,k)*qc(i,k)   ! FC is finite volume flux;
              ksource(i,k)=k              ! ksource(:,k) is index of
            enddo                         ! vertical grid box which
          enddo                           ! contains the departure
          do k=1,N                        ! point (restricted by N);
            do ks=k,N-1                   ! During the search: also
              do i=istr,iend
                if (WL(i,k) > z_w(i,j,ks)) then
                  ksource(i,k)=ks+1
                  FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                endif
              enddo                       ! add in content of whole
            enddo                         ! grid boxes participating
          enddo       !--> discard WR     ! in FC.

          do k=1,N                        ! Finalize computation of
            do i=istr,iend                ! flux: add fractional part
              ks=ksource(i,k)
              cu=min(1.,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
              FC(i,k-1)=FC(i,k-1) + Hz(i,j,ks)*cu*( qL(i,ks)
     &                             +cu*( 0.5*(qR(i,ks)-qL(i,ks))
     &             -(1.5-cu)*(qR(i,ks)+qL(i,ks)-2.*qc(i,ks))  ))
            enddo
          enddo              !--> discard WL
          do k=1,N,+1
            do i=istr,iend
              qc(i,k)=qc(i,k) + (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
            enddo
          enddo    !--> discard everything, except qc,FC(:,0)
!
!
! Deposition and resuspension near the bottom:  Update thickness of
!============ === =========== ==== === =======  top sediment layer on
! the bottom "bed" (meters), which changes due to precepitation flux
! FC(:,0) which is already computed and resuspension (erosion, hence
! called ero_flux), which is also applied to the bottom-most grid box
! value qc(:,1) (hence total mass is concerved; also restrict
! "ero_flux" so that "bed" cannot go negative after both fluxes are
! applied).
!
          cff=dt*Erate(ised)*(1.-poros(ised))
          cff1=1./tau_ce(ised)
          cff2=Srho(ised)*(1.-poros(ised))
          cff3=1./cff2
          do i=istr,iend
# ifdef BBL
            tau_water=0.5*sqrt( (bustrw(i,j)+bustrw(i+1,j))**2
     &                         +(bvstrw(i,j)+bvstrw(i,j+1))**2)
# else
            tau_water=0.5*sqrt( (bustr(i,j)+bustr(i+1,j))**2
     &                         +(bvstr(i,j)+bvstr(i,j+1))**2)
# endif
            ero_flux=min( max(0., cff*(cff1*tau_water-1.)),
     &                  cff2*bed(i,j,1,ithck)+0.001*FC(i,0) )

            bed(i,j,1,ithck)=( bed(i,j,1,ithck) -cff3*( ero_flux
     &                                          -0.001*FC(i,0) ))
# ifdef MASKING
     &                                                *rmask(i,j)
# endif
            qc(i,1)=qc(i,1) + 1000.*ero_flux*Hz_inv(i,1)
          enddo

          do k=1,N                         ! Write concentration of
            do i=istr,iend                 ! suspended sediment back
                t(i,j,k,nnew,indx)=qc(i,k) ! to its global array.
# ifdef MASKING
     &                         *rmask(i,j)
# endif
            enddo
          enddo
        enddo    ! <-- ised
      enddo    ! <-- j
      return
      end
#else
      subroutine sediment_empty
      end
#endif

