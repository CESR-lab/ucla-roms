#ifdef UV_ADV

! This module "compute_vert_rhs_w_terms.h" computes vertical
! advective fluxes for the w equation and adds to rw.
! This should be included in a giant j-loop. Be careful about re-using
! work arrays.
! Used work arrays:  curv, Fz. Flxw is a scalar

! The only thing I've done is UPSTREAM_W (so no cpp flags :) )
!  w_t =  Fz(k+1/2) - Fz(k-1/2)
!  Fz(k+1/2) = w(k+1/2)*0.5(We(i,j,k+1)+We(i,j,k))
!  w(k+1/2) =  0.5*(w(i,j,k+1)+w(i,j,k)) - 0.166666*curv
!  at k+1/2: add curv(k..k+2) if w<0 otherw add curv(k-1...k+1)
!  curv = w(k+1) - 2*w(k) + w(k-1)

# define curv CF

          do k=1,N-1
            do i=istr,iend
              curv(i,k) = w(i,j,k+1,nrhs) - 2*w(i,j,k,nrhs) + w(i,j,k-1,nrhs)
            enddo
          enddo
          do i=istr,iend !! extrapolate to k=N
            curv(i,N) = 2*curv(i,N-1) - curv(i,N-2)
          enddo
          do k=1,N
            do i=istr,iend
              Flxw = 0.5*(We(i,j,k)+We(i,j,k-1))
              FC(i,k) = 0.5*(w(i,j,k,nrhs)+w(i,j,k-1,nrhs))*Flxw
     &          -0.1666666666666666*( curv(i,k-1)*max(Flxw,0.)
     &                               +curv(i,k  )*min(Flxw,0.))
            enddo
          enddo

          !! These guys (FC) are fluxes (m3/s) of velocity (m/s) -> m4/s2
          do k=1,N-1
            do i=istr,iend
              rw(i,j,k) = rw(i,j,k) - FC(i,k+1) + FC(i,k)
            enddo
          enddo

          ! flux at surface is zero because it moves with omega.
          ! this also takes care of the half-volume of w(N)
          do i=istr,iend
            rw(i,j,N) = rw(i,j,N)   + FC(i,N)
          enddo
#endif
