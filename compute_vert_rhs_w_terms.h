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
!  Extrapolate for advw(N) using 2*advw(N-1)- advw(N-2)

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

          do k=1,N-1
            do i=istr,iend
              rw(i,j,k) = rw(i,j,k) - FC(i,k+1) + FC(i,k)
            enddo
          enddo

	  ! Extrapolate adv(w): adv(N) = 2*adv(N-1) -adv(N-2)
	  !                          = -2*FC(N) + 3*FC(N-1) -FC(N-2)
	  ! Multiply by 0.5 because of half volume at surface
          do i=istr,iend
            rw(i,j,N) = rw(i,j,N)-FC(i,N)+1.5*FC(i,N-1)-0.5*FC(i,N-2)
          enddo
#endif
