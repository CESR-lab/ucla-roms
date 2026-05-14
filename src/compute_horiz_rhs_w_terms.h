! Compute horizontal fluxes for vertical flux w. Essentially interpolate
! tracer values from their native locations on C grid to horizontal
! velocity points with simultaneous translation from grid-box-averages
! to instantaneous values at interface location.  2 options an be
! selected: 3-point upstream-biased parabolic interpolation (UPSTREAM_W);
!  4-point symmetric fourth-order method (undefined state of both CPP
! switches); 

! This code is extracted into a special module because it is used
! twice, in predictor and corrector substeps for w

#ifdef UV_ADV
# define UPSTREAM_W

#ifdef UPSTREAM_W
# define curv wrk1
#else 
# define grad wrk1
#endif
          do j=1,ny
            do i=0,nx+2
              FX(i,j)=(w(i,j,k,nrhs)-w(i-1,j,k,nrhs))*umask(i,j)
            enddo
          enddo

          if (west_bnd) then
            do j=1,ny
              FX(0,j)=FX(1,j)
            enddo
          endif
          if (east_bnd) then
            do j=1,ny
              FX(nx+2,j)=FX(nx+1,j)
            enddo
          endif

          do j=1,ny
            do i=0,nx+1
#if defined UPSTREAM_W
              curv(i,j)=FX(i+1,j)-FX(i,j)
#else
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
#endif
            enddo
          enddo             !--> discard FX
          do j=1,ny
            do i=1,nx+1
              if (k<nz) then
                Uflxw= 0.5*(FlxU(i,j,k)+FlxU(i,j,k+1))
              else
                Uflxw= 0.5*(FlxU(i,j,k)              )
              endif
#ifdef UPSTREAM_W
              FX(i,j)=0.5*(w(i,j,k,nrhs)+w(i-1,j,k,nrhs))
     &                                       *Uflxw
     &          -0.1666666666666666*( curv(i-1,j)*max(Uflxw,0.)
     &                               +curv(i  ,j)*min(Uflxw,0.))
#else
              FX(i,j)=0.5*( w(i,j,k,nrhs)+w(i-1,j,k,nrhs)
     &                   -0.3333333333333333*(grad(i,j)-grad(i-1,j))
     &                                                 )*Uflxw
#endif
            enddo           !--> discard curv,grad, keep FX
          enddo

          do j=0,ny+2
            do i=1,nx
              FE(i,j)=(w(i,j,k,nrhs)-w(i,j-1,k,nrhs))*vmask(i,j)
            enddo
          enddo

          if (south_bnd) then
            do i=1,nx
              FE(i,0)=FE(i,1)
            enddo
          endif
          if (north_bnd) then
            do i=1,nx
              FE(i,ny+2)=FE(i,ny+1)
            enddo
          endif

          do j=0,ny+1
            do i=1,nx
#if defined UPSTREAM_W
              curv(i,j)=FE(i,j+1)-FE(i,j)
#else
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
#endif
            enddo
          enddo            !--> discard FE

          do j=1,ny+1
            do i=1,nx
              if (k<nz) then
                Vflxw= 0.5*(FlxV(i,j,k)+FlxV(i,j,k+1))
              else
                Vflxw= 0.5*(FlxV(i,j,k)              )
              endif
#ifdef UPSTREAM_W
              FE(i,j)=0.5*(w(i,j,k,nrhs)+w(i,j-1,k,nrhs))
     &                                                  *Vflxw
     &          -0.1666666666666666*( curv(i,j-1)*max(Vflxw,0.)
     &                               +curv(i,j  )*min(Vflxw,0.))
#else
              FE(i,j)=0.5*( w(i,j,k,nrhs)+w(i,j-1,k,nrhs)
     &                   -0.3333333333333333*(grad(i,j)-grad(i,j-1))
     &                                                 )*Vflxw
#endif
            enddo
          enddo             !--> discard curv,grad, keep FE

        do j=1,ny
          do i=1,nx
            rw(i,j,k)=rw(i,j,k)-(FX(i+1,j)-FX(i,j))
     &                         -(FE(i,j+1)-FE(i,j))
#ifdef DIAGNOSTICS_NHMG
            Wdiag(i,j,k,iwhoriadv)=rw(i,j,k)-Wdiag(i,j,k,iwprsgr) ! loop also 1:nz. Only called if NHMG anyway
#endif
          enddo
        enddo
#endif

#undef curv
