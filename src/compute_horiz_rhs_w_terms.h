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
#ifndef EW_PERIODIC
          if (WESTERN_EDGE) then       ! Determine extended index
            imin=istr                  ! range for computation of
          else                         ! elementary differences: it
            imin=istr-1                ! needs to be restricted
          endif                        ! because in the vicinity of
          if (EASTERN_EDGE) then       ! physical boundary the extra
            imax=iend                  ! point may be not available,
          else                         ! and extrapolation of slope
            imax=iend+1                ! is used instead.
          endif
#else
          imin=istr-1
          imax=iend+1
#endif
          do j=jstr,jend
            do i=imin,imax+1
              FX(i,j)=(w(i,j,k,nrhs)-w(i-1,j,k,nrhs))
#ifdef MASKING
     &                                               *umask(i,j)
#endif
            enddo
          enddo
#ifndef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=jstr,jend
              FX(istr-1,j)=FX(istr,j)
            enddo
          endif
          if (EASTERN_EDGE) then
            do j=jstr,jend
              FX(iend+2,j)=FX(iend+1,j)
            enddo
          endif
#endif
          do j=jstr,jend
            do i=istr-1,iend+1
#if defined UPSTREAM_W
              curv(i,j)=FX(i+1,j)-FX(i,j)
#else
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
#endif
            enddo
          enddo             !--> discard FX
          do j=jstr,jend
            do i=istr,iend+1
              if (k<N) then
                Uflxw= 0.5*(FlxU(i,j,k)+FlxU(i,j,k+1))
              else
                Uflxw= 0.5*(FlxU(i,j,k)              )
              endif
#ifdef UPSTREAM_W
              FX(i,j)=0.5*(w(i,j,k,nrhs)+w(i-1,j,k,nrhs))
     &                                       *Uflxw
     &          -0.1666666666666666*( curv(i-1,j)*max(Uflxw,0.)
     &                               +curv(i  ,j)*min(Uflxw,0.))
!# ifdef DIAGNOSTICS_W
!              TruncFX(i,j)=0.04166666666666667*(curv(i,j)-curv(i-1,j))
!     &                                               *abs(Uflxw)
!# endif
#else
              FX(i,j)=0.5*( w(i,j,k,nrhs)+w(i-1,j,k,nrhs)
     &                   -0.3333333333333333*(grad(i,j)-grad(i-1,j))
     &                                                 )*Uflxw
#endif
            enddo           !--> discard curv,grad, keep FX
          enddo

#ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            jmin=jstr
          else
            jmin=jstr-1
          endif
          if (NORTHERN_EDGE) then
            jmax=jend
          else
            jmax=jend+1
          endif
#else
          jmin=jstr-1
          jmax=jend+1
#endif
          do j=jmin,jmax+1
            do i=istr,iend
              FE(i,j)=(w(i,j,k,nrhs)-w(i,j-1,k,nrhs))
#ifdef MASKING
     &                                               *vmask(i,j)
#endif
            enddo
          enddo
#ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=istr,iend
              FE(i,jstr-1)=FE(i,jstr)
            enddo
          endif
          if (NORTHERN_EDGE) then
            do i=istr,iend
              FE(i,jend+2)=FE(i,jend+1)
            enddo
          endif
#endif
          do j=jstr-1,jend+1
            do i=istr,iend
#if defined UPSTREAM_W
              curv(i,j)=FE(i,j+1)-FE(i,j)
#else
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
#endif
            enddo
          enddo            !--> discard FE

          do j=jstr,jend+1
            do i=istr,iend
              if (k<N) then
                Vflxw= 0.5*(FlxV(i,j,k)+FlxV(i,j,k+1))
              else
                Vflxw= 0.5*(FlxV(i,j,k)              )
              endif
#ifdef UPSTREAM_W
              FE(i,j)=0.5*(w(i,j,k,nrhs)+w(i,j-1,k,nrhs))
     &                                                  *Vflxw
     &          -0.1666666666666666*( curv(i,j-1)*max(Vflxw,0.)
     &                               +curv(i,j  )*min(Vflxw,0.))
!# ifdef DIAGNOSTICS_W
!              TruncFE(i,j)=0.04166666666666667*(curv(i,j)-curv(i-1,j))
!     &                                               *abs(Vflxw)
!# endif
#else
              FE(i,j)=0.5*( w(i,j,k,nrhs)+w(i,j-1,k,nrhs)
     &                   -0.3333333333333333*(grad(i,j)-grad(i,j-1))
     &                                                 )*Vflxw
#endif
            enddo
          enddo             !--> discard curv,grad, keep FE

        do j=jstr,jend
          do i=istr,iend
            rw(i,j,k)=rw(i,j,k)-(FX(i+1,j)-FX(i,j))
     &                         -(FE(i,j+1)-FE(i,j))
#ifdef DIAGNOSTICS_NHMG
            Wdiag(i,j,k,iwhoriadv)=rw(i,j,k)-Wdiag(i,j,k,iwprsgr) ! loop also 1:N. Only called if NHMG anyway
#endif
          enddo
        enddo
#endif

#undef curv
