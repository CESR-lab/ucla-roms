! Compute horizontal fluxes for tracers.  Essentially interpolate
! tracer values from their native locations on C grid to horizontal
! velocity points with simultaneous translation from grid-box-averages
! to instantaneous values at interface location.  Three options an be
! selected: 3-point upstream-biased parabolic interpolation (UPSTREAM_TS);
!  4-point symmetric fourth-order method (undefined state of both CPP
! switches); 4-point scheme where arithmetic averaging of elementary
! differences is replaced by harmonic averaging (AKIMA), resulting in
! mid-point values be bounded by two nearest values at native location,
! regardless of grid-scale roughness of the interpolated field, while
! still retaining asymptotic fourth-order behavior for smooth fields.
! This code is extracted into a special module because it is used
! twice, in predictor and corrector substeps for tracer variables.


c---#define BIO_1ST_USTREAM_TEST
#ifdef  BIO_1ST_USTREAM_TEST
        if (itrc>isalt) then       ! biological tracer components:
!       if (itrc>0) then       ! biological tracer components: 
          if (nrhs==3) then         ! compute fluxes during corrector
            do j=jstr,jend            ! stage only
              do i=istr,iend+1
                FX(i,j)=t(i-1,j,k,nstp,itrc)*max(FlxU(i,j,k),0.)
     &                 +t(i  ,j,k,nstp,itrc)*min(FlxU(i,j,k),0.)
              enddo
            enddo
            do j=jstr,jend+1
              do i=istr,iend
                FE(i,j)=t(i,j-1,k,nstp,itrc)*max(FlxV(i,j,k),0.)
     &                 +t(i,j  ,k,nstp,itrc)*min(FlxV(i,j,k),0.)
              enddo
            enddo
          else                         ! there is no need to compute
            do j=jstr,jend+1           ! fluxes during predictor stage
              do i=istr,iend+1         ! because there is no use for
                FX(i,j)=0.             ! t(:,:,:,:,n+1/2) in the case
                FE(i,j)=0.             ! of 1st-order upsteam (note
              enddo                    ! index "nstp" instead of
            enddo                      ! "nrhs" above.
          endif
        else       !--> standard code applies for T,S
#endif


#ifdef UPSTREAM_TS
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
              FX(i,j)=(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
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
#if defined UPSTREAM_TS
              curv(i,j)=FX(i+1,j)-FX(i,j)
#elif defined AKIMA
              cff=2.*FX(i+1,j)*FX(i,j)
              if (cff>epsil) then
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              else
                grad(i,j)=0.
              endif
#else
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
#endif
            enddo
          enddo             !--> discard FX
          do j=jstr,jend
            do i=istr,iend+1
#ifdef UPSTREAM_TS
              FX(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc))
     &                                                  *FlxU(i,j,k)
     &          -0.1666666666666666*( curv(i-1,j)*max(FlxU(i,j,k),0.)
     &                               +curv(i  ,j)*min(FlxU(i,j,k),0.))
#else
              FX(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc)
     &                   -0.3333333333333333*(grad(i,j)-grad(i-1,j))
     &                                                 )*FlxU(i,j,k)
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
              FE(i,j)=(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
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
#if defined UPSTREAM_TS
              curv(i,j)=FE(i,j+1)-FE(i,j)
#elif defined AKIMA
              cff=2.*FE(i,j+1)*FE(i,j)
              if (cff>epsil) then
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              else
                grad(i,j)=0.
              endif
#else
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
#endif
            enddo
          enddo            !--> discard FE

          do j=jstr,jend+1
            do i=istr,iend
#ifdef UPSTREAM_TS
              FE(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc))
     &                                                  *FlxV(i,j,k)
     &          -0.1666666666666666*( curv(i,j-1)*max(FlxV(i,j,k),0.)
     &                               +curv(i,j  )*min(FlxV(i,j,k),0.))
# ifdef DIAGNOSTICS_TS
              TruncFE(i,j)=0.04166666666666667*(curv(i,j)-curv(i-1,j))
     &                                               *abs(FlxV(i,j,k))
# endif
#else
              FE(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc)
     &                   -0.3333333333333333*(grad(i,j)-grad(i,j-1))
     &                                                 )*FlxV(i,j,k)
#endif
            enddo
          enddo             !--> discard curv,grad, keep FE
#ifdef BIO_1ST_USTREAM_TEST
        endif  !<-- itrc>isalt, bio-components only.
#endif

        if (river_source) then
          !! inefficient because this is inside a k-loop
          !! we could try to compute riv_uvel(i,j) somewhere else
          do j=jstr,jend
            do i=istr,iend+1
              if (abs(riv_uflx(i,j)).gt.1e-3) then
                riv_depth = 0.5*( z_w(i-1,j,N)-z_w(i-1,j,0)
     &                      + z_w(i  ,j,N)-z_w(i  ,j,0) )
                iriver = nint(riv_uflx(i,j)/10)
                riv_uvel = riv_vol(iriver)*(riv_uflx(i,j)-10*iriver)/riv_depth
                FX(i,j)= riv_trc(iriver,itrc)*
     &            0.5*(Hz(i-1,j,k)+Hz(i,j,k))*riv_uvel

              endif
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              if (abs(riv_vflx(i,j)).gt.1e-3) then
                riv_depth = 0.5*( z_w(i,j-1,N)-z_w(i,j-1,0)
     &                      + z_w(i  ,j,N)-z_w(i  ,j,0) )
                iriver = nint(riv_vflx(i,j)/10)
                riv_vvel = riv_vol(iriver)*(riv_vflx(i,j)-10*iriver)/riv_depth
                FE(i,j)= riv_trc(iriver,itrc)*
     &            0.5*(Hz(i,j-1,k)+Hz(i,j,k))*riv_vvel

              endif
            enddo
          enddo
        endif  !<-- river_source
