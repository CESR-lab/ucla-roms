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
            do j=1,ny            ! stage only
              do i=1,nx+1
                UFx(i,j)=t(i-1,j,k,nstp,itrc)*max(FlxU(i,j,k),0.)
     &                 +t(i  ,j,k,nstp,itrc)*min(FlxU(i,j,k),0.)
              enddo
            enddo
            do j=1,ny+1
              do i=1,nx
                VFe(i,j)=t(i,j-1,k,nstp,itrc)*max(FlxV(i,j,k),0.)
     &                 +t(i,j  ,k,nstp,itrc)*min(FlxV(i,j,k),0.)
              enddo
            enddo
          else                         ! there is no need to compute
            do j=1,ny+1           ! fluxes during predictor stage
              do i=1,nx+1         ! because there is no use for
                UFx(i,j)=0.             ! t(:,:,:,:,n+1/2) in the case
                VFe(i,j)=0.             ! of 1st-order upsteam (note
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
          if (.not.(ew_periodic)) then
            if (inode==0) then       ! Determine extended index
              imin=1                  ! range for computation of
            else                         ! elementary differences: it
              imin=0                ! needs to be restricted
            endif                        ! because in the vicinity of
            if (inode==np_xi-1) then       ! physical boundary the extra
              imax=nx                  ! point may be not available,
            else                         ! and extrapolation of slope
              imax=nx+1                ! is used instead.
            endif
          else
            imin=0
            imax=nx+1
          endif
          wrk3(1-bf:nx+bf,1-bf:ny+bf) => scratch_2d(:,7)
          do j=1,ny
            do i=imin,imax+1
              wrk3(i,j)=(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
#ifdef MASKING
     &                                               *umask(i,j)
#endif
            enddo
          enddo
          if (.not.(ew_periodic)) then
            if (inode==0) then
              do j=1,ny
                wrk3(0,j)=wrk3(1,j) ! issue
              enddo
            endif
            if (inode==np_xi-1) then
              do j=1,ny
                wrk3(nx+2,j)=wrk3(nx+1,j)
              enddo
            endif
          endif
          do j=1,ny
            do i=0,nx+1
#if defined UPSTREAM_TS
              curv(i,j)=wrk3(i+1,j)-wrk3(i,j)
#elif defined AKIMA
              cff=2.*wrk3(i+1,j)*wrk3(i,j)
              if (cff>epsil) then
                grad(i,j)=cff/(wrk3(i+1,j)+wrk3(i,j))
              else
                grad(i,j)=0.
              endif
#else
              grad(i,j)=0.5*(wrk3(i+1,j)+wrk3(i,j))
#endif
            enddo
          enddo             !--> discard UFx
          nullify(wrk3)
          do j=1,ny
            do i=1,nx+1
#ifdef UPSTREAM_TS
              UFx(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc))
     &                                                  *FlxU(i,j,k)
     &          -0.1666666666666666*( curv(i-1,j)*max(FlxU(i,j,k),0.)
     &                               +curv(i  ,j)*min(FlxU(i,j,k),0.))
#else
              UFx(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc)
     &                   -0.3333333333333333*(grad(i,j)-grad(i-1,j))
     &                                                 )*FlxU(i,j,k)
#endif
            enddo           !--> discard curv,grad, keep UFx
          enddo

          if (.not.(ns_periodic)) then
            if (jnode==0) then
              jmin=1
            else
              jmin=0
            endif
            if (jnode==np_eta-1) then
              jmax=ny
            else
              jmax=ny+1
            endif
          else
            jmin=0
            jmax=ny+1
          endif
          wrk3(1-bf:nx+bf,1-bf:ny+bf) => scratch_2d(:,7)
          do j=jmin,jmax+1
            do i=1,nx
              wrk3(i,j)=(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
#ifdef MASKING
     &                                               *vmask(i,j)
#endif
            enddo
          enddo

          if (.not.(ns_periodic)) then
            if (jnode==0) then
              do i=1,nx
                wrk3(i,0)=wrk3(i,1) ! issue
              enddo
            endif
            if (jnode==np_eta-1) then
              do i=1,nx
                wrk3(i,ny+2)=wrk3(i,ny+1)
              enddo
            endif
          endif
          do j=0,ny+1
            do i=1,nx
#if defined UPSTREAM_TS
              curv(i,j)=wrk3(i,j+1)-wrk3(i,j)
#elif defined AKIMA
              cff=2.*wrk3(i,j+1)*wrk3(i,j)
              if (cff>epsil) then
                grad(i,j)=cff/(wrk3(i,j+1)+wrk3(i,j))
              else
                grad(i,j)=0.
              endif
#else
              grad(i,j)=0.5*(wrk3(i,j+1)+wrk3(i,j))
#endif
            enddo
          enddo            !--> discard wrk3/VFe

	  nullify(wrk3)

          do j=1,ny+1
            do i=1,nx
#ifdef UPSTREAM_TS
              VFe(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc))
     &                                                  *FlxV(i,j,k)
     &          -0.1666666666666666*( curv(i,j-1)*max(FlxV(i,j,k),0.)
     &                               +curv(i,j  )*min(FlxV(i,j,k),0.))
# ifdef DIAGNOSTICS_TS
              TruncFE(i,j)=0.04166666666666667*(curv(i,j)-curv(i-1,j))
     &                                               *abs(FlxV(i,j,k))
# endif
#else
              VFe(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc)
     &                   -0.3333333333333333*(grad(i,j)-grad(i,j-1))
     &                                                 )*FlxV(i,j,k)
#endif
            enddo
          enddo             !--> discard curv,grad, keep VFe
#ifdef BIO_1ST_USTREAM_TEST
        endif  !<-- itrc>isalt, bio-components only.
#endif

        if (river_source) then
          !! inefficient because this is inside a k-loop
          !! we could try to compute riv_uvel(i,j) somewhere else
          do j=1,ny
            do i=1,nx+1
              if (abs(riv_uflx(i,j)).gt.1e-3) then
                riv_depth = 0.5*( z_w(i-1,j,N)-z_w(i-1,j,0)
     &                      + z_w(i  ,j,N)-z_w(i  ,j,0) )
                iriver = nint(riv_uflx(i,j)/10)
                riv_uvel = riv_vol(iriver)*(riv_uflx(i,j)-10*iriver)/riv_depth
                UFx(i,j)= riv_trc(iriver,itrc)*
     &            0.5*(Hz(i-1,j,k)+Hz(i,j,k))*riv_uvel

              endif
            enddo
          enddo
          do j=1,ny+1
            do i=1,nx
              if (abs(riv_vflx(i,j)).gt.1e-3) then
                riv_depth = 0.5*( z_w(i,j-1,N)-z_w(i,j-1,0)
     &                      + z_w(i  ,j,N)-z_w(i  ,j,0) )
                iriver = nint(riv_vflx(i,j)/10)
                riv_vvel = riv_vol(iriver)*(riv_vflx(i,j)-10*iriver)/riv_depth
                VFe(i,j)= riv_trc(iriver,itrc)*
     &            0.5*(Hz(i,j-1,k)+Hz(i,j,k))*riv_vvel

              endif
            enddo
          enddo
        endif  !<-- river_source
