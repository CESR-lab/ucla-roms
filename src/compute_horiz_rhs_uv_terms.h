#if defined UV_COR || (defined CURVGRID && defined UV_ADV)
        do j=jstrV-1,jend                ! Add Coriolis terms and
          do i=istrU-1,iend              ! contribution to advection
            cff=0.5*Hz(i,j,k)*(          ! associated with curvilinear
# ifdef UV_COR
     &              fomn(i,j)            ! horizontal coordinates.
# endif
# if (defined CURVGRID && defined UV_ADV)
     &             +0.5*( (v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)
     &                   -(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
# endif
     &                                                             )
# ifdef WEC
#  if (defined CURVGRID && defined UV_ADV)
     ! Define a cff1 that is added to UFx, VFe for stokes terms
           cff1=0.5*Hz(i,j,k)*(
     &          0.5*( dndx(i,j)*(vst(i,j,k)+vst(i,j+1,k))
     &               -dmde(i,j)*(ust(i,j,k)+ust(i+1,j,k)) ))
#  else
           cff1 = 0.0
#  endif
           UFx(i,j)=(cff+cff1)*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
           UFe(i,j)=cff*(vst(i,j,k)+vst(i,j+1,k))
           VFe(i,j)=(cff+cff1)*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
           VFx(i,j)=cff*(ust(i,j,k)+ust(i+1,j,k))
# else /* no WEC */

            UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
# endif
          enddo
        enddo
        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)+0.5*(UFx(i,j)+UFx(i-1,j))
# ifdef WEC
     &               + 0.5*(UFe(i,j)+UFe(i-1,j))
# endif
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-0.5*(VFe(i,j)+VFe(i,j-1))
# ifdef WEC
     &               -0.5*(VFx(i,j)+VFx(i,j-1))
# endif
          enddo
        enddo

# ifdef DIAGNOSTICS
	! JM It would be better to not do the u,v boundary points
        if (CORR_STAGE) then
          if (diag_uv.and.calc_diag) then
            Udiag(:,:,k,icori) = 0.5*(UFx(1:nx,1:ny)+UFx(0:nx-1,1:ny))*dxdyi_u
            Vdiag(:,:,k,icori) =-0.5*(VFe(1:nx,1:ny)+VFe(1:nx,0:ny-1))*dxdyi_v
          endif
        endif
# endif /* DIAGNOSTICS */

#endif  /* UV_COR */



#ifdef UV_ADV

! Add horizontal advection of momentum: compute diagonal [UFx,VFe]
! and off-diagonal [UFe,VFx] components of momentum flux tensor due
! to horizontal advection; after that add their divergence to r.h.s.

# define uxx wrk1
# define Huxx wrk2
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then        ! Sort out bounding indices for
          imin=istrU                  ! the extended ranges: note that
        else                          ! in the vicinity of physical
          imin=istrU-1                ! boundaries values at the
        endif                         ! extremal points of stencil
        if (EASTERN_EDGE) then        ! are not available, so an
          imax=iend                   ! extrapolation rule needs to
        else                          ! be applied.  Also note that
          imax=iend+1                 ! for this purpose periodic
        endif                         ! ghost points and MPI margins
# else
        imin=istr-1                   ! are not considered as
        imax=iend+1                   ! physical boundaries.
# endif
        do j=jstr,jend
          do i=imin,imax
            uxx(i,j)=u(i-1,j,k,nrhs)-2.*u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
            Huxx(i,j)=FlxU(i-1,j,k) -2.*FlxU(i,j,k) +FlxU(i+1,j,k)
          enddo
        enddo
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr,jend
            uxx(istrU-1,j) =uxx(istrU,j)
            Huxx(istrU-1,j)=Huxx(istrU,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
            uxx(iend+1,j) =uxx(iend,j)
            Huxx(iend+1,j)=Huxx(iend,j)
          enddo
        endif
# endif
        do j=jstr,jend
          do i=istrU-1,iend
# ifdef UPSTREAM_UV
            cff=FlxU(i,j,k)+FlxU(i+1,j,k)-delta*( Huxx(i  ,j)
     &                                            +Huxx(i+1,j))
            UFx(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*uxx(i  ,j)
     &                                  +min(cff,0.)*uxx(i+1,j)
     &                                                      ))
# else
            UFx(i,j)=0.25*( u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
     &                         -delta*(uxx(i,j)+uxx(i+1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i+1,j,k)
     &                      -delta*(Huxx(i,j)+Huxx(i+1,j)))
# endif
          enddo
        enddo
# undef Huxx
# undef uxx

# define vee wrk1
# define Hvee wrk2
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=jstrV
        else
          jmin=jstrV-1
        endif
        if (NORTHERN_EDGE) then
          jmax=jend
        else
          jmax=jend+1
        endif
# else
        jmin=jstr-1
        jmax=jend+1
# endif
        do j=jmin,jmax
          do i=istr,iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.*v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
            Hvee(i,j)=FlxV(i,j-1,k) -2.*FlxV(i,j,k)  +FlxV(i,j+1,k)
          enddo
        enddo
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr,iend
            vee(i,jstrV-1)=vee(i,jstrV)
            Hvee(i,jstrV-1)=Hvee(i,jstrV)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
            vee(i,jend+1)=vee(i,jend)
            Hvee(i,jend+1)=Hvee(i,jend)
          enddo
        endif
# endif
        do j=jstrV-1,jend
          do i=istr,iend
# ifdef UPSTREAM_UV
            cff=FlxV(i,j,k)+FlxV(i,j+1,k)-delta*( Hvee(i,j  )
     &                                           +Hvee(i,j+1))
            VFe(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
     &                          -gamma*( max(cff,0.)*vee(i,j  )
     &                                  +min(cff,0.)*vee(i,j+1)
     &                                                      ))
# else
            VFe(i,j)=0.25*( v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
     &                        -delta*(vee(i,j)+vee(i,j+1))
     &                   )*( FlxV(i,j,k)+FlxV(i,j+1,k)
     &                      -delta*(Hvee(i,j)+Hvee(i,j+1)))
# endif
          enddo
        enddo
# undef Hvee
# undef vee

# define uee wrk1
# define Hvxx wrk2
# ifndef NS_PERIODIC
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
# else
        jmin=jstr-1
        jmax=jend+1
# endif
        do j=jmin,jmax
          do i=istrU,iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.*u(i,j,k,nrhs)+u(i,j+1,k,nrhs)
          enddo
        enddo
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istrU,iend
            uee(i,jstr-1)=uee(i,jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istrU,iend
            uee(i,jend+1)=uee(i,jend)
          enddo
        endif
# endif
        do j=jstr,jend+1
          do i=istrU-1,iend
           Hvxx(i,j)=FlxV(i-1,j,k)-2.*FlxV(i,j,k)+FlxV(i+1,j,k)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istrU,iend
# ifdef UPSTREAM_UV
            cff=FlxV(i,j,k)+FlxV(i-1,j,k)-delta*( Hvxx(i  ,j)
     &                                           +Hvxx(i-1,j))
            UFe(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i,j-1,k,nrhs))
     &                          -gamma*( max(cff,0.)*uee(i,j-1)
     &                                  +min(cff,0.)*uee(i,j  )
     &                                                      ))
# else
            UFe(i,j)=0.25*( u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
     &                        -delta*(uee(i,j)+uee(i,j-1))
     &                  )*( FlxV(i,j,k)+FlxV(i-1,j,k)
     &                     -delta*(Hvxx(i,j)+Hvxx(i-1,j)))
# endif
          enddo
        enddo
# undef Hvxx
# undef uee

# define vxx wrk1
# define Huee wrk2
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          imin=istr
        else
          imin=istr-1
        endif
        if (EASTERN_EDGE) then
          imax=iend
        else
          imax=iend+1
        endif
# else
        imin=istr-1
        imax=iend+1
# endif
        do j=jstrV,jend
          do i=imin,imax
            vxx(i,j)=v(i-1,j,k,nrhs)-2.*v(i,j,k,nrhs)+v(i+1,j,k,nrhs)
          enddo
        enddo
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstrV,jend
            vxx(istr-1,j)=vxx(istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstrV,jend
            vxx(iend+1,j)=vxx(iend,j)
          enddo
        endif
# endif
        do j=jstrV-1,jend
          do i=istr,iend+1
           Huee(i,j)=FlxU(i,j-1,k)-2.*FlxU(i,j,k)+FlxU(i,j+1,k)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend+1
# ifdef UPSTREAM_UV
            cff=FlxU(i,j,k)+FlxU(i,j-1,k)-delta*( Huee(i,j  )
     &                                           +Huee(i,j-1))
            VFx(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i-1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*vxx(i-1,j)
     &                                  +min(cff,0.)*vxx(i  ,j)
     &                                                      ))
# else
            VFx(i,j)=0.25*( v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
     &                        -delta*(vxx(i,j)+vxx(i-1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i,j-1,k)
     &                     -delta*(Huee(i,j)+Huee(i,j-1)))
# endif
          enddo
        enddo
# undef Huee
# undef vxx
        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)-UFx(i,j  )+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i  ,j)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j  )
     &                         -VFe(i  ,j)+VFe(i,j-1)
          enddo
        enddo
#endif /* UV_ADV */
