! Content of this package:
!------------------------- ! initialization of vertical interpolation
! r2r_init_vertint_thread  ! by inverse mapping by cubic splines: given
! r2r_check_vertint_thread ! parent-grid zp_r(kp=1,Np) and child-grid
! r2r_init_vrtint_tile     ! z_r(k=1,N) compute parent grid index-space
! r2r_check_vrtint_tile    ! coordinate kprnt(k=1,N) such that
! comp_zspline_tile        !      zspline[kp=kprnt(k)] = z_r(k)
!
!                          ! interpolation routines for fields:
! r2r_vrtint_thread        ! compute vertical derivatives in parent
! r2r_vertint_tile          ! grid index-coordinate space, then knowing
! r2r_vsplnint_tile        ! kp=kprnt(k=1,N), interpolate parent-grid
!                          ! data onto child grid.
!
! r2r_set_depth_tile       ! compute z_r=z_r[h(j,j), Cs(k)]
!
! compute_uvbar_thread     ! vertical integration of u,v-velocity
! compute_uvbar_tile       ! components to compute barotripic mode

#define MASK_LAND_DATA
#define SPLINE
#define VERBOSE
#undef NEWTON

#ifndef MASK_LAND_DATA
# define FillValue 0.
#endif


      module r2r_vertint_vars
       integer, save :: trd_count=0
       real(kind=8), save :: vert_int_error=0.
       integer, save :: alloc_zc_size=0, allc_zpr_size=0
C$OMP THREADPRIVATE(allc_zpr_size, alloc_zc_size)
       real(kind=8), allocatable, dimension(:) :: zc,zp_r,drv
C$OMP THREADPRIVATE(zc,zp_r,drv)
      end module r2r_vertint_vars



      subroutine r2r_init_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     &                                            h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer ncx,ncy, Np,N
      real(kind=8) hprnt(ncx,ncy),h(ncx,ncy),hcp,hc,Csp_r(Np),Cs_r(N),
     &                                                kprnt(ncx,ncy,N)
      integer ntrds,trd, nsub_x,nsub_y, ntls, isize,
     &        my_first,my_last,tile, istr,iend,jstr,jend
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zc_size < isize*(N+1) .or.
     &    allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) deallocate(zp_r,drv,zc)
        allocate( zc(alloc_zc_size), zp_r(allc_zpr_size),
     &                                drv(allc_zpr_size) )
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)')  'r2r_init_vertint_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     &                           'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
#ifdef VERBOSE
      elseif (trd == 0) then
        write(*,'(2x,A)',advance='no')
     & 'entering r2r_init_vertint_thread ...'
#endif
      endif

      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     &                                    istr,iend,jstr,jend )
        call r2r_init_vrtint_tile( istr,iend,jstr,jend, ncx,ncy,
     &                              hprnt,Np,hcp,Csp_r,zp_r,drv,
     &                                 h, N,hc,Cs_r,zc, kprnt)
      enddo
#ifdef VERBOSE
      if (trd == 0) write(*,'(2x,A)')
     & 'leaving r2r_init_vertint_thread'
#endif
C$OMP BARRIER
      end

      subroutine r2r_check_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     &                                             h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer ncx,ncy, Np,N
      real(kind=8) hprnt(ncx,ncy),h(ncx,ncy), hcp,hc,Csp_r(Np),Cs_r(N),
     &                                        kprnt(ncx,ncy,N),my_error
      integer ntrds,trd, nsub_x,nsub_y, ntls, isize,
     &        my_first,my_last, tile, istr,iend,jstr,jend
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zc_size < isize*(N+1) .or.
     &    allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) deallocate(zp_r,drv,zc)
        allocate(zc(alloc_zc_size), zp_r(allc_zpr_size),
     &                               drv(allc_zpr_size))
      endif
      my_error=0.D0 !<-- initialize
      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     &                                    istr,iend,jstr,jend )
        call r2r_check_vrtint_tile(istr,iend,jstr,jend, ncx,ncy,
     &                             hprnt, Np,hcp,Csp_r,zp_r,drv,
     &                        h, N,hc,Cs_r,zc, kprnt, my_error)
      enddo
C$OMP CRITICAL(r2r_vert_crgn)
      if (trd_count == 0) vert_int_error=0.D0
      trd_count=trd_count+1
      vert_int_error=max(vert_int_error, my_error)
      if (trd_count == ntrds) then
        trd_count=0
        write(*,*)
        write(*,*) 'maximum vert_int_error =', vert_int_error
        write(*,*)
      endif
C$OMP END CRITICAL(r2r_vert_crgn)
      end



      subroutine r2r_init_vrtint_tile( istr,iend,jstr,jend, ncx,ncy,
     &                                 hprnt,Np,hcp,Csp_r, zp_r,drv,
     &                                     h,N, hc, Cs_r,z_r, kprnt)

! Given a set of parent-grid z-levels zp_r=z_r(kp), kp=1,..,Np defined
! at vertical RHO-points, as well as knowing that zp=zeta at surface
! kp=Np+1/2, (always presume zeta=0 for this purpose), and zp=-h_parent
! at bottom, kp=1/2, construct

! (i) a set of spine derivatives drv=drv(kp) = d zp_r/dkp (vertical
!     derivative in index space of parent grid) defined at the same
!     RHO-points, then

! (ii) construct an inverse mapping function
!
!           kprnt=kprnt(z), 1/2 < kprnt < Np+1/2
!
!      which maps child-grid z-levels z_r(k), k=1,..,N  onto
!      parent-grid vertical index kp in such a way that index-space
!      spline-interpolated value
!
!         z* = spline[p, zp_r(kp), zp_r(kp+1), drv(kp), drv(kp+1)]
!
!      where index "kp" in the integer part of kprnt(z) and "p" is
!      the fractional reminder, is equal to "z" itself.
!
! Once this mapping  kp = kp(k) = kprnt[z_r(k)] is established, any
! field "qsrc" defined on the parent grid, hence qsrc(kp), kp=1,...,Np
! can be spline interpolated in index space,
!
!      q(k)= spline[p, qsrc(kp), src(kp+1), drvq(kp), drvq(kp+1)]
!
! where  drvq = d qsrc/d kp is vertical derivative of qsrc computed
! in index space of parent grid.

! Task (ii) above is achieved by first using linear interpolation to
! map child-grid z-levels z_r(k) onto parent-grid zp_r(kp). This leads
! to an an initial approximation for kprnt(k). Subsequently use Newton
! iterations to compute fractional offsets "p" in such a way that spline
! version of z_r(k) = spline[p, zp_r(kp),...] item (ii) above, holds.

! The linear stage constructs "continuous index" mapping coordinate
! kp=kprnt(k) such that interpolation of field "qsrc" defined on the
! parent grid z=zp_r(kp), kp=1,...,Np into child grid z=z_r(k),
! k=1,...,N is approximated as
!
!           kp=int(kprnt(k)) ; p = kprnt(k)-float(kp)
!           q(k)=p*qsrc(kp+1) + (1-p)*qsrc(kp)
!
! where, "kp" and "p" are the integer and fractional parts of "kprnt"
! defined in such a way that zp_r(kp), kp=1,...,Np interpolated the
! same way into the location z_r(k) should yield z_r(k) itself, that
! is
!           z_r(k)=p*zp_r(kp+1) + (1-p)*zp_r(kp)
!
! where index "kp" is chosen to make z_r(k) bounded as
!
!           zp_r(kp) <=  z_r(k) < zp_r(kp+1)
!
! and the fractional part
!
!           p = [z_r(k)-zp_r(kp)]/[zp_r(kp+1)-zp_r(kp)]
!
! The caveats are:
!
! (iii) it is possible that some child-grid "z_r" points may be above
!       the uppermost available zp_r(Np); in this case kp=Np while
!       fraction "p" is set using unperturbed free surface z=0 instead
!       of nonexisting zp_r(kp+1); interpolation becomes extrapolation
!       toward surface;
!
! (iv) for a similar reason some "z_r" may be below the lowest
!      available zp_r(kp=1), and even worse, because of inconsistency
!      between parent and child topography some "z_r" points may be
!      even below bottom of the parent grid.

! To address both (or at least mitigate in the case of (iv)) spline
! interpolation algorithm for zp_r (see "comp_zspline_tile" below with
! CPP-switches STAGG_BCS and/or STAGG_NOT_A_KNOT_BCS activated) is
! designed to utilize boundary conditions
!
!               zp(kprnt=Np+1/2) = zeta     !<-- surface
!               zp(kprnt=1/2) = -h_parent   !<-- bottom
!
! in addition to the usual Np equations for spline derivatives: i.e.,
! the extremal kp=1 and kp=Np are no-longer treated as boundaries, but
! are formed from the same assumption of second-derivative continuity
! as for all the interior points. Furthermore, to simplify the searches
! and interpolations "comp_zspline_tile" produces ghost point values
!
!            zp_r(kp=0), drv(0), zp_r(Np+1),  drv(Np+1)
!
! which are at RHO-points half-grid-interval below the bottom and
! above the surface. These ghost points are constructed in such a way
! that spline polynomial of zp_r taken at kprnt=1/2 yields exactly
! -h_parent (bottom) and at kprnt=N+1/2 yields zeta (=0 in all cases
! here, surface).  Then, finally,
!
!        kprnt > Np+1/2  above surface, should never occur
!   Np < kprnt < Np+1/2  extrapolation toward surface; accepted
!   1  < kprnt < Np      interpolation within the vertical column
!  1/2 < kprnt < 1       extrapolation toward bottom;
!        kprnt < 1/2     below bottom of parent grid; use filling


      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, Np,N
      real(kind=8) hprnt(ncx,ncy), h(ncx,ncy), Csp_r(Np), Cs_r(N),
     &             zp_r(istr:iend,0:Np+1), z_r(istr:iend,N), hcp,hc,
     &              drv(istr:iend,0:Np+1), kprnt(ncx,ncy,N), p
      integer i,j,k, kp
#ifdef SPLINE
      integer iter
      real(kind=8) q,pq, zerr, dZds
# ifndef NEWTON
     &                 , dZs2, d2Zds2
# endif
#endif
      do j=jstr,jend
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &                               h, N, hc,Cs_r, z_r)
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &               hprnt, Np,hcp,Csp_r, zp_r(istr,1) )
        do i=istr,iend
          zp_r(i,Np+1)=0.D0     !<-- surface
          zp_r(i,0)=-hprnt(i,j) !<-- bottom
        enddo
#ifdef SPLINE
        call comp_zspline_tile(istr,iend, Np, zp_r,drv)
#endif
        do i=istr,iend          !--> search loop
          kp=Np
          do k=N,1,-1  !--> recursive because of kp
            do while(z_r(i,k) < zp_r(i,kp) .and. kp > 0)
              kp=kp-1
            enddo
            kprnt(i,j,k)=dble(kp)+0.5D0  !<-- temporarily
          enddo
        enddo                            ! set initial approximation
        do k=1,N                         ! for fractional distance "p"
          do i=istr,iend                 ! by linear interpolation then
            kp=int(kprnt(i,j,k))
            p=(z_r(i,k)-zp_r(i,kp))/(zp_r(i,kp+1)-zp_r(i,kp))
            kprnt(i,j,k)=dble(kp)+p
          enddo                          ! apply Newton iterations
#ifdef SPLINE
          do i=istr,iend
            if (kprnt(i,j,k) > 0.499999999999D0) then
              kp=int(kprnt(i,j,k)) ; p=kprnt(i,j,k)-dble(kp)
              do iter=1,8
                q=1.D0-p ; pq=p*q

                zerr=p*((p+2.D0*pq)*zp_r(i,kp+1) -pq*drv(i,kp+1))
     &              +q*((q+2.D0*pq)*zp_r(i,kp  ) +pq*drv(i,kp  ))
     &                                                  -z_r(i,k)
                dZds=6.D0*pq*(zp_r(i,kp+1)-zp_r(i,kp))
     &                           +(p*p-2.D0*pq)*drv(i,kp+1)
     &                           +(q*q-2.D0*pq)*drv(i,kp  )

# ifdef NEWTON
                p=p-zerr/dZds !-> d2Zds2 below is (1/2)*d^2 Z/ds^2
# else
                d2Zds2=(2.D0*p-q)*drv(i,kp+1) -(2.D0*q-p)*drv(i,kp)
     &                        -3.D0*(p-q)*(zp_r(i,kp+1)-zp_r(i,kp))

                dZs2=dZds*dZds
                p=p-zerr*(dZs2+d2Zds2*zerr)/(dZs2*dZds)
# endif
              enddo
              kprnt(i,j,k)=dble(kp)+p !<-- corrected value
            else
              kprnt(i,j,k)=0.D0
            endif
          enddo
#endif
        enddo !<-- k
      enddo !<-- j
      end

      subroutine r2r_check_vrtint_tile(istr,iend,jstr,jend, ncx,ncy,
     &                                 hprnt, Np,hcp,Csp_r, zp_r,drv,
     &                            h, N,hc,Cs_r,z_r, kprnt, my_error)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, Np,N, i,j,k, kp
      real(kind=8) hprnt(ncx,ncy), h(ncx,ncy), Csp_r(Np), Cs_r(N),
     &             zp_r(istr:iend,0:Np+1), z_r(istr:iend,N), hcp,hc,
     &             drv(istr:iend,0:Np+1), kprnt(ncx,ncy,N), my_error
     &           , p
#ifdef SPLINE
     &           , q, pq
#endif
      do j=jstr,jend
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &                               h, N, hc,Cs_r, z_r)
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &               hprnt, Np,hcp,Csp_r, zp_r(istr,1) )
        do i=istr,iend
          zp_r(i,Np+1)=0.D0     !<-- surface
          zp_r(i,0)=-hprnt(i,j) !<-- bottom
        enddo
#ifdef SPLINE
        call comp_zspline_tile(istr,iend, Np, zp_r,drv)
#endif
        do i=istr,iend
          do k=1,N
            kp=max(int(kprnt(i,j,k)), 0) ; p=kprnt(i,j,k)-float(kp)
#ifdef SPLINE
            if (kprnt(i,j,k) > 0.499999999999D0) then
              q=1.D0-p ; pq=p*q
              my_error=max( my_error, abs(
     &             p*((p+2.D0*pq)*zp_r(i,kp+1) -pq*drv(i,kp+1))
     &            +q*((q+2.D0*pq)*zp_r(i,kp  ) +pq*drv(i,kp  ))
     &                                            -z_r(i,k) ) )
            endif
#else
            if (0 < kp .and. kp < Np) then
              my_error=max( my_error, abs(       p*zp_r(i,kp+1)
     &                                     +(1.D0-p)*zp_r(i,kp)
     &                                            -z_r(i,k) ) )
            endif
#endif
          enddo
        enddo
      enddo !<-- j
      end




#ifdef SPLINE
# define STAGG_NOT_A_KNOT_BCS

      subroutine comp_zspline_tile(istr,iend, Np, zp_r,drv)

! Compute spline derivatives of drv=d zp_r/d kp for the parent grid
! vertical coordinate.  The incoming argument zp_r=zp_r(kp) contains
! z-levels at RHO-point for indices k=1,..,Np, AND, IN ADDITION TO
! THAT indices kp=0 and kp=Np+1 contain bottom zp_r(0)=-h_parent and
! surface zp_r(Np+1)=zeta value.  On output zp_r(0) and zp_r(Np+1)
! are converted into EXTRAPOLATED ghost-point values, kind of half-
! grid interval below the bottom and above the surface. Derivatives
! drv(0) and drv(Np+1) are also at ghost points, while all others
! for kp=1,...,Np are for regular RHO-points.

      implicit none
      integer istr,iend, Np
      real(kind=8) zp_r(istr:iend,0:Np+1), drv(istr:iend,0:Np+1)
      integer i,k
      real(kind=8) CF(1:Np-1), cff,cff1

# ifdef STAGG_NOT_A_KNOT_BCS
      CF(1)=1.D0/3.D0 ; cff=8.D0/9.D0
      do i=istr,iend
        drv(i,1)=cff*(zp_r(i,2) -zp_r(i,0))
      enddo
# elif defined STAGG_BCS
      CF(1)=0.25D0 ; cff=CF(1)
      do i=istr,iend
        drv(i,1)=cff*( 3.D0*zp_r(i,2) +zp_r(i,1) -4.D0*zp_r(i,0) )
      enddo
# else
      CF(1)=0.5D0
      do i=istr,iend
        drv(i,1)=1.5D0*(zp_r(i,2)-zp_r(i,1))
      enddo
# endif

      do k=2,Np-1,+1   !--> forward elimination
        cff=1.D0/(4.D0-CF(k-1)) ; CF(k)=cff
        do i=istr,iend
          drv(i,k)=cff*( 3.D0*(zp_r(i,k+1)-zp_r(i,k-1)) -drv(i,k-1))
        enddo
      enddo

# ifdef STAGG_NOT_A_KNOT_BCS
      cff=1.D0/(3.D0-CF(Np-1)) ; cff1=8.D0/3.D0
      do i=istr,iend
        drv(i,Np)=cff*( cff1*(zp_r(i,Np+1)-zp_r(i,Np-1)) -drv(i,Np-1) )
      enddo
# elif defined STAGG_BCS
      cff=1.D0/(4.D0-CF(Np-1))
      do i=istr,iend
        drv(i,Np)=cff*( 4.D0*zp_r(i,Np+1)-zp_r(i,Np)-3.D0*zp_r(i,Np-1)
     &                                                  -drv(i,Np-1) )
      enddo
# else
      cff=1.D0/(2.D0-CF(Np-1))
      do i=istr,iend
        drv(i,Np)=cff*( 3.D0*(zp_r(i,Np)-zp_r(i,Np-1)) -drv(i,Np-1))
      enddo
# endif

      do k=Np-1,1,-1   !--> backsubstitution
        cff=CF(k)
        do i=istr,iend
          drv(i,k)=drv(i,k)-cff*drv(i,k+1)
        enddo
      enddo

      do i=istr,iend
# ifdef STAGG_NOT_A_KNOT_BCS
        drv(i,0)=2.D0*drv(i,2)-drv(i,1) -8.D0*zp_r(i,0)+12.D0*zp_r(i,1)
     &                                                  -4.D0*zp_r(i,2)
        zp_r(i,0)=4.D0*zp_r(i,0)-4.D0*zp_r(i,1)+zp_r(i,2)
     &                         -0.5D0*(drv(i,2)-drv(i,1))

        drv(i,Np+1)=2.D0*drv(i,Np-1)-drv(i,Np)  +8.D0*zp_r(i,Np+1)
     &                          -12.D0*zp_r(i,Np)+4.D0*zp_r(i,Np-1)
        zp_r(i,Np+1)=4.D0*zp_r(i,Np+1)-4.D0*zp_r(i,Np)+zp_r(i,Np-1)
     &                               -0.5D0*(drv(i,Np)-drv(i,Np-1))
# elif defined STAGG_BCS
        drv(i,0)=8.D0*(zp_r(i,1)-zp_r(i,0))-3.D0*drv(i,1)
        zp_r(i,0)=4.D0*zp_r(i,0)-3.D0*zp_r(i,1) +drv(i,1)

        drv(i,Np+1)=8.D0*(zp_r(i,Np+1)-zp_r(i,Np))-3.D0*drv(i,Np)
        zp_r(i,Np+1)=4.D0*zp_r(i,Np+1)-3.D0*zp_r(i,Np) -drv(i,Np)
# endif
      enddo
      end
#endif




      subroutine r2r_vrtint_thread(ncx,ncy, lmsk,mask, btm_bc,
     &                                  Np,qsrc, N,kprnt,qtr)
      use r2r_vertint_vars
      implicit none
      integer ncx,ncy, lmsk, btm_bc, Np,N
      integer(kind=2) mask(ncx,ncy)
      real(kind=8) kprnt(ncx,ncy,N)
      real(kind=4) qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
!>
      integer ntrds,trd, nsub_x,nsub_y, ntls, my_first,my_last,
     &                                tile, istr,iend,jstr,jend
C$    integer omp_get_thread_num, omp_get_num_threads

      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      do tile=my_first,my_last,+1
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                  istr,iend,jstr,jend)
#ifdef SPLINE
        call r2r_vsplnint_tile( istr,iend,jstr,jend, ncx,ncy,
     &           lmsk,mask, btm_bc, Np,qsrc,drv, N,kprnt,qtr)
#else
        call r2r_vertint_tile(  istr,iend,jstr,jend, ncx,ncy,
     &           lmsk,mask,         Np,qsrc,     N,kprnt,qtr)
#endif
      enddo
      end

      subroutine r2r_vretint_tile(istr,iend,jstr,jend, ncx,ncy,
     &                         lmsk,mask, Np,qsrc, N,kprnt,qtr)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, lmsk,Np,N, i,j,k,kp
      integer(kind=2) mask(ncx,ncy)
      real(kind=8) kprnt(ncx,ncy,N), p
      real(kind=4) qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
#ifdef MASK_LAND_DATA
# include "spval.h"
#endif
      do j=jstr,jend
        do k=1,N
          do i=istr,iend
            kp=int(kprnt(i,j,k)) ; p=kprnt(i,j,k)-float(kp)
            if (kp > Np-1) then
              qtr(i,j,k)=qsrc(i,j,Np)
            elseif (kp > 0) then
              qtr(i,j,k)=p*qsrc(i,j,kp+1)+(1.D0-p)*qsrc(i,j,kp)
            else
              qtr(i,j,k)=qsrc(i,j,1)
            endif
          enddo
        enddo
        if (lmsk == 1) then
          do i=istr,iend
            if (mask(i,j) == 0) then
              do k=1,N
                qtr(i,j,k)=FillValue
              enddo
            endif
          enddo
        endif
      enddo
      end


#ifdef SPLINE
      subroutine r2r_vsplnint_tile(istr,iend,jstr,jend, ncx,ncy, lmsk,
     &                          mask, btm_bc,Np,qsrc,drv, N,kprnt,qtr)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, lmsk,btm_bc,Np,N, i,j,k,kp
      integer(kind=2) mask(ncx,ncy)
      real(kind=4) qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
      real(kind=8) kprnt(ncx,ncy,N), CF(Np),  drv(istr:iend,Np),
     &                                              cff, p,q,pq
#ifdef MASK_LAND_DATA
# include "spval.h"
#endif
      do j=jstr,jend
        if (btm_bc == -1) then             ! staggered no-slip b.c.
          CF(1)=1.D0/3.D0 ; cff=8.D0/9.D0  ! combined with not-a-knot
          do i=istr,iend                   ! at botton
            drv(i,1)=cff*qsrc(i,j,2)
          enddo
        elseif (btm_bc == +1) then         ! staggered Neumann b.c.
          CF(1)=7.D0/15.D0 ; cff=1.2D0     ! combined with not-a-knot
          do i=istr,iend
            drv(i,1)=cff*(qsrc(i,j,2)-qsrc(i,j,1))
          enddo
        else
          CF(1)=0.5D0                      ! "natural" lower b.c.
          do i=istr,iend
            drv(i,1)=1.5D0*(qsrc(i,j,2)-qsrc(i,j,1))
          enddo
        endif

        do k=2,Np-1,+1            !--> forward elimination
          CF(k)=1.D0/(4.D0-CF(k-1)) ; cff=CF(k)
          do i=istr,iend
            drv(i,k)=cff*( 3.D0*(qsrc(i,j,k+1)-qsrc(i,j,k-1))
     &                                           -drv(i,k-1) )
          enddo
        enddo

        cff=1.D0/(2.D0-CF(Np-1))  !<-- upper b.c.
        do i=istr,iend
          drv(i,Np)=cff*( 3.D0*(qsrc(i,j,Np)-qsrc(i,j,Np-1))
     &                                         -drv(i,Np-1) )
        enddo
        do k=Np-1,1,-1            !--> backsubstitution
          cff=CF(k)
          do i=istr,iend
            drv(i,k)=drv(i,k)-cff*drv(i,k+1)
          enddo
        enddo
        do k=1,N
          do i=istr,iend
            kp=int(kprnt(i,j,k))
            kp=max( 1, int(kprnt(i,j,k)) )
            p=max(-0.5D0, kprnt(i,j,k)-float(kp) )
            if (kp > Np-1) then                   ! extrapolate within
              qtr(i,j,k)=qsrc(i,j,Np)+p*drv(i,Np) ! upper half of upper
            elseif (kp > 0) then                  ! most grid box.
              q=1.D0-p ; pq=p*q
              qtr(i,j,k)=p*((p+2.D0*pq)*qsrc(i,j,kp+1) -pq*drv(i,kp+1))
     &                  +q*((q+2.D0*pq)*qsrc(i,j,kp  ) +pq*drv(i,kp  ))
c**         else
c**           qtr(i,j,k)=qsrc(i,j,1)
            endif
          enddo
        enddo
        if (lmsk == 1) then
          do i=istr,iend
            if (mask(i,j) == 0) then
              do k=1,N
                qtr(i,j,k)=FillValue
              enddo
            endif
          enddo
        endif
      enddo
      end
#endif


      subroutine r2r_set_depth_tile(istr,iend,jstr,jend, ncx,ncy,
     &                                             h, N,hc,Cs, z)

! Compute
!                                  hc*s(k) + Cs(k)*h(i,j)
!             z(i,j,k) = h(i,j) * ------------------------
!                                       hc + h(i,j)
!
! at vertical RHO- or W-points which is determined by looking at the
! first value of Cs(k): for for W-type -1 <= Cs <= 0  reaching both -1
! and 0 at the ends; for RHO-type it stays fully inside the interval.
! Place the outcome into private array.

      implicit none
      integer istr,iend,jstr,jend, ncx,ncy,N, i,j,k
      real(kind=8) h(ncx,ncy), hc, Cs(N), z(istr:iend,jstr:jend,N),
     &                                           hcds,bias,cf1,cf2
      if (Cs(1) > -0.999999D0) then
        hcds=hc/dble(N) ; bias=0.5D0  !<-- RHO-points
      else
        hcds=hc/dble(N-1) ; bias=0.D0  !<-- W-points
      endif
      do j=jstr,jend
        do i=istr,iend
          z(i,j,N)=h(i,j)/(hc+h(i,j)) !<-- temporarily
        enddo
        do k=1,N
          cf1=hcds*(dble(k-N)-bias) ; cf2=Cs(k)
          do i=istr,iend
            z(i,j,k)= z(i,j,N)*(cf1 + cf2*h(i,j))
          enddo
        enddo
      enddo
      end


      module comp_uvbar_vars
      integer, save :: alloc_zw_size=0
C$OMP THREADPRIVATE(alloc_zw_size)
      real(kind=8), allocatable, dimension(:) :: D,z_w,Hz
C$OMP THREADPRIVATE(D,z_w,Hz)
      end module comp_uvbar_vars


      subroutine compute_uvbar_thread(ncx,ncy, N, hc,Cs_w, h, u,v,
     &                                                  ubar,vbar)

! Compute barotropic velocities "ubar" and "vbar" from "u" and "v" --
! basically it is just vertical integration, except that the use of 3D
! arrays z_w Hz is avoided for memory reason: using vertical slice of
! "z_w" instead, resulting in somewhat awkward-looking code.

      use comp_uvbar_vars
      implicit none
      integer ncx,ncy, N
      real(kind=8) hc, Cs_w(0:N),  h(ncx,ncy)
      real(kind=4) u(2:ncx,ncy,N), ubar(2:ncx,ncy),
     &             v(ncx,2:ncy,N), vbar(ncx,2:ncy)
      integer ntrds,trd, nsub_x,nsub_y, tile,my_first,my_last,
     &        ntls,  isize,  istr,iend,jstr,jend, istrU,jstrV
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zw_size < (isize+2)*(N+1)) then
        alloc_zw_size=(isize+2)*(N+1)
        if (allocated(z_w)) deallocate(Hz,z_w,D)
        allocate(D(isize+2), z_w(alloc_zw_size), Hz(2*alloc_zw_size))
C$OMP CRITICAL(comp_uvbar_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)')  'compute_uvbar_thread: ',
     &   'allocated', (2*(isize+1)+6*alloc_zw_size)/float(262144),
     &                        'MB private workspace, trd =', trd
C$OMP END CRITICAL(comp_uvbar_crgn)
      endif

      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     &                                    istr,iend,jstr,jend)
        istrU=max(istr,2) ; jstrV=max(jstr,2)
        call compute_uvbar_tile(istr,iend,jstr,jend, istrU,jstrV,
     &         ncx,ncy, N, hc,Cs_w, h, u,v, ubar,vbar, D,z_w,Hz)
      enddo
      end

! Note that istrU,jstrV are computed above by the driver rather than
! inside the working routine below because the first dimension of "z_w"
! must be set from "istrU-1" rather than based on "istr" (say "istr-1")
! because private array dimensions inside "r2r_set_depth_tile" are tied
! to its loop ranges, leaving no other choice but to have horizontal
! dimensions of "z_w" be the same in both "compute_uvbar_tile" and
! "r2r_set_depth_tile".

      subroutine compute_uvbar_tile(istr,iend,jstr,jend, istrU,jstrV,
     &              ncx,ncy, N, hc,Cs_w, h, u,v, ubar,vbar, D,z_w,Hz)
      implicit none
      integer istr,iend,jstr,jend, istrU,jstrV, ncx,ncy, N, i,j,k,j1
      real(kind=8) hc, Cs_w(0:N),  h(ncx,ncy),  cff,  D(istrU-1:iend),
     &                   z_w(istrU-1:iend,0:N),  Hz(istrU-1:iend,2,N)
      real(kind=4) u(2:ncx,ncy,N), ubar(2:ncx,ncy),
     &             v(ncx,2:ncy,N), vbar(ncx,2:ncy)

      j1=2
      do j=jstrV-1,jend,+1 !--> recursive because of Hz.
        call r2r_set_depth_tile(istrU-1,iend, j,j, ncx,ncy,
     &                                h, N+1,hc,Cs_w, z_w)
        j1=3-j1 !<-- rotate index
        do k=1,N
          do i=istrU-1,iend
            Hz(i,j1,k)=z_w(i,k)-z_w(i,k-1)
          enddo
        enddo
        if (j > jstr-1) then
          do i=istrU,iend
            ubar(i,j)=0. ; D(i)=0.D0
          enddo
          do k=N,1,-1
            do i=istrU,iend
              cff=Hz(i,j1,k)+Hz(i-1,j1,k)
              D(i)=D(i)+cff
              ubar(i,j)=ubar(i,j)+cff*u(i,j,k)
            enddo
          enddo
          do i=istrU,iend
            ubar(i,j)=ubar(i,j)/D(i)
          enddo
        endif
        if (j > jstrV-1) then
          do i=istr,iend
            vbar(i,j)=0. ; D(i)=0.D0
          enddo
          do k=N,1,-1                        ! Here Hz(i,3-j1,k) is
            do i=istr,iend                   ! the j-slice computed
              cff=Hz(i,j1,k)+Hz(i,3-j1,k)    ! during the previous j.
              D(i)=D(i)+cff
              vbar(i,j)=vbar(i,j)+cff*v(i,j,k)
            enddo
          enddo
          do i=istr,iend
            vbar(i,j)=vbar(i,j)/D(i)
          enddo
        endif
      enddo !<-- j
      end
