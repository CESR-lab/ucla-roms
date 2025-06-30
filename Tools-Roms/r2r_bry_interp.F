! Same functionality as "r2r_interp_init.F", except that now it is
! just for a single 1D line rather than a 2D array;  all algorithms
! are the same, with the exception of not having alternating ic-index
! sweeps consecutive jcs.

      subroutine bry_init_interp(nx,ny, x,y, ncx, xc,yc, ip,jp, xi,eta)
      implicit none
      integer nx,ny, ncx
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
C$OMP PARALLEL SHARED(nx,ny,x,y, ncx, xc,yc, ip,jp, xi,eta)
      call bry_init_interp_thread(nx,ny,x,y, ncx,xc,yc,
     &                                       ip,jp, xi,eta)
C$OMP END PARALLEL
      call check_search_line(nx,ny, x,y, ncx, xc,yc, ip,jp)
      call check_offsts_line(nx,ny, x,y, ncx, xc,yc, ip,jp,
     &                                              xi,eta)
      end

      subroutine bry_init_interp_thread(nx,ny,x,y, ncx, xc,yc, ip,jp,
     &                                                         xi,eta)
      implicit none
      integer nx,ny, ncx,    icmin,icmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$    icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      call bry_search_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc, ip,jp)
      call compute_offsts_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc,
     &                                                  ip,jp, xi,eta)
      end


      subroutine r2r_bry_search(nx,ny,x,y, ncx,xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc
      integer(kind=4), dimension(ncx) :: ip,jp
C$OMP PARALLEL SHARED(nx,ny,x,y, ncx,xc,yc, ip,jp)
      call r2r_bry_search_thread(nx,ny,x,y, ncx,xc,yc, ip,jp)
C$OMP END PARALLEL
      call check_search_line(    nx,ny,x,y, ncx,xc,yc, ip,jp)
      end

      subroutine r2r_bry_search_thread(nx,ny,x,y, ncx,xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx,  icmin,icmax             ! Search for parent
      real(kind=8), dimension(nx,ny) :: x,y        ! grid indices only,
      real(kind=8), dimension(ncx) :: xc,yc        ! do not compute
      integer(kind=4), dimension(ncx) :: ip,jp     ! xi,eta offsets.
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$    icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      call bry_search_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc, ip,jp)
      end


      subroutine bry_search_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc,
     &                                                         ip,jp)
      implicit none
      integer icmin,icmax, nx,ny, ncx,    ic,i,j,inew,jnew
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx)  ::  xc,yc
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) dx4,dy4, dx3,dy3, r41,r34,
     &             dx1,dy1, dx2,dy2, r12,r23
#ifdef VERBOSE
      write(*,*) 'enter search_line, icmin =',icmin, ' icmax =', icmax
#endif
      i=nx/2 ; j=ny/2 ; ic=icmin                       !   3 <---- 2
      dx1=x(i+1,j  )-xc(ic) ; dy1=y(i+1,j  )-yc(ic)    !   !       ^
      dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic)    !   !       !
      dx3=x(i  ,j+1)-xc(ic) ; dy3=y(i  ,j+1)-yc(ic)    !   v       !
      dx4=x(i  ,j  )-xc(ic) ; dy4=y(i  ,j  )-yc(ic)    !   4 ----> 1

      r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
      r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      do while(.true.)
        do while( r12 >= 0.D0 .and. r23 >= 0.D0
     &      .and. r34 >= 0.D0 .and. r41 >= 0.D0 )
#ifdef VERBOSE
          write(*,*) 'ic =', ic
#endif
          ip(ic)=i ; jp(ic)=j    !<-- record bounding indices
          ic=ic+1                !--> proceed to the next target point
          if (ic > icmax) return

          dx1=x(i+1,j  )-xc(ic) ; dy1=y(i+1,j  )-yc(ic)
          dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic)
          dx3=x(i  ,j+1)-xc(ic) ; dy3=y(i  ,j+1)-yc(ic)
          dx4=x(i  ,j  )-xc(ic) ; dy4=y(i  ,j  )-yc(ic)

          r12=dx1*dy2-dx2*dy1   ; r23=dx2*dy3-dx3*dy2
          r34=dx3*dy4-dx4*dy3   ; r41=dx4*dy1-dx1*dy4
        enddo

        inew=i ; jnew=j
        if (r12 < 0.D0) inew=min(inew+1, nx-1)
        if (r23 < 0.D0) jnew=min(jnew+1, ny-1)
        if (r34 < 0.D0) inew=max(inew-1,    1)
        if (r41 < 0.D0) jnew=max(jnew-1,    1)

        if (inew == i .and. jnew == j) then
#ifdef VERBOSE
          write(*,*) 'ic =', ic
#endif
          ip(ic)=-1 ; jp(ic)=-1         !<-- cannot be bounded
          ic=ic+1
          if (ic > icmax) return
        else
          i=inew ; j=jnew      !--> accept move of source location.
        endif

        dx1=x(i+1,j  )-xc(ic) ; dy1=y(i+1,j  )-yc(ic)
        dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic)
        dx3=x(i  ,j+1)-xc(ic) ; dy3=y(i  ,j+1)-yc(ic)
        dx4=x(i  ,j  )-xc(ic) ; dy4=y(i  ,j  )-yc(ic)

        r12=dx1*dy2-dx2*dy1   ; r23=dx2*dy3-dx3*dy2
        r34=dx3*dy4-dx4*dy3   ; r41=dx4*dy1-dx1*dy4

      enddo !<-- while(.true.)
      end

      subroutine compute_offsts_line(icmin,icmax, nx,ny, x,y, ncx,
     &                                        xc,yc,ip,jp, xi,eta)
      implicit none
      integer icmin,icmax,  nx,ny, ncx,  ic, i,j, iter
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) a11,a12,a21,a22, det, dX,dY, p,p1,q,q1

      do ic=icmin,icmax
        i=ip(ic) ; j=jp(ic)
        if (i > 0 .and. j > 0) then
          a11=0.5D0*(x(i+1,j+1)-x(i,j+1) +x(i+1,j)-x(i,j))
          a12=0.5D0*(x(i+1,j+1)-x(i+1,j) +x(i,j+1)-x(i,j))
          a21=0.5D0*(y(i+1,j+1)-y(i,j+1) +y(i+1,j)-y(i,j))
          a22=0.5D0*(y(i+1,j+1)-y(i+1,j) +y(i,j+1)-y(i,j))

          dX=xc(ic)-0.25D0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))
          dY=yc(ic)-0.25D0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))

          det=1.D0/(a11*a22-a12*a21)
          xi(ic) =0.5D0 + det*(a22*dX-a12*dY)
          eta(ic)=0.5D0 + det*(a11*dY-a21*dX)
        else
          xi(ic)=-1.D0 ; eta(ic)=-1.D0  !<-- special values
        endif
      enddo
      do iter=1,10
        do ic=icmin,icmax
          i=ip(ic) ; j=jp(ic)
          if (i > 0 .and. j > 0) then
            p=xi(ic)  ; p1=1.D0-p
            q=eta(ic) ; q1=1.D0-q

            a11=q*(x(i+1,j+1)-x(i,j+1)) +q1*(x(i+1,j)-x(i,j))
            a12=p*(x(i+1,j+1)-x(i+1,j)) +p1*(x(i,j+1)-x(i,j))
            a21=q*(y(i+1,j+1)-y(i,j+1)) +q1*(y(i+1,j)-y(i,j))
            a22=p*(y(i+1,j+1)-y(i+1,j)) +p1*(y(i,j+1)-y(i,j))

            dX=xc(ic) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     &                 -p*q1*x(i+1,j)  -p1*q1*x(i,j)
            dY=yc(ic) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     &                 -p*q1*y(i+1,j)  -p1*q1*y(i,j)

            det=1.D0/(a11*a22-a12*a21)

            xi(ic) =p + det*(a22*dX-a12*dY)
            eta(ic)=q + det*(a11*dY-a21*dX)
          endif
        enddo
      enddo
      end

      subroutine check_search_line(nx,ny,x,y, ncx, xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx,       ic, i,j
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) dx4,dy4, dx3,dy3, r41,r34,
     &             dx1,dy1, dx2,dy2, r12,r23

      write(*,'(2x,A)',advance='no') 'check_search_line ...'
      do ic=1,ncx
        i=ip(ic) ; j=jp(ic)
        if (i > 0 .and. j > 0) then
          dx1=x(i+1,j  )-xc(ic) ; dy1=y(i+1,j  )-yc(ic)
          dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic)
          dx3=x(i  ,j+1)-xc(ic) ; dy3=y(i  ,j+1)-yc(ic)
          dx4=x(i  ,j  )-xc(ic) ; dy4=y(i  ,j  )-yc(ic)

          r12=dx1*dy2-dx2*dy1   ; r23=dx2*dy3-dx3*dy2
          r34=dx3*dy4-dx4*dy3   ; r41=dx4*dy1-dx1*dy4

          if (r12<0.D0 .or. r23<0.D0 .or. r34<0.D0 .or. r41<0.D0)
     &     write(*,*) '### ERROR: Search algorithm failure at ',
     &                                     'ic =', ic
        endif
      enddo
      write(*,'(2x,A)') '...done'
      end

      subroutine check_offsts_line(nx,ny, x,y, ncx, xc,yc, ip,jp,
     &                                                     xi,eta)
      implicit none
      integer nx,ny, ncx,  ic, i,j
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) dX,dY, p,p1,q,q1, errX,errY

      write(*,'(2x,A)',advance='no') 'checking offsets...'
      errX=0.D0 ; errY=0.D0
      do ic=1,ncx                   ! Notice simple semantic rule
        i=ip(ic) ; j=jp(ic)         ! here: one and only one "1" is
        if (i > 0 .and. j > 0) then ! always present:  either as +1
          p=xi(ic)  ; p1=1.D0-p     ! in index, or as p1,q1 in the
          q=eta(ic) ; q1=1.D0-q     ! corresponding coefficient.

          dX=xc(ic) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     &               -p*q1*x(i+1,j)  -p1*q1*x(i,j)
          dY=yc(ic) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     &               -p*q1*y(i+1,j)  -p1*q1*y(i,j)

          if (abs(dX) > errX) errX=abs(dX)
          if (abs(dY) > errY) errY=abs(dY)
        endif
      enddo
      write(*,'(2x,A,2ES22.15)') 'max errors =', errX,errY
      end



      subroutine bry_interp(nx,ny,N, src, ncx, ip,jp,xi,eta, msk, targ)
      implicit none
      integer nx,ny,N, ncx
      integer(kind=4) ip(ncx),jp(ncx)
      real(kind=8) xi(ncx), eta(ncx)
      integer(kind=2) msk(ncx)
      real(kind=4) src(nx,ny,N), targ(ncx,N)
C$OMP PARALLEL SHARED(nx,ny,N, src, ncx, ip,jp,xi,eta, msk, targ)
      call bry_interp_thread(nx,ny,N, src, ncx, ip,jp,xi,eta,msk, targ)
C$OMP END PARALLEL
      end

      subroutine bry_interp_thread(nx,ny,N, src, ncx, ip,jp,xi,eta,
     &                                                     msk,targ)
      implicit none
      integer nx,ny,N, ncx
      integer(kind=4) ip(ncx), jp(ncx)
      real(kind=8) xi(ncx), eta(ncx)
      integer(kind=2) msk(ncx)
      real(kind=4) src(nx,ny,N), targ(ncx,N)
      integer icmin,icmax,isize, istr,iend,tile
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$    icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call bry_interp_line( istr,iend, nx,ny,N, src, ncx, ip,jp,
     &                                         xi,eta, msk, targ)
      enddo
      end

/*
      subroutine bry_interp_line(istr,iend, nx,ny,N, src, ncx, ip,jp,
     &                                             xi,eta, msk, targ)
      implicit none
      integer istr,iend, nx,ny,N, ncx, ic, i,j,k
      integer(kind=4) ip(ncx), jp(ncx)
      real(kind=8) xi(ncx), eta(ncx), p,p1,q,q1
      integer(kind=2) msk(ncx)
      real(kind=4) src(nx,ny,N), targ(ncx,N)

      do ic=istr,iend
        if (ip(ic) > 0) then
          i=ip(ic) ; p=xi(ic)  ; p1=1.D0-p
          j=jp(ic) ; q=eta(ic) ; q1=1.D0-q
          do k=1,N
            targ(ic,k)=p*q*src(i+1,j+1,k) +p1*q*src(i,j+1,k)
     &                 +p*q1*src(i+1,j,k)  +p1*q1*src(i,j,k)
          enddo
        else
          do k=1,N
            targ(ic,k)=0.
          enddo
        endif
      enddo
      end
*/




      subroutine bry_interp_line(istr,iend, nx,ny,N, src, ncx, ip,jp,
     &                                             xi,eta, msk, targ)
      implicit none
      integer istr,iend, nx,ny,N, ncx, ic, i,j,k
      integer(kind=4) ip(ncx), jp(ncx)
      real(kind=8) xi(ncx), eta(ncx)
      integer(kind=2) msk(ncx)
      real(kind=4) src(nx,ny,N), targ(ncx,N)
      real(kind=8), parameter :: TwoThird=2.D0/3.D0,
     &                         FourNineth=4.D0/9.D0
      real(kind=8) px,qx,pqx, HxL,HxR,GxL,GxR, FxLL,FxRL,FxLR,FxRR,
     &             py,qy,pqy, HyL,HyR,GyL,GyR, FyLL,FyRL,FyLR,FyRR,
     &                                     FxyLL,FxyRL,FxyLR,FxyRR

      do ic=istr,iend
        i=ip(ic) ; j=jp(ic) ;
        if ( msk(ic) == 1 .and. 2 < i .and. i < nx-2 .and.
     &                          2 < j .and. j < ny-2 ) then
          px=xi(ic) ; qx=1.D0-px; pqx=px*qx
          GxR=-px*pqx ; HxR=px*(px+2.D0*pqx)
          GxL= qx*pqx ; HxL=qx*(qx+2.D0*pqx)

          py=eta(ic); qy=1.D0-py; pqy=py*qy
          GyR=-py*pqy ; HyR=py*(py+2.D0*pqy)
          GyL= qy*pqy ; HyL=qy*(qy+2.D0*pqy)

          do k=1,N
            FxLL=TwoThird*( src(i+1,j,k)-src(i-1,j,k)
     &            -0.125D0*(src(i+2,j,k)-src(i-2,j,k)))
            FxRL=TwoThird*( src(i+2,j,k)-src(i  ,j,k)
     &            -0.125D0*(src(i+3,j,k)-src(i-1,j,k)))
            FxLR=TwoThird*( src(i+1,j+1,k)-src(i-1,j+1,k)
     &            -0.125D0*(src(i+2,j+1,k)-src(i-2,j+1,k)))
            FxRR=TwoThird*( src(i+2,j+1,k)-src(i  ,j+1,k)
     &            -0.125D0*(src(i+3,j+1,k)-src(i-1,j+1,k)))

            FyLL=TwoThird*( src(i,j+1,k)-src(i,j-1,k)
     &            -0.125D0*(src(i,j+2,k)-src(i,j-2,k)))
            FyRL=TwoThird*( src(i+1,j+1,k)-src(i+1,j-1,k)
     &            -0.125D0*(src(i+1,j+2,k)-src(i+1,j-2,k)))
            FyLR=TwoThird*( src(i,j+2,k)-src(i,j  ,k)
     &            -0.125D0*(src(i,j+3,k)-src(i,j-1,k)))
            FyRR=TwoThird*( src(i+1,j+2,k)-src(i+1,j  ,k)
     &            -0.125D0*(src(i+1,j+3,k)-src(i+1,j-1,k)))

            FxyLL=FourNineth*( -src(i-1,j+1,k) +src(i+1,j+1,k)
     &                         +src(i-1,j-1,k) -src(i+1,j-1,k)
     &            +0.125D0*(  src(i-1,j+2,k)     -src(i+1,j+2,k)
     &                      +src(i-2,j+1,k)       -src(i+2,j+1,k)
     &                      -src(i-2,j-1,k)       +src(i+2,j-1,k)
     &                       -src(i-1,j-2,k)     +src(i+1,j-2,k)
     &            +0.125D0*(   -src(i-2,j+2,k) +src(i+2,j+2,k)
     &                         +src(i-2,j-2,k) -src(i+2,j-2,k)
     &                                                   )))
            FxyRL=FourNineth*( -src(i  ,j+1,k) +src(i+2,j+1,k)
     &                         +src(i  ,j-1,k) -src(i+2,j-1,k)
     &            +0.125D0*(  src(i  ,j+2,k)     -src(i+2,j+2,k)
     &                      +src(i-1,j+1,k)       -src(i+3,j+1,k)
     &                      -src(i-1,j-1,k)       +src(i+3,j-1,k)
     &                       -src(i  ,j-2,k)     +src(i+2,j-2,k)
     &            +0.125D0*(   -src(i-1,j+2,k) +src(i+3,j+2,k)
     &                         +src(i-1,j-2,k) -src(i+3,j-2,k)
     &                                                   )))
            FxyLR=FourNineth*( -src(i-1,j+2,k) +src(i+1,j+2,k)
     &                         +src(i-1,j  ,k) -src(i+1,j  ,k)
     &            +0.125D0*(  src(i-1,j+3,k)     -src(i+1,j+3,k)
     &                      +src(i-2,j+2,k)       -src(i+2,j+2,k)
     &                      -src(i-2,j  ,k)       +src(i+2,j  ,k)
     &                       -src(i-1,j-1,k)     +src(i+1,j-1,k)
     &            +0.125D0*(   -src(i-2,j+3,k) +src(i+2,j+3,k)
     &                         +src(i-2,j-1,k) -src(i+2,j-1,k)
     &                                                    )))
            FxyRR=FourNineth*( -src(i  ,j+2,k) +src(i+2,j+2,k)
     &                         +src(i  ,j  ,k) -src(i+2,j  ,k)
     &            +0.125D0*(  src(i  ,j+3,k)     -src(i+2,j+3,k)
     &                      +src(i-1,j+2,k)       -src(i+3,j+2,k)
     &                      -src(i-1,j  ,k)       +src(i+3,j  ,k)
     &                       -src(i  ,j-1,k)     +src(i+2,j-1,k)
     &            +0.125D0*(   -src(i-1,j+3,k) +src(i+3,j+3,k)
     &                         +src(i-1,j-1,k) -src(i+3,j-1,k)
     &                                                    )))

            targ(ic,k) = HxL*HyL*src(i,j  ,k) + HxR*HyL*src(i+1,j  ,k)
     &                 + HxL*HyR*src(i,j+1,k) + HxR*HyR*src(i+1,j+1,k)

     &                         + GxL*HyL*FxLL + GxR*HyL*FxRL
     &                         + GxL*HyR*FxLR + GxR*HyR*FxRR

     &                         + HxL*GyL*FyLL + HxR*GyL*FyRL
     &                         + HxL*GyR*FyLR + HxR*GyR*FyRR

     &                        + GxL*GyL*FxyLL + GxR*GyL*FxyRL
     &                        + GxL*GyR*FxyLR + GxR*GyR*FxyRR
          enddo
        else
          do k=1,N
            targ(ic,k)=0.  !<-- no special value masking here yet.
          enddo
        endif
      enddo
      end





! An alternative driver for "r2r_init_vrtint_tile" (r2r_vert_interp.F)
! designed for 1D line as needed by r2r_bry.   The differences are:
! (i) the second horizontal dimension no longer exists; (ii) instead of
! using 2D tiling, it cuts "ncx" into a set of chunks according to the
! number of theads; (iii) external layer with parallel region is added;
! (iv) there is no option VERBOSE mode.


      subroutine bry_init_vertinterp(ncx, hprnt, Np,hcp,Csp_r,
     &                                    h, N,hc,Cs_r, kprnt)
      implicit none
      integer ncx, Np,N
      real(kind=8) hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     &                                              kprnt(ncx,N)
C$OMP PARALLEL SHARED(ncx,hprnt, Np,hcp,Csp_r,h,N,hc,Cs_r,kprnt)
      call bry_init_vinterp_thread(ncx, hprnt, Np,hcp,Csp_r,
     &                                      h, N,hc,Cs_r, kprnt)
C$OMP END PARALLEL
      end

      subroutine bry_init_vinterp_thread(ncx, hprnt, Np,hcp,Csp_r,
     &                                        h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer ncx, Np,N
      real(kind=8) hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     &                                              kprnt(ncx,N)
      integer icmin,icmax,isize, istr,iend,tile
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$    icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      if (alloc_zc_size < isize*(N+1) .or.
     &    allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) deallocate(zp_r,drv,zc)
        allocate( zc(alloc_zc_size), zp_r(allc_zpr_size),
     &                                drv(allc_zpr_size) )
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)')  'bry_init_vinterp_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     &                           'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
      endif
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call r2r_init_vrtint_tile( istr,iend, 1, 1, ncx,1,
     &                           hprnt,Np,hcp,Csp_r, zp_r,drv,
     &                               h,N, hc, Cs_r, zc, kprnt)
      enddo
      end

! Same, but for checking routine.

      subroutine bry_check_init_vertinterp(ncx, hprnt, Np,hcp,Csp_r,
     &                                         h, N,hc,Cs_r, kprnt)
      implicit none
      integer ncx, Np,N
      real(kind=8) hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     &                                              kprnt(ncx,N)
C$OMP PARALLEL SHARED(ncx,hprnt, Np,hcp,Csp_r,h,N,hc,Cs_r,kprnt)
      call bry_check_vrtint_thread(ncx, hprnt, Np,hcp,Csp_r,
     &                                      h, N,hc,Cs_r, kprnt)
C$OMP END PARALLEL
      end

      subroutine bry_check_vrtint_thread(ncx, hprnt, Np,hcp,Csp_r,
     &                                        h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer ncx, Np,N
      real(kind=8) hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     &                                    kprnt(ncx,N), my_error
      integer icmin,icmax,isize, istr,iend,tile
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$    icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      if (alloc_zc_size < isize*(N+1) .or.
     &    allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) deallocate(zp_r,drv,zc)
        allocate( zc(alloc_zc_size), zp_r(allc_zpr_size),
     &                                drv(allc_zpr_size) )
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)')  'bry_check_vrtint_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     &                           'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
      endif
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call r2r_check_vrtint_tile( istr,iend, 1,1, ncx,1,
     &                          hprnt, Np,hcp,Csp_r,zp_r,drv,
     &                       h, N,hc,Cs_r,zc, kprnt, my_error)
      enddo
C$OMP CRITICAL(r2r_vert_crgn)
      if (trd_count == 0) vert_int_error=0.D0
      trd_count=trd_count+1
      vert_int_error=max(vert_int_error, my_error)
C$    if (trd_count == numthreads) then
        trd_count=0
        write(*,*) '          maximum vert_int_error =', vert_int_error
C$    endif
C$OMP END CRITICAL(r2r_vert_crgn)
      end



#define SPLINE


      subroutine bry_vertinterp(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
      implicit none
      integer ncx, btm_bc, Np,N
      real(kind=4) qsrc(ncx,Np), qtr(ncx,N)
      real(kind=8) kprnt(ncx,N)
C$OMP PARALLEL SHARED(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
      call bry_vinterp_thread(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
C$OMP END PARALLEL
      end

      subroutine bry_vinterp_thread(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
      use r2r_vertint_vars !<-- needed for "drv"
      implicit none
      integer ncx, btm_bc, Np,N
      real(kind=4) qsrc(ncx,Np), qtr(ncx,N)
      real(kind=8) kprnt(ncx,N)
      integer, parameter :: lmsk=0   !<-- suppress land masking
      integer(kind=2) mask(8)  !<-- for compatibility; not used
      integer icmin,icmax,isize, istr,iend,tile
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$    icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1,icmax)
#ifdef SPLINE
        call r2r_vsplnint_tile( istr,iend, 1,1, ncx,1,
     &              lmsk,mask, btm_bc, Np,qsrc,drv, N,kprnt,qtr)
#else
        call r2r_vertint_tile(  iistr,iend, 1,1, ncx,1,
     &              lmsk,mask,         Np,qsrc,     N,kprnt,qtr)
#endif
      enddo
      end



      subroutine bry_vert_average(ncx, h, N,hc,Cs_w, qsrc,qbar)
      implicit none
      integer ncx,N
      real(kind=8) h(ncx), hc, Cs_w(0:N)
      real(kind=4) qsrc(ncx,N),qbar(ncx)
C$OMP PARALLEL SHARED(ncx,N, h, hc,Cs_w, qsrc,qbar)
      call bry_vert_average_thread(ncx, h, N, hc,Cs_w, qsrc,qbar)
C$OMP END PARALLEL
      end

      subroutine bry_vert_average_thread(ncx, h, N,hc,Cs_w, qsrc,qbar)
      use r2r_vertint_vars
      implicit none
      integer ncx,N
      real(kind=8) h(ncx), hc, Cs_w(0:N)
      real(kind=4) qsrc(ncx,N), qbar(ncx)
      integer icmin,icmax,isize, istr,iend,tile
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$    icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      if (alloc_zc_size < isize*(N+1)) then
        alloc_zc_size=isize*(N+1)
        if (allocated(zc)) deallocate(zc)
        allocate(zc(alloc_zc_size))
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)')  'bry_vert_average_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     &                           'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
      endif
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call bry_vert_average_line(istr,iend, ncx,h, N,hc,Cs_w,
     &                                           qsrc,qbar, zc)
      enddo
      end

      subroutine bry_vert_average_line(istr,iend, ncx, h, N,hc,Cs_w,
     &                                               qsrc,qbar, z_w)
      implicit none
      integer istr,iend, ncx,N, i,k
      real(kind=8) h(ncx), hc, Cs_w(0:N), z_w(istr:iend,0:N)
      real(kind=4) qsrc(ncx,N),qbar(ncx)

      call r2r_set_depth_tile(istr,iend,1,1, ncx,1,h, N+1,hc,Cs_w,z_w)
      do i=istr,iend
        qbar(i)=qsrc(i,N)*(z_w(i,N)-z_w(i,N-1))
      enddo
      do k=N-1,1,-1
        do i=istr,iend
          qbar(i)=qbar(i)+qsrc(i,k)*(z_w(i,k)-z_w(i,k-1))
        enddo
      enddo
      do i=istr,iend
        qbar(i)=qbar(i)/(z_w(i,N)-z_w(i,0))
      enddo
      end
