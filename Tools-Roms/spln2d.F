#ifndef KIND
# define KIND 4
# define spln2d_TYPE             spln2d_real
# define spln2d_interp_TYPE      spln2d_interp_real
# define spln2d_interp_TYPE_tile spln2d_interp_real_tile
c--# define VERBOSE
#elif KIND == 8
# define spln2d_TYPE             spln2d_double
# define spln2d_interp_TYPE      spln2d_interp_double
# define spln2d_interp_TYPE_tile spln2d_interp_double_tile
#endif


! The following routine takes a two-dimensional array "src" and
! computes its derivatives
!
!        d src           d src           d^2 src           d^2 src
! srX = -------   srY = -------   sXY = ---------   sYX = ---------
!          dx             dy              dx dy             dy dx
!
! to construct  bi-cubic spline interpolation.  All the derivatives
! are defined at the same grid locations as the field itself and are
! obtained from the condition of continuity of the second derivative
! in the same direction.  The spline algorithm is one-dimensional
! Gaussian elimination procedure applied individually in each
! direction, so the mixed derivatives are computed as derivative of
! previously computed derivative in transversal direction.
! It is expected that sXY == sYX, however this routine computes both
! of them (just a matter of which direction goes first and which
! follows) and compares the outcome.
!
! All the derivatives are defined as "with respect to array index",
! i.e., dx=dy=1 in both directions.

! "Natural" boundary conditions are assumed on the sides (the second
! derivative vanishes at the ends).

! Array "src" which can be be either single or double precision,
! matching TYPE in the name of the subroutine, however all computed
! derivatives srX,srY,sXY,sYX are stored in double precision arrays
! regardless of type of "src".


      subroutine spln2d_TYPE(nx,ny, src, srX,srY,sXY,sYX)
      implicit none
      integer nx,ny, i,j
      real(kind=KIND), dimension(nx,ny) :: src
      real(kind=8), dimension(nx,ny) :: srX,srY,sXY,sYX
      real(kind=8) cff, error
      real(kind=8), dimension(:), allocatable :: CF,FC

#ifdef VERBOSE
      write(*,*) 'enter spln2d_TYPE...'
#endif
      allocate(CF(nx), FC(ny), stat=i)
      if (i /= 0) then
        write(*,'(/1x,2A/)') '### ERROR: spln2d_TYPE :: ',
     &                       'Memory allocation failure.'
        stop
      endif

      CF(1)=0.5D0
      do j=1,ny
        srX(1,j)=1.5D0*(src(2,j)-src(1,j))
      enddo
      do i=2,nx-1,+1
        CF(i)=1.D0/(4.D0-CF(i-1))
        cff=CF(i)
        do j=1,ny
          srX(i,j)=cff*(3.D0*(src(i+1,j)-src(i-1,j)) -srX(i-1,j))
        enddo
      enddo
      cff=1.D0/(2.D0-CF(nx-1))
      do j=1,ny
        srX(nx,j)=cff*(3.D0*(src(nx,j)-src(nx-1,j))-srX(nx-1,j))
      enddo
      do i=nx-1,1,-1
        cff=CF(i)
        do j=1,ny
          srX(i,j)=srX(i,j)-cff*srX(i+1,j)
        enddo
      enddo

      FC(1)=0.5D0
      do i=1,nx
        srY(i,1)=1.5D0*(src(i,2)-src(i,1))
        sXY(i,1)=1.5D0*(srX(i,2)-srX(i,1))
      enddo
      do j=2,ny-1,+1
        FC(j)=1.D0/(4.D0-FC(j-1))
        cff=FC(j)
        do i=1,nx
          srY(i,j)=cff*( 3.D0*(src(i,j+1)-src(i,j-1)) -srY(i,j-1))
          sXY(i,j)=cff*( 3.D0*(srX(i,j+1)-srX(i,j-1)) -sXY(i,j-1))
        enddo
      enddo
      cff=1.D0/(2.D0-FC(ny-1))
      do i=1,nx
        srY(i,ny)=cff*( 3.D0*(src(i,ny)-src(i,ny-1)) -srY(i,ny-1))
        sXY(i,ny)=cff*( 3.D0*(srX(i,ny)-srX(i,ny-1)) -sXY(i,ny-1))
      enddo
      do j=ny-1,1,-1
        cff=FC(j)
        do i=1,nx
          srY(i,j)=srY(i,j)-cff*srY(i,j+1)
          sXY(i,j)=sXY(i,j)-cff*sXY(i,j+1)
        enddo
      enddo

      CF(1)=0.5D0
      do j=1,ny
        sYX(1,j)=1.5D0*(srY(2,j)-srY(1,j))
      enddo
      do i=2,nx-1,+1
        CF(i)=1.D0/(4.D0-CF(i-1))
        cff=CF(i)
        do j=1,ny
          sYX(i,j)=cff*(3.D0*(srY(i+1,j)-srY(i-1,j)) -sYX(i-1,j))
        enddo
      enddo
      cff=1.D0/(2.D0-CF(nx-1))
      do j=1,ny
        sYX(nx,j)=cff*(3.D0*(srY(nx,j)-srY(nx-1,j))-sYX(nx-1,j))
      enddo
      do i=nx-1,1,-1
        cff=CF(i)
        do j=1,ny
          sYX(i,j)=sYX(i,j)-cff*sYX(i+1,j)
        enddo
      enddo

      error=0.D0
      do j=1,ny
        do i=1,nx
# ifdef VERBOSE
          cff=abs((sXY(i,j)-sYX(i,j))/(sXY(i,j)+sYX(i,j)))
          if (cff > error) error=cff
# endif
          sXY(1,j)=0.5D0*(sXY(i,j)+sYX(i,j))
        enddo
      enddo
# ifdef VERBOSE
      write(*,'(1x,2A/25x,A,ES22.15)') '...leaving spln2d_TYPE, ',
     &        'maximum normalized error max{|dXY-dYX|/|dXY+dYX|}',
     &        'in computing mixed derivative', error
# endif
      end



      subroutine spln2d_interp_TYPE(nx,ny, src, srX, srY, sXY,
     &                           ncx,ncy, ip,jp, xi,eta, targ)
      implicit none
      integer nx,ny, ncx,ncy
      real(kind=KIND), dimension(nx,ny) :: src
      real(kind=8),    dimension(nx,ny) :: srX,srY,sXY
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      real(kind=8),    dimension(ncx,ncy) :: xi,eta
      real(kind=KIND), dimension(ncx,ncy) :: targ
C$    integer numthreads, trd, omp_get_num_threads, omp_get_thread_num
      integer icmin,icmax,jcmin,jcmax, chunk_size
      icmin=1 ; icmax=ncx ; jcmin=1 ; jcmax=ncy ; chunk_size=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncy+numthreads-1)/numthreads
C$    jcmin=1+trd*chunk_size ; jcmax=min(jcmin+chunk_size-1,ncy)

      call spln2d_interp_TYPE_tile(nx,ny, src, srX, srY, sXY,
     &                           ncx,ncy, ip,jp, xi,eta, targ,
     &                               icmin,icmax,jcmin,jcmax)
      end

      subroutine spln2d_interp_TYPE_tile(nx,ny, src, srX, srY, sXY,
     &                                 ncx,ncy, ip,jp, xi,eta, targ,
     &                                     icmin,icmax,jcmin,jcmax)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=KIND), dimension(nx,ny) :: src
      real(kind=8),    dimension(nx,ny) :: srX,srY,sXY
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      real(kind=8),    dimension(ncx,ncy) :: xi,eta
      real(kind=KIND), dimension(ncx,ncy) :: targ
!-->  internal variables
      integer ic,jc, iL,iR, jL,jR
      real(kind=8) px,qx,pqx, HxR,HxL,GxR,GxL,
     &             py,qy,pqy, HyR,HyL,GyR,GyL

# ifdef VERBOSE
      write(*,*) 'enter spln2d_interp_TYPE_tile...'
# endif

      do jc=jcmin,jcmax
        do ic=icmin,icmax
          if (ip(ic,jc) > 0 .and. jp(ic,jc) > 0) then
            iL=ip(ic,jc)  ;  px=xi(ic,jc)
            iR=iL+1       ;  qx=1.D0-px ; pqx=px*qx
            GxR=-px*pqx   ;  HxR=px*(px+2.D0*pqx)
            GxL= qx*pqx   ;  HxL=qx*(qx+2.D0*pqx)

            jL=jp(ic,jc)  ;  py=eta(ic,jc)
            jR=jL+1       ;  qy=1.D0-py ; pqy=py*qy
            GyR=-py*pqy   ;  HyR=py*(py+2.D0*pqy)
            GyL= qy*pqy   ;  HyL=qy*(qy+2.D0*pqy)

            targ(ic,jc) = HxL*HyL*src(iL,jL) + HxR*HyL*src(iR,jL)
     &                  + HxL*HyR*src(iL,jR) + HxR*HyR*src(iR,jR)

     &                  + GxL*HyL*srX(iL,jL) + GxR*HyL*srX(iR,jL)
     &                  + GxL*HyR*srX(iL,jR) + GxR*HyR*srX(iR,jR)

     &                  + HxL*GyL*srY(iL,jL) + HxR*GyL*srY(iR,jL)
     &                  + HxL*GyR*srY(iL,jR) + HxR*GyR*srY(iR,jR)

     &                  + GxL*GyL*sXY(iL,jL) + GxR*GyL*sXY(iR,jL)
     &                  + GxL*GyR*sXY(iL,jR) + GxR*GyR*sXY(iR,jR)
          else
# if KIND == 8
            targ(ic,jc) = 0.D0
# else
            targ(ic,jc) = 0.
# endif
          endif
        enddo
      enddo
# ifdef VERBOSE
      write(*,*) '...leaving spln2d_interp_TYPE_tile'
# endif
      end

#if KIND == 4
# undef KIND
# define KIND 8
# undef spln2d_TYPE
# undef spln2d_interp_TYPE
# undef spln2d_interp_TYPE_tile
# include "spln2d.F"
#endif
