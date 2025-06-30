      subroutine r2r_interp_init_thread( nx,ny,x,y, ncx,ncy, xc,yc,
     &                                                ip,jp, xi,eta)

! Initialize interpolation between two arbitrary oriented curvilinear
! grids: given arrays of coordinates [x(nx,ny),y(nx,ny)] for points of
! a non-Cartesian "source" grid, and similar [xc(ncx,ncy),yc(ncx,ncy)]
! for the "target" grid (hereafter "c" stands for "child") find arrays
! of indices [ip(ncx,ncy),jp(ncx,ncy)] and arrays fractional distances
! [xi(ncx,ncy),eta(ncx,ncy)] such that each point of the target grid
! (xc,yc) is surrounded by 4 points of the source,
!
!              [x,y](ip,jp+1) --__
!                     /             --__ [x,y](ip+1,jp+1)
!                    /   [xc,yc]            /
!                   /                      /
!              [x,y](ip,jp) -- __         /
!                               --__ [x,y](ip+1,jp)
!
! such that bi-linear interpolation of [x,y] into location of [xc,yc]
! yields [xc,yc] themselves,
!
!             xc = (1-xi)*(1-eta)*x(i,j) + xi*(1-eta)*x(i+1,j)
!                 +(1-xi)* eta *x(i,j+1) + xi* eta *x(i+1,j+1)
!
!             yc = (1-xi)*(1-eta)*y(i,j) + xi*(1-eta)*y(i+1,j)
!                 +(1-xi)* eta *y(i,j+1) + xi* eta *y(i+1,j+1)
!
! in other words, (ip+xi,jp+eta) are coordinates of point (xc,yc) in
! "continuous index" space (i,j) of the source grid.

! Normally the target grid is expected to be entirely within the area
! covered by the source, however it is not strictly required by the
! algorithm below: if some portions of the target grid are outside the
! source (hence interpolation is impossible for these locations) then
! ip=-1 and jp=-1 are set as special values.

      implicit none
      integer nx,ny, ncx,ncy
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer icmin,icmax,jcmin,jcmax
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncy+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx ; jcmin=1 ; jcmax=ncy
C$    jcmin=1+trd*chunk_size ; jcmax=min(jcmin+chunk_size-1,ncy)

      call search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     &                                       ncx,ncy, xc,yc, ip,jp)

      call comp_offsets_tile(  icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     &                               ncx,ncy, xc,yc, ip,jp, xi,eta)
      end



      subroutine r2r_interp_search_thread( nx,ny,x,y, ncx,ncy, xc,yc,
     &                                                          ip,jp)
      implicit none
      integer nx,ny, ncx,ncy                       ! Search for parent
      real(kind=8), dimension(nx,ny) :: x,y        ! grid indices only,
      real(kind=8), dimension(ncx,ncy) :: xc,yc    ! do not compute
      integer(kind=4), dimension(ncx,ncy) :: ip,jp ! xi,eta offsets.
      integer icmin,icmax,jcmin,jcmax
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncy+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx ; jcmin=1 ; jcmax=ncy
C$    jcmin=1+trd*chunk_size ; jcmax=min(jcmin+chunk_size-1,ncy)

      call search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     &                                       ncx,ncy, xc,yc, ip,jp)
      end



      subroutine search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny,
     &                                  x,y, ncx,ncy, xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc         !   4 <---- 3
      integer(kind=4), dimension(ncx,ncy) :: ip,jp      !   !       ^
      integer i,j, ic,jc, inew,jnew, iinc               !   !       !
      real(kind=8) dx4,dy4, dx3,dy3, r41,r34,           !   v       !
     &             dx1,dy1, dx2,dy2, r12,r23            !   1 ----> 2
#ifdef VERBOSE
      write(*,*) 'enter search_tile, jcmin =',jcmin, ' jcmax =', jcmax
#endif

! The search algorithm is organized as follows:
!
!    initialization: select starting point ic,jc
!    and compute cyclic vector products r12...r41
!
!    do while(.true.)
!      do while( r12...r41 >= 0 -- meaning that point
!                 ic,jc is inside of [i:i+1]x[j:j+1] )
!
!        record ip,jp for this point, and proceed to the next point
!        of target grid using by incrementing only ONE index, either
!        "ic" or "jc" -- if "ic" already reached its bound, while
!        simultaneously reversing the direction of "ic" sweep -- the
!        reversing is done to make sure that the next target point is
!        always nearby, to maximize the probability that it is still
!        inside of [i:i+1]x[j:j+1] cell, so the inner while-loop would
!        not break off too often. Recompute r12...r41 for the next
!        ic,jc point to check logical condition during next iteration
!        of while loop; [Note that there is no attempt to change
!        source-grid indices i,j inside this loop.]
!
!      enddo
!
!      once the while(r12...r41 >= 0) loop breaks off, it means that
!      the point ic,jc is no longer inside [i:i+1]x[j:j+1], so i and/or
!      j index of source grid must be incremented, depending which of
!      the r12...r41 is negative.  The indices i,j are incremented by
!      one and are restricted to be within the range of source grid.
!      This leads to two possibilities:
!      either
!             (i) at least one of them, i,j is incremented, then keep
!                 ic,jc the same and proceed with recomputing r12...r41
!                 to resume while(r12...r41 >= 0) loop;
!      or
!            (ii) both increments of i and j are canceled by the
!                 restriction, so i,j get "stuck"  which is detected
!                 by inew,jnew both having the same values. This means
!                 that this ic,jc point cannot be bounded properly
!                 because its location is outside the source grid, so
!                 mark it by special values of ip,jp (essentially skip
!                 this point) and proceed to the next one by
!                 incrementing  ip,jp;
!
!    enddo !<-- while(.true.)
!
! Note (a) that not having special care for (i) leads to an infinite
!          loop in the algorithm if some of the target points cannot be
!          bounded; and
!      (b) assuming that the target grid has finer resolution than the
!          source, it is expected that most of the computing time is
!          spent inside while(r12...r41 >= 0) loop, hence its body is
!          minimized at the expense of somewhat awkwardness of what is
!          around it.

      i=nx/2 ; j=ny/2 ; ic=icmin ; jc=jcmin ; iinc=+1

      dx1=x(i+1,j  )-xc(ic,jc) ; dy1=y(i+1,j  )-yc(ic,jc)   !  3 <-- 2
      dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)   !  !     !
      dx3=x(i  ,j+1)-xc(ic,jc) ; dy3=y(i  ,j+1)-yc(ic,jc)   !  !     !
      dx4=x(i  ,j  )-xc(ic,jc) ; dy4=y(i  ,j  )-yc(ic,jc)   !  4 --> 1

      r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
      r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      do while(.true.)
        do while( r12 >= 0.D0 .and. r23 >= 0.D0
     &      .and. r34 >= 0.D0 .and. r41 >= 0.D0 )

          ip(ic,jc)=i ; jp(ic,jc)=j        !<-- record bounding indices

          if (iinc > 0 .and. ic < icmax) then    !--> proceed to the
            ic=ic+1                              !    next target point
          elseif (iinc < 0 .and. ic > icmin) then
            ic=ic-1
          else
            jc=jc+1 ; iinc=-iinc ; if (jc > jcmax) return
#ifdef VERBOSE
            write(*,*) 'jc =', jc
#endif
          endif

          dx1=x(i+1,j  )-xc(ic,jc) ; dy1=y(i+1,j  )-yc(ic,jc)
          dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
          dx3=x(i  ,j+1)-xc(ic,jc) ; dy3=y(i  ,j+1)-yc(ic,jc)
          dx4=x(i  ,j  )-xc(ic,jc) ; dy4=y(i  ,j  )-yc(ic,jc)

          r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
          r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4
        enddo

        inew=i ; jnew=j
        if (r12 < 0.D0) inew=min(inew+1, nx-1)
        if (r23 < 0.D0) jnew=min(jnew+1, ny-1)
        if (r34 < 0.D0) inew=max(inew-1,    1)
        if (r41 < 0.D0) jnew=max(jnew-1,    1)

        if (inew == i .and. jnew == j) then
          ip(ic,jc)=-1 ; jp(ic,jc)=-1       !<-- cannot be bounded

          if (iinc > 0 .and. ic < icmax) then  !--> proceed to the
            ic=ic+1                            !    next target point
          elseif (iinc < 0 .and. ic > icmin) then
            ic=ic-1
          else
            jc=jc+1 ; iinc=-iinc ; if (jc > jcmax) return
#ifdef VERBOSE
            write(*,*) 'jc =', jc
#endif
          endif
        else
          i=inew ; j=jnew       !<-- accept move of source location.
        endif

        dx1=x(i+1,j  )-xc(ic,jc) ; dy1=y(i+1,j  )-yc(ic,jc)
        dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
        dx3=x(i  ,j+1)-xc(ic,jc) ; dy3=y(i  ,j+1)-yc(ic,jc)
        dx4=x(i  ,j  )-xc(ic,jc) ; dy4=y(i  ,j  )-yc(ic,jc)

        r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
        r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      enddo
      end



      subroutine comp_offsets_tile(icmin,icmax,jcmin,jcmax, nx,ny,
     &                         x,y, ncx,ncy, xc,yc, ip,jp, xi,eta)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer ic,jc, i,j, iter
      real(kind=8) a11,a12,a21,a22, det, dX,dY, p,p1,q,q1

      do jc=jcmin,jcmax
        do ic=icmin,icmax
          i=ip(ic,jc) ; j=jp(ic,jc)
          if (i > 0 .and. j > 0) then
            a11=0.5D0*(x(i+1,j+1)-x(i,j+1) +x(i+1,j)-x(i,j))
            a12=0.5D0*(x(i+1,j+1)-x(i+1,j) +x(i,j+1)-x(i,j))
            a21=0.5D0*(y(i+1,j+1)-y(i,j+1) +y(i+1,j)-y(i,j))
            a22=0.5D0*(y(i+1,j+1)-y(i+1,j) +y(i,j+1)-y(i,j))

            dX=xc(ic,jc)-0.25D0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))
            dY=yc(ic,jc)-0.25D0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))

            det=1.D0/(a11*a22-a12*a21)
            xi(ic,jc) =0.5D0 + det*(a22*dX-a12*dY)
            eta(ic,jc)=0.5D0 + det*(a11*dY-a21*dX)
          else
            xi(ic,jc)=-1.D0 ; eta(ic,jc)=-1.D0  !<-- special values
          endif
        enddo
        do iter=1,10
          do ic=icmin,icmax
            i=ip(ic,jc) ; j=jp(ic,jc)
            if (i > 0 .and. j > 0) then
              p=xi(ic,jc)  ; p1=1.D0-p
              q=eta(ic,jc) ; q1=1.D0-q

              a11=q*(x(i+1,j+1)-x(i,j+1)) +q1*(x(i+1,j)-x(i,j))
              a12=p*(x(i+1,j+1)-x(i+1,j)) +p1*(x(i,j+1)-x(i,j))
              a21=q*(y(i+1,j+1)-y(i,j+1)) +q1*(y(i+1,j)-y(i,j))
              a22=p*(y(i+1,j+1)-y(i+1,j)) +p1*(y(i,j+1)-y(i,j))

              dX=xc(ic,jc) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     &                      -p*q1*x(i+1,j)  -p1*q1*x(i,j)
              dY=yc(ic,jc) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     &                      -p*q1*y(i+1,j)  -p1*q1*y(i,j)

              det=1.D0/(a11*a22-a12*a21)

              xi(ic,jc) =p + det*(a22*dX-a12*dY)
              eta(ic,jc)=q + det*(a11*dY-a21*dX)
            endif
          enddo
        enddo
      enddo
      end


      subroutine check_search_indices(nx,ny,x,y, ncx,ncy, xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx,ncy
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer ic,jc, i,j
      real(kind=8) dx1,dy1, dx2,dy2, dx3,dy3, dx4,dy4, r12,r23,r34,r41

      write(*,'(2x,A)',advance='no') 'check_search_indices ...'
      do jc=1,ncy
        do ic=1,ncx
          i=ip(ic,jc) ; j=jp(ic,jc)
          if (i > 0 .and. j > 0) then
            dx1=x(i+1,j  )-xc(ic,jc) ; dy1=y(i+1,j  )-yc(ic,jc)
            dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
            dx3=x(i  ,j+1)-xc(ic,jc) ; dy3=y(i  ,j+1)-yc(ic,jc)
            dx4=x(i  ,j  )-xc(ic,jc) ; dy4=y(i  ,j  )-yc(ic,jc)

            r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
            r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

            if (r12<0.D0 .or. r23<0.D0 .or. r34<0.D0 .or. r41<0.D0)
     &        write(*,*) '### ERROR: Search algorithm failure at ',
     &                                    'ic =', ic, ' jc =', jc
          endif
        enddo
      enddo
      write(*,'(2x,A)') '...done'
      end



      subroutine check_offsets(nx,ny, x,y, ncx,ncy, xc,yc,
     &                                      ip,jp, xi,eta)
      implicit none
      integer nx,ny, ncx,ncy
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer ic,jc, i,j
      real(kind=8) dX,dY, p,p1,q,q1, errX,errY

      write(*,'(2x,A)',advance='no') 'checking offsets...'
      errX=0.D0 ; errY=0.D0
      do jc=1,ncy
        do ic=1,ncx                     ! Notice simple semantic rule
          i=ip(ic,jc) ; j=jp(ic,jc)     ! here: one and only one "1" is
          if (i > 0 .and. j > 0) then   ! always present:  either as +1
            p=xi(ic,jc)  ; p1=1.D0-p    ! in index, or as p1,q1 in the
            q=eta(ic,jc) ; q1=1.D0-q    ! corresponding coefficient.

            dX=xc(ic,jc) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     &                    -p*q1*x(i+1,j)  -p1*q1*x(i,j)
            dY=yc(ic,jc) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     &                    -p*q1*y(i+1,j)  -p1*q1*y(i,j)

            if (abs(dX) > errX) errX=abs(dX)
            if (abs(dY) > errY) errY=abs(dY)

c**         xi(ic,jc)=dX ; eta(ic,jc)=dY
          endif
        enddo
      enddo
      write(*,'(2x,A,2ES22.15)') 'max errors =', errX,errY
      end
