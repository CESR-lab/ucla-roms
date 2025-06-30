! A set of generic tiling tools for 2D subdomain decomposition needed
! for a coarse-grained shared-memory OpenMP parallelization approach
! similar to that of the actual ROMS code with the exception that now
! it is adapted for dynamically changing dimensions, so the numbers of
! partitions are selected automatically from the grid array dimensions
! and number of threads, with an effort to minimize the miss-balance
! (if the dimensions are not evenly divisible by the desired number of
! tiles), as well as to bring subdomain sizes close to the optimum for
! cache performance, subject to judicial compromise in the case of
! conflict between the two goals.
!
! Content:
!
!   set_tiles(nx,ny, nsub_x,nsub_y)
!
!         input: nx,ny -- grid dimensions;
!                number of threads (determined internally);
!
!        output: nsub_x  -- numbers of tiles in each direction (their
!                nsub_y     product is guaranteed to be divisible by
!                           the number of threads;
!
!   comp_tile_bounds(tile, nx,ny,nsub_x,nsub_y, istr,iend,jstr,jend)
!
!         input: tile -- tile number ranging from 0 to nsub_x*nsub_y-1;
!                nx,ny, nsub_x,nsub_y -- same as above;
!
!        output: istr,iend,jstr,jend -- starting and ending indices
!                                       of subdomain "tile";

c--#define TEST     /* For self-testing only */
#ifdef TEST
      program test_tiling
      implicit none
      integer nx,ny, nsub_x,nsub_y
  1    write(*,*) 'enter nx,ny:'
       read(*,*) nx,ny
C$OMP PARALLEL SHARED(nx,ny)
         call set_tiles(nx,ny, nsub_x,nsub_y)
C$OMP END PARALLEL
       goto 1
      end
#endif

      subroutine set_tiles(nx,ny, nsub_x,nsub_y)
      implicit none
      integer nx,ny, nsub_x,nsub_y
      integer ntrds, ntx,nty, nsb, i,size,excess,max_exc
!C$    integer omp_get_num_threads
      integer,parameter :: min_vect_lenght=96, targ_length=128

      ntrds=1
!C$    ntrds=omp_get_num_threads()

      if (nx>8*min_vect_lenght-1 .and. mod(ntrds,8)==0) then
        ntx=8
      elseif (nx>6*min_vect_lenght-1 .and. mod(ntrds,6)==0) then
        ntx=6
      elseif (nx>4*min_vect_lenght-1 .and. mod(ntrds,4)==0) then
        ntx=4
      elseif (nx>3*min_vect_lenght-1 .and. mod(ntrds,3)==0) then
        ntx=3
      elseif (nx>2*min_vect_lenght-1 .and. mod(ntrds,2)==0) then
        ntx=2
      else                        ! Scan possible tile sizes within
        ntx=1                     ! the range of "i" index and select
      endif                       ! one which yields the least number
      nty=ntrds/ntx               ! of excess points, subject to the
                                  ! constraint that "nsb" can be
      max_exc=ny                  ! evenly divided by "nty"...
      do i=9,25
        nsb=(ny+i-1)/i            !<-- prospective tile size
        nsb=nsb-mod(nsb,nty)      !<-- make "nsb" divisible by "nty"
        if (nsb == 0) nsb=nty
        size=(ny+nsb-1)/nsb       !--> actual tile size
        excess=nsb*size-ny
        if (excess<max_exc) then
          max_exc=excess
          nsub_y=nsb
        endif
#ifdef TEST
        write(*,*) nsb,size,excess
#endif
      enddo

      if (mod(nsub_y,ntrds) > 0) then
        nsb=nsub_y/nty
        if (mod(nsb,3)==0 .and. mod(ntx,3)==0) then
          nsb=nsb/3
          ntx=ntx/3
        endif
        do while(mod(nsb,2)==0 .and. mod(ntx,2)==0)
          nsb=nsb/2
          ntx=ntx/2
        enddo
      else
        ntx=1
      endif

      nsub_x=(nx+targ_length-1)/targ_length
      if (ntx > 1) then
        i=mod(nsub_x,ntx)           !--> make sure that "nsub_x"
        if (i > ntx/2) then         !    is divisible by "ntx"
          nsub_x=nsub_x+ntx-i
        else
          nsub_x=nsub_x-i
        endif
      endif

#ifdef TEST
C$OMP CRITICAL (tiling_cr_rgn)
      size=(nx+nsub_x-1)/nsub_x
      excess=nsub_x*size-nx
      write(*,'(4(1x,A,I4))') 'nx =',nx, 'nsub_x=',nsub_x,
     &                  'size =',size, 'excess =', excess
      size=(ny+nsub_y-1)/nsub_y
      excess=nsub_y*size-ny
      write(*,'(4(1x,A,I4))') 'ny =',ny, 'nsub_y=',nsub_y,
     &                  'size =',size, 'excess =', excess
      write(*,*)
C$OMP END CRITICAL (tiling_cr_rgn)
#elif defined VERBOSE
C$OMP CRITICAL (tiling_cr_rgn)
      write(*,'(2(1x,A,I4),3(1x,A,I3))') 'set_tiles: nx =', nx,
     &        'ny =',ny, 'numthreads =', ntrds, 'nsub_x =', nsub_x,
     &                                          'nsub_y =', nsub_y
C$OMP END CRITICAL (tiling_cr_rgn)
#endif
      end

      subroutine comp_tile_bounds(tile, nx,ny,nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
      implicit none
      integer tile, nx,ny,nsub_x,nsub_y, istr,iend,jstr,jend,
     &                                              i,j,size
      j=tile/nsub_x
      i=tile-nsub_x*j                    ! This tiling algorithm
                                         ! generally follows ROMS,
      size=(nx+nsub_x-1)/nsub_x          ! where tile sizes in each
      istr=1+i*size -(nsub_x*size-nx)/2  ! direction are computed by
      iend=min(istr+size-1,nx)           ! integer division rounded
      istr=max(istr,1)                   ! upward, while the excess
                                         ! (in the case when the
      size=(ny+nsub_y-1)/nsub_y          ! domain cannot be divided
      jstr=1+j*size -(nsub_y*size-ny)/2  ! evenly) is half-and-half
      jend=min(jstr+size-1,ny)           ! split between the ends.
      jstr=max(jstr,1)                   ! As the result tiles the
                                         ! first and the last in each
      end                                ! direction may be smaller.
