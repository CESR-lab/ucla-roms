      module mrg_zone_comm_vars
      integer :: trd_count=0, npts=0, npts_bak=-1
      end module mrg_zone_comm_vars

#define WESTERN_EDGE istr==1
#define EASTERN_EDGE iend==ncx
#define SOUTHERN_EDGE jstr==1
#define NORTHERN_EDGE jend==ncy

      subroutine etch_mgz_weights_thread(ncx,ncy, mask, mgz,wrk, wdth)

! Assuming that "mgz" is initialized as mgz=1 just on the row of
! perimeter points along the unmasked part of open boundary etch it
! into the interior of the domain by increasing "mgz" values at
! unmasked points by 1 at every step in the location where the value
! iself or at least one of its 4 immediate neighbors is already
! positive.  After completion this procedure creates a constant-slope
! shape function function connected to the boundary by water points.
! Argument "mask" is input-only land mask; "mgz" is input-output; 
! "wrk" is used as work array;   

      implicit none
      integer ncx,ncy, wdth
      integer(kind=2), dimension(ncx,ncy) :: mask, mgz,wrk
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     &             range, istr,iend,jstr,jend, iter, max_iters, i,j
C$    integer omp_get_num_threads, omp_get_thread_num
      numthreads=1                     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c***  nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

! Note that "mgz" and "wrk" switch roles in the two code segments
! below: output of the first becomes input of the second, and vice
! versa. Other than that and tile reversal the two are identical.

      max_iters=wdth/2
      if (wdth > 2*max_iters) max_iters=max_iters+1
      do iter=1,max_iters
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          call etch_weights_tile(istr,iend,jstr,jend, ncx,ncy,
     &                                         mask, mgz,wrk)
        enddo
C$OMP BARRIER
        do tile=my_first,my_last+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)

          if (iter < max_iters  .or.  wdth == 2*max_iters) then
            call etch_weights_tile(istr,iend,jstr,jend, ncx,ncy,
     &                                           mask, wrk,mgz)
          else
            do j=jstr,jend        ! Copy the outcome from the previous
              do i=istr,iend      ! call to etch_weights_tile back into
                mgz(i,j)=wrk(i,j) ! array "mgz" in the case when "wdth"
              enddo               ! is an odd number.  No need to do so
            enddo                 ! for even "wdth" as the final state
          endif                   ! is naturally there.
        enddo
C$OMP BARRIER
      enddo
      end


      subroutine etch_weights_tile(istr,iend,jstr,jend, ncx,ncy,
     &                                            mask, ms1,ms2)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, i,j, iw,ie,js,jn
      integer(kind=2), dimension(ncx,ncy) :: mask, ms1, ms2
      integer(kind=2), parameter :: i1=1
      do j=jstr,jend
        js=max(j-1,1) ; jn=min(j+1,ncy)
        do i=istr,iend
          iw=max(i-1,1) ; ie=min(i+1,ncx)
          if ( mask(i,j) > 0 .and. ( ms1(i,j) > 0 .or.
     &           ms1(iw,j) > 0 .or. ms1(ie,j) > 0   .or.
     &           ms1(i,js) > 0 .or. ms1(i,jn) > 0 )) then
            ms2(i,j)=ms1(i,j)+i1
          else
            ms2(i,j)=0
          endif
        enddo
      enddo
      end



      subroutine mrg_zone_cont_thread(ncx,ncy, mgz,ms1,ms2)

! Enforce the property that every water point within the merging
! zone identified as mgz=1 can be reached from the boundary by water
! (unreachable points mgz=1 are reset to mgz=0).   The procedure is
! essentially the same as in "single_connect" with the exception that
! the initial points are set on the perimeter, and, in principle,
! there is a possibility of multiple unconnected merging zones
! separated from each other by land.

      use mrg_zone_comm_vars
      implicit none
      integer ncx,ncy
      integer(kind=2), dimension(ncx,ncy) :: mgz, ms1, ms2

      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     &                 range, istr,iend,jstr,jend, iter, my_sum, i,j
C$    integer omp_get_num_threads, omp_get_thread_num
      numthreads=1                     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c***  nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      if (trd == 0) then
        write(*,*) 'Enforcing continuity of merging zone, ',
     &                           'numthreads =', numthreads
        write(*,*) 'total number of points in grid',
     &         (ncx-2)*(ncy-2), ' excluding perimeter rows'
        trd_count=0 ; npts=0 ; npts_bak=-1   !<-- initialize
      endif

      do tile=my_first,my_last
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                istr,iend,jstr,jend)
        do j=jstr,jend
          do i=iend,istr
            ms1(i,j)=0 ; ms2(i,j)=0 !<-- initialize
          enddo
        enddo
        if (WESTERN_EDGE) then                   ! initialize search
          do j=jstr,jend                         ! by setting merging
            if (mgz(istr,j) > 0) then            ! zone masks along
              ms1(istr,j)=1 ; ms2(istr,j)=1      ! the perimeter
            endif
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
            if (mgz(iend,j) > 0) then
              ms1(iend,j)=1 ; ms2(iend,j)=1
            endif
          enddo
        endif
        if (SOUTHERN_EDGE) then
          do i=istr,iend
            if (mgz(i,jstr) > 0) then
              ms1(i,jstr)=1 ; ms2(i,jstr)=1
            endif
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
            if (mgz(i,jend) > 0) then
              ms1(i,jend)=1 ; ms2(i,jend)=1
            endif
          enddo
        endif
      enddo  !<-- tile
C$OMP BARRIER

      iter=0
      do while(npts /= npts_bak)
        my_sum=0 ; iter=iter+1
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          if (WESTERN_EDGE) istr=istr+1
          if (EASTERN_EDGE) iend=iend-1
          if (SOUTHERN_EDGE) jstr=jstr+1
          if (NORTHERN_EDGE) jend=jend-1

          do j=jstr,jend
            do i=istr,iend
              if ( mgz(i,j) > 0 .and. ( ms1(i,j ) > 0 .or.
     &             ms1(i-1,j) > 0  .or. ms1(i+1,j) > 0   .or.
     &             ms1(i,j-1) > 0  .or. ms1(i,j+1) > 0 )) then
                ms2(i,j)=1
                my_sum=my_sum+1
              endif
            enddo
          enddo
        enddo  !<-- tile
C$OMP CRITICAL(cr_region)
        if (trd_count == 0) then
          npts_bak=npts ; npts=0
        endif
        npts=npts+my_sum
        trd_count=trd_count+1
        if (trd_count == numthreads) then
          trd_count=0
          if (mod(iter,20) == 0 .or. npts == npts_bak) then
            write(*,'(8x,A,I7,2(2x,A,I10))') 'iter =', iter,
     &           'npts =', npts, 'changes =', npts-npts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER

        my_sum=0 ; iter=iter+1
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          if (WESTERN_EDGE) istr=istr+1
          if (EASTERN_EDGE) iend=iend-1
          if (SOUTHERN_EDGE) jstr=jstr+1
          if (NORTHERN_EDGE) jend=jend-1

          do j=jstr,jend
            do i=istr,iend
              if ( mgz(i,j) > 0 .and. ( ms2(i,j ) > 0 .or.
     &             ms2(i-1,j) > 0  .or. ms2(i+1,j) > 0   .or.
     &             ms2(i,j-1) > 0  .or. ms2(i,j+1) > 0 )) then
                ms1(i,j)=1
                my_sum=my_sum+1
              endif
            enddo
          enddo
        enddo
C$OMP CRITICAL(cr_region)
        if (trd_count == 0) then
          npts_bak=npts ; npts=0
        endif
        npts=npts+my_sum
        trd_count=trd_count+1
        if (trd_count == numthreads) then
          trd_count=0
          if (mod(iter,20) == 0 .or. npts == npts_bak) then
            write(*,'(8x,A,I7,2(2x,A,I10))') 'iter =', iter,
     &           'npts =', npts, 'changes =', npts-npts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo  !<-- while

      do tile=my_last,my_first,-1
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                istr,iend,jstr,jend)
        do j=jstr,jend
          do i=istr,iend                     ! cancel points
             if (ms1(i,j) == 0) mgz(i,j)=0   ! which cannot be
          enddo                              ! reached by water
        enddo
      enddo
C$OMP BARRIER
      end


      subroutine etch_mgz_into_land_thread(ncx,ncy, mask, mgz,ms1,ms2)
      use mrg_zone_comm_vars
      implicit none
      integer ncx,ncy
      integer(kind=2), dimension(ncx,ncy) :: mask, mgz,ms1,ms2
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     &    range, istr,iend,jstr,jend, iter, my_sum, i,j, iw,ie,js,jn
C$    integer omp_get_num_threads, omp_get_thread_num
      numthreads=1                     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c**   nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first+range-1, nsub_x*nsub_y-1)

      if (trd == 0) then
        write(*,*) 'Etching merging zone area into land...'
        trd_count=0 ; npts=0 ; npts_bak=-1   !<-- initialize
      endif

      do tile=my_first,my_last
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)
        do j=jstr,jend
          do i=istr,iend
            if (mask(i,j) > 0) then            ! initialize etching
              if (mgz(i,j) == 0) mgz(i,j)=-1   ! procedure: -1
            endif                              !  -1 = water
            ms1(i,j)=mgz(i,j)                  !   0 = land
            ms2(i,j)=mgz(i,j)                  !  +1 = merging zone
          enddo
        enddo
      enddo
C$OMP BARRIER

      iter=0
      do while(npts /= npts_bak)
        iter=iter+1 ; my_sum=0
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              ms2(i,j)=ms1(i,j)
              if (ms1(i,j) == 0) then
                if (ms1(iw,j)+ms1(ie,j)+ms1(i,js)+ms1(i,jn) < 0) then
                  if ( ms1(iw,j) <= 0 .and. ms1(ie,j) <= 0 .and.
     &                 ms1(i,js) <= 0 .and. ms1(i,jn) <= 0) then
                    ms2(i,j)=-1
                    my_sum=my_sum+1
                  endif
                elseif (ms1(iw,j)+ms1(ie,j)+ms1(i,js)+ms1(i,jn) > 0
     &                                                         ) then
                  if ( ms1(iw,j) >= 0 .and. ms1(ie,j) >= 0 .and.
     &                 ms1(i,js) >= 0 .and. ms1(i,jn) >= 0) then
                    ms2(i,j)=+1
                    my_sum=my_sum+1
                  endif
                endif
              endif
            enddo
          enddo
        enddo !<-- tile
C$OMP BARRIER

        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              ms1(i,j)=ms2(i,j)
              if (ms2(i,j) == 0) then
                if (ms2(iw,j)+ms2(ie,j)+ms2(i,js)+ms2(i,jn) < 0) then
                  if ( ms2(iw,j) <= 0 .and. ms2(ie,j) <= 0 .and.
     &                 ms2(i,js) <= 0 .and. ms2(i,jn) <= 0) then
                    ms1(i,j)=-1
                    my_sum=my_sum+1
                  endif
                elseif (ms2(iw,j)+ms2(ie,j)+ms2(i,js)+ms2(i,jn) > 0
     &                                                         ) then
                  if ( ms2(iw,j) >= 0 .and. ms2(ie,j) >= 0 .and.
     &                 ms2(i,js) >= 0 .and. ms2(i,jn) >= 0) then
                    ms1(i,j)=+1
                    my_sum=my_sum+1
                  endif
                endif
              endif
            enddo
          enddo
        enddo !<-- tile

C$OMP CRITICAL(cr_region)
        if (trd_count == 0) then
          npts_bak=npts ; npts=my_sum
        else
          npts=npts+my_sum
        endif
        trd_count=trd_count+1
        if (trd_count == numthreads) then
          trd_count=0
          if (mod(iter,20) == 0 .or. npts == npts_bak) then
            write(*,'(8x,A,I7,2(2x,A,I10))') 'iter =', iter,
     &           'npts =', npts, 'changes =', npts-npts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo !<-- while
      end




      subroutine etch_weights_into_land_thread(wdth, ncx,ncy, mgz,
     &             ms1,ms2, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
      implicit none
      integer ncx,ncy, wdth
      integer(kind=2), dimension(ncx,ncy) :: mgz, ms1,ms2
      logical OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     &                  range, istr,iend,jstr,jend, iter, max_iters
C$    integer omp_get_num_threads, omp_get_thread_num
      numthreads=1                     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c***  nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)
      if (trd == 0) then
        write(*,*) 'Etching merging weight function into land...'
      endif
      max_iters=(wdth+1)/2
      do iter=1,max_iters
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          call mgz_step_into_land_tile( istr,iend,jstr,jend,
     &                                 ncx,ncy, mgz, ms1,ms2,
     &           wdth, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
        enddo
C$OMP BARRIER
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          call mgz_step_into_land_tile( istr,iend,jstr,jend,
     &                                 ncx,ncy, mgz, ms2,ms1,
     &           wdth, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
        enddo
C$OMP BARRIER
        if (trd == 0 .and. mod(iter,20) == 0) then
          write(*,'(8x,A,I7)') 'iter =', iter
        endif
      enddo
      end

      subroutine mgz_step_into_land_tile(istr,iend,jstr,jend, ncx,ncy,
     &      mgz, ms1,ms2, wdth, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)

! Perform one step of etching "ms2" into allowed land area using "ms1"
! as input.  Note that in the calling sequence above "ms1" and "ms2"
! alternate, so the algorithm here is quasi time stepping.
! At the beginning of this procedure "ms1" and "ms2" are initialized
! as integer constant-slope functions: starting with the maximum value
! of "width" at the open boundary and decreasing by 1 for every row of
! points as proceeding into the interior of the domain. However, this
! applies only to water points; on all the land points both "ms1" and
! "ms2" are set to zero;   Condition "mgz > 0" specifies merging zone
! which already extended into some land (etched), so now "ms1" and
! "ms2" will be allowed to intrude into land areas where "mgz > 0" in
! such a way that their values will be decreasing by 1 every step
! moving away from the perimeter of the grid and from the coastline
! into the land, however, the expansion is not be allowed to go
! beyond "mgz > 0".

      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, wdth
      integer(kind=2), dimension(ncx,ncy) :: mgz, ms1,ms2
      logical OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH
      integer i,j, iw,ie,js,jn, imgz,jmgz
      integer(kind=2), parameter :: i1=1

      if (OBC_WEST) then
        imgz=min(iend, wdth)
        if (istr <= imgz) then
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,imgz
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max(  ms2(i,j), ms1(i,j), ms1(ie,j),
     &                max(ms1(iw,j), ms1(i,jn), ms1(i,js))-i1
     &                                                      )
              endif
            enddo
          enddo
        endif
      endif
      if (OBC_SOUTH) then
        jmgz=min(jend, wdth)
        if (jstr <= jmgz) then
          do j=jstr,jmgz
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max( ms2(i,j), ms1(i,j), ms1(i,jn),
     &               max(ms1(i,js), ms1(ie,j), ms1(iw,j))-i1
     &                                                     )
              endif
            enddo
          enddo
        endif
      endif
      if (OBC_EAST) then
        imgz=max(istr, ncx-wdth+1)
        if (imgz <= iend) then
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=imgz,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max( ms2(i,j), ms1(i,j), ms1(iw,j),
     &               max(ms1(ie,j), ms1(i,jn), ms1(i,js))-i1
     &                                                     )
              endif
            enddo
          enddo
        endif
      endif
      if (OBC_NORTH) then
        jmgz=max(jstr, ncy-wdth+1)
        if (jmgz <= jend) then
          do j=jmgz,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max( ms2(i,j), ms1(i,j), ms1(i,js),
     &               max(ms1(i,jn), ms1(ie,j), ms1(iw,j))-i1
     &                                                     )
              endif
            enddo
          enddo
        endif
      endif
      end


      subroutine smooth_wgt_thread(max_iters, ncx,ncy, mgz,wgt,hwg)
      implicit none
      integer max_iters, ncx,ncy
      integer(kind=2) mgz(ncx,ncy)
      real(kind=8) wgt(ncx,ncy), hwg (ncx,ncy)

      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last,
     &                  tile, range, istr,iend,jstr,jend, iter
C$    integer omp_get_num_threads, omp_get_thread_num
      numthreads=1                     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)

c***  nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      if (trd == 0) then
        write(*,*) 'Smoothing weight function, iters =',
     &          max_iters, '  numthreads =', numthreads
      endif

      do iter=1,max_iters
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          call gauss_sm_wgt_tile(istr,iend,jstr,jend,
     &                                 ncx,ncy, mgz, wgt,hwg)
        enddo
C$OMP BARRIER
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          call smooth_wgt_tile(istr,iend,jstr,jend,
     &                                 ncx,ncy, mgz, hwg,wgt)
        enddo
C$OMP BARRIER
        if (mod(iter,20) == 0 .and. trd == 0) then
          write(*,'(8x,A,I7)') 'iter =', iter
        endif
      enddo
      end

      subroutine smooth_wgt_tile(istr,iend,jstr,jend,
     &                          ncx,ncy, mgz, src,targ)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, i,j
      integer(kind=2) mgz(ncx,ncy)
      real(kind=8) src(ncx,ncy), targ(ncx,ncy)

      if (WESTERN_EDGE) istr=istr+1   ! isotropic smoothing operator
      if (EASTERN_EDGE) iend=iend-1   !
      if (SOUTHERN_EDGE) jstr=jstr+1  !        1/32  1/8  1/32
      if (NORTHERN_EDGE) jend=jend-1  !
                                      !        1/8   3/8   1/8
      do j=jstr,jend                  !
        do i=istr,iend                !        1/32  1/8  1/32
          if (mgz(i,j) > 0) then
            targ(i,j)=0.125D0*( 3.D0*src(i,j)+src(i-1,j)+src(i,j-1)
     &                                        +src(i+1,j)+src(i,j+1)
     &                           +0.25D0*( src(i-1,j+1)+src(i+1,j+1)
     &                                    +src(i-1,j-1)+src(i+1,j-1)
     &                                                          ) )
          else
            targ(i,j)=0.D0
          endif
        enddo
      enddo
      end


      subroutine gauss_sm_wgt_tile(istr,iend,jstr,jend,
     &                          ncx,ncy, mgz, src,targ)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, i,j
      integer(kind=2) mgz(ncx,ncy)
      real(kind=8) src(ncx,ncy), targ(ncx,ncy)

      if (WESTERN_EDGE) istr=istr+1 ; if (SOUTHERN_EDGE) jstr=jstr+1
      if (EASTERN_EDGE) iend=iend-1 ; if (NORTHERN_EDGE) jend=jend-1

      do j=jstr,jend
        do i=istr,iend
          if (mgz(i,j) > 0) then
            targ(i,j)=0.2D0*(                 src(i-1,j)+src(i,j-1)
     &                                        +src(i+1,j)+src(i,j+1)
     &                           +0.25D0*( src(i-1,j+1)+src(i+1,j+1)
     &                                    +src(i-1,j-1)+src(i+1,j-1)
     &                                                          ) )
          else
            targ(i,j)=0.D0
          endif
        enddo
      enddo
      end
