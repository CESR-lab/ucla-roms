      module r2r_rotate_scratch
       integer, save :: alloc_FX_size=0
C$OMP THREADPRIVATE(alloc_FX_size)
       real(kind=4), allocatable, dimension(:):: FX,FE
C$OMP THREADPRIVATE(FX,FE)
      end module r2r_rotate_scratch 


      subroutine r2r_rotate_shift_thread(nx,ny,N, csA,snA, ur,vr, u,v)

! Rotate vector components (ur,vr) colocated ar RHO-points while
! simultaneously shifting them into U- and V-locations on C-grid. 
! inputs csA,snA, ur,vr (all remain unchanged); outpits u,v

      use r2r_rotate_scratch
      implicit none
      integer nx,ny,N
      real(kind=8) csA(nx,ny),snA(nx,ny)
      real(kind=4) ur(nx,ny,N),vr(nx,ny,N), u(nx-1,ny,N),v(nx,ny-1,N)
!>
      integer ntrds,trd, nsub_x,nsub_y, ntls,isize,jsize,
     &        my_first,my_last, tile, istr,iend,jstr,jend
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()

      call set_tiles(nx,ny, nsub_x,nsub_y)

      isize=(nx+nsub_x-1)/nsub_x ; jsize=(ny+nsub_y-1)/nsub_y
      if (alloc_FX_size < (isize+4)*(jsize+4)) then
        alloc_FX_size=(isize+4)*(jsize+4)
        if (allocated(FX)) deallocate(FX,FE)
        allocate(FX(alloc_FX_size),FE(alloc_FX_size))
C$OMP CRITICAL(ext_wr_rgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)') 'r2r_rotate_shift_thread: ',
     &           'allocated',  float(2*alloc_FX_size)/float(262144),
     &                        'MB of private workspace, trd =', trd
#if VERBOSE > 1
        write(*,'(4x,2(A,I5,1x),2(1x,A,I4),1x,2(1x,A,I4))')
     &  'nx =',nx, 'ny =',ny, 'nsub_x =',nsub_x, 'nsub_y =',nsub_y,
     &                         'isize =',isize,   'jsize =',jsize
#endif
C$OMP END CRITICAL(ext_wr_rgn)
#ifdef VERBOSE
      elseif (trd == 0) then
        write(*,'(1x,A)',advance='no')
     &  'enter r2r_rotate_shift_thread...'
#endif
      endif

      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)

      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, nx,ny, nsub_x,nsub_y,
     &                               istr,iend,jstr,jend )
        call r2r_rotate_shift_tile( istr,iend,jstr,jend,
     &               nx,ny,N, csA,snA, ur,vr, u,v, FX,FE )
      enddo
#ifdef VERBOSE
      if (trd == 0) write(*,'(2x,A)') 'leaving r2r_rotate_shift_thread'
#endif
      end


      subroutine r2r_rotate_shift_tile(istr,iend,jstr,jend, nx,ny,N,
     &                                   csA,snA, ur,vr, u,v, FX,FE)
      implicit none
      integer istr,iend,jstr,jend, nx,ny,N, istrU,jstrV, i,j,k
      real(kind=8) csA(nx,ny),snA(nx,ny)
      real(kind=4) ur(nx,ny,N), u(2:nx,ny,N),
     &             vr(nx,ny,N), v(nx,2:ny,N),
     &           FX(istr-1:iend,jstr-1:jend),
     &           FE(istr-1:iend,jstr-1:jend)
      istrU=max(istr,2) ; jstrV=max(jstr,2)
      do k=1,N
        do j=jstrV-1,jend
          do i=istrU-1,iend
            FX(i,j)=ur(i,j,k)*csA(i,j)+vr(i,j,k)*snA(i,j)
            FE(i,j)=vr(i,j,k)*csA(i,j)-ur(i,j,k)*snA(i,j)
          enddo
        enddo
        do j=jstr,jend
          do i=istrU,iend
            u(i,j,k)=0.5*(FX(i,j)+FX(i-1,j))
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            v(i,j,k)=0.5*(FE(i,j)+FE(i,j-1))
          enddo
        enddo
      enddo  !<-- k
      end



      subroutine r2r_rotate_in_place_thread(nx,ny,N, csA,snA, ur,vr)

! Rotate vector components (ur,vr) colocated ar RHO-points and place
! the outcome back into the same arrays.  

      implicit none
      integer nx,ny,N
      real(kind=8) csA(nx,ny),  snA(nx,ny)
      real(kind=4) ur(nx,ny,N), vr(nx,ny,N)
!>
      real(kind=4) cff1,cff2
      integer ntrds,trd, nsub_x,nsub_y, istr,iend,jstr,jend,
     &                   ntls,tile, my_first,my_last, i,j,k
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      if (trd == 0) write(*,'(2x,A)',advance='no')
     &             'entering r2r_rotate_r_thread...'

      call set_tiles(nx,ny, nsub_x,nsub_y)

      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)

      do tile=my_first,my_last,+1
        call comp_tile_bounds(tile, nx,ny, nsub_x,nsub_y,
     &                               istr,iend,jstr,jend)
        do j=jstr,jend
          do k=1,N
            do i=istr,iend
              cff1=ur(i,j,k)*csA(i,j)+vr(i,j,k)*snA(i,j)
              cff2=vr(i,j,k)*csA(i,j)-ur(i,j,k)*snA(i,j)
              ur(i,j,k)=cff1 ; vr(i,j,k)=cff2
            enddo
          enddo
        enddo
      enddo
      if (trd == 0) write(*,'(2x,A)') '...done'
      end
