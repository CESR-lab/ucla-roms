! Content of this package is user-callable "etch_into_land" and/or
!-------- -- ---- -------  "etch_into_land_thread",

      subroutine etch_into_land(nx,ny, mask, qfld)
      implicit none
      integer nx,ny
      integer(kind=2) mask(nx,ny)
      real(kind=4) qfld(nx,ny)
C$OMP PARALLEL SHARED(nx,ny, mask, qfld)
      call etch_into_land_thread(nx,ny, mask, qfld)
C$OMP END PARALLEL
      end

! which progressively fill in land masked and/or special/missing value
! areas of a given 2D field "qfld" in order to subsequently interpolate
! it by an algorithm not designed to handle missing values.  All other
! routines in this package are for its own internal use and not meant
! to be user callable.  This package replaces "extend_on_land.F" by
! with the same functionality at a significantly reduced computational
! cost.  While mathematically equivalent, this package uses radically
! different approach to code organization (and optimization) by
! replacing 2D-i,j-index sweeps with progressively reduction of a
! precomputed list of indices special-valued points:

!   (1) make list of indices of all points with special-valued
!   (2) select the ones from list (1) above which have immediate
!       neighbors with valid values - these are the points which
!       potentially can be filled in;
!
!   (3) attempt to fill in all points from the list (2) - in order
!       to be filled in a point must pass certain threshold of having
!       enough valid neighbors - not all from list (2) are eventually
!       be selected; compute values for those which will, an place
!       these values into a special array (not directly into "qfld" at
!       this moment) in order to avoid interference); mark positions
!       of these finally selected points within list (1) by zeroing
!       out the first index;
!
!   (4) once step (3) is complete for the entire array and by all
!       the threads, fill in selected values into "qfld";
!
!   (5) shorten the list (1) by excluding zero-index points which
!       were marked (3) and filled in by step (4);
!
!   keep repeating steps (2)-(5) until nothing left in the list.
!

! The roles of the individual routines in the above algorithm are:
!---- -----
!        call copy_extend_tile  ! preliminary step
!        call set_qext_bc_tile
!        call init_ijmsk_tile   !<-- list (1) of all sp.val. points
!        do while(mskd_pts>0)
!          call select_coast_pts  !<-- list (2) as subset of (1)
!          call comp_patch_pts     !<-- list (3) as subset of (2)
!          call apply_patch_pts    !<-- (4)
!          call set_qext_bc_tile
!          call shortlist         !<--(5)
!        enddo
!        call strip_halo_tile

! Meaning of the variables: [dimensional variables "nmsk", "ncst", and
!-------- -- --- ---------- "npths" indicate meaningful portions of the
!                            arrays as rather than allocated sizes;
! qext(0:nx+1,0:ny+1) working array with one row of ghost points all
!           around, initially copy of "qfld" with Neumann b.c. applied;
! ijmsk(2,nmsk) list of i,j-indices of all special-valued points of
!                                 array "qext", this is list (1) above;
! nijmsk(ncst) subset of indices within ijmsk(2,*) identifying points
!               which can be potentially filled during the current
!                                            iteration (list (2) above;
! ijptch(2:npths) list of i,j-indices of "qext" which are about to be
!                    filled; this is list (3) which is a subset nijmsk;
! ptch(npths) a set of fill values computed by step (3) and applied by
!                                                          step(4);

! Parallelization:  Because of the mathematical nature of the problem,
!-----------------  the standard tiling approach is would not be most
! efficient because the work done by each thread depends on the number
! of special valued points rather than total number of points within a
! subdomain, and therefore not predictable. This may cause significant
! load miss-balance.  Note that the algorithm within each thread is
! fundamentally sequential as it involves counting to form index lists.
! So the approach below is to one-dimensionally divide all the domain
! into set of narrow stripes (parameter "jsize" below) and make them
! interleaved among the threads to, hopefully, make all of them get
! approximately the same amount of land.  All the lists of indices are
! strictly PRIVATE in the code below, so each thread is responsible for
! processing only points within its own pre-determined set of stripes.
! Because the first list "ijmsk" is formed by a "_tile" (i.e. striped)
! routine, this establishes AFFINITY between the stripes and the
! subsequent going-by-the-list loops in sense that all the indices
! from the lists are guaranteed to be within the set of stripes of
! the thread.  This eliminates the necessity to have barrier copying
! fill values into "qext" by a going-along-index-list routine and
! applying Neumann b.c. by tiled. Thus, the only operations needed
! to be barrier-separated are the ones where "qext" is input (forming
! lists and calculating fill values) and where "qext" is output
! (applying fill values into "qext").

#define spv -9.9E+9
#define spv_set -1.D+10

#define WESTERN_EDGE istr==1
#define EASTERN_EDGE iend==nx
#define SOUTHERN_EDGE jstr==1
#define NORTHERN_EDGE jend==ny

      module etch_on_land_vars
        real(kind=4), allocatable, dimension(:) :: qext
        integer, save :: allc_ext_size=0, trd_count=0, mskd_pts=0

        integer, save :: alloc_msk_size=0
C$OMP THREADPRIVATE(alloc_msk_size)
        real(kind=4), allocatable, dimension(:) :: ptch
C$OMP THREADPRIVATE(ptch)
        integer(kind=2), allocatable, dimension(:,:) :: ijmsk
C$OMP THREADPRIVATE(ijmsk)
        integer(kind=2), allocatable, dimension(:,:) :: ijptch
C$OMP THREADPRIVATE(ijptch)
        integer(kind=4), allocatable, dimension(:) :: nijmsk
C$OMP THREADPRIVATE(nijmsk)
      end module etch_on_land_vars




      subroutine etch_into_land_thread(nx,ny, mask, qfld)
      use etch_on_land_vars
      implicit none
      integer nx,ny
      integer(kind=2) mask(nx,ny)
      real(kind=4) qfld(nx,ny)
!>
      integer ntrds,trd, nmsk,nmsk_new, ncst, npths,
     &        istr,iend,jstr,jend, j,j0, jskip, ierr
      integer, parameter :: jsize=3
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      j0=trd*jsize ; jskip=ntrds*jsize ; istr=1 ; iend=nx


C$OMP MASTER
      mskd_pts=1!<-- set to start while loop
      if ((nx+2)*(ny+2) > allc_ext_size) then
        allc_ext_size=(nx+2)*(ny+2)
        if (allocated(qext)) deallocate(qext)
        allocate(qext(allc_ext_size))
        write(*,'(1x,2A,F16.8,1x,A)') 'etch_into_land_thread :: ',
     &             'allocated',  dble((nx+2)*(ny+2))/dble(262144),
     &                                'MB shared workspace array'
      endif
C$OMP END MASTER
C$OMP BARRIER

      nmsk=0
      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call copy_extend_tile(istr,iend,jstr,jend,  nx,ny, mask,
     &                                          qfld,qext, nmsk)
         call set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
      enddo
C$OMP BARRIER

      if (nmsk > alloc_msk_size) then
        alloc_msk_size=nmsk
        if (allocated(ptch)) deallocate(ptch,ijptch, nijmsk,ijmsk)
        allocate(  ijmsk(2,alloc_msk_size),  nijmsk(alloc_msk_size),
     &            ijptch(2,alloc_msk_size),    ptch(alloc_msk_size),
     &                                                 stat=ierr )
        if (ierr == 0) then
C$OMP CRITICAL(etch_cr_rgn)
          write(*,'(1x,2A,F16.8,1x,A,I3)') 'etch_into_land_thread :: ',
     &               'allocated',  dble(4*alloc_msk_size)/dble(262144),
     &                             'MB private workspace, trd =', trd
C$OMP END CRITICAL(etch_cr_rgn)
        else
          write(*,*) '### ERROR: etch_into_land_thread :: ',
     &                           'memory allocation error.'
        endif
      endif


      nmsk=0 !<-- number of masked points
      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call init_ijmsk_tile(istr,iend,jstr,jend, nx,ny, qext,
     &                            alloc_msk_size, ijmsk, nmsk)
      enddo
C$OMP BARRIER

      do while (mskd_pts > 0)
        ncst=0 !<-- number of coastal points
        call select_coast_pts(nx,ny, qext, nmsk,ijmsk, ncst,nijmsk)

        npths=0 !<-- number of points to be patched
        call comp_patch_pts(nx,ny, qext, nmsk,ijmsk, ncst,nijmsk,
     &                                       npths, ijptch,ptch)
C$OMP BARRIER
        call apply_patch_pts(npths, ijptch,ptch,   nx,ny, qext)
        do j=j0,ny,jskip
          jstr=max(1,j) ; jend=min(j+jsize-1,ny)
          call set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
        enddo
        nmsk_new=0
        call shortlist(nmsk,ijmsk, nmsk_new)
        nmsk=nmsk_new
C$OMP CRITICAL(etch_cr_rgn)
        if (trd_count == 0) mskd_pts=0
        trd_count=trd_count+1
        mskd_pts=mskd_pts+nmsk
        if (trd_count == ntrds) then
          trd_count=0
#ifdef VERBOSE
          write(*,*) 'rmaining masked points =', mskd_pts
#endif
        endif
C$OMP END CRITICAL(etch_cr_rgn)
C$OMP BARRIER
      enddo  !<-- while

      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call strip_halo_tile(istr,iend,jstr,jend, nx,ny, qext,qfld)
      enddo
C$OMP BARRIER
      end

      subroutine copy_extend_tile(istr,iend,jstr,jend, nx,ny, mask,
     &                                             qsrc,qext, nmsk)

! Copy array "qsrc" into "qext" with has one row of ghost points all
! around, while setting land points to special values. The logic here
! is designed to work both ways: either there a non-trivial land mask
! array, or masked data is already set to some special value, while
! mask(i,j) is trivially set to all-water mask(:,:)=1 status and makes
! no effect. The secondary goal is to determine the total number of
! special-value points encountered by this thread so an appropiate
! sized arrays can be allocated to hold list of indices.

      implicit none
      integer istr,iend,jstr,jend, nx,ny, nmsk, i,j
      integer(kind=2) mask(nx,ny)
      real(kind=4) qsrc(nx,ny), qext(0:nx+1,0:ny+1)
!>    write(*,*) 'enter copy_extend_tile'
      do j=jstr,jend
        do i=istr,iend
          if (mask(i,j) > 0 .and. abs(qsrc(i,j)) < abs(spv)) then
            qext(i,j)=qsrc(i,j)
          else
            qext(i,j)=spv_set ; nmsk=nmsk+1
          endif
        enddo
      enddo
      end

      subroutine strip_halo_tile(istr,iend,jstr,jend, nx,ny, qext,qsrc)
      implicit none
      integer istr,iend,jstr,jend, nx,ny, i,j
      real(kind=4) qext(0:nx+1,0:ny+1), qsrc(nx,ny)
!>    write(*,*) 'enter strip_halo_tile'
      do j=jstr,jend
        do i=istr,iend
          qsrc(i,j)=qext(i,j)
        enddo
      enddo
      end

      subroutine set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
      implicit none
      integer istr,iend,jstr,jend, nx,ny, i,j
      real(kind=4) qext(0:nx+1,0:ny+1)
!>    write(*,*) 'enter set_qext_bc_tile'
      if (WESTERN_EDGE) then
        do j=jstr,jend
          qext(istr-1,j)=qext(istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          qext(iend+1,j)=qext(iend,j)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          qext(i,jstr-1)=qext(i,jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          qext(i,jend+1)=qext(i,jend)
        enddo
      endif
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        qext(istr-1,jstr-1)=qext(istr,jstr)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        qext(istr-1,jend+1)=qext(istr,jend)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        qext(iend+1,jstr-1)=qext(iend,jstr)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        qext(iend+1,jend+1)=qext(iend,jend)
      endif
      end

      subroutine init_ijmsk_tile(istr,iend,jstr,jend, nx,ny, src,
     &                                      max_pts, ijmsk, nmsk)
      implicit none
      integer istr,iend,jstr,jend, nx,ny,  max_pts, nmsk, i,j
      real(kind=4) src(0:nx+1,0:ny+1)
      integer(kind=2) ijmsk(2,max_pts)
!>    write(*,*) 'enter init_ijmsk_tile'
      do j=jstr,jend                       ! Form list of indices of
        do i=istr,iend                     ! points with special value.
          if (src(i,j) < spv) then
            nmsk=nmsk+1
            ijmsk(1,nmsk)=i
            ijmsk(2,nmsk)=j
          endif
        enddo
      enddo
      end

      subroutine select_coast_pts(nx,ny, src, nmsk,ijmsk, ncst,nijmsk)
      implicit none
      integer nx,ny, nmsk, ncst, i,j,n     ! Take list of previously
      real(kind=4) src(0:nx+1,0:ny+1)      ! identified masked points
      integer(kind=2) ijmsk(2,nmsk)        ! and select the ones among
      integer(kind=4) nijmsk(nmsk)         ! them which have at least
!>    write(*,*) 'enter select_coast_pts'  ! one water neighbor.
      do n=1,nmsk
        i=ijmsk(1,n) ; j=ijmsk(2,n)
        if (src(i,j) < spv) then
          if (src(i+1,j) > spv .or. src(i,j+1) > spv .or.
     &        src(i-1,j) > spv .or. src(i,j-1) > spv ) then
            ncst=ncst+1 ; nijmsk(ncst)=n
          endif
        endif
      enddo
      end

      subroutine comp_patch_pts(nx,ny, src, nmsk,ijmsk, ncst,nijmsk,
     &                                          npths, ijptch, ptch)
      implicit none
      integer nx,ny, nmsk, ncst, npths, i,j,n
      real(kind=4) src(0:nx+1,0:ny+1),   ptch(ncst), wgt,vlu
      integer(kind=2) ijmsk(2,nmsk), ijptch(2,ncst)
      integer(kind=4)  nijmsk(nmsk)

      real(kind=4), parameter :: grad=1./3.,  corn=0.707106781186547,
     &                          corngrad=0.5*corn*grad, threshold=2.4

!>    write(*,*) 'enter comp_patch_pts, ncst =',ncst

      do n=1,ncst                         ! check surrounding points:
        i=ijmsk(1,nijmsk(n)) ; wgt=0.     ! counterclockwise direction
        j=ijmsk(2,nijmsk(n)) ; vlu=0.     ! starting from the east.

        if (src(i+1,j) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i+1,j)
          if (i < nx) then
            if (src(i+2,j) > spv) vlu=vlu+grad*(src(i+1,j)-src(i+2,j))
          endif
        endif
#define DIAGONAL
#ifdef DIAGONAL
        if (src(i+1,j+1) > spv) then
          wgt=wgt+corn ;  vlu=vlu + corn*src(i+1,j+1)
          if (i < nx) then
            if (src(i+2,j+1) > spv) vlu=vlu + corngrad*( src(i+1,j+1)
     &                                                  -src(i+2,j+1))
          endif
          if (j < ny) then
            if (src(i+1,j+2) > spv) vlu=vlu  + corngrad*(src(i+1,j+1)
     &                                                  -src(i+1,j+2))
          endif
        endif
#endif

        if (src(i,j+1) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i,j+1)
          if (j < ny) then
            if (src(i,j+2) > spv) vlu=vlu+grad*(src(i,j+1)-src(i,j+2))
          endif
        endif

#ifdef DIAGONAL
        if (src(i-1,j+1) > spv) then
          wgt=wgt+corn ;   vlu=vlu + corn*src(i-1,j+1)
          if (j < ny) then
            if (src(i-1,j+2) > spv) vlu=vlu  + corngrad*(src(i-1,j+1)
     &                                                  -src(i-1,j+2))
          endif
          if (i > 1) then
            if (src(i-2,j+1) > spv) vlu=vlu  + corngrad*(src(i-1,j+1)
     &                                                  -src(i-2,j+1))
          endif
        endif
#endif

        if (src(i-1,j) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i-1,j)
          if (i > 1) then
            if (src(i-2,j) > spv) vlu=vlu+grad*(src(i-1,j)-src(i-2,j))
          endif
        endif

#ifdef DIAGONAL
        if (src(i-1,j-1) > spv) then
          wgt=wgt+corn ;   vlu=vlu + corn*src(i-1,j-1)
          if (i > 1) then
            if (src(i-2,j-1) > spv) vlu=vlu  + corngrad*(src(i-1,j-1)
     &                                                  -src(i-2,j-1))
          endif
          if (j > 1) then
            if (src(i-1,j-2) > spv) vlu=vlu  + corngrad*(src(i-1,j-1)
     &                                                  -src(i-1,j-2))
          endif
        endif
#endif

        if (src(i,j-1) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i,j-1)
          if (j > 1) then
            if (src(i,j-2) > spv) vlu=vlu+grad*(src(i,j-1)-src(i,j-2))
          endif
        endif

#ifdef DIAGONAL
        if (src(i+1,j-1) > spv) then
          wgt=wgt+corn ;  vlu=vlu + corn*src(i+1,j-1)
          if (j > 1) then
            if (src(i+1,j-2) > spv) vlu=vlu  + corngrad*(src(i+1,j-1)
     &                                                  -src(i+1,j-2))
          endif
          if (i < nx) then
            if (src(i+2,j-1) > spv) vlu=vlu  + corngrad*(src(i+1,j-1)
     &                                                  -src(i+2,j-1))
          endif
        endif
#endif
        if (wgt > threshold) then
          npths=npths+1              ! At the end set "ijmsk" i-index
          ptch(npths)=vlu/wgt        ! to zero to signal that the point
          ijptch(1,npths)=i          ! is no longer a special value.
          ijptch(2,npths)=j
          ijmsk(1,nijmsk(n))=0
        endif
      enddo
      end

      subroutine apply_patch_pts(npths, ijptch,ptch,  nx,ny, src)
      implicit none
      integer npths, nx,ny, i,j,n
      integer(kind=2) ijptch(2,npths)
      real(kind=4) ptch(npths), src(0:nx+1,0:ny+1)
!>    write(*,*) 'enter apply_patch_pts'
      do n=1,npths
        i=ijptch(1,n) ; j=ijptch(2,n)
        src(i,j)=ptch(n)
      enddo
      end


      subroutine shortlist(nmsk,ijmsk, nmsk_new)
      implicit none
      integer nmsk,nmsk_new,n
      integer(kind=2) ijmsk(2,nmsk)
      do n=1,nmsk
        if (ijmsk(1,n) > 0) then
          nmsk_new=nmsk_new+1
          ijmsk(1,nmsk_new)=ijmsk(1,n)
          ijmsk(2,nmsk_new)=ijmsk(2,n)
        endif
      enddo
!>    write(*,*) 'shortlist, nmsk=', nmsk, ' -->', nmsk_new
      end
