      program r2r_match_topo

! Modifies child-grid bottom topography in the vicinity of forced
! open boundaries to exactly match topography interpolated from the
! parent grid right at the boundary, while at the same time, time
! making smooth transition into the interior of child-grid domain
! where the topography is expected to be less smooth because of
! finer resolution.

! The method is essentially generating a 2D weighting function,
! wgt=wgt(i,j), 0 <= wgt <= 1, hence
!
!             h_matched = (1-wgt) * h + wgt * h_parent
!
! where wgt=1 at open boundaries, and wgt=0 in the interior beyond
! merging zone, with smooth transition within - sounds quite simple,
! however is should be noted that merging should be avoided in the
! areas closed to the edge of the child domain, but blocked by the
! land mask. This in its turn leads to a rather complicated logical
! "etching" to determine the exact coverage of merging zone.
! Interpolation of parent-grid topography is by bi-cubic spline
! using Hermite polynomials.

! Created/maintained by Alexander Shchepetkin, old_galaxy@yahoo.com

c--#define LEGACY_METHOD
#define SMOOTH_JUNCTIONS

      implicit none
      character(len=4) merging_flags
      character(len=64) prntgrd, chldgrd, str
      real(kind=8), allocatable, dimension(:,:) :: hp,xp,yp, hc,xc,yc,
     &                                xi,eta, wgt,hwg, srX,srY,sXY,sYX
      real(kind=4), allocatable, dimension(:,:) :: mwgt
      integer(kind=4), allocatable, dimension(:,:) :: ip,jp
      integer(kind=2), allocatable, dimension(:,:) :: mask,mgz,ms1,ms2
      integer nargs, wdth, net_alloc_mem, nx,ny, ncx,ncy, ncsrc,nctarg,
     &        ncid, ierr, ipvar,jpvar, xivar,etavar, hpvar,hvar, dhvar,
     &        wgtvar,mwgvar, r2dgrd(2), i,j,k, iter, lprnt,lchld,
     &                                      lstr, lenstr, read_str_int
      real(kind=8) xmin,xmax, xcmin,xcmax, cff
      logical OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH, rename_hvar

       integer rad,rad2, ict,jct

#include "netcdf.inc"

      net_alloc_mem=0 !<-- allocated memory expressed in 4-byte numbers
      OBC_WEST=.false. ; OBC_SOUTH=.false. ; merging_flags='    '
      OBC_EAST=.false. ; OBC_NORTH=.false. ; rename_hvar=.false.

! Check whether all the arguments are specified correctly and both
! netCDF files can be opened (at least in read-only mode); diagnose
! errors and quite if an errors occurs; write help page and quit if
! the program is called with insufficient number of arguments.

      nargs=iargc()
      if (nargs == 4) then
        call getarg(1,str) ; lstr=lenstr(str)
        do i=1,lstr
          if (str(i:i) == 'W'  .or.  str(i:i) == 'w') then
            OBC_WEST=.true.
          elseif (str(i:i) == 'E' .or. str(i:i) == 'e') then
            OBC_EAST=.true.
          elseif (str(i:i) == 'S' .or. str(i:i) == 's') then
            OBC_SOUTH=.true.
          elseif (str(i:i) == 'N' .or. str(i:i) == 'n') then
            OBC_NORTH=.true.
          else
            write(*,'(/1x,7A/12x,A)')       '### ERROR: Illegal ',
     &       'argument: letter ''',  str(i:i),  ''' should not ',
     &       'be present in the first argument ''', str(1:lstr),
     &       '''.',   'Only "W", "E", "S",  and "N" are allowed.'
            stop
          endif
        enddo
        k=0
        if (OBC_WEST) then                    ! prepare signature
          k=k+1 ; merging_flags(k:k)='W'      ! string to be saved
        endif                                 ! as an attribute in
        if (OBC_EAST) then                    ! netCDF file
         k=k+1 ; merging_flags(k:k)='E'
        endif
        if (OBC_SOUTH) then
          k=k+1 ; merging_flags(k:k)='S'
        endif
        if (OBC_NORTH) then
          k=k+1 ; merging_flags(k:k)='N'
        endif
        call getarg(2,str) ; lstr=lenstr(str)
        ierr=read_str_int(str, wdth)
        if (ierr == 0) then
          call getarg(3,prntgrd) ; lprnt=lenstr(prntgrd)
          ierr=nf_open(prntgrd, nf_nowrite, ncsrc)
          if (ierr == nf_noerr) then
            call getarg(4,chldgrd) ; lchld=lenstr(chldgrd)
            ierr=nf_open(chldgrd, nf_nowrite, nctarg)
            if (ierr == nf_noerr) then
              write(*,'(/1x,A,4(1x,L1),4x,A,I4/1x,3A,4x,3A/)')
     &              'OBC(W,E,S,N) =', OBC_WEST,OBC_EAST,OBC_SOUTH,
     &               OBC_NORTH,   'width of merging zone =', wdth,
     &                  'Parent grid =''', prntgrd(1:lprnt), '''',
     &                   'Child grid =''', chldgrd(1:lchld), ''''
            else
              write(*,'(/1x,4A/12x,A)')  '### ERROR: Cannot open ',
     &                       'file ''', chldgrd(1:lenstr(chldgrd)),
     &                                    '''.', nf_strerror(ierr)
              stop
            endif
          else
            write(*,'(/1x,4A/12x,A)')  '### ERROR: Cannot open ''',
     &                       'file ''', prntgrd(1:lenstr(prntgrd)),
     &                                    '''.', nf_strerror(ierr)
            stop
          endif
        else
          write(*,'(/1x,4A/)')   '### ERROR: The second argument ',
     &                 'should be an integer number instead of ''',
     &                                         str(1:lstr),  '''.'
          stop
        endif
      else
        write(*,'(3(/1x,2A)/1x,A/,/8x,A)')      'This command-line ',
     &   'operator modifies child-grid topography in the vicinity of',
     &        'its open boundaries to ensure exact matching between ',
     &     'the parent and the child', 'at the boundary followed by ',
     &            'smooth transition into the interior child domain.',
     &     'It should be used as',
     &       'r2r_match_topo WESN width parent_grid.nc child_grid.nc'
        write(*,'(4(/1x,2A))')     'where argument WESN should have ',
     &                 'up to 4 letters (W, E, S, N only, upper- or',
     &           'lowercase letters for Western, Eastern, Southern, ',
     &           'Northern sides) indicating','which ones have open ',
     &     'boundaries; argument "width" (an integer number) is the',
     &       'desired width of merging zone expressed as the number ',
     &                                        'of child-grid points.'
        write(*,'(9(/1x,2A)/1x,A/)') 'NOTE: The original child-grid ',
     &                  'topography stored as variable "h" in netCDF',
     &   'file "child_grid.nc" will be renamed into "hsmth" without ',
     &   'any modification',  '(hence backed up), while the matched ',
     & 'topography will be placed into a newly',  'created variable ',
     &  'also named "h" -- as expected by ROMS model.  On the other',
     &      'hand if "hsmth" already exists in the file, then it is ',
     &  'presumed that it is',  'the original version of topography ',
     & '(smoothed, but not matched), so in this','case "hsmth" will ',
     &   'be the one used as input, while no backup of "h" needs to',
     &         'be made.  This mechanism ensures that applying this ',
     &       'operator to the same', 'file again and again does not ',
     &               'modify topography which already been modified',
     &        'but rather results in consistent outcome every time.'
        stop
      endif

! Allocate arrays and read coordinates first for the target grid,
! then for the source.  This leads to a more optimal memory use as
! xp,yp arrays for the source grid can be deallocated after computing
! ip,jp-indices and fractional offsets xi,eta, the xp,yp are allocated
! at the end to be deallocated to free memory for arrays associated
! with spline interpolation.

      call roms_find_dims(nctarg, chldgrd, i,j,k)
      ncx=i+2 ;  ncy=j+2
      allocate( xi(ncx,ncy),eta(ncx,ncy), xc(ncx,ncy),yc(ncx,ncy),
     &                       wgt(ncx,ncy),hwg(ncx,ncy), stat=ierr )
      if (ierr == 0) allocate( mwgt(ncx,ncy), stat=ierr )
      if (ierr == 0) allocate( ip(ncx,ncy),jp(ncx,ncy), stat=ierr )
      if (ierr == 0) allocate( mask(ncx,ncy), mgz(ncx,ncy),
     &                       ms1(ncx,ncy),ms2(ncx,ncy), stat=ierr )
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem +ncx*ncy*(6*2 + 1*1 + 2*1)
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 1.'
        stop
      endif
      call get_var_by_name_double(nctarg, 'lon_rho',  xc)
      call get_var_by_name_double(nctarg, 'lat_rho',  yc)
      call get_var_by_name_double(nctarg, 'mask_rho', xi)
      do j=1,ncy
        do i=1,ncx
          if (xi(i,j) > 0.5D0) then   ! temporarily use array "xi"
            mask(i,j)=1               ! to read mask from the file;
          else                        ! thereafter this program uses
            mask(i,j)=0               ! only integer(kind=2) version
          endif                       ! of mask
        enddo
      enddo
      ierr=nf_close(nctarg)

      xcmin=xc(1,1) ; xcmax=xc(1,1)
      do j=1,ncy
        do i=1,ncx
          if (xc(i,j) < xcmin) then
            xcmin=xc(i,j)
          elseif (xc(i,j) > xcmax) then
            xcmax=xc(i,j)
          endif
        enddo
      enddo
      write(*,*) 'Child grid xcmin,xcmax =', xcmin,xcmax

! The same for parent grid

      call roms_find_dims(ncsrc, prntgrd, i,j,k)
      nx=i+2 ; ny=j+2
      allocate(hp(nx,ny), xp(nx,ny), yp(nx,ny), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem+(3*2)*nx*ny
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 2.'
      stop
      endif
      call get_var_by_name_double(ncsrc, 'lon_rho', xp)
      call get_var_by_name_double(ncsrc, 'lat_rho', yp)
      call get_var_by_name_double(ncsrc, 'h',       hp)
      ierr=nf_close(ncsrc)

! Interpolate parent-grid topography onto child grid...

      xmin=xp(1,1) ; xmax=xp(1,1)
      do j=1,ny
        do i=1,nx
          if (xp(i,j) < xmin) then
            xmin=xp(i,j)
          elseif (xp(i,j) > xmax) then
            xmax=xp(i,j)
          endif
        enddo
      enddo
      write(*,*) 'Parent grid xmin,xmax =', xmin,xmax
      if (xmin > xcmax) then
        do j=1,ny
          do i=1,nx
            xp(i,j)=xp(i,j)-360.D0
          enddo
        enddo
        write(*,*) 'Adjusted to',xmin-360.D0,xmax-360.D0
      elseif (xmax < xcmin) then
        do j=1,ny
          do i=1,nx
            xp(i,j)=xp(i,j)+360.D0
          enddo
        enddo
        write(*,*) 'Adjusted to',xmin+360.D0, xmax+360.D0
      endif

C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp, xi,eta)
        call r2r_interp_init_thread(nx,ny, xp,yp, ncx,ncy, xc,yc,
     &                                              ip,jp, xi,eta)
C$OMP END PARALLEL

      call check_search_indices(nx,ny,xp,yp, ncx,ncy, xc,yc, ip,jp)
      call check_offsets(nx,ny,xp,yp, ncx,ncy,xc,yc, ip,jp, xi,eta)

      deallocate(xp,yp) ; net_alloc_mem=net_alloc_mem-4*nx*ny

      allocate(srX(nx,ny),srY(nx,ny), sXY(nx,ny),sYX(nx,ny), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + (4*2)*nx*ny
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 3.'
        stop
      endif

      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)

      deallocate(sYX) ; net_alloc_mem=net_alloc_mem-2*nx*ny
      allocate(hc(ncx,ncy), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + (1*2)*ncx*ncy
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 4.'
        stop
      endif

C$OMP PARALLEL SHARED(nx,ny, hp, srX,srY,sXY, ncx,ncy,
C$OMP&                               ip,jp, xi,eta, hc)
      call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY,
     &                      ncx,ncy, ip,jp, xi,eta, hc)
C$OMP END PARALLEL


! Start forming merging zone and weights: the two methods below differ
!------ --------------- -------- -------- only by handling land mask
! and produce identical results in the case of no land.  Either way,
! the outcome is "mgz" assuming values 0 or 1 to define merging zone
! (unmasked water points only) and "wgt" initialized as constant-slope
! function wgt=wdth at each open boundary water points from which it
! descends linearly into the interior of the domain by decreasing by
! 1 for each row of grid points until it vanishes to 0.  In the case
! of non-LEGACY_METHOD "wgt" also decreases by 1 every step when it
! goes around the corner of the coastline -- the minimal path from
! the costline to the given point is counted as the number of steps
! connecting unmasked water points by C-grid masking rules rather than
! along a straight line to the nearest boundary point, so the resultant
! weights and the extent of merging zone expected to be somewhat
! smaller in areas behind land.

#ifdef LEGACY_METHOD

! Initialize constant-slope function "wgt" and the initial version of
! merging zone "mgz" where wgt > 0 in a straightforward analytical way,
! then apply land mask which potentially may break up "mgz" into
! several unconnected areas some of which cannot be reached by water
! points from the open boundary, then exclude such points by applying
! algorithm similar to enforcing connectivity of unmasked water area
! by single connect.

      do j=1,ncy
        do i=1,ncx
          k=0
          if (OBC_WEST)  k=max(k, wdth-i+1)
          if (OBC_EAST)  k=max(k, wdth-ncx+i)
          if (OBC_SOUTH) k=max(k, wdth-j+1)
          if (OBC_NORTH) k=max(k, wdth-ncy+j)

          if (k > 0 .and. mask(i,j) > 0) then
            mgz(i,j)=1 ; wgt(i,j)=dble(k)
          else
            mgz(i,j)=0 ; wgt(i,j)=0.D0
          endif
        enddo
      enddo

C$OMP PARALLEL SHARED(ncx,ncy, mgz,ms1,ms2)
      call mrg_zone_cont_thread(ncx,ncy, mgz,ms1,ms2)
C$OMP END PARALLEL

#else
      do j=1,ncy                     ! Initialize "mgz" as a single
        do i=1,ncx                   ! row of perimeter points along
          mgz(i,j)=0                 ! the unmasked parts of open
        enddo                        ! boundaries, then "etch" inward,
      enddo                          ! while obeying the connectivity
      if (OBC_WEST) then             ! rules of land-water masking:
        do j=1,ncy                   ! the by the construction "mgz"
          mgz(1,j)=mask(1,j)         ! has the property that each of
        enddo                        ! its points can be reached from
      endif                          ! the open boundary by water.
      if (OBC_SOUTH) then
        do i=1,ncx
          mgz(i,1)=mask(i,1)
        enddo
      endif
      if (OBC_EAST) then
        do j=1,ncy
          mgz(ncx,j)=mask(ncx,j)
        enddo
      endif
      if (OBC_NORTH) then
        do i=1,ncx
          mgz(i,ncy)=mask(i,ncy)
        enddo
      endif

C$OMP PARALLEL SHARED(ncx,ncy, mask,ms1,ms2, wdth)
      call etch_mgz_weights_thread(ncx,ncy, mask,mgz,ms2, wdth)
C$OMP END PARALLEL

      do j=1,ncy
        do i=1,ncx
          if (mgz(i,j) > 0) then
            wgt(i,j)=dble(mgz(i,j)) ; mgz(i,j)=1
          else
            wgt(i,j)=0.D0
          endif
        enddo
      enddo
#endif /* LEGACY_METHOD */

! The following part makes round junctions between two adjacent open
! boundaries by extending merging zone inward. This is useful to avoid
! steep gradient near inner corner of merging zone when weighting
! finction "wgt" is subjected to Laplacian smoothing.  Note that all
! loops below are non-reversible as they process points in row-by-row
! manner (both in i and j directions) with checking that the previous
! row is set to 1.

#ifdef SMOOTH_JUNCTIONS
      rad=(wdth+1)/3 ; rad2=rad**2
      if (OBC_NORTH .and. OBC_WEST) then
        ict=wdth+rad+1 ; jct=ncy-wdth-rad  !<-- center of the circle
        do j=jct+rad,jct,-1
          do i=ict-rad,ict,+1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i-1,j) > 0 .and. mgz(i,j+1) > 0
     &                         .and. mgz(i-1,j+1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_WEST .and. OBC_SOUTH) then
        ict=wdth+rad+1 ; jct=wdth+rad+1
        do j=jct-rad,jct,+1
          do i=ict-rad,ict,+1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i-1,j) > 0 .and. mgz(i,j-1) > 0
     &                         .and. mgz(i-1,j-1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_SOUTH .and. OBC_EAST) then
        ict=ncx-wdth-rad ; jct=wdth+rad+1  !<-- center of the circle
        do j=jct-rad,jct,+1
          do i=ict+rad,ict,-1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i+1,j) > 0 .and. mgz(i,j-1) > 0
     &                         .and. mgz(i+1,j-1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_EAST .and. OBC_NORTH) then
        ict=ncx-wdth-rad ; jct=ncy-wdth-rad
        do j=jct+rad,jct,-1
          do i=ict+rad,ict,-1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i+1,j) > 0 .and. mgz(i,j+1) > 0
     &                         .and. mgz(i+1,j+1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
#endif

! Etch the merging zone into the land mask, in such a way that
! any possibility of reaching water behind mask is excluded.

C$OMP PARALLEL SHARED(ncx,ncy, mask, mgz,ms1,ms2)
      call etch_mgz_into_land_thread(ncx,ncy, mask, mgz,ms1,ms2)
C$OMP END PARALLEL

! At this stage "ms1 > 0" defines the area to which the merging zone
! should be allowed to grow. This area consists of all water points
! of the merging zone which can be reached by water starting from
! the boundary row of non-masked points, and the adjacent land area
! to which to merging zone can be expanded to ensure smoothness of
! the merged topography, but without interfering with water interior.
! The next stage is to etch the actual merging function into the
! land-masked area.

      do j=1,ncy                            ! Set "mgz" to be merging
        do i=1,ncx                          ! zone mask (positive "ms1"
          if (ms1(i,j) > 0) then            ! area), which will remain
            if (mask(i,j) > 0) then         ! unchanged from now on,
              ms1(i,j)=int(wgt(i,j)+0.5D0)  ! then initialize both
            else                            ! "ms1" and "ms2" to the
              ms1(i,j)=0                    ! constant-slope merging
            endif                           ! function (which later
            mgz(i,j)=1                      ! be used as the argument
          else                              ! for the actual merging
            ms1(i,j)=0                      ! function), using integer
            mgz(i,j)=0                      ! numbers from the range
          endif                             ! from 0 to "wdth".
          ms2(i,j)=ms1(i,j)
        enddo
      enddo

C$OMP PARALLEL SHARED(wdth, ncx,ncy, mgz,ms1,ms2, OBC_WEST,OBC_EAST,
C$OMP&                                           OBC_SOUTH,OBC_NORTH)
      call etch_weights_into_land_thread(wdth, ncx,ncy, mgz,ms1,ms2,
     &                         OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
C$OMP END PARALLEL

! At this moment if there would be no land "ms1" is constant-slope
! function: the same maximum value at the outermost row of points and
! descending linearly to zero within the band "wdth" into the domain -
! so smoothing iterations below (subject to Dirichlet B.C) would not
! actually cause any change to it with the exception of (i) in the
! vicinity of corners between adjacent open boundaries (resulting in
! sharp bend along the 45-degree row points where two side slopes
! come together), and (ii) irregular shapes caused by "etching" of
! "ms1" into the land.

      do j=1,ncy
        do i=1,ncx
          wgt(i,j)=dble(ms1(i,j)) ; hwg(i,j)=wgt(i,j)
        enddo
      enddo
      iter=4*wdth
C$OMP PARALLEL SHARED(iter, ncx,ncy, mgz,wgt,hwg)
      call smooth_wgt_thread(iter, ncx,ncy, mgz,wgt,hwg)
C$OMP END PARALLEL

! Read in child grid topography:  To avoid repeated modification of
!----- -- ----- ---- -----------  the same topography stored in netCDF
! file as variable "h" if this program is applied again and again, the
! algorithm below is designed to rename "h" --> "hsmth" in the netCDF
! file without modifying it (hence making a back up copy of the original
! topography), then create a brand new variable also named "h" for the
! merged version of topography.  Therefore, in anticipation that the
! file may already be processed by this program, check first if "hsmth"
! already exists:
!         if yes, read it in and use it as input;
!         if not, read "h", rename, and create new "h".
! As the result, if applied again and again to the same file with
! the same input arguments, the program creates exactly the same
! merged "h" (not successively modified).

      ierr=nf_open(chldgrd, nf_nowrite, ncid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_varid(ncid, 'hsmth', hvar)
        if (ierr == nf_noerr) then
          write(*,'(1x,4A/8x,A)') 'Found ''hsmth'' - smoothed but ',
     &         'not yet merged topo in ''', chldgrd(1:lchld), '''.',
     &                                  'Using ''hsmth'' as input.'
          rename_hvar=.false.
        else
          ierr=nf_inq_varid(ncid, 'h', hvar)
          if (ierr == nf_noerr) then
            write(*,*) 'No variable ''hsmth'' is present in ''',
     &             chldgrd(1:lchld), '''. Using ''h'' as input.'
            rename_hvar=.true.
          else
            write(*,*) '### ERROR: Neither ''hsmth'', no ''h'' is ',
     &                     'present in ''', chldgrd(1:lchld), '''.'
          endif
        endif
        if (ierr == nf_noerr) then
          ierr=nf_get_var_double(ncid, hvar, hwg)
          if (ierr == nf_noerr) then
            ierr=nf_close(ncid)
            if (ierr == nf_noerr) then
              write(*,'(8x,A)') 'read child grid topography'
            else
              write(*,*)   '### ERROR: Cannot close ''',
     &                          chldgrd(1:lchld), '''.'
            endif
          else
            write(*,*) '### ERROR: Cannot read topography from ''',
     &                                     chldgrd(1:lchld), '''.'
          endif
        endif
      else
        write(*,*) '### ERROR: Cannot open netCDF file ''',
     &                chldgrd(1:lchld), ''' for reading.'
      endif
      if (ierr /= nf_noerr) stop

! Convert the originally constant-slope --> etched --> smoothed "wgt"
! into merging function, then perform the actual merging of topography.
! Note that the incoming "wgt" is within the range of 0 <= wgt <= wdth
! nominally changing by 1 from one point to the next (nominally because
! masking-etching-smoothing may modify this). A slightly different
! rescaling,
!               -1/(2*(wdth-1)) <= wgt <= (wdth-1/2)(wdth-1)
!
! subject to limiting 0 <= wgt <= 1 is used below with the rationale
! to take out discretely jagged inner border of merging zone (obtained
! by etching), and to make cos() function a bit closer to 1 at the
! second row of points near the open boundary.

      write(*,*) 'Merging topography...'
      cff=1.D0/dble(wdth-1)
      do j=1,ncy
        do i=1,ncx
          if (mgz(i,j) > 0) then
            wgt(i,j)=min(1.D0, cff*(wgt(i,j)-0.5D0))
            if (wgt(i,j) < 0.D0) then
              wgt(i,j)=0.D0 ; mgz(i,j)=0
            endif
            wgt(i,j)=0.5D0-0.5D0*cos(3.14159265358979323D0*wgt(i,j))
            hwg(i,j)=hwg(i,j)*(1.D0-wgt(i,j))+hc(i,j)*wgt(i,j)
          endif

          if (mgz(i,j) > 0) then
            if (mask(i,j) > 0) then
              mwgt(i,j)=1.
            else                          ! Array "mwgt" is merely
              mwgt(i,j)=0.8               ! for 4-color illustration
            endif                         ! of the layout of merging
          else                            ! zone relatively to land
            wgt(i,j)=0.D0                 ! mask:
            if (mask(i,j) > 0) then       !
              mwgt(i,j)=0.                !    mwgt=1 or 0.8 merge
            else                          !    mwgt=0 or 0.2 intact
              mwgt(i,j)=0.2               !
            endif                         !    mwgt=0  or  1 water
          endif                           !    mwgt= 0.8 or 0.2 land
        enddo
      enddo

! Write back into the file, also save all the input arguments given
! to the program, so, in principle, all what is going on here can be
! reproduced exactly relying solely on info stored in netCDF file.

      write(*,*) 'Writing merged topography into child grid file ''',
     &                                     chldgrd(1:lchld), '''...'
      ierr=nf_open(chldgrd, nf_write, ncid)
      if (ierr == nf_noerr) then
        ierr=nf_redef(ncid)
        if (ierr == nf_noerr) then
          ierr=nf_inq_dimid(ncid, 'xi_rho', r2dgrd(1))
          if (ierr == nf_noerr) then
            ierr=nf_inq_dimid(ncid, 'eta_rho', r2dgrd(2))
            if (ierr == nf_noerr) then
              ierr=nf_inq_varid(ncid, 'h', hvar)
              if (ierr == nf_noerr) then
                if (rename_hvar) then
                  ierr=nf_rename_var(ncid, hvar, 'hsmth')
                  if (ierr == nf_noerr) then
                    write(*,'(8x,A)') 'renamed "h" --> "hsmth"'
                    ierr=nf_def_var(ncid, 'h', nf_double, 2,
     &                                           r2dgrd, hvar)
                    if (ierr == nf_noerr) then
                      write(*,'(8x,A)') 'created new variable "h"'
                    else
                      write(*,*) '### ERROR: Cannot create netCDF ',
     &                           'variable "h".', nf_strerror(ierr)
                    endif
                  else
                    write(*,*) '### ERROR: Cannot rename "h" --> ',
     &                                                  '"hsmth".'
                  endif
                endif  !<-- rename_hvar
                ierr=nf_inq_varid(ncid, 'wgt', wgtvar)
                if (ierr == nf_noerr) then
                   write(*,'(8x,A)') 'variable "wgt" aready exists'
                else
                  ierr=nf_def_var(ncid, 'wgt', nf_double, 2, r2dgrd,
     &                                                       wgtvar)
                  if (ierr == nf_noerr) then
                    write(*,'(8x,A)') 'created new variable "wgt"'
                  else
                    write(*,*) '### ERROR: Cannot create netCDF ',
     &                         'variable "wgt".', nf_strerror(ierr)
                  endif
                endif
              else
                write(*,*) '### ERROR: Cannot find variable "h" in ''',
     &                                         chldgrd(1:lchld), '''.'
              endif
            else
              write(*,*) '### ERROR: Cannot find dimension "eta_rho" ',
     &             'in ''', chldgrd(1:lchld), '''.', nf_strerror(ierr)
            endif
          else
            write(*,*)    '### ERROR: Cannot find dimension "xi_rho" ',
     &             'in ''', chldgrd(1:lchld), '''.', nf_strerror(ierr)
          endif


          if (ierr == nf_noerr) then
            ierr=nf_put_att_text(ncid, hvar, 'long_name', 60,
     &  'smooth topography merged with parent grid at open boundaries')
            if (ierr == nf_noerr) then
              ierr=nf_put_att_text(ncid, hvar, 'parent_grid', lprnt,
     &                                            prntgrd(1:lprnt))
              if (ierr == nf_noerr) then
                lstr=lenstr(merging_flags)
                ierr=nf_put_att_text(ncid, hvar, 'merging_flags', lstr,
     &                                          merging_flags(1:lstr))
                if (ierr == nf_noerr) then
                  ierr=nf_put_att_int(ncid, hvar, 'merging_width',
     &                                            nf_int, 1,wdth)
                  if (ierr == nf_noerr) then
                    ierr=nf_put_att_text(ncid, wgtvar, 'long_name', 42,
     &                   'parent-to-child topography merging weights')
                    if (ierr == nf_noerr) then
                      ierr=nf_enddef(ncid)
                      if (ierr == nf_noerr) then
                        ierr=nf_put_var_double(ncid, hvar, hwg)
                        if (ierr == nf_noerr) then
                          write(*,'(8x,A)') 'wrote merged topography'
                          ierr=nf_put_var_double(ncid, wgtvar, wgt)
                          if (ierr == nf_noerr) then
                            write(*,'(8x,A)') 'wrote merging weights'
                            ierr=nf_close(ncid)
                            if (ierr /= nf_noerr) then
                              write(*,*)    '### ERROR: Cannot close ',
     &                       'netCDF file ''', chldgrd(1:lchld), '''.'
                            endif
                          else
                            write(*,*) '### ERROR: Cannot write "wgt"',
     &                              ' into ''', chldgrd(1:lchld), '''.'
                          endif
                        else
                          write(*,*)    '### ERROR: Cannot write "h"',
     &                             ' into ''', chldgrd(1:lchld), '''.'
                        endif
                      else
                        write(*,*) '### ERROR: Cannot switch netCDF ',
     &                                    'file ''', chldgrd(1:lchld),
     &                                          ''' into input mode.'
                      endif
                    else
                      write(*,*)  '### ERROR: Cannot put attribute ',
     &                    '''long_name'' for merging weights "wgt".'
                    endif
                  else
                    write(*,*)  '### ERROR: Cannot put attribute ',
     &                                        '''merging_width''.'
                  endif
                else
                  write(*,*)  '### ERROR: Cannot put attribute ',
     &                                      '''merging_flags''.'
                endif
              else
                write(*,*)  '### ERROR: Cannot put attribute ',
     &                                      '''parent_grid''.'
              endif
            else
              write(*,*)  '### ERROR: Cannot put attribute ',
     &                    '''long name'' for topography "h".'
            endif
          endif
        else
          write(*,*)  '### ERROR: Cannot switch netCDF file ''',
     &           chldgrd(1:lchld), ''' into redefinition mode.'
        endif
      else
        write(*,*)   '### ERROR: Cannot open netCDF file ''',
     &                   chldgrd(1:lchld), ''' for writing.'
      endif
      if (ierr /= nf_noerr) stop

! Create diagnostic file: the following part - all the way to the very
!------- ---------- ----- end is merely to illustrate what this program
! is doing, and does not participate in any computation by ROMS model.

      write(*,*) 'Creating diagnostic file "r2r_match_topo_diag.nc"...'

      ierr=nf_create('r2r_match_topo_diag.nc', nf_64bit_offset, ncid)
      ierr=nf_def_dim(ncid, 'xi_rho',  ncx, r2dgrd(1))
      ierr=nf_def_dim(ncid, 'eta_rho', ncy, r2dgrd(2))

      ierr=nf_def_var(ncid, 'ip', nf_int, 2, r2dgrd, ipvar)
      ierr=nf_put_att_text(ncid, ipvar, 'long_name', 37,
     &                    'parent-to-child interpolation i-index')
      ierr=nf_put_att_int(ncid, ipvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'jp', nf_int, 2, r2dgrd, jpvar)
      ierr=nf_put_att_text(ncid, jpvar, 'long_name', 37,
     &                    'parent-to-child interpolation j-index')
      ierr=nf_put_att_int(ncid, jpvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'xi',  nf_double, 2, r2dgrd,  xivar)
      ierr=nf_put_att_text(ncid, xivar, 'long_name', 38,
     &                  'parent-to-child interpolation offset X')
      ierr=nf_put_att_double(ncid, xivar, '_FillValue', nf_double,
     &                                                  1, -1.D0)

      ierr=nf_def_var(ncid, 'eta', nf_double, 2, r2dgrd, etavar)
      ierr=nf_put_att_text(ncid, etavar, 'long_name', 38,
     &                  'parent-to-child interpolation offset Y')
      ierr=nf_put_att_double(ncid, etavar,'_FillValue', nf_double,
     &                                                  1, -1.D0)

      ierr=nf_def_var(ncid, 'hprnt', nf_double, 2, r2dgrd, hpvar)
      ierr=nf_put_att_text(ncid, hpvar, 'long_name', 35,
     &                     'interpolated parent grid topography')

      ierr=nf_def_var(ncid, 'h', nf_double, 2, r2dgrd, hvar)
      ierr=nf_put_att_text(ncid, hvar, 'long_name', 24,
     &                           'merged bottom topography')
      ierr=nf_put_att_text(ncid, hvar, 'parent_grid', lprnt,
     &                                     prntgrd(1:lprnt))
      lstr=lenstr(merging_flags)
      ierr=nf_put_att_text(ncid, hvar, 'merging_flags', lstr,
     &                                  merging_flags(1:lstr))
      ierr=nf_put_att_int(ncid,  hvar, 'merging_width', nf_int,
     &                                                 1,wdth)

      ierr=nf_def_var(ncid, 'dh', nf_double, 2, r2dgrd, dhvar)
      ierr=nf_put_att_text(ncid, dhvar, 'long_name', 37,
     &                  'merged child - parent topo difference')

      ierr=nf_def_var(ncid, 'wgt', nf_double, 2, r2dgrd, wgtvar)
      ierr=nf_put_att_text(ncid, wgtvar, 'long_name', 39,
     &                'parent-to-child merging weight function')

      ierr=nf_def_var(ncid, 'mwgt', nf_real, 2, r2dgrd,  mwgvar)
      ierr=nf_put_att_text(ncid, wgtvar, 'long_name', 28,
     &                          'parent-to-child merging zone')

      ierr=nf_enddef(ncid)
      ierr=nf_put_var_int(ncid, ipvar, ip)
      ierr=nf_put_var_int(ncid, jpvar, jp)
      ierr=nf_put_var_double(ncid, xivar,  xi)
      ierr=nf_put_var_double(ncid, etavar, eta)
      ierr=nf_put_var_double(ncid, hpvar,  hc)
      ierr=nf_put_var_double(ncid, hvar,   hwg)
      do j=1,ncy
        do i=1,ncx
          hwg(i,j)=hwg(i,j)-hc(i,j)
        enddo
      enddo
      ierr=nf_put_var_double(ncid, dhvar, hwg)
      ierr=nf_put_var_double(ncid, wgtvar, wgt)
      ierr=nf_put_var_real(ncid, mwgvar, mwgt)
      ierr=nf_close(ncid)
      end
