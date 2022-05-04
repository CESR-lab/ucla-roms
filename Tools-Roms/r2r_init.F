      program r2r_init

! Creates initial condition file for child-grid model by interpolating
! user-specified record of parent-grid solution.  Should be used as
!
!        r2r_init chld_grd.nc theta_s theta_b hc N prnt_grd.nc \
!                                              prnt_data.nc rec
!
! where theta_s theta_b hc N are S-coordinate parameters and number of
! vertical levels for the intended child-grid model;  chld_grd.nc and
! prnt_grd.nc are netCDF grid files for child and parent respectively;
! prnt_data.nc is file name of parent-grid history file containing data
! for zeta,u,v,T,S (mandatory variables for the initial conditions),
! "rec" is record number within this file (only one record will be
! read); upon completion the interpolated fields are written into
! netCDF file "roms_init.nc" (always named this way, specified by
! parameter roms_init below) containing just a single record and an
! auxiliary file "r2r_init_diag.nc" for diagnostic purposes for this
! program itself -- there is no use for it other than to ncview it
! and verify sanity of the algorithms below.

! Method: 2D horizontal bi-cubic spline interpolation via Hermite basis
! functions followed by vertical splines by spline inversion (z-levels
! of child grid are translated into continuous k-index space of parent
! grid in such a way spline interpolation of the resultant kp=kp(z_r)
! values back into z-space of child grid yields the original child-grid
! z_r=z_r(k) exactly. in doing so vertical spline derivatives of both
! z_r and the field to be interpolated are constructed in k-index
! space).  Prior to the horizontal interpolation data is extended to
! land-masked areas by etching algorithm.

! CPP-switch PARENT_GRID_SUBREGION make this program read only relevant
! subdomain within the parent grid as opposite to reading the whole
! data. The pros and cons for this are as follows: (1) all parent-grid
! arrays are allocated with smaller size; (2) computational savings
! mainly due to less data etching to land-masked areas (other savings
! as well, but not so dramatic); (3) as for the reading netCDF files
! it may or may not improve depending on the specific situation:
! reading sub-array from netCDF-3 file is often slower than reading
! the whole thing because it effectively reading many small records
! instead of reading in one single touch; reading compressed netCDF-4
! files may be faster because fewer blocks needs to be read (highly
! dependent on block structure of the source file).

! CPP-switch CONTOUR_CHILD_MASK reduces masking in the resultant
! file "roms_init.nc" to just one row of points along the coastline,
! while the internal pints inside land remain unmasked. This is useful
! to see what the result of etching procedure, but overall selecting
! this switch does not affect the usability of the file as the initial
! condition because ROMS code will mask all the land points on its own
! any way.

! CPP-switch WITH_REC_DIM makes all the variables the resultant file
! have time dimension which is unlimited dimension, but still having
! only one record written.  Both versions can be used as the initial
! to start ROMS model.

! Parameter setting "btm_slp" and "btm_trc" control type of boundary
! condition at bottom for the purpose of computing spline derivatives:
! btm_slp=-1 is no-slip, meaning that u,v=0 exactly at the bottom, i.e.
! half-grid-interval below the k=1 value; This applies to velocities
! only; for tracers the meaningful b.c. are natural on Neumann (0 and
! +1 respectively). Top b.c. is always natural.

! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com

#define PARENT_GRID_SUBREGION
#define CONTOUR_CHILD_MASK
c-#define WITH_REC_DIM
#define TIMING
c--#define VERBOSE

      use mod_io_size_acct
      implicit none
      character(len=16) VertCoordType, str
      character(len=80) chldgrd, prntgrd, prnt_data
      character(len=12), parameter :: roms_init='roms_init.nc'
      character(len=10), parameter :: time_var_name='ocean_time'
      integer, parameter :: ntrc=2
      character(len=4) tracer(ntrc)
      real(kind=8) theta_s,theta_b, hc,hcp, xcmin,xcmax, cff,cff1,cff2
      real(kind=8), allocatable, dimension(:) :: Cs_w,Cs_r,Csp_w,Csp_r,
     &                                           srX,srY,sXY,sYX, kprnt
      real(kind=8), allocatable, dimension(:,:) :: h,xc,yc,  hp,xp,yp,
     &      xi,eta, xpu,ypu,xiu,etau, xpv,ypv,xiv,etav, hprnt, csA,snA
      integer(kind=4), allocatable, dimension(:,:) :: ip,jp, ipu,jpu,
     &                                                       ipv,jpv
      integer(kind=2), allocatable, dimension(:,:) :: mskp, umsp, vmsp,
     &                                                mask, umask,vmask
      real(kind=4), allocatable, dimension(:) :: wrk1,wrk2,wrk3,wrmx

      integer net_alloc_mem, nargs, rec, ncpgrd,nccgrd, ncsrc,nctarg,
     &        nx,ny,Np, ncx,ncy,N, i,j,k, itrc,isrc,itrg, size, varid,
     &        tvar_in,tvar_out, natts, ierr, lstr, lpgrd,lprnt,lcgrd
      integer read_str_int, read_str_dble, lenstr !<-- functions
      integer, parameter :: btm_slp=-1 ! no-slip=-1; +1 free; 0 natural
     &                    , btm_trc=+1 ! +1 Neumann; 0 natural;
#ifdef PARENT_GRID_SUBREGION
      integer iwest,jsouth, imin,imax,jmin,jmax
      integer, parameter :: margn=8
#else
      integer, parameter :: iwest=1, jsouth=1
#endif
#ifdef TIMING
      real(kind=4) tstart,tend
      integer nclk, iclk(2), clk_rate, clk_max
      integer(kind=8) :: inc_clk, net_clk=0
#endif
#include "phys_const.h"
#include "netcdf.inc"
#ifdef TIMING
      call cpu_time(tstart)
      nclk=1
      call system_clock(iclk(nclk), clk_rate, clk_max)
# ifdef VERBOSE
      write(*,*) 'clk_rate=', clk_rate, ' clk_max =', clk_max
# endif
#endif
      tracer(1)='temp' ; tracer(2)='salt'
      nargs=iargc()
      if (nargs == 8) then
        call getarg(1,chldgrd) ; lcgrd=lenstr(chldgrd)
        ierr=nf_open(chldgrd(1:lcgrd), nf_nowrite, nccgrd)
        if (ierr == nf_noerr) then
          call getarg(2,str) ; ierr=read_str_dble(str, theta_s)
          if (ierr == 0) then
            call getarg(3,str) ; ierr=read_str_dble(str,theta_b)
            if (ierr == 0) then
              call getarg(4,str) ; ierr=read_str_dble(str, hc)
              if (ierr == 0) then
                call getarg(5,str) ; ierr=read_str_int(str, N)
                if (ierr == 0) then
                  call getarg(6,prntgrd) ; lpgrd=lenstr(prntgrd)
                  ierr=nf_open(prntgrd(1:lpgrd), nf_nowrite, ncpgrd)
                  if (ierr == nf_noerr) then
                    call getarg(7,prnt_data) ; lprnt=lenstr(prnt_data)
                    ierr=nf_open(prnt_data(1:lprnt), nf_nowrite, ncsrc)
                    if (ierr == nf_noerr) then
                      call getarg(8,str) ; ierr=read_str_int(str, rec)
                      if (ierr == 0) then
                        write(*,'(/1x,3A,T52,A,I4,T62,A,F9.5)')
     &                               'r2r_init args: child grid = ''',
     &                                         chldgrd(1:lcgrd), '''',
     &                                 'N =',N,  'theta_s =', theta_s
                        write(*,'(15x,3A,T62,A,F9.5)')
     &                           'parent grid = ''', prntgrd(1:lpgrd),
     &                                    '''',  'theta_b =', theta_b
                        write(*,'(15x,3A,T52,A,I4,T67,A,F9.3/)')
     &                        'parent data = ''',  prnt_data(1:lprnt),
     &                                '''', 'rec =', rec,  'hc =', hc
                        goto 1 !--> SUCCESS
                      else
                        write(*,'(/1x,3A/)')     '### ERROR: Wrong ',
     &                     'argument #8, rec = ', str(1:lenstr(str))
                      endif
                    else
                      write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #7 ',
     &                        ':: Cannot open ''', prnt_data(1:lprnt),
     &                                       '''.', nf_strerror(ierr)
                    endif
                  else
                    write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #6 ::',
     &                     ' Cannot open ''', prntgrd(1:lpgrd), '''.',
     &                                              nf_strerror(ierr)
                  endif
                else
                  write(*,'(/1x,3A/)') '### ERROR: Wrong argument #5',
     &                                   ', N = ', str(1:lenstr(str))
                endif
              else
                write(*,'(/1x,3A/)') '### ERROR: Wrong argument #4, ',
     &                                   'hc = ',  str(1:lenstr(str))
              endif
            else
              write(*,'(/1x,3A/)')   '### ERROR: Wrong argument #3, ',
     &                               'theta_b = ', str(1:lenstr(str))
            endif
          else
            write(*,'(/1x,3A/)')     '### ERROR: Wrong argument #2, ',
     &                               'theta_s = ', str(1:lenstr(str))
          endif
        else
          write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #1 :: Cannot ',
     &        'open ''', chldgrd(1:lcgrd), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(2(/1x,2A)/1x,A/)')    'This command interpolates ',
     &           'user-specified record of the parent-grid solution', 
     &  'into the child grid and creates an initial condition file ',
     &  'to start the child', 'grid model.  It should be invoked as'
        write(*,'(/2x,2A/)')       'r2r_init child_grid.nc theta_s ',
     &              'theta_b hc N parent_grid.nc parent_data.nc rec'
        write(*,'(3(/1x,2A)/)')    'where theta_s theta_b hc N are ',
     &              'for vertical coordinate of the intended child-',
     &     'grid configuration (parent-grid vertical coordinate is ',
     &     'determined automatically',  'from the parent-grid data ',
     &         'file);  rec is record number within parent_data.nc.'
      endif
      stop

  1   call roms_find_dims(nccgrd, chldgrd, i,j,k)
      ncx=i+2 ;  ncy=j+2
      allocate(Cs_r(N), Cs_w(0:N))
      call set_scoord(theta_s,theta_b, N, Cs_r,Cs_w)
      call roms_find_dims(ncsrc,  prnt_data, i,j,k)
      call roms_check_dims(ncpgrd, prntgrd, i,j,0)
      nx=i+2 ; ny=j+2 ; Np=k
      allocate(Csp_w(0:Np),Csp_r(Np))
      call read_scoord(ncsrc, Np, Csp_r,Csp_w,hcp, VertCoordType)

! Allocate arrays and read coordinates first for the child grid, then
! for the parent.  This leads to a more optimal memory use because
! coordinate arrays xc,yc, xp,yp, xpu,ypu, xpv,ypv  for both child and
! parent grid can be deallocated after computing ip,jp-indices and
! fractional offsets xi,eta.  Keep track of the total amount of
! allocated memory "net_alloc_mem" expressed in 4-byte numbers.

      allocate( csA(ncx,ncy),snA(ncx,ncy), h(ncx,ncy), hprnt(ncx,ncy),
     &               mask(ncx,ncy), umask(ncx-1,ncy), vmask(ncx,ncy-1),
     &                                                     stat=ierr )
      if (ierr == 0) then
        net_alloc_mem=2*(2*N+1) +2*(2*Np+1) +4*2*ncx*ncy +3*ncx*ncy/2
        allocate( xi(ncx,ncy),eta(ncx,ncy), xiu(ncx,ncy),etau(ncx,ncy),
     &                                      xiv(ncx,ncy),etav(ncx,ncy),
     &                                                     stat=ierr )
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem +6*2*ncx*ncy
          allocate( ip(ncx,ncy),jp(ncx,ncy), ipu(ncx,ncy),jpu(ncx,ncy),
     &                                       ipv(ncx,ncy),jpv(ncx,ncy),
     &                                                     stat=ierr )
          if (ierr == 0) then
            net_alloc_mem = net_alloc_mem + 6*ncx*ncy
            allocate( xc(ncx,ncy),yc(ncx,ncy), stat=ierr)
            if (ierr == 0) then
              net_alloc_mem = net_alloc_mem + 2*2*ncx*ncy
            endif
          endif
        endif
      endif
      if (ierr == 0) then
        write(*,'(10x,A,F9.2,1x,A)') 'allocated', dble(net_alloc_mem)
     &                                          /262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)')  '### ERROR: Memory allocation failure 1.'
        stop
      endif

! Read child grid land mask and topography. Because land mask may or
! may not present check for it first using netCDF native functions; if
! there is no variable "mask_rho" in netCDF file, set it to water
! everywhere. Temporarily use array "xc" to read mask from the file,
! thereafter this program uses only integer(kind=2) version of mask.

      ierr=nf_inq_varid(nccgrd, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        ierr=nf_get_var_double(nccgrd, varid, xc) !<-- termorarily
        if (ierr == nf_noerr) then
#ifdef CONTOUR_CHILD_MASK
          call set_contour_mask(ncx,ncy, xc, mask,umask,vmask)
#else
          call         set_mask(ncx,ncy, xc, mask,umask,vmask)
#endif
        else
          write(*,'(/1x,5A/)') '### ERROR: Cannot read ''mask_rho'' ',
     &         'from ''', chldgrd(1:lcgrd), ''', ', nf_strerror(ierr)
          stop
        endif
      else
        mask=1 ; umask=1 ; vmask=1
        write(*,'(9x,4A)')   'No land mask ''mask_rho'' is present ',
     &  'in ''', chldgrd(1:lcgrd), ''', assuming mask=1 everywhere.'
      endif
      call get_var_by_name_double(nccgrd, 'h',  h)

! Read horizontal coordinates for the child grid...

      call get_var_by_name_double(nccgrd, 'lon_rho',  xc)
      call get_var_by_name_double(nccgrd, 'lat_rho',  yc)
      write(*,'(1x,A)',advance='no') 'child grid longitude '
      call compute_min_max(ncx,ncy, xc, xcmin,xcmax)

#ifdef PARENT_GRID_SUBREGION

! Preliminary step: read horizontal coordinates for the parent
!------------ ----- grid, then pretend initializing parent --> child
! interpolation, but actually all what we need is (ip,jp)-indices:
! use them to find the smallest logically rectangular subdomain within
! the parent grid which encloses the child (the unmasked portion of it,
! to be more specific).   This region is characterized by bounds
! imin,imax,jmin,jmax defined within the parent grid index space.
! Expand the region in all four directions by making "margn"-wide halo
! to allow more points to reduce the influence of artificial boundary
! conditions associated with spline interpolations, then save
! coordinates of south-west corner of the expanded region and redefine
! nx,ny consistently with its size.  Thereafter deallocate coordinate
! arrays for both parent and child -- basically all what this step is
! needed for is just 4 integer numbers: iwest,jsouth,nx,ny.

      allocate(xp(nx,ny), yp(nx,ny), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 2*2*nx*ny
      endif

      call get_var_by_name_double(ncpgrd, 'lon_rho',  xp)
      call get_var_by_name_double(ncpgrd, 'lat_rho',  yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)

C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp)
      call r2r_interp_search_thread(  nx,ny, xp,yp, ncx,ncy, xc,yc,
     &                                                        ip,jp)
C$OMP END PARALLEL
      call check_search_indices(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp)

      call compute_index_bounds( ncx,ncy, ip,jp, mask, imin,imax,
     &                                                 jmin,jmax)

      write(*,'(1x,2A/2(4x,A,2I5,1x,A,I5))') 'minimal parent-grid ',
     &                    'index bounds to accommodate child grid:',
     &                     'imin,imax =', imin,imax, 'of nx =', nx,
     &                     'jmin,jmax =', jmin,jmax, 'of ny =', ny
      imin=max(imin-margn,1)  ; jmin=max(jmin-margn,1)
      imax=min(imax+margn,nx) ; jmax=min(jmax+margn,ny)
      write(*,'(1x,A/4x,A,2I5,17x,A,2I5)')  'adjusted to',
     &   'imin,imax =',imin,imax, 'jmin,jmax =',jmin,jmax

      iwest=imin ; nx=imax-imin+1 ; jsouth=jmin ; ny=jmax-jmin+1
      write(*,'(/2(2x,A,I5)/)') 'setting subdomain sizes  nx =', nx,
     &                                                   'ny =', ny
      deallocate(xp,yp)
      net_alloc_mem=net_alloc_mem -2*2*nx*ny
      write(*,'(8x,A,F9.2,1x,A)') 'deallocated',
     &           dble(2*2*nx*ny)/262144.D0,'MB'
#endif

! Note that the only outcome of the code segment above is 4 integer
! numbers,  iwest,jsouth, nx,ny,  which are the indices of south-west
! corner of the subdomain within the parent grid, and the sizes of the
! subdomain.   Everything else is discarded.  Allocate subdomain-sized
! (within the parent grid) 2D arrays: note that xp(nx,ny),yp(nx,ny)
! were deallocated above and now are allocated again with a different
! (expected to be smaller) size).

      allocate( srX(nx*ny), srY(nx*ny), sXY(nx*ny), sYX(nx*ny),
     &                                              stat=ierr )
      if (ierr == 0) then
        allocate(hp(nx,ny), mskp(nx,ny), umsp(nx-1,ny), vmsp(nx,ny-1),
     &           xp(nx,ny), yp(nx,ny),    xpu(nx-1,ny),  xpv(nx,ny-1),
     &                                    ypu(nx-1,ny),  ypv(nx,ny-1),
     &                                                    stat=ierr )
      endif
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem +(4+7)*2*nx*ny +3*(nx*ny)/2
        write(*,'(10x,A,F9.2,1x,A)') 'allocated',dble(net_alloc_mem)
     &                                         /262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 2.'
        stop
      endif

! Re-read (read first time if  CPP PARENT_GRID_SUBREGION is undefined)
! horizontal coordinates for parent grid. however this time only within
! the subdomain of parent grid. Add/subtract 360 degrees to/from
! longitude if necessary to be consistent with child-grid xcmin,xcmax
! range determined above.

      call get_patch_by_name_double(ncpgrd, prntgrd, 'lon_rho',
     &                            iwest,jsouth, nx,ny,0,0, xp)
      call get_patch_by_name_double(ncpgrd, prntgrd, 'lat_rho',
     &                            iwest,jsouth, nx,ny,0,0, yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)

! Initialize horizontal interpolation:

#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(/F8.2,1x,A)') dble(net_clk)/dble(clk_rate),
#else
      write(*,'(/1x,A)')
#endif
     &           'initializing horizontal interpolation...'

C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp, xi,eta)
      call r2r_interp_init_thread(nx,ny, xp,yp, ncx,ncy, xc,yc,
     &                                               ip,jp, xi,eta)
C$OMP END PARALLEL
      call check_search_indices(nx,ny,xp,yp, ncx,ncy, xc,yc, ip,jp)
      call check_offsets(nx,ny,xp,yp, ncx,ncy,xc,yc, ip,jp, xi,eta)

      do j=1,ny
        do i=1,nx-1
          xpu(i,j)=0.5D0*(xp(i,j)+xp(i+1,j))
          ypu(i,j)=0.5D0*(yp(i,j)+yp(i+1,j))
        enddo
      enddo
C$OMP PARALLEL SHARED(nx,ny, xpu,ypu, ncx,ncy,xc,yc, ipu,jpu,xiu,etau)
      call r2r_interp_init_thread(nx-1,ny, xpu,ypu, ncx,ncy, xc,yc,
     &                                               ipu,jpu, xiu,etau)
C$OMP END PARALLEL
      call check_search_indices(nx-1,ny,xpu,ypu, ncx,ncy,xc,yc,ipu,jpu)
      call check_offsets(nx-1,ny,xpu,ypu, ncx,ncy,xc,yc,
     &                                               ipu,jpu,xiu,etau)

      do j=1,ny-1
        do i=1,nx
          xpv(i,j)=0.5D0*(xp(i,j)+xp(i,j+1))
          ypv(i,j)=0.5D0*(yp(i,j)+yp(i,j+1))
        enddo
      enddo
C$OMP PARALLEL SHARED(nx,ny, xpv,ypv, ncx,ncy,xc,yc, ipv,jpv,xiv,etav)
      call r2r_interp_init_thread(nx,ny-1, xpv,ypv, ncx,ncy, xc,yc,
     &                                               ipv,jpv, xiv,etav)
C$OMP END PARALLEL
      call check_search_indices(nx,ny-1,xpv,ypv, ncx,ncy,xc,yc,ipv,jpv)
      call check_offsets(nx,ny-1,xpv,ypv, ncx,ncy,xc,yc,
     &                                               ipv,jpv,xiv,etav)

#ifdef PARENT_GRID_SUBREGION
# ifdef VERBOSE
      call compute_index_bounds( ncx,ncy, ip,jp, mask, imin,imax,
     &                                                 jmin,jmax)

      write(*,'(1x,2A/2(4x,A,2I5,1x,A,I5))') 're-checking parent-',
     &              'grid index bounds to accommodate child grid:',
     &                     'imin,imax =', imin,imax, 'of nx =', nx,
     &                     'jmin,jmax =', jmin,jmax, 'of ny =', ny
      imin=max(imin-margn,1)  ; jmin=max(jmin-margn,1)
      imax=min(imax+margn,nx) ; jmax=min(jmax+margn,ny)
      write(*,'(1x,A/4x,A,2I5,17x,A,2I5)')  'adjusted to',
     &   'imin,imax =',imin,imax, 'jmin,jmax =',jmin,jmax
      if ( imin == 1 .and. imax == nx .and.
     &     jmin == 1 .and. jmax == ny ) then
        write(*,*) 'parent-child bounding check passed'
      else
        write(*,*) '### ERROR: Algorithm failure.' ; stop
      endif
# endif
#endif
#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &    'initialization complete in', dble(inc_clk)/dble(clk_rate),
     &    'sec'
#endif

! Read parent-grid land mask and topography. Similarly to above,
! temporarily use use array "hp" to read "mask_rho" from the file,
! thereafter this program uses only integer(kind=2) version of mask.

      ierr=nf_inq_varid(ncpgrd, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        call get_patch_by_name_double(ncpgrd, prntgrd, 'mask_rho',
     &                               iwest,jsouth, nx,ny,0,0, hp)
        call set_mask(nx,ny, hp, mskp,umsp,vmsp)
      else
        mskp=1 ; umsp=1 ; vmsp=1
        write(*,'(9x,4A)')   'No land mask ''mask_rho'' is present ',
     &  'in ''', prntgrd(1:lpgrd), ''', assuming mask=1 everywhere.'
      endif
      call get_patch_by_name_double(ncpgrd, prntgrd, 'h',
     &                      iwest,jsouth, nx,ny,0,0, hp)

! Interpolate parent grid topography onto child grid...

      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, hp,  srX,srY,sXY, ncx,ncy,
C$OMP&                            ip,jp, xi,eta, hprnt)
      call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY,
     &                   ncx,ncy, ip,jp, xi,eta, hprnt)
C$OMP END PARALLEL
                                                       ! to prevent
      do j=1,ncy                                       ! zeros left by
        do i=1,ncx                                     ! horizontal
          if (hprnt(i,j) < 0.0001D0) hprnt(i,j)=h(i,j) ! interpolation
        enddo                                          ! where parent
      enddo                                            ! and child do
                                                       ! not overlap

! Read angle between true East and local direction of XI-coordinate
! of ROMS grid, then compute csA=cos(alpha) and snA=sin(alpha), first
! for child grid then for the parent. For the latter temporarily place
! the outcome into arrays xp,yp -- after all horizontal interpolations
! have been initialized the content arrays is no longer needed. Then
! interpolate cos(A) and sin(A) of the parent into child grid, then
! compute cos and sin of the child-parent difference of angles; these
! will be used to rotate velocity vector components.

      call read_angle(nccgrd, chldgrd, 1,1,    ncx,ncy,   csA,snA)
      call read_angle(ncpgrd, prntgrd, iwest,jsouth, nx,ny, xp,yp)

      call spln2d_double(nx,ny, xp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, xp, srX,srY,sXY, ncx,ncy,
C$OMP&                              ip,jp, xi,eta, xc)
      call spln2d_interp_double(nx,ny, xp, srX,srY,sXY,
     &                     ncx,ncy, ip,jp, xi,eta, xc)
C$OMP END PARALLEL

      call spln2d_double(nx,ny, yp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, yp, srX,srY,sXY, ncx,ncy,
C$OMP&                              ip,jp, xi,eta, yc)
      call spln2d_interp_double(nx,ny, yp, srX,srY,sXY,
     &                     ncx,ncy, ip,jp, xi,eta, yc)
C$OMP END PARALLEL

! Convert child-grid angles cosA,sinA into cosA=cos(child-parent) and
! snA=sin(child-parent) where the parent-grid angle data is coming from
! above as arrays xc=cos(parent) and yc=sin(parent) interpolated into
! child grid.  Because the intepolation does not preserve the property
! of having cos^2+sin^2=1 exactly, re-normalize them in the process.

      do j=1,ncy
        do i=1,ncx
          cff=1.D0/sqrt(xc(i,j)*xc(i,j)+yc(i,j)*yc(i,j))
          cff1=csA(i,j)*xc(i,j) +snA(i,j)*yc(i,j)
          cff2=snA(i,j)*xc(i,j) -csA(i,j)*yc(i,j)
          csA(i,j)=cff*cff1  ;  snA(i,j)=cff*cff2
        enddo
      enddo

      deallocate(hp, xp,yp, xpu,ypu, xpv,ypv, xc,yc)
      net_alloc_mem=net_alloc_mem-(7*2*nx*ny+2*2*ncy*ncy)
      write(*,'(8x,A,F10.2,1x,A)') 'deallocated',
     &      dble(7*2*nx*ny+2*2*ncy*ncy)/262144.D0,'MB'
      ierr=nf_close(nccgrd) ; ierr=nf_close(ncpgrd)


! Initialize vertical interpolation

      size=ncx*ncy*N ; allocate(kprnt(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 2*size
        write(*,'(/8x,A,F9.2,1x,A/)')    'allocated kprnt, reaching',
     &                      dble(net_alloc_mem)/262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 3.'
        stop
      endif

#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(/F8.2,1x,A)') dble(net_clk)/dble(clk_rate),
#else
      write(*,'(/1x,A)')
#endif
     &          'initializing vertical interpolation...'


C$OMP PARALLEL SHARED(ncx,ncy, hprnt,Np,hcp,Csp_r, h,N,hc,Cs_r, kprnt)
      call r2r_init_vertint_thread(ncx,ncy,  hprnt, Np,hcp,Csp_r,
     &                                          h, N, hc, Cs_r, kprnt)

      call r2r_check_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     &                                          h, N, hc, Cs_r, kprnt)
C$OMP END PARALLEL

#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &    'initialization complete in', dble(inc_clk)/dble(clk_rate),
     &    'sec'
#endif

      call r2r_init_diag_file(ncx,ncy, N, ip,jp,xi,eta,
     &               ipu,jpu,xiu,etau, ipv,jpv,xiv,etav,
     &                        csA,snA, h, hprnt, kprnt)

      write(*,'(1x,A/)') 'initialization complete'

c**      stop !<-- to test initialization


! Allocate large 3D arrays: wrk1, wrk2, wrk3 must be of sufficient
! size to hold one time record of the largest-possible 3D field which
! may be either parent- or child-grid-size variable.  Furthermore,
! horizontal interpolation creates an intermediate field with mixed
! dimensions: horizontally on child grid (ncx,ncy), but still of parent
! size vertically and therefore its total size may be bigger than
! either parent or child [say, if parent has fewer horizontal points
! but more vertical levels than the child].  So the forth array, wrmx,
! is allocated to a size sufficient for all three cases, parent, child,
! and intermediate.

      size=max(nx*ny*(Np+1), ncx*ncy*(N+1))
      allocate(wrk1(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem+size
        write(*,'(8x,A,F10.2,1x,A)')    'allocated wrk1, reaching',
     &                    dble(net_alloc_mem)/262144.D0, 'MB total'
        allocate(wrk2(size), stat=ierr)
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem+size
          write(*,'(8x,A,F10.2,1x,A)')   'allocated wrk2, reaching',
     &                     dble(net_alloc_mem)/262144.D0, 'MB total'
          allocate(wrk3(size), stat=ierr)
          if (ierr == 0) then
            net_alloc_mem=net_alloc_mem+size
            write(*,'(8x,A,F10.2,1x,A)')  'allocated wrk3, reaching',
     &                      dble(net_alloc_mem)/262144.D0, 'MB total'

            size=max(size, ncx*ncy*(Np+1)) !<-- mixed size array
            allocate(wrmx(size), stat=ierr)
            if (ierr == 0) then
              net_alloc_mem=net_alloc_mem+size
              write(*,'(8x,A,F10.2,1x,A/)')'allocated wrmx, reaching',
     &                       dble(net_alloc_mem)/262144.D0, 'MB total'
            endif
          endif
        endif
      endif
      if (ierr /= 0) then
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 4.'
        stop
      endif


!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     ***
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***

! Create target file, save names of the parent grid and data source
! files, record number as global attributes so this target file can
! be reproduced later using exactly the same conditions...

      ierr=nf_create(roms_init, nf_64bit_offset, nctarg)
      if (ierr == nf_noerr) then
        call def_roms_file(nctarg, roms_init, ncx,ncy,N,
     &                     theta_s,theta_b, hc, Cs_w,Cs_r)
        ierr=nf_put_att_text(nctarg, nf_global, 'memo', 19,
     &                              'created by r2r_init')
        if (ierr == nf_noerr) then
          ierr=nf_put_att_text(nctarg, nf_global, 'parent_grid',
     &                                 lpgrd, prntgrd(1:lpgrd))
          if (ierr == nf_noerr) then
            ierr=nf_put_att_text(nctarg, nf_global,
     &            'parent_data_file', lprnt, prnt_data(1:lprnt))
            if (ierr == nf_noerr) then
              ierr=nf_put_att_int(nctarg, nf_global,
     &             'parent_data_file_record', nf_int, 1, rec)
              if (ierr == nf_noerr) then
                write(*,*) 'added r2r_init attributes'
              else
                write(*,*) '### ERROR 5: put_att_int'
              endif
            else
              write(*,*) '### ERROR 4: put_att_text'
            endif
          else
            write(*,*) '### ERROR 3: put_att_text'
          endif
        else
          write(*,*) '### ERROR 2: put_att_text'
        endif
      else
        write(*,*) '### ERROR 1: nf_create'
      endif
      if (ierr /= nf_noerr) stop

! Copy all attributes for time variable while the target file is still
! in redefinition mode, then switch it into input and copy time itself.

      ierr=nf_inq_varid(nctarg, time_var_name, tvar_out)
      if (ierr == nf_noerr) then
        ierr=nf_inq_varid(ncsrc, time_var_name, tvar_in)
        if (ierr == nf_noerr) then
          ierr=nf_inq_varnatts(ncsrc, tvar_in, natts)
          if (ierr == nf_noerr) then
            do i=1,natts
              ierr=nf_inq_attname(ncsrc, tvar_in, i, str)
              if (ierr == nf_noerr) then
                lstr=lenstr(str)
                ierr=nf_copy_att(ncsrc, tvar_in, str(1:lstr),
     &                                     nctarg, tvar_out)
                if (ierr == nf_noerr) then
                  write(*,*) 'copied attribute ''',str(1:lstr),''''
                else
                  write(*,*) '### ERROR 10: copy_att' ; stop
                endif
              else
                write(*,*) '### ERROR 9: inq_attname' ; stop
              endif
            enddo
          else
            write(*,*) '### ERROR 8: inq_varnatts'
          endif
        else
          write(*,*) '### ERROR 7: inq_varid'
        endif
      else
        write(*,*) '### ERROR 6: inq_varid'
      endif
      if (ierr /= nf_noerr) stop

      ierr=nf_enddef(nctarg)  !<-- set to input mode
      if (ierr == nf_noerr) then
        ierr=nf_get_vara_double(ncsrc, tvar_in, rec,1, cff1)
        if (ierr == nf_noerr) then
          ierr=nf_put_var_double(nctarg, tvar_out, cff1)
          if (ierr == nf_noerr) then
            write(*,*) 'wrote ''', time_var_name, ''''
          else
            write(*,*) '### ERROR 13: put_var'
          endif
        else
          write(*,*) '### ERROR 12: get_vara'
        endif
      else
        write(*,*) '### ERROR 11: nf_enddef'
      endif
      if (ierr /= nf_noerr) stop


! Process variables

      call get_patch_by_name_real(ncsrc, prnt_data, 'zeta',
     &                  iwest,jsouth, nx,ny,0,rec,  wrk1)
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1)
      call etch_into_land_thread(nx,ny, mskp, wrk1)
C$OMP END PARALLEL
      call spln2d_real(nx,ny, wrk1, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, wrk1, srX,srY,sXY, ncx,ncy,
C$OMP&                              ip,jp, xi,eta, wrk2)
      call spln2d_interp_real(nx,ny, wrk1,  srX,srY,sXY,
     &                     ncx,ncy, ip,jp, xi,eta, wrk2)
C$OMP END PARALLEL
!     call apply_mask(ncx,ncy,1, mask,wrk2)
      call put_rec_by_name_real(nctarg,roms_init, 'zeta',
     &                                ncx,ncy,0,1, wrk2)
#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
#endif

      call get_patch_by_name_real(ncsrc, prnt_data, 'u', iwest,jsouth,
     &                                         nx-1,ny,Np, rec, wrk1)
      do k=1,Np
        isrc=1+(k-1)*(nx-1)*ny ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, umsp, wrk1, isrc)
        call etch_into_land_thread(nx-1,ny, umsp, wrk1(isrc))
C$OMP END PARALLEL
        call spln2d_real(nx-1,ny, wrk1(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc, wrk1, srX,srY,sXY,
C$OMP&                   ncx,ncy, itrg, ipu,jpu, xiu,etau, wrmx)
        call spln2d_interp_real(nx-1,ny, wrk1(isrc), srX,srY,sXY,
     &                   ncx,ncy, ipu,jpu, xiu,etau, wrmx(itrg))
C$OMP END PARALLEL
        write(6,'(A)',advance='no') '.' ; flush(6)
      enddo
      write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask,  Np,wrmx,  N,kprnt, wrk1)
      call r2r_vrtint_thread(ncx,ncy, 0,mask, btm_slp, Np,wrmx,
     &                                          N,kprnt, wrk1)
C$OMP END PARALLEL

#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
#endif


      call get_patch_by_name_real(ncsrc, prnt_data, 'v', iwest,jsouth,
     &                                         nx,ny-1,Np, rec, wrk2)
      do k=1,Np
        isrc=1+(k-1)*nx*(ny-1) ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, vmsp, wrk2, isrc)
        call etch_into_land_thread(nx,ny-1, vmsp, wrk2(isrc))
C$OMP END PARALLEL
        call spln2d_real(nx,ny-1, wrk2(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc,  wrk2, srX,srY,sXY,
C$OMP&                   ncx,ncy, itrg, ipv,jpv, xiv,etav, wrmx)
        call spln2d_interp_real(nx,ny-1, wrk2(isrc), srX,srY,sXY,
     &                   ncx,ncy, ipv,jpv, xiv,etav, wrmx(itrg))
C$OMP END PARALLEL
        write(6,'(A)',advance='no') '.' ; flush(6)
      enddo
      write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask,  Np,wrmx,  N,kprnt, wrk2)
      call r2r_vrtint_thread(ncx,ncy, 0,mask, btm_slp, Np,wrmx,
     &                                          N,kprnt, wrk2)
C$OMP END PARALLEL

#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
#endif


C$OMP PARALLEL SHARED(ncx,ncy, N, csA,snA, wrk1,wrk2, wrk3,wrmx)
      call r2r_rotate_shift_thread(ncx,ncy,N, csA,snA, wrk1,wrk2,
     &                                                 wrk3,wrmx)
!$OMP END PARALLEL
!     call apply_mask(ncx-1,ncy,N, umask, wrk3)
!     call apply_mask(ncx,ncy-1,N, vmask, wrmx)

      call put_rec_by_name_real(nctarg, roms_init, 'u',
     &                              ncx-1,ncy,N,1,  wrk3)
      call put_rec_by_name_real(nctarg, roms_init, 'v',
     &                              ncx,ncy-1,N,1,  wrmx)


C$OMP PARALLEL SHARED(ncx,ncy, N, hc,Cs_w, h, wrk3,wrmx, wrk1,wrk2)
      call compute_uvbar_thread(ncx,ncy, N, hc,Cs_w, h,  wrk3,wrmx,
     &                                                   wrk1,wrk2)
C$OMP END PARALLEL
!     call apply_mask(ncx-1,ncy,1, umask, wrk1)
!     call apply_mask(ncx,ncy-1,1, vmask, wrk2)

      call put_rec_by_name_real(nctarg, roms_init, 'ubar',
     &                              ncx-1,ncy,0,1,  wrk1)
      call put_rec_by_name_real(nctarg, roms_init, 'vbar',
     &                              ncx,ncy-1,0,1,  wrk2)
#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
#endif

      do itrc=1,ntrc
        call get_patch_by_name_real(ncsrc, prnt_data, tracer(itrc),
     &                          iwest,jsouth, nx,ny,Np,rec, wrk1)
        do k=1,Np
          isrc=1+(k-1)*nx*ny ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, wrk1, isrc)
          call etch_into_land_thread(nx,ny, mskp, wrk1(isrc))
C$OMP END PARALLEL
          call spln2d_real(nx,ny, wrk1(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc,  wrk1, srX,srY,sXY,
C$OMP&                    ncx,ncy, itrg, ip,jp, xi,eta, wrmx)
          call spln2d_interp_real(nx,ny, wrk1(isrc), srX,srY,sXY,
     &                      ncx,ncy, ip,jp, xi,eta, wrmx(itrg))
C$OMP END PARALLEL
          write(6,'(A)',advance='no') '.' ; flush(6)
        enddo
        write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask, Np,wrmx, N,kprnt,wrk1)
      call r2r_vrtint_thread(ncx,ncy, 1,mask, btm_trc, Np,wrmx,
     &                                        N,kprnt,wrk1)
C$OMP END PARALLEL
        call put_rec_by_name_real(nctarg,roms_init, tracer(itrc),
     &                                        ncx,ncy,N,1, wrk1)
#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
#endif
      enddo !<--itrc

      ierr=nf_close(nctarg)
#ifdef TIMING
      write(*,'(/4x,A,F12.3,1x,A,F10.1,1x,A)')    'total data read ',
     &                  dble(sz_read_acc)/dble(1024**2), 'MBytes in',
     &                          dble(read_clk)/dble(clk_rate), 'sec'
      write(*,'(1x,A,F12.3,1x,A,F10.1,1x,A/)') 'total data written ',
     &                 dble(sz_write_acc)/dble(1024**2), 'MBytes in',
     &                         dble(write_clk)/dble(clk_rate), 'sec'
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate, clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      call cpu_time(tend)
      write(*,'(1x,A,F8.2,1x,A,4x,A,F8.2,1x,A,F8.1,1x,A)')
     &   'Wall Clock time:', dble(net_clk)/dble(clk_rate), 'sec',
     &                           'CPU time:', tend-tstart, 'sec',
     & (tend-tstart)*dble(clk_rate)/dble(net_clk)*100.D0,'% CPU'
#endif
      lstr=lenstr(roms_init)
      write(*,'(/1x,5A/)')  'Files ''', roms_init(1:lstr),
     &     ''' and ''', 'r2r_init_diag.nc', ''' are ready.'
      end

#include "def_roms_file.F"
