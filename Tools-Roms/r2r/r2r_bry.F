      program r2r_bry

! Creates open boundary forcing file by extracting and interpolating
! relevant data from parent-grid solution. Should be used as
!
!    r2r_bry child_grid.nc theta_s theta_b hc N WESN parent_grid.nc \
!                                 parent_his1.nc parent_his2.nc ...
!
! where theta_s theta_b hc N specify vertical coordinate of child-grid
! configuration; argument WESN (stands for West, East, North, South)
! specifies which sides are open; arguments beyond the parent grid file
! name are parent-grid history files, may be one or more, may be named
! individually or via wildcards.

! CPP-switch PARENT_GRID_SUBREGION make this program read only relevant
! subdomain within the parent grid as opposite to reading the whole
! data.  Refer to "r2r_init.F" for pros and cons of selecting it.

! Method: 2D horizontal interpolation is by bi-cubic pseudo-splines
! via Hermite basis functions followed by vertical splines by spline
! inversion (z-levels of child grid are translated into continuous
! k-index space of parent grid in such a way spline interpolation of
! the resultant kp=kp(z_r) values back into z-space of child grid
! yields the initial z_r=z_r(k) exactly. In doing so verticall spline
! derivatives for both z_r of parent grid and for the 3D fields to be
! interpolated are constructed in k-index space of horizontally-child
! but still vertically-parent grid). All horizontal interpolations are
! direct from source to destination locations on C-grids:
!              RHO_parent --> RHO_child, U-->U, V-->V,
!                V-->U, V-->U,  RHO-->U, and RHO-->V
! overall a total of 7 permutations. Cross-interpolations, V-->U and
! V-->U, are needed for the rotation of vector components; RHO-->U and
! RHO-->V are is for computing barotropic mode (vertical averaging).
! Before horizontal interpolation data is extended to land-masked areas
! by etching.

! Parameters setting "btm_slp" and "btm_trc" control type of boundary
! condition at bottom for the purpose of computing spline derivatives:
! btm_slp=-1 is no-slip, meaning that u,v=0 exactly at the bottom, i.e.
! half-grid-interval below the k=1 value; This applies to velocities
! only; for tracers the meaningful b.c. are natural on Neumann (0 and
! +1 respectively). Top b.c. is always natural.

! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com

#define PARENT_GRID_SUBREGION
c--#define VERBOSE
#define TIMING
      use mod_io_size_acct
      implicit none
      character(len=80) chldgrd, prntgrd
      character(len=80), allocatable, dimension(:) :: prnt_data
      character(len=11), parameter :: roms_bry='roms_bry.nc'
      character(len=10), parameter :: time_var_name='ocean_time'
      character(len=16) VertCoordType, str, vname, trgname, time_units
      integer, parameter :: ntrc=2
      character(len=4) tracer(ntrc)
      integer, parameter :: btm_slp=-1 ! no-slip=-1; +1 free; 0 natural
     &                    , btm_trc=+1 ! +1 Neumann; 0 natural;
      logical OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH
      real(kind=8) theta_s,theta_b, hc,hcp, xcmin,xcmax, time,
     &                                         cff, cs,csP,sn,snP
      real(kind=8), allocatable, dimension(:,:) :: csA,snA, h,hprnt,
     &                             hp, xc,yc, xp,yp, xpu,ypu, xpv,ypv
      integer(kind=2), allocatable, dimension(:,:) :: mskp, umsp,vmsp
      real(kind=8), allocatable, dimension(:) :: Cs_w,Cs_r,Csp_w,Csp_r,
     &  xr_west,yr_west,    xu_west,yu_west,     xiur_west,etaur_west,
     &                     xiu_west,etau_west,   xiuv_west,etauv_west,
     &  xi_west,eta_west,   xv_west,yv_west,     xivr_west,etavr_west,
     &                     xiv_west,etav_west,   xivu_west,etavu_west,
     &                    csAu_west,snAu_west,   csAv_west,snAv_west,
     &   h_west,hp_west,    hu_west,hpu_west,      hv_west,hpv_west,
     &          kpr_west,           kpu_west,              kpv_west,

     &  xr_east,yr_east,    xu_east,yu_east,     xiur_east,etaur_east,
     &                     xiu_east,etau_east,   xiuv_east,etauv_east,
     &  xi_east,eta_east,   xv_east,yv_east,     xivr_east,etavr_east,
     &                     xiv_east,etav_east,   xivu_east,etavu_east,
     &                    csAu_east,snAu_east,   csAv_east,snAv_east,
     &   h_east,hp_east,    hu_east,hpu_east,      hv_east,hpv_east,
     &          kpr_east,           kpu_east,              kpv_east

       real(kind=8), allocatable, dimension(:) ::
     & xr_south,yr_south,  xu_south,yu_south,   xiur_south,etaur_south,
     &                    xiu_south,etau_south, xiuv_south,etauv_south,
     & xi_south,eta_south, xv_south,yv_south,   xivr_south,etavr_south,
     &                    xiv_south,etav_south, xivu_south,etavu_south,
     &                   csAu_south,snAu_south, csAv_south,snAv_south,
     &   h_south,hp_south, hu_south,hpu_south,   hv_south,hpv_south,
     &          kpr_south,          kpu_south,            kpv_south,

     & xr_north,yr_north,  xu_north,yu_north,   xiur_north,etaur_north,
     &                    xiu_north,etau_north, xiuv_north,etauv_north,
     & xi_north,eta_north, xv_north,yv_north,   xivr_north,etavr_north,
     &                    xiv_north,etav_north, xivu_north,etavu_north,
     &                   csAu_north,snAu_north, csAv_north,snAv_north,
     &   h_north,hp_north, hu_north,hpu_north,   hv_north,hpv_north,
     &          kpr_north,          kpu_north,            kpv_north

      integer(kind=4), allocatable, dimension(:) ::
     &     ir_west,jr_west,    iu_west,ju_west,      iv_west,jv_west,
     &                        iur_west,jur_west,    ivr_west,jvr_west,
     &                        iuv_west,juv_west,    ivu_west,jvu_west,

     &     ir_east,jr_east,    iu_east,ju_east,      iv_east,jv_east,
     &                        iur_east,jur_east,    ivr_east,jvr_east,
     &                        iuv_east,juv_east,    ivu_east,jvu_east,

     &    ir_south,jr_south,  iu_south,ju_south,    iv_south,jv_south,
     &                       iur_south,jur_south,  ivr_south,jvr_south,
     &                       iuv_south,juv_south,  ivu_south,jvu_south,

     &    ir_north,jr_north,  iu_north,ju_north,    iv_north,jv_north,
     &                       iur_north,jur_north,  ivr_north,jvr_north,
     &                       iuv_north,juv_north,  ivu_north,jvu_north
      integer(kind=2), allocatable, dimension(:) ::
     &                    mskr_west, mskr_east, mskr_south, mskr_north,
     &                    msku_west, msku_east, msku_south, msku_north,
     &                    mskv_west, mskv_east, mskv_south, mskv_north

      real(kind=8), allocatable, dimension(:) :: srX, srY, sXY, sYX
      real(kind=4), allocatable, dimension(:) :: wrk1,wrk2, wrk3,wrk4,
     &                                           wrk5,wrk6, wrk7,wrk8
      integer net_alloc_mem, nargs, nx,ny,Np, ncx,ncy,N, nhists,ihis,
     &        nccgrd, ncpgrd, ncsrc, nctarg, nrecs,rec, recout, varid,
     &        vtype, ndims, vdimids(4), natts, tvar_in,tvar_out, ierr,
     &        itrc,size,isrc, i,j,k, lstr,lvar,ltgv,lcgrd,lpgrd,lprnt

      integer lenstr, read_str_int, read_str_dble !<-- functions

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
      ihis=1 ; rec=1; recout=0 ; time_units='                  '
      theta_s=0.D0 ; hc=0.D0 ; OBC_WEST=.false. ; OBC_SOUTH=.false.
      theta_b=0.D0 ; N=0     ; OBC_EAST=.false. ; OBC_NORTH=.false.

      nargs=iargc()
      if (nargs > 7) then
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
                  call getarg(6,str) ; lvar=lenstr(str)
                  do i=1,lvar
                    if (str(i:i) == 'W'  .or.  str(i:i) == 'w') then
                      OBC_WEST=.true.
                    elseif (str(i:i) == 'E' .or. str(i:i) == 'e') then
                      OBC_EAST=.true.
                    elseif (str(i:i) == 'S' .or. str(i:i) == 's') then
                      OBC_SOUTH=.true.
                    elseif (str(i:i) == 'N' .or. str(i:i) == 'n') then
                      OBC_NORTH=.true.
                    else
                      write(*,'(/1x,6A/12x,A/)')  '### ERROR: Wrong ',
     &                  'argument # 6 ''',str(1:lvar),''': letter ''',
     &                          str(i:i), ''' should not be present.',
     &                     'Only "W", "E", "S", and "N" are allowed.'
                      stop
                    endif
                  enddo
                  call getarg(7,prntgrd) ; lpgrd=lenstr(prntgrd)
                  ierr=nf_open(prntgrd(1:lpgrd), nf_nowrite, ncpgrd)
                  if (ierr == nf_noerr) then
                    nhists=nargs-7                 ! Open and close
                    allocate(prnt_data(nhists))    ! the files to check
                    do i=1,nhists                  ! that they are OK,
                      call getarg(i+7,prnt_data(i))! however leave the
                      lprnt=lenstr(prnt_data(i))   ! first one open.
                      if (i == 1) then
                        ierr=nf_open(prnt_data(i)(1:lprnt),
     &                                 nf_nowrite, ncsrc)
                      else                                  ! using "j"
                        ierr=nf_open(prnt_data(i)(1:lprnt), ! as netCDF
     &                                      nf_nowrite, j)  ! file ID
                      endif                                 ! here is
                      if (ierr == nf_noerr) then            ! just for
                        if (i > 1) ierr=nf_close(j)         ! testing
                      else
                        write(*,'(/1x,A,I3,1x,3A/24x,A/)')
     &                  '### ERROR: arg #', i+7,   ':: Cannot open ''',
     &                   prnt_data(i)(1:lprnt),'''.', nf_strerror(ierr)
                        stop
                      endif
                    enddo

                    write(*,'(/1x,3A,T52,A,I4,T62,A,F9.5)')
     &               'r2r_bry args: child grid = ''',chldgrd(1:lcgrd),
     &                           '''', 'N =',N,  'theta_s =', theta_s
                    write(*,'(15x,3A,T62,A,F9.5)') 'parent grid = ''',
     &                  prntgrd(1:lpgrd), '''',  'theta_b =', theta_b
                    write(*,'(12x,4(3x,A,1x,L1),T67,A,F9.3)')
     &                'W =',OBC_WEST, 'E =',OBC_EAST, 'S =',OBC_SOUTH,
     &                                   'N =', OBC_NORTH, 'hc =', hc
#ifdef VERBOSE
                    do i=1,nhists
                      write(*,'(12x,I4,1x,3A)')  i, 'src file = ''',
     &                   prnt_data(i)(1:lenstr(prnt_data(i))), ''''
                    enddo
#endif
                    write(*,*)
                    goto 1 !--> SUCCESS
                  else
                    write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #7 ::',
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
        write(*,'(2(/1x,2A)/,/8x,2A/)')   'This command extracts ',
     &    'side-boundary forcing data from the parent-grid solution',
     &    'and interpolates it into perimeter of the child grid. It',
     &                                       ' should be invoked as',
     &            'r2r_bry child_grid.nc theta_s theta_b hc N WESN ',
     &            'parent_grid.nc parent_his1.nc parent_his2.nc ...'

        write(*,'(4(/1x,2A)/)')    'where theta_s theta_b hc N are ',
     &              'for vertical coordinate of the intended child-',
     &     'grid configuration (parent-grid vertical coordinate is ',
     &     'determined automatically',  'from the parent-grid data ',
     &        'file); argument WESN is to indicate which among the',
     &       'west, east, south, north sides of the child domain ',
     &                                     'have open boundaries.'
      endif
      stop

  1   call roms_find_dims(nccgrd, chldgrd, i,j,k)
      ncx=i+2 ;  ncy=j+2
      call roms_find_dims(ncsrc, prnt_data(1), i,j,k)
      call roms_check_dims(ncpgrd, prntgrd, i,j,0)
      nx=i+2 ; ny=j+2 ; Np=k
      allocate(Cs_r(N), Cs_w(0:N), Csp_w(0:Np),Csp_r(Np))
      call set_scoord(theta_s,theta_b, N, Cs_r,Cs_w)
      call read_scoord(ncsrc, Np, Csp_r,Csp_w,hcp, VertCoordType)
      net_alloc_mem=2*(2*N+1)+ 2*(2*Np+1) ; ierr=0

      if (OBC_WEST) then
        allocate( kpr_west(ncy*N), kpu_west(ncy*N), kpv_west(ncy*N-N),
     &       xr_west(ncy),yr_west(ncy), xu_west(ncy),yu_west(ncy),
     &                                csAu_west(ncy),snAu_west(ncy),
     &     h_west(ncy),hp_west(ncy),    hu_west(ncy),hpu_west(ncy),
     &                                xv_west(ncy-1),yv_west(ncy-1),
     &                              csAv_west(ncy-1),snAv_west(ncy-1),
     &                                hv_west(ncy-1),hpv_west(ncy-1),

     &     xi_west(ncy),eta_west(ncy), xiu_west(ncy),etau_west(ncy),
     &                                xiuv_west(ncy),etauv_west(ncy),
     &                                xiur_west(ncy),etaur_west(ncy),
     &                              xivr_west(ncy-1),etavr_west(ncy-1),
     &                              xivu_west(ncy-1),etavu_west(ncy-1),
     &                               xiv_west(ncy-1),etav_west(ncy-1),

     &       ir_west(ncy),jr_west(ncy), iu_west(ncy),ju_west(ncy),
     &                                 iuv_west(ncy),juv_west(ncy),
     &                                 iur_west(ncy),jur_west(ncy),
     &                               ivr_west(ncy-1),jvr_west(ncy-1),
     &                               ivu_west(ncy-1),jvu_west(ncy-1),
     &                                iv_west(ncy-1),jv_west(ncy-1),
     &             mskr_west(ncy), msku_west(ncy), mskv_west(ncy-1),
     &                                                stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncy*N + 2*(N*ncy-1)
     &   +18*2*ncy +12*2*(ncy-1) +8*ncy +6*(ncy-1) +(3*ncy-1)/2
        write(*,*) 'allocated western  side coordinate arrays'
      endif

      if (OBC_EAST  .and. ierr == 0 ) then
        allocate(kpr_east(ncy*N), kpu_east(ncy*N), kpv_east(ncy*N-N),
     &       xr_east(ncy),yr_east(ncy), xu_east(ncy),yu_east(ncy),
     &                                csAu_east(ncy),snAu_east(ncy),
     &     h_east(ncy),hp_east(ncy),    hu_east(ncy),hpu_east(ncy),
     &                                xv_east(ncy-1),yv_east(ncy-1),
     &                              csAv_east(ncy-1),snAv_east(ncy-1),
     &                                hv_east(ncy-1),hpv_east(ncy-1),

     &     xi_east(ncy),eta_east(ncy), xiu_east(ncy),etau_east(ncy),
     &                                xiuv_east(ncy),etauv_east(ncy),
     &                                xiur_east(ncy),etaur_east(ncy),
     &                              xivr_east(ncy-1),etavr_east(ncy-1),
     &                              xivu_east(ncy-1),etavu_east(ncy-1),
     &                               xiv_east(ncy-1),etav_east(ncy-1),

     &       ir_east(ncy),jr_east(ncy), iu_east(ncy),ju_east(ncy),
     &                                 iuv_east(ncy),juv_east(ncy),
     &                                 iur_east(ncy),jur_east(ncy),
     &                               ivr_east(ncy-1),jvr_east(ncy-1),
     &                               ivu_east(ncy-1),jvu_east(ncy-1),
     &                                iv_east(ncy-1),jv_east(ncy-1),
     &             mskr_east(ncy), msku_east(ncy), mskv_east(ncy-1),
     &                                                stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncy*N + 2*(N*ncy-1)
     &    +18*2*ncy +12*2*(ncy-1) +8*ncy +6*(ncy-1) +(3*ncy-1)/2
        write(*,*) 'allocated eastern  side coordinate arrays'
      endif

      if (OBC_SOUTH  .and. ierr == 0 ) then
        allocate(kpr_south(ncx*N),kpu_south(ncx*N-N),kpv_south(ncx*N),
     &   xr_south(ncx),yr_south(ncx), xv_south(ncx),yv_south(ncx),
     &                              csAv_south(ncx),snAv_south(ncx),
     &    h_south(ncx),hp_south(ncx), hv_south(ncx),hpv_south(ncx),
     &                              xu_south(ncx-1),yu_south(ncx-1),
     &                            csAu_south(ncx-1),snAu_south(ncx-1),
     &                              hu_south(ncx-1),hpu_south(ncx-1),
     &       xi_south(ncx),eta_south(ncx),
     &                               xiv_south(ncx),etav_south(ncx),
     &                              xivu_south(ncx),etavu_south(ncx),
     &                              xivr_south(ncx),etavr_south(ncx),
     &                            xiur_south(ncx-1),etaur_south(ncx-1),
     &                            xiuv_south(ncx-1),etauv_south(ncx-1),
     &                             xiu_south(ncx-1),etau_south(ncx-1),
     &      ir_south(ncx),jr_south(ncx),
     &                                iv_south(ncx),jv_south(ncx),
     &                               ivu_south(ncx),jvu_south(ncx),
     &                               ivr_south(ncx),jvr_south(ncx),
     &                             iur_south(ncx-1),jur_south(ncx-1),
     &                             iuv_south(ncx-1),juv_south(ncx-1),
     &                              iu_south(ncx-1),ju_south(ncx-1),
     &        mskr_south(ncx), msku_south(ncx-1), mskv_south(ncx),
     &                                               stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncx*N + 2*(N*ncx-1)
     &   +18*2*ncx +12*2*(ncx-1) +8*ncx +6*(ncx-1)+(3*ncx-1)/2
        write(*,*) 'allocated southern side coordinate arrays'
      endif

      if (OBC_NORTH  .and. ierr == 0 ) then
        allocate(kpr_north(ncx*N),kpu_north(ncx*N-N),kpv_north(ncx*N),
     &   xr_north(ncx),yr_north(ncx), xv_north(ncx),yv_north(ncx),
     &                              csAv_north(ncx),snAv_north(ncx),
     &    h_north(ncx),hp_north(ncx), hv_north(ncx),hpv_north(ncx),
     &                              xu_north(ncx-1),yu_north(ncx-1),
     &                            csAu_north(ncx-1),snAu_north(ncx-1),
     &                              hu_north(ncx-1),hpu_north(ncx-1),
     &       xi_north(ncx),eta_north(ncx),
     &                               xiv_north(ncx),etav_north(ncx),
     &                              xivu_north(ncx),etavu_north(ncx),
     &                              xivr_north(ncx),etavr_north(ncx),
     &                            xiur_north(ncx-1),etaur_north(ncx-1),
     &                            xiuv_north(ncx-1),etauv_north(ncx-1),
     &                             xiu_north(ncx-1),etau_north(ncx-1),
     &      ir_north(ncx),jr_north(ncx),
     &                                iv_north(ncx),jv_north(ncx),
     &                               ivu_north(ncx),jvu_north(ncx),
     &                               ivr_north(ncx),jvr_north(ncx),
     &                             iur_north(ncx-1),jur_north(ncx-1),
     &                             iuv_north(ncx-1),juv_north(ncx-1),
     &                              iu_north(ncx-1),ju_north(ncx-1),
     &        mskr_north(ncx), msku_north(ncx-1), mskv_north(ncx),
     &                                               stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncx*N + 2*(N*ncx-1)
     &   +18*2*ncx +12*2*(ncx-1) +8*ncx +6*(ncx-1) +(3*ncx-1)/2
        write(*,*) 'allocated northern side coordinate arrays'
      endif
      write(*,'(8x,A,F9.2,1x,A)') 'allocated', dble(net_alloc_mem)
     &                                        /262144.D0,'MB total'

      allocate( csA(ncx,ncy),snA(ncx,ncy), h(ncx,ncy), hprnt(ncx,ncy),
     &                                                     stat=ierr)
      if (ierr == 0 ) then
        net_alloc_mem=net_alloc_mem +4*2*ncx*ncy
        allocate(xc(ncx,ncy),yc(ncx,ncy), stat=ierr)
        if (ierr == 0 ) then
          net_alloc_mem=net_alloc_mem +2*2*ncx*ncy
#ifdef PARENT_GRID_SUBREGION
          allocate(xp(nx,ny),yp(nx,ny), stat=ierr)
          if (ierr == 0) then
            net_alloc_mem=net_alloc_mem +2*2*nx*ny
          endif
#endif
        endif
      endif
      if (ierr == 0) then
        write(*,'(8x,A,F9.2,1x,A)') 'allocated', dble(net_alloc_mem)
     &                                          /262144.D0,'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 1.'
        stop
      endif

! Setup child-grid land mask: ! This program does not use 2D arrays
!------ ---------- ---- ----- ! for child-grid mask, but rather sets up
                              ! individual 1D arrays for RHO- U- and V-
      if (OBC_WEST) then      ! points at each open boundary.  Because
        mskr_west(:)=1        ! land mask may or may not be present for
        msku_west(:)=1        ! for the particular configuration, set
        mskv_west(:)=1        ! all masks to "water" everywhere status
      endif                   ! by assigning mask=1 (on the left), then
      if (OBC_EAST) then      ! using native netCDF functions check
        mskr_east(:)=1        ! whether variable "mask_rho" is present
        msku_east(:)=1        ! in netCDF file: if yes, proceed reading
        mskv_east(:)=1        ! and set mask=0 where it is land; if no
      endif                   ! then all-water status on the left will
      if (OBC_SOUTH) then     ! remain unchanged.
        mskr_south(:)=1
        mskv_south(:)=1       ! Temporarily use array "h" to read mask
        msku_south(:)=1       ! from the file as it is expected to be
      endif                   ! of real type there -- thereafter this
      if (OBC_NORTH) then     ! program uses only 2-byte integer
        mskr_north(:)=1       ! version of mask.
        mskv_north(:)=1
        msku_north(:)=1
      endif                                         ! Note that code
                                                    ! below relies on
      ierr=nf_inq_varid(nccgrd, 'mask_rho', varid)  ! having all mask
      if (ierr == nf_noerr) then                    ! arrays mskX_side
        ierr=nf_get_var_double(nccgrd, varid, h)    ! initialized to 1
        if (ierr == nf_noerr) then                  ! everywhere.
#define ORAND .and.
          if (OBC_WEST) then
            do j=1,ncy
              if (h(1,j) < 0.5D0) mskr_west(j)=0
              if (h(1,j) < 0.5D0 ORAND  h(2,j) < 0.5D0) msku_west(j)=0
            enddo
            do j=1,ncy-1
              if (h(1,j) < 0.5D0 ORAND h(1,j+1) < 0.5D0) mskv_west(j)=0
            enddo
          endif
          if (OBC_EAST) then
            do j=1,ncy
              if (h(ncx,j) < 0.5D0) mskr_east(j)=0
              if (h(ncx,j) < 0.5D0 ORAND h(ncx-1,j) < 0.5D0)
     &                                       msku_east(j)=0
            enddo
            do j=1,ncy-1
              if (h(ncx,j) < 0.5D0 ORAND h(ncx,j+1) < 0.5D0)
     &                                       mskv_east(j)=0
            enddo
          endif
          if (OBC_SOUTH) then
            do i=1,ncx
              if (h(i,1) < 0.5D0) mskr_south(i)=0
              if (h(i,1) < 0.5D0 ORAND h(i,2) < 0.5D0) mskv_south(i)=0
            enddo
            do i=1,ncx-1
              if (h(i,1) < 0.5D0 ORAND h(i+1,1)< 0.5D0) msku_south(i)=0
            enddo
          endif
          if (OBC_NORTH) then
            do i=1,ncx
              if (h(i,ncy) < 0.5D0) mskr_north(i)=0
              if (h(i,ncy) < 0.5D0 ORAND h(i,ncy-1) < 0.5D0)
     &                                      mskv_north(i)=0
            enddo
            do i=1,ncx-1
              if (h(i,ncy) < 0.5D0 ORAND h(i+1,ncy) < 0.5D0)
     &                                      msku_north(i)=0
            enddo
          endif
#undef ORAND
        else
          write(*,'(/1x,5A/)') '### ERROR: Cannot read ''mask_rho'' ',
     &        'from ''', chldgrd(1:lcgrd), ''', ', nf_strerror(ierr)
          stop
        endif
      else
        write(*,'(9x,4A)')   'Variable ''mask_rho'' does not exist ',
     &  'in ''', chldgrd(1:lcgrd), ''', assuming mask=1 everywhere.'
      endif

! Read child-grid horizontal coordinates, extract and save them along
!------ ---------- --------- -----------  the open boundary lines. Only
! "r"-versions, xr_side,yr_side, are used initially for the parent grid
! index search in the following segment below.  The other two, "u"- and
! "v"-types will used during the final search and initialization of the
! horizontal interpolation after identifying the relevant subdomain
! within the parent-grid.

      call get_var_by_name_double(nccgrd, 'lon_rho',  xc)
      call get_var_by_name_double(nccgrd, 'lat_rho',  yc)

      if (OBC_WEST) then
        do j=1,ncy
          xr_west(j)=xc(1,j); yr_west(j)=yc(1,j)
          xu_west(j)=0.5D0*(xc(1,j)+xc(2,j))
          yu_west(j)=0.5D0*(yc(1,j)+yc(2,j))
        enddo
        do j=1,ncy-1
          xv_west(j)=0.5D0*(xc(1,j)+xc(1,j+1))
          yv_west(j)=0.5D0*(yc(1,j)+yc(1,j+1))
        enddo
      endif

      if (OBC_EAST) then
        do j=1,ncy
          xr_east(j)=xc(ncx,j); yr_east(j)=yc(ncx,j)
          xu_east(j)=0.5D0*(xc(ncx,j)+xc(ncx-1,j))
          yu_east(j)=0.5D0*(yc(ncx,j)+yc(ncx-1,j))
        enddo
        do j=1,ncy-1
          xv_east(j)=0.5D0*(xc(ncx,j)+xc(ncx,j+1))
          yv_east(j)=0.5D0*(yc(ncx,j)+yc(ncx,j+1))
        enddo
      endif

      if (OBC_SOUTH) then
        do i=1,ncx
          xr_south(i)=xc(i,1); yr_south(i)=yc(i,1)
          xv_south(i)=0.5D0*(xc(i,1)+xc(i,2))
          yv_south(i)=0.5D0*(yc(i,1)+yc(i,2))
        enddo
        do i=1,ncx-1
          xu_south(i)=0.5D0*(xc(i,1)+xc(i+1,1))
          yu_south(i)=0.5D0*(yc(i,1)+yc(i+1,1))
        enddo
      endif

      if (OBC_NORTH) then
        do i=1,ncx
          xr_north(i)=xc(i,ncy) ; yr_north(i)=yc(i,ncy)
          xv_north(i)=0.5D0*(xc(i,ncy)+xc(i,ncy-1))
          yv_north(i)=0.5D0*(yc(i,ncy)+yc(i,ncy-1))
        enddo
        do i=1,ncx-1
          xu_north(i)=0.5D0*(xc(i,ncy)+xc(i+1,ncy))
          yu_north(i)=0.5D0*(yc(i,ncy)+yc(i+1,ncy))
        enddo
      endif

      write(*,'(1x,A)',advance='no') 'child grid longitude '
      call compute_min_max(ncx,ncy, xc, xcmin,xcmax)

#ifdef PARENT_GRID_SUBREGION

! Preliminary step: read horizontal coordinates for the parent
!------------ ----- grid. then pretend initializing parent --> child
! interpolation, but actually all what we need is (ip,jp)-indices: use
! them to find the smallest logically rectangular subdomain within the
! parent grid which encloses the child (to be more specifically the
! unmasked portion of it).  This region is characterized by bounds
! imin,imax,jmin,jmax defined within the parent grid index space.
! Expand the region in all four directions by making "margn"-wide halo
! to allow more points to reduce the influence of artificial boundary
! conditions associated with spline interpolations, then save
! coordinates of south-west corner of the expanded region and redefine
! nx,ny consistently with its size.  Thereafter deallocate coordinate
! arrays for both parent and child -- basically all what this step is
! needed for is just 4 integer numbers: iwest,jsouth,nx,ny.

      call get_var_by_name_double(ncpgrd, 'lon_rho',  xp)
      call get_var_by_name_double(ncpgrd, 'lat_rho',  yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)

      imin=nx+1 ; imax=-1  !<-- initialize to unrealistic
      jmin=ny+1 ; jmax=-1  !<-- values outside the range

      if (OBC_WEST) then
        call r2r_bry_search( nx,ny, xp,yp,  ncy, xr_west,yr_west,
     &                                             ir_west,jr_west)
        call r2r_bry_index_bounds(ncy,  ir_west,jr_west,  mskr_west,
     &                                         imin,imax,jmin,jmax)
      endif

      if (OBC_EAST) then
        call r2r_bry_search( nx,ny, xp,yp,  ncy, xr_east,yr_east,
     &                                             ir_east,jr_east)
        call r2r_bry_index_bounds(ncy,  ir_east,jr_east,  mskr_east,
     &                                         imin,imax,jmin,jmax)
      endif

      if (OBC_SOUTH) then
        call r2r_bry_search( nx,ny, xp,yp,  ncx, xr_south,yr_south,
     &                                            ir_south,jr_south)
        call r2r_bry_index_bounds(ncx, ir_south,jr_south, mskr_south,
     &                                          imin,imax,jmin,jmax)
      endif

      if (OBC_NORTH) then
        call r2r_bry_search( nx,ny, xp,yp,  ncx, xr_north,yr_north,
     &                                            ir_north,jr_north)
        call r2r_bry_index_bounds(ncx, ir_north,jr_north, mskr_north,
     &                                          imin,imax,jmin,jmax)
      endif

      imax=imax+1 ; jmax=jmax+1
      write(*,'(1x,2A/2(4x,A,2I5,1x,A,I5))')  'minimal parent-grid ',
     &                     'index bounds to accommodate child grid:',
     &                      'imin,imax =', imin,imax, 'of nx =', nx,
     &                      'jmin,jmax =', jmin,jmax, 'of ny =', ny
      imin=max(imin-margn,1)  ; jmin=max(jmin-margn,1)
      imax=min(imax+margn,nx) ; jmax=min(jmax+margn,ny)
      write(*,'(1x,A/4x,A,2I5,17x,A,2I5)')  'adjusted to',
     &   'imin,imax =',imin,imax, 'jmin,jmax =',jmin,jmax

      iwest=imin ; nx=imax-imin+1 ; jsouth=jmin ; ny=jmax-jmin+1
      write(*,'(/2(2x,A,I5)/)') 'setting subdomain sizes  nx =', nx,
     &                                                   'ny =', ny
      deallocate(xp,yp)
      net_alloc_mem=net_alloc_mem -2*2*nx*ny
      write(*,'(6x,A,F9.2,1x,A)') 'deallocated',
     &           dble(2*2*nx*ny)/262144.D0,'MB'
#endif

! Note that the only outcome of the segment above is 4 integer numbers,
!                 iwest,jsouth, nx,ny
! which are the indices of south-west corner of the subdomain within
! the parent grid, and the sizes of the subdomain.  Everything else is
! discarded.  The next step is to allocate subdomain-sized 2D arrays
! for parent-grid variables: note that xp(nx,ny),yp(nx,ny) were
! deallocated above and now are allocated again with a different
! (expected to be smaller) size).

      allocate( srX(nx*ny), srY(nx*ny), sXY(nx*ny), sYX(nx*ny),
     &                                               stat=ierr )
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

! Re-read (read first time if CPP PARENT_GRID_SUBREGION is undefined)
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
      do j=1,ny
        do i=1,nx-1
          xpu(i,j)=0.5D0*(xp(i,j)+xp(i+1,j))
          ypu(i,j)=0.5D0*(yp(i,j)+yp(i+1,j))
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          xpv(i,j)=0.5D0*(xp(i,j)+xp(i,j+1))
          ypv(i,j)=0.5D0*(yp(i,j)+yp(i,j+1))
        enddo
      enddo

! Initialize horizontal interpolation:  The purpose of this stage is
!----------- ---------- --------------  to find (i) parent grid indices
! ip=ip(ic),jp=jp(ic) enclosing child-grid point "ic" into the parent
! grid element defined by 4 vertices, (ip,jp), (ip+1,jp), (ip,jp+1),
! (ip,jp+1)), and (ip+1,jp+1), and (ii) fractions xi,eta such that
! bi-linear interpolation of coordinates of the 4 vertices of parent
! grid yields coordinates (xc,yc) of the child grid.

! The semantic rules for variable names are:
!
!       role        placement           side
!      [i,j,xi,eta][r,u,v.ur,vr,uv.vu]_[west,east,south,north]
!
! allearing with all possible permutations (there are 4*7*4=102 arrays
! computed as sets of 4 in 28 calls to "bry_init_interp").
! Single-letter placement "r" indicates RHO-points of child grid
! interpolated from RHO-points of parent; similarly for "u" and "v".
! Dual-letter placements use first letter to identify location on
! the child grid and while second letter is type of source on parent.

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

#ifdef PARENT_GRID_SUBREGION
      imin=nx+1 ; imax=-1  !<-- initialize to unrealistic
      jmin=ny+1 ; jmax=-1  !<-- values outside the range
#endif

      if (OBC_WEST) then
        call bry_init_interp(nx,ny, xp,yp,  ncy,  xr_west,yr_west,
     &                           ir_west,jr_west, xi_west,eta_west)
#ifdef PARENT_GRID_SUBREGION
        call r2r_bry_index_bounds( ncy, ir_west,jr_west, mskr_west,
     &                                         imin,imax,jmin,jmax)
#endif
        call bry_init_interp(nx,ny, xp,yp,   ncy,  xu_west,yu_west,
     &                     iur_west,jur_west, xiur_west,etaur_west)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncy, xu_west,yu_west,
     &                          iu_west,ju_west, xiu_west,etau_west)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncy, xu_west,yu_west,
     &                      iuv_west,juv_west, xiuv_west,etauv_west)

        call bry_init_interp(nx,ny, xp,yp,  ncy-1,  xv_west,yv_west,
     &                      ivr_west,jvr_west, xivr_west,etavr_west)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncy-1,xv_west,yv_west,
     &                      ivu_west,jvu_west, xivu_west,etavu_west)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncy-1,xv_west,yv_west,
     &                          iv_west,jv_west, xiv_west,etav_west)
      endif

      if (OBC_EAST) then
        call bry_init_interp(nx,ny, xp,yp,  ncy,  xr_east,yr_east,
     &                           ir_east,jr_east, xi_east,eta_east)
#ifdef PARENT_GRID_SUBREGION
        call r2r_bry_index_bounds( ncy, ir_east,jr_east, mskr_east,
     &                                         imin,imax,jmin,jmax)
#endif
        call bry_init_interp(nx,ny,  xp,yp,   ncy, xu_east,yu_east,
     &                     iur_east,jur_east, xiur_east,etaur_east)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncy, xu_east,yu_east,
     &                          iu_east,ju_east, xiu_east,etau_east)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncy, xu_east,yu_east,
     &                      iuv_east,juv_east, xiuv_east,etauv_east)

        call bry_init_interp(nx,ny, xp,yp,   ncy-1, xv_east,yv_east,
     &                      ivr_east,jvr_east, xivr_east,etavr_east)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncy-1,xv_east,yv_east,
     &                      ivu_east,jvu_east, xivu_east,etavu_east)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncy-1,xv_east,yv_east,
     &                          iv_east,jv_east, xiv_east,etav_east)
      endif

      if (OBC_SOUTH) then
        call bry_init_interp(nx,ny,  xp,yp, ncx,  xr_south,yr_south,
     &                         ir_south,jr_south, xi_south,eta_south)
#ifdef PARENT_GRID_SUBREGION
        call r2r_bry_index_bounds( ncx, ir_south,jr_south, mskr_south,
     &                                            imin,imax,jmin,jmax)
#endif
        call bry_init_interp(nx,ny,  xp,yp, ncx-1, xu_south,yu_south,
     &                    iur_south,jur_south, xiur_south,etaur_south)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncx-1,xu_south,yu_south,
     &                        iu_south,ju_south, xiu_south,etau_south)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncx-1,xu_south,yu_south,
     &                    iuv_south,juv_south, xiuv_south,etauv_south)

        call bry_init_interp(nx,ny,  xp,yp,   ncx,  xv_south,yv_south,
     &                    ivr_south,jvr_south, xivr_south,etavr_south)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncx, xv_south,yv_south,
     &                    ivu_south,jvu_south, xivu_south,etavu_south)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncx,xv_south,yv_south,
     &                       iv_south,jv_south, xiv_south,etav_south)
      endif

      if (OBC_NORTH) then
        call bry_init_interp(nx,ny,  xp,yp,  ncx,  xr_north,yr_north,
     &                          ir_north,jr_north, xi_north,eta_north)
#ifdef PARENT_GRID_SUBREGION
        call r2r_bry_index_bounds( ncx, ir_north,jr_north, mskr_north,
     &                                            imin,imax,jmin,jmax)
#endif
        call bry_init_interp(nx,ny,  xp,yp,  ncx-1, xu_north,yu_north,
     &                    iur_north,jur_north, xiur_north,etaur_north)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncx-1,xu_north,yu_north,
     &                        iu_north,ju_north, xiu_north,etau_north)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncx-1,xu_north,yu_north,
     &                    iuv_north,juv_north, xiuv_north,etauv_north)

        call bry_init_interp(nx,ny,  xp,yp,   ncx,  xv_north,yv_north,
     &                    ivr_north,jvr_north, xivr_north,etavr_north)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncx, xv_north,yv_north,
     &                    ivu_north,jvu_north, xivu_north,etavu_north)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncx, xv_north,yv_north,
     &                        iv_north,jv_north, xiv_north,etav_north)
      endif

#ifdef PARENT_GRID_SUBREGION
      imax=imax+1 ; jmax=jmax+1
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
#endif
#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &                       'horizontal initialization complete in',
     &                           dble(inc_clk)/dble(clk_rate), 'sec'
#endif

! Setup (cosA,sinA) arrays to rotate velocity components at each open
!------ ----------- -----  boundary. Variable "angle" stored in netCDF
! grid files the angle between true East and local direction of ROMS
! XI-coordinate of the grid.  Read it then compute csA=cos(alpha) and
! snA=sin(alpha), first for child grid then for the parent. For the
! latter temporarily place the outcome into arrays xp,yp -- after all
! horizontal interpolations have been initialized above the content
! arrays is no longer needed.  Then interpolate cos(A) and sin(A) of
! the parent into child grid and compute cos and sin of the child-
! parent difference of angles; these will be used to rotate velocity
! vector components.

      call read_angle(nccgrd, chldgrd, 1,1,    ncx,ncy,   csA,snA)
      call read_angle(ncpgrd, prntgrd, iwest,jsouth, nx,ny, xp,yp)

      call spln2d_double(nx,ny, xp, srX,srY,sXY,sYX)
      if (OBC_WEST) then
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncy,1,
     &         iur_west,jur_west, xiur_west,etaur_west, csAu_west)
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncy-1,1,
     &         ivr_west,jvr_west, xivr_west,etavr_west, csAv_west)
      endif
      if (OBC_EAST) then
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncy,1,
     &         iur_east,jur_east, xiur_east,etaur_east, csAu_east)
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncy-1,1,
     &         ivr_east,jvr_east, xivr_east,etavr_east, csAv_east)
      endif
      if (OBC_SOUTH) then
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncx-1,1,
     &    iur_south,jur_south, xiur_south,etaur_south, csAu_south)
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncx,1,
     &    ivr_south,jvr_south, xivr_south,etavr_south, csAv_south)
      endif
      if (OBC_NORTH) then
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncx-1,1,
     &    iur_north,jur_north, xiur_north,etaur_north, csAu_north)
        call spln2d_interp_double(nx,ny, xp,  srX,srY,sXY, ncx,1,
     &    ivr_north,jvr_north, xivr_north,etavr_north, csAv_north)
      endif

      call spln2d_double(nx,ny, yp, srX,srY,sXY,sYX)
      if (OBC_WEST) then
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncy,1,
     &         iur_west,jur_west, xiur_west,etaur_west, snAu_west)
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncy-1,1,
     &         ivr_west,jvr_west, xivr_west,etavr_west, snAv_west)
      endif
      if (OBC_EAST) then
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncy,1,
     &         iur_east,jur_east, xiur_east,etaur_east, snAu_east)
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncy-1,1,
     &         ivr_east,jvr_east, xivr_east,etavr_east, snAv_east)
      endif
      if (OBC_SOUTH) then
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncx-1,1,
     &    iur_south,jur_south, xiur_south,etaur_south, snAu_south)
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncx,1,
     &    ivr_south,jvr_south, xivr_south,etavr_south, snAv_south)
      endif
      if (OBC_NORTH) then
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncx-1,1,
     &    iur_north,jur_north, xiur_north,etaur_north, snAu_north)
        call spln2d_interp_double(nx,ny, yp,  srX,srY,sXY, ncx,1,
     &    ivr_north,jvr_north, xivr_north,etavr_north, snAv_north)
      endif
                                                    ! Thus far arrays
      if (OBC_WEST) then                            ! csAu_west and
        do j=1,ncy                                  ! snAu_west contain
          csP=csAu_west(j) ; cs=csA(1,j)+csA(2,j)   ! cos(parent) and
          snP=snAu_west(j) ; sn=snA(1,j)+snA(2,j)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_west(j)=cff*(cs*csP+sn*snP)
          snAu_west(j)=cff*(sn*csP-cs*snP)          ! sin(parent)
        enddo                                       ! interpolated
        do j=1,ncy-1                                ! from parent
          csP=csAv_west(j) ; cs=csA(1,j)+csA(1,j+1) ! to child grid.
          snP=snAv_west(j) ; sn=snA(1,j)+snA(1,j+1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_west(j)=cff*(cs*csP+sn*snP)
          snAv_west(j)=cff*(sn*csP-cs*snP)          ! Convert them into
        enddo                                       ! cos(child-parent)
      endif                                         !       and
                                                    ! sin(child-parent)
      if (OBC_EAST) then
        do j=1,ncy
          csP=csAu_east(j) ;  cs=csA(ncx,j)+csA(ncx-1,j)
          snP=snAu_east(j) ;  sn=snA(ncx,j)+snA(ncx-1,j)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_east(j)=cff*(cs*csP+sn*snP)
          snAu_east(j)=cff*(sn*csP-cs*snP)
        enddo
        do j=1,ncy-1
          csP=csAv_east(j) ;  cs=csA(ncx,j)+csA(ncx,j+1)
          snP=snAv_east(j) ;  sn=snA(ncx,j)+snA(ncx,j+1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_east(j)=cff*(cs*csP+sn*snP)
          snAv_east(j)=cff*(sn*csP-cs*snP)
        enddo
      endif

      if (OBC_SOUTH) then
        do i=1,ncx
          csP=csAv_south(i) ; cs=csA(i,1)+csA(i,2)
          snP=snAv_south(i) ; sn=snA(i,1)+snA(i,2)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_south(i)=cff*(cs*csP+sn*snP)
          snAv_south(i)=cff*(sn*csP-cs*snP)
        enddo
        do i=1,ncx-1
          csP=csAu_south(i) ; cs=csA(i,1)+csA(i+1,1)
          snP=snAu_south(i) ; sn=snA(i,1)+snA(i+1,1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_south(i)=cff*(cs*csP+sn*snP)
          snAu_south(i)=cff*(sn*csP-cs*snP)
        enddo
      endif

      if (OBC_NORTH) then
        do i=1,ncx
          csP=csAv_north(i) ; cs=csA(i,ncy)+csA(i,ncy-1)
          snP=snAv_north(i) ; sn=snA(i,ncy)+snA(i,ncy-1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_north(i)=cff*(cs*csP+sn*snP)
          snAv_north(i)=cff*(sn*csP-cs*snP)
        enddo
        do i=1,ncx-1
          csP=csAu_north(i) ; cs=csA(i,ncy)+csA(i+1,ncy)
          snP=snAu_north(i) ; sn=snA(i,ncy)+snA(i+1,ncy)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_north(i)=cff*(cs*csP+sn*snP)
          snAu_north(i)=cff*(sn*csP-cs*snP)
        enddo
      endif

! Read parent-grid land mask. Similarly to the above: temporarily
! use array "hp" to read "mask_rho" from the file, thereafter this
! program uses only integer(kind=2) version of mask.

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

! Read parent-grid topography and interpolate it to child grid, then
! initialize vertical interpolation for each line along each boundary
! (there are 4 boundaries, each having 3 lines: RHO-, U-, and V-points
! all of which have distinct type of vertical interpolations, so the
! ultimate goal is to initialize 12 arrays, kpX_size, where X={r,y,v}
! and _side={west,east,north,south}). Note that it is possible that
! some portion of child grid is outside the parent, so in that areas
! horizontal interpolation results to zero values.  Even though in any
! meaningful model configuration these areas are under land mask, and
! therefore make no effect, leaving zero-valued "hp" there confuses
! vertical search algorithm inside "bry_init_vertinterp" resulting
! in out-of bound array index.  For performance reason it is better
! not to add extra protective logic there, but rather fill-in these
! areas with child-grid topography.

      call get_var_by_name_double(nccgrd, 'h',  h)
      call get_patch_by_name_double(ncpgrd, prntgrd, 'h',
     &                      iwest,jsouth, nx,ny,0,0, hp)

      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)

      if (OBC_WEST) then
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncy,1,
     &                 ir_west,jr_west, xi_west,eta_west, hp_west)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncy,1,
     &          iur_west,jur_west, xiur_west,etaur_west, hpu_west)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncy-1,1,
     &          ivr_west,jvr_west, xivr_west,etavr_west, hpv_west)
        do j=1,ncy
          h_west(j)=h(1,j)
          if (hp_west(j) < 0.0001D0)  hp_west(j)=h_west(j)
          hu_west(j)=0.5D0*(h(1,j)+h(2,j))
          if (hpu_west(j) < 0.0001D0) hpu_west(j)=hu_west(j)
        enddo
        do j=1,ncy-1
          hv_west(j)=0.5D0*(h(1,j)+h(1,j+1))
          if (hpv_west(j) < 0.0001D0) hpv_west(j)=hv_west(j)
        enddo
        call bry_init_vertinterp(ncy, hp_west, Np,hcp,Csp_r,
     &                            h_west, N,hc,Cs_r, kpr_west)
        call bry_init_vertinterp(ncy, hpu_west, Np,hcp,Csp_r,
     &                           hu_west, N,hc,Cs_r, kpu_west)
        call bry_init_vertinterp(ncy-1, hpv_west, Np,hcp,Csp_r,
     &                           hv_west, N,hc,Cs_r, kpv_west)
      endif

      if (OBC_EAST) then
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncy,1,
     &                 ir_east,jr_east, xi_east,eta_east, hp_east)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncy,1,
     &          iur_east,jur_east, xiur_east,etaur_east, hpu_east)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncy-1,1,
     &          ivr_east,jvr_east, xivr_east,etavr_east, hpv_east)
        do j=1,ncy
          h_east(j)=h(ncx,j)
          hu_east(j)=0.5D0*( h(ncx,j)+h(ncx-1,j))
          if (hp_east(j) < 0.0001D0)  hp_east(j)=h_east(j)
          hu_east(j)=0.5D0*( h(ncx,j)+h(ncx-1,j))
          if (hpu_east(j) < 0.0001D0) hpu_east(j)=hu_east(j)
        enddo
        do j=1,ncy-1
          hv_east(j)=0.5D0*(h(ncx,j)+h(ncx,j+1))
          if (hpv_east(j) < 0.0001D0) hpv_east(j)=hv_east(j)
        enddo
        call bry_init_vertinterp(ncy, hp_east, Np,hcp,Csp_r,
     &                            h_east, N,hc,Cs_r, kpr_east)
        call bry_init_vertinterp(ncy, hpu_east, Np,hcp,Csp_r,
     &                           hu_east, N,hc,Cs_r, kpu_east)
        call bry_init_vertinterp(ncy-1, hpv_east, Np,hcp,Csp_r,
     &                           hv_east, N,hc,Cs_r, kpv_east)
      endif

      if (OBC_SOUTH) then
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncx,1,
     &            ir_south,jr_south, xi_south,eta_south, hp_south)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncx-1,1,
     &     iur_south,jur_south, xiur_south,etaur_south, hpu_south)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncx,1,
     &     ivr_south,jvr_south, xivr_south,etavr_south, hpv_south)
        do i=1,ncx
          h_south(i)=h(i,1)
          if (hp_south(i) < 0.0001D0)  hp_south(i)=h_south(i)
          hv_south(i)=0.5D0*(h(i,1)+h(i,2))
          if (hpv_south(i) < 0.0001D0) hpv_south(i)=hv_south(i)
        enddo
        do i=1,ncx-1
          hu_south(i)=0.5D0*(h(i,1)+h(i+1,1))
          if (hpu_south(i) < 0.0001D0) hpu_south(i)=hu_south(i)
        enddo
        call bry_init_vertinterp(ncx, hp_south, Np,hcp,Csp_r,
     &                            h_south, N,hc,Cs_r, kpr_south)
        call bry_init_vertinterp(ncx-1, hpu_south, Np,hcp,Csp_r,
     &                           hu_south, N,hc,Cs_r, kpu_south)
        call bry_init_vertinterp(ncx, hpv_south, Np,hcp,Csp_r,
     &                           hv_south, N,hc,Cs_r, kpv_south)
      endif

      if (OBC_NORTH) then
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncx,1,
     &            ir_north,jr_north, xi_north,eta_north, hp_north)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncx-1,1,
     &     iur_north,jur_north, xiur_north,etaur_north, hpu_north)
        call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY, ncx,1,
     &     ivr_north,jvr_north, xivr_north,etavr_north, hpv_north)
        do i=1,ncx
          h_north(i)=h(i,ncy)
          if (hp_north(i) < 0.0001D0)  hp_north(i)=h_north(i)
          hv_north(i)=0.5D0*(h(i,ncy) +h(i,ncy-1))
          if (hpv_north(i) < 0.0001D0) hpv_north(i)=hv_north(i)
        enddo
        do i=1,ncx-1
          hu_north(i)=0.5D0*(h(i,ncy)+h(i+1,ncy))
          if (hpu_north(i) < 0.0001D0) hpu_north(i)=hu_north(i)
        enddo
        call bry_init_vertinterp(ncx, hp_north, Np,hcp,Csp_r,
     &                            h_north, N,hc,Cs_r, kpr_north)
        call bry_init_vertinterp(ncx-1, hpu_north, Np,hcp,Csp_r,
     &                           hu_north, N,hc,Cs_r, kpu_north)
        call bry_init_vertinterp(ncx, hpv_north, Np,hcp,Csp_r,
     &                           hv_north, N,hc,Cs_r, kpv_north)
      endif

      write(*,'(/1x,A)') 'Checking...'
      if (OBC_WEST) then
        write(*,'(2x,A)') 'western'
        call bry_check_init_vertinterp(ncy, hp_west, Np,hcp,Csp_r,
     &                                 h_west, N,hc,Cs_r, kpr_west)
        call bry_check_init_vertinterp(ncy, hpu_west, Np,hcp,Csp_r,
     &                                 hu_west, N,hc,Cs_r, kpu_west)
        call bry_check_init_vertinterp(ncy-1, hpv_west, Np,hcp,Csp_r,
     &                                 hv_west, N,hc,Cs_r, kpv_west)
      endif
      if (OBC_EAST) then
        write(*,'(2x,A)') 'eastern'
        call bry_check_init_vertinterp(ncy, hp_east, Np,hcp,Csp_r,
     &                                 h_east, N,hc,Cs_r, kpr_east)
        call bry_check_init_vertinterp(ncy, hpu_east, Np,hcp,Csp_r,
     &                                 hu_east, N,hc,Cs_r, kpu_east)
        call bry_check_init_vertinterp(ncy-1, hpv_east, Np,hcp,Csp_r,
     &                                 hv_east, N,hc,Cs_r, kpv_east)
      endif
      if (OBC_SOUTH) then
        write(*,'(2x,A)') 'southern'
        call bry_check_init_vertinterp(ncx, hp_south, Np,hcp,Csp_r,
     &                                h_south, N,hc,Cs_r, kpr_south)
        call bry_check_init_vertinterp(ncx-1, hpu_south, Np,hcp,Csp_r,
     &                                hu_south, N,hc,Cs_r, kpu_south)
        call bry_check_init_vertinterp(ncx, hpv_south, Np,hcp,Csp_r,
     &                                hv_south, N,hc,Cs_r, kpv_south)
      endif
      if (OBC_NORTH) then
        write(*,'(2x,A)') 'northern'
        call bry_check_init_vertinterp(ncx, hp_north, Np,hcp,Csp_r,
     &                                h_north, N,hc,Cs_r, kpr_north)
        call bry_check_init_vertinterp(ncx-1, hpu_north, Np,hcp,Csp_r,
     &                                hu_north, N,hc,Cs_r, kpu_north)
        call bry_check_init_vertinterp(ncx, hpv_north, Np,hcp,Csp_r,
     &                              hv_north, N,hc,Cs_r, kpv_north)
      endif
#ifdef TIMING
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(/F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &                          'vertical initialization complete in',
     &                            dble(inc_clk)/dble(clk_rate), 'sec'
#endif

! Create boundary forcing file.

      ierr=nf_create(roms_bry, nf_64bit_offset, nctarg)
      if (ierr == nf_noerr) then
        call def_bry_file(roms_bry, nctarg, ncx,ncy,N,
     &         OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH,
     &                 theta_s, theta_b, hc, Cs_w,Cs_r)
        ierr=nf_inq_varid(nctarg, 'bry_time', tvar_out)
        if (ierr /= nf_noerr) then
          write(*,*) '### ERROR 2: nf_inq_varid(. tvar_out .)'
        endif
      else
        write(*,*) '### ERROR 1: nf_create'
      endif
      if (ierr /= nf_noerr) stop
      ierr=nf_enddef(nctarg)  !<-- set to input mode
      ierr=nf_sync(nctarg)


! Finally clean up memory which is no longer needed and allocate
! big buffer arrays to read the actual data.


      deallocate(hp, xp,yp, xpu,ypu, xpv,ypv, xc,yc)
      net_alloc_mem=net_alloc_mem -(6*2*nx*ny+2*2*ncx*ncy)
      write(*,'(8x,A,F9.2,1x,A)') 'deallocated',
     &           dble(6*2*nx*ny+2*2*ncx*ncy)/262144.D0,'MB'


      size=max(ncx,ncy)*max(Np,N)
      allocate(wrk3(size), wrk4(size), wrk5(size), wrk6(size),
     &                        wrk7(size), wrk8(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 6*size
        write(*,'(8x,A,F10.2,1x,A)')  'allocated wrk3,3.4, reaching',
     &                     dble(net_alloc_mem)/262144.D0, 'MB total'
        size=nx*ny*(Np+1)
        allocate(wrk1(size), wrk2(size), stat=ierr)
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem + 2*size
          write(*,'(8x,A,F10.2,1x,A)')    'allocated wrk1,2, reaching',
     &                       dble(net_alloc_mem)/262144.D0, 'MB total'
        endif
      endif
      write(*,'(2/7(/4x,A)/)')

     &'      *****    *********    ******   *******    *********  ',
     &'    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *  ',
     &'    ***           ***     **   ***   ***   ***     ***     ',
     &'      *****       ***    ***   ***   ***   **      ***     ',
     &'          ***     ***    *********   ******        ***     ',
     &'    ***   ***     ***    ***   ***   ***  **       ***     ',
     &'      *****       ***    ***   ***   ***   ***     ***     '


! At this moment the first history file is already open, while all
! others need to be open when their time comes.  Determine the number
! of records in the file, netCDF ID and units for timing variable (time
! in seconds needs to be converted ito days).

      do ihis=1,nhists
        if (ihis > 1) then
          lprnt=lenstr(prnt_data(ihis))
          ierr=nf_open(prnt_data(ihis)(1:lprnt), nf_nowrite, ncsrc)
        endif
        ierr=nf_inq_varid(ncsrc, time_var_name, tvar_in)
        if (ierr == nf_noerr) then
          ierr=nf_inq_var(ncsrc, tvar_in, vname,   vtype,
     &                           ndims,   vdimids, natts)
          if (ierr == nf_noerr) then
            lvar=lenstr(vname)
            if (ndims == 1) then
              ierr=nf_inq_dimlen(ncsrc, vdimids, nrecs)
              if (ierr == nf_noerr) then
#ifdef VERBOSE
                write(*,'(4(2x,A,I3),3I3)')  'vtype =', vtype,
     &                     'natts =', natts, 'ndims =', ndims,
     &                     'vdimids =', (vdimids(i),i=1,ndims)
#endif
                do i=1,natts
                  ierr=nf_inq_attname(ncsrc, tvar_in, i, str)
                  if (ierr == nf_noerr) then
                    lstr=lenstr(str)
                    if (str(1:lenstr(str)) == 'units') then
                      ierr=nf_get_att_text(ncsrc,tvar_in, 'units',str)
                      if (ierr == nf_noerr) then
                        lstr=lenstr(str) ; time_units=str(1:lstr)
                      else
                        write(*,'(/1x,7A/12x,A/)')     '### ERROR: ',
     &                 'Found, but cannot read attribute ''units'' ',
     &                 'for variable ''', vname(1:lvar),  ''' in ''',
     &                               prnt_data(ihis)(1:lprnt), '''.',
     &                                             nf_strerror(ierr)
                        stop
                      endif
                    endif
                  else
                    write(*,'(/1x,2A,I3,1x,5A/12x,A/)') '### ERROR: ',
     &                       'Cannot make inquiry for attribute #', i,
     &                 'for variable ''',  vname(1:lvar),  ''' in ''',
     &              prnt_data(ihis)(1:lprnt),'''.', nf_strerror(ierr)
                    stop
                  endif
                enddo
              else
                write(*,'(1x,4A/12x,A/)')    '### ERROR: Cannot ',
     &                   'determine number of time records in ''',
     &          prnt_data(ihis)(1:lprnt),'''.', nf_strerror(ierr)
              endif
            elseif (ndims == 0) then
              write(*,'(2(/1x,3A)/)')         'Time variable ''',
     &               vname(1:lvar), ''' does not have dimensions.',
     &               'Presuming that ''', prnt_data(ihis)(1:lprnt),
     &                               ''' is a single-record file.'
              nrecs=1
            else
              write(*,'(/1x,5A/)')   '### ERROR: Time variable ''',
     &         vname(1:lvar), ''' in ''', prnt_data(ihis)(1:lprnt),
     &                           ''' has more than one dimension.'
            endif
          else
            write(*,'(/1x,6A/12x,A/)')  '### ERROR: Cannot make ',
     &       'inquiry for variable ''', vname(1:lvar), ''' in ''',
     &         prnt_data(ihis)(1:lprnt), '''.', nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,6A/12x,A/)')  '### ERROR: Cannot find ',
     &                 'variable ''', vname(1:lvar), ''' in ''',
     &       prnt_data(ihis)(1:lprnt), '''.', nf_strerror(ierr)
        endif
        if (ierr == nf_noerr) then
          lstr=lenstr(time_units)
          if (lstr > 0) then
            write(*,'(/1x,A,I5,1x,6A)') 'Found', nrecs, 'time ',
     &        'records in ''', prnt_data(ihis)(1:lprnt), ''', ',
     &                      'time units = ', time_units(1:lstr)
          else
            write(*,'(/1x,5A/)')  '### ERROR: Time variable ''',
     &      vname(1:lvar), ''' in ''', prnt_data(ihis)(1:lprnt),
     &                  ''' does not have attribute ''units''.'
            stop
          endif
        else
          stop
        endif

        do rec=1,nrecs
          recout=recout+1
#ifdef TIMING
          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate,clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          if (inc_clk < 0) inc_clk=inc_clk+clk_max
          net_clk=net_clk+inc_clk
          write(*,'(/F10.2,1x,A,I6/)') dble(net_clk)/dble(clk_rate),
     &                            'sec  Processing record ', recout
#endif

! Time
          ierr=nf_get_vara_double(ncsrc, tvar_in, rec,1, time)
          if (ierr == nf_noerr) then
            if (time_units(1:6) == 'second') time=time/86400.D0
            ierr=nf_put_vara_double(nctarg, tvar_out, recout,1, time)
            if (ierr == nf_noerr) then
              write(*,*) '    bry_time =', time, ' days'
            else
              write(*,'(/1x,3A/)')  '### ERROR: Cannot write ''',
     &                         'bry_time'', ', nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,5A/)')     '### ERROR: Cannot read ''',
     &                                time_var_name,''' from ''',
     &      prnt_data(ihis)(1:lprnt), ''' ,', nf_strerror(ierr)
          endif

! Free surface ...

          call get_patch_by_name_real(ncsrc, prnt_data(ihis),'zeta',
     &                             iwest,jsouth, nx,ny,0,rec, wrk1)
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1)
          call etch_into_land_thread(nx,ny, mskp, wrk1)
C$OMP END PARALLEL
          call spln2d_real(nx,ny, wrk1, srX,srY,sXY,sYX)

          if (OBC_WEST) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncy,1,
     &                  ir_west,jr_west,   xi_west,eta_west,   wrk2)
            call bry_apply_mask(ncy,1, mskr_west, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_west',
     &                                         ncy,0,0,recout, wrk2)
          endif
          if (OBC_EAST) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncy,1,
     &                  ir_east,jr_east,   xi_east,eta_east,   wrk2)
            call bry_apply_mask(ncy,1, mskr_east, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_east',
     &                                         ncy,0,0,recout, wrk2)
          endif
          if (OBC_SOUTH) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncx,1,
     &                  ir_south,jr_south, xi_south,eta_south, wrk2)
            call bry_apply_mask(ncx,1, mskr_south, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_south',
     &                                         ncx,0,0,recout, wrk2)
          endif
          if (OBC_NORTH) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncx,1,
     &                  ir_north,jr_north, xi_north,eta_north, wrk2)
            call bry_apply_mask(ncx,1, mskr_north, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_north',
     &                                         ncx,0,0,recout, wrk2)
          endif

! Horizontal velocities, both u,v and barotropic...

          call get_patch_by_name_real(ncsrc, prnt_data(ihis), 'u',
     &                        iwest,jsouth, nx-1,ny,Np,rec, wrk1)
          call get_patch_by_name_real(ncsrc, prnt_data(ihis), 'v',
     &                        iwest,jsouth, nx,ny-1,Np,rec, wrk2)
          do k=1,Np
            isrc=1+(k-1)*(nx-1)*ny
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1, isrc)
            call etch_into_land_thread(nx-1,ny, umsp, wrk1(isrc))
C$OMP END PARALLEL
          enddo
          do k=1,Np
            isrc=1+(k-1)*nx*(ny-1)
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1, isrc)
            call etch_into_land_thread(nx,ny-1, vmsp, wrk2(isrc))
C$OMP END PARALLEL
          enddo

          if (OBC_WEST) then
            call bry_interp(nx-1,ny,Np, wrk1, ncy,  iu_west,ju_west,
     &                           xiu_west,etau_west, msku_west, wrk3)
            call bry_interp(nx,ny-1,Np, wrk2, ncy, iuv_west,juv_west,
     &                         xiuv_west,etauv_west, msku_west, wrk4)
            call bry_rotate_u_in_place(ncy,Np, csAu_west,snAu_west,
     &                                                     wrk3,wrk4)
            call bry_vertinterp(ncy, btm_slp, Np, wrk3,  N, kpu_west,
     &                                                          wrk5)
            call bry_vert_average(ncy, hu_west, N,hc,Cs_w, wrk5,wrk7)
            call bry_apply_mask(ncy,N, msku_west, wrk5)
            call bry_apply_mask(ncy,1, msku_west, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_west',
     &                                       ncy,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry, 'ubar_west',
     &                                       ncy,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np, wrk2, ncy-1, iv_west,jv_west,
     &                            xiv_west,etav_west, mskv_west, wrk4)
            call bry_interp(nx-1,ny,Np, wrk1, ncy-1, ivu_west,jvu_west,
     &                          xivu_west,etavu_west, mskv_west, wrk3)
            call bry_rotate_v_in_place(ncy-1,Np, csAv_west,snAv_west,
     &                                                      wrk3,wrk4)
            call bry_vertinterp(ncy-1,btm_slp, Np, wrk4,  N, kpv_west,
     &                                                           wrk6)
            call bry_vert_average(ncy-1,hv_west, N,hc,Cs_w, wrk6,wrk8)
            call bry_apply_mask(ncy-1,N, mskv_west, wrk6)
            call bry_apply_mask(ncy-1,1, mskv_west, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_west',
     &                                     ncy-1,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry, 'vbar_west',
     &                                     ncy-1,0,0,recout, wrk8)
          endif

          if (OBC_EAST) then
            call bry_interp(nx-1,ny,Np, wrk1, ncy,  iu_east,ju_east,
     &                           xiu_east,etau_east, msku_east, wrk3)
            call bry_interp(nx,ny-1,Np, wrk2, ncy, iuv_east,juv_east,
     &                         xiuv_east,etauv_east, msku_east, wrk4)
            call bry_rotate_u_in_place(ncy,Np, csAu_east,snAu_east,
     &                                                     wrk3,wrk4)
            call bry_vertinterp(ncy, btm_slp, Np,wrk3, N,kpu_east,
     &                                                          wrk5)
            call bry_vert_average(ncy, hu_east, N,hc,Cs_w, wrk5,wrk7)
            call bry_apply_mask(ncy,N, msku_east, wrk5)
            call bry_apply_mask(ncy,1, msku_east, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_east',
     &                                       ncy,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry, 'ubar_east',
     &                                       ncy,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np,wrk2, ncy-1, iv_east,jv_east,
     &                           xiv_east,etav_east, mskv_east, wrk4)
            call bry_interp(nx-1,ny,Np,wrk1, ncy-1, ivu_east,jvu_east,
     &                         xivu_east,etavu_east, mskv_east, wrk3)
            call bry_rotate_v_in_place(ncy-1,Np, csAv_east,snAv_east,
     &                                                     wrk3,wrk4)
            call bry_vertinterp(ncy-1, btm_slp, Np,wrk4, N,kpv_east,
     &                                                          wrk6)
            call bry_vert_average(ncy-1,hv_east, N,hc,Cs_w,wrk6,wrk8)
            call bry_apply_mask(ncy-1,N, mskv_east, wrk6)
            call bry_apply_mask(ncy-1,1, mskv_east, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_east',
     &                                   ncy-1,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry, 'vbar_east',
     &                                   ncy-1,0,0,recout, wrk8)
          endif

          if (OBC_SOUTH) then
            call bry_interp(nx-1,ny,Np,wrk1, ncx-1,iu_south,ju_south,
     &                         xiu_south,etau_south, msku_south, wrk3)
            call bry_interp(nx,ny-1,Np,wrk2, ncx-1,iuv_south,juv_south,
     &                       xiuv_south,etauv_south, msku_south, wrk4)
            call bry_rotate_u_in_place(ncx-1,Np, csAu_south,snAu_south,
     &                                                      wrk3,wrk4)
            call bry_vertinterp(ncx-1,btm_slp, Np,wrk3, N,kpu_south,
     &                                                           wrk5)
            call bry_vert_average(ncx-1,hu_south, N,hc,Cs_w,wrk5,wrk7)
            call bry_apply_mask(ncx-1,N, msku_south, wrk5)
            call bry_apply_mask(ncx-1,1, msku_south, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_south',
     &                                      ncx-1,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry,'ubar_south',
     &                                      ncx-1,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np, wrk2, ncx, iv_south,jv_south,
     &                         xiv_south,etav_south, mskv_south, wrk4)
            call bry_interp(nx-1,ny,Np, wrk1, ncx, ivu_south,jvu_south,
     &                       xivu_south,etavu_south, mskv_south, wrk3)
            call bry_rotate_v_in_place(ncx,Np, csAv_south,snAv_south,
     &                                                      wrk3,wrk4)
            call bry_vertinterp(ncx, btm_slp, Np,wrk4, N,kpv_south,
     &                                                           wrk6)
            call bry_vert_average(ncx, hv_south, N,hc,Cs_w, wrk6,wrk8)
            call bry_apply_mask(ncx,N, mskv_south, wrk6)
            call bry_apply_mask(ncx,1, mskv_south, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_south',
     &                                       ncx,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry,'vbar_south',
     &                                       ncx,0,0,recout, wrk8)
          endif

          if (OBC_NORTH) then
            call bry_interp(nx-1,ny,Np,wrk1, ncx-1, iu_north,ju_north,
     &                         xiu_north,etau_north, msku_north, wrk3)
            call bry_interp(nx,ny-1,Np,wrk2, ncx-1,iuv_north,juv_north,
     &                       xiuv_north,etauv_north, msku_north, wrk4)
            call bry_rotate_u_in_place(ncx-1,Np, csAu_north,snAu_north,
     &                                                      wrk3,wrk4)
            call bry_vertinterp(ncx-1, btm_slp, Np, wrk3, N,kpu_north,
     &                                                           wrk5)
            call bry_vert_average(ncx-1,hu_north, N,hc,Cs_w,wrk5,wrk7)
            call bry_apply_mask(ncx-1,N, msku_north, wrk5)
            call bry_apply_mask(ncx-1,1, msku_north, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_north',
     &                                      ncx-1,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry, 'ubar_north',
     &                                      ncx-1,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np, wrk2, ncx, iv_north,jv_north,
     &                         xiv_north,etav_north, mskv_north, wrk4)
            call bry_interp(nx-1,ny,Np, wrk1, ncx, ivu_north,jvu_north,
     &                       xivu_north,etavu_north, mskv_north, wrk3)
            call bry_rotate_v_in_place(ncx,Np, csAv_north,snAv_north,
     &                                                      wrk3,wrk4)
            call bry_vertinterp(ncx, btm_slp, Np,wrk4, N,kpv_north,
     &                                                           wrk6)
            call bry_vert_average(ncx, hv_north, N,hc,Cs_w, wrk6,wrk8)
            call bry_apply_mask(ncx,N, mskv_north, wrk6)
            call bry_apply_mask(ncx,1, mskv_north, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_north',
     &                                      ncx,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry,'vbar_north',
     &                                      ncx,0,0,recout, wrk8)
          endif

! Tracers...

          do itrc=1,ntrc
            vname=tracer(itrc) ; lvar=lenstr(vname)
            call get_patch_by_name_real(ncsrc, prnt_data(ihis),vname,
     &                              iwest,jsouth, nx,ny,Np,rec, wrk1)
            do k=1,Np
              isrc=1+(k-1)*nx*ny
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1, isrc)
              call etch_into_land_thread(nx,ny, mskp, wrk1(isrc))
C$OMP END PARALLEL
            enddo

            if (OBC_WEST) then
              trgname=vname(1:lvar)/ /'_west'; ltgv=lenstr(trgname)
              call bry_interp(nx,ny,Np, wrk1,  ncy,  ir_west,jr_west,
     &                            xi_west,eta_west, mskr_west, wrk3)
              call bry_vertinterp(ncy, btm_trc, Np,wrk3, N,kpr_west,
     &                                                         wrk5)
              call bry_apply_mask(ncy,N, mskr_west, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     &                                      0,ncy,N,recout, wrk5)
            endif

            if (OBC_EAST) then
              trgname=vname(1:lvar)/ /'_east'; ltgv=lenstr(trgname)
              call bry_interp(nx,ny,Np, wrk1,  ncy,  ir_east,jr_east,
     &                            xi_east,eta_east, mskr_east, wrk3)
              call bry_vertinterp(ncy, btm_trc, Np,wrk3, N,kpr_east,
     &                                                         wrk5)
              call bry_apply_mask(ncy,N, mskr_east, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     &                                      0,ncy,N,recout, wrk5)
            endif

            if (OBC_SOUTH) then
              trgname=vname(1:lvar)/ /'_south'; ltgv=lenstr(trgname)
              call bry_interp(nx,ny,Np, wrk1,  ncx, ir_south,jr_south,
     &                          xi_south,eta_south, mskr_south, wrk3)
              call bry_vertinterp(ncx, btm_trc, Np,wrk3, N,kpr_south,
     &                                                          wrk5)
              call bry_apply_mask(ncx,N, mskr_south, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     &                                       ncx,0,N,recout, wrk5)
            endif

            if (OBC_NORTH) then
              trgname=vname(1:lvar)/ /'_north'; ltgv=lenstr(trgname)
              call bry_interp(nx,ny,Np, wrk1,  ncx, ir_north,jr_north,
     &                          xi_north,eta_north, mskr_north, wrk3)
              call bry_vertinterp(ncx, btm_trc, Np,wrk3, N,kpr_north,
     &                                                          wrk5)
              call bry_apply_mask(ncx,N, mskr_north, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     &                                      ncx,0,N,recout, wrk5)
            endif
          enddo  !<-- itrc        ! forcefully sync target file once
        enddo  !<-- rec           ! in a while to make it readable if
        ierr=nf_close(ncsrc)
        if (mod(recout,32) == 0) ierr=nf_sync(nctarg)
      enddo  !<-- ihis
      ierr=nf_close(nctarg)       ! the program is interrupted.
#ifdef TIMING
      write(*,'(/4x,A,F12.3,1x,A,F10.1,1x,A)')    'total data read ',
     &                  dble(sz_read_acc)/dble(1024**2), 'MBytes in',
     &                          dble(read_clk)/dble(clk_rate), 'sec'
      write(*,'(1x,A,F12.3,1x,A,F10.1,1x,A/)') 'total data written ',
     &                 dble(sz_write_acc)/dble(1024**2), 'MBytes in',
     &                         dble(write_clk)/dble(clk_rate), 'sec'
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      call cpu_time(tend)
      write(*,'(/1x,A,F8.2,1x,A,4x,A,F8.2,1x,A,F8.1,1x,A)')
     &   'Wall Clock time:', dble(net_clk)/dble(clk_rate), 'sec',
     &                           'CPU time:', tend-tstart, 'sec',
     & (tend-tstart)*dble(clk_rate)/dble(net_clk)*100.D0,'% CPU'
#endif
      end
