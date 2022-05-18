
      module ocean_vars
      ! ocean variables and output routines for those variables

#include "cppdefs.opt"

      use param
      use netcdf
      use nc_read_write
      use roms_read_write
      use scalars
      implicit none
      private

#include "ocean_vars.opt"

      ! module preamble:  ![
      ! 2D - taken from what was ocean2d.F
      real,public,allocatable,dimension(:,:,:) :: zeta            ! free surface elevation [m] and barotropic
      real,public,allocatable,dimension(:,:,:) :: ubar            ! velocity components in XI-directions
      real,public,allocatable,dimension(:,:,:) :: vbar
      real,public,allocatable,dimension(:,:) :: zeta_avg
      real,public,allocatable,dimension(:,:) :: ubar_avg
      real,public,allocatable,dimension(:,:) :: vbar_avg
#ifdef NHMG
      real,public,allocatable,dimension(:,:) :: nh_ubar
      real,public,allocatable,dimension(:,:) :: nh_vbar
      real,public,allocatable,dimension(:,:) :: nh_wcor
#endif

      ! 3D - taken from what was ocean2d.F
#ifdef SOLVE3D
      real,public,allocatable,dimension(:,:,:,:) :: u
      real,public,allocatable,dimension(:,:,:,:) :: v
      real,public,allocatable,dimension(:,:,:)   :: u_avg
      real,public,allocatable,dimension(:,:,:)   :: v_avg
# if defined NHMG
      real,public,allocatable,dimension(:,:,:,:) :: w
      real,public,allocatable,dimension(:,:,:,:) :: nhdu
      real,public,allocatable,dimension(:,:,:,:) :: nhdv
      real,public,allocatable,dimension(:,:,:,:) :: nhdw
# endif

      real,public,allocatable,dimension(:,:,:) :: FlxU
      real,public,allocatable,dimension(:,:,:) :: FlxV
      real,public,allocatable,dimension(:,:,:) :: We      ! explicit
      real,public,allocatable,dimension(:,:,:) :: Wi      ! implicit

      real,public,allocatable,dimension(:,:,:) :: w_avg
      real,public,allocatable,dimension(:,:,:) :: wvl_avg

      real,public,allocatable,dimension(:,:,:) :: Hz      ! height of rho-cell
      real,public,allocatable,dimension(:,:,:) :: z_r     ! depth at rho-points
      real,public,allocatable,dimension(:,:,:) :: z_w     ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      real,public,allocatable,dimension(:,:,:) :: dzdxi
      real,public,allocatable,dimension(:,:,:) :: dzdeta
# endif
#endif  /* SOLVE3D */

      ! netcdf outputting:
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_ovars=0
      integer :: navg_ovars = 0                              ! number of samples in average

      public :: init_arrays_ocean2d
      public :: init_arrays_ocean3d
      public :: wrt_rst
      public :: wrt_his_ocean_vars
      public :: wrt_avg_ocean_vars
      public :: calc_avg_ocean_vars

      contains  !]

!----------------------------------------------------------------------
      subroutine init_arrays_ocean2d  ![
      use ncvars
      use scalars

      implicit none

      allocate( zeta(GLOBAL_2D_ARRAY,4) ); zeta=0.         ! zeta(:,:,knew) needs to be =0. for set_depth_tile
      allocate( ubar(GLOBAL_2D_ARRAY,4) ); ubar=0.         ! since knew can change if exact restart or not
      allocate( vbar(GLOBAL_2D_ARRAY,4) ); vbar=0.         ! set all zeta = 0.

      if (wrt_file_avg) then
        if (wrt_avg_Z)  then
          allocate( zeta_avg(GLOBAL_2D_ARRAY) )
          zeta_avg=0.                                        ! avg needs to be 0. because Nan x 0. = Nan in set_avg.F
        endif
        if (wrt_avg_Ub) then
          allocate( ubar_avg(GLOBAL_2D_ARRAY) )
          ubar_avg=0.
        endif
        if (wrt_avg_Vb) then
          allocate( vbar_avg(GLOBAL_2D_ARRAY) )
          vbar_avg=0.
        endif
      endif

#ifdef NHMG
      allocate( nh_ubar(GLOBAL_2D_ARRAY) )
      allocate( nh_vbar(GLOBAL_2D_ARRAY) )
      allocate( nh_wcor(GLOBAL_2D_ARRAY) )
#endif

      ! NOTE: call this here to prevent output in between timesteps which would break code
      ! checking script.
      if (wrt_file_his) call display_output_settings_to_terminal
      if (wrt_file_avg) call display_output_settings_to_terminal_avg

      end subroutine init_arrays_ocean2d  !]

!----------------------------------------------------------------------
      subroutine init_arrays_ocean3d  ![
      use ncvars
      use scalars
      use eos_vars
      use mixing
      implicit none

      ! DevinD - again lots of extra first touched variables that weren't in init_arrays originally?

#ifdef SOLVE3D
      allocate( u(GLOBAL_2D_ARRAY,N,3) )
      u(:,:,:,1)=init
      u(:,:,:,2)=0.0                             ! index 2 used on rhs u(indx) in pre_step for DC.
      u(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.
      allocate( v(GLOBAL_2D_ARRAY,N,3) )
      v(:,:,:,1)=init
      v(:,:,:,2)=0.0                             ! index 2 used on rhs v(indx) in pre_step for DC.
      v(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.

# if defined NHMG
      allocate( w(GLOBAL_2D_ARRAY,0:N,3) )
      allocate( nhdu(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdv(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdw(GLOBAL_2D_ARRAY,0:N,2) )
# endif

      allocate( FlxU(GLOBAL_2D_ARRAY,N) ) ; FlxU=init
      allocate( FlxV(GLOBAL_2D_ARRAY,N) ) ; FlxV=init
      allocate( We(GLOBAL_2D_ARRAY,0:N) ) ; We=init        ! explicit
      allocate( Wi(GLOBAL_2D_ARRAY,0:N) ) ; Wi=init        ! implicit

      allocate( Hz(GLOBAL_2D_ARRAY,N) )    ; Hz=init       ! height of rho-cell
      allocate( z_r(GLOBAL_2D_ARRAY,N) )   ; z_r=init      ! depth at rho-points
      allocate( z_w(GLOBAL_2D_ARRAY,0:N) ) ; z_w=init      ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      allocate( dzdxi(GLOBAL_2D_ARRAY,1:N)  )
      allocate( dzdeta(GLOBAL_2D_ARRAY,1:N) )
# endif
#endif  /* SOLVE3D */

      if (wrt_file_avg) then                                         ! wrt_file_avg prevents allocation when
                                                                     ! when variable flags are true but avg not wanted.
        if (wrt_avg_U) then
          allocate( u_avg(GLOBAL_2D_ARRAY,N) )
          u_avg=0.
        endif
        if (wrt_avg_V) then
          allocate( v_avg(GLOBAL_2D_ARRAY,N) )
          v_avg=0.
        endif
        if (wrt_avg_O) then
          allocate( w_avg(GLOBAL_2D_ARRAY,0:N) )
          w_avg=0
        endif
        if (wrt_avg_W) then
          allocate( wvl_avg(GLOBAL_2D_ARRAY,0:N) )
          wvl_avg=0.
        endif
        ! the following variables do not 'live' in this module, but in order to prevent
        ! circular reference cause by e.g. wrt_R logical, they are allocated here.
        if (wrt_avg_R) then
          allocate( rho_avg(GLOBAL_2D_ARRAY,N) )
          rho_avg=0.
        endif
        if (wrt_avg_Akv)  allocate( akv_avg(GLOBAL_2D_ARRAY,0:N) )
        if (wrt_avg_Akt)  allocate( akt_avg(GLOBAL_2D_ARRAY,0:N) )
# ifdef SALINITY
        if (wrt_avg_Aks)  allocate( aks_avg(GLOBAL_2D_ARRAY,0:N) )
# endif
# ifdef LMD_KPP
        if (wrt_avg_Hbls) allocate( hbl_avg(GLOBAL_2D_ARRAY) )
# endif
# ifdef LMD_BKPP
        if (wrt_avg_Hbbl) allocate( hbbl_avg(GLOBAL_2D_ARRAY) )
# endif

      endif  ! <-- wrt_file_avg

      end subroutine init_arrays_ocean3d  !]

! ----------------------------------------------------------------------
      subroutine wrt_his_ocean_vars  ![
      ! write ocean_vars variables to output netcdf file
      ! ocean_vars variables are calculated for t=n in timestep t=n
      ! (unlike u/v/etc which are calculated for t=n+1 in timestep t=n)

      use ncvars  ! REMOVE THIS LATER AND MOVED WHAT IS NEEDED HERE!
      use mixing
      use eos_vars
      use work_mod
      use grid
      use tracers    !, only: wrt_his_trc

      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_his                              ! current file output record
      integer,save           :: total_rec_his=0                      ! total his output records so far
      real,save              :: output_time_his=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_his
      integer                :: tile, ierr, i, j, k

      if (wrt_file_his) then

      if (.not. first_step) output_time_his = output_time_his + dt   ! only start count after first timestep
      first_step=.false.                                             ! as first step the ocean_vars values are for t=0

      if (output_time_his>=output_period_his .or.
     &    output_time_his==0                      ) then             ! catch to get t=0 output

        if (mod(total_rec_his,nrpf_his)==0) then
          call create_file_ocean_vars(total_rec_his,fname_his,.false.)
          rec_his = 0
        endif
        total_rec_his = total_rec_his +1
        rec_his = rec_his + 1

        ierr=nf90_open(fname_his,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/rec_his/))

        call write_time_step( ncid, rec_his, total_rec_his )

        start=1; start(3)=rec_his                                    ! back to 2D vars
        if (wrt_Z)  call ncwrite(ncid,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start)
        if (wrt_Ub) call ncwrite(ncid,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start)
        if (wrt_Vb) call ncwrite(ncid,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start)

#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_his called at the middle of next timestep as some variables only calculated there for t=n.

        start(3)=1; start(4)=rec_his
        if (wrt_U) call ncwrite(ncid,vname(1,indxU),u( 1:i1,j0:j1,:,nstp),start)
        if (wrt_V) call ncwrite(ncid,vname(1,indxV),v(i0:i1, 1:j1,:,nstp),start)

        call wrt_his_trc(ncid,start)                              ! tracer variables

        if (wrt_R) then
# ifdef SPLIT_EOS
          call ncwrite(ncid,vname(1,indxR),rho1(i0:i1,j0:j1,:),start)
# else
          call ncwrite(ncid,vname(1,indxR), rho(i0:i1,j0:j1,:),start)
# endif
        endif
        if (wrt_O) then                              ! s-coordinate omega vertical velocity (m/s).
          do k=0,N
            do j=0,Mm+1
              do i=0,Lm+1
                work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
              enddo
            enddo
          enddo
          call ncwrite(ncid,vname(1,indxO),work(i0:i1,j0:j1,:),start)
        endif
        if (wrt_W) then                              ! true vertical velocity (m/s).
          do tile=0,NSUB_X*NSUB_E-1
            call wvlcty (tile, work)
          enddo
          ! need 1:N as work is 1:N+1 in size
          ! put as 0:N-1 because z-dimension array is declared in correctly in wvlcty_tile.
          call ncwrite(ncid,vname(1,indxW),work(i0:i1,j0:j1,0:N-1),start)
        endif
        if (wrt_Akv) call ncwrite(ncid,vname(1,indxAkv),Akv(i0:i1,j0:j1,:),start)
        if (wrt_Akt) call ncwrite(ncid,vname(1,indxAkt),Akt(i0:i1,j0:j1,:,itemp),start)
# ifdef SALINITY
        if (wrt_Aks) call ncwrite(ncid,vname(1,indxAks),Akt(i0:i1,j0:j1,:,isalt),start)
# endif
        start(3)=rec_his                                      ! back to 2D vars
# ifdef LMD_KPP
        if (wrt_Hbls) call ncwrite(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1,nstp),start)
# endif
# ifdef LMD_BKPP
        if (wrt_Hbbl) call ncwrite(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1,nstp),start)
# endif
#endif

        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'ocean_vars :: wrote history, tdays =', tdays,
     &     'step =', iic-1, 'rec =', rec_his, '/', total_rec_his     ! -1 since iic=t(n+1)
        endif

        output_time_his=0
      endif

      endif  ! <-- wrt_file_his

      end subroutine wrt_his_ocean_vars  !]

! ----------------------------------------------------------------------
      subroutine wrt_avg_ocean_vars  ![
      ! write averaged variables to output netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      use ncvars  ! remove later
      use eos_vars
      use mixing
      use tracers
      use grid
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_avg                              ! current file output record
      integer,save           :: total_rec_avg=0                      ! total avg output records so far
      real,save              :: output_time_avg=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_avg
      integer :: tile, tn, ierr, k

      if (wrt_file_avg) then

      if (.not. first_step) then
!        call calc_avg_ocean_vars ! don't include t=0 in averaging. can't calc_avg here as rho, omega and w change
        output_time_avg = output_time_avg + dt                       ! only start count after first timestep
      endif
      first_step=.false.                                             ! as first step the ocean_vars values are for t=0

      if (output_time_avg>=output_period_avg) then                   ! catch to get t=0 output

        if (mod(total_rec_avg,nrpf_avg)==0) then
          call create_file_ocean_vars(total_rec_avg,fname_avg,.true.)
          rec_avg = 0
        endif
        total_rec_avg = total_rec_avg +1
        rec_avg = rec_avg + 1

        ierr=nf90_open(fname_avg,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/t_avg_ovars/),(/rec_avg/))

        start=1; start(3)=rec_avg                                    ! back to 2D vars
        if (wrt_avg_Z)   call ncwrite(ncid, vname(1,indxZ),  zeta_avg(i0:i1,j0:j1), start)
        if (wrt_avg_Ub)  call ncwrite(ncid, vname(1,indxUb), ubar_avg( 1:i1,j0:j1), start)
        if (wrt_avg_Vb)  call ncwrite(ncid, vname(1,indxVb), vbar_avg(i0:i1, 1:j1), start)
# ifdef SOLVE3D
        start(3)=1; start(4)=rec_avg
        if (wrt_avg_U)   call ncwrite(ncid, vname(1,indxU),     u_avg( 1:i1,j0:j1,:), start)
        if (wrt_avg_V)   call ncwrite(ncid, vname(1,indxV),     v_avg(i0:i1, 1:j1,:), start)

        call wrt_avg_trc(ncid,start)     ! STILL NEED TO UPDATE                                   ! tracer variables

        if (wrt_avg_R)   call ncwrite(ncid, vname(1,indxR),   rho_avg(i0:i1,j0:j1,:), start)
        if (wrt_avg_O) then
          do k=0,N
            w_avg(i0:i1,j0:j1,k) = w_avg(i0:i1,j0:j1,k)*pm(i0:i1,j0:j1)*pn(i0:i1,j0:j1)  ! convert before write.
          enddo
          call ncwrite(ncid, vname(1,indxO),     w_avg(i0:i1,j0:j1,:), start)            ! here rather than calc_avg
        endif                                                                            ! for efficiency
        ! need 1:N as work is 1:N+1 in size
        ! put as 0:N-1 because z-dimension array is declared in correctly in wvlcty_tile.
        if (wrt_avg_W)   call ncwrite(ncid, vname(1,indxW),   wvl_avg(i0:i1,j0:j1,0:N-1), start)
        if (wrt_avg_Akv) call ncwrite(ncid, vname(1,indxAkv), akv_avg(i0:i1,j0:j1,:), start)
        if (wrt_avg_Akt) call ncwrite(ncid, vname(1,indxAkt), akt_avg(i0:i1,j0:j1,:), start)
#  ifdef SALINITY
        if (wrt_avg_Aks) call ncwrite(ncid, vname(1,indxAks), aks_avg(i0:i1,j0:j1,:), start)
#  endif
        start(3)=rec_avg                                                ! back to 2D vars
#  ifdef LMD_KPP
        if (wrt_avg_Hbls) call ncwrite(ncid, vname(1,indxHbls),  hbl_avg(i0:i1,j0:j1), start)
#  endif
#  ifdef LMD_BKPP
        if (wrt_avg_Hbbl) call ncwrite(ncid, vname(1,indxHbbl), hbbl_avg(i0:i1,j0:j1), start)
#  endif
# endif /* SOLVE3D */

        navg_ovars=0
        output_time_avg=0

        ierr=nf90_close(ncid)
        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'ocean_vars :: wrote averages, tdays =', tdays,
     &        'step =', iic-1, 'rec =', rec_avg, '/', total_rec_avg  ! -1 since iic=t(n+1)
     &         MYID
        endif
      endif

      endif  ! <-- wrt_file_avg

      end subroutine wrt_avg_ocean_vars  !]

!----------------------------------------------------------------------
      subroutine create_file_ocean_vars(total_rec,fname,avg)  ![
      implicit none

      !input/output
      integer          ,intent(in)  :: total_rec           ! total netcdf records so far for these variables
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ierr,varid
      character(len=10),dimension(4) :: dimnames           ! dimension names
      integer,          dimension(4) :: dimsizes

      integer :: ndims ! debug

      if (      avg) fname=trim(output_root_name) / / '_avg.nc'
      if (.not. avg) fname=trim(output_root_name) / / '_his.nc'
      call create_filename(fname,total_rec)

      ierr=nf90_create(fname,nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) call handle_ierr(ierr,'create_file_ocean_vars: ')

      ! ADD 4 DEFAULT DIMENSIONS SO NCJOIN WILL WORK!!!

      varid = nccreate(ncid,'ocean_time',(/'time'/),(/0/),nf90_double)
      if (avg) then
        ierr=nf90_put_att(ncid,varid,'long_name',
     &                     'Average time since origin')
      else
        ierr=nf90_put_att(ncid,varid,'long_name',
     &                     'Time since origin')
      endif
      ierr=nf90_put_att(ncid,varid,'units','second' )

      call put_global_atts(ncid, ierr)                     ! put global attributes in file

      if (avg) then
        call def_vars_ocean_vars(  .true.  )
        ierr=nf90_put_att(ncid,nf90_global,'type','ROMS averages file')
        ierr=nf90_put_att(ncid,nf90_global,
     &         'output_period_avg_seconds',output_period_avg)
      else
        call def_vars_ocean_vars(  .false. )
        ierr=nf90_put_att(ncid,nf90_global,'type','ROMS history file')
        ierr=nf90_put_att(ncid,nf90_global,
     &         'output_period_his_seconds',output_period_his)
      endif

      ierr=nf90_enddef(ncid)
      if (mynode == 0) then
        write(*,'(7x,2A)')
     &    'ocean_vars :: created new netCDF file ', trim(fname)
      endif

      end subroutine create_file_ocean_vars !]

! ----------------------------------------------------------------------
      subroutine def_vars_ocean_vars( avg )  ![
      ! define output variable & attributes in netcdf results file
      ! for instantaneous or averaged variables
      use ncvars ! vname and indxZ
      use tracers
      implicit none

      ! input
      logical,intent(in) :: avg
      ! local
      integer           :: ierr, varid
      character(len=64) :: text_lname

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
      varid = nccreate(ncid,'time_step',(/'auxil','time'/),(/iaux,0/),nf90_int)
      ierr=nf90_put_att (ncid, varid, 'long_name',
     &       'time step and record numbers from initialization')

      if ( (wrt_Z .and. .not. avg) .or. (wrt_avg_Z .and. avg) ) then ! .or. needed for his or avg output selection
        if (.not. avg) text_lname=vname(2,indxZ)
        if (      avg) text_lname='averaged '/ /vname(2,indxZ)
        varid = nccreate(ncid,vname(1,indxZ),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxZ))
      endif
      if ( (wrt_Ub .and. .not. avg) .or. (wrt_avg_Ub .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxUb)
        if (      avg) text_lname='averaged '/ /vname(2,indxUb)
        varid = nccreate(ncid,vname(1,indxUb),(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxUb))
      endif
      if ( (wrt_Vb .and. .not. avg) .or. (wrt_avg_Vb .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxVb)
        if (      avg) text_lname='averaged '/ /vname(2,indxVb)
        varid = nccreate(ncid,vname(1,indxVb),(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxVb))
      endif
#ifdef SOLVE3D
      if ( (wrt_U .and. .not. avg) .or. (wrt_avg_U .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxU)
        if (      avg) text_lname='averaged '/ /vname(2,indxU)
        varid = nccreate(ncid,vname(1,indxU),(/'xi_u','eta_rho','s_rho','time'/),(/xi_u,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxU))
      endif
      if ( (wrt_V .and. .not. avg) .or. (wrt_avg_V .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxV)
        if (      avg) text_lname='averaged '/ /vname(2,indxV)
        varid = nccreate(ncid,vname(1,indxV),(/'xi_rho','eta_v','s_rho','time'/),(/xi_rho,eta_v,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxV))
      endif

      if (avg) then
        call def_avg_trc( ncid )
      else
        call def_his_trc( ncid )
      endif

      if ( (wrt_R .and. .not. avg) .or. (wrt_avg_R .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxR)
        if (      avg) text_lname='averaged '/ /vname(2,indxR)
        varid = nccreate(ncid,vname(1,indxR),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxR))
      endif
      if ( (wrt_O .and. .not. avg) .or. (wrt_avg_O .and. avg) ) then                                      ! s-coordinate "omega" vertical velocity.
        if (.not. avg) text_lname=vname(2,indxO)
        if (      avg) text_lname='averaged '/ /vname(2,indxO)
        varid = nccreate(ncid,vname(1,indxO),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxO))
      endif
      if ( (wrt_W .and. .not. avg) .or. (wrt_avg_W .and. avg) ) then                                      ! true W-vertical velocity.
        if (.not. avg) text_lname=vname(2,indxW)
        if (      avg) text_lname='averaged '/ /vname(2,indxW)
        varid = nccreate(ncid,vname(1,indxW),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxW))
      endif
      if ( (wrt_Akv .and. .not. avg) .or. (wrt_avg_Akv .and. avg) ) then                                    ! vertical viscosity coefficient.
        if (.not. avg) text_lname=vname(2,indxAkv)
        if (      avg) text_lname='averaged '/ /vname(2,indxAkv)
        varid = nccreate(ncid,vname(1,indxAkv),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkv))
      endif
      if ( (wrt_Akt .and. .not. avg) .or. (wrt_avg_Akt .and. avg) ) then                                    ! vertical thermal conductivity coefficient.
        if (.not. avg) text_lname=vname(2,indxAkt)
        if (      avg) text_lname='averaged '/ /vname(2,indxAkt)
        varid = nccreate(ncid,vname(1,indxAkt),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkt))
      endif
# ifdef SALINITY
      if ( (wrt_Aks .and. .not. avg) .or. (wrt_avg_Aks .and. avg) ) then                                    ! vertical diffusion coefficient for salinity.
        if (.not. avg) text_lname=vname(2,indxAks)
        if (      avg) text_lname='averaged '/ /vname(2,indxAks)
        varid = nccreate(ncid,vname(1,indxAks),(/'xi_rho','eta_rho','s_w','time'/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAks))
      endif
# endif /* SALINITY */
# ifdef LMD_KPP
      if ( (wrt_Hbls .and. .not. avg) .or. (wrt_avg_Hbls .and. avg) ) then                                   ! depth of surface planetary boundary layer.
        if (.not. avg) text_lname=vname(2,indxHbls)
        if (      avg) text_lname='averaged '/ /vname(2,indxHbls)
        varid = nccreate(ncid,vname(1,indxHbls),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbls))
      endif
# endif
# ifdef LMD_BKPP
      if ( (wrt_Hbbl .and. .not. avg) .or. (wrt_avg_Hbbl .and. avg) ) then                                   ! thickness of bottom boundary layer.
        if (.not. avg) text_lname=vname(2,indxHbbl)
        if (      avg) text_lname='averaged '/ /vname(2,indxHbbl)
        varid = nccreate(ncid,vname(1,indxHbbl),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbbl))
      endif
# endif
#endif /* SOLVE3D */

      end subroutine def_vars_ocean_vars  !]

! ----------------------------------------------------------------------
      subroutine calc_avg_ocean_vars  ![
      ! set averaged variables
      use eos_vars
      use work_mod
      use mixing
      use tracers
      use private_scratch
      use dimensions
      implicit none

      ! local
      real         :: coef
      logical,save :: first_step=.true.

      if (.not. first_step) then                           ! don't include t=0 in average

        navg_ovars = navg_ovars +1
        coef = 1./navg_ovars

        if (coef==1) then                                  ! this refreshes average (1-coef)=0
          if (mynode==0) write(*,'(7x,2A,F9.1)')
     &      'ocean_vars :: started averaging. ',
     &      'output_period_avg (s) =', output_period_avg
        endif

        t_avg_ovars = t_avg_ovars*(1-coef) + time*coef

        ! need i0:i1 indices because arrays still GLOBAL_2D therefore wasted margin
        if (wrt_avg_Z)  zeta_avg(i0:i1,j0:j1) = zeta_avg(i0:i1,j0:j1)*(1-coef) + zeta(i0:i1,j0:j1,knew)*coef
        if (wrt_avg_Ub) ubar_avg( 1:i1,j0:j1) = ubar_avg( 1:i1,j0:j1)*(1-coef) + ubar( 1:i1,j0:j1,knew)*coef
        if (wrt_avg_Vb) vbar_avg(i0:i1, 1:j1) = vbar_avg(i0:i1, 1:j1)*(1-coef) + vbar(i0:i1, 1:j1,knew)*coef
# ifdef SOLVE3D
        if (wrt_avg_U)   u_avg( 1:i1,j0:j1,:) = u_avg( 1:i1,j0:j1,:)*(1-coef) + u( 1:i1,j0:j1,:,nstp)*coef
        if (wrt_avg_V)   v_avg(i0:i1, 1:j1,:) = v_avg(i0:i1, 1:j1,:)*(1-coef) + v(i0:i1, 1:j1,:,nstp)*coef

        call set_avg_trc(coef)                                        ! tracer variables

        if (wrt_avg_R) then
#  ifdef SPLIT_EOS
          rho_avg(i0:i1,j0:j1,:) = rho_avg(i0:i1,j0:j1,:)*(1-coef) +
     &     ( rho1(i0:i1,j0:j1,:) - qp1(i0:i1,j0:j1,:) * z_r(i0:i1,j0:j1,:) )*coef
#  else
          rho_avg(i0:i1,j0:j1,:) = rho_avg(i0:i1,j0:j1,:)*(1-coef) + rho(i0:i1,j0:j1,:)*coef
#  endif
        endif
        if (wrt_avg_O) w_avg(i0:i1,j0:j1,:) = w_avg(i0:i1,j0:j1,:)*(1-coef)
     &                                      + ( We(i0:i1,j0:j1,:)+Wi(i0:i1,j0:j1,:) ) *coef
        ! need 1:N as work is 1:N+1 in size
        ! put as 0:N-1 because z-dimension array is declared in correctly in wvlcty_tile.
        if (wrt_avg_W) then
          ! the w averaging, calculation and output should be checked at some point.
          ! things seem inconsistent and buggy.
          call wvlcty_tile(1,nx,1,ny, work, A2d(1,1), A2d(1,2), A2d(1,3))  ! as per set_avg

          wvl_avg(i0:i1,j0:j1,1:N) = wvl_avg(i0:i1,j0:j1,1:N)*(1-coef)
     &                               + work(i0:i1,j0:j1,1:N)*coef
        endif
        if (wrt_avg_Akv) akv_avg(i0:i1,j0:j1,:) = akv_avg(i0:i1,j0:j1,:)*(1-coef) + akv(i0:i1,j0:j1,:)*coef
        if (wrt_avg_Akt) akt_avg(i0:i1,j0:j1,:) = akt_avg(i0:i1,j0:j1,:)*(1-coef) + akt(i0:i1,j0:j1,:,itemp)*coef
#  ifdef SALINITY
        if (wrt_avg_Aks) aks_avg(i0:i1,j0:j1,:) = aks_avg(i0:i1,j0:j1,:)*(1-coef) + akt(i0:i1,j0:j1,:,isalt)*coef
#  endif
#  ifdef LMD_KPP
        if (wrt_avg_Hbls) hbl_avg(i0:i1,j0:j1)  =  hbl_avg(i0:i1,j0:j1)*(1-coef) + hbls(i0:i1,j0:j1,nstp)*coef
#  endif
#  ifdef LMD_BKPP
        if (wrt_avg_Hbbl) hbbl_avg(i0:i1,j0:j1) = hbbl_avg(i0:i1,j0:j1)*(1-coef) + hbbl(i0:i1,j0:j1,nstp)*coef
#  endif
# endif /* SOLVE3D */

      endif  ! <-- .not. first_step

      first_step=.false.

      end subroutine calc_avg_ocean_vars  !]

! ----------------------------------------------------------------------
      subroutine write_time_step( ncid, record, total_recs )  ![
      use ncvars
      implicit none

      ! inputs
      integer, intent(in) :: ncid
      integer, intent(in) :: record                        ! current file record number
      integer, intent(in) :: total_recs                    ! total records for variable

      ! local
      integer :: ibuff(iaux), start(2), count(2)           ! iaux = 6 from wrt_his.F
      integer :: var_id_tmp, ierr

      ibuff(1)=iic-1  ; ibuff(2)=nrecrst                   ! time step and
      ibuff(4:iaux)=0 ; ibuff(3)=total_recs                ! rechis numbers.

      start(1)=1      ; count(1)=6                         ! iaux = 6 in ncvars
      start(2)=record ; count(2)=1

      ierr=nf90_inq_varid(ncid, 'time_step', var_id_tmp)
      if (ierr/=0) call handle_ierr(ierr,'ocean_vars write_time_step inq_varid')
      ierr=nf90_put_var(ncid, var_id_tmp, ibuff, start, count)  ! & record time step info
      if (ierr/=0) call handle_ierr(ierr,'ocean_vars write_time_step put_var')
      ! can't use ncwrite() as ibuff is an integer array not real

      end subroutine write_time_step  !]

!-----------------------------------------------------------------------
      subroutine display_output_settings_to_terminal  ![
      use tracers
      use ncvars
      implicit none

      ! local
      integer :: itrc

      if (mynode==0) write(*,'(/7x,2A,F6.1,2x,A,I4)')                ! write to terminal output in simulation pre-amble text which
     &  'ocean_vars :: history file ',                               ! result variables are being stored
     &  'ouput_period =', output_period_his,
     &  'recs/file =', nrpf_his

      if (mynode==0) write(*,'(9x,A,13(/11x,A,T20,L1,4x,A),/)')      ! 13(....) repeats formatting 11 times.
     &               'his fields to be saved: (T/F)'                 ! T20 moves to the 20th character on line
     &                 , 'zeta',   wrt_Z,    vname(2,indxZ)
     &                 , 'ubar',   wrt_Ub,   vname(2,indxUb)
     &                 , 'vbar',   wrt_Vb,   vname(2,indxVb)
#ifdef SOLVE3D
     &                 , 'u',      wrt_U,    vname(2,indxU)
     &                 , 'v',      wrt_V,    vname(2,indxV)

     &                 , 'rho',    wrt_R,    vname(2,indxR)
     &                 , 'Omega',  wrt_O,    vname(2,indxO)
     &                 , 'W',      wrt_W,    vname(2,indxW)

     &                 , 'Akv',    wrt_Akv,  vname(2,indxAkv)
     &                 , 'Akt',    wrt_Akt,  vname(2,indxAkt)
# ifdef SALINITY
     &                 , 'Aks',    wrt_Aks,  vname(2,indxAks)
# endif
# ifdef LMD_KPP
     &                 , 'hbls',   wrt_Hbls, vname(2,indxHbls)
# endif
# ifdef LMD_BKPP
     &                 , 'hbbl',   wrt_Hbbl, vname(2,indxHbbl)
# endif
          do itrc=1,iTandS+nt_passive
            if (mynode==0) write(*,'(11x,A,I2,A,T16,L1,T20,A)') 't(',
     &        itrc, ')', wrt_t(itrc), t_vname(itrc)
          enddo
#endif /* SOLVE3D */

      end subroutine display_output_settings_to_terminal  !]

!-----------------------------------------------------------------------
      subroutine display_output_settings_to_terminal_avg  ![
      use tracers
      use ncvars
      implicit none

      ! local
      integer :: itrc

      if (mynode==0) write(*,'(/7x,2A,F6.1,2x,A,I4)')                ! write to terminal output in simulation pre-amble text which
     &  'ocean_vars :: average file ',                               ! result variables are being stored
     &  'ouput_period =', output_period_avg,
     &  'recs/file =', nrpf_avg

      if (mynode==0) write(*,'(9x,A,13(/11x,A,T20,L1,4x,A),/)')      ! 13(....) repeats formatting 11 times.
     &               'avg fields to be saved: (T/F)'                 ! T20 moves to the 20th character on line
     &                 , 'zeta',   wrt_avg_Z,    vname(2,indxZ)
     &                 , 'ubar',   wrt_avg_Ub,   vname(2,indxUb)
     &                 , 'vbar',   wrt_avg_Vb,   vname(2,indxVb)
#ifdef SOLVE3D
     &                 , 'u',      wrt_avg_U,    vname(2,indxU)
     &                 , 'v',      wrt_avg_V,    vname(2,indxV)

     &                 , 'rho',    wrt_avg_R,    vname(2,indxR)
     &                 , 'Omega',  wrt_avg_O,    vname(2,indxO)
     &                 , 'W',      wrt_avg_W,    vname(2,indxW)

     &                 , 'Akv',    wrt_avg_Akv,  vname(2,indxAkv)
     &                 , 'Akt',    wrt_avg_Akt,  vname(2,indxAkt)
# ifdef SALINITY
     &                 , 'Aks',    wrt_avg_Aks,  vname(2,indxAks)
# endif
# ifdef LMD_KPP
     &                 , 'hbls',   wrt_avg_Hbls, vname(2,indxHbls)
# endif
# ifdef LMD_BKPP
     &                 , 'hbbl',   wrt_avg_Hbbl, vname(2,indxHbbl)
# endif
          do itrc=1,iTandS+nt_passive
            if (mynode==0) write(*,'(11x,A,I2,A,T16,L1,T20,A)') 't(',
     &        itrc, ')', wrt_t_avg(itrc), t_vname(itrc)
          enddo
#endif /* SOLVE3D */

      end subroutine display_output_settings_to_terminal_avg  !]

!----------------------------------------------------------------------
      subroutine def_rst(total_rec, ierr)  ![
      ! create or openrestart netCDF file.

      use tracers,         only: t_vname, t_lname, t_units, wrt_t
      use roms_read_write, only: output_root_name
      use ncvars
      use scalars
      use grid                                                       ! def_grid

      implicit none

      ! input
      integer :: total_rec, ierr
      ! local
      logical :: create_new_file
      integer :: lenstr, checkdims, my_nf_def_dim, varid
      integer :: rec, lrst, lvar, timedim, prev_fill_mode, itrc
#ifdef MASK_LAND_DATA
# ifdef DBLEPREC
      real(kind=8),parameter :: spv_set=1.D+33
# else
      real(kind=4),parameter :: spv_set=1.E+33
# endif
#endif

! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'rst_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.

      ierr=0
      rstname=trim(output_root_name) / / '_rst.nc'
      lrst=lenstr(rstname)
#if defined MPI && defined PARALLEL_FILES
        ! Insert MPI node numbers to file name
        call insert_node(rstname, lrst, mynode, NNODES, ierr)  ! previously done in ncdf_create_file:
#endif
      lrst=lenstr(rstname)
      if (nrpfrst > 0) then
        lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
        call insert_time_index (rstname, lrst, lvar, ierr)
        if (ierr /= 0) goto 99
      endif

! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncrst==-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just reopen, no checking, all ids are assumed
!        .... to be already known (this branch is only for
!        .... output into single file from MPI code).
!     endif
!
! which is designed to implement a flexible opening policy: if
!  (i) ldefhis=.true., it forces creation of a new file [if the file
!      already exists, it will be overwritten]; on the other hand, if
! (ii) ldefhis=.false., it is assumed that the file already exists
!      and an attempt to open it is made first; if the attempt is
!      successful, the file is prepared for appending hew data;
!      if it fails, a new file is created.

!      create_new_file=ldefhis
      create_new_file=.true.  ! DevinD edit to work
      if (ncrst /= -1) create_new_file=.false.
#if defined MPI && !defined PARALLEL_FILES
      if (mynode > 0) create_new_file=.false.
#endif

  10  if (create_new_file) then
        ierr=nf90_create(rstname(1:lrst), nf90_netcdf4, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf90_set_fill(ncrst, nf90_nofill, prev_fill_mode)
          if (ierr /= nf90_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_rst :: Cannot ',
     &               'create restart NetCDF file ''', rstname(1:lrst),
     &                                       '''',  nf90_strerror(ierr)
          goto 99              !--> ERROR
        endif

        call put_global_atts (ncrst, ierr)

#ifdef PUT_GRID_INTO_RESTART
        if (total_rec == 1) call def_grid(ncrst)
#endif

! Define evolving model variables:
! ------ -------- ----- ----------
! Time step number and time record numbers:
        varid = nccreate(ncrst,'time_step',(/'auxil','time'/),(/iaux,0/),nf90_int)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &             'time step and record numbers from initialization')

! Time.
        varid = nccreate(ncrst,vname(1,indxTime),(/'time'/),(/0/),nf90_double)
        ierr=nf90_put_att(ncrst, varid, 'long_name', vname(2,indxTime))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxTime))

! Free-surface.
        varid = nccreate(ncrst,vname(1,indxZ),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name', vname(2,indxZ))
        ierr=nf90_put_att (ncrst, varid, 'units', vname(3,indxZ))
        ierr = nf90_put_att(ncrst,varid,'_FillValue',spv_set)

! 2D momenta in XI- and ETA-components
        varid = nccreate(ncrst,vname(1,indxUb),(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                                   vname(2,indxUb))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxUb))
        ierr = nf90_put_att(ncrst,varid,'_FillValue',spv_set)
        varid = nccreate(ncrst,vname(1,indxVb),(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                                   vname(2,indxVb))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxVb))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        varid = nccreate(ncrst,'DU_avg2',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &           '<<fast-time averaged ubar(:,:,n+1/2)>>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DV_avg2',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &           '<<fast-time-averaged vbar(:,:,n+1/2)>>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DU_avg_bak',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged ubar(:,:,n-1)>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DV_avg_bak',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged vbar(:,:,n-1)>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
#  elif defined PRED_COUPLED_MODE
        varid = nccreate(ncrst,'rufrc_bak',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                        '3D-to-2D forcing, XI-component')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
        varid = nccreate(ncrst,'rvfrc_bak',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                       '3D-to-2D forcing, ETA-component')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
#  endif
# endif  /* EXACT_RESTART */

! 3D momenta in XI- and ETA-directions.
        varid = nccreate(ncrst, vname(1,indxU),(/'xi_u','eta_rho','s_rho','time'/),(/xi_u,eta_rho,N,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                     vname(2,indxU))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                     vname(3,indxU))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,vname(1,indxV),(/'xi_rho','eta_v','s_rho','time'/),(/xi_rho,eta_v,N,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                     vname(2,indxV))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                     vname(3,indxV))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

! Tracer variables.
        do itrc=1,NT
          varid = nccreate(ncrst,t_vname(itrc),(/'xi_rho','eta_rho','s_rho','time'/),
     &                                         (/xi_rho,eta_rho,N,0/), NF_FTYPE)
          ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                       t_lname(itrc))
          ierr=nf90_put_att (ncrst, varid, 'units', t_units(itrc))
          ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.
!        ierr=nf_def_var (ncrst, vname(1,indxHbls)(1:lvar),
!     &                            NF_FTYPE, 3, r2dgrd, rstHbls)
        varid = nccreate(ncrst,vname(1,indxHbls),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                 vname(2,indxHbls))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                 vname(3,indxHbls))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.
!        ierr=nf_def_var (ncrst, vname(1,indxHbbl),
!     &                            NF_FTYPE, 3, r2dgrd, rstHbbl)
        varid = nccreate(ncrst,vname(1,indxHbbl),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                 vname(2,indxHbbl))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                 vname(3,indxHbbl))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
# endif
#endif /* SOLVE3D */

        ierr=nf90_enddef(ncrst)
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
          write(*,'(7x,4A,1x,A,i4)') 'def_rst :: Created new ',
     &          'netCDF file ''', rstname(1:lrst), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Check consistency of the dimensions of fields from the file with
! model dimensions. Determine the current size of unlimited dimension
! and set initial record [in the case of MPI serialized output, at
! this moment the last time record is assumed to be **partially**
! written by MPI processes with lower rank. Thus the next write is
! expected to be into the same record rather than next one (except
! MPI-master, who initializes the record).
!
! In the case when file is rejected (whether it cannot be opened, or
! something is wrong with its dimensions, create new file.

      elseif (ncrst == -1) then
        ierr=nf90_open(rstname, nf90_write, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf90_set_fill(ncrst, nf90_nofill, prev_fill_mode)
          if (ierr /= nf90_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif

          ierr=checkdims(ncrst, rstname, rec)
          if (ierr == nf90_noerr) then
            if (nrpfrst==0) then
              ierr=rec+1 - nrecrst
            else
              ierr=rec+1 - (1+mod(nrecrst-1, abs(nrpfrst)))
            endif
            if (ierr > 0) then
              mpi_master_only write(*,
     &           '(/1x,A,I5,1x,3A/21x,A,I5,1x,A,I5,1x,A/21x,A/)')
     &           'WARNING: def_rst :: The actual number of records',
     &            rec, 'in netCDF file ''',  rstname(1:lrst), '''',
     &           'exceeds record', rec+1-ierr, '/', total_rec,
     &           'for restart initial conditions.',
     &           'All records beyond this will be overwriding.'
              rec=rec-ierr
            elseif (nrpfrst==0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI && !defined PARALLEL_FILES
              if (mynode > 0) total_rec=total_rec-1
#endif
            endif
            ierr=nf90_noerr
          endif
        endif

        if (ierr /= nf90_noerr) then
#if defined MPI && !defined PARALLEL_FILES
          if (mynode==0) then
            create_new_file=.true.
            goto 10
          else
            write(*,'(/1x,4A, 1x,A,I4/)')  '### ERROR: def_rst :: ',
     &        'Cannot open restart netCDF file ''', rstname(1:lrst),
     &        '''.' MYID
            goto 99                                     !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif

! Find netCDF variable IDs of evolving model variables:
! ---- ------ -------- --- -- -------- ----- ----------

! Time step indices and time:

        ierr=nf90_inq_varid(ncrst, 'time_step', rstTstep)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxTime))
          ierr=nf90_inq_varid(ncrst, vname(1,indxTime)(1:lvar), rstTime)
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxTime)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) 'time_step', rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

! Free-surface and 2D momenta in XI- and ETA-directions

        lvar=lenstr(vname(1,indxZ))
        ierr=nf90_inq_varid(ncrst, vname(1,indxZ)(1:lvar), rstZ)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf90_inq_varid(ncrst, vname(1,indxUb)(1:lvar), rstUb)
          if (ierr == nf90_noerr) then
            lvar=lenstr(vname(1,indxVb))
            ierr=nf90_inq_varid(ncrst, vname(1,indxVb)(1:lvar), rstVb)
            if (ierr /= nf90_noerr) then
              write(*,1) vname(1,indxVb)(1:lvar), rstname(1:lrst)
            endif
          else
            write(*,1) vname(1,indxUb)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxZ)(1:lvar),  rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        ierr=nf90_inq_varid(ncrst,   'DU_avg2',   rst_DU_avg2)
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncrst,  'DV_avg2',     rst_DV_avg2)
          if (ierr == nf90_noerr) then
            ierr=nf90_inq_varid(ncrst, 'DU_avg_bak', rst_DU_avg_bak)
            if (ierr == nf90_noerr) then
              ierr=nf90_inq_varid(ncrst, 'DV_avg_bak', rst_DV_avg_bak)
              if (ierr /= nf90_noerr) then
                write(*,1) 'DV_avg_bak', rstname(1:lrst)
              endif
            else
              write(*,1) 'DU_avg_bak', rstname(1:lrst)
            endif
          else
            write(*,1) 'DV_avg2', rstname(1:lrst)
          endif
        else
          write(*,1) 'DU_avg2', rstname(1:lrst)
        endif
#  elif defined PRED_COUPLED_MODE
        ierr=nf90_inq_varid(ncrst, 'rufrc_bak', rst_rufrc)
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncrst, 'rvfrc_bak', rst_rvfrc)
          if (ierr /= nf90_noerr) then
            write(*,1) 'rvfrc_bak', rstname(1:lrst)
          endif
        else
          write(*,1) 'rufrc_bak', rstname(1:lrst)
        endif
#  endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR
# endif

! 3D momenta n XI- and ETA-directions.

        lvar=lenstr(vname(1,indxU))
        ierr=nf90_inq_varid (ncrst, vname(1,indxU)(1:lvar), rstU)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf90_inq_varid (ncrst, vname(1,indxV)(1:lvar), rstV)
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxV)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxU)(1:lvar), rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

! Tracer variables.

        do itrc=1,NT
          lvar=lenstr(vname(1,indxT+itrc-1))
          ierr=nf90_inq_varid(ncrst, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  rstT(itrc))
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxT+itrc-1)(1:lvar), rstname(1:lrst)
            goto 99                                       !--> ERROR
          endif
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        lvar=lenstr(vname(1,indxHbls))
        ierr=nf90_inq_varid(ncrst, vname(1,indxHbls)(1:lvar), rstHbls)
        if (ierr /= nf90_noerr) then
          write(*,1) vname(1,indxHbls)(1:lvar), rstname(1:lrst)
          goto 99                                        !--> ERROR
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.
        ierr=nf90_inq_varid(ncrst, vname(1,indxHbbl), rstHbbl)
        if (ierr /= nf90_noerr) then
        write(*,1) vname(1,indxHbbl), rstname
          goto 99                                       !--> ERROR
        endif
# endif
#endif
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(*,'(7x,4A,I4,2x,A,I4)') 'def_rst :: open existing ',
     &            'file ''', rstname(1:lrst), ''' from rec =', rec
     &             MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

#if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf90_open(rstname(1:lrst), nf_write, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf_set_fill(ncrst, nf_nofill, prev_fill_mode)
          if (ierr == nf90_noerr) then
# ifndef MPI_SILENT_MODE
            write(*,'(7x,4A,I4)') 'def_rst :: reopen existing file ''',
     &              rstname(1:lrst), '''.' MYID
# endif
          else
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A,1x,A,I4/)') '### ERROR: def_rst :: Cannot ',
     &            'reopen file ''', rstname(1:lrst),'''.'
     &             MYID
          goto 99                                         !--> ERROR
        endif
#endif
      endif              !<-- create_new_file

#ifdef PUT_GRID_INTO_RESTART

! Write grid variables.
! ----- ---- ----------

      if (total_rec == 1) call wrt_grid(ncrst, rstname, lrst)
#endif
  99  return                                              !--> ERROR
   1  format(/1x,'### ERROR: def_rst :: Cannot find variable ''', A,
     &                            ''' in netCDF file ''', A, '''.'/)

      end subroutine def_rst  !]

!----------------------------------------------------------------------
      subroutine wrt_rst  ![
      ! write model prognostic variables into restart netCDF file.

      use tracers, only: wrt_rst_trc
      use coupling
      use mixing
      use ncvars
      use scalars
      use mpi
      use read_write                                       ! temporary addtion until writing fixed

      implicit none

      integer :: record, start(2), count(2), ibuff(iaux),
     &           i, ierr, lstr, lvar, lenstr, ncdf_write, varid, itrc
      character(len=18) :: tstring

#if defined MPI && !defined PARALLEL_FILES
      integer :: status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 1, ocean_grid_comm, status, ierr)
#endif

! Create/open restart and set record index within the file.

      call def_rst(nrecrst, ierr)
      if (ierr /= nf90_noerr) goto 99
      lstr=lenstr(rstname)
      nrecrst=max(nrecrst,1)                  !!! WARNING: Here it
      if (nrpfrst == 0) then                  !!! is assumed that
        record=nrecrst                        !!! global restart
      else                                    !!! record index
        record=1+mod(nrecrst-1, abs(nrpfrst)) !!! nrecrst is already
      endif                                   !!! advanced by main.

! Write out evolving model variables:
! ----- --- -------- ----- ----------

      ibuff(1)=iic-1  ; ibuff(2)=nrecrst        ! Time step and
      ibuff(4:iaux)=0 ; ibuff(3)=1 !nrechis     ! record numbers. DevinD put in 1 as don't need this anymore
!#ifdef AVERAGES
!      ibuff(4)=nrecavg  ! DevinD
!#endif
      start(1)=1      ; count(1)=iaux
      start(2)=record ; count(2)=1

      ierr=nf90_inq_varid(ncrst,'time_step',varid)
      ierr=nf90_put_var(ncrst, varid, ibuff, start, count)
      if (ierr /= nf90_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_rst :: ',
     &       'Cannot write variable ''time_step'' into restart file, ',
     &       'rec =', record, nf90_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif

! Time: Note that (1) time is always double precision; and (2) which
! one to write into the file, time or tdays, is decided by attribute
! "units" stored as vname(3,indxTime).  When computng time is always
! in seconds, however it is more convenient to have it in days in all
! the files; (3) this code is designed for politeness - it reads time
! units from the initial condition file (see "get_init.F") and uses
! it throughout the entire run.

      ierr=nf90_inq_varid(ncrst, 'ocean_time', varid)
      if (vname(3,indxTime)(1:6) == 'second') then
        ierr=nf90_put_var(ncrst, varid, time,(/record/))
      else
        ierr=nf90_put_var(ncrst, varid, tdays,(/record/))
      endif
      if (ierr /= nf90_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,'(/1x,4A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_rst :: ',
     &        'Cannot write variable ''',    vname(1,indxTime)(1:lvar),
     &        ''' into restart file, rec =', record, nf90_strerror(ierr)
     &            MYID
        goto 99                                           !--> ERROR
      endif

! Free-surface and 2D momentum components in XI- and ETA-directions.
      call nc_write_var(ncrst, zeta(START_2D_ARRAY,knew), 1,
     &                    vname(1,indxZ), rp_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, ubar(START_2D_ARRAY,knew), 1,
     &                    vname(1,indxUb), up_var, record, ierr)
        if (ierr == nf90_noerr) then
          call nc_write_var(ncrst, vbar(START_2D_ARRAY,knew), 1,
     &                      vname(1,indxVb), vp_var, record, ierr)
          if (ierr /= nf90_noerr) then
            lvar=lenstr(vname(1,indxVb))
            write(*,1) vname(1,indxVb)(1:lvar), record MYID
          endif
        else
          lvar=lenstr(vname(1,indxUb))
          write(*,1) vname(1,indxUb)(1:lvar), record MYID
        endif
      else
        lvar=lenstr(vname(1,indxZ))
        write(*,1) vname(1,indxZ)(1:lvar), record MYID
      endif
      if (ierr /= nf90_noerr) goto 99                !--> ERROR

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
      call nc_write_var(ncrst, DU_avg2, 1, 'DU_avg2', up_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, DV_avg2, 1, 'DV_avg2', vp_var, record, ierr)
        if (ierr == nf90_noerr) then
          call nc_write_var(ncrst, DU_avg_bak, 1, 'DU_avg_bak', up_var, record, ierr)
          if (ierr == nf90_noerr) then
            call nc_write_var(ncrst, DV_avg_bak, 1, 'DV_avg_bak', vp_var, record, ierr)
            if (ierr /= nf90_noerr) then
              write(*,1) 'DV_avg_bak', record MYID
            endif
          else
            write(*,1) 'DU_avg_bak', record MYID
          endif
        else
          write(*,1) 'DV_avg2', record MYID
        endif
      else
        write(*,1) 'DU_avg2', record MYID
      endif
#  elif defined PRED_COUPLED_MODE
      call nc_write_var(ncrst, rufrc_bak(START_2D_ARRAY,nstp), 1, 'rufrc', up_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, rvfrc_bak(START_2D_ARRAY,nstp), 1, 'rvfrc', vp_var, record, ierr)

        if (ierr /= nf90_noerr) then
          write(*,1) 'rvfrc_bak', record MYID
        endif
      else
        write(*,1) 'rufrc_bak', record MYID
      endif
#  endif
      if (ierr /= nf90_noerr) goto 99                !--> ERROR
# endif

! 3D momentum components in XI- and ETA-directions.
      call nc_write_var(ncrst, u(START_2D_ARRAY,1,nstp), N,
     &                    vname(1,indxU), up_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, v(START_2D_ARRAY,1,nstp), N,
     &                    vname(1,indxV), vp_var, record, ierr)
        if (ierr /= nf90_noerr) then
          lvar=lenstr(vname(1,indxV))
          write(*,1) vname(1,indxV)(1:lvar), record MYID
        endif
      else
        lvar=lenstr(vname(1,indxU))
        write(*,1) vname(1,indxU)(1:lvar), record MYID
      endif
      if (ierr /= nf90_noerr) goto 99                !--> ERROR

! Tracer variables.
      call wrt_rst_trc(record)

# ifdef LMD_KPP
! Depth of planetary boundary layer [m].
      call nc_write_var(ncrst, hbls(START_2D_ARRAY,nstp), 1,
     &                    vname(1,indxHbls), rp_var, record, ierr)
      if (ierr /= nf90_noerr) then
        lvar=lenstr(vname(1,indxHbls))
        write(*,1) vname(1,indxHbls)(1:lvar), record MYID
        goto 99                                         !--> ERROR
      endif
# endif
# ifdef LMD_BKPP
! Thickness of KPP bottom boundary layer [m].
      call nc_write_var(ncrst, hbbl(START_2D_ARRAY,nstp), 1,
     &                    vname(1,indxHbbl), rp_var, record, ierr)
      if (ierr /= nf90_noerr) then
        lvar=lenstr(vname(1,indxHbbl))
        write(*,1) vname(1,indxHbbl)(1:lvar), record MYID
        goto 99                                         !--> ERROR
      endif
# endif
#endif /* SOLVE3D */


  1   format(/1x, '### ERROR: wrt_rst :: Cannot write variable ''',
     &             A, ''' into restart file, rec =', i6, 3x,A,i4/)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to disk to the allow other
! processes to access data immediately after it is written.

#if defined MPI && !defined PARALLEL_FILES
      ierr=nf90_close (ncrst)
      if (nrpfrst > 0 .and. record >= nrpfrst) ncrst=-1
#else
      if (nrpfrst > 0 .and. record >= nrpfrst) then
        ierr=nf90_close (ncrst)
        ncrst=-1
      else
        ierr=nf90_sync(ncrst)
      endif
#endif
      if (ierr == nf90_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(tstring,'(F18.8)') tdays
          i=1
          do while (i<18 .and. tstring(i:i)==' ')
            i=i+1
          enddo
          write(*,'(7x,A,1x,A,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &      'wrt_rst :: wrote restart, tdays =', tstring(i:i+8),
     &      'step =', ibuff(1),  'rec =', record, '/',  nrecrst
     &       MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(/1x,2A/)')    '### ERROR: wrt_rst :: Cannot ',
     &                          'synchronize/close restart file.'
        if (may_day_flag == 0) may_day_flag=3
      endif

#if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                        mynode+1, 1, ocean_grid_comm, ierr)
#endif

      end subroutine wrt_rst  !]

!----------------------------------------------------------------------

      end module ocean_vars
