
      module basic_output
      ! Basic output routines for ocean variables
      ! his, avg, and rst files

#include "cppdefs.opt"

      use param
      use dimensions
      use netcdf
      use nc_read_write
      use roms_read_write
      use obc_tune
      use diagnostics
      use scalars
      use ocean_vars
      use eos_vars
      use grid
      use coupling
      use mixing
      use tracers

      implicit none
      private

#include "ocean_vars.opt"

      real,allocatable,dimension(:,:,:) :: akv_avg
      real,allocatable,dimension(:,:,:) :: akt_avg

      real,allocatable,dimension(:,:)   :: hbls_avg
      real,allocatable,dimension(:,:)   :: hbbl_avg
      real,allocatable,dimension(:,:,:) :: rho_avg

      integer, parameter, public :: iaux=6           ! Length of netCDF variable "time_step"


      ! netcdf outputting:
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_ovars=0
      integer :: navg_ovars = 0                             ! number of samples in average
      character(len=99),public :: fname_rst                 ! restart filename needed by diagnostics
      integer,public           :: rec_rst=nrpf_rst          ! current file output record needed by diags
      integer                :: prev_month
      logical                :: write_another_step = .false.

      public :: wrt_his_ocean_vars
      public :: wrt_avg_ocean_vars
      public :: calc_avg_ocean_vars
      public :: wrt_rst_ocean_vars
      public :: init_avg_arrays

      contains  !]
!----------------------------------------------------------------------
      subroutine init_restarts  ![

      implicit none

      integer,dimension(6)   :: date

      if (wrt_file_rst .and. monthly_restarts) then
          call sec2date(time,date)
          prev_month = date(2)
      endif

      end subroutine init_restarts !]
!----------------------------------------------------------------------
      subroutine init_avg_arrays  ![

      implicit none

      if (wrt_file_avg) then
        if (wrt_avg_Z)  then
          allocate( zeta_avg(GLOBAL_2D_ARRAY) )
          zeta_avg=0.                            ! avg needs to be 0. because Nan x 0. = Nan in set_avg.F
        endif
        if (wrt_avg_Ub) then
          allocate( ubar_avg(GLOBAL_2D_ARRAY) )
          ubar_avg=0.
        endif
        if (wrt_avg_Vb) then
          allocate( vbar_avg(GLOBAL_2D_ARRAY) )
          vbar_avg=0.
        endif

        if (wrt_avg_U) then
          allocate( u_avg(GLOBAL_2D_ARRAY,N) )
          u_avg=0.
        endif
        if (wrt_avg_V) then
          allocate( v_avg(GLOBAL_2D_ARRAY,N) )
          v_avg=0.
        endif

        if (wrt_avg_Tem) then
          allocate( tem_avg(GLOBAL_2D_ARRAY,N) )
          u_avg=0.
        endif
        if (wrt_avg_Sal) then
          allocate( sal_avg(GLOBAL_2D_ARRAY,N) )
          sal_avg=0.
        endif

        if (wrt_avg_O) then
          allocate( w_avg(GLOBAL_2D_ARRAY,0:N) )
          w_avg=0
        endif
        if (wrt_avg_W) then
          allocate( wvl_avg(GLOBAL_2D_ARRAY,nz) )
          wvl_avg=0.
        endif

        ! the following variables do not 'live' in this module, but in order to prevent
        ! circular reference cause by e.g. wrt_R logical, they are allocated here.

        if (wrt_avg_R) then
          allocate( rho_avg(GLOBAL_2D_ARRAY,nz) )
          rho_avg=0.
        endif
        if (wrt_avg_Akv) then
          allocate( akv_avg(GLOBAL_2D_ARRAY,0:nz) )
          Akv_avg = 0
        endif
        if (wrt_avg_Akt) then
          allocate( akt_avg(GLOBAL_2D_ARRAY,0:nz) )
          Akt_avg = 0
        endif

# ifdef LMD_KPP
        if (wrt_avg_hbls) allocate( hbls_avg(GLOBAL_2D_ARRAY) )
# endif
# ifdef LMD_BKPP
        if (wrt_avg_hbbl) allocate( hbbl_avg(GLOBAL_2D_ARRAY) )
# endif

      endif  ! <-- wrt_file_avg

      ! NOTE: call this here to prevent output in between timesteps which would break code
      ! checking script.
      if (wrt_file_his) call display_output_settings_to_terminal
      if (wrt_file_avg) call display_output_settings_to_terminal_avg

      end subroutine init_avg_arrays !]
! ----------------------------------------------------------------------
      subroutine wrt_his_ocean_vars(special)  ![
      ! write ocean_vars variables to output netcdf file
      ! ocean_vars variables are calculated for t=n in timestep t=n
      ! (unlike u/v/etc which are calculated for t=n+1 in timestep t=n)

      implicit none

      ! import/export
      logical,optional       :: special !! for blowup and initial history file
      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_his=nrpf_his                     ! current file output record
      integer,save           :: total_rec_his=0                      ! total his output records so far
      real,save              :: output_time_his=output_period_his+10
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_his
      integer                :: tile, ierr, i, j, k
      logical,save           :: supress=.false.
      logical                :: blowup,initial

      real,dimension(:,:,:),allocatable  :: work


      blowup = .false.
      initial= .false.
      if (present(special)) then
        if (special) then
          blowup = .true.
        else
          initial = .true.           
        endif
      endif


      if (wrt_file_his.or.blowup) then

        allocate( work(1-bf:nx+bf,1-bf:ny+bf,0:nz))

        if (.not.present(special)) then 
         output_time_his = output_time_his + dt  
        endif


      if (output_time_his>=output_period_his.or.present(special) ) then

         if (initial) then
           if (mynode==0) print *,'Writing initial history file'
         elseif (blowup) then
           if (mynode==0) print *,'Writing emergency history file'
         endif

        if (rec_his==nrpf_his) then
          call create_file_ocean_vars(fname_his,.false.)
          rec_his = 0
        endif
        rec_his = rec_his + 1

        ierr=nf90_open(fname_his,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/rec_his/))

        call write_time_step( ncid, rec_his, total_rec_his )

        start=1; start(3)=rec_his                                    ! back to 2D vars
        if (wrt_Z)  call ncwrite(ncid,'zeta',zeta(i0:i1,j0:j1,knew),start)
        if (wrt_Ub) call ncwrite(ncid,'ubar',ubar( 1:i1,j0:j1,knew),start)
        if (wrt_Vb) call ncwrite(ncid,'vbar',vbar(i0:i1, 1:j1,knew),start)

! 3D momentum components in x- and y-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_his called at the middle of next timestep as some variables only calculated there for t=n.

        start(3)=1; start(4)=rec_his
        if (wrt_U) call ncwrite(ncid,'u',u( 1:i1,j0:j1,:,nnew),start)
        if (wrt_V) call ncwrite(ncid,'v',v(i0:i1, 1:j1,:,nnew),start)

        if (wrt_Tem) call ncwrite(ncid,'temp',t(i0:i1,j0:j1,:,nnew,1),start)
        if (wrt_Sal) call ncwrite(ncid,'salt',t(i0:i1,j0:j1,:,nnew,2),start)

        call wrt_his_trc(ncid,start)                              ! tracer variables

        if (wrt_R) then
# ifdef SPLIT_EOS
          call ncwrite(ncid,'rho', rho1(i0:i1,j0:j1,:),start)
# else
          call ncwrite(ncid,'rho', rho(i0:i1,j0:j1,:),start)
# endif
        endif
        if (wrt_O) then                              ! s-coordinate omega vertical velocity (m/s).
          do k=0,N
            do j=0,Mm+1
              do i=0,Lm+1
                work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
                ! Implicit omega only
!               work(i,j,k)=pm(i,j)*pn(i,j)*Wi(i,j,k)
              enddo
            enddo
          enddo
          call ncwrite(ncid,'omega',work(i0:i1,j0:j1,:),start)
        endif
        if (wrt_W) then                              ! true vertical velocity (m/s).
          call wvlcty(work(:,:,1:nz))
          call ncwrite(ncid,'w',work(i0:i1,j0:j1,1:nz),start)
        endif
        if (wrt_Akv) call ncwrite(ncid,'Akv',Akv(i0:i1,j0:j1,:),start)
        if (wrt_Akt) call ncwrite(ncid,'Akt',Akt(i0:i1,j0:j1,:),start)

        start(3)=rec_his                                      ! back to 2D vars
        if (wrt_hbls) call ncwrite(ncid,'hbls',hbls(i0:i1,j0:j1),start)
        if (wrt_hbbl) call ncwrite(ncid,'hbbl',hbbl(i0:i1,j0:j1),start)

        ierr=nf90_close(ncid)
        if (mynode == 0 .and. (.not.code_check)) then                 ! supress output for code_check
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'ocean_vars :: wrote history, tdays =', tdays,
     &     'step =', iic, 'rec =', rec_his
        endif
        code_check=.false.

        output_time_his=0
      endif

      deallocate(work)

      endif  ! <-- wrt_file_his


      end subroutine wrt_his_ocean_vars  !]
! ----------------------------------------------------------------------
      subroutine wrt_avg_ocean_vars  ![
      ! write averaged variables to output netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.

      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_avg=nrpf_avg
      integer,save           :: total_rec_avg=0                      ! total avg output records so far
      real,save              :: output_time_avg=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_avg
      character(len=99)      :: output_time_string
      character(len=99)      :: formatted_string
      integer :: tile, tn, ierr, k

      if (wrt_file_avg) then

        output_time_avg = output_time_avg + dt

        if (output_time_avg>=output_period_avg) then

          ! save info to netcdf
          avg_output = ''
          write(output_time_string,'(F12.1)') output_time_avg
          write(formatted_string,'(A,A)') trim(adjustl(output_time_string)),
     &                                    ' seconds'
          call store_string_att(avg_output,formatted_string)

          if (rec_avg == nrpf_avg) then
            call create_file_ocean_vars(fname_avg,.true.)
            rec_avg = 0
          endif
        rec_avg = rec_avg + 1

        ierr=nf90_open(fname_avg,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/t_avg_ovars/),(/rec_avg/))

        start=1; start(3)=rec_avg                                    ! back to 2D vars
        if (wrt_avg_Z)   call ncwrite(ncid,'zeta', zeta_avg(i0:i1,j0:j1), start)
        if (wrt_avg_Ub)  call ncwrite(ncid,'ubar', ubar_avg( 1:i1,j0:j1), start)
        if (wrt_avg_Vb)  call ncwrite(ncid,'vbar', vbar_avg(i0:i1, 1:j1), start)

        start(3)=1; start(4)=rec_avg
        if (wrt_avg_U)   call ncwrite(ncid,'u',u_avg( 1:i1,j0:j1,:), start)
        if (wrt_avg_V)   call ncwrite(ncid,'v',v_avg(i0:i1, 1:j1,:), start)

        if (wrt_avg_Tem)   call ncwrite(ncid,'temp',tem_avg( 1:i1,j0:j1,:), start)
        if (wrt_avg_Sal)   call ncwrite(ncid,'salt',sal_avg(i0:i1, 1:j1,:), start)

        call wrt_avg_trc(ncid,start)     ! STILL NEED TO UPDATE                                   ! tracer variables

        if (wrt_avg_R)   call ncwrite(ncid,'rho',rho_avg(i0:i1,j0:j1,:), start)

        if (wrt_avg_O) then
          do k=0,nz
            w_avg(i0:i1,j0:j1,k) = w_avg(i0:i1,j0:j1,k)*pm(i0:i1,j0:j1)*pn(i0:i1,j0:j1)  ! convert before write.
          enddo
          call ncwrite(ncid,'omega',w_avg(i0:i1,j0:j1,:), start)            ! here rather than calc_avg
        endif                                                                            ! for efficiency
        if (wrt_avg_W)   call ncwrite(ncid,'w',   wvl_avg(i0:i1,j0:j1,1:N), start)
        if (wrt_avg_Akv) call ncwrite(ncid,'Akv', akv_avg(i0:i1,j0:j1,:), start)
        if (wrt_avg_Akt) call ncwrite(ncid,'Akt', akt_avg(i0:i1,j0:j1,:), start)

        start(3)=rec_avg                                                ! back to 2D vars
#  ifdef LMD_KPP
        if (wrt_avg_Hbls) call ncwrite(ncid,'hbls',  hbls_avg(i0:i1,j0:j1), start)
#  endif
#  ifdef LMD_BKPP
        if (wrt_avg_Hbbl) call ncwrite(ncid,'hbbl', hbbl_avg(i0:i1,j0:j1), start)
#  endif

        navg_ovars=0
        output_time_avg=0

        ierr=nf90_close(ncid)
        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'ocean_vars :: wrote averages, tdays =', tdays,
     &        'step =', iic, 'rec =', rec_avg, '/', total_rec_avg  ! 
     &         MYID
        endif
      endif

      endif  ! <-- wrt_file_avg

      end subroutine wrt_avg_ocean_vars  !]
!----------------------------------------------------------------------
      subroutine wrt_rst_ocean_vars  ![
      ! write ocean_vars variables to output netcdf restart file
      ! ocean_vars variables are calculated for t=n in timestep t=n
      ! (unlike u/v/etc which are calculated for t=n+1 in timestep t=n)

      implicit none

      ! local
      real,save              :: output_time_rst=0                    ! time since last output
      integer,dimension(6)   :: date
      character(len=15)      :: datestr

      if (wrt_file_rst) then

      if (monthly_restarts) then
          ! We need to save the last two timesteps of the month, so that when the new simulation
          ! loads the restart file the first timestep it takes will be the first timestep of the
          ! new month.
          ! time + 2*dt here so that it looks two timesteps ahead to see if we're changing months.
          call sec2date(time+2*dt,date)
          if ((date(2) - prev_month) /= 0) then
            call wrt_restart_file
            write_another_step = .true.
          endif

          if (write_another_step) then
            call wrt_restart_file
            write_another_step = .false.
          endif

          prev_month = date(2)
      else
!     if (.not. first_step) output_time_rst = output_time_rst + dt   ! only start count after first timestep
!     first_step=.false.                                             ! as first step the ocean_vars values are for t=0
          output_time_rst = output_time_rst + dt

          if (output_time_rst    >=output_period_rst
#ifdef EXACT_RESTART
     &        .or.
     &        output_time_rst+dt >=output_period_rst         ! timestep before
#endif
     &        ) then
             call wrt_restart_file
             if (output_time_rst>=output_period_rst) then
                 output_time_rst=0                             ! catch for exact restart to not
             endif                                             ! reset the time
          endif ! it output_time > output_period
      endif ! monthly_restarts

      endif ! wrt_file_rst

      end subroutine wrt_rst_ocean_vars  !]
!----------------------------------------------------------------------
      subroutine wrt_restart_file  ![
      ! write ocean_vars variables to output netcdf restart file
      ! ocean_vars variables are calculated for t=n in timestep t=n
      ! (unlike u/v/etc which are calculated for t=n+1 in timestep t=n)

#ifdef MARBL
!     add MARBL saved state to restart to restart file
      use marbl_driver, only: marbldrv_write_ss_vars_to_rst
#endif       
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: total_rec_rst=0                      ! total his output records so far
      real,save              :: output_time_rst=0                    ! time since last output
      logical,save           :: first_step=.true.
      integer                :: tile, ierr, i, j, k


        if (rec_rst==nrpf_rst) then
#ifdef EXACT_RESTART
          time = time + dt
#endif
          call create_file_rst_ocean_vars(fname_rst)
#ifdef EXACT_RESTART
          time = time - dt
#endif
          rec_rst = 0
        endif
        total_rec_rst = total_rec_rst +1
        rec_rst = rec_rst + 1

        ierr=nf90_open(fname_rst,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/rec_rst/))

        call write_time_step( ncid, rec_rst, total_rec_rst )

        start=1; start(3)=rec_rst                                    ! back to 2D vars
        call ncwrite(ncid,'zeta',zeta(i0:i1,j0:j1,knew),start)
        call ncwrite(ncid,'ubar',ubar( 1:i1,j0:j1,knew),start)
        call ncwrite(ncid,'vbar',vbar(i0:i1, 1:j1,knew),start)

! 3D momentum components in x- and y-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_rst called at the middle of next timestep as some variables only calculated there for t=n.

        start(3)=1; start(4)=rec_rst
        call ncwrite(ncid,'u',u( 1:i1,j0:j1,:,nnew),start)
        call ncwrite(ncid,'v',v(i0:i1, 1:j1,:,nnew),start)

        call ncwrite(ncid,'temp',t(i0:i1,j0:j1,:,nnew,1),start)
        call ncwrite(ncid,'salt',t(i0:i1,j0:j1,:,nnew,2),start)

        call wrt_rst_trc(ncid,start)  ! beyond T and S tracer variables

        start(3)=rec_rst                             ! back to 2D vars

# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        call ncwrite(ncid,'DU_avg1',    DU_avg1( 1:i1,j0:j1),start)
        call ncwrite(ncid,'DV_avg1',    DV_avg1(i0:i1, 1:j1),start)
        call ncwrite(ncid,'DU_avg2',    DU_avg2( 1:i1,j0:j1),start)
        call ncwrite(ncid,'DV_avg2',    DV_avg2(i0:i1, 1:j1),start)
        call ncwrite(ncid,'DU_avg_bak', DU_avg_bak( 1:i1,j0:j1),start)
        call ncwrite(ncid,'DV_avg_bak', DV_avg_bak(i0:i1, 1:j1),start)

#  elif defined PRED_COUPLED_MODE
        call ncwrite(ncid,'rufrc',rufrc_bak( 1:i1,j0:j1,nstp),start)
        call ncwrite(ncid,'rvfrc',rvfrc_bak(i0:i1, 1:j1,nstp),start)
#  endif
# endif


# ifdef LMD_KPP
        call ncwrite(ncid,'hbls',hbls(i0:i1,j0:j1),start)
# endif
# ifdef LMD_BKPP
        call ncwrite(ncid,'hbbl',hbbl(i0:i1,j0:j1),start)
# endif

        if (obc_tuning) call wrt_rst_tune(ncid,rec_rst)
        if (diag_pflx)  call wrt_rst_diag_slow(ncid,rec_rst)

#ifdef MARBL
        call marbldrv_write_ss_vars_to_rst(ncid,rec_rst)
#endif         

        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'ocean_vars :: wrote restart, tdays =', tdays,
     &     'step =', iic, 'rec =', rec_rst, '/', total_rec_rst
        endif

      end subroutine wrt_restart_file  !]
!----------------------------------------------------------------------
      subroutine create_file_ocean_vars(fname,avg)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ierr,varid

      if (avg) then
        call create_file('_avg',fname)
      else
        call create_file('_his',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      if (avg) then
        call def_vars_ocean_vars(  .true.  )
        ierr=nf90_put_att(ncid,nf90_global,'type','ROMS averages file')
        ierr=nf90_put_att(ncid,nf90_global,'averaging_timescale',adjustl(avg_output))
      else
        call def_vars_ocean_vars(  .false. )
        ierr=nf90_put_att(ncid,nf90_global,'type','ROMS history file')
      endif

      ierr=nf90_close(ncid)

      end subroutine create_file_ocean_vars !]
!----------------------------------------------------------------------
      subroutine create_file_rst_ocean_vars(fname)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      ! local
      integer :: ierr,varid

      call create_file('_rst',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      ierr=nf90_put_att(ncid,nf90_global,'type','ROMS restart file')

      call def_vars_rst_ocean_vars

      ierr=nf90_close(ncid)

      end subroutine create_file_rst_ocean_vars !]
! ----------------------------------------------------------------------
      subroutine def_vars_ocean_vars( avg )  ![
      ! define output variable & attributes in netcdf results file
      ! for instantaneous or averaged variables

      use tracers
      implicit none

      ! input
      logical,intent(in) :: avg
      ! local
      integer           :: ierr, varid
      character(len=64) :: text_lname
      character(len=7)  :: dn_aux = 'auxil'

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
      varid = nccreate(ncid,'time_step',(/dn_aux,dn_tm/),(/iaux,0/),nf90_int)
      ierr=nf90_put_att (ncid, varid, 'long_name',
     &       'time step and record numbers from initialization')

      if ( (wrt_Z .and. .not. avg) .or. (wrt_avg_Z .and. avg) ) then ! .or. needed for his or avg output selection
        varid = nccreate(ncid,'zeta',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','free surface elevation')
        ierr = nf90_put_att(ncid,varid,'units','meter')
      endif
      if ( (wrt_Ub .and. .not. avg) .or. (wrt_avg_Ub .and. avg) ) then
        varid = nccreate(ncid,'ubar',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,
     &                  'long_name','vertically averaged u-momentum')
        ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif
      if ( (wrt_Vb .and. .not. avg) .or. (wrt_avg_Vb .and. avg) ) then
        varid = nccreate(ncid,'vbar',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/))
        ierr = nf90_put_att(ncid,varid,
     &                  'long_name','vertically averaged v-momentum')
        ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif

      if ( (wrt_U .and. .not. avg) .or. (wrt_avg_U .and. avg) ) then
        varid = nccreate(ncid,'u',(/dn_xu,dn_yr,dn_zr,dn_tm/),(/xi_u,eta_rho,nz,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','u-momentum')
        ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif

      if ( (wrt_V .and. .not. avg) .or. (wrt_avg_V .and. avg) ) then
        varid = nccreate(ncid,'v',(/dn_xr,dn_yv,dn_zr,dn_tm/),(/xi_rho,eta_v,nz,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','v-momentum')
        ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif

      if ( (wrt_Tem .and. .not. avg) .or. (wrt_avg_tem .and. avg) ) then
        varid = nccreate(ncid,'temp',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,nz,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','potential temperature')
        ierr = nf90_put_att(ncid,varid,'units','Degrees Celsius')
      endif

      if ( (wrt_Sal .and. .not. avg) .or. (wrt_avg_sal .and. avg) ) then
        varid = nccreate(ncid,'salt',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,nz,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','salinity')
        ierr = nf90_put_att(ncid,varid,'units','PSU')
      endif

      if ( (wrt_R .and. .not. avg) .or. (wrt_avg_R .and. avg) ) then
        varid = nccreate(ncid,'rho',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,nz,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','density anomaly')
        ierr = nf90_put_att(ncid,varid,'units','kilogram meter-3')
      endif

      if ( (wrt_O .and. .not. avg) .or. (wrt_avg_O .and. avg) ) then
        varid = nccreate(ncid,'omega',(/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,nz+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','S-coordinate vertical flux')
        ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif
      if ( (wrt_W .and. .not. avg) .or. (wrt_avg_W .and. avg) ) then
        varid = nccreate(ncid,'w',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,nz,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','vertical velocity')
        ierr = nf90_put_att(ncid,varid,'units','meter second-1')
      endif
      if ( (wrt_Akv .and. .not. avg) .or. (wrt_avg_Akv .and. avg)) then
        varid = nccreate(ncid,'Akv',(/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,nz+1,0/))
        ierr = nf90_put_att(ncid,varid,
     &      'long_name','vertical momentum mixing coefficient')
        ierr = nf90_put_att(ncid,varid,'units','meter2 second-1')
      endif
      if ( (wrt_Akt .and. .not. avg) .or. (wrt_avg_Akt .and. avg) ) then
        varid = nccreate(ncid,'Akt',
     &      (/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,nz+1,0/))
        ierr = nf90_put_att(ncid,varid,
     &      'long_name','vertical tracer mixing coefficient')
        ierr = nf90_put_att(ncid,varid,'units','meter2 second-1')
      endif

# ifdef LMD_KPP
      if ( (wrt_Hbls .and. .not. avg) .or. (wrt_avg_Hbls .and. avg) ) then 
        varid = nccreate(ncid,'hbls',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','surface BL thickness')
        ierr = nf90_put_att(ncid,varid,'units','meter')
      endif
# endif
# ifdef LMD_BKPP
      if ( (wrt_Hbbl .and. .not. avg) .or. (wrt_avg_Hbbl .and. avg) ) then
        varid = nccreate(ncid,'hbbl',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name','bottom BL thickness')
        ierr = nf90_put_att(ncid,varid,'units','meter')
      endif
# endif

      end subroutine def_vars_ocean_vars  !]
! ----------------------------------------------------------------------
      subroutine def_vars_rst_ocean_vars  ![
      ! define output variable & attributes in netcdf restart file
      use tracers
#ifdef MARBL
!     add MARBL saved state to restart to restart file
      use marbl_driver, only: marbldrv_create_ss_vars_in_rst
#endif       
      implicit none

      ! local
      integer :: ierr, varid, itrc
      character(len=7)  :: dn_aux = 'auxil'
      
! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);

      varid = nccreate(ncid,'time_step',(/dn_aux,dn_tm/),(/iaux,0/),nf90_int)

      ierr=nf90_put_att (ncid, varid, 'long_name',
     &       'time step and record numbers from initialization')

      varid = nccreate(ncid,'zeta',(/dn_xr,dn_yr,dn_tm/),
     &                                     (/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','free surface elevation')
      ierr = nf90_put_att(ncid,varid,'units','meter')

      varid = nccreate(ncid,'ubar',(/dn_xu,dn_yr,dn_tm/),
     &                                      (/xi_u,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,
     &           'long_name','vertically averaged u-momentum')
      ierr = nf90_put_att(ncid,varid,'units','meter second -1')

      varid = nccreate(ncid,'vbar',(/dn_xr,dn_yv,dn_tm/),
     &                                      (/xi_rho,eta_v,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,
     &           'long_name','vertically averaged v-momentum')
      ierr = nf90_put_att(ncid,varid,'units','meter second -1')

#ifdef MARBL
      call marbldrv_create_ss_vars_in_rst(ncid)
#endif      
      varid = nccreate(ncid,'u',(/dn_xu,dn_yr,dn_zr,dn_tm/),
     &                                     (/xi_u,eta_rho,nz,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','u-momentum')
      ierr = nf90_put_att(ncid,varid,'units','meter second -1')

      varid = nccreate(ncid,'v',(/dn_xr,dn_yv,dn_zr,dn_tm/),
     &                                     (/xi_rho,eta_v,nz,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','v-momentum')
      ierr = nf90_put_att(ncid,varid,'units','meter second -1')

      varid = nccreate(ncid,'temp',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &                              (/xi_rho,eta_rho,nz,0/), nf90_double)
      ierr=nf90_put_att (ncid,varid,'long_name', 'potential temperature')
      ierr=nf90_put_att (ncid,varid,'units', 'degrees Celsius')

      varid = nccreate(ncid,'salt',(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &                              (/xi_rho,eta_rho,nz,0/), nf90_double)
      ierr=nf90_put_att (ncid,varid,'long_name', 'salinity')
      ierr=nf90_put_att (ncid,varid,'units', 'PSU')

      ! tracers beyond T and S
      do itrc=3,nt
        varid = nccreate(ncid,t_vname(itrc),(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &                                         (/xi_rho,eta_rho,N,0/), nf90_double)
        ierr=nf90_put_att (ncid, varid, 'long_name', t_lname(itrc))
        ierr=nf90_put_att (ncid, varid, 'units', t_units(itrc))
      enddo

# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
      varid = nccreate(ncid,'DU_avg1',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/),nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           '<<fast-time averaged uflx>>')

      varid = nccreate(ncid,'DV_avg1',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/),nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           '<<fast-time-averaged vflx>>')

      varid = nccreate(ncid,'DU_avg2',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/),nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           '<<fast-time averaged ubar(:,:,n+1/2)>>')

      varid = nccreate(ncid,'DV_avg2',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/),nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           '<<fast-time-averaged vbar(:,:,n+1/2)>>')

      varid = nccreate(ncid,'DU_avg_bak',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/),nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged ubar(:,:,n-1)>')

      varid = nccreate(ncid,'DV_avg_bak',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/),nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged vbar(:,:,n-1)>')
#  elif defined PRED_COUPLED_MODE
        varid = nccreate(ncid,'rufrc_bak',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/),nf90_double)
        ierr=nf90_put_att(ncid, varid, 'long_name',
     &                        '3D-to-2D forcing, XI-component')
        varid = nccreate(ncid,'rvfrc_bak',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/),nf90_double)
        ierr=nf90_put_att(ncid, varid, 'long_name',
     &                       '3D-to-2D forcing, ETA-component')
#  endif
# endif  /* EXACT_RESTART */

# ifdef LMD_KPP
      varid = nccreate(ncid,'hbls',
     &   (/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','surface BL thickness')
      ierr = nf90_put_att(ncid,varid,'units','meter')
# endif
# ifdef LMD_BKPP
      varid = nccreate(ncid,'hbbl',
     &   (/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','bottom BL thickness')
      ierr = nf90_put_att(ncid,varid,'units','meter')
# endif

      if (diag_pflx) then
        varid = nccreate(ncid,'u_slow',(/dn_xu,dn_yr,dn_zr,dn_tm/),
     &                               (/xi_u,eta_rho,N,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name','time filtered u')
        ierr = nf90_put_att(ncid,varid,'units','m/s')

        varid = nccreate(ncid,'v_slow',(/dn_xr,dn_yv,dn_zr,dn_tm/),
     &                               (/xi_rho,eta_v,N,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name','time filtered v')
        ierr = nf90_put_att(ncid,varid,'units','m/s')

        varid = nccreate(ncid,'p_slow',
     &                       (/dn_xr,dn_yr,dn_zr,dn_tm/),
     &                       (/xi_rho,eta_rho,n,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &                               'time filtered pressure')
        ierr = nf90_put_att(ncid,varid,'units','Pa??')
      endif

      if (obc_tuning) then
        varid = nccreate(ncid,'ub_west',
     &                       (/dn_yr,dn_tm/),
     &                       (/eta_rho,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &                               'west bc tuning coeff')
        varid = nccreate(ncid,'ub_east',
     &                       (/dn_yr,dn_tm/),
     &                       (/eta_rho,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &                               'east bc tuning coeff')
        varid = nccreate(ncid,'ub_north',
     &                       (/dn_xr,dn_tm/),
     &                       (/xi_rho,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &                               'north bc tuning coeff')
        varid = nccreate(ncid,'ub_south',
     &                       (/dn_xr,dn_tm/),
     &                       (/xi_rho,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &                               'south bc tuning coeff')
      endif

      end subroutine def_vars_rst_ocean_vars  !]
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
!     logical,save :: first_step=.true.

      if (wrt_file_avg) then

!     if (.not. first_step) then                           ! don't include t=0 in average

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
        if (wrt_avg_U)   u_avg( 1:i1,j0:j1,:) = u_avg( 1:i1,j0:j1,:)*(1-coef) + u( 1:i1,j0:j1,:,nstp)*coef
        if (wrt_avg_V)   v_avg(i0:i1, 1:j1,:) = v_avg(i0:i1, 1:j1,:)*(1-coef) + v(i0:i1, 1:j1,:,nstp)*coef

        if (wrt_avg_Tem) tem_avg( 1:i1,j0:j1,:) = tem_avg(i0:i1,j0:j1,:)*(1-coef) + t(i0:i1,j0:j1,:,nstp,1)*coef
        if (wrt_avg_Sal) sal_avg(i0:i1, 1:j1,:) = sal_avg(i0:i1,j0:j1,:)*(1-coef) + t(i0:i1,j0:j1,:,nstp,2)*coef

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
        if (wrt_avg_W) then
          ! the w averaging, calculation and output should be checked at some point.
          call wvlcty(work)

          wvl_avg(i0:i1,j0:j1,1:N) = wvl_avg(i0:i1,j0:j1,1:N)*(1-coef)
     &                               + work(i0:i1,j0:j1,1:N)*coef
        endif
        if (wrt_avg_Akv) akv_avg(i0:i1,j0:j1,:) = akv_avg(i0:i1,j0:j1,:)*(1-coef) + akv(i0:i1,j0:j1,:)*coef
        if (wrt_avg_Akt) akt_avg(i0:i1,j0:j1,:) = akt_avg(i0:i1,j0:j1,:)*(1-coef) + akt(i0:i1,j0:j1,:)*coef
#  ifdef LMD_KPP
        if (wrt_avg_Hbls) hbls_avg(i0:i1,j0:j1)  =  hbls_avg(i0:i1,j0:j1)*(1-coef) + hbls(i0:i1,j0:j1)*coef
#  endif
#  ifdef LMD_BKPP
        if (wrt_avg_Hbbl) hbbl_avg(i0:i1,j0:j1) = hbbl_avg(i0:i1,j0:j1)*(1-coef) + hbbl(i0:i1,j0:j1)*coef
#  endif

!     endif  ! <-- .not. first_step

!     first_step=.false.

      endif ! <-- wrt_file_avg

      end subroutine calc_avg_ocean_vars  !]
! ----------------------------------------------------------------------
      subroutine write_time_step( ncid, record, total_recs )  ![
      implicit none

      ! inputs
      integer, intent(in) :: ncid
      integer, intent(in) :: record                        ! current file record number
      integer, intent(in) :: total_recs                    ! total records for variable

      ! local
      integer :: ibuff(iaux), start(2), count(2)           ! iaux = 6 from wrt_his.F
      integer :: var_id_tmp, ierr

      ibuff(1)=iic    ; ibuff(2)=record                    ! time step and
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
      implicit none

      ! local
      integer :: itrc

      if (mynode==0) write(*,'(/7x,2A,F6.1,2x,A,I4)')                ! write to terminal output in simulation pre-amble text which
     &  'ocean_vars :: history file ',                               ! result variables are being stored
     &  'ouput_period =', output_period_his,
     &  'recs/file =', nrpf_his

      if (mynode==0) write(*,'(9x,A,13(/11x,A,T20,L1,4x,A),/)')      ! 13(....) repeats formatting 11 times.
     &               'his fields to be saved: (T/F)'                 ! T20 moves to the 20th character on line
     &                 , 'zeta',   wrt_Z,    'free surface elevation'
     &                 , 'ubar',   wrt_Ub,   'vertically averaged u-momentum'
     &                 , 'vbar',   wrt_Vb,   'vertically averaged v-momentum'
     &                 , 'u',      wrt_U,    'u-momentum'
     &                 , 'v',      wrt_V,    'v-momentum'
     &                 , 'temp',   wrt_Tem,  'Temperature'
     &                 , 'salt',   wrt_Sal,  'Salinity'
     &                 , 'rho',    wrt_R,    'density anomaly'
     &                 , 'Omega',  wrt_O,    'S-coordinate vertical flux'
     &                 , 'W',      wrt_W,    'vertical velocity'

     &                 , 'Akv',    wrt_Akv,  'vertical momentum mixing'
     &                 , 'Akt',    wrt_Akt,  'vertical tracer mixing'
# ifdef LMD_KPP
     &                 , 'hbls',   wrt_Hbls, 'surface bl thickness'
# endif
# ifdef LMD_BKPP
     &                 , 'hbbl',   wrt_Hbbl, 'bottom bl thickness'
# endif
          do itrc=3,2+nt_passive
            if (mynode==0) write(*,'(11x,A,I2,A,T16,L1,T20,A)') 't(',
     &        itrc, ')', wrt_t(itrc), t_vname(itrc)
          enddo

      end subroutine display_output_settings_to_terminal  !]
!-----------------------------------------------------------------------
      subroutine display_output_settings_to_terminal_avg  ![
      use tracers
      implicit none

      ! local
      integer :: itrc

      if (mynode==0) write(*,'(/7x,2A,F6.1,2x,A,I4)')                ! write to terminal output in simulation pre-amble text which
     &  'ocean_vars :: average file ',                               ! result variables are being stored
     &  'ouput_period =', output_period_avg,
     &  'recs/file =', nrpf_avg

      if (mynode==0) write(*,'(9x,A,13(/11x,A,T20,L1,4x,A),/)')      ! 13(....) repeats formatting 11 times.
     &               'avg fields to be saved: (T/F)'                 ! T20 moves to the 20th character on line
     &                 , 'zeta',   wrt_Z,    'free surface elevation'
     &                 , 'ubar',   wrt_Ub,   'vertically averaged u-momentum'
     &                 , 'vbar',   wrt_Vb,   'vertically averaged v-momentum'
     &                 , 'u',      wrt_U,    'u-momentum'
     &                 , 'v',      wrt_V,    'v-momentum'
     &                 , 'temp',   wrt_Tem,  'Temperature'
     &                 , 'salt',   wrt_Sal,  'Salinity'
     &                 , 'rho',    wrt_R,    'density anomaly'
     &                 , 'Omega',  wrt_O,    'S-coordinate vertical flux'
     &                 , 'W',      wrt_W,    'vertical velocity'

     &                 , 'Akv',    wrt_Akv,  'vertical momentum mixing'
     &                 , 'Akt',    wrt_Akt,  'vertical tracer mixing'
# ifdef LMD_KPP
     &                 , 'hbls',   wrt_hbls, 'surface bl thickness'
# endif
# ifdef LMD_BKPP
     &                 , 'hbbl',   wrt_hbbl, 'bottom bl thickness'
# endif
          do itrc=3,2+nt_passive
            if (mynode==0) write(*,'(11x,A,I2,A,T16,L1,T20,A)') 't(',
     &        itrc, ')', wrt_t_avg(itrc), t_vname(itrc)
          enddo

      end subroutine display_output_settings_to_terminal_avg  !]
! ----------------------------------------------------------------------

      end module basic_output
