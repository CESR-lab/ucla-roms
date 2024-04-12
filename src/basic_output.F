
      module basic_output
      ! Basic output routines for ocean variables
      ! his, avg, and rst files

#include "cppdefs.opt"

      use param
      use netcdf
      use nc_read_write
      use roms_read_write
      use sponge_tune
      use diagnostics
      use scalars
      use ocean_vars
      implicit none
      private

#include "ocean_vars.opt"

# ifdef SOLVE3D
      real,allocatable,dimension(:,:,:) :: akv_avg
      real,allocatable,dimension(:,:,:) :: akt_avg
#  ifdef SALINITY
      real,allocatable,dimension(:,:,:) :: aks_avg
#  endif
#  ifdef LMD_KPP
      real,allocatable,dimension(:,:) :: hbl_avg
#  endif
#  ifdef LMD_BKPP
      real,allocatable,dimension(:,:) :: hbbl_avg
#  endif
      real,allocatable,dimension(:,:,:) :: rho_avg
# endif

      ![   Taken from old ncvars (should be improved eventually):
      integer, parameter, public :: indxTime=1, indxZ=2, indxUb=3, indxVb=4
#ifdef SOLVE3D
     &                    , indxU=5, indxV=6, indxO=7, indxW=8
     &                    , indxR=9, indxT=10
# ifdef SALINITY
     &                    , indxS=indxT+1
# endif
     &                    , indxAkv=indxT+NT,   indxAkt=indxAkv+1
# ifdef SALINITY
     &                    , indxAks=indxAkt+1
# endif
# ifdef LMD_KPP
#  ifdef SALINITY
     &                    , indxHbls=indxAks+1
#  else
     &                    , indxHbls=indxAkt+1
#  endif
# ifdef LMD_BKPP
     &                    , indxHbbl=indxHbls+1
# endif
# else
# ifdef LMD_BKPP
     &                    , indxHbbl=indxAks+1
# endif
# endif

#endif
      !]

      character*42, public ::  vname(3, 180) ! DevinD 180 copied from old code only
      integer, parameter, public :: iaux=6           ! Length of netCDF variable "time_step"

      ! end old ncvars

      ! netcdf outputting:
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_ovars=0
      integer :: navg_ovars = 0                             ! number of samples in average
      character(len=99),public :: fname_rst                 ! restart filename needed by diagnostics
      integer,public           :: rec_rst=nrpf_rst          ! current file output record needed by diags

      public :: wrt_his_ocean_vars
      public :: wrt_avg_ocean_vars
      public :: calc_avg_ocean_vars
      public :: wrt_rst_ocean_vars
      public :: init_avg_arrays

      contains  !]
!----------------------------------------------------------------------
      subroutine init_avg_arrays  ![
      use scalars

      implicit none

      if (wrt_file_his .and. mod(output_period_his,dt) /= 0) then
        stop 'ERROR :: ocean_vars.opt - dt not a factor of output_period_his!'
      endif
      if (wrt_file_avg .and. mod(output_period_avg,dt) /= 0) then
        stop 'ERROR :: ocean_vars.opt - dt not a factor of output_period_avg!'
      endif

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

! Names of variables in NetCDF output files. The first element
! is the name of the variable; the other two are are attributes.

      vname(1,indxTime)='ocean_time'
      vname(2,indxTime)='Time since initialization'
      vname(3,indxTime)='second'

      vname(1,indxZ)='zeta'
      vname(2,indxZ)='free-surface elevation'
      vname(3,indxZ)='meter'

      vname(1,indxUb)='ubar'
      vname(2,indxUb)='vertically averaged u-momentum component'
      vname(3,indxUb)='meter second-1'

      vname(1,indxVb)='vbar'
      vname(2,indxVb)='vertically averaged v-momentum component'
      vname(3,indxVb)='meter second-1'

#ifdef SOLVE3D
      vname(1,indxU)='u'
      vname(2,indxU)='u-momentum component'
      vname(3,indxU)='meter second-1'

      vname(1,indxV)='v'
      vname(2,indxV)='v-momentum component'
      vname(3,indxV)='meter second-1'

      vname(1,indxO)='omega'
      vname(2,indxO)='S-coordinate vertical velocity'
      vname(3,indxO)='meter second-1'

      vname(1,indxW)='w'
      vname(2,indxW)='vertical velocity'
      vname(3,indxW)='meter second-1'

      vname(1,indxR)='rho'
      vname(2,indxR)='density anomaly'
      vname(3,indxR)='kilogram meter-3'

      vname(1,indxT)='temp'
      vname(2,indxT)='potential temperature'
      vname(3,indxT)='Celsius'

# ifdef SALINITY
      vname(1,indxS)='salt'
      vname(2,indxS)='salinity'
      vname(3,indxS)='PSU'
# endif

      vname(1,indxAkv)='Akv'
      vname(2,indxAkv)='vertical viscosity coefficient'
      vname(3,indxAkv)='meter2 second-1'

      vname(1,indxAkt)='Akt'
      vname(2,indxAkt)='vertical thermal conductivity coefficient'
      vname(3,indxAkt)='meter2 second-1'
# ifdef SALINITY
      vname(1,indxAks)='AKs'
      vname(2,indxAks)='salinity vertical diffusion coefficient'
      vname(3,indxAks)='meter2 second-1'
# endif
# ifdef LMD_KPP
      vname(1,indxHbls)='hbls'
      vname(2,indxHbls)='Thickness of KPP surface boundary layer'
      vname(3,indxHbls)='meter'
# endif
# ifdef LMD_BKPP
      vname(1,indxHbbl)='hbbl'
      vname(2,indxHbbl)='Thickness of KPP bottom boundary layer'
      vname(3,indxHbbl)='meter'
# endif
#endif

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

      use mixing
      use eos_vars
      use work_mod
      use grid
      use tracers    !, only: wrt_his_trc

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
        if (wrt_Z)  call ncwrite(ncid,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start)
        if (wrt_Ub) call ncwrite(ncid,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start)
        if (wrt_Vb) call ncwrite(ncid,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start)

#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_his called at the middle of next timestep as some variables only calculated there for t=n.

        start(3)=1; start(4)=rec_his
        if (wrt_U) call ncwrite(ncid,vname(1,indxU),u( 1:i1,j0:j1,:,nnew),start)
        if (wrt_V) call ncwrite(ncid,vname(1,indxV),v(i0:i1, 1:j1,:,nnew),start)

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
                ! Implicit omega only
!               work(i,j,k)=pm(i,j)*pn(i,j)*Wi(i,j,k)
              enddo
            enddo
          enddo
          call ncwrite(ncid,vname(1,indxO),work(i0:i1,j0:j1,:),start)
        endif
        if (wrt_W) then                              ! true vertical velocity (m/s).
          do tile=0,NSUB_X*NSUB_E-1
            call wvlcty (tile, work)
          enddo
          call ncwrite(ncid,vname(1,indxW),work(i0:i1,j0:j1,1:N),start)
        endif
        if (wrt_Akv) call ncwrite(ncid,vname(1,indxAkv),Akv(i0:i1,j0:j1,:),start)
        if (wrt_Akt) call ncwrite(ncid,vname(1,indxAkt),Akt(i0:i1,j0:j1,:,itemp),start)
# ifdef SALINITY
        if (wrt_Aks) call ncwrite(ncid,vname(1,indxAks),Akt(i0:i1,j0:j1,:,isalt),start)
# endif
        start(3)=rec_his                                      ! back to 2D vars
# ifdef LMD_KPP
        if (wrt_Hbls) call ncwrite(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1),start)
# endif
# ifdef LMD_BKPP
        if (wrt_Hbbl) call ncwrite(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1),start)
# endif
#endif
        ierr=nf90_close(ncid)
        if (mynode == 0 .and. (.not.code_check)) then                 ! supress output for code_check
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'ocean_vars :: wrote history, tdays =', tdays,
     &     'step =', iic, 'rec =', rec_his
        endif
        code_check=.false.

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

      use eos_vars
      use mixing
      use tracers
      use grid
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_avg=nrpf_avg
      integer,save           :: total_rec_avg=0                      ! total avg output records so far
      real,save              :: output_time_avg=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_avg
      integer :: tile, tn, ierr, k

      if (wrt_file_avg) then

        output_time_avg = output_time_avg + dt

        if (output_time_avg>=output_period_avg) then

          if (rec_avg == nrpf_avg) then
            call create_file_ocean_vars(fname_avg,.true.)
            rec_avg = 0
          endif
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
        if (wrt_avg_W)   call ncwrite(ncid, vname(1,indxW),   wvl_avg(i0:i1,j0:j1,1:N), start)
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

      use mixing
      use eos_vars
      use work_mod
      use grid
      use tracers 
      use coupling
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

      if (wrt_file_rst) then

!     if (.not. first_step) output_time_rst = output_time_rst + dt   ! only start count after first timestep
!     first_step=.false.                                             ! as first step the ocean_vars values are for t=0
      output_time_rst = output_time_rst + dt

      if (output_time_rst    >=output_period_rst
#ifdef EXACT_RESTART
     &    .or.
     &    output_time_rst+dt >=output_period_rst                     ! timestep before
#endif
     &                    ) then


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
        call ncwrite(ncid,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start)
        call ncwrite(ncid,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start)
        call ncwrite(ncid,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start)
#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_rst called at the middle of next timestep as some variables only calculated there for t=n.

        start(3)=1; start(4)=rec_rst
        call ncwrite(ncid,vname(1,indxU),u( 1:i1,j0:j1,:,nnew),start)
        call ncwrite(ncid,vname(1,indxV),v(i0:i1, 1:j1,:,nnew),start)

        call wrt_rst_trc(ncid,start)                              ! tracer variables

        start(3)=rec_rst                                      ! back to 2D vars

# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
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
        call ncwrite(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1),start)
# endif
# ifdef LMD_BKPP
        call ncwrite(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1),start)
# endif
#endif
        ! Added for perfect restart reasons
!       start(3)=1; start(4)=rec_rst
!       call ncwrite(ncid,'Akv',Akv(i0:i1,j0:j1,:),start)
!       call ncwrite(ncid,'Akt',Akt(i0:i1,j0:j1,:,itemp),start)
!       call ncwrite(ncid,'Aks',Akt(i0:i1,j0:j1,:,isalt),start)

        if (ub_tune)   call wrt_rst_ub(ncid,rec_rst)
        if (diag_pflx) call wrt_rst_diag_slow(ncid,rec_rst)

#ifdef MARBL
        call marbldrv_write_ss_vars_to_rst(ncid,rec_rst)
#endif         

        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'ocean_vars :: wrote restart, tdays =', tdays,
     &     'step =', iic, 'rec =', rec_rst, '/', total_rec_rst
        endif

        if (output_time_rst>=output_period_rst) then
          output_time_rst=0                                          ! catch for exact restart to not
        endif                                                        ! reset the time
      endif

      endif  ! <-- wrt_file_rst

      end subroutine wrt_rst_ocean_vars  !]
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
        if (.not. avg) text_lname=vname(2,indxZ)
        if (      avg) text_lname='averaged '/ /vname(2,indxZ)
        varid = nccreate(ncid,vname(1,indxZ),(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxZ))
      endif
      if ( (wrt_Ub .and. .not. avg) .or. (wrt_avg_Ub .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxUb)
        if (      avg) text_lname='averaged '/ /vname(2,indxUb)
        varid = nccreate(ncid,vname(1,indxUb),(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxUb))
      endif
      if ( (wrt_Vb .and. .not. avg) .or. (wrt_avg_Vb .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxVb)
        if (      avg) text_lname='averaged '/ /vname(2,indxVb)
        varid = nccreate(ncid,vname(1,indxVb),(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxVb))
      endif
#ifdef SOLVE3D
      if ( (wrt_U .and. .not. avg) .or. (wrt_avg_U .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxU)
        if (      avg) text_lname='averaged '/ /vname(2,indxU)
        varid = nccreate(ncid,vname(1,indxU),(/dn_xu,dn_yr,dn_zr,dn_tm/),(/xi_u,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxU))
      endif
      if ( (wrt_V .and. .not. avg) .or. (wrt_avg_V .and. avg) ) then
        if (.not. avg) text_lname=vname(2,indxV)
        if (      avg) text_lname='averaged '/ /vname(2,indxV)
        varid = nccreate(ncid,vname(1,indxV),(/dn_xr,dn_yv,dn_zr,dn_tm/),(/xi_rho,eta_v,N,0/))
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
        varid = 
     &   nccreate(ncid,vname(1,indxR),(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,n,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxR))
      endif
      if ( (wrt_O .and. .not. avg) .or. (wrt_avg_O .and. avg) ) then    ! s-coordinate "omega" vertical velocity.
        if (.not. avg) text_lname=vname(2,indxO)
        if (      avg) text_lname='averaged '/ /vname(2,indxO)
        varid = 
     &    nccreate(ncid,vname(1,indxO),(/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,n+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxO))
      endif
      if ( (wrt_W .and. .not. avg) .or. (wrt_avg_W .and. avg) ) then   ! true W-vertical velocity.
        if (.not. avg) text_lname=vname(2,indxW)
        if (      avg) text_lname='averaged '/ /vname(2,indxW)
        varid = 
     &    nccreate(ncid,vname(1,indxW),(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,n,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxW))
      endif
      if ( (wrt_Akv .and. .not. avg) .or. (wrt_avg_Akv .and. avg)) then ! vertical viscosity coefficient.
        if (.not. avg) text_lname=vname(2,indxAkv)
        if (      avg) text_lname='averaged '/ /vname(2,indxAkv)
        varid = nccreate(ncid,vname(1,indxAkv),
     &    (/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,n+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkv))
      endif
      if ( (wrt_Akt .and. .not. avg) .or. (wrt_avg_Akt .and. avg) ) then  ! vertical thermal conductivity coefficient.
        if (.not. avg) text_lname=vname(2,indxAkt)
        if (      avg) text_lname='averaged '/ /vname(2,indxAkt)
        varid = nccreate(ncid,vname(1,indxAkt),
     &      (/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAkt))
      endif
# ifdef SALINITY
      if ( (wrt_Aks .and. .not. avg) .or. (wrt_avg_Aks .and. avg) ) then! vertical diffusion coefficient for salinity.
        if (.not. avg) text_lname=vname(2,indxAks)
        if (      avg) text_lname='averaged '/ /vname(2,indxAks)
        varid = 
     &    nccreate(ncid,vname(1,indxAks),(/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxAks))
      endif
# endif /* SALINITY */
# ifdef LMD_KPP
      if ( (wrt_Hbls .and. .not. avg) .or. (wrt_avg_Hbls .and. avg) ) then                                   ! depth of surface planetary boundary layer.
        if (.not. avg) text_lname=vname(2,indxHbls)
        if (      avg) text_lname='averaged '/ /vname(2,indxHbls)
        varid = nccreate(ncid,vname(1,indxHbls),(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbls))
      endif
# endif
# ifdef LMD_BKPP
      if ( (wrt_Hbbl .and. .not. avg) .or. (wrt_avg_Hbbl .and. avg) ) then       ! thickness of bottom boundary layer.
        if (.not. avg) text_lname=vname(2,indxHbbl)
        if (      avg) text_lname='averaged '/ /vname(2,indxHbbl)
        varid = nccreate(ncid,vname(1,indxHbbl),
     &         (/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
        ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbbl))
      endif
# endif
#endif /* SOLVE3D */

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

      varid = nccreate(ncid,vname(1,indxZ),(/dn_xr,dn_yr,dn_tm/),
     &                                     (/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxZ))
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxZ))

      varid = nccreate(ncid,vname(1,indxUb),(/dn_xu,dn_yr,dn_tm/),
     &                                      (/xi_u,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxUb))
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxUb))

      varid = nccreate(ncid,vname(1,indxVb),(/dn_xr,dn_yv,dn_tm/),
     &                                      (/xi_rho,eta_v,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxVb))
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxVb))
#ifdef MARBL
      call marbldrv_create_ss_vars_in_rst(ncid)
#endif      
#ifdef SOLVE3D
      varid = nccreate(ncid,vname(1,indxU),(/dn_xu,dn_yr,dn_zr,dn_tm/),
     &                                     (/xi_u,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxU))
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxU))

      varid = nccreate(ncid,vname(1,indxV),(/dn_xr,dn_yv,dn_zr,dn_tm/),
     &                                     (/xi_rho,eta_v,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxV))
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxV))

      do itrc=1,NT
        varid = nccreate(ncid,t_vname(itrc),(/dn_xr,dn_yr,dn_zr,dn_tm/),
     &                                         (/xi_rho,eta_rho,N,0/), nf90_double)
        ierr=nf90_put_att (ncid, varid, 'long_name', t_lname(itrc))
        ierr=nf90_put_att (ncid, varid, 'units', t_units(itrc))
      enddo

# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
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
      varid = nccreate(ncid,vname(1,indxHbls),
     &   (/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxHbls))
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbls))
# endif
# ifdef LMD_BKPP
      varid = nccreate(ncid,vname(1,indxHbbl),(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',vname(2,indxHbbl))
      ierr = nf90_put_att(ncid,varid,'units',vname(3,indxHbbl))
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
#ifdef SPONGE_TUNE
       if (ub_tune) then
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
#endif

#endif /* SOLVE3D */

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
        if (wrt_avg_W) then
          ! the w averaging, calculation and output should be checked at some point.
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
        if (wrt_avg_Hbls) hbl_avg(i0:i1,j0:j1)  =  hbl_avg(i0:i1,j0:j1)*(1-coef) + hbls(i0:i1,j0:j1)*coef
#  endif
#  ifdef LMD_BKPP
        if (wrt_avg_Hbbl) hbbl_avg(i0:i1,j0:j1) = hbbl_avg(i0:i1,j0:j1)*(1-coef) + hbbl(i0:i1,j0:j1)*coef
#  endif
# endif /* SOLVE3D */

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
! ----------------------------------------------------------------------

      end module basic_output
