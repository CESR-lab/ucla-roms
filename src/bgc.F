      module bgc

      ! BGC - bio-geo-chemical module
      ! -----------------------------

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Dec)
      ! Refactoring of ETH Zurich roms code with BEC2, code which has been used
      ! by Pierre Damien to run full pacific models.

#include "cppdefs.opt"
#ifdef BIOLOGY_BEC2

      ! param needed for GLOBAL_2D_array to work. NT = number tracer from param
      use param
      use read_write, only: set_frc_var_tile, ncvs
      use tracers     ! for iPO4, etc, indices of bgc tracers (formerly in param.h)
      use ncvars
      use bgc_ecosys_vars
      use nc_read_write
      use roms_read_write
      use netcdf
      use scalars

      implicit none

      private ! Make all variable private to this module unless public specified

      ! Includes:
!#include "bgc_ecosys_bec2.h"

#include "bgc.opt"


#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
      real, public    :: srflx_dailyavg(GLOBAL_2D_ARRAY)   ! daily average short-wave rad
#endif

      ! netcdf outputting:
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_bgc=0
      integer :: navg_bgc = 0                              ! number of samples in average

      ! misc from porting ETH code:
      character(len=40) tclm_name(NT) ! Taken from ETH code ncvars

      public set_bgc_surf_frc
#ifdef BEC2_DIAG
!      public def_his_bgc_diag
!      public wrt_his_bgc_diag
# ifdef AVERAGES
!      public def_avg_bgc_diag
!      public wrt_avg_bgc_diag
      public set_avg_bgc_diag
# endif
#endif
      public init_arrays_bgc_frc
      public init_scalars_bec2
      public wrt_his_bgc_tracers
      public wrt_avg_bgc_tracers

      contains

! ----------------------------------------------------------------------
      subroutine set_bgc_surf_frc(istr,iend,jstr,jend)  ![

      ! read in bgc surface flux and interpolate to model time.
      ! Taken from get_smth.F & set_smth.F of ETH code.

      ! Since bgc is set up such that iron and dust have flux at surface but
      ! are not directly bgc tracers, I have separated them into set_bgc_frc

      ! this routine is called by # if defined BIOLOGY_BEC2 in set_forces.F of ETH code

      use scalars
      use bgc_forces

      implicit none

! Needed for iron and dust variables, not yet ported to this module. Should do.
!#include "bgc_forces.h"

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      
      ! local
      integer :: i, j

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
# include "compute_extended_bounds.h"
#endif


      call set_frc_var_tile( istr, iend, jstr, jend, nc_dust,    dust,    interp_frc )  ! dust
      call set_frc_var_tile( istr, iend, jstr, jend, nc_iron,    iron,    interp_frc )  ! iron


#ifdef PCO2AIR_FORCING
      call set_frc_var_tile( istr, iend, jstr, jend, nc_pco2air, pco2air, interp_frc )
#endif /* PCO2AIR_FORCING */   


#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
      ! Daily average Short wave radiation
      call set_frc_var_tile( istr, iend, jstr, jend, nc_swrad_avg, srflx_dailyavg, interp_frc )
     
      do j=jstrR,jendR
        do i=istrR,iendR
          srflx_dailyavg(i,j)= srflx_dailyavg(i,j)/(rho0*Cp)
          
#ifdef SEA_ICE_NOFLUX
          if( t(i,j,N,nrhs,itemp) .le. -1.8 ) then    ! Restrict stflx to prevent surface temperature to go below -2 degrees C.
#   if defined LMD_KPP
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
            srflx_dailyavg(i,j)=0.
#endif                   
#   endif
          endif
#endif

        enddo !<- istrR
      enddo   !<- jstrR
     
#endif /* DAILYPAR_PHOTOINHIBITION */ 

      end subroutine set_bgc_surf_frc  !]

! ----------------------------------------------------------------------
      subroutine wrt_his_bgc_tracers  ![
      ! write variables to output netcdf file
      ! bgc variables are in fact tracers, but we choose to output
      ! them to their own file.
      ! this file also include bgc diagnostic variables.
      ! selection of variables to output is still done in tracers.opt
      ! otherwise would cause circular reference as needed in both directions.
      !
      ! TIME: since bgc is calculated in step3d_t_ISO.F and uses 'nnew'
      ! for the tracers (i.e. t=n+1). The bgc tracers are for t=n+1.
      ! since we now write them out immediately, we need to set the time to
      ! t=n+1, i.e. time+dt.
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_his                              ! current file output record
      integer,save           :: total_rec_his=0                      ! total his output records so far
      real,save              :: output_time_his=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_his
      integer                :: itrc, ierr

      ! check this as different from WEC first_step and time. Called from step3d_t_ISO
      if (.not. first_step) output_time_his = output_time_his + dt   ! only start count after first timestep
      first_step=.false.                                             ! as first step the wec values are for t=0

      if (output_time_his>=output_period_his .or.
     &    output_time_his==0                      ) then             ! catch to get t=0 output

        if (mod(total_rec_his,nrpf_his)==0) then
          call create_bgc_file(total_rec_his,fname_his,.false.)
          rec_his = 0
        endif
        total_rec_his = total_rec_his +1
        rec_his = rec_his + 1

        ierr=nf90_open(fname_his,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time+dt/),(/rec_his/))      ! t=n+1 thus time+dt

        start=1; start(4)=rec_his                                    ! 3D vars
        do itrc=iTandS+nt_passive+1, NT
          if (wrt_t(itrc)) call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nstp,itrc),start)
        enddo

        call wrt_his_bgc_diag(ncid,rec_his)

        ierr=nf90_close(ncid)
        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'bgc :: wrote history, tdays =', tdays,
     &     'step =', iic, 'rec =', rec_his, '/', total_rec_his
        endif

        output_time_his=0
      endif

      end subroutine wrt_his_bgc_tracers  !]

! ----------------------------------------------------------------------
      subroutine wrt_avg_bgc_tracers  ![
      ! write averaged variables to output netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      implicit none

      ! local
      integer,dimension(4)   :: start
      integer,save           :: rec_avg                              ! current file output record
      integer,save           :: total_rec_avg=0                      ! total avg output records so far
      real,save              :: output_time_avg=0                    ! time since last output
      logical,save           :: first_step=.true.
      character(len=99),save :: fname_avg
      integer                :: itrc, ierr, itavg

      if (.not. first_step) then
        call calc_avg_bgc                                            ! don't include t=0 in averaging
        output_time_avg = output_time_avg + dt                       ! only start count after first timestep
      endif
      first_step=.false.                                             ! as first step the bgc values are for t=0

      if (output_time_avg>=output_period_avg) then                   ! catch to get t=0 output

        if (mod(total_rec_avg,nrpf_avg)==0) then
          call create_bgc_file(total_rec_avg,fname_avg,.true.)
          rec_avg = 0
        endif
        total_rec_avg = total_rec_avg +1
        rec_avg = rec_avg + 1

        ierr=nf90_open(fname_avg,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/t_avg_bgc/),(/rec_avg/))

        start=1; start(4)=rec_avg                                    ! back to 3D vars

        do itrc=iTandS+nt_passive+1, NT
          if (wrt_t_avg(itrc)) then
            itavg = NT_2_t_avg(itrc)                                 ! get respective index for t_avg(itavg) -> t(itrc)
            call ncwrite(ncid,t_vname(itrc),t_avg(i0:i1,j0:j1,:,itavg),start)
          endif
        enddo

        call wrt_avg_bgc_diag(ncid,rec_avg)

        navg_bgc=0
        output_time_avg=0

        ierr=nf90_close(ncid)
        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'bgc :: wrote averages, tdays =', tdays,
     &        'step =', iic-1, 'rec =', rec_avg, '/', total_rec_avg  ! -1 since iic=t(n+1)
     &         MYID
        endif
      endif

      end subroutine wrt_avg_bgc_tracers  !]

!----------------------------------------------------------------------
      subroutine create_bgc_file(total_rec,fname,avg)  ![
      implicit none

      !input/output
      integer          ,intent(in)  :: total_rec           ! total netcdf records so far for these variables
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ierr,varid
      character(len=10),dimension(4) :: dimnames           ! dimension names
      integer,          dimension(4) :: dimsizes

      if (      avg) fname=trim(output_root_name) / / '_bgc_avg.nc'
      if (.not. avg) fname=trim(output_root_name) / / '_bgc.nc'
      call create_filename(fname,total_rec)

      ierr=nf90_create(fname,nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) call handle_ierr(ierr,'create_bgc_file: ')

      varid = nccreate(ncid,'ocean_time',(/'time'/),(/0/),nf90_double)
      if (      avg) ierr = nf90_put_att(ncid,varid,'long_name','Average time since origin')
      if (.not. avg) ierr = nf90_put_att(ncid,varid,'long_name','Time since origin')
      ierr = nf90_put_att(ncid,varid,'units','second' )

      call put_global_atts(ncid, ierr)                     ! put global attributes in file
      if (ierr/=nf90_noerr)
     &  call handle_ierr(ierr,'create_bgc_file: ')

      if (      avg) call def_vars_bgc( ncid, .true.  )
      if (.not. avg) call def_vars_bgc( ncid, .false. )

      ierr = nf90_enddef(ncid)
      if (mynode == 0) then
        write(*,'(7x,2A)')
     &    'bgc :: created new netCDF file ', trim(fname)
      endif

      end subroutine create_bgc_file !]

! ----------------------------------------------------------------------
      subroutine def_vars_bgc( ncid, avg )  ![
      ! define output variable & attributes in netcdf results file
      ! for actual or averaged variables
      implicit none

      ! input
      integer,intent(in) :: ncid
      logical,intent(in) :: avg
      ! local
      integer           :: ierr, varid, itrc
      character(len=64) :: text_lname

      do itrc=iTandS+nt_passive+1, NT
        if (wrt_t(itrc)) then
          if (avg) then
            text_lname='averaged '/ /t_lname(itrc)
          else
            text_lname=t_lname(itrc)
          endif
          varid = nccreate(ncid,t_vname(itrc),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
          ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
        endif
      enddo

      if (avg) then
        call def_avg_bgc_diag(ncid)
      else
        call def_his_bgc_diag(ncid)
      endif


      end subroutine def_vars_bgc  !]

! ----------------------------------------------------------------------
#ifdef AVERAGES
      subroutine calc_avg_bgc  ![
      implicit none

      ! local
      real    :: coef
      integer :: itrc, itavg

      navg_bgc = navg_bgc +1
      coef = 1./navg_bgc

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'bgc :: started averaging. ',
     &   'output_period_avg (s) =', output_period_avg
      endif

      t_avg_bgc = t_avg_bgc*(1-coef) + time*coef

      ! need i0:i1 indices because arrays still GLOBAL_2D therefore wasted margin
      do itrc=iTandS+nt_passive+1, NT
        if (wrt_t_avg(itrc)) then
          itavg = NT_2_t_avg(itrc)                         ! get respective index for t_avg(itavg) -> t(itrc)
          t_avg(i0:i1,j0:j1,:,itavg) = t_avg(i0:i1,j0:j1,:,itavg)    *(1-coef)
     &                               +     t(i0:i1,j0:j1,:,nstp,itrc)*coef     ! CONFIRM NSTP OR NNEW!!!
!                t_avg(i,j,k,itavg)=cff1*t_avg(i,j,k,itavg)
!     &                             +cff*t(i,j,k,nstp,itrc)
        endif
      enddo

      end subroutine calc_avg_bgc  !]
#endif /* AVERAGES */

! ----------------------------------------------------------------------
#ifdef BEC2_DIAG
      subroutine def_his_bgc_diag(ncid)  ![

      ! Define history/avg file variables for BEC2_DIAG
      ! Taken from def_his.F of ETH code.
      ! Called from def_his.F

      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer :: idiag, ierr, varid

      ! 2d diagnostics:
      do idiag=1,nr_bec2_diag_2d
        if (wrtavg_bec2_diag_2d(idiag)) then
          varid = nccreate(ncid,vname_bec2_diag_2d(1,idiag),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',vname_bec2_diag_2d(2,idiag))
          ierr = nf90_put_att(ncid,varid,'units',vname_bec2_diag_2d(3,idiag))
        endif
      enddo

      ! 3d diagnostics:
      do idiag=1,nr_bec2_diag_3d
        if (wrtavg_bec2_diag_3d(idiag)) then
          varid = nccreate(ncid,vname_bec2_diag_3d(1,idiag),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',vname_bec2_diag_3d(2,idiag))
          ierr = nf90_put_att(ncid,varid,'units',vname_bec2_diag_3d(3,idiag))
        endif
      enddo

      end subroutine def_his_bgc_diag  !]
#endif /* BEC2_DIAG */

! ----------------------------------------------------------------------
! Need AVERAGES as 'ncavg' only exists in AVERAGES is defined.
#ifdef AVERAGES
#ifdef BEC2_DIAG
      subroutine def_avg_bgc_diag(ncid)  ![
      ! Define history/avg file variables for BEC2_DIAG
      ! Taken from def_his.F of ETH code.
      ! Called from def_his.F
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer idiag, ierr, varid
      character(len=64) long_name

      ! 2d diagnostics:
      do idiag=1,nr_bec2_diag_2d
        if (wrtavg_bec2_diag_2d(idiag)) then
          long_name='averaged '/ /vname_bec2_diag_2d(2,idiag) ! Add averaged to long names
          varid = nccreate(ncid,vname_bec2_diag_2d(1,idiag),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',long_name)
          ierr = nf90_put_att(ncid,varid,'units',vname_bec2_diag_2d(3,idiag))
        endif
      enddo

      ! 3d diagnostics:
      do idiag=1,nr_bec2_diag_3d
        if (wrtavg_bec2_diag_3d(idiag)) then
          long_name='averaged '/ /vname_bec2_diag_3d(2,idiag) ! Add averaged to long name
!          call nc_define_var(ncid, vname_bec2_diag_3d(1,idiag),
!     &    long_name, vname_bec2_diag_3d(3,idiag), r3dgrd, ierr )
          varid = nccreate(ncid,vname_bec2_diag_3d(1,idiag),(/'xi_rho','eta_rho','s_rho','time'/),(/xi_rho,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',long_name)
          ierr = nf90_put_att(ncid,varid,'units',vname_bec2_diag_3d(3,idiag))
        endif
      enddo

      end subroutine def_avg_bgc_diag  !]
#endif /* BEC2_DIAG */
#endif /* AVERAGES */

! ----------------------------------------------------------------------
#ifdef BEC2_DIAG
      subroutine wrt_his_bgc_diag (ncid,record)  ![
      ! Write BEC2_DIAG to history file
      implicit none

      ! inputs
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      ! local
      integer :: itrc, ierr

      ! 2d diagnostics:
      do itrc=1,nr_bec2_diag_2d
         if (wrtavg_bec2_diag_2d(itrc)) then
           call nc_write_var(ncid, bec2_diag_2d(START_2D_ARRAY,itrc), 1,
     &             vname_bec2_diag_2d(1,itrc), rp_var, record, ierr )
         endif
      enddo

      ! 3d diagnostics:
      do itrc=1,nr_bec2_diag_3d
         if (wrtavg_bec2_diag_3d(itrc)) then
           call nc_write_var(ncid, bec2_diag_3d(START_2D_ARRAY,1,itrc), N,
     &             vname_bec2_diag_3d(1,itrc), rp_var, record, ierr )
         endif
      enddo

      end subroutine wrt_his_bgc_diag  !]
#endif /* BEC2_DIAG */

! ----------------------------------------------------------------------
#ifdef AVERAGES
#ifdef BEC2_DIAG
      subroutine wrt_avg_bgc_diag (ncid,record)  ![
      ! write BEC2_DIAG to avg file
      implicit none

      ! inputs
      integer,intent(in) :: ncid
      integer,intent(in) :: record
      ! local
      integer :: itrc, ierr

      ! 2d diagnostics:
      do itrc=1,nr_bec2_diag_2d
         if (wrtavg_bec2_diag_2d(itrc)) then
           call nc_write_var(ncid, bec2_diag_2d_avg(START_2D_ARRAY,itrc),
     &          1, vname_bec2_diag_2d(1,itrc), rp_var, record, ierr )
         endif
      end do

      ! 3d diagnostics:
      do itrc=1,nr_bec2_diag_3d
         if (wrtavg_bec2_diag_3d(itrc)) then
           call nc_write_var(ncid, bec2_diag_3d_avg(START_2D_ARRAY,1,itrc),
     &          N, vname_bec2_diag_3d(1,itrc), rp_var, record, ierr )
         endif
      end do

      end subroutine wrt_avg_bgc_diag  !]
#endif /* BEC2_DIAG */
#endif /* AVERAGES */

! ----------------------------------------------------------------------
#ifdef AVERAGES
#ifdef BEC2_DIAG
      subroutine set_avg_bgc_diag(istrR,iendR,jstrR,jendR,cff,cff1)  ![

      implicit none

      ! inputs
      integer istrR,iendR,jstrR,jendR
      real cff,cff1
      ! local
      integer i,j,k,itrc

      do itrc=1,nr_bec2_diag_2d
         do j=jstrR,jendR
            do i=istrR,iendR
               bec2_diag_2d_avg(i,j,itrc) = cff1*bec2_diag_2d_avg(i,j,itrc)
     &                                    + cff*bec2_diag_2d(i,j,itrc)
            end do
         end do
      end do
      do itrc=1,nr_bec2_diag_3d
        do k=1,N
          do j=jstrR,jendR
             do i=istrR,iendR
                bec2_diag_3d_avg(i,j,k,itrc) = cff1*bec2_diag_3d_avg(i,j,k,itrc)
     &                                       + cff*bec2_diag_3d(i,j,k,itrc)
             end do
          end do
        end do
      end do

      end subroutine set_avg_bgc_diag  !]
#endif /* BEC2_DIAG */
#endif /* AVERAGES */

! ----------------------------------------------------------------------
      subroutine init_arrays_bgc_frc  ![
      implicit none

      allocate( nc_dust%data( GLOBAL_2D_ARRAY,2) )
      allocate( nc_iron%data( GLOBAL_2D_ARRAY,2) )

#ifdef PCO2AIR_FORCING
      allocate( nc_pco2air%data( GLOBAL_2D_ARRAY,2) )
#endif

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
      allocate( nc_swrad_avg%data( GLOBAL_2D_ARRAY,2) )
#endif

      end subroutine  !]

! ----------------------------------------------------------------------
      subroutine init_scalars_bec2()  ![
!
! Set initial values for  globally accessible (stored in common
! blocks) scalar variables of the BEC model.
!
      use comm_vars
      use scalars

      implicit none
      integer i,j, itrc ,lvar,lenstr ! ierr, DevinD removed ierr as arguement
      integer omp_get_num_threads
!#include "param.h"
!#include "ncvars"
!#include "ncroms_bgcflux.h"

! CN: take this part out of the NEVER_DEF part
!#ifdef BEC2_DIAG
!# include "ecosys_bec2.h"
!#endif

! DevinD commented out 20201215
!#ifdef NEVER_DEF
!
!!
!! Find out how many threads are created and check whether the number
!! of threads exceeds maximum allowed, as well as check if the number
!! of tiles is divisible by the number of threads, so the job can be
!! evenly partitioned, complain and signal to terminate, if something
!! is wrong.
!!
!C$OMP PARALLEL
!C$OMP CRITICAL (isca_cr_rgn)
!      numthreads=1
!C$      numthreads=omp_get_num_threads()
!C$OMP END CRITICAL (isca_cr_rgn)
!C$OMP END PARALLEL
!      MPI_master_only write(*,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
!     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'
!
!#ifdef MPI
!# define NUMPROC NNODES
!#else
!# define NUMPROC 1
!#endif
!
!      if (numthreads.gt.NUMPROC) then
!        MPI_master_only write(*,'(/1x,A,I3/)')
!     &    'ERROR: Requested number of threads exceeds setting: ',NUMPROC
!        ierr=ierr+1                                      !--> ERROR
!      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
!        MPI_master_only write(*,
!     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
!     &                'wrong choice of numthreads =', numthreads,
!     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
!        ierr=ierr+1                                      !--> ERROR
!      endif
!!
!      do i=0,2*NT+1
!        global_sum(i)=QuadZero
!      enddo
!#endif /* NEVER_DEF */

!
! Define names of variables in NetCDF output files.
! Inner dimension is for variable type and outer is for variable
! attributes (name, long-name, units, field).
!

      ! DevinD: NOTE: 'vname' was previously used to store all
      ! variable names (primary vars, tracers, etc).
      ! However, tracers variable names, etc are now stored in tracers.F
      ! in t_vname, t_lname, etc. Hence populating vname for bgc should not
      ! be necessary.
      ! Further, ETH vname(4,NT) but roms vname(3,NT), see below for PO4 where
      ! row '4' is populated with '  '. This overide data as only 3 rows!
      ! Hence I have removed.
      ! HOWEVER, vname_bec2_diag_2d & vname_bec2_diag_3d are STILL used and must
      ! be correctly populated.

#ifdef SOLVE3D
!      vname(1,indxPO4)='PO4                                       '
!      vname(2,indxpO4)='Phosphate                                 '
!      vname(3,indxpO4)='mMol P m-3                                '
!      vname(4,indxpO4)='  '

! DevinD left this in for now and not sure what it does but removed offending
! entry to row '4' of vname.
# ifndef BEC2_DIAG
      vname(1,indxPH)='PH'
      vname(2,indxPH)='surf pH value'
      vname(3,indxPH)=' '
!      vname(4,indxPH)=' '
      vname(1,indxPCO2)='pCO2                                 '
      vname(2,indxPCO2)='Surface water pCO2                   '
      vname(3,indxPCO2)='ppm                                  '
!      vname(4,indxPCO2)='  '
      vname(1,indxPCO2air)='pCO2air                           '
      vname(2,indxPCO2air)='Atmospheric pCO2                  '
      vname(3,indxPCO2air)='ppm                               '
!      vname(4,indxPCO2air)='           '
      vname(1,indxPARinc)='PARinc                                '
      vname(2,indxPARinc)=
     &  'Incoming photosynth. active radiation '
      vname(3,indxPARinc)='W m-2                                 '
!      vname(4,indxPARinc)='  '
      vname(1,indxPAR)='PAR'
      vname(2,indxPAR)=
     &  'Photosynthetically active radiation'
      vname(3,indxPAR)='W m-2'
!      vname(4,indxPAR)='  '
# endif /*BEC2_DIAG */



# ifdef BEC2_DIAG
      !
      ! 2d diagnotics:
      !
      vname_bec2_diag_2d(1,indxPH)='PH'
      vname_bec2_diag_2d(2,indxPH)='surf pH value'
      vname_bec2_diag_2d(3,indxPH)=' '
      vname_bec2_diag_2d(4,indxPH)=' '
      vname_bec2_diag_2d(1,indxPCO2)='pCO2'
      vname_bec2_diag_2d(2,indxPCO2)='Surface water pCO2'
      vname_bec2_diag_2d(3,indxPCO2)='ppm'
      vname_bec2_diag_2d(4,indxPCO2)='  '
      vname_bec2_diag_2d(1,indxPCO2air)='pCO2air'
      vname_bec2_diag_2d(2,indxPCO2air)='Atmospheric pCO2'
      vname_bec2_diag_2d(3,indxPCO2air)='ppm'
      vname_bec2_diag_2d(4,indxPCO2air)='           '
      vname_bec2_diag_2d(1,indxPARinc)='PARinc'
      vname_bec2_diag_2d(2,indxPARinc)=
     &  'Incoming photosynth. active radiation '
      vname_bec2_diag_2d(3,indxPARinc)='W m-2'
      vname_bec2_diag_2d(4,indxPARinc)='  '
      vname_bec2_diag_2d(1,indxFGO2)='FG_O2'
      vname_bec2_diag_2d(2,indxFGO2)='Air-sea flux of O2'
      vname_bec2_diag_2d(3,indxFGO2)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxFGO2)=' '
      vname_bec2_diag_2d(1,indxFGCO2)='FG_CO2'
      vname_bec2_diag_2d(2,indxFGCO2)='Air-sea flux of CO2'
      vname_bec2_diag_2d(3,indxFGCO2)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxFGCO2)=' '
      vname_bec2_diag_2d(1,indxWS10m)='WS10m'
      vname_bec2_diag_2d(2,indxWS10m)='Wind speed at 10m'
      vname_bec2_diag_2d(3,indxWS10m)='m/s '
      vname_bec2_diag_2d(4,indxWS10m)=' '
      vname_bec2_diag_2d(1,indxXKW)='XKW'
      vname_bec2_diag_2d(2,indxXKW)='XKW'
      vname_bec2_diag_2d(3,indxXKW)='m/s '
      vname_bec2_diag_2d(4,indxXKW)=' '
      vname_bec2_diag_2d(1,indxATMPRESS)='ATM_PRESS'
      vname_bec2_diag_2d(2,indxATMPRESS)='Atmospheric pressure'
      vname_bec2_diag_2d(3,indxATMPRESS)='atm '
      vname_bec2_diag_2d(4,indxATMPRESS)=' '

      vname_bec2_diag_2d(1,indxSCHMIDTO2)='SCHMIDT_O2'
      vname_bec2_diag_2d(2,indxSCHMIDTO2)='Schmidt number for O2'
      vname_bec2_diag_2d(3,indxSCHMIDTO2)=' '
      vname_bec2_diag_2d(4,indxSCHMIDTO2)=' '

      vname_bec2_diag_2d(1,indxO2SAT)='O2SAT'
      vname_bec2_diag_2d(2,indxO2SAT)='O2 saturation concentration'
      vname_bec2_diag_2d(3,indxO2SAT)='mmol/m3 '
      vname_bec2_diag_2d(4,indxO2SAT)=' '
      vname_bec2_diag_2d(1,indxSCHMIDTCO2)='SCHMIDT_CO2'
      vname_bec2_diag_2d(2,indxSCHMIDTCO2)='Schmidt number for CO2'
      vname_bec2_diag_2d(3,indxSCHMIDTCO2)=' '
      vname_bec2_diag_2d(4,indxSCHMIDTCO2)=' '

      vname_bec2_diag_2d(1,indxPVO2)='PV_O2'
      vname_bec2_diag_2d(2,indxPVO2)='Piston velocity for O2'
      vname_bec2_diag_2d(3,indxPVO2)='m/s '
      vname_bec2_diag_2d(4,indxPVO2)=' '
      vname_bec2_diag_2d(1,indxPVCO2)='PV_CO2'
      vname_bec2_diag_2d(2,indxPVCO2)='Piston velocity for CO2'
      vname_bec2_diag_2d(3,indxPVCO2)='m/s '
      vname_bec2_diag_2d(4,indxPVCO2)=' '
# ifdef Ncycle_SY
      vname_bec2_diag_2d(1,indxfgn2o)='FG_N2O'
      vname_bec2_diag_2d(2,indxfgn2o)=
     &   'Air-sea flux of N2O '
      vname_bec2_diag_2d(3,indxfgn2o)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2o)=' '

      vname_bec2_diag_2d(1,indxschmidt_n2o)='SCHMIDT_N2O'
      vname_bec2_diag_2d(2,indxschmidt_n2o)='Schmidt number for N2O'
      vname_bec2_diag_2d(3,indxschmidt_n2o)=' '
      vname_bec2_diag_2d(4,indxschmidt_n2o)=' '

      vname_bec2_diag_2d(1,indxn2osat)='N2OSAT'
      vname_bec2_diag_2d(2,indxn2osat)='N2O saturation concentration'
      vname_bec2_diag_2d(3,indxn2osat)='mmol/m3 '
      vname_bec2_diag_2d(4,indxn2osat)=' '

      vname_bec2_diag_2d(1,indxpvn2o)='PV_N2O'
      vname_bec2_diag_2d(2,indxpvn2o)='Piston velocity for N2O'
      vname_bec2_diag_2d(3,indxpvn2o)='m/s '
      vname_bec2_diag_2d(4,indxpvn2o)=' '

      vname_bec2_diag_2d(1,indxfgn2)='FG_N2'
      vname_bec2_diag_2d(2,indxfgn2)=
     &   'Air-sea flux of N2 '
      vname_bec2_diag_2d(3,indxfgn2)='mmol/m2/s '
      vname_bec2_diag_2d(4,indxfgn2)=' '

      vname_bec2_diag_2d(1,indxschmidt_n2)='SCHMIDT_N2'
      vname_bec2_diag_2d(2,indxschmidt_n2)='Schmidt number for N2'
      vname_bec2_diag_2d(3,indxschmidt_n2)=' '
      vname_bec2_diag_2d(4,indxschmidt_n2)=' '

      vname_bec2_diag_2d(1,indxn2sat)='N2SAT'
      vname_bec2_diag_2d(2,indxn2sat)='N2 saturation concentration'
      vname_bec2_diag_2d(3,indxn2sat)='mmol/m3 '
      vname_bec2_diag_2d(4,indxn2sat)=' '

      vname_bec2_diag_2d(1,indxpvn2)='PV_N2'
      vname_bec2_diag_2d(2,indxpvn2)='Piston velocity for N2'
      vname_bec2_diag_2d(3,indxpvn2)='m/s '
      vname_bec2_diag_2d(4,indxpvn2)=' '
#endif
      vname_bec2_diag_2d(1,indxCO2STAR)='CO2STAR'
      vname_bec2_diag_2d(2,indxCO2STAR)='CO2STAR'
      vname_bec2_diag_2d(3,indxCO2STAR)='mmol/m3 '
      vname_bec2_diag_2d(4,indxCO2STAR)=' '

      vname_bec2_diag_2d(1,indxDCO2STAR)='DCO2STAR'
      vname_bec2_diag_2d(2,indxDCO2STAR)='DCO2STAR'
      vname_bec2_diag_2d(3,indxDCO2STAR)='mmol/m3'
      vname_bec2_diag_2d(4,indxDCO2STAR)=' '

      vname_bec2_diag_2d(1,indxPCO2OC)='PCO2OC'
      vname_bec2_diag_2d(2,indxPCO2OC)='PCO2OC'
      vname_bec2_diag_2d(3,indxPCO2OC)='not looked up yet '
      vname_bec2_diag_2d(4,indxPCO2OC)=' '

      vname_bec2_diag_2d(1,indxIRONFLUX)='IRON_FLUX'
      vname_bec2_diag_2d(2,indxIRONFLUX)='Iron surface flux'
      vname_bec2_diag_2d(3,indxIRONFLUX)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxIRONFLUX)=' '

      vname_bec2_diag_2d(1,indxSEDDENITRIF)='SED_DENITRIF'
      vname_bec2_diag_2d(2,indxSEDDENITRIF)='Sediment denitrification'
      vname_bec2_diag_2d(3,indxSEDDENITRIF)='mmol/m3/s'
      vname_bec2_diag_2d(4,indxSEDDENITRIF)=' '

      vname_bec2_diag_2d(1,indxFESEDFLUX)='SED_FE_FLUX'
      vname_bec2_diag_2d(2,indxFESEDFLUX)='Benthic Iron Flux'
      vname_bec2_diag_2d(3,indxFESEDFLUX)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxFESEDFLUX)=' '

      vname_bec2_diag_2d(1,indxFLUXTOSED)='FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxFLUXTOSED)='POC reaching the sediments'
      vname_bec2_diag_2d(3,indxFLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxFLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxCACO3FLUXTOSED)='CACO3_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxCACO3FLUXTOSED)='CaCO3 reaching the sed'
      vname_bec2_diag_2d(3,indxCACO3FLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxCACO3FLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxSIO2FLUXTOSED)='SIO2_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxSIO2FLUXTOSED)='Opal reaching the sed'
      vname_bec2_diag_2d(3,indxSIO2FLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxSIO2FLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxPIRONFLUXTOSED)='PIRON_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxPIRONFLUXTOSED)='P iron reaching the sed'
      vname_bec2_diag_2d(3,indxPIRONFLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxPIRONFLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxDUSTFLUXTOSED)='DUST_FLUX_TO_SED'
      vname_bec2_diag_2d(2,indxDUSTFLUXTOSED)='Dust reaching the sed'
      vname_bec2_diag_2d(3,indxDUSTFLUXTOSED)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxDUSTFLUXTOSED)=' '

      vname_bec2_diag_2d(1,indxPOCSEDLOSS)='POC_SED_LOSS'
      vname_bec2_diag_2d(2,indxPOCSEDLOSS)='POC lost to sed at bottom'
      vname_bec2_diag_2d(3,indxPOCSEDLOSS)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxPOCSEDLOSS)=' '

      vname_bec2_diag_2d(1,indxOTHERREMIN)='OTHER_REMIN'
      vname_bec2_diag_2d(2,indxOTHERREMIN)='other remin pathways in sed'
      vname_bec2_diag_2d(3,indxOTHERREMIN)='mmol/m3/s'
      vname_bec2_diag_2d(4,indxOTHERREMIN)=' '

      vname_bec2_diag_2d(1,indxCACO3SEDLOSS)='P_CACO3_SED_LOSS'
      vname_bec2_diag_2d(2,indxCACO3SEDLOSS)='CaCO3 lost to sed'
      vname_bec2_diag_2d(3,indxCACO3SEDLOSS)='mmol/m2/s'
      vname_bec2_diag_2d(4,indxCACO3SEDLOSS)=' '

      vname_bec2_diag_2d(1,indxSIO2SEDLOSS)='SIO2_SED_LOSS'
      vname_bec2_diag_2d(2,indxSIO2SEDLOSS)='SiO2 lost to sed'
      vname_bec2_diag_2d(3,indxSIO2SEDLOSS)='mmol Si/m2/s'
      vname_bec2_diag_2d(4,indxSIO2SEDLOSS)=' '


      ! MF: Control BEC2_DIAG Output vars
      wrtavg_bec2_diag_2d=(/.false.,.false.,.false.,.true.,.true.,   ! pCO2air, PARinc,FG_O2,FG_CO2, WS10m
     &                     .false.,.false.,.false.,.false.,.false.,   ! XKW, ATM_PRESS,SCHMIDT_O2,O2SAT, SCHMIDT_CO2
     &                     .false.,.false.,.false.,.false.,.false.,    ! PVO2, PV_CO2,IRON_FLUX, SED_DENITRI, PH
     &                     .false.,.false.,.false.,.false.,          ! PCO2, CO2STAR, PCO2OC, DCO2STAR
#ifdef Ncycle_SY
     &                     .false.,.false.,.false.,.false., !SCHMIDT_N2O, PV_N2O, N2OSAT, FG_N2O
     &                     .false.,.false.,.false.,.false., !SCHMIDT_N2,PV_N2, N2SAT, FG_N2
# endif
     &                     .true.,.false.,.false.,.false.,.false.,   ! FE_SED_FLUX,FLUX_TO_SED,CACO3_FLUX_TO_SED,SIO2_FLUX_TO_SED,PIRON_FLUX_TO_SED
     &                     .false.,.false.,.false.,.false.,.false./)  ! DUST_FLUX_SED, POC_SED_LOSS, OTHER_REMIN,CaCO3_sed_loss,P_SiO2_sed_loss

      !
      ! More 3d diagnostics:
      !
      vname_bec2_diag_3d(1,indxPAR)='PAR'
      vname_bec2_diag_3d(2,indxPAR)=
     &  'Photosynthetically active radiation'
      vname_bec2_diag_3d(3,indxPAR)='W m-2'
      vname_bec2_diag_3d(4,indxPAR)='  '

      vname_bec2_diag_3d(1,indxPOCFLUXIN)='POC_FLUX_IN'
      vname_bec2_diag_3d(2,indxPOCFLUXIN)='POC flux into cell'
      vname_bec2_diag_3d(3,indxPOCFLUXIN)='mmol C/m2/s '
      vname_bec2_diag_3d(4,indxPOCFLUXIN)=' '

      vname_bec2_diag_3d(1,indxPOCPROD)='POC_PROD'
      vname_bec2_diag_3d(2,indxPOCPROD)='POC production'
      vname_bec2_diag_3d(3,indxPOCPROD)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPOCPROD)=' '
      vname_bec2_diag_3d(1,indxPOCREMIN)='POC_REMIN'
      vname_bec2_diag_3d(2,indxPOCREMIN)='POC remineralization'
      vname_bec2_diag_3d(3,indxPOCREMIN)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxPOCREMIN)=' '

      vname_bec2_diag_3d(1,indxCACO3FLUXIN)='CACO3_FLUX_IN'
      vname_bec2_diag_3d(2,indxCACO3FLUXIN)=
     &  'Incoming flux of large detrital CaCO3'
      vname_bec2_diag_3d(3,indxCACO3FLUXIN)='mmol CaCO3/m2/s '
      vname_bec2_diag_3d(4,indxCACO3FLUXIN)=' '
      vname_bec2_diag_3d(1,indxPCACO3PROD)='PCACO3_PROD'
      vname_bec2_diag_3d(2,indxPCACO3PROD)=
     &  'Production of large detrital CaCO3'
      vname_bec2_diag_3d(3,indxPCACO3PROD)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxPCACO3PROD)=' '
      vname_bec2_diag_3d(1,indxCACO3REMIN)='CACO3_REMIN'
      vname_bec2_diag_3d(2,indxCACO3REMIN)=
     &  'Remineralization of large detrital CaCO3'
      vname_bec2_diag_3d(3,indxCACO3REMIN)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxCACO3REMIN)=' '

      vname_bec2_diag_3d(1,indxSIO2FLUXIN)='SIO2_FLUX_IN'
      vname_bec2_diag_3d(2,indxSIO2FLUXIN)=
     &  'Incoming flux of large detritus SiO2'
      vname_bec2_diag_3d(3,indxSIO2FLUXIN)='mmol Si/m2/s '
      vname_bec2_diag_3d(4,indxSIO2FLUXIN)=' '
      vname_bec2_diag_3d(1,indxSIO2PROD)='SIO2_PROD'
      vname_bec2_diag_3d(2,indxSIO2PROD)=
     &  'Production of large detritus SiO2'
      vname_bec2_diag_3d(3,indxSIO2PROD)='mmol Si/m3/s '
      vname_bec2_diag_3d(4,indxSIO2PROD)=' '
      vname_bec2_diag_3d(1,indxSIO2REMIN)='SIO2_REMIN'
      vname_bec2_diag_3d(2,indxSIO2REMIN)=
     &  'Remineralization of large detritus SiO2'
      vname_bec2_diag_3d(3,indxSIO2REMIN)='mmol Si/m3/s'
      vname_bec2_diag_3d(4,indxSIO2REMIN)=' '

      vname_bec2_diag_3d(1,indxDUSTFLUXIN)='DUST_FLUX_IN'
      vname_bec2_diag_3d(2,indxDUSTFLUXIN)=
     &  'Incoming dust flux'
      vname_bec2_diag_3d(3,indxDUSTFLUXIN)='kg dust/m2/s '
      vname_bec2_diag_3d(4,indxDUSTFLUXIN)=' '
      vname_bec2_diag_3d(1,indxDUSTREMIN)='DUST_REMIN'
      vname_bec2_diag_3d(2,indxDUSTREMIN)=
     &  'Remineralization of dust'
      vname_bec2_diag_3d(3,indxDUSTREMIN)='kg dust/m3/s '
      vname_bec2_diag_3d(4,indxDUSTREMIN)=' '

      vname_bec2_diag_3d(1,indxPIRONFLUXIN)='P_IRON_FLUX_IN'
      vname_bec2_diag_3d(2,indxPIRONFLUXIN)=
     &  'Incoming flux of large detritus iron'
      vname_bec2_diag_3d(3,indxPIRONFLUXIN)='mmol Fe/m2/s '
      vname_bec2_diag_3d(4,indxPIRONFLUXIN)=' '
      vname_bec2_diag_3d(1,indxPIRONPROD)='P_IRON_PROD'
      vname_bec2_diag_3d(2,indxPIRONPROD)=
     &  'Production of large detritus iron'
      vname_bec2_diag_3d(3,indxPIRONPROD)='mmol Fe/m3/s '
      vname_bec2_diag_3d(4,indxPIRONPROD)=' '
      vname_bec2_diag_3d(1,indxPIRONREMIN)='P_IRON_REMIN'
      vname_bec2_diag_3d(2,indxPIRONREMIN)=
     &  'Remineralization of large detritus iron'
      vname_bec2_diag_3d(3,indxPIRONREMIN)='mmol Fe/m3/s '
      vname_bec2_diag_3d(4,indxPIRONREMIN)=' '

      vname_bec2_diag_3d(1,indxGRAZESP)='GRAZE_SP'
      vname_bec2_diag_3d(2,indxGRAZESP)=
     &  'Grazing rate on small phytoplankton'
      vname_bec2_diag_3d(3,indxGRAZESP)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxGRAZESP)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIAT)='GRAZE_DIAT'
      vname_bec2_diag_3d(2,indxGRAZEDIAT)=
     &  'Grazing rate on diatoms'
      vname_bec2_diag_3d(3,indxGRAZEDIAT)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIAT)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIAZ)='GRAZE_DIAZ'
      vname_bec2_diag_3d(2,indxGRAZEDIAZ)=
     &  'Grazing rate on diazotrophs'
      vname_bec2_diag_3d(3,indxGRAZEDIAZ)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIAZ)=' '

      vname_bec2_diag_3d(1,indxSPLOSS)='SP_LOSS'
      vname_bec2_diag_3d(2,indxSPLOSS)=
     &  'Small phytoplankton non-grazing mortality'
      vname_bec2_diag_3d(3,indxSPLOSS)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxSPLOSS)=' '
      vname_bec2_diag_3d(1,indxDIATLOSS)='DIAT_LOSS'
      vname_bec2_diag_3d(2,indxDIATLOSS)=
     &  'Diatom non-grazing mortality'
      vname_bec2_diag_3d(3,indxDIATLOSS)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxDIATLOSS)=' '
      vname_bec2_diag_3d(1,indxZOOLOSS)='ZOO_LOSS'
      vname_bec2_diag_3d(2,indxZOOLOSS)=
     &  'Mortality due to higher trophic grazing on zooplankton'
      vname_bec2_diag_3d(3,indxZOOLOSS)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxZOOLOSS)=' '

      vname_bec2_diag_3d(1,indxSPAGG)='SP_AGG'
      vname_bec2_diag_3d(2,indxSPAGG)=
     &  'Aggregation of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxSPAGG)=' '
      vname_bec2_diag_3d(1,indxDIATAGG)='DIAT_AGG'
      vname_bec2_diag_3d(2,indxDIATAGG)=
     &  'Aggregation of diatoms'
      vname_bec2_diag_3d(3,indxDIATAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxDIATAGG)=' '

      vname_bec2_diag_3d(1,indxPHOTOCSP)='PHOTOC_SP'
      vname_bec2_diag_3d(2,indxPHOTOCSP)=
     &  'C fixation rate by small phytoplankton'
      vname_bec2_diag_3d(3,indxPHOTOCSP)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPHOTOCSP)=' '
      vname_bec2_diag_3d(1,indxPHOTOCDIAT)='PHOTOC_DIAT'
      vname_bec2_diag_3d(2,indxPHOTOCDIAT)=
     &  'C fixation rate by diatoms'
      vname_bec2_diag_3d(3,indxPHOTOCDIAT)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPHOTOCDIAT)=' '

      vname_bec2_diag_3d(1,indxTOTPROD)='TOT_PROD'
      vname_bec2_diag_3d(2,indxTOTPROD)=
     &  'Total autotroph production'
      vname_bec2_diag_3d(3,indxTOTPROD)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxTOTPROD)=' '
      vname_bec2_diag_3d(1,indxDOCPROD)='DOC_PROD'
      vname_bec2_diag_3d(2,indxDOCPROD)='DOC production'
      vname_bec2_diag_3d(3,indxDOCPROD)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxDOCPROD)=' '

      vname_bec2_diag_3d(1,indxDOCREMIN)='DOC_REMIN'
      vname_bec2_diag_3d(2,indxDOCREMIN)=
     &  'Remineralization rate of DOC'
      vname_bec2_diag_3d(3,indxDOCREMIN)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxDOCREMIN)=' '
      vname_bec2_diag_3d(1,indxFESCAVENGE)='FE_SCAVENGE'
      vname_bec2_diag_3d(2,indxFESCAVENGE)=
     &  'Loss of dissolved iron, scavenging'
      vname_bec2_diag_3d(3,indxFESCAVENGE)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxFESCAVENGE)=' '

      vname_bec2_diag_3d(1,indxSPNLIM)='SP_N_LIM'
      vname_bec2_diag_3d(2,indxSPNLIM)=
     &  'N limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPNLIM)=' '
      vname_bec2_diag_3d(4,indxSPNLIM)=' '
      vname_bec2_diag_3d(1,indxSPFEUPTAKE)='SP_FE_LIM'
      vname_bec2_diag_3d(2,indxSPFEUPTAKE)=
     &  'Fe limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxSPFEUPTAKE)=' '
      vname_bec2_diag_3d(1,indxSPPO4UPTAKE)='SP_PO4_LIM'
      vname_bec2_diag_3d(2,indxSPPO4UPTAKE)=
     &  'PO4 limitation factor of small phytoplankton'
      vname_bec2_diag_3d(3,indxSPPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxSPPO4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxSPLIGHTLIM)='SP_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxSPLIGHTLIM)=
     &  'Small phytoplankton light limitation factor'
      vname_bec2_diag_3d(3,indxSPLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxSPLIGHTLIM)=' '

      vname_bec2_diag_3d(1,indxDIATNLIM)='DIAT_N_LIM'
      vname_bec2_diag_3d(2,indxDIATNLIM)=
     &  'N limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATNLIM)=' '
      vname_bec2_diag_3d(4,indxDIATNLIM)=' '
      vname_bec2_diag_3d(1,indxDIATFEUPTAKE)='DIAT_FE_LIM'
      vname_bec2_diag_3d(2,indxDIATFEUPTAKE)=
     &  'Fe limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIATFEUPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATPO4UPTAKE)='DIAT_PO4_LIM'
      vname_bec2_diag_3d(2,indxDIATPO4UPTAKE)=
     &  'PO4 limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIATPO4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATSIO3UPTAKE)='DIAT_SIO3_LIM'
      vname_bec2_diag_3d(2,indxDIATSIO3UPTAKE)=
     &  'SiO3 limitation factor of diatoms'
      vname_bec2_diag_3d(3,indxDIATSIO3UPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIATSIO3UPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIATLIGHTLIM)='DIAT_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxDIATLIGHTLIM)=
     &  'Diatom light limitation factor'
      vname_bec2_diag_3d(3,indxDIATLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxDIATLIGHTLIM)=' '

      vname_bec2_diag_3d(1,indxCACO3PROD)='CACO3_PROD'
      vname_bec2_diag_3d(2,indxCACO3PROD)=
     &  'Production of CaCO3 by small phytoplankton'
      vname_bec2_diag_3d(3,indxCACO3PROD)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxCACO3PROD)=' '
      vname_bec2_diag_3d(1,indxDIAZNFIX)='DIAZ_NFIX'
      vname_bec2_diag_3d(2,indxDIAZNFIX)=
     &  'Total nitrogen fixation by diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNFIX)='mmol/m3/s   '
      vname_bec2_diag_3d(4,indxDIAZNFIX)=' '

      vname_bec2_diag_3d(1,indxDIAZLOSS)='DIAZ_LOSS'
      vname_bec2_diag_3d(2,indxDIAZLOSS)=
     &  'Diazotroph non-grazing mortality'
      vname_bec2_diag_3d(3,indxDIAZLOSS)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZLOSS)=' '
      vname_bec2_diag_3d(1,indxPHOTOCDIAZ)='PHOTOC_DIAZ'
      vname_bec2_diag_3d(2,indxPHOTOCDIAZ)=
     &   'Diazotroph C-fixation'
      vname_bec2_diag_3d(3,indxPHOTOCDIAZ)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxPHOTOCDIAZ)=' '

      vname_bec2_diag_3d(1,indxDIAZPO4UPTAKE)='DIAZ_PO4_LIM'
      vname_bec2_diag_3d(2,indxDIAZPO4UPTAKE)=
     &  'PO4 limitation factor of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIAZPO4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIAZFEUPTAKE)='DIAZ_FE_LIM'
      vname_bec2_diag_3d(2,indxDIAZFEUPTAKE)=
     &  'Fe limitation factor of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxDIAZFEUPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIAZLIGHTLIM)='DIAZ_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxDIAZLIGHTLIM)=
     &   'Diazotroph light limitation'
      vname_bec2_diag_3d(3,indxDIAZLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxDIAZLIGHTLIM)=' '
      vname_bec2_diag_3d(1,indxFESCAVENGERATE)='FE_SCAVENGE_RATE'
      vname_bec2_diag_3d(2,indxFESCAVENGERATE)=
     &  'Annual scavenging rate of iron as % of ambient'
      vname_bec2_diag_3d(3,indxFESCAVENGERATE)=' '
      vname_bec2_diag_3d(4,indxFESCAVENGERATE)=' '

      vname_bec2_diag_3d(1,indxDONPROD)='DON_PROD'
      vname_bec2_diag_3d(2,indxDONPROD)=
     &  'Production of dissolved organic N'
      vname_bec2_diag_3d(3,indxDONPROD)=' '
      vname_bec2_diag_3d(4,indxDONPROD)=' '
      vname_bec2_diag_3d(1,indxDONREMIN)='DON_REMIN'
      vname_bec2_diag_3d(2,indxDONREMIN)=
     &  'Remineralization rate of DON'
      vname_bec2_diag_3d(3,indxDONREMIN)='mmol N/m3/s  '
      vname_bec2_diag_3d(4,indxDONREMIN)=' '

      vname_bec2_diag_3d(1,indxDOFEPROD)='DOFE_PROD'
      vname_bec2_diag_3d(2,indxDOFEPROD)=
     &  'Production of dissolved organic iron'
      vname_bec2_diag_3d(3,indxDOFEPROD)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDOFEPROD)=' '
      vname_bec2_diag_3d(1,indxDOFEREMIN)='DOFE_REMIN'
      vname_bec2_diag_3d(2,indxDOFEREMIN)=
     &  'Remineralization rate of DOFE'
      vname_bec2_diag_3d(3,indxDOFEREMIN)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDOFEREMIN)=' '

      vname_bec2_diag_3d(1,indxDOPPROD)='DOP_PROD'
      vname_bec2_diag_3d(2,indxDOPPROD)=
     &  'Production of dissolved organic P'
      vname_bec2_diag_3d(3,indxDOPPROD)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDOPPROD)=' '
      vname_bec2_diag_3d(1,indxDOPREMIN)='DOP_REMIN'
      vname_bec2_diag_3d(2,indxDOPREMIN)=
     &  'Remineralization rate of DOP'
      vname_bec2_diag_3d(3,indxDOPREMIN)='mmol/m3/s '
      vname_bec2_diag_3d(4,indxDOPREMIN)=' '

      vname_bec2_diag_3d(1,indxDIATSIUPTAKE)='DIAT_SI_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATSIUPTAKE)=
     &  'Silicon uptake rate by diatoms'
      vname_bec2_diag_3d(3,indxDIATSIUPTAKE)='mmol Si/m3/s  '
      vname_bec2_diag_3d(4,indxDIATSIUPTAKE)=' '

      vname_bec2_diag_3d(1,indxIRONUPTAKESP)='PHOTOFE_SP'
      vname_bec2_diag_3d(2,indxIRONUPTAKESP)=
     &  'Iron uptake rate by small phytoplankton'
      vname_bec2_diag_3d(3,indxIRONUPTAKESP)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKESP)=' '
      vname_bec2_diag_3d(1,indxIRONUPTAKEDIAT)='PHOTOFE_DIAT'
      vname_bec2_diag_3d(2,indxIRONUPTAKEDIAT)=
     &  'Iron uptake rate by diatoms'
      vname_bec2_diag_3d(3,indxIRONUPTAKEDIAT)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKEDIAT)=' '
      vname_bec2_diag_3d(1,indxIRONUPTAKEDIAZ)='PHOTOFE_DIAZ'
      vname_bec2_diag_3d(2,indxIRONUPTAKEDIAZ)=
     &  'Iron uptake rate by diatotrophs'
      vname_bec2_diag_3d(3,indxIRONUPTAKEDIAZ)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKEDIAZ)=' '

      vname_bec2_diag_3d(1,indxNITRIF)='NITRIF'
      vname_bec2_diag_3d(2,indxNITRIF)='Nitrification'
      vname_bec2_diag_3d(3,indxNITRIF)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxNITRIF)=' '
      vname_bec2_diag_3d(1,indxDENITRIF)='DENITRIF'
      vname_bec2_diag_3d(2,indxDENITRIF)='Denitrification'
      vname_bec2_diag_3d(3,indxDENITRIF)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDENITRIF)=' '

      vname_bec2_diag_3d(1,indxSPNO3UPTAKE)='SP_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxSPNO3UPTAKE)=
     &   'NO3 uptake rate of small phyto'
      vname_bec2_diag_3d(3,indxSPNO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxSPNO3UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATNO3UPTAKE)='DIAT_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATNO3UPTAKE)=
     &   'NO3 uptake rate of diatoms'
      vname_bec2_diag_3d(3,indxDIATNO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIATNO3UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIAZNO3UPTAKE)='DIAZ_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxDIAZNO3UPTAKE)=
     &   'NO3 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZNO3UPTAKE)=' '

      vname_bec2_diag_3d(1,indxSPNH4UPTAKE)='SP_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxSPNH4UPTAKE)=
     &   'NH4 uptake rate of small phyto'
      vname_bec2_diag_3d(3,indxSPNH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxSPNH4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIATNH4UPTAKE)='DIAT_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATNH4UPTAKE)=
     &   'NH4 uptake rate of diatoms'
      vname_bec2_diag_3d(3,indxDIATNH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIATNH4UPTAKE)=' '
      vname_bec2_diag_3d(1,indxDIAZNH4UPTAKE)='DIAZ_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxDIAZNH4UPTAKE)=
     &   'NH4 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZNH4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxGRAZEDICSP)='SP_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxGRAZEDICSP)=
     &   'Small phyto grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxGRAZEDICSP)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxGRAZEDICSP)=' '
      vname_bec2_diag_3d(1,indxGRAZEDICDIAT)='DIAT_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxGRAZEDICDIAT)=
     &   'Diatom grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxGRAZEDICDIAT)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxGRAZEDICDIAT)=' '
      vname_bec2_diag_3d(1,indxGRAZEDICDIAZ)='DIAZ_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxGRAZEDICDIAZ)=
     &   'Diazotroph grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxGRAZEDICDIAZ)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxGRAZEDICDIAZ)=' '

      vname_bec2_diag_3d(1,indxLOSSDICSP)='SP_LOSS_DIC'
      vname_bec2_diag_3d(2,indxLOSSDICSP)=
     &   'Small phyto non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxLOSSDICSP)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxLOSSDICSP)=' '
      vname_bec2_diag_3d(1,indxLOSSDICDIAT)='DIAT_LOSS_DIC'
      vname_bec2_diag_3d(2,indxLOSSDICDIAT)=
     &   'Diatom non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxLOSSDICDIAT)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxLOSSDICDIAT)=' '
      vname_bec2_diag_3d(1,indxLOSSDICDIAZ)='DIAZ_LOSS_DIC'
      vname_bec2_diag_3d(2,indxLOSSDICDIAZ)=
     &   'Diazotroph non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxLOSSDICDIAZ)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxLOSSDICDIAZ)=' '

      vname_bec2_diag_3d(1,indxZOOLOSSDIC)='ZOO_LOSS_DIC'
      vname_bec2_diag_3d(2,indxZOOLOSSDIC)=
     &   'Zooplankton loss routed to DIC'
      vname_bec2_diag_3d(3,indxZOOLOSSDIC)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxZOOLOSSDIC)=' '

      vname_bec2_diag_3d(1,indxDIAZAGG)='DIAZ_AGG'
      vname_bec2_diag_3d(2,indxDIAZAGG)=
     &  'Aggregation of diatoms'
      vname_bec2_diag_3d(3,indxDIAZAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxDIAZAGG)=' '

      vname_bec2_diag_3d(1,indxGRAZESPZOO)='GRAZE_SP_ZOO'
      vname_bec2_diag_3d(2,indxGRAZESPZOO)=
     &  'Grazing rate on small phytoplankton routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZESPZOO)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxGRAZESPZOO)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIATZOO)='GRAZE_DIAT_ZOO'
      vname_bec2_diag_3d(2,indxGRAZEDIATZOO)=
     &  'Grazing rate on diatoms routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZEDIATZOO)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIATZOO)=' '
      vname_bec2_diag_3d(1,indxGRAZEDIAZZOO)='GRAZE_DIAZ_ZOO'
      vname_bec2_diag_3d(2,indxGRAZEDIAZZOO)=
     &  'Grazing rate on diazotrophs routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZEDIAZZOO)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZEDIAZZOO)=' '

      vname_bec2_diag_3d(1,indxSPQCACO3)='QCACO3'
      vname_bec2_diag_3d(2,indxSPQCACO3)=
     &  'Small phyto CaCO3/C ratio'
      vname_bec2_diag_3d(3,indxSPQCACO3)='mmol CaCO3/mmol C'
      vname_bec2_diag_3d(4,indxSPQCACO3)=' '

      vname_bec2_diag_3d(1,indxSPPHOTOACC)='PHOTOACC_SP'
      vname_bec2_diag_3d(2,indxSPPHOTOACC)=
     &  'Chl synthesis in photoadaptation for small phyto'
      vname_bec2_diag_3d(3,indxSPPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxSPPHOTOACC)=' '
      vname_bec2_diag_3d(1,indxDIATPHOTOACC)='PHOTOACC_DIAT'
      vname_bec2_diag_3d(2,indxDIATPHOTOACC)=
     &  'Chl synthesis in photoadaptation for diatoms'
      vname_bec2_diag_3d(3,indxDIATPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxDIATPHOTOACC)=' '
      vname_bec2_diag_3d(1,indxDIAZPHOTOACC)='PHOTOACC_DIAZ'
      vname_bec2_diag_3d(2,indxDIAZPHOTOACC)=
     &  'Chl synthesis in photoadaptation for diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxDIAZPHOTOACC)=' '

      vname_bec2_diag_3d(1,indxSPCZERO)='SPC_ZERO'
      vname_bec2_diag_3d(2,indxSPCZERO)=
     &  'Change caused by setting negative SPC values to zero'
      vname_bec2_diag_3d(3,indxSPCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxSPCZERO)=' '
      vname_bec2_diag_3d(1,indxDIATCZERO)='DIATC_ZERO'
      vname_bec2_diag_3d(2,indxDIATCZERO)=
     &  'Change caused by setting negative DIATC values to zero'
      vname_bec2_diag_3d(3,indxDIATCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxDIATCZERO)=' '
      vname_bec2_diag_3d(1,indxDIAZCZERO)='DIAZC_ZERO'
      vname_bec2_diag_3d(2,indxDIAZCZERO)=
     &  'Change caused by setting negative DIAZC values to zero'
      vname_bec2_diag_3d(3,indxDIAZCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxDIAZCZERO)=' '
      vname_bec2_diag_3d(1,indxDOCZERO)='DOC_ZERO'
      vname_bec2_diag_3d(2,indxDOCZERO)=
     &  'Change caused by setting negative DOC values to zero'
      vname_bec2_diag_3d(3,indxDOCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxDOCZERO)=' '
      vname_bec2_diag_3d(1,indxZOOCZERO)='ZOOC_ZERO'
      vname_bec2_diag_3d(2,indxZOOCZERO)=
     &  'Change caused by setting negative ZOOC values to zero'
      vname_bec2_diag_3d(3,indxZOOCZERO)='mMol C/m3/s'
      vname_bec2_diag_3d(4,indxZOOCZERO)=' '
      vname_bec2_diag_3d(1,indxSPCACO3ZERO)='SPCACO3_ZERO'
      vname_bec2_diag_3d(2,indxSPCACO3ZERO)=
     &  'Change caused by setting negative SPCACO3 values to zero'
      vname_bec2_diag_3d(3,indxSPCACO3ZERO)='mMol CaCO3/m3/s'
      vname_bec2_diag_3d(4,indxSPCACO3ZERO)=' '

      vname_bec2_diag_3d(1,indxDONRREMIN) = 'DONr_REMIN'
      vname_bec2_diag_3d(2,indxDONRREMIN)=
     &  'Portion of refractory DON remineralized'
      vname_bec2_diag_3d(3,indxDONRREMIN)='mMol N/m3/s'
      vname_bec2_diag_3d(4,indxDONRREMIN)=' '

      vname_bec2_diag_3d(1,indxTOTCHL)='TOT_CHL'
      vname_bec2_diag_3d(2,indxTOTCHL)=
     &  'Total Chlorophyll'
      vname_bec2_diag_3d(3,indxTOTCHL)='mg Chl/m3'
      vname_bec2_diag_3d(4,indxTOTCHL)=' '

      vname_bec2_diag_3d(1,indxTOTPHYTOC)='TOT_PHYTOC'
      vname_bec2_diag_3d(2,indxTOTPHYTOC)=
     &  'Total Phytoplankton Carbon'
      vname_bec2_diag_3d(3,indxTOTPHYTOC)='mMol C/m3'
      vname_bec2_diag_3d(4,indxTOTPHYTOC)=' '

      vname_bec2_diag_3d(1,indxSPPLIM)='SP_P_LIM'
      vname_bec2_diag_3d(2,indxSPPLIM)=
     &  'Small phytoplankton P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxSPPLIM)=' '
      vname_bec2_diag_3d(4,indxSPPLIM)=' '

      vname_bec2_diag_3d(1,indxDIATPLIM)='DIAT_P_LIM'
      vname_bec2_diag_3d(2,indxDIATPLIM)=
     &  'Diatom P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxDIATPLIM)=' '
      vname_bec2_diag_3d(4,indxDIATPLIM)=' '

      vname_bec2_diag_3d(1,indxDIAZPLIM)='DIAZ_P_LIM'
      vname_bec2_diag_3d(2,indxDIAZPLIM)=
     &  'Diazotroph P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxDIAZPLIM)=' '
      vname_bec2_diag_3d(4,indxDIAZPLIM)=' '

      vname_bec2_diag_3d(1,indxo2cons)='O2_CONSUMPTION'
      vname_bec2_diag_3d(2,indxo2cons)=
     &  'O2 consumption rate'
      vname_bec2_diag_3d(3,indxo2cons)='mmol O2/m3/s'
      vname_bec2_diag_3d(4,indxo2cons)=' '

      vname_bec2_diag_3d(1,indxo2prod)='O2_PRODUCTION'
      vname_bec2_diag_3d(2,indxo2prod)=
     &  'O2 production rate'
      vname_bec2_diag_3d(3,indxo2prod)='mmol O2/m3/s'
      vname_bec2_diag_3d(4,indxo2prod)=' '

# ifdef BEC_COCCO
      vname_bec2_diag_3d(1,indxGRAZECOCCO)='GRAZE_COCCO'
      vname_bec2_diag_3d(2,indxGRAZECOCCO)=
     &  'Grazing rate on Coccos'
      vname_bec2_diag_3d(3,indxGRAZECOCCO)='mmol C/m3/s'
      vname_bec2_diag_3d(4,indxGRAZECOCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCOLOSS)='COCCO_LOSS'
      vname_bec2_diag_3d(2,indxCOCCOLOSS)=
     &  'Cocco non-grazing mortality'
      vname_bec2_diag_3d(3,indxCOCCOLOSS)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxCOCCOLOSS)=' '

      vname_bec2_diag_3d(1,indxCOCCOAGG)='COCCO_AGG'
      vname_bec2_diag_3d(2,indxCOCCOAGG)=
     &  'Aggregation of coccos'
      vname_bec2_diag_3d(3,indxCOCCOAGG)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxCOCCOAGG)=' '

      vname_bec2_diag_3d(1,indxPHOTOCCOCCO)='PHOTOC_COCCO'
      vname_bec2_diag_3d(2,indxPHOTOCCOCCO)=
     &  'C fixation rate by Coccos'
      vname_bec2_diag_3d(3,indxPHOTOCCOCCO)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxPHOTOCCOCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCONLIM)='COCCO_N_LIM'
      vname_bec2_diag_3d(2,indxCOCCONLIM)=
     &  'N limitation factor of Coccos'
      vname_bec2_diag_3d(3,indxCOCCONLIM)=' '
      vname_bec2_diag_3d(4,indxCOCCONLIM)=' '
      vname_bec2_diag_3d(1,indxCOCCOFEUPTAKE)='COCCO_FE_LIM'
      vname_bec2_diag_3d(2,indxCOCCOFEUPTAKE)=
     &  'Fe limitation factor of Coccos'
      vname_bec2_diag_3d(3,indxCOCCOFEUPTAKE)=' '
      vname_bec2_diag_3d(4,indxCOCCOFEUPTAKE)=' '
      vname_bec2_diag_3d(1,indxCOCCOPO4UPTAKE)='COCCO_PO4_LIM'
      vname_bec2_diag_3d(2,indxCOCCOPO4UPTAKE)=
     &  'PO4 limitation factor of Coccos'
      vname_bec2_diag_3d(3,indxCOCCOPO4UPTAKE)=' '
      vname_bec2_diag_3d(4,indxCOCCOPO4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxCOCCOLIGHTLIM)='COCCO_LIGHT_LIM'
      vname_bec2_diag_3d(2,indxCOCCOLIGHTLIM)=
     &  'Cocco light limitation factor'
      vname_bec2_diag_3d(3,indxCOCCOLIGHTLIM)=' '
      vname_bec2_diag_3d(4,indxCOCCOLIGHTLIM)=' '

      vname_bec2_diag_3d(1,indxCACO3PRODCOCCO)='CACO3_PROD_COCCO'
      vname_bec2_diag_3d(2,indxCACO3PRODCOCCO)=
     &  'Production of CaCO3 by Coccos'
      vname_bec2_diag_3d(3,indxCACO3PRODCOCCO)='mmol CaCO3/m3/s '
      vname_bec2_diag_3d(4,indxCACO3PRODCOCCO)=' '

      vname_bec2_diag_3d(1,indxIRONUPTAKECOCCO)='PHOTOFE_COCCO'
      vname_bec2_diag_3d(2,indxIRONUPTAKECOCCO)=
     &  'Iron uptake rate by Coccos'
      vname_bec2_diag_3d(3,indxIRONUPTAKECOCCO)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxIRONUPTAKECOCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCONO3UPTAKE)='COCCO_NO3_UPTAKE'
      vname_bec2_diag_3d(2,indxCOCCONO3UPTAKE)=
     &   'NO3 uptake rate of Coccos'
      vname_bec2_diag_3d(3,indxCOCCONO3UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCONO3UPTAKE)=' '

      vname_bec2_diag_3d(1,indxCOCCONH4UPTAKE)='COCCO_NH4_UPTAKE'
      vname_bec2_diag_3d(2,indxCOCCONH4UPTAKE)=
     &   'NH4 uptake rate of Coccos'
      vname_bec2_diag_3d(3,indxCOCCONH4UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCONH4UPTAKE)=' '

      vname_bec2_diag_3d(1,indxCOCCOGRAZEDIC)='COCCO_GRAZE_DIC'
      vname_bec2_diag_3d(2,indxCOCCOGRAZEDIC)=
     &   'Cocco grazing rate routed to DIC'
      vname_bec2_diag_3d(3,indxCOCCOGRAZEDIC)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCOGRAZEDIC)=' '

      vname_bec2_diag_3d(1,indxCOCCOLOSSDIC)='COCCO_LOSS_DIC'
      vname_bec2_diag_3d(2,indxCOCCOLOSSDIC)=
     &   'Cocco non-grazing mortality routed to DIC'
      vname_bec2_diag_3d(3,indxCOCCOLOSSDIC)='mmol C/m3/s  '
      vname_bec2_diag_3d(4,indxCOCCOLOSSDIC)=' '

      vname_bec2_diag_3d(1,indxGRAZECOCCOZOO)='GRAZE_COCCO_ZOO'
      vname_bec2_diag_3d(2,indxGRAZECOCCOZOO)=
     &  'Grazing rate on Coccos routed to zoo'
      vname_bec2_diag_3d(3,indxGRAZECOCCOZOO)='mmol C/m3/s '
      vname_bec2_diag_3d(4,indxGRAZECOCCOZOO)=' '

      vname_bec2_diag_3d(1,indxQCACO3COCCO)='QCACO3_COCCO'
      vname_bec2_diag_3d(2,indxQCACO3COCCO)=
     &  'Cocco CaCO3/C ratio'
      vname_bec2_diag_3d(3,indxQCACO3COCCO)='mmol CaCO3/mmol C'
      vname_bec2_diag_3d(4,indxQCACO3COCCO)=' '

      vname_bec2_diag_3d(1,indxCOCCOPHOTOACC)='PHOTOACC_COCCO'
      vname_bec2_diag_3d(2,indxCOCCOPHOTOACC)=
     &  'Chl synthesis in photoadaptation for Coccos'
      vname_bec2_diag_3d(3,indxCOCCOPHOTOACC)='mg Chl/m3/s'
      vname_bec2_diag_3d(4,indxCOCCOPHOTOACC)=' '

      vname_bec2_diag_3d(1,indxCOCCOPLIM)='COCCO_P_LIM'
      vname_bec2_diag_3d(2,indxCOCCOPLIM)=
     &  'Cocco P limitation (PO4 + DOP)'
      vname_bec2_diag_3d(3,indxCOCCOPLIM)=' '
      vname_bec2_diag_3d(4,indxCOCCOPLIM)=' '

# endif /* BEC_COCCO */
# ifdef Ncycle_SY
      vname_bec2_diag_3d(1,indxammox)='AMMOX'
      vname_bec2_diag_3d(2,indxammox)=
     &  'Ammonium oxidation (NH4 -> NO2)'
      vname_bec2_diag_3d(3,indxammox)='mmol N/m3/s '
      vname_bec2_diag_3d(4,indxammox)=' '

      vname_bec2_diag_3d(1,indxnitrox)='NITROX'
      vname_bec2_diag_3d(2,indxnitrox)=
     &  'Nitrite oxidation (NO2 -> NO3)'
      vname_bec2_diag_3d(3,indxnitrox)='mmol N/m3/s'
      vname_bec2_diag_3d(4,indxnitrox)=' '

      vname_bec2_diag_3d(1,indxanammox)='ANAMMOX'
      vname_bec2_diag_3d(2,indxanammox)=
     &  'Anammox (NO2+NH4 -> N2)'
      vname_bec2_diag_3d(3,indxanammox)='mmol N/m3/s'
      vname_bec2_diag_3d(4,indxanammox)=' '

      vname_bec2_diag_3d(1,indxDENITRIF1)='DENITRIF1'
      vname_bec2_diag_3d(2,indxDENITRIF1)='Denitrification1'
      vname_bec2_diag_3d(3,indxDENITRIF1)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxDENITRIF1)=' '

      vname_bec2_diag_3d(1,indxDENITRIF2)='DENITRIF2'
      vname_bec2_diag_3d(2,indxDENITRIF2)='Denitrification2'
      vname_bec2_diag_3d(3,indxDENITRIF2)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxDENITRIF2)=' '

      vname_bec2_diag_3d(1,indxDENITRIF3)='DENITRIF3'
      vname_bec2_diag_3d(2,indxDENITRIF3)='Denitrification3'
      vname_bec2_diag_3d(3,indxDENITRIF3)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxDENITRIF3)=' '

      vname_bec2_diag_3d(1,indxN2OAMMOX)='N2OAMMOX'
      vname_bec2_diag_3d(2,indxN2OAMMOX)='N2O produced by Ammox'
      vname_bec2_diag_3d(3,indxN2OAMMOX)='mmol/m3/s'
      vname_bec2_diag_3d(4,indxN2OAMMOX)=''

      vname_bec2_diag_3d(1,indxSPNO2UPTAKE)='SP_NO2_UPTAKE'
      vname_bec2_diag_3d(2,indxSPNO2UPTAKE)=
     &   'NO2 uptake rate of small phyto'
      vname_bec2_diag_3d(3,indxSPNO2UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxSPNO2UPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIATNO2UPTAKE)='DIAT_NO2_UPTAKE'
      vname_bec2_diag_3d(2,indxDIATNO2UPTAKE)=
     &   'NO3 uptake rate of diatoms'
      vname_bec2_diag_3d(3,indxDIATNO2UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIATNO2UPTAKE)=' '

      vname_bec2_diag_3d(1,indxDIAZNO2UPTAKE)='DIAZ_NO2_UPTAKE'
      vname_bec2_diag_3d(2,indxDIAZNO2UPTAKE)=
     &   'NO3 uptake rate of diazotrophs'
      vname_bec2_diag_3d(3,indxDIAZNO2UPTAKE)='mmol/m3/s  '
      vname_bec2_diag_3d(4,indxDIAZNO2UPTAKE)=' '

# endif
# endif /* BEC2_DIAG */

# ifdef BEC2_DIAG

        ! MF: Control BEC2_DIAG Output vars
      wrtavg_bec2_diag_3d=(/.true.,.true.,.false.,.true.,.false.,      !PAR, POC_FLUX_IN, POC_Prod,POC_REMIN, CACO3_FLUX_IN,
     &                     .false.,.false.,.false.,.false.,.false.,       !PCACO3_PROD, CACO3_REMIN,SIO2_FLUX_IN, SIO2_PROD, SIO2_REMIN,
     &                     .false.,.false.,.false.,.false.,.false.,       !DUST_FLUX_IN, DUST_REMIN,P_IRON_FLUX_IN, P_IRON_PROD, P_IRON_REMIN,
     &                     .false.,.false.,.false.,.false.,.false.,       !GRAZE_SP, GRAZE_DIAT,GRAZE_DIAZ, SP_LOSS, DIAT_LOSS,
     &                     .false.,.false.,.false.,.false.,.false.,       !ZOO_LOSS, SP_AGG, DIAT_AGG,PHOTOC_SP, PHOTOC_DIA,
     &                     .true.,.false.,.true.,.true.,.false.,       !TOT_PROD, DOC_PROD,DOC_REMIN, FE_SCAVENGE, SP_N_LIM,
     &                     .false.,.false.,.false.,.false.,.false.,       !SP_FE_LIM, SP_PO4_LIM,SP_LIGHT_LIM, DIAT_N_LIM, DIAT_FE_LIM,
     &                     .false.,.false.,.false.,.false.,.false.,       !DIAT_PO4_LIM, DIAT_SIO3_LIM,DIAT_LIGHT_LIM, CACO3_PROD, DIAZ_NFIX,
     &                     .false.,.false.,.false.,.false.,.false.,       !DIAZ_LOSS, PHOTOC_DIAZ,DIAZ_PO4_LIM, DIAZ_FE_LIM, DIAZ_LIGHT_LIM,
     &                     .false.,.false.,.false.,.false.,.false.,       !FE_SCAVENGE_RATE, DON_PROD,DON_REMIN, DOFE_PROD, DOFE_REMIN,
     &                     .false.,.false.,.false.,.false.,.false.,     !DOP_PROD, DOP_REMIN,DIAT_SI_UPTAKE, PHOTOFE_SP, PHOTOFE_DIA,
     &                     .false.,.false.,.false.,.false.,.false.,    !PHOTOFE_DIAZ, NITRI,DENITRI,SP_NO3_UPTAKE,DIAT_NO3_UPTAKE,
     &                     .false.,.false.,.false.,.false.,.false.,  !DIAZ_NO3_UPTAKE,SP_NH4_UPTAKE, DIAT_NH4_UPTAKE,DIAZ_NH4_UPTAKE,SP_GRAZE_DIC,
     &                     .false.,.false.,.false.,.false.,.false.,  !DIAT_GRAZE_DIC,DIAZ_GRAZE_DIC, SP_LOSS_DIC,DIAT_LOSS_DIC,DIAZ_LOSS_DIC,
     &                     .false.,.false.,.false.,.false.,.false.,  !ZOO_LOSS_DIC, DIAZ_AGG,GRAZE_SP_ZOO, GRAZE_DIAT_ZOO, GRAZE_DIAZ_ZOO,
     &                     .false.,.false.,.false.,.false.,          !SPQCACO3,PHOTOACC_SP,PHOTOACC_DIAT, PHOTOACC_DIAZ,
     &                     .false.,.false.,.false.,.false.,.false.,  !SPC_ZERO, DIATC_ZERO, DIAZC_ZERO, DOC_ZERO, ZOOC_ZERO
     &                     .false.,.false.,.true.,.false.,           !SPCACO3_ZERO, DONR_REMIN, TOTCHL, SP_P_LIM
     &                     .false.,.false.,.false.,.false.,.false.      !DIAT_P_LIM,DIAZ_P_LIM, TOT_PHYTO_C, O2_CONSUMPTION, O2_PRODUCTION
     ! 93 tracer
# ifdef BEC_COCCO
     ! + 18 -> 119 tracer (w/ expl sinking) or 109 tracer (w/ impl sinking)
     &                     ,.true., ! GRAZE_COCCO
     &                     .true.,.true.,.true.,.true.,     ! COCCO_LOSS,COCCO_AGG, PHOTOC_COCCO, COCCO_N_LIM
     &                     .true.,.true.,.true.,.true.,     ! COCCO_PO4_LIM, COCCO_FE_LIM, COCCO_LIGHT_LIM, CACO3_PROD_COCCO
     &                     .true.,.true.,.true.,.true., ! IRON_UPTAKE_COCCO, COCCO_NO3_UPTAKE, COCCO_NH4_UPTAKE, COCCO_GRAZE_DIC
     &                     .true.,.true.,.true.,.true.,    ! COCCO_LOSS_DIC, GRAZE_COCCO_ZOO, QCACO3_COCCO, COCCO_PHOTOACC
     &                     .true.   ! COCCO_P_LIM
# endif
# ifdef Ncycle_SY
     &                     ,.true.,.true.,.false.,.false.,.false.,.false.!AMMOX, NITROX,ANAMMOX,DENITRIF1,DENITRIF2,DENITRIF3
     &                     ,.false.,.false.,.false.,.false. !SP_NO2_UPTAKE,DIAT_NO2_UPTAKE, DIAZ_NO2_UPTAKE, N2OAMMOX,
# endif
     &                     /)

# endif /* BEC2_DIAG */

      tclm_name(iPO4) = 'po4_time'
      tclm_name(iNO3) = 'no3_time'
      tclm_name(iSIO3) = 'sio3_time'
      tclm_name(iO2) = 'o2_time'
      tclm_name(iFE) = 'fe_time'
      tclm_name(iDIC) = 'dic_time'
      tclm_name(iALK) = 'alk_time'
      tclm_name(iNH4) = 'nh4_time'
      tclm_name(iDOC) = 'doc_time'
      tclm_name(iDOP) = 'dop_time'
      tclm_name(iDOP) = 'dop_time'
      tclm_name(iDOPR) = 'dop_time'
      tclm_name(iDON) = 'don_time'
      tclm_name(iDONR) = 'don_time'
      tclm_name(iDOFE) = 'dofe_time'
      tclm_name(iSPC) = 'spc_time'
      tclm_name(iSPCHL) = 'spchl_time'
      tclm_name(iSPCACO3) = 'spcaco3_time'
      tclm_name(iSPFE) = 'spfe_time'
      tclm_name(iDIATC) = 'diatc_time'
      tclm_name(iDIATCHL) = 'diatchl_time'
      tclm_name(iDIATSI) = 'diatsi_time'
      tclm_name(iDIATFE) = 'diatfe_time'
      tclm_name(iZOOC) = 'zooc_time'
      tclm_name(iDIAZC) = 'diazc_time'
      tclm_name(iDIAZCHL) = 'diazchl_time'
      tclm_name(iDIAZFE) = 'diazfe_time'
# ifdef BEC_COCCO
      tclm_name(iCOCCOC) = 'coccoc_time'
      tclm_name(iCOCCOCHL) = 'coccochl_time'
      tclm_name(iCOCCOCAL) = 'coccocal_time'
      tclm_name(iCOCCOFE) = 'coccofe_time'
      tclm_name(iCAL) = 'cal_time'
# endif
# ifdef Ncycle_SY
      tclm_name(iNO2) = 'no2_time'
      tclm_name(iN2) = 'n2_time'
      tclm_name(iN2O) = 'n2o_time'
# endif
#endif /* SOLVE3D */

      end subroutine init_scalars_bec2  !]

! ----------------------------------------------------------------------

#endif /* BIOLOGY_BEC2 - entire module */
      end module bgc
