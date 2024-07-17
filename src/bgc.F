      module bgc

      ! BGC - bio-geo-chemical module
      ! -----------------------------

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Dec)
      ! Refactoring of ETH Zurich roms code with BEC2, code which has been used
      ! by Pierre Damien to run full pacific models.
#include "cppdefs.opt"

#if defined BIOLOGY_BEC2 || defined MARBL

      ! param needed for GLOBAL_2D_array to work. NT = number tracer from param
      use param
      use tracers               ! for iPO4, etc, indices of bgc tracers (formerly in param.h)
#ifdef MARBL      
      use marbl_driver 
#endif      
      use bgc_ecosys_vars
      use nc_read_write
      use roms_read_write
      use netcdf
      use scalars
      use dimensions !
      use ocean_vars
      use dimensions

      implicit none

      private ! Make all variable private to this module unless public specified

      ! Includes:

      ! netcdf outputting:
      integer :: ncid=-1, prev_fill_mode
      real    :: t_avg_bgc=0
      real    :: t_avg_dia_bgc=0
      integer,save :: navg_bgc = 0               ! number of samples in average
      integer,save :: navg_dia_bgc = 0           ! number of samples in average
      integer :: record_avg = nrpf_avg           ! Triggers making of initial file
      integer :: record_his = nrpf_his           
      integer :: record_dia_avg = nrpf_avg_dia   ! Triggers making of initial file
      integer :: record_dia_his = nrpf_his_dia
      real,save :: output_time_his = 0
      real,save :: output_time_avg = 0
      real,save :: output_time_dia_his = 0
      real,save :: output_time_dia_avg = 0

      ! misc from porting ETH code:
      character(len=40) tclm_name(NT) ! Taken from ETH code ncvars

      public set_bgc_surf_frc
      public init_arrays_bgc_frc
      public init_scalars_bec2
      public wrt_bgc_diags
      public wrt_bgc

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


!     call set_frc_var_tile( istr, iend, jstr, jend, nc_dust,    dust,    interp_frc )  ! dust
!     call set_frc_var_tile( istr, iend, jstr, jend, nc_iron,    iron,    interp_frc )  ! iron
      call set_frc_data(nc_dust,dust,'r')
      call set_frc_data(nc_iron,iron,'r')


#ifdef PCO2AIR_FORCING
!     call set_frc_var_tile( istr, iend, jstr, jend, nc_pco2air, pco2air, interp_frc )
      call set_frc_data(nc_pco2air,pco2air,'r')
#ifdef MARBL
      call set_frc_data(nc_pco2air_alt,pco2air_alt,'r')
#endif      
#endif /* PCO2AIR_FORCING */
      
#ifdef NOX_FORCING
      call set_frc_data(nc_nox,nox,'r')
#endif /* NOX_FORCING */   
#ifdef NHY_FORCING
      call set_frc_data(nc_nhy,nhy,'r')
#endif /* NHY_FORCING */


#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
      ! Daily average Short wave radiation
!     call set_frc_var_tile( istr, iend, jstr, jend, nc_swrad_avg, swrad_avg, interp_frc )
      call set_frc_data(nc_swrad_avg,swrad_avg,'r')
     
      do j=jstrR,jendR
        do i=istrR,iendR
          swrad_avg(i,j)= swrad_avg(i,j)/(rho0*Cp)
          
#ifdef SEA_ICE_NOFLUX
          if( t(i,j,N,nrhs,itemp) .le. -1.8 ) then    ! Restrict stflx to prevent surface temperature to go below -2 degrees C.
#   if defined LMD_KPP
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
            swrad_avg(i,j)=0.
#endif                   
#   endif
          endif
#endif

        enddo !<- istrR
      enddo   !<- jstrR
     
#endif /* DAILYPAR_PHOTOINHIBITION */ 

      end subroutine set_bgc_surf_frc  !]
! ----------------------------------------------------------------------
      subroutine wrt_bgc_tracers(avg)  ![
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

      ! import/export
      logical, intent(in) :: avg

      ! local
      character(len=99),save :: fname_his,fname_avg
      integer                :: itrc, itavg, ierr
      logical,save           :: first_step_avg=.true.
      logical,save           :: first_step_his=.true.


      if (avg) then

!        if (.not. first_step_avg) then
           call calc_avg_bgc                          ! don't include t=0 in averaging
           output_time_avg = output_time_avg + dt     ! only start count after first timestep
!        endif
!        first_step_avg=.false.
        

!        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,F11.4,1x,A,I4)')
!     &       'wrt_bgc_tracers :: average, time =', output_time_avg,
!     &       'period =', output_period_avg, 'rec =', record_avg
!        endif

        if (output_time_avg>=output_period_avg) then  ! time for an output
          if (mod(record_avg,nrpf_avg)==0) then
             call create_bgc_file(fname_avg,avg)
             record_avg = 0
          endif
          record_avg = record_avg + 1
          output_time_avg = 0
          navg_bgc = 0
          ierr=nf90_open(fname_avg,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_avg/))
          do itrc=iTandS+nt_passive+1, nt
             if (wrt_t_avg(itrc)) then
                if (mynode == 0) then
                   write(*,'(7x,A,1x,A)')
     &             'wrt_bgc_tracers :: average, trc = ', t_vname(itrc) 
                endif
             itavg = NT_2_t_avg(itrc)                         ! get respective index for t_avg(itavg) -> t(itrc)
             call ncwrite(ncid,t_vname(itrc),t_avg(i0:i1,j0:j1,:,itavg),(/1,1,1,record_avg/))
             endif
          enddo
          if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'bgc :: wrote average, tdays =', tdays,
     &       'step =', iic, 'rec =', record_avg
          endif
        endif  

      else

        if (.not. first_step_his) then
           output_time_his = output_time_his + dt
        endif
        first_step_his=.false.

!        if (mynode == 0) then
!           write(*,'(7x,A,1x,F11.4,2x,A,F11.4,1x,A,I4)')
!     &       'wrt_bgc_tracers :: history, time =', output_time_his,
!     &       'period =', output_period_his, 'rec =', record_his
!        endif

        if (output_time_his>=output_period_his .or.
     &      output_time_his==0                 ) then  ! time for an output
          if (mod(record_his,nrpf_his)==0) then
             call create_bgc_file(fname_his,avg)
             record_his = 0
          endif
          record_his = record_his + 1
          output_time_his=0
          ierr=nf90_open(fname_his,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_his/))
          do itrc=iTandS+nt_passive+1, nt
            if (wrt_t(itrc)) then
                if (mynode == 0) then
                   write(*,'(7x,A,1x,A)')
     &             'wrt_bgc_tracers :: history , trc = ', t_vname(itrc)
                endif
!            call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nstp,itrc),(/1,1,1,record_his/))
            call ncwrite(ncid,t_vname(itrc),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record_his/))
            endif
          enddo
          if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'bgc :: wrote history, tdays =', tdays,
     &       'step =', iic-1, 'rec =', record_his
          endif
        endif  ! <-- wrt_file_his

      endif

      ierr=nf90_close(ncid)

      end subroutine wrt_bgc_tracers  !]
! ----------------------------------------------------------------------
      subroutine wrt_bgc_diags ![
      implicit none
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
      ! local
      character(len=99),save :: fname_his,fname_avg
      integer                :: itrc, ierr, idiag
      logical,save           :: first_step_dia_avg=.true.
      logical,save           :: first_step_dia_his=.true.

      if (wrt_avg_dia) then

!         if (.not. first_step_dia_avg) then
            call calc_avg_dia_bgc                              ! don't include t=0 in averaging
            output_time_dia_avg = output_time_dia_avg + dt     ! only start count after first timestep
!         endif
!         first_step_dia_avg=.false.

        if (output_time_dia_avg>=output_period_avg_dia) then  ! time for an output
          if (mod(record_dia_avg,nrpf_avg_dia)==0) then
             call create_bgc_dia_file(fname_avg,.true.)
             record_dia_avg = 0
          endif
          record_dia_avg = record_dia_avg + 1
          output_time_dia_avg = 0
          navg_dia_bgc = 0
          ierr=nf90_open(fname_avg,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_dia_avg/))
          do itrc=1,nr_bec2_diag_2d
            if (wrt_bec2_diag_2d(itrc)) then
            idiag = idx_bec2_diag_2d(itrc)
            call ncwrite(ncid,vname_bec2_diag_2d(1,itrc),
     &           bec2_diag_2d_avg(i0:i1,j0:j1,idiag),(/1,1,record_dia_avg/))
            endif
          enddo
          do itrc=1,nr_bec2_diag_3d
            if (wrt_bec2_diag_3d(itrc)) then
            idiag = idx_bec2_diag_3d(itrc)
            call ncwrite(ncid,vname_bec2_diag_3d(1,itrc),
     &           bec2_diag_3d_avg(i0:i1,j0:j1,:,idiag),(/1,1,1,record_dia_avg/))
            endif
          enddo
          if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'bgc diag :: wrote average, tdays =', tdays,
     &       'step =', iic, 'rec =', record_dia_avg
          endif
        endif  

      endif  

      if (wrt_his_dia) then

        if (.not. first_step_dia_his) then
           output_time_dia_his = output_time_dia_his + dt
        endif
        first_step_dia_his=.false.

        if (output_time_dia_his>=output_period_his_dia .or.
     &      output_time_dia_his==0                     ) then  ! time for an output
          if (mod(record_dia_his,nrpf_his_dia)==0) then
             call create_bgc_dia_file(fname_his,.false.)
             record_dia_his = 0
          endif
          record_dia_his = record_dia_his + 1
          output_time_dia_his = 0
          ierr=nf90_open(fname_his,nf90_write,ncid)
          call ncwrite(ncid,'ocean_time',(/time/),(/record_dia_his/))
          do itrc=1,nr_bec2_diag_2d
            if (wrt_bec2_diag_2d(itrc)) then
!               if (mynode == 0) then
!                   write(*,'(7x,A,1x,A)')
!     &            'wrt_bgc_dia :: history, diag = ', vname_bec2_diag_2d(1,itrc)
!               endif
            idiag = idx_bec2_diag_2d(itrc)
            call ncwrite(ncid,vname_bec2_diag_2d(1,itrc),
     &           bec2_diag_2d(i0:i1,j0:j1,idiag),(/1,1,record_dia_his/))
            endif
          enddo
          do itrc=1,nr_bec2_diag_3d
            if (wrt_bec2_diag_3d(itrc)) then
!               if (mynode == 0) then
!                   write(*,'(7x,A,1x,A)')
!     &            'wrt_bgc_dia :: history, diag = ', vname_bec2_diag_3d(1,itrc)
!               endif
            idiag = idx_bec2_diag_3d(itrc)
            call ncwrite(ncid,vname_bec2_diag_3d(1,itrc),
     &           bec2_diag_3d(i0:i1,j0:j1,:,idiag),(/1,1,1,record_dia_his/))
            endif
          enddo
          if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')
     &       'bgc diag :: wrote history, tdays =', tdays,
     &       'step =', iic, 'rec =', record_dia_his
          endif
        endif  ! <-- wrt_file_his

      endif

      ierr=nf90_close(ncid)
#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */

      end subroutine wrt_bgc_diags  !]
!----------------------------------------------------------------------
      subroutine create_bgc_file(fname,avg)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ncid,ierr,varid

      if (avg) then
        call create_file('_bgc_avg',fname)
      else
        call create_file('_bgc',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      call def_vars_bgc( ncid,avg )

      ierr = nf90_enddef(ncid)

      end subroutine create_bgc_file !]
!----------------------------------------------------------------------
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
      subroutine create_bgc_dia_file(fname,avg)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname
      logical,          intent(in)  :: avg                 ! his or average file

      ! local
      integer :: ncid,ierr,varid

      if (avg) then
        call create_file('_bgc_dia_avg',fname)
      else
        call create_file('_bgc_dia',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      call def_bgc_diag(ncid,avg)

      ierr = nf90_enddef(ncid)

      end subroutine create_bgc_dia_file !]
#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */
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

      if (avg) then
        do itrc=iTandS+nt_passive+1,NT
          if (wrt_t_avg(itrc)) then
             text_lname='avg_'/ /t_lname(itrc)
             varid = nccreate(ncid,t_vname(itrc),
     &         (/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))

             ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
             ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
          endif
        enddo
      else
        do itrc=iTandS+nt_passive+1,NT
          if (wrt_t(itrc)) then
             text_lname='avg_'/ /t_lname(itrc)
             varid = nccreate(ncid,t_vname(itrc),
     &         (/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
             ierr = nf90_put_att(ncid,varid,'long_name',text_lname)
             ierr = nf90_put_att(ncid,varid,'units',t_units(itrc))
          endif
        enddo
      endif

      end subroutine def_vars_bgc  !]
! ----------------------------------------------------------------------
      subroutine calc_avg_bgc  ![
      implicit none

      ! local
      real    :: coef
      integer :: itrc, itavg, k

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
        endif
      enddo

      end subroutine calc_avg_bgc  !]
! ----------------------------------------------------------------------
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
      subroutine calc_avg_dia_bgc  ![
      implicit none

      ! local
      real    :: coef
      integer :: itrc, itavg, k

      navg_dia_bgc = navg_dia_bgc +1
      coef = 1./navg_dia_bgc

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'bgc diag :: started averaging. ',
     &   'output_period_avg_dia (s) =', output_period_avg_dia
      endif

      t_avg_dia_bgc = t_avg_dia_bgc*(1-coef) + time*coef

      do itrc=1,nr_bec2_wrdiag_2d
        bec2_diag_2d_avg(i0:i1,j0:j1,itrc) = bec2_diag_2d_avg(i0:i1,j0:j1,itrc) *(1-coef)
     &                                     + bec2_diag_2d(i0:i1,j0:j1,itrc)     *coef
      end do
      do itrc=1,nr_bec2_wrdiag_3d
        do k=1,N
          bec2_diag_3d_avg(i0:i1,j0:j1,k,itrc) = bec2_diag_3d_avg(i0:i1,j0:j1,k,itrc) *(1-coef)
     &                                         + bec2_diag_3d(i0:i1,j0:j1,k,itrc)     *coef
        end do
      end do

      end subroutine calc_avg_dia_bgc  !]
#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */
! ----------------------------------------------------------------------
#if defined (BEC2_DIAG) || defined (MARBL_DIAGS)
      subroutine def_bgc_diag(ncid,avg)  ![

      ! Define history/avg file variables for BEC2_DIAG
      ! Taken from def_his.F of ETH code.
      ! Called from def_his.F

      implicit none

      ! input
      integer,intent(in) :: ncid
      logical,intent(in) :: avg
      ! local
      integer :: idiag, ierr, varid

      ! 2d diagnostics:
      do idiag=1,nr_bec2_diag_2d
        if (wrt_bec2_diag_2d(idiag)) then
          varid = nccreate(ncid,vname_bec2_diag_2d(1,idiag),
     &            (/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',vname_bec2_diag_2d(2,idiag))
          ierr = nf90_put_att(ncid,varid,'units',vname_bec2_diag_2d(3,idiag))
        endif
      enddo

      ! 3d diagnostics:
      do idiag=1,nr_bec2_diag_3d
        if (wrt_bec2_diag_3d(idiag)) then
          varid = nccreate(ncid,vname_bec2_diag_3d(1,idiag),
     &            (/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',vname_bec2_diag_3d(2,idiag))
          ierr = nf90_put_att(ncid,varid,'units',vname_bec2_diag_3d(3,idiag))
        endif
      enddo

      end subroutine def_bgc_diag  !]
#endif /* (BEC2_DIAG) || defined (MARBL_DIAGS) */
! ----------------------------------------------------------------------
      subroutine init_arrays_bgc_frc  ![
      implicit none

      allocate( nc_dust%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_iron%vdata( GLOBAL_2D_ARRAY,2) )

#ifdef PCO2AIR_FORCING
      allocate( nc_pco2air%vdata( GLOBAL_2D_ARRAY,2) )
#ifdef MARBL      
      allocate( nc_pco2air_alt%vdata( GLOBAL_2D_ARRAY,2) )
#endif      
#endif

#ifdef NOX_FORCING
      allocate( nc_nox%vdata( GLOBAL_2D_ARRAY,2) )
#endif
#ifdef NHY_FORCING
      allocate( nc_nhy%vdata( GLOBAL_2D_ARRAY,2) )
#endif

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
      allocate( nc_swrad_avg%vdata( GLOBAL_2D_ARRAY,2) )
#endif

      end subroutine  init_arrays_bgc_frc !]
! ----------------------------------------------------------------------
      subroutine init_scalars_bec2()  ![
#ifndef MARBL
!     
! Set initial values for  globally accessible (stored in common
! blocks) scalar variables of the BEC model.
!
      use comm_vars
      use scalars

      implicit none
      integer i,j, itrc ,lvar,lenstr ! ierr, DevinD removed ierr as arguement
      integer omp_get_num_threads

      if (wrt_his .and. mod(output_period_his,dt) /= 0) then
        stop 'ERROR :: bgc.opt - dt not a factor of output_period_his!'
      endif
      if (wrt_avg .and. mod(output_period_avg,dt) /= 0) then
        stop 'ERROR :: bgc.opt - dt not a factor of output_period_avg!'
      endif
#ifdef BEC2_DIAG
      if (wrt_his_dia .and. mod(output_period_his_dia,dt) /= 0) then
        stop 'ERROR :: bgc.opt - dt not a factor ' //
     &        'of output_period_his_dia for diagnostics!'
      endif
      if (wrt_avg_dia .and. mod(output_period_avg_dia,dt) /= 0) then
        stop 'ERROR :: bgc.opt - dt not a factor' //
     &        'of output_period_avg_dia for diagnostics!'
      endif
#endif /* BEC2_DIA */

#ifdef SOLVE3D

! set the bgc tracer index for the bgc routines
      po4_ind_t     = iPO4-iTandS-nt_passive     
      no3_ind_t     = iNO3-iTandS-nt_passive
      sio3_ind_t    = iSIO3-iTandS-nt_passive
      nh4_ind_t     = iNH4-iTandS-nt_passive
      fe_ind_t      = iFE-iTandS-nt_passive
      o2_ind_t      = iO2-iTandS-nt_passive
      dic_ind_t     = iDIC-iTandS-nt_passive
      alk_ind_t     = iALK-iTandS-nt_passive
      doc_ind_t     = iDOC-iTandS-nt_passive
      don_ind_t     = iDon-iTandS-nt_passive
      dofe_ind_t    = iDofe-iTandS-nt_passive
      dop_ind_t     = iDop-iTandS-nt_passive
      dopr_ind_t    = iDopr-iTandS-nt_passive
      donr_ind_t    = iDonr-iTandS-nt_passive
      zooc_ind_t    = iZOOC-iTandS-nt_passive
      spc_ind_t     = iSPC-iTandS-nt_passive
      spchl_ind_t   = iSPCHL-iTandS-nt_passive
      spfe_ind_t    = iSPFE-iTandS-nt_passive
      spcaco3_ind_t = iSPCACO3-iTandS-nt_passive
      diatc_ind_t   = iDIATC-iTandS-nt_passive
      diatchl_ind_t = iDIATCHL-iTandS-nt_passive
      diatfe_ind_t  = iDIATFE-iTandS-nt_passive
      diatsi_ind_t  = iDIATSI-iTandS-nt_passive
      diazc_ind_t   = iDiazc-iTandS-nt_passive
      diazchl_ind_t = iDiazchl-iTandS-nt_passive
      diazfe_ind_t  = iDiazfe-iTandS-nt_passive
#ifdef Ncycle_SY
      no2_ind_t     = iNO2-iTandS-nt_passive
      n2_ind_t      = iN2-iTandS-nt_passive
      n2o_ind_t     = iN2O-iTandS-nt_passive
#endif

# ifdef BEC2_DIAG
#include "bgc_2Ddiagnostics.opt"
#include "bgc_3Ddiagnostics.opt"
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
# ifdef Ncycle_SY
      tclm_name(iNO2) = 'no2_time'
      tclm_name(iN2) = 'n2_time'
      tclm_name(iN2O) = 'n2o_time'
# endif
#endif /* SOLVE3D */

#endif /* ndef MARBL */      
      end subroutine init_scalars_bec2  !]
! ----------------------------------------------------------------------
      subroutine wrt_bgc(init)  ![
      implicit none

! import/export
      logical,optional       :: init

      ! initialisation writting
      if (wrt_his.and.present(init)) then
         call wrt_bgc_tracers(.false.) ! write his data
      else
      ! local
         if (wrt_his)     call wrt_bgc_tracers(.false.) ! write his data
         if (wrt_avg)     call wrt_bgc_tracers(.true.) ! write avg data
      endif

#if defined BEC2_DIAG || defined MARBL_DIAGS
      call wrt_bgc_diags
#endif 


      end subroutine wrt_bgc  !]
! ----------------------------------------------------------------------
#endif /*(BIOLOGY_BEC2 || MARBL)*/

      end module bgc
