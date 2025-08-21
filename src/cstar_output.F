      module cstar_output
      ! Collection of variables dedicated to the C-star system
      ! Jeroen Molemaker Feb 2025

#include "cppdefs.opt"

#if defined MARBL && defined MARBL_DIAGS
      use param
      use tracers
      use marbl_driver
      use bgc_shared_vars
      use dimensions
      use roms_read_write
      use nc_read_write
      use netcdf
      use scalars
      use ocean_vars
      use diagnostics
      use cdr_frc
      implicit none

      private

#include "cstar_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation

      integer,dimension(6) :: date
      character(len=15)  :: datestr
      integer :: prev_month

      integer :: navg = 0
      integer :: iALK, iDIC, iALK_alt, iDIC_alt
      integer :: iPH, iPH_alt, iFG, iFG_alt
      real,allocatable,dimension(:,:,:) :: hALK_tmp
      real,allocatable,dimension(:,:,:) :: hDIC_tmp
      real,allocatable,dimension(:,:,:) :: hALK_alt_tmp
      real,allocatable,dimension(:,:,:) :: hDIC_alt_tmp

      real,allocatable,dimension(:,:,:) :: ALK_source
      real,allocatable,dimension(:,:,:) :: ALK_alt_source
      real,allocatable,dimension(:,:,:) :: DIC_source
      real,allocatable,dimension(:,:,:) :: DIC_alt_source

      ! Needed for averaging
      real,allocatable,dimension(:,:) :: zeta__avg
      real,allocatable,dimension(:,:,:) :: temp_avg
      real,allocatable,dimension(:,:,:) :: salt_avg
      real,allocatable,dimension(:,:,:) :: ALK_avg
      real,allocatable,dimension(:,:,:) :: hALK_avg
      real,allocatable,dimension(:,:,:) :: DIC_avg
      real,allocatable,dimension(:,:,:) :: hDIC_avg
      real,allocatable,dimension(:,:,:) :: ALK_alt_avg
      real,allocatable,dimension(:,:,:) :: hALK_alt_avg
      real,allocatable,dimension(:,:,:) :: DIC_alt_avg
      real,allocatable,dimension(:,:,:) :: hDIC_alt_avg
      real,allocatable,dimension(:,:,:) :: pH_avg
      real,allocatable,dimension(:,:,:) :: pH_alt_avg
      real,allocatable,dimension(:,:) :: FG_CO2_avg
      real,allocatable,dimension(:,:) :: FG_ALT_CO2_avg
      real,allocatable,dimension(:,:,:) :: ALK_source_avg
      real,allocatable,dimension(:,:,:) :: ALK_alt_source_avg
      real,allocatable,dimension(:,:,:) :: DIC_source_avg
      real,allocatable,dimension(:,:,:) :: DIC_alt_source_avg

      ! Public functions
      public wrt_cstar,init_cstar

      contains
!----------------------------------------------------------------------
      subroutine init_cstar ![
      ! Allocate and initialize arrays.
      implicit none

      ! local
      logical,save :: done=.false.
      integer :: itot=0
      integer :: idx

      if (done) then
        return
      else
#ifndef MARBL
        error stop 'ERROR: cstar_output must have MARBL enabled.'
#endif
        done = .true.
      endif

      if (monthly_averages .and. .not. do_avg) then
         error stop 'ERROR: cstar_output:: monthly avgs=t, but do_avg=f'
      endif

      if (monthly_averages) then
          call sec2date(time,date)
          prev_month = date(2)
      endif

      if (mynode==0) print *,'init random'

      ! Loop over MARBL tracers...
      ! ...but use itot as index (to include other tracers)
      do idx=1,ntrc_bio
         itot=itot+1
         if (t_vname(itot)=='ALK') then
           iALK = itot
         endif
         if (t_vname(itot)=='DIC') then
           iDIC = itot
         endif
         if (t_vname(itot)=='ALK_ALT_CO2') then
           iALK_alt = itot
         endif
         if (t_vname(itot)=='DIC_ALT_CO2') then
           iDIC_alt = itot
         endif
      enddo

      itot = 0
      ! Loop over 2D BGC diagnostics...
      ! ...but use itot as index (to include other diagnostics)
      do idx=1,nr_bgc_diag_2d
         itot=itot+1
         if (vname_bgc_diag_2d(1,idx)=='FG_CO2') then
           iFG = itot
         endif
         if (vname_bgc_diag_2d(1,idx)=='FG_ALT_CO2') then
           iFG_alt = itot
         endif
      enddo

      itot = 0
      ! Loop over 3D BGC diagnostics...
      ! ...but use itot as index (to include other diagnostics)
      do idx=1,nr_bgc_diag_3d
         itot=itot+1
         if (vname_bgc_diag_3d(1,idx)=='pH_3D') then
           iPH = itot
         endif
         if (vname_bgc_diag_3d(1,idx)=='pH_3D_ALT_CO2') then
           iPH_alt = itot
         endif
      enddo

      ! put the relevant part of your code here

        if (cdr_source) then
          allocate(ALK_source(GLOBAL_2D_ARRAY,1:N) )
          allocate(ALK_alt_source(GLOBAL_2D_ARRAY,1:N) )
          allocate(DIC_source(GLOBAL_2D_ARRAY,1:N) )
          allocate(DIC_alt_source(GLOBAL_2D_ARRAY,1:N) )
        endif

      if (do_avg) then
        allocate(zeta__avg(GLOBAL_2D_ARRAY) )
        allocate(temp_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(salt_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(ALK_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(hALK_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(DIC_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(hDIC_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(ALK_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(hALK_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(DIC_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(hDIC_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(pH_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(pH_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        allocate(FG_CO2_avg(GLOBAL_2D_ARRAY) )
        allocate(FG_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        if (cdr_source) then
          allocate(ALK_source_avg(GLOBAL_2D_ARRAY,1:N) )
          allocate(ALK_alt_source_avg(GLOBAL_2D_ARRAY,1:N) )
          allocate(DIC_source_avg(GLOBAL_2D_ARRAY,1:N) )
          allocate(DIC_alt_source_avg(GLOBAL_2D_ARRAY,1:N) )
        endif
      else
        allocate(hALK_tmp(GLOBAL_2D_ARRAY,1:N) )
        allocate(hDIC_tmp(GLOBAL_2D_ARRAY,1:N) )
        allocate(hALK_alt_tmp(GLOBAL_2D_ARRAY,1:N) )
        allocate(hDIC_alt_tmp(GLOBAL_2D_ARRAY,1:N) )
      endif

      if (mynode==0) print *,'init cstar'

      end subroutine init_cstar  !]
!----------------------------------------------------------------------
      subroutine calc_average ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      use param
      implicit none

      ! local
      real :: coef

      navg = navg+1

      coef = 1./navg

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'cstar :: started averaging. ',
     &   'output_period (s) =', output_period
      endif

      zeta__avg(:,:) = zeta__avg(:,:)*(1-coef) + zeta(:,:,knew)*coef

      temp_avg(:,:,:) = temp_avg(:,:,:)*(1-coef) + t(:,:,:,knew,itemp)*coef

      salt_avg(:,:,:) = salt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,isalt)*coef

      ALK_avg(:,:,:) = ALK_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK)*coef

      hALK_avg(:,:,:) = hALK_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK)*Hz(:,:,:)*coef

      DIC_avg(:,:,:) = DIC_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC)*coef

      hDIC_avg(:,:,:) = hDIC_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC)*Hz(:,:,:)*coef

      ALK_alt_avg(:,:,:) = ALK_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK_alt)*coef

      hALK_alt_avg(:,:,:) = hALK_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK_alt)*Hz(:,:,:)*coef

      DIC_alt_avg(:,:,:) = DIC_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC_alt)*coef

      hDIC_alt_avg(:,:,:) = hDIC_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC_alt)*Hz(:,:,:)*coef

      pH_avg(:,:,:) = pH_avg(:,:,:)*(1-coef) + bgc_diag_3d(:,:,:,iPH)*coef

      pH_alt_avg(:,:,:) = pH_alt_avg(:,:,:)*(1-coef) + bgc_diag_3d(:,:,:,iPH_alt)*coef

      FG_CO2_avg(:,:) = FG_CO2_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,iFG)*coef

      FG_ALT_CO2_avg(:,:) = FG_ALT_CO2_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,iFG_alt)*coef

      if (cdr_source) then
        ALK_source_avg(:,:,:) = ALK_source_avg(:,:,:)*(1-coef) + ALK_source(:,:,:)*coef

        ALK_alt_source_avg(:,:,:) = ALK_alt_source_avg(:,:,:)*(1-coef) + ALK_alt_source(:,:,:)*coef

        DIC_source_avg(:,:,:) = DIC_source_avg(:,:,:)*(1-coef) + DIC_source(:,:,:)*coef

        DIC_alt_source_avg(:,:,:) = DIC_alt_source_avg(:,:,:)*(1-coef) + DIC_alt_source(:,:,:)*coef
      endif

      end subroutine calc_average !]
! ----------------------------------------------------------------------
      subroutine multiply_by_thickness ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      use param
      implicit none

      hALK_tmp(:,:,:) = t(i0:i1,j0:j1,:,knew,iALK)*Hz(i0:i1,j0:j1,:)

      hDIC_tmp(:,:,:) = t(i0:i1,j0:j1,:,knew,iDIC)*Hz(i0:i1,j0:j1,:)

      hALK_alt_tmp(:,:,:) = t(i0:i1,j0:j1,:,knew,iALK_alt)*Hz(i0:i1,j0:j1,:)

      hDIC_alt_tmp(:,:,:) = t(i0:i1,j0:j1,:,knew,iDIC_alt)*Hz(i0:i1,j0:j1,:)

      end subroutine multiply_by_thickness !]
!----------------------------------------------------------------------
      subroutine calc_cdr_source ![
      ! Update source terms from the CDR module
      use param
      implicit none

      integer :: i,j,k,icdr,cidx

      ALK_source(:,:,:) = 0
      ALK_alt_source(:,:,:) = 0
      DIC_source(:,:,:) = 0
      DIC_alt_source(:,:,:) = 0

      ! Loop over cdr release locations in this subdomain
      do cidx=1,cdr_nprf
        icdr = cdr_icdr(cidx)
        i = cdr_iloc(cidx)
        j = cdr_jloc(cidx)
        do k=1,nz
          ALK_source(i,j,k) = ALK_source(i,j,k)
     &      +cdr_prf(cidx,k)*cdr_flx(icdr,iALK)
          ALK_alt_source(i,j,k) = ALK_alt_source(i,j,k)
     &      + cdr_prf(cidx,k)*cdr_flx(icdr,iALK_alt)
          DIC_source(i,j,k) = DIC_source(i,j,k)
     &      + cdr_prf(cidx,k)*cdr_flx(icdr,iDIC)
          DIC_alt_source(i,j,k) = DIC_alt_source(i,j,k)
     &      + cdr_prf(cidx,k)*cdr_flx(icdr,iDIC_alt)
        enddo
      enddo

      end subroutine calc_cdr_source !]

! ----------------------------------------------------------------------

      subroutine def_vars_cstar(ncid)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer                        :: ierr, varid

      varid = nccreate(ncid,'zeta',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','sea surface heigth')
      ierr = nf90_put_att(ncid,varid,'units','m')

      varid = nccreate(ncid,'temp',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','temperature')
      ierr = nf90_put_att(ncid,varid,'units','degrees C')

      varid = nccreate(ncid,'salt',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','salinity')
      ierr = nf90_put_att(ncid,varid,'units','PSU')

      varid = nccreate(ncid,'ALK',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','alkalinity')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'DIC',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','DIC')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'ALK_ALT_CO2',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','alt alkalinity')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'DIC_ALT_CO2',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','alt DIC')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'hALK',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','thickness-weighted alkalinity')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'hDIC',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','thickness-weighted DIC')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'hALK_ALT_CO2',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','alt thickness-weighted alkalinity')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'hDIC_ALT_CO2',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','alt thickness-weighted DIC')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'pH',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','pH')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'pH_ALT_CO2',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','alt pH')
      ierr = nf90_put_att(ncid,varid,'units','nondimensional')

      varid = nccreate(ncid,'FG_CO2',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','surface flux of CO2')
      ierr = nf90_put_att(ncid,varid,'units','m')

      varid = nccreate(ncid,'FG_ALT_CO2',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','alt surface flux of CO2')
      ierr = nf90_put_att(ncid,varid,'units','m')

      if (cdr_source) then
        varid = nccreate(ncid,'ALK_source',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name','ALK source from CDR module')
        ierr = nf90_put_att(ncid,varid,'units','meq/s')

        varid = nccreate(ncid,'ALK_ALT_source',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name','alt ALK source from CDR module')
        ierr = nf90_put_att(ncid,varid,'units','meq/s')

        varid = nccreate(ncid,'DIC_source',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name','DIC source from CDR module')
        ierr = nf90_put_att(ncid,varid,'units','mmol/s')

        varid = nccreate(ncid,'DIC_ALT_source',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'long_name','alt DIC source from CDR module')
        ierr = nf90_put_att(ncid,varid,'units','mmol/s')
      endif



      if (do_avg) then
        ierr=nf90_put_att(ncid,nf90_global,'averaging',
     &      'All variables are averaged in time')
      endif

      end subroutine def_vars_cstar  !]
!----------------------------------------------------------------------
      subroutine wrt_cstar  ![
      ! Check whether it is time to write to file
      implicit none

      if (cdr_source) call calc_cdr_source

      if (do_avg) call calc_average

      if (monthly_averages) then
        call sec2date(time+dt,date)

        if ((date(2) - prev_month) /= 0) call wrt_output

        prev_month = date(2)
      else

        output_time = output_time + dt

        if (output_time>=output_period) then
          call wrt_output
          output_time = 0
        endif

      endif

      end subroutine wrt_cstar  !]
!----------------------------------------------------------------------
      subroutine wrt_output  ![
      ! Call wrt after completion of the time-step
      implicit none

      ! local
      character(len=99),save :: fname
      integer,dimension(3)   :: start
      integer                :: ncid,ierr

      if (record==nrpf) then
          call create_file('_cstar',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_cstar(ncid)
          ierr = nf90_close(ncid)
          record = 0
      endif

        record = record+1

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr)
     &      call handle_ierr(ierr,'wrt_cstar, opening: ',fname)

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        if (do_avg) then
          call ncwrite(ncid,'zeta'  ,zeta__avg(i0:i1,j0:j1),(/1,1,record/))
          call ncwrite(ncid,'temp',temp_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'salt',salt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'ALK',ALK_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'DIC',DIC_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'ALK_ALT_CO2',ALK_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'DIC_ALT_CO2',DIC_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'hALK',hALK_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'hDIC',hDIC_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'hALK_ALT_CO2',hALK_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'hDIC_ALT_CO2',hDIC_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'pH',pH_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'pH_ALT_CO2',pH_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'FG_CO2'  ,FG_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
          call ncwrite(ncid,'FG_ALT_CO2'  ,FG_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
          if (cdr_source) then
            call ncwrite(ncid,'ALK_source',ALK_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'ALK_ALT_source',ALK_alt_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'DIC_source',DIC_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'DIC_ALT_source',DIC_alt_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          endif
        else
          call multiply_by_thickness

          call ncwrite(ncid,'zeta'  ,zeta(i0:i1,j0:j1,knew),(/1,1,record/))
          call ncwrite(ncid,'temp',t(i0:i1,j0:j1,:,nnew,itemp),(/1,1,1,record/))
          call ncwrite(ncid,'salt',t(i0:i1,j0:j1,:,nnew,isalt),(/1,1,1,record/))
          call ncwrite(ncid,'ALK',t(i0:i1,j0:j1,:,nnew,iALK),(/1,1,1,record/))
          call ncwrite(ncid,'DIC',t(i0:i1,j0:j1,:,nnew,iDIC),(/1,1,1,record/))
          call ncwrite(ncid,'ALK_ALT_CO2',t(i0:i1,j0:j1,:,nnew,iALK_alt),(/1,1,1,record/))
          call ncwrite(ncid,'DIC_ALT_CO2',t(i0:i1,j0:j1,:,nnew,iDIC_alt),(/1,1,1,record/))
          call ncwrite(ncid,'hALK',hALK_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'hDIC',hDIC_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'hALK_ALT_CO2',hALK_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'hDIC_ALT_CO2',hDIC_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
          call ncwrite(ncid,'pH',bgc_diag_3d(i0:i1,j0:j1,:,iPH),(/1,1,1,record/))
          call ncwrite(ncid,'pH_ALT_CO2',bgc_diag_3d(i0:i1,j0:j1,:,iPH_alt),(/1,1,1,record/))
          call ncwrite(ncid,'FG_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG),(/1,1,record/))
          call ncwrite(ncid,'FG_ALT_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG_alt),(/1,1,record/))
          if (cdr_source) then
            call ncwrite(ncid,'ALK_source',ALK_source(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'ALK_ALT_source',ALK_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'DIC_source',DIC_source(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'DIC_ALT_source',DIC_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/))
          endif
        endif

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_cstar :: wrote cstar, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        navg = 0

!
      end subroutine wrt_output !]
!----------------------------------------------------------------------

#else /* MARBL && MARBL_DIAGS */

!----------------------------------------------------------------------
      use param
      use tracers
      use dimensions
      use roms_read_write
      use scalars
      use ocean_vars

      implicit none

      private

#include "cstar_output.opt"

      ! Public functions
      public init_cstar

      contains

      subroutine init_cstar ![
      ! Allocate and initialize arrays.
      implicit none

#ifndef MARBL
        error stop 'ERROR: cstar_output must have MARBL enabled.'
#endif

#ifndef MARBL_DIAGS
        error stop 'ERROR: cstar_output must have MARBL_DIAGS enabled.'
#endif

      end subroutine init_cstar !]
!----------------------------------------------------------------------

#endif /* MARBL && MARBL_DIAGS */

      end module cstar_output
